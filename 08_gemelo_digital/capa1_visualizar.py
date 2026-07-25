#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa1_visualizar.py  --  Bio-Kidney AI 2026
============================================
Visualizacion de las NEFRONAS (Capa 1) superpuestas sobre el DOMINIO (Capa 0).

  >>> Se ejecuta con bpy (Blender), NO en env_biokidney <<<

Carga capa0_dominio.npz (contexto tenue de fondo) y capa1_nefronas.npz
(protagonista) y dibuja:
  - Dominio de Capa 0 como nube de puntos muy submuestreada y gris tenue.
  - Hilio / seno renal como esfera roja en (0,-32,0).
  - 1300 glomerulos como puntos coloreados por tipo
        cortical -> azul,  yuxtamedular -> naranja/rojo.
  - Tubulos (polilineas de 6 vertices) como tubos coloreados por tipo (tenues).

Todo con objetos de malla unicos + atributo de color de vertice + Geometry
Nodes (Mesh to Points / Mesh to Curve), igual que en capa0_visualizar.py.

Uso:
    blender --python 08_gemelo_digital/capa1_visualizar.py
"""

import os
import bpy
import numpy as np

# ============================================================================
#  PARAMETROS
# ============================================================================
DOM_SUBSAMPLE = 20          # dominio de fondo: 1 de cada N puntos (muy tenue)
DOM_POINT_RADIUS = 0.22
DOM_COLOR = (0.55, 0.55, 0.55)   # gris
DOM_EMISSION = 0.12              # muy bajo -> no compite con las nefronas

GLOM_POINT_RADIUS = 0.40
TUBE_RADIUS = 0.12

COLOR_CORTICAL = (0.12, 0.35, 0.95)   # azul
COLOR_JUXTA = (1.00, 0.40, 0.05)      # naranja/rojo intenso
TUBE_DIM = 0.55                       # tubulos un poco mas tenues que glomerulos

HILIO = (0.0, -32.0, 0.0)
HILIO_RADIO = 2.5

# --- CORTE (rebanada delgada) -----------------------------------------------
CORTE_ACTIVO = True
CORTE_EJE = 'Y'            # 'X' (axial) | 'Y' (sagital) | 'Z' (coronal)
CORTE_CENTRO = 0.0        # posicion del plano de corte sobre el eje [mm]
CORTE_GROSOR = 6.0        # grosor de la rebanada [mm]; |coord - centro| < grosor/2
CORTE_SUBSAMPLE = 4      # en corte: dominio mas denso para ver la seccion
CORTE_DOM_EMISSION = 0.6  # en corte: dominio mas visible (estratificacion)

# Color del dominio por region EN EL CORTE (ver estratificacion)
COLOR_CORTEX_SLICE = (0.78, 0.80, 0.88)    # corteza: gris-azulado claro
COLOR_MEDULLA_SLICE = (0.45, 0.40, 0.42)   # medula: gris medio oscuro
COLOR_PIRAMIDE_SLICE = (0.62, 0.45, 0.30)  # piramides: tono calido

_EJE_IDX = {'X': 0, 'Y': 1, 'Z': 2}
_EJE_VIEW = {'X': 'RIGHT', 'Y': 'FRONT', 'Z': 'TOP'}  # vista perpendicular

DOM_NPZ = "capa0_dominio.npz"
NEF_NPZ = "capa1_nefronas.npz"


# ============================================================================
#  UTILIDADES (misma logica robusta que capa0_visualizar.py)
# ============================================================================
def find_npz(name):
    """Localiza un .npz priorizando la raiz del repo (funciona desde el editor
    de Scripting de Blender, donde __file__ no esta definido)."""
    here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() \
        else os.getcwd()
    candidates = [
        os.path.expanduser(os.path.join("~", "Escritorio", "BioKidney-AI", name)),
        os.path.join(here, name),
        os.path.join(here, "..", name),
        os.path.join(os.getcwd(), name),
    ]
    for c in candidates:
        if os.path.isfile(c):
            path = os.path.abspath(c)
            print(f"  [find_npz] {name} -> {path}")
            return path
    raise FileNotFoundError(f"No se encontro {name}. Buscado en: {candidates}")


def clear_scene():
    """Borra todos los objetos de malla existentes y mallas huerfanas."""
    bpy.ops.object.select_all(action='DESELECT')
    for obj in list(bpy.data.objects):
        if obj.type == 'MESH':
            obj.select_set(True)
    bpy.ops.object.delete()
    for m in list(bpy.data.meshes):
        if m.users == 0:
            bpy.data.meshes.remove(m)


def make_emission_material(name, color_attr_name, strength=1.0):
    """Material Emission cuyo color = atributo de color de vertice."""
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    attr = nt.nodes.new("ShaderNodeAttribute")
    attr.attribute_type = 'GEOMETRY'
    attr.attribute_name = color_attr_name
    attr.location = (-400, 0)
    emi = nt.nodes.new("ShaderNodeEmission")
    emi.location = (-150, 0)
    emi.inputs["Strength"].default_value = strength
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    out.location = (100, 0)
    nt.links.new(attr.outputs["Color"], emi.inputs["Color"])
    nt.links.new(emi.outputs["Emission"], out.inputs["Surface"])
    return mat


def make_solid_material(name, color, strength=1.0):
    """Material Emission de color plano (no afectado por sombras)."""
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    emi = nt.nodes.new("ShaderNodeEmission")
    emi.inputs["Color"].default_value = (*color, 1.0)
    emi.inputs["Strength"].default_value = strength
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    out.location = (200, 0)
    nt.links.new(emi.outputs["Emission"], out.inputs["Surface"])
    return mat


def _gn_add_io(ng, name, socket_type):
    """Sockets de E/S del node group (compatible Blender 3.x / 4.x / 5.x)."""
    if hasattr(ng, "interface"):
        ng.interface.new_socket(name, in_out='INPUT', socket_type=socket_type)
        ng.interface.new_socket(name, in_out='OUTPUT', socket_type=socket_type)
    else:
        ng.inputs.new(socket_type, name)
        ng.outputs.new(socket_type, name)


def make_mesh_object(name, verts, edges, rgba, color_attr_name):
    """Crea un objeto de malla (verts + edges opcionales) con atributo de color
    de vertice (FLOAT_COLOR, dominio POINT)."""
    mesh = bpy.data.meshes.new(name + "_mesh")
    mesh.from_pydata(verts.tolist(), edges, [])
    mesh.update()
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.collection.objects.link(obj)
    col = mesh.color_attributes.new(name=color_attr_name,
                                    type='FLOAT_COLOR', domain='POINT')
    col.data.foreach_set("color", rgba.reshape(-1))
    try:
        mesh.color_attributes.active_color = col
    except Exception:
        pass
    mesh.update()
    return obj


def add_points_geonodes(obj, material, radius):
    """Geometry Nodes: Mesh -> Points (vertices) + Set Material."""
    ng = bpy.data.node_groups.new(obj.name + "_GNpts", 'GeometryNodeTree')
    _gn_add_io(ng, "Geometry", 'NodeSocketGeometry')
    nodes, links = ng.nodes, ng.links
    n_in = nodes.new("NodeGroupInput"); n_in.location = (-400, 0)
    n_out = nodes.new("NodeGroupOutput"); n_out.location = (400, 0)
    m2p = nodes.new("GeometryNodeMeshToPoints"); m2p.location = (-150, 0)
    m2p.mode = 'VERTICES'
    m2p.inputs["Radius"].default_value = radius
    setmat = nodes.new("GeometryNodeSetMaterial"); setmat.location = (150, 0)
    setmat.inputs["Material"].default_value = material
    links.new(n_in.outputs[0], m2p.inputs["Mesh"])
    links.new(m2p.outputs["Points"], setmat.inputs["Geometry"])
    links.new(setmat.outputs["Geometry"], n_out.inputs[0])
    mod = obj.modifiers.new(obj.name + "_GNpts", 'NODES')
    mod.node_group = ng
    return mod


def add_tube_geonodes(obj, material, radius):
    """Geometry Nodes: Mesh(edges) -> Curve -> Mesh (tubos) + Set Material.

    Convierte las aristas en curvas y les da grosor con un perfil circular,
    para que los tubulos se vean como tubos coloreados en el viewport
    manteniendo un unico objeto de malla.
    """
    ng = bpy.data.node_groups.new(obj.name + "_GNtube", 'GeometryNodeTree')
    _gn_add_io(ng, "Geometry", 'NodeSocketGeometry')
    nodes, links = ng.nodes, ng.links
    n_in = nodes.new("NodeGroupInput"); n_in.location = (-600, 0)
    n_out = nodes.new("NodeGroupOutput"); n_out.location = (400, 0)

    m2c = nodes.new("GeometryNodeMeshToCurve"); m2c.location = (-400, 100)

    circle = nodes.new("GeometryNodeCurvePrimitiveCircle"); circle.location = (-400, -150)
    circle.mode = 'RADIUS'
    circle.inputs["Resolution"].default_value = 6
    circle.inputs["Radius"].default_value = radius

    c2m = nodes.new("GeometryNodeCurveToMesh"); c2m.location = (-120, 0)

    setmat = nodes.new("GeometryNodeSetMaterial"); setmat.location = (150, 0)
    setmat.inputs["Material"].default_value = material

    links.new(n_in.outputs[0], m2c.inputs["Mesh"])
    links.new(m2c.outputs["Curve"], c2m.inputs["Curve"])
    links.new(circle.outputs["Curve"], c2m.inputs["Profile Curve"])
    links.new(c2m.outputs["Mesh"], setmat.inputs["Geometry"])
    links.new(setmat.outputs["Geometry"], n_out.inputs[0])
    mod = obj.modifiers.new(obj.name + "_GNtube", 'NODES')
    mod.node_group = ng
    return mod


def setup_viewport(view_axis=None):
    """Viewport en modo Material, ejes visibles, encuadre (guarda si headless).
    Si `view_axis` (p.ej. 'FRONT'/'RIGHT'/'TOP') se da, orienta la camara
    perpendicular al plano correspondiente."""
    try:
        for area in bpy.context.screen.areas:
            if area.type == 'VIEW_3D':
                for space in area.spaces:
                    if space.type == 'VIEW_3D':
                        space.shading.type = 'MATERIAL'
                        space.overlay.show_axis_x = True
                        space.overlay.show_axis_y = True
                        space.overlay.show_axis_z = True
                for region in area.regions:
                    if region.type == 'WINDOW':
                        try:
                            with bpy.context.temp_override(area=area, region=region):
                                if view_axis:
                                    bpy.ops.view3d.view_axis(type=view_axis)
                                bpy.ops.view3d.view_all()
                        except Exception:
                            pass
    except Exception:
        pass


def slice_mask(coords, eje, centro, grosor):
    """Mascara de puntos dentro de la rebanada |coord_eje - centro| < grosor/2."""
    return np.abs(coords[:, _EJE_IDX[eje]] - centro) < (grosor * 0.5)


def domain_region_colors(region_label):
    """rgba (n,4) del dominio por region: cortex / medulla / piramide."""
    n = len(region_label)
    rgba = np.empty((n, 4), dtype=np.float32)
    rgba[:, 3] = 1.0
    is_cortex = region_label == "cortex"
    is_pir = np.char.startswith(region_label, "piramide")
    is_med = ~is_cortex & ~is_pir
    rgba[is_cortex, :3] = COLOR_CORTEX_SLICE
    rgba[is_med, :3] = COLOR_MEDULLA_SLICE
    rgba[is_pir, :3] = COLOR_PIRAMIDE_SLICE
    return rgba


def colors_by_tipo(tipo, dim=1.0):
    """Devuelve rgba (n,4) segun tipo (cortical=azul, yuxtamedular=naranja)."""
    n = len(tipo)
    rgba = np.empty((n, 4), dtype=np.float32)
    rgba[:, 3] = 1.0
    is_juxta = tipo == "yuxtamedular"
    rgba[~is_juxta, :3] = np.array(COLOR_CORTICAL) * dim
    rgba[is_juxta, :3] = np.array(COLOR_JUXTA) * dim
    return rgba


# ============================================================================
#  MAIN
# ============================================================================
def main():
    print("=" * 70)
    print("  CAPA 1 - VISUALIZACION DE NEFRONAS (sobre dominio Capa 0)")
    print("=" * 70)

    dom_path = find_npz(DOM_NPZ)
    nef_path = find_npz(NEF_NPZ)
    dom = np.load(dom_path, allow_pickle=False)
    nef = np.load(nef_path, allow_pickle=False)
    print("  [OK] cargados ambos .npz (dominio + nefronas)")

    clear_scene()

    if CORTE_ACTIVO:
        print(f"  MODO CORTE: eje {CORTE_EJE}, centro {CORTE_CENTRO} mm, "
              f"grosor {CORTE_GROSOR} mm")
    dom_subsample = CORTE_SUBSAMPLE if CORTE_ACTIVO else DOM_SUBSAMPLE
    dom_emission = CORTE_DOM_EMISSION if CORTE_ACTIVO else DOM_EMISSION

    # ---------------------------------------------------------------- dominio
    dcoords = dom["coords"][::dom_subsample]
    dlabels = dom["region_label"][::dom_subsample]
    if CORTE_ACTIVO:
        m = slice_mask(dcoords, CORTE_EJE, CORTE_CENTRO, CORTE_GROSOR)
        dcoords, dlabels = dcoords[m], dlabels[m]
        rgba_dom = domain_region_colors(dlabels)       # color por region
    else:
        rgba_dom = np.empty((len(dcoords), 4), dtype=np.float32)
        rgba_dom[:, :3] = DOM_COLOR
        rgba_dom[:, 3] = 1.0
    nd = len(dcoords)
    obj_dom = make_mesh_object("capa0_dominio_bg", dcoords, [], rgba_dom, "col_dom")
    mat_dom = make_emission_material("mat_dom_bg", "col_dom", strength=dom_emission)
    add_points_geonodes(obj_dom, mat_dom, DOM_POINT_RADIUS)

    # ------------------------------------------------------------------ hilio
    bpy.ops.mesh.primitive_uv_sphere_add(radius=HILIO_RADIO, location=HILIO)
    hobj = bpy.context.active_object
    hobj.name = "hilio_seno_renal"
    hobj.data.materials.append(make_solid_material("mat_hilio", (1.0, 0.0, 0.0),
                                                    strength=3.0))

    # --- filtro de nefronas: glomerulo dentro de la rebanada ----------------
    glom = nef["glomerulos"].astype(np.float64)
    tipo = nef["tipo"]
    tub = nef["tubulos"].astype(np.float64)   # (N, NPTS, 3)
    if CORTE_ACTIVO:
        gm = slice_mask(glom, CORTE_EJE, CORTE_CENTRO, CORTE_GROSOR)
        glom, tipo, tub = glom[gm], tipo[gm], tub[gm]

    # ------------------------------------------------------------- glomerulos
    ng_ = len(glom)
    rgba_glom = colors_by_tipo(tipo, dim=1.0)
    obj_glom = make_mesh_object("capa1_glomerulos", glom, [], rgba_glom, "col_glom")
    mat_glom = make_emission_material("mat_glom", "col_glom", strength=1.6)
    add_points_geonodes(obj_glom, mat_glom, GLOM_POINT_RADIUS)

    # ---------------------------------------------------------------- tubulos
    n_tub, npts, _ = tub.shape
    verts = tub.reshape(-1, 3)                # (n_tub*NPTS, 3)
    # aristas consecutivas dentro de cada tubulo
    base = (np.arange(n_tub) * npts)[:, None]
    seg = np.arange(npts - 1)[None, :]
    a = (base + seg).reshape(-1)
    edges = np.stack([a, a + 1], axis=1).tolist()
    # color por vertice = color del tipo de su nefrona (tenue)
    rgba_tub = np.repeat(colors_by_tipo(tipo, dim=TUBE_DIM), npts, axis=0)
    obj_tub = make_mesh_object("capa1_tubulos", verts, edges, rgba_tub, "col_tub")
    mat_tub = make_emission_material("mat_tub", "col_tub", strength=1.0)
    add_tube_geonodes(obj_tub, mat_tub, TUBE_RADIUS)

    # ----------------------------------------------------------------- vista
    view = _EJE_VIEW.get(CORTE_EJE) if CORTE_ACTIVO else None
    setup_viewport(view_axis=view)

    # --------------------------------------------------------------- resumen
    n_cort = int(np.count_nonzero(tipo == "cortical"))
    n_jux = int(np.count_nonzero(tipo == "yuxtamedular"))
    print("-" * 70)
    if CORTE_ACTIVO:
        print(f"  Dominio (rebanada): {nd} puntos (SUBSAMPLE={dom_subsample}) "
              f"coloreado por region")
        print(f"  NEFRONAS DENTRO DEL CORTE: {ng_}  de {len(nef['tipo'])} totales")
    else:
        print(f"  Dominio de fondo (Capa 0): {nd} puntos (SUBSAMPLE={dom_subsample})")
        print(f"  Glomerulos dibujados: {ng_}")
    print(f"     corticales (azul)        : {n_cort}")
    print(f"     yuxtamedulares (naranja) : {n_jux}")
    print(f"  Tubulos dibujados: {n_tub}  ({npts} vertices c/u)")
    print("  Hilio / seno renal en:", HILIO)
    print(f"  Vista: {'perpendicular ' + str(view) if view else 'libre'}  |  "
          f"Eje X = superior-inferior")
    print("=" * 70)


if __name__ == "__main__":
    main()
