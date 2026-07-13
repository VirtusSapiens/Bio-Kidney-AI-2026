#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa0_visualizar.py  --  Bio-Kidney AI 2026
============================================
Visualizacion del DOMINIO renal (capa 0) DENTRO de Blender.

  >>> Se ejecuta con bpy (Blender), NO en env_biokidney <<<

Lee capa0_dominio.npz y dibuja el dominio renal como una nube de puntos
coloreada por region (corteza / medula / piramides medulares) usando un
unico objeto de malla (solo vertices) + atributo de color de vertice,
renderizado como puntos mediante Geometry Nodes (Mesh to Points).

Uso:
    blender --python 08_gemelo_digital/capa0_visualizar.py
    # o desde el editor de scripting de Blender: abrir y "Run Script"
"""

import os
import bpy
import numpy as np

# ============================================================================
#  PARAMETROS
# ============================================================================
SUBSAMPLE = 8          # dibujar 1 de cada SUBSAMPLE puntos (~25000 con N=200k)
POINT_RADIUS = 0.35    # radio de los puntos en viewport [mm]
# Hilio / seno renal: desplazado ligeramente hacia afuera (-Y) para que
# sobresalga de la capsula (la superficie medial esta en y=-30).
HILIO = (0.0, -32.0, 0.0)
HILIO_RADIO = 2.5                  # radio del marcador esferico [mm]
# Flecha medial->lateral: del hilio hacia el centro del organo
HILIO_FLECHA_FIN = (0.0, -15.0, 0.0)

# Nombre del archivo de dominio (se busca junto al script y en el padre)
NPZ_NAME = "capa0_dominio.npz"

# Colores de region (RGB 0-1)
COLOR_CORTEX = (0.80, 0.80, 0.80)      # gris claro
COLOR_MEDULLA = (0.42, 0.42, 0.42)     # gris medio (columnas de Bertin)

# Paleta cualitativa de 10 colores bien diferenciados (tipo tab10, sin grises)
PALETTE10 = np.array([
    (31, 119, 180),    # azul
    (255, 127, 14),    # naranja
    (44, 160, 44),     # verde
    (214, 39, 40),     # rojo
    (148, 103, 189),   # violeta
    (140, 86, 75),     # marron
    (227, 119, 194),   # rosa
    (188, 189, 34),    # oliva
    (23, 190, 207),    # cian
    (255, 215, 0),     # dorado
], dtype=np.float64) / 255.0


# ============================================================================
#  UTILIDADES
# ============================================================================
def find_npz():
    """Localiza capa0_dominio.npz.

    Prioriza la ruta absoluta de la raiz del repo (~/Escritorio/BioKidney-AI),
    para que funcione tambien desde el editor de Scripting de Blender, donde
    __file__ no esta definido y las rutas relativas al script fallan. Las
    rutas relativas quedan como respaldo.
    """
    here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() \
        else os.getcwd()
    candidates = [
        os.path.expanduser(os.path.join("~", "Escritorio", "BioKidney-AI", NPZ_NAME)),
        os.path.join(here, NPZ_NAME),
        os.path.join(here, "..", NPZ_NAME),
        os.path.join(os.getcwd(), NPZ_NAME),
    ]
    for c in candidates:
        if os.path.isfile(c):
            path = os.path.abspath(c)
            print("  [find_npz] usando:", path)
            return path
    raise FileNotFoundError(
        f"No se encontro {NPZ_NAME}. Buscado en: {candidates}")


def clear_scene():
    """Borra todos los objetos de malla existentes."""
    bpy.ops.object.select_all(action='DESELECT')
    for obj in list(bpy.data.objects):
        if obj.type == 'MESH':
            obj.select_set(True)
    bpy.ops.object.delete()
    # limpiar mallas huerfanas
    for m in list(bpy.data.meshes):
        if m.users == 0:
            bpy.data.meshes.remove(m)


def make_emission_material(name, color_attr_name):
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
    emi.inputs["Strength"].default_value = 1.0

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
    """Anade sockets de E/S al node group (compatible Blender 3.x / 4.x)."""
    if hasattr(ng, "interface"):   # Blender 4.x
        ng.interface.new_socket(name, in_out='INPUT', socket_type=socket_type)
        ng.interface.new_socket(name, in_out='OUTPUT', socket_type=socket_type)
    else:                           # Blender 3.x
        ng.inputs.new(socket_type, name)
        ng.outputs.new(socket_type, name)


def add_points_geonodes(obj, material, radius):
    """Modificador Geometry Nodes: Mesh -> Points (vertices) + Set Material.

    Permite ver los vertices como puntos coloreados en viewport/render
    manteniendo un unico objeto de malla.
    """
    ng = bpy.data.node_groups.new("GN_capa0_puntos", 'GeometryNodeTree')
    _gn_add_io(ng, "Geometry", 'NodeSocketGeometry')

    nodes, links = ng.nodes, ng.links
    n_in = nodes.new("NodeGroupInput")
    n_in.location = (-400, 0)
    n_out = nodes.new("NodeGroupOutput")
    n_out.location = (400, 0)

    m2p = nodes.new("GeometryNodeMeshToPoints")
    m2p.location = (-150, 0)
    m2p.mode = 'VERTICES'
    m2p.inputs["Radius"].default_value = radius

    setmat = nodes.new("GeometryNodeSetMaterial")
    setmat.location = (150, 0)
    setmat.inputs["Material"].default_value = material

    links.new(n_in.outputs[0], m2p.inputs["Mesh"])
    links.new(m2p.outputs["Points"], setmat.inputs["Geometry"])
    links.new(setmat.outputs["Geometry"], n_out.inputs[0])

    mod = obj.modifiers.new("GN_capa0_puntos", 'NODES')
    mod.node_group = ng
    return mod


def setup_viewport():
    """Pone el viewport en modo Material y encuadra la escena (guarda si headless)."""
    try:
        for area in bpy.context.screen.areas:
            if area.type == 'VIEW_3D':
                for space in area.spaces:
                    if space.type == 'VIEW_3D':
                        space.shading.type = 'MATERIAL'
                        space.overlay.show_axis_x = True
                        space.overlay.show_axis_y = True
                        space.overlay.show_axis_z = True
                # encuadrar todo
                for region in area.regions:
                    if region.type == 'WINDOW':
                        override = {'area': area, 'region': region}
                        try:
                            with bpy.context.temp_override(**override):
                                bpy.ops.view3d.view_all()
                        except Exception:
                            pass
    except Exception:
        pass


# ============================================================================
#  MAIN
# ============================================================================
def main():
    npz_path = find_npz()
    print("=" * 70)
    print("  CAPA 0 - VISUALIZACION DEL DOMINIO RENAL (Blender)")
    print("=" * 70)
    print("  Leyendo:", npz_path)

    data = np.load(npz_path, allow_pickle=False)
    coords = data["coords"]
    labels = data["region_label"]
    # depth = data["depth"]  # disponible si se quiere colorear por profundidad

    # --- submuestreo ---
    coords = coords[::SUBSAMPLE]
    labels = labels[::SUBSAMPLE]
    n = len(coords)
    print(f"  SUBSAMPLE = {SUBSAMPLE}  ->  {n} puntos dibujados")

    # --- color por region (vectorizado) ---
    rgba = np.empty((n, 4), dtype=np.float32)
    rgba[:, 3] = 1.0

    counts = {}
    cortex_mask = labels == "cortex"
    medulla_mask = labels == "medulla"
    rgba[cortex_mask, :3] = COLOR_CORTEX
    rgba[medulla_mask, :3] = COLOR_MEDULLA
    counts["cortex"] = int(cortex_mask.sum())
    counts["medulla"] = int(medulla_mask.sum())

    n_pir = len(PALETTE10)
    for i in range(n_pir):
        name = f"piramide_{i:02d}"
        mask = labels == name
        c = int(mask.sum())
        if c:
            rgba[mask, :3] = PALETTE10[i % n_pir]
            counts[name] = c

    # --- crear malla (solo vertices) ---
    clear_scene()
    mesh = bpy.data.meshes.new("capa0_dominio_mesh")
    mesh.from_pydata(coords.tolist(), [], [])
    mesh.update()

    obj = bpy.data.objects.new("capa0_dominio", mesh)
    bpy.context.collection.objects.link(obj)

    # --- atributo de color de vertice (dominio POINT) ---
    col = mesh.color_attributes.new(name="region_color",
                                    type='FLOAT_COLOR', domain='POINT')
    col.data.foreach_set("color", rgba.reshape(-1))
    try:
        mesh.color_attributes.active_color = col
    except Exception:
        pass
    mesh.update()

    # --- material + geometry nodes (mesh -> points) ---
    mat = make_emission_material("mat_capa0_dominio", "region_color")
    add_points_geonodes(obj, mat, POINT_RADIUS)

    # --- marcar el hilio / seno renal ---
    # Esfera roja pura e intensa, desplazada hacia afuera para sobresalir.
    mat_hilio = make_solid_material("mat_hilio", (1.0, 0.0, 0.0), strength=3.0)
    bpy.ops.mesh.primitive_uv_sphere_add(radius=HILIO_RADIO, location=HILIO)
    hilio_obj = bpy.context.active_object
    hilio_obj.name = "hilio_seno_renal"
    hilio_obj.data.materials.append(mat_hilio)

    # Flecha medial->lateral: cilindro delgado del hilio hacia el centro.
    import mathutils
    start = mathutils.Vector(HILIO)
    end = mathutils.Vector(HILIO_FLECHA_FIN)
    vec = end - start
    length = vec.length
    bpy.ops.mesh.primitive_cylinder_add(radius=0.5, depth=length,
                                        location=(start + end) * 0.5)
    flecha = bpy.context.active_object
    flecha.name = "hilio_eje_medial_lateral"
    # alinear el eje Z (por defecto) del cilindro con la direccion hilio->centro
    flecha.rotation_mode = 'QUATERNION'
    flecha.rotation_quaternion = mathutils.Vector((0, 0, 1)).rotation_difference(vec)
    flecha.data.materials.append(mat_hilio)

    # --- viewport ---
    setup_viewport()

    # --- resumen ---
    print("-" * 70)
    print("  Conteo por region (puntos dibujados):")
    for k in ["cortex", "medulla"]:
        print(f"     {k:14s}: {counts.get(k, 0):>7d}")
    for i in range(n_pir):
        name = f"piramide_{i:02d}"
        if name in counts:
            print(f"     {name:14s}: {counts[name]:>7d}")
    print(f"  TOTAL puntos dibujados: {n}")
    print(f"  Hilio / seno renal (esfera roja) en: {HILIO}  radio={HILIO_RADIO} mm")
    print(f"  Flecha medial->lateral: {HILIO} -> {HILIO_FLECHA_FIN}")
    print("  Eje X = superior-inferior  (eje largo del rinon)")
    print("=" * 70)


if __name__ == "__main__":
    main()
