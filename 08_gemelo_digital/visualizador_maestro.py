#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
visualizador_maestro.py  --  Bio-Kidney AI 2026
================================================
Visualizador UNIFICADO del gemelo digital renal. Reemplaza la logica de
capa0_visualizar.py + capa1_visualizar.py en un solo script con CAPAS
ACTIVABLES por flags booleanos.

  >>> Se ejecuta con bpy (Blender), NO en env_biokidney <<<

- Carga PEREZOSA: intenta cargar cada .npz; si falta, desactiva su capa con
  un aviso y continua (sirve aunque falten capas futuras).
- Una funcion por capa (draw_*), cada una respeta su flag y el CORTE.
- Patron probado: objeto unico de malla + atributo de color + Geometry Nodes
  (Mesh to Points para puntos, Curve to Mesh para tubulos).
- Cada capa va a su propia coleccion de Blender (visible/ocultable en outliner).

Uso:
    blender --python 08_gemelo_digital/visualizador_maestro.py
"""

import os
import bpy
import numpy as np

# ============================================================================
#  FLAGS DE CONTROL  (editar aqui)
# ============================================================================
VER_DOMINIO    = False   # Capa 0: nube del parenquima (corteza + medula)
VER_PIRAMIDES  = False   # Capa 0: piramides medulares coloreadas
VER_HILIO      = True    # Capa 0: esfera + flecha del hilio/seno
VER_GLOMERULOS = False   # Capa 1: glomerulos azul/naranja
VER_TUBULOS    = False   # Capa 1: tubulos
VER_DEMANDA    = False   # Capa 2: mapa de calor de demanda metabolica de O2
VER_ARTERIAL   = True    # Capa 3a: arbol arterial (space colonization + Murray)
VER_VENOSO     = True    # Capa 3b: arbol venoso (reparado, canonico) -> azul
VER_COLECTOR   = True    # Capa 3c: arbol colector urinario (bosque 10 subarboles) -> amarillo
VER_PERITUBULAR = True   # Capa 3a-bis: red peritubular + vasa recta (puntos de drenaje)
VER_CALICES    = True    # Capa 4: sistema calicial alto (papila->caliz->pelvis->ureter)

# DEBUG: colorear el colector por piramide_id (10 tonos) en vez de amarillo plano,
# para verificar de un vistazo que son 10 territorios segregados y no una maraña.
COLECTOR_POR_PIRAMIDE = False
# INSPECCION CAPA 3a-bis: peritubular + arterial ON (para ver la relacion entre la
# segunda red capilar y el territorio arterial que la alimenta); dominio/glomerulos/
# demanda OFF para leer limpio como la VASA RECTA violeta se hunde en la medula
# mientras la red peritubular verde tapiza la corteza. El corte sagital
# (CORTE_ACTIVO=True, CORTE_EJE='Y') es especialmente revelador aqui.
# NOTA: al inspeccionar la demanda (VER_DEMANDA=True) conviene poner
# VER_DOMINIO=False, para que la nube gris plana del dominio no tape el mapa
# de calor (ambos pintan los mismos puntos del parenquima). No es obligatorio,
# pero el mapa de calor se lee mucho mejor sin el dominio gris encima.

CORTE_ACTIVO = False     # corte sagital (rebanada delgada)
CORTE_EJE = 'Y'          # 'X' (axial) | 'Y' (sagital) | 'Z' (coronal)
CORTE_CENTRO = 0.0       # posicion del plano [mm]
CORTE_GROSOR = 6.0       # grosor de la rebanada [mm]

# ============================================================================
#  PARAMETROS DE DIBUJO
# ============================================================================
DOM_SUBSAMPLE = 12       # nube de dominio (1 de cada N)
DOM_SUBSAMPLE_CORTE = 4  # mas denso en corte para ver la seccion
DOM_POINT_RADIUS = 0.25
PIR_POINT_RADIUS = 0.30
GLOM_POINT_RADIUS = 0.40
TUBE_RADIUS = 0.12
SUBSAMPLE_DEMANDA = 4    # mapa de calor de demanda (1 de cada N)
DEMANDA_POINT_RADIUS = 0.30
DEMANDA_EMISSION = 1.3   # mapa de calor algo emisivo para "brillar"

DOM_EMISSION = 0.25      # tenue en vista completa
DOM_EMISSION_CORTE = 0.60
PIR_EMISSION = 0.8
GLOM_EMISSION = 1.6
TUBE_EMISSION = 1.0
TUBE_DIM = 0.55          # tubulos algo mas tenues que glomerulos

# Capa 3a - arbol arterial
# Los radios reales de Murray van de ~0.012 mm (terminal) a ~0.53 mm (raiz):
# invisibles a escala real. Se remapean al rango VISIBLE [MIN, MAX] conservando
# el orden, para LEER el calibre (tronco grueso -> terminales finas) sin que las
# arteriolas desaparezcan. El reporte de consola da los radios REALES en mm.
ART_TUBE_MIN = 0.06      # radio de tubo visible para la arteriola terminal [mm display]
ART_TUBE_MAX = 0.70      # radio de tubo visible para el tronco/raiz   [mm display]
ART_PROFILE_RES = 6      # resolucion del perfil circular del tubo
ART_ROOT_RADIUS = 2.0    # esfera marcadora de la RAIZ (hilio)
ART_TERMINAL_POINT_RADIUS = 0.16   # puntos de los nodos TERMINALES (arteriolas)
ART_EMISSION = 1.2

# Capa 3b - arbol venoso (mismo patron que el arterial; remapeo de calibre analogo)
VEN_TUBE_MIN = 0.06      # radio de tubo visible para la venula terminal [mm display]
VEN_TUBE_MAX = 0.70      # radio de tubo visible para el colector/raiz    [mm display]
VEN_ROOT_RADIUS = 2.0    # esfera marcadora de la RAIZ (colector en el hilio)
VEN_TERMINAL_POINT_RADIUS = 0.16   # puntos de los nodos TERMINALES (venulas)
VEN_EMISSION = 1.2

# Capa 3c - arbol colector urinario (grafo nodos+parent+radios, mismo patron de tubo).
# RECONCILIACION DE UNIDADES: 'radios' del .npz esta en mm (0.015-0.125 mm = 15-125 um),
# la MISMA unidad que arterial/venoso (coords en mm). A diferencia del arterial/venoso
# (que remapean min-max a un rango display), aqui el radio del tubo = radios * escala:
# se conserva el TAMAÑO FISICO relativo (el colector es fino de verdad), y el factor
# es SOLO-VISUAL para que no quede invisible junto a los grandes vasos.
RADIO_ESCALA_COLECTOR = 4.0   # factor COSMETICO (no altera el .npz): 15-125um -> tubo
                              # 0.06-0.50 mm (0.06 = mismo suelo visible que el arterial).
                              # Subir si el colector queda invisible junto al arterial.
COLECTOR_ROOT_RADIUS = 1.2    # esfera marcadora de cada PAPILA (apex, 10 raices)
COLECTOR_TERMINAL_POINT_RADIUS = 0.12   # puntos de los nodos TERMINALES (bocas colectoras)
COLECTOR_EMISSION = 1.2

# Capa 4 - sistema calicial alto (centerline + radio; radios ABSOLUTOS, sin escala cosmetica:
# el gate es ver si la pelvis-como-cilindro basta a tamaño fisico real 0.2-4.5 mm)
CALIZ_CONE_VERTS = 20         # segmentos del cono tapered (cosmetico, no altera el .npz)
CALIZ_BELLINI_FRAC = 0.35     # fraccion del segmento Bellini->copa que es el stub FINO r=0.2
                              # (arista_tipo==1): stub fino + copa r=1.5 tope-a-tope = ESCALON visible
CALIZ_NODE_MARKER = 0.30      # radio de la esferita marcadora por nodo (color por nivel)

# Capa 3a-bis - red peritubular + vasa recta (puntos de drenaje venoso)
PERI_POINT_RADIUS = 0.25   # puntos de drenaje (mismo patron Mesh to Points)
PERI_EMISSION = 1.4        # algo emisivo para que la vasa recta "brille" en la medula

# Colores
COLOR_CORTEX = (0.78, 0.80, 0.88)     # corteza gris-azulado claro
COLOR_MEDULLA = (0.45, 0.40, 0.42)    # medula gris medio
COLOR_CORTICAL = (0.12, 0.35, 0.95)   # glomerulo cortical -> azul
COLOR_JUXTA = (1.00, 0.40, 0.05)      # glomerulo yuxtamedular -> naranja
COLOR_ARTERIAL = (0.85, 0.06, 0.06)   # arbol arterial -> rojo
COLOR_VENOSO = (0.08, 0.22, 0.92)     # arbol venoso -> azul
COLOR_VENOSO_TERMINAL = (0.35, 0.65, 1.00)   # nodos terminales (venulas) -> azul claro
COLOR_COLECTOR = (0.95, 0.82, 0.05)          # arbol colector urinario -> amarillo
COLOR_COLECTOR_TERMINAL = (1.00, 0.92, 0.45)   # bocas colectoras -> amarillo claro
COLOR_TERMINAL = (1.00, 0.55, 0.35)   # nodos terminales (arteriolas) -> salmon claro
COLOR_PERITUBULAR = (0.10, 0.85, 0.25)  # red peritubular cortical -> verde
COLOR_VASA_RECTA = (0.80, 0.10, 0.90)   # vasa recta (medula) -> magenta/violeta intenso
# Capa 4: color por 'nivel' (enum 0..4)
COLOR_CALIZ = [
    (0.55, 0.55, 0.55),   # 0 papila_junction -> gris
    (1.00, 0.75, 0.45),   # 1 caliz_menor     -> naranja claro
    (1.00, 0.55, 0.10),   # 2 caliz_mayor     -> naranja
    (0.90, 0.10, 0.10),   # 3 pelvis          -> rojo
    (0.90, 0.10, 0.90),   # 4 ureter          -> magenta
]
PALETTE10 = np.array([
    (31, 119, 180), (255, 127, 14), (44, 160, 44), (214, 39, 40),
    (148, 103, 189), (140, 86, 75), (227, 119, 194), (188, 189, 34),
    (23, 190, 207), (255, 215, 0),
], dtype=np.float64) / 255.0

# Hilio
HILIO = (0.0, -32.0, 0.0)
HILIO_RADIO = 2.5
HILIO_FLECHA_FIN = (0.0, -15.0, 0.0)

# Archivos
DOM_NPZ = "capa0_dominio.npz"
NEF_NPZ = "capa1_nefronas.npz"
DEM_NPZ = "capa2_demanda.npz"
ART_NPZ = "capa3a_arterial.npz"
VEN_NPZ = "capa3b_venoso.npz"
COLECTOR_NPZ = "capa3c_colector.npz"
PERI_NPZ = "capa3ab_peritubular.npz"
CALICES_NPZ = "capa4_colector_alto.npz"

# Colecciones
COL_DOMINIO = "Coleccion_Dominio"
COL_NEFRONAS = "Coleccion_Nefronas"
COL_DEMANDA = "Coleccion_Demanda"
COL_ARTERIAL = "Coleccion_Arterial"
COL_VENOSO = "Coleccion_Venoso"
COL_COLECTOR = "Coleccion_Colector"
COL_PERITUBULAR = "Coleccion_Peritubular"
COL_CALICES = "Coleccion_Calices"

_EJE_IDX = {'X': 0, 'Y': 1, 'Z': 2}
_EJE_VIEW = {'X': 'RIGHT', 'Y': 'FRONT', 'Z': 'TOP'}


# ============================================================================
#  UTILIDADES DE ARCHIVO / ESCENA
# ============================================================================
def find_npz(name):
    """Localiza un .npz priorizando la raiz del repo. Devuelve ruta o None."""
    here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() \
        else os.getcwd()
    for c in (
        os.path.expanduser(os.path.join("~", "Escritorio", "BioKidney-AI", name)),
        os.path.join(here, name),
        os.path.join(here, "..", name),
        os.path.join(os.getcwd(), name),
    ):
        if os.path.isfile(c):
            return os.path.abspath(c)
    return None


def try_load(name):
    """Carga perezosa: devuelve (data|None, ruta|None)."""
    path = find_npz(name)
    if path is None:
        return None, None
    return np.load(path, allow_pickle=False), path


def clear_scene():
    """Borra objetos de malla, mallas huerfanas y colecciones de capa previas.

    Borra directamente desde bpy.data (NO usa bpy.ops.object.select_all /
    .delete, que dependen del View Layer y fallan al re-ejecutar el script con
    'Object X cannot be selected because it is not in View Layer'). Cada borrado
    va en su propio try/except para que un objeto problematico no aborte todo.
    """
    # 1) objetos de malla (cubre tambien hilio/flecha y nubes de puntos)
    for obj in list(bpy.data.objects):
        if obj.type == 'MESH':
            try:
                bpy.data.objects.remove(obj, do_unlink=True)
            except Exception as e:
                print(f"  [clear_scene] no se pudo borrar objeto {obj.name!r}: {e}")

    # 2) mallas huerfanas (sin usuarios)
    for m in list(bpy.data.meshes):
        if m.users == 0:
            try:
                bpy.data.meshes.remove(m)
            except Exception:
                pass

    # 3) colecciones de capa creadas por el script (no acumular vacias)
    for cname in (COL_DOMINIO, COL_NEFRONAS, COL_DEMANDA, COL_ARTERIAL,
                  COL_VENOSO, COL_COLECTOR, COL_PERITUBULAR):
        coll = bpy.data.collections.get(cname)
        if coll is not None:
            try:
                bpy.data.collections.remove(coll)
            except Exception:
                pass


def get_collection(name):
    """Crea (o recupera) una coleccion enlazada a la escena."""
    if name in bpy.data.collections:
        return bpy.data.collections[name]
    coll = bpy.data.collections.new(name)
    bpy.context.scene.collection.children.link(coll)
    return coll


def link_to_collection(obj, coll):
    for c in list(obj.users_collection):
        c.objects.unlink(obj)
    coll.objects.link(obj)


# ============================================================================
#  MATERIALES Y GEOMETRY NODES (patron probado)
# ============================================================================
def make_emission_material(name, color_attr_name, strength=1.0):
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
    if hasattr(ng, "interface"):
        ng.interface.new_socket(name, in_out='INPUT', socket_type=socket_type)
        ng.interface.new_socket(name, in_out='OUTPUT', socket_type=socket_type)
    else:
        ng.inputs.new(socket_type, name)
        ng.outputs.new(socket_type, name)


def make_mesh_object(name, verts, edges, rgba, color_attr_name, coll):
    mesh = bpy.data.meshes.new(name + "_mesh")
    mesh.from_pydata(verts.tolist(), edges, [])
    mesh.update()
    obj = bpy.data.objects.new(name, mesh)
    coll.objects.link(obj)
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


def add_point_float_attr(obj, name, values):
    """Anade un atributo FLOAT por punto (p.ej. radio por nodo para modular tubo)."""
    mesh = obj.data
    attr = mesh.attributes.new(name=name, type='FLOAT', domain='POINT')
    attr.data.foreach_set("value", np.asarray(values, dtype=np.float32).reshape(-1))
    mesh.update()
    return attr


def add_tree_tube_geonodes(obj, material, radius_attr="radius",
                           profile_res=ART_PROFILE_RES):
    """Curve to Mesh con RADIO MODULADO por un atributo de punto.

    Igual que add_tube_geonodes pero el radio del tubo se toma del atributo
    'radius_attr' (por nodo) via Set Curve Radius: el perfil circular es de radio
    1.0 y la curva lo escala punto a punto -> tronco grueso, terminales finas.
    """
    ng = bpy.data.node_groups.new(obj.name + "_GNtree", 'GeometryNodeTree')
    _gn_add_io(ng, "Geometry", 'NodeSocketGeometry')
    nodes, links = ng.nodes, ng.links
    n_in = nodes.new("NodeGroupInput"); n_in.location = (-800, 0)
    n_out = nodes.new("NodeGroupOutput"); n_out.location = (600, 0)
    m2c = nodes.new("GeometryNodeMeshToCurve"); m2c.location = (-600, 120)
    na = nodes.new("GeometryNodeInputNamedAttribute"); na.location = (-600, -220)
    na.data_type = 'FLOAT'
    try:
        na.inputs["Name"].default_value = radius_attr
    except Exception:
        pass
    setr = nodes.new("GeometryNodeSetCurveRadius"); setr.location = (-360, 60)
    circle = nodes.new("GeometryNodeCurvePrimitiveCircle"); circle.location = (-360, -220)
    circle.mode = 'RADIUS'
    circle.inputs["Resolution"].default_value = profile_res
    circle.inputs["Radius"].default_value = 1.0
    c2m = nodes.new("GeometryNodeCurveToMesh"); c2m.location = (-100, 0)
    setmat = nodes.new("GeometryNodeSetMaterial"); setmat.location = (220, 0)
    setmat.inputs["Material"].default_value = material
    links.new(n_in.outputs[0], m2c.inputs["Mesh"])
    links.new(m2c.outputs["Curve"], setr.inputs["Curve"])
    links.new(na.outputs["Attribute"], setr.inputs["Radius"])
    links.new(setr.outputs["Curve"], c2m.inputs["Curve"])
    links.new(circle.outputs["Curve"], c2m.inputs["Profile Curve"])
    links.new(c2m.outputs["Mesh"], setmat.inputs["Geometry"])
    links.new(setmat.outputs["Geometry"], n_out.inputs[0])
    mod = obj.modifiers.new(obj.name + "_GNtree", 'NODES')
    mod.node_group = ng
    return mod


# ============================================================================
#  HELPERS DE DATOS
# ============================================================================
def slice_mask(coords, eje=None, centro=None, grosor=None):
    """Mascara de la rebanada. Si CORTE inactivo -> todo True."""
    if not CORTE_ACTIVO:
        return np.ones(len(coords), dtype=bool)
    eje = CORTE_EJE if eje is None else eje
    centro = CORTE_CENTRO if centro is None else centro
    grosor = CORTE_GROSOR if grosor is None else grosor
    return np.abs(coords[:, _EJE_IDX[eje]] - centro) < (grosor * 0.5)


def colors_by_tipo(tipo, dim=1.0):
    n = len(tipo)
    rgba = np.empty((n, 4), dtype=np.float32)
    rgba[:, 3] = 1.0
    is_juxta = tipo == "yuxtamedular"
    rgba[~is_juxta, :3] = np.array(COLOR_CORTICAL) * dim
    rgba[is_juxta, :3] = np.array(COLOR_JUXTA) * dim
    return rgba


def dom_subsample():
    return DOM_SUBSAMPLE_CORTE if CORTE_ACTIVO else DOM_SUBSAMPLE


# Colormap tipo inferno (frio->caliente) por interpolacion lineal sobre paradas.
# Sin matplotlib: np.interp canal a canal.
_HEAT_POS = np.array([0.00, 0.35, 0.60, 0.80, 1.00])
_HEAT_RGB = np.array([
    (0.08, 0.02, 0.25),   # demanda baja (medula ~0.3): morado/azul oscuro frio
    (0.55, 0.08, 0.45),   # magenta/purpura
    (0.95, 0.45, 0.10),   # naranja (demanda media)
    (1.00, 0.85, 0.20),   # amarillo
    (1.00, 1.00, 0.95),   # demanda alta (~1.0): blanco caliente
])


def heatmap_colors(values):
    """rgba (n,4) mapeando values en [0,1] a un colormap frio->caliente."""
    v = np.clip(np.asarray(values, dtype=np.float64), 0.0, 1.0)
    rgba = np.empty((len(v), 4), dtype=np.float32)
    rgba[:, 3] = 1.0
    for ch in range(3):
        rgba[:, ch] = np.interp(v, _HEAT_POS, _HEAT_RGB[:, ch])
    return rgba


# ============================================================================
#  FUNCIONES POR CAPA
# ============================================================================
def draw_dominio(dom, coll):
    """Capa 0: nube del parenquima (corteza + medula generica, sin piramides)."""
    sub = dom_subsample()
    coords = dom["coords"][::sub]
    labels = dom["region_label"][::sub]
    keep = ~np.char.startswith(labels, "piramide")
    keep &= slice_mask(coords)
    coords, labels = coords[keep], labels[keep]
    rgba = np.empty((len(coords), 4), dtype=np.float32)
    rgba[:, 3] = 1.0
    is_cortex = labels == "cortex"
    rgba[is_cortex, :3] = COLOR_CORTEX
    rgba[~is_cortex, :3] = COLOR_MEDULLA
    emi = DOM_EMISSION_CORTE if CORTE_ACTIVO else DOM_EMISSION
    obj = make_mesh_object("capa0_dominio", coords, [], rgba, "col_dom", coll)
    add_points_geonodes(obj, make_emission_material("mat_dom", "col_dom", emi),
                        DOM_POINT_RADIUS)
    return len(coords)


def draw_piramides(dom, coll):
    """Capa 0: puntos de las piramides medulares, color por id de piramide."""
    sub = dom_subsample()
    coords = dom["coords"][::sub]
    labels = dom["region_label"][::sub]
    keep = np.char.startswith(labels, "piramide")
    keep &= slice_mask(coords)
    coords, labels = coords[keep], labels[keep]
    rgba = np.empty((len(coords), 4), dtype=np.float32)
    rgba[:, 3] = 1.0
    n_pir = len(PALETTE10)
    for i in range(n_pir):
        m = labels == f"piramide_{i:02d}"
        if np.any(m):
            rgba[m, :3] = PALETTE10[i % n_pir]
    emi = DOM_EMISSION_CORTE if CORTE_ACTIVO else PIR_EMISSION
    obj = make_mesh_object("capa0_piramides", coords, [], rgba, "col_pir", coll)
    add_points_geonodes(obj, make_emission_material("mat_pir", "col_pir", emi),
                        PIR_POINT_RADIUS)
    return len(coords)


def draw_hilio(dom, coll):
    """Capa 0: esfera + flecha medial->lateral del hilio (referencia, sin corte)."""
    import mathutils
    mat_h = make_solid_material("mat_hilio", (1.0, 0.0, 0.0), strength=3.0)
    bpy.ops.mesh.primitive_uv_sphere_add(radius=HILIO_RADIO, location=HILIO)
    sph = bpy.context.active_object
    sph.name = "hilio_seno_renal"
    sph.data.materials.append(mat_h)
    link_to_collection(sph, coll)

    start = mathutils.Vector(HILIO)
    end = mathutils.Vector(HILIO_FLECHA_FIN)
    vec = end - start
    bpy.ops.mesh.primitive_cylinder_add(radius=0.5, depth=vec.length,
                                        location=(start + end) * 0.5)
    arr = bpy.context.active_object
    arr.name = "hilio_eje_medial_lateral"
    arr.rotation_mode = 'QUATERNION'
    arr.rotation_quaternion = mathutils.Vector((0, 0, 1)).rotation_difference(vec)
    arr.data.materials.append(mat_h)
    link_to_collection(arr, coll)
    return 2


def draw_glomerulos(nef, coll):
    """Capa 1: glomerulos azul (cortical) / naranja (yuxtamedular)."""
    glom = nef["glomerulos"].astype(np.float64)
    tipo = nef["tipo"]
    m = slice_mask(glom)
    glom, tipo = glom[m], tipo[m]
    rgba = colors_by_tipo(tipo, dim=1.0)
    obj = make_mesh_object("capa1_glomerulos", glom, [], rgba, "col_glom", coll)
    add_points_geonodes(obj, make_emission_material("mat_glom", "col_glom",
                        GLOM_EMISSION), GLOM_POINT_RADIUS)
    n_c = int(np.count_nonzero(tipo == "cortical"))
    n_j = int(np.count_nonzero(tipo == "yuxtamedular"))
    return len(glom), n_c, n_j


def draw_tubulos(nef, coll):
    """Capa 1: tubulos como tubos; en corte se incluyen los de glomerulo en rebanada."""
    glom = nef["glomerulos"].astype(np.float64)
    tipo = nef["tipo"]
    tub = nef["tubulos"].astype(np.float64)   # (N, NPTS, 3)
    m = slice_mask(glom)
    tipo, tub = tipo[m], tub[m]
    n_tub, npts, _ = tub.shape
    if n_tub == 0:
        return 0
    verts = tub.reshape(-1, 3)
    base = (np.arange(n_tub) * npts)[:, None]
    seg = np.arange(npts - 1)[None, :]
    a = (base + seg).reshape(-1)
    edges = np.stack([a, a + 1], axis=1).tolist()
    rgba = np.repeat(colors_by_tipo(tipo, dim=TUBE_DIM), npts, axis=0)
    obj = make_mesh_object("capa1_tubulos", verts, edges, rgba, "col_tub", coll)
    add_tube_geonodes(obj, make_emission_material("mat_tub", "col_tub",
                      TUBE_EMISSION), TUBE_RADIUS)
    return n_tub


def draw_demanda(dem, coll):
    """Capa 2: mapa de calor de la demanda metabolica de O2 sobre el parenquima.

    Pinta los puntos del dominio con un colormap frio->caliente segun la demanda
    [0,1]: medula (~0.3) fria, corteza+nefronas (~0.8-1.0) caliente. Respeta el
    CORTE (en rebanada el mapa de calor es especialmente revelador).

    Devuelve (n_puntos, demanda_min, demanda_max).
    """
    coords = dem["coords"][::SUBSAMPLE_DEMANDA]
    demanda = dem["demanda"][::SUBSAMPLE_DEMANDA]
    m = slice_mask(coords)
    coords, demanda = coords[m], demanda[m]
    rgba = heatmap_colors(demanda)
    obj = make_mesh_object("capa2_demanda", coords, [], rgba, "col_dem", coll)
    add_points_geonodes(obj, make_emission_material("mat_dem", "col_dem",
                        DEMANDA_EMISSION), DEMANDA_POINT_RADIUS)
    if len(demanda) == 0:
        return 0, 0.0, 0.0
    return len(coords), float(demanda.min()), float(demanda.max())


def draw_arterial(art, coll):
    """Capa 3a: arbol arterial como GRAFO -> tubos con calibre de Murray.

    - Malla unica (vertices=nodos, edges=aristas), color rojo arterial.
    - Geometry Nodes Curve to Mesh con radio modulado por nodo (atributo
      'radius'): a cada nodo se le asigna el radio de su arista ENTRANTE (la
      raiz toma el radio maximo del tronco). Los radios reales de Murray se
      remapean al rango visible [ART_TUBE_MIN, ART_TUBE_MAX] conservando el
      orden para LEER el calibre sin que las terminales desaparezcan.
    - Esfera marcadora en la RAIZ (hilio) + puntos en los nodos TERMINALES.
    - Respeta CORTE_ACTIVO: conserva las aristas con AMBOS extremos en la
      rebanada (la raiz, como referencia, se dibuja siempre).

    Devuelve (n_segmentos, (radio_min_mm, radio_max_mm), n_terminales).
    """
    nodos = art["nodos"].astype(np.float64)
    aristas = art["aristas"].astype(np.int64)
    radios = art["radios"].astype(np.float64)          # por ARISTA, en mm reales
    terminales = art["terminales"].astype(np.int64)
    M = len(nodos)
    root = 0                                            # capa3a enraiza en el indice 0

    # radio por NODO = radio de su arista entrante (hijo). La raiz no tiene
    # arista entrante -> le damos el radio maximo (es el tronco mas grueso).
    node_r = np.zeros(M, dtype=np.float64)
    if len(aristas):
        node_r[aristas[:, 1]] = radios
    node_r[root] = float(radios.max()) if len(radios) else ART_TUBE_MAX

    # --- CORTE: conservar aristas con ambos extremos en la rebanada ---------
    nm = slice_mask(nodos)
    if len(aristas):
        em = nm[aristas[:, 0]] & nm[aristas[:, 1]]
    else:
        em = np.zeros(0, dtype=bool)
    arist_kept = aristas[em]
    radios_kept = radios[em]
    keep_idx = np.where(nm)[0]
    remap = -np.ones(M, dtype=np.int64)
    remap[keep_idx] = np.arange(len(keep_idx))
    verts = nodos[nm]
    edges2 = remap[arist_kept] if len(arist_kept) else np.empty((0, 2), dtype=np.int64)
    node_r_kept = node_r[nm]
    n_seg = len(edges2)

    # --- remapeo de radio real -> radio visible (display) -------------------
    if node_r_kept.size:
        rmn, rmx = node_r_kept.min(), node_r_kept.max()
        if rmx > rmn:
            disp = ART_TUBE_MIN + (node_r_kept - rmn) / (rmx - rmn) * \
                (ART_TUBE_MAX - ART_TUBE_MIN)
        else:
            disp = np.full_like(node_r_kept, ART_TUBE_MAX)
    else:
        disp = node_r_kept

    # --- malla del arbol (tubos con calibre Murray) -------------------------
    if n_seg > 0:
        rgba = np.empty((len(verts), 4), dtype=np.float32)
        rgba[:, 3] = 1.0
        rgba[:, :3] = COLOR_ARTERIAL
        obj = make_mesh_object("capa3a_arterial", verts, edges2.tolist(),
                               rgba, "col_art", coll)
        add_point_float_attr(obj, "radius", disp)
        add_tree_tube_geonodes(obj, make_emission_material("mat_art", "col_art",
                               ART_EMISSION), "radius")

    # --- RAIZ (hilio): esfera marcadora (referencia, siempre visible) -------
    raiz_pos = tuple(nodos[root])
    mat_r = make_solid_material("mat_art_raiz", COLOR_ARTERIAL, strength=3.0)
    bpy.ops.mesh.primitive_uv_sphere_add(radius=ART_ROOT_RADIUS, location=raiz_pos)
    sph = bpy.context.active_object
    sph.name = "capa3a_raiz_hilio"
    sph.data.materials.append(mat_r)
    link_to_collection(sph, coll)

    # --- nodos TERMINALES: puntos pequenos (arteriolas aferentes) -----------
    term_pos = nodos[terminales]
    tmask = slice_mask(term_pos)
    term_pos = term_pos[tmask]
    n_term = len(term_pos)
    if n_term > 0:
        rgba_t = np.empty((n_term, 4), dtype=np.float32)
        rgba_t[:, 3] = 1.0
        rgba_t[:, :3] = COLOR_TERMINAL
        objt = make_mesh_object("capa3a_terminales", term_pos, [], rgba_t,
                                "col_term", coll)
        add_points_geonodes(objt, make_emission_material("mat_art_term", "col_term",
                            ART_EMISSION), ART_TERMINAL_POINT_RADIUS)

    rr = (float(radios_kept.min()), float(radios_kept.max())) if len(radios_kept) \
        else (0.0, 0.0)
    return n_seg, rr, n_term


def draw_venoso(ven, coll):
    """Capa 3b: arbol venoso (reparado, canonico) como GRAFO -> tubos con calibre.

    ESPEJA draw_arterial: malla unica (vertices=nodos, edges=aristas) en AZUL,
    Geometry Nodes Curve to Mesh con radio modulado por nodo (atributo 'radius'),
    esfera marcadora en la RAIZ (colector) + puntos en los nodos TERMINALES
    (venulas). Unica diferencia estructural: la raiz se lee del campo 'raiz' del
    .npz (el venoso enraiza en su colector), no del indice 0 fijo.

    Devuelve (n_segmentos, (radio_min_mm, radio_max_mm), n_terminales).
    """
    nodos = ven["nodos"].astype(np.float64)
    aristas = ven["aristas"].astype(np.int64)
    radios = ven["radios"].astype(np.float64)          # por ARISTA, en mm reales
    terminales = ven["terminales"].astype(np.int64)
    M = len(nodos)
    root = int(ven["raiz"]) if "raiz" in ven.files else 0

    # radio por NODO = radio de su arista entrante (hijo). La raiz no tiene
    # arista entrante -> le damos el radio maximo (es el colector mas grueso).
    node_r = np.zeros(M, dtype=np.float64)
    if len(aristas):
        node_r[aristas[:, 1]] = radios
    node_r[root] = float(radios.max()) if len(radios) else VEN_TUBE_MAX

    # --- CORTE: conservar aristas con ambos extremos en la rebanada ---------
    nm = slice_mask(nodos)
    if len(aristas):
        em = nm[aristas[:, 0]] & nm[aristas[:, 1]]
    else:
        em = np.zeros(0, dtype=bool)
    arist_kept = aristas[em]
    radios_kept = radios[em]
    keep_idx = np.where(nm)[0]
    remap = -np.ones(M, dtype=np.int64)
    remap[keep_idx] = np.arange(len(keep_idx))
    verts = nodos[nm]
    edges2 = remap[arist_kept] if len(arist_kept) else np.empty((0, 2), dtype=np.int64)
    node_r_kept = node_r[nm]
    n_seg = len(edges2)

    # --- remapeo de radio real -> radio visible (display) -------------------
    if node_r_kept.size:
        rmn, rmx = node_r_kept.min(), node_r_kept.max()
        if rmx > rmn:
            disp = VEN_TUBE_MIN + (node_r_kept - rmn) / (rmx - rmn) * \
                (VEN_TUBE_MAX - VEN_TUBE_MIN)
        else:
            disp = np.full_like(node_r_kept, VEN_TUBE_MAX)
    else:
        disp = node_r_kept

    # --- malla del arbol (tubos con calibre) --------------------------------
    if n_seg > 0:
        rgba = np.empty((len(verts), 4), dtype=np.float32)
        rgba[:, 3] = 1.0
        rgba[:, :3] = COLOR_VENOSO
        obj = make_mesh_object("capa3b_venoso", verts, edges2.tolist(),
                               rgba, "col_ven", coll)
        add_point_float_attr(obj, "radius", disp)
        add_tree_tube_geonodes(obj, make_emission_material("mat_ven", "col_ven",
                               VEN_EMISSION), "radius")

    # --- RAIZ (colector): esfera marcadora (referencia, siempre visible) -----
    raiz_pos = tuple(nodos[root])
    mat_r = make_solid_material("mat_ven_raiz", COLOR_VENOSO, strength=3.0)
    bpy.ops.mesh.primitive_uv_sphere_add(radius=VEN_ROOT_RADIUS, location=raiz_pos)
    sph = bpy.context.active_object
    sph.name = "capa3b_raiz_colector"
    sph.data.materials.append(mat_r)
    link_to_collection(sph, coll)

    # --- nodos TERMINALES: puntos pequenos (venulas) ------------------------
    term_pos = nodos[terminales]
    tmask = slice_mask(term_pos)
    term_pos = term_pos[tmask]
    n_term = len(term_pos)
    if n_term > 0:
        rgba_t = np.empty((n_term, 4), dtype=np.float32)
        rgba_t[:, 3] = 1.0
        rgba_t[:, :3] = COLOR_VENOSO_TERMINAL
        objt = make_mesh_object("capa3b_terminales", term_pos, [], rgba_t,
                                "col_ven_term", coll)
        add_points_geonodes(objt, make_emission_material("mat_ven_term", "col_ven_term",
                            VEN_EMISSION), VEN_TERMINAL_POINT_RADIUS)

    rr = (float(radios_kept.min()), float(radios_kept.max())) if len(radios_kept) \
        else (0.0, 0.0)
    return n_seg, rr, n_term


def draw_colector(col_npz, coll):
    """Capa 3c: arbol colector urinario (bosque de 10 subarboles) como GRAFO -> tubos.

    ESPEJA draw_venoso, con dos diferencias:
      1. La topologia viene de 'parent' (no 'aristas'): las aristas se derivan como
         (parent[i], i) para parent[i] >= 0; las 10 raices (parent = -1) son las papilas.
      2. Los radios ya son POR-NODO en el .npz (mm, 15-125 um) y NO se remapean min-max:
         el radio del tubo = radios * RADIO_ESCALA_COLECTOR (tamaño fisico * factor
         cosmetico solo-visual). Color AMARILLO plano, o por piramide_id si
         COLECTOR_POR_PIRAMIDE (10 tonos, para ver los territorios segregados).

    Respeta CORTE_ACTIVO (conserva aristas con ambos extremos en la rebanada).
    Devuelve (n_segmentos, (radio_min_mm, radio_max_mm) REALES, n_terminales).
    """
    nodos = col_npz["nodos"].astype(np.float64)
    parent = col_npz["parent"].astype(np.int64)
    radios = col_npz["radios"].astype(np.float64)         # POR-NODO, en mm reales
    piramide_id = col_npz["piramide_id"].astype(np.int64)
    papila_nodo = col_npz["papila_nodo"].astype(np.int64)
    terminales = col_npz["terminales"].astype(np.int64)
    M = len(nodos)

    # aristas derivadas de parent (parent[i] = -1 en las 10 raices/papilas)
    child = np.where(parent >= 0)[0]
    aristas = np.stack([parent[child], child], axis=1) if len(child) \
        else np.empty((0, 2), dtype=np.int64)

    def _col_por_nodo(idx):
        """rgb (n,3) por nodo: amarillo plano o paleta por piramide_id."""
        if COLECTOR_POR_PIRAMIDE:
            return PALETTE10[piramide_id[idx] % len(PALETTE10)]
        out = np.empty((len(idx), 3), dtype=np.float64)
        out[:] = COLOR_COLECTOR
        return out

    # --- CORTE: conservar aristas con ambos extremos en la rebanada ---------
    nm = slice_mask(nodos)
    if len(aristas):
        em = nm[aristas[:, 0]] & nm[aristas[:, 1]]
    else:
        em = np.zeros(0, dtype=bool)
    arist_kept = aristas[em]
    keep_idx = np.where(nm)[0]
    remap = -np.ones(M, dtype=np.int64)
    remap[keep_idx] = np.arange(len(keep_idx))
    verts = nodos[nm]
    edges2 = remap[arist_kept] if len(arist_kept) else np.empty((0, 2), dtype=np.int64)
    node_r_kept = radios[nm]
    n_seg = len(edges2)

    # --- radio del tubo = radio FISICO por nodo * factor cosmetico -----------
    disp = node_r_kept * RADIO_ESCALA_COLECTOR

    # --- malla del arbol (tubos con calibre por-nodo) -----------------------
    if n_seg > 0:
        rgba = np.empty((len(verts), 4), dtype=np.float32)
        rgba[:, 3] = 1.0
        rgba[:, :3] = _col_por_nodo(keep_idx)
        obj = make_mesh_object("capa3c_colector", verts, edges2.tolist(),
                               rgba, "col_col", coll)
        add_point_float_attr(obj, "radius", disp)
        add_tree_tube_geonodes(obj, make_emission_material("mat_col", "col_col",
                               COLECTOR_EMISSION), "radius")

    # --- PAPILAS (10 raices): esferas marcadoras (respetan el corte) ---------
    pap_pos = nodos[papila_nodo]
    pmask = slice_mask(pap_pos)
    for k in np.where(pmask)[0]:
        cpap = tuple(PALETTE10[k % len(PALETTE10)]) if COLECTOR_POR_PIRAMIDE else COLOR_COLECTOR
        mat_p = make_solid_material(f"mat_col_papila_{k}", cpap, strength=3.0)
        bpy.ops.mesh.primitive_uv_sphere_add(radius=COLECTOR_ROOT_RADIUS,
                                             location=tuple(pap_pos[k]))
        sph = bpy.context.active_object
        sph.name = f"capa3c_papila_{k}"
        sph.data.materials.append(mat_p)
        link_to_collection(sph, coll)

    # --- nodos TERMINALES: puntos pequenos (bocas colectoras que reciben tubulos) ---
    tmask = slice_mask(nodos[terminales])
    term_idx = terminales[tmask]
    term_pos = nodos[term_idx]
    n_term = len(term_pos)
    if n_term > 0:
        rgba_t = np.empty((n_term, 4), dtype=np.float32)
        rgba_t[:, 3] = 1.0
        rgba_t[:, :3] = _col_por_nodo(term_idx) if COLECTOR_POR_PIRAMIDE \
            else COLOR_COLECTOR_TERMINAL
        objt = make_mesh_object("capa3c_terminales", term_pos, [], rgba_t,
                                "col_col_term", coll)
        add_points_geonodes(objt, make_emission_material("mat_col_term", "col_col_term",
                            COLECTOR_EMISSION), COLECTOR_TERMINAL_POINT_RADIUS)

    rr = (float(node_r_kept.min()), float(node_r_kept.max())) if node_r_kept.size \
        else (0.0, 0.0)
    return n_seg, rr, n_term


def draw_peritubular(peri, coll):
    """Capa 3a-bis: SEGUNDA red capilar (peritubular + vasa recta) como PUNTOS.

    Dibuja los puntos de DRENAJE de la Capa 3ab (puente arterial->venoso) con el
    patron probado de Mesh to Points. Color por TIPO DE RED:
      - peritubular_cortical -> VERDE  (tapiza los tubulos en la corteza)
      - vasa_recta           -> MAGENTA/VIOLETA (se hunde en la medula siguiendo
        el asa de Henle; en CORTE sagital se ve penetrar la medula).

    Respeta CORTE_ACTIVO. Devuelve (n_total, n_peritubular, n_vasa_recta).
    """
    pts = peri["puntos_drenaje"].astype(np.float64)
    tipo_red = peri["tipo_red"]
    m = slice_mask(pts)
    pts, tipo_red = pts[m], tipo_red[m]
    if len(pts) == 0:
        return 0, 0, 0

    is_vasa = tipo_red == "vasa_recta"
    rgba = np.empty((len(pts), 4), dtype=np.float32)
    rgba[:, 3] = 1.0
    rgba[~is_vasa, :3] = COLOR_PERITUBULAR
    rgba[is_vasa, :3] = COLOR_VASA_RECTA

    obj = make_mesh_object("capa3ab_peritubular", pts, [], rgba, "col_peri", coll)
    add_points_geonodes(obj, make_emission_material("mat_peri", "col_peri",
                        PERI_EMISSION), PERI_POINT_RADIUS)

    n_vasa = int(np.count_nonzero(is_vasa))
    n_peri = int(len(pts) - n_vasa)
    return len(pts), n_peri, n_vasa


def draw_calices(cal, coll):
    """Capa 4: sistema calicial alto (papila->caliz menor->caliz mayor->pelvis->ureter)
    como CENTERLINE + RADIO. Cada arista = cilindro CONICO (tapered) entre nodos[parent]
    y nodos[child], con radio_inicio=radios[parent] y radio_fin=radios[child]. NO usa
    primitivas de camara: el gate es ver si la pelvis-como-cilindro basta a tamaño fisico
    real (radios ABSOLUTOS 0.2-4.5 mm, sin escala cosmetica).

    ESPEJA el patron de draw_colector, con dos particularidades:
      - Color por 'nivel' (0..4): gris / naranja claro / naranja / rojo / magenta (COLOR_CALIZ).
      - arista_tipo==1 (JUNCTION_DISCRETA, Bellini->copa): NO interpola el radio; dibuja un
        stub FINO r=radios[parent] tope-a-tope con un cilindro r=radios[child] -> ESCALON de
        luz visible (para inspeccionar la discontinuidad Bellini->copa).

    Respeta CORTE_ACTIVO. Devuelve (n_aristas_dibujadas, (r_min, r_max) mm, n_nodos).
    """
    import mathutils
    nodos = cal["nodos"].astype(np.float64)
    radios = cal["radios"].astype(np.float64)          # mm ABSOLUTOS (0.2 - 4.5)
    nivel = cal["nivel"].astype(np.int64)
    aristas = cal["aristas"].astype(np.int64)          # [parent, child]
    arista_tipo = cal["arista_tipo"].astype(np.int64)

    # materiales por nivel (una vez; se reutilizan en aristas y nodos)
    mats = [make_solid_material(f"mat_caliz_nivel{L}", COLOR_CALIZ[L], strength=1.5)
            for L in range(len(COLOR_CALIZ))]

    def _tubo(p_ini, p_fin, r_ini, r_fin, mat, nombre):
        start = mathutils.Vector(tuple(p_ini)); end = mathutils.Vector(tuple(p_fin))
        vec = end - start
        if vec.length < 1e-9:
            return
        # cono tapered: radius1 (base, -Z) en p_ini=parent ; radius2 (+Z) en p_fin=child
        bpy.ops.mesh.primitive_cone_add(
            vertices=CALIZ_CONE_VERTS, radius1=float(r_ini), radius2=float(r_fin),
            depth=vec.length, location=tuple((start + end) * 0.5))
        obj = bpy.context.active_object
        obj.name = nombre
        obj.rotation_mode = 'QUATERNION'
        obj.rotation_quaternion = mathutils.Vector((0, 0, 1)).rotation_difference(vec)
        obj.data.materials.append(mat)
        link_to_collection(obj, coll)

    nm = slice_mask(nodos)
    n_seg = 0
    for e in range(len(aristas)):
        pa, ch = int(aristas[e, 0]), int(aristas[e, 1])
        if not (nm[pa] and nm[ch]):
            continue
        p0, p1 = nodos[pa], nodos[ch]
        if int(arista_tipo[e]) == 1:
            # ESCALON DISCRETO (Bellini->copa): stub fino r cte + copa r cte, tope-a-tope
            mid = p0 + CALIZ_BELLINI_FRAC * (p1 - p0)
            _tubo(p0, mid, radios[pa], radios[pa], mats[int(nivel[pa])],
                  f"capa4_bellini_{e}")     # stub fino r=radios[parent] (nivel papila=gris)
            _tubo(mid, p1, radios[ch], radios[ch], mats[int(nivel[ch])],
                  f"capa4_copa_{e}")        # copa r=radios[child] -> el salto es el ESCALON
        else:
            # TAPER continuo: cono radios[parent] -> radios[child]; color por nivel del hijo
            _tubo(p0, p1, radios[pa], radios[ch], mats[int(nivel[ch])], f"capa4_arista_{e}")
        n_seg += 1

    # marcadores de NODO (esferitas por nivel) para leer la jerarquia calicial
    n_nod = 0
    for i in np.where(nm)[0]:
        bpy.ops.mesh.primitive_uv_sphere_add(radius=CALIZ_NODE_MARKER, location=tuple(nodos[i]))
        s = bpy.context.active_object
        s.name = f"capa4_nodo_{i}_nivel{int(nivel[i])}"
        s.data.materials.append(mats[int(nivel[i])])
        link_to_collection(s, coll)
        n_nod += 1

    rr = (float(radios[nm].min()), float(radios[nm].max())) if nm.any() else (0.0, 0.0)
    return n_seg, rr, n_nod


# ============================================================================
#  VIEWPORT
# ============================================================================
def setup_viewport(view_axis=None):
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


# ============================================================================
#  MAIN
# ============================================================================
def main():
    print("=" * 70)
    print("  VISUALIZADOR MAESTRO - Gemelo Digital Bio-Kidney AI")
    print("=" * 70)

    # --- carga perezosa -----------------------------------------------------
    dom, dom_path = try_load(DOM_NPZ)
    nef, nef_path = try_load(NEF_NPZ)
    dem, dem_path = (None, None)
    if VER_DEMANDA:
        dem, dem_path = try_load(DEM_NPZ)
    art, art_path = (None, None)
    if VER_ARTERIAL:
        art, art_path = try_load(ART_NPZ)
    ven, ven_path = (None, None)
    if VER_VENOSO:
        ven, ven_path = try_load(VEN_NPZ)
    colec, colec_path = (None, None)
    if VER_COLECTOR:
        colec, colec_path = try_load(COLECTOR_NPZ)
    peri, peri_path = (None, None)
    if VER_PERITUBULAR:
        peri, peri_path = try_load(PERI_NPZ)
    cal, cal_path = (None, None)
    if VER_CALICES:
        cal, cal_path = try_load(CALICES_NPZ)

    # flags efectivos segun disponibilidad
    ver_dominio = VER_DOMINIO
    ver_piramides = VER_PIRAMIDES
    ver_hilio = VER_HILIO
    ver_glom = VER_GLOMERULOS
    ver_tub = VER_TUBULOS
    ver_dem = VER_DEMANDA
    ver_art = VER_ARTERIAL
    ver_ven = VER_VENOSO
    ver_colec = VER_COLECTOR
    ver_peri = VER_PERITUBULAR
    ver_cal = VER_CALICES

    cargadas, desactivadas = [], []
    if dom is not None:
        cargadas.append(f"{DOM_NPZ}")
    else:
        for fl, nm in ((VER_DOMINIO, "VER_DOMINIO"), (VER_PIRAMIDES, "VER_PIRAMIDES"),
                       (VER_HILIO, "VER_HILIO")):
            if fl:
                desactivadas.append(f"{nm} (falta {DOM_NPZ})")
        ver_dominio = ver_piramides = ver_hilio = False

    if nef is not None:
        cargadas.append(f"{NEF_NPZ}")
    else:
        for fl, nm in ((VER_GLOMERULOS, "VER_GLOMERULOS"), (VER_TUBULOS, "VER_TUBULOS")):
            if fl:
                desactivadas.append(f"{nm} (falta {NEF_NPZ})")
        ver_glom = ver_tub = False

    if VER_DEMANDA:
        if dem is not None:
            cargadas.append(f"{DEM_NPZ}")
        else:
            desactivadas.append(f"VER_DEMANDA (falta {DEM_NPZ})")
            ver_dem = False

    if VER_ARTERIAL:
        if art is not None:
            cargadas.append(f"{ART_NPZ}")
        else:
            desactivadas.append(f"VER_ARTERIAL (falta {ART_NPZ})")
            ver_art = False

    if VER_VENOSO:
        if ven is not None:
            cargadas.append(f"{VEN_NPZ}")
        else:
            desactivadas.append(f"VER_VENOSO (falta {VEN_NPZ})")
            ver_ven = False

    if VER_COLECTOR:
        if colec is not None:
            cargadas.append(f"{COLECTOR_NPZ}")
        else:
            desactivadas.append(f"VER_COLECTOR (falta {COLECTOR_NPZ})")
            ver_colec = False

    if VER_PERITUBULAR:
        if peri is not None:
            cargadas.append(f"{PERI_NPZ}")
        else:
            desactivadas.append(f"VER_PERITUBULAR (falta {PERI_NPZ})")
            ver_peri = False

    if VER_CALICES:
        if cal is not None:
            cargadas.append(f"{CALICES_NPZ}")
        else:
            desactivadas.append(f"VER_CALICES (falta {CALICES_NPZ})")
            ver_cal = False

    print("  Archivos cargados :", ", ".join(cargadas) if cargadas else "(ninguno)")
    if desactivadas:
        print("  Capas desactivadas:", "; ".join(desactivadas))
    if CORTE_ACTIVO:
        print(f"  MODO CORTE: eje {CORTE_EJE}, centro {CORTE_CENTRO} mm, "
              f"grosor {CORTE_GROSOR} mm")
    print("-" * 70)

    # --- escena -------------------------------------------------------------
    clear_scene()
    col_dom = get_collection(COL_DOMINIO)
    col_nef = get_collection(COL_NEFRONAS)

    counts = {}
    if ver_dominio:
        counts["dominio (pts)"] = draw_dominio(dom, col_dom)
    if ver_piramides:
        counts["piramides (pts)"] = draw_piramides(dom, col_dom)
    if ver_hilio:
        counts["hilio (objs)"] = draw_hilio(dom, col_dom)
    if ver_glom:
        ng, nc, nj = draw_glomerulos(nef, col_nef)
        counts["glomerulos"] = ng
        counts["  corticales"] = nc
        counts["  yuxtamedulares"] = nj
    if ver_tub:
        counts["tubulos"] = draw_tubulos(nef, col_nef)
    dem_range = None
    if ver_dem:
        col_dem = get_collection(COL_DEMANDA)
        nd_dem, dmin, dmax = draw_demanda(dem, col_dem)
        counts["demanda (pts)"] = nd_dem
        dem_range = (dmin, dmax)
    art_radio_range = None
    if ver_art:
        col_art = get_collection(COL_ARTERIAL)
        n_seg, art_radio_range, n_term = draw_arterial(art, col_art)
        counts["arterial (segmentos)"] = n_seg
        counts["  terminales (arteriolas)"] = n_term
    ven_radio_range = None
    if ver_ven:
        col_ven = get_collection(COL_VENOSO)
        n_seg_v, ven_radio_range, n_term_v = draw_venoso(ven, col_ven)
        counts["venoso (segmentos)"] = n_seg_v
        counts["  terminales (venulas)"] = n_term_v
    colec_radio_range = None
    if ver_colec:
        col_colec = get_collection(COL_COLECTOR)
        n_seg_c, colec_radio_range, n_term_c = draw_colector(colec, col_colec)
        counts["colector (segmentos)"] = n_seg_c
        counts["  bocas colectoras (term.)"] = n_term_c
    if ver_peri:
        col_peri = get_collection(COL_PERITUBULAR)
        n_peri_tot, n_peri_cort, n_peri_vasa = draw_peritubular(peri, col_peri)
        counts["peritubular (pts drenaje)"] = n_peri_tot
        counts["  peritubular_cortical (verde)"] = n_peri_cort
        counts["  vasa_recta (violeta)"] = n_peri_vasa
    cal_radio_range = None
    if ver_cal:
        col_cal = get_collection(COL_CALICES)
        n_cal_seg, cal_radio_range, n_cal_nod = draw_calices(cal, col_cal)
        counts["calices (aristas)"] = n_cal_seg
        counts["  nodos calices (por nivel)"] = n_cal_nod

    # --- vista --------------------------------------------------------------
    view = _EJE_VIEW.get(CORTE_EJE) if CORTE_ACTIVO else None
    setup_viewport(view_axis=view)

    # --- reporte ------------------------------------------------------------
    print("  Capas dibujadas:")
    for k, v in counts.items():
        print(f"     {k:18s}: {v}")
    if dem_range is not None:
        print(f"  Mapa de calor demanda: rango representado "
              f"[{dem_range[0]:.3f}, {dem_range[1]:.3f}] (frio->caliente)")
    if art_radio_range is not None:
        print(f"  Arbol arterial: rango de radios REALES (Murray) "
              f"[{art_radio_range[0]:.4f}, {art_radio_range[1]:.4f}] mm "
              f"-> remapeado a tubo visible [{ART_TUBE_MIN}, {ART_TUBE_MAX}] mm")
    if ven_radio_range is not None:
        print(f"  Arbol venoso: rango de radios REALES "
              f"[{ven_radio_range[0]:.4f}, {ven_radio_range[1]:.4f}] mm "
              f"-> remapeado a tubo visible [{VEN_TUBE_MIN}, {VEN_TUBE_MAX}] mm")
    if colec_radio_range is not None:
        print(f"  Arbol colector: radios REALES (mm, POR-NODO) "
              f"[{colec_radio_range[0]:.4f}, {colec_radio_range[1]:.4f}] mm "
              f"-> tubo = radios x {RADIO_ESCALA_COLECTOR} (cosmetico) = "
              f"[{colec_radio_range[0]*RADIO_ESCALA_COLECTOR:.4f}, "
              f"{colec_radio_range[1]*RADIO_ESCALA_COLECTOR:.4f}] mm"
              f"{'   [color por piramide]' if COLECTOR_POR_PIRAMIDE else '   [amarillo plano]'}")
    if cal_radio_range is not None:
        print(f"  Sistema calicial (Capa 4): radios ABSOLUTOS (mm) "
              f"[{cal_radio_range[0]:.3f}, {cal_radio_range[1]:.3f}] mm "
              f"(sin escala cosmetica)   [color por nivel 0..4]")
    print(f"  Vista: {'perpendicular ' + str(view) if view else 'libre'}")
    cols = [COL_DOMINIO, COL_NEFRONAS]
    if ver_dem:
        cols.append(COL_DEMANDA)
    if ver_art:
        cols.append(COL_ARTERIAL)
    if ver_ven:
        cols.append(COL_VENOSO)
    if ver_colec:
        cols.append(COL_COLECTOR)
    if ver_peri:
        cols.append(COL_PERITUBULAR)
    if ver_cal:
        cols.append(COL_CALICES)
    print("  Colecciones:", " | ".join(cols))
    print("=" * 70)


if __name__ == "__main__":
    main()
