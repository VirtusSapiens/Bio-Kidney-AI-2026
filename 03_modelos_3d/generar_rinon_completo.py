import bpy
import csv
import math
import os
import numpy as np

CSV_ARBOL    = os.path.expanduser(
    "~/Escritorio/BioKidney-AI/02_vascular_cco/arbol_vascular_cco_v7.csv")
SALIDA_BLEND = os.path.expanduser(
    "~/Escritorio/BioKidney-AI/03_modelos_3d/rinon_biokidney_ai_v7.blend")

PUNTOS_CURVA = 6
SEMILLA      = 42
np.random.seed(SEMILLA)

def spline_bezier(p0, p1, curvatura=0.18, n=6):
    eje = p1 - p0
    L   = np.linalg.norm(eje)
    if L < 1e-10:
        return np.array([p0, p1])
    en    = eje / L
    perp  = np.array([-en[1], en[0], 0.0])
    if np.linalg.norm(perp) < 1e-10:
        perp = np.array([0., -en[2], en[1]])
    perp  = perp / np.linalg.norm(perp)
    perp2 = np.cross(en, perp)
    pc = ((p0 + p1) / 2
          + perp  * L * curvatura * np.random.uniform(-1, 1)
          + perp2 * L * curvatura * 0.5 * np.random.uniform(-1, 1))
    t = np.linspace(0, 1, n)
    return (np.outer((1-t)**2, p0) +
            np.outer(2*(1-t)*t, pc) +
            np.outer(t**2, p1))

def limpiar_escena():
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False)
    for col in bpy.data.collections:
        bpy.data.collections.remove(col)

def crear_material(nombre, color, alpha=1.0):
    mat = bpy.data.materials.new(name=nombre)
    mat.use_nodes = True
    if alpha < 1.0:
        mat.blend_method = 'BLEND'
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links
    nodes.clear()
    out   = nodes.new('ShaderNodeOutputMaterial')
    princ = nodes.new('ShaderNodeBsdfPrincipled')
    princ.inputs['Base Color'].default_value = color
    princ.inputs['Alpha'].default_value      = alpha
    princ.inputs['Roughness'].default_value  = 0.25
    links.new(princ.outputs['BSDF'], out.inputs['Surface'])
    out.location   = (300, 0)
    princ.location = (0, 0)
    return mat

def crear_capsula():
    bpy.ops.mesh.primitive_uv_sphere_add(
        radius=27.5, segments=64, ring_count=32,
        location=(0, 0, 0))
    rinon = bpy.context.active_object
    rinon.name = "Capsula_Renal"
    rinon.scale = (1.0, 0.545, 0.454)
    bpy.ops.object.transform_apply(scale=True)
    bpy.ops.mesh.primitive_uv_sphere_add(
        radius=10.0, segments=32, ring_count=16,
        location=(0.0, -25.0, 0.0))
    hilio = bpy.context.active_object
    hilio.name = "Hilio_Cutter"
    rinon.select_set(True)
    bpy.context.view_layer.objects.active = rinon
    mod = rinon.modifiers.new(name="Hilio", type='BOOLEAN')
    mod.operation = 'DIFFERENCE'
    mod.object    = hilio
    bpy.ops.object.modifier_apply(modifier="Hilio")
    hilio.hide_set(True)
    hilio.hide_render = True
    mat = crear_material("Capsula_Renal",
                         (0.85, 0.55, 0.40, 1.0), alpha=0.20)
    rinon.data.materials.append(mat)
    return rinon

def cargar_segmentos(ruta):
    segs = []
    with open(ruta, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            segs.append({
                'inicio'  : np.array([float(row['x1_mm']),
                                      float(row['y1_mm']),
                                      float(row['z1_mm'])]),
                'fin'     : np.array([float(row['x2_mm']),
                                      float(row['y2_mm']),
                                      float(row['z2_mm'])]),
                'radio'   : float(row['radio_um']) / 1000,
                'sistema' : row['sistema'],
            })
    return segs

def crear_tubo_curvo(seg, mat):
    curva = spline_bezier(
        seg['inicio'], seg['fin'],
        curvatura=0.20, n=PUNTOS_CURVA)
    radio_vis = max(seg['radio'], 0.15)
    objetos   = []
    for i in range(len(curva) - 1):
        p0 = curva[i]
        p1 = curva[i+1]
        eje  = p1 - p0
        long = np.linalg.norm(eje)
        if long < 0.01:
            continue
        centro = (p0 + p1) / 2
        bpy.ops.mesh.primitive_cylinder_add(
            radius=radio_vis,
            depth=long,
            vertices=5,
            location=(float(centro[0]),
                      float(centro[1]),
                      float(centro[2])))
        obj = bpy.context.active_object
        eje_n = eje / long
        z     = np.array([0.0, 0.0, 1.0])
        if abs(np.dot(eje_n, z)) < 0.9999:
            rot_eje    = np.cross(z, eje_n)
            rot_eje    = rot_eje / np.linalg.norm(rot_eje)
            rot_angulo = math.acos(np.clip(np.dot(z, eje_n), -1, 1))
            obj.rotation_mode = 'AXIS_ANGLE'
            obj.rotation_axis_angle = (
                rot_angulo,
                float(rot_eje[0]),
                float(rot_eje[1]),
                float(rot_eje[2]))
        elif np.dot(eje_n, z) < 0:
            obj.rotation_euler[0] = math.pi
        obj.data.materials.append(mat)
        objetos.append(obj)
    return objetos

def crear_arbol_vascular(segs):
    print("  Creando " + str(len(segs)) + " segmentos con curvas Bezier...")
    mat_art = crear_material("Arteria",  (0.85, 0.08, 0.08, 1.0))
    mat_ven = crear_material("Vena",     (0.08, 0.08, 0.90, 1.0))
    mat_col = crear_material("Colector", (0.90, 0.78, 0.10, 1.0))
    todos_art = []
    todos_ven = []
    todos_col = []
    for i, seg in enumerate(segs):
        if seg['sistema'] == 'art':
            mat = mat_art
        elif seg['sistema'] == 'ven':
            mat = mat_ven
        else:
            mat = mat_col
        objs = crear_tubo_curvo(seg, mat)
        if seg['sistema'] == 'art':
            todos_art.extend(objs)
        elif seg['sistema'] == 'ven':
            todos_ven.extend(objs)
        else:
            todos_col.extend(objs)
        if (i + 1) % 100 == 0:
            print("  -> " + str(i+1) + "/" + str(len(segs)) + " segmentos")
    for nombre, lista in [("ArbolArterial", todos_art),
                           ("ArbolVenoso",   todos_ven),
                           ("SistemaColector", todos_col)]:
        if lista:
            bpy.ops.object.select_all(action='DESELECT')
            for o in lista:
                o.select_set(True)
            bpy.context.view_layer.objects.active = lista[0]
            bpy.ops.object.join()
            bpy.context.active_object.name = nombre
            print("  " + nombre + ": " + str(len(lista)) + " objetos unidos")

def configurar_iluminacion():
    bpy.ops.object.light_add(type='AREA', location=(80, -80, 100))
    l1 = bpy.context.active_object
    l1.data.energy = 50000
    l1.data.size   = 80
    l1.name        = "Luz_Principal"
    bpy.ops.object.light_add(type='AREA', location=(-60, 50, 50))
    l2 = bpy.context.active_object
    l2.data.energy = 15000
    l2.data.size   = 50
    l2.name        = "Luz_Relleno"
    bpy.data.worlds["World"].node_tree.nodes[
        "Background"].inputs[0].default_value = (0.02, 0.04, 0.10, 1)
    bpy.data.worlds["World"].node_tree.nodes[
        "Background"].inputs[1].default_value = 0.2

def configurar_camara():
    bpy.ops.object.camera_add(location=(90, -90, 60))
    cam = bpy.context.active_object
    cam.name = "Camara_BioKidney"
    cam.rotation_euler = (
        math.radians(55),
        math.radians(0),
        math.radians(45))
    cam.data.lens = 85
    bpy.context.scene.camera = cam

def main():
    print("BIO-KIDNEY AI 2026 - GENERADOR 3D v7 CORREGIDO")
    print("Curvas Bezier - Hilio correcto - Escala mm")
    print("[1/5] Limpiando escena...")
    limpiar_escena()
    print("[2/5] Creando capsula renal...")
    crear_capsula()
    print("[3/5] Cargando segmentos CCO v7...")
    segs = cargar_segmentos(os.path.expanduser(CSV_ARBOL))
    print(str(len(segs)) + " segmentos cargados")
    print("[4/5] Construyendo arbol vascular con curvas...")
    crear_arbol_vascular(segs)
    print("[5/5] Iluminacion y camara...")
    configurar_iluminacion()
    configurar_camara()
    bpy.ops.wm.save_as_mainfile(
        filepath=os.path.expanduser(SALIDA_BLEND))
    print("Guardado -> " + SALIDA_BLEND)

main()
