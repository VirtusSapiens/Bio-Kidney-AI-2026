import bpy
import math
import os

def limpiar():
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete(use_global=False)

def material(nombre, color, alpha=1.0):
    mat = bpy.data.materials.new(name=nombre)
    mat.use_nodes = True
    if alpha < 1.0:
        mat.blend_method = "BLEND"
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links
    nodes.clear()
    out = nodes.new("ShaderNodeOutputMaterial")
    princ = nodes.new("ShaderNodeBsdfPrincipled")
    princ.inputs["Base Color"].default_value = color
    princ.inputs["Alpha"].default_value = alpha
    princ.inputs["Roughness"].default_value = 0.25
    links.new(princ.outputs["BSDF"], out.inputs["Surface"])
    return mat

def main():
    print("Creando capsula renal...")
    limpiar()

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
    mod = rinon.modifiers.new(name="Hilio", type="BOOLEAN")
    mod.operation = "DIFFERENCE"
    mod.object = hilio
    bpy.ops.object.modifier_apply(modifier="Hilio")
    hilio.hide_set(True)
    hilio.hide_render = True

    mat = material("Capsula_Renal", (0.85, 0.55, 0.40, 1.0), alpha=0.20)
    rinon.data.materials.append(mat)

    bpy.ops.object.light_add(type="AREA", location=(80, -80, 100))
    l1 = bpy.context.active_object
    l1.data.energy = 50000
    l1.data.size = 80

    bpy.ops.object.light_add(type="AREA", location=(-60, 50, 50))
    l2 = bpy.context.active_object
    l2.data.energy = 15000
    l2.data.size = 50

    bpy.data.worlds["World"].node_tree.nodes["Background"].inputs[0].default_value = (0.02, 0.04, 0.10, 1)
    bpy.data.worlds["World"].node_tree.nodes["Background"].inputs[1].default_value = 0.2

    bpy.ops.object.camera_add(location=(90, -90, 60))
    cam = bpy.context.active_object
    cam.name = "Camara_BioKidney"
    cam.rotation_euler = (math.radians(55), 0, math.radians(45))
    cam.data.lens = 85
    bpy.context.scene.camera = cam

    salida = os.path.expanduser("~/Escritorio/BioKidney-AI/03_modelos_3d/rinon_capsula.blend")
    bpy.ops.wm.save_as_mainfile(filepath=salida)
    print("Listo -> " + salida)

main()
