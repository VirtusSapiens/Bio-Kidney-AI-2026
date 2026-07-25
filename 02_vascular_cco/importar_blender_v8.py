"""
╔══════════════════════════════════════════════════════════════════════════════╗
║  Bio-Kidney AI 2026 — Importador Blender (BPy) optimizado                  ║
║  Entrada: renal_data_v1.json (CCO v8)                                      ║
║                                                                            ║
║  Estrategia de rendimiento:                                                ║
║    - UNA sola malla (1 objeto Blender, no 1902)                            ║
║    - Tubos de 4 lados (mínima geometría por segmento)                      ║
║    - Pre-allocación NumPy → from_pydata() en batch                         ║
║    - Sin Bezier, sin subdivisión — segmentos rectos                        ║
║    - 3 materiales (arterial/venoso/colector) vía face-material-index       ║
║    - Pico RAM estimado: ~50 MB (seguro en 16 GB)                           ║
║                                                                            ║
║  Uso en Blender:                                                           ║
║    1. Abrir Blender → Edit → Preferences → Scripting:                      ║
║       desmarcar "Auto Run Python Scripts" si da error de seguridad         ║
║    2. Scripting workspace → Open → importar_blender_v8.py → Run Script     ║
║    O desde terminal:                                                       ║
║       blender --python importar_blender_v8.py                              ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import bpy
import bmesh
import json
import os
import time
import math
from mathutils import Vector, Matrix

# ── CONFIGURACIÓN ────────────────────────────────────────────────────────────
JSON_PATH = os.path.expanduser(
    "~/Escritorio/BioKidney-AI/02_vascular_cco/renal_data_v1.json"
)

LADOS            = 4       # 4 caras por tubo (cuadrado → ligero)
ESCALA_RADIO     = 1.0     # multiplicador de radio visual
RADIO_MIN_VIS    = 0.04    # mm — radio mínimo visible en viewport
UNIDAD           = "mm"    # coordenadas del JSON ya están en mm

# Colores materiales (RGBA)
COLORES = {
    "art": (0.92, 0.12, 0.08, 1.0),   # rojo arterial
    "ven": (0.10, 0.18, 0.90, 1.0),   # azul venoso
    "col": (0.95, 0.82, 0.12, 1.0),   # amarillo colector
}


def log(msg):
    """Print con timestamp para la System Console de Blender."""
    t = time.strftime("%H:%M:%S")
    print(f"  [{t}] {msg}", flush=True)


# ── CREAR MATERIALES ─────────────────────────────────────────────────────────
def crear_materiales():
    mats = {}
    nombres = {
        "art": "BK_Arterial",
        "ven": "BK_Venoso",
        "col": "BK_Colector",
    }
    for sistema, nombre in nombres.items():
        mat = bpy.data.materials.get(nombre)
        if mat is None:
            mat = bpy.data.materials.new(name=nombre)
        mat.use_nodes = True
        bsdf = mat.node_tree.nodes.get("Principled BSDF")
        if bsdf:
            color = COLORES[sistema]
            bsdf.inputs["Base Color"].default_value = color
            bsdf.inputs["Roughness"].default_value = 0.55
            # Subsurface para look orgánico
            bsdf.inputs["Subsurface Weight"].default_value = 0.15
            bsdf.inputs["Subsurface Radius"].default_value = (
                color[0] * 0.5, color[1] * 0.5, color[2] * 0.5
            )
        mats[sistema] = mat
    return mats


# ── GENERAR ANILLO DE VÉRTICES ───────────────────────────────────────────────
def anillo(centro, eje, radio, n_lados):
    """
    Genera n_lados vértices en un anillo perpendicular al eje.
    Retorna lista de Vector.
    """
    eje_n = eje.normalized()

    # Buscar vector perpendicular estable
    if abs(eje_n.z) < 0.9:
        up = Vector((0, 0, 1))
    else:
        up = Vector((1, 0, 0))

    perp1 = eje_n.cross(up).normalized()
    perp2 = eje_n.cross(perp1).normalized()

    verts = []
    for i in range(n_lados):
        ang = 2.0 * math.pi * i / n_lados
        offset = (perp1 * math.cos(ang) + perp2 * math.sin(ang)) * radio
        verts.append(centro + offset)
    return verts


# ── CONSTRUIR MALLA BATCH ───────────────────────────────────────────────────
def construir_malla(segmentos, mats):
    """
    Construye UNA sola malla con todos los tubos vasculares.
    Cada segmento = 2 anillos de LADOS vértices + LADOS caras quad.
    """
    n_segs = len(segmentos)
    n_verts_total = n_segs * LADOS * 2
    n_faces_total = n_segs * LADOS

    log(f"Construyendo malla: {n_segs} segmentos -> "
        f"{n_verts_total} verts, {n_faces_total} quads")

    # Pre-alocar listas
    all_verts = []
    all_faces = []
    face_mat_indices = []

    # Mapa sistema → índice de material
    mat_keys = list(mats.keys())
    mat_idx = {s: i for i, s in enumerate(mat_keys)}

    t0 = time.time()
    v_offset = 0

    for seg_i, seg in enumerate(segmentos):
        # Coordenadas (ya en mm)
        p0 = Vector((seg["x1_mm"], seg["y1_mm"], seg["z1_mm"]))
        p1 = Vector((seg["x2_mm"], seg["y2_mm"], seg["z2_mm"]))

        eje = p1 - p0
        if eje.length < 1e-6:
            continue

        # Radio con escala y mínimo visible
        radio = max(seg["radio_um"] / 1000.0 * ESCALA_RADIO, RADIO_MIN_VIS)

        # Generar anillos
        ring0 = anillo(p0, eje, radio, LADOS)
        ring1 = anillo(p1, eje, radio, LADOS)

        # Añadir vértices
        all_verts.extend(ring0)
        all_verts.extend(ring1)

        # Caras quad: conectar anillo 0 con anillo 1
        sistema = seg.get("sistema", "art")
        m_idx = mat_idx.get(sistema, 0)

        for k in range(LADOS):
            k_next = (k + 1) % LADOS
            a = v_offset + k
            b = v_offset + k_next
            c = v_offset + LADOS + k_next
            d = v_offset + LADOS + k
            all_faces.append((a, b, c, d))
            face_mat_indices.append(m_idx)

        v_offset += LADOS * 2

        # Progreso cada 200 segmentos
        if (seg_i + 1) % 200 == 0 or seg_i == n_segs - 1:
            elapsed = time.time() - t0
            pct = 100 * (seg_i + 1) / n_segs
            rate = (seg_i + 1) / max(elapsed, 0.01)
            remaining = (n_segs - seg_i - 1) / max(rate, 0.01)
            log(f"  Tubos: {seg_i+1:5d}/{n_segs} ({pct:5.1f}%) | "
                f"{elapsed:.1f}s elapsed | ~{remaining:.0f}s restante")

    log(f"Geometria lista: {len(all_verts)} vertices, {len(all_faces)} quads")

    # ── Crear mesh Blender con from_pydata ──
    log("Creando objeto Blender...")
    mesh = bpy.data.meshes.new("BK_VascularTree_v8")

    # Convertir Vectors a tuplas para from_pydata
    verts_tuples = [(v.x, v.y, v.z) for v in all_verts]
    mesh.from_pydata(verts_tuples, [], all_faces)

    # Asignar materiales
    for mat in mats.values():
        mesh.materials.append(mat)

    # Asignar material index por cara
    mesh.update()
    for i, poly in enumerate(mesh.polygons):
        if i < len(face_mat_indices):
            poly.material_index = face_mat_indices[i]

    # Smooth shading
    for poly in mesh.polygons:
        poly.use_smooth = True

    mesh.update()
    mesh.validate()

    # Crear objeto y linkar a escena
    obj = bpy.data.objects.new("BioKidney_Vascular_v8", mesh)
    bpy.context.collection.objects.link(obj)
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)

    log(f"Objeto '{obj.name}' creado en escena")
    return obj


# ── MAIN ─────────────────────────────────────────────────────────────────────
def main():
    print("\n" + "=" * 60)
    print("  Bio-Kidney AI 2026 — Importador Blender v8")
    print("  1 objeto | 4-sided tubes | batch geometry")
    print("=" * 60)

    # ── Cargar JSON ──
    log(f"Cargando {JSON_PATH}...")
    if not os.path.exists(JSON_PATH):
        log(f"ERROR: no se encontro {JSON_PATH}")
        return

    with open(JSON_PATH, "r") as f:
        data = json.load(f)

    segmentos = data["segments"]
    stats = data.get("statistics", {})
    log(f"JSON cargado: {len(segmentos)} segmentos | "
        f"v{data.get('version', '?')}")
    log(f"  Vasculares: {stats.get('n_vascular', '?')} | "
        f"Terminales art: {stats.get('n_terminals_art', '?')} | "
        f"Cobertura: {stats.get('coverage_pct', '?')}%")

    # ── Limpiar escena (solo objetos mesh previos de BioKidney) ──
    log("Limpiando objetos BioKidney previos...")
    for obj in list(bpy.data.objects):
        if obj.name.startswith("BioKidney_") or obj.name.startswith("BK_"):
            bpy.data.objects.remove(obj, do_unlink=True)

    # Purgar meshes huerfanas
    for mesh in list(bpy.data.meshes):
        if mesh.users == 0:
            bpy.data.meshes.remove(mesh)

    # ── Crear materiales ──
    log("Creando materiales (arterial/venoso/colector)...")
    mats = crear_materiales()

    # ── Construir malla ──
    t_total = time.time()
    obj = construir_malla(segmentos, mats)

    # ── Estadísticas finales ──
    n_verts = len(obj.data.vertices)
    n_polys = len(obj.data.polygons)
    elapsed = time.time() - t_total

    # Contar por sistema
    n_art = sum(1 for s in segmentos if s["sistema"] == "art")
    n_ven = sum(1 for s in segmentos if s["sistema"] == "ven")
    n_col = sum(1 for s in segmentos if s["sistema"] == "col")

    print("\n" + "=" * 60)
    print("  IMPORTACION COMPLETADA")
    print("=" * 60)
    log(f"  Segmentos arteriales  : {n_art}")
    log(f"  Segmentos venosos     : {n_ven}")
    log(f"  Segmentos colectores  : {n_col}")
    log(f"  Vertices totales      : {n_verts:,}")
    log(f"  Quads totales         : {n_polys:,}")
    log(f"  Materiales            : {len(obj.data.materials)}")
    log(f"  Tiempo total          : {elapsed:.1f}s")
    log(f"  RAM estimada malla    : ~{n_verts * 12 / 1024 / 1024:.1f} MB")
    print("=" * 60)
    log("Tip: activa Material Preview (Z) para ver colores")
    log("Tip: si necesitas mas grosor, selecciona el objeto y")
    log("     añade modifier Solidify con grosor 0.1-0.5 mm\n")


# ── EJECUCIÓN ────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    main()
else:
    # Cuando se ejecuta desde Blender Text Editor (Run Script)
    main()
