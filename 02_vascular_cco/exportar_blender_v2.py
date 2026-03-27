import csv
import math
import os
import numpy as np

# ── CONFIGURACIÓN ──────────────────────────────────────────────────────────────
ENTRADA = os.path.expanduser(
    "~/Escritorio/BioKidney-AI/02_vascular_cco/arbol_vascular_cco_v7.csv")
SALIDA  = os.path.expanduser(
    "~/Escritorio/BioKidney-AI/02_vascular_cco/arbol_vascular_cco_v7_curvo.obj")

LADOS         = 6    # caras por cilindro
PUNTOS_CURVA  = 8    # puntos de curvatura por segmento
SEMILLA       = 42

np.random.seed(SEMILLA)

# ── SPLINE BEZIER ──────────────────────────────────────────────────────────────
def spline_bezier(p0, p1, curvatura=0.18, n=8):
    """Genera curva Bezier entre dos puntos"""
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

    # Punto de control con curvatura anatómica
    pc = ((p0 + p1) / 2
          + perp  * L * curvatura * np.random.uniform(-1, 1)
          + perp2 * L * curvatura * 0.5 * np.random.uniform(-1, 1))

    t = np.linspace(0, 1, n)
    return (np.outer((1-t)**2, p0) +
            np.outer(2*(1-t)*t, pc) +
            np.outer(t**2, p1))

# ── CARGAR SEGMENTOS ───────────────────────────────────────────────────────────
def cargar_segmentos(ruta):
    segmentos = []
    with open(ruta, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            segmentos.append({
                'inicio'  : np.array([float(row['x1_mm']) / 1000,
                                      float(row['y1_mm']) / 1000,
                                      float(row['z1_mm']) / 1000]),
                'fin'     : np.array([float(row['x2_mm']) / 1000,
                                      float(row['y2_mm']) / 1000,
                                      float(row['z2_mm']) / 1000]),
                'radio'   : float(row['radio_um']) / 1_000_000,
                'sistema' : row['sistema'],
                'nivel'   : int(row['nivel']),
            })
    return segmentos

# ── CILINDRO CURVO ─────────────────────────────────────────────────────────────
def cilindro_curvo(curva, radio, lados):
    """
    Genera un tubo siguiendo la curva Bezier.
    Cada par de puntos consecutivos genera un segmento de cilindro.
    """
    vertices = []
    caras    = []
    offset   = 0

    for i in range(len(curva) - 1):
        p0 = curva[i]
        p1 = curva[i + 1]
        eje = p1 - p0
        L   = np.linalg.norm(eje)
        if L < 1e-10:
            continue

        eje_n = eje / L
        perp  = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(perp, eje_n)) > 0.9:
            perp = np.array([0.0, 1.0, 0.0])
        perp  = perp - np.dot(perp, eje_n) * eje_n
        perp  = perp / np.linalg.norm(perp)
        perp2 = np.cross(eje_n, perp)

        # Radio varía ligeramente a lo largo del segmento
        for j, punto in enumerate([p0, p1]):
            for k in range(lados):
                angulo = 2 * math.pi * k / lados
                offset_r = radio * (math.cos(angulo) * perp +
                                    math.sin(angulo) * perp2)
                vertices.append(punto * 1000 + offset_r * 1000)

        # Caras del cilindro
        base = offset
        for k in range(lados):
            a = base + k
            b = base + k + lados
            c = base + ((k + 1) % lados) + lados
            d = base + ((k + 1) % lados)
            caras.append((a, b, c, d))

        offset += lados * 2

    return vertices, caras

# ── EXPORTAR OBJ ──────────────────────────────────────────────────────────────
def exportar_obj(segmentos, ruta):
    vertices_totales = []
    caras_totales    = []
    v_offset         = 0

    print(f"\n  Procesando {len(segmentos)} segmentos...")

    for i, seg in enumerate(segmentos):
        # Generar curva Bezier para este segmento
        curva = spline_bezier(
            seg['inicio'], seg['fin'],
            curvatura=0.20,
            n=PUNTOS_CURVA)

        # Ajustar radio mínimo visible en Blender
        radio_vis = max(seg['radio'], 0.0003)

        verts, caras = cilindro_curvo(curva, radio_vis, LADOS)

        if not verts:
            continue

        vertices_totales.extend(verts)
        for cara in caras:
            caras_totales.append(
                tuple(c + v_offset for c in cara))
        v_offset += len(verts)

        if (i + 1) % 200 == 0:
            print(f"    → {i+1}/{len(segmentos)} segmentos procesados")

    # Escribir archivo OBJ
    print(f"\n  Escribiendo {len(vertices_totales)} vértices "
          f"y {len(caras_totales)} caras...")

    with open(ruta, 'w') as f:
        f.write("# Bio-Kidney AI 2026 — Árbol Vascular CCO v7\n")
        f.write("# Curvas Bezier anatómicas · Ley de Murray α=3.0\n")
        f.write(f"# Segmentos: {len(segmentos)}\n")
        f.write("# Unidades: mm\n\n")
        f.write("o ArbolVascular_CCO_v7\n\n")

        for v in vertices_totales:
            f.write(f"v {v[0]:.4f} {v[1]:.4f} {v[2]:.4f}\n")

        f.write("\n")
        for cara in caras_totales:
            indices = " ".join(str(c + 1) for c in cara)
            f.write(f"f {indices}\n")

    print(f"\n{'═'*54}")
    print(f"  EXPORTACIÓN COMPLETADA")
    print(f"{'═'*54}")
    print(f"  Segmentos exportados : {len(segmentos)}")
    print(f"  Vértices totales     : {len(vertices_totales)}")
    print(f"  Caras totales        : {len(caras_totales)}")
    print(f"  Archivo              : {ruta}")
    print(f"{'═'*54}")
    print(f"\n  Para abrir en Blender:")
    print(f"  File → Import → Wavefront (.obj)")
    print(f"  Selecciona: arbol_vascular_cco_v7_curvo.obj\n")

# ── MAIN ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("\n" + "═"*54)
    print("  BIO-KIDNEY AI 2026 — EXPORTADOR BLENDER v2")
    print("  Curvas Bezier · Cilindros curvos · 3 sistemas")
    print("═"*54)

    print("\n  Cargando segmentos vasculares...")
    segs = cargar_segmentos(os.path.expanduser(ENTRADA))
    print(f"  {len(segs)} segmentos cargados")

    exportar_obj(segs, os.path.expanduser(SALIDA))

