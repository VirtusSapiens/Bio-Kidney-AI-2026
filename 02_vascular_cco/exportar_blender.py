import csv
import math
import os
import numpy as np

# ── CONFIGURACIÓN ──────────────────────────────────────────────────────────────
ENTRADA  = os.path.expanduser(
    "~/Escritorio/BioKidney-AI/02_vascular_cco/arbol_vascular_cco_v7.csv")
SALIDA   = os.path.expanduser(
    "~/Escritorio/BioKidney-AI/02_vascular_cco/arbol_vascular_cco_v7.obj")
LADOS    = 8      # caras por cilindro (8 = octágono, suficiente para Blender)
ESCALA   = 1.0    # 1 mm en CSV = 1 unidad Blender

def cargar_segmentos(ruta):
    segmentos = []
    with open(ruta, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            segmentos.append({
                'inicio': np.array([float(row['x1_mm']) / 1000,
                                    float(row['y1_mm']) / 1000,
                                    float(row['z1_mm']) / 1000]),
                'fin'   : np.array([float(row['x2_mm']) / 1000,
                                    float(row['y2_mm']) / 1000,
                                    float(row['z2_mm']) / 1000]),
                'radio' : float(row['radio_um']) / 1000000.0,
            })
    return segmentos

def base_cilindro(inicio, fin, radio, lados):
    """
    Genera vértices e índices de un cilindro entre dos puntos.
    Retorna (vertices, caras) donde caras son índices locales.
    """
    eje = fin - inicio
    long = np.linalg.norm(eje)
    if long < 1e-10:
        return [], []
    eje_n = eje / long

    # Vector perpendicular
    perp = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(perp, eje_n)) > 0.9:
        perp = np.array([0.0, 1.0, 0.0])
    perp = perp - np.dot(perp, eje_n) * eje_n
    perp = perp / np.linalg.norm(perp)
    perp2 = np.cross(eje_n, perp)

    vertices = []
    for i in range(lados):
        angulo = 2 * math.pi * i / lados
        offset = radio * (math.cos(angulo) * perp +
                          math.sin(angulo) * perp2)
        vertices.append(inicio * ESCALA + offset)
        vertices.append(fin    * ESCALA + offset)

    caras = []
    for i in range(lados):
        a = i * 2
        b = i * 2 + 1
        c = ((i + 1) % lados) * 2 + 1
        d = ((i + 1) % lados) * 2
        caras.append((a, b, c, d))

    return vertices, caras

def exportar_obj(segmentos, ruta):
    vertices_totales = []
    caras_totales    = []
    offset           = 0

    for seg in segmentos:
        verts, caras = base_cilindro(
            seg['inicio'], seg['fin'], seg['radio'], LADOS)
        if not verts:
            continue
        vertices_totales.extend(verts)
        for cara in caras:
            caras_totales.append(tuple(c + offset for c in cara))
        offset += len(verts)

    with open(ruta, 'w') as f:
        f.write("# Bio-Kidney AI 2026 — Árbol Vascular CCO\n")
        f.write("# Generado con Ley de Murray (α=3.0)\n")
        f.write(f"# Segmentos: {len(segmentos)}\n")
        f.write(f"# Unidades: mm\n\n")
        f.write("o ArbolVascular_CCO\n\n")

        for v in vertices_totales:
            f.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")

        f.write("\n")
        for cara in caras_totales:
            # OBJ usa índices base 1
            indices = " ".join(str(c + 1) for c in cara)
            f.write(f"f {indices}\n")

    print(f"\n{'═'*50}")
    print(f"  EXPORTACIÓN A BLENDER COMPLETADA")
    print(f"{'═'*50}")
    print(f"  Segmentos exportados : {len(segmentos)}")
    print(f"  Vértices totales     : {len(vertices_totales)}")
    print(f"  Caras totales        : {len(caras_totales)}")
    print(f"  Archivo              : {ruta}")
    print(f"{'═'*50}\n")
    print("  Para abrir en Blender:")
    print("  File → Import → Wavefront (.obj)")
    print(f"  Selecciona: arbol_vascular_cco_v7.obj\n")

if __name__ == "__main__":
    print("\n  Cargando segmentos vasculares...")
    segmentos = cargar_segmentos(os.path.expanduser(ENTRADA))
    print(f"  Segmentos cargados: {len(segmentos)}")
    print("  Generando cilindros 3D...")
    exportar_obj(segmentos, os.path.expanduser(SALIDA))
