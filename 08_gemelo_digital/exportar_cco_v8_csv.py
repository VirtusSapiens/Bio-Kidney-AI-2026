#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
exportar_cco_v8_csv.py  -  ADAPTADOR gemelo -> CSV de segmentos

Lee las capas validadas del gemelo digital (.npz) y escribe un CSV de segmentos
con EXACTAMENTE el mismo esquema de columnas que arbol_vascular_cco_v7.csv, de
modo que los simuladores (oxigeno, glomerular) puedan consumir el gemelo nuevo
sin modificar su codigo.

  - NO modifica las capas .npz.
  - NO modifica los simuladores.
  - NO sobrescribe el v7.

Esquema v7 (canonico, 10 columnas):
    id,sistema,nivel,x1_mm,y1_mm,z1_mm,x2_mm,y2_mm,z2_mm,radio_um

Mapeo capa -> sistema:
    capa3a_arterial     -> 'art'
    capa3b_venoso       -> 'ven'
    capa3c_colector     -> 'col'
    capa4_colector_alto -> 'col'

Convenciones verificadas contra el v7 y contra los simuladores:
    - Coordenadas: las capas ya estan en mm  -> se escriben en mm (columnas *_mm),
      SIN conversion. (El simulador de oxigeno hace *0.1 mm->cm al cargar.)
    - Radios: las capas guardan radio en mm -> se convierte a um (x1000) para la
      columna radio_um. (El simulador hace *1e-4 um->cm al cargar.)
    - sistema: se usan los codigos cortos 'art'/'ven'/'col' identicos a los del v7
      (el loader del oxigeno deriva la presion de este string; NO se lee del CSV).
    - nivel: profundidad topologica (saltos) desde la raiz del sistema. Dato real
      derivado del grafo, no inventado. (El loader del oxigeno lo ignora.)

Uso:
    python exportar_cco_v8_csv.py [--salida RUTA_CSV] [--raiz DIR_NPZ]
"""
import argparse
import csv
from collections import deque
from pathlib import Path

import numpy as np

# Ubicaciones por defecto ------------------------------------------------------
DIR_NPZ_DEFAULT = Path(__file__).resolve().parent.parent  # raiz del repo, donde viven los .npz
CSV_SALIDA_DEFAULT = DIR_NPZ_DEFAULT / "02_vascular_cco" / "arbol_vascular_cco_v8_gemelo.csv"

MM_A_UM = 1000.0  # radios de las capas estan en mm -> um

HEADER_V7 = ["id", "sistema", "nivel",
             "x1_mm", "y1_mm", "z1_mm",
             "x2_mm", "y2_mm", "z2_mm", "radio_um"]


def _profundidad_bfs(n_nodos, aristas, raiz_idx):
    """Profundidad (saltos) de cada nodo desde raiz_idx via BFS sobre el grafo
    no dirigido. Nodos inalcanzables -> profundidad 0."""
    ady = [[] for _ in range(n_nodos)]
    for a, b in aristas:
        ady[a].append(b)
        ady[b].append(a)
    prof = [0] * n_nodos
    visto = [False] * n_nodos
    visto[raiz_idx] = True
    q = deque([raiz_idx])
    while q:
        u = q.popleft()
        for v in ady[u]:
            if not visto[v]:
                visto[v] = True
                prof[v] = prof[u] + 1
                q.append(v)
    return prof


def segmentos_arista_por_arista(npz_path, sistema, radio_por_arista=True):
    """Capas con arrays 'nodos','aristas','radios' (arterial, venoso, capa4).
    Devuelve lista de dicts de segmento. radios en mm -> um."""
    d = np.load(npz_path, allow_pickle=True)
    nodos = np.asarray(d["nodos"], dtype=float)      # (N,3) en mm
    aristas = np.asarray(d["aristas"], dtype=int)    # (E,2) indices [padre, hijo]
    radios = np.asarray(d["radios"], dtype=float)    # (E,) mm  o  (N,) mm

    # Determinar raiz para nivel BFS
    if "raiz" in d.files and np.ndim(d["raiz"]) == 0:
        raiz_idx = int(d["raiz"])
    elif "raiz_pos" in d.files:
        raiz_idx = int(np.argmin(np.linalg.norm(nodos - np.asarray(d["raiz_pos"], float), axis=1)))
    elif "raiz" in d.files and np.asarray(d["raiz"]).shape == (3,):
        raiz_idx = int(np.argmin(np.linalg.norm(nodos - np.asarray(d["raiz"], float), axis=1)))
    else:
        raiz_idx = 0
    prof = _profundidad_bfs(len(nodos), aristas, raiz_idx)

    segs = []
    for i, (a, b) in enumerate(aristas):
        p1 = nodos[a]
        p2 = nodos[b]
        if radio_por_arista and len(radios) == len(aristas):
            r_mm = float(radios[i])
        else:  # radios por nodo -> usar nodo hijo
            r_mm = float(radios[b])
        segs.append({
            "sistema": sistema,
            "nivel": int(prof[b]),
            "x1": p1[0], "y1": p1[1], "z1": p1[2],
            "x2": p2[0], "y2": p2[1], "z2": p2[2],
            "radio_um": r_mm * MM_A_UM,
        })
    return segs


def segmentos_parent(npz_path, sistema):
    """Capa colector 3c: arrays 'nodos','parent','radios' (parent-pointer tree).
    Una arista por nodo con parent>=0. radio = radio del nodo hijo (mm->um)."""
    d = np.load(npz_path, allow_pickle=True)
    nodos = np.asarray(d["nodos"], dtype=float)
    parent = np.asarray(d["parent"], dtype=int)
    radios = np.asarray(d["radios"], dtype=float)  # (N,) mm por nodo

    # profundidad via cadena de parents
    prof = np.zeros(len(nodos), dtype=int)
    orden = np.argsort(parent)  # padres antes que hijos aprox; recomputamos robusto
    # calculo robusto de profundidad con memoizacion
    memo = {}

    def depth(i):
        if parent[i] < 0:
            return 0
        if i in memo:
            return memo[i]
        dv = depth(int(parent[i])) + 1
        memo[i] = dv
        return dv

    segs = []
    for i in range(len(nodos)):
        p = int(parent[i])
        if p < 0:
            continue
        p1 = nodos[p]   # padre = inicio
        p2 = nodos[i]   # hijo  = fin
        segs.append({
            "sistema": sistema,
            "nivel": int(depth(i)),
            "x1": p1[0], "y1": p1[1], "z1": p1[2],
            "x2": p2[0], "y2": p2[1], "z2": p2[2],
            "radio_um": float(radios[i]) * MM_A_UM,
        })
    return segs


def exportar(dir_npz, csv_salida):
    dir_npz = Path(dir_npz)
    todos = []
    resumen = {}

    fuentes = [
        ("capa3a_arterial.npz", "art", segmentos_arista_por_arista),
        ("capa3b_venoso.npz",   "ven", segmentos_arista_por_arista),
        ("capa3c_colector.npz", "col", segmentos_parent),
        ("capa4_colector_alto.npz", "col", segmentos_arista_por_arista),
    ]

    for nombre, sistema, fn in fuentes:
        ruta = dir_npz / nombre
        if not ruta.exists():
            print(f"  [AVISO] no existe {ruta} - se omite")
            continue
        segs = fn(ruta, sistema)
        todos.extend(segs)
        resumen[nombre] = (sistema, len(segs))
        print(f"  {nombre:28s} sistema={sistema}  segmentos={len(segs)}")

    # Escritura CSV con esquema v7 exacto
    csv_salida = Path(csv_salida)
    csv_salida.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_salida, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(HEADER_V7)
        for i, s in enumerate(todos):
            w.writerow([
                i, s["sistema"], s["nivel"],
                f"{s['x1']:.4f}", f"{s['y1']:.4f}", f"{s['z1']:.4f}",
                f"{s['x2']:.4f}", f"{s['y2']:.4f}", f"{s['z2']:.4f}",
                f"{s['radio_um']:.2f}",
            ])

    print(f"\n  Escrito: {csv_salida}  ({len(todos)} segmentos)")
    return todos, resumen, csv_salida


def verificar(csv_salida):
    """Recarga el CSV escrito y reporta rangos fisicos y particion por sistema."""
    import pandas as pd
    df = pd.read_csv(csv_salida)
    print("\n  --- VERIFICACION (recarga) ---")
    print(f"  Header leido: {list(df.columns)}")
    print(f"  Filas: {len(df)}")
    print("  Particion por sistema:")
    for s, n in df["sistema"].value_counts().items():
        sub = df[df["sistema"] == s]
        print(f"    {s}: {n}  radio_um [{sub['radio_um'].min():.2f}, {sub['radio_um'].max():.2f}]")
    coords = df[["x1_mm", "y1_mm", "z1_mm", "x2_mm", "y2_mm", "z2_mm"]].values
    print(f"  Coords mm  min={coords.min():.3f}  max={coords.max():.3f}")
    print(f"  radio_um   min={df['radio_um'].min():.2f}  max={df['radio_um'].max():.2f}")
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--raiz", default=str(DIR_NPZ_DEFAULT), help="directorio con los .npz")
    ap.add_argument("--salida", default=str(CSV_SALIDA_DEFAULT), help="ruta CSV de salida")
    args = ap.parse_args()

    print("== EXPORTADOR gemelo -> CSV segmentos (esquema v7) ==")
    todos, resumen, ruta = exportar(args.raiz, args.salida)
    verificar(ruta)
