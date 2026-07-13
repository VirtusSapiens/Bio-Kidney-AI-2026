#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa3b_severidad_vv.py  --  Bio-Kidney AI 2026
==============================================
AUDITORIA DE SOLO LECTURA. NO escribe .npz, NO promueve, NO toca el reparador ni
el visualizador. Solo extiende el analisis e imprime a stdout.

Barrido de SEVERIDAD sobre el cuadrante critico VV (no-exacta ∩ >=1 interno) del
venoso reparado, con el arterial canonico como control.

Para cada colision critica: penetracion p = thr - dd (thr = r_i + r_j).
  - Barre grosores de pared imprimible w ∈ {20,30,50,75,100,150} µm y cuenta
    cuantas criticas tienen p > w (candidatas a FUNDIR dos lumenes a esa resolucion
    = "shunt").
  - Cola severa (p > COLA_MM): nº de segmentos distintos implicados (concentracion
    vs dispersion), grado/profundidad de los nodos, ubicacion cortex/medula y banda
    al hilio, y los top-20 segmentos por nº de colisiones (reparacion dirigida).
  - Control: mismo barrido sobre el arterial; colas lado a lado.

Reutiliza clasificar()/find_npz/HILIO/BANDAS_HILIO de capa3b_clasificar_vv.

Solo NumPy / SciPy. SIN bpy. Entorno: env_biokidney.
"""

import os
import sys
import numpy as np
from scipy.spatial import cKDTree

# ============================================================================
#  PARAMETROS
# ============================================================================
W_UM = [20, 30, 50, 75, 100, 150]    # grosores de pared imprimible (micras)
COLA_MM = 0.20                        # umbral de "cola severa": p > 0.20 mm
W_PLACEHOLDER_UM = 50                 # grosor de referencia para el veredicto
TOP_N = 20                            # segmentos mas problematicos a listar

IN_ART = "capa3a_arterial.npz"
IN_VEN = "capa3b_venoso_reparado.npz"
IN_DOM = "capa0_dominio.npz"

# ============================================================================
#  REUTILIZACION
# ============================================================================
_here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
if _here not in sys.path:
    sys.path.insert(0, _here)
from capa3b_clasificar_vv import clasificar, find_npz, HILIO, BANDAS_HILIO


def pct(n, tot):
    return f"{(100.0*n/tot):.1f}%" if tot else "  - "


def profundidades(arist, n_nodos):
    """Profundidad (n.º de aristas desde la raiz) por nodo. parent[i] < i."""
    parent = -np.ones(n_nodos, dtype=np.int64)
    parent[arist[:, 1]] = arist[:, 0]
    depth = np.zeros(n_nodos, dtype=np.int64)
    for i in range(1, n_nodos):
        if parent[i] >= 0:
            depth[i] = depth[parent[i]] + 1
    return depth


def critico(cl):
    """Extrae el cuadrante critico (no-exacta ∩ >=1 interno) de un dict clasificar()."""
    if cl["n_col"] == 0:
        z = np.empty(0)
        return dict(i=z.astype(int), j=z.astype(int), p=z, cont=np.empty((0, 3)), n=0)
    m = (~cl["exacta"]) & (~cl["ambos_term"])
    return dict(i=cl["i"][m], j=cl["j"][m],
                p=(cl["thr"] - cl["dd"])[m], cont=cl["cont"][m],
                n=int(np.count_nonzero(m)))


# ============================================================================
#  MAIN
# ============================================================================
def main():
    art_path = find_npz(IN_ART)
    ven_path = find_npz(IN_VEN)
    dom_path = find_npz(IN_DOM)

    art = np.load(art_path, allow_pickle=False)
    ven = np.load(ven_path, allow_pickle=False)

    print("=" * 80)
    print("  CAPA 3b - BARRIDO DE SEVERIDAD DEL CUADRANTE CRITICO VV (solo lectura)")
    print("=" * 80)
    print(f"  Venoso (candidato) : {ven_path}")
    print(f"  Arterial (control) : {art_path}")
    print(f"  p = thr - dd   thr = r_i + r_j   cola severa: p > {COLA_MM} mm")
    print("-" * 80)

    datos = {}
    for nombre, npz in (("VENOSO", ven), ("ARTERIAL", art)):
        nodos = npz["nodos"].astype(np.float64)
        arist = npz["aristas"].astype(np.int64)
        radios = npz["radios"].astype(np.float64)
        cl = clasificar(nodos, arist, radios)
        cr = critico(cl)
        datos[nombre] = dict(nodos=nodos, arist=arist, radios=radios, cl=cl, cr=cr,
                             grado=cl["grado"], n_seg=cl["n_seg"])

    n_crit_ven = datos["VENOSO"]["cr"]["n"]
    n_crit_art = datos["ARTERIAL"]["cr"]["n"]

    # ------------------------------------------------------------------
    #  BARRIDO DE GROSOR DE PARED (shunts = colisiones criticas con p > w)
    # ------------------------------------------------------------------
    print("  BARRIDO DE GROSOR DE PARED w  (shunts = criticas con p > w)")
    print("  " + "-" * 76)
    print(f"  {'w (um)':>7} | {'VENOSO shunts':>22} | {'ARTERIAL shunts':>24}")
    print(f"  {'':>7} | {'n   %crit   %seg':>22} | {'n   %crit   %seg':>24}")
    print("  " + "-" * 76)
    for w in W_UM:
        wm = w / 1000.0
        nv = int(np.count_nonzero(datos["VENOSO"]["cr"]["p"] > wm))
        na = int(np.count_nonzero(datos["ARTERIAL"]["cr"]["p"] > wm))
        cv = f"{nv:>5d} {pct(nv,n_crit_ven):>6} {pct(nv,datos['VENOSO']['n_seg']):>6}"
        ca = f"{na:>5d} {pct(na,n_crit_art):>6} {pct(na,datos['ARTERIAL']['n_seg']):>6}"
        print(f"  {w:>7d} | {cv:>22} | {ca:>24}")
    print("  " + "-" * 76)
    print(f"  {'(total critico)':>7} | {n_crit_ven:>5d}{'':>17}| {n_crit_art:>5d}")
    print("-" * 80)

    # capa0 para ubicacion (una sola vez)
    d0 = np.load(dom_path, allow_pickle=False)
    coords0 = d0["coords"].astype(np.float64)
    region0 = d0["region_label"]
    tree0 = cKDTree(coords0)

    # ------------------------------------------------------------------
    #  COLA SEVERA  (p > COLA_MM)  -- detalle por arbol
    # ------------------------------------------------------------------
    resumen_cola = {}
    for nombre in ("VENOSO", "ARTERIAL"):
        D = datos[nombre]; cr = D["cr"]
        mcola = cr["p"] > COLA_MM
        ic = cr["i"][mcola]; jc = cr["j"][mcola]
        pc = cr["p"][mcola]; cc = cr["cont"][mcola]
        n_cola = int(mcola.sum())
        resumen_cola[nombre] = dict(n=n_cola, i=ic, j=jc, p=pc, cont=cc)

        print(f"  COLA SEVERA (p > {COLA_MM} mm) -- {nombre}")
        print("  " + "-" * 76)
        if n_cola == 0:
            print("  (cola vacia)")
            print("-" * 80)
            continue

        segs_unicos = np.unique(np.concatenate([ic, jc]))
        print(f"  colisiones en la cola : {n_cola}   "
              f"({pct(n_cola, cr['n'])} de las criticas, {pct(n_cola, D['n_seg'])} de segmentos)")
        print(f"  segmentos DISTINTOS implicados : {len(segs_unicos)}   "
              f"-> {n_cola/max(len(segs_unicos),1):.2f} colisiones por segmento "
              f"({'concentrado (plegado local)' if n_cola > 1.5*len(segs_unicos) else 'disperso'})")

        # grado y profundidad de los nodos implicados
        arist = D["arist"]; grado = D["grado"]
        depth = profundidades(arist, len(D["nodos"]))
        nodos_impl = np.unique(arist[np.concatenate([ic, jc])].ravel())
        gr = grado[nodos_impl]
        dp = depth[nodos_impl]
        dmax = int(depth.max())
        hist_g = {int(g): int(np.count_nonzero(gr == g)) for g in np.unique(gr)}
        print(f"  grado de nodos implicados : {hist_g}   "
              f"(1=hoja, 2=pasante, 3=bifurcacion)")
        print(f"  profundidad nodos implicados (de {dmax} max): "
              f"min {dp.min()}  mediana {int(np.median(dp))}  p90 {int(np.percentile(dp,90))}  max {dp.max()}")
        tercio = dmax / 3.0
        n_sup = int(np.count_nonzero(dp < tercio))
        n_med = int(np.count_nonzero((dp >= tercio) & (dp < 2 * tercio)))
        n_prof = int(np.count_nonzero(dp >= 2 * tercio))
        print(f"     tercio troncal {n_sup} ({pct(n_sup,len(nodos_impl))})  "
              f"medio {n_med} ({pct(n_med,len(nodos_impl))})  "
              f"distal {n_prof} ({pct(n_prof,len(nodos_impl))})")

        # ubicacion cortex/medula + banda al hilio
        _, idx0 = tree0.query(cc, k=1)
        es_ctx = region0[idx0] == "cortex"
        n_ctx = int(es_ctx.sum())
        print(f"  ubicacion: cortex {n_ctx} ({pct(n_ctx,n_cola)})   "
              f"medula {n_cola-n_ctx} ({pct(n_cola-n_ctx,n_cola)})")
        dist_h = np.linalg.norm(cc - HILIO, axis=1)
        bandas = []
        for lo, hi in BANDAS_HILIO:
            mb = (dist_h >= lo) & (dist_h < hi)
            etiq = f"[{lo:.0f},{'inf' if np.isinf(hi) else f'{hi:.0f}'})"
            bandas.append(f"{etiq}:{int(mb.sum())}")
        print(f"  banda al hilio (mm): " + "  ".join(bandas))

        # top-N segmentos por nº de colisiones de la cola
        todos = np.concatenate([ic, jc])
        segs, cuenta = np.unique(todos, return_counts=True)
        orden = np.argsort(-cuenta)
        seg_term = D["cl"]["seg_terminal"]
        print(f"  top {TOP_N} segmentos por nº de colisiones en la cola "
              f"(reparar estos limpia mas cruces):")
        topk = orden[:TOP_N]
        partes = []
        for k in topk:
            s = int(segs[k]); c = int(cuenta[k])
            tipo = "T" if seg_term[s] else "I"
            partes.append(f"#{s}({tipo}):{c}")
        # imprime en filas de 5
        for r0 in range(0, len(partes), 5):
            print("     " + "  ".join(partes[r0:r0 + 5]))
        # cobertura del top-N: cuantas colisiones de la cola tocan al menos uno
        top_set = set(int(segs[k]) for k in topk)
        cubiertas = int(np.count_nonzero([(a in top_set) or (b in top_set)
                                          for a, b in zip(ic, jc)]))
        print(f"  cobertura: reparar esos {len(top_set)} segmentos resuelve "
              f"{cubiertas}/{n_cola} colisiones de la cola ({pct(cubiertas,n_cola)})")
        print("-" * 80)

    # ------------------------------------------------------------------
    #  TABLA COMPACTA de colas lado a lado + VEREDICTO
    # ------------------------------------------------------------------
    print("  COLAS LADO A LADO")
    print("  " + "-" * 76)
    print(f"  {'metrica':<40}{'VENOSO':>16}{'ARTERIAL':>16}")
    print("  " + "-" * 76)
    rv, ra = resumen_cola["VENOSO"], resumen_cola["ARTERIAL"]
    sv = len(np.unique(np.concatenate([rv['i'], rv['j']]))) if rv['n'] else 0
    sa = len(np.unique(np.concatenate([ra['i'], ra['j']]))) if ra['n'] else 0
    print(f"  {'criticas totales':<40}{n_crit_ven:>16d}{n_crit_art:>16d}")
    print(f"  {f'cola severa (p>{COLA_MM}mm)':<40}{rv['n']:>16d}{ra['n']:>16d}")
    print(f"  {'segmentos distintos en la cola':<40}{sv:>16d}{sa:>16d}")
    wm = W_PLACEHOLDER_UM / 1000.0
    shunt_v = int(np.count_nonzero(datos['VENOSO']['cr']['p'] > wm))
    shunt_a = int(np.count_nonzero(datos['ARTERIAL']['cr']['p'] > wm))
    print(f"  {f'shunts a w={W_PLACEHOLDER_UM}um (p>w)':<40}{shunt_v:>16d}{shunt_a:>16d}")
    print("  " + "-" * 76)
    print("-" * 80)
    print("  VEREDICTO")
    print(f"  A un grosor de pared de {W_PLACEHOLDER_UM} um (placeholder): "
          f"{shunt_v} shunts reales en el VENOSO vs {shunt_a} en el ARTERIAL.")
    print("=" * 80)


if __name__ == "__main__":
    main()
