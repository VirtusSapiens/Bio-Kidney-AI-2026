#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa3b_clasificar_vv.py  --  Bio-Kidney AI 2026
===============================================
AUDITORIA DE SOLO LECTURA. NO escribe ningun .npz, NO promueve nada, NO toca el
visualizador. Solo imprime a stdout.

Clasifica las AUTO-COLISIONES venoso-venoso (VV) de un arbol en una matriz 2x2
para decidir si BLOQUEAN el cierre de la Capa 3b. Se corre sobre dos arboles:
  - capa3a_arterial.npz          (referencia ya canonica)
  - capa3b_venoso_reparado.npz   (candidato a promover)

Una colision VV = dos segmentos NO ADYACENTES (que no comparten extremo) cuya
distancia segmento-segmento dd < thr, con thr = r_i + r_j.

Matriz 2x2:
  eje A: EXACTA (dd < 1e-4)  vs  NO-EXACTA (1e-4 <= dd < thr)
  eje B: AMBOS TERMINALES    vs  >=1 INTERNO
  (segmento TERMINAL = al menos un extremo de grado 1; INTERNO = ambos grado >=2)

El cuadrante critico es NO-EXACTA ∩ >=1 INTERNO: cruces reales entre vasos con
recorrido, que si serian un defecto imprimible. El cuadrante EXACTA ∩ AMBOS
TERMINALES suele ser convergencia benigna (dos hojas que cayeron en el mismo
atractor compartido).

Reutiliza seg_seg_dist (y find_npz) de capa3b_reparar_colisiones.

Solo NumPy / SciPy. SIN bpy. Entorno: env_biokidney.
"""

import os
import sys
import numpy as np
from scipy.spatial import cKDTree

# ============================================================================
#  PARAMETROS
# ============================================================================
UMBRAL_EXACTA = 1e-4     # mm  dd por debajo de esto -> colision "exacta" (coincidente)
UMBRAL_CAND   = 3.0      # mm  candidatos: pares con mid-mid < esto (cubre cualquier
                         #     colision: thr_max ~ 2*r_max ~ 1.1 mm + longitud de seg)
CHUNK         = 2_000_000
HILIO = np.array([0.0, -30.0, 0.0])
BANDAS_HILIO = [(0.0, 5.0), (5.0, 15.0), (15.0, 30.0), (30.0, np.inf)]

IN_ART = "capa3a_arterial.npz"
IN_VEN = "capa3b_venoso_reparado.npz"
IN_DOM = "capa0_dominio.npz"

# ============================================================================
#  REUTILIZACION de la maquinaria de reparacion/auditoria
# ============================================================================
_here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
if _here not in sys.path:
    sys.path.insert(0, _here)
from capa3b_reparar_colisiones import seg_seg_dist, find_npz   # reuso de la logica


# ============================================================================
#  DETECCION DE COLISIONES VV NO ADYACENTES
# ============================================================================
def detectar_vv(nodos, arist, radios):
    """Devuelve (i, j, dd, contacto) de las colisiones VV no adyacentes (dd < r_i+r_j).

    i, j : indices de segmento (no adyacentes, i<j)
    dd   : distancia segmento-segmento
    contacto : (n,3) punto medio del par mas cercano (para ubicar la colision)
    """
    P = nodos[arist[:, 0]]; Q = nodos[arist[:, 1]]
    mid = 0.5 * (P + Q)
    tree = cKDTree(mid)
    pairs = tree.query_pairs(r=UMBRAL_CAND, output_type='ndarray')   # (M,2), i<j
    if len(pairs) == 0:
        return (np.empty(0, np.int64), np.empty(0, np.int64),
                np.empty(0), np.empty((0, 3)))
    i, j = pairs[:, 0], pairs[:, 1]
    # excluir ADYACENTES (comparten algun nodo extremo)
    comparten = (
        (arist[i, 0] == arist[j, 0]) | (arist[i, 0] == arist[j, 1]) |
        (arist[i, 1] == arist[j, 0]) | (arist[i, 1] == arist[j, 1])
    )
    i = i[~comparten]; j = j[~comparten]
    ci, cj, cd, cc = [], [], [], []
    for s0 in range(0, len(i), CHUNK):
        s1 = min(s0 + CHUNK, len(i))
        ii = i[s0:s1]; jj = j[s0:s1]
        dd, cont = seg_seg_dist(P[ii], Q[ii], P[jj], Q[jj])
        thr = radios[ii] + radios[jj]
        hit = dd < thr
        if np.any(hit):
            ci.append(ii[hit]); cj.append(jj[hit]); cd.append(dd[hit]); cc.append(cont[hit])
    if ci:
        return (np.concatenate(ci), np.concatenate(cj),
                np.concatenate(cd), np.concatenate(cc))
    return (np.empty(0, np.int64), np.empty(0, np.int64), np.empty(0), np.empty((0, 3)))


def clasificar(nodos, arist, radios):
    """Detecta y clasifica en la 2x2. Devuelve un dict con todo lo necesario."""
    i, j, dd, cont = detectar_vv(nodos, arist, radios)
    n_nodos = len(nodos)
    n_seg = len(arist)
    grado = np.bincount(arist.ravel(), minlength=n_nodos)
    # segmento terminal si AL MENOS un extremo tiene grado 1
    seg_terminal = (grado[arist[:, 0]] == 1) | (grado[arist[:, 1]] == 1)

    n_col = len(dd)
    if n_col == 0:
        return dict(n_seg=n_seg, n_col=0, grado=grado, seg_terminal=seg_terminal,
                    i=i, j=j, dd=dd, cont=cont, thr=np.empty(0),
                    exacta=np.empty(0, bool), ambos_term=np.empty(0, bool))

    thr = radios[i] + radios[j]
    exacta = dd < UMBRAL_EXACTA
    ambos_term = seg_terminal[i] & seg_terminal[j]        # eje B: ambos terminales
    return dict(n_seg=n_seg, n_col=n_col, grado=grado, seg_terminal=seg_terminal,
                i=i, j=j, dd=dd, cont=cont, thr=thr,
                exacta=exacta, ambos_term=ambos_term)


def matriz_2x2(cl):
    """Devuelve los 4 conteos en orden:
       (exacta&ambos_term, exacta&>=1int, noexacta&ambos_term, noexacta&>=1int)."""
    ex, bt = cl["exacta"], cl["ambos_term"]
    q_ex_bt = int(np.count_nonzero(ex & bt))
    q_ex_in = int(np.count_nonzero(ex & ~bt))
    q_ne_bt = int(np.count_nonzero(~ex & bt))
    q_ne_in = int(np.count_nonzero(~ex & ~bt))
    return q_ex_bt, q_ex_in, q_ne_bt, q_ne_in


def leaf_node(e, arist, grado):
    """Nodo de grado 1 del segmento terminal e (si ambos, el primero)."""
    a, b = int(arist[e, 0]), int(arist[e, 1])
    return a if grado[a] == 1 else b


# ============================================================================
#  MAIN
# ============================================================================
def main():
    art_path = find_npz(IN_ART)
    ven_path = find_npz(IN_VEN)
    dom_path = find_npz(IN_DOM)

    art = np.load(art_path, allow_pickle=False)
    ven = np.load(ven_path, allow_pickle=False)

    print("=" * 78)
    print("  CAPA 3b - CLASIFICACION 2x2 DE AUTO-COLISIONES VENOSO-VENOSO (solo lectura)")
    print("=" * 78)
    print(f"  Arterial (referencia): {art_path}")
    print(f"  Venoso  (candidato)  : {ven_path}")
    print(f"  thr = r_i + r_j      umbral exacta < {UMBRAL_EXACTA}   candidatos mid<{UMBRAL_CAND} mm")
    print("-" * 78)

    cl_art = clasificar(art["nodos"].astype(np.float64),
                        art["aristas"].astype(np.int64),
                        art["radios"].astype(np.float64))
    cl_ven = clasificar(ven["nodos"].astype(np.float64),
                        ven["aristas"].astype(np.int64),
                        ven["radios"].astype(np.float64))

    q_art = matriz_2x2(cl_art)
    q_ven = matriz_2x2(cl_ven)

    def pct(n, tot):
        return f"{(100.0*n/tot):.1f}%" if tot else "  - "

    # --- tabla compacta arterial vs venoso ----------------------------------
    print("  MATRIZ 2x2  (conteo | % sobre total de colisiones del arbol)")
    print("  " + "-" * 74)
    print(f"  {'cuadrante':<34}{'ARTERIAL':>18}{'VENOSO':>18}")
    print("  " + "-" * 74)
    filas = [
        ("EXACTA    & ambos terminales", q_art[0], q_ven[0]),
        ("EXACTA    & >=1 interno",      q_art[1], q_ven[1]),
        ("NO-EXACTA & ambos terminales", q_art[2], q_ven[2]),
        ("NO-EXACTA & >=1 interno  <== CRITICO", q_art[3], q_ven[3]),
    ]
    for nombre, na, nv in filas:
        ca = f"{na:>6d} {pct(na, cl_art['n_col']):>6}"
        cv = f"{nv:>6d} {pct(nv, cl_ven['n_col']):>6}"
        print(f"  {nombre:<34}{ca:>18}{cv:>18}")
    print("  " + "-" * 74)
    print(f"  {'TOTAL colisiones VV':<34}{cl_art['n_col']:>18d}{cl_ven['n_col']:>18d}")
    print(f"  {'(segmentos en el arbol)':<34}{cl_art['n_seg']:>18d}{cl_ven['n_seg']:>18d}")
    print("-" * 78)

    # ========================================================================
    #  CUADRANTE CRITICO (NO-EXACTA ∩ >=1 INTERNO) -- solo VENOSO
    # ========================================================================
    print("  CUADRANTE CRITICO (NO-EXACTA ∩ >=1 INTERNO) -- detalle VENOSO")
    print("  " + "-" * 74)
    crit = (~cl_ven["exacta"]) & (~cl_ven["ambos_term"])
    n_crit = int(np.count_nonzero(crit))
    n_seg_ven = cl_ven["n_seg"]
    print(f"  conteo: {n_crit}   ({100.0*n_crit/n_seg_ven:.3f} % sobre {n_seg_ven} segmentos venosos)")

    if n_crit > 0:
        dd_c = cl_ven["dd"][crit]
        thr_c = cl_ven["thr"][crit]
        pen = thr_c - dd_c                          # profundidad de penetracion
        cont_c = cl_ven["cont"][crit]
        print(f"  profundidad de penetracion thr-dd (mm): "
              f"min {pen.min():.4f}  mediana {np.median(pen):.4f}  "
              f"p90 {np.percentile(pen,90):.4f}  max {pen.max():.4f}")

        # ubicacion: corteza vs medula (capa0, vecino mas cercano del contacto)
        d0 = np.load(dom_path, allow_pickle=False)
        coords = d0["coords"].astype(np.float64)
        region = d0["region_label"]
        tree0 = cKDTree(coords)
        _, idx0 = tree0.query(cont_c, k=1)
        es_cortex = region[idx0] == "cortex"
        n_ctx = int(np.count_nonzero(es_cortex)); n_med = n_crit - n_ctx
        print(f"  ubicacion: cortex {n_ctx} ({pct(n_ctx,n_crit)})   "
              f"medula {n_med} ({pct(n_med,n_crit)})")

        # banda de distancia al hilio
        dist_hil = np.linalg.norm(cont_c - HILIO, axis=1)
        print(f"  banda distancia al hilio (hilio={HILIO.tolist()}):")
        for lo, hi in BANDAS_HILIO:
            m = (dist_hil >= lo) & (dist_hil < hi)
            etiq = f"[{lo:.0f},{'inf' if np.isinf(hi) else f'{hi:.0f}'})"
            n_b = int(np.count_nonzero(m))
            print(f"     {etiq:>10} mm : {n_b:>4d}  ({pct(n_b,n_crit)})")
    else:
        print("  (cuadrante critico vacio en el venoso)")
    print("-" * 78)

    # ========================================================================
    #  CUADRANTE EXACTA ∩ AMBOS TERMINALES -- confirmacion de convergencia benigna
    # ========================================================================
    print("  CUADRANTE EXACTA ∩ AMBOS TERMINALES -- confirmacion (VENOSO)")
    print("  " + "-" * 74)
    ben = cl_ven["exacta"] & cl_ven["ambos_term"]
    n_ben = int(np.count_nonzero(ben))
    if n_ben > 0:
        arist_v = ven["aristas"].astype(np.int64)
        nodos_v = ven["nodos"].astype(np.float64)
        grado_v = cl_ven["grado"]
        ii = cl_ven["i"][ben]; jj = cl_ven["j"][ben]
        coincide = 0
        for a_seg, b_seg in zip(ii, jj):
            la = leaf_node(a_seg, arist_v, grado_v)
            lb = leaf_node(b_seg, arist_v, grado_v)
            if np.linalg.norm(nodos_v[la] - nodos_v[lb]) < UMBRAL_EXACTA:
                coincide += 1
        print(f"  total exacta∩ambos-terminales: {n_ben}")
        print(f"  con los DOS extremos hoja (grado 1) COINCIDENTES: {coincide}  "
              f"({pct(coincide,n_ben)})  -> convergencia genuina por atractor compartido (benigna)")
    else:
        print("  (cuadrante exacta∩ambos-terminales vacio en el venoso)")
    print("-" * 78)

    # ========================================================================
    #  VEREDICTO
    # ========================================================================
    print("  VEREDICTO")
    print("  " + "-" * 74)
    print(f"  Cuadrante CRITICO del venoso (NO-EXACTA ∩ >=1 interno): {n_crit} colisiones "
          f"de {cl_ven['n_col']} VV totales  ({pct(n_crit, cl_ven['n_col'])}).")
    if n_crit == 0:
        print("  -> Ninguna colision critica: las VV son convergencias/tangencias benignas; "
              "no bloquean el cierre de Capa 3b.")
    else:
        print(f"  -> Hay {n_crit} cruce(s) real(es) entre vasos con recorrido: REVISAR antes de "
              "promover (estos si serian defecto imprimible).")
    print("=" * 78)


if __name__ == "__main__":
    main()
