#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa3b_anastomosis_vv.py  --  Bio-Kidney AI 2026
================================================
AUDITORIA DE SOLO LECTURA. NO escribe .npz, NO promueve, NO toca generador /
reparador / visualizador. Solo el .py y stdout.

Clasifica los contactos VV del cuadrante CRITICO (no-exacta ∩ >=1 interno, los
7106 del venoso) en "ANASTOMOSIS-PLAUSIBLE" vs "ROCE ACCIDENTAL", para decidir
cuales soldar (union fisiologica) y cuales separar (defecto imprimible).

Por cada colision critica:
  - RATIO de calibre = min(r_i,r_j)/max(r_i,r_j). Anastomosis = vasos comparables
    (ratio->1); roce twig-vs-tronco = ratio->0.
  - AMBOS INTERNOS: los dos segmentos con sus dos extremos en grado >=2 (una
    anastomosis une dos vasos CON RECORRIDO, no una hoja contra un vaso).
  - ANGULO entre segmentos (P->Q): cruce/confluencia (alto) vs paralelas que rozan (bajo).
  - PENETRACION p = thr - dd.
  - UBICACION: cortex / union corticomedular (CMJ) / medula profunda (capa0) y banda al hilio.

Criterio (PROVISIONAL, fijado tras ver las distribuciones crudas):
  anastomosis-plausible = AMBOS INTERNOS  &  ratio >= RATIO_MIN  &
                          ubicacion ∈ {cortex, CMJ}  &  angulo >= ANG_CRUCE.
El resto del cuadrante critico = roce accidental. De los accidentales, los que
superan w=50 µm son shunts reales a reparar/regenerar.

Control: mismo conteo sobre el arterial (esperado ~0 anastomosis: las arterias
renales son terminales y no anastomosan; si salen muchas, el criterio esta laxo).

Reutiliza clasificar()/find_npz/HILIO/BANDAS_HILIO (capa3b_clasificar_vv) y
critico() (capa3b_severidad_vv).

Solo NumPy / SciPy. SIN bpy. Entorno: env_biokidney.
"""

import os
import sys
import numpy as np
from scipy.spatial import cKDTree

# ============================================================================
#  PARAMETROS  (provisionales; se justifican con las distribuciones crudas)
# ============================================================================
RATIO_MIN      = 0.5     # calibre comparable: min/max r >= esto
ANG_CRUCE_DEG  = 30.0    # cruce/confluencia: angulo >= esto (no paralelas)
CMJ_BAND_MM    = 2.0     # banda de la union corticomedular alrededor del limite cortex/medula
W_SHUNT_UM     = 50      # grosor de pared: accidental con p > esto = shunt real
SUB_DOM        = 2       # submuestreo de capa0 para los KD-tree de region (1 = sin)

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
from capa3b_severidad_vv import critico


def pct(n, tot):
    return f"{(100.0*n/tot):.1f}%" if tot else "  - "


def deciles(v):
    if len(v) == 0:
        return "(vacio)"
    qs = np.percentile(v, np.arange(0, 101, 10))
    return "  ".join(f"{q:.3f}" for q in qs)


# ============================================================================
#  CLASIFICADOR DE UBICACION cortex / CMJ / medula profunda
# ============================================================================
def construir_region_trees(dom_path):
    d0 = np.load(dom_path, allow_pickle=False)
    coords = d0["coords"].astype(np.float64)[::SUB_DOM]
    region = d0["region_label"][::SUB_DOM]
    es_ctx = region == "cortex"
    tree_ctx = cKDTree(coords[es_ctx])
    tree_med = cKDTree(coords[~es_ctx])      # medulla + piramides = medula
    return tree_ctx, tree_med


def clasificar_ubicacion(cont, tree_ctx, tree_med, band=CMJ_BAND_MM):
    """Devuelve un array de etiquetas: 'cortex' | 'CMJ' | 'medula_prof'.

    CMJ = punto a < band del limite cortex/medula (su region opuesta esta cerca).
    """
    if len(cont) == 0:
        return np.empty(0, dtype="<U12")
    d_ctx, _ = tree_ctx.query(cont, k=1)
    d_med, _ = tree_med.query(cont, k=1)
    en_cortex = d_ctx <= d_med
    d_opuesta = np.where(en_cortex, d_med, d_ctx)     # distancia a la otra region
    loc = np.where(en_cortex, "cortex", "medula_prof").astype("<U12")
    loc[d_opuesta < band] = "CMJ"
    return loc


# ============================================================================
#  ANALISIS POR ARBOL
# ============================================================================
def analizar(npz, tree_ctx, tree_med):
    nodos = npz["nodos"].astype(np.float64)
    arist = npz["aristas"].astype(np.int64)
    radios = npz["radios"].astype(np.float64)
    cl = clasificar(nodos, arist, radios)
    cr = critico(cl)                                   # cuadrante critico
    n = cr["n"]
    if n == 0:
        return dict(n=0, n_seg=cl["n_seg"])

    i, j = cr["i"], cr["j"]
    ri, rj = radios[i], radios[j]
    ratio = np.minimum(ri, rj) / np.maximum(ri, rj)
    p = cr["p"]
    cont = cr["cont"]

    seg_term = cl["seg_terminal"]
    ambos_int = (~seg_term[i]) & (~seg_term[j])        # los dos segmentos internos

    # angulo entre segmentos (no orientado, 0..90)
    vi = nodos[arist[i, 1]] - nodos[arist[i, 0]]
    vj = nodos[arist[j, 1]] - nodos[arist[j, 0]]
    ni = np.linalg.norm(vi, axis=1); nj = np.linalg.norm(vj, axis=1)
    cosang = np.abs(np.einsum("ij,ij->i", vi, vj)) / np.maximum(ni * nj, 1e-12)
    ang = np.degrees(np.arccos(np.clip(cosang, 0.0, 1.0)))

    loc = clasificar_ubicacion(cont, tree_ctx, tree_med)
    dist_h = np.linalg.norm(cont - HILIO, axis=1)

    return dict(n=n, n_seg=cl["n_seg"], i=i, j=j, ratio=ratio, p=p, cont=cont,
                ambos_int=ambos_int, ang=ang, loc=loc, dist_h=dist_h)


def hist2d_ratio_loc(ratio, loc):
    """Matriz (4 bins de ratio) x (cortex, CMJ, medula_prof)."""
    edges = [0.0, 0.25, 0.5, 0.75, 1.0001]
    rb = np.digitize(ratio, edges) - 1
    cats = ["cortex", "CMJ", "medula_prof"]
    M = np.zeros((4, 3), dtype=int)
    for b in range(4):
        for c, cat in enumerate(cats):
            M[b, c] = int(np.count_nonzero((rb == b) & (loc == cat)))
    return M, cats


# ============================================================================
#  MAIN
# ============================================================================
def main():
    art_path = find_npz(IN_ART)
    ven_path = find_npz(IN_VEN)
    dom_path = find_npz(IN_DOM)

    print("=" * 82)
    print("  CAPA 3b - ANASTOMOSIS-PLAUSIBLE vs ROCE ACCIDENTAL (VV critico, solo lectura)")
    print("=" * 82)
    print(f"  Venoso (candidato): {ven_path}")
    print(f"  Arterial (control): {art_path}")
    print(f"  Criterio provisional: ambos-internos & ratio>={RATIO_MIN} & loc∈(cortex,CMJ) "
          f"& ang>={ANG_CRUCE_DEG}deg")
    print("-" * 82)

    tree_ctx, tree_med = construir_region_trees(dom_path)
    ven = analizar(np.load(ven_path, allow_pickle=False), tree_ctx, tree_med)
    art = analizar(np.load(art_path, allow_pickle=False), tree_ctx, tree_med)

    # ------------------------------------------------------------------
    #  DISTRIBUCIONES CRUDAS (sin aplicar cortes todavia) -- VENOSO
    # ------------------------------------------------------------------
    print("  DISTRIBUCIONES CRUDAS del cuadrante critico VENOSO (deciles 0..100%)")
    print("  " + "-" * 78)
    print(f"  ratio calibre min/max : {deciles(ven['ratio'])}")
    print(f"  angulo segmentos (deg): {deciles(ven['ang'])}")
    print(f"  penetracion p (mm)    : {deciles(ven['p'])}")
    print("  " + "-" * 78)
    a_int = int(np.count_nonzero(ven["ambos_int"]))
    print(f"  ambos-internos : {a_int} ({pct(a_int, ven['n'])})   "
          f"solo-uno-interno : {ven['n']-a_int} ({pct(ven['n']-a_int, ven['n'])})  "
          f"de {ven['n']} criticas")
    print("-" * 82)

    # ------------------------------------------------------------------
    #  HISTOGRAMA 2D (ratio x ubicacion) -- VENOSO
    # ------------------------------------------------------------------
    M, cats = hist2d_ratio_loc(ven["ratio"], ven["loc"])
    print("  HISTOGRAMA 2D  (bins de ratio de calibre)  x  (ubicacion)  -- VENOSO")
    print("  " + "-" * 78)
    print(f"  {'ratio \\ loc':<14}{'cortex':>12}{'CMJ':>10}{'medula_prof':>14}{'fila':>9}")
    et = ["[0.00-0.25)", "[0.25-0.50)", "[0.50-0.75)", "[0.75-1.00]"]
    for b in range(4):
        fila = int(M[b].sum())
        print(f"  {et[b]:<14}{M[b,0]:>12d}{M[b,1]:>10d}{M[b,2]:>14d}{fila:>9d}")
    col = M.sum(axis=0)
    print(f"  {'col':<14}{col[0]:>12d}{col[1]:>10d}{col[2]:>14d}{int(M.sum()):>9d}")
    print(f"  -> calibre comparable (ratio>={RATIO_MIN}) en cortex+CMJ: "
          f"{int(M[2:,0:2].sum())}   en medula_prof: {int(M[2:,2].sum())}")
    print("-" * 82)

    # ------------------------------------------------------------------
    #  CLASIFICACION FINAL  -- VENOSO y ARTERIAL (control)
    # ------------------------------------------------------------------
    def clasifica_final(D):
        loc_ok = (D["loc"] == "cortex") | (D["loc"] == "CMJ")
        plausible = (D["ambos_int"] & (D["ratio"] >= RATIO_MIN) &
                     loc_ok & (D["ang"] >= ANG_CRUCE_DEG))
        accidental = ~plausible
        shunt = accidental & (D["p"] > W_SHUNT_UM / 1000.0)
        return plausible, accidental, shunt

    pv, av, sv = clasifica_final(ven)
    npl_v = int(pv.sum()); nac_v = int(av.sum()); nsh_v = int(sv.sum())

    if art["n"] > 0:
        pa, aa, sa = clasifica_final(art)
        npl_a = int(pa.sum()); nac_a = int(aa.sum()); nsh_a = int(sa.sum())
    else:
        npl_a = nac_a = nsh_a = 0

    print("  TABLA FINAL  (de las criticas del cuadrante)")
    print("  " + "-" * 78)
    print(f"  {'metrica':<46}{'VENOSO':>14}{'ARTERIAL':>14}")
    print("  " + "-" * 78)
    print(f"  {'criticas (no-exacta ∩ >=1 interno)':<46}{ven['n']:>14d}{art['n']:>14d}")
    print(f"  {'ANASTOMOSIS-PLAUSIBLE':<46}{npl_v:>14d}{npl_a:>14d}")
    print(f"  {'  (% sobre criticas)':<46}{pct(npl_v,ven['n']):>14}{pct(npl_a,art['n']):>14}")
    print(f"  {'ROCE ACCIDENTAL':<46}{nac_v:>14d}{nac_a:>14d}")
    print(f"  {f'  de ellos shunt (p>{W_SHUNT_UM}um)':<46}{nsh_v:>14d}{nsh_a:>14d}")
    print("  " + "-" * 78)

    # detalle de las plausibles del venoso (para auditar el criterio)
    if npl_v > 0:
        print(f"  plausibles VENOSO: ratio[mediana {np.median(ven['ratio'][pv]):.2f}] "
              f"ang[mediana {np.median(ven['ang'][pv]):.0f}deg] "
              f"p[mediana {np.median(ven['p'][pv]):.3f}mm]")
        loc_pv = ven["loc"][pv]
        print(f"     ubicacion: cortex {int((loc_pv=='cortex').sum())}  "
              f"CMJ {int((loc_pv=='CMJ').sum())}")
    print("-" * 82)

    # ------------------------------------------------------------------
    #  VEREDICTO
    # ------------------------------------------------------------------
    print("  VEREDICTO")
    control_ok = npl_a <= max(5, int(0.01 * max(art["n"], 1)))
    print(f"  Venoso: {npl_v} anastomosis-plausibles y {nac_v} roces accidentales "
          f"(de ellos {nsh_v} shunts > {W_SHUNT_UM}um a separar/regenerar).")
    print(f"  Control arterial: {npl_a} anastomosis-plausibles de {art['n']} criticas  "
          f"-> {'OK criterio (arterial ~0, no anastomosa)' if control_ok else 'CRITERIO LAXO: aprietalo (arterial deberia dar ~0)'}")
    print("=" * 82)


if __name__ == "__main__":
    main()
