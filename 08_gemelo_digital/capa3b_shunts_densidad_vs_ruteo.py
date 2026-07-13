#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa3b_shunts_densidad_vs_ruteo.py  --  Bio-Kidney AI 2026
==========================================================
AUDITORIA DE SOLO LECTURA. NO escribe ningun .npz, NO toca el generador ni el
visualizador. Solo imprime a stdout.

PREGUNTA: el empaquetamiento INTRA-LOBULAR de los shunts VV de calibre que funden
a 100um, ¿es DENSIDAD GENUINA (limite de escala: no cabe mas) o RUTEO INEFICIENTE
(defecto corregible: habia hueco y la generacion no los separo)?

Trabaja sobre el canonico capa3b_venoso.npz, subconjunto = shunts de calibre que
funden a 100um (p = (r_i+r_j)-dd >= 0.100 mm  y  >=1 segmento de calibre), que ya
sabemos 99.5% intra-zona.

TRES DISCRIMINADORES
--------------------
1) CORRELACION CARGA-vs-SHUNTS por lobulo. Para cada uno de los 10 lobulos:
   (a) nº de vasos de calibre asignados al lobulo (carga), (b) nº de shunts de
   calibre intra-lobulares en el lobulo. Pearson r + dispersion del ratio
   shunts/vaso (coef. de variacion). r alto y CV bajo -> densidad; CV alto
   (lobulos con carga similar y conteos muy distintos) -> ruteo.
2) ESPACIO LIBRE local. Ocupacion = nº de vasos (midpoints de segmento) en una
   esfera R=1mm alrededor del cruce, COMPARADA con la ocupacion de fondo
   (mismo R alrededor de puntos de vaso aleatorios). Y d_next = distancia al
   vaso mas cercano que NO es del par. Si la ocupacion del cruce ~ fondo y
   d_next deja hueco para separarse -> ruteo; si los cruces viven en los puntos
   mas densos y sin hueco -> densidad.
3) PARALELISMO del par. Angulo entre las direcciones de los dos segmentos del
   par (0=paralelos, 90=perpendiculares). Casi-paralelos = vasos que corren
   juntos y se rozan (ruteo); angulo alto = cruces por destinos opuestos
   (inevitables).

Zona de drenaje = lobulo: piramide cuyo eje (apice->base) esta mas cerca del
punto medio del vaso. Reutiliza clasificar()/critico() (cuadrante critico VV).

Solo NumPy / SciPy. SIN bpy. Entorno: env_biokidney.
"""

import os
import sys
import numpy as np
from scipy.spatial import cKDTree

# ============================================================================
#  PARAMETROS
# ============================================================================
CALIBRE_FACTOR = 2.0     # "de calibre": radio > 2x terminal
V_MM           = 0.100   # 100 um: umbral de "funde lumenes"
R_LOCAL        = 1.0     # mm  radio de la esfera de ocupacion local
N_FONDO        = 4000    # nº de puntos aleatorios para la ocupacion de fondo
SEED           = 2026
GAP_FACTOR     = 2.0     # "hay hueco" si d_next > GAP_FACTOR * p (sitio para separar)
PARALELO_DEG   = 20.0    # umbral de "casi paralelo"

IN_VEN = "capa3b_venoso.npz"   # CANONICO (solo lectura)
IN_DOM = "capa0_dominio.npz"

_here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
if _here not in sys.path:
    sys.path.insert(0, _here)
from capa3b_clasificar_vv import clasificar, find_npz   # usa seg_seg_dist por dentro
from capa3b_severidad_vv import critico                  # cuadrante critico VV


def lobulo_de(P, apex, base):
    """Indice de piramide cuyo eje (segmento apice-base) esta mas cerca de cada P."""
    n = len(P); NP = len(apex); dist = np.empty((n, NP))
    for k in range(NP):
        A = apex[k]; B = base[k]; AB = B - A; den = float(AB @ AB)
        t = np.clip(((P - A) @ AB) / den, 0.0, 1.0)
        proj = A + t[:, None] * AB
        dist[:, k] = np.linalg.norm(P - proj, axis=1)
    return dist.argmin(axis=1)


def pctl(x, q):
    return np.percentile(x, q) if len(x) else float('nan')


def main():
    rng = np.random.default_rng(SEED)

    ven = np.load(find_npz(IN_VEN), allow_pickle=False)
    nodos = ven["nodos"].astype(np.float64)
    arist = ven["aristas"].astype(np.int64)
    radios = ven["radios"].astype(np.float64)
    P = nodos[arist[:, 0]]; Q = nodos[arist[:, 1]]
    mid = 0.5 * (P + Q)
    seg_dir = Q - P
    seg_len = np.linalg.norm(seg_dir, axis=1)
    u = seg_dir / np.maximum(seg_len[:, None], 1e-12)

    # --- cuadrante critico + calibre + funde a 100um ---
    cl = clasificar(nodos, arist, radios); cr = critico(cl)
    terminal_r = float(radios.min()); thr_cal = CALIBRE_FACTOR * terminal_r
    interno = ~cl["seg_terminal"]
    seg_cal = interno & (radios > thr_cal)
    i, j, p, cont = cr["i"], cr["j"], cr["p"], cr["cont"]
    sel = (p >= V_MM) & (seg_cal[i] | seg_cal[j])
    i, j, p, cont = i[sel], j[sel], p[sel], cont[sel]
    N = len(i)

    # --- piramides / lobulos ---
    d0 = np.load(find_npz(IN_DOM), allow_pickle=True)
    apex = d0["pyramid_apex"].astype(np.float64)
    axis = d0["pyramid_axis"].astype(np.float64)
    L = d0["pyramid_length"].astype(np.float64)
    axu = axis / np.linalg.norm(axis, axis=1, keepdims=True)
    base = apex + axu * L[:, None]
    NP = len(apex)

    zi = lobulo_de(mid[i], apex, base)
    zj = lobulo_de(mid[j], apex, base)
    intra = zi == zj

    print("=" * 86)
    print("  CAPA 3b - SHUNTS VV INTRA-LOBULARES: DENSIDAD GENUINA vs RUTEO INEFICIENTE")
    print("=" * 86)
    print(f"  capa3b_venoso.npz CANONICO (solo lectura)")
    print(f"  Subconjunto: shunts de calibre que funden a {V_MM*1000:.0f}um = {N}   "
          f"(intra-zona {int(intra.sum())} = {100*intra.mean():.1f}%)")
    print(f"  calibre: interno Y radio>{thr_cal*1000:.1f}um   R_local={R_LOCAL}mm   "
          f"paralelo<{PARALELO_DEG:.0f}deg")
    print("-" * 86)

    # ======================================================================
    #  (1) CORRELACION CARGA (vasos de calibre) vs SHUNTS por lobulo
    # ======================================================================
    lob_seg = lobulo_de(mid, apex, base)             # lobulo de TODO segmento
    carga = np.array([int(np.count_nonzero(seg_cal & (lob_seg == k))) for k in range(NP)])
    sel_intra = intra
    shunts_lob = np.array([int(np.count_nonzero(sel_intra & (zi == k))) for k in range(NP)])

    print("  (1) CORRELACION CARGA-vs-SHUNTS POR LOBULO")
    print("  " + "-" * 82)
    print(f"  {'lobulo':>7} | {'carga (vasos calibre)':>22} | {'shunts intra':>13} | {'shunts/vaso':>12}")
    print("  " + "-" * 82)
    ratios = []
    for k in range(NP):
        r = shunts_lob[k] / carga[k] if carga[k] else 0.0
        ratios.append(r)
        print(f"  {('P%02d'%k):>7} | {carga[k]:>22d} | {shunts_lob[k]:>13d} | {r:>12.4f}")
    ratios = np.array(ratios)
    print("  " + "-" * 82)
    # Pearson r
    if carga.std() > 0 and shunts_lob.std() > 0:
        pear = float(np.corrcoef(carga, shunts_lob)[0, 1])
    else:
        pear = float('nan')
    cv = float(ratios.std() / ratios.mean()) if ratios.mean() > 0 else float('nan')
    print(f"  Pearson r(carga, shunts) = {pear:.3f}   |   "
          f"ratio shunts/vaso: media {ratios.mean():.4f}  CV {cv:.2f}  "
          f"(min {ratios.min():.4f} max {ratios.max():.4f})")
    print("-" * 86)

    # ======================================================================
    #  (2) ESPACIO LIBRE local: ocupacion del cruce vs fondo, y d_next
    # ======================================================================
    tree_mid = cKDTree(mid)
    # ocupacion local en los cruces (excluye los 2 segmentos del propio par)
    listas = tree_mid.query_ball_point(cont, R_LOCAL)
    occ_cruce = np.array([len(s) for s in listas]) - 2     # quita el par
    occ_cruce = np.clip(occ_cruce, 0, None)
    # ocupacion de fondo: misma esfera alrededor de midpoints aleatorios
    idx_f = rng.choice(len(mid), size=min(N_FONDO, len(mid)), replace=False)
    listas_f = tree_mid.query_ball_point(mid[idx_f], R_LOCAL)
    occ_fondo = np.array([len(s) for s in listas_f]) - 1    # quita el propio
    occ_fondo = np.clip(occ_fondo, 0, None)
    # d_next: vaso mas cercano al cruce que NO es del par
    dd, ii = tree_mid.query(cont, k=6)
    d_next = np.full(N, np.inf)
    for a in range(N):
        for col in range(ii.shape[1]):
            s = int(ii[a, col])
            if s != int(i[a]) and s != int(j[a]):
                d_next[a] = dd[a, col]; break
    hay_hueco = d_next > (GAP_FACTOR * p)
    finite = np.isfinite(d_next)

    print("  (2) ESPACIO LIBRE LOCAL (ocupacion en esfera R=%.1fmm)" % R_LOCAL)
    print("  " + "-" * 82)
    print(f"  ocupacion (nº vasos vecinos)   {'p25':>6} {'mediana':>8} {'p75':>6} {'p90':>6} {'max':>6}")
    print(f"     en los CRUCES de calibre :  {pctl(occ_cruce,25):>6.0f} {pctl(occ_cruce,50):>8.0f} "
          f"{pctl(occ_cruce,75):>6.0f} {pctl(occ_cruce,90):>6.0f} {occ_cruce.max():>6.0f}")
    print(f"     FONDO (puntos aleatorios):  {pctl(occ_fondo,25):>6.0f} {pctl(occ_fondo,50):>8.0f} "
          f"{pctl(occ_fondo,75):>6.0f} {pctl(occ_fondo,90):>6.0f} {occ_fondo.max():>6.0f}")
    sobredens = (np.median(occ_cruce) / max(np.median(occ_fondo), 1e-9))
    print(f"     -> sobredensidad cruce/fondo (mediana): {sobredens:.2f}x")
    print(f"  d_next (vaso mas cercano fuera del par, mm): "
          f"p25 {pctl(d_next[finite],25):.3f}  mediana {pctl(d_next[finite],50):.3f}  "
          f"p75 {pctl(d_next[finite],75):.3f}")
    print(f"  shunts con HUECO para separarse (d_next > {GAP_FACTOR:g}*p): "
          f"{int(hay_hueco.sum())} / {N}  ({100*hay_hueco.mean():.1f}%)")
    print("-" * 86)

    # ======================================================================
    #  (3) PARALELISMO del par (solo este subconjunto)
    # ======================================================================
    cosang = np.abs(np.einsum('ij,ij->i', u[i], u[j]))
    cosang = np.clip(cosang, 0.0, 1.0)
    ang = np.degrees(np.arccos(cosang))               # 0=paralelo, 90=perp
    print("  (3) ANGULO DEL PAR (solo shunts de calibre que funden a 100um)")
    print("  " + "-" * 82)
    print(f"  angulo (deg): media {ang.mean():.1f}  mediana {np.median(ang):.1f}  "
          f"p25 {pctl(ang,25):.1f}  p75 {pctl(ang,75):.1f}")
    bins = [0, 10, 20, 30, 45, 60, 90.0001]
    print("  histograma:")
    for b0, b1 in zip(bins[:-1], bins[1:]):
        m = (ang >= b0) & (ang < b1)
        n = int(m.sum())
        barra = "#" * int(round(40 * n / max(N, 1)))
        print(f"     [{b0:>2.0f},{b1:>2.0f})deg : {n:>5d} ({100*n/N:>4.1f}%) {barra}")
    n_par = int((ang < PARALELO_DEG).sum())
    print(f"  casi-paralelos (<{PARALELO_DEG:.0f}deg): {n_par} / {N}  ({100*n_par/N:.1f}%)")
    print("-" * 86)

    # ======================================================================
    #  VEREDICTO (combina los 3)
    # ======================================================================
    s_corr_density = (not np.isnan(pear) and pear > 0.8 and cv < 0.5)
    s_space_routing = (sobredens < 1.5) or (hay_hueco.mean() > 0.5)
    s_parallel_routing = (n_par / N) > 0.5
    votos_ruteo = int(not s_corr_density) + int(s_space_routing) + int(s_parallel_routing)

    print("  VEREDICTO")
    print("  " + "-" * 82)
    print(f"    (1) correlacion: Pearson r={pear:.3f}, CV ratio={cv:.2f} -> "
          f"{'DENSIDAD (carga explica los shunts)' if s_corr_density else 'RUTEO (dispersion alta / no lineal)'}")
    print(f"    (2) espacio:     sobredensidad {sobredens:.2f}x, con hueco {100*hay_hueco.mean():.0f}% -> "
          f"{'RUTEO (hay sitio libre)' if s_space_routing else 'DENSIDAD (saturado)'}")
    print(f"    (3) paralelismo: {100*n_par/N:.0f}% casi-paralelos -> "
          f"{'RUTEO (vasos que corren juntos)' if s_parallel_routing else 'DENSIDAD/inevitable (angulo alto)'}")
    print("  " + "-" * 82)
    if votos_ruteo >= 2:
        print("  >> El empaquetamiento intra-lobular es principalmente RUTEO INEFICIENTE "
              f"({votos_ruteo}/3 señales): defecto corregible, no limite de escala.")
    else:
        print("  >> El empaquetamiento intra-lobular es principalmente DENSIDAD GENUINA "
              f"({3-votos_ruteo}/3 señales): limite de escala, no corregible por ruteo.")
    print("=" * 86)


if __name__ == "__main__":
    main()
