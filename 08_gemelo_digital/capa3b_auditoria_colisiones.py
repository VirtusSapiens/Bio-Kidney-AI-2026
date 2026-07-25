#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa3b_auditoria_colisiones.py  --  Bio-Kidney AI 2026
======================================================
Auditoria de COLISIONES entre el arbol VENOSO (Capa 3b) y el ARTERIAL (Capa 3a).

Una colision = un segmento arterial y uno venoso que se acercan a menos de la SUMA
de sus radios (se solaparian fisicamente; serian imposibles de imprimir / de
coexistir como dos vasos reales). Este script NO modifica nada: solo evalua y
reporta para validar el sesgo satelital de la Capa 3b.

METODO (distancia segmento-segmento en 3D):
  1. Cada arista es un segmento 3D (par de nodos) con su radio (Murray).
  2. CANDIDATOS por cKDTree: se comparan solo los pares arterial-venoso cuyos
     puntos medios distan < UMBRAL_CANDIDATO (evita el O(N*M) de todos contra
     todos). Como los segmentos miden ~PASO_CRECIMIENTO (~0.6 mm), un umbral de
     5 mm entre puntos medios cubre con holgura cualquier par que pueda solapar.
  3. Para cada candidato, distancia MINIMA real segmento-segmento (no solo entre
     puntos medios): algoritmo de Ericson (Real-Time Collision Detection).
  4. COLISION si  dist_min < (radio_art + radio_ven) * FACTOR_SEGURIDAD.

Ademas, para caracterizar "cuan cerca corren los arboles en general", se calcula
para CADA segmento venoso su distancia al segmento arterial mas cercano (via los
K vecinos mas proximos por punto medio, exacta), sin truncar por el umbral.

Solo NumPy / SciPy. SIN bpy. Entorno: env_biokidney.

Entrada : capa3a_arterial.npz, capa3b_venoso.npz
Salida  : capa3b_auditoria_colisiones.npz  +  capa3b_auditoria_colisiones.txt
"""

import os
import itertools
import numpy as np
from scipy.spatial import cKDTree

# ============================================================================
#  PARAMETROS  (editar aqui)
# ============================================================================
FACTOR_SEGURIDAD  = 1.0    # colision si dist < (r_a + r_v) * FACTOR. Subirlo exige
                           # MAS separacion (margen de pared/impresion).
UMBRAL_CANDIDATO  = 5.0    # mm  pares con puntos medios mas cerca que esto -> candidatos
K_VECINOS         = 8      # vecinos mas cercanos por punto medio para la distancia
                           # exacta al segmento arterial mas proximo (distribucion global)
RADIO_HILIO_MM    = 10.0   # mm  una colision se considera "cerca del hilio" si su punto
                           # de contacto esta a < esto del hilio (ambos troncos convergen ahi)
CHUNK             = 2_000_000   # pares por bloque (limita memoria en seg-seg)

IN_ART = "capa3a_arterial.npz"
IN_VEN = "capa3b_venoso.npz"
OUT_NPZ = "capa3b_auditoria_colisiones.npz"
OUT_TXT = "capa3b_auditoria_colisiones.txt"


# ============================================================================
#  UTILIDADES
# ============================================================================
def find_npz(name):
    here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() \
        else os.getcwd()
    for c in (
        os.path.expanduser(os.path.join("~", "Escritorio", "BioKidney-AI", name)),
        os.path.join(here, name),
        os.path.join(here, "..", name),
        os.path.join(os.getcwd(), name),
    ):
        if os.path.isfile(c):
            return os.path.abspath(c)
    raise FileNotFoundError(f"No se encontro {name}")


def seg_seg_dist(P1, Q1, P2, Q2):
    """Distancia minima entre segmentos P1-Q1 y P2-Q2 (vectorizado, N pares).

    Algoritmo de closest-point segment-segment de Ericson (Real-Time Collision
    Detection, 2005), generalizado a lotes con NumPy y guardas para segmentos
    degenerados (longitud ~0).

    Devuelve (dist (N,), contacto (N,3))  donde 'contacto' es el punto medio del
    par de puntos mas cercanos (util para localizar la colision en el espacio).
    """
    d1 = Q1 - P1
    d2 = Q2 - P2
    r = P1 - P2
    a = np.einsum("ij,ij->i", d1, d1)     # |d1|^2
    e = np.einsum("ij,ij->i", d2, d2)     # |d2|^2
    f = np.einsum("ij,ij->i", d2, r)
    c = np.einsum("ij,ij->i", d1, r)
    b = np.einsum("ij,ij->i", d1, d2)
    EPS = 1e-12
    a_safe = np.where(a > EPS, a, 1.0)
    e_safe = np.where(e > EPS, e, 1.0)

    denom = a * e - b * b
    # caso general (segmentos no paralelos): s optimo en la recta, acotado a [0,1]
    s = np.where(denom > EPS, np.clip((b * f - c * e) / np.where(denom > EPS, denom, 1.0),
                                      0.0, 1.0), 0.0)
    t = (b * s + f) / e_safe
    # re-acotar t a [0,1] y recomputar s en consecuencia
    t_lt = t < 0.0
    t_gt = t > 1.0
    t = np.clip(t, 0.0, 1.0)
    s = np.where(t_lt, np.clip(-c / a_safe, 0.0, 1.0), s)
    s = np.where(t_gt, np.clip((b - c) / a_safe, 0.0, 1.0), s)
    # guardas de degeneracion
    a_deg = a <= EPS                       # segmento arterial = punto
    s = np.where(a_deg, 0.0, s)
    t = np.where(a_deg, np.clip(f / e_safe, 0.0, 1.0), t)
    e_deg = e <= EPS                       # segmento venoso = punto
    t = np.where(e_deg, 0.0, t)
    s = np.where(e_deg, np.clip(-c / a_safe, 0.0, 1.0), s)

    c1 = P1 + d1 * s[:, None]
    c2 = P2 + d2 * t[:, None]
    diff = c1 - c2
    dist = np.sqrt(np.einsum("ij,ij->i", diff, diff))
    contacto = 0.5 * (c1 + c2)
    return dist, contacto


def deciles_str(vals):
    """Cadena con los deciles (0,10,...,100 %) de 'vals'."""
    qs = np.percentile(vals, np.arange(0, 101, 10))
    return "  ".join(f"{q:.2f}" for q in qs)


# ============================================================================
#  MAIN
# ============================================================================
def main():
    art_path = find_npz(IN_ART)
    ven_path = find_npz(IN_VEN)

    art = np.load(art_path, allow_pickle=False)
    ven = np.load(ven_path, allow_pickle=False)

    nodos_a = art["nodos"].astype(np.float64)
    arist_a = art["aristas"].astype(np.int64)
    radios_a = art["radios"].astype(np.float64)           # por arista
    nodos_v = ven["nodos"].astype(np.float64)
    arist_v = ven["aristas"].astype(np.int64)
    radios_v = ven["radios"].astype(np.float64)
    hilio = ven["raiz_pos"].astype(np.float64) if "raiz_pos" in ven.files \
        else art["raiz"].astype(np.float64)

    # endpoints y puntos medios de cada segmento
    P1 = nodos_a[arist_a[:, 0]]; Q1 = nodos_a[arist_a[:, 1]]
    P2 = nodos_v[arist_v[:, 0]]; Q2 = nodos_v[arist_v[:, 1]]
    mid_a = 0.5 * (P1 + Q1)
    mid_v = 0.5 * (P2 + Q2)
    na, nv = len(arist_a), len(arist_v)

    R = []
    def out(line=""):
        print(line)
        R.append(line)

    out("=" * 70)
    out("  CAPA 3b - AUDITORIA DE COLISIONES  (venoso vs arterial)")
    out("=" * 70)
    out(f"  Arterial : {art_path}")
    out(f"  Venoso   : {ven_path}")
    out("-" * 70)
    out(f"  Segmentos arteriales : {na}")
    out(f"  Segmentos venosos    : {nv}")
    out(f"  FACTOR_SEGURIDAD = {FACTOR_SEGURIDAD}   UMBRAL_CANDIDATO = {UMBRAL_CANDIDATO} mm")
    out(f"  Radio arterial [{radios_a.min():.4f}, {radios_a.max():.4f}] mm   "
        f"venoso [{radios_v.min():.4f}, {radios_v.max():.4f}] mm")
    out("-" * 70)

    tree_a = cKDTree(mid_a)

    # ------------------------------------------------------------------
    #  A) DISTRIBUCION GLOBAL: distancia de cada segmento venoso al
    #     segmento arterial mas cercano (exacta, via K vecinos por mid).
    # ------------------------------------------------------------------
    k = min(K_VECINOS, na)
    _, idx_knn = tree_a.query(mid_v, k=k)
    if idx_knn.ndim == 1:
        idx_knn = idx_knn[:, None]
    dmin_v = np.full(nv, np.inf)
    contacto_v = np.zeros((nv, 3))
    for j in range(idx_knn.shape[1]):
        aj = idx_knn[:, j]
        dj, cj = seg_seg_dist(P1[aj], Q1[aj], P2, Q2)
        mask = dj < dmin_v
        dmin_v[mask] = dj[mask]
        contacto_v[mask] = cj[mask]

    dist_min_global = float(dmin_v.min())
    dist_media_segseg = float(dmin_v.mean())

    # ------------------------------------------------------------------
    #  B) CANDIDATOS y COLISIONES: todos los pares con mid a < UMBRAL.
    # ------------------------------------------------------------------
    neigh = tree_a.query_ball_point(mid_v, r=UMBRAL_CANDIDATO, workers=-1)
    counts = np.fromiter((len(x) for x in neigh), dtype=np.int64, count=nv)
    n_cand = int(counts.sum())
    vidx = np.repeat(np.arange(nv, dtype=np.int64), counts)
    aidx = np.fromiter(itertools.chain.from_iterable(neigh), dtype=np.int64,
                       count=n_cand)

    col_a, col_v, col_d, col_c = [], [], [], []
    for s0 in range(0, n_cand, CHUNK):
        s1 = min(s0 + CHUNK, n_cand)
        ai = aidx[s0:s1]; vi = vidx[s0:s1]
        dd, cc = seg_seg_dist(P1[ai], Q1[ai], P2[vi], Q2[vi])
        umbral = (radios_a[ai] + radios_v[vi]) * FACTOR_SEGURIDAD
        hit = dd < umbral
        if np.any(hit):
            col_a.append(ai[hit]); col_v.append(vi[hit])
            col_d.append(dd[hit]); col_c.append(cc[hit])

    if col_a:
        col_a = np.concatenate(col_a); col_v = np.concatenate(col_v)
        col_d = np.concatenate(col_d); col_c = np.concatenate(col_c)
    else:
        col_a = np.empty(0, np.int64); col_v = np.empty(0, np.int64)
        col_d = np.empty(0); col_c = np.empty((0, 3))

    n_col = len(col_a)
    n_ven_col = int(np.unique(col_v).size) if n_col else 0
    pct_col = 100.0 * n_col / nv
    pct_ven_col = 100.0 * n_ven_col / nv

    # desglose: ¿colisiones cerca del hilio (troncos convergentes) o repartidas?
    if n_col:
        d_hil = np.linalg.norm(col_c - hilio, axis=1)
        n_cerca = int(np.count_nonzero(d_hil < RADIO_HILIO_MM))
        pct_cerca = 100.0 * n_cerca / n_col
    else:
        d_hil = np.empty(0)
        n_cerca = 0; pct_cerca = 0.0

    # ========================================================================
    #  REPORTE
    # ========================================================================
    out(f"  Pares de segmentos CANDIDATOS evaluados : {n_cand}")
    out(f"     (mid arterial-venoso < {UMBRAL_CANDIDATO} mm; evita {na}x{nv} = {na*nv} comparaciones)")
    out("")
    out(f"  COLISIONES detectadas : {n_col}   "
        f"(dist < (r_art + r_ven) x {FACTOR_SEGURIDAD})")
    out(f"     % sobre segmentos venosos (pares)      : {pct_col:.3f} %")
    out(f"     segmentos venosos distintos en colision: {n_ven_col} / {nv}  ({pct_ven_col:.3f} %)")
    out("")
    out(f"  Distancia MINIMA arterial-venosa global (peor caso): {dist_min_global:.4f} mm")
    out(f"  Distancia MEDIA segmento-segmento (al arterial mas cercano): {dist_media_segseg:.3f} mm")
    out(f"     -> referencia: offset medio nodo-arterial reportado en 3b ~2.9 mm  "
        f"({'consistente' if 2.0 <= dist_media_segseg <= 3.8 else 'revisar'})")
    out("")
    out(f"  HISTOGRAMA de distancias minimas venoso->arterial (deciles 0..100%, mm):")
    out(f"     {deciles_str(dmin_v)}")
    out("")
    out(f"  DESGLOSE de colisiones (cerca del hilio vs repartidas):")
    if n_col:
        out(f"     contacto a < {RADIO_HILIO_MM} mm del hilio : {n_cerca} / {n_col}  ({pct_cerca:.1f} %)")
        out(f"     dist. de la colision al hilio: min {d_hil.min():.2f}  "
            f"mediana {np.median(d_hil):.2f}  max {d_hil.max():.2f} mm")
        if pct_cerca >= 80.0:
            out(f"     -> CONCENTRADAS cerca del hilio (artefacto: ambos troncos "
                f"comparten el mismo puerto/raiz)")
        elif pct_cerca <= 20.0:
            out(f"     -> REPARTIDAS por el parenquima (el sesgo satelital no separa lo suficiente)")
        else:
            out(f"     -> MIXTAS (parte en el hilio, parte repartidas)")
    else:
        out(f"     (sin colisiones: nada que desglosar)")
    out("-" * 70)

    # --- guardar npz --------------------------------------------------------
    out_npz = os.path.join(os.path.dirname(ven_path), OUT_NPZ)
    np.savez_compressed(
        out_npz,
        col_arista_arterial=col_a.astype(np.int32),
        col_arista_venosa=col_v.astype(np.int32),
        col_distancia=col_d.astype(np.float32),
        col_contacto=col_c.astype(np.float32),
        col_dist_hilio=d_hil.astype(np.float32),
        dist_min_venoso=dmin_v.astype(np.float32),       # por segmento venoso
        # metadata
        factor_seguridad=np.float64(FACTOR_SEGURIDAD),
        umbral_candidato=np.float64(UMBRAL_CANDIDATO),
        radio_hilio_mm=np.float64(RADIO_HILIO_MM),
        hilio=hilio.astype(np.float64),
        n_seg_arterial=np.int32(na),
        n_seg_venoso=np.int32(nv),
        n_candidatos=np.int32(n_cand),
        n_colisiones=np.int32(n_col),
        n_venosos_en_colision=np.int32(n_ven_col),
        dist_min_global=np.float64(dist_min_global),
        dist_media_segseg=np.float64(dist_media_segseg),
    )

    # --- VEREDICTO (una linea) ----------------------------------------------
    veredicto = (f"VEREDICTO: {n_col} colisiones de {nv} segmentos venosos "
                 f"({pct_col:.3f}%); separacion minima arterial-venosa = "
                 f"{dist_min_global:.4f} mm.")
    out("")
    out(veredicto)
    out("=" * 70)
    out(f"  -> npz : {out_npz}")
    out_txt = os.path.join(os.path.dirname(ven_path), OUT_TXT)
    out(f"  -> txt : {out_txt}")

    with open(out_txt, "w", encoding="utf-8") as fh:
        fh.write("\n".join(R) + "\n")


if __name__ == "__main__":
    main()
