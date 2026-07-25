#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa3b_shunts_cola_corregible.py  --  Bio-Kidney AI 2026
========================================================
AUDITORIA DE SOLO LECTURA. NO escribe ningun .npz, NO toca el generador ni el
visualizador. Solo imprime a stdout. Trabaja sobre el canonico capa3b_venoso.npz.

CONTEXTO: el discriminador previo (capa3b_shunts_densidad_vs_ruteo.py) hallo que los
1927 shunts VV de calibre (funden a 100um) son PISO DE DENSIDAD GENUINA + COLA DE
RUTEO CORREGIBLE. PERO su test de paralelismo (54.5% casi-paralelos <20deg) esta
SESGADO: todos los segmentos miden 0.6mm (paso_crecimiento fijo), asi que dos vasos
que GLOBALMENTE divergen tienen tramos locales casi-paralelos por cuantizacion.

TAREA: medir el tamaño REAL de la cola de ruteo corregible, descontando el artefacto
de malla de 0.6mm.

OJO -- HALLAZGO DE ESCALA: definir "vaso" como cadena de segmentos entre nodos de
union (grado != 2) NO funciona en este arbol: es un CCO que bifurca en casi cada
nodo, asi que el 94% de esas cadenas es UN SOLO segmento de 0.6mm -> la cadena inter-
bifurcacion ES el segmento de malla y el angulo "vaso-a-vaso" coincide con el seg-a-
seg (no deflacta nada; es tautologia). Se reporta esa version como BASELINE con su
diagnostico, y se mide la direccion global de verdad por la LINEA PRINCIPAL del vaso:
continuacion por el hijo de mayor radio (Murray) en una ventana de +-L a lo largo del
arbol (atraviesa las bifurcaciones). Ese es el "vaso anatomico" real.

  1) PARALELISMO a tres escalas: seg-a-seg (0.6mm), cadena-inter-bifurcacion
     (degenerada), y LINEA-PRINCIPAL +-L (real). ¿Cuanto del casi-paralelo sobrevive
     al mirar a escala de vaso de verdad? Barrido en L.
  2) COLA CORREGIBLE DEPURADA: par casi-paralelo REAL = angulo linea-principal < 20deg.
     nº absoluto, % sobre 1927, % sobre el total de shunts venosos.
  3) SEPARABILIDAD: perpendicular a la direccion global del par, ¿hay hueco lateral
     (> 2x radio) en al menos un lado para separarlos sin chocar con un TERCER vaso?
     Distancia REAL punto-a-segmento (no el d_next de midpoints, contaminado por el
     suelo de 0.3mm). Se excluyen los DOS vasos del par. Se reporta que fraccion de la
     holgura viene de un sector lateral VACIO (saturacion del radio de busqueda).
  4) DESCOMPOSICION de los 1927 en 3 cubetas:
       (a) densidad irreducible : no casi-paralelos (linea-principal), o encajonados
       (b) ruteo corregible separable : casi-paralelos reales con hueco lateral > 2r
       (c) ambiguos : casi-paralelos con hueco lateral marginal (entre 1r y 2r)

Solo NumPy / SciPy. SIN bpy. Entorno: env_biokidney.
"""

import os
import sys
from collections import defaultdict
import numpy as np
from scipy.spatial import cKDTree

# ============================================================================
#  PARAMETROS
# ============================================================================
CALIBRE_FACTOR = 2.0       # "de calibre": radio > 2x terminal
V_MM           = 0.100     # 100 um: umbral de "funde lumenes"
PARALELO_DEG   = 20.0      # umbral de "casi paralelo"
L_MAIN         = 1.5       # mm  semiventana de la linea principal (ventana ~2L ~ 3mm)
L_BARRIDO      = [0.6, 1.0, 1.5, 2.0, 3.0]   # sensibilidad al tamaño de ventana

R_SEP   = 1.0              # mm  radio de busqueda de terceros vasos alrededor del cruce
AX_WIN  = 0.5              # mm  semiventana axial (misma seccion transversal del cruce)
K_SECT  = 12              # nº de sectores angulares laterales (30 deg c/u)

IN_VEN = "capa3b_venoso.npz"    # CANONICO (solo lectura)
IN_DOM = "capa0_dominio.npz"

_here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
if _here not in sys.path:
    sys.path.insert(0, _here)
from capa3b_clasificar_vv import clasificar, find_npz
from capa3b_severidad_vv import critico


def pctl(x, q):
    return np.percentile(x, q) if len(x) else float('nan')


# ============================================================================
#  BASELINE: vaso = cadena de segmentos entre nodos de union (grado != 2)
# ============================================================================
def cadenas_interbifurcacion(nodos, arist):
    """seg_vessel[e] = id de la cadena; dirs[v] = direccion (PCA/cuerda); cnt = nº seg."""
    n_nodos = len(nodos); n_seg = len(arist)
    grado = np.bincount(arist.ravel(), minlength=n_nodos)
    adj = defaultdict(list)
    for e in range(n_seg):
        a, b = int(arist[e, 0]), int(arist[e, 1])
        adj[a].append((b, e)); adj[b].append((a, e))
    union = grado != 2
    seg_vessel = -np.ones(n_seg, dtype=np.int64)
    dirs = []; cnt = []
    vid = 0
    for jn in np.where(union)[0]:
        for (nb, e0) in adj[jn]:
            if seg_vessel[e0] != -1:
                continue
            seg_vessel[e0] = vid; node_seq = [jn, nb]
            arrival = e0; cur = nb
            while not union[cur]:
                nxt = next((nn, ee) for (nn, ee) in adj[cur] if ee != arrival)
                nn, ee = nxt
                seg_vessel[ee] = vid; node_seq.append(nn)
                arrival = ee; cur = nn
            V = nodos[np.array(node_seq)]; chord = V[-1] - V[0]
            if len(V) >= 3:
                _, _, vt = np.linalg.svd(V - V.mean(0), full_matrices=False)
                d = vt[0]
                if d @ chord < 0:
                    d = -d
            else:
                d = chord
            nn_ = np.linalg.norm(d)
            dirs.append(d / nn_ if nn_ > 1e-12 else np.array([1.0, 0.0, 0.0]))
            cnt.append(len(node_seq) - 1)
            vid += 1
    assert not np.any(seg_vessel < 0)
    return grado, seg_vessel, np.array(dirs), np.array(cnt)


# ============================================================================
#  REAL: linea principal del vaso (continuacion por el hijo de mayor radio)
# ============================================================================
def arbol_padres_e_hijo_principal(nodos, arist, radios):
    """parent[nodo] y main_child[nodo] (hijo por el que continua la linea = mayor radio)."""
    n_nodos = len(nodos)
    parent = -np.ones(n_nodos, dtype=np.int64)
    parent[arist[:, 1]] = arist[:, 0]                 # arista = (padre, hijo)
    main_child = -np.ones(n_nodos, dtype=np.int64)
    best_r = np.full(n_nodos, -1.0)
    for e in range(len(arist)):
        par = int(arist[e, 0])
        if radios[e] > best_r[par]:
            best_r[par] = radios[e]; main_child[par] = int(arist[e, 1])
    return parent, main_child


def anclas_linea_principal(edges, arist, nodos, parent, main_child, L):
    """Para cada segmento devuelve (p_up, p_dn): el ancla a distancia ~L hacia la raiz
       (por el padre) y hacia la hoja (por la linea principal = hijo de mayor radio),
       atravesando bifurcaciones. Tambien los nodos extremos propios del segmento."""
    n = len(edges)
    p_up = np.zeros((n, 3)); p_dn = np.zeros((n, 3))
    n_par = nodos[arist[edges, 0]]; n_hij = nodos[arist[edges, 1]]
    for k, e in enumerate(edges):
        cur = int(arist[e, 0]); d = 0.0
        while d < L and parent[cur] >= 0:
            nxt = int(parent[cur]); d += float(np.linalg.norm(nodos[nxt] - nodos[cur])); cur = nxt
        p_up[k] = nodos[cur]
        cur = int(arist[e, 1]); d = 0.0
        while d < L and main_child[cur] >= 0:
            nxt = int(main_child[cur]); d += float(np.linalg.norm(nodos[nxt] - nodos[cur])); cur = nxt
        p_dn[k] = nodos[cur]
    return p_up, p_dn, n_par, n_hij


def _norm(v):
    return v / np.maximum(np.linalg.norm(v, axis=1, keepdims=True), 1e-12)


def heading_linea_principal(edges, arist, nodos, parent, main_child, L):
    """Direccion global (cuerda p_up->p_dn) unitaria de cada segmento."""
    p_up, p_dn, _, _ = anclas_linea_principal(edges, arist, nodos, parent, main_child, L)
    return _norm(p_dn - p_up)


def punto_a_segmentos(c, A, B):
    AB = B - A
    den = np.maximum(np.einsum('ij,ij->i', AB, AB), 1e-12)
    t = np.clip(np.einsum('ij,ij->i', c - A, AB) / den, 0.0, 1.0)
    foot = A + t[:, None] * AB
    return foot, np.linalg.norm(foot - c, axis=1)


def hueco_lateral(c, a_axis, cand_A, cand_B):
    """Mayor holgura lateral libre (mm) perpendicular al eje a_axis alrededor de c.
       Devuelve (best_clearance, sector_vacio_bool)."""
    tmp = np.array([1.0, 0.0, 0.0]) if abs(a_axis[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = tmp - (tmp @ a_axis) * a_axis; e1 /= np.linalg.norm(e1)
    e2 = np.cross(a_axis, e1)
    if len(cand_A) == 0:
        return R_SEP, True
    foot, _ = punto_a_segmentos(c, cand_A, cand_B)
    v = foot - c
    axial = v @ a_axis
    vperp = v - np.outer(axial, a_axis)
    rho = np.linalg.norm(vperp, axis=1)
    m = (np.abs(axial) < AX_WIN) & (rho < R_SEP) & (rho > 1e-9)
    if not np.any(m):
        return R_SEP, True
    x = vperp[m] @ e1; y = vperp[m] @ e2
    sect = ((np.arctan2(y, x) + np.pi) / (2 * np.pi) * K_SECT).astype(int) % K_SECT
    rr = rho[m]
    clear = np.full(K_SECT, R_SEP)
    for s in range(K_SECT):
        seln = sect == s
        if np.any(seln):
            clear[s] = rr[seln].min()
    best = float(clear.max())
    return best, bool(np.any(clear >= R_SEP - 1e-9))


def angulo(dA, dB):
    return np.degrees(np.arccos(np.clip(np.abs(np.einsum('ij,ij->i', dA, dB)), 0, 1)))


def main():
    ven = np.load(find_npz(IN_VEN), allow_pickle=False)
    nodos = ven["nodos"].astype(np.float64)
    arist = ven["aristas"].astype(np.int64)
    radios = ven["radios"].astype(np.float64)
    P = nodos[arist[:, 0]]; Q = nodos[arist[:, 1]]
    mid = 0.5 * (P + Q)
    u = (Q - P) / np.maximum(np.linalg.norm(Q - P, axis=1)[:, None], 1e-12)

    # --- subconjunto CANONICO: shunts de calibre que funden a 100um ---
    cl = clasificar(nodos, arist, radios); cr = critico(cl)
    terminal_r = float(radios.min()); thr_cal = CALIBRE_FACTOR * terminal_r
    seg_cal = (~cl["seg_terminal"]) & (radios > thr_cal)
    i0, j0, p0, cont0 = cr["i"], cr["j"], cr["p"], cr["cont"]
    n_shunt_total = int(np.count_nonzero(p0 >= V_MM))
    sel = (p0 >= V_MM) & (seg_cal[i0] | seg_cal[j0])
    i, j, p, cont = i0[sel], j0[sel], p0[sel], cont0[sel]
    N = len(i)
    r_pair = np.maximum(radios[i], radios[j])

    # --- lobulos (reportar intra) ---
    d0 = np.load(find_npz(IN_DOM), allow_pickle=True)
    apex = d0["pyramid_apex"].astype(np.float64)
    axu = d0["pyramid_axis"].astype(np.float64)
    axu = axu / np.linalg.norm(axu, axis=1, keepdims=True)
    base = apex + axu * d0["pyramid_length"].astype(np.float64)[:, None]

    def lobulo_de(Pm):
        dist = np.empty((len(Pm), len(apex)))
        for k in range(len(apex)):
            AB = base[k] - apex[k]; den = float(AB @ AB)
            t = np.clip(((Pm - apex[k]) @ AB) / den, 0.0, 1.0)
            dist[:, k] = np.linalg.norm(Pm - (apex[k] + t[:, None] * AB), axis=1)
        return dist.argmin(axis=1)
    intra = lobulo_de(mid[i]) == lobulo_de(mid[j])

    print("=" * 88)
    print("  CAPA 3b - COLA DE RUTEO CORREGIBLE REAL (paralelismo a escala de vaso)")
    print("=" * 88)
    print(f"  capa3b_venoso.npz CANONICO (solo lectura)")
    print(f"  Shunts de calibre que funden a {V_MM*1000:.0f}um = {N}   "
          f"(intra-lobulares {int(intra.sum())} = {100*intra.mean():.1f}%)")
    print(f"  Total shunts venosos a {V_MM*1000:.0f}um = {n_shunt_total}")
    print("-" * 88)

    # ======================================================================
    #  (1) PARALELISMO a tres escalas
    # ======================================================================
    grado, seg_vessel, dirs_ch, cnt_ch = cadenas_interbifurcacion(nodos, arist)
    vi, vj = seg_vessel[i], seg_vessel[j]
    parent, main_child = arbol_padres_e_hijo_principal(nodos, arist, radios)

    ang_seg = angulo(u[i], u[j])                                    # seg-a-seg (sesgado)
    ang_ch = angulo(dirs_ch[vi], dirs_ch[vj]); ang_ch[vi == vj] = 0.0   # cadena (degenerada)

    print("  DIAGNOSTICO DE ESCALA:")
    print(f"     segmentos por cadena inter-bifurcacion: mediana {int(np.median(cnt_ch))}, "
          f"{100*np.mean(cnt_ch==1):.1f}% son UN solo segmento de 0.6mm")
    print(f"     -> 'cadena inter-bifurcacion' == segmento de malla: la definicion literal "
          f"de vaso es degenerada aqui")
    print("  " + "-" * 84)

    # linea principal a varias L
    edges_i = i; edges_j = j
    print(f"  {'escala':<28}{'mediana':>9}{'<20deg':>16}")
    n_lt = int((ang_seg < PARALELO_DEG).sum())
    print(f"  {'seg-a-seg (0.6mm)':<28}{np.median(ang_seg):>9.1f}{n_lt:>8d} ({100*n_lt/N:>4.1f}%)")
    n_lt = int((ang_ch < PARALELO_DEG).sum())
    print(f"  {'cadena inter-bifurc.':<28}{np.median(ang_ch):>9.1f}{n_lt:>8d} ({100*n_lt/N:>4.1f}%)  <== degenerada")
    ang_main = None
    for L in L_BARRIDO:
        hi = heading_linea_principal(edges_i, arist, nodos, parent, main_child, L)
        hj = heading_linea_principal(edges_j, arist, nodos, parent, main_child, L)
        a = angulo(hi, hj)
        n_lt = int((a < PARALELO_DEG).sum())
        etiq = f"linea-principal +-{L:.1f}mm"
        marca = "  <== REAL (primaria)" if abs(L - L_MAIN) < 1e-9 else ""
        print(f"  {etiq:<28}{np.median(a):>9.1f}{n_lt:>8d} ({100*n_lt/N:>4.1f}%){marca}")
        if abs(L - L_MAIN) < 1e-9:
            ang_main = a
    print("  " + "-" * 84)
    print(f"  histograma linea-principal +-{L_MAIN}mm (direccion global REAL del vaso):")
    for b0, b1 in zip([0, 10, 20, 30, 45, 60], [10, 20, 30, 45, 60, 90.0001]):
        m = (ang_main >= b0) & (ang_main < b1); n = int(m.sum())
        print(f"     [{b0:>2.0f},{b1:>2.0f})deg : {n:>5d} ({100*n/N:>4.1f}%) {'#'*int(round(40*n/max(N,1)))}")
    n_seg_par = int((ang_seg < PARALELO_DEG).sum()); n_main_par = int((ang_main < PARALELO_DEG).sum())
    print(f"  PREMISA INVERTIDA: al quitar el ruido de malla el paralelismo SUBE, no baja;")
    print(f"  el casi-paralelo <20deg crece x{n_main_par/max(n_seg_par,1):.2f} "
          f"(de {n_seg_par} seg-a-seg a {n_main_par} a escala de vaso +-{L_MAIN}mm).")
    # --- confound: paralelismo por convergencia al colector (raiz) vs correr juntos ---
    pu_i, pd_i, np_i, nh_i = anclas_linea_principal(i, arist, nodos, parent, main_child, L_MAIN)
    pu_j, pd_j, np_j, nh_j = anclas_linea_principal(j, arist, nodos, parent, main_child, L_MAIN)
    ang_root = angulo(_norm(np_i - pu_i), _norm(np_j - pu_j))   # tramo HACIA RAIZ (colector)
    ang_leaf = angulo(_norm(pd_i - nh_i), _norm(pd_j - nh_j))   # tramo HACIA HOJA (territorio)
    print(f"  descomposicion del heading (mediana angulo entre los dos vasos):")
    print(f"     tramo HACIA RAIZ  (convergencia al colector): {np.median(ang_root):.1f}deg  "
          f"<20deg {100*(ang_root<PARALELO_DEG).mean():.0f}%")
    print(f"     tramo HACIA HOJA  (territorios de drenaje)  : {np.median(ang_leaf):.1f}deg  "
          f"<20deg {100*(ang_leaf<PARALELO_DEG).mean():.0f}%")
    print("-" * 88)

    # ======================================================================
    #  (2) COLA CORREGIBLE DEPURADA
    #  El paralelismo "full" (80%) es convergencia al colector (no corregible).
    #  CORRIDA REDUNDANTE = paralelo en AMBOS sentidos (hacia raiz Y hacia hoja):
    #  vasos que de verdad corren juntos, no una confluencia en Y con territorios
    #  divergentes. Esa es la unica cola que un offset lateral podria arreglar sin
    #  romper la topologia de drenaje.
    # ======================================================================
    conflu = (ang_root < PARALELO_DEG) & (ang_leaf >= PARALELO_DEG)   # Y de tributarios
    redund = (ang_root < PARALELO_DEG) & (ang_leaf < PARALELO_DEG)    # corrida redundante
    real_par = redund
    n_full = int((ang_main < PARALELO_DEG).sum())
    n_real = int(real_par.sum())
    print("  (2) COLA CORREGIBLE DEPURADA")
    print("  " + "-" * 84)
    print(f"  paralelo 'full' +-{L_MAIN}mm (INFLADO por convergencia al colector): "
          f"{n_full} ({100*n_full/N:.1f}%)  -> NO es la cola corregible")
    print(f"  confluencia en Y (paralelo hacia raiz, DIVERGENTE hacia hoja = tributarios "
          f"de distinto territorio): {int(conflu.sum())} ({100*conflu.sum()/N:.1f}%)  "
          f"-> anatomia, no ruteo")
    print(f"  CORRIDA REDUNDANTE (paralelo en AMBOS sentidos = vasos que corren juntos): "
          f"{n_real}   ({100*n_real/N:.1f}% de {N}, {100*n_real/max(n_shunt_total,1):.1f}% "
          f"de {n_shunt_total} venosos)")
    print("-" * 88)

    # ======================================================================
    #  (3) SEPARABILIDAD de la cola real (distancia REAL punto-a-segmento)
    # ======================================================================
    hi_main = heading_linea_principal(edges_i, arist, nodos, parent, main_child, L_MAIN)
    hj_main = heading_linea_principal(edges_j, arist, nodos, parent, main_child, L_MAIN)
    tree_mid = cKDTree(mid)
    clear = np.full(N, np.nan); vacio = np.zeros(N, bool)
    for a in np.where(real_par)[0]:
        c = cont[a]
        ax = hi_main[a] + (hj_main[a] if (hi_main[a] @ hj_main[a]) >= 0 else -hj_main[a])
        na = np.linalg.norm(ax); ax = ax / na if na > 1e-12 else hi_main[a]
        cand = np.array(tree_mid.query_ball_point(c, R_SEP + 0.4), dtype=np.int64)
        if len(cand):
            cand = cand[(seg_vessel[cand] != vi[a]) & (seg_vessel[cand] != vj[a])]
        clear[a], vacio[a] = hueco_lateral(c, ax, P[cand], Q[cand])

    gap_need = 2.0 * r_pair
    cl_real = clear[real_par]
    sep = real_par & (clear > gap_need)
    box = real_par & (clear <= r_pair)
    amb = real_par & (clear > r_pair) & (clear <= gap_need)

    print("  (3) SEPARABILIDAD de la cola real (holgura lateral punto-a-segmento REAL)")
    print("  " + "-" * 84)
    print(f"  holgura lateral (mm): p25 {pctl(cl_real,25):.3f}  mediana {pctl(cl_real,50):.3f}  "
          f"p75 {pctl(cl_real,75):.3f}   (R_busqueda={R_SEP}mm)")
    print(f"  con sector lateral VACIO (holgura llega al tope {R_SEP}mm, no medida): "
          f"{int(vacio[real_par].sum())}/{n_real} ({100*vacio[real_par].mean():.0f}%)")
    print(f"  umbral 2r: mediana {2*np.median(r_pair[real_par]):.4f} mm")
    print(f"  separable (holgura > 2r)      : {int(sep.sum())}  ({100*sep.sum()/max(n_real,1):.1f}% de la cola)")
    print(f"  encajonado (holgura <= 1r)    : {int(box.sum())}  ({100*box.sum()/max(n_real,1):.1f}% de la cola)")
    print(f"  marginal  (1r < holgura <= 2r): {int(amb.sum())}  ({100*amb.sum()/max(n_real,1):.1f}% de la cola)")
    print("-" * 88)

    # ======================================================================
    #  (4) DESCOMPOSICION FINAL de los 1927
    # ======================================================================
    b_sep = int(sep.sum()); c_amb = int(amb.sum()); a_dens = N - b_sep - c_amb
    print("  (4) DESCOMPOSICION DE LOS %d SHUNTS DE CALIBRE" % N)
    print("  " + "-" * 84)
    print(f"  (a) densidad/topologia irreducible : {a_dens:>5d}  ({100*a_dens/N:>4.1f}%)")
    print(f"        confluencias en Y (tributarios divergentes) {int(conflu.sum())} + "
          f"no-paralelos {int((~(conflu|redund)).sum())} + encajonados {int(box.sum())}")
    print(f"  (b) ruteo corregible SEPARABLE     : {b_sep:>5d}  ({100*b_sep/N:>4.1f}%)")
    print(f"        corrida redundante (paralela en ambos sentidos) con hueco lateral > 2r")
    print(f"  (c) ambiguos (hueco marginal)      : {c_amb:>5d}  ({100*c_amb/N:>4.1f}%)")
    print(f"        corrida redundante con hueco lateral entre 1r y 2r")
    print(f"  suma: {a_dens + b_sep + c_amb} == {N}")
    print("-" * 88)

    # ======================================================================
    #  VEREDICTO
    # ======================================================================
    print("  VEREDICTO")
    print("  " + "-" * 84)
    print(f"  De los {N} shunts de calibre, la cubeta (b) = ruteo-corregible-separable REAL "
          f"es {b_sep}")
    print(f"  ({100*b_sep/N:.1f}% de {N}; {100*b_sep/max(n_shunt_total,1):.1f}% de los "
          f"{n_shunt_total} venosos): ese es el maximo que un offset lateral podria eliminar.")
    print(f"  OJO: el test seg-a-seg previo (54.5% paralelo) NO estaba inflado por la malla "
          f"-> al corregirla")
    print(f"  el paralelismo SUBE a {100*n_full/N:.0f}%, pero es convergencia OBLIGATORIA al "
          f"colector (raiz 2deg,")
    print(f"  hoja 35deg): son confluencias en Y de tributarios, no vasos redundantes. La cola "
          f"corregible")
    print(f"  de verdad (paralelos en AMBOS sentidos) es solo {n_real} ({100*n_real/N:.0f}%), y "
          f"con hueco lateral -> (b)={b_sep}.")
    print("=" * 88)


if __name__ == "__main__":
    main()
