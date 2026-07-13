#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa3b_reparar_colisiones.py  --  Bio-Kidney AI 2026
====================================================
Reparacion LOCAL de las colisiones arterial-venosas detectadas por
capa3b_auditoria_colisiones.py, SIN re-generar el arbol venoso.

Solo se mueven NODOS venosos implicados; se preservan EXACTAMENTE la topologia
(mismas aristas), la conectividad (arbol enraizado) y los radios de Murray (que
dependen de la estructura, no de las posiciones). El arbol arterial es FIJO.

PASO 1 - SEPARAR EL PUERTO DEL HILIO
    El venoso nace en el mismo punto que el arterial -> solape 0.0000 mm en la raiz.
    Se desplaza la RAIZ venosa a un PUERTO separado (hilio + offset, ~3 mm en
    direccion anterior; anatomicamente la vena va por delante de la arteria) y se
    PROPAGA el desplazamiento a los primeros nodos del tronco con peso decreciente
    con la profundidad (graph-depth), para no romper la convergencia del tronco.

PASO 2 - REPARACION LOCAL DE CRUCES EN PARENQUIMA
    Para cada segmento venoso en colision: se halla el segmento arterial que lo
    cruza, la direccion de separacion (del punto mas cercano arterial al venoso) y
    se empuja el segmento venoso lo MINIMO para que la distancia supere
    (r_art + r_ven) * MARGEN (MARGEN=1.5: holgura de pared, no tangencia). Los
    nodos vecinos (1-2 saltos) se desplazan con peso decreciente (suavizado, sin
    quiebres). Topologia y radios intactos.

PASO 3 - RE-VERIFICACION
    Se re-corre la logica de auditoria (reutiliza seg_seg_dist de la auditoria):
    colisiones arterial-venosas restantes y, ademas, colisiones NUEVAS
    venoso-venoso introducidas por los empujes. Maximo MAX_PASADAS reparaciones;
    el residual se reporta (sin bucle infinito).

Solo NumPy / SciPy. SIN bpy. Entorno: env_biokidney.

Entrada : capa3a_arterial.npz (fijo), capa3b_venoso.npz,
          capa3b_auditoria_colisiones.npz (referencia de partida).
Salida  : capa3b_venoso_reparado.npz   (NO sobrescribe el original).
"""

import os
import sys
import itertools
import numpy as np
from scipy.spatial import cKDTree

# ============================================================================
#  PARAMETROS  (editar aqui)
# ============================================================================
# --- PASO 1: puerto del hilio ---
OFFSET_PUERTO_MM   = 3.0     # separacion del puerto venoso respecto al arterial
EJE_PUERTO         = 'Z'     # eje anterior-posterior (el mas fino, semieje ~18 mm):
SIGNO_PUERTO       = +1.0    # la vena va por delante (anterior) de la arteria
PROF_FADE_PUERTO   = 12      # profundidad (n.º de aristas) en que se desvanece el offset

# --- PASO 2: reparacion local ---
FACTOR_DETECCION   = 1.0     # criterio de COLISION (reproduce el 67 de la auditoria)
MARGEN             = 1.5     # objetivo de separacion al reparar: (r_a+r_v)*MARGEN
SUAVE_PESOS        = (0.5, 0.25)   # pesos del suavizado a 1 y 2 saltos del segmento
# MAX_PASADAS 3->4 (jul 2026, cascada corticalseno): la reparacion local OSCILA
# (empujar un segmento + su suavizado a 1-2 saltos crea colisiones transitorias en los
# vecinos, que se reparan y el segmento re-colisiona una vez mas). Con la geometria
# regenerada (cortex 6.6mm) esa oscilacion converge en la 4a pasada; el tope de 3
# cortaba justo UNA pasada antes (residual 1). Se sube a 4. NO se tocan MARGEN/SUAVE/
# OFFSET ni el factor de deteccion: la fisica del reparador es identica. Cada pasada
# solo empuja segmentos AUN en colision, asi que las reparaciones ya buenas no se
# re-tocan -> 0 regresiones (diagnostico: 33->3->6->1->0).
MAX_PASADAS        = 4       # tope de pasadas de reparacion
MAX_EMPUJE_MM      = 3.0     # clamp de seguridad por nodo y pasada

# --- candidatos (cKDTree) ---
UMBRAL_CAND_AV     = 5.0     # mm  candidatos arterial-venoso (igual que la auditoria)
UMBRAL_CAND_VV     = 2.5     # mm  candidatos venoso-venoso (colision vv exige <(2*rmax)*MARGEN)
CHUNK              = 2_000_000

EPS = 1e-9

IN_ART = "capa3a_arterial.npz"
IN_VEN = "capa3b_venoso.npz"
IN_AUD = "capa3b_auditoria_colisiones.npz"
OUT_NPZ = "capa3b_venoso_reparado.npz"


# ============================================================================
#  REUTILIZACION de la auditoria
# ============================================================================
_here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
if _here not in sys.path:
    sys.path.insert(0, _here)
try:
    from capa3b_auditoria_colisiones import seg_seg_dist, find_npz   # reuso de la logica
except Exception:  # respaldo si no se puede importar
    def find_npz(name):
        for c in (
            os.path.expanduser(os.path.join("~", "Escritorio", "BioKidney-AI", name)),
            os.path.join(_here, name), os.path.join(_here, "..", name),
            os.path.join(os.getcwd(), name),
        ):
            if os.path.isfile(c):
                return os.path.abspath(c)
        raise FileNotFoundError(name)
    seg_seg_dist = None


def seg_seg_closest(P1, Q1, P2, Q2):
    """Como seg_seg_dist pero devuelve tambien los puntos mas cercanos.

    Devuelve (dist (N,), c1 (N,3) en el segmento arterial, c2 (N,3) en el venoso).
    Mismo algoritmo de Ericson que la auditoria; aqui necesitamos c1,c2 para la
    DIRECCION de separacion del empuje.
    """
    d1 = Q1 - P1; d2 = Q2 - P2; r = P1 - P2
    a = np.einsum("ij,ij->i", d1, d1)
    e = np.einsum("ij,ij->i", d2, d2)
    f = np.einsum("ij,ij->i", d2, r)
    c = np.einsum("ij,ij->i", d1, r)
    b = np.einsum("ij,ij->i", d1, d2)
    a_safe = np.where(a > EPS, a, 1.0)
    e_safe = np.where(e > EPS, e, 1.0)
    denom = a * e - b * b
    s = np.where(denom > EPS, np.clip((b * f - c * e) / np.where(denom > EPS, denom, 1.0),
                                      0.0, 1.0), 0.0)
    t = (b * s + f) / e_safe
    t_lt = t < 0.0; t_gt = t > 1.0
    t = np.clip(t, 0.0, 1.0)
    s = np.where(t_lt, np.clip(-c / a_safe, 0.0, 1.0), s)
    s = np.where(t_gt, np.clip((b - c) / a_safe, 0.0, 1.0), s)
    a_deg = a <= EPS
    s = np.where(a_deg, 0.0, s); t = np.where(a_deg, np.clip(f / e_safe, 0.0, 1.0), t)
    e_deg = e <= EPS
    t = np.where(e_deg, 0.0, t); s = np.where(e_deg, np.clip(-c / a_safe, 0.0, 1.0), s)
    c1 = P1 + d1 * s[:, None]
    c2 = P2 + d2 * t[:, None]
    diff = c1 - c2
    dist = np.sqrt(np.einsum("ij,ij->i", diff, diff))
    return dist, c1, c2


if seg_seg_dist is None:   # respaldo: distancia a partir de seg_seg_closest
    def seg_seg_dist(P1, Q1, P2, Q2):
        d, c1, c2 = seg_seg_closest(P1, Q1, P2, Q2)
        return d, 0.5 * (c1 + c2)


# ============================================================================
#  AUDITORIAS (arterial-venoso y venoso-venoso)
# ============================================================================
def auditar_av(Pa, Qa, ra, Pv, Qv, rv, factor, umbral=UMBRAL_CAND_AV):
    """Colisiones arterial-venoso. Devuelve dict con pares, distancias y min global."""
    mid_a = 0.5 * (Pa + Qa); mid_v = 0.5 * (Pv + Qv)
    tree = cKDTree(mid_a)
    neigh = tree.query_ball_point(mid_v, r=umbral, workers=-1)
    counts = np.fromiter((len(x) for x in neigh), np.int64, len(neigh))
    n_cand = int(counts.sum())
    vidx = np.repeat(np.arange(len(mid_v), dtype=np.int64), counts)
    aidx = np.fromiter(itertools.chain.from_iterable(neigh), np.int64, n_cand)
    col_a, col_v, col_d = [], [], []
    dmin = np.inf
    for s0 in range(0, n_cand, CHUNK):
        s1 = min(s0 + CHUNK, n_cand)
        ai = aidx[s0:s1]; vi = vidx[s0:s1]
        dd, _ = seg_seg_dist(Pa[ai], Qa[ai], Pv[vi], Qv[vi])
        if dd.size:
            dmin = min(dmin, float(dd.min()))
        thr = (ra[ai] + rv[vi]) * factor
        hit = dd < thr
        if np.any(hit):
            col_a.append(ai[hit]); col_v.append(vi[hit]); col_d.append(dd[hit])
    if col_a:
        col_a = np.concatenate(col_a); col_v = np.concatenate(col_v); col_d = np.concatenate(col_d)
    else:
        col_a = np.empty(0, np.int64); col_v = np.empty(0, np.int64); col_d = np.empty(0)
    return dict(n_cand=n_cand, col_a=col_a, col_v=col_v, col_d=col_d,
                n_col=len(col_a), dmin=dmin)


def auditar_vv(P, Q, r, arist, factor, umbral=UMBRAL_CAND_VV):
    """Colisiones venoso-venoso, excluyendo segmentos ADYACENTES (comparten nodo).

    Devuelve (n_col, dmin_no_adyacente).
    """
    mid = 0.5 * (P + Q)
    tree = cKDTree(mid)
    pairs = tree.query_pairs(r=umbral, output_type='ndarray')   # (M,2) con i<j
    if len(pairs) == 0:
        return 0, np.inf
    i, j = pairs[:, 0], pairs[:, 1]
    # excluir adyacentes: comparten algun nodo extremo
    comparten = (
        (arist[i, 0] == arist[j, 0]) | (arist[i, 0] == arist[j, 1]) |
        (arist[i, 1] == arist[j, 0]) | (arist[i, 1] == arist[j, 1])
    )
    i = i[~comparten]; j = j[~comparten]
    if len(i) == 0:
        return 0, np.inf
    n_col = 0
    dmin = np.inf
    for s0 in range(0, len(i), CHUNK):
        s1 = min(s0 + CHUNK, len(i))
        ii = i[s0:s1]; jj = j[s0:s1]
        dd, _ = seg_seg_dist(P[ii], Q[ii], P[jj], Q[jj])
        dmin = min(dmin, float(dd.min()))
        thr = (r[ii] + r[jj]) * factor
        n_col += int(np.count_nonzero(dd < thr))
    return n_col, dmin


# ============================================================================
#  MAIN
# ============================================================================
def main():
    art_path = find_npz(IN_ART)
    ven_path = find_npz(IN_VEN)
    aud_path = find_npz(IN_AUD)

    art = np.load(art_path, allow_pickle=False)
    ven = dict(np.load(ven_path, allow_pickle=False))
    aud = np.load(aud_path, allow_pickle=False)

    nodos_a = art["nodos"].astype(np.float64)
    arist_a = art["aristas"].astype(np.int64)
    radios_a = art["radios"].astype(np.float64)
    Pa = nodos_a[arist_a[:, 0]]; Qa = nodos_a[arist_a[:, 1]]

    nodos_v0 = ven["nodos"].astype(np.float64)          # ORIGINAL (no se toca)
    arist_v = ven["aristas"].astype(np.int64)
    radios_v = ven["radios"].astype(np.float64)
    K = len(nodos_v0)
    hilio = ven["raiz_pos"].astype(np.float64) if "raiz_pos" in ven else nodos_v0[0].copy()

    nodos_v = nodos_v0.copy()                            # copia de trabajo

    print("=" * 70)
    print("  CAPA 3b - REPARACION LOCAL DE COLISIONES (venoso)")
    print("=" * 70)
    print(f"  Arterial (fijo): {art_path}")
    print(f"  Venoso         : {ven_path}")
    print(f"  Auditoria ref. : {aud_path}")
    print(f"  Nodos venosos {K}   aristas {len(arist_v)}   (arbol: {len(arist_v)==K-1})")
    print(f"  Colisiones de referencia (auditoria): {int(aud['n_colisiones'])}  "
          f"(sep. min {float(aud['dist_min_global']):.4f} mm)")
    print("-" * 70)

    def segs_v():
        return nodos_v[arist_v[:, 0]], nodos_v[arist_v[:, 1]]

    # --- baseline propio (debe reproducir la auditoria) ---------------------
    Pv0, Qv0 = segs_v()
    base = auditar_av(Pa, Qa, radios_a, Pv0, Qv0, radios_v, FACTOR_DETECCION)
    vv_base, vv_base_dmin = auditar_vv(Pv0, Qv0, radios_v, arist_v, FACTOR_DETECCION)
    n_col_antes = base["n_col"]; dmin_antes = base["dmin"]
    print(f"  [baseline] colisiones AV = {n_col_antes}   sep.min AV = {dmin_antes:.4f} mm")
    print(f"  [baseline] colisiones VV (venoso-venoso) = {vv_base}   "
          f"sep.min VV = {vv_base_dmin:.4f} mm")
    print("-" * 70)

    # ========================================================================
    #  PASO 1 - separar el puerto del hilio
    # ========================================================================
    # profundidad (n.º de aristas desde la raiz). parent[i] < i en capa3b.
    parent = -np.ones(K, dtype=np.int64)
    parent[arist_v[:, 1]] = arist_v[:, 0]
    depth = np.zeros(K, dtype=np.int64)
    for i in range(1, K):
        depth[i] = depth[parent[i]] + 1

    eje = {'X': 0, 'Y': 1, 'Z': 2}[EJE_PUERTO]
    offset = np.zeros(3); offset[eje] = SIGNO_PUERTO * OFFSET_PUERTO_MM
    w_port = np.clip(1.0 - depth / float(PROF_FADE_PUERTO), 0.0, 1.0)   # 1 en raiz -> 0
    disp_port = w_port[:, None] * offset[None, :]
    nodos_v = nodos_v + disp_port
    puerto_venoso = nodos_v[0].copy()
    n_port = int(np.count_nonzero(w_port > 0))
    print(f"  PASO 1: puerto venoso desplazado a {puerto_venoso}  "
          f"(offset {OFFSET_PUERTO_MM} mm en {EJE_PUERTO})")
    print(f"          nodos del tronco desplazados (fade depth<{PROF_FADE_PUERTO}): {n_port}")

    Pv1, Qv1 = segs_v()
    a1 = auditar_av(Pa, Qa, radios_a, Pv1, Qv1, radios_v, FACTOR_DETECCION)
    print(f"          tras PASO 1: colisiones AV = {a1['n_col']}   sep.min AV = {a1['dmin']:.4f} mm")
    print("-" * 70)

    # ========================================================================
    #  PASO 2 - reparacion local (hasta MAX_PASADAS)
    # ========================================================================
    # adyacencia para el suavizado
    vecinos = [[] for _ in range(K)]
    for (u, v) in arist_v:
        vecinos[u].append(v); vecinos[v].append(u)

    nodos_reparados = set()
    print(f"  PASO 2: reparacion local (MARGEN={MARGEN}, max {MAX_PASADAS} pasadas)")
    for passo in range(1, MAX_PASADAS + 1):
        Pv, Qv = segs_v()
        a = auditar_av(Pa, Qa, radios_a, Pv, Qv, radios_v, FACTOR_DETECCION)
        if a["n_col"] == 0:
            print(f"     pasada {passo}: 0 colisiones AV -> nada que reparar, fin.")
            break

        # peor socio arterial (min distancia) por cada segmento venoso en colision
        orden = np.argsort(a["col_d"])
        peor = {}
        for idx in orden:
            v = int(a["col_v"][idx])
            if v not in peor:
                peor[v] = (int(a["col_a"][idx]), float(a["col_d"][idx]))

        disp = np.zeros((K, 3), dtype=np.float64)
        for v, (ai, _) in peor.items():
            u0, u1 = int(arist_v[v, 0]), int(arist_v[v, 1])
            dd, c1, c2 = seg_seg_closest(
                Pa[ai][None, :], Qa[ai][None, :],
                nodos_v[u0][None, :], nodos_v[u1][None, :])
            d = float(dd[0])
            direccion = (c2[0] - c1[0])
            nrm = np.linalg.norm(direccion)
            if nrm < EPS:   # solape exacto: usa la direccion entre puntos medios
                direccion = (0.5 * (nodos_v[u0] + nodos_v[u1]) - 0.5 * (Pa[ai] + Qa[ai]))
                nrm = np.linalg.norm(direccion)
                if nrm < EPS:
                    direccion = np.array([0.0, 0.0, 1.0]); nrm = 1.0
            direccion = direccion / nrm
            need = (radios_a[ai] + radios_v[v]) * MARGEN
            empuje = need - d
            if empuje <= 0:
                continue
            empuje = min(empuje + 1e-3, MAX_EMPUJE_MM)
            vec = direccion * empuje
            # empuje pleno a los dos nodos del segmento
            for nodo in (u0, u1):
                disp[nodo] += vec
                nodos_reparados.add(nodo)
            # suavizado: vecinos a 1 y 2 saltos con peso decreciente
            seed = {u0, u1}
            hop1 = {w for s in seed for w in vecinos[s]} - seed
            hop2 = {w for s in hop1 for w in vecinos[s]} - seed - hop1
            for nodo in hop1:
                disp[nodo] += vec * SUAVE_PESOS[0]
            for nodo in hop2:
                disp[nodo] += vec * SUAVE_PESOS[1]

        nodos_v = nodos_v + disp
        Pv, Qv = segs_v()
        a_post = auditar_av(Pa, Qa, radios_a, Pv, Qv, radios_v, FACTOR_DETECCION)
        mov = int(np.count_nonzero(np.linalg.norm(disp, axis=1) > EPS))
        print(f"     pasada {passo}: reparados {len(peor)} segmentos venosos, "
              f"{mov} nodos movidos -> colisiones AV restantes = {a_post['n_col']}  "
              f"(sep.min {a_post['dmin']:.4f} mm)")
        if a_post["n_col"] == 0:
            break

    # ========================================================================
    #  PASO 3 - re-verificacion final
    # ========================================================================
    Pv, Qv = segs_v()
    fin_av = auditar_av(Pa, Qa, radios_a, Pv, Qv, radios_v, FACTOR_DETECCION)
    vv_fin, vv_fin_dmin = auditar_vv(Pv, Qv, radios_v, arist_v, FACTOR_DETECCION)

    # desplazamientos totales aplicados (PASO1 + PASO2) respecto al original
    disp_total = nodos_v - nodos_v0
    mag = np.linalg.norm(disp_total, axis=1)
    movidos = mag > EPS
    n_movidos = int(np.count_nonzero(movidos))
    desp_medio = float(mag[movidos].mean()) if n_movidos else 0.0
    desp_max = float(mag.max()) if n_movidos else 0.0

    # Murray intacto: los radios no se tocaron
    murray_ok = np.array_equal(radios_v, ven["radios"].astype(np.float64))
    conect_ok = (len(arist_v) == K - 1)

    # --- guardar (mismo esquema que capa3b_venoso, nodos reparados) ---------
    out = {k: ven[k] for k in ven}     # copia todos los campos originales
    out["nodos"] = nodos_v.astype(np.float32)
    out["raiz_pos"] = puerto_venoso.astype(np.float64)
    # metadata de reparacion
    out["rep_offset_puerto_mm"] = np.float64(OFFSET_PUERTO_MM)
    out["rep_eje_puerto"] = np.str_(EJE_PUERTO)
    out["rep_margen"] = np.float64(MARGEN)
    out["rep_n_nodos_movidos"] = np.int32(n_movidos)
    out["rep_desp_medio_mm"] = np.float64(desp_medio)
    out["rep_desp_max_mm"] = np.float64(desp_max)
    out["rep_colisiones_antes"] = np.int32(n_col_antes)
    out["rep_colisiones_despues"] = np.int32(fin_av["n_col"])
    out["rep_sepmin_antes_mm"] = np.float64(dmin_antes)
    out["rep_sepmin_despues_mm"] = np.float64(fin_av["dmin"])
    out["rep_vv_antes"] = np.int32(vv_base)
    out["rep_vv_despues"] = np.int32(vv_fin)
    out_path = os.path.join(os.path.dirname(ven_path), OUT_NPZ)
    np.savez_compressed(out_path, **out)

    # ========================================================================
    #  VERIFICACION (consola)
    # ========================================================================
    print("-" * 70)
    print("  VERIFICACION FINAL")
    print("-" * 70)
    print(f"  Colisiones AV  antes -> despues : {n_col_antes}  ->  {fin_av['n_col']}")
    print(f"  Separacion min AV antes -> despues: {dmin_antes:.4f}  ->  {fin_av['dmin']:.4f} mm  "
          f"{'(sana, >0)' if fin_av['dmin'] > 0 else '(REVISAR: sigue en 0)'}")
    print(f"  Colisiones VV  antes -> despues : {vv_base}  ->  {vv_fin}   "
          f"{'(sin nuevas vv)' if vv_fin <= vv_base else 'REVISAR: aparecieron vv'}")
    print(f"  Nodos movidos (PASO1 tronco + PASO2 reparacion): {n_movidos} / {K}  "
          f"({100.0*n_movidos/K:.2f} %)")
    print(f"     PASO 2 (reparacion local) movio: {len(nodos_reparados)} nodos")
    print(f"  Desplazamiento aplicado: medio {desp_medio:.3f} mm   max {desp_max:.3f} mm  "
          f"{'(pequeno/local)' if desp_max <= MAX_EMPUJE_MM + OFFSET_PUERTO_MM + 1e-6 else 'REVISAR'}")
    print(f"  Conectividad intacta (aristas = nodos-1): {conect_ok}")
    print(f"  Radios Murray sin cambios: {murray_ok}")
    print("-" * 70)
    residual = fin_av["n_col"]
    if residual == 0 and vv_fin <= vv_base and fin_av["dmin"] > 0:
        print(f"  RESULTADO: OK - 0 colisiones AV, separacion sana "
              f"({fin_av['dmin']:.4f} mm), sin nuevas vv, Murray y topologia intactos.")
    else:
        print(f"  RESULTADO: residual = {residual} colisiones AV tras {MAX_PASADAS} pasadas; "
              f"vv {vv_base}->{vv_fin}. Revisar parametros (MARGEN/SUAVE/OFFSET).")
    print(f"  -> guardado en: {out_path}  (NO se sobrescribio el original)")
    print("=" * 70)


if __name__ == "__main__":
    main()
