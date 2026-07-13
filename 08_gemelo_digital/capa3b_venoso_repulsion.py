#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa3b_venoso_repulsion.py  --  Bio-Kidney AI 2026
==================================================
Re-genera el arbol VENOSO con AUTO-SEPARACION por penalizacion de proximidad
durante el crecimiento (termino repulsivo contra los vasos ya colocados), y
barre el peso de repulsion W_REP contra un contrato de aceptacion.

ESCRITURA con respaldo, SIN promover a canonico: cada corrida sale a
capa3b_venoso_wrep{W}.npz. El capa3b_venoso.npz canonico NO se toca (ademas se
respaldo capa3b_venoso.py -> capa3b_venoso_v1_sinrepulsion_backup.py).

Pipeline por W_REP (identico salvo el peso de repulsion):
  1. Space colonization CONVERGENTE hacia los puntos de drenaje (atractores), con:
       - atraccion (media hacia atractores)  [original]
       - REPULSION nueva: cKDTree sobre puntos medios de aristas ya colocadas,
         radio de busqueda acotado; direccion final = norm(atraccion + W_REP*repulsion)
       - sesgo satelital vs arterial          [original, intacto]
       - binario GRADO_MAX_HIJOS=2, radios de Murray (FACTOR_RADIO_VENOSO=1.3,
         RADIO_TERMINAL_VENOSO=0.018), convergencia al hilio   [original, intacto]
  2. Separar el PUERTO venoso del arterial  (reusa parametros de capa3b_reparar_colisiones).
  3. Reparacion local AV  (reusa auditar_av/seg_seg_closest de capa3b_reparar_colisiones)
     -> el estado base (W_REP=0) reproduce capa3b_venoso_reparado.npz (AV~0, ~4070 shunts VV).
  4. Auditoria con los auditores existentes (cobertura, Murray, VV severidad, AV).

W_REP=0 es el CONTROL (repulsion off -> reproduce el arbol actual).

Solo NumPy / SciPy. SIN bpy. Entorno: env_biokidney.
"""

import os
import sys
import time
import numpy as np
from scipy.spatial import cKDTree

_here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
if _here not in sys.path:
    sys.path.insert(0, _here)

import capa3b_venoso as V                       # generador original (constantes + helpers)
import capa3b_reparar_colisiones as R           # puerto + reparacion AV (constantes + helpers)
from capa3b_clasificar_vv import clasificar
from capa3b_severidad_vv import critico

# ============================================================================
#  PARAMETROS DEL BARRIDO
# ============================================================================
W_REP_SWEEP = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0]
RADIO_REP   = 1.5      # mm  radio de busqueda de vasos cercanos (acotado, barato)
MIN_REP     = 0.45     # mm  ignora midpoints mas cercanos que esto (auto-repulsion del propio vaso)
W_SHUNT_MM  = 0.050    # umbral de shunt VV (50 um)

IN_PERI = "capa3ab_peritubular.npz"
IN_DOM  = "capa0_dominio.npz"
IN_ART  = "capa3a_arterial.npz"

# alias de constantes del generador original (intactas)
DIST_INFLUENCIA = V.DIST_INFLUENCIA
DIST_MATAR      = V.DIST_MATAR
PASO            = V.PASO_CRECIMIENTO
MAX_ITER        = V.MAX_ITERACIONES
GRADO_MAX       = V.GRADO_MAX_HIJOS
RADIO_TERMINAL  = V.RADIO_TERMINAL_VENOSO
FACTOR_RADIO    = V.FACTOR_RADIO_VENOSO
MURRAY_EXP      = V.MURRAY_EXP
REBUILD_CADA_K  = V.REBUILD_CADA_K
REBUILD_NUEVOS  = V.REBUILD_NUEVOS


# ============================================================================
#  REPULSION
# ============================================================================
def rep_dirs(pos, mid_tree, mids_arr, radio=RADIO_REP, min_rep=MIN_REP):
    """Direccion unitaria que ALEJA cada posicion de los puntos medios cercanos.

    pos (n,3). Para cada nodo: suma ponderada (mas cerca pesa mas) de los vectores
    (pos - midpoint) dentro de 'radio', ignorando los < min_rep (vaso propio).
    Devuelve (n,3) unitarios (0 donde no hay vecinos).
    """
    out = np.zeros((len(pos), 3))
    if mid_tree is None:
        return out
    listas = mid_tree.query_ball_point(pos, radio)
    for a, ids in enumerate(listas):
        if not ids:
            continue
        diff = pos[a] - mids_arr[ids]
        dd = np.linalg.norm(diff, axis=1)
        keep = dd > min_rep
        if not np.any(keep):
            continue
        diff = diff[keep]; dd = dd[keep]
        w = 1.0 - dd / radio                  # (0,1], mas cerca -> mas fuerte
        v = (w[:, None] * (diff / dd[:, None])).sum(axis=0)
        nv = np.linalg.norm(v)
        if nv > 1e-9:
            out[a] = v / nv
    return out


# ============================================================================
#  SPACE COLONIZATION CON REPULSION  (copia del original + termino repulsivo)
# ============================================================================
def space_colonization_rep(raiz, atractores, art_tree, w_rep):
    """Igual que capa3b_venoso.space_colonization, mas la REPULSION (gateada por
    w_rep>0, de modo que w_rep=0 reproduce EXACTAMENTE el original)."""
    atrac = np.asarray(atractores, dtype=np.float64).copy()

    cap = 1024
    pos_buf = np.empty((cap, 3)); parent_buf = np.empty(cap, np.int64)
    nhijos_buf = np.zeros(cap, np.int64); abierto_buf = np.zeros(cap, bool)
    pos_buf[0] = np.asarray(raiz, float); parent_buf[0] = -1; abierto_buf[0] = True
    M = 1
    mids = []                                   # puntos medios de aristas colocadas
    state = {"tree": None, "ids": None, "iters": 0, "nuevos": 0,
             "mid_tree": None, "mid_arr": None}

    def _ensure(n):
        nonlocal cap, pos_buf, parent_buf, nhijos_buf, abierto_buf
        if n <= cap:
            return
        while cap < n:
            cap *= 2
        pos_buf = np.resize(pos_buf, (cap, 3)); parent_buf = np.resize(parent_buf, cap)
        nhijos_buf = np.resize(nhijos_buf, cap); abierto_buf = np.resize(abierto_buf, cap)

    def agregar_hijo(p, pos):
        nonlocal M
        _ensure(M + 1)
        pos_buf[M] = pos; parent_buf[M] = p; nhijos_buf[M] = 0; abierto_buf[M] = True
        if w_rep > 0:
            mids.append(0.5 * (pos_buf[p] + pos))
        M += 1
        nhijos_buf[p] += 1
        if nhijos_buf[p] >= GRADO_MAX:
            abierto_buf[p] = False
        state["nuevos"] += 1

    def reconstruir():
        ids = np.nonzero(abierto_buf[:M])[0]
        state["tree"] = cKDTree(pos_buf[ids]); state["ids"] = ids
        state["iters"] = 0; state["nuevos"] = 0
        if w_rep > 0 and mids:
            state["mid_arr"] = np.asarray(mids)
            state["mid_tree"] = cKDTree(state["mid_arr"])

    def matar(atrac, nuevos_pos):
        if len(atrac) == 0 or len(nuevos_pos) == 0:
            return atrac
        dk, _ = cKDTree(np.asarray(nuevos_pos)).query(atrac, k=1)
        return atrac[dk > DIST_MATAR]

    atrac = matar(atrac, [pos_buf[0]])
    reconstruir()

    n_iter = 0
    for it in range(MAX_ITER):
        n_iter = it + 1
        if len(atrac) == 0:
            break
        if state["iters"] >= REBUILD_CADA_K or state["nuevos"] >= REBUILD_NUEVOS:
            reconstruir()
        state["iters"] += 1

        dist, rows = state["tree"].query(atrac, k=1)
        near = state["ids"][rows]
        en_rango = dist <= DIST_INFLUENCIA
        if not np.any(en_rango):
            a = int(np.argmin(dist)); nd = int(near[a])
            d = atrac[a] - pos_buf[nd]; dn = np.linalg.norm(d)
            if dn < 1e-9:
                break
            base = (d / dn)[None, :]
            sesg = V.sesgo_satelital(pos_buf[nd][None, :], base, art_tree)[0]
            nuevo = pos_buf[nd] + PASO * sesg
            agregar_hijo(nd, nuevo); atrac = matar(atrac, [nuevo]); reconstruir()
            continue

        activo = en_rango & abierto_buf[near]
        if not np.any(activo):
            reconstruir(); continue
        near_act = near[activo]
        vecs = atrac[activo] - pos_buf[near_act]
        norm = np.linalg.norm(vecs, axis=1, keepdims=True); norm[norm == 0] = 1.0
        vecs /= norm
        uniq, inv = np.unique(near_act, return_inverse=True)
        acc = np.zeros((len(uniq), 3)); np.add.at(acc, inv, vecs)
        gn = np.linalg.norm(acc, axis=1); keep = gn > 1e-9
        uniq = uniq[keep]; g = acc[keep] / gn[keep, None]
        if len(uniq) == 0:
            reconstruir(); continue

        # --- REPULSION: direccion final = norm(atraccion + W_REP * repulsion) ---
        if w_rep > 0:
            rep = rep_dirs(pos_buf[uniq], state["mid_tree"], state["mid_arr"])
            g = g + w_rep * rep
            gn2 = np.linalg.norm(g, axis=1, keepdims=True); gn2[gn2 < 1e-9] = 1.0
            g = g / gn2

        # --- sesgo satelital vs arterial (intacto) ---
        g = V.sesgo_satelital(pos_buf[uniq], g, art_tree)

        nuevos = pos_buf[uniq] + PASO * g
        nuevos_pos = []
        for k in range(len(uniq)):
            agregar_hijo(int(uniq[k]), nuevos[k]); nuevos_pos.append(nuevos[k])
        atrac = matar(atrac, nuevos_pos)

    return pos_buf[:M].copy(), parent_buf[:M].copy(), len(atrac), n_iter


# ============================================================================
#  PUERTO + REPARACION AV  (replica fiel de capa3b_reparar_colisiones, reusando sus helpers)
# ============================================================================
def separar_puerto(nodos, arist):
    K = len(nodos)
    parent = -np.ones(K, np.int64); parent[arist[:, 1]] = arist[:, 0]
    depth = np.zeros(K, np.int64)
    for i in range(1, K):
        if parent[i] >= 0:
            depth[i] = depth[parent[i]] + 1
    eje = {'X': 0, 'Y': 1, 'Z': 2}[R.EJE_PUERTO]
    offset = np.zeros(3); offset[eje] = R.SIGNO_PUERTO * R.OFFSET_PUERTO_MM
    w = np.clip(1.0 - depth / float(R.PROF_FADE_PUERTO), 0.0, 1.0)
    return nodos + w[:, None] * offset[None, :]


def reparar_av(nodos, arist, radios, Pa, Qa, radios_a):
    K = len(nodos)
    vecinos = [[] for _ in range(K)]
    for (u, v) in arist:
        vecinos[u].append(v); vecinos[v].append(u)
    for _ in range(R.MAX_PASADAS):
        Pv = nodos[arist[:, 0]]; Qv = nodos[arist[:, 1]]
        a = R.auditar_av(Pa, Qa, radios_a, Pv, Qv, radios, R.FACTOR_DETECCION)
        if a["n_col"] == 0:
            break
        orden = np.argsort(a["col_d"]); peor = {}
        for idx in orden:
            v = int(a["col_v"][idx])
            if v not in peor:
                peor[v] = int(a["col_a"][idx])
        disp = np.zeros((K, 3))
        for v, ai in peor.items():
            u0, u1 = int(arist[v, 0]), int(arist[v, 1])
            dd, c1, c2 = R.seg_seg_closest(Pa[ai][None, :], Qa[ai][None, :],
                                           nodos[u0][None, :], nodos[u1][None, :])
            d = float(dd[0]); direccion = c2[0] - c1[0]; nrm = np.linalg.norm(direccion)
            if nrm < R.EPS:
                direccion = 0.5 * (nodos[u0] + nodos[u1]) - 0.5 * (Pa[ai] + Qa[ai])
                nrm = np.linalg.norm(direccion)
                if nrm < R.EPS:
                    direccion = np.array([0.0, 0.0, 1.0]); nrm = 1.0
            direccion /= nrm
            need = (radios_a[ai] + radios[v]) * R.MARGEN
            empuje = need - d
            if empuje <= 0:
                continue
            vec = direccion * min(empuje + 1e-3, R.MAX_EMPUJE_MM)
            for nodo in (u0, u1):
                disp[nodo] += vec
            seed = {u0, u1}
            hop1 = {w2 for s in seed for w2 in vecinos[s]} - seed
            hop2 = {w2 for s in hop1 for w2 in vecinos[s]} - seed - hop1
            for nodo in hop1:
                disp[nodo] += vec * R.SUAVE_PESOS[0]
            for nodo in hop2:
                disp[nodo] += vec * R.SUAVE_PESOS[1]
        nodos = nodos + disp
    return nodos


# ============================================================================
#  AUDITORIA DE UNA CORRIDA
# ============================================================================
def murray_pct(parent, r_in):
    K = len(parent)
    hijos = [[] for _ in range(K)]
    for i in range(K):
        if parent[i] >= 0:
            hijos[parent[i]].append(i)
    ok = tot = 0
    for i in range(K):
        if hijos[i]:
            tot += 1
            lhs = r_in[i] ** MURRAY_EXP
            rhs = sum(r_in[c] ** MURRAY_EXP for c in hijos[i])
            if abs(lhs - rhs) <= 1e-3 * max(lhs, 1e-12):
                ok += 1
    return (100.0 * ok / tot if tot else 100.0), hijos


def generar(w_rep, puntos_drenaje, hilio, art_tree, Pa, Qa, radios_a, out_dir):
    """Genera + puerto + reparacion AV; guarda npz; devuelve dict del arbol."""
    nodos, parent, n_no_alc, n_iter = space_colonization_rep(
        hilio, puntos_drenaje, art_tree, w_rep)
    r_in, _ = V.asignar_radios_murray(parent, RADIO_TERMINAL)
    r_in = r_in * FACTOR_RADIO

    edges = [(parent[i], i) for i in range(len(nodos)) if parent[i] >= 0]
    arist = np.asarray(edges, np.int64)
    radios = np.asarray([r_in[i] for (_, i) in edges], float)
    hijos_cnt = np.bincount(arist[:, 0], minlength=len(nodos))
    terminales = np.where(hijos_cnt == 0)[0].astype(np.int64)

    # puerto + reparacion AV
    nodos = separar_puerto(nodos, arist)
    nodos = reparar_av(nodos, arist, radios, Pa, Qa, radios_a)

    # guardar (NO canonico)
    tag = f"{w_rep:g}".replace(".", "p")
    out_path = os.path.join(out_dir, f"capa3b_venoso_wrep{tag}.npz")
    np.savez_compressed(
        out_path,
        nodos=nodos.astype(np.float32), aristas=arist.astype(np.int32),
        radios=radios.astype(np.float32), terminales=terminales.astype(np.int32),
        raiz=np.int32(0), raiz_pos=nodos[0].astype(np.float64),
        w_rep=np.float64(w_rep), radio_rep=np.float64(RADIO_REP),
        factor_radio_venoso=np.float64(FACTOR_RADIO),
        radio_terminal_venoso=np.float64(RADIO_TERMINAL),
        n_iteraciones=np.int32(n_iter),
    )
    return dict(nodos=nodos, arist=arist, radios=radios, terminales=terminales,
                parent=parent, r_in=r_in, n_no_alc=n_no_alc, out=out_path)


def auditar_run(tree, puntos_drenaje, Pa, Qa, radios_a):
    nodos, arist, radios = tree["nodos"], tree["arist"], tree["radios"]
    # cobertura
    kd = cKDTree(nodos)
    d_at, _ = kd.query(puntos_drenaje, k=1)
    cov = 100.0 * np.count_nonzero(d_at < DIST_MATAR) / len(puntos_drenaje)
    # Murray
    mur, hijos = murray_pct(tree["parent"], tree["r_in"])
    # grados de salida
    out_deg = np.array([len(h) for h in hijos])
    n_branch = int(np.count_nonzero(out_deg >= 1))
    n_bin = int(np.count_nonzero(out_deg == 2))
    pct_bin = 100.0 * n_bin / n_branch if n_branch else 100.0
    megahubs = int(np.count_nonzero(out_deg > GRADO_MAX))
    # VV shunts >50um (cuadrante critico)
    cl = clasificar(nodos, arist, radios)
    cr = critico(cl)
    vv_shunt = int(np.count_nonzero(cr["p"] > W_SHUNT_MM)) if cr["n"] else 0
    # AV colisiones vs arterial
    Pv = nodos[arist[:, 0]]; Qv = nodos[arist[:, 1]]
    av = R.auditar_av(Pa, Qa, radios_a, Pv, Qv, radios, R.FACTOR_DETECCION)
    return dict(cov=cov, mur=mur, n_nodos=len(nodos), n_term=len(tree["terminales"]),
                pct_bin=pct_bin, megahubs=megahubs, vv_shunt=vv_shunt, av=av["n_col"])


# ============================================================================
#  MAIN
# ============================================================================
def main():
    base = os.path.dirname(V.find_npz(IN_PERI))
    peri = np.load(V.find_npz(IN_PERI), allow_pickle=False)
    puntos_drenaje = peri["puntos_drenaje"].astype(np.float64)
    d0 = np.load(V.find_npz(IN_DOM), allow_pickle=False)
    hilio = d0["hilio"].astype(np.float64)
    art = np.load(V.find_npz(IN_ART), allow_pickle=False)
    nodos_a = art["nodos"].astype(np.float64); arist_a = art["aristas"].astype(np.int64)
    radios_a = art["radios"].astype(np.float64)
    Pa = nodos_a[arist_a[:, 0]]; Qa = nodos_a[arist_a[:, 1]]
    art_tree = cKDTree(nodos_a)

    # referencia fija: shunts VV >50um del ARTERIAL
    cl_a = clasificar(nodos_a, arist_a, radios_a)
    cr_a = critico(cl_a)
    vv_art_ref = int(np.count_nonzero(cr_a["p"] > W_SHUNT_MM)) if cr_a["n"] else 0

    print("=" * 100)
    print("  CAPA 3b - BARRIDO DE REPULSION (auto-separacion durante el crecimiento)")
    print("=" * 100)
    print(f"  Atractores (drenaje): {len(puntos_drenaje)}   RADIO_REP={RADIO_REP} mm   "
          f"MIN_REP={MIN_REP} mm   shunt>{W_SHUNT_MM*1000:.0f}um")
    print(f"  Referencia fija - shunts VV>50um del ARTERIAL: {vv_art_ref}")
    print("-" * 100)
    hdr = (f"  {'W_REP':>6} | {'cov%':>6} {'Murray%':>8} | {'nodos':>7} {'term':>6} "
           f"{'bin%':>6} {'mhub':>5} | {'VV>50um':>8} {'(art ref)':>9} | {'AV':>4} | {'t(s)':>6}")
    print(hdr)
    print("  " + "-" * 96)

    filas = []
    for w in W_REP_SWEEP:
        t0 = time.perf_counter()
        tree = generar(w, puntos_drenaje, hilio, art_tree, Pa, Qa, radios_a, base)
        m = auditar_run(tree, puntos_drenaje, Pa, Qa, radios_a)
        dt = time.perf_counter() - t0
        filas.append((w, m, dt, tree["out"]))
        print(f"  {w:>6.2f} | {m['cov']:>6.1f} {m['mur']:>8.1f} | "
              f"{m['n_nodos']:>7d} {m['n_term']:>6d} {m['pct_bin']:>6.1f} {m['megahubs']:>5d} | "
              f"{m['vv_shunt']:>8d} {vv_art_ref:>9d} | {m['av']:>4d} | {dt:>6.1f}")
    print("  " + "-" * 96)

    # ------------------------------------------------------------------
    #  VEREDICTO
    # ------------------------------------------------------------------
    print("-" * 100)
    print("  VEREDICTO")
    validos = [(w, m, dt, p) for (w, m, dt, p) in filas
               if m["cov"] >= 99.99 and m["mur"] >= 99.99 and m["av"] <= 2
               and m["megahubs"] == 0]
    base_vv = filas[0][1]["vv_shunt"]
    print(f"  Base (W_REP=0): cobertura {filas[0][1]['cov']:.1f}%, Murray {filas[0][1]['mur']:.1f}%, "
          f"AV {filas[0][1]['av']}, shunts VV>50um {base_vv}")
    if validos:
        # entre los que cumplen contrato, el de menor VV; desempate por menor W
        mejor = min(validos, key=lambda x: (x[1]["vv_shunt"], x[0]))
        w, m, dt, p = mejor
        baja = base_vv - m["vv_shunt"]
        print(f"  Candidatos que mantienen cov=100% & Murray=100% & AV~0 & sin mega-hubs: "
              f"{[f'{v[0]:g}' for v in validos]}")
        print(f"  -> MEJOR candidato a canonico: W_REP={w:g}  "
              f"(shunts VV>50um {m['vv_shunt']} vs base {base_vv}, baja {baja}; "
              f"ref arterial {vv_art_ref})")
        if m["vv_shunt"] > 2 * vv_art_ref:
            print(f"     NOTA: aun lejos del nivel arterial ({vv_art_ref}); mejora parcial, "
                  f"no resuelve del todo.")
    else:
        print("  NINGUN W_REP mantiene cobertura 100% + Murray 100% + AV~0 simultaneamente.")
        # mostrar el trade-off: cobertura vs shunts
        for (w, m, dt, p) in filas:
            print(f"     W={w:g}: cov {m['cov']:.1f}%  VV>50um {m['vv_shunt']}  AV {m['av']}  "
                  f"mhub {m['megahubs']}")
        print("  -> TENSION ESTRUCTURAL: la repulsion necesaria para separar tambien impide "
              "alcanzar atractores (cae la cobertura). Hay que discutir ARQUITECTURA, no solo el peso.")
    print(f"  (todos los .npz por W_REP quedan en {base} para inspeccion; nada promovido)")
    print("=" * 100)


if __name__ == "__main__":
    main()
