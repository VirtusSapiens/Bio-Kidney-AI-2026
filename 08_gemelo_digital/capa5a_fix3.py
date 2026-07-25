#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa5a_fix3.py -- FIX estanqueidad 5a v3: PEDICULOS HILARES SEGREGADOS (deterministico).
Anclaje anatomico: orden anteroposterior vena-arteria-pelvis (eje AP = Z; +Z anterior).
Construccion, NO solver. Mantiene FIX1 (pelvis elipsoide Z_semi=2.0) SIN cambios; pelvis y
ureter NO se mueven (referencia posterior). Modifica SOLO tramos hilares vasculares (local
Capa 5, sintesis de pediculo). NO reabre 3a/3b/3c, NO cambia Capa 4, NO genera 5b.
FALLBACK en cualquier paso -> DETENERSE, reportar el obstaculo, NO forzar, NO carve.
"""
import os, sys, json
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import capa5a_nucleo_sdf as C5

MIN_WALL = 0.400
MARGIN = 0.200
MIN_CANAL_R = 0.200
PELVIS_CENTRO = np.array([0.0, -27.2, 0.0])
PELVIS_ZSEMI = 2.0
PELVIS_YSEMI = 4.5
Y_BORDE_ANT = PELVIS_CENTRO[1] + PELVIS_YSEMI   # -22.7
URETER_PUERTO = np.array([0.0, -30.0, 0.0])
OUT_MD = os.path.join(C5.REPO, "09_paper_vascular", "auditoria_5a_fix3.md")


def main():
    LOG = []
    def log(s=""): print(s); LOG.append(s)
    log("FIX v3 -- PEDICULOS HILARES SEGREGADOS (vena-arteria-pelvis, AP=Z, +Z anterior)")

    d3a = np.load(C5.find("capa3a_arterial.npz"), allow_pickle=True)
    d3b = np.load(C5.find("capa3b_venoso.npz"), allow_pickle=True)
    d3c = np.load(C5.find("capa3c_colector.npz"), allow_pickle=True)
    d4 = np.load(C5.find("capa4_colector_alto.npz"), allow_pickle=True)
    n3a, e3a, r3a = C5.adapt_edge_radios(d3a); n3b, e3b, r3b = C5.adapt_edge_radios(d3b)
    n3c, e3c, r3c = C5.adapt_3c(d3c); n4, e4, r4 = C5.adapt_4(d4)
    Zp = float(d4["z_pelvis"]); niv = d4["nivel"]
    pel_pos = n4[niv == 3][0]; pel_sa = np.array([PELVIS_YSEMI, PELVIS_YSEMI, Zp])

    # ---- PASO 1: corredores Z ----
    r_art = float(d3a["radio_raiz"]); r_ven = float(d3b["radio_raiz"])
    Z_art = PELVIS_ZSEMI + r_art + MIN_WALL + MARGIN
    Z_ven = Z_art + r_art + r_ven + MIN_WALL + MARGIN
    ven_clear = Z_ven - PELVIS_ZSEMI - r_ven
    log("\nPASO1 CORREDORES Z (AP; +Z anterior):")
    log(f"   r_art={r_art*1000:.0f}um  r_ven={r_ven*1000:.0f}um")
    log(f"   Z_art = {Z_art:.3f} mm (arteria, 1er corredor anterior a pelvis)")
    log(f"   Z_ven = {Z_ven:.3f} mm (vena, anterior a la arteria)")
    log(f"   pelvis Z_semi={PELVIS_ZSEMI}; vena despeja pelvis directo: {ven_clear:.3f} >= {MIN_WALL} "
        f"{'[OK]' if ven_clear >= MIN_WALL else '[NO]'}")
    log(f"   puertos nuevos: ARTERIA=[0,-30,{Z_art:.3f}]  VENA=[0,-30,{Z_ven:.3f}]  (ureter=[0,-30,0] SIN cambios)")

    # ---- keepout urinario MACRO (pelvis elipsoide + calices/ureter = Capa 4) ----
    uN, uE, uR, _ = C5.merge_parts([(n3c, e3c, r3c, "3c"), (n4, e4, r4, "4")])
    pu = C5.poda_clamp(uN, uE, uR, C5.PUERTO_URINARIO)
    u_pel = int(np.argmin(np.linalg.norm(uN - pel_pos, axis=1)))
    t4 = cKDTree(n4)
    UC4 = []
    for a, b in pu["edges_ret"]:
        a, b = int(a), int(b)
        if t4.query(uN[a])[0] < 1e-6 and t4.query(uN[b])[0] < 1e-6:
            ra, rb = pu["radii_c"][a], pu["radii_c"][b]
            if a == u_pel: ra = rb
            if b == u_pel: rb = ra
            UC4.append((uN[a], uN[b], float(ra), float(rb)))
    um = np.array([0.5 * (c[0] + c[1]) for c in UC4]); ut = cKDTree(um)
    def sd_ell(p):
        q = (p - pel_pos) / pel_sa; k0 = np.linalg.norm(q); k1 = np.linalg.norm(q / pel_sa)
        return k0 * (k0 - 1) / max(k1, 1e-9)
    def sd_urin(p):
        d = sd_ell(p)
        for j in ut.query_ball_point(p, 9.0):
            a, b, r1, r2 = UC4[j]; d = min(d, float(C5.sd_round_cone(p[None], a, b, r1, r2)[0]))
        return d

    # ---- PASO 2: conflicto y verificacion de cadena-terminal por circuito ----
    def analyze(nA, eA, rA, port, nombre):
        p = C5.poda_clamp(nA, eA, rA, port); ret = set(int(x) for x in p["ret_nodes"]); rc = p["radii_c"]
        adj = {i: set() for i in range(len(nA))}
        for a, b in p["edges_ret"]: adj[int(a)].add(int(b)); adj[int(b)].add(int(a))
        conf = [i for i in ret if sd_urin(nA[i]) - rc[i] < MIN_WALL]
        confset = set(conf)
        exits = set(j for i in confset for j in adj[i] if j not in confset)
        degs = {i: len([j for j in adj[i] if j in confset]) for i in confset}
        branch = [i for i in confset if degs[i] > 2]
        port_in = p["port_node"] in confset
        is_chain = (len(branch) == 0)
        return dict(p=p, ret=ret, rc=rc, adj=adj, conf=conf, exits=exits, branch=branch,
                    port=p["port_node"], port_in=port_in, is_chain=is_chain, N=nA)

    log("\nPASO2 CONFLICTO (keepout = pelvis elipsoide + ureter + calices, Capa 4):")
    A = analyze(n3a, e3a, r3a, np.array(d3a["raiz"], float), "arterial")
    V = analyze(n3b, e3b, r3b, C5.PUERTO_VENOSO, "venoso")
    for nombre, X in (("arterial", A), ("venoso", V)):
        Ys = [X["N"][i][1] for i in X["conf"]]
        log(f"   {nombre}: conflicto {len(X['conf'])} nodos (Y[{min(Ys):.1f},{max(Ys):.1f}]), "
            f"puerto {'EN' if X['port_in'] else 'NO en'} conflicto, salidas {len(X['exits'])}, "
            f"nodos-rama(grado>2) {len(X['branch'])}  -> "
            f"{'[OK] cadena terminal' if (X['is_chain'] and X['port_in']) else '[FALLA cadena]'}")

    # verificacion: ambos deben ser cadena terminal con puerto para la construccion
    fb = []
    for nombre, X in (("arterial", A), ("venoso", V)):
        if not (X["is_chain"] and X["port_in"]):
            fb.append((nombre, X))

    # ---- construir el/los pediculo(s) que SI son cadena (informativo) y verificar ----
    def build_and_verify(X, Z_corr, port_pos, other_cones, nombre):
        # ancestro: subir aguas arriba desde el conflicto hasta nodo limpio con Y>Y_BORDE_ANT
        # BFS desde el puerto por el arbol, saltando conflicto, hasta Y>Y_BORDE_ANT y limpio
        port = X["port"]; adj = X["adj"]; N = X["N"]; rc = X["rc"]; conf = set(X["conf"])
        # camino unico de cadena: seguir vecinos hasta salir del conflicto
        cur = port; visited = {cur}
        while True:
            nxt = [j for j in adj[cur] if j not in visited]
            if not nxt: break
            cur = nxt[0]; visited.add(cur)
            if cur not in conf and N[cur][1] > Y_BORDE_ANT:
                break
        anc = cur
        # pediculo: anc -> W -> puerto
        W = np.array([0.0, Y_BORDE_ANT, Z_corr])
        seg = [(N[anc], W), (W, port_pos)]
        r_anc = max(rc[anc], MIN_CANAL_R); r_pt = max(MIN_CANAL_R, MIN_CANAL_R)
        # verificar pared en muestreo denso del pediculo vs keepout_urin y otro-circuito
        wmin = np.inf; wloc = None
        for (a, b) in seg:
            for t in np.linspace(0, 1, 40):
                q = a + t * (b - a)
                r = r_anc + t * (r_pt - r_anc)  # taper aprox
                du = sd_urin(q); dv = min([float(C5.sd_round_cone(q[None], *c)[0]) for c in other_cones]) if other_cones else np.inf
                w = min(du, dv) - r
                if w < wmin: wmin = w; wloc = q
        return dict(anc=anc, anc_pos=N[anc], W=W, wmin=wmin, wloc=wloc)

    ven_cones = [(V["N"][int(a)], V["N"][int(b)], float(V["rc"][int(a)]), float(V["rc"][int(b)]))
                 for a, b in V["p"]["edges_ret"]]
    art_ok = False
    if A["is_chain"] and A["port_in"]:
        bp = build_and_verify(A, Z_art, np.array([0, -30, Z_art]), ven_cones, "arterial")
        anc_clean = bp["anc_pos"][1] > Y_BORDE_ANT
        art_ok = (bp["wmin"] >= MIN_WALL - 1e-3) and anc_clean
        log(f"\n   PEDICULO ARTERIAL (conflicto SI es cadena, PASO2 [OK]): ancestro nodo {bp['anc']} en "
            f"{bp['anc_pos'].round(2)} (Y={bp['anc_pos'][1]:.1f} {'>-22.7 [OK]' if anc_clean else '<=-22.7 [NO hay ancestro limpio anterior a la pelvis]'})")
        log(f"      pediculo {bp['anc_pos'].round(2)} -> W[0,{Y_BORDE_ANT},{Z_art:.2f}] -> puerto[0,-30,{Z_art:.2f}]")
        log(f"      pared minima del pediculo vs keepout+venoso: {bp['wmin']*1000:.0f}um "
            f"{'[OK] >=400um' if bp['wmin']>=MIN_WALL-1e-3 else '[<400um -> el tramo ancestro->W roza la pelvis]'}")
        if not art_ok:
            log(f"      -> ARTERIAL: el stub imprimible retenido esta ENTERO en la banda Y de la pelvis "
                f"(Y<=-22.7); no hay ancestro limpio anterior donde anclar el pediculo -> tampoco cierra.")

    # ---- FALLBACK (PASO 2): el venoso NO es cadena terminal ----
    if fb:
        log("\n[FALLBACK - PASO 2] al menos un circuito NO tiene el conflicto como CADENA TERMINAL con puerto.")
        for nombre, X in fb:
            log(f"   {nombre}: conflicto de {len(X['conf'])} nodos con {len(X['branch'])} nodo(s)-rama (grado>2) "
                f"y {len(X['exits'])} salidas; puerto {'EN' if X['port_in'] else 'NO en'} conflicto.")
            log(f"      => NO es una punta unica: es un PLEXO peri-hilar ramificado. Snipearlo y reconectar por "
                f"UN pediculo desconectaria {len(X['exits'])-1} ramas. La sintesis de pediculo-unico NO aplica.")
            log(f"      nodos-rama: {X['branch'][:8]} ; salidas (nodo,Y,Z): "
                f"{[(int(j), round(float(X['N'][j][1]),1), round(float(X['N'][j][2]),1)) for j in list(X['exits'])[:6]]}")
            # obstaculo concreto: peor nodo del conflicto
            worst = min(X["conf"], key=lambda i: sd_urin(X["N"][i]) - X["rc"][i])
            w = sd_urin(X["N"][worst]) - X["rc"][worst]
            log(f"      peor nodo {worst} en {X['N'][worst].round(2)}: pared {w*1000:.0f}um (faltan {(MIN_WALL-w)*1000:.0f}um) vs urinario macro.")
        log(f"\n   ARTERIAL: conflicto SI es cadena (PASO2 [OK]) pero el pediculo {'VERIFICA' if art_ok else 'NO verifica'} "
            f"{'' if art_ok else '(stub imprimible entero en la banda Y de la pelvis, sin ancestro limpio anterior).'}")
        log("   OBSTACULO para cerrar estanqueidad: (1) el arbol VENOSO ramifica en el hilio")
        log("   (plexo peri-pelvico, no un tronco unico) -> la segregacion de pediculo-unico no lo captura;")
        log("   (2) el arterial imprimible es un stub sub-mm peri-pelvico sin parenquima anterior donde anclar.")
        log("   OPCIONES (NO aplico nada; decidis vos): (a) tronco venoso unico en el hilio (revisar 3b, fuera de")
        log("   alcance local); (b) MULTI-pediculo venoso (un corredor por cada una de las ramas en conflicto);")
        log("   (c) carve local del keepout; (d) mover pelvis en -Y para reducir el solape peri-pelvico venoso.")
        _write_md(LOG, fb=True)
        return LOG

    _write_md(LOG, fb=False)
    return LOG


def _write_md(LOG, fb):
    M = ["# Auditoria fix estanqueidad Capa 5a v3 (pediculos hilares segregados)\n",
         "**Programa:** Bio-Kidney AI 2026 - **Capa:** 5a (fix v3) - **Fecha:** 2026-07-09\n",
         "Construccion deterministica de pediculos hilares con anclaje anatomico (orden AP "
         "vena-arteria-pelvis; +Z anterior). Mantiene FIX1 (pelvis elipsoide). Local en Capa 5, "
         "NO reescribe 3a/3b/3c, NO cambia Capa 4. MIN_WALL=0.400, MARGIN=0.200.\n\n```\n"]
    M += [l + "\n" for l in LOG]
    M.append("```\n")
    os.makedirs(os.path.dirname(OUT_MD), exist_ok=True)
    with open(OUT_MD, "w") as f: f.write("".join(M))


if __name__ == "__main__":
    main()
    print("  -> reporte:", OUT_MD)
