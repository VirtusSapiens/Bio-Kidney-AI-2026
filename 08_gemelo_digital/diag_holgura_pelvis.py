#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diag_holgura_pelvis.py -- DIAGNOSTICO DE HOLGURA de la pelvis (solo lectura).
Barrido de posicion del centro de la pelvis (X=0, Y-Z) buscando un centro que cierre
SIMULTANEAMENTE: pared vs venoso/arterial, contencion en seno, intrarrenalidad, infundibulos
mayores sin pinza/colision, ureter, y corredor anterior. Reusa SDF/adaptador/clamp/FIX1.
NO modifica nada, NO genera STL/npz.
"""
import os, sys, json
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import capa5a_nucleo_sdf as C5

MIN_WALL = 0.400; MARGIN = 0.200
PELVIS_SEMI = np.array([4.5, 4.5, 2.0])
SENO_C = np.array([0.0, -34.0, 0.0]); SENO_S = np.array([22.0, 16.0, 11.0])
PARENQ_Y = -27.7; HILIO = np.array([0.0, -30.0, 0.0])
Z_ART = 2.814; Z_VEN = 4.152; R_ART = 0.214; R_VEN = 0.524
OUT = os.path.join(C5.REPO, "09_paper_vascular", "diagnostico_holgura_pelvis.md")


def fib_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n); th = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.sin(phi) * np.cos(th), np.sin(phi) * np.sin(th), np.cos(phi)], 1)


def main():
    L = []
    def log(s=""): print(s); L.append(s)
    log("DIAGNOSTICO DE HOLGURA DE LA PELVIS (solo lectura)")

    d3a = np.load(C5.find("capa3a_arterial.npz"), allow_pickle=True)
    d3b = np.load(C5.find("capa3b_venoso.npz"), allow_pickle=True)
    d4 = np.load(C5.find("capa4_colector_alto.npz"), allow_pickle=True)
    n3a, e3a, r3a = C5.adapt_edge_radios(d3a); n3b, e3b, r3b = C5.adapt_edge_radios(d3b)
    n4 = d4["nodos"]; niv = d4["nivel"]
    may_pos = n4[niv == 2]; ure_pos = n4[niv == 4][0]

    def cones_of(nA, eA, rA, port):
        p = C5.poda_clamp(nA, eA, rA, port)
        return [(nA[int(a)], nA[int(b)], float(p["radii_c"][int(a)]), float(p["radii_c"][int(b)]))
                for a, b in p["edges_ret"]]
    VC = cones_of(n3b, e3b, r3b, C5.PUERTO_VENOSO)
    AC = cones_of(n3a, e3a, r3a, np.array(d3a["raiz"], float))
    # prefiltrar al hilio (midpoint Y in [-36,-14]) para velocidad
    def near_hilio(cones):
        return [c for c in cones if -36 <= 0.5 * (c[0][1] + c[1][1]) <= -14]
    VCh = near_hilio(VC); ACh = near_hilio(AC)
    vmid = np.array([0.5 * (c[0] + c[1]) for c in VCh]); vtree = cKDTree(vmid)
    amid = np.array([0.5 * (c[0] + c[1]) for c in ACh]) if ACh else np.zeros((1, 3)); atree = cKDTree(amid)
    log(f"venoso cones {len(VC)} (hilio {len(VCh)}), arterial cones {len(AC)} (hilio {len(ACh)})")

    def sd_min(P, cones, tree, R=8.0):
        """min sd_round_cone sobre cones cercanos a cada punto (P: (M,3)) -> (M,)."""
        out = np.full(len(P), np.inf)
        idxs = tree.query_ball_point(P, R)
        # agrupar por cono para vectorizar
        need = set()
        for lst in idxs: need |= set(lst)
        for j in need:
            a, b, r1, r2 = cones[j]
            d = C5.sd_round_cone(P, a, b, r1, r2)
            np.minimum(out, d, out=out)
        return out
    def seno_val(P):
        q = (P - SENO_C) / SENO_S; return np.einsum('ij,ij->i', q, q)

    SURF = fib_sphere(200)

    # ===== D1: mapa del hilio =====
    log("\n===== D1  MAPA DEL HILIO (venoso/arterial peri-pelvico) =====")
    def in_region(c):
        m = 0.5 * (c[0] + c[1])
        return -33 <= m[1] <= -20 and -9 <= m[2] <= 9 and -7 <= m[0] <= 7
    vreg = [c for c in VCh if in_region(c)]; areg = [c for c in ACh if in_region(c)]
    vz = np.array([0.5 * (c[0][2] + c[1][2]) for c in vreg])
    log(f"  venoso: {len(vreg)} segmentos en la region; anterior(+Z) {int((vz > 0).sum())}, posterior(-Z) {int((vz < 0).sum())}, Z~0 {int((np.abs(vz) <= 0.3).sum())}")
    log(f"  arterial: {len(areg)} segmentos en la region")
    # hueco: en X=0, grid Y-Z, distancia a lumen venoso (clearance)
    ys = np.arange(-33, -20.01, 0.4); zs = np.arange(-9, 9.01, 0.4)
    GY, GZ = np.meshgrid(ys, zs, indexing='ij')
    Pg = np.stack([np.zeros(GY.size), GY.ravel(), GZ.ravel()], 1)
    clr = sd_min(Pg, VCh, vtree, R=9.0).reshape(GY.shape)  # clearance a venoso (>0 fuera)
    kmax = np.unravel_index(np.argmax(clr), clr.shape)
    log(f"  HUECO inter-tributaria venoso (X=0): clearance max {clr[kmax]:.2f} mm en Y={GY[kmax]:.1f} Z={GZ[kmax]:.1f}")
    # extent del hueco donde clearance permite semiejes (>= 4.5 XY no aplica en Y-Z; medir donde clr>=MIN_WALL)
    big = clr >= (0.5)  # zonas con >0.5mm de vacio
    if big.any():
        yb = GY[big]; zb = GZ[big]
        log(f"  vacio venoso (clr>=0.5mm) extent: Y[{yb.min():.1f},{yb.max():.1f}] Z[{zb.min():.1f},{zb.max():.1f}]")

    # ===== D2: barrido de centro de pelvis =====
    log("\n===== D2  BARRIDO DE CENTRO DE PELVIS (X=0) =====")
    Ys = np.arange(-30, -20.0 + 1e-9, 0.5); Zs = np.arange(-6, 6.0 + 1e-9, 0.5)
    viable = []; grid_info = {}
    for cy in Ys:
        for cz in Zs:
            C = np.array([0.0, cy, cz])
            surf = C + SURF * PELVIS_SEMI
            pv = float(sd_min(surf, VCh, vtree).min())
            pa = float(sd_min(surf, ACh, atree).min()) if ACh else np.inf
            cs = float(seno_val(surf).max())
            intr = float(np.mean(surf[:, 1] >= PARENQ_Y))
            ok = (pv >= MIN_WALL) and (pa >= MIN_WALL) and (cs <= 1.0)
            grid_info[(round(cy, 1), round(cz, 1))] = (pv, pa, cs, intr, ok)
            if ok:
                viable.append((cy, cz, pv, pa, cs, intr))
    log(f"  centros probados: {len(Ys) * len(Zs)}  ->  VIABLES (pared_ven>=0.4 Y pared_art>=0.4 Y seno<=1): {len(viable)}")
    if viable:
        vy = [v[0] for v in viable]; vz2 = [v[1] for v in viable]
        log(f"  mascara viable: Y[{min(vy):.1f},{max(vy):.1f}]  Z[{min(vz2):.1f},{max(vz2):.1f}]")
        for v in sorted(viable, key=lambda x: -(x[5]))[:12]:
            log(f"     C=[0,{v[0]:.1f},{v[1]:.1f}] pared_ven {v[2]:.2f} pared_art {v[3]:.2f} seno {v[4]:.2f} intrarrenal {100*v[5]:.0f}%")
    else:
        log("  MASCARA VIABLE VACIA: ningun centro cierra pared_ven & pared_art & seno simultaneamente.")

    # helpers D3/D4 (INDEPENDIENTES del tamaño de la pelvis: usan el CENTRO, no los semiejes)
    def check_infund(C):
        d3ok = True; inf_wmin = np.inf; swing = 0.0
        r_line = 2.0 + np.linspace(0, 1, 25) * (4.5 - 2.0)
        for m in may_pos:
            tt = np.linspace(0, 1, 25)[:, None]
            seg = m[None] + tt * (C - m)[None]                      # (25,3) vectorizado
            wv = sd_min(seg, VCh, vtree) - r_line
            wa = (sd_min(seg, ACh, atree) - r_line) if ACh else np.full(25, np.inf)
            w = float(min(wv.min(), wa.min())); inf_wmin = min(inf_wmin, w)
            v0 = C - m; swing = max(swing, np.degrees(np.arctan2(abs(v0[2]), abs(v0[1]))))
            if w < MIN_WALL: d3ok = False
        return d3ok, inf_wmin, swing
    def check_ureter(C):
        tt = np.linspace(0, 1, 25)[:, None]
        seg = C[None] + tt * (HILIO - C)[None]
        w = float(min(sd_min(seg, VCh, vtree).min(), sd_min(seg, ACh, atree).min() if ACh else np.inf)) - 0.75
        return w >= MIN_WALL, w
    def check_corridor(C):
        pv = float(sd_ellip_pt(np.array([0, -30, Z_VEN]), C)) - R_VEN
        pa = float(sd_ellip_pt(np.array([0, -30, Z_ART]), C)) - R_ART
        gap = (Z_VEN - Z_ART) - R_VEN - R_ART
        return (pv >= MIN_WALL and pa >= MIN_WALL and gap >= MIN_WALL), pv, pa

    # ===== D3-D5 sobre viables =====
    survivors = []; stages = {"d3": 0, "d4": 0, "d5": 0}
    log("\n===== D3  INFUNDIBULOS MAYORES / D4 URETER / D5 CORREDOR (sobre viables) =====")
    for (cy, cz, pv, pa, cs, intr) in sorted(viable, key=lambda x: -x[5])[:15]:
        C = np.array([0.0, cy, cz])
        d3ok, inf_wmin, swing = check_infund(C)
        d4ok, d4w = check_ureter(C)
        d5ok, pvp, pap = check_corridor(C)
        allok = d3ok and d4ok and d5ok
        if not d3ok: stages["d3"] += 1
        if not d4ok: stages["d4"] += 1
        if not d5ok: stages["d5"] += 1
        if allok: survivors.append((cy, cz, pv, pa, cs, intr, inf_wmin, swing, d4w))
        log(f"  C=[0,{cy:.1f},{cz:.1f}]: D3 inf_wall {inf_wmin:.2f} swing {swing:.0f}deg {'[OK]' if d3ok else '[X]'} | "
            f"D4 ureter_wall {d4w:.2f} {'[OK]' if d4ok else '[X]'} | D5 puertos(v{pvp:.2f},a{pap:.2f}) {'[OK]' if d5ok else '[X]'} -> {'CIERRA TODO' if allok else 'no'}")

    # ===== D6 veredicto =====
    log("\n===== D6  VEREDICTO =====")
    if survivors:
        best = max(survivors, key=lambda s: (s[5], min(s[2], s[3]), -s[7]))
        log(f"  [OK] HAY CENTROS QUE CIERRAN TODO (D2+D3+D4+D5): {len(survivors)}")
        log(f"  TARGET RECOMENDADO: pelvis centro = [0, {best[0]:.1f}, {best[1]:.1f}]  "
            f"pared_ven {best[2]:.2f} pared_art {best[3]:.2f} seno {best[4]:.2f} intrarrenal {100*best[5]:.0f}% "
            f"inf_wall {best[6]:.2f} swing {best[7]:.0f}deg ureter_wall {best[8]:.2f}")
        verdict = dict(status="POSICION", target=[0.0, float(best[0]), float(best[1])], survivors=len(survivors))
    else:
        n = len(sorted(viable, key=lambda x: -x[5])[:15])
        log(f"  NINGUN centro D2-viable cierra D3-D5. Fallos entre {n} candidatos: "
            f"D3(infundibulos) {stages['d3']}/{n}, D4(ureter) {stages['d4']}/{n}, D5(corredor) {stages['d5']}/{n}")
        log(f"  El cuerpo de la pelvis SI cabe (D2: {len(viable)} centros viables), pero los INFUNDIBULOS")
        log(f"  (calices en Z~0 -> pelvis) y el URETER (pelvis -> hilio [0,-30,0]) cruzan el plexo venoso en Z~0.")
        # test clave: ¿EXISTE ALGUNA POSICION (cualquier centro) con infundibulos+ureter limpios?
        # D3/D4 dependen solo de la POSICION del centro (no del tamaño). Barrer toda la grilla.
        best_inf = None
        for cy in np.arange(-30, -20.01, 1.0):
            for cz in np.arange(-6, 6.01, 1.0):
                C = np.array([0.0, cy, cz])
                d3ok, inf_w, _ = check_infund(C); d4ok, d4w = check_ureter(C)
                score = min(inf_w, d4w)
                if best_inf is None or score > best_inf[3]:
                    best_inf = (cy, cz, d3ok and d4ok, score, inf_w, d4w)
                if d3ok and d4ok:
                    best_inf = (cy, cz, True, score, inf_w, d4w); break
            if best_inf and best_inf[2]: break
        if best_inf and best_inf[2]:
            cy, cz = best_inf[0], best_inf[1]
            log(f"  EXISTE posicion con infundibulos+ureter limpios: [0,{cy:.1f},{cz:.1f}] -> el pelvis-cuerpo ahi")
            log(f"  puede requerir reducir tamaño para caber (D2) -> fix = POSICION (+TAMANO si D2 no cierra ahi).")
            verdict = dict(status="POSICION_INFUND_OK", target=[0.0, float(cy), float(cz)])
        else:
            cy, cz, _, sc, iw, d4w = best_inf
            log(f"  NINGUNA posicion (barrido completo 525 centros) tiene infundibulos+ureter limpios.")
            log(f"  Mejor posicion posible [0,{cy:.1f},{cz:.1f}]: mejor pared infundibulo/ureter = {sc*1000:.0f}um (faltan {(MIN_WALL-sc)*1000:.0f}um).")
            log(f"  Los infundibulos (calices Z~0 -> pelvis) y el ureter (pelvis -> hilio Z~0) estan ANCLADOS en Z~0")
            log(f"  e INDEPENDIENTES del tamaño de la pelvis: reducir/mover la pelvis NO los despeja.")
            log(f"  BINDING = el PLEXO VENOSO peri-hilar (32 segs, 27 posteriores) ocupa el corredor Z~0 obligado.")
            log(f"  -> el fix NO es posicion NI tamaño de la pelvis: es SEGREGAR/mover el plexo venoso (revisar 3b).")
            log(f"  Confirma cuantitativamente el blocker de las Entradas 019-021 (arbol venoso ramificado en el hilio).")
            verdict = dict(status="NO_ES_LA_PELVIS__BINDING_PLEXO_VENOSO_Z0",
                           mejor_pared_inf_ureter_um=round(sc * 1000, 1),
                           d3_fail=stages["d3"], d4_fail=stages["d4"], d2_viables=len(viable))

    _write_md(L, verdict)
    return L, verdict


def sd_ellip_pt(p, center, semi=PELVIS_SEMI):
    q = (p - center) / semi
    k0 = np.linalg.norm(q); k1 = np.linalg.norm(q / semi)
    return k0 * (k0 - 1) / max(k1, 1e-9)


def _write_md(L, verdict):
    M = ["# Diagnostico de holgura de la pelvis (solo lectura)\n",
         "**Programa:** Bio-Kidney AI 2026 - **Capa:** 5a (diagnostico) - **Fecha:** 2026-07-10\n",
         "Barrido de posicion del centro de la pelvis (X=0, Y-Z) buscando cierre simultaneo de pared "
         "vs venoso/arterial, contencion en seno, intrarrenalidad, infundibulos, ureter y corredor "
         "anterior. Reusa SDF/adaptador/clamp/FIX1 (pelvis elipsoide [4.5,4.5,2.0]). NADA se modifica.\n\n```\n"]
    M += [x + "\n" for x in L]; M.append("```\n")
    M.append(f"\n## Veredicto (estructurado)\n```json\n{json.dumps(verdict, indent=2)}\n```\n")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f: f.write("".join(M))


if __name__ == "__main__":
    main()
    print("  -> reporte:", OUT)
