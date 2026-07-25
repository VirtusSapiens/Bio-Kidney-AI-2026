#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diag_estanqueidad_5a.py -- DIAGNOSTICO (solo lectura) de la falla de estanqueidad 5a.
Reusa la normalizacion/adaptadores/poda de capa5a_nucleo_sdf (mismo clamp). NO re-malla,
NO modifica 5a ni genera 5b. Clasifica los 1118 voxeles urin<->vascular compartidos.
"""
import os, sys, json
import numpy as np
from scipy.spatial import cKDTree

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _here)
import capa5a_nucleo_sdf as C5   # REUSO: adaptadores, merge, poda_clamp, eval_sdf, cones, SDF

VX = C5.VOXEL_5A
MINR = C5.MIN_CANAL_R
REPO = C5.REPO
LEVEL_NAME = {0: "papila_junction", 1: "caliz_menor", 2: "caliz_mayor", 3: "pelvis", 4: "ureter"}
OUT = os.path.join(REPO, "09_paper_vascular", "diagnostico_estanqueidad_5a.md")


def seg_seg(p1, q1, p2, q2):
    d1 = q1 - p1; d2 = q2 - p2; r = p1 - p2
    a = d1 @ d1; e = d2 @ d2; f = d2 @ r; EPS = 1e-12
    if a <= EPS and e <= EPS:
        return np.linalg.norm(p1 - p2), 0.0, 0.0
    if a <= EPS:
        s = 0.0; t = np.clip(f / e, 0, 1)
    else:
        c = d1 @ r
        if e <= EPS:
            t = 0.0; s = np.clip(-c / a, 0, 1)
        else:
            b = d1 @ d2; den = a * e - b * b
            s = np.clip((b * f - c * e) / den, 0, 1) if den > EPS else 0.0
            t = (b * s + f) / e
            if t < 0: t = 0.0; s = np.clip(-c / a, 0, 1)
            elif t > 1: t = 1.0; s = np.clip((b - c) / a, 0, 1)
    c1 = p1 + d1 * s; c2 = p2 + d2 * t
    return float(np.linalg.norm(c1 - c2)), float(s), float(t)


def argmin_cone(cones, P):
    best = np.full(len(P), np.inf); arg = np.full(len(P), -1, int)
    for k, (a, b, r1, r2) in enumerate(cones):
        d = C5.sd_round_cone(P, a, b, r1, r2)
        m = d < best; best[m] = d[m]; arg[m] = k
    return best, arg


def main():
    L = []
    def log(s=""): print(s); L.append(s)
    log("DIAGNOSTICO ESTANQUEIDAD 5a")

    # ---- reconstruir circuitos EXACTOS como 5a ----
    d3a = np.load(C5.find("capa3a_arterial.npz"), allow_pickle=True)
    d3b = np.load(C5.find("capa3b_venoso.npz"), allow_pickle=True)
    d3c = np.load(C5.find("capa3c_colector.npz"), allow_pickle=True)
    d4 = np.load(C5.find("capa4_colector_alto.npz"), allow_pickle=True)
    n3a, e3a, r3a = C5.adapt_edge_radios(d3a)
    n3b, e3b, r3b = C5.adapt_edge_radios(d3b)
    n3c, e3c, r3c = C5.adapt_3c(d3c)
    n4, e4, r4 = C5.adapt_4(d4)
    n3ab = np.load(C5.find("capa3ab_peritubular.npz"), allow_pickle=True)["puntos_drenaje"].astype(float)

    # URINARIO
    uN, uE, uR, _ = C5.merge_parts([(n3c, e3c, r3c, "3c"), (n4, e4, r4, "4")])
    pu = C5.poda_clamp(uN, uE, uR, C5.PUERTO_URINARIO)
    # etiqueta por nodo urinario
    tree4 = cKDTree(n4); dd4, jj4 = tree4.query(uN, k=1)
    ulabel = np.array(["colector_3c"] * len(uN), dtype=object)
    m4 = dd4 < 1e-3
    for i in np.where(m4)[0]:
        ulabel[i] = LEVEL_NAME[int(d4["nivel"][jj4[i]])]
    u_clamp_node = uR < MINR                        # nativo < piso -> clampeado
    u_cones = C5.build_cones(uN, pu["edges_ret"], pu["radii_c"])
    u_edges = pu["edges_ret"]

    # VASCULAR por componente
    vN, vE, vR, _ = C5.merge_parts([(n3a, e3a, r3a, "3a"),
                                    (n3ab, np.empty((0, 2), int), np.zeros(len(n3ab)), "3ab"),
                                    (n3b, e3b, r3b, "3b")])
    ncv, labv = C5.comps(len(vN), vE)
    lab_art = labv[int(np.argmin(np.linalg.norm(vN - np.array(d3a["raiz"], float), axis=1)))]
    lab_ven = labv[int(np.argmin(np.linalg.norm(vN - C5.PUERTO_VENOSO, axis=1)))]
    vasc = {}
    for lab, nombre, port in [(lab_art, "arterial", np.array(d3a["raiz"], float)),
                              (lab_ven, "venoso", C5.PUERTO_VENOSO)]:
        mask = labv == lab; idx = np.where(mask)[0]
        remap = -np.ones(len(vN), int); remap[idx] = np.arange(len(idx))
        em = mask[vE[:, 0]] & mask[vE[:, 1]]
        sube = remap[vE[em]]; subn = vN[idx]; subr = vR[idx]
        p = C5.poda_clamp(subn, sube, subr, port)
        vasc[nombre] = dict(N=subn, edges=p["edges_ret"], rc=p["radii_c"], rnat=subr,
                            cones=C5.build_cones(subn, p["edges_ret"], p["radii_c"]))
    v_all_cones = vasc["arterial"]["cones"] + vasc["venoso"]["cones"]
    # etiqueta vascular por cono: (nombre, edge_local)
    v_cone_meta = [("arterial", e) for e in vasc["arterial"]["edges"]] + \
                  [("venoso", e) for e in vasc["venoso"]["edges"]]

    # ---- (1) grid comun identico a 5a -> voxeles compartidos ----
    allpts = np.vstack([np.vstack([c[0], c[1]]) for c in u_cones + v_all_cones])
    rmax = max(max(c[2], c[3]) for c in u_cones + v_all_cones)
    lo = allpts.min(0) - (rmax + C5.BBOX_MARGIN); hi = allpts.max(0) + (rmax + C5.BBOX_MARGIN)
    shp = tuple(int(np.ceil((hi[i] - lo[i]) / VX)) + 1 for i in range(3))
    su = C5.eval_sdf(u_cones, lo, shp, VX)
    sv = C5.eval_sdf(v_all_cones, lo, shp, VX)
    both = (su < 0) & (sv < 0)
    idx = np.argwhere(both)
    P = lo + idx * VX
    log(f"(1) voxeles compartidos: {len(P)} (5a reporto 1118)")

    # ---- (2) segmento responsable urin/vasc por voxel ----
    _, uarg = argmin_cone(u_cones, P)
    _, varg = argmin_cone(v_all_cones, P)

    def urin_class(vi, edge):
        a, b = int(edge[0]), int(edge[1])
        na = ulabel[a] if np.linalg.norm(P[vi] - uN[a]) < np.linalg.norm(P[vi] - uN[b]) else ulabel[b]
        cl = bool(u_clamp_node[a] or u_clamp_node[b])
        return na, cl

    def vasc_class(vi, k):
        nombre, edge = v_cone_meta[k]
        a, b = int(edge[0]), int(edge[1])
        rnat = vasc[nombre]["rnat"]
        cl = bool(rnat[a] < MINR or rnat[b] < MINR)
        return nombre, cl

    # ---- (3) tabla nivel_urinario x origen_vascular ----
    from collections import defaultdict
    tab = defaultdict(lambda: [0, 0])   # key (ulevel, vorigin) -> [both_native, any_clamped]
    per_side_clamp = dict(u_clamp=0, v_clamp=0, both_native=0)
    for vi in range(len(P)):
        ul, ucl = urin_class(vi, u_edges[uarg[vi]])
        vo, vcl = vasc_class(vi, varg[vi])
        any_clamp = ucl or vcl
        tab[(ul, vo)][1 if any_clamp else 0] += 1
        if ucl: per_side_clamp["u_clamp"] += 1
        if vcl: per_side_clamp["v_clamp"] += 1
        if not any_clamp: per_side_clamp["both_native"] += 1
    log("(3) tabla (nivel_urinario x vascular): [both_native | any_clamped]")
    for (ul, vo), (bn, ac) in sorted(tab.items(), key=lambda x: -(x[1][0] + x[1][1])):
        log(f"    {ul:16s} x {vo:8s}: both_native={bn:5d}  any_clamped={ac:5d}")
    log(f"    -> voxeles con >=1 lado clampeado: {sum(v[1] for v in tab.values())} ; ambos nativos: {sum(v[0] for v in tab.values())}")
    log(f"       (urin clampeado en {per_side_clamp['u_clamp']}, vasc clampeado en {per_side_clamp['v_clamp']})")

    # ---- (4) pared inter-circuito: NATIVA (revela separacion real) por voxel ----
    walls_nat = np.empty(len(P)); walls_clamp = np.empty(len(P))
    for vi in range(len(P)):
        ue = u_edges[uarg[vi]]; ua, ub = int(ue[0]), int(ue[1])
        vname, ve = v_cone_meta[varg[vi]]; va, vb = int(ve[0]), int(ve[1])
        d, s, t = seg_seg(uN[ua], uN[ub], vasc[vname]["N"][va], vasc[vname]["N"][vb])
        # nativo
        run = uR[ua] + s * (uR[ub] - uR[ua]); rvn = vasc[vname]["rnat"][va] + t * (vasc[vname]["rnat"][vb] - vasc[vname]["rnat"][va])
        walls_nat[vi] = d - run - rvn
        # clampeado
        ruc = pu["radii_c"][ua] + s * (pu["radii_c"][ub] - pu["radii_c"][ua])
        rvc = vasc[vname]["rc"][va] + t * (vasc[vname]["rc"][vb] - vasc[vname]["rc"][va])
        walls_clamp[vi] = d - ruc - rvc
    b_neg = int((walls_nat < 0).sum())
    b_100 = int(((walls_nat >= 0) & (walls_nat < 0.100)).sum())
    b_100_400 = int(((walls_nat >= 0.100) & (walls_nat < 0.400)).sum())
    b_400 = int((walls_nat >= 0.400).sum())
    log(f"(4) PARED NATIVA inter-circuito por voxel compartido (histograma):")
    log(f"    <0 (solape real, radios nativos): {b_neg}")
    log(f"    0-100um : {b_100}")
    log(f"    100-400um : {b_100_400}")
    log(f"    >=400um (pared real imprimible-separable -> SUB-MUESTREO/clamp): {b_400}")
    # min superficie-superficie (clampeado, geometria actual de 5a) con localizacion
    imin = int(np.argmin(walls_clamp))
    ue = u_edges[uarg[imin]]; vname, ve = v_cone_meta[varg[imin]]
    log(f"    min pared CLAMPEADA (geometria 5a): {walls_clamp[imin]*1000:.1f}um  en {P[imin].round(2)}")
    log(f"       urin edge {tuple(int(x) for x in ue)} ({ulabel[int(ue[0])]}/{ulabel[int(ue[1])]}) vs {vname} edge {tuple(int(x) for x in ve)}")
    log(f"    min pared NATIVA: {walls_nat.min()*1000:.1f}um")

    # ---- (5) zoom hilio: pelvis/ureter vs tronco venoso en Y[-31,-27] ----
    # cones urinarios cuyo nivel es pelvis o ureter
    pel_ure = []
    for k, e in enumerate(u_edges):
        a, b = int(e[0]), int(e[1])
        if ulabel[a] in ("pelvis", "ureter") or ulabel[b] in ("pelvis", "ureter"):
            pel_ure.append(u_cones[k])
    ven_hilio = []
    vc = vasc["venoso"]
    for k, e in enumerate(vc["edges"]):
        a, b = int(e[0]), int(e[1])
        y = 0.5 * (vc["N"][a][1] + vc["N"][b][1])
        if -31 <= y <= -27:
            ven_hilio.append(vc["cones"][k])
    hil_min = np.inf; hil_loc = None
    for (a1, b1, r1a, r1b) in pel_ure:
        for (a2, b2, r2a, r2b) in ven_hilio:
            d, s, t = seg_seg(a1, b1, a2, b2)
            wall = d - (r1a + s * (r1b - r1a)) - (r2a + t * (r2b - r2a))
            if wall < hil_min:
                hil_min = wall; hil_loc = (0.5 * (a1 + b1), 0.5 * (a2 + b2))
    log(f"(5) ZOOM HILIO (pelvis/ureter vs venoso Y[-31,-27]): pares {len(pel_ure)}x{len(ven_hilio)}")
    if hil_loc is not None:
        log(f"    pared lumen-lumen minima: {hil_min*1000:.1f}um  "
            f"{'[<400um: componente hilio]' if hil_min < 0.400 else '[>=400um: NO hay colision de hilio]'}")
        log(f"    urin~{hil_loc[0].round(2)}  venoso~{hil_loc[1].round(2)}")
    else:
        log("    sin segmentos venosos en la banda de hilio")

    # ---- (6) genero: localizar handle = fusiones urin no-adyacentes ----
    mids = np.array([0.5 * (uN[int(e[0])] + uN[int(e[1])]) for e in u_edges])
    tm = cKDTree(mids)
    fus = []
    for i in range(len(u_edges)):
        for j in tm.query_ball_point(mids[i], 4.0):
            if j <= i: continue
            if len(set(u_edges[i].tolist()) & set(u_edges[j].tolist())): continue  # adyacentes
            a1, b1 = uN[int(u_edges[i][0])], uN[int(u_edges[i][1])]
            a2, b2 = uN[int(u_edges[j][0])], uN[int(u_edges[j][1])]
            d, s, t = seg_seg(a1, b1, a2, b2)
            ru = pu["radii_c"][int(u_edges[i][0])] + s * (pu["radii_c"][int(u_edges[i][1])] - pu["radii_c"][int(u_edges[i][0])])
            rv = pu["radii_c"][int(u_edges[j][0])] + t * (pu["radii_c"][int(u_edges[j][1])] - pu["radii_c"][int(u_edges[j][0])])
            if d < ru + rv:   # lumenes fusionan -> posible handle
                fus.append((i, j, d, ru + rv, 0.5 * (a1 + b1)))
    log(f"(6) GENERO/handle: fusiones urin no-adyacentes (lumen solapa, clampeado): {len(fus)}")
    for (i, j, d, s, m) in fus[:8]:
        a, b = int(u_edges[i][0]), int(u_edges[i][1]); c, e = int(u_edges[j][0]), int(u_edges[j][1])
        log(f"    edge{i}({ulabel[a]}/{ulabel[b]}) x edge{j}({ulabel[c]}/{ulabel[e]}): d={d:.3f} solape={s-d:.3f}mm en {m.round(2)}")

    # ---- (7) numeros cortados de 5a ----
    meta = np.load(C5.find("capa5a_meta.npz"), allow_pickle=True)
    vol_lumen = json.loads(str(meta["vol_lumen"]))
    log(f"(7) NUMEROS CORTADOS de 5a:")
    log(f"    porosidad: {float(meta['porosidad_pct']):.4f} %")
    log(f"    lumen por circuito (mm3): " + ", ".join(f"{k}={v:.2f}" for k, v in vol_lumen.items()))
    log(f"    gap arterial<->venoso (nodos): {float(meta['gap_av_mm']):.3f} mm")

    # ---- VEREDICTO ----
    frac_sub = 100.0 * b_400 / len(P)
    frac_col = 100.0 * (b_neg + b_100 + b_100_400) / len(P)
    log("\nVEREDICTO:")
    log(f"    SUB-MUESTREO (pared nativa >=400um): {b_400}/{len(P)} ({frac_sub:.1f}%)")
    log(f"    COLISION REAL (pared nativa <400um): {b_neg+b_100+b_100_400}/{len(P)} ({frac_col:.1f}%)")
    log(f"    componente hilio: {'SI' if hil_min < 0.400 else 'NO'} (pared pelvis/ureter<->venoso {hil_min*1000:.1f}um)")

    return dict(L=L, nP=len(P), tab=dict(tab), per=per_side_clamp,
                b_neg=b_neg, b_100=b_100, b_100_400=b_100_400, b_400=b_400,
                walls_nat=walls_nat, walls_clamp=walls_clamp, imin=imin, P=P,
                uarg=uarg, varg=varg, u_edges=u_edges, ulabel=ulabel, uN=uN,
                v_cone_meta=v_cone_meta, vasc=vasc, hil_min=hil_min, hil_loc=hil_loc,
                pel_ure=len(pel_ure), ven_hilio=len(ven_hilio), fus=fus,
                meta=meta, vol_lumen=vol_lumen, frac_sub=frac_sub, frac_col=frac_col)


if __name__ == "__main__":
    R = main()
    # ---- escribir markdown ----
    b_neg, b_100, b_100_400, b_400, nP = R["b_neg"], R["b_100"], R["b_100_400"], R["b_400"], R["nP"]
    M = []
    M.append("# Diagnostico estanqueidad Capa 5a (solo lectura)\n")
    M.append("**Programa:** Bio-Kidney AI 2026 - **Capa:** 5a (diagnostico) - **Fecha:** 2026-07-09\n")
    M.append("Localiza y clasifica los voxeles compartidos urinario<->vascular de la auditoria 5a. "
             "Reusa centerlines/radios ya normalizados (mismo adaptador y clamp de 5a). NO re-malla, "
             "NO modifica 5a, NO genera 5b.\n")
    M.append(f"## (1) Voxeles compartidos\n{nP} voxeles urinario & vascular (5a reporto 1118).\n")
    M.append("## (2)-(3) Clasificacion nivel_urinario x origen_vascular\n")
    M.append("| nivel urinario | vascular | both_native | any_clamped |\n|---|---|---|---|\n")
    for (ul, vo), (bn, ac) in sorted(R["tab"].items(), key=lambda x: -(x[1][0] + x[1][1])):
        M.append(f"| {ul} | {vo} | {bn} | {ac} |\n")
    tot_clamp = sum(v[1] for v in R["tab"].values()); tot_nat = sum(v[0] for v in R["tab"].values())
    M.append(f"\n- Voxeles con **>=1 lado clampeado**: **{tot_clamp}** ; **ambos nativos**: **{tot_nat}**.\n")
    M.append(f"- (urin clampeado en {R['per']['u_clamp']}, vascular clampeado en {R['per']['v_clamp']} voxeles)\n")
    M.append("## (4) Pared inter-circuito NATIVA (radios originales) por voxel\n")
    M.append("Distingue sub-muestreo (pared real >=400um, solo el voxel 200um la puentea) de colision real (<400um).\n\n")
    M.append("| pared nativa | voxeles |\n|---|---|\n")
    M.append(f"| <0 (solape real) | {b_neg} |\n| 0-100um | {b_100} |\n| 100-400um | {b_100_400} |\n| >=400um (sub-muestreo) | {b_400} |\n")
    imin = R["imin"]; P = R["P"]; ue = R["u_edges"][R["uarg"][imin]]; vname, ve = R["v_cone_meta"][R["varg"][imin]]
    M.append(f"\n- Min pared CLAMPEADA (geometria actual de 5a): **{R['walls_clamp'][imin]*1000:.1f} um** en {P[imin].round(2)} "
             f"(urin {R['ulabel'][int(ue[0])]}/{R['ulabel'][int(ue[1])]} vs {vname}).\n")
    M.append(f"- Min pared NATIVA: **{R['walls_nat'].min()*1000:.1f} um**.\n")
    M.append("## (5) Zoom hilio (pelvis/ureter vs tronco venoso, Y en [-31,-27])\n")
    if R["hil_loc"] is not None:
        M.append(f"- Pared lumen-lumen minima pelvis/ureter <-> venoso: **{R['hil_min']*1000:.1f} um** "
                 f"{'(<400um: HAY componente de hilio)' if R['hil_min'] < 0.400 else '(>=400um: NO hay colision de hilio; los puertos a 3mm no se tocan)'}\n")
        M.append(f"- urin ~{R['hil_loc'][0].round(2)}  venoso ~{R['hil_loc'][1].round(2)} ({R['pel_ure']} segs pelvis/ureter x {R['ven_hilio']} segs venosos en banda)\n")
    M.append("## (6) Genero / handle de la malla urinaria\n")
    M.append(f"Fusiones urinarias no-adyacentes (lumen clampeado solapa -> crean asa): **{len(R['fus'])}**.\n")
    for (i, j, d, s, m) in R["fus"][:8]:
        ul = R["ulabel"]; ue2 = R["u_edges"]
        M.append(f"- edge{i} ({ul[int(ue2[i][0])]}/{ul[int(ue2[i][1])]}) x edge{j} "
                 f"({ul[int(ue2[j][0])]}/{ul[int(ue2[j][1])]}): d={d:.3f} solape={s-d:.3f} mm en {m.round(2)}\n")
    M.append("Coincide con la fusion co-mayor (calices menores del mismo mayor, solape 85um de la Entrada 015).\n")
    M.append("## (7) Numeros que quedaron cortados en 5a\n")
    M.append(f"- Porosidad: **{float(R['meta']['porosidad_pct']):.4f} %**\n")
    for k, v in R["vol_lumen"].items():
        M.append(f"- Lumen {k}: {v:.2f} mm3\n")
    M.append(f"- Gap arterial<->venoso (nodos): **{float(R['meta']['gap_av_mm']):.3f} mm**\n")
    M.append("## VEREDICTO\n")
    M.append(f"- **Sub-muestreo** (pared nativa >=400um, se separa a 100um/5b): **{b_400}/{nP} ({R['frac_sub']:.1f}%)**\n")
    M.append(f"- **Colision real** (pared nativa <400um, requiere fix de geometria): "
             f"**{b_neg+b_100+b_100_400}/{nP} ({R['frac_col']:.1f}%)**\n")
    M.append(f"- **Componente de hilio (pelvis/ureter vs venoso <400um):** "
             f"**{'SI' if R['hil_min'] < 0.400 else 'NO'}** (pared {R['hil_min']*1000:.1f} um).\n")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        f.write("".join(M))
    print("  -> reporte:", OUT)
