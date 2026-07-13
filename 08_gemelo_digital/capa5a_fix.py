#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa5a_fix.py -- FIX de estanqueidad de 5a (Entrada 019).
FIX1: pelvis elipsoide aplanado (ya en capa4_colector_alto.npz: tipo_primitiva/semieje_z).
FIX2: reruteo LOCAL del tronco venoso terminal a +Z (fuera del keep-out pelvis/ureter).
FIX3: separacion del puerto/tramo arterial a -Z respecto al ureter.
Los reruteos son transformaciones LOCALES en Capa 5 (NO reescriben 3a/3b/3c). Re-malla
grueso (200um) y re-audita. Si el solver no logra MIN_WALL sin mover puertos>5mm o romper
conectividad -> FALLBACK (reporta, NO fuerza, NO carve).
"""
import os, sys, json, time
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import capa5a_nucleo_sdf as C5

MIN_WALL = 0.400
VX = C5.VOXEL_5A
REPO = C5.REPO
DZMAX = 5.0
MAXPASS = 60
OUT_MD = os.path.join(REPO, "09_paper_vascular", "auditoria_5a_fix.md")
EZ = np.array([0.0, 0.0, 1.0])


def sd_ellipsoid(P, c, r):
    p = P - c
    k0 = np.linalg.norm(p / r, axis=1)
    k1 = np.linalg.norm(p / (r * r), axis=1)
    return k0 * (k0 - 1.0) / np.maximum(k1, 1e-9)


def eval_sdf_mixed(cones, ellips, origin, shape, vx):
    sdf = C5.eval_sdf(cones, origin, shape, vx)
    nx, ny, nz = shape
    for (c, r) in ellips:
        rmax = float(r.max()) + 1.5 * vx
        lo = c - (r + 1.5 * vx); hi = c + (r + 1.5 * vx)
        i0 = np.clip(np.floor((lo - origin) / vx).astype(int), 0, [nx, ny, nz])
        i1 = np.clip(np.ceil((hi - origin) / vx).astype(int) + 1, 0, [nx, ny, nz])
        if np.any(i1 <= i0):
            continue
        xs = origin[0] + np.arange(i0[0], i1[0]) * vx
        ys = origin[1] + np.arange(i0[1], i1[1]) * vx
        zs = origin[2] + np.arange(i0[2], i1[2]) * vx
        gx, gy, gz = np.meshgrid(xs, ys, zs, indexing='ij')
        Pg = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
        d = sd_ellipsoid(Pg, c, r).reshape(gx.shape).astype(np.float32)
        sub = sdf[i0[0]:i1[0], i0[1]:i1[1], i0[2]:i1[2]]
        np.minimum(sub, d, out=sub)
    return sdf


def main():
    LOG = []
    def log(s=""): print(s); LOG.append(s)
    log("FIX ESTANQUEIDAD 5a (pelvis elipsoide + reruteo venoso/arterial)")
    fallback = None

    d3a = np.load(C5.find("capa3a_arterial.npz"), allow_pickle=True)
    d3b = np.load(C5.find("capa3b_venoso.npz"), allow_pickle=True)
    d3c = np.load(C5.find("capa3c_colector.npz"), allow_pickle=True)
    d4 = np.load(C5.find("capa4_colector_alto.npz"), allow_pickle=True)
    n3a, e3a, r3a = C5.adapt_edge_radios(d3a)
    n3b, e3b, r3b = C5.adapt_edge_radios(d3b)
    n3c, e3c, r3c = C5.adapt_3c(d3c)
    n4, e4, r4 = C5.adapt_4(d4)
    n3ab = np.load(C5.find("capa3ab_peritubular.npz"), allow_pickle=True)["puntos_drenaje"].astype(float)
    Zp = float(d4["z_pelvis"]); niv4 = d4["nivel"]
    pel_pos = n4[niv4 == 3][0]; ure_pos = n4[niv4 == 4][0]; may_pos = n4[niv4 == 2]
    pel_sa = np.array([float(d4["radios"][niv4 == 3][0]), float(d4["radios"][niv4 == 3][0]), Zp])
    log(f"FIX1 pelvis elipsoide {pel_sa} en {pel_pos.round(2)}  (Z_pelvis={Zp})")

    # ---- URINARIO con pelvis elipsoide ----
    uN, uE, uR, _ = C5.merge_parts([(n3c, e3c, r3c, "3c"), (n4, e4, r4, "4")])
    pu = C5.poda_clamp(uN, uE, uR, C5.PUERTO_URINARIO)
    # localizar nodos macro en el merge (pelvis, ureter, mayores) por posicion
    def nearest(P, q): return int(np.argmin(np.linalg.norm(P - q, axis=1)))
    u_pel = nearest(uN, pel_pos); u_ure = nearest(uN, ure_pos)
    # cones urinarios: los incidentes a la pelvis usan el radio del OTRO extremo (no 4.5)
    def urin_cones_ellip(edges, N, rc):
        cones = []
        for a, b in edges:
            a, b = int(a), int(b)
            ra, rb = rc[a], rc[b]
            if a == u_pel: ra = rb
            if b == u_pel: rb = ra
            cones.append((N[a], N[b], float(ra), float(rb)))
        return cones
    u_cones = urin_cones_ellip(pu["edges_ret"], uN, pu["radii_c"])
    u_ellip = [(pel_pos, pel_sa)]
    # keep-out macro urinario para el solver (elipsoide + infundibulos + ureter)
    keep_cones = []
    for a, b in pu["edges_ret"]:
        a, b = int(a), int(b)
        # niveles: mapear a nivel4 por posicion
        lo_pel = (a == u_pel or b == u_pel)
        # solo macro: mayor(r2)->pelvis y pelvis->ureter
        if lo_pel:
            other = b if a == u_pel else a
            r_other = float(pu["radii_c"][other])
            keep_cones.append((uN[a], uN[b], r_other, r_other))
    log(f"URINARIO: {len(uN)} nodos, {len(pu['edges_ret'])} aristas retenidas, pelvis nodo {u_pel}, keep-out cones {len(keep_cones)}")

    # ---- circuitos vasculares ----
    def build_vasc(nA, eA, rA, port):
        p = C5.poda_clamp(nA, eA, rA, port)
        return dict(N=nA.copy(), N0=nA.copy(), edges=p["edges_ret"], rc=p["radii_c"],
                    ret=p["ret_nodes"], port=p["port_node"])
    ven = build_vasc(n3b, e3b, r3b, C5.PUERTO_VENOSO)
    art = build_vasc(n3a, e3a, r3a, np.array(d3a["raiz"], float))

    # ---- SOLVER de reruteo (wall vs pelvis elipsoide + keep_cones ; empuje en signo dado) ----
    def wall_at(p, rv):
        d = sd_ellipsoid(p[None, :], pel_pos, pel_sa)[0]
        for (a, b, r1, r2) in keep_cones:
            d = min(d, float(C5.sd_round_cone(p[None, :], a, b, r1, r2)[0]))
        return d - rv
    def solve_dz(p, rv, sign):
        if wall_at(p, rv) >= MIN_WALL: return 0.0
        hi = DZMAX
        if wall_at(p + sign * hi * EZ, rv) < MIN_WALL: return DZMAX + 99  # infeasible
        loo = 0.0
        for _ in range(28):
            mid = 0.5 * (loo + hi)
            if wall_at(p + sign * mid * EZ, rv) >= MIN_WALL: hi = mid
            else: loo = mid
        return hi

    def dz_needed(p, rv, sign, cap=20.0):
        if wall_at(p, rv) >= MIN_WALL: return 0.0
        hi = cap
        if wall_at(p + sign * hi * EZ, rv) < MIN_WALL: return np.inf
        loo = 0.0
        for _ in range(28):
            mid = 0.5 * (loo + hi)
            if wall_at(p + sign * mid * EZ, rv) >= MIN_WALL: hi = mid
            else: loo = mid
        return hi

    def reroute(circ, sign, nombre):
        """(1) MIDE req por nodo desde ORIGINAL (reporte). (2) APLICA moves iterativos con
        desplazamiento CUMULATIVO capado a 5mm y suavizado ligero (evita overshoot). (3) verifica."""
        nonlocal fallback
        N0 = circ["N0"]; ret = circ["ret"]; rc = circ["rc"]; edges = circ["edges"]
        # medicion desde original
        req = {}; infeas = []
        for i in ret:
            i = int(i); d = dz_needed(N0[i], rc[i], sign)
            if d <= 0: continue
            req[i] = d
            if not np.isfinite(d) or d > 5.0: infeas.append(i)
        vals = [v for v in req.values() if np.isfinite(v)]
        port_req = req.get(int(circ["port"]), 0.0)
        log(f"REROUTE {nombre} ({sign:+.0f}Z): {len(req)} nodos a mover  max {max(vals) if vals else 0:.2f}mm  "
            f"mediana {np.median(vals) if vals else 0:.2f}mm  puerto {port_req if np.isfinite(port_req) else 999:.2f}mm  "
            f"infeasibles(>5mm): {len(infeas)}")
        circ["req"] = req; circ["infeas"] = infeas; circ["port_req"] = float(port_req) if np.isfinite(port_req) else 99.0
        if len(infeas) or (np.isfinite(port_req) and port_req > 5.0):
            fallback = fallback or []
            fallback.append(dict(circuito=nombre, n_mover=len(req), infeas=len(infeas),
                                 max_req=max(vals) if vals else 0.0, port_req=float(port_req) if np.isfinite(port_req) else 99.0))
            return circ
        # APLICAR (bounded): cumul capado a 5mm, suavizado (0.3,0.15), max 25 pasadas
        neigh = {i: set() for i in range(len(N0))}
        for a, b in edges: neigh[int(a)].add(int(b)); neigh[int(b)].add(int(a))
        N = circ["N"]; cumul = np.zeros(len(N))
        for passo in range(25):
            viol = [int(i) for i in ret if wall_at(N[int(i)], rc[int(i)]) < MIN_WALL - 1e-4]
            if not viol: break
            disp = np.zeros(len(N))
            for i in viol:
                d = dz_needed(N[i], rc[i], sign, cap=5.0)
                if not np.isfinite(d): d = 5.0 - cumul[i]
                d = max(0.0, min(d, 5.0 - cumul[i]))
                disp[i] = max(disp[i], d)
                h1 = neigh[i]; h2 = {w for s in h1 for w in neigh[s]} - h1 - {i}
                for w in h1: disp[w] = max(disp[w], 0.3 * d)
                for w in h2: disp[w] = max(disp[w], 0.15 * d)
            disp = np.minimum(disp, 5.0 - cumul)
            N[:, 2] += sign * disp; cumul += disp
        remain = [int(i) for i in ret if wall_at(N[int(i)], rc[int(i)]) < MIN_WALL - 1e-3]
        moved = set(int(i) for i in np.where(cumul > 1e-6)[0])
        circ["remain_nodes"] = remain
        circ["max_move"] = float(cumul.max()); circ["port_move"] = float(cumul[circ["port"]]); circ["n_moved"] = len(moved)
        log(f"   aplicado: {len(moved)} nodos movidos, max cumul {cumul.max():.2f}mm, puerto {cumul[circ['port']]:.2f}mm, restan violando {len(remain)}")
        if len(remain) or cumul[circ["port"]] > 5.0:
            fallback = fallback or []
            fallback.append(dict(circuito=nombre, n_mover=len(req), infeas=len(remain),
                                 max_req=float(cumul.max()), port_req=float(cumul[circ["port"]])))
        return circ

    reroute(ven, +1.0, "venoso")
    reroute(art, -1.0, "arterial")

    if fallback:
        log("\n[FALLBACK] la reruteo de direccion-fija NO logra MIN_WALL con movimientos <=5mm.")
        log("CAUSA GEOMETRICA: los arboles vascular(es) CRUZAN la pelvis por AMBOS lados (Z+ y Z-);")
        log("empujar todo un arbol en un solo sentido obliga a nodos anteriores a atravesar la pelvis.")
        seen = set()
        for fb in fallback:
            if fb["circuito"] in seen: continue
            seen.add(fb["circuito"])
            circ = ven if fb["circuito"] == "venoso" else art
            rem = circ.get("remain_nodes", [])
            zp = [int(circ["N0"][i][2] > 0) for i in rem]   # 1=anterior Z+, 0=posterior Z-
            n_ant = sum(zp); n_pos = len(zp) - n_ant
            log(f"   {fb['circuito']}: tras aplicar (cap 5mm) restan {len(rem)} nodos violando MIN_WALL "
                f"(puerto movido {circ.get('port_move',0):.2f}mm, max cumul {circ.get('max_move',0):.2f}mm)")
            log(f"      STRADDLE: de los {len(rem)} residuales, {n_ant} ANTERIORES (Z+) y {n_pos} POSTERIORES (Z-) "
                f"de la pelvis -> empujarlos todos a un lado obliga a cruzarla.")
            for i in rem[:4]:
                p = circ["N0"][i]; de = float(sd_ellipsoid(p[None], pel_pos, pel_sa)[0])
                w = wall_at(circ["N"][i], circ["rc"][i])
                log(f"      nodo {i} orig {p.round(2)} (Z={p[2]:+.2f}): pared final {w*1000:.0f}um "
                    f"(faltan {(MIN_WALL-w)*1000:.0f}um); a pelvis-elipsoide {de:.2f}mm")
        log("\n   OPCIONES (NO aplico ninguna; las decidis vos):")
        log("   (a) mas aplanamiento de pelvis (Z_pelvis<2.0) -> reduce el keep-out en Z.")
        log("   (b) mover el nodo pelvis en -Y (hundir hacia el seno) -> la aleja del cruce vascular.")
        log("   (c) reruteo POR-GRADIENTE (cada nodo se aparta de la superficie mas cercana, no direccion fija).")
        log("   (d) carve keep-out en la matriz (NO por defecto; preferis decidirlo).")
        log("   FIX1 (pelvis elipsoide en capa4) QUEDA aplicado y verificado (byte-identico); solo el reruteo espera decision.")
        _write_md(LOG, fallback=fallback)
        return LOG, dict(fallback=fallback)

    # ---- RE-MALLAR + RE-AUDITAR (grueso 200um) ----
    t0 = time.perf_counter()
    circuits = {"urinaria": dict(cones=u_cones, ellip=u_ellip, N=uN, ret=pu["ret_nodes"], edges=pu["edges_ret"],
                                 rc=pu["radii_c"], port=pu["port_node"], portpos=C5.PUERTO_URINARIO),
                "vascular_0_arterial": dict(cones=C5.build_cones(art["N"], art["edges"], art["rc"]), ellip=[],
                                            N=art["N"], ret=art["ret"], edges=art["edges"], rc=art["rc"],
                                            port=art["port"], portpos=np.array(d3a["raiz"], float)),
                "vascular_1_venoso": dict(cones=C5.build_cones(ven["N"], ven["edges"], ven["rc"]), ellip=[],
                                          N=ven["N"], ret=ven["ret"], edges=ven["edges"], rc=ven["rc"],
                                          port=ven["port"], portpos=C5.PUERTO_VENOSO)}
    meshes = {}
    for nombre, c in circuits.items():
        allp = np.vstack([np.vstack([k[0], k[1]]) for k in c["cones"]] + ([e[0][None] for e in []]))
        pts = [c["N"][c["ret"]]]
        for (cc, rr) in c["ellip"]: pts.append(cc[None] + np.array([[rr[0], rr[1], rr[2]]])); pts.append(cc[None] - np.array([[rr[0], rr[1], rr[2]]]))
        P = np.vstack(pts)
        rmax = max([max(k[2], k[3]) for k in c["cones"]] + [float(rr.max()) for (_, rr) in c["ellip"]] + [0.3])
        lo = P.min(0) - (rmax + C5.BBOX_MARGIN); hi = P.max(0) + (rmax + C5.BBOX_MARGIN)
        shp = tuple(int(np.ceil((hi[i] - lo[i]) / VX)) + 1 for i in range(3))
        sdf = eval_sdf_mixed(c["cones"], c["ellip"], lo, shp, VX)
        vol = int((sdf < 0).sum()) * VX ** 3
        from skimage import measure
        v, f, _, _ = measure.marching_cubes(sdf, level=0.0, spacing=(VX, VX, VX)); v = v + lo
        topo = C5.mesh_topology(v, f)
        vc, fc, nout = C5.clip_to_domain(v, f)
        C5.write_stl(os.path.join(REPO, f"capa5a_fix_lumen_{nombre}.stl"), vc, fc)
        meshes[nombre] = dict(v=v, f=f, vol=vol, topo=topo, sdf=sdf, lo=lo, shp=shp)
        log(f"MALLA {nombre}: verts {len(v)} faces {len(f)} vol {vol:.1f}mm3 genero {topo['genus']} "
            f"manifold {topo['manifold']} bordes {topo['boundary']}")

    # estanqueidad en grid comun
    up = np.vstack([np.vstack([k[0], k[1]]) for k in u_cones]); up = np.vstack([up, pel_pos + pel_sa, pel_pos - pel_sa])
    vp = np.vstack([np.vstack([k[0], k[1]]) for k in circuits["vascular_0_arterial"]["cones"] + circuits["vascular_1_venoso"]["cones"]])
    allp = np.vstack([up, vp]); rmax = 4.5
    lo = allp.min(0) - (rmax + C5.BBOX_MARGIN); hi = allp.max(0) + (rmax + C5.BBOX_MARGIN)
    shp = tuple(int(np.ceil((hi[i] - lo[i]) / VX)) + 1 for i in range(3))
    su = eval_sdf_mixed(u_cones, u_ellip, lo, shp, VX)
    sv = C5.eval_sdf(circuits["vascular_0_arterial"]["cones"] + circuits["vascular_1_venoso"]["cones"], lo, shp, VX)
    both = int(((su < 0) & (sv < 0)).sum())
    # pared minima inter-circuito con localizacion (seg-seg centerline - radios)
    # min sobre pares urin(macro)xvasc de (dist - r_u - r_v); ellipsoide aprox por su bbox de cones
    def min_wall_loc():
        best = np.inf; loc = None
        vcones = circuits["vascular_0_arterial"]["cones"] + circuits["vascular_1_venoso"]["cones"]
        # muestreo: verts de la malla urinaria vs superficie vascular via SDF
        vu = meshes["urinaria"]["v"]
        dv = C5.eval_sdf(vcones, lo, shp, VX)
        # distancia de cada vert urinario a la superficie vascular = interp SDF ~ nearest grid
        ijk = np.clip(((vu - lo) / VX).round().astype(int), 0, np.array(shp) - 1)
        dvals = dv[ijk[:, 0], ijk[:, 1], ijk[:, 2]]
        k = int(np.argmin(dvals))
        return float(dvals[k]), vu[k]
    wmin, wloc = min_wall_loc()
    dt = time.perf_counter() - t0
    log(f"\nESTANQUEIDAD: voxeles en AMBOS = {both}  (5a: 1118)  {'[OK] meta 0' if both == 0 else '[REVISAR]'}")
    log(f"  pared minima inter-circuito (superficie urin -> superficie vascular): {wmin*1000:.0f}um "
        f"en {wloc.round(2)}  {'[OK] >=400um' if wmin >= MIN_WALL - 1e-3 else '[REVISAR]'}")

    # auditorias restantes
    log("WATERTIGHT:")
    for nombre, m in meshes.items():
        t = m["topo"]; ok = t["manifold"] and t["boundary"] == 0 and t["nonmanifold"] == 0
        log(f"  {nombre}: manifold {t['manifold']} bordes {t['boundary']} genero {t['genus']} {'[OK]' if ok else '[REVISAR]'}")
    log("CONECTIVIDAD:")
    for nombre, c in circuits.items():
        ncc, ll = C5.comps(len(c["N"]), c["edges"])
        big = np.bincount(ll[c["ret"]]).argmax()
        reach = ll[c["port"]] == big
        log(f"  {nombre}: puerto alcanzado {'[OK]' if reach else '[REVISAR]'}")
    log("CONTENCION en elipsoide Capa 0:")
    for nombre, c in circuits.items():
        val = C5.in_ellipsoid(c["N"][c["ret"]]); viol = int((val > 1 + 1e-9).sum())
        log(f"  {nombre}: violaciones {viol}/{len(c['ret'])} (max {val.max():.3f})"
            f"{' (puerto salida)' if viol else ' [OK]'}")
    vmat = 4 / 3 * np.pi * np.prod(C5.DOMINIO_SEMIEJES)
    vt = sum(m["vol"] for m in meshes.values())
    log("VOLUMENES:")
    for nombre, m in meshes.items(): log(f"  lumen {nombre}: {m['vol']:.2f} mm3")
    log(f"  matriz {vmat:.1f} mm3  lumen total {vt:.2f}  porosidad {100*vt/vmat:.4f} %")
    log(f"\nTIEMPO re-mallado+audit: {dt:.1f}s")

    np.savez_compressed(os.path.join(REPO, "capa5a_fix_meta.npz"),
                        min_wall=MIN_WALL, z_pelvis=Zp, both_voxels=both, pared_min_um=wmin * 1000,
                        ven_max_move=ven["max_move"], ven_port_move=ven["port_move"],
                        art_max_move=art["max_move"], art_port_move=art["port_move"],
                        porosidad=100 * vt / vmat, vol_total=vt)
    _write_md(LOG, ok=True, both=both, wmin=wmin, ven=ven, art=art, meshes=meshes, Zp=Zp,
              vmat=vmat, vt=vt, val_may=[float(sd_ellipsoid(m[None], pel_pos, pel_sa)[0]) for m in may_pos])
    return LOG, dict(both=both, wmin=wmin)


def _write_md(LOG, ok=False, fallback=None, **kw):
    M = ["# Auditoria fix estanqueidad Capa 5a\n",
         "**Programa:** Bio-Kidney AI 2026 - **Capa:** 5a (fix) - **Fecha:** 2026-07-09\n",
         "FIX1 pelvis elipsoide aplanado (Capa 4) + FIX2 reruteo venoso +Z + FIX3 arterial -Z "
         "(local en Capa 5, sin reescribir 3a/3b/3c). MIN_WALL=0.400 mm.\n\n```\n"]
    M += [l + "\n" for l in LOG]
    M.append("```\n")
    os.makedirs(os.path.dirname(OUT_MD), exist_ok=True)
    with open(OUT_MD, "w") as f:
        f.write("".join(M))


if __name__ == "__main__":
    main()
    print("  -> reporte:", OUT_MD)
