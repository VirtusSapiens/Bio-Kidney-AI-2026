#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa5a_fix2.py -- FIX estanqueidad 5a v2: reruteo POR-GRADIENTE (Entrada 020).
Reemplaza el solver de direccion-fija (que cayo en FALLBACK por straddle en Z). Cada nodo en
conflicto se aparta de la superficie de OTRO circuito mas cercana por +grad(SDF_keepout): las
ramas posteriores van -Z y las anteriores +Z, cada una por el camino corto.
Mantiene FIX1 (pelvis elipsoide Z_pelvis=2.0 en capa4) SIN cambios. Reruteo LOCAL en Capa 5
(NO reescribe 3a/3b/3c). NO genera 5b. Si algun nodo no converge / >MAX_DESPL / rompe
conectividad -> FALLBACK (reporta, no fuerza, no carve).
"""
import os, sys, json, time
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import capa5a_nucleo_sdf as C5

MIN_WALL = 0.400
MAX_DESPL = 5.0
MAX_ITER = 200
VX = C5.VOXEL_5A
REPO = C5.REPO
OUT_MD = os.path.join(REPO, "09_paper_vascular", "auditoria_5a_fix2.md")


def sd_ellipsoid(P, c, r):
    p = P - c
    k0 = np.linalg.norm(p / r, axis=1)
    k1 = np.linalg.norm(p / (r * r), axis=1)
    return k0 * (k0 - 1.0) / np.maximum(k1, 1e-9)


def eval_sdf_mixed(cones, ellips, origin, shape, vx):
    sdf = C5.eval_sdf(cones, origin, shape, vx)
    nx, ny, nz = shape
    for (c, r) in ellips:
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
    log("FIX v2 ESTANQUEIDAD 5a -- reruteo POR-GRADIENTE (mantiene FIX1 pelvis elipsoide)")
    fallback = []

    d3a = np.load(C5.find("capa3a_arterial.npz"), allow_pickle=True)
    d3b = np.load(C5.find("capa3b_venoso.npz"), allow_pickle=True)
    d3c = np.load(C5.find("capa3c_colector.npz"), allow_pickle=True)
    d4 = np.load(C5.find("capa4_colector_alto.npz"), allow_pickle=True)
    n3a, e3a, r3a = C5.adapt_edge_radios(d3a)
    n3b, e3b, r3b = C5.adapt_edge_radios(d3b)
    n3c, e3c, r3c = C5.adapt_3c(d3c)
    n4, e4, r4 = C5.adapt_4(d4)
    Zp = float(d4["z_pelvis"]); niv4 = d4["nivel"]
    pel_pos = n4[niv4 == 3][0]; ure_pos = n4[niv4 == 4][0]
    pel_sa = np.array([float(d4["radios"][niv4 == 3][0])] * 2 + [Zp])
    log(f"FIX1 (sin cambios): pelvis elipsoide {pel_sa} en {pel_pos.round(2)}")

    # ---- URINARIO con pelvis elipsoide (fijo) ----
    uN, uE, uR, _ = C5.merge_parts([(n3c, e3c, r3c, "3c"), (n4, e4, r4, "4")])
    pu = C5.poda_clamp(uN, uE, uR, C5.PUERTO_URINARIO)
    u_pel = int(np.argmin(np.linalg.norm(uN - pel_pos, axis=1)))
    def urin_cones():
        cs = []
        for a, b in pu["edges_ret"]:
            a, b = int(a), int(b); ra, rb = pu["radii_c"][a], pu["radii_c"][b]
            if a == u_pel: ra = rb
            if b == u_pel: rb = ra
            cs.append((uN[a], uN[b], float(ra), float(rb)))
        return cs
    UCONES = urin_cones(); UELL = [(pel_pos, pel_sa)]
    umid = np.array([0.5 * (c[0] + c[1]) for c in UCONES]); utree = cKDTree(umid)

    def build_vasc(nA, eA, rA, port):
        p = C5.poda_clamp(nA, eA, rA, port)
        parent = -np.ones(len(nA), int)
        for a, b in p["edges_ret"]: parent[int(b)] = int(a)
        return dict(N=nA.copy(), N0=nA.copy(), edges=p["edges_ret"], rc=p["radii_c"],
                    ret=set(int(x) for x in p["ret_nodes"]), port=p["port_node"], parent=parent)
    ven = build_vasc(n3b, e3b, r3b, C5.PUERTO_VENOSO)
    art = build_vasc(n3a, e3a, r3a, np.array(d3a["raiz"], float))
    for c in (ven, art):
        ch = {i: [] for i in range(len(c["N"]))}
        for a, b in c["edges"]: ch[int(a)].append(int(b))
        c["children"] = ch
        c["origlen"] = {int(b): float(np.linalg.norm(c["N0"][int(b)] - c["N0"][int(a)])) for a, b in c["edges"]}

    # ---- keep-out (SDF de "todo lo no-propio") ----
    def sd_urin(p):
        d = float(sd_ellipsoid(p[None], pel_pos, pel_sa)[0])
        for j in utree.query_ball_point(p, 9.0):
            a, b, r1, r2 = UCONES[j]
            d = min(d, float(C5.sd_round_cone(p[None], a, b, r1, r2)[0]))
        return d
    def vasc_cones(c):
        return [(c["N"][int(a)], c["N"][int(b)], float(c["rc"][int(a)]), float(c["rc"][int(b)])) for a, b in c["edges"]]
    def sd_vasc(cones, p):
        d = np.inf
        for (a, b, r1, r2) in cones:
            d = min(d, float(C5.sd_round_cone(p[None], a, b, r1, r2)[0]))
        return d
    def keepout(p, other_cones):
        return min(sd_urin(p), sd_vasc(other_cones, p))
    def wall(p, r, other_cones):
        return keepout(p, other_cones) - r
    def grad(p, other_cones):
        h = 0.05; g = np.zeros(3)
        for k in range(3):
            e = np.zeros(3); e[k] = h
            g[k] = (keepout(p + e, other_cones) - keepout(p - e, other_cones)) / (2 * h)
        n = np.linalg.norm(g)
        return g / n if n > 1e-9 else np.array([0.0, 0.0, 1.0])

    # ---- PASO 1: identificar conflictos iniciales (terminales Y in [-31,-26]) ----
    def conflicts(c, other_cones):
        out = []
        for i in c["ret"]:
            p = c["N"][i]
            if -31 <= p[1] <= -26 and wall(p, c["rc"][i], other_cones) < MIN_WALL - 1e-4:
                out.append(i)
        return out
    art_cones0 = vasc_cones(art); ven_cones0 = vasc_cones(ven)
    v_conf0 = conflicts(ven, art_cones0); a_conf0 = conflicts(art, ven_cones0)
    log(f"PASO1 keep-out: venoso {len(v_conf0)} nodos en conflicto; arterial {len(a_conf0)} nodos en conflicto")

    # ---- PASO 2: solver por-gradiente ----
    def conn_ok(c, i, newp):
        pa = c["parent"][i]
        if pa >= 0:
            L0 = c["origlen"].get(i, 1.0)
            d = newp - c["N"][pa]; L = np.linalg.norm(d)
            if L > 1.5 * max(L0, 1e-6): return False
            if np.dot(d, c["N0"][i] - c["N0"][pa]) <= 0: return False  # invertida
        for b in c["children"][i]:
            L0 = c["origlen"].get(b, 1.0)
            if np.linalg.norm(c["N"][b] - newp) > 1.5 * max(L0, 1e-6): return False
        return True

    def step_node(c, i, other_cones):
        p = c["N"][i]; r = c["rc"][i]
        w = wall(p, r, other_cones)
        if w >= MIN_WALL: return False
        g = grad(p, other_cones)
        deficit = MIN_WALL - w
        st = min(deficit + 0.05, 0.4)
        while st > 1e-3:
            newp = p + g * st
            if np.linalg.norm(newp - c["N0"][i]) <= MAX_DESPL and conn_ok(c, i, newp):
                c["N"][i] = newp
                return True
            st *= 0.5
        return False

    t0 = time.perf_counter()
    it = 0
    for it in range(MAX_ITER):
        art_cones = vasc_cones(art); ven_cones = vasc_cones(ven)
        moved = 0
        for i in list(ven["ret"]):
            if step_node(ven, i, art_cones): moved += 1
        for i in list(art["ret"]):
            if step_node(art, i, ven_cones): moved += 1
        if moved == 0:
            break
    dt_solve = time.perf_counter() - t0

    # verificar convergencia final (vs TODAS las superficies de otro circuito)
    art_cones = vasc_cones(art); ven_cones = vasc_cones(ven)
    def residuals(c, other_cones):
        return [i for i in c["ret"] if wall(c["N"][i], c["rc"][i], other_cones) < MIN_WALL - 1e-3]
    v_res = residuals(ven, art_cones); a_res = residuals(art, ven_cones)
    def stats(c):
        disp = np.linalg.norm(c["N"] - c["N0"], axis=1)
        mv = np.where(disp > 1e-4)[0]
        return len(mv), float(disp.max()) if len(mv) else 0.0, float(disp[list(c["ret"])].mean()), float(disp[c["port"]])
    vn, vmax, vmean, vport = stats(ven); an, amax, amean, aport = stats(art)
    log(f"PASO2 solver ({it+1} iter, {dt_solve:.1f}s):")
    log(f"   venoso  : {vn} nodos movidos, despl max {vmax:.2f}mm medio {vmean:.2f}mm puerto {vport:.2f}mm; residual {len(v_res)}")
    log(f"   arterial: {an} nodos movidos, despl max {amax:.2f}mm medio {amean:.2f}mm puerto {aport:.2f}mm; residual {len(a_res)}")

    if v_res or a_res or vmax > MAX_DESPL or amax > MAX_DESPL:
        log("\n[FALLBACK] el reruteo por-gradiente NO convergio a MIN_WALL para todos los nodos.")
        for nombre, c, res, other in (("venoso", ven, v_res, art_cones), ("arterial", art, a_res, ven_cones)):
            for i in res[:5]:
                p = c["N"][i]; du = sd_urin(p); dv = sd_vasc(other, p)
                w = min(du, dv) - c["rc"][i]
                sup = "urinario(pelvis/ureter/calices)" if du < dv else "otro-arbol-vascular"
                log(f"   {nombre} nodo {i} en {p.round(2)}: pared {w*1000:.0f}um (faltan {(MIN_WALL-w)*1000:.0f}um) "
                    f"vs {sup}; grad estancado (atrapado entre superficies).")
        log("   (decision del usuario: carve local vs revisar el hilio; NO aplico nada)")
        _write_md(LOG, None); _save(None); return LOG

    # ---- PASO 3: re-mallar + re-auditar ----
    def mesh_circ(cones, ellip, nodes_pts):
        pts = [nodes_pts] + [c[None] + np.array([[r[0], r[1], r[2]]]) for (c, r) in ellip] + \
              [c[None] - np.array([[r[0], r[1], r[2]]]) for (c, r) in ellip]
        P = np.vstack(pts)
        rmax = max([max(k[2], k[3]) for k in cones] + [float(r.max()) for (_, r) in ellip] + [0.3])
        lo = P.min(0) - (rmax + C5.BBOX_MARGIN); hi = P.max(0) + (rmax + C5.BBOX_MARGIN)
        shp = tuple(int(np.ceil((hi[i] - lo[i]) / VX)) + 1 for i in range(3))
        sdf = eval_sdf_mixed(cones, ellip, lo, shp, VX)
        vol = int((sdf < 0).sum()) * VX ** 3
        from skimage import measure
        v, f, _, _ = measure.marching_cubes(sdf, level=0.0, spacing=(VX, VX, VX)); v = v + lo
        return v, f, vol, sdf, lo, shp
    circuits = {
        "urinaria": dict(cones=UCONES, ellip=UELL, N=uN, ret=pu["ret_nodes"], edges=pu["edges_ret"], port=pu["port_node"]),
        "vascular_0_arterial": dict(cones=vasc_cones(art), ellip=[], N=art["N"], ret=np.array(list(art["ret"])), edges=art["edges"], port=art["port"]),
        "vascular_1_venoso": dict(cones=vasc_cones(ven), ellip=[], N=ven["N"], ret=np.array(list(ven["ret"])), edges=ven["edges"], port=ven["port"]),
    }
    meshes = {}
    for nombre, c in circuits.items():
        v, f, vol, sdf, lo, shp = mesh_circ(c["cones"], c["ellip"], c["N"][np.array(list(c["ret"]))] if not isinstance(c["ret"], np.ndarray) else c["N"][c["ret"]])
        topo = C5.mesh_topology(v, f)
        vc, fc, nout = C5.clip_to_domain(v, f)
        C5.write_stl(os.path.join(REPO, f"capa5a_fix2_lumen_{nombre}.stl"), vc, fc)
        meshes[nombre] = dict(v=v, f=f, vol=vol, topo=topo, sdf=sdf, lo=lo, shp=shp)

    # estanqueidad grid comun
    up = np.vstack([np.vstack([k[0], k[1]]) for k in UCONES] + [pel_pos + pel_sa, pel_pos - pel_sa])
    vcs = circuits["vascular_0_arterial"]["cones"] + circuits["vascular_1_venoso"]["cones"]
    vp = np.vstack([np.vstack([k[0], k[1]]) for k in vcs])
    allp = np.vstack([up, vp]); lo = allp.min(0) - 6.5; hi = allp.max(0) + 6.5
    shp = tuple(int(np.ceil((hi[i] - lo[i]) / VX)) + 1 for i in range(3))
    su = eval_sdf_mixed(UCONES, UELL, lo, shp, VX)
    sa = C5.eval_sdf(circuits["vascular_0_arterial"]["cones"], lo, shp, VX)
    sv = C5.eval_sdf(circuits["vascular_1_venoso"]["cones"], lo, shp, VX)
    both = int(((su < 0) & ((sa < 0) | (sv < 0))).sum())
    both_a = int(((su < 0) & (sa < 0)).sum()); both_v = int(((su < 0) & (sv < 0)).sum())
    # pared minima superficie-superficie: verts urinario -> SDF vascular
    vu = meshes["urinaria"]["v"]
    svasc = np.minimum(sa, sv)
    ijk = np.clip(((vu - lo) / VX).round().astype(int), 0, np.array(shp) - 1)
    dvals = svasc[ijk[:, 0], ijk[:, 1], ijk[:, 2]]
    kmin = int(np.argmin(dvals)); wmin = float(dvals[kmin]); wloc = vu[kmin]

    log("\nPASO3 AUDITORIAS:")
    log(f" (1) ESTANQUEIDAD: voxeles AMBOS = {both} (urin&art {both_a}, urin&ven {both_v})  ANTES 1118  {'[OK] meta 0' if both == 0 else '[REVISAR]'}")
    log(f"     pared minima inter-circuito (superf urin->vascular): {wmin*1000:.0f}um en {wloc.round(2)}  {'[OK] >=400um' if wmin >= MIN_WALL - 0.02 else '[REVISAR]'}")
    log(" (2) WATERTIGHT:")
    for nombre, m in meshes.items():
        t = m["topo"]; ok = t["manifold"] and t["boundary"] == 0 and t["nonmanifold"] == 0
        log(f"     {nombre}: manifold {t['manifold']} bordes {t['boundary']} genero {t['genus']} {'[OK]' if ok else '[REVISAR]'}"
            + ("  (genero 1 co-mayor 85um -> 5b, no tocar)" if nombre == "urinaria" and t["genus"] > 0 else ""))
    log(" (3) CONECTIVIDAD:")
    for nombre, c in circuits.items():
        ncc, ll = C5.comps(len(c["N"]), c["edges"]); ret = c["ret"] if isinstance(c["ret"], np.ndarray) else np.array(list(c["ret"]))
        big = np.bincount(ll[ret]).argmax(); reach = ll[c["port"]] == big
        log(f"     {nombre}: puerto alcanzado {'[OK]' if reach else '[REVISAR]'}")
    log(" (4) CONTENCION en elipsoide Capa 0 (incl. nodos movidos):")
    for nombre, c in circuits.items():
        ret = c["ret"] if isinstance(c["ret"], np.ndarray) else np.array(list(c["ret"]))
        val = C5.in_ellipsoid(c["N"][ret]); viol = int((val > 1 + 1e-9).sum())
        log(f"     {nombre}: violaciones {viol}/{len(ret)} (max {val.max():.3f}){' (puerto salida)' if viol else ' [OK]'}")
    vmat = 4 / 3 * np.pi * np.prod(C5.DOMINIO_SEMIEJES); vt = sum(m["vol"] for m in meshes.values())
    log(" (5) VOLUMENES:")
    for nombre, m in meshes.items(): log(f"     lumen {nombre}: {m['vol']:.2f} mm3")
    log(f"     matriz {vmat:.1f}  lumen total {vt:.2f}  porosidad {100*vt/vmat:.4f} %")
    log(f"\n desplazamientos: venoso {vn} nodos (max {vmax:.2f}mm), arterial {an} nodos (max {amax:.2f}mm)")

    res = dict(both=both, both_a=both_a, both_v=both_v, wmin=wmin, wloc=wloc, meshes=meshes,
               vn=vn, vmax=vmax, vmean=vmean, vport=vport, an=an, amax=amax, amean=amean, aport=aport,
               vmat=vmat, vt=vt, v_conf0=len(v_conf0), a_conf0=len(a_conf0), it=it + 1)
    _write_md(LOG, res); _save(res)
    return LOG


def _save(res):
    if res is None: return
    np.savez_compressed(os.path.join(REPO, "capa5a_fix2_meta.npz"),
                        min_wall=MIN_WALL, both_voxels=res["both"], pared_min_um=res["wmin"] * 1000,
                        ven_nodos=res["vn"], ven_despl_max=res["vmax"], ven_puerto=res["vport"],
                        art_nodos=res["an"], art_despl_max=res["amax"], art_puerto=res["aport"],
                        porosidad=100 * res["vt"] / res["vmat"], vol_total=res["vt"])


def _write_md(LOG, res):
    M = ["# Auditoria fix estanqueidad Capa 5a v2 (reruteo por-gradiente)\n",
         "**Programa:** Bio-Kidney AI 2026 - **Capa:** 5a (fix v2) - **Fecha:** 2026-07-09\n",
         "Reruteo POR-GRADIENTE (cada nodo se aparta de la superficie de otro circuito mas cercana, "
         "no direccion fija) manteniendo FIX1 (pelvis elipsoide). Local en Capa 5, NO reescribe 3a/3b/3c. "
         "MIN_WALL=0.400 mm, MAX_DESPL=5.0 mm.\n\n```\n"]
    M += [l + "\n" for l in LOG]
    M.append("```\n")
    os.makedirs(os.path.dirname(OUT_MD), exist_ok=True)
    with open(OUT_MD, "w") as f: f.write("".join(M))


if __name__ == "__main__":
    main()
    print("  -> reporte:", OUT_MD)
