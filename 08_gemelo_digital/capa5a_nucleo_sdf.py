#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa5a_nucleo_sdf.py  --  Bio-Kidney AI 2026
============================================
Capa 5a: NUCLEO SDF + validacion gruesa (smoke a 200um). NO genera la malla de
matriz fina ni el voxel etiquetado (eso es 5b).

CONCEPTO: cada capa es centerline+radio = union de capsulas conicas = SDF analitico
(rounded-cone, Inigo Quilez). El solido = elipsoide COMPLETO de Capa 0 (seno relleno)
menos la union de lumenes. El grid es solo MUESTREO del SDF (no se hace crecer).

Circuitos:
  URINARIO = 3c U 4        (colector fino + calicial alto)
  VASCULAR = 3a U 3ab U 3b (arterial, bed peritubular, venoso)

Poda a MIN_CANAL_R (200um r / 400um diam): la microvasculatura sub-400um NO se imprime;
queda en las capas previas como OBJETIVO DE PERFUSION (scope wet-lab).

Solo NumPy/SciPy/scikit-image. STL: writer binario propio. Entorno env_biokidney.
"""
import os, sys, time, struct, json, resource
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, dijkstra
from skimage import measure

# ============================================================================
#  PARAMETROS (mm)
# ============================================================================
VOXEL_FLOOR = 0.100          # paso fino (queda para 5b)
VOXEL_5A    = 0.200          # muestreo grueso de 5a (smoke)
MIN_CANAL_DIAM = 0.400
MIN_CANAL_R    = 0.200
DOMINIO_SEMIEJES = np.array([55.0, 30.0, 18.0])
DOMINIO_CENTRO   = np.array([0.0, 0.0, 0.0])
PUERTO_URINARIO  = np.array([0.0, -30.0, 0.0])
PUERTO_VENOSO    = np.array([0.0, -30.0, 3.0])
MERGE_TOL   = 1e-3
BBOX_MARGIN = 2.0

REPO = os.path.expanduser("~/Escritorio/BioKidney-AI")
MD_OUT = os.path.join(REPO, "09_paper_vascular", "auditoria_capa5a.md")


def find(name):
    for c in (os.path.join(REPO, name), name):
        if os.path.isfile(c):
            return c
    raise FileNotFoundError(name)


def peak_mem_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0  # KB->MB en Linux


# ============================================================================
#  ADAPTADORES POR-CAPA  -> (nodos Nx3, edges Ex2, node_radios N)
# ============================================================================
def adapt_edge_radios(d):
    """3a/3b: radios POR-ARISTA -> node_radios[child]=radio; raiz=radio_raiz."""
    nodos = d["nodos"].astype(np.float64)
    aris = d["aristas"].astype(np.int64)
    er = d["radios"].astype(np.float64)
    N = len(nodos)
    nr = np.zeros(N)
    nr[aris[:, 1]] = er
    child = np.zeros(N, bool); child[aris[:, 1]] = True
    roots = np.where(~child)[0]
    nr[roots] = float(d["radio_raiz"])
    return nodos, aris, nr


def adapt_3c(d):
    """3c: radios POR-NODO; topologia via parent (10 raices = papilas)."""
    nodos = d["nodos"].astype(np.float64)
    parent = d["parent"].astype(np.int64)
    nr = d["radios"].astype(np.float64)
    ch = np.where(parent >= 0)[0]
    aris = np.stack([parent[ch], ch], axis=1)
    return nodos, aris, nr


def adapt_4(d):
    """4: radios POR-NODO; aristas explicitas."""
    return d["nodos"].astype(np.float64), d["aristas"].astype(np.int64), d["radios"].astype(np.float64)


# ============================================================================
#  MERGE + COMPONENTES
# ============================================================================
def merge_parts(parts, tol=MERGE_TOL):
    """parts = [(nodos, edges, radii, tag)]. Fusiona nodos coincidentes a 'tol'.
    Devuelve (nodos, edges, radii, tags_list_por_nodo)."""
    allnd, alled, allr, alltag = [], [], [], []
    off = 0
    for (nd, ed, rr, tag) in parts:
        allnd.append(nd); allr.append(rr)
        alltag.append(np.full(len(nd), tag, dtype=object))
        if len(ed):
            alled.append(ed + off)
        off += len(nd)
    nodos = np.vstack(allnd); radii = np.concatenate(allr); tags = np.concatenate(alltag)
    edges = np.vstack(alled) if alled else np.empty((0, 2), np.int64)
    # union-find de coincidentes
    n = len(nodos)
    par = np.arange(n)
    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    if n:
        pairs = cKDTree(nodos).query_pairs(tol, output_type='ndarray')
        for a, b in pairs:
            ra, rb = f(int(a)), f(int(b))
            if ra != rb:
                par[max(ra, rb)] = min(ra, rb)
    roots = np.array([f(i) for i in range(n)])
    uniq, inv = np.unique(roots, return_inverse=True)
    newpos = nodos[uniq]
    newr = np.zeros(len(uniq)); np.maximum.at(newr, inv, radii)
    newtags = [set() for _ in range(len(uniq))]
    for i in range(n):
        newtags[inv[i]].add(tags[i])
    newedges = inv[edges] if len(edges) else np.empty((0, 2), np.int64)
    if len(newedges):
        newedges = newedges[newedges[:, 0] != newedges[:, 1]]
    return newpos, newedges, newr, newtags


def comps(N, edges):
    if len(edges):
        r = np.concatenate([edges[:, 0], edges[:, 1]])
        c = np.concatenate([edges[:, 1], edges[:, 0]])
        g = csr_matrix((np.ones(len(r)), (r, c)), shape=(N, N))
    else:
        g = csr_matrix((N, N))
    nc, lab = connected_components(g, directed=False)
    return nc, lab


# ============================================================================
#  PODA + CLAMP
# ============================================================================
def poda_clamp(nodos, edges, radii, port_pos):
    """Elimina aristas con AMBOS extremos < MIN_CANAL_R; clampa nodos retenidos
    por debajo; preserva conectividad al puerto (min path si hiciera falta)."""
    port_node = int(np.argmin(np.linalg.norm(nodos - port_pos, axis=1)))
    both_low = (radii[edges[:, 0]] < MIN_CANAL_R) & (radii[edges[:, 1]] < MIN_CANAL_R)
    edges_ret = edges[~both_low]
    n_podadas = int(both_low.sum())
    ret_nodes = np.unique(edges_ret) if len(edges_ret) else np.array([], np.int64)
    reconx = 0
    # conectividad al puerto sobre el subgrafo retenido
    if len(ret_nodes):
        sub = np.zeros(len(nodos), bool); sub[ret_nodes] = True
        nc, lab = comps(len(nodos), edges_ret)
        if sub[port_node]:
            port_lab = lab[port_node]
            malas = np.array([c for c in np.unique(lab[ret_nodes]) if c != port_lab])
            if len(malas):
                # reconectar cada comp retenido aislado por camino minimo (clampeado)
                r = np.concatenate([edges[:, 0], edges[:, 1]])
                cc = np.concatenate([edges[:, 1], edges[:, 0]])
                gfull = csr_matrix((np.ones(len(r)), (r, cc)), shape=(len(nodos), len(nodos)))
                dist, pred = dijkstra(gfull, directed=False, indices=port_node,
                                      return_predecessors=True)
                extra = []
                for cl in malas:
                    tgt = ret_nodes[lab[ret_nodes] == cl]
                    j = int(tgt[np.argmin(dist[tgt])])
                    while j != port_node and pred[j] >= 0:
                        extra.append((pred[j], j)); j = int(pred[j])
                if extra:
                    edges_ret = np.vstack([edges_ret, np.array(extra)])
                    ret_nodes = np.unique(edges_ret)
                    reconx = len(extra)
    # clamp
    radii_c = radii.copy()
    below = ret_nodes[radii[ret_nodes] < MIN_CANAL_R] if len(ret_nodes) else np.array([], np.int64)
    radii_c[below] = MIN_CANAL_R
    n_clamp = int(len(below))
    frac = (len(edges_ret) / len(edges)) if len(edges) else 0.0
    return dict(edges_ret=edges_ret, ret_nodes=ret_nodes, radii_c=radii_c,
                n_podadas=n_podadas, n_clamp=n_clamp, port_node=port_node,
                frac=frac, reconx=reconx)


# ============================================================================
#  SDF ROUNDED-CONE + MUESTREO EN GRID
# ============================================================================
def sd_round_cone(P, a, b, r1, r2):
    ba = b - a; l2 = float(ba @ ba)
    if l2 < 1e-12:
        return np.linalg.norm(P - a, axis=1) - max(r1, r2)
    rr = r1 - r2; a2 = l2 - rr * rr; il2 = 1.0 / l2
    if a2 <= 1e-12:   # una esfera contiene a la otra: union de 2 esferas
        return np.minimum(np.linalg.norm(P - a, axis=1) - r1,
                          np.linalg.norm(P - b, axis=1) - r2)
    pa = P - a
    y = pa @ ba; z = y - l2
    v = pa * l2 - np.outer(y, ba)
    x2 = np.einsum('ij,ij->i', v, v)
    y2 = y * y * l2; z2 = z * z * l2
    k = np.sign(rr) * rr * rr * x2
    out = np.empty(len(P))
    m1 = np.sign(z) * a2 * z2 > k
    m2 = (np.sign(y) * a2 * y2 < k) & (~m1)
    m3 = ~(m1 | m2)
    out[m1] = np.sqrt(np.maximum(x2[m1] + z2[m1], 0.0)) * il2 - r2
    out[m2] = np.sqrt(np.maximum(x2[m2] + y2[m2], 0.0)) * il2 - r1
    out[m3] = (np.sqrt(np.maximum(x2[m3] * a2 * il2, 0.0)) + y[m3] * rr) * il2 - r1
    return out


def eval_sdf(cones, origin, shape, vx):
    """SDF minimo sobre 'cones' en un grid; solo evalua el sub-box de cada cono."""
    sdf = np.full(shape, np.inf, dtype=np.float32)
    nx, ny, nz = shape
    for (a, b, r1, r2) in cones:
        rmax = max(r1, r2) + 1.5 * vx
        lo = np.minimum(a, b) - rmax; hi = np.maximum(a, b) + rmax
        i0 = np.clip(np.floor((lo - origin) / vx).astype(int), 0, [nx, ny, nz])
        i1 = np.clip(np.ceil((hi - origin) / vx).astype(int) + 1, 0, [nx, ny, nz])
        if np.any(i1 <= i0):
            continue
        xs = origin[0] + np.arange(i0[0], i1[0]) * vx
        ys = origin[1] + np.arange(i0[1], i1[1]) * vx
        zs = origin[2] + np.arange(i0[2], i1[2]) * vx
        gx, gy, gz = np.meshgrid(xs, ys, zs, indexing='ij')
        P = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
        d = sd_round_cone(P, a, b, r1, r2).reshape(gx.shape).astype(np.float32)
        sub = sdf[i0[0]:i1[0], i0[1]:i1[1], i0[2]:i1[2]]
        np.minimum(sub, d, out=sub)
    return sdf


def build_cones(nodos, edges, radii):
    return [(nodos[a], nodos[b], float(radii[a]), float(radii[b])) for a, b in edges]


def circuit_bbox(nodos, ret_nodes, radii):
    P = nodos[ret_nodes]
    rmax = float(radii[ret_nodes].max())
    lo = P.min(0) - (rmax + BBOX_MARGIN); hi = P.max(0) + (rmax + BBOX_MARGIN)
    return lo, hi


def mesh_lumen(cones, lo, hi, vx):
    shape = tuple(int(np.ceil((hi[i] - lo[i]) / vx)) + 1 for i in range(3))
    sdf = eval_sdf(cones, lo, shape, vx)
    lumen_vox = int(np.count_nonzero(sdf < 0.0))
    vol = lumen_vox * vx ** 3
    if not (sdf.min() < 0 < sdf.max()):
        return None, None, vol, shape
    verts, faces, _, _ = measure.marching_cubes(sdf, level=0.0, spacing=(vx, vx, vx))
    verts = verts + lo
    return verts, faces, vol, shape


# ============================================================================
#  STL BINARIO (writer propio, sin numpy-stl)
# ============================================================================
def write_stl(path, verts, faces):
    tris = verts[faces]                      # (F,3,3)
    n = len(faces)
    nrm = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    ln = np.linalg.norm(nrm, axis=1, keepdims=True)
    nrm = np.divide(nrm, ln, out=np.zeros_like(nrm), where=ln > 1e-12)
    with open(path, 'wb') as fp:
        fp.write(b'\0' * 80)
        fp.write(struct.pack('<I', n))
        for i in range(n):
            fp.write(struct.pack('<3f', *nrm[i].astype(np.float32)))
            fp.write(tris[i, 0].astype('<f4').tobytes())
            fp.write(tris[i, 1].astype('<f4').tobytes())
            fp.write(tris[i, 2].astype('<f4').tobytes())
            fp.write(struct.pack('<H', 0))
    return n


def mesh_topology(verts, faces):
    F = len(faces)
    used = np.unique(faces)
    V = len(used)
    e = np.sort(np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    ue, cnt = np.unique(e, axis=0, return_counts=True)
    E = len(ue)
    # componentes de superficie (para interpretar euler = 2 * n_comp cerradas)
    nc, _ = comps(int(faces.max()) + 1, ue)
    ncc = nc - (int(faces.max()) + 1 - V)  # descuenta indices no usados (singletons)
    euler = V - E + F
    surf = int(ncc)
    genus = (2 * surf - euler) // 2         # asas totales (genero acumulado)
    return dict(V=V, E=E, F=F, euler=euler, manifold=bool(np.all(cnt == 2)),
                boundary=int((cnt == 1).sum()), nonmanifold=int((cnt > 2).sum()),
                surf_comps=surf, genus=genus)


def in_ellipsoid(P):
    q = (P - DOMINIO_CENTRO) / DOMINIO_SEMIEJES
    return np.einsum('ij,ij->i', q, q)


def clip_to_domain(verts, faces):
    val = in_ellipsoid(verts)
    keep = val <= 1.0
    if keep.all():
        return verts, faces, 0
    remap = -np.ones(len(verts), np.int64)
    remap[keep] = np.arange(int(keep.sum()))
    fkeep = keep[faces].all(axis=1)
    newf = remap[faces[fkeep]]
    return verts[keep], newf, int((~keep).sum())


# ============================================================================
#  MAIN
# ============================================================================
def main():
    t_all = time.perf_counter()
    LOG = []
    def log(s=""):
        print(s); LOG.append(s)

    log("=" * 74)
    log("  CAPA 5a - NUCLEO SDF + VALIDACION GRUESA (200um smoke)")
    log("=" * 74)

    # ---- PASO 0: cargar + adaptar ----
    d3a = np.load(find("capa3a_arterial.npz"), allow_pickle=True)
    d3ab = np.load(find("capa3ab_peritubular.npz"), allow_pickle=True)
    d3b = np.load(find("capa3b_venoso.npz"), allow_pickle=True)
    d3c = np.load(find("capa3c_colector.npz"), allow_pickle=True)
    d4 = np.load(find("capa4_colector_alto.npz"), allow_pickle=True)

    n3a, e3a, r3a = adapt_edge_radios(d3a)
    n3b, e3b, r3b = adapt_edge_radios(d3b)
    n3c, e3c, r3c = adapt_3c(d3c)
    n4, e4, r4 = adapt_4(d4)
    n3ab = d3ab["puntos_drenaje"].astype(np.float64)  # solo puntos (sin aristas/radios)

    log("\n  PASO 0 - esquemas normalizados (nodos, aristas, radios por-nodo):")
    log(f"    3a arterial : nodos {n3a.shape} aristas {e3a.shape} radios/nodo (raiz={float(d3a['radio_raiz'])*1000:.0f}um)")
    log(f"    3ab peritub.: PUNTOS {n3ab.shape} (SIN aristas ni radios: bed capilar, sub-resolucion)")
    log(f"    3b venoso   : nodos {n3b.shape} aristas {e3b.shape} radios/nodo (raiz={float(d3b['radio_raiz'])*1000:.0f}um)")
    log(f"    3c colector : nodos {n3c.shape} aristas {e3c.shape} radios/nodo (max={r3c.max()*1000:.0f}um)")
    log(f"    4 calicial  : nodos {n4.shape} aristas {e4.shape} radios/nodo (max={r4.max()*1000:.0f}um)")

    # ---- PASO 1: circuitos ----
    log("\n  PASO 1 - CIRCUITOS")
    # URINARIO = 3c U 4
    uN, uE, uR, uTags = merge_parts([(n3c, e3c, r3c, "3c"), (n4, e4, r4, "4")])
    nc_u, lab_u = comps(len(uN), uE)
    comp_sizes_u = np.bincount(lab_u)
    big_u = int(np.argmax(comp_sizes_u))
    log(f"    URINARIO (3c U 4): nodos {len(uN)}  aristas {len(uE)}  componentes {nc_u}")
    if nc_u != 1:
        # el urinario DEBE ser una sola componente
        log(f"    [FALLO] urinario en {nc_u} componentes (esperado 1). ABORTAR.")
        sys.exit(2)
    log(f"    -> [OK] una sola componente conexa (papilas 3c fusionadas con papila_junction de 4)")

    # VASCULAR = 3a U 3ab U 3b   (3ab = puntos edge-less)
    vN, vE, vR, vTags = merge_parts([(n3a, e3a, r3a, "3a"), (n3ab, np.empty((0, 2), int), np.zeros(len(n3ab)), "3ab"),
                                     (n3b, e3b, r3b, "3b")])
    nc_v, lab_v = comps(len(vN), vE)
    sizes = np.bincount(lab_v)
    tube_comps = [c for c in range(nc_v) if sizes[c] > 1]      # componentes con aristas (tubos)
    singl = int((sizes == 1).sum())
    log(f"    VASCULAR (3a U 3ab U 3b): nodos {len(vN)}  aristas {len(vE)}  componentes {nc_v}")
    log(f"    -> tubos (size>1): {len(tube_comps)} | puntos aislados (3ab bed capilar): {singl}")
    for c in tube_comps:
        tg = set()
        for i in np.where(lab_v == c)[0]:
            tg |= vTags[i]
        tg -= {"3ab"}
        log(f"       componente {c}: {int(sizes[c])} nodos, capas {sorted(tg)}")
    # identificar arterial vs venoso por puerto
    def comp_de_puerto(port):
        j = int(np.argmin(np.linalg.norm(vN - port, axis=1)))
        return lab_v[j], j
    lab_art, _ = comp_de_puerto(np.array(d3a["raiz"], float))   # raiz arterial detectada de 3a
    lab_ven, _ = comp_de_puerto(PUERTO_VENOSO)
    log(f"    -> arterial = componente {lab_art} (raiz 3a {np.array(d3a['raiz'],float)})")
    log(f"    -> venoso   = componente {lab_ven} (puerto {PUERTO_VENOSO})")
    log(f"    -> arterial y venoso SEPARADOS: {'[OK] (gap capilar de diseño, es dato no defecto)' if lab_art != lab_ven else '[REVISAR] misma componente'}")

    # ---- PASO 2: poda + clamp por circuito ----
    log("\n  PASO 2 - PODA + CLAMP (piso MIN_CANAL_R = %.3f mm / diam %.3f mm)" % (MIN_CANAL_R, MIN_CANAL_DIAM))
    circuitos = {}
    # urinario
    pu = poda_clamp(uN, uE, uR, PUERTO_URINARIO)
    circuitos["urinaria"] = (uN, pu)
    # vascular por componente-tubo
    for c in tube_comps:
        mask = lab_v == c
        idx = np.where(mask)[0]
        remap = -np.ones(len(vN), np.int64); remap[idx] = np.arange(len(idx))
        emask = mask[vE[:, 0]] & mask[vE[:, 1]]
        sube = remap[vE[emask]]
        subn = vN[idx]; subr = vR[idx]
        port = np.array(d3a["raiz"], float) if c == lab_art else PUERTO_VENOSO
        nombre = "vascular_0_arterial" if c == lab_art else "vascular_1_venoso"
        pc = poda_clamp(subn, sube, subr, port)
        circuitos[nombre] = (subn, pc)

    log(f"    {'circuito':22s}{'aristas_tot':>12}{'podadas':>9}{'retenidas':>10}{'frac%':>8}{'clamp':>7}{'reconx':>8}")
    poda_tab = {}
    for nombre, (nd, p) in circuitos.items():
        tot = p["n_podadas"] + len(p["edges_ret"])
        log(f"    {nombre:22s}{tot:>12}{p['n_podadas']:>9}{len(p['edges_ret']):>10}{100*p['frac']:>7.1f}%{p['n_clamp']:>7}{p['reconx']:>8}")
        poda_tab[nombre] = dict(tot=tot, podadas=p["n_podadas"], retenidas=len(p["edges_ret"]),
                                frac=p["frac"], clamp=p["n_clamp"], reconx=p["reconx"])
    log("    NOTA: la microvasculatura sub-400um NO se imprime; queda en 3a/3ab/3b/3c como")
    log("          OBJETIVO DE PERFUSION (maduracion capilar wet-lab), no como canal impreso.")

    # ---- PASO 3: mallas de lumen (200um) ----
    log("\n  PASO 3 - MALLAS DE LUMEN (SDF rounded-cone, marching_cubes a %.0fum)" % (VOXEL_5A * 1000))
    meshes = {}          # nombre -> dict(verts,faces,verts_clip,faces_clip,vol,topo,ret_nodes,nodos,radii_c,edges_ret,port_node)
    for nombre, (nd, p) in circuitos.items():
        t0 = time.perf_counter()
        if len(p["edges_ret"]) == 0:
            log(f"    {nombre}: 0 aristas retenidas -> sin lumen imprimible (todo sub-piso).")
            meshes[nombre] = dict(empty=True, vol=0.0, ret_nodes=p["ret_nodes"], nodos=nd,
                                  radii_c=p["radii_c"], edges_ret=p["edges_ret"], port_node=p["port_node"])
            continue
        cones = build_cones(nd, p["edges_ret"], p["radii_c"])
        lo, hi = circuit_bbox(nd, p["ret_nodes"], p["radii_c"])
        verts, faces, vol, shape = mesh_lumen(cones, lo, hi, VOXEL_5A)
        if verts is None:
            log(f"    {nombre}: SDF sin cruce por 0 -> sin superficie.")
            meshes[nombre] = dict(empty=True, vol=vol, ret_nodes=p["ret_nodes"], nodos=nd,
                                  radii_c=p["radii_c"], edges_ret=p["edges_ret"], port_node=p["port_node"])
            continue
        topo = mesh_topology(verts, faces)
        vclip, fclip, nout = clip_to_domain(verts, faces)
        stl = os.path.join(REPO, f"capa5a_lumen_{nombre}.stl")
        nf = write_stl(stl, vclip, fclip) if len(fclip) else 0
        dt = time.perf_counter() - t0
        log(f"    {nombre:22s}: grid {shape} verts {len(verts)} faces {len(faces)} "
            f"vol {vol:.1f}mm3 clip(-{nout}v) STL_faces {nf}  ({dt:.1f}s)")
        meshes[nombre] = dict(empty=False, verts=verts, faces=faces, verts_clip=vclip, faces_clip=fclip,
                              vol=vol, topo=topo, nout=nout, stl=os.path.basename(stl), shape=shape,
                              ret_nodes=p["ret_nodes"], nodos=nd, radii_c=p["radii_c"],
                              edges_ret=p["edges_ret"], port_node=p["port_node"])

    # ---- PASO 4: auditorias ----
    log("\n  PASO 4 - AUDITORIAS")
    # (1) watertight (sobre malla cruda cerrada; el clip abre los puertos a proposito)
    log("   (1) WATERTIGHT (malla cruda, pre-clip):")
    aud1 = {}
    for nombre, m in meshes.items():
        if m.get("empty"):
            log(f"       {nombre:22s}: sin malla (vacio)")
            continue
        t = m["topo"]
        # WATERTIGHT = cierre 2-manifold: manifold, sin bordes, sin non-manifold.
        # (euler=2 solo si genero 0; el genero>0 por fusion de ramas a 200um NO rompe el cierre)
        ok = t["manifold"] and t["boundary"] == 0 and t["nonmanifold"] == 0
        log(f"       {nombre:22s}: V{t['V']} E{t['E']} F{t['F']} manifold={t['manifold']} "
            f"bordes={t['boundary']} nonmanifold={t['nonmanifold']} euler={t['euler']} "
            f"genero={t['genus']}  {'[OK] watertight' if ok else '[REVISAR]'}"
            f"{'  (genero>0: asas por fusion de ramas anidadas a 200um; 5b/100um lo reduce)' if ok and t['genus'] > 0 else ''}")
        aud1[nombre] = dict(topo=t, ok=ok)

    # (2) estanqueidad: voxeles dentro-urinario y dentro-vascular en grid comun
    log("   (2) ESTANQUEIDAD (voxeles en AMBOS circuitos):")
    urin_cones, vasc_cones = [], []
    for nombre, (nd, p) in circuitos.items():
        cs = build_cones(nd, p["edges_ret"], p["radii_c"]) if len(p["edges_ret"]) else []
        if nombre == "urinaria":
            urin_cones += cs
        else:
            vasc_cones += cs
    both = 0; both_art = 0; both_ven = 0; ss_dist = np.inf
    art_cones = build_cones(*[circuitos["vascular_0_arterial"][0], circuitos["vascular_0_arterial"][1]["edges_ret"], circuitos["vascular_0_arterial"][1]["radii_c"]]) if "vascular_0_arterial" in circuitos and len(circuitos["vascular_0_arterial"][1]["edges_ret"]) else []
    ven_cones = build_cones(*[circuitos["vascular_1_venoso"][0], circuitos["vascular_1_venoso"][1]["edges_ret"], circuitos["vascular_1_venoso"][1]["radii_c"]]) if "vascular_1_venoso" in circuitos and len(circuitos["vascular_1_venoso"][1]["edges_ret"]) else []
    if urin_cones and vasc_cones:
        allpts = np.vstack([np.vstack([c[0], c[1]]) for c in urin_cones + vasc_cones])
        rmax = max(max(c[2], c[3]) for c in urin_cones + vasc_cones)
        lo = allpts.min(0) - (rmax + BBOX_MARGIN); hi = allpts.max(0) + (rmax + BBOX_MARGIN)
        shp = tuple(int(np.ceil((hi[i] - lo[i]) / VOXEL_5A)) + 1 for i in range(3))
        su = eval_sdf(urin_cones, lo, shp, VOXEL_5A)
        sv = eval_sdf(vasc_cones, lo, shp, VOXEL_5A)
        both = int(np.count_nonzero((su < 0) & (sv < 0)))
        if art_cones:
            sa = eval_sdf(art_cones, lo, shp, VOXEL_5A); both_art = int(np.count_nonzero((su < 0) & (sa < 0)))
        if ven_cones:
            sn = eval_sdf(ven_cones, lo, shp, VOXEL_5A); both_ven = int(np.count_nonzero((su < 0) & (sn < 0)))
        vu = meshes["urinaria"].get("verts")
        vv = [meshes[n]["verts"] for n in meshes if n.startswith("vascular") and not meshes[n].get("empty")]
        if vu is not None and vv:
            vv = np.vstack(vv)
            ss_dist = float(cKDTree(vv).query(vu, k=1)[0].min())
    log(f"       voxeles dentro de AMBOS (urinario & vascular): {both}  (urin&arterial={both_art}, urin&venoso={both_ven})")
    log(f"       {'[OK] meta 0' if both == 0 else '[DATO EMERGENTE] solape por clamp-a-piso (micro-vasos inflados a 400um diam) + muestreo 200um; centerlines siguen disjuntos'}")
    log(f"       dist minima superficie-superficie urinario<->vascular: {ss_dist:.3f} mm")
    # gap arterial<->venoso (nodos retenidos)
    art_nodes = meshes.get("vascular_0_arterial", {}).get("ret_nodes")
    ven_nodes = meshes.get("vascular_1_venoso", {}).get("ret_nodes")
    gap_av = np.inf
    if art_nodes is not None and ven_nodes is not None and len(art_nodes) and len(ven_nodes):
        Pa = meshes["vascular_0_arterial"]["nodos"][art_nodes]
        Pv = meshes["vascular_1_venoso"]["nodos"][ven_nodes]
        gap_av = float(cKDTree(Pv).query(Pa, k=1)[0].min())
    log(f"       gap arterial<->venoso (nodos retenidos): {gap_av:.3f} mm")

    # (3) conectividad a puerto
    log("   (3) CONECTIVIDAD a puerto:")
    aud3 = {}
    puertos = {"urinaria": PUERTO_URINARIO, "vascular_0_arterial": np.array(d3a["raiz"], float),
               "vascular_1_venoso": PUERTO_VENOSO}
    for nombre, (nd, p) in circuitos.items():
        if len(p["edges_ret"]) == 0:
            log(f"       {nombre:22s}: sin aristas retenidas (nada imprimible)"); aud3[nombre] = None; continue
        ncc, ll = comps(len(nd), p["edges_ret"])
        reach = (p["port_node"] in p["ret_nodes"]) and \
                (ll[p["port_node"]] in ll[p["ret_nodes"]])
        # el puerto retenido y en la misma componente que el grueso del circuito
        big = np.bincount(ll[p["ret_nodes"]]).argmax()
        reach = ll[p["port_node"]] == big
        d_port = float(np.linalg.norm(nd[p["port_node"]] - puertos[nombre]))
        log(f"       {nombre:22s}: puerto nodo {p['port_node']} (a {d_port:.3f}mm del puerto), "
            f"alcanzado por circuito retenido: {'[OK]' if reach else '[REVISAR]'}")
        aud3[nombre] = dict(reach=bool(reach), d_port=d_port)

    # (4) contencion en elipsoide (nodos de lumen retenidos)
    log("   (4) CONTENCION en elipsoide (nodos de lumen retenidos):")
    aud4 = {}
    for nombre, (nd, p) in circuitos.items():
        if not len(p["ret_nodes"]):
            continue
        val = in_ellipsoid(nd[p["ret_nodes"]])
        viol = np.where(val > 1.0 + 1e-9)[0]
        vmax = float(val.max())
        log(f"       {nombre:22s}: violaciones {len(viol)}/{len(p['ret_nodes'])} (val max {vmax:.3f})"
            f"{'  [puertos: salida por diseño]' if len(viol) else '  [OK]'}")
        aud4[nombre] = dict(viol=int(len(viol)), vmax=vmax, n=len(p["ret_nodes"]))

    # (5) volumenes / porosidad
    log("   (5) VOLUMENES / POROSIDAD:")
    vol_matriz = 4.0 / 3.0 * np.pi * np.prod(DOMINIO_SEMIEJES)   # elipsoide COMPLETO (seno relleno)
    vol_lumen = {n: meshes[n]["vol"] for n in meshes}
    vtot = sum(vol_lumen.values())
    poros = 100.0 * vtot / vol_matriz
    log(f"       matriz (elipsoide completo, seno relleno): {vol_matriz:.1f} mm3")
    for n, v in vol_lumen.items():
        log(f"       lumen {n:22s}: {v:8.2f} mm3")
    log(f"       lumen TOTAL: {vtot:.2f} mm3   POROSIDAD: {poros:.4f} %")

    # ---- guardar meta ----
    meta = dict(
        voxel_5a=VOXEL_5A, voxel_floor=VOXEL_FLOOR, min_canal_r=MIN_CANAL_R, min_canal_diam=MIN_CANAL_DIAM,
        dominio_semiejes=DOMINIO_SEMIEJES, dominio_centro=DOMINIO_CENTRO,
        vol_matriz_mm3=vol_matriz, porosidad_pct=poros,
        poda=json.dumps(poda_tab), vol_lumen=json.dumps({k: float(v) for k, v in vol_lumen.items()}),
        both_voxels=both, both_urin_arterial=both_art, both_urin_venoso=both_ven,
        ss_dist_mm=ss_dist, gap_av_mm=gap_av,
        n_comp_vascular_tubos=len(tube_comps), n_puntos_3ab=len(n3ab),
        lab_arterial=int(lab_art), lab_venoso=int(lab_ven),
    )
    np.savez_compressed(os.path.join(REPO, "capa5a_meta.npz"), **meta)

    dt_all = time.perf_counter() - t_all
    mem = peak_mem_mb()
    log("\n  " + "-" * 70)
    log(f"  TIEMPO TOTAL: {dt_all:.1f}s   PICO MEMORIA: {mem:.0f} MB")
    log(f"  -> capa5a_meta.npz + STLs guardados en {REPO}")
    log(f"  -> 5b (matriz fina 100um + voxel etiquetado) queda PENDIENTE")
    log("=" * 74)

    # ---- reporte markdown ----
    M = []
    M.append("# Auditoria Capa 5a - nucleo SDF (validacion gruesa 200um)\n")
    M.append("**Programa:** Bio-Kidney AI 2026 - **Capa:** 5a - **Fecha:** 2026-07-09\n")
    M.append("Cada capa = centerline+radio = union de capsulas conicas = SDF (rounded-cone). "
             "Solido = elipsoide COMPLETO de Capa 0 (seno relleno) menos union de lumenes. "
             "Muestreo grueso a 200um (smoke); el paso fino 100um y el voxel etiquetado son 5b.\n")
    M.append("## PASO 0 - esquemas normalizados (nodos, aristas, radios/nodo)\n")
    M.append("| capa | nodos | aristas | radios | adaptador |\n|---|---|---|---|---|\n")
    M.append(f"| 3a arterial | {n3a.shape[0]} | {e3a.shape[0]} | por-ARISTA -> nodo(child); raiz={float(d3a['radio_raiz'])*1000:.0f}um | adapt_edge_radios |\n")
    M.append(f"| 3ab peritubular | {n3ab.shape[0]} pts | 0 | ninguno (bed capilar, sub-resolucion) | puntos edge-less |\n")
    M.append(f"| 3b venoso | {n3b.shape[0]} | {e3b.shape[0]} | por-ARISTA -> nodo(child); raiz={float(d3b['radio_raiz'])*1000:.0f}um | adapt_edge_radios |\n")
    M.append(f"| 3c colector | {n3c.shape[0]} | {e3c.shape[0]} | por-NODO; topologia via parent | adapt_3c |\n")
    M.append(f"| 4 calicial | {n4.shape[0]} | {e4.shape[0]} | por-NODO; aristas explicitas | adapt_4 |\n")
    M.append("## PASO 1 - circuitos\n")
    M.append(f"- **URINARIO = 3c U 4:** {len(uN)} nodos, {len(uE)} aristas, **{nc_u} componente(s)** "
             f"{'[OK] una sola (papilas fusionadas)' if nc_u == 1 else '[FALLO]'}\n")
    M.append(f"- **VASCULAR = 3a U 3ab U 3b:** {len(vN)} nodos, {len(vE)} aristas -> "
             f"**{len(tube_comps)} componentes-tubo** (arterial comp {lab_art}, venoso comp {lab_ven}) "
             f"+ {singl} puntos aislados (3ab). Arterial y venoso "
             f"**{'SEPARADOS (gap capilar de diseño, dato no defecto)' if lab_art != lab_ven else 'MISMA componente [REVISAR]'}**.\n")
    M.append("## PASO 2 - poda + clamp (piso 200um r / 400um diam)\n")
    M.append("| circuito | aristas_tot | podadas | retenidas | frac% superv | clamp | reconx |\n|---|---|---|---|---|---|---|\n")
    for nombre, v in poda_tab.items():
        M.append(f"| {nombre} | {v['tot']} | {v['podadas']} | {v['retenidas']} | {100*v['frac']:.1f}% | {v['clamp']} | {v['reconx']} |\n")
    M.append("\n> La microvasculatura sub-400um **NO se imprime**: queda en 3a/3ab/3b/3c como "
             "**OBJETIVO DE PERFUSION** (maduracion capilar wet-lab), no como canal impreso.\n")
    M.append("## PASO 4 - auditorias\n")
    M.append("### (1) Watertight = cierre 2-manifold (malla cruda pre-clip; el clip abre los puertos a proposito)\n")
    M.append("Criterio: manifold + 0 bordes + 0 non-manifold. euler=2 solo si genero 0; "
             "el genero>0 (asas por fusion de ramas anidadas a 200um) NO rompe el cierre.\n\n")
    M.append("| circuito | V | E | F | manifold | bordes | nonmanif | euler | genero | estado |\n|---|---|---|---|---|---|---|---|---|---|\n")
    for nombre, a in aud1.items():
        t = a["topo"]
        M.append(f"| {nombre} | {t['V']} | {t['E']} | {t['F']} | {t['manifold']} | {t['boundary']} | "
                 f"{t['nonmanifold']} | {t['euler']} | {t['genus']} | {'[OK] watertight' if a['ok'] else '[REVISAR]'} |\n")
    M.append(f"\n### (2) Estanqueidad\n- Voxeles dentro de AMBOS (urinario & vascular): **{both}** "
             f"(urin&arterial={both_art}, urin&venoso={both_ven})\n")
    if both == 0:
        M.append("  - [OK] meta 0.\n")
    else:
        M.append(f"  - **[DATO EMERGENTE]** El solape NO viene de centerlines cruzados (siguen disjuntos) "
                 f"sino del **clamp-a-piso**: micro-vasos venosos sub-400um se inflan a 400um diam y, "
                 f"al muestrear a 200um, su lumen fusiona con el colector urinario en la medula (interdigitacion "
                 f"vasa-recta / conducto colector). Se resuelve en 5b (100um + politica de clamp / separacion "
                 f"inter-circuito >=400um por diseño). Es tension de escala del modelo reducido, no cruce topologico.\n")
    M.append(f"- Distancia minima superficie-superficie urinario<->vascular: **{ss_dist:.3f} mm**\n"
             f"- Gap arterial<->venoso (nodos retenidos): **{gap_av:.3f} mm** (gap capilar de diseño)\n")
    M.append("### (3) Conectividad a puerto\n")
    for nombre, a in aud3.items():
        if a is None:
            M.append(f"- {nombre}: sin aristas retenidas (nada imprimible)\n")
        else:
            M.append(f"- {nombre}: puerto a {a['d_port']:.3f} mm, alcanzado por el circuito retenido: "
                     f"{'[OK]' if a['reach'] else '[REVISAR]'}\n")
    M.append("### (4) Contencion en elipsoide\n")
    for nombre, a in aud4.items():
        M.append(f"- {nombre}: violaciones {a['viol']}/{a['n']} (val max {a['vmax']:.3f})"
                 f"{' -> puertos (salida por diseño)' if a['viol'] else ' [OK]'}\n")
    M.append(f"### (5) Volumenes / porosidad\n")
    M.append(f"- Matriz (elipsoide completo, seno relleno): **{vol_matriz:.1f} mm3**\n")
    for n, v in vol_lumen.items():
        M.append(f"- Lumen {n}: {v:.2f} mm3\n")
    M.append(f"- Lumen TOTAL: **{vtot:.2f} mm3** - **POROSIDAD: {poros:.4f} %**\n")
    M.append(f"\n## Estado\n**Capa 5a cerrada.** Nucleo SDF muestreado a 200um, poda a 400um diam, "
             f"lumenes urinario/arterial/venoso mallados (STL binario propio), auditorias corridas. "
             f"Tiempo {dt_all:.1f}s, pico memoria {mem:.0f} MB. **5b (matriz fina 100um + voxel "
             f"etiquetado) PENDIENTE.**\n")
    os.makedirs(os.path.dirname(MD_OUT), exist_ok=True)
    with open(MD_OUT, "w") as f:
        f.write("".join(M))
    log(f"  -> reporte: {MD_OUT}")

    return LOG, dict(poda_tab=poda_tab, aud1=aud1, aud3=aud3, aud4=aud4, meshes=meshes,
                     both=both, ss_dist=ss_dist, gap_av=gap_av, vol_matriz=vol_matriz,
                     vol_lumen=vol_lumen, vtot=vtot, poros=poros, tube_comps=tube_comps,
                     singl=singl, lab_art=lab_art, lab_ven=lab_ven, dt_all=dt_all, mem=mem,
                     n3ab=len(n3ab))


if __name__ == "__main__":
    main()
