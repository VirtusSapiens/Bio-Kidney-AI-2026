#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa3a_arterial.py  --  Bio-Kidney AI 2026
==========================================
Capa 3a del gemelo digital: ARBOL ARTERIAL renal.

Genera el arbol arterial por SPACE COLONIZATION (Runions et al. 2005) desde el
hilio (entrada de la arteria renal) hacia los 1300 glomerulos de la Capa 1, que
actuan como PUNTOS DE ATRACCION. En cada bifurcacion se respeta la LEY DE MURRAY
(r^3 = sum r_hijas^3) para asignar radios coherentes de la raiz a las hojas.

La salida es un GRAFO (nodos + aristas con radio), NO voxeles. Las hojas del
arbol son las arteriolas aferentes terminales: se usan como capilares-proxy para
medir la cobertura del tejido y compararla con la linea base de la Capa 2.

Solo NumPy / SciPy. SIN bpy. Entorno: env_biokidney.

Entrada : capa0_dominio.npz, capa1_nefronas.npz, capa2_demanda.npz (cobertura)
Salida  : capa3a_arterial.npz
"""

import os
import sys
import numpy as np
from scipy.spatial import cKDTree

# ============================================================================
#  PARAMETROS  (editar aqui)
# ============================================================================
# RAIZ se toma del hilio de la Capa 0 (no se fija aqui; ver main()).
DIST_INFLUENCIA = 8.0     # mm  radio en que un atractor influye en el nodo mas cercano
DIST_MATAR      = 1.2     # mm  un atractor a esta distancia de un nodo se da por alcanzado
PASO_CRECIMIENTO = 0.6    # mm  longitud de cada segmento nuevo por iteracion
MAX_ITERACIONES = 2000    # tope de seguridad
MURRAY_EXP      = 3.0     # exponente de la ley de Murray: r^3 = sum(r_hijas^3)
RADIO_TERMINAL  = 0.012   # mm  ~12 micras, arteriola aferente terminal

SEED = 2026

IN_DOM = "capa0_dominio.npz"
IN_NEF = "capa1_nefronas.npz"
IN_DEM = "capa2_demanda.npz"
OUT_NPZ = "capa3a_arterial.npz"


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


def get_cobertura():
    """Devuelve la funcion cobertura() de la Capa 2.

    Intenta importarla del modulo capa2_demanda (misma carpeta). Si no se puede,
    re-implementa la metrica identica: % de puntos de dominio a < umbral del
    capilar mas cercano, distancia media/maxima y deficit ponderado por demanda.
    """
    here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() \
        else os.getcwd()
    if here not in sys.path:
        sys.path.insert(0, here)
    try:
        from capa2_demanda import cobertura as _cob  # noqa: import local
        print("  cobertura(): importada de capa2_demanda.py")
        return _cob
    except Exception as e:                            # re-implementacion de respaldo
        print(f"  cobertura(): re-implementada (no se pudo importar capa2_demanda: {e})")

        def cobertura(puntos_capilares, coords, demanda, umbral_mm=0.150):
            n = len(coords)
            cap = np.asarray(puntos_capilares, dtype=np.float64).reshape(-1, 3) \
                if puntos_capilares is not None else np.empty((0, 3))
            if len(cap) == 0:
                return dict(pct_cubierto=0.0, dist_media=np.inf, dist_max=np.inf,
                            deficit=float(demanda.sum()), deficit_frac=1.0,
                            dist=np.full(n, np.inf))
            tree = cKDTree(cap)
            dist, _ = tree.query(coords, k=1)
            cubierto = dist < umbral_mm
            dem_total = float(demanda.sum())
            deficit = float(demanda[~cubierto].sum())
            return dict(
                pct_cubierto=100.0 * np.count_nonzero(cubierto) / n,
                dist_media=float(dist.mean()),
                dist_max=float(dist.max()),
                deficit=deficit,
                deficit_frac=(deficit / dem_total) if dem_total > 0 else 0.0,
                dist=dist,
            )
        return cobertura


# ============================================================================
#  SPACE COLONIZATION  (genera la topologia: nodos + parentesco)
# ============================================================================
def space_colonization(raiz, atractores):
    """Crece el arbol desde 'raiz' (3,) hacia 'atractores' (A,3).

    Algoritmo clasico de Runions et al. 2005:
      a) cada atractor vota por el NODO mas cercano (cKDTree),
      b) solo votan los atractores dentro de DIST_INFLUENCIA,
      c) cada nodo votado crece PASO_CRECIMIENTO en la direccion media
         (normalizada) hacia sus atractores,
      d) se eliminan los atractores a < DIST_MATAR de cualquier nodo.

    Como cada nodo nuevo nace de uno ya existente, parent[hijo] < hijo siempre:
    el grafo es un arbol enraizado en el indice 0.

    Devuelve (nodos (M,3), parent (M,) int  con parent[0]=-1, n_iter).
    """
    nodos = [np.asarray(raiz, dtype=np.float64)]
    parent = [-1]
    atrac = np.asarray(atractores, dtype=np.float64).copy()

    n_iter = 0
    for it in range(MAX_ITERACIONES):
        n_iter = it + 1
        if len(atrac) == 0:
            break

        nodos_arr = np.asarray(nodos)
        tree = cKDTree(nodos_arr)
        dist, idx = tree.query(atrac, k=1)            # nodo mas cercano por atractor

        activo = dist <= DIST_INFLUENCIA              # solo atractores con influencia
        if not np.any(activo):
            # BOOTSTRAP del tronco: el hilio esta a >DIST_INFLUENCIA del glomerulo
            # mas cercano (la arteria renal entra "lejos" del parenquima). Crecemos
            # un segmento recto desde el nodo mas cercano hacia el atractor mas
            # cercano hasta que la copa entre en rango y arranque el algoritmo.
            a = int(np.argmin(dist))
            nd = int(idx[a])
            d = atrac[a] - nodos_arr[nd]
            dn = np.linalg.norm(d)
            if dn < 1e-9:
                break
            nodos.append(nodos_arr[nd] + PASO_CRECIMIENTO * (d / dn))
            parent.append(nd)
            nodos_arr = np.asarray(nodos)
            dk, _ = cKDTree(nodos_arr).query(atrac, k=1)
            atrac = atrac[dk > DIST_MATAR]
            continue

        near = idx[activo]                            # nodo que cada atractor activo atrae
        vecs = atrac[activo] - nodos_arr[near]
        norm = np.linalg.norm(vecs, axis=1, keepdims=True)
        norm[norm == 0] = 1.0
        vecs /= norm                                  # direcciones unitarias hacia atractores

        # suma de direcciones por nodo votado
        acc = np.zeros_like(nodos_arr)
        np.add.at(acc, near, vecs)
        uniq = np.unique(near)
        g = acc[uniq]
        gn = np.linalg.norm(g, axis=1)
        keep = gn > 1e-9                              # descarta direcciones que se cancelan
        uniq = uniq[keep]
        g = g[keep] / gn[keep, None]

        if len(uniq) == 0:
            break

        nuevos = nodos_arr[uniq] + PASO_CRECIMIENTO * g
        base = len(nodos)
        for k in range(len(uniq)):
            nodos.append(nuevos[k])
            parent.append(int(uniq[k]))

        # d) matar atractores alcanzados por CUALQUIER nodo (incluidos los nuevos)
        nodos_arr = np.asarray(nodos)
        tree2 = cKDTree(nodos_arr)
        dk, _ = tree2.query(atrac, k=1)
        atrac = atrac[dk > DIST_MATAR]

    return np.asarray(nodos), np.asarray(parent, dtype=np.int64), len(atrac), n_iter


# ============================================================================
#  RADIOS POR LEY DE MURRAY  (post-proceso, de hojas a raiz)
# ============================================================================
def asignar_radios_murray(parent):
    """Asigna a cada nodo el radio del segmento que lo ALIMENTA (su arista al padre).

    - Hoja (sin hijos)        -> RADIO_TERMINAL.
    - Nodo con hijos          -> (sum r_hijo^MURRAY_EXP)^(1/MURRAY_EXP).
      (con un solo hijo el radio se conserva; con varios, bifurcacion de Murray.)

    Como parent[i] < i, recorrer i de M-1 a 0 procesa los hijos antes que el padre.

    Devuelve r_in (M,) radio entrante por nodo; r_in[raiz] es el radio de la raiz.
    """
    M = len(parent)
    hijos = [[] for _ in range(M)]
    for i in range(M):
        p = parent[i]
        if p >= 0:
            hijos[p].append(i)

    r_in = np.zeros(M, dtype=np.float64)
    for i in range(M - 1, -1, -1):
        h = hijos[i]
        if not h:
            r_in[i] = RADIO_TERMINAL
        else:
            r_in[i] = sum(r_in[c] ** MURRAY_EXP for c in h) ** (1.0 / MURRAY_EXP)
    return r_in, hijos


# ============================================================================
#  MAIN
# ============================================================================
def main():
    np.random.seed(SEED)

    dom_path = find_npz(IN_DOM)
    nef_path = find_npz(IN_NEF)
    dem_path = find_npz(IN_DEM)

    print("=" * 70)
    print("  CAPA 3a - ARBOL ARTERIAL (Space Colonization + Ley de Murray)")
    print("=" * 70)
    print("  Entrada dominio :", dom_path)
    print("  Entrada nefronas:", nef_path)
    print("  Entrada demanda :", dem_path)
    cobertura = get_cobertura()
    print("-" * 70)
    print(f"  DIST_INFLUENCIA = {DIST_INFLUENCIA} mm   DIST_MATAR = {DIST_MATAR} mm")
    print(f"  PASO_CRECIMIENTO = {PASO_CRECIMIENTO} mm   MAX_ITERACIONES = {MAX_ITERACIONES}")
    print(f"  MURRAY_EXP = {MURRAY_EXP}   RADIO_TERMINAL = {RADIO_TERMINAL} mm")
    print("-" * 70)

    # --- cargar datos -------------------------------------------------------
    d0 = np.load(dom_path, allow_pickle=False)
    coords = d0["coords"].astype(np.float64)
    region = d0["region_label"]
    hilio = d0["hilio"].astype(np.float64)

    d1 = np.load(nef_path, allow_pickle=False)
    glomerulos = d1["glomerulos"].astype(np.float64)

    d2 = np.load(dem_path, allow_pickle=False)
    demanda = d2["demanda"].astype(np.float64)
    umbral_mm = float(d2["umbral_difusion_mm"])
    base_cov_pct = float(d2["cov_proxy_pct"])

    raiz = hilio.copy()
    n_glom = len(glomerulos)
    print(f"  RAIZ (hilio): {raiz}")
    print(f"  PUNTOS_ATRACCION (glomerulos): {n_glom}")
    print(f"  Puntos de dominio: {len(coords)}   umbral cobertura: {umbral_mm} mm")
    print("-" * 70)

    # --- 1) topologia por space colonization --------------------------------
    print("  Creciendo arbol por space colonization ...")
    nodos, parent, n_no_alcanzados, n_iter = space_colonization(raiz, glomerulos)
    n_alcanzados = n_glom - n_no_alcanzados

    # --- 2) radios por Murray -----------------------------------------------
    r_in, hijos = asignar_radios_murray(parent)

    # aristas (padre -> hijo) y su radio = radio entrante del hijo
    edge_list = [(parent[i], i) for i in range(len(nodos)) if parent[i] >= 0]
    aristas = np.asarray(edge_list, dtype=np.int64)
    radios = np.asarray([r_in[i] for (_, i) in edge_list], dtype=np.float64)

    # nodos terminales (hojas, sin hijos) -> arteriolas aferentes
    terminales = np.asarray([i for i in range(len(nodos)) if not hijos[i]], dtype=np.int64)
    term_pos = nodos[terminales]
    radio_raiz = float(r_in[0])

    # --- 3) cobertura usando terminales como capilares-proxy ----------------
    cov = cobertura(term_pos, coords, demanda, umbral_mm)

    # --- guardar ------------------------------------------------------------
    out_path = os.path.join(os.path.dirname(dom_path), OUT_NPZ)
    np.savez_compressed(
        out_path,
        nodos=nodos.astype(np.float32),
        aristas=aristas.astype(np.int32),
        radios=radios.astype(np.float32),
        terminales=terminales.astype(np.int32),
        # metadata: parametros
        raiz=raiz.astype(np.float64),
        dist_influencia=np.float64(DIST_INFLUENCIA),
        dist_matar=np.float64(DIST_MATAR),
        paso_crecimiento=np.float64(PASO_CRECIMIENTO),
        max_iteraciones=np.int32(MAX_ITERACIONES),
        murray_exp=np.float64(MURRAY_EXP),
        radio_terminal=np.float64(RADIO_TERMINAL),
        radio_raiz=np.float64(radio_raiz),
        seed=np.int32(SEED),
        # metadata: resultados
        n_glomerulos=np.int32(n_glom),
        n_glomerulos_alcanzados=np.int32(n_alcanzados),
        n_glomerulos_no_alcanzados=np.int32(n_no_alcanzados),
        n_nodos=np.int32(len(nodos)),
        n_aristas=np.int32(len(aristas)),
        n_iteraciones=np.int32(n_iter),
        cov_pct=np.float64(cov["pct_cubierto"]),
        cov_dist_media=np.float64(cov["dist_media"]),
        cov_dist_max=np.float64(cov["dist_max"]),
        cov_deficit=np.float64(cov["deficit"]),
        cov_deficit_frac=np.float64(cov["deficit_frac"]),
    )

    # ========================================================================
    #  VERIFICACION
    # ========================================================================
    print("\n  VERIFICACION")
    print("-" * 70)
    print(f"  Iteraciones ejecutadas : {n_iter} / {MAX_ITERACIONES}")
    print(f"  Nodos generados        : {len(nodos)}")
    print(f"  Aristas generadas      : {len(aristas)}  (arbol: aristas = nodos-1 = {len(nodos)-1})")
    print(f"  Nodos terminales (hojas): {len(terminales)}")

    pct_alcanzados = 100.0 * n_alcanzados / n_glom
    print(f"\n  Glomerulos alcanzados  : {n_alcanzados} / {n_glom}  ({pct_alcanzados:.2f} %)")
    print(f"  Glomerulos NO alcanzados: {n_no_alcanzados}  "
          f"({'OK >95%' if pct_alcanzados > 95 else 'REVISAR <95%'})")

    print(f"\n  Radios (Murray):")
    print(f"     radio en la RAIZ     : {radio_raiz:.4f} mm")
    print(f"     radio TERMINAL (hoja): {RADIO_TERMINAL:.4f} mm")
    print(f"     rango de aristas     : [{radios.min():.4f}, {radios.max():.4f}] mm")
    print(f"     -> arteria renal real ~2-3 mm; orden de magnitud: "
          f"{'razonable' if 0.5 <= radio_raiz <= 6.0 else 'REVISAR'}")
    print(f"     -> la raiz es el radio maximo: "
          f"{'OK' if abs(radio_raiz - radios.max()) < 1e-9 else 'REVISAR'}")

    # --- cobertura global ---------------------------------------------------
    print(f"\n  COBERTURA  (terminales del arbol como capilares-proxy)")
    print("-" * 70)
    print(f"     % tejido cubierto (< {umbral_mm} mm): {cov['pct_cubierto']:.3f} %")
    print(f"     distancia media al capilar : {cov['dist_media']:.3f} mm")
    print(f"     distancia maxima (peor caso): {cov['dist_max']:.3f} mm")
    print(f"     deficit (demanda sin servir): {cov['deficit']:.1f}  "
          f"({100.0 * cov['deficit_frac']:.2f} % de la demanda total)")
    delta = cov["pct_cubierto"] - base_cov_pct
    print(f"     LINEA BASE Capa 2 (glomerulos): {base_cov_pct:.3f} %")
    print(f"     -> mejora: {delta:+.3f} puntos  (x{cov['pct_cubierto']/base_cov_pct:.1f} "
          f"sobre la base)" if base_cov_pct > 0 else "")

    # --- cobertura desglosada por region ------------------------------------
    print(f"\n  COBERTURA POR REGION  (cortex vs medula)")
    print("-" * 70)
    dist = cov["dist"]
    is_cortex = region == "cortex"
    is_medula = ~is_cortex                       # medulla + piramides (columnas de Bertin incl.)
    is_pir = np.char.startswith(region, "piramide")
    for nombre, mask in (("cortex", is_cortex),
                         ("medula (medulla+piramides)", is_medula),
                         ("  piramides", is_pir)):
        nm = int(np.count_nonzero(mask))
        if nm == 0:
            continue
        dsub = dist[mask]
        pct = 100.0 * np.count_nonzero(dsub < umbral_mm) / nm
        print(f"     {nombre:<28s}: n={nm:>7d}  cubierto={pct:6.3f} %  "
              f"dist_media={dsub.mean():.3f} mm  peor={dsub.max():.3f} mm")
    print("     (deficit medular esperado: el arbol arterial llega a glomerulos")
    print("      corticales; la medula se irriga via vasa recta, Capa 3 posterior.)")

    print("\n  -> guardado en:", out_path)
    print("=" * 70)


if __name__ == "__main__":
    main()
