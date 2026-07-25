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

# --- topologia (corrige los "mega-hubs" de la v1) -----------------------------
# DEFINITIVO = 2 (arbol BINARIO). Validado por comparacion empirica frente a la
# variante trinaria (GRADO_MAX_HIJOS=3): el binario IGUALA el 100% de glomerulos
# alcanzados del trinario, pero con la mitad de terminales, es mas fisiologico
# (un vaso real se bifurca en 2) y mas imprimible. No hay coste de cobertura, asi
# que el binario es la topologia canonica de la Capa 3a. La trinaria se conserva
# como respaldo (capa3a_arterial_v2_trinario.npz). Subir a 3 reactiva trifurcaciones.
GRADO_MAX_HIJOS = 2        # max. hijos por nodo en TODA la generacion (bifurcacion)
# --- rendimiento del cKDTree (evita reconstruirlo en cada iteracion) ----------
REBUILD_CADA_K  = 10       # reconstruye el arbol cada K iteraciones ...
REBUILD_NUEVOS  = 500      # ... o antes si se acumulan tantos nodos nuevos

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

    Algoritmo de Runions et al. 2005 con dos correcciones respecto a la v1:

    TOPOLOGIA (corrige los "mega-hubs"):
      Ningun nodo puede tener mas de GRADO_MAX_HIJOS hijos en TODA la generacion
      (se lleva un contador por nodo). Al alcanzar el tope el nodo queda CERRADO
      y sale del conjunto consultable: los atractores que lo tenian como nodo mas
      cercano se reasignan automaticamente al NODO ABIERTO mas cercano -- que
      suele ser uno de sus hijos --, forzando que el crecimiento continue DESDE
      las ramas hijas hacia abajo (ramificacion recursiva real). El cKDTree se
      construye SOLO sobre nodos abiertos, de modo que esa reasignacion es
      implicita.

    RENDIMIENTO (evita reconstruir el cKDTree de ~90k nodos cada iteracion):
      El arbol de nodos abiertos se reconstruye solo cada REBUILD_CADA_K
      iteraciones, o antes si se acumulan REBUILD_NUEVOS nodos nuevos. Entre
      reconstrucciones se consulta el arbol existente (ligeramente desactualizado
      es aceptable: los atractores afectados esperan, como mucho, a la proxima
      reconstruccion). Los atractores se "matan" solo contra los nodos NUEVOS de
      cada iteracion, que es equivalente a hacerlo contra todos (un atractor solo
      puede quedar a <DIST_MATAR por el avance de una rama recien creada).

    Pasos por iteracion (Runions):
      a) cada atractor vota por el NODO ABIERTO mas cercano (cKDTree),
      b) solo votan los atractores dentro de DIST_INFLUENCIA,
      c) cada nodo votado crece PASO_CRECIMIENTO en la direccion MEDIA
         (normalizada) hacia sus atractores  -> un hijo por iteracion,
      d) se eliminan los atractores a < DIST_MATAR de cualquier nodo nuevo.

    Como cada nodo nuevo nace de uno ya existente, parent[hijo] < hijo siempre:
    el grafo es un arbol enraizado en el indice 0.

    Devuelve (nodos (M,3), parent (M,) int  con parent[0]=-1, n_no_alcanzados,
              n_iter).
    """
    atrac = np.asarray(atractores, dtype=np.float64).copy()

    # --- buffers de nodos (crecimiento amortizado, evita reconstruir listas) ---
    cap = 1024
    pos_buf = np.empty((cap, 3), dtype=np.float64)
    parent_buf = np.empty(cap, dtype=np.int64)
    nhijos_buf = np.zeros(cap, dtype=np.int64)
    abierto_buf = np.zeros(cap, dtype=bool)
    pos_buf[0] = np.asarray(raiz, dtype=np.float64)
    parent_buf[0] = -1
    abierto_buf[0] = True
    M = 1  # n.º de nodos vivos en los buffers

    # --- estado del cKDTree incremental ---
    state = {"tree": None, "ids": None, "iters": 0, "nuevos": 0}

    def _ensure(n):
        nonlocal cap, pos_buf, parent_buf, nhijos_buf, abierto_buf
        if n <= cap:
            return
        while cap < n:
            cap *= 2
        pos_buf = np.resize(pos_buf, (cap, 3))
        parent_buf = np.resize(parent_buf, cap)
        nhijos_buf = np.resize(nhijos_buf, cap)
        abierto_buf = np.resize(abierto_buf, cap)

    def agregar_hijo(p, pos):
        """Anade un nodo hijo de 'p'; cierra 'p' si alcanza GRADO_MAX_HIJOS."""
        nonlocal M
        _ensure(M + 1)
        pos_buf[M] = pos
        parent_buf[M] = p
        nhijos_buf[M] = 0
        abierto_buf[M] = True
        M += 1
        nhijos_buf[p] += 1
        if nhijos_buf[p] >= GRADO_MAX_HIJOS:
            abierto_buf[p] = False          # nodo CERRADO: ya no emite mas ramas
        state["nuevos"] += 1

    def reconstruir():
        ids = np.nonzero(abierto_buf[:M])[0]            # solo nodos ABIERTOS
        state["tree"] = cKDTree(pos_buf[ids])
        state["ids"] = ids
        state["iters"] = 0
        state["nuevos"] = 0

    def matar(atrac, nuevos_pos):
        """Elimina atractores a < DIST_MATAR de algun nodo recien creado."""
        if len(atrac) == 0 or len(nuevos_pos) == 0:
            return atrac
        dk, _ = cKDTree(np.asarray(nuevos_pos)).query(atrac, k=1)
        return atrac[dk > DIST_MATAR]

    # matar atractores que ya nacen pegados a la raiz (normalmente ninguno)
    atrac = matar(atrac, [pos_buf[0]])
    reconstruir()

    n_iter = 0
    for it in range(MAX_ITERACIONES):
        n_iter = it + 1
        if len(atrac) == 0:
            break

        # reconstruir el arbol solo cada K iteraciones o si hay muchos nodos nuevos
        if state["iters"] >= REBUILD_CADA_K or state["nuevos"] >= REBUILD_NUEVOS:
            reconstruir()
        state["iters"] += 1

        dist, rows = state["tree"].query(atrac, k=1)    # nodo abierto mas cercano
        near = state["ids"][rows]                        # indice real de nodo

        en_rango = dist <= DIST_INFLUENCIA
        if not np.any(en_rango):
            # BOOTSTRAP del tronco: el hilio esta a >DIST_INFLUENCIA del glomerulo
            # mas cercano (la arteria renal entra "lejos" del parenquima). Crecemos
            # un segmento recto desde el nodo abierto mas cercano hacia el atractor
            # mas cercano hasta que la copa entre en rango y arranque el algoritmo.
            a = int(np.argmin(dist))
            nd = int(near[a])
            d = atrac[a] - pos_buf[nd]
            dn = np.linalg.norm(d)
            if dn < 1e-9:
                break
            nuevo = pos_buf[nd] + PASO_CRECIMIENTO * (d / dn)
            agregar_hijo(nd, nuevo)
            atrac = matar(atrac, [nuevo])
            reconstruir()                                # el tronco es corto: lo incluimos ya
            continue

        # un atractor cuyo nodo mas cercano se CERRO durante la ventana se descarta
        # esta iteracion; se reasignara a un nodo abierto tras la proxima reconstruccion
        activo = en_rango & abierto_buf[near]
        if not np.any(activo):
            reconstruir()                                # reasigna a nodos abiertos ya
            continue

        near_act = near[activo]
        vecs = atrac[activo] - pos_buf[near_act]
        norm = np.linalg.norm(vecs, axis=1, keepdims=True)
        norm[norm == 0] = 1.0
        vecs /= norm                                     # direcciones unitarias a atractores

        # direccion MEDIA por nodo votado (logica de Runions intacta)
        uniq, inv = np.unique(near_act, return_inverse=True)
        acc = np.zeros((len(uniq), 3), dtype=np.float64)
        np.add.at(acc, inv, vecs)
        gn = np.linalg.norm(acc, axis=1)
        keep = gn > 1e-9                                 # descarta direcciones que se cancelan
        uniq = uniq[keep]
        g = acc[keep] / gn[keep, None]
        if len(uniq) == 0:
            reconstruir()
            continue

        nuevos = pos_buf[uniq] + PASO_CRECIMIENTO * g     # un hijo por nodo votado
        nuevos_pos = []
        for k in range(len(uniq)):
            agregar_hijo(int(uniq[k]), nuevos[k])
            nuevos_pos.append(nuevos[k])

        # d) matar atractores alcanzados por los nodos NUEVOS
        atrac = matar(atrac, nuevos_pos)

    return (pos_buf[:M].copy(), parent_buf[:M].copy(), len(atrac), n_iter)


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
    print("-" * 70)
    print(f"  DIST_INFLUENCIA = {DIST_INFLUENCIA} mm   DIST_MATAR = {DIST_MATAR} mm")
    print(f"  PASO_CRECIMIENTO = {PASO_CRECIMIENTO} mm   MAX_ITERACIONES = {MAX_ITERACIONES}")
    print(f"  MURRAY_EXP = {MURRAY_EXP}   RADIO_TERMINAL = {RADIO_TERMINAL} mm")
    print("-" * 70)

    # --- cargar datos -------------------------------------------------------
    d0 = np.load(dom_path, allow_pickle=False)
    coords = d0["coords"].astype(np.float64)
    hilio = d0["hilio"].astype(np.float64)

    d1 = np.load(nef_path, allow_pickle=False)
    glomerulos = d1["glomerulos"].astype(np.float64)

    d2 = np.load(dem_path, allow_pickle=False)
    umbral_mm = float(d2["umbral_difusion_mm"])   # solo informativo en consola

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
    radio_raiz = float(r_in[0])

    # NOTA: la cobertura tisular (terminales como capilares-proxy) NO se calcula:
    # ya se auditó que no aplica al arbol ARTERIAL (la perfusion es trabajo de la
    # capa capilar). La metrica de exito aqui es el % de glomerulos ALCANZADOS.

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
        grado_max_hijos=np.int32(GRADO_MAX_HIJOS),
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

    # --- DISTRIBUCION DE GRADOS DE SALIDA (CRITICO: sin mega-hubs) -----------
    grados = np.array([len(hijos[i]) for i in range(len(nodos))], dtype=np.int64)
    grado_max = int(grados.max())
    print(f"\n  DISTRIBUCION DE GRADOS DE SALIDA (n.º de hijos por nodo)")
    print("-" * 70)
    for g in range(0, grado_max + 1):
        etiqueta = {0: "terminales", 1: "pasantes", 2: "bifurcaciones",
                    3: "trifurcaciones"}.get(g, "")
        print(f"     {g} hijos: {int(np.count_nonzero(grados == g)):>7d}   {etiqueta}")
    n_exceso = int(np.count_nonzero(grados > GRADO_MAX_HIJOS))
    print(f"     grado MAXIMO de salida : {grado_max}   "
          f"(limite GRADO_MAX_HIJOS = {GRADO_MAX_HIJOS})")
    print(f"     nodos con > {GRADO_MAX_HIJOS} hijos    : {n_exceso}   "
          f"{'OK (sin mega-hubs)' if n_exceso == 0 else 'FALLO: hay mega-hubs'}")

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

    # --- metrica de exito: % de glomerulos ALCANZADOS -----------------------
    # (NO se reporta cobertura tisular: no aplica al arbol arterial; la perfusion
    #  del tejido corresponde a la capa capilar posterior.)
    print(f"\n  METRICA DE EXITO: % de glomerulos alcanzados = {pct_alcanzados:.2f} %")
    print("-" * 70)
    print(f"     {n_alcanzados} / {n_glom} glomerulos conectados al arbol arterial")
    print(f"     {'OK: mantiene/supera la v1 (92.92%)' if pct_alcanzados >= 92.92 else 'REVISAR: por debajo de la v1 (92.92%)'}")

    print("\n  -> guardado en:", out_path)
    print("=" * 70)


if __name__ == "__main__":
    main()
