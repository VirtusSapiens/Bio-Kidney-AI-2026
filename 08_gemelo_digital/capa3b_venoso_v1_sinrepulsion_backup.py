#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa3b_venoso.py  --  Bio-Kidney AI 2026
========================================
Capa 3b del gemelo digital: ARBOL VENOSO renal.

Recoge la sangre desde la SEGUNDA red capilar (peritubular + vasa recta, Capa 3ab)
y la devuelve al hilio por la vena renal. A diferencia del arterial (que DIVERGE
del hilio hacia los glomerulos), el venoso CONVERGE: nace en muchos puntos de
drenaje finos y se une en troncos cada vez mas gruesos hacia el hilio.

  glomerulo -> eferente -> capilares peritubulares / vasa recta (Capa 3ab)
            -> [PUNTOS DE DRENAJE]  -> VENULAS -> vena renal (esta capa)

Algoritmicamente se reutiliza la MISMA maquinaria depurada de la Capa 3a (space
colonization de Runions 2005, bifurcacion BINARIA con GRADO_MAX_HIJOS, reasignacion
de atractores a los hijos al cerrarse un nodo, reconstruccion del cKDTree cada K
iteraciones y bootstrapping del tronco desde el hilio). El arbol resultante esta
ENRAIZADO en el hilio y alcanza todos los puntos de drenaje: el caracter
"convergente" es la interpretacion del FLUJO (de las hojas hacia la raiz) y de los
RADIOS de Murray (raiz = el tronco mas grueso).

Dos diferencias fisiologicas respecto al arterial:
  1. RADIOS MAYORES para el mismo caudal: la vena opera a baja presion, su lumen es
     mas ancho. Se modela con RADIO_TERMINAL_VENOSO (> arterial) y un
     FACTOR_RADIO_VENOSO global opcional.
  2. TRAYECTO SATELITAL: la vena corre junto a la arteria pero DESPLAZADA
     (~DESPLAZAMIENTO_SATELITAL), entrelazada sin montarse encima. Se implementa
     como un SESGO SUAVE en la direccion de crecimiento (no evasion dura, no
     recalcula rutas): un termino que empuja cada rama a mantener esa separacion
     del segmento arterial mas cercano (KD-tree contra los nodos arteriales),
     ponderado por PESO_SATELITAL.

Solo NumPy / SciPy. SIN bpy. Entorno: env_biokidney.

Entrada : capa3ab_peritubular.npz (puntos_drenaje = atractores), capa0_dominio.npz
          (hilio = raiz), capa3a_arterial.npz (arterial, para el sesgo satelital).
Salida  : capa3b_venoso.npz
"""

import os
import time
import numpy as np
from scipy.spatial import cKDTree

# ============================================================================
#  PARAMETROS  (editar aqui)
# ============================================================================
# RAIZ = hilio (se toma de capa0; la vena renal sale por el MISMO puerto que entra
#               la arteria). PUNTOS_ATRACCION = puntos_drenaje de la Capa 3ab.
DIST_INFLUENCIA = 8.0     # mm  radio en que un atractor influye en el nodo mas cercano
DIST_MATAR      = 1.2     # mm  un atractor a esta distancia de un nodo se da por alcanzado
PASO_CRECIMIENTO = 0.6    # mm  longitud de cada segmento nuevo por iteracion
MAX_ITERACIONES = 2000    # tope de seguridad
GRADO_MAX_HIJOS = 2       # max. hijos por nodo (bifurcacion BINARIA, como el arterial)

# --- sesgo SATELITAL respecto al arbol arterial -------------------------------
DESPLAZAMIENTO_SATELITAL = 1.5   # mm  separacion preferente respecto al arterial
PESO_SATELITAL = 0.3             # cuanto pesa el sesgo satelital vs la atraccion a
                                 # los puntos (0 = ignora arterial; 1 = domina).
                                 # Suave a proposito: crecer cerca pero desplazado.

# --- radios (vena = lumen mas ancho que arteria por baja presion) -------------
MURRAY_EXP             = 3.0     # r^3 = sum(r_hijas^3)
RADIO_TERMINAL_VENOSO  = 0.018   # mm  venula terminal. MAYOR que arterial (0.012)
FACTOR_RADIO_VENOSO    = 1.3     # escala global del lumen venoso sobre Murray. Se
                                 # aplica a TODO el arbol tras Murray: refleja que,
                                 # a igual caudal, la vena es mas ancha que la
                                 # arteria. Documentado: el radio terminal/raiz que
                                 # se reporta es POST-factor.

# --- rendimiento del cKDTree (igual que la Capa 3a) ---------------------------
REBUILD_CADA_K  = 10
REBUILD_NUEVOS  = 500

SEED = 2026

IN_PERI = "capa3ab_peritubular.npz"
IN_DOM = "capa0_dominio.npz"
IN_ART = "capa3a_arterial.npz"
OUT_NPZ = "capa3b_venoso.npz"


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


def sesgo_satelital(pos, base_dir, art_tree):
    """Aplica el SESGO SATELITAL a unas direcciones de crecimiento.

    pos      (n,3)  posiciones de los nodos que crecen
    base_dir (n,3)  direccion unitaria de crecimiento por atraccion a los puntos
    art_tree cKDTree de los nodos arteriales (o None -> no hay sesgo)

    Para cada nodo se busca el nodo arterial mas cercano y se mide la distancia d.
    El sesgo empuja a mantener d ~ DESPLAZAMIENTO_SATELITAL:
      - d < objetivo  -> empuja ALEJANDOSE del arterial (no montarse encima)
      - d > objetivo  -> tira suavemente HACIA el arterial (correr a su lado)
      - d ~ objetivo  -> sesgo nulo
    El termino es radial (en la direccion arterial->nodo), de magnitud proporcional
    al error normalizado y acotada a 1, y se mezcla con base_dir ponderado por
    PESO_SATELITAL. NO es evasion dura: solo inclina la direccion, no recalcula rutas.

    Devuelve la direccion unitaria resultante (n,3).
    """
    if art_tree is None or PESO_SATELITAL <= 0.0:
        return base_dir
    d, idx = art_tree.query(pos, k=1)
    cercano = art_tree.data[idx]                 # nodo arterial mas cercano
    away = pos - cercano                          # arterial -> nodo (radial)
    dn = np.linalg.norm(away, axis=1)
    safe = dn > 1e-9
    away_u = np.zeros_like(away)
    away_u[safe] = away[safe] / dn[safe, None]
    # error normalizado: >0 si estamos demasiado cerca (hay que alejarse)
    err = np.clip((DESPLAZAMIENTO_SATELITAL - dn) / DESPLAZAMIENTO_SATELITAL, -1.0, 1.0)
    bias = err[:, None] * away_u                  # +away (alejar) o -away (acercar)
    nuevo = base_dir + PESO_SATELITAL * bias
    nn = np.linalg.norm(nuevo, axis=1, keepdims=True)
    nn[nn < 1e-9] = 1.0
    return nuevo / nn


# ============================================================================
#  SPACE COLONIZATION CONVERGENTE  (con sesgo satelital)
# ============================================================================
def space_colonization(raiz, atractores, art_tree=None):
    """Crece el arbol desde 'raiz' (3,) hacia 'atractores' (A,3).

    Identico al de la Capa 3a (Runions 2005 + correccion de mega-hubs por
    GRADO_MAX_HIJOS + cKDTree incremental + bootstrapping de tronco), con UNA
    adicion: la direccion de cada rama nueva se pasa por sesgo_satelital() para
    mantener el trayecto DESPLAZADO del arterial. El sesgo se aplica tanto en el
    crecimiento normal como en el bootstrap del tronco.

    Devuelve (nodos (M,3), parent (M,) int con parent[0]=-1, n_no_alcanzados, n_iter).
    """
    atrac = np.asarray(atractores, dtype=np.float64).copy()

    # --- buffers de nodos (crecimiento amortizado) --------------------------
    cap = 1024
    pos_buf = np.empty((cap, 3), dtype=np.float64)
    parent_buf = np.empty(cap, dtype=np.int64)
    nhijos_buf = np.zeros(cap, dtype=np.int64)
    abierto_buf = np.zeros(cap, dtype=bool)
    pos_buf[0] = np.asarray(raiz, dtype=np.float64)
    parent_buf[0] = -1
    abierto_buf[0] = True
    M = 1

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
        nonlocal M
        _ensure(M + 1)
        pos_buf[M] = pos
        parent_buf[M] = p
        nhijos_buf[M] = 0
        abierto_buf[M] = True
        M += 1
        nhijos_buf[p] += 1
        if nhijos_buf[p] >= GRADO_MAX_HIJOS:
            abierto_buf[p] = False
        state["nuevos"] += 1

    def reconstruir():
        ids = np.nonzero(abierto_buf[:M])[0]
        state["tree"] = cKDTree(pos_buf[ids])
        state["ids"] = ids
        state["iters"] = 0
        state["nuevos"] = 0

    def matar(atrac, nuevos_pos):
        if len(atrac) == 0 or len(nuevos_pos) == 0:
            return atrac
        dk, _ = cKDTree(np.asarray(nuevos_pos)).query(atrac, k=1)
        return atrac[dk > DIST_MATAR]

    atrac = matar(atrac, [pos_buf[0]])
    reconstruir()

    n_iter = 0
    for it in range(MAX_ITERACIONES):
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
            # BOOTSTRAP del tronco desde el hilio hasta entrar en el parenquima.
            a = int(np.argmin(dist))
            nd = int(near[a])
            d = atrac[a] - pos_buf[nd]
            dn = np.linalg.norm(d)
            if dn < 1e-9:
                break
            base = (d / dn)[None, :]
            sesg = sesgo_satelital(pos_buf[nd][None, :], base, art_tree)[0]
            nuevo = pos_buf[nd] + PASO_CRECIMIENTO * sesg
            agregar_hijo(nd, nuevo)
            atrac = matar(atrac, [nuevo])
            reconstruir()
            continue

        activo = en_rango & abierto_buf[near]
        if not np.any(activo):
            reconstruir()
            continue

        near_act = near[activo]
        vecs = atrac[activo] - pos_buf[near_act]
        norm = np.linalg.norm(vecs, axis=1, keepdims=True)
        norm[norm == 0] = 1.0
        vecs /= norm

        uniq, inv = np.unique(near_act, return_inverse=True)
        acc = np.zeros((len(uniq), 3), dtype=np.float64)
        np.add.at(acc, inv, vecs)
        gn = np.linalg.norm(acc, axis=1)
        keep = gn > 1e-9
        uniq = uniq[keep]
        g = acc[keep] / gn[keep, None]            # direccion media por atraccion
        if len(uniq) == 0:
            reconstruir()
            continue

        # SESGO SATELITAL: inclina la direccion para mantener el offset arterial
        g = sesgo_satelital(pos_buf[uniq], g, art_tree)

        nuevos = pos_buf[uniq] + PASO_CRECIMIENTO * g
        nuevos_pos = []
        for k in range(len(uniq)):
            agregar_hijo(int(uniq[k]), nuevos[k])
            nuevos_pos.append(nuevos[k])

        atrac = matar(atrac, nuevos_pos)

    return (pos_buf[:M].copy(), parent_buf[:M].copy(), len(atrac), n_iter)


# ============================================================================
#  RADIOS POR LEY DE MURRAY  (de hojas a raiz)
# ============================================================================
def asignar_radios_murray(parent, radio_terminal):
    """Radio del segmento que ALIMENTA a cada nodo (su arista al padre).

    - Hoja (sin hijos) -> radio_terminal.
    - Nodo con hijos   -> (sum r_hijo^MURRAY_EXP)^(1/MURRAY_EXP).
    parent[i] < i, asi que recorrer i de M-1 a 0 procesa hijos antes que padres.

    Devuelve (r_in (M,), hijos (lista)). NOTA: el FACTOR_RADIO_VENOSO global se
    aplica fuera, en main(), sobre el resultado completo.
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
            r_in[i] = radio_terminal
        else:
            r_in[i] = sum(r_in[c] ** MURRAY_EXP for c in h) ** (1.0 / MURRAY_EXP)
    return r_in, hijos


# ============================================================================
#  MAIN
# ============================================================================
def main():
    np.random.seed(SEED)

    peri_path = find_npz(IN_PERI)
    dom_path = find_npz(IN_DOM)
    art_path = find_npz(IN_ART)

    print("=" * 70)
    print("  CAPA 3b - ARBOL VENOSO (Space Colonization CONVERGENTE + satelital)")
    print("=" * 70)
    print("  Entrada drenaje (atractores):", peri_path)
    print("  Entrada dominio (hilio/raiz):", dom_path)
    print("  Entrada arterial (satelital):", art_path)
    print("-" * 70)
    print(f"  DIST_INFLUENCIA = {DIST_INFLUENCIA} mm   DIST_MATAR = {DIST_MATAR} mm")
    print(f"  PASO_CRECIMIENTO = {PASO_CRECIMIENTO} mm   MAX_ITERACIONES = {MAX_ITERACIONES}")
    print(f"  GRADO_MAX_HIJOS = {GRADO_MAX_HIJOS} (binario)")
    print(f"  DESPLAZAMIENTO_SATELITAL = {DESPLAZAMIENTO_SATELITAL} mm   "
          f"PESO_SATELITAL = {PESO_SATELITAL}")
    print(f"  MURRAY_EXP = {MURRAY_EXP}   RADIO_TERMINAL_VENOSO = {RADIO_TERMINAL_VENOSO} mm   "
          f"FACTOR_RADIO_VENOSO = {FACTOR_RADIO_VENOSO}")
    print("-" * 70)

    # --- cargar datos -------------------------------------------------------
    peri = np.load(peri_path, allow_pickle=False)
    puntos_drenaje = peri["puntos_drenaje"].astype(np.float64)   # atractores
    region_dren = peri["region"]                                  # cortex / medula

    d0 = np.load(dom_path, allow_pickle=False)
    hilio = d0["hilio"].astype(np.float64)

    art = np.load(art_path, allow_pickle=False)
    nodos_art = art["nodos"].astype(np.float64)
    radio_raiz_art = float(art["radio_raiz"]) if "radio_raiz" in art.files else None

    raiz = hilio.copy()
    n_atrac = len(puntos_drenaje)
    art_tree = cKDTree(nodos_art)                # KD-tree arterial para el sesgo

    print(f"  RAIZ (hilio): {raiz}")
    print(f"  PUNTOS_ATRACCION (drenaje): {n_atrac}")
    print(f"  Nodos arteriales (referencia satelital): {len(nodos_art)}")
    print("-" * 70)

    # --- 1) topologia por space colonization convergente --------------------
    print("  Creciendo arbol venoso por space colonization (satelital) ...")
    t0 = time.perf_counter()
    nodos, parent, n_no_alcanzados, n_iter = space_colonization(
        raiz, puntos_drenaje, art_tree=art_tree)
    t_run = time.perf_counter() - t0

    # --- 2) radios por Murray + factor venoso -------------------------------
    r_in, hijos = asignar_radios_murray(parent, RADIO_TERMINAL_VENOSO)
    r_in *= FACTOR_RADIO_VENOSO                  # lumen venoso mas ancho (global)

    edge_list = [(parent[i], i) for i in range(len(nodos)) if parent[i] >= 0]
    aristas = np.asarray(edge_list, dtype=np.int64)
    radios = np.asarray([r_in[i] for (_, i) in edge_list], dtype=np.float64)

    terminales = np.asarray([i for i in range(len(nodos)) if not hijos[i]],
                            dtype=np.int64)
    radio_raiz = float(r_in[0])
    radio_terminal_efectivo = RADIO_TERMINAL_VENOSO * FACTOR_RADIO_VENOSO

    # --- cobertura de atractores y desglose por region ----------------------
    # Recomputamos "alcanzado" contra el arbol final (un atractor esta drenado si
    # hay un nodo venoso a < DIST_MATAR). Permite el desglose cortex/medula.
    tree_fin = cKDTree(nodos)
    d_at, _ = tree_fin.query(puntos_drenaje, k=1)
    alcanzado = d_at < DIST_MATAR
    n_alcanzados = int(np.count_nonzero(alcanzado))

    es_medula = region_dren == "medula"
    es_cortex = ~es_medula
    n_ctx = int(np.count_nonzero(es_cortex))
    n_med = int(np.count_nonzero(es_medula))
    ctx_alc = int(np.count_nonzero(alcanzado & es_cortex))
    med_alc = int(np.count_nonzero(alcanzado & es_medula))

    # --- guardar ------------------------------------------------------------
    out_path = os.path.join(os.path.dirname(peri_path), OUT_NPZ)
    np.savez_compressed(
        out_path,
        nodos=nodos.astype(np.float32),
        aristas=aristas.astype(np.int32),
        radios=radios.astype(np.float32),
        terminales=terminales.astype(np.int32),
        raiz=np.int32(0),                        # el hilio es el nodo 0
        raiz_pos=raiz.astype(np.float64),
        # metadata: parametros
        dist_influencia=np.float64(DIST_INFLUENCIA),
        dist_matar=np.float64(DIST_MATAR),
        paso_crecimiento=np.float64(PASO_CRECIMIENTO),
        max_iteraciones=np.int32(MAX_ITERACIONES),
        grado_max_hijos=np.int32(GRADO_MAX_HIJOS),
        desplazamiento_satelital=np.float64(DESPLAZAMIENTO_SATELITAL),
        peso_satelital=np.float64(PESO_SATELITAL),
        murray_exp=np.float64(MURRAY_EXP),
        radio_terminal_venoso=np.float64(RADIO_TERMINAL_VENOSO),
        factor_radio_venoso=np.float64(FACTOR_RADIO_VENOSO),
        radio_terminal_efectivo=np.float64(radio_terminal_efectivo),
        radio_raiz=np.float64(radio_raiz),
        seed=np.int32(SEED),
        # metadata: resultados
        n_puntos_drenaje=np.int32(n_atrac),
        n_puntos_alcanzados=np.int32(n_alcanzados),
        n_no_alcanzados=np.int32(n_atrac - n_alcanzados),
        n_nodos=np.int32(len(nodos)),
        n_aristas=np.int32(len(aristas)),
        n_terminales=np.int32(len(terminales)),
        n_iteraciones=np.int32(n_iter),
        t_run_s=np.float64(t_run),
    )

    # ========================================================================
    #  VERIFICACION
    # ========================================================================
    print("\n  VERIFICACION")
    print("-" * 70)
    print(f"  Tiempo de corrida      : {t_run:.2f} s   ({n_iter} iteraciones)")
    print(f"  Nodos generados        : {len(nodos)}")
    print(f"  Aristas generadas      : {len(aristas)}  (arbol: aristas = nodos-1 = {len(nodos)-1})")
    print(f"  Nodos terminales (venulas): {len(terminales)}")

    # --- distribucion de grados de salida (sin mega-hubs) -------------------
    grados = np.array([len(hijos[i]) for i in range(len(nodos))], dtype=np.int64)
    grado_max = int(grados.max())
    print(f"\n  DISTRIBUCION DE GRADOS DE SALIDA (n.º de hijos por nodo)")
    print("-" * 70)
    for g in range(0, grado_max + 1):
        etiqueta = {0: "terminales", 1: "pasantes", 2: "bifurcaciones",
                    3: "trifurcaciones"}.get(g, "")
        print(f"     {g} hijos: {int(np.count_nonzero(grados == g)):>7d}   {etiqueta}")
    n_exceso = int(np.count_nonzero(grados > GRADO_MAX_HIJOS))
    print(f"     grado MAXIMO de salida : {grado_max}   (limite GRADO_MAX_HIJOS = {GRADO_MAX_HIJOS})")
    print(f"     nodos con > {GRADO_MAX_HIJOS} hijos    : {n_exceso}   "
          f"{'OK (sin mega-hubs)' if n_exceso == 0 else 'FALLO: hay mega-hubs'}")

    # --- cobertura de puntos de drenaje -------------------------------------
    pct = 100.0 * n_alcanzados / n_atrac
    print(f"\n  Puntos de drenaje alcanzados: {n_alcanzados} / {n_atrac}  ({pct:.2f} %)")
    print(f"     {'OK >95%' if pct > 95 else 'REVISAR <95%'}")

    # --- COBERTURA POR REGION (critico: el venoso DEBE drenar la medula) ----
    pct_ctx = 100.0 * ctx_alc / n_ctx if n_ctx else 0.0
    pct_med = 100.0 * med_alc / n_med if n_med else 0.0
    print(f"\n  COBERTURA POR REGION (puntos de drenaje alcanzados)")
    print("-" * 70)
    print(f"     cortex : {ctx_alc:>5d} / {n_ctx:<5d}  ({pct_ctx:5.1f} %)")
    print(f"     medula : {med_alc:>5d} / {n_med:<5d}  ({pct_med:5.1f} %)")
    print(f"     -> {'OK: el venoso drena la MEDULA, no solo la corteza' if pct_med > 90 else 'REVISAR: drenaje medular bajo'}")

    # --- radios (raiz venosa debe SUPERAR la raiz arterial) -----------------
    print(f"\n  Radios (Murray x FACTOR_RADIO_VENOSO):")
    print(f"     radio en la RAIZ venosa : {radio_raiz:.4f} mm")
    print(f"     radio TERMINAL efectivo : {radio_terminal_efectivo:.4f} mm  "
          f"(base {RADIO_TERMINAL_VENOSO} x {FACTOR_RADIO_VENOSO})")
    print(f"     rango de aristas        : [{radios.min():.4f}, {radios.max():.4f}] mm")
    if radio_raiz_art is not None:
        ok_v = radio_raiz > radio_raiz_art
        print(f"     raiz ARTERIAL (capa3a)  : {radio_raiz_art:.4f} mm")
        print(f"     -> raiz venosa > raiz arterial: "
              f"{'OK (vena mas ancha por baja presion)' if ok_v else 'REVISAR: vena no mas ancha'}")

    # --- offset satelital logrado (confirma que corre AL LADO, no encima) ---
    d_va, _ = art_tree.query(nodos, k=1)
    print(f"\n  Offset satelital al arterial (distancia nodo venoso -> arterial):")
    print(f"     objetivo : {DESPLAZAMIENTO_SATELITAL:.2f} mm")
    print(f"     media    : {d_va.mean():.3f} mm   mediana: {np.median(d_va):.3f} mm")
    print(f"     p5-p95   : [{np.percentile(d_va,5):.3f}, {np.percentile(d_va,95):.3f}] mm")
    sep_min = float(d_va[1:].min()) if len(d_va) > 1 else 0.0   # excluye la raiz (puerto comun)
    print(f"     -> minima (excl. hilio): {sep_min:.3f} mm   "
          f"{'corre satelital, no encima' if sep_min > 0.05 else 'REVISAR: solapa el arterial'}")

    # NOTA: NO se calcula cobertura tisular (no aplica al arbol venoso, igual que
    # al arterial; la perfusion del tejido es trabajo de la capa capilar).
    print(f"\n  (cobertura tisular NO aplica al arbol venoso, omitida como en 3a)")
    print("\n  -> guardado en:", out_path)
    print("=" * 70)


if __name__ == "__main__":
    main()
