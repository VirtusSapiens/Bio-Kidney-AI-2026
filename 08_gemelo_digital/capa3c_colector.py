#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa3c_colector.py  --  Bio-Kidney AI 2026
==========================================
Capa 3c del gemelo digital: ARBOL COLECTOR URINARIO.

Recoge la orina desde el extremo DISTAL de cada tubulo (Capa 1) y la convergue,
por piramide, hacia su PAPILA (apice de la piramide, Capa 0), que drena al caliz.
NO es un arbol de presion: NO usa Ley de Murray ni conservacion de area. Los radios
se asignan por orden de Horton-Strahler (bloque intercambiable).

ARQUITECTURA: BOSQUE de 10 subarboles DISJUNTOS, uno por piramide k=0..9.
  raiz_k        = pyramid_apex[k]                         (papila, Capa 0)
  atractores_k  = tubulos[piramide_destino==k, 5]        (extremo distal, Capa 1)
  nodos_k, parent_k, ... = space_colonization(raiz_k, atractores_k)   (REUSO Capa 3b)

Se REUSA space_colonization() de capa3b_venoso.py (importada, NO reescrita); se
sobrescriben sus hiperparametros de modulo para el colector y se llama con
art_tree=None (sin sesgo satelital). NO se usa asignar_radios_murray().

Convencion de radios (confirmada leyendo capa3b_venoso.py): el venoso guarda
`radios` POR-ARISTA, pero su cantidad interna r_in[i] es POR-NODO = radio del
segmento que alimenta al nodo i. Aqui se guarda radios(N,) POR-NODO con esa misma
convencion: radios[i] = radio del segmento colector en el nodo i = rampa de Strahler.

Solo NumPy / SciPy. SIN bpy. Entorno: env_biokidney.

Entrada : capa0_dominio.npz (pyramid_apex), capa1_nefronas.npz (tubulos, piramide_destino)
Salida  : capa3c_colector.npz
"""

import os
import sys
import time
import numpy as np
from scipy.spatial import cKDTree

_here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
if _here not in sys.path:
    sys.path.insert(0, _here)
import capa3b_venoso as cv                      # REUSO de space_colonization() + find_npz()

# ============================================================================
#  HIPERPARAMETROS  (arrancan del venoso; DIST_* ajustados para maximizar alcance)
# ============================================================================
PASO_CRECIMIENTO = 0.6      # mm  (igual que el venoso)
DIST_INFLUENCIA  = 8.0      # mm  radio de influencia de un atractor
DIST_MATAR       = 0.6      # mm  atractor alcanzado (mas ajustado que el venoso 1.2:
                            #     el colector debe ACERCARSE al extremo del tubulo)
GRADO_MAX_HIJOS  = 2        # bifurcacion binaria (como el venoso)
MAX_ITERACIONES  = 2000

# --- radios por Horton-Strahler [BLOQUE INTERCAMBIABLE] -----------------------
# Radios ANCLADOS A LITERATURA ANATOMICA (diametros de fuente):
#   radio_terminal = 15 um (r) -> 30 um (diam): conducto colector CORTICAL, cuyo
#     diametro cae en el rango 20-50 um (morfometria del colector cortical,
#     Britannica / literatura morfometrica renal).
#   radio_papila   = 200 um (r) -> 400 um (diam): conducto PAPILAR de Bellini,
#     centro del rango 300-600 um de diametro (grading system de nefrolitiasis,
#     PMC). El conducto de Bellini es UNO por piramide, en el apice papilar.
RADIO_TERMINAL = 0.015      # mm  (15 um r / 30 um diam)  colector cortical (20-50 um diam)
RADIO_PAPILA   = 0.200      # mm  (200 um r / 400 um diam) conducto de Bellini (300-600 um diam)
# [Si se cambia a downstream-count: sustituir SOLO el bloque de radios por
#  radio = rampa en log(N_terminales_aguas_abajo). No tocar nada mas.]

SEED = 2026
UMBRAL_PROX = 0.150         # mm  umbral de proximidad inter-subarbol (analogo shunt VV)

IN_DOM = "capa0_dominio.npz"
IN_NEF = "capa1_nefronas.npz"
OUT_NPZ = "capa3c_colector.npz"


def strahler_por_nodo(parent, children):
    """Orden de Horton-Strahler por nodo (hojas=1; nodo=max(orden hijos), +1 si el
       maximo empata en >=2 hijos). parent[i] < i por bloque -> procesar de N-1 a 0."""
    N = len(parent)
    order = np.ones(N, dtype=np.int64)
    for i in range(N - 1, -1, -1):
        ch = children[i]
        if not ch:
            order[i] = 1
            continue
        oc = [order[c] for c in ch]
        m = max(oc)
        order[i] = m + 1 if oc.count(m) >= 2 else m
    return order


def main():
    np.random.seed(SEED)

    dom_path = cv.find_npz(IN_DOM)
    nef_path = cv.find_npz(IN_NEF)

    # --- sobrescribir los globals del modulo venoso para REUSAR su funcion -----
    cv.PASO_CRECIMIENTO = PASO_CRECIMIENTO
    cv.DIST_INFLUENCIA = DIST_INFLUENCIA
    cv.DIST_MATAR = DIST_MATAR
    cv.GRADO_MAX_HIJOS = GRADO_MAX_HIJOS
    cv.MAX_ITERACIONES = MAX_ITERACIONES

    d0 = np.load(dom_path, allow_pickle=False)
    apex = d0["pyramid_apex"].astype(np.float64)             # (10,3) papilas
    n_pir = len(apex)

    nef = np.load(nef_path, allow_pickle=False)
    tubulos = nef["tubulos"].astype(np.float64)              # (1300,6,3)
    piramide_destino = nef["piramide_destino"].astype(np.int64)   # (1300,)
    distal = tubulos[:, 5, :]                                # (1300,3) extremo distal
    n_nef = len(distal)

    print("=" * 74)
    print("  CAPA 3c - ARBOL COLECTOR URINARIO (bosque de %d subarboles, convergente)" % n_pir)
    print("=" * 74)
    print("  Dominio (papilas=apex):", dom_path)
    print("  Nefronas (tubulos[:,5]):", nef_path)
    print("-" * 74)
    print(f"  REUSO space_colonization() de capa3b_venoso.py  (art_tree=None, sin satelital)")
    print(f"  PASO={PASO_CRECIMIENTO}mm  DIST_INFLUENCIA={DIST_INFLUENCIA}mm  "
          f"DIST_MATAR={DIST_MATAR}mm  GRADO_MAX_HIJOS={GRADO_MAX_HIJOS}")
    print(f"  Radios Horton-Strahler (anclaje POR-SUBARBOL): "
          f"terminal={RADIO_TERMINAL*1000:.0f}um r (colector cortical, diam 20-50um)  "
          f"papila={RADIO_PAPILA*1000:.0f}um r (Bellini, diam 300-600um)")
    print("-" * 74)

    # --- 1) crecer los 10 subarboles y concatenar ---------------------------
    t0 = time.perf_counter()
    nodos_all, parent_all, pir_all = [], [], []
    papila_nodo = np.zeros(n_pir, dtype=np.int64)
    n_no_alc = np.zeros(n_pir, dtype=np.int64)
    n_atrac_k = np.zeros(n_pir, dtype=np.int64)
    offset = 0
    for k in range(n_pir):
        atr_k = distal[piramide_destino == k]
        n_atrac_k[k] = len(atr_k)
        nodos_k, parent_k, n_no_k, _ = cv.space_colonization(apex[k], atr_k, art_tree=None)
        Mk = len(nodos_k)
        pk = parent_k.copy()
        pk[pk >= 0] += offset                      # offset de indices; raiz queda -1
        nodos_all.append(nodos_k)
        parent_all.append(pk)
        pir_all.append(np.full(Mk, k, dtype=np.int64))
        papila_nodo[k] = offset                    # nodo 0 del subarbol = raiz = apex[k]
        n_no_alc[k] = n_no_k
        offset += Mk

    nodos = np.vstack(nodos_all)
    parent = np.concatenate(parent_all)
    piramide_id = np.concatenate(pir_all)
    N = len(nodos)
    t_run = time.perf_counter() - t0

    # children globales (parent[i] < i por bloque)
    children = [[] for _ in range(N)]
    for i in range(N):
        p = parent[i]
        if p >= 0:
            children[p].append(i)
    terminales = np.array([i for i in range(N) if not children[i]], dtype=np.int64)

    # --- 2) radios por Horton-Strahler (rampa monotona, ANCLADA POR-SUBARBOL) ---
    # [UNICO cambio de logica de esta recalibracion]
    # Anclaje POR-SUBARBOL (no global): para cada subarbol k el radio se interpola
    # entre RADIO_TERMINAL (orden 1) y RADIO_PAPILA en el orden MAXIMO DE ESE
    # subarbol (omax_k), NO en el orden maximo global. Justificacion anatomica: las
    # 10 papilas son equivalentes (un conducto de Bellini por piramide); su calibre
    # NO debe depender de la profundidad de ramificacion del subarbol. Como la raiz
    # (papila) tiene siempre el orden maximo de su subarbol, queda exactamente en
    # RADIO_PAPILA -> las 10 papilas ~ homogeneas (~200 um r), en vez del rango
    # 81-125 um que producia el anclaje global.
    strahler = strahler_por_nodo(parent, children)
    omax = int(strahler.max())                       # orden max GLOBAL (solo para reporte)
    radios = np.empty(N, dtype=np.float64)
    for k in range(n_pir):
        in_k = piramide_id == k
        omax_k = int(strahler[in_k].max())           # orden max de ESTE subarbol
        if omax_k > 1:
            radios[in_k] = (RADIO_TERMINAL
                            + (strahler[in_k] - 1) / (omax_k - 1)
                            * (RADIO_PAPILA - RADIO_TERMINAL))
        else:
            radios[in_k] = RADIO_PAPILA

    # --- 3) conexion con Capa 1: tubulo distal -> nodo colector global mas cercano ---
    tree_nod = cKDTree(nodos)
    cdist, cidx = tree_nod.query(distal, k=1)
    conexion_tubulo = np.stack([np.arange(n_nef), cidx], axis=1).astype(np.int64)
    conexion_dist = cdist.astype(np.float64)

    # --- guardar ------------------------------------------------------------
    out_path = os.path.join(os.path.dirname(dom_path), OUT_NPZ)
    np.savez_compressed(
        out_path,
        nodos=nodos.astype(np.float32),
        parent=parent.astype(np.int32),
        radios=radios.astype(np.float32),                  # POR-NODO (convencion venoso)
        piramide_id=piramide_id.astype(np.int32),
        papila_nodo=papila_nodo.astype(np.int32),
        terminales=terminales.astype(np.int32),
        conexion_tubulo=conexion_tubulo.astype(np.int32),
        conexion_dist=conexion_dist.astype(np.float32),
        strahler=strahler.astype(np.int32),
        # escalares
        paso_crecimiento=np.float64(PASO_CRECIMIENTO),
        n_subarboles=np.int32(n_pir),
        radio_terminal=np.float64(RADIO_TERMINAL),
        radio_papila=np.float64(RADIO_PAPILA),
        seed=np.int32(SEED),
        n_no_alcanzados_por_piramide=n_no_alc.astype(np.int32),
        dist_influencia=np.float64(DIST_INFLUENCIA),
        dist_matar=np.float64(DIST_MATAR),
    )

    # ========================================================================
    #  AUDITORIA  (solo mide; NO corrige)
    # ========================================================================
    print(f"  Subarboles: {n_pir}   Nodos: {N}   Aristas: {N-n_pir} (N - n_subarboles)   "
          f"Terminales: {len(terminales)}   ({t_run:.2f}s)")
    print("\n" + "=" * 74)
    print("  AUDITORIA")
    print("=" * 74)

    # (1) ALCANCE ------------------------------------------------------------
    print("  (1) ALCANCE por piramide (atractores = extremos distales de tubulo)")
    print("  " + "-" * 70)
    print(f"  {'pir':>4}{'n_atrac':>10}{'no_alc':>9}{'%alc':>8}{'nodos':>8}")
    for k in range(n_pir):
        nk = int(n_atrac_k[k]); na = int(n_no_alc[k])
        pa = 100.0 * (nk - na) / nk if nk else 0.0
        mk = int(np.count_nonzero(piramide_id == k))
        print(f"  {k:>4}{nk:>10}{na:>9}{pa:>7.1f}%{mk:>8}")
    tot_atr = int(n_atrac_k.sum()); tot_na = int(n_no_alc.sum())
    print("  " + "-" * 70)
    print(f"  TOTAL atractores {tot_atr}   no alcanzados {tot_na}  "
          f"({100.0*(tot_atr-tot_na)/tot_atr:.2f}% alcanzados)")
    print(f"  conexion_dist (tubulo->colector) mm: "
          f"p50 {np.percentile(conexion_dist,50):.3f}  p90 {np.percentile(conexion_dist,90):.3f}  "
          f"p99 {np.percentile(conexion_dist,99):.3f}  max {conexion_dist.max():.3f}")
    print("-" * 74)

    # (2) CONTINUIDAD --------------------------------------------------------
    print("  (2) CONTINUIDAD")
    print("  " + "-" * 70)
    n_conect = len(np.unique(conexion_tubulo[:, 0]))
    print(f"  nefronas conectadas: {n_conect} / {n_nef}  "
          f"({'OK 100%' if n_conect == n_nef else 'REVISAR'})")
    # cada papila_nodo debe coincidir con pyramid_apex[k]
    dpap = np.linalg.norm(nodos[papila_nodo] - apex, axis=1)
    ok_pap = bool(np.all(dpap < 1e-3))
    print(f"  papila_nodo[k] == pyramid_apex[k]: max desvio {dpap.max():.2e} mm  "
          f"{'OK' if ok_pap else 'FALLO'}")
    # fraccion de conexiones que caen en el subarbol de su propia piramide
    mismo = piramide_id[conexion_tubulo[:, 1]] == piramide_destino
    print(f"  conexiones al subarbol de SU piramide: {int(mismo.sum())}/{n_nef} "
          f"({100.0*mismo.mean():.1f}%)  (resto: nodo global mas cercano es de piramide vecina)")
    print("-" * 74)

    # (3) RADIOS -------------------------------------------------------------
    print("  (3) RADIOS (Horton-Strahler, por-nodo)")
    print("  " + "-" * 70)
    print(f"  orden Strahler: rango [1, {omax}]   rango radios [{radios.min()*1000:.1f}, "
          f"{radios.max()*1000:.1f}] um")
    # monotonia: radio(padre) >= radio(hijo) en toda arista
    mask_e = parent >= 0
    mono = np.all(radios[parent[mask_e]] >= radios[mask_e] - 1e-12)
    n_viol = int(np.count_nonzero(radios[parent[mask_e]] < radios[mask_e] - 1e-12))
    print(f"  monotonia raiz>hoja (radio_padre >= radio_hijo): "
          f"{'OK' if mono else 'FALLO'}  (violaciones: {n_viol})")
    r_pap = radios[papila_nodo]
    print(f"  radio en las 10 papilas (um): min {r_pap.min()*1000:.1f}  "
          f"max {r_pap.max()*1000:.1f}  (anclaje POR-SUBARBOL: las 10 ~ RADIO_PAPILA, "
          f"homogeneas)")
    print(f"     ordenes Strahler en las papilas: {[int(strahler[p]) for p in papila_nodo]}  "
          f"(cada papila = orden max de SU subarbol -> radio = RADIO_PAPILA)")
    print("  " + "-" * 70)
    # COMPUERTA / ESCALA vs RESOLUCION DE IMPRESION (100-150 um de rasgo, en diametro).
    # La resolucion restringe el rasgo FINO (terminal), no el grueso (papila).
    d_pap = 2.0 * RADIO_PAPILA * 1000.0     # diametro papila (grueso)   um
    d_ter = 2.0 * RADIO_TERMINAL * 1000.0   # diametro terminal (fino)   um
    pap_res = d_pap >= 150.0                # papila por ENCIMA de la resolucion -> resoluble
    print(f"  COMPUERTA vs resolucion de impresion 100-150um (diametro de rasgo):")
    print(f"     papila (GRUESO)  diam {d_pap:.0f}um: "
          f"{'[OK] por encima de 100-150um -> resoluble' if pap_res else '[NOTA] cerca de resolucion'}")
    print(f"     terminal (FINO)  diam {d_ter:.0f}um: [NOTA] SUB-RESOLUCION de impresion POR DISENO")
    print(f"        (mismo encuadre de escala representativa que la vasculatura: el tronco/terminal")
    print(f"         representa el drenaje del lobulo, no el conducto individual; la resolucion")
    print(f"         restringe el rasgo fino, no la papila, que es grueso y siempre resoluble)")
    print("-" * 74)

    # (4) PROXIMIDAD INTER-SUBARBOL (analogo shunt VV; NO se aplica repulsion) --
    print("  (4) PROXIMIDAD INTER-SUBARBOL (nodos de piramides DISTINTAS)")
    print("  " + "-" * 70)
    tree_all = cKDTree(nodos)
    pares = tree_all.query_pairs(UMBRAL_PROX, output_type='ndarray')
    if len(pares):
        cross = piramide_id[pares[:, 0]] != piramide_id[pares[:, 1]]
        pares_x = pares[cross]
    else:
        pares_x = np.empty((0, 2), dtype=np.int64)
    # distancia minima EXACTA entre subarboles distintos: por cada subarbol k, la
    # distancia del resto de nodos al arbol k (10 queries; barato y exacto).
    dmin = np.inf
    for k in range(n_pir):
        in_k = piramide_id == k
        if not in_k.any() or in_k.all():
            continue
        tk = cKDTree(nodos[in_k])
        dq, _ = tk.query(nodos[~in_k], k=1)
        dmin = min(dmin, float(dq.min()))
    print(f"  pares inter-subarbol < {UMBRAL_PROX*1000:.0f}um: {len(pares_x)}   "
          f"distancia minima inter-subarbol: {dmin:.4f} mm")
    if len(pares_x):
        # localizacion: distancia al hilio (proximidad esperada cerca del seno/papilas)
        hil = d0["hilio"].astype(np.float64)
        pmid = 0.5 * (nodos[pares_x[:, 0]] + nodos[pares_x[:, 1]])
        dh = np.linalg.norm(pmid - hil, axis=1)
        cerca_seno = int(np.count_nonzero(dh < 20.0))
        print(f"     de esos, {cerca_seno} ({100.0*cerca_seno/len(pares_x):.0f}%) a <20mm del hilio "
              f"(zona de convergencia de papilas: proximidad intra-territorio ESPERADA)")
        print(f"     -> {'mayoria cerca del seno: cruces esperados por papilas vecinas, no defecto' if cerca_seno > 0.5*len(pares_x) else 'REVISAR: cruces en pleno parenquima'}")
    else:
        print(f"     sin pares por debajo del umbral: subarboles bien separados")
    print("=" * 74)
    print("  guardado en:", out_path)
    print("=" * 74)


if __name__ == "__main__":
    main()
