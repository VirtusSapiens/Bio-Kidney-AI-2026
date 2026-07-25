#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa3a_auditoria_murray.py  --  Bio-Kidney AI 2026
==================================================
AUDITORIA numerica de la Ley de Murray sobre el arbol arterial de
capa3a_arterial.npz. NO modifica el arbol: solo lo evalua y produce evidencia
cuantitativa para revision cientifica.

Ley de Murray en una bifurcacion:  r_padre^k = sum(r_hijos^k),  con k = 3.

Pasos:
  1) Reconstruye la estructura padre-hijo (arbol dirigido desde el hilio) y
     reporta la distribucion de grados de salida (0/1/2/3+ hijos).
  2) Audita Murray SOLO en nodos de bifurcacion (>=2 hijos): error relativo
     entre r_padre y (sum r_hijos^3)^(1/3); estadisticos y % bajo tolerancia.
  3) Trata las multifurcaciones (3+ hijos) por separado y honestamente.
  4) Recupera el exponente de Murray EMPIRICO 'k' (verificacion cruzada
     independiente): si la construccion uso k=3, debe salir ~3.

Solo NumPy / SciPy. SIN bpy. Entorno: env_biokidney.

Entrada : capa3a_arterial.npz
Salida  : capa3a_auditoria_murray.npz  (errores por bifurcacion, para histograma)
          capa3a_auditoria_murray.txt  (reporte de texto)
"""

import os
import numpy as np
from collections import defaultdict
from scipy.optimize import brentq, minimize_scalar

# ============================================================================
#  PARAMETROS
# ============================================================================
MURRAY_K = 3.0                      # exponente de diseno (lo que uso la construccion)
TOLERANCIAS = (0.01, 0.05, 0.10)    # 1%, 5%, 10%
TOL_VEREDICTO = 0.05                # tolerancia de la linea de veredicto final

IN_NPZ = "capa3a_arterial.npz"
OUT_NPZ = "capa3a_auditoria_murray.npz"
OUT_TXT = "capa3a_auditoria_murray.txt"


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


def k_empirico_bifurcacion(ratios):
    """Resuelve k tal que sum_i (r_i/r_padre)^k = 1 para UNA bifurcacion.

    'ratios' = r_hijos / r_padre (todos < 1 en un arbol bien formado). La funcion
    h(k)=sum(ratios^k)-1 es estrictamente decreciente: h(0)=n_hijos-1>0 y
    h(+inf)=-1<0, asi que existe una unica raiz. Devuelve k o np.nan si degenera
    (algun hijo >= padre, que rompe la monotonia).
    """
    ratios = np.asarray(ratios, dtype=np.float64)
    if ratios.size < 2 or np.any(ratios >= 1.0) or np.any(ratios <= 0.0):
        return np.nan
    h = lambda k: np.sum(ratios ** k) - 1.0
    try:
        return float(brentq(h, 1e-3, 60.0, xtol=1e-8, rtol=1e-10))
    except Exception:
        return np.nan


# ============================================================================
#  MAIN
# ============================================================================
def main():
    in_path = find_npz(IN_NPZ)
    d = np.load(in_path, allow_pickle=False)
    nodos = d["nodos"].astype(np.float64)
    aristas = d["aristas"].astype(np.int64)          # (E,2) = (padre, hijo)
    radios = d["radios"].astype(np.float64)          # (E,) radio de la arista (= r entrante del hijo)
    terminales_in = d["terminales"].astype(np.int64)
    M = len(nodos)
    E = len(aristas)

    lines = []                                        # buffer del reporte de texto

    def P(s=""):
        print(s)
        lines.append(s)

    P("=" * 72)
    P("  CAPA 3a - AUDITORIA DE LA LEY DE MURRAY (r_padre^k = sum r_hijos^k)")
    P("=" * 72)
    P(f"  Entrada: {in_path}")
    P(f"  Nodos: {M}    Aristas: {E}    k de diseno: {MURRAY_K}")
    P("  (Auditoria de solo lectura: NO modifica el arbol.)")

    # ------------------------------------------------------------------------
    #  PASO 1 - reconstruir estructura padre-hijo y grados de salida
    # ------------------------------------------------------------------------
    # Las aristas ya estan dirigidas (padre, hijo); el radio de la arista es el
    # radio del segmento ENTRANTE al hijo -> r de cada nodo (no-raiz) = node_r[hijo].
    parent = np.full(M, -1, dtype=np.int64)
    node_r = np.full(M, np.nan, dtype=np.float64)     # radio del segmento entrante por nodo
    parent[aristas[:, 1]] = aristas[:, 0]
    node_r[aristas[:, 1]] = radios

    outdeg = np.bincount(aristas[:, 0], minlength=M)  # n.o de hijos por nodo
    roots = np.where(parent < 0)[0]

    n_term = int(np.count_nonzero(outdeg == 0))
    n_cont = int(np.count_nonzero(outdeg == 1))
    n_bin = int(np.count_nonzero(outdeg == 2))
    n_multi = int(np.count_nonzero(outdeg >= 3))
    n_bif = n_bin + n_multi                            # bifurcaciones (>=2 hijos)

    P("\n" + "-" * 72)
    P("  PASO 1 - ESTRUCTURA DEL ARBOL DIRIGIDO")
    P("-" * 72)
    P(f"  Raices detectadas (sin padre): {len(roots)}  -> nodo(s) {roots.tolist()}")
    if len(roots) == 1:
        P(f"    OK: arbol con raiz unica en el indice {int(roots[0])} (hilio).")
    else:
        P("    AVISO: se esperaba exactamente 1 raiz (el hilio).")
    P(f"  Distribucion de grados de salida (n.o de hijos):")
    P(f"     0 hijos  (terminales)        : {n_term:>7d}  ({100.0*n_term/M:5.2f} %)")
    P(f"     1 hijo   (continuacion)      : {n_cont:>7d}  ({100.0*n_cont/M:5.2f} %)")
    P(f"     2 hijos  (bifurcacion binaria): {n_bin:>7d}  ({100.0*n_bin/M:5.2f} %)")
    P(f"     3+ hijos (multifurcacion)    : {n_multi:>7d}  ({100.0*n_multi/M:5.2f} %)")
    P(f"     -> total nodos de bifurcacion (>=2 hijos): {n_bif}")
    # coherencia con 'terminales' guardado en el .npz
    coincide = (n_term == len(terminales_in))
    P(f"  Terminales por grado=0: {n_term}  vs  'terminales' del .npz: "
      f"{len(terminales_in)}  -> {'coinciden' if coincide else 'DIFIEREN'}")

    # ------------------------------------------------------------------------
    #  PASO 2 + 3 - Murray en bifurcaciones (suma de TODOS los hijos^3)
    # ------------------------------------------------------------------------
    # Suma de r_hijo^k por nodo padre (vectorizado). Vale para Murray binario y
    # generalizado: r_pred = (sum_hijos r^3)^(1/3).
    S3 = np.bincount(aristas[:, 0], weights=radios ** MURRAY_K, minlength=M)
    r_pred_all = np.where(S3 > 0, S3 ** (1.0 / MURRAY_K), 0.0)

    # nodos de bifurcacion, EXCLUIDA la raiz (no tiene segmento entrante medible)
    is_bif = (outdeg >= 2) & (parent >= 0)
    bif_ids = np.where(is_bif)[0]
    r_padre = node_r[bif_ids]
    r_pred = r_pred_all[bif_ids]
    err_rel = np.abs(r_padre - r_pred) / r_padre
    nh = outdeg[bif_ids]                               # n.o de hijos por bifurcacion
    es_binaria = nh == 2

    n_raiz_bif = int(np.count_nonzero((outdeg >= 2) & (parent < 0)))

    def stats_block(titulo, mask):
        e = err_rel[mask]
        P(f"  {titulo}  (n={e.size})")
        if e.size == 0:
            P("     (sin casos)")
            return
        P(f"     error rel.  medio={e.mean()*100:.4f} %   mediana={np.median(e)*100:.4f} %")
        P(f"                 maximo={e.max()*100:.4f} %   p95={np.percentile(e,95)*100:.4f} %")
        for t in TOLERANCIAS:
            pct = 100.0 * np.count_nonzero(e <= t) / e.size
            P(f"     cumplen <= {t*100:>4.0f} % : {pct:6.2f} %")

    P("\n" + "-" * 72)
    P("  PASO 2 - AUDITORIA DE MURRAY EN BIFURCACIONES (>=2 hijos, excl. raiz)")
    P("-" * 72)
    P(f"  Comparacion: r_padre  vs  (sum_hijos r^{MURRAY_K:g})^(1/{MURRAY_K:g})")
    if n_raiz_bif:
        P(f"  (La raiz tiene >=2 hijos pero NO se audita: carece de segmento")
        P(f"   entrante medible. radio_raiz se fija por construccion = Murray.)")
    stats_block("TODAS las bifurcaciones:", np.ones(bif_ids.size, dtype=bool))

    P("\n" + "-" * 72)
    P("  PASO 3 - BINARIAS LIMPIAS vs MULTIFURCACIONES")
    P("-" * 72)
    frac_bin = 100.0 * n_bin / n_bif if n_bif else 0.0
    frac_multi = 100.0 * n_multi / n_bif if n_bif else 0.0
    P(f"  Bifurcaciones binarias (2 hijos): {n_bin}  ({frac_bin:.2f} % de las bifurcaciones)")
    P(f"  Multifurcaciones (3+ hijos)     : {n_multi}  ({frac_multi:.2f} % de las bifurcaciones)")

    # Distribucion del n.o de hijos por buckets (las cuentas crudas son ilegibles
    # porque hay 'hubs' con ~1900 hijos). Separamos multifurcaciones pequenas de
    # MEGA-HUBS, claramente no fisiologicos (un vaso no emite cientos de ramas).
    MEGA_HUB = 50                                     # umbral de hub degenerado
    es_mega = nh > MEGA_HUB
    es_multi_peq = (nh >= 3) & (~es_mega)
    n_mega = int(np.count_nonzero(es_mega))
    n_multi_peq = int(np.count_nonzero(es_multi_peq))
    if n_multi:
        P(f"  Distribucion del grado en multifurcaciones (buckets de n.o de hijos):")
        edges_b = [3, 4, 5, 6, 11, 51, 200, 1000, 10000]
        nhm = nh[~es_binaria]
        for lo, hi in zip(edges_b[:-1], edges_b[1:]):
            c = int(np.count_nonzero((nhm >= lo) & (nhm < hi)))
            if c:
                rng = f"{lo}" if hi - lo == 1 else f"{lo}-{hi-1}"
                P(f"     {rng:>9s} hijos: {c}")
        P(f"     grado maximo observado: {int(nh.max())} hijos")

    if n_mega:
        hijos_en_mega = int(nh[es_mega].sum())
        P(f"\n  *** HALLAZGO ESTRUCTURAL (para revision) ***")
        P(f"  {n_mega} nodos son MEGA-HUBS (> {MEGA_HUB} hijos): concentran "
          f"{hijos_en_mega} aristas ({100.0*hijos_en_mega/E:.1f} % del arbol).")
        P(f"  No son bifurcaciones fisiologicas: un vaso real se divide en 2 (rara")
        P(f"  vez 3) ramas. Estos 'abanicos' son un ARTEFACTO del space colonization")
        P(f"  (un mismo nodo fue el mas cercano a un cluster y emitio una rama nueva")
        P(f"  cada iteracion). La Ley de Murray se cumple igual (los radios se")
        P(f"  asignaron con esa formula), pero la TOPOLOGIA no es ramificacion")
        P(f"  recursiva realista. Recomendacion: revisar la generacion en capa3a")
        P(f"  (limitar grado de salida / re-bifurcar los hubs) antes de conclusiones")
        P(f"  morfometricas. Murray OK no implica arbol anatomicamente correcto.")

    stats_block("Solo BINARIAS (2 hijos):", es_binaria)
    if n_multi_peq:
        stats_block(f"Multifurcaciones PEQUENAS (3-{MEGA_HUB} hijos):", es_multi_peq)
    if n_mega:
        stats_block("MEGA-HUBS (> %d hijos, Murray generalizado):" % MEGA_HUB, es_mega)

    # ------------------------------------------------------------------------
    #  PASO 4 - exponente de Murray empirico (verificacion cruzada)
    # ------------------------------------------------------------------------
    # (a) ajuste GLOBAL de un unico k que minimiza el residual log al cuadrado:
    #     residual(k) = sum_bif [ ln(sum_hijos r^k) - k*ln(r_padre) ]^2.
    #     Vectorizado: para cada k, Sk = bincount(padre, weights=r^k).
    padres_e = aristas[:, 0]
    ln_rp = np.log(node_r)                             # ln r_padre por nodo (nan en raiz)

    def residual_global(k):
        Sk = np.bincount(padres_e, weights=radios ** k, minlength=M)
        res = np.log(Sk[bif_ids]) - k * ln_rp[bif_ids]
        return float(np.sum(res ** 2))

    opt = minimize_scalar(residual_global, bounds=(0.5, 8.0), method="bounded")
    k_global = float(opt.x)

    # (b) k EMPIRICO por bifurcacion (resolviendo sum (r_hijo/r_padre)^k = 1) y
    #     su distribucion. Requiere los radios de los hijos agrupados por padre.
    hijos_r = defaultdict(list)
    for p, r in zip(padres_e, radios):
        hijos_r[int(p)].append(float(r))

    k_bif = np.full(bif_ids.size, np.nan, dtype=np.float64)
    for j, nd in enumerate(bif_ids):
        rp = node_r[nd]
        ratios = np.asarray(hijos_r[int(nd)], dtype=np.float64) / rp
        k_bif[j] = k_empirico_bifurcacion(ratios)
    k_validos = k_bif[np.isfinite(k_bif)]

    P("\n" + "-" * 72)
    P("  PASO 4 - EXPONENTE DE MURRAY EMPIRICO (verificacion cruzada)")
    P("-" * 72)
    P(f"  (a) Ajuste GLOBAL (min. residual log^2 sobre todas las bifurcaciones):")
    P(f"        k_global = {k_global:.6f}    (residual = {opt.fun:.3e})")
    if k_validos.size:
        P(f"  (b) k EMPIRICO por bifurcacion (n validos={k_validos.size}/{k_bif.size}):")
        P(f"        media={k_validos.mean():.6f}  mediana={np.median(k_validos):.6f}  "
          f"std={k_validos.std():.3e}")
        P(f"        min={k_validos.min():.6f}  max={k_validos.max():.6f}")
    else:
        P("  (b) k por bifurcacion: sin casos validos.")
    desv = abs(k_global - MURRAY_K)
    P(f"  Desviacion de k_global respecto al k de diseno ({MURRAY_K:g}): {desv:.2e}")
    P(f"  -> {'OK: el exponente empirico recupera el de diseno.' if desv < 0.05 else 'REVISAR: el exponente empirico NO coincide con el de diseno (posible bug).'}")

    # ------------------------------------------------------------------------
    #  VEREDICTO
    # ------------------------------------------------------------------------
    pct_veredicto = 100.0 * np.count_nonzero(err_rel <= TOL_VEREDICTO) / err_rel.size \
        if err_rel.size else 0.0
    k_report = np.median(k_validos) if k_validos.size else k_global

    P("\n" + "=" * 72)
    P("  VEREDICTO")
    P("-" * 72)
    P(f"  Murray se cumple en {pct_veredicto:.2f}% de {err_rel.size} bifurcaciones "
      f"dentro del {TOL_VEREDICTO*100:.0f}% de tolerancia; "
      f"exponente empirico k={k_report:.4f}.")
    P("=" * 72)

    # ------------------------------------------------------------------------
    #  GUARDAR SALIDAS
    # ------------------------------------------------------------------------
    out_npz = os.path.join(os.path.dirname(in_path), OUT_NPZ)
    np.savez_compressed(
        out_npz,
        bif_node_ids=bif_ids.astype(np.int32),
        n_hijos=nh.astype(np.int32),
        es_binaria=es_binaria,
        r_padre=r_padre.astype(np.float64),
        r_murray_pred=r_pred.astype(np.float64),
        err_rel=err_rel.astype(np.float64),
        k_empirico_bif=k_bif.astype(np.float64),
        # resumen
        k_global=np.float64(k_global),
        murray_k_diseno=np.float64(MURRAY_K),
        n_terminales=np.int32(n_term),
        n_continuacion=np.int32(n_cont),
        n_binarias=np.int32(n_bin),
        n_multifurcaciones=np.int32(n_multi),
        n_mega_hubs=np.int32(n_mega),
        umbral_mega_hub=np.int32(MEGA_HUB),
        grado_maximo=np.int32(int(nh.max()) if nh.size else 0),
        n_bifurcaciones=np.int32(n_bif),
        pct_veredicto=np.float64(pct_veredicto),
        tol_veredicto=np.float64(TOL_VEREDICTO),
    )

    out_txt = os.path.join(os.path.dirname(in_path), OUT_TXT)
    with open(out_txt, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"\n  -> errores por bifurcacion: {out_npz}")
    print(f"  -> reporte de texto       : {out_txt}")


if __name__ == "__main__":
    main()
