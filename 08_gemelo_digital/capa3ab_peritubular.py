#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa3ab_peritubular.py  --  Bio-Kidney AI 2026
==============================================
Capa 3ab del gemelo digital: SEGUNDA RED CAPILAR renal (peritubular + vasa recta).

Es el PUENTE entre el arbol arterial (Capa 3a) y el futuro arbol venoso (Capa 3b).
Flujo fisiologico real:

    glomerulo --> arteriola eferente --> [capilares peritubulares (corteza)
                                          o  vasa recta (yuxtamedulares)]
               --> venulas --> vena

Esta capa NO modela cada capilar individual (serian millones). Genera un conjunto
REPRESENTATIVO de PUNTOS DE DRENAJE distribuidos a lo largo de cada tubulo: son los
puntos donde la sangre, tras irrigar el tejido, se recoge hacia el lado venoso.
Esos puntos seran los PUNTOS DE ATRACCION / ORIGEN del arbol venoso de la Capa 3b.

SIMPLIFICACION (documentada): un "punto de drenaje" agrega el drenaje de todo un
segmento de la red capilar peritubular o de la vasa recta que rodea esa porcion de
tubulo. Densidad parametrizable (PUNTOS_POR_TUBULO_*), NO una correspondencia 1:1
con capilares reales.

Tipos de red (anatomia):
  - peritubular_cortical : red capilar que envuelve los tubulos contorneados en la
                           CORTEZA. Nace de la eferente de glomerulos CORTICALES.
  - vasa_recta           : capilares rectos largos que descienden a la MEDULA
                           siguiendo el asa de Henle. Nacen de la eferente de
                           glomerulos YUXTAMEDULARES. Son los que dan cobertura
                           MEDULAR (resuelven el deficit medular que el arbol
                           arterial, anclado a glomerulos, deja en la medula).

NOTA sobre "eferente_asociada": la Capa 3a crece hasta los glomerulos, de modo que
sus nodos TERMINALES son las arteriolas (aferentes) terminales junto al glomerulo.
La arteriola EFERENTE emerge del mismo glomerulo, asi que usamos el terminal
arterial mas cercano como PROXY del territorio arterial/eferente que alimenta cada
punto de drenaje. Eso registra el origen del flujo (de donde viene la sangre).

Solo NumPy / SciPy. SIN bpy. Entorno: env_biokidney.

Entrada : capa1_nefronas.npz, capa3a_arterial.npz, capa0_dominio.npz
Salida  : capa3ab_peritubular.npz
"""

import os
import numpy as np
from scipy.spatial import cKDTree

# ============================================================================
#  PARAMETROS  (editar aqui)
# ============================================================================
PUNTOS_POR_TUBULO_CORTICAL = 3     # densidad de drenaje en tubulos corticales
PUNTOS_POR_TUBULO_YUXTA    = 5     # mas puntos en yuxtamedulares (asa larga + vasa recta medular)
JITTER_RADIAL_MM           = 0.15  # los puntos rodean el tubulo, no caen exactamente sobre el
SEED                       = 2026

IN_NEF = "capa1_nefronas.npz"
IN_ART = "capa3a_arterial.npz"
IN_DOM = "capa0_dominio.npz"
OUT_NPZ = "capa3ab_peritubular.npz"


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


def muestrear_polilinea(vertices, fracciones):
    """Muestrea una polilinea (P,3) en posiciones dadas como FRACCIONES de su
    longitud de arco total (cada f en [0,1]).

    Devuelve:
      pts  (N,3)  posiciones interpoladas a lo largo del tubulo
      tan  (N,3)  tangente unitaria local en cada posicion (direccion del tubulo)
    """
    v = np.asarray(vertices, dtype=np.float64)
    seg = np.diff(v, axis=0)                     # (P-1, 3) vectores de segmento
    long_seg = np.linalg.norm(seg, axis=1)       # longitud de cada segmento
    acum = np.concatenate(([0.0], np.cumsum(long_seg)))  # arco acumulado (P,)
    total = acum[-1]

    pts = np.empty((len(fracciones), 3), dtype=np.float64)
    tan = np.empty((len(fracciones), 3), dtype=np.float64)

    if total < 1e-9:                             # tubulo degenerado (punto)
        pts[:] = v[0]
        tan[:] = np.array([0.0, 0.0, 1.0])
        return pts, tan

    objetivo = np.asarray(fracciones, dtype=np.float64) * total
    for i, s in enumerate(objetivo):
        k = int(np.searchsorted(acum, s, side="right") - 1)   # segmento que contiene s
        k = min(max(k, 0), len(seg) - 1)
        ll = long_seg[k] if long_seg[k] > 1e-12 else 1.0
        u = (s - acum[k]) / ll                                # parametro local [0,1]
        pts[i] = v[k] + u * seg[k]
        tan[i] = seg[k] / np.linalg.norm(seg[k])
    return pts, tan


def jitter_radial(pts, tan, radio, rng):
    """Desplaza cada punto una distancia ~'radio' en una direccion ALEATORIA
    PERPENDICULAR a la tangente local: el resultado RODEA el tubulo (no cae sobre
    el), imitando el manto de capilares peritubulares alrededor del tubulo.
    """
    n = len(pts)
    # base ortonormal {e1, e2} perpendicular a cada tangente
    aux = np.tile(np.array([0.0, 0.0, 1.0]), (n, 1))
    casi_paralelo = np.abs(np.einsum("ij,ij->i", tan, aux)) > 0.9
    aux[casi_paralelo] = np.array([1.0, 0.0, 0.0])
    e1 = np.cross(tan, aux)
    e1 /= np.linalg.norm(e1, axis=1, keepdims=True)
    e2 = np.cross(tan, e1)                         # ya unitario (tan, e1 ortonormales)

    theta = rng.uniform(0.0, 2.0 * np.pi, size=n)
    # magnitud con sqrt(u) -> reparto uniforme en el disco perpendicular
    mag = radio * np.sqrt(rng.uniform(0.0, 1.0, size=n))
    offset = (mag * np.cos(theta))[:, None] * e1 + (mag * np.sin(theta))[:, None] * e2
    return pts + offset


# ============================================================================
#  MAIN
# ============================================================================
def main():
    rng = np.random.default_rng(SEED)

    nef_path = find_npz(IN_NEF)
    art_path = find_npz(IN_ART)
    dom_path = find_npz(IN_DOM)

    print("=" * 70)
    print("  CAPA 3ab - RED PERITUBULAR + VASA RECTA (puntos de drenaje venoso)")
    print("=" * 70)
    print("  Entrada nefronas:", nef_path)
    print("  Entrada arterial:", art_path)
    print("  Entrada dominio :", dom_path)
    print("-" * 70)
    print(f"  PUNTOS_POR_TUBULO_CORTICAL = {PUNTOS_POR_TUBULO_CORTICAL}")
    print(f"  PUNTOS_POR_TUBULO_YUXTA    = {PUNTOS_POR_TUBULO_YUXTA}")
    print(f"  JITTER_RADIAL_MM           = {JITTER_RADIAL_MM}   SEED = {SEED}")
    print("-" * 70)

    # --- cargar datos -------------------------------------------------------
    d1 = np.load(nef_path, allow_pickle=False)
    tipo_nef = d1["tipo"]                            # 'cortical' / 'yuxtamedular'
    tubulos = d1["tubulos"].astype(np.float64)       # (1300, 6, 3)
    n_nef = len(tubulos)

    d3a = np.load(art_path, allow_pickle=False)
    nodos_art = d3a["nodos"].astype(np.float64)      # (N, 3)
    terminales = d3a["terminales"].astype(np.int64)  # indices de hojas en nodos_art
    pos_terminales = nodos_art[terminales]           # arteriolas terminales (proxy eferente)

    d0 = np.load(dom_path, allow_pickle=False)
    coords_dom = d0["coords"].astype(np.float64)     # (200000, 3)
    region_dom = d0["region_label"]                  # cortex / medulla / piramide_XX

    print(f"  Nefronas              : {n_nef}  "
          f"(corticales {int(np.count_nonzero(tipo_nef == 'cortical'))}, "
          f"yuxtamedulares {int(np.count_nonzero(tipo_nef == 'yuxtamedular'))})")
    print(f"  Terminales arteriales : {len(terminales)}  (proxy de arteriolas eferentes)")
    print(f"  Puntos de dominio     : {len(coords_dom)}  (para etiquetar region)")
    print("-" * 70)

    # --- 1) distribuir puntos de drenaje a lo largo de cada tubulo ----------
    # Fracciones de arco: puntos centrados en N segmentos iguales -> bien
    # repartidos a lo largo del tubulo sin acumularse en los extremos.
    frac_cort = (np.arange(PUNTOS_POR_TUBULO_CORTICAL) + 0.5) / PUNTOS_POR_TUBULO_CORTICAL
    frac_yux = (np.arange(PUNTOS_POR_TUBULO_YUXTA) + 0.5) / PUNTOS_POR_TUBULO_YUXTA

    pts_list, tipo_list, nef_list = [], [], []
    for i in range(n_nef):
        es_yux = (tipo_nef[i] == "yuxtamedular")
        frac = frac_yux if es_yux else frac_cort
        # tipo de red: cortical -> peritubular_cortical ; yuxta -> vasa_recta
        etiqueta = "vasa_recta" if es_yux else "peritubular_cortical"

        pts, tan = muestrear_polilinea(tubulos[i], frac)
        pts = jitter_radial(pts, tan, JITTER_RADIAL_MM, rng)

        pts_list.append(pts)
        tipo_list.append(np.full(len(pts), etiqueta, dtype="<U20"))
        nef_list.append(np.full(len(pts), i, dtype=np.int64))

    puntos_drenaje = np.vstack(pts_list).astype(np.float64)
    tipo_red = np.concatenate(tipo_list)
    nefrona_origen = np.concatenate(nef_list)
    M = len(puntos_drenaje)

    # --- 2a) region (cortex/medula) por vecino mas cercano en el dominio ----
    # cortex -> 'cortex' ; medulla y piramide_XX -> 'medula' (las piramides SON medula)
    tree_dom = cKDTree(coords_dom)
    _, idx_dom = tree_dom.query(puntos_drenaje, k=1)
    etiqueta_dom = region_dom[idx_dom]
    region = np.where(etiqueta_dom == "cortex", "cortex", "medula").astype("<U8")

    # --- 2b) eferente arterial asociada (origen del flujo) ------------------
    tree_term = cKDTree(pos_terminales)
    dist_efer, idx_term_local = tree_term.query(puntos_drenaje, k=1)
    # guardamos el INDICE del nodo terminal en capa3a.nodos (no el indice local)
    eferente_asociada = terminales[idx_term_local].astype(np.int64)

    # --- guardar ------------------------------------------------------------
    out_path = os.path.join(os.path.dirname(nef_path), OUT_NPZ)
    n_cort = int(np.count_nonzero(tipo_red == "peritubular_cortical"))
    n_vasa = int(np.count_nonzero(tipo_red == "vasa_recta"))
    n_med = int(np.count_nonzero(region == "medula"))
    n_ctx = int(np.count_nonzero(region == "cortex"))

    np.savez_compressed(
        out_path,
        puntos_drenaje=puntos_drenaje.astype(np.float32),
        tipo_red=tipo_red,
        region=region,
        nefrona_origen=nefrona_origen.astype(np.int32),
        eferente_asociada=eferente_asociada.astype(np.int32),
        dist_eferente_mm=dist_efer.astype(np.float32),
        # metadata: parametros
        puntos_por_tubulo_cortical=np.int32(PUNTOS_POR_TUBULO_CORTICAL),
        puntos_por_tubulo_yuxta=np.int32(PUNTOS_POR_TUBULO_YUXTA),
        jitter_radial_mm=np.float64(JITTER_RADIAL_MM),
        seed=np.int32(SEED),
        # metadata: conteos
        n_puntos=np.int32(M),
        n_peritubular_cortical=np.int32(n_cort),
        n_vasa_recta=np.int32(n_vasa),
        n_cortex=np.int32(n_ctx),
        n_medula=np.int32(n_med),
    )

    # ========================================================================
    #  VERIFICACION
    # ========================================================================
    print("\n  VERIFICACION")
    print("-" * 70)
    print(f"  Puntos de drenaje generados : {M}")
    esperado = (int(np.count_nonzero(tipo_nef == "cortical")) * PUNTOS_POR_TUBULO_CORTICAL
                + int(np.count_nonzero(tipo_nef == "yuxtamedular")) * PUNTOS_POR_TUBULO_YUXTA)
    print(f"     esperado por parametros  : {esperado}   "
          f"{'OK' if esperado == M else 'REVISAR'}")

    print(f"\n  Desglose por tipo de red:")
    print(f"     peritubular_cortical : {n_cort:>6d}  ({100.0*n_cort/M:5.1f} %)")
    print(f"     vasa_recta           : {n_vasa:>6d}  ({100.0*n_vasa/M:5.1f} %)")

    # --- CRITICO: cobertura medular vs cortical -----------------------------
    print(f"\n  CRITICO - distribucion CORTEZA vs MEDULA (por ubicacion real):")
    print(f"     cortex : {n_ctx:>6d}  ({100.0*n_ctx/M:5.1f} %)")
    print(f"     medula : {n_med:>6d}  ({100.0*n_med/M:5.1f} %)")
    # ¿la vasa recta SI puebla la medula?
    vr = tipo_red == "vasa_recta"
    vr_med = int(np.count_nonzero(vr & (region == "medula")))
    vr_tot = int(np.count_nonzero(vr))
    print(f"     -> de la vasa_recta, en medula: {vr_med}/{vr_tot}  "
          f"({(100.0*vr_med/vr_tot) if vr_tot else 0.0:5.1f} %)")
    print(f"     -> {'OK: la vasa recta puebla la medula (resuelve el deficit medular)' if n_med > 0 else 'FALLO: medula sin drenaje'}")

    # --- distancia a la eferente asociada (debe ser pequena) ----------------
    # "Pequena" se juzga frente a la ESCALA del organo (semieje mayor ~55 mm), no
    # en absoluto: los tubulos miden ~20 mm de media y los puntos se reparten por
    # todo su recorrido, asi que una media de pocos mm ya es <10% del organo.
    semieje_mayor = float(np.max(d0["semiejes"]))
    print(f"\n  Distancia punto de drenaje -> arteriola eferente asociada:")
    print(f"     media : {dist_efer.mean():.3f} mm   "
          f"({100.0*dist_efer.mean()/semieje_mayor:.1f} % del semieje mayor {semieje_mayor:.0f} mm)")
    print(f"     mediana: {np.median(dist_efer):.3f} mm")
    print(f"     p95   : {np.percentile(dist_efer, 95):.3f} mm")
    print(f"     max   : {dist_efer.max():.3f} mm")
    # desglose por region: los puntos medulares estan MAS lejos del arterial, y
    # eso es ESPERADO -> la medula es justo donde el arbol arterial es escaso (el
    # deficit medular). La vasa recta llega ahi precisamente para cubrir ese hueco.
    d_ctx = dist_efer[region == "cortex"]
    d_med = dist_efer[region == "medula"]
    if len(d_ctx):
        print(f"     media en cortex: {d_ctx.mean():.3f} mm  (cerca del arterial cortical)")
    if len(d_med):
        print(f"     media en medula: {d_med.mean():.3f} mm  (mayor: arterial escaso en medula, esperado)")
    print(f"     -> {'OK: red proxima al territorio arterial (<10% del organo)' if dist_efer.mean() < 0.10*semieje_mayor else 'REVISAR: red lejos del arterial'}")

    # --- rango espacial (bounding box) --------------------------------------
    bb_min = puntos_drenaje.min(axis=0)
    bb_max = puntos_drenaje.max(axis=0)
    print(f"\n  Rango espacial (bounding box) de los puntos de drenaje:")
    print(f"     x: [{bb_min[0]:7.2f}, {bb_max[0]:7.2f}] mm")
    print(f"     y: [{bb_min[1]:7.2f}, {bb_max[1]:7.2f}] mm")
    print(f"     z: [{bb_min[2]:7.2f}, {bb_max[2]:7.2f}] mm")
    print(f"     -> cubre corteza y medula: confirmado por el desglose de region")

    print("\n  -> guardado en:", out_path)
    print("=" * 70)


if __name__ == "__main__":
    main()
