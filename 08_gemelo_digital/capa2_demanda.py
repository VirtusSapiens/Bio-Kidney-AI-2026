#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa2_demanda.py  --  Bio-Kidney AI 2026
=========================================
Capa 2 del gemelo digital: CAMPO ESCALAR DE DEMANDA METABOLICA de O2 sobre el
parenquima, en coordenadas reales (mm).

La demanda es "el problema" que el arbol vascular (Capa 3) debera resolver:
cada punto del tejido necesita estar a < UMBRAL_DIFUSION de un capilar. Este
script construye el mapa de cuanta O2 demanda cada region, concentrando la
demanda donde hay unidades funcionales (glomerulos + tubulos de la Capa 1).

Incluye una funcion cobertura() lista para la Capa 3: dado un array de
posiciones de capilares, mide que tan bien el tejido queda servido. Por ahora
se ejecuta con los GLOMERULOS como capilares-proxy (resultado malo esperado:
es el punto de partida que el arbol vendra a mejorar).

Solo NumPy / SciPy. SIN bpy. Entorno: env_biokidney.

Entrada : capa0_dominio.npz, capa1_nefronas.npz
Salida  : capa2_demanda.npz
"""

import os
import numpy as np
from scipy.spatial import cKDTree

# ============================================================================
#  PARAMETROS  (editar aqui)
# ============================================================================
# Limite fisiologico de difusion de O2 ~200 um (hipoxia). Disenamos con margen.
UMBRAL_DIFUSION_UM = 150        # distancia max de diseno tejido-capilar [micras]
UMBRAL_DIFUSION_MM = 0.150      # mismo valor en mm  (1 um = 1e-3 mm -> 150 um = 0.150 mm)

DEMANDA_CORTEX = 1.0            # demanda relativa alta (glomerulos, tub. contorneados)
DEMANDA_MEDULA = 0.4           # demanda relativa media-baja

RADIO_INFLUENCIA_NEFRONA_MM = 0.5   # radio de realce de demanda alrededor de nefronas
REALCE_PESO = 0.3              # peso del realce local frente a la demanda base

SEED = 2026

IN_DOM = "capa0_dominio.npz"
IN_NEF = "capa1_nefronas.npz"
OUT_NPZ = "capa2_demanda.npz"


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


# ============================================================================
#  CAMPO DE DEMANDA
# ============================================================================
def build_demand(coords, region_label, nefrona_pts):
    """Construye el campo de demanda en [0,1] sobre los puntos del dominio.

    1) Demanda base por region: cortex -> DEMANDA_CORTEX, resto -> DEMANDA_MEDULA.
    2) Realce local: kernel lineal decreciente alrededor de cada punto-nefrona
       dentro de RADIO_INFLUENCIA_NEFRONA_MM (concentra demanda en unidades activas).
    3) Normaliza min-max a [0,1].

    Devuelve (demanda_norm (N,), base (N,), realce (N,)).
    """
    n = len(coords)
    is_cortex = region_label == "cortex"

    # 1) base por region
    base = np.where(is_cortex, DEMANDA_CORTEX, DEMANDA_MEDULA).astype(np.float64)

    # 2) realce local por cercania de nefronas (scatter desde puntos-nefrona)
    realce = np.zeros(n, dtype=np.float64)
    r = RADIO_INFLUENCIA_NEFRONA_MM
    dom_tree = cKDTree(coords)
    vecinos = dom_tree.query_ball_point(nefrona_pts, r)   # lista por punto-nefrona
    for k, idxs in enumerate(vecinos):
        if not idxs:
            continue
        idxs = np.asarray(idxs)
        d = np.linalg.norm(coords[idxs] - nefrona_pts[k], axis=1)
        w = 1.0 - d / r                                  # kernel lineal (0..1)
        np.add.at(realce, idxs, w)

    # Normaliza el realce a [0,1] ANTES de sumarlo, usando el percentil 99 como
    # referencia (no el maximo): donde muchos tubulos convergen el realce crudo
    # se dispara y un solo outlier aplastaria todo el campo al normalizar.
    pos = realce[realce > 0]
    ref = np.percentile(pos, 99) if pos.size else 1.0
    realce_norm = np.clip(realce / ref, 0.0, 1.0) if ref > 0 else realce

    # 3) campo final + normalizacion DIVIDE-POR-MAXIMO (no min-max).
    #    Justificacion: la medula NO es avascular -- requiere vasa recta para
    #    irrigar las asas de Henle yuxtamedulares (sembradas en Capa 1). Con
    #    min-max la base medular (0.4) caeria al suelo del rango -> demanda = 0,
    #    y entonces el arbol de Capa 3 podria dejar la medula sin irrigar y aun
    #    asi reportar exito, lo cual seria anatomicamente incoherente. Dividir
    #    por el maximo conserva la demanda medular residual (~0.31): la medula
    #    demanda MENOS que la corteza, pero no cero.
    final = base + REALCE_PESO * realce_norm
    hi = float(final.max())
    demanda = final / hi if hi > 0 else np.zeros_like(final)
    return demanda, base, realce


# ============================================================================
#  COBERTURA  (lista para Capa 3; ejecutable ya, incluso con input vacio)
# ============================================================================
def cobertura(puntos_capilares, coords, demanda, umbral_mm=UMBRAL_DIFUSION_MM):
    """Dado un array (M,3) de posiciones de capilares [mm], mide la cobertura
    del tejido: para cada punto del dominio, distancia al capilar mas cercano.

    Devuelve un dict con:
        pct_cubierto   : % de puntos con dist < umbral
        dist_media     : distancia media al capilar mas cercano [mm]
        dist_max       : peor caso [mm]
        deficit        : demanda total de los puntos NO cubiertos (dist > umbral)
        deficit_frac   : deficit / demanda_total  (fraccion de demanda sin servir)
        dist           : (N,) distancias por punto [mm]
    """
    n = len(coords)
    cap = np.asarray(puntos_capilares, dtype=np.float64).reshape(-1, 3) \
        if puntos_capilares is not None else np.empty((0, 3))

    if len(cap) == 0:
        # sin capilares: nada cubierto, distancia infinita, deficit = toda la demanda
        return dict(pct_cubierto=0.0, dist_media=np.inf, dist_max=np.inf,
                    deficit=float(demanda.sum()), deficit_frac=1.0,
                    dist=np.full(n, np.inf))

    tree = cKDTree(cap)
    dist, _ = tree.query(coords, k=1)
    cubierto = dist < umbral_mm
    no_cub = ~cubierto
    dem_total = float(demanda.sum())
    deficit = float(demanda[no_cub].sum())
    return dict(
        pct_cubierto=100.0 * np.count_nonzero(cubierto) / n,
        dist_media=float(dist.mean()),
        dist_max=float(dist.max()),
        deficit=deficit,
        deficit_frac=(deficit / dem_total) if dem_total > 0 else 0.0,
        dist=dist,
    )


def print_cobertura(nombre, m):
    print(f"  Cobertura [{nombre}]:")
    print(f"     % tejido cubierto (< {UMBRAL_DIFUSION_MM} mm): {m['pct_cubierto']:.3f} %")
    print(f"     distancia media al capilar : {m['dist_media']:.3f} mm")
    print(f"     distancia maxima (peor caso): {m['dist_max']:.3f} mm")
    print(f"     deficit (demanda sin servir): {m['deficit']:.1f}  "
          f"({100.0 * m['deficit_frac']:.2f} % de la demanda total)")


# ============================================================================
#  MAIN
# ============================================================================
def main():
    np.random.default_rng(SEED)  # reproducibilidad (no hay aleatoriedad esencial aqui)

    dom_path = find_npz(IN_DOM)
    nef_path = find_npz(IN_NEF)

    print("=" * 70)
    print("  CAPA 2 - CAMPO DE DEMANDA METABOLICA DE O2")
    print("=" * 70)
    print("  Entrada dominio :", dom_path)
    print("  Entrada nefronas:", nef_path)
    print(f"  Umbral de difusion: {UMBRAL_DIFUSION_UM} um = {UMBRAL_DIFUSION_MM} mm "
          f"(margen bajo el limite de hipoxia ~200 um)")
    print(f"  Demanda base: cortex={DEMANDA_CORTEX}  medula/piramides={DEMANDA_MEDULA}")
    print(f"  Radio de influencia nefrona: {RADIO_INFLUENCIA_NEFRONA_MM} mm")
    print("-" * 70)

    d0 = np.load(dom_path, allow_pickle=False)
    coords = d0["coords"].astype(np.float64)
    region = d0["region_label"]
    n = len(coords)

    d1 = np.load(nef_path, allow_pickle=False)
    glomerulos = d1["glomerulos"].astype(np.float64)
    tubulos = d1["tubulos"].astype(np.float64)            # (M, NPTS, 3)
    tub_pts = tubulos.reshape(-1, 3)
    nefrona_pts = np.concatenate([glomerulos, tub_pts], axis=0)
    print(f"  Puntos de dominio: {n}")
    print(f"  Puntos-nefrona (glomerulos + vertices de tubulo): {len(nefrona_pts)}")

    # --- campo de demanda ---------------------------------------------------
    demanda, base, realce = build_demand(coords, region, nefrona_pts)

    # --- cobertura-proxy con glomerulos (punto de partida) ------------------
    cov_proxy = cobertura(glomerulos, coords, demanda)

    # --- guardar ------------------------------------------------------------
    out_path = os.path.join(os.path.dirname(dom_path), OUT_NPZ)
    np.savez_compressed(
        out_path,
        coords=coords.astype(np.float32),
        demanda=demanda.astype(np.float32),
        # metadata
        umbral_difusion_um=np.float64(UMBRAL_DIFUSION_UM),
        umbral_difusion_mm=np.float64(UMBRAL_DIFUSION_MM),
        demanda_cortex=np.float64(DEMANDA_CORTEX),
        demanda_medula=np.float64(DEMANDA_MEDULA),
        radio_influencia_mm=np.float64(RADIO_INFLUENCIA_NEFRONA_MM),
        realce_peso=np.float64(REALCE_PESO),
        seed=np.int32(SEED),
        # referencia: cobertura-proxy inicial con glomerulos
        cov_proxy_pct=np.float64(cov_proxy["pct_cubierto"]),
        cov_proxy_dist_media=np.float64(cov_proxy["dist_media"]),
        cov_proxy_dist_max=np.float64(cov_proxy["dist_max"]),
        cov_proxy_deficit=np.float64(cov_proxy["deficit"]),
        cov_proxy_deficit_frac=np.float64(cov_proxy["deficit_frac"]),
    )

    # ========================================================================
    #  VERIFICACION
    # ========================================================================
    is_cortex = region == "cortex"
    is_pir = np.char.startswith(region, "piramide")
    is_med = ~is_cortex & ~is_pir            # medula generica (columnas de Bertin)
    dc = demanda[is_cortex]
    dm = demanda[is_med]
    dp = demanda[is_pir]
    dmed_all = demanda[~is_cortex]           # toda la medula (generica + piramides)

    print("\n  VERIFICACION")
    print("-" * 70)
    print("  Demanda (normalizada, divide-por-maximo) por region:")
    print(f"     cortex          : n={dc.size:>7d}  min={dc.min():.3f}  "
          f"media={dc.mean():.3f}  max={dc.max():.3f}")
    print(f"     medula generica : n={dm.size:>7d}  min={dm.min():.3f}  "
          f"media={dm.mean():.3f}  max={dm.max():.3f}")
    print(f"     piramides       : n={dp.size:>7d}  min={dp.min():.3f}  "
          f"media={dp.mean():.3f}  max={dp.max():.3f}")
    print(f"     -> media cortex ({dc.mean():.3f}) > media medula "
          f"({dmed_all.mean():.3f}): "
          f"{'OK' if dc.mean() > dmed_all.mean() else 'REVISAR'}")
    print(f"     -> medula base conservada (~0.31, no aplastada a 0): "
          f"min medula = {dmed_all.min():.3f}")
    print(f"     puntos con realce local > 0: {int(np.count_nonzero(realce > 0))} "
          f"({100.0 * np.count_nonzero(realce > 0) / n:.1f} %)")

    # histograma de texto (deciles)
    print("\n  Histograma de demanda (deciles):")
    edges = np.linspace(0.0, 1.0, 11)
    hist, _ = np.histogram(demanda, bins=edges)
    hmax = hist.max() if hist.max() > 0 else 1
    for i in range(10):
        bar = "#" * int(40 * hist[i] / hmax)
        print(f"     [{edges[i]:.1f}-{edges[i+1]:.1f}) {hist[i]:>7d} |{bar}")

    # cobertura-proxy
    print("\n  COBERTURA DE PARTIDA (proxy: glomerulos como capilares)")
    print("-" * 70)
    print_cobertura("glomerulos-proxy", cov_proxy)
    print("  (cobertura BAJA esperada: es el punto de partida que la Capa 3,")
    print("   el arbol vascular, vendra a mejorar.)")

    print("\n  -> guardado en:", out_path)
    print("=" * 70)


if __name__ == "__main__":
    main()
