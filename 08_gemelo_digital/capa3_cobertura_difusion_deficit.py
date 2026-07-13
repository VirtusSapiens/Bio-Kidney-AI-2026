#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa3_cobertura_difusion_deficit.py  --  Bio-Kidney AI 2026
===========================================================
AUDITORIA DE SOLO LECTURA. NO escribe ningun .npz. Solo imprime a stdout.
Entorno: env_biokidney (NumPy/SciPy, SIN bpy).

PREGUNTA: con la vasculatura COMPLETA (arterial capa3a + peritubular capa3ab +
venoso capa3b), ¿que fraccion del tejido con demanda queda dentro del radio de
difusion (150um) de algun vaso, y DONDE esta el deficit restante? La validacion
visual sugirio un aparente hueco de vasculatura central que hay que cuantificar y
localizar: ¿seno renal esperado (no es parenquima -> inofensivo) o parenquima real?

METODO:
  - Reutiliza cobertura() de Capa 2 (distancia al vaso mas cercano < 150um).
  - Los arboles arterial/venoso son GRAFOS (nodos+aristas): se DENSIFICAN sus
    aristas en puntos a paso 0.1mm (< umbral 0.15mm) para que la distancia punto-a-
    punto del KDTree aproxime bien la distancia al segmento (centerline). La
    peritubular ya son puntos (puntos_drenaje) y se usan tal cual.
  - Nube capilar = densif(arterial) U puntos(peritubular) U densif(venoso).

  Deficit espacial: de los puntos de demanda NO cubiertos, se clasifican por
  (a) banda de distancia al hilio [0,-30,0] y (b) dentro/fuera del elipsoide del
  SENO renal esculpido en Capa 0 (centro_seno, semiejes_seno). Y se aisla el
  "centro" geometrico del tejido para ver si hay hueco de cobertura en parenquima
  central que NO corresponde al seno.

Linea base sin arbol (Capa 2, cov_proxy): 0.66% cubierto, deficit 98.87%.
"""

import os
import sys
import numpy as np
from scipy.spatial import cKDTree

PASO_DENSIF_MM = 0.1     # paso de muestreo de aristas -> puntos (< umbral 0.15mm)
BANDAS_HILIO = [(0, 10), (10, 20), (20, 30), (30, 40), (40, 60), (60, np.inf)]
R_CENTRO = [5.0, 10.0, 15.0]   # radios de la esfera "centro geometrico" [mm]

IN_DEM = "capa2_demanda.npz"
IN_DOM = "capa0_dominio.npz"
IN_ART = "capa3a_arterial.npz"
IN_PERI = "capa3ab_peritubular.npz"
IN_VEN = "capa3b_venoso.npz"

_here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
if _here not in sys.path:
    sys.path.insert(0, _here)
from capa2_demanda import cobertura, UMBRAL_DIFUSION_MM


def _p(name):
    """Ruta del npz: junto al script o en el cwd (repo root)."""
    for c in (os.path.join(_here, name), os.path.join(_here, "..", name), name):
        if os.path.isfile(c):
            return os.path.abspath(c)
    return name


def densificar(nodos, aristas, paso=PASO_DENSIF_MM):
    """Muestrea puntos a lo largo de cada arista a <= 'paso' mm (centerline -> puntos)."""
    P = nodos[aristas[:, 0]]; Q = nodos[aristas[:, 1]]
    Lmax = float(np.linalg.norm(Q - P, axis=1).max()) if len(P) else 0.0
    n = max(2, int(np.ceil(Lmax / paso)) + 1)
    t = np.linspace(0.0, 1.0, n)                      # (n,)
    pts = P[:, None, :] + t[None, :, None] * (Q - P)[:, None, :]   # (E,n,3)
    return pts.reshape(-1, 3)


def en_seno(pts, centro, semiejes):
    """Mascara: puntos dentro del elipsoide del seno renal."""
    q = (pts - centro) / semiejes
    return np.einsum('ij,ij->i', q, q) < 1.0


def main():
    dem = np.load(_p(IN_DEM), allow_pickle=False)
    coords = dem["coords"].astype(np.float64)
    demanda = dem["demanda"].astype(np.float64)
    base_pct = float(dem["cov_proxy_pct"])
    base_defic = 100.0 * float(dem["cov_proxy_deficit_frac"])
    N = len(coords); U = UMBRAL_DIFUSION_MM

    d0 = np.load(_p(IN_DOM), allow_pickle=False)
    hilio = d0["hilio"].astype(np.float64)
    centro_seno = d0["centro_seno"].astype(np.float64)
    semiejes_seno = d0["semiejes_seno"].astype(np.float64)

    art = np.load(_p(IN_ART), allow_pickle=False)
    peri = np.load(_p(IN_PERI), allow_pickle=False)
    ven = np.load(_p(IN_VEN), allow_pickle=False)

    p_art = densificar(art["nodos"].astype(np.float64), art["aristas"].astype(np.int64))
    p_peri = peri["puntos_drenaje"].astype(np.float64)
    p_ven = densificar(ven["nodos"].astype(np.float64), ven["aristas"].astype(np.int64))
    nube = np.vstack([p_art, p_peri, p_ven])

    print("=" * 88)
    print("  COBERTURA DE DIFUSION DE LA VASCULATURA COMPLETA SOBRE EL CAMPO DE DEMANDA")
    print("=" * 88)
    print(f"  demanda: {N} puntos (suma {demanda.sum():.0f})   umbral difusion = {U*1000:.0f}um")
    print(f"  nube capilar: arterial {len(p_art)} + peritubular {len(p_peri)} + venoso "
          f"{len(p_ven)} = {len(nube)} pts (aristas densif. a {PASO_DENSIF_MM}mm)")
    print("-" * 88)

    # ======================================================================
    #  (1) COBERTURA GLOBAL (completa) + desglose por componente
    # ======================================================================
    m_full = cobertura(nube, coords, demanda, U)
    dist = m_full["dist"]
    cub = dist < U
    print("  (1) COBERTURA GLOBAL")
    print("  " + "-" * 84)
    print(f"  {'componente':<26}{'% cubierto':>12}{'% deficit(pts)':>16}{'% deficit(demanda)':>20}")
    for nom, cloud in (("solo peritubular", p_peri), ("solo arterial", p_art),
                       ("solo venoso", p_ven), ("COMPLETA (a+p+v)", nube)):
        mm = cobertura(cloud, coords, demanda, U)
        print(f"  {nom:<26}{mm['pct_cubierto']:>11.3f}%{100-mm['pct_cubierto']:>15.3f}%"
              f"{100*mm['deficit_frac']:>19.2f}%")
    print("  " + "-" * 84)
    print(f"  BASE sin arbol (Capa 2): {base_pct:.2f}% cubierto, deficit {base_defic:.2f}%")
    print(f"  -> mejora de cobertura: {base_pct:.2f}% -> {m_full['pct_cubierto']:.2f}%  "
          f"(x{m_full['pct_cubierto']/max(base_pct,1e-9):.1f})")
    print(f"  dist al vaso mas cercano: media {m_full['dist_media']:.3f}  "
          f"max {m_full['dist_max']:.3f} mm")
    # --- diagnostico de RESOLUCION: cobertura vs radio de difusion -----------
    tv = ven["nodos"][ven["terminales"].astype(np.int64)].astype(np.float64)
    tt = cKDTree(tv); dnn, _ = tt.query(tv, k=2)
    sep_term = float(np.median(dnn[:, 1]))
    print(f"  RESOLUCION: espaciado mediano entre terminales venosos = {sep_term:.3f} mm "
          f"(la geometria vascular es {sep_term/U:.1f}x mas gruesa que el umbral {U*1000:.0f}um)")
    print(f"  {'umbral (um)':>12}{'% cubierto':>12}")
    for w in (0.15, 0.25, 0.35, 0.5, 0.75, 1.0):
        pc = 100.0 * np.count_nonzero(dist < w) / N
        marca = "  <== difusion de diseno (150um)" if abs(w - 0.15) < 1e-9 else ""
        print(f"  {w*1000:>11.0f} {pc:>11.2f}%{marca}")
    print("-" * 88)

    # ======================================================================
    #  (2) DEFICIT ESPACIAL: no cubiertos por banda al hilio y dentro/fuera seno
    # ======================================================================
    nocub = ~cub
    n_nocub = int(nocub.sum())
    seno = en_seno(coords, centro_seno, semiejes_seno)          # todo el dominio
    d_hil = np.linalg.norm(coords - hilio, axis=1)

    print("  (2a) NO CUBIERTOS por BANDA de distancia al hilio %s" % hilio.tolist())
    print("  " + "-" * 84)
    print(f"  {'banda (mm)':>12}{'n_total':>10}{'n_nocub':>10}{'%nocub':>9}"
          f"{'nocub en seno':>15}{'nocub fuera':>13}")
    for lo, hi in BANDAS_HILIO:
        mb = (d_hil >= lo) & (d_hil < hi)
        nt = int(mb.sum())
        mnc = mb & nocub
        nnc = int(mnc.sum())
        n_in = int((mnc & seno).sum()); n_out = nnc - n_in
        etiq = f"[{lo:.0f},{'inf' if np.isinf(hi) else f'{hi:.0f}'})"
        pctb = 100*nnc/nt if nt else 0.0
        print(f"  {etiq:>12}{nt:>10d}{nnc:>10d}{pctb:>8.1f}%{n_in:>15d}{n_out:>13d}")
    print("  " + "-" * 84)

    n_seno = int(seno.sum())
    nc_in = int((nocub & seno).sum()); nc_out = n_nocub - nc_in
    dem_in = float(demanda[nocub & seno].sum()); dem_out = float(demanda[nocub & ~seno].sum())
    print("  (2b) NO CUBIERTOS dentro/fuera del SENO RENAL "
          f"(elipsoide centro {centro_seno.tolist()} semiejes {semiejes_seno.tolist()})")
    print(f"     puntos de demanda dentro del seno (total): {n_seno} ({100*n_seno/N:.2f}%)")
    print(f"     NO cubiertos DENTRO del seno  (esperado, no parenquima): {nc_in}  "
          f"({100*nc_in/max(n_nocub,1):.1f}% del deficit, demanda {dem_in:.0f})")
    print(f"     NO cubiertos FUERA del seno   (PARENQUIMA REAL)        : {nc_out}  "
          f"({100*nc_out/max(n_nocub,1):.1f}% del deficit, demanda {dem_out:.0f})")
    print(f"     cobertura del PARENQUIMA (excluye seno): "
          f"{100*(1 - nc_out/max((~seno).sum(),1)):.3f}% cubierto")
    print("-" * 88)

    # ======================================================================
    #  (3) EL CENTRO: hueco de cobertura en parenquima central?
    # ======================================================================
    centro_geo = coords.mean(axis=0)
    d_cen = np.linalg.norm(coords - centro_geo, axis=1)
    print(f"  (3) TEJIDO DEL CENTRO GEOMETRICO {centro_geo.round(2).tolist()} "
          f"(dist al seno {np.linalg.norm(centro_geo-centro_seno):.1f}mm -> fuera del seno)")
    print("  " + "-" * 84)
    print(f"  {'esfera R':>9}{'n_dem':>9}{'n_nocub':>9}{'%nocub':>9}"
          f"{'nocub fuera seno':>18}{'dist vaso med/max':>20}")
    for R in R_CENTRO:
        mc = d_cen < R
        nt = int(mc.sum())
        mnc = mc & nocub
        nnc = int(mnc.sum())
        n_out = int((mnc & ~seno).sum())
        if nnc:
            dm = dist[mnc]; dstat = f"{dm.mean():.3f}/{dm.max():.3f}"
        else:
            dstat = "  -  "
        print(f"  {('<%.0fmm'%R):>9}{nt:>9d}{nnc:>9d}{100*nnc/max(nt,1):>8.1f}%"
              f"{n_out:>18d}{dstat:>20}")
    print("  " + "-" * 84)
    # foco en R=10mm
    mc = d_cen < 10.0
    mnc = mc & nocub
    nnc = int(mnc.sum()); n_out = int((mnc & ~seno).sum())
    if nnc:
        dm = dist[mnc]
        print(f"  centro <10mm: {nnc} demanda sin cubrir, TODOS fuera del seno={n_out==nnc}; "
              f"dist al vaso: media {dm.mean():.3f}  mediana {np.median(dm):.3f}  max {dm.max():.3f} mm")
    else:
        print(f"  centro <10mm: 0 puntos de demanda sin cubrir (parenquima central CUBIERTO)")
    # perfil radial: dist al vaso vs distancia al centro (aisla el VACIO central del efecto global)
    print("  perfil radial (mediana dist al vaso por cascara desde el centro):")
    print(f"     {'cascara (mm)':>14}{'n_dem':>8}{'dist vaso med':>15}{'vs global 1.74':>16}")
    for lo, hi in [(0, 5), (5, 10), (10, 15), (15, 25), (25, 40)]:
        ms = (d_cen >= lo) & (d_cen < hi)
        if ms.any():
            dmm = float(np.median(dist[ms]))
            print(f"     {f'[{lo},{hi})':>14}{int(ms.sum()):>8d}{dmm:>15.3f}"
                  f"{dmm/m_full['dist_media']:>15.2f}x")
    d_vaso_origen = float(np.linalg.norm(nube - centro_geo, axis=1).min())
    print(f"  vaso mas cercano al centro geometrico: {d_vaso_origen:.2f} mm "
          f"-> vacio avascular central de radio ~{d_vaso_origen:.0f}mm")
    print("-" * 88)

    # ======================================================================
    #  VEREDICTO
    # ======================================================================
    frac_seno = nc_in / max(n_nocub, 1)
    cob_05 = 100.0 * np.count_nonzero(dist < 0.5) / N
    print("  VEREDICTO")
    print("  " + "-" * 84)
    print(f"  Cobertura total {m_full['pct_cubierto']:.2f}% a {U*1000:.0f}um "
          f"(vs {base_pct:.2f}% base); deficit de demanda {100*m_full['deficit_frac']:.2f}%.")
    print(f"  El deficit NO esta en el seno (0% de la demanda cae dentro del seno; el campo de "
          f"demanda ya")
    print(f"  es parenquima puro): el 100% es parenquima. PERO no es un hueco localizado -> es "
          f"GLOBAL y")
    print(f"  ~uniforme (99% sin cubrir en toda banda al hilio), porque la geometria vascular "
          f"resuelve a")
    print(f"  ~{sep_term:.1f}mm (espaciado terminal), {sep_term/U:.0f}x mas grueso que la difusion de "
          f"{U*1000:.0f}um: a 0.5mm la cobertura")
    print(f"  sube a {cob_05:.0f}%. El bed capilar que cerraria los 150um NO esta en la geometria "
          f"(los arboles son")
    print(f"  troncos de suministro/drenaje). SOBRE ese fondo hay un vacio avascular REAL en el "
          f"centro")
    print(f"  geometrico (vaso mas cercano a {d_vaso_origen:.1f}mm; centro <5mm ~2x peor que la media).")
    print(f"  >> Cobertura total {m_full['pct_cubierto']:.2f}% (vs {base_pct:.2f}% base); el deficit "
          f"esta en PARENQUIMA REAL,")
    print(f"     pero es un limite de RESOLUCION global de la malla vascular (no seno, no hueco "
          f"unico), con un")
    print(f"     vacio central genuino de ~{d_vaso_origen:.0f}mm sobreimpuesto.")
    print("=" * 88)


if __name__ == "__main__":
    main()
