#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa1_nefronas.py  --  Bio-Kidney AI 2026
==========================================
Capa 1 del gemelo digital: siembra de NEFRONAS sobre el dominio de la Capa 0.

Cada nefrona = un GLOMERULO (punto en la corteza) + un TUBULO (polilinea que
desciende desde el glomerulo hacia una piramide medular). Los glomerulos son
los puntos de atraccion anatomicos para el arbol vascular de capas futuras.

  NOTA DE ESCALA  --  Se siembran N_NEFRONAS = 1300 nefronas, NO el ~1.000.000
  del rinon humano real. Es una version REDUCIDA y REPRESENTATIVA, coherente
  con la demanda del generador CCO v8 previo. Densidades y proporciones
  (corticales/yuxtamedulares, profundidad, drenaje por piramide) se conservan;
  el conteo absoluto se escala. Documentado a proposito.

Solo NumPy / SciPy. SIN bpy. Entorno: env_biokidney.

Entrada : capa0_dominio.npz
Salida  : capa1_nefronas.npz
"""

import os
import numpy as np
from scipy.spatial import cKDTree

# ============================================================================
#  PARAMETROS  (editar aqui)
# ============================================================================
N_NEFRONAS = 1300          # version reducida representativa (real ~1e6)
FRAC_YUXTAMEDULAR = 0.15   # 15% yuxtamedulares, 85% corticales
SEED = 2026

DIST_MIN = 0.8             # distancia minima entre glomerulos [mm] (Poisson-disk)

# REANCLAJE A MM (jul 2026): con la Capa 0 corregida, la profundidad cortical viene en mm
# absolutos (depth_cortical_mm = distancia SOLO a la capsula). Los pesos y la banda
# yuxtamedular se anclan en mm, NO a la fraccion del viejo depth_norm (ya no es un
# normalizador estable: cambio de 28->52 mm al separar el seno; ver diagnostico).
# Preferencia de los glomerulos corticales por la mitad mas externa (depth_mm bajo):
CORTICAL_DEPTH_SCALE_MM = 3.3   # mm  escala de exp(-depth_mm/scale) (~= GROSOR/2, mitad externa)

# Banda yuxtamedular: tercio PROFUNDO del cortex, junto a la union cortico-medular.
#   juxta_min_depth_mm = GROSOR_CORTICAL_MM * JUXTA_BAND_FRAC   [mm]
JUXTA_BAND_FRAC = 0.75          # -> banda [0.75*GROSOR, GROSOR] mm (tercio profundo del cortex)

# Penetracion del tubulo en la piramide (fraccion de su longitud, base->apice)
#   0 = se queda en la union corticomedular (base);  1 = llega al apice (papila)
CORTICAL_PENETRACION = 0.25   # corticales: asa corta, medula externa
JUXTA_PENETRACION = 0.85      # yuxtamedulares: asa larga, hacia la papila
PENETRACION_JITTER = 0.08     # variabilidad aleatoria de la penetracion

TUBULE_NPTS = 6            # vertices por tubulo (version simple; se refinara)

IN_NPZ = "capa0_dominio.npz"
OUT_NPZ = "capa1_nefronas.npz"


# ============================================================================
#  UTILIDADES
# ============================================================================
def find_input():
    """Localiza capa0_dominio.npz (raiz del repo o junto al script)."""
    here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() \
        else os.getcwd()
    candidates = [
        os.path.expanduser(os.path.join("~", "Escritorio", "BioKidney-AI", IN_NPZ)),
        os.path.join(here, IN_NPZ),
        os.path.join(here, "..", IN_NPZ),
        os.path.join(os.getcwd(), IN_NPZ),
    ]
    for c in candidates:
        if os.path.isfile(c):
            return os.path.abspath(c)
    raise FileNotFoundError(f"No se encontro {IN_NPZ}. Buscado en: {candidates}")


def weighted_order(weights, rng):
    """Orden aleatorio SIN reemplazo ponderado por `weights`
    (Efraimidis-Spirakis): devuelve indices en orden de prioridad decreciente.
    """
    u = rng.random(len(weights))
    keys = u ** (1.0 / np.maximum(weights, 1e-12))
    return np.argsort(-keys)


def poisson_pick(coords_pool, order, n_target, accepted, count, dist_min, label):
    """Acepta puntos de `coords_pool` en el `order` dado respetando distancia
    minima `dist_min` contra TODO lo ya aceptado (array preasignado `accepted`,
    con `count` filas validas). Devuelve (count_nuevo, indices_locales, n_aceptados).
    """
    picked_local = []
    got = 0
    for li in order:
        if got >= n_target:
            break
        p = coords_pool[li]
        if count > 0:
            dmin = np.min(np.linalg.norm(accepted[:count] - p, axis=1))
            if dmin < dist_min:
                continue
        accepted[count] = p
        count += 1
        picked_local.append(li)
        got += 1
    if got < n_target:
        print(f"  [AVISO] solo se colocaron {got}/{n_target} glomerulos {label} "
              f"(candidatos/dist_min agotados)")
    return count, np.array(picked_local, dtype=np.int64), got


def assign_pyramids(glom, apex, axis, length):
    """Asigna cada glomerulo a la piramide cuyo EJE (segmento apice->base) esta
    mas cerca. Devuelve (pir_id (N,), base (P,3))."""
    base = apex + length[:, None] * axis
    n = len(glom)
    best_d = np.full(n, np.inf)
    best_i = np.full(n, -1, dtype=np.int32)
    for i in range(len(apex)):
        A, B = apex[i], base[i]
        AB = B - A
        denom = float(AB @ AB)
        t = np.clip(((glom - A) @ AB) / denom, 0.0, 1.0)
        proj = A + t[:, None] * AB
        d = np.linalg.norm(glom - proj, axis=1)
        better = d < best_d
        best_d[better] = d[better]
        best_i[better] = i
    return best_i, base


# ============================================================================
#  MAIN
# ============================================================================
def main():
    rng = np.random.default_rng(SEED)
    in_path = find_input()

    print("=" * 70)
    print("  CAPA 1 - SIEMBRA DE NEFRONAS  (glomerulos + tubulos)")
    print("=" * 70)
    print("  Entrada:", in_path)

    d0 = np.load(in_path, allow_pickle=False)
    coords = d0["coords"].astype(np.float64)
    region = d0["region_label"]
    depth = d0["depth"].astype(np.float64)                    # normalizada [0,1] (compat/reporte)
    depth_mm = d0["depth_cortical_mm"].astype(np.float64)     # ABSOLUTA [mm] (solo capsula) - criterio
    grosor = float(d0["grosor_cortical_mm"])                  # umbral cortico-medular [mm]
    umbral = float(d0["umbral_cm"])                           # fraccion equivalente (compat lectura)
    apex = d0["pyramid_apex"].astype(np.float64)
    axis = d0["pyramid_axis"].astype(np.float64)
    length = d0["pyramid_length"].astype(np.float64)
    n_pir = int(d0["n_piramides"])

    print(f"  N_NEFRONAS = {N_NEFRONAS}  (version reducida; real ~1e6)")
    print(f"  frac_yuxtamedular = {FRAC_YUXTAMEDULAR}  | dist_min = {DIST_MIN} mm")
    print(f"  grosor cortical = {grosor} mm (umbral CM en mm)  | banda yuxta = "
          f"[{grosor*JUXTA_BAND_FRAC:.2f}, {grosor:.2f}] mm (tercio profundo)")
    print("-" * 70)

    # --- pool de corteza ----------------------------------------------------
    cortex_idx = np.where(region == "cortex")[0]
    cx_coords = coords[cortex_idx]
    cx_depth = depth[cortex_idx]              # normalizada (reporte/guardado)
    cx_depth_mm = depth_mm[cortex_idx]        # ABSOLUTA [mm] (criterio de siembra)
    print(f"  Puntos de corteza disponibles: {len(cortex_idx)}")

    # cuotas
    n_juxta = int(round(FRAC_YUXTAMEDULAR * N_NEFRONAS))
    n_cortical = N_NEFRONAS - n_juxta

    # --- ordenes de prioridad ----------------------------------------------
    # Yuxtamedulares: tercio PROFUNDO del cortex (mm), junto a la union cortico-medular;
    # preferencia por el depth_mm mas alto (mas cerca del umbral de 6.6 mm).
    juxta_min_depth_mm = grosor * JUXTA_BAND_FRAC
    juxta_pool = np.where(cx_depth_mm >= juxta_min_depth_mm)[0]
    w_juxta = np.exp((cx_depth_mm[juxta_pool] - juxta_min_depth_mm) / CORTICAL_DEPTH_SCALE_MM)
    juxta_order = juxta_pool[weighted_order(w_juxta, rng)]

    # Corticales: toda la corteza, preferencia por la mitad mas externa (depth_mm bajo).
    w_cort = np.exp(-cx_depth_mm / CORTICAL_DEPTH_SCALE_MM)
    cortical_order = weighted_order(w_cort, rng)

    # --- siembra con distancia minima (Poisson-disk greedy) -----------------
    accepted = np.empty((N_NEFRONAS, 3), dtype=np.float64)
    count = 0
    # 1) yuxtamedulares primero (region mas restringida)
    count, juxta_local, got_j = poisson_pick(
        cx_coords, juxta_order, n_juxta, accepted, count, DIST_MIN, "yuxtamedulares")
    # 2) corticales despues (rellenan el resto), excluyendo solapes con los anteriores
    count, cort_local, got_c = poisson_pick(
        cx_coords, cortical_order, n_cortical, accepted, count, DIST_MIN, "corticales")

    n_total = count
    glomerulos = accepted[:n_total]

    # etiquetas, depth e indice-corteza por glomerulo (en orden de insercion)
    local_idx = np.concatenate([juxta_local, cort_local])
    tipo = np.array(["yuxtamedular"] * got_j + ["cortical"] * got_c, dtype="<U12")
    glom_depth = cx_depth[local_idx]              # normalizada [0,1] (compat/reporte)
    glom_depth_mm = cx_depth_mm[local_idx]        # ABSOLUTA [mm] (solo capsula)

    # --- destino piramidal + generacion de tubulos --------------------------
    pir_dest, base = assign_pyramids(glomerulos, apex, axis, length)

    penetracion = np.where(tipo == "yuxtamedular",
                           JUXTA_PENETRACION, CORTICAL_PENETRACION)
    penetracion = np.clip(
        penetracion + rng.uniform(-PENETRACION_JITTER, PENETRACION_JITTER, n_total),
        0.0, 1.0)

    tubulos = np.empty((n_total, TUBULE_NPTS, 3), dtype=np.float64)
    for k in range(n_total):
        i = pir_dest[k]
        A, B = apex[i], base[i]
        # objetivo dentro de la piramide: base (union) -> apice (papila)
        T = B + penetracion[k] * (A - B)
        tubulos[k] = np.linspace(glomerulos[k], T, TUBULE_NPTS)

    # --- guardar ------------------------------------------------------------
    out_path = os.path.join(os.path.dirname(in_path), OUT_NPZ)
    np.savez_compressed(
        out_path,
        glomerulos=glomerulos.astype(np.float32),
        tipo=tipo,
        tubulos=tubulos.astype(np.float32),       # (N, TUBULE_NPTS, 3)
        tubulos_npts=np.int32(TUBULE_NPTS),
        piramide_destino=pir_dest.astype(np.int32),
        depth_glomerulo=glom_depth.astype(np.float32),
        depth_glomerulo_mm=glom_depth_mm.astype(np.float32),   # profundidad cortical ABSOLUTA [mm]
        # metadata
        n_nefronas=np.int32(n_total),
        frac_yuxtamedular=np.float64(FRAC_YUXTAMEDULAR),
        seed=np.int32(SEED),
        distancia_minima=np.float64(DIST_MIN),
        cortical_penetracion=np.float64(CORTICAL_PENETRACION),
        juxta_penetracion=np.float64(JUXTA_PENETRACION),
        grosor_cortical_mm=np.float64(grosor),                 # umbral cortico-medular [mm]
        umbral_cm=np.float64(umbral),
        n_piramides=np.int32(n_pir),
    )

    # ========================================================================
    #  VERIFICACION
    # ========================================================================
    print("\n  VERIFICACION")
    print("-" * 70)

    n_c = int((tipo == "cortical").sum())
    n_j = int((tipo == "yuxtamedular").sum())
    print(f"  Total nefronas sembradas: {n_total}")
    print(f"     corticales    : {n_c:>5d}  ({100.0 * n_c / n_total:5.2f} %)  "
          f"[objetivo ~85%]")
    print(f"     yuxtamedulares: {n_j:>5d}  ({100.0 * n_j / n_total:5.2f} %)  "
          f"[objetivo ~15%]")

    # 0 glomerulos fuera de corteza (criterio en mm ABSOLUTOS: depth_mm < grosor)
    fuera_por_depth = int(np.count_nonzero(glom_depth_mm >= grosor))
    # confirmar que cada glomerulo coincide con un punto de corteza de Capa 0
    tree_cx = cKDTree(cx_coords)
    dd, _ = tree_cx.query(glomerulos, k=1)
    fuera_por_match = int(np.count_nonzero(dd > 1e-6))
    print(f"\n  Glomerulos fuera del cortex (depth_mm >= {grosor} mm): {fuera_por_depth}")
    print(f"  Glomerulos que no son punto-corteza de Capa 0       : {fuera_por_match}")
    print(f"  Profundidad cortical de glomerulos [mm]: min {glom_depth_mm.min():.2f}  "
          f"media {glom_depth_mm.mean():.2f}  max {glom_depth_mm.max():.2f}  (todos < {grosor})")
    print("  -> ambos deben ser 0")

    # rango de profundidad por tipo (yuxtamedulares deben tener depth mayor)
    dc = glom_depth[tipo == "cortical"]
    dj = glom_depth[tipo == "yuxtamedular"]
    print("\n  Profundidad (depth normalizada) de glomerulos:")
    print(f"     corticales    : min={dc.min():.3f}  media={dc.mean():.3f}  "
          f"max={dc.max():.3f}")
    print(f"     yuxtamedulares: min={dj.min():.3f}  media={dj.mean():.3f}  "
          f"max={dj.max():.3f}")
    print(f"     -> media yuxtamedular ({dj.mean():.3f}) > "
          f"media cortical ({dc.mean():.3f}): "
          f"{'OK' if dj.mean() > dc.mean() else 'REVISAR'}")

    # distribucion por piramide de destino
    print("\n  Nefronas por piramide de destino:")
    for i in range(n_pir):
        cnt = int(np.count_nonzero(pir_dest == i))
        print(f"     piramide_{i:02d} : {cnt:>4d}  ({100.0 * cnt / n_total:5.2f} %)")

    bbmin = glomerulos.min(axis=0)
    bbmax = glomerulos.max(axis=0)
    print("\n  Bounding box glomerulos [mm]:")
    print(f"     X:[{bbmin[0]:7.2f},{bbmax[0]:7.2f}]  "
          f"Y:[{bbmin[1]:7.2f},{bbmax[1]:7.2f}]  "
          f"Z:[{bbmin[2]:7.2f},{bbmax[2]:7.2f}]")

    print("\n  -> guardado en:", out_path)
    print("=" * 70)


if __name__ == "__main__":
    main()
