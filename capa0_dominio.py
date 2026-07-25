#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa0_dominio.py  --  Bio-Kidney AI 2026
=========================================
Capa 0 del pipeline: definicion del DOMINIO geometrico del rinon.

Construye un dominio renal con FORMA DE FRIJOL (no un elipsoide convexo):
un elipsoide principal al que se le RESTA un elipsoide de exclusion en la
cara medial, esculpiendo la hendidura del hilio / seno renal.

  PARENQUIMA = dentro(elipsoide principal) Y fuera(elipsoide del seno)

Define un campo de profundidad corteza->medula que trata tanto la capsula
externa COMO la pared concava del seno como "superficie" (los puntos
cercanos a cualquiera de las dos son corticales/superficiales), particiona
el volumen en corteza / medula y construye las piramides medulares
(conos cuyo apice -papila- APUNTA HACIA EL SENO, donde drenan los calices).

Solo NumPy / SciPy. SIN bpy (no depende de Blender).
Entorno: env_biokidney

Salida:
    capa0_dominio.npz
        coords        (N,3) float32   coordenadas de los puntos de parenquima
        region_label  (N,)  str       'cortex' | 'medulla' | 'piramide_XX'
        depth         (N,)  float32   profundidad normalizada [0,1]
        + metadata del dominio (semiejes, hilio, seno, umbral, n_piramides...)
"""

import numpy as np

# ============================================================================
#  PARAMETROS DEL DOMINIO  (editar aqui)
# ============================================================================
A_SEMI = 55.0          # semieje X  [mm]  -> eje largo (superior-inferior)
B_SEMI = 30.0          # semieje Y  [mm]  -> medial-lateral
C_SEMI = 18.0          # semieje Z  [mm]  -> anterior-posterior
SEMIEJES = np.array([A_SEMI, B_SEMI, C_SEMI], dtype=np.float64)

# Hilio / seno renal: cara MEDIAL del elipsoide -> punto (0, -b, 0)
HILIO = np.array([0.0, -B_SEMI, 0.0], dtype=np.float64)

# --- SENO RENAL (elipsoide de EXCLUSION que esculpe la hendidura medial) ---
# Centrado mas hacia afuera en -Y; al restarse del dominio crea la
# invaginacion del frijol. Ajustado para ocupar la zona media de la cara
# medial sin partir el organo en dos.
CENTRO_SENO = np.array([0.0, -34.0, 0.0], dtype=np.float64)
SEMIEJES_SENO = np.array([22.0, 16.0, 11.0], dtype=np.float64)

# --- UNION CORTICO-MEDULAR: umbral en MM ABSOLUTOS sobre la profundidad CORTICAL ---
# CORRECCION DE LOGICA RAIZ (jul 2026): la profundidad CORTICAL se mide SOLO contra la
# capsula externa. La pared del seno es superficie INTERNA (bajo ella hay medula / grasa
# sinusal, NO cortex): ya NO debe generar "cortex" peri-sinusal falso. El seno CONSERVA
# su rol geometrico (exclusion de parenquima e build_pyramids), que no usan este campo.
# El umbral es ABSOLUTO en mm (no una fraccion de depth_norm) para ser INVARIANTE ante
# cambios del normalizador (evita la trampa depth_norm 28->52mm detectada en diagnostico).
GROSOR_CORTICAL_MM = 6.6   # mm  espesor de la corteza renal (cortical width, MDCT n=2068)
#   depth_cortical_mm < GROSOR_CORTICAL_MM  -> corteza ;  >= -> medula
# UMBRAL_CM (fraccion) QUEDA DEPRECADO como criterio de corte: se conserva/guarda en el
# .npz solo como fraccion EQUIVALENTE (= GROSOR_CORTICAL_MM / depth_norm) por
# compatibilidad con lectores que hagan "depth < umbral_cm" sobre la profundidad normal.
UMBRAL_CM = 0.30           # [DEPRECADO] antigua fraccion; ya NO se usa para cortar cortex

N_PIRAMIDES = 10       # piramides medulares (rango humano fisiologico: 8-18)
CONE_HALF_ANGLE_DEG = 22.0   # semiangulo de apertura de cada cono piramidal

N_POINTS = 200_000     # numero de puntos VALIDOS de parenquima a generar
SEED = 2026            # semilla reproducible

OUT_NPZ = "capa0_dominio.npz"


# ============================================================================
#  GEOMETRIA: NIVELES ELIPSOIDALES Y PERTENENCIA AL PARENQUIMA
# ============================================================================
def ellipsoid_level(coords, center, semiejes):
    """Nivel radial normalizado de un elipsoide centrado en `center`.

        level = sqrt(sum(((p-center)/semiejes)^2))
        level < 1 -> dentro,  = 1 -> superficie,  > 1 -> fuera
    """
    norm = (coords - center) / semiejes
    return np.sqrt(np.sum(norm * norm, axis=-1))


def is_parenchyma(coords):
    """Mascara de parenquima: dentro del elipsoide principal Y fuera del seno."""
    inside_main = ellipsoid_level(coords, 0.0, SEMIEJES) < 1.0
    outside_seno = ellipsoid_level(coords, CENTRO_SENO, SEMIEJES_SENO) >= 1.0
    return inside_main & outside_seno


# ============================================================================
#  CAMPO DE PROFUNDIDAD CORTICAL  (SOLO capsula externa; el seno ya no cuenta)
# ============================================================================
def _surface_radius(coords, center, semiejes):
    """Radio del elipsoide a lo largo del rayo centro->punto (distancia del
    centro a la superficie en esa direccion). Estable para r>0."""
    q = coords - center
    r = np.linalg.norm(q, axis=-1)
    r_safe = np.where(r > 1e-12, r, 1.0)
    unit = q / r_safe[:, None]
    inv = np.sqrt(np.sum((unit / semiejes) ** 2, axis=-1))
    return r, 1.0 / inv  # (distancia centro->punto, radio de superficie)


def nearest_surface_distance(coords):
    """[CONSERVADA - ya NO define la profundidad cortical]  Distancia a la frontera
    MAS CERCANA del parenquima: capsula externa O pared del seno (el minimo).

    Se conserva por si algun consumidor la necesita, pero la CORRECCION DE LOGICA
    RAIZ (jul 2026) la retira del calculo de la profundidad cortical: incluir la
    pared del seno generaba "cortex" peri-sinusal falso (puntos interiores marcados
    corticales por cercania al seno). La etiqueta cortex/medula usa ahora
    capsule_distance() (solo capsula externa).
    """
    r_main, rsurf_main = _surface_radius(coords, np.zeros(3), SEMIEJES)
    dist_main = rsurf_main - r_main
    r_seno, rsurf_seno = _surface_radius(coords, CENTRO_SENO, SEMIEJES_SENO)
    dist_seno = r_seno - rsurf_seno
    return np.minimum(dist_main, dist_seno)


def capsule_distance(coords):
    """PROFUNDIDAD CORTICAL en mm: distancia (a lo largo del rayo) a la CAPSULA
    EXTERNA (elipsoide principal), positiva dentro del parenquima.

    SOLO la capsula externa define cuan 'cortical' es un punto. La pared del seno
    es superficie INTERNA (bajo ella hay medula / grasa sinusal): NO la contamos.
    """
    r_main, rsurf_main = _surface_radius(coords, np.zeros(3), SEMIEJES)
    return np.clip(rsurf_main - r_main, 0.0, None)


def compute_depth(coords):
    """Profundidad CORTICAL: en mm absolutos y su version normalizada [0,1].

    depth_mm = capsule_distance(coords)  (SOLO capsula externa).
      0 mm  = sobre la capsula externa -> maxima corticalidad.
      crece hacia el interior del parenquima.
    depth_norm  = max(depth_mm) [mm]  (normalizador).
    depth_normalizada = depth_mm / depth_norm  en [0,1]  (para visualizacion).

    La etiqueta cortex/medula se decide con umbral ABSOLUTO en mm
    (GROSOR_CORTICAL_MM), NO con la fraccion normalizada -> invariante al
    normalizador (no reaparece la trampa depth_norm 28->52mm del diagnostico).
    Devuelve (depth_mm, depth_norm, depth_normalizada).
    """
    depth_mm = capsule_distance(coords)
    depth_norm = float(depth_mm.max()) if depth_mm.size else 1.0
    depth_normalizada = np.clip(depth_mm / depth_norm, 0.0, 1.0)
    return depth_mm, depth_norm, depth_normalizada


# ============================================================================
#  MUESTREO DEL PARENQUIMA  (con rechazo del seno)
# ============================================================================
def sample_parenchyma(n, rng):
    """Genera n puntos uniformes en el PARENQUIMA (dentro del elipsoide
    principal y fuera del seno) por muestreo con rechazo.

    Devuelve (coords (n,3), n_tried, n_accept) donde la fraccion de
    aceptacion estima el volumen relativo del parenquima.
    """
    out = []
    n_tried = 0
    n_have = 0
    batch = max(n, 50_000)
    while n_have < n:
        # muestreo uniforme dentro del elipsoide principal (metodo de la bola)
        dirs = rng.normal(size=(batch, 3))
        dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
        radii = rng.uniform(0.0, 1.0, size=(batch, 1)) ** (1.0 / 3.0)
        pts = (dirs * radii) * SEMIEJES
        n_tried += batch

        # rechazar los que caen en el seno
        keep = ellipsoid_level(pts, CENTRO_SENO, SEMIEJES_SENO) >= 1.0
        pts = pts[keep]
        out.append(pts)
        n_have += len(pts)

    coords = np.concatenate(out, axis=0)[:n]
    n_accept = sum(len(p) for p in out)  # aceptados en el total muestreado
    return coords, n_tried, n_accept


# ============================================================================
#  PIRAMIDES MEDULARES (conos con apice hacia el seno)
# ============================================================================
def build_pyramids(n_piramides=N_PIRAMIDES):
    """Define la geometria de las piramides medulares.

    Cada piramide es un cono con:
        - APICE (papila) sobre la pared del seno, APUNTANDO HACIA LA CAVIDAD
          (hacia el seno), que es donde drena hacia los calices.
        - BASE en la union corticomedular, abriendose hacia la corteza.

    Los apices se situan en la cara +Y del elipsoide del seno (la pared que
    bordea el parenquima), repartidos a lo largo del eje largo (X) y en dos
    hileras (anterior/posterior, +-Z), imitando la disposicion humana.

    Devuelve apex (n,3), axis (n,3 unit), length (n,), r_base (n,).
    """
    a, b, c = SEMIEJES
    sx, sy, sz = SEMIEJES_SENO
    cx, cy, cz = CENTRO_SENO
    half = np.deg2rad(CONE_HALF_ANGLE_DEG)

    apex = np.zeros((n_piramides, 3))
    base = np.zeros((n_piramides, 3))

    for i in range(n_piramides):
        f = (i + 0.5) / n_piramides           # (0,1)
        xf = (f - 0.5) * 2.0                   # (-1,1) reparto a lo largo de X
        zrow = 1.0 if (i % 2 == 0) else -1.0   # hilera anterior / posterior

        # Coordenadas normalizadas en la seccion transversal del seno
        ux = xf * 0.60
        uz = zrow * 0.40
        rad2 = ux * ux + uz * uz               # < 1 para estar en la pared +Y

        # Apice sobre la pared +Y del seno (cara que mira al parenquima),
        # es decir el punto del elipsoide del seno con la mayor Y.
        uy = np.sqrt(max(0.0, 1.0 - rad2))
        apex[i] = np.array([cx + ux * sx,
                            cy + uy * sy,        # cara superior (+Y) del seno
                            cz + uz * sz])

        # Base hacia la union corticomedular (mas lateral, +Y y abierta en X,Z)
        base[i] = np.array([xf * 0.50 * a,
                            +0.40 * b,
                            zrow * 0.62 * c])

    vec = base - apex                          # apice -> base (apice mira al seno)
    length = np.linalg.norm(vec, axis=1)
    axis = vec / length[:, None]
    r_base = length * np.tan(half)
    return apex, axis, length, r_base


def assign_pyramids(coords, medulla_mask, apex, axis, length, r_base):
    """Asigna cada punto de la medula a una piramide (o a 'medula' generica).

    Un punto pertenece al cono i si su proyeccion axial t cae en [0, L_i]
    y su distancia radial al eje es menor que el radio del cono en t
    (radio crece linealmente de 0 en el apice a r_base en la base).
    Si cae en varios conos se asigna al de mejor ajuste (menor radio
    normalizado).
    """
    n = len(coords)
    best_score = np.full(n, np.inf)
    best_pyr = np.full(n, -1, dtype=np.int32)

    for i in range(len(apex)):
        v = coords - apex[i]
        t = v @ axis[i]                        # proyeccion axial
        perp = v - np.outer(t, axis[i])
        radial = np.linalg.norm(perp, axis=1)
        cone_r = (t / length[i]) * r_base[i]   # radio del cono en t
        inside = (medulla_mask & (t >= 0.0) & (t <= length[i])
                  & (radial <= cone_r))
        score = radial / np.maximum(cone_r, 1e-9)
        better = inside & (score < best_score)
        best_score[better] = score[better]
        best_pyr[better] = i

    return best_pyr


# ============================================================================
#  MAIN
# ============================================================================
def main():
    rng = np.random.default_rng(SEED)

    print("=" * 70)
    print("  CAPA 0 - DOMINIO RENAL FRIJOL  (elipsoide - seno + piramides)")
    print("=" * 70)
    print(f"  Semiejes principal (a,b,c) : {SEMIEJES} mm   (eje largo = X)")
    print(f"  Hilio / seno               : {HILIO}")
    print(f"  Seno: centro {CENTRO_SENO}  semiejes {SEMIEJES_SENO}")
    print(f"  Grosor cortical (umbral CM): {GROSOR_CORTICAL_MM} mm ABSOLUTO "
          f"(solo capsula; seno NO genera cortex)")
    print(f"  Piramides medulares        : {N_PIRAMIDES}")
    print(f"  Puntos de parenquima       : {N_POINTS}")
    print("-" * 70)

    # 1) Muestreo del parenquima (rechazando el seno)
    coords, n_tried, n_accept = sample_parenchyma(N_POINTS, rng)

    # 2) Campo de profundidad CORTICAL (SOLO capsula externa; el seno ya no cuenta)
    depth_mm, depth_norm, d = compute_depth(coords)

    # 3) Particion corteza / medula por umbral ABSOLUTO en mm (invariante al normalizador)
    umbral_cm_equiv = GROSOR_CORTICAL_MM / depth_norm   # fraccion equivalente (compat lectura)
    cortex_mask = depth_mm < GROSOR_CORTICAL_MM
    medulla_mask = ~cortex_mask

    # 4) Piramides medulares (apice hacia el seno)
    apex, axis, length, r_base = build_pyramids(N_PIRAMIDES)
    pyr_id = assign_pyramids(coords, medulla_mask, apex, axis, length, r_base)

    # 5) Etiquetas
    region_label = np.empty(N_POINTS, dtype="<U12")
    region_label[cortex_mask] = "cortex"
    region_label[medulla_mask] = "medulla"
    in_pyr = pyr_id >= 0
    for i in range(N_PIRAMIDES):
        region_label[pyr_id == i] = f"piramide_{i:02d}"

    # 6) Volumenes (Monte Carlo)
    vol_main = 4.0 / 3.0 * np.pi * np.prod(SEMIEJES)
    accept_frac = n_accept / n_tried           # parenquima / elipsoide principal
    vol_parenquima = vol_main * accept_frac
    vol_seno_excluido = vol_main - vol_parenquima

    # 7) Guardar
    np.savez_compressed(
        OUT_NPZ,
        coords=coords.astype(np.float32),
        region_label=region_label,
        depth=d.astype(np.float32),                    # profundidad cortical NORMALIZADA [0,1]
        depth_cortical_mm=depth_mm.astype(np.float32), # profundidad cortical ABSOLUTA [mm] (solo capsula)
        # --- metadata del dominio ---
        semiejes=SEMIEJES,
        hilio=HILIO,
        centro_seno=CENTRO_SENO,
        semiejes_seno=SEMIEJES_SENO,
        grosor_cortical_mm=np.float64(GROSOR_CORTICAL_MM),  # umbral cortico-medular ABSOLUTO [mm]
        umbral_cm=np.float64(umbral_cm_equiv),         # fraccion EQUIVALENTE (compat: depth<umbral_cm)
        depth_norm=np.float64(depth_norm),
        n_piramides=np.int32(N_PIRAMIDES),
        cone_half_angle_deg=np.float64(CONE_HALF_ANGLE_DEG),
        n_points=np.int32(N_POINTS),
        seed=np.int32(SEED),
        vol_parenquima_mm3=np.float64(vol_parenquima),
        vol_elipsoide_mm3=np.float64(vol_main),
        pyramid_apex=apex.astype(np.float64),
        pyramid_axis=axis.astype(np.float64),
        pyramid_length=length.astype(np.float64),
        pyramid_r_base=r_base.astype(np.float64),
    )

    # 8) Resumen
    n_cortex = int(cortex_mask.sum())
    n_medulla = int(medulla_mask.sum())
    n_in_pyr = int(in_pyr.sum())
    n_med_generic = n_medulla - n_in_pyr

    print("\n  VOLUMENES")
    print("-" * 70)
    print(f"  Elipsoide principal : {vol_main:11.1f} mm^3")
    print(f"  Parenquima (frijol) : {vol_parenquima:11.1f} mm^3  "
          f"({100.0 * accept_frac:5.2f} % del elipsoide)")
    print(f"  Seno excluido       : {vol_seno_excluido:11.1f} mm^3  "
          f"({100.0 * (1 - accept_frac):5.2f} % del elipsoide)")

    print("\n  PARTICION  (sobre el parenquima)")
    print("-" * 70)
    print(f"  Corteza : {n_cortex:>8d}  ({100.0 * n_cortex / N_POINTS:6.2f} %)")
    print(f"  Medula  : {n_medulla:>8d}  ({100.0 * n_medulla / N_POINTS:6.2f} %)")
    print(f"     - en piramides   : {n_in_pyr:>8d}  "
          f"({100.0 * n_in_pyr / N_POINTS:6.2f} %)")
    print(f"     - medula generica: {n_med_generic:>8d}  "
          f"({100.0 * n_med_generic / N_POINTS:6.2f} %)")

    print("\n  Conteo por piramide:")
    for i in range(N_PIRAMIDES):
        cnt = int((pyr_id == i).sum())
        print(f"     piramide_{i:02d} : {cnt:>7d}  "
              f"({100.0 * cnt / N_POINTS:5.2f} %)")

    bb_min = coords.min(axis=0)
    bb_max = coords.max(axis=0)
    print("\n  Bounding box [mm]:")
    print(f"     X: [{bb_min[0]:8.3f}, {bb_max[0]:8.3f}]")
    print(f"     Y: [{bb_min[1]:8.3f}, {bb_max[1]:8.3f}]")
    print(f"     Z: [{bb_min[2]:8.3f}, {bb_max[2]:8.3f}]")
    print(f"     depth cortical [mm]: [{depth_mm.min():.3f}, {depth_mm.max():.3f}]  "
          f"(norm={depth_norm:.3f} mm)   grosor cortical={GROSOR_CORTICAL_MM} mm")
    print(f"     cortex mas profundo (capsula): {depth_mm[cortex_mask].max():.3f} mm "
          f"(debe ser < {GROSOR_CORTICAL_MM} mm)")

    # --- INVARIANTE DE LA CORRECCION: cortex ya NO peri-sinusal (nearest surface = seno) ---
    r_main2, rsurf_main2 = _surface_radius(coords, np.zeros(3), SEMIEJES)
    d_cap = rsurf_main2 - r_main2
    r_seno2, rsurf_seno2 = _surface_radius(coords, CENTRO_SENO, SEMIEJES_SENO)
    d_sen = r_seno2 - rsurf_seno2
    peri_cortex = int(np.count_nonzero(cortex_mask & (d_sen < d_cap)))
    print(f"\n  INVARIANTE seno: puntos cortex cuya superficie mas cercana es el SENO: "
          f"{peri_cortex}  ({100.0*peri_cortex/max(1,n_cortex):.3f} % del cortex)")
    print(f"     (todo cortex esta < {GROSOR_CORTICAL_MM} mm de la CAPSULA; los residuales son")
    print(f"      esquina capsula/seno, corticales legitimos, NO interior profundo)")

    print("\n  -> guardado en:", OUT_NPZ)
    print("=" * 70)


if __name__ == "__main__":
    main()
