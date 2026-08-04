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
# ANCLA: Glodny et al. 2009, MDCT n=2068 (NO Beland; ver correccion de atribucion).
# La fuente mide ESPESOR CORTICAL PERPENDICULAR a la capsula -> la magnitud que se
# compara contra este umbral debe ser tambien perpendicular (ver capsule_distance).
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


def capsule_distance_radial(coords):
    """[DEPRECADA - conservada para auditoria]  Antigua profundidad cortical:
    distancia A LO LARGO DEL RAYO centroide->punto hasta la capsula.

    NO es un espesor: mide una cuerda radial desde el centroide, que sobreestima
    el espesor perpendicular (mas cuanto mas oblicuo es el rayo a la normal de la
    superficie). Reemplazada por capsule_distance() en jul 2026 (ver nota alli).
    """
    r_main, rsurf_main = _surface_radius(coords, np.zeros(3), SEMIEJES)
    return np.clip(rsurf_main - r_main, 0.0, None)


def _nearest_point_ellipsoid(coords, semiejes, n_iter=100):
    """Punto MAS CERCANO de la superficie del elipsoide (centrado en el origen,
    semiejes `semiejes`) para cada punto de `coords`. Devuelve (n,3).

    Metodo: multiplicador de Lagrange. El pie de la perpendicular cumple

        x_i = a_i^2 p_i / (a_i^2 + lam)     con   F(lam) = sum_i (a_i p_i)^2
                                                          / (a_i^2 + lam)^2 - 1 = 0

    F es estrictamente decreciente en el intervalo valido -> raiz unica por biseccion.

    GUARDA NUMERICA (critica).  Se reparametriza  mu = lam + min(a_i^2)  y se
    bisecta en mu >= 0, NO en lam. Motivo: el limite inferior del bracket es
    -min(a_i^2) sobre TODOS los ejes, y en los puntos que caen sobre un plano
    coordenado (p_i = 0 en el eje corto) la raiz se DEGENERA justo en ese borde:
    ahi el pie de la perpendicular sale FUERA del plano (x_i != 0 aunque p_i = 0),
    y la formula de Lagrange se vuelve 0/0 en ese eje. Dos formas de equivocarse:

      (a) evaluar 0/0 como 0  -> devuelve un punto que NO esta en la superficie.
          Para (0,-18,0) daria 10.125 mm (el punto (0,-28.125,0), interior).
      (b) subir el bracket a -min(a_i^2) SOLO sobre los ejes con p_i != 0
          -> se salta la rama degenerada y converge a OTRA raiz de F, que es
          punto estacionario pero NO el minimo global. Para (0,-18,0) daria
          12.000 mm (el polo (0,-30,0)), en vez de los 11.906 correctos.

    Tratamiento correcto: si en mu = 0 no hay raiz (F(0) <= 0), el minimo ESTA en
    mu = 0 y las coordenadas de los ejes degenerados (a_i^2 = min, p_i = 0) se
    recuperan de la ECUACION DEL ELIPSOIDE (el residuo), no de la formula 0/0.
    """
    p = np.atleast_2d(np.asarray(coords, dtype=np.float64))
    a = np.asarray(semiejes, dtype=np.float64)
    a2 = a * a
    m = a2.min()
    d2 = a2 - m                      # >= 0 ; vale 0 en el/los eje(s) mas corto(s)

    q = np.abs(p)                    # simetria por octante; los signos se restauran al final
    num = (a * q) ** 2               # (n,3)  numerador de cada termino de F

    def F(mu):
        """F(mu) con la convencion 0/0 -> 0 y x>0 / 0 -> +inf (limite correcto)."""
        den = d2 + mu[:, None]       # = a_i^2 + lam
        out = np.zeros_like(num)
        pos = num > 0.0
        with np.errstate(divide="ignore", invalid="ignore"):
            out[pos] = num[pos] / (den[pos] ** 2)
        return out.sum(axis=1) - 1.0

    n = len(q)
    mu = np.zeros(n)
    hi = np.linalg.norm(a * q, axis=1) + 1.0     # F(hi) < 0 garantizado
    tiene_raiz = F(np.zeros(n)) > 0.0            # False -> minimo degenerado en mu = 0

    # --- biseccion vectorizada solo donde hay raiz (mu > 0) ---
    lo_b = np.zeros(n)
    hi_b = hi.copy()
    for _ in range(n_iter):
        mid = 0.5 * (lo_b + hi_b)
        f = F(mid)
        arriba = f > 0.0                          # la raiz esta por encima de mid
        lo_b = np.where(arriba, mid, lo_b)
        hi_b = np.where(arriba, hi_b, mid)
    mu = np.where(tiene_raiz, 0.5 * (lo_b + hi_b), 0.0)

    # --- pie de la perpendicular ---
    den = d2 + mu[:, None]
    x = np.zeros_like(q)
    ok = den > 0.0
    x[ok] = (np.broadcast_to(a2, q.shape)[ok] * q[ok]) / den[ok]

    # --- rama degenerada: ejes con a_i^2 = min y p_i = 0 -> del residuo del elipsoide ---
    deg = (d2 == 0.0)                             # eje(s) mas corto(s)
    if np.any(~tiene_raiz):
        sel = ~tiene_raiz
        resid = 1.0 - np.sum((x[sel] ** 2) / a2, axis=1)
        resid = np.clip(resid, 0.0, None)
        j = int(np.argmax(deg))                   # primer eje degenerado
        x[sel, j] = a[j] * np.sqrt(resid)

    return np.sign(np.where(p == 0.0, 1.0, p)) * x


def capsule_distance(coords):
    """PROFUNDIDAD CORTICAL en mm: distancia al PUNTO MAS CERCANO de la CAPSULA
    externa (superficie del elipsoide principal 55/30/18).

    CAMBIO DE MAGNITUD (jul 2026, metodo B).  Antes se media la distancia RADIAL
    desde el centroide (capsule_distance_radial, ahora deprecada). Esa magnitud no
    es un espesor: es una cuerda desde el centro, y sobreestima el espesor real en
    todo punto cuyo rayo no sea normal a la capsula. El umbral de corte
    cortico-medular (GROSOR_CORTICAL_MM = 6.6 mm, Glodny 2009, MDCT n=2068) esta
    definido sobre el espesor cortical PERPENDICULAR a la capsula, de modo que
    comparar 6.6 mm contra una distancia radial mezclaba dos magnitudes distintas.
    La distancia al punto mas cercano de la capsula es normal a la superficie por
    construccion, y por tanto APROXIMA el espesor perpendicular medido en MDCT.

    SOLO la capsula externa define cuan 'cortical' es un punto. La pared del seno
    es superficie INTERNA (bajo ella hay medula / grasa sinusal): NO la contamos.

    LIMITACION CONOCIDA (capsula fantasma): el elipsoide principal es una superficie
    CERRADA, pero la porcion de el que queda dentro del elipsoide del seno fue
    excavada y no existe como capsula real. Los puntos vecinos a la pared del seno
    pueden por tanto medir su distancia contra un tramo de capsula fantasma. El
    bloque VERIFICACION de main() cuantifica cuantos puntos estan afectados.
    """
    x = _nearest_point_ellipsoid(coords, SEMIEJES)
    return np.linalg.norm(np.atleast_2d(coords) - x, axis=-1)


def compute_depth(coords):
    """Profundidad CORTICAL: en mm absolutos y su version normalizada [0,1].

    depth_mm = capsule_distance(coords)  (SOLO capsula externa, distancia al punto
      MAS CERCANO de la capsula = espesor PERPENDICULAR aproximado; jul 2026).
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

    # ========================================================================
    #  VERIFICACION
    # ========================================================================
    print("\n  VERIFICACION")
    print("-" * 70)

    # --- (1) GUARDA NUMERICA: puntos SOBRE ejes coordenados (caso degenerado) ---
    # En estos puntos algun p_i = 0 y el pie de la perpendicular sale FUERA de ese
    # plano: es donde la formula de Lagrange se vuelve 0/0 y donde un bracket mal
    # elegido devuelve basura SILENCIOSA (ver docstring de _nearest_point_ellipsoid).
    # Valores esperados verificados de forma independiente (minimizacion directa
    # sobre la parametrizacion de la superficie).
    casos = [
        ((0.0, -18.0,   0.0), 11.905881, "eje Y - CASO CRITICO (0,-18,0)"),
        ((0.0,   0.0,   0.0), 18.000000, "centroide (todos los ejes degenerados)"),
        ((0.0, -29.0,   0.0),  1.000000, "eje Y, cerca de la capsula"),
        ((0.0,  25.0,   0.0),  5.000000, "eje Y, lado lateral"),
        ((0.0,   0.0, -15.0),  3.000000, "eje Z (eje corto = eje degenerado)"),
        ((50.0,  0.0,   0.0),  5.000000, "eje X, cerca del polo"),
        ((40.0,  0.0,   0.0), 11.492218, "eje X, pie perpendicular fuera del eje"),
    ]
    pts = np.array([c[0] for c in casos], dtype=np.float64)
    got = capsule_distance(pts)
    foot = _nearest_point_ellipsoid(pts, SEMIEJES)
    resid = np.abs(ellipsoid_level(foot, np.zeros(3), SEMIEJES) - 1.0)
    radial = capsule_distance_radial(pts)

    print("  Guarda numerica - distancia perpendicular en puntos sobre ejes coordenados:")
    print(f"     {'punto':>18s} {'esperado':>10s} {'obtenido':>10s} {'err':>9s} "
          f"{'|lvl-1|':>9s} {'radial(dep)':>11s}")
    ok_guarda = True
    for k, (p, exp, nota) in enumerate(casos):
        err = abs(got[k] - exp)
        bien = (err < 1e-5) and (resid[k] < 1e-9)
        ok_guarda &= bien
        print(f"     {str(p):>18s} {exp:10.6f} {got[k]:10.6f} {err:9.2e} "
              f"{resid[k]:9.2e} {radial[k]:11.3f}   "
              f"{'OK' if bien else 'FALLO'}  <- {nota}")
    print(f"     -> el pie de la perpendicular cae SOBRE la superficie "
          f"(|level-1| < 1e-9) en los {len(casos)} casos: "
          f"{'OK' if np.all(resid < 1e-9) else 'FALLO'}")
    print(f"     -> GUARDA NUMERICA: {'OK' if ok_guarda else 'FALLO'}")
    print("        (bracket global mal aplicado daria 10.125 mm en (0,-18,0);")
    print("         bracket solo sobre ejes con p_i != 0 daria 12.000 mm -> ambos MAL)")

    # --- (2) IMPACTO DEL CAMBIO DE MAGNITUD (radial -> perpendicular) ---
    d_rad = capsule_distance_radial(coords)
    cortex_rad = d_rad < GROSOR_CORTICAL_MM
    n_reetiq = int(np.count_nonzero(cortex_rad != cortex_mask))
    print("\n  Impacto del cambio de magnitud (radial deprecada -> perpendicular):")
    print(f"     cortex con magnitud RADIAL (deprecada) : {int(cortex_rad.sum()):>8d}  "
          f"({100.0 * cortex_rad.mean():6.2f} %)")
    print(f"     cortex con magnitud PERPENDICULAR      : {n_cortex:>8d}  "
          f"({100.0 * n_cortex / N_POINTS:6.2f} %)")
    print(f"     puntos REETIQUETADOS                   : {n_reetiq:>8d}  "
          f"({100.0 * n_reetiq / N_POINTS:6.2f} % del parenquima)")
    print(f"     depth perpendicular <= radial en todo punto (la radial "
          f"sobreestima): "
          f"{'OK' if np.all(depth_mm <= d_rad + 1e-6) else 'FALLO'}")

    # --- (3) LIMITACION: CAPSULA FANTASMA ---
    # El tramo del elipsoide principal excavado por el seno no existe como capsula
    # real. Un punto mide contra capsula fantasma si su pie de perpendicular cae
    # DENTRO del elipsoide del seno.
    foot_all = _nearest_point_ellipsoid(coords, SEMIEJES)
    fantasma = ellipsoid_level(foot_all, CENTRO_SENO, SEMIEJES_SENO) < 1.0
    n_fant = int(fantasma.sum())
    print("\n  Limitacion conocida - capsula fantasma (tramo excavado por el seno):")
    print(f"     puntos cuyo pie de perpendicular cae en el seno: {n_fant:>8d}  "
          f"({100.0 * n_fant / N_POINTS:.3f} % del parenquima)")
    if n_fant:
        # distancia PERPENDICULAR a la pared del seno (misma magnitud que depth)
        q_seno = coords[fantasma].astype(np.float64) - CENTRO_SENO
        foot_seno = _nearest_point_ellipsoid(q_seno, SEMIEJES_SENO)
        d_pared_seno = np.linalg.norm(q_seno - foot_seno, axis=1)
        print(f"     todos ellos a <= {d_pared_seno.max():.1f} mm de la pared del seno "
              f"(mediana {np.median(d_pared_seno):.2f} mm) -> banda delgada peri-sinusal")
        print(f"     de ellos etiquetados cortex: "
              f"{int(np.count_nonzero(fantasma & cortex_mask))}")

    # --- (4) SVD / RANGO DE LA NUBE REGENERADA ---
    centro = coords.mean(axis=0)
    sv = np.linalg.svd(coords - centro, compute_uv=False)
    rango = int(np.linalg.matrix_rank(coords - centro))
    print("\n  SVD / rango de la nube regenerada:")
    print(f"     centroide [mm]        : "
          f"[{centro[0]:8.4f}, {centro[1]:8.4f}, {centro[2]:8.4f}]")
    print(f"     valores singulares    : "
          f"[{sv[0]:.4f}, {sv[1]:.4f}, {sv[2]:.4f}]")
    print(f"     sv normalizados (/sv0): "
          f"[{sv[0]/sv[0]:.4f}, {sv[1]/sv[0]:.4f}, {sv[2]/sv[0]:.4f}]")
    print(f"     rango                 : {rango}  "
          f"{'OK (nube 3D, no degenerada)' if rango == 3 else 'FALLO'}")
    print(f"     desv. tipica por eje  : "
          f"[{coords[:,0].std():.4f}, {coords[:,1].std():.4f}, {coords[:,2].std():.4f}] mm")

    print("\n  -> guardado en:", OUT_NPZ)
    print("=" * 70)


if __name__ == "__main__":
    main()
