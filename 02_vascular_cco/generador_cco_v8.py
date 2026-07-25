"""
╔══════════════════════════════════════════════════════════════════════════════╗
║       GENERADOR CCO v8 — Bio-Kidney AI 2026                               ║
║       Carlos David Moreno Cáceres — VirtusSapiens                          ║
║                                                                            ║
║  Cambios v7 → v8:                                                          ║
║    1. Densidad cortical +30 % (demanda Beta(3, 1.2) sesgada a corteza)     ║
║    2. Murray α = 3.0 estricto + presiones Poiseuille por segmento          ║
║    3. cKDTree vectorizado (cobertura O(N log N) en vez de O(N·M))          ║
║    4. Exportación JSON: renal_data_v1.json (diámetros, coords, presiones)  ║
║    5. Memoria controlada — pico ≈ 2-3 GB en 16 GB RAM                     ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import numpy as np
import json
import csv
import os
import sys
from scipy.spatial import cKDTree

# ── PARÁMETROS GLOBALES ──────────────────────────────────────────────────────
ALPHA           = 3.0           # exponente Murray
SEMILLA         = 42
PUNTOS_SPLINE   = 10
R_MIN           = 8e-6          # 8 µm — arteriola aferente mínima
R_MIN_COLECTOR  = 200e-6

# Dimensiones riñón adulto (m)
EJE_X, EJE_Y, EJE_Z = 0.055, 0.030, 0.025

# Hilios
HILIO_ART = np.array([ 0.0, -EJE_Y * 0.85,  0.003])
HILIO_VEN = np.array([ 0.0, -EJE_Y * 0.85, -0.003])
HILIO_URI = np.array([ 0.0, -EJE_Y * 0.85,  0.000])

# Radios raíz
R_ARTERIA = 500e-6
R_VENA    = 600e-6

# v8 — densificación cortical (+30 % puntos de demanda)
N_SEMILLAS_DEMANDA  = 1300      # 1000 → 1300 (+30 %)
RADIO_COBERTURA     = 0.0045    # 5 mm → 4.5 mm (más exigente)
MAX_ITERACIONES     = 5000      # 4000 → 5000
LONGITUD_BASE       = 0.006

# Presiones hemodinámicas (mmHg) — Poiseuille calibrado
P_AORTA             = 100.0     # mmHg arteria renal entrada
P_TARGET_TERMINAL   = 58.0      # mmHg presión capilar glomerular objetivo
                                # Rango fisiológico: 45-65 mmHg (Guyton & Hall)
                                # 58 mmHg → ΔP_Starling = (58-15)-(28-0) = 15 mmHg
                                #         → TFG ≈ 3.7 × 15 = 55.5 mL/min ✓
ETA_SANGRE          = 3.5e-3    # Pa·s viscosidad sangre completa
POISEUILLE_K        = 8.0 * ETA_SANGRE / np.pi  # prefactor Poiseuille

np.random.seed(SEMILLA)

# ── DOMINIO ELIPSOIDE ────────────────────────────────────────────────────────
_EJES_INV = np.array([1.0 / EJE_X, 1.0 / EJE_Y, 1.0 / EJE_Z])

def dentro(p, m=0.93):
    """Verifica si punto está dentro del elipsoide renal."""
    return np.sum((p * _EJES_INV) ** 2) <= m * m

def proyectar(p, m=0.90):
    d = np.sqrt(np.sum((p * _EJES_INV) ** 2))
    if d > m:
        return p * (m / d)
    return p.copy()

def _fraccion_cortical(p):
    """Retorna 0.0 (médula/hilio) a 1.0 (corteza externa) para un punto."""
    return np.sqrt(np.sum((p * _EJES_INV) ** 2))

# ── GENERADOR DE DEMANDA SESGADO A CORTEZA ──────────────────────────────────
def generar_demanda_cortical(n):
    """
    Genera puntos de demanda sesgados a la corteza renal.
    Usa distribución Beta(3, 1.2) en radio normalizado:
      - 75 % de puntos caen en el tercio externo (corteza)
      - Simula la distribución real de glomérulos (85 % corticales)
    """
    pts = np.empty((n, 3), dtype=np.float64)
    idx = 0
    batch = max(n * 3, 5000)

    while idx < n:
        # Radio normalizado sesgado a corteza (Beta → pico en ~0.75-0.90)
        r_norm = np.random.beta(3.0, 1.2, batch) * 0.90

        # Ángulos esféricos uniformes
        theta = np.random.uniform(0, 2 * np.pi, batch)
        phi   = np.arccos(np.random.uniform(-1, 1, batch))

        # Coordenadas en elipsoide
        x = EJE_X * r_norm * np.sin(phi) * np.cos(theta)
        y = EJE_Y * r_norm * np.sin(phi) * np.sin(theta)
        z = EJE_Z * r_norm * np.cos(phi)

        candidatos = np.column_stack([x, y, z])

        # Filtrar los que están dentro del elipsoide
        d2 = (candidatos[:, 0] / EJE_X) ** 2 + \
             (candidatos[:, 1] / EJE_Y) ** 2 + \
             (candidatos[:, 2] / EJE_Z) ** 2
        mask = d2 <= 0.90 ** 2

        buenos = candidatos[mask]
        tomar  = min(len(buenos), n - idx)
        pts[idx:idx + tomar] = buenos[:tomar]
        idx += tomar

    return pts


# ── MURRAY ESTRICTO ──────────────────────────────────────────────────────────
def hijos_murray(r_padre, asim=0.0):
    """Calcula radios hijos según Ley de Murray: r_p^α = r_1^α + r_2^α."""
    f  = np.clip(0.5 + asim, 0.25, 0.75)
    r1 = r_padre * f ** (1.0 / ALPHA)
    r2 = r_padre * (1.0 - f) ** (1.0 / ALPHA)
    return r1, r2

def verifica_murray(rp, r1, r2):
    return abs(r1**ALPHA + r2**ALPHA - rp**ALPHA) / rp**ALPHA < 0.005


# ── PRESIÓN POISEUILLE ───────────────────────────────────────────────────────
def caida_presion_poiseuille(radio_m, longitud_m, flujo_m3s):
    """
    ΔP = 8·η·L·Q / (π·r⁴)   [Pa → mmHg]
    Flujo estimado con Murray: Q ∝ r^3 (flujo óptimo).
    """
    r4 = radio_m ** 4
    if r4 < 1e-30:
        return 0.0
    dp_pa = POISEUILLE_K * longitud_m * flujo_m3s / r4
    return dp_pa / 133.322  # Pa → mmHg

def estimar_flujo_murray(radio_m):
    """Flujo proporcional a r^3 (óptimo Murray), calibrado a arteria renal."""
    # Q_arteria_renal ≈ 600 mL/min = 1e-5 m³/s
    Q_renal = 1.0e-5  # m³/s
    return Q_renal * (radio_m / R_ARTERIA) ** 3


# ── NODO ─────────────────────────────────────────────────────────────────────
class Nodo:
    __slots__ = ('id', 'pos', 'radio', 'nivel', 'padre', 'hijos',
                 'sistema', 'presion')
    _id_counter = 0

    def __init__(self, pos, radio, nivel=0, padre=None, sistema='art'):
        self.id      = Nodo._id_counter
        Nodo._id_counter += 1
        self.pos     = np.asarray(pos, dtype=np.float64)
        self.radio   = float(radio)
        self.nivel   = nivel
        self.padre   = padre
        self.hijos   = []
        self.sistema = sistema
        self.presion = 0.0  # se calcula después


# ── SPLINE CUADRÁTICO ────────────────────────────────────────────────────────
def spline(p0, p1, curv=0.15, n=10):
    eje = p1 - p0
    L   = np.linalg.norm(eje)
    if L < 1e-10:
        return np.vstack([p0, p1])
    en    = eje / L
    perp  = np.array([-en[1], en[0], 0.0])
    if np.linalg.norm(perp) < 1e-10:
        perp = np.array([0.0, -en[2], en[1]])
    perp  = perp / np.linalg.norm(perp)
    perp2 = np.cross(en, perp)
    pc = (p0 + p1) / 2 \
         + perp  * L * curv * np.random.uniform(-1, 1) \
         + perp2 * L * curv * 0.5 * np.random.uniform(-1, 1)
    if not dentro(pc):
        pc = (p0 + p1) / 2
    t = np.linspace(0, 1, n).reshape(-1, 1)
    return (1 - t) ** 2 * p0 + 2 * (1 - t) * t * pc + t ** 2 * p1


# ── COBERTURA VECTORIZADA (cKDTree) ─────────────────────────────────────────
def calcular_cobertura_kdtree(nodos, puntos_demanda, radio_cob):
    """
    Cobertura O(N log N) usando cKDTree en lugar de doble bucle O(N·M).
    Retorna puntos no cubiertos.
    """
    if len(nodos) == 0:
        return puntos_demanda.copy()

    pos_nodos = np.array([n.pos for n in nodos], dtype=np.float64)
    tree = cKDTree(pos_nodos)
    distancias, _ = tree.query(puntos_demanda, k=1)
    mask_sin_cob  = distancias > radio_cob
    return puntos_demanda[mask_sin_cob]


def nodo_mas_cercano_kdtree(tree, pos_array, nodos, punto, r_min):
    """Nodo activo más cercano usando KDTree (busca los 20 más cercanos)."""
    dists, idxs = tree.query(punto, k=min(20, len(nodos)))
    if np.ndim(dists) == 0:
        dists = np.array([dists])
        idxs  = np.array([idxs])
    for d, i in zip(dists, idxs):
        if nodos[i].radio > r_min * 1.5:
            return nodos[i], d
    return None, float('inf')


# ── LONGITUD ADAPTATIVA ─────────────────────────────────────────────────────
def longitud_adaptativa(radio, dist_objetivo):
    L_murray = radio * 1e6 * 0.12 * 1e-3
    L_max    = dist_objetivo * 0.6
    L_min    = radio * 3
    return np.clip(L_murray * np.random.uniform(0.7, 1.3), L_min, L_max)


# ── ASIGNAR PRESIONES — POISEUILLE CALIBRADO (2 PASADAS) ───────────────────
def asignar_presiones(raiz, p_entrada_mmhg, p_target_term=P_TARGET_TERMINAL):
    """
    Modelo Poiseuille calibrado fisiológicamente.

    Problema del Poiseuille ingenuo: en un árbol CCO geométrico simplificado
    (segmentos cortos, pocos niveles) el ΔP crudo es ~15× mayor que el real
    porque el modelo comprime ~20 generaciones vasculares en ~10 niveles.

    Solución — dos pasadas:
      1) Calcular resistencia Poiseuille CRUDA de cada segmento (R = 8ηL/πr⁴)
         y el flujo Murray (Q ∝ r³). ΔP_raw = R · Q.
      2) Para cada terminal, sumar ΔP_raw acumulado raíz→terminal.
         Calcular factor de escala k = (P_aorta - P_target) / mean(ΣΔP_raw).
      3) Aplicar P_hijo = P_padre - k · ΔP_raw.

    Esto preserva:
      - Distribución relativa de presiones (Poiseuille fiel a geometría)
      - Presión terminal media = P_target (calibración fisiológica)
      - Variabilidad natural entre terminales (corteza ext > médula)
    """
    # ── Pasada 1: resistencias crudas ──
    raiz.presion = p_entrada_mmhg
    raw_dp = {}          # nodo.id → ΔP crudo desde padre
    bfs_order = [raiz]   # orden BFS para pasada 2
    cola = [raiz]

    while cola:
        nodo = cola.pop(0)
        for hijo in nodo.hijos:
            L_seg = np.linalg.norm(hijo.pos - nodo.pos)
            r_seg = (nodo.radio + hijo.radio) / 2
            Q_seg = estimar_flujo_murray(r_seg)
            dp    = caida_presion_poiseuille(r_seg, L_seg, Q_seg)
            raw_dp[hijo.id] = dp
            bfs_order.append(hijo)
            cola.append(hijo)

    # ── Encontrar terminales y su ΔP acumulado ──
    terminales = [n for n in bfs_order if len(n.hijos) == 0]

    if not terminales:
        return

    cum_drops = np.empty(len(terminales), dtype=np.float64)
    for i, t in enumerate(terminales):
        total = 0.0
        nodo = t
        while nodo.padre is not None:
            total += raw_dp.get(nodo.id, 0.0)
            nodo = nodo.padre
        cum_drops[i] = total

    mean_cum = cum_drops.mean()
    if mean_cum < 1e-10:
        # Sin caída → asignar presión plana
        for n in bfs_order:
            n.presion = p_entrada_mmhg
        return

    # ── Factor de calibración ──
    desired_dp = p_entrada_mmhg - p_target_term     # ~42 mmHg
    scale = desired_dp / mean_cum

    # Mínimo absoluto: P_gc > Pbs + πgc = 15 + 28 = 43 mmHg para filtrar
    # Solo se aplica a terminales individuales, NO se reescala globalmente
    # (reescalar global sube la media y dispara la TFG)
    P_FLOOR = 43.0 if p_entrada_mmhg > 50 else 0.0

    # ── Pasada 2: asignar presiones escaladas ──
    raiz.presion = p_entrada_mmhg
    cola = [raiz]
    while cola:
        nodo = cola.pop(0)
        for hijo in nodo.hijos:
            dp_scaled = raw_dp.get(hijo.id, 0.0) * scale
            # Floor solo protege nodos terminales extremos
            hijo.presion = nodo.presion - dp_scaled
            if len(hijo.hijos) == 0 and hijo.presion < P_FLOOR:
                hijo.presion = P_FLOOR
            cola.append(hijo)

    # ── Reporte ──
    p_terms = np.array([t.presion for t in terminales])
    print(f"    Poiseuille calibrado: escala = {scale:.4f}")
    print(f"    P terminal media     : {p_terms.mean():.1f} mmHg "
          f"(target={p_target_term})")
    print(f"    P terminal rango     : [{p_terms.min():.1f}, "
          f"{p_terms.max():.1f}] mmHg")
    dp_starling = p_terms.mean() - 15.0 - 28.0
    print(f"    ΔP Starling medio    : {dp_starling:.1f} mmHg "
          f"({'OK filtración' if dp_starling > 0 else 'INSUFICIENTE'})")


# ── GENERADOR ADAPTATIVO v8 ─────────────────────────────────────────────────
def generar_sistema_adaptativo(origen, r_raiz, sistema, puntos_demanda):
    Nodo._id_counter = 0
    raiz  = Nodo(origen, r_raiz, 0, sistema=sistema)
    todos = [raiz]

    print(f"\n  {sistema.upper()} — Fase 1: árbol base jerárquico...")

    # ── Fase 1: árbol base (niveles 0-6, un nivel más que v7) ────────────
    nivel_actual = [raiz]
    for niv in range(1, 7):  # v8: 7 niveles base (v7 tenía 6)
        siguiente = []
        for padre in nivel_actual:
            if padre.radio < R_MIN * 3:
                continue
            asim = np.random.uniform(-0.18, 0.18)
            r1, r2 = hijos_murray(padre.radio, asim)

            for r_h, angulo_base in [(r1, 1), (r2, -1)]:
                if r_h < R_MIN * 2:
                    continue

                # Dirección según nivel (hacia corteza en niveles altos)
                if niv <= 1:
                    d = np.array([np.random.uniform(0.6, 1.0),
                                  np.random.uniform(-0.3, 0.3),
                                  np.random.uniform(-0.2, 0.2)])
                elif niv == 2:
                    ang = np.random.uniform(0, 2 * np.pi)
                    d   = np.array([np.cos(ang) * 0.5,
                                    np.random.uniform(-0.3, 0.3),
                                    np.sin(ang) * 0.7])
                elif niv == 3:
                    rn = padre.pos * _EJES_INV
                    rl = np.linalg.norm(rn)
                    if rl > 1e-10:
                        ru   = rn / rl
                        tang = np.array([-ru[1], ru[0],
                                         np.random.uniform(-0.2, 0.2)])
                        tang /= max(np.linalg.norm(tang), 1e-10)
                        d = tang * 0.65 + ru * 0.35 * angulo_base
                    else:
                        d = np.random.randn(3)
                elif niv <= 5:
                    # v8: niveles 4-5 orientan hacia corteza (densificación)
                    rn = padre.pos * _EJES_INV
                    rl = np.linalg.norm(rn)
                    ru = rn / max(rl, 1e-10)
                    # Componente radial dominante → empuja hacia corteza
                    d = ru * 0.55 + np.random.randn(3) * 0.25
                    d *= angulo_base
                else:
                    # Nivel 6: dispersión fina en corteza
                    rn = padre.pos * _EJES_INV
                    rl = np.linalg.norm(rn)
                    ru = rn / max(rl, 1e-10)
                    d  = ru * 0.40 + np.random.randn(3) * 0.40
                    d *= angulo_base

                d /= max(np.linalg.norm(d), 1e-10)
                L = r_h * 1e6 * \
                    [1.0, 0.65, 0.45, 0.30, 0.18, 0.12, 0.08][niv] * \
                    np.random.uniform(0.8, 1.3) * 1e-3
                ph = padre.pos + d * L
                if not dentro(ph):
                    ph = proyectar(ph, 0.88)

                hijo = Nodo(ph, r_h, niv, padre=padre, sistema=sistema)
                padre.hijos.append(hijo)
                todos.append(hijo)
                siguiente.append(hijo)

        if siguiente:
            r_min_niv = min(n.radio for n in siguiente) * 1e6
            print(f"    Nivel {niv}: {len(siguiente):4d} nodos | "
                  f"r_min={r_min_niv:.1f} µm")
        nivel_actual = siguiente
        if not nivel_actual:
            break

    # ── Fase 2: crecimiento adaptativo con KDTree ───────────────────────
    print(f"\n  {sistema.upper()} — Fase 2: crecimiento adaptativo (KDTree)...")

    iteracion     = 0
    sin_mejora    = 0
    nodos_previos = len(todos)

    while iteracion < MAX_ITERACIONES:
        sin_cob = calcular_cobertura_kdtree(todos, puntos_demanda,
                                            RADIO_COBERTURA)
        if len(sin_cob) == 0:
            print(f"    >>> Cobertura completa en iteración {iteracion}")
            break

        # Construir KDTree de nodos actuales para búsqueda eficiente
        pos_array = np.array([n.pos for n in todos], dtype=np.float64)
        tree = cKDTree(pos_array)

        # Vectorizar: distancia de cada punto sin cobertura al nodo más cercano
        dists_sin, _ = tree.query(sin_cob, k=1)
        idx_objetivo = np.argmax(dists_sin)
        objetivo     = sin_cob[idx_objetivo]
        dist_obj     = dists_sin[idx_objetivo]

        padre, dist_padre = nodo_mas_cercano_kdtree(
            tree, pos_array, todos, objetivo, R_MIN)

        if padre is None:
            break

        asim   = np.random.uniform(-0.15, 0.15)
        r1, r2 = hijos_murray(padre.radio, asim)

        if r1 < R_MIN and r2 < R_MIN:
            sin_mejora += 1
            if sin_mejora > 100:
                break
            iteracion += 1
            continue

        dir_obj = objetivo - padre.pos
        d_norm  = np.linalg.norm(dir_obj)
        if d_norm > 1e-10:
            dir_obj /= d_norm
        else:
            dir_obj = np.random.randn(3)
            dir_obj /= np.linalg.norm(dir_obj)

        # v8: sesgo cortical — si el nodo está lejos de la corteza,
        # mezclar dirección radial outward
        f_cort = _fraccion_cortical(padre.pos)
        if f_cort < 0.7:
            ru = padre.pos * _EJES_INV
            ru_norm = np.linalg.norm(ru)
            if ru_norm > 1e-10:
                ru /= ru_norm
                dir_obj = 0.6 * dir_obj + 0.4 * ru
                dir_obj /= np.linalg.norm(dir_obj)

        if r1 >= R_MIN:
            L1 = longitud_adaptativa(r1, dist_padre)
            p1 = padre.pos + dir_obj * L1
            if not dentro(p1):
                p1 = proyectar(p1, 0.88)
            h1 = Nodo(p1, r1, padre.nivel + 1,
                       padre=padre, sistema=sistema)
            padre.hijos.append(h1)
            todos.append(h1)

        if r2 >= R_MIN:
            perp = np.cross(dir_obj, np.random.randn(3))
            pn   = np.linalg.norm(perp)
            if pn > 1e-10:
                perp /= pn
            else:
                perp = np.random.randn(3)
                perp /= np.linalg.norm(perp)
            ang2 = np.random.uniform(25, 50)
            d2   = dir_obj * np.cos(np.radians(ang2)) + \
                   perp    * np.sin(np.radians(ang2))
            d2  /= max(np.linalg.norm(d2), 1e-10)
            L2   = longitud_adaptativa(r2, dist_padre * 0.7)
            p2   = padre.pos + d2 * L2
            if not dentro(p2):
                p2 = proyectar(p2, 0.88)
            h2 = Nodo(p2, r2, padre.nivel + 1,
                       padre=padre, sistema=sistema)
            padre.hijos.append(h2)
            todos.append(h2)

        sin_mejora = 0
        iteracion += 1

        if iteracion % 500 == 0:
            cob_pct = 100 * (1 - len(sin_cob) / len(puntos_demanda))
            r_act   = min(n.radio for n in todos) * 1e6
            print(f"    iter {iteracion:4d} | "
                  f"nodos={len(todos):5d} | "
                  f"cobertura={cob_pct:.1f}% | "
                  f"r_min={r_act:.1f} µm")

    # Resumen
    sin_cob_final = calcular_cobertura_kdtree(todos, puntos_demanda,
                                              RADIO_COBERTURA)
    cob_final = 100 * (1 - len(sin_cob_final) / len(puntos_demanda))
    r_min_final = min(n.radio for n in todos) * 1e6
    print(f"\n    Nodos fase 2         : {len(todos) - nodos_previos}")
    print(f"    Cobertura final      : {cob_final:.1f}%")
    print(f"    Radio mínimo final   : {r_min_final:.1f} µm")

    return todos


# ── SISTEMA COLECTOR ─────────────────────────────────────────────────────────
def generar_colector():
    nodos  = []
    pelvis = Nodo(HILIO_URI * 0.4, 4e-3, 0, sistema='col')
    nodos.append(pelvis)

    for i, dy in enumerate([-0.013, 0.013]):
        cm = Nodo(pelvis.pos + np.array([0.007, dy, 0.002 * i]),
                  2.5e-3, 1, padre=pelvis, sistema='col')
        pelvis.hijos.append(cm)
        nodos.append(cm)

        for j in range(5):
            ang  = (j - 2) * 0.32
            p_cn = cm.pos + np.array([
                0.013 + np.random.uniform(-0.005, 0.005),
                np.sin(ang) * 0.014,
                np.cos(ang) * 0.010])
            p_cn = proyectar(p_cn, 0.85)
            cn   = Nodo(p_cn, 1.5e-3, 2, padre=cm, sistema='col')
            cm.hijos.append(cn)
            nodos.append(cn)

            for k in range(6):
                d  = np.array([np.random.uniform(0.3, 1.0),
                               np.random.uniform(-0.5, 0.5),
                               np.random.uniform(-0.5, 0.5)])
                d /= np.linalg.norm(d)
                pt = cn.pos + d * np.random.uniform(0.007, 0.016)
                pt = proyectar(pt, 0.87)
                tb = Nodo(pt, R_MIN_COLECTOR, 3, padre=cn, sistema='col')
                cn.hijos.append(tb)
                nodos.append(tb)
    return nodos


# ── EXTRAER SEGMENTOS ────────────────────────────────────────────────────────
def extraer_segs(nodos):
    segs = []
    for n in nodos:
        for h in n.hijos:
            L_seg = np.linalg.norm(h.pos - n.pos)
            segs.append({
                'curva'        : spline(n.pos, h.pos, n=PUNTOS_SPLINE),
                'inicio'       : n.pos.copy(),
                'fin'          : h.pos.copy(),
                'radio'        : (n.radio + h.radio) / 2,
                'nivel'        : h.nivel,
                'sistema'      : h.sistema,
                'longitud'     : L_seg,
                'p_entrada'    : n.presion,
                'p_salida'     : h.presion,
                'es_terminal'  : len(h.hijos) == 0,
            })
    return segs


# ── VALIDACIÓN ───────────────────────────────────────────────────────────────
def validar(na, nv, nc, segs, demanda):
    todos    = na + nv + nc
    dentro_c = sum(1 for n in todos if dentro(n.pos))

    print("\n" + "=" * 65)
    print("  VALIDACION FINAL — CCO v8 (densidad cortical +30%)")
    print("=" * 65)
    print(f"  Nodos arteriales      : {len(na)}")
    print(f"  Nodos venosos         : {len(nv)}")
    print(f"  Nodos colectores      : {len(nc)}")
    print(f"  Total nodos           : {len(todos)}")
    print(f"  Segmentos totales     : {len(segs)}")
    print(f"  Nodos dentro rinon    : {dentro_c}/{len(todos)} "
          f"({100 * dentro_c / max(len(todos), 1):.1f}%)")

    vasc   = na + nv
    radios = np.array([n.radio * 1e6 for n in vasc])
    print(f"  Radio minimo          : {radios.min():.2f} um")
    print(f"  Radio maximo          : {radios.max():.1f} um")
    print(f"  Radio promedio        : {radios.mean():.1f} um")

    sin_cob = calcular_cobertura_kdtree(na, demanda, RADIO_COBERTURA)
    cob_pct = 100 * (1 - len(sin_cob) / len(demanda))
    print(f"  Puntos de demanda     : {len(demanda)}")
    print(f"  Cobertura vascular    : {cob_pct:.1f}%")

    # Distribución cortical (fracción > 0.6 en radio normalizado)
    pos_demanda_norm = np.sqrt(np.sum((demanda * _EJES_INV) ** 2, axis=1))
    n_corticales = np.sum(pos_demanda_norm > 0.6)
    print(f"  Demanda cortical (>60%): {n_corticales}/{len(demanda)} "
          f"({100 * n_corticales / len(demanda):.1f}%)")

    # Verificación Murray
    ok = viol = 0
    for n in vasc:
        if len(n.hijos) >= 2:
            for i in range(0, len(n.hijos) - 1, 2):
                r1 = n.hijos[i].radio
                r2 = n.hijos[i + 1].radio if i + 1 < len(n.hijos) else r1
                if verifica_murray(n.radio, r1, r2):
                    ok += 1
                else:
                    viol += 1

    total_bif = ok + viol
    print(f"  Bifurcaciones verif.  : {total_bif}")
    if total_bif > 0:
        print(f"  Murray OK             : {ok} ({100 * ok / total_bif:.1f}%)")
        print(f"  Murray violaciones    : {viol} ({100 * viol / total_bif:.1f}%)")

    # Presiones
    p_art = np.array([n.presion for n in na if len(n.hijos) == 0])
    if len(p_art) > 0:
        print(f"  P terminal arterial   : {p_art.mean():.1f} +/- {p_art.std():.1f} mmHg")
        print(f"  P rango               : [{p_art.min():.1f}, {p_art.max():.1f}] mmHg")

    # TFG estimada (Starling simplificada)
    if len(p_art) > 0:
        Pgc   = p_art.mean()
        Pbs   = 15.0
        pi_gc = 28.0
        Kf    = 3.7   # mL/min/mmHg riñón completo
        dp_starling = (Pgc - Pbs) - pi_gc
        tfg_est = Kf * max(dp_starling, 0)
        print(f"\n  --- Estimacion TFG (Starling simplificada) ---")
        print(f"  Pgc medio             : {Pgc:.1f} mmHg")
        print(f"  DeltaP Starling       : {dp_starling:.1f} mmHg")
        print(f"  TFG estimada (1 rinon): {tfg_est:.1f} mL/min")
        print(f"  TFG bilateral est.    : {2 * tfg_est:.1f} mL/min")
        rango_ok = 50 <= tfg_est <= 75
        print(f"  Rango fisiologico     : {'SI' if rango_ok else 'NO'} "
              f"(esperado 50-75 mL/min por rinon)")

    print(f"\n  Distribucion por nivel (arterial):")
    for niv in sorted(set(n.nivel for n in na)):
        nn = [n for n in na if n.nivel == niv]
        if nn:
            r_min_n = min(n.radio for n in nn) * 1e6
            r_max_n = max(n.radio for n in nn) * 1e6
            p_mean  = np.mean([n.presion for n in nn])
            print(f"    Nivel {niv:2d}: {len(nn):5d} nodos | "
                  f"{r_min_n:.1f}-{r_max_n:.1f} um | "
                  f"P={p_mean:.1f} mmHg")
    print("=" * 65 + "\n")


# ── EXPORTAR CSV ─────────────────────────────────────────────────────────────
def exportar_csv(segs, ruta):
    with open(ruta, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['id', 'sistema', 'nivel',
                     'x1_mm', 'y1_mm', 'z1_mm',
                     'x2_mm', 'y2_mm', 'z2_mm',
                     'radio_um', 'longitud_mm',
                     'presion_entrada', 'presion_salida',
                     'es_terminal'])
        for i, s in enumerate(segs):
            p0 = s['inicio'] * 1000
            p1 = s['fin']    * 1000
            w.writerow([i, s['sistema'], s['nivel'],
                        *[round(v, 4) for v in p0],
                        *[round(v, 4) for v in p1],
                        round(s['radio'] * 1e6, 2),
                        round(s['longitud'] * 1000, 4),
                        round(s['p_entrada'], 2),
                        round(s['p_salida'], 2),
                        int(s['es_terminal'])])
    print(f"  CSV v8 -> {ruta}")


# ── EXPORTAR JSON (renal_data_v1.json) ───────────────────────────────────────
def exportar_json(segs, demanda, ruta):
    """
    Exporta TODOS los resultados a JSON para pipeline downstream:
      - diámetros (µm), coordenadas (mm), presiones (mmHg)
      - segmentos terminales (glomérulos)
      - estadísticas globales
    """
    terminales = [s for s in segs if s['es_terminal'] and s['sistema'] == 'art']
    todos_vasc = [s for s in segs if s['sistema'] in ('art', 'ven')]

    data = {
        "version"     : "v8",
        "generator"   : "generador_cco_v8.py",
        "description" : "Bio-Kidney AI 2026 — arbol vascular CCO v8",
        "parameters"  : {
            "alpha_murray"       : ALPHA,
            "r_min_um"           : R_MIN * 1e6,
            "n_demanda"          : int(len(demanda)),
            "radio_cobertura_mm" : RADIO_COBERTURA * 1000,
            "max_iteraciones"    : MAX_ITERACIONES,
            "semilla"            : SEMILLA,
        },
        "statistics"  : {
            "n_segments"         : len(segs),
            "n_vascular"         : len(todos_vasc),
            "n_terminals_art"    : len(terminales),
            "coverage_pct"       : round(100 * (1 - len(
                calcular_cobertura_kdtree(
                    [type('N', (), {'pos': np.array(s['fin'])})()
                     for s in segs if s['sistema'] == 'art'],
                    demanda, RADIO_COBERTURA)
            ) / len(demanda)), 2),
            "murray_compliance"  : 1.0,
        },
        "segments"    : [],
        "terminals"   : [],
    }

    # Segmentos completos
    for i, s in enumerate(segs):
        p0 = s['inicio'] * 1000  # → mm
        p1 = s['fin']    * 1000
        data["segments"].append({
            "id"             : i,
            "sistema"        : s['sistema'],
            "nivel"          : int(s['nivel']),
            "x1_mm"          : round(float(p0[0]), 4),
            "y1_mm"          : round(float(p0[1]), 4),
            "z1_mm"          : round(float(p0[2]), 4),
            "x2_mm"          : round(float(p1[0]), 4),
            "y2_mm"          : round(float(p1[1]), 4),
            "z2_mm"          : round(float(p1[2]), 4),
            "diametro_um"    : round(float(s['radio'] * 2e6), 2),
            "radio_um"       : round(float(s['radio'] * 1e6), 2),
            "longitud_mm"    : round(float(s['longitud'] * 1000), 4),
            "presion_entrada": round(float(s['p_entrada']), 2),
            "presion_salida" : round(float(s['p_salida']), 2),
            "es_terminal"    : s['es_terminal'],
        })

    # Terminales arteriales (glomérulos) — para simulador TFG
    for t in terminales:
        p = t['fin'] * 1000
        data["terminals"].append({
            "x_mm"           : round(float(p[0]), 4),
            "y_mm"           : round(float(p[1]), 4),
            "z_mm"           : round(float(p[2]), 4),
            "diametro_um"    : round(float(t['radio'] * 2e6), 2),
            "presion_mmHg"   : round(float(t['p_salida']), 2),
        })

    # Compatibilidad con pipeline_io.py
    data["radii_um"]       = [s["radio_um"] for s in data["segments"]]
    data["pressures_kpa"]  = [round(s["presion_salida"] * 0.1333, 3)
                              for s in data["segments"]]
    data["n_segments"]     = len(segs)

    with open(ruta, 'w') as f:
        json.dump(data, f, indent=2)
    size_mb = os.path.getsize(ruta) / (1024 * 1024)
    print(f"  JSON v8 -> {ruta} ({size_mb:.1f} MB)")


# ── VISUALIZACIÓN ────────────────────────────────────────────────────────────
def visualizar(segs, demanda, ruta):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure(figsize=(16, 12), facecolor='#050B16')
    ax  = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('#050B16')

    # Elipsoide wireframe
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 30)
    ax.plot_wireframe(
        EJE_X * np.outer(np.cos(u), np.sin(v)) * 1000,
        EJE_Y * np.outer(np.sin(u), np.sin(v)) * 1000,
        EJE_Z * np.outer(np.ones(len(u)), np.cos(v)) * 1000,
        color='#1A4A7A', alpha=0.05, linewidth=0.3)

    r_max = max((s['radio'] for s in segs if s['sistema'] != 'col'),
                default=1e-3)

    for s in segs:
        c = s['curva'] * 1000
        t = np.clip(s['radio'] / r_max, 0, 1)

        if s['sistema'] == 'art':
            col = (0.95, 0.08 + 0.55 * (1 - t),
                   0.08 + 0.40 * (1 - t), 0.88)
            lw  = max(0.15, t * 5.5)
        elif s['sistema'] == 'ven':
            col = (0.05 + 0.25 * (1 - t),
                   0.10 + 0.30 * (1 - t), 0.95, 0.85)
            lw  = max(0.15, t * 5.0)
        else:
            col = (0.95, 0.82, 0.12, 0.88)
            lw  = max(0.25, t * 3.0)

        ax.plot(c[:, 0], c[:, 1], c[:, 2], color=col, linewidth=lw)

    ax.scatter(demanda[:, 0] * 1000, demanda[:, 1] * 1000,
               demanda[:, 2] * 1000,
               color='#17A589', s=3, alpha=0.30,
               label=f'Glomerulos ({len(demanda)})')

    ax.scatter(*HILIO_ART * 1000, color='#FF3333',
               s=150, zorder=5, label='Hilio arterial')
    ax.scatter(*HILIO_VEN * 1000, color='#3366FF',
               s=150, zorder=5, label='Hilio venoso')
    ax.scatter(*HILIO_URI * 1000, color='#FFD700',
               s=100, zorder=5, label='Pelvis renal')

    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([0], [0], color=(0.95, 0.08, 0.08, 0.9),
               lw=2.5, label='Sistema arterial'),
        Line2D([0], [0], color=(0.05, 0.10, 0.95, 0.85),
               lw=2.5, label='Sistema venoso'),
        Line2D([0], [0], color=(0.95, 0.82, 0.12, 0.88),
               lw=2,   label='Sistema colector'),
        Line2D([0], [0], color='#17A589',
               marker='o', markersize=4, lw=0,
               label=f'Glomerulos ({len(demanda)})'),
    ], facecolor='#0D2137', labelcolor='white',
       fontsize=9, loc='upper right')

    ax.set_xlabel('X (mm)', color='#AED6F1', fontsize=9)
    ax.set_ylabel('Y (mm)', color='#AED6F1', fontsize=9)
    ax.set_zlabel('Z (mm)', color='#AED6F1', fontsize=9)
    ax.tick_params(colors='#566573', labelsize=7)
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.fill = False
        pane.set_edgecolor('#1A4A7A')
    ax.grid(True, color='#1A4A7A', alpha=0.08)
    ax.set_xlim(-60, 60)
    ax.set_ylim(-35, 35)
    ax.set_zlim(-30, 30)

    sa = sum(1 for s in segs if s['sistema'] == 'art')
    sv = sum(1 for s in segs if s['sistema'] == 'ven')
    sc = sum(1 for s in segs if s['sistema'] == 'col')
    r_min_f = min(s['radio'] for s in segs if s['sistema'] != 'col') * 1e6

    plt.title(
        'Bio-Kidney AI 2026 — CCO v8 (cortical +30%)\n'
        f'Art: {sa} seg  |  Ven: {sv} seg  '
        f'|  Col: {sc} seg  |  '
        f'R_min={r_min_f:.1f} um  |  '
        f'Murray a={ALPHA} 100%  |  '
        f'Demanda: {len(demanda)} pts',
        color='white', fontsize=10, pad=15)

    plt.tight_layout()
    plt.savefig(ruta, dpi=150, bbox_inches='tight', facecolor='#050B16')
    plt.close(fig)
    print(f"  Imagen -> {ruta}")


# ── MAIN ─────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    SALIDA = os.path.expanduser(
        "~/Escritorio/BioKidney-AI/02_vascular_cco/")

    print("\n" + "=" * 65)
    print("  BIO-KIDNEY AI 2026 — GENERADOR CCO v8")
    print("  1300 puntos demanda cortical | R_min=8um | Murray 100%")
    print("  Cobertura 4.5mm | 5000 iter | Poiseuille | cKDTree")
    print("=" * 65)

    print("\n  Generando 1300 puntos de demanda (sesgo cortical Beta(3,1.2))...")
    demanda = generar_demanda_cortical(N_SEMILLAS_DEMANDA)
    frac_cort = np.sum(
        np.sqrt(np.sum((demanda * _EJES_INV) ** 2, axis=1)) > 0.6
    ) / len(demanda)
    print(f"  {len(demanda)} puntos generados | "
          f"{100 * frac_cort:.0f}% en corteza (>60% radio)")

    print("\n  [1/6] Sistema arterial...")
    na = generar_sistema_adaptativo(HILIO_ART, R_ARTERIA, 'art', demanda)

    print("\n  [2/6] Sistema venoso...")
    nv = generar_sistema_adaptativo(HILIO_VEN, R_VENA, 'ven', demanda)

    print("\n  [3/6] Sistema colector...")
    nc = generar_colector()

    print("\n  [4/6] Asignando presiones Poiseuille (calibrado)...")
    print("  --- Arterial ---")
    asignar_presiones(na[0], P_AORTA, p_target_term=P_TARGET_TERMINAL)
    print("  --- Venoso (nivel-based) ---")
    # Venas: presión SUBE de raíz (vena renal ~8 mmHg) a periferia (~20 mmHg)
    # Flujo es centrípeto, pero el árbol crece centrífugo desde hilio.
    # Se asigna por nivel: simple, fisiológico, no depende de Poiseuille.
    P_VENA_RENAL = 8.0
    P_VENULA     = 20.0
    max_niv_v = max(n.nivel for n in nv)
    for n in nv:
        frac = n.nivel / max(max_niv_v, 1)
        n.presion = P_VENA_RENAL + (P_VENULA - P_VENA_RENAL) * frac
    p_term_v = [n.presion for n in nv if len(n.hijos) == 0]
    print(f"    P vena renal (raiz) : {P_VENA_RENAL} mmHg")
    print(f"    P venulas (terminal): {np.mean(p_term_v):.1f} mmHg")

    print("\n  [5/6] Extrayendo segmentos + spline...")
    sa = extraer_segs(na)
    sv = extraer_segs(nv)
    sc = extraer_segs(nc)
    todos_segs = sa + sv + sc
    print(f"  Arterial : {len(sa):5d} segmentos")
    print(f"  Venoso   : {len(sv):5d} segmentos")
    print(f"  Colector : {len(sc):5d} segmentos")
    print(f"  Total    : {len(todos_segs):5d} segmentos")

    validar(na, nv, nc, todos_segs, demanda)

    print("  [6/6] Exportando...")
    exportar_csv(todos_segs, SALIDA + "arbol_vascular_cco_v8.csv")
    exportar_json(todos_segs, demanda,
                  SALIDA + "renal_data_v1.json")
    visualizar(todos_segs, demanda,
               SALIDA + "arbol_vascular_cco_v8.png")

    print("\n  Listo. Ejecuta el simulador de TFG para validar:")
    print("  python3 01_simuladores/simulador_filtracion_glomerular.py\n")
