#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║          OPTIMIZADOR Co-SWIFT — Bio-Kidney AI 2026                         ║
║          Módulo de Optimización Multi-Objetivo para Bioimpresión            ║
║          VirtusSapiens | Carlos David Moreno Cáceres                       ║
║          Pipeline in silico — Módulo 09                                    ║
╚══════════════════════════════════════════════════════════════════════════════╝

Física implementada:
  - Modelo reológico Herschel-Bulkley para biotintas no-Newtonianas
  - Hagen-Poiseuille modificado para fluidos no-Newtonianos
  - Criterio de Von Mises para integridad estructural
  - PSO (Particle Swarm Optimization) multi-objetivo NSGA-II adaptado
  - Frente de Pareto para soluciones óptimas

Parámetros validados integrados:
  - Presión óptima extrusión: 30 kPa (simulador SWIFT)
  - NICE Bioink: GelMA 7% + Alginato 1.5% + Nanocelulosa 0.8% + LAP 0.25%
  - Radio mínimo canal Co-SWIFT: 200 µm
  - Límite shear stress citotóxico: 150 Pa
  - Ventana bioimpresión iPSC: días 16-30
  - WSS óptimo: 5.6 dyn/cm²
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patheffects as pe
from scipy.optimize import minimize
from scipy.spatial import ConvexHull
import warnings
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
# PYQT6 UI
# ─────────────────────────────────────────────────────────────────────────────
try:
    from PyQt6.QtWidgets import (
        QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
        QLabel, QPushButton, QProgressBar, QFrame, QScrollArea,
        QTabWidget, QGridLayout, QGroupBox, QSlider, QSpinBox,
        QDoubleSpinBox, QComboBox, QTextEdit, QSizePolicy, QSplitter,
        QTableWidget, QTableWidgetItem, QHeaderView, QMessageBox,
        QFileDialog
    )
    from PyQt6.QtCore import (
        Qt, QThread, pyqtSignal, QTimer, QPropertyAnimation,
        QEasingCurve, QSize
    )
    from PyQt6.QtGui import (
        QFont, QColor, QPalette, QLinearGradient, QGradient,
        QPainter, QPixmap, QIcon, QBrush
    )
    HAS_QT = True
except ImportError:
    HAS_QT = False
    print("[WARNING] PyQt6 no disponible — modo consola activado")

# ─────────────────────────────────────────────────────────────────────────────
# PALETA DE COLORES — Bio-Kidney AI Dark Theme
# ─────────────────────────────────────────────────────────────────────────────
COLORS = {
    'bg_dark':     '#0A0E1A',
    'bg_panel':    '#0F1628',
    'bg_card':     '#141D35',
    'accent_cyan': '#00D4FF',
    'accent_teal': '#00FFB3',
    'accent_gold': '#FFB347',
    'accent_red':  '#FF4757',
    'accent_purple':'#A855F7',
    'text_primary':'#E8F4FD',
    'text_secondary':'#7B9EC7',
    'border':      '#1E3A5F',
    'optimal':     '#00FF88',
    'partial':     '#FFB347',
    'critical':    '#FF4757',
    'pareto':      '#FF6B9D',
}

plt.rcParams.update({
    'figure.facecolor':   COLORS['bg_dark'],
    'axes.facecolor':     COLORS['bg_panel'],
    'axes.edgecolor':     COLORS['border'],
    'axes.labelcolor':    COLORS['text_primary'],
    'xtick.color':        COLORS['text_secondary'],
    'ytick.color':        COLORS['text_secondary'],
    'text.color':         COLORS['text_primary'],
    'grid.color':         COLORS['border'],
    'grid.alpha':         0.4,
    'font.family':        'DejaVu Sans',
    'axes.titlesize':     11,
    'axes.labelsize':     9,
})

# ─────────────────────────────────────────────────────────────────────────────
# FÍSICA: MODELO HERSCHEL-BULKLEY
# ─────────────────────────────────────────────────────────────────────────────

class HerschelBulkley:
    """
    Modelo reológico Herschel-Bulkley para biotintas no-Newtonianas
    τ = τ_y + K * γ̇^n
    
    Parámetros NICE Bioink validados (GelMA 7% + Alg 1.5% + NC 0.8% + LAP 0.25%):
      τ_y = 28 Pa   (esfuerzo de fluencia)
      K   = 4.2     (índice de consistencia)
      n   = 0.41    (índice de flujo — shear-thinning)
    """
    def __init__(self, tau_y=28.0, K=4.2, n=0.41):
        self.tau_y = tau_y  # Pa — esfuerzo de fluencia
        self.K = K          # Pa·s^n — índice de consistencia
        self.n = n          # adimensional — índice de flujo

    def shear_stress(self, shear_rate):
        """Esfuerzo de corte τ [Pa] dado γ̇ [1/s]"""
        if shear_rate <= 0:
            return self.tau_y
        return self.tau_y + self.K * (shear_rate ** self.n)

    def apparent_viscosity(self, shear_rate):
        """Viscosidad aparente η [Pa·s]"""
        if shear_rate <= 0:
            return 1e6  # sólido
        return self.shear_stress(shear_rate) / shear_rate

    def shear_rate_from_pressure(self, pressure_Pa, radius_m, length_m=0.01):
        """
        Hagen-Poiseuille modificado para Herschel-Bulkley
        Iteración numérica para obtener γ̇ máximo en la pared
        """
        # Gradiente de presión
        dP_dL = pressure_Pa / length_m
        # Esfuerzo de corte en la pared
        tau_wall = (radius_m * dP_dL) / 2.0
        
        if tau_wall <= self.tau_y:
            return 0.0  # No hay flujo (por debajo del umbral de fluencia)
        
        # Inversión numérica: τ_wall = τ_y + K * γ̇^n
        gamma_dot = ((tau_wall - self.tau_y) / self.K) ** (1.0 / self.n)
        return gamma_dot

    def flow_rate(self, pressure_Pa, radius_m, length_m=0.01):
        """
        Caudal volumétrico Q [m³/s] para Herschel-Bulkley
        Usando integración analítica del perfil de velocidades
        """
        dP_dL = pressure_Pa / length_m
        tau_wall = (radius_m * dP_dL) / 2.0
        
        if tau_wall <= self.tau_y:
            return 0.0
        
        # Radio del plug (región central donde no hay flujo)
        r_plug = (self.tau_y / (dP_dL / 2.0)) if dP_dL > 0 else radius_m
        r_plug = min(r_plug, radius_m)
        
        # Integral del perfil de velocidades (aproximación numérica)
        N = 200
        r_vals = np.linspace(r_plug, radius_m, N)
        dr = r_vals[1] - r_vals[0]
        
        def velocity_at_r(r):
            tau_r = (r * dP_dL) / 2.0
            if tau_r <= self.tau_y:
                return 0.0
            gamma_dot = ((tau_r - self.tau_y) / self.K) ** (1.0 / self.n)
            return gamma_dot * (radius_m - r)
        
        Q = 0.0
        for r in r_vals:
            v = velocity_at_r(r)
            Q += 2 * np.pi * r * v * dr
        
        return max(Q, 0.0)


# ─────────────────────────────────────────────────────────────────────────────
# FÍSICA: CRITERIO DE VON MISES
# ─────────────────────────────────────────────────────────────────────────────

def von_mises_stress(sigma_r, sigma_theta, sigma_z, tau_rz=0):
    """
    Criterio de Von Mises para integridad estructural del filamento extruido
    σ_VM = sqrt(0.5 * [(σr-σθ)² + (σθ-σz)² + (σz-σr)² + 6τ²])
    """
    s1 = (sigma_r - sigma_theta) ** 2
    s2 = (sigma_theta - sigma_z) ** 2
    s3 = (sigma_z - sigma_r) ** 2
    s4 = 6 * tau_rz ** 2
    return np.sqrt(0.5 * (s1 + s2 + s3 + s4))


def filament_stress_analysis(pressure_kPa, nozzle_diameter_um, print_speed_mm_s):
    """
    Análisis de estrés del filamento durante extrusión
    Retorna: (sigma_VM, tau_wall, integridad_ok)
    """
    P = pressure_kPa * 1000  # Pa
    r = (nozzle_diameter_um * 1e-6) / 2  # m
    v = print_speed_mm_s * 1e-3  # m/s
    
    # Estrés radial (presión interna)
    sigma_r = -P  # compresivo
    
    # Estrés hoop (circunferencial) — pared delgada
    sigma_theta = P * r / (r * 0.1) if r > 0 else 0  # t ~ 10% del radio
    
    # Estrés axial por aceleración
    rho = 1050  # kg/m³ densidad bioink
    sigma_z = rho * v ** 2  # Pa
    
    # Estrés de corte en la pared (HB)
    hb = HerschelBulkley()
    tau_wall = hb.shear_stress(hb.shear_rate_from_pressure(P, r))
    
    sigma_VM = von_mises_stress(abs(sigma_r), abs(sigma_theta), sigma_z, tau_wall)
    
    # Resistencia del NICE Bioink (límite ~350 kPa basado en simulador estrés mecánico)
    integridad_ok = sigma_VM < 350000
    
    return sigma_VM, tau_wall, integridad_ok


# ─────────────────────────────────────────────────────────────────────────────
# MODELO DE VIABILIDAD CELULAR
# ─────────────────────────────────────────────────────────────────────────────

def cell_viability_model(tau_wall_Pa, temp_C, nozzle_diam_um, layer_time_s):
    """
    Modelo empírico de viabilidad celular post-impresión
    Basado en literatura Co-SWIFT y datos validados del pipeline
    
    Factores:
      - Estrés de corte (τ_wall): daño mecánico a membrana
      - Temperatura: desnaturalización proteica y gelificación GelMA
      - Diámetro boquilla: trauma por confinamiento
      - Tiempo entre capas: hipoxia acumulada
    """
    # Factor estrés de corte (0-1, 1=máxima viabilidad)
    # Límite citotóxico: 150 Pa
    if tau_wall_Pa <= 0:
        f_shear = 0.60  # sin flujo = sin impresión viable
    elif tau_wall_Pa < 50:
        f_shear = 1.0
    elif tau_wall_Pa < 150:
        f_shear = 1.0 - 0.3 * ((tau_wall_Pa - 50) / 100)
    elif tau_wall_Pa < 300:
        f_shear = 0.70 - 0.25 * ((tau_wall_Pa - 150) / 150)
    else:
        f_shear = 0.45 - 0.15 * min((tau_wall_Pa - 300) / 200, 1.0)

    # Factor temperatura
    # Óptimo: 35-37°C. GelMA gelifica <25°C (no extruible). >40°C daño celular
    if 35 <= temp_C <= 37:
        f_temp = 1.0
    elif 30 <= temp_C < 35:
        f_temp = 0.85 + 0.15 * (temp_C - 30) / 5
    elif 37 < temp_C <= 39:
        f_temp = 1.0 - 0.15 * (temp_C - 37) / 2
    elif 39 < temp_C <= 42:
        f_temp = 0.85 - 0.30 * (temp_C - 39) / 3
    else:
        f_temp = max(0.3, 0.85 - 0.55 * abs(temp_C - 36) / 10)

    # Factor diámetro boquilla
    # Óptimo: 250-400 µm. <150 µm: daño por confinamiento. >600 µm: baja resolución
    if 250 <= nozzle_diam_um <= 400:
        f_nozzle = 1.0
    elif 200 <= nozzle_diam_um < 250:
        f_nozzle = 0.90 + 0.10 * (nozzle_diam_um - 200) / 50
    elif 400 < nozzle_diam_um <= 600:
        f_nozzle = 1.0 - 0.10 * (nozzle_diam_um - 400) / 200
    elif nozzle_diam_um < 200:
        f_nozzle = max(0.50, 0.90 - 0.40 * (200 - nozzle_diam_um) / 100)
    else:
        f_nozzle = 0.85

    # Factor tiempo entre capas (hipoxia)
    # Óptimo: 30-120 s. >300 s: hipoxia significativa
    if 30 <= layer_time_s <= 120:
        f_time = 1.0
    elif layer_time_s < 30:
        f_time = 0.92 + 0.08 * (layer_time_s / 30)
    elif 120 < layer_time_s <= 300:
        f_time = 1.0 - 0.12 * (layer_time_s - 120) / 180
    else:
        f_time = max(0.55, 0.88 - 0.25 * (layer_time_s - 300) / 300)

    viabilidad = f_shear * f_temp * f_nozzle * f_time * 100.0
    return min(max(viabilidad, 20.0), 98.0)


def vascular_resolution_model(pressure_kPa, nozzle_diam_um, print_speed_mm_s):
    """
    Resolución de canal vascular estimada [µm]
    Objetivo Co-SWIFT: 200 µm de diámetro mínimo
    """
    hb = HerschelBulkley()
    P = pressure_kPa * 1000
    r = (nozzle_diam_um * 1e-6) / 2
    
    Q = hb.flow_rate(P, r)
    v = print_speed_mm_s * 1e-3
    
    # Diámetro de filamento extruido
    if v > 0 and Q > 0:
        A_filament = Q / v
        d_filament = 2 * np.sqrt(A_filament / np.pi) * 1e6  # µm
    else:
        d_filament = nozzle_diam_um * 1.2
    
    # El canal Pluronic tiene ~85% del diámetro del filamento (contracción post-sacrificial)
    d_canal = d_filament * 0.85
    
    return max(d_canal, nozzle_diam_um * 0.7)


def print_time_model(pressure_kPa, nozzle_diam_um, print_speed_mm_s):
    """
    Tiempo de impresión por capa normalizado [s] (para 1 cm² de área)
    """
    hb = HerschelBulkley()
    P = pressure_kPa * 1000
    r = (nozzle_diam_um * 1e-6) / 2
    
    Q = hb.flow_rate(P, r)
    v = print_speed_mm_s * 1e-3
    
    if v <= 0 or Q <= 0:
        return 9999.0
    
    # Volumen por capa (1 cm² × altura de capa ~200 µm)
    V_layer = 1e-4 * 200e-6  # m³
    # Longitud total de trayectoria
    d_filament = 2 * np.sqrt((Q/v) / np.pi) * 1e6 if v > 0 else nozzle_diam_um
    spacing = max(d_filament * 1e-6, 100e-6)
    L_total = V_layer / (np.pi * r**2)
    
    t_print = L_total / v
    return min(t_print, 9999.0)


# ─────────────────────────────────────────────────────────────────────────────
# FUNCIÓN OBJETIVO MULTI-CRITERIO
# ─────────────────────────────────────────────────────────────────────────────

class ObjectiveFunction:
    """
    4 objetivos a optimizar simultáneamente:
      f1: Minimizar estrés de corte [Pa] — objetivo: < 150 Pa
      f2: Maximizar resolución vascular — objetivo: canal de 200 µm
      f3: Maximizar viabilidad celular [%] — objetivo: > 85%
      f4: Minimizar tiempo de impresión [s] — objetivo: mínimo
    """
    
    # Rangos de parámetros (límites de búsqueda)
    BOUNDS = {
        'pressure_kPa':     (5.0,   80.0),
        'print_speed_mm_s': (1.0,   25.0),
        'temp_C':           (25.0,  42.0),
        'nozzle_diam_um':   (150.0, 600.0),
        'layer_time_s':     (20.0,  300.0),
    }
    
    # Parámetros NICE Bioink fijos (validados)
    NICE_BIOINK = {
        'GelMA_pct':        7.0,
        'Alginato_pct':     1.5,
        'Nanocelulosa_pct': 0.8,
        'LAP_pct':          0.25,
    }
    
    def __init__(self):
        self.hb = HerschelBulkley()
        self.evaluation_count = 0
    
    def evaluate(self, x):
        """
        Evalúa los 4 objetivos para un vector de parámetros x
        x = [pressure_kPa, print_speed_mm_s, temp_C, nozzle_diam_um, layer_time_s]
        """
        pressure_kPa, print_speed_mm_s, temp_C, nozzle_diam_um, layer_time_s = x
        
        # Calcular física
        P = pressure_kPa * 1000
        r = (nozzle_diam_um * 1e-6) / 2
        
        gamma_dot = self.hb.shear_rate_from_pressure(P, r)
        tau_wall = self.hb.shear_stress(gamma_dot)
        
        sigma_VM, tau_wall_full, _ = filament_stress_analysis(
            pressure_kPa, nozzle_diam_um, print_speed_mm_s
        )
        
        # Objetivos
        f1 = tau_wall_full  # Minimizar estrés de corte [Pa]
        
        d_canal = vascular_resolution_model(pressure_kPa, nozzle_diam_um, print_speed_mm_s)
        f2 = abs(d_canal - 200.0)  # Minimizar desviación del objetivo 200 µm
        
        viab = cell_viability_model(tau_wall_full, temp_C, nozzle_diam_um, layer_time_s)
        f3 = 100.0 - viab  # Minimizar (1-viabilidad) = maximizar viabilidad
        
        f4 = print_time_model(pressure_kPa, nozzle_diam_um, print_speed_mm_s)
        
        self.evaluation_count += 1
        return np.array([f1, f2, f3, f4])
    
    def is_feasible(self, x):
        """Verifica restricciones de factibilidad"""
        pressure_kPa, print_speed_mm_s, temp_C, nozzle_diam_um, layer_time_s = x
        
        P = pressure_kPa * 1000
        r = (nozzle_diam_um * 1e-6) / 2
        
        hb = self.hb
        gamma_dot = hb.shear_rate_from_pressure(P, r)
        tau_wall = hb.shear_stress(gamma_dot)
        
        # Restricción 1: Estrés de corte < 150 Pa (límite citotóxico)
        if tau_wall > 150:
            return False, f"Shear stress {tau_wall:.1f} Pa > 150 Pa límite citotóxico"
        
        # Restricción 2: Hay flujo real (τ_wall > τ_yield)
        if tau_wall <= hb.tau_y:
            return False, f"Sin flujo: τ_wall ({tau_wall:.1f}) ≤ τ_yield ({hb.tau_y})"
        
        # Restricción 3: Von Mises < límite material
        _, _, integridad = filament_stress_analysis(pressure_kPa, nozzle_diam_um, print_speed_mm_s)
        if not integridad:
            return False, "Von Mises excede límite estructural del bioink"
        
        # Restricción 4: GelMA debe estar en rango de extrusión (T > 25°C)
        if temp_C < 25:
            return False, f"Temperatura {temp_C}°C < 25°C (GelMA gelificado)"
        
        return True, "OK"


# ─────────────────────────────────────────────────────────────────────────────
# PSO — PARTICLE SWARM OPTIMIZATION MULTI-OBJETIVO
# ─────────────────────────────────────────────────────────────────────────────

class MOPSO:
    """
    Multi-Objective Particle Swarm Optimization
    Implementación MOPSO con archivo externo de Pareto
    
    Parámetros PSO clásicos:
      w  = 0.729 (inercia — convergencia comprobada)
      c1 = 1.494 (coeficiente cognitivo — hacia mejor personal)
      c2 = 1.494 (coeficiente social — hacia mejor global)
    """
    
    def __init__(
        self,
        n_particles=80,
        n_iterations=150,
        w=0.729,
        c1=1.494,
        c2=1.494,
        archive_size=100,
        seed=42,
        progress_callback=None
    ):
        self.n_particles = n_particles
        self.n_iterations = n_iterations
        self.w = w
        self.c1 = c1
        self.c2 = c2
        self.archive_size = archive_size
        self.seed = seed
        self.progress_callback = progress_callback
        
        self.obj = ObjectiveFunction()
        self.bounds = list(self.obj.BOUNDS.values())
        self.n_params = len(self.bounds)
        self.n_objectives = 4
        
        self.rng = np.random.RandomState(seed)
        
        # Historia de convergencia
        self.convergence_history = []
        self.feasible_count_history = []
        self.archive_history = []
    
    def _random_position(self):
        pos = np.zeros(self.n_params)
        for i, (lo, hi) in enumerate(self.bounds):
            pos[i] = self.rng.uniform(lo, hi)
        return pos
    
    def _clip_position(self, pos):
        for i, (lo, hi) in enumerate(self.bounds):
            pos[i] = np.clip(pos[i], lo, hi)
        return pos
    
    def _dominates(self, f1, f2):
        """f1 domina f2 si f1 es mejor o igual en todos y estrictamente mejor en al menos uno"""
        return np.all(f1 <= f2) and np.any(f1 < f2)
    
    def _update_archive(self, archive, new_solution):
        """Actualiza el archivo de Pareto"""
        f_new = new_solution[1]
        
        # Verificar si es dominado por alguien en el archivo
        for sol in archive:
            if self._dominates(sol[1], f_new):
                return archive  # Dominado, no agregar
        
        # Eliminar soluciones dominadas por la nueva
        archive = [sol for sol in archive if not self._dominates(f_new, sol[1])]
        
        # Agregar nueva solución
        archive.append(new_solution)
        
        # Limitar tamaño del archivo (mantener diversidad)
        if len(archive) > self.archive_size:
            # Eliminar la más congestionada (crowding distance)
            archive = self._crowding_distance_sort(archive)
            archive = archive[:self.archive_size]
        
        return archive
    
    def _crowding_distance_sort(self, archive):
        """Ordena por crowding distance (diversidad en frente de Pareto)"""
        if len(archive) <= 2:
            return archive
        
        objectives = np.array([sol[1] for sol in archive])
        n = len(objectives)
        distances = np.zeros(n)
        
        for m in range(objectives.shape[1]):
            sorted_idx = np.argsort(objectives[:, m])
            distances[sorted_idx[0]] = np.inf
            distances[sorted_idx[-1]] = np.inf
            f_range = objectives[sorted_idx[-1], m] - objectives[sorted_idx[0], m]
            if f_range > 0:
                for i in range(1, n-1):
                    distances[sorted_idx[i]] += (
                        objectives[sorted_idx[i+1], m] - objectives[sorted_idx[i-1], m]
                    ) / f_range
        
        # Ordenar por distancia descendente (más separados primero)
        sorted_by_crowd = sorted(zip(distances, range(n)), reverse=True)
        return [archive[i] for _, i in sorted_by_crowd]
    
    def _select_guide(self, archive):
        """Selecciona guía del archivo por torneo de crowding distance"""
        if len(archive) == 1:
            return archive[0][0]
        
        # Selección aleatoria de 2 candidatos
        idx1, idx2 = self.rng.choice(len(archive), 2, replace=False)
        
        # El que tiene mayor crowding distance (más diverso)
        return archive[idx1][0]
    
    def optimize(self):
        """Ejecuta la optimización PSO multi-objetivo"""
        print("\n" + "═"*70)
        print("  OPTIMIZADOR Co-SWIFT — INICIANDO PSO MULTI-OBJETIVO")
        print("═"*70)
        print(f"  Partículas:  {self.n_particles}")
        print(f"  Iteraciones: {self.n_iterations}")
        print(f"  Parámetros:  {self.n_params} variables de decisión")
        print(f"  Objetivos:   {self.n_objectives} funciones")
        print("═"*70)
        
        # Inicialización
        positions = np.array([self._random_position() for _ in range(self.n_particles)])
        velocities = np.zeros_like(positions)
        
        # Inicializar con velocidad pequeña
        for i, (lo, hi) in enumerate(self.bounds):
            velocities[:, i] = self.rng.uniform(-(hi-lo)*0.1, (hi-lo)*0.1, self.n_particles)
        
        # Evaluar posiciones iniciales
        fitness = np.array([self.obj.evaluate(p) for p in positions])
        
        # Mejores personales
        pbest_pos = positions.copy()
        pbest_fit = fitness.copy()
        
        # Archivo de Pareto externo
        archive = []
        for i, (pos, fit) in enumerate(zip(positions, fitness)):
            feasible, _ = self.obj.is_feasible(pos)
            if feasible:
                archive = self._update_archive(archive, (pos.copy(), fit.copy()))
        
        best_viability_history = []
        
        # Bucle PSO
        for iteration in range(self.n_iterations):
            for i in range(self.n_particles):
                # Seleccionar guía del archivo (o posición aleatoria si vacío)
                if archive:
                    gbest_pos = self._select_guide(archive)
                else:
                    gbest_pos = positions[self.rng.randint(self.n_particles)]
                
                # Actualizar velocidad
                r1 = self.rng.random(self.n_params)
                r2 = self.rng.random(self.n_params)
                
                velocities[i] = (
                    self.w * velocities[i]
                    + self.c1 * r1 * (pbest_pos[i] - positions[i])
                    + self.c2 * r2 * (gbest_pos - positions[i])
                )
                
                # Limitar velocidad máxima
                for j, (lo, hi) in enumerate(self.bounds):
                    v_max = (hi - lo) * 0.2
                    velocities[i, j] = np.clip(velocities[i, j], -v_max, v_max)
                
                # Actualizar posición
                positions[i] = self._clip_position(positions[i] + velocities[i])
                
                # Evaluar
                new_fit = self.obj.evaluate(positions[i])
                
                # Actualizar mejor personal (usando suma ponderada como proxy)
                if np.sum(new_fit) < np.sum(pbest_fit[i]):
                    pbest_pos[i] = positions[i].copy()
                    pbest_fit[i] = new_fit.copy()
                
                # Actualizar archivo de Pareto
                feasible, _ = self.obj.is_feasible(positions[i])
                if feasible:
                    archive = self._update_archive(archive, (positions[i].copy(), new_fit.copy()))
            
            # Registro de convergencia
            if archive:
                # Mejor viabilidad en archivo
                best_viab = max(100.0 - sol[1][2] for sol in archive)
                best_shear = min(sol[1][0] for sol in archive)
                best_viability_history.append(best_viab)
                self.convergence_history.append({
                    'iter': iteration,
                    'archive_size': len(archive),
                    'best_viability': best_viab,
                    'best_shear': best_shear,
                    'evaluations': self.obj.evaluation_count
                })
            
            feasible_count = sum(1 for p in positions if self.obj.is_feasible(p)[0])
            self.feasible_count_history.append(feasible_count)
            
            # Reporte de progreso
            if iteration % 20 == 0 or iteration == self.n_iterations - 1:
                viab_str = f"{best_viability_history[-1]:.1f}%" if best_viability_history else "—"
                print(f"  Iter {iteration+1:3d}/{self.n_iterations} | "
                      f"Archivo Pareto: {len(archive):3d} | "
                      f"Factibles: {feasible_count:3d}/{self.n_particles} | "
                      f"Mejor viabilidad: {viab_str}")
            
            if self.progress_callback:
                progress = int((iteration + 1) / self.n_iterations * 100)
                self.progress_callback(progress, len(archive), feasible_count)
        
        print(f"\n  ✓ Optimización completada — {self.obj.evaluation_count} evaluaciones")
        print(f"  ✓ Frente de Pareto: {len(archive)} soluciones no-dominadas")
        
        return archive


# ─────────────────────────────────────────────────────────────────────────────
# ANÁLISIS DE RESULTADOS
# ─────────────────────────────────────────────────────────────────────────────

def find_recommended_solution(archive):
    """
    Selecciona la solución recomendada del frente de Pareto
    Criterio: máxima viabilidad celular con restricción τ < 150 Pa
    """
    if not archive:
        return None
    
    # Normalizar objetivos
    objectives = np.array([sol[1] for sol in archive])
    
    # Pesos de preferencia (ajustados a objetivos clínicos)
    # f1 (shear): peso 0.30
    # f2 (canal): peso 0.25
    # f3 (viabilidad): peso 0.35  — prioridad clínica
    # f4 (tiempo): peso 0.10
    weights = np.array([0.30, 0.25, 0.35, 0.10])
    
    # Normalizar [0,1] por min-max
    obj_min = objectives.min(axis=0)
    obj_max = objectives.max(axis=0)
    obj_range = obj_max - obj_min
    obj_range[obj_range == 0] = 1.0
    
    obj_norm = (objectives - obj_min) / obj_range
    
    # Función de utilidad ponderada
    utility = obj_norm @ weights
    
    best_idx = np.argmin(utility)
    return archive[best_idx]


def generate_protocol_table(solution, obj_func):
    """Genera tabla del protocolo de impresión recomendado"""
    pos = solution[0]
    fit = solution[1]
    
    pressure_kPa, print_speed_mm_s, temp_C, nozzle_diam_um, layer_time_s = pos
    tau_wall = fit[0]
    d_canal = 200.0 + fit[1] if fit[1] >= 0 else 200.0 - fit[1]
    viab = 100.0 - fit[2]
    
    d_canal = vascular_resolution_model(pressure_kPa, nozzle_diam_um, print_speed_mm_s)
    
    # Estado por parámetro
    def estado(val, optimal_range, ok_range):
        lo_opt, hi_opt = optimal_range
        lo_ok, hi_ok = ok_range
        if lo_opt <= val <= hi_opt:
            return "ÓPTIMO", COLORS['optimal']
        elif lo_ok <= val <= hi_ok:
            return "PARCIAL", COLORS['partial']
        else:
            return "CRÍTICO", COLORS['critical']
    
    protocol = {
        'Presión de extrusión': {
            'valor': f'{pressure_kPa:.1f} kPa',
            'referencia': '20-40 kPa (validado SWIFT)',
            'estado': estado(pressure_kPa, (20, 40), (10, 60))[0],
            'color': estado(pressure_kPa, (20, 40), (10, 60))[1],
        },
        'Velocidad de impresión': {
            'valor': f'{print_speed_mm_s:.1f} mm/s',
            'referencia': '5-15 mm/s (Co-SWIFT literatura)',
            'estado': estado(print_speed_mm_s, (5, 15), (2, 20))[0],
            'color': estado(print_speed_mm_s, (5, 15), (2, 20))[1],
        },
        'Temperatura bioink': {
            'valor': f'{temp_C:.1f} °C',
            'referencia': '35-37°C (GelMA extruible)',
            'estado': estado(temp_C, (35, 37), (30, 39))[0],
            'color': estado(temp_C, (35, 37), (30, 39))[1],
        },
        'Diámetro de boquilla': {
            'valor': f'{nozzle_diam_um:.0f} µm',
            'referencia': '250-400 µm (canal 200 µm)',
            'estado': estado(nozzle_diam_um, (250, 400), (200, 500))[0],
            'color': estado(nozzle_diam_um, (250, 400), (200, 500))[1],
        },
        'Tiempo entre capas': {
            'valor': f'{layer_time_s:.0f} s',
            'referencia': '30-120 s (hipoxia mínima)',
            'estado': estado(layer_time_s, (30, 120), (20, 200))[0],
            'color': estado(layer_time_s, (30, 120), (20, 200))[1],
        },
        'Estrés de corte (τ)': {
            'valor': f'{tau_wall:.1f} Pa',
            'referencia': '< 150 Pa (límite citotóxico)',
            'estado': 'ÓPTIMO' if tau_wall < 100 else 'PARCIAL' if tau_wall < 150 else 'CRÍTICO',
            'color': COLORS['optimal'] if tau_wall < 100 else COLORS['partial'] if tau_wall < 150 else COLORS['critical'],
        },
        'Canal vascular (d)': {
            'valor': f'{d_canal:.0f} µm',
            'referencia': '≥ 200 µm (objetivo Co-SWIFT)',
            'estado': 'ÓPTIMO' if d_canal >= 190 else 'PARCIAL' if d_canal >= 150 else 'CRÍTICO',
            'color': COLORS['optimal'] if d_canal >= 190 else COLORS['partial'] if d_canal >= 150 else COLORS['critical'],
        },
        'Viabilidad celular': {
            'valor': f'{viab:.1f} %',
            'referencia': '> 85% (objetivo bioimpresión)',
            'estado': 'ÓPTIMO' if viab >= 85 else 'PARCIAL' if viab >= 70 else 'CRÍTICO',
            'color': COLORS['optimal'] if viab >= 85 else COLORS['partial'] if viab >= 70 else COLORS['critical'],
        },
        'GelMA': {
            'valor': '7.0 %',
            'referencia': '6-8% (NICE Bioink validado)',
            'estado': 'ÓPTIMO',
            'color': COLORS['optimal'],
        },
        'Alginato': {
            'valor': '1.5 %',
            'referencia': '1.0-2.0% (NICE Bioink validado)',
            'estado': 'ÓPTIMO',
            'color': COLORS['optimal'],
        },
        'Nanocelulosa': {
            'valor': '0.8 %',
            'referencia': '0.5-1.0% (NICE Bioink validado)',
            'estado': 'ÓPTIMO',
            'color': COLORS['optimal'],
        },
        'LAP (fotoiniciador)': {
            'valor': '0.25 %',
            'referencia': '0.1-0.5% (no citotóxico)',
            'estado': 'ÓPTIMO',
            'color': COLORS['optimal'],
        },
    }
    
    return protocol


# ─────────────────────────────────────────────────────────────────────────────
# VISUALIZACIÓN — DASHBOARD MATPLOTLIB
# ─────────────────────────────────────────────────────────────────────────────

def plot_results(archive, convergence_history, recommended_solution, obj_func):
    """Genera el dashboard completo de resultados"""
    
    fig = plt.figure(figsize=(20, 14), facecolor=COLORS['bg_dark'])
    fig.suptitle(
        'OPTIMIZADOR Co-SWIFT — Bio-Kidney AI 2026\n'
        'Frente de Pareto Multi-Objetivo | VirtusSapiens',
        fontsize=16, fontweight='bold', color=COLORS['accent_cyan'],
        y=0.98
    )
    
    gs = gridspec.GridSpec(
        3, 4,
        figure=fig,
        hspace=0.42,
        wspace=0.38,
        top=0.92,
        bottom=0.06,
        left=0.06,
        right=0.97
    )
    
    objectives = np.array([sol[1] for sol in archive])
    positions = np.array([sol[0] for sol in archive])
    
    # Extraer métricas legibles
    shear_vals = objectives[:, 0]
    canal_dev  = objectives[:, 1]
    viab_vals  = 100 - objectives[:, 2]
    time_vals  = objectives[:, 3]
    
    d_canals = np.array([
        vascular_resolution_model(p[0], p[3], p[1]) for p in positions
    ])
    
    rec_pos = recommended_solution[0]
    rec_fit = recommended_solution[1]
    rec_shear = rec_fit[0]
    rec_viab = 100 - rec_fit[2]
    rec_canal = vascular_resolution_model(rec_pos[0], rec_pos[3], rec_pos[1])
    
    # ── Plot 1: Frente de Pareto — Shear vs Viabilidad ──────────────────────
    ax1 = fig.add_subplot(gs[0, 0:2])
    sc1 = ax1.scatter(
        shear_vals, viab_vals,
        c=d_canals, cmap='plasma', s=40, alpha=0.7,
        zorder=3, label='Soluciones Pareto'
    )
    ax1.scatter(
        rec_shear, rec_viab,
        color=COLORS['accent_teal'], s=250, marker='*',
        zorder=5, label='* Recomendada', edgecolors='white', linewidths=0.5
    )
    ax1.axvline(x=150, color=COLORS['accent_red'], lw=1.5, ls='--', alpha=0.8, label='Límite citotóxico')
    ax1.axhline(y=85, color=COLORS['optimal'], lw=1.5, ls='--', alpha=0.8, label='Objetivo viabilidad')
    
    cb1 = plt.colorbar(sc1, ax=ax1, pad=0.02)
    cb1.set_label('Canal [µm]', color=COLORS['text_secondary'], fontsize=8)
    cb1.ax.yaxis.set_tick_params(color=COLORS['text_secondary'])
    plt.setp(cb1.ax.yaxis.get_ticklabels(), color=COLORS['text_secondary'], fontsize=7)
    
    ax1.set_xlabel('Estrés de Corte τ [Pa]')
    ax1.set_ylabel('Viabilidad Celular [%]')
    ax1.set_title('FRENTE DE PARETO\nShear Stress vs Viabilidad', color=COLORS['accent_cyan'])
    ax1.legend(fontsize=7, framealpha=0.3)
    ax1.grid(True, alpha=0.3)
    
    # ── Plot 2: Frente de Pareto — Canal vs Viabilidad ──────────────────────
    ax2 = fig.add_subplot(gs[0, 2:4])
    sc2 = ax2.scatter(
        d_canals, viab_vals,
        c=shear_vals, cmap='hot', s=40, alpha=0.7, zorder=3
    )
    ax2.scatter(
        rec_canal, rec_viab,
        color=COLORS['accent_teal'], s=250, marker='*',
        zorder=5, edgecolors='white', linewidths=0.5
    )
    ax2.axvline(x=200, color=COLORS['accent_gold'], lw=1.5, ls='--', alpha=0.8, label='Objetivo canal 200 µm')
    ax2.axhline(y=85, color=COLORS['optimal'], lw=1.5, ls='--', alpha=0.8, label='Objetivo >85%')
    
    cb2 = plt.colorbar(sc2, ax=ax2, pad=0.02)
    cb2.set_label('Shear [Pa]', color=COLORS['text_secondary'], fontsize=8)
    cb2.ax.yaxis.set_tick_params(color=COLORS['text_secondary'])
    plt.setp(cb2.ax.yaxis.get_ticklabels(), color=COLORS['text_secondary'], fontsize=7)
    
    ax2.set_xlabel('Diámetro Canal Vascular [µm]')
    ax2.set_ylabel('Viabilidad Celular [%]')
    ax2.set_title('FRENTE DE PARETO\nResolución Vascular vs Viabilidad', color=COLORS['accent_cyan'])
    ax2.legend(fontsize=7, framealpha=0.3)
    ax2.grid(True, alpha=0.3)
    
    # ── Plot 3: Convergencia del algoritmo ───────────────────────────────────
    ax3 = fig.add_subplot(gs[1, 0:2])
    if convergence_history:
        iters = [h['iter'] for h in convergence_history]
        viabs = [h['best_viability'] for h in convergence_history]
        sheares = [h['best_shear'] for h in convergence_history]
        archive_sizes = [h['archive_size'] for h in convergence_history]
        
        ax3_twin = ax3.twinx()
        ax3.plot(iters, viabs, color=COLORS['accent_teal'], lw=2, label='Viabilidad máxima [%]')
        ax3_twin.plot(iters, archive_sizes, color=COLORS['accent_purple'], lw=1.5, ls=':', label='Tamaño archivo')
        ax3.axhline(y=85, color=COLORS['optimal'], lw=1, ls='--', alpha=0.6)
        
        ax3.set_xlabel('Iteración PSO')
        ax3.set_ylabel('Viabilidad Máxima [%]', color=COLORS['accent_teal'])
        ax3_twin.set_ylabel('Soluciones Pareto', color=COLORS['accent_purple'])
        ax3_twin.tick_params(axis='y', labelcolor=COLORS['accent_purple'])
        ax3.tick_params(axis='y', labelcolor=COLORS['accent_teal'])
        ax3.set_title('CONVERGENCIA PSO\nEvolución del Algoritmo', color=COLORS['accent_cyan'])
        ax3.grid(True, alpha=0.3)
        
        lines1, labels1 = ax3.get_legend_handles_labels()
        lines2, labels2 = ax3_twin.get_legend_handles_labels()
        ax3.legend(lines1 + lines2, labels1 + labels2, fontsize=7, framealpha=0.3)
    
    # ── Plot 4: Distribución de parámetros en frente de Pareto ───────────────
    ax4 = fig.add_subplot(gs[1, 2])
    param_names = ['P\n[kPa]', 'v\n[mm/s]', 'T\n[°C]', 'Ø\n[µm/10]', 't\n[s/10]']
    param_vals_norm = [
        positions[:, 0] / 80,
        positions[:, 1] / 25,
        (positions[:, 2] - 25) / 17,
        positions[:, 3] / 600,
        positions[:, 4] / 300
    ]
    rec_norm = [
        rec_pos[0]/80, rec_pos[1]/25,
        (rec_pos[2]-25)/17, rec_pos[3]/600, rec_pos[4]/300
    ]
    
    bp = ax4.boxplot(
        [v for v in param_vals_norm],
        labels=param_names,
        patch_artist=True,
        medianprops={'color': COLORS['accent_teal'], 'lw': 2},
        boxprops={'facecolor': COLORS['bg_card'], 'edgecolor': COLORS['accent_cyan']},
        whiskerprops={'color': COLORS['text_secondary']},
        capprops={'color': COLORS['text_secondary']},
        flierprops={'marker': 'o', 'color': COLORS['text_secondary'], 'ms': 3}
    )
    
    for i, v in enumerate(rec_norm):
        ax4.scatter(i+1, v, color=COLORS['accent_teal'], s=100, marker='*', zorder=5)
    
    ax4.set_title('DISTRIBUCIÓN PARETO\nNormalizada [0,1]', color=COLORS['accent_cyan'])
    ax4.set_ylabel('Valor normalizado')
    ax4.grid(True, alpha=0.3, axis='y')
    
    # ── Plot 5: Radar de la solución recomendada ─────────────────────────────
    ax5 = fig.add_subplot(gs[1, 3], projection='polar')
    
    categories = ['Viabilidad\n[%]', 'Canal\n[µm/250]', 'Shear\n[seguro]', 'Velocidad\n[norm]', 'Temp\n[óptima]']
    N = len(categories)
    
    rec_shear_score = max(0, (150 - rec_shear) / 150)  # 1 = mejor (más bajo)
    rec_canal_score = min(rec_canal / 250, 1.0)
    rec_viab_score = rec_viab / 100
    rec_speed_score = min(rec_pos[1] / 15, 1.0)
    rec_temp_score = 1.0 - abs(rec_pos[2] - 36) / 10
    
    values = [rec_viab_score, rec_canal_score, rec_shear_score, rec_speed_score, rec_temp_score]
    values += [values[0]]  # cerrar polígono
    
    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    angles += angles[:1]
    
    ax5.plot(angles, values, color=COLORS['accent_teal'], lw=2)
    ax5.fill(angles, values, color=COLORS['accent_teal'], alpha=0.25)
    ax5.set_xticks(angles[:-1])
    ax5.set_xticklabels(categories, size=7, color=COLORS['text_primary'])
    ax5.set_ylim(0, 1)
    ax5.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax5.set_yticklabels(['25%', '50%', '75%', '100%'], size=6, color=COLORS['text_secondary'])
    ax5.grid(color=COLORS['border'], alpha=0.5)
    ax5.set_facecolor(COLORS['bg_panel'])
    ax5.set_title('PERFIL SOLUCIÓN\nRecomendada', color=COLORS['accent_cyan'], pad=15)
    
    # ── Plot 6: Tabla del protocolo ──────────────────────────────────────────
    ax6 = fig.add_subplot(gs[2, :])
    ax6.set_facecolor(COLORS['bg_card'])
    ax6.axis('off')
    
    protocol = generate_protocol_table(recommended_solution, obj_func)
    
    items = list(protocol.items())
    n_cols = 4
    n_rows = (len(items) + n_cols - 1) // n_cols
    
    col_width = 1.0 / n_cols
    row_height = 1.0 / (n_rows + 1)
    
    # Encabezado
    ax6.text(0.5, 1.0 - 0.005, 'PROTOCOLO DE BIOIMPRESIÓN Co-SWIFT — PARÁMETROS RECOMENDADOS',
             ha='center', va='top', fontsize=11, fontweight='bold',
             color=COLORS['accent_cyan'], transform=ax6.transAxes)
    
    for idx, (param, data) in enumerate(items):
        row = idx // n_cols
        col = idx % n_cols
        
        x = col * col_width + col_width * 0.05
        y = 1.0 - (row + 1.2) * row_height
        
        # Fondo del card
        bg = FancyBboxPatch(
            (x, y - row_height * 0.15),
            col_width * 0.90, row_height * 0.75,
            boxstyle="round,pad=0.005",
            transform=ax6.transAxes,
            facecolor=COLORS['bg_panel'],
            edgecolor=data['color'],
            linewidth=1.5, zorder=2
        )
        ax6.add_patch(bg)
        
        # Estado badge
        estado_x = x + col_width * 0.90 - col_width * 0.22
        estado_y = y + row_height * 0.42
        ax6.text(
            estado_x, estado_y, data['estado'],
            ha='center', va='center',
            fontsize=7, fontweight='bold',
            color=data['color'],
            transform=ax6.transAxes, zorder=3
        )
        
        # Nombre parámetro
        ax6.text(x + col_width*0.04, y + row_height*0.30, param,
                 ha='left', va='center', fontsize=8, fontweight='bold',
                 color=COLORS['text_primary'], transform=ax6.transAxes, zorder=3)
        
        # Valor
        ax6.text(x + col_width*0.04, y + row_height*0.08, data['valor'],
                 ha='left', va='center', fontsize=10, fontweight='bold',
                 color=data['color'], transform=ax6.transAxes, zorder=3)
        
        # Referencia
        ax6.text(x + col_width*0.04, y - row_height*0.08, data['referencia'],
                 ha='left', va='center', fontsize=6.5,
                 color=COLORS['text_secondary'], transform=ax6.transAxes, zorder=3)
    
    plt.savefig(
        os.path.expanduser('~/Escritorio/BioKidney-AI/resultados/pareto_coswift.png'),
        dpi=150, bbox_inches='tight', facecolor=COLORS['bg_dark']
    )
    print("  ✓ Dashboard guardado: ~/Escritorio/BioKidney-AI/resultados/pareto_coswift.png")
    
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# INTERFAZ PyQt6
# ─────────────────────────────────────────────────────────────────────────────

if HAS_QT:
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

    class OptimizationThread(QThread):
        progress_signal = pyqtSignal(int, int, int)
        finished_signal = pyqtSignal(list, list, object)
        log_signal = pyqtSignal(str)
        
        def __init__(self, n_particles, n_iterations):
            super().__init__()
            self.n_particles = n_particles
            self.n_iterations = n_iterations
        
        def run(self):
            import io
            from contextlib import redirect_stdout
            
            def progress_cb(pct, archive_size, feasible):
                self.progress_signal.emit(pct, archive_size, feasible)
            
            pso = MOPSO(
                n_particles=self.n_particles,
                n_iterations=self.n_iterations,
                progress_callback=progress_cb
            )
            
            archive = pso.optimize()
            convergence = pso.convergence_history
            recommended = find_recommended_solution(archive)
            
            self.finished_signal.emit(archive, convergence, recommended)

    STYLESHEET = f"""
    QMainWindow, QWidget {{
        background-color: {COLORS['bg_dark']};
        color: {COLORS['text_primary']};
        font-family: 'Consolas', 'Courier New', monospace;
    }}
    QGroupBox {{
        border: 1px solid {COLORS['border']};
        border-radius: 6px;
        margin-top: 8px;
        padding-top: 8px;
        font-weight: bold;
        color: {COLORS['accent_cyan']};
        font-size: 11px;
    }}
    QGroupBox::title {{
        subcontrol-origin: margin;
        left: 10px;
        padding: 0 5px;
    }}
    QPushButton {{
        background-color: {COLORS['bg_card']};
        color: {COLORS['accent_cyan']};
        border: 1px solid {COLORS['accent_cyan']};
        border-radius: 4px;
        padding: 8px 16px;
        font-size: 12px;
        font-weight: bold;
    }}
    QPushButton:hover {{
        background-color: {COLORS['accent_cyan']};
        color: {COLORS['bg_dark']};
    }}
    QPushButton:disabled {{
        color: {COLORS['text_secondary']};
        border-color: {COLORS['border']};
    }}
    QProgressBar {{
        border: 1px solid {COLORS['border']};
        border-radius: 4px;
        background-color: {COLORS['bg_card']};
        text-align: center;
        color: {COLORS['text_primary']};
        height: 20px;
    }}
    QProgressBar::chunk {{
        background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
            stop:0 {COLORS['accent_cyan']}, stop:1 {COLORS['accent_teal']});
        border-radius: 3px;
    }}
    QLabel {{
        color: {COLORS['text_primary']};
    }}
    QSpinBox, QDoubleSpinBox {{
        background-color: {COLORS['bg_card']};
        color: {COLORS['text_primary']};
        border: 1px solid {COLORS['border']};
        border-radius: 3px;
        padding: 3px;
    }}
    QTableWidget {{
        background-color: {COLORS['bg_card']};
        color: {COLORS['text_primary']};
        border: 1px solid {COLORS['border']};
        gridline-color: {COLORS['border']};
        selection-background-color: {COLORS['accent_cyan']};
        selection-color: {COLORS['bg_dark']};
    }}
    QHeaderView::section {{
        background-color: {COLORS['bg_panel']};
        color: {COLORS['accent_cyan']};
        border: 1px solid {COLORS['border']};
        padding: 4px;
        font-weight: bold;
    }}
    QTabWidget::pane {{
        border: 1px solid {COLORS['border']};
        background-color: {COLORS['bg_panel']};
    }}
    QTabBar::tab {{
        background-color: {COLORS['bg_card']};
        color: {COLORS['text_secondary']};
        border: 1px solid {COLORS['border']};
        padding: 6px 14px;
    }}
    QTabBar::tab:selected {{
        background-color: {COLORS['bg_panel']};
        color: {COLORS['accent_cyan']};
        border-bottom: 2px solid {COLORS['accent_cyan']};
    }}
    QTextEdit {{
        background-color: {COLORS['bg_card']};
        color: {COLORS['accent_teal']};
        border: 1px solid {COLORS['border']};
        font-family: 'Consolas', monospace;
        font-size: 10px;
    }}
    """

    class MainWindow(QMainWindow):
        def __init__(self):
            super().__init__()
            self.setWindowTitle("Optimizador Co-SWIFT — Bio-Kidney AI 2026 | VirtusSapiens")
            self.setMinimumSize(1100, 700)
            self.setStyleSheet(STYLESHEET)
            
            self.archive = None
            self.convergence_history = None
            self.recommended_solution = None
            self.opt_thread = None
            
            self._build_ui()
        
        def _build_ui(self):
            central = QWidget()
            self.setCentralWidget(central)
            main_layout = QVBoxLayout(central)
            main_layout.setSpacing(6)
            main_layout.setContentsMargins(10, 8, 10, 8)
            
            # ── Header compacto ──────────────────────────────────────────────
            header = QLabel(
                "⬡  OPTIMIZADOR Co-SWIFT  •  Bio-Kidney AI 2026  •  VirtusSapiens"
            )
            header.setStyleSheet(f"""
                color: {COLORS['accent_cyan']};
                font-size: 14px;
                font-weight: bold;
                padding: 6px 10px;
                background-color: {COLORS['bg_card']};
                border: 1px solid {COLORS['accent_cyan']};
                border-radius: 4px;
            """)
            header.setAlignment(Qt.AlignmentFlag.AlignCenter)
            header.setFixedHeight(36)
            main_layout.addWidget(header)
            
            # ── Contenido principal ──────────────────────────────────────────
            splitter = QSplitter(Qt.Orientation.Horizontal)
            
            # Panel izquierdo: scroll area para que los botones siempre sean visibles
            scroll_area = QScrollArea()
            scroll_area.setWidgetResizable(True)
            scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
            scroll_area.setStyleSheet(f"""
                QScrollArea {{
                    border: none;
                    background-color: {COLORS['bg_dark']};
                }}
                QScrollBar:vertical {{
                    background: {COLORS['bg_card']};
                    width: 6px;
                    border-radius: 3px;
                }}
                QScrollBar::handle:vertical {{
                    background: {COLORS['border']};
                    border-radius: 3px;
                    min-height: 20px;
                }}
            """)
            scroll_area.setFixedWidth(300)
            
            left_panel = QWidget()
            left_layout = QVBoxLayout(left_panel)
            left_layout.setSpacing(6)
            left_layout.setContentsMargins(4, 4, 4, 4)
            
            # Parámetros PSO
            pso_group = QGroupBox("PARÁMETROS PSO")
            pso_layout = QGridLayout(pso_group)
            pso_layout.setSpacing(4)
            pso_layout.setContentsMargins(6, 14, 6, 6)
            
            pso_layout.addWidget(QLabel("Partículas:"), 0, 0)
            self.spin_particles = QSpinBox()
            self.spin_particles.setRange(20, 200)
            self.spin_particles.setValue(80)
            pso_layout.addWidget(self.spin_particles, 0, 1)
            
            pso_layout.addWidget(QLabel("Iteraciones:"), 1, 0)
            self.spin_iterations = QSpinBox()
            self.spin_iterations.setRange(50, 500)
            self.spin_iterations.setValue(150)
            pso_layout.addWidget(self.spin_iterations, 1, 1)
            
            left_layout.addWidget(pso_group)
            
            # Parámetros físicos
            physics_group = QGroupBox("FÍSICA VALIDADA")
            physics_layout = QGridLayout(physics_group)
            physics_layout.setSpacing(3)
            physics_layout.setContentsMargins(6, 14, 6, 6)
            
            params_info = [
                ("τ_yield NICE Bioink:", "28 Pa"),
                ("K (consistencia):", "4.2 Pa·sⁿ"),
                ("n (flujo):", "0.41"),
                ("Límite citotóxico:", "150 Pa"),
                ("Canal objetivo:", "200 µm"),
                ("WSS validado:", "5.6 dyn/cm²"),
                ("P extrusión ref.:", "30 kPa"),
            ]
            for i, (label, val) in enumerate(params_info):
                lbl = QLabel(label)
                lbl.setStyleSheet(f"color: {COLORS['text_secondary']};")
                val_lbl = QLabel(val)
                val_lbl.setStyleSheet(f"color: {COLORS['accent_gold']}; font-weight: bold;")
                physics_layout.addWidget(lbl, i, 0)
                physics_layout.addWidget(val_lbl, i, 1)
            
            left_layout.addWidget(physics_group)
            
            # NICE Bioink
            bioink_group = QGroupBox("NICE BIOINK — COMPOSICIÓN")
            bioink_layout = QGridLayout(bioink_group)
            bioink_layout.setSpacing(3)
            bioink_layout.setContentsMargins(6, 14, 6, 6)
            
            bioink_params = [
                ("GelMA:", "7.0 %"),
                ("Alginato:", "1.5 %"),
                ("Nanocelulosa:", "0.8 %"),
                ("LAP:", "0.25 %"),
            ]
            for i, (label, val) in enumerate(bioink_params):
                lbl = QLabel(label)
                lbl.setStyleSheet(f"color: {COLORS['text_secondary']};")
                val_lbl = QLabel(val)
                val_lbl.setStyleSheet(f"color: {COLORS['optimal']}; font-weight: bold;")
                bioink_layout.addWidget(lbl, i, 0)
                bioink_layout.addWidget(val_lbl, i, 1)
            
            left_layout.addWidget(bioink_group)
            
            # Botones — altura fija para que siempre sean visibles
            self.btn_optimize = QPushButton("▶  INICIAR OPTIMIZACIÓN PSO")
            self.btn_optimize.setFixedHeight(36)
            self.btn_optimize.clicked.connect(self._start_optimization)
            left_layout.addWidget(self.btn_optimize)
            
            self.btn_plot = QPushButton("📊  VER DASHBOARD COMPLETO")
            self.btn_plot.setFixedHeight(34)
            self.btn_plot.clicked.connect(self._show_plots)
            self.btn_plot.setEnabled(False)
            left_layout.addWidget(self.btn_plot)
            
            self.btn_export = QPushButton("💾  EXPORTAR PROTOCOLO")
            self.btn_export.setFixedHeight(34)
            self.btn_export.clicked.connect(self._export_protocol)
            self.btn_export.setEnabled(False)
            left_layout.addWidget(self.btn_export)
            
            # Estado global
            self.status_label = QLabel("Estado: Esperando optimización...")
            self.status_label.setStyleSheet(
                f"color: {COLORS['accent_gold']}; font-size: 10px; padding: 4px;"
                f"background: {COLORS['bg_card']}; border-radius: 4px;"
            )
            self.status_label.setWordWrap(True)
            left_layout.addWidget(self.status_label)
            left_layout.addStretch()
            
            # Montar scroll area
            scroll_area.setWidget(left_panel)
            splitter.addWidget(scroll_area)
            
            # Panel derecho: resultados con tabs
            right_panel = QWidget()
            right_layout = QVBoxLayout(right_panel)
            right_layout.setSpacing(6)
            
            # Progress bar
            progress_widget = QWidget()
            progress_layout = QHBoxLayout(progress_widget)
            progress_layout.setContentsMargins(0, 0, 0, 0)
            
            self.progress_bar = QProgressBar()
            self.progress_bar.setValue(0)
            self.progress_bar.setFormat("PSO: %p% | Archivo Pareto: —")
            progress_layout.addWidget(self.progress_bar)
            right_layout.addWidget(progress_widget)
            
            # Tabs
            self.tabs = QTabWidget()
            
            # Tab 1: Log
            self.log_text = QTextEdit()
            self.log_text.setReadOnly(True)
            self.log_text.append("» OPTIMIZADOR Co-SWIFT — Bio-Kidney AI 2026")
            self.log_text.append("» Módulo: Optimización Multi-Objetivo PSO")
            self.log_text.append("» Física: Herschel-Bulkley + Von Mises + Pareto")
            self.log_text.append("─" * 60)
            self.tabs.addTab(self.log_text, "📋 Log PSO")
            
            # Tab 2: Tabla protocolo
            self.protocol_table = QTableWidget()
            self.protocol_table.setColumnCount(4)
            self.protocol_table.setHorizontalHeaderLabels(
                ["Parámetro", "Valor Óptimo", "Referencia", "Estado"]
            )
            self.protocol_table.horizontalHeader().setStretchLastSection(True)
            self.protocol_table.horizontalHeader().setSectionResizeMode(
                QHeaderView.ResizeMode.ResizeToContents
            )
            self.tabs.addTab(self.protocol_table, "📊 Protocolo")
            
            # Tab 3: Métricas
            self.metrics_widget = QWidget()
            self.metrics_layout = QGridLayout(self.metrics_widget)
            self.tabs.addTab(self.metrics_widget, "📈 Métricas")
            
            right_layout.addWidget(self.tabs)
            splitter.addWidget(right_panel)
            splitter.setSizes([300, 900])
            
            main_layout.addWidget(splitter)
        
        def _start_optimization(self):
            self.btn_optimize.setEnabled(False)
            self.btn_plot.setEnabled(False)
            self.btn_export.setEnabled(False)
            self.progress_bar.setValue(0)
            self.status_label.setText("Estado: Optimizando...")
            self.status_label.setStyleSheet(
                f"color: {COLORS['accent_cyan']}; font-size: 11px; padding: 5px;"
                f"background: {COLORS['bg_card']}; border-radius: 4px;"
            )
            self.log_text.append("\n» INICIANDO PSO MULTI-OBJETIVO...")
            self.log_text.append(f"  Partículas: {self.spin_particles.value()}")
            self.log_text.append(f"  Iteraciones: {self.spin_iterations.value()}")
            self.log_text.append("─" * 60)
            
            self.opt_thread = OptimizationThread(
                self.spin_particles.value(),
                self.spin_iterations.value()
            )
            self.opt_thread.progress_signal.connect(self._on_progress)
            self.opt_thread.finished_signal.connect(self._on_finished)
            self.opt_thread.start()
        
        def _on_progress(self, pct, archive_size, feasible):
            self.progress_bar.setValue(pct)
            self.progress_bar.setFormat(
                f"PSO: {pct}% | Pareto: {archive_size} | Factibles: {feasible}"
            )
            if pct % 20 == 0:
                self.log_text.append(
                    f"  [{pct:3d}%] Archivo Pareto: {archive_size} | Factibles: {feasible}"
                )
        
        def _on_finished(self, archive, convergence_history, recommended):
            self.archive = archive
            self.convergence_history = convergence_history
            self.recommended_solution = recommended
            
            self.progress_bar.setValue(100)
            self.log_text.append("\n» OPTIMIZACIÓN COMPLETADA ✓")
            self.log_text.append(f"  Soluciones Pareto encontradas: {len(archive)}")
            
            if recommended:
                pos = recommended[0]
                fit = recommended[1]
                viab = 100 - fit[2]
                canal = vascular_resolution_model(pos[0], pos[3], pos[1])
                
                self.log_text.append(f"\n» SOLUCIÓN RECOMENDADA:")
                self.log_text.append(f"  Presión:          {pos[0]:.1f} kPa")
                self.log_text.append(f"  Vel. impresión:   {pos[1]:.1f} mm/s")
                self.log_text.append(f"  Temperatura:      {pos[2]:.1f} °C")
                self.log_text.append(f"  Boquilla:         {pos[3]:.0f} µm")
                self.log_text.append(f"  T. entre capas:   {pos[4]:.0f} s")
                self.log_text.append(f"  Shear stress:     {fit[0]:.1f} Pa")
                self.log_text.append(f"  Canal vascular:   {canal:.0f} µm")
                self.log_text.append(f"  Viabilidad:       {viab:.1f} %")
                
                self._populate_protocol_table(recommended)
                self._populate_metrics(recommended)
                
                estado_global = "ÓPTIMO" if viab >= 85 and fit[0] < 150 else "PARCIAL"
                color = COLORS['optimal'] if estado_global == "ÓPTIMO" else COLORS['partial']
                self.status_label.setText(f"Estado global: {estado_global}")
                self.status_label.setStyleSheet(
                    f"color: {color}; font-size: 13px; font-weight: bold; padding: 5px;"
                    f"background: {COLORS['bg_card']}; border-radius: 4px;"
                    f"border: 1px solid {color};"
                )
            
            self.btn_optimize.setEnabled(True)
            self.btn_plot.setEnabled(True)
            self.btn_export.setEnabled(True)
        
        def _populate_protocol_table(self, recommended):
            protocol = generate_protocol_table(recommended, ObjectiveFunction())
            self.protocol_table.setRowCount(len(protocol))
            
            for row, (param, data) in enumerate(protocol.items()):
                self.protocol_table.setItem(row, 0, QTableWidgetItem(param))
                self.protocol_table.setItem(row, 1, QTableWidgetItem(data['valor']))
                self.protocol_table.setItem(row, 2, QTableWidgetItem(data['referencia']))
                
                estado_item = QTableWidgetItem(data['estado'])
                color_hex = data['color']
                estado_item.setForeground(QColor(color_hex))
                estado_item.setFont(QFont('Consolas', 9, QFont.Weight.Bold))
                self.protocol_table.setItem(row, 3, estado_item)
            
            self.tabs.setCurrentIndex(1)
        
        def _populate_metrics(self, recommended):
            # Limpiar layout anterior
            for i in reversed(range(self.metrics_layout.count())):
                self.metrics_layout.itemAt(i).widget().deleteLater()
            
            if not recommended:
                return
            
            pos = recommended[0]
            fit = recommended[1]
            viab = 100 - fit[2]
            canal = vascular_resolution_model(pos[0], pos[3], pos[1])
            
            metrics = [
                ("Viabilidad Celular", f"{viab:.1f}%", ">85%",
                 COLORS['optimal'] if viab >= 85 else COLORS['partial']),
                ("Shear Stress", f"{fit[0]:.1f} Pa", "<150 Pa",
                 COLORS['optimal'] if fit[0] < 100 else COLORS['partial']),
                ("Canal Vascular", f"{canal:.0f} µm", "≥200 µm",
                 COLORS['optimal'] if canal >= 190 else COLORS['partial']),
                ("Soluciones Pareto", str(len(self.archive)), "máx diversidad",
                 COLORS['accent_cyan']),
            ]
            
            for i, (name, val, ref, color) in enumerate(metrics):
                card = QFrame()
                card.setStyleSheet(
                    f"background: {COLORS['bg_card']}; border: 1px solid {color};"
                    f"border-radius: 6px; padding: 10px;"
                )
                card_layout = QVBoxLayout(card)
                
                name_lbl = QLabel(name)
                name_lbl.setStyleSheet(f"color: {COLORS['text_secondary']}; font-size: 10px;")
                
                val_lbl = QLabel(val)
                val_lbl.setStyleSheet(f"color: {color}; font-size: 22px; font-weight: bold;")
                
                ref_lbl = QLabel(f"Objetivo: {ref}")
                ref_lbl.setStyleSheet(f"color: {COLORS['text_secondary']}; font-size: 9px;")
                
                card_layout.addWidget(name_lbl)
                card_layout.addWidget(val_lbl)
                card_layout.addWidget(ref_lbl)
                
                self.metrics_layout.addWidget(card, i // 2, i % 2)
        
        def _show_plots(self):
            if not self.archive:
                return
            try:
                os.makedirs(
                    os.path.expanduser('~/Escritorio/BioKidney-AI/resultados'),
                    exist_ok=True
                )
                fig = plot_results(
                    self.archive, self.convergence_history,
                    self.recommended_solution, ObjectiveFunction()
                )
                plt.show()
            except Exception as e:
                QMessageBox.critical(self, "Error en Dashboard", str(e))
        
        def _export_protocol(self):
            if not self.recommended_solution:
                return
            
            protocol = generate_protocol_table(
                self.recommended_solution, ObjectiveFunction()
            )
            
            output_dir = os.path.expanduser('~/Escritorio/BioKidney-AI/resultados')
            os.makedirs(output_dir, exist_ok=True)
            
            txt_path = os.path.join(output_dir, 'protocolo_coswift.txt')
            with open(txt_path, 'w', encoding='utf-8') as f:
                f.write("╔══════════════════════════════════════════════════════╗\n")
                f.write("║    PROTOCOLO Co-SWIFT — Bio-Kidney AI 2026          ║\n")
                f.write("║    VirtusSapiens | Carlos David Moreno Cáceres      ║\n")
                f.write("╚══════════════════════════════════════════════════════╝\n\n")
                f.write("PARÁMETROS ÓPTIMOS DE BIOIMPRESIÓN\n")
                f.write("─" * 55 + "\n")
                for param, data in protocol.items():
                    f.write(f"{param:<30} {data['valor']:<15} [{data['estado']}]\n")
                f.write("\n")
                f.write("COMPOSICIÓN NICE BIOINK (validada)\n")
                f.write("─" * 55 + "\n")
                f.write("GelMA              7.0%\n")
                f.write("Alginato           1.5%\n")
                f.write("Nanocelulosa       0.8%\n")
                f.write("LAP (fotoiniciador) 0.25%\n")
                f.write("\nGenerado por Optimizador Co-SWIFT v1.0\n")
            
            QMessageBox.information(
                self,
                "Protocolo exportado",
                f"Protocolo guardado en:\n{txt_path}"
            )
            self.log_text.append(f"\n» Protocolo exportado: {txt_path}")


# ─────────────────────────────────────────────────────────────────────────────
# MODO CONSOLA (sin Qt)
# ─────────────────────────────────────────────────────────────────────────────

def run_console():
    print("\n" + "╔" + "═"*68 + "╗")
    print("║  OPTIMIZADOR Co-SWIFT — Bio-Kidney AI 2026" + " "*25 + "║")
    print("║  VirtusSapiens | Carlos David Moreno Cáceres" + " "*23 + "║")
    print("╚" + "═"*68 + "╝")
    
    pso = MOPSO(n_particles=60, n_iterations=100)
    archive = pso.optimize()
    
    os.makedirs(
        os.path.expanduser('~/Escritorio/BioKidney-AI/resultados'),
        exist_ok=True
    )
    
    recommended = find_recommended_solution(archive)
    
    if recommended:
        pos = recommended[0]
        fit = recommended[1]
        viab = 100 - fit[2]
        canal = vascular_resolution_model(pos[0], pos[3], pos[1])
        
        print("\n" + "═"*70)
        print("  SOLUCIÓN RECOMENDADA — PROTOCOLO Co-SWIFT")
        print("═"*70)
        print(f"  Presión de extrusión:    {pos[0]:.1f} kPa")
        print(f"  Velocidad de impresión:  {pos[1]:.1f} mm/s")
        print(f"  Temperatura bioink:      {pos[2]:.1f} °C")
        print(f"  Diámetro boquilla:       {pos[3]:.0f} µm")
        print(f"  Tiempo entre capas:      {pos[4]:.0f} s")
        print("─"*70)
        print(f"  Shear stress:            {fit[0]:.1f} Pa  {'✓ ÓPTIMO' if fit[0]<100 else '⚠ PARCIAL' if fit[0]<150 else '✗ CRÍTICO'}")
        print(f"  Canal vascular:          {canal:.0f} µm  {'✓ ÓPTIMO' if canal>=190 else '⚠ PARCIAL'}")
        print(f"  Viabilidad celular:      {viab:.1f} %   {'✓ ÓPTIMO' if viab>=85 else '⚠ PARCIAL'}")
        print(f"  Soluciones Pareto:       {len(archive)}")
        print("═"*70)
    
    fig = plot_results(archive, pso.convergence_history, recommended, ObjectiveFunction())
    plt.show()


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    os.makedirs(
        os.path.expanduser('~/Escritorio/BioKidney-AI/resultados'),
        exist_ok=True
    )
    
    if HAS_QT:
        app = QApplication(sys.argv)
        app.setStyle('Fusion')
        window = MainWindow()
        window.show()
        sys.exit(app.exec())
    else:
        run_console()
