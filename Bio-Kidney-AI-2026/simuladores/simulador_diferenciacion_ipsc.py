"""
╔══════════════════════════════════════════════════════════════════════════════╗
║       Bio-Kidney AI 2026 — Simulador de Diferenciación iPSC                ║
║       Módulo: Conversión iPSC → Linajes Renales Específicos                ║
║       Referencia: Takasato et al. 2015 | Freedman et al. 2015              ║
║       Autor: Carlos David Moreno Cáceres — VirtusSapiens                   ║
╚══════════════════════════════════════════════════════════════════════════════╝

Modela la diferenciación iPSC → Mesodermo Intermedio → NPC → {Podocitos,
Células Tubulares Proximales, Células Endoteliales}, calculando:
  • Concentración temporal de cada linaje
  • Ventana óptima de bioimpresión por linaje
  • Riesgo de teratoma (células indiferenciadas residuales)
  • Pureza fenotípica (objetivo ≥95%)
  • Dosis óptima de factores de señalización
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
import warnings
warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────
#  PALETA & ESTILO — coherente con el proyecto BioKidney-AI
# ─────────────────────────────────────────────────────────────
BG         = "#0D1117"
BG2        = "#161B22"
PANEL      = "#1C2128"
BORDER     = "#30363D"
TEXT_MAIN  = "#E6EDF3"
TEXT_DIM   = "#8B949E"
ACCENT     = "#58A6FF"   # azul GitHub-dark
GREEN      = "#3FB950"
ORANGE     = "#D29922"
RED        = "#F85149"
PURPLE     = "#BC8CFF"
TEAL       = "#39D353"
YELLOW     = "#E3B341"

# Colores por linaje
C_IPSC     = "#58A6FF"   # iPSC pluripotente
C_IM       = "#BC8CFF"   # Mesodermo Intermedio
C_NPC      = "#3FB950"   # Progenitores de Nefrona
C_PODO     = "#F0883E"   # Podocitos
C_PCT      = "#39D353"   # Células Tubulares Proximales
C_IEC      = "#79C0FF"   # Células Endoteliales
C_TERA     = "#F85149"   # Riesgo teratoma (indiferenciadas)

matplotlib.rcParams.update({
    "figure.facecolor":  BG,
    "axes.facecolor":    PANEL,
    "axes.edgecolor":    BORDER,
    "axes.labelcolor":   TEXT_MAIN,
    "xtick.color":       TEXT_DIM,
    "ytick.color":       TEXT_DIM,
    "text.color":        TEXT_MAIN,
    "legend.facecolor":  BG2,
    "legend.edgecolor":  BORDER,
    "grid.color":        BORDER,
    "grid.linestyle":    "--",
    "grid.alpha":        0.5,
    "font.family":       "monospace",
    "font.size":         9,
})

# ─────────────────────────────────────────────────────────────
#  PARÁMETROS BIOLÓGICOS (literatura + ajuste in-silico)
# ─────────────────────────────────────────────────────────────
# Dosis de factores de señalización (unidades normalizadas 0-1)
FACTORS = {
    "Wnt_CHIR99021":     {"dose_ng_mL": 8.0,  "window": (0,  5),  "Kd": 2.0},
    "BMP7":              {"dose_ng_mL": 50.0, "window": (5,  12), "Kd": 25.0},
    "AcidoRetinoico":    {"dose_ng_mL": 0.1,  "window": (5,  12), "Kd": 0.05},
    "GDNF":              {"dose_ng_mL": 100.0,"window": (12, 21), "Kd": 50.0},
    "FGF9":              {"dose_ng_mL": 200.0,"window": (12, 21), "Kd": 100.0},
}

# Constantes cinéticas del modelo ODE
k = {
    # Tasas de diferenciación (día⁻¹)
    "k_iPSC_to_IM":   0.55,   # iPSC → Mesodermo Intermedio
    "k_IM_to_NPC":    0.40,   # IM → NPC
    "k_NPC_to_PODO":  0.18,   # NPC → Podocitos
    "k_NPC_to_PCT":   0.22,   # NPC → Células Tubulares Prox.
    "k_NPC_to_IEC":   0.14,   # NPC → Células Endoteliales
    # Proliferación
    "p_iPSC":         0.30,   # tasa basal iPSC
    "p_NPC":          0.25,   # proliferación NPC
    "p_PODO":         0.05,
    "p_PCT":          0.08,
    "p_IEC":          0.10,
    # Apoptosis / muerte
    "d_iPSC":         0.10,
    "d_IM":           0.08,
    "d_NPC":          0.06,
    "d_PODO":         0.04,
    "d_PCT":          0.04,
    "d_IEC":          0.04,
    # Capacidad de carga (fracción del total)
    "K_total":        1.0,
}

# Densidad objetivo para bioimpresión: 200M células/mL (normalizada a 1.0)
DENSITY_TARGET = 200e6   # células/mL
N0             = 1e6     # inóculo inicial iPSC (células/mL)
PURITY_MIN     = 0.95    # 95%

# ─────────────────────────────────────────────────────────────
#  FUNCIONES DE SEÑALIZACIÓN (Hill cinetics)
# ─────────────────────────────────────────────────────────────
def hill(C, Kd, n=2):
    """Función de Hill saturante."""
    return (C**n) / (Kd**n + C**n)

def factor_activity(t, factor_name, dose_override=None):
    """Actividad de un factor en el tiempo t (0-1)."""
    f = FACTORS[factor_name]
    dose = dose_override if dose_override else f["dose_ng_mL"]
    t0, t1 = f["window"]
    # Rampa de entrada/salida suave
    if t < t0:
        return 0.0
    elif t > t1:
        return max(0.0, 1.0 - 0.3*(t - t1))
    else:
        ramp_in  = min(1.0, (t - t0) / 0.5)
        return ramp_in * hill(dose, f["Kd"])

# ─────────────────────────────────────────────────────────────
#  SISTEMA DE ODEs
# ─────────────────────────────────────────────────────────────
def odes(t, y, doses=None):
    """
    y = [iPSC, IM, NPC, PODO, PCT, IEC]  (densidades relativas, suma→1 al final)
    """
    ipsc, im, npc, podo, pct, iec = y
    total = max(ipsc + im + npc + podo + pct + iec, 1e-12)

    # Presión logística global
    logistic = max(0, 1.0 - total / k["K_total"])

    # Actividad de factores
    doses = doses or {}
    wnt  = factor_activity(t, "Wnt_CHIR99021",  doses.get("Wnt_CHIR99021"))
    bmp  = factor_activity(t, "BMP7",           doses.get("BMP7"))
    ra   = factor_activity(t, "AcidoRetinoico", doses.get("AcidoRetinoico"))
    gdnf = factor_activity(t, "GDNF",           doses.get("GDNF"))
    fgf  = factor_activity(t, "FGF9",           doses.get("FGF9"))

    # Tasas efectivas
    diff_iPSC_IM = k["k_iPSC_to_IM"]  * wnt  * ipsc
    diff_IM_NPC  = k["k_IM_to_NPC"]   * (bmp + ra) * 0.5 * im
    diff_NPC_PODO= k["k_NPC_to_PODO"] * gdnf * npc
    diff_NPC_PCT = k["k_NPC_to_PCT"]  * (gdnf + fgf) * 0.5 * npc
    diff_NPC_IEC = k["k_NPC_to_IEC"]  * fgf  * npc

    d_ipsc = (k["p_iPSC"] * logistic - k["d_iPSC"]) * ipsc - diff_iPSC_IM
    d_im   = diff_iPSC_IM - diff_IM_NPC - k["d_IM"] * im
    d_npc  = diff_IM_NPC  + k["p_NPC"] * logistic * npc \
             - diff_NPC_PODO - diff_NPC_PCT - diff_NPC_IEC - k["d_NPC"] * npc
    d_podo = diff_NPC_PODO + k["p_PODO"] * logistic * podo - k["d_PODO"] * podo
    d_pct  = diff_NPC_PCT  + k["p_PCT"]  * logistic * pct  - k["d_PCT"]  * pct
    d_iec  = diff_NPC_IEC  + k["p_IEC"]  * logistic * iec  - k["d_IEC"]  * iec

    return [d_ipsc, d_im, d_npc, d_podo, d_pct, d_iec]

# ─────────────────────────────────────────────────────────────
#  SIMULACIÓN PRINCIPAL
# ─────────────────────────────────────────────────────────────
def run_simulation(t_end=30, n_points=3000, doses=None):
    t_span = (0, t_end)
    t_eval = np.linspace(0, t_end, n_points)
    # Condiciones iniciales: 100% iPSC
    y0 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sol = solve_ivp(
        lambda t, y: odes(t, y, doses),
        t_span, y0, t_eval=t_eval,
        method="RK45", rtol=1e-8, atol=1e-10, max_step=0.05
    )
    return sol.t, sol.y  # y shape: (6, n_points)

# ─────────────────────────────────────────────────────────────
#  ANÁLISIS: ventanas, pureza, riesgo teratoma
# ─────────────────────────────────────────────────────────────
def compute_metrics(t, y):
    ipsc, im, npc, podo, pct, iec = y
    total = ipsc + im + npc + podo + pct + iec
    total = np.where(total < 1e-12, 1e-12, total)

    fracs = {
        "iPSC":  ipsc / total,
        "IM":    im   / total,
        "NPC":   npc  / total,
        "PODO":  podo / total,
        "PCT":   pct  / total,
        "IEC":   iec  / total,
    }

    # Teratoma ≈ fracción de iPSC residual (log-riesgo)
    teratoma_risk = fracs["iPSC"]  # [0,1]

    # Ventana óptima: fracción ≥ 0.95 de un linaje terminal
    windows = {}
    for linaje in ["PODO", "PCT", "IEC"]:
        mask = fracs[linaje] >= PURITY_MIN
        idx  = np.where(mask)[0]
        if len(idx) > 0:
            windows[linaje] = (round(t[idx[0]], 1), round(t[idx[-1]], 1))
        else:
            # encontrar máximo aunque no llegue al 95%
            peak_idx = np.argmax(fracs[linaje])
            windows[linaje] = (None, round(t[peak_idx], 1))

    # Índice de pureza fenotípica (área bajo curva >95% / área total)
    purity_indices = {}
    for linaje in ["PODO", "PCT", "IEC"]:
        f = fracs[linaje]
        area_above = np.trapezoid(np.where(f >= PURITY_MIN, f, 0), t)
        area_total = np.trapezoid(f, t)
        purity_indices[linaje] = area_above / (area_total + 1e-12)

    return fracs, teratoma_risk, windows, purity_indices

# ─────────────────────────────────────────────────────────────
#  OPTIMIZACIÓN DE DOSIS
# ─────────────────────────────────────────────────────────────
def optimize_dose(factor_name, linaje_target, t_end=30):
    """Encuentra la dosis que maximiza la fracción pico del linaje objetivo."""
    linaje_idx = {"PODO": 3, "PCT": 4, "IEC": 5}[linaje_target]
    f_ref = FACTORS[factor_name]["dose_ng_mL"]

    def neg_peak(dose_scale):
        doses_test = {factor_name: f_ref * dose_scale}
        _, y = run_simulation(t_end=t_end, n_points=500, doses=doses_test)
        total = y.sum(axis=0)
        total = np.where(total < 1e-12, 1e-12, total)
        frac = y[linaje_idx] / total
        return -np.max(frac)

    result = minimize_scalar(neg_peak, bounds=(0.1, 5.0), method="bounded")
    return f_ref * result.x, -result.fun  # (dosis_óptima, pureza_máx)

# ─────────────────────────────────────────────────────────────
#  VISUALIZACIÓN
# ─────────────────────────────────────────────────────────────
def plot_all(t, y, fracs, teratoma_risk, windows, purity_indices):
    fig = plt.figure(figsize=(20, 14), facecolor=BG)
    fig.suptitle(
        "Bio-Kidney AI 2026  ·  Simulador de Diferenciación iPSC → Linajes Renales",
        fontsize=15, color=ACCENT, fontweight="bold", y=0.98,
        fontfamily="monospace"
    )

    gs = gridspec.GridSpec(
        3, 3, figure=fig,
        hspace=0.52, wspace=0.38,
        top=0.93, bottom=0.07, left=0.06, right=0.97
    )

    ipsc_arr, im_arr, npc_arr, podo_arr, pct_arr, iec_arr = y
    total = y.sum(axis=0)
    total = np.where(total < 1e-12, 1e-12, total)

    # ── Panel 1: Densidades absolutas (escala relativa) ──────────────────
    ax1 = fig.add_subplot(gs[0, :2])
    ax1.set_facecolor(PANEL)
    ax1.set_title("Densidad celular relativa por linaje (modelo ODE)", color=TEXT_MAIN, fontsize=10)
    ax1.plot(t, ipsc_arr, color=C_IPSC,  lw=1.8, label="iPSC")
    ax1.plot(t, im_arr,   color=C_IM,    lw=1.8, label="Mesodermo Int.")
    ax1.plot(t, npc_arr,  color=C_NPC,   lw=1.8, label="NPC")
    ax1.plot(t, podo_arr, color=C_PODO,  lw=2.2, label="Podocitos", ls="--")
    ax1.plot(t, pct_arr,  color=C_PCT,   lw=2.2, label="Tub. Prox.", ls="--")
    ax1.plot(t, iec_arr,  color=C_IEC,   lw=2.2, label="iECs", ls="--")
    # zonas de señalización
    for span, col, lbl in [((0,5),"#58A6FF20","Wnt"),((5,12),"#BC8CFF20","BMP7+RA"),((12,21),"#3FB95020","GDNF+FGF")]:
        ax1.axvspan(*span, color=col, alpha=0.25)
        ax1.text((span[0]+span[1])/2, ax1.get_ylim()[1]*0.02, lbl,
                 color=TEXT_DIM, ha="center", fontsize=7.5)
    ax1.set_xlabel("Tiempo (días)", color=TEXT_DIM)
    ax1.set_ylabel("Densidad relativa", color=TEXT_DIM)
    ax1.legend(loc="upper right", fontsize=8)
    ax1.grid(True)
    ax1.xaxis.set_minor_locator(MultipleLocator(1))

    # ── Panel 2: Fracción fenotípica (pureza) ────────────────────────────
    ax2 = fig.add_subplot(gs[0, 2])
    ax2.set_facecolor(PANEL)
    ax2.set_title("Pureza fenotípica (%)", color=TEXT_MAIN, fontsize=10)
    for linaje, col, arr in [
        ("PODO", C_PODO, fracs["PODO"]),
        ("PCT",  C_PCT,  fracs["PCT"]),
        ("IEC",  C_IEC,  fracs["IEC"]),
        ("iPSC", C_TERA, fracs["iPSC"]),
    ]:
        ax2.plot(t, arr*100, color=col, lw=1.8, label=linaje)
    ax2.axhline(95, color=ORANGE, lw=1.2, ls=":", label="95% objetivo")
    ax2.set_xlabel("Tiempo (días)", color=TEXT_DIM)
    ax2.set_ylabel("Fracción (%)", color=TEXT_DIM)
    ax2.set_ylim(0, 105)
    ax2.legend(fontsize=8)
    ax2.grid(True)

    # ── Panel 3: Riesgo de teratoma ───────────────────────────────────────
    ax3 = fig.add_subplot(gs[1, :2])
    ax3.set_facecolor(PANEL)
    ax3.set_title("Riesgo de Teratoma — iPSC indiferenciadas residuales", color=TEXT_MAIN, fontsize=10)
    risk_pct = teratoma_risk * 100
    # gradiente de color por nivel de riesgo
    for i in range(len(t)-1):
        r = risk_pct[i]
        col = RED if r > 5 else ORANGE if r > 1 else GREEN
        ax3.fill_between(t[i:i+2], risk_pct[i:i+2], color=col, alpha=0.7)
    ax3.plot(t, risk_pct, color=RED, lw=0.8, alpha=0.4)
    ax3.axhline(5.0, color=RED,    lw=1.2, ls="--", label="Riesgo alto (>5%)")
    ax3.axhline(1.0, color=ORANGE, lw=1.2, ls="--", label="Riesgo moderado (>1%)")
    # anotación ventanas
    for linaje, (t0, t1) in windows.items():
        if t0 is not None:
            ax3.axvspan(t0, t1, alpha=0.1, color=GREEN)
            ax3.text((t0+t1)/2, ax3.get_ylim()[1]*0.8 if i==0 else ax3.get_ylim()[1]*0.6,
                     f"✓{linaje}", color=GREEN, ha="center", fontsize=8)
    ax3.set_xlabel("Tiempo (días)", color=TEXT_DIM)
    ax3.set_ylabel("iPSC residual (%)", color=TEXT_DIM)
    ax3.legend(fontsize=8)
    ax3.grid(True)
    ax3.set_yscale("log")

    # ── Panel 4: Actividad de factores de señalización ───────────────────
    ax4 = fig.add_subplot(gs[1, 2])
    ax4.set_facecolor(PANEL)
    ax4.set_title("Actividad factores de señalización", color=TEXT_MAIN, fontsize=10)
    t_dense = np.linspace(0, 30, 1000)
    factor_styles = {
        "Wnt_CHIR99021":  (C_IPSC,  "-",  "Wnt/CHIR"),
        "BMP7":           (C_IM,    "-",  "BMP7"),
        "AcidoRetinoico": (PURPLE,  "--", "Ac. Retin."),
        "GDNF":           (C_NPC,   "-",  "GDNF"),
        "FGF9":           (C_IEC,   "--", "FGF9"),
    }
    for fname, (col, ls, lbl) in factor_styles.items():
        act = [factor_activity(ti, fname) for ti in t_dense]
        ax4.plot(t_dense, act, color=col, lw=1.8, ls=ls, label=lbl)
    ax4.set_xlabel("Tiempo (días)", color=TEXT_DIM)
    ax4.set_ylabel("Actividad (norm.)", color=TEXT_DIM)
    ax4.legend(fontsize=7.5)
    ax4.grid(True)

    # ── Panel 5: Mapa de calor de pureza fenotípica ───────────────────────
    ax5 = fig.add_subplot(gs[2, :2])
    ax5.set_facecolor(PANEL)
    ax5.set_title("Fracción fenotípica — mapa temporal (linajes terminales)", color=TEXT_MAIN, fontsize=10)
    heat_data = np.vstack([fracs["PODO"], fracs["PCT"], fracs["IEC"], fracs["iPSC"]])
    im_plot = ax5.imshow(
        heat_data, aspect="auto",
        extent=[t[0], t[-1], -0.5, 3.5],
        cmap="plasma", vmin=0, vmax=1, origin="lower"
    )
    ax5.set_yticks([0, 1, 2, 3])
    ax5.set_yticklabels(["iPSC (riesgo)", "iECs", "Tub. Prox.", "Podocitos"], color=TEXT_MAIN)
    ax5.set_xlabel("Tiempo (días)", color=TEXT_DIM)
    cbar = fig.colorbar(im_plot, ax=ax5, fraction=0.02, pad=0.01)
    cbar.set_label("Fracción fenotípica", color=TEXT_DIM)
    cbar.ax.yaxis.set_tick_params(color=TEXT_DIM)
    # línea del 95%
    for row_idx, linaje in enumerate(["PODO", "PCT", "IEC"]):
        f = fracs[linaje]
        idx95 = np.where(f >= PURITY_MIN)[0]
        if len(idx95):
            ax5.axvline(t[idx95[0]], color=GREEN, lw=0.8, ls=":", alpha=0.7)

    # ── Panel 6: Tabla de resumen ─────────────────────────────────────────
    ax6 = fig.add_subplot(gs[2, 2])
    ax6.set_facecolor(BG2)
    ax6.axis("off")
    ax6.set_title("Resumen de Ventanas Óptimas", color=TEXT_MAIN, fontsize=10)

    linaje_info = [
        ("Podocitos",   "PODO", C_PODO,  "Filtración glomerular"),
        ("Tub. Prox.",  "PCT",  C_PCT,   "Reabsorción tubular"),
        ("iECs",        "IEC",  C_IEC,   "Vascularización"),
    ]

    y_pos = 0.92
    for full, key, col, func in linaje_info:
        w = windows.get(key, (None, None))
        pi = purity_indices.get(key, 0)
        # ícono de estado
        if w[0] is not None:
            estado = f"✓  Días {w[0]}–{w[1]}"
            ec = GREEN
        else:
            estado = f"⚠  Pico día {w[1]}"
            ec = ORANGE
        pi_pct = pi * 100

        ax6.text(0.02, y_pos,       f"▌ {full}", transform=ax6.transAxes,
                 color=col, fontsize=9, fontweight="bold")
        ax6.text(0.02, y_pos-0.055, f"   {func}", transform=ax6.transAxes,
                 color=TEXT_DIM, fontsize=7.5)
        ax6.text(0.02, y_pos-0.11,  f"   {estado}", transform=ax6.transAxes,
                 color=ec, fontsize=8)
        ax6.text(0.02, y_pos-0.165, f"   Índice pureza: {pi_pct:.1f}%",
                 transform=ax6.transAxes, color=TEXT_DIM, fontsize=7.5)
        ax6.plot([0.02, 0.98], [y_pos-0.20, y_pos-0.20],
                    color=BORDER, lw=0.5, transform=ax6.transAxes)
        y_pos -= 0.28

    # Riesgo teratoma al día 21
    idx21 = np.argmin(np.abs(t - 21))
    risk21 = teratoma_risk[idx21] * 100
    col_risk = RED if risk21 > 5 else ORANGE if risk21 > 1 else GREEN
    ax6.text(0.02, y_pos, f"⚠ Teratoma día 21:", transform=ax6.transAxes,
             color=TEXT_DIM, fontsize=8)
    ax6.text(0.50, y_pos, f"{risk21:.2f}%", transform=ax6.transAxes,
             color=col_risk, fontsize=9, fontweight="bold")

    # ── Footer ────────────────────────────────────────────────────────────
    fig.text(0.5, 0.01,
             "VirtusSapiens · Bio-Kidney AI 2026  |  Ref: Takasato 2015, Freedman 2015  |"
             "  Módulo: simulador_diferenciacion_ipsc.py",
             ha="center", color=TEXT_DIM, fontsize=7.5)

    return fig

# ─────────────────────────────────────────────────────────────
#  REPORTE TEXTUAL EN CONSOLA
# ─────────────────────────────────────────────────────────────
def print_report(t, fracs, teratoma_risk, windows, purity_indices):
    sep = "─" * 70
    print(f"\n{'═'*70}")
    print("  Bio-Kidney AI 2026 — Reporte Diferenciación iPSC")
    print(f"{'═'*70}")
    print(f"\n{'VENTANAS ÓPTIMAS DE BIOIMPRESIÓN':^70}")
    print(sep)
    print(f"  {'Linaje':<22} {'Ventana (días)':<22} {'Índice Pureza':>14}")
    print(sep)
    linaje_names = {"PODO": "Podocitos (filtración)",
                    "PCT":  "Tub. Proximales",
                    "IEC":  "Endoteliales (iECs)"}
    for key, name in linaje_names.items():
        w  = windows.get(key, (None, None))
        pi = purity_indices.get(key, 0) * 100
        if w[0] is not None:
            w_str = f"Días {w[0]:.1f} – {w[1]:.1f}"
            estado = "✓ APTO"
        else:
            w_str = f"Pico en día {w[1]:.1f} (sub-95%)"
            estado = "⚠ AJUSTAR"
        print(f"  {name:<22} {w_str:<22} {pi:>9.1f}%  {estado}")
    print(sep)

    print(f"\n{'RIESGO DE TERATOMA':^70}")
    print(sep)
    for dia in [7, 14, 21, 28, 30]:
        idx = np.argmin(np.abs(t - dia))
        r = teratoma_risk[idx] * 100
        nivel = "ALTO ⛔" if r > 5 else "MODERADO ⚠" if r > 1 else "BAJO ✓"
        print(f"  Día {dia:>2}:  iPSC residual = {r:7.3f}%   →  {nivel}")
    print(sep)

    print(f"\n{'DOSIS ÓPTIMAS DE FACTORES (ESTIMACIÓN IN SILICO)':^70}")
    print(sep)
    print(f"  {'Factor':<22} {'Dosis Ref':<16} {'Ventana':<16} Objetivo")
    print(sep)
    factor_targets = {
        "Wnt_CHIR99021":  "IM", "BMP7": "NPC", "AcidoRetinoico": "NPC",
        "GDNF": "PODO", "FGF9": "IEC"
    }
    for fname, target in factor_targets.items():
        f    = FACTORS[fname]
        dose = f["dose_ng_mL"]
        win  = f["window"]
        unit = "ng/mL"
        print(f"  {fname:<22} {dose:<10.1f} {unit:<6} Días {win[0]}-{win[1]:<8} → {target}")
    print(sep)

    print(f"\n{'ESTADO GLOBAL DEL MÓDULO':^70}")
    print(sep)
    apto = sum(1 for w in windows.values() if w[0] is not None)
    print(f"  Linajes con pureza ≥95%:      {apto}/3")
    idx21 = np.argmin(np.abs(t - 21))
    r21 = teratoma_risk[idx21] * 100
    print(f"  Riesgo teratoma día 21:        {r21:.3f}%")
    print(f"  Densidad objetivo:             {DENSITY_TARGET/1e6:.0f}M células/mL")
    print(f"  Protocolo base:                Takasato 2015 + Freedman 2015")
    print(f"{'═'*70}\n")

# ─────────────────────────────────────────────────────────────
#  MAIN
# ─────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("\n[Bio-Kidney AI 2026] Iniciando simulación de diferenciación iPSC...")
    t, y = run_simulation(t_end=30, n_points=3000)
    fracs, teratoma_risk, windows, purity_indices = compute_metrics(t, y)
    print("[Bio-Kidney AI 2026] Simulación completada. Generando análisis...")

    print_report(t, fracs, teratoma_risk, windows, purity_indices)

    fig = plot_all(t, y, fracs, teratoma_risk, windows, purity_indices)

    # ── Guardar resultado ────────────────────────────────────────────────
    output_path = "simulador_diferenciacion_ipsc_output.png"
    fig.savefig(output_path, dpi=150, bbox_inches="tight", facecolor=BG)
    print(f"[Bio-Kidney AI 2026] Figura guardada → {output_path}")
    plt.tight_layout(rect=[0, 0.02, 1, 0.97])
    plt.show()
