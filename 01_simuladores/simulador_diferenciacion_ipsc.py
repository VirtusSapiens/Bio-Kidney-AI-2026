"""
╔══════════════════════════════════════════════════════════════════════════════╗
║       Bio-Kidney AI 2026 — Simulador de Diferenciación iPSC                ║
║       Módulo: Conversión iPSC → Linajes Renales Específicos                ║
║       Referencia: Takasato et al. 2015 | Freedman et al. 2015              ║
║       Autor: Carlos David Moreno Cáceres — VirtusSapiens                   ║
╚══════════════════════════════════════════════════════════════════════════════╝

Arquitectura: 3 protocolos de cultivo independientes (como en laboratorio real)
  • Protocolo A → Podocitos        (GDNF + BMP4 dominante)
  • Protocolo B → Tub. Proximales  (FGF9 + EGF)
  • Protocolo C → iECs             (FGF9 alto + VEGFA)

Cada protocolo comparte la fase común iPSC→IM→NPC (días 0-12),
luego diverge con sus factores específicos (días 12-30).
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
from scipy.integrate import solve_ivp
import warnings
warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────
#  PALETA & ESTILO
# ─────────────────────────────────────────────────────────────
BG        = "#0D1117"
BG2       = "#161B22"
PANEL     = "#1C2128"
BORDER    = "#30363D"
TEXT_MAIN = "#E6EDF3"
TEXT_DIM  = "#8B949E"
ACCENT    = "#58A6FF"
GREEN     = "#3FB950"
ORANGE    = "#D29922"
RED       = "#F85149"
PURPLE    = "#BC8CFF"

C_IPSC  = "#58A6FF"
C_IM    = "#BC8CFF"
C_NPC   = "#3FB950"
C_PODO  = "#F0883E"
C_PCT   = "#39D353"
C_IEC   = "#79C0FF"

matplotlib.rcParams.update({
    "figure.facecolor": BG,    "axes.facecolor":   PANEL,
    "axes.edgecolor":   BORDER,"axes.labelcolor":  TEXT_MAIN,
    "xtick.color":      TEXT_DIM,"ytick.color":    TEXT_DIM,
    "text.color":       TEXT_MAIN,"legend.facecolor": BG2,
    "legend.edgecolor": BORDER,  "grid.color":     BORDER,
    "grid.linestyle":   "--",    "grid.alpha":     0.5,
    "font.family":      "monospace","font.size":   9,
})

# ─────────────────────────────────────────────────────────────
#  PARÁMETROS
# ─────────────────────────────────────────────────────────────
PURITY_MIN     = 0.95
DENSITY_TARGET = 200e6

FACTORS_COMMON = {
    "Wnt_CHIR99021":  {"dose_ng_mL": 8.0,  "window": (0,  5),  "Kd": 2.0},
    "BMP7":           {"dose_ng_mL": 50.0, "window": (5,  12), "Kd": 25.0},
    "AcidoRetinoico": {"dose_ng_mL": 0.1,  "window": (5,  12), "Kd": 0.05},
}

PROTOCOLS = {
    "PODO": {
        "label":   "Podocitos",
        "funcion": "Filtración glomerular",
        "color":   C_PODO,
        "factors": {
            "GDNF": {"dose_ng_mL": 100.0, "window": (12, 28), "Kd": 30.0},
            "BMP4": {"dose_ng_mL": 10.0,  "window": (14, 28), "Kd": 5.0},
        },
        "k_NPC_to_target": 1.40,
        "p_target": 0.06,
        "d_target": 0.02,
    },
    "PCT": {
        "label":   "Tub. Proximales",
        "funcion": "Reabsorción tubular",
        "color":   C_PCT,
        "factors": {
            "FGF9": {"dose_ng_mL": 200.0, "window": (12, 26), "Kd": 80.0},
            "EGF":  {"dose_ng_mL": 10.0,  "window": (15, 26), "Kd": 4.0},
        },
        "k_NPC_to_target": 1.20,
        "p_target": 0.08,
        "d_target": 0.02,
    },
    "IEC": {
        "label":   "Endoteliales (iECs)",
        "funcion": "Vascularización",
        "color":   C_IEC,
        "factors": {
            "FGF9":  {"dose_ng_mL": 300.0, "window": (12, 28), "Kd": 100.0},
            "VEGFA": {"dose_ng_mL": 50.0,  "window": (14, 28), "Kd": 20.0},
        },
        "k_NPC_to_target": 1.00,
        "p_target": 0.10,
        "d_target": 0.02,
    },
}

K_COMMON = {
    "k_iPSC_to_IM": 1.20,
    "k_IM_to_NPC":  0.90,
    "p_iPSC":       0.08,
    "p_NPC":        0.18,
    "d_iPSC":       0.60,
    "d_IM":         0.10,
    "d_NPC":        0.04,
    "K_total":      1.0,
}

# ─────────────────────────────────────────────────────────────
#  SEÑALIZACIÓN
# ─────────────────────────────────────────────────────────────
def hill(C, Kd, n=2):
    return (C**n) / (Kd**n + C**n)

def factor_activity(t, fdict):
    dose = fdict["dose_ng_mL"]
    t0, t1 = fdict["window"]
    if t < t0:
        return 0.0
    elif t > t1:
        return max(0.0, 1.0 - 0.4*(t - t1))
    else:
        ramp = min(1.0, (t - t0) / 0.8)
        return ramp * hill(dose, fdict["Kd"])

def common_signal(t):
    wnt = factor_activity(t, FACTORS_COMMON["Wnt_CHIR99021"])
    bmp = factor_activity(t, FACTORS_COMMON["BMP7"])
    ra  = factor_activity(t, FACTORS_COMMON["AcidoRetinoico"])
    return wnt, bmp, ra

# ─────────────────────────────────────────────────────────────
#  ODE POR PROTOCOLO: y = [iPSC, IM, NPC, TARGET]
# ─────────────────────────────────────────────────────────────
def ode_protocol(t, y, proto_key):
    ipsc, im, npc, target = y
    proto    = PROTOCOLS[proto_key]
    kc       = K_COMMON
    total    = max(ipsc + im + npc + target, 1e-12)
    logistic = max(0.0, 1.0 - total / kc["K_total"])

    wnt, bmp, ra = common_signal(t)
    spec_signals = [factor_activity(t, fd) for fd in proto["factors"].values()]
    spec = sum(spec_signals) / len(spec_signals)

    diff_iPSC_IM    = kc["k_iPSC_to_IM"] * wnt * ipsc
    diff_IM_NPC     = kc["k_IM_to_NPC"]  * (bmp + ra) * 0.5 * im
    diff_NPC_target = proto["k_NPC_to_target"] * spec * npc

    d_ipsc  = (kc["p_iPSC"] * logistic - kc["d_iPSC"]) * ipsc - diff_iPSC_IM
    d_im    = diff_iPSC_IM - diff_IM_NPC - kc["d_IM"] * im
    d_npc   = diff_IM_NPC + kc["p_NPC"] * logistic * npc \
              - diff_NPC_target - kc["d_NPC"] * npc
    d_target = diff_NPC_target \
               + proto["p_target"] * logistic * target \
               - proto["d_target"] * target

    return [d_ipsc, d_im, d_npc, d_target]

# ─────────────────────────────────────────────────────────────
#  SIMULACIÓN
# ─────────────────────────────────────────────────────────────
def run_all_protocols(t_end=30, n_points=3000):
    t_eval = np.linspace(0, t_end, n_points)
    y0 = [1.0, 0.0, 0.0, 0.0]
    results = {}
    for pk in PROTOCOLS:
        sol = solve_ivp(
            lambda t, y, pk=pk: ode_protocol(t, y, pk),
            (0, t_end), y0, t_eval=t_eval,
            method="RK45", rtol=1e-8, atol=1e-10, max_step=0.05
        )
        results[pk] = {"t": sol.t, "y": sol.y}
    return results

# ─────────────────────────────────────────────────────────────
#  MÉTRICAS
# ─────────────────────────────────────────────────────────────
def compute_metrics(results):
    metrics = {}
    for pk, res in results.items():
        t = res["t"]
        ipsc, im, npc, target = res["y"]
        total    = ipsc + im + npc + target
        total    = np.where(total < 1e-12, 1e-12, total)
        purity   = target / total
        teratoma = ipsc   / total

        mask = purity >= PURITY_MIN
        idx  = np.where(mask)[0]
        if len(idx):
            window = (round(t[idx[0]], 1), round(t[idx[-1]], 1))
        else:
            peak_idx = np.argmax(purity)
            window   = (None, round(t[peak_idx], 1))

        area_above   = np.trapezoid(np.where(purity >= PURITY_MIN, purity, 0), t)
        area_total   = np.trapezoid(purity, t)
        purity_index = area_above / (area_total + 1e-12)
        peak_purity  = float(np.max(purity))

        metrics[pk] = {
            "t": t, "purity": purity, "teratoma": teratoma,
            "total": total, "window": window,
            "purity_index": purity_index, "peak_purity": peak_purity,
        }
    return metrics

# ─────────────────────────────────────────────────────────────
#  REPORTE CONSOLA
# ─────────────────────────────────────────────────────────────
def print_report(metrics):
    sep = "─" * 70
    print(f"\n{'═'*70}")
    print("  Bio-Kidney AI 2026 — Reporte Diferenciación iPSC")
    print(f"{'═'*70}")
    print(f"\n{'VENTANAS ÓPTIMAS DE BIOIMPRESIÓN':^70}")
    print(sep)
    print(f"  {'Protocolo':<24} {'Ventana (días)':<24} {'Pureza pico':>10}  Estado")
    print(sep)
    for pk, m in metrics.items():
        proto = PROTOCOLS[pk]
        w  = m["window"]
        pp = m["peak_purity"] * 100
        if w[0] is not None:
            w_str  = f"Días {w[0]:.1f} – {w[1]:.1f}"
            estado = "✓ APTO"
        else:
            w_str  = f"Pico día {w[1]:.1f} ({pp:.1f}%)"
            estado = "⚠ AJUSTAR"
        print(f"  {proto['label']:<24} {w_str:<24} {pp:>8.1f}%  {estado}")
    print(sep)

    print(f"\n{'RIESGO DE TERATOMA POR PROTOCOLO':^70}")
    print(sep)
    for pk, m in metrics.items():
        proto = PROTOCOLS[pk]
        print(f"\n  [{proto['label']}]")
        for dia in [7, 14, 21, 28, 30]:
            idx   = np.argmin(np.abs(m["t"] - dia))
            r     = m["teratoma"][idx] * 100
            nivel = "ALTO ⛔" if r > 5 else "MODERADO ⚠" if r > 1 else "BAJO ✓"
            print(f"    Día {dia:>2}: iPSC residual = {r:7.4f}%  → {nivel}")
    print(sep)

    aptos = sum(1 for m in metrics.values() if m["window"][0] is not None)
    print(f"\n{'ESTADO GLOBAL':^70}")
    print(sep)
    print(f"  Protocolos con pureza ≥95%:   {aptos}/3")
    print(f"  Densidad objetivo:             {DENSITY_TARGET/1e6:.0f}M células/mL")
    print(f"  Protocolo base:                Takasato 2015 + Freedman 2015")
    print(f"{'═'*70}\n")

# ─────────────────────────────────────────────────────────────
#  VISUALIZACIÓN
# ─────────────────────────────────────────────────────────────
def plot_all(results, metrics):
    fig = plt.figure(figsize=(20, 15), facecolor=BG)
    fig.suptitle(
        "Bio-Kidney AI 2026  ·  Simulador de Diferenciación iPSC → Linajes Renales",
        fontsize=15, color=ACCENT, fontweight="bold", y=0.98, fontfamily="monospace"
    )
    gs = gridspec.GridSpec(3, 3, figure=fig,
        hspace=0.52, wspace=0.38, top=0.93, bottom=0.07, left=0.06, right=0.97)

    proto_list = list(PROTOCOLS.keys())

    # ── Panel 1: Densidades protocolo PODO ──────────────────────────────
    ax1 = fig.add_subplot(gs[0, :2])
    ax1.set_facecolor(PANEL)
    ax1.set_title("Densidad celular relativa — Protocolo Podocitos (referencia)", color=TEXT_MAIN, fontsize=10)
    m_ref = metrics["PODO"]
    t_ref = m_ref["t"]
    y_ref = results["PODO"]["y"]
    total_ref = m_ref["total"]
    ax1.plot(t_ref, y_ref[0]/total_ref, color=C_IPSC, lw=1.8, label="iPSC")
    ax1.plot(t_ref, y_ref[1]/total_ref, color=C_IM,   lw=1.8, label="Mesodermo Int.")
    ax1.plot(t_ref, y_ref[2]/total_ref, color=C_NPC,  lw=1.8, label="NPC")
    ax1.plot(t_ref, m_ref["purity"],    color=C_PODO, lw=2.4, ls="--", label="Podocitos")
    ax1.plot(metrics["PCT"]["t"], metrics["PCT"]["purity"], color=C_PCT, lw=1.2, ls=":", alpha=0.6, label="PCT (ref)")
    ax1.plot(metrics["IEC"]["t"], metrics["IEC"]["purity"], color=C_IEC, lw=1.2, ls=":", alpha=0.6, label="IEC (ref)")
    for span, col, lbl in [((0,5),"#58A6FF20","Wnt"),((5,12),"#BC8CFF20","BMP7+RA"),((12,28),"#F0883E20","GDNF+BMP4")]:
        ax1.axvspan(*span, color=col, alpha=0.22)
        ax1.text((span[0]+span[1])/2, 0.97, lbl, color=TEXT_DIM,
                 ha="center", fontsize=7.5, transform=ax1.get_xaxis_transform(), va="top")
    ax1.set_xlabel("Tiempo (días)", color=TEXT_DIM)
    ax1.set_ylabel("Fracción total", color=TEXT_DIM)
    ax1.legend(loc="upper right", fontsize=7.5, ncol=2)
    ax1.grid(True)

    # ── Panel 2: Pureza 3 protocolos ────────────────────────────────────
    ax2 = fig.add_subplot(gs[0, 2])
    ax2.set_facecolor(PANEL)
    ax2.set_title("Pureza fenotípica por protocolo (%)", color=TEXT_MAIN, fontsize=10)
    for pk in proto_list:
        m   = metrics[pk]
        col = PROTOCOLS[pk]["color"]
        lbl = PROTOCOLS[pk]["label"]
        ax2.plot(m["t"], m["purity"]*100, color=col, lw=2.0, label=lbl)
        w = m["window"]
        if w[0] is not None:
            ax2.axvspan(w[0], w[1], color=col, alpha=0.10)
    ax2.axhline(95, color=ORANGE, lw=1.2, ls=":", label="95% objetivo")
    ax2.set_xlabel("Tiempo (días)", color=TEXT_DIM)
    ax2.set_ylabel("Pureza (%)", color=TEXT_DIM)
    ax2.set_ylim(0, 105)
    ax2.legend(fontsize=8)
    ax2.grid(True)

    # ── Panel 3: Riesgo teratoma ─────────────────────────────────────────
    ax3 = fig.add_subplot(gs[1, :2])
    ax3.set_facecolor(PANEL)
    ax3.set_title("Riesgo de Teratoma — iPSC residuales por protocolo", color=TEXT_MAIN, fontsize=10)
    for pk in proto_list:
        m        = metrics[pk]
        col      = PROTOCOLS[pk]["color"]
        lbl      = PROTOCOLS[pk]["label"]
        risk_pct = m["teratoma"] * 100
        ax3.plot(m["t"], risk_pct, color=col, lw=1.8, label=lbl)
        for mask, fc in [
            (risk_pct > 5.0,                       RED),
            ((risk_pct > 1.0) & (risk_pct <= 5.0), ORANGE),
            (risk_pct <= 1.0,                      GREEN),
        ]:
            ax3.fill_between(m["t"], risk_pct, where=mask, color=fc, alpha=0.10)
    ax3.axhline(5.0, color=RED,    lw=1.0, ls="--", label="Alto (>5%)")
    ax3.axhline(1.0, color=ORANGE, lw=1.0, ls="--", label="Mod. (>1%)")
    ax3.set_yscale("log")
    ax3.set_xlabel("Tiempo (días)", color=TEXT_DIM)
    ax3.set_ylabel("iPSC residual (%)", color=TEXT_DIM)
    ax3.legend(fontsize=8, ncol=2)
    ax3.grid(True)

    # ── Panel 4: Actividad factores ──────────────────────────────────────
    ax4 = fig.add_subplot(gs[1, 2])
    ax4.set_facecolor(PANEL)
    ax4.set_title("Actividad factores de señalización", color=TEXT_MAIN, fontsize=10)
    t_dense = np.linspace(0, 30, 1000)
    for fname, col, ls, lbl in [
        ("Wnt_CHIR99021",  C_IPSC, "-",  "Wnt/CHIR"),
        ("BMP7",           C_IM,   "-",  "BMP7"),
        ("AcidoRetinoico", PURPLE, "--", "Ac. Retin."),
    ]:
        act = [factor_activity(ti, FACTORS_COMMON[fname]) for ti in t_dense]
        ax4.plot(t_dense, act, color=col, lw=1.6, ls=ls, label=lbl)
    for pk, fname, col, ls, lbl in [
        ("PODO", "GDNF",  C_PODO, "--", "GDNF→Podo"),
        ("PCT",  "FGF9",  C_PCT,  "--", "FGF9→PCT"),
        ("IEC",  "VEGFA", C_IEC,  ":",  "VEGFA→iEC"),
    ]:
        act = [factor_activity(ti, PROTOCOLS[pk]["factors"][fname]) for ti in t_dense]
        ax4.plot(t_dense, act, color=col, lw=1.6, ls=ls, label=lbl)
    ax4.set_xlabel("Tiempo (días)", color=TEXT_DIM)
    ax4.set_ylabel("Actividad (norm.)", color=TEXT_DIM)
    ax4.legend(fontsize=7, ncol=2)
    ax4.grid(True)

    # ── Panel 5: Heatmap ─────────────────────────────────────────────────
    ax5 = fig.add_subplot(gs[2, :2])
    ax5.set_facecolor(PANEL)
    ax5.set_title("Mapa de pureza fenotípica — 3 protocolos independientes", color=TEXT_MAIN, fontsize=10)
    heat_data = np.vstack([
        metrics["PODO"]["purity"],
        metrics["PCT"]["purity"],
        metrics["IEC"]["purity"],
        metrics["PODO"]["teratoma"],
    ])
    im_plot = ax5.imshow(heat_data, aspect="auto",
        extent=[t_ref[0], t_ref[-1], -0.5, 3.5],
        cmap="plasma", vmin=0, vmax=1, origin="lower")
    ax5.set_yticks([0, 1, 2, 3])
    ax5.set_yticklabels(["iPSC riesgo", "iECs", "Tub. Prox.", "Podocitos"], color=TEXT_MAIN)
    ax5.set_xlabel("Tiempo (días)", color=TEXT_DIM)
    cbar = fig.colorbar(im_plot, ax=ax5, fraction=0.02, pad=0.01)
    cbar.set_label("Pureza fenotípica", color=TEXT_DIM)
    for pk in ["PODO", "PCT", "IEC"]:
        w = metrics[pk]["window"]
        if w[0] is not None:
            for xv in [w[0], w[1]]:
                ax5.axvline(xv, color=GREEN, lw=0.9, ls=":", alpha=0.8)

    # ── Panel 6: Resumen ─────────────────────────────────────────────────
    ax6 = fig.add_subplot(gs[2, 2])
    ax6.set_facecolor(BG2)
    ax6.axis("off")
    ax6.set_title("Resumen Ventanas Óptimas", color=TEXT_MAIN, fontsize=10)
    y_pos = 0.92
    for pk in proto_list:
        proto = PROTOCOLS[pk]
        m     = metrics[pk]
        w     = m["window"]
        pp    = m["peak_purity"] * 100
        pi    = m["purity_index"] * 100
        col   = proto["color"]
        if w[0] is not None:
            estado = f"✓  Días {w[0]}–{w[1]}"
            ec = GREEN
        else:
            estado = f"⚠  Pico día {w[1]} ({pp:.1f}%)"
            ec = ORANGE
        ax6.text(0.02, y_pos,       f"▌ {proto['label']}", transform=ax6.transAxes,
                 color=col, fontsize=9, fontweight="bold")
        ax6.text(0.02, y_pos-0.055, f"   {proto['funcion']}", transform=ax6.transAxes,
                 color=TEXT_DIM, fontsize=7.5)
        ax6.text(0.02, y_pos-0.11,  f"   {estado}", transform=ax6.transAxes,
                 color=ec, fontsize=8)
        ax6.text(0.02, y_pos-0.165, f"   Pico: {pp:.1f}%  |  Índice: {pi:.1f}%",
                 transform=ax6.transAxes, color=TEXT_DIM, fontsize=7.5)
        ax6.plot([0.02, 0.98], [y_pos-0.20, y_pos-0.20],
                 color=BORDER, lw=0.5, transform=ax6.transAxes)
        y_pos -= 0.28

    idx21    = np.argmin(np.abs(metrics["PODO"]["t"] - 21))
    risk21   = metrics["PODO"]["teratoma"][idx21] * 100
    col_risk = RED if risk21 > 5 else ORANGE if risk21 > 1 else GREEN
    ax6.text(0.02, y_pos,  "⚠ Teratoma día 21:", transform=ax6.transAxes,
             color=TEXT_DIM, fontsize=8)
    ax6.text(0.55, y_pos, f"{risk21:.4f}%", transform=ax6.transAxes,
             color=col_risk, fontsize=9, fontweight="bold")

    fig.text(0.5, 0.01,
             "VirtusSapiens · Bio-Kidney AI 2026  |  Ref: Takasato 2015, Freedman 2015  |"
             "  simulador_diferenciacion_ipsc.py",
             ha="center", color=TEXT_DIM, fontsize=7.5)
    return fig

# ─────────────────────────────────────────────────────────────
#  MAIN
# ─────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("\n[Bio-Kidney AI 2026] Iniciando simulación — 3 protocolos independientes...")
    results = run_all_protocols(t_end=30, n_points=3000)
    metrics = compute_metrics(results)
    print("[Bio-Kidney AI 2026] Simulación completada. Generando análisis...")
    print_report(metrics)
    fig = plot_all(results, metrics)
    output_path = "simulador_diferenciacion_ipsc_output.png"
    plt.tight_layout(rect=[0, 0.02, 1, 0.97])
    fig.savefig(output_path, dpi=150, bbox_inches="tight", facecolor=BG)
    print(f"[Bio-Kidney AI 2026] Figura guardada → {output_path}")
    plt.show()
