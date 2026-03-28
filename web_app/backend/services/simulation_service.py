"""
Servicio de Simulación BioKidney-AI — Capa de orquestación completa.
Expone cada módulo del pipeline como un método independiente con resultados estructurados.
"""

import time
import numpy as np
from typing import Dict, Any, List
from scipy.integrate import solve_ivp

from web_app.backend.utils.logger import bk_logger
from web_app.backend.database.models import SessionLocal, SimulationRecord
from biokidney.aggregator import BioKidneyEngine
from biokidney.core.config import cfg_physio, cfg_vasc, cfg_sim
from biokidney.experts.fluids import FluidDynamicsExpert
from biokidney.experts.cellular import CellularExpert


class SimulationService:
    """Orquesta las simulaciones del core biokidney y persiste resultados."""

    def __init__(self):
        self.engine = BioKidneyEngine()
        self.fluids = FluidDynamicsExpert()
        self.cellular = CellularExpert()

    # ──────────────────────────────────────────────────────────────
    # 1. VASCULAR (CCO)
    # ──────────────────────────────────────────────────────────────
    async def run_vascular(self, params: Dict[str, Any]) -> Dict[str, Any]:
        start = time.time()
        bk_logger.info("Simulando sistema vascular CCO...")

        n_seeds = int(params.get("n_seeds", 1000))
        murray_exp = float(params.get("murray_exponent", 3.0))

        np.random.seed(42)
        n_segments = int(n_seeds * 1.448)

        # Generate tree via recursive bifurcation with Murray's Law enforcement
        # Track parent-child triplets for compliance verification
        radii = []
        levels = []
        triplets = []  # (r_parent, r_child1, r_child2)

        def bifurcate(r_parent, level, max_level=11):
            if level >= max_level or r_parent < cfg_vasc.RADIO_MIN_UM:
                radii.append(r_parent)
                levels.append(level)
                return
            radii.append(r_parent)
            levels.append(level)
            # Murray's Law: r_p^alpha = r_1^alpha + r_2^alpha
            asym = np.clip(0.5 + np.random.normal(0, 0.08), 0.3, 0.7)
            r1 = r_parent * asym**(1.0 / murray_exp)
            r2 = r_parent * (1.0 - asym)**(1.0 / murray_exp)
            triplets.append((r_parent, r1, r2))
            bifurcate(r1, level + 1, max_level)
            bifurcate(r2, level + 1, max_level)

        # Generate for 3 systems (arterial, venous, collector)
        for root_r in [500, 600, 200]:
            bifurcate(float(root_r), 0)

        radii = np.array(radii[:n_segments])
        levels = np.array(levels[:n_segments])
        n_segments = len(radii)
        radii = np.clip(radii, cfg_vasc.RADIO_MIN_UM, cfg_vasc.RADIO_MAX_UM)

        # Murray compliance: verify actual parent-child triplets
        murray_ok = 0
        for rp, r1, r2 in triplets:
            ratio = (r1**murray_exp + r2**murray_exp) / (rp**murray_exp + 1e-12)
            if 0.85 <= ratio <= 1.15:
                murray_ok += 1
        murray_pct = round((murray_ok / max(len(triplets), 1)) * 100, 1)

        systems = {"arterial": int(n_segments * 0.334),
                    "venoso": int(n_segments * 0.333),
                    "colector": n_segments - int(n_segments * 0.334) - int(n_segments * 0.333)}

        radius_dist = {
            "labels": [f"{int(b)}" for b in np.linspace(cfg_vasc.RADIO_MIN_UM, cfg_vasc.RADIO_MAX_UM, 20)],
            "counts": np.histogram(radii, bins=20)[0].tolist()
        }

        level_dist = {
            "labels": [str(i) for i in range(1, 12)],
            "counts": [int(np.sum(levels == i)) for i in range(1, 12)]
        }

        duration = (time.time() - start) * 1000
        result = {
            "module": "vascular",
            "n_segments": n_segments,
            "n_seeds": n_seeds,
            "murray_compliance_pct": round(murray_pct, 1),
            "murray_exponent": murray_exp,
            "systems": systems,
            "domain_cm": {"x": cfg_vasc.DOMINIO_X, "y": cfg_vasc.DOMINIO_Y, "z": cfg_vasc.DOMINIO_Z},
            "radius_range_um": {"min": round(float(radii.min()), 1), "max": round(float(radii.max()), 1)},
            "radius_distribution": radius_dist,
            "level_distribution": level_dist,
            "status": "OPTIMO" if murray_pct >= 95 else "FUNCIONAL" if murray_pct >= 80 else "REVISAR",
            "duration_ms": round(duration, 2)
        }
        self._persist(result, "Vascular CCO", params, duration)
        return result

    # ──────────────────────────────────────────────────────────────
    # 2. FILTRACION GLOMERULAR
    # ──────────────────────────────────────────────────────────────
    async def run_filtration(self, params: Dict[str, Any]) -> Dict[str, Any]:
        start = time.time()
        bk_logger.info("Simulando filtración glomerular...")

        pgc = float(params.get("pgc", cfg_physio.PGC_ENTRY))
        pbs = float(params.get("pbs", cfg_physio.PBS))
        pi_gc = float(params.get("pi_gc", cfg_physio.PI_GC_ENTRY))
        n_glom = int(params.get("n_glomerulos", cfg_physio.N_GLOMERULOS))
        kf = float(params.get("kf", cfg_physio.KF_GLOM_NL))

        # Generate glomerular pressure distribution
        np.random.seed(42)
        n_sample = min(n_glom, 1000)
        frac_cortex = np.random.beta(2, 1, n_sample)
        p_terminals = np.clip(pgc + (frac_cortex - 0.5) * 9.0 + np.random.normal(0, 2, n_sample), 48, 75)

        tfg_per_glom = np.zeros(n_sample)
        ff_per_glom = np.zeros(n_sample)
        dp_per_glom = np.zeros(n_sample)

        # Starling pressure profile along capillary
        x_norm = np.linspace(0, 1, 100)
        pressure_profile = {"x": (x_norm * 100).tolist()}

        for i in range(n_sample):
            tfg_g, ff_g, dp_g = self.fluids.integrate_capillary_filtration(p_terminals[i], pi_gc)
            tfg_per_glom[i] = tfg_g
            ff_per_glom[i] = ff_g
            dp_per_glom[i] = dp_g

        tfg_total = float(np.mean(tfg_per_glom) * 1e-6 * n_glom)
        ff_mean = float(np.mean(ff_per_glom))
        dp_mean = float(np.mean(dp_per_glom))
        kf_eff = tfg_total / dp_mean if dp_mean > 0 else kf

        # Starling profile for reference case
        pi_profile = self.fluids.starling_deen_model(x_norm, pi_gc, 0.20)
        dp_profile = np.maximum((pgc - pbs) - (pi_profile - 0.0), 0)
        pressure_profile["p_hidrostatic"] = (np.full(100, pgc - pbs)).tolist()
        pressure_profile["pi_oncotic"] = pi_profile.tolist()
        pressure_profile["dp_starling"] = dp_profile.tolist()

        # TFG vs arterial pressure curve
        p_art_range = np.linspace(60, 160, 40)
        tfg_curve_vals = []
        for p_art in p_art_range:
            pgc_est = np.clip(p_art * 0.60, 30, 80)
            tfg_i, _, _ = self.fluids.integrate_capillary_filtration(pgc_est, pi_gc)
            tfg_curve_vals.append(round(float(tfg_i * 1e-6 * n_glom), 2))

        tfg_curve = {"p_arterial": p_art_range.tolist(), "tfg": tfg_curve_vals}

        # TFG histogram
        tfg_hist_counts, tfg_hist_edges = np.histogram(tfg_per_glom, bins=30)
        tfg_histogram = {
            "centers": ((tfg_hist_edges[:-1] + tfg_hist_edges[1:]) / 2).tolist(),
            "counts": tfg_hist_counts.tolist()
        }

        pct_native = (tfg_total / 62.5) * 100
        if tfg_total >= 90:
            status = "OPTIMO"
        elif tfg_total >= 60:
            status = "FUNCIONAL"
        else:
            status = "INSUFICIENTE"

        duration = (time.time() - start) * 1000
        result = {
            "module": "filtration",
            "tfg_mL_min": round(tfg_total, 2),
            "ff_mean_pct": round(ff_mean * 100, 1),
            "dp_starling_mean": round(dp_mean, 2),
            "kf_effective": round(kf_eff, 3),
            "pct_native": round(pct_native, 1),
            "n_glomerulos": n_glom,
            "pgc_entry": pgc,
            "pi_gc_entry": pi_gc,
            "pressure_profile": pressure_profile,
            "tfg_curve": tfg_curve,
            "tfg_histogram": tfg_histogram,
            "threshold_dialysis": 60.0,
            "threshold_optimal": 90.0,
            "status": status,
            "duration_ms": round(duration, 2)
        }
        self._persist(result, "Filtración Glomerular", params, duration)
        return result

    # ──────────────────────────────────────────────────────────────
    # 3. REABSORCION TUBULAR
    # ──────────────────────────────────────────────────────────────
    async def run_reabsorption(self, params: Dict[str, Any]) -> Dict[str, Any]:
        start = time.time()
        bk_logger.info("Simulando reabsorción tubular...")

        gfr_input = float(params.get("gfr_input", 82.0))
        aldosterone = float(params.get("aldosterone_factor", 1.0))
        adh = float(params.get("adh_factor", 1.0))

        vol = gfr_input
        segments = []

        # Túbulo Proximal
        frac_sglt2 = self.fluids.michaelis_menten(cfg_physio.TM_SGLT2, cfg_physio.KM_SGLT2, cfg_physio.PLASMA_GLUCOSA / 18.0)
        frac_nhe3 = self.fluids.michaelis_menten(cfg_physio.TM_NHE3, cfg_physio.KM_NHE3, cfg_physio.PLASMA_NA)
        tp_vol_out = vol * 0.33
        segments.append({
            "name": "Tubulo Proximal", "abbr": "TP",
            "vol_in": round(vol, 2), "vol_out": round(tp_vol_out, 2),
            "vol_reabs": round(vol - tp_vol_out, 2), "frac_reabs": 67.0,
            "transporters": {"SGLT2": round(frac_sglt2 * 100, 1), "NHE3": round(frac_nhe3 * 100, 1)},
            "osm_out": 296, "na_out": 139
        })

        # Asa Descendente
        ahd_vol_out = tp_vol_out * 0.27
        segments.append({
            "name": "Asa Descendente", "abbr": "AHD",
            "vol_in": round(tp_vol_out, 2), "vol_out": round(ahd_vol_out, 2),
            "vol_reabs": round(tp_vol_out - ahd_vol_out, 2), "frac_reabs": 73.0,
            "transporters": {"AQP1": 90.0},
            "osm_out": 1200, "na_out": 350
        })

        # Asa Ascendente
        frac_nkcc2 = self.fluids.michaelis_menten(cfg_physio.TM_NKCC2, cfg_physio.KM_NKCC2, 350)
        aha_vol_out = ahd_vol_out  # impermeable
        segments.append({
            "name": "Asa Ascendente", "abbr": "AHA",
            "vol_in": round(ahd_vol_out, 2), "vol_out": round(aha_vol_out, 2),
            "vol_reabs": 0.0, "frac_reabs": 0.0,
            "transporters": {"NKCC2": round(frac_nkcc2 * 100, 1)},
            "osm_out": 100, "na_out": 65
        })

        # Túbulo Distal
        frac_enac = self.fluids.michaelis_menten(cfg_physio.TM_ENAC, cfg_physio.KM_ENAC, 65) * aldosterone
        td_reabs = self.fluids.kedem_katchalsky_water(cfg_physio.LP_TD, 20.0, aha_vol_out)
        td_vol_out = aha_vol_out - td_reabs
        segments.append({
            "name": "Tubulo Distal", "abbr": "TD",
            "vol_in": round(aha_vol_out, 2), "vol_out": round(td_vol_out, 2),
            "vol_reabs": round(td_reabs, 2), "frac_reabs": round((td_reabs / aha_vol_out) * 100, 1) if aha_vol_out > 0 else 0,
            "transporters": {"ENaC": round(min(frac_enac, 1.0) * 100, 1)},
            "osm_out": 120, "na_out": 45
        })

        # Túbulo Colector
        lp_tc = cfg_physio.LP_TC * adh
        tc_reabs = self.fluids.kedem_katchalsky_water(lp_tc, 800.0, td_vol_out)
        tc_vol_out = max(td_vol_out - tc_reabs, 0.8)
        osm_orina = 100 + (adh * 1100)
        segments.append({
            "name": "Tubulo Colector", "abbr": "TC",
            "vol_in": round(td_vol_out, 2), "vol_out": round(tc_vol_out, 3),
            "vol_reabs": round(td_vol_out - tc_vol_out, 2),
            "frac_reabs": round(((td_vol_out - tc_vol_out) / td_vol_out) * 100, 1) if td_vol_out > 0 else 0,
            "transporters": {"AQP2": round(min(adh * 85, 100), 1)},
            "osm_out": round(osm_orina), "na_out": 80
        })

        total_reabs = gfr_input - tc_vol_out
        reabs_pct = (total_reabs / gfr_input) * 100
        urine_L_day = tc_vol_out * 1440 / 1000

        # Osmotic profile
        osm_profile = {
            "labels": ["Filtrado", "TP", "AHD (punta)", "AHA", "TD", "TC (orina)"],
            "values": [290, 296, 1200, 100, 120, round(osm_orina)]
        }

        # Volume waterfall
        vol_waterfall = {
            "labels": ["Filtrado Primario"] + [s["name"] for s in segments],
            "values": [round(gfr_input, 2)] + [s["vol_out"] for s in segments]
        }

        status = "OPTIMO" if 1.0 <= tc_vol_out <= 2.0 else "FUNCIONAL" if tc_vol_out < 2.5 else "POLIURIA"

        duration = (time.time() - start) * 1000
        result = {
            "module": "reabsorption",
            "gfr_input": gfr_input,
            "urine_vol_mL_min": round(tc_vol_out, 3),
            "urine_vol_L_day": round(urine_L_day, 2),
            "urine_osm": round(osm_orina),
            "total_reabs_mL_min": round(total_reabs, 2),
            "reabs_pct": round(reabs_pct, 1),
            "segments": segments,
            "osm_profile": osm_profile,
            "vol_waterfall": vol_waterfall,
            "aldosterone_factor": aldosterone,
            "adh_factor": adh,
            "status": status,
            "duration_ms": round(duration, 2)
        }
        self._persist(result, "Reabsorción Tubular", params, duration)
        return result

    # ──────────────────────────────────────────────────────────────
    # 4. OXIGENO
    # ──────────────────────────────────────────────────────────────
    async def run_oxygen(self, params: Dict[str, Any]) -> Dict[str, Any]:
        start = time.time()
        bk_logger.info("Simulando difusión de oxígeno...")

        grid_size = int(params.get("grid_size", 40))
        p_art = float(params.get("p_art_o2", cfg_physio.P_ART_O2))
        p_ven = float(params.get("p_ven_o2", cfg_physio.P_VEN_O2))
        d_o2 = cfg_physio.D_O2

        np.random.seed(42)
        nx, ny = grid_size, grid_size

        # Analytical steady-state model based on vascular source proximity
        # Models Krogh cylinder diffusion from distributed vascular sources
        x = np.linspace(0, 1, nx)
        y = np.linspace(0, 1, ny)
        xx, yy = np.meshgrid(x, y, indexing='ij')

        # Generate vascular source positions (matching CCO v7 density)
        n_sources = max(60, nx * ny // 10)
        src_x = np.random.uniform(0.05, 0.95, n_sources)
        src_y = np.random.uniform(0.05, 0.95, n_sources)
        src_p = np.random.choice([p_art, p_ven, (p_art + p_ven) / 2],
                                  n_sources, p=[0.34, 0.33, 0.33])

        # Compute PO2 field as superposition of Krogh cylinders
        # Each source creates a radial diffusion zone
        krogh_radius = 0.15  # normalized coverage radius
        po2 = np.full((nx, ny), 2.0)  # baseline low PO2

        for sx, sy, sp in zip(src_x, src_y, src_p):
            dist = np.sqrt((xx - sx)**2 + (yy - sy)**2)
            contribution = sp * np.exp(-dist**2 / (2 * krogh_radius**2))
            po2 = np.maximum(po2, contribution)

        # Apply metabolic consumption using Michaelis-Menten
        consumption = self.cellular.oxygen_consumption_rate(po2)
        po2 = np.maximum(po2 - consumption * 0.3, 1.0)

        # Smooth with Gaussian-like averaging (3 passes)
        for _ in range(3):
            po2_smooth = po2.copy()
            po2_smooth[1:-1, 1:-1] = (
                po2[:-2, 1:-1] + po2[2:, 1:-1] +
                po2[1:-1, :-2] + po2[1:-1, 2:] +
                4 * po2[1:-1, 1:-1]
            ) / 8
            po2 = po2_smooth

        po2_min = float(po2.min())
        po2_max = float(po2.max())
        po2_mean = float(po2.mean())
        hypoxia_pct = float(np.sum(po2 < cfg_physio.P_HIPOXIA) / po2.size * 100)

        # Heatmap data
        step = max(1, grid_size // 30)
        heatmap = po2[::step, ::step].round(2).tolist()

        # PO2 histogram
        po2_flat = po2.flatten()
        hist_counts, hist_edges = np.histogram(po2_flat, bins=25)
        po2_histogram = {
            "centers": ((hist_edges[:-1] + hist_edges[1:]) / 2).round(1).tolist(),
            "counts": hist_counts.tolist()
        }

        # Cross-section profile
        mid = nx // 2
        cross_section = {
            "x": list(range(ny)),
            "po2": po2[mid, :].round(2).tolist()
        }

        status = "OPTIMO" if hypoxia_pct == 0 else "FUNCIONAL" if hypoxia_pct < 5 else "HIPOXIA"

        duration = (time.time() - start) * 1000
        result = {
            "module": "oxygen",
            "grid_size": grid_size,
            "po2_min": round(po2_min, 2),
            "po2_max": round(po2_max, 2),
            "po2_mean": round(po2_mean, 2),
            "hypoxia_pct": round(hypoxia_pct, 2),
            "hypoxia_threshold": cfg_physio.P_HIPOXIA,
            "heatmap": heatmap,
            "po2_histogram": po2_histogram,
            "cross_section": cross_section,
            "convergence": [],
            "iterations": 3,
            "p_art_o2": p_art,
            "p_ven_o2": p_ven,
            "status": status,
            "duration_ms": round(duration, 2)
        }
        self._persist(result, "Difusión O2", params, duration)
        return result

    # ──────────────────────────────────────────────────────────────
    # 5. iPSC DIFERENCIACION
    # ──────────────────────────────────────────────────────────────
    async def run_ipsc(self, params: Dict[str, Any]) -> Dict[str, Any]:
        start = time.time()
        bk_logger.info("Simulando diferenciación iPSC...")

        t_end = float(params.get("t_end", 30))
        n_points = 500

        PROTOCOLS = {
            "PODO": {"label": "Podocitos", "function": "Filtracion glomerular",
                     "color": "#F0883E", "k_npc": 1.40, "p_t": 0.06, "d_t": 0.02},
            "PCT":  {"label": "Tub. Proximales", "function": "Reabsorcion tubular",
                     "color": "#39D353", "k_npc": 1.20, "p_t": 0.08, "d_t": 0.02},
            "IEC":  {"label": "Endoteliales (iECs)", "function": "Vascularizacion",
                     "color": "#79C0FF", "k_npc": 1.00, "p_t": 0.10, "d_t": 0.02},
        }

        K_COMMON = {
            "k_iPSC_to_IM": 1.20, "k_IM_to_NPC": 0.90,
            "p_iPSC": 0.08, "p_NPC": 0.18,
            "d_iPSC": 0.60, "d_IM": 0.10, "d_NPC": 0.04, "K_total": 1.0,
        }

        def hill(c, kd, n=2):
            return (c**n) / (kd**n + c**n)

        def wnt_signal(t):
            if t < 0: return 0.0
            if t > 5: return max(0.0, 1.0 - 0.4 * (t - 5))
            return min(1.0, t / 0.8) * hill(8.0, 2.0)

        def bmp_ra_signal(t):
            if t < 5: return 0.0
            if t > 12: return max(0.0, 1.0 - 0.4 * (t - 12))
            return min(1.0, (t - 5) / 0.8) * 0.9

        def spec_signal(t, start=12, end=28):
            if t < start: return 0.0
            if t > end: return max(0.0, 1.0 - 0.4 * (t - end))
            return min(1.0, (t - start) / 0.8) * 0.85

        def ode(t, y, kc, proto):
            ipsc, im, npc, target = y
            total = max(ipsc + im + npc + target, 1e-12)
            logistic = max(0.0, 1.0 - total / kc["K_total"])
            wnt = wnt_signal(t)
            bmp_ra = bmp_ra_signal(t)
            spec = spec_signal(t)

            d_iPSC_IM = kc["k_iPSC_to_IM"] * wnt * ipsc
            d_IM_NPC = kc["k_IM_to_NPC"] * bmp_ra * 0.5 * im
            d_NPC_T = proto["k_npc"] * spec * npc

            return [
                (kc["p_iPSC"] * logistic - kc["d_iPSC"]) * ipsc - d_iPSC_IM,
                d_iPSC_IM - d_IM_NPC - kc["d_IM"] * im,
                d_IM_NPC + kc["p_NPC"] * logistic * npc - d_NPC_T - kc["d_NPC"] * npc,
                d_NPC_T + proto["p_t"] * logistic * target - proto["d_t"] * target,
            ]

        t_eval = np.linspace(0, t_end, n_points)
        protocols_result = []

        for pk, proto in PROTOCOLS.items():
            sol = solve_ivp(
                lambda t, y, p=proto: ode(t, y, K_COMMON, p),
                (0, t_end), [1.0, 0.0, 0.0, 0.0],
                t_eval=t_eval, method="RK45", rtol=1e-8, atol=1e-10, max_step=0.1
            )
            ipsc, im, npc, target = sol.y
            total = np.maximum(ipsc + im + npc + target, 1e-12)
            purity = target / total
            teratoma = ipsc / total

            peak_purity = float(np.max(purity))
            mask_95 = purity >= 0.95
            idx_95 = np.where(mask_95)[0]
            if len(idx_95):
                window = (round(float(sol.t[idx_95[0]]), 1), round(float(sol.t[idx_95[-1]]), 1))
            else:
                window = None

            # Downsample for JSON
            step = max(1, n_points // 150)
            protocols_result.append({
                "key": pk,
                "label": proto["label"],
                "function": proto["function"],
                "color": proto["color"],
                "peak_purity_pct": round(peak_purity * 100, 1),
                "window_95": window,
                "time": sol.t[::step].round(1).tolist(),
                "purity": (purity[::step] * 100).round(1).tolist(),
                "teratoma": (teratoma[::step] * 100).round(4).tolist(),
                "ipsc": (ipsc[::step] / total[::step] * 100).round(1).tolist(),
                "npc": (npc[::step] / total[::step] * 100).round(1).tolist(),
                "target": (target[::step] / total[::step] * 100).round(1).tolist(),
            })

        all_apt = all(p["window_95"] is not None for p in protocols_result)

        # Teratoma at key days
        teratoma_summary = {}
        for p in protocols_result:
            t_arr = np.array(p["time"])
            ter_arr = np.array(p["teratoma"])
            for day in [7, 14, 21, 28]:
                idx = int(np.argmin(np.abs(t_arr - day)))
                teratoma_summary[f"{p['key']}_day{day}"] = round(float(ter_arr[idx]), 4)

        status = "OPTIMO" if all_apt else "PARCIAL"

        duration = (time.time() - start) * 1000
        result = {
            "module": "ipsc",
            "t_end_days": t_end,
            "protocols": protocols_result,
            "teratoma_summary": teratoma_summary,
            "all_protocols_apt": all_apt,
            "status": status,
            "duration_ms": round(duration, 2)
        }
        self._persist(result, "Diferenciación iPSC", params, duration)
        return result

    # ──────────────────────────────────────────────────────────────
    # 6. CO-SWIFT BIOPRINTING
    # ──────────────────────────────────────────────────────────────
    async def run_bioprinting(self, params: Dict[str, Any]) -> Dict[str, Any]:
        start = time.time()
        bk_logger.info("Simulando optimización Co-SWIFT bioprinting...")

        n_particles = int(params.get("n_particles", 100))
        np.random.seed(42)

        # Herschel-Bulkley rheology parameters
        tau_y = 5.0     # yield stress Pa
        K_hb = 0.8      # consistency index
        n_hb = 0.45     # flow behavior index
        r_nozzle = 200e-6  # 200 um

        # Generate Pareto-optimal solutions
        pressures = np.random.uniform(20, 100, n_particles)  # kPa
        viabilities = np.zeros(n_particles)
        shear_stresses = np.zeros(n_particles)
        flow_rates = np.zeros(n_particles)

        for i in range(n_particles):
            p = pressures[i] * 1000  # Pa
            gamma_dot = (p * r_nozzle) / (2 * K_hb * 0.01)  # shear rate estimate
            tau = tau_y + K_hb * (gamma_dot ** n_hb)
            shear_stresses[i] = tau
            viabilities[i] = max(0, min(100, 100 - 0.12 * tau))
            flow_rates[i] = 3.14159 * r_nozzle**2 * gamma_dot * 1e9  # nL/s

        # Find Pareto front
        pareto_mask = np.ones(n_particles, dtype=bool)
        for i in range(n_particles):
            for j in range(n_particles):
                if i != j:
                    if pressures[j] <= pressures[i] and viabilities[j] >= viabilities[i] and \
                       (pressures[j] < pressures[i] or viabilities[j] > viabilities[i]):
                        pareto_mask[i] = False
                        break

        optimal_idx = np.argmax(viabilities[pareto_mask] * (1 - np.abs(pressures[pareto_mask] - 60) / 100))
        if np.sum(pareto_mask) > 0:
            pareto_pressures = pressures[pareto_mask]
            pareto_viabilities = viabilities[pareto_mask]
            best_pressure = float(pareto_pressures[min(optimal_idx, len(pareto_pressures) - 1)])
            best_viability = float(pareto_viabilities[min(optimal_idx, len(pareto_viabilities) - 1)])
        else:
            best_pressure = 60.0
            best_viability = 98.0

        best_idx = np.argmin(np.abs(pressures - best_pressure))
        best_wss = float(shear_stresses[best_idx])

        # WSS in dyn/cm2
        wss_dyncm2 = best_wss * 0.1

        scatter_data = {
            "pressures": pressures.round(1).tolist(),
            "viabilities": viabilities.round(1).tolist(),
            "shear_stresses": shear_stresses.round(1).tolist(),
            "is_pareto": pareto_mask.tolist()
        }

        # Pressure vs viability curve
        p_range = np.linspace(20, 100, 50)
        via_curve = [round(float(max(0, min(100, 100 - 0.12 * (tau_y + K_hb * ((p * 1000 * r_nozzle / (2 * K_hb * 0.01)) ** n_hb))))), 1)
                     for p in p_range]
        pv_curve = {"pressure": p_range.round(1).tolist(), "viability": via_curve}

        status = "OPTIMO" if best_viability >= 95 else "FUNCIONAL" if best_viability >= 85 else "REVISAR"

        duration = (time.time() - start) * 1000
        result = {
            "module": "bioprinting",
            "n_particles": n_particles,
            "optimal_pressure_kPa": round(best_pressure, 1),
            "optimal_viability_pct": round(best_viability, 1),
            "optimal_wss_Pa": round(best_wss, 1),
            "wss_dyncm2": round(wss_dyncm2, 1),
            "nozzle_radius_um": 200,
            "bioink": "GelMA 7% + Alginato 1.5% + Nanocelulosa 0.8%",
            "scatter": scatter_data,
            "pv_curve": pv_curve,
            "pareto_count": int(np.sum(pareto_mask)),
            "pressure_range": {"min": 20, "max": 100, "unit": "kPa"},
            "viability_threshold": 85.0,
            "status": status,
            "duration_ms": round(duration, 2)
        }
        self._persist(result, "Co-SWIFT Bioprinting", params, duration)
        return result

    # ──────────────────────────────────────────────────────────────
    # FULL PIPELINE
    # ──────────────────────────────────────────────────────────────
    async def run_full_pipeline(self, user_params: Dict[str, Any]) -> Dict[str, Any]:
        start = time.time()
        bk_logger.info(f"Iniciando pipeline completo con params: {user_params}")

        vasc = await self.run_vascular(user_params.get("vascular", {}))
        filt = await self.run_filtration(user_params.get("filtration", {}))
        reabs = await self.run_reabsorption({"gfr_input": filt["tfg_mL_min"],
                                              **user_params.get("reabsorption", {})})
        oxy = await self.run_oxygen(user_params.get("oxygen", {}))
        ipsc = await self.run_ipsc(user_params.get("ipsc", {}))
        bio = await self.run_bioprinting(user_params.get("bioprinting", {}))

        duration = (time.time() - start) * 1000

        modules_status = {
            "vascular": vasc["status"],
            "filtration": filt["status"],
            "reabsorption": reabs["status"],
            "oxygen": oxy["status"],
            "ipsc": ipsc["status"],
            "bioprinting": bio["status"],
        }
        n_optimal = sum(1 for s in modules_status.values() if s == "OPTIMO")
        global_status = "OPTIMO" if n_optimal >= 5 else "FUNCIONAL" if n_optimal >= 3 else "REVISAR"

        return {
            "pipeline": "complete",
            "global_status": global_status,
            "modules_status": modules_status,
            "n_optimal": n_optimal,
            "summary": {
                "tfg": filt["tfg_mL_min"],
                "urine": reabs["urine_vol_mL_min"],
                "reabs_pct": reabs["reabs_pct"],
                "hypoxia_pct": oxy["hypoxia_pct"],
                "viability_pct": bio["optimal_viability_pct"],
                "murray_pct": vasc["murray_compliance_pct"],
                "ipsc_apt": ipsc["all_protocols_apt"],
            },
            "vascular": vasc,
            "filtration": filt,
            "reabsorption": reabs,
            "oxygen": oxy,
            "ipsc": ipsc,
            "bioprinting": bio,
            "duration_ms": round(duration, 2)
        }

    # ──────────────────────────────────────────────────────────────
    # HISTORY
    # ──────────────────────────────────────────────────────────────
    def get_history(self, limit: int = 20):
        db = SessionLocal()
        try:
            records = db.query(SimulationRecord).order_by(
                SimulationRecord.timestamp.desc()
            ).limit(limit).all()
            return [
                {
                    "id": r.id,
                    "timestamp": r.timestamp.isoformat() if r.timestamp else None,
                    "type": r.tipo_simulacion,
                    "status": r.estado,
                    "duration_ms": r.duracion_ms,
                    "results": r.resultados,
                }
                for r in records
            ]
        finally:
            db.close()

    # ──────────────────────────────────────────────────────────────
    # PERSISTENCE
    # ──────────────────────────────────────────────────────────────
    def _persist(self, result: dict, sim_type: str, params: dict, duration: float):
        try:
            db = SessionLocal()
            # Remove large data before persisting
            persist_result = {k: v for k, v in result.items()
                             if k not in ("heatmap", "pressure_profile", "scatter",
                                          "protocols", "tfg_histogram", "po2_histogram",
                                          "cross_section", "convergence", "pv_curve",
                                          "tfg_curve", "radius_distribution", "level_distribution",
                                          "vol_waterfall", "osm_profile")}
            record = SimulationRecord(
                tipo_simulacion=sim_type,
                parametros_entrada=params,
                resultados=persist_result,
                duracion_ms=duration,
                estado=result.get("status", "COMPLETED")
            )
            db.add(record)
            db.commit()
        except Exception as e:
            bk_logger.error(f"Error persistiendo resultado: {e}")
        finally:
            db.close()
