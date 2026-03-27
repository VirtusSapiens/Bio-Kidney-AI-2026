"""
Módulo de Experto en Dinámica de Fluidos BioKidney-AI
----------------------------------------------------
Modela la fenomenología física de la filtración glomerular y el transporte
tubular utilizando las ecuaciones de Starling y Deen.
"""

import numpy as np
from typing import Tuple, List
from biokidney.core.config import cfg_physio, cfg_sim

class FluidDynamicsExpert:
    """
    Motor experto en transporte de masa y hemodinámica renal.
    """
    
    @staticmethod
    def starling_deen_model(x_norm: np.ndarray, pi_0: float, ff_local: float) -> np.ndarray:
        """Modelo Deen para la presión oncótica a lo largo del capilar."""
        denominador = 1.0 - ff_local * x_norm
        denominador = np.maximum(denominador, 0.01)
        return pi_0 / denominador

    @staticmethod
    def calculate_gfr_glom(pgc: float, pbs: float, pi_gc: float, pi_bs: float, kf: float) -> float:
        """Ecuación de Starling: TFG = Kf * [(Pgc - Pbs) - (pi_gc - pi_bs)]"""
        delta_p = (pgc - pbs) - (pi_gc - pi_bs)
        return kf * max(delta_p, 0.0)

    @staticmethod
    def michaelis_menten(tm: float, km: float, concentracion: float) -> float:
        """
        Calcula la tasa de transporte activo: J = Tm * C / (Km + C)
        Retorna la fracción de capacidad utilizada (0-1).
        """
        if tm <= 0: return 0.0
        j = tm * concentracion / (km + concentracion)
        return min(j / tm, 1.0)

    @staticmethod
    def kedem_katchalsky_water(lp: float, delta_osm: float, vol_in: float) -> float:
        """
        Calcula el flujo de agua pasivo: Jv = Lp * delta_osm.
        Retorna el volumen reabsorbido (mL/min).
        """
        jv = lp * abs(delta_osm) * vol_in
        return min(jv, vol_in * 0.85)

    def integrate_capillary_filtration(self, pgc: float, pi_0: float) -> Tuple[float, float, float]:
        """
        Realiza la integración numérica (RK4 o Trapecio) a lo largo del capilar.
        Retorna: (TFG_glom, FF_final, dP_medio)
        """
        x_norm = np.linspace(0, 1, cfg_sim.N_PASOS_INTEGRACION)
        tfg_glom = 0.0
        ff_est = 0.20
        
        for xn in x_norm:
            pi_local = self.starling_deen_model(xn, pi_0, ff_est)
            dq = self.calculate_gfr_glom(pgc, cfg_physio.PBS, pi_local, cfg_physio.PI_BS, cfg_physio.KF_GLOM_NL)
            tfg_glom += dq * (1.0 / cfg_sim.N_PASOS_INTEGRACION)
            
            # Autoregresión de Fracción de Filtración
            q_plasma_nl = (cfg_physio.RPF_RIÑON / cfg_physio.N_GLOMERULOS) * 1e6
            ff_est = min(tfg_glom / q_plasma_nl, 0.9) if q_plasma_nl > 0 else 0.2
            
        return tfg_glom, ff_est, pgc - cfg_physio.PBS - pi_0
