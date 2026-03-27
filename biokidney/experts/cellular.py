"""
Módulo de Experto Celular BioKidney-AI
--------------------------------------
Modelado de procesos biológicos a nivel celular: consumo metabólico, 
difusión de gases y cinética de diferenciación (iPSC).
"""

import numpy as np
from biokidney.core.config import cfg_physio

class CellularExpert:
    """
    Motor experto en fenomenología celular y metabolismo.
    """
    
    @staticmethod
    def oxygen_consumption_rate(po2: np.ndarray) -> np.ndarray:
        """
        Calcula la tasa de consumo de O2 usando Michaelis-Menten.
        M = M_max * P / (P + P50)
        """
        pp = np.maximum(po2, 0.0)
        return cfg_physio.M_MAX_O2 * pp / (pp + cfg_physio.P50_O2 + 1e-12)

    def evaluate_hypoxia(self, po2_grid: np.ndarray) -> float:
        """Calcula el porcentaje de volumen en estado de hipoxia."""
        vóxeles_hipoxicos = np.sum(po2_grid < cfg_physio.P_HIPOXIA)
        return (vóxeles_hipoxicos / po2_grid.size) * 100.0

    def simulate_ipsc_maturation(self, days: np.ndarray) -> np.ndarray:
        """Modelado de la curva de maduración de iPSCs."""
        # Lógica basada en simulador_ipsc_biokidney.py
        return 1.0 / (1.0 + np.exp(-0.3 * (days - 15)))
