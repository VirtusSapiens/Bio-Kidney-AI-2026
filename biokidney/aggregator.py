"""
Orquestador Principal BioKidney-AI (Patrón Agregación)
------------------------------------------------------
Este módulo actúa como el punto de integración central del framework,
agregando a los expertos especializados en cada dominio renal.
"""

from biokidney.experts.vascular import VascularExpert
from biokidney.experts.fluids import FluidDynamicsExpert
from biokidney.core.config import cfg_physio, cfg_vasc
import numpy as np

class BioKidneyEngine:
    """
    Agregador Maestro que orquesta la simulación completa del riñón artificial.
    """
    def __init__(self):
        # Agregación de Expertos (MoE)
        self.vascular = VascularExpert()
        self.fluids = FluidDynamicsExpert()
        
        # Estado Global del Órgano
        self.tfg_actual = 0.0
        self.orina_actual = 0.0
        self.viabilidad = 1.0

    def ejecutar_pipeline_completo(self):
        """
        Ejecuta la secuencia lógica de simulación:
        Vascularización -> Filtración -> Reabsorción
        """
        print("🚀 BioKidney Engine: Iniciando pipeline multiescala...")
        
        # 1. Generar Árboles
        hilio = np.array([0.0, -cfg_vasc.DOMINIO_Y * 0.85, 0.0])
        self.vascular.sintetizar_sistema(hilio, cfg_vasc.RADIO_MAX_UM * 1e-6, 'art')
        
        # 2. Simular Filtración (usando presiones del árbol)
        # Por ahora usamos la presión de configuración
        tfg_glom, _, _ = self.fluids.integrate_capillary_filtration(cfg_physio.PGC_ENTRY, cfg_physio.PI_GC_ENTRY)
        self.tfg_actual = tfg_glom * 1e-6 * cfg_physio.N_GLOMERULOS
        
        print(f"✅ Pipeline Completado. TFG Resultante: {self.tfg_actual:.2f} mL/min")
        return self.tfg_actual

if __name__ == "__main__":
    engine = BioKidneyEngine()
    engine.ejecutar_pipeline_completo()
