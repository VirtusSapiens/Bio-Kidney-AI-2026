"""
Módulo de Configuración Centralizada BioKidney-AI
------------------------------------------------
Este módulo define las constantes fisiológicas y parámetros de simulación
para asegurar la consistencia en todo el ecosistema.
"""

from dataclasses import dataclass, field
import numpy as np

@dataclass(frozen=True)
class PhysiologicalConfig:
    # Presiones Glomerulares (mmHg)
    PGC_ENTRY: float = 60.0
    PBS: float = 15.0
    PI_GC_ENTRY: float = 28.0
    PI_BS: float = 0.0
    
    # Coeficientes de Filtración
    KF_REF: float = 3.7  # mL/min/mmHg (Riñón completo)
    KF_GLOM_NL: float = 3.7  # nL/min/mmHg (Por glomérulo)
    
    # Flujos y Volúmenes
    N_GLOMERULOS: int = 1_000_000
    RPF_RIÑON: float = 625.0  # mL/min (Flujo Plasmático Renal)
    GFR_TARGET: float = 62.5  # mL/min por riñón
    
    # Umbrales Clínicos
    UMBRAL_DIALISIS: float = 60.0
    UMBRAL_OPTIMO: float = 90.0

@dataclass(frozen=True)
class VascularConfig:
    # Ley de Murray
    MURRAY_EXPONENT: float = 3.0
    
    # Geometría
    RADIO_MAX_UM: float = 600.0
    RADIO_MIN_UM: float = 15.0
    
    # Dominio Renal (cm)
    DOMINIO_X: float = 5.5
    DOMINIO_Y: float = 3.0
    DOMINIO_Z: float = 2.5

@dataclass(frozen=True)
class SimulationConfig:
    # Resolución
    N_PASOS_INTEGRACION: int = 100
    N_SEMILLAS_DEMANDA: int = 1000
    
    # Colores Estándar
    COLOR_ARTERIAL: str = "#FF4444"
    COLOR_VENOSO: str = "#00D4FF"
    COLOR_COLECTOR: str = "#FFD700"
    COLOR_PARENQUIMA: str = "#161B22"

# Instancia Global de Configuración
cfg_physio = PhysiologicalConfig()
cfg_vasc = VascularConfig()
cfg_sim = SimulationConfig()
