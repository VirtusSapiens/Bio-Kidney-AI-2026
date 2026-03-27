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
    
    # Composición Plasmática (mEq/L o mg/dL)
    PLASMA_NA: float = 140.0
    PLASMA_K: float = 4.5
    PLASMA_GLUCOSA: float = 90.0
    PLASMA_CREATININA: float = 1.0
    PLASMA_OSM: float = 290.0
    
    # Capacidades de Transporte (Tm en mg/min o µEq/min)
    TM_SGLT2: float = 375.0
    TM_NHE3: float = 200.0
    TM_NKCC2: float = 180.0
    TM_ENAC: float = 30.0
    
    # Constantes de Afinidad (Km)
    KM_SGLT2: float = 2.0
    KM_NHE3: float = 15.0
    KM_NKCC2: float = 20.0
    KM_ENAC: float = 5.0
    
    # Permeabilidades Hidráulicas (Lp en mL/min/mOsm)
    LP_TP: float = 0.0045
    LP_AHD: float = 0.0035
    LP_AHA: float = 0.0001
    LP_TD: float = 0.0008
    LP_TC: float = 0.0030
    
    # Flujos y Volúmenes
    N_GLOMERULOS: int = 1_000_000
    RPF_RIÑON: float = 625.0
    GFR_TARGET: float = 62.5
    
    # Umbrales Clínicos
    UMBRAL_DIALISIS: float = 60.0
    UMBRAL_OPTIMO: float = 90.0
    
    # Parámetros de Difusión y Metabolismo (O2)
    D_O2: float = 2.0e-5
    M_MAX_O2: float = 5.0e-4
    P50_O2: float = 1.0
    P_HIPOXIA: float = 1.0
    P_ART_O2: float = 40.0
    P_VEN_O2: float = 20.0

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
