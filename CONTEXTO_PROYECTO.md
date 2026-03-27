# Informe de Contexto: BioKidney-AI
**Fecha de generación:** 2026-03-26 21:17:02
**Estado del Sistema:** Pipeline 100% Completado (Hito Marzo 2026)

## 1. Estructura de Directorios
```text
00_bitacora/
        BITACORA.md
        Bitacora_BioKidney_AI_2026.pdf
01_simuladores/
        BioKidney_AI_2026_Informe_Sesion_TFG.pdf
        BiokidneyAI_Informe_Avance_Mar2026.docx
        Code_Generated_Image.png
        Informe_Modulo09_BioKidney_AI_2026.docx
        dashboard_maestro (1).png
        dashboard_maestro.png
        informe_biokidney_ai_2026_Simulador_de_Reabsorción_Tubular.docx
        informe_biokidney_coswift.docx
        informe_sesion_biokidney_marzo2026.pdf
        optimizador_coswift.py
        reabsorcion_tubular_dashboard (1).png
        reabsorcion_tubular_dashboard.png
        reabsorcion_tubular_informe.pdf
        resultado_oxigeno_biokidney.png
        simulador_diferenciacion_ipsc.py
        simulador_diferenciacion_ipsc_output.png
        simulador_estres_mecanico_dECM.png
        simulador_estres_mecanico_dECM.py
        simulador_filtracion_glomerular (1).png
        simulador_filtracion_glomerular.png
        simulador_filtracion_glomerular.py
        simulador_filtracion_glomerular_G.py
        simulador_filtracion_glomerular_reporte.pdf
        simulador_oxigeno_biokidney.py
        simulador_reabsorcion_tubular.py
    v1_originales/
            simulador_ipsc_biokidney.py
            simulador_swift_biokidney.py
            simulador_wss_biokidney.py
    v2_dashboard/
02_vascular_cco/
        arbol_vascular_cco.csv
        arbol_vascular_cco.obj
        arbol_vascular_cco.png
        arbol_vascular_cco_v2.csv
        arbol_vascular_cco_v2.png
        arbol_vascular_cco_v3.csv
        arbol_vascular_cco_v3.png
        arbol_vascular_cco_v4.csv
        arbol_vascular_cco_v4.png
        arbol_vascular_cco_v5.csv
        arbol_vascular_cco_v5.png
        arbol_vascular_cco_v6.csv
        arbol_vascular_cco_v6.png
        arbol_vascular_cco_v7.csv
        arbol_vascular_cco_v7.obj
        arbol_vascular_cco_v7.png
        arbol_vascular_cco_v7_curvo.obj
        exportar_blender.py
        exportar_blender_v2.py
        generador_cco.py
        generador_cco_v2.py
        generador_cco_v3.py
        generador_cco_v4.py
        generador_cco_v5.py
        generador_cco_v6.py
        generador_cco_v7.py
        senyal_hipoxia_para_cco.csv
03_modelos_3d/
        V7_recto_untitled.blend
        arbol_vascular_cco.blend
        arbol_vascular_cco.blend1
        generar_rinon_completo.py
        render_biokidney_v7_definitivo.png
        rinon_biokidney_ai.blend
        rinon_biokidney_ai_v7.blend
        rinon_biokidney_ai_v7.blend1
        rinon_capsula.blend
        rinon_capsula.blend1
        solo_capsula.py
        untitled.png
        v7_curvountitled.blend
04_literatura/
05_presentacion/
    pitch_inversores/
06_app/
        BioKidney_AI_Informe_Desarrollo_dashboard_maestro.docx
        biokidney_app.py
        dashboard_maestro_app.py
        filtracion_glomerular_gui.py
07_Videos Luma/
        Arbol vascular CCO v7.mp4
        Vascular_Tree_Video_An_abstract_glowing_animation_displays_two_YtuUY7DV.mp4
07_presentacion_final/
        dashboard_maestro.png
        reporte_validacion_completo.pdf
        resumen_ejecutivo_una_pagina.pdf
        scorecard_biokidney_ai_2026.pdf
    BLUEPRINT_INGENIERIA.md
    CONTEXTO_PROYECTO.md
    Dockerfile
    README.md
    analizador_proyecto_biokidney.py
biokidney/
        __init__.py
        aggregator.py
    core/
            __init__.py
            config.py
    experts/
            __init__.py
            cellular.py
            fluids.py
            vascular.py
    utils/
            __init__.py
    biokidney_architect.py
    docker-compose.yml
    inicio.sh
    requirements.txt
resultados/
        pareto_coswift.png
web_app/
    admin/
    backend/
        core/
        database/
                models.py
            main.py
        services/
                simulation_service.py
        utils/
                logger.py
    frontend/
        components/
            dashboard.html
```

## 2. Resumen de Componentes Clave

### Archivo: `00_bitacora/BITACORA.md`
- **Tipo:** Documentación Markdown

---

### Archivo: `01_simuladores/optimizador_coswift.py`
- **Tipo:** Script de Python
- **Dependencias:** `from matplotlib.colors import LinearSegmentedColormap, from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, from scipy.optimize import minimize, from scipy.spatial import ConvexHull, import matplotlib, import matplotlib.gridspec as gridspec, import matplotlib.patheffects as pe, import matplotlib.pyplot as plt...`
- **Propósito:** ╔══════════════════════════════════════════════════════════════════════════════╗ ║          OPTIMIZADOR Co-SWIFT — Bio-Kidney AI 2026                         ║ ║          Módulo de Optimización Multi-Objetivo para Bioimpresión            ║ ║         ...

---

### Archivo: `01_simuladores/simulador_diferenciacion_ipsc.py`
- **Tipo:** Script de Python
- **Dependencias:** `from matplotlib.ticker import MultipleLocator, from scipy.integrate import solve_ivp, import matplotlib, import matplotlib.gridspec as gridspec, import matplotlib.pyplot as plt, import numpy as np, import warnings`
- **Propósito:** ╔══════════════════════════════════════════════════════════════════════════════╗ ║       Bio-Kidney AI 2026 — Simulador de Diferenciación iPSC                ║ ║       Módulo: Conversión iPSC → Linajes Renales Específicos                ║ ║       Ref...

---

### Archivo: `01_simuladores/simulador_estres_mecanico_dECM.py`
- **Tipo:** Script de Python
- **Dependencias:** `from matplotlib.colors import LinearSegmentedColormap, from matplotlib.patches import Circle, import matplotlib, import matplotlib.gridspec as gridspec, import matplotlib.pyplot as plt, import numpy as np, import warnings`
- **Propósito:** BIO-KIDNEY AI 2026 - VirtusSapiens Carlos David Moreno Caceres SIMULADOR DE ESTRES MECANICO EN dECM - Modulo 09 v5  FISICA IMPLEMENTADA: 1. Elasticidad lineal   : sigma = E * epsilon 2. Kelvin-Voigt         : epsilon(t) = (sigma/E)*(1 - exp(-E*t/eta)...

---

### Archivo: `01_simuladores/simulador_filtracion_glomerular.py`
- **Tipo:** Script de Python
- **Dependencias:** `from datetime import datetime, from matplotlib.colors import LinearSegmentedColormap, from scipy.integrate import solve_ivp, import matplotlib, import matplotlib.gridspec as gridspec, import matplotlib.pyplot as plt, import numpy as np, import os, sys, csv, warnings`
- **Propósito:** ╔══════════════════════════════════════════════════════════════════════════════════╗ ║          SIMULADOR DE FILTRACIÓN GLOMERULAR — Bio-Kidney AI 2026               ║ ║          Carlos David Moreno Cáceres — VirtusSapiens                           ║...

---

### Archivo: `01_simuladores/simulador_filtracion_glomerular_G.py`
- **Tipo:** Script de Python
- **Dependencias:** `from datetime import datetime, from matplotlib.colors import LinearSegmentedColormap, from scipy.integrate import solve_ivp, import matplotlib, import matplotlib.gridspec as gridspec, import matplotlib.pyplot as plt, import numpy as np, import os, sys, csv, warnings`
- **Propósito:** ╔══════════════════════════════════════════════════════════════════════════════════╗ ║          SIMULADOR DE FILTRACIÓN GLOMERULAR — Bio-Kidney AI 2026               ║ ║          Carlos David Moreno Cáceres — VirtusSapiens                           ║...

---

### Archivo: `01_simuladores/simulador_oxigeno_biokidney.py`
- **Tipo:** Script de Python
- **Dependencias:** `from biokidney.core.config import cfg_physio, cfg_sim, cfg_vasc, from biokidney.experts.cellular import CellularExpert, from matplotlib.colors import LinearSegmentedColormap, from pathlib import Path, import matplotlib, import matplotlib.gridspec as gridspec, import matplotlib.pyplot as plt, import numpy as np...`

---

### Archivo: `01_simuladores/simulador_reabsorcion_tubular.py`
- **Tipo:** Script de Python
- **Dependencias:** `from biokidney.aggregator import BioKidneyEngine, from biokidney.core.config import cfg_physio, cfg_sim, from biokidney.experts.fluids import FluidDynamicsExpert, from datetime import datetime, from matplotlib.gridspec import GridSpec, from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, from reportlab.lib import colors, from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT...`
- **Propósito:** ╔══════════════════════════════════════════════════════════════════════════════╗ ║         BIO-KIDNEY AI 2026 — VirtusSapiens © Carlos David Moreno Cáceres   ║ ║              SIMULADOR DE REABSORCIÓN TUBULAR — Módulo 12                   ║ ║         ...

---

### Archivo: `01_simuladores/v1_originales/simulador_ipsc_biokidney.py`
- **Tipo:** Script de Python
- **Dependencias:** `import matplotlib.pyplot as plt, import numpy as np`

---

### Archivo: `01_simuladores/v1_originales/simulador_swift_biokidney.py`
- **Tipo:** Script de Python
- **Dependencias:** `import matplotlib.pyplot as plt, import numpy as np`

---

### Archivo: `01_simuladores/v1_originales/simulador_wss_biokidney.py`
- **Tipo:** Script de Python
- **Dependencias:** `import matplotlib.pyplot as plt, import numpy as np`

---

### Archivo: `02_vascular_cco/arbol_vascular_cco.csv`
- **Tipo:** Dataset / Datos Estructurales
- **Columnas:** `id,x1,y1,z1,x2,y2,z2,radio_um,longitud_mm,reynolds`
- **Registros:** 161

---

### Archivo: `02_vascular_cco/arbol_vascular_cco_v2.csv`
- **Tipo:** Dataset / Datos Estructurales
- **Columnas:** `id,x1_mm,y1_mm,z1_mm,x2_mm,y2_mm,z2_mm,radio_um,longitud_mm`
- **Registros:** 362

---

### Archivo: `02_vascular_cco/arbol_vascular_cco_v3.csv`
- **Tipo:** Dataset / Datos Estructurales
- **Columnas:** `id,sistema,nivel,x1_mm,y1_mm,z1_mm,x2_mm,y2_mm,z2_mm,radio_um`
- **Registros:** 401

---

### Archivo: `02_vascular_cco/arbol_vascular_cco_v4.csv`
- **Tipo:** Dataset / Datos Estructurales
- **Columnas:** `id,sistema,nivel,x1_mm,y1_mm,z1_mm,x2_mm,y2_mm,z2_mm,radio_um`
- **Registros:** 302

---

### Archivo: `02_vascular_cco/arbol_vascular_cco_v5.csv`
- **Tipo:** Dataset / Datos Estructurales
- **Columnas:** `id,sistema,nivel,x1_mm,y1_mm,z1_mm,x2_mm,y2_mm,z2_mm,radio_um`
- **Registros:** 1078

---

### Archivo: `02_vascular_cco/arbol_vascular_cco_v6.csv`
- **Tipo:** Dataset / Datos Estructurales
- **Columnas:** `id,sistema,nivel,x1_mm,y1_mm,z1_mm,x2_mm,y2_mm,z2_mm,radio_um`
- **Registros:** 488

---

### Archivo: `02_vascular_cco/arbol_vascular_cco_v7.csv`
- **Tipo:** Dataset / Datos Estructurales
- **Columnas:** `id,sistema,nivel,x1_mm,y1_mm,z1_mm,x2_mm,y2_mm,z2_mm,radio_um`
- **Registros:** 1448

---

### Archivo: `02_vascular_cco/exportar_blender.py`
- **Tipo:** Script de Python
- **Dependencias:** `import csv, import math, import numpy as np, import os`
- **Propósito:** Genera vértices e índices de un cilindro entre dos puntos. Retorna (vertices, caras) donde caras son índices locales....

---

### Archivo: `02_vascular_cco/exportar_blender_v2.py`
- **Tipo:** Script de Python
- **Dependencias:** `import csv, import math, import numpy as np, import os`
- **Propósito:** Genera curva Bezier entre dos puntos...

---

### Archivo: `02_vascular_cco/generador_cco.py`
- **Tipo:** Script de Python
- **Dependencias:** `from mpl_toolkits.mplot3d import Axes3D, from mpl_toolkits.mplot3d.art3d import Line3DCollection, import csv, import matplotlib.pyplot as plt, import numpy as np, import os`

---

### Archivo: `02_vascular_cco/generador_cco_v2.py`
- **Tipo:** Script de Python
- **Dependencias:** `from mpl_toolkits.mplot3d import Axes3D, from mpl_toolkits.mplot3d.art3d import Line3DCollection, import csv, import matplotlib.pyplot as plt, import numpy as np, import os`
- **Propósito:** Verifica si un punto está dentro del elipsoide renal...

---

### Archivo: `02_vascular_cco/generador_cco_v3.py`
- **Tipo:** Script de Python
- **Dependencias:** `from mpl_toolkits.mplot3d import Axes3D, from mpl_toolkits.mplot3d.art3d import Line3DCollection, import csv, import matplotlib.pyplot as plt, import numpy as np, import os`

---

### Archivo: `02_vascular_cco/generador_cco_v4.py`
- **Tipo:** Script de Python
- **Dependencias:** `from mpl_toolkits.mplot3d import Axes3D, import csv, os, import matplotlib.pyplot as plt, import numpy as np`

---

### Archivo: `02_vascular_cco/generador_cco_v5.py`
- **Tipo:** Script de Python
- **Dependencias:** `from mpl_toolkits.mplot3d import Axes3D, import csv, os, import matplotlib.pyplot as plt, import numpy as np`

---

### Archivo: `02_vascular_cco/generador_cco_v6.py`
- **Tipo:** Script de Python
- **Dependencias:** `from mpl_toolkits.mplot3d import Axes3D, import csv, os, import matplotlib.pyplot as plt, import numpy as np`

---

### Archivo: `02_vascular_cco/generador_cco_v7.py`
- **Tipo:** Script de Python
- **Dependencias:** `from mpl_toolkits.mplot3d import Axes3D, import csv, os, import matplotlib.pyplot as plt, import numpy as np`

---

### Archivo: `02_vascular_cco/senyal_hipoxia_para_cco.csv`
- **Tipo:** Dataset / Datos Estructurales
- **Columnas:** `ix,iy,iz,x_cm,y_cm,z_cm,prioridad`
- **Registros:** 520

---

### Archivo: `03_modelos_3d/generar_rinon_completo.py`
- **Tipo:** Script de Python
- **Dependencias:** `import bpy, import csv, import math, import numpy as np, import os`

---

### Archivo: `03_modelos_3d/solo_capsula.py`
- **Tipo:** Script de Python
- **Dependencias:** `import bpy, import math, import os`

---

### Archivo: `06_app/biokidney_app.py`
- **Tipo:** Script de Python
- **Dependencias:** `from PyQt6.QtCharts import (, from PyQt6.QtCore import Qt, QTimer, from PyQt6.QtGui import QFont, QColor, from PyQt6.QtGui import QPainter, QBrush, QPen, from PyQt6.QtWidgets import (, import numpy as np, import sys`
- **Propósito:** QMainWindow { background-color: #F4F6F8; } QTabWidget::pane { border: 1px solid #D5D8DC; background: #FFFFFF; border-radius: 8px; } QTabBar::tab { background: #D5D8DC; color: #2C3E50; padding: 10px 24px; font-size: 13px; font-weight: bold; border-top...

---

### Archivo: `06_app/dashboard_maestro_app.py`
- **Tipo:** Script de Python
- **Dependencias:** `from PyQt6.QtCore import Qt,QTimer,QRect,QPropertyAnimation,QEasingCurve,pyqtSignal, from PyQt6.QtGui import QFont,QColor,QPalette,QLinearGradient,QPainter,QBrush,QPen,QRadialGradient,QCursor, from PyQt6.QtWidgets import QApplication,QMainWindow,QWidget,QVBoxLayout,QHBoxLayout,QLabel,QPushButton,QGroupBox,QGridLayout,QFrame,QSizePolicy,QFileDialog,QMessageBox,QStatusBar,QScrollArea, from datetime import datetime, import sys,os,subprocess`

---

### Archivo: `06_app/filtracion_glomerular_gui.py`
- **Tipo:** Script de Python
- **Dependencias:** `from PyQt6.QtCharts import QChart, QChartView, QLineSeries, QValueAxis, QAreaSeries, from PyQt6.QtCore import Qt, QTimer, from PyQt6.QtGui import QFont, QColor, QPainter, QBrush, QPen, from PyQt6.QtWidgets import (, import numpy as np, import os, import sys`
- **Propósito:** QMainWindow {{ background-color: {AZUL_OSCURO}; }} QFrame#MainFrame {{ background-color: {AZUL_OSCURO}; }} QGroupBox {{ font-size: 13px; font-weight: bold; color: {AZUL_CLARO}; border: 1.5px solid #30363D; border-radius: 8px; margin-top: 15px; paddin...

---

### Archivo: `BLUEPRINT_INGENIERIA.md`
- **Tipo:** Documentación Markdown

---

### Archivo: `CONTEXTO_PROYECTO.md`
- **Tipo:** Documentación Markdown

---

### Archivo: `README.md`
- **Tipo:** Documentación Markdown

---

### Archivo: `analizador_proyecto_biokidney.py`
- **Tipo:** Script de Python
- **Dependencias:** `from pathlib import Path, from typing import List, Dict, Optional, import datetime, import logging, import os`
- **Propósito:** Genera una representación visual de la estructura de directorios....

---

### Archivo: `biokidney/__init__.py`
- **Tipo:** Script de Python

---

### Archivo: `biokidney/aggregator.py`
- **Tipo:** Script de Python
- **Dependencias:** `from biokidney.core.config import cfg_physio, cfg_vasc, from biokidney.experts.fluids import FluidDynamicsExpert, from biokidney.experts.vascular import VascularExpert, import numpy as np`
- **Propósito:** Orquestador Principal BioKidney-AI (Patrón Agregación) ------------------------------------------------------ Este módulo actúa como el punto de integración central del framework, agregando a los expertos especializados en cada dominio renal. Agregad...

---

### Archivo: `biokidney/core/__init__.py`
- **Tipo:** Script de Python

---

### Archivo: `biokidney/core/config.py`
- **Tipo:** Script de Python
- **Dependencias:** `from dataclasses import dataclass, field, import numpy as np`
- **Propósito:** Módulo de Configuración Centralizada BioKidney-AI ------------------------------------------------ Este módulo define las constantes fisiológicas y parámetros de simulación para asegurar la consistencia en todo el ecosistema....

---

### Archivo: `biokidney/experts/__init__.py`
- **Tipo:** Script de Python

---

### Archivo: `biokidney/experts/cellular.py`
- **Tipo:** Script de Python
- **Dependencias:** `from biokidney.core.config import cfg_physio, import numpy as np`
- **Propósito:** Calcula el porcentaje de volumen en estado de hipoxia....

---

### Archivo: `biokidney/experts/fluids.py`
- **Tipo:** Script de Python
- **Dependencias:** `from biokidney.core.config import cfg_physio, cfg_sim, from typing import Tuple, List, import numpy as np`
- **Propósito:** Modelo Deen para la presión oncótica a lo largo del capilar....

---

### Archivo: `biokidney/experts/vascular.py`
- **Tipo:** Script de Python
- **Dependencias:** `from biokidney.core.config import cfg_vasc, cfg_sim, from pathlib import Path, from typing import List, Tuple, Optional, import numpy as np`
- **Propósito:** Verifica si un punto está dentro del elipsoide renal....

---

### Archivo: `biokidney/utils/__init__.py`
- **Tipo:** Script de Python

---

### Archivo: `biokidney_architect.py`
- **Tipo:** Script de Python
- **Dependencias:** `from datetime import datetime, from pathlib import Path, import logging, import os, import sys`
- **Propósito:** Verifica si la base de conocimiento existe....

---

### Archivo: `requirements.txt`

---

### Archivo: `web_app/backend/database/models.py`
- **Tipo:** Script de Python
- **Dependencias:** `from datetime import datetime, from sqlalchemy import Column, Integer, String, Float, DateTime, JSON, ForeignKey, create_engine, from sqlalchemy.ext.declarative import declarative_base, from sqlalchemy.orm import sessionmaker, relationship`
- **Propósito:** Historial de simulaciones clínicas. Esencial para trazabilidad médica y revisión de inversionistas. Métricas de performance del sistema. Trazas de auditoría para el módulo de administración....

---

### Archivo: `web_app/backend/main.py`
- **Tipo:** Script de Python
- **Dependencias:** `from fastapi import FastAPI, HTTPException, Request, from fastapi.middleware.cors import CORSMiddleware, from fastapi.responses import HTMLResponse, from fastapi.staticfiles import StaticFiles, from pydantic import BaseModel, from typing import List, Dict, Any, from web_app.backend.database.models import init_db, from web_app.backend.services.simulation_service import SimulationService...`
- **Propósito:** Sirve la SPA del dashboard....

---

### Archivo: `web_app/backend/services/simulation_service.py`
- **Tipo:** Script de Python
- **Dependencias:** `from biokidney.aggregator import BioKidneyEngine, from biokidney.core.config import cfg_physio, from typing import Dict, Any, from web_app.backend.database.models import SessionLocal, SimulationRecord, from web_app.backend.utils.logger import bk_logger, import time`
- **Propósito:** Capa de servicio que orquesta las simulaciones del core biokidney y persiste los resultados para auditoría....

---

### Archivo: `web_app/backend/utils/logger.py`
- **Tipo:** Script de Python
- **Dependencias:** `from loguru import logger, from pathlib import Path, import logging, import sys`
- **Propósito:** Configura un sistema de logging profesional con Loguru. - Consola con colores para desarrollo. - Archivo rotativo para auditoría médica. - Estructura JSON para fácil integración con ELK/CloudWatch....

---

## 3. Inventario de Activos Binarios

#### Modelos Blender (.blend)
- **Total:** 6
- **Muestras:** `V7_recto_untitled.blend, arbol_vascular_cco.blend, rinon_biokidney_ai.blend, rinon_capsula.blend, v7_curvountitled.blend...`
#### Imágenes/Visualizaciones (.png)
- **Total:** 21
- **Muestras:** `dashboard_maestro.png, pareto_coswift.png, render_biokidney_v7_definitivo.png, simulador_estres_mecanico_dECM.png, untitled.png...`
#### Informes (.pdf)
- **Total:** 8
- **Muestras:** `Bitacora_BioKidney_AI_2026.pdf, reporte_validacion_completo.pdf, resumen_ejecutivo_una_pagina.pdf, scorecard_biokidney_ai_2026.pdf, simulador_filtracion_glomerular_reporte.pdf...`
#### Videos (.mp4)
- **Total:** 2
- **Muestras:** `Arbol vascular CCO v7.mp4, Vascular_Tree_Video_An_abstract_glowing_animation_displays_two_YtuUY7DV.mp4`
#### Mallas 3D (.obj)
- **Total:** 3
- **Muestras:** `arbol_vascular_cco.obj, arbol_vascular_cco_v7.obj, arbol_vascular_cco_v7_curvo.obj`
#### Documentos Word (.docx)
- **Total:** 5
- **Muestras:** `BioKidney_AI_Informe_Desarrollo_dashboard_maestro.docx, BiokidneyAI_Informe_Avance_Mar2026.docx, Informe_Modulo09_BioKidney_AI_2026.docx, informe_biokidney_ai_2026_Simulador_de_Reabsorción_Tubular.docx, informe_biokidney_coswift.docx`