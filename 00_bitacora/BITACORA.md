# Bitácora de Proyecto: Bio-Kidney AI
**Investigador:** [Tu nombre]
**Inicio del proyecto:** Enero 2026
**Objetivo:** Desarrollar el modelo computacional completo de un riñón bioimpreso funcional, desde la síntesis vascular hasta la validación hemodinámica, como base científica para el primer riñón artificial bioimpreso viable.

---

## ENTRADA 001 — Enero 2026
**Estado:** Investigación inicial

### Lo que hice
- Investigué el estado del arte en bioimpresión renal
- Estudié iPSCs, factores Yamanaka, uso de Wnt para redireccionamiento celular
- Aprendí anatomía y fisiología renal aplicada al proyecto
- Comencé a aprender inglés de forma empírica (nivel intermedio alcanzado)

### Conceptos clave aprendidos
- iPSCs como fuente celular autóloga
- Técnica Co-SWIFT para bioimpresión coaxial
- Ley de Murray para geometría vascular
- OBBs (Organ Building Blocks) de organoides renales

### Decisiones tomadas
- Usar Ubuntu + Python como entorno principal de desarrollo
- Enfoque in silico antes de cualquier trabajo húmedo de laboratorio

---

## ENTRADA 002 — Febrero 2026
**Estado:** Primeros simuladores funcionales

### Lo que hice
- Desarrollé 3 simuladores en Python independientes:
  1. `simulador_wss_biokidney.py` — Hemodinámica (Wall Shear Stress)
  2. `simulador_ipsc_biokidney.py` — Cinética de maduración celular
  3. `simulador_swift_biokidney.py` — Optimización de extrusión

### Resultados obtenidos (v1)
- Presión óptima de extrusión: 30 kPa
- Ventana de bioimpresión: día 21–30 de cultivo
- Radio mínimo de vaso funcional: 50 µm

### Problemas identificados en v1
- WSS v1 usaba zona de seguridad incorrecta (10–70 dyn/cm²)
  → El rango fisiológico renal real es 1–10 dyn/cm²
- Factor de conversión incorrecto inflaba los valores de WSS artificialmente
- Los simuladores no estaban conectados entre sí

---

## ENTRADA 003 — Marzo 2026
**Estado:** Dashboard integrado v2.0

### Lo que hice
- Corrección del error de WSS: rango fisiológico real 1–10 dyn/cm²
- Integración de los 3 simuladores en un dashboard HTML interactivo
- Pipeline conectado: iPSC → SWIFT → WSS → Estado global del órgano
- Inicio del modelo 3D en Blender (forma exterior del riñón)
- Diseñé arquitectura del generador vascular CCO con Ley de Murray

### Archivos generados
- `dashboard_v2.html` — Dashboard integrado interactivo
- `untitled.blend` — Modelo 3D inicial (esfera + boolean difference para hilio)

### Correcciones científicas aplicadas
- WSS fisiológico renal en arteriolas: 1–10 dyn/cm² (no 10–70)
- Estrés mecánico límite para viabilidad celular: <150 Pa (confirmado)
- Zona de bioimpresión: día 21–30 validada con modelo cinético corregido

### Próximo paso
- Generador vascular CCO (Constructive Constrained Optimization)
- Implementar Ley de Murray: r_padre³ = r_hijo1³ + r_hijo2³
- Añadir número de Reynolds para validar régimen laminar

---

## ENTRADA 004 — [FECHA]
**Estado:** 

### Lo que hice


### Resultados


### Problemas / dudas


### Próximo paso


---

## IDEAS PENDIENTES
- [ ] Generador CCO con Ley de Murray en Python
- [ ] Integrar CFD con OpenFOAM
- [ ] Modelo de difusión de oxígeno (gradiente en 200M células/ml)
- [ ] CNN para optimización del árbol vascular
- [ ] Exportar geometría vascular a Blender para render
- [ ] Preparar pitch para inversores

## CUELLOS DE BOTELLA IDENTIFICADOS
1. Micro-vascularización < 200 µm (angiogénesis no controlable aún)
2. Densidad celular 200M/ml → muerte por hipoxia sin perfusión inmediata
3. Maduración funcional completa del nefrón (túbulos + glomérulo)
4. Anastomosis sin trombosis aguda

## REFERENCIAS CLAVE
- Lewis Lab (Harvard) — SWIFT en tejido cardiaco
- Takahashi et al. — iPSCs y organoides renales
- Murray CW (1926) — Ley de Murray, geometría vascular óptima
- Finet et al. — Regla de bifurcación vascular


---
## Entrada: 21 Marzo 2026

### App de escritorio v2.0 — Bio-Kidney AI Simulador
- Desarrollada en PyQt6
- Tres módulos: WSS, iPSC, SWIFT
- Archivo: 06_app/biokidney_app.py
- Acceso directo: BioKidney.desktop en el escritorio
- Parámetros validados con valores ÓPTIMO por defecto
- Simuladores originales preservados en: 01_simuladores/v1_originales/

---
## Entrada: 21 Marzo 2026 — Generador Vascular CCO

### Resultados de validación
- Segmentos generados : 161
- Bifurcaciones       : 80
- Cumplimiento Murray : 100.0%
- Flujo laminar       : 95.0% (153/161)
- Radio mínimo        : 32.5 µm
- Radio máximo        : 500.0 µm
- Reynolds promedio   : 479.0

### Archivos generados
- 02_vascular_cco/generador_cco.py
- 02_vascular_cco/arbol_vascular_cco.csv
- 02_vascular_cco/arbol_vascular_cco.png

### Pendiente
- Reducir Reynolds máximo en segmentos raíz
- Exportar geometría a Blender (.obj)
- Aumentar resolución a 500+ segmentos

---
## Entrada: 21 Marzo 2026 — Modelo 3D Riñón Completo

### Archivos generados
- 03_modelos_3d/generar_rinon_completo.py
- 03_modelos_3d/rinon_biokidney_ai.blend
- 03_modelos_3d/render_biokidney_v1.png

### Componentes del modelo
- Cápsula renal: elipsoide con proporción anatómica (1.0 x 0.6 x 0.55)
- Hilio renal: boolean difference aplicado
- Árbol vascular CCO: 161 segmentos integrados
- Materiales: cápsula translúcida (alpha 0.25) + vasos rojos
- Iluminación: área principal 800W + relleno 200W
- Render: EEVEE, 2 segundos, fondo azul oscuro

### Pendientes modelo 3D
- Contener árbol vascular dentro de la cápsula (ajuste dominio CCO)
- Aumentar lados cilindros de 8 a 16 para mejor resolución
- Cambiar render a Cycles para mayor realismo
- Agregar material subsurface scattering en cápsula

### Pendiente general
- Versión web (Dash + Render.com) para acceso móvil sin WiFi

---
## Entrada: 21 Marzo 2026 — Generador CCO v2

### Resultados validación v2
- Nodos totales        : 363
- Segmentos con spline : 362
- Terminales (glomér.) : 120
- Radio mínimo         : 111.6 µm
- Radio máximo         : 500.0 µm
- Nodos dentro riñón   : 100.0%
- Cumplimiento Murray  : 87.8%

### Mejoras sobre v1
- Dominio elipsoide anatómico (5.5 x 3.0 x 2.5 cm)
- Curvas spline cúbico en cada segmento
- 120 glomérulos distribuidos uniformemente
- Hilio renal como punto de origen

### Pendiente v3
- Jerarquía anatómica de 6 niveles
- Radio mínimo hasta 15-20 µm
- Cumplimiento Murray >95%

---
## Entrada: 21 Marzo 2026 — Generador CCO v4

### Resultados validación v4
- Nodos arteriales    : 127
- Nodos venosos       : 127
- Nodos colectores    : 51
- Total nodos         : 305
- Nodos dentro riñón  : 100.0%
- Radio mínimo        : 81.4 µm
- Radio máximo        : 600.0 µm
- Murray OK           : 100.0% (126/126 bifurcaciones)
- Murray violaciones  : 0

### Mejoras sobre v3
- Murray estricto con verificación por assertion
- Radios propagados exclusivamente desde Murray
- Jerarquía anatómica de 6 niveles con direcciones por nivel
- Arterial y venoso perfectamente simétricos

### Pendientes v5
- Bajar radio mínimo a 15-20 µm (agregar 2-3 niveles)
- Mejorar distribución hacia polos del riñón
- Exportar v4 a Blender para render completo

---
## Entrada: 22 Marzo 2026 — Generador CCO v5

### Resultados validación v5
- Nodos arteriales    : 511
- Nodos venosos       : 511
- Nodos colectores    : 59
- Total nodos         : 1081
- Segmentos totales   : 1078
- Nodos dentro riñón  : 100.0%
- Radio mínimo        : 41.2 µm
- Radio máximo        : 600.0 µm
- Radio promedio      : 110.8 µm
- Murray OK           : 510/510 (100.0%)
- Murray violaciones  : 0

### Distribución por nivel (arterial)
- Nivel 0:   1 nodo  | 500.0 µm
- Nivel 1:   2 nodos | 383–410 µm
- Nivel 2:   4 nodos | 286–327 µm
- Nivel 3:   8 nodos | 220–280 µm
- Nivel 4:  16 nodos | 156–233 µm
- Nivel 5:  32 nodos | 112–191 µm
- Nivel 6:  64 nodos | 75–166 µm
- Nivel 7: 128 nodos | 53–142 µm
- Nivel 8: 256 nodos | 41–125 µm

### Evolución del proyecto CCO
| Versión | Murray | R_min   | Nodos | Sistemas      |
|---------|--------|---------|-------|---------------|
| v1      | 100%   | 32.5 µm |   161 | Solo arterial |
| v2      | 87.8%  | 111.6 µm|   363 | Solo arterial |
| v3      | 57.4%  | 103.9 µm|   404 | 3 sistemas    |
| v4      | 100%   | 81.4 µm |   305 | 3 sistemas    |
| v5      | 100%   | 41.2 µm | 1081  | 3 sistemas    |

### Objetivo siguiente: v6
- Llegar a R_min = 10-15 µm (arteriolas aferentes reales)
- Mejorar distribución polar del árbol
- Estrategia: cambiar bifurcación binaria pura por
  crecimiento adaptativo hacia zonas sin cobertura vascular
- Estimado: 3000-5000 nodos, nivel 10-12

---
## Entrada: 22 Marzo 2026 — Generador CCO v6

### Resultados validación v6
- Nodos arteriales    : 203
- Nodos venosos       : 215
- Nodos colectores    : 73
- Total nodos         : 491
- Segmentos totales   : 488
- Nodos dentro riñón  : 100.0%
- Radio mínimo        : 25.7 µm
- Radio máximo        : 600.0 µm
- Cobertura vascular  : 100.0% (400/400 glomérulos)
- Murray OK           : 208/208 (100.0%)
- Niveles alcanzados  : 12

### Hito crítico
Primera versión con cobertura vascular completa del dominio.
El algoritmo adaptativo creció hacia zonas sin irrigación
replicando el comportamiento de la angiogénesis real.

### Evolución completa CCO
| Versión | Murray | R_min   | Nodos | Cobertura |
|---------|--------|---------|-------|-----------|
| v1      | 100%   | 32.5 µm |   161 | N/A       |
| v2      | 87.8%  | 111.6 µm|   363 | N/A       |
| v3      | 57.4%  | 103.9 µm|   404 | N/A       |
| v4      | 100%   | 81.4 µm |   305 | N/A       |
| v5      | 100%   | 41.2 µm | 1081  | N/A       |
| v6      | 100%   | 25.7 µm |  491  | 100%      |

### Pendiente v7
- Bajar R_min a 10-15 µm
- Aumentar N_SEMILLAS_DEMANDA a 800
- Aumentar MAX_ITERACIONES a 4000

---
## Entrada: 24 Marzo 2026 — Render Definitivo Bio-Kidney v7

### Modelo 3D completado
- Cápsula renal: elipsoide anatómico con hilio boolean
- Árbol vascular: CCO v7 con curvas Bezier reales
- Tres sistemas: arterial (rojo) + venoso + colector
- Motor: EEVEE, 2.43 segundos de render
- Archivo: 03_modelos_3d/render_biokidney_v7_definitivo.png
- Archivo Blender: 03_modelos_3d/rinon_capsula.blend

### Estado actual del proyecto
- Simuladores Python: COMPLETADOS (WSS, iPSC, SWIFT)
- App de escritorio PyQt6: COMPLETADA
- Generador CCO v7: COMPLETADO (Murray 100%, cobertura 100%)
- Modelo 3D Blender: COMPLETADO
- @ComunidadRenal: 2 posts publicados, Post 3 programado

---
## Entrada: 25 Marzo 2026 — Simulador Difusión O₂

### OBJETIVO ALCANZADO — 100% tejido oxigenado
- PO₂ mínima    : 20.000 mmHg (umbral crítico: 1.0 mmHg)
- PO₂ media     : 28.323 mmHg (rango fisiológico normal)
- PO₂ máxima    : 40.000 mmHg
- Vóxeles hipóxicos: 0 de 24 (0.00%)
- Densidad celular: 200M células/mL

### Física implementada
- Ley de Fick 3D con consumo Michaelis-Menten
- Solver SOR (omega=1.6), tolerancia 1e-4 mmHg
- Grid: 60x60x40 vóxeles, dominio 95x54x45 mm

### Archivos generados
- 01_simuladores/simulador_oxigeno_biokidney.py
- 01_simuladores/resultado_oxigeno_biokidney.png
- 02_vascular_cco/senyal_hipoxia_para_cco.csv
- Informe PDF: informe_sesion_biokidney_marzo2026.pdf

---
## Entrada: 25 Marzo 2026 — Simulador Diferenciación iPSC

### Resultados — 3 protocolos independientes
- Podocitos        : ventana días 15.9-30 | pureza 100% | APTO
- Tub. Proximales  : ventana días 16.9-30 | pureza 100% | APTO
- Endoteliales iECs: ventana días 16.9-30 | pureza 100% | APTO
- Riesgo teratoma día 21: 0.0000% — BAJO en todos
- Protocolos ≥95% pureza: 3/3
- Base científica: Takasato 2015 + Freedman 2015

### Convergencia del pipeline
La ventana días 16-30 confirma el simulador de cinética iPSC v1.
El riesgo de teratoma cero valida el protocolo de purificación.

### Archivos generados
- 01_simuladores/simulador_diferenciacion_ipsc.py
- 01_simuladores/simulador_diferenciacion_ipsc_output.png
- BiokidneyAI_Informe_Avance_Mar2026.docx

---
## Entrada: 25 Marzo 2026 — Simulador Estrés Mecánico dECM

### Resultados — 4 materiales validados
- GelMA 7%      : ÓPTIMO | ventana 0-47 kPa
- Alginato 1.5% : ÓPTIMO | ventana 0-28 kPa
- NICE Bioink   : ÓPTIMO | ventana 0-56 kPa (MEJOR)
- dECM Espinaca : ÓPTIMO | ventana 0-36 kPa
- P_ext SWIFT 30 kPa: dentro de ventana segura en TODOS
- Shear canal Q=0.5nL/s: 0.004 Pa (límite 150 Pa — factor 37,500x)
- Estado global: ÓPTIMO MECÁNICO dECM

### Física validada
- Elasticidad lineal (Hooke)
- Viscoelástico Kelvin-Voigt
- Hagen-Poiseuille en canales vasculares
- Von Mises criterio de fallo (Lamé pared gruesa)
- Deformación radial de canales Co-SWIFT

### Archivos generados
- 01_simuladores/simulador_estres_mecanico_dECM.py
- 01_simuladores/simulador_estres_mecanico_dECM.png
- Informe_Modulo09_BioKidney_AI_2026.docx

---
## Entrada: 25 Marzo 2026 — Optimizador Co-SWIFT PSO

### Resultados PSO Multi-Objetivo
- Algoritmo: PSO 80 partículas x 150 iteraciones
- Evaluaciones totales: 12,080
- Frente de Pareto: 100 soluciones no-dominadas
- Convergencia: iteración 21

### Protocolo óptimo encontrado
- Presión extrusión : 5.6 kPa
- Velocidad impresión: 21.5 mm/s
- Temperatura bioink : 36.5°C
- Diámetro boquilla  : 283 µm
- Tiempo entre capas : 43 s
- Shear stress       : 39.6 Pa (límite 150 Pa)
- Canal vascular     : 198 µm (objetivo 200 µm)
- Viabilidad celular : 98.0% (objetivo >85%)

### Archivos generados
- 01_simuladores/optimizador_coswift.py
- resultados/pareto_coswift.png
- informe_biokidney_coswift.docx

---
## Entrada: 26 Marzo 2026 — Simulador Filtración Glomerular

### HITO ALCANZADO — Pipeline al 80%
- TFG calculado     : 82.02 mL/min
- TFG objetivo      : >= 60 mL/min (umbral diálisis)
- Presión neta dP   : 20.01 mmHg
- Fracción filtración: 13.1%
- Estado global     : FUNCIONAL
- TFG nativo ref.   : 62.5 mL/min → Bio-Kidney supera en 31%

### Física implementada
- Ecuación de Starling (4 presiones)
- Modelo Deen-Robertson-Brenner (concentración oncótica)
- Kf calibrado: 3.7 mL/min/mmHg
- 1000 glomérulos escalados desde árbol CCO v7

### Pendientes de este módulo
- GUI PyQt6 interactiva con sliders en tiempo real
- Corrección error PDF (imagen demasiado grande)

### Próximo módulo
- Reabsorción Tubular (+10%) → pipeline al 90%
- Submodulos: Túbulo proximal, Asa de Henle, Túbulo distal, Colector
- Salida esperada: ~1.5 mL/min orina funcional

### Archivos generados
- 01_simuladores/simulador_filtracion_glomerular.py
- 01_simuladores/simulador_filtracion_glomerular.png
- BioKidney_AI_2026_Informe_Sesion_TFG.pdf

---
## Entrada: 26 Marzo 2026 — Simulador Reabsorción Tubular

### HITO ALCANZADO — Pipeline al 90%
- Estado global     : ÓPTIMO — 6/6 criterios
- Orina producida   : 1.520 mL/min = 2.19 L/día
- Reabsorción total : 80.48 mL/min (98.1% del filtrado)
- Función tubular   : 78.2% vs riñón nativo
- Osmolaridad orina : 1200 mOsm/kg (máxima concentración)
- Glucosa orina     : 0.18 mg/dL (reabsorción casi completa)

### Segmentos tubulares validados
- Túbulo Proximal   : reabsorbe 67% (27.06 mL/min salida)
- Asa Descendente   : concentra hasta 1200 mOsm/kg
- Asa Ascendente    : Na⁺/K⁺/Cl⁻ sin agua, 100 mOsm/kg
- Túbulo Distal     : ajuste fino aldosterona
- Túbulo Colector   : concentración final ADH

### Transportadores validados
- SGLT2 (TP)  : 71.6% saturación
- NHE3 (TP)   : 90.3% saturación
- AQP1 (Asc)  : 90.9% saturación
- NKCC2 (AHA) : 94.3% saturación
- ENaC (TD)   : 96.3% saturación
- AQP2 (TC)   : 95.3% saturación

### Archivos generados
- 01_simuladores/simulador_reabsorcion_tubular.py
- 01_simuladores/reabsorcion_tubular_dashboard.png
- 01_simuladores/reabsorcion_tubular_informe.pdf
- informe_biokidney_ai_2026_Simulador_de_Reabsorcion_Tubular.docx

---
## Entrada: 26 Marzo 2026 — PIPELINE 100% COMPLETO

### HITO HISTÓRICO
Bio-Kidney AI 2026 — Validación In Silico Completa
12/12 módulos ÓPTIMOS
7/12 KPIs superan el riñón nativo humano

### Estado final del pipeline
- CCO v7 Vascularización  : ✓ Murray 100%, 1448 seg
- Difusión O₂             : ✓ 100% oxigenado, 0% hipóxico
- Diferenciación iPSC     : ✓ 3/3 linajes, pureza 100%
- Bioimpresión Co-SWIFT   : ✓ 98% viabilidad, Pareto 100
- Filtración Glomerular   : ✓ TFG 82 mL/min
- Reabsorción Tubular     : ✓ 2.19 L/día, 6/6 criterios
- Dashboard Maestro       : ✓ Integración completa

### Documentos finales generados
- 07_presentacion_final/dashboard_maestro.png
- 07_presentacion_final/scorecard_biokidney_ai_2026.pdf
- 07_presentacion_final/reporte_validacion_completo.pdf
- 07_presentacion_final/resumen_ejecutivo_una_pagina.pdf

### Próxima fase
- Enviar resumen ejecutivo al Dr. Hincapié (UDEA)
- Preparar preprint para bioRxiv
- Contactar Harvard Wyss Institute
- Contactar Oxford IBME
- Fase experimental con laboratorio aliado

---
---
## Entrada: 26 Marzo 2026 — Ajustes Finales y GUI Interactiva

### Hitos técnicos alcanzados
- **Bug Fix PDF:** Corregido el tamaño de imagen en el reporte del Simulador de Filtración Glomerular (ajustado a 22x13 cm para evitar saltos de página).
- **Nueva GUI Filtración:** Desarrollada aplicación PyQt6 `06_app/filtracion_glomerular_gui.py` con sliders interactivos para Pgc, Pbs, π y Kf.
- **Sincronización 100%:** Actualizados los pies de página y estados internos en los simuladores de Filtración y Reabsorción Tubular para reflejar el estado final del pipeline.

### Resultados de validación final
- Estado global: **RIÑÓN FUNCIONAL ✓**
- Pipeline in silico: **100% COMPLETADO**
- Consistencia documental: Verificada en todos los módulos.

### Próxima fase (Inmediata)
- Ejecutar `analizador_proyecto_biokidney.py` para generar el `CONTEXTO_PROYECTO.md` final.
- Preparar el envío formal al Dr. Hincapié (UDEA).
- Iniciar redacción de abstract para bioRxiv.

---
## Entrada: 27-28 Marzo 2026 — Identidad Científica y Primera Publicación
### Hitos alcanzados
- **ORCID creado:** 0009-0005-3933-5072 (https://orcid.org/0009-0005-3933-5072)
- **GitHub público:** https://github.com/VirtusSapiens/Bio-Kidney-AI-2026
- **4 tags de versión:** v1.0.0 → v1.1.0 → v1.2.0 → v2.0.0
- **Preprint v1 escrito:** Abstract, Introduction, Methods, Results, Discussion, Conclusion, 10 referencias
- **bioRxiv:** Sometido como BIORXIV/2026/714957 — rechazado por falta de afiliación institucional (no por contenido)
- **Cuenta bioRxiv creada:** david.moreno.159cm@gmail.com
- **LinkedIn actualizado:** Publicación agregada con DOI
### Archivos generados
- 00_bitacora/preprint_biokidney_2026.md
- 00_bitacora/preprint_biokidney_2026.pdf

---
## Entrada: 29-30 Marzo 2026 — CCO v8 y Preprint v2
### Hitos técnicos
- **CCO v8 implementado:** 1,902 segmentos vasculares, 915 bifurcaciones, 100% Murray
- **TFG mejorada:** 115.2 mL/min (rango normal adulto sano, +40% vs v7)
- **Modelo Poiseuille calibrado:** Dos pasos, presiones terminales 58.6 +/- 13.4 mmHg
- **Distribución Beta(3,1.2):** 63% demanda glomerular hacia cortex
- **Blender v8:** Modelo 3D renderizado (sab_4_abril_1200_am_v8_1.blend)
- **Preprint v2 EN/ES:** Titulo mejorado, CCO v8 integrado, material suplementario
- **Segunda sumisión bioRxiv:** BIORXIV/2026/715287 — rechazado por afiliación institucional
### Archivos generados
- 02_vascular_cco/generador_cco_v8.py
- 02_vascular_cco/arbol_vascular_cco_v8.csv
- 02_vascular_cco/arbol_vascular_cco_v8.png
- 02_vascular_cco/importar_blender_v8.py
- 03_modelos_3d/sab_4_abril_1200_am_v8_1.blend
- 00_bitacora/preprint_biokidney_2026_EN.pdf
- 00_bitacora/preprint_biokidney_2026_ES.md
- 00_bitacora/supplementary_material_v8.md
- biokidney/simulation/fractal_vascularizer.py
- biokidney/simulation/v8_to_fractal.py

---
## Entrada: 11 Abril 2026 — Zenodo y Correo Dr. Hincapié
### Hitos alcanzados
- **Zenodo publicado:** DOI 10.5281/zenodo.19508077
  URL: https://zenodo.org/records/19508077
  Preprint abierto, indexado en OpenAIRE
- **LinkedIn actualizado:** Publicación con DOI de Zenodo
- **Correo Dr. Hincapié programado:** Lunes 14 abril 7:00 AM
  Asunto: Propuesta de colaboracion tecnica e institucional: Framework Bio-Kidney AI 2026 (TFG 115.2 mL/min)
- **Rama fractal creada:** feature/cco-v8-fractal
- **Git push exitoso:** Todos los archivos v8 en GitHub
### Segunda fase definida
- Pasar de validacion matematica a factibilidad tecnica de fabricacion
- Especificacion de biotintas por modulo
- Protocolo de impresion 3D capa por capa
- Integracion angiogenesis post-impresion
- Manual de implementacion para laboratorios
- Automatizacion: Python genera coordenadas para 19,000+ glomérulos en Blender

---
## Entrada: 27-28 Marzo 2026 — Identidad Cientifica y Primera Publicacion
### Hitos alcanzados
- ORCID creado: 0009-0005-3933-5072 (https://orcid.org/0009-0005-3933-5072)
- GitHub publico: https://github.com/VirtusSapiens/Bio-Kidney-AI-2026
- 4 tags de version: v1.0.0 a v2.0.0
- Preprint v1 escrito: Abstract, Introduction, Methods, Results, Discussion, Conclusion, 10 referencias
- bioRxiv: Sometido como BIORXIV/2026/714957 — rechazado por falta de afiliacion institucional
- Cuenta bioRxiv creada: david.moreno.159cm@gmail.com
- LinkedIn actualizado: Publicacion agregada con DOI
### Archivos generados
- 00_bitacora/preprint_biokidney_2026.md
- 00_bitacora/preprint_biokidney_2026.pdf

---
## Entrada: 29-30 Marzo 2026 — CCO v8 y Preprint v2
### Hitos tecnicos
- CCO v8 implementado: 1,902 segmentos vasculares, 915 bifurcaciones, 100% Murray
- TFG mejorada: 115.2 mL/min (rango normal adulto sano, +40% vs v7)
- Modelo Poiseuille calibrado de dos pasos: presiones terminales 58.6 +/- 13.4 mmHg
- Distribucion Beta(3,1.2): 63% demanda glomerular hacia cortex
- Blender v8: Modelo 3D renderizado (sab_4_abril_1200_am_v8_1.blend)
- Preprint v2 EN/ES: Titulo mejorado, CCO v8 integrado, material suplementario
- Segunda sumision bioRxiv: BIORXIV/2026/715287 — rechazado por afiliacion institucional
- Video para Instagram generado: video_para_instagram.mp4
### Archivos generados
- 02_vascular_cco/generador_cco_v8.py
- 02_vascular_cco/arbol_vascular_cco_v8.csv
- 02_vascular_cco/arbol_vascular_cco_v8.png
- 02_vascular_cco/importar_blender_v8.py
- 03_modelos_3d/sab_4_abril_1200_am_v8_1.blend
- 00_bitacora/preprint_biokidney_2026_EN.pdf
- 00_bitacora/preprint_biokidney_2026_ES.md
- 00_bitacora/supplementary_material_v8.md
- biokidney/simulation/fractal_vascularizer.py
- biokidney/simulation/v8_to_fractal.py
- renal_data_v8_fractal.json

---
## Entrada: 11 Abril 2026 — Zenodo, Bitacora y Correo Dr. Hincapie
### Hitos alcanzados
- Zenodo publicado: DOI 10.5281/zenodo.19508077
  URL: https://zenodo.org/records/19508077
  Preprint abierto, indexado en OpenAIRE
- LinkedIn actualizado: Publicacion con DOI de Zenodo
- Correo Dr. Hincapie programado: Lunes 14 abril 7:00 AM
  Asunto: Propuesta de colaboracion tecnica e institucional: Framework Bio-Kidney AI 2026 (TFG 115.2 mL/min)
- Rama fractal creada: feature/cco-v8-fractal
- Git push exitoso: Todos los archivos v8 en GitHub
- Bitacora actualizada con entradas 27 marzo al 11 abril
### Segunda fase definida
- Pasar de validacion matematica a factibilidad tecnica de fabricacion
- Especificacion de biotintas por modulo
- Protocolo de impresion 3D capa por capa
- Integracion angiogenesis post-impresion
- Manual de implementacion para laboratorios
- Automatizacion: Python genera coordenadas para 19,000+ glomérulos en Blender
### Pendientes
- Respuesta Dr. Hincapie (afiliacion GIB-UDEA para resometer a bioRxiv)
- ResearchGate: abrir con DOI de Zenodo como evidencia
- CCO v8 fractal: implementar L-systems + Murray + CFD
- Deployment publico SPA web (Railway/Render)
- Implementar KD-Trees en CCO para escalabilidad mayor a 5000 nodos
