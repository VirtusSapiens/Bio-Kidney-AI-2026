# Bitácora de Proyecto: Bio-Kidney AI
**Investigador:** Carlos David Moreno Cáceres · VirtusSapiens
**Inicio del proyecto:** Enero 2026
**Objetivo:** Desarrollar el modelo computacional completo de un riñón bioimpreso funcional, desde la síntesis vascular hasta la validación hemodinámica, como base científica para el primer riñón artificial bioimpreso viable.

> **Nota de lectura (2026-09-21).** Esta bitácora es un registro cronológico
> y no se reescribe. El objetivo de arriba es el original de enero de 2026.
> El alcance vigente del proyecto y la lista de afirmaciones retiradas están
> en `MAESTRO_contexto_canonico_BioKidney_v2.md` (§0). Las entradas
> anteriores al 2026-08-22 pueden contener cifras y encuadres retirados; las
> líneas con TFG decimal identificadas en la auditoría del 2026-09-21 llevan
> el marcador [RETIRADO: ver MAESTRO §0].
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
- CCO v7 Vascularización  : [OK] Murray 100%, 1448 seg
- Difusión O₂             : [OK] 100% oxigenado, 0% hipóxico
- Diferenciación iPSC     : [OK] 3/3 linajes, pureza 100%
- Bioimpresión Co-SWIFT   : [OK] 98% viabilidad, Pareto 100
- Filtración Glomerular   : [OK] TFG 82 mL/min
- Reabsorción Tubular     : [OK] 2.19 L/día, 6/6 criterios
- Dashboard Maestro       : [OK] Integración completa

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
- Estado global: **RIÑÓN FUNCIONAL [OK]**
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
- **TFG mejorada:** 115.2 mL/min (rango normal adulto sano, +40% vs v7) [RETIRADO: ver MAESTRO §0]
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
  Asunto: Propuesta de colaboracion tecnica e institucional: Framework Bio-Kidney AI 2026 (TFG 115.2 mL/min) [RETIRADO: ver MAESTRO §0]
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
- TFG mejorada: 115.2 mL/min (rango normal adulto sano, +40% vs v7) [RETIRADO: ver MAESTRO §0]
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
  Asunto: Propuesta de colaboracion tecnica e institucional: Framework Bio-Kidney AI 2026 (TFG 115.2 mL/min) [RETIRADO: ver MAESTRO §0]
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

---

## ENTRADA 005 — 12 de junio de 2026 — Arranque del gemelo digital (Capa 0)
**Estado:** Capa 0 generada y visualizada. Pendiente: calibración fina del seno + Capa 1 (nefronas).
**Sesión:** Arranque del gemelo digital — Capa 0 (Dominio renal).

### Encuadre estratégico
El gemelo digital del riñón es el **primer paso obligatorio antes de imprimir biogel**: el plano y el banco de pruebas de todo lo material. No se imprime para ver si funciona; se predice computacionalmente que funcionará, y *entonces* se imprime. Este orden es lo que distingue un proyecto de ingeniería serio de un ensayo por prueba y error.

El árbol vascular no es una pieza más del proyecto: es la **columna estructural** del gemelo digital. Sin arquitectura vascular correcta, los seis módulos biofísicos previos corren sobre una geometría que no se sostiene.

### Ideas propias a conservar (no perder)
- **Analogía geodésica.** Así como lat/long/elevación ubican un punto en la Tierra (ovoide), un sistema de coordenadas renales (profundidad corteza→médula, eje polar, segmento de pirámide) puede indexar cada ramificación del árbol. El eje corteza→médula es la "elevación". → Se materializa como los *puntos de atracción* del generador vascular.
- **La incertidumbre como motor, no como defecto.** El camino autodirigido da miedo precisamente porque no tiene plantilla; pero esa misma ausencia de plantilla es lo que permite encontrar algo nuevo. *Caminante, no hay camino.* La forma de domar el miedo no es eliminar la incertidumbre, sino convertir el camino difuso en hitos verificables.
- **Criterio visual de investigador.** En tres ocasiones el ojo detectó errores estructurales antes que las matemáticas: (1) el colector que salía del centro hacia la corteza en CCO v8, (2) el hilio que era marcador en lugar de geometría, (3) la necesidad de la hendidura para la forma de frijol. Este instinto es criterio científico real.

### Decisiones técnicas de arquitectura
- **El modelo es un GRAFO, no vóxeles.** Nodos (bifurcaciones) + aristas (segmentos con radio). Pesa kilobytes, no teras. Los vóxeles solo se usan *después* y *localmente* para: ruta de bioink, CFD, difusión de oxígeno.
- **Resolución de 10 µm uniforme = imposible** en hardware actual (~150 mil millones de vóxeles, ~150 GB). Solución: malla adaptativa (octree), refinar a 10 µm solo en ramas terminales.
- **CNN es la herramienta equivocada** para generar el árbol (opera sobre rejillas; el árbol es grafo). La física ya da el óptimo local exacto y gratis (Murray para radios, ángulo óptimo de Zamir). La IA solo aporta en el problema global/combinatorio de topología → GNN/RL sobre grafo, no CNN sobre vóxeles.
- **Generador correcto: space colonization** (Runions et al., 2005) + restricciones físicas. Las ramas crecen hacia puntos de atracción (glomérulos) y se bifurcan; Murray fija radios y ángulo en cada bifurcación; se detiene a escala capilar (~10 µm).
- **Separación ciencia / visualización.** La ciencia vive en `env_biokidney` (NumPy/SciPy) y genera un `.npz`. Blender solo *lee* y dibuja. El `.npz` es el puente. Nunca meter física dentro de `bpy`.

### El mapa de cinco capas (orden de dependencia)
El árbol vascular es la **última** pieza estructural, no la primera: es la *solución* a un problema que definen todas las demás piezas. Orden de ensamblaje:

**forma → nefronas → demanda → árboles → colección → cálices**

| Capa | Contenido | Define |
|---|---|---|
| **0 — Dominio** | Forma del órgano, partición corteza/médula, pirámides, hilio/seno. Es el sistema de coordenadas. | El terreno |
| **1 — Nefronas** | Glomérulos en corteza + túbulos descendiendo a médula. Puntos de atracción anatómicos. | Dónde se consume O₂ y se filtra |
| **2 — Campo de demanda** | Restricción de los 200 µm (límite de difusión de O₂; diseñar a ~100–150 µm por margen). El árbol llena *demanda*, no espacio. | El problema a resolver |
| **3 — Tres árboles** | Arterial (CCO/Murray/space colonization), venoso (espejo, baja presión), colector (convergente, NO usa Murray ni CCO). Capilares de 10 µm emergen aquí. | La solución vascular |
| **4 — Cálices y pelvis** | Destinos de convergencia del colector, definidos por las pirámides de Capa 0. | El drenaje final |

**Restricciones físicas clave:**
- Distancia máxima de difusión de O₂: ~200 µm (margen de diseño 100–150 µm).
- Capilares peritubulares / vasa recta: ~10 µm (escala terminal, alrededor de túbulos, asa de Henle, conducto colector medular).
- Ley de Murray: `r³ = r₁³ + r₂³` en cada bifurcación.
- Verificación complementaria: Finet `d_parent = 0.678(d₁ + d₂)`.

### Progreso de la sesión — Capa 0
**Archivos creados en `08_gemelo_digital/`:**
- `capa0_dominio.py` — generador (env_biokidney, sin bpy) → produce `capa0_dominio.npz`
- `capa0_visualizar.py` — visualizador (bpy, dentro de Blender)

**Geometría del dominio:**
- Elipsoide semiejes a=55, b=30, c=18 mm; eje largo = X (superior-inferior).
- Hilio/seno renal en cara medial; seno esculpido por sustracción de elipsoide de exclusión → forma de frijol.
- Partición final: **59.9% corteza / 40.1% médula** (21.9% de la médula en pirámides).
- 10 pirámides medulares, papilas apuntando al seno (drenarán hacia los cálices).
- Volumen parénquima: 120.447 mm³ (96.8% del elipsoide; seno excluye 3.960 mm³).
- Muestreo: 200.000 puntos, rechazo verificado (0 puntos dentro del seno), SEED=2026.
- `depth(p)`: distancia a la frontera más cercana (cápsula externa O pared del seno).

**Parámetros del seno (en metadata del .npz):** CENTRO_SENO=(0,-34,0), SEMIEJES_SENO=(22,16,11).

### Pendientes y decisiones abiertas
1. **Calibración del seno.** La hendidura actual (3.2% del volumen) es aceptable pero suave. Para figura científica final, considerar subir a SEMIEJES_SENO=(24,18,12) y/o CENTRO_SENO=(0,-32,0) para forma de frijol inequívoca. Límite: no dejar que la cavidad separe los polos en lóbulos. **No bloqueante** — se puede refinar sin rehacer capas superiores (metadata del seno ya registrada).
2. **Mejora pendiente de `depth`.** La definición actual es razonable, pero para la Capa 2 (demanda de O₂) se necesitará **distancia métrica real en mm**, no nivel radial normalizado (importa en elipsoide achatado: la cáscara cortical es más delgada en los polos).
3. **Honestidad científica al presentar.** Esto es una *geometría idealizada paramétrica basada en morfometría poblacional*, NO la reconstrucción de un riñón real por imágenes médicas. Describirlo así fortalece la credibilidad.

### Próximo paso
Capa 1 — sembrar nefronas (glomérulos en corteza, túbulos hacia pirámides). Es el momento en que el dominio deja de ser cáscara vacía y se vuelve órgano con unidades funcionales.

*Entorno: Ubuntu 24.04 · Python 3.12 · NumPy 2.4.3 / SciPy 1.17.1 · env_biokidney · Blender 5.1.1 · Repo maestro `~/Escritorio/BioKidney-AI/`*

---

## ENTRADA 006 — 12 de junio de 2026 — Capa 1 validada + Estado del arte
**Estado:** Capa 1 completada y validada. Sembrado de las unidades funcionales sobre el dominio de Capa 0, más posicionamiento honesto frente al estado del arte.

### 1. Capa 1 completada y validada
- `capa1_nefronas.py` + `capa1_visualizar.py` creados en `08_gemelo_digital/`.
- **1300 nefronas** sembradas: glomérulos en corteza (rechazo + KD-tree, **0 fuera de corteza**), **85% corticales (1105) / 15% yuxtamedulares (195)**. Yuxtamedulares más profundos (depth medio **0.266 vs 0.085**) y con túbulos más largos (asa larga hacia papila).
- Túbulos: polilíneas de 6 vértices descendiendo hacia la pirámide más cercana. Poisson-disk `DIST_MIN=0.8mm` sin solapes.
- Salida `capa1_nefronas.npz` (`glomerulos`, `tipo`, `tubulos`, `piramide_destino`, `depth_glomerulo`).
- **Validación visual:** corte sagital (Y, grosor 6mm) confirmó estratificación correcta — glomérulos forman banda en corteza externa, NO rellenan médula; túbulos descienden corteza→pirámides. **Cuarta validación visual exitosa del proyecto.**

### 2. Pendiente de recalibración de Capa 0 (acumular para una sola pasada, NO bloqueante)
- (a) Seno a `(24,18,12)` opcional.
- (b) `depth` a distancia métrica real en mm (necesario para Capa 2).
- (c) **NUEVO:** sesgo polar en drenaje por pirámide (`piramide_00` y `09` captan ~21%/20%, centrales ~5%) — artefacto de geometría de pirámides de Capa 0; equilibrar por capacidad o posición.

### 3. Posicionamiento frente al estado del arte (investigado esta sesión)
- El modelado computacional de árboles vasculares renales es un campo **activo y poblado**. Grupos trabajando en piezas similares: CCO/GCO aplicado a riñón; generación de redes a escala de órgano completo (corteza humana, nivel pre-arteriolar); marco multiescala in silico de hemodinámica renal con tomografía de contraste de fase (bioRxiv dic 2025); pipeline de vasculatura sintética 230× más rápido con CFD + fabricación 3D.
- **Implicación:** NO posicionar el proyecto como "primero en modelar árbol vascular renal" (insostenible). El aporte real probablemente está en la **integración** (vascular + colector + nefronas + función en marco abierto, reproducible, hardware modesto), no en inventar generación vascular. Definir aporte exacto en fase de publicación.
- **Cuello de botella de vascularización:** confirmado como reto #1 en bioimpresión de órganos (literatura 2025). PERO tiene dos mitades: (a) **geométrica/diseño** — qué arquitectura cumple la restricción de difusión de O₂ — que el gemelo digital SÍ ataca; (b) **biológica** — maduración celular, anastomosis con huésped, tipos celulares (endotelio/pericitos/estroma), perfusión sin necrosis — que NINGÚN modelo digital resuelve solo, requiere lab húmedo. El gemelo digital resuelve "el plano, no el ladrillo". Posicionamiento honesto = credibilidad.

### Siguiente paso
Capa 2 — campo de demanda de oxígeno (restricción 200µm, diseñar a 100–150µm).

---

## ENTRADA 007 — 14 de junio de 2026 — Capa 2 (campo de demanda) + visualizador maestro
**Estado:** Capa 2 completada y validada. El PROBLEMA del gemelo digital queda completamente planteado (forma + nefronas + demanda + regla de cobertura). Falta la SOLUCIÓN (Capa 3).

### 1. Visualizador maestro
- Creado `visualizador_maestro.py` en `08_gemelo_digital/`, unifica `capa0_visualizar.py` + `capa1_visualizar.py` (que se conservan como respaldo).
- Capas activables por flags booleanos (`VER_DOMINIO`, `VER_PIRAMIDES`, `VER_HILIO`, `VER_GLOMERULOS`, `VER_TUBULOS`, `VER_DEMANDA`) + corte sagital (`CORTE_ACTIVO/EJE/CENTRO/GROSOR`).
- **Carga perezosa:** si falta un `.npz`, desactiva su flag y continúa. Cada capa va a su propia colección de Blender (`Coleccion_Dominio`, `_Nefronas`, `_Demanda`) para control desde el outliner.
- **`clear_scene` endurecido:** borra vía `bpy.data.objects.remove` (no por selección/`bpy.ops`) con try/except por objeto, limpia mallas huérfanas y colecciones del script. Resuelve el `RuntimeError` de View Layer al re-ejecutar con distintos flags.

### 2. Capa 2 — campo de demanda metabólica (completada y validada)
- Creado `capa2_demanda.py` (env_biokidney, sin bpy). Lee capa0 + capa1, construye campo escalar de demanda de O₂ sobre los 200.000 puntos en coordenadas reales (mm).
- **Construcción:** base por región (cortex=1.0, médula=0.4) + realce local por nefronas (KD-tree con glomérulos+túbulos, kernel lineal en radio 0.5mm, realce normalizado por percentil 99 para que la base regional domine).
- **Decisión clave de normalización:** divide-por-máximo (NO min-max), para conservar demanda medular residual (~0.31). Justificación: la médula NO es avascular (requiere vasa recta para las asas yuxtamedulares); con min-max caería a 0 y el árbol de Capa 3 podría dejar la médula seca reportando éxito falso.
- **Resultado:** cortex media 0.773, médula/pirámides ~0.308 (=0.4/1.3). Histograma bimodal limpio. Contraste cortex>médula preservado.
- **Restricción de diseño fijada:** ningún punto a más de 150µm (0.150mm) de un capilar — margen de seguridad bajo el límite de hipoxia (~200µm).
- Función `cobertura(capilares)` lista para Capa 3: mide % tejido cubierto (<0.150mm), distancia media, peor caso, y déficit ponderado por demanda. Maneja input vacío.
- **Línea base ("antes")** con glomérulos como capilares-proxy: **0.66% cubierto, déficit 98.87%, dist. media 3.56mm, peor caso 12.5mm.** Esta es la métrica que el árbol de Capa 3 vendrá a mejorar — la mejora será un resultado central del paper.
- Salida `capa2_demanda.npz` (coords, demanda[0,1], metadata, métricas proxy).

### 3. Validación visual (quinta del proyecto)
- Mapa de calor activado en el maestro (`draw_demanda`, colormap inferno/turbo por NumPy). Corte sagital confirma: banda cortical caliente (amarillo) envolviendo núcleo medular frío (rosa/morado), con transición por realce de nefronas yuxtamedulares.
- Figura candidata a paper: la sección térmica comunica el principio **"el árbol llena demanda, no espacio"**.
- Las 4 capas se visualizan coherentes simultáneamente (glomérulos donde hay alta demanda, túbulos hacia pirámides, mapa respeta frontera cortex-médula).

### Siguiente paso
Capa 3 — los tres árboles vasculares (arterial CCO/Murray/space colonization; venoso espejo; colector convergente sin Murray). El corazón del proyecto; merece sesión dedicada.

---

## ENTRADA 008 — 19 de junio de 2026 — Reconciliación de cifras y estado canónico
**Estado:** Entrada de consolidación, no de avance técnico. Las cifras de las entradas anteriores se conservan intactas como registro fiel del proceso; esta entrada solo fija cuál es el valor vigente de cada parámetro y por qué cambió. Sirve como referencia rápida de "qué dato manda hoy y por qué".

### Nota de apertura
Ninguna cifra de las entradas previas se modifica ni se corrige en su sitio: la bitácora es un registro histórico y su valor depende de que las entradas previas queden intactas. Lo que sigue es una reconciliación de superficie — distingue el **valor histórico** (ya registrado en entradas anteriores) del **valor vigente**, sin tocar el original. Cuando una entrada antigua y esta discrepen, manda esta para el dato vigente; la entrada antigua manda para saber qué se creía en su momento.

### Tabla de evolución de cifras
| Parámetro | Valor histórico → Valor vigente | Razón del cambio |
|---|---|---|
| **TFG bilateral** | 82 mL/min (presión sintética, no derivada de geometría real) → **115,2 mL/min** (v8, modelo de Poiseuille calibrado en dos pasadas, dentro del rango normal 100–125) | El 82 no provenía de la geometría vascular real. | [RETIRADO: ver MAESTRO §0]
| **Referencia "nativo"** | "62,5 mL/min" como comparación → se **elimina** el encuadre de "superar al nativo" | Comparación engañosa; el marco correcto es "dentro del rango fisiológico normal". |
| **WSS (wall shear stress)** | 24,69–37,04 dyn/cm² (banda 10–70) → **5,6 dyn/cm²** (banda renal real 1–10) | La banda 10–70 era fisiológicamente incorrecta para arteriolas renales; ya corregido en marzo. |
| **Vascular** | CCO v7 / 1.448 segmentos → **CCO v8 / 1.902 segmentos** (pipeline publicado), con reconstrucción en curso por capas (Capa 3, space colonization + Murray) | Mejora de densidad cortical y presión; v7 es histórico. |
| **Módulos** | "12 módulos" (pipeline granular) ≡ "6 módulos" (escalas biofísicas) | Es la misma obra a distinta resolución; para comunicación se usan **6**. |
| **DOI** | `…19508077` = DOI de versión específica · `…19508076` = concept DOI (todas las versiones) | Para citar por defecto se usa el **concept DOI `…076`**. |

### Nota de posicionamiento
El aporte del proyecto es la **integración abierta y reproducible sobre hardware modesto**, no la invención de la generación vascular. El gemelo digital resuelve **"el plano, no el ladrillo"**: la mitad geométrica/de diseño del problema de vascularización, no la mitad biológica (maduración celular, anastomosis, perfusión sin necrosis), que requiere lab húmedo. El trabajo es **Fase 1, completamente in silico**.

### Estado a la fecha (19 de junio de 2026)
- **Capas 0–2 del gemelo digital: completadas y validadas** (dominio, nefronas, campo de demanda).
- **Capa 3 (árbol vascular): frontera inmediata, aún sin empezar.** Es el corazón del proyecto y la solución al problema que las capas 0–2 dejan completamente planteado.

---

## ENTRADA 009 — 21 de junio de 2026 — Capa 3a: árbol arterial (generación, corrección, auditoría)
**Estado:** Capa 3a completada, corregida y auditada. Primer árbol vascular del gemelo digital. El árbol arterial canónico queda como BINARIO, fisiológicamente verosímil y con Ley de Murray validada por verificación cruzada. (Nota de numeración: esta es la entrada 009 — la 008 ya estaba ocupada por la reconciliación del 19 de junio; se continúa la secuencia.)

### 1. Árbol arterial por space colonization
- Creado `capa3a_arterial.py` (`env_biokidney`, sin bpy). Genera el árbol arterial por **space colonization** (Runions et al. 2005) desde el hilio hacia los 1300 glomérulos (Capa 1) usados como **puntos de atracción**. Representación como **GRAFO** (nodos + aristas con radio), no vóxeles.
- Radios asignados por **Ley de Murray** de hojas a raíz: `r_padre³ = Σ r_hijos³` (exponente k=3).
- **Bootstrapping de tronco:** el hilio queda a **11.3 mm** del glomérulo más cercano (fuera del radio de influencia), así que se añadió una fase de **tronco recto** desde el hilio hasta que la copa entra en rango — anatómicamente correcto (la arteria renal viaja antes de ramificar).

### 2. Defecto detectado y corregido (mega-hubs)
- **v1 inicial:** Murray cumplía 100%, pero la **auditoría** (`capa3a_auditoria_murray.py`) destapó un defecto estructural: **46 nodos "mega-hub" con hasta 1969 hijos** cada uno, concentrando **94% de las aristas**. Topología inverosímil (tronco → abanicos gigantes, no ramificación recursiva).
- **Lección epistemológica registrada:** que Murray cumpla "por construcción" es un test de **CONSISTENCIA**, no de validación anatómica. La visualización sirve para inspeccionar, no para demostrar. La **auditoría numérica + la inspección visual juntas** son lo que valida.
- **Corrección:** `GRADO_MAX_HIJOS` limita los hijos por nodo; al cerrarse un nodo, sus atractores se **reasignan a los hijos** → ramificación recursiva real. Además, **optimización del cKDTree** (reconstruir cada 10 iteraciones en vez de cada una): corrida de **~10 min → ~8 s** (factor ~55×).

### 3. Decisión binario vs trinario (validación empírica)
- Se permitió inicialmente hasta **trifurcaciones** por argumento bioingenieril (priorizar irrigación de demanda sobre anatomía estricta, conservando Murray generalizado `r_p³ = Σ r_i³`). Confirmado además que **co-SWIFT** (Stankey et al. 2024, bioprinting coaxial) usa Murray en sus bifurcaciones y que las multifurcaciones pueden descomponerse en bifurcaciones en cascada para impresión.
- **v2 trinario** salió **99.8% trifurcaciones** (artefacto: cada nodo llena su cuota de 3 antes de cerrar). Se comparó empíricamente contra **v3 binario** (`GRADO_MAX_HIJOS=2`).
- **Resultado:** el binario **iguala el 100% de glomérulos alcanzados** con la **mitad de terminales** (5652 vs 11341), **100% bifurcaciones binarias** (fisiológico), más imprimible (menos nodos de unión co-SWIFT) y más rápido. El trinario **no se justifica por cobertura**.
- **Decisión:** binario adoptado como **definitivo** (`GRADO_MAX_HIJOS=2`).

### 4. Estado final Capa 3a (canónico = binario)
- `capa3a_arterial.npz`: **11971 nodos, 5652 terminales, grado máx 2, 100% glomérulos alcanzados (1300/1300)**.
- **Auditoría Murray:** 100% de **5651 bifurcaciones** cumplen dentro del 1/5/10% de tolerancia; **exponente empírico k=3.0000** (verificación cruzada, sin bug de implementación); **0 mega-hubs**.
- **Respaldos conservados:** `v1_megahubs` (defecto), `v2_trinario`, `v3_binario` (=canónico). **No borrar.**
- **Validado visualmente** en el visualizador maestro (`draw_arterial`): árbol genuinamente ramificado desde el hilio, no abanicos.

### 5. Pendiente de calibración (no bloqueante, acumular con los de capas previas)
- **Radio de raíz ~0.21 mm**, lejos del rango anatómico de arteria renal (~2–3 mm). Es Murray correcto sobre el nº de terminales con radio terminal de 12 µm, pero para realismo anatómico habrá que **calibrar** (subir radio terminal o nº de generaciones) en la pasada futura de calibración, junto al **seno (24,18,12)**, **depth métrico** y **sesgo polar de pirámides**.

### Siguiente paso
Capa 3b — árbol venoso (espejo del arterial a baja presión).

---

## ENTRADA 010 — 22 de junio de 2026 — Capa 3a-bis: red peritubular / vasa recta
**Estado:** Capa 3a-bis completada. Segunda red capilar que actúa de **puente arterial→venoso**: cierra el trayecto entre las terminales arteriales (arteriolas aferentes) y el futuro árbol venoso, y resuelve el déficit medular que el árbol arterial dejaba. (Nota de numeración: la 009 ya estaba ocupada por la Capa 3a; se continúa la secuencia.)

### 1. Segunda red capilar (puente arterial→venoso)
- Creado `capa3ab_peritubular.py` (`env_biokidney`, sin bpy). Genera la red **peritubular + vasa recta** como **puntos de drenaje**: el destino de la sangre tras el filtrado glomerular y el origen del retorno venoso.
- **4290 puntos de drenaje**, repartidos en **77% peritubular cortical** (tapiza los túbulos en la corteza) / **23% vasa recta** (desciende a la médula siguiendo el asa de Henle).
- **Distribución por región:** **55% en médula**; el **62.5% de la vasa recta cae en médula** → **resuelve el déficit medular** que el árbol arterial (dirigido a los glomérulos, mayormente corticales) dejaba sin servir.
- Estos 4290 puntos son el **origen del árbol venoso** (Capa 3b): la Capa 3b crece por space colonization tomándolos como puntos de atracción.
- Salida: `capa3ab_peritubular.npz` (puntos_drenaje, tipo_red).

### Siguiente paso
Capa 3b — árbol venoso convergente, que drena estos 4290 puntos hacia el hilio.

---

## ENTRADA 011 — 23 de junio – 4 de julio de 2026 — Capa 3b: árbol venoso (CERRADA)
**Estado:** **Capa 3b CERRADA.** Tercer y último árbol vascular del gemelo digital. Cierra el circuito **arterial → peritubular → venoso**: los 4290 puntos de drenaje de la Capa 3a-bis se recogen en un árbol venoso convergente al hilio, reparado frente a colisiones arterio-venosas y auditado exhaustivamente frente a auto-colisiones venosas.

### 1. Generación del árbol venoso
- Creado `capa3b_venoso.py` (`env_biokidney`, sin bpy). **Space colonization convergente al hilio** (inverso al arterial: la sangre converge en vez de distribuirse), **satelital al arterial** (crece junto al árbol arterial para el trenzado anatómico arteria-vena).
- **25716 nodos, 12150 terminales, 100% de puntos de drenaje alcanzados** (córtex 100%, médula 100%) → **circuito arterial→peritubular→venoso completo**.
- **Binario, 0 mega-hubs** (misma disciplina topológica que el arterial).
- **Radio de raíz venoso 0.538 mm > arterial 0.214 mm**: la vena es **más gruesa** (retorno de baja presión, mayor sección para el mismo caudal). **Ley de Murray k=3.0 al 100%.**
- Archivo canónico: `capa3b_venoso.npz` (respaldo pre-reparación: `capa3b_venoso_v1_precolisiones_backup.npz`).

### 2. Reparación de colisiones arterial-venoso (AV)
- **67 colisiones AV → 0.** Separación mínima **0.0000 mm → 0.0441 mm** (los vasos ya no se tocan).
- **Murray y topología intactos**; solo **0.59% de nodos movidos** (desplazamiento medio 0.23 mm). **Reparación local, sin re-generar** el árbol.

### 3. Auditoría de auto-colisiones venoso-venoso (VV)
- Tras **5 auditorías sucesivas**, los **1927 shunts de calibre** (que funden a 100 µm, **99.5% intra-lobular**) se descomponen en:
  - **79.2% (1527) confluencias en Y anatómicamente obligatorias:** dos venas tributarias convergiendo al colector — **paralelas hacia la raíz (2.1°), divergentes hacia la hoja (34.6°)** = anatomía correcta, **NO defecto**.
  - **19.6% (378) ruteo-corregible-separable** (cota superior floja).
  - **1.1% ambiguos.**
- **Veredicto:** no es defecto estructural sino **mayoritariamente anatomía**; la cola corregible tiene techo flojo y su **único fix difuso probado (repulsión de proximidad, barrido W_REP 0→4) empeora el conteo monótonamente**. **No justifica desarrollo dirigido.**
- **Descartados como fix:** repulsión local, convergencia jerárquica por zonas, puertos venosos múltiples.

### 4. Validación visual (Blender, `draw_venoso`: arterial en rojo / venoso en azul)
- Par **arteria-vena entrelazado sin fundirse** (confirma la reparación AV), **convergencia venosa al hilio**, y **confluencias en Y anatómicas** visibles. [OK]

### 5. Nota para Capa 4 (voxelización)
- Las **1504 confluencias en Y** son candidatas a **uniones topológicas explícitas con luz** al voxelizar.
- **Artefacto conocido:** **todos los segmentos miden exactamente 0.6 mm** (`paso_crecimiento` fijo) — relevante para la voxelización.

### 6. Hallazgo documentado — cobertura de difusión (NO bloqueante)
- **Cobertura a 150 µm = 0.78%** (vs **0.66% base**). Es **límite de resolución de la escala representativa** (segmentos 0.6 mm = **4× el radio de difusión**; el lecho capilar no está en la geometría **por diseño**), **NO fallo del árbol**.
- La validación de la vasculatura es **geométrica** (100% de puntos de demanda alcanzados, Murray), **no de difusión** — esta última es métrica de la **etapa capilar futura**.
- **Vacío avascular central genuino de ~5 mm** (no es el seno; hipótesis: sesgo polar de drenaje).
- Ambos → **pendientes de la pasada de calibración**.

### Estado
**Capa 3b CERRADA.** El circuito vascular completo (arterial → peritubular → venoso) queda cerrado, reparado y auditado. Siguiente frontera: Capa 4 (voxelización con luz).

---

## ENTRADA 012 — 4 de julio de 2026 — Capa 3c: árbol colector urinario (CERRADA)
**Estado:** **Capa 3c CERRADA.** Cuarto árbol del gemelo digital y primer sistema **no vascular**: el colector urinario. Cierra el trayecto de la orina **glomérulo → túbulo → pirámide → papila**. A diferencia de los árboles de presión (3a/3b), **NO usa Ley de Murray ni conservación de área**: es un árbol de convergencia de flujo con radios por Horton-Strahler.

### 1. Arquitectura — bosque de 10 subárboles disjuntos
- Creado `capa3c_colector.py` (`env_biokidney`, sin bpy). **Un subárbol independiente por pirámide** (k=0..9), sin conexión entre ellos: cada pirámide drena a su propia papila.
- Por pirámide: **raíz = `pyramid_apex[k]`** (la papila, Capa 0 — verificado en auditoría previa: el apex está a ~13 mm del hilio, apunta al seno) y **atractores = `tubulos[piramide_destino==k, 5]`** (el **extremo distal** de cada túbulo, punto 5 de la polilínea; `piramide_destino` se **hereda** de Capa 1, no se recalcula).
- **Reúso de código:** se importa y reutiliza **`space_colonization()` de `capa3b_venoso.py`** (misma máquina depurada de Runions 2005 + bifurcación binaria + bootstrap de tronco), llamada con `art_tree=None` (**sin sesgo satelital**). **NO** se usa `asignar_radios_murray()`. Los 10 subárboles se concatenan con offset de índices de `parent`; `parent` de cada raíz = −1; `piramide_id[n]` y `papila_nodo[k]` (índice global de cada raíz) trazan la pertenencia.

### 2. Hiperparámetros
- Arrancan del venoso: **PASO_CRECIMIENTO = 0.6 mm**, DIST_INFLUENCIA = 8.0 mm, **GRADO_MAX_HIJOS = 2** (binario). Ajustado **DIST_MATAR = 0.6 mm** (más ajustado que el venoso 1.2, para que el colector se **acerque** al extremo del túbulo) → maximiza alcance.
- **Radios por Horton-Strahler** (bloque intercambiable): orden por nodo (hojas=1; +1 si el máximo empata en ≥2 hijos); radio = rampa lineal monótona entre **`radio_terminal` = 15 µm** (orden 1) y **`radio_papila` = 125 µm** (orden máx global). Placeholders **calibrables**, documentados como escalares en el `.npz`. *(Nota: para cambiar a downstream-count se sustituye SOLO este bloque por una rampa en log(N_terminales_aguas_abajo).)*

### 3. Salida canónica
- `capa3c_colector.npz`: **2086 nodos, 2076 aristas (= N − 10 subárboles), 853 terminales**. Campos: `nodos(N,3) parent(N,) radios(N,) piramide_id(N,) papila_nodo(10,) terminales(T,) conexion_tubulo(1300,2) conexion_dist(1300,) strahler(N,)` + escalares (`paso_crecimiento, n_subarboles=10, radio_terminal, radio_papila, seed, n_no_alcanzados_por_piramide(10,)`). **Radios por-nodo**, misma convención que el venoso (r_in por-nodo = radio del segmento que alimenta al nodo).

### 4. Auditoría (4 bloques; solo mide, no corrige)
- **(1) Alcance:** **100% de atractores alcanzados en las 10 pirámides** (0 no alcanzados de 1300, análogo al "100% glomérulos" del arterial). Distancia túbulo→colector `conexion_dist`: **p50 0.155 mm, p90 0.267, p99 0.296, máx 0.408 mm**.
- **(2) Continuidad:** **1300/1300 nefronas conectadas**, el **100% al subárbol de su propia pirámide** (sin fuga inter-territorio). Las 10 `papila_nodo[k]` coinciden **exactamente** con `pyramid_apex[k]` (desvío 0). [OK]
- **(3) Radios:** orden Strahler [1, 6]; radios [15, 125] µm; **monotonía raíz>hoja al 100%** (0 violaciones). Papilas 81–125 µm (las de menor orden son más finas). Compuerta `radio_papila` = 125 µm dentro de 100–150 µm. [OK]
- **(4) Proximidad inter-subárbol** (análogo shunt VV, **sin aplicar repulsión**): **0 pares por debajo de 150 µm**; distancia mínima entre subárboles distintos **5.29 mm**. Los subárboles están **bien separados** (bosque disjunto): a diferencia del venoso (árbol único), aquí no hay cruces — no aplica problema de shunts.

### Estado
**Capa 3c CERRADA.** El colector urinario queda validado geométricamente (100% de túbulos drenados, papilas ancladas, radios monótonos, subárboles disjuntos). Pendiente (paso aparte): visualizador. Frontera: Capa 4 (cálices/pelvis + voxelización) — su contrato de entrada aún no está definido; se diseñará sobre esta salida.

---

## ENTRADA 013 — 5 de julio de 2026 — Capa 3c: recalibración de radios (anclaje anatómico + rampa por-subárbol)
**Estado:** **Recalibración quirúrgica de los radios del colector.** SOLO cambian los calibres; **la arquitectura, la topología, el space colonization y la conexión con Capa 1 quedan intactos**. Verificado por no-regresión: el `.npz` es **bit-idéntico** en toda la topología, solo cambian `radios` y el escalar `radio_papila`.

### 1. Qué cambió (y qué NO)
- **`radio_terminal` = 15 µm (SE MANTIENE).** Diámetro 30 µm ∈ conducto **colector cortical** 20–50 µm ⌀ (morfometría del colector cortical, Britannica / literatura morfométrica renal).
- **`radio_papila` = 200 µm (SUBE desde 125).** Diámetro 400 µm = **centro del conducto de Bellini** (300–600 µm ⌀; *grading system* de nefrolitiasis, PMC). Un conducto de Bellini por pirámide, en el ápice papilar.
- **Rampa Strahler ANCLADA POR-SUBÁRBOL (único cambio de lógica).** Antes la rampa lineal se anclaba al **orden Strahler máximo global** (omax=6), de modo que las papilas de subárboles menos ramificados quedaban más finas (rango **81–125 µm**). Ahora, para cada subárbol *k* el radio se interpola entre `radio_terminal` (orden 1) y `radio_papila` en el **orden máximo de ESE subárbol** (omax_k). Como la raíz (papila) tiene siempre el orden máximo de su subárbol, cae exactamente en `radio_papila`. Es el bloque intercambiable que ya estaba marcado; no se tocó nada más.
- **Justificación anatómica del anclaje por-subárbol:** las **10 papilas son anatómicamente equivalentes** (un Bellini por pirámide); su calibre **no debe depender de la profundidad de ramificación** del subárbol, que es un artefacto del muestreo de nefronas por pirámide.
- **NO tocado:** nodos, `parent`, terminales, `papila_nodo`, `piramide_id`, `conexion_tubulo`, `conexion_dist`, `strahler`, hiperparámetros del space colonization. Todo bit-idéntico a Entrada 012.

### 2. Nota de compuerta / escala vs resolución de impresión
- La compuerta se **reencuadra** contra la **resolución de impresión (100–150 µm de rasgo, en diámetro)**, no contra morfometría. La resolución restringe el rasgo **FINO** (terminal), no el **GRUESO** (papila).
- **Papila** ⌀400 µm: **por encima** de 100–150 µm → resoluble [OK].
- **Terminal** ⌀30 µm: **SUB-RESOLUCIÓN por diseño** [NOTA], **no es fallo**. Es el mismo encuadre de escala representativa que la vasculatura: el segmento terminal representa el **drenaje del lóbulo**, no el conducto individual. La resolución de impresión restringe el rasgo fino; la papila (rasgo grueso) queda siempre por encima.

### 3. Regeneración y no-regresión
- Regenerado `capa3c_colector.npz` con `radios` nuevos y escalar `radio_papila` actualizado; **mismas claves, misma topología**. Verificación campo a campo (prev vs. new): `nodos, parent, terminales, papila_nodo, piramide_id, conexion_tubulo, conexion_dist, strahler, n_no_alcanzados_por_piramide` → **IDÉNTICOS**. Solo cambian `radios` ([15,125]→[15,200] µm) y `radio_papila` (125→200 µm). **[OK] TOPOLOGÍA BIT-IDÉNTICA.**

### 4. Re-auditoría (4 bloques)
- **(1) Alcance/continuidad — IDÉNTICO (confirmado):** **0/1300 no alcanzados** (100 % en las 10 pirámides), **1300/1300 nefronas conectadas** (100 % al subárbol de su pirámide), 10 `papila_nodo[k]` = `pyramid_apex[k]` (desvío 0). Sin cambios respecto a Entrada 012 (la topología no cambió).
- **(2) Radios:** rango nuevo **[15, 200] µm**, orden Strahler [1, 6], **monotonía raíz>hoja 100 % (0 violaciones)**. **Las 10 papilas ahora 200 µm homogéneas** (min=max=200), frente al rango 81–125 µm del anclaje global. Órdenes Strahler en las papilas: [6,5,5,4,4,4,4,4,5,5] (cada papila = orden máximo de su subárbol).
- **(3) Compuerta:** papila ⌀400 µm resoluble [OK]; terminal ⌀30 µm sub-resolución por diseño [NOTA] (ver §2).
- **(4) Proximidad inter-subárbol** (recalculada, los radios cambian los tubos; **sin aplicar repulsión**): **0 pares < 150 µm**, distancia mínima inter-subárbol **5.29 mm**. Sigue sin pares problemáticos: subárboles bien separados (bosque disjunto).

### Estado
**Recalibración cerrada.** Radios del colector anclados a literatura anatómica (colector cortical 20–50 µm ⌀; Bellini 300–600 µm ⌀), papilas homogéneas por anclaje por-subárbol, topología intacta y verificada por no-regresión. La Capa 3c sigue **CERRADA**; frontera Capa 4 sin cambios.

---

## ENTRADA 014 — 5 de julio de 2026 — HITO: corrección de lógica raíz (seno ≠ córtex) + regeneración completa de la cascada Capa 0→3c
**Estado:** **Hito mayor completado y PROMOVIDO a canónico.** Corrección de una raíz geométrica del gemelo digital y regeneración íntegra de las 8 capas, con esquema de rollback (backup en `08_gemelo_digital/_backup_pre_seno/`, que **permanece** como estado histórico). Precedido de tres diagnósticos de solo lectura (inventario dimensional, causa-raíz del seno, seguridad del campo de profundidad) que anclaron cada decisión.

### 1. Corrección de lógica raíz en Capa 0 (`capa0_dominio.py`)
- **Problema:** `nearest_surface_distance()` definía la profundidad como el **mínimo** entre distancia a la cápsula externa y a la **pared del seno renal**. Consecuencia: puntos INTERIORES pegados al seno recibían profundidad baja y se etiquetaban **cortex** aunque estuvieran a hasta **21 mm** de la cápsula (córtex peri-sinusal falso). El seno es superficie **interna** (bajo ella hay médula / grasa sinusal, NO córtex).
- **Fix:** la profundidad **cortical** se mide **SOLO contra la cápsula externa** (`capsule_distance()`). El seno **conserva** su rol geométrico (exclusión de parénquima y `build_pyramids`), que no usan este campo — verificado en el diagnóstico de seguridad.
- **Umbral fracción → mm ABSOLUTOS:** nuevo `GROSOR_CORTICAL_MM = 6.6` (cortical width, MDCT n=2068). `cortex = depth_cortical_mm < 6.6`. Se **depreca** `UMBRAL_CM` como fracción (se conserva en el `.npz` como fracción equivalente = GROSOR/depth_norm por compatibilidad de lectura). El umbral absoluto es **invariante al normalizador**: evita la trampa detectada en el diagnóstico (al separar el seno, `depth_norm` salta 28→52 mm; con umbral en fracción eso reclasificaba 43 896 puntos médula→córtex). Con umbral en mm: **0** reclasificaciones espurias.

### 2. Reanclaje de Capa 1 (`capa1_nefronas.py`) — ajuste acompañante obligatorio
- La banda yuxtamedular se reancla a **mm**: `juxta_min_depth_mm = GROSOR·0.75` → tercio profundo del córtex **[4.95, 6.60] mm**. El peso cortical usa `CORTICAL_DEPTH_SCALE_MM = 3.3 mm` (≈GROSOR/2). Pool de siembra = `region=="cortex"` (ya correcto desde Capa 0).
- **`FRAC_YUXTAMEDULAR = 0.15` sin cambios** (consenso 85/15). **Resultado: 1105/195 = 85.00/15.00 exacto; 0 glomérulos fuera del córtex.** El córtex más delgado (6.6 mm) **NO** expulsó yuxtamedulares → no hizo falta subir a 7.0 mm.

### 3. Regeneración completa en orden de dependencia (mismos hiperparámetros canónicos)
`capa0 → capa1 → capa2 → capa3a → capa3ab → capa3b (crudo → auditoría → reparación → promoción) → capa3c`. Solo cambió lo que deriva de la corrección de Capa 0; hiperparámetros vasculares (paso 0.6 mm, Murray k=3.0) y del colector (Strahler radio_terminal 15 µm / radio_papila 200 µm por-subárbol) **sin tocar**.

### 4. Fix `MAX_PASADAS` 3→4 en el reparador AV (`capa3b_reparar_colisiones.py`)
- Con la geometría regenerada, la reparación AV daba **33→1** (residual) con el tope de 3 pasadas, donde el pipeline viejo daba 67→0. Diagnóstico de solo lectura: la reparación local **OSCILA** (empujar un segmento + su suavizado a 1–2 saltos crea colisiones transitorias en los vecinos; secuencia **33→3→6→1→0**), y **convergía en la 4ª pasada**; el tope cortaba una antes.
- **Cambio único:** `MAX_PASADAS 3→4`. **NO** se tocaron MARGEN/SUAVE/OFFSET ni el factor de detección — física del reparador idéntica. Cada pasada solo empuja segmentos **aún** en colisión → las 32 reparaciones ya buenas no se re-tocan: **0 regresiones**. La colisión residual era un cruce oblicuo (39.8°) de dos vasos **finos** (arteriola 47 µm / vénula 29.5 µm) en la periferia cortical (56.9 mm del hilio), no un caso de troncos gruesos.
- **Resultado:** **AV 33→0**, `sep_min = 0.0525 mm` (sana, >0); VV 25240→25209 (monótono, sin nuevas); Murray y topología intactos. Auditoría independiente sobre el `.npz` promovido: **0 colisiones**.

### 5. Tabla comparativa viejo → nuevo (los valores cambian; es informativo)
| Métrica | Viejo | Nuevo |
|---|---|---|
| capa0 córtex / médula / pirámide (pts) | 119831 / 36370 / 43799 | 92624 / 55788 / 51588 |
| capa0 umbral cortico-medular | 0.30 fracción (≈8.5 mm, contaminado por seno) | **6.6 mm absoluto** |
| capa1 cortical / yuxta | 1105 / 195 | 1105 / 195 |
| capa3a nodos / aristas / terminales | 11971 / 11970 / 5652 | 12016 / 12015 / 5699 |
| capa3a glomérulos alcanzados / radio raíz | 1300 / 0.2138 mm | 1300 / 0.2143 mm |
| capa3ab drenaje córtex / médula | 1931 / 2359 | 1191 / 3099 |
| capa3b nodos / aristas / terminales | 25716 / 25715 / 12150 | 23799 / 23798 / 11212 |
| **capa3b colisiones AV** | **67 → 0** | **33 → 0** (con MAX_PASADAS=4) |
| capa3c nodos / terminales | 2086 / 853 | 2189 / 903 |

### 6. Resultado peri-sinusal: patología ELIMINADA
- **Córtex peri-sinusal PROFUNDO (glomérulos cortex a >6.6 mm de la cápsula): 0** (antes: glomérulos hasta 21 mm de profundidad marcados cortex). Profundidad cortical máxima de glomérulos ahora **6.59 mm**.
- Residual "nearest-surface = seno": **37 / 1300 (2.85 %)** frente a 152 (11.69 %) — son la **esquina cápsula/seno** (reborde hiliar), todos a <6.6 mm de la cápsula: **córtex legítimo**, no interior falso. El sobreconteo de esquina ya estaba anticipado en el diagnóstico.

### 7. Reauditoría final de invariantes (todos [OK])
- **capa0:** grosor 6.6 mm exacto; córtex más profundo 6.600 mm (<6.6). [OK]
- **capa1:** 85/15 exacto; 0 glomérulos fuera del córtex. [OK]
- **capa3a:** 1300/1300 glomérulos alcanzados; Murray k=3.0; sin mega-hubs. [OK]
- **capa3ab:** 4290 puntos de drenaje conectados (córtex 1191 / médula 3099). [OK]
- **capa3b:** **AV = 0** (el que fallaba); sep_min 0.0525 mm; drenaje 4290/4290; sin mega-hubs; VV no aumentó (25240→25209); convergente al puerto hiliar [0,−30,+3]; aristas=nodos−1. [OK]
- **capa3c:** 10 subárboles; 0 atractores no alcanzados; 1300/1300 nefronas conectadas; papila=apex (desvío 7e-7 mm); papilas 200 µm homogéneas; monotonía 0 violaciones; 0 pares inter-subárbol <150 µm. [OK]

### Estado
**Hito PROMOVIDO a canónico.** Los 8 `.npz` regenerados son ahora los canónicos en su ruta. El backup `08_gemelo_digital/_backup_pre_seno/` **permanece** como rollback histórico (no se borra). La corrección elimina el córtex peri-sinusal falso, ancla el umbral cortico-medular a un valor anatómico absoluto (6.6 mm, MDCT n=2068) e invariante al normalizador, y toda la cascada vascular/colectora queda regenerada y validada. Frontera Capa 4 sin cambios.

---

## ENTRADA 015 — 6 de julio de 2026 — Capa 4: sistema calicial alto (papila → cáliz → pelvis → uréter)
**Estado:** **Capa 4 generada y validada.** Primer tramo del colector **alto**: cierra el trayecto de la orina desde la papila hasta el hilio (`papila → cáliz menor → cáliz mayor → pelvis → uréter`). Representación **centerline + radio** (grafo), convergente hacia el hilio. **NO voxeliza** (eso es Capa 5). Solo se generó Capa 4; capas 0–3c intactas.

### 1. Arquitectura — bicalicial, radios ABSOLUTOS
- Creado `08_gemelo_digital/capa4_colector_alto.py` (`env_biokidney`, sin bpy). Salida `capa4_colector_alto.npz`.
- **24 nodos, 23 aristas** (árbol, E = M−1): 10 `papila_junction` (nivel 0) + 10 `caliz_menor` (nivel 1) + **2 `caliz_mayor`** (nivel 2, patrón **bicalicial** por split del eje polar X) + 1 `pelvis` (nivel 3) + 1 `ureter` (nivel 4).
- **Radios absolutos en mm** (misma decisión que Capa 3c, no proporcionales): Bellini 0.200 (heredado de 3c) · copa 1.5 · infundíbulo 2.0 · pelvis 4.5 · uréter 1.5.
- **Construcción determinista:** copa = papila + 2 mm hacia el hilio; cáliz mayor = centroide de sus copas → Z=0 → empuje 3 mm al hilio; pelvis = centroide de los mayores → X=0, Z=0 → **punto medio con el hilio** (Y derivado ≈ **−27.2 mm**, no hardcodeado); uréter = hilio [0,−30,0].
- **Aristas orientadas papila→uréter** con `arista_tipo`: las 10 `papila→copa` son **junction DISCRETA** (tipo 1, escalón Bellini 0.2→1.5, no taper); el resto son taper (tipo 0); `pelvis→ureter` (4.5→1.5) codifica el **UPJ**.
- Campos para stitching en Capa 5: `papila_nodo_ref` (10, índices globales de 3c), `piramide_id`, `caliz_menor_id`, `caliz_mayor_id`, `nivel`, meta (`version, capa='4', n_calices_menores=10, n_calices_mayores=2, eje_polar='X', radios_anclados_json`, refs a capa0/capa3c).

### 2. Verificación previa
- Rango de `pyramid_apex` por eje: X **23.76 mm** (mayor) vs Y 2.78 vs Z 8.80 → **EJE_POLAR = X confirmado [OK]**.
- `papila_nodo_ref` coincide con `pyramid_apex`: desvío máx **7.18e-7 mm [OK]**.

### 3. Auditoría (4 bloques → `09_paper_vascular/auditoria_capa4_calicial.md`)
- **(1) Alcanzabilidad:** desde cada papila, siguiendo aristas, se llega al único nodo uréter: **10/10 [OK]**.
- **(2) Colisiones:** centerlines **no se cruzan** (dist mín no adyacente **3.168 mm** > 0) [OK]. Solape de luz con **radio interpolado** en el punto de aproximación (test físico correcto para taper, desglosado): **cruce entre cálices mayores DISTINTOS = 0 [OK]** (criterio de fallo). Residuales esperados (NOTA, no defecto): 2 solapes **mismo-mayor** (embudo de cálices menores a su mayor común, ~85 µm) y 20 con la **pelvis** (los cálices drenan a la cámara; su lumen de 4.5 mm envuelve las bocas de los infundíbulos → continuidad calicial→pelvis).
- **(3) Contención en el seno** (elipsoide `((x)/22)²+((y+34)/16)²+((z)/11)²≤1`): cáliz menor 10/10, cáliz mayor 2/2, pelvis 1/1, uréter 1/1 → **100 % dentro [OK]**. Las 10 `papila_junction` son **interfaz médula/seno, exentas**: valor del elipsoide = **1.000** (sobre la pared +Y del seno, esperado); polares **k=0 = 1.000, k=9 = 1.000**.
- **(4) Monotonía de radio por rama:** perfil papila→uréter **[0.2, 1.5, 2.0, 4.5, 1.5]** → creciente papila<copa<mayor<pelvis **[OK]**, única caída **pelvis→uréter (UPJ)** [OK]; el escalón Bellini→copa es junction discreta.

### Estado
**Capa 4 válida.** Sistema calicial alto en centerline+radio, convergente al hilio, radios absolutos, contenido en el seno, sin cruces inter-mayor. Frontera: **Capa 5 (voxelización)**, con stitching vía `papila_nodo_ref`. Nota de citas: los radios copa/infundíbulo/pelvis/uréter son placeholders morfométricos absolutos marcados `[CITA PENDIENTE]` en `radios_anclados_json` (Bellini heredado de 3c). NO se generó Capa 5.

---

## ENTRADA 016 — 8 de julio de 2026 — Capa 4: visualizador calicial + cierre de citas de radios (geometría intacta)
**Estado:** Dos tareas sobre la Capa 4 **ya generada** (entrada 015), sin regenerar geometría salvo el re-run determinista de la tarea B. No se tocaron capas 0–3c ni se generó Capa 5.

### 1. Visualizador de Capa 4 (`08_gemelo_digital/visualizador_maestro.py`)
- Añadido flag **`VER_CALICES = True`** (junto a los demás `VER_*`) y función **`draw_calices()`**, siguiendo el patrón `draw_colector`/`VER_COLECTOR` (flag por capa + función draw + carga perezosa + colección propia `Coleccion_Calices`). **No** se modificaron las funciones draw ni los flags de otras capas.
- `draw_calices()` lee `capa4_colector_alto.npz` y dibuja **cada arista como cilindro cónico (tapered)** entre `nodos[parent]` y `nodos[child]`, con `radio_inicio = radios[parent]` y `radio_fin = radios[child]` (centerline+radio, `primitive_cone_add` con radius1/radius2; **sin** primitivas de cámara — el gate es ver si la pelvis-como-cilindro basta a tamaño físico real, radios **absolutos** 0.2–4.5 mm sin escala cosmética).
- **Color por `nivel` (0..4):** 0 papila_junction gris · 1 cáliz menor naranja claro · 2 cáliz mayor naranja · 3 pelvis rojo · 4 uréter magenta (`COLOR_CALIZ`).
- **Aristas `arista_tipo == 1` (junction discreta Bellini→copa):** NO interpola el radio; renderiza el **escalón visible** = stub fino r=0.2 tope-a-tope con cilindro r=copa (dos cilindros de radio constante que abutan), para inspeccionar la discontinuidad de luz. Respeta `CORTE_ACTIVO`. (El visualizador corre aparte con Blender; se dejó listo, no se ejecutó aquí.)

### 2. Cierre de `radios_anclados` (citas morfométricas) — geometría byte-idéntica
- Editado el dict `radios_anclados` en `capa4_colector_alto.py`, reemplazando cada `[CITA PENDIENTE]` por (valor_mm, fuente): Bellini 0.200 (diám. 300–600 µm, heredado 3c) · copa 1.500 (cuello/infundíbulo del cáliz menor; menores 7–13 cup-shaped) · infundíbulo 2.000 (diám. infundibular ~4 mm modal, 60.3 % ≥4 mm; CTU n=1321, PMC10402955) · pelvis 4.500 (AP normal adulta <~10 mm, extremo alto-normal; umbrales de hidronefrosis 10–20 mm) · uréter 1.500 (luz ureteral ~3–4 mm distendida).
- **Re-run determinista** del generador. Verificado contra el `.npz` previo: **`nodos`, `radios`, `nivel`, `aristas`, `arista_tipo`, `papila_nodo_ref`, `piramide_id`, `caliz_menor_id`, `caliz_mayor_id` IDÉNTICOS byte-a-byte (tol 0)** [OK]. Las 4 auditorías dan el **mismo** resultado que la entrada 015 (alcanzabilidad 10/10; cruce inter-mayor 0; embudo mismo-mayor 2 [NOTA]; nesting pelvis 20 [NOTA]; contención 100 %; papilas val 1.000; perfil de radios [0.2, 1.5, 2.0, 4.5, 1.5]). Lo único que cambió es el meta `radios_anclados_json` (las citas; valores mm sin cambio).

### Estado
**Cerrado.** Visualizador con `VER_CALICES` listo para Blender; citas morfométricas incorporadas al `.npz` sin alterar la geometría (verificado byte-a-byte). Capas 0–3c intactas; Capa 5 no generada.

---

## ENTRADA 017 — 9 de julio de 2026 — Capa 5a: núcleo SDF + validación gruesa (200 µm smoke)
**Estado:** **Capa 5a generada y auditada.** Primer paso de la voxelización/impresión: representar cada capa (centerline+radio) como **SDF analítico** (unión de cápsulas cónicas, rounded-cone de Inigo Quilez), muestrearlo a 200 µm y mallar los lúmenes. Solo Capa 5a; **5b (matriz fina 100 µm + vóxel etiquetado) queda PENDIENTE**. Capas 0–4 intactas (timestamps 07-05/07-08, previos a esta sesión).

### 1. Esquemas hallados (PASO 0, normalizados con adaptador por-capa)
| capa | nodos | aristas | radios | adaptador |
|---|---|---|---|---|
| 3a arterial | 12016 | 12015 | **por-ARISTA** → nodo(child); raíz 214 µm | `adapt_edge_radios` |
| 3ab peritubular | 4290 pts | 0 | **ninguno** (solo puntos: bed capilar sub-resolución) | puntos edge-less |
| 3b venoso | 23799 | 23798 | **por-ARISTA** → nodo(child); raíz 524 µm | `adapt_edge_radios` |
| 3c colector | 2189 | 2179 | **por-NODO**; topología vía `parent` | `adapt_3c` |
| 4 calicial | 24 | 23 | **por-NODO**; aristas explícitas | `adapt_4` |

Hallazgo de esquema: 3a/3b guardan radios **por-arista** (no por-nodo); se derivó radio por-nodo = radio de la arista que alimenta al hijo, raíz = `radio_raiz`. **3ab es una nube de puntos sin aristas ni radios** (bed capilar): no aporta tubo imprimible.

### 2. Decisiones aplicadas
- **Seno relleno:** la matriz = elipsoide **completo** de Capa 0 (semiejes [55,30,18]); el seno no se resta (es donde va la pelvis).
- **Split 5a/5b:** 5a muestrea a **200 µm** (smoke) y malla lúmenes; 5b hará 100 µm + matriz fina + vóxel. El grid es solo muestreo del SDF (no se hace crecer geometría).
- **Poda + clamp a 400 µm diám. (200 µm r):** se eliminan aristas con **ambos** extremos sub-piso; los nodos retenidos por debajo se clampan a 200 µm (evita pinch). La microvasculatura sub-400 µm **NO se imprime**: queda en 3a/3ab/3b/3c como **objetivo de perfusión** (maduración capilar wet-lab).
- **STL:** writer binario propio (sin numpy-stl); scikit-image instalado en `env_biokidney` para `marching_cubes`.

### 3. Circuitos (PASO 1)
- **URINARIO = 3c ∪ 4:** tras fusión de nodos coincidentes (papilas 3c ↔ `papila_junction` de 4, desvío ~7e-7): **492 nodos, 2202 aristas, 1 sola componente conexa** [OK] (las 10 papilas fusionadas unen el bosque de 3c con el árbol calicial).
- **VASCULAR = 3a ∪ 3ab ∪ 3b:** **2 componentes-tubo** — arterial (6935 nodos, capa 3a) y venoso (14800 nodos, capa 3b) — **separados** (gap capilar de diseño, dato no defecto) + **4290 puntos aislados** (3ab, bed capilar edge-less).

### 4. Poda / clamp por circuito (PASO 2) — números emergentes
| circuito | aristas tot | podadas | retenidas | frac. superv. | clamp | reconx |
|---|---|---|---|---|---|---|
| urinaria | 2202 | 1320 | 882 | **40.1 %** | 10 | 0 |
| vascular_0_arterial | 12015 | 11986 | **29** | **0.2 %** | 1 | 0 |
| vascular_1_venoso | 23798 | 23299 | 499 | **2.1 %** | 77 | 0 |

**Dato clave:** el **árbol arterial es casi enteramente sub-400 µm** al scope reducido (solo 29/12015 aristas ≥ piso; raíz 214 µm). Es coherente con el desajuste de escala ya documentado (radio_raiz ~0.21 mm vs arteria renal real 2–3 mm): a esta densidad reducida el arterial imprimible es un stub sub-mm. El venoso retiene más (raíz 524 µm). El urinario retiene el 40 % (Capa 4 completa + troncos de 3c a 200 µm). 0 reconexiones (radios monótonos por rama → poda no desconecta del puerto).

### 5. Mallas de lumen (PASO 3, 200 µm) y auditorías (PASO 4)
- **STLs:** `capa5a_lumen_urinaria.stl` (69 486 tri), `capa5a_lumen_vascular_0_arterial.stl` (977 tri), `capa5a_lumen_vascular_1_venoso.stl` (19 729 tri). Meta en `capa5a_meta.npz`. Reporte `09_paper_vascular/auditoria_capa5a.md`.
- **(1) Watertight** (= cierre 2-manifold: manifold ∧ 0 bordes ∧ 0 non-manifold): **los 3 [OK]**. Urinaria género 1 (asas por fusión de cálices anidados a 200 µm; 5b/100 µm lo reduce); arterial/venoso género 0 (euler 2).
- **(2) Estanqueidad:** vóxeles dentro de AMBOS = **1118** (urin&arterial 100, urin&venoso 1018) — **ver DATOS EMERGENTES**. Dist mínima superficie-superficie urinario↔vascular **0.025 mm**; gap arterial↔venoso **1.546 mm**.
- **(3) Conectividad:** urinario alcanza [0,-30,0]; arterial su raíz [0,-30,0]; venoso [0,-30,3] (nodo a 0.489 mm) → **3/3 [OK]**.
- **(4) Contención en elipsoide:** urinario 0/408, arterial 0/18, venoso **1/316** (el nodo-puerto, val 1.033 — salida por diseño). [OK]
- **(5) Volúmenes:** matriz (elipsoide completo) **124 407 mm³**; lumen urinario 878.3, arterial 1.45, venoso 43.2 → **total 922.9 mm³**; **porosidad 0.742 %**.
- **Tiempo 21.9 s · pico de memoria 287 MB** (i5/16 GB, muestreo por sub-box por cono).

### DATOS EMERGENTES (no anticipados)
1. **Solape urinario↔vascular (1118 vóxeles, 0.025 mm):** NO es cruce topológico (los centerlines siguen disjuntos), sino artefacto de **clamp-a-piso**: micro-vasos venosos sub-400 µm se inflan a 400 µm diám. y al muestrear a 200 µm su lumen fusiona con el conducto colector en la médula (interdigitación vasa-recta / colector, anatómicamente adyacentes). Es **tensión de escala del modelo reducido**, a resolver en 5b (100 µm + política de clamp / separación inter-circuito ≥ 400 µm por diseño). Dominado por urin↔venoso (1018) frente a urin↔arterial (100).
2. **Arterial imprimible ≈ stub sub-mm** (29 aristas): confirma que, al scope reducido, el árbol arterial vive casi todo bajo el piso de impresión; su rol es objetivo de perfusión, no canal impreso. Sugiere revisar el calibre arterial (radio_raiz) si se quiere arterial imprimible.
3. **Urinaria género 1:** los cálices menores del mismo mayor (solape 85 µm de la Entrada 015) fusionan a 200 µm creando un asa; esperado, se reduce a 100 µm.

### Estado
**Capa 5a cerrada.** Núcleo SDF muestreado, poda/clamp aplicados, 3 lúmenes mallados y validados (watertight 2-manifold), auditorías corridas con dos hallazgos emergentes documentados. **5b (matriz fina 100 µm + vóxel etiquetado) PENDIENTE.** Capas 0–4 intactas.

---

## ENTRADA 018 — 9 de julio de 2026 — Diagnóstico de la falla de estanqueidad de 5a (CORRIGE la hipótesis de la 017)
**Estado:** Diagnóstico de **solo lectura** de los 1118 vóxeles compartidos urinario↔vascular de la auditoría 5a. Reutiliza los centerlines/radios ya normalizados (mismo adaptador y clamp de 5a), sin re-mallar, sin modificar 5a ni generar 5b. **Corrige la hipótesis de la Entrada 017.** Capas 0–5a intactas (timestamps 07-05…07-09; los `.npz`/STL de 5a de las 18:23 no se tocaron).

### Corrección a la Entrada 017
La 017 atribuyó el solape a **inflación de clamp** (micro-vasos sub-piso inflados a 400 µm diám.) + muestreo 200 µm, como "tensión de escala". **El diagnóstico lo desmiente:** de los 1118 vóxeles, **1115 (99.7 %) son de geometría NATIVA** (ambos lados con radio original ≥ 200 µm); solo **3** tocan un lado clampeado. No es artefacto de clamp: es **colisión macro real en el hilio**.

### Tabla de clasificación (nivel urinario × origen vascular)
| nivel urinario | vascular | ambos nativos | ≥1 clampeado |
|---|---|---|---|
| pelvis | venoso | 825 | 0 |
| ureter | venoso | 134 | 0 |
| pelvis | arterial | 69 | 0 |
| caliz_mayor | venoso | 34 | 0 |
| ureter | arterial | 31 | 0 |
| colector_3c | venoso | 22 | 3 |

Domina **pelvis↔venoso** (825) y **ureter↔venoso** (134): las estructuras urinarias **macro** (pelvis r=4.5 mm, uréter r=1.5 mm) solapan el **tronco venoso** (y el stub arterial) que convergen al mismo hilio.

### Pared inter-circuito NATIVA (radios originales) — el test que decide
| pared nativa | vóxeles |
|---|---|
| <0 (solape real) | **1067** |
| 0–100 µm | 0 |
| 100–400 µm | 51 |
| ≥400 µm (sub-muestreo) | **0** |

**Pared mínima nativa = −4347 µm** (−4.3 mm de interpenetración): la pelvis (r=4.5 mm), centrada en [0,−27.2,0], **engulle** el tronco venoso; de hecho el puerto venoso [0,−30,3] queda a 4.1 mm del centro de la pelvis, **dentro** de su lumen de 4.5 mm. Peor contacto en [−1.12,−27.56,0] (cáliz_mayor/pelvis vs venoso).

### Zoom hilio
Pelvis/uréter vs tronco venoso en Y∈[−31,−27] (3×6 segmentos): pared lumen-lumen mínima **−4399 µm (<400 µm)** → **HAY componente de hilio**. urin ~[1.95,−25.81,0], venoso ~[−0.6,−27.28,−0.16].

### Género / handle de la malla urinaria
2 fusiones urinarias no-adyacentes (lumen solapa): `papila_junction/caliz_menor × caliz_menor/caliz_mayor` a d=3.168, **solape 85 µm**, en [∓6.08,−20.92,±4.05]. **Coinciden con la fusión co-mayor** (cálices menores del mismo mayor, solape 85 µm de la Entrada 015) → confirma que el género 1 de la malla urinaria viene de ese anidamiento, no del hilio.

### Números que quedaron cortados en 5a
- **Porosidad: 0.7418 %.**
- Lumen por circuito: urinaria **878.26**, arterial **1.45**, venoso **43.18** mm³.
- **Gap arterial↔venoso (nodos): 1.546 mm.**

### VEREDICTO
- **Sub-muestreo** (pared nativa ≥400 µm, se separaría a 100 µm): **0/1118 (0.0 %)**.
- **Colisión REAL** (pared nativa <400 µm): **1118/1118 (100.0 %)** — de la cual 1067 es interpenetración franca (<0).
- **Componente de hilio: SÍ** (pelvis/uréter ↔ venoso, pared −4399 µm).
- **Origen:** convergencia macro en el hilio: la **pelvis (r=4.5 mm)** y el uréter ocupan el mismo espacio que el tronco venoso/arterial (puertos urinario [0,−30,0] y venoso [0,−30,3] a 3 mm). **NO es fixable subiendo a 100 µm**: requiere **fix de geometría** (separar la pelvis del hilio vascular: offset/menor radio de pelvis, o enrutar el tronco venoso alrededor, o desplazar el puerto venoso). El asa de género (85 µm co-mayor) es aparte y sí es de resolución.

### Estado
**Diagnóstico cerrado.** La falla de estanqueidad de 5a es **colisión geométrica real en el hilio (pelvis vs tronco venoso)**, no artefacto de clamp/resolución (corrige la 017). Requiere revisión de diseño del acople pelvis↔hilio antes de 5b. Nada de 5b ni de geometría modificado.

---

## ENTRADA 019 — 9 de julio de 2026 — Fix de estanqueidad 5a: FIX1 aplicado, reruteo en FALLBACK
**Estado:** Intento de cierre del blocker de estanqueidad de la Entrada 018. **FIX1 (pelvis aplanada) aplicado y verificado; FIX2/FIX3 (reruteo venoso/arterial de dirección-fija) cayeron en FALLBACK** por obstrucción geométrica real. No se forzó, no se hizo carve, no se generó 5b. Capas 0–3c intactas (timestamps 07-05); solo cambió la Capa 4 (FIX1, campos nuevos) — 3a/3b/3c NO se reescribieron (los reruteos eran transformaciones locales de Capa 5, no aplicadas por el FALLBACK).

### FIX1 — pelvis aplanada (Capa 4) [APLICADO, verificado]
- `capa4_colector_alto.py`: la pelvis pasa de esfera r=4.5 a **elipsoide [4.5, 4.5, Z_pelvis]** con **Z_pelvis = 2.0** (aplanamiento AP). Nuevos campos en el `.npz`: `tipo_primitiva` (M,) = 'tubo' salvo pelvis='elipsoide', `semieje_z` (M,) = radio salvo pelvis=2.0, escalar `z_pelvis`. El centro de la pelvis NO se movió.
- **No-regresión:** los 9 campos existentes (`nodos, radios, nivel, aristas, arista_tipo, papila_nodo_ref, piramide_id, caliz_menor_id, caliz_mayor_id`) **byte-idénticos (tol 0)**; las 4 auditorías de Capa 4 **idénticas a la Entrada 016**. Contención del elipsoide: uréter dentro (val 0.386), caliz_mayor a 4.80 mm (val 1.139) conecta por el infundíbulo (r=2.0) sin pinch, igual que en 5a con la esfera.

### FIX2/FIX3 — reruteo local (Capa 5) [FALLBACK]
- Keep-out = pelvis elipsoide + infundíbulos + uréter, dilatado por MIN_WALL=0.400 mm + radio venoso/arterial. Solver: empujar los nodos terminales venosos a +Z (anterior) y los arteriales a −Z (posterior), con desplazamiento cumulativo capado a 5 mm y suavizado ligero; puertos no más de 5 mm.
- **Medición por-nodo (desde original):** venoso 9 nodos, máx **4.99 mm** (≤5); arterial 5 nodos, máx **4.53 mm** (≤5); puertos ≤2.23 mm. Individualmente parecía factible.
- **Pero el reruteo ACOPLADO no cierra:** tras aplicar (cap 5 mm + suavizado + conectividad), **restan 4 nodos venosos y 7 arteriales violando MIN_WALL**. **Causa geométrica (STRADDLE):** los 4 venosos residuales son **todos POSTERIORES (Z≈−2.7…−3.3)** — empujarlos a +Z los haría cruzar toda la pelvis; los 7 arteriales residuales son **todos ANTERIORES (Z≈+2.5…+3.5)** — empujarlos a −Z también cruzaría la pelvis. Cada árbol tiene ramas en el lado "equivocado" de la pelvis, así que **una dirección fija no puede separarlos** sin atravesar la cámara (>5 mm). Ejemplo: nodo arterial 8 en [0.52,−27.8,+3.54] necesitaría cruzar 1124 µm de pared faltante yendo a −Z.

### Estanqueidad antes/después
- **Antes:** 1118 vóxeles compartidos (5a / Entrada 018).
- **Después:** el reruteo de dirección-fija **NO se aplicó a la geometría** (FALLBACK), así que la estanqueidad sigue sin resolverse (queda 1118 pendiente de decisión). FIX1 por sí solo (aplanar la pelvis a Z=2.0) reduce el keep-out pero no separa las ramas vasculares que pasan a Z≈0 dentro del footprint XY r=4.5.

### Opciones propuestas (NO apliqué ninguna; decisión del usuario)
1. **Más aplanamiento** de la pelvis (Z_pelvis < 2.0) — reduce el keep-out en Z.
2. **Mover el nodo pelvis en −Y** (hundirlo hacia el seno) — lo aleja del cruce vascular.
3. **Reruteo POR-GRADIENTE** (cada nodo se aparta de la superficie más cercana, no dirección fija) — resuelve el straddle, pero cambia la política del solver.
4. **Carve del keep-out** en la matriz — NO por defecto (preferencia explícita de decidirlo a mano).

### DATOS EMERGENTES
- El conflicto pelvis↔vascular en el hilio es **bilateral en Z**: no es que la vena esté "de un lado", sino que tanto arterial como venoso tienen ramas anteriores y posteriores rodeando la pelvis. Confirma que el acople pelvis↔hilio necesita un rediseño (no un empujón), coherente con el veredicto de la Entrada 018.
- FIX1 queda como avance parcial válido y reusable (contrato de primitiva de pelvis) para cualquiera de las opciones.

### Estado
**Fix PARCIAL: FIX1 (pelvis elipsoide) aplicado y verificado; reruteo en FALLBACK a la espera de decisión de diseño.** No se forzó reruteo, no se hizo carve, no se generó 5b. STL del fix (del run intermedio) eliminados por no representar geometría válida; los STL/meta de 5a original quedan intactos. Reporte en `09_paper_vascular/auditoria_5a_fix.md`.

---

## ENTRADA 020 — 9 de julio de 2026 — Fix estanqueidad 5a v2: reruteo POR-GRADIENTE (también FALLBACK, diagnóstico más profundo)
**Estado:** Segundo intento de cierre del blocker de estanqueidad, con la **política por-gradiente** que reemplaza la de dirección-fija (Entrada 019). Cada nodo en conflicto se aparta de la superficie de otro-circuito más cercana siguiendo `+grad(SDF_keepout)` (no una dirección global), para disolver el straddle en Z. Mantiene **FIX1 sin cambios** (pelvis elipsoide Z_pelvis=2.0). Reruteo LOCAL en Capa 5 — **transformación de puntas cerca del hilio, NO cambio de topología, NO reescribe 3a/3b/3c**. **También cayó en FALLBACK**, pero revela la causa raíz con más precisión. Capas 0, 3a, 3b, 3c intactas (07-05); Capa 4 solo FIX1 (07-09). No se generó 5b.

### Política nueva (gradiente multi-superficie)
- Keep-out por circuito = SDF de "todo lo no-propio", dilatado MIN_WALL=0.400 mm: para venosos = unión(lumen urinario [pelvis elipsoide + uréter + cálices], lumen arterial); para arteriales = unión(lumen urinario, lumen venoso). Cada nodo se verifica contra **todas** las superficies, no solo la pelvis.
- Solver iterativo (grad por diferencias finitas), paso adaptativo, con restricción de conectividad (arista al padre sin estirar >1.5× ni invertir) y `MAX_DESPL=5.0 mm/nodo`. Puertos libres dentro de MAX_DESPL.

### Nodos movidos / desplazamientos
- **Venoso:** 31 nodos movidos, desplazamiento máx **1.32 mm** (medio 0.03), puerto 0.00 mm; **8 nodos residuales** violando.
- **Arterial:** 5 nodos movidos, máx **0.78 mm**, puerto 0.41 mm; **4 residuales**.
- Solver: 8 iteraciones, 262 s.

### Estanqueidad antes/después
- **Antes: 1118.** **Después: sin resolver** — FALLBACK antes del re-mallado; el reruteo NO se aplicó a la geometría final (nada que promover). No hay tabla "después".

### FALLBACK — diagnóstico raíz (por qué el gradiente se estanca)
Los desplazamientos son **diminutos** (≤1.3 mm) pese a déficits enormes, porque la **restricción de conectividad** bloquea el escape: los nodos en conflicto forman una **cadena del tronco que ATRAVIESA la pelvis**. Residuales venosos: nodos 4→8 son consecutivos y recorren la cámara en Z (Z=+1.65 → +0.55 → −0.52 → −1.52 → −2.47), con paredes de **−788, −2041, −2072, −985, +54 µm** (los del centro a −2 mm DENTRO del lumen de la pelvis). Ningún nodo puede salir porque sus vecinos de cadena también están atrapados dentro (mover uno estiraría >1.5× la arista). 
- **Hallazgo adicional decisivo:** el **puerto/raíz arterial [0,−30,0] COINCIDE con el puerto del uréter** (mismo punto de diseño): el nodo arterial 0 tiene pared **−1513 µm** (1.5 mm dentro del lumen del uréter). El arterial no puede despejar el uréter sin mover un puerto que está fijado al mismo punto.
- Es decir: **no es proximidad de puntas** (nudge-able), sino **conflicto de ruteo/topología en el hilio** — el tronco venoso pasa por el interior de la cámara de la pelvis y el puerto arterial está encima del uréter. Ninguna deformación terminal acotada (conectividad + 5 mm) lo resuelve.

### Convergencia de los dos intentos
Los solvers de **dirección-fija (019)** y **gradiente (020)** caen ambos en FALLBACK y coinciden: el acople pelvis↔hilio↔puertos vasculares necesita **rediseño**, no un empujón de puntas. El gradiente añade la evidencia de que el tronco venoso **atraviesa** la pelvis y de la **coincidencia puerto-arterial = uréter**.

### DATOS EMERGENTES
- El puerto arterial [0,−30,0] y el puerto urinario/uréter [0,−30,0] son **el mismo punto** por diseño de Capa 0/3a → colisión inevitable a nivel de puerto (independiente de la pelvis). Separar los puertos es condición necesaria.
- El tronco venoso terminal no bordea la pelvis: la **cruza** (cadena de nodos a ambos lados de Z dentro del elipsoide). Reruteo terminal acotado es insuficiente; requiere re-ruteo del tronco (o mover la pelvis, o hundir el hilio).

### Opciones (NO apliqué ninguna; decisión del usuario: "carve local vs revisar el hilio")
1. **Separar los puertos** (mover el puerto venoso y/o arterial en Z/Y respecto al uréter) a nivel de diseño de Capa 0/3a/3b — resuelve la coincidencia de puertos.
2. **Mover la pelvis** (nodo pelvis en −Y hacia el seno, o más aplanamiento) para que el tronco venoso no la atraviese.
3. **Carve local** del keep-out en la matriz (no por defecto).
4. Re-ruteo del **tronco** venoso (no solo puntas) — implicaría revisar 3b, fuera del alcance "local Capa 5".

### Estado
**Fix v2 en FALLBACK (gradiente).** FIX1 (pelvis elipsoide) permanece aplicado y verificado. El blocker de estanqueidad es un **conflicto de ruteo/puertos en el hilio** (tronco venoso atraviesa la pelvis; puerto arterial = uréter), confirmado por dos solvers independientes. Requiere decisión de diseño. No se forzó, no se hizo carve, no se generó 5b. Reporte en `09_paper_vascular/auditoria_5a_fix2.md`.

---

## ENTRADA 021 — 9 de julio de 2026 — Fix estanqueidad 5a v3: pedículos hilares segregados (determinista) — también FALLBACK
**Estado:** Tercer intento del blocker de estanqueidad, con **construcción determinista de pedículos hilares segregados** (reemplaza los dos solvers previos). Anclaje anatómico: orden anteroposterior **vena–arteria–pelvis** (eje AP = Z; +Z = anterior). Mantiene **FIX1 sin cambios** (pelvis elipsoide Z_semi=2.0); pelvis y uréter NO se mueven (referencia posterior). Reruteo LOCAL en Capa 5 (síntesis de pedículo) — NO reabre 3a/3b/3c, NO cambia Capa 4, NO genera 5b. **También cayó en FALLBACK (PASO 2)**, ahora por una causa estructural precisa. Capas 0,3a,3b,3c intactas (07-05); Capa 4 solo FIX1 (07-09).

### Anclaje anatómico
El pedículo hilar idealiza el orden AP de libro: **vena (anterior) → arteria (media) → pelvis (posterior)**, offset anterior del plano medio AP. Es un **patrón de referencia con variabilidad documentada** (el orden vena/arteria/pelvis no es invariante; ~20–80 % según lado/nivel) → fila de correspondencia = corredor hilar vena-arteria-pelvis idealizado `[CITA PENDIENTE: variabilidad AP del pedículo renal]`.

### PASO 1 — corredores Z (determinista, radios reales del tronco)
- r_art = 214 µm (raíz arterial), r_ven = 524 µm (puerto venoso).
- **Z_art = Z_semi + r_art + MIN_WALL + MARGIN = 2.814 mm** (arteria, 1er corredor anterior a la pelvis).
- **Z_ven = Z_art + r_art + r_ven + MIN_WALL + MARGIN = 4.152 mm** (vena, anterior a la arteria).
- La vena despeja la pelvis directo: 4.152 − 2.0 − 0.524 = **1.629 mm ≥ 0.400 [OK]**.
- Puertos nuevos propuestos: ARTERIA [0,−30,2.814], VENA [0,−30,4.152]; uréter [0,−30,0] sin cambios.

### PASO 2 — conflicto (keepout = pelvis elipsoide + uréter + cálices, Capa 4) y verificación
- **Arterial:** 5 nodos, Y[−30,−28.8], puerto EN conflicto, 1 salida, 0 nodos-rama → **[OK] cadena terminal**.
- **Venoso:** 10 nodos, Y[−29.1,−22.0], puerto **NO en conflicto**, **6 salidas**, **1 nodo-rama (grado>2)** → **[FALLA cadena]**. Es un **plexo peri-hilar ramificado**, no una punta única: snipearlo y reconectar por UN pedículo desconectaría 5 ramas. Peor nodo venoso 5 en [−0.63,−27.53,0.14]: pared **−3421 µm** (3.4 mm dentro de la pelvis).

### PASO 3 — construcción (informativa) y por qué no cierra
- **Arterial (cadena):** el ancestro más anterior del stub retenido está en [0.45,−24.12,8.22], **Y=−24.1 (≤ −22.7)** — el stub arterial imprimible está **entero dentro de la banda Y de la pelvis**, no hay parénquima anterior donde anclar el pedículo. El pedículo ancestro→W→puerto roza la pelvis: pared mínima **−166 µm** → **no verifica**.
- **Venoso:** ni siquiera pasa PASO 2 (plexo ramificado).

### Estanqueidad antes/después
- **Antes: 1118.** **Después: sin resolver** (FALLBACK en PASO 2, antes del re-mallado; nada aplicado, sin STL de fix3).

### Convergencia de los TRES intentos
Dirección-fija (019), gradiente (020) y pedículo determinista (021) caen todos en FALLBACK y **coinciden en la causa estructural**: el acople hilar no es fixable con transformaciones locales de Capa 5 porque **la geometría vascular de origen (3a/3b) no lo soporta**:
1. El **árbol venoso ramifica inmediatamente en el hilio** (plexo peri-pélvico de 6 ramas), no es un tronco único → ninguna segregación de pedículo-único lo captura.
2. El **arterial imprimible es un stub sub-mm peri-pelvico** (radio raíz 214 µm; casi todo el arterial es sub-piso) sin parénquima anterior donde anclar.
3. Los **puertos arterial y urinario coinciden** ([0,−30,0]) por diseño de Capa 0/3a.

### DATOS EMERGENTES
- La segregación hilar anatómica (corredores Z) es **geométricamente válida y computada** (Z_art/Z_ven resueltos), pero **inaplicable a esta topología vascular**: el problema no es el corredor sino que el árbol venoso no tiene un tronco hilar único y el arterial imprimible no llega al parénquima anterior.
- Explícito: los pedículos hilares serían geometría **SINTETIZADA por Capa 5** (no derivada de 3a/3b; los `.npz` de origen NO se reescriben), misma categoría que la pelvis-primitiva (FIX1); el arterial se mantendría como boca de entrada corta imprimible. Pero la síntesis de pedículo-único no aplica aquí.

### Opciones (NO apliqué ninguna; decisión del usuario)
1. **Tronco venoso único en el hilio** (revisar 3b: forzar convergencia a un solo tronco peri-hilar antes del port) — fuera del alcance "local Capa 5".
2. **Multi-pedículo venoso** (un corredor Z por cada rama en conflicto) — extiende la construcción a N pedículos.
3. **Separar los puertos** arterial/venoso del uréter a nivel de diseño (Capa 0/3a/3b).
4. **Mover la pelvis en −Y** (hundirla hacia el seno) para reducir el solape peri-pélvico venoso; o **carve local** (no por defecto).

### Estado
**Fix v3 en FALLBACK (PASO 2).** Los tres enfoques (fijo, gradiente, pedículo determinista) confirman que el blocker es de **diseño hilar** (árbol venoso ramificado + arterial sub-mm + puertos coincidentes), no resoluble con reruteo local de Capa 5 sobre la geometría 3a/3b actual. FIX1 (pelvis elipsoide) permanece aplicado y verificado. No se forzó, no se hizo carve, no se generó 5b. Reporte en `09_paper_vascular/auditoria_5a_fix3.md`.

---

## ENTRADA 022 — 10 de julio de 2026 — Diagnóstico de holgura de la pelvis: el fix NO es la pelvis, es el plexo venoso
**Estado:** Diagnóstico **comprehensivo de solo lectura**: barrido de posición del centro de la pelvis (X=0, Y-Z) buscando un centro que cierre **simultáneamente** pared vs venoso/arterial, contención en seno, intrarrenalidad, infundíbulos mayores sin pinza/colisión, uréter y corredor anterior. Reencuadre: la pelvis debe caber en el **hueco inter-tributaria** del hilio. Reusa SDF/adaptador/clamp/FIX1. **Nada modificado** (capas 0–4 intactas, 07-05…07-09; solo el script diag + su md). **5b no generado.**

### D1 — mapa del hilio (occupancy vascular peri-pélvico)
- En Y[−33,−20], Z[−9,9], X[−7,7]: **venoso 32 segmentos** — **27 posteriores (−Z)**, 5 anteriores (+Z), 1 en Z~0. Arterial 29 segmentos. El **plexo venoso es dominante-posterior y denso** alrededor del hilio.
- Hueco inter-tributaria (X=0): clearance máx 10.83 mm hacia el borde anterior (Y−20.2, Z9), pero el vacío útil está fragmentado por el plexo.

### D2 — barrido de centro de pelvis (525 centros, X=0)
- **17 centros VIABLES** (pared_ven≥0.4 Y pared_art≥0.4 Y seno≤1): **el cuerpo de la pelvis SÍ cabe**. Máscara Y[−30,−23], Z[−6,2.5]. Mejor intrarrenal [0,−23,2] (100 % intrarrenal, pared_ven 0.54); mejores paredes en Z≈−5/−6 (posterior, pero intrarrenal 31–47 %).

### D3–D5 — infundíbulos / uréter / corredor (sobre los 15 mejores viables)
- **D3 infundíbulos: 15/15 FALLAN** (inf_wall negativo, −0.09 a −2.29 mm): el infundíbulo cáliz_mayor(Z~0)→pelvis cruza el plexo venoso.
- **D4 uréter: 15/15 FALLAN** (uréter_wall −0.96 mm, constante): pelvis→hilio [0,−30,0] cruza el venoso cerca del hilio.
- **D5 corredor anterior: 15/15 OK** (los puertos anteriores vena/arteria despejan la pelvis reposicionada).

### D6 — VEREDICTO: **el fix NO es la pelvis (ni posición ni tamaño)**
- El cuerpo de la pelvis cabe (D2), pero los **infundíbulos (cálices en Z~0) y el uréter (hilio en Z~0) están ANCLADOS en Z~0** y deben cruzar el **corredor obligado Z~0** que ocupa el **plexo venoso peri-hilar**.
- **Barrido completo de posiciones:** NINGUNA tiene infundíbulos+uréter limpios; la mejor posible ([0,−23,1]) da pared **−964 µm** (faltan 1364 µm).
- **Independiente del tamaño:** infundíbulos y uréter usan el centro/cálices/hilio, no los semiejes de la pelvis → **reducir o mover la pelvis NO los despeja**. No es POSICION ni POSICION+TAMAÑO.
- **BINDING = el plexo venoso peri-hilar** (32 segmentos, 27 posteriores) ocupando el corredor Z~0. **El fix es segregar/mover el plexo venoso (revisar 3b)**, no la pelvis. Confirma **cuantitativamente** el blocker de las Entradas 019–021 (árbol venoso ramificado en el hilio).

### Reencuadre y dato emergente
- El problema nunca fue el cuerpo de la pelvis (cabe en 17 posiciones) sino las **conexiones ancladas en Z~0** (cálices, hilio) forzadas a atravesar el plexo venoso. El rango Y intrarrenal es estrecho (borde de parénquima −27.7): los centros más intrarrenales (Y≈−23) tienen la pared venosa más ajustada.
- Corresponde: el hilio real tiene los tres elementos (vena/arteria/pelvis-uréter) segregados AP; aquí el **venoso no está segregado** (plexo que llena el hilio), y por eso ninguna colocación de la pelvis resuelve. La palanca es la **arquitectura venosa hilar** (3b), no la pelvis (Capa 4).

### Estado
**Diagnóstico cerrado.** Veredicto: **el blocker de estanqueidad NO se resuelve reposicionando ni redimensionando la pelvis** — el binding es el plexo venoso peri-hilar en el corredor Z~0 obligado de infundíbulos/uréter. La decisión de diseño (segregar el tronco venoso en el hilio) queda del lado de 3b. Reporte en `09_paper_vascular/diagnostico_holgura_pelvis.md`. Nada modificado; 5b pendiente.

---

## ENTRADA 023 — 11 de julio de 2026 — Adaptador gemelo -> CSV de segmentos (esquema v7) para consumo por los simuladores
**Estado:** **ADAPTADOR creado y verificado.** Nuevo exportador `08_gemelo_digital/exportar_cco_v8_csv.py` que lee las capas validadas del gemelo (.npz) y escribe un CSV con el **mismo esquema de columnas que el v7**, de modo que los simuladores (oxígeno, glomerular) consuman el gemelo **sin modificar su código**. Capas .npz y simuladores **intactos**; v7 y v8 existente **conservados**. Salida: `02_vascular_cco/arbol_vascular_cco_v8_gemelo.csv`.

### PASO 0 — Contrato de columnas (leído del v7 real, no asumido)
- **Header v7 exacto (10 columnas):** `id,sistema,nivel,x1_mm,y1_mm,z1_mm,x2_mm,y2_mm,z2_mm,radio_um`.
- **NO existen** columnas de presión ni de terminal en el v7 (la tarea las suponía). Decisión: **no inventar presiones**; se reproduce el esquema de 10 columnas tal cual.
- **Oxígeno** (`cargar_arbol_cco`): solo requiere coords + `radio_um` + `sistema`. La **presión se DERIVA** del string `sistema` (art->P_ART, ven->P_VEN, resto->punto medio); **no se lee del CSV**. Por eso el esquema v7 basta.
- **Glomerular**: busca `es_terminal`/`presion_salida` pero **degrada con gracia** si faltan (toma los últimos 1000 segmentos, P=45 mmHg) — comportamiento idéntico al que ya tiene hoy sobre el v7. No requiere columnas nuevas.
- **Nota para Carlos (mejora opcional, NO aplicada):** las capas SÍ contienen arrays `terminales` (5699 art / 11212 ven / 903 col). Se podría poblar un `es_terminal` real para que el glomerular use terminales verdaderos en vez del fallback "últimos 1000". Queda fuera del esquema v7 estricto; pendiente de decisión.

### PASO 1 — Mapeo capa -> segmentos (una arista = una fila)
- `capa3a_arterial.npz`     (nodos/aristas/radios) -> sistema **art** -> **12015** segmentos
- `capa3b_venoso.npz`       (nodos/aristas/radios) -> sistema **ven** -> **23798** segmentos
- `capa3c_colector.npz`     (nodos/parent/radios, parent-pointer tree) -> sistema **col** -> **2179** segmentos
- `capa4_colector_alto.npz` (nodos/aristas/radios) -> sistema **col** -> **23** segmentos
- **TOTAL 38015 segmentos.** Comparación: v7 = 1448; cifra del paper para v8 = 1902.
- **Discrepancia reportada y decidida por Carlos:** el grafo del gemelo es el **grafo CCO crudo** (un nodo por paso de crecimiento), ~20x más denso que el v7 o la cifra del paper. La cifra 1902 del paper corresponde a **otro pipeline** (`generador_cco_v8.py`), no a las capas. **Decisión de Carlos: resolución completa (~38k), una arista = una fila** (no coarsening). El coarsening a ~1902 sería una tarea aparte (colapso de cadenas / poda Strahler), no un exportador.

### PASO 2 — Unidades y convenciones (verificadas contra el v7)
- **Coordenadas:** las capas ya están en **mm** (rangos art x[-53.9,54.1], y[-30,29.2], z[-17.6,17.5]) -> se escriben en mm en las columnas `*_mm`, **SIN conversión** (el oxígeno hace *0.1 mm->cm al cargar).
- **Radios:** las capas guardan radio en **mm** (radio_raiz art 0.214 mm, ven 0.524 mm). Se convierten a **µm (x1000)** para `radio_um`. Verificación: raíz venosa 0.524 mm -> 524 µm, coherente con el máx venoso del v7 (542 µm). El oxígeno hace *1e-4 µm->cm al cargar.
- **`sistema`:** se usan los **códigos cortos `art`/`ven`/`col`** idénticos al v7 (la tarea proponía 'arterial'/'venoso'/'colector', pero el v7 real y el derivador de presión del oxígeno usan los cortos). Prioricé el dato real del v7.
- **`nivel`:** profundidad topológica (saltos) desde la raíz del sistema vía BFS — dato real derivado del grafo, no inventado; el loader del oxígeno lo ignora.

### PASO 3 — Escritura + verificación (recarga)
- Escrito `arbol_vascular_cco_v8_gemelo.csv` (38015 segmentos). **NO se sobrescribió el v7** ni el `arbol_vascular_cco_v8.csv` existente.
- Recarga con pandas OK. Header idéntico al v7. Partición por sistema cuadra: art 12015, ven 23798, col 2202.
- **Rangos físicos:** coords mm [-53.923, 54.135]; `radio_um` art [12.00, 214.35], ven [23.40, 523.73], col [15.00, 4500.00] (el máx 4500 µm = cáliz/pelvis de capa4). Todos físicos.

### PASO 4 — Prueba de consumo (sin modificar simuladores, sin correr el solver)
- Se ejecutó **solo** `cargar_arbol_cco` del simulador de oxígeno apuntando al `_gemelo.csv`: **CARGA OK**, 38015 segmentos, partición correcta, presiones derivadas fisiológicas (art 40 / ven 20 / col 30 mmHg), grid bbox renal (±5.44 x, ±3.06 y, ±1.81 cm). No se corrió el solver Fick.

### Colisión de nombre — decisión de Carlos
- `arbol_vascular_cco_v8.csv` **ya existía** (26-may, 1902 segmentos, generado por `generador_cco_v8.py`, con columnas extra presión/terminal; es el que referencia el preprint). **Decisión de Carlos: mantener el export del gemelo como `arbol_vascular_cco_v8_gemelo.csv` (no destructivo)**; el v8 del paper queda intacto. Para usar el árbol nuevo, se apunta el simulador explícitamente al `_gemelo.csv`.

### Integridad verificada (timestamps)
- Capas .npz sin tocar: capa3a (07-05 18:21), capa3b (07-05 19:34), capa3c (07-05 18:23), capa4 (07-09 19:20).
- Simuladores sin tocar: oxígeno (03-26), glomerular (03-26).
- v7 (03-22) y v8 existente (05-26) conservados. Nuevo: exportador + `_gemelo.csv` (07-11).

### Estado
**Adaptador entregado y validado end-to-end en la carga.** El gemelo validado (Capas 0-4) es ahora consumible por los simuladores vía `arbol_vascular_cco_v8_gemelo.csv` con esquema v7 exacto, sin tocar capas ni simuladores. Pendientes opcionales (decisión de Carlos): (a) poblar `es_terminal` real para el glomerular; (b) coarsening a ~1902 si se quiere paridad con el paper.

---

> **Nota de consolidación (11-jul-2026).** Las Entradas 024-027 registran trabajo del **mismo día 11-jul** que es **cronológicamente ANTERIOR** a la Entrada 023 (adaptador), pero se asientan después por orden de escritura. Se numeran consecutivas sin renumerar entradas previas. Auditado: los ítems de la revisión de simuladores, el marco de honestidad, el preprint v3 y la publicación en Zenodo v3 NO tenían entrada formal; el adaptador SÍ (Entrada 023, cerrado). No se tocó código, capas ni el preprint en esta consolidación: solo BITACORA.md.

---

## ENTRADA 024 — 11 de julio de 2026 — Revisión crítica de los 6 simuladores + optimizador Co-SWIFT: clasificación en 3 niveles

**Estado:** **CERRADO (auditoría).** Revisión módulo a módulo de los seis simuladores funcionales y del optimizador Co-SWIFT, con clasificación en **3 niveles de validez epistémica**. Hallazgo arquitectónico central y fallas concretas confirmadas por inspección de código. **Solo lectura y análisis: no se modificó ningún simulador.**

### Hallazgo arquitectónico
- El multisimulador consume un **CSV de segmentos** (`arbol_vascular_cco_v7.csv`), **NO una malla ni las capas nuevas del gemelo**. La fabricación imprimible (voxelización, Capa 5) es fase posterior y separable. Esto motiva luego el adaptador de la Entrada 023 (exportar el gemelo al esquema v7).

### Clasificación en 3 niveles
- **Nivel 1 — Solver de campo real (valor científico genuino):** difusión de O₂ (`simulador_oxigeno_biokidney.py`). Resuelve EDP reacción-difusión (Fick 3D) por SOR sobre el árbol real, con consumo Michaelis-Menten y realimentación de hipoxia. Física real.
- **Nivel 2 — Orden reducido correcto pero calibrado (estimaciones de diseño, no predicciones):** estrés mecánico dECM, Co-SWIFT/optimizador, filtración glomerular y reabsorción tubular (esta última en cuarentena).
- **Nivel 3 — Fenomenológico/decorativo (ilustrativo, sin validez predictiva):** iPSC simple (`simulador_ipsc_biokidney.py`, tasas "Vibe Coding" inventadas por comentario del autor), WSS simple (`simulador_wss_biokidney.py`). Caso aparte: diferenciación iPSC con EDO (mejor, EDOs tipo Hill referenciadas a Takasato/Freedman, pero constantes cinéticas elegidas, no ajustadas).

### Fallas concretas confirmadas
- **Reabsorción tubular:** NO parsea — `IndentationError` en la línea 208; salidas `n_ok=6` y `funcion_global=98.0` **hardcodeadas** (el propio comentario dice "simulado para brevedad del refactor"). Cifras no reproducibles.
- **Filtración glomerular:** dos `Kf_glom_nL` distintos entre archivos — **3.7** (no-G, "calibrado exacto") vs **4.1** (`_G`, "ajuste técnico"); el resultado de TFG depende de cuál se corra. La versión `_G` inyecta **+15 mmHg** de "vasodilatación aferente" dentro de una función de auto-ajuste; el "mapa de presión CCO" son posiciones `np.random`.
- **GFR_INPUT=82 mmHg hardcodeado** como constante de entrada, que contradice el ~62.5 de los sub-modelos.
- **Co-SWIFT:** la viabilidad "98%" es el **techo de un clamp [20, 98]** de `cell_viability_model`, no un óptimo calculado; funciones objetivo = curvas empíricas a trozos con puntos de quiebre elegidos.
- **Estrés dECM:** el mapa de Von Mises del andamio (Panel 7) es **un campo pintado con gaussianas**, no un cálculo espacial; el factor de poro `P_poro = P_ext x 0.85 x 0.03` son constantes a mano.
- **WSS simple:** usa presión absoluta en vez de caída de presión y define viscosidad sin usarla — **dimensionalmente mal aplicado**.
- **Teratoma y pureza:** salen de una EDO con constantes inventadas; no tienen calibración empírica y no pueden usarse como resultado de seguridad biológica.

### Archivos tocados
- Ninguno de código. Insumo para el marco de honestidad (Entrada 025).

### Estado
**Cerrado como auditoría.** Deuda derivada: (a) auditar el paquete `biokidney` (`cfg_physio.P_HIPOXIA`, `D_O2`, tasa de consumo de `CellularExpert`) del que cuelga el número del solver O₂; (b) reparar la reabsorción tubular (queda en cuarentena); (c) elegir y justificar un único Kf glomerular.

---

## ENTRADA 025 — 11 de julio de 2026 — Marco de honestidad epistémica (documento canónico, 3 niveles)

**Estado:** **CERRADO.** Creación del documento canónico `MARCO_honestidad_epistemica_BioKidney.md`, que formaliza la clasificación en 3 niveles de la Entrada 024 y fija la afirmación central defendible. Es el cuerpo único del que se derivan la sección de honestidad de Zenodo v3, la reescritura del estado del Dashboard Maestro y la espina del material de fundraising. Reemplaza toda afirmación de tipo "12/12 ÓPTIMO · VALIDACIÓN COMPLETA · 100%".

### Contenido
- **Afirmación central:** la contribución validada es la **fidelidad anatómica y geométrica** del gemelo digital multicapa (morfometría poblacional, ley de Murray k=3.0, cobertura geométrica) + **un único solver de campo real** (difusión de O₂). Los módulos funcionales son verificaciones de factibilidad de orden reducido, explícitamente calibradas — **no predicciones fisiológicas**. Se defiende geometría, no función.
- Clasificación por niveles (§2), correcciones transversales obligatorias (§3: teratoma/pureza fuera como resultados; calibración declarada, no oculta; Kf; reabsorción en cuarentena) y qué NO cambia (§5: Capa 5 separable).

### Archivos tocados
- Nuevo: `00_bitacora/MARCO_honestidad_epistemica_BioKidney.md`.

### Estado
**Cerrado.** Documento canónico de referencia para todo reencuadre posterior (preprint v3, Zenodo, deck).

---

## ENTRADA 026 — 11 de julio de 2026 — Preprint v3: reencuadre de honestidad sobre v2 (20 cambios + fila 21 del caption Fig 1)

**Estado:** **CERRADO.** Redacción de la versión v3 del preprint sobre la base de v2, aplicando el reencuadre de honestidad del marco canónico (Entrada 025). Se conserva intacta la ciencia vascular real (Murray k=3.0, tabla de presiones, morfometría, calibración Poiseuille de 2 pasos, solver O₂); se reposiciona todo lo calibrado/ilustrativo. **No se inventó ningún dato nuevo.**

### Cambios (20 filas + fila 21)
- **Título adoptado:** "A Multi-Layer Geometric Digital Twin of the Human Kidney…" (se retira el framing predictivo "Six-Module Integration Predicts Physiological Renal Output…").
- Abstract, Introducción y Conclusión reencuadrados por niveles (real-field solver / calibrated reduced-order / illustrative); se retiran como logros "full iPSC purity", "98% viability", "98.1% reabsorption".
- 115.2 mL/min reetiquetado como *feasibility check under a calibrated scale, not a geometric prediction*; 82 mL/min = estimación standalone con presión sintética. Se **fija Kf = 3.7** (se nota la variante 4.1 no usada).
- Reabsorción degradada a orden reducido y **retirada de resultados** (2.19 L/día, 98.1%, 6/6): módulo en cuarentena, no ejecuta end-to-end.
- Solver O₂ mantenido como campo real, con las 4 salvedades del marco (umbral = anoxia; longitud de difusión escala-órgano; "0% hipoxia" condicionado por resolución de grilla; dependencia de `cfg_physio.P_HIPOXIA` pendiente de auditar).
- Teratoma/pureza reformulados como objetivos de diseño de literatura (Takasato/Freedman), nunca como resultado de seguridad.
- §4.7 Limitations (párrafo nuevo): el multisimulador consume un grafo de segmentos (CSV/JSON), no una malla; el gemelo geométrico multicapa (Capas 0-4) es el resultado validado; el cruce hilar es límite de escala, no bloqueo del in-silico.
- **Fila 21:** Fig. 1 caption — frase puente de procedencia (CCO v8 -> fuente de Capas 0-4).

### Archivos tocados
- Nuevos: `00_bitacora/preprint_biokidney_2026_EN_v3.md`, `00_bitacora/preprint_biokidney_2026_EN_v3.pdf`, `00_bitacora/CHANGELOG_v2_a_v3.md`.
- El v2 (`preprint_biokidney_2026.md`) se conserva como base.

### Estado
**Cerrado.** Preprint v3 redactado y compilado; changelog v2->v3 documentado con las 21 filas y su razón (marco §).

---

## ENTRADA 027 — 11 de julio de 2026 — Zenodo v3 PUBLICADO (DOI de versión)

**Estado:** **CERRADO.** Publicación en Zenodo de la versión v3 (preprint con reencuadre de honestidad, Entrada 026), bajo el DOI de concepto ya existente del proyecto.

### DOIs
- **DOI de versión (v3):** 10.5281/zenodo.21314073
- **DOI de concepto (todas las versiones):** 10.5281/zenodo.19508076
- Antecedente: la primera publicación en Zenodo (11-abr-2026, ver entrada de esa fecha) usó DOI 10.5281/zenodo.19508077.

### Archivos tocados
- Ninguno en el repo (acción externa de publicación). Insumo: preprint v3 (Entrada 026) y sección de honestidad derivada del marco (Entrada 025).

### Estado
**Cerrado.** Registro público de v3 disponible con DOI citable.

## ENTRADA 028 — 13 de julio de 2026 — Auditoría del solver de O₂: diverge en producción (regresión sin línea base), contaminación `col`, blast radius acotado

**Estado:** **AUDITORÍA CERRADA (solo lectura, ningún fix aplicado).** Rastreo forense del único "solver de campo real" del framework (difusión de O₂). Hallazgo mayor: el solver **no converge** — diverge a NaN en producción. El efecto en cadena está acotado: ningún otro módulo consume el campo. Consecuencia epistémica sobre el documento canónico (Entrada 025): el framework queda con **cero módulos de Nivel 1**.

### 1. El solver de O₂ diverge (evidencia dura)
- **Causa:** `resolver_fick_3d` actualiza todo el interior en una sola asignación vectorizada leyendo el array viejo → es **Jacobi**, no Gauss-Seidel. Con **sobre-relajación ω=1.6** (`relajacion = 1.6`), Jacobi-sobre-relajado (JOR) es **inestable para ω>1**: el campo crece exponencialmente (40 → 9e99 → 4e200 → overflow) y colapsa a NaN hacia it≈1500–2000.
- **Producción, árbol por defecto (v7, 1448 seg), grilla 60×60×40:** **98.05% de NaN** en el campo P (141 197 / 144 000); **99.98% del tejido** (141 197 / 141 221); solo **24 vóxeles finitos**. Medido sobre la fuente exacta, sin editar el original.
- **Original intacto:** md5 de `senyal_hipoxia_para_cco.csv` sin cambios (`02221d8946b3cbeefbff4e4c0dc8092c`) antes y después.

### 2. Fallo silencioso (el reporte se autodelata)
`analizar_oxigenacion` usa `nanmean`/`nanmin`, que **ignoran los NaN**. Con 24 vóxeles finitos, el pipeline imprime `🎯 OBJETIVO ALCANZADO: 100% del tejido oxigenado`, `Hipóxicos: 0 (0.00%)`, `PO₂ mínima 20.000 / media 28.323` — pero la línea inmediatamente anterior dice `Oxigenados: 24 (0.02%)`. El 99.98% restante no es ni hipóxico ni oxigenado: es NaN. El reporte contradice su propio banner.

### 3. Regresión sin línea base recuperable (git)
- `git blame`: `relajacion = 1.6` existe desde el **commit raíz `^37f958c` (2026-03-26 19:59)**. No hay **ningún** commit con ω≤1 ni Gauss-Seidel. En git, el módulo **nunca convergió**.
- El único artefacto **finito** —`senyal_hipoxia_para_cco.csv`, **520 vóxeles hipóxicos**, 2026-03-24 01:50— es **anterior por ~2.7 días** al primer commit del módulo. La versión que convergía existió **pre-commit** y **nunca se versionó**: regresión sin línea base recuperable. (Nota: incluso esa versión reportaba 520 hipóxicos, no "cero".)

### 4. Contaminación `col` (sistema colector urinario en el árbol vascular)
- En v7/v8, `sistema=col` es el **sistema colector urinario** (pelvis/cálices; radios **850–3250 µm** vs máx venoso 542 µm; 72 segmentos), no vena colectora.
- El loader lo admite al término fuente vía el `return` por defecto de `pres()` → Dirichlet a **30 mmHg** (punto medio art-ven). Al rasterizar pinta **~1.4% del grid = 2.6× más volumen que arterias+venas juntas** (col 1.42% vs art+ven 0.54%, v7).
- **A/B estabilizado (ω=1.0, v8, 8000 iter):** quitar `col` **sube la hipoxia de 14.48% a 16.35%** (+2 952 vóxeles) y baja la media 1.26 mmHg → `col` **enmascaraba hipoxia real** como fuente espuria de O₂. |A−B|: media ~1.4, máx ~18–24 mmHg; 47k vóxeles >1 mmHg, 13k >5 mmHg. **Caveat:** ESTAB no convergió en 8000 iter; A y B reciben idéntico tratamiento (el delta es comparable), pero los absolutos no son definitivos.

### 5. Blast radius: contenido
- `resolver_fick_3d` tiene **un solo consumidor: él mismo**. Nada más lo importa.
- `senyal_hipoxia_para_cco.csv` se **escribe, no se lee** por ningún `.py`. El bucle "hipoxia→CCO" no existe: **el generador CCO (v2–v8) nunca la usó para crecer el árbol**.
- El campo P **no se persiste** (no hay `.npz`).
- `CellularExpert.oxygen_consumption_rate` lo consumen el solver roto y el **backend web** (`simulation_service.py:352`), que tiene su **propio** modelo O₂ 2D (Krogh + suavizado, acotado, **no diverge**; comparte `P_HIPOXIA=1.0` y el consumo).
- **NO dependen del O₂:** CCO v2–v8, filtración glomerular, viabilidad SWIFT/MOPSO (sus "hipoxia" son heurística de tiempo de impresión), iPSC, reabsorción.
- **Geometría limpia:** capas 0–4 (dominio, nefronas, arterial, venoso, colector, calicial) tienen **cero** consumo del campo O₂ / señal de hipoxia. Única "influencia" del O₂ sobre la geometría: la constante de diseño **150 µm** (cobertura por distancia KDTree en capa2/capa3), nunca el campo.

### Archivos tocados
- **Ninguno de producción.** Solo lectura + `01_simuladores/_ab_test.py` (arnés A/B; importa el original sin modificarlo) y un runner en scratchpad. Documentación: esta entrada + actualización de `MARCO_honestidad_epistemica_BioKidney.md` (§2, §5, §6).

### Nota — código committeado que nunca compiló (patrón)
Al limpiar el working tree (13-jul) se confirmó que **`analizador_proyecto_biokidney.py` no compila ni en HEAD**: la línea 9 (`def obtener_estructura_texto`) está indentada como método pero **no hay `class` que la contenga** — nunca fue Python válido (independiente de una corrupción de sesión posterior que se descartó). Es el **tercer módulo con código committeado que nunca corrió**, junto a: (i) **reabsorción tubular** (`simulador_reabsorcion_tubular.py`, IndentationError línea 208 — marco §3.4); (ii) el **solver de O₂** (corre pero **diverge** a NaN, §1 de esta entrada). Se deja **sin reparar** (no es ciencia, no citado por el paper, no bloquea ningún commit); registrado para dimensionar el código muerto committeado (ver barrido `py_compile` del 13-jul).

### Estado
**Cerrado (auditoría).** Deuda derivada, sin aplicar hasta decisión: (a) reparar el solver (Gauss-Seidel red-black o ω≤1) + reemplazar `nanmean`/`nanmin` por detección de NaN que falle ruidosamente; (b) corregir las cifras de O₂ del preprint v3 (§2.2/§3.2, líneas 26/106/108/210); (c) excluir `col` del término fuente.

## ENTRADA 029 — 16 de julio de 2026 — Cierre de la sesión de corrección: preprint v4 (O₂ retirado), reproducibilidad versionada, sincronización repo↔v4

**Estado:** **CERRADO.** Sesión de corrección que implementa las decisiones derivadas de la auditoría (Entrada 028). Se redactó el preprint v4 retirando el módulo de O₂ (no se arregló), se versionó la contribución central del paper que no estaba en git, y se sincronizó todo el repo con el nuevo registro de honestidad. Sin reescritura de historial; sin publicación en Zenodo (paso irreversible, aparte).

### 1. Preprint v4 redactado (`preprint_biokidney_2026_EN_v4.md`, commit `7870574`)
- **Módulo de O₂ retirado, no arreglado.** §2.2 (Krogh) y §3.4 (Oxygen Diffusion) **eliminadas**; añadida **§4.8 "Oxygen Transport — Withdrawn from This Version"** (divergencia Jacobi ω=1.6, contaminación `col`, límite de escala 0.6 mm vs 150 µm); añadida **§3.7 Reproducibility**.
- **Título** sin la cláusula "a Real-Field Oxygen Solver"; encuadre epistémico → el framework **no contiene módulos de nivel solver-real**.
- **GFR corregido:** 115.2→**115.4** bilateral, 57.6→**57.7** por riñón (aritmética 3.7×15.6; **target-vs-logrado aclarado** en §2.1.3: 55.5 = presión objetivo 58.0, 57.7 = lograda 58.6). Propagado a los 9 sitios de 115.2.
- **Referencias renumeradas a [1]–[12]:** se eliminaron 3 huérfanas reales (Krogh, Michaelis, Herschel) y se **restauró Grebenyuk** tras detectar que estaba citada en la **cita agrupada `[3,4]`** que el grep simple `\[N\]` no capturaba. Verificación group-aware: cero huérfanas.
- **Aclaración segmento-vs-nodo** en §3.3: 1.902 segmentos = 1.905 nodos − 3 (un nodo raíz por sistema art/ven/col).
- Verificación: grep de cifras de O₂ = **0** en el v4.

### 2. Reproducibilidad versionada (commit `141466b`)
- **Las Capas 0–4 (`08_gemelo_digital/`) y `renal_data_v1.json` NO estaban en git** — la contribución central del paper era irreproducible por clonación. Ahora versionadas.
- **Determinismo confirmado:** el generador v7 con semilla fija produce un CSV **idéntico bit-a-bit** (md5 `5999d564…`, **1.448 segmentos**: 692 art / 684 ven / 72 col).
- **Excluidos** (`.gitignore`, para no publicar geometría transitoria): `08_gemelo_digital/_backup_pre_seno/` (geometría pre-corrección del defecto cortical 152/1300) y **13 `.npz` experimentales** de la raíz (variantes `v1/v2/v3`, sweeps `wrep*`, auditorías, backups).

### 3. Sincronización repo↔v4 (commit `9571d79`)
- Alineados con el registro v4: `supplementary_material_v8.md` (115.4/57.7, O₂ fuera, six→five módulos), `README.md` (**retirado "Pipeline 100% Completado"** y "todas las fases de validación" → contribución = gemelo geométrico + feasibility checks + cuarentenas), `dashboard_maestro_app.py` (82→115.4 feasibility, reabsorción→cuarentena, viabilidad/WSS etiquetadas clamp/operating-point), `MAESTRO_contexto_canonico_BioKidney.md` (cinco módulos, 115,4, O₂/reabsorción en cuarentena, viabilidad/GFR/WSS recalificadas; guard de cifras prohibidas intacto).
- **Manuscritos pre-honestidad archivados** (commit `4988ca0`): `preprint_biokidney_2026.md` (v2 EN) y `..._ES.md` → `99_archivo/SUPERSEDED_*` (solo movidos, contenido intacto).

### 4. Auditoría de salud del repo
- **85 de 87 módulos versionados compilan** (`py_compile`).
- **Tres módulos en cuarentena confirmados:** (i) O₂ (`simulador_oxigeno_biokidney.py`, diverge en runtime); (ii) reabsorción tubular (`simulador_reabsorcion_tubular.py`, IndentationError L208); (iii) `analizador_proyecto_biokidney.py` (IndentationError L9, roto en HEAD — nunca compiló).
- **Deuda anotada:** el MVP web no arranca en un clon fresco por `init_db()` apuntando a `web_app/database/` sin `mkdir` (falta bootstrap, no falta código); **tres modelos de O₂ coexistiendo** — Fick 3D roto, Krogh 2D del backend web (acotado, vivo en la app), Krogh 30×30 del texto (retirado).

### 5. Preservado sin tocar (md5 verificado a lo largo de la sesión)
BITACORA previa (entradas ≤028), `MARCO_honestidad_epistemica_BioKidney.md`, `CHANGELOG_v2_a_v3.md`, `preprint_biokidney_2026_EN_v3.md` (DOI 10.5281/zenodo.21314073), `implementation_protocol_v1.md` (DOI 10.5281/zenodo.19508077). **No se reescribió historial.** (`INDICE_BioKidneyAI.md` y los `CONTEXTO_*.md` quedaron untracked, fuera de los commits.)

### 6. Pendiente
- **Publicar el v4 en Zenodo** — NO hecho (paso irreversible, decisión aparte); implica sincronizar también `implementation_protocol_v1.md` (tiene DOI propio, aún con cifras viejas).
- **Retomar el hilo del adaptador** (Entrada 023): poblar `es_terminal` real, **excluir `col`** del término fuente de O₂, y **repuntar los simuladores al gemelo** (`arbol_vascular_cco_v8_gemelo.csv`).

### Estado
**Cerrado.** Preprint v4 + reproducibilidad + sincronización commiteados (`7870574`, `141466b`/`03964c8`/`0153294`/`6565b32`, `9571d79`, `4988ca0`). El repo es clonable y su puerta de entrada (README/MAESTRO) dice lo mismo que el paper. Próximo paso de ciencia: cerrar el adaptador del gemelo.

---

## ENTRADA 030 — 25 de julio de 2026 — Corrección de reproducibilidad: el commit `141466b` versionó un stub, no el árbol

**Estado:** **CORREGIDO.** Durante la auditoría de radios del Apéndice A del preprint v4 se descubrió que la "contribución central versionada" en Entrada 029 (§2) era, para el vaso vascular, un placeholder — no el dato real. Corregido y verificado desde clon limpio. Solo staging/commit; ningún contenido de archivo modificado en esta corrección de datos.

### 1. El hallazgo
- El commit `141466b` ("versionar `renal_data_v1.json` — contribución central") había commiteado por error un **stub sintético de 3 segmentos (473 bytes)**, no el **export real de 1.902 segmentos** (861 KB). El diff de ese commit sobre el JSON era `1 +` (una sola línea); el stub es un JSON de una línea, sin campo `sistema`, con `nivel` 0/1 y radios 2500/1200 µm de juguete.
- **Causa:** `git add` del archivo equivocado — el disco tenía el stub en esa ruta en el momento de `141466b`; el export real de 1.902 segmentos **sobreescribió** el archivo en el working tree **después**, sin recommitear. Por eso aparecía como modificado (`M`, +40.220 líneas, sin stagear) hasta hoy.
- **NO fue `.gitignore`:** `git check-ignore -v 02_vascular_cco/renal_data_v1.json` → no ignorado. Las reglas del `.gitignore` solo tapan variantes transitorias (`capa3a_arterial_v*.npz`, `capa3b_venoso_v1_*`, `capa*_auditoria_*`, `_backup_pre_seno/`).
- **Los 10 `.npz` de Capas 0–4 del mismo commit SÍ eran reales:** verificado tamaño `git show HEAD:<path>` == working tree para los 10 (capa0 3.8 MB, capa2 2.3 MB, etc.). El único stub era el JSON.

### 2. La consecuencia (por qué pasó desapercibido toda la sesión)
- Todas las verificaciones de "reproducible por clonación" previas se hicieron contra el **working tree**, no contra un clon. Un `git clone` limpio habría entregado el **stub de 3 segmentos**, y la Tabla A1 del paper **no** habría sido reproducible desde el repo público.
- El clone-test de `141466b` (Entrada 029) verificó **presencia** del archivo en el clon, no su **contenido** — y el stub pasó esa prueba.

### 3. La corrección (commit `859f13e`)
- `git add 02_vascular_cco/renal_data_v1.json` + commit único: `fix(datos): versionar renal_data_v1.json real (1902 segmentos) — 141466b habia commiteado un stub de 3 por error` (`1 file changed, 40220 insertions(+), 1 deletion(-)`).
- **Verificado desde clon limpio** (`git clone .` → `/tmp/clone_test2`, no el working tree): **1.902 segmentos** en el JSON del clon; verificador (`verify_a1.py` con `ROOT` reapuntado al clon) → **Tabla A1 11/11 filas MATCH** (R_min/R_max, tol 0.1 µm) y **§3.1 8/8 cifras** (`452.4, 163.0, 265.3, 120.1, 222.1, 47.7, 137.5, 556.8`).
- Contexto: esta corrección cierra la coherencia con la edición de L108 del v4 (commit `1582831`), que afirma que el JSON "holds the complete export … versioned in the git repository" — ahora **verdadera** respecto de lo commiteado.

### 4. Regla de método derivada (extiende la regla de verificación-por-archivo de esta sesión)
- Verificar **"reproducible por clonación"** exige un `git clone` real a un **directorio limpio** y correr el verificador **ahí**, no en el working tree.
- **"El archivo está presente en el clon" ≠ "el archivo es correcto".** El clone-test de `141466b` verificó presencia, no contenido, y el stub pasó.
- **Todo artefacto de datos versionado se valida por tamaño/contenido desde el clon**, no por presencia. (Para binarios: `git show HEAD:<path> | wc -c` vs working tree; para datos parseables: contar registros desde el clon.)

### Estado
**Corregido y verificado desde clon.** El JSON real (1.902 seg) está en `859f13e`. Pendiente sin ejecutar (decisión aparte): push de `feature/cco-v8-fractal` → `origin` y merge `feature` → `main` (el default público sigue en `b550317`). No se reescribió historial; ningún contenido de archivo científico fue modificado en esta corrección.

---

## ENTRADA 031 — 25 de julio de 2026 — Publicación del preprint v4 en Zenodo + cierre de la sesión de corrección

**Estado:** **CERRADO — PUBLICADO.** El preprint v4 (honesto, con O₂ retirado y las capas geométricas versionadas y reproducibles) está publicado en Zenodo y todo el trabajo de la sesión está en `origin/main`. Cierre de la línea de corrección abierta tras la auditoría (Entradas 028–030).

### 1. Preprint v4 PUBLICADO en Zenodo
- **Version DOI:** `10.5281/zenodo.21576269` (v4, 25 jul 2026). **Concept DOI:** `10.5281/zenodo.19508076` (resuelve siempre a la última versión — **usar por defecto para citar**).
- Publicado como *new version* del registro existente; **cadena de versiones preservada:** v1 (`19508077`) → v2 (`20390315`) → v3 (`21314073`) → v4 (`21576269`). Historial intacto, **sin reescritura**.
- **Título final:** *"A Multi-Layer Geometric Digital Twin of the Human Kidney: Anatomical Fidelity and Calibrated Reduced-Order Feasibility Checks"* — sin la cláusula "Real-Field Oxygen Solver" del v3, coherente con la retirada del módulo O₂. **Licencia:** CC-BY 4.0. **Tipo:** Preprint.
- Nota de versión pública registrada: retirada del módulo O₂, corrección de §3.7 y Apéndice A, versionado de las capas geométricas.
- **Trampa evitada en el último paso:** Zenodo **precargó el título y el abstract deshonestos del v1/v2**; se corrigieron **manualmente** antes de *Publish*, copiando título y abstract de la portada del PDF v4 verificado.

### 2. Correcciones finales del v4 (esta sesión — commiteadas y pusheadas a `origin/main`)
- **§3.7 Reproducibility:** md5/conteo corregidos del árbol **v7** (`5999d564…`, 1.448) al **v8 real** (`4f881a05…`, 1.902 seg — 904 art / 926 ven / 72 col). **Determinismo del v8 VERIFICADO:** regenerando `generador_cco_v8.py` (semilla = 42) en `/tmp` aislado, el md5 del `renal_data_v1.json` producido coincide **bit-a-bit** con el versionado. Commit `7b51009`.
- **§4.7:** retirado el overclaim residual *"whole-organ functional prediction"* → *"geometric fidelity for pre-fabrication analysis"*.
- **Fig 1 caption:** aclarado que el banner *"R_min=45.3 µm"* es el **mínimo global art+ven**, mientras el texto reporta **47.7 µm** (mínimo arterial). Ambos verificados.
- **§2.1.3:** verificado que el **55.5** (target de calibración, P<sub>gc</sub>=58.0) vs **~58** (logrado, 58.6) ya estaba correctamente distinguido; **sin cambios**.
- **PDF regenerado** con las 4 correcciones y **verificado por extracción de texto** (`pdftotext`: `4f881a05` presente; `5999d564`, `whole-organ functional`, `1,448 seg—692` ausentes). Commit del PDF ya en `origin/main`.

### 3. Estado de git al cierre
- **`main` = `origin/main`, sincronizada.** Contiene todo el trabajo de la sesión (v4 corregido, PDF corregido, `renal_data_v1.json` real de 1.902 seg, README/supplementary/dashboard/MAESTRO en registro de honestidad, bitácora).
- **`feature/cco-v8-fractal` quedó un paso atrás** (en `63304f3`, sin las 4 correcciones finales) — **divergencia aceptada deliberadamente**; `main` es la rama canónica de aquí en adelante.
- **Merge `feature`→`main` (`e72c7aa`)** ejecutado y pusheado.
- **Autenticación GitHub migrada a SSH** (clave ed25519 registrada). `recovery-codes.txt` movido fuera de la carpeta del repo (a `~/`, **nunca estuvo trackeado**).

### 4. DEUDA PENDIENTE (próximas sesiones — no urgente, no bloquea nada)
- **Higiene del repo:** working tree sucio con `.pyc` trackeados (agregar a `.gitignore` + `git rm --cached`); carpeta anidada `Bio-Kidney-AI-2026/` con material viejo (decidir archivar); tres `simulador_..._ipsc (2)/(3)/(4).py` duplicados; `requirements.txt` con `stripe==8.10.0` sin commitear; PDFs generados trackeados. Stashes `wip-pre-merge-25jul` sin resolver. `main`→`feature` opcional si se quiere reconciliar ramas.
- **Ruta científica** (el hilo que abrió esta sesión): el adaptador `.npz`→CSV; poblar `es_terminal` real; decidir exclusión de `col` (urinario) del CSV vascular; repuntar simuladores al gemelo. El solver de O₂ queda **EN CUARENTENA, no se arregla** (límite de escala: terminales 0.6 mm vs difusión 150 µm; sin capilares por diseño).
- **Módulos en cuarentena confirmados:** O₂ (diverge), reabsorción tubular (`IndentationError` L208), `analizador_proyecto_biokidney.py` (`IndentationError` L9). Deuda menor: bootstrap `init_db()` del MVP sin `mkdir`; tres modelos de O₂ coexistiendo.

### 5. HALLAZGO DE MÉTODO (el aprendizaje central de la sesión — más importante que cualquier bug técnico)
Un **ejecutor puede reportar verificaciones que no ejecutó**, y una **plataforma (Zenodo) puede precargar metadata obsoleta**. **Regla instituida:** ningún reporte de estado cuenta como verificación; **todo cambio numérico o publicación se valida por** (i) **artefacto reproducible** (script versionado, md5, extracción de texto del PDF), (ii) **ejecutado/inspeccionado directamente por el operador**, (iii) **cruzado contra la fuente real** — nunca contra el resumen del agente ni contra campos precargados. **La fuente de verdad es el estado real, verificado en el momento, sin excepción** — ni para el código, ni para Claude Code, ni para Zenodo.

### Estado
**Cerrado y publicado.** Preprint v4 en Zenodo (Concept DOI `10.5281/zenodo.19508076`; version DOI `10.5281/zenodo.21576269`). `main` = `origin/main` con todo el trabajo de honestidad. Próximo paso de ciencia: cerrar el adaptador del gemelo y repuntar los simuladores.

---

## ENTRADA 032 — 3 de agosto de 2026 (registrada el 6 de agosto) — Capa 0: corrección del campo de profundidad cortical, de cuerda radial a distancia perpendicular (método B)

**Estado:** **CORREGIDO Y VERIFICADO.** El campo `depth_cortical_mm` de Capa 0 medía una magnitud que no era un espesor. Corregido en el commit `fbb6614`. Impacto recomputado desde `capa0_dominio.npz` y validado desde terminal. **Pasivo abierto: Capas 1-4 no regeneradas.**

### 1. Qué se corrigió
- **Antes:** la profundidad cortical era la **cuerda radial** desde el centroide del elipsoide hasta la cápsula, a lo largo del rayo centroide→punto. Esa función se conserva deprecada como `capsule_distance_radial` (`capa0_dominio.py:128`), con cuerpo `return np.clip(rsurf_main - r_main, 0.0, None)`.
- **Ahora:** la profundidad es la **distancia euclídea al punto más cercano de la cápsula externa** (elipsoide 55/30/18). Implementada en `capsule_distance` (`capa0_dominio.py:221`), cuerpo en `:244-245`:
  ```python
  x = _nearest_point_ellipsoid(coords, SEMIEJES)
  return np.linalg.norm(np.atleast_2d(coords) - x, axis=-1)
  ```
- El pie de la perpendicular lo resuelve `_nearest_point_ellipsoid` (`capa0_dominio.py:140`) por multiplicador de Lagrange con bisección vectorizada, reparametrizada en `mu = lam + min(a_i^2)`. Esa reparametrización **no es cosmética**: en los puntos que caen sobre un plano coordenado la raíz degenera justo en el borde del bracket, y el docstring (`:156-163`) registra las dos formas documentadas de equivocarse en silencio (10.125 mm y 12.000 mm en `(0,-18,0)`, contra los 11.906 correctos). El bloque `VERIFICACION` de `main()` incluye 7 casos de guarda numérica para esto (`capa0_dominio.py:513-553`; la lista de casos en `:522-530`).
- **Cadena vigente:** `main()` (`:388`) → `compute_depth()` (`:248`) → `capsule_distance()` (`:221`) → `_nearest_point_ellipsoid()` (`:140`); persistencia en `:437-438` (`depth`, `depth_cortical_mm`).

### 2. Por qué: se comparaban dos magnitudes distintas
- El umbral córtico-medular es `GROSOR_CORTICAL_MM = 6.6` (`capa0_dominio.py:61`), y la fuente lo define como **espesor cortical PERPENDICULAR a la cápsula** (*cortical width*).
- La cuerda radial **no es un espesor**: sobreestima el espesor real en todo punto cuyo rayo no sea normal a la superficie, y tanto más cuanto más oblicuo es. En un elipsoide achatado (55/30/18) la oblicuidad es grande en amplias zonas.
- Comparar 6.6 mm (perpendicular) contra una distancia radial mezclaba unidades conceptuales distintas. El docstring de `capsule_distance` (`:225-233`) deja la justificación anclada en el propio código.
- Propiedad que se sigue del cambio y que sirve de control: la perpendicular es **≤** la radial en todo punto. Recomputado: `perp <= radial en todo punto: True`. Por eso el reetiquetado sólo puede ir en un sentido.

### 3. Impacto verificado (recomputado desde `capa0_dominio.npz`, N = 200 000)
| Magnitud | Córtex | % |
|---|---|---|
| Radial (deprecada) | 92 624 | **46.31 %** |
| Perpendicular (vigente, la del `.npz`) | 117 703 | **58.85 %** |

- **Reetiquetados: 25 079 puntos = 12.54 % del parénquima.**
- **Dirección: 25 079 médula→córtex, 0 córtex→médula.** Estrictamente unidireccional, como exige `perp <= radial`.
- **Reparto: 46.31 / 53.69 → 58.85 / 41.15.**
- Coincide con lo declarado en el mensaje de `fbb6614` (`reetiqueta 25079 pts medula->cortex; reparto 46.3/53.7 -> 58.85/41.15`).
- Comando de reproducción (no escribe en disco; el guard `if __name__ == "__main__":` de `:610` impide que el import dispare `main()`):
  ```bash
  .venv/bin/python -c "import numpy as np, capa0_dominio as C; d=np.load('capa0_dominio.npz',allow_pickle=True); co=d['coords'].astype(np.float64); dm=d['depth_cortical_mm'].astype(np.float64); rad=C.capsule_distance_radial(co); G=float(C.GROSOR_CORTICAL_MM); cr,cp=rad<G,dm<G; print('cortex radial',int(cr.sum()),round(100*cr.mean(),2),'% | cortex perp',int(cp.sum()),round(100*cp.mean(),2),'%'); print('reetiquetados',int((cr!=cp).sum()),round(100*(cr!=cp).mean(),2),'% | med->cortex',int((~cr&cp).sum()),'cortex->med',int((cr&~cp).sum()))"
  ```

### 4. El contrafactual: la medida de lo que se evitó
La corrección del método B se apoya sobre una corrección anterior — que la profundidad dependa **sólo** de la cápsula, no del mínimo con la pared del seno. Cuantificado hoy sobre el `.npz` vigente:

- **31 311 puntos (15.656 %)** están más cerca de la pared del seno que de la cápsula (`dist_seno < dist_main`).
- De ellos, sólo **4 097** son córtex bajo el criterio actual — corticales legítimos de la esquina cápsula/seno.
- Si el campo persistido usara `np.minimum(dist_capsula, dist_seno)`, el córtex sería **99 004 = 49.50 %** y **29 307 etiquetas diferirían** de las actuales.

Esos 29 307 puntos son la magnitud del *córtex peri-sinusal falso* que la separación cápsula/seno eliminó: puntos interiores marcados corticales por vecindad a una superficie **interna**, bajo la cual hay médula y grasa sinusal, no corteza.

Nota de trazabilidad: `nearest_surface_distance` (`capa0_dominio.py:111`) sigue en el archivo, con su `np.minimum` intacto en `:125`, pero **sin un solo llamador en todo el repositorio**. Es código muerto declarado — su docstring (`:112-119`) anuncia el retiro y la razón. Su permanencia es deliberada, no un residuo.

### 5. Corrección de atribución (separada del cambio de código)
Esta corrección es **independiente** del método B y no altera ninguna cifra del gemelo. Afecta a de dónde procede el número 6.6 mm.

- **El NÚMERO 6.6 mm procede de Glodny B, Unterholzner V, Taferner B, Hofmann KJ, Rehder P, Strasak A, Petersen J. BMC Urology 2009;9:19** (DOI 10.1186/1471-2490-9-19): *cortical width* CW = 6.6 ± 1.9 mm (dcho), 6.6 ± 2.0 mm (izdo), MDCT 64-cortes, **n = 2068 riñones / 1040 adultos asintomáticos**, ICC 0.96.
- **El MÉTODO procede de Beland MD, Walle NL, Machan JT, Cronan JJ. AJR 2010;195(2):W146-149**, que define el protocolo explícito de medición **perpendicular a la cápsula** (plano sagital, sobre pirámide medular) — orientación que Glodny no detalla.
- **La atribución previa del NÚMERO a Beland era errónea.** Beland es n=25 pacientes con ERC (edad media 73) y su media es 5.9 mm: contiene el 6.6 en su rango 3.2–11.0 mm, lo que es validación *post hoc*, **no procedencia**.
- **Alcance de la propagación (contado con grep, no de memoria):** `grep -rn "Beland"` sobre el repo (excluyendo `.git`, `.venv`, `__pycache__`) devuelve **10 ocurrencias en 2 archivos**. Desglose:
  - `capa0_dominio.py:58` — **1 ocurrencia, ya correcta**: `# ANCLA: Glodny et al. 2009, MDCT n=2068 (NO Beland; ver correccion de atribucion).` El código **nunca** atribuyó el número a Beland.
  - `09_paper_vascular/auditoria_correspondencia_anatomica.md` — **9 ocurrencias**. De ellas, **un único sitio preserva el texto erróneo literalmente**: `:142-143`, la línea del listado de citas colocadas el 2026-07-06 (`GROSOR_CORTICAL_MM (6.6 mm, grosor cortical) → Beland MD, Walle NL, Machan JT, Cronan JJ. AJR 2010;195(2):W146-149`). Se conserva **a propósito, sin modificar**, con el bloque de corrección inmediatamente debajo (`:144-152`). Las 8 ocurrencias restantes (`:9`, `:12`, `:15`, `:47`, `:114`, `:116`, `:150`, `:151`) son **texto de la corrección**, no propagación del error.
- **Conclusión del alcance:** la atribución errónea **no está viva en ningún punto del repo**. Está preservada en un sitio, marcada como errónea en ese mismo sitio. No queda nada que corregir por este concepto.
- Salvedad de método: el grep cubre archivos de texto. **No cubre el contenido de los PDF** del repo (`00_bitacora/*.pdf`, preprints). Si la atribución a Beland viajó a un PDF publicado, ese conteo no lo detecta y queda fuera de esta entrada.

### 6. Limitación declarada: cápsula fantasma
El elipsoide principal es una superficie **cerrada**, pero el tramo excavado por el seno no existe como cápsula real. Un punto mide contra cápsula fantasma si el pie de su perpendicular cae dentro del elipsoide del seno.

- **596 puntos = 0.298 % del parénquima** miden contra cápsula fantasma.
- De ellos, **209 están etiquetados córtex**.
- Todos a **≤ 4.1 mm** de la pared del seno (mediana **0.75 mm**): es una banda delgada peri-sinusal, no un volumen difuso.
- La zona es la misma que Glodny excluye explícitamente en su protocolo (evitar columnas de Bertin y seno renal), de modo que la limitación cae donde la fuente tampoco mide.
- Cuantificado en vivo por el bloque `VERIFICACION` de `main()` (`capa0_dominio.py:570-588`).

### 7. La trampa evitada — mecanismo, no anécdota
Durante la auditoría del **2026-08-06** se emitió un misdiagnóstico: se reportó como "inconsistencia interna abierta" que `capa0_dominio.py:53-55` declarase que la pared del seno ya no genera córtex mientras `:123-125` seguía devolviendo `np.minimum(dist_main, dist_seno)`. **Se creyó reabierto un bug que estaba cerrado desde julio.**

**Mecanismo del error, paso a paso:**
1. Los documentos `09_paper_vascular/diagnostico_causa_raiz_seno.md` y `09_paper_vascular/diagnostico_seguridad_campo_profundidad.md` (ambos del **2026-07-05**, previos a `fbb6614`) afirman que `compute_depth` **usa** `nearest_surface_distance`, y lo clasifican `[DEPENDE-AMBAS]` (`diagnostico_seguridad_campo_profundidad.md:22`).
2. Esos documentos citan **rangos de línea caducos** (`:98-111`, `:114-128`, `:268`, `:271-272`). En el archivo actual esos rangos apuntan a **otras funciones**. La cita numérica aparenta precisión y por eso no se cuestiona.
3. La búsqueda por nombre (`grep "nearest_surface_distance"`) devuelve **7 aciertos**, de los cuales **sólo 1 es la definición** y **ninguno es una llamada**. El resto son prosa de los `.md` obsoletos y la copia de respaldo. Un conteo de aciertos sin distinguir *definición / llamada / prosa* sugiere una función viva.
4. Se leyó el cuerpo de la función (`:121-125`) pero **no su docstring** (`:112-119`), que dice literalmente `[CONSERVADA - ya NO define la profundidad cortical]` y explica que la etiqueta usa ahora `capsule_distance()`. La respuesta estaba a cuatro líneas de donde se miró.

**Reglas derivadas:**
- **Un `grep` por nombre no es un grafo de llamadas.** Separar siempre definición, invocación y mención en prosa antes de afirmar que algo está vivo.
- **La documentación forense caduca y su precisión aparente la hace peligrosa.** Un `.md` con `archivo:línea` envejece peor que uno sin ellos, porque invita a confiar sin recomprobar. → Acción tomada: ambos documentos llevan ahora un bloque de aviso al inicio (2026-08-06) que los declara superados por `fbb6614` y remite a esta entrada. Se conservan íntegros por su valor forense.
- **Leer el docstring antes de leer el cuerpo.** El estado declarado de una función vive en su docstring; el cuerpo sólo dice qué calcularía si alguien la llamara.
- **Antes de reportar un bug reabierto, comprobar contra el artefacto**, no contra el código fuente. Bastaba comparar `depth_cortical_mm` del `.npz` con ambas funciones candidatas: `max |depth_npz − capsule_distance| = 1.9e-06 mm` (redondeo `float32`) frente a `max |depth_npz − nearest_surface_distance| = 14.03 mm`. Eso zanja la pregunta en un comando.
- Esta entrada extiende la regla instituida en la Entrada 031 (§5): **ningún reporte cuenta como verificación**. Aquí el reporte defectuoso lo produjo el agente auditor, y lo que lo desmintió fue el artefacto.

### 8. Commit de referencia
- **`fbb6614`** — *"Capa 0: corrige depth cortical de radial a perpendicular (Glodny 2009 CW, 6.6mm)"*, 2026-08-03 23:58:54 -0500.
- Alcance: `capa0_dominio.py` (+219 líneas), `capa0_dominio.npz` (Bin 3875256 → 3859371 bytes), `09_paper_vascular/auditoria_correspondencia_anatomica.md` (+33). Total: 3 archivos, 242 inserciones, 10 supresiones.
- `coords` bit-idéntico respecto del estado anterior: **sólo cambian etiquetas, no geometría**. SVD de la nube regenerada: rango 3.
- Backup del estado previo en `_backup_pre_depth_2026-07-30/`.

### Estado
**Corregido en Capa 0; verificado desde terminal.** El `.npz` vigente contiene la profundidad perpendicular y el reparto 58.85/41.15.

**PASIVO ABIERTO — Capas 1-4 NO regeneradas.** Declarado en el propio mensaje de `fbb6614` (`Capas 1-4 NO regeneradas, quedan desincronizadas hasta la cascada`) y vigente desde el 2026-08-03. Las capas 1-4 se construyeron sobre el etiquetado córtex/médula **anterior**, en el que 25 079 puntos hoy corticales figuraban como médula. Toda capa que siembre desde el pool `region == "cortex"` — Capa 1 siembra los 1300 glomérulos así — está desincronizada respecto de Capa 0. **Hasta que corra la cascada, Capa 0 y Capas 1-4 no describen el mismo riñón.** Ninguna cifra que cruce ambos niveles debe publicarse en ese estado.

---

## ENTRADA 033 — 6 de agosto de 2026 — El seno renal: supuesto geométrico sin ancla morfométrica, y la distinción de magnitudes del volumen

**Estado:** **DECLARADO, NO ANCLADO.** El seno renal del gemelo es un supuesto de diseño geométricamente consistente, sin cita morfométrica que lo sostenga. Esta entrada fija su estado epistémico y separa dos magnitudes que se venían confundiendo.

### 1. Clasificación epistémica: SUPUESTO GEOMÉTRICAMENTE CONSISTENTE, no anclado
- Los dos parámetros que definen el seno son `CENTRO_SENO = [0,−34,0]` (`capa0_dominio.py:48`) y `SEMIEJES_SENO = [22,16,11]` (`capa0_dominio.py:49`).
- **Ninguno tiene ancla en literatura.** Ambos figuran como pendientes en la auditoría de correspondencia anatómica: `09_paper_vascular/auditoria_correspondencia_anatomica.md:128-129`, puestos 1 y 2 de los **10 parámetros clase [A] sin cita** (`1. CENTRO_SENO = [0,−34,0] mm (capa0_dominio.py:48) — posición del seno renal` / `2. SEMIEJES_SENO = [22,16,11] mm (capa0_dominio.py:49) — dimensiones del seno renal`).
- Lo que sí está anclado es el elipsoide **principal**: semiejes 55/30/18 → Emamian SA, Nielsen MB, Pedersen JF, Ytte L. AJR 1993;160(1):83-86. El seno se talló **por dentro** de esa forma anclada, calibrado por criterio geométrico ("forma de frijol sin partir el órgano en lóbulos"), no por morfometría del seno.
- El propio código lo declara en el comentario de `capa0_dominio.py:45-47`: `Centrado mas hacia afuera en -Y; al restarse del dominio crea la invaginacion del frijol. Ajustado para ocupar la zona media de la cara medial sin partir el organo en dos.`
- **Consecuencia para publicación:** el seno es un supuesto de diseño declarado. Es legítimo como geometría idealizada paramétrica; **no** es reconstrucción anatómica ni está validado contra población. Cualquier figura o texto que lo presente debe decirlo.

### 2. Distinción de magnitudes — el punto central de esta entrada
Se venían usando indistintamente dos números que **no son la misma cosa**:

| Magnitud | Valor | Qué es |
|---|---|---|
| **Elipsoide del seno completo** | **16.22 mL** (16 218.996 mm³) | Volumen analítico del elipsoide de exclusión: `4/3·π·22·16·11` |
| **Intersección con el parénquima** | **3.9599 mL** (3 959.877 mm³) | Volumen **excluido** del elipsoide principal, lo que `capa0_dominio.py:430` calcula |

- **3.96 mL NO es el volumen del seno renal.** Es el volumen **excluido por intersección**: la parte del elipsoide del seno que cae dentro del elipsoide principal y por tanto se resta del parénquima. Los otros ~12.26 mL del elipsoide del seno quedan **fuera** del elipsoide principal y no recortan nada.
- Expresión exacta que lo produce, `capa0_dominio.py:426-430`:
  ```python
  # 6) Volumenes (Monte Carlo)
  vol_main = 4.0 / 3.0 * np.pi * np.prod(SEMIEJES)
  accept_frac = n_accept / n_tried           # parenquima / elipsoide principal
  vol_parenquima = vol_main * accept_frac
  vol_seno_excluido = vol_main - vol_parenquima
  ```
- Valores en el `.npz` vigente: `vol_elipsoide_mm3 = 124407.0690821558`, `vol_parenquima_mm3 = 120447.19207327077`, diferencia **3959.8770088850288 mm³**.

**Etiquetado incorrecto localizado en el repo (grep, no memoria):**
- `00_bitacora/BITACORA.md:690` — `- Volumen parénquima: 120.447 mm³ (96.8% del elipsoide; seno excluye 3.960 mm³).` Esta línea es **correcta**: dice *"seno excluye"*, no *"volumen del seno"*. Se deja como está.
- `00_bitacora/_archivo/2026-06-12_bitacora_gemelo_digital.md:68` — texto idéntico, igualmente correcto.
- **`grep -rn "3\.96\|3,96"` sobre el repo no devuelve ningún sitio que etiquete 3.96 como "volumen del seno renal".** El único otro acierto, `09_paper_vascular/diagnostico_holgura_pelvis.md:20` (`pared_art 3.96`), es una **holgura de pared arterial en mm**, magnitud distinta y sin relación.
- **La cadena `3.96 mL` no existe en ningún archivo del repo.** El número nunca se ha escrito en mililitros en disco. Si aparece en material externo (deck, preprint, correspondencia) etiquetado como volumen del seno, ese material está mal y esta entrada es la referencia para corregirlo.

### 2-bis. COTA INFERIOR DE VOLUMEN VIOLADA (añadido 2026-08-11, tras el anclaje de Caglar 2014)

Cuando se redactó esta entrada no existía referencia de volumen del seno y §7 lo declaraba pendiente. **Ya existe una cota, y el gemelo la viola.**

**El argumento es de contención, no de ajuste.** El seno renal contiene grasa **más** el sistema colector (pelvis, cálices), la arteria y vena renales y sus ramas, linfáticos y tejido conectivo. Por tanto, para un mismo riñón:

> **V_seno > V_grasa_sinusal**, siempre y por definición.

Caglar 2014 mide **sólo la grasa** del seno por estereología sobre TC (n = 240, 21–80 años). Sus medias por grupo:

| grupo | grasa sinusal (cm³) |
|---|---|
| hombres, riñón izquierdo | **5.70 ± 2.87** |
| hombres, riñón derecho | 4.15 ± 2.39 |
| mujeres, riñón izquierdo | 3.51 ± 2.67 |
| mujeres, riñón derecho | 2.49 ± 2.16 |

**Aritmética de la violación.** El volumen **excavado** del gemelo (la intersección seno ∩ parénquima, que es el hueco real que el modelo abre) es **3.9599 mL**:

| comparación | cociente seno/grasa | ¿cumple V_seno > V_grasa? |
|---|---|---|
| 3.9599 vs 5.70 (H izq) | **0.700** | **NO — violada por 1.740 mL** |
| 3.9599 vs 4.15 (H der) | **0.962** | **NO — violada por 0.190 mL** |
| 3.9599 vs 3.51 (M izq) | 1.128 | sí, con margen de sólo 0.450 mL |
| 3.9599 vs 2.49 (M der) | 1.590 | sí |

**Lectura:** el hueco que el gemelo excava es *más pequeño que la grasa sola* de un riñón izquierdo masculino medio, y apenas mayor que la grasa sola de un riñón derecho masculino medio. Como en ese hueco debe caber además todo el sistema colector y el pedículo vascular — que es exactamente lo que Capas 3b/4/5a intentan meter ahí —, **la cota inferior está violada en los dos grupos masculinos**.

- El gemelo **no declara lateralidad ni sexo**: **[PENDIENTE DE ANCLA]** cuál de los cuatro grupos le corresponde. Por eso la cota se reporta contra los cuatro.
- La magnitud comparable es el **volumen excavado (3.9599 mL)**, no el elipsoide completo (16.22 mL). El elipsoide no es un volumen anatómico: el 76 % de él cae fuera del parénquima y no excava nada.
- Reproducible:
  ```bash
  .venv/bin/python -c "import numpy as np; d=np.load('capa0_dominio.npz'); v=(float(d['vol_elipsoide_mm3'])-float(d['vol_parenquima_mm3']))/1000; print('excavado %.4f mL'%v); [print('  vs %-6s %.2f -> %.3f %s'%(k,g,v/g,'VIOLA' if v<g else 'ok')) for k,g in [('H izq',5.70),('H der',4.15),('M izq',3.51),('M der',2.49)]]"
  ```

**Esto invierte el diagnóstico de §6.** Ver la revisión de esa sección más abajo.

### 3. Fragilidad: el número no se persiste
- `vol_seno_excluido` se calcula en `capa0_dominio.py:430` y **nunca se guarda**. El bloque `np.savez_compressed` (`:433-457`) persiste `vol_parenquima_mm3` (`:451`) y `vol_elipsoide_mm3` (`:452`), pero no el volumen excluido. Las claves de volumen del `.npz` son exactamente `['vol_parenquima_mm3', 'vol_elipsoide_mm3']`.
- El valor sólo existe como salida de consola, `capa0_dominio.py:470-471`:
  ```python
  print(f"  Seno excluido       : {vol_seno_excluido:11.1f} mm^3  "
        f"({100.0 * (1 - accept_frac):5.2f} % del elipsoide)")
  ```
- **Por tanto, toda cita del número hoy es transcripción de consola o resta manual**, no lectura de un campo almacenado. Es el mismo patrón de fragilidad que la Entrada 030 (el stub de `renal_data_v1.json`): un dato que se afirma pero cuyo soporte versionado no existe.
- Mitigación mientras no se corrija — re-derivarlo del `.npz`, sin escribir en disco:
  ```bash
  .venv/bin/python -c "import numpy as np; d=np.load('capa0_dominio.npz'); print((float(d['vol_elipsoide_mm3'])-float(d['vol_parenquima_mm3']))/1000,'mL')"
  ```
- **Corrección propuesta y DIFERIDA:** añadir en `capa0_dominio.py` la línea `vol_seno_excluido_mm3=np.float64(vol_seno_excluido),` inmediatamente después de `:452`. **No aplicada.** Aplicarla obliga a re-ejecutar `capa0_dominio.py`, que sobrescribe `capa0_dominio.npz`; con las Capas 1-4 ya desincronizadas (Entrada 032, §Estado), esa regeneración debe decidirse **dentro de la cascada**, no de forma aislada.

### 4. Naturaleza estocástica del número
- `vol_seno_excluido` **no es analítico**: es un estimador **Monte Carlo** por muestreo con rechazo. `sample_parenchyma` (`capa0_dominio.py:272`) genera puntos uniformes en el elipsoide principal y rechaza los que caen dentro del seno (`:292`); la fracción de aceptación estima el volumen relativo.
- Parámetros que lo determinan: **`SEED = 2026`** (`capa0_dominio.py:71`) y **`N_POINTS = 200_000`** (`capa0_dominio.py:70`).
- Es reproducible bit a bit con esa semilla, pero **cambiar la semilla cambia el último dígito**. Los 4 decimales de 3.9599 mL no son precisión física: son precisión de un muestreo concreto. Reportar el número sin declarar semilla y N es reportar de más.
- Reproducción desde cero, re-ejecutando el muestreo sin escribir en disco:
  ```bash
  .venv/bin/python -c "import numpy as np, capa0_dominio as C; rng=np.random.default_rng(int(C.SEED)); c,t,a=C.sample_parenchyma(int(C.N_POINTS),rng); vm=4/3*np.pi*np.prod(C.SEMIEJES); print('seno excluido mm3 =',vm-vm*a/t,' -> mL =',(vm-vm*a/t)/1000)"
  ```

### 5. Orientación y rotación: simplificaciones de diseño declaradas
- El seno es un **elipsoide alineado con los ejes**, sin rotación: `SEMIEJES_SENO = [22,16,11]` se aplica directamente sobre X/Y/Z. En el riñón real el eje del seno no es paralelo al eje del órgano.
- Está **centrado en X = 0** (`CENTRO_SENO = [0,−34,0]`): perfectamente simétrico respecto del ecuador del órgano, sin la asimetría supero-inferior anatómica.
- Está **centrado en Z = 0**: sin inclinación anterior-posterior.
- El hilio se define como el punto medial puro `HILIO = [0,−30,0]` (`capa0_dominio.py:42`), no como una región.
- **Las tres son simplificaciones de diseño, no hallazgos.** Se declaran aquí para que ninguna figura las presente como anatomía. El seno del gemelo reproduce la *topología* (invaginación medial, forma de frijol, papilas drenando hacia la cavidad) sin reproducir la *orientación*.
- Consistencia interna sí verificada: 0 puntos de parénquima dentro del seno (rechazo por construcción, `:292`) y las 10 pirámides medulares apuntan sus papilas a la pared del seno (`build_pyramids`, `capa0_dominio.py:305`).

### 6. ¿Está el seno subdimensionado y aprieta a las capas superiores?
**No hay evidencia en el repo que lo sostenga, y sí evidencia que lo contradice.** El 3.96 mL **no** es argumento para esta afirmación en ningún sentido: es un volumen de intersección, no una holgura.

La única evidencia localizable sobre holgura en la cavidad es `09_paper_vascular/diagnostico_holgura_pelvis.md` (Capa **5a**, 2026-07-10; barrido de 525 posiciones de la pelvis). Sus números propios:

- `:16` — `centros probados: 525  ->  VIABLES (pared_ven>=0.4 Y pared_art>=0.4 Y seno<=1): 17`. **El criterio de contención en el seno (`seno<=1`) lo cumplen 17 centros**: el seno **sí** contiene la pelvis. Los valores de contención de los candidatos van de `seno 0.52` a `seno 1.00` (`:17-28`).
- `:49` — `NINGUN centro D2-viable cierra D3-D5. Fallos entre 15 candidatos: D3(infundibulos) 15/15, D4(ureter) 15/15, D5(corredor) 0/15`.
- `:56` — `BINDING = el PLEXO VENOSO peri-hilar (32 segs, 27 posteriores) ocupa el corredor Z~0 obligado.`
- `:57` — `-> el fix NO es posicion NI tamaño de la pelvis: es SEGREGAR/mover el plexo venoso (revisar 3b).`
- `:64` — veredicto estructurado: `"status": "NO_ES_LA_PELVIS__BINDING_PLEXO_VENOSO_Z0"`, con `"mejor_pared_inf_ureter_um": -964.3`.

**Lectura:** el documento identifica el cuello de botella en el **plexo venoso peri-hilar**, no en las dimensiones del seno, y descarta explícitamente el tamaño como la variable a tocar. Confirma cuantitativamente el blocker de las Entradas 019-021 (árbol venoso ramificado en el hilio).

**Salvedades que impiden cerrar la pregunta con este documento:**
- Es de **Capa 5a (pelvis)**, no de Capa 4. **Para Capa 4 específicamente: [PENDIENTE DE ANCLA — sin respaldo en repo].**
- Es del **2026-07-10**, anterior a `fbb6614`. Se apoya en geometría de capas que hoy están desincronizadas (Entrada 032, §Estado). Debe recomputarse tras la cascada antes de tratarse como vigente.

**Redacción autorizada:** *"el seno contiene la pelvis en 17 de 525 posiciones barridas; el binding declarado es el plexo venoso peri-hilar, no el tamaño del seno (diagnostico_holgura_pelvis.md:56-57), pendiente de recomputar tras la cascada."* **Redacción NO autorizada:** cualquier variante de *"el seno está subdimensionado"*, con o sin el 3.96 como apoyo.

> **REVISIÓN 2026-08-11 — la prohibición de arriba queda LEVANTADA.** El párrafo anterior se conserva sin modificar como registro de lo que se sabía el 2026-08-06. Su razonamiento era correcto **con la evidencia disponible entonces**: el 3.96 mL no era argumento porque no había con qué compararlo, y `diagnostico_holgura_pelvis.md` sólo hablaba de la *pelvis*, no del *volumen*.
>
> **Caglar 2014 aporta la cota que faltaba** (§2-bis). La afirmación pasa a estar sostenida, pero **por una vía distinta de la que se había intentado**: no por holgura de Capa 5a, sino por violación de la cota inferior de volumen.
>
> **REDACCIÓN AHORA AUTORIZADA** (usar literalmente, o una paráfrasis que conserve las tres precisiones marcadas):
>
> > *"El volumen **excavado** por el seno del gemelo — **3.9599 mL**, la intersección del elipsoide de exclusión con el parénquima, **no** el elipsoide completo de 16.22 mL — es **menor que el volumen de grasa sinusal sola** reportado por Caglar et al. (Folia Morphol. 2014;73(3):302-308) para riñón izquierdo masculino (**5.70 ± 2.87 cm³**) y riñón derecho masculino (**4.15 ± 2.39 cm³**), medido por estereología sobre TC en n = 240 sujetos. Como el seno renal contiene, además de grasa, el sistema colector, el pedículo vascular, linfáticos y tejido conectivo, se cumple necesariamente V_seno > V_grasa. **El seno del gemelo está por tanto subdimensionado** frente a esa cota inferior en los dos grupos masculinos."*
>
> **Tres precisiones que la redacción debe conservar siempre:**
> 1. **El número es el excavado (3.9599 mL), nunca el elipsoide (16.22 mL).** El elipsoide no es un volumen anatómico.
> 2. **Caglar mide grasa, no seno.** La comparación es contra una **cota inferior**, no contra un valor de referencia del seno. Nunca escribir "el seno debería medir 5.70 mL".
> 3. **La cota se viola en los grupos masculinos** (izq y der); frente a los femeninos el gemelo queda por encima, con margen estrecho. Si no se declara lateralidad ni sexo, decir "en los dos grupos masculinos", no "siempre".
>
> **Sigue NO AUTORIZADO:** (a) presentar 16.22 mL como volumen del seno; (b) afirmar que el seno "aprieta a Capa 4" o a Capa 5a — eso sigue siendo `[PENDIENTE DE ANCLA]`, y `diagnostico_holgura_pelvis.md:56-57` sigue señalando el plexo venoso como binding; (c) derivar un tamaño de seno cruzando Caglar con Zhang 2023 (ver §7).

### 7. Qué haría falta para anclar el seno
- Una fuente morfométrica del **seno renal** (volumen, ejes o proporción respecto del parénquima) sobre población adulta sana, medida por MDCT o RM. Glodny 2009 mide *parenchymal width* (PW ≈ 15.5 mm, cápsula → seno) pero **no** dimensiona el seno como cavidad.
- Hasta entonces, `CENTRO_SENO` y `SEMIEJES_SENO` permanecen en la lista de 10 pendientes de `auditoria_correspondencia_anatomica.md:126-137`, puestos 1 y 2.
- **Sin esa fuente, el seno no debe aparecer en ninguna tabla de parámetros anclados del preprint.**

> **ACTUALIZACIÓN 2026-08-11 — dos anclas PARCIALES incorporadas.** Registro completo en `04_literatura/anclas_seno_renal.md`.
>
> **[ANCLA PARCIAL 1] Caglar V, Kucuk A, Aktas S, et al.** *"Volumetric evaluation of fat in the renal sinus in normal subjects using stereological method on computed tomography images and its relationship with body composition."* **Folia Morphologica. 2014;73(3):302-308.** DOI **10.5603/FM.2014.0016**. PMID **25242158**.
> n = 240 sujetos, 21–80 años, TC, método estereológico. **Grasa del seno renal:** hombres izq **5.70 ± 2.87** / der **4.15 ± 2.39** cm³; mujeres izq **3.51 ± 2.67** / der **2.49 ± 2.16** cm³.
> **Qué ancla:** una **cota inferior** de volumen del seno (§2-bis). **Qué NO ancla:** el volumen del seno, sus ejes ni su posición. Mide grasa, que es un subconjunto propio del contenido sinusal.
>
> **[ANCLA PARCIAL 2] Zhang QH, et al.** *Frontiers in Endocrinology.* **2023;14:1187781.** DOI **10.3389/fendo.2023.1187781**.
> n = 126, RM, mapeo de fracción grasa. **Fracción grasa del seno:** hombres **28.33 ± 6.73** (der) / **31.21 ± 6.29** (izq); mujeres **23.82 ± 7.74** / **27.92 ± 8.15**.
> **Qué ancla:** la **proporción** de grasa dentro del seno. **Qué NO ancla:** ningún volumen absoluto.
>
> **[NO ANCLADO — HIPÓTESIS CONVERGENTE, PROHIBIDO USAR COMO PARÁMETRO]** Cruzar la fracción de Zhang con el volumen de Caglar sugeriría un seno total de orden ~14–15 cm³, cercano al elipsoide de exclusión de 16.22 mL, lo que apuntaría a que **el tamaño del elipsoide es aproximadamente correcto y el error está en `CENTRO_SENO`** — coherente con el diagnóstico independiente de `08_gemelo_digital/_auditoria/2026-08-09` y `2026-08-10`, donde `CENTRO_SENO[1]` resulta ser la variable dominante y `SEMIEJES_SENO` no puede corregir el defecto de los ápices.
> **Por qué NO puede usarse como valor:** los volúmenes absolutos de ambos estudios difieren en casi un orden de magnitud por metodologías distintas (estereología sobre TC vs mapeo de fracción grasa por RM, poblaciones y criterios de segmentación distintos). **Un número derivado cruzando dos estudios no comparables no es un anclaje.** Se registra como convergencia cualitativa entre dos vías independientes, jamás como parámetro.
>
> **SIGUE PENDIENTE:** las dimensiones del área ecogénica central de **Emamian 1993** no aparecen en ninguna fuente secundaria accesible. **[PENDIENTE DE ANCLA]** para `SEMIEJES_SENO` y `CENTRO_SENO` como valores. Lo que Caglar y Zhang permiten hoy es **acotar y falsar**, no fijar.
>
> **Aviso de integridad (2026-08-11):** circuló un documento externo con dimensiones *reconstruidas* de la CEA de Emamian presentadas como datos del paper (L = 4.5, W = 1.3, T = 1.2 cm; V = 3.68 cm³). **Esas cifras no están en la fuente y no deben usarse.** Búsqueda en el repo: **no aparecen** (los únicos aciertos de `3.68` son coordenadas del árbol vascular en `02_vascular_cco/renal_data_v1.json`).

### Estado
**Declarado y acotado, no anclado.** Las dos magnitudes quedan separadas: elipsoide del seno completo **16.22 mL**, volumen excluido por intersección **3.9599 mL** (Monte Carlo, SEED 2026, N 200 000). El repo no contiene hoy ningún etiquetado incorrecto del 3.96; la cadena `3.96 mL` no existe en disco.

**Pendientes que esta entrada abre:**
1. Persistir `vol_seno_excluido_mm3` en el `.npz` — línea propuesta en §3, **diferida a la cascada**.
2. Anclar `CENTRO_SENO` y `SEMIEJES_SENO`, o mantenerlos declarados como supuesto en toda publicación.
3. Recomputar `diagnostico_holgura_pelvis.md` tras la cascada, y sólo entonces pronunciarse sobre la holgura de Capa 4.
4. **Bloqueante heredado de la Entrada 032:** Capas 1-4 sin regenerar.

> **ADENDA 2026-08-11 — el estado epistémico cambia de NO ANCLADO a ACOTADO.**
> El encabezado de esta entrada (`DECLARADO, NO ANCLADO`) se conserva sin modificar: era exacto el 2026-08-06. Hoy la situación es distinta y esta adenda la fija.
>
> - **Lo que cambia:** existe una **cota inferior anclada** (Caglar 2014) y el gemelo la **viola** en los dos grupos masculinos (§2-bis). La redacción *"el seno está subdimensionado"* pasa de prohibida a autorizada, con las tres precisiones de §6.
> - **Lo que NO cambia:** `SEMIEJES_SENO` y `CENTRO_SENO` **siguen sin ancla de valor** y siguen en los puestos 1 y 2 de los 10 pendientes. Caglar y Zhang permiten **falsar**, no **fijar**. El seno **sigue sin poder aparecer** en una tabla de parámetros anclados del preprint.
> - **Estado resultante:** **SUPUESTO DECLARADO, ACOTADO POR ABAJO Y FALSADO.** Un supuesto que ahora se sabe incorrecto en una dirección conocida — que es más que "no anclado" y menos que "anclado".
>
> **Pendiente nuevo que esta adenda abre:**
> 5. La corrección de `CENTRO_SENO` deja de ser opcional: hoy corrige **cuatro** defectos independientes — ápices corticales, PW frente a Glodny, espesor medular, y ahora la cota de volumen de Caglar. Intervalo admisible bajo las tres condiciones ancladas: **`CENTRO_SENO[1] ∈ [−31.37, −27.80]`** (0 ápices corticales · PW dentro de ±1 s.d. de Glodny · volumen excavado > 5.70 mL). Diagnóstico y barridos en `08_gemelo_digital/_auditoria/2026-08-09_seno_apices_particion.md` y `_auditoria/2026-08-10_centro_seno_particion.md`; fuentes en `04_literatura/anclas_seno_renal.md`.

> **REVISIÓN 2026-08-12 — el pendiente 5, tal como está redactado arriba, contiene DOS ERRORES METODOLÓGICOS.** El texto anterior se conserva sin modificar como registro de lo que se afirmó el 2026-08-11. Auditoría completa en `08_gemelo_digital/_auditoria/2026-08-12_reservas_pw_region2d.md`.
>
> **Error 1 — «cuatro defectos independientes» es incorrecto.**
> Las tres condiciones geométricas son **estrictamente monótonas crecientes** en `CENTRO_SENO[1]`. Barrido con aserto de monotonía en 13 pasos (`sy = 16`):
>
> | `cy` | `depth` mín de ápice | PW_max | V excavado (mL) |
> |---|---|---|---|
> | −34.0 | 5.295 | 12.00 | 3.994 |
> | −32.0 | 6.757 | 14.00 | 5.275 |
> | −30.0 | 8.058 | 16.00 | 6.635 |
> | −28.0 | 9.190 | 18.00 | 8.128 |
>
> Y la monotonía es **estructural, no numérica**: las tres miden, con métricas distintas, **cuánto se adentra el seno en el parénquima**. Por tanto **no podían discrepar en dirección**, y su acuerdo **no es corroboración mutua**. La formulación correcta es: *cuatro síntomas del mismo defecto que coinciden en la **DIRECCIÓN**, no cuatro anclas independientes que convergen en un **VALOR**.* Lo que sí aportan por separado son **umbrales distintos**, y con ellos los dos bordes de la región.
>
> **Error 2 — el intervalo `[−31.37, −27.80]` está mal planteado y queda RETIRADO.**
> Se cumple la identidad `PW_max = |−B_SEMI − (cy + sy)|` (contrastada contra malla en 5 configuraciones, diferencia ≤ 1.8 × 10⁻⁴ mm). Luego **PW depende únicamente de la SUMA `cy + sy`**, no de `cy` por separado. El intervalo publicado es un **corte en `sy = 16`**, y `sy` está tan sin anclar como `cy` (puestos 1 y 2 de los 10 pendientes). Además `sx` y `sz` tampoco están anclados: **el espacio real es 4D en `(cy, sy, sx, sz)`**.
>
> **Formulación correcta de lo que las anclas restringen hoy:**
> - **Glodny (PW ±1 s.d.):** restringe la **suma**, `cy + sy ∈ [−17.4, −11.8]`. Es una banda diagonal del plano `(cy, sy)`, no un intervalo en `cy`.
> - **Caglar (volumen excavado > 5.70 mL):** una región curva en `(cy, sy, sx, sz)`; cota inferior.
> - **Región admisible bajo ambas**, con `sx = 22` y `sz = 11` fijos *(rebanada 2D de un espacio 4D)*: banda diagonal de anchura 1.8–4.4 mm en `cy` según `sy` — p. ej. `sy = 12` → `cy ∈ [−28.50, −24.10]`; `sy = 16` → `cy ∈ [−31.30, −27.90]`; `sy = 22` → `cy ∈ [−35.70, −33.90]`. Tabla completa en `_auditoria/2026-08-12`, Reserva 2c.
>
> **Error 3 (derivado) — los 4 ápices corticales NO son vinculantes.**
> Con `sy = 16`, los umbrales son: ápices `cy ≥ −32.50`; PW ≥ 12.6 `cy ≥ −33.40`; **volumen > 5.70 `cy ≥ −31.30`**; **PW ≤ 18.2 `cy ≤ −27.80`**. Es decir: **Caglar fija el borde inferior, Glodny el superior, y la condición de ápices es REDUNDANTE en ambos bordes** — no fija nada, se satisface automáticamente en cuanto se cumple Caglar. **Siguen siendo un defecto real** (4 papilas en territorio córtex, `depth` 5.295 y 6.428 mm frente al umbral 6.6), y su diagnóstico sigue en pie; **pero no son un pilar del intervalo** y presentarlos como uno de los cuatro criterios de calibración era inexacto.
>
> **Premisa incorrecta sobre el protocolo de Glodny.**
> Se había asumido que Glodny 2009 mide *"una medida por riñón en localización estandarizada"*. **La fuente no lo dice.** Consultado el texto completo (BMC Urology 2009;9:19, PMC2813848): especifica `width of the parenchyma (PW) and the cortex (CW) in the arterial phase` y el pie de la Figura 1 lo sitúa en **corte axial**; el `measurements were performed twice in a random sample of 50 data sets` es control de reproducibilidad, **no** el protocolo de rutina. Quedan **[PENDIENTE DE ANCLA]**: (i) el **nivel anatómico** al que se mide el PW (polo superior / tercio medio / hilio / polo inferior), (ii) el **número de medidas por riñón**, (iii) si es **localización estandarizada o promedio**.
>
> **Consecuencia:** no se sabe qué estadístico del gemelo es comparable con el de Glodny. La **dirección** del hallazgo es robusta —máximo 12.000, p95 11.528, mediana 6.855 y media 6.590 caen **todos** por debajo de 12.6 mm—, pero la **magnitud** del déficit no lo es (3.4 mm con el máximo, 8.5 mm con la mediana), y la magnitud es justo lo que fijaría el borde superior de la región.
>
> **Qué queda FIRME tras esta revisión:** la **cota inferior de volumen violada** (§2-bis) — no depende de PW, ni del estadístico elegido, ni de `sy`: compara dos volúmenes bajo el mismo criterio de contención. Y que **`CENTRO_SENO` es la variable dominante** (barrido de `_auditoria/2026-08-10` §1c, independiente de estas reservas).
>
> **Qué NO debe commitearse ni publicarse:** ninguna cifra de intervalo para `CENTRO_SENO[1]`, en ninguna forma, hasta (i) determinar el estadístico de PW de Glodny y (ii) reformular la región sobre `(cy, sy)` o sobre la suma `cy + sy`.

## ENTRADA — 2026-08-22 · Higiene documental y MAESTRO v2

Se detectó deriva entre los documentos canónicos y el estado real del proyecto
(preprint v4, Capas 0–4 cerradas). Acciones:

- Creado MAESTRO_contexto_canonico_BioKidney_v2.md como fuente única.
  Verificado contra preprint_biokidney_2026_EN_v4.md y BITÁCORA desde terminal.
- Fijada la TFG canónica en ~115 mL/min sin decimal (115,2 y 115,4 superadas;
  el v4 eliminó el decimal por tratarse de un chequeo aritmético manual).
- Registrado como sobreclaim el término "pipeline integrado predictivo": el
  módulo de filtración corre sobre el árbol CCO v7 legacy y no consume el
  campo de presión v8 (v4 §2.4, §3.3).
- Archivados con cabecera SUPERSEDED: implementation_protocol_v1.md,
  BITACORA_sesion_resometimiento_biorxiv.md, preprint_EN_v3.md,
  MAESTRO v1, informe de desarrollo del dashboard (mar 2026).
- INDICE_BioKidneyAI.md marcado como superado, pendiente de regeneración.
- Parcheadas las instrucciones de los cuatro agentes de marca.

Hallazgo relevante: el informe del dashboard describía el estado previo a la
corrección de dashboard_maestro_app.py (entrada 1458). El código fue alineado
al v4; su documentación no. Verificar en futuras correcciones que la
documentación de un artefacto se corrija junto con el artefacto.

## ENTRADA — 2026-09-21 · Cierre de higiene documental (frente A)

Continuación de la entrada del 22 ago 2026. Criterio: ningún archivo activo
contiene cifras retiradas (MAESTRO §0) ni encuadre predictivo; ante
ambigüedad, se archiva.

Commits:
- 51a3944: MAESTRO v2 y archivado del 22 ago (trabajo de esa fecha,
  commiteado hoy).
- 508e97c: archiva preprint EN v1 (PDF), 01_simuladores/Informes/,
  01_simuladores/screenshots/, 07_presentacion_final/ y dos imágenes de
  preprints superados (a 99_archivo/images/).
- 8dc4190: archiva dashboard_maestro_app.py, simulador_reabsorcion_tubular.py
  y run_pipeline.sh; inicio.sh L54 pasa a comentario.
- 4eb7029: archiva simulador_oxigeno_biokidney.py; _ab_test.py con ROOT
  relativo, apuntando al archivo.

Movidos a 99_archivo/ sin commitear: material que nunca estuvo en git
(preprints PDF superados, documentos de contexto, un manuscrito de abril,
un deck de presentación externa, diffs de auditoría). No se publica
material contaminado que nunca fue público.

Verificado desde terminal:
- Las 7 rutas del repo citadas por preprint_biokidney_2026_EN_v4.md siguen
  intactas.
- _ab_test.py ejecutado con .venv: con omega=1.6 diverge a NaN (it=1524);
  con omega=1.0 no converge en 8000 iteraciones. El md5 de
  senyal_hipoxia_para_cco.csv no cambia tras la ejecución.

Hallazgos (pendientes, fuera del alcance de esta sesión):
- Figura 3 del v4 (image_f02d63.png): la imagen no corresponde a su caption
  (sin paneles A/B, sin líneas de 90 y 60 mL/min, sin inset), y su título y
  la etiqueta "Calibrated" contradicen el encuadre de feasibility check. Se
  queda en su ruta por ser la versión depositada; se corrige en una v5.
- Figura 2 del v4: el mínimo de presión terminal coincide con el umbral de
  43 mmHg; revisar si es un piso impuesto por la calibración.
- Módulos iPSC y dECM: sus veredictos contradicen sus propios datos.
- web_app/ sigue calculando métricas de O2, de flujo urinario y el
  porcentaje de TFG nativa.
- Simuladores de filtración glomerular: comparaciones con el riñón nativo
  en código activo.

Lecciones de proceso:
- Un proceso automatizado ejecutó 15 movimientos de archivos (ctime
  11:14:08, en ~150 ms) antes de la ejecución manual, pese a reportes de
  "no ejecutado". Origen pendiente de confirmar. Regla: mientras se ejecuta
  en terminal, el agente no recibe instrucciones, y se comprueba que el
  índice esté vacío (git diff --cached) antes de cada commit.
- Un mensaje de commit afirmó una verificación antes de ejecutarla (la
  primera corrida falló por el entorno); se confirmó después. Regla:
  verificar, leer la salida, y solo entonces redactar.
