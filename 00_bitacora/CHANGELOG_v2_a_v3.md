# CHANGELOG — preprint Bio-Kidney AI 2026, v2 → v3

**Fecha:** 2026-07-11 · **Autor:** Carlos David Moreno Cáceres (VirtusSapiens)
**Base:** `preprint_biokidney_2026_EN_v2` (= `preprint_biokidney_2026.md`) → `preprint_biokidney_2026_EN_v3.md`
**Fuente de las decisiones:** `MARCO_honestidad_epistemica_BioKidney.md` (documento canónico).

**Regla rectora:** la contribución validada es la **fidelidad anatómica y geométrica** (gemelo digital
multicapa) + un **solver de campo real** (O₂). Los módulos funcionales son **verificaciones de
factibilidad de orden reducido, calibradas — no predicciones**. Se conserva intacta la ciencia
vascular real (Murray k=3.0, tabla de presiones, morfometría, calibración Poiseuille de 2 pasos,
solver O₂); se reposiciona todo lo calibrado/ilustrativo. No se inventó ningún dato nuevo.

---

## Cambios por sección

| # | Sección (v3) | Cambio | Razón (Marco) |
|---|---|---|---|
| 1 | **Título** | Se quitó "Six-Module Integration Predicts Physiological Renal Output…". Nuevo título centrado en "geometric digital twin + reduced-order feasibility checks". | Marco §1: se defiende geometría, no función; el framing predictivo es indefendible. |
| 2 | **Abstract** | Reformulado: contribución validada = fidelidad anatómica/geométrica + solver O₂ real. 115.2 mL/min reetiquetado como *feasibility check under calibrated scale, not a geometric prediction*. Se retiran del abstract, como logros, "full iPSC purity", "98% viability", "98.1% reabsorption"; se mencionan sólo como objetivos de literatura. Se conservan Murray k=3.0, Poiseuille 2-pasos, O₂. | Marco §1, §3.1, §3.3, §3.5. |
| 3 | **1. Introduction** | "complete in-silico validation of a functional bioprinted human kidney" → framing por niveles (real-field solver / calibrated reduced-order / illustrative). | Marco §3.5 (quitar sobreafirmación), §2 (tiering). |
| 4 | **2.3 iPSC (métodos)** | Reformulado a módulo *ilustrativo*; constantes cinéticas **elegidas, no ajustadas**; pureza/teratoma = **objetivos de diseño** (Takasato/Freedman), nunca resultado de seguridad. Se quitó "low teratoma risk" como hallazgo. | Marco §3.1 (teratoma/pureza fuera como resultados), §2 Nivel 3. |
| 5 | **2.4 Co-SWIFT (métodos)** | Se mantienen Herschel-Bulkley + MOPSO como reales; se declara que las funciones objetivo son **curvas empíricas a trozos** y que el 98% viabilidad es el **techo (clamp [20,98])**, no un óptimo; 5.6 dyn/cm² = operating point de entrada. | Marco §2 Nivel 2 (Co-SWIFT), §3.3. |
| 6 | **2.5 Filtración (métodos)** | Se declara: 115.2 = K_f(lit) × ΔP_Starling con presión calibrada (k=0.0132) → **calibración, no predicción emergente**. 82 mL/min = estimación de módulo standalone con presión **sintética** (constante de entrada). Se **fija K_f = 3.7**; se nota que existía una variante 4.1 (`_G`) **no usada**. | Marco §3.2 (calibración declarada), §3.3 (82 como entrada), §2 Nivel 2 (elección de K_f). |
| 7 | **2.6 Reabsorción (métodos)** | Degradado a orden reducido con fracciones **fijadas a valores de libro**; salidas **hardcodeadas** (comentario "simulado para brevedad"); módulo **no ejecuta end-to-end**; cifras (2.19 L/día, 98.1%, 6/6) **retiradas**; requiere reparación antes de reportar. | Marco §3.4 (cuarentena), §2 Nivel 2. |
| 8 | **Figura 4 (caption, panel F)** | "Functional validation: GFR=115.2" → "Reduced-order feasibility check under calibrated scale (… not a geometric GFR prediction)". | Marco §3.3, §3.5. |
| 9 | **3.3 GFR (resultados)** | Se conserva la descomposición de Starling (matemática real). Se reetiqueta "central functional outcome / functionally viable output" → **feasibility/consistency check**; 82 = constante de entrada, inconsistente con sub-modelos (~62.5). Se explica que cae en rango **por construcción** (calibración). | Marco §3.2, §3.3. |
| 10 | **3.4 O₂ (resultados)** | Se mantiene como **solver de campo real**. Se añaden las 4 salvedades del Marco: umbral = **anoxia** no hipoxia; longitud de difusión **escala-órgano**; "0% hipoxia" condicionado por **resolución de grilla**; dependencia crítica de `cfg_physio.P_HIPOXIA` **pendiente de auditar**. | Marco §2 Nivel 1 (salvedades O₂). |
| 11 | **3.5 iPSC (resultados)** | Se quitan "100% purity by day 21" y "teratoma low" como resultados; se dejan como objetivos de diseño de literatura. | Marco §3.1. |
| 12 | **3.6 Bioprinting (resultados)** | 98% viabilidad = **clamp/techo** de la función; 60 Pa / 5.6 dyn/cm² = **valores de entrada**, no predicciones convergidas. | Marco §2 Nivel 2, §3.3. |
| 13 | **3.7 Reabsorción (resultados)** | **Sin resultados reportados**: cifras 2.19 L/día, 1.52 mL/min, 98.1%, 1200 mOsm/kg, 6/6 **retiradas**; módulo en **cuarentena**. Sólo objetivo de literatura ~97–99%. | Marco §3.4. |
| 14 | **4.1 Principal Findings** | Se elimina "complete in-silico validation" y "convergence of all six modules toward physiologically plausible outputs"; el outcome principal se reencuadra como el **campo hemodinámico calibrado** (método), no un número de GFR. | Marco §3.5, §1. |
| 15 | **4.2 (cadena causal)** | "densidad cortical → nº de glomérulos perfundidos → superficie de filtración → GFR" → **correlato morfológico**, no mecanismo computado (el pipeline no calcula esa cadena). | Corrección señalada para v3; Marco §1. |
| 16 | **4.6 (iPSC + Co-SWIFT)** | iPSC "predicts full purity / teratoma low" → objetivos de literatura; Co-SWIFT 98%/5.6 = clamp/entrada. | Marco §3.1, §2 Nivel 2. |
| 17 | **4.7 Limitations** | **Párrafo nuevo:** el multisimulador consume un **grafo de segmentos (CSV/JSON)**, no una malla; la fabricación imprimible (voxelización, Capa 5) es fase posterior/separable; el **gemelo geométrico multicapa (Capas 0–4) es el resultado validado**; el cruce hilar es límite de escala, no bloqueo del in-silico. Párrafo de estado de módulos (reducido/cuarentena/ilustrativo). Se corrige el cierre ("predict whole-organ filtration performance" → check de factibilidad autoconsistente). | Marco §5 (qué no cambia / Capa 5 separable). |
| 18 | **5. Conclusion** | Se elimina "complete in-silico validation" y la lista de logros funcionales (full purity, 98% viab., 98.1% reabs.); resultados defendibles = geométricos y de campo; 115.2 = feasibility check; reabsorción en cuarentena. | Marco §3.5, §1. |
| 19 | **Tabla v7→v8 (4.5) + Fig 3 caption** | Fila GFR: "Within normal range" → "**Feasibility check, in-range by calibration**" con nota al pie; caption Fig 3 (B) añade "feasibility check under a calibrated scale — not a predicted GFR". | Consistencia con §3.3. |
| 20 | **Encabezado** | Nota de versión "v3 revision (epistemic tiering / honesty framing), July 2026". | Trazabilidad Zenodo. |
| 21 | **Fig. 1 caption** | Figura 1 caption: frase puente de procedencia (CCO v8 → fuente de Capas 0–4). | Consistencia con §4.7 (gemelo multicapa validado). |

## Consistencia de números de entrada (verificado)
Se marcaron explícitamente como **constante de entrada / literatura**, no como convergencia de pipeline:
`82 mL/min` (standalone, presión sintética), `5.6 dyn/cm²` (operating point), `1.52 mL/min` y `2.19 L/día`
(retirados de reabsorción), `98%` (clamp de la función de viabilidad). `115.2 mL/min` = feasibility check
bajo escala calibrada (K_f=3.7 × ΔP_Starling calibrada, k=0.0132).

## Lo que NO cambió (ciencia sólida conservada)
Murray k=3.0000 y 100% cumplimiento (915 bifurcaciones); tabla de presiones por nivel y Apéndice A;
morfometría vascular (calibres, radios); calibración Poiseuille de dos pasos (k=0.0132) y su justificación
metodológica (§4.4); solver de difusión de O₂ como único solver de campo real; figuras vasculares (1, 2) y
sus datos numéricos; referencias.

## Títulos alternativos propuestos (para tu elección final)
- **A (adoptado):** *A Multi-Layer Geometric Digital Twin of the Human Kidney: Anatomical Fidelity, a Real-Field Oxygen Solver, and Calibrated Reduced-Order Feasibility Checks.*
- **B:** *Bio-Kidney AI 2026: An Open, Reproducible Geometric Digital Twin of the Kidney — Murray-Compliant Vasculature, a Real Oxygen-Field Solver, and Honestly-Tiered Reduced-Order Functional Checks.*
