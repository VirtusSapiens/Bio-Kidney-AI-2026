# Marco de honestidad epistémica — Bio-Kidney AI 2026
### Clasificación por niveles de los módulos funcionales del multisimulador

**Autor:** Carlos David Moreno Cáceres (VirtusSapiens) · ORCID 0009-0005-3933-5072
**Propósito del documento:** cuerpo único y verdadero del que se derivan (i) la sección de honestidad de Zenodo v3, (ii) la reescritura del estado del Dashboard Maestro, y (iii) la espina del material de fundraising. Reemplaza toda afirmación de tipo "12/12 ÓPTIMO · VALIDACIÓN COMPLETA · 100%".

---

## 1. Afirmación central (la única defendible)

> La contribución validada de Bio-Kidney AI es la **fidelidad anatómica y geométrica** de un gemelo digital multicapa, reproducible en hardware modesto, con integración abierta: morfometría poblacional, ley de Murray (k=3.0), cobertura geométrica y un único solver de campo real (difusión de O₂). Los módulos funcionales son **verificaciones de factibilidad de orden reducido y exploración del espacio de diseño, explícitamente calibradas — no predicciones fisiológicas.**

Se defiende **geometría, no función**. Un solver de campo real más varias estimaciones de ingeniería honestamente etiquetadas es una contribución sólida y publicable. Seis "validaciones" verdes uniformes son un blanco fácil para el primer revisor que abra el módulo más débil.

> **Actualización 2026-07-13 (auditoría — BITACORA Entrada 028):** la cláusula "un único solver de campo real (difusión de O₂)" de la afirmación central **ya no se sostiene**. El solver de O₂ **diverge en producción** (98.05% de NaN en el campo P) y queda **EN CUARENTENA**. El framework tiene ahora **cero módulos de Nivel 1**. La afirmación defendible se reduce a **fidelidad geométrica** (Capas 0–4), que la auditoría confirmó **limpia** de cualquier dependencia del O₂. Ver §2 (Nivel 1 vacío), §5 y §6.

---

## 2. Clasificación por niveles

### Nivel 1 — Solver numérico de campo real

**VACÍO a partir del 2026-07-13.** El framework queda con **cero módulos de Nivel 1.**

El único candidato, la **difusión de oxígeno** (`simulador_oxigeno_biokidney.py`), fue **degradado y puesto EN CUARENTENA** tras la auditoría del 2026-07-13 (BITACORA Entrada 028). Evidencia:
- **Diverge numéricamente.** `resolver_fick_3d` es un update **Jacobi** vectorizado con **sobre-relajación ω=1.6** (`relajacion = 1.6`); JOR es inestable para ω>1. El campo crece exponencialmente hasta overflow → **98.05% de NaN** en el campo P (24 vóxeles finitos de 141 221 de tejido) en la corrida de producción con el árbol por defecto.
- **Fallo silencioso.** `analizar_oxigenacion` usa `nanmean`/`nanmin`, que enmascaran los NaN: el pipeline imprime "0% hipoxia / 100% del tejido oxigenado" mientras su propia línea reporta "Oxigenados: 24 (0.02%)".
- **Regresión sin línea base.** `git blame`: ω=1.6 desde el commit raíz `^37f958c` (2026-03-26). El único artefacto finito (`senyal_hipoxia_para_cco.csv`, 520 vóxeles hipóxicos, 2026-03-24) es anterior al primer commit → la versión que convergía nunca se versionó.

Por tanto **no es un solver operativo** y no puede sostener ninguna afirmación de oxigenación. Queda reclasificado como **NO OPERATIVO / EN CUARENTENA**, junto a reabsorción tubular (§3.4). Las salvedades de interpretación que antes lo acompañaban (anoxia ≠ hipoxia; longitud de difusión de escala-órgano; "0%" condicionado por resolución de grilla; config `cfg_physio.P_HIPOXIA` / `D_O2` / consumo de `CellularExpert` sin auditar) siguen vigentes pero son ahora **secundarias**: el módulo ni siquiera converge. Ver §6 para los hallazgos completos y la deuda derivada.

### Nivel 2 — Orden reducido correcto pero calibrado (estimaciones de diseño, no predicciones)

**Estrés mecánico dECM** (`simulador_estres_mecanico_dECM.py`).
Real: Kelvin-Voigt para creep, Von Mises por Lamé de pared gruesa, deformación radial. Declarar: el factor de transmisión de poro `P_poro = P_ext × 0.85 × 0.03` son constantes puestas a mano; el mapa de Von Mises del andamio (Panel 7) es **un campo pintado con gaussianas** (`S = 15·(1−exp) + 35·exp(−d²/0.25)`, escalado al valor de Lamé), no un cálculo espacial.

**Co-SWIFT / optimizador** (`optimizador_coswift.py`).
Real: reología Herschel-Bulkley, Hagen-Poiseuille modificado, MOPSO con archivo de Pareto externo y crowding distance. Declarar: las funciones objetivo (`cell_viability_model`, `vascular_resolution_model`) son **curvas empíricas a trozos con puntos de quiebre elegidos, no ajustadas a datos**. Además `cell_viability_model` está **clampeada a [20, 98]** — el "98% de viabilidad" del dashboard es el techo de la función, no un óptimo calculado. La maquinaria es real; lo que optimiza es semi-fenomenológico.

**Filtración glomerular** (`simulador_filtracion_glomerular.py` / `_G`).
Real: ecuación de Starling-Deen bien implementada. Declarar: está calibrada para dar en el blanco. La versión `_G` aplica **+15 mmHg de "vasodilatación aferente"** dentro de una función de auto-ajuste cuando la presión entra baja; la curva de "autorregulación" es un mapeo lineal prescrito `Pgc = 0.60 × P_art`, no autorregulación miogénica/TGF; el "mapa de presión CCO" son posiciones `np.random`. **Inconsistencia interna a corregir:** los dos archivos usan un `Kf_glom_nL` distinto — 3.7 (no-G, "calibrado exacto") vs 4.1 (`_G`, "ajuste técnico"). Hay que **elegir uno y justificarlo**; el resultado de TFG depende de cuál se corra.

**Reabsorción tubular** (`simulador_reabsorcion_tubular.py`). **EN CUARENTENA — ver §3.4.**
Real (en diseño): cascada de 5 segmentos anatómicamente fiel, con formas de transporte reales (Michaelis-Menten, Kedem-Katchalsky, contracorriente). Declarar: las fracciones de reabsorción están casi todas fijadas a mano a valores de libro (`0.67` agua/Na proximal, `vol *= 0.926` por paso "para ~73%"), y `n_ok = 6`, `funcion_global = 98.0` están **hardcodeados** (el propio comentario dice "Simulado para brevedad del refactor").

### Nivel 3 — Fenomenológico o decorativo (ilustrativo, sin validez predictiva)

- **iPSC simple** (`simulador_ipsc_biokidney.py`): tasas "basadas en Vibe Coding" (comentario del propio autor: inventadas).
- **WSS simple** (`simulador_wss_biokidney.py`): usa presión absoluta en vez de caída de presión, con viscosidad definida pero sin usar. Dimensionalmente mal aplicado.

**Caso aparte — diferenciación iPSC con EDO** (`simulador_diferenciacion_ipsc.py`).
Mejor que el simple: sistema de EDOs real con señalización tipo Hill, referenciado a Takasato 2015 y Freedman 2015. Pero **todas** las constantes cinéticas (`k_iPSC_to_IM`, etc.) están elegidas, no ajustadas a datos. Produce curvas con forma de diferenciación; las cifras de pureza y teratoma **no tienen calibración empírica**.

---

## 3. Correcciones transversales obligatorias

### 3.1 Teratoma y pureza: eliminar como resultados
`teratoma = ipsc/total` y `purity = target/total` salen de una EDO con parámetros inventados. Ninguna ecuación in silico puede "validar" que el riesgo de teratoma es bajo ni la pureza del 100%. En laboratorio real los organoides renales de iPSC no se acercan al 100% y el teratoma es un riesgo abierto. **Reformular como objetivos de diseño del protocolo tomados de la literatura (Takasato/Freedman), nunca como resultado de seguridad biológica.** Es la afirmación más peligrosa del deck.

### 3.2 Calibración declarada, no oculta
Varios módulos se afinan para superar umbrales clínicos (TFG>60, viabilidad>85%, orina ~1.5 mL/min) y luego se reporta el haberlos superado como hallazgo — lo cual es circular. La calibración es legítima en modelos de orden reducido; **ocultarla es lo letal**. Un `+15 mmHg` escondido en una función llamada "auto-ajuste" se lee como fraude cuando un revisor lo encuentra. Declarar explícitamente y justificar: el +15 mmHg, la elección de Kf (3.7 vs 4.1), las fracciones de reabsorción, el 0.926/paso.

### 3.3 Constantes de entrada, no convergencias
`GFR_INPUT = 82.0` ("RESULTADO PREVIO"), `WSS = 5.6 dyn/cm²` ("validado") y la orina ~1.5 mL/min **no emergen de un pipeline acoplado; están fijados**. Peor: el 82 contradice a sus dos supuestos productores (ambos glomerulares apuntan a ~62.5) y el módulo que lo consumiría no corre. Dejar de presentarlos como si el sistema los produjera de forma consistente; etiquetarlos como constantes de entrada.

### 3.4 Reabsorción: reparar o poner en cuarentena antes de citar
El archivo **no parsea** (IndentationError, línea 208: bloque de diccionario huérfano duplicado tras el `return` de la 207 — verificado con `py_compile`). Como no se ejecuta, **la imagen del Dashboard Maestro con cifras de reabsorción (2.19 L/día, 98.1%, 6/6) no pudo salir de este código**: la procedencia imagen↔código está rota. Antes de que cualquier cifra de reabsorción entre a un artefacto: (a) reparar el syntax; (b) al hacerlo aparecen KeyErrors porque el dashboard referencia claves que `simular()` nunca crea (`funcion_nativa`, `checks`, `orina['Na']/['K']/['glucosa']/...`); (c) dado que las fracciones son de libro puestas a mano, **la jugada honesta es cuarentenar y reetiquetar como orden-reducido-calibrado, no parchear claves para pintar de verde**; (d) confirmar qué código generó realmente la imagen del dashboard.

### 3.5 Lenguaje de sobreafirmación: quitar de todo artefacto
Eliminar "VALIDACIÓN COMPLETA / 100% completado / VALIDACIÓN FINAL" y cualquier invocación a instituciones (Harvard Wyss / Oxford, etc.) que no correspondan a una colaboración real. Es justo el tipo de sobreafirmación que termina una conversación con un revisor serio y que erosiona la credibilidad, único capital de un investigador independiente.

---

## 4. Reemplazo del estado del Dashboard

Sustituir el bloque "OK · RIÑÓN FUNCIONAL · VALIDACIÓN IN SILICO COMPLETA · Pipeline 100% · 12/12 módulos ÓPTIMOS" por un **estado por niveles**, que etiquete cada módulo como:

- `SOLVER` — solver de campo real (solo O₂), con sus salvedades de interpretación.
- `ORDEN REDUCIDO (CALIBRADO)` — dECM, Co-SWIFT, filtración, reabsorción (esta última: `EN CUARENTENA`).
- `ILUSTRATIVO` — iPSC simple, WSS simple.

La honestidad de nivel **es** lo que da credibilidad. El estado global deja de ser "100%" y pasa a describir qué está validado geométricamente vs qué es estimación de diseño.

---

## 5. Qué NO cambia (lo validado de verdad)

El cuerpo científico permanece intacto y es la contribución real:
- Gemelo geométrico multicapa (Capas 0–4): dominio anclado a morfometría (cortical width MDCT n=2068), 1300 nefronas representativas, árbol arterial con Murray k=3.0000 y 100% de glomérulos alcanzados, sistema colector y calicial anclado a literatura (infundíbulo ~4 mm CTU n=1321, pelvis normal, conducto de Bellini), contención y conectividad auditadas.
- ~~Difusión de O₂ como único solver de campo real.~~ **Retirado 2026-07-13:** el solver de O₂ diverge y está EN CUARENTENA (§2, §6). La geometría (Capas 0–4) permanece válida y la auditoría confirmó que **no depende del campo de O₂ en ningún punto** — el bug no toca ni un vóxel del gemelo geométrico.
- Escala representativa y sus límites (radio raíz arterial, lecho capilar ausente por diseño, cruce hilar a escala idealizada) documentados como filas de correspondencia, no como defectos.

La fabricación imprimible (Capa 5, malla/voxelización) es un entregable posterior y separable; su cruce de estanqueidad hilar es un límite de escala del hilio idealizado, **no un bloqueo del in silico** — el multisimulador consume un grafo de segmentos, no una malla sólida.

---

## 6. Hallazgos de la auditoría del solver de O₂ (2026-07-13)

Registro derivado de BITACORA Entrada 028. Solo lectura; ningún fix aplicado al código.

### 6.1 Tres modelos de O₂ coexistiendo (deuda arquitectónica)
El framework contiene **tres** implementaciones distintas de difusión de O₂, que no coinciden entre sí ni con el texto publicado:
1. **Fick 3D 60×60×40** (`simulador_oxigeno_biokidney.py`) — el "solver real"; **roto** (diverge a NaN — §2).
2. **Krogh 2D** del backend web (`web_app/backend/services/simulation_service.py`) — superposición de cilindros + suavizado, **acotado** (no diverge); es el que usa la app. Comparte `P_HIPOXIA=1.0` y el consumo, pero es otro modelo.
3. **Krogh 30×30** del texto del preprint — descrito en §2.2 de Zenodo v3; **sin código versionado que lo produzca**.
Deben unificarse a una sola fuente de verdad antes de citar cualquier cifra de O₂.

### 6.2 Contradicción interna del preprint v3 (Zenodo) — a corregir
- **§2.2 (Métodos)** describe un **cilindro de Krogh en grilla 30×30** (2D, capilar único, 3 iteraciones, umbral 4 mmHg).
- **§3.2 (Resultados)** describe un **Fick 3D reacción-difusión sobre la geometría real, resuelto "a convergencia"** (el "único solver de campo real").
- **Ambas citan las mismas cifras:** min 5.6 mmHg, media 31.76 mmHg, 0% hipoxia.
- Esas cifras **no son trazables a ningún código versionado:** el Fick 3D diverge (98% NaN); el único artefacto finito daba 520 vóxeles hipóxicos (no cero); el umbral del código es `P_HIPOXIA=1.0` (no 4); la grilla es 60×60×40 (no 30×30).
- El número está **hardcodeado como string** en el dashboard (`06_app/dashboard_maestro_app.py:8-9`: "PO2 minima 5.6 mmHg sobre umbral critico 4 mmHg. 0% hipoxia").
- **Líneas a corregir en v3:** 26 (abstract, "solved to convergence"), 106 (título "Krogh Cylinder Model"), 108 (30×30 / 3 iter / 5.6 / umbral 4 / zero hypoxic), 210 (5.6 / 31.76 / "to convergence" / "no voxels below threshold").

### 6.3 Hallazgo de método: etiquetas que prometen más de lo que hace el código
Patrón sistemático — nombres/etiquetas que afirman una semántica que el código no cumple:
- **`sistema=col`** en un CSV *vascular* es en realidad el **sistema colector urinario** (pelvis/cálices), no una estructura vascular; se colaba al campo de O₂ como fuente sanguínea espuria a 30 mmHg.
- **`senyal_hipoxia_para_cco.csv`** promete un bucle de realimentación hipoxia→CCO que **nadie lee**; el árbol nunca creció con esa señal.
- **"v8"** nombra **dos geometrías distintas** (`arbol_vascular_cco_v8.csv`, 1902 seg, vs `arbol_vascular_cco_v8_gemelo.csv`, 38 015 seg).
- El `return` por defecto de `pres()` **asigna 30 mmHg a cualquier etiqueta desconocida** en vez de rechazarla — así entró `col` al término fuente.

**Regla derivada (canónica):**
1. **Los nombres deben describir lo que el código hace,** no lo que se aspira a que haga. Un archivo `_para_cco` que nadie consume, un `sistema=col` que no es vascular y un "v8" ambiguo son deuda de trazabilidad que hay que saldar.
2. **Todo `return` por defecto ante una etiqueta desconocida debe fallar ruidosamente** (excepción / log de error), **nunca asignar un valor silenciosamente**. El mismo principio aplica a `nanmean`/`nanmin`: enmascarar NaN es la versión numérica del mismo pecado — un cálculo que colapsó debe gritar, no imprimir "100% oxigenado".

---

*Documento canónico. Cualquier cifra que entre a Zenodo v3, al Dashboard o a material de campaña debe ser trazable a este marco: solver / orden-reducido-calibrado / ilustrativo, con toda calibración declarada.*
