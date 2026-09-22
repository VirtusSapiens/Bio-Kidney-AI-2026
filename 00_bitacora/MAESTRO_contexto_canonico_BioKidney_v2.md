# BIO-KIDNEY AI — DOCUMENTO MAESTRO DE CONTEXTO (CANÓNICO) · v2

**Fuente única de verdad para GENESIS y los tres agentes Nexus.**
**VirtusSapiens · Carlos David Moreno Cáceres · Medellín, Colombia**

**Estado de referencia: 22 de agosto de 2026**
**Supersede: MAESTRO v1 (19 jun 2026) e INDICE_BioKidneyAI.md (mayo 2026).**
**Cifras y título verificados por Carlos contra `preprint_biokidney_2026_EN_v4.md`
y la BITÁCORA desde su propia terminal (22 ago 2026), conforme a §7.**

> Si cualquier otro documento, instrucción de agente o memoria de proyecto
> contradice este archivo, **prevalece este**. El MAESTRO v1 debe archivarse,
> no consultarse.

---

## 0. AFIRMACIONES RETIRADAS — NO REPETIR NUNCA

Estas afirmaciones circularon en versiones anteriores del MAESTRO, del INDICE o
de las instrucciones de los agentes. **Todas están retiradas.** Si aparecen en
un borrador, es un error de fuente, no una omisión a rellenar.

| Afirmación retirada | Estado real |
|---|---|
| "Oxigenación tisular simulada / solver de O₂ validado" | **Retirado del preprint v4.** El solver diverge a NaN (Jacobi + ω=1,6 inestable; ~98% del dominio NaN enmascarado por nanmean/nanmin). No se comunica ningún resultado de O₂. |
| "Producción de orina 2,19 L/día" / "módulo de reabsorción validado" | **Retirado.** El módulo no compila (IndentationError). Cifra no computada. |
| "Capa 3 es la frontera, aún sin empezar" | Falso desde junio. Capas 0–4 **cerradas**. |
| "La mejora de cobertura que traerá Capa 3 será un resultado central del paper" | **Falso.** Ver §6: la cobertura a 150 µm NO es métrica de validación de Capa 3. |
| "Figura 4 / arquitectura MoE / TFG asociada a esa figura" | **Eliminada** del preprint v4. |
| "7/12 KPIs superan al nativo" · "supera al riñón nativo" | Superadas. Nunca. |
| **"TFG 115,2 mL/min" y "TFG 115,4 mL/min"** | **Ambas superadas.** El v4 reporta `~115 mL/min` sin decimal, deliberadamente (§4). Nunca citar un decimal. |
| **"Pipeline integrado predictivo" / "los cinco módulos integrados en un solo pipeline"** | **Sobreclaim.** El módulo de filtración corre sobre el árbol CCO **v7 legacy**; la geometría v8 aún no está conectada a los simuladores (v4 §2.4). El ~115 es un chequeo aritmético a mano, no salida de una cadena integrada. |
| "TFG 82 mL/min" como resultado actual | El 82 **sí existe** en el v4, pero solo como **constante de entrada** del módulo standalone bajo presión sintética, retenida para documentar el cambio metodológico v7→v8. Nunca presentarlo como TFG del proyecto. |
| "115 ± 13,4 mL/min" o cualquier ± aplicado a la TFG | El ±13,4 mmHg es la **desviación espacial de la presión terminal**, no dispersión de la TFG (v4, Fig. 3). |
| "WSS 24,69–37,04 dyn/cm²" (rango 10–70) | Superada. Rango renal real 1–10. |
| "CCO v7 / 1.448 segmentos" como cifra actual | Histórica. |
| "Primero en modelar vasculatura renal" | Insostenible. Campo activo y poblado. |
| "6 simuladores validados" (INDICE, mayo) | Falso. Ver §4: tres niveles epistémicos, no "validado" uniforme. |
| "renal_data_v1.json aún no exportado a disco" (INDICE) | Resuelto. Ver §5. |

---

## 1. INVESTIGADOR (identidad base neutra)

**Carlos David Moreno Cáceres** — venezolano, residente en Medellín, Colombia.
Investigador independiente en bioingeniería computacional, autodidacta, sin
afiliación institucional. Marca científica: **VirtusSapiens** (*Virtus et
Sapientia*).

Formación previa en la Fuerza Aérea Venezolana (ex oficial, guerra electrónica
y aeronáutica, 2010–2018). Esa formación es la **base metodológica** declarada:
análisis de señales, modelado de sistemas complejos, decisión rigurosa bajo
incertidumbre.

En registro científico esta es la única identidad que se usa: **investigador
independiente**. La narrativa personal de salud se activa solo en Divulgación y
Fundraising (§9).

Colaborador técnico: **John Tapias** (VECANOVA) — ingeniero de sistemas; co-desarrolló la SPA web.
ORCID: **0009-0005-3933-5072**.

---

## 2. MARCO DE HONESTIDAD EPISTÉMICA (regla que gobierna todo lo demás)

Este marco (`MARCO_honestidad_epistemica_BioKidney.md`) es la columna vertebral
de cómo se comunica el proyecto. **Ningún agente puede describir un resultado
sin declarar a qué nivel pertenece.**

**Nivel 1 — Solver de campo real.** Resuelve ecuaciones sobre un dominio.
Actualmente: **ninguno activo.** El único que existía (O₂) está retirado
pendiente de reparación.

**Nivel 2 — Verificación de factibilidad de orden reducido, calibrada.** No
predice; comprueba que un valor cae en rango fisiológico bajo una calibración
declarada. Aquí viven TFG, Pgc, WSS, viabilidad Co-SWIFT. **Nunca llamarlos
"predicción" ni "resultado validado".**

**Nivel 3 — Ilustrativo.** Sirve para visualizar o comunicar. No es evidencia.

### El retiro de O₂ no daña la contribución (dato comunicable)

El v4 documenta un **rastreo de dependencias** que confirma que ningún módulo
geométrico consume el campo de oxígeno: el generador vascular, la morfometría
poblacional y las Capas 0–4 se computan de forma independiente. La
retroalimentación hipoxia→CCO que sugería el nombre de archivo
`senyal_hipoxia_para_cco.csv` **no existe en el código** — esa señal se escribe
pero nunca se lee.

Esto se comunica así: *"retiramos el módulo de oxígeno y verificamos, con un
rastreo de dependencias, que la contribución geométrica queda intacta."* Es un
ejemplo de autocorrección rigurosa, y funciona mejor como argumento de
credibilidad que cualquier cifra que el solver hubiera producido.

### El reencuadre central (usar esta frase como base)

> **La fidelidad anatómica y geométrica del gemelo digital es el resultado
> científico validado. Los módulos funcionales son estimaciones de
> factibilidad explícitamente calibradas.**

Esta es la afirmación defendible del proyecto. Es más modesta que la de marzo y
mucho más sólida.

---

## 3. ESTADO CANÓNICO DEL PROYECTO (a 22 ago 2026)

Bio-Kidney AI es un **gemelo digital geométrico multicapa del riñón humano**,
enteramente *in silico*. El riñón físico **no existe**. El propósito declarado
es la **validación in-silico previa a fabricación** para bioimpresión.

**Hito actual:** las Capas 0 a 4 del gemelo están construidas, auditadas y
cerradas. La Capa 5 (mallado de lumen por SDF) está **formalmente diferida**
como entregable de fase de fabricación — es una decisión, no un pendiente.

**Cómo describir el momento actual (texto aprobado):**

> "El gemelo digital geométrico está completo desde el dominio del órgano hasta
> el sistema colector: nefronas, campo de demanda, árbol arterial, red
> peritubular, árbol venoso, árbol colector urinario y cálices. Cada capa está
> anclada a morfometría poblacional publicada y auditada contra literatura. La
> siguiente fase es una pasada de calibración, no una capa nueva."

**Auditoría de correspondencia anatómica:** 47 filas — 24 Clase A (afirma
anatomía) y 23 Clase B (simplificación declarada). Que existan 23
simplificaciones declaradas es una fortaleza del trabajo, no una debilidad;
comunicarlo así.

---

## 4. CIFRAS CANÓNICAS (usar SOLO estas, con su nivel epistémico)

### Geometría y morfometría — Nivel 1 de confianza (medido sobre el modelo)

| Parámetro | Valor | Nota |
|---|---|---|
| Nefronas | 1.300 | 85% corticales / 15% yuxtamedulares |
| Pirámides medulares | 10 | Ancla: Bonsib / Heptinstall |
| Corrección de profundidad cortical | método B (perpendicular a cápsula) | Reclasificó 25.079 puntos (12,54%) de médula a córtex |
| Volumen excavado del seno renal | ~3,96 mL | Por debajo del límite inferior de Caglar 2014 → **seno declarado subdimensionado** |
| Cobertura del eje largo por pirámides | 45,1% | 47,1% de la médula no pertenece a ninguna pirámide |
| Ápices piramidales en territorio cortical | 4 | Hallazgo de auditoría, declarado |
| Anclas de espesor | CW 6,6 mm · PW 15,4 ± 2,8 mm | Glodny et al. 2009 |
| Árbol vascular v8 | 1.902 segmentos (904 art / 926 ven / 72 col) | Determinista, seed 42, md5 `4f881a05cbbf8f5713d416c61bd1ceaa` |
| Puntos de demanda alcanzados | 100% | Métrica de validación real de Capa 3 |
| Ley de Murray | k = 3,0 | Cumplida |
| Puntos de drenaje peritubular | 4.290 | Capa 3a-bis |
| Confluencias Y venosas | 1.527 (79,2% de 1.927) | Anatomía correcta, no defecto |
| Cálices / pelvis (Capa 4) | 24 nodos, 23 aristas | Pelvis elipsoide achatada |

### Módulos funcionales — Nivel 2 (factibilidad calibrada)

| Parámetro | Valor | Encuadre obligatorio |
|---|---|---|
| TFG bilateral | **~115 mL/min** (sin decimal) | *Feasibility check* bajo escala calibrada, **NO** predicción geométrica de TFG. Es K_f de literatura (3,7) × ΔP Starling (15,6 mmHg) **calculado a mano**. En rango normal **por construcción de la calibración**. |
| Presión capilar glomerular (Pgc) | 58,6 ± 13,4 mmHg | Por encima del umbral de filtración (~43 mmHg). |
| WSS hemodinámico | 5,6 dyn/cm² | **Punto de operación de entrada** (Co-SWIFT), no resultado. Rango renal real 1–10. |
| Viabilidad Co-SWIFT | 98% | **Techo (clamp) de la función objetivo [20,98]**, no óptimo computado. |
| Ventana de bioimpresión iPSC | Días 21–30 de cultivo | — |
| Oxigenación tisular | **RETIRADA** | Solver diverge. No comunicar cifra alguna. |
| Producción de orina | **RETIRADA** | Módulo no compila. |

> **Por qué "~115" y no un decimal.** El 115,2 fue el valor del v3; el 115,4 fue
> una corrección aritmética posterior (BITÁCORA 1447). El v4 abandonó ambos y
> adoptó `~115` porque la cifra es un **chequeo aritmético hecho a mano** que
> toma 58,6 mmHg como constante de entrada — un decimal le daría falsa
> apariencia de salida de solver. **Un agente que escriba 115,2 o 115,4 está
> citando una versión retirada.**
>
> **Cadena real (decirlo así si alguien pregunta):** el generador de geometría
> v8 produce el campo de presión terminal; el módulo de filtración corre
> aparte sobre el árbol v7 legacy y no consume ese campo. El ~115 es un
> chequeo de consistencia entre ambos, no una salida acoplada.

### Hardware
Dell Inspiron 14-3467, i5, 16 GB RAM, sin GPU. Ubuntu 24.04, Python 3.12,
entorno `env_biokidney`.

---

## 5. ARQUITECTURA DEL GEMELO DIGITAL (Capas 0–5)

Orden de dependencia: forma → nefronas → demanda → árboles → colección →
cálices → (lumen).

- **Capa 0 — Dominio** ✅ **auditada.** Riñón en forma de frijol, seno esculpido, hilio, partición córtex/médula, 10 pirámides. Sistema de coordenadas del proyecto.
- **Capa 1 — Nefronas** ✅ 1.300 nefronas, split 85/15. *Nota: el sesgo polar de 4,8× reside en la función de asignación de esta capa, no en Capa 0.*
- **Capa 2 — Campo de demanda** ✅ Campo escalar de demanda de O₂.
- **Capa 3a — Árbol arterial** ✅ CCO binario, Murray k=3,0, 100% de glomérulos alcanzados.
- **Capa 3a-bis — Red peritubular / vasa recta** ✅ 4.290 puntos de drenaje.
- **Capa 3b — Árbol venoso** ✅ **CERRADA tras cinco auditorías.** Los 1.927 shunts VV a nivel de calibre se descompusieron en 79,2% confluencias Y anatómicamente obligatorias, 19,6% reenrutables (cota superior laxa) y 1,1% ambiguos. **Veredicto: no es un defecto estructural.** La vascularización venosa **pasa** la compuerta de factibilidad anatómica. *Artefacto conocido: todos los segmentos venosos miden exactamente 0,6 mm (paso de crecimiento fijo).*
- **Capa 3c — Árbol colector urinario** ✅ Bosque de 10 subárboles, uno por pirámide. Radios de conductos de Bellini anclados. Sin Murray ni CCO (convergente).
- **Capa 4 — Cálices y pelvis** ✅ 24 nodos, 23 aristas.
- **Capa 5 — Mallado de lumen (SDF)** ⏸️ **DIFERIDA por decisión**, entregable de fase de fabricación.

**Principios de arquitectura:** el modelo es un **grafo** (nodos + aristas), no
vóxeles. Regla de escritor único: la ciencia vive en `env_biokidney`
(NumPy/SciPy) y produce `.npz`; **Blender solo lee y dibuja**. Nunca física
dentro de `bpy`.

**Datos vasculares en disco:** `renal_data_v1.json` corregido (861 KB) ya
reemplazó al stub de 3 segmentos que se había commiteado por error. El INDICE de
mayo, que marca ese archivo como "no exportado", está obsoleto en ese punto.

**Límite de rol:** el trabajo técnico (código, pipeline, paper, BITÁCORA) se
hace en Claude Code, **no** dentro de los agentes de marca. GENESIS explica esta
ciencia; no la desarrolla.

### ✅ Higiene documental de `00_bitacora/` (cerrada 2026-09-21)

Los dos archivos que esta sección marcaba como la contaminación más seria
del repo ya no están en rutas activas:

- `implementation_protocol_v1.md` → `99_archivo/SUPERSEDED_implementation_protocol_v1.md`
  (commit 51a3944).
- `BITACORA_sesion_resometimiento_biorxiv.md` → `99_archivo/SUPERSEDED_BITACORA_sesion_resometimiento_biorxiv.md`
  (fuera de rutas activas; no se commitea porque nunca estuvo en git).

La limpieza se extendió al resto del repo: PDFs de preprints superados,
informes y capturas de simuladores, materiales de presentación, dashboard
maestro y simuladores de reabsorción y de O₂ (commits 508e97c, 8dc4190 y
4eb7029). Detalle en BITÁCORA, entrada del 2026-09-21.

**Verificar fuera del repo:** si el conocimiento de GENESIS es una copia
subida aparte y no una lectura directa de `00_bitacora/`, estos archivos
deben retirarse también de esa copia. Archivar en el repo no los borra de
ahí.

### Frentes abiertos tras la higiene documental

- **v5 del depósito en Zenodo:** Figura 3 (la imagen no corresponde a su
  caption y tiene encuadre de mejora), Figura 2 (el mínimo de presión
  coincide con el umbral), título retirado en el registro Zenodo y en
  `ORCID/works.bib`, y decisión sobre `supplementary_material_v8.md`.
- **`web_app/`:** sigue calculando métricas de O₂, de flujo urinario y el
  porcentaje de TFG nativa.
- **Simuladores de filtración glomerular**
  (`simulador_filtracion_glomerular{,_G}.py`): comparaciones con el riñón
  nativo en código activo; el v4 cita el módulo standalone (§2.4).
- **Módulos iPSC y dECM:** veredictos que contradicen sus propios datos.
- **`INDICE_BioKidneyAI.md`:** superado, pendiente de regeneración.
- **Fuera del repo:** ubicar `biokidney_fixed.jpg` (exportación del
  dashboard maestro) y verificar si se publicó.

---

## 6. HALLAZGO DE COBERTURA DE DIFUSIÓN (leer antes de comunicar cualquier %)

Con vasculatura completa (arterial + peritubular + venosa), la cobertura de
difusión a 150 µm es **0,78%** (línea base sin árbol: 0,66%).

**Esto NO es un fallo del árbol.** Es un límite de escala de resolución: los
segmentos terminales miden 0,6 mm — cuatro veces el radio de difusión de 150 µm
— y el lecho capilar que cerraría esa brecha **no está en la geometría por
diseño** (los árboles son troncos de suministro y drenaje hacia 4.290 puntos
representativos, no capilares individuales). Barrido de umbral: 150 µm → 0,8%;
500 µm → 9%; 1000 µm → 28%. Escala como problema de resolución, no como brecha.

### Reglas duras para agentes

1. **La cobertura a 150 µm NO es métrica de validación de Capa 3.** Es una
   métrica de la futura etapa capilar. No usarla como titular, ni positivo ni
   negativo.
2. **La validación de Capa 3 es geométrica:** 100% de puntos de demanda
   alcanzados, Murray k=3,0, morfometría dentro de rango.
3. Existe además un **vacío avascular central genuino** de ~5 mm de radio (vaso
   más cercano a 4,9 mm). Está localizado, es distinto del déficit global de
   resolución, y **no es el seno renal** (0 demanda cae en el seno; 100% del
   déficit está en parénquima). Hipótesis en estudio: vinculado al sesgo polar
   de drenaje conocido. Pendiente para la pasada de calibración.
4. Ningún agente comunica el vacío central como resuelto. Está **abierto**.

---

## 7. REGLA DE VERIFICACIÓN (no negociable)

> **Ningún reporte de un agente cuenta como verificación. Un resultado existe
> solo cuando Carlos lo confirma desde su propia terminal.**

Origen: la BITÁCORA documenta (entrada 030 y posteriores) episodios de
fabricación por parte de asistentes de IA — propagación inventada de números de
línea, afirmaciones falsas de "salida cruda pegada arriba", escalamiento de
falsos positivos sin verificar docstrings, fuentes inventadas y confusión
repetida de métricas (radial vs. perpendicular).

**Corolario para contenido:** ninguna referencia bibliográfica sale de un agente
a una pieza publicable sin que Carlos la verifique en la fuente. Un agente que
no tiene la cita **dice que no la tiene**; no la reconstruye.

**Anclas de literatura confirmadas y utilizables:**
Glodny et al. 2009 (CW 6,6 mm; PW 15,4 ± 2,8 mm) · Caglar 2014 (volúmenes de
grasa sinusal) · Bonsib / Heptinstall (N = 10 pirámides) · Takasato (cinética
iPSC) · Starling-Deen (filtración) · Murray (ley de ramificación).

---

## 8. POSICIONAMIENTO HONESTO

- El aporte es **integración y fidelidad geométrica auditada**, no invención de métodos.
- El gemelo resuelve **"el plano, no el ladrillo"**: la mitad geométrica y de diseño del cuello de botella de vascularización. La mitad biológica —maduración, anastomosis, perfusión sin necrosis— requiere laboratorio húmedo y **ningún modelo digital la resuelve solo**.
- Limitación explícita (Fase 1): trabajo **completamente in silico**, sin experimentos in vitro ni in vivo.
- La geometría del dominio es **idealizada paramétrica** (morfometría poblacional), **no** reconstrucción de un riñón real por imágenes médicas. Decirlo así fortalece la credibilidad.
- Estrategia de tallas S/M/L: envolvente del órgano y número de nefronas escalan; la anatomía fina se mantiene fija; la TFG escala con superficie corporal. **El modelo M debe estar completamente validado antes de parametrizar S y L.**

### Título del preprint (verificado contra v4)

**Título vigente:**
> *A Multi-Layer Geometric Digital Twin of the Human Kidney: Anatomical Fidelity
> and Calibrated Reduced-Order Feasibility Checks*

El título anterior — *"…Six-Module Integration Predicts Physiological Renal
Output…"* — fue **retirado en el v3** (CHANGELOG, fila 1), con el razonamiento
correcto: se defiende geometría, no función; el encuadre predictivo es
indefendible. **Ningún agente cita el título viejo.** El MAESTRO v1 lo
arrastraba en su §7; ese es uno de los errores que este v2 corrige.

**Conteo de módulos:** cinco es el canónico del v4 (generador vascular +
filtración glomerular + Co-SWIFT + cinética iPSC + reabsorción tubular). El
paso `six→five` fue deliberado al retirar O₂. No mezclar ambos números en una
misma pieza.

---

## 9. NARRATIVA PERSONAL (solo Divulgación y Fundraising)

*No se usa en registro científico ni en alianzas técnicas.* Esta separación es
un estándar editorial deliberado de Carlos.

**Frase núcleo:** *"No es mi enfermedad la que me define. Soy yo, con mis
acciones, a pesar de ella."*

**Tono obligatorio:** admiración e inspiración, nunca lástima. Protagonista, no
víctima.

**Material estrictamente privado:** hay vivencias que Carlos compartió como
contexto y que **no** forman parte de la narrativa pública. No se usan, no se
referencian, no se convierten en contenido. Carlos es el único que decide qué
cruza a lo público.

> ⚠️ **Inconsistencia estructural pendiente de resolver.** Nexus Copy y Nexus
> Creator tienen la narrativa personal escrita dentro de su sección de identidad
> base, lo que contradice esta regla de activación selectiva. Debe resolverse en
> una de dos direcciones (decisión de Carlos): (a) mover la narrativa a una
> sección de registro activable en esos dos agentes, o (b) modificar esta §9
> para reconocer que en los canales de marca personal la narrativa es base. Hoy
> los documentos se contradicen.

---

## 10. CREDENCIALES VERIFICABLES

- **ORCID:** 0009-0005-3933-5072 — https://orcid.org/0009-0005-3933-5072
- **DOI (concept, todas las versiones):** 10.5281/zenodo.19508076 — **usar este por defecto**; nunca queda obsoleto.
- **DOI de versión:** 10.5281/zenodo.19508077 — solo para referir una versión exacta.
- **Repositorio:** github.com/VirtusSapiens/Bio-Kidney-AI-2026
- **Repo local:** `~/Escritorio/BioKidney-AI/` · BITÁCORA en `00_bitacora/BITACORA.md`
- Licencia CC-BY 4.0.

---

## 11. CONTACTOS ACADÉMICOS

- **Dr. José David Hincapié** (GIB — UDEA): abierto a asesorar; posible endoso/afiliación para resometer a bioRxiv. Contacto prioritario.
- **Dr. José Nelson Carvajal Quiroz**: referente de nefrología en Colombia/LATAM; primer contacto pendiente.
- **Dr. Bustamante** (UPB): escéptico.
- **Prof. Jhon Freddy Ochoa** (UDEA): pendiente vía Hincapié.
- Objetivos internacionales (largo plazo): Wyss Institute, IBME Oxford, Wake Forest, McGowan.

---

## 12. QUÉ NO AFIRMAR NUNCA (guardrails)

1. Que el riñón físico existe o está "construido". Es un gemelo digital.
2. Que el proyecto "supera al riñón nativo" o fue "primero en modelar vasculatura renal".
3. Cualquier cifra o afirmación del §0.
4. Ningún resultado de oxigenación tisular ni de producción de orina.
5. Cobertura de difusión a 150 µm como validación o como fallo (§6).
6. Que un resultado está verificado sin confirmación de Carlos desde su terminal (§7).
7. Cualquier referencia bibliográfica no verificada por Carlos.
8. Consejo médico a otros pacientes renales (dietas, dosis, valores, tratamientos). Carlos comparte experiencia y señala recursos; el manejo clínico es del equipo médico de cada persona.
9. Que un grant o programa existe o aplica sin verificar su vigencia al momento de postular.
10. Que la Capa 5 está pendiente por falta de tiempo. Está **diferida por decisión**.

---

**VirtusSapiens — Virtus et Sapientia**
Carlos David Moreno Cáceres · Medellín, Colombia · Documento canónico v2, agosto 2026
