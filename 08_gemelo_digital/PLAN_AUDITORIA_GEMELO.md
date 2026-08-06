# Plan de Auditoría del Gemelo Digital — Bio-Kidney AI

**Propósito.** Recorrer el gemelo digital de capas (Capa 0 → 4), capa por capa,
verificando que cada parámetro geométrico esté (a) anclado a literatura, (b)
declarado explícitamente como supuesto, o (c) marcado como no verificado. El
objetivo es que el gemelo salga como preprint propio —separado de la rama
CCO v8— solo cuando cada capa haya pasado esta puerta.

**Regla de oro del régimen nuevo.** Verificar ANTES de generar/simular, no
después. Un número sin ancla no entra, aunque produzca una figura bonita. Un
resultado real y malo es un hallazgo; uno falso y bonito es deuda.

**Qué NO es este documento.** No es una auditoría ya ejecutada. Es el mapa de
qué verificar. Cada celda de estado se llena con evidencia producida por Carlos
desde su propia terminal (regla: ningún reporte de agente cuenta como
verificación sin confirmación independiente).

---

## Leyenda de estados

| Símbolo | Significado |
|---|---|
| ✅ ANCLADO | Parámetro con cita de literatura verificada y valor coherente |
| 📐 SUPUESTO | Decisión de diseño declarada explícitamente; sin pretensión de exactitud anatómica, pero razonable y documentada |
| ⚠️ SIN VERIFICAR | Aún no auditado; estado desconocido |
| 🔴 DEFECTO | Problema confirmado que requiere corrección |
| 🟡 FIDELIDAD | Funciona en rango/topología pero la geometría no refleja bien la anatomía |
| 🚫 N/A | No aplica a esta capa |

**Tier de cita** (heredado de `auditoria_correspondencia_anatomica.md`):
A = requiere ancla cuantitativa de literatura; B = cualitativo/estructural.

---

## Separación de ramas (contexto crítico)

Antes de auditar, fijar la frontera que motivó todo esto:

- **Rama CCO** (`02_vascular_cco/generador_cco_v8.py` + multisimuladores):
  objeto propio. Genera su demanda con muestreo Beta(3,1.2) sintético. **NUNCA
  abre `capa1_nefronas.npz`.** Sus cifras (1902 seg, Murray α=3.0, presiones,
  GFR) son independientes del gemelo. → Se documenta y cierra por separado.
- **Rama Gemelo Digital** (`08_gemelo_digital/`, Capas 0–4): objeto de este
  plan. Aún en construcción/verificación. → Sale como preprint propio cuando
  pase esta auditoría.
- **Prohibición de aquí en adelante:** no volver a mezclar ambos objetos en un
  mismo documento/versión. El gemelo no es "v(n+1)" de la rama CCO.

---

## Capa 0 — Dominio (parénquima)

Archivo: `capa0_dominio.py` → `capa0_dominio.npz` (coords, region_label)

| Parámetro | Valor actual | Estado | Ancla / Nota |
|---|---|---|---|
| Forma envolvente (elipsoide) | semiejes X/Y/Z | ⚠️ SIN VERIFICAR | ¿Dimensiones ancladas a morfometría poblacional? |
| Umbral córtex/médula | 6.6 mm (corrección Entrada 014) | ⚠️ SIN VERIFICAR | Anclado a MDCT n=2068 según bitácora — reconfirmar cita |
| Definición del seno renal | (24,18,12)? | 🔴 DEFECTO | Deuda pendiente desde junio (ver bitácora) |
| `depth` (métrico vs normalizado) | — | ⚠️ SIN VERIFICAR | Deuda pendiente: ¿está en mm reales? |
| region_label (particiones) | — | ⚠️ SIN VERIFICAR | ¿Córtex/médula/seno bien etiquetados tras fix 014? |
| Sesgo polar de pirámides | pol ~21% vs central ~5% | 🟡 FIDELIDAD | Hipótesis viva para el vacío avascular central |

**Preguntas de auditoría para Code:**
1. ¿Los semiejes del elipsoide corresponden a un riñón de referencia citado?
2. ¿El umbral 6.6 mm sigue anclado a la fuente MDCT, y esa fuente es verificable?
3. Estado de las 4 deudas de Capa 0 (seno, depth, sesgo polar, radio raíz).

---

## Capa 1 — Nefronas

Archivo: `capa1_nefronas.py` → `capa1_nefronas.npz`
(glomerulos, tubulos, tipo, piramide_destino, depth_*)

| Parámetro | Valor actual | Estado | Ancla / Nota |
|---|---|---|---|
| N glomérulos | 1300 | 📐 SUPUESTO | Representativos, no 1:1 con ~1M reales. Declarado. |
| Reparto cortical/yuxta | 85% / 15% | ⚠️ SIN VERIFICAR | Code midió 85.0/15.0 exacto — reconfirmar cita del 85% |
| Poisson-disk d_min | 0.8 mm (real 0.807) | 📐 SUPUESTO | Separación de diseño; ¿tiene base o es solo espaciado? |
| CORTICAL_PENETRACION | — | 🔴 DEFECTO | Tier A "cita pendiente" + es parte del bug distal |
| JUXTA_PENETRACION | — | 🔴 DEFECTO | Tier A "cita pendiente" + parte del bug distal |
| **tubulos[:,5] (distal)** | sobre eje pirámide | 🔴 DEFECTO | **Bug confirmado L219: sin dispersión radial/azimutal → colineal por pirámide** |
| tubulos como recta | linspace, 6 pts | 📐 SUPUESTO | Declarado "versión simple" L54; cosmético p/colector, importa p/Capa 2 |
| depth cortical vs yuxta | 2.08 vs 5.85 mm | ⚠️ SIN VERIFICAR | Orden correcto; ¿valores anclados? |

**Defecto raíz confirmado (este diagnóstico):**
`T = B + p·(A−B)` coloca el extremo distal sobre el eje de la pirámide. Corrección
propuesta: dispersión radial/azimutal dentro del CONO de la pirámide.
→ **Requisito del régimen nuevo:** el ángulo del cono y la distribución radial
deben anclarse a literatura (geometría de la papila/médula), NO inventarse.

**Preguntas de auditoría para Code:**
1. ¿Qué dice la literatura sobre la geometría de convergencia tubular hacia la
   papila? ¿Qué ángulo/dispersión es anatómicamente defendible?
2. ¿CORTICAL/JUXTA_PENETRACION pueden anclarse a una fuente, o quedan como
   supuesto declarado?

---

## Capa 2 — Demanda

Archivo: `capa2_demanda.py` → `capa2_demanda.npz`

| Parámetro | Valor actual | Estado | Ancla / Nota |
|---|---|---|---|
| Fuente de puntos | 7800 vértices de tubulos | 🟡 FIDELIDAD | Rango 3 OK, pero hereda geometría "recta al eje" |
| Campo de demanda | realce por túbulo | ⚠️ SIN VERIFICAR | ¿Ponderación con base fisiológica? |
| Cobertura difusión 150µm | 0.78% | ✅ ANCLADO | Confirmado como límite de resolución, NO métrica de validación de Capa 3 |
| Umbral validación Capa 3 | geométrico (100% alcance) | ✅ ANCLADO | Decisión correcta ya tomada en bitácora |

**Nota:** Capa 2 se regenera tras corregir Capa 1 (hereda los tubulos).

---

## Capa 3a — Arterial

Archivo: `capa3a_arterial.py` → `capa3a_arterial.npz` (12016 nodos, 12015 aristas)

| Parámetro | Valor actual | Estado | Ancla / Nota |
|---|---|---|---|
| Atractores | glomerulos (sano) | ✅ ANCLADO | Campo sano; topología limpia |
| Ley de Murray | α=3.0 impuesto | 📐 SUPUESTO | Impuesto por construcción → test de consistencia, NO resultado |
| Radio de raíz | 0.214 mm | 🔴 DEFECTO | Orden de magnitud vs arteria renal real (2-3mm) — deuda desde junio |
| 100% glomérulos alcanzados | sí | ✅ ANCLADO | Métrica geométrica válida |
| Radios por-arista | — | ⚠️ SIN VERIFICAR | ¿Coherentes con Murray tras remapeo? |

---

## Capa 3b — Venoso

Archivo: `capa3b_venoso.py` → `capa3b_venoso.npz` (23799 nodos, 23798 aristas)

| Parámetro | Valor actual | Estado | Ancla / Nota |
|---|---|---|---|
| Shunts VV (1927) | 79% Y-confluencias | ✅ ANCLADO | Cerrado en bitácora: no es defecto estructural |
| Segmentos a 0.6mm exacto | artefacto paso fijo | 🟡 FIDELIDAD | Conocido; relevante para voxelización |
| Feasibility gate venoso | PASA | ✅ ANCLADO | No en límite duro de imprimibilidad |

---

## Capa 3ab — Peritubular

Archivo: `capa3ab_peritubular.py`

| Parámetro | Valor actual | Estado | Ancla / Nota |
|---|---|---|---|
| Muestreo de tubulos | — | 🟡 FIDELIDAD | Hereda geometría recta-al-eje de Capa 1 |
| 4290 puntos representativos | — | ⚠️ SIN VERIFICAR | ¿Sin aristas/radios? (nube, no árbol) |

---

## Capa 3c — Colector urinario

Archivo: `capa3c_colector.py` → `capa3c_colector.npz`

| Parámetro | Valor actual | Estado | Ancla / Nota |
|---|---|---|---|
| **Geometría del árbol** | 10 filamentos colineales | 🔴 DEFECTO | **Consecuencia del bug de Capa 1. 492 posiciones únicas de 2189** |
| Radios Horton-Strahler | terminal 15µm, papila 200µm | ✅ ANCLADO | Anclado a ducto de Bellini según bitácora |
| Forest de 10 subárboles | 1 por pirámide | 📐 SUPUESTO | Estructura correcta; falla la geometría espacial |
| Auditoría alcance/continuidad | PASA | ⚠️ SIN VERIFICAR | **Pasa sobre una recta — la auditoría no comprueba dimensionalidad** |

**Lección metodológica (para BITÁCORA):** la auditoría de capa3c midió alcance,
continuidad, radios y proximidad, y todo "pasó" sobre geometría colineal, porque
**nunca comprobó dimensionalidad (rango/SVD)**. Añadir chequeo de dimensionalidad
a la batería de auditoría de toda capa con geometría espacial.

---

## Capa 4 — Calicial (cálices + pelvis)

Archivo: `capa4_colector_alto.py` → `capa4_colector_alto.npz` (24 nodos, 23 aristas)

| Parámetro | Valor actual | Estado | Ancla / Nota |
|---|---|---|---|
| Geometría | hereda de Capa 3c | 🔴 DEFECTO | Se construye sobre papilas/nodos del colector roto |
| Infundíbulo | ~4 mm | ✅ ANCLADO | CTU n=1321 según bitácora — reconfirmar |
| Pelvis (elipsoide aplanado) | Z_semi=2.0mm FIX1 | ⚠️ SIN VERIFICAR | ¿Anclado o supuesto? |
| Ducto de Bellini | — | ✅ ANCLADO | Literatura, según bitácora |

---

## Orden de ejecución propuesto (cascada de dependencias)

La cascada obliga a corregir de arriba a abajo. Regenerar Capa 1 invalida 2→4.

```
1. Capa 0  — cerrar las 4 deudas (seno, depth, sesgo polar, radio raíz)
             [algunas afectan la cascada entera → hacer PRIMERO]
2. Capa 1  — corregir bug tubulos[:,5] con dispersión anclada a literatura
3. Capa 2  — regenerar (hereda tubulos)
4. Capa 3ab— regenerar (hereda tubulos)
5. Capa 3a — regenerar métricas (radio raíz si se corrige en paso 1)
6. Capa 3c — regenerar (ahora crece hacia atractores 3D sanos)
7. Capa 4  — regenerar (hereda colector sano)
8. Re-auditar dimensionalidad (SVD rango 3) en TODA capa espacial
```

**Backup obligatorio:** antes de sobrescribir cualquier `.npz`, copia fechada.

---

## Protocolo de verificación (regla del proyecto)

Para cada celda que pase a ✅ o 📐:
1. Code propone el ancla/valor y muestra la evidencia.
2. Carlos ejecuta la verificación desde su propia terminal.
3. Solo entonces se marca el estado y se registra la cita.
4. Ningún reporte de agente cuenta como verificación por sí solo.

---

## Registro de cambios de este plan

| Fecha | Cambio |
|---|---|
| (hoy) | Documento base creado. Estados iniciales desde diagnóstico del colector y bitácora. Todo lo no confirmado marcado ⚠️. |
