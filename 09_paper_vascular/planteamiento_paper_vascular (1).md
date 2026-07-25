# Arquitectura vascular renal *in silico*: un marco computacional multicapa guiado por restricciones físicas y demanda metabólica

**Título (EN):** *In silico renal vascular architecture: a multilayer computational framework driven by physical constraints and metabolic demand*

**Programa:** Bio-Kidney AI 2026
**Autor:** Carlos David Moreno Cáceres · VirtusSapiens
**ORCID:** 0009-0005-3933-5072
**Repositorio:** github.com/VirtusSapiens/Bio-Kidney-AI-2026
**Estado del documento:** PLANTEAMIENTO (esqueleto) — se completa a medida que avanzan las capas.

> **Nota de alcance.** Este trabajo aborda la mitad *geométrica/de diseño* del cuello de botella de vascularización en ingeniería de tejidos renales — el diseño de una arquitectura vascular que satisfaga la restricción de difusión de oxígeno. NO aborda la mitad biológica (maduración celular, anastomosis con el huésped, perfusión sin necrosis), que requiere trabajo experimental. El modelo es una geometría idealizada paramétrica basada en morfometría poblacional, no una reconstrucción de un riñón real por imágenes médicas.

---

## Resumen (Abstract) — [POR COMPLETAR al cerrar Capa 3]

*Estructura objetivo (~250 palabras):*
- **Contexto:** la vascularización es el cuello de botella #1 en bioimpresión de órganos; el riñón es el órgano sólido más complejo.
- **Brecha:** la mayoría de los enfoques tratan la generación vascular de forma aislada o priorizan precisión geométrica sobre acoplamiento funcional; faltan marcos abiertos, integrados y reproducibles en hardware modesto.
- **Aporte:** un marco multicapa que integra dominio anatómico, unidades funcionales (nefronas), campo de demanda metabólica y generación vascular bajo restricciones físicas (Murray, límite de difusión de O₂).
- **Métodos:** [resumen de las 5 capas].
- **Resultados:** [métricas — cobertura de demanda, % de tejido dentro del umbral de difusión, conformidad con Murray, etc.].
- **Significancia:** marco abierto y reproducible que sirve de plano para etapas experimentales posteriores.

---

## 1. Introducción

- El problema del trasplante renal y la escasez de órganos (contexto motivacional, sin sobreventa).
- La vascularización como cuello de botella central en ingeniería de tejidos (con referencias 2025).
- Las dos mitades del problema: geométrica/diseño vs. biológica. Delimitación honesta del alcance de este trabajo a la primera.
- Estado del arte: CCO/GCO aplicado a riñón, generación a escala de órgano, marcos multiescala recientes. **Posicionamiento: el aporte es la INTEGRACIÓN abierta y reproducible, no la invención de la generación vascular.**
- Objetivo y contribuciones del trabajo.

## 2. Métodos — El marco multicapa

> Principio rector: orden de dependencia **forma → nefronas → demanda → árboles → colección**. El árbol vascular es la *solución* a un problema definido por las capas previas; la representación nativa es un GRAFO (no vóxeles).

### 2.1 Capa 0 — Dominio renal [HECHO]
- Elipsoide paramétrico (semiejes 55×30×18 mm), seno renal esculpido por sustracción (forma de frijol).
- Partición corteza/médula (~60/40), pirámides medulares, hilio.
- Campo de profundidad. *(Pendiente: migrar a distancia métrica real en mm.)*

### 2.2 Capa 1 — Nefronas [HECHO]
- Siembra de 1300 glomérulos en corteza (modelo reducido representativo; real ~10⁶).
- Dos poblaciones: corticales (85%) y yuxtamedulares (15%).
- Túbulos como polilíneas descendentes hacia pirámides.
- Validación: corte sagital confirma estratificación cortical.

### 2.3 Capa 2 — Campo de demanda metabólica [EN CURSO]
- Campo escalar de demanda de O₂ sobre el dominio (corteza > médula).
- Restricción de diseño: ningún punto a más de **100–150 µm** de un capilar (margen de seguridad bajo el límite de hipoxia de ~200 µm).
- Métrica de distancia: KD-tree sobre puntos (coordenadas reales en mm) para iteración barata.

### 2.4 Capa 3 — Árboles vasculares [PENDIENTE]
- Arterial: generación guiada por demanda (space colonization + CCO), Ley de Murray (r³=r₁³+r₂³) en cada bifurcación, ángulo óptimo.
- Venoso: espejo a baja presión.
- Colector urinario: árbol convergente (NO usa Murray ni CCO), drena glomérulo→pirámide→cáliz.
- Acoplamiento: cada glomérulo como nodo de transferencia.

### 2.5 Capa 4 — Cálices, pelvis y voxelización [PENDIENTE]
- Convergencia del colector hacia el seno.
- Voxelización adaptativa (octree) para imprimibilidad y campos de difusión.

## 3. Resultados [POR COMPLETAR]

- Métricas de cobertura de demanda (% de tejido dentro del umbral de difusión).
- Conformidad con la Ley de Murray (distribución de error por bifurcación).
- Estadísticas morfométricas del árbol (n.º de segmentos, órdenes, radios terminales) y comparación con rangos de la literatura.
- [Eventual] Recuperación de función fisiológica (TFG) emergente de la geometría.

## 4. Discusión [POR COMPLETAR]

- Qué resuelve el marco y qué no (recordar las dos mitades).
- Comparación honesta con el estado del arte.
- Limitaciones: geometría idealizada, modelo reducido de nefronas, validación funcional pendiente, falta de la dimensión biológica.
- Trabajo futuro: aumento de densidad, parametrización paciente-específica, acoplamiento con modelos de maduración.

## 5. Disponibilidad de datos y código

- Código abierto en el repositorio; datos generados reproducibles vía semilla (SEED=2026).
- Entorno: Ubuntu 24.04, Python 3.12 (NumPy/SciPy), hardware de consumo (i5, 16 GB).

---

## Hoja de ruta de redacción

| Sección | Depende de | Estado |
|---|---|---|
| Métodos 2.1–2.2 | Capas 0–1 | Listo para redactar |
| Métodos 2.3 | Capa 2 | En curso |
| Métodos 2.4–2.5 | Capas 3–4 | Pendiente |
| Resultados | Capa 3 cerrada | Pendiente |
| Abstract + Discusión | Todo lo anterior | Pendiente |

**Decisión de publicación:** redactar el preprint completo cuando se cierre la Capa 3 (árbol vascular funcional), que es el resultado que justifica la publicación. La bitácora (`BITACORA.md`) es el borrador natural de las secciones de métodos.
