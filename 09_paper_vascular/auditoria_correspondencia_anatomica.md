# Auditoría de correspondencia anatómica — Inventario de parámetros dimensionales

**Programa:** Bio-Kidney AI 2026
**Alcance:** `08_gemelo_digital/` (capa0 → capa3c) + `09_paper_vascular/`. El generador
canónico de Capa 0 (`capa0_dominio.py`) vive en la raíz del repo, no en `08_gemelo_digital/`;
se incluye por ser la Capa 0.
**Naturaleza:** SOLO LECTURA. No se cambió ningún valor. Fecha inventario: 2026-07-05.
**Actualización 2026-07-06:** colocadas **2 citas verificadas** (Emamian 1993 → semiejes;
Beland 2010 → grosor cortical). El resto de filas [A] queda como `[CITA PENDIENTE: <dato>]`.
No se inventó ninguna otra referencia. Fila `UMBRAL_CM` actualizada a `GROSOR_CORTICAL_MM`
(el código cambió en la entrada 014 de bitácora: umbral fracción → mm absolutos, 6.6 mm).

## Propósito
Inventario de todo parámetro que **afirma una magnitud física o una proporción anatómica**,
para auditar su correspondencia con la literatura. Cada fila queda clasificada en:

- **[A] AFIRMA-ANATOMÍA** — pretende ser un valor anatómico/físico real (radios, longitudes,
  profundidades, diámetros, proporciones de población, exponentes de leyes físicas, ángulos).
  **Necesita fuente.**
- **[B] ESCALA-DECLARADA** — reducción representativa explícita, no pretende ser el valor real
  (conteos reducidos, pasos de crecimiento, umbrales de algoritmo, discretización). **No necesita
  fuente**: ya viene encuadrada.

Las columnas `valor de literatura` y `veredicto` se completan a medida que se **verifican**
citas. Las filas [A] sin cita verificada llevan `[CITA PENDIENTE: <dato>]`; las [B] llevan
`— (no requiere, [B])`.

---

## Tabla de inventario

| parámetro | archivo:línea | valor (unidad si declarada) | unidad declarada (sí/no) | categoría [A/B] | origen conocido (comentario) o SIN DOCUMENTAR | valor de literatura | veredicto |
|---|---|---|---|---|---|---|---|
| `A_SEMI` (semieje X, eje largo) | capa0_dominio.py:36 | 55.0 mm | sí (mm) | **A** | SIN DOCUMENTAR (comentario da eje, no fuente; paper repite 55×30×18) | **Emamian SA et al. AJR 1993;160(1):83-86** (DOI 10.2214/ajr.160.1.8416654); dimensiones renales por ecografía, n=665 adultos | 2×semieje = 110 mm (longitud del órgano) → dentro del rango de longitud renal adulta; afinar contra la media exacta de Emamian |
| `B_SEMI` (semieje Y, medial-lateral) | capa0_dominio.py:37 | 30.0 mm | sí (mm) | **A** | SIN DOCUMENTAR | **Emamian SA et al. AJR 1993;160(1):83-86** (DOI 10.2214/ajr.160.1.8416654) | 2×semieje = 60 mm (anchura del órgano) → dentro del rango adulto (Emamian); afinar con media |
| `C_SEMI` (semieje Z, antero-posterior) | capa0_dominio.py:38 | 18.0 mm | sí (mm) | **A** | SIN DOCUMENTAR | **Emamian SA et al. AJR 1993;160(1):83-86** (DOI 10.2214/ajr.160.1.8416654) | 2×semieje = 36 mm (espesor AP) → dentro del rango adulto (Emamian); afinar con media |
| `HILIO` (posición del hilio) | capa0_dominio.py:42 | [0, −30, 0] mm | sí (mm, vía B_SEMI) | **A** | Construcción: cara medial del elipsoide (= −B_SEMI); sin fuente independiente | `[CITA PENDIENTE: posición del hilio en la cara medial renal]` |  |
| `CENTRO_SENO` (centro elipsoide de exclusión) | capa0_dominio.py:48 | [0, −34, 0] mm | sí (mm) | **A** | SIN DOCUMENTAR (comentario: "ajustado" para no partir el órgano) | `[CITA PENDIENTE: geometría/posición del seno renal]` |  |
| `SEMIEJES_SENO` (dimensiones del seno renal) | capa0_dominio.py:49 | [22, 16, 11] mm | sí (mm) | **A** | SIN DOCUMENTAR (comentario: "ajustado") | `[CITA PENDIENTE: dimensiones del seno renal]` |  |
| `GROSOR_CORTICAL_MM` (unión cortico-medular; reemplazó a `UMBRAL_CM` frac. en entrada 014) | capa0_dominio.py:58 | 6.6 mm (absoluto) | sí (mm) | **A** | DOCUMENTADO con cita (ver →); el antiguo `UMBRAL_CM=0.30` (fracción, línea 63) queda **deprecado** | **Beland MD, Walle NL, Machan JT, Cronan JJ. AJR 2010;195(2):W146-149**; espesor cortical renal por ecografía: **media 5.9 mm, rango 3.2–11.0 mm** | **[OK]** 6.6 mm cae **dentro del rango 3.2–11.0 mm** (por encima de la media 5.9 mm). Nota: el comentario en código citaba "MDCT n=2068" (no verificado); la fuente **verificada** es Beland (ecografía) |
| `N_PIRAMIDES` (pirámides medulares) | capa0_dominio.py:55 | 10 | sí (conteo) | **A** | DOCUMENTADO en comentario: "rango humano fisiológico 8–18" (sin cita formal) | `[CITA PENDIENTE: nº de pirámides medulares (rango 8–18)]` |  |
| `CONE_HALF_ANGLE_DEG` (semiángulo cono piramidal) | capa0_dominio.py:56 | 22.0 ° | sí (grados) | **A** | SIN DOCUMENTAR | `[CITA PENDIENTE: geometría/ángulo de la pirámide medular]` |  |
| `N_POINTS` (puntos de parénquima muestreados) | capa0_dominio.py:58 | 200 000 | sí (conteo) | **B** | Escala computacional (resolución de muestreo, no anatómica) | — (no requiere, [B]) | — |
| `N_NEFRONAS` (nefronas sembradas) | capa1_nefronas.py:31 | 1300 | sí (conteo) | **B** | DECLARADA: "versión reducida representativa (real ~1e6)" | — (no requiere, [B]) | — |
| `FRAC_YUXTAMEDULAR` (proporción de población) | capa1_nefronas.py:32 | 0.15 (15% yuxta / 85% cortical) | sí (fracción) | **A** | SIN DOCUMENTAR (comentario afirma la proporción, sin fuente) | `[CITA PENDIENTE: proporción yuxtamedular/cortical ~15/85]` |  |
| `DIST_MIN` (separación mín. glomérulos, Poisson-disk) | capa1_nefronas.py:35 | 0.8 mm | sí (mm) | **B** | Escala: espaciado ligado a la densidad reducida (1300), no distancia inter-glomerular real | — (no requiere, [B]) | — |
| `CORTICAL_DEPTH_SCALE` (peso exp(−depth/scale)) | capa1_nefronas.py:38 | 0.12 | no (peso adimensional) | **B** | Modelado (sesgo de muestreo, no afirma anatomía) | — (no requiere, [B]) | — |
| `JUXTA_BAND_FRAC` (banda yuxtamedular = umbral×frac) | capa1_nefronas.py:42 | 0.75 | no (fracción de umbral) | **B** | Modelado (define banda de siembra; no valor medido) | — (no requiere, [B]) | — |
| `CORTICAL_PENETRACION` (penetración túbulo cortical) | capa1_nefronas.py:46 | 0.25 (fracción base→ápice) | sí (fracción de long. pirámide) | **A** | SIN DOCUMENTAR (comentario: "asa corta, médula externa", sin fuente) | `[CITA PENDIENTE: profundidad del asa de Henle cortical]` |  |
| `JUXTA_PENETRACION` (penetración túbulo yuxta) | capa1_nefronas.py:47 | 0.85 (fracción base→ápice) | sí (fracción de long. pirámide) | **A** | SIN DOCUMENTAR (comentario: "asa larga, hacia la papila", sin fuente) | `[CITA PENDIENTE: profundidad del asa de Henle yuxtamedular]` |  |
| `PENETRACION_JITTER` (variabilidad penetración) | capa1_nefronas.py:48 | 0.08 | sí (fracción) | **B** | Modelado (ruido aleatorio) | — (no requiere, [B]) | — |
| `TUBULE_NPTS` (vértices por túbulo) | capa1_nefronas.py:50 | 6 | sí (conteo) | **B** | Discretización: "versión simple; se refinará" | — (no requiere, [B]) | — |
| `UMBRAL_DIFUSION_UM` (umbral difusión O₂) | capa2_demanda.py:33 | 150 µm | sí (µm) — **FLAG: es DISTANCIA tejido→capilar, NO radio ni diámetro de vaso** | **A** | DOCUMENTADO en comentario: "límite hipoxia ~200 µm, diseño con margen" (sin cita formal) | `[CITA PENDIENTE: límite de difusión de O₂ ~100–200 µm]` |  |
| `UMBRAL_DIFUSION_MM` (mismo umbral en mm) | capa2_demanda.py:34 | 0.150 mm | sí (mm) — misma FLAG (distancia de cobertura) | **A** | DOCUMENTADO (igual que arriba; el paper lo cita como 100–150 µm) | `[CITA PENDIENTE: límite de difusión de O₂ ~100–200 µm]` |  |
| `DEMANDA_CORTEX` (demanda relativa córtex) | capa2_demanda.py:36 | 1.0 | no (relativo/normalización) | **B** | Referencia de normalización (no magnitud medida) | — (no requiere, [B]) | — |
| `DEMANDA_MEDULA` (demanda relativa médula) | capa2_demanda.py:37 | 0.4 (médula = 40% del córtex) | no (relativo) | **A** | SIN DOCUMENTAR (afirma proporción metabólica córtex/médula, sin fuente) | `[CITA PENDIENTE: proporción de consumo de O₂ médula/córtex]` |  |
| `RADIO_INFLUENCIA_NEFRONA_MM` (radio de realce) | capa2_demanda.py:39 | 0.5 mm | sí (mm) | **B** | Modelado (kernel de realce de demanda) | — (no requiere, [B]) | — |
| `REALCE_PESO` (peso del realce local) | capa2_demanda.py:40 | 0.3 | no (peso) | **B** | Modelado | — (no requiere, [B]) | — |
| `DIST_INFLUENCIA` (space colonization) | capa3a_arterial.py:32 | 8.0 mm | sí (mm) | **B** | Parámetro de algoritmo (Runions 2005), no anatómico | — (no requiere, [B]) | — |
| `DIST_MATAR` (atractor alcanzado) | capa3a_arterial.py:33 | 1.2 mm | sí (mm) | **B** | Parámetro de algoritmo | — (no requiere, [B]) | — |
| `PASO_CRECIMIENTO` (segmento por iteración) | capa3a_arterial.py:34 | 0.6 mm | sí (mm) | **B** | Escala de crecimiento (declarada) | — (no requiere, [B]) | — |
| `MURRAY_EXP` (exponente Ley de Murray, arterial) | capa3a_arterial.py:36 | 3.0 | sí (exponente adimensional) | **A** | DOCUMENTADO: "Ley de Murray r³=Σr_hijas³" (ley nombrada; sin cita formal) | `[CITA PENDIENTE: exponente de la Ley de Murray k=3 (Murray 1926)]` |  |
| `RADIO_TERMINAL` (arteriola aferente terminal) | capa3a_arterial.py:37 | 0.012 mm (~12 µm r) | sí (mm) | **A** | SIN DOCUMENTAR (comentario nombra "arteriola aferente terminal", sin medición/cita) | `[CITA PENDIENTE: radio/diámetro de la arteriola aferente]` |  |
| `GRADO_MAX_HIJOS` (bifurcación binaria, arterial) | capa3a_arterial.py:46 | 2 | sí (conteo) | **A** | DOCUMENTADO: "un vaso real se bifurca en 2" (justificado en comentario) | `[CITA PENDIENTE: topología binaria de bifurcación vascular]` |  |
| `radio_raiz` **(DERIVADO, no constante)** — radio raíz arterial | capa3a_arterial.py:372,394 | ~0.21 mm (r_in[0], sale de Murray) | sí (mm) | **A** | Comentario compara con "arteria renal real ~2–3 mm"; gate propio 0.5–6.0 mm → el valor derivado queda MUY por debajo (posible desajuste de escala) | `[CITA PENDIENTE: radio de la arteria renal ~2–3 mm]` |  |
| `PUNTOS_POR_TUBULO_CORTICAL` (densidad drenaje) | capa3ab_peritubular.py:53 | 3 | sí (conteo) | **B** | Escala de drenaje representativa | — (no requiere, [B]) | — |
| `PUNTOS_POR_TUBULO_YUXTA` (densidad drenaje) | capa3ab_peritubular.py:54 | 5 | sí (conteo) | **B** | Escala (asa larga + vasa recta) | — (no requiere, [B]) | — |
| puntos de drenaje totales **(DERIVADO)** | capa3ab_peritubular.py:53–54 | ~4290 (3×corticales + 5×yuxta) | sí (conteo) | **B** | Escala representativa vs. capilares reales (~10⁶); derivado de los dos anteriores | — (no requiere, [B]) | — |
| `JITTER_RADIAL_MM` (dispersión radial de puntos) | capa3ab_peritubular.py:55 | 0.15 mm | sí (mm) | **B** | Modelado (los puntos rodean el túbulo) | — (no requiere, [B]) | — |
| `DESPLAZAMIENTO_SATELITAL` (separación vena-arteria) | capa3b_venoso.py:59 | 1.5 mm | sí (mm) | **B** | Modelado (sesgo satelital; la relación vena-satélite es real, la magnitud es de diseño) | — (no requiere, [B]) | — |
| `PESO_SATELITAL` (peso del sesgo satelital) | capa3b_venoso.py:60 | 0.3 | no (peso) | **B** | Modelado | — (no requiere, [B]) | — |
| `MURRAY_EXP` (exponente Ley de Murray, venoso) | capa3b_venoso.py:65 | 3.0 | sí (exponente) | **A** | DOCUMENTADO: "Ley de Murray" (nombrada; sin cita formal) | `[CITA PENDIENTE: exponente de la Ley de Murray k=3 (Murray 1926)]` |  |
| `RADIO_TERMINAL_VENOSO` (vénula terminal) | capa3b_venoso.py:66 | 0.018 mm (~18 µm r) | sí (mm) | **A** | SIN DOCUMENTAR (comentario: "MAYOR que arterial 0.012", sin cita) | `[CITA PENDIENTE: radio/diámetro de la vénula terminal]` |  |
| `FACTOR_RADIO_VENOSO` (escala lumen venoso) | capa3b_venoso.py:67 | 1.3 (× Murray) | no (factor) | **A** | SIN DOCUMENTAR (afirma proporción vena/arteria por baja presión, sin fuente) | `[CITA PENDIENTE: proporción de calibre lumen vena/arteria]` |  |
| `RADIO_TERMINAL` (colector cortical) **[FILA YA VERIFICADA — formato]** | capa3c_colector.py:62 | 0.015 mm (15 µm r / 30 µm ⌀) | sí (mm; r y ⌀ explícitos) | **A** | DOCUMENTADO con cita informal en código: colector cortical 20–50 µm ⌀ (Britannica / morfometría renal) | `[CITA PENDIENTE: formalizar ref — colector cortical 20–50 µm ⌀]` |  |
| `RADIO_PAPILA` (conducto de Bellini) **[FILA YA VERIFICADA — formato]** | capa3c_colector.py:63 | 0.200 mm (200 µm r / 400 µm ⌀) | sí (mm; r y ⌀ explícitos) | **A** | DOCUMENTADO con cita informal en código: Bellini 300–600 µm ⌀ (grading nefrolitiasis, PMC) | `[CITA PENDIENTE: formalizar ref — conducto de Bellini 300–600 µm ⌀]` |  |
| `DIST_MATAR` (colector; más ajustado que arterial) | capa3c_colector.py:49 | 0.6 mm | sí (mm) | **B** | Parámetro de algoritmo (0.6 vs 1.2 arterial, para acercarse al distal) | — (no requiere, [B]) | — |
| `UMBRAL_PROX` (proximidad inter-subárbol) | capa3c_colector.py:68 | 0.150 mm | sí (mm) | **B** | Umbral de auditoría (análogo shunt VV), no anatómico | — (no requiere, [B]) | — |
| `PASO_DENSIF_MM` (paso de densificación de aristas) | capa3_cobertura_difusion_deficit.py:37 | 0.1 mm | sí (mm) | **B** | Discretización (< umbral 0.15 mm) | — (no requiere, [B]) | — |
| `R_CENTRO` (radios de sonda "centro geométrico") | capa3_cobertura_difusion_deficit.py:39 | [5, 10, 15] mm | sí (mm) | **B** | Radios de análisis (sonda), no anatómicos | — (no requiere, [B]) | — |

---

## Parámetros de modelado adimensionales (no afirman anatomía — fuera de la clase [A])
Incluidos arriba como [B] por completitud, no requieren fuente: `CORTICAL_DEPTH_SCALE`,
`JUXTA_BAND_FRAC`, `PENETRACION_JITTER`, `REALCE_PESO`, `PESO_SATELITAL`, `DEMANDA_CORTEX`.
Los `SEED = 2026` (todas las capas) no son dimensionales y se omiten.

## Notas de encuadre relevantes para la auditoría
- **Umbral de difusión O₂ (150 µm):** es una **distancia máxima tejido→capilar** (radio de
  cobertura de suministro), **NO** un radio ni un diámetro de vaso. Al contrastar con literatura,
  compárese contra la distancia de difusión de O₂ (~100–200 µm), no contra calibres vasculares.
- **`radio_raiz` arterial (~0.21 mm):** derivado de Murray desde `RADIO_TERMINAL=0.012 mm` y 1300
  terminales. El propio código lo compara con "arteria renal real ~2–3 mm" y su gate acepta
  0.5–6.0 mm; el valor derivado queda por debajo. Es un **desajuste de escala esperable** del
  modelo reducido (menos terminales ⇒ raíz más fina), a resolver en la columna `veredicto`.
- **Radios de Capa 3c** (15 µm / 200 µm) son la referencia de formato: ambas filas [A] con cita
  **informal en el código** (colector cortical 20–50 µm ⌀; Bellini 300–600 µm ⌀), que aún debe
  **formalizarse** con referencia bibliográfica → marcadas `[CITA PENDIENTE]`. El resto de radios
  ([A]) del árbol arterial y venoso tampoco tienen ese anclaje todavía.
- **`GROSOR_CORTICAL_MM = 6.6 mm`** (antes `UMBRAL_CM = 0.30` fracción): tras la entrada 014 el
  umbral cortico-medular es **absoluto en mm** e invariante al normalizador. Cita **verificada**
  (Beland 2010): 6.6 mm dentro del rango 3.2–11.0 mm. **Discrepancia a resolver:** el comentario
  en `capa0_dominio.py:58` cita "MDCT n=2068", que **no** es la fuente verificada (Beland es
  ecografía) — conviene corregir el comentario del código a la referencia real.

---

## Resultado: lista de trabajo real (clase [A] SIN cita verificada)

De los **14 parámetros [A]** originalmente SIN DOCUMENTAR, **4 quedan ahora con cita verificada**
(3 semiejes → Emamian 1993; grosor cortical → Beland 2010). **Quedan 10 pendientes** de fuente:

1. `CENTRO_SENO = [0,−34,0] mm` (capa0_dominio.py:48) — posición del seno renal
2. `SEMIEJES_SENO = [22,16,11] mm` (capa0_dominio.py:49) — dimensiones del seno renal
3. `CONE_HALF_ANGLE_DEG = 22.0°` (capa0_dominio.py:56) — apertura de la pirámide medular
4. `FRAC_YUXTAMEDULAR = 0.15` (capa1_nefronas.py:32) — proporción yuxta/cortical
5. `CORTICAL_PENETRACION = 0.25` (capa1_nefronas.py:46) — profundidad de asa cortical
6. `JUXTA_PENETRACION = 0.85` (capa1_nefronas.py:47) — profundidad de asa yuxtamedular
7. `RADIO_TERMINAL = 0.012 mm` (capa3a_arterial.py:37) — radio arteriola aferente
8. `DEMANDA_MEDULA = 0.4` (capa2_demanda.py:37) — proporción demanda O₂ médula/córtex
9. `RADIO_TERMINAL_VENOSO = 0.018 mm` (capa3b_venoso.py:66) — radio vénula terminal
10. `FACTOR_RADIO_VENOSO = 1.3` (capa3b_venoso.py:67) — proporción lumen vena/arteria

**Citas verificadas colocadas (2026-07-06):**
- `A_SEMI`, `B_SEMI`, `C_SEMI` (semiejes / dimensiones renales) → **Emamian SA, Nielsen MB,
  Pedersen JF, Ytte L. AJR 1993;160(1):83-86** (DOI 10.2214/ajr.160.1.8416654).
- `GROSOR_CORTICAL_MM` (6.6 mm, grosor cortical) → **Beland MD, Walle NL, Machan JT, Cronan JJ.
  AJR 2010;195(2):W146-149** (media 5.9 mm, rango 3.2–11.0 mm; 6.6 mm dentro del rango).

**Clase [A] con justificación en código pero cita aún pendiente de formalizar (no en la lista
de 10, pero conviene formalizar):** `N_PIRAMIDES=10` (rango 8–18), `MURRAY_EXP=3.0` (Ley de
Murray, ×2 archivos), `UMBRAL_DIFUSION=150 µm` (límite hipoxia ~200 µm), `GRADO_MAX_HIJOS=2`
(bifurcación binaria), `HILIO` (construcción medial), `radio_raiz` (comparado con arteria renal
2–3 mm), y las dos filas de formato de Capa 3c (`RADIO_TERMINAL`/`RADIO_PAPILA`, citas informales
Britannica/PMC por formalizar).
