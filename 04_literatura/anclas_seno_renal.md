# Anclas de literatura — seno renal

**Programa:** Bio-Kidney AI 2026 · **Parámetros afectados:** `SEMIEJES_SENO`, `CENTRO_SENO`
(`capa0_dominio.py:48-49`) · **Última actualización:** 2026-08-11

Registro de las fuentes disponibles hoy para acotar el seno renal del gemelo, qué mide cada una,
qué permite y qué **no** permite. Documento de trazabilidad: si una cifra del seno aparece en el
preprint, debe poder rastrearse hasta aquí.

**Estado del parámetro a día de hoy: SUPUESTO DECLARADO, ACOTADO POR ABAJO Y FALSADO.**
Ninguna de las fuentes de abajo fija un valor. Las dos permiten **falsar**, no **anclar**.

---

## 1. Caglar 2014 — grasa del seno renal por estereología sobre TC

**Referencia completa**

> Caglar V, Kucuk A, Aktas S, et al. *"Volumetric evaluation of fat in the renal sinus in normal
> subjects using stereological method on computed tomography images and its relationship with body
> composition."* **Folia Morphologica. 2014;73(3):302-308.**
> DOI **10.5603/FM.2014.0016** · PMID **25242158**

**Qué mide**

Volumen **de la grasa** contenida en el seno renal, en sujetos normales, por **método
estereológico** aplicado a imágenes de **tomografía computarizada**. n = **240** sujetos,
**21–80 años**. El estudio relaciona además ese volumen con la composición corporal.

**Valores**

| grupo | grasa sinusal (cm³) |
|---|---|
| hombres, riñón **izquierdo** | **5.70 ± 2.87** |
| hombres, riñón **derecho** | **4.15 ± 2.39** |
| mujeres, riñón **izquierdo** | **3.51 ± 2.67** |
| mujeres, riñón **derecho** | **2.49 ± 2.16** |

**Qué ancla — una COTA INFERIOR, no un valor**

El seno renal contiene grasa **más** el sistema colector (pelvis, cálices mayores y menores), la
arteria y la vena renales con sus ramas de primer orden, vasos linfáticos y tejido conectivo laxo.
Por tanto, para un mismo riñón:

> **V_seno > V_grasa_sinusal**, necesariamente y por definición de contención.

Eso convierte los valores de Caglar en una **cota inferior falsable** del volumen del seno. Es todo
lo que aportan, y es suficiente para detectar un seno demasiado pequeño.

**Aplicación al gemelo (2026-08-11)**

El volumen **excavado** por el seno del gemelo — la intersección del elipsoide de exclusión con el
parénquima — es **3.9599 mL**:

```bash
.venv/bin/python -c "import numpy as np; d=np.load('capa0_dominio.npz'); v=(float(d['vol_elipsoide_mm3'])-float(d['vol_parenquima_mm3']))/1000; print('excavado %.4f mL'%v); [print('  vs %-6s %.2f -> %.3f %s'%(k,g,v/g,'VIOLA' if v<g else 'ok')) for k,g in [('H izq',5.70),('H der',4.15),('M izq',3.51),('M der',2.49)]]"
```

| comparación | cociente | cumple V_seno > V_grasa |
|---|---|---|
| 3.9599 vs **5.70** (H izq) | 0.700 | **NO — violada por 1.740 mL** |
| 3.9599 vs **4.15** (H der) | 0.962 | **NO — violada por 0.190 mL** |
| 3.9599 vs 3.51 (M izq) | 1.128 | sí, margen 0.450 mL |
| 3.9599 vs 2.49 (M der) | 1.590 | sí |

**La cota está violada en los dos grupos masculinos.** Registro completo en `00_bitacora/BITACORA.md`,
ENTRADA 033 §2-bis.

**Limitaciones para este uso**

1. **Mide grasa, no seno.** Es un subconjunto propio del contenido sinusal. **Nunca escribir "el seno
   debería medir 5.70 mL".** El único uso legítimo es como cota inferior.
2. **No da ejes ni posición.** No permite fijar `SEMIEJES_SENO` ni `CENTRO_SENO`; sólo detectar que el
   volumen resultante es demasiado pequeño.
3. **El gemelo no declara lateralidad ni sexo.** **[PENDIENTE DE ANCLA]** cuál de los cuatro grupos le
   corresponde. Por eso la comparación se reporta contra los cuatro y la conclusión se limita a "en
   los dos grupos masculinos".
4. **Las desviaciones típicas son grandes** (hasta ±2.87 sobre una media de 5.70, un 50 %). La grasa
   sinusal correlaciona con composición corporal, que el gemelo no modela.
5. **Estereología sobre TC** define el borde graso por umbral de densidad; no es la misma frontera
   geométrica que el elipsoide de exclusión del gemelo.

---

## 2. Zhang 2023 — fracción grasa del seno por RM

**Referencia completa**

> Zhang QH, et al. **Frontiers in Endocrinology. 2023;14:1187781.**
> DOI **10.3389/fendo.2023.1187781**

**Qué mide**

**Fracción grasa** (proporción, no volumen) del seno renal, por **mapeo de fracción grasa en
resonancia magnética**. n = **126**.

**Valores**

| grupo | fracción grasa del seno (%) |
|---|---|
| hombres, riñón **derecho** | **28.33 ± 6.73** |
| hombres, riñón **izquierdo** | **31.21 ± 6.29** |
| mujeres, riñón **derecho** | **23.82 ± 7.74** |
| mujeres, riñón **izquierdo** | **27.92 ± 8.15** |

**Qué ancla**

La **proporción** de grasa dentro del seno: en torno al 24–31 % según grupo. Es decir, la grasa es
**menos de un tercio** del contenido sinusal; el resto es colector, vasos y conectivo.

Eso refuerza cualitativamente el argumento de contención de Caglar: si la grasa es ~30 % del seno, la
cota inferior real es bastante más alta que el valor de Caglar. **Pero ese refuerzo es cualitativo y
no debe cuantificarse** (ver §3).

**Limitaciones para este uso**

1. **No da ningún volumen absoluto.** Sola, no acota nada dimensionalmente.
2. **No es comparable en absoluto con Caglar.** Metodología distinta (RM vs TC), población distinta
   (n = 126 vs 240), criterio de segmentación distinto. Los volúmenes absolutos derivables de ambos
   estudios difieren **en casi un orden de magnitud**.
3. **No da ejes ni posición.** Igual que Caglar, no permite fijar los parámetros.

---

## 2-bis. Glodny 2009 — ancho parenquimatoso (PW), ancla INDIRECTA del seno

**Referencia completa**

> Glodny B, Unterholzner V, Taferner B, Hofmann KJ, Rehder P, Strasak A, Petersen J.
> *"Normal kidney size and its influencing factors — a 64-slice MDCT study of 1.040 asymptomatic
> patients."* **BMC Urology. 2009;9:19.** DOI **10.1186/1471-2490-9-19**.
> Texto completo: `https://pmc.ncbi.nlm.nih.gov/articles/PMC2813848/`

**Por qué es un ancla del seno pese a no medirlo**

Glodny no dimensiona el seno como cavidad. Pero el **ancho parenquimatoso** (PW) es la distancia
cápsula → seno: acota **dónde empieza el seno**, y por tanto restringe su posición y tamaño de forma
indirecta. En la parametrización del gemelo esa restricción es exacta y cerrada:

> `PW_max = |−B_SEMI − (CENTRO_SENO[1] + SEMIEJES_SENO[1])|`

es decir, **PW restringe la SUMA `cy + sy`**, no ninguno de los dos por separado. Contrastada contra
malla en 5 configuraciones, diferencia ≤ 1.8 × 10⁻⁴ mm (`_auditoria/2026-08-12`, Reserva 2a).

**Valores**

| magnitud | derecho | izquierdo |
|---|---|---|
| **PW** (parenchymal width) | **15.4 ± 2.8 mm** | **15.9 ± 2.7 mm** |
| **CW** (cortical width) | 6.6 ± 1.9 mm | 6.6 ± 1.9 mm |
| LPP (longitud polo-polo) | 108.5 ± 12.2 mm | 111.3 ± 12.6 mm |

n = 2068 riñones / 1040 adultos asintomáticos, MDCT 64 cortes.

### Qué ESPECIFICA el paper sobre el protocolo de medida del PW

- Fase de adquisición: `width of the parenchyma (PW) and the cortex (CW) in the arterial phase`.
- **Plano: axial.** Pie de la Figura 1: `Axial 0.625 mm collimated slice of the kidney in an arterial
  phase, with the strongly contrasted kidney cortex`, mostrando `Cortical width (CW), and parenchymal
  width (PW)`.
- Control de reproducibilidad: `The measurements were performed twice in a random sample of 50 data
  sets`.
- Para la **longitud** (no para PW) sí detalla: `axes were adjusted individually in double oblique
  planes using 3D software`.

### Qué NO especifica — [PENDIENTE DE ANCLA]

| pregunta | estado |
|---|---|
| **Nivel anatómico** al que se mide el PW (polo superior / tercio medio / hilio / polo inferior) | **[PENDIENTE DE ANCLA]** — el texto de Métodos no lo indica |
| **Número de medidas por riñón** | **[PENDIENTE DE ANCLA]** — no consta. El "twice in a random sample of 50" es control inter-observador, no el protocolo de rutina |
| Si es **localización estandarizada** o **promedio** de varias | **[PENDIENTE DE ANCLA]** |
| Si se mide en el punto de **máximo espesor** o en uno anatómicamente definido | **[PENDIENTE DE ANCLA]** |

**Advertencia de uso — es el punto crítico de esta ficha.** Sin saber a qué nivel y con qué
agregación mide Glodny, **no se sabe qué estadístico del gemelo es comparable con su media
poblacional**. Comparar un **máximo espacial sobre una geometría** contra una **media ± s.d. sobre
2068 riñones** es un error de magnitud de la misma clase que el corregido por el método B
(ENTRADA 032 §2): las unidades coinciden (mm) y la comparación parece legítima, y no lo es.

**Lo que sí se puede afirmar hoy:** que el parénquima del gemelo es más delgado que la referencia,
**en dirección**, porque **todos** los estadísticos candidatos caen bajo el límite inferior de Glodny
a −1 s.d. (12.6 mm): máximo 12.000 · p95 11.528 · mediana 6.855 · media 6.590.
**Lo que NO se puede afirmar:** la **magnitud** del déficit, que va de 3.4 mm (máximo) a 8.5 mm
(mediana) según el estadístico. Distribución completa en `_auditoria/2026-08-12`, Reserva 1d.

**Nota de trazabilidad:** «PW del gemelo» **no existe como código en el repo**
(`grep -rn "parenchymal\|PW\b" --include=*.py` → vacío). Es una construcción de auditoría introducida
el 2026-08-09. **No confundir con `depth_cortical_mm`**, que sí es un campo del `.npz` y mide otra
cosa (cápsula → punto, en todo el parénquima, incluidos los polos donde no hay seno). Definición
operativa exacta en `_auditoria/2026-08-12`, Reserva 1a.

---

## 3. NO ANCLADO — la hipótesis del cruce, y por qué no puede usarse

**La hipótesis.** Cruzar la fracción grasa de Zhang (~28–31 % en hombres) con el volumen de grasa de
Caglar (5.70 cm³ en hombre izquierdo) sugeriría un seno total del orden de **~14–15 cm³**, valor
cercano al elipsoide de exclusión del gemelo (**16.22 mL**). Eso apuntaría a que **el tamaño del
elipsoide es aproximadamente correcto y el error está en `CENTRO_SENO`**, no en `SEMIEJES_SENO`.

**Esa convergencia es interesante porque coincide con un diagnóstico independiente.** Los barridos de
`08_gemelo_digital/_auditoria/2026-08-09_seno_apices_particion.md` (Tarea 1c) y
`_auditoria/2026-08-10_centro_seno_particion.md` (Tarea 1a) muestran, **por vía puramente geométrica y
sin usar ninguna de estas dos fuentes**, que:

- `SEMIEJES_SENO` **no puede** corregir el defecto de los cuatro ápices corticales (ni escalando ×1.3
  baja de 2 ápices en territorio córtex);
- `CENTRO_SENO[1]` **sí** lo corrige, y es la variable dominante.

Dos vías independientes apuntan al mismo parámetro. **Eso es lo que se registra: una convergencia
cualitativa.**

> **Salvedad (2026-08-12) — qué significa aquí «independientes».**
> Se refiere al **tipo de evidencia**: una vía es literatura (cruce Caglar × Zhang), la otra es
> geometría (barridos de sensibilidad sobre el `.npz`). En ese sentido sí son independientes, y la
> convergencia es informativa.
>
> **No se refiere a las condiciones geométricas entre sí.** Las tres que se usan para acotar el seno
> —ápices corticales, PW y volumen excavado— son **estrictamente monótonas crecientes** en
> `CENTRO_SENO[1]`, porque las tres miden, con métricas distintas, cuánto se adentra el seno en el
> parénquima. **No pueden discrepar en dirección**, de modo que su acuerdo **no es corroboración
> mutua**: coinciden en la **dirección**, no convergen independientemente en un **valor**.
>
> Consecuencia práctica: la condición de ápices corticales resulta **redundante** (no fija ningún
> borde de la región admisible); los bordes los fijan Caglar por abajo y Glodny por arriba.
> Demostración en `08_gemelo_digital/_auditoria/2026-08-12_reservas_pw_region2d.md`, Reserva 3.

### Por qué el número ~14–15 cm³ NO puede usarse

> **PROHIBIDO usarlo como parámetro, como valor de referencia, o como justificación de `SEMIEJES_SENO`.**

1. **Es un cociente entre dos estudios no comparables.** Caglar (TC, estereología, n = 240) y Zhang
   (RM, mapeo de fracción grasa, n = 126) miden sobre poblaciones, técnicas y criterios de
   segmentación distintos. Sus volúmenes absolutos difieren en casi un orden de magnitud.
2. **Un número derivado cruzando dos fuentes no comparables no es un anclaje.** Presentarlo como tal
   sería exactamente el error de atribución corregido en la ENTRADA 032 §5 (el 6.6 mm atribuido a
   Beland cuando procedía de Glodny): confundir *validación post hoc* con *procedencia*.
3. **La coincidencia con 16.22 mL puede ser casual.** Y además compara con la magnitud equivocada: el
   elipsoide completo **no es un volumen anatómico** — el 76 % de él cae fuera del parénquima y no
   excava nada. La magnitud comparable es el **volumen excavado, 3.9599 mL**.

**Uso legítimo:** citarlo como *"dos vías independientes convergen en señalar `CENTRO_SENO` y no
`SEMIEJES_SENO`"*. **Uso ilegítimo:** cualquier frase que contenga «~14 cm³» o «~15 cm³» como cifra
del seno.

---

## 4. Lo que sigue PENDIENTE DE ANCLA

**Emamian 1993 — dimensiones del área ecogénica central.**

> Emamian SA, Nielsen MB, Pedersen JF, Ytte L. AJR 1993;160(1):83-86. DOI 10.2214/ajr.160.1.8416654

El repo cita esta fuente **sólo** para los tres semiejes del órgano (55/30/18 mm):
`09_paper_vascular/auditoria_correspondencia_anatomica.md:41-43` y `:140-141`. Se ha indicado que el
mismo estudio midió longitud, anchura, grosor y volumen del **área ecogénica central** (el seno), pero
**esos valores no aparecen en ninguna fuente secundaria accesible**. **[PENDIENTE DE ANCLA].**

Ese dato es el que permitiría **fijar** `SEMIEJES_SENO` y `CENTRO_SENO` en lugar de sólo acotarlos.

### ⚠ Aviso de integridad — cifras fabricadas en circulación

Circuló un documento externo que presenta dimensiones **reconstruidas** de la CEA de Emamian 1993
**como si fueran datos del paper**:

> **L = 4.5 cm · W = 1.3 cm · T = 1.2 cm · V = 3.68 cm³ — ESTAS CIFRAS NO ESTÁN EN LA FUENTE.**

**No usarlas, no citarlas, no derivar nada de ellas.** Búsqueda en el repo (2026-08-11): **no
aparecen**. Los únicos aciertos de la cadena `3.68` son coordenadas del árbol vascular sin relación
con el seno (`02_vascular_cco/renal_data_v1.json:2409,3377,4276,4293,5098`). La cadena `CEA` no
aparece en ningún archivo.

Si esas cifras se encuentran en algún material del proyecto, hay que retirarlas y anotar su
procedencia, no integrarlas.

---

## 5. Resumen operativo

| fuente | qué mide | qué permite | qué NO permite |
|---|---|---|---|
| **Caglar 2014** | volumen de **grasa** sinusal (TC, estereología, n=240) | **cota inferior** falsable de V_seno | fijar volumen, ejes o posición |
| **Zhang 2023** | **fracción** grasa del seno (RM, n=126) | saber que la grasa es ~24–31 % del seno | cualquier volumen absoluto |
| **Glodny 2009 (PW)** | cápsula → seno (MDCT, n=2068) | restringir la **suma** `cy + sy` | fijar `cy` o `sy` por separado; y su **nivel de medida** es [PENDIENTE DE ANCLA] |
| **cruce Caglar × Zhang** | — | registrar convergencia cualitativa hacia `CENTRO_SENO` | **cualquier cifra** («~14–15 cm³») |
| **Emamian 1993 (CEA)** | dimensiones y volumen del seno | **fijaría** los parámetros | **[PENDIENTE DE ANCLA]** — no disponible |

**Consecuencia para el gemelo:** el volumen excavado (3.9599 mL) **viola la cota inferior de Caglar**
en los dos grupos masculinos. Eso autoriza la afirmación *"el seno del gemelo está subdimensionado"*
con las precisiones de `00_bitacora/BITACORA.md`, ENTRADA 033 §6, y hace que la corrección de
`CENTRO_SENO` deje de ser opcional.

### Lo que las anclas restringen, y lo que NO

> **AVISO (2026-08-12): este apartado sustituye a un intervalo numérico retirado.**
> Una versión anterior de esta sección publicaba `CENTRO_SENO[1] ∈ [−31.37, −27.80]` como
> "intervalo admisible bajo las tres condiciones ancladas". **Esa cifra estaba mal planteada y queda
> retirada:** era un corte en `SEMIEJES_SENO[1] = 16`, valor que está tan sin anclar como `cy`.
> Auditoría en `08_gemelo_digital/_auditoria/2026-08-12_reservas_pw_region2d.md`.

| ancla | qué restringe realmente |
|---|---|
| **Glodny (PW ±1 s.d.)** | la **suma** `cy + sy ∈ [−17.4, −11.8]` — banda diagonal, no un intervalo en `cy` |
| **Caglar (V excavado > 5.70 mL)** | región curva en `(cy, sy, sx, sz)`; **cota inferior**, fija el borde inferior |
| **ápices corticales** | **nada** — redundante, se satisface automáticamente bajo Caglar |

**El espacio completo es 4D en `(cy, sy, sx, sz)`**: `SEMIEJES_SENO` entero está sin anclar, no sólo
`sy`. Cualquier región 2D que se reporte es una **rebanada** con `sx` y `sz` congelados.

Región admisible con `sx = 22`, `sz = 11` fijos — **rebanada, no resultado**:

| `sy` | `cy` admisible | anchura |
|---|---|---|
| 12 | [−28.50, −24.10] | 4.40 mm |
| 16 | [−31.30, −27.90] | 3.40 mm |
| 22 | [−35.70, −33.90] | 1.80 mm |

Tabla completa y comandos en `_auditoria/2026-08-12`, Reserva 2c.
**Ningún valor está elegido ni aplicado, y ninguna cifra de esta sección debe publicarse** hasta
determinar el nivel de medida del PW de Glodny (ver §2-bis, [PENDIENTE DE ANCLA]).

**Qué reduciría la región a un punto:** **las tres dimensiones de la CEA** de Emamian (`sx`, `sy`,
`sz`), tras lo cual Glodny fijaría `cy` por la suma. **No basta con su volumen**: es una ecuación
sobre cuatro incógnitas.
