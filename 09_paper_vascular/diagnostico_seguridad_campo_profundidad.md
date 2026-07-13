# Diagnóstico de seguridad — separar profundidad-cortical del rol del seno renal

**Programa:** Bio-Kidney AI 2026 · Fecha: 2026-07-05
**Naturaleza:** SOLO LECTURA. No se modificó ningún archivo, `.npz`, parámetro ni lógica. La
simulación del fix es **en memoria**; nada se guardó. El objetivo es decidir si es seguro que la
profundidad-cortical dependa SOLO de la cápsula externa, separándola del rol geométrico del seno.

## Contexto verificado en el código (leído, no asumido)
- `nearest_surface_distance()` (`capa0_dominio.py:98-111`) devuelve `np.minimum(dist_main, dist_seno)`
  — el **mínimo** entre distancia a la cápsula externa y distancia a la pared del seno. **Confirmado.**
- `compute_depth()` (`:114-128`) normaliza ese mínimo: `depth = raw / raw.max()`, con
  `depth_norm = raw.max()`. **Confirmado.**
- `cortex_mask = d < UMBRAL_CM` (`:271`) parte cortex/medula. **Confirmado.** `UMBRAL_CM=0.30` es
  una **fracción del `depth_norm`**, no un valor absoluto en mm — dato clave para la Pregunta 3.

---

## PREGUNTA 1 — Inventario de consumidores del campo combinado

| símbolo | archivo:línea | qué hace | clase |
|---|---|---|---|
| `nearest_surface_distance` | capa0_dominio.py:98-111 | `min(dist_capsula, dist_seno)`; único llamador es `compute_depth` | **[DEPENDE-AMBAS]** (es el mínimo, por definición) |
| `compute_depth` / `depth` | capa0_dominio.py:114-128, :268 | normaliza el mínimo a [0,1] | **[DEPENDE-AMBAS]** (hereda el min) |
| `cortex_mask = d < UMBRAL_CM` | capa0_dominio.py:271 | define cortex vs medula | **[DEPENDE-CÁPSULA]** en intención (el córtex es relativo a la cápsula), hoy contaminado por el seno |
| `medulla_mask` → `assign_pyramids` | capa0_dominio.py:272,276,218-244 | restringe la asignación piramidal a puntos de médula; el test de cono es **geometría pura** | [DEPENDE-CÁPSULA] (solo necesita "¿es médula?") |
| `depth_norm` (guardado) | capa0_dominio.py:304 | escalar informativo | **[NEUTRO]** — **ningún** módulo aguas abajo lo lee (verificado por grep) |
| `region_label` (guardado) | capa0_dominio.py:279-284 | cortex/medulla/piramide_XX | consumido aguas abajo (ver filas siguientes) |
| `cortex_idx = region=="cortex"` | capa1_nefronas.py:156 | **pool de siembra de TODOS los glomérulos** | **[DEPENDE-CÁPSULA]** (los glomérulos deben estar en córtex real) |
| `cx_depth`, `w_cort=exp(-cx_depth/scale)` | capa1_nefronas.py:158,174 | pesa corticales hacia depth bajo (externo) | **[DEPENDE-CÁPSULA]** (usa depth como "distancia desde la cápsula"; hoy el seno lo falsea) |
| `juxta_min_depth = umbral*JUXTA_BAND_FRAC` | capa1_nefronas.py:168-169 | banda yuxtamedular = depth alto (cerca de unión CM) | **[DEPENDE-CÁPSULA]** (asume depth = profundidad desde cápsula) |
| `glom_depth` (guardado/verif.) | capa1_nefronas.py:193,248 | verificación `depth>=umbral` | **[NEUTRO]** (diagnóstico) |
| `assign_pyramids` (por eje) | capa1_nefronas.py:109-126,196 | asigna glomérulo a pirámide por eje apex→base | **[DEPENDE-SENO vía elipsoide]**, NO vía campo depth |
| `is_cortex = region=="cortex"` → demanda | capa2_demanda.py:80,83,233 | demanda base cortex=1.0 vs medula=0.4 | **[DEPENDE-CÁPSULA]** (demanda alta debe ser córtex real) |
| `region_dom` (etiqueta por vecino) | capa3ab_peritubular.py:170,205-209 | cortex/medula para desglose de cobertura | **[DEPENDE-CÁPSULA]** (solo split cortex/medula) |
| `region_dren` (heredada de peritubular) | capa3b_venoso.py:332 | desglose cortex/medula de alcance | [NEUTRO] (heredada) |
| `pyramid_apex` (papilas) | capa3c_colector.py:105 / capa3a—no | raíces del colector (papilas) | **[DEPENDE-SENO vía elipsoide]**, NO vía campo depth |
| `hilio` | capa3a:342, capa3b:335, capa3c:316 | raíz vascular / distancia de auditoría | [NEUTRO] respecto al campo depth |

**Ningún consumidor lee la distancia-al-seno a través del campo `depth`.** Todos los usos del campo
combinado son o [DEPENDE-CÁPSULA] (quieren profundidad relativa a la cápsula) o [DEPENDE-AMBAS] (el
propio `compute_depth`). Los roles del seno que sí importan viajan por otra ruta (elipsoide del seno).

---

## PREGUNTA 2 — Papilas y médula interna (el riesgo principal)

**Cómo se calcula la geometría piramidal:** `build_pyramids()` (`capa0_dominio.py:167-215`) construye
`apex/axis/length/r_base` **solo** a partir de `SEMIEJES`, `SEMIEJES_SENO`, `CENTRO_SENO` y
`CONE_HALF_ANGLE_DEG` (constantes de módulo). Los ápices se colocan sobre la pared +Y del **elipsoide
del seno** (`:199-204`: `apex = CENTRO_SENO + u·SEMIEJES_SENO`). **No** llama a
`nearest_surface_distance`, **no** usa `depth` ni `region_label`.

**Orden real de cómputo en `main()`:** `compute_depth` (`:268`) → `cortex/medulla_mask` (`:271-272`)
→ `build_pyramids` (`:275`) → `assign_pyramids` (`:276`). Aunque `build_pyramids` se llama *después*
de `compute_depth`, sus **entradas no incluyen el campo depth** — es independiente del orden. La
afirmación previa ("depth/depth_norm y la geometría de pirámides no cambian, son previos/aparte del
umbral") **queda verificada**: la geometría piramidal es función pura de los elipsoides.

**Único acoplamiento depth↔pirámides:** `assign_pyramids` (`:218-244`) recibe `medulla_mask` y solo
asigna a un cono los puntos que son médula **y** caen dentro del cono. El **test de cono es
geometría pura**; lo único que depende de depth es la máscara "¿es médula?".

**¿Algún módulo define la médula interna o la convergencia del colector usando la distancia AL SENO
específicamente?** No a través del campo depth. El colector (`capa3c`) converge hacia
`pyramid_apex` (posición geométrica sobre el seno) — es [DEPENDE-SENO] **legítimo**, pero por la vía
del elipsoide, que **no se toca** al separar cápsula/seno en `compute_depth`. → No hay dependencia
[DEPENDE-SENO] que se rompa vía el campo combinado.

---

## PREGUNTA 3 — Simulación en memoria del fix (nada guardado)

Se recomputó una profundidad-cortical alternativa = distancia SOLO a la cápsula externa. **Aparece
un efecto no anticipado que depende de CÓMO se aplique el umbral:**

### Variante A — literal del encargo: solo-cápsula, MISMO `UMBRAL_CM=0.30` (fracción)
Al quitar el `min` con el seno, el punto más profundo pasa a medirse contra la cápsula, así que
`depth_norm` **cambia**: `28.44 mm → 51.58 mm`. Como 0.30 es *fracción* de `depth_norm`, el umbral
absoluto salta de 8.53 mm a **15.47 mm**.

| efecto | valor |
|---|---|
| grosor cortical implícito (0.30 × depth_norm) | **15.47 mm** (hoy 8.53 mm) — se **aleja** de los ~6.6 mm anatómicos |
| puntos cortex → medula | 1556 (0.78 %) |
| puntos **medula → cortex** (INESPERADO) | **43 896** |
| partición nueva (cortex / medula / piramide) | 162 171 / 14 015 / 23 814 (hoy 119 831 / 36 370 / 43 799) |
| glomérulos peri-sinusales que pasan a medula | **11 / 152** (el síntoma NO se resuelve) |

→ La variante literal es **insegura**: engorda el córtex en ~44 k puntos y apenas mueve 11 de los
152 glomérulos. Causa: `UMBRAL_CM` normalizado está acoplado a `depth_norm`, que cambia al separar.

### Variante B — fix limpio: solo-cápsula + umbral ABSOLUTO en mm
Fijando el umbral en mm sobre la distancia-a-cápsula (probado a 8.53 mm = grosor actual):

| efecto | valor |
|---|---|
| puntos cortex → medula | 7 305 (3.65 %) |
| puntos **medula → cortex** | **0** (sin efecto colateral) |
| glomérulos peri-sinusales que pasan a medula | **83 / 152** (no-peri-sinusales que flipan: **0**) |

Barrido del umbral absoluto (solo-cápsula), para ver el acercamiento a ~6.6 mm anatómicos:

| umbral mm | cortex % | peri-sinusales→medula | no-peri que flipan |
|---|---|---|---|
| 5.00 | 37.0 % | 119 / 152 | 286 |
| 6.00 | 42.9 % | 110 / 152 | 233 |
| **6.60** | **46.3 %** | **104 / 152** | 203 |
| 7.00 | 48.5 % | 102 / 152 | 168 |
| 8.00 | 53.7 % | 90 / 152 | 62 |
| 8.53 (actual) | 56.3 % | 83 / 152 | 0 |

**Matiz importante:** ni la variante limpia mueve los 152 a médula al grosor actual (solo 83);
~69 de los "152 peri-sinusales" están **también** genuinamente cerca de la cápsula (región de
esquina donde cápsula y seno se juntan; `d_capsula` desde 0.43 mm). Es decir, "152 peri-sinusal"
(superficie más cercana = seno) **sobrecuenta** el conjunto "falso córtex que debería ser médula":
el verdadero falso-córtex interior es ~83 al grosor actual. Apretar hacia 6.6 mm resuelve más
peri-sinusales (104) pero reclasifica también 203 glomérulos corticales normales a médula
(cambio de diseño mayor, no solo un fix del seno).

Efecto colateral geométrico correcto: de los puntos cortex→medula, **162 entran en un cono
piramidal** (`assign_pyramids`) → pasarían de "cortex" a "piramide_XX". Es coherente (son médula
profunda junto a la papila), no una rotura.

---

## PREGUNTA 4 — Efectos colaterales aguas abajo (predicción, sin ejecutar)

Con el fix limpio (Variante B), dado el inventario de la Pregunta 1:

**Se benefician (correcto) — [DEPENDE-CÁPSULA]:**
- Capa 1, pool de siembra `region=="cortex"` (`:156`): deja de admitir puntos interiores peri-sinusales.
- Capa 1, gradiente `w_cort` y banda yuxta (`:168-174`): `depth` pasa a significar de verdad
  "profundidad desde la cápsula"; el peso exp(−depth) y la banda yuxta se vuelven físicamente correctos.
- Capa 2, demanda `cortex=1.0` (`:83`): deja de asignar demanda cortical alta a tejido interior.
- Capa 3ab, desglose cortex/medula (`:205-209`): etiquetas más fieles.

**No se rompe nada [DEPENDE-SENO] legítimo:** papilas (`build_pyramids`), asignación de cono
(`assign_pyramids`), raíces del colector (`pyramid_apex`) y exclusión de parénquima
(`is_parenchyma`) usan el **elipsoide del seno directamente**, no el campo `depth`. Separar
cápsula/seno en `compute_depth` **no** toca esa ruta.

**Usos que REQUIEREN ajuste acompañante (por eso no es VERDE):**
1. **`UMBRAL_CM` debe dejar de ser fracción normalizada y pasar a mm absolutos** (o re-derivarse
   contra el nuevo `depth_norm`). Sin esto, la renormalización engorda el córtex en 43 896 puntos
   (Variante A). Es el ajuste que sostiene todo el fix.
2. **`CORTICAL_DEPTH_SCALE=0.12` y `JUXTA_BAND_FRAC=0.75`** (capa1) están calibrados contra la
   normalización actual (`depth_norm=28.44`). Con el campo redefinido siguen teniendo sentido pero
   **hay que retunearlos**; en particular `juxta_min_depth = umbral*JUXTA_BAND_FRAC` cambia de escala
   si `umbral` pasa a mm.
3. Regeneración en cascada de los 8 `.npz` (capa0→capa3c), ya conocida: `region_label` cambia →
   glomérulos, demanda, árboles y colector se regeneran (esperado, no es un daño).

---

## VEREDICTO DE SEGURIDAD: **[AMARILLO]**

Separar la profundidad-cortical del rol del seno es **viable y correcto a nivel de consumidores**:
ningún uso legítimo depende de la distancia-al-seno **a través del campo `depth`** (los roles reales
del seno —papilas, conos, exclusión, colector— viajan por el elipsoide, intactos), así que **no es
ROJO**. Pero **no es VERDE**: el fix literal (solo-cápsula con el mismo `UMBRAL_CM=0.30`) es inseguro
—renormaliza `depth_norm` (28→52 mm), engorda el córtex en 43 896 puntos y solo mueve 11/152
glomérulos—. El fix seguro exige **ajustes acompañantes**: (1) convertir el umbral cortico-medular a
mm absolutos, (2) retunear `CORTICAL_DEPTH_SCALE`/`JUXTA_BAND_FRAC`, (3) regenerar la cascada. Con
esos ajustes, `medula→cortex = 0` y el síntoma peri-sinusal se resuelve limpiamente.

---

## Números clave
- **Puntos que cambian de etiqueta:** Variante A (literal, insegura): 1 556 cortex→medula **pero
  43 896 medula→cortex**. Variante B (limpia, umbral absoluto 8.53 mm): **7 305 cortex→medula, 0
  medula→cortex**.
- **Nuevo grosor cortical:** Variante A: **15.47 mm** (se aleja de 6.6 mm). Variante B: el que se fije
  en mm (8.53 mm mantiene el actual; 6.6 mm → córtex 46 %, pero reclasifica 203 glomérulos normales).
- **¿Los 152 glomérulos peri-sinusales quedan en médula?** No del todo: Variante A **11/152**;
  Variante B a grosor actual **83/152** (los ~69 restantes están también genuinamente cerca de la
  cápsula, no son falso córtex).
- **Veredicto:** **[AMARILLO]** — viable con 3 ajustes acompañantes (umbral en mm, retuneo de
  constantes de capa1, regeneración de cascada).
