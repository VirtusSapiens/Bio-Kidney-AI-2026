# Diagnóstico de causa raíz — el "matiz del seno" en el campo de profundidad

**Programa:** Bio-Kidney AI 2026 · Fecha: 2026-07-05
**Naturaleza:** SOLO LECTURA / forense. No se modificó ningún `.npz`, parámetro ni lógica.
Todas las cifras se recomputan **en memoria** desde los `.npz` existentes con la geometría
exacta de `capa0_dominio.py` (`surface_radius`, elipsoide principal y elipsoide del seno).

> **Presento evidencia numérica, no conclusiones causales.** Las tres hipótesis del encargo
> (banda peri-sinusal falsa → captura glomérulos / coincide con el vacío avascular / explica el
> sesgo polar) se dejan al juicio del lector con los números delante.

## Verificación del contexto (confirmado leyendo el código)
- `compute_depth` (`capa0_dominio.py:114-128`) usa `nearest_surface_distance` (`:98-111`), que
  toma el **mínimo** entre distancia a la cápsula externa (`rsurf_main − r_main`) y distancia a la
  pared del seno (`r_seno − rsurf_seno`). **Confirmado.**
- Consecuencia confirmada: un punto pegado a la pared del seno recibe `depth` baja → se etiqueta
  `cortex` aunque esté lejos de la cápsula externa.
- **Nota metodológica:** en Capa 1 *todos* los glomérulos se siembran desde el pool
  `region == "cortex"`, así que la condición "Y está etiquetado cortex" (Pregunta 1) se cumple por
  construcción para los 1300. Por eso "peri-sinusal falso" se reduce a: **superficie más cercana =
  seno**, es decir `d_seno < d_capsula`.
- **El vacío avascular NO está guardado en ningún `.npz`** (el `capa3b_auditoria_colisiones.npz` es
  de colisiones arterio-venosas, no de cobertura). Se **recomputa** (Pregunta 2).

---

## PREGUNTA 1 — MAGNITUD de la banda peri-sinusal

Para cada uno de los 1300 glomérulos se recomputaron por separado `d_capsula` (distancia al
elipsoide externo) y `d_seno` (distancia al elipsoide de exclusión del seno). "Peri-sinusal" =
`d_seno < d_capsula`.

- **Peri-sinusales: 152 / 1300 = 11.69 %**
- "Interior real" de esos 152 (distancia a la cápsula externa, `d_capsula` [mm]):
  **min 0.43 · media 8.90 · max 21.11 mm** → una parte está a >8 mm de la cápsula externa (es decir,
  puntos interiores que quedan "cortex" solo por la cercanía al seno).

Distribución por `piramide_destino`:

| pir k | total nefronas | peri-sinusales | % peri-sinusal de la pirámide |
|---|---|---|---|
| 0 | 274 | 25 | 9.1 % |
| 1 | 164 | 15 | 9.1 % |
| 2 | 70 | 12 | 17.1 % |
| 3 | 85 | 17 | 20.0 % |
| 4 | 77 | 6 | 7.8 % |
| 5 | 73 | 8 | 11.0 % |
| 6 | 66 | 10 | 15.2 % |
| 7 | 53 | 11 | 20.8 % |
| 8 | 180 | 24 | 13.3 % |
| 9 | 258 | 24 | 9.3 % |

Umbral interpretativo declarado por el encargo (número crudo, sin veredicto): **<1 % cosmético /
>10 % estructural. El valor crudo es 11.69 %.** En *fracción* por pirámide, las centrales pequeñas
(k=2,3,6,7) muestran mayor % peri-sinusal (17–21 %) que las extremas (k=0,9: ~9 %).

---

## PREGUNTA 2 — COINCIDENCIA con el vacío avascular central

Recomputado: nube de vasos = aristas arteriales densificadas a 0.1 mm (como
`capa3_cobertura_difusion_deficit.py`) + puntos de drenaje peritubular. Para cada uno de los 200 000
puntos de parénquima (Capa 0) se calculó la distancia al vaso más cercano; el máximo es el "vacío".

| nube de vasos | n vasos | posición del vacío | dist al vaso más cercano |
|---|---|---|---|
| arterial + peritubular | 82 135 | **[−8.4, 6.3, 1.1] mm** | **6.64 mm** |
| arterial + peritubular + venoso | 249 424 | **[−8.4, 6.3, 1.1] mm** (idéntico) | **6.64 mm** |

Distancias del vacío a las referencias del seno:

| referencia | posición | distancia al vacío |
|---|---|---|
| (a) centro del seno | [0, −34, 0] | **41.17 mm** |
| (b) pared del seno (a lo largo del rayo) | — | **25.02 mm** |
| (c) hilio | [0, −30, 0] | **37.27 mm** |

Top-5 puntos peor cubiertos (posición = distancia): `[−8.4, 6.3, 1.1]=6.6`; `[−8.9, 5.6, 1.5]=6.6`;
`[10.1, −0.3, −0.6]=6.4`; `[−8.8, 11.7, 1.5]=6.4`; `[−8.8, 5.7, 1.1]=6.4`.

Observaciones numéricas (sin veredicto): el vacío está en **Y ≈ +6.3** (lado opuesto al seno, que
está en Y ≈ −34/−18); la magnitud (6.64 mm) es mayor que el "~4.9 mm" mencionado en el encargo;
añadir el árbol venoso **no cambia** la posición ni la distancia. El vacío está a >25 mm de la pared
del seno.

---

## PREGUNTA 3 — CORRELACIÓN con el sesgo polar de drenaje

Fracción de drenaje = nefronas asignadas a cada pirámide / 1300. `d_papila_seno` = distancia de
`pyramid_apex[k]` al centro del seno. `n_perisinusal` = glomérulos peri-sinusales de esa pirámide.

| k | frac drenaje | d_papila_seno (mm) | n_perisinusal |
|---|---|---|---|
| 0 | 21.1 % | 17.35 | 25 |
| 1 | 12.6 % | 16.57 | 15 |
| 2 | 5.4 % | 15.97 | 12 |
| 3 | 6.5 % | 15.55 | 17 |
| 4 | 5.9 % | 15.34 | 6 |
| 5 | 5.6 % | 15.34 | 8 |
| 6 | 5.1 % | 15.55 | 10 |
| 7 | 4.1 % | 15.97 | 11 |
| 8 | 13.8 % | 16.57 | 24 |
| 9 | 19.8 % | 17.35 | 24 |

- **corr(frac_drenaje, d_papila_seno) = 0.937** (positiva: más drenaje ↔ papila *más lejos* del
  centro del seno; las pirámides polares k=0,9 tienen la papila más alejada, 17.35 mm, y el mayor
  drenaje).
- Posición axial de las papilas (X): `[−11.9, −9.2, −6.6, −4.0, −1.3, 1.3, 4.0, 6.6, 9.2, 11.9]`
  → k=0 y k=9 son los **polos axiales**; k=4,5 son las centrales.

---

## PREGUNTA 4 — CONTROL: ¿explicación del sesgo INDEPENDIENTE del seno?

El drenaje por pirámide **es** la asignación de siembra de Capa 1 (`piramide_destino`), así que el
sesgo ya está presente en la siembra. Para separar "el seno causa el sesgo" de "el sesgo es
geometría piramidal", se comparó con dos medidas geométricas independientes del campo de profundidad
del seno: volumen del cono piramidal (de `pyramid_r_base`, `pyramid_length`) y catchment cortical
(cada punto de corteza asignado al eje piramidal más cercano).

| k | frac drenaje | vol_cono (mm³) | n_med_pts (territorio) | cortex_catchment | catch % |
|---|---|---|---|---|---|
| 0 | 21.1 % | 8740.9 | 6428 | 23751 | 19.8 % |
| 1 | 12.6 % | 7398.0 | 5548 | 16207 | 13.5 % |
| 2 | 5.4 % | 6493.6 | 3627 | 6340 | 5.3 % |
| 3 | 6.5 % | 5932.9 | 3307 | 6739 | 5.6 % |
| 4 | 5.9 % | 5664.3 | 2974 | 6976 | 5.8 % |
| 5 | 5.6 % | 5664.3 | 2984 | 7153 | 6.0 % |
| 6 | 5.1 % | 5932.9 | 3242 | 6602 | 5.5 % |
| 7 | 4.1 % | 6493.6 | 3601 | 6227 | 5.2 % |
| 8 | 13.8 % | 7398.0 | 5617 | 16081 | 13.4 % |
| 9 | 19.8 % | 8740.9 | 6471 | 23755 | 19.8 % |

- **corr(frac_drenaje, cortex_catchment) = 0.994**
- **corr(frac_drenaje, vol_cono) = 0.950**

Observaciones numéricas (sin veredicto): las pirámides polares (k=0,9) tienen ~1.5× el volumen de
cono y ~3.7× el catchment cortical de las centrales (k=4,5); el reparto de drenaje sigue el
catchment cortical casi 1:1 (r=0.994). La geometría del cono es simétrica respecto al eje X (k y
9−k comparten `vol_cono` y `d_papila_seno`), reflejando el reparto `xf=(f−0.5)·2` de `build_pyramids`.

---

## Los 4 números-clave

1. **% peri-sinusal:** **11.69 %** (152 / 1300 glomérulos con superficie más cercana = pared del seno).
2. **d_vacío_a_seno:** vacío en [−8.4, 6.3, 1.1]; **25.02 mm a la pared del seno**, 41.17 mm al centro
   del seno, 37.27 mm al hilio (invariante a añadir el árbol venoso; distancia al vaso 6.64 mm).
3. **¿drenaje correlaciona con d_papila_seno?** Sí, **r = +0.937** (más drenaje ↔ papila más lejos del
   seno); pero el drenaje correlaciona **más fuerte con el catchment cortical (r = 0.994)** y con el
   volumen del cono (r = 0.950), ambos independientes del campo de profundidad del seno.
4. **¿el sesgo preexiste en la siembra?** Sí: el drenaje por pirámide es la asignación de Capa 1, y
   queda predicho casi 1:1 por la geometría del cono/catchment (r = 0.994), simétrica en el eje X.

*(Todos los datos anteriores son recomputables desde los `.npz` presentes; el único ítem no
almacenado —el vacío avascular— se recomputó explícitamente y se indicó como tal.)*
