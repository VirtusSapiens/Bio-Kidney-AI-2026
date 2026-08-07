> # ⚠️ DOCUMENTO PARCIALMENTE SUPERADO — LEER CON EL FILTRO DE ABAJO
>
> **Revisado por `_auditoria/2026-08-12_reservas_pw_region2d.md` (2026-08-12).** Ese informe somete
> las conclusiones de éste a tres reservas metodológicas y **retira algunas**. Este documento **no es
> el estado actual**; se conserva íntegro y sin modificar por su valor forense.
>
> ## SUPERADO — no citar de este documento
>
> - **El intervalo `CENTRO_SENO[1] ∈ [−32.5, −27.8]`** (Tarea 1b) y toda cifra derivada de él.
>   **Mal planteado:** se cumple la identidad `PW_max = |−B_SEMI − (cy + sy)|`, luego PW restringe
>   sólo la **suma** `cy + sy ∈ [−17.4, −11.8]`. El intervalo publicado es un **corte en `sy = 16`**,
>   valor tan sin anclar como `cy` (puestos 1 y 2 de los 10 pendientes). Y `sx`, `sz` tampoco están
>   anclados: **el espacio real es 4D**. La formulación correcta es una **región 2D** (rebanada de
>   una 4D), tabulada en `_auditoria/2026-08-12`, Reserva 2c.
> - **La comparación «PW del gemelo vs Glodny» sin reserva de estadístico** (Tarea 2). Un **máximo
>   espacial sobre una geometría** no es comparable con una **media ± s.d. sobre n = 2068 riñones**.
>   Error de magnitud de la misma clase que el corregido por el método B (ENTRADA 032 §2).
> - **La premisa sobre el protocolo de medida de Glodny.** El paper (BMC Urology 2009;9:19,
>   PMC2813848) especifica fase arterial y plano axial, pero **no** el nivel anatómico, **ni** el
>   número de medidas por riñón, **ni** si es localización estandarizada o promedio.
>   **[PENDIENTE DE ANCLA]** los tres.
> - **La condición de ápices corticales como criterio de calibración** (Tareas 1a-1b). Es
>   **redundante**: con `sy = 16` su umbral es `cy ≥ −32.50`, menos restrictivo que el de Caglar
>   (`cy ≥ −31.30`). **No fija ningún borde.** Sigue siendo un defecto real; no es un pilar.
> - **⚠ AUTOCORRECCIÓN — «el volumen de Emamian es exactamente el desempate» (Tarea 3, B3) es
>   INEXACTO.** Lo escribí yo en este documento y lo corrijo yo en `_auditoria/2026-08-12`
>   (Tarea Final). **El volumen es UNA ecuación sobre CUATRO incógnitas** `(cy, sy, sx, sz)`: reduce
>   el espacio de 4D a 3D, no a un punto. **Lo que reduce la región a un punto son las TRES
>   DIMENSIONES de la CEA** (`sx`, `sy`, `sz`); con ellas ancladas, Glodny fija `cy` por la suma. Si
>   Emamian aportase además la posición, el sistema quedaría sobredeterminado — comprobación de
>   consistencia, no grado de libertad.
>
> ## VIGENTE — sigue siendo la referencia
>
> - **`CENTRO_SENO` es la variable dominante frente a `SEMIEJES_SENO`** (Tarea 1a). El barrido
>   completo cy ∈ [−34, −28] es geométrico y no depende de ninguna reserva.
> - **El defecto de `capa4_colector_alto.py:449`** (Tarea 1c): la fórmula del seno escrita como cadena
>   literal mientras `:340` computa desde el `.npz`. Hoy es correcta **por coincidencia, no por
>   construcción**; cualquier cambio de seno la vuelve falsa en silencio. El texto ya está versionado
>   en `09_paper_vascular/auditoria_capa4_calicial.md:19` y `00_bitacora/BITACORA.md:998`. **Debe
>   corregirse aunque no se toque `CENTRO_SENO`.** El f-string propuesto sigue siendo válido.
> - **La especificación funcional de la partición** (Tarea 2) y su resultado: **0.0000 % de médula
>   huérfana (82 297 / 82 297)** con el seno actual, frente al 47.050 % de los conos. `pyramid_r_base`
>   y `cone_half_angle_deg` pierden su único consumo interno (`capa0_dominio.py:374`) y tienen **cero
>   consumidores externos**. Independiente de todas las reservas.
> - **El inventario de decisiones bloqueadas** (Tarea 3), **con la corrección de B3 indicada arriba**.
> - **Las cuatro afirmaciones de papilas declaradas medulares** (Tarea 4) y las notas propuestas.
>
> **Referencia vigente:** `08_gemelo_digital/_auditoria/2026-08-12_reservas_pw_region2d.md`;
> fuentes en `04_literatura/anclas_seno_renal.md`; estado en `00_bitacora/BITACORA.md`, ENTRADA 033.
>
> **Aviso de trazabilidad:** «PW» **no existe como código en el repo** — es una construcción de
> auditoría introducida en `_auditoria/2026-08-09`. **No confundir con `depth_cortical_mm`**, que sí
> es un campo del `.npz` y mide otra cosa (cápsula → punto, en todo el parénquima, incluidos los polos
> donde no hay seno). Definición operativa en `_auditoria/2026-08-12`.

---

# Auditoría — `CENTRO_SENO` consolidado, especificación de la partición, y papilas declaradas medulares

**Programa:** Bio-Kidney AI 2026 · **Capa:** 0 · **Fecha:** 2026-08-10

> **NATURALEZA: SOLO LECTURA / DIAGNÓSTICO.**
> No se modificó ningún `.py`, ningún `.npz`, ningún parámetro ni ninguna lógica. No se ejecutó
> `capa0_dominio.py`. Todas las cifras se recomputan **en memoria** desde `capa0_dominio.npz` y
> funciones puras importadas de `capa0_dominio` (el guard `if __name__ == "__main__":` de
> `capa0_dominio.py:610` impide que el import dispare `main()`).
>
> **Ninguna configuración de seno y ningún esquema de partición está aplicado.** La Tarea 2 es una
> especificación funcional, no código. Los textos de la Tarea 4 son propuestas de redacción no
> escritas en ningún archivo.
>
> Continúa `_auditoria/2026-08-07`, `_auditoria/2026-08-08` y `_auditoria/2026-08-09`.

---

## 0. Chequeo de integridad — cifras externas no fiables

Se advierte de un documento externo que presenta dimensiones **reconstruidas** del área ecogénica
central de Emamian 1993 como si fueran datos del paper (L = 4.5, W = 1.3, T = 1.2 cm; V = 3.68 cm³).
**Esos números no están en la fuente y no se usan en este informe.**

Búsqueda en el repo:

```bash
grep -rn "3\.68\|3,68" --include=*.md --include=*.py --include=*.txt --include=*.json . --exclude-dir=.git --exclude-dir=.venv
grep -rn "\bCEA\b" --include=*.md --include=*.py . --exclude-dir=.git
grep -rn "ecogénic\|ecogenic\|echogenic" --include=*.md --include=*.py . --exclude-dir=.git
```

**Resultado: las cifras NO están en el repo.**

| patrón | resultado |
|---|---|
| `3.68` como volumen de seno | **NO ENCONTRADO.** Los únicos aciertos de la cadena `3.68` son coordenadas del árbol vascular sin relación con el seno: `02_vascular_cco/renal_data_v1.json:2409` (`"y2_mm": 3.6887`), `:3377` (`"x2_mm": 3.6804`), `:4276`, `:4293`, `:5098` (`"radio_um": 123.68`) |
| `CEA` | **NO ENCONTRADO** (cero apariciones en todo el repo) |
| `ecogénica` / `echogenic` | Aparece **sólo** en `08_gemelo_digital/_auditoria/2026-08-09_seno_apices_particion.md:23,432,433,435,452`, y **sin ningún número**: describe qué dato haría falta, no lo aporta |
| Menciones de Emamian | 7, todas en `00_bitacora/BITACORA.md:1647` y `09_paper_vascular/auditoria_correspondencia_anatomica.md:8,41,42,43,126,140`, y **todas referidas exclusivamente a los tres semiejes del órgano** |

**El repo está limpio.** Las dimensiones y el volumen de la CEA siguen **[PENDIENTE DE ANCLA]** y este
informe no los estima por ninguna vía.

---

## TAREA 1 — `CENTRO_SENO[1]` consolidado

### (a) Barrido completo, cy ∈ [−34, −28] en pasos de 0.5

```bash
.venv/bin/python -c "import numpy as np, capa0_dominio as C
SS=np.array([22.,16.,11.]); A=C.SEMIEJES
def ad(cyv):
    o=[]
    for i in range(10):
        xf=((i+0.5)/10-0.5)*2.0; z=1.0 if i%2==0 else -1.0
        ux=xf*0.60; uz=z*0.40; uy=np.sqrt(max(0,1-ux*ux-uz*uz))
        o.append(float(C.capsule_distance(np.array([[ux*SS[0], cyv+uy*SS[1], uz*SS[2]]]))[0]))
    return np.array(o)
rng=np.random.default_rng(2026); n=600000
u=rng.normal(size=(n,3)); u/=np.linalg.norm(u,axis=1,keepdims=True)
P=(u*rng.uniform(0,1,(n,1))**(1/3))*A; VM=4/3*np.pi*np.prod(A)
def vol(cyv): return VM*np.mean(C.ellipsoid_level(P,np.array([0.,cyv,0.]),SS)<1.0)/1000
print('  cy   apices<6.6  depth_min  PW_max  vol_seno_mL  PW en [12.6,18.2]?')
for cy in np.arange(-34.0,-27.99,0.5):
    a=ad(cy); pw=cy+SS[1]+30
    print('%7.1f %10d %11.3f %8.2f %12.3f   %s'%(cy,(a<6.6).sum(),a.min(),pw,vol(cy),'SI' if 12.6<=pw<=18.2 else 'no'))"
```

| `cy` | ápices < 6.6 | `depth` mín | PW_max | volumen seno | PW ∈ [12.6, 18.2] |
|---|---|---|---|---|---|
| **−34.0 (actual)** | **4** | 5.295 | 12.00 | 3.949 mL | no |
| −33.5 | 2 | 5.674 | 12.50 | 4.259 mL | no |
| −33.0 | 2 | 6.045 | 13.00 | 4.572 mL | sí |
| −32.5 | 2 | 6.406 | 13.50 | 4.914 mL | sí |
| **−32.0** | **0** | 6.757 | 14.00 | 5.246 mL | sí |
| −31.5 | 0 | 7.098 | 14.50 | 5.587 mL | sí |
| −31.0 | 0 | 7.428 | 15.00 | 5.926 mL | sí |
| −30.5 | 0 | 7.748 | 15.50 | 6.280 mL | sí |
| −30.0 | 0 | 8.058 | 16.00 | 6.635 mL | sí |
| −29.5 | 0 | 8.357 | 16.50 | 6.982 mL | sí |
| −29.0 | 0 | 8.645 | 17.00 | 7.351 mL | sí |
| −28.5 | 0 | 8.923 | 17.50 | 7.732 mL | sí |
| −28.0 | 0 | 9.190 | 18.00 | 8.105 mL | sí |

Relación cerrada: `PW_max = CENTRO_SENO[1] + SEMIEJES_SENO[1] + B_SEMI = cy + 46`. El volumen es
estimación Monte Carlo con `SEED = 2026`, N = 600 000; el valor en cy = −34 (3.949 mL) coincide con
el canónico 3.9599 mL de ENTRADA 033 §2 dentro del ruido de muestreo.

### (b) ¿Existe un `cy` que cumpla las tres condiciones?

Frontera exacta del cruce 2 → 0 ápices corticales, por bisección:

```bash
.venv/bin/python -c "import numpy as np, capa0_dominio as C
SS=np.array([22.,16.,11.])
def nc(cyv):
    o=[]
    for i in range(10):
        xf=((i+0.5)/10-0.5)*2.0; z=1.0 if i%2==0 else -1.0
        ux=xf*0.60; uz=z*0.40; uy=np.sqrt(max(0,1-ux*ux-uz*uz))
        o.append(float(C.capsule_distance(np.array([[ux*SS[0], cyv+uy*SS[1], uz*SS[2]]]))[0]))
    return int((np.array(o)<6.6).sum())
lo,hi=-32.5,-32.0
for _ in range(30):
    m=(lo+hi)/2
    if nc(m)==0: lo=m
    else: hi=m
print('frontera 2->0 apices: cy = %.4f  (PW_max %.4f)'%(hi,hi+46))
print('interseccion: cy in [%.4f, %.4f]'%(max(hi,12.6-46),18.2-46))
print('PW exacto 15.4 -> cy = %.4f'%(15.4-46))"
```

```
frontera 2->0 apices: cy = -32.5000  (PW_max 13.5000)
interseccion: cy in [-32.5000, -27.8000]
PW exacto 15.4 -> cy = -30.6000
```

**Sí. Intervalo conjunto: `cy ∈ [−32.5, −27.8]`.**

| condición | restricción |
|---|---|
| 0 ápices corticales | cy ≥ **−32.5** |
| PW_max ≥ 12.6 (Glodny −1 s.d.) | cy ≥ −33.4 |
| PW_max ≤ 18.2 (Glodny +1 s.d.) | cy ≤ **−27.8** |
| **intersección** | **cy ∈ [−32.5, −27.8]**, anchura **4.7 mm** |
| PW_max = 15.4 exacto | **cy = −30.6** — cae dentro, aproximadamente en el centro del intervalo |

**Sobre "rango razonable" de volumen: no puedo evaluarlo.** El volumen del seno recorre 4.914 →
8.105 mL a lo largo del intervalo, pero **no existe en el repo ni entre las fuentes a la vista ningún
valor de referencia del volumen del seno renal**. Marcar un subintervalo como "razonable" exigiría
precisamente el dato de Emamian que falta. **[PENDIENTE DE ANCLA]** — el intervalo de arriba está
acotado por los ápices y por Glodny, **no** por volumen.

Lo único que sí puedo decir con respaldo: cy = −30.6 multiplica el volumen actual por **≈ 1.57**
(3.949 → ~6.21 mL), y eso reabre ENTRADA 033, que hoy declara 3.9599 mL como volumen excluido.

### (c) ¿Rompe algo aguas arriba?

Dos consumidores de `centro_seno`, con comportamiento **opuesto**.

**`capa3_cobertura_difusion_deficit.py` — paramétrico, seguro.**

```bash
grep -n "centro_seno\|semiejes_seno" 08_gemelo_digital/capa3_cobertura_difusion_deficit.py
```
```
87:    centro_seno = d0["centro_seno"].astype(np.float64)
88:    semiejes_seno = d0["semiejes_seno"].astype(np.float64)
145:    seno = en_seno(coords, centro_seno, semiejes_seno)          # todo el dominio
167:          f"(elipsoide centro {centro_seno.tolist()} semiejes {semiejes_seno.tolist()})")
183:          f"(dist al seno {np.linalg.norm(centro_geo-centro_seno):.1f}mm -> fuera del seno)")
```

Lee del `.npz`, calcula con lo leído (`:145`) e **imprime los valores reales** (`:167`). Cambiar
`CENTRO_SENO` se propaga correctamente. **No rompe nada.**

**`capa4_colector_alto.py` — DEFECTO: fórmula hardcodeada en el informe generado.**

Calcula con los valores del `.npz`, `:340`:

```python
    val = np.array([elipsoide_val(nodos[i], CENTRO_SENO, SEMIEJES_SENO) for i in range(M)])
```

pero **escribe una cadena literal** con los valores viejos, `:449`:

```python
    L.append("Test elipsoide: ((x)/22)^2 + ((y+34)/16)^2 + ((z)/11)^2 <= 1.0\n\n")
```

**Si `CENTRO_SENO` cambia, Capa 4 computa con el seno nuevo y declara por escrito el viejo.** El
informe afirmaría haber aplicado un test que no aplicó. Es un fallo silencioso: no hay excepción, no
hay `[REVISAR]`, sólo una línea de documentación que deja de ser cierta.

**Y ese texto ya está versionado en dos sitios**, de modo que el error se propagaría a documentos del
repo:

- `09_paper_vascular/auditoria_capa4_calicial.md:19` — `Test elipsoide: ((x)/22)^2 + ((y+34)/16)^2 + ((z)/11)^2 <= 1.0`
- `00_bitacora/BITACORA.md:998` — `- **(3) Contención en el seno** (elipsoide `((x)/22)²+((y+34)/16)²+((z)/11)²≤1`): cáliz menor 10/10, cáliz mayor 2/2, pelvis 1/1, uréter 1/1 → **100 % dentro [OK]**.`

**Corrección propuesta (texto, sin aplicar):** sustituir `capa4_colector_alto.py:449` por una cadena
formateada a partir de las variables leídas, de forma que la fórmula impresa no pueda divergir de la
computada:

```python
    L.append(f"Test elipsoide: ((x-{CENTRO_SENO[0]:g})/{SEMIEJES_SENO[0]:g})^2 + "
             f"((y-{CENTRO_SENO[1]:g})/{SEMIEJES_SENO[1]:g})^2 + "
             f"((z-{CENTRO_SENO[2]:g})/{SEMIEJES_SENO[2]:g})^2 <= 1.0\n\n")
```

Esta corrección es **independiente** de que se cambie o no `CENTRO_SENO`: hoy la línea es correcta por
coincidencia, no por construcción.

---

## TAREA 2 — Especificación funcional de la partición (sin implementar)

### Qué reemplaza exactamente

**Función objetivo:** `assign_pyramids`, declarada en `capa0_dominio.py:356`.

Firma actual (`:356`):

```python
def assign_pyramids(coords, medulla_mask, apex, axis, length, r_base):
```

Criterio actual, docstring `:359-363`:

```
    Un punto pertenece al cono i si su proyeccion axial t cae en [0, L_i]
    y su distancia radial al eje es menor que el radio del cono en t
    (radio crece linealmente de 0 en el apice a r_base en la base).
    Si cae en varios conos se asigna al de mejor ajuste (menor radio
    normalizado).
```

**Especificación propuesta.** Sustituir la contención en cono por asignación al **ápice más cercano**:

| aspecto | actual | propuesto |
|---|---|---|
| **entradas** | `coords`, `medulla_mask`, `apex`, `axis`, `length`, `r_base` | `coords`, `medulla_mask`, `apex` |
| **criterio** | contención en cono + desempate por radio normalizado (`:377`) | distancia euclídea mínima a `apex[i]` |
| **salida** | `best_pyr` (N,) int32, `−1` si no cae en ningún cono | `pyr_id` (N,) int32, `−1` **sólo** fuera de `medulla_mask` |
| **escape** | sí — genera la etiqueta `'medulla'` genérica | **no existe dentro de la médula** |
| **líneas afectadas** | `:356-382` (cuerpo completo), `:415-416` (llamada) | mismas |

`build_pyramids` (`:305`) **no se toca**: sigue produciendo `apex`, `axis`, `length`. Sólo deja de ser
necesario su cuarto retorno.

### Qué pasa con `pyramid_r_base` y `cone_half_angle_deg`

Trazado completo de sus usos:

```bash
grep -n "r_base\|half\b\|CONE_HALF_ANGLE_DEG" capa0_dominio.py
```

| línea | uso | bajo partición |
|---|---|---|
| `:69` | `CONE_HALF_ANGLE_DEG = 22.0` (declaración) | **queda sin ningún uso** |
| `:322` | `half = np.deg2rad(CONE_HALF_ANGLE_DEG)` | muere |
| `:352` | `r_base = length * np.tan(half)` | muere |
| `:353` | `return apex, axis, length, r_base` | el cuarto retorno sobra |
| `:374` | `cone_r = (t / length[i]) * r_base[i]` — **único consumo real** | **desaparece con el criterio de cono** |
| `:448` | `cone_half_angle_deg=np.float64(CONE_HALF_ANGLE_DEG)` (persistencia) | persiste un valor sin significado |
| `:456` | `pyramid_r_base=r_base.astype(np.float64)` (persistencia) | ídem |

Consumidores externos:

```bash
cd 08_gemelo_digital && grep -rn '\["pyramid_r_base"\]\|\["cone_half_angle_deg"\]' *.py | wc -l
```
→ **`0`**

**Ambas claves pueden retirarse del `.npz` sin romper a nadie.** Y `CONE_HALF_ANGLE_DEG` es hoy la
única constante **[SIN DECLARAR]** de `capa0_dominio.py` (informe del 2026-08-07, Tarea 4b): la
partición la **elimina** en lugar de obligar a anclarla. Es la única opción evaluada que reduce el
número de parámetros libres en vez de aumentarlo.

### Confirmación numérica bajo la geometría ACTUAL, sin tocar el seno

Éste es el punto de la tarea: la partición debe funcionar **con el seno tal como está hoy**, es decir
sin depender del ancla pendiente.

```bash
.venv/bin/python -c "import numpy as np; d=np.load('capa0_dominio.npz',allow_pickle=True); co=d['coords'].astype(np.float64); dm=d['depth_cortical_mm'].astype(np.float64); lab=d['region_label']; ap=d['pyramid_apex']; med=dm>=6.6; M=co[med]; lb=np.linalg.norm(M[:,None,:]-ap[None,:,:],axis=2).argmin(axis=1); c=np.array([int((lb==i).sum()) for i in range(10)]); print('medula (depth>=6.6):',int(med.sum())); print('asignados:',int((lb>=0).sum()),'| HUERFANOS:',int((lb<0).sum()),'= %.4f %%'%(100*(lb<0).mean())); print('celdas no vacias: %d/10'%len(np.unique(lb))); print('conteo',list(c)); print('%share',[round(100*x/c.sum(),2) for x in c]); print('min %d max %d ratio %.3f'%(c.min(),c.max(),c.max()/c.min())); print('conos actuales, huerfana: %.3f %%'%(100*(lab=='medulla').sum()/med.sum()))"
```

```
medula (depth>=6.6): 82297
asignados: 82297 | HUERFANOS: 0 = 0.0000 %
celdas no vacias: 10/10
conteo [13439, 9251, 6226, 6100, 6020, 6038, 6004, 6161, 9372, 13686]
%share [16.33, 11.24, 7.57, 7.41, 7.31, 7.34, 7.3, 7.49, 11.39, 16.63]
min 6004 max 13686 ratio 2.279
conos actuales, huerfana: 47.050 %
```

| métrica | conos (actual) | **partición** |
|---|---|---|
| médula huérfana | **47.050 %** | **0.0000 %** |
| puntos asignados / médula total | 43 576 / 82 297 | **82 297 / 82 297** |
| celdas no vacías | 10/10 | 10/10 |
| solape (≥2 asignaciones) | 21.78 % | **imposible por construcción** |
| ratio max/min | 1.609 | 2.279 |
| ratio polar/central | ~1.45 | **1.878** |

**La partición da 0 % de médula huérfana con el seno en `[0, −34, 0] / [22, 16, 11]`, es decir sin
tocar nada de lo que está bloqueado por Emamian.** Ése era el punto a demostrar.

**Dos salvedades, dichas claramente:**

1. El ratio polar/central 1.878 se obtiene con los ápices **actuales**, cuatro de los cuales están en
   territorio córtex (Tarea 4). Si el seno se corrige, ese número cambiará. Es el estado de hoy, no
   una predicción.
2. La partición **no resuelve** el déficit de médula del informe del 2026-08-09 (5.4 mm en el eje
   hiliar frente a los 8.8 que implica Glodny). Repartir mejor una médula delgada no la engrosa. Ese
   defecto sólo lo toca el seno, y sigue bloqueado.

---

## TAREA 3 — Inventario de lo bloqueado por el volumen de la CEA

| # | decisión bloqueada | qué se decidiría con el dato | si nunca se consigue |
|---|---|---|---|
| **B1** | Anclar `SEMIEJES_SENO` | Las tres dimensiones pasarían de [SUPUESTO DECLARADO] a [ANCLADA]. Hoy son los puestos 1 y 2 de los 10 pendientes (`09_paper_vascular/auditoria_correspondencia_anatomica.md:128-129`) | **Se declara supuesto y se sigue.** Es lo que ya hace ENTRADA 033. **No bloquea la cascada** |
| **B2** | Anclar `CENTRO_SENO` | Ídem, más la posición dentro del contorno renal | **Se declara supuesto y se sigue. No bloquea la cascada** |
| **B3** | **Elegir el punto dentro de `cy ∈ [−32.5, −27.8]`** | El volumen de la CEA sería el desempate: el intervalo recorre 4.914 → 8.105 mL y todas sus posiciones cumplen Glodny igual de bien | **Bloquea de facto la calibración fina, pero NO la cascada.** Ver abajo |
| **B4** | Cerrar ENTRADA 033 declarando el volumen de exclusión como magnitud anclada | El 3.9599 mL (o su sucesor) pasaría a ser comparable con una referencia | **Se declara supuesto y se sigue.** ENTRADA 033 ya está redactada así |
| **B5** | Decidir si la asimetría lobar (Opción F) exige seno asimétrico | Si Emamian reporta asimetría supero-inferior de la CEA | **Se declara supuesto.** [PENDIENTE DE ANCLA] incluso sobre si la fuente lo contiene |

### Lo que cambia respecto del informe del 2026-08-09

En aquel informe escribí que B3 era **el nudo** y que sin él no se podía avanzar. **El barrido de la
Tarea 1a lo matiza:** existe un intervalo de 4.7 mm de anchura donde **todas** las posiciones cumplen
las dos condiciones anclables (0 ápices corticales, PW dentro del rango de Glodny). El dato de
Emamian elegiría **un punto dentro** de ese intervalo; su ausencia no impide elegir **alguno**.

**Regla propuesta si el dato nunca llega** (texto, no aplicada): adoptar `cy = −30.6` como
**supuesto declarado**, con dos justificaciones que sí tienen respaldo:

- es el valor que hace `PW_max = 15.4` exacto, la media de Glodny 2009;
- cae aproximadamente en el centro del intervalo admisible [−32.5, −27.8], que es la elección de
  máxima robustez frente al error de la propia estimación.

Y declararlo en la bitácora como **[SUPUESTO DECLARADO, calibrado contra Glodny PW]**, nunca como
anclado a Emamian.

### Qué bloquea realmente la cascada, y qué no

| | bloquea la cascada | no la bloquea |
|---|---|---|
| **D5** — adoptar la partición | | **✓** — 0 % huérfana con el seno actual (Tarea 2) |
| **D7** — corregir la declaración papila-medular | | **✓** — es documental (Tarea 4) |
| Corregir `capa4_colector_alto.py:449` | | **✓** — independiente del valor del seno (Tarea 1c) |
| **B3** — elegir `cy` | **✗ sólo si se exige anclaje**; con supuesto declarado, no | |
| Regenerar Capa 0 y las 9 `.npz` | depende de B3 **o** de aceptar el supuesto | |

**Ninguna de las cinco decisiones bloqueadas impide avanzar si se acepta declararlas como supuesto.**
Lo que el dato de Emamian compraría es convertir un supuesto declarado en un anclaje — que es
exactamente lo que el proyecto ha estado distinguiendo desde ENTRADA 032.

---

## TAREA 4 — Papilas declaradas medulares

### Las cuatro afirmaciones localizadas

| # | ruta:línea | texto literal |
|---|---|---|
| **1** | `08_gemelo_digital/capa4_colector_alto.py:351` | `    pap_val = val[idx_papila]  # papilas: interfaz medula/seno, exentas` |
| **2** | `08_gemelo_digital/capa4_colector_alto.py:395` | `    print(f"        papila_junction (EXENTA, interfaz medula/seno): "` |
| **3** | `08_gemelo_digital/capa4_colector_alto.py:453-456` | `` L.append(f"`papila_junction` (10) es nodo-interfaz medula/seno, **EXENTA** del test " f"(se sienta sobre la pared +Y del seno). …") `` |
| **4** | `09_paper_vascular/auditoria_capa4_calicial.md:28` | `` `papila_junction` (10) es nodo-interfaz medula/seno, **EXENTA** del test (se sienta sobre la pared +Y del seno). Valor del elipsoide en las papilas: [1.000, 1.000]; polares k=0=1.000, k=9=1.000 (~1.0 = sobre la pared del seno, esperado). `` |

Una quinta, relacionada, en la bitácora: `00_bitacora/BITACORA.md:998` — `Las 10 `papila_junction` son **interfaz médula/seno, exentas**`.

**#4 y #5 son salida generada por #3 y ya están versionadas.** Es exactamente el tipo de afirmación
que llega al preprint: está redactada, tiene números, y se lee como resultado.

### El hecho que las contradice

```bash
.venv/bin/python -c "import numpy as np, capa0_dominio as C; d=np.load('capa0_dominio.npz',allow_pickle=True); da=C.capsule_distance(d['pyramid_apex']); print('depth_cortical_mm de las 10 papilas:',np.round(da,3)); print('umbral GROSOR_CORTICAL_MM =',float(C.GROSOR_CORTICAL_MM)); print('papilas con depth<6.6:',int((da<6.6).sum()),'->',list(np.where(da<6.6)[0]))"
```
```
depth_cortical_mm de las 10 papilas: [5.295 6.428 7.187 7.655 7.879 7.879 7.655 7.187 6.428 5.295]
umbral GROSOR_CORTICAL_MM = 6.6
papilas con depth<6.6: 4 -> [0, 1, 8, 9]
```

El test de contención en el seno da `val ≈ 1.000` para las diez —están sobre la pared del seno, como
dice el texto— pero **eso no las hace medulares**. Son dos condiciones distintas y sólo se comprueba
una. Las papilas k = 0, 1, 8, 9 están simultáneamente **sobre la pared del seno** y **en territorio
córtex** por el umbral de Capa 0.

### Notas propuestas (texto, SIN APLICAR)

**Para #1 — `capa4_colector_alto.py:351`.** Ampliar el comentario:

```python
    pap_val = val[idx_papila]  # papilas: interfaz medula/seno, exentas
    # NOTA (2026-08-10): "interfaz medula/seno" describe la INTENCION geometrica.
    # El test de contencion (val ~ 1.0) comprueba que la papila esta SOBRE la pared
    # del seno, NO que este en territorio medular. Bajo la geometria actual
    # (CENTRO_SENO=[0,-34,0]) las papilas k=0,1,8,9 tienen depth_cortical_mm de
    # 5.295 y 6.428 mm, por DEBAJO del umbral GROSOR_CORTICAL_MM=6.6 de Capa 0:
    # 4 de 10 caen hoy en territorio CORTEX. Ver _auditoria/2026-08-09 y 2026-08-10.
```

**Para #2 — `capa4_colector_alto.py:395`.** Añadir a la salida de consola, tras la línea existente:

```python
    print(f"        AVISO: 'interfaz medula/seno' NO se comprueba. 4 de 10 papilas "
          f"(k=0,1,8,9) estan hoy a <6.6 mm de la capsula = territorio CORTEX.")
```

**Para #3 — `capa4_colector_alto.py:453-456`.** Añadir un párrafo al informe generado, inmediatamente
después del bloque actual:

```python
    L.append("\n> **AVISO DE ALCANCE (2026-08-10).** El test de contencion comprueba que las "
             "papilas esten SOBRE la pared del seno (val ~ 1.0), **no** que esten en territorio "
             "medular. Son condiciones distintas. Bajo la geometria actual, las papilas "
             "**k=0, 1, 8, 9** tienen `depth_cortical_mm` de **5.295** y **6.428 mm**, por debajo "
             "del umbral `GROSOR_CORTICAL_MM = 6.6` de Capa 0: **4 de 10 caen en territorio "
             "cortex**. La expresion 'interfaz medula/seno' describe la intencion de diseno, no un "
             "hecho comprobado. Diagnostico y causa en "
             "`08_gemelo_digital/_auditoria/2026-08-09_seno_apices_particion.md`.\n\n")
```

**Para #4 — `09_paper_vascular/auditoria_capa4_calicial.md:28`.** Como es salida generada, la
corrección de fondo es #3 y regenerar. Mientras tanto, insertar inmediatamente después de la línea 28,
**sin borrar el texto existente**:

```markdown
> **⚠ AVISO AÑADIDO 2026-08-10 — el texto anterior es de alcance mayor que lo comprobado.**
> El test comprueba `val ≈ 1.0` (papila **sobre la pared del seno**), no que la papila sea
> **medular**. Bajo la geometría vigente, las papilas **k = 0, 1, 8, 9** tienen
> `depth_cortical_mm` = **5.295** y **6.428 mm**, por debajo del umbral
> `GROSOR_CORTICAL_MM = 6.6` (`capa0_dominio.py:61`): **4 de 10 están hoy en territorio córtex**.
> Reproducible con:
> ```bash
> .venv/bin/python -c "import numpy as np, capa0_dominio as C; d=np.load('capa0_dominio.npz',allow_pickle=True); print(np.round(C.capsule_distance(d['pyramid_apex']),3))"
> ```
> Causa y barrido de la variable dominante en
> `08_gemelo_digital/_auditoria/2026-08-09_seno_apices_particion.md` (Tarea 1) y
> `_auditoria/2026-08-10_centro_seno_particion.md` (Tarea 1a).
> **No usar la afirmación "interfaz médula/seno" en publicación mientras el defecto siga abierto.**
```

**Para #5 — `00_bitacora/BITACORA.md:998`.** No editar la entrada histórica. La corrección
corresponde a una entrada nueva de bitácora, en la línea de lo ya hecho con el sesgo polar mal
atribuido (ENTRADA 034 propuesta).

**Nota de método:** las cuatro notas dicen lo mismo y **ninguna borra el texto original**. El patrón es
el mismo que se usó con la atribución errónea a Beland (`auditoria_correspondencia_anatomica.md:142-152`):
preservar lo que se afirmó, marcar por qué no se sostiene.

---

## Resumen

| tarea | resultado |
|---|---|
| **0** | Las cifras externas (L=4.5, W=1.3, T=1.2 cm; V=3.68 cm³) **NO están en el repo**. Los aciertos de `3.68` son coordenadas vasculares en `renal_data_v1.json`. `CEA` no aparece. `ecogénica` sólo en mi informe previo, sin números. **No se usan ni se citan** |
| **1a** | Barrido cy ∈ [−34, −28] paso 0.5: 4 ápices corticales en −34, **0 desde −32.0**; PW_max = cy + 46; volumen 3.949 → 8.105 mL |
| **1b** | **Sí existe: `cy ∈ [−32.5, −27.8]`** (anchura 4.7 mm) cumple 0 ápices + PW en [12.6, 18.2]. `cy = −30.6` da PW = 15.4 exacto y cae cerca del centro. **El "rango razonable" de volumen es [PENDIENTE DE ANCLA]** |
| **1c** | `capa3_cobertura…py:87,88,145,167` es paramétrico, **seguro**. **`capa4_colector_alto.py:449` hardcodea la fórmula** en el informe generado: computaría con el seno nuevo y declararía el viejo. El texto ya está versionado en `auditoria_capa4_calicial.md:19` y `BITACORA.md:998`. Corrección propuesta (f-string) |
| **2** | Especificación completa. Reemplaza `capa0_dominio.py:356-382`; entradas `coords`, `medulla_mask`, `apex`; sin escape dentro de la médula. `pyramid_r_base` y `cone_half_angle_deg` pierden su único consumo interno (`:374`) y tienen **0 consumidores externos** → retirables. **0.0000 % de médula huérfana (82 297/82 297) con el seno ACTUAL**, frente a 47.050 % |
| **3** | 5 decisiones bloqueadas (B1–B5). **Ninguna bloquea la cascada** si se acepta declararlas supuesto. B3 se relaja: el barrido muestra un intervalo de 4.7 mm donde todo cumple lo anclable; Emamian elegiría un punto, su ausencia no impide elegir alguno. Regla propuesta: `cy = −30.6` como [SUPUESTO DECLARADO, calibrado contra Glodny PW] |
| **4** | 4 afirmaciones localizadas (+1 en bitácora), dos de ellas ya versionadas. El test comprueba `val ≈ 1.0` (sobre la pared), **no** medularidad. Notas propuestas para las cuatro, ninguna borra el texto original |

### Lo que este informe cambia

1. **B3 deja de ser bloqueante absoluto.** Hay un intervalo admisible de 4.7 mm, no un punto. La
   ausencia de Emamian impide *anclar*, no *decidir*.
2. **Un defecto nuevo, independiente del seno:** `capa4_colector_alto.py:449` imprime una fórmula
   literal. Hoy es correcta por coincidencia; cualquier cambio de seno la vuelve falsa en silencio.
   **Debe corregirse aunque no se toque `CENTRO_SENO`.**
3. **La partición queda especificada y numéricamente sostenida sobre la geometría actual.** Es la
   única decisión de fondo lista para tomarse sin fuentes nuevas, y la única que **reduce**
   parámetros libres: elimina `CONE_HALF_ANGLE_DEG`, la última constante [SIN DECLARAR] del archivo.

**Sin cambios en disco fuera de este archivo.** No se ejecutó `capa0_dominio.py`. No se hizo
`git add`, `git commit`, `git stash` ni `git restore`. Ningún `.py`, ningún `.npz` y ningún documento
existente fue modificado.
