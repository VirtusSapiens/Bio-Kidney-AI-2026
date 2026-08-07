# Auditoría — Tres reservas metodológicas sobre el intervalo de `CENTRO_SENO`

**Programa:** Bio-Kidney AI 2026 · **Capa:** 0 · **Fecha:** 2026-08-12

> **NATURALEZA: SOLO LECTURA / DIAGNÓSTICO.**
> No se modificó ningún `.py`, ningún `.npz`, ningún parámetro ni ninguna lógica. No se ejecutó
> `capa0_dominio.py`. Cálculo en memoria sobre `capa0_dominio.npz` y funciones puras importadas de
> `capa0_dominio` (guard `if __name__ == "__main__":` en `:610`).
>
> **Este informe REVISA CRÍTICAMENTE conclusiones propias** de `_auditoria/2026-08-09`,
> `_auditoria/2026-08-10` y de las ediciones a ENTRADA 033. Las correcciones de la Reserva 3b son
> propuestas de texto: **no se han aplicado**.
>
> **El resultado de la cota inferior violada (3.9599 mL vs 5.70/4.15 cm³) NO se toca:** sobrevive
> intacto a las tres reservas y se explica por qué en la Tarea Final.

---

## AVISO DE TRAZABILIDAD — «PW» no es un campo del repo

> **Léase antes que nada lo demás. Añadido 2026-08-13.**
>
> **`PW` no existe como código, como campo del `.npz` ni como variable de ninguna capa.**
>
> ```bash
> grep -rn "parenchymal\|PW\b\|ancho_parenq\|parenq_width" --include=*.py . --exclude-dir=.venv
> ```
> → **vacío. NO ENCONTRADO.**
>
> Es una **construcción de auditoría** introducida por mí el **2026-08-09**
> (`_auditoria/2026-08-09_seno_apices_particion.md`, Tarea 2a) y reutilizada en los informes del
> 08-10, 08-11 y 08-12. No lo calcula `capa0_dominio.py`, no se persiste y ninguna capa lo consume.
>
> ### Definición operativa exacta
>
> Sobre una **malla uniforme de 800 × 800 en el plano (X, Z)**, para cada celda donde coexisten la
> cápsula medial y la pared del seno:
>
> - `y_capsula(X,Z) = −B_SEMI · sqrt(1 − (X/A_SEMI)² − (Z/C_SEMI)²)` — rama **medial** (−Y) del
>   elipsoide principal 55/30/18;
> - `y_seno(X,Z) = CENTRO_SENO[1] + SEMIEJES_SENO[1] · sqrt(1 − (X/sx)² − (Z/sz)²)` — cara **+Y** del
>   elipsoide del seno;
> - `PW(X,Z) = y_seno − y_capsula`, definido sólo donde ambas existen y `y_seno > y_capsula`.
>
> El **«PW máximo»** es `max PW(X,Z)`, que se alcanza en `X = 0, Z = 0` (eje hiliar).
>
> ### Tres propiedades que hay que tener delante
>
> 1. **Es una distancia en Y, no una distancia perpendicular a la cápsula.** Coinciden sólo en el eje
>    hiliar, donde la normal de la cápsula es paralela a Y. Fuera de ahí, **no**.
> 2. **Los percentiles dependen del muestreo** (malla uniforme en (X, Z), sin ponderar por área ni
>    por volumen). Cambiar la ponderación mueve mediana y percentiles; **sólo el máximo es
>    invariante**.
> 3. **Es un valor sobre una geometría determinista**, sin dispersión. No es un estadístico
>    poblacional.
>
> ### NO CONFUNDIR con `depth_cortical_mm`
>
> | | `depth_cortical_mm` | «PW» (esta construcción) |
> |---|---|---|
> | ¿existe en el repo? | **Sí** — campo del `.npz`, `capa0_dominio.py:438` | **No** — sólo en informes de `_auditoria/` |
> | qué mide | distancia **perpendicular** de un punto del parénquima a la cápsula externa | separación **en Y** entre cápsula medial y pared del seno |
> | dónde está definido | en **todo** el parénquima, incluidos los polos | **sólo** donde el seno está detrás de la cápsula |
> | cómo se calcula | `capsule_distance()` (`:221`) vía `_nearest_point_ellipsoid` (`:140`) | malla analítica en (X, Z), fuera del código de producción |
> | contra qué se compara | `GROSOR_CORTICAL_MM = 6.6` (Glodny CW) | PW de Glodny — **con la reserva de estadístico de §1b–1c** |
>
> **Comparar percentiles de `depth_cortical_mm` contra el PW de Glodny es un error de magnitud**:
> `depth_cortical_mm` cubre los polos, donde no hay seno detrás y por tanto no hay PW que medir. Los
> dos campos coinciden en unidades (mm) y en nada más.

---

## RESERVA 1 — ¿Es mi PW la misma magnitud que la de Glodny?

### (a) Cómo definí el "PW máximo (eje hiliar)"

**Primer hecho que hay que dejar sentado: no hay ningún cálculo de PW en el repo.**

```bash
grep -rn "parenchymal\|PW\b\|ancho_parenq\|parenq_width" --include=*.py . --exclude-dir=.venv
```
→ **vacío. NO ENCONTRADO.**

`capa0_dominio.py` no calcula PW, no lo persiste y no lo compara con nada. **El PW es una construcción
mía**, introducida en `_auditoria/2026-08-09` (Tarea 2a) y reutilizada después. Su definición vive
sólo en el comando que la produce:

```bash
.venv/bin/python -c "import numpy as np
a,b,c=55.,30.,18.; sx,sy,sz=22.,16.,11.; cy=-34.
N=800; X=np.linspace(-a,a,N); Z=np.linspace(-c,c,N); XX,ZZ=np.meshgrid(X,Z,indexing='ij')
r=1-(XX/a)**2-(ZZ/c)**2; okc=r>0; ycap=np.where(okc,-b*np.sqrt(np.clip(r,0,None)),np.nan)
rs=1-(XX/sx)**2-(ZZ/sz)**2; oks=rs>0; ysen=np.where(oks,cy+sy*np.sqrt(np.clip(rs,0,None)),np.nan)
m=okc&oks&(ysen>ycap); w=(ysen-ycap)[m]
print('max %.3f  mediana %.3f'%(w.max(),np.median(w)))"
```

En palabras: sobre una malla uniforme en el plano (X, Z), para cada celda donde coexisten cápsula
medial y pared del seno, la **diferencia de coordenada Y** entre `ysen` y `ycap`. El "PW máximo" es
el máximo de ese campo, que se alcanza en X = 0, Z = 0.

**Tres propiedades de esa definición que conviene tener delante:**

1. Es una **distancia en Y**, no una distancia perpendicular a la cápsula. En X = 0, Z = 0 coinciden
   (la normal de la cápsula ahí es paralela a Y), pero fuera del eje hiliar **no**.
2. El estadístico depende del **muestreo**: malla uniforme en (X, Z), no ponderada por área de
   superficie ni por volumen. Cambiar la ponderación cambia la mediana, no el máximo.
3. Es un valor sobre **una** geometría determinista, sin dispersión.

### (b) ¿Son el mismo estadístico?

**NO.**

Glodny reporta **media ± desviación típica sobre n = 2068 riñones**: un estadístico poblacional que
resume la variabilidad entre individuos. Yo reporto un **máximo sobre la superficie de una única
geometría**: un extremo espacial dentro de un solo objeto, sin dispersión de ningún tipo.

Comparar un máximo espacial contra una media poblacional **es un error de magnitud de la misma clase
que el que corrigió el método B** (ENTRADA 032 §2): allí se comparaba una cuerda radial contra un
espesor perpendicular; aquí se compara un extremo espacial contra una media poblacional. En ambos
casos las unidades coinciden (mm) y la comparación parece legítima, y no lo es.

**Lo escribí como si lo fuera.** En `_auditoria/2026-08-09` (Tarea 2b) presenté "PW_max = 11.999 mm
frente a Glodny 15.4 ± 2.8" y lo llamé "el análogo más directo a una medida de corte axial único".
**Esa justificación no basta**: que ambos sean "una medida en un corte" no los hace el mismo
estadístico si uno es el extremo del campo y el otro una media entre sujetos.

### (c) ¿Cuál es la medida del gemelo comparable con la de Glodny?

Consultada la fuente primaria:

> Glodny B, Unterholzner V, Taferner B, Hofmann KJ, Rehder P, Strasak A, Petersen J.
> *"Normal kidney size and its influencing factors — a 64-slice MDCT study of 1.040 asymptomatic
> patients."* **BMC Urology. 2009;9:19.** DOI 10.1186/1471-2490-9-19.
> Texto completo: `https://pmc.ncbi.nlm.nih.gov/articles/PMC2813848/`

**Lo que la fuente sí especifica:**

- Los parámetros medidos incluyen `width of the parenchyma (PW) and the cortex (CW) in the arterial phase`.
- El pie de la Figura 1 sitúa la medida en un **corte axial**: `Axial 0.625 mm collimated slice of the kidney in an arterial phase, with the strongly contrasted kidney cortex`, mostrando `Cortical width (CW), and parenchymal width (PW)`.
- Control de calidad: `The measurements were performed twice in a random sample of 50 data sets`.
- Valores: **PW derecho 15.4 ± 2.8 mm**, **PW izquierdo 15.9 ± 2.7 mm**; **CW 6.6 ± 1.9 mm en ambos lados**. Longitud polo-polo 108.5 ± 12.2 (dcho) / 111.3 ± 12.6 mm (izdo).

**Lo que la fuente NO especifica:**

| pregunta | respuesta |
|---|---|
| ¿A qué **nivel** del riñón se mide el PW (polo superior, tercio medio, hilio, polo inferior)? | **[PENDIENTE DE ANCLA]** — el texto de Métodos no lo indica; el pie de figura dice sólo "axial" |
| ¿**Cuántas** medidas por riñón? | **[PENDIENTE DE ANCLA]** — no consta. El "twice in a random sample of 50" es reproducibilidad inter-observador, no el protocolo de rutina |
| ¿Es una localización estandarizada o un promedio de varias? | **[PENDIENTE DE ANCLA]** |
| ¿Se mide en el punto de máximo espesor, o en uno anatómicamente definido? | **[PENDIENTE DE ANCLA]** |

**Corrijo también la premisa del encargo.** La reserva dice *"Glodny reporta [...] UNA medida por
riñón en localización estandarizada"*. **La fuente no afirma eso.** No consta ni el número de medidas
ni la estandarización de la localización. Que Glodny obtenga una media ± s.d. sobre 2068 riñones no
implica que haya una sola medida por riñón: podría ser el promedio de varias, y el paper no lo dice.

**Conclusión de (c): no puedo determinar cuál es la medida del gemelo comparable con la de Glodny.**
Sin saber a qué nivel se mide el PW, no sé si el análogo correcto es el máximo (12.000), la mediana
(6.855), la media (6.590) o el valor en una localización concreta que aún no está definida.
**[PENDIENTE DE ANCLA].** No lo supongo.

### (d) Distribución completa del ancho parenquimatoso del gemelo

Para que la comparación pueda hacerse con el estadístico correcto cuando se sepa cuál es:

```bash
.venv/bin/python -c "import numpy as np
a,b,c=55.,30.,18.; sx,sy,sz=22.,16.,11.; cy=-34.
N=800; X=np.linspace(-a,a,N); Z=np.linspace(-c,c,N); XX,ZZ=np.meshgrid(X,Z,indexing='ij')
r=1-(XX/a)**2-(ZZ/c)**2; okc=r>0; ycap=np.where(okc,-b*np.sqrt(np.clip(r,0,None)),np.nan)
rs=1-(XX/sx)**2-(ZZ/sz)**2; oks=rs>0; ysen=np.where(oks,cy+sy*np.sqrt(np.clip(rs,0,None)),np.nan)
m=okc&oks&(ysen>ycap); w=(ysen-ycap)[m]
print('n=%d'%len(w)); [print('p%-3d %8.3f'%(q,np.percentile(w,q))) for q in [5,25,50,75,95]]
print('min %.3f max %.3f media %.3f sd %.3f'%(w.min(),w.max(),w.mean(),w.std()))"
```

| estadístico | valor (mm) |
|---|---|
| p5 | 0.846 |
| p25 | 3.781 |
| **p50 (mediana)** | **6.855** |
| p75 | 9.553 |
| p95 | 11.528 |
| mínimo | 0.001 |
| **máximo** | **12.000** |
| media ± sd espacial | 6.590 ± 3.400 |

n = 96 760 celdas, malla uniforme en (X, Z), `cy = −34`, `sy = 16`.

**Advertencia sobre la ponderación:** estos percentiles son sobre malla uniforme en (X, Z), **no**
ponderados por área de superficie. Una ponderación distinta desplaza mediana y percentiles; **sólo el
máximo es invariante**. Cuando se determine el estadístico de Glodny habrá que redefinir también la
ponderación del gemelo.

**Lo que sobrevive a la reserva y lo que no:**

- **Sobrevive la DIRECCIÓN.** Sea cual sea el estadístico elegido —máximo 12.000, mediana 6.855 o
  media 6.590— **todos caen por debajo de 12.6 mm**, el límite inferior de Glodny a −1 s.d. El
  parénquima del gemelo es más delgado que la referencia bajo cualquier lectura razonable. El p95
  (11.528) tampoco llega. **Esa conclusión es robusta.**
- **NO sobrevive la MAGNITUD del déficit.** Con el máximo el déficit es 3.4 mm; con la mediana, 8.5 mm.
  Y la magnitud es exactamente lo que fija el intervalo de `cy`. **El intervalo hereda una elección de
  estadístico que no está justificada.**

---

## RESERVA 2 — El intervalo es un corte 1D de un espacio sin anclar

### (a) La identidad `PW_max = |−B_SEMI − (cy + sy)|`

**Se cumple.** Contrastada contra la malla en cinco configuraciones:

| cy | sy | PW por malla | identidad | diferencia |
|---|---|---|---|---|
| −34.000 | 16.00 | 11.9998 | 12.0000 | 1.6 × 10⁻⁴ |
| −30.600 | 16.00 | 15.3998 | 15.4000 | 1.6 × 10⁻⁴ |
| −34.000 | 19.40 | 15.3998 | 15.4000 | 1.8 × 10⁻⁴ |
| −28.000 | 16.00 | 17.9998 | 18.0000 | 1.6 × 10⁻⁴ |
| −31.372 | 16.00 | 14.6278 | 14.6280 | 1.6 × 10⁻⁴ |

La diferencia residual es la resolución de la malla (110 mm / 800 celdas ≈ 0.1375 mm de paso, con
interpolación implícita en el máximo). **La identidad es exacta analíticamente:** en X = 0, Z = 0,
`ycap = −B_SEMI` y `ysen = cy + sy`.

**Consecuencia estructural: `PW_max` depende ÚNICAMENTE de la suma `cy + sy`.** No de cada uno por
separado. La condición PW no es una restricción sobre `cy`: es una restricción sobre una recta
diagonal del plano (cy, sy).

### (b) Reformulación sobre la suma

Sea **s = cy + sy**. Entonces `PW_max = s + B_SEMI = s + 30`.

| condición | restricción sobre `s` |
|---|---|
| PW_max = 15.4 (media de Glodny, riñón derecho) | **s = −14.6** |
| PW_max ∈ [12.6, 18.2] (±1 s.d.) | **s ∈ [−17.4, −11.8]** |

**El intervalo `cy ∈ [−31.37, −27.80]` que publiqué en ENTRADA 033 es válido SOLO para `sy = 16`.**
Y `sy = 16` está tan sin anclar como `cy`: ambos son los puestos 1 y 2 de los 10 pendientes en
`09_paper_vascular/auditoria_correspondencia_anatomica.md:128-129`. **La reserva es correcta y el
intervalo, tal como lo escribí, está mal planteado.**

### (c) Región 2D bajo las tres condiciones

Las otras dos condiciones **sí** dependen de `sy` por separado, no sólo de la suma:

- **Ápices:** la posición del ápice es `cy + uy·sy` con `uy` distinto para cada pirámide (0.741 a
  0.915). Sólo el ápice ecuatorial (uy → 1) dependería de la suma; los polares no.
- **Volumen excavado:** depende del solapamiento de dos elipsoides, función de `cy`, `sy`, `sx`, `sz`.

Región admisible bajo las tres, con `sx = 22` y `sz = 11` fijos:

```bash
.venv/bin/python -c "import numpy as np, capa0_dominio as C
A=C.SEMIEJES; sx,sz=22.,11.; b=30.
rng=np.random.default_rng(2026); n=300000
u=rng.normal(size=(n,3)); u/=np.linalg.norm(u,axis=1,keepdims=True)
P=(u*rng.uniform(0,1,(n,1))**(1/3))*A; VM=4/3*np.pi*np.prod(A)
Rq=1-(P[:,0]/sx)**2-(P[:,2]/sz)**2; Yq=P[Rq>0,1]; Sq=np.sqrt(Rq[Rq>0])
def vol(cy,sy): return VM*(np.abs(Yq-cy)<sy*Sq).sum()/n/1000
XF=np.array([((i+0.5)/10-0.5)*2.0 for i in range(10)]); ZR=np.array([1.0 if i%2==0 else -1.0 for i in range(10)])
UX=XF*0.60; UZ=ZR*0.40; UY=np.sqrt(np.clip(1-UX**2-UZ**2,0,None))
def nap(cy,sy): return int((C.capsule_distance(np.stack([UX*sx, cy+UY*sy, UZ*sz],1))<6.6).sum())
for sy in [12,14,16,18,20,22]:
    ok=[cy for cy in np.arange(-36,-24,0.1) if 12.6<=abs(-b-(cy+sy))<=18.2 and nap(cy,sy)==0 and vol(cy,sy)>5.70]
    Vell=4/3*np.pi*sx*sy*sz/1000
    lo,hi=min(ok),max(ok)
    print('sy=%2.0f cy in [%7.2f,%7.2f] anchura %.2f | s in [%6.2f,%6.2f] | dentro %.1f%%-%.1f%%'%(sy,lo,hi,hi-lo,lo+sy,hi+sy,100*vol(lo,sy)/Vell,100*vol(hi,sy)/Vell))"
```

| `sy` | `cy` admisible | anchura | `s = cy+sy` | % del elipsoide DENTRO del parénquima |
|---|---|---|---|---|
| 12 | [−28.50, −24.10] | 4.40 | [−16.50, −12.10] | 47.1 % → 74.4 % |
| 14 | [−29.90, −25.90] | 4.00 | [−15.90, −11.90] | 40.3 % → 61.4 % |
| **16** | **[−31.30, −27.90]** | 3.40 | [−15.30, −11.90] | 35.3 % → 50.6 % |
| 18 | [−32.70, −29.90] | 2.80 | [−14.70, −11.90] | 31.5 % → 42.4 % |
| 20 | [−34.20, −31.90] | 2.30 | [−14.20, −11.90] | 28.3 % → 36.0 % |
| 22 | [−35.70, −33.90] | 1.80 | [−13.70, −11.90] | 25.7 % → 31.1 % |

*(Verificación cruzada con un segundo cálculo, paso 0.05 y otra semilla: bordes dentro de ±0.05 mm.)*

**La región admisible es una banda diagonal**, no un intervalo. Se estrecha al crecer `sy`. La fila
`sy = 16` reproduce mi intervalo previo — **como un corte, no como el resultado**.

**Caveat que la propia reserva induce y hay que declarar:** `sx = 22` y `sz = 11` **también están sin
anclar**. La región verdadera es **4D en (cy, sy, sx, sz)**; lo de arriba es una rebanada 2D de ella.
La crítica de la Reserva 2 se aplica recursivamente a mi propia respuesta. Se declara, no se disimula.

### (d) ¿Alguna combinación deja el elipsoide mayoritariamente DENTRO?

**Sí, en la esquina de `sy` pequeño y `cy` interior.**

Estado actual, para referencia:

```
cy=-34, sy=16 -> Vexc 3.99 mL | elipsoide 16.22 mL -> DENTRO 24.6 %  FUERA 75.4 %
```

Sobre la región admisible la fracción interior recorre **25.7 % → 74.4 %**. Supera el 50 % en:

- `sy = 12`, todo el intervalo admisible (47.1 % → 74.4 %; cruza el 50 % en torno a `cy ≈ −27.9`)
- `sy = 14`, la mitad superior del intervalo (hasta 61.4 %)
- `sy = 16`, sólo el extremo superior (50.6 % en `cy = −27.90`)

**Nunca en `sy ≥ 18`.** El máximo alcanzable en la región es **74.4 %** (`sy = 12`, `cy = −24.10`).

**Lectura:** un seno **más pequeño y más interior** deja mucha mayor fracción de su elipsoide
haciendo trabajo real de excavación, y sigue cumpliendo las tres condiciones. La configuración actual
—grande y muy exterior— es la peor de la familia en ese sentido: tres cuartas partes del elipsoide no
excavan nada.

**Lo que esto NO autoriza a concluir:** que `sy = 12` sea "mejor". No hay ancla que diga que el
elipsoide deba estar mayoritariamente dentro; es una propiedad de parsimonia del modelo, no un hecho
anatómico. **[PENDIENTE DE ANCLA].**

---

## RESERVA 3 — Independencia de las condiciones

### (a) Monotonía

```bash
.venv/bin/python -c "import numpy as np, capa0_dominio as C
A=C.SEMIEJES; sx,sz=22.,11.; b=30.; sy=16.
rng=np.random.default_rng(2026); n=300000
u=rng.normal(size=(n,3)); u/=np.linalg.norm(u,axis=1,keepdims=True)
P=(u*rng.uniform(0,1,(n,1))**(1/3))*A; VM=4/3*np.pi*np.prod(A)
Rq=1-(P[:,0]/sx)**2-(P[:,2]/sz)**2; Yq=P[Rq>0,1]; Sq=np.sqrt(Rq[Rq>0])
XF=np.array([((i+0.5)/10-0.5)*2.0 for i in range(10)]); ZR=np.array([1.0 if i%2==0 else -1.0 for i in range(10)])
UX=XF*0.60; UZ=ZR*0.40; UY=np.sqrt(np.clip(1-UX**2-UZ**2,0,None))
prev=None
for cy in np.arange(-34,-27.9,0.5):
    v=(float(C.capsule_distance(np.stack([UX*sx, cy+UY*sy, UZ*sz],1)).min()), abs(-b-(cy+sy)), VM*(np.abs(Yq-cy)<sy*Sq).sum()/n/1000)
    print('%7.1f %14.3f %9.2f %10.3f'%(cy,*v))
    if prev: assert all(x>y for x,y in zip(v,prev)),'NO monotona'
    prev=v
print('-> las tres estrictamente CRECIENTES en cy')"
```

| `cy` | `depth` mín de ápice | PW_max | V excavado (mL) |
|---|---|---|---|
| −34.0 | 5.295 | 12.00 | 3.994 |
| −33.0 | 6.045 | 13.00 | 4.622 |
| −32.0 | 6.757 | 14.00 | 5.275 |
| −31.0 | 7.428 | 15.00 | 5.939 |
| −30.0 | 8.058 | 16.00 | 6.635 |
| −29.0 | 8.645 | 17.00 | 7.371 |
| −28.0 | 9.190 | 18.00 | 8.128 |

**Las tres son estrictamente crecientes en `cy`.** El aserto de monotonía pasa en los 13 pasos. Y es
estructural, no numérico: las tres miden, con métricas distintas, **cuánto se adentra el seno en el
parénquima**.

### (b) La formulación correcta — y la corrección que debo a ENTRADA 033

**Tres condiciones monótonas en la misma dirección no son tres anclas independientes.** Son **tres
lecturas de un mismo eje**. Cada una aporta un umbral sobre la misma variable latente; su acuerdo no
es corroboración mutua, porque no podían discrepar en dirección.

Lo que sí aportan por separado: **umbrales distintos**, y por tanto **bordes distintos** del
intervalo. Eso es real y es lo que hay que decir.

**Frase incorrecta, localizada:**

- `00_bitacora/BITACORA.md:1817` — `> 5. La corrección de `CENTRO_SENO` deja de ser opcional: hoy corrige **cuatro** defectos independientes — ápices corticales, PW frente a Glodny, espesor medular, y ahora la cota de volumen de Caglar.`

**Corrección propuesta (texto, SIN APLICAR)** para esa línea:

> `> 5. La corrección de CENTRO_SENO deja de ser opcional. **Cuatro síntomas del mismo defecto** —ápices corticales, PW frente a Glodny, espesor medular y la cota de volumen de Caglar— **apuntan en la misma dirección: el seno está demasiado afuera**. No son cuatro anclas independientes: las cuatro magnitudes son **estrictamente monótonas crecientes en `CENTRO_SENO[1]`** (barrido en `_auditoria/2026-08-12`), de modo que su acuerdo en dirección era necesario y no constituye corroboración mutua. Lo que sí aportan por separado son **umbrales distintos**, y por tanto los dos bordes del intervalo admisible: Caglar fija el borde inferior, Glodny el superior, y la condición de ápices es **redundante** (nunca vinculante). Región admisible en `_auditoria/2026-08-12`, §Reserva 2c — es **2D en (cy, sy)**, no un intervalo en `cy`.`

**Segunda frase a matizar, no a corregir:**

- `04_literatura/anclas_seno_renal.md:141` — `Dos vías independientes apuntan al mismo parámetro. **Eso es lo que se registra: una convergencia cualitativa.**`
- `04_literatura/anclas_seno_renal.md:158` — `**Uso legítimo:** citarlo como *"dos vías independientes convergen en señalar `CENTRO_SENO` y no `SEMIEJES_SENO`"*.`

Aquí **"independientes" sí es defendible**: una vía es literatura (cruce Caglar × Zhang), la otra es
geometría (barridos de sensibilidad). Son tipos de evidencia distintos. **Pero conviene añadir la
salvedad** de que el lado geométrico es una familia monótona:

> **Nota propuesta a añadir tras `anclas_seno_renal.md:141`:**
> *"Salvedad (2026-08-12): «independientes» se refiere al tipo de evidencia —literatura frente a
> geometría—, no a las condiciones geométricas entre sí. Las tres condiciones geométricas (ápices,
> PW, volumen) son estrictamente monótonas en `CENTRO_SENO[1]` y por tanto no pueden discrepar en
> dirección; su acuerdo no es corroboración mutua. Ver `_auditoria/2026-08-12`, Reserva 3."*

### (c) Qué fija cada borde

Con `sy = 16`:

| condición | umbral sobre `cy` | papel |
|---|---|---|
| PW ≥ 12.6 (Glodny −1 s.d.) | cy ≥ −33.40 | **redundante** (el menos restrictivo) |
| 0 ápices corticales | cy ≥ −32.50 | **redundante** |
| **V excavado > 5.70 (Caglar)** | **cy ≥ −31.30** | **VINCULANTE — borde inferior** |
| **PW ≤ 18.2 (Glodny +1 s.d.)** | **cy ≤ −27.80** | **VINCULANTE — borde superior, único** |

**Resultados:**

1. **Borde inferior: Caglar.** Es la única condición que lo fija. Sin Caglar, el borde lo pondría la
   condición de ápices en −32.50, y el intervalo sería 1.2 mm más ancho.
2. **Borde superior: Glodny.** Es la **única** condición con cota superior. Ni el volumen ni los
   ápices penalizan mover el seno más adentro; los dos siguen mejorando indefinidamente.
3. **La condición de ápices corticales es redundante en los dos bordes.** No fija nada.

**Consecuencia para ENTRADA 033.** El defecto de los cuatro ápices corticales —que descubrí en
`_auditoria/2026-08-08` y traté como hallazgo principal— **no aporta ninguna restricción al
intervalo**. Se corrige solo en cuanto se satisface Caglar. Sigue siendo un defecto real y su
diagnóstico sigue en pie; **pero no es un criterio de calibración**, y presentarlo como uno de los
cuatro pilares del intervalo era inexacto.

---

## TAREA FINAL — Estado de la decisión

### ¿Queda región admisible no vacía?

**Sí.** Ninguna de las tres reservas la vacía. Sobre `(cy, sy)` con `sx = 22`, `sz = 11`, es una banda
diagonal de anchura 1.8–4.4 mm en `cy` según `sy`, existente al menos para `sy ∈ [12, 22]`.

### ¿Sobre qué parámetros exactamente?

**Sobre `(CENTRO_SENO[1], SEMIEJES_SENO[1])` = `(cy, sy)`** — y con la salvedad de que `sx` y `sz`
también están sin anclar, de modo que el espacio real es **4D**. Lo reportado es una rebanada.

**Lo que las anclas actuales restringen realmente:**

| ancla | restringe |
|---|---|
| Glodny (PW) | la **suma** `cy + sy ∈ [−17.4, −11.8]` — una banda diagonal, no un punto |
| Caglar (volumen) | una región curva en `(cy, sy, sx, sz)`, cota inferior |
| ápices corticales | **nada** (redundante) |

### ¿Qué haría falta para reducirla a un punto?

**No basta con Emamian, y no basta con un volumen.**

| dato | qué fijaría | ¿reduce a punto? |
|---|---|---|
| **Volumen de la CEA** (Emamian) | una ecuación en 4 incógnitas | **No** — reduce 4D a 3D |
| **L, W, T de la CEA** (Emamian) | `sx`, `sy`, `sz` (3 ecuaciones) | **Casi** — quedaría `cy` libre, y Glodny lo fijaría por `s = cy + sy`. **Sí, en combinación** |
| **Posición de la CEA** dentro del contorno renal | `cy` directamente | **Sí, en combinación con las dimensiones** |

**Respuesta directa: lo que reduciría la región a un punto son las TRES DIMENSIONES de la CEA, no su
volumen.** Con `sx`, `sy`, `sz` anclados, `cy` queda determinado por la condición de Glodny sobre la
suma. Si Emamian aporta además la posición, el sistema queda sobredeterminado y eso sería una
comprobación de consistencia, no un grado de libertad.

En `_auditoria/2026-08-10` (Tarea 3, B3) escribí que *"el volumen de Emamian es exactamente el
desempate"*. **Eso era inexacto por la misma razón que la Reserva 2:** el volumen es **una** ecuación
sobre un espacio de **cuatro** incógnitas. El desempate son las dimensiones.

### Qué queda firme y qué no

| resultado | estado tras las tres reservas |
|---|---|
| **Cota inferior de volumen violada** (3.9599 mL excavados < 5.70 y < 4.15 cm³ de grasa sola) | **FIRME.** No depende de PW, ni de la elección de estadístico, ni de `sy`: compara dos volúmenes con el mismo criterio de contención. Es la única conclusión de esta línea que sobrevive intacta a las tres reservas |
| **El parénquima del gemelo es más delgado que Glodny** | **FIRME EN DIRECCIÓN.** Máximo (12.000), p95 (11.528), mediana (6.855) y media (6.590) caen todos bajo 12.6. La **magnitud** del déficit, no |
| **Los 4 ápices corticales son un defecto** | **FIRME como defecto.** Pero **redundante como criterio**: no fija ningún borde |
| **`CENTRO_SENO` es la variable dominante** | **FIRME** — sostenido por el barrido de `_auditoria/2026-08-10` §1c, que no depende de estas reservas |
| **Intervalo `cy ∈ [−31.37, −27.80]`** | **NO FIRME.** Es un corte en `sy = 16` de una región 2D (realmente 4D), y su borde superior depende de un estadístico de PW no determinado. **No debe commitearse** |
| *"cuatro defectos independientes"* | **INCORRECTO.** Son cuatro síntomas monótonos de un eje. Corrección propuesta en Reserva 3b |

### Recomendación de alcance para el commit

**Commitear:** el resultado de la cota violada (ya aceptado), y la corrección de la frase
"cuatro defectos independientes".

**No commitear todavía:** el intervalo numérico de `CENTRO_SENO[1]`, en ninguna de sus formas. Está
pendiente de (i) determinar el estadístico de PW de Glodny —**[PENDIENTE DE ANCLA]**— y (ii)
reformularse como región sobre `(cy, sy)` o sobre la suma `cy + sy`, nunca como intervalo en `cy`
solo.

**Sin cambios en disco fuera de este archivo.** No se ejecutó `capa0_dominio.py`. No se hizo
`git add`, `git commit`, `git stash` ni `git restore`. Ningún `.py`, ningún `.npz` y ningún documento
existente fue modificado; las correcciones de la Reserva 3b son texto propuesto.

---

**Fuente primaria consultada para la Reserva 1c:**
[Glodny B et al., *Normal kidney size and its influencing factors — a 64-slice MDCT study of 1.040 asymptomatic patients*, BMC Urology 2009;9:19 (PMC2813848)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2813848/)
