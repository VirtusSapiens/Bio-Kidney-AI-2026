# Auditoría — Médula huérfana y simulación del paquete A+C+E

**Programa:** Bio-Kidney AI 2026 · **Capa:** 0 · **Fecha:** 2026-08-08

> **NATURALEZA: SOLO LECTURA / DIAGNÓSTICO.**
> No se modificó ningún `.py`, ningún `.npz`, ningún parámetro ni ninguna lógica. No se ejecutó
> `capa0_dominio.py`. Todas las cifras se recomputan **en memoria** desde `capa0_dominio.npz`,
> `capa1_nefronas.npz` y `capa3c_colector.npz` mediante `np.load` y funciones puras importadas de
> `capa0_dominio` (el guard `if __name__ == "__main__":` de `capa0_dominio.py:610` impide que el
> import dispare `main()`).
>
> **LA SIMULACIÓN DE LA TAREA 2 NO ESTÁ APLICADA.** Es cálculo en memoria a partir de un script
> desechable en scratchpad; **no se escribió código de producción** y ningún `.py` del repo fue
> tocado. Los valores simulados no existen en disco fuera de este documento.
>
> Continúa `_auditoria/2026-08-07_reparametrizacion_piramides.md`, del que hereda la numeración de
> opciones A–F.

---

## Punto de partida — las tres precisiones aceptadas

1. **La Opción C no es alternativa a B: es la definición.** `GROSOR_CORTICAL_MM` está anclado a
   Glodny 2009 como *cortical width* = distancia **perpendicular cápsula → base de la pirámide**.
   "Base" y "punto a 6.6 mm perpendicular" son el mismo objeto. B aproxima con **2.331 mm** de error
   en los polos lo que C calcula exacto. **B queda subsumida por C y no vuelve a considerarse.**
2. **Hipótesis de trabajo:** el confinamiento de ápices a ±18 mm es **anatomía, no defecto** — en el
   riñón real las papilas convergen a un seno central y las pirámides polares tienen ejes largos y
   oblicuos. **`SEMIEJES_SENO` no se reabre.**
3. **Consecuencia:** `r_base = length·tan(22°)` (`capa0_dominio.py:352`) es el problema real. La
   Opción E tiene ancla propia — compuestas (polares) vs simples (mediopolares) son categorías
   **discretas** en Bonsib, no un gradiente proporcional a la longitud.

---

## TAREA 1 — Médula huérfana

### (a) Cuantificación

```bash
.venv/bin/python -c "import numpy as np; d=np.load('capa0_dominio.npz',allow_pickle=True); lab=d['region_label']; cx=lab=='cortex'; gen=lab=='medulla'; pir=np.char.startswith(lab,'piramide'); med=gen|pir; N=len(lab); print('cortex',int(cx.sum()),round(100*cx.mean(),3),'%'); print('medulla generica',int(gen.sum()),round(100*gen.mean(),3),'%'); print('piramide_XX',int(pir.sum()),round(100*pir.mean(),3),'%'); print('MEDULA total',int(med.sum()),round(100*med.mean(),3),'%'); print('HUERFANA = %d/%d = %.3f %% de la medula'%(gen.sum(),med.sum(),100*gen.sum()/med.sum()))"
```

| etiqueta | puntos | % del dominio |
|---|---|---|
| `cortex` | 117 703 | 58.852 % |
| `medulla` (genérica, sin pirámide) | **38 721** | 19.360 % |
| `piramide_XX` | 43 576 | 21.788 % |
| **médula total** | **82 297** | 41.148 % |

**Fracción de médula huérfana: 38 721 / 82 297 = 47.050 %.**

Casi la mitad de la médula del gemelo no pertenece a ningún cono piramidal.

### (b) Dónde está esa médula

```bash
.venv/bin/python -c "import numpy as np, capa0_dominio as C; d=np.load('capa0_dominio.npz',allow_pickle=True); lab=d['region_label']; co=d['coords'].astype(np.float64); dm=d['depth_cortical_mm'].astype(np.float64); gen=lab=='medulla'; pir=np.char.startswith(lab,'piramide'); [ (print(nm), [print('  %s min %8.2f p25 %7.2f mediana %7.2f p75 %7.2f max %7.2f |%s| mediana %6.2f'%(ax,S[:,i].min(),np.percentile(S[:,i],25),np.median(S[:,i]),np.percentile(S[:,i],75),S[:,i].max(),ax,np.median(np.abs(S[:,i])))) for i,ax in enumerate('XYZ')], print('  depth_mm min %.2f mediana %.2f max %.2f'%(D.min(),np.median(D),D.max())), print('  dist perp pared seno: mediana %.2f p10 %.2f max %.2f'%tuple(np.percentile(np.linalg.norm((S-C.CENTRO_SENO)-C._nearest_point_ellipsoid(S-C.CENTRO_SENO,C.SEMIEJES_SENO),axis=1),[50,10,100])[[0,1,2]]))) for nm,S,D in [('HUERFANA',co[gen],dm[gen]),('EN CONO',co[pir],dm[pir])]]"
```

| | HUÉRFANA | EN CONO |
|---|---|---|
| X: mediana / \|X\| mediana / max | 0.01 / **24.27** / 47.91 | 0.14 / **12.76** / 37.57 |
| Y: mediana / \|Y\| mediana | −4.08 / 10.92 | 2.75 / 6.41 |
| Z: mediana / \|Z\| mediana | −0.01 / **2.44** | 0.01 / **5.39** |
| `depth_mm`: min / mediana / max | 6.60 / 9.43 / 17.98 | 6.60 / 10.69 / 17.33 |
| dist. perp. a pared del seno: p10 / mediana / max | 5.66 / 24.48 / 42.93 | 11.67 / 23.80 / 38.30 |

**No es peri-sinusal.** Su distancia mediana a la pared del seno (24.48 mm) es incluso mayor que la de
la médula en cono (23.80 mm).

Descomposición:

```bash
.venv/bin/python -c "import numpy as np; d=np.load('capa0_dominio.npz',allow_pickle=True); lab=d['region_label']; co=d['coords'].astype(np.float64); gen=lab=='medulla'; G=co[gen]; pol=np.abs(G[:,0])>24.75; inter=(~pol)&(np.abs(G[:,2])<2.2); print('huerfana total',len(G)); print('  POLAR (|X|>24.75)      %6d  %5.1f %%'%(pol.sum(),100*pol.mean())); print('  INTER-HILERA (|Z|<2.2) %6d  %5.1f %%'%(inter.sum(),100*inter.mean())); print('  resto                  %6d  %5.1f %%'%((~pol&~inter).sum(),100*(~pol&~inter).mean()))"
```

| zona | puntos | % de la huérfana |
|---|---|---|
| **Polar** (\|X\| > 24.75, más allá de la base piramidal actual) | 18 883 | **48.8 %** |
| **Inter-hilera** (\|Z\| < 2.2, entre las dos filas Z = ±4.4) | 9 431 | 24.4 % |
| resto | 10 407 | 26.9 % |

Perfil radial — fracción de la médula **local** que queda huérfana:

```bash
.venv/bin/python -c "import numpy as np; d=np.load('capa0_dominio.npz',allow_pickle=True); lab=d['region_label']; co=d['coords'].astype(np.float64); gen=lab=='medulla'; pir=np.char.startswith(lab,'piramide'); [print('  |X| [%5.2f,%5.2f)  medula %6d  huerfana %6d = %5.1f %%'%(lo,hi,(mg:=gen&(np.abs(co[:,0])>=lo)&(np.abs(co[:,0])<hi)).sum()+(mp:=pir&(np.abs(co[:,0])>=lo)&(np.abs(co[:,0])<hi)).sum(),mg.sum(),100*mg.sum()/max(mg.sum()+mp.sum(),1))) for lo,hi in [(0,10),(10,20),(20,24.75),(24.75,30),(30,40),(40,55)]]"
```

| banda \|X\| | médula local | huérfana | % |
|---|---|---|---|
| [0, 10) | 25 713 | 8 861 | 34.5 % |
| [10, 20) | 23 791 | 6 752 | 28.4 % |
| [20, 24.75) | 9 561 | 4 225 | 44.2 % |
| [24.75, 30) | 9 157 | 5 790 | 63.2 % |
| [30, 40) | 11 438 | 10 456 | **91.4 %** |
| [40, 55) | 2 637 | 2 637 | **100.0 %** |

**Más allá de X = ±40 mm, el 100 % de la médula es huérfana.** No hay ni una pirámide que la reclame.

### (b-bis) HALLAZGO COLATERAL — cuatro ápices están en territorio córtex

```bash
.venv/bin/python -c "import numpy as np, capa0_dominio as C; d=np.load('capa0_dominio.npz',allow_pickle=True); da=C.capsule_distance(d['pyramid_apex']); print('depth_cortical_mm de los 10 apices:',np.round(da,3)); print('umbral',float(C.GROSOR_CORTICAL_MM)); print('apices con depth<6.6:',int((da<6.6).sum()),list(np.where(da<6.6)[0]))"
```

```
depth_cortical_mm de los 10 apices: [5.295 6.428 7.187 7.655 7.879 7.879 7.655 7.187 6.428 5.295]
umbral 6.6
apices con depth<6.6: 4 [0, 1, 8, 9]
```

**Los ápices k = 0, 1, 8, 9 están a 5.295 y 6.428 mm de la cápsula — por debajo del umbral 6.6 mm
de Capa 0.** Es decir: **cuatro papilas se encuentran en territorio etiquetado córtex.**
Anatómicamente la papila es el punto más profundo de la pirámide, en la interfaz médula/seno. Este
defecto **no aparecía en ninguna auditoría previa** y tiene consecuencia directa sobre la Tarea 2.

### (c) ¿Sostiene el repo la médula genérica?

Lectura anatómica propuesta: la médula **es** el conjunto de las pirámides; lo intercalado entre
ellas son **columnas de Bertin**, que son **córtex**. Bajo esa lectura, una etiqueta "médula que no
pertenece a ninguna pirámide" no debería existir.

**Búsqueda de justificación en el código:**

```bash
grep -n "medula generica\|generic\|n_med_generic" capa0_dominio.py
```

Sólo cuatro apariciones, **todas descriptivas, ninguna justificativa**:

- `capa0_dominio.py:26` — `region_label  (N,)  str       'cortex' | 'medulla' | 'piramide_XX'`
- `capa0_dominio.py:357` — `"""Asigna cada punto de la medula a una piramide (o a 'medula' generica).`
- `capa0_dominio.py:463` — `n_med_generic = n_medulla - n_in_pyr`
- `capa0_dominio.py:479-480` — `print(f"     - medula generica: {n_med_generic:>8d}  "` …

Es un **residuo aritmético** (`n_medulla − n_in_pyr`), no una categoría con intención declarada.

**Búsqueda en la bitácora y documentos:**

```bash
grep -rn "médula genérica\|medula generica\|columnas de Bertin\|Bertin" --include=*.md . --exclude-dir=.git
```

**NO ENCONTRADO** ningún pasaje que justifique la médula genérica como decisión de diseño. Las dos
únicas apariciones de "columnas de Bertin" son parte de la **definición de CW de Glodny**, y apuntan
en dirección contraria:

- `09_paper_vascular/auditoria_correspondencia_anatomica.md:47` — `**DEFINICIÓN DE CW (*cortical width*):** distancia **perpendicular desde la cápsula renal hasta la base de la pirámide medular** — es decir, **corteza sola** —, medida sobre **cortes axiales** del MDCT, **evitando las columnas de Bertin y el seno renal**.`
- `00_bitacora/BITACORA.md:1608` — `- La zona es la misma que Glodny excluye explícitamente en su protocolo (evitar columnas de Bertin y seno renal), de modo que la limitación cae donde la fuente tampoco mide.`

**Respuesta: el repo NO sostiene la médula genérica.** La fuente que el propio gemelo cita para su
umbral córtico-medular trata las columnas de Bertin como **córtex** y define la base de la pirámide
como el final del córtex. La categoría `'medulla'` sin pirámide **no tiene razón de diseño
localizable**: es lo que queda cuando los conos no cubren.

### (d) La fracción huérfana como criterio de aceptación

**Propuesta:** adoptar `fraccion_huerfana = n('medulla') / (n('medulla') + n(piramide_XX))` como
métrica de aceptación cuantitativa de cualquier reparametrización, calculada sobre el `.npz`
regenerado.

**Valor de referencia actual: 47.050 %.**

**Ancla para el umbral.** La literatura no da un número directo de "fracción de médula no piramidal"
— **NO ENCONTRADO**. Lo que sí da es la estructura: bajo Bonsib, la médula **es** el conjunto de
pirámides y lo intercalado es córtex (Bertin). El objetivo teórico es por tanto **0 %**, y todo
residuo es error de teselación del modelo cónico, no anatomía.

Como el modelo es de conos rectos sobre un elipsoide achatado, 0 % no es alcanzable. Un umbral
defendible se puede anclar en lo que la propia geometría permite: **la fracción huérfana debería ser
menor que la fracción de médula que queda fuera del alcance de cualquier cono por razones de
contorno**. Con eso, tres tramos:

| tramo | lectura |
|---|---|
| **< 10 %** | la teselación es buena; el residuo es esquina de contorno |
| **10–25 %** | aceptable con limitación declarada, indicando dónde queda el residuo |
| **> 25 %** | el modelo cónico no describe la médula; el `region_label` no es informativo |

**El estado actual (47.05 %) cae en el tercer tramo.** Esa es la lectura honesta: hoy la etiqueta
`piramide_XX` describe **poco más de la mitad** de la médula del gemelo.

**Salvedad:** este umbral es **[SUPUESTO DECLARADO]**, no anclado. Lo propongo como criterio
operativo, no como cifra con respaldo bibliográfico.

---

## TAREA 2 — Simulación del paquete A+C+E (NO APLICADA)

### Montaje

- **A:** `base_objetivo = (xf·k·a, 0.40·b, zrow·0.62·c)` con k ∈ {0.65, 0.75, 0.85}
- **C:** base real = punto donde `capsule_distance` = 6.6 mm marchando desde el ápice en la dirección
  ápice→objetivo
- **E:** `r_base` categórico: **15.00 mm** para polares (k ∈ {0,1,8,9}, compuestas) y **12.98 mm**
  para centrales (k ∈ {2..7}, simples), en lugar de `length·tan(22°)`

Ápices **sin cambio** respecto del código actual (hipótesis de trabajo 2). Asignación reproducida
según `capa0_dominio.py:369-380`; médula definida como `depth_cortical_mm ≥ 6.6`.

### Línea base (código actual)

| métrica | valor |
|---|---|
| huérfana | **47.05 %** |
| solape (≥2 conos) | 21.78 % de la médula |
| `length` | 37.12, 35.11, 33.62, 32.62, 32.12, 32.12, 32.62, 33.62, 35.11, 37.12 |
| `r_base` | 15.00, 14.19, 13.58, 13.18, 12.98, 12.98, 13.18, 13.58, 14.19, 15.00 |

### Resultado A+C+E, dirección ápice→objetivo

| k | ápices resueltos | base X (resueltos) | `length` | oblicuidad vs X | **huérfana** | solape |
|---|---|---|---|---|---|---|
| 0.65 | **6 / 10** | ±16.00 … ±3.34 | 28.8–29.0 | 71.1°–86.0° | **68.50 %** | 7.76 % |
| 0.75 | **6 / 10** | ±18.02 … ±3.83 | 28.8–29.1 | 66.9°–85.0° | **66.79 %** | 6.07 % |
| 0.85 | **6 / 10** | ±19.90 … ±4.32 | 28.8–29.3 | 62.9°–84.0° | **65.21 %** | 4.59 % |

En los tres casos, `depth` de las bases resueltas = **6.600 mm exacto** — la Opción C cumple su
definición por construcción.

**Dos resultados negativos, ambos firmes:**

1. **Los ápices k = 0, 1, 8, 9 no admiten solución.** Su profundidad (5.295 y 6.428 mm) ya es
   **menor** que 6.6, de modo que marchando hacia fuera no existe cruce con la UCM: el ápice está del
   lado cortical de la superficie que la base debería tocar. **No es un problema de dirección ni de
   k: es el hallazgo (b-bis).** Mientras esas cuatro papilas estén donde están, la Opción C es
   inaplicable a las cuatro pirámides polares — justo las que la literatura señala como compuestas.
2. **La huérfana EMPEORA** (47.05 % → 65–69 %). Causa: la Opción E fija `r_base` = 12.98 mm en las
   seis centrales, mientras el código actual les daba 12.98–13.58 derivado, y sobre todo las bases se
   acortan (`length` 32–34 → 28.8–29.3), de modo que los conos cubren menos volumen. El solape cae
   (21.78 % → 4.59 %), que es el efecto buscado, pero a costa de dejar más médula sin reclamar.

### Cota superior: A con dirección de máximo alcance polar

La dirección ápice→objetivo no es la que más lejos llega. Un barrido de dirección en el plano XY
muestra que desde los ápices centrales la UCM es alcanzable hasta **\|X\| = 42.14 mm**:

```bash
.venv/bin/python -c "import numpy as np, capa0_dominio as C; sx,sy,sz=C.SEMIEJES_SENO; cx,cy,cz=C.CENTRO_SENO; A=np.array([[cx+((i+0.5)/10-0.5)*2*0.60*sx, cy+np.sqrt(max(0,1-(((i+0.5)/10-0.5)*2*0.60)**2-0.16))*sy, cz+(1.0 if i%2==0 else -1.0)*0.40*sz] for i in range(10)]); th=np.linspace(-np.pi/2,np.pi/2,241); U=np.stack([np.sin(th),np.cos(th),np.zeros_like(th)],1); ts=np.linspace(0,70,220); [print('%2d apexX %7.2f depth %7.3f maxBaseX %s'%(i,A[i,0],float(C.capsule_distance(A[i:i+1])[0]),(lambda b: '%.2f'%b if b>-1e8 else 'SIN SOLUCION')(max([abs((A[i]+np.interp(0.0,[(v:=C.capsule_distance((A[i][None,None,:]+ts[None,:,None]*U[:,None,:]).reshape(-1,3)).reshape(241,220)[j]-6.6)[(k:=np.where(v<0)[0][0])],v[k-1]],[ts[k],ts[k-1]])*U[j])[0]) for j in range(241) if (vv:=C.capsule_distance((A[i][None,None,:]+ts[None,:,None]*U[:,None,:]).reshape(-1,3)).reshape(241,220)[j]-6.6)[0]>=0 and len(np.where(vv<0)[0])>0]+[-1e9])))) for i in range(10)]"
```

```
 0 apexX  -11.88 depth   5.295 maxBaseX SIN SOLUCION
 1 apexX   -9.24 depth   6.428 maxBaseX SIN SOLUCION
 2 apexX   -6.60 depth   7.187 maxBaseX 42.14
 …
 9 apexX   11.88 depth   5.295 maxBaseX SIN SOLUCION
```

**Corrección a mi lectura del informe anterior:** lo que topa las bases en ~19.9 mm no es la UCM,
sino la **dirección elegida** (ápice→objetivo tiene una componente +Y grande y sale por arriba antes
de avanzar en X). Con dirección optimizada, la UCM permite llegar a **42.14 mm = 76.6 % del semieje**.

Mejor caso resultante, con las cuatro polares heredando su base actual:

| variante | base X | `length` | oblicuidad vs X | **huérfana** | solape |
|---|---|---|---|---|---|
| A(máx polar)+C+**E** | ±42.14 (centrales), ±24.75/±19.25 (polares heredadas) | 35.1–45.7 | 25.5°–73.4° | **44.17 %** | 22.25 % |
| A(máx polar)+C, **sin E** (`r_base = L·tan22`) | ídem | ídem | ídem | **41.13 %** | 33.42 % |

**Aun en el mejor caso alcanzable con los ápices actuales, la huérfana sólo baja de 47.05 % a
44.17 % (con E) o 41.13 % (sin E).** Sigue en el tramo "> 25 %" del criterio de la Tarea 1d.

### Solape de conos

En todas las variantes hay solape sustancial. La Opción E lo reduce con fuerza (21.78 % → 4.59 % con
dirección ápice→objetivo) porque desacopla el radio de la longitud e impide que los conos largos se
inflen. **Ése es el efecto que E promete y lo cumple.** Lo que E no hace —ni pretende— es cubrir la
médula polar.

### Lectura de la Tarea 2

El paquete A+C+E **no resuelve la médula huérfana**, y en la parametrización literal solicitada la
empeora. La razón no está en ninguna de las tres opciones: está en que **cuatro de los diez ápices
están del lado cortical de la UCM**, lo que (i) impide aplicarles C y (ii) deja las cuatro pirámides
que deben cubrir los polos ancladas donde están. Mientras eso no se trate, A, C y E operan sólo sobre
las seis centrales.

---

## TAREA 3 — ¿Contradice algo del repo que los ápices sean centrales?

### Ni Capa 3c ni Capa 4 leen `pyramid_axis`

```bash
grep -c "pyramid_axis" 08_gemelo_digital/capa4_colector_alto.py 08_gemelo_digital/capa3c_colector.py
```
```
08_gemelo_digital/capa4_colector_alto.py:0
08_gemelo_digital/capa3c_colector.py:0
```

**Ambas consumen únicamente `pyramid_apex`.** La orientación del eje piramidal les es invisible. Ejes
muy oblicuos **no rompen nada** aguas arriba por esa vía. Esto **apoya la hipótesis de trabajo 2.**

### Capa 3c — crecimiento desde la papila

`08_gemelo_digital/capa3c_colector.py:105` — `apex = d0["pyramid_apex"].astype(np.float64)             # (10,3) papilas`

y `:143` — `nodos_k, parent_k, n_no_k, _ = cv.space_colonization(apex[k], atr_k, art_tree=None)`

El árbol colector crece por *space colonization* desde `apex[k]` hacia los extremos distales de las
nefronas asignadas a la pirámide k. **No asume ninguna orientación**: sólo una raíz y un conjunto de
atractores. Parámetros en `:48-49`: `DIST_INFLUENCIA = 8.0` mm, `DIST_MATAR = 0.6` mm.

El riesgo teórico —que atractores muy lejanos no se alcancen desde una raíz central— **no se
materializa hoy**:

```bash
.venv/bin/python -c "import numpy as np; d=np.load('capa3c_colector.npz',allow_pickle=True); print('n_no_alcanzados_por_piramide =',d['n_no_alcanzados_por_piramide']); print('dist_influencia',d['dist_influencia'],'dist_matar',d['dist_matar'])"
```
```
n_no_alcanzados_por_piramide = [0 0 0 0 0 0 0 0 0 0]
dist_influencia 8.0 dist_matar 0.6
```

**Cero atractores no alcanzados en las 10 pirámides**, incluida k=0 que hoy recibe 288 glomérulos
extendidos hasta X = −54.35 mm. El crecimiento tolera alcances largos desde una raíz central.

### Capa 4 — dos supuestos sobre los ÁPICES (no sobre los ejes)

**Supuesto 1 — el eje polar debe ser X, con fallo duro.** `08_gemelo_digital/capa4_colector_alto.py`:

```python
    rangos = [float(apex[:, i].max() - apex[:, i].min()) for i in range(3)]
    eje_max = int(np.argmax(rangos))
    ver_eje = (eje_max == EJE_POLAR)
```
```python
    if not ver_eje:
        print("  [FALLO DURO] el eje polar no es X; abortando.")
        sys.exit(1)
```

Rango actual de ápices: X = 23.76 mm, Y = 2.78 mm, Z = 8.80 mm. X domina con holgura. **Mover ápices
MÁS polares es seguro** para este test; **comprimirlos** por debajo de 8.80 mm en X lo rompería con
`sys.exit(1)`.

**Supuesto 2 — patrón bicalicial por signo de X.**

```python
    m_de_k = np.where(apex[:, EJE_POLAR] < 0.0, 0, 1).astype(np.int64)  # X<0 ->0, X>=0 ->1
```

Exige ápices repartidos a ambos lados de X = 0 (hoy 5/5). Cualquier reparametrización que los
concentre en un lado produce un cáliz mayor vacío.

**Supuesto 3 — la copa calicial se orienta hacia el hilio.**

```python
        pos = apex[k] + D_COPA * unit(HILIO - apex[k])
```

Geométricamente robusto para cualquier posición de ápice, pero desplaza las posiciones de los cálices
si los ápices se mueven.

### Respuesta

**NO ENCONTRADO nada en el repo que contradiga que los ápices permanezcan centrales.** Al contrario:
ambas capas consumen sólo `apex`, ignoran `pyramid_axis`, capa3c alcanza el 100 % de sus atractores
desde raíces centrales, y el único fallo duro de capa4 penaliza **comprimir** los ápices, no
mantenerlos. La hipótesis de trabajo 2 queda **sin contradicción localizable**.

**Con una salvedad que la matiza:** el hallazgo (b-bis) muestra que cuatro ápices están a 5.295 y
6.428 mm de la cápsula, del lado cortical de la UCM. Eso no es "central": es **superficial**. La
hipótesis de que los ápices son centrales por anatomía se sostiene; la posición concreta de esos
cuatro no se sigue de ella.

---

## TAREA 4 — `CONE_HALF_ANGLE_DEG` bajo criterio de teselación

Manteniendo ápices, ejes y longitudes actuales y variando sólo el semiángulo
(`r_base = length·tan(semiángulo)`):

```bash
.venv/bin/python -c "import numpy as np, capa0_dominio as C; d=np.load('capa0_dominio.npz',allow_pickle=True); co=d['coords'].astype(np.float64); dm=d['depth_cortical_mm'].astype(np.float64); med=dm>=float(C.GROSOR_CORTICAL_MM); ap=d['pyramid_apex']; ax=d['pyramid_axis']; L=d['pyramid_length']; T=[];R=[]; [ (T.append(v@ax[i]),R.append(np.linalg.norm((v:=co-ap[i])-np.outer(v@ax[i],ax[i]),axis=1))) for i in range(10)]; T=np.array([ (co-ap[i])@ax[i] for i in range(10)]); R=np.array([np.linalg.norm((co-ap[i])-np.outer((co-ap[i])@ax[i],ax[i]),axis=1) for i in range(10)]); print('semiang r_base_pol huerfana% solape%'); [print('%5d° %10.2f %9.2f %8.2f'%(h,L[0]*np.tan(np.radians(h)),100*((nc:=sum((med&(T[i]>=0)&(T[i]<=L[i])&(R[i]<=(T[i]/L[i])*(L[i]*np.tan(np.radians(h))))).astype(np.int32) for i in range(10)))==0)[med].sum()/med.sum(),100*(nc>=2)[med].sum()/med.sum())) for h in [22,25,30,35,40,45,50,55,60,70,80]]"
```

| semiángulo | `r_base` polar | **huérfana %** | solape (≥2 conos) % | cubre médula % |
|---|---|---|---|---|
| **22° (actual)** | 15.00 | **47.05** | 21.78 | 52.95 |
| 25° | 17.31 | 40.05 | 32.09 | 59.95 |
| 30° | 21.43 | 30.21 | 46.04 | 69.79 |
| 35° | 25.99 | 22.48 | 57.39 | 77.52 |
| 40° | 31.14 | 16.62 | 67.22 | 83.38 |
| 45° | 37.12 | 11.19 | 76.04 | 88.81 |
| 50° | 44.23 | 5.83 | 83.34 | 94.17 |
| 55° | 53.01 | 2.68 | 89.91 | 97.32 |
| 60° | 64.29 | 1.49 | 94.55 | 98.51 |
| 70° | 101.98 | 0.55 | 97.95 | 99.45 |
| 80° | 210.50 | 0.13 | 99.22 | 99.87 |

### ¿22° es plausible o está lejos?

**Está lejos, y la curva muestra que el problema no es el semiángulo.**

- Para bajar la huérfana al tramo "< 10 %" del criterio de la Tarea 1d haría falta **≈ 46–47°**, más
  del doble del valor actual.
- A 45°, `r_base` polar = **37.12 mm**, mayor que el semieje medial-lateral del riñón entero
  (`B_SEMI = 30.0`, `capa0_dominio.py:37`). El cono deja de ser una pirámide y pasa a ser medio
  órgano.
- Y el precio es el solape: a 45°, **el 76.04 % de la médula está dentro de dos o más conos**. La
  etiqueta `piramide_XX` deja de particionar y pasa a depender enteramente del desempate por `score`
  (`capa0_dominio.py:377`).

**No existe ningún semiángulo que cubra la médula sin destruir la partición.** La curva es monótona en
ambos sentidos: cobertura y solape suben juntos. Eso identifica el cuello de botella con precisión —
**no es `CONE_HALF_ANGLE_DEG`, es que 10 conos con ápices confinados a ±11.88 mm no pueden teselar una
médula que se extiende hasta ±48 mm.** El déficit es de *cobertura del dominio*, no de apertura.

Esto **no es una recomendación de valor**. Es la constatación de que 22° no está "mal calibrado": el
modelo cónico con ápices centrales no admite calibración que satisfaga el criterio de teselación.

---

## Resumen

| tarea | resultado |
|---|---|
| **1a** | Médula huérfana **38 721 / 82 297 = 47.050 %** |
| **1b** | 48.8 % polar (\|X\|>24.75), 24.4 % inter-hilera (\|Z\|<2.2). Más allá de \|X\|=40 mm, **100 %** huérfana. No es peri-sinusal |
| **1b-bis** | **Hallazgo nuevo:** 4 ápices (k=0,1,8,9) a **5.295 / 6.428 mm** de la cápsula → **en territorio córtex** |
| **1c** | **El repo NO sostiene la médula genérica.** Es residuo aritmético (`n_medulla − n_in_pyr`), sin justificación en código ni bitácora. La definición de CW de Glodny que el gemelo cita trata Bertin como córtex |
| **1d** | Criterio propuesto: `fraccion_huerfana`. Objetivo teórico 0 %; tramos <10 % / 10–25 % / >25 %. **Actual 47.05 % → tercer tramo.** Umbral [SUPUESTO DECLARADO], sin ancla bibliográfica |
| **2** | A+C+E **no resuelve**: 65–69 % con dirección ápice→objetivo (empeora), **44.17 %** en el mejor caso. **4/10 ápices sin solución para C.** E sí cumple lo suyo: solape 21.78 % → 4.59 % |
| **3** | **Nada contradice ápices centrales.** capa3c y capa4 leen sólo `apex`, ignoran `pyramid_axis` (grep = 0 en ambos); capa3c alcanza 100 % de atractores (`n_no_alcanzados = [0]×10`). Restricción real: capa4 aborta con `sys.exit(1)` si el rango X de ápices deja de dominar |
| **4** | 22° → 47.05 % huérfana. Haría falta **≈46–47°** para bajar de 10 %, con `r_base` = 37 mm (> `B_SEMI`=30) y **76 % de solape**. **No hay semiángulo que funcione** |

### Lo que este informe cambia respecto del anterior

1. **La Opción C queda bloqueada para las 4 pirámides polares** por el hallazgo (b-bis), no por la
   geometría del seno. Es un defecto nuevo, anterior a cualquier reparametrización.
2. **Corrijo mi lectura del tope de bases:** lo que las limitaba a ~19.9 mm era la dirección
   ápice→objetivo, no la UCM. Con dirección optimizada se alcanzan **42.14 mm**.
3. **La hipótesis de trabajo 2 sale reforzada** — nada aguas arriba depende de la orientación del eje.
4. **La Tarea 4 reencuadra el problema:** el déficit es de cobertura del dominio por 10 conos de
   ápice central, no de apertura angular. Ninguna de las opciones A–F, sola o combinada, ataca eso.

**Sin cambios en disco fuera de este archivo.** No se ejecutó `capa0_dominio.py`. No se hizo
`git add`, `git commit`, `git stash` ni `git restore`. Ningún `.py`, ningún `.npz` y ningún documento
existente fue modificado.
