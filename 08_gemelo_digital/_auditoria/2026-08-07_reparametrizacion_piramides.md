# Auditoría — Reparametrización de `build_pyramids` frente a la literatura lobar

**Programa:** Bio-Kidney AI 2026 · **Capa:** 0 · **Fecha:** 2026-08-07

> **NATURALEZA: SOLO LECTURA / DIAGNÓSTICO.**
> No se modificó ningún `.py`, ningún `.npz`, ningún parámetro ni ninguna lógica. No se ejecutó
> `capa0_dominio.py` (habría sobrescrito `capa0_dominio.npz`). Todas las cifras se recomputan **en
> memoria** desde `capa0_dominio.npz` y `capa1_nefronas.npz` mediante `np.load` y funciones puras
> importadas de `capa0_dominio` (el guard `if __name__ == "__main__":` de `capa0_dominio.py:610`
> impide que el import dispare `main()`).
>
> **NINGUNA DE LAS OPCIONES DE LA TAREA 2 ESTÁ APLICADA.** Las opciones A–F son un inventario de
> alternativas con sus consecuencias, no una recomendación ni una implementación. Los textos de la
> Tarea 5 son propuestas de redacción: **no se han escrito en `BITACORA.md` ni en
> `auditoria_correspondencia_anatomica.md`.**
>
> Todos los comandos de este documento son de lectura y **no escriben en disco**. Se ejecutan desde
> la raíz del repo (`~/Escritorio/BioKidney-AI`).

---

## Contexto: los anclajes de literatura que motivan esta auditoría

- **Bonsib SM**, *Heptinstall's Pathology of the Kidney*: 11–14 lóbulos; **~6 en polo superior, 4 en
  zona media, 4 en polo inferior**; **9–11 pirámides** en el adulto tras la fusión lobar; pirámides
  **compuestas** principalmente polares, **simples** en las mediopolares.
- **Nesterenko et al., Морфологія 2016** (150 riñones, 634 pirámides): **3–8 pirámides sólo en el
  extremo superior**, media **4.22 ± 0.15**.
- Múltiples fuentes: las **arterias arcuatas** siguen un curso en arco **paralelo a la cápsula** a lo
  largo de la unión córtico-medular. **La UCM no es un plano.**

### Consecuencia sobre `build_pyramids`

Dos literales quedan **contradichos**:

- `capa0_dominio.py:345` — `xf * 0.50 * a` confina las bases a X = ±24.75 mm (**45.0 %** del eje
  largo); la literatura sitúa pirámides en **ambos polos**.
- `capa0_dominio.py:346` — `+0.40 * b` pone las 10 bases en el plano **Y = 12.00**; medidas, están a
  **3.014–5.113 mm** de la cápsula, **todas dentro del córtex** bajo el propio umbral 6.6 mm de
  Capa 0.

Dos elementos quedan **validados**:

- `N_PIRAMIDES = 10` (`capa0_dominio.py:68`) cae en el rango 9–11 de Bonsib.
- Que las polares sean mayores: `r_base` 15.00 vs 12.98 mm, `length` 37.12 vs 32.12 mm — coherente
  con pirámides compuestas vs simples.

---

## TAREA 1 — La superficie UCM real de Capa 0

### (a) Puntos en la banda 6.6 ± 0.2 mm

```bash
.venv/bin/python -c "import numpy as np; d=np.load('capa0_dominio.npz',allow_pickle=True); dm=d['depth_cortical_mm'].astype(np.float64); m=np.abs(dm-6.6)<=0.2; print('n =',int(m.sum()),'=',round(100*m.mean(),3),'% del dominio')"
```

**5 412 puntos = 2.706 % del dominio.**

### (b) Caracterización de la superficie

```bash
.venv/bin/python -c "import numpy as np; d=np.load('capa0_dominio.npz',allow_pickle=True); co=d['coords'].astype(np.float64); dm=d['depth_cortical_mm'].astype(np.float64); P=co[np.abs(dm-6.6)<=0.2]; n=len(P); [print('%s [%.3f, %.3f] span %.3f'%(ax,P[:,i].min(),P[:,i].max(),P[:,i].ptp())) for i,ax in enumerate('XYZ')]; print('centroide',np.round(P.mean(axis=0),4)); w=np.linalg.lstsq(P**2,np.ones(n),rcond=None)[0]; s=1/np.sqrt(w); print('semiejes ajustados a=%.4f b=%.4f c=%.4f'%tuple(s)); lv=(P**2/s**2).sum(axis=1); print('|level-1| mediana %.5f p95 %.5f max %.5f'%(np.median(abs(lv-1)),np.percentile(abs(lv-1),95),abs(lv-1).max()))"
```

| eje | extensión | span |
|---|---|---|
| X | [−48.493, 48.488] | 96.982 |
| Y | [−22.286, 23.589] | 45.875 |
| Z | [−11.540, 11.563] | 23.104 |

Centroide `[−0.0008, 0.8564, 0.1220]`. El desplazamiento en +Y (0.856 mm) es la firma del seno: la
banda no existe en la cara medial excavada.

**Ajuste elipsoidal alineado a ejes: a = 46.2027, b = 23.1695, c = 11.2021 mm.**

**El ajuste es mediocre y la superficie NO es un elipsoide.** Residuo `|level−1|`: mediana 0.01744,
p95 0.04862, **max 0.10369**. En mm, distancia perpendicular al elipsoide ajustado: mediana 0.148,
p95 0.517, **max 2.331 mm**.

La razón es estructural: la UCM es la **superficie desplazada (offset) 6.6 mm hacia dentro** de un
elipsoide, y el offset de un elipsoide no es un elipsoide.

```bash
.venv/bin/python -c "import numpy as np, capa0_dominio as C; d=np.load('capa0_dominio.npz',allow_pickle=True); co=d['coords'].astype(np.float64); dm=d['depth_cortical_mm'].astype(np.float64); P=co[np.abs(dm-6.6)<=0.2]; [print('%-9s'%nm,'desv mm: mediana %.3f p95 %.3f max %.3f'%tuple(np.percentile(np.linalg.norm(P-C._nearest_point_ellipsoid(P,s),axis=1),[50,95,100]))) for nm,s in [('ajustado',np.array([46.2027,23.1695,11.2021])),('ingenuo',np.array([48.4,23.4,11.4]))]]"
```

| elipsoide candidato | desviación mediana | p95 | max |
|---|---|---|---|
| ajustado 46.20/23.17/11.20 | 0.148 mm | 0.517 mm | **2.331 mm** |
| ingenuo 48.4/23.4/11.4 (= cápsula − 6.6) | 0.388 mm | 1.013 mm | **1.312 mm** |

El elipsoide **ingenuo** tiene peor mediana pero **mejor máximo**. El ajustado gana en el grueso de
los puntos y pierde en los polos: `55 − 46.20 = 8.80` mm de retranqueo en X frente a 6.83 y 6.80 en
Y/Z. En el polo real, `X_max = 48.488 ≈ 55 − 6.6` — el offset ahí es 6.6, no 8.8. **El ajuste
elipsoidal sacrifica los polos**, que es justamente la región que la literatura obliga a poblar.

**Conclusión de (b): la UCM es aproximable por un elipsoide sólo con error de ~2.3 mm en los polos.
Cualquier reparametrización que coloque bases sobre "un elipsoide UCM" hereda ese error donde más
importa. El campo `depth_cortical_mm` sí da la UCM exacta, punto a punto.**

### (c) El plano Y = 12.00 frente a la UCM

Las 10 bases comparten un único valor de Y:

```bash
.venv/bin/python -c "import numpy as np, capa0_dominio as C; d=np.load('capa0_dominio.npz',allow_pickle=True); b=d['pyramid_apex']+d['pyramid_length'][:,None]*d['pyramid_axis']; print('Y unico:',np.unique(np.round(b[:,1],6))); db=C.capsule_distance(b); print('depth_mm bases:',np.round(db,3)); print('min %.3f max %.3f | deficit hasta 6.6: min %.3f max %.3f mediana %.3f'%(db.min(),db.max(),(6.6-db).min(),(6.6-db).max(),np.median(6.6-db)))"
```

```
Y unico: [12.]
depth_mm bases: [3.014 3.879 4.506 4.913 5.113 5.113 4.913 4.506 3.879 3.014]
min 3.014 max 5.113 | deficit hasta 6.6: min 1.487 max 3.586 mediana 2.094
```

Las 10 bases están dentro del córtex bajo el propio umbral de Capa 0. **Déficit hasta la UCM:
1.487 mm (polares) a 3.586 mm (centrales), mediana 2.094 mm.**

Y el plano completo, no sólo los 10 puntos:

```bash
.venv/bin/python -c "import numpy as np; d=np.load('capa0_dominio.npz',allow_pickle=True); co=d['coords'].astype(np.float64); dm=d['depth_cortical_mm'].astype(np.float64); pl=np.abs(co[:,1]-12.0)<=0.3; dp=dm[pl]; print('n en plano Y=12+-0.3:',int(pl.sum())); print('depth min %.3f mediana %.3f max %.3f'%(dp.min(),np.median(dp),dp.max())); print('|depth-6.6| mediana %.3f max %.3f'%(np.median(abs(dp-6.6)),abs(dp-6.6).max())); print('fraccion dentro del cortex: %.1f%%'%(100*(dp<6.6).mean()))"
```

```
n en plano Y=12+-0.3: 2675
depth min 0.005 mediana 6.047 max 15.526
|depth-6.6| mediana 3.429 max 8.926
fraccion dentro del cortex: 53.8%
```

**El plano Y = 12 barre profundidades de 0.005 a 15.526 mm.** Desviación respecto de la UCM:
**mediana 3.429 mm, máxima 8.926 mm**. Sólo el 53.8 % del plano cae en córtex. **Ningún plano puede
ser la UCM** — coincide con la fuente citada: la UCM no es un plano y las arcuatas la recorren en
arco paralelo a la cápsula.

---

## TAREA 2 — Opciones de reparametrización (NO IMPLEMENTADAS)

Los cuatro literales en juego, `capa0_dominio.py:333-334` y `:345-347`:

```python
        ux = xf * 0.60
        uz = zrow * 0.40
```

```python
        base[i] = np.array([xf * 0.50 * a,
                            +0.40 * b,
                            zrow * 0.62 * c])
```

### Restricción de acoplamiento que condiciona todas las opciones

`capa0_dominio.py:349-352`:

```python
    vec = base - apex                          # apice -> base (apice mira al seno)
    length = np.linalg.norm(vec, axis=1)
    axis = vec / length[:, None]
    r_base = length * np.tan(half)
```

**`r_base` no es libre: es `length · tan(22°)`.** Alejar una base aumenta la longitud y, con ella, el
radio del cono. Toda opción que empuje bases hacia los polos infla los conos polares salvo que se
toque también `CONE_HALF_ANGLE_DEG`.

---

### Opción A — Ampliar sólo el factor X de la base

| aspecto | detalle |
|---|---|
| **Literal que toca** | `capa0_dominio.py:345`, `xf * 0.50 * a` → `xf * k * a` con k > 0.50 (p. ej. k = 0.85 → X = ±46.75 mm) |
| **Ancla** | Bonsib (lóbulos en ambos polos: ~6 superior / 4 medio / 4 inferior); Nesterenko 2016 (3–8 pirámides sólo en el extremo superior, media 4.22 ± 0.15) |
| **Qué NO resuelve** | Las bases siguen en el plano Y = 12 (déficit 1.487–3.586 mm). No toca la UCM en absoluto |
| **Qué rompe aguas arriba** | `pyramid_length` crece; con k = 0.85 y ápice fijo la longitud polar pasa de 37.12 a ~48–49 mm y `r_base` a ~19.7 mm. Los conos polares se solapan con los mediopolares y el desempate por `score` (`capa0_dominio.py:377`) pasa a gobernar el reparto en lugar de la contención. Cambia `region_label` → invalida Capas 1–4 completas |
| **Compatible con r_base polar > central** | **Sí, y lo exagera** — pasaría de 1.16× a ~1.5×+. Puede sobrepasar lo que la categoría "pirámide compuesta" justifica |

### Opción B — Sustituir el plano por un elipsoide UCM

| aspecto | detalle |
|---|---|
| **Literal que toca** | `capa0_dominio.py:346`, `+0.40 * b` → componente Y derivada de un elipsoide UCM; en la práctica se reescribe `base[i]` entero como punto sobre ese elipsoide en la dirección elegida |
| **Ancla** | "La UCM no es un plano; las arcuatas la siguen en arco paralelo a la cápsula" |
| **Qué NO resuelve** | La extensión X, si no se combina con A. Y arrastra el error de ajuste de la Tarea 1b: **hasta 2.331 mm en los polos** con el elipsoide ajustado, 1.312 mm con el ingenuo |
| **Qué rompe aguas arriba** | Ídem A: cambia `region_label` y `pyramid_*` → invalida Capas 1–4 |
| **Compatible con r_base polar > central** | **Sí** — la relación polar/central la sigue fijando la longitud |

### Opción C — Colocar bases por inversión del campo `depth_cortical_mm`

| aspecto | detalle |
|---|---|
| **Literal que toca** | `capa0_dominio.py:345-347` completos, sustituidos por: elegir la dirección ápice→base y avanzar hasta que `capsule_distance` = `GROSOR_CORTICAL_MM` |
| **Ancla** | La definición misma de *cortical width* en Glodny 2009 (perpendicular cápsula → base de la pirámide medular) — la base **es** el punto a 6.6 mm |
| **Qué NO resuelve** | La extensión X (depende de la dirección elegida, o sea de A). Introduce una raíz numérica por pirámide dentro de `build_pyramids`, que hoy es analítica y cerrada |
| **Qué rompe aguas arriba** | Ídem A/B, más un coste estructural: `build_pyramids` pasa a depender de `capsule_distance`. El orden actual de `main()` ya lo permite — campo de profundidad en `:407`, pirámides en `:414` |
| **Compatible con r_base polar > central** | **Sí.** Es la única opción que garantiza déficit 0 por construcción, sin error de ajuste elipsoidal |

### Opción D — Reparametrizar también los ápices

| aspecto | detalle |
|---|---|
| **Literal que toca** | `capa0_dominio.py:333`, `ux = xf * 0.60` → factor mayor; y/o `:334`, `uz = zrow * 0.40` |
| **Ancla** | **Ninguna directa.** Bonsib describe distribución **lobar**, no posición de las papilas sobre la pared sinusal |
| **Qué NO resuelve** | Nada por sí sola; y **choca con un techo duro** — ver Tarea 3 |
| **Qué rompe aguas arriba** | `pyramid_apex` es la clave más consumida de la geometría piramidal: `capa1_nefronas.py:152`, `capa3c_colector.py:105`, `capa4_colector_alto.py:123`, más `capa3b_shunts_cola_corregible.py:228` y `capa3b_shunts_densidad_vs_ruteo.py:107`. Mover ápices reubica las 10 papilas y con ellas todo el árbol calicial de Capa 4 |
| **Compatible con r_base polar > central** | **Sí**, es ortogonal a `r_base` |

### Opción E — Desacoplar `r_base` de `length`

| aspecto | detalle |
|---|---|
| **Literal que toca** | `capa0_dominio.py:352`, `r_base = length * np.tan(half)` → `r_base` por pirámide, o `CONE_HALF_ANGLE_DEG` vectorizado |
| **Ancla** | Bonsib — pirámides **compuestas** (polares) vs **simples** (mediopolares) son categorías discretas, no un gradiente proporcional a la longitud. Sostiene fijar dos radios, no derivar uno |
| **Qué NO resuelve** | Ni la extensión X ni la UCM |
| **Qué rompe aguas arriba** | Poco directo: `pyramid_r_base` **no tiene consumidor** (defecto D-2 de la auditoría de contrato). Pero cambia `region_label` → invalida Capas 1–4 igualmente |
| **Compatible con r_base polar > central** | **Es la única que da control explícito** sobre la relación polar/central en lugar de heredarla de la longitud |

### Opción F — Distribución lobar no uniforme

| aspecto | detalle |
|---|---|
| **Literal que toca** | `capa0_dominio.py:328-329`, `f = (i + 0.5) / n_piramides` y `xf = (f - 0.5) * 2.0` → reparto no equiespaciado (p. ej. densidad ~6/4/4 superior/medio/inferior) |
| **Ancla** | **La más directa de todas** — Bonsib (11–14 lóbulos: ~6 superior, 4 medio, 4 inferior; 9–11 pirámides tras fusión) y Nesterenko 2016 (media 4.22 ± 0.15 en el extremo superior) |
| **Qué NO resuelve** | Ni la UCM ni el techo de los ápices |
| **Qué rompe aguas arriba** | Ídem; además **rompe la simetría X actual**, que hoy es exacta (`pyramid_length` es un palíndromo: 37.12, 35.11, 33.62, 32.62, 32.12, 32.12, 32.62, 33.62, 35.11, 37.12). Toda verificación aguas arriba que asuma simetría superior/inferior dejaría de valer |
| **Compatible con r_base polar > central** | Sí, pero **reasigna cuáles son "las polares"** — con 6 arriba y 4 abajo, el par extremo deja de ser simétrico |

**Nota transversal:** A, B, C y F cambian `region_label` y por tanto invalidan lo mismo. La diferencia
de coste entre aplicar una o las cuatro es marginal; la diferencia de riesgo, no.

---

## TAREA 3 — Interacción con el seno

```bash
.venv/bin/python -c "import numpy as np, capa0_dominio as C; from scipy.optimize import brentq; sx,sy,sz=C.SEMIEJES_SENO; cx,cy,cz=C.CENTRO_SENO; f=lambda ux,uz:(float(C.ellipsoid_level(np.array([[cx+ux*sx,cy+np.sqrt(max(0,1-ux*ux-uz*uz))*sy,cz+uz*sz]]),np.zeros(3),C.SEMIEJES)[0])-1.0); [print('uz=%.2f -> ux_max=%.4f apexX_max=%.3f mm (%.1f%% del semieje)'%(uz,r:=brentq(lambda u:f(u,uz),0.1,np.sqrt(1-uz*uz)-1e-9),r*sx,100*r*sx/55)) for uz in [0.0,0.20,0.40]]"
```

| `uz` | `ux_max` | apex X máx | % del semieje 55 |
|---|---|---|---|
| 0.00 | 0.9239 | **20.325 mm** | 37.0 % |
| 0.20 | 0.8987 | 19.771 mm | 35.9 % |
| **0.40 (actual)** | 0.8179 | **17.995 mm** | 32.7 % |

Estado actual: `ux = 0.54` → apex X = **±11.880 mm = 21.6 %** del semieje. Los 10 ápices están dentro
del parénquima con holgura (`level` ∈ [0.6907, 0.8072]).

```bash
.venv/bin/python -c "import numpy as np, capa0_dominio as C; d=np.load('capa0_dominio.npz',allow_pickle=True); ap=d['pyramid_apex']; lv=C.ellipsoid_level(ap,np.zeros(3),C.SEMIEJES); print('apex X [%.3f, %.3f] | level principal min %.4f max %.4f (<1 = dentro)'%(ap[:,0].min(),ap[:,0].max(),lv.min(),lv.max())); print('apex Y:',np.round(ap[:,1],3)); print('apex Z:',np.round(np.unique(ap[:,2]),3))"
```

```
apex X [-11.880, 11.880] | level principal min 0.6907 max 0.8072 (<1 = dentro)
apex Y: [-22.151 -20.966 -20.144 -19.621 -19.367 -19.367 -19.621 -20.144 -20.966 -22.151]
apex Z: [-4.4  4.4]
```

### El techo es duro y su causa es doble

1. **`SEMIEJES_SENO[0] = 22` acota geométricamente:** ningún punto del elipsoide del seno tiene
   |X| > 22. **Cobertura máxima absoluta del eje largo: 40.0 %** (2·22/110).
2. **Antes de llegar a 22, el ápice sale del parénquima.** Al aumentar `ux`,
   `uy = sqrt(1 − ux² − uz²)` cae y el ápice se desliza desde la pared medial (+Y del seno) hacia el
   ecuador del seno, es decir hacia Y = −34. Con `uz = 0.40`, en `ux = 0.85` el ápice ya está fuera
   (`level` = 1.0387).

Barrido literal que lo muestra:

```
 ux    apexX    uy     apexY    level_principal   dentro?
0.540   11.880  0.741  -22.151         0.8072   SI
0.600   13.200  0.693  -22.915         0.8371   SI
0.650   14.300  0.646  -23.662         0.8657   SI
0.700   15.400  0.592  -24.534         0.8983   SI
0.750   16.500  0.527  -25.571         0.9361   SI
0.800   17.600  0.447  -26.845         0.9813   SI
0.850   18.700  0.343  -28.515         1.0387   NO
0.916   20.163  0.005  -33.916         1.2134   NO
```

### Respuesta directa

Con `SEMIEJES_SENO = [22,16,11]` los ápices **no** pueden acompañar a bases polares. Su techo real es
**±18.0 mm (32.7 %)** manteniendo `uz = 0.40`, o **±20.3 mm (37.0 %)** si se aplana la hilera
anterior/posterior a `uz = 0`. Frente a bases que la literatura pide llevar hacia ±46.75 mm (85 % del
semieje), el desfase ápice/base sería de ~29 mm.

**Consecuencias, todas cuantificables hoy:**

- Los ejes piramidales polares se vuelven **muy oblicuos** respecto de X. Anatómicamente eso no es
  falso —las pirámides polares sí apuntan oblicuamente al seno—, pero
- `length` polar crece de 37.12 a ~48–49 mm y, por `capa0_dominio.py:352`, `r_base` de 15.00 a
  ~19.7 mm. Los conos polares se solapan con los mediopolares y el desempate por `score`
  (`capa0_dominio.py:377`) pasa a gobernar el reparto en lugar de la contención.

**¿Obliga a reabrir `SEMIEJES_SENO`? Sí, si se quiere cobertura polar de ápices por encima del 37 %.**
No, si se acepta que los ápices queden confinados al tercio central y sean los ejes los que se
alarguen oblicuamente.

**Advertencia de alcance:** `SEMIEJES_SENO` está **[SUPUESTO DECLARADO], no anclado** — puesto 2 de
los 10 pendientes en `09_paper_vascular/auditoria_correspondencia_anatomica.md:129`. Ya hay una
propuesta previa de cambiarlo, en `00_bitacora/BITACORA.md:691`: *"considerar subir a
SEMIEJES_SENO=(24,18,12) y/o CENTRO_SENO=(0,-32,0)"*. Subir a 24 sólo mueve el techo absoluto del
40.0 % al 43.6 %. **No alcanza para cobertura polar; el seno tendría que crecer mucho más, y eso
reabre el volumen de exclusión (ENTRADA 033 §2) y la contención de Capa 5a.**

---

## TAREA 4 — Alcance de la regeneración

Grafo de dependencias (`IN_*` de cada script en `08_gemelo_digital/`):

```bash
cd 08_gemelo_digital && for f in capa1_nefronas capa2_demanda capa3a_arterial capa3ab_peritubular capa3b_venoso capa3c_colector capa4_colector_alto; do echo "--- $f"; grep -n 'IN_[A-Z]* *= *"' $f.py; done
```

| script | lee | produce |
|---|---|---|
| `capa1_nefronas.py` | capa0 | `capa1_nefronas.npz` |
| `capa2_demanda.py` | capa0, capa1 | `capa2_demanda.npz` |
| `capa3a_arterial.py` | capa0, capa1, capa2 | `capa3a_arterial.npz` |
| `capa3ab_peritubular.py` | capa1, capa3a, capa0 | `capa3ab_peritubular.npz` |
| `capa3b_venoso.py` | capa3ab, capa0, capa3a | `capa3b_venoso.npz` |
| `capa3c_colector.py` | capa0, capa1 | `capa3c_colector.npz` |
| `capa4_colector_alto.py` | capa0 | `capa4_colector_alto.npz` |

### Vía de invalidación 1 — `pyramid_apex` / `pyramid_axis` / `pyramid_length` (lectura directa)

- `08_gemelo_digital/capa1_nefronas.py:152-154` → `capa1_nefronas.npz`
- `08_gemelo_digital/capa3c_colector.py:105` → `capa3c_colector.npz`
- `08_gemelo_digital/capa4_colector_alto.py:123` → `capa4_colector_alto.npz`

### Vía de invalidación 2 — `region_label` (cambia porque cambia `pyr_id`)

- `08_gemelo_digital/capa2_demanda.py:191`
- `08_gemelo_digital/capa3ab_peritubular.py:170`
- `08_gemelo_digital/capa3b_anastomosis_vv.py:82`
- `08_gemelo_digital/capa3b_clasificar_vv.py:213`
- `08_gemelo_digital/capa3b_severidad_vv.py:132`
- `08_gemelo_digital/capa3a_arterial_v1_backup.py:248`

### Orden de regeneración (topológico)

```
0. capa0_dominio.py             -> capa0_dominio.npz          [raíz]
1. capa1_nefronas.py            -> capa1_nefronas.npz         (apex/axis/length + region_label + depth_mm)
2. capa2_demanda.py             -> capa2_demanda.npz          (region_label)
3. capa3a_arterial.py           -> capa3a_arterial.npz
4. capa3ab_peritubular.py       -> capa3ab_peritubular.npz
5. capa3b_venoso.py             -> capa3b_venoso.npz
6. capa3b_reparar_colisiones.py -> capa3b_venoso_reparado.npz
   capa3c_colector.py           -> capa3c_colector.npz        (sólo necesita capa0+capa1; puede ir tras el paso 1)
7. capa4_colector_alto.py       -> capa4_colector_alto.npz    (sólo necesita capa0; puede ir tras el paso 0)
8. capa5a_*                     -> capa5a_meta.npz + 3 .stl
```

`capa3c_colector` y `capa4_colector_alto` son **paralelizables** (dependen de capa0/capa1, no de la
rama vascular). El camino crítico es **0→1→2→3a→3ab→3b→reparar**.

### `.npz` invalidados (9 canónicos)

`capa1_nefronas.npz`, `capa2_demanda.npz`, `capa3a_arterial.npz`, `capa3ab_peritubular.npz`,
`capa3b_venoso.npz`, `capa3b_venoso_reparado.npz`, `capa3c_colector.npz`, `capa4_colector_alto.npz`,
`capa5a_meta.npz`.

### Además quedan invalidados

Los diagnósticos que leen `pyramid_*`: `capa3b_shunts_cola_corregible.py:228-231`,
`capa3b_shunts_densidad_vs_ruteo.py:107-109`, `capa3b_auditoria_colisiones.npz`, y el documento
`09_paper_vascular/diagnostico_holgura_pelvis.md` (ya señalado como pendiente de recomputar en
ENTRADA 033 §6).

### NO invalidado

`02_vascular_cco/renal_data_v1.json` y el árbol CCO v8 — rama independiente, no lee capa0.

---

## TAREA 5 — Propuestas documentales (texto, SIN APLICAR)

> Nada de esta sección se ha escrito en `BITACORA.md` ni en
> `auditoria_correspondencia_anatomica.md`. Es redacción propuesta.

### (a) ENTRADA 034 — propuesta

> ## ENTRADA 034 — [fecha] — Sesgo polar mal atribuido, y dos literales de `build_pyramids` contradichos por la literatura lobar
>
> **Estado:** **DIAGNOSTICADO, NO CORREGIDO.** El sesgo polar no es un defecto de Capa 0. Y, por
> separado, la búsqueda de literatura contradice dos literales de `build_pyramids` que hasta ahora no
> estaban ni anclados ni declarados.
>
> ### 1. Corrección de atribución del sesgo polar
> - `00_bitacora/BITACORA.md:721` afirma: `sesgo polar en drenaje por pirámide (piramide_00 y 09 captan ~21%/20%, centrales ~5%) — artefacto de geometría de pirámides de Capa 0`.
>   `08_gemelo_digital/PLAN_AUDITORIA_GEMELO.md:63` lo lista igual, como deuda de Capa 0.
> - **La atribución a Capa 0 es incorrecta.** Ambas líneas se conservan sin modificar como registro
>   histórico; esta entrada las corrige.
> - Reparto real en Capa 0: min 3543, max 5702, **ratio 1.609**, y sigue el volumen del cono
>   (correlación de Pearson **0.957**, ratio de volumen 1.543, desviación máxima 1.225 puntos
>   porcentuales). **Es geometría, no defecto.**
> - El 21 %/5 % vive en `capa1_nefronas.npz`: `[288,168,60,73,73,74,65,67,174,258]` → 22.15 % /
>   4.62 %, **ratio 4.8**.
> - **Mecanismo:** `08_gemelo_digital/capa1_nefronas.py:113` reimplementa `assign_pyramids` como
>   Voronoi de eje más cercano, **sin contención ni escape**, incompatible con la contención en cono
>   de `capa0_dominio.py:356`. Los ejes piramidales cubren X ∈ [−24.75, 24.75]; la corteza llega a
>   ±54.35. **El 44.1 % de los glomérulos vive más allá del último eje, y el 69.5 % de esos cae en
>   k=0/k=9.**
> - **Causa próxima:** `pyramid_r_base` y `cone_half_angle_deg` se persisten y **no tiene consumidor
>   ninguna de las dos**. Sin el radio, Capa 1 no puede hacer contención.
> - **Consecuencia operativa: corregir Capa 0 no resuelve el 21 %/5 %.** El parche va en
>   `capa1_nefronas.py:113`.
>
> ### 2. Hallazgo de literatura
> - **Bonsib** (*Heptinstall's Pathology of the Kidney*): 11–14 lóbulos, ~6 polo superior / 4 zona
>   media / 4 polo inferior; 9–11 pirámides en adulto tras fusión lobar; pirámides **compuestas**
>   principalmente polares, **simples** en mediopolares.
> - **Nesterenko et al., Морфологія 2016** (150 riñones, 634 pirámides): 3–8 pirámides sólo en el
>   extremo superior, media **4.22 ± 0.15**.
> - Múltiples fuentes: las **arterias arcuatas** siguen un curso en arco **paralelo a la cápsula** a
>   lo largo de la UCM. **La UCM no es un plano.**
>
> ### 3. Lo que la literatura VALIDA
> - `N_PIRAMIDES = 10` (`capa0_dominio.py:68`) cae en el rango 9–11 de Bonsib → pasa de
>   [SUPUESTO DECLARADO] a **[ANCLADA]**.
> - Que las polares sean mayores: `r_base` 15.00 vs 12.98, `length` 37.12 vs 32.12 → coherente con
>   compuestas vs simples.
>
> ### 4. Lo que la literatura CONTRADICE
> - **`capa0_dominio.py:345`, `xf * 0.50 * a`:** confina las bases a X = ±24.75 mm, **45.0 % del eje
>   largo**. Bonsib y Nesterenko sitúan pirámides en **ambos polos**.
> - **`capa0_dominio.py:346`, `+0.40 * b`:** pone las 10 bases en el plano **Y = 12.00**. Medido,
>   están a **3.014–5.113 mm** de la cápsula — **todas dentro del córtex** bajo el propio umbral
>   6.6 mm de Capa 0. Déficit hasta la UCM: 1.487–3.586 mm, mediana 2.094 mm. El plano completo
>   Y = 12 barre profundidades de 0.005 a 15.526 mm (desviación mediana a la UCM 3.429 mm, máxima
>   8.926 mm): **ningún plano puede ser la UCM.**
>
> ### 5. Techo del seno
> Los ápices no pueden acompañar bases polares: `SEMIEJES_SENO = [22,16,11]` acota la cobertura
> absoluta al **40.0 %** del eje largo, y el ápice sale del parénquima antes, en **±17.995 mm
> (32.7 %)** con `uz = 0.40`. Subir a `[24,18,12]` (propuesta abierta en `BITACORA.md:691`) sólo lleva
> el techo al 43.6 %. **Cobertura polar de ápices exigiría reabrir `SEMIEJES_SENO` en magnitud, y con
> ello el volumen de exclusión (ENTRADA 033 §2) y la contención de Capa 5a.**
>
> ### 6. Estado epistémico de los cuatro literales
> `0.60` (`:333`), `0.40` (`:334`), `0.50` (`:345`), `0.40` (`:346`), más `0.62` (`:347`):
> **[SIN DECLARAR]** — ni ancla ni declaración de supuesto. Dos de ellos ahora **contradichos**.
> Entran al inventario de `auditoria_correspondencia_anatomica.md`.
>
> ### Estado
> **Capa 0 NO cumple el criterio de cierre.** Falla (a) por `pyramid_r_base`/`cone_half_angle_deg`
> huérfanas y (c) por los literales sin declarar. La cascada **no debe lanzarse** hasta decidir la
> reparametrización: congelaría los dos literales contradichos en las nueve capas regeneradas.

### (b) Actualización de `auditoria_correspondencia_anatomica.md` — propuesta

**Fila a MODIFICAR** (la de `N_PIRAMIDES`, hoy en `:48`):

> | `N_PIRAMIDES` (pirámides medulares) | capa0_dominio.py:**68** | 10 | sí (conteo) | **A** | **[ANCLADA — 2026-08-07]** El comentario "rango humano fisiológico 8–18" queda **superado**: la fuente da 9–11. | **Bonsib SM, *Heptinstall's Pathology of the Kidney*** — 11–14 lóbulos (~6 polo superior, 4 zona media, 4 polo inferior); **9–11 pirámides** en el adulto tras fusión lobar. Pirámides **compuestas** principalmente polares, **simples** en mediopolares. Corroborado por **Nesterenko et al., Морфологія 2016** (150 riñones, 634 pirámides): 3–8 en el extremo superior, media **4.22 ± 0.15**. | **[OK]** 10 dentro de 9–11 |

**Fila a CONSERVAR sin cambio de clase:**

> | `CONE_HALF_ANGLE_DEG` | capa0_dominio.py:**69** | 22.0 ° | sí (grados) | **A** | **[SIN DECLARAR]** — el comentario del código describe qué es, no de dónde sale. Gobierna `r_base = length·tan(half)` (`:352`) y con ello el 21.79 % del dominio que cae en pirámides. | `[CITA PENDIENTE]` | — |

**Cinco filas NUEVAS** (los literales de `build_pyramids`, hoy fuera del inventario):

> | literal `0.60` (factor X del ápice) | capa0_dominio.py:333 | — | no (fracción) | **A** | **[SIN DECLARAR]**. Fija apex X = ±11.88 mm (21.6 % del semieje). Techo geométrico ±17.995 mm con `uz=0.40`. | `[CITA PENDIENTE]` | — |
> | literal `0.40` (factor Z del ápice) | capa0_dominio.py:334 | — | no | **A** | **[SIN DECLARAR]**. Dos hileras ±4.4 mm. Reducirlo a 0 sube el techo de apex X a ±20.325 mm. | `[CITA PENDIENTE]` | — |
> | literal `0.50` (factor X de la base) | capa0_dominio.py:345 | — | no | **A** | **⚠ [SIN DECLARAR — CONTRADICHO 2026-08-07]** Confina bases a X = ±24.75 mm (45.0 % del eje largo). Bonsib/Nesterenko sitúan pirámides en ambos polos. | Bonsib; Nesterenko 2016 | **[CONTRADICHO]** |
> | literal `0.40` (factor Y de la base) | capa0_dominio.py:346 | — | no | **A** | **⚠ [SIN DECLARAR — CONTRADICHO 2026-08-07]** Coloca las 10 bases en el plano Y = 12.00; medidas a 3.014–5.113 mm de la cápsula, todas dentro del córtex (umbral 6.6 mm). La UCM no es un plano. | arcuatas en arco paralelo a la cápsula | **[CONTRADICHO]** |
> | literal `0.62` (factor Z de la base) | capa0_dominio.py:347 | — | no | **A** | **[SIN DECLARAR]** | `[CITA PENDIENTE]` | — |

### (c) Las 5 referencias `ruta:línea` caducas (defecto D-5) — corrección propuesta

Todas en `09_paper_vascular/auditoria_correspondencia_anatomica.md`:

| línea del doc | dice hoy | debe decir | constante |
|---|---|---|---|
| `:47` | `capa0_dominio.py:58` | `capa0_dominio.py:61` | `GROSOR_CORTICAL_MM` |
| `:47` | `(fracción, línea 63)` | `(fracción, línea 66)` | `UMBRAL_CM` |
| `:48` | `capa0_dominio.py:55` | `capa0_dominio.py:68` | `N_PIRAMIDES` |
| `:49` | `capa0_dominio.py:56` | `capa0_dominio.py:69` | `CONE_HALF_ANGLE_DEG` |
| `:50` | `capa0_dominio.py:58` | `capa0_dominio.py:71` | `N_POINTS` — **además resuelve la colisión**: `:58` estaba asignada a dos constantes distintas |
| `:113`, `:116` | `capa0_dominio.py:58` | `capa0_dominio.py:61` | `GROSOR_CORTICAL_MM` (prosa) |
| `:130` | `capa0_dominio.py:56` | `capa0_dominio.py:69` | `CONE_HALF_ANGLE_DEG` (lista de 10) |

Vigentes, sin tocar: `:41`→`:36`, `:42`→`:37`, `:43`→`:38`, `:44`→`:42`, `:45`/`:128`→`:48`,
`:46`/`:129`→`:49`.

**Sugerencia de método, no de contenido:** este es el segundo ciclo en que las referencias
`ruta:línea` de este documento caducan tras un commit que desplaza líneas. Un verificador de una
línea, ejecutable antes de cada cierre de capa, evitaría el tercero. No está escrito.

---

## Resumen

| tarea | resultado |
|---|---|
| **1a** | **5 412 pts** en 6.6 ± 0.2 mm (2.706 % del dominio) |
| **1b** | X ±48.49, Y [−22.29, 23.59], Z ±11.55. Ajuste 46.20/23.17/11.20 pero **no es un elipsoide**: hasta **2.331 mm** de error, concentrado en los polos |
| **1c** | Bases a **3.014–5.113 mm** (déficit 1.487–3.586). Plano Y=12 barre 0.005–15.526 mm; desviación a la UCM mediana **3.429**, máx **8.926 mm** |
| **2** | 6 opciones (A–F) con literal, ancla, límite, rotura aguas arriba y compatibilidad. **Sin recomendación** |
| **3** | Techo de ápices **±17.995 mm (32.7 %)**; máximo absoluto del seno **40.0 %**. **Sí obliga a reabrir `SEMIEJES_SENO`** si se quiere cobertura polar |
| **4** | **9 `.npz` canónicos** invalidados; camino crítico 0→1→2→3a→3ab→3b→reparar; 3c y 4 paralelizables |
| **5** | Textos propuestos para ENTRADA 034, 7 filas de `auditoria_correspondencia_anatomica.md` y las 7 correcciones de referencia |

**Veredicto heredado de la auditoría de contrato:** Capa 0 **NO** cumple el criterio de cierre.
Falla (a) —consumidor identificado— por `pyramid_r_base` y `cone_half_angle_deg` huérfanas, y falla
(c) —estado epistémico declarado— por `CONE_HALF_ANGLE_DEG` y los cinco literales de
`build_pyramids`. Esta auditoría añade que **dos de esos literales están además contradichos por la
literatura**.

**Sin cambios en disco.** No se ejecutó `capa0_dominio.py`. No se hizo `git add`, `git commit`,
`git stash` ni `git restore`. Ningún `.py`, ningún `.npz` y ningún documento existente fue
modificado para producir este informe.
