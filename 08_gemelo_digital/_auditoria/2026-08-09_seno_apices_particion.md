> # ⚠️ DOCUMENTO PARCIALMENTE SUPERADO — LEER CON EL FILTRO DE ABAJO
>
> **Revisado por `_auditoria/2026-08-12_reservas_pw_region2d.md` (2026-08-12).** Ese informe somete
> las conclusiones de éste a tres reservas metodológicas y **retira algunas**. Este documento **no es
> el estado actual**; se conserva íntegro y sin modificar por su valor forense.
>
> ## SUPERADO — no citar de este documento
>
> - **Cualquier intervalo o valor de `CENTRO_SENO[1]`**, en cualquier forma. La Tarea 2c de este
>   informe da tres vías a PW 15.4 (`cy = −30.6`, `sy = 19.4`, escala ×1.2125). **Están mal
>   planteadas:** se cumple la identidad `PW_max = |−B_SEMI − (cy + sy)|`, luego PW restringe sólo la
>   **suma** `cy + sy`, y el espacio real es **4D en `(cy, sy, sx, sz)`** — `SEMIEJES_SENO` entero
>   está sin anclar, no sólo `sy`. Lo que aquí aparece como "tres vías" son tres puntos de una banda
>   diagonal continua.
> - **La comparación «PW del gemelo vs Glodny» tal como está formulada en la Tarea 2b.** Se compara un
>   **máximo espacial sobre una geometría** (11.999 mm) contra una **media ± s.d. sobre n = 2068
>   riñones** (15.4 ± 2.8). **No son el mismo estadístico.** Es un error de magnitud de la misma clase
>   que el corregido por el método B (ENTRADA 032 §2): las unidades coinciden (mm) y la comparación
>   parece legítima, y no lo es. La justificación que se da aquí —"el análogo más directo a una medida
>   de corte axial único"— **no basta**.
> - **La premisa sobre el protocolo de medida de Glodny.** Se asumió una localización estandarizada.
>   Consultada la fuente primaria (BMC Urology 2009;9:19, PMC2813848), el paper especifica fase
>   arterial y plano axial, pero **no** el nivel anatómico, **ni** el número de medidas por riñón,
>   **ni** si es localización estandarizada o promedio. Los tres quedan **[PENDIENTE DE ANCLA]**.
> - **La afirmación de que las condiciones geométricas son anclas independientes.** Ápices, PW y
>   volumen excavado son **estrictamente monótonas crecientes** en `CENTRO_SENO[1]` y miden lo mismo
>   por métricas distintas: coinciden en la **dirección**, no convergen independientemente en un
>   **valor**. Y la condición de ápices resulta **redundante** — no fija ningún borde.
>
> ## VIGENTE — sigue siendo la referencia
>
> - **`CENTRO_SENO` es la variable dominante frente a `SEMIEJES_SENO`** (Tarea 1c). El barrido que lo
>   sostiene —`SEMIEJES_SENO` no baja de 2 ápices corticales ni escalando ×1.3; `CENTRO_SENO[1] = −32`
>   los lleva a 0— es geométrico y no depende de ninguna de las tres reservas.
> - **Los 4 ápices corticales (k = 0, 1, 8, 9) son un defecto real** (Tarea 1a-1b): `depth` 5.295 y
>   6.428 mm frente al umbral `GROSOR_CORTICAL_MM = 6.6`, con la contradicción documental de
>   `capa4_colector_alto.py:351,395,453`. **Matiz añadido el 08-12:** son defecto, **pero no criterio
>   de calibración** — resultan redundantes bajo la cota de Caglar.
> - **La especificación de la partición por vecindad a la papila** (Tarea 3) y su resultado numérico:
>   **0 % de médula huérfana** frente al 47.050 % de los conos, con la geometría actual y sin tocar el
>   seno. Independiente de todas las reservas.
> - **La dirección del déficit de parénquima** (Tarea 2b): el parénquima del gemelo es más delgado que
>   la referencia de Glodny. **Todos** los estadísticos candidatos caen bajo 12.6 mm (máximo 12.000 ·
>   p95 11.528 · mediana 6.855 · media 6.590). **La dirección es robusta; la magnitud del déficit
>   no** — va de 3.4 a 8.5 mm según el estadístico.
> - **La descomposición de la médula huérfana** y el techo de los ápices en el seno (Tarea 3, Tarea 1).
>
> **Referencia vigente:** `08_gemelo_digital/_auditoria/2026-08-12_reservas_pw_region2d.md`;
> fuentes en `04_literatura/anclas_seno_renal.md`; estado en `00_bitacora/BITACORA.md`, ENTRADA 033.
>
> **Aviso de trazabilidad:** «PW» **no existe como código en el repo** — es una construcción de
> auditoría introducida en este mismo documento (Tarea 2a). **No confundir con `depth_cortical_mm`**,
> que sí es un campo del `.npz` y mide otra cosa. Definición operativa en `_auditoria/2026-08-12`.

---

# Auditoría — Ápices corticales, control PW contra Glodny, y partición frente a conos

**Programa:** Bio-Kidney AI 2026 · **Capa:** 0 · **Fecha:** 2026-08-09

> **NATURALEZA: SOLO LECTURA / DIAGNÓSTICO.**
> No se modificó ningún `.py`, ningún `.npz`, ningún parámetro ni ninguna lógica. No se ejecutó
> `capa0_dominio.py`. Todas las cifras se recomputan **en memoria** desde `capa0_dominio.npz`,
> `capa3c_colector.npz` y funciones puras importadas de `capa0_dominio` (el guard
> `if __name__ == "__main__":` de `capa0_dominio.py:610` impide que el import dispare `main()`).
>
> **Ninguna configuración de seno, ningún esquema de partición y ninguna opción A–F está aplicada.**
> Los valores de las Tareas 2c y 3 son simulación; no existen en disco fuera de este documento.
>
> Continúa `_auditoria/2026-08-07_reparametrizacion_piramides.md` y
> `_auditoria/2026-08-08_paquete_ACE_simulacion.md`, de los que hereda la numeración de opciones A–F.

---

## Estado de los tres datos nuevos frente al repo

| dato aportado | estado en el repo |
|---|---|
| Emamian 1993 midió longitud, anchura, grosor **y volumen del área ecogénica central** (= seno renal) | **[PENDIENTE DE ANCLA]** — el repo cita Emamian 1993 **sólo** para los tres semiejes del órgano (`09_paper_vascular/auditoria_correspondencia_anatomica.md:41-43`, `:140-141`). **NO ENCONTRADO** ningún número del área ecogénica central en el repo. La afirmación de que la fuente los contiene queda registrada, pero **sus valores no están disponibles aquí** |
| Glodny 2009 da **PW = 15.4 ± 2.8 mm (dcho) / 15.9 ± 2.7 (izdo)**; espesor medular implícito **8.8 mm** | Parcialmente en el repo: `09_paper_vascular/auditoria_correspondencia_anatomica.md:47` dice `**NO CONFUNDIR con *parenchymal width* (PW) ≈ 15.5 mm**, que es cápsula → **seno renal** e incluye **corteza + médula**`. El valor exacto 15.4 ± 2.8 y la descomposición 6.6 + 8.8 son **nuevos**; los uso como ancla declarada |
| La huérfana se descompone 48.8 % polar / 51.2 % local | Recomputado y coincide (ver `_auditoria/2026-08-08`, Tarea 1b) |

**Consecuencia operativa inmediata:** la Tarea 2c se puede resolver **hoy**, porque sólo necesita PW
de Glodny. La calibración del seno contra Emamian **no**, porque sus números no están en el repo.

---

## TAREA 1 — El defecto de los cuatro ápices corticales

### (a) Cuáles, a qué profundidad, y por qué

```bash
.venv/bin/python -c "import numpy as np, capa0_dominio as C; d=np.load('capa0_dominio.npz',allow_pickle=True); ap=d['pyramid_apex']; da=C.capsule_distance(ap); [print('%2d apexX %8.3f apexY %8.3f depth %7.3f ux %6.3f uy %6.3f uz %6.3f %s'%(i,ap[i,0],ap[i,1],da[i],(u:=((i+0.5)/10-0.5)*2*0.60),np.sqrt(max(0,1-u*u-0.16)),(1.0 if i%2==0 else -1.0)*0.40,'<-- CORTEX' if da[i]<6.6 else '')) for i in range(10)]"
```

| i | apex X | apex Y | `depth` | `ux` | `uy` | |
|---|---|---|---|---|---|---|
| 0 | −11.880 | −22.151 | **5.295** | −0.540 | 0.741 | **CÓRTEX** |
| 1 | −9.240 | −20.966 | **6.428** | −0.420 | 0.815 | **CÓRTEX** |
| 2 | −6.600 | −20.144 | 7.187 | −0.300 | 0.866 | |
| 3 | −3.960 | −19.621 | 7.655 | −0.180 | 0.899 | |
| 4 | −1.320 | −19.367 | 7.879 | −0.060 | 0.915 | |
| 5..7 | simétricos | | | | | |
| 8 | 9.240 | −20.966 | **6.428** | 0.420 | 0.815 | **CÓRTEX** |
| 9 | 11.880 | −22.151 | **5.295** | 0.540 | 0.741 | **CÓRTEX** |

**Cuatro ápices (k = 0, 1, 8, 9) están por debajo del umbral `GROSOR_CORTICAL_MM = 6.6`
(`capa0_dominio.py:61`).**

**Causa, línea a línea.** Es un acoplamiento forzado por la esfera unidad. `capa0_dominio.py:333-339`:

```python
        ux = xf * 0.60
        uz = zrow * 0.40
        rad2 = ux * ux + uz * uz               # < 1 para estar en la pared +Y
        ...
        uy = np.sqrt(max(0.0, 1.0 - rad2))
```

y `capa0_dominio.py:340-342`:

```python
        apex[i] = np.array([cx + ux * sx,
                            cy + uy * sy,        # cara superior (+Y) del seno
                            cz + uz * sz])
```

El mecanismo tiene **dos brazos que se suman**, ambos actuando al alejarse del ecuador:

1. **El ápice baja hacia la cápsula medial.** `uy = sqrt(1 − ux² − uz²)` decrece monótonamente con
   |ux|. Como `apex_Y = cy + uy·sy = −34 + 16·uy`, al pasar de uy = 0.915 (k=4) a uy = 0.741 (k=0) el
   ápice cae de Y = −19.367 a Y = **−22.151**: se desplaza **2.78 mm hacia la cápsula medial**.
2. **La cápsula medial sube hacia el centro.** El elipsoide principal en (X, Z) tiene su superficie
   medial en `Y = −30·sqrt(1 − (X/55)² − (Z/18)²)`. En k=4 (X=−1.32) eso da Y = −29.08; en k=0
   (X=−11.88) da Y = −28.36: la cápsula se acerca **0.72 mm**.

Los dos brazos convergen y el hueco perpendicular pasa de 7.879 a 5.295 mm.

**El comentario del código describe la intención correcta pero no la garantiza.** `capa0_dominio.py:335`:

```python
        rad2 = ux * ux + uz * uz               # < 1 para estar en la pared +Y
```

`rad2 < 1` asegura que el punto esté **sobre el elipsoide del seno**; **no** asegura que esté a más
de 6.6 mm de la cápsula. No existe en `build_pyramids` ninguna comprobación de profundidad.
**NO ENCONTRADO** ningún test que lo verifique en `capa0_dominio.py` ni en el bloque `VERIFICACION`
(`:513-608`).

### (b) ¿Se propaga aguas arriba?

**Sí, como contradicción documental — no como fallo numérico.**

**Capa 3c** trata el ápice como conducto de Bellini, que es una estructura **medular** por
definición. `08_gemelo_digital/capa3c_colector.py:59-61`:

```python
#   radio_papila   = 200 um (r) -> 400 um (diam): conducto PAPILAR de Bellini,
#     ...
#     PMC). El conducto de Bellini es UNO por piramide, en el apice papilar.
```

y `:14` — `  raiz_k        = pyramid_apex[k]                         (papila, Capa 0)`

**Capa 4 lo declara explícitamente médula/seno, y lo escribe en su informe generado.**
`08_gemelo_digital/capa4_colector_alto.py:351`:

```python
    pap_val = val[idx_papila]  # papilas: interfaz medula/seno, exentas
```

`:395`:

```python
    print(f"        papila_junction (EXENTA, interfaz medula/seno): "
```

`:453-454`:

```python
    L.append(f"`papila_junction` (10) es nodo-interfaz medula/seno, **EXENTA** del test "
             f"(se sienta sobre la pared +Y del seno). Valor del elipsoide en las papilas: "
```

Ese texto ya viajó al repo: `09_paper_vascular/auditoria_capa4_calicial.md:28` reproduce
`` `papila_junction` (10) es nodo-interfaz medula/seno, **EXENTA** del test ``.

**Lectura:** Capa 4 afirma en un documento versionado que las 10 papilas son interfaz **médula**/seno.
Para cuatro de ellas, Capa 0 dice que están a 5.295 y 6.428 mm de la cápsula, es decir en territorio
**córtex** bajo su propio umbral. **Las dos capas se contradicen.**

Ninguna de las dos rompe numéricamente: Capa 4 **exime** a las papilas del test de contención y
ninguna lee `region_label` en el ápice. El defecto es de coherencia, no de ejecución — pero es
exactamente el tipo de afirmación que acaba en una figura o en el preprint.

### (c) ¿Qué variable domina? Barrido de una en una

```bash
.venv/bin/python -c "import numpy as np, capa0_dominio as C
def ad(uxf=0.60,uzf=0.40,cyv=-34.0,SS=None):
    SS=C.SEMIEJES_SENO if SS is None else np.array(SS,float); o=[]
    for i in range(10):
        xf=((i+0.5)/10-0.5)*2.0; z=1.0 if i%2==0 else -1.0
        ux=xf*uxf; uz=z*uzf; uy=np.sqrt(max(0,1-ux*ux-uz*uz))
        o.append(float(C.capsule_distance(np.array([[ux*SS[0], cyv+uy*SS[1], uz*SS[2]]]))[0]))
    return np.array(o)
print('ux ',[(v,round(ad(uxf=v).min(),3),int((ad(uxf=v)<6.6).sum())) for v in [0.30,0.40,0.50,0.60,0.70,0.80]])
print('uz ',[(v,round(ad(uzf=v).min(),3),int((ad(uzf=v)<6.6).sum())) for v in [0.0,0.20,0.40,0.60]])
print('cy ',[(v,round(ad(cyv=v).min(),3),int((ad(cyv=v)<6.6).sum())) for v in [-38,-36,-34,-32,-30,-28]])
print('esc',[(s,round(ad(SS=[22*s,16*s,11*s]).min(),3),int((ad(SS=[22*s,16*s,11*s])<6.6).sum())) for s in [0.8,0.9,1.0,1.1,1.2,1.3]])"
```

| variable | valores → (depth mín, nº ápices < 6.6) |
|---|---|
| factor `ux` (`:333`) | 0.30 → (7.329, **0**) · 0.40 → (6.848, **0**) · 0.50 → (6.184, 2) · **0.60 → (5.295, 4)** · 0.70 → (4.111, 4) · 0.80 → (2.515, 6) |
| factor `uz` (`:334`) | **0.00 → (8.689, 0)** · 0.20 → (7.574, **0**) · **0.40 → (5.295, 4)** · 0.60 → (2.058, **10**) |
| `CENTRO_SENO[1]` (`:48`) | −38 → (1.979, **10**) · −36 → (3.691, 10) · **−34 → (5.295, 4)** · **−32 → (6.757, 0)** · −30 → (8.058, 0) · −28 → (9.190, 0) |
| escala de `SEMIEJES_SENO` (`:49`) | ×0.8 → (4.065, **10**) · ×0.9 → (4.735, 6) · **×1.0 → (5.295, 4)** · ×1.1 → (5.745, 2) · ×1.2 → (6.091, **2**) · ×1.3 → (6.340, **2**) |

**Variable dominante: `CENTRO_SENO[1]`.** Es la que más rápido cruza el umbral y la única que lo hace
con un cambio pequeño y en la dirección que además corrige el PW (Tarea 2). Pasar de −34 a −32 lleva
los 10 ápices por encima de 6.6.

**Resultado que invierte la intuición previa: `SEMIEJES_SENO` NO puede resolver el defecto.** Ni
siquiera escalando ×1.3 baja de **2** ápices corticales, y la curva se aplana (5.295 → 5.745 → 6.091
→ 6.340: rendimientos decrecientes). Agrandar el seno mueve el ápice hacia +Y, pero también lo
empuja hacia fuera en X, donde la cápsula está más cerca. **Los dos efectos casi se cancelan.**

`uz = 0.40` es la segunda variable más sensible: llevarla a 0.60 pone **los 10** ápices en córtex, y
llevarla a 0 los saca todos. Es el factor de las dos hileras anterior/posterior.

**Conclusión de la Tarea 1:** el defecto **no** es que el seno sea pequeño. Es que su **centro está
demasiado afuera** (−34, es decir 4 mm fuera del elipsoide principal, cuyo `B_SEMI` es 30), lo que
obliga a la pared +Y a acercarse a la cápsula en cuanto el ápice se aparta del ecuador.

---

## TAREA 2 — Control PW contra Glodny

### (a) Ancho parenquimatoso del gemelo

PW = espesor del parénquima desde la cápsula hasta la pared del seno. Medido como en un corte axial:
para cada (X, Z) donde ambos existen, distancia entre la cápsula medial y la pared +Y del seno.

```bash
.venv/bin/python -c "import numpy as np
a,b,c=55.,30.,18.; sx,sy,sz=22.,16.,11.; cy=-34.
N=400; X=np.linspace(-a,a,N); Z=np.linspace(-c,c,N); XX,ZZ=np.meshgrid(X,Z,indexing='ij')
r=1-(XX/a)**2-(ZZ/c)**2; okc=r>0; ycap=np.where(okc,-b*np.sqrt(np.clip(r,0,None)),np.nan)
rs=1-(XX/sx)**2-(ZZ/sz)**2; oks=rs>0; ysen=np.where(oks,cy+sy*np.sqrt(np.clip(rs,0,None)),np.nan)
m=okc&oks&(ysen>ycap); w=(ysen-ycap)[m]
print('n %d mediana %.3f p5 %.3f p95 %.3f min %.3f MAX %.3f'%(len(w),np.median(w),np.percentile(w,5),np.percentile(w,95),w.min(),w.max()))
print('dentro de [12.6,18.2]: %.1f %% | por DEBAJO de 12.6: %.1f %%'%(100*((w>=12.6)&(w<=18.2)).mean(),100*(w<12.6).mean()))"
```

| estadístico | valor |
|---|---|
| mediana | **6.853 mm** |
| p5 | 0.841 mm |
| p95 | 11.531 mm |
| mínimo | 0.001 mm |
| **máximo** (eje hiliar, X=0, Z=0) | **11.999 mm** |

### (b) Comparación con Glodny — el parénquima del gemelo es MÁS DELGADO

**Glodny PW = 15.4 ± 2.8 mm → rango ±1 s.d. = [12.6, 18.2].**

**El 100 % del PW del gemelo cae por debajo de 12.6 mm. Fracción dentro del rango de Glodny: 0.0 %.**

**Responde a la pregunta en sentido contrario al planteado: el parénquima del gemelo no es más
grueso, es sistemáticamente más delgado.** Ni su valor máximo (11.999 mm, en el eje hiliar, que es el
análogo más directo a una medida de corte axial único) alcanza el límite inferior del rango.

La descomposición lo hace nítido. En el eje hiliar (X = 0, Z = 0): cápsula medial en Y = −30.000,
pared del seno en Y = −34 + 16 = −18.000.

| | PW | córtex (CW) | médula |
|---|---|---|---|
| **Glodny 2009** | 15.4 mm | 6.6 mm | **8.8 mm** |
| **gemelo (eje hiliar)** | **12.0 mm** | 6.6 mm | **5.4 mm** |

**El déficit está enteramente en la médula: 5.4 mm donde la fuente implica 8.8 — un 39 % menos.** El
córtex es correcto por construcción (el umbral es 6.6 mm), así que todo el error de PW se descarga
sobre la médula.

Esto **explica la médula huérfana desde otro ángulo**: no sólo los conos no la cubren (informe del
2026-08-08); es que además hay **menos médula de la que debería haber**. Son dos defectos
independientes que apuntan al mismo sitio.

### (c) ¿Qué `SEMIEJES_SENO` da PW = 15.4?

La relación es cerrada en el eje hiliar: `PW_max = CENTRO_SENO[1] + SEMIEJES_SENO[1] + B_SEMI`.

```bash
.venv/bin/python -c "import numpy as np
a,b,c=55.,30.,18.
def pw(sxv,syv,szv,cyv=-34.,N=400):
    X=np.linspace(-a,a,N); Z=np.linspace(-c,c,N); XX,ZZ=np.meshgrid(X,Z,indexing='ij')
    r=1-(XX/a)**2-(ZZ/c)**2; okc=r>0; ycap=np.where(okc,-b*np.sqrt(np.clip(r,0,None)),np.nan)
    rs=1-(XX/sxv)**2-(ZZ/szv)**2; oks=rs>0; ysen=np.where(oks,cyv+syv*np.sqrt(np.clip(rs,0,None)),np.nan)
    m=okc&oks&(ysen>ycap); return (ysen-ycap)[m]
print('sy  ',[(v,round(np.median(pw(22,v,11)),3),round(pw(22,v,11).max(),3)) for v in [16,18,19.4,20,22,24]])
print('cy  ',[(v,round(np.median(pw(22,16,11,v)),3),round(pw(22,16,11,v).max(),3)) for v in [-34,-32,-30.6,-30,-28]])
print('esc ',[(s,round(np.median(pw(22*s,16*s,11*s)),3),round(pw(22*s,16*s,11*s).max(),3)) for s in [1.0,1.1,1.2,1.3]])"
```

Tres formas de llegar a **PW_max = 15.4**, con `SEMIEJES = 55/30/18` intacto:

| vía | valor | PW_max |
|---|---|---|
| sólo `SEMIEJES_SENO[1]` | **sy = 19.4** (seno 22 / **19.4** / 11) | 15.399 |
| sólo `CENTRO_SENO[1]` | **cy = −30.6** | 15.400 |
| escala uniforme | **×1.2125** (seno 26.7 / 19.4 / 13.3) | 15.400 |

**Aviso sobre la mediana.** Llevar la **mediana** a 15.4 no es alcanzable con un seno razonable:
sy = 24 sólo da mediana 12.028. La mediana se toma sobre toda la huella del seno, incluidas sus
orillas donde PW → 0, y **no es el análogo de la medida de Glodny**, que es una medida en un corte
concreto a la altura del hilio. **El estadístico comparable es PW_max (eje hiliar) = 11.999 mm.**
Doy los valores para PW_max y lo declaro explícitamente.

### (d) Compatibilidad con el reparto 6.6 / 8.8 — y la convergencia con la Tarea 1

Por construcción, sí: 6.6 + 8.8 = 15.4. Como el umbral cortical es absoluto (`GROSOR_CORTICAL_MM`,
`capa0_dominio.py:61`), llevar PW_max a 15.4 deja automáticamente 8.8 mm de médula en el eje hiliar.

Lo interesante es que **las tres vías no son equivalentes**. Evaluadas contra los tres criterios a la
vez:

```bash
.venv/bin/python -c "import numpy as np, capa0_dominio as C
def ad(SS,cyv):
    o=[]
    for i in range(10):
        xf=((i+0.5)/10-0.5)*2.0; z=1.0 if i%2==0 else -1.0
        ux=xf*0.60; uz=z*0.40; uy=np.sqrt(max(0,1-ux*ux-uz*uz))
        o.append(float(C.capsule_distance(np.array([[ux*SS[0], cyv+uy*SS[1], uz*SS[2]]]))[0]))
    return np.array(o)
def ve(SS,cyv,n=400000):
    rng=np.random.default_rng(3); u=rng.normal(size=(n,3)); u/=np.linalg.norm(u,axis=1,keepdims=True)
    P=(u*rng.uniform(0,1,(n,1))**(1/3))*C.SEMIEJES
    return 4/3*np.pi*np.prod(C.SEMIEJES)*np.mean(C.ellipsoid_level(P,np.array([0.,cyv,0.]),np.array(SS,float))<1.0)/1000
for nm,SS,cy in [('ACTUAL',[22,16,11],-34.0),('cy=-30.6',[22,16,11],-30.6),('sy=19.4',[22,19.4,11],-34.0),('escala x1.2125',[22*1.2125,16*1.2125,11*1.2125],-34.0)]:
    a=ad(SS,cy); print('%-15s PWmax %6.2f  apex_min %7.3f  n<6.6 %2d  vol_excl %6.3f mL'%(nm,cy+SS[1]+30,a.min(),(a<6.6).sum(),ve(SS,cy)))"
```

| configuración | PW_max | `depth` mín de ápice | ápices < 6.6 | volumen excluido |
|---|---|---|---|---|
| **ACTUAL** (22/16/11, cy=−34) | 12.00 | 5.295 | **4** | 3.993 mL |
| **cy = −30.6** (22/16/11) | **15.40** | **7.685** | **0** | 6.214 mL |
| **sy = 19.4** (22/19.4/11) | **15.40** | 7.110 | **0** | 5.591 mL |
| escala ×1.2125 | 15.40 | 6.127 | **2** | 7.413 mL |

*(El 3.993 mL del caso ACTUAL es una estimación Monte Carlo con semilla y N distintos del `.npz`
canónico, que da 3.9599 mL según ENTRADA 033 §2. La diferencia es ruido de muestreo, no discrepancia.)*

**Convergencia:** las dos vías que llevan PW_max a 15.4 **también eliminan los cuatro ápices
corticales**, y la que menos ayuda a los ápices (escala uniforme, deja 2) es además la que más
volumen de exclusión añade. **Mover `CENTRO_SENO[1]` de −34 a −30.6 corrige simultáneamente el
defecto de la Tarea 1 y el déficit de PW de la Tarea 2**, y es la variable que la Tarea 1c ya había
señalado como dominante por una vía independiente.

**Coste declarado, no evaluado:** cualquiera de las tres sube el volumen excluido de ~4 mL a
5.6–7.4 mL. Eso **reabre ENTRADA 033** (la distinción 16.22 mL / 3.9599 mL) y toca la contención de
Capa 5a. **Y ninguna de las tres está anclada a Emamian:** todas se derivan de ajustar PW a Glodny,
que es un ancla sobre el *parénquima*, no sobre el *seno*. Marco los tres valores
**[PENDIENTE DE ANCLA]** en cuanto a las dimensiones del seno propiamente dichas.

---

## TAREA 3 — Partición por vecindad a la papila, frente a los conos

Esquema evaluado: cada punto de médula se asigna a la papila (`pyramid_apex`) **más cercana**.
Diagrama de Voronoi de 10 semillas, restringido a la médula.

```bash
.venv/bin/python -c "import numpy as np; d=np.load('capa0_dominio.npz',allow_pickle=True); co=d['coords'].astype(np.float64); dm=d['depth_cortical_mm'].astype(np.float64); med=dm>=6.6; ap=d['pyramid_apex']; M=co[med]; lab=np.linalg.norm(M[:,None,:]-ap[None,:,:],axis=2).argmin(axis=1); cnt=np.array([int((lab==i).sum()) for i in range(10)]); print('conteo',list(cnt)); print('%share',[round(100*x/cnt.sum(),2) for x in cnt]); print('min %d max %d ratio %.3f'%(cnt.min(),cnt.max(),cnt.max()/cnt.min())); print('polares(0,1,8,9) %.2f %% | centrales %.2f %% | ratio medias %.3f'%(100*cnt[[0,1,8,9]].sum()/cnt.sum(),100*cnt[2:8].sum()/cnt.sum(),cnt[[0,1,8,9]].mean()/cnt[2:8].mean()))"
```

```
conteo [13439, 9251, 6226, 6100, 6020, 6038, 6004, 6161, 9372, 13686]
%share [16.33, 11.24, 7.57, 7.41, 7.31, 7.34, 7.3, 7.49, 11.39, 16.63]
min 6004 max 13686 ratio 2.279
polares(0,1,8,9) 55.59 % | centrales 44.41 % | ratio medias 1.878
```

### Respuestas punto por punto

**¿0 % de médula huérfana por construcción?** **Sí.** Todo punto de médula tiene una papila más
cercana; no hay test de contención que pueda fallar. La huérfana pasa de **47.05 % → 0 %** sin tocar
ningún parámetro geométrico.

**¿Qué ancla lo sostiene o lo contradice?**

- **A favor — el lóbulo como unidad del desarrollo.** Bonsib describe el riñón como 11–14 lóbulos
  fusionados, cada uno una pirámide con su córtex suprayacente. Una partición del parénquima en
  territorios, uno por papila, **es** la descripción lobar. Un cono independiente no lo es: los conos
  pueden solaparse (hoy el 21.78 % de la médula está en ≥2 conos) y dejar huecos, cosa que los
  lóbulos no hacen.
- **A favor — las columnas de Bertin como frontera.** Anatómicamente, lo que separa dos pirámides
  adyacentes es una columna de Bertin. En un Voronoi sembrado por papilas, la frontera entre celdas
  cae exactamente donde equidistan dos papilas: es el lugar geométrico de las columnas. El repo ya
  trata Bertin como frontera y como córtex, en la definición de CW que cita
  (`09_paper_vascular/auditoria_correspondencia_anatomica.md:47`): `**evitando las columnas de Bertin y el seno renal**`.
- **En contra — la pirámide es un cono, no un poliedro.** La pirámide medular tiene forma cónica con
  base en la UCM y ápice en la papila. Las celdas de Voronoi son poliedros de caras planas, y su
  "base" no queda a profundidad constante. La partición **acierta el reparto y pierde la forma**.
- **[PENDIENTE DE ANCLA]:** no localizo ninguna fuente, en el repo ni entre las tres citadas, que
  cuantifique cuánto se aparta la pirámide real de un territorio de Voronoi. La objeción es
  cualitativa.

**¿Qué se pierde?** `pyramid_r_base` y `cone_half_angle_deg` dejan de tener significado.

**¿Qué consumidor los necesitaría? Ninguno.** Ambas claves tienen **cero consumidores** en todo el
repo (defecto D-2 de `_auditoria/2026-08-07`, Tarea 1b):

```bash
cd 08_gemelo_digital && grep -rn '\["pyramid_r_base"\]\|\["cone_half_angle_deg"\]' *.py | wc -l
```
→ `0`

Retirarlas **no rompe nada aguas arriba**. Y `CONE_HALF_ANGLE_DEG` es hoy la única constante
**[SIN DECLARAR]** de `capa0_dominio.py` (`:69`): la partición la elimina en lugar de anclarla.

**¿Reproduce que las polares sean mayores?** **Sí, y con más fuerza que los conos.**

| esquema | %share polar (k=0,9) | %share central | ratio max/min | ratio medias polar/central |
|---|---|---|---|---|
| conos (actual) | 12.87 / 13.09 | 8.13–8.63 | 1.609 | ~1.45 |
| **Voronoi** | **16.33 / 16.63** | 7.30–7.57 | **2.279** | **1.878** |

Las polares captan el **55.59 %** de la médula entre cuatro celdas. Eso es coherente con Bonsib
(compuestas polares, simples mediopolares) **sin ningún parámetro que lo imponga**: emerge de que las
papilas polares tienen más parénquima a su alrededor. **Es el hecho anclado reproducido por
construcción, no por calibración** — a diferencia de la Opción E, que lo impone fijando dos radios.

**Salvedad honesta:** ese 2.279 se obtiene con los ápices **actuales**, cuatro de los cuales están en
territorio córtex (Tarea 1). Si los ápices se corrigen, el ratio cambiará. El número de arriba es del
estado de hoy, no una predicción.

### Comparación contra A–F

| | A (X base) | C (base sobre UCM) | E (r_base categórico) | F (lobar no uniforme) | **Partición (nueva, "G")** |
|---|---|---|---|---|---|
| huérfana resultante | 65–69 % | — (no aplicable a 4/10) | — | — | **0 % por construcción** |
| resuelve cobertura polar | parcial | no | no | no | **sí** |
| resuelve solape | no | no | **sí** (21.8→4.6 %) | no | **sí** (no existe) |
| reproduce polares mayores | exagera (~1.5×+) | — | lo **impone** | reasigna cuáles | **lo produce** (1.878) |
| ancla | Bonsib, Nesterenko | Glodny (CW = cápsula→base) | Bonsib (compuestas/simples) | Bonsib, Nesterenko | Bonsib (lóbulo, Bertin) |
| toca | `:345` | `:345-347` | `:352` | `:328-329` | **`:356-382`** (`assign_pyramids`) |
| bloqueada por los 4 ápices corticales | no | **sí, 4/10** | no | no | **no** |
| `pyramid_r_base` / `cone_half_angle_deg` | siguen | siguen | E los reescribe | siguen | **quedan sin sentido (0 consumidores)** |
| cambia `region_label` | sí | sí | sí | sí | **sí** |
| invalida Capas 1–4 | sí | sí | sí | sí | **sí** |

**La partición no compite con A, C, E ni F: opera en otra función.** A/C/F cambian
`build_pyramids` (`:305`); E cambia `:352`; la partición cambia `assign_pyramids` (`:356`). Es
combinable con cualquiera de ellas — y también es la única que hace innecesarias las Opciones A y E,
porque cubre sin ampliar bases y reparte sin fijar radios.

**Lo que la partición NO resuelve:** el déficit de médula de la Tarea 2 (5.4 mm frente a 8.8). Repartir
mejor una médula demasiado delgada no la engrosa. **Ese defecto sólo lo toca el seno.**

*(No recomiendo ninguna opción. La tabla compara; la decisión no está en este informe.)*

---

## TAREA 4 — Orden de decisión

### Decisiones que pueden tomarse YA (ancla disponible en el repo o en Glodny)

| # | decisión | ancla | por qué no está bloqueada |
|---|---|---|---|
| **D1** | Reconocer los 4 ápices corticales como defecto y su variable dominante (`CENTRO_SENO[1]`) | `GROSOR_CORTICAL_MM = 6.6`, Glodny — ya en el repo (`capa0_dominio.py:61`) | El umbral que viola ya está anclado. No necesita Emamian |
| **D2** | Fijar el objetivo **PW_max = 15.4 mm** como criterio de aceptación del seno | Glodny 2009, PW 15.4 ± 2.8 | Es un ancla sobre el **parénquima**, que el gemelo controla con `SEMIEJES` (anclado) + seno |
| **D3** | Fijar el objetivo **médula = 8.8 mm** en el eje hiliar | Glodny, PW − CW = 15.4 − 6.6 | Derivado de D2 y del umbral ya anclado |
| **D4** | Adoptar `fraccion_huerfana` como criterio cuantitativo | criterio operativo, [SUPUESTO DECLARADO] | No necesita fuente nueva; es métrica interna |
| **D5** | Sustituir contención-en-cono por **partición** en `assign_pyramids` | Bonsib (lóbulo, Bertin) | Independiente de las dimensiones del seno. Da 0 % huérfana con cualquier seno |
| **D6** | Retirar `pyramid_r_base` y `cone_half_angle_deg` **si** se adopta D5 | 0 consumidores en el repo | Hecho comprobable hoy |
| **D7** | Corregir la contradicción documental papila-medular (Tarea 1b) | interno | Es coherencia entre capas, no morfometría |

**D2 y D3 son el hallazgo aprovechable de esta sesión:** el repo ya tenía todo para detectar que su
parénquima es 3.4 mm más delgado de lo que su propia fuente indica. No hacía falta Emamian.

### Decisiones BLOQUEADAS hasta tener los números del seno de Emamian 1993

| # | decisión | qué falta exactamente |
|---|---|---|
| **B1** | Anclar `SEMIEJES_SENO` | Longitud / anchura / grosor del área ecogénica central. Hoy es puesto 2 de los 10 pendientes (`09_paper_vascular/auditoria_correspondencia_anatomica.md:129`) |
| **B2** | Anclar `CENTRO_SENO` | Posición del área ecogénica dentro del contorno renal. Emamian da dimensiones; **[PENDIENTE DE ANCLA]** si da posición |
| **B3** | **Elegir entre las tres vías a PW 15.4** (cy=−30.6 / sy=19.4 / escala ×1.2125) | Las tres dan el mismo PW y **difieren en el volumen del seno**: 6.214 / 5.591 / 7.413 mL. **El volumen de Emamian es exactamente el desempate** |
| **B4** | Declarar el volumen de exclusión como magnitud anclada | Volumen del área ecogénica central. Cerraría ENTRADA 033, que hoy declara el seno **SUPUESTO NO ANCLADO** |
| **B5** | Decidir si la asimetría lobar (Opción F) exige seno asimétrico | Emamian mide un seno simétrico; **[PENDIENTE DE ANCLA]** cualquier asimetría supero-inferior |

**B3 es el nudo.** Las tres configuraciones son indistinguibles bajo Glodny —las tres dan PW_max
exactamente 15.4— y difieren en un 33 % de volumen de seno. **Elegir sin el volumen de Emamian sería
ajustar un parámetro libre a un solo dato y llamarlo anclaje.**

### Secuencia sugerida

```
AHORA (sin Emamian):
  D1 → registrar el defecto de los 4 ápices y su causa
  D2, D3 → fijar PW_max 15.4 y médula 8.8 como criterios de aceptación
  D4 → fijar fraccion_huerfana como métrica
  D5, D6 → decidir partición vs conos   [independiente del seno]
  D7 → corregir la contradicción documental

BLOQUEADO (requiere Emamian 1993 — dimensiones y volumen del área ecogénica central):
  B1, B2 → anclar SEMIEJES_SENO y CENTRO_SENO
  B3 → elegir entre cy=-30.6 / sy=19.4 / escala x1.2125
  B4 → cerrar ENTRADA 033
  B5 → asimetría lobar

DESPUÉS de B3:
  regenerar Capa 0 y la cascada 1→4 (9 .npz, orden en _auditoria/2026-08-07 Tarea 4)
```

**Nota de secuencia:** D5 (partición) **no depende** de B3. Puede decidirse y quedar escrito antes de
tener Emamian, porque da 0 % de huérfana con cualquier configuración de seno. Es la única decisión de
fondo que esta auditoría deja lista para tomar sin fuentes nuevas.

---

## Resumen

| tarea | resultado |
|---|---|
| **1a** | 4 ápices (k=0,1,8,9) a **5.295 / 6.428 mm** < 6.6. Causa: `uy = sqrt(1−ux²−uz²)` (`:339`) baja el ápice 2.78 mm hacia la cápsula mientras la cápsula sube 0.72 mm. `:335` sólo garantiza estar sobre el elipsoide, **no** la profundidad; **no existe test** |
| **1b** | Sí, como contradicción documental. `capa4_colector_alto.py:351,395,453` declara las papilas `interfaz medula/seno`; el texto ya está versionado en `auditoria_capa4_calicial.md:28`. `capa3c_colector.py:59-61` las trata como conducto de Bellini (medular). Ninguna falla numéricamente |
| **1c** | Dominante: **`CENTRO_SENO[1]`** (−32 → 0 ápices corticales). **`SEMIEJES_SENO` NO puede resolverlo**: ni ×1.3 baja de 2 |
| **2a** | PW: mediana **6.853**, p5 0.841, p95 11.531, **max 11.999 mm** |
| **2b** | **Más delgado, no más grueso.** 0.0 % dentro de [12.6, 18.2]; 100 % por debajo. Eje hiliar: **12.0 = 6.6 córtex + 5.4 médula** vs Glodny **15.4 = 6.6 + 8.8**. El déficit es **todo medular (−39 %)** |
| **2c** | Tres vías a PW_max 15.4: **sy = 19.4**, **cy = −30.6**, **escala ×1.2125**. La mediana no es el estadístico comparable; lo declaro |
| **2d** | Compatible por construcción. **Convergencia:** cy=−30.6 y sy=19.4 eliminan además los 4 ápices corticales. Coste: volumen excluido 4 → 5.6–7.4 mL, reabre ENTRADA 033 |
| **3** | Partición: **0 % huérfana por construcción**; ratio polar/central **1.878** (los conos dan 1.609) — reproduce Bonsib **sin calibrar**; `pyramid_r_base` y `cone_half_angle_deg` quedan sin sentido y **no los necesita nadie** (0 consumidores). Opera en `assign_pyramids` (`:356`), no en `build_pyramids` — **combinable con A–F** |
| **4** | **7 decisiones tomables ya** (D1–D7), **5 bloqueadas** (B1–B5). El nudo es **B3**: las tres vías a PW 15.4 son indistinguibles bajo Glodny y difieren 33 % en volumen de seno; **el volumen de Emamian es el desempate** |

### Lo que este informe cambia

1. **`SEMIEJES_SENO` deja de ser sospechoso principal.** No puede corregir los ápices y no es la
   variable dominante. El sospechoso es **`CENTRO_SENO[1] = −34`**, un centro situado 4 mm **fuera**
   del elipsoide principal.
2. **Un defecto nuevo, detectable con anclas que el repo ya tenía:** el parénquima del gemelo es
   **3.4 mm más delgado** que Glodny, y el déficit es enteramente medular (5.4 vs 8.8 mm).
3. **Aparece una opción que no estaba en A–F** y que no compite con ellas: partición en
   `assign_pyramids`, con 0 % de huérfana por construcción y el ratio polar/central emergente.
4. **La decisión de partición (D5) no está bloqueada por Emamian.** Es la única de fondo que puede
   tomarse hoy.

**Sin cambios en disco fuera de este archivo.** No se ejecutó `capa0_dominio.py`. No se hizo
`git add`, `git commit`, `git stash` ni `git restore`. Ningún `.py`, ningún `.npz` y ningún documento
existente fue modificado.
