# Auditoria Capa 5a - nucleo SDF (validacion gruesa 200um)
**Programa:** Bio-Kidney AI 2026 - **Capa:** 5a - **Fecha:** 2026-07-09
Cada capa = centerline+radio = union de capsulas conicas = SDF (rounded-cone). Solido = elipsoide COMPLETO de Capa 0 (seno relleno) menos union de lumenes. Muestreo grueso a 200um (smoke); el paso fino 100um y el voxel etiquetado son 5b.
## PASO 0 - esquemas normalizados (nodos, aristas, radios/nodo)
| capa | nodos | aristas | radios | adaptador |
|---|---|---|---|---|
| 3a arterial | 12016 | 12015 | por-ARISTA -> nodo(child); raiz=214um | adapt_edge_radios |
| 3ab peritubular | 4290 pts | 0 | ninguno (bed capilar, sub-resolucion) | puntos edge-less |
| 3b venoso | 23799 | 23798 | por-ARISTA -> nodo(child); raiz=524um | adapt_edge_radios |
| 3c colector | 2189 | 2179 | por-NODO; topologia via parent | adapt_3c |
| 4 calicial | 24 | 23 | por-NODO; aristas explicitas | adapt_4 |
## PASO 1 - circuitos
- **URINARIO = 3c U 4:** 492 nodos, 2202 aristas, **1 componente(s)** [OK] una sola (papilas fusionadas)
- **VASCULAR = 3a U 3ab U 3b:** 26025 nodos, 35813 aristas -> **2 componentes-tubo** (arterial comp 0, venoso comp 4291) + 4290 puntos aislados (3ab). Arterial y venoso **SEPARADOS (gap capilar de diseño, dato no defecto)**.
## PASO 2 - poda + clamp (piso 200um r / 400um diam)
| circuito | aristas_tot | podadas | retenidas | frac% superv | clamp | reconx |
|---|---|---|---|---|---|---|
| urinaria | 2202 | 1320 | 882 | 40.1% | 10 | 0 |
| vascular_0_arterial | 12015 | 11986 | 29 | 0.2% | 1 | 0 |
| vascular_1_venoso | 23798 | 23299 | 499 | 2.1% | 77 | 0 |

> La microvasculatura sub-400um **NO se imprime**: queda en 3a/3ab/3b/3c como **OBJETIVO DE PERFUSION** (maduracion capilar wet-lab), no como canal impreso.
## PASO 4 - auditorias
### (1) Watertight = cierre 2-manifold (malla cruda pre-clip; el clip abre los puertos a proposito)
Criterio: manifold + 0 bordes + 0 non-manifold. euler=2 solo si genero 0; el genero>0 (asas por fusion de ramas anidadas a 200um) NO rompe el cierre.

| circuito | V | E | F | manifold | bordes | nonmanif | euler | genero | estado |
|---|---|---|---|---|---|---|---|---|---|
| urinaria | 37072 | 111216 | 74144 | True | 0 | 0 | 0 | 1 | [OK] watertight |
| vascular_0_arterial | 508 | 1518 | 1012 | True | 0 | 0 | 2 | 0 | [OK] watertight |
| vascular_1_venoso | 10042 | 30120 | 20080 | True | 0 | 0 | 2 | 0 | [OK] watertight |

### (2) Estanqueidad
- Voxeles dentro de AMBOS (urinario & vascular): **1118** (urin&arterial=100, urin&venoso=1018)
  - **[DATO EMERGENTE]** El solape NO viene de centerlines cruzados (siguen disjuntos) sino del **clamp-a-piso**: micro-vasos venosos sub-400um se inflan a 400um diam y, al muestrear a 200um, su lumen fusiona con el colector urinario en la medula (interdigitacion vasa-recta / conducto colector). Se resuelve en 5b (100um + politica de clamp / separacion inter-circuito >=400um por diseño). Es tension de escala del modelo reducido, no cruce topologico.
- Distancia minima superficie-superficie urinario<->vascular: **0.025 mm**
- Gap arterial<->venoso (nodos retenidos): **1.546 mm** (gap capilar de diseño)
### (3) Conectividad a puerto
- urinaria: puerto a 0.000 mm, alcanzado por el circuito retenido: [OK]
- vascular_0_arterial: puerto a 0.000 mm, alcanzado por el circuito retenido: [OK]
- vascular_1_venoso: puerto a 0.489 mm, alcanzado por el circuito retenido: [OK]
### (4) Contencion en elipsoide
- urinaria: violaciones 0/408 (val max 1.000) [OK]
- vascular_0_arterial: violaciones 0/18 (val max 1.000) [OK]
- vascular_1_venoso: violaciones 1/316 (val max 1.033) -> puertos (salida por diseño)
### (5) Volumenes / porosidad
- Matriz (elipsoide completo, seno relleno): **124407.1 mm3**
- Lumen urinaria: 878.26 mm3
- Lumen vascular_0_arterial: 1.45 mm3
- Lumen vascular_1_venoso: 43.18 mm3
- Lumen TOTAL: **922.88 mm3** - **POROSIDAD: 0.7418 %**

## Estado
**Capa 5a cerrada.** Nucleo SDF muestreado a 200um, poda a 400um diam, lumenes urinario/arterial/venoso mallados (STL binario propio), auditorias corridas. Tiempo 21.9s, pico memoria 287 MB. **5b (matriz fina 100um + voxel etiquetado) PENDIENTE.**
