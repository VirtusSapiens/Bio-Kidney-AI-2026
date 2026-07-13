# Diagnostico estanqueidad Capa 5a (solo lectura)
**Programa:** Bio-Kidney AI 2026 - **Capa:** 5a (diagnostico) - **Fecha:** 2026-07-09
Localiza y clasifica los voxeles compartidos urinario<->vascular de la auditoria 5a. Reusa centerlines/radios ya normalizados (mismo adaptador y clamp de 5a). NO re-malla, NO modifica 5a, NO genera 5b.
## (1) Voxeles compartidos
1118 voxeles urinario & vascular (5a reporto 1118).
## (2)-(3) Clasificacion nivel_urinario x origen_vascular
| nivel urinario | vascular | both_native | any_clamped |
|---|---|---|---|
| pelvis | venoso | 825 | 0 |
| ureter | venoso | 134 | 0 |
| pelvis | arterial | 69 | 0 |
| caliz_mayor | venoso | 34 | 0 |
| ureter | arterial | 31 | 0 |
| colector_3c | venoso | 22 | 3 |

- Voxeles con **>=1 lado clampeado**: **3** ; **ambos nativos**: **1115**.
- (urin clampeado en 0, vascular clampeado en 3 voxeles)
## (4) Pared inter-circuito NATIVA (radios originales) por voxel
Distingue sub-muestreo (pared real >=400um, solo el voxel 200um la puentea) de colision real (<400um).

| pared nativa | voxeles |
|---|---|
| <0 (solape real) | 1067 |
| 0-100um | 0 |
| 100-400um | 51 |
| >=400um (sub-muestreo) | 0 |

- Min pared CLAMPEADA (geometria actual de 5a): **-4347.3 um** en [ -1.12 -27.56   0.  ] (urin caliz_mayor/pelvis vs venoso).
- Min pared NATIVA: **-4347.3 um**.
## (5) Zoom hilio (pelvis/ureter vs tronco venoso, Y en [-31,-27])
- Pared lumen-lumen minima pelvis/ureter <-> venoso: **-4398.9 um** (<400um: HAY componente de hilio)
- urin ~[  1.95 -25.81   0.  ]  venoso ~[ -0.6  -27.28  -0.16] (3 segs pelvis/ureter x 6 segs venosos en banda)
## (6) Genero / handle de la malla urinaria
Fusiones urinarias no-adyacentes (lumen clampeado solapa -> crean asa): **2**.
- edge861 (papila_junction/caliz_menor) x edge869 (caliz_menor/caliz_mayor): d=3.168 solape=0.085 mm en [ -6.08 -20.92   4.05]
- edge866 (papila_junction/caliz_menor) x edge878 (caliz_menor/caliz_mayor): d=3.168 solape=0.085 mm en [  6.08 -20.92  -4.05]
Coincide con la fusion co-mayor (calices menores del mismo mayor, solape 85um de la Entrada 015).
## (7) Numeros que quedaron cortados en 5a
- Porosidad: **0.7418 %**
- Lumen urinaria: 878.26 mm3
- Lumen vascular_0_arterial: 1.45 mm3
- Lumen vascular_1_venoso: 43.18 mm3
- Gap arterial<->venoso (nodos): **1.546 mm**
## VEREDICTO
- **Sub-muestreo** (pared nativa >=400um, se separa a 100um/5b): **0/1118 (0.0%)**
- **Colision real** (pared nativa <400um, requiere fix de geometria): **1118/1118 (100.0%)**
- **Componente de hilio (pelvis/ureter vs venoso <400um):** **SI** (pared -4398.9 um).
