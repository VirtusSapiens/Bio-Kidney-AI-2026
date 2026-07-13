# Auditoria Capa 4 — sistema calicial alto
**Programa:** Bio-Kidney AI 2026 · **Capa:** 4 (colector alto) · **Fecha:** 2026-07-06
Representacion centerline+radio, convergente papila -> ureter. Radios ABSOLUTOS (mm). Bicalicial (10 calices menores -> 2 mayores -> pelvis -> ureter). NO voxeliza (Capa 5).
- Nodos: **24** · Aristas: **23** (arbol: E=M-1=23)
## Verificacion previa
- Rango de `pyramid_apex` por eje (X,Y,Z): [23.760, 2.784, 8.800] mm -> eje de mayor rango = **X** (EJE_POLAR=X): [OK]
- `papila_nodo_ref` vs `pyramid_apex`: desvio max **7.18e-07 mm** [OK]
## (1) Alcanzabilidad
Desde cada papila_junction, siguiendo aristas, se llega al unico nodo ureter: **10/10** [OK]
## (2) Colisiones (aristas no adyacentes)
- **Centerline (test de cruce):** dist minima entre centerlines no adyacentes = **3.168 mm** (>0 -> las ramas NO se cruzan) [OK]
- **Luz (solape de lumen):** radio INTERPOLADO en el punto de aproximacion (r(s)=r_parent+s*(r_child-r_parent)); flag si dist < r_i(s)+r_j(t). Test fisico correcto para segmentos con taper (evita el sobre-conteo de usar el radio grande de la pelvis a lo largo de toda su arista).
  - **Mayores DISTINTOS (CRITERIO DE FALLO, meta 0):** **0** [OK] — ramas de calices mayores diferentes NO deben solaparse.
  - **Mismo mayor (embudo esperado):** 2 [NOTA, no defecto]: calices menores del mismo caliz mayor que convergen hacia su infundibulo comun (near-tangencia; magnitud del solape sub-0.1 mm).
    - arista 2 (2, 12) vs 10 (10, 20) (mayor 0): d=3.168 mm, solape 0.085 mm
    - arista 7 (7, 17) vs 19 (19, 21) (mayor 1): d=3.168 mm, solape 0.085 mm
  - **Con la pelvis:** 20 [NOTA, no defecto]: la pelvis (r=4.5 mm) es la camara donde drenan los calices; su lumen envuelve las bocas de los infundibulos que la alimentan -> continuidad anatomica calices->pelvis.
## (3) Contencion en el seno renal
Test elipsoide: ((x)/22)^2 + ((y+34)/16)^2 + ((z)/11)^2 <= 1.0

| grupo | dentro/total | val max | estado |
|---|---|---|---|
| caliz_menor | 10/10 | 0.794 | [OK] |
| caliz_mayor | 2/2 | 0.391 | [OK] |
| pelvis | 1/1 | 0.180 | [OK] |
| ureter | 1/1 | 0.062 | [OK] |

`papila_junction` (10) es nodo-interfaz medula/seno, **EXENTA** del test (se sienta sobre la pared +Y del seno). Valor del elipsoide en las papilas: [1.000, 1.000]; polares k=0=1.000, k=9=1.000 (~1.0 = sobre la pared del seno, esperado).
## (4) Monotonia de radio por rama
Perfil papila->ureter (mm): **[0.2, 1.5, 2.0, 4.5, 1.5]**. Creciente papila<copa<mayor<pelvis: [OK]. Caida unica pelvis->ureter (UPJ 4.5->1.5): [OK]. El escalon Bellini->copa (0.2->1.5) es junction DISCRETA (arista_tipo=1), no taper.
## Estado
**[OK] Capa 4 valida**: alcanzabilidad, colisiones, contencion y monotonia. Radios absolutos. Frontera: Capa 5 (voxelizacion), stitching via `papila_nodo_ref`.
