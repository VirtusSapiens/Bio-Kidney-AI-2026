# Diagnostico de holgura de la pelvis (solo lectura)
**Programa:** Bio-Kidney AI 2026 - **Capa:** 5a (diagnostico) - **Fecha:** 2026-07-10
Barrido de posicion del centro de la pelvis (X=0, Y-Z) buscando cierre simultaneo de pared vs venoso/arterial, contencion en seno, intrarrenalidad, infundibulos, ureter y corredor anterior. Reusa SDF/adaptador/clamp/FIX1 (pelvis elipsoide [4.5,4.5,2.0]). NADA se modifica.

```
DIAGNOSTICO DE HOLGURA DE LA PELVIS (solo lectura)
venoso cones 472 (hilio 306), arterial cones 29 (hilio 29)

===== D1  MAPA DEL HILIO (venoso/arterial peri-pelvico) =====
  venoso: 32 segmentos en la region; anterior(+Z) 5, posterior(-Z) 27, Z~0 1
  arterial: 29 segmentos en la region
  HUECO inter-tributaria venoso (X=0): clearance max 10.83 mm en Y=-20.2 Z=9.0
  vacio venoso (clr>=0.5mm) extent: Y[-33.0,-20.2] Z[-9.0,9.0]

===== D2  BARRIDO DE CENTRO DE PELVIS (X=0) =====
  centros probados: 525  ->  VIABLES (pared_ven>=0.4 Y pared_art>=0.4 Y seno<=1): 17
  mascara viable: Y[-30.0,-23.0]  Z[-6.0,2.5]
     C=[0,-23.0,2.0] pared_ven 0.54 pared_art 0.87 seno 0.98 intrarrenal 100%
     C=[0,-23.0,2.5] pared_ven 0.85 pared_art 0.59 seno 1.00 intrarrenal 100%
     C=[0,-28.0,-6.0] pared_ven 0.62 pared_art 3.96 seno 0.79 intrarrenal 47%
     C=[0,-28.5,-6.0] pared_ven 0.93 pared_art 3.93 seno 0.75 intrarrenal 42%
     C=[0,-28.5,-5.5] pared_ven 0.60 pared_art 3.44 seno 0.69 intrarrenal 42%
     C=[0,-29.0,-6.0] pared_ven 1.24 pared_art 3.84 seno 0.72 intrarrenal 35%
     C=[0,-29.0,-5.5] pared_ven 0.88 pared_art 3.35 seno 0.66 intrarrenal 35%
     C=[0,-29.0,-5.0] pared_ven 0.56 pared_art 2.86 seno 0.60 intrarrenal 35%
     C=[0,-29.5,-6.0] pared_ven 1.54 pared_art 3.80 seno 0.68 intrarrenal 31%
     C=[0,-29.5,-5.5] pared_ven 1.21 pared_art 3.30 seno 0.63 intrarrenal 31%
     C=[0,-29.5,-5.0] pared_ven 0.88 pared_art 2.80 seno 0.57 intrarrenal 31%
     C=[0,-29.5,-4.5] pared_ven 0.52 pared_art 2.30 seno 0.52 intrarrenal 31%

===== D3  INFUNDIBULOS MAYORES / D4 URETER / D5 CORREDOR (sobre viables) =====
  C=[0,-23.0,2.0]: D3 inf_wall -0.34 swing 55deg [X] | D4 ureter_wall -0.96 [X] | D5 puertos(v2.11,a2.22) [OK] -> no
  C=[0,-23.0,2.5]: D3 inf_wall -0.09 swing 61deg [X] | D4 ureter_wall -0.96 [X] | D5 puertos(v1.97,a2.27) [OK] -> no
  C=[0,-28.0,-6.0]: D3 inf_wall -2.29 swing 59deg [X] | D4 ureter_wall -0.96 [X] | D5 puertos(v7.69,a6.67) [OK] -> no
  C=[0,-28.5,-6.0]: D3 inf_wall -2.14 swing 56deg [X] | D4 ureter_wall -0.96 [X] | D5 puertos(v7.66,a6.64) [OK] -> no
  C=[0,-28.5,-5.5]: D3 inf_wall -2.25 swing 53deg [X] | D4 ureter_wall -0.96 [X] | D5 puertos(v7.17,a6.14) [OK] -> no
  C=[0,-29.0,-6.0]: D3 inf_wall -1.99 swing 53deg [X] | D4 ureter_wall -0.96 [X] | D5 puertos(v7.64,a6.62) [OK] -> no
  C=[0,-29.0,-5.5]: D3 inf_wall -2.09 swing 50deg [X] | D4 ureter_wall -0.96 [X] | D5 puertos(v7.14,a6.12) [OK] -> no
  C=[0,-29.0,-5.0]: D3 inf_wall -2.20 swing 47deg [X] | D4 ureter_wall -0.96 [X] | D5 puertos(v6.65,a5.62) [OK] -> no
  C=[0,-29.5,-6.0]: D3 inf_wall -1.85 swing 50deg [X] | D4 ureter_wall -0.96 [X] | D5 puertos(v7.63,a6.60) [OK] -> no
  C=[0,-29.5,-5.5]: D3 inf_wall -1.94 swing 47deg [X] | D4 ureter_wall -0.96 [X] | D5 puertos(v7.13,a6.10) [OK] -> no
  C=[0,-29.5,-5.0]: D3 inf_wall -2.04 swing 44deg [X] | D4 ureter_wall -0.96 [X] | D5 puertos(v6.63,a5.61) [OK] -> no
  C=[0,-29.5,-4.5]: D3 inf_wall -2.15 swing 41deg [X] | D4 ureter_wall -0.96 [X] | D5 puertos(v6.13,a5.11) [OK] -> no
  C=[0,-30.0,-6.0]: D3 inf_wall -1.72 swing 47deg [X] | D4 ureter_wall -0.96 [X] | D5 puertos(v7.63,a6.60) [OK] -> no
  C=[0,-30.0,-5.5]: D3 inf_wall -1.81 swing 45deg [X] | D4 ureter_wall -0.96 [X] | D5 puertos(v7.13,a6.10) [OK] -> no
  C=[0,-30.0,-5.0]: D3 inf_wall -1.89 swing 42deg [X] | D4 ureter_wall -0.96 [X] | D5 puertos(v6.63,a5.60) [OK] -> no

===== D6  VEREDICTO =====
  NINGUN centro D2-viable cierra D3-D5. Fallos entre 15 candidatos: D3(infundibulos) 15/15, D4(ureter) 15/15, D5(corredor) 0/15
  El cuerpo de la pelvis SI cabe (D2: 17 centros viables), pero los INFUNDIBULOS
  (calices en Z~0 -> pelvis) y el URETER (pelvis -> hilio [0,-30,0]) cruzan el plexo venoso en Z~0.
  NINGUNA posicion (barrido completo 525 centros) tiene infundibulos+ureter limpios.
  Mejor posicion posible [0,-23.0,1.0]: mejor pared infundibulo/ureter = -964um (faltan 1364um).
  Los infundibulos (calices Z~0 -> pelvis) y el ureter (pelvis -> hilio Z~0) estan ANCLADOS en Z~0
  e INDEPENDIENTES del tamaño de la pelvis: reducir/mover la pelvis NO los despeja.
  BINDING = el PLEXO VENOSO peri-hilar (32 segs, 27 posteriores) ocupa el corredor Z~0 obligado.
  -> el fix NO es posicion NI tamaño de la pelvis: es SEGREGAR/mover el plexo venoso (revisar 3b).
  Confirma cuantitativamente el blocker de las Entradas 019-021 (arbol venoso ramificado en el hilio).
```

## Veredicto (estructurado)
```json
{
  "status": "NO_ES_LA_PELVIS__BINDING_PLEXO_VENOSO_Z0",
  "mejor_pared_inf_ureter_um": -964.3,
  "d3_fail": 15,
  "d4_fail": 15,
  "d2_viables": 17
}
```
