# Auditoria fix estanqueidad Capa 5a
**Programa:** Bio-Kidney AI 2026 - **Capa:** 5a (fix) - **Fecha:** 2026-07-09
FIX1 pelvis elipsoide aplanado (Capa 4) + FIX2 reruteo venoso +Z + FIX3 arterial -Z (local en Capa 5, sin reescribir 3a/3b/3c). MIN_WALL=0.400 mm.

```
FIX ESTANQUEIDAD 5a (pelvis elipsoide + reruteo venoso/arterial)
FIX1 pelvis elipsoide [4.5 4.5 2. ] en [  0.   -27.21   0.  ]  (Z_pelvis=2.0)
URINARIO: 492 nodos, 882 aristas retenidas, pelvis nodo 490, keep-out cones 3
REROUTE venoso (+1Z): 9 nodos a mover  max 4.99mm  mediana 3.38mm  puerto 0.00mm  infeasibles(>5mm): 0
   aplicado: 19 nodos movidos, max cumul 5.00mm, puerto 0.08mm, restan violando 4
REROUTE arterial (-1Z): 5 nodos a mover  max 4.53mm  mediana 3.41mm  puerto 2.23mm  infeasibles(>5mm): 0
   aplicado: 16 nodos movidos, max cumul 5.00mm, puerto 2.23mm, restan violando 7

[FALLBACK] la reruteo de direccion-fija NO logra MIN_WALL con movimientos <=5mm.
CAUSA GEOMETRICA: los arboles vascular(es) CRUZAN la pelvis por AMBOS lados (Z+ y Z-);
empujar todo un arbol en un solo sentido obliga a nodos anteriores a atravesar la pelvis.
   venoso: tras aplicar (cap 5mm) restan 4 nodos violando MIN_WALL (puerto movido 0.08mm, max cumul 5.00mm)
      STRADDLE: de los 4 residuales, 0 ANTERIORES (Z+) y 4 POSTERIORES (Z-) de la pelvis -> empujarlos todos a un lado obliga a cruzarla.
      nodo 11 orig [ -0.96 -25.04  -2.74] (Z=-2.74): pared final 11um (faltan 389um); a pelvis-elipsoide 0.98mm
      nodo 12 orig [ -0.96 -25.04  -2.74] (Z=-2.74): pared final 335um (faltan 65um); a pelvis-elipsoide 0.98mm
      nodo 13 orig [ -1.1  -24.52  -3.25] (Z=-3.25): pared final -314um (faltan 714um); a pelvis-elipsoide 1.58mm
      nodo 14 orig [ -1.1  -24.52  -3.25] (Z=-3.25): pared final 9um (faltan 391um); a pelvis-elipsoide 1.58mm
   arterial: tras aplicar (cap 5mm) restan 7 nodos violando MIN_WALL (puerto movido 2.23mm, max cumul 5.00mm)
      STRADDLE: de los 7 residuales, 7 ANTERIORES (Z+) y 0 POSTERIORES (Z-) de la pelvis -> empujarlos todos a un lado obliga a cruzarla.
      nodo 5 orig [  0.36 -28.45   2.54] (Z=+2.54): pared final 322um (faltan 78um); a pelvis-elipsoide 0.62mm
      nodo 6 orig [  0.44 -28.14   3.05] (Z=+3.05): pared final -211um (faltan 611um); a pelvis-elipsoide 1.09mm
      nodo 7 orig [  0.44 -28.14   3.05] (Z=+3.05): pared final -196um (faltan 596um); a pelvis-elipsoide 1.09mm
      nodo 8 orig [  0.52 -27.8    3.54] (Z=+3.54): pared final -724um (faltan 1124um); a pelvis-elipsoide 1.56mm

   OPCIONES (NO aplico ninguna; las decidis vos):
   (a) mas aplanamiento de pelvis (Z_pelvis<2.0) -> reduce el keep-out en Z.
   (b) mover el nodo pelvis en -Y (hundir hacia el seno) -> la aleja del cruce vascular.
   (c) reruteo POR-GRADIENTE (cada nodo se aparta de la superficie mas cercana, no direccion fija).
   (d) carve keep-out en la matriz (NO por defecto; preferis decidirlo).
   FIX1 (pelvis elipsoide en capa4) QUEDA aplicado y verificado (byte-identico); solo el reruteo espera decision.
```
