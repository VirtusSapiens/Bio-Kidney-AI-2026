# Auditoria fix estanqueidad Capa 5a v2 (reruteo por-gradiente)
**Programa:** Bio-Kidney AI 2026 - **Capa:** 5a (fix v2) - **Fecha:** 2026-07-09
Reruteo POR-GRADIENTE (cada nodo se aparta de la superficie de otro circuito mas cercana, no direccion fija) manteniendo FIX1 (pelvis elipsoide). Local en Capa 5, NO reescribe 3a/3b/3c. MIN_WALL=0.400 mm, MAX_DESPL=5.0 mm.

```
FIX v2 ESTANQUEIDAD 5a -- reruteo POR-GRADIENTE (mantiene FIX1 pelvis elipsoide)
FIX1 (sin cambios): pelvis elipsoide [4.5 4.5 2. ] en [  0.   -27.21   0.  ]
PASO1 keep-out: venoso 7 nodos en conflicto; arterial 5 nodos en conflicto
PASO2 solver (8 iter, 262.4s):
   venoso  : 31 nodos movidos, despl max 1.32mm medio 0.03mm puerto 0.00mm; residual 8
   arterial: 5 nodos movidos, despl max 0.78mm medio 0.07mm puerto 0.41mm; residual 4

[FALLBACK] el reruteo por-gradiente NO convergio a MIN_WALL para todos los nodos.
   venoso nodo 4 en [ -0.91 -28.1    1.65]: pared -788um (faltan 1188um) vs urinario(pelvis/ureter/calices); grad estancado (atrapado entre superficies).
   venoso nodo 5 en [ -0.6  -27.52   0.55]: pared -2041um (faltan 2441um) vs urinario(pelvis/ureter/calices); grad estancado (atrapado entre superficies).
   venoso nodo 6 en [ -0.56 -27.04  -0.52]: pared -2072um (faltan 2472um) vs urinario(pelvis/ureter/calices); grad estancado (atrapado entre superficies).
   venoso nodo 7 en [ -0.62 -26.48  -1.52]: pared -985um (faltan 1385um) vs urinario(pelvis/ureter/calices); grad estancado (atrapado entre superficies).
   venoso nodo 8 en [ -0.63 -25.84  -2.47]: pared 54um (faltan 346um) vs urinario(pelvis/ureter/calices); grad estancado (atrapado entre superficies).
   arterial nodo 0 en [  0.   -30.41   0.  ]: pared -1513um (faltan 1913um) vs urinario(pelvis/ureter/calices); grad estancado (atrapado entre superficies).
   arterial nodo 1 en [  0.07 -29.71   0.56]: pared -1475um (faltan 1875um) vs urinario(pelvis/ureter/calices); grad estancado (atrapado entre superficies).
   arterial nodo 2 en [  0.15 -29.48   1.43]: pared -509um (faltan 909um) vs urinario(pelvis/ureter/calices); grad estancado (atrapado entre superficies).
   arterial nodo 3 en [  0.24 -29.23   2.28]: pared 270um (faltan 130um) vs urinario(pelvis/ureter/calices); grad estancado (atrapado entre superficies).
   (decision del usuario: carve local vs revisar el hilio; NO aplico nada)
```
