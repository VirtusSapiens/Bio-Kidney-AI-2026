# Auditoria fix estanqueidad Capa 5a v3 (pediculos hilares segregados)
**Programa:** Bio-Kidney AI 2026 - **Capa:** 5a (fix v3) - **Fecha:** 2026-07-09
Construccion deterministica de pediculos hilares con anclaje anatomico (orden AP vena-arteria-pelvis; +Z anterior). Mantiene FIX1 (pelvis elipsoide). Local en Capa 5, NO reescribe 3a/3b/3c, NO cambia Capa 4. MIN_WALL=0.400, MARGIN=0.200.

```
FIX v3 -- PEDICULOS HILARES SEGREGADOS (vena-arteria-pelvis, AP=Z, +Z anterior)

PASO1 CORREDORES Z (AP; +Z anterior):
   r_art=214um  r_ven=524um
   Z_art = 2.814 mm (arteria, 1er corredor anterior a pelvis)
   Z_ven = 4.152 mm (vena, anterior a la arteria)
   pelvis Z_semi=2.0; vena despeja pelvis directo: 1.629 >= 0.4 [OK]
   puertos nuevos: ARTERIA=[0,-30,2.814]  VENA=[0,-30,4.152]  (ureter=[0,-30,0] SIN cambios)

PASO2 CONFLICTO (keepout = pelvis elipsoide + ureter + calices, Capa 4):
   arterial: conflicto 5 nodos (Y[-30.0,-28.8]), puerto EN conflicto, salidas 1, nodos-rama(grado>2) 0  -> [OK] cadena terminal
   venoso: conflicto 10 nodos (Y[-29.1,-22.0]), puerto NO en conflicto, salidas 6, nodos-rama(grado>2) 1  -> [FALLA cadena]

   PEDICULO ARTERIAL (conflicto SI es cadena, PASO2 [OK]): ancestro nodo 28 en [  0.45 -24.12   8.22] (Y=-24.1 <=-22.7 [NO hay ancestro limpio anterior a la pelvis])
      pediculo [  0.45 -24.12   8.22] -> W[0,-22.7,2.81] -> puerto[0,-30,2.81]
      pared minima del pediculo vs keepout+venoso: -166um [<400um -> el tramo ancestro->W roza la pelvis]
      -> ARTERIAL: el stub imprimible retenido esta ENTERO en la banda Y de la pelvis (Y<=-22.7); no hay ancestro limpio anterior donde anclar el pediculo -> tampoco cierra.

[FALLBACK - PASO 2] al menos un circuito NO tiene el conflicto como CADENA TERMINAL con puerto.
   venoso: conflicto de 10 nodos con 1 nodo(s)-rama (grado>2) y 6 salidas; puerto NO en conflicto.
      => NO es una punta unica: es un PLEXO peri-hilar ramificado. Snipearlo y reconectar por UN pediculo desconectaria 5 ramas. La sintesis de pediculo-unico NO aplica.
      nodos-rama: [8] ; salidas (nodo,Y,Z): [(1, -29.6, 2.6), (11, -25.0, -2.7), (12, -25.0, -2.7), (21, -22.5, -4.8), (25, -21.5, -5.3), (26, -21.5, -5.3)]
      peor nodo 5 en [ -0.63 -27.53   0.14]: pared -3421um (faltan 3821um) vs urinario macro.

   ARTERIAL: conflicto SI es cadena (PASO2 [OK]) pero el pediculo NO verifica (stub imprimible entero en la banda Y de la pelvis, sin ancestro limpio anterior).
   OBSTACULO para cerrar estanqueidad: (1) el arbol VENOSO ramifica en el hilio
   (plexo peri-pelvico, no un tronco unico) -> la segregacion de pediculo-unico no lo captura;
   (2) el arterial imprimible es un stub sub-mm peri-pelvico sin parenquima anterior donde anclar.
   OPCIONES (NO aplico nada; decidis vos): (a) tronco venoso unico en el hilio (revisar 3b, fuera de
   alcance local); (b) MULTI-pediculo venoso (un corredor por cada una de las ramas en conflicto);
   (c) carve local del keepout; (d) mover pelvis en -Y para reducir el solape peri-pelvico venoso.
```
