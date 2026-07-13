#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
capa3b_shunts_vs_voxel.py  --  Bio-Kidney AI 2026
=================================================
AUDITORIA DE SOLO LECTURA. NO escribe ningun .npz, NO promueve nada, NO toca el
generador ni el visualizador. Solo imprime a stdout.

PREGUNTA: los shunts venoso-venoso (VV) del cuadrante critico, ¿funden lumenes
de VERDAD a la resolucion de voxelizacion prevista para la Capa 4, o caen por
debajo de resolucion (y entonces son un artefacto sub-vóxel, no un defecto real)?

Cruza la distribucion de PENETRACION de los shunts VV con la resolucion de
voxelizacion candidata de Capa 4.  Trabaja sobre:
  - capa3b_venoso.npz    (canonico actual; NO el reparado, NO los wrep)
  - capa3a_arterial.npz  (control: nivel ya aceptado)

CRITERIO DE "FUNDIR LUMENES" a un tamaño de vóxel v
---------------------------------------------------
Para cada shunt critico, la penetracion  p = thr - dd = (r_i + r_j) - dd  ES el
SOLAPAMIENTO FISICO de los dos tubos (ya considera los radios): es cuanto se
interpenetran los dos lumenes cilindricos cuyos ejes distan dd.
  - Si p >= v  -> los dos lumenes se solapan MAS que un vóxel: al voxelizar caen
    en el mismo vacio y se FUNDEN en un cortocircuito (shunt real a esa escala).
  - Si p <  v  -> el solape es sub-vóxel: caen en vóxeles separados o en la pared
    y NO se funden (artefacto que la voxelizacion borra).

CALIBRE (la columna critica)
----------------------------
Un shunt entre vasos de CALIBRE (vasos con recorrido y lumen grueso) seria un
cortocircuito hemodinamico real; un cruce entre venulas terminales finas es un
roce de ramitas sin consecuencia. Segmento "de calibre" = ambos nodos de grado
>= 2 (interno, no hoja) Y radio > 2x el radio terminal del arbol. Un shunt
"involucra calibre" si AL MENOS uno de sus dos segmentos es de calibre.

VOXELIZACION DE CAPA 4
----------------------
Capa 4 (calices/pelvis + voxelizacion adaptativa por octree) esta PENDIENTE: el
repo NO define aun un tamaño de vóxel objetivo unico (solo la aspiracion de
refinar a ~10 um en ramas terminales via octree; 10 um uniforme = inviable,
~150 GB). Por tanto se BARRE un rango de tamaños candidatos y se reporta como
SUPUESTO EXPLICITO.

Reutiliza clasificar() (que internamente usa seg_seg_dist) de
capa3b_clasificar_vv y la logica de cuadrante critico critico() de
capa3b_severidad_vv.

Solo NumPy / SciPy. SIN bpy. Entorno: env_biokidney.
"""

import os
import sys
import numpy as np

# ============================================================================
#  PARAMETROS
# ============================================================================
# Tamaños de vóxel candidatos (Capa 4 no define uno todavia -> SUPUESTO).
VOXELS_UM = [50, 100, 150, 200, 300]
PRINTABLE_UM = (100, 150)      # ventana de resolucion imprimible realista
CALIBRE_FACTOR = 2.0           # "de calibre": radio > CALIBRE_FACTOR x terminal

IN_VEN = "capa3b_venoso.npz"   # CANONICO actual (NO reparado, NO wrep)
IN_ART = "capa3a_arterial.npz" # control

# ============================================================================
#  REUTILIZACION (solo lectura)
# ============================================================================
_here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
if _here not in sys.path:
    sys.path.insert(0, _here)
from capa3b_clasificar_vv import clasificar, find_npz   # usa seg_seg_dist por dentro
from capa3b_severidad_vv import critico                  # cuadrante critico VV


def pct(n, tot):
    return f"{(100.0*n/tot):.3f}%" if tot else "  -  "


# ============================================================================
#  ANALISIS DE UN ARBOL
# ============================================================================
def analizar(npz):
    """Devuelve, para el cuadrante critico VV de un arbol:
       p (penetracion=solape, mm), involucra_calibre (bool), n_seg, terminal_r,
       calibre_thr, n_seg_calibre."""
    nodos = npz["nodos"].astype(np.float64)
    arist = npz["aristas"].astype(np.int64)
    radios = npz["radios"].astype(np.float64)

    cl = clasificar(nodos, arist, radios)
    cr = critico(cl)                      # i, j (idx de segmento), p = thr - dd, cont, n

    n_seg = len(arist)
    terminal_r = float(radios.min())     # los terminales son los mas finos (Murray)
    calibre_thr = CALIBRE_FACTOR * terminal_r

    # segmento de calibre: interno (ningun extremo de grado 1) Y radio grueso
    interno = ~cl["seg_terminal"]
    seg_calibre = interno & (radios > calibre_thr)
    n_seg_calibre = int(np.count_nonzero(seg_calibre))

    if cr["n"] > 0:
        involucra_calibre = seg_calibre[cr["i"]] | seg_calibre[cr["j"]]
        p = cr["p"]
    else:
        involucra_calibre = np.empty(0, bool)
        p = np.empty(0)

    return dict(p=p, calibre=involucra_calibre, n_crit=cr["n"], n_seg=n_seg,
                terminal_r=terminal_r, calibre_thr=calibre_thr,
                n_seg_calibre=n_seg_calibre)


def tabla_arbol(nombre, A):
    """Imprime la tabla por tamaño de vóxel de un arbol. Devuelve dict
       {v_um: (n_funde, n_funde_calibre)} para el resumen comparativo."""
    print(f"  ARBOL: {nombre}")
    print(f"    segmentos={A['n_seg']}   criticos VV (no-exacta ∩ >=1 interno)={A['n_crit']}   "
          f"r_terminal={A['terminal_r']*1000:.1f}um   calibre>{A['calibre_thr']*1000:.1f}um   "
          f"segmentos de calibre={A['n_seg_calibre']}")
    print("    " + "-" * 78)
    print(f"    {'vóxel v':>8} | {'shunts funden':>13} | {'% sobre seg':>11} | "
          f"{'de esos, >=1 calibre':>20}")
    print(f"    {'(um)':>8} | {'(p >= v)':>13} | {'(n/n_seg)':>11} | "
          f"{'(cortocircuito real)':>20}")
    print("    " + "-" * 78)
    res = {}
    for v_um in VOXELS_UM:
        v_mm = v_um / 1000.0
        funde = A["p"] >= v_mm
        n_funde = int(np.count_nonzero(funde))
        n_cal = int(np.count_nonzero(funde & A["calibre"])) if A["n_crit"] else 0
        res[v_um] = (n_funde, n_cal)
        print(f"    {v_um:>8d} | {n_funde:>13d} | {pct(n_funde, A['n_seg']):>11} | "
              f"{n_cal:>20d}")
    print("    " + "-" * 78)
    return res


# ============================================================================
#  MAIN
# ============================================================================
def main():
    ven_path = find_npz(IN_VEN)
    art_path = find_npz(IN_ART)
    ven = np.load(ven_path, allow_pickle=False)
    art = np.load(art_path, allow_pickle=False)

    print("=" * 84)
    print("  CAPA 3b - SHUNTS VV vs RESOLUCION DE VOXELIZACION DE CAPA 4 (solo lectura)")
    print("=" * 84)
    print(f"  Venoso (canonico) : {ven_path}")
    print(f"  Arterial (control): {art_path}")
    print(f"  p = thr - dd = (r_i+r_j) - dd  = SOLAPAMIENTO FISICO de los dos lumenes (mm).")
    print(f"  Funde a vóxel v  <=>  p >= v.   Calibre: interno Y radio > {CALIBRE_FACTOR:g}x terminal.")
    print("-" * 84)
    print("  SUPUESTO DE VOXELIZACION (Capa 4 PENDIENTE, sin tamaño de vóxel canonico):")
    print(f"     Se barre el rango candidato {VOXELS_UM} um.")
    print(f"     (Bitacora: voxelizacion adaptativa por octree, aspiracion ~10um SOLO en")
    print(f"      ramas terminales; 10um uniforme inviable ~150GB. Ventana imprimible ~{PRINTABLE_UM[0]}-{PRINTABLE_UM[1]}um.)")
    print("-" * 84)

    A_ven = analizar(ven)
    A_art = analizar(art)

    res_ven = tabla_arbol("VENOSO (canonico capa3b_venoso.npz)", A_ven)
    res_art = tabla_arbol("ARTERIAL (control capa3a_arterial.npz)", A_art)

    # ------------------------------------------------------------------
    #  COMPARATIVA COMPACTA venoso vs arterial
    # ------------------------------------------------------------------
    print("  COMPARATIVA VENOSO vs ARTERIAL  (shunts que FUNDEN lumenes)")
    print("  " + "-" * 80)
    print(f"  {'vóxel':>6} | {'VEN funde':>9} {'VEN calibre':>11} | "
          f"{'ART funde':>9} {'ART calibre':>11} | {'ratio calibre V/A':>17}")
    print("  " + "-" * 80)
    for v_um in VOXELS_UM:
        nv, nvc = res_ven[v_um]
        na, nac = res_art[v_um]
        if nac > 0:
            ratio = f"{nvc/nac:>13.1f}x"
        elif nvc > 0:
            ratio = f"{'inf (A=0)':>14}"
        else:
            ratio = f"{'-- (0/0)':>14}"
        print(f"  {v_um:>6d} | {nv:>9d} {nvc:>11d} | {na:>9d} {nac:>11d} | {ratio:>17}")
    print("  " + "-" * 80)
    print("-" * 84)

    # ------------------------------------------------------------------
    #  VEREDICTO (una linea por escala + conclusion)
    # ------------------------------------------------------------------
    print("  VEREDICTO -- shunts de CALIBRE que funden lumenes (venoso vs arterial control)")
    print("  " + "-" * 80)
    for v_um in VOXELS_UM:
        nvc = res_ven[v_um][1]
        nac = res_art[v_um][1]
        marca = "  <-- ventana imprimible" if v_um in PRINTABLE_UM else ""
        print(f"  a v={v_um:>3d}um : VENOSO {nvc:>5d} cortocircuitos de calibre  |  "
              f"ARTERIAL {nac:>5d} (control){marca}")
    print("  " + "-" * 80)

    # conclusion data-driven sobre la ventana imprimible
    ven_win = [res_ven[v][1] for v in PRINTABLE_UM]
    art_win = [res_art[v][1] for v in PRINTABLE_UM]
    ven_max = max(ven_win); art_ref = max(art_win)
    print("  CONCLUSION")
    if art_ref == 0 and ven_max == 0:
        print(f"  A resolucion imprimible ({PRINTABLE_UM[0]}-{PRINTABLE_UM[1]}um) NI el venoso NI el arterial "
              "funden lumenes de calibre.")
        print("  -> Los shunts VV son artefactos SUB-VOXEL: la voxelizacion de Capa 4 los borra. "
              "VV NO es un defecto real a la escala de impresion.")
    else:
        if art_ref == 0:
            orden = float('inf')
            rel = "el arterial control no funde ninguno de calibre"
        else:
            orden = ven_max / art_ref
            rel = f"~{orden:.1f}x el arterial control"
        if orden <= 3.0:
            print(f"  A resolucion imprimible ({PRINTABLE_UM[0]}-{PRINTABLE_UM[1]}um) el venoso funde "
                  f"{ven_max} cortocircuitos de calibre, {rel} ({art_ref}).")
            print("  -> Mismo ORDEN que el nivel arterial ya aceptado: VV es una LIMITACION DE ESCALA "
                  "consistente con lo ya documentado, no un defecto nuevo.")
        else:
            print(f"  A resolucion imprimible ({PRINTABLE_UM[0]}-{PRINTABLE_UM[1]}um) el venoso funde "
                  f"{ven_max} cortocircuitos de calibre, {rel} ({art_ref}).")
            print("  -> EXCEDE al arterial por encima de su orden incluso filtrando por calibre: hay un "
                  "PROBLEMA ESTRUCTURAL real que SOBREVIVE a la voxelizacion. Discutir arquitectura.")
    print("=" * 84)


if __name__ == "__main__":
    main()
