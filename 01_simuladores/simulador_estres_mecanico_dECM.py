#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BIO-KIDNEY AI 2026 - VirtusSapiens
Carlos David Moreno Caceres
SIMULADOR DE ESTRES MECANICO EN dECM - Modulo 09 v5

FISICA IMPLEMENTADA:
  1. Elasticidad lineal   : sigma = E * epsilon
  2. Kelvin-Voigt         : epsilon(t) = (sigma/E)*(1 - exp(-E*t/eta))
  3. Hagen-Poiseuille     : tau_wall = 4*eta*Q / (pi*r^3)
  4. Von Mises (Lame)     : cilindro pared gruesa, presion de poro externa
  5. Deformacion radial   : u_r/r_int  (Lame, colapso funcional del canal)

MODELO DE PRESION Co-SWIFT:
  P_poro = P_ext * 0.85 * 0.03  (~765 Pa @ 30kPa)
  Los sigma_fallo de hidrogeles (20-55 kPa) >> P_poro => ventana OPTIMA
"""

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle
import warnings
warnings.filterwarnings('ignore')

BG_DARK  = '#0A0E1A'
BG_PANEL = '#0F1623'
BG_CARD  = '#141B2D'
CYAN     = '#00D4FF'
GREEN    = '#00FF94'
GOLD     = '#FFB800'
RED      = '#FF3366'
VIOLET   = '#A855F7'
WHITE    = '#F0F4FF'
DIM      = '#6B7A99'
GRID_C   = '#1E2A40'

P_SWIFT     = 30_000.0
TAU_LIMITE  = 150.0
R_INT       = 200e-6
R_EXT       = 600e-6
ETA_BIOINK  = 0.05
Q_CAPILAR   = 0.5e-12
Q_ARTERIOLA = 2.0e-12
Q_VENULA    = 5.0e-12
POROSIDAD   = 0.85
F_TRANS     = 0.03
LIMITE_DR   = 40.0

MATERIALES = {
    'GelMA 7%':      {'E_min':2000,'E_max':8000,'E_mean':5000,'eta':120, 'sigma_fallo':35000,'nu':0.49,'color':CYAN,  'limite_dr':40.0},
    'Alginato 1.5%': {'E_min':1000,'E_max':3000,'E_mean':2000,'eta':80,  'sigma_fallo':20000,'nu':0.48,'color':GREEN, 'limite_dr':60.0},
    'NICE Bioink':   {'E_min':3000,'E_max':9000,'E_mean':6000,'eta':200, 'sigma_fallo':55000,'nu':0.49,'color':GOLD,  'limite_dr':40.0},
    'dECM Espinaca': {'E_min':1000,'E_max':5000,'E_mean':3000,'eta':60,  'sigma_fallo':28000,'nu':0.45,'color':VIOLET,'limite_dr':50.0},
}

def p_poro(P_ext):
    return P_ext * POROSIDAD * F_TRANS

def von_mises_lame(Pp, r_int=R_INT, r_ext=R_EXT, nu=0.49):
    A  = r_ext**2 - r_int**2
    st = -2.0 * Pp * r_ext**2 / A
    sr = 0.0
    sz = nu * (sr + st)
    vm = np.sqrt(0.5 * ((sr-st)**2 + (st-sz)**2 + (sz-sr)**2))
    return abs(st), abs(sz), vm

def deformacion_radial(Pp, E, nu, r_int=R_INT, r_ext=R_EXT):
    A   = r_ext**2 - r_int**2
    u_r = Pp * r_ext**2 * r_int * (1.0 + nu) / (E * A)
    return (u_r / r_int) * 100.0

def shear_wall(Q, eta=ETA_BIOINK, r=R_INT):
    return (4.0 * eta * Q) / (np.pi * r**3)

def kelvin_voigt(E, eta, sigma, t):
    eps_t  = (sigma / E) * (1.0 - np.exp(-t * E / eta))
    eps_eq = sigma / E
    return eps_t, eps_eq

def ventana_segura(props):
    E       = props['E_mean']
    nu      = props['nu']
    sf      = props['sigma_fallo']
    lim_dr  = props['limite_dr']   # limite especifico por material
    tau_ref = shear_wall(Q_CAPILAR)
    crit_C  = tau_ref < TAU_LIMITE
    Ps  = np.linspace(0, 60000, 800)
    res = []
    for P_ext in Ps:
        Pp       = p_poro(P_ext)
        _,_,vm   = von_mises_lame(Pp, nu=nu)
        dr       = deformacion_radial(Pp, E, nu)
        seguro   = (vm < sf) and (dr < lim_dr) and crit_C
        res.append({'P_ext':P_ext,'Pp':Pp,'vm':vm,'dr':dr,'tau':tau_ref,'seguro':seguro})
    ps = [r['P_ext'] for r in res if r['seguro']]
    return Ps, res, (min(ps) if ps else 0.0), (max(ps) if ps else 0.0)

def evaluar(props):
    Ps, res, p0, p1 = ventana_segura(props)
    vent     = p1 - p0
    swift_ok = any(r['seguro'] and abs(r['P_ext'] - P_SWIFT) < 2000 for r in res)
    if   swift_ok and vent >= 25000: estado,col = 'OPTIMO', GREEN
    elif swift_ok and vent >= 5000:  estado,col = 'PARCIAL',GOLD
    elif swift_ok:                   estado,col = 'PARCIAL',GOLD
    else:                            estado,col = 'CRITICO',RED
    return {'estado':estado,'color':col,'p0_kPa':p0/1000,'p1_kPa':p1/1000,'vent_Pa':vent,'swift_ok':swift_ok}

def ax_style(ax, titulo, xl, yl, xlim=None, ylim=None):
    ax.set_facecolor(BG_CARD)
    for sp in ax.spines.values(): sp.set_color(GRID_C)
    ax.tick_params(colors=DIM, labelsize=8)
    ax.xaxis.label.set_color(DIM); ax.yaxis.label.set_color(DIM)
    ax.set_title(titulo, color=WHITE, fontsize=9.5, fontweight='bold', pad=8)
    ax.set_xlabel(xl, fontsize=8); ax.set_ylabel(yl, fontsize=8)
    ax.grid(True, color=GRID_C, lw=0.5, alpha=0.7); ax.set_axisbelow(True)
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)

def construir_figura():
    EVS = {nm: evaluar(pr) for nm, pr in MATERIALES.items()}
    Pp_30 = p_poro(P_SWIFT)

    fig = plt.figure(figsize=(22, 26), facecolor=BG_DARK)
    gs  = gridspec.GridSpec(6, 4, figure=fig, hspace=0.52, wspace=0.38,
                            top=0.93, bottom=0.04, left=0.06, right=0.97)

    # HEADER
    ah = fig.add_subplot(gs[0,:])
    ah.set_facecolor(BG_PANEL); ah.axis('off')
    ah.axhline(0.95, color=CYAN, lw=1.5, alpha=0.6)
    ah.text(0.5,0.82,'BIO-KIDNEY AI 2026  |  VirtusSapiens', transform=ah.transAxes,
            ha='center', fontsize=10, color=DIM, fontfamily='monospace')
    ah.text(0.5,0.52,'SIMULADOR DE ESTRES MECANICO EN dECM', transform=ah.transAxes,
            ha='center', fontsize=20, color=WHITE, fontweight='bold')
    ah.text(0.5,0.20,'Modulo 09  |  Validacion estructural Co-SWIFT  |  Bioimpresion renal  |  v5',
            transform=ah.transAxes, ha='center', fontsize=9.5, color=CYAN, fontfamily='monospace')
    ah.axhline(0.05, color=CYAN, lw=0.5, alpha=0.3)

    # P1 Estres vs Deformacion
    a1 = fig.add_subplot(gs[1,:2])
    ax_style(a1,'Estres vs Deformacion - Elasticidad Lineal','Deformacion epsilon','Estres sigma (kPa)',xlim=(0,0.15),ylim=(0,12))
    eps_r = np.linspace(0,0.15,300)
    for nm,pr in MATERIALES.items():
        for Ev,ls,lb in [(pr['E_min'],'--',None),(pr['E_max'],'-',nm)]:
            a1.plot(eps_r, Ev*eps_r/1000, color=pr['color'], ls=ls,
                    lw=1.8 if ls=='-' else 0.9, alpha=0.92 if ls=='-' else 0.4, label=lb)
        a1.axhline(pr['sigma_fallo']/1000, color=pr['color'], lw=0.7, ls=':', alpha=0.5)
    a1.axhline(Pp_30/1000, color=RED, lw=1.8, ls='-.', label=f'P poro @ 30kPa = {Pp_30:.0f} Pa')
    a1.legend(loc='upper left', fontsize=7.5, facecolor=BG_CARD, edgecolor=GRID_C, labelcolor=WHITE, framealpha=0.95)

    # P2 Kelvin-Voigt
    a2 = fig.add_subplot(gs[1,2:])
    ax_style(a2,'Respuesta Viscoelastica - Kelvin-Voigt','Tiempo (s)','Deformacion epsilon',xlim=(0,60))
    t = np.linspace(0,60,500)
    for nm,pr in MATERIALES.items():
        eps_t,eps_eq = kelvin_voigt(pr['E_mean'],pr['eta'],Pp_30,t)
        a2.plot(t, eps_t, color=pr['color'], lw=2.0, label=nm, alpha=0.92)
        a2.axhline(eps_eq, color=pr['color'], lw=0.7, ls=':', alpha=0.4)
    a2.axhline(0.10, color=RED, lw=1.5, ls='--', label='Limite deformacion 10%', alpha=0.85)
    a2.legend(loc='lower right', fontsize=7.5, facecolor=BG_CARD, edgecolor=GRID_C, labelcolor=WHITE, framealpha=0.95)

    # P3 Hagen-Poiseuille
    a3 = fig.add_subplot(gs[2,:2])
    ax_style(a3,'Shear Stress en Canal - Hagen-Poiseuille','Radio canal (um)','Shear stress tau (Pa)',xlim=(150,600),ylim=(0,400))
    radios_um = np.linspace(150,600,300)
    radios_m  = radios_um*1e-6
    for lbl,Q,ls,col in [('Q=0.5 nL/s  capilar ref',Q_CAPILAR,'-',CYAN),
                          ('Q=2.0 nL/s  arteriola',Q_ARTERIOLA,'--',DIM),
                          ('Q=5.0 nL/s  venula',Q_VENULA,':',DIM)]:
        tau_v = [shear_wall(Q, r=r) for r in radios_m]
        a3.plot(radios_um, tau_v, color=col, lw=2.0 if Q==Q_CAPILAR else 1.3, ls=ls, label=lbl, alpha=0.9)
    a3.axhline(TAU_LIMITE, color=RED, lw=1.8, ls='-.', label=f'Limite citotoxico {TAU_LIMITE} Pa')
    a3.axvline(200, color=GOLD, lw=1.5, ls='--', label='Radio min Co-SWIFT 200 um', alpha=0.8)
    a3.fill_between(radios_um, 0, TAU_LIMITE,   color=GREEN, alpha=0.06)
    a3.fill_between(radios_um, TAU_LIMITE, 400, color=RED,   alpha=0.06)
    a3.legend(loc='upper right', fontsize=7.5, facecolor=BG_CARD, edgecolor=GRID_C, labelcolor=WHITE, framealpha=0.95)

    # P4 Von Mises vs P_ext
    a4 = fig.add_subplot(gs[2,2:])
    ax_style(a4,'Von Mises vs Presion Extrusion - Lame Pared Gruesa','Presion extrusion (kPa)','Von Mises (Pa)')
    Pext_arr = np.linspace(0,60000,400)
    Pp_arr   = np.array([p_poro(P) for P in Pext_arr])
    vm_all   = []
    for nm,pr in MATERIALES.items():
        vm_v = [von_mises_lame(Pp, nu=pr['nu'])[2] for Pp in Pp_arr]
        vm_all.extend(vm_v)
        a4.plot(Pext_arr/1000, vm_v, color=pr['color'], lw=2.0, label=nm, alpha=0.92)
        a4.axhline(pr['sigma_fallo'], color=pr['color'], lw=0.8, ls=':', alpha=0.55)
    a4.axvline(30, color=RED, lw=1.8, ls='--', label='P SWIFT 30 kPa', alpha=0.9)
    ylim4 = max(vm_all)*1.25
    a4.fill_betweenx([0,ylim4], 0, 30, color=GREEN, alpha=0.04)
    a4.set_ylim(0, ylim4)
    a4.text(0.98,0.05, f'P poro @ 30kPa = {Pp_30:.0f} Pa\n(Pext x {POROSIDAD} x {F_TRANS})',
            transform=a4.transAxes, ha='right', va='bottom', fontsize=7.5, color=DIM, fontfamily='monospace')
    a4.legend(loc='upper left', fontsize=7.5, facecolor=BG_CARD, edgecolor=GRID_C, labelcolor=WHITE, framealpha=0.95)

    # P5 Deformacion radial
    a5 = fig.add_subplot(gs[3,:2])
    ax_style(a5,'Deformacion Radial Canal - Colapso Funcional','Presion extrusion (kPa)','Delta r / r_int (%)',xlim=(0,60))
    dr_max_all = []
    for nm,pr in MATERIALES.items():
        dr_v = [deformacion_radial(p_poro(P), pr['E_mean'], pr['nu']) for P in Pext_arr]
        dr_max_all.append(max(dr_v))
        a5.plot(Pext_arr/1000, dr_v, color=pr['color'], lw=2.0, label=nm, alpha=0.92)
        # Limite especifico por material (linea punteada del color del material)
        a5.axhline(pr['limite_dr'], color=pr['color'], lw=1.0, ls=':', alpha=0.7)
    # Limite referencia global 40%
    a5.axhline(40.0, color=RED, lw=1.8, ls='--', label='Limite colapso GelMA/NICE 40%', alpha=0.9)
    a5.axvline(30, color=RED, lw=1.3, ls=':', alpha=0.6)
    a5.fill_between(Pext_arr/1000, 0, 40.0, color=GREEN, alpha=0.05)
    a5.set_ylim(0, max(dr_max_all)*1.4)
    a5.legend(loc='upper left', fontsize=7.5, facecolor=BG_CARD, edgecolor=GRID_C, labelcolor=WHITE, framealpha=0.95)

    # P6 Ventana segura
    a6 = fig.add_subplot(gs[3,2:])
    ax_style(a6,'Ventana Segura de Presion de Extrusion','Presion (kPa)','Material',xlim=(0,65))
    nms = list(MATERIALES.keys())
    for i,(nm,pr) in enumerate(MATERIALES.items()):
        ev = EVS[nm]; p0=ev['p0_kPa']; p1=ev['p1_kPa']
        a6.barh(i, 60, left=0, height=0.55, color=BG_CARD, alpha=0.5, zorder=1)
        a6.barh(i, max(p1-p0,0.5), left=p0, height=0.55, color=pr['color'], alpha=0.82, zorder=2)
        a6.text(min(p1+1,62), i, f'{p0:.0f}-{p1:.0f} kPa  [{ev["estado"]}]',
                va='center', fontsize=8.5, color=ev['color'], fontfamily='monospace', fontweight='bold')
    a6.axvline(30, color=RED, lw=2.0, ls='--', label='P SWIFT 30 kPa', zorder=5)
    a6.set_yticks(np.arange(len(nms))); a6.set_yticklabels(nms, fontsize=9, color=WHITE)
    a6.set_ylim(-0.5, len(nms)-0.5)
    a6.legend(loc='lower right', fontsize=7.5, facecolor=BG_CARD, edgecolor=GRID_C, labelcolor=WHITE, framealpha=0.95)

    # P7 Mapa Von Mises
    a7 = fig.add_subplot(gs[4,:2])
    ax_style(a7,'Mapa Concentracion Von Mises - Andamio dECM','Posicion X (mm)','Posicion Y (mm)',xlim=(0,10),ylim=(0,10))
    xg,yg = np.linspace(0,10,200),np.linspace(0,10,200)
    X,Y   = np.meshgrid(xg,yg)
    chs   = [(2,2),(5,5),(8,8),(3.5,7.5),(6.5,3.0)]
    S     = 15*(1-np.exp(-0.3*np.abs(X-5)))+12*(1-np.exp(-0.3*np.abs(Y-5)))
    for cx,cy in chs:
        d = np.sqrt((X-cx)**2+(Y-cy)**2)
        S += 35*np.exp(-d**2/0.25)+15*np.exp(-d**2/1.0)
    S = S/S.max()*von_mises_lame(Pp_30)[2]
    cmap_s = LinearSegmentedColormap.from_list('bk',[BG_DARK,'#0A2A4A',CYAN,GOLD,RED],N=256)
    im = a7.contourf(X,Y,S,levels=30,cmap=cmap_s,alpha=0.92)
    for cx,cy in chs:
        a7.add_patch(Circle((cx,cy),0.2,color=WHITE,fill=False,lw=1.5,alpha=0.8))
        a7.plot(cx,cy,'o',color=WHITE,ms=2,zorder=6)
    cb = plt.colorbar(im, ax=a7, fraction=0.035, pad=0.02)
    cb.set_label('Von Mises (Pa)', color=DIM, fontsize=8)
    cb.ax.yaxis.set_tick_params(labelcolor=DIM)
    a7.text(0.5,0.02,'Canales Co-SWIFT (circulos)  |  Von Mises @ P_ext=30kPa',
            transform=a7.transAxes, ha='center', fontsize=7.5, color=DIM, fontfamily='monospace')

    # P8 Comparativa
    a8 = fig.add_subplot(gs[4,2:])
    ax_style(a8,'Comparativa de Materiales - Indices de Rendimiento','Categoria','Puntuacion (0-10)',ylim=(0,11))
    cats = ['Modulo\nElastico','Resist.\nCompresion','Ventana\nSegura','Estab.\nViscoelast.']
    xc   = np.arange(len(cats)); nmat=len(MATERIALES); bw=0.18
    offs = np.linspace(-(nmat-1)/2,(nmat-1)/2,nmat)*bw
    for i,(nm,pr) in enumerate(MATERIALES.items()):
        ev = EVS[nm]
        sc = [min(pr['E_mean']/9000*10,10), min(pr['sigma_fallo']/55000*10,10),
              min(ev['vent_Pa']/60000*10,10), min(pr['eta']/200*10,10)]
        a8.bar(xc+offs[i], sc, bw, color=pr['color'], alpha=0.82, label=nm, zorder=3)
    a8.set_xticks(xc); a8.set_xticklabels(cats, fontsize=8, color=WHITE)
    a8.axhline(7, color=DIM, lw=0.7, ls=':', alpha=0.5)
    a8.legend(loc='upper right', fontsize=7.5, facecolor=BG_CARD, edgecolor=GRID_C, labelcolor=WHITE, framealpha=0.95)

    # ESTADO GLOBAL
    ae = fig.add_subplot(gs[5,:])
    ae.set_facecolor(BG_PANEL); ae.axis('off')
    ae.axhline(0.98, color=CYAN, lw=0.7, alpha=0.4)
    ae.text(0.5,0.88,'ESTADO GLOBAL - SIMULADOR ESTRES MECANICO dECM',
            transform=ae.transAxes, ha='center', fontsize=11, color=DIM, fontfamily='monospace')
    cnt = {'OPTIMO':0,'PARCIAL':0,'CRITICO':0}
    for ev in EVS.values(): cnt[ev['estado']] += 1
    if   cnt['CRITICO']==0 and cnt['OPTIMO']>=3: eg,cg='OPTIMO',GREEN; desc='Andamio dECM validado - Co-SWIFT seguro a 30 kPa en todos los materiales'
    elif cnt['CRITICO']==0:                       eg,cg='OPTIMO',GREEN; desc='Andamio dECM validado - Co-SWIFT seguro a 30 kPa'
    elif cnt['CRITICO']<=1:                       eg,cg='PARCIAL',GOLD; desc='Mayoria de materiales validados - revisar material con ventana reducida'
    else:                                          eg,cg='CRITICO',RED;  desc='Revision urgente de parametros de bioimpresion'
    ae.text(0.5,0.64,f'[ {eg} ]', transform=ae.transAxes, ha='center', fontsize=26,
            color=cg, fontweight='bold', fontfamily='monospace')
    ae.text(0.5,0.44,desc, transform=ae.transAxes, ha='center', fontsize=10.5, color=WHITE)
    for j,(nm,ev) in enumerate(EVS.items()):
        xj = 0.13+j*0.24
        ae.text(xj,0.28,nm, transform=ae.transAxes, ha='center', fontsize=8, color=DIM, fontfamily='monospace')
        ae.text(xj,0.16,f'[ {ev["estado"]} ]', transform=ae.transAxes, ha='center', fontsize=10,
                color=ev['color'], fontweight='bold', fontfamily='monospace')
        ae.text(xj,0.06,f'{ev["p0_kPa"]:.0f}-{ev["p1_kPa"]:.0f} kPa', transform=ae.transAxes,
                ha='center', fontsize=8.5, color=ev['color'], fontfamily='monospace')
    ae.axhline(0.02, color=CYAN, lw=0.5, alpha=0.3)
    fig.text(0.97,0.005,'Bio-Kidney AI 2026  |  Carlos David Moreno Caceres  |  VirtusSapiens',
             ha='right', va='bottom', fontsize=7, color=DIM, fontfamily='monospace', alpha=0.6)
    return fig, EVS

def reporte(EVS):
    sep   = '='*72
    Pp_30 = p_poro(P_SWIFT)
    tau_r = shear_wall(Q_CAPILAR)
    vm_30 = von_mises_lame(Pp_30)[2]
    print(f'\n{sep}')
    print('  BIO-KIDNEY AI 2026 | VirtusSapiens')
    print('  SIMULADOR ESTRES MECANICO dECM v5 - REPORTE')
    print(sep)
    print(f'  P_ext SWIFT          : {P_SWIFT/1000:.0f} kPa')
    print(f'  P_poro (Pext x {POROSIDAD} x {F_TRANS}): {Pp_30:.1f} Pa')
    print(f'  Von Mises @ P_poro   : {vm_30:.4f} Pa')
    print(f'  Shear (Q=0.5nL/s)   : {tau_r:.5f} Pa  (limite {TAU_LIMITE} Pa)')
    print(f'  Canal: r_int=200um, r_ext=600um  |  Limite colapso: por material (40-60%)')
    print(f'{"-"*72}')
    print(f'{"Material":<20} {"Estado":<10} {"Ventana kPa":<22} {"SWIFT":<6} {"Vent total":<14} {"Lim DR"}')
    print(f'{"-"*72}')
    for nm,ev in EVS.items():
        rng    = f'{ev["p0_kPa"]:.0f}-{ev["p1_kPa"]:.0f} kPa'
        sw     = 'SI' if ev['swift_ok'] else 'NO'
        lim_dr = MATERIALES[nm]['limite_dr']
        print(f'{nm:<20} {ev["estado"]:<10} {rng:<22} {sw:<6} {ev["vent_Pa"]/1000:.1f} kPa       {lim_dr:.0f}%')
    print(f'{"-"*72}')
    print('  OK Elasticidad lineal  : sigma = E*epsilon')
    print('  OK Kelvin-Voigt        : epsilon(t) = (sigma/E)*(1-exp(-Et/eta))')
    print('  OK Hagen-Poiseuille    : tau = 4*eta*Q/(pi*r^3)  Q=0.5nL/s fisiologico')
    print('  OK Von Mises (Lame)    : cilindro pared gruesa, presion externa')
    print('  OK Deformacion radial  : u_r = Pp*r_ext^2*r_int*(1+nu)/(E*A)')
    print(f'{sep}\n')

def main():
    print('\n  Bio-Kidney AI 2026 | Simulador Estres Mecanico dECM v5\n')
    for msg in ['  Calculando elasticidad lineal...',
                '  Modelando respuesta Kelvin-Voigt...',
                '  Evaluando Hagen-Poiseuille (caudal fisiologico)...',
                '  Aplicando Von Mises - Lame pared gruesa...',
                '  Calculando deformacion radial de canales...',
                '  Construyendo mapa de zonas criticas...',
                '  Generando ventana segura de presion...']:
        print(msg)
    fig, EVS = construir_figura()
    reporte(EVS)
    import os
    out = os.path.expanduser('~/Escritorio/BioKidney-AI/01_simuladores/')
    os.makedirs(out, exist_ok=True)
    png = os.path.join(out, 'simulador_estres_mecanico_dECM.png')
    fig.savefig(png, dpi=150, bbox_inches='tight', facecolor=BG_DARK, edgecolor='none')
    print(f'\n  Figura guardada: {png}')
    plt.tight_layout()
    plt.show()
    print('\n  Completado | VirtusSapiens Bio-Kidney AI 2026\n')

if __name__ == '__main__':
    main()
