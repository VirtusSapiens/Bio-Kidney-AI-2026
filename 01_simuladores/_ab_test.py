#!/usr/bin/env python3
# _ab_test.py  —  A/B de contaminacion por sistema=col en el solver de O2 (arbol v8).
# NO modifica el original: importa simulador_oxigeno_biokidney.py y reutiliza
# rasterizador + consumo. Corre en dos regimenes:
#   PROD  : solver de produccion tal cual (Jacobi + omega=1.6)  -> DIVERGE a NaN.
#   ESTAB : mismo esquema pero omega=1.0 (Jacobi estable) para obtener el campo real.
# En cada regimen: A = con col, B = sin col. Misma grilla (arbol completo) en todo.
import matplotlib; matplotlib.use('Agg')
import sys, importlib.util
from pathlib import Path
import numpy as np

ROOT = Path.home() / 'Escritorio/BioKidney-AI'
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / '01_simuladores'))
spec = importlib.util.spec_from_file_location('o2mod', str(ROOT / '01_simuladores/simulador_oxigeno_biokidney.py'))
o2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(o2)
from biokidney.experts.cellular import CellularExpert as CE
cfgp = o2.cfg_physio
U = cfgp.P_HIPOXIA
print(f"P_HIPOXIA (umbral) = {U} mmHg")

V8 = ROOT / '02_vascular_cco/arbol_vascular_cco_v8.csv'
segs, cfg = o2.cargar_arbol_cco(str(V8))

def solver(mask, pv, omega, max_iter=8000, tol=1e-4):
    dx2,dy2,dz2 = cfg.dx**2, cfg.dy**2, cfg.dz**2; D = cfgp.D_O2
    coef = 2/dx2 + 2/dy2 + 2/dz2
    P = np.full((cfg.nx,cfg.ny,cfg.nz), 20.); P[mask] = pv[mask]; tejido = ~mask
    for it in range(max_iter):
        Po = P.copy(); i=slice(1,cfg.nx-1);j=slice(1,cfg.ny-1);k=slice(1,cfg.nz-1); Pi=P[i,j,k]
        lap=((P[2:,1:-1,1:-1]+P[:-2,1:-1,1:-1])/dx2+(P[1:-1,2:,1:-1]+P[1:-1,:-2,1:-1])/dy2+(P[1:-1,1:-1,2:]+P[1:-1,1:-1,:-2])/dz2)
        Pgs=(lap - CE.oxygen_consumption_rate(Pi)/D)/coef
        P[i,j,k]=np.where(tejido[i,j,k],(1-omega)*Pi+omega*Pgs,Pi)
        P[0,:,:]=P[1,:,:];P[-1,:,:]=P[-2,:,:];P[:,0,:]=P[:,1,:];P[:,-1,:]=P[:,-2,:];P[:,:,0]=P[:,:,1];P[:,:,-1]=P[:,:,-2]
        P[mask]=pv[mask]; P=np.maximum(P,0.)
        if np.isnan(P).any(): return P, it, 'DIVERGE(NaN)'
        r=float(np.max(np.abs(P-Po)))
        if r<tol: return P, it, f'conv r={r:.2e}'
    return P, max_iter, 'sin converger'

def stats(P, mask):
    t=~mask; Pt=P[t]; fin=np.isfinite(Pt); Ptf=Pt[fin]; u=U
    n=int(t.sum()); nfin=int(fin.sum())
    nh=int((Ptf<u).sum())
    return dict(n=n, nfin=nfin, pmed=float(Ptf.mean()) if nfin else float('nan'),
               pmin=float(Ptf.min()) if nfin else float('nan'),
               hip=100*nh/max(nfin,1), nh=nh)

def corre(regimen, omega):
    print(f"\n=========== REGIMEN {regimen}  (omega={omega}) ===========")
    ma,pva = o2.mapear_vasculatura(segs, cfg)
    segs_b=[s for s in segs if not str(s['sistema']).lower().strip().startswith('col')]
    mb,pvb = o2.mapear_vasculatura(segs_b, cfg)
    Pa,ita,sa = solver(ma,pva,omega); Pb,itb,sb = solver(mb,pvb,omega)
    ra,rb = stats(Pa,ma), stats(Pb,mb)
    print(f"A(con col): it={ita} {sa}")
    print(f"B(sin col): it={itb} {sb}")
    for lab,r in (('A',ra),('B',rb)):
        print(f"  {lab}: tejido={r['n']:,} finitos={r['nfin']:,}  "
              f"pO2_med={r['pmed']:.3f}  pO2_min={r['pmin']:.3f}  hipoxia={r['hip']:.3f}% ({r['nh']} vox)")
    fin=np.isfinite(Pa)&np.isfinite(Pb); diff=np.abs(np.where(fin,Pa-Pb,0.0))
    comun=(~ma)&(~mb)&fin
    print(f"  |A-B| grid-finito: media={diff[fin].mean():.4f} max={diff.max():.4f} mmHg "
          f"(>1mmHg:{(diff>1).sum()}, >5mmHg:{(diff>5).sum()})")
    dc=np.abs(Pa-Pb)[comun]
    print(f"  |A-B| tejido comun finito (n={comun.sum():,}): media={dc.mean():.4f} max={dc.max():.4f} mmHg")

corre("PROD (produccion tal cual)", 1.6)
corre("ESTAB (omega=1.0, diagnostico)", 1.0)
