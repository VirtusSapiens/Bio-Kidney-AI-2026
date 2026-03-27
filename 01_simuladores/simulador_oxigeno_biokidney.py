#!/usr/bin/env python3
import os, sys, time, warnings
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path
warnings.filterwarnings("ignore")

class ParametrosBiofisicos:
    D_O2             = 2.0e-5
    M_max            = 5.0e-4
    P_50             = 1.0
    P_arteriolar     = 40.0
    P_venular        = 20.0
    radio_influencia = 150.0
    P_hipoxia        = 1.0

class ConfiguracionGrid:
    nx, ny, nz = 60, 60, 40
    max_iter   = 5000
    tolerancia = 1e-4
    relajacion = 1.6
    def __init__(self, xmin, xmax, ymin, ymax, zmin, zmax):
        m = 0.05
        self.Lx0=xmin-m; self.Lx1=xmax+m
        self.Ly0=ymin-m; self.Ly1=ymax+m
        self.Lz0=zmin-m; self.Lz1=zmax+m
        self.Lx=self.Lx1-self.Lx0
        self.Ly=self.Ly1-self.Ly0
        self.Lz=self.Lz1-self.Lz0
    @property
    def dx(self): return self.Lx/(self.nx-1)
    @property
    def dy(self): return self.Ly/(self.ny-1)
    @property
    def dz(self): return self.Lz/(self.nz-1)
    def xs(self): return np.linspace(self.Lx0,self.Lx1,self.nx)
    def ys(self): return np.linspace(self.Ly0,self.Ly1,self.ny)
    def zs(self): return np.linspace(self.Lz0,self.Lz1,self.nz)

def cargar_arbol_cco(ruta_csv, params):
    ruta = Path(ruta_csv)
    if not ruta.exists():
        print(f"\n⚠️  CSV no encontrado. Usando árbol sintético...")
        return _generar_sintetico(params)
    print(f"\n📂 Cargando árbol CCO v7: {ruta}")
    df = pd.read_csv(ruta)
    print(f"   ✓ {len(df)} segmentos | columnas: {list(df.columns)}")
    x0=df['x1_mm'].values*0.1; y0=df['y1_mm'].values*0.1; z0=df['z1_mm'].values*0.1
    x1=df['x2_mm'].values*0.1; y1=df['y2_mm'].values*0.1; z1=df['z2_mm'].values*0.1
    radios_um=df['radio_um'].values
    sistemas=df['sistema'].astype(str).values
    def pres(s):
        s=s.lower().strip()
        if s.startswith('art'): return params.P_arteriolar
        elif s.startswith('ven'): return params.P_venular
        return (params.P_arteriolar+params.P_venular)/2.0
    segs=[{'p0':np.array([x0[i],y0[i],z0[i]],dtype=float),
           'p1':np.array([x1[i],y1[i],z1[i]],dtype=float),
           'radio_um':float(radios_um[i]),
           'presion':pres(sistemas[i]),
           'sistema':sistemas[i]} for i in range(len(df))]
    tx=np.concatenate([x0,x1]); ty=np.concatenate([y0,y1]); tz=np.concatenate([z0,z1])
    cfg=ConfiguracionGrid(tx.min(),tx.max(),ty.min(),ty.max(),tz.min(),tz.max())
    print(f"   ✓ Dominio: {cfg.Lx*10:.1f}x{cfg.Ly*10:.1f}x{cfg.Lz*10:.1f} mm")
    return segs, cfg

def _generar_sintetico(params):
    np.random.seed(42)
    Lx0,Lx1=-6.,6.; Ly0,Ly1=-3.,3.; Lz0,Lz1=-2.,2.
    defs=[('art',params.P_arteriolar,483),('ven',params.P_venular,482),
          ('col',(params.P_arteriolar+params.P_venular)/2,483)]
    segs=[]
    for sist,pres,n in defs:
        for _ in range(n):
            cx=np.random.uniform(Lx0*.9,Lx1*.9); cy=np.random.uniform(Ly0*.9,Ly1*.9)
            cz=np.random.uniform(Lz0*.9,Lz1*.9); a=np.random.uniform(0,2*np.pi)
            b=np.random.uniform(-.3,.3); L=np.random.exponential(.4)
            ex=np.clip(cx+L*np.cos(a)*np.cos(b),Lx0,Lx1)
            ey=np.clip(cy+L*np.sin(a)*np.cos(b),Ly0,Ly1)
            ez=np.clip(cz+L*np.sin(b),Lz0,Lz1)
            r=np.clip(np.random.lognormal(np.log(50),.5),5,500)
            segs.append({'p0':np.array([cx,cy,cz]),'p1':np.array([ex,ey,ez]),
                         'radio_um':r,'presion':pres,'sistema':sist})
    cfg=ConfiguracionGrid(Lx0,Lx1,Ly0,Ly1,Lz0,Lz1)
    print(f"   ✓ Árbol sintético: {len(segs)} segmentos")
    return segs,cfg

def mapear_vasculatura(segmentos, config, params):
    cfg=config; r_inf=params.radio_influencia*1e-4
    xs=cfg.xs(); ys=cfg.ys(); zs=cfg.zs()
    mascara=np.zeros((cfg.nx,cfg.ny,cfg.nz),dtype=bool)
    presion=np.zeros((cfg.nx,cfg.ny,cfg.nz))
    peso_acc=np.zeros((cfg.nx,cfg.ny,cfg.nz))
    print(f"\n🗺️  Mapeando {len(segmentos)} segmentos → grid {cfg.nx}x{cfg.ny}x{cfg.nz}...")
    t0=time.time()
    for seg in segmentos:
        p0=seg['p0']; p1=seg['p1']; vp=p1-p0
        long=np.linalg.norm(vp)
        if long<1e-12: continue
        r_seg=max(r_inf,seg['radio_um']*1e-4)
        xlo=max(xs[0],min(p0[0],p1[0])-r_seg); xhi=min(xs[-1],max(p0[0],p1[0])+r_seg)
        ylo=max(ys[0],min(p0[1],p1[1])-r_seg); yhi=min(ys[-1],max(p0[1],p1[1])+r_seg)
        zlo=max(zs[0],min(p0[2],p1[2])-r_seg); zhi=min(zs[-1],max(p0[2],p1[2])+r_seg)
        ix0=max(0,int(np.searchsorted(xs,xlo))-1); ix1=min(cfg.nx,int(np.searchsorted(xs,xhi))+2)
        iy0=max(0,int(np.searchsorted(ys,ylo))-1); iy1=min(cfg.ny,int(np.searchsorted(ys,yhi))+2)
        iz0=max(0,int(np.searchsorted(zs,zlo))-1); iz1=min(cfg.nz,int(np.searchsorted(zs,zhi))+2)
        if ix0>=ix1 or iy0>=iy1 or iz0>=iz1: continue
        # FIX: meshgrid indexing='ij' → mismo shape (N,M,K) para los tres arrays
        gx,gy,gz=np.meshgrid(xs[ix0:ix1],ys[iy0:iy1],zs[iz0:iz1],indexing='ij')
        wx=gx-p0[0]; wy=gy-p0[1]; wz=gz-p0[2]
        t_=np.clip((wx*vp[0]+wy*vp[1]+wz*vp[2])/(long*long+1e-24),0.,1.)
        cx_=p0[0]+t_*vp[0]; cy_=p0[1]+t_*vp[1]; cz_=p0[2]+t_*vp[2]
        dist=np.sqrt((gx-cx_)**2+(gy-cy_)**2+(gz-cz_)**2)
        en=dist<=r_seg
        sigma=r_seg*.4
        peso=np.where(en,np.exp(-.5*(dist/sigma)**2),0.)
        mascara[ix0:ix1,iy0:iy1,iz0:iz1]|=en
        presion[ix0:ix1,iy0:iy1,iz0:iz1]+=peso*seg['presion']
        peso_acc[ix0:ix1,iy0:iy1,iz0:iz1]+=peso
    ok=peso_acc>1e-30
    presion[ok]=presion[ok]/peso_acc[ok]
    presion[~mascara]=0.
    n_v=int(mascara.sum())
    print(f"   ✓ Voxeles vasculares: {n_v:,}/{mascara.size:,} ({100.*n_v/mascara.size:.1f}%)  {time.time()-t0:.1f}s")
    return mascara,presion

def resolver_fick_3d(mascara, presion_vasc, config, params):
    cfg=config; dx2=cfg.dx**2; dy2=cfg.dy**2; dz2=cfg.dz**2
    omega=cfg.relajacion; D=params.D_O2; Mm=params.M_max; P50=params.P_50
    coef=2./dx2+2./dy2+2./dz2
    P=np.full((cfg.nx,cfg.ny,cfg.nz),20.)
    P[mascara]=presion_vasc[mascara]
    tejido=~mascara
    print(f"\n🔬 Solver Fick 3D — SOR (ω={omega})  grid {cfg.nx}x{cfg.ny}x{cfg.nz}")
    t0=time.time()
    for it in range(cfg.max_iter):
        P_old=P.copy()
        i=slice(1,cfg.nx-1); j=slice(1,cfg.ny-1); k=slice(1,cfg.nz-1)
        Pi=P[i,j,k]
        lap=((P[2:,1:-1,1:-1]+P[:-2,1:-1,1:-1])/dx2+
             (P[1:-1,2:,1:-1]+P[1:-1,:-2,1:-1])/dy2+
             (P[1:-1,1:-1,2:]+P[1:-1,1:-1,:-2])/dz2)
        Pp=np.maximum(Pi,0.)
        Mv=Mm*Pp/(Pp+P50+1e-12)
        Pgs=(lap-Mv/D)/coef
        P[i,j,k]=np.where(tejido[i,j,k],(1.-omega)*Pi+omega*Pgs,Pi)
        P[0,:,:]=P[1,:,:]; P[-1,:,:]=P[-2,:,:]
        P[:,0,:]=P[:,1,:]; P[:,-1,:]=P[:,-2,:]
        P[:,:,0]=P[:,:,1]; P[:,:,-1]=P[:,:,-2]
        P[mascara]=presion_vasc[mascara]
        P=np.maximum(P,0.)
        res=float(np.max(np.abs(P-P_old)))
        if it%500==0: print(f"   iter {it:5d}  residuo={res:.5f} mmHg")
        if res<cfg.tolerancia:
            print(f"\n   ✅ Convergido en {it+1} iter ({time.time()-t0:.1f}s) residuo={res:.2e}")
            break
    else:
        print(f"\n   ⚠️  Max iter. Residuo={res:.4f} mmHg")
    return P

def analizar_oxigenacion(P, mascara, params):
    tejido=~mascara; Pt=P[tejido].astype(float); u=params.P_hipoxia
    n=int(tejido.sum()); nh=int((Pt<u).sum()); no=int((Pt>=u).sum())
    return {'n_total':n,'n_hipoxicos':nh,'n_oxigenados':no,
            'pct_hipoxia':100.*nh/max(n,1),'pct_oxigenado':100.*no/max(n,1),
            'P_min':float(np.nanmin(Pt)),'P_max':float(np.nanmax(Pt)),
            'P_media':float(np.nanmean(Pt)),
            'P_p5':float(np.nanpercentile(Pt,5)),'P_p95':float(np.nanpercentile(Pt,95)),
            'coords_hipo':np.argwhere((P<u)&tejido),'objetivo':nh==0}

def imprimir_reporte(r, params):
    print("\n"+"═"*62)
    print("  📊  REPORTE DE OXIGENACIÓN — Bio-Kidney AI 2026")
    print("═"*62)
    print(f"  Voxeles tejido:   {r['n_total']:>10,}")
    print(f"  Oxigenados:       {r['n_oxigenados']:>10,}  ({r['pct_oxigenado']:6.2f}%)")
    print(f"  Hipóxicos:        {r['n_hipoxicos']:>10,}  ({r['pct_hipoxia']:6.2f}%)")
    print("─"*62)
    print(f"  PO₂ mínima:    {r['P_min']:>8.3f} mmHg")
    print(f"  PO₂ media:     {r['P_media']:>8.3f} mmHg")
    print(f"  PO₂ máxima:    {r['P_max']:>8.3f} mmHg")
    print(f"  PO₂ perc.  5%: {r['P_p5']:>8.3f} mmHg")
    print(f"  PO₂ perc. 95%: {r['P_p95']:>8.3f} mmHg")
    print("─"*62)
    if r['objetivo']: print("  🎯  OBJETIVO ALCANZADO: 100% del tejido oxigenado")
    else: print(f"  ⚠️  {r['n_hipoxicos']} voxeles hipóxicos → señal enviada a CCO")
    print("═"*62)

def exportar_senyal_cco(r, config):
    coords=r['coords_hipo']
    if len(coords)==0:
        print("\n  ✅ Sin zonas hipóxicas"); return pd.DataFrame()
    xs=config.xs(); ys=config.ys(); zs=config.zs()
    filas=[{'x_mm':round(xs[ix]*10,4),'y_mm':round(ys[iy]*10,4),
             'z_mm':round(zs[iz]*10,4),'prioridad':'CRITICA'} for ix,iy,iz in coords]
    df=pd.DataFrame(filas)
    ruta=Path("~/Escritorio/BioKidney-AI/02_vascular_cco/senyal_hipoxia_para_cco.csv").expanduser()
    ruta.parent.mkdir(parents=True,exist_ok=True)
    df.to_csv(ruta,index=False)
    print(f"\n  💾 Señal CCO: {ruta}  ({len(df)} zonas hipóxicas)")
    return df

def visualizar(P, mascara, r, config, params):
    matplotlib.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
        'axes.titlesize':11,'figure.facecolor':'#0f1117','axes.facecolor':'#1a1d2e',
        'text.color':'#e0e0e0','axes.labelcolor':'#e0e0e0',
        'xtick.color':'#aaaaaa','ytick.color':'#aaaaaa'})
    cfg=config; u=params.P_hipoxia; vmax=params.P_arteriolar
    cmap=LinearSegmentedColormap.from_list('O2',['#0d1b7a','#1a4baf','#2196f3',
        '#26c6da','#66bb6a','#ffee58','#ff7043','#b71c1c'])
    fig=plt.figure(figsize=(18,12),facecolor='#0f1117')
    fig.suptitle('BIO-KIDNEY AI 2026 · Simulador de Difusión de Oxígeno\n'
        'VirtusSapiens · Carlos David Moreno Cáceres · Medellín',
        fontsize=14,color='#64b5f6',fontweight='bold',y=.98)
    gs=gridspec.GridSpec(2,3,figure=fig,hspace=.38,wspace=.32)
    xm=cfg.xs()*10; ym=cfg.ys()*10; zm=cfg.zs()*10

    ax1=fig.add_subplot(gs[0,0]); iz=cfg.nz//2; Pxy=P[:,:,iz].T
    im1=ax1.imshow(Pxy,origin='lower',cmap=cmap,vmin=0,vmax=vmax,
        extent=[xm[0],xm[-1],ym[0],ym[-1]],aspect='auto')
    ax1.contour(xm,ym,Pxy,levels=[u],colors='#ff4444',linewidths=1.5,linestyles='--')
    ax1.set_title(f'PO₂ Plano XY (z={zm[iz]:.1f}mm)'); ax1.set_xlabel('X(mm)'); ax1.set_ylabel('Y(mm)')
    plt.colorbar(im1,ax=ax1,label='PO₂(mmHg)')

    ax2=fig.add_subplot(gs[0,1]); iy=cfg.ny//2; Pxz=P[:,iy,:].T
    im2=ax2.imshow(Pxz,origin='lower',cmap=cmap,vmin=0,vmax=vmax,
        extent=[xm[0],xm[-1],zm[0],zm[-1]],aspect='auto')
    ax2.contour(xm,zm,Pxz,levels=[u],colors='#ff4444',linewidths=1.5,linestyles='--')
    ax2.set_title(f'PO₂ Plano XZ (y={ym[iy]:.1f}mm)'); ax2.set_xlabel('X(mm)'); ax2.set_ylabel('Z(mm)')
    plt.colorbar(im2,ax=ax2,label='PO₂(mmHg)')

    ax3=fig.add_subplot(gs[0,2]); tej=~mascara; Pt=P[tej].flatten()
    bins=np.linspace(0,vmax,60)
    ax3.hist(Pt[Pt<u],bins=bins,color='#ef5350',alpha=.85,label=f'Hipóxico({r["n_hipoxicos"]:,})')
    ax3.hist(Pt[Pt>=u],bins=bins,color='#42a5f5',alpha=.75,label=f'Oxigenado({r["n_oxigenados"]:,})')
    ax3.axvline(u,color='#ff4444',lw=2,ls='--',label=f'Umbral {u}mmHg')
    ax3.set_xlabel('PO₂(mmHg)'); ax3.set_ylabel('Voxeles'); ax3.set_title('Distribución PO₂')
    ax3.legend(fontsize=8)

    ax4=fig.add_subplot(gs[1,0],projection='3d'); ax4.set_facecolor('#1a1d2e')
    idx=np.argwhere(tej); ns=min(3000,len(idx)); p=np.random.permutation(len(idx))[:ns]; m=idx[p]
    sc=ax4.scatter(cfg.xs()[m[:,0]]*10,cfg.ys()[m[:,1]]*10,cfg.zs()[m[:,2]]*10,
        c=P[m[:,0],m[:,1],m[:,2]],cmap=cmap,vmin=0,vmax=vmax,s=2,alpha=.4)
    ax4.set_xlabel('X(mm)'); ax4.set_ylabel('Y(mm)'); ax4.set_zlabel('Z(mm)'); ax4.set_title('PO₂ 3D')
    plt.colorbar(sc,ax=ax4,label='PO₂(mmHg)',shrink=.6)

    ax5=fig.add_subplot(gs[1,1])
    rk=np.linspace(0,300e-4,300); Rc=params.radio_influencia*1e-4
    Pan=np.maximum(params.P_arteriolar-(params.M_max/(4*params.D_O2))*(rk**2-Rc**2),0)
    ax5.plot(rk*1e4,Pan,color='#42a5f5',lw=2.5,label='Krogh analítico')
    ax5.axhline(u,color='#ff4444',lw=1.5,ls='--',label=f'Umbral {u}mmHg')
    ax5.fill_between(rk*1e4,0,u,color='#ef5350',alpha=.15,label='Zona hipóxica')
    ax5.set_xlabel('Distancia capilar(µm)'); ax5.set_ylabel('PO₂(mmHg)')
    ax5.set_title('Perfil Krogh'); ax5.legend(fontsize=8)

    ax6=fig.add_subplot(gs[1,2]); ax6.axis('off')
    est="✅ OBJETIVO ALCANZADO" if r['objetivo'] else "⚠️  ZONAS HIPÓXICAS DETECTADAS"
    col='#66bb6a' if r['objetivo'] else '#ef5350'
    txt=(f"  SIMULACIÓN PO₂ — CCO v7\n  {'─'*28}\n"
         f"  Árbol vascular: 1448 segmentos\n  Grid: {cfg.nx}x{cfg.ny}x{cfg.nz}\n"
         f"  Dominio: {cfg.Lx*10:.0f}x{cfg.Ly*10:.0f}x{cfg.Lz*10:.0f} mm\n\n"
         f"  PO₂ mínima:    {r['P_min']:>7.3f} mmHg\n"
         f"  PO₂ media:     {r['P_media']:>7.3f} mmHg\n"
         f"  PO₂ máxima:    {r['P_max']:>7.3f} mmHg\n"
         f"  PO₂ perc.  5%: {r['P_p5']:>7.3f} mmHg\n"
         f"  PO₂ perc. 95%: {r['P_p95']:>7.3f} mmHg\n\n"
         f"  Oxigenado: {r['pct_oxigenado']:>6.2f}%\n"
         f"  Hipóxico:  {r['pct_hipoxia']:>6.2f}%\n\n"
         f"  Densidad: 200M cel/mL\n  Umbral: {params.P_hipoxia} mmHg\n")
    ax6.text(.05,.95,txt,transform=ax6.transAxes,fontsize=9,va='top',
        fontfamily='monospace',color='#b0bec5',
        bbox=dict(boxstyle='round',facecolor='#1a1d2e',alpha=.8))
    ax6.text(.5,.04,est,transform=ax6.transAxes,fontsize=11,ha='center',
        va='bottom',color=col,fontweight='bold')
    ruta=Path("~/Escritorio/BioKidney-AI/01_simuladores/resultado_oxigeno_biokidney.png").expanduser()
    ruta.parent.mkdir(parents=True,exist_ok=True)
    plt.savefig(ruta,dpi=150,bbox_inches='tight',facecolor='#0f1117',edgecolor='none')
    print(f"\n  🖼️  Figura: {ruta}")
    plt.show()

def main():
    print("\n╔"+"═"*64+"╗")
    print("║  BIO-KIDNEY AI 2026 · Simulador de Difusión de Oxígeno      ║")
    print("║  VirtusSapiens · Carlos David Moreno Cáceres · Medellín     ║")
    print("╚"+"═"*64+"╝\n")
    params=ParametrosBiofisicos()
    print(f"⚙️  D_O2={params.D_O2:.2e}  M_max={params.M_max:.2e}  P_50={params.P_50}mmHg  umbral={params.P_hipoxia}mmHg  r_inf={params.radio_influencia:.0f}µm")
    ruta=Path("~/Escritorio/BioKidney-AI/02_vascular_cco/arbol_vascular_cco_v7.csv").expanduser()
    segs,cfg=cargar_arbol_cco(str(ruta),params)
    mascara,pvasc=mapear_vasculatura(segs,cfg,params)
    P=resolver_fick_3d(mascara,pvasc,cfg,params)
    r=analizar_oxigenacion(P,mascara,params)
    imprimir_reporte(r,params)
    exportar_senyal_cco(r,cfg)
    print("\n📈 Generando visualizaciones...")
    visualizar(P,mascara,r,cfg,params)
    print("\n✅ Simulación completada.\n")

if __name__=="__main__":
    main()
