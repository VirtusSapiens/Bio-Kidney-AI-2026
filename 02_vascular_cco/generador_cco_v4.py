import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv, os

# ── PARÁMETROS ─────────────────────────────────────────────────────────────────
ALPHA          = 3.0
SEMILLA        = 42
PUNTOS_SPLINE  = 10
N_NIVELES      = 6      # niveles de ramificación por sistema
RAMAS_POR_NODO = 2      # bifurcación binaria estricta

# Dimensiones riñón (m)
EJE_X, EJE_Y, EJE_Z = 0.055, 0.030, 0.025

# Radios raíz (m)
R_ARTERIA = 500e-6
R_VENA    = 600e-6

# Entradas en el hilio
HILIO_ART = np.array([ 0.0, -EJE_Y*0.85,  0.003])
HILIO_VEN = np.array([ 0.0, -EJE_Y*0.85, -0.003])
HILIO_URI = np.array([ 0.0, -EJE_Y*0.85,  0.000])

np.random.seed(SEMILLA)

# ── DOMINIO ────────────────────────────────────────────────────────────────────
def dentro(p, m=0.93):
    return (p[0]/EJE_X)**2+(p[1]/EJE_Y)**2+(p[2]/EJE_Z)**2 <= m**2

def proyectar(p, m=0.90):
    d = np.sqrt((p[0]/EJE_X)**2+(p[1]/EJE_Y)**2+(p[2]/EJE_Z)**2)
    if d > m:
        # Escalar cada componente respetando la forma elipsoide
        return np.array([p[0]/d*m*EJE_X/EJE_X,
                         p[1]/d*m*EJE_Y/EJE_Y,
                         p[2]/d*m*EJE_Z/EJE_Z])
    return p.copy()

# ── MURRAY ESTRICTO ────────────────────────────────────────────────────────────
def hijos_murray(r_padre, asim=0.0):
    """
    Calcula radios hijos EXCLUSIVAMENTE desde Murray.
    r³_padre = r³_h1 + r³_h2
    asim ∈ [-0.2, 0.2] para variación natural
    """
    f  = np.clip(0.5 + asim, 0.3, 0.7)
    r1 = r_padre * f**(1.0/ALPHA)
    r2 = r_padre * (1.0-f)**(1.0/ALPHA)
    # Verificación explícita
    assert abs(r1**ALPHA + r2**ALPHA - r_padre**ALPHA) / r_padre**ALPHA < 0.001
    return r1, r2

def verificar_murray(r_p, r1, r2):
    lhs = r_p**ALPHA
    rhs = r1**ALPHA + r2**ALPHA
    return abs(lhs - rhs) / lhs < 0.01

# ── NODO ───────────────────────────────────────────────────────────────────────
class Nodo:
    _id = 0
    def __init__(self, pos, radio, nivel=0, padre=None, sistema='art'):
        self.id      = Nodo._id; Nodo._id += 1
        self.pos     = np.array(pos, dtype=float)
        self.radio   = float(radio)
        self.nivel   = nivel
        self.padre   = padre
        self.hijos   = []
        self.sistema = sistema

# ── SPLINE ─────────────────────────────────────────────────────────────────────
def spline(p0, p1, curv=0.18, n=10):
    eje = p1-p0; L = np.linalg.norm(eje)
    if L < 1e-10: return np.array([p0,p1])
    en = eje/L
    perp = np.array([-en[1],en[0],0.0])
    if np.linalg.norm(perp)<1e-10: perp=np.array([0.,-en[2],en[1]])
    perp /= np.linalg.norm(perp)
    perp2 = np.cross(en, perp)
    pc = (p0+p1)/2 + perp*L*curv*np.random.uniform(-1,1) \
                   + perp2*L*curv*0.5*np.random.uniform(-1,1)
    if not dentro(pc): pc=(p0+p1)/2
    t = np.linspace(0,1,n)
    return np.outer((1-t)**2,p0)+np.outer(2*(1-t)*t,pc)+np.outer(t**2,p1)

# ── DIRECCIÓN ANATÓMICA ────────────────────────────────────────────────────────
def dir_anatomica(nivel, sistema):
    """
    Cada nivel tiene una tendencia de crecimiento anatómica:
    Nivel 0-1: desde hilio hacia interior (dirección +X)
    Nivel 2:   ramificación interlobar (arcos en Y/Z)
    Nivel 3:   arcuatas (horizontal, paralelas a cápsula)
    Nivel 4-5: interlobulillares (radial hacia cápsula)
    """
    if nivel <= 1:
        base = np.array([0.7, np.random.uniform(0.2,0.6),
                         np.random.uniform(-0.3,0.3)])
    elif nivel == 2:
        base = np.array([np.random.uniform(0.2,0.5),
                         np.random.uniform(-0.7,0.7),
                         np.random.uniform(-0.5,0.5)])
    elif nivel == 3:
        # Arcuatas — crecen en arco horizontal
        angulo = np.random.uniform(0, 2*np.pi)
        base   = np.array([np.cos(angulo)*0.3,
                           np.sin(angulo)*0.7,
                           np.random.uniform(-0.2,0.2)])
    else:
        # Interlobulillares — radial hacia la cápsula
        base = np.array([np.random.uniform(-0.3,0.3),
                         np.random.uniform(-0.3,0.3),
                         np.random.choice([-1,1]) *
                         np.random.uniform(0.6,1.0)])

    if sistema == 'ven':
        base = base * np.array([1, -1, 1])

    n = np.linalg.norm(base)
    return base/n if n > 1e-10 else np.array([1.,0.,0.])

def longitud_nivel(nivel, radio):
    """Longitud proporcional al radio y al nivel — escala alométrica"""
    factores = [1.0, 0.65, 0.45, 0.30, 0.18, 0.10]
    f = factores[min(nivel, len(factores)-1)]
    return radio * np.random.uniform(80, 160) * f * 1000

# ── GENERADOR JERÁRQUICO CON MURRAY ESTRICTO ───────────────────────────────────
def generar_sistema_murray(origen, r_raiz, sistema, n_niveles):
    """
    Árbol vascular donde CADA bifurcación satisface Murray exactamente.
    Los radios se propagan desde la raíz hacia los terminales.
    """
    Nodo._id = 0
    raiz  = Nodo(origen, r_raiz, 0, sistema=sistema)
    todos = [raiz]
    nivel_actual = [raiz]

    print(f"\n  Sistema {sistema.upper()} — propagación Murray:")
    print(f"    Nivel 0: radio raíz = {r_raiz*1e6:.0f} µm")

    for niv in range(1, n_niveles + 1):
        siguiente = []

        for padre in nivel_actual:
            # Generar exactamente 2 hijos por Murray
            asim = np.random.uniform(-0.18, 0.18)
            try:
                r1, r2 = hijos_murray(padre.radio, asim)
            except AssertionError:
                # Si falla la verificación, recalcular sin asimetría
                r1, r2 = hijos_murray(padre.radio, 0.0)

            for r_hijo, signo in [(r1, 1), (r2, -1)]:
                # Radio mínimo fisiológico
                if r_hijo < 8e-6:
                    continue

                d    = dir_anatomica(niv, sistema)
                d   *= signo if np.random.random() > 0.5 else 1
                L    = longitud_nivel(niv, r_hijo)
                p_h  = padre.pos + d * L

                if not dentro(p_h):
                    p_h = proyectar(p_h)

                hijo = Nodo(p_h, r_hijo, niv,
                            padre=padre, sistema=sistema)
                padre.hijos.append(hijo)
                todos.append(hijo)
                siguiente.append(hijo)

        if siguiente:
            r_prom = np.mean([n.radio for n in siguiente])*1e6
            r_min  = min(n.radio for n in siguiente)*1e6
            print(f"    Nivel {niv}: {len(siguiente)} nodos | "
                  f"radio prom={r_prom:.1f} µm | "
                  f"radio min={r_min:.1f} µm")

        nivel_actual = siguiente
        if not nivel_actual:
            break

    return todos

# ── SISTEMA COLECTOR ───────────────────────────────────────────────────────────
def generar_colector():
    nodos = []
    pelvis = Nodo(HILIO_URI * 0.4, 4e-3, 0, sistema='col')
    nodos.append(pelvis)

    for i, offset_y in enumerate([-0.010, 0.010]):
        cm = Nodo(pelvis.pos + np.array([0.008, offset_y, 0.002*i]),
                  2.5e-3, 1, padre=pelvis, sistema='col')
        pelvis.hijos.append(cm); nodos.append(cm)

        for j in range(4):
            ang  = (j - 1.5) * 0.35
            p_cn = cm.pos + np.array([
                0.010 + np.random.uniform(-0.004, 0.004),
                np.sin(ang) * 0.012,
                np.cos(ang) * 0.008])
            p_cn = proyectar(p_cn, 0.86)
            cn   = Nodo(p_cn, 1.5e-3, 2, padre=cm, sistema='col')
            cm.hijos.append(cn); nodos.append(cn)

            for k in range(5):
                d = np.array([np.random.uniform(0.4,1.0),
                              np.random.uniform(-0.4,0.4),
                              np.random.uniform(-0.4,0.4)])
                d /= np.linalg.norm(d)
                pt = cn.pos + d * np.random.uniform(0.006, 0.014)
                pt = proyectar(pt, 0.87)
                tb = Nodo(pt, 0.6e-3, 3, padre=cn, sistema='col')
                cn.hijos.append(tb); nodos.append(tb)
    return nodos

# ── EXTRAER SEGMENTOS ──────────────────────────────────────────────────────────
def extraer_segs(nodos):
    segs = []
    for n in nodos:
        for h in n.hijos:
            segs.append({
                'curva'  : spline(n.pos, h.pos, n=PUNTOS_SPLINE),
                'inicio' : n.pos.copy(),
                'fin'    : h.pos.copy(),
                'radio'  : (n.radio + h.radio) / 2,
                'nivel'  : h.nivel,
                'sistema': h.sistema,
            })
    return segs

# ── VALIDACIÓN ─────────────────────────────────────────────────────────────────
def validar(na, nv, nc, segs):
    todos = na + nv + nc
    dentro_count = sum(1 for n in todos if dentro(n.pos))

    print("\n" + "═"*56)
    print("  VALIDACIÓN FINAL — CCO v4")
    print("═"*56)
    print(f"  Nodos arteriales    : {len(na)}")
    print(f"  Nodos venosos       : {len(nv)}")
    print(f"  Nodos colectores    : {len(nc)}")
    print(f"  Total nodos         : {len(todos)}")
    print(f"  Segmentos totales   : {len(segs)}")
    print(f"  Nodos dentro riñón  : {dentro_count}/{len(todos)} "
          f"({100*dentro_count/len(todos):.1f}%)")

    vasc = [n for n in na+nv]
    if vasc:
        radios = [n.radio*1e6 for n in vasc]
        print(f"  Radio mínimo        : {min(radios):.1f} µm")
        print(f"  Radio máximo        : {max(radios):.1f} µm")
        print(f"  Radio promedio      : {np.mean(radios):.1f} µm")

    # Verificación Murray estricta
    ok = viol = 0
    for n in na + nv:
        if len(n.hijos) >= 2:
            for i in range(0, len(n.hijos)-1, 2):
                r_p  = n.radio
                r_h1 = n.hijos[i].radio
                r_h2 = n.hijos[i+1].radio if i+1 < len(n.hijos) else r_h1
                if verificar_murray(r_p, r_h1, r_h2):
                    ok += 1
                else:
                    viol += 1

    total_bif = ok + viol
    if total_bif > 0:
        print(f"  Bifurcaciones verif.: {total_bif}")
        print(f"  Murray OK           : {ok} ({100*ok/total_bif:.1f}%)")
        print(f"  Murray violaciones  : {viol} ({100*viol/total_bif:.1f}%)")
    print("═"*56 + "\n")

# ── EXPORTAR CSV ───────────────────────────────────────────────────────────────
def exportar(segs, ruta):
    with open(ruta, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['id','sistema','nivel',
                    'x1_mm','y1_mm','z1_mm',
                    'x2_mm','y2_mm','z2_mm','radio_um'])
        for i,s in enumerate(segs):
            p0=s['inicio']*1000; p1=s['fin']*1000
            w.writerow([i,s['sistema'],s['nivel'],
                        *[round(v,4) for v in p0],
                        *[round(v,4) for v in p1],
                        round(s['radio']*1e6,2)])
    print(f"  CSV v4 → {ruta}")

# ── VISUALIZACIÓN ──────────────────────────────────────────────────────────────
def visualizar(segs, ruta):
    fig = plt.figure(figsize=(15,11), facecolor='#080F1E')
    ax  = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('#080F1E')

    # Elipsoide de referencia
    u=np.linspace(0,2*np.pi,40); v=np.linspace(0,np.pi,25)
    ax.plot_wireframe(
        EJE_X*np.outer(np.cos(u),np.sin(v))*1000,
        EJE_Y*np.outer(np.sin(u),np.sin(v))*1000,
        EJE_Z*np.outer(np.ones(len(u)),np.cos(v))*1000,
        color='#1A4A7A', alpha=0.05, linewidth=0.3)

    r_max = max((s['radio'] for s in segs
                 if s['sistema']!='col'), default=1e-3)

    for s in segs:
        c = s['curva']*1000
        t = s['radio']/r_max
        if s['sistema']=='art':
            col=(0.95, 0.08+0.35*(1-t), 0.08+0.25*(1-t), 0.90)
            lw=max(0.3, t*4.5)
        elif s['sistema']=='ven':
            col=(0.05+0.2*(1-t), 0.15+0.2*(1-t), 0.90, 0.85)
            lw=max(0.3, t*4.0)
        else:
            col=(0.95, 0.80, 0.10, 0.85)
            lw=max(0.3, t*2.5)
        ax.plot(c[:,0],c[:,1],c[:,2], color=col, linewidth=lw)

    # Puntos hilio
    ax.scatter(*HILIO_ART*1000,color='#FF4444',s=120,
               zorder=5,label='Hilio arterial')
    ax.scatter(*HILIO_VEN*1000,color='#4488FF',s=120,
               zorder=5,label='Hilio venoso')
    ax.scatter(*HILIO_URI*1000,color='#FFD700',s=80,
               zorder=5,label='Pelvis renal')

    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([0],[0],color=(0.95,0.08,0.08,0.9),lw=2,
               label='Sistema arterial'),
        Line2D([0],[0],color=(0.05,0.15,0.90,0.85),lw=2,
               label='Sistema venoso'),
        Line2D([0],[0],color=(0.95,0.80,0.10,0.85),lw=2,
               label='Sistema colector'),
    ], facecolor='#1A3A5C', labelcolor='white',
       fontsize=9, loc='upper right')

    ax.set_xlabel('X — Longitud (mm)',color='#AED6F1',fontsize=8)
    ax.set_ylabel('Y — Ancho (mm)',   color='#AED6F1',fontsize=8)
    ax.set_zlabel('Z — Grosor (mm)',  color='#AED6F1',fontsize=8)
    ax.tick_params(colors='#566573',labelsize=7)
    for pane in [ax.xaxis.pane,ax.yaxis.pane,ax.zaxis.pane]:
        pane.fill=False; pane.set_edgecolor('#1A4A7A')
    ax.grid(True,color='#1A4A7A',alpha=0.12)
    ax.set_xlim(-60,60); ax.set_ylim(-35,35); ax.set_zlim(-30,30)

    sa=sum(1 for s in segs if s['sistema']=='art')
    sv=sum(1 for s in segs if s['sistema']=='ven')
    sc=sum(1 for s in segs if s['sistema']=='col')
    plt.title(
        'Bio-Kidney AI 2026 — Vascularización Completa CCO v4\n'
        f'Arterial: {sa} seg  ·  Venoso: {sv} seg  '
        f'·  Colector: {sc} seg  ·  Murray α={ALPHA} estricto',
        color='white',fontsize=10,pad=15)
    plt.tight_layout()
    plt.savefig(ruta,dpi=150,bbox_inches='tight',facecolor='#080F1E')
    plt.show()
    print(f"  Imagen → {ruta}")

# ── MAIN ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    SALIDA = os.path.expanduser(
        "~/Escritorio/BioKidney-AI/02_vascular_cco/")

    print("\n"+"═"*56)
    print("  BIO-KIDNEY AI 2026 — GENERADOR CCO v4")
    print("  Murray estricto · Jerarquía anatómica · 3 sistemas")
    print("═"*56)

    print("\n  [1/5] Sistema arterial...")
    na = generar_sistema_murray(
        HILIO_ART, R_ARTERIA, 'art', N_NIVELES)

    print("\n  [2/5] Sistema venoso...")
    nv = generar_sistema_murray(
        HILIO_VEN, R_VENA, 'ven', N_NIVELES)

    print("\n  [3/5] Sistema colector...")
    nc = generar_colector()

    print("\n  [4/5] Extrayendo segmentos con spline...")
    sa = extraer_segs(na)
    sv = extraer_segs(nv)
    sc = extraer_segs(nc)
    todos = sa + sv + sc
    print(f"  Arterial: {len(sa)} · Venoso: {len(sv)} · "
          f"Colector: {len(sc)} · Total: {len(todos)}")

    validar(na, nv, nc, todos)

    print("  [5/5] Exportando...")
    exportar(todos, SALIDA+"arbol_vascular_cco_v4.csv")
    visualizar(todos, SALIDA+"arbol_vascular_cco_v4.png")

