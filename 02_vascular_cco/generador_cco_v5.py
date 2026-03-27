import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv, os

# ── PARÁMETROS ─────────────────────────────────────────────────────────────────
ALPHA         = 3.0
SEMILLA       = 42
PUNTOS_SPLINE = 10
N_NIVELES     = 8       # 8 niveles — llega hasta arteriolas finas

# Dimensiones riñón (m)
EJE_X, EJE_Y, EJE_Z = 0.055, 0.030, 0.025

# Radios raíz
R_ARTERIA = 500e-6
R_VENA    = 600e-6

# Hilio
HILIO_ART = np.array([ 0.0, -EJE_Y*0.85,  0.003])
HILIO_VEN = np.array([ 0.0, -EJE_Y*0.85, -0.003])
HILIO_URI = np.array([ 0.0, -EJE_Y*0.85,  0.000])

# Radio mínimo fisiológico (arteriola aferente)
R_MIN = 10e-6

np.random.seed(SEMILLA)

# ── DOMINIO ELIPSOIDE ──────────────────────────────────────────────────────────
def dentro(p, m=0.93):
    return (p[0]/EJE_X)**2+(p[1]/EJE_Y)**2+(p[2]/EJE_Z)**2 <= m**2

def proyectar(p, m=0.90):
    d = np.sqrt((p[0]/EJE_X)**2+(p[1]/EJE_Y)**2+(p[2]/EJE_Z)**2)
    if d > m:
        f = m / d
        return np.array([p[0]*f, p[1]*f, p[2]*f])
    return p.copy()

# ── MURRAY ESTRICTO ────────────────────────────────────────────────────────────
def hijos_murray(r_padre, asim=0.0):
    f  = np.clip(0.5 + asim, 0.25, 0.75)
    r1 = r_padre * f**(1.0/ALPHA)
    r2 = r_padre * (1.0-f)**(1.0/ALPHA)
    return r1, r2

def verifica_murray(rp, r1, r2):
    return abs(r1**ALPHA + r2**ALPHA - rp**ALPHA) / rp**ALPHA < 0.005

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

# ── SPLINE BEZIER ──────────────────────────────────────────────────────────────
def spline(p0, p1, curv=0.15, n=10):
    eje = p1-p0; L = np.linalg.norm(eje)
    if L < 1e-10:
        return np.array([p0, p1])
    en = eje/L
    perp = np.array([-en[1], en[0], 0.0])
    if np.linalg.norm(perp) < 1e-10:
        perp = np.array([0., -en[2], en[1]])
    perp  = perp / np.linalg.norm(perp)
    perp2 = np.cross(en, perp)
    pc = ((p0+p1)/2
          + perp  * L * curv * np.random.uniform(-1, 1)
          + perp2 * L * curv * 0.5 * np.random.uniform(-1, 1))
    if not dentro(pc):
        pc = (p0+p1)/2
    t = np.linspace(0, 1, n)
    return (np.outer((1-t)**2, p0) +
            np.outer(2*(1-t)*t, pc) +
            np.outer(t**2, p1))

# ── DIRECCIÓN ANATÓMICA v5 ────────────────────────────────────────────────────
def dir_anatomica_v5(nivel, sistema, padre_pos):
    """
    Direcciones anatómicas corregidas para distribución polar:

    Niv 0-1 : desde hilio hacia interior (+X dominante)
    Niv 2   : arterias interlobares — divergen hacia polos
              (distribución en abanico en X y Z)
    Niv 3   : arcuatas — arcos paralelos a la cápsula
              (tangentes al elipsoide)
    Niv 4   : interlobulillares — radiales hacia cápsula
              (dirección desde centro hacia superficie)
    Niv 5-6 : arteriolas — locales, pequeña desviación
    Niv 7-8 : capilares terminales — muy locales
    """
    if nivel <= 1:
        # Hacia el interior desde el hilio
        base = np.array([
            np.random.uniform(0.5, 1.0),
            np.random.uniform(-0.3, 0.3),
            np.random.uniform(-0.2, 0.2)
        ])

    elif nivel == 2:
        # Distribución en abanico hacia los polos
        # Los polos del riñón están en ±X y ±Z
        polo = np.random.choice(['X+', 'X-', 'Z+', 'Z-',
                                 'XZ1', 'XZ2', 'XZ3', 'XZ4'])
        if polo == 'X+':
            base = np.array([0.8,  np.random.uniform(-0.3,0.3),
                             np.random.uniform(-0.2,0.2)])
        elif polo == 'X-':
            base = np.array([-0.7, np.random.uniform(-0.3,0.3),
                             np.random.uniform(-0.2,0.2)])
        elif polo == 'Z+':
            base = np.array([np.random.uniform(0.2,0.5),
                             np.random.uniform(-0.2,0.2), 0.8])
        elif polo == 'Z-':
            base = np.array([np.random.uniform(0.2,0.5),
                             np.random.uniform(-0.2,0.2), -0.8])
        else:
            ang  = np.random.uniform(0, 2*np.pi)
            base = np.array([np.cos(ang)*0.6,
                             np.random.uniform(-0.2,0.2),
                             np.sin(ang)*0.6])

    elif nivel == 3:
        # Arcuatas — tangentes a la superficie del elipsoide
        # Dirección perpendicular al radio desde el centro
        r_norm = padre_pos / np.array([EJE_X, EJE_Y, EJE_Z])
        r_norm_len = np.linalg.norm(r_norm)
        if r_norm_len > 1e-10:
            r_unit = r_norm / r_norm_len
            # Vector tangente — perpendicular al radio
            tang = np.array([-r_unit[1], r_unit[0],
                             np.random.uniform(-0.3, 0.3)])
            tang = tang / max(np.linalg.norm(tang), 1e-10)
            # Mezcla tangente + algo radial para no ser perfectamente plano
            base = tang * 0.7 + r_unit * 0.3
        else:
            ang  = np.random.uniform(0, 2*np.pi)
            base = np.array([np.cos(ang), np.sin(ang), 0.1])

    elif nivel == 4:
        # Interlobulillares — radiales desde el centro hacia la cápsula
        # Dirección desde padre hacia la superficie más cercana
        r_norm = padre_pos / np.array([EJE_X, EJE_Y, EJE_Z])
        r_len  = np.linalg.norm(r_norm)
        if r_len > 1e-10:
            # Hacia afuera (cápsula) con algo de desviación
            base = (r_norm/r_len +
                    np.random.randn(3)*0.25)
        else:
            base = np.random.randn(3)

    else:
        # Niveles 5-8: ramificación local muy fina
        # Pequeña desviación desde la dirección padre
        if padre_pos is not None and nivel >= 2:
            dir_local = np.random.randn(3)
            dir_local = dir_local / np.linalg.norm(dir_local)
            base = dir_local
        else:
            base = np.random.randn(3)

    if sistema == 'ven':
        # Venoso es casi espejo del arterial en Y
        base = base * np.array([1.0, 1.0, -1.0])

    n = np.linalg.norm(base)
    return base/n if n > 1e-10 else np.array([1., 0., 0.])

def longitud_nivel_v5(nivel, radio):
    """
    Longitud por nivel — escala alométrica estricta.
    Longitud ∝ radio (relación Murray-Poiseuille)
    """
    factores = [1.0, 0.70, 0.50, 0.35, 0.22, 0.14, 0.08, 0.05, 0.03]
    f   = factores[min(nivel, len(factores)-1)]
    # Longitud base en mm, convertida a metros
    L_mm = radio * 1e6 * f * np.random.uniform(0.8, 1.3)
    return max(L_mm * 1e-3, radio * 5)

# ── GENERADOR CON MURRAY ESTRICTO ─────────────────────────────────────────────
def generar_sistema_v5(origen, r_raiz, sistema, n_niveles):
    Nodo._id = 0
    raiz  = Nodo(origen, r_raiz, 0, sistema=sistema)
    todos = [raiz]
    nivel_actual = [raiz]

    print(f"\n  {sistema.upper()} — propagación Murray v5:")
    print(f"    Nivel 0: 1 nodo | radio = {r_raiz*1e6:.0f} µm")

    for niv in range(1, n_niveles+1):
        siguiente = []

        for padre in nivel_actual:
            if padre.radio < R_MIN * 1.5:
                continue

            asim = np.random.uniform(-0.20, 0.20)
            r1, r2 = hijos_murray(padre.radio, asim)

            for r_hijo in [r1, r2]:
                if r_hijo < R_MIN:
                    continue

                d = dir_anatomica_v5(niv, sistema, padre.pos)
                L = longitud_nivel_v5(niv, r_hijo)
                p_h = padre.pos + d * L

                if not dentro(p_h):
                    p_h = proyectar(p_h, 0.89)

                hijo = Nodo(p_h, r_hijo, niv,
                            padre=padre, sistema=sistema)
                padre.hijos.append(hijo)
                todos.append(hijo)
                siguiente.append(hijo)

        if siguiente:
            r_min_niv  = min(n.radio for n in siguiente)*1e6
            r_prom_niv = np.mean([n.radio for n in siguiente])*1e6
            print(f"    Nivel {niv}: {len(siguiente):4d} nodos | "
                  f"r_min={r_min_niv:6.1f} µm | "
                  f"r_prom={r_prom_niv:6.1f} µm")

        nivel_actual = siguiente
        if not nivel_actual:
            print(f"    → Radio mínimo alcanzado en nivel {niv}")
            break

    return todos

# ── SISTEMA COLECTOR ───────────────────────────────────────────────────────────
def generar_colector():
    nodos  = []
    pelvis = Nodo(HILIO_URI * 0.4, 4e-3, 0, sistema='col')
    nodos.append(pelvis)

    for i, dy in enumerate([-0.012, 0.012]):
        cm = Nodo(pelvis.pos + np.array([0.006, dy, 0.002*i]),
                  2.5e-3, 1, padre=pelvis, sistema='col')
        pelvis.hijos.append(cm); nodos.append(cm)

        for j in range(4):
            ang  = (j - 1.5) * 0.38
            p_cn = cm.pos + np.array([
                0.012 + np.random.uniform(-0.005, 0.005),
                np.sin(ang) * 0.013,
                np.cos(ang) * 0.009])
            p_cn = proyectar(p_cn, 0.85)
            cn   = Nodo(p_cn, 1.5e-3, 2, padre=cm, sistema='col')
            cm.hijos.append(cn); nodos.append(cn)

            for k in range(6):
                d  = np.array([np.random.uniform(0.3, 1.0),
                               np.random.uniform(-0.5, 0.5),
                               np.random.uniform(-0.5, 0.5)])
                d /= np.linalg.norm(d)
                pt  = cn.pos + d * np.random.uniform(0.006, 0.015)
                pt  = proyectar(pt, 0.87)
                tb  = Nodo(pt, 0.5e-3, 3,
                           padre=cn, sistema='col')
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
    todos  = na + nv + nc
    dentro_c = sum(1 for n in todos if dentro(n.pos))

    print("\n" + "═"*58)
    print("  VALIDACIÓN FINAL — CCO v5")
    print("═"*58)
    print(f"  Nodos arteriales    : {len(na)}")
    print(f"  Nodos venosos       : {len(nv)}")
    print(f"  Nodos colectores    : {len(nc)}")
    print(f"  Total nodos         : {len(todos)}")
    print(f"  Segmentos totales   : {len(segs)}")
    print(f"  Nodos dentro riñón  : {dentro_c}/{len(todos)} "
          f"({100*dentro_c/len(todos):.1f}%)")

    vasc   = na + nv
    radios = [n.radio*1e6 for n in vasc]
    print(f"  Radio mínimo        : {min(radios):.1f} µm")
    print(f"  Radio máximo        : {max(radios):.1f} µm")
    print(f"  Radio promedio      : {np.mean(radios):.1f} µm")

    # Verificación Murray
    ok = viol = 0
    for n in vasc:
        if len(n.hijos) >= 2:
            for i in range(0, len(n.hijos)-1, 2):
                r1 = n.hijos[i].radio
                r2 = n.hijos[i+1].radio if i+1 < len(n.hijos) else r1
                if verifica_murray(n.radio, r1, r2):
                    ok += 1
                else:
                    viol += 1

    total = ok + viol
    print(f"  Bifurcaciones verif.: {total}")
    if total > 0:
        print(f"  Murray OK           : {ok} ({100*ok/total:.1f}%)")
        print(f"  Murray violaciones  : {viol} ({100*viol/total:.1f}%)")

    # Distribución por nivel
    print(f"\n  Distribución por nivel (arterial):")
    for niv in range(N_NIVELES+1):
        nodos_niv = [n for n in na if n.nivel == niv]
        if nodos_niv:
            r_min = min(n.radio for n in nodos_niv)*1e6
            r_max = max(n.radio for n in nodos_niv)*1e6
            print(f"    Nivel {niv}: {len(nodos_niv):4d} nodos | "
                  f"{r_min:.1f}–{r_max:.1f} µm")
    print("═"*58 + "\n")

# ── EXPORTAR CSV ───────────────────────────────────────────────────────────────
def exportar(segs, ruta):
    with open(ruta, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['id','sistema','nivel',
                    'x1_mm','y1_mm','z1_mm',
                    'x2_mm','y2_mm','z2_mm','radio_um'])
        for i, s in enumerate(segs):
            p0 = s['inicio']*1000
            p1 = s['fin']*1000
            w.writerow([i, s['sistema'], s['nivel'],
                        *[round(v,4) for v in p0],
                        *[round(v,4) for v in p1],
                        round(s['radio']*1e6, 2)])
    print(f"  CSV v5 → {ruta}")

# ── VISUALIZACIÓN ──────────────────────────────────────────────────────────────
def visualizar(segs, ruta):
    fig = plt.figure(figsize=(16, 12), facecolor='#060D1A')
    ax  = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('#060D1A')

    # Elipsoide de referencia
    u = np.linspace(0, 2*np.pi, 50)
    v = np.linspace(0, np.pi,   30)
    ax.plot_wireframe(
        EJE_X * np.outer(np.cos(u), np.sin(v)) * 1000,
        EJE_Y * np.outer(np.sin(u), np.sin(v)) * 1000,
        EJE_Z * np.outer(np.ones(len(u)), np.cos(v)) * 1000,
        color='#1A4A7A', alpha=0.05, linewidth=0.3)

    r_max = max((s['radio'] for s in segs
                 if s['sistema'] != 'col'), default=1e-3)

    for s in segs:
        c  = s['curva'] * 1000
        t  = s['radio'] / r_max
        lv = s['nivel']

        if s['sistema'] == 'art':
            # Rojo oscuro (arterias gruesas) → rosa (capilares)
            col = (0.95, 0.08+0.55*(1-t), 0.08+0.40*(1-t), 0.88)
            lw  = max(0.2, t * 5.0)
        elif s['sistema'] == 'ven':
            # Azul oscuro (venas gruesas) → celeste (vénulas)
            col = (0.05+0.30*(1-t), 0.15+0.35*(1-t), 0.92, 0.85)
            lw  = max(0.2, t * 4.5)
        else:
            # Amarillo — sistema colector
            col = (0.95, 0.82, 0.12, 0.88)
            lw  = max(0.3, t * 3.0)

        ax.plot(c[:,0], c[:,1], c[:,2],
                color=col, linewidth=lw)

    # Puntos hilio
    ax.scatter(*HILIO_ART*1000, color='#FF3333',
               s=150, zorder=5, label='Hilio arterial')
    ax.scatter(*HILIO_VEN*1000, color='#3366FF',
               s=150, zorder=5, label='Hilio venoso')
    ax.scatter(*HILIO_URI*1000, color='#FFD700',
               s=100, zorder=5, label='Pelvis renal')

    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([0],[0], color=(0.95,0.08,0.08,0.9), lw=2.5,
               label='Sistema arterial'),
        Line2D([0],[0], color=(0.05,0.15,0.92,0.85), lw=2.5,
               label='Sistema venoso'),
        Line2D([0],[0], color=(0.95,0.82,0.12,0.88), lw=2,
               label='Sistema colector'),
    ], facecolor='#0D2137', labelcolor='white',
       fontsize=10, loc='upper right')

    ax.set_xlabel('X — Longitud (mm)', color='#AED6F1', fontsize=9)
    ax.set_ylabel('Y — Ancho (mm)',    color='#AED6F1', fontsize=9)
    ax.set_zlabel('Z — Grosor (mm)',   color='#AED6F1', fontsize=9)
    ax.tick_params(colors='#566573', labelsize=7)
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.fill = False
        pane.set_edgecolor('#1A4A7A')
    ax.grid(True, color='#1A4A7A', alpha=0.10)
    ax.set_xlim(-60, 60)
    ax.set_ylim(-35, 35)
    ax.set_zlim(-30, 30)

    sa = sum(1 for s in segs if s['sistema']=='art')
    sv = sum(1 for s in segs if s['sistema']=='ven')
    sc = sum(1 for s in segs if s['sistema']=='col')
    r_min_final = min(s['radio'] for s in segs
                      if s['sistema']!='col')*1e6

    plt.title(
        'Bio-Kidney AI 2026 — Vascularización Completa CCO v5\n'
        f'Arterial: {sa} seg  ·  Venoso: {sv} seg  '
        f'·  Colector: {sc} seg  ·  '
        f'R_min={r_min_final:.1f} µm  ·  Murray α={ALPHA} 100%',
        color='white', fontsize=10, pad=15)

    plt.tight_layout()
    plt.savefig(ruta, dpi=150, bbox_inches='tight',
                facecolor='#060D1A')
    plt.show()
    print(f"  Imagen → {ruta}")

# ── MAIN ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    SALIDA = os.path.expanduser(
        "~/Escritorio/BioKidney-AI/02_vascular_cco/")

    print("\n" + "═"*58)
    print("  BIO-KIDNEY AI 2026 — GENERADOR CCO v5")
    print("  8 niveles anatómicos · R_min=10µm · Murray 100%")
    print("  Distribución polar corregida · 3 sistemas")
    print("═"*58)

    print("\n  [1/5] Sistema arterial (8 niveles)...")
    na = generar_sistema_v5(HILIO_ART, R_ARTERIA, 'art', N_NIVELES)

    print("\n  [2/5] Sistema venoso (8 niveles)...")
    nv = generar_sistema_v5(HILIO_VEN, R_VENA, 'ven', N_NIVELES)

    print("\n  [3/5] Sistema colector...")
    nc = generar_colector()

    print("\n  [4/5] Extrayendo segmentos con spline...")
    sa = extraer_segs(na)
    sv = extraer_segs(nv)
    sc = extraer_segs(nc)
    todos = sa + sv + sc
    print(f"  Arterial : {len(sa):4d} segmentos")
    print(f"  Venoso   : {len(sv):4d} segmentos")
    print(f"  Colector : {len(sc):4d} segmentos")
    print(f"  Total    : {len(todos):4d} segmentos")

    validar(na, nv, nc, todos)

    print("  [5/5] Exportando...")
    exportar(todos, SALIDA + "arbol_vascular_cco_v5.csv")
    visualizar(todos, SALIDA + "arbol_vascular_cco_v5.png")
