import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv, os

# ── PARÁMETROS ─────────────────────────────────────────────────────────────────
ALPHA         = 3.0
SEMILLA       = 42
PUNTOS_SPLINE = 10
R_MIN         = 10e-6    # 10 µm — arteriola aferente real
R_MIN_COLECTOR= 200e-6   # 200 µm — túbulo colector

# Dimensiones riñón (m)
EJE_X, EJE_Y, EJE_Z = 0.055, 0.030, 0.025

# Hilio
HILIO_ART = np.array([ 0.0, -EJE_Y*0.85,  0.003])
HILIO_VEN = np.array([ 0.0, -EJE_Y*0.85, -0.003])
HILIO_URI = np.array([ 0.0, -EJE_Y*0.85,  0.000])

# Radios raíz
R_ARTERIA = 500e-6
R_VENA    = 600e-6

# Parámetros de crecimiento adaptativo
N_SEMILLAS_DEMANDA   = 400   # puntos de demanda (glomérulos simulados)
RADIO_COBERTURA      = 0.008 # 8 mm — radio de influencia por terminal
MAX_ITERACIONES      = 2000  # máximo de segmentos por sistema
LONGITUD_BASE        = 0.006 # 6 mm longitud base nivel 0

np.random.seed(SEMILLA)

# ── DOMINIO ELIPSOIDE ──────────────────────────────────────────────────────────
def dentro(p, m=0.93):
    return ((p[0]/EJE_X)**2 +
            (p[1]/EJE_Y)**2 +
            (p[2]/EJE_Z)**2) <= m**2

def proyectar(p, m=0.90):
    d = np.sqrt((p[0]/EJE_X)**2 +
                (p[1]/EJE_Y)**2 +
                (p[2]/EJE_Z)**2)
    if d > m:
        f = m / d
        return np.array([p[0]*f, p[1]*f, p[2]*f])
    return p.copy()

def generar_demanda(n):
    """
    Genera N puntos de demanda vascular distribuidos
    uniformemente en el elipsoide renal.
    Representan los glomérulos y tejido que necesita irrigación.
    """
    pts = []
    intentos = 0
    while len(pts) < n and intentos < n * 50:
        p = np.array([
            np.random.uniform(-EJE_X*0.90, EJE_X*0.90),
            np.random.uniform(-EJE_Y*0.90, EJE_Y*0.90),
            np.random.uniform(-EJE_Z*0.90, EJE_Z*0.90),
        ])
        if dentro(p, 0.90):
            pts.append(p)
        intentos += 1
    return np.array(pts)

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

# ── SPLINE ─────────────────────────────────────────────────────────────────────
def spline(p0, p1, curv=0.15, n=10):
    eje = p1 - p0
    L   = np.linalg.norm(eje)
    if L < 1e-10:
        return np.array([p0, p1])
    en    = eje / L
    perp  = np.array([-en[1], en[0], 0.0])
    if np.linalg.norm(perp) < 1e-10:
        perp = np.array([0., -en[2], en[1]])
    perp  = perp / np.linalg.norm(perp)
    perp2 = np.cross(en, perp)
    pc = ((p0 + p1) / 2
          + perp  * L * curv * np.random.uniform(-1, 1)
          + perp2 * L * curv * 0.5 * np.random.uniform(-1, 1))
    if not dentro(pc):
        pc = (p0 + p1) / 2
    t = np.linspace(0, 1, n)
    return (np.outer((1-t)**2, p0) +
            np.outer(2*(1-t)*t, pc) +
            np.outer(t**2, p1))

# ── COBERTURA VASCULAR ────────────────────────────────────────────────────────
def calcular_cobertura(nodos, puntos_demanda):
    """
    Para cada punto de demanda, verifica si ya tiene
    un nodo vascular dentro del radio de cobertura.
    Retorna los puntos SIN cobertura (necesitan irrigación).
    """
    sin_cobertura = []
    for pd in puntos_demanda:
        cubierto = False
        for n in nodos:
            if np.linalg.norm(n.pos - pd) < RADIO_COBERTURA:
                cubierto = True
                break
        if not cubierto:
            sin_cobertura.append(pd)
    return np.array(sin_cobertura) if sin_cobertura else np.array([]).reshape(0,3)

def nodo_mas_cercano_activo(nodos, punto, r_min):
    """
    Encuentra el nodo más cercano al punto de demanda
    que todavía tenga radio suficiente para bifurcar.
    """
    mejor     = None
    min_dist  = float('inf')
    for n in nodos:
        if n.radio > r_min * 1.5:
            d = np.linalg.norm(n.pos - punto)
            if d < min_dist:
                min_dist = d
                mejor    = n
    return mejor, min_dist

# ── LONGITUD ADAPTATIVA ───────────────────────────────────────────────────────
def longitud_adaptativa(radio, dist_objetivo):
    """
    La longitud del segmento se adapta a la distancia
    al punto de demanda y al radio del vaso.
    Nunca más larga que la distancia al objetivo.
    """
    L_murray = radio * 1e6 * 0.12 * 1e-3
    L_max    = dist_objetivo * 0.6
    L_min    = radio * 3
    return np.clip(L_murray * np.random.uniform(0.7, 1.3),
                   L_min, L_max)

# ── GENERADOR ADAPTATIVO POR DEMANDA ─────────────────────────────────────────
def generar_sistema_adaptativo(origen, r_raiz, sistema,
                                puntos_demanda):
    """
    Algoritmo CCO adaptativo:
    1. Genera árbol base jerárquico (niveles 0-4)
    2. Detecta zonas sin cobertura vascular
    3. Hace crecer nuevos vasos hacia esas zonas
    4. Repite hasta cubrir todo el dominio o llegar a R_min
    """
    Nodo._id = 0
    raiz  = Nodo(origen, r_raiz, 0, sistema=sistema)
    todos = [raiz]

    print(f"\n  {sistema.upper()} — Fase 1: árbol base...")

    # ── FASE 1: Árbol base jerárquico (niveles 0-5) ────────────────
    nivel_actual = [raiz]
    for niv in range(1, 6):
        siguiente = []
        for padre in nivel_actual:
            if padre.radio < R_MIN * 3:
                continue
            asim = np.random.uniform(-0.18, 0.18)
            r1, r2 = hijos_murray(padre.radio, asim)

            for r_h, angulo_base in [(r1, 1), (r2, -1)]:
                if r_h < R_MIN * 2:
                    continue

                # Dirección anatómica por nivel
                if niv <= 1:
                    d = np.array([np.random.uniform(0.6,1.0),
                                  np.random.uniform(-0.3,0.3),
                                  np.random.uniform(-0.2,0.2)])
                elif niv == 2:
                    ang = np.random.uniform(0, 2*np.pi)
                    d   = np.array([np.cos(ang)*0.5,
                                    np.random.uniform(-0.3,0.3),
                                    np.sin(ang)*0.7])
                elif niv == 3:
                    # Arcuatas — tangentes al elipsoide
                    rn = padre.pos / np.array([EJE_X, EJE_Y, EJE_Z])
                    rl = np.linalg.norm(rn)
                    if rl > 1e-10:
                        ru   = rn / rl
                        tang = np.array([-ru[1], ru[0],
                                         np.random.uniform(-0.2,0.2)])
                        tang = tang / max(np.linalg.norm(tang), 1e-10)
                        d    = tang * 0.65 + ru * 0.35 * angulo_base
                    else:
                        d = np.random.randn(3)
                elif niv == 4:
                    # Interlobulillares — radiales
                    rn = padre.pos / np.array([EJE_X, EJE_Y, EJE_Z])
                    rl = np.linalg.norm(rn)
                    d  = (rn/max(rl,1e-10) +
                          np.random.randn(3)*0.3) * angulo_base
                else:
                    d = np.random.randn(3)

                d = d / max(np.linalg.norm(d), 1e-10)
                L = (r_h * 1e6 *
                     [1.0,0.65,0.45,0.30,0.18,0.12][niv] *
                     np.random.uniform(0.8,1.3) * 1e-3)
                ph = padre.pos + d * L
                if not dentro(ph):
                    ph = proyectar(ph, 0.88)

                hijo = Nodo(ph, r_h, niv,
                            padre=padre, sistema=sistema)
                padre.hijos.append(hijo)
                todos.append(hijo)
                siguiente.append(hijo)

        r_min_niv = min(n.radio for n in siguiente)*1e6 if siguiente else 0
        print(f"    Nivel {niv}: {len(siguiente):4d} nodos | "
              f"r_min={r_min_niv:.1f} µm")
        nivel_actual = siguiente
        if not nivel_actual:
            break

    # ── FASE 2: Crecimiento adaptativo hacia zonas sin cobertura ──
    print(f"\n  {sistema.upper()} — Fase 2: crecimiento adaptativo...")

    iteracion    = 0
    sin_mejora   = 0
    nodos_previos = len(todos)

    while iteracion < MAX_ITERACIONES:
        # Detectar zonas sin cobertura
        sin_cob = calcular_cobertura(todos, puntos_demanda)

        if len(sin_cob) == 0:
            print(f"    ✓ Cobertura completa en iteración {iteracion}")
            break

        # Tomar el punto sin cobertura más alejado de cualquier vaso
        distancias = []
        for pc in sin_cob:
            d_min = min(np.linalg.norm(n.pos - pc) for n in todos)
            distancias.append(d_min)

        idx_objetivo = np.argmax(distancias)
        objetivo     = sin_cob[idx_objetivo]
        dist_obj     = distancias[idx_objetivo]

        # Encontrar nodo activo más cercano
        padre, dist_padre = nodo_mas_cercano_activo(
            todos, objetivo, R_MIN)

        if padre is None:
            break

        # Calcular radios hijos con Murray
        asim   = np.random.uniform(-0.15, 0.15)
        r1, r2 = hijos_murray(padre.radio, asim)

        if r1 < R_MIN and r2 < R_MIN:
            sin_mejora += 1
            if sin_mejora > 50:
                break
            iteracion += 1
            continue

        # Dirección hacia el objetivo
        dir_obj = objetivo - padre.pos
        d_norm  = np.linalg.norm(dir_obj)
        if d_norm > 1e-10:
            dir_obj = dir_obj / d_norm
        else:
            dir_obj = np.random.randn(3)
            dir_obj = dir_obj / np.linalg.norm(dir_obj)

        # Hijo 1 — hacia el objetivo
        if r1 >= R_MIN:
            L1  = longitud_adaptativa(r1, dist_padre)
            p1  = padre.pos + dir_obj * L1
            if not dentro(p1):
                p1 = proyectar(p1, 0.88)
            h1  = Nodo(p1, r1, padre.nivel+1,
                       padre=padre, sistema=sistema)
            padre.hijos.append(h1)
            todos.append(h1)

        # Hijo 2 — dirección complementaria
        if r2 >= R_MIN:
            perp = np.cross(dir_obj,
                            np.random.randn(3))
            if np.linalg.norm(perp) > 1e-10:
                perp = perp / np.linalg.norm(perp)
            else:
                perp = np.random.randn(3)
                perp = perp / np.linalg.norm(perp)

            ang2 = np.random.uniform(25, 50)
            d2   = (dir_obj * np.cos(np.radians(ang2)) +
                    perp    * np.sin(np.radians(ang2)))
            d2   = d2 / max(np.linalg.norm(d2), 1e-10)
            L2   = longitud_adaptativa(r2, dist_padre * 0.7)
            p2   = padre.pos + d2 * L2
            if not dentro(p2):
                p2 = proyectar(p2, 0.88)
            h2   = Nodo(p2, r2, padre.nivel+1,
                        padre=padre, sistema=sistema)
            padre.hijos.append(h2)
            todos.append(h2)

        sin_mejora = 0
        iteracion += 1

        if iteracion % 200 == 0:
            cob_pct = 100*(1 - len(sin_cob)/len(puntos_demanda))
            r_act   = min(n.radio for n in todos)*1e6
            print(f"    iter {iteracion:4d} | "
                  f"nodos={len(todos):5d} | "
                  f"cobertura={cob_pct:.1f}% | "
                  f"r_min={r_act:.1f} µm")

    # Reporte final fase 2
    sin_cob_final = calcular_cobertura(todos, puntos_demanda)
    cob_final = 100*(1 - len(sin_cob_final)/len(puntos_demanda))
    r_min_final = min(n.radio for n in todos)*1e6
    print(f"\n    Nodos añadidos fase 2 : {len(todos)-nodos_previos}")
    print(f"    Cobertura final       : {cob_final:.1f}%")
    print(f"    Radio mínimo final    : {r_min_final:.1f} µm")

    return todos

# ── SISTEMA COLECTOR ───────────────────────────────────────────────────────────
def generar_colector():
    nodos  = []
    pelvis = Nodo(HILIO_URI * 0.4, 4e-3, 0, sistema='col')
    nodos.append(pelvis)

    for i, dy in enumerate([-0.013, 0.013]):
        cm = Nodo(pelvis.pos + np.array([0.007, dy, 0.002*i]),
                  2.5e-3, 1, padre=pelvis, sistema='col')
        pelvis.hijos.append(cm)
        nodos.append(cm)

        for j in range(5):
            ang  = (j - 2) * 0.32
            p_cn = cm.pos + np.array([
                0.013 + np.random.uniform(-0.005, 0.005),
                np.sin(ang) * 0.014,
                np.cos(ang) * 0.010])
            p_cn = proyectar(p_cn, 0.85)
            cn   = Nodo(p_cn, 1.5e-3, 2,
                        padre=cm, sistema='col')
            cm.hijos.append(cn)
            nodos.append(cn)

            for k in range(6):
                d  = np.array([np.random.uniform(0.3,1.0),
                               np.random.uniform(-0.5,0.5),
                               np.random.uniform(-0.5,0.5)])
                d /= np.linalg.norm(d)
                pt  = cn.pos + d * np.random.uniform(0.007, 0.016)
                pt  = proyectar(pt, 0.87)
                tb  = Nodo(pt, R_MIN_COLECTOR, 3,
                           padre=cn, sistema='col')
                cn.hijos.append(tb)
                nodos.append(tb)
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
def validar(na, nv, nc, segs, demanda):
    todos    = na + nv + nc
    dentro_c = sum(1 for n in todos if dentro(n.pos))

    print("\n" + "═"*60)
    print("  VALIDACIÓN FINAL — CCO v6 (Adaptativo)")
    print("═"*60)
    print(f"  Nodos arteriales      : {len(na)}")
    print(f"  Nodos venosos         : {len(nv)}")
    print(f"  Nodos colectores      : {len(nc)}")
    print(f"  Total nodos           : {len(todos)}")
    print(f"  Segmentos totales     : {len(segs)}")
    print(f"  Nodos dentro riñón    : {dentro_c}/{len(todos)} "
          f"({100*dentro_c/len(todos):.1f}%)")

    vasc   = na + nv
    radios = [n.radio*1e6 for n in vasc]
    print(f"  Radio mínimo          : {min(radios):.2f} µm")
    print(f"  Radio máximo          : {max(radios):.1f} µm")
    print(f"  Radio promedio        : {np.mean(radios):.1f} µm")

    # Cobertura vascular
    sin_cob = calcular_cobertura(na, demanda)
    cob_pct = 100*(1 - len(sin_cob)/len(demanda))
    print(f"  Puntos de demanda     : {len(demanda)}")
    print(f"  Cobertura vascular    : {cob_pct:.1f}%")

    # Murray
    ok = viol = 0
    for n in vasc:
        if len(n.hijos) >= 2:
            for i in range(0, len(n.hijos)-1, 2):
                r1 = n.hijos[i].radio
                r2 = (n.hijos[i+1].radio
                      if i+1 < len(n.hijos) else r1)
                if verifica_murray(n.radio, r1, r2):
                    ok += 1
                else:
                    viol += 1

    total = ok + viol
    print(f"  Bifurcaciones verif.  : {total}")
    if total > 0:
        print(f"  Murray OK             : {ok} ({100*ok/total:.1f}%)")
        print(f"  Murray violaciones    : {viol} ({100*viol/total:.1f}%)")

    # Distribución por nivel
    print(f"\n  Distribución por nivel (arterial):")
    niveles = sorted(set(n.nivel for n in na))
    for niv in niveles:
        nn = [n for n in na if n.nivel == niv]
        if nn:
            r_min = min(n.radio for n in nn)*1e6
            r_max = max(n.radio for n in nn)*1e6
            print(f"    Nivel {niv:2d}: {len(nn):5d} nodos | "
                  f"{r_min:.1f}–{r_max:.1f} µm")
    print("═"*60 + "\n")

# ── EXPORTAR CSV ───────────────────────────────────────────────────────────────
def exportar(segs, ruta):
    with open(ruta, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['id','sistema','nivel',
                    'x1_mm','y1_mm','z1_mm',
                    'x2_mm','y2_mm','z2_mm',
                    'radio_um'])
        for i, s in enumerate(segs):
            p0 = s['inicio'] * 1000
            p1 = s['fin']    * 1000
            w.writerow([i, s['sistema'], s['nivel'],
                        *[round(v,4) for v in p0],
                        *[round(v,4) for v in p1],
                        round(s['radio']*1e6, 2)])
    print(f"  CSV v6 → {ruta}")

# ── VISUALIZACIÓN ──────────────────────────────────────────────────────────────
def visualizar(segs, demanda, ruta):
    fig = plt.figure(figsize=(16, 12), facecolor='#050B16')
    ax  = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('#050B16')

    # Elipsoide
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
        t  = np.clip(s['radio'] / r_max, 0, 1)

        if s['sistema'] == 'art':
            col = (0.95, 0.08+0.55*(1-t),
                   0.08+0.40*(1-t), 0.88)
            lw  = max(0.15, t * 5.5)
        elif s['sistema'] == 'ven':
            col = (0.05+0.25*(1-t),
                   0.10+0.30*(1-t), 0.95, 0.85)
            lw  = max(0.15, t * 5.0)
        else:
            col = (0.95, 0.82, 0.12, 0.88)
            lw  = max(0.25, t * 3.0)

        ax.plot(c[:,0], c[:,1], c[:,2],
                color=col, linewidth=lw)

    # Puntos de demanda (glomérulos)
    ax.scatter(demanda[:,0]*1000, demanda[:,1]*1000,
               demanda[:,2]*1000,
               color='#17A589', s=4, alpha=0.35,
               label=f'Glomérulos ({len(demanda)})')

    # Hilio
    ax.scatter(*HILIO_ART*1000, color='#FF3333',
               s=150, zorder=5, label='Hilio arterial')
    ax.scatter(*HILIO_VEN*1000, color='#3366FF',
               s=150, zorder=5, label='Hilio venoso')
    ax.scatter(*HILIO_URI*1000, color='#FFD700',
               s=100, zorder=5, label='Pelvis renal')

    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([0],[0], color=(0.95,0.08,0.08,0.9),
               lw=2.5, label='Sistema arterial'),
        Line2D([0],[0], color=(0.05,0.10,0.95,0.85),
               lw=2.5, label='Sistema venoso'),
        Line2D([0],[0], color=(0.95,0.82,0.12,0.88),
               lw=2,   label='Sistema colector'),
        Line2D([0],[0], color='#17A589',
               marker='o', markersize=4, lw=0,
               label=f'Glomérulos ({len(demanda)})'),
    ], facecolor='#0D2137', labelcolor='white',
       fontsize=9, loc='upper right')

    ax.set_xlabel('X — Longitud (mm)', color='#AED6F1', fontsize=9)
    ax.set_ylabel('Y — Ancho (mm)',    color='#AED6F1', fontsize=9)
    ax.set_zlabel('Z — Grosor (mm)',   color='#AED6F1', fontsize=9)
    ax.tick_params(colors='#566573', labelsize=7)
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.fill = False
        pane.set_edgecolor('#1A4A7A')
    ax.grid(True, color='#1A4A7A', alpha=0.08)
    ax.set_xlim(-60, 60)
    ax.set_ylim(-35, 35)
    ax.set_zlim(-30, 30)

    sa = sum(1 for s in segs if s['sistema']=='art')
    sv = sum(1 for s in segs if s['sistema']=='ven')
    sc = sum(1 for s in segs if s['sistema']=='col')
    r_min_f = min(s['radio'] for s in segs
                  if s['sistema']!='col') * 1e6

    plt.title(
        'Bio-Kidney AI 2026 — Vascularización Completa CCO v6\n'
        f'Arterial: {sa} seg  ·  Venoso: {sv} seg  '
        f'·  Colector: {sc} seg  ·  '
        f'R_min={r_min_f:.1f} µm  ·  '
        f'Murray α={ALPHA} estricto  ·  Adaptativo',
        color='white', fontsize=10, pad=15)

    plt.tight_layout()
    plt.savefig(ruta, dpi=150, bbox_inches='tight',
                facecolor='#050B16')
    plt.show()
    print(f"  Imagen → {ruta}")

# ── MAIN ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    SALIDA = os.path.expanduser(
        "~/Escritorio/BioKidney-AI/02_vascular_cco/")

    print("\n" + "═"*60)
    print("  BIO-KIDNEY AI 2026 — GENERADOR CCO v6")
    print("  Crecimiento adaptativo por demanda vascular")
    print(f"  R_min objetivo: {R_MIN*1e6:.0f} µm  ·  "
          f"Semillas demanda: {N_SEMILLAS_DEMANDA}")
    print("═"*60)

    print("\n  Generando puntos de demanda (glomérulos)...")
    demanda = generar_demanda(N_SEMILLAS_DEMANDA)
    print(f"  {len(demanda)} puntos generados dentro del elipsoide")

    print("\n  [1/5] Sistema arterial...")
    na = generar_sistema_adaptativo(
        HILIO_ART, R_ARTERIA, 'art', demanda)

    print("\n  [2/5] Sistema venoso...")
    nv = generar_sistema_adaptativo(
        HILIO_VEN, R_VENA, 'ven', demanda)

    print("\n  [3/5] Sistema colector...")
    nc = generar_colector()

    print("\n  [4/5] Extrayendo segmentos con spline...")
    sa = extraer_segs(na)
    sv = extraer_segs(nv)
    sc = extraer_segs(nc)
    todos = sa + sv + sc
    print(f"  Arterial : {len(sa):5d} segmentos")
    print(f"  Venoso   : {len(sv):5d} segmentos")
    print(f"  Colector : {len(sc):5d} segmentos")
    print(f"  Total    : {len(todos):5d} segmentos")

    validar(na, nv, nc, todos, demanda)

    print("  [5/5] Exportando...")
    exportar(todos, SALIDA + "arbol_vascular_cco_v6.csv")
    visualizar(todos, demanda,
               SALIDA + "arbol_vascular_cco_v6.png")
