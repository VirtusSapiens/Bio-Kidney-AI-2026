import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import csv
import os

# ── PARÁMETROS ─────────────────────────────────────────────────────────────────
RADIO_RAIZ      = 500e-6   # 500 µm — arteria renal principal
VISCOSIDAD      = 3.5e-3   # Pa·s
PRESION         = 13300    # Pa (~100 mmHg)
N_TERMINALES    = 120      # puntos terminales a alcanzar
RADIO_MIN       = 15e-6    # 15 µm — límite capilar
ALPHA_MURRAY    = 3.0
PUNTOS_SPLINE   = 8        # puntos de curvatura por segmento
SEMILLA         = 42

# Dimensiones elipsoide renal (metros)
EJE_X = 0.055   # 5.5 cm — eje mayor (longitud)
EJE_Y = 0.030   # 3.0 cm — eje medio (ancho)
EJE_Z = 0.025   # 2.5 cm — eje menor (grosor)

# Hilio renal — entrada vascular (lado medial)
HILIO = np.array([0.0, -EJE_Y * 0.85, 0.0])

np.random.seed(SEMILLA)

# ── DOMINIO ELIPSOIDE ──────────────────────────────────────────────────────────
def dentro_rinon(punto):
    """Verifica si un punto está dentro del elipsoide renal"""
    x, y, z = punto
    return (x/EJE_X)**2 + (y/EJE_Y)**2 + (z/EJE_Z)**2 <= 1.0

def proyectar_al_rinon(punto):
    """Si el punto está fuera, lo proyecta al interior del elipsoide"""
    x, y, z = punto
    d = np.sqrt((x/EJE_X)**2 + (y/EJE_Y)**2 + (z/EJE_Z)**2)
    if d > 1.0:
        factor = 0.95 / d
        return np.array([x * factor * EJE_X / abs(x) if x != 0 else 0,
                         y * factor * EJE_Y / abs(y) if y != 0 else 0,
                         z * factor * EJE_Z / abs(z) if z != 0 else 0])
    return punto

def generar_terminales_en_rinon(n):
    """
    Genera N puntos terminales distribuidos uniformemente
    dentro del elipsoide renal usando rechazo de muestras.
    Simula los glomérulos (puntos de filtración).
    """
    terminales = []
    intentos   = 0
    while len(terminales) < n and intentos < n * 50:
        # Muestra aleatoria en el cubo contenedor
        p = np.array([
            np.random.uniform(-EJE_X * 0.92, EJE_X * 0.92),
            np.random.uniform(-EJE_Y * 0.92, EJE_Y * 0.92),
            np.random.uniform(-EJE_Z * 0.92, EJE_Z * 0.92)
        ])
        if dentro_rinon(p):
            # Evitar zona del hilio
            if np.linalg.norm(p - HILIO) > 0.008:
                terminales.append(p)
        intentos += 1
    return np.array(terminales)

# ── SPLINE CÚBICO ──────────────────────────────────────────────────────────────
def spline_segmento(p0, p1, curvatura=0.3, n_puntos=8):
    """
    Genera una curva suave entre dos puntos con desviación lateral.
    Simula la curvatura natural de los vasos sanguíneos.
    """
    # Vector perpendicular aleatorio para la curvatura
    eje = p1 - p0
    long = np.linalg.norm(eje)
    if long < 1e-10:
        return np.array([p0, p1])

    eje_n = eje / long
    perp  = np.array([-eje_n[1], eje_n[0], 0])
    if np.linalg.norm(perp) < 1e-10:
        perp = np.array([0, -eje_n[2], eje_n[1]])
    perp = perp / np.linalg.norm(perp)

    # Punto de control para la curva (Bezier cuadrático)
    desv    = long * curvatura * np.random.uniform(-1, 1)
    desv2   = long * curvatura * 0.5 * np.random.uniform(-1, 1)
    p_ctrl  = (p0 + p1) / 2 + perp * desv
    p_ctrl += np.cross(eje_n, perp) * desv2

    # Mantener dentro del riñón
    if not dentro_rinon(p_ctrl):
        p_ctrl = (p0 + p1) / 2

    # Interpolar curva de Bezier
    t      = np.linspace(0, 1, n_puntos)
    curva  = np.outer((1-t)**2, p0) + \
             np.outer(2*(1-t)*t, p_ctrl) + \
             np.outer(t**2, p1)
    return curva

# ── NODO DEL ÁRBOL ────────────────────────────────────────────────────────────
class Nodo:
    _id = 0
    def __init__(self, posicion, radio, padre=None):
        self.id       = Nodo._id
        Nodo._id     += 1
        self.pos      = np.array(posicion)
        self.radio    = radio
        self.padre    = padre
        self.hijos    = []
        self.terminal = False

# ── GENERADOR CCO v2 ───────────────────────────────────────────────────────────
def radio_hijos_murray(r_padre, asimetria=0.0):
    f  = 0.5 + asimetria
    r1 = r_padre * f**(1.0/ALPHA_MURRAY)
    r2 = r_padre * (1-f)**(1.0/ALPHA_MURRAY)
    return r1, r2

def nodo_mas_cercano(arbol, punto):
    """Encuentra el nodo del árbol más cercano a un punto terminal"""
    min_dist = float('inf')
    mejor    = None
    for nodo in arbol:
        if not nodo.terminal and nodo.radio > RADIO_MIN * 2:
            d = np.linalg.norm(nodo.pos - punto)
            if d < min_dist:
                min_dist = d
                mejor    = nodo
    return mejor

def generar_arbol_cco_v2():
    Nodo._id = 0

    # Raíz en el hilio renal
    raiz = Nodo(HILIO, RADIO_RAIZ)
    arbol = [raiz]

    # Generar puntos terminales (glomérulos)
    terminales = generar_terminales_en_rinon(N_TERMINALES)
    print(f"  Terminales (glomérulos) generados: {len(terminales)}")

    # Primer segmento hacia el centro del riñón
    centro = np.array([0.0, 0.0, 0.0])
    r1, r2 = radio_hijos_murray(raiz.radio, 0.1)
    hijo_izq = Nodo(centro + np.array([0.01, 0.005, 0.002]), r1, raiz)
    hijo_der = Nodo(centro + np.array([0.01, -0.005, -0.002]), r2, raiz)
    raiz.hijos = [hijo_izq, hijo_der]
    arbol.extend([hijo_izq, hijo_der])

    # Conectar cada terminal al nodo más cercano
    np.random.shuffle(terminales)
    for i, terminal in enumerate(terminales):
        padre = nodo_mas_cercano(arbol, terminal)
        if padre is None:
            continue

        r_padre = padre.radio
        asim    = np.random.uniform(-0.15, 0.15)
        r1, r2  = radio_hijos_murray(r_padre, asim)

        if r1 < RADIO_MIN or r2 < RADIO_MIN:
            # Crear solo un hijo terminal
            if r1 >= RADIO_MIN:
                hijo = Nodo(terminal, r1, padre)
                hijo.terminal = True
                padre.hijos.append(hijo)
                arbol.append(hijo)
            continue

        # Punto de bifurcación entre padre y terminal
        t_bif  = np.random.uniform(0.3, 0.7)
        p_bif  = padre.pos + t_bif * (terminal - padre.pos)

        # Mantener bifurcación dentro del riñón
        if not dentro_rinon(p_bif):
            p_bif = proyectar_al_rinon(p_bif)

        nodo_bif = Nodo(p_bif, r_padre, padre)
        padre.hijos.append(nodo_bif)
        arbol.append(nodo_bif)

        # Hijo 1 — hacia el terminal
        h1 = Nodo(terminal, r1, nodo_bif)
        h1.terminal = True

        # Hijo 2 — dirección aleatoria dentro del riñón
        dir_alt = np.random.randn(3)
        dir_alt = dir_alt / np.linalg.norm(dir_alt)
        long2   = r2 * np.random.uniform(50, 120)
        p2      = p_bif + dir_alt * long2

        # Proyectar al riñón si sale
        if not dentro_rinon(p2):
            p2 = p_bif + dir_alt * long2 * 0.5
            if not dentro_rinon(p2):
                p2 = p_bif

        h2 = Nodo(p2, r2, nodo_bif)
        nodo_bif.hijos = [h1, h2]
        arbol.extend([h1, h2])

        if (i + 1) % 20 == 0:
            print(f"  → Procesando terminal {i+1}/{len(terminales)}...")

    return arbol, terminales

# ── EXTRAER SEGMENTOS ──────────────────────────────────────────────────────────
def extraer_segmentos_con_spline(arbol):
    """
    Extrae pares (padre→hijo) con curvas spline para visualización.
    """
    segmentos = []
    for nodo in arbol:
        for hijo in nodo.hijos:
            curva = spline_segmento(
                nodo.pos, hijo.pos,
                curvatura=0.25,
                n_puntos=PUNTOS_SPLINE)
            segmentos.append({
                'curva'  : curva,
                'radio'  : (nodo.radio + hijo.radio) / 2,
                'inicio' : nodo.pos,
                'fin'    : hijo.pos,
                'radio_inicio': nodo.radio,
                'radio_fin'   : hijo.radio,
            })
    return segmentos

# ── VALIDACIÓN ─────────────────────────────────────────────────────────────────
def validar_arbol(arbol, segmentos):
    print("\n" + "═"*52)
    print("  VALIDACIÓN — ÁRBOL VASCULAR CCO v2")
    print("═"*52)
    print(f"  Nodos totales       : {len(arbol)}")
    print(f"  Segmentos con spline: {len(segmentos)}")

    terminales = [n for n in arbol if n.terminal]
    print(f"  Terminales (glomér.): {len(terminales)}")

    radios = [n.radio * 1e6 for n in arbol]
    print(f"  Radio mínimo        : {min(radios):.1f} µm")
    print(f"  Radio máximo        : {max(radios):.1f} µm")
    print(f"  Radio promedio      : {np.mean(radios):.1f} µm")

    # Verificar que los nodos estén dentro del riñón
    dentro = sum(1 for n in arbol if dentro_rinon(n.pos))
    print(f"  Nodos dentro riñón  : {dentro}/{len(arbol)} "
          f"({100*dentro/len(arbol):.1f}%)")

    # Verificar Murray
    viol = 0
    bif  = 0
    for n in arbol:
        if len(n.hijos) == 2:
            bif += 1
            lhs  = n.radio**ALPHA_MURRAY
            rhs  = sum(h.radio**ALPHA_MURRAY for h in n.hijos)
            if abs(lhs - rhs) / lhs > 0.02:
                viol += 1
    print(f"  Bifurcaciones       : {bif}")
    print(f"  Cumplimiento Murray : "
          f"{100*(1-viol/max(bif,1)):.1f}%")
    print("═"*52 + "\n")

# ── EXPORTAR CSV ───────────────────────────────────────────────────────────────
def exportar_csv_v2(segmentos, ruta):
    with open(ruta, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['id', 'x1_mm', 'y1_mm', 'z1_mm',
                         'x2_mm', 'y2_mm', 'z2_mm',
                         'radio_um', 'longitud_mm'])
        for i, seg in enumerate(segmentos):
            p0   = seg['inicio'] * 1000
            p1   = seg['fin']    * 1000
            long = np.linalg.norm(seg['fin'] - seg['inicio']) * 1000
            writer.writerow([
                i,
                *[round(v, 4) for v in p0],
                *[round(v, 4) for v in p1],
                round(seg['radio'] * 1e6, 2),
                round(long, 3)
            ])
    print(f"  CSV v2 exportado → {ruta}")

# ── VISUALIZACIÓN ──────────────────────────────────────────────────────────────
def visualizar_v2(arbol, segmentos, terminales, ruta_img):
    fig = plt.figure(figsize=(14, 10), facecolor='#0D2137')
    ax  = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('#0D2137')

    radio_max = max(s['radio'] for s in segmentos)
    radio_min = min(s['radio'] for s in segmentos)

    # Dibujar elipsoide renal de referencia (wireframe)
    u = np.linspace(0, 2*np.pi, 30)
    v = np.linspace(0, np.pi, 20)
    xe = EJE_X * np.outer(np.cos(u), np.sin(v))
    ye = EJE_Y * np.outer(np.sin(u), np.sin(v))
    ze = EJE_Z * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_wireframe(xe*1000, ye*1000, ze*1000,
                      color='#1A4A7A', alpha=0.08,
                      linewidth=0.3)

    # Dibujar segmentos con spline y color por radio
    for seg in segmentos:
        curva = seg['curva'] * 1000
        t     = (seg['radio'] - radio_min) / max(radio_max - radio_min, 1e-10)
        # Rojo intenso = arteria grande, rosa = capilar
        color = (0.9, 0.15 + 0.5*(1-t), 0.15 + 0.4*(1-t), 0.85)
        lw    = max(0.4, t * 3.5)
        ax.plot(curva[:,0], curva[:,1], curva[:,2],
                color=color, linewidth=lw, alpha=0.85)

    # Hilio renal
    ax.scatter(*HILIO*1000, color='#F1C40F',
               s=120, zorder=5, label='Hilio renal')

    # Terminales (glomérulos)
    if len(terminales) > 0:
        pts = terminales * 1000
        ax.scatter(pts[:,0], pts[:,1], pts[:,2],
                   color='#17A589', s=12, alpha=0.7,
                   label=f'Glomérulos ({len(terminales)})')

    # Formato
    ax.set_xlabel('X — Longitud (mm)', color='#AED6F1', fontsize=9)
    ax.set_ylabel('Y — Ancho (mm)',    color='#AED6F1', fontsize=9)
    ax.set_zlabel('Z — Grosor (mm)',   color='#AED6F1', fontsize=9)
    ax.tick_params(colors='#566573', labelsize=7)
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.fill = False
        pane.set_edgecolor('#1A4A7A')
    ax.grid(True, color='#1A4A7A', alpha=0.2)
    ax.set_xlim(-60, 60)
    ax.set_ylim(-35, 35)
    ax.set_zlim(-30, 30)

    plt.title(
        'Bio-Kidney AI 2026 — Árbol Vascular CCO v2\n'
        f'Dominio elipsoide renal ({EJE_X*100:.0f}×'
        f'{EJE_Y*100:.0f}×{EJE_Z*100:.0f} cm) · '
        f'Ley de Murray α={ALPHA_MURRAY} · '
        f'{len(segmentos)} segmentos · Spline cúbico',
        color='white', fontsize=10, pad=15)
    plt.legend(facecolor='#1A4A7A', labelcolor='white',
               fontsize=9, loc='upper right')
    plt.tight_layout()
    plt.savefig(ruta_img, dpi=150, bbox_inches='tight',
                facecolor='#0D2137')
    plt.show()
    print(f"  Imagen guardada → {ruta_img}")

# ── MAIN ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    SALIDA = os.path.expanduser(
        "~/Escritorio/BioKidney-AI/02_vascular_cco/")

    print("\n" + "═"*52)
    print("  BIO-KIDNEY AI 2026 — GENERADOR CCO v2")
    print("  Dominio anatómico elipsoide renal")
    print("  Spline cúbico + distribución de glomérulos")
    print("═"*52)

    print("\n  [1/4] Generando árbol vascular...")
    arbol, terminales = generar_arbol_cco_v2()

    print("\n  [2/4] Extrayendo segmentos con spline...")
    segmentos = extraer_segmentos_con_spline(arbol)

    print("\n  [3/4] Validando...")
    validar_arbol(arbol, segmentos)

    print("  [4/4] Exportando y visualizando...")
    exportar_csv_v2(segmentos,
                    SALIDA + "arbol_vascular_cco_v2.csv")
    visualizar_v2(arbol, segmentos, terminales,
                  SALIDA + "arbol_vascular_cco_v2.png")
