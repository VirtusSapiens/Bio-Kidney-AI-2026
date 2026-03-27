import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import csv
import os

# ── PARÁMETROS GLOBALES ────────────────────────────────────────────────────────
ALPHA_MURRAY   = 3.0
VISCOSIDAD     = 3.5e-3
PUNTOS_SPLINE  = 10
SEMILLA        = 42

# Dimensiones elipsoide renal (metros)
EJE_X = 0.055
EJE_Y = 0.030
EJE_Z = 0.025

# Puntos de entrada/salida vasculares en el hilio
HILIO_ARTERIA = np.array([0.0, -EJE_Y * 0.85,  0.003])
HILIO_VENA    = np.array([0.0, -EJE_Y * 0.85, -0.003])
HILIO_URETER  = np.array([0.0, -EJE_Y * 0.85,  0.0  ])

# Jerarquía anatómica — 6 niveles por sistema
JERARQUIA_ARTERIAL = [
    # (nombre,                radio_mm,  n_ramas, nivel)
    ("Arteria renal",          0.500,     1,       0),
    ("Arterias segmentarias",  0.300,     5,       1),
    ("Arterias interlobares",  0.180,    10,       2),
    ("Arterias arcuatas",      0.100,    20,       3),
    ("Interlobulillares",      0.055,    40,       4),
    ("Arteriolas aferentes",   0.020,   120,       5),
]

JERARQUIA_VENOSA = [
    ("Vena renal",             0.600,     1,       0),
    ("Venas segmentarias",     0.350,     5,       1),
    ("Venas interlobares",     0.200,    10,       2),
    ("Venas arcuatas",         0.120,    20,       3),
    ("Interlobulillares v.",   0.065,    40,       4),
    ("Vénulas eferentes",      0.025,   120,       5),
]

np.random.seed(SEMILLA)

# ── DOMINIO ELIPSOIDE ──────────────────────────────────────────────────────────
def dentro_rinon(p, margen=0.95):
    return (p[0]/EJE_X)**2 + (p[1]/EJE_Y)**2 + (p[2]/EJE_Z)**2 <= margen**2

def proyectar_interior(p, margen=0.92):
    d = np.sqrt((p[0]/EJE_X)**2 + (p[1]/EJE_Y)**2 + (p[2]/EJE_Z)**2)
    if d > margen:
        escala = margen / d
        return np.array([p[0]*escala, p[1]*escala, p[2]*escala])
    return p.copy()

def generar_puntos_en_rinon(n, zona_y=None, margen=0.88):
    puntos  = []
    intentos = 0
    while len(puntos) < n and intentos < n * 100:
        p = np.array([
            np.random.uniform(-EJE_X * margen, EJE_X * margen),
            np.random.uniform(-EJE_Y * margen, EJE_Y * margen),
            np.random.uniform(-EJE_Z * margen, EJE_Z * margen),
        ])
        if zona_y is not None:
            p[1] = abs(p[1]) * zona_y
        if dentro_rinon(p, margen):
            puntos.append(p)
        intentos += 1
    return puntos

# ── SPLINE BEZIER ──────────────────────────────────────────────────────────────
def spline_bezier(p0, p1, curvatura=0.2, n=10):
    eje  = p1 - p0
    long = np.linalg.norm(eje)
    if long < 1e-10:
        return np.array([p0, p1])
    eje_n = eje / long

    perp = np.array([-eje_n[1], eje_n[0], 0.0])
    if np.linalg.norm(perp) < 1e-10:
        perp = np.array([0.0, -eje_n[2], eje_n[1]])
    perp  = perp / np.linalg.norm(perp)
    perp2 = np.cross(eje_n, perp)

    desv1  = long * curvatura * np.random.uniform(-1, 1)
    desv2  = long * curvatura * 0.5 * np.random.uniform(-1, 1)
    p_ctrl = (p0 + p1) / 2 + perp * desv1 + perp2 * desv2

    if not dentro_rinon(p_ctrl):
        p_ctrl = (p0 + p1) / 2

    t     = np.linspace(0, 1, n)
    curva = (np.outer((1-t)**2, p0) +
             np.outer(2*(1-t)*t, p_ctrl) +
             np.outer(t**2, p1))
    return curva

# ── NODO VASCULAR ──────────────────────────────────────────────────────────────
class Nodo:
    _id = 0
    def __init__(self, pos, radio, nivel=0, padre=None, sistema='arterial'):
        self.id      = Nodo._id
        Nodo._id    += 1
        self.pos     = np.array(pos, dtype=float)
        self.radio   = radio
        self.nivel   = nivel
        self.padre   = padre
        self.hijos   = []
        self.sistema = sistema

def radio_murray(r_padre, asim=0.0):
    f  = 0.5 + asim
    r1 = r_padre * max(f, 0.1)**(1.0/ALPHA_MURRAY)
    r2 = r_padre * max(1-f, 0.1)**(1.0/ALPHA_MURRAY)
    return r1, r2

# ── GENERADOR JERÁRQUICO ───────────────────────────────────────────────────────
def generar_sistema(origen, jerarquia, sistema='arterial', invertir_y=False):
    """
    Genera un árbol vascular jerárquico con 6 niveles anatómicos.
    """
    Nodo._id = 0
    raiz     = Nodo(origen, jerarquia[0][1] * 1e-3, 0, sistema=sistema)
    nodos    = [raiz]
    hojas    = [raiz]  # nodos sin hijos en el nivel actual

    for nivel in range(1, len(jerarquia)):
        nombre, radio_mm, n_ramas, _ = jerarquia[nivel]
        radio   = radio_mm * 1e-3
        nuevas_hojas = []

        # Distribuir los n_ramas entre las hojas del nivel anterior
        ramas_por_hoja = max(1, n_ramas // len(hojas))

        for hoja in hojas:
            for _ in range(ramas_por_hoja):
                # Dirección de crecimiento — hacia el interior del riñón
                dir_base = np.array([
                    np.random.uniform(-1, 1),
                    np.random.uniform(0.2, 1.0) * (1 if not invertir_y else -1),
                    np.random.uniform(-1, 1)
                ])
                dir_base = dir_base / np.linalg.norm(dir_base)

                # Longitud proporcional al nivel
                factor_long = [1.0, 0.7, 0.5, 0.35, 0.2, 0.12]
                long = (EJE_X * factor_long[nivel] *
                        np.random.uniform(0.6, 1.4))

                p_hijo = hoja.pos + dir_base * long
                p_hijo = proyectar_interior(p_hijo, 0.90)

                asim  = np.random.uniform(-0.15, 0.15)
                r1, r2 = radio_murray(hoja.radio, asim)
                r_usar = max(r1, radio * 0.8)

                hijo = Nodo(p_hijo, r_usar, nivel,
                            padre=hoja, sistema=sistema)
                hoja.hijos.append(hijo)
                nodos.append(hijo)
                nuevas_hojas.append(hijo)

                # Segunda rama (bifurcación)
                if nivel < 5 and np.random.random() > 0.3:
                    dir2 = np.array([
                        np.random.uniform(-1, 1),
                        np.random.uniform(0.1, 0.8) * (1 if not invertir_y else -1),
                        np.random.uniform(-1, 1)
                    ])
                    dir2   = dir2 / np.linalg.norm(dir2)
                    p_hij2 = hoja.pos + dir2 * long * 0.8
                    p_hij2 = proyectar_interior(p_hij2, 0.90)

                    hijo2  = Nodo(p_hij2, max(r2, radio * 0.7),
                                  nivel, padre=hoja, sistema=sistema)
                    hoja.hijos.append(hijo2)
                    nodos.append(hijo2)
                    nuevas_hojas.append(hijo2)

        hojas = nuevas_hojas
        print(f"    Nivel {nivel} ({nombre}): "
              f"{len(nuevas_hojas)} ramas, "
              f"radio ~{radio*1e6:.0f} µm")

    return nodos

# ── SISTEMA COLECTOR ───────────────────────────────────────────────────────────
def generar_sistema_colector():
    """
    Sistema colector urinario:
    Pelvis renal → Cálices mayores → Cálices menores → Tubos colectores
    """
    nodos = []

    # Pelvis renal — zona central-medial
    pelvis = Nodo(HILIO_URETER * 0.5, 0.004, 0, sistema='colector')
    nodos.append(pelvis)

    # 2 Cálices mayores
    for i, dy in enumerate([-0.012, 0.012]):
        caliz_mayor = Nodo(
            pelvis.pos + np.array([0.005, dy, 0.002*i]),
            0.003, 1, padre=pelvis, sistema='colector')
        pelvis.hijos.append(caliz_mayor)
        nodos.append(caliz_mayor)

        # 3-4 Cálices menores por cáliz mayor
        for j in range(3):
            angulo = (j - 1) * 0.4
            p_cm   = caliz_mayor.pos + np.array([
                0.008 + np.random.uniform(-0.003, 0.003),
                np.sin(angulo) * 0.010,
                np.cos(angulo) * 0.008
            ])
            p_cm = proyectar_interior(p_cm, 0.85)
            caliz_menor = Nodo(p_cm, 0.002, 2,
                               padre=caliz_mayor, sistema='colector')
            caliz_mayor.hijos.append(caliz_menor)
            nodos.append(caliz_menor)

            # 4-5 Tubos colectores por cáliz menor
            for k in range(4):
                dir_t = np.array([
                    np.random.uniform(0.3, 1.0),
                    np.random.uniform(-0.5, 0.5),
                    np.random.uniform(-0.5, 0.5)
                ])
                dir_t  = dir_t / np.linalg.norm(dir_t)
                p_tubo = caliz_menor.pos + dir_t * np.random.uniform(0.005, 0.012)
                p_tubo = proyectar_interior(p_tubo, 0.88)
                tubo   = Nodo(p_tubo, 0.0008, 3,
                              padre=caliz_menor, sistema='colector')
                caliz_menor.hijos.append(tubo)
                nodos.append(tubo)

    return nodos

# ── EXTRAER SEGMENTOS ──────────────────────────────────────────────────────────
def extraer_segmentos(nodos):
    segs = []
    for n in nodos:
        for h in n.hijos:
            curva = spline_bezier(n.pos, h.pos,
                                  curvatura=0.18,
                                  n=PUNTOS_SPLINE)
            segs.append({
                'curva'  : curva,
                'inicio' : n.pos.copy(),
                'fin'    : h.pos.copy(),
                'radio'  : (n.radio + h.radio) / 2,
                'nivel'  : h.nivel,
                'sistema': h.sistema,
            })
    return segs

# ── VALIDACIÓN ─────────────────────────────────────────────────────────────────
def validar(nodos_art, nodos_ven, nodos_col, segs):
    print("\n" + "═"*54)
    print("  VALIDACIÓN — ÁRBOL VASCULAR CCO v3")
    print("═"*54)

    todos = nodos_art + nodos_ven + nodos_col
    dentro = sum(1 for n in todos if dentro_rinon(n.pos))
    print(f"  Nodos arteriales    : {len(nodos_art)}")
    print(f"  Nodos venosos       : {len(nodos_ven)}")
    print(f"  Nodos colectores    : {len(nodos_col)}")
    print(f"  Total nodos         : {len(todos)}")
    print(f"  Segmentos totales   : {len(segs)}")
    print(f"  Nodos dentro riñón  : {dentro}/{len(todos)} "
          f"({100*dentro/len(todos):.1f}%)")

    radios = [n.radio*1e6 for n in todos if n.sistema != 'colector']
    if radios:
        print(f"  Radio mínimo vasc.  : {min(radios):.1f} µm")
        print(f"  Radio máximo vasc.  : {max(radios):.1f} µm")
        print(f"  Radio promedio      : {np.mean(radios):.1f} µm")

    # Murray
    viol = bif = 0
    for n in nodos_art + nodos_ven:
        if len(n.hijos) >= 2:
            bif += 1
            lhs  = n.radio**ALPHA_MURRAY
            rhs  = sum(h.radio**ALPHA_MURRAY for h in n.hijos[:2])
            if abs(lhs - rhs)/max(lhs, 1e-20) > 0.05:
                viol += 1
    print(f"  Bifurcaciones       : {bif}")
    if bif > 0:
        print(f"  Cumplimiento Murray : {100*(1-viol/bif):.1f}%")
    print("═"*54 + "\n")

# ── EXPORTAR CSV ───────────────────────────────────────────────────────────────
def exportar_csv(segs, ruta):
    with open(ruta, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['id','sistema','nivel',
                    'x1_mm','y1_mm','z1_mm',
                    'x2_mm','y2_mm','z2_mm',
                    'radio_um'])
        for i, s in enumerate(segs):
            p0 = s['inicio']*1000
            p1 = s['fin']*1000
            w.writerow([i, s['sistema'], s['nivel'],
                        *[round(v,4) for v in p0],
                        *[round(v,4) for v in p1],
                        round(s['radio']*1e6, 2)])
    print(f"  CSV v3 exportado → {ruta}")

# ── VISUALIZACIÓN ──────────────────────────────────────────────────────────────
def visualizar(segs, ruta_img):
    fig = plt.figure(figsize=(15, 11), facecolor='#0A1628')
    ax  = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('#0A1628')

    # Elipsoide de referencia
    u = np.linspace(0, 2*np.pi, 40)
    v = np.linspace(0, np.pi,   25)
    ax.plot_wireframe(
        EJE_X * np.outer(np.cos(u), np.sin(v)) * 1000,
        EJE_Y * np.outer(np.sin(u), np.sin(v)) * 1000,
        EJE_Z * np.outer(np.ones(len(u)), np.cos(v)) * 1000,
        color='#1A4A7A', alpha=0.06, linewidth=0.3)

    radio_max = max((s['radio'] for s in segs
                     if s['sistema'] != 'colector'), default=1)

    for s in segs:
        curva = s['curva'] * 1000
        t     = s['radio'] / radio_max

        if s['sistema'] == 'arterial':
            color = (0.9, 0.1 + 0.3*(1-t), 0.1 + 0.2*(1-t), 0.9)
            lw    = max(0.4, t * 4.0)
        elif s['sistema'] == 'venoso':
            color = (0.1 + 0.2*(1-t), 0.2 + 0.2*(1-t), 0.85, 0.85)
            lw    = max(0.4, t * 3.5)
        else:  # colector
            color = (0.9, 0.75, 0.1, 0.80)
            lw    = max(0.3, t * 2.0)

        ax.plot(curva[:,0], curva[:,1], curva[:,2],
                color=color, linewidth=lw)

    # Puntos de hilio
    ax.scatter(*HILIO_ARTERIA*1000, color='#E74C3C',
               s=100, zorder=5, label='Hilio arterial')
    ax.scatter(*HILIO_VENA*1000,    color='#3498DB',
               s=100, zorder=5, label='Hilio venoso')
    ax.scatter(*HILIO_URETER*1000,  color='#F1C40F',
               s=80,  zorder=5, label='Uréter / Pelvis')

    # Leyenda de sistemas
    from matplotlib.lines import Line2D
    leyenda = [
        Line2D([0],[0], color=(0.9,0.1,0.1,0.9), lw=2,
               label='Sistema arterial'),
        Line2D([0],[0], color=(0.1,0.2,0.85,0.85), lw=2,
               label='Sistema venoso'),
        Line2D([0],[0], color=(0.9,0.75,0.1,0.8), lw=2,
               label='Sistema colector'),
    ]
    ax.legend(handles=leyenda, facecolor='#1A3A5C',
              labelcolor='white', fontsize=9,
              loc='upper right')

    ax.set_xlabel('X — Longitud (mm)', color='#AED6F1', fontsize=8)
    ax.set_ylabel('Y — Ancho (mm)',    color='#AED6F1', fontsize=8)
    ax.set_zlabel('Z — Grosor (mm)',   color='#AED6F1', fontsize=8)
    ax.tick_params(colors='#566573', labelsize=7)
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.fill = False
        pane.set_edgecolor('#1A4A7A')
    ax.grid(True, color='#1A4A7A', alpha=0.15)
    ax.set_xlim(-60, 60)
    ax.set_ylim(-35, 35)
    ax.set_zlim(-30, 30)

    seg_art = sum(1 for s in segs if s['sistema']=='arterial')
    seg_ven = sum(1 for s in segs if s['sistema']=='venoso')
    seg_col = sum(1 for s in segs if s['sistema']=='colector')

    plt.title(
        'Bio-Kidney AI 2026 — Vascularización Completa CCO v3\n'
        f'Arterial: {seg_art} seg  ·  Venoso: {seg_ven} seg  '
        f'·  Colector: {seg_col} seg  ·  '
        f'Ley de Murray α={ALPHA_MURRAY}',
        color='white', fontsize=10, pad=15)

    plt.tight_layout()
    plt.savefig(ruta_img, dpi=150, bbox_inches='tight',
                facecolor='#0A1628')
    plt.show()
    print(f"  Imagen guardada → {ruta_img}")

# ── MAIN ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    SALIDA = os.path.expanduser(
        "~/Escritorio/BioKidney-AI/02_vascular_cco/")

    print("\n" + "═"*54)
    print("  BIO-KIDNEY AI 2026 — GENERADOR CCO v3")
    print("  Tres sistemas: Arterial · Venoso · Colector")
    print("  Jerarquía anatómica de 6 niveles")
    print("═"*54)

    print("\n  [1/5] Sistema arterial...")
    nodos_art = generar_sistema(
        HILIO_ARTERIA, JERARQUIA_ARTERIAL,
        sistema='arterial', invertir_y=False)

    print("\n  [2/5] Sistema venoso...")
    nodos_ven = generar_sistema(
        HILIO_VENA, JERARQUIA_VENOSA,
        sistema='venoso', invertir_y=False)

    print("\n  [3/5] Sistema colector...")
    nodos_col = generar_sistema_colector()

    print("\n  [4/5] Extrayendo segmentos con spline...")
    segs_art = extraer_segmentos(nodos_art)
    segs_ven = extraer_segmentos(nodos_ven)
    segs_col = extraer_segmentos(nodos_col)
    todos_segs = segs_art + segs_ven + segs_col

    print(f"  Segmentos arteriales : {len(segs_art)}")
    print(f"  Segmentos venosos    : {len(segs_ven)}")
    print(f"  Segmentos colectores : {len(segs_col)}")

    validar(nodos_art, nodos_ven, nodos_col, todos_segs)

    print("  [5/5] Exportando...")
    exportar_csv(todos_segs,
                 SALIDA + "arbol_vascular_cco_v3.csv")
    visualizar(todos_segs,
               SALIDA + "arbol_vascular_cco_v3.png")
