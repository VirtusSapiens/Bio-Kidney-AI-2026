import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import csv
import os

# ── PARÁMETROS GLOBALES ────────────────────────────────────────────────────────
RADIO_RAIZ      = 500e-6   # 500 µm — arteria renal principal
PRESION_ENTRADA = 13300    # Pa (~100 mmHg)
VISCOSIDAD      = 3.5e-3   # Pa·s (sangre)
N_SEGMENTOS     = 80       # número de segmentos vasculares a generar
RADIO_MIN       = 10e-6    # 10 µm — límite capilar
ALPHA_MURRAY    = 3.0      # exponente Ley de Murray
DOMINIO         = 0.05     # cubo de 5 cm de lado (m)
SEMILLA         = 42       # reproducibilidad

np.random.seed(SEMILLA)

# ── ESTRUCTURA DE SEGMENTO ─────────────────────────────────────────────────────
class Segmento:
    _id = 0
    def __init__(self, inicio, fin, radio, padre=None):
        self.id     = Segmento._id
        Segmento._id += 1
        self.inicio = np.array(inicio, dtype=float)
        self.fin    = np.array(fin,    dtype=float)
        self.radio  = radio
        self.padre  = padre
        self.hijos  = []
        self.longitud = np.linalg.norm(self.fin - self.inicio)

    def reynolds(self):
        rho = 1060  # kg/m³ densidad sangre
        v   = PRESION_ENTRADA * self.radio**2 / (8 * VISCOSIDAD * max(self.longitud, 1e-6))
        return 2 * rho * v * self.radio / VISCOSIDAD

# ── LEY DE MURRAY ──────────────────────────────────────────────────────────────
def radio_hijos_murray(r_padre, asimetria=0.0):
    """
    r³_padre = r³_hijo1 + r³_hijo2
    asimetria ∈ [-0.3, 0.3] para variación natural
    """
    fraccion = 0.5 + asimetria
    r1 = r_padre * (fraccion) ** (1.0 / ALPHA_MURRAY)
    r2 = r_padre * (1.0 - fraccion) ** (1.0 / ALPHA_MURRAY)
    return r1, r2

# ── PUNTO DE BIFURCACIÓN ───────────────────────────────────────────────────────
def punto_bifurcacion(seg, t=0.6):
    """Punto a lo largo del segmento donde ocurre la bifurcación"""
    return seg.inicio + t * (seg.fin - seg.inicio)

# ── DIRECCIÓN DE HIJOS ─────────────────────────────────────────────────────────
def direcciones_hijos(dir_padre):
    """
    Genera dos direcciones hijas con ángulo de bifurcación basado en Murray.
    Ángulo óptimo ~37° para bifurcación simétrica.
    """
    angulo_base = np.radians(np.random.uniform(25, 45))
    rotacion    = np.random.uniform(0, 2 * np.pi)

    # Vector perpendicular al padre
    perp = np.array([-dir_padre[1], dir_padre[0], 0])
    if np.linalg.norm(perp) < 1e-10:
        perp = np.array([0, -dir_padre[2], dir_padre[1]])
    perp = perp / np.linalg.norm(perp)

    # Rotar perp alrededor de dir_padre
    c, s   = np.cos(rotacion), np.sin(rotacion)
    eje    = dir_padre / np.linalg.norm(dir_padre)
    perp_r = (perp * c +
              np.cross(eje, perp) * s +
              eje * np.dot(eje, perp) * (1 - c))
    perp_r = perp_r / np.linalg.norm(perp_r)

    dir1 = (dir_padre * np.cos(angulo_base) +
            perp_r    * np.sin(angulo_base))
    dir2 = (dir_padre * np.cos(angulo_base) -
            perp_r    * np.sin(angulo_base))

    return (dir1 / np.linalg.norm(dir1),
            dir2 / np.linalg.norm(dir2))

# ── GENERADOR CCO PRINCIPAL ────────────────────────────────────────────────────
def generar_arbol_cco():
    Segmento._id = 0
    segmentos = []

    # Segmento raíz — arteria renal entrando al hilio
    raiz = Segmento(
        inicio = [0.0, 0.0, 0.0],
        fin    = [0.01, 0.0, 0.0],
        radio  = RADIO_RAIZ
    )
    segmentos.append(raiz)
    cola = [raiz]

    iteracion = 0
    while iteracion < N_SEGMENTOS and cola:
        # Seleccionar segmento a bifurcar (preferir los más grandes)
        pesos = np.array([s.radio for s in cola])
        pesos = pesos / pesos.sum()
        idx   = np.random.choice(len(cola), p=pesos)
        padre = cola[idx]

        # Verificar radio mínimo
        r1, r2 = radio_hijos_murray(
            padre.radio,
            asimetria=np.random.uniform(-0.2, 0.2)
        )
        if r1 < RADIO_MIN or r2 < RADIO_MIN:
            cola.pop(idx)
            continue

        # Punto y dirección de bifurcación
        t      = np.random.uniform(0.4, 0.8)
        p_bif  = punto_bifurcacion(padre, t)
        dir_p  = padre.fin - padre.inicio
        dir_p  = dir_p / np.linalg.norm(dir_p)
        d1, d2 = direcciones_hijos(dir_p)

        # Longitud hija proporcional al radio (escala alométrica)
        long1 = max(r1 * np.random.uniform(8, 15) * 100, 2e-3)
        long2 = max(r2 * np.random.uniform(8, 15) * 100, 2e-3)

        fin1  = p_bif + d1 * long1
        fin2  = p_bif + d2 * long2

        # Mantener dentro del dominio
        fin1  = np.clip(fin1, -DOMINIO, DOMINIO)
        fin2  = np.clip(fin2, -DOMINIO, DOMINIO)

        h1 = Segmento(p_bif, fin1, r1, padre)
        h2 = Segmento(p_bif, fin2, r2, padre)
        padre.hijos = [h1, h2]

        segmentos.append(h1)
        segmentos.append(h2)
        cola.pop(idx)
        cola.append(h1)
        cola.append(h2)
        iteracion += 1

    return segmentos

# ── VALIDACIÓN ─────────────────────────────────────────────────────────────────
def validar_arbol(segmentos):
    print("\n" + "═"*50)
    print("  VALIDACIÓN DEL ÁRBOL VASCULAR CCO")
    print("═"*50)
    print(f"  Segmentos totales   : {len(segmentos)}")

    radios   = [s.radio * 1e6 for s in segmentos]
    reynolds = [s.reynolds()  for s in segmentos]

    print(f"  Radio mínimo        : {min(radios):.1f} µm")
    print(f"  Radio máximo        : {max(radios):.1f} µm")
    print(f"  Radio promedio      : {np.mean(radios):.1f} µm")
    print(f"  Reynolds máximo     : {max(reynolds):.1f}")
    print(f"  Reynolds promedio   : {np.mean(reynolds):.1f}")

    laminar = sum(1 for r in reynolds if r < 2100)
    print(f"  Flujo laminar       : {laminar}/{len(segmentos)} "
          f"({100*laminar/len(segmentos):.1f}%)")

    # Verificar Ley de Murray en bifurcaciones
    violaciones = 0
    for s in segmentos:
        if len(s.hijos) == 2:
            lhs = s.radio ** ALPHA_MURRAY
            rhs = (s.hijos[0].radio ** ALPHA_MURRAY +
                   s.hijos[1].radio ** ALPHA_MURRAY)
            if abs(lhs - rhs) / lhs > 0.01:
                violaciones += 1

    bifurcaciones = sum(1 for s in segmentos if len(s.hijos) == 2)
    print(f"  Bifurcaciones       : {bifurcaciones}")
    print(f"  Violaciones Murray  : {violaciones}")
    print(f"  Cumplimiento Murray : {100*(1-violaciones/max(bifurcaciones,1)):.1f}%")
    print("═"*50 + "\n")

# ── EXPORTAR CSV ───────────────────────────────────────────────────────────────
def exportar_csv(segmentos, ruta):
    with open(ruta, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['id', 'x1', 'y1', 'z1',
                         'x2', 'y2', 'z2',
                         'radio_um', 'longitud_mm', 'reynolds'])
        for s in segmentos:
            writer.writerow([
                s.id,
                *s.inicio * 1000,
                *s.fin    * 1000,
                round(s.radio * 1e6, 2),
                round(s.longitud * 1000, 3),
                round(s.reynolds(), 2)
            ])
    print(f"  CSV exportado → {ruta}")

# ── VISUALIZACIÓN 3D ───────────────────────────────────────────────────────────
def visualizar(segmentos, ruta_img):
    fig = plt.figure(figsize=(12, 9), facecolor='#0D2137')
    ax  = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('#0D2137')

    radio_max = max(s.radio for s in segmentos)
    radio_min = min(s.radio for s in segmentos)

    lineas  = []
    colores = []

    for s in segmentos:
        lineas.append([s.inicio * 1000, s.fin * 1000])
        # Color por radio: rojo=grande, azul=pequeño
        t = (s.radio - radio_min) / max(radio_max - radio_min, 1e-10)
        colores.append((t, 0.3, 1.0 - t, 0.85))

    col = Line3DCollection(lineas, colors=colores,
                           linewidths=[max(0.5, s.radio/radio_max*4)
                                       for s in segmentos])
    ax.add_collection3d(col)

    # Punto raíz
    ax.scatter(*segmentos[0].inicio * 1000,
               color='#F1C40F', s=80, zorder=5, label='Hilio renal')

    # Terminales
    terminales = [s for s in segmentos if not s.hijos]
    if terminales:
        pts = np.array([s.fin * 1000 for s in terminales])
        ax.scatter(pts[:,0], pts[:,1], pts[:,2],
                   color='#17A589', s=8, alpha=0.6, label='Terminales')

    # Formato
    todos = np.array([[s.inicio, s.fin] for s in segmentos]).reshape(-1, 3) * 1000
    for set_lim, datos in zip(
        [ax.set_xlim, ax.set_ylim, ax.set_zlim],
        [todos[:,0], todos[:,1], todos[:,2]]
    ):
        mn, mx = datos.min(), datos.max()
        pad = (mx - mn) * 0.1 or 1
        set_lim(mn - pad, mx + pad)

    ax.set_xlabel('X (mm)', color='white', fontsize=9)
    ax.set_ylabel('Y (mm)', color='white', fontsize=9)
    ax.set_zlabel('Z (mm)', color='white', fontsize=9)
    ax.tick_params(colors='white', labelsize=7)
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.fill = False
        pane.set_edgecolor('#1A4A7A')
    ax.grid(True, color='#1A4A7A', alpha=0.3)

    plt.title('Bio-Kidney AI 2026 — Árbol Vascular CCO\n'
              f'Ley de Murray (α={ALPHA_MURRAY}) · '
              f'{len(segmentos)} segmentos · '
              f'Raíz: {RADIO_RAIZ*1e6:.0f} µm',
              color='white', fontsize=11, pad=15)
    plt.legend(facecolor='#1A4A7A', labelcolor='white', fontsize=9)
    plt.tight_layout()
    plt.savefig(ruta_img, dpi=150, bbox_inches='tight',
                facecolor='#0D2137')
    plt.show()
    print(f"  Imagen guardada  → {ruta_img}")

# ── MAIN ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    SALIDA = os.path.expanduser(
        "~/Escritorio/BioKidney-AI/02_vascular_cco/"
    )
    os.makedirs(SALIDA, exist_ok=True)

    print("\n  Generando árbol vascular CCO...")
    print(f"  Segmentos objetivo : {N_SEGMENTOS}")
    print(f"  Radio raíz        : {RADIO_RAIZ*1e6:.0f} µm")
    print(f"  Radio mínimo      : {RADIO_MIN*1e6:.0f} µm")
    print(f"  Ley de Murray α   : {ALPHA_MURRAY}")

    segmentos = generar_arbol_cco()
    validar_arbol(segmentos)
    exportar_csv(segmentos,
                 SALIDA + "arbol_vascular_cco.csv")
    visualizar(segmentos,
               SALIDA + "arbol_vascular_cco.png")
