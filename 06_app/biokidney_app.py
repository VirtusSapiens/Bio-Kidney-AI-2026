import sys
import numpy as np
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QTabWidget, QWidget,
    QVBoxLayout, QHBoxLayout, QGridLayout, QLabel,
    QDoubleSpinBox, QPushButton, QFrame, QSlider,
    QGroupBox, QSizePolicy
)
from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QFont, QColor
from PyQt6.QtCharts import (
    QChart, QChartView, QLineSeries, QValueAxis, QAreaSeries
)
from PyQt6.QtGui import QPainter, QBrush, QPen

# ── PALETA DE COLORES ──────────────────────────────────────────────────────────
AZUL_OSCURO  = "#0D2137"
AZUL_MEDIO   = "#1A4A7A"
AZUL_CLARO   = "#2E86C1"
CIAN         = "#17A589"
ROJO         = "#C0392B"
AMARILLO     = "#F1C40F"
GRIS_FONDO   = "#F4F6F8"
BLANCO       = "#FFFFFF"

ESTILO_GLOBAL = """
QMainWindow { background-color: #F4F6F8; }
QTabWidget::pane { border: 1px solid #D5D8DC; background: #FFFFFF; border-radius: 8px; }
QTabBar::tab {
    background: #D5D8DC; color: #2C3E50;
    padding: 10px 24px; font-size: 13px; font-weight: bold;
    border-top-left-radius: 6px; border-top-right-radius: 6px;
    margin-right: 2px;
}
QTabBar::tab:selected { background: #1A4A7A; color: white; }
QTabBar::tab:hover:!selected { background: #AEB6BF; }
QGroupBox {
    font-size: 12px; font-weight: bold; color: #1A4A7A;
    border: 1.5px solid #D5D8DC; border-radius: 8px;
    margin-top: 12px; padding-top: 10px;
}
QGroupBox::title { subcontrol-origin: margin; left: 12px; padding: 0 6px; }
QDoubleSpinBox {
    border: 1.5px solid #AEB6BF; border-radius: 6px;
    padding: 6px 10px; font-size: 13px; background: white;
    min-height: 32px;
}
QDoubleSpinBox:focus { border: 1.5px solid #2E86C1; }
QPushButton {
    background: #1A4A7A; color: white;
    border-radius: 8px; padding: 10px 24px;
    font-size: 13px; font-weight: bold;
}
QPushButton:hover { background: #2E86C1; }
QPushButton:pressed { background: #0D2137; }
QSlider::groove:horizontal {
    height: 6px; background: #D5D8DC; border-radius: 3px;
}
QSlider::handle:horizontal {
    background: #2E86C1; width: 18px; height: 18px;
    margin: -6px 0; border-radius: 9px;
}
QSlider::sub-page:horizontal { background: #2E86C1; border-radius: 3px; }
QLabel { color: #2C3E50; }
"""

# ── WIDGET DE RESULTADO ────────────────────────────────────────────────────────
class ResultadoWidget(QFrame):
    def __init__(self, titulo, unidad):
        super().__init__()
        self.setStyleSheet("""
            QFrame { background: #EAF4FB; border-radius: 10px;
                     border: 1.5px solid #AED6F1; }
        """)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 12, 16, 12)

        self.lbl_titulo = QLabel(titulo)
        self.lbl_titulo.setStyleSheet("font-size: 11px; color: #566573; border: none;")
        self.lbl_titulo.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.lbl_valor = QLabel("—")
        self.lbl_valor.setStyleSheet("font-size: 26px; font-weight: bold; color: #1A4A7A; border: none;")
        self.lbl_valor.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.lbl_unidad = QLabel(unidad)
        self.lbl_unidad.setStyleSheet("font-size: 11px; color: #566573; border: none;")
        self.lbl_unidad.setAlignment(Qt.AlignmentFlag.AlignCenter)

        layout.addWidget(self.lbl_titulo)
        layout.addWidget(self.lbl_valor)
        layout.addWidget(self.lbl_unidad)

    def actualizar(self, valor, estado=None):
        self.lbl_valor.setText(f"{valor:.2f}")
        if estado == "ÓPTIMO":
            self.setStyleSheet("QFrame { background: #EAFAF1; border-radius: 10px; border: 1.5px solid #27AE60; }")
            self.lbl_valor.setStyleSheet("font-size: 26px; font-weight: bold; color: #1E8449; border: none;")
        elif estado == "CRÍTICO":
            self.setStyleSheet("QFrame { background: #FDEDEC; border-radius: 10px; border: 1.5px solid #E74C3C; }")
            self.lbl_valor.setStyleSheet("font-size: 26px; font-weight: bold; color: #C0392B; border: none;")
        else:
            self.setStyleSheet("QFrame { background: #FEF9E7; border-radius: 10px; border: 1.5px solid #F1C40F; }")
            self.lbl_valor.setStyleSheet("font-size: 26px; font-weight: bold; color: #B7950B; border: none;")

# ── BADGE DE ESTADO ────────────────────────────────────────────────────────────
class EstadoBadge(QLabel):
    def __init__(self):
        super().__init__("— Calcular para ver estado —")
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setMinimumHeight(44)
        self.setStyleSheet("""
            background: #D5D8DC; color: #566573;
            border-radius: 8px; font-size: 14px; font-weight: bold;
            padding: 8px 20px;
        """)

    def set_estado(self, estado):
        if estado == "ÓPTIMO":
            self.setText("✔  ESTADO: ÓPTIMO")
            self.setStyleSheet("background: #1E8449; color: white; border-radius: 8px; font-size: 14px; font-weight: bold; padding: 8px 20px;")
        elif estado == "CRÍTICO":
            self.setText("✖  ESTADO: CRÍTICO")
            self.setStyleSheet("background: #C0392B; color: white; border-radius: 8px; font-size: 14px; font-weight: bold; padding: 8px 20px;")
        else:
            self.setText("⚠  ESTADO: PARCIAL")
            self.setStyleSheet("background: #B7950B; color: white; border-radius: 8px; font-size: 14px; font-weight: bold; padding: 8px 20px;")

# ── TAB 1: WSS ─────────────────────────────────────────────────────────────────
class TabWSS(QWidget):
    def __init__(self):
        super().__init__()
        layout = QHBoxLayout(self)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(20)

        # Panel izquierdo — entradas
        panel_izq = QVBoxLayout()

        gb_params = QGroupBox("Parámetros de entrada")
        grid = QGridLayout(gb_params)
        grid.setSpacing(12)

        grid.addWidget(QLabel("Radio arteriola (µm):"), 0, 0)
        self.spin_radio = QDoubleSpinBox()
        self.spin_radio.setRange(10, 200)
        self.spin_radio.setValue(50)
        self.spin_radio.setSuffix(" µm")
        grid.addWidget(self.spin_radio, 0, 1)

        grid.addWidget(QLabel("Viscosidad sangre (cP):"), 1, 0)
        self.spin_visc = QDoubleSpinBox()
        self.spin_visc.setRange(1, 10)
        self.spin_visc.setValue(3.5)
        self.spin_visc.setSingleStep(0.1)
        self.spin_visc.setSuffix(" cP")
        grid.addWidget(self.spin_visc, 1, 1)

        grid.addWidget(QLabel("Velocidad flujo (mm/s):"), 2, 0)
        self.spin_vel = QDoubleSpinBox()
        self.spin_vel.setRange(0.1, 50)
        self.spin_vel.setValue(2.0)
        self.spin_vel.setSingleStep(0.5)
        self.spin_vel.setSuffix(" mm/s")
        grid.addWidget(self.spin_vel, 2, 1)

        panel_izq.addWidget(gb_params)

        gb_zona = QGroupBox("Zona fisiológica de referencia")
        zona_layout = QVBoxLayout(gb_zona)
        lbl_ref = QLabel("Rango seguro: 1 – 10 dyn/cm²\nArteriolas renales (Hagen-Poiseuille)")
        lbl_ref.setStyleSheet("font-size: 12px; color: #566573; line-height: 1.6;")
        zona_layout.addWidget(lbl_ref)
        panel_izq.addWidget(gb_zona)

        btn_calcular = QPushButton("Calcular WSS")
        btn_calcular.clicked.connect(self.calcular)
        panel_izq.addWidget(btn_calcular)
        panel_izq.addStretch()

        layout.addLayout(panel_izq, 1)

        # Panel derecho — resultados
        panel_der = QVBoxLayout()

        gb_resultados = QGroupBox("Resultados")
        res_layout = QVBoxLayout(gb_resultados)

        fila_res = QHBoxLayout()
        self.res_wss = ResultadoWidget("Wall Shear Stress", "dyn/cm²")
        self.res_flujo = ResultadoWidget("Flujo volumétrico", "µL/min")
        fila_res.addWidget(self.res_wss)
        fila_res.addWidget(self.res_flujo)
        res_layout.addLayout(fila_res)

        self.badge = EstadoBadge()
        res_layout.addWidget(self.badge)
        panel_der.addWidget(gb_resultados)

        # Gráfico
        gb_grafico = QGroupBox("WSS vs Radio arteriola")
        grafico_layout = QVBoxLayout(gb_grafico)
        self.chart_view = self._crear_grafico()
        grafico_layout.addWidget(self.chart_view)
        panel_der.addWidget(gb_grafico)

        layout.addLayout(panel_der, 2)

    def _crear_grafico(self):
        self.chart = QChart()
        self.chart.setTitle("")
        self.chart.setBackgroundBrush(QBrush(QColor("#FFFFFF")))
        self.chart.legend().setVisible(False)
        self.chart.setMargins(__import__('PyQt6.QtCore', fromlist=['QMargins']).QMargins(10, 10, 10, 10))

        self.serie_wss = QLineSeries()
        self.serie_wss.setColor(QColor(AZUL_CLARO))
        pen = QPen(QColor(AZUL_CLARO))
        pen.setWidth(3)
        self.serie_wss.setPen(pen)

        self.chart.addSeries(self.serie_wss)

        self.axis_x = QValueAxis()
        self.axis_x.setTitleText("Radio (µm)")
        self.axis_x.setRange(10, 200)
        self.axis_x.setLabelFormat("%.0f")

        self.axis_y = QValueAxis()
        self.axis_y.setTitleText("WSS (dyn/cm²)")
        self.axis_y.setRange(0, 50)
        self.axis_y.setLabelFormat("%.1f")

        self.chart.addAxis(self.axis_x, Qt.AlignmentFlag.AlignBottom)
        self.chart.addAxis(self.axis_y, Qt.AlignmentFlag.AlignLeft)
        self.serie_wss.attachAxis(self.axis_x)
        self.serie_wss.attachAxis(self.axis_y)

        view = QChartView(self.chart)
        view.setRenderHint(QPainter.RenderHint.Antialiasing)
        view.setMinimumHeight(220)
        return view

    def calcular(self):
        r_um = self.spin_radio.value()
        r_m = r_um * 1e-6
        mu = self.spin_visc.value() * 1e-3
        v = self.spin_vel.value() * 1e-3

        wss_pa = 4 * mu * v / r_m
        wss_dyn = wss_pa * 10
        flujo = np.pi * r_m**2 * v * 1e9

        estado = "ÓPTIMO" if 1 <= wss_dyn <= 10 else ("PARCIAL" if wss_dyn <= 20 else "CRÍTICO")

        self.res_wss.actualizar(wss_dyn, estado)
        self.res_flujo.actualizar(flujo, estado)
        self.badge.set_estado(estado)

        # Actualizar gráfico
        self.serie_wss.clear()
        radios = np.linspace(10, 200, 100)
        for r in radios:
            r_m2 = r * 1e-6
            w = 4 * mu * v / r_m2 * 10
            self.serie_wss.append(r, w)

        max_wss = 4 * mu * v / (10e-6) * 10
        self.axis_y.setRange(0, min(max_wss * 1.2, 100))

# ── TAB 2: iPSC ────────────────────────────────────────────────────────────────
class TabIPSC(QWidget):
    def __init__(self):
        super().__init__()
        layout = QHBoxLayout(self)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(20)

        panel_izq = QVBoxLayout()

        gb_params = QGroupBox("Parámetros de maduración")
        grid = QGridLayout(gb_params)
        grid.setSpacing(12)

        grid.addWidget(QLabel("Día de cultivo:"), 0, 0)
        self.spin_dia = QDoubleSpinBox()
        self.spin_dia.setRange(0, 60)
        self.spin_dia.setValue(21)
        self.spin_dia.setSuffix(" días")
        grid.addWidget(self.spin_dia, 0, 1)

        grid.addWidget(QLabel("Tasa de diferenciación (%):"), 1, 0)
        self.spin_tasa = QDoubleSpinBox()
        self.spin_tasa.setRange(50, 100)
        self.spin_tasa.setValue(85)
        self.spin_tasa.setSuffix(" %")
        grid.addWidget(self.spin_tasa, 1, 1)

        grid.addWidget(QLabel("Células iniciales (M):"), 2, 0)
        self.spin_cel = QDoubleSpinBox()
        self.spin_cel.setRange(1, 500)
        self.spin_cel.setValue(200)
        self.spin_cel.setSuffix(" M")
        grid.addWidget(self.spin_cel, 2, 1)

        panel_izq.addWidget(gb_params)

        gb_ventana = QGroupBox("Ventana de bioimpresión")
        v_layout = QVBoxLayout(gb_ventana)
        lbl_v = QLabel("Óptimo: Día 21 – 30\nFuncionalidad renal: ~79.8%\nRiesgo teratoma: <5% con purificación")
        lbl_v.setStyleSheet("font-size: 12px; color: #566573; line-height: 1.6;")
        v_layout.addWidget(lbl_v)
        panel_izq.addWidget(gb_ventana)

        btn = QPushButton("Calcular maduración")
        btn.clicked.connect(self.calcular)
        panel_izq.addWidget(btn)
        panel_izq.addStretch()
        layout.addLayout(panel_izq, 1)

        panel_der = QVBoxLayout()
        gb_res = QGroupBox("Resultados")
        res_l = QVBoxLayout(gb_res)
        fila = QHBoxLayout()
        self.res_func = ResultadoWidget("Funcionalidad renal", "%")
        self.res_viables = ResultadoWidget("Células viables", "M")
        fila.addWidget(self.res_func)
        fila.addWidget(self.res_viables)
        res_l.addLayout(fila)
        self.badge = EstadoBadge()
        res_l.addWidget(self.badge)
        panel_der.addWidget(gb_res)

        gb_g = QGroupBox("Curva de maduración iPSC")
        g_l = QVBoxLayout(gb_g)
        self.chart_view = self._crear_grafico()
        g_l.addWidget(self.chart_view)
        panel_der.addWidget(gb_g)
        layout.addLayout(panel_der, 2)

    def _crear_grafico(self):
        self.chart = QChart()
        self.chart.setBackgroundBrush(QBrush(QColor("#FFFFFF")))
        self.chart.legend().setVisible(False)

        self.serie = QLineSeries()
        pen = QPen(QColor(CIAN))
        pen.setWidth(3)
        self.serie.setPen(pen)
        self.chart.addSeries(self.serie)

        self.ax = QValueAxis()
        self.ax.setTitleText("Día de cultivo")
        self.ax.setRange(0, 60)
        self.ax.setLabelFormat("%.0f")

        self.ay = QValueAxis()
        self.ay.setTitleText("Funcionalidad (%)")
        self.ay.setRange(0, 100)
        self.ay.setLabelFormat("%.0f")

        self.chart.addAxis(self.ax, Qt.AlignmentFlag.AlignBottom)
        self.chart.addAxis(self.ay, Qt.AlignmentFlag.AlignLeft)
        self.serie.attachAxis(self.ax)
        self.serie.attachAxis(self.ay)

        view = QChartView(self.chart)
        view.setRenderHint(QPainter.RenderHint.Antialiasing)
        view.setMinimumHeight(220)
        return view

    def calcular(self):
        dia = self.spin_dia.value()
        tasa = self.spin_tasa.value() / 100
        cel_init = self.spin_cel.value()

        k = 0.15
        func = 100 * (1 - np.exp(-k * dia)) * tasa
        viables = cel_init * tasa * (1 - np.exp(-k * dia))

        if 21 <= dia <= 30:
            estado = "ÓPTIMO"
        elif 15 <= dia < 21 or 30 < dia <= 40:
            estado = "PARCIAL"
        else:
            estado = "CRÍTICO"

        self.res_func.actualizar(func, estado)
        self.res_viables.actualizar(viables, estado)
        self.badge.set_estado(estado)

        self.serie.clear()
        for d in np.linspace(0, 60, 200):
            f = 100 * (1 - np.exp(-k * d)) * tasa
            self.serie.append(d, f)

# ── TAB 3: SWIFT ───────────────────────────────────────────────────────────────
class TabSWIFT(QWidget):
    def __init__(self):
        super().__init__()
        layout = QHBoxLayout(self)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(20)

        panel_izq = QVBoxLayout()
        gb = QGroupBox("Parámetros de extrusión")
        grid = QGridLayout(gb)
        grid.setSpacing(12)

        grid.addWidget(QLabel("Presión (kPa):"), 0, 0)
        self.spin_pres = QDoubleSpinBox()
        self.spin_pres.setRange(5, 200)
        self.spin_pres.setValue(30)
        self.spin_pres.setSuffix(" kPa")
        grid.addWidget(self.spin_pres, 0, 1)

        grid.addWidget(QLabel("Radio boquilla (µm):"), 1, 0)
        self.spin_boquilla = QDoubleSpinBox()
        self.spin_boquilla.setRange(50, 500)
        self.spin_boquilla.setValue(200)
        self.spin_boquilla.setSuffix(" µm")
        grid.addWidget(self.spin_boquilla, 1, 1)

        grid.addWidget(QLabel("Viscosidad bioink (Pa·s):"), 2, 0)
        self.spin_visc = QDoubleSpinBox()
        self.spin_visc.setRange(0.1, 50)
        self.spin_visc.setValue(5.0)
        self.spin_visc.setSingleStep(0.5)
        self.spin_visc.setSuffix(" Pa·s")
        grid.addWidget(self.spin_visc, 2, 1)

        panel_izq.addWidget(gb)

        gb_ref = QGroupBox("Límites del protocolo SWIFT")
        ref_l = QVBoxLayout(gb_ref)
        lbl = QLabel("Presión óptima: 30 kPa\nEstrés máximo (citotoxicidad): 150 Pa\nBioink NICE: GelMA 7% + Alginato 1.5%")
        lbl.setStyleSheet("font-size: 12px; color: #566573; line-height: 1.6;")
        ref_l.addWidget(lbl)
        panel_izq.addWidget(gb_ref)

        btn = QPushButton("Calcular extrusión")
        btn.clicked.connect(self.calcular)
        panel_izq.addWidget(btn)
        panel_izq.addStretch()
        layout.addLayout(panel_izq, 1)

        panel_der = QVBoxLayout()
        gb_res = QGroupBox("Resultados")
        res_l = QVBoxLayout(gb_res)
        fila = QHBoxLayout()
        self.res_estres = ResultadoWidget("Estrés en boquilla", "Pa")
        self.res_flujo = ResultadoWidget("Flujo de extrusión", "µL/s")
        fila.addWidget(self.res_estres)
        fila.addWidget(self.res_flujo)
        res_l.addLayout(fila)
        self.badge = EstadoBadge()
        res_l.addWidget(self.badge)
        panel_der.addWidget(gb_res)

        gb_g = QGroupBox("Estrés vs Presión de extrusión")
        g_l = QVBoxLayout(gb_g)
        self.chart_view = self._crear_grafico()
        g_l.addWidget(self.chart_view)
        panel_der.addWidget(gb_g)
        layout.addLayout(panel_der, 2)

    def _crear_grafico(self):
        self.chart = QChart()
        self.chart.setBackgroundBrush(QBrush(QColor("#FFFFFF")))
        self.chart.legend().setVisible(False)

        self.serie = QLineSeries()
        pen = QPen(QColor(AMARILLO))
        pen.setWidth(3)
        self.serie.setPen(pen)

        self.serie_limite = QLineSeries()
        pen2 = QPen(QColor(ROJO))
        pen2.setWidth(2)
        pen2.setStyle(Qt.PenStyle.DashLine)
        self.serie_limite.setPen(pen2)

        self.chart.addSeries(self.serie)
        self.chart.addSeries(self.serie_limite)

        self.ax = QValueAxis()
        self.ax.setTitleText("Presión (kPa)")
        self.ax.setRange(5, 200)
        self.ax.setLabelFormat("%.0f")

        self.ay = QValueAxis()
        self.ay.setTitleText("Estrés (Pa)")
        self.ay.setRange(0, 300)
        self.ay.setLabelFormat("%.0f")

        self.chart.addAxis(self.ax, Qt.AlignmentFlag.AlignBottom)
        self.chart.addAxis(self.ay, Qt.AlignmentFlag.AlignLeft)
        self.serie.attachAxis(self.ax)
        self.serie.attachAxis(self.ay)
        self.serie_limite.attachAxis(self.ax)
        self.serie_limite.attachAxis(self.ay)

        # Línea de límite citotoxicidad
        self.serie_limite.append(5, 150)
        self.serie_limite.append(200, 150)

        view = QChartView(self.chart)
        view.setRenderHint(QPainter.RenderHint.Antialiasing)
        view.setMinimumHeight(220)
        return view

    def calcular(self):
        P = self.spin_pres.value() * 1e3
        r = self.spin_boquilla.value() * 1e-6
        mu = self.spin_visc.value()

        L = 0.05
        estres = (P * r) / (2 * L)
        flujo = (np.pi * r**4 * P) / (8 * mu * L) * 1e9

        estado = "ÓPTIMO" if estres <= 120 else ("PARCIAL" if estres <= 150 else "CRÍTICO")

        self.res_estres.actualizar(estres, estado)
        self.res_flujo.actualizar(flujo, estado)
        self.badge.set_estado(estado)

        self.serie.clear()
        for p in np.linspace(5, 200, 150):
            e = (p * 1e3 * r) / (2 * 0.01)
            self.serie.append(p, e)

# ── VENTANA PRINCIPAL ──────────────────────────────────────────────────────────
class BioKidneyApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Bio-Kidney AI 2026 — Simulador de Validación In Silico")
        self.setMinimumSize(1000, 680)
        self.setStyleSheet(ESTILO_GLOBAL)

        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QVBoxLayout(central)
        main_layout.setContentsMargins(16, 12, 16, 12)
        main_layout.setSpacing(10)

        # Header
        header = QLabel("Bio-Kidney AI 2026  ·  Simulador de Validación In Silico")
        header.setStyleSheet(f"""
            background: {AZUL_OSCURO}; color: white;
            font-size: 16px; font-weight: bold;
            padding: 14px 20px; border-radius: 8px;
        """)
        main_layout.addWidget(header)

        # Tabs
        tabs = QTabWidget()
        tabs.addTab(TabWSS(),   "  WSS — Hemodinámica  ")
        tabs.addTab(TabIPSC(),  "  iPSC — Maduración Celular  ")
        tabs.addTab(TabSWIFT(), "  SWIFT — Extrusión  ")
        main_layout.addWidget(tabs)

        # Footer
        footer = QLabel("Carlos David Moreno Cáceres · VirtusSapiens · Medellín, Colombia · 2026")
        footer.setAlignment(Qt.AlignmentFlag.AlignCenter)
        footer.setStyleSheet("font-size: 11px; color: #AEB6BF; padding: 4px;")
        main_layout.addWidget(footer)


# ── MAIN ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setApplicationName("Bio-Kidney AI 2026")
    window = BioKidneyApp()
    window.show()
    sys.exit(app.exec())
