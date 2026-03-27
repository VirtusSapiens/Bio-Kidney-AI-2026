import sys
import os
import numpy as np
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, 
    QHBoxLayout, QGridLayout, QLabel, QSlider, 
    QFrame, QGroupBox, QSizePolicy
)
from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QFont, QColor, QPainter, QBrush, QPen
from PyQt6.QtCharts import QChart, QChartView, QLineSeries, QValueAxis, QAreaSeries

# ── PALETA DE COLORES BIOKIDNEY AI ─────────────────────────────────────────────
AZUL_OSCURO  = "#0D1117"
PANEL        = "#161B22"
AZUL_CLARO   = "#00D4FF"
VERDE        = "#00FF88"
ROJO         = "#FF4444"
AMARILLO     = "#FFD700"
BLANCO       = "#E6EDF3"
GRIS         = "#8B949E"

ESTILO_GUI = f"""
QMainWindow {{ background-color: {AZUL_OSCURO}; }}
QFrame#MainFrame {{ background-color: {AZUL_OSCURO}; }}
QGroupBox {{
    font-size: 13px; font-weight: bold; color: {AZUL_CLARO};
    border: 1.5px solid #30363D; border-radius: 8px;
    margin-top: 15px; padding-top: 10px;
}}
QGroupBox::title {{ subcontrol-origin: margin; left: 12px; padding: 0 6px; }}
QLabel {{ color: {BLANCO}; font-size: 12px; }}
QSlider::groove:horizontal {{
    height: 6px; background: #30363D; border-radius: 3px;
}}
QSlider::handle:horizontal {{
    background: {AZUL_CLARO}; width: 16px; height: 16px;
    margin: -5px 0; border-radius: 8px;
}}
"""

# ── WIDGET DE RESULTADO ESTILIZADO ─────────────────────────────────────────────
class ResultadoPanel(QFrame):
    def __init__(self, titulo, unidad, color_val=AZUL_CLARO):
        super().__init__()
        self.color_val = color_val
        self.setStyleSheet(f"""
            QFrame {{ 
                background: {PANEL}; 
                border-radius: 10px; 
                border: 1.5px solid #30363D; 
            }}
        """)
        layout = QVBoxLayout(self)
        self.lbl_titulo = QLabel(titulo)
        self.lbl_titulo.setStyleSheet(f"color: {GRIS}; font-size: 11px; border: none;")
        self.lbl_titulo.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.lbl_valor = QLabel("0.00")
        self.lbl_valor.setStyleSheet(f"color: {color_val}; font-size: 28px; font-weight: bold; border: none;")
        self.lbl_valor.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.lbl_unidad = QLabel(unidad)
        self.lbl_unidad.setStyleSheet(f"color: {GRIS}; font-size: 10px; border: none;")
        self.lbl_unidad.setAlignment(Qt.AlignmentFlag.AlignCenter)

        layout.addWidget(self.lbl_titulo)
        layout.addWidget(self.lbl_valor)
        layout.addWidget(self.lbl_unidad)

    def actualizar(self, valor, color=None):
        self.lbl_valor.setText(f"{valor:.2f}")
        if color:
            self.lbl_valor.setStyleSheet(f"color: {color}; font-size: 28px; font-weight: bold; border: none;")

# ── APLICACIÓN PRINCIPAL ───────────────────────────────────────────────────────
class FiltracionGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Bio-Kidney AI 2026 — Simulador de Filtración Glomerular")
        self.resize(1100, 750)
        self.setStyleSheet(ESTILO_GUI)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)
        main_layout.setContentsMargins(20, 20, 20, 20)
        main_layout.setSpacing(20)

        # ── COLUMNA IZQUIERDA: CONTROLES ──────────────────────────────────────
        controles_panel = QFrame()
        controles_panel.setFixedWidth(320)
        controles_layout = QVBoxLayout(controles_panel)
        controles_layout.setContentsMargins(0, 0, 0, 0)

        header = QLabel("LABORATORIO DE FILTRACIÓN")
        header.setStyleSheet(f"font-size: 18px; font-weight: bold; color: {AZUL_CLARO};")
        controles_layout.addWidget(header)

        # Grupo Presiones Hidrostáticas
        grp_presion = QGroupBox("Presiones Hidrostáticas (mmHg)")
        l_pres = QVBoxLayout(grp_presion)
        
        self.sld_pgc = self.crear_control(l_pres, "P_capilar (Pgc):", 40, 80, 60)
        self.sld_pbs = self.crear_control(l_pres, "P_bowman (Pbs):", 5, 25, 15)
        controles_layout.addWidget(grp_presion)

        # Grupo Presión Oncótica
        grp_oncotica = QGroupBox("Presión Oncótica (mmHg)")
        l_onc = QVBoxLayout(grp_oncotica)
        self.sld_pi = self.crear_control(l_onc, "π_entrada (πgc):", 20, 35, 28)
        controles_layout.addWidget(grp_oncotica)

        # Grupo Permeabilidad
        grp_perm = QGroupBox("Permeabilidad Glomerular")
        l_perm = QVBoxLayout(grp_perm)
        self.sld_kf = self.crear_control(l_perm, "Kf (mL/min/mmHg):", 1, 6, 3) # escala x10 en slider
        controles_layout.addWidget(grp_perm)

        controles_layout.addStretch()

        # Branding
        footer = QLabel("Bio-Kidney AI 2026\nVirtusSapiens")
        footer.setStyleSheet(f"color: {GRIS}; font-size: 10px;")
        controles_layout.addWidget(footer)

        main_layout.addWidget(controles_panel)

        # ── COLUMNA DERECHA: DASHBOARD ───────────────────────────────────────
        dashboard_layout = QVBoxLayout()
        
        # Resultados Superiores
        res_layout = QHBoxLayout()
        self.res_tfg = ResultadoPanel("TFG TOTAL", "mL / min", color_val=VERDE)
        self.res_ff  = ResultadoPanel("FRACCIÓN FILTRACIÓN", "%", color_val=AMARILLO)
        self.res_dp  = ResultadoPanel("ΔP STARLING MEDIO", "mmHg", color_val=AZUL_CLARO)
        res_layout.addWidget(self.res_tfg)
        res_layout.addWidget(self.res_ff)
        res_layout.addWidget(self.res_dp)
        dashboard_layout.addLayout(res_layout)

        # Gráfico Starling
        self.chart_view = self.crear_grafico_starling()
        dashboard_layout.addWidget(self.chart_view)

        main_layout.addLayout(dashboard_layout)

        # Conectar señales
        self.sld_pgc.valueChanged.connect(self.actualizar_simulacion)
        self.sld_pbs.valueChanged.connect(self.actualizar_simulacion)
        self.sld_pi.valueChanged.connect(self.actualizar_simulacion)
        self.sld_kf.valueChanged.connect(self.actualizar_simulacion)

        # Iniciar
        self.actualizar_simulacion()

    def crear_control(self, layout, texto, min_v, max_v, def_v):
        l_label = QHBoxLayout()
        lbl_name = QLabel(texto)
        lbl_val  = QLabel(str(def_v))
        lbl_val.setFixedWidth(40)
        lbl_val.setAlignment(Qt.AlignmentFlag.AlignRight)
        lbl_val.setStyleSheet(f"color: {AZUL_CLARO}; font-weight: bold;")
        l_label.addWidget(lbl_name)
        l_label.addWidget(lbl_val)
        layout.addLayout(l_label)

        slider = QSlider(Qt.Orientation.Horizontal)
        slider.setRange(min_v * 10, max_v * 10)
        slider.setValue(def_v * 10)
        slider.valueChanged.connect(lambda v: lbl_val.setText(f"{v/10:.1f}"))
        layout.addWidget(slider)
        return slider

    def crear_grafico_starling(self):
        chart = QChart()
        chart.setBackgroundBrush(QBrush(QColor(PANEL)))
        chart.setTitle("DINÁMICA DE STARLING — MODELO DEEN")
        chart.setTitleFont(QFont("monospace", 10, QFont.Weight.Bold))
        chart.setTitleBrush(QBrush(QColor(AZUL_CLARO)))
        
        self.series_dp = QLineSeries()
        self.series_dp.setName("ΔP Starling Neto")
        pen = QPen(QColor(AZUL_CLARO))
        pen.setWidth(3)
        self.series_dp.setPen(pen)

        self.series_pi = QLineSeries()
        self.series_pi.setName("π Oncótica (Deen)")
        pen_pi = QPen(QColor("#A855F7"))
        pen_pi.setWidth(2)
        pen_pi.setStyle(Qt.PenStyle.DashLine)
        self.series_pi.setPen(pen_pi)

        chart.addSeries(self.series_dp)
        chart.addSeries(self.series_pi)

        axis_x = QValueAxis()
        axis_x.setRange(0, 100)
        axis_x.setTitleText("Posición en Capilar (%)")
        axis_x.setLabelsColor(QColor(GRIS))
        axis_x.setTitleBrush(QBrush(QColor(GRIS)))
        chart.addAxis(axis_x, Qt.AlignmentFlag.AlignBottom)
        self.series_dp.attachAxis(axis_x)
        self.series_pi.attachAxis(axis_x)

        axis_y = QValueAxis()
        axis_y.setRange(0, 50)
        axis_y.setTitleText("Presión (mmHg)")
        axis_y.setLabelsColor(QColor(GRIS))
        axis_y.setTitleBrush(QBrush(QColor(GRIS)))
        chart.addAxis(axis_y, Qt.AlignmentFlag.AlignLeft)
        self.series_dp.attachAxis(axis_y)
        self.series_pi.attachAxis(axis_y)

        chart.legend().setVisible(True)
        chart.legend().setAlignment(Qt.AlignmentFlag.AlignTop)
        chart.legend().setLabelColor(QColor(BLANCO))

        view = QChartView(chart)
        view.setRenderHint(QPainter.RenderHint.Antialiasing)
        view.setStyleSheet(f"background-color: {PANEL}; border: none;")
        return view

    def actualizar_simulacion(self):
        # Obtener valores de sliders
        Pgc = self.sld_pgc.value() / 10.0
        Pbs = self.sld_pbs.value() / 10.0
        pi_0 = self.sld_pi.value() / 10.0
        Kf_total = self.sld_kf.value() / 10.0 # mL/min/mmHg

        # Modelo Físico Simplificado (Basado en el simulador completo)
        x_norm = np.linspace(0, 1, 50)
        FF_est = 0.20
        delta_p_list = []
        pi_list = []
        
        # Integración local
        kf_glom = Kf_total / 1_000_000 * 1e6 # nL/min/mmHg
        tfg_glom = 0
        
        for xn in x_norm:
            # pi(x) = pi_0 / (1 - FF*x)
            pi_x = pi_0 / (1 - FF_est * xn)
            dp_x = max((Pgc - Pbs) - (pi_x - 0), 0)
            
            delta_p_list.append(dp_x)
            pi_list.append(pi_x)
            
            dq = kf_glom * dp_x * (1.0 / len(x_norm))
            tfg_glom += dq
            
            # Ajustar FF para la siguiente iteración
            q_plasma = 0.625 # nL/min normal
            FF_est = min(tfg_glom / q_plasma, 0.8)

        # Escalar a riñón completo
        tfg_total = tfg_glom * 1e-6 * 1_000_000
        ff_final = FF_est * 100
        dp_medio = np.mean(delta_p_list)

        # Actualizar UI
        color_tfg = VERDE if tfg_total >= 60 else (AMARILLO if tfg_total >= 40 else ROJO)
        self.res_tfg.actualizar(tfg_total, color=color_tfg)
        self.res_ff.actualizar(ff_final)
        self.res_dp.actualizar(dp_medio)

        # Actualizar Gráfico
        points_dp = []
        points_pi = []
        for i, (dp, pi) in enumerate(zip(delta_p_list, pi_list)):
            points_dp.append(Qt.QtCore.QPointF(i * 2, dp))
            points_pi.append(Qt.QtCore.QPointF(i * 2, pi))
        
        self.series_dp.replace(points_dp)
        self.series_pi.replace(points_pi)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = FiltracionGUI()
    window.show()
    sys.exit(app.exec())
