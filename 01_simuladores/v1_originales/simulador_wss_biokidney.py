import numpy as np
import matplotlib.pyplot as plt

def simular_hemodinamica():
    viscosidad_sangre = 0.0035
    radio_vaso = 50e-6
    presiones = np.linspace(80, 120, 100) * 133.322
    longitud = 0.05
    wss_pascal = (radio_vaso * presiones) / (2 * longitud)
    wss_dynas = wss_pascal * 10
    plt.figure(figsize=(10, 6))
    plt.plot(presiones / 133.322, wss_dynas, color='red', label='WSS Calculado')
    plt.axhspan(10, 70, color='green', alpha=0.2, label='Zona de Seguridad')
    plt.title("Validación Hemodinámica: Bio-Kidney AI 2026")
    plt.xlabel("Presión (mmHg)")
    plt.ylabel("Wall Shear Stress (dyn/cm²)")
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    simular_hemodinamica()
