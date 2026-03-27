import numpy as np
import matplotlib.pyplot as plt

def simular_reprogramacion_celular():
    print("--- BIO-KIDNEY AI 2026: SIMULADOR DE CINÉTICA iPSC V2.0 ---")
    
    # Tiempo de simulación (30 días de cultivo)
    dias = np.linspace(0, 30, 300)
    
    # Parámetros de transformación (Tasas de éxito basadas en Vibe Coding)
    tasa_activacion = 0.15   # Velocidad de activación de Oct4/Sox2
    tasa_diferenciacion = 0.1 # Velocidad de especialización renal
    
    # Modelado de Marcadores
    # 1. Pérdida de identidad original (Somatic Memory)
    identidad_original = 100 * np.exp(-tasa_activacion * dias)
    
    # 2. Ascenso de Pluripotencia (Estado iPSC)
    pluripotencia = 100 * (1 - np.exp(-tasa_activacion * dias)) * np.exp(-tasa_diferenciacion * (dias - 10).clip(min=0))
    
    # 3. Especialización Renal (Podocitos/Túbulos)
    especializacion_renal = 100 * (1 - np.exp(-tasa_diferenciacion * (dias - 15).clip(min=0)))

    # Graficación de la maduración celular
    plt.figure(figsize=(12, 7))
    plt.plot(dias, identidad_original, label='Memoria Somática (Célula Original)', color='gray', linestyle='--')
    plt.plot(dias, pluripotencia, label='Estado Pluripotente (iPSC)', color='gold', linewidth=2)
    plt.plot(dias, especializacion_renal, label='Funcionalidad Renal (Podocitos)', color='blue', linewidth=3)
    
    # Hitos del Proyecto
    plt.axvline(x=15, color='red', linestyle=':', label='Inicio de Diferenciación')
    plt.axvspan(21, 30, color='green', alpha=0.1, label='Ventana de Bioimpresión')
    
    plt.title("Cinética de Maduración Celular: Bio-Kidney AI 2026", fontsize=14, fontweight='bold')
    plt.xlabel("Días de Cultivo en Biorreactor", fontsize=12)
    plt.ylabel("Porcentaje de Expresión Fenotípica (%)", fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Veredicto de la IA
    eficiencia_final = especializacion_renal[-1]
    print(f"[*] Eficiencia de diferenciación proyectada a día 30: {eficiencia_final:.2f}%")
    
    if eficiencia_final > 80:
        print("[VERDICTO]: Cultivo VALIDADO. Células aptas para carga en biotinta SWIFT.")
    else:
        print("[ALERTA]: Baja eficiencia. Ajustar cóctel de citoquinas.")
        
    plt.show()

if __name__ == "__main__":
    simular_reprogramacion_celular()
