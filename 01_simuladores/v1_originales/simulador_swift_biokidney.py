import numpy as np
import matplotlib.pyplot as plt

def simular_bioimpresion_swift():
    print("--- BIO-KIDNEY AI 2026: SIMULADOR DE EXTRUSIÓN SWIFT 2.0 ---")
    
    # Parámetros de la Boquilla (Nozzle)
    diametro_boquilla = 200e-6  # 200 micras
    radio = diametro_boquilla / 2
    
    # Rango de presiones de aire aplicadas a la jeringa (en kPa)
    presiones_kpa = np.linspace(10, 100, 50)
    presiones_pa = presiones_kpa * 1000
    
    # Reología de la Biotinta (Gelatina/Fibrina + Células)
    viscosidad_aparente = 0.5  # Pa·s (Consistencia tipo miel espesa)
    longitud_boquilla = 0.012  # 12mm de aguja
    
    # 1. Cálculo del Flujo Volumétrico (Q) - Ley de Hagen-Poiseuille
    flujo_q = (np.pi * (radio**4) * presiones_pa) / (8 * viscosidad_aparente * longitud_boquilla)
    
    # 2. Velocidad de Extrusión (m/s)
    velocidad_extrusion = flujo_q / (np.pi * radio**2)
    
    # 3. Estrés Mecánico sobre las Células (Shear Stress en la pared de la boquilla)
    # Este es el valor crítico: si supera los 150 Pa, las células mueren.
    estres_celular_pa = (4 * viscosidad_aparente * velocidad_extrusion) / radio

    # Graficación
    fig, ax1 = plt.subplots(figsize=(10, 6))

    color = 'tab:blue'
    ax1.set_xlabel('Presión de la Jeringa (kPa)', fontsize=12)
    ax1.set_ylabel('Velocidad de Salida (mm/s)', color=color, fontsize=12)
    ax1.plot(presiones_kpa, velocidad_extrusion * 1000, color=color, linewidth=2, label='Velocidad de Impresión')
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # Segundo eje para el estrés celular
    color = 'tab:red'
    ax2.set_ylabel('Estrés Mecánico en Células (Pa)', color=color, fontsize=12)
    ax2.plot(presiones_kpa, estres_celular_pa, color=color, linestyle='--', linewidth=2, label='Estrés Mecánico')
    ax2.tick_params(axis='y', labelcolor=color)

    # Límite de Viabilidad Celular (150 Pa es el estándar de seguridad)
    plt.axhspan(0, 150, color='green', alpha=0.15, label='Zona de Viabilidad (Células Vivas)')
    plt.axhline(y=150, color='darkred', linestyle=':', label='Límite de Citotoxicidad Mecánica')

    plt.title("Optimización de Extrusión SWIFT 2.0: Bio-Kidney AI 2026", fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    fig.tight_layout()
    
    # Veredicto
    presion_optima = presiones_kpa[np.argmin(np.abs(estres_celular_pa - 100))]
    print(f"[*] Presión de trabajo recomendada: {presion_optima:.1f} kPa")
    print(f"[*] Estrés máximo detectado: {max(estres_celular_pa):.2f} Pa")
    
    if max(estres_celular_pa) < 150:
        print("[VERDICTO]: FLUJO LAMINAR SEGURO. Viabilidad celular garantizada.")
    else:
        print("[ALERTA]: EXCESO DE PRESIÓN. Riesgo de lisis celular por cizallamiento.")
        
    plt.show()

if __name__ == "__main__":
    simular_bioimpresion_swift()
