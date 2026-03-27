# BioKidney-AI 2026

## 🚀 Propósito del Sistema
**BioKidney-AI** es un ecosistema de simulación multiescala diseñado para el desarrollo y validación *in-silico* de riñones artificiales bioimpresos funcionales. El sistema integra desde la generación sintética de árboles vasculares basados en leyes biológicas (Murray) hasta la modelación compleja de la filtración glomerular y reabsorción tubular, permitiendo predecir la viabilidad fisiológica de un órgano antes de su biofabricación.

---

## 🏛️ Arquitectura del Proyecto
El software está organizado en capas de dominio científico y técnico:

### 1. Capa Geométrica y Vascular (`02_vascular_cco`)
Contiene algoritmos de **Optimización Constructiva Restringida (CCO)** para generar redes vasculares (arteriales, venosas y colectoras) que cumplen con la Ley de Murray y aseguran una cobertura total del parénquima renal.

### 2. Capa de Fenomenología Física (`01_simuladores`)
Módulos que ejecutan simulaciones de transporte de masa:
- **Filtración Glomerular:** Ecuaciones de Starling + Modelo Deen.
- **Reabsorción Tubular:** Cinética de transportadores (Michaelis-Menten) y balance hídrico.
- **Difusión de Oxígeno:** Modelado 3D de gradientes de presión parcial de O2 (Fick).
- **Cinética Celular:** Maduración de iPSCs y protocolos de diferenciación.

### 3. Capa de Visualización y CAD (`03_modelos_3d`)
Integración con **Blender** mediante scripts de Python (`bpy`) para la generación de mallas 3D realistas, renderizado de alta fidelidad y preparación de archivos para bioimpresoras 3D.

### 4. Capa de Interfaz y Control (`06_app`)
Aplicaciones de escritorio desarrolladas en **PyQt6** que actúan como dashboards maestros para visualizar KPIs en tiempo real y realizar ajustes paramétricos.

---

## 🛠️ Herramientas de Gestión (Root Scripts)

En la raíz del proyecto se encuentran las herramientas de orquestación mejoradas bajo estándares de ingeniería de software:

### 1. `inicio.sh`
**Propósito:** Punto de entrada único para el desarrollador.
- Configura automáticamente el entorno virtual (`env_biokidney`).
- Verifica el estado actual del sistema.
- Muestra los últimos archivos modificados y comandos rápidos de ejecución.
- **Uso:** `bash inicio.sh`

### 2. `analizador_proyecto_biokidney.py`
**Propósito:** Generador de la Base de Conocimiento.
- Realiza un escaneo recursivo del proyecto extrayendo metadatos, dependencias y propósitos de cada módulo.
- Genera el archivo `CONTEXTO_PROYECTO.md` utilizado por agentes de IA y arquitectos de sistema.
- **Uso:** `python3 analizador_proyecto_biokidney.py`

### 3. `biokidney_architect.py`
**Propósito:** Planificador Estratégico.
- Analiza el `CONTEXTO_PROYECTO.md` para proponer mejoras arquitectónicas (patrones de diseño, validaciones de datos).
- Genera un `BLUEPRINT_INGENIERIA.md` para guiar el desarrollo de nuevas funcionalidades.
- **Uso:** `python3 biokidney_architect.py "Refactorización de simuladores"`

---

## 📊 Estado Actual: Pipeline 100% Completado
A fecha de Marzo 2026, el sistema ha superado con éxito todas las fases de validación *in-silico*, demostrando la capacidad teórica de un riñón bioimpreso para alcanzar un TFG de ~82 mL/min y una producción de orina homeostática de ~1.5 mL/min.

**Investigador Principal:** Carlos David Moreno Cáceres (VirtusSapiens)  
**Ubicación:** Medellín, Colombia
