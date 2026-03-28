# 🧬 Guía de Despliegue y Uso: BioKidney-AI 2026

Esta guía detalla los pasos necesarios para configurar, desplegar y operar la plataforma de simulación de bioingeniería renal **BioKidney-AI**.

---

## 1. Requisitos Previos

Antes de comenzar, asegúrate de tener instalado:

*   **Python 3.10+**
*   **uv**: Herramienta de gestión de paquetes ultra rápida.
    *   *Instalación:* `curl -LsSf https://astral.sh/uv/install.sh | sh`
*   **Docker & Docker Compose**: Para el despliegue de contenedores.
*   **Blender 4.0+** (Opcional): Solo si deseas visualizar o renderizar los modelos 3D generados (`.blend`).

---

## 2. Configuración del Entorno de Desarrollo

La plataforma utiliza `uv` para una gestión eficiente de dependencias.

### Paso 1: Clonar el repositorio
```bash
git clone <url-del-repositorio>
cd BioKidney-AI
```

### Paso 2: Crear el Entorno Virtual e Instalar Dependencias
```bash
# Crear el entorno virtual
uv venv

# Activar el entorno
source .venv/bin/activate  # Linux/macOS
# .venv\Scripts\activate   # Windows

# Instalar dependencias desde el archivo de requerimientos
uv pip install -r requirements.txt
```

---

## 3. Ejecución de Simuladores (CLI)

Los simuladores biofísicos se encuentran en la carpeta `01_simuladores/`. Puedes ejecutarlos individualmente para generar datos y reportes.

### Ejemplo: Simulador de Filtración Glomerular
```bash
python 01_simuladores/simulador_filtracion_glomerular.py
```
*   **Resultado:** Generará archivos `.png` con gráficas de presión y un reporte PDF en la misma carpeta.

### Ejemplo: Generador de Árbol Vascular (CCO)
```bash
python 02_vascular_cco/generador_cco_v7.py
```
*   **Resultado:** Crea archivos `.csv` y `.obj` que representan la estructura vascular del riñón bioartificial.

---

## 4. Despliegue de la Aplicación Web

La plataforma incluye un Dashboard Maestro para visualizar las simulaciones de forma integrada.

### Opción A: Ejecución Local (Desarrollo)
1.  **Iniciar el Backend (FastAPI):**
    ```bash
    python web_app/backend/main.py
    ```
2.  **Abrir el Frontend:**
    Abre el archivo `web_app/frontend/dashboard.html` en tu navegador.

### Opción B: Despliegue con Docker (Producción/Testing)
```bash
docker-compose up --build
```
*   La aplicación estará disponible en `http://localhost:8000`.

---

## 5. Integración con Blender (Modelado 3D)

Para convertir los datos vasculares en modelos 3D volumétricos:
1.  Abre Blender.
2.  Ve a la pestaña `Scripting`.
3.  Carga y ejecuta el script `02_vascular_cco/exportar_blender_v2.py`.
4.  Importa el archivo `.csv` generado por los simuladores CCO.

---

## 6. Estructura de Resultados

*   **`/resultados`**: Gráficas de optimización (Pareto).
*   **`/01_simuladores`**: Informes técnicos y capturas de dashboards.
*   **`/07_presentacion_final`**: Scorecards y resúmenes ejecutivos para inversores o revisores científicos.

---

## 7. Solución de Problemas

*   **Error de dependencias**: Ejecuta `uv pip install -r requirements.txt --upgrade` para asegurar que tienes las versiones correctas.
*   **Falta de permisos en Linux**: Si usas `inicio.sh`, asegúrate de darle permisos de ejecución: `chmod +x inicio.sh`.
