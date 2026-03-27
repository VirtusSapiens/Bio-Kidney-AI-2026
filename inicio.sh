#!/bin/bash
# ------------------------------------------------------------------------------
# Script de Inicio: Bio-Kidney AI 2026
# Propósito: Configurar el entorno de desarrollo y proporcionar accesos rápidos.
# ------------------------------------------------------------------------------

# Configuración de Rutas
PROYECTO_DIR="$HOME/Escritorio/BioKidney-AI"
VENV_DIR="$HOME/env_biokidney"

# Colores para la interfaz
CYAN='\033[0;36m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

clear
echo -e "${CYAN}==================================================${NC}"
echo -e "${CYAN}       BIO-KIDNEY AI 2026 - VirtusSapiens         ${NC}"
echo -e "${CYAN}==================================================${NC}"

# 1. Verificación de Entorno
if [ -d "$VENV_DIR" ]; then
    source "$VENV_DIR/bin/activate"
    echo -e "${GREEN}[OK]${NC} Entorno Virtual: env_biokidney activado."
else
    echo -e "${YELLOW}[!]${NC} Advertencia: No se encontró el entorno en $VENV_DIR"
fi

# 2. Navegación
if [ -d "$PROYECTO_DIR" ]; then
    cd "$PROYECTO_DIR"
    echo -e "${GREEN}[OK]${NC} Directorio de trabajo: $(pwd)"
else
    echo -e "${YELLOW}[!]${NC} Error: No se pudo encontrar la carpeta del proyecto."
fi

echo ""
echo -e "${CYAN}-- ESTADO DEL SISTEMA --${NC}"
# Mostrar hito actual desde el analizador si existe
if [ -f "CONTEXTO_PROYECTO.md" ]; then
    grep "Estado del Sistema" CONTEXTO_PROYECTO.md | cut -d':' -f2
else
    echo "Analizando sistema por primera vez..."
    python3 analizador_proyecto_biokidney.py > /dev/null
fi

echo ""
echo -e "${CYAN}-- ARCHIVOS RECIENTES (Python) --${NC}"
find . -maxdepth 3 -name "*.py" -not -path '*/.*' -printf '%T@ %p\n' | sort -nr | head -5 | cut -d' ' -f2

echo ""
echo -e "${CYAN}-- COMANDOS RÁPIDOS --${NC}"
echo -e " ${GREEN}Dashboard   :${NC} python3 06_app/dashboard_maestro_app.py"
echo -e " ${GREEN}Filtración  :${NC} python3 06_app/filtracion_glomerular_gui.py"
echo -e " ${GREEN}Analizador  :${NC} python3 analizador_proyecto_biokidney.py"
echo -e " ${GREEN}Arquitecto  :${NC} python3 biokidney_architect.py"
echo ""
echo -e "${CYAN}==================================================${NC}"
