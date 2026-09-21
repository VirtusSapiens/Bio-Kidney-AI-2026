#!/bin/bash
set -e
source /home/virtussapiens/Escritorio/BioKidney-AI/.venv/bin/activate
cd /home/virtussapiens/Escritorio/BioKidney-AI
echo "=== [1/4] Generando red vascular CCO v7 ==="
python 02_vascular_cco/generador_cco_v7.py
echo "=== [2/4] Simulando filtración glomerular ==="
python 01_simuladores/simulador_filtracion_glomerular.py
echo "=== [3/4] Simulando difusión O2 ==="
python 01_simuladores/simulador_oxigeno_biokidney.py
echo "=== [4/4] Simulando reabsorción tubular ==="
python 01_simuladores/simulador_reabsorcion_tubular.py
echo "=== PIPELINE COMPLETO ==="
