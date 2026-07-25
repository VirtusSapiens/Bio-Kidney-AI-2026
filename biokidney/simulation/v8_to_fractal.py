import json
import numpy as np
import os
from biokidney.simulation.fractal_vascularizer import FractalBrancher

def integrate_fractal_layer(input_json, output_json):
    if not os.path.exists(input_json):
        print(f'Error: No se encuentra {input_json}')
        return
    with open(input_json, 'r') as f:
        data = json.load(f)
    segments = data.get('segments', [])
    
    # Filtro de terminales reales
    terminal_nodes = [s for s in segments if s.get('nivel', 0) >= 10 or s.get('es_terminal') == True]

    brancher = FractalBrancher(murray_alpha=3.0, min_radius_um=8.0)
    all_fractal_segments = []
    print(f'Procesando {len(terminal_nodes)} terminales con flujo orientado...')
    
    for term in terminal_nodes:
        # Coordenadas actuales
        p1 = np.array([term['x2_mm'], term['y2_mm'], term['z2_mm']])
        p0 = np.array([term['x1_mm'], term['y1_mm'], term['z1_mm']])
        
        # DIRECCIÓN: Calculamos el vector de salida para que el fractal siga la inercia
        v_dir = p1 - p0
        r_mm = term['radio_um'] / 1000.0
        
        brancher.segments = []
        # Crecimiento fractal orientado (Inercia centrífuga)
        brancher.grow_tree(p1, v_dir, r_mm, 0, 3, np.linalg.norm(v_dir)*0.6)
        all_fractal_segments.extend(brancher.segments)
        
    extended_data = {
        'metadata': {'version': 'BioKidney v8.2 Orientada', 'hilio_lateral': True},
        'segments': segments + all_fractal_segments
    }
    
    with open(output_json, 'w') as f:
        json.dump(extended_data, f, indent=4)
    print(f'¡Éxito! Red coherente generada en: {output_json}')

if __name__ == '__main__':
    integrate_fractal_layer('02_vascular_cco/renal_data_v1.json', 'renal_data_v8_fractal.json')
