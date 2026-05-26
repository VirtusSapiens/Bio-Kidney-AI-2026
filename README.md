# Bio-Kidney AI 2026

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19508077.svg)](https://doi.org/10.5281/zenodo.19508077)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/release/python-3120/)
[![Platform](https://img.shields.io/badge/platform-Ubuntu%2024.04-orange.svg)](https://ubuntu.com/)

**Multi-Scale Biophysical Simulation of Functional Viability in an iPSC-Derived Bioprinted Human Kidney**

> *Six-Module Integration Predicts Physiological Renal Output Across Vascular, Cellular, and Functional Scales*

**Carlos David Moreno Cáceres** — Independent Researcher, VirtusSapiens, Medellín, Colombia
ORCID: [0009-0005-3933-5072](https://orcid.org/0009-0005-3933-5072)

---

## Overview

Bio-Kidney AI 2026 is a multi-scale computational framework for the complete *in-silico* validation of a functional bioprinted human kidney. The framework integrates six interdependent biophysical modules into a unified predictive pipeline, enabling pre-fabrication viability assessment of bioengineered organ constructs.

The complete pipeline executes in **under 60 seconds** on consumer hardware (Intel Core i5, 16 GB RAM, Ubuntu Linux) and is implemented as an open-source Python package with a **Mixture of Experts (MoE)** architecture.

### Key Results (Phase 1)

| Module | Output | Physiological Range |
|--------|--------|-------------------|
| Vascular CCO v8 | 1,902 segments, 100% Murray's Law compliance | ✅ |
| Hemodynamic pressure | 58.6 ± 13.4 mmHg terminal capillary pressure | ✅ (threshold: 43 mmHg) |
| Glomerular filtration | **115.2 mL/min bilateral GFR** | ✅ (normal: 100–125 mL/min) |
| Oxygen diffusion | 0 hypoxic voxels, min PO₂ = 5.6 mmHg | ✅ (threshold: 4 mmHg) |
| iPSC differentiation | 100% phenotypic purity by day 21 | ✅ |
| Co-SWIFT bioprinting | 98% post-print cell viability | ✅ |
| Tubular reabsorption | 98.1% reabsorption efficiency, 2.19 L/day urine | ✅ |

---

## Architecture

The system is organized into scientific domain layers:

```
BioKidney-AI/
├── biokidney/                  # Core MoE package
│   ├── core/                   # Pipeline orchestration (config, I/O)
│   ├── experts/                # Domain modules (vascular, cellular, fluids)
│   └── simulation/             # Fractal vascularizer, v8 converter
├── 01_simuladores/             # Standalone simulation modules
├── 02_vascular_cco/            # CCO v1→v8 generators, data exports
├── 03_modelos_3d/              # Blender scenes, STL files
├── 04_literatura/              # Reference papers
├── 06_app/                     # Desktop GUI applications (PyQt6)
├── 07_presentacion_final/      # Academic presentations, reports
└── web_app/                    # FastAPI backend + dashboard frontend
```

### Module 1 — Vascular Tree Generation (CCO v8)
Constrained Constructive Optimization algorithm generating hierarchical arterial, venous, and collecting trees. Key innovations in v8:
- **Beta(3, 1.2) cortical demand distribution** — biases 63.1% of glomerular demand points toward the cortex, replicating human nephron anatomy
- **Two-pass calibrated Poiseuille model** — resolves pressure computation failure in geometrically simplified trees; delivers physiological terminal pressures of 58.6 ± 13.4 mmHg
- **cKDTree spatial indexing** — reduces coverage queries from O(N·M) to O(N log N); ~100× speedup over v7
- Full Murray's Law compliance (α = 3.0) across all 915 bifurcations

### Module 2 — Oxygen Diffusion (Krogh Cylinder)
Steady-state diffusion solved via Successive Over-Relaxation (SOR) on a 30×30 voxel grid. Zero hypoxic zones confirmed across the entire tissue volume.

### Module 3 — iPSC Differentiation Kinetics
First-order kinetic model following the Takasato 2015 protocol. Three renal lineages (podocytes, proximal tubule cells, loop of Henle cells) reach >95% phenotypic purity by day 15. Residual OCT4 < 0.1% at day 21 (low teratoma risk).

### Module 4 — Co-SWIFT Bioprinting Optimization
Pareto multi-objective optimization of bioprinting parameters using Herschel-Bulkley rheology, calibrated for NICE bioink (GelMA 7% + Alginate 1.5% + Nanocellulose 0.8% + LAP 0.25%). 100 Pareto-optimal solutions identified.

### Module 5 — Glomerular Filtration (Starling-Deen)
Starling equation integrated along the glomerular capillary, coupled with the Deen-Robertson-Brenner oncotic pressure model. In integrated pipeline mode, receives calibrated pressure field from CCO v8, yielding bilateral GFR = 115.2 mL/min.

### Module 6 — Tubular Reabsorption (Michaelis-Menten / Kedem-Katchalsky)
Five-segment nephron model covering proximal tubule, descending and ascending loop of Henle, distal tubule, and collecting duct. Key transporters: SGLT2, NHE3, AQP1, NKCC2, ENaC, AQP2.

---

## Installation

```bash
# Clone the repository
git clone https://github.com/VirtusSapiens/Bio-Kidney-AI-2026.git
cd Bio-Kidney-AI-2026

# Create and activate virtual environment
python3 -m venv env_biokidney
source env_biokidney/bin/activate

# Install dependencies
pip install -r requirements.txt
```

**Scientific computing environment:**
- Python 3.12
- NumPy 2.4.3 · SciPy 1.17.1 · Matplotlib 3.10.8

**Web platform environment (isolated):**
- FastAPI 0.104.1 · SQLAlchemy 2.0.23 · Loguru 0.7.2 · Docker

---

## Usage

### Run the complete pipeline
```bash
bash run_pipeline.sh
```

### Run individual simulation modules
```bash
source env_biokidney/bin/activate

# Vascular tree generation (CCO v8)
python3 02_vascular_cco/generador_cco_v8.py

# Glomerular filtration simulator
python3 01_simuladores/simulador_filtracion_glomerular.py

# iPSC differentiation kinetics
python3 01_simuladores/simulador_diferenciacion_ipsc.py

# Oxygen diffusion
python3 01_simuladores/simulador_oxigeno_biokidney.py

# Tubular reabsorption
python3 01_simuladores/simulador_reabsorcion_tubular.py

# Co-SWIFT bioprinting optimizer
python3 01_simuladores/optimizador_coswift.py
```

### Launch web dashboard
```bash
bash inicio.sh
# Dashboard: http://localhost:8000/
# API: http://localhost:8000/api/simulate
```

### Docker deployment
```bash
docker-compose build
docker-compose up
```

---

## Utility Scripts

| Script | Purpose |
|--------|---------|
| `inicio.sh` | Single entry point — configures environment, shows system status |
| `analizador_proyecto_biokidney.py` | Recursive project scan → generates `CONTEXTO_PROYECTO.md` for AI agents |
| `biokidney_architect.py` | Strategic planner → generates `BLUEPRINT_INGENIERIA.md` |
| `blender_vascular_tree.py` | Blender integration for 3D vascular tree rendering |

---

## Preprint

The Phase 1 preprint is available on Zenodo:

> Moreno Cáceres, C.D. (2026). *Multi-Scale Biophysical Simulation of Functional Viability in an iPSC-Derived Bioprinted Human Kidney: Six-Module Integration Predicts Physiological Renal Output Across Vascular, Cellular, and Functional Scales*. Zenodo. https://doi.org/10.5281/zenodo.19508077

---

## Roadmap

### Phase 1 — Complete ✅
- CCO v8 vascular tree with calibrated hemodynamics
- Six-module integrated pipeline
- Bilateral GFR = 115.2 mL/min (normal physiological range)
- Open-source Python package with MoE architecture
- Preprint published (Zenodo, May 2026 v2)

### Phase 2 — Planned
- [ ] Single-nephron hemodynamics resolution
- [ ] Tubuloglomerular feedback modeling
- [ ] Myogenic autoregulation and renin-angiotensin-aldosterone axis
- [ ] Full CFD analysis on generated vascular tree (OpenFOAM)
- [ ] Patient-specific parameterization
- [ ] CNN-based vascular tree optimization
- [ ] *In-vitro* validation protocol

---

## Citation

If you use this framework in your research, please cite:

```bibtex
@misc{morenoCaceres2026biokidney,
  author       = {Moreno Cáceres, Carlos David},
  title        = {Multi-Scale Biophysical Simulation of Functional Viability
                  in an iPSC-Derived Bioprinted Human Kidney},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.19508077},
  url          = {https://doi.org/10.5281/zenodo.19508077}
}
```

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

## Acknowledgements

The author thanks John Tapias (VECANOVA) for software engineering contributions.

---

*Independent research — VirtusSapiens, Medellín, Colombia, 2026*
