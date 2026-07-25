"""
Pipeline I/O — carga renal_data_v1.json generado por CCO v8.
Compatible hacia atrás con el formato plano de v7.
"""

import json
import numpy as np
import os

PIPELINE_PATH = os.path.expanduser(
    "~/Escritorio/BioKidney-AI/02_vascular_cco/renal_data_v1.json"
)


def load_cco_data(path=PIPELINE_PATH):
    """
    Carga datos del arbol vascular CCO.
    Soporta:
      - renal_data_v1.json (v8, esquema completo)
      - pipeline_data.json (v7, esquema plano)
    """
    with open(path, "r") as f:
        raw = json.load(f)

    # v8 format (tiene "segments" y "terminals")
    if "segments" in raw:
        segments = raw["segments"]
        terminals = raw.get("terminals", [])
        return {
            "version"          : raw.get("version", "unknown"),
            "n_segments"       : raw["n_segments"],
            "radii_um"         : np.array(raw["radii_um"], dtype=np.float32),
            "pressures_kpa"    : np.array(raw["pressures_kpa"], dtype=np.float32),
            "coverage_pct"     : raw.get("statistics", {}).get("coverage_pct", 0),
            "murray_compliance": raw.get("statistics", {}).get("murray_compliance", 0),
            "n_terminals"      : len(terminals),
            "terminal_pressures_mmHg": np.array(
                [t["presion_mmHg"] for t in terminals], dtype=np.float32
            ) if terminals else np.array([], dtype=np.float32),
            "terminal_positions_mm": np.array(
                [[t["x_mm"], t["y_mm"], t["z_mm"]] for t in terminals],
                dtype=np.float32
            ) if terminals else np.empty((0, 3), dtype=np.float32),
        }

    # v7 fallback (esquema plano)
    return {
        "version"          : "v7",
        "n_segments"       : raw["n_segments"],
        "radii_um"         : np.array(raw["radii_um"], dtype=np.float32),
        "pressures_kpa"    : np.array(raw["pressures_kpa"], dtype=np.float32),
        "coverage_pct"     : raw.get("coverage_pct", 0),
        "murray_compliance": raw.get("murray_compliance", 0),
        "n_terminals"      : 0,
        "terminal_pressures_mmHg": np.array([], dtype=np.float32),
        "terminal_positions_mm"  : np.empty((0, 3), dtype=np.float32),
    }
