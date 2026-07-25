#!/usr/bin/env python3
"""Verificacion cruzada Tabla A1 + parrafo morfologia (3.1) vs renal_data_v1.json.
Escribe el resultado a /tmp/verify_result.txt (no a stdout)."""
import json, re, os

ROOT = "/home/virtussapiens/Escritorio/BioKidney-AI"
JSON_PATH = os.path.join(ROOT, "02_vascular_cco/renal_data_v1.json")
MD_PATH = os.path.join(ROOT, "00_bitacora/preprint_biokidney_2026_EN_v4.md")
OUT = "/tmp/verify_result.txt"
TOL = 0.1

# ---- 1. recomputar desde JSON (misma logica que audit_radii.py) ----
d = json.load(open(JSON_PATH))
art = [s for s in d["segments"] if s["sistema"] == "art"]
ven = [s for s in d["segments"] if s["sistema"] == "ven"]
by = {}
for s in art:
    by.setdefault(s["nivel"], []).append(s["radio_um"])
json_lvl = {lv: (min(rs), max(rs)) for lv, rs in by.items()}  # lv -> (rmin, rmax)

all_art_r = [s["radio_um"] for s in art]
json_mean = sum(all_art_r) / len(all_art_r)
json_art_min = min(all_art_r)
json_ven_max = max(s["radio_um"] for s in ven)

# ---- 2. parsear la tabla A1 (Appendix A, tabla de 6 columnas) ----
lines = open(MD_PATH, encoding="utf-8").read().splitlines()
ap_start = next(i for i, l in enumerate(lines) if l.startswith("## Appendix A"))
table_rows = {}  # lv -> (rmin, rmax) segun el archivo
for l in lines[ap_start:]:
    m = re.match(r"\|\s*(\d+)\s*\|[^|]*\|\s*(\d+)\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|", l)
    if m:
        lv = int(m.group(1))
        table_rows[lv] = (float(m.group(3)), float(m.group(4)))

# ---- 3. comparar celda por celda (niveles 1..11) ----
out = []
out.append("=== TABLA A1: R_min / R_max  publicado vs JSON (tol %.1f um) ===" % TOL)
out.append("LVL | tabla_Rmin | json_Rmin | MATCH | tabla_Rmax | json_Rmax | MATCH")
ok = 0
for lv in range(1, 12):
    trmin, trmax = table_rows.get(lv, (None, None))
    jrmin, jrmax = json_lvl.get(lv, (None, None))
    mmin = (trmin is not None and abs(trmin - jrmin) <= TOL)
    mmax = (trmax is not None and abs(trmax - jrmax) <= TOL)
    if mmin and mmax:
        ok += 1
    out.append("%3d | %10s | %9.1f | %5s | %10s | %9.1f | %5s" % (
        lv, trmin, jrmin, "MATCH" if mmin else "MISMATCH",
        trmax, jrmax, "MATCH" if mmax else "MISMATCH"))
out.append("RESULTADO: %d/11 filas coinciden (R_min y R_max ambos)" % ok)

# ---- 4. parrafo morfologia 3.1: extraer los numeros um y comparar ----
para = next(l for l in lines if l.startswith("**Vascular morphology.**"))
nums = [float(x) for x in re.findall(r"\d+\.\d+", para)]  # solo decimales (um)
expected = {
    "root_seg_max(L1)": round(json_lvl[1][1], 1),
    "arcuate_min(L4)": round(json_lvl[4][0], 1),
    "arcuate_max(L4)": round(json_lvl[4][1], 1),
    "interlob_min(L5)": round(json_lvl[5][0], 1),
    "interlob_max(L5)": round(json_lvl[5][1], 1),
    "finest_terminal": round(json_art_min, 1),
    "mean_art_radius": round(json_mean, 1),
    "venous_max": round(json_ven_max, 1),
}
out.append("")
out.append("=== 3.1 morfologia: cifra JSON presente en el parrafo? ===")
p_ok = 0
for name, val in expected.items():
    present = any(abs(val - n) <= TOL for n in nums)
    if present:
        p_ok += 1
    out.append("%-18s json=%8.1f | %s" % (name, val, "PRESENTE" if present else "AUSENTE"))
out.append("RESULTADO: %d/%d cifras del JSON presentes en 3.1" % (p_ok, len(expected)))
out.append("(numeros um hallados en el parrafo: %s)" % nums)

open(OUT, "w", encoding="utf-8").write("\n".join(out) + "\n")
