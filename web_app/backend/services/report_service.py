"""
Servicio de Reportes BioKidney-AI (Pantalla 5 — Centro de Descarga)
-------------------------------------------------------------------
Genera los artefactos descargables tras un procesamiento exitoso:
  - Reporte analítico PDF (reportlab)
  - Malla estructural .OBJ y .STL derivada del árbol vascular fractal
"""

import json
import math
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, List, Tuple

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
)

from web_app.backend.services import storage_service
from web_app.backend.utils.logger import bk_logger

# Paleta de marca (Spec Maestra)
BRAND_BLUE = colors.HexColor("#0B2545")
BRAND_GREEN = colors.HexColor("#2E7D32")

FRACTAL_JSON = Path("renal_data_v8_fractal.json")


# ──────────────────────────────────────────────────────────────
# PDF
# ──────────────────────────────────────────────────────────────
def generate_pdf(project: Dict[str, Any], results: Dict[str, Any]) -> str:
    """Genera el Reporte Completo PDF y devuelve la ruta local."""
    path = storage_service.report_path(project["id"], "pdf")
    doc = SimpleDocTemplate(
        str(path), pagesize=A4,
        topMargin=2 * cm, bottomMargin=2 * cm,
        leftMargin=2 * cm, rightMargin=2 * cm,
    )

    styles = getSampleStyleSheet()
    h1 = ParagraphStyle("h1", parent=styles["Title"], textColor=BRAND_BLUE, fontSize=20)
    h2 = ParagraphStyle("h2", parent=styles["Heading2"], textColor=BRAND_GREEN)
    body = styles["BodyText"]

    story: List = []
    story.append(Paragraph("Bio-Kidney AI — Reporte de Simulación In-Sílico", h1))
    story.append(Spacer(1, 6))
    story.append(Paragraph(
        f"Validación hemodinámica renal · Generado {datetime.utcnow():%Y-%m-%d %H:%M UTC}",
        body))
    story.append(Spacer(1, 16))

    # Ficha del proyecto
    story.append(Paragraph("1. Ficha del Proyecto", h2))
    meta = [
        ["ID de proyecto", project["id"]],
        ["Nombre", project["project_name"]],
        ["Institución", project["institution"]],
        ["Email institucional", project["user_institution_email"]],
        ["TFG objetivo (mL/min)", str(project["gfr_constant"])],
        ["Dataset", project.get("dataset_filename") or "—"],
        ["Créditos de cómputo", str(project["simulation_credits_cost"])],
        ["Facturado (USD)", f"${project['total_billed_usd']:.2f}"],
    ]
    story.append(_kv_table(meta))
    story.append(Spacer(1, 16))

    # Resultados clave del pipeline
    story.append(Paragraph("2. Resultados del Pipeline Multiescala", h2))
    summary = results.get("summary", {})
    ms = results.get("modules_status", {})
    res_rows = [
        ["Estado global", results.get("global_status", "—")],
        ["TFG (mL/min)", str(summary.get("tfg", "—"))],
        ["Volumen orina (mL/min)", str(summary.get("urine", "—"))],
        ["Reabsorción (%)", str(summary.get("reabs_pct", "—"))],
        ["Hipoxia (%)", str(summary.get("hypoxia_pct", "—"))],
        ["Viabilidad biofabricación (%)", str(summary.get("viability_pct", "—"))],
        ["Cumplimiento Ley de Murray (%)", str(summary.get("murray_pct", "—"))],
    ]
    story.append(_kv_table(res_rows))
    story.append(Spacer(1, 16))

    # Estado por módulo
    story.append(Paragraph("3. Estado por Módulo", h2))
    mod_rows = [["Módulo", "Estado"]] + [[k.capitalize(), v] for k, v in ms.items()]
    t = Table(mod_rows, colWidths=[8 * cm, 6 * cm])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), BRAND_BLUE),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F2F6FA")]),
        ("PADDING", (0, 0), (-1, -1), 6),
    ]))
    story.append(t)
    story.append(Spacer(1, 20))

    story.append(Paragraph(
        "<i>Documento generado automáticamente por la plataforma Bio-Kidney AI. "
        "Uso científico e investigador. Enfocado en la dignidad humana.</i>",
        ParagraphStyle("foot", parent=body, textColor=colors.grey, fontSize=8)))

    doc.build(story)
    bk_logger.success(f"PDF generado: {path}")
    return str(path)


def _kv_table(rows: List[List[str]]) -> Table:
    t = Table(rows, colWidths=[7 * cm, 9 * cm])
    t.setStyle(TableStyle([
        ("FONTNAME", (0, 0), (0, -1), "Helvetica-Bold"),
        ("TEXTCOLOR", (0, 0), (0, -1), BRAND_BLUE),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#D0D7DE")),
        ("PADDING", (0, 0), (-1, -1), 5),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
    ]))
    return t


# ──────────────────────────────────────────────────────────────
# MALLA 3D (.OBJ / .STL) desde el árbol fractal
# ──────────────────────────────────────────────────────────────
def _load_segments() -> List[Dict[str, Any]]:
    """
    Carga y NORMALIZA los segmentos del árbol fractal a un esquema único:
      {id, nivel, x1_mm, y1_mm, z1_mm, x2_mm, y2_mm, z2_mm, radio_um}

    El JSON v8 mezcla dos formatos: uno explícito (x1_mm…, radio_um en µm) y
    otro compacto (p0=[x,y,z], p1=[x,y,z], r en mm, lvl).
    """
    if not FRACTAL_JSON.exists():
        return []
    with open(FRACTAL_JSON, "r", encoding="utf-8") as f:
        raw = json.load(f).get("segments", [])

    out: List[Dict[str, Any]] = []
    for i, s in enumerate(raw):
        if "x1_mm" in s:  # formato explícito
            out.append({
                "id": s.get("id", i), "nivel": s.get("nivel", 0),
                "x1_mm": s["x1_mm"], "y1_mm": s["y1_mm"], "z1_mm": s["z1_mm"],
                "x2_mm": s["x2_mm"], "y2_mm": s["y2_mm"], "z2_mm": s["z2_mm"],
                "radio_um": s.get("radio_um", 200.0),
            })
        elif "p0" in s and "p1" in s:  # formato compacto (r en mm)
            p0, p1 = s["p0"], s["p1"]
            out.append({
                "id": s.get("id", i), "nivel": s.get("lvl", 0),
                "x1_mm": p0[0], "y1_mm": p0[1], "z1_mm": p0[2],
                "x2_mm": p1[0], "y2_mm": p1[1], "z2_mm": p1[2],
                "radio_um": s.get("r", 0.2) * 1000.0,
            })
    return out


def _cylinder(p1: Tuple[float, float, float], p2: Tuple[float, float, float],
              radius: float, sides: int = 6):
    """Devuelve (vertices, triangles) de un cilindro low-poly entre p1 y p2."""
    ax, ay, az = p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]
    length = math.sqrt(ax * ax + ay * ay + az * az)
    if length < 1e-9:
        return [], []
    axis = (ax / length, ay / length, az / length)
    # vector perpendicular
    ref = (0.0, 0.0, 1.0) if abs(axis[2]) < 0.9 else (1.0, 0.0, 0.0)
    u = _cross(axis, ref)
    u = _norm(u)
    v = _norm(_cross(axis, u))

    verts = []
    for ring, center in ((0, p1), (1, p2)):
        for s in range(sides):
            ang = 2 * math.pi * s / sides
            offx = radius * (math.cos(ang) * u[0] + math.sin(ang) * v[0])
            offy = radius * (math.cos(ang) * u[1] + math.sin(ang) * v[1])
            offz = radius * (math.cos(ang) * u[2] + math.sin(ang) * v[2])
            verts.append((center[0] + offx, center[1] + offy, center[2] + offz))

    tris = []
    for s in range(sides):
        a = s
        b = (s + 1) % sides
        c = sides + s
        d = sides + (s + 1) % sides
        tris.append((a, c, b))
        tris.append((b, c, d))
    return verts, tris


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _norm(a):
    m = math.sqrt(a[0] ** 2 + a[1] ** 2 + a[2] ** 2) or 1.0
    return (a[0] / m, a[1] / m, a[2] / m)


def _build_mesh():
    """Ensambla vértices/triángulos globales del árbol vascular."""
    all_verts: List[Tuple[float, float, float]] = []
    all_tris: List[Tuple[int, int, int]] = []
    for seg in _load_segments():
        p1 = (seg["x1_mm"], seg["y1_mm"], seg["z1_mm"])
        p2 = (seg["x2_mm"], seg["y2_mm"], seg["z2_mm"])
        r = max(seg.get("radio_um", 200) / 1000.0, 0.15)  # µm -> mm, mínimo visible
        verts, tris = _cylinder(p1, p2, r)
        base = len(all_verts)
        all_verts.extend(verts)
        all_tris.extend((a + base, b + base, c + base) for a, b, c in tris)
    return all_verts, all_tris


def generate_obj(project_id: str) -> str:
    verts, tris = _build_mesh()
    path = storage_service.report_path(project_id, "obj")
    lines = ["# Bio-Kidney AI — arbol vascular renal (mm)", "o vascular_tree"]
    for x, y, z in verts:
        lines.append(f"v {x:.4f} {y:.4f} {z:.4f}")
    for a, b, c in tris:
        lines.append(f"f {a+1} {b+1} {c+1}")  # OBJ es 1-indexado
    path.write_text("\n".join(lines), encoding="utf-8")
    bk_logger.success(f"OBJ generado: {path} ({len(verts)} v, {len(tris)} f)")
    return str(path)


def generate_stl(project_id: str) -> str:
    verts, tris = _build_mesh()
    path = storage_service.report_path(project_id, "stl")
    out = ["solid vascular_tree"]
    for a, b, c in tris:
        v1, v2, v3 = verts[a], verts[b], verts[c]
        n = _norm(_cross(
            (v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]),
            (v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2]),
        ))
        out.append(f"  facet normal {n[0]:.5f} {n[1]:.5f} {n[2]:.5f}")
        out.append("    outer loop")
        for vx in (v1, v2, v3):
            out.append(f"      vertex {vx[0]:.4f} {vx[1]:.4f} {vx[2]:.4f}")
        out.append("    endloop")
        out.append("  endfacet")
    out.append("endsolid vascular_tree")
    path.write_text("\n".join(out), encoding="utf-8")
    bk_logger.success(f"STL generado: {path} ({len(tris)} triángulos)")
    return str(path)


def preview_segments() -> List[Dict[str, Any]]:
    """
    Segmentos low-res para el preview Three.js (Pantalla 3).
    Se entregan tal cual (ya son de baja resolución en el JSON fractal).
    """
    return _load_segments()
