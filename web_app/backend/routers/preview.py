"""
Router de Preview — Pantalla 3 (Preview de Malla 3D)
Entrega los segmentos vasculares low-res para Three.js y el cálculo de
recursos requeridos.
"""

from fastapi import APIRouter, HTTPException

from web_app.backend.database.models import SessionLocal, Project
from web_app.backend.services import report_service, pricing_service

router = APIRouter(prefix="/api/projects", tags=["preview"])


@router.get("/{project_id}/preview")
async def get_preview(project_id: str):
    db = SessionLocal()
    try:
        proj = db.query(Project).get(project_id)
        if not proj:
            raise HTTPException(status_code=404, detail="Proyecto no encontrado")

        quote = pricing_service.quote(proj.dataset_size_bytes, proj.gfr_constant)
        if proj.status == "DRAFT":
            proj.status = "PREVIEWED"
            db.commit()

        return {
            "project_id": project_id,
            "segments": report_service.preview_segments(),
            "quote": quote,
        }
    finally:
        db.close()
