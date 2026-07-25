"""
Router de Proyectos — Pantalla 2 (Configuración / Input)
Crea el proyecto y recibe el dataset biológico (Drag & Drop, ≤500MB).
"""

from fastapi import APIRouter, UploadFile, File, Form, HTTPException

from web_app.backend.database.models import SessionLocal, Project
from web_app.backend.services import storage_service, pricing_service
from web_app.backend.services.storage_service import UploadTooLarge
from web_app.backend.utils.logger import bk_logger

router = APIRouter(prefix="/api/projects", tags=["projects"])


@router.post("")
async def create_project(
    project_name: str = Form(...),
    institution: str = Form(...),
    user_institution_email: str = Form(...),
    gfr_constant: float = Form(62.5),
    dataset: UploadFile = File(...),
):
    """Crea un proyecto en estado DRAFT y guarda el dataset biológico."""
    db = SessionLocal()
    try:
        proj = Project(
            project_name=project_name,
            institution=institution,
            user_institution_email=user_institution_email,
            gfr_constant=gfr_constant,
            status="DRAFT",
        )
        db.add(proj)
        db.commit()
        db.refresh(proj)
        project_id = proj.id

        # Guardar dataset con límite de 500MB (regla de negocio)
        try:
            saved = storage_service.save_dataset(
                project_id, dataset.filename, dataset.file
            )
        except UploadTooLarge as e:
            db.delete(proj)
            db.commit()
            raise HTTPException(status_code=413, detail=str(e))

        # Cotización de recursos (para la Pantalla 3)
        quote = pricing_service.quote(saved["size_bytes"], gfr_constant)

        proj.biological_dataset_path = saved["path"]
        proj.dataset_filename = saved["filename"]
        proj.dataset_size_bytes = saved["size_bytes"]
        proj.simulation_credits_cost = quote["simulation_credits_cost"]
        proj.total_billed_usd = quote["total_billed_usd"]
        proj.estimated_nodes = quote["estimated_nodes"]
        proj.status = "PREVIEWED"
        db.commit()
        db.refresh(proj)

        bk_logger.info(f"Proyecto creado: {project_id} ({project_name})")
        return {"project": proj.to_dict(), "quote": quote}
    finally:
        db.close()


@router.get("/{project_id}")
async def get_project(project_id: str):
    db = SessionLocal()
    try:
        proj = db.query(Project).get(project_id)
        if not proj:
            raise HTTPException(status_code=404, detail="Proyecto no encontrado")
        return proj.to_dict()
    finally:
        db.close()


@router.get("/{project_id}/status")
async def project_status(project_id: str):
    """Polling de la barra de progreso del clúster (Pantalla 5)."""
    db = SessionLocal()
    try:
        proj = db.query(Project).get(project_id)
        if not proj:
            raise HTTPException(status_code=404, detail="Proyecto no encontrado")
        return {
            "id": proj.id,
            "status": proj.status,
            "progress": proj.progress,
            "report_url": proj.report_url,
        }
    finally:
        db.close()
