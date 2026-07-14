"""
Router de Reportes — Pantalla 5 (Centro de Descarga)
Sirve el PDF analítico y los archivos estructurales .obj / .stl.
El acceso exige que el proyecto esté COMPLETED (gating de descarga).
"""

from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse

from web_app.backend.database.models import SessionLocal, Project
from web_app.backend.services import storage_service

router = APIRouter(prefix="/api/reports", tags=["reports"])

_MEDIA = {
    "pdf": "application/pdf",
    "obj": "text/plain",
    "stl": "application/vnd.ms-pki.stl",
}


def _guard_completed(project_id: str) -> Project:
    db = SessionLocal()
    try:
        proj = db.query(Project).get(project_id)
        if not proj:
            raise HTTPException(status_code=404, detail="Proyecto no encontrado")
        if proj.status != "COMPLETED":
            raise HTTPException(status_code=409, detail="La simulación aún no ha finalizado")
        return proj
    finally:
        db.close()


@router.get("/{project_id}/report.{ext}")
async def download(project_id: str, ext: str):
    if ext not in _MEDIA:
        raise HTTPException(status_code=404, detail="Formato no disponible")
    _guard_completed(project_id)
    path = storage_service.report_path(project_id, ext)
    if not path.exists():
        raise HTTPException(status_code=404, detail=f"Archivo .{ext} no generado")
    return FileResponse(
        str(path), media_type=_MEDIA[ext],
        filename=f"biokidney_{project_id}.{ext}",
    )
