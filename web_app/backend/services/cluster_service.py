"""
Orquestador del Clúster BioKidney-AI
------------------------------------
Ejecuta el pipeline científico completo (el "clúster HPC") y genera los
artefactos descargables. IMPORTANTE (regla de negocio de la spec):

    Esta función SOLO debe invocarse desde el handler del webhook de pago,
    una vez confirmado Payment_Success. Nunca desde un endpoint accesible
    directamente por el frontend.
"""

import asyncio
from typing import Dict, Any

from web_app.backend.database.models import SessionLocal, Project
from web_app.backend.services.simulation_service import SimulationService
from web_app.backend.services import report_service, storage_service
from web_app.backend.utils.logger import bk_logger

_sim = SimulationService()


def _set(project_id: str, **fields) -> None:
    db = SessionLocal()
    try:
        proj = db.query(Project).get(project_id)
        if not proj:
            return
        for k, v in fields.items():
            setattr(proj, k, v)
        db.commit()
    finally:
        db.close()


async def process_project(project_id: str) -> None:
    """Corre el pipeline y produce PDF/OBJ/STL. Actualiza progreso y estado."""
    bk_logger.info(f"[CLÚSTER] Iniciando procesamiento de {project_id}")
    db = SessionLocal()
    try:
        proj = db.query(Project).get(project_id)
        if not proj:
            bk_logger.error(f"[CLÚSTER] Proyecto {project_id} no existe")
            return
        project_dict = proj.to_dict()
        params = {
            "vascular": {"n_seeds": max(500, proj.estimated_nodes // 2)},
            "filtration": {"n_glomerulos": 1_000_000},
            **(proj.sim_params or {}),
        }
    finally:
        db.close()

    _set(project_id, status="PROCESSING", progress=10)

    # 1. Pipeline científico completo
    results = await _sim.run_full_pipeline(params)
    _set(project_id, progress=70, results=results)

    # 2. Artefactos descargables
    await asyncio.to_thread(report_service.generate_pdf, project_dict, results)
    _set(project_id, progress=85)
    await asyncio.to_thread(report_service.generate_obj, project_id)
    await asyncio.to_thread(report_service.generate_stl, project_id)

    # 3. Finalizar — poblar Report_Generated_URL (crítico #5)
    _set(
        project_id,
        progress=100,
        status="COMPLETED",
        report_url=storage_service.public_url(project_id),
    )
    bk_logger.success(f"[CLÚSTER] Proyecto {project_id} COMPLETADO")
