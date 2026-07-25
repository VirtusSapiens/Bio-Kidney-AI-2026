"""
Servicio de Almacenamiento BioKidney-AI
---------------------------------------
Abstrae el almacenamiento de datasets biológicos subidos y de los reportes
generados. En el MVP usa disco local; la interfaz está lista para migrar a
S3 (Report_Generated_URL) sin tocar los routers.
"""

import os
import shutil
from pathlib import Path
from typing import BinaryIO

from web_app.backend.utils.logger import bk_logger

# Regla de negocio de la spec: máximo 500MB por proyecto
MAX_UPLOAD_BYTES = 500 * 1024 * 1024

_BASE = Path(os.environ.get("BIOKIDNEY_STORAGE", "web_app/storage"))
UPLOAD_DIR = _BASE / "uploads"
REPORT_DIR = _BASE / "reports"


def _ensure_dirs() -> None:
    UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)


class UploadTooLarge(Exception):
    """Se supera el límite de 500MB permitido por proyecto."""


def save_dataset(project_id: str, filename: str, fileobj: BinaryIO) -> dict:
    """
    Guarda el dataset biológico en streaming, aplicando el límite de 500MB.
    Devuelve ruta y tamaño en bytes.
    """
    _ensure_dirs()
    safe_name = os.path.basename(filename or "dataset.bin")
    dest = UPLOAD_DIR / f"{project_id}__{safe_name}"

    written = 0
    with open(dest, "wb") as out:
        while True:
            chunk = fileobj.read(1024 * 1024)
            if not chunk:
                break
            written += len(chunk)
            if written > MAX_UPLOAD_BYTES:
                out.close()
                dest.unlink(missing_ok=True)
                raise UploadTooLarge(
                    f"El archivo supera el límite de {MAX_UPLOAD_BYTES // (1024*1024)}MB."
                )
            out.write(chunk)

    bk_logger.info(f"Dataset guardado: {dest} ({written} bytes)")
    return {"path": str(dest), "filename": safe_name, "size_bytes": written}


def report_path(project_id: str, ext: str) -> Path:
    """Ruta local de un artefacto de reporte (pdf/obj/stl)."""
    _ensure_dirs()
    return REPORT_DIR / f"{project_id}.{ext}"


def public_url(project_id: str) -> str:
    """
    URL segura de descarga (Report_Generated_URL, crítico #5).
    En el MVP apunta al endpoint local; en producción sería una URL S3 firmada.
    """
    return f"/api/reports/{project_id}/report.pdf"
