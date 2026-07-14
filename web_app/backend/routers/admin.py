"""
Router de Administración — Pantalla 6 (Backoffice)
Tabla de proyectos con: ID, usuario, fecha, tipo de modelo, estado y USD.
"""

from fastapi import APIRouter

from web_app.backend.database.models import SessionLocal, Project, Transaction
from sqlalchemy import func

router = APIRouter(prefix="/api/admin", tags=["admin"])


@router.get("/projects")
async def list_projects():
    """Filas del backoffice + métricas de recaudación."""
    db = SessionLocal()
    try:
        projects = db.query(Project).order_by(Project.created_at.desc()).all()
        rows = []
        for p in projects:
            paid = (
                db.query(Transaction)
                .filter(Transaction.project_id == p.id,
                        Transaction.status == "SUCCESS")
                .first()
            )
            rows.append({
                "id": p.id,
                "user": p.user_institution_email,
                "institution": p.institution,
                "date": p.created_at.isoformat() if p.created_at else None,
                "model_type": p.project_name,
                "status": p.status,
                "credits": p.simulation_credits_cost,
                "billed_usd": p.total_billed_usd if paid else 0.0,
                "paid": bool(paid),
            })

        total_revenue = (
            db.query(func.coalesce(func.sum(Transaction.amount_usd), 0.0))
            .filter(Transaction.status == "SUCCESS")
            .scalar()
        )
        return {
            "projects": rows,
            "metrics": {
                "total_projects": len(rows),
                "total_revenue_usd": round(float(total_revenue), 2),
                "completed": sum(1 for r in rows if r["status"] == "COMPLETED"),
                "processing": sum(1 for r in rows if r["status"] == "PROCESSING"),
            },
        }
    finally:
        db.close()
