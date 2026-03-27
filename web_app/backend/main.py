from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from typing import List, Dict, Any
from web_app.backend.services.simulation_service import SimulationService
from web_app.backend.database.models import init_db
from web_app.backend.utils.logger import bk_logger

app = FastAPI(title="BioKidney API - Medical & Investor Platform")

# Configurar CORS para frontend moderno
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Inicializar Base de Datos
init_db()

# Instanciar servicios
sim_service = SimulationService()

# --- MODELOS API ---
class SimulationRequest(BaseModel):
    params: Dict[str, Any]

# --- ENDPOINTS ---

@app.get("/", response_class=HTMLResponse)
async def read_root():
    """Sirve la SPA del dashboard."""
    with open("web_app/frontend/dashboard.html", "r") as f:
        return f.read()

@app.post("/api/simulate")
async def run_simulation(req: SimulationRequest):
    bk_logger.info("Recibida petición de simulación web")
    try:
        results = await sim_service.run_full_pipeline(req.params)
        return results
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/history")
async def get_history():
    return sim_service.get_history()

@app.get("/api/system/health")
async def system_health():
    return {
        "status": "Healthy",
        "version": "2026.1.0",
        "engine": "BioKidney-Aggregator-v1",
        "database": "SQLite-Active"
    }

# --- MÓDULO ADMIN ---
@app.get("/api/admin/metrics")
async def get_admin_metrics():
    """Métricas de uso del sistema para administradores."""
    return {
        "total_simulations": 154, # Ejemplo
        "avg_duration_ms": 12.4,
        "active_users": 3,
        "logs_path": "web_app/logs/audit.log"
    }

if __name__ == "__main__":
    import uvicorn
    bk_logger.success("Iniciando Servidor Web BioKidney AI...")
    uvicorn.run(app, host="0.0.0.0", port=8000)
