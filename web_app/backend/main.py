"""
BioKidney AI — API Backend (FastAPI)
Plataforma médica e inversora para validación in-silico de riñón bioimpreso.
"""

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse
from pydantic import BaseModel
from typing import Dict, Any, Optional
from web_app.backend.services.simulation_service import SimulationService
from web_app.backend.database.models import init_db
from web_app.backend.utils.logger import bk_logger

app = FastAPI(title="BioKidney AI — Platform API", version="2026.2.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

init_db()
sim_service = SimulationService()


# --- Request Models ---
class SimulationRequest(BaseModel):
    params: Dict[str, Any] = {}


# --- Frontend SPA ---
@app.get("/", response_class=HTMLResponse)
async def serve_spa():
    with open("web_app/frontend/dashboard.html", "r", encoding="utf-8") as f:
        return f.read()


# --- Individual Module Endpoints ---
@app.post("/api/simulate/vascular")
async def simulate_vascular(req: SimulationRequest):
    try:
        return await sim_service.run_vascular(req.params)
    except Exception as e:
        bk_logger.error(f"Vascular error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/simulate/filtration")
async def simulate_filtration(req: SimulationRequest):
    try:
        return await sim_service.run_filtration(req.params)
    except Exception as e:
        bk_logger.error(f"Filtration error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/simulate/reabsorption")
async def simulate_reabsorption(req: SimulationRequest):
    try:
        return await sim_service.run_reabsorption(req.params)
    except Exception as e:
        bk_logger.error(f"Reabsorption error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/simulate/oxygen")
async def simulate_oxygen(req: SimulationRequest):
    try:
        return await sim_service.run_oxygen(req.params)
    except Exception as e:
        bk_logger.error(f"Oxygen error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/simulate/ipsc")
async def simulate_ipsc(req: SimulationRequest):
    try:
        return await sim_service.run_ipsc(req.params)
    except Exception as e:
        bk_logger.error(f"iPSC error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/simulate/bioprinting")
async def simulate_bioprinting(req: SimulationRequest):
    try:
        return await sim_service.run_bioprinting(req.params)
    except Exception as e:
        bk_logger.error(f"Bioprinting error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# --- Full Pipeline ---
@app.post("/api/simulate")
async def run_full_pipeline(req: SimulationRequest):
    bk_logger.info("Petición de pipeline completo recibida")
    try:
        return await sim_service.run_full_pipeline(req.params)
    except Exception as e:
        bk_logger.error(f"Pipeline error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# --- Utility Endpoints ---
@app.get("/api/history")
async def get_history():
    return sim_service.get_history()


@app.get("/api/system/health")
async def system_health():
    return {
        "status": "Healthy",
        "version": "2026.2.0",
        "engine": "BioKidney-Aggregator-MoE-v2",
        "database": "SQLite-Active",
        "modules": ["vascular", "filtration", "reabsorption", "oxygen", "ipsc", "bioprinting"]
    }


if __name__ == "__main__":
    import uvicorn
    bk_logger.success("Iniciando BioKidney AI Platform v2026.2.0...")
    uvicorn.run(app, host="0.0.0.0", port=8000)
