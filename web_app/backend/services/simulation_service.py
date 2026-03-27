import time
from typing import Dict, Any
from web_app.backend.utils.logger import bk_logger
from web_app.backend.database.models import SessionLocal, SimulationRecord
from biokidney.aggregator import BioKidneyEngine
from biokidney.core.config import cfg_physio

class SimulationService:
    """
    Capa de servicio que orquesta las simulaciones del core biokidney
    y persiste los resultados para auditoría.
    """
    def __init__(self):
        self.engine = BioKidneyEngine()

    async def run_full_pipeline(self, user_params: Dict[str, Any]) -> Dict[str, Any]:
        start_time = time.time()
        bk_logger.info(f"Iniciando simulación solicitada con parámetros: {user_params}")
        
        try:
            # Aquí se podrían sobreescribir los parámetros de config si el usuario envía cambios
            result_tfg = self.engine.ejecutar_pipeline_completo()
            
            duration = (time.time() - start_time) * 1000
            
            resultados = {
                "tfg_final": result_tfg,
                "viabilidad": 98.5, # Ejemplo
                "status": "ÓPTIMO",
                "duration_ms": duration
            }
            
            # Persistir en BD
            db = SessionLocal()
            record = SimulationRecord(
                tipo_simulacion="Full Pipeline",
                parametros_entrada=user_params,
                resultados=resultados,
                duracion_ms=duration,
                estado="COMPLETED"
            )
            db.add(record)
            db.commit()
            db.close()
            
            bk_logger.success(f"Simulación completada con éxito en {duration:.2f}ms")
            return resultados
            
        except Exception as e:
            bk_logger.error(f"Error en la ejecución del pipeline: {str(e)}")
            raise e

    def get_history(self, limit: int = 10):
        db = SessionLocal()
        history = db.query(SimulationRecord).order_by(SimulationRecord.timestamp.desc()).limit(limit).all()
        db.close()
        return history
