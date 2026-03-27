import logging
import sys
from pathlib import Path
from loguru import logger

# Configuración de niveles de log
LOG_LEVELS = {
    "DEBUG": "DEBUG",
    "INFO": "INFO",
    "SUCCESS": "SUCCESS",
    "WARNING": "WARNING",
    "ERROR": "ERROR",
    "CRITICAL": "CRITICAL"
}

def setup_advanced_logging():
    """
    Configura un sistema de logging profesional con Loguru.
    - Consola con colores para desarrollo.
    - Archivo rotativo para auditoría médica.
    - Estructura JSON para fácil integración con ELK/CloudWatch.
    """
    # Eliminar handlers por defecto
    logger.remove()

    # Formato Vercel/V0 Style para consola
    fmt = "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>"
    
    # Handler Consola
    logger.add(sys.stdout, format=fmt, level="DEBUG", colorize=True)

    # Handler Archivo (Auditoría)
    log_path = Path("web_app/logs/audit.log")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger.add(log_path, rotation="10 MB", retention="1 month", level="INFO", compression="zip")

    logger.success("🚀 Sistema de Logging BioKidney Activado (Nivel: DEBUG)")

# Exportar logger configurado
setup_advanced_logging()
bk_logger = logger
