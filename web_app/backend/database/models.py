import os
from datetime import datetime
from sqlalchemy import Column, Integer, String, Float, DateTime, JSON, ForeignKey, create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship

Base = declarative_base()

class SimulationRecord(Base):
    """
    Historial de simulaciones clínicas.
    Esencial para trazabilidad médica y revisión de inversionistas.
    """
    __tablename__ = "simulaciones"
    
    id = Column(Integer, primary_key=True, index=True)
    timestamp = Column(DateTime, default=datetime.utcnow)
    tipo_simulacion = Column(String) # 'Filtración', 'Oxígeno', 'Tubular'
    parametros_entrada = Column(JSON)
    resultados = Column(JSON)
    duracion_ms = Column(Float)
    estado = Column(String) # 'COMPLETED', 'FAILED'

class SystemMetric(Base):
    """
    Métricas de performance del sistema.
    """
    __tablename__ = "metricas_sistema"
    
    id = Column(Integer, primary_key=True, index=True)
    timestamp = Column(DateTime, default=datetime.utcnow)
    nombre_metrica = Column(String) # 'cpu_usage', 'mem_usage', 'integration_steps'
    valor = Column(Float)

class AdminLog(Base):
    """
    Trazas de auditoría para el módulo de administración.
    """
    __tablename__ = "admin_logs"
    
    id = Column(Integer, primary_key=True, index=True)
    timestamp = Column(DateTime, default=datetime.utcnow)
    usuario = Column(String)
    accion = Column(String)
    detalle = Column(String)

# Configuración de base de datos
DATABASE_URL = os.environ.get("DATABASE_URL", "sqlite:///web_app/database/biokidney.db")
engine = create_engine(DATABASE_URL, connect_args={"check_same_thread": False})
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

def init_db():
    Base.metadata.create_all(bind=engine)
