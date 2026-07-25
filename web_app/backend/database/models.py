import os
import uuid
from datetime import datetime
from sqlalchemy import Column, Integer, String, Float, DateTime, JSON, ForeignKey, create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship

Base = declarative_base()


def _new_id() -> str:
    """Identificador público corto para proyectos (trazabilidad clínica)."""
    return uuid.uuid4().hex[:12]


# ──────────────────────────────────────────────────────────────
# MODELO DE NEGOCIO — Funnel del MVP (Spec Maestra V1)
# ──────────────────────────────────────────────────────────────

class Project(Base):
    """
    Proyecto de simulación in-silico subido por un investigador.

    Contiene los 5 datos críticos exigidos por la spec:
      #1 User_Institution_Email   -> user_institution_email
      #2 Biological_Dataset_File  -> biological_dataset_path
      #4 Simulation_Credits_Cost  -> simulation_credits_cost
      #5 Report_Generated_URL     -> report_url
    (#3 Transaction_Stripe_ID vive en la tabla Transaction relacionada.)

    Máquina de estados (regla de negocio: no se procesa sin pago):
      DRAFT -> PREVIEWED -> AWAITING_PAYMENT -> PAID -> PROCESSING -> COMPLETED
                                                     \\-> FAILED
    """
    __tablename__ = "proyectos"

    id = Column(String, primary_key=True, default=_new_id, index=True)
    created_at = Column(DateTime, default=datetime.utcnow)

    # --- Pantalla 2: Configuración (Input) ---
    project_name = Column(String, nullable=False)
    institution = Column(String, nullable=False)
    user_institution_email = Column(String, nullable=False)  # CRÍTICO #1
    gfr_constant = Column(Float, default=62.5)               # TFG objetivo (mL/min)
    sim_params = Column(JSON, default=dict)                   # constantes extra del modelo

    # CRÍTICO #2 — archivo biológico (≤500MB, validado en el router)
    biological_dataset_path = Column(String, nullable=True)
    dataset_filename = Column(String, nullable=True)
    dataset_size_bytes = Column(Integer, default=0)

    # --- Pantalla 3: Preview + cálculo de recursos ---
    simulation_credits_cost = Column(Float, default=0.0)     # CRÍTICO #4
    total_billed_usd = Column(Float, default=0.0)
    estimated_nodes = Column(Integer, default=0)

    # --- Pantalla 5: Resultados ---
    report_url = Column(String, nullable=True)               # CRÍTICO #5
    results = Column(JSON, nullable=True)
    progress = Column(Integer, default=0)                    # 0-100 barra de progreso

    # Estado del funnel / clúster
    status = Column(String, default="DRAFT", index=True)

    transactions = relationship("Transaction", back_populates="project")

    def to_dict(self):
        return {
            "id": self.id,
            "created_at": self.created_at.isoformat() if self.created_at else None,
            "project_name": self.project_name,
            "institution": self.institution,
            "user_institution_email": self.user_institution_email,
            "gfr_constant": self.gfr_constant,
            "dataset_filename": self.dataset_filename,
            "dataset_size_bytes": self.dataset_size_bytes,
            "simulation_credits_cost": self.simulation_credits_cost,
            "total_billed_usd": self.total_billed_usd,
            "estimated_nodes": self.estimated_nodes,
            "report_url": self.report_url,
            "progress": self.progress,
            "status": self.status,
        }


class Transaction(Base):
    """
    Registro de cobro Stripe (Pantalla 4 — MOMENTO DE COBRO).

    CRÍTICO #3: Transaction_Stripe_ID -> stripe_transaction_id
    El estado SUCCESS SOLO lo fija el webhook de Stripe (o el bypass
    local en modo test). Ese evento es el único que habilita el
    procesamiento en el clúster.
    """
    __tablename__ = "transacciones"

    id = Column(Integer, primary_key=True, index=True)
    created_at = Column(DateTime, default=datetime.utcnow)
    project_id = Column(String, ForeignKey("proyectos.id"), index=True)

    stripe_transaction_id = Column(String, index=True)  # CRÍTICO #3 (payment_intent id)
    amount_usd = Column(Float, default=0.0)
    currency = Column(String, default="usd")

    # Facturación de subvención / laboratorio
    billing_name = Column(String)         # Nombre en tarjeta
    billing_address = Column(String)
    tax_id = Column(String)               # NIT / Tax ID

    status = Column(String, default="PENDING")  # PENDING -> SUCCESS / FAILED

    project = relationship("Project", back_populates="transactions")

    def to_dict(self):
        return {
            "id": self.id,
            "created_at": self.created_at.isoformat() if self.created_at else None,
            "project_id": self.project_id,
            "stripe_transaction_id": self.stripe_transaction_id,
            "amount_usd": self.amount_usd,
            "currency": self.currency,
            "billing_name": self.billing_name,
            "tax_id": self.tax_id,
            "status": self.status,
        }

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
