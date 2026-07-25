"""
Router de Pagos — Pantalla 4 (MOMENTO DE COBRO ⚡)
--------------------------------------------------
Flujo:
  1. POST /api/payments/intent   -> crea PaymentIntent (Stripe o bypass)
  2. Frontend confirma la tarjeta (Stripe Elements) o, en bypass,
     POST /api/payments/confirm-bypass simula la confirmación.
  3. POST /api/payments/webhook  -> ÚNICO punto que marca SUCCESS y
     dispara el procesamiento del clúster.

REGLA DE NEGOCIO: sin evento 'payment_intent.succeeded' NO hay simulación.
"""

from fastapi import APIRouter, Request, BackgroundTasks, HTTPException
from pydantic import BaseModel

from web_app.backend.database.models import SessionLocal, Project, Transaction
from web_app.backend.services import payment_service, cluster_service
from web_app.backend.utils.logger import bk_logger

router = APIRouter(prefix="/api/payments", tags=["payments"])


class IntentRequest(BaseModel):
    project_id: str
    billing_name: str
    billing_address: str = ""
    tax_id: str = ""  # NIT / Tax ID


class BypassConfirm(BaseModel):
    payment_intent_id: str


@router.get("/config")
async def payment_config():
    """Clave pública y modo, para inicializar Stripe Elements en el frontend."""
    return {
        "publishable_key": payment_service.publishable_key(),
        "bypass": payment_service.is_bypass(),
    }


@router.post("/intent")
async def create_intent(req: IntentRequest):
    """Crea el PaymentIntent y una Transaction PENDING con datos de facturación."""
    db = SessionLocal()
    try:
        proj = db.query(Project).get(req.project_id)
        if not proj:
            raise HTTPException(status_code=404, detail="Proyecto no encontrado")

        intent = payment_service.create_payment_intent(
            proj.total_billed_usd, proj.id,
            metadata={"institution": proj.institution},
        )

        tx = Transaction(
            project_id=proj.id,
            stripe_transaction_id=intent["payment_intent_id"],
            amount_usd=proj.total_billed_usd,
            billing_name=req.billing_name,
            billing_address=req.billing_address,
            tax_id=req.tax_id,
            status="PENDING",
        )
        db.add(tx)
        proj.status = "AWAITING_PAYMENT"
        db.commit()

        return {
            "payment_intent_id": intent["payment_intent_id"],
            "client_secret": intent["client_secret"],
            "mode": intent["mode"],
            "amount_usd": proj.total_billed_usd,
        }
    finally:
        db.close()


@router.post("/confirm-bypass")
async def confirm_bypass(req: BypassConfirm, request: Request,
                         background: BackgroundTasks):
    """
    Modo test: simula la confirmación de tarjeta y reenvía el evento
    al mismo pipeline de webhook para respetar el gating real.
    """
    if not payment_service.is_bypass():
        raise HTTPException(status_code=400, detail="Bypass deshabilitado (modo Stripe real).")
    event = payment_service.confirm_bypass_payment(req.payment_intent_id)
    return await _handle_event(event, background)


@router.post("/webhook")
async def stripe_webhook(request: Request, background: BackgroundTasks):
    """Endpoint que Stripe invoca. Verifica firma y procesa el evento."""
    payload = await request.body()
    sig = request.headers.get("stripe-signature")
    try:
        event = payment_service.verify_webhook(payload, sig)
    except Exception as e:
        bk_logger.error(f"Webhook inválido: {e}")
        raise HTTPException(status_code=400, detail="Firma de webhook inválida")
    return await _handle_event(event, background)


async def _handle_event(event: dict, background: BackgroundTasks):
    """Marca la transacción como SUCCESS y ARRANCA el clúster (gating)."""
    if event.get("type") != "payment_intent.succeeded":
        return {"received": True, "handled": False}

    intent_id = event["data"]["object"]["id"]
    db = SessionLocal()
    try:
        tx = db.query(Transaction).filter(
            Transaction.stripe_transaction_id == intent_id
        ).first()
        if not tx:
            bk_logger.warning(f"Webhook: transacción {intent_id} no encontrada")
            return {"received": True, "handled": False}

        if tx.status == "SUCCESS":
            return {"received": True, "handled": True, "duplicate": True}

        tx.status = "SUCCESS"
        proj = db.query(Project).get(tx.project_id)
        if proj:
            proj.status = "PAID"
            proj.progress = 5
        db.commit()
        project_id = tx.project_id
    finally:
        db.close()

    # ⚡ Único disparador del clúster, tras Payment_Success confirmado.
    background.add_task(cluster_service.process_project, project_id)
    bk_logger.success(f"Pago confirmado {intent_id} → clúster disparado para {project_id}")
    return {"received": True, "handled": True, "project_id": project_id}
