"""
Servicio de Pago BioKidney-AI (Pantalla 4 — MOMENTO DE COBRO ⚡)
---------------------------------------------------------------
Integración estructural con Stripe lista para producción, con un BYPASS
local para modo test (sin claves ni librería instalada).

Modo de operación (env BIOKIDNEY_PAYMENTS):
  - "stripe"  -> usa el SDK real de Stripe (requiere STRIPE_SECRET_KEY).
  - "bypass"  -> (por defecto) simula PaymentIntent + webhook localmente,
                 permitiendo probar la simulación real de inmediato.

En AMBOS modos el flujo de gating es idéntico: el procesamiento del clúster
solo se dispara cuando el webhook (real o simulado) confirma Payment_Success.
"""

import os
import time
import uuid
from typing import Dict, Any, Optional

from web_app.backend.utils.logger import bk_logger

PAYMENT_MODE = os.environ.get("BIOKIDNEY_PAYMENTS", "bypass").lower()
STRIPE_SECRET_KEY = os.environ.get("STRIPE_SECRET_KEY", "")
STRIPE_PUBLISHABLE_KEY = os.environ.get("STRIPE_PUBLISHABLE_KEY", "pk_test_PLACEHOLDER")
STRIPE_WEBHOOK_SECRET = os.environ.get("STRIPE_WEBHOOK_SECRET", "")

# Import opcional: el bypass funciona aunque 'stripe' no esté instalado.
_stripe = None
if PAYMENT_MODE == "stripe":
    try:
        import stripe as _stripe  # type: ignore
        _stripe.api_key = STRIPE_SECRET_KEY
        bk_logger.success("PaymentService: modo STRIPE real activado.")
    except ImportError:
        bk_logger.warning(
            "PaymentService: 'stripe' no instalado; cayendo a modo BYPASS."
        )
        _stripe = None
        PAYMENT_MODE = "bypass"
else:
    bk_logger.info("PaymentService: modo BYPASS local (test) activado.")


def publishable_key() -> str:
    return STRIPE_PUBLISHABLE_KEY


def is_bypass() -> bool:
    return PAYMENT_MODE != "stripe" or _stripe is None


def create_payment_intent(amount_usd: float, project_id: str,
                          metadata: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Crea el PaymentIntent (Pantalla 4). Devuelve el client_secret que el
    frontend (Stripe Elements) usa para confirmar el pago.
    """
    amount_cents = int(round(amount_usd * 100))
    metadata = {"project_id": project_id, **(metadata or {})}

    if not is_bypass():
        intent = _stripe.PaymentIntent.create(
            amount=amount_cents,
            currency="usd",
            metadata=metadata,
            automatic_payment_methods={"enabled": True},
        )
        return {
            "payment_intent_id": intent.id,
            "client_secret": intent.client_secret,
            "mode": "stripe",
        }

    # --- BYPASS local ---
    intent_id = f"pi_bypass_{uuid.uuid4().hex[:16]}"
    bk_logger.info(f"[BYPASS] PaymentIntent simulado {intent_id} por ${amount_usd}")
    return {
        "payment_intent_id": intent_id,
        "client_secret": f"{intent_id}_secret_test",
        "mode": "bypass",
    }


def confirm_bypass_payment(payment_intent_id: str) -> Dict[str, Any]:
    """
    Simula la confirmación de tarjeta en modo bypass. Emula el evento
    'payment_intent.succeeded' que en producción llegaría al webhook.
    """
    time.sleep(0.4)  # latencia bancaria simulada
    return {
        "type": "payment_intent.succeeded",
        "data": {"object": {"id": payment_intent_id, "status": "succeeded"}},
    }


def verify_webhook(payload: bytes, sig_header: Optional[str]) -> Dict[str, Any]:
    """
    Verifica y parsea un evento de webhook de Stripe.
    En modo bypass acepta el payload JSON tal cual (sin firma).
    """
    if not is_bypass() and STRIPE_WEBHOOK_SECRET:
        event = _stripe.Webhook.construct_event(
            payload, sig_header, STRIPE_WEBHOOK_SECRET
        )
        return event
    # BYPASS / sin secret configurado: parseo directo.
    import json
    return json.loads(payload.decode("utf-8"))
