"""
Servicio de Precios BioKidney-AI
--------------------------------
Traduce la complejidad del modelo a créditos de cómputo y a USD.
Usado en la Pantalla 3 (Preview) para mostrar el "cálculo de recursos
requeridos" antes del cobro (Pantalla 4).
"""

from typing import Dict, Any

# Económica del clúster (MVP)
USD_PER_CREDIT = 2.50          # precio por crédito de cómputo
BASE_CREDITS = 40.0           # coste fijo de arrancar un job
CREDITS_PER_MB = 0.20         # coste por MB del dataset biológico
CREDITS_PER_KNODE = 6.0       # coste por cada 1000 nodos vasculares estimados
GFR_COMPLEXITY_REF = 62.5     # TFG de referencia (riñón nativo)


def estimate_nodes(dataset_size_bytes: int, gfr_constant: float) -> int:
    """
    Estima los nodos del árbol vascular a sintetizar según el tamaño del
    dataset y la TFG objetivo (a mayor filtración, mayor densidad capilar).
    """
    size_mb = max(dataset_size_bytes, 0) / (1024 * 1024)
    gfr_factor = max(gfr_constant, 1.0) / GFR_COMPLEXITY_REF
    nodes = int((2000 + size_mb * 800) * gfr_factor)
    return max(nodes, 500)


def quote(dataset_size_bytes: int, gfr_constant: float) -> Dict[str, Any]:
    """
    Devuelve la cotización de recursos para un proyecto.
    """
    size_mb = max(dataset_size_bytes, 0) / (1024 * 1024)
    nodes = estimate_nodes(dataset_size_bytes, gfr_constant)

    credits = (
        BASE_CREDITS
        + size_mb * CREDITS_PER_MB
        + (nodes / 1000.0) * CREDITS_PER_KNODE
    )
    credits = round(credits, 1)
    usd = round(credits * USD_PER_CREDIT, 2)

    # Estimación de tiempo de clúster (informativa)
    est_minutes = round(2 + nodes / 1500.0, 1)

    return {
        "estimated_nodes": nodes,
        "dataset_size_mb": round(size_mb, 2),
        "simulation_credits_cost": credits,
        "usd_per_credit": USD_PER_CREDIT,
        "total_billed_usd": usd,
        "estimated_cluster_minutes": est_minutes,
        "breakdown": {
            "base": BASE_CREDITS,
            "dataset": round(size_mb * CREDITS_PER_MB, 1),
            "vascular_nodes": round((nodes / 1000.0) * CREDITS_PER_KNODE, 1),
        },
    }
