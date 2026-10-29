"""Blocks 5–6 — Rate constants: barrier models, family parameters, named models.

    barriers.py    ΔG‡ as a function of ΔG_rxn (capped BEP, Marcus)
    parameters.py  intrinsic barrier per reaction family (parameter sets)
    rates.py       Eyring TST + detailed balance → k_f, k_r
    models.py      ModelSpec and the registry of named models to compare and select

Source: Y. Alcaraz Galván
"""

from kinetics.microkinetics.barriers import (
    BARRIER_MODELS, calculate_barrier, reversed_params,
)
from kinetics.microkinetics.parameters import (
    LEVEL1_PARAMETERS, REFERENCE_BEP_PARAMETERS, normalize_family_params,
)
from kinetics.microkinetics.rates import (
    KINETIC_MODELS, calculate_network_rates, calculate_rate_constants,
)
from kinetics.microkinetics.models import MODELS, ModelSpec, get_model

__all__ = [
    "BARRIER_MODELS", "calculate_barrier", "reversed_params",
    "LEVEL1_PARAMETERS", "REFERENCE_BEP_PARAMETERS", "normalize_family_params",
    "KINETIC_MODELS", "calculate_network_rates", "calculate_rate_constants",
    "MODELS", "ModelSpec", "get_model",
]
