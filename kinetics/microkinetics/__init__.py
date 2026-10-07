"""Blocks 5–6 — Rate constants: barrier models, family parameters, named models and Arrhenius fits.

    barriers.py    ΔG‡ as a function of ΔG_rxn (capped BEP, Marcus, Agmon-Levine, Blowers-Masel, two-parabola)
    parameters.py  intrinsic barrier per reaction family (parameter sets)
    rates.py       Eyring TST + detailed balance → k_f, k_r
    models.py      ModelSpec and the registry of named models to compare and select
    arrhenius.py   modified-Arrhenius regression of k(T)

Source: Y. Alcaraz Galván
"""

from kinetics.microkinetics.barriers import (
    BARRIER_MODELS, REVERSAL_INVARIANT, calculate_barrier, invert_marcus_barrier, reversed_params,
)
from kinetics.microkinetics.parameters import (
    FAMILY_BEP_PARAMETERS, LEVEL1_PARAMETERS, REFERENCE_BEP_PARAMETERS, normalize_family_params,
)
from kinetics.microkinetics.rates import (
    KINETIC_MODELS, calculate_eyring_rate, calculate_network_rates, calculate_rate_constants,
)
from kinetics.microkinetics.models import MODELS, ModelSpec, describe_models, get_model, tabulate_family_barriers
from kinetics.microkinetics.arrhenius import build_arrhenius_table, fit_modified_arrhenius

__all__ = [
    "BARRIER_MODELS", "REVERSAL_INVARIANT", "calculate_barrier", "invert_marcus_barrier", "reversed_params",
    "FAMILY_BEP_PARAMETERS", "LEVEL1_PARAMETERS", "REFERENCE_BEP_PARAMETERS", "normalize_family_params",
    "KINETIC_MODELS", "calculate_eyring_rate", "calculate_network_rates", "calculate_rate_constants",
    "MODELS", "ModelSpec", "describe_models", "get_model", "tabulate_family_barriers",
    "build_arrhenius_table", "fit_modified_arrhenius",
]
