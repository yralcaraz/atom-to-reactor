"""Blocks 7 and 10 — Reactors: mass-action ODEs integrated in time, and checks against experiment.

    engine.py       mass-action right-hand side, analytic Jacobian and element balances
    protocol.py     recipe, temperature program, multi-stage protocol reactor and sample histories
    validation.py   the Gogoi 2024 observables as functions of a state, and the water-series recipe

Source: Y. Alcaraz Galván
"""

from kinetics.reactor.engine import IntegrationBudgetExceeded, MassActionSystem, calculate_element_totals
from kinetics.reactor.protocol import (
    Stage, build_protocol_schedule, calculate_recipe_molarities, simulate_history, simulate_protocol, tabulate_acquisitions,
)
from kinetics.reactor.validation import OBSERVABLES, calculate_water_series_c0

__all__ = [
    "IntegrationBudgetExceeded", "MassActionSystem", "calculate_element_totals",
    "Stage", "build_protocol_schedule", "calculate_recipe_molarities", "simulate_history", "simulate_protocol",
    "tabulate_acquisitions",
    "OBSERVABLES", "calculate_water_series_c0",
]
