"""Blocks 7 and 10 — Reactors: mass-action ODEs integrated in time, and checks against experiment.

    engine.py       mass-action right-hand side, analytic Jacobian and element balances
    batch.py        isothermal batch reactor
    protocol.py     recipe, temperature program and multi-stage protocol reactor
    observables.py  characteristic times and tables extracted from a run
    validation.py   re-simulation of the Gogoi 2024 control experiments and water series for any model
"""

from kinetics.reactor.engine import MassActionSystem, calculate_element_totals
from kinetics.reactor.batch import simulate_batch_reactor
from kinetics.reactor.protocol import Stage, build_protocol_schedule, calculate_recipe_molarities, simulate_protocol
from kinetics.reactor.observables import (
    calculate_remaining_fraction, calculate_worst_case_pressure_bar, find_crossing_time, summarize_batch_runs, tabulate_acquisitions,
    tabulate_trajectory,
)
from kinetics.reactor.validation import (
    OBSERVABLE_REACTIONS, OBSERVABLES, calculate_phosphate_fractions, calculate_water_series_c0,
    evaluate_control_experiments, evaluate_heating_observation, evaluate_water_series, is_within_window,
    simulate_control_experiment, simulate_phosphate_path,
)

__all__ = [
    "MassActionSystem", "calculate_element_totals", "simulate_batch_reactor",
    "Stage", "build_protocol_schedule", "calculate_recipe_molarities", "simulate_protocol",
    "calculate_remaining_fraction", "calculate_worst_case_pressure_bar", "find_crossing_time", "summarize_batch_runs", "tabulate_acquisitions",
    "tabulate_trajectory",
    "OBSERVABLE_REACTIONS", "OBSERVABLES", "calculate_phosphate_fractions", "calculate_water_series_c0",
    "evaluate_control_experiments", "evaluate_heating_observation", "evaluate_water_series", "is_within_window",
    "simulate_control_experiment", "simulate_phosphate_path",
]
