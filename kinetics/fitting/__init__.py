"""Block 14 — Fitting: the models against the lab NMR observables.

    residuals.py    lab observables against a model under a scenario of the unknown ages: residuals and tables
    candidates.py   candidate structures of the first fit: which barriers and energies are free
    estimation.py   best fit and profile intervals
    reporting.py    a stored fit against the measurements, spectrum by spectrum and in time
    config.py       the setup of one fit run (solvation medium, readings, paper observations), stored with its results

Source: Y. Alcaraz Galván
"""

from kinetics.fitting.residuals import (
    FREE_AGE_BOUNDS_H, REPORT_SOLVER, SCENARIO_AGES_H, SCENARIOS, SEARCH_SOLVER, build_sample_history,
    build_sample_set, calculate_predicted_shares, calculate_residuals, get_sample_overrides,
    merge_sample_overrides, set_sample_overrides, tabulate_residuals,
)
from kinetics.fitting.candidates import (
    FIT_STRUCTURES, FitParameter, FitStructure, get_structure, tabulate_reaction_barriers,
)
from kinetics.fitting.config import (
    DEFAULT_FIT_CONFIG, FIT_CONFIG_FILE, apply_fit_config, check_fit_medium, check_fit_setup, load_fit_config,
)
from kinetics.fitting.estimation import (
    DELTA_CHI2_95, FIT_RULE_MAX_Z, FitProblem, build_profile_grid, calculate_profile, find_confidence_interval,
    fit_structure, load_fit_result, write_fit_result,
)

__all__ = [
    "DEFAULT_FIT_CONFIG", "FIT_CONFIG_FILE", "apply_fit_config", "check_fit_medium", "check_fit_setup", "load_fit_config",
    "FIT_STRUCTURES", "FitParameter", "FitStructure", "get_structure", "tabulate_reaction_barriers",
    "DELTA_CHI2_95", "FIT_RULE_MAX_Z", "FitProblem", "build_profile_grid", "calculate_profile", "find_confidence_interval",
    "fit_structure", "load_fit_result", "write_fit_result",
    "FREE_AGE_BOUNDS_H", "REPORT_SOLVER", "SCENARIO_AGES_H", "SCENARIOS", "SEARCH_SOLVER", "build_sample_history",
    "build_sample_set", "calculate_predicted_shares", "calculate_residuals", "get_sample_overrides",
    "merge_sample_overrides", "set_sample_overrides", "tabulate_residuals",
]
