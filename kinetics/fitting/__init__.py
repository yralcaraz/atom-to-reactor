"""Block 14 — Fitting: the models against the lab NMR observables.

    residuals.py    lab observables against a model under a scenario of the unknown ages: residuals and tables
    candidates.py   candidate structures of the first fit: which barriers and energies are free
    estimation.py   best fit, profile intervals, leave-one-out and what the data cannot determine
    reporting.py    a stored fit against the measurements, spectrum by spectrum and in time
    config.py       the setup of one fit run (solvation medium, readings, paper observations), stored with its results

Source: Y. Alcaraz Galván
"""

from kinetics.fitting.residuals import (
    FREE_AGE_BOUNDS_H, REPORT_SOLVER, SCENARIO_AGES_H, SCENARIOS, SEARCH_SOLVER, build_sample_history,
    build_sample_set, calculate_predicted_shares, calculate_residuals, find_free_age, get_sample_overrides,
    merge_sample_overrides, set_sample_overrides, summarize_residuals, tabulate_residuals,
)
from kinetics.fitting.candidates import (
    FIT_STRUCTURES, FitParameter, FitStructure, build_fitted_model, embed_parent_theta, get_structure,
    tabulate_reaction_barriers,
)
from kinetics.fitting.config import (
    DEFAULT_FIT_CONFIG, FIT_CONFIG_FILE, apply_fit_config, check_fit_medium, check_fit_setup, load_fit_config,
)
from kinetics.fitting.estimation import (
    DELTA_CHI2_95, FIT_RULE_MAX_Z, FitProblem, add_profile_points, build_profile_grid, calculate_parameter_directions,
    calculate_prediction_band, calculate_profile, compare_structures, extend_profile, find_confidence_interval,
    fit_leave_one_out,
    fit_structure, get_profile_theta, load_fit_result, refine_profile_edges, simulate_synthetic_shares,
    write_fit_result,
)

__all__ = [
    "DEFAULT_FIT_CONFIG", "FIT_CONFIG_FILE", "apply_fit_config", "check_fit_medium", "check_fit_setup", "load_fit_config",
    "FIT_STRUCTURES", "FitParameter", "FitStructure", "build_fitted_model", "embed_parent_theta", "get_structure",
    "tabulate_reaction_barriers",
    "DELTA_CHI2_95", "FIT_RULE_MAX_Z", "FitProblem", "add_profile_points", "build_profile_grid",
    "calculate_parameter_directions", "extend_profile",
    "calculate_prediction_band", "calculate_profile", "compare_structures", "find_confidence_interval",
    "fit_leave_one_out", "fit_structure", "get_profile_theta", "load_fit_result", "refine_profile_edges",
    "simulate_synthetic_shares", "write_fit_result",
    "FREE_AGE_BOUNDS_H", "REPORT_SOLVER", "SCENARIO_AGES_H", "SCENARIOS", "SEARCH_SOLVER", "build_sample_history",
    "build_sample_set", "calculate_predicted_shares", "calculate_residuals", "find_free_age", "get_sample_overrides",
    "merge_sample_overrides", "set_sample_overrides", "summarize_residuals", "tabulate_residuals",
]
