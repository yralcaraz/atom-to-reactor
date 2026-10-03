"""Block 14 — Fitting and experiment design: what the data constrain now and what to measure next.

    feasibility.py  bounds on the family barriers from windowed observations (Gogoi 2024 controls)
    constraints.py  observations with a time interval (lab and paper) and the barrier values each allows
    structure.py    what a measurement implies whatever the barriers are (equilibrium, heating, water balance)
    readouts.py     NMR peaks a fit can integrate, derived from the measured shifts (Gogoi 2024)
    design.py       NMR readouts, candidate experiments and their expected parameter precision (Fisher information)
    residuals.py    lab observables against a model under a scenario of the unknown ages: residuals and tables

Classification: PUBLIC (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván
"""

from kinetics.fitting.feasibility import (
    DEFAULT_G_GRID_EV, FAMILIES, find_barrier_bounds, get_family_parameter, project_bounds_to_entropy,
    scan_family_barriers, summarize_feasible_intervals,
)
from kinetics.fitting.constraints import (
    DEFAULT_SCAN_GRID_EV, SHARED_BARRIER, build_observation, build_share_window, describe_observations,
    evaluate_observation, find_allowed_intervals, intersect_intervals, observations_from_controls,
    replace_observation, scan_barrier, select_observations, summarize_region, trace_boundaries,
)
from kinetics.fitting.structure import (
    SILYL_RELEASED, calculate_equilibrium_locus, calculate_path_distance, calculate_released_silyl_M,
    calculate_required_water_M, calculate_rt_equivalent_hours, find_closest_on_locus,
)
from kinetics.fitting.readouts import (
    READOUT_NUCLEI, build_nmr_readouts, calculate_site_dft_shift, find_reachable_species, find_site_atoms,
    find_unread_species, tabulate_readout_evidence,
)
from kinetics.fitting.design import (
    DEFAULT_ACQUISITION, DEFAULT_PRIOR_SIGMA, NMR_READOUTS, ExperimentDesign, build_composition,
    add_readout_noise, build_in_situ_design, build_isothermal_design, calculate_barrier_error_budget, build_step_design, calculate_design_information, calculate_eyring_precision,
    calculate_half_life_map, calculate_parameter_precision, calculate_readouts, compare_designs, describe_designs,
    scan_design_precision, simulate_design, split_family_by_reaction,
)

from kinetics.fitting.residuals import (
    FREE_AGE_BOUNDS_H, REPORT_SOLVER, SCENARIO_AGES_H, SCENARIOS, SEARCH_SOLVER, build_sample_history,
    build_sample_set, calculate_predicted_shares, calculate_residuals, find_free_age, summarize_residuals,
    tabulate_residuals,
)
__all__ = [
    "FREE_AGE_BOUNDS_H", "REPORT_SOLVER", "SCENARIO_AGES_H", "SCENARIOS", "SEARCH_SOLVER", "build_sample_history",
    "build_sample_set", "calculate_predicted_shares", "calculate_residuals", "find_free_age", "summarize_residuals",
    "tabulate_residuals",
    "DEFAULT_G_GRID_EV", "FAMILIES", "find_barrier_bounds", "get_family_parameter", "project_bounds_to_entropy",
    "scan_family_barriers", "summarize_feasible_intervals",
    "DEFAULT_SCAN_GRID_EV", "SHARED_BARRIER", "build_observation", "build_share_window", "describe_observations",
    "evaluate_observation", "find_allowed_intervals", "intersect_intervals", "observations_from_controls",
    "replace_observation", "scan_barrier", "select_observations", "summarize_region", "trace_boundaries",
    "SILYL_RELEASED", "calculate_equilibrium_locus", "calculate_path_distance", "calculate_released_silyl_M",
    "calculate_required_water_M", "calculate_rt_equivalent_hours", "find_closest_on_locus",
    "READOUT_NUCLEI", "build_nmr_readouts", "calculate_site_dft_shift", "find_reachable_species", "find_site_atoms",
    "find_unread_species", "tabulate_readout_evidence",
    "DEFAULT_ACQUISITION", "DEFAULT_PRIOR_SIGMA", "NMR_READOUTS", "ExperimentDesign", "build_composition",
    "add_readout_noise", "build_in_situ_design", "build_isothermal_design", "calculate_barrier_error_budget", "build_step_design", "calculate_design_information", "calculate_eyring_precision",
    "calculate_half_life_map", "calculate_parameter_precision", "calculate_readouts", "compare_designs",
    "describe_designs", "scan_design_precision", "simulate_design", "split_family_by_reaction",
]
