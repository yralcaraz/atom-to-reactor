"""Block 14 — Fitting and experiment design: what the data constrain now and what to measure next.

    feasibility.py  bounds on the family barriers from windowed observations (Gogoi 2024 controls)
    readouts.py     NMR peaks a fit can integrate, derived from the measured shifts (Gogoi 2024)
    design.py       NMR readouts, candidate experiments and their expected parameter precision (Fisher information)
"""

from kinetics.fitting.feasibility import (
    DEFAULT_G_GRID_EV, FAMILIES, find_barrier_bounds, get_family_parameter, project_bounds_to_entropy,
    scan_family_barriers, summarize_feasible_intervals,
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

__all__ = [
    "DEFAULT_G_GRID_EV", "FAMILIES", "find_barrier_bounds", "get_family_parameter", "project_bounds_to_entropy",
    "scan_family_barriers", "summarize_feasible_intervals",
    "READOUT_NUCLEI", "build_nmr_readouts", "calculate_site_dft_shift", "find_reachable_species", "find_site_atoms",
    "find_unread_species", "tabulate_readout_evidence",
    "DEFAULT_ACQUISITION", "DEFAULT_PRIOR_SIGMA", "NMR_READOUTS", "ExperimentDesign", "build_composition",
    "add_readout_noise", "build_in_situ_design", "build_isothermal_design", "calculate_barrier_error_budget", "build_step_design", "calculate_design_information", "calculate_eyring_precision",
    "calculate_half_life_map", "calculate_parameter_precision", "calculate_readouts", "compare_designs",
    "describe_designs", "scan_design_precision", "simulate_design", "split_family_by_reaction",
]
