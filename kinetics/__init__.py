"""atom-to-reactor: multiscale microkinetics, operando NMR and chemometrics for TMSPA in wet EC.

Department of Chemistry – Ångström Laboratory, Uppsala University. Author: Yeray Alcaraz Galván.

The package follows the pipeline, one folder per stage:

    kinetics/constants.py   physical constants and unit conversions
    kinetics/data/          Block 1     inputs: snapshot, species, network, experimental data
    kinetics/thermo/        Blocks 2–4  free energies in solution, reaction thermodynamics
    kinetics/microkinetics/ Blocks 5–6  barrier models, rate constants, named models, Arrhenius
    kinetics/reactor/       Blocks 7,10 batch and protocol reactors, checks against experiment
    kinetics/spectroscopy/  Blocks 8,11–13 NMR shifts, spectra, fingerprints
    kinetics/fitting/       Block 14    bounds from windowed data, experiment design

Figures and display tables for the notebooks live in demo/ (no science there). The most used functions
are re-exported here.
"""

from kinetics.data import (
    NETWORK, NETWORK_SPECIES, format_equation, get_barrier_windows, get_measured_shifts, load_experimental_data,
    load_species_database,
)
from kinetics.thermo import (
    calculate_cycle_residuals, calculate_network_thermo, calculate_reaction_thermo, calculate_solution_gibbs,
)
from kinetics.microkinetics import (
    MODELS, ModelSpec, build_arrhenius_table, calculate_network_rates, calculate_rate_constants, describe_models,
    get_model,
)
from kinetics.reactor import (
    build_protocol_schedule, calculate_recipe_molarities, calculate_remaining_fraction,
    evaluate_control_experiments, simulate_batch_reactor, simulate_protocol, summarize_batch_runs,
    tabulate_acquisitions, tabulate_trajectory,
)
from kinetics.spectroscopy import (
    build_nmr_catalog, calculate_nmr_peaks, calculate_water_mass_balance, recover_reaction_extents,
    run_fingerprint_analysis, simulate_acquisition_spectra, simulate_nmr_spectra,
)

__all__ = [
    "NETWORK", "NETWORK_SPECIES", "format_equation", "get_barrier_windows", "get_measured_shifts",
    "load_experimental_data", "load_species_database",
    "calculate_cycle_residuals", "calculate_network_thermo", "calculate_reaction_thermo", "calculate_solution_gibbs",
    "MODELS", "ModelSpec", "build_arrhenius_table", "calculate_network_rates", "calculate_rate_constants",
    "describe_models", "get_model",
    "build_protocol_schedule", "calculate_recipe_molarities", "calculate_remaining_fraction",
    "evaluate_control_experiments", "simulate_batch_reactor", "simulate_protocol", "summarize_batch_runs",
    "tabulate_acquisitions", "tabulate_trajectory",
    "build_nmr_catalog", "calculate_nmr_peaks", "calculate_water_mass_balance", "recover_reaction_extents",
    "run_fingerprint_analysis", "simulate_acquisition_spectra", "simulate_nmr_spectra",
]
