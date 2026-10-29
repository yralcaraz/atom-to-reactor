"""atom-to-reactor: multiscale microkinetics, operando NMR and chemometrics for TMSPA in wet EC.

Department of Chemistry – Ångström Laboratory, Uppsala University. Author: Yeray Alcaraz Galván.

The package follows the pipeline, one folder per stage:

    kinetics/constants.py   physical constants and unit conversions
    kinetics/data/          Block 1     inputs: snapshot, species, network, experimental data
    kinetics/thermo/        Blocks 2–4  free energies in solution, reaction thermodynamics
    kinetics/microkinetics/ Blocks 5–6  barrier models, rate constants, named models, Arrhenius
    kinetics/reactor/       Blocks 7,10 batch and protocol reactors, checks against experiment
    kinetics/fitting/       Block 14    the fit of the models to the lab NMR observables

Figures and display tables for the notebooks live in demo/ (no science there). The most used functions
are re-exported here.

Source: Y. Alcaraz Galván
"""

from kinetics.data import (
    NETWORK, NETWORK_SPECIES, format_equation, load_experimental_data,
)
from kinetics.thermo import (
    calculate_network_thermo, calculate_reaction_thermo, calculate_solution_gibbs,
)
from kinetics.microkinetics import (
    MODELS, ModelSpec, calculate_network_rates, calculate_rate_constants, get_model,
)
from kinetics.reactor import (
    build_protocol_schedule, calculate_recipe_molarities, simulate_protocol, tabulate_acquisitions,
)

__all__ = [
    "NETWORK", "NETWORK_SPECIES", "format_equation", "load_experimental_data",
    "calculate_network_thermo", "calculate_reaction_thermo", "calculate_solution_gibbs",
    "MODELS", "ModelSpec", "calculate_network_rates", "calculate_rate_constants", "get_model",
    "build_protocol_schedule", "calculate_recipe_molarities", "simulate_protocol", "tabulate_acquisitions",
]
