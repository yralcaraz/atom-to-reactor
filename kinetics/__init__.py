"""Atom-to-Reactor: Multiscale Microkinetics, Operando Spectroscopy & Chemometrics Package.

Department of Chemistry – Ångström Laboratory, Uppsala University.
Author: Yeray Alcaraz Galván
"""

from kinetics.thermo.gas_thermo import calculate_gas_thermo
from kinetics.thermo.standard_state import calculate_standard_state_shift
from kinetics.thermo.solution_gibbs import calculate_solution_gibbs
from kinetics.thermo.reaction_thermo import calculate_reaction_thermo
from kinetics.thermo.species_data import load_default_species_database

from kinetics.microkinetics.rate_constants import calculate_rate_constants, bep_eyring, marcus_eyring, level1_eyring
from kinetics.microkinetics.barrier_models import BARRIER_MODELS, barrier, invert_marcus
from kinetics.microkinetics.kinetic_parameters import (
    DEFAULT_FAMILY_BEP_PARAMETERS,
    PETER_REFERENCE_PARAMETERS,
    LEVEL1_FAMILY_PARAMETERS
)
from kinetics.microkinetics.arrhenius import fit_modified_arrhenius, generate_arrhenius_summary

from kinetics.reactor.batch_reactor import simulate_tank_reactor, build_stoichiometric_matrix
from kinetics.reactor.protocol_reactor import (
    simulate_protocol_reactor,
    compute_recipe_molarities,
    build_default_protocol_schedule
)

from kinetics.spectroscopy.molecular_symmetry import (
    build_referenced_nmr_sites,
    build_bond_graph,
    find_equivalence_classes,
    extract_nmr_sites,
    DEFAULT_NUCLEI,
    DEFAULT_SPECIES
)
from kinetics.spectroscopy.multinuclear_nmr import (
    simulate_multinuclear_spectra,
    simulate_virtual_nmr,
    lorentzian,
    DEFAULT_NMR_29SI,
    site_peaks,
    auto_regions,
    compute_water_mass_balance
)
from kinetics.spectroscopy.reaction_fingerprints import (
    build_multinuclear_feature_space,
    build_pure_component_matrix,
    build_reaction_fingerprints,
    analyze_reaction_identifiability,
    recover_reaction_extents
)

__all__ = [
    "calculate_gas_thermo",
    "calculate_standard_state_shift",
    "calculate_solution_gibbs",
    "calculate_reaction_thermo",
    "load_default_species_database",
    "calculate_rate_constants",
    "bep_eyring",
    "marcus_eyring",
    "level1_eyring",
    "BARRIER_MODELS",
    "barrier",
    "invert_marcus",
    "DEFAULT_FAMILY_BEP_PARAMETERS",
    "PETER_REFERENCE_PARAMETERS",
    "LEVEL1_FAMILY_PARAMETERS",
    "fit_modified_arrhenius",
    "generate_arrhenius_summary",
    "simulate_tank_reactor",
    "build_stoichiometric_matrix",
    "simulate_protocol_reactor",
    "compute_recipe_molarities",
    "build_default_protocol_schedule",
    "build_referenced_nmr_sites",
    "build_bond_graph",
    "find_equivalence_classes",
    "extract_nmr_sites",
    "DEFAULT_NUCLEI",
    "DEFAULT_SPECIES",
    "simulate_virtual_nmr",
    "lorentzian",
    "DEFAULT_NMR_29SI",
    "simulate_multinuclear_spectra",
    "site_peaks",
    "auto_regions",
    "compute_water_mass_balance",
    "build_multinuclear_feature_space",
    "build_pure_component_matrix",
    "build_reaction_fingerprints",
    "analyze_reaction_identifiability",
    "recover_reaction_extents"
]
