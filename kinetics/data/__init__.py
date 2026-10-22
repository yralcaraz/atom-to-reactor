"""Block 1 — Inputs: the computed-data snapshot, the species database, the reaction network and the
experimental reference data.

    snapshot.py      read-only access to the Tank snapshot (DFT, MD solvation, NMR, geometries; read from the data folder)
    species.py       species database built from the snapshot (energies converted to eV)
    medium.py        the solvation medium the species database takes (neat solvent or EC/DEC mixture)
    network.py       the 9-reaction network, its species and stoichiometry helpers
    experimental.py  Gogoi et al. 2024: measured shifts, barrier windows, control experiments, water series
    lab_nmr.py       measured NMR spectra (JEOL .jdf, read from LAB_NMR_DIR): inventory and processing
    observables.py   curated lab observables (shares, uncertainties, sample histories) that a fit reads

Source: Y. Alcaraz Galván
"""

from kinetics.data.snapshot import (
    calculate_molar_mass, count_elements, get_dataset_records, get_solvation_records, get_solvent,
    list_solvation_sets, load_snapshot,
)
from kinetics.data.medium import DEFAULT_MEDIUM, SolvationMedium, get_active_medium, set_active_medium
from kinetics.data.species import load_medium_species_database, load_species_database
from kinetics.data.network import (
    EXTRA_REACTIONS, NETWORK, NETWORK_SPECIES, NETWORK_VARIANTS, build_network, build_stoichiometric_matrix,
    find_reaction_cycles, format_equation, get_active_network, get_active_network_variant, list_species,
    set_active_network_variant,
)
from kinetics.data.experimental import (
    get_barrier_windows, get_measured_shifts, get_water_series, load_experimental_data,
)
from kinetics.data.lab_nmr import (
    DEFAULT_LAB_NMR_DIR, build_lab_nmr_inventory, calculate_area_shares, calculate_window_integrals, estimate_noise,
    find_heated_windows, load_spectrum, load_text_spectrum, read_jdf, refine_window_phase, tabulate_area_shares,
)
from kinetics.data.observables import (
    DEFAULT_ERROR_FLOOR, DEFAULT_OBSERVABLES_PATH, build_lab_inventory_for_observables, build_lab_observables,
    calculate_replicate_scatter, describe_lab_observables, list_spectrum_sets, load_lab_observables,
    load_lab_sample_folders, load_share_tables, tabulate_observable_shares, write_lab_observables,
)

__all__ = [
    "load_snapshot", "get_dataset_records", "get_solvation_records", "get_solvent", "list_solvation_sets", "count_elements", "calculate_molar_mass",
    "DEFAULT_MEDIUM", "SolvationMedium", "get_active_medium", "set_active_medium",
    "load_medium_species_database", "load_species_database",
    "EXTRA_REACTIONS", "NETWORK", "NETWORK_SPECIES", "NETWORK_VARIANTS", "build_network", "build_stoichiometric_matrix",
    "find_reaction_cycles", "format_equation", "get_active_network", "get_active_network_variant", "list_species",
    "set_active_network_variant",
    "load_experimental_data", "get_measured_shifts", "get_barrier_windows", "get_water_series",
    "DEFAULT_LAB_NMR_DIR", "build_lab_nmr_inventory", "read_jdf", "load_spectrum", "load_text_spectrum",
    "estimate_noise", "calculate_window_integrals", "refine_window_phase", "calculate_area_shares",
    "tabulate_area_shares", "find_heated_windows",
    "DEFAULT_ERROR_FLOOR", "DEFAULT_OBSERVABLES_PATH", "build_lab_inventory_for_observables",
    "build_lab_observables", "calculate_replicate_scatter", "list_spectrum_sets",
    "describe_lab_observables", "load_lab_observables", "load_lab_sample_folders", "load_share_tables", "tabulate_observable_shares",
    "write_lab_observables",
]
