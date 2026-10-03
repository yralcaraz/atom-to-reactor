"""Block 1 — Inputs: the computed-data snapshot, the species database, the reaction network and the
experimental reference data.

    snapshot.py      read-only access to data/tank_api_snapshot.json (DFT, MD solvation, NMR, geometries)
    species.py       species database built from the snapshot (energies converted to eV)
    network.py       the 9-reaction network, its species and stoichiometry helpers
    experimental.py  Gogoi et al. 2024: measured shifts, barrier windows, control experiments, water series
    lab_nmr.py       measured lab spectra (JEOL .jdf, private, read from LAB_NMR_DIR): inventory and processing

Classification: PUBLIC (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván
"""

from kinetics.data.snapshot import (
    calculate_molar_mass, count_elements, get_dataset_records, get_solvation_records, load_snapshot,
)
from kinetics.data.species import load_species_database
from kinetics.data.network import (
    NETWORK, NETWORK_SPECIES, build_stoichiometric_matrix, find_reaction_cycles, format_equation, list_species,
)
from kinetics.data.experimental import (
    get_barrier_windows, get_measured_shifts, get_water_series, load_experimental_data,
)
from kinetics.data.lab_nmr import (
    DEFAULT_LAB_NMR_DIR, build_lab_nmr_inventory, calculate_area_shares, calculate_window_integrals, estimate_noise,
    find_heated_windows, load_spectrum, load_text_spectrum, read_jdf, refine_window_phase, tabulate_area_shares,
)

__all__ = [
    "load_snapshot", "get_dataset_records", "get_solvation_records", "count_elements", "calculate_molar_mass",
    "load_species_database",
    "NETWORK", "NETWORK_SPECIES", "build_stoichiometric_matrix", "find_reaction_cycles", "format_equation",
    "list_species",
    "load_experimental_data", "get_measured_shifts", "get_barrier_windows", "get_water_series",
    "DEFAULT_LAB_NMR_DIR", "build_lab_nmr_inventory", "read_jdf", "load_spectrum", "load_text_spectrum",
    "estimate_noise", "calculate_window_integrals", "refine_window_phase", "calculate_area_shares",
    "tabulate_area_shares", "find_heated_windows",
]
