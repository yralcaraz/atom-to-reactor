"""Blocks 8, 11–13 — Spectroscopy: NMR shifts, synthetic spectra, water balance and reaction fingerprints.

    symmetry.py      chemical shifts from DFT shieldings, averaged over symmetry-equivalent atoms
    spectra.py       Lorentzian spectra from concentrations, OH exchange, water mass balance
    fingerprints.py  reaction fingerprints, identifiability and recovery of reaction extents

Classification: PUBLIC (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván
"""

from kinetics.spectroscopy.symmetry import (
    DEFAULT_NUCLEI, build_bond_graph, build_nmr_catalog, extract_shielding_sites, find_equivalence_classes,
)
from kinetics.spectroscopy.spectra import (
    EXCHANGE_PEAK, LINE_SHAPES, calculate_nmr_peaks, calculate_water_mass_balance, find_shift_windows,
    lorentzian, simulate_nmr_spectra,
)
from kinetics.spectroscopy.fingerprints import (
    FeatureSpace, analyze_reaction_identifiability, build_feature_space, build_pure_spectra,
    build_reaction_fingerprints, find_visible_species, recover_reaction_extents, run_fingerprint_analysis,
    simulate_acquisition_spectra,
)

__all__ = [
    "DEFAULT_NUCLEI",
    "build_bond_graph",
    "build_nmr_catalog",
    "extract_shielding_sites",
    "find_equivalence_classes",
    "EXCHANGE_PEAK",
    "LINE_SHAPES",
    "calculate_nmr_peaks",
    "calculate_water_mass_balance",
    "find_shift_windows",
    "lorentzian",
    "simulate_nmr_spectra",
    "FeatureSpace",
    "analyze_reaction_identifiability",
    "build_feature_space",
    "build_pure_spectra",
    "build_reaction_fingerprints",
    "find_visible_species",
    "recover_reaction_extents",
    "run_fingerprint_analysis",
    "simulate_acquisition_spectra",
]
