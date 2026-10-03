"""Block 11 — NMR chemical-shift catalog from DFT shieldings and molecular symmetry.

Atoms are grouped into topological equivalence classes (covalent bond graph + 1-WL refinement); their
isotropic shieldings are averaged and referenced, δ = σ(reference) − σ(sample). Protons bonded to O or N
are flagged as labile (fast exchange). A catalog maps {element: {species: [(δ_ppm, n_atoms, labile)]}}.

Classification: PUBLIC (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván
"""

import warnings
from functools import lru_cache

import numpy as np

from kinetics.data.network import NETWORK_SPECIES
from kinetics.data.snapshot import get_nmr_shieldings, get_structure

DEFAULT_COVALENT_RADII_A = {
    'H': 0.31, 'C': 0.76, 'O': 0.66, 'Si': 1.11, 'P': 1.07,
    'F': 0.57, 'S': 1.05, 'Li': 1.28, 'N': 0.71, 'Cl': 1.02,
}

# element → (reference species, isotope label)
DEFAULT_NUCLEI = {
    'Si': ('TMS', '29Si'),
    'P':  ('H3PO4', '31P'),
    'C':  ('TMS', '13C'),
    'H':  ('TMS', '1H'),
}


def build_bond_graph(elements: list, coordinates_A, covalent_radii_A: dict = None,
                     scale_factor: float = 1.2) -> np.ndarray:
    """Boolean adjacency matrix: bonded if d < scale_factor · (r_i + r_j)."""
    radii = covalent_radii_A or DEFAULT_COVALENT_RADII_A
    xyz = np.asarray(coordinates_A, dtype=float)
    dist = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    r = np.array([radii.get(e, 0.70) for e in elements])
    return (dist < scale_factor * (r[:, None] + r[None, :])) & ~np.eye(len(elements), dtype=bool)


def find_equivalence_classes(elements: list, bonded: np.ndarray) -> list:
    """Class label per atom from 1-WL colour refinement of the bond graph."""
    n = len(elements)
    labels = list(elements)
    for _ in range(n):
        keys = [(labels[i], tuple(sorted(labels[j] for j in np.flatnonzero(bonded[i])))) for i in range(n)]
        lookup = {k: f'class_{i}' for i, k in enumerate(sorted(set(keys)))}
        new_labels = [lookup[k] for k in keys]
        if len(set(new_labels)) == len(set(labels)):
            break
        labels = new_labels
    return labels


def extract_shielding_sites(species: str, element: str) -> list:
    """[(mean isotropic shielding σ [ppm], n_atoms, labile)] per equivalence class; [] without data."""
    shieldings = get_nmr_shieldings(species)
    if not shieldings:
        return []
    elements, coords = get_structure(species)
    bonded = build_bond_graph(elements, coords)
    classes = find_equivalence_classes(elements, bonded)

    grouped = {}
    for entry in shieldings:
        if entry.get('element') != element:
            continue
        i = entry['atom_index']
        labile = element == 'H' and any(elements[j] in ('O', 'N') for j in np.flatnonzero(bonded[i]))
        grouped.setdefault((classes[i], labile), []).append(entry['isotropic_ppm'])
    return [(float(np.mean(v)), len(v), labile) for (_, labile), v in grouped.items()]


def build_nmr_catalog(species: list = None, nuclei: dict = None, labile_shift_ppm: dict = None) -> dict:
    """Referenced shifts {element: {species: [(δ_ppm, n_atoms, labile)]}}, sorted by δ (descending).

    labile_shift_ppm: measured OH shifts in EC ({species: δ}) replacing the gas-phase labile values.
    The result is cached and shared: do not modify it.
    """
    return _cached_catalog(
        tuple(species or NETWORK_SPECIES),
        tuple(sorted((nuclei or DEFAULT_NUCLEI).items())),
        tuple(sorted((labile_shift_ppm or {}).items())),
    )


@lru_cache(maxsize=16)
def _cached_catalog(species: tuple, nuclei: tuple, labile_shift_ppm: tuple) -> dict:
    overrides = dict(labile_shift_ppm)
    catalog = {}
    for element, (reference, isotope) in nuclei:
        ref_sites = extract_shielding_sites(reference, element)
        if not ref_sites:
            warnings.warn(f"Reference '{reference}' has no {isotope} shieldings; {element} skipped")
            continue
        sigma_ref = ref_sites[0][0]
        catalog[element] = {}
        for sp in species:
            sites = [(float(overrides[sp]) if labile and sp in overrides else sigma_ref - sigma, n, labile)
                     for sigma, n, labile in extract_shielding_sites(sp, element)]
            if sites:
                catalog[element][sp] = sorted(sites, key=lambda s: s[0], reverse=True)
    return catalog
