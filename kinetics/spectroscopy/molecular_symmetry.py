"""Molecular Symmetry and Equivalence Class Partitioning for Multi-nuclear NMR.

Computes topological equivalence classes from 3D Cartesian coordinates using
covalent bond connectivity graphs and Weisfeiler-Lehman multiset refinement.
Averages DFT GIAO isotropic magnetic shieldings across equivalent nuclei,
identifies exchangeable labile protons (-OH, -NH), and computes chemical shifts
referenced to experimental secondary standards (TMS, H3PO4).
"""

import os
import json
import numpy as np
import pandas as pd

DEFAULT_COVALENT_RADII = {
    "H": 0.31, "C": 0.76, "O": 0.66, "Si": 1.11, "P": 1.07,
    "F": 0.57, "S": 1.05, "Li": 1.28, "N": 0.71, "Cl": 1.02
}

DEFAULT_NUCLEI = {
    "Si": ("TMS",   "29Si"),
    "P":  ("H3PO4", "31P"),
    "C":  ("TMS",   "13C"),
    "H":  ("TMS",   "1H"),
}

DEFAULT_SPECIES = [
    "TMSPA", "BMSPA", "MMSPA", "H3PO4", "H2O", "TMSOH", "HMDSO",
    "EC", "TMSOEG", "TMSOdiEG", "CO2"
]


def load_dataset_snapshot(snapshot_path: str = None) -> dict:
    """Loads dataset API snapshot JSON containing datasets, structures, and NMR shieldings."""
    if snapshot_path is None:
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        candidate_paths = [
            os.path.join(base_dir, "data", "tank_api_snapshot.json"),
            os.path.join(base_dir, "tank_api_snapshot.json"),
        ]
        for p in candidate_paths:
            if os.path.exists(p):
                snapshot_path = p
                break
        if snapshot_path is None:
            raise FileNotFoundError("Could not locate tank_api_snapshot.json in data/ or root directory.")

    with open(snapshot_path, "r", encoding="utf-8") as f:
        return json.load(f)


def build_bond_graph(elements: list[str], coordinates_angstrom: list[list[float]] | np.ndarray,
                     covalent_radii: dict = None, scale_factor: float = 1.2) -> np.ndarray:
    """Constructs covalent bond adjacency matrix from Cartesian coordinates."""
    radii_map = covalent_radii or DEFAULT_COVALENT_RADII
    xyz = np.asarray(coordinates_angstrom, dtype=float)
    n_atoms = len(elements)
    
    dist_mat = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    r_arr = np.array([radii_map.get(e, 0.70) for e in elements], dtype=float)
    r_cutoff = scale_factor * (r_arr[:, None] + r_arr[None, :])
    
    bonded = (dist_mat < r_cutoff) & (~np.eye(n_atoms, dtype=bool))
    return bonded


def find_equivalence_classes(elements: list[str], bond_matrix: np.ndarray) -> list[str]:
    """Partitions atoms into topological equivalence classes via 1-WL multiset refinement."""
    n = len(elements)
    labels = list(elements)
    
    for _ in range(n):
        keys = [
            (labels[i], tuple(sorted(labels[j] for j in np.flatnonzero(bond_matrix[i]))))
            for i in range(n)
        ]
        unique_keys = sorted(set(keys))
        lookup = {k: f"class_{idx}" for idx, k in enumerate(unique_keys)}
        new_labels = [lookup[k] for k in keys]
        
        if len(set(new_labels)) == len(set(labels)):
            break
        labels = new_labels
        
    return labels


def extract_nmr_sites(pipeline_id: str, element: str, snapshot_data: dict = None) -> list[tuple[float, int, bool]] | None:
    """Extracts symmetry-averaged isotropic magnetic shieldings for a target element in a molecule."""
    data = snapshot_data if snapshot_data is not None else load_dataset_snapshot()
    
    datasets = data.get("datasets", [])
    matched_entry = None
    for entry in datasets:
        if entry.get("pipeline_id") == pipeline_id and entry.get("has_nmr", False):
            matched_entry = entry
            break
            
    if matched_entry is None:
        return None
        
    uuid = matched_entry["dataset_uuid"]
    nmr_shieldings = data.get("nmr", {}).get(uuid, {}).get("nmr_shieldings", [])
    structure = data.get("structure", {}).get(uuid, {})
    
    if not nmr_shieldings or not structure:
        return None
        
    elements = structure["elements"]
    coords = structure["coordinates_angstrom"]
    
    bonded = build_bond_graph(elements, coords)
    classes = find_equivalence_classes(elements, bonded)
    
    grouped_shieldings = {}
    for entry in nmr_shieldings:
        if entry.get("element") != element:
            continue
        atom_idx = entry["atom_index"]
        cls_id = classes[atom_idx]
        
        is_labile = False
        if element == "H":
            neighbor_elements = [elements[j] for j in np.flatnonzero(bonded[atom_idx])]
            is_labile = any(nbr in ("O", "N") for nbr in neighbor_elements)
            
        key = (cls_id, is_labile)
        grouped_shieldings.setdefault(key, []).append(entry["isotropic_ppm"])
        
    results = [
        (float(np.mean(vals)), len(vals), is_labile)
        for (_, is_labile), vals in grouped_shieldings.items()
    ]
    return results


def build_referenced_nmr_sites(
    species_list: list[str] = None,
    nuclei_dict: dict = None,
    labile_shift_overrides: dict = None,
    snapshot_data: dict = None
) -> dict:
    """Builds the comprehensive referenced chemical shift catalog for multi-nuclear NMR."""
    species = species_list or DEFAULT_SPECIES
    nuclei = nuclei_dict or DEFAULT_NUCLEI
    overrides = labile_shift_overrides or {}
    snapshot = snapshot_data if snapshot_data is not None else load_dataset_snapshot()
    
    nmr_sites_catalog = {}
    
    for el, (ref_species, isotope_label) in nuclei.items():
        ref_data = extract_nmr_sites(ref_species, el, snapshot_data=snapshot)
        if not ref_data:
            print(f"Warning: Reference species '{ref_species}' for {isotope_label} has no NMR shielding data.")
            continue
            
        sigma_ref = ref_data[0][0]
        nmr_sites_catalog[el] = {}
        
        for sp in species:
            sites = extract_nmr_sites(sp, el, snapshot_data=snapshot)
            if not sites:
                continue
                
            processed_sites = []
            for s_mean, mult, is_lab in sites:
                if is_lab and sp in overrides:
                    delta = float(overrides[sp])
                else:
                    delta = float(sigma_ref - s_mean)
                processed_sites.append((delta, mult, is_lab))
                
            processed_sites.sort(key=lambda item: item[0], reverse=True)
            nmr_sites_catalog[el][sp] = processed_sites
            
    return nmr_sites_catalog
