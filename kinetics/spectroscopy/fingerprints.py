"""Block 13 — Reaction fingerprints, identifiability and extent recovery.

The fingerprint of reaction r is F_r = Σ_i ν_ir P_i, with P_i the pure spectrum of species i (per mM)
concatenated over the nuclei. Each nucleus block is scaled to unit norm. A greedy SVD-rank basis gives
the identifiable (lumped) reactions; the Net Analyte Signal (NAS) gives, per nucleus, the share of each
fingerprint no other reaction can mimic; the pseudoinverse recovers per-step extents Δξ from spectra.

Source: Y. Alcaraz Galván; adapted from P. Broqvist
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd

from kinetics.data.network import NETWORK, NETWORK_SPECIES, build_stoichiometric_matrix
from kinetics.spectroscopy.symmetry import DEFAULT_NUCLEI, build_nmr_catalog
from kinetics.spectroscopy.spectra import LINE_SHAPES, find_shift_windows, lorentzian


@dataclass
class FeatureSpace:
    """Concatenated multinuclear ppm axis: per element its windows, full axis and slice in the feature vector."""
    windows: dict       # element → [(high_ppm, low_ppm, x_ppm)]
    x_ppm: dict         # element → concatenated axis
    blocks: dict        # element → slice
    n_features: int


def find_visible_species(species: list = None, nmr_catalog: dict = None) -> list:
    """Species with at least one non-labile site (H2O, whose only protons exchange, is invisible)."""
    catalog = nmr_catalog if nmr_catalog is not None else build_nmr_catalog()
    return [sp for sp in (species or NETWORK_SPECIES) if sp != 'EC'
            and any(not labile for el in catalog for _, _, labile in catalog[el].get(sp, []))]


def build_feature_space(visible_species: list, nmr_catalog: dict = None, nuclei=('Si', 'P', 'C', 'H'),
                        line_shapes: dict = None) -> FeatureSpace:
    catalog = nmr_catalog if nmr_catalog is not None else build_nmr_catalog()
    shapes = line_shapes or LINE_SHAPES
    windows = {}
    for el in nuclei:
        if el not in catalog:
            continue
        shifts = [d for sp in visible_species for d, _, labile in catalog[el].get(sp, []) if not labile]
        if not shifts:
            continue
        windows[el] = []
        for hi, lo in find_shift_windows(shifts, shapes[el]['pad_ppm'], shapes[el]['max_gap_ppm']):
            n_pts = int(np.clip(8.0 * (hi - lo) / shapes[el]['fwhm_ppm'], 600, 6000))
            windows[el].append((hi, lo, np.linspace(hi, lo, n_pts)))
    x_ppm = {el: np.concatenate([x for *_, x in w]) for el, w in windows.items()}
    blocks, start = {}, 0
    for el, x in x_ppm.items():
        blocks[el] = slice(start, start + len(x))
        start += len(x)
    return FeatureSpace(windows, x_ppm, blocks, start)


def build_pure_spectra(visible_species: list, feature_space: FeatureSpace, nmr_catalog: dict = None,
                       line_shapes: dict = None) -> np.ndarray:
    """Pure-species spectra per mM, shape (n_species, n_features)."""
    catalog = nmr_catalog if nmr_catalog is not None else build_nmr_catalog()
    shapes = line_shapes or LINE_SHAPES
    P = np.zeros((len(visible_species), feature_space.n_features))
    for k, sp in enumerate(visible_species):
        for el, sl in feature_space.blocks.items():
            for d, n, labile in catalog.get(el, {}).get(sp, []):
                if not labile:
                    P[k, sl] += n * lorentzian(feature_space.x_ppm[el], d, shapes[el]['fwhm_ppm'])
    return P


def build_reaction_fingerprints(network: dict, visible_species: list, pure_spectra: np.ndarray,
                                feature_space: FeatureSpace) -> dict:
    """Raw and block-normalized fingerprints: {'f_raw', 'f_norm', 'weights', 'S_visible'}."""
    S_visible = build_stoichiometric_matrix(network, visible_species).T
    f_raw = S_visible @ pure_spectra
    weights = np.ones(feature_space.n_features)
    for sl in feature_space.blocks.values():
        norm = np.linalg.norm(f_raw[:, sl])
        if norm > 1e-12:
            weights[sl] = 1.0 / norm
    return {'f_raw': f_raw, 'f_norm': f_raw * weights, 'weights': weights, 'S_visible': S_visible}


def analyze_reaction_identifiability(fingerprints: dict, rxn_labels: list, feature_space: FeatureSpace,
                                     nuclei: dict = None) -> dict:
    """Greedy rank basis, lumped labels (e.g. 'R1 ⊕ R5'), projection and NAS selectivity per nucleus."""
    nuclei = nuclei or DEFAULT_NUCLEI
    f_norm, f_raw = fingerprints['f_norm'], fingerprints['f_raw']

    def rank(m, tol=1e-8):
        s = np.linalg.svd(m, compute_uv=False)
        return int((s > tol * s[0]).sum())

    basis = []
    for i in range(len(f_norm)):
        if rank(f_norm[basis + [i]]) > len(basis):
            basis.append(i)
    fb_norm = f_norm[basis]
    projection = np.linalg.lstsq(fb_norm.T, f_norm.T, rcond=None)[0].T
    projection[np.abs(projection) < 1e-8] = 0.0
    projection = np.round(projection, 8)

    lumped = []
    for j, b in enumerate(basis):
        members = [rxn_labels[r] for r in range(len(f_norm)) if r != b and abs(projection[r, j]) > 1e-6]
        lumped.append(f"{rxn_labels[b]} ⊕ {','.join(members)}" if members else rxn_labels[b])

    def net_analyte_signal(f, j):
        others = np.delete(f, j, axis=0)
        return f[j] - others.T @ np.linalg.lstsq(others.T, f[j], rcond=None)[0]

    columns = list(feature_space.blocks) + ['all']
    selectivity = np.full((len(basis), len(columns)), np.nan)
    for c, col in enumerate(columns):
        fb = fb_norm if col == 'all' else fb_norm[:, feature_space.blocks[col]]
        scale = np.abs(fb).max()
        for j in range(len(basis)):
            norm_j = np.linalg.norm(fb[j])
            if norm_j > 1e-6 * scale:
                selectivity[j, c] = np.linalg.norm(net_analyte_signal(fb, j)) / norm_j
    names = [nuclei[el][1] if el != 'all' else 'all nuclei' for el in columns]
    return {
        'basis_indices': basis,
        'lumped_labels': lumped,
        'projection': projection,
        'F_basis_norm': fb_norm,
        'F_basis_raw': f_raw[basis],
        'selectivity': pd.DataFrame(selectivity, index=lumped, columns=names),
        'effective_rank': len(basis),
    }


def recover_reaction_extents(dD: np.ndarray, f_basis_norm: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Per-step extents Δξ [mM] of the basis reactions from spectral differences, by pseudoinverse."""
    return (dD * weights) @ np.linalg.pinv(f_basis_norm)


def run_fingerprint_analysis(network: dict = None, species: list = None, nmr_catalog: dict = None) -> dict:
    """Feature space, pure spectra, fingerprints and identifiability of a network in one call."""
    net = network or NETWORK
    catalog = nmr_catalog if nmr_catalog is not None else build_nmr_catalog()
    visible = find_visible_species(species, catalog)
    space = build_feature_space(visible, catalog)
    pure = build_pure_spectra(visible, space, catalog)
    fingerprints = build_reaction_fingerprints(net, visible, pure, space)
    identifiability = analyze_reaction_identifiability(fingerprints, list(net), space)
    return {'visible_species': visible, 'feature_space': space, 'pure_spectra': pure,
            **fingerprints, **identifiability}


def simulate_acquisition_spectra(sim: dict, analysis: dict, noise_rel: float = 0.0, seed: int = 42) -> np.ndarray:
    """Concatenated spectra at injection and at every acquisition of a protocol run, with optional noise.

    Noise is Gaussian with σ = noise_rel × the block maximum, independently per nucleus.
    """
    steps = [sim['injection_idx']] + [a['idx'] for a in sim['acquisitions']]
    rows = [sim['idx'][sp] for sp in analysis['visible_species']]
    D = (sim['C_M'][rows][:, steps].T * 1000.0) @ analysis['pure_spectra']
    if noise_rel > 0:
        rng = np.random.default_rng(seed)
        for sl in analysis['feature_space'].blocks.values():
            D[:, sl] += rng.normal(0.0, noise_rel * np.abs(D[:, sl]).max(), D[:, sl].shape)
    return D
