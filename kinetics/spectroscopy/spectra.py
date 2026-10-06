"""Block 8 / 12 — Synthetic multinuclear NMR (29Si, 31P, 13C, 1H) from concentrations.

Every resolved site contributes a Lorentzian of area ∝ n_atoms · C. Labile OH protons exchange faster
than the NMR time scale and coalesce into one population-weighted peak. The solvent (EC) is not summed
into the spectra. Intensities are in mM of nuclei.

Source: Y. Alcaraz Galván; adapted from P. Broqvist
"""

import numpy as np
import pandas as pd

from kinetics.spectroscopy.symmetry import build_nmr_catalog

# Default linewidths (FWHM) and window clustering per nucleus [ppm]
LINE_SHAPES = {
    'Si': {'fwhm_ppm': 0.4, 'pad_ppm': 4.0, 'max_gap_ppm': 40.0},
    'P':  {'fwhm_ppm': 0.5, 'pad_ppm': 4.0, 'max_gap_ppm': 40.0},
    'C':  {'fwhm_ppm': 0.4, 'pad_ppm': 4.0, 'max_gap_ppm': 40.0},
    'H':  {'fwhm_ppm': 0.02, 'pad_ppm': 0.4, 'max_gap_ppm': 2.0},
}
EXCHANGE_PEAK = 'OH (exch.)'


def lorentzian(x, x0, fwhm: float):
    """Lorentzian of unit height centred at x0."""
    gamma = fwhm / 2.0
    return gamma**2 / ((x - x0)**2 + gamma**2)


def calculate_nmr_peaks(element: str, C_M: np.ndarray, species_idx: dict, nmr_catalog: dict = None,
                        solvent=('EC',)) -> list:
    """[(label, shift_ppm[t], intensity_mM[t])] for every resolved site plus the exchange-averaged OH peak.

    C_M has one column per time point.
    """
    catalog = nmr_catalog if nmr_catalog is not None else build_nmr_catalog()
    n_t = C_M.shape[1]
    peaks = []
    oh_weight = np.zeros(n_t)
    oh_weighted_shift = np.zeros(n_t)
    for sp, sites in catalog.get(element, {}).items():
        if sp not in species_idx or sp in solvent:
            continue
        conc_mM = C_M[species_idx[sp]] * 1000.0
        for shift, n_atoms, labile in sites:
            intensity = conc_mM * n_atoms
            if labile:
                oh_weight += intensity
                oh_weighted_shift += intensity * shift
            else:
                peaks.append((sp, np.full(n_t, shift), intensity))
    if np.any(oh_weight > 0):
        peaks.append((EXCHANGE_PEAK, oh_weighted_shift / np.where(oh_weight > 0, oh_weight, 1.0), oh_weight))
    return peaks


def simulate_nmr_spectra(element: str, C_M: np.ndarray, species_idx: dict, x_ppm: np.ndarray, *,
                         fwhm_ppm: float = None, nmr_catalog: dict = None, solvent=('EC',)) -> np.ndarray:
    """Spectra (n_times × len(x_ppm)) as a sum of Lorentzians."""
    fwhm = fwhm_ppm if fwhm_ppm is not None else LINE_SHAPES[element]['fwhm_ppm']
    spectra = np.zeros((C_M.shape[1], len(x_ppm)))
    for _, shift, intensity in calculate_nmr_peaks(element, C_M, species_idx, nmr_catalog, solvent):
        spectra += intensity[:, None] * lorentzian(x_ppm[None, :], shift[:, None], fwhm)
    return spectra


def find_shift_windows(shifts_ppm, pad_ppm: float = 4.0, max_gap_ppm: float = 40.0) -> list:
    """Cluster shifts into disjoint display windows [(high_ppm, low_ppm)], splitting at gaps > max_gap_ppm."""
    if len(shifts_ppm) == 0:
        return [(10.0, -10.0)]
    s = sorted(shifts_ppm, reverse=True)
    windows, high = [], s[0]
    for a, b in zip(s, s[1:]):
        if a - b > max_gap_ppm:
            windows.append((high + pad_ppm, a - pad_ppm))
            high = b
    windows.append((high + pad_ppm, s[-1] - pad_ppm))
    return windows


def calculate_water_mass_balance(C_M: np.ndarray, species_idx: dict, h2o0_M: float,
                                 nmr_catalog: dict = None) -> dict:
    """Water from the OH mass balance, [H2O] = [H2O]0 − ½ Σ_{i≠H2O} n_OH,i C_i, and the OH carriers.

    The exchange-averaged OH integral stays at 2·[H2O]0 (no reaction destroys an OH proton), so water
    must be read from the other species' resolved peaks.
    """
    catalog = nmr_catalog if nmr_catalog is not None else build_nmr_catalog()
    n_t = C_M.shape[1]
    oh_total_M = np.zeros(n_t)
    oh_non_water_M = np.zeros(n_t)
    carriers = {}
    for sp, sites in catalog.get('H', {}).items():
        n_oh = sum(n for _, n, labile in sites if labile)
        if sp not in species_idx or n_oh == 0:
            continue
        oh_M = n_oh * C_M[species_idx[sp]]
        oh_total_M += oh_M
        carriers[f'{sp} ({n_oh} OH)'] = oh_M * 1000.0
        if sp != 'H2O':
            oh_non_water_M += oh_M
    return {
        'h2o_mass_balance_M': h2o0_M - 0.5 * oh_non_water_M,
        'h2o_direct_M': C_M[species_idx['H2O']] if 'H2O' in species_idx else np.zeros(n_t),
        'oh_total_M': oh_total_M,
        'oh_carriers_mM': pd.DataFrame(carriers),
    }
