"""Synthetic share tables in the format of notebook 03's exports, for tests that must not need the lab data.

The numbers are invented: they have the layout of the lab set (same samples, same number of spectra, same
kind of history) but none of its values.

Classification: PUBLIC (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván
"""

import numpy as np
import pandas as pd

from kinetics.data.observables import LAB_SAMPLES, PHOSPHATE_WINDOWS, SILYL_WINDOWS, build_lab_observables

# sample → [(hours since the first spectrum, temperature °C, scans, shares in window order)]
SYNTHETIC_SPECTRA = {
    '0.5 % H2O': [(19.0, 22.0, 16, (0.90, 0.08, 0.01, 0.01))],
    '2 % H2O': [(0.1, 22.5, 16, (0.00, 0.05, 0.30, 0.65)), (33.0, 22.5, 256, (0.00, 0.05, 0.30, 0.65)),
                (72.0, 22.5, 256, (0.00, 0.05, 0.30, 0.65)), (96.0, 22.5, 256, (0.00, 0.05, 0.30, 0.65)),
                (120.0, 22.5, 256, (0.00, 0.05, 0.30, 0.65)), (145.0, 22.5, 256, (0.00, 0.05, 0.30, 0.65))],
    'TMSPa + TMSOH (A)': [(0.1, 22.0, 16, (0.10, 0.70, 0.19, 0.01))],
    'TMSPa + TMSOH (B)': [(2.0, 22.0, 16, (0.10, 0.85, 0.04, 0.01))],
    'TMSPa alone': [(0.05, 22.0, 16, (0.80, 0.15, 0.03, 0.02)), (240.0, 22.0, 16, (0.25, 0.70, 0.03, 0.02))],
    'TMSOH, probe': [(0.2, 23.0, 64, (0.97, 0.02, 0.01)), (138.7, 80.0, 64, (0.95, 0.03, 0.02)),
                     (144.0, 22.0, 64, (0.95, 0.03, 0.02))],
    'TMSOH, glovebox': [(0.1, 22.0, 64, (0.45, 0.15, 0.40))],
}
SYNTHETIC_SIGMA = {'31P': {16: (0.010, 0.012), 256: (0.003, 0.012)}, '13C': {64: (0.015, 0.015)}}


def build_synthetic_share_tables(spectra: dict = None) -> tuple:
    """(³¹P table, ¹³C table) with one row per (spectrum, window), like `tabulate_area_shares` + notebook 03."""
    spectra = spectra if spectra is not None else SYNTHETIC_SPECTRA
    rows = {'31P': [], '13C': []}
    t0 = pd.Timestamp('2020-01-01 12:00')
    for k, (sample, records) in enumerate(spectra.items()):
        nucleus = LAB_SAMPLES[sample]['nucleus']
        windows = PHOSPHATE_WINDOWS if nucleus == '31P' else SILYL_WINDOWS
        for j, (hours, T_C, scans, shares) in enumerate(records):
            noise, baseline = SYNTHETIC_SIGMA[nucleus][scans]
            for window, share in zip(windows, shares):
                rows[nucleus].append({
                    'sample': sample, 'file': f'synthetic_{k}_{j}.jdf', 'nucleus': nucleus,
                    'acquired_at': t0 + pd.Timedelta(days=30 * k, hours=hours), 'hours_since_first': hours,
                    'temperature_C': T_C, 'scans': scans, 'relaxation_delay_s': 5.0 if nucleus == '31P' else 2.0,
                    'pulse_angle_deg': 45.0 if nucleus == '31P' else 30.0, 'noe': nucleus == '13C',
                    'lock_solvent': 'DMSO-D6', 'data_md5': f'{k:02d}{j:02d}', 'window': window, 'share': share,
                    'sigma_noise': noise, 'sigma_baseline': baseline, 'sigma': float(np.hypot(noise, baseline)),
                })
    return pd.DataFrame(rows['31P']), pd.DataFrame(rows['13C'])


def build_synthetic_observables(spectra: dict = None) -> dict:
    """Observables dict with the layout of the lab set and invented shares."""
    shares_P, shares_C = build_synthetic_share_tables(spectra)
    return build_lab_observables(shares_P, shares_C, built_from='synthetic numbers for tests (no lab data)')
