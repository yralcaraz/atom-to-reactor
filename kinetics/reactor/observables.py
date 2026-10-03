"""Observables extracted from reactor results: characteristic times and composition at acquisitions.

Classification: PUBLIC (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván
"""

import numpy as np
import pandas as pd

from kinetics.constants import R_SI


def find_crossing_time(t: np.ndarray, y: np.ndarray, level: float) -> float:
    """First time y falls to `level` (linear interpolation); NaN if it never does."""
    below = np.nonzero(y <= level)[0]
    if len(below) == 0:
        return np.nan
    i = below[0]
    if i == 0:
        return float(t[0])
    return float(t[i - 1] + (level - y[i - 1]) * (t[i] - t[i - 1]) / (y[i] - y[i - 1]))


def calculate_worst_case_pressure_bar(C_gas_M: float, liquid_mL: float, gas_mL: float, T_K: float) -> float:
    """Headspace pressure [bar] if all of a dissolved gas left the liquid (ideal gas); an upper bound for safety.

    p = C·V_liquid·R·T / V_gas. Henry partitioning keeps part of the gas dissolved, so the real pressure is lower.
    """
    return C_gas_M * liquid_mL * R_SI * T_K / gas_mL / 100.0     # mol/L·mL = mmol; Pa → bar


def summarize_batch_runs(sims: dict) -> pd.DataFrame:
    """Scorecard per run: TMSPA half-life, 90 % water consumption, peak TMSOH, final CO2 and HMDSO."""
    rows = {}
    for name, sim in sims.items():
        C, idx, t_h = sim['C_mM'], sim['idx'], sim['t_h']
        rows[name] = {
            't½ TMSPA (h)': find_crossing_time(t_h, C[idx['TMSPA']], 0.5 * C[idx['TMSPA'], 0]),
            't 90% H2O used (h)': find_crossing_time(t_h, C[idx['H2O']], 0.1 * C[idx['H2O'], 0]),
            'peak TMSOH (mM)': C[idx['TMSOH']].max(),
            'final CO2 (mM)': C[idx['CO2'], -1],
            'final HMDSO (mM)': C[idx['HMDSO'], -1],
        }
    return pd.DataFrame.from_dict(rows, orient='index')


def tabulate_acquisitions(sim: dict, species: list = None) -> pd.DataFrame:
    """Composition [mM] at each NMR acquisition, indexed by the preceding hold temperature (°C)."""
    species = species or [sp for sp in sim['species'] if sp != 'EC']
    return pd.DataFrame(
        {sp: [a['C_M'][sim['idx'][sp]] * 1000.0 for a in sim['acquisitions']] for sp in species},
        index=pd.Index([a['T_C'] for a in sim['acquisitions']], name='hold T (°C)'),
    )


def calculate_remaining_fraction(sim: dict, species: str = 'TMSPA') -> pd.Series:
    """Fraction of `species` left at each acquisition, relative to the injected amount."""
    c0 = sim['C_M'][sim['idx'][species], sim['injection_idx']]
    return tabulate_acquisitions(sim, [species])[species] / (c0 * 1000.0)


def tabulate_trajectory(sim: dict) -> pd.DataFrame:
    """Concentrations [mM] of every species over time, with the temperature [°C], for export."""
    df = pd.DataFrame(sim['C_mM'].T, columns=[f'{sp} (mM)' for sp in sim['species']])
    df.insert(0, 'T (°C)', sim['T_K'] - 273.15)
    df.index = pd.Index(sim['t_h'], name='t (h)')
    return df
