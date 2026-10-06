"""Experimental reference data: Gogoi et al., J. Phys. Chem. C 2024, 128, 1654.

`data/experimental_gogoi2024.json` holds four kinds of data:
    nmr_shifts           measured ³¹P, ²⁹Si, ¹³C and ¹H shifts, each with its site, peak and figure
    barrier_constraints  ΔG‡ windows derived from the reported observations (R4, R5, R8)
    control_experiments  the paper's simple mixtures, with a window on one observable each
    water_series         ³¹P of TMSPA with 0.5–5 vol% water, and of the 2 vol% sample heated to 80 °C
Entries marked 'assumed' rest on a detection limit or reaction time the paper does not state; 'reading'
marks our quantitative reading of a qualitative description.

Source: Y. Alcaraz Galván; data from Gogoi et al., J. Phys. Chem. C 2024, 128, 1654
"""

import json
import os
from functools import lru_cache

import pandas as pd

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_EXPERIMENTAL_PATH = os.path.join(_REPO_ROOT, 'data', 'experimental_gogoi2024.json')


@lru_cache(maxsize=4)
def load_experimental_data(path: str = DEFAULT_EXPERIMENTAL_PATH) -> dict:
    """Parsed reference file (cached; do not modify)."""
    with open(path, 'r', encoding='utf-8') as f:
        return json.load(f)


def get_measured_shifts(nucleus: str = None) -> pd.DataFrame:
    """Measured shifts [ppm]: one row per (nucleus, species), with 'mid_ppm' the centre of the reported range."""
    df = pd.DataFrame(load_experimental_data()['nmr_shifts'])
    df['mid_ppm'] = 0.5 * (df['low_ppm'] + df['high_ppm'])
    return df[df['nucleus'] == nucleus].reset_index(drop=True) if nucleus else df


def get_barrier_windows() -> pd.DataFrame:
    """ΔG‡ windows per reaction [eV] (NaN = open side), with the temperature and basis of each."""
    return pd.DataFrame(load_experimental_data()['barrier_constraints']).set_index('rxn_id').astype(
        {'low_eV': float, 'high_eV': float})


def get_water_series() -> dict:
    """Water series of TMSPA in EC/DEC: 'medium', 'samples' (windows on P fractions) and 'heating'."""
    return load_experimental_data()['water_series']
