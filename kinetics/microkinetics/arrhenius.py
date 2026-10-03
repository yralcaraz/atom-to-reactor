"""Block 6 — Modified Arrhenius regression of k(T).

    ln k = ln A + β ln(T/T0) - Ea/(RT),   fitted by least squares.

With a temperature-independent barrier the Eyring rate already has this form (A = kB·T0/h, β = 1,
Ea = ΔG‡), so the fit returns its inputs with R² = 1 (Finding 9). It becomes informative with g(T),
a diffusion ceiling, or measured k(T).

Classification: PUBLIC (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván
"""

import numpy as np
import pandas as pd

from kinetics.constants import EV_TO_KJ_MOL, R_SI, T_STD_K
from kinetics.microkinetics.rates import calculate_rate_constants
from kinetics.data.network import NETWORK


def fit_modified_arrhenius(rxn_id: str, T_grid_K: np.ndarray, *, direction: str = 'f', T0_K: float = T_STD_K,
                           **rate_options) -> dict:
    """A, β, Ea and R² of k_f (direction='f') or k_r ('r') over T_grid_K.

    rate_options are passed to calculate_rate_constants (kinetic_model, family_params, network, ...).
    """
    T = np.asarray(T_grid_K, dtype=float)
    key = {'f': 'k_f', 'r': 'k_r'}[direction]
    k = np.array([calculate_rate_constants(rxn_id, T_K, **rate_options)[key] for T_K in T])

    y = np.log(k)
    X = np.column_stack([np.ones_like(T), np.log(T / T0_K), -1.0 / (R_SI * T)])
    params, *_ = np.linalg.lstsq(X, y, rcond=None)
    ln_A, beta, Ea_J_mol = params
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1.0 - np.sum((y - X @ params) ** 2) / ss_tot if ss_tot > 0 else 1.0
    return {
        'rxn_id': rxn_id,
        'direction': direction,
        'A': np.exp(ln_A),
        'beta': beta,
        'Ea_kJ_mol': Ea_J_mol / 1000.0,
        'Ea_eV': Ea_J_mol / 1000.0 / EV_TO_KJ_MOL,
        'R2': r2,
    }


def build_arrhenius_table(T_grid_K: np.ndarray, *, network: dict = None, **rate_options) -> pd.DataFrame:
    """Forward and reverse modified-Arrhenius parameters of every reaction (numeric columns)."""
    net = network or NETWORK
    rows = {}
    for rxn_id, rxn in net.items():
        row = {'class': rxn['class']}
        for direction in ('f', 'r'):
            fit = fit_modified_arrhenius(rxn_id, T_grid_K, direction=direction, network=net, **rate_options)
            row.update({f'A_{direction}': fit['A'], f'beta_{direction}': fit['beta'],
                        f'Ea_{direction}_eV': fit['Ea_eV'], f'R2_{direction}': fit['R2']})
        rows[rxn_id] = row
    return pd.DataFrame.from_dict(rows, orient='index')
