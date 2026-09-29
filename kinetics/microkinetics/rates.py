"""Block 5 — Rate constants: barrier model + Eyring TST + detailed balance.

    ΔG‡_f = wR + F(ΔG_rxn - wR + wP; g(T)),   g(T) = g(T_ref) - (T - T_ref)·ΔS‡
    k_f   = (kB·T/h) · exp(-ΔG‡_f / RT)        [optionally in series with the diffusion limit k_D]
    k_r   = k_f / K_eq

Kinetic models (the `kinetic_model` argument):
    'bep_eyring'     legacy capped BEP, ΔG‡_f = max(g, g + α·ΔG_rxn); direction-dependent
    'marcus_eyring'  Marcus, ΔG‡_f = g (1 + ΔG_rxn / 4g)²
    'agmon_levine', 'blowers_masel', 'two_parabola'   smooth Level 1 shapes with the family g
    'level1'         per-family shape and LEVEL1_PARAMETERS
Each model has a fallback parameter set used when `family_params` is None or lacks a family.
Reactions flagged 'canonical': False in the network use the reversed family parameters.
"""

import numpy as np
import pandas as pd

from kinetics.constants import EV_TO_KJ_MOL, H_SI, KB_EV, KB_SI, R_SI
from kinetics.microkinetics.barriers import calculate_barrier, reversed_params
from kinetics.microkinetics.parameters import (
    FAMILY_BEP_PARAMETERS, LEVEL1_PARAMETERS, resolve_family_params,
)
from kinetics.data.network import NETWORK
from kinetics.thermo.reaction import calculate_reaction_thermo

# kinetic_model → (barrier shape, or None to read each family's 'shape'; fallback parameter set)
KINETIC_MODELS = {
    'bep_eyring':    ('bep_cap', FAMILY_BEP_PARAMETERS),
    'marcus_eyring': ('marcus', FAMILY_BEP_PARAMETERS),
    'agmon_levine':  ('agmon_levine', FAMILY_BEP_PARAMETERS),
    'blowers_masel': ('blowers_masel', FAMILY_BEP_PARAMETERS),
    'two_parabola':  ('two_parabola', FAMILY_BEP_PARAMETERS),
    'level1':        (None, LEVEL1_PARAMETERS),
}


def calculate_rate_constants(rxn_id: str, T_K: float, *, kinetic_model: str = 'bep_eyring',
                             family_params: dict = None, network: dict = None, species_db: dict = None,
                             thermo_mode: str = 'wb97mv', viscosity_Pa_s: float = None) -> dict:
    """k_f, k_r and the barriers of one reaction at T_K.

    viscosity_Pa_s enables the Collins-Kimball ceiling k_D = 8RT/3η for bimolecular steps.
    """
    if kinetic_model not in KINETIC_MODELS:
        raise ValueError(f"Unknown kinetic_model '{kinetic_model}'. Supported: {list(KINETIC_MODELS)}")
    shape, fallback_set = KINETIC_MODELS[kinetic_model]
    net = network or NETWORK
    rxn = net[rxn_id]

    thermo = calculate_reaction_thermo(rxn_id, T_K, network=net, species_db=species_db, thermo_mode=thermo_mode)
    fam = resolve_family_params(thermo['class'], family_params if family_params is not None else fallback_set,
                                fallback_set)
    active_shape = shape or fam.get('shape', 'marcus')
    g_eV = fam['g_eV'] - (T_K - fam['T_ref_K']) * fam['dS_act_J_mol_K'] / 1000.0 / EV_TO_KJ_MOL
    if not rxn.get('canonical', True):
        fam = reversed_params(fam)

    shape_params = {'wR_eV': fam['wR_eV'], 'wP_eV': fam['wP_eV']}
    if active_shape == 'bep_cap':
        shape_params['alpha'] = fam['alpha']
    elif active_shape == 'blowers_masel':
        shape_params['w'] = fam['w_eV']
    elif active_shape == 'two_parabola':
        shape_params['alpha0'] = fam['alpha0']

    dG_rxn_kJ_mol = thermo['dG_rxn_kJ_mol']
    dG_barrier_f_eV, alpha_eff = calculate_barrier(thermo['dG_rxn_eV'], active_shape, g_eV, **shape_params)
    dG_barrier_f_kJ_mol = dG_barrier_f_eV * EV_TO_KJ_MOL

    eyring_prefactor = KB_SI * T_K / H_SI
    k_f = eyring_prefactor * np.exp(-(dG_barrier_f_kJ_mol * 1000.0) / (R_SI * T_K))

    # Collins-Kimball: bimolecular steps cannot exceed the Smoluchowski limit k_D = 8RT / 3η
    k_D = None
    if viscosity_Pa_s is not None and sum(rxn['reactants'].values()) == 2:
        k_D = 8.0 * R_SI * T_K / (3.0 * viscosity_Pa_s) * 1000.0   # m³/(mol·s) → M⁻¹·s⁻¹
        k_f = 1.0 / (1.0 / k_f + 1.0 / k_D)
        dG_barrier_f_kJ_mol = -R_SI * T_K * np.log(k_f / eyring_prefactor) / 1000.0
        dG_barrier_f_eV = dG_barrier_f_kJ_mol / EV_TO_KJ_MOL

    # Detailed balance, applied after any ceiling
    K_eq = thermo['K_eq']
    k_r = k_f / K_eq if K_eq > 1e-300 else 0.0
    dG_barrier_r_kJ_mol = dG_barrier_f_kJ_mol - dG_rxn_kJ_mol

    result = {
        'rxn_id': rxn_id,
        'T_K': T_K,
        'kinetic_model': kinetic_model,
        'barrier_model': active_shape,
        'thermo_mode': thermo_mode,
        'class': thermo['class'],
        'g_eV': g_eV,
        'alpha_eff': alpha_eff,
        'dG_rxn_eV': thermo['dG_rxn_eV'],
        'dG_rxn_kJ_mol': dG_rxn_kJ_mol,
        'dG_barrier_f_kJ_mol': dG_barrier_f_kJ_mol,
        'dG_barrier_f_eV': dG_barrier_f_eV,
        'dG_barrier_r_kJ_mol': dG_barrier_r_kJ_mol,
        'dG_barrier_r_eV': dG_barrier_r_kJ_mol / EV_TO_KJ_MOL,
        'k_f': k_f,
        'k_r': k_r,
        'K_eq': K_eq,
    }
    if active_shape == 'marcus':
        result['lambda_eV'] = 4.0 * g_eV
    if k_D is not None:
        result['k_D'] = k_D
    return result


def calculate_eyring_rate(dG_barrier_eV: float, T_K: float) -> float:
    """Eyring rate constant k = (kB·T/h) exp(-ΔG‡/kB·T) [s⁻¹ or M⁻¹s⁻¹] for a barrier in eV."""
    return KB_SI * T_K / H_SI * np.exp(-dG_barrier_eV / (KB_EV * T_K))


def calculate_network_rates(T_K: float, *, network: dict = None, **rate_options) -> pd.DataFrame:
    """calculate_rate_constants for every reaction, one row per reaction (network order).

    rate_options: kinetic_model, family_params, species_db, thermo_mode, viscosity_Pa_s.
    """
    net = network or NETWORK
    columns = ['class', 'barrier_model', 'g_eV', 'dG_rxn_eV', 'dG_barrier_f_eV', 'dG_barrier_r_eV',
               'alpha_eff', 'k_f', 'k_r', 'K_eq']
    rows = {rxn_id: calculate_rate_constants(rxn_id, T_K, network=net, **rate_options) for rxn_id in net}
    return pd.DataFrame.from_dict(rows, orient='index')[columns]
