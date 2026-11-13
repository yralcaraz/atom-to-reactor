"""Block 2 — Gas-phase thermochemistry: stored DFT Gibbs energies or Grimme quasi-RRHO.

Thermo modes (the `thermo_mode` argument used across the pipeline):

    'wb97mv'    stored ωB97M-V Gibbs energy at 298.15 K; G does not change with T.
    'wb97mv_T'  stored ωB97M-V enthalpy and Gibbs energy at 298.15 K give the entropy, S = (H − G)/298.15 K, and
                G(T) = H − T·S: G follows the temperature (H and S are kept at their 298.15 K values). Equal to
                'wb97mv' at 298.15 K. Nothing is fitted (notebooks/H_DFT_Implementation.ipynb).
    'qRRHO'   statistical mechanics from frequencies and moments of inertia (Grimme 2012 damping).
              The snapshot has no frequencies, so for its species this falls back to the stored G
              (with a warning) and only the standard-state shift differs from 'wb97mv'.

Source: Y. Alcaraz Galván
"""

import warnings

import numpy as np

from kinetics.constants import (
    AMU_TO_KG, ANGSTROM_TO_M, C_CM_S, EV_TO_KJ_MOL, H_SI, HBAR_SI, KB_SI, NA, P_STD_PA, R_SI, T_STD_K,
)
from kinetics.data.species import resolve_species_database

THERMO_MODES = ('wb97mv', 'wb97mv_T', 'qRRHO')
STORED_MODES = ('wb97mv', 'wb97mv_T')        # read the stored DFT values; no standard-state shift (solution.py)
DEFAULT_THERMO_MODE = 'wb97mv'

_active_thermo_mode = DEFAULT_THERMO_MODE


def check_thermo_mode(thermo_mode: str) -> None:
    if thermo_mode not in THERMO_MODES:
        raise ValueError(f"Unknown thermo_mode '{thermo_mode}'. Supported: {THERMO_MODES}")


def set_active_thermo_mode(thermo_mode: str = None) -> None:
    """The thermo mode every calculation uses when none is passed (None: the default, 'wb97mv')."""
    global _active_thermo_mode
    thermo_mode = thermo_mode or DEFAULT_THERMO_MODE
    check_thermo_mode(thermo_mode)
    _active_thermo_mode = thermo_mode


def get_active_thermo_mode() -> str:
    return _active_thermo_mode


def resolve_thermo_mode(thermo_mode: str = None) -> str:
    """thermo_mode, or the active one when it is None."""
    thermo_mode = thermo_mode or _active_thermo_mode
    check_thermo_mode(thermo_mode)
    return thermo_mode


def calculate_gas_thermo(species: str, T_K: float, *, species_db: dict = None,
                         thermo_mode: str = None, nu0_cm1: float = 100.0) -> dict:
    """H°, S° and G° of one gas-phase species at T_K [kJ/mol, J/(mol·K)].

    thermo_mode None: the active one (set_active_thermo_mode). nu0_cm1 is the quasi-RRHO damping frequency (qRRHO only).
    """
    thermo_mode = resolve_thermo_mode(thermo_mode)
    data = resolve_species_database(species_db)[species]

    has_freqs = len(data.get('frequencies_cm1', ())) > 0
    if thermo_mode == 'qRRHO' and not has_freqs:
        warnings.warn(f"qRRHO requested but '{species}' has no frequencies; using the stored Gibbs energy",
                      stacklevel=2)
    if thermo_mode == 'wb97mv_T':
        return _stored_gibbs(data, T_K)
    if thermo_mode == 'wb97mv' or not has_freqs:
        return _stored_gibbs(data)
    return _quasi_rrho(data, T_K, nu0_cm1)


def _stored_gibbs(data: dict, T_K: float = None) -> dict:
    """Stored RRHO values at 298.15 K. T_K None: G stays at its 298.15 K value ('wb97mv'); else G = H − T_K·S with
    S = (H − G)/298.15 K ('wb97mv_T'). A species without a stored H has S = 0 and its G does not move."""
    G298_kJ_mol = data['G_gas_eV'] * EV_TO_KJ_MOL
    H_kJ_mol = data.get('H_gas_eV', data['G_gas_eV']) * EV_TO_KJ_MOL
    S_J_mol_K = (H_kJ_mol - G298_kJ_mol) * 1000.0 / T_STD_K
    G_kJ_mol = G298_kJ_mol if T_K is None else H_kJ_mol - T_K * S_J_mol_K / 1000.0
    return {
        'H_gas_kJ_mol': H_kJ_mol,
        'S_gas_J_mol_K': S_J_mol_K,
        'G_gas_kJ_mol': G_kJ_mol,
        'G_gas_eV': data['G_gas_eV'] if T_K is None else G_kJ_mol / EV_TO_KJ_MOL,
    }


def _quasi_rrho(data: dict, T_K: float, nu0_cm1: float) -> dict:
    mass_kg = (data['mass_g_mol'] / 1000.0) / NA
    E_scf_J_mol = data['E_scf_eV'] * EV_TO_KJ_MOL * 1000.0
    sigma_rot = data['sigma_rot']

    # Translation: particle in a box at the 1 bar ideal-gas volume (Sackur-Tetrode)
    lambda_th = np.sqrt(H_SI**2 / (2.0 * np.pi * mass_kg * KB_SI * T_K))
    q_trans = (KB_SI * T_K / P_STD_PA) / lambda_th**3
    H_trans = 2.5 * R_SI * T_K
    S_trans = R_SI * (np.log(q_trans) + 2.5)

    # Rotation: rigid rotor (a zero moment marks a linear molecule)
    I_kg_m2 = [I * AMU_TO_KG * ANGSTROM_TO_M**2 for I in data['moments_amu_A2']]
    if any(I == 0.0 for I in I_kg_m2):
        q_rot = (8.0 * np.pi**2 * max(I_kg_m2) * KB_SI * T_K) / (sigma_rot * H_SI**2)
        H_rot = R_SI * T_K
        S_rot = R_SI * (np.log(q_rot) + 1.0)
    else:
        q_rot = (np.sqrt(np.pi) / sigma_rot) * (8.0 * np.pi**2 * KB_SI * T_K / H_SI**2)**1.5 * np.sqrt(np.prod(I_kg_m2))
        H_rot = 1.5 * R_SI * T_K
        S_rot = R_SI * (np.log(q_rot) + 1.5)

    # Vibration: harmonic enthalpy; entropy interpolated to a free rotor below nu0 (Grimme 2012)
    freqs = np.array(data['frequencies_cm1'], dtype=float)
    ZPE_J_mol = np.sum(0.5 * H_SI * C_CM_S * freqs) * NA
    theta_v = H_SI * C_CM_S * freqs / KB_SI
    x = theta_v / T_K
    exp_x = np.exp(np.clip(x, 0, 100))
    H_vib = ZPE_J_mol + R_SI * np.sum(theta_v / (exp_x - 1.0))
    S_vib_ho = R_SI * (x / (exp_x - 1.0) - np.log(1.0 - np.exp(-np.clip(x, 0, 100))))
    I_eff = HBAR_SI / (4.0 * np.pi * C_CM_S * freqs)
    S_free_rot = R_SI * (0.5 + np.log(np.sqrt(8.0 * np.pi**3 * I_eff * KB_SI * T_K / H_SI**2)))
    w = 1.0 / (1.0 + (nu0_cm1 / freqs)**4)
    S_vib = np.sum(w * S_vib_ho + (1.0 - w) * S_free_rot)

    H_kJ_mol = (E_scf_J_mol + H_trans + H_rot + H_vib) / 1000.0
    S_J_mol_K = S_trans + S_rot + S_vib
    G_kJ_mol = H_kJ_mol - T_K * S_J_mol_K / 1000.0
    return {
        'H_gas_kJ_mol': H_kJ_mol,
        'S_gas_J_mol_K': S_J_mol_K,
        'G_gas_kJ_mol': G_kJ_mol,
        'G_gas_eV': G_kJ_mol / EV_TO_KJ_MOL,
    }
