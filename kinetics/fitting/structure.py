"""Block 14C — Tests of the model structure: what a measurement implies whatever the barriers are.

A barrier sets how fast a run moves; it cannot move where a reaction stops (equilibrium) or which
compositions a run passes through (that is set by the rate law and the ratios of rate constants). These
functions turn measurements into statements of that kind, using stoichiometry and thermodynamics only.

    calculate_released_silyl_M     TMS groups released from TMSPA, from the measured ³¹P shares
    calculate_path_distance        how close a simulated composition path comes to a measured composition
    calculate_equilibrium_locus    ΔG_rxn of R2, R3 (and R4) implied if a measured composition is an equilibrium
    find_closest_on_locus          the point of the locus closest to the computed ΔG_rxn, in standard errors
    calculate_rt_equivalent_hours  room-temperature time that a heating history is worth, for a given barrier
    calculate_required_water_M     water needed to convert a measured amount of TMSPA

Classification: PUBLIC (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván
"""

import numpy as np
import pandas as pd

from kinetics.constants import KB_EV, ZERO_CELSIUS_K
from kinetics.data.snapshot import calculate_molar_mass

# TMS groups each phosphate has lost relative to TMSPA
SILYL_RELEASED = {'TMSPA': 0, 'BMSPA': 1, 'MMSPA': 2, 'H3PO4': 3}


def calculate_released_silyl_M(shares: dict, tmspa0_M: float) -> float:
    """Concentration [M] of TMS groups released from TMSPA, given the share of all P in each phosphate."""
    return tmspa0_M * sum(SILYL_RELEASED[sp] * shares.get(sp, 0.0) for sp in SILYL_RELEASED)


def calculate_path_distance(path: pd.DataFrame, shares: dict, sigmas: dict) -> dict:
    """Closest approach of a composition path (rows of shares) to a measured composition.

    The distance of a path point is the largest |share − measured| / σ over the measured species, so a
    distance ≤ 2 means every species is within 2 standard errors at once. Returns the smallest distance and
    the path row where it occurs.
    """
    species = list(shares)
    z = np.max([np.abs(path[sp].to_numpy() - shares[sp]) / sigmas[sp] for sp in species], axis=0)
    best = int(np.argmin(z))
    return {'distance': float(z[best]), 'row': path.iloc[best]}


def calculate_equilibrium_locus(shares: dict, *, tmspa0_M: float, h2o0_M: float, T_K: float,
                                n_points: int = 60) -> pd.DataFrame:
    """ΔG_rxn of R2, R3 and R4 [eV] for which a measured phosphate composition would be an equilibrium.

    The released TMS groups end either as TMSOH or, two at a time, as HMDSO; each HMDSO gives one water back
    (R4). With h the HMDSO formed [M]:
        [TMSOH] = released − 2h,   [H2O] = H2O₀ − released + h
        K2 = ([MMSPA]/[BMSPA])·[TMSOH]/[H2O],   K3 = ([H3PO4]/[MMSPA])·[TMSOH]/[H2O],
        K4 = [HMDSO][H2O]/[TMSOH]²,   ΔG = −k_B·T·ln K.
    h runs from 0 (no condensation: the most TMSOH, so the most negative ΔG2 and ΔG3) to just below
    released/2. Each row is one possible split; measuring TMSOH or water would pick one row.
    """
    released = calculate_released_silyl_M(shares, tmspa0_M)
    if released >= h2o0_M:
        raise ValueError('more TMS released than water added: the composition cannot come from hydrolysis alone')
    h = np.linspace(0.0, 0.5 * released, n_points + 1)[:-1]
    tmsoh = released - 2.0 * h
    h2o = h2o0_M - released + h
    ratio = tmsoh / h2o
    kT = KB_EV * T_K
    with np.errstate(divide='ignore'):
        dG4 = -kT * np.log(h * h2o / tmsoh ** 2)
    return pd.DataFrame({
        'HMDSO_M': h, 'TMSOH_M': tmsoh, 'H2O_M': h2o, 'TMSOH_to_H2O': ratio,
        'dG_R2_eV': -kT * np.log(shares['MMSPA'] / shares['BMSPA'] * ratio),
        'dG_R3_eV': -kT * np.log(shares['H3PO4'] / shares['MMSPA'] * ratio),
        'dG_R4_eV': dG4,
    })


def find_closest_on_locus(locus: pd.DataFrame, center_eV: dict, cov_eV2: pd.DataFrame) -> dict:
    """Point of the locus closest to the computed ΔG_rxn in Mahalanobis distance (standard errors).

    center_eV: computed ΔG_rxn per reaction id; cov_eV2: their covariance (kinetics.thermo
    calculate_solvation_covariance). Uses the reactions present in both (e.g. R2 and R3).
    """
    rxns = [r for r in center_eV if f'dG_{r}_eV' in locus.columns and r in cov_eV2.index]
    inv = np.linalg.inv(cov_eV2.loc[rxns, rxns].to_numpy())
    d = locus[[f'dG_{r}_eV' for r in rxns]].to_numpy() - np.array([center_eV[r] for r in rxns])
    dist = np.sqrt(np.einsum('ij,jk,ik->i', d, inv, d))
    best = int(np.argmin(dist))
    return {'distance': float(dist[best]), 'row': locus.iloc[best]}


def calculate_rt_equivalent_hours(heated: pd.DataFrame, dG_barrier_eV, *, T_ref_K: float) -> np.ndarray:
    """Hours at T_ref_K that give the same progress as the heating history, for each barrier ΔG‡ [eV].

    heated: rows with 'T_C' and 'minutes'. Each interval counts as minutes × k(T)/k(T_ref), with the Eyring
    ratio (T/T_ref)·exp[ΔG‡/k_B·(1/T_ref − 1/T)] at ΔS‡ = 0, as in the model.
    """
    dG = np.atleast_1d(np.asarray(dG_barrier_eV, dtype=float))
    total = np.zeros_like(dG)
    for _, seg in heated.iterrows():
        T = seg['T_C'] + ZERO_CELSIUS_K
        total += seg['minutes'] / 60.0 * (T / T_ref_K) * np.exp(dG / KB_EV * (1.0 / T_ref_K - 1.0 / T))
    return total


def calculate_required_water_M(tmspa_converted_M: float, *, density_g_mL: float) -> dict:
    """Water needed to convert tmspa_converted_M of TMSPA into BMSPA, as M and as ppm by mass of the liquid.

    One water hydrolyses one TMSPA (R1) and releases one TMSOH, which can convert a second TMSPA by transfer
    (R5). So the water needed lies between half the converted TMSPA (every TMSOH transfers) and all of it
    (none does).
    """
    to_ppm = calculate_molar_mass('H2O') / (density_g_mL * 1000.0) * 1e6
    low, high = 0.5 * tmspa_converted_M, tmspa_converted_M
    return {'low_M': low, 'high_M': high, 'low_ppm': low * to_ppm, 'high_ppm': high * to_ppm}
