"""Block 4 — Reaction thermodynamics: ΔG_rxn, K_eq and Wegscheider cycle closure.

Classification: PUBLIC (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván
"""

import numpy as np
import pandas as pd

from kinetics.constants import EV_TO_KJ_MOL, R_SI, T_STD_K
from kinetics.data.network import NETWORK, find_reaction_cycles, format_equation
from kinetics.thermo.gas import calculate_gas_thermo
from kinetics.thermo.solution import calculate_solution_gibbs
from kinetics.data.species import resolve_species_database
from kinetics.thermo.uncertainty import calculate_solvation_sigma


def calculate_reaction_thermo(rxn_id: str, 
                              T_K: float, 
                              *, 
                              network: dict = None, 
                              species_db: dict = None,
                              thermo_mode: str = 'wb97mv') -> dict:
    
    """
    [Checked - YA]
    ΔH, ΔS, ΔG and K_eq = exp(-ΔG/RT) of one reaction in EC solution at T_K.

    Computes reaction thermodynamics for the whole reaction network at T

    """
    
    rxn = (network or NETWORK)[rxn_id]
    dG_kJ_mol = dH_kJ_mol = dS_J_mol_K = 0.0
    for side, sign in (('products', 1.0), ('reactants', -1.0)):
        for sp, nu in rxn[side].items():
            th = calculate_solution_gibbs(sp, T_K, species_db=species_db, thermo_mode=thermo_mode)
            dG_kJ_mol += sign * nu * th['G_sol_kJ_mol']
            dH_kJ_mol += sign * nu * th['H_sol_kJ_mol']
            dS_J_mol_K += sign * nu * th['S_sol_J_mol_K']
    return {
        'rxn_id': rxn_id,
        'class': rxn['class'],
        'dH_rxn_kJ_mol': dH_kJ_mol,
        'dS_rxn_J_mol_K': dS_J_mol_K,
        'dG_rxn_kJ_mol': dG_kJ_mol,
        'dG_rxn_eV': dG_kJ_mol / EV_TO_KJ_MOL,
        'K_eq': np.exp(-dG_kJ_mol * 1000.0 / (R_SI * T_K)),
    }


def calculate_network_thermo(T_K: float = T_STD_K, *, network: dict = None, species_db: dict = None,
                             thermo_mode: str = 'wb97mv') -> pd.DataFrame:
    """One row per reaction: gas-phase ΔG, solvation contribution ΔΔE_solv ± σ, solution ΔG_rxn and K_eq.

    sigma_solv_eV is the ±1 standard error of ΔΔE_solv (NaN when the snapshot has no MD record).
    """
    net = network or NETWORK
    db = resolve_species_database(species_db)
    rows = {}
    for rxn_id, rxn in net.items():
        th = calculate_reaction_thermo(rxn_id, T_K, network=net, species_db=db, thermo_mode=thermo_mode)
        dG_gas_eV = ddE_solv_eV = 0.0
        for side, sign in (('products', 1.0), ('reactants', -1.0)):
            for sp, nu in rxn[side].items():
                dG_gas_eV += sign * nu * calculate_gas_thermo(sp, T_K, species_db=db, thermo_mode=thermo_mode)['G_gas_eV']
                ddE_solv_eV += sign * nu * db[sp]['dE_solv_eV']
        try:
            sigma_eV = calculate_solvation_sigma(rxn['reactants'], rxn['products'])
        except KeyError:
            sigma_eV = np.nan
        rows[rxn_id] = {
            'class': rxn['class'],
            'equation': format_equation(rxn),
            'dG_gas_eV': dG_gas_eV,
            'ddE_solv_eV': ddE_solv_eV,
            'sigma_solv_eV': sigma_eV,
            'dG_rxn_eV': th['dG_rxn_eV'],
            'dG_rxn_kJ_mol': th['dG_rxn_kJ_mol'],
            'K_eq': th['K_eq'],
        }
    return pd.DataFrame.from_dict(rows, orient='index')


def calculate_cycle_residuals(dG_rxn: dict, network: dict = None) -> list[tuple[str, float]]:
    """Wegscheider check: [(cycle label, Σ c_j ΔG_j)] for every independent cycle of the network.

    dG_rxn maps rxn_id → ΔG_rxn (any unit; the residual has the same unit).
    """
    net = network or NETWORK
    ids = list(net)
    out = []
    for cycle in find_reaction_cycles(net):
        terms = [(c, ids[j]) for j, c in enumerate(cycle) if c != 0]
        label = ''
        for k, (c, rxn_id) in enumerate(terms):
            coef = '' if abs(c) == 1 else f'{abs(c):g} '
            sign = ('' if c > 0 else '−') if k == 0 else (' + ' if c > 0 else ' − ')
            label += f'{sign}{coef}{rxn_id}'
        out.append((label, float(sum(c * dG_rxn[rxn_id] for c, rxn_id in terms))))
    return out
