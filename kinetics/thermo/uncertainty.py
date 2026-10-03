"""Block 3C — Propagation of MD solvation uncertainties into reaction energies.

Each snapshot solvation record reports three standard errors of the mean MD energies from the
MACE runs (`raw_metadata`, stored under the `*_std_eV` keys): the gas-phase solute (E_gas_std_eV),
the solute-in-EC box (E_solution_std_eV) and the pure-EC reference box (E_solvent_std_eV). The
snapshot's `uncertainty_kjmol` is their quadrature sum.

The pure-EC box is the same run for every species, so its term is fully correlated:
it enters a reaction only through the net change in the number of solutes, and cancels
in the 2 → 2 steps of the network. Adding the per-species totals in quadrature would
count it several times.

The values are standard errors already, so every σ returned here is the ±1 standard error of
the solvation contribution and is used as stored (no further division by √n).

Classification: PUBLIC (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván
"""

import math
from collections import Counter

import pandas as pd

from kinetics.constants import KJ_MOL_TO_EV
from kinetics.data.network import NETWORK
from kinetics.data.snapshot import get_solvation_records


def calculate_solvation_sigma_components(solvation_records: dict = None) -> dict:
    """{species: (sigma_solute_eV, sigma_solvent_eV)}.

    sigma_solute combines the gas-phase and solution-box standard errors (independent per species);
    sigma_solvent is the standard error of the shared pure-EC reference term.
    """
    records = solvation_records if solvation_records is not None else get_solvation_records()
    components = {}
    for name, record in records.items():
        raw = record['raw_metadata']
        components[name] = (math.hypot(raw['E_gas_std_eV'], raw['E_solution_std_eV']), raw['E_solvent_std_eV'])
    return components


def calculate_solvation_sigma(reactants, products, solvation_records: dict = None) -> float:
    """Standard error σ of ΔΔE_solv [eV] for one reaction, treating the pure-EC reference term as shared.

    reactants, products: lists of species (repeated for stoichiometry, e.g. ['TMSOH', 'TMSOH'])
    or {species: coefficient} dicts.
    """
    components = calculate_solvation_sigma_components(solvation_records)
    var_solute = 0.0
    net_solutes = 0
    sigma_solvent = None
    for name, nu in _net_stoichiometry(reactants, products).items():
        if name not in components:
            raise KeyError(f"No solvation record for '{name}'")
        sigma_solute, sigma_solvent = components[name]
        var_solute += (nu * sigma_solute) ** 2
        net_solutes += nu
    # Each solute subtracts one pure-EC reference energy, so the shared term scales with
    # the net change in the number of solutes
    var_solvent = (net_solutes * sigma_solvent) ** 2 if sigma_solvent is not None else 0.0
    return math.sqrt(var_solute + var_solvent)


def calculate_solvation_covariance(rxn_ids, network: dict = None, solvation_records: dict = None) -> pd.DataFrame:
    """Covariance [eV²] of the solvation contributions ΔΔE_solv of several reactions (DataFrame).

    Two reactions are correlated through the species they share: Cov(i, j) = Σ_s ν_si ν_sj σ_solute,s²,
    plus the shared pure-EC term scaled by each reaction's net change in solutes (zero for 2 → 2 steps).
    The diagonal equals calculate_solvation_sigma².
    """
    net = network or NETWORK
    components = calculate_solvation_sigma_components(solvation_records)
    nets = {r: _net_stoichiometry(net[r]['reactants'], net[r]['products']) for r in rxn_ids}
    sigma_solvent = next(iter(components.values()))[1]
    cov = pd.DataFrame(0.0, index=list(rxn_ids), columns=list(rxn_ids))
    for i in rxn_ids:
        for j in rxn_ids:
            shared = set(nets[i]) & set(nets[j])
            value = sum(nets[i][sp] * nets[j][sp] * components[sp][0] ** 2 for sp in shared)
            value += sum(nets[i].values()) * sum(nets[j].values()) * sigma_solvent ** 2
            cov.loc[i, j] = value
    return cov


def calculate_solvation_sigma_naive(reactants, products, solvation_records: dict = None) -> float:
    """σ from adding each species' snapshot `uncertainty_kjmol` in quadrature (overestimate)."""
    records = solvation_records if solvation_records is not None else get_solvation_records()
    var = 0.0
    for name, nu in _net_stoichiometry(reactants, products).items():
        var += (nu * records[name]['uncertainty_kjmol'] * KJ_MOL_TO_EV) ** 2
    return math.sqrt(var)


def _net_stoichiometry(reactants, products) -> dict:
    """{species: net coefficient} (products − reactants), zero entries dropped."""
    net = Counter()
    for sp, nu in Counter(products).items():
        net[sp] += nu
    for sp, nu in Counter(reactants).items():
        net[sp] -= nu
    return {sp: nu for sp, nu in net.items() if nu != 0}
