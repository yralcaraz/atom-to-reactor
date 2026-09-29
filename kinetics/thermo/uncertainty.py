"""Block 3C — Propagation of MD solvation uncertainties into reaction energies.

Each snapshot solvation record reports three standard deviations from the MACE MD runs
(`raw_metadata`): the gas-phase solute (E_gas_std_eV), the solute-in-EC box
(E_solution_std_eV) and the pure-EC reference box (E_solvent_std_eV). The snapshot's
`uncertainty_kjmol` is their quadrature sum.

The pure-EC box is the same run for every species, so its term is fully correlated:
it enters a reaction only through the net change in the number of solutes, and cancels
in the 2 → 2 steps of the network. Adding the per-species totals in quadrature would
count it several times.

Open question: whether these stds are per-frame deviations or standard errors of the
mean is not documented, so the σ values returned here are provisional.
"""

import math
from collections import Counter

from kinetics.constants import KJ_MOL_TO_EV
from kinetics.data.snapshot import get_solvation_records


def calculate_solvation_sigma_components(solvation_records: dict = None) -> dict:
    """{species: (sigma_solute_eV, sigma_solvent_eV)}.

    sigma_solute combines the gas-phase and solution-box terms (independent per species);
    sigma_solvent is the shared pure-EC reference term.
    """
    records = solvation_records if solvation_records is not None else get_solvation_records()
    components = {}
    for name, record in records.items():
        raw = record['raw_metadata']
        components[name] = (math.hypot(raw['E_gas_std_eV'], raw['E_solution_std_eV']), raw['E_solvent_std_eV'])
    return components


def calculate_solvation_sigma(reactants, products, solvation_records: dict = None) -> float:
    """σ of ΔΔE_solv [eV] for one reaction, treating the pure-EC reference term as shared.

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
