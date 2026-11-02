"""Checks of the thermodynamic fit: the equilibrium solver and a recovery on invented spectra.

    python tests/test_thermo_fit.py      (needs the data folder for the species energies; not the raw spectra)

1. Solver: for the recipes of the tubes, every reaction in the span of R1–R4 (N10, N12–N15 included) has Q = K at the
   solution, and phosphorus and silyl groups are conserved.
2. Recovery: shares computed at equilibrium with known solvation shifts, given the error model of a real block, are
   fitted back; the reaction energies the tube sees (R2, R3, R4 in the 2 % tube) come back within 0.02 eV.

Source: Y. Alcaraz Galván
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from kinetics.data.medium import SolvationMedium, set_active_medium
from kinetics.data.network import get_active_network, set_active_network_variant
from kinetics.fitting import fit_structure
from kinetics.fitting.residuals import _lab_block
from kinetics.fitting.thermo_fit import ThermoProblem, tabulate_reaction_energies
from kinetics.reactor.validation import calculate_water_series_c0
from kinetics.thermo.equilibrium import calculate_equilibrium, list_equilibrium_reactions
from kinetics.thermo.reaction import calculate_reaction_thermo

set_active_medium(SolvationMedium.from_dict({'snapshot': 'tank_api_snapshot_2026-10-08.json',
                                             'volume_fractions': {'EC': 0.5, 'DEC': 0.5}, 'rule': 'linear'}))
set_active_network_variant('R1-R9+N10-N15')
NET = get_active_network()
T_K = 295.65
RECIPES = {'TMSOH': {'TMSOH': 0.45, 'EC': 7.1},
           '2 % H2O': {sp: c for sp, c in calculate_water_series_c0(2.0).items() if c > 0},
           'TMSPa + TMSOH': {'TMSPA': 0.149, 'TMSOH': 0.18, 'EC': 7.1},
           'TMSPa + water': {'TMSPA': 0.149, 'H2O': 0.0656, 'EC': 7.1}}


def silyl(c):
    return 3 * c.get('TMSPA', 0) + 2 * c.get('BMSPA', 0) + c.get('MMSPA', 0) + c.get('TMSOH', 0) + 2 * c.get('HMDSO', 0)


def phosphorus(c):
    return sum(c.get(sp, 0) for sp in ('TMSPA', 'BMSPA', 'MMSPA', 'H3PO4'))


def test_solver():
    span = list_equilibrium_reactions(('R1', 'R2', 'R3', 'R4'), NET)
    assert span == ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'N10', 'N12', 'N13', 'N14', 'N15'], span
    for name, c0 in RECIPES.items():
        C = calculate_equilibrium(c0, T_K, network=NET)
        for rxn_id in span:
            rxn = NET[rxn_id]
            if any(C.get(sp, 0) <= 0 for side in ('reactants', 'products') for sp in rxn[side]):
                continue
            lnQ = sum(nu * np.log(C[sp]) for sp, nu in rxn['products'].items()) - sum(nu * np.log(C[sp]) for sp, nu in rxn['reactants'].items())
            lnK = np.log(calculate_reaction_thermo(rxn_id, T_K, network=NET)['K_eq'])
            assert abs(lnQ - lnK) < 1e-6, f'{name} {rxn_id}: ln Q {lnQ:.6f}, ln K {lnK:.6f}'
        assert abs(phosphorus(C) - phosphorus(c0)) < 1e-7 and abs(silyl(C) - silyl(c0)) < 1e-7, name   # the stop of the solver: 1e-8 of the largest conserved amount


def test_recovery():
    truth = {'TMSPA': 0.0, 'BMSPA': -0.06, 'MMSPA': 0.03, 'H3PO4': 0.05, 'H2O': -0.03, 'TMSOH': 0.05, 'HMDSO': -0.05}
    sample = {'name': '2 % H2O', 'c0_M': RECIPES['2 % H2O'], 'T_K': T_K, 'unknown_M': {}, 'blocks': []}
    probe = ThermoProblem([sample], free_species=truth, xi=3.0, fixed={f'dE_solv {sp}': v for sp, v in truth.items()})
    C = calculate_equilibrium(sample['c0_M'], T_K, network=NET, species_shifts_eV=truth)
    shares = {'31P': {w: C[w] for w in ('TMSPA', 'BMSPA', 'MMSPA', 'H3PO4')},
              '13C': {'TMSOH': C['TMSOH'], 'HMDSO': 2 * C['HMDSO'], 'TMSOEG': 0.0,
                      'P-silyl': 3 * C['TMSPA'] + 2 * C['BMSPA'] + C['MMSPA']}}
    floor = {'31P': {'absolute': 0.01, 'relative': 0.05}, '13C': {'absolute': 0.03, 'relative': 0.15}}
    for nucleus, amounts in shares.items():
        total = sum(amounts.values())
        spectra = [{'t_since_first_h': 24.0 * k, 'temperature_C': 22.5,
                    'shares': {w: {'share': a / total, 'sigma_noise': 0.005, 'sigma_baseline': 0.005} for w, a in amounts.items()}}
                   for k in range(2)]
        sample['blocks'].append(_lab_block('2 % H2O', nucleus, 'fit', spectra, {}, floor, 1.0))
    assert probe.chi2(np.zeros(0)) < 1e-12
    problem = ThermoProblem([sample], free_species=truth, xi=3.0)
    fit = fit_structure(problem, n_screen=256, n_starts=8, max_nfev=60)
    fitted = tabulate_reaction_energies(problem.species_shifts(fit['theta']), T_K)['ΔG (eV)']
    true = tabulate_reaction_energies(truth, T_K)['ΔG (eV)']
    for rxn_id in ('R2', 'R3', 'R4'):
        assert abs(fitted[rxn_id] - true[rxn_id]) < 0.02, f'{rxn_id}: fitted {fitted[rxn_id]:+.3f}, true {true[rxn_id]:+.3f} eV'
    assert fit['chi2'] < 1e-3, fit['chi2']


if __name__ == '__main__':
    test_solver()
    test_recovery()
    print('test_thermo_fit: solver and recovery pass')
