"""Tests for the propagation of MD solvation uncertainties into reaction energies."""

import math
import os
import sys

repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

from kinetics.constants import EV_TO_KJ_MOL
from kinetics.data.snapshot import get_solvation_records
from kinetics.data.network import NETWORK
from kinetics.thermo import (
    calculate_solvation_covariance, calculate_solvation_sigma, calculate_solvation_sigma_naive,
)

SOLVATION = get_solvation_records()


def test_snapshot_uncertainty_is_quadrature_of_md_stds():
    for record in SOLVATION.values():
        raw = record['raw_metadata']
        total = math.sqrt(raw['E_gas_std_eV'] ** 2 + raw['E_solution_std_eV'] ** 2 + raw['E_solvent_std_eV'] ** 2)
        assert math.isclose(total * EV_TO_KJ_MOL, record['uncertainty_kjmol'], rel_tol=1e-6)


def test_shared_solvent_term_cancels_in_two_to_two_reactions():
    reactants, products = ['TMSPA', 'H2O'], ['BMSPA', 'TMSOH']
    expected = math.sqrt(sum(
        SOLVATION[s]['raw_metadata']['E_gas_std_eV'] ** 2 + SOLVATION[s]['raw_metadata']['E_solution_std_eV'] ** 2
        for s in reactants + products
    ))
    sigma = calculate_solvation_sigma(reactants, products)
    assert math.isclose(sigma, expected, rel_tol=1e-12)
    assert sigma < calculate_solvation_sigma_naive(reactants, products)


def test_list_and_dict_stoichiometry_agree():
    as_list = calculate_solvation_sigma(['TMSOH', 'TMSOH'], ['HMDSO', 'H2O'])
    as_dict = calculate_solvation_sigma({'TMSOH': 2}, {'HMDSO': 1, 'H2O': 1})
    assert math.isclose(as_list, as_dict, rel_tol=1e-12)


def test_covariance_matches_sigmas_and_shared_species():
    cov = calculate_solvation_covariance(['R2', 'R3', 'R4'])
    for r in cov.index:
        sigma = calculate_solvation_sigma(NETWORK[r]['reactants'], NETWORK[r]['products'])
        assert math.isclose(math.sqrt(cov.loc[r, r]), sigma, rel_tol=1e-12)
    assert math.isclose(cov.loc['R2', 'R3'], cov.loc['R3', 'R2'])
    # R2 and R3 share H2O and TMSOH on the same sides (+) and MMSPA on opposite sides (−)
    def solute_var(sp):
        raw = SOLVATION[sp]['raw_metadata']
        return raw['E_gas_std_eV'] ** 2 + raw['E_solution_std_eV'] ** 2
    expected = solute_var('H2O') + solute_var('TMSOH') - solute_var('MMSPA')
    assert math.isclose(cov.loc['R2', 'R3'], expected, rel_tol=1e-12)


if __name__ == "__main__":
    test_snapshot_uncertainty_is_quadrature_of_md_stds()
    test_shared_solvent_term_cancels_in_two_to_two_reactions()
    test_list_and_dict_stoichiometry_agree()
    test_covariance_matches_sigmas_and_shared_species()
    print("All solvation-uncertainty tests passed.")
