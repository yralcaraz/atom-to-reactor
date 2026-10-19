"""Tests for the solvation medium (kinetics/data/medium.py) and the setup of a fit run (kinetics/fitting/config.py).

Source: Y. Alcaraz Galván; expected mixture energies as in P. Broqvist's tmspa_reaction_network_model notebook
(EC:DEC linear, 2026-10-08 snapshot: R1 −5.8, R2 −4.2, R3 −1.4, R4 +2.3 kcal/mol)
"""

import json
import math
import os
import sys
import tempfile

repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

import numpy as np

from kinetics.data.medium import DEFAULT_MEDIUM, SolvationMedium, get_active_medium, set_active_medium
from kinetics.data.snapshot import SNAPSHOT_2026_10_08_PATH, get_solvation_records
from kinetics.data.species import load_species_database, resolve_species_database
from kinetics.fitting.config import apply_fit_config, check_fit_medium, load_fit_config
from kinetics.thermo.reaction import calculate_network_thermo
from kinetics.thermo.uncertainty import calculate_solvation_covariance, calculate_solvation_sigma

EC_DEC = SolvationMedium('tank_api_snapshot_2026-10-08.json', (('EC', 0.5), ('DEC', 0.5)))
KCAL_TO_EV = 1.0 / 23.0605
HAVE_NEW = os.path.exists(SNAPSHOT_2026_10_08_PATH)


def test_default_medium_is_the_old_database():
    set_active_medium(None)
    assert get_active_medium() == DEFAULT_MEDIUM
    assert resolve_species_database(None) == load_species_database()


def test_linear_ec_dec_mixture():
    if not HAVE_NEW:
        print("skipped: no 2026-10-08 snapshot")
        return
    x = EC_DEC.mole_fractions()
    assert math.isclose(x['EC'], 0.645, abs_tol=0.002) and math.isclose(sum(x.values()), 1.0)
    ec, dec = get_solvation_records(EC_DEC.path, 'EC'), get_solvation_records(EC_DEC.path, 'DEC')
    try:
        set_active_medium(EC_DEC)
        db = resolve_species_database(None)
        expected = x['EC'] * ec['TMSPA']['raw_metadata']['dE_solv_eV'] + x['DEC'] * dec['TMSPA']['raw_metadata']['dE_solv_eV']
        assert math.isclose(db['TMSPA']['dE_solv_eV'], expected, rel_tol=1e-9)
        assert db['TMSPA']['solvent'] == 'EC:DEC'
        thermo = calculate_network_thermo(298.15)
        for rxn, kcal in {'R1': -5.8, 'R2': -4.2, 'R3': -1.4, 'R4': 2.3}.items():
            assert abs(thermo.loc[rxn, 'dG_rxn_eV'] - kcal * KCAL_TO_EV) < 0.01, (rxn, thermo.loc[rxn, 'dG_rxn_eV'])
        cov = calculate_solvation_covariance(['R1', 'R2', 'R3', 'R4'])
        for rxn in cov.index:
            assert math.isclose(math.sqrt(cov.loc[rxn, rxn]), thermo.loc[rxn, 'sigma_solv_eV'], rel_tol=1e-12)
        # weights x_s: the mixture's σ lies below that of either neat solvent
        sigma_neat = [calculate_solvation_sigma({'TMSPA': 1, 'H2O': 1}, {'BMSPA': 1, 'TMSOH': 1}, records)
                      for records in (ec, dec)]
        assert thermo.loc['R1', 'sigma_solv_eV'] < min(sigma_neat)
    finally:
        set_active_medium(None)
    assert np.isclose(calculate_network_thermo(298.15).loc['R1', 'dG_rxn_eV'],
                      calculate_network_thermo(298.15, species_db=load_species_database()).loc['R1', 'dG_rxn_eV'])


def test_fit_config_round_trip():
    with tempfile.TemporaryDirectory() as folder:
        default = load_fit_config(folder)
        assert default['medium'] == DEFAULT_MEDIUM and 'E1' in default['paper_observations']
        with open(os.path.join(folder, 'fit_config.json'), 'w', encoding='utf-8') as f:
            json.dump({'solvation': EC_DEC.to_dict(), 'scenarios': ['middle', 'free'], 'paper_observations': []}, f)
        try:
            config = apply_fit_config(folder)
            assert config['medium'] == EC_DEC == get_active_medium()
            assert config['scenarios'] == ['middle', 'free'] and config['paper_observations'] == []
            try:
                check_fit_medium({'theta': {}})            # a fit of the first three runs: default medium
            except ValueError:
                pass
            else:
                raise AssertionError('a default-medium fit must not run in the EC:DEC medium')
            check_fit_medium({'medium': EC_DEC.to_dict()})
        finally:
            set_active_medium(None)
        with open(os.path.join(folder, 'fit_config.json'), 'w', encoding='utf-8') as f:
            json.dump({'scenarios': ['soon']}, f)
        try:
            load_fit_config(folder)
        except ValueError:
            pass
        else:
            raise AssertionError('an unknown reading must be refused')


def test_pre_history_override_fixes_the_age():
    from kinetics.data.observables import DEFAULT_OBSERVABLES_PATH
    if not os.path.exists(DEFAULT_OBSERVABLES_PATH):
        print("skipped: no lab observables")
        return
    from kinetics.data import load_lab_observables
    from kinetics.fitting import build_sample_set
    observables = load_lab_observables()
    changed = {'2 % H2O': {'pre_history': [(None, 8.0), (30.0, 8.0)], 'heated_duration_h': 8.0}}
    for scenario in ('middle', 'free'):
        sample = next(s for s in build_sample_set(observables, scenario, overrides=changed) if s['name'] == '2 % H2O')
        assert sample['age_mode'] == 'none'
        assert [round(d / 3600.0, 6) for _, d in sample['pre_segments']] == [8.0, 8.0]
        assert math.isclose(sample['pre_segments'][0][0], sample['T_rt_K']) and math.isclose(sample['pre_segments'][1][0], 303.15)
        assert all(d <= 8.0 * 3600.0 + 1e-6 for _, _, d in sample['heated'])


if __name__ == "__main__":
    test_default_medium_is_the_old_database()
    test_linear_ec_dec_mixture()
    test_fit_config_round_trip()
    test_pre_history_override_fixes_the_age()
    print("All medium and fit-setup tests passed.")
