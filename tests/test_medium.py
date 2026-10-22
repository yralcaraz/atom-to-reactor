"""Tests for the solvation medium (kinetics/data/medium.py), the network variants (kinetics/data/network.py) and the
setup of a fit run (kinetics/fitting/config.py, with the run-wide changes to the samples of kinetics/fitting/residuals.py).

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
from dataclasses import replace

from kinetics.data.snapshot import SNAPSHOT_2026_10_08_PATH, get_cosmors_records, get_solvation_records
from kinetics.data.species import load_species_database, resolve_species_database
from kinetics.fitting.config import apply_fit_config, check_fit_medium, load_fit_config
from kinetics.thermo.reaction import calculate_network_thermo
from kinetics.fitting.candidates import _energy_tools
from kinetics.thermo.uncertainty import calculate_method_gap, calculate_solvation_covariance, calculate_solvation_sigma

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


def test_cosmors_medium():
    """Fit 6.2: COSMO-RS solvation mixed like the MD; MD standard errors kept; the gap to the MD enters the prior."""
    if not HAVE_NEW:
        print("skipped: no 2026-10-08 snapshot")
        return
    cosmo = replace(EC_DEC, method='cosmors')
    assert SolvationMedium.from_dict(cosmo.to_dict()) == cosmo and 'method' not in EC_DEC.to_dict()
    x = cosmo.mole_fractions()
    expected = x['EC'] * get_cosmors_records(cosmo.path, 'EC')['TMSPA'] + x['DEC'] * get_cosmors_records(cosmo.path, 'DEC')['TMSPA']
    freed = ('R1', 'R2', 'R3', 'R4')
    try:
        set_active_medium(EC_DEC)
        md_thermo = calculate_network_thermo(298.15)
        md_cov = calculate_solvation_covariance(list(freed)).to_numpy()
        assert np.allclose(calculate_method_gap(freed), 0.0)
        set_active_medium(cosmo)
        assert math.isclose(resolve_species_database(None)['TMSPA']['dE_solv_eV'], expected, rel_tol=1e-12)
        thermo = calculate_network_thermo(298.15)
        gap = calculate_method_gap(freed)
        for rxn in freed:          # computed energy moves by the gap; same MD standard error
            assert math.isclose(thermo.loc[rxn, 'dG_rxn_eV'] - md_thermo.loc[rxn, 'dG_rxn_eV'], gap[rxn], abs_tol=1e-9)
        assert abs(gap['R1'] - 0.387) < 0.002 and abs(gap['R4'] + 0.233) < 0.002
        assert np.allclose(calculate_solvation_covariance(list(freed)).to_numpy(), md_cov)
        sigma = _energy_tools(freed)[3]
        assert np.allclose(sigma, np.sqrt(np.diag(md_cov) + gap.to_numpy() ** 2))
    finally:
        set_active_medium(None)


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


def test_network_variants():
    from kinetics.data.network import (EXTRA_REACTIONS, build_network, get_active_network_variant,
                                       set_active_network_variant)
    from kinetics.fitting.candidates import FIXED_BARRIERS_EV, get_structure
    assert get_active_network_variant() == 'R1-R9' and set(build_network()) == {f'R{i}' for i in range(1, 10)}
    thermo = calculate_network_thermo(298.15, network=build_network('R1-R9+N10-N15'))['dG_rxn_eV']
    # each added step is a combination of R1-R9 (Peter's notebook), so its ΔG is that combination exactly
    for rxn, combination in {'N10': {'R2': 1, 'R1': -1}, 'N11': {'R9': 1, 'R8': -1}, 'N12': {'R2': 1, 'R7': 1},
                             'N13': {'R2': 1, 'R3': -1}, 'N14': {'R3': 1, 'R1': -1}, 'N15': {'R1': 1, 'R6': 1}}.items():
        assert abs(thermo[rxn] - sum(c * thermo[r] for r, c in combination.items())) < 1e-9, rxn
    try:
        set_active_network_variant('R1-R9+N10-N15')
        split, plain = get_structure('M3-split'), get_structure('M3')
        assert split.network()['N15']['class'] == 'hydrolysis_R1' and split.network()['N12']['class'] == 'hydrolysis_R23'
        assert plain.network()['N15']['class'] == plain.network()['N12']['class'] == 'hydrolysis'
        names = [p.name for p in split.parameters()]
        assert names == ['g hydrolysis_R1', 'g hydrolysis_R23', 'g transfer', 'g condensation', 'g solvent_attack',
                         'dG R1', 'dG R2', 'dG R3', 'dG R4'], names          # no new free parameter
        theta = {'g hydrolysis_R1': 1.4, 'g hydrolysis_R23': 0.95, 'g transfer': 0.7, 'g condensation': 1.3,
                 'g solvent_attack': 1.23, 'dG R1': 0.1, 'dG R2': 0.1, 'dG R3': 0.1, 'dG R4': -0.1}
        model = split.build(theta)['model']
        assert model.family_params['exchange']['g_eV'] == FIXED_BARRIERS_EV['exchange']
        rates = model.calculate_rates(295.65)
        assert set(EXTRA_REACTIONS) <= set(rates.index)
        # with the shifts of the freed energies, the combinations still hold
        assert abs(rates.loc['N15', 'dG_rxn_eV'] - rates.loc['R1', 'dG_rxn_eV'] - rates.loc['R6', 'dG_rxn_eV']) < 1e-9
        shared = get_structure('M0-Marcus').build({'g all reactions': 1.0})['model']
        assert shared.family_params['exchange']['g_eV'] == FIXED_BARRIERS_EV['exchange']
        set_active_network_variant('R1-R9+N10-N15, Peter barriers')
        assert all(get_structure('M3-split').network()[r]['class'] == 'exchange' for r in EXTRA_REACTIONS)
    finally:
        set_active_network_variant(None)
    assert get_structure('M3-split').network() == {**build_network(), 'R1': {**build_network()['R1'], 'class': 'hydrolysis_R1'},
                                                   'R2': {**build_network()['R2'], 'class': 'hydrolysis_R23'},
                                                   'R3': {**build_network()['R3'], 'class': 'hydrolysis_R23'}}


def test_heating_before_the_first_spectrum():
    from kinetics.fitting.residuals import build_sample_history, merge_sample_overrides
    hour = 3600.0
    sample = {'kind': 'lab', 'name': 'tube', 'pre_segments': [], 'age_h': 5.0, 'T_rt_K': 296.0,
              'heated': [(-9 * hour, 303.15, 8 * hour), (hour, 313.15, hour)], 'blocks': [{'t_spec_s': np.array([0.0, 3 * hour])}]}
    segments, times = build_sample_history(sample)
    # the 30 °C step is cut at mixing (age 5 h): 4 h at 30 °C, 1 h at room temperature, then the step after the spectrum
    assert [(round(T, 2), round(d / hour, 6)) for T, d in segments] == [(303.15, 4.0), (296.0, 1.0), (296.0, 1.0),
                                                                        (313.15, 1.0), (296.0, 1.0)]
    assert np.allclose(times[0] / hour, [5.0, 8.0])
    segments, _ = build_sample_history(sample, age_h=24.0)
    assert [(round(T, 2), round(d / hour, 6)) for T, d in segments][:3] == [(296.0, 15.0), (303.15, 8.0), (296.0, 1.0)]
    base = {'2 % H2O': {'heated_episodes': [[1.0, 30.0, 8.0]]}}
    merged = merge_sample_overrides(base, {'2 % H2O': {'added_M': {'TMSOH': 1e-6}}, 'other': {'hold_h': 2.0}})
    assert merged == {'2 % H2O': {'heated_episodes': [[1.0, 30.0, 8.0]], 'added_M': {'TMSOH': 1e-6}}, 'other': {'hold_h': 2.0}}
    assert merge_sample_overrides(base, {'2 % H2O': {'heated_episodes': None}})['2 % H2O']['heated_episodes'] is None
    assert base == {'2 % H2O': {'heated_episodes': [[1.0, 30.0, 8.0]]}}                 # the run's own are not changed
    cooled = {'kind': 'lab', 'name': 'tube', 'blocks': [{'t_spec_s': np.array([0.0, 20 * hour, 40 * hour])}],
              'heated': [(5 * hour, 323.15, 8 * hour), (30 * hour, 343.15, 8 * hour)]}
    from kinetics.fitting.residuals import _cool_before_spectra
    # 50 °C ends at 13 h, 7 h before the spectrum at 20 h: kept; 70 °C ends at 38 h, 2 h before 40 h: moved 1 h earlier
    assert [(round(s / hour, 6), round(d / hour, 6)) for s, _, d in _cool_before_spectra(cooled, 3.0)] == [(5.0, 8.0), (29.0, 8.0)]


def test_fit_config_network_and_samples():
    from kinetics.data.network import get_active_network_variant
    from kinetics.fitting.config import check_fit_setup
    from kinetics.fitting.residuals import get_sample_overrides
    steps = {'2 % H2O': {'heated_episodes': [[13.3, 30.0, 8.0]]}}
    with tempfile.TemporaryDirectory() as folder:
        assert load_fit_config(folder)['network'] == 'R1-R9' and load_fit_config(folder)['sensitivities'] is None
        with open(os.path.join(folder, 'fit_config.json'), 'w', encoding='utf-8') as f:
            json.dump({'network': 'R1-R9+N10-N15', 'sample_overrides': steps}, f)
        try:
            apply_fit_config(folder)
            assert get_active_network_variant() == 'R1-R9+N10-N15' and get_sample_overrides() == steps
            check_fit_setup({'network': 'R1-R9+N10-N15', 'sample_overrides': steps})
            for stored in ({}, {'network': 'R1-R9+N10-N15', 'sample_overrides': {}}):
                try:
                    check_fit_setup(stored)          # a fit of fits 1-4, or one made without the steps
                except ValueError:
                    pass
                else:
                    raise AssertionError(f'{stored} must not run in this setup')
        finally:
            apply_fit_config(os.path.join(folder, 'nowhere'))    # back to the default setup
        assert get_active_network_variant() == 'R1-R9' and get_sample_overrides() == {}
        with open(os.path.join(folder, 'fit_config.json'), 'w', encoding='utf-8') as f:
            json.dump({'network': 'R1-R99'}, f)
        try:
            load_fit_config(folder)
        except ValueError:
            pass
        else:
            raise AssertionError('an unknown network must be refused')


if __name__ == "__main__":
    test_default_medium_is_the_old_database()
    test_linear_ec_dec_mixture()
    test_cosmors_medium()
    test_fit_config_round_trip()
    test_pre_history_override_fixes_the_age()
    test_network_variants()
    test_heating_before_the_first_spectrum()
    test_fit_config_network_and_samples()
    print("All medium and fit-setup tests passed.")
