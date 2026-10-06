"""Tests for Block 14: feasibility bounds from the control experiments, the water-series checks and experiment design.

Source: Y. Alcaraz Galván; expected bounds computed from the Tank dataset snapshot (P. Broqvist) and Gogoi et al., J. Phys. Chem. C 2024, 128, 1654
"""

import os
import sys

import numpy as np
import pandas as pd

repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

from kinetics import calculate_recipe_molarities, get_model
from kinetics.constants import EV_TO_KJ_MOL, KB_EV
from kinetics.fitting import (
    build_composition, build_in_situ_design, calculate_barrier_error_budget, calculate_design_information,
    calculate_eyring_precision, calculate_parameter_precision, find_barrier_bounds, get_family_parameter,
    project_bounds_to_entropy, simulate_design, split_family_by_reaction,
)
from kinetics.fitting.readouts import (
    build_nmr_readouts, find_reachable_species, find_site_atoms, find_unread_species,
)
from kinetics.reactor import calculate_water_series_c0, evaluate_water_series, simulate_phosphate_path
from kinetics.fitting import (
    SHARED_BARRIER, build_observation, build_share_window, calculate_equilibrium_locus, calculate_path_distance,
    calculate_required_water_M, calculate_rt_equivalent_hours, evaluate_observation, find_allowed_intervals,
    find_closest_on_locus, intersect_intervals, observations_from_controls, replace_observation, summarize_region,
)
from kinetics.reactor import simulate_control_experiment
from kinetics.data import load_experimental_data

COARSE_GRID_EV = np.round(np.arange(0.60, 1.71, 0.10), 2)


def test_family_parameter_updates():
    model = get_model('level1')
    variant = model.with_family_params('hydrolysis', g_eV=1.2, dS_act_J_mol_K=-50.0)
    assert get_family_parameter(variant, 'hydrolysis') == 1.2
    assert get_family_parameter(variant, 'hydrolysis', 'dS_act_J_mol_K') == -50.0
    assert get_family_parameter(variant, 'transfer') == get_family_parameter(model, 'transfer')
    assert get_family_parameter(model, 'hydrolysis') == 0.80, 'the registered model must not change'
    # A family without its own entry starts from 'default'
    reference = get_model('peter_reference').with_family_params('condensation', g_eV=1.4)
    assert get_family_parameter(reference, 'condensation') == 1.4
    assert get_family_parameter(reference, 'hydrolysis') == 1.15


def test_feasibility_bounds_level1():
    model = get_model('level1')
    bounds = find_barrier_bounds(model, ['solvent_attack', 'transfer'], g_grid_eV=COARSE_GRID_EV)
    direct = bounds[bounds['acts'] == 'directly']
    sa = direct[direct['family'] == 'solvent_attack']
    assert abs(sa[sa['bound'] == 'g ≥']['g_eV'].max() - 1.273) < 0.005, sa
    assert abs(sa[sa['bound'] == 'g ≤']['g_eV'].min() - 1.387) < 0.005, sa
    tr = direct[direct['family'] == 'transfer']
    assert len(tr) == 1 and tr.iloc[0]['bound'] == 'g ≤' and abs(tr.iloc[0]['g_eV'] - 1.129) < 0.005, tr
    assert find_barrier_bounds(model, ['hydrolysis'], g_grid_eV=(0.9, 1.3, 1.7)).empty, 'no control contains water'
    # A bound on g(T_exp) moves along ΔS‡ with slope (T_exp − T_ref)
    proj = project_bounds_to_entropy(direct, 'solvent_attack', [0.0, -100.0])
    e2 = proj[(proj['experiment'] == 'E2') & (proj['bound'] == 'g ≤')].set_index('dS_J_mol_K')['g_ref_eV']
    assert abs((e2[0.0] - e2[-100.0]) - 55.0 * 100.0 / 1000.0 / EV_TO_KJ_MOL) < 1e-9


def test_water_series_contradiction():
    water = evaluate_water_series(['level1'], method='BDF')
    assert not water['all ok'].any(), 'level1 should not reproduce the water series at any common time'
    assert not (water['W1 ok'] & water['W2 ok']).any()


def test_readouts_follow_measured_shifts():
    readouts = build_nmr_readouts()
    # ³¹P: four separate signals; ²⁹Si: the three silyl phosphates are one reported signal (xMSPA)
    assert readouts['P'] == {sp: {sp: 1} for sp in ('TMSPA', 'BMSPA', 'MMSPA', 'H3PO4')}
    assert readouts['Si']['xMSPA'] == {'TMSPA': 3, 'BMSPA': 2, 'MMSPA': 1}
    assert readouts['Si']['TMSOH'] == {'TMSOH': 1} and readouts['Si']['TMSOEG'] == {'TMSOEG': 1}
    # Methyl sites: 3 C and 9 H per TMS group
    assert readouts['C']['HMDSO'] == {'HMDSO': 6} and readouts['H']['xMSPA']['TMSPA'] == 27
    assert len(find_site_atoms('EC', 'C', 'Si-CH3')) == 0
    # Species with no measured shift are not read
    assert find_unread_species('P') == [] and find_unread_species('Si') == ['TMSOdiEG']
    # TMSOH alone in EC can form the ring-opened products but no phosphate
    reachable = find_reachable_species({'TMSOH': 0.45, 'EC': 15.0})
    assert 'TMSOdiEG' in reachable and 'BMSPA' not in reachable


def test_design_readouts_and_information():
    wet = calculate_recipe_molarities(0.02, 0.05)['after']
    c0 = build_composition(wet['EC'], TMSPA=wet['TMSPA'], H2O=wet['H2O'])
    design = build_in_situ_design('t', 'test', c0, segments=((40, 2), (60, 2)), nuclei=('P', 'Si'), dead_time_h=0.25)
    model = get_model('level1').with_family_params('hydrolysis', g_eV=1.25)
    run = simulate_design(design, model)
    total_P = sum(run['readouts']['P'].values())
    total_Si = sum(run['readouts']['Si'].values())
    assert np.allclose(total_P, wet['TMSPA'], rtol=1e-6), 'P readouts must close the phosphorus balance'
    t_Si = np.asarray(run['sampling_h']['Si'])
    unread_Si = sum(np.interp(t_Si, run['t_h'], run['C_M'][run['idx'][sp]]) for sp in find_unread_species('Si'))
    assert np.allclose(total_Si + unread_Si, 3.0 * wet['TMSPA'], rtol=1e-6), \
        'Si readouts plus the unread Si species must close the silicon balance'
    # Splitting a family per reaction leaves the kinetics unchanged
    split, network, steps = split_family_by_reaction(model, 'hydrolysis')
    assert steps == ['hydrolysis_R1', 'hydrolysis_R2', 'hydrolysis_R3']
    assert np.allclose(split.calculate_rates(313.15, network=network)['k_f'].values,
                       model.calculate_rates(313.15)['k_f'].values)
    # At 40 °C for 2 h EC ring-opening does not happen, so the data cannot inform solvent attack
    design_40C = build_in_situ_design('t40', 'test', c0, segments=((40, 2),), nuclei=('P',), dead_time_h=0.25)
    info = calculate_design_information(design_40C, model, ['hydrolysis', 'solvent_attack'])
    sigma, corr = calculate_parameter_precision(info, ['hydrolysis', 'solvent_attack'])
    assert sigma['hydrolysis'] < 0.01, sigma
    assert abs(sigma['solvent_attack'] - 0.30) < 0.01, 'a parameter the data do not see keeps the prior σ'
    assert abs(corr.loc['hydrolysis', 'hydrolysis'] - 1.0) < 1e-12


def test_error_budget_and_eyring():
    budget = calculate_barrier_error_budget(1.0, 40.0, dT_K=0.5, rel_conc_error=0.05, sigma_dG_rxn_eV={'x': 0.2})
    kT = KB_EV * 313.15
    assert abs(budget.iloc[1]['on measured ΔG‡ (meV)'] - 1000.0 * kT * np.log(1.05)) < 1e-6
    assert budget.iloc[2]['on measured ΔG‡ (meV)'] == 0.0 and abs(budget.iloc[2]['on family g (meV)'] - 100.0) < 1e-9
    eyring = calculate_eyring_precision({'two': (40, 60)}, sigma_g_eV=0.003)
    expected = 0.003 * EV_TO_KJ_MOL * 1000.0 / np.sqrt(2 * 10.0 ** 2)
    assert abs(eyring.loc['two', 'σ ΔS‡ (J/mol/K)'] - expected) < 1e-9


def test_observations_reproduce_control_checks():
    controls = load_experimental_data()['control_experiments']
    observations = observations_from_controls()
    model = get_model('level1')
    for exp, obs in zip(controls, observations):
        value = simulate_control_experiment(exp, model, method='BDF')
        check = evaluate_observation(obs, model)
        assert abs(check['values'][exp['observable']] - value) < 1e-9
    # Same bounds as feasibility.find_barrier_bounds for solvent attack (E1, E2)
    sa = [o for o in observations if o['id'] in ('E1', 'E2')]
    region = intersect_intervals(find_allowed_intervals(model, 'solvent_attack', sa, g_grid_eV=COARSE_GRID_EV)).iloc[0]
    assert abs(region['g_low_eV'] - 1.273) < 0.005 and abs(region['g_high_eV'] - 1.387) < 0.005, region
    assert summarize_region(find_allowed_intervals(model, 'solvent_attack', sa, g_grid_eV=COARSE_GRID_EV)).startswith('1.27')


def test_time_interval_and_windows():
    assert build_share_window(0.93, 0.02, 'at least') == (0.89, None)
    assert build_share_window(0.0, 0.03, 'at most', k_sigma=2.0) == (None, 0.06)
    assert build_share_window(0.5, 0.4, 'both') == (0.0, 1.0), 'windows are clipped to [0, 1]'
    c0 = calculate_water_series_c0(2.0)
    model = get_model('level1').with_family_params('hydrolysis', g_eV=1.25)
    gone = build_observation('x', 'TMSPA gone', c0_M=c0, T_K=298.15, t_min_s=60.0, t_max_s=3600.0,
                             windows={'tmspa_share': (None, 0.05)}, status='measured', source='test',
                             time_basis='test', probes='hydrolysis')
    assert not evaluate_observation(gone, model)['consistent'], 'TMSPA cannot be gone within 1 h at g = 1.25'
    # A wider time interval can only help: some time in it may satisfy the window
    longer = replace_observation(gone, t_max_s=3600.0 * 24 * 365)
    assert evaluate_observation(longer, model)['consistent']
    # One shared barrier (M0): every family moves together
    shared = find_allowed_intervals('level1', SHARED_BARRIER, [gone], g_grid_eV=COARSE_GRID_EV)
    assert shared.iloc[0]['verdict'] == 'upper bound only'


def test_structure_tests():
    shares = {'BMSPA': 0.04, 'MMSPA': 0.31, 'H3PO4': 0.65}
    locus = calculate_equilibrium_locus(shares, tmspa0_M=0.149, h2o0_M=1.05, T_K=298.15)
    row = locus.iloc[0]          # no condensation: TMSOH = released TMS
    released = 0.149 * (0.04 + 2 * 0.31 + 3 * 0.65)
    assert abs(row['TMSOH_M'] - released) < 1e-12 and abs(row['H2O_M'] - (1.05 - released)) < 1e-12
    K2 = 0.31 / 0.04 * row['TMSOH_M'] / row['H2O_M']
    assert abs(row['dG_R2_eV'] + KB_EV * 298.15 * np.log(K2)) < 1e-12
    assert (np.diff(locus['dG_R3_eV']) > 0).all(), 'more condensation means less TMSOH and a less negative ΔG3'
    cov = pd.DataFrame([[0.01, 0.0], [0.0, 0.04]], index=['R2', 'R3'], columns=['R2', 'R3'])
    closest = find_closest_on_locus(locus, {'R2': row['dG_R2_eV'] - 0.1, 'R3': row['dG_R3_eV']}, cov)
    assert abs(closest['distance'] - 1.0) < 0.05
    # Heating at the reference temperature counts one to one
    heated = pd.DataFrame({'T_C': [25.0, 80.0], 'minutes': [60.0, 60.0]})
    hours = calculate_rt_equivalent_hours(heated, [0.0, 1.0], T_ref_K=298.15)
    eyring_ratio = 353.15 / 298.15 * np.exp(1.0 / KB_EV * (1.0 / 298.15 - 1.0 / 353.15))
    assert abs(hours[0] - (1.0 + 353.15 / 298.15)) < 1e-9 and abs(hours[1] - (1.0 + eyring_ratio)) < 1e-6
    water = calculate_required_water_M(0.08, density_g_mL=1.0)
    assert water['low_M'] == 0.04 and abs(water['high_ppm'] - 0.08 * 18.015 * 1000.0) < 1e-6
    # A barrier moves a run along its path but does not change it
    c0 = calculate_water_series_c0(2.0)
    paths = [simulate_phosphate_path(c0, get_model('level1').with_family_params('hydrolysis', g_eV=g), t_end_s=3.2e10)
             for g in (1.1, 1.3)]
    targets = {'BMSPA': 0.04, 'MMSPA': 0.31, 'H3PO4': 0.65}
    sigmas = dict.fromkeys(targets, 0.02)
    d = [calculate_path_distance(path, targets, sigmas)['distance'] for path in paths]
    assert min(d) > 3.0 and abs(d[0] - d[1]) < 1.0, d


if __name__ == "__main__":
    test_family_parameter_updates()
    test_feasibility_bounds_level1()
    test_water_series_contradiction()
    test_readouts_follow_measured_shifts()
    test_design_readouts_and_information()
    test_error_budget_and_eyring()
    test_observations_reproduce_control_checks()
    test_time_interval_and_windows()
    test_structure_tests()
    print("All fitting and design tests passed.")
