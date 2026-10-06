"""Tests for the first fit: sample histories, species shifts, residuals, structures, optimiser and profiles.

Covers:
1. simulate_history against the batch and design reactors, the call budget, a model with its own network.
2. Species shifts: Wegscheider cycles stay closed, the requested reaction shifts come out, registered models unchanged.
3. Residuals: scenarios, the replicate block, the one-sided paper term, ages found inside an evaluation.
4. Solver tolerance: search and report settings agree far inside the measurement error.
5. Structures, the optimiser and profiles: known parameters are recovered from synthetic shares, a parameter the
   samples cannot see is reported as not determined, nesting never raises χ².
6. With the lab file present: the windows of notebook 03 give its barrier intervals through the new forward path.

Source: Y. Alcaraz Galván; expected intervals read from notebooks/results/03 (computed from the Tank
dataset snapshot and the experimental NMR spectra)
"""

import os
import sys

import numpy as np
import pandas as pd

repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for path in (repo_root, os.path.dirname(os.path.abspath(__file__))):
    if path not in sys.path:
        sys.path.insert(0, path)

from kinetics.constants import ZERO_CELSIUS_K
from kinetics.data import DEFAULT_OBSERVABLES_PATH, NETWORK, NETWORK_SPECIES, load_experimental_data, load_lab_observables
from kinetics.fitting import (
    ExperimentDesign, FitProblem, calculate_parameter_directions, calculate_prediction_band, calculate_profile,
    compare_structures, embed_parent_theta, extend_profile, find_confidence_interval, fit_leave_one_out, fit_structure,
    get_structure,
    simulate_design, simulate_synthetic_shares, split_family_by_reaction, tabulate_reaction_barriers,
)
from kinetics.fitting.residuals import (
    FREE_AGE_BOUNDS_H, PENALTY, REPORT_SOLVER, SEARCH_SOLVER, _shares_from_state, build_sample_history,
    build_sample_set, calculate_predicted_shares, calculate_residuals, summarize_residuals, tabulate_residuals,
)
from kinetics.microkinetics.models import MODELS, get_model
from kinetics.reactor import simulate_batch_reactor, simulate_history
from kinetics.thermo import (
    build_shifted_species_database, calculate_cycle_residuals, calculate_network_thermo, calculate_solvation_covariance,
    calculate_solvation_sigma_components, calculate_species_shifts,
)
from synthetic_observables import build_synthetic_observables

RESULTS_03 = os.path.join(repo_root, 'notebooks', 'results', '03')
MID = (get_model('level1').with_family_params('hydrolysis', g_eV=1.28).with_family_params('transfer', g_eV=1.10)
       .with_family_params('solvent_attack', g_eV=1.315))


def _full(c0):
    return {**dict.fromkeys(NETWORK_SPECIES, 0.0), **c0}


# ------------------------------------------------------------------------------
# 1. Sample histories
# ------------------------------------------------------------------------------
def test_history_matches_existing_reactors():
    c0 = {'TMSPA': 0.149, 'H2O': 1.05, 'EC': 7.1}
    # One segment: the batch reactor with the same solver takes the same steps
    batch = simulate_batch_reactor(_full(c0), T_K=295.65, t_end_s=86400.0, model=MID, n_points=20, t_start_s=60.0, method='BDF')
    history = simulate_history(c0, [(295.65, 86400.0)], batch['t_s'], model=MID)
    assert history['success'] and history['species'] == list(NETWORK_SPECIES)
    assert np.max(np.abs(history['C_M'] - batch['C_M'])) < 1e-10, np.max(np.abs(history['C_M'] - batch['C_M']))
    # Several segments: the design reactor (Radau) agrees within the solver tolerances
    segments_h = ((22.5, 24.0), (60.0, 0.5), (22.5, 12.0), (80.0, 0.3), (22.5, 6.0))
    times_h = (1.0, 24.0, 24.25, 30.0, 36.6, 42.8)
    design = ExperimentDesign('test', 'test', _full(c0), segments_h, {'P': times_h})
    reference = simulate_design(design, MID, points_per_segment=5)
    history = simulate_history(c0, [(T + ZERO_CELSIUS_K, 3600.0 * h) for T, h in segments_h],
                               3600.0 * np.array(times_h), model=MID)
    for species in ('TMSPA', 'BMSPA', 'MMSPA', 'H3PO4'):
        expected = reference['readouts']['P'][species]
        assert np.max(np.abs(history['C_M'][history['idx'][species]] - expected)) < 2e-7, species
    # Sample times in any order, t = 0 included; times outside the history are refused
    shuffled = simulate_history(c0, [(295.65, 3600.0), (313.15, 3600.0)], [7200.0, 0.0, 1800.0], model=MID)
    assert shuffled['C_M'][shuffled['idx']['TMSPA'], 1] == 0.149
    assert shuffled['C_M'][shuffled['idx']['TMSPA'], 0] < shuffled['C_M'][shuffled['idx']['TMSPA'], 2] < 0.149
    try:
        simulate_history(c0, [(295.65, 3600.0)], [7200.0], model=MID)
        raise AssertionError('a sample time beyond the history must be refused')
    except ValueError:
        pass
    # An exhausted call budget is reported, not raised
    starved = simulate_history(c0, [(295.65, 86400.0)], [86400.0], model=MID, max_rhs_calls=5)
    assert not starved['success'] and np.isnan(starved['C_M']).all()
    print('  ✓ simulate_history: equals the batch reactor (1e-10 M) and the design reactor (2e-7 M); budget reported')


def test_model_with_own_network():
    split, network, steps = split_family_by_reaction('level1', 'hydrolysis')
    assert steps == ['hydrolysis_R1', 'hydrolysis_R2', 'hydrolysis_R3']
    variant = split.with_family_params('hydrolysis_R1', g_eV=1.30).with_family_params('hydrolysis_R2', g_eV=1.05)
    carried = variant.__class__(**{**variant.__dict__, 'network': network})
    c0 = {'TMSPA': 0.149, 'H2O': 1.05, 'EC': 7.1}
    explicit = simulate_history(c0, [(295.65, 86400.0)], [86400.0], model=variant, network=network)
    own = simulate_history(c0, [(295.65, 86400.0)], [86400.0], model=carried)
    assert np.array_equal(explicit['C_M'], own['C_M']), 'a model that carries its network needs none passed'
    batch = simulate_batch_reactor(_full(c0), T_K=295.65, t_end_s=86400.0, model=carried, n_points=5, t_start_s=60.0, method='BDF')
    assert np.max(np.abs(batch['C_M'][:, -1] - own['C_M'][:, 0])) < 1e-10
    rates = carried.calculate_rates(295.65)
    assert rates.loc['R1', 'g_eV'] == 1.30 and rates.loc['R2', 'g_eV'] == 1.05 and rates.loc['R3', 'g_eV'] == 0.80
    print('  ✓ a model with its own network is used as such by the history, batch and rate functions')


# ------------------------------------------------------------------------------
# 2. Species shifts
# ------------------------------------------------------------------------------
def test_species_shifts():
    x_hat = calculate_network_thermo(298.15)['dG_rxn_eV']
    wanted = {'R1': 0.05, 'R2': 0.12, 'R3': 0.17, 'R4': -0.08}
    shifts = calculate_species_shifts(wanted)
    db = build_shifted_species_database(shifts)
    moved = calculate_network_thermo(298.15, species_db=db)['dG_rxn_eV'] - x_hat
    for rxn_id, d in wanted.items():
        assert abs(moved[rxn_id] - d) < 1e-9, (rxn_id, moved[rxn_id])   # species energies are ~5e4 eV
    for rxn_id, hydrolysis in (('R5', 'R1'), ('R6', 'R2'), ('R7', 'R3')):     # transfer = hydrolysis + condensation
        assert abs(moved[rxn_id] - (wanted[hydrolysis] + wanted["R4"])) < 1e-9
    assert abs(moved["R8"]) < 1e-9 and abs(moved["R9"]) < 1e-9, 'the other reactions keep their energy'
    thermo = calculate_network_thermo(298.15, species_db=db)
    assert all(abs(res) < 1e-9 for _, res in calculate_cycle_residuals(thermo['dG_rxn_kJ_mol'].to_dict()))
    try:
        calculate_species_shifts({'R1': 0.1, 'R4': 0.0, 'R5': 0.3})
        raise AssertionError('a shift that breaks a Wegscheider cycle must be refused')
    except ValueError:
        pass
    # The most probable shift behind "R2 and R3 free": 0.75 standard errors (notebook 03 §5.3)
    sigma = {sp: s for sp, (s, _) in calculate_solvation_sigma_components().items()}
    d = np.array([0.120, 0.170])
    probable = calculate_species_shifts({'R2': d[0], 'R3': d[1]}, species_sigma_eV=sigma, keep_others=False)
    distance = np.sqrt(sum((v / sigma[sp]) ** 2 for sp, v in probable.items()))
    cov = calculate_solvation_covariance(['R2', 'R3']).to_numpy()
    assert abs(distance - np.sqrt(d @ np.linalg.solve(cov, d))) < 1e-9 and abs(distance - 0.75) < 0.01, distance
    assert abs(probable['TMSOH'] - 0.091) < 0.002 and abs(probable['H2O'] + 0.036) < 0.002, probable
    # A model that carries the shifts gives the rates of the shifted database; registered models carry none
    level1 = get_model('level1')
    carried = level1.__class__(**{**level1.__dict__, 'species_shifts_eV': shifts})
    assert np.array_equal(carried.calculate_rates(300.0)['k_f'].values, level1.calculate_rates(300.0, species_db=db)['k_f'].values)
    assert all(m.network is None and m.species_shifts_eV == {} for m in MODELS.values())
    assert np.array_equal(level1.calculate_rates(300.0)['k_f'].values,
                          level1.calculate_rates(300.0, network=NETWORK, species_db=None)['k_f'].values)
    print(f'  ✓ species shifts: reactions moved as asked, cycles closed, others unchanged; most probable shift {distance:.2f} SE')


# ------------------------------------------------------------------------------
# 3. Residuals
# ------------------------------------------------------------------------------
def test_sample_set_and_histories():
    observables = build_synthetic_observables()
    middle = {s['name']: s for s in build_sample_set(observables, 'middle')}
    assert list(middle) == ['0.5 % H2O', '2 % H2O', 'TMSPa + TMSOH (A)', 'TMSOH, probe', 'TMSOH, glovebox', 'E1']
    assert all(s['age_h'] == 24.0 for s in middle.values() if s['age_mode'] == 'fixed')
    assert middle['TMSOH, glovebox']['age_mode'] == 'none' and middle['2 % H2O']['replicate']
    assert not middle['TMSOH, probe']['replicate'], 'a spectrum at 80 °C is not a repeat of the room-temperature ones'
    segments, times = build_sample_history(middle['2 % H2O'])
    assert len(segments) == 14 and abs(sum(d for _, d in segments) / 3600.0 - (24.0 + 145.0)) < 1e-9
    assert [round(T - ZERO_CELSIUS_K) for T, _ in segments[2::2]] == [30, 40, 50, 60, 70, 80]
    assert np.allclose(times / 3600.0, 24.0 + np.array([0.1, 33.0, 72.0, 96.0, 120.0, 145.0]))
    segments, times = build_sample_history(middle['TMSOH, glovebox'])
    assert segments[0] == (353.15, 8 * 3600.0) and abs(times[0] / 3600.0 - 8.1) < 1e-9
    free = {s['name']: s for s in build_sample_set(observables, 'free')}
    assert [free[n]['age_mode'] for n in ('0.5 % H2O', '2 % H2O', 'TMSPa + TMSOH (A)', 'TMSOH, probe')] == \
        ['profile', 'parameter', 'profile', 'fixed']
    held = build_sample_set(observables, 'short', roles=('hold_out',), paper=())
    assert [s['name'] for s in held] == ['TMSPa + TMSOH (B)']
    assert 'E1' not in [s['name'] for s in build_sample_set(observables, 'short', exclude=('E1', '2 % H2O'))]
    longer = build_sample_set(observables, 'middle', overrides={'TMSOH, glovebox': {'hold_h': 16.0},
                                                               '2 % H2O': {'heated_duration_h': 8.0}})
    by_name = {s['name']: s for s in longer}
    assert by_name['TMSOH, glovebox']['pre_segments'][0][1] == 16 * 3600.0
    assert all(d == 8 * 3600.0 for _, _, d in by_name['2 % H2O']['heated'])
    build_sample_history(by_name['2 % H2O'])                 # 8 h holds still fit between the spectra
    seeded = {s['name']: s for s in build_sample_set(observables, 'middle', overrides={'2 % H2O': {'added_M': {'TMSOH': 1e-6}}})}
    assert seeded['2 % H2O']['c0_M']['TMSOH'] == middle['2 % H2O']['c0_M'].get('TMSOH', 0.0) + 1e-6
    assert seeded['0.5 % H2O']['c0_M'] == middle['0.5 % H2O']['c0_M'], 'an addition applies to the named sample only'
    assert 'TMSOH' not in observables['samples']['2 % H2O']['c0_M'], 'the observables record must not be changed'
    print('  ✓ sample sets: scenarios, known pre-history, heated segments, hold-out, exclusions and overrides')


def test_replicate_block_and_limits():
    observables = build_synthetic_observables()
    samples = build_sample_set(observables, 'middle', exclude=('0.5 % H2O', 'TMSPa + TMSOH (A)', 'TMSOH, probe',
                                                               'TMSOH, glovebox', 'E1'))
    block = samples[0]
    assert block['replicate'] and block['birge'] == 1.0
    from kinetics.fitting.residuals import _rows_of_sample, _whiten
    # A model that predicts the same composition in every spectrum, off by a common amount in one window
    offset = np.zeros_like(block['measured'])
    offset[:, 2] = 0.02
    rows = pd.DataFrame(_rows_of_sample(block, block['measured'] + offset, 24.0))
    drift = rows[rows['kind'] == 'drift']
    assert np.max(np.abs(drift['z'])) < 1e-9, 'six identical spectra and a constant prediction: no drift residual'
    mean = rows[rows['kind'] == 'mean'].set_index('quantity')
    assert abs(mean.loc['MMSPA (mean of 6 spectra)', 'z'] - 0.02 / mean.loc['MMSPA (mean of 6 spectra)', 'sigma']) < 1e-9
    # χ² of a common offset δ in one window: δ² / (σ_c² + 1/Σ(1/σ_i²)), times (K − 1)/K
    sigma_mean = mean.loc['MMSPA (mean of 6 spectra)', 'sigma']
    chi2 = np.sum(_whiten(block, block['measured'] + offset) ** 2)
    assert abs(chi2 - 0.75 * (0.02 / sigma_mean) ** 2) < 1e-9, (chi2, 0.75 * (0.02 / sigma_mean) ** 2)
    # A drift in one spectrum is weighed with the independent error only
    bump = np.zeros_like(block['measured'])
    bump[3, 1] = 0.01
    rows = pd.DataFrame(_rows_of_sample(block, block['measured'] + bump, 24.0))
    z = rows[(rows['kind'] == 'drift') & (rows['quantity'] == 'BMSPA drift, spectrum 4')]['z'].iloc[0]
    assert 2.0 < z < 0.01 / block['sigma_ind'][3, 1] + 1e-9, z
    # One-sided paper term: zero inside the limit, linear outside; a failed simulation costs PENALTY per value
    e1 = build_sample_set(observables, 'middle', exclude=tuple(observables['samples']))
    assert [s['name'] for s in e1] == ['E1']
    slow = get_model('level1')
    fast = slow.with_family_params('solvent_attack', g_eV=1.15)
    assert calculate_residuals(slow, e1)[0] == 0.0
    r_fast = calculate_residuals(fast, e1)[0]
    predicted = calculate_predicted_shares(e1[0], fast)
    assert r_fast > 3.0 and abs(r_fast - (predicted - 0.02) / 0.01) < 1e-12
    failed = calculate_residuals(MID, samples, solver={**SEARCH_SOLVER, 'max_rhs_calls': 3})
    assert failed.size == 24 and np.all(failed == PENALTY)
    print('  ✓ residuals: no drift for identical spectra, common offset weighed once, one-sided limit, failure penalty')


def test_free_ages_and_tables():
    observables = build_synthetic_observables()
    free = build_sample_set(observables, 'free')
    residuals, details = calculate_residuals(MID, free, ages={'2 % H2O': 100.0}, return_details=True)
    assert residuals.size == 45 and all(d['ok'] for d in details.values())
    # The age found inside the evaluation is the best of a brute-force scan
    tube = next(s for s in free if s['name'] == 'TMSPa + TMSOH (A)')
    from kinetics.fitting.residuals import _whiten
    ages = np.geomspace(*FREE_AGE_BOUNDS_H, 200)
    scan = [np.sum(_whiten(tube, calculate_predicted_shares(tube, MID, age_h=a)) ** 2) for a in ages]
    found = details['TMSPa + TMSOH (A)']['age_h']
    chi2_found = np.sum(_whiten(tube, calculate_predicted_shares(tube, MID, age_h=found)) ** 2)
    assert FREE_AGE_BOUNDS_H[0] <= found <= FREE_AGE_BOUNDS_H[1]
    assert chi2_found <= min(scan) + 0.02, (chi2_found, min(scan), found)
    # A fixed scenario is the free one evaluated at those ages
    middle = build_sample_set(observables, 'middle')
    same = calculate_residuals(MID, free, ages={'2 % H2O': 24.0, '0.5 % H2O': 24.0, 'TMSPa + TMSOH (A)': 24.0})
    fixed = calculate_residuals(MID, middle)
    probe = slice(32, 41)                                    # the probe sample sits at the lower bound in 'free'
    assert np.allclose(np.delete(same, probe), np.delete(fixed, probe), atol=1e-9)
    table = tabulate_residuals(MID, middle)
    summary = summarize_residuals(calculate_residuals(MID, middle, solver=REPORT_SOLVER), table)
    assert set(table['kind']) == {'share', 'mean', 'drift', 'limit'} and len(table) == 49
    assert summary['n_residuals'] == 45 and summary['max_abs_z'] == table['z'].abs().max()
    # Available water scales the water of the samples mixed with water only
    dry = calculate_residuals(MID, middle, water_fraction=0.5)
    assert not np.allclose(dry[:28], fixed[:28]) and np.allclose(dry[32:], fixed[32:])
    print('  ✓ free ages: the age found in one simulation matches a 200-point scan; tables and water fraction')


# ------------------------------------------------------------------------------
# 4. Solver tolerance
# ------------------------------------------------------------------------------
def test_search_and_report_tolerances_agree():
    observables = build_synthetic_observables()
    samples = build_sample_set(observables, 'middle')
    worst = 0.0
    for g_H, g_T, g_SA in ((1.28, 1.10, 1.315), (1.10, 0.90, 1.30), (1.40, 1.30, 1.40)):
        model = (get_model('level1').with_family_params('hydrolysis', g_eV=g_H).with_family_params('transfer', g_eV=g_T)
                 .with_family_params('solvent_attack', g_eV=g_SA))
        for sample in samples:
            if sample['kind'] == 'paper':
                continue
            coarse = calculate_predicted_shares(sample, model, solver=SEARCH_SOLVER)
            fine = calculate_predicted_shares(sample, model, solver=REPORT_SOLVER)
            worst = max(worst, np.max(np.abs(coarse - fine) / sample['sigma_ind']))
    assert worst < 0.1, worst
    print(f'  ✓ rtol 1e-6 and 1e-8 agree within {worst:.1e} standard errors on every share')


# ------------------------------------------------------------------------------
# 5. Structures, optimiser and profiles
# ------------------------------------------------------------------------------
def test_structures_build_models():
    observables = build_synthetic_observables()
    free = build_sample_set(observables, 'free')
    names = {name: [p.name for p in get_structure(name).parameters(free)] for name in
             ('M0-BEP', 'M0-Marcus', 'M1', 'M1-split', 'M3', 'M3-split', 'M3-split+W')}
    assert names['M0-BEP'] == ['g all reactions', 'log10 age 2 % H2O']
    assert names['M1-split'][:5] == ['g hydrolysis_R1', 'g hydrolysis_R23', 'g transfer', 'g condensation', 'g solvent_attack']
    assert names['M3-split+W'][5:] == ['dG R1', 'dG R2', 'dG R3', 'dG R4', 'water fraction', 'log10 age 2 % H2O']
    # M0 reproduces the registered models at their own barrier; M1 reproduces level1
    for structure, registered, g in (('M0-BEP', 'peter_reference', 1.15), ('M0-Marcus', 'level1', 1.15)):
        built = get_structure(structure).build({'g all reactions': g})
        expected = get_model(registered).with_barrier(g).calculate_rates(300.0)
        assert np.array_equal(built['model'].calculate_rates(300.0)['k_f'].values, expected['k_f'].values), structure
    level1_theta = {'g hydrolysis': 0.80, 'g transfer': 0.80, 'g condensation': 1.30, 'g solvent_attack': 1.32}
    built = get_structure('M1').build(level1_theta)
    assert np.array_equal(built['model'].calculate_rates(300.0)['k_f'].values, get_model('level1').calculate_rates(300.0)['k_f'].values)
    assert built['prior'].size == 0 and built['ages'] == {} and built['water_fraction'] == 1.0
    # The split gives R2 and R3 one barrier; freed energies move the ladder and leave R8, R9 alone
    theta = embed_parent_theta('M3-split', level1_theta)
    assert theta['g hydrolysis_R1'] == theta['g hydrolysis_R23'] == 0.80 and theta['dG R2'] == 0.0
    same = get_structure('M3-split').build(theta)
    assert np.allclose(same['model'].calculate_rates(300.0)['k_f'].values, built['model'].calculate_rates(300.0)['k_f'].values,
                       rtol=1e-12), 'a child at its parent\'s values is the parent'
    theta.update({'g hydrolysis_R23': 1.0, 'dG R2': 0.12, 'dG R3': 0.17, 'log10 age 2 % H2O': 2.0, 'water fraction': 0.5})
    moved = get_structure('M3-split+W').build(theta)
    rates, base = moved['model'].calculate_rates(300.0), built['model'].calculate_rates(300.0)
    assert rates.loc['R2', 'g_eV'] == rates.loc['R3', 'g_eV'] == 1.0 and rates.loc['R1', 'g_eV'] == 0.80
    assert abs(rates.loc['R2', 'dG_rxn_eV'] - base.loc['R2', 'dG_rxn_eV'] - 0.12) < 1e-9
    assert abs(rates.loc['R6', 'dG_rxn_eV'] - base.loc['R6', 'dG_rxn_eV'] - 0.12) < 1e-9, 'R6 = R2 + R4 follows'
    assert abs(rates.loc['R8', 'dG_rxn_eV'] - base.loc['R8', 'dG_rxn_eV']) < 1e-9
    assert moved['ages'] == {'2 % H2O': 100.0} and moved['water_fraction'] == 0.5
    cov = calculate_solvation_covariance(['R1', 'R2', 'R3', 'R4']).to_numpy()
    d = np.array([0.0, 0.12, 0.17, 0.0])
    assert abs(moved['prior'] @ moved['prior'] - d @ np.linalg.solve(cov, d)) < 1e-9, 'the constraint is the Mahalanobis distance'
    assert [round(row['z'], 6) for row in moved['prior_rows']] == [round(v, 6) for v in d / np.sqrt(np.diag(cov))]
    barriers = tabulate_reaction_barriers(get_structure('M3-split'), theta)
    assert barriers['R1']['T_K'] == 295.65 and barriers['R8']['T_K'] == 353.15 and len(barriers) == 9
    print('  ✓ structures: parameter lists, registered models reproduced, split and freed energies, constraint')


def test_interval_finder():
    v = np.linspace(-3.0, 3.0, 25)
    parabola = pd.DataFrame({'parameter': 'p', 'value': 1.0 + 0.1 * v, 'chi2': 5.0 + v ** 2})
    parabola['delta'] = parabola['chi2'] - 5.0
    ci = find_confidence_interval(parabola, 'p')
    assert ci['status'] == 'interval' and abs(ci['low'] - (1.0 - 0.196)) < 0.002 and abs(ci['high'] - (1.0 + 0.196)) < 0.002
    # A grid coarser than the interval (step 2.5 σ): the edges must not be pulled towards the best value
    coarse = pd.DataFrame({'parameter': 'p', 'value': 1.0 + 0.1 * np.array([-5.0, -2.5, 0.0, 2.5, 5.0]),
                           'chi2': 5.0 + np.array([-5.0, -2.5, 0.0, 2.5, 5.0]) ** 2})
    coarse['delta'] = coarse['chi2'] - 5.0
    ci = find_confidence_interval(coarse, 'p')
    half = 0.1 * np.sqrt(3.84)
    assert abs(ci['low'] - (1.0 - half)) < 1e-9 and abs(ci['high'] - (1.0 + half)) < 1e-9, ci
    one_sided = parabola.assign(chi2=np.where(v > 0, 5.0, 5.0 + v ** 2))
    ci = find_confidence_interval(one_sided, 'p')
    assert ci['status'] == 'lower bound only' and np.isnan(ci['high']) and abs(ci['low'] - (1.0 - 0.196)) < 0.002
    ci_low = ci['low']
    flat = parabola.assign(chi2=5.0)
    assert find_confidence_interval(flat, 'p')['status'] == 'not determined'
    # Two basins: the interval is their envelope and says so; an accepted local optimum widens a profile
    w = np.linspace(-5.0, 5.0, 41)
    two = pd.DataFrame({'parameter': 'p', 'value': 1.0 + 0.1 * w, 'chi2': 5.0 + np.minimum((w + 2.0) ** 2, (w - 2.0) ** 2 + 1.0)})
    two['delta'] = two['chi2'] - 5.0
    ci = find_confidence_interval(two, 'p')
    assert ci['separate ranges'] and abs(ci['low'] - 0.604) < 0.003 and abs(ci['high'] - 1.3685) < 0.003, ci
    assert not find_confidence_interval(parabola, 'p')['separate ranges']
    from kinetics.fitting import add_profile_points
    narrow = parabola[np.abs(v) <= 3.0].copy()
    widened = add_profile_points(narrow, [{'theta': {'p': 1.5}, 'chi2': 6.0}, {'theta': {'p': 1.7}, 'chi2': 30.0}])
    assert len(widened) == len(narrow) + 1, 'a local optimum outside the accepted range of χ² is not a profile point'
    ci = find_confidence_interval(widened, 'p')
    assert ci['separate ranges'] and np.isnan(ci['high']) and ci['status'] == 'lower bound only', ci
    band = calculate_prediction_band(lambda theta: theta['a'] ** 2, [{'a': 1.0}, {'a': -3.0}, {'a': 2.0}])
    assert band['low'] == 1.0 and band['high'] == 9.0 and band['theta_high'] == {'a': -3.0}
    # An open side survives a result file: JSON has no NaN, and it must come back as one
    import tempfile
    from kinetics.fitting import load_fit_result, write_fit_result
    with tempfile.TemporaryDirectory() as folder:
        path = write_fit_result({'intervals': pd.DataFrame([find_confidence_interval(one_sided, 'p')]), 'n': np.int64(3)},
                                os.path.join(folder, 'result.json'))
        stored = load_fit_result(path)
    row = stored['intervals'].iloc[0]
    assert stored['n'] == 3 and np.isnan(row['high']) and f"{row['high']:.2f}" == 'nan' and abs(row['low'] - ci_low) < 1e-12
    print('  ✓ intervals: ±1.96 σ on a parabola even on a coarse grid, an open side stays open, a flat profile is '
          '"not determined", an open side survives a result file')


def test_fit_recovers_and_reports_blind_spots(workers=4):
    """Synthetic shares from a known parameter set, on the TMSOH samples only: the barriers these samples see
    come back inside their intervals (an open side counts as inside); the hydrolysis barrier, which they cannot
    see, is 'not determined'."""
    observables = build_synthetic_observables()
    tmsoh = build_sample_set(observables, 'middle', exclude=('0.5 % H2O', '2 % H2O', 'TMSPa + TMSOH (A)'))
    truth = {'g hydrolysis': 1.25, 'g transfer': 1.00, 'g condensation': 1.22, 'g solvent_attack': 1.315}
    exact = simulate_synthetic_shares('M1', truth, tmsoh)
    assert np.max(np.abs(FitProblem('M1', exact).residuals([truth[k] for k in truth]))) < 1e-4, 'no noise: no residual'
    data = simulate_synthetic_shares('M1', truth, tmsoh, seed=1)
    problem = FitProblem('M1', data, fixed={'g transfer': 1.00})
    at_truth = problem.chi2(problem.vector(truth))
    fit = fit_structure(problem, n_screen=64, n_starts=4, max_nfev=15, workers=workers)
    assert fit['chi2'] <= at_truth + 1e-6, (fit['chi2'], at_truth)
    assert fit['n_residuals'] == 13 and fit['n_parameters'] == 3 and fit['fits']
    profile = calculate_profile(problem, fit, ['g solvent_attack', 'g condensation', 'g hydrolysis'], max_nfev=10, workers=workers)
    other = {**fit['theta'], 'g hydrolysis': 1.60, 'g condensation': min(fit['theta']['g condensation'] + 0.05, 1.69)}
    extended = extend_profile(problem, profile, [other], ['g solvent_attack'], max_nfev=6, workers=workers)
    assert len(extended) == len(profile), 'extending a profile keeps its grid'
    both = extended.merge(profile, on=['parameter', 'value'], suffixes=('', '_first'))
    assert (both['chi2'] <= both['chi2_first'] + 1e-9).all(), 'a second start can only lower a profile'
    for name in ('g solvent_attack', 'g condensation'):
        ci = find_confidence_interval(profile, name)
        assert not np.isfinite(ci['low']) or ci['low'] - 0.002 <= truth[name], ci
        assert not np.isfinite(ci['high']) or truth[name] <= ci['high'] + 0.002, ci
    assert find_confidence_interval(profile, 'g solvent_attack')['status'] == 'interval'
    assert find_confidence_interval(profile, 'g hydrolysis')['status'] == 'not determined'
    assert profile[profile['parameter'] == 'g hydrolysis']['delta'].abs().max() < 0.5
    directions = calculate_parameter_directions(problem, fit['theta'])
    assert np.isinf(directions['sigma_local']['g hydrolysis']), 'a parameter the data do not see has no standard error'
    assert directions['sigma_local']['g solvent_attack'] < 0.02
    table = compare_structures({'M1': fit})
    assert bool(table.loc['M1', 'fits']) and table.loc['M1', 'residuals'] == 13
    # Without the glovebox sample nothing bounds solvent attack from above
    without = {'TMSOH, glovebox': [s for s in data if s['name'] != 'TMSOH, glovebox']}
    loo = fit_leave_one_out('M1', without, fit, n_screen=32, n_starts=2, max_nfev=10, workers=workers)
    assert loo['TMSOH, glovebox']['n_residuals'] == 10
    print(f"  ✓ recovery: χ² {fit['chi2']:.1f} ≤ {at_truth:.1f} at the truth; truth inside the intervals; "
          "the unseen barrier is 'not determined'")


def test_nesting_never_raises_chi2(workers=4):
    observables = build_synthetic_observables()
    samples = build_sample_set(observables, 'middle', exclude=('2 % H2O', 'TMSOH, probe', 'E1'))
    parent = fit_structure(FitProblem('M1', samples), n_screen=64, n_starts=4, max_nfev=15, workers=workers, report=False)
    child = fit_structure(FitProblem('M1-split', samples), n_screen=0, max_nfev=15, workers=workers, report=False,
                          start_points=[embed_parent_theta('M1-split', parent['theta'])])
    assert child['chi2'] <= parent['chi2'] + 1e-6, (child['chi2'], parent['chi2'])
    freed = fit_structure(FitProblem('M3-split', samples), n_screen=0, max_nfev=15, workers=workers, report=False,
                          start_points=[embed_parent_theta('M3-split', child['theta'])])
    assert freed['chi2'] <= child['chi2'] + 1e-6, (freed['chi2'], child['chi2'])
    print(f"  ✓ nesting: χ² {parent['chi2']:.1f} (M1) ≥ {child['chi2']:.1f} (M1-split) ≥ {freed['chi2']:.1f} (M3-split)")


# ------------------------------------------------------------------------------
# 6. Notebook 03's intervals through the new forward path (needs the lab file)
# ------------------------------------------------------------------------------
GRID_EV = np.round(np.arange(0.60, 1.701, 0.025), 3)        # the scan grid of notebook 03


def _violation(c0, T_K, t_min_s, t_max_s, windows, model):
    """Smallest distance outside the windows over 25 log-spaced times, as kinetics.fitting.evaluate_observation,
    but through simulate_history and the share functions of the fit."""
    single = t_max_s <= t_min_s * (1.0 + 1e-9)
    times = np.array([t_max_s]) if single else np.geomspace(max(t_min_s, 1.0), t_max_s, 25)
    sim = simulate_history(c0, [(T_K, t_max_s)], times, model=model, max_rhs_calls=2_000_000)
    assert sim['success'], sim['message']
    C, idx = np.maximum(sim['C_M'], 0.0), sim['idx']
    values = {
        'tmspa_share': lambda: _shares_from_state(C, idx, '31P')[:, 0],
        'ring_opened': lambda: _shares_from_state(C, idx, '13C')[:, 2],
        'hmdso': lambda: _shares_from_state(C, idx, '13C')[:, 1],
        'tmspa_conversion': lambda: 1.0 - C[idx['TMSPA']] / c0['TMSPA'],
    }
    out = np.zeros(len(times))
    for name, (low, high) in windows.items():
        v = values[name]()
        out += (np.maximum(0.0, low - v) if low is not None else 0.0) + (np.maximum(0.0, v - high) if high is not None else 0.0)
    return float(out.min())


def _interval(violation, tol_eV=0.002):
    """(g_low, g_high) allowed along the scan grid, edges refined by bisection; NaN = open at the grid limit.

    As kinetics.fitting.find_allowed_intervals: when no grid point is allowed, a window narrower than the grid
    step is looked for next to the point that comes closest.
    """
    from scipy.optimize import minimize_scalar
    grid = GRID_EV
    v = np.array([violation(g) for g in grid])
    ok = v <= 0.0
    if not ok.any():
        i = int(np.argmin(v))
        res = minimize_scalar(violation, bounds=(grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]),
                              method='bounded', options={'xatol': tol_eV})
        if res.fun > 0.0:
            return None
        order = np.argsort(np.append(grid, res.x))
        grid, ok = np.append(grid, res.x)[order], np.append(ok, True)[order]
    first, last = np.flatnonzero(ok)[[0, -1]]

    def edge(g_out, g_in):
        while abs(g_in - g_out) > tol_eV:
            g_mid = 0.5 * (g_in + g_out)
            g_in, g_out = (g_mid, g_out) if violation(g_mid) <= 0.0 else (g_in, g_mid)
        return 0.5 * (g_in + g_out)

    return (edge(grid[first - 1], grid[first]) if first > 0 else np.nan,
            edge(grid[last + 1], grid[last]) if last < len(grid) - 1 else np.nan)


def _notebook_03_observations(observables):
    """L1–L5, E1 and E4 as (c0, T_K, t_min_s, t_max_s, windows), built from the observables file as in notebook 03 §2.2."""
    S = observables['samples']

    def window(sample, spectrum, name, side, k=2.0):
        v = S[sample]['spectra'][spectrum]['shares'][name]
        sigma = np.hypot(v['sigma_noise'], v['sigma_baseline'])
        return (max(0.0, v['share'] - k * sigma) if side in ('at least', 'both') else None,
                min(1.0, v['share'] + k * sigma) if side in ('at most', 'both') else None)

    def lab(sample, spectrum, name, side):
        sp = S[sample]['spectra'][spectrum]
        t = 3600.0 * sp['t_since_first_h']
        return (S[sample]['c0_M'], sp['temperature_C'] + ZERO_CELSIUS_K, t, t + 86400.0,
                {'tmspa_share': window(sample, spectrum, name, side)})

    heated = S['TMSOH, probe']['heated'][0]
    controls = {e['id']: e for e in load_experimental_data()['control_experiments']}
    return {
        'L1': lab('0.5 % H2O', 0, 'TMSPA', 'at least'),
        'L2': lab('2 % H2O', 0, 'TMSPA', 'at most'),
        'L3': lab('TMSPa + TMSOH (A)', 0, 'TMSPA', 'at most'),
        'L4': (S['TMSOH, probe']['c0_M'], 353.15, 3600.0 * heated['duration_h'], 5.6 * 86400.0,
               {'ring_opened': window('TMSOH, probe', -1, 'TMSOEG', 'at most')}),
        'L5': (S['TMSOH, glovebox']['c0_M'], 353.15, 8 * 3600.0, 8 * 3600.0,
               {'ring_opened': window('TMSOH, glovebox', 0, 'TMSOEG', 'both'),
                'hmdso': window('TMSOH, glovebox', 0, 'HMDSO', 'at most')}),
        'E1': (controls['E1']['c0_M'], 298.15, 604800.0, 604800.0, {'ring_opened': (None, 0.02)}),
        'E4': (controls['E4']['c0_M'], 298.15, 86400.0, 86400.0, {'tmspa_conversion': (0.95, None)}),
    }


def test_windows_reproduce_notebook_03_if_available(tolerance_eV=0.005):
    if not os.path.exists(DEFAULT_OBSERVABLES_PATH):
        print('  - skipped: lab_observables.json not found in the data folder (set ATOM_DATA_DIR)')
        return None
    observations = _notebook_03_observations(load_lab_observables())
    level1, peter = get_model('level1'), get_model('peter_reference')
    structures = {
        'hydrolysis': lambda g: level1.with_family_params('hydrolysis', g_eV=g),
        'transfer': lambda g: level1.with_family_params('transfer', g_eV=g),
        'condensation': lambda g: level1.with_family_params('condensation', g_eV=g),
        'solvent_attack': lambda g: level1.with_family_params('solvent_attack', g_eV=g),
        'M0 · capped BEP, one E0': peter.with_barrier,
        'M0 · Marcus, one g': level1.with_barrier,
    }
    expected = pd.concat([pd.read_csv(os.path.join(RESULTS_03, name), index_col=0)
                          for name in ('m1_intervals.csv', 'm0_intervals.csv')])
    rows = []
    for structure, grp in expected.groupby(level=0, sort=False):
        for _, row in grp.iterrows():
            obs = observations[row['observation']]
            found = _interval(lambda g: _violation(*obs, structures[structure](float(g))))
            new_low, new_high = (np.nan, np.nan) if found is None else found
            none_expected = str(row['verdict']).startswith('no value')
            rows.append({'structure': structure, 'observation': row['observation'],
                         'nb03 low': row['g_low_eV'], 'new low': new_low, 'nb03 high': row['g_high_eV'], 'new high': new_high,
                         'ok': (found is None) == none_expected and all(
                             (np.isnan(a) and np.isnan(b)) or abs(a - b) <= tolerance_eV
                             for a, b in ((row['g_low_eV'], new_low), (row['g_high_eV'], new_high)))})
    table = pd.DataFrame(rows)
    worst = np.nanmax(np.abs(np.concatenate([table['nb03 low'] - table['new low'], table['nb03 high'] - table['new high']])))
    print(table.round(4).to_string(index=False))
    assert table['ok'].all(), table[~table['ok']]
    print(f'  ✓ notebook 03: {len(table)} intervals reproduced through simulate_history, largest difference {worst:.4f} eV')
    return table


if __name__ == "__main__":
    test_history_matches_existing_reactors()
    test_model_with_own_network()
    test_species_shifts()
    test_sample_set_and_histories()
    test_replicate_block_and_limits()
    test_free_ages_and_tables()
    test_search_and_report_tolerances_agree()
    test_structures_build_models()
    test_interval_finder()
    test_fit_recovers_and_reports_blind_spots()
    test_nesting_never_raises_chi2()
    test_windows_reproduce_notebook_03_if_available()
    print("All estimation tests passed.")
