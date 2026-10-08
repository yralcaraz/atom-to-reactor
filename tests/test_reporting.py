"""Tests for the reporting of a stored fit: shares per spectrum, trajectories, window amounts and correlations.

Covers:
1. evaluate_fit and tabulate_fit_shares: what a stored parameter set predicts, one row per spectrum and window.
2. simulate_fit_trajectories and tabulate_window_amounts: the trajectory passes through the predicted shares,
   phosphorus and silyl groups are conserved, the history of each tube is followed and can be run on.
3. Free ages: the age of an unheated tube is found again on the trajectory itself; given ages are used as they are.
4. calculate_parameter_correlation: a correlation matrix with a finite standard error for what the samples see.

Source: Y. Alcaraz Galván (synthetic numbers, no lab data)
"""

import os
import sys

import numpy as np

repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for path in (repo_root, os.path.dirname(os.path.abspath(__file__))):
    if path not in sys.path:
        sys.path.insert(0, path)

from kinetics.fitting.reporting import (
    calculate_parameter_correlation, evaluate_fit, simulate_fit_trajectories, tabulate_fit_shares, tabulate_window_amounts,
)
from kinetics.fitting.residuals import REPORT_SOLVER, _whiten_sample, calculate_predicted_shares
from synthetic_observables import build_synthetic_observables

THETA = {'g hydrolysis_R1': 1.28, 'g hydrolysis_R23': 1.05, 'g transfer': 1.10, 'g condensation': 1.25,
         'g solvent_attack': 1.315, 'c0 H2O TMSPa alone': 0.06}
FIT = {'structure': 'M1-split', 'scenario': 'middle', 'theta': THETA, 'ages': {}}


def _chi2(sample, evaluation, age_h):
    predicted = calculate_predicted_shares(sample, evaluation['built']['model'], age_h=age_h, solver=REPORT_SOLVER,
                                           added_M=evaluation['built']['added_M'].get(sample['name']))
    return float(np.sum(_whiten_sample(sample, predicted) ** 2))


def test_shares_of_a_stored_fit():
    evaluation = evaluate_fit(FIT, build_synthetic_observables())
    shares = tabulate_fit_shares(evaluation)
    lab = [s for s in evaluation['samples'] if s['kind'] == 'lab']
    assert len(shares) == sum(block['measured'].size for s in lab for block in s['blocks'])
    assert np.allclose(shares.groupby(['sample', 'nucleus', 'spectrum'])['predicted'].sum(), 1.0)
    for sample in lab:
        direct = calculate_predicted_shares(sample, evaluation['built']['model'], age_h=evaluation['ages'][sample['name']],
                                            added_M=evaluation['built']['added_M'].get(sample['name']), solver=REPORT_SOLVER)
        rows = shares[shares['sample'] == sample['name']]
        assert np.allclose(rows['predicted'].to_numpy(), np.concatenate([p.ravel() for p in direct]), atol=1e-12)
        assert np.allclose(rows['measured'].to_numpy(), np.concatenate([b['measured'].ravel() for b in sample['blocks']]))
    assert evaluation['ages']['2 % H2O'] == 24.0 and evaluation['ages']['TMSOH, glovebox'] is None
    totals = shares.groupby('sample')['total_M'].first()
    assert np.isclose(totals['TMSOH, probe'], 0.45) and np.isclose(totals['2 % H2O'], evaluation['samples'][1]['c0_M']['TMSPA'])
    print(f'  ✓ shares of a stored fit: {len(shares)} rows, one per spectrum and window, equal to the direct prediction')


def test_trajectories_and_window_amounts():
    evaluation = evaluate_fit(FIT, build_synthetic_observables())
    shares = tabulate_fit_shares(evaluation)
    trajectories = simulate_fit_trajectories(evaluation)
    phosphate, silyl = tabulate_window_amounts(trajectories, '31P'), tabulate_window_amounts(trajectories, '13C')
    for name in ('0.5 % H2O', '2 % H2O', 'TMSPa alone', 'TMSOH, probe', 'TMSOH, glovebox'):
        rows = shares[shares['sample'] == name]
        last = rows[rows['t_h'] == rows['t_h'].max()]
        table = phosphate if rows['nucleus'].iloc[0] == '31P' else silyl
        end = table[table['sample'] == name].iloc[-1]
        assert np.isclose(end['t_h'], last['t_h'].iloc[0]), 'a trajectory ends at the last spectrum'
        amounts = np.array([end[w] for w in last['window']])
        assert np.allclose(amounts / amounts.sum(), last['predicted'].to_numpy(), atol=1e-8), name
        assert np.allclose(amounts.sum(), last['total_M'].iloc[0], rtol=1e-6), 'the windows count every nucleus of the tube'
    # Phosphorus and silyl groups are conserved along a trajectory; HMDSO carries two silyl groups
    tube = trajectories['sample'] == '2 % H2O'
    tmspa0 = next(s for s in evaluation['samples'] if s['name'] == '2 % H2O')['c0_M']['TMSPA']
    assert np.allclose(phosphate.loc[tube, ['TMSPA', 'BMSPA', 'MMSPA', 'H3PO4']].sum(axis=1), tmspa0, rtol=1e-6)
    assert np.allclose(silyl.loc[tube].drop(columns=['sample', 't_h', 'T_C']).sum(axis=1), 3.0 * tmspa0, rtol=1e-6)
    assert np.allclose(silyl['HMDSO'], 2.0 * trajectories['HMDSO'].clip(lower=0.0))
    # Each tube follows its own history, and can be run on at room temperature
    box = trajectories[trajectories['sample'] == 'TMSOH, glovebox']
    assert np.isclose(box['T_C'].max(), 80.0) and np.isclose(box.loc[box['T_C'] > 25.0, 't_h'].max(), 8.0)
    assert trajectories.loc[tube, 'T_C'].max() > 25.0, 'the heating steps of the 2 % tube are in its trajectory'
    longer = simulate_fit_trajectories(evaluation, until_h=500.0)
    assert np.isclose(longer.groupby('sample')['t_h'].max(), 500.0).all()
    assert np.isclose(longer[longer['sample'] == 'TMSOH, glovebox'].iloc[-1]['T_C'], evaluation['samples'][-2]['T_rt_K'] - 273.15)
    print('  ✓ trajectories: through the predicted shares, P and silyl groups conserved, histories followed and run on')


def test_free_ages_are_found_on_the_trajectory():
    observables = build_synthetic_observables()
    stored = {'0.5 % H2O': 30.0, 'TMSPa + TMSOH (A)': 5.0, 'TMSPa alone': 60.0}
    fit = {**FIT, 'scenario': 'free', 'theta': {**THETA, 'log10 age 2 % H2O': 2.0}, 'ages': dict(stored)}
    evaluation = evaluate_fit(fit, observables)
    assert np.isclose(evaluation['ages']['2 % H2O'], 100.0) and evaluation['ages']['TMSOH, probe'] == 1.0
    step = 720.0 ** (1.0 / 52.0)
    for sample in evaluation['samples']:
        name = sample['name']
        if name not in stored:
            continue
        found = evaluation['ages'][name]
        assert stored[name] / step - 1e-9 <= found <= stored[name] * step + 1e-9, 'inside one step of the stored age'
        scan = [_chi2(sample, evaluation, age) for age in np.geomspace(stored[name] / step, stored[name] * step, 25)]
        assert _chi2(sample, evaluation, found) <= min(scan) + 1e-3, (name, found)
    kept = evaluate_fit(fit, observables, ages={**stored, '2 % H2O': 100.0})
    assert all(kept['ages'][name] == age for name, age in stored.items()), 'given ages are used as they are'
    print('  ✓ free ages: found again on the trajectory within one step of the stored age; given ages are kept')


def test_parameter_correlation():
    found = calculate_parameter_correlation(FIT, build_synthetic_observables())
    matrix, sigma = found['correlation'].to_numpy(), found['sigma_local']
    assert list(found['correlation'].index) == list(THETA)
    seen = np.isfinite(sigma.to_numpy())
    assert seen[list(THETA).index('g solvent_attack')] and sigma['g solvent_attack'] < 0.05
    block = matrix[np.ix_(seen, seen)]
    assert np.allclose(block, block.T) and np.allclose(np.diag(block), 1.0) and np.all(np.abs(block) <= 1.0 + 1e-9)
    assert np.isnan(matrix[~seen]).all(), 'a parameter the samples do not see has no correlation'
    print(f"  ✓ correlation: {int(seen.sum())} of {len(seen)} parameters seen; solvent attack to ±{sigma['g solvent_attack']:.3f} eV")


if __name__ == "__main__":
    test_shares_of_a_stored_fit()
    test_trajectories_and_window_amounts()
    test_free_ages_are_found_on_the_trajectory()
    test_parameter_correlation()
    print("All reporting tests passed.")
