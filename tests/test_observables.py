"""Tests for the curated lab observables (kinetics/data/observables.py).

Covers:
1. Building the observables from share tables, the round trip through JSON and the display table.
2. Heating episodes: transcribed ones by default, `find_heated_windows` tables when given.
3. Replicate scatter against the stated noise.
4. The real file, when it is present: it reproduces the share tables exported by notebook 03 exactly.

Source: Y. Alcaraz Galván
"""

import os
import sys
import tempfile

import numpy as np
import pandas as pd

repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for path in (repo_root, os.path.dirname(os.path.abspath(__file__))):
    if path not in sys.path:
        sys.path.insert(0, path)

from kinetics.data import (
    DEFAULT_OBSERVABLES_PATH, build_lab_observables, calculate_replicate_scatter, describe_lab_observables,
    load_lab_observables, load_share_tables, tabulate_observable_shares, write_lab_observables,
)
from kinetics.data.observables import DEFAULT_SHARES_DIR, LAB_SAMPLES, RECORDED_HISTORIES
from synthetic_observables import SYNTHETIC_SPECTRA, build_synthetic_observables, build_synthetic_share_tables


def test_build_and_round_trip():
    observables = build_synthetic_observables()
    assert set(observables['samples']) == set(LAB_SAMPLES)
    with tempfile.TemporaryDirectory() as tmp:
        path = write_lab_observables(observables, os.path.join(tmp, 'observables.json'))
        assert load_lab_observables(path) == observables, 'JSON round trip must not change the content'
    shares_P, shares_C = build_synthetic_share_tables()
    table = tabulate_observable_shares(observables)
    assert len(table) == len(shares_P) + len(shares_C)
    for name, sample in observables['samples'].items():
        assert len(sample['spectra']) == len(SYNTHETIC_SPECTRA[name])
        for spectrum, (hours, T_C, scans, shares) in zip(sample['spectra'], SYNTHETIC_SPECTRA[name]):
            assert spectrum['t_since_first_h'] == hours and spectrum['scans'] == scans
            assert [v['share'] for v in spectrum['shares'].values()] == list(shares)
        assert sample['c0_M']['EC'] == 7.1 and all(c > 0 for c in sample['c0_M'].values())
    assert observables['samples']['TMSOH, glovebox']['age']['known'], 'the glovebox hold is the known pre-history'
    assert not observables['samples']['2 % H2O']['age']['known']
    described = describe_lab_observables(observables)
    assert described.loc['2 % H2O', 'spectra'] == '6 31P' and described.loc['TMSOH, probe', 'role'] == 'fit'
    print('  ✓ observables: built from share tables, JSON round trip exact, 7 samples described')


def test_heating_episodes():
    observables = build_synthetic_observables()
    ramp = observables['samples']['2 % H2O']['heated']
    assert [h['T_C'] for h in ramp] == [30.0, 40.0, 50.0, 60.0, 70.0, 80.0]
    assert abs(sum(h['duration_h'] for h in ramp) - 1.8) < 1e-9, 'notebook 03: 1.8 h of heating in all'
    assert ramp == RECORDED_HISTORIES['2 % H2O']['heated']
    # A table of find_heated_windows replaces the transcribed episodes
    shares_P, shares_C = build_synthetic_share_tables()
    first = shares_P.loc[shares_P['sample'] == '2 % H2O', 'acquired_at'].min() - pd.Timedelta(hours=0.1)
    table = pd.DataFrame({'start': [first + pd.Timedelta(hours=10)], 'end': [first + pd.Timedelta(hours=10.5)],
                          'T_C': [60.0], 'n_spectra': [3], 'minutes': [30.0]})
    rebuilt = build_lab_observables(shares_P, shares_C, built_from='test', heated_windows={'2 % H2O': table})
    heated = rebuilt['samples']['2 % H2O']['heated']
    assert len(heated) == 1 and abs(heated[0]['start_h'] - 10.0) < 1e-9 and abs(heated[0]['duration_h'] - 0.5) < 1e-9
    assert rebuilt['samples']['TMSOH, probe']['heated'] == RECORDED_HISTORIES['TMSOH, probe']['heated']
    print('  ✓ heating episodes: transcribed by default, replaced by find_heated_windows tables when given')


def test_replicate_scatter():
    scatter = calculate_replicate_scatter(build_synthetic_observables(), '2 % H2O')
    assert scatter['chi2'] < 1e-20 and scatter['n_spectra'] == 6 and scatter['dof'] == 15
    rng = np.random.default_rng(3)
    spectra = {k: list(v) for k, v in SYNTHETIC_SPECTRA.items()}
    ratios = []
    for _ in range(200):
        noisy = []
        for hours, T_C, scans, shares in SYNTHETIC_SPECTRA['2 % H2O']:
            sigma = 0.010 if scans == 16 else 0.003
            noisy.append((hours, T_C, scans, tuple(s + rng.normal(0.0, sigma) for s in shares)))
        spectra['2 % H2O'] = noisy
        ratios.append(calculate_replicate_scatter(build_synthetic_observables(spectra), '2 % H2O')['birge'])
    # independent noise on four windows over 15 degrees of freedom: sqrt(chi2 / dof) with chi2 ~ 4 × chi2(5)
    assert abs(np.mean(ratios) - np.sqrt(20.0 / 15.0)) < 0.05, np.mean(ratios)
    print('  ✓ replicate scatter: zero for identical spectra; Birge ratio as expected for pure noise')


def test_real_file_reproduces_share_tables_if_available():
    csv = os.path.join(DEFAULT_SHARES_DIR, 'lab_shares_31P.csv')
    if not (os.path.exists(DEFAULT_OBSERVABLES_PATH) and os.path.exists(csv)):
        print('  - skipped: lab_observables.json or the notebook 03 share tables not found in the data folder (set ATOM_DATA_DIR)')
        return
    observables = load_lab_observables()
    table = tabulate_observable_shares(observables)
    source = pd.concat(load_share_tables(), ignore_index=True)
    # The tables of notebook 03 hold one nucleus per sample: the spectra of a second nucleus (raw route) are extra rows
    extra = sum(len(spectrum['shares']) for sample in observables['samples'].values()
                for more in sample.get('more_nuclei', {}).values() for spectrum in more['spectra'])
    assert len(source) == 56 and len(table) == 56 + extra, (len(table), len(source), extra)
    merged = source.assign(acquired_at=source['acquired_at'].map(lambda t: pd.Timestamp(t).isoformat())).merge(
        table, on=['sample', 'acquired_at', 'window'], suffixes=('_csv', ''), validate='one_to_one')
    assert len(merged) == 56
    # A file built from the tables holds their numbers; one built from the raw spectra recomputes them (rounding)
    tolerance = 1e-12 if observables['_meta']['built_from'].startswith('raw spectra') else 0.0
    for column in ('share', 'sigma_noise', 'sigma_baseline'):
        assert (merged[column] - merged[f'{column}_csv']).abs().max() <= tolerance, f'{column} differs from the CSV'
    assert (merged['t_since_first_h'] - merged['hours_since_first']).abs().max() <= tolerance
    scatter = calculate_replicate_scatter(observables, '2 % H2O')
    print(f"  ✓ real file: 56 shares of 15 spectra {'identical' if tolerance == 0.0 else 'equal within 1e-12'} to the "
          f"notebook 03 tables, {extra} more shares of a second nucleus "
          f"(built from: {observables['_meta']['built_from'][:40]}…)")
    print(f"    replicate scatter of the six 2 % spectra: chi2 {scatter['chi2']:.1f} on {scatter['dof']} dof, "
          f"Birge ratio {scatter['birge']:.2f}")
    print(scatter['per_window'].round(4).to_string())


if __name__ == "__main__":
    test_build_and_round_trip()
    test_heating_episodes()
    test_replicate_scatter()
    test_real_file_reproduces_share_tables_if_available()
    print("All observables tests passed.")
