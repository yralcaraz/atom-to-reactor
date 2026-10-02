"""Lab NMR reader tests.

Covers:
1. FID → spectrum on a synthetic signal: peak position (Delta axis convention), digital-filter delay, phasing.
2. The real data directory, when it is available: every file parses and FIDs match the exported spectrum.
3. Area shares with their uncertainty on a synthetic spectrum, and heated intervals from acquisition times.
"""

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

from kinetics.data.lab_nmr import (
    DEFAULT_LAB_NMR_DIR, build_lab_nmr_inventory, calculate_area_shares, calculate_spectrum, find_heated_windows,
    load_spectrum, load_text_spectrum, phase_spectrum,
)


def _synthetic_fid(peaks_ppm, *, n=8192, sweep_Hz=20000.0, freq_MHz=161.8, offset_ppm=0.0, delay=20,
                   phase_rad=1.1):
    """Delta convention: a component at +f Hz sits at offset_ppm − f/freq_MHz."""
    t = np.arange(n) / sweep_Hz
    fid = sum(np.exp(-2j * np.pi * (ppm - offset_ppm) * freq_MHz * t) for ppm in peaks_ppm)
    fid = fid * np.exp(-t / 0.2) * np.exp(1j * phase_rad)
    rng = np.random.default_rng(0)
    fid = fid + 0.002 * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    return np.concatenate([np.zeros(delay), fid[:-delay]]) + 0j


def test_synthetic_fid():
    peaks = [-24.6, -6.4]
    fid = _synthetic_fid(peaks)
    ppm, spectrum = calculate_spectrum(fid, sweep_Hz=20000.0, freq_MHz=161.8, offset_ppm=0.0)
    assert ppm[0] > ppm[-1], 'ppm axis must be descending'
    real, _ = phase_spectrum(spectrum)
    for target in peaks:
        window = np.abs(ppm - target) < 1.0
        assert abs(ppm[window][np.argmax(real[window])] - target) < 0.02, f'peak at {target} ppm misplaced'
    assert real.max() > 20 * abs(real.min()), 'phased spectrum should be absorptive (no deep negative lobes)'
    print('  ✓ synthetic FID: peaks within 0.02 ppm, absorptive after automatic phasing')


def test_real_data_if_available():
    root = Path(DEFAULT_LAB_NMR_DIR)
    if not root.is_dir():
        print(f'  - skipped: {root} not found (set LAB_NMR_DIR)')
        return
    inventory = build_lab_nmr_inventory(root)
    assert len(inventory) > 0 and inventory['acquired_at'].notna().all()
    folder = root / '5% TMSPa + 0.5v% H2O in EC-DEC'
    ppm, real = load_spectrum(folder / 'WW-NG-S1 TMSPa+H2O in EC-DEC_31P-1-1.jdf')
    ppm_ref, real_ref = load_text_spectrum(folder / 'WW-NG-S1 TMSPa+H2O in EC-DEC_31P-1.asc')
    window, window_ref = np.abs(ppm + 24.6) < 1.0, np.abs(ppm_ref + 24.6) < 1.0
    top = ppm[window][np.argmax(real[window])]
    top_ref = ppm_ref[window_ref][np.argmax(real_ref[window_ref])]
    assert abs(top - top_ref) < 0.03, f'TMSPA peak {top:.3f} vs export {top_ref:.3f} ppm'
    print(f'  ✓ {len(inventory)} files parsed; ³¹P TMSPA peak {top:.2f} ppm vs exported {top_ref:.2f} ppm')


def test_area_shares_on_synthetic_spectrum():
    ppm = np.linspace(6.0, -32.0, 8000)
    def lorentzian(center, area, fwhm=0.3):
        return area * (fwhm / 2) / np.pi / ((ppm - center) ** 2 + (fwhm / 2) ** 2)
    rng = np.random.default_rng(1)
    baseline = 0.02 + 0.001 * ppm                      # sloped baseline
    y = lorentzian(-24.6, 3.0) + lorentzian(-14.0, 1.0) + baseline + 0.002 * rng.standard_normal(len(ppm))
    windows = {'TMSPA': (-27.0, -22.0), 'BMSPA': (-16.5, -12.5), 'MMSPA': (-8.5, -5.0)}
    shares = calculate_area_shares(ppm, y, windows, region_ppm=(-32.0, 6.0))
    # Lorentzian tails outside the windows are lost equally from both peaks, so the ratio survives
    assert abs(shares.loc['TMSPA', 'share'] - 0.75) < 0.01, shares
    assert abs(shares.loc['BMSPA', 'share'] - 0.25) < 0.01, shares
    assert abs(shares.loc['MMSPA', 'share']) < 3 * shares.loc['MMSPA', 'sigma'] + 0.005, 'empty window must be ~0'
    assert (shares['sigma'] >= shares['sigma_noise']).all() and (shares['sigma'] > 0).all()
    print('  ✓ area shares: 3:1 peaks recovered within 0.01 on a sloped, noisy baseline; empty window ≈ 0')


def test_heated_windows_from_acquisition_times():
    t0 = pd.Timestamp('2022-06-23 10:00')
    spectra = pd.DataFrame({
        'acquired_at': [t0, t0 + pd.Timedelta(minutes=10), t0 + pd.Timedelta(hours=5), t0 + pd.Timedelta(hours=6),
                        t0 + pd.Timedelta(hours=7)],
        'duration_min': [3.0, 6.0, 3.0, 6.0, 3.0],
        'temperature_C': [30.0, 30.0, 22.5, 40.0, 22.4],
    })
    heated = find_heated_windows(spectra)
    assert list(heated['T_C']) == [30.0, 40.0] and list(heated['n_spectra']) == [2, 1]
    assert abs(heated.loc[0, 'minutes'] - 16.0) < 1e-9 and abs(heated.loc[1, 'minutes'] - 6.0) < 1e-9
    print('  ✓ heated windows: consecutive spectra at one temperature merge; room-temperature spectra split them')


if __name__ == "__main__":
    print("================================================================================")
    print("RUNNING LAB NMR READER TESTS")
    print("================================================================================")
    test_synthetic_fid()
    test_real_data_if_available()
    test_area_shares_on_synthetic_spectrum()
    test_heated_windows_from_acquisition_times()
    print("\n================================================================================")
    print("ALL LAB NMR READER TESTS PASSED")
    print("================================================================================")
