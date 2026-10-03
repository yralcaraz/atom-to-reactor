"""Measured NMR spectra from the lab (JEOL Delta .jdf files and their text exports).

The raw files are private and are not stored in the repository. They are read from the directory given by
the environment variable LAB_NMR_DIR, or from the thesis OneDrive copy when it is not set.

    read_jdf                 header, acquisition parameters and data of one .jdf file
    build_lab_nmr_inventory  one row per .jdf file with the acquisition metadata (no spectral processing)
    calculate_spectrum       FID → spectrum: digital-filter delay, exponential broadening, FFT, ppm axis
    phase_spectrum           automatic zero- and first-order phasing (ACME entropy minimisation)
    refine_window_phase      zero-order phase touch-up for one ppm window
    load_spectrum            ppm axis and real spectrum of one .jdf file, processed or FID
    load_text_spectrum       ppm axis and intensity of a Delta text export (.txt or .asc)
    calculate_window_integrals  area, height and position of the signal in named ppm windows
    calculate_area_shares    share of the total area in each window, with a noise and a baseline uncertainty
    tabulate_area_shares     area shares of several spectra, one row per (spectrum, window)
    find_heated_windows      clock intervals in which a sample was recorded above room temperature

Chemical-shift axes are as acquired (lock-based, not referenced to an internal standard).

Classification: CONFIDENTIAL (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván; reads N. Gogoi, raw lab NMR spectra 2022–2023 (unpublished)
"""

import hashlib
import os
import re
import struct
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

DEFAULT_LAB_NMR_DIR = os.environ.get('LAB_NMR_DIR', str(
    Path.home() / 'Library/CloudStorage/OneDrive-Uppsalauniversitet/TFM - Master Thesis'
    / '03 - Validation & writing/VAL - Experimental validation/VAL - NMR raw data - Gogoi 2023'))

# Delta stores times as seconds since this epoch
_JEOL_EPOCH = datetime(1990, 1, 1)
_UNIT_CODES = {13: 'Hz', 26: 'ppm', 28: 's'}
_NUCLEUS_LABELS = {'Proton': '1H', 'Carbon13': '13C', 'Phosphorus31': '31P', 'Silicon29': '29Si',
                   'Fluorine19': '19F'}
# Nucleus implied by the experiment token at the end of a file name
_NAME_TOKENS = {'proton': '1H', 'carbon': '13C', '31p': '31P', '31pcpd': '31P', 'phosphorus31': '31P',
                '29si': '29Si', 'silicon': '29Si', 'inept_dec': '29Si', '19f': '19F'}
# Default exponential line broadening per nucleus [Hz]
LINE_BROADENING_HZ = {'1H': 0.3, '13C': 1.0, '31P': 1.0, '29Si': 2.0, '19F': 1.0}


def _cstr(raw: bytes) -> str:
    return raw.split(b'\0')[0].decode(errors='replace').strip()


def read_jdf(path) -> dict:
    """One 1D .jdf file: 'title', 'nucleus', 'unit' ('s' for an FID, 'ppm' for a processed spectrum),
    'axis_start', 'axis_stop', 'data' (complex, trimmed to the stored region) and 'params' (lower-case names)."""
    raw = Path(path).read_bytes()
    if raw[:8] != b'JEOL.NMR':
        raise ValueError(f'not a JEOL Delta file: {path}')
    n_points = struct.unpack('>I', raw[176:180])[0]
    i_start = struct.unpack('>I', raw[208:212])[0]
    i_stop = struct.unpack('>I', raw[240:244])[0]
    param_start = struct.unpack('>I', raw[1212:1216])[0]
    data_start = struct.unpack('>I', raw[1284:1288])[0]
    data_length = struct.unpack('>Q', raw[1288:1296])[0]

    # Parameter section: 16-byte header, then 64-byte little-endian records
    params = {}
    n_records = struct.unpack('<I', raw[param_start + 8:param_start + 12])[0] + 1
    for i in range(n_records):
        rec = raw[param_start + 16 + 64 * i: param_start + 80 + 64 * i]
        if len(rec) < 64:
            break
        value_type = struct.unpack('<I', rec[32:36])[0]
        name = _cstr(rec[36:64]).lower()
        if not name:
            continue
        if value_type == 0:
            params[name] = _cstr(rec[16:32])
        elif value_type == 1:
            params[name] = struct.unpack('<i', rec[16:20])[0]
        elif value_type == 2:
            params[name] = struct.unpack('<d', rec[16:24])[0]

    dtype = '<f8' if raw[8] == 1 else '>f8'
    values = np.frombuffer(raw[data_start:data_start + data_length], dtype=dtype)
    data = values[:n_points] + 1j * values[n_points:2 * n_points] if raw[24] == 3 else values[:n_points]
    nucleus = _cstr(raw[808:840])
    return {
        'title': _cstr(raw[48:172]),
        'nucleus': _NUCLEUS_LABELS.get(nucleus, nucleus),
        'unit': _UNIT_CODES.get(raw[33], str(raw[33])),
        'axis_start': struct.unpack('>d', raw[272:280])[0],
        'axis_stop': struct.unpack('>d', raw[336:344])[0],
        'data': data[i_start:i_stop + 1],
        'params': params,
    }


def _acquired_at(params: dict):
    seconds = params.get('actual_start_time')
    return _JEOL_EPOCH + timedelta(seconds=float(seconds)) if seconds else pd.NaT


def _nucleus_from_name(path: Path):
    stem = re.sub(r'(-\d+)+(\[\d+\])?$', '', path.stem).lower()   # drop the '-1-2' / '-1-2[1]' suffix
    for token, nucleus in sorted(_NAME_TOKENS.items(), key=lambda kv: -len(kv[0])):
        if stem.endswith('_' + token):
            return nucleus
    return None


def build_lab_nmr_inventory(root=DEFAULT_LAB_NMR_DIR) -> pd.DataFrame:
    """One row per .jdf file under root, sorted by acquisition time, with the acquisition settings that decide
    whether a spectrum is quantitative (scans, relaxation delay, pulse angle, NOE) and an MD5 of the data."""
    root = Path(root)
    rows = []
    for path in sorted(root.rglob('*.jdf')):
        jdf = read_jdf(path)
        p = jdf['params']
        rows.append({
            'folder': str(path.parent.relative_to(root)),
            'file': path.name,
            'title': jdf['title'],
            'nucleus': jdf['nucleus'],
            'nucleus_in_name': _nucleus_from_name(path),
            'data_kind': 'FID' if jdf['unit'] == 's' else 'spectrum',
            'acquired_at': _acquired_at(p),
            'duration_min': (p['end_time'] - p['actual_start_time']) / 60 if p.get('end_time') else np.nan,
            'temperature_C': p.get('temp_get'),
            'temperature_control': p.get('temp_state'),
            'experiment': p.get('experiment'),
            'scans': p.get('scans'),
            'relaxation_delay_s': p.get('relaxation_delay'),
            'pulse_angle_deg': p.get('x_angle'),
            'noe': p.get('irr_noe'),
            'lock_solvent': p.get('solvent'),
            'spectrometer_MHz': p.get('x_freq', np.nan) / 1e6,
            'data_md5': hashlib.md5(jdf['data'].tobytes()).hexdigest(),
            'path': str(path),
        })
    return pd.DataFrame(rows).sort_values('acquired_at').reset_index(drop=True)


def calculate_spectrum(fid: np.ndarray, *, sweep_Hz: float, freq_MHz: float, offset_ppm: float,
                       line_broadening_Hz: float = 1.0, clip_fraction: float = 0.1):
    """Unphased complex spectrum of a Delta FID and its ppm axis (descending).

    The digital-filter delay is removed as an integer shift to the FID maximum; the remaining fraction of a
    point is left to the first-order phase. The outer clip_fraction of the window on each side (filter
    roll-off) is dropped, as Delta does.
    """
    n = len(fid)
    delay = int(np.argmax(np.abs(fid[:min(n, 256)])))
    signal = np.roll(fid, -delay) * np.exp(-np.pi * line_broadening_Hz * np.arange(n) / sweep_Hz)
    signal[0] *= 0.5
    spectrum = np.fft.fftshift(np.fft.fft(signal))
    ppm = offset_ppm - (np.arange(n) - n // 2) * sweep_Hz / n / freq_MHz
    keep = slice(int(clip_fraction * n), int((1 - clip_fraction) * n))
    return ppm[keep], spectrum[keep]


def phase_spectrum(spectrum: np.ndarray, *, block: int = 2, n_starts: int = 6):
    """Real spectrum after automatic phasing, and (p0, p1) in rad.

    ACME criterion (Chen et al., J. Magn. Reson. 2002, 158, 164): entropy of the first derivative plus a
    penalty on negative intensity, evaluated on a block-summed copy of the spectrum for speed.
    """
    n = len(spectrum)
    u = np.linspace(-0.5, 0.5, n)
    m = n // block * block
    s_block = spectrum[:m].reshape(-1, block).sum(axis=1)
    u_block = u[:m].reshape(-1, block).mean(axis=1)

    def cost(phase):
        s = (s_block * np.exp(1j * (phase[0] + phase[1] * u_block))).real
        d = np.abs(np.diff(s))
        h = d / d.sum()
        h = h[h > 0]
        negative = s[s < 0]
        return -(h * np.log(h)).sum() + 1000.0 * (negative ** 2).sum() / (s ** 2).sum()

    fits = [minimize(cost, [p0, 0.0], method='Nelder-Mead', options={'xatol': 1e-4, 'fatol': 1e-9, 'maxiter': 3000})
            for p0 in np.linspace(-np.pi, np.pi, n_starts, endpoint=False)]
    p0, p1 = min(fits, key=lambda r: r.fun).x
    return (spectrum * np.exp(1j * (p0 + p1 * u))).real, (p0, p1)


def refine_window_phase(ppm: np.ndarray, spectrum: np.ndarray, window_ppm, *, peak_fraction: float = 0.2) -> np.ndarray:
    """Complex spectrum rotated by the zero-order phase that maximises the real intensity of the peaks inside
    window_ppm (points whose magnitude exceeds peak_fraction of the window maximum). The first-order phase
    barely changes across a narrow window, so one rotation suffices there."""
    lo, hi = sorted(window_ppm)
    inside = spectrum[(ppm >= lo) & (ppm <= hi)]
    peaks = inside[np.abs(inside) >= peak_fraction * np.abs(inside).max()]
    angle = -np.angle(peaks.sum())
    return spectrum * np.exp(1j * angle)


def load_spectrum(path, *, line_broadening_Hz: float = None, window_ppm=None):
    """ppm axis (descending) and real spectrum of one .jdf file.

    FIDs are transformed and phased automatically; spectra already processed on the spectrometer keep their
    phase. With window_ppm, the zero-order phase is refined for that window (refine_window_phase), which is the
    version to use for plots or areas of that window.
    """
    jdf = read_jdf(path)
    if jdf['unit'] == 'ppm':
        ppm = np.linspace(jdf['axis_start'], jdf['axis_stop'], len(jdf['data']))
        spectrum = jdf['data']
    else:
        p = jdf['params']
        lb = LINE_BROADENING_HZ.get(jdf['nucleus'], 1.0) if line_broadening_Hz is None else line_broadening_Hz
        ppm, spectrum = calculate_spectrum(jdf['data'], sweep_Hz=p['x_sweep'], freq_MHz=p['x_freq'] / 1e6,
                                           offset_ppm=p['x_offset'], line_broadening_Hz=lb)
        _, (p0, p1) = phase_spectrum(spectrum)
        spectrum = spectrum * np.exp(1j * (p0 + p1 * np.linspace(-0.5, 0.5, len(spectrum))))
    if window_ppm is not None:
        spectrum = refine_window_phase(ppm, spectrum, window_ppm)
    return ppm, spectrum.real


def load_text_spectrum(path):
    """ppm axis and real intensity of a Delta text export: two columns (.txt) or X/Real/Imaginary (.asc)."""
    path = Path(path)
    data = np.loadtxt(path, skiprows=1 if path.suffix == '.asc' else 0, usecols=(0, 1))
    return data[:, 0], data[:, 1]


def estimate_noise(intensity: np.ndarray, *, edge_fraction: float = 0.05) -> float:
    """Standard deviation of the first edge_fraction of the spectrum after removing a linear baseline."""
    edge = intensity[:max(10, int(edge_fraction * len(intensity)))]
    x = np.arange(len(edge))
    return float(np.std(edge - np.polyval(np.polyfit(x, edge, 1), x)))


def calculate_window_integrals(ppm: np.ndarray, intensity: np.ndarray, windows_ppm: dict, *,
                               baseline_pad_ppm: float = 0.5) -> dict:
    """Area and height of the signal in each named window {name: (lo, hi)} above a local baseline.

    The baseline is the median intensity in the two pads of baseline_pad_ppm just outside the window. Areas
    are in intensity × ppm and comparable only within one spectrum (or spectra with equal settings).
    """
    step = abs(ppm[1] - ppm[0])
    out = {}
    for name, (lo, hi) in windows_ppm.items():
        inside = (ppm >= lo) & (ppm <= hi)
        pads = ((ppm >= lo - baseline_pad_ppm) & (ppm < lo)) | ((ppm > hi) & (ppm <= hi + baseline_pad_ppm))
        y = intensity[inside] - np.median(intensity[pads])
        out[f'{name}_area'] = float(y.sum() * step)
        out[f'{name}_height'] = float(y.max())
        out[f'{name}_ppm'] = float(ppm[inside][np.argmax(y)])
    return out


def _polynomial_baseline(ppm: np.ndarray, intensity: np.ndarray, signal_free: np.ndarray, order: int) -> np.ndarray:
    return np.polyval(np.polyfit(ppm[signal_free], intensity[signal_free], order), ppm)


def calculate_area_shares(ppm: np.ndarray, intensity: np.ndarray, windows_ppm: dict, *, region_ppm,
                          baseline_orders=(1, 3), baseline_pad_ppm: float = 0.5,
                          exclusion_pad_ppm: float = 0.3) -> pd.DataFrame:
    """Share of the summed area in each named window {name: (lo, hi)}, with its uncertainty.

    The areas are computed with three baselines: the local one of calculate_window_integrals (median of two
    pads beside each window) and polynomials of each order in baseline_orders, fitted to the points of
    region_ppm outside every window (± exclusion_pad_ppm). The share is the median over the three.
    - sigma_noise: white noise of the signal-free points, propagated to the share.
    - sigma_baseline: half the spread of the share between the three baselines.
    - sigma: the two in quadrature.
    A share can come out slightly negative for an empty window; that is the size of the error, not a signal.
    """
    in_region = (ppm >= region_ppm[0]) & (ppm <= region_ppm[1])
    x, y = ppm[in_region], intensity[in_region]
    step = abs(x[1] - x[0])
    signal_free = np.ones(len(x), dtype=bool)
    for lo, hi in windows_ppm.values():
        signal_free &= ~((x >= lo - exclusion_pad_ppm) & (x <= hi + exclusion_pad_ppm))
    masks = {name: (x >= lo) & (x <= hi) for name, (lo, hi) in windows_ppm.items()}

    local = calculate_window_integrals(x, y, windows_ppm, baseline_pad_ppm=baseline_pad_ppm)
    areas = {'local': {name: local[f'{name}_area'] for name in windows_ppm}}
    for order in baseline_orders:
        corrected = y - _polynomial_baseline(x, y, signal_free, order)
        areas[f'poly{order}'] = {name: float(corrected[m].sum() * step) for name, m in masks.items()}
    noise = float(np.std(corrected[signal_free]))     # residual of the highest-order baseline

    shares = pd.DataFrame({method: pd.Series(a) / sum(a.values()) for method, a in areas.items()})
    total = sum(areas[f'poly{baseline_orders[-1]}'].values())
    sigma_area = {name: noise * step * np.sqrt(m.sum()) for name, m in masks.items()}
    out = pd.DataFrame(index=list(windows_ppm))
    out['share'] = shares.median(axis=1)
    sigma_noise = []
    for name in windows_ppm:
        s = out.loc[name, 'share']
        others = sum(sigma_area[k] ** 2 for k in windows_ppm if k != name)
        sigma_noise.append(np.sqrt((1.0 - s) ** 2 * sigma_area[name] ** 2 + s ** 2 * others) / abs(total))
    out['sigma_noise'] = sigma_noise
    out['sigma_baseline'] = 0.5 * (shares.max(axis=1) - shares.min(axis=1))
    out['sigma'] = np.hypot(out['sigma_noise'], out['sigma_baseline'])
    for method in shares.columns:
        out[f'share_{method}'] = shares[method]
    return out


def tabulate_area_shares(spectra: pd.DataFrame, windows_ppm: dict, *, region_ppm, **options) -> pd.DataFrame:
    """calculate_area_shares for every row of an inventory slice (needs 'path'), one row per (spectrum, window).

    The phase of each spectrum is refined for region_ppm. The inventory columns of each spectrum are kept.
    """
    rows = []
    for _, spectrum in spectra.iterrows():
        ppm, intensity = load_spectrum(spectrum['path'], window_ppm=region_ppm)
        shares = calculate_area_shares(ppm, intensity, windows_ppm, region_ppm=region_ppm, **options)
        for window, values in shares.iterrows():
            rows.append({**spectrum.to_dict(), 'window': window, **values.to_dict()})
    return pd.DataFrame(rows)


def find_heated_windows(spectra: pd.DataFrame, *, rt_max_C: float = 25.0) -> pd.DataFrame:
    """Clock intervals in which one sample was recorded above rt_max_C (inventory rows of that sample).

    Consecutive spectra at the same temperature form one interval, from the start of the first to the end of
    the last. This is the shortest time the sample can have spent at that temperature: it may have been held
    there longer before or after the spectra, which the files do not record.
    """
    rows = spectra.sort_values('acquired_at')
    out = []
    for _, r in rows.iterrows():
        end = r['acquired_at'] + pd.to_timedelta(r['duration_min'], unit='min')
        heated = r['temperature_C'] > rt_max_C
        if heated and out and out[-1]['open'] and abs(out[-1]['T_C'] - r['temperature_C']) < 0.5:
            out[-1].update(end=end, n_spectra=out[-1]['n_spectra'] + 1)
        elif heated:
            out.append({'start': r['acquired_at'], 'end': end, 'T_C': round(float(r['temperature_C']), 1),
                        'n_spectra': 1, 'open': True})
        elif out:
            out[-1]['open'] = False
    table = pd.DataFrame(out, columns=['start', 'end', 'T_C', 'n_spectra', 'open']).drop(columns='open')
    table['minutes'] = (table['end'] - table['start']).dt.total_seconds() / 60.0
    return table
