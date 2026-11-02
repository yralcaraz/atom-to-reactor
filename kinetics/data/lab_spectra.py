"""Every spectrum of every tube: its class (used in the fit, usable, noise) and its shares.

The fits read 22 of the 89 .jdf files of the raw folders (kinetics.data.observables). Notebook 05 rev08 §2 sorted
every file of every tube into three classes and read the usable ones as the fitted ones are read; this module is
that sorting and reading, so that a fit run can use it (the thermodynamic fit, kinetics.fitting.thermo_fit). The
rules are those of rev08 §2, unchanged:

- 'used in the fit': a spectrum of a block of the observables with the role 'fit';
- 'usable, not fitted': a block of another role (check, hold-out), or a line at least LINE_SNR times the noise of
  the region where lines are looked for (LINE_REGIONS);
- 'noise': the rest.

A file whose data are those of another one (same md5) is a copy: it takes the class of the original.

Shares: ³¹P and ¹³C as in the observables (area of the windows, moved by the axis offset of the tube); ²⁹Si by height
(its lines are one or two points wide on a rolling baseline); ¹H in two windows of the Si–CH3 region (silyl still on
phosphate, free silyl). A spectrum of a block of the observables keeps the shares stored there.

    list_tube_folders        the raw folders of each tube (lab_sample_folders.json, check folders included)
    calculate_peak_to_noise  tallest point of the line region over its noise
    classify_tube_files      one row per .jdf of a tube: header, class, time on the clock of the tube
    get_share_reading        windows, region and options of the shares of one nucleus in one tube
    calculate_height_shares  ²⁹Si: shares of the summed height of the windows
    read_spectrum_shares     shares of one spectrum, with their noise and baseline errors
    tabulate_tube_shares     shares of every usable spectrum of a tube, one row per (spectrum, window)

Source: Y. Alcaraz Galván; reads NMR data from experiments; rules of notebook 05 rev08 §2 (2026-10-28)
"""

import functools

import numpy as np
import pandas as pd

from kinetics.data.lab_nmr import calculate_area_shares, load_spectrum
from kinetics.data.observables import (
    C_PHOSPHATE_WINDOWS_PPM, C_REGION_PPM, C_WINDOWS_PPM, P_WINDOWS_PPM, list_spectrum_sets, load_lab_sample_folders,
)

FILE_CLASSES = ('used in the fit', 'usable, not fitted', 'noise')
P_REGION_PPM = (-32.0, 6.0)              # region of the ³¹P shares (build_lab_inventory_for_observables)
LINE_REGIONS = {'31P': (-32.0, 6.0), '13C': (-3.0, 2.5), '1H': (-0.35, 0.10), '29Si': (0.0, 35.0), '19F': (-300.0, 100.0)}   # ppm: where a line is looked for
LINE_SNR = 10.0       # a line: its tallest point at least 10 × the noise of the region. Found on 28 Oct: lines 14 to 600, noise 3 to 7
# ²⁹Si, ppm as recorded (not referenced): lines at 7.8 (HMDSO), 14.2 (TMSOH), 18.0 (TMS-EG) and 21.6–23.0 (silyl on
# phosphate). One Si per silyl group, so a share is a share of the silyl groups, as in ¹³C
SI_WINDOWS = {'TMSOH': (13.0, 16.0), 'HMDSO': (6.5, 9.0), 'TMSOEG': (17.0, 20.0), 'P-silyl': (20.5, 25.0)}
SI_REGION = (-5.0, 35.0)
# ¹H, Si–CH3 region: one line per form of the silyl group, 9 H each, so an area share is a share of the silyl groups.
# TMSOH and HMDSO are only 0.004–0.011 ppm apart, so one window takes all the free silyl
H1_WINDOWS = {'P-bound silyl': (-0.100, 0.060), 'free silyl': (-0.300, -0.150)}
H1_REGION = (-1.0, 0.6)
H1_OPTIONS = {'baseline_pad_ppm': 0.1, 'exclusion_pad_ppm': 0.05}
# Error floor of a share that the observables do not hold: ²⁹Si takes that of ¹³C (NOE, like it), ¹H that of ³¹P (no NOE). DECLARED
FLOOR_OF = {'31P': '31P', '13C': '13C', '29Si': '13C', '1H': '31P'}



@functools.lru_cache(maxsize=1)
def list_tube_folders() -> dict:
    """{tube: [raw folder, ...]}: lab_sample_folders.json of the data folder, its check folders included."""
    listed = load_lab_sample_folders()
    return {tube: folders + listed.get('raw_check_folders', {}).get(tube, [])
            for tube, folders in listed['raw_sample_folders'].items()}


@functools.lru_cache(maxsize=None)
def _load_phased(path: str, window=None):
    # one file is read for its class and for its shares
    return load_spectrum(path, window_ppm=window)


def get_axis_offset(tube: str, nucleus: str, observables: dict) -> float:
    """The shift of the axis of a tube for one nucleus (LAB_SAMPLES: ¹³C of the tubes with phosphate); 0 elsewhere."""
    return observables['samples'][tube].get('more_nuclei', {}).get(nucleus, {}).get('axis_offset_ppm', 0.0)


def calculate_peak_to_noise(path: str, nucleus: str, offset: float = 0.0) -> float:
    """Tallest point of the region where a line is looked for, over the noise of that region: the spread (MAD) of the
    differences between neighbouring points, which a sloping baseline or a broad hump (the glass, near −110 ppm) does not raise."""
    low, high = LINE_REGIONS[nucleus]
    ppm, intensity = _load_phased(path, (low + offset, high + offset))
    y = intensity[(ppm >= low + offset) & (ppm <= high + offset)]
    step = np.diff(y)
    return float((y.max() - np.median(y)) / (1.4826 * np.median(np.abs(step - np.median(step))) / np.sqrt(2)))


def classify_tube_files(tube: str, observables: dict, inventory: pd.DataFrame) -> pd.DataFrame:
    """One row per .jdf of the folders of a tube, in the order acquired: the inventory columns, 'role' (that of its block
    in the observables, NaN for a file the observables do not hold), 'copy of', 'peak/noise', 'class' and 't_h' (hours
    since the first file of the tube, the clock of the observables).

    inventory: build_lab_nmr_inventory of the raw folder, copies included.
    """
    roles = {spectrum['data_md5']: role for _, role, spectra in list_spectrum_sets(observables['samples'][tube])
             for spectrum in spectra}
    files = inventory[inventory['folder'].isin(list_tube_folders()[tube])].sort_values(['acquired_at', 'file']).reset_index(drop=True)
    files['role'] = files['data_md5'].map(roles)
    original = files.drop_duplicates('data_md5').set_index('data_md5')['file']
    files['copy of'] = [original[md5] if copy else None for md5, copy in zip(files['data_md5'], files.duplicated('data_md5'))]
    files['peak/noise'] = [calculate_peak_to_noise(r['path'], r['nucleus'], get_axis_offset(tube, r['nucleus'], observables))
                           for _, r in files.iterrows()]
    files['class'] = [FILE_CLASSES[0] if role == 'fit' else FILE_CLASSES[1] if isinstance(role, str) or snr >= LINE_SNR
                      else FILE_CLASSES[2] for role, snr in zip(files['role'], files['peak/noise'])]
    files['t_h'] = (files['acquired_at'] - files['acquired_at'].min()).dt.total_seconds() / 3600.0
    return files


def get_share_reading(nucleus: str, tube: str, observables: dict):
    """(windows, region, options) of the shares of one nucleus in one tube: those of the fitted spectra for ³¹P and ¹³C
    (moved by the axis offset of the tube), SI_WINDOWS for ²⁹Si, H1_WINDOWS for ¹H; None for a nucleus that no species
    of the model holds. A tube without phosphate has no window for the silyl on phosphate."""
    offset = get_axis_offset(tube, nucleus, observables)
    phosphate = 'TMSPA' in observables['samples'][tube]['c0_M']
    if nucleus == '31P':
        return P_WINDOWS_PPM, P_REGION_PPM, {}
    if nucleus == '13C':
        windows = C_PHOSPHATE_WINDOWS_PPM if phosphate else C_WINDOWS_PPM
        return ({window: (low + offset, high + offset) for window, (low, high) in windows.items()},
                (C_REGION_PPM[0] + offset, C_REGION_PPM[1] + offset), {})
    if nucleus == '29Si':
        return {window: limits for window, limits in SI_WINDOWS.items() if window != 'P-silyl' or phosphate}, SI_REGION, {}
    if nucleus == '1H':
        return ({window: limits for window, limits in H1_WINDOWS.items() if window != 'P-bound silyl' or phosphate},
                H1_REGION, H1_OPTIONS)
    return None


def calculate_height_shares(ppm: np.ndarray, intensity: np.ndarray, windows_ppm: dict, *, region_ppm, pad_ppm: float = 0.5) -> pd.DataFrame:
    """²⁹Si: share of the summed height of the windows. Its lines are one or two points wide on a rolling baseline, so the
    area of a window is mostly baseline; every line is an Si(CH3)3 group of about the same width. Each height is the
    tallest point of the window over the median of two pads beside it; an empty window reads the top of its noise,
    about 4 σ. sigma: the noise (MAD of neighbouring differences) propagated to the share."""
    inside = (ppm >= region_ppm[0]) & (ppm <= region_ppm[1])
    x, y = ppm[inside], intensity[inside]
    step = np.diff(y)
    noise = 1.4826 * np.median(np.abs(step - np.median(step))) / np.sqrt(2)
    heights = {}
    for name, (low, high) in windows_ppm.items():
        pads = ((x >= low - pad_ppm) & (x < low)) | ((x > high) & (x <= high + pad_ppm))
        heights[name] = float(y[(x >= low) & (x <= high)].max() - np.median(y[pads]))
    total = sum(heights.values())
    shares = pd.DataFrame({'share': {name: h / total for name, h in heights.items()}, 'height': heights})
    shares['sigma'] = [noise * np.sqrt((1.0 - s) ** 2 + s ** 2 * (len(heights) - 1)) / total for s in shares['share']]
    return shares


def read_spectrum_shares(path: str, nucleus: str, tube: str, observables: dict, *, data_kind: str = 'FID') -> pd.DataFrame:
    """Shares of one spectrum ('share', 'sigma_noise', 'sigma_baseline'; one row per window), read as rev08 §2 reads them.
    The phase of a ¹H FID is refined on the Si–CH3 region; a ¹H spectrum processed on the spectrometer keeps its phase."""
    windows, region, options = get_share_reading(nucleus, tube, observables)
    window = (LINE_REGIONS['1H'] if data_kind == 'FID' else None) if nucleus == '1H' else region
    ppm, intensity = _load_phased(path, window)
    if nucleus == '29Si':
        shares = calculate_height_shares(ppm, intensity, windows, region_ppm=region)
        return pd.DataFrame({'share': shares['share'], 'sigma_noise': shares['sigma'], 'sigma_baseline': 0.0})
    return calculate_area_shares(ppm, intensity, windows, region_ppm=region, **options)[['share', 'sigma_noise', 'sigma_baseline']]


def tabulate_tube_shares(tube: str, observables: dict, files: pd.DataFrame) -> pd.DataFrame:
    """Shares of every usable spectrum of a tube (classify_tube_files), copies once, one row per (spectrum, window): 'tube',
    'nucleus', 'file', 'path', 'data_md5', 'acquired_at', 't_h', 'T_C', 'role', 'scans', 'noe', 'relaxation_delay_s',
    'experiment', 'data_kind', 'window', 'share', 'sigma_noise', 'sigma_baseline', 'source' ('observables' or 'read here').
    A spectrum of a block of the observables keeps its stored shares; a nucleus no species holds (¹⁹F) is left out."""
    stored = {spectrum['data_md5']: spectrum['shares'] for _, _, spectra in list_spectrum_sets(observables['samples'][tube])
              for spectrum in spectra}
    usable = files[(files['class'] != FILE_CLASSES[2]) & files['copy of'].isna()]
    rows = []
    for _, r in usable.iterrows():
        if get_share_reading(r['nucleus'], tube, observables) is None:
            continue
        if r['data_md5'] in stored:
            shares, source = pd.DataFrame(stored[r['data_md5']]).T, 'observables'
        else:
            shares, source = read_spectrum_shares(r['path'], r['nucleus'], tube, observables, data_kind=r['data_kind']), 'read here'
        head = {'tube': tube, 'nucleus': r['nucleus'], 'file': r['file'], 'path': r['path'], 'data_md5': r['data_md5'],
                'acquired_at': r['acquired_at'], 't_h': r['t_h'], 'T_C': r['temperature_C'], 'role': r['role'],
                'scans': r['scans'], 'noe': r['noe'], 'relaxation_delay_s': r['relaxation_delay_s'], 'experiment': r['experiment'],
                'data_kind': r['data_kind']}
        rows += [{**head, 'window': window, 'share': float(s['share']), 'sigma_noise': float(s['sigma_noise']),
                  'sigma_baseline': float(s['sigma_baseline']), 'source': source} for window, s in shares.iterrows()]
    return pd.DataFrame(rows)
