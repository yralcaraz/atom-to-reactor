"""Curated lab observables: everything a fit reads from the lab NMR spectra, in one file.

The file (`lab_observables.json`, in the data folder, see `kinetics/paths.py`) holds, per sample: the recipe, what is
known of its temperature history, and the area shares of every quantitative spectrum with their two measured
uncertainties (noise, baseline spread). It is built from two share tables in the format of
`kinetics.data.lab_nmr.tabulate_area_shares` (plus the columns 'sample' and 'hours_since_first'), so the same
builder runs from the raw spectra or from the tables exported by notebook 03.

    load_share_tables          the two share tables exported by notebook 03 (CSV)
    build_lab_inventory_for_observables  the same two tables, and the heating episodes, from the raw spectra
    build_lab_observables      share tables + recorded histories → the observables dict
    write_lab_observables      save it as JSON
    load_lab_observables       read it back
    describe_lab_observables   one row per sample: spectra, temperatures, times, shares
    calculate_replicate_scatter  scatter of repeated spectra of one sample against their stated noise

What the files do not record is kept explicit: the time from mixing to a sample's first spectrum ('age') is
unknown for every sample, and the heating episodes carry a 'basis' saying where each number comes from.

Source: Y. Alcaraz Galván; sample histories transcribed from the acquisition times of the experimental NMR
spectra, as recorded in notebooks 03 and 04
"""

import json
import os
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from kinetics.paths import DATA_DIR_VAR, data_path

DEFAULT_OBSERVABLES_PATH = data_path('data', 'lab_observables.json')
DEFAULT_SHARES_DIR = data_path('results', '03')
DEFAULT_SAMPLE_FOLDERS_PATH = data_path('data', 'lab_sample_folders.json')

PHOSPHATE_WINDOWS = ('TMSPA', 'BMSPA', 'MMSPA', 'H3PO4')        # ³¹P: one P per molecule
SILYL_WINDOWS = ('TMSOH', 'HMDSO', 'TMSOEG')                    # ¹³C, Si–CH3 region: three methyl C per Si
RT_MAX_C = 25.0

# Error floor on a share, per nucleus: the larger of an absolute value and a fraction of the share. It stands
# for what the two measured terms cannot see (relaxation, NOE differences between peaks). DECLARED, not measured.
DEFAULT_ERROR_FLOOR = {'31P': {'absolute': 0.01, 'relative': 0.05},
                       '13C': {'absolute': 0.03, 'relative': 0.15}}

# Samples, named as in notebook 03 (§2.1). 'recipe' is resolved by _recipe_molarities.
# role: 'fit' enters the fit; 'hold_out' is predicted, never fitted; 'check' is compared after the fit.
LAB_SAMPLES = {
    '0.5 % H2O': {
        'label': '5 vol% TMSPa + 0.5 vol% H2O in EC/DEC', 'nucleus': '31P', 'role': 'fit',
        'recipe': ('water_series', 0.5),
        'flags': ['lock solvent D2O; every other sample is locked on DMSO-d6 (notebook 04 §5)',
                  '31P lines near -140.1, -144.5 and -148.9 ppm, 4.4 ppm apart, the pattern of PF6- '
                  '(notebook 04, Appendix A): if the sample holds LiPF6 it is outside the network. AUDIT OPEN'],
    },
    '2 % H2O': {
        'label': '5 vol% TMSPa + 2 vol% H2O in EC/DEC', 'nucleus': '31P', 'role': 'fit',
        'recipe': ('water_series', 2.0), 'flags': [],
    },
    'TMSPa + TMSOH (A)': {
        'label': '5 vol% TMSPa + 2 vol% TMSOH in EC/DEC, tube A', 'nucleus': '31P', 'role': 'fit',
        'recipe': ('tmspa_tmsoh',),
        'flags': ['later spectra of this tube sit in a folder whose files are titled 5 v% TMSOH (notebook 04 §5): '
                  'not used'],
    },
    'TMSPa + TMSOH (B)': {
        'label': '5 vol% TMSPa + 2 vol% TMSOH in EC/DEC, tube B', 'nucleus': '31P', 'role': 'hold_out',
        'recipe': ('tmspa_tmsoh',), 'flags': ['16 scans: standard errors of 0.08 to 0.18 on the shares'],
    },
    'TMSPa alone': {
        'label': '5 vol% TMSPa in EC/DEC, no water added', 'nucleus': '31P', 'role': 'check',
        'recipe': ('tmspa',),
        'flags': ['TMSPA converts to BMSPA in 10 days without added water: the water content is unknown '
                  '(notebook 03 §5.4)'],
    },
    'TMSOH, probe': {
        'label': '5 vol% TMSOH in EC/DEC, heated in the NMR probe', 'nucleus': '13C', 'role': 'fit',
        'recipe': ('tmsoh',), 'flags': [],
    },
    'TMSOH, glovebox': {
        'label': '5 vol% TMSOH in EC/DEC, stirred at 80 °C in a glovebox', 'nucleus': '13C', 'role': 'fit',
        'recipe': ('tmsoh',),
        'flags': ['a 19F spectrum of this sample exists and has not been inspected. AUDIT OPEN',
                  'an unassigned 13C line at -1.75 ppm falls in the TMSOEG window and is counted with it '
                  '(notebook 03 §2.1)'],
    },
}

# What notebooks 03 and 04 record of each sample's history besides its room-temperature spectra.
# 'heated': episodes after the sample's first spectrum [start_h and duration_h on the file clock].
# 'pre_history': what happened between mixing and the first spectrum, where something is known.
RECORDED_HISTORIES = {
    '2 % H2O': {'heated': [
        {'T_C': 30.0, 'start_h': 16.0, 'duration_h': 0.30,
         'basis': 'notebook 04 §5: the 30 °C spectrum follows the first room-temperature one by 16 h'},
        {'T_C': 40.0, 'start_h': 24.0, 'duration_h': 0.30,
         'basis': 'date from notebook 04 §4.1 (23 June, after 30 °C, before the second room-temperature spectrum); '
                  'clock time ASSUMED'},
        {'T_C': 50.0, 'start_h': 64.0, 'duration_h': 0.30,
         'basis': 'date from notebook 04 §4.1 (25 June); clock time ASSUMED from the 30 °C step'},
        {'T_C': 60.0, 'start_h': 88.0, 'duration_h': 0.30,
         'basis': 'date from notebook 04 §4.1 (26 June); clock time ASSUMED from the 30 °C step'},
        {'T_C': 70.0, 'start_h': 112.0, 'duration_h': 0.30,
         'basis': 'date from notebook 04 §4.1 (27 June); clock time ASSUMED from the 30 °C step'},
        {'T_C': 80.0, 'start_h': 136.0, 'duration_h': 0.30,
         'basis': 'date from notebook 04 §4.1 (28 June); clock time ASSUMED from the 30 °C step'},
    ], 'heated_note': 'durations: notebook 03 §5.2, 1.8 h in all shared equally (about 17 min per step). They are '
                      'the time covered by the spectra at each temperature, so a MINIMUM; the hold may have been '
                      'longer (the paper protocol has 8 h holds)'},
    'TMSOH, probe': {'heated': [
        {'T_C': 80.0, 'start_h': 138.66, 'duration_h': 0.145,
         'basis': 'notebook 03 §2.2: the spectra at 80 °C span 8.7 min (MINIMUM heating; at most the 5.6 days '
                  'between the room-temperature spectra); placed around the 80 °C 13C spectrum'},
    ]},
    'TMSOH, glovebox': {'pre_history': [
        {'T_C': 80.0, 'duration_h': 8.0,
         'basis': 'Gogoi 2024: 8 h hold at 80 °C; ASSUMED to be this sample (notebook 03 §2.2)'},
    ]},
}

OPEN_AUDIT = [
    'Raw spectra: check the 31P region near -144 ppm (PF6-) in every TMSPa sample, not only the 0.5 % one.',
    'Raw spectra: inspect the 19F spectrum of the glovebox TMSOH sample (PF6- near -72 ppm, fluorosilanes).',
    'Raw spectra: replace the ASSUMED clock times and the minimum durations of the heating steps by '
    'find_heated_windows on the inventory.',
]


def load_share_tables(directory: str = DEFAULT_SHARES_DIR) -> tuple:
    """(³¹P table, ¹³C table) as exported by notebook 03 (`lab_shares_31P.csv`, `lab_shares_13C.csv`, in the data folder)."""
    tables = []
    for name in ('lab_shares_31P.csv', 'lab_shares_13C.csv'):
        path = os.path.join(directory, name)
        if not os.path.exists(path):
            raise FileNotFoundError(f'{path} not found: the share tables are read from the data folder (set {DATA_DIR_VAR}, '
                                    f'see .env.example)')
        tables.append(pd.read_csv(path, index_col=0, parse_dates=['acquired_at']))
    return tuple(tables)


def load_lab_sample_folders(path: str = DEFAULT_SAMPLE_FOLDERS_PATH) -> dict:
    """Names of the raw lab folders (and file titles) that hold each sample, read from the data folder.

    {'raw_sample_folders': {sample: [folder, ...]}  (notebook 03 §2.1),
     'sample_of_folder': {folder: sample}, 'sample_of_title': {title: sample}  (notebook 04 §1)}.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f'{path} not found: the folder names are read from the data folder (set {DATA_DIR_VAR}, '
                                f'see .env.example)')
    with open(path, 'r', encoding='utf-8') as f:
        return json.load(f)


# Integration windows used for the area shares (notebook 03 §2.1)
P_WINDOWS_PPM = {'TMSPA': (-27.0, -22.0), 'BMSPA': (-16.5, -12.5), 'MMSPA': (-8.5, -5.0), 'H3PO4': (-1.0, 3.5)}
C_WINDOWS_PPM = {'TMSOH': (0.05, 0.45), 'HMDSO': (0.6, 1.1), 'TMSOEG': (-1.95, -1.45)}


def build_lab_inventory_for_observables(root) -> tuple:
    """(³¹P share table, ¹³C share table, {sample: heated windows}) from the raw spectra under `root`.

    The selection is that of notebook 03 §2.1: ³¹P at room temperature with the NOE off and a relaxation delay
    of at least 5 s; ¹³C (Si–CH3 region) of the two TMSOH samples. Needs the raw spectra (LAB_NMR_DIR).
    """
    from kinetics.data.lab_nmr import build_lab_nmr_inventory, find_heated_windows, tabulate_area_shares
    inventory = build_lab_nmr_inventory(root).drop_duplicates('data_md5')
    folders = load_lab_sample_folders()['raw_sample_folders']
    inventory['sample'] = inventory['folder'].map({f: s for s, names in folders.items() for f in names})
    inventory = inventory[inventory['sample'].notna()].copy()
    first = inventory.groupby('sample')['acquired_at'].transform('min')
    inventory['hours_since_first'] = (inventory['acquired_at'] - first).dt.total_seconds() / 3600.0
    p31 = inventory[(inventory['nucleus'] == '31P') & inventory['sample'].str.contains('TMSPa|H2O')
                    & (inventory['temperature_C'] < RT_MAX_C) & (inventory['noe'] == 'FALSE')
                    & (inventory['relaxation_delay_s'] >= 5)]
    c13 = inventory[(inventory['nucleus'] == '13C') & inventory['sample'].str.startswith('TMSOH')]
    shares_P = tabulate_area_shares(p31, P_WINDOWS_PPM, region_ppm=(-32.0, 6.0))
    shares_C = tabulate_area_shares(c13, C_WINDOWS_PPM, region_ppm=(-3.0, 2.0))
    heated = {}
    for sample in ('2 % H2O', 'TMSOH, probe'):
        table = find_heated_windows(inventory[inventory['sample'] == sample])
        if len(table):
            heated[sample] = table
    return shares_P, shares_C, heated


def _recipe_molarities(recipe: tuple) -> tuple:
    """(c0_M of the solutes and EC, basis) for a recipe key; concentrations as in notebook 03 §2.2."""
    from kinetics.reactor.validation import calculate_water_series_c0   # reactor imports data: keep this local
    tmspa_M = calculate_water_series_c0(0.5)['TMSPA']
    kind = recipe[0]
    if kind == 'water_series':
        c0 = {sp: c for sp, c in calculate_water_series_c0(recipe[1]).items() if c > 0}
        return c0, f'Gogoi 2024 water series: {recipe[1]:g} vol% H2O stock diluted by 5 vol% TMSPa; EC 7.1 M'
    if kind == 'tmspa_tmsoh':
        return {'TMSPA': tmspa_M, 'TMSOH': 0.18, 'EC': 7.1}, 'recipe of control E4: 5 vol% TMSPa, 2 vol% TMSOH = 0.18 M'
    if kind == 'tmsoh':
        return {'TMSOH': 0.45, 'EC': 7.1}, 'recipe of controls E1–E3: 5 vol% TMSOH = 0.45 M'
    if kind == 'tmspa':
        return {'TMSPA': tmspa_M, 'EC': 7.1}, '5 vol% TMSPa; water not added and not measured'
    raise KeyError(f'Unknown recipe {recipe!r}')


def _spectrum_record(rows: pd.DataFrame, windows: tuple) -> dict:
    first = rows.iloc[0]
    by_window = rows.set_index('window')
    missing = [w for w in windows if w not in by_window.index]
    if missing:
        raise KeyError(f"spectrum {first['file']}: no share for window(s) {missing}")
    return {
        't_since_first_h': float(first['hours_since_first']),
        'acquired_at': pd.Timestamp(first['acquired_at']).isoformat(),
        'temperature_C': float(first['temperature_C']),
        'scans': int(first['scans']),
        'relaxation_delay_s': float(first['relaxation_delay_s']),
        'pulse_angle_deg': float(first['pulse_angle_deg']),
        'noe': str(first['noe']).upper(),
        'lock_solvent': str(first['lock_solvent']),
        'data_md5': str(first['data_md5']),
        'shares': {w: {'share': float(by_window.loc[w, 'share']),
                       'sigma_noise': float(by_window.loc[w, 'sigma_noise']),
                       'sigma_baseline': float(by_window.loc[w, 'sigma_baseline'])} for w in windows},
    }


def build_lab_observables(shares_P: pd.DataFrame, shares_C: pd.DataFrame, *, built_from: str,
                          heated_windows: dict = None, samples: dict = None,
                          error_floor: dict = None) -> dict:
    """The observables dict from the two share tables and the recorded histories.

    shares_P, shares_C: one row per (spectrum, window) with the columns of `tabulate_area_shares` plus 'sample'
    (a key of `samples`) and 'hours_since_first' (time since the first spectrum of that sample).
    built_from: what the tables were computed from; stored in the file.
    heated_windows: {sample: table of `find_heated_windows`} when the raw inventory is at hand; the episodes
    of RECORDED_HISTORIES are used for every sample without an entry.
    """
    samples = samples if samples is not None else LAB_SAMPLES
    out = {}
    for name, spec in samples.items():
        table, windows = (shares_P, PHOSPHATE_WINDOWS) if spec['nucleus'] == '31P' else (shares_C, SILYL_WINDOWS)
        rows = table[table['sample'] == name]
        if rows.empty:
            raise KeyError(f"no spectrum of sample '{name}' in the {spec['nucleus']} table")
        spectra = [_spectrum_record(grp, windows) for _, grp in rows.groupby('acquired_at', sort=True)]
        first_at = pd.Timestamp(spectra[0]['acquired_at']) - pd.Timedelta(hours=spectra[0]['t_since_first_h'])
        recorded = RECORDED_HISTORIES.get(name, {})
        if heated_windows is not None and name in heated_windows:
            heated = [{'T_C': float(r['T_C']), 'start_h': (pd.Timestamp(r['start']) - first_at).total_seconds() / 3600.0,
                       'duration_h': float(r['minutes']) / 60.0,
                       'basis': 'find_heated_windows on the raw inventory (MINIMUM duration: time covered by the spectra)'}
                      for _, r in heated_windows[name].iterrows()]
        else:
            heated = [dict(h) for h in recorded.get('heated', [])]
        rt = [s['temperature_C'] for s in spectra if s['temperature_C'] < RT_MAX_C]
        c0_M, recipe_basis = _recipe_molarities(spec['recipe'])
        pre_history = [dict(h) for h in recorded.get('pre_history', [])]
        out[name] = {
            'label': spec['label'],
            'campaign': int(first_at.year),
            'nucleus': spec['nucleus'],
            'role': spec['role'],
            'c0_M': c0_M,
            'recipe_basis': recipe_basis,
            'first_spectrum_at': first_at.isoformat(),
            'T_rt_C': float(np.mean(rt)),
            'age': {'known': bool(pre_history),
                    'note': 'time from mixing to the first spectrum is not recorded' if not pre_history
                    else 'history before the first spectrum taken from the paper'},
            'pre_history': pre_history,
            'heated': heated,
            'heated_note': recorded.get('heated_note', ''),
            'spectra': spectra,
            'flags': list(spec['flags']),
        }
        if spec['nucleus'] == '13C':
            out[name]['flags'].append('13C recorded with NOE on and a 2 s relaxation delay: shares assume the same '
                                      'NOE and relaxation for the Si-CH3 carbons of every species')
    return {
        '_meta': {
            'source': 'Y. Alcaraz Galván; NMR data from experiments',
            'built_from': built_from,
            'built_at_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
            'shares': 'share of the summed area of the windows of one spectrum: median over three baselines; '
                      'sigma_noise from the white noise, sigma_baseline half the spread between the baselines '
                      '(kinetics.data.lab_nmr.calculate_area_shares)',
            'windows': {'31P': list(PHOSPHATE_WINDOWS), '13C': list(SILYL_WINDOWS)},
            'error_floor': error_floor if error_floor is not None else DEFAULT_ERROR_FLOOR,
            'open_audit': list(OPEN_AUDIT),
        },
        'samples': out,
    }


def write_lab_observables(observables: dict, path: str = DEFAULT_OBSERVABLES_PATH) -> str:
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(observables, f, indent=1, ensure_ascii=False)
    return path


def load_lab_observables(path: str = DEFAULT_OBSERVABLES_PATH) -> dict:
    """The observables dict from `lab_observables.json` in the data folder (a fresh copy on every call)."""
    if not os.path.exists(path):
        raise FileNotFoundError(f'{path} not found: build it with scripts/build_observables.py '
                                f'(set {DATA_DIR_VAR}, see .env.example)')
    with open(path, 'r', encoding='utf-8') as f:
        return json.load(f)


def tabulate_observable_shares(observables: dict) -> pd.DataFrame:
    """One row per (sample, spectrum, window): share, sigma_noise, sigma_baseline and the spectrum's time and T."""
    rows = []
    for name, sample in observables['samples'].items():
        for k, spectrum in enumerate(sample['spectra']):
            for window, v in spectrum['shares'].items():
                rows.append({'sample': name, 'spectrum': k, 'acquired_at': spectrum['acquired_at'],
                             't_since_first_h': spectrum['t_since_first_h'], 'temperature_C': spectrum['temperature_C'],
                             'window': window, **v})
    return pd.DataFrame(rows)


def describe_lab_observables(observables: dict) -> pd.DataFrame:
    """One row per sample: role, spectra, temperatures, times on the file clock, what is known of its history."""
    rows = {}
    for name, s in observables['samples'].items():
        t = [sp['t_since_first_h'] for sp in s['spectra']]
        T = sorted({round(sp['temperature_C']) for sp in s['spectra']})
        heated = ', '.join(f"{h['T_C']:g} °C × {60 * h['duration_h']:.0f} min" for h in s['heated'])
        pre = ', '.join(f"{h['T_C']:g} °C × {h['duration_h']:g} h" for h in s['pre_history'])
        rows[name] = {
            'role': s['role'], 'campaign': s['campaign'], 'spectra': f"{len(s['spectra'])} {s['nucleus']}",
            'solutes at t = 0': ', '.join(f'{sp} {1000 * c:.0f} mM' for sp, c in s['c0_M'].items() if sp != 'EC'),
            'T of the spectra (°C)': ', '.join(str(v) for v in T),
            'time since first spectrum (h)': f'{min(t):.2f}' if len(t) == 1 else f'{min(t):.2f} to {max(t):.1f}',
            'before the first spectrum': pre or 'unknown (age)',
            'heated after the first spectrum': heated or '—',
            'flags': len(s['flags']),
        }
    return pd.DataFrame.from_dict(rows, orient='index')


def calculate_replicate_scatter(observables: dict, sample: str) -> dict:
    """Scatter of the repeated room-temperature spectra of one sample against their stated white noise.

    Per window: χ² of the shares about their noise-weighted mean (n − 1 degrees of freedom). Pooled: the sum
    over windows on (K − 1)(n − 1) degrees of freedom, because the K shares of a spectrum sum to one. The Birge
    ratio sqrt(χ²/dof) says by how much the spectrum-to-spectrum scatter exceeds the noise; a ratio near 1 means
    the baseline error is common to the spectra.
    """
    spectra = [sp for sp in observables['samples'][sample]['spectra'] if sp['temperature_C'] < RT_MAX_C]
    if len(spectra) < 2:
        raise ValueError(f"'{sample}' has no repeated room-temperature spectra")
    windows = list(spectra[0]['shares'])
    rows, chi2_sum = {}, 0.0
    for w in windows:
        y = np.array([sp['shares'][w]['share'] for sp in spectra])
        s = np.array([sp['shares'][w]['sigma_noise'] for sp in spectra])
        mean = np.sum(y / s ** 2) / np.sum(1.0 / s ** 2)
        chi2 = float(np.sum(((y - mean) / s) ** 2))
        chi2_sum += chi2
        rows[w] = {'mean': float(mean), 'std': float(np.std(y, ddof=1)), 'median sigma_noise': float(np.median(s)),
                   'chi2': chi2, 'dof': len(y) - 1, 'birge': float(np.sqrt(chi2 / (len(y) - 1)))}
    dof = (len(windows) - 1) * (len(spectra) - 1)
    return {'per_window': pd.DataFrame.from_dict(rows, orient='index'), 'chi2': chi2_sum, 'dof': dof,
            'birge': float(np.sqrt(chi2_sum / dof)), 'n_spectra': len(spectra)}
