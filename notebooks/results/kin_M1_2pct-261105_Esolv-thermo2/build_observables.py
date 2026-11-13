"""The observables of this run: ONE tube, the 2 % H2O of 5 Nov 2026 (our spectra, fit_thermo_rev01 §2). No other tube.

    python notebooks/results/kin_M1_2pct-261105_Esolv-thermo2/build_observables.py

- ³¹P only: the five quantitative spectra (NOE off, d1 5 s), read as every ³¹P of the fit (tabulate_area_shares,
  P_WINDOWS_PPM, region −32…6 ppm). 31P-2-1 holds the data of 31P-1-1 (same md5): read once.
- ¹H not used: its Si–CH3 region is noise (S/N 3–7, §2.2.2); no ¹³C was recorded.
- Times known: water added at 16:33:04 (6 min 50 s before the end of the first ³¹P, as in §2) = t 0. The time to the start
  of the first ³¹P is a known room-temperature history ('pre_history'). Spectrum times: the start of each acquisition.
- Composition at mixing: the recipe 5 vol% TMSPa + 2 vol% H2O in EC/DEC (calculate_water_series_c0(2.0), EC 7.1 M).
  Room temperature: the mean of the spectra.
- 'replicates': False — the tube reacts while recorded: each spectrum has its own error (noise ⊕ baseline ⊕ floor).
- '_meta': that of the data folder's observables (error floor), the samples dropped.
"""

import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))

from kinetics.data.lab_nmr import build_lab_nmr_inventory, tabulate_area_shares   # noqa: E402
from kinetics.data.observables import P_WINDOWS_PPM, WINDOWS_OF_NUCLEUS, _spectrum_record   # noqa: E402
from kinetics.reactor.validation import calculate_water_series_c0   # noqa: E402

META_SOURCE = Path('/home/yeray/data/atom_to_reactor_private/data/lab_observables.json')
TARGET = HERE / 'data_dir' / 'data' / 'lab_observables.json'
NMR_ROOT = Path('/home/yeray/data/nmr_261105')
NAME = '2 % H2O, 05/11/2026'

meta = json.loads(META_SOURCE.read_text(encoding='utf-8'))['_meta']
files = build_lab_nmr_inventory(NMR_ROOT)
p31 = files[files['nucleus'] == '31P'].sort_values(['acquired_at', 'file'])
first = p31.iloc[0]
water_added = first['acquired_at'] + pd.Timedelta(minutes=first['duration_min']) - pd.Timedelta(minutes=6, seconds=50)
p31 = p31[(p31['noe'] == 'FALSE') & (p31['relaxation_delay_s'] >= 5)].drop_duplicates('data_md5')
p31 = p31.assign(hours_since_first=(p31['acquired_at'] - p31['acquired_at'].min()).dt.total_seconds() / 3600.0)
shares = tabulate_area_shares(p31, P_WINDOWS_PPM, region_ppm=(-32.0, 6.0))
spectra = [_spectrum_record(rows, WINDOWS_OF_NUCLEUS['31P']) for _, rows in shares.groupby('acquired_at', sort=True)]
age_h = (p31['acquired_at'].min() - water_added).total_seconds() / 3600.0
T_rt = float(p31['temperature_C'].mean())
c0 = {sp: float(c) for sp, c in calculate_water_series_c0(2.0).items() if c > 0}
observables = {
    '_meta': {**meta, 'built_from': f'{NMR_ROOT} (kin_M1_2pct-261105_Esolv-thermo2/build_observables.py): the 5 Nov 2026 tube only'},
    'samples': {NAME: {
        'label': '5 vol% TMSPa + 2 vol% H2O in EC/DEC, 5 Nov 2026 (our spectra)',
        'campaign': 2026, 'nucleus': '31P', 'role': 'fit',
        'c0_M': c0,
        'recipe_basis': 'water series recipe: 2 vol% H2O stock diluted by 5 vol% TMSPa; EC 7.1 M',
        'first_spectrum_at': p31['acquired_at'].min().isoformat(),
        'T_rt_C': T_rt,
        'age': {'known': True, 'note': f'water added {water_added:%Y-%m-%d %H:%M:%S} (6 min 50 s before the end of the first 31P); '
                                       f'first 31P starts {60 * age_h:.2f} min later'},
        'pre_history': [{'T_C': T_rt, 'duration_h': age_h}],
        'heated': [], 'heated_note': 'none: every spectrum at room temperature',
        'replicates': False,
        'spectra': spectra,
        'flags': ['31P-2-1 holds the data of 31P-1-1 (same md5): read once',
                  '1H recorded (5 spectra), Si-CH3 region noise (S/N 3-7): not used; no 13C recorded'],
        'more_nuclei': {},
    }},
}
TARGET.write_text(json.dumps(observables, indent=2), encoding='utf-8')
print(f'{TARGET}: {NAME} only; c0 ' + ', '.join(f'{sp} {1000 * c:.0f} mM' for sp, c in c0.items()) + f'; {len(spectra)} 31P at '
      + ', '.join(f"{60 * (age_h + s['t_since_first_h']):.1f}" for s in spectra) + ' min since the water')
