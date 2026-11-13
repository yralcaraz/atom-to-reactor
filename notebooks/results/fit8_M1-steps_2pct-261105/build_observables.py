"""The observables of this run: those of fit8_M1-steps_S4-at-age (the last M1-steps run) plus the 2 % H2O tube of
5 Nov 2026 (our spectra, fit_thermo_rev01 §2) as a tube of its own, '2 % H2O, 05/11/2026'.

    python notebooks/results/fit8_M1-steps_2pct-261105/build_observables.py

- ³¹P only: the five quantitative spectra (NOE off, d1 5 s), read as every ³¹P of the observables (tabulate_area_shares,
  P_WINDOWS_PPM, region −32…6 ppm). 31P-2-1 holds the data of 31P-1-1 (same md5): read once.
- ¹H not used: its Si–CH3 region is noise (S/N 3–7, §2.2.2); no ¹³C was recorded.
- Times known: water added at 16:33:04 (6 min 50 s before the end of the first ³¹P, as in §2); the time between it
  and the start of the first ³¹P is a known room-temperature history ('pre_history'), so no reading of the age.
  Spectrum times: the start of each acquisition, as for every tube of the observables.
- Composition at mixing: that of the 2 % H2O tube (same recipe). Room temperature: the mean of the spectra.
- 'replicates': False — the tube reacts while recorded, so its spectra are not repeats of one composition.
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

SOURCE = REPO / 'notebooks' / 'results' / 'fit8_M1-steps_S4-at-age' / 'data_dir' / 'data' / 'lab_observables.json'
TARGET = HERE / 'data_dir' / 'data' / 'lab_observables.json'
NMR_ROOT = Path('/home/yeray/data/nmr_261105')
NAME = '2 % H2O, 05/11/2026'

observables = json.loads(SOURCE.read_text(encoding='utf-8'))
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
recipe = observables['samples']['2 % H2O']
observables['samples'][NAME] = {
    'label': '5 vol% TMSPa + 2 vol% H2O in EC/DEC, 5 Nov 2026 (our spectra)',
    'campaign': 2026, 'nucleus': '31P', 'role': 'fit',
    'c0_M': dict(recipe['c0_M']),
    'recipe_basis': recipe['recipe_basis'] + ' (same recipe as the 2 % H2O tube of 2022)',
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
}
observables['_meta']['built_from'] = (observables['_meta'].get('built_from', '') + f'; plus {NAME} from {NMR_ROOT} '
                                      '(fit8_M1-steps_2pct-261105/build_observables.py)')
TARGET.write_text(json.dumps(observables, indent=2), encoding='utf-8')
print(f'{TARGET}: {NAME}, {len(spectra)} 31P spectra at '
      + ', '.join(f"{60 * (age_h + s['t_since_first_h']):.1f}" for s in spectra) + ' min since the water')
