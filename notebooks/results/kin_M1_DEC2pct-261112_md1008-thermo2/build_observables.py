"""The observables of this run: ONE tube, the 2 % H2O in neat DEC of 12 Nov 2026 (our spectra, fit_thermo_rev01 §3).

    python notebooks/results/kin_M1_DEC2pct-261112_md1008-thermo2/build_observables.py

- ³¹P only: the five spectra with the title 'TMSPa 2%water in DEC' (4 scans, NOE off, d1 5 s), read with the windows of
  neat DEC of fit_thermo_rev01 §3 (BMSPA and MMSPA 1.5 ppm above their EC/DEC place; region −32…6 ppm).
- ¹H not used (ringing, silyl windows read only roughly, §3); no ¹³C recorded. The spectra of 'TMSPa with DEC' (before the
  water) are another tube: not used.
- Times known: water added at 14:21:23, 1 min before the start of the first ³¹P (§3) = t 0; the time to the first ³¹P is a
  known room-temperature history ('pre_history'). Spectrum times: the start of each acquisition.
- Composition at mixing: the volumes of §3 (25 µL TMSPa + 465 µL DEC + 10 µL H2O; densities of the model): TMSPA 149 mM,
  H2O 1107 mM, no EC. The trace of BMSPA before the water (3.6 ± 1.8 mM, §3) is left out.
- 'replicates': False — the tube reacts while recorded.
- TMSPA window of the spectra from 5 min on: σ 1e3 (no weight). The TMSPA line turns dispersive after the water (no tuning or
  shimming) and its area dips below the baseline: 15, −15, −5 mM after 5 mM at 3 min.
"""

import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))

from kinetics.data.lab_nmr import build_lab_nmr_inventory, tabulate_area_shares   # noqa: E402
from kinetics.data.observables import WINDOWS_OF_NUCLEUS, _spectrum_record   # noqa: E402
from kinetics.data.snapshot import calculate_molar_mass   # noqa: E402
from kinetics.reactor.protocol import DEFAULT_DENSITY_G_ML   # noqa: E402

META_SOURCE = Path('/home/yeray/data/atom_to_reactor_private/data/lab_observables.json')
TARGET = HERE / 'data_dir' / 'data' / 'lab_observables.json'
NMR_ROOT = Path('/home/yeray/data/nmr_261112')
TITLE = 'TMSPa 2%water in DEC'
NAME = '2 % H2O in DEC, 12/11/2026'
DEC_P_WINDOWS = {'TMSPA': (-27.0, -22.0), 'BMSPA': (-15.5, -11.5), 'MMSPA': (-7.0, -3.5), 'H3PO4': (-1.0, 5.0)}
V_UL = {'TMSPA': 25.0, 'DEC': 465.0, 'H2O': 10.0}
TMSPA_OUT_FROM_MIN = 4.0          # TMSPA windows of the spectra after this time: no weight (dispersive line)

meta = json.loads(META_SOURCE.read_text(encoding='utf-8'))['_meta']
files = build_lab_nmr_inventory(NMR_ROOT)
p31 = files[(files['nucleus'] == '31P') & (files['title'] == TITLE)].sort_values(['acquired_at', 'file'])
water_added = p31['acquired_at'].min() - pd.Timedelta(minutes=1)
p31 = p31[(p31['noe'] == 'FALSE') & (p31['relaxation_delay_s'] >= 5)].drop_duplicates('data_md5')
p31 = p31.assign(hours_since_first=(p31['acquired_at'] - p31['acquired_at'].min()).dt.total_seconds() / 3600.0)
shares = tabulate_area_shares(p31, DEC_P_WINDOWS, region_ppm=(-32.0, 6.0))
spectra = [_spectrum_record(rows, WINDOWS_OF_NUCLEUS['31P']) for _, rows in shares.groupby('acquired_at', sort=True)]
age_h = (p31['acquired_at'].min() - water_added).total_seconds() / 3600.0
dropped = []
for k, s in enumerate(spectra):
    t_min = 60.0 * (age_h + s['t_since_first_h'])
    if t_min > TMSPA_OUT_FROM_MIN:
        v = s['shares']['TMSPA']
        v['sigma_noise_measured'] = v['sigma_noise']
        v['sigma_noise'] = 1.0e3
        dropped.append(f"TMSPA at {t_min:.1f} min (share {v['share']:+.3f})")
T_rt = float(p31['temperature_C'].mean())
total_uL = sum(V_UL.values())
c0 = {sp: 1000.0 * V_UL[sp] * DEFAULT_DENSITY_G_ML[sp] / calculate_molar_mass(sp) / total_uL for sp in ('TMSPA', 'H2O')}   # mol/L
observables = {
    '_meta': {**meta, 'built_from': f'{NMR_ROOT} ({HERE.name}/build_observables.py): the 12 Nov 2026 neat-DEC tube only'},
    'samples': {NAME: {
        'label': '5 vol% TMSPa + 2 % H2O (of the sample) in neat DEC, 12 Nov 2026 (our spectra)',
        'campaign': 2026, 'nucleus': '31P', 'role': 'fit',
        'c0_M': c0,
        'recipe_basis': '25 µL TMSPa + 465 µL DEC + 10 µL H2O (2 % of the 500 µL sample, to confirm), volumes adding up, no EC',
        'first_spectrum_at': p31['acquired_at'].min().isoformat(),
        'T_rt_C': T_rt,
        'age': {'known': True, 'note': f'water added {water_added:%Y-%m-%d %H:%M:%S} (1 min before the first 31P); '
                                       f'first 31P starts {60 * age_h:.2f} min later'},
        'pre_history': [{'T_C': T_rt, 'duration_h': age_h}],
        'heated': [], 'heated_note': 'none: every spectrum at room temperature',
        'replicates': False,
        'spectra': spectra,
        'flags': ['31P windows of neat DEC (fit_thermo_rev01 §3): BMSPA −15.5…−11.5, MMSPA −7.0…−3.5, H3PO4 −1.0…5.0 ppm',
                  'no weight (σ 1e3), dispersive TMSPA line: ' + '; '.join(dropped),
                  '1H recorded (5 spectra), ringing: not used; no 13C recorded'],
        'more_nuclei': {},
    }},
}
TARGET.write_text(json.dumps(observables, indent=2), encoding='utf-8')
print(f'{TARGET}: {NAME} only; c0 ' + ', '.join(f'{sp} {1000 * c:.0f} mM' for sp, c in c0.items()) + f'; {len(spectra)} 31P at '
      + ', '.join(f"{60 * (age_h + s['t_since_first_h']):.1f}" for s in spectra) + ' min since the water; ' + '; '.join(dropped))
