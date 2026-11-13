"""Error model of this run (2026-11-09, the user: weight the precise points more, forget the noisy ones), applied to the
observables written by build_observables.py:

    python notebooks/results/<run>/build_observables.py && python notebooks/results/<run>/weights.py

- no declared floor: σ of a share = noise ⊕ baseline of its spectrum (the objective of fit_thermo_rev01 §2.1). The floor
  of the 2022 fits (³¹P: 0.01 absolute or 5 % of the share) set H3PO4 at 5–7 % where its spectra give 1–4 %.
- NOISE_OUT (run …_nofloor-noiseout): a window whose share is below 2 × (noise ⊕ baseline) (S/N < 2) is noise: its σ is
  set to 1e3, so it has no weight in χ² (it stays in the table, z ≈ 0).
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
NOISE_OUT = HERE.name.endswith('noiseout')
path = HERE / 'data_dir' / 'data' / 'lab_observables.json'
obs = json.loads(path.read_text(encoding='utf-8'))
obs['_meta']['error_floor'] = {**obs['_meta']['error_floor'], '31P': {'absolute': 0.0, 'relative': 0.0}}
dropped = []
if NOISE_OUT:
    for name, record in obs['samples'].items():
        for k, spectrum in enumerate(record['spectra']):
            for window, v in spectrum['shares'].items():
                sigma = float(np.hypot(v['sigma_noise'], v['sigma_baseline']))
                if abs(v['share']) < 2.0 * sigma:
                    v['sigma_noise_measured'] = v['sigma_noise']
                    v['sigma_noise'] = 1.0e3
                    dropped.append(f"{window} spectrum {k} (share {v['share']:+.3f}, σ {sigma:.3f})")
        record.setdefault('flags', []).append('error model: no floor' + ('; windows with S/N < 2 weighted 0: ' + '; '.join(dropped) if dropped else ''))
path.write_text(json.dumps(obs, indent=2), encoding='utf-8')
print(f'{path}: no 31P floor' + (f"; weight 0 for {len(dropped)} windows: " + ', '.join(dropped) if NOISE_OUT else ''))
