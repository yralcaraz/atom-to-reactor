#!/usr/bin/env python3
"""Run every test script and report. Every test needs the Tank snapshot (the species database and the
network are built from it, even for synthetic data), so all are skipped when it is absent.

    python tests/run_all.py

The data folder is given by ATOM_DATA_DIR (see kinetics/paths.py and .env.example). Inside a test,
the checks that need further data files (the legacy baseline, the lab observables, the raw spectra) skip
themselves when the file is missing.

Source: Y. Alcaraz Galván
"""

import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from kinetics.data.snapshot import DEFAULT_SNAPSHOT_PATH
from kinetics.paths import DATA_DIR_VAR


def main() -> int:
    have_snapshot = os.path.exists(DEFAULT_SNAPSHOT_PATH)
    results = {}
    for test in sorted(Path(__file__).parent.glob('test_*.py')):
        if not have_snapshot:
            results[test.name] = 'skipped (needs the Tank snapshot)'
            continue
        done = subprocess.run([sys.executable, str(test)], cwd=REPO, capture_output=True, text=True)
        results[test.name] = 'ok' if done.returncode == 0 else f'FAILED (exit {done.returncode})\n{done.stdout[-600:]}{done.stderr[-1200:]}'
    width = max(len(name) for name in results)
    for name, status in results.items():
        print(f'{name:<{width}}  {status}')
    if not have_snapshot:
        print(f'\n{DATA_DIR_VAR} is not set or has no data/tank_api_snapshot.json: see .env.example')
    return int(any(status.startswith('FAILED') for status in results.values()))


if __name__ == '__main__':
    sys.exit(main())
