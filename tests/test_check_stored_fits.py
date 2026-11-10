"""Regression check: scripts/check_stored_fits.py recomputes the χ² of small stored fits exactly.

    python -m pytest tests/test_check_stored_fits.py      (needs the data folder, ATOM_DATA_DIR; skipped without it)

One kinetic run (fit8_M1-steps_S1-tube, 21 fits) and one thermodynamic run (thermo2_2pctH2O, 3 steps): a χ² that
moves means the code no longer runs a stored fit as it was fitted.

Source: Y. Alcaraz Galván
"""

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from kinetics.paths import data_dir

RUNS = ['fit8_M1-steps_S1-tube', 'thermo2_2pctH2O']


@pytest.mark.parametrize('run', RUNS)
def test_stored_fits_unchanged(run):
    results = REPO / 'notebooks' / 'results' / run
    if data_dir() is None or not data_dir().is_dir():
        pytest.skip('no data folder (ATOM_DATA_DIR)')
    if not results.is_dir():
        pytest.skip(f'{run} not on disk')
    done = subprocess.run([sys.executable, str(REPO / 'scripts' / 'check_stored_fits.py'), str(results)],
                          capture_output=True, text=True)
    assert done.returncode == 0, done.stdout + done.stderr
    assert 'largest relative change 0.0e+00' in done.stdout, done.stdout
