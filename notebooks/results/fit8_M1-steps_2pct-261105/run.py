"""scripts/run_fit.py steps for this run, with the tube of 5 Nov 2026 among the lab tubes a step can fit.

    ATOM_DATA_DIR=notebooks/results/fit8_M1-steps_2pct-261105/data_dir \
        python notebooks/results/fit8_M1-steps_2pct-261105/run.py steps --results notebooks/results/fit8_M1-steps_2pct-261105
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / 'scripts'))
sys.path.insert(0, str(REPO))
import run_fit  # noqa: E402

run_fit.LAB_FITTED_SAMPLES = run_fit.LAB_FITTED_SAMPLES + ('2 % H2O, 05/11/2026',)
run_fit.main()
