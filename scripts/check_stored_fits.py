#!/usr/bin/env python3
"""Recompute the χ² of the stored fits of a results directory, without fitting, and compare it with the stored one.

    python scripts/check_stored_fits.py notebooks/results/fit8_M1-steps

Each fit (fits/<structure>__<reading>.json) and each step (steps/*.json) is run forward with the report solver,
in the setup of its run (fit_config.json), on the samples it was fitted to. A χ² that changed means the code no
longer runs that fit as it was fitted, or the lab observables changed: run it after any change to kinetics/.
The fits of a variant (a third part in the file name) are not checked. The directory is only read.

Source: Y. Alcaraz Galván
"""

import argparse
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from kinetics.data import load_lab_observables
from kinetics.fitting import REPORT_SOLVER, FitProblem, apply_fit_config, build_sample_set, load_fit_result

LAB_SAMPLES = ('0.5 % H2O', '2 % H2O', 'TMSPa + TMSOH (A)', 'TMSPa alone', 'TMSOH, probe', 'TMSOH, glovebox')


def check_fit(path: Path, observables: dict) -> float:
    """Relative change of the χ² of one stored fit (None: its structure cannot be rebuilt by name)."""
    fit = load_fit_result(str(path))
    exclude = tuple(s for s in LAB_SAMPLES if s not in fit['step_fit']) if 'step_fit' in fit else ()
    samples = build_sample_set(observables, fit['scenario'], exclude=exclude)
    try:
        problem = FitProblem(fit['structure'], samples, solver=REPORT_SOLVER, fixed=fit.get('fixed') or None)
    except KeyError:
        return None
    chi2 = problem.chi2(problem.vector(fit['theta']))
    return abs(chi2 - fit['chi2']) / max(1.0, abs(fit['chi2']))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('results', help='results directory, with its fit_config.json')
    parser.add_argument('--tolerance', type=float, default=1e-9, help='largest relative change accepted')
    args = parser.parse_args()
    results = Path(args.results).resolve()
    apply_fit_config(results)
    observables = load_lab_observables()
    paths = sorted(p for p in (results / 'fits').glob('*.json') if p.stem.count('__') == 1)
    paths += sorted((results / 'steps').glob('*.json'))
    t0, worst, changed, skipped = time.perf_counter(), 0.0, [], []
    for path in paths:
        change = check_fit(path, observables)
        if change is None:
            skipped.append(path.name)
            continue
        worst = max(worst, change)
        if change > args.tolerance:
            changed.append(path.name)
            print(f'  CHANGED {path.name}: relative change {change:.1e}')
    for name in skipped:
        print(f'  not checked {name}: its structure cannot be rebuilt by name')
    print(f'{results.name}: {len(paths) - len(skipped)} stored fits recomputed, largest relative change {worst:.1e}, '
          f'{len(changed)} changed ({time.perf_counter() - t0:.0f} s)')
    sys.exit(1 if changed else 0)


if __name__ == '__main__':
    main()
