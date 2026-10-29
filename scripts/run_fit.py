#!/usr/bin/env python3
"""Staged, resumable runner of M1 fitted to the lab NMR observables one tube at a time.

    python scripts/run_fit.py steps --results DIR   the main structure fitted one tube at a time, each step scoring every tube

Every job writes one file under the results directory (--results, required) and is skipped when that file exists,
so a stage can be interrupted and started again. A new fit needs a new directory: in one that holds results the old
fits are returned as they are.

Setup of a run (kinetics/fitting/config.py): `fit_config.json` in the results directory, read at the start of every
stage, fixes the solvation medium, the network, the run-wide changes to the samples, the readings, the main
structure, the search effort and the steps. Each fit records them, so a notebook runs it forward as it was fitted.
A structure fits a reading if no standardised residual exceeds 3.

Steps (eighth fit, notebooks/results/fit8_M1-steps/fit_config.json; MET 2026-10-23: one experiment at a time): the
steps stage fits the main structure to one tube at a time ('steps': the tubes, the free barriers and the step the
others come from; a barrier no earlier step set is held off), and after every step scores every tube with that
step's parameters, so the tubes not used yet are predicted ('step_variants': other readings of a tube's history,
scored the same way).

The stages of the fits of every tube at once (throughput, recovery, baselines, extended, refit, night, traces,
predictions) are in the branch M5-full.

Source: Y. Alcaraz Galván
"""

import os

for _var in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_var, '1')            # one BLAS thread per worker process

import argparse
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from kinetics.data import load_lab_observables
from kinetics.data.medium import get_active_medium
from kinetics.data.network import get_active_network_variant
from kinetics.fitting import (
    REPORT_SOLVER, SCENARIOS, FitProblem, apply_fit_config, build_sample_set, calculate_profile,
    find_confidence_interval, fit_structure, get_structure, load_fit_result, tabulate_reaction_barriers,
    write_fit_result,
)
from kinetics.fitting.candidates import get_energy_box_sigma
from kinetics.fitting.residuals import get_sample_overrides

RESULTS = None                           # the results directory (--results); set in main

# Search effort per structure: (Sobol screen points, starts polished, least-squares iterations per start); the
# setup can change it ('effort')
EFFORT = {'M1': (128, 4, 30), 'M1-split': (128, 4, 30)}


def log(message: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {message}"
    print(line, flush=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    with open(RESULTS / 'run.log', 'a', encoding='utf-8') as f:
        f.write(line + '\n')


def cached(path: Path, compute):
    """Result stored at `path`, computed and written when it is not there yet."""
    if path.exists():
        return load_fit_result(str(path))
    result = compute()
    write_fit_result(result, str(path))
    return result


def sample_set(scenario: str, **options) -> list:
    return build_sample_set(load_lab_observables(), scenario, **options)


def short(theta: dict) -> str:
    return ', '.join(f"{k.replace('g ', '').replace('hydrolysis_', '')} {v:.3f}" for k, v in theta.items())


# ------------------------------------------------------------------------------
# steps: the main structure fitted one tube at a time
# ------------------------------------------------------------------------------
SCENARIO_ORDER = ('10 min', '1 h', '8 h')
MAIN_STRUCTURE = 'M1'                    # the structure the steps fit; set in main
LAB_FITTED_SAMPLES = ('0.5 % H2O', '2 % H2O', 'TMSPa + TMSOH (A)', 'TMSPa alone', 'TMSOH, probe', 'TMSOH, glovebox')


def base_name(structure: str) -> str:
    return structure.split('+')[0].split(' (')[0]


STEPS = ()                               # the steps of the run (fit_config.json 'steps'); set in main
STEP_VARIANTS = {}                       # other readings of a tube's history, scored after every step ('step_variants')


def step_path(step: str, scenario: str, folder: str = 'steps') -> Path:
    return RESULTS / folder / f'{MAIN_STRUCTURE}__{step}__{scenario}.json'.replace(' ', '_')


def score_step(problem: FitProblem, theta: dict) -> dict:
    """A step without free parameters: its tubes scored with the parameters it takes over, reported as
    fit_structure reports a fit (report solver)."""
    final = problem.replace(solver=REPORT_SOLVER)
    x = final.vector(theta)
    table, residuals = final.table(x), final.residuals(x)
    worst = table.loc[table['z'].abs().idxmax()]
    return {'structure': problem.structure.name, 'parameters': [], 'fixed': dict(problem.fixed), 'theta': final.theta(x),
            'x': [], 'chi2': float(residuals @ residuals), 'n_residuals': int(residuals.size), 'n_parameters': 0,
            'max_abs_z': float(abs(worst['z'])), 'worst': f"{worst['sample']}: {worst['quantity']}",
            'fits': bool(abs(worst['z']) <= 3.0), 'table': table, 'ages': {}, 'local_optima': [], 'at_bounds': [],
            'n_evaluations': 1, 'seconds': 0.0}


def score_tubes(theta: dict, scenario: str, tubes, overrides: dict = None) -> list:
    """χ² (whitened as in the fit, so the tubes add up to the χ² of all of them) and worst miss of each tube under
    one parameter set; nothing is fitted."""
    rows = []
    for sample in sample_set(scenario, overrides=overrides or {}):
        if sample['name'] not in tubes:
            continue
        problem = FitProblem(MAIN_STRUCTURE, [sample], fixed=theta, solver=REPORT_SOLVER)
        table = problem.table(np.zeros(0))
        worst = table.loc[table['z'].abs().idxmax()]
        rows.append({'tube': sample['name'], 'chi2': problem.chi2(np.zeros(0)), 'max |z|': float(abs(worst['z'])),
                     'worst': worst['quantity'], 'fits': bool(abs(worst['z']) <= 3.0)})
    return rows


def stage_steps(args) -> None:
    """The main structure fitted one tube at a time, in the order of the setup's steps (eighth fit). A step fits
    the parameters of 'free' to the tubes of 'fit' alone; every other parameter comes from the step named in 'from'
    (same reading), or is held off when no earlier step set it (the top of its box: the reaction does not run in
    these tubes). A step without free parameters scores its tubes with the parameters it takes over. After every
    step each lab tube is scored with the step's parameters, and the tubes of 'step_variants' also under their
    other history: the tubes not used yet are predicted, so each step tests the ones before it."""
    structure = get_structure(MAIN_STRUCTURE)
    off = {p.name: p.upper for p in structure.parameters(sample_set(SCENARIO_ORDER[0]))}
    n_screen, n_starts, max_nfev = EFFORT[base_name(MAIN_STRUCTURE)]
    by_id = {step['id']: step for step in STEPS}
    summary, scores, intervals = [], [], []
    for scenario in SCENARIO_ORDER:
        thetas = {}
        for step in STEPS:
            name, free = step['id'], list(step['free'])
            base = thetas[step['from']] if step.get('from') else off
            data = sample_set(scenario, exclude=tuple(s for s in LAB_FITTED_SAMPLES if s not in step['fit']))
            problem = FitProblem(structure, data, fixed={k: v for k, v in base.items() if k not in free})

            def compute():
                if free:
                    fit = fit_structure(problem, n_screen=max(int(n_screen * args.scale), 0), n_starts=n_starts,
                                        max_nfev=max_nfev, seed=args.seed, workers=args.workers, start_points=[base])
                else:
                    fit = score_step(problem, base)
                fit.update(scenario=scenario, tag=f'__{name}', step=name, step_fit=list(step['fit']), step_free=free,
                           step_from=step.get('from'), label=structure.label, energy_box_sigma=get_energy_box_sigma(),
                           barriers=tabulate_reaction_barriers(structure, fit['theta']), medium=get_active_medium().to_dict(),
                           network=get_active_network_variant(), sample_overrides=get_sample_overrides(),
                           variant_overrides={})
                return fit

            fresh = not step_path(name, scenario).exists()
            fit = cached(step_path(name, scenario), compute)
            thetas[name] = fit['theta']
            log(f"step {name:4s} {scenario:6s} {', '.join(step['fit']):20s} chi2 {fit['chi2']:9.2f}  max |z| {fit['max_abs_z']:6.2f}  "
                f"{'FITS' if fit['fits'] else 'fails'}  [{short({k: fit['theta'][k] for k in free}) or 'nothing free'}]"
                + (f"  ({fit['n_evaluations']} evaluations, {fit['seconds']:.0f} s)" if fresh else '  (cached)'))
            summary.append({'step': name, 'scenario': scenario, 'tubes fitted': ', '.join(step['fit']), 'free': ', '.join(free),
                            'from': step.get('from') or 'off', 'chi2': fit['chi2'], 'max |z|': fit['max_abs_z'],
                            'worst': fit['worst'], 'fits': fit['fits'], 'at bounds': ', '.join(fit['at_bounds']), **fit['theta']})
            if free:
                bundle = cached(step_path(name, scenario, 'step_profiles'), lambda: {
                    'profile': calculate_profile(problem, fit, free, max_nfev=args.profile_nfev, workers=args.workers)})
                intervals += [{**find_confidence_interval(bundle['profile'], p), 'step': name, 'scenario': scenario} for p in free]
            # the tubes this step's parameters rest on: its own and those of the steps it takes them from
            used, link = set(step['fit']), step.get('from')
            while link:
                used |= set(by_id[link]['fit'])
                link = by_id[link].get('from')
            for variant, overrides in {'as fitted': {}, **STEP_VARIANTS}.items():
                tubes = LAB_FITTED_SAMPLES if not overrides else tuple(overrides)
                for row in score_tubes(fit['theta'], scenario, tubes, overrides):
                    role = 'fitted here' if row['tube'] in step['fit'] else 'fitted before' if row['tube'] in used else 'not used'
                    scores.append({'step': name, 'scenario': scenario, 'variant': variant, 'role': role, **row})
    pd.DataFrame(summary).to_csv(RESULTS / 'steps_summary.csv', index=False)
    pd.DataFrame(scores).to_csv(RESULTS / 'steps_scores.csv', index=False)
    pd.DataFrame(intervals).to_csv(RESULTS / 'steps_intervals.csv', index=False)
    predicted = pd.DataFrame(scores).query("role == 'not used' and variant == 'as fitted'")
    log('predicted tubes, largest |z| per step and reading: ' + '; '.join(
        f"{st} {sc} {g['max |z|'].max():.1f}" for (st, sc), g in predicted.groupby(['step', 'scenario'], sort=False)))


STAGES = {'steps': stage_steps}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('stage', choices=list(STAGES))
    parser.add_argument('--results', required=True, help='results directory, with its fit_config.json')
    parser.add_argument('--workers', type=int, default=12)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--scale', type=float, default=1.0, help='multiplies every Sobol screen (smoke tests: 0.1)')
    parser.add_argument('--profile-nfev', type=int, default=15, help='least-squares iterations per profile point')
    parser.add_argument('--quick', action='store_true', help='smoke test: tiny screens, two starts, few iterations')
    args = parser.parse_args()
    global RESULTS, SCENARIOS, SCENARIO_ORDER, MAIN_STRUCTURE, STEPS, STEP_VARIANTS
    RESULTS = Path(args.results).resolve()
    config = apply_fit_config(RESULTS)
    SCENARIOS = tuple(sc for sc in SCENARIOS if sc in config['scenarios'])
    SCENARIO_ORDER = tuple(sc for sc in SCENARIO_ORDER if sc in config['scenarios'])
    if config['structures']:
        MAIN_STRUCTURE = config['structures']['main']
    EFFORT.update({name: tuple(effort) for name, effort in config['effort'].items()})
    STEPS, STEP_VARIANTS = tuple(config['steps']), dict(config['step_variants'])
    for i, step in enumerate(STEPS):                 # a step takes its parameters from an earlier one
        if step.get('from') and step['from'] not in [earlier['id'] for earlier in STEPS[:i]]:
            raise SystemExit(f"{config['source']}: step {step['id']} takes its parameters from '{step['from']}', not an earlier step")
        unknown = set(step['fit']) - set(LAB_FITTED_SAMPLES)
        if unknown:
            raise SystemExit(f"{config['source']}: step {step['id']} fits unknown tubes {sorted(unknown)}")
    if args.quick:
        args.scale, args.profile_nfev = 0.05, 3
        for name in EFFORT:
            EFFORT[name] = (EFFORT[name][0], 2, 5)
    warnings.filterwarnings('ignore', category=RuntimeWarning)
    RESULTS.mkdir(parents=True, exist_ok=True)
    log(f'--- stage {args.stage} (workers {args.workers}) ---')
    log(f"setup ({config['source']}): {config['medium'].label()}; network {config['network']}; energy box "
        f"{'±%g standard errors' % config['energy_box_sigma'] if config['energy_box_sigma'] else 'fixed'}; main structure "
        f"{MAIN_STRUCTURE}; readings "
        f"{', '.join(SCENARIOS)}; samples changed for the run: "
        f"{', '.join(config['sample_overrides']) or 'none'}")
    t0 = time.perf_counter()
    STAGES[args.stage](args)
    log(f'--- stage {args.stage} done in {(time.perf_counter() - t0) / 60:.1f} min ---')


if __name__ == '__main__':
    main()
