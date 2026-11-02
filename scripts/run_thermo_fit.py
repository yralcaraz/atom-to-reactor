#!/usr/bin/env python3
"""Staged, resumable runner of the thermodynamic fit: the solvation energies from the spectra at equilibrium, one tube at a time.

    python scripts/run_thermo_fit.py points --results DIR   every usable spectrum read, and which are at equilibrium (raw spectra)
    python scripts/run_thermo_fit.py steps  --results DIR   the steps of the setup at every ξ, each step scoring every tube
    python scripts/run_thermo_fit.py all    --results DIR   both

The model, the free parameters and the rule that picks the spectra are in kinetics/fitting/thermo_fit.py (notebook 06
§1 and §5: E_solv,i free within ±ξ·SD_i, χ² of the shares, no other constraint).

Setup of a run: `fit_config.json` in the results directory, as for scripts/run_fit.py (kinetics/fitting/config.py): the
solvation medium, the network and the run-wide changes to the samples, and under 'thermo' the part of this fit
(thermo_fit.DEFAULT_THERMO_SETUP: the ξ to fit, the free species, the reactions at equilibrium, the nuclei fitted and
checked, the rule, the files the setup decides, the search effort and the steps).

points: reads every .jdf of the raw folders (LAB_NMR_DIR) as notebook 05 rev08 §2 does and writes `spectra_shares.csv`
(one row per usable spectrum and window) and `equilibrium_points.csv` (one row per spectrum: at equilibrium or not,
and why). Both hold lab data: they stay in the results directory, never in git. Run again (after deleting them) when
the rule or 'overrides' change; the steps read them.

steps 'cumulative' (setup 'order'): one tube added at a time, every tube added so far fitted together, each step from the
one before; tubes of the same recipe that end apart (thermo_fit.compare_replicates, written to `replicates.csv`) split the
chain into two solutions (thermo_fit.build_cumulative_steps). Or a list of steps (as stage_steps of scripts/run_fit.py): a step fits the parameters of 'free' to the tubes of 'fit' alone ('all':
every tube with a spectrum at equilibrium), with every other parameter from the step named in 'from' (same ξ) or at its
computed value when no earlier step set it. 'free' is 'data' (the default): every species present in those tubes at
equilibrium and every amount at mixing they leave open, so nothing the spectra see is held; or a list of species names
and parameter names such as 'c0 H2O TMSPa alone'. After every step each tube at equilibrium is scored
with the step's parameters, on the fitted nuclei and on the checked ones. Every job writes one file and is skipped
when that file exists, so a stage can be interrupted and started again; a new setup needs a new directory.

Parallel run (run_parallel.sh in the results directory): the ξ and the chains are independent, so one process per ξ and
chain (--xi, --chain) fits its part and writes its fits; a last run without them finds every fit on disk and writes the
tables. Parts that share a step (the steps before a chain splits) are run by the first chain only.

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

from kinetics.data import build_lab_nmr_inventory, load_lab_observables
from kinetics.data.medium import get_active_medium
from kinetics.data.network import get_active_network_variant
from kinetics.fitting import apply_fit_config, calculate_profile, find_confidence_interval, fit_structure, load_fit_result, write_fit_result
from kinetics.fitting.candidates import RECIPE_PREFIX
from kinetics.fitting.residuals import get_sample_overrides
from kinetics.fitting.thermo_fit import (
    SOLVATION_PREFIX, ThermoProblem, build_cumulative_steps, build_equilibrium_samples, calculate_energy_ranges,
    compare_replicates, find_equilibrium_points,
    load_thermo_setup, tabulate_lab_spectra, tabulate_reaction_energies,
)
from kinetics.paths import lab_nmr_dir

RESULTS = None                           # the results directory (--results); set in main
SETUP = None                             # the 'thermo' part of its fit_config.json; set in main


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


# ------------------------------------------------------------------------------
# points: every usable spectrum, and which are at equilibrium
# ------------------------------------------------------------------------------
def stage_points(args) -> None:
    observables = load_lab_observables()
    if not (RESULTS / 'spectra_shares.csv').exists():
        log(f'reading every .jdf of {lab_nmr_dir()}')
        tabulate_lab_spectra(observables, build_lab_nmr_inventory(lab_nmr_dir())).to_csv(RESULTS / 'spectra_shares.csv', index=False)
    shares = read_shares()
    points = find_equilibrium_points(shares, rule=SETUP['rule'], overrides=SETUP['overrides'],
                                     error_floor=observables['_meta'].get('error_floor'))
    points.to_csv(RESULTS / 'equilibrium_points.csv', index=False)
    for tube, group in points.groupby('tube', sort=False):
        log(f"{tube}: " + '; '.join(f"{r['nucleus']} {r['t_h']:.1f} h {'AT EQUILIBRIUM' if r['at equilibrium'] else 'no'}"
                                    f"{' (setup)' if r['decided by'] == 'setup' else ''}" for _, r in group.iterrows()))


def read_shares() -> pd.DataFrame:
    return pd.read_csv(RESULTS / 'spectra_shares.csv', parse_dates=['acquired_at'])


def read_points() -> pd.DataFrame:
    path = RESULTS / 'equilibrium_points.csv'
    if not path.exists():
        raise SystemExit(f'{path} not found: run the stage points first')
    return pd.read_csv(path)


# ------------------------------------------------------------------------------
# steps: the free solvation energies fitted one tube at a time, at every ξ
# ------------------------------------------------------------------------------
def parameter_name(item: str) -> str:
    return f'{SOLVATION_PREFIX}{item}' if item in SETUP['free_species'] else item


def data_parameters(problem: ThermoProblem) -> list:
    """The parameters the tubes of a problem have data on: the solvation energy of every species present in one of them
    at equilibrium (computed energies), and their amounts at mixing that the recipe does not fix."""
    theta = start_values(problem)
    present = {sp for sample in problem.samples for sp, c in problem.equilibrium(sample, theta).items() if c > 0}
    return [p.name for p in problem.all_parameters if p.kind == 'recipe' or p.target in present]


def step_path(step: str, xi: float, folder: str = 'steps') -> Path:
    return RESULTS / folder / f'thermo__{step}__xi{xi:g}.json'


def start_values(problem: ThermoProblem) -> dict:
    """Every parameter at its computed value (shift 0) or, for an amount at mixing, in the middle of its range."""
    return {p.name: 0.0 if p.kind == 'solvation' else 0.5 * (p.lower + p.upper) for p in problem.all_parameters}


def score(problem: ThermoProblem, theta: dict) -> dict:
    """A problem without free parameters, reported as fit_structure reports a fit."""
    x = problem.vector(theta)
    table, residuals = problem.table(x), problem.residuals(x)
    worst = table.loc[table['z'].abs().idxmax()]
    return {'structure': 'thermo', 'parameters': [], 'fixed': dict(problem.fixed), 'theta': problem.theta(x), 'x': [],
            'chi2': float(residuals @ residuals), 'n_residuals': int(residuals.size), 'n_parameters': 0,
            'max_abs_z': float(abs(worst['z'])), 'worst': f"{worst['sample']}: {worst['nucleus']} {worst['quantity']}",
            'fits': bool(abs(worst['z']) <= 3.0), 'table': table, 'ages': {}, 'local_optima': [], 'at_bounds': [],
            'n_evaluations': 1, 'seconds': 0.0}


def score_tubes(samples: list, theta: dict, xi: float) -> list:
    """χ² and worst miss of each tube under one parameter set, on the fitted nuclei and on the checked ones."""
    rows = []
    for sample in samples:
        problem = ThermoProblem([sample], free_species=SETUP['free_species'], xi=xi, reactions=SETUP['reactions'])
        problem = problem.replace(fixed={**start_values(problem), **theta})     # an amount no step has fitted: mid-range
        for role in ('fit', 'check'):
            table = problem.table(np.zeros(0), roles=(role,))
            if table.empty:
                continue
            residuals = problem.residuals(np.zeros(0), roles=(role,))
            worst = table.loc[table['z'].abs().idxmax()]
            rows.append({'tube': sample['name'], 'nuclei': role, 'chi2': float(residuals @ residuals), 'n_residuals': int(residuals.size),
                         'max |z|': float(abs(worst['z'])), 'worst': f"{worst['nucleus']} {worst['quantity']}",
                         'fits': bool(abs(worst['z']) <= 3.0)})
    return rows


def stage_steps(args) -> None:
    """The steps of the setup at every ξ (as stage_steps of scripts/run_fit.py)."""
    observables = load_lab_observables()
    samples = build_equilibrium_samples(read_shares(), read_points(), observables, SETUP['nuclei'])
    tubes = [s['name'] for s in samples]
    contents = {s['name']: ', '.join(f"{b['nucleus']} × {len(b['t_spec_s'])} {b['role']}" for b in s['blocks']) for s in samples}
    log('tubes at equilibrium: ' + '; '.join(f'{tube} ({blocks})' for tube, blocks in contents.items()))
    steps = SETUP['steps']
    if steps == 'cumulative':
        replicates = compare_replicates(read_shares(), read_points(), observables)
        replicates.to_csv(RESULTS / 'replicates.csv', index=False)
        steps = build_cumulative_steps(SETUP['order'], replicates, tubes)
        log('replicates (same recipe): ' + ('; '.join(f"{r['tube']} / {r['other']} {r['nucleus']} max |z| {r['max |z|']:.1f} "
                                                      f"{'together' if r['compatible'] else 'APART: two solutions'}"
                                                      for _, r in replicates.iterrows()) or 'none'))
        log('steps: ' + '; '.join(f"{s['id']} [{', '.join(s['fit'])}]" + (f" from {s['from']}" if s.get('from') else '') for s in steps))
    if args.chain is not None:                         # one chain: the steps of that chain and the steps they come from
        by_id = {step['id']: step for step in steps}
        keep, todo = set(), [s['id'] for s in steps if s['id'].lstrip('T0123456789') == args.chain]
        if not todo:
            raise SystemExit(f"no step of chain '{args.chain}' (steps: {', '.join(by_id)})")
        while todo:
            step_id = todo.pop()
            keep.add(step_id)
            if by_id[step_id].get('from') and by_id[step_id]['from'] not in keep:
                todo.append(by_id[step_id]['from'])
        steps = [s for s in steps if s['id'] in keep]
    for step in steps:
        missing = [] if step['fit'] == 'all' else [t for t in step['fit'] if t not in tubes]
        if missing:
            raise SystemExit(f"step {step['id']} fits {missing}, which has no spectrum at equilibrium (equilibrium_points.csv)")
    n_screen, n_starts, max_nfev = SETUP['effort']
    by_id = {step['id']: step for step in steps}
    summary, scores, intervals, energies = [], [], [], []
    partial = args.xi is not None or args.chain is not None
    for xi in [x for x in SETUP['xi'] if args.xi is None or abs(x - args.xi) < 1e-12]:
        thetas = {}
        for step in steps:
            name, fitted = step['id'], tubes if step['fit'] == 'all' else list(step['fit'])
            data = [s for s in samples if s['name'] in fitted]
            problem = ThermoProblem(data, free_species=SETUP['free_species'], xi=xi, reactions=SETUP['reactions'])
            free = data_parameters(problem) if step.get('free', 'data') == 'data' else [parameter_name(p) for p in step['free']]
            base = {**start_values(problem), **(thetas[step['from']] if step.get('from') else {})}
            problem = problem.replace(fixed={k: v for k, v in base.items() if k not in free})

            def compute():
                if problem.n:
                    fit = fit_structure(problem, n_screen=max(int(n_screen * args.scale), 0), n_starts=n_starts,
                                        max_nfev=max_nfev, seed=args.seed, workers=args.workers, start_points=[base])
                else:
                    fit = score(problem, base)
                fit.update(xi=xi, step=name, step_fit=fitted, step_free=free, step_from=step.get('from'),
                           species_shifts_eV=problem.species_shifts(fit['theta']), reactions=list(SETUP['reactions']),
                           nuclei=SETUP['nuclei'], medium=get_active_medium().to_dict(), network=get_active_network_variant(),
                           sample_overrides=get_sample_overrides())
                return fit

            fresh = not step_path(name, xi).exists()
            fit = cached(step_path(name, xi), compute)
            thetas[name] = fit['theta']
            log(f"step {name:4s} ξ {xi:g}  {', '.join(fitted):30s} chi2 {fit['chi2']:9.2f}  max |z| {fit['max_abs_z']:6.2f}  "
                f"{'FITS' if fit['fits'] else 'fails'}  at bounds [{', '.join(fit['at_bounds'])}]"
                + (f"  ({fit['n_evaluations']} evaluations, {fit['seconds']:.0f} s)" if fresh else '  (cached)'))
            summary.append({'step': name, 'xi': xi, 'tubes fitted': ', '.join(fitted), 'free': ', '.join(free),
                            'from': step.get('from') or 'computed', 'chi2': fit['chi2'], 'n_residuals': fit['n_residuals'],
                            'max |z|': fit['max_abs_z'], 'worst': fit['worst'], 'fits': fit['fits'],
                            'at bounds': ', '.join(fit['at_bounds']), **fit['theta']})
            table = tabulate_reaction_energies(fit['species_shifts_eV'])
            if problem.n:
                bundle = cached(step_path(name, xi, 'step_profiles'), lambda: {
                    'profile': calculate_profile(problem, fit, problem.names, max_nfev=args.profile_nfev or SETUP['profile_nfev'],
                                                 workers=args.workers)})
                intervals += [{**find_confidence_interval(bundle['profile'], p), 'step': name, 'xi': xi} for p in problem.names]
                table = table.join(calculate_energy_ranges(bundle['profile'], problem))
            energies += [{'step': name, 'xi': xi, 'reaction': rxn, **row} for rxn, row in table.iterrows()]
            # the tubes this step's parameters rest on: its own and those of the steps it takes them from
            used, link = set(fitted), step.get('from')
            while link:
                used |= set(tubes if by_id[link]['fit'] == 'all' else by_id[link]['fit'])
                link = by_id[link].get('from')
            for row in score_tubes(samples, fit['theta'], xi):
                role = 'fitted here' if row['tube'] in fitted else 'fitted before' if row['tube'] in used else 'not used'
                scores.append({'step': name, 'xi': xi, 'role': role, **row})
    if partial:                                        # one part of a parallel run: its fits are on disk, the tables come
        log(f'part ξ {args.xi} chain {args.chain} done (tables: run steps without --xi and --chain)')   # from the whole run
        return
    pd.DataFrame(summary).to_csv(RESULTS / 'steps_summary.csv', index=False)
    pd.DataFrame(scores).to_csv(RESULTS / 'steps_scores.csv', index=False)
    pd.DataFrame(intervals).to_csv(RESULTS / 'steps_intervals.csv', index=False)
    pd.DataFrame(energies).to_csv(RESULTS / 'steps_energies.csv', index=False)


STAGES = {'points': stage_points, 'steps': stage_steps}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('stage', choices=list(STAGES) + ['all'])
    parser.add_argument('--results', required=True, help='results directory, with its fit_config.json')
    parser.add_argument('--workers', type=int, default=12)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--scale', type=float, default=1.0, help='multiplies every Sobol screen (smoke tests: 0.1)')
    parser.add_argument('--profile-nfev', type=int, default=None, help='least-squares iterations per profile point (default: the setup)')
    parser.add_argument('--xi', type=float, default=None, help='fit this ξ only (one part of a parallel run)')
    parser.add_argument('--chain', default=None, help="fit this chain only ('a', 'b': one part of a parallel run)")
    parser.add_argument('--quick', action='store_true', help='smoke test: tiny screens, two starts, few iterations')
    args = parser.parse_args()
    global RESULTS, SETUP
    RESULTS = Path(args.results).resolve()
    config = apply_fit_config(RESULTS)
    SETUP = load_thermo_setup(config)
    for i, step in enumerate([] if SETUP['steps'] == 'cumulative' else SETUP['steps']):   # a step takes its parameters from an earlier one
        if step.get('from') and step['from'] not in [earlier['id'] for earlier in SETUP['steps'][:i]]:
            raise SystemExit(f"{config['source']}: step {step['id']} takes its parameters from '{step['from']}', not an earlier step")
        free = step.get('free', 'data')
        unknown = [] if free == 'data' else [p for p in free if p not in SETUP['free_species'] and not p.startswith(RECIPE_PREFIX)]
        if unknown:
            raise SystemExit(f"{config['source']}: step {step['id']} frees {unknown}: neither a free species nor an amount at mixing")
    if args.quick:
        args.scale, args.profile_nfev = 0.05, 3
        SETUP['effort'] = [SETUP['effort'][0], 2, 5]
    warnings.filterwarnings('ignore', category=RuntimeWarning)
    RESULTS.mkdir(parents=True, exist_ok=True)
    log(f'--- stage {args.stage} (workers {args.workers}) ---')
    log(f"setup ({config['source']}): {config['medium'].label()}; network {config['network']}; reactions at equilibrium "
        f"{', '.join(SETUP['reactions'])}; free {', '.join(SETUP['free_species'])} within ±ξ standard errors, ξ "
        f"{', '.join(f'{x:g}' for x in SETUP['xi'])}; nuclei fitted {', '.join(SETUP['nuclei']['fit'])}, checked "
        f"{', '.join(SETUP['nuclei']['check'])}; rule {SETUP['rule']}; samples changed for the run: "
        f"{', '.join(config['sample_overrides']) or 'none'}")
    t0 = time.perf_counter()
    for stage in (list(STAGES) if args.stage == 'all' else [args.stage]):
        STAGES[stage](args)
    log(f'--- stage {args.stage} done in {(time.perf_counter() - t0) / 60:.1f} min ---')


if __name__ == '__main__':
    main()
