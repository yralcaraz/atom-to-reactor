#!/usr/bin/env python3
"""Staged, resumable runner of the first fit of the kinetic models to the lab NMR observables.

    python scripts/run_fit.py registered     the four registered models as they are, four scenarios (seconds)
    python scripts/run_fit.py throughput     evaluations per second of the full objective on this machine
    python scripts/run_fit.py recovery       known parameters refitted from synthetic shares (the phase 3 check)
    python scripts/run_fit.py baselines      M0 (capped BEP, Marcus) and M1, four scenarios
    python scripts/run_fit.py extended       M1-split, M3 and M3-split, four scenarios
    python scripts/run_fit.py refit          fits of M3 and M3-split checked at the tightened search tolerance, redone if needed
    python scripts/run_fit.py night          profiles, leave-one-sample-out, indistinguishability, sensitivities
    python scripts/run_fit.py traces         the main structure with a trace of TMSOH at mixing: as fitted, and fitted again
    python scripts/run_fit.py predictions    hold-out and checks, storage at 25 °C, Peter's protocol with its ΔS‡ band

Every job writes one file under notebooks/results/05/ and is skipped when that file exists, so a stage can be
interrupted and started again. The results derive from confidential lab data and are not stored in the repository.

Rules fixed before any fit (see the plan of the first fit, an internal document): scenarios 'short' (1 h), 'middle' (1 d), 'long' (7 d)
and 'free' (each unknown age between 1 h and 7 d); a structure fits a scenario if no standardised residual
exceeds 3; the first consistent parameter set is that of the smallest structure that fits, a common-age scenario
taking precedence over free ages.

Classification: PUBLIC (see CLASSIFICATION.md)
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

from kinetics.data import load_experimental_data, load_lab_observables
from kinetics.fitting import (
    DELTA_CHI2_95, REPORT_SOLVER, SCENARIOS, FitProblem, FitStructure, add_profile_points, build_sample_set,
    calculate_parameter_directions, extend_profile,
    calculate_profile, calculate_residuals, compare_structures, embed_parent_theta, find_confidence_interval,
    find_free_age, fit_leave_one_out, fit_structure, get_profile_theta, get_structure, load_fit_result,
    calculate_prediction_band, refine_profile_edges, simulate_synthetic_shares, summarize_residuals,
    tabulate_reaction_barriers, tabulate_residuals, write_fit_result,
)
from kinetics.microkinetics.models import MODELS, get_model
from kinetics.reactor import (
    build_protocol_schedule, calculate_recipe_molarities, calculate_remaining_fraction, find_crossing_time,
    simulate_history, simulate_protocol,
)

RESULTS = REPO / 'notebooks' / 'results' / '05'

# Search effort per structure: (Sobol screen points, starts polished, least-squares iterations per start)
EFFORT = {'M0-BEP': (128, 6, 30), 'M0-Marcus': (128, 6, 30), 'M1': (512, 12, 40), 'M1-split': (512, 12, 40),
          'M3': (2048, 16, 40), 'M3-split': (2048, 16, 40)}


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
# registered: the models of the registry as they are
# ------------------------------------------------------------------------------
def stage_registered(args) -> None:
    rows, tables = [], []
    for name in MODELS:
        for scenario in SCENARIOS:
            samples, model = sample_set(scenario), get_model(name)
            ages = find_free_age(model, samples, solver=REPORT_SOLVER)
            residuals = calculate_residuals(model, samples, ages=ages, solver=REPORT_SOLVER)
            table = tabulate_residuals(model, samples, ages=ages)
            rows.append({'model': name, 'scenario': scenario, **summarize_residuals(residuals, table)})
            tables.append(table.assign(model=name, scenario=scenario))
            log(f"registered {name:16s} {scenario:7s} chi2 {rows[-1]['chi2']:12.1f}  max |z| {rows[-1]['max_abs_z']:8.1f}  "
                f"({rows[-1]['worst']})")
    pd.DataFrame(rows).to_csv(RESULTS / 'registered.csv', index=False)
    pd.concat(tables, ignore_index=True).to_csv(RESULTS / 'registered_residuals.csv', index=False)


# ------------------------------------------------------------------------------
# throughput: what one night can hold
# ------------------------------------------------------------------------------
def stage_throughput(args) -> None:
    from scipy.stats import qmc
    from kinetics.fitting.estimation import _chi2_unit, _map
    rows = []
    for structure, scenario in (('M3-split', 'middle'), ('M3-split', 'free'), ('M1', 'middle')):
        problem = FitProblem(structure, sample_set(scenario))
        points = list(1.0 + qmc.Sobol(problem.n, scramble=True, seed=11).random(256))
        t0 = time.perf_counter()
        [problem.chi2(problem.from_unit(u)) for u in points[:24]]
        single = 24 / (time.perf_counter() - t0)
        t0 = time.perf_counter()
        _map(_chi2_unit, points, problem, args.workers)
        parallel = len(points) / (time.perf_counter() - t0)
        rows.append({'structure': structure, 'scenario': scenario, 'parameters': problem.n,
                     'evaluations per s, 1 process': single, f'evaluations per s, {args.workers} processes': parallel})
        log(f'throughput {structure} {scenario}: {single:.1f} /s single, {parallel:.1f} /s with {args.workers} processes')
    pd.DataFrame(rows).to_csv(RESULTS / 'throughput.csv', index=False)


# ------------------------------------------------------------------------------
# recovery: can this pipeline find known parameters in data like ours?
# ------------------------------------------------------------------------------
RECOVERY_TRUTHS = {
    'A': ('M1-split', 'middle', {'g hydrolysis_R1': 1.28, 'g hydrolysis_R23': 1.05, 'g transfer': 1.05,
                                 'g condensation': 1.25, 'g solvent_attack': 1.315}, (1, 2, 3),
          ['g hydrolysis_R1', 'g hydrolysis_R23', 'g transfer', 'g condensation', 'g solvent_attack'], (512, 12, 40)),
    'B': ('M3-split', 'middle', {'g hydrolysis_R1': 1.30, 'g hydrolysis_R23': 1.00, 'g transfer': 0.90,
                                 'g condensation': 1.25, 'g solvent_attack': 1.315,
                                 'dG R1': 0.0, 'dG R2': 0.12, 'dG R3': 0.17, 'dG R4': 0.0}, (1, 2),
          ['g hydrolysis_R1', 'g hydrolysis_R23', 'g transfer', 'g condensation', 'g solvent_attack', 'dG R2', 'dG R3'],
          (1024, 8, 40)),
}


def stage_recovery(args) -> None:
    """The design, error model and scenario are those of the lab set; only the shares are synthetic."""
    rows = []
    for label, (structure, scenario, truth, seeds, profiled, effort) in RECOVERY_TRUTHS.items():
        for seed in seeds:
            def compute():
                samples = sample_set(scenario)
                data = simulate_synthetic_shares(structure, truth, samples, seed=seed)
                problem = FitProblem(structure, data)
                report = problem.replace(solver=REPORT_SOLVER)
                chi2_truth = report.chi2(report.vector(truth))
                fit = fit_structure(problem, n_screen=effort[0], n_starts=effort[1], max_nfev=effort[2], seed=seed,
                                    workers=args.workers)
                profile = calculate_profile(problem, fit, profiled, max_nfev=12, workers=args.workers)
                return {'truth': truth, 'chi2_truth': chi2_truth, 'fit': fit, 'profile': profile,
                        'evaluations': fit['n_evaluations'] + int(profile['nfev'].sum())}
            t0 = time.perf_counter()
            result = cached(RESULTS / 'recovery' / f'{label}_seed{seed}.json', compute)
            fit = result['fit']
            intervals = pd.DataFrame([find_confidence_interval(result['profile'], name) for name in profiled])
            intervals['truth'] = [truth[name] for name in profiled]
            log(f"recovery {label} seed {seed}: chi2 fit {fit['chi2']:.2f} vs truth {result['chi2_truth']:.2f}, "
                f"{result['evaluations']} evaluations, {time.perf_counter() - t0:.0f} s")
            for _, ci in intervals.iterrows():
                inside = ((not np.isfinite(ci['low'])) or ci['low'] - 1e-3 <= ci['truth']) and \
                         ((not np.isfinite(ci['high'])) or ci['truth'] <= ci['high'] + 1e-3)
                rows.append({'case': f'{label}, seed {seed}', 'structure': structure, 'parameter': ci['parameter'],
                             'truth': ci['truth'], 'fit': fit['theta'][ci['parameter']], 'low': ci['low'], 'high': ci['high'],
                             'status': ci['status'], 'truth inside': bool(inside),
                             'chi2 fit': fit['chi2'], 'chi2 truth': result['chi2_truth'],
                             'optimiser ok': bool(fit['chi2'] <= result['chi2_truth'] + 0.5)})
    summary = pd.DataFrame(rows)
    summary.to_csv(RESULTS / 'recovery_summary.csv', index=False)
    misses = int((~summary['truth inside']).sum())
    optimiser = bool(summary['optimiser ok'].all())
    log(f"recovery: {len(summary)} intervals, truth outside in {misses}; optimiser at least as good as the truth in "
        f"every case: {optimiser}")
    log('recovery status by parameter: ' + '; '.join(
        f"{name}: {', '.join(sorted(set(grp['status'])))}" for name, grp in summary.groupby('parameter', sort=False)))
    passed = optimiser and misses <= max(1, round(0.15 * len(summary)))
    log(f"RECOVERY {'PASSED' if passed else 'FAILED'} (criterion: optimiser ok in every case and at most 15 % of the "
        'truths outside their 95 % interval)')
    return passed


# ------------------------------------------------------------------------------
# Fits of the candidate structures
# ------------------------------------------------------------------------------
FIRST_PASS_SOLVER = {'method': 'BDF', 'rtol': 1e-6, 'atol': 1e-10, 'max_rhs_calls': 30000}   # search tolerance of the first pass
PARENTS = {'M1': ('M0-Marcus',), 'M1-split': ('M1',), 'M3': ('M1',), 'M3-split': ('M1-split', 'M3')}
SCENARIO_ORDER = ('middle', 'short', 'long', 'free')
SIZE_ORDER = ('M0-BEP', 'M0-Marcus', 'M1', 'M1-split', 'M3', 'M3-split')
MAIN_STRUCTURE = 'M3-split'          # the largest structure: every other one is a special case of it


def base_name(structure: str) -> str:
    return structure.split('+')[0].split(' (')[0]


def fit_path(structure: str, scenario: str, tag: str = '') -> Path:
    return RESULTS / 'fits' / f'{structure}__{scenario}{tag}.json'.replace(' ', '_')


def complete_start(problem: FitProblem, theta: dict):
    """A start for this problem from a parameter set that may lack its age or water fraction; None if it lacks more."""
    out = dict(theta)
    for p in problem.parameters:
        if p.name not in out:
            if p.kind == 'age':
                out[p.name] = p.upper
            elif p.kind == 'water':
                out[p.name] = 1.0
            else:
                return None
    return out


TMSOH_BLOCK = ('TMSOH, probe', 'TMSOH, glovebox', 'E1')        # samples that see condensation and solvent attack alone
BLOCK_PARAMETERS = ('g condensation', 'g solvent_attack')


def fit_tmsoh_block(structure: FitStructure, data: list, args) -> dict:
    """{barrier: value} of condensation and solvent attack from the TMSOH-only samples of a sample set.

    These samples hold no phosphate, so the other barriers do not act in them; they are held at mid-box (freed
    energies at their computed values). Used only to seed the screen of the full fit, where every parameter is free.
    """
    block = [sample for sample in data if sample['name'] in TMSOH_BLOCK]
    if structure.shared_barrier or len(block) < 2:
        return {}
    fixed = {p.name: (0.0 if p.kind == 'energy' else 1.0 if p.kind == 'water' else 0.5 * (p.lower + p.upper))
             for p in structure.parameters(block) if p.name not in BLOCK_PARAMETERS}
    fit = fit_structure(FitProblem(structure, block, fixed=fixed), n_screen=max(int(256 * args.scale), 16), n_starts=4,
                        max_nfev=30, seed=args.seed, workers=args.workers, report=False)
    return {name: fit['theta'][name] for name in BLOCK_PARAMETERS}


def run_one(structure, scenario: str, args, *, tag: str = '', starts=(), sample_options: dict = None,
            prior_scale: float = 1.0, effort: tuple = None, samples: list = None) -> dict:
    """One cached fit. `samples` replaces the lab sample set (synthetic shares)."""
    structure = get_structure(structure)

    def compute():
        data = samples if samples is not None else sample_set(scenario, **(sample_options or {}))
        problem = FitProblem(structure, data, prior_scale=prior_scale)
        n_screen, n_starts, max_nfev = effort or EFFORT[base_name(structure.name)]
        usable = [c for c in (complete_start(problem, s) for s in starts) if c is not None]
        block = fit_tmsoh_block(structure, data, args) if n_screen > 0 else {}
        fit = fit_structure(problem, n_screen=max(int(n_screen * args.scale), 0), n_starts=n_starts, max_nfev=max_nfev,
                            seed=args.seed, workers=args.workers, start_points=usable, seed_values=block)
        fit['tmsoh_block'] = block
        fit.update(scenario=scenario, tag=tag, label=structure.label,
                   barriers=tabulate_reaction_barriers(structure, fit['theta']))
        return fit

    path = fit_path(structure.name, scenario, tag)
    fresh = not path.exists()
    fit = cached(path, compute)
    log(f"fit {structure.name:10s} {scenario:7s}{tag:24s} chi2 {fit['chi2']:10.2f}  max |z| {fit['max_abs_z']:6.2f}  "
        f"{'FITS' if fit['fits'] else 'fails'}  [{short(fit['theta'])}]"
        + (f"  ({fit['n_evaluations']} evaluations, {fit['seconds']:.0f} s)" if fresh else '  (cached)'))
    return fit


def run_structure(structure: str, scenario: str, args, extra_starts=()) -> dict:
    """Fit of a registered structure, started also from its parents in this scenario and from its own other scenarios."""
    starts = list(extra_starts)
    for parent in PARENTS.get(structure, ()):
        if fit_path(parent, scenario).exists():
            parent_fit = load_fit_result(str(fit_path(parent, scenario)))
            starts.append(embed_parent_theta(structure, parent_fit['theta']))
    for other in SCENARIO_ORDER:
        if other != scenario and fit_path(structure, other).exists():
            starts.append({k: v for k, v in load_fit_result(str(fit_path(structure, other)))['theta'].items()
                           if not k.startswith('log10 age')})
    return run_one(structure, scenario, args, starts=starts)


def load_main_fits() -> dict:
    """{(structure, scenario): fit} of every untagged fit on disk."""
    out = {}
    for structure in SIZE_ORDER:
        for scenario in SCENARIOS:
            if fit_path(structure, scenario).exists():
                out[(structure, scenario)] = load_fit_result(str(fit_path(structure, scenario)))
    return out


def select_consistent(fits: dict):
    """(structure, scenario, kind) by the rule fixed before fitting, or None: the smallest structure that fits a
    common-age scenario; failing that, the smallest that fits with free ages."""
    for kind, scenarios in (('common age', ('short', 'middle', 'long')), ('free ages', ('free',))):
        for structure in SIZE_ORDER:
            fitting = [sc for sc in scenarios if (structure, sc) in fits and fits[(structure, sc)]['fits']]
            if fitting:
                best = min(fitting, key=lambda sc: fits[(structure, sc)]['chi2'])
                return structure, best, kind
    return None


def summarize_fits() -> pd.DataFrame:
    fits = load_main_fits()
    if not fits:
        return pd.DataFrame()
    table = compare_structures({f'{st} | {sc}': fit for (st, sc), fit in fits.items()})
    table.insert(1, 'scenario', [sc for _, sc in fits])
    theta = pd.DataFrame([{'structure': st, 'scenario': sc, **fit['theta']} for (st, sc), fit in fits.items()])
    table.reset_index(drop=True).to_csv(RESULTS / 'fits_summary.csv', index=False)
    theta.to_csv(RESULTS / 'fits_parameters.csv', index=False)
    pd.concat([fit['table'].assign(structure=st, scenario=sc) for (st, sc), fit in fits.items()],
              ignore_index=True).to_csv(RESULTS / 'fits_residuals.csv', index=False)
    pd.DataFrame([{'structure': st, 'scenario': sc, 'reaction': r, **b} for (st, sc), fit in fits.items()
                  for r, b in fit['barriers'].items()]).to_csv(RESULTS / 'fits_barriers.csv', index=False)
    selected = select_consistent(fits)
    log('selection rule: ' + ('no structure fits any scenario' if selected is None else
                              f'{selected[0]} in the {selected[1]} scenario ({selected[2]})'))
    return table


def stage_baselines(args) -> None:
    for scenario in SCENARIO_ORDER:
        for structure in ('M0-BEP', 'M0-Marcus', 'M1'):
            run_structure(structure, scenario, args)
    summarize_fits()
    # Check: the single barriers against the intervals of notebook 03 (one per observation, 24 h assumed age)
    nb03 = pd.read_csv(REPO / 'notebooks' / 'results' / '03' / 'm0_intervals.csv', index_col=0)
    for structure, label in (('M0-BEP', 'M0 · capped BEP, one E0'), ('M0-Marcus', 'M0 · Marcus, one g')):
        rows = nb03.loc[label]
        for scenario in SCENARIOS:
            g = load_fit_result(str(fit_path(structure, scenario)))['theta']['g all reactions']
            inside = [row['observation'] for _, row in rows.iterrows() if not str(row['verdict']).startswith('no value')
                      and (np.isnan(row['g_low_eV']) or g >= row['g_low_eV'])
                      and (np.isnan(row['g_high_eV']) or g <= row['g_high_eV'])]
            log(f"check {structure} {scenario}: best barrier {g:.3f} eV lies in the notebook 03 interval of "
                f"{', '.join(inside) or 'no observation'}")


def stage_extended(args) -> None:
    for scenario in SCENARIO_ORDER:
        for structure in ('M1-split', 'M3', 'M3-split'):
            run_structure(structure, scenario, args)
    summarize_fits()
    fits, ok = load_main_fits(), True
    for scenario in SCENARIOS:
        for child, parents in PARENTS.items():
            for parent in parents:
                if (child, scenario) in fits and (parent, scenario) in fits:
                    c, p = fits[(child, scenario)]['chi2'], fits[(parent, scenario)]['chi2']
                    if c > p + 0.05:
                        ok = False
                        log(f'NESTING VIOLATED in {scenario}: {child} chi2 {c:.2f} > {parent} chi2 {p:.2f}')
    log(f"nesting check: {'chi2 never rises from a structure to the one that extends it' if ok else 'FAILED'}")


def stage_refit(args) -> None:
    """Fits of M3 and M3-split checked against the tightened absolute tolerance of the search, and redone where needed.

    The first pass searched with atol = 1e-10 M. A fit whose χ² at that tolerance differs from its reported χ²
    by more than 0.5 was found on a noisy surface: it is fitted again from scratch, started also from the old
    solution and from the lowest point of its first-pass profile. The others are polished again from the same
    points. An improved fit replaces the old one, which is kept under fits_first_pass/.
    """
    (RESULTS / 'fits_first_pass').mkdir(exist_ok=True)
    for structure in ('M3', MAIN_STRUCTURE):
        for scenario in SCENARIO_ORDER:
            path = fit_path(structure, scenario)
            old = load_fit_result(str(path))
            coarse = FitProblem(structure, sample_set(scenario), solver=FIRST_PASS_SOLVER)
            at_search = coarse.chi2(coarse.vector(old['theta']))
            starts = [old['theta']]
            for profile_path in (result_path('profiles', structure, scenario),
                                 RESULTS / 'first_pass_atol1e-10' / f'profile_{structure}__{scenario}.json'):
                if profile_path.exists():
                    profile = load_fit_result(str(profile_path))['profile']
                    starts.append(get_profile_theta(profile.loc[profile['chi2'].idxmin()]))
            noisy = abs(at_search - old['chi2']) > 0.5
            keep = RESULTS / 'fits_first_pass' / path.name
            if not keep.exists():
                path.rename(keep)
            else:
                path.unlink()
            if noisy:
                new = run_structure(structure, scenario, args, extra_starts=starts)
            else:
                new = run_one(structure, scenario, args, starts=starts, effort=(0, 1, 40))
            if new['chi2'] > old['chi2'] - 0.05:
                write_fit_result(old, str(path))
                log(f'refit {structure} {scenario}: kept (chi2 {old["chi2"]:.2f}; second pass {new["chi2"]:.2f})')
                continue
            log(f'refit {structure} {scenario}: chi2 {old["chi2"]:.2f} -> {new["chi2"]:.2f}'
                + (' (the first pass was tolerance-sensitive: full refit)' if noisy else ' (polished from the profile minimum)'))
    summarize_fits()


# ------------------------------------------------------------------------------
# night: blind spots
# ------------------------------------------------------------------------------
def result_path(kind: str, structure: str, scenario: str, tag: str = '') -> Path:
    return RESULTS / kind / f'{structure}__{scenario}{tag}.json'.replace(' ', '_')


NIGHT_BARRIER_OFFSETS_EV = (-0.30, -0.15, -0.10, -0.05, -0.025, 0.025, 0.05, 0.10, 0.15, 0.30)
NIGHT_ENERGY_OFFSETS_SE = (-2.0, -1.0, -0.5, -0.25, 0.25, 0.5, 1.0, 2.0)


def night_grids(problem: FitProblem, theta: dict, names: list) -> dict:
    """Profile grids of the night run: 10 offsets for a barrier, 8 for a freed energy (the library default has
    14 and 12), the whole box in 13 steps for an age or the water fraction. Interval edges are refined afterwards."""
    sigma = {row['quantity']: row['sigma'] for row in problem.build(problem.vector(theta))['prior_rows']}
    grids = {}
    for p in problem.all_parameters:
        if p.name not in names:
            continue
        if p.kind == 'barrier':
            grid = theta[p.name] + np.array(NIGHT_BARRIER_OFFSETS_EV)
        elif p.kind == 'energy':
            grid = theta[p.name] + next(v for k, v in sigma.items() if p.target in k) * np.array(NIGHT_ENERGY_OFFSETS_SE)
        else:
            grid = np.linspace(p.lower, p.upper, 13)
        grid = grid[(grid >= p.lower - 1e-12) & (grid <= p.upper + 1e-12)]
        grids[p.name] = np.unique(np.round(grid[np.abs(grid - theta[p.name]) > 1e-9], 9))
    return grids


def select_alternatives(optima: list, theta: dict, max_n: int = 1, max_delta: float = 10.0) -> list:
    """Other local optima worth following in a profile: within max_delta of the best χ² and a different solution
    (some barrier or energy more than 0.03 eV from the best fit and from each other)."""
    def differs(a, b):
        return any(abs(a[k] - b[k]) > 0.03 for k in a if k.startswith(('g ', 'dG ')))

    lowest, chosen = min(o['chi2'] for o in optima), []
    for optimum in sorted(optima, key=lambda o: o['chi2']):
        if optimum['chi2'] - lowest > max_delta or len(chosen) == max_n:
            break
        if all(differs(optimum['theta'], other) for other in [theta] + chosen):
            chosen.append(optimum['theta'])
    return chosen


def run_profiles(structure: str, scenario: str, args, *, refine: bool, names=None, tag: str = '', solver: dict = None,
                 fit_tag: str = '', follow_optima: bool = True) -> dict:
    """Profiles of one fit: a grid walk from the best fit, the same from its nearest other local optimum, the
    accepted local optima as extra points, then more points where the intervals end.

    The recovery check showed why the other optima matter: a chain that only follows the best fit can stay in
    its basin and end an interval too early.
    """
    def compute():
        fit = load_fit_result(str(fit_path(structure, scenario, fit_tag)))
        problem = FitProblem(structure, sample_set(scenario), solver=solver)
        which = list(names) if names is not None else list(problem.names)
        profile = calculate_profile(problem, fit, which, grids=night_grids(problem, fit['theta'], which),
                                    max_nfev=args.profile_nfev, workers=args.workers)
        n_followed = 0
        if follow_optima and fit['fits']:
            # χ² of the stored optima again, at this problem's tolerance
            optima = [{'theta': o['theta'], 'chi2': problem.chi2(problem.vector(o['theta']))} for o in fit['local_optima'][:12]]
            alternatives = select_alternatives(optima, fit['theta'])
            n_followed = len(alternatives)
            if alternatives:
                profile = extend_profile(problem, profile, alternatives, which, max_nfev=args.profile_nfev, workers=args.workers)
            profile = add_profile_points(profile, optima, which)
        if refine:
            profile = refine_profile_edges(problem, profile, which, max_nfev=args.profile_nfev, workers=args.workers)
        intervals = pd.DataFrame([find_confidence_interval(profile, name) for name in which])
        directions = calculate_parameter_directions(problem.replace(solver=REPORT_SOLVER), fit['theta'])
        return {'structure': structure, 'scenario': scenario, 'profile': profile, 'intervals': intervals,
                'directions': directions['directions'], 'sigma_local': directions['sigma_local'].to_dict(),
                'other_optima_followed': n_followed, 'evaluations': int(profile['nfev'].sum())}

    path = result_path('profiles', structure, scenario, tag)
    fresh = not path.exists()
    t0 = time.perf_counter()
    result = cached(path, compute)
    log(f'profiles {structure} {scenario}{tag}: ' + '; '.join(
        f"{r['parameter'].replace('g ', '')} " + (f"{r['low']:.3f}–{r['high']:.3f}" if r['status'] == 'interval' else
                                                   f"≥ {r['low']:.3f}" if r['status'] == 'lower bound only' else
                                                   f"≤ {r['high']:.3f}" if r['status'] == 'upper bound only' else 'not determined')
        for _, r in result['intervals'].iterrows())
        + (f"  ({result['evaluations']} evaluations, {time.perf_counter() - t0:.0f} s)" if fresh else '  (cached)'))
    return result


def adopt_profile_minimum(structure: str, scenario: str, args, *, refine: bool) -> bool:
    """A profile can find a lower χ² than the fit it started from. The fit is then polished from that point; if
    it improves by more than 0.05, it replaces the fit (the old one is kept under fits_before_profiles/), the
    profile is walked again from the new best fit where the structure fits, and Δχ² is measured from the new χ².
    Returns True when the fit was replaced."""
    path, profile_path = fit_path(structure, scenario), result_path('profiles', structure, scenario)
    fit, bundle = load_fit_result(str(path)), load_fit_result(str(profile_path))
    profile = bundle['profile']
    if profile['chi2'].min() > fit['chi2'] - 0.05:
        return False
    lowest = get_profile_theta(profile.loc[profile['chi2'].idxmin()])
    new = run_one(structure, scenario, args, tag='__from_profile_minimum', starts=[lowest, fit['theta']], effort=(0, 2, 60))
    if new['chi2'] > fit['chi2'] - 0.05:
        log(f"profile minimum {structure} {scenario}: fit kept (chi2 {fit['chi2']:.2f}; polished from the profile minimum {new['chi2']:.2f})")
        return False
    (RESULTS / 'fits_before_profiles').mkdir(exist_ok=True)
    path.rename(RESULTS / 'fits_before_profiles' / path.name)
    write_fit_result({**new, 'tag': '', 'local_optima': new['local_optima'] + fit['local_optima']}, str(path))
    problem = FitProblem(structure, sample_set(scenario))
    names = list(dict.fromkeys(profile['parameter']))
    reference = problem.chi2(problem.vector(new['theta']))
    if new['fits']:
        profile = extend_profile(problem, profile, [new['theta']], names, max_nfev=args.profile_nfev, workers=args.workers)
    profile = add_profile_points(profile.assign(delta=profile['chi2'] - reference), [{'theta': new['theta'], 'chi2': reference}], names)
    if refine and new['fits']:
        profile = refine_profile_edges(problem, profile, names, max_nfev=args.profile_nfev, workers=args.workers)
    profile['delta'] = profile['chi2'] - reference
    directions = calculate_parameter_directions(problem.replace(solver=REPORT_SOLVER), new['theta'])
    bundle.update(profile=profile, intervals=pd.DataFrame([find_confidence_interval(profile, name) for name in names]),
                  directions=directions['directions'], sigma_local=directions['sigma_local'].to_dict(),
                  evaluations=int(profile['nfev'].sum()), fit_replaced_from_profile_minimum=True)
    write_fit_result(bundle, str(profile_path))
    log(f"profile minimum {structure} {scenario}: fit replaced, chi2 {fit['chi2']:.2f} -> {new['chi2']:.2f}, "
        f"max |z| {new['max_abs_z']:.2f}, {'FITS' if new['fits'] else 'fails'}  [{short(new['theta'])}]")
    return True


FITTED_SAMPLES = ('0.5 % H2O', '2 % H2O', 'TMSPa + TMSOH (A)', 'TMSOH, probe', 'TMSOH, glovebox', 'E1')


def run_leave_one_out(structure: str, scenario: str, args) -> dict:
    def compute():
        fit = load_fit_result(str(fit_path(structure, scenario)))
        sets = {name: sample_set(scenario, exclude=(name,)) for name in FITTED_SAMPLES}
        results = fit_leave_one_out(structure, sets, fit, n_screen=max(int(128 * args.scale), 0), n_starts=2, max_nfev=30,
                                    seed=args.seed, workers=args.workers)
        rows = []
        for name, r in results.items():
            barriers = tabulate_reaction_barriers(get_structure(structure), {**fit['theta'], **r['theta']})
            rows.append({'left out': name, 'chi2': r['chi2'], 'max |z|': r['max_abs_z'], 'worst': r['worst'],
                         'fits': r['fits'], 'evaluations': r['n_evaluations'], **r['theta'],
                         **{f'dG‡ {rxn}': b['dG_barrier_eV'] for rxn, b in barriers.items()}})
        return {'structure': structure, 'scenario': scenario, 'table': pd.DataFrame(rows),
                'evaluations': int(sum(r['n_evaluations'] for r in results.values()))}

    path = result_path('loo', structure, scenario)
    fresh = not path.exists()
    t0 = time.perf_counter()
    result = cached(path, compute)
    log(f"leave-one-out {structure} {scenario}: " + '; '.join(
        f"without {r['left out']}: chi2 {r['chi2']:.1f}, max |z| {r['max |z|']:.1f}" for _, r in result['table'].iterrows())
        + (f"  ({result['evaluations']} evaluations, {time.perf_counter() - t0:.0f} s)" if fresh else '  (cached)'))
    return result


SENSITIVITIES = {
    'floor x0': ({'floor_scale': 0.0}, 1.0),
    'floor x2': ({'floor_scale': 2.0}, 1.0),
    'covariance x0.5': ({}, 0.5),
    'covariance x2': ({}, 2.0),
    'ramp held 8 h': ({'overrides': {'2 % H2O': {'heated_duration_h': 8.0}}}, 1.0),
    'glovebox 4 h': ({'overrides': {'TMSOH, glovebox': {'hold_h': 4.0}}}, 1.0),
    'glovebox 16 h': ({'overrides': {'TMSOH, glovebox': {'hold_h': 16.0}}}, 1.0),
}


def run_sensitivities(structure: str, scenario: str, args) -> pd.DataFrame:
    main = load_fit_result(str(fit_path(structure, scenario)))
    rows = [{'variant': 'as fitted', 'scenario': scenario, 'search': 'full', 'chi2': main['chi2'], 'max |z|': main['max_abs_z'],
             'worst': main['worst'], 'fits': main['fits'], **main['theta'],
             **{f'dG‡ {rxn}': b['dG_barrier_eV'] for rxn, b in main['barriers'].items()}}]
    for variant, (options, prior_scale) in SENSITIVITIES.items():
        tag = f'__{variant}'.replace(' ', '_')
        fit = run_one(structure, scenario, args, tag=tag, starts=[main['theta']],
                      sample_options=options, prior_scale=prior_scale, effort=(128, 3, 30))
        search = 'reduced'
        if not fit['fits']:
            # A variant is not called a failure on a reduced search: it is fitted again with the full one
            full = run_one(structure, scenario, args, tag=tag + '__full_search', starts=[main['theta'], fit['theta']],
                           sample_options=options, prior_scale=prior_scale)
            fit, search = (full, 'full') if full['chi2'] < fit['chi2'] else (fit, 'full, no better')
        rows.append({'variant': variant, 'scenario': scenario, 'search': search, 'chi2': fit['chi2'],
                     'max |z|': fit['max_abs_z'], 'worst': fit['worst'], 'fits': fit['fits'], **fit['theta'],
                     **{f'dG‡ {rxn}': b['dG_barrier_eV'] for rxn, b in fit['barriers'].items()}})
    return pd.DataFrame(rows)


def run_indistinguishability(scenario: str, args) -> list:
    """Can the data separate (a) equilibrium from exhausted water, (b) Marcus from Agmon–Levine, (c) one hydrolysis
    barrier from two? Each row: one structure fitted to the lab shares, or to noise-free shares made by another."""
    main = load_fit_result(str(fit_path(MAIN_STRUCTURE, scenario)))
    lab = sample_set(scenario)
    rows = []

    def record(question, fitted, data, fit, reference_chi2=np.nan):
        rows.append({'scenario': scenario, 'question': question, 'fitted': fitted, 'data': data, 'chi2': fit['chi2'],
                     'max |z|': fit['max_abs_z'], 'worst': fit['worst'], 'within 3 SE': fit['fits'],
                     'chi2 of the reference fit': reference_chi2,
                     'water fraction': fit['theta'].get('water fraction', np.nan)})

    # (a) equilibrium or exhausted water
    water = run_one('M3-split+W', scenario, args, starts=[main['theta']], effort=(256, 4, 40))
    record('equilibrium or exhausted water', 'M3-split+W', 'lab', water, main['chi2'])
    split = load_fit_result(str(fit_path('M1-split', scenario)))
    exhausted = run_one('M1-split+W', scenario, args, starts=[split['theta']], effort=(256, 6, 40))
    record('equilibrium or exhausted water', 'M1-split+W', 'lab', exhausted, main['chi2'])
    made_by_main = simulate_synthetic_shares(MAIN_STRUCTURE, main['theta'], lab)
    mimic = run_one('M1-split+W', scenario, args, tag='__on_M3-split_shares', starts=[exhausted['theta']],
                    effort=(256, 6, 40), samples=made_by_main)
    record('equilibrium or exhausted water', 'M1-split+W', 'made by M3-split', mimic, 0.0)
    made_by_water = simulate_synthetic_shares('M1-split+W', exhausted['theta'], lab)
    mimic = run_one(MAIN_STRUCTURE, scenario, args, tag='__on_M1-split+W_shares', starts=[main['theta']],
                    effort=(256, 6, 40), samples=made_by_water)
    record('equilibrium or exhausted water', MAIN_STRUCTURE, 'made by M1-split+W', mimic, 0.0)
    # (b) barrier shape
    shaped = get_structure(MAIN_STRUCTURE).with_shape('agmon_levine')
    other_shape = run_one(shaped, scenario, args, starts=[main['theta']], effort=(0, 1, 40))
    record('Marcus or Agmon–Levine', shaped.name, 'lab', other_shape, main['chi2'])
    # (c) one hydrolysis barrier or two
    mimic = run_one('M3', scenario, args, tag='__on_M3-split_shares',
                    starts=[load_fit_result(str(fit_path('M3', scenario)))['theta']], effort=(256, 6, 40), samples=made_by_main)
    record('one hydrolysis barrier or two', 'M3', 'made by M3-split', mimic, 0.0)
    record('one hydrolysis barrier or two', 'M3', 'lab', load_fit_result(str(fit_path('M3', scenario))), main['chi2'])
    return rows


def run_water_profile(scenario: str, args) -> dict:
    structure = 'M3-split+W'

    def compute():
        fit = load_fit_result(str(fit_path(structure, scenario)))
        problem = FitProblem(structure, sample_set(scenario))
        profile = calculate_profile(problem, fit, ['water fraction'], max_nfev=args.profile_nfev, workers=args.workers)
        return {'structure': structure, 'scenario': scenario, 'profile': profile,
                'intervals': pd.DataFrame([find_confidence_interval(profile, 'water fraction')]),
                'evaluations': int(profile['nfev'].sum())}

    result = cached(result_path('profiles', structure, scenario), compute)
    r = result['intervals'].iloc[0]
    log(f"water fraction {scenario}: best {r['best']:.2f}, 95 % range {r['low']:.2f} to {r['high']:.2f} ({r['status']})")
    return result


def stage_night(args) -> None:
    # 0. The recovery check once more, with the code as it now stands (the first run is kept in recovery_first_run/)
    if not stage_recovery(args):
        raise SystemExit('the synthetic recovery failed: the night run stops here')
    fits = load_main_fits()
    selected = select_consistent(fits)
    profiled = [MAIN_STRUCTURE, 'M1'] + ([selected[0]] if selected and selected[0] not in (MAIN_STRUCTURE, 'M1') else [])
    log(f'night: profiles of {profiled}; leave-one-out, sensitivities and indistinguishability on {MAIN_STRUCTURE}')
    # 1. Profiles of the main structure, then of the baseline structure
    for structure in profiled:
        for scenario in SCENARIO_ORDER:
            run_profiles(structure, scenario, args, refine=scenario in ('middle', 'free'))
    replaced = [adopt_profile_minimum(structure, scenario, args, refine=scenario in ('middle', 'free'))
                for structure in profiled for scenario in SCENARIO_ORDER]
    if any(replaced):
        summarize_fits()
        fits = load_main_fits()
    # 2. Which sample drives what, in the scenarios the structure fits
    for scenario in SCENARIO_ORDER:
        if fits[(MAIN_STRUCTURE, scenario)]['fits']:
            run_leave_one_out(MAIN_STRUCTURE, scenario, args)
        else:
            log(f'leave-one-out {MAIN_STRUCTURE} {scenario}: skipped, the structure does not fit this scenario')
    # 3. Structures the data cannot tell apart
    rows = []
    for scenario in ('middle', 'free'):
        rows += run_indistinguishability(scenario, args)
        run_water_profile(scenario, args)
    pd.DataFrame(rows).to_csv(RESULTS / 'indistinguishability.csv', index=False)
    # 4. Sensitivity to the declared error floor, the covariance and the times the files do not fix
    pd.concat([run_sensitivities(MAIN_STRUCTURE, scenario, args) for scenario in ('middle', 'free')],
              ignore_index=True).to_csv(RESULTS / 'sensitivities.csv', index=False)
    # 5. Solver tolerance: two profiles repeated at rtol 1e-8
    names = ['g solvent_attack', 'g hydrolysis_R23']
    loose = run_profiles(MAIN_STRUCTURE, 'middle', args, refine=True)['intervals'].set_index('parameter')
    tight = run_profiles(MAIN_STRUCTURE, 'middle', args, refine=True, names=names, tag='__rtol1e-8',
                         solver=REPORT_SOLVER)['intervals'].set_index('parameter')
    check = pd.DataFrame({'low, rtol 1e-6': loose.loc[names, 'low'], 'low, rtol 1e-8': tight.loc[names, 'low'],
                          'high, rtol 1e-6': loose.loc[names, 'high'], 'high, rtol 1e-8': tight.loc[names, 'high'],
                          'status, rtol 1e-6': loose.loc[names, 'status'], 'status, rtol 1e-8': tight.loc[names, 'status']})
    check.to_csv(RESULTS / 'tolerance_check.csv')
    worst = np.nanmax(np.abs(np.concatenate([check['low, rtol 1e-6'] - check['low, rtol 1e-8'],
                                             check['high, rtol 1e-6'] - check['high, rtol 1e-8'], [0.0]])))
    same = bool((check['status, rtol 1e-6'] == check['status, rtol 1e-8']).all())
    log(f"tolerance check: interval edges move by at most {worst:.4f} eV from rtol 1e-6 to 1e-8; same status: {same}")
    summarize_night()


def summarize_night() -> None:
    """Flat tables of the night's results for the notebook. Intervals are read off the stored profiles here, so
    every table uses the same interval rule."""
    intervals, bands, loo = [], [], []
    for path in sorted((RESULTS / 'profiles').glob('*.json')):
        result = load_fit_result(str(path))
        tag = path.stem.split('__')[2] if path.stem.count('__') > 1 else ''
        names = list(dict.fromkeys(result['profile']['parameter']))
        current = pd.DataFrame([find_confidence_interval(result['profile'], name) for name in names])
        intervals.append(current.assign(structure=result['structure'], scenario=result['scenario'], tag=tag))
        if tag or '+W' in result['structure']:
            continue
        structure, accepted = get_structure(result['structure']), result['profile'][result['profile']['delta'] <= DELTA_CHI2_95]
        barriers = pd.DataFrame([{rxn: b['dG_barrier_eV'] for rxn, b in tabulate_reaction_barriers(structure, get_profile_theta(row)).items()}
                                 for _, row in accepted.iterrows()])
        energies = pd.DataFrame([{rxn: b['dG_rxn_eV'] for rxn, b in tabulate_reaction_barriers(structure, get_profile_theta(row)).items()}
                                 for _, row in accepted.iterrows()])
        for rxn in barriers:
            bands.append({'structure': result['structure'], 'scenario': result['scenario'], 'reaction': rxn,
                          'dG‡ low': barriers[rxn].min(), 'dG‡ high': barriers[rxn].max(),
                          'dG_rxn low': energies[rxn].min(), 'dG_rxn high': energies[rxn].max(),
                          'accepted sets': len(accepted)})
    for path in sorted((RESULTS / 'loo').glob('*.json')):
        result = load_fit_result(str(path))
        loo.append(result['table'].assign(structure=result['structure'], scenario=result['scenario']))
    if intervals:
        pd.concat(intervals, ignore_index=True).to_csv(RESULTS / 'intervals.csv', index=False)
        pd.DataFrame(bands).to_csv(RESULTS / 'barrier_bands.csv', index=False)
    if loo:
        pd.concat(loo, ignore_index=True).to_csv(RESULTS / 'leave_one_out.csv', index=False)


# ------------------------------------------------------------------------------
# traces: does a fit survive a trace of TMSOH at mixing?
# ------------------------------------------------------------------------------
WATER_SAMPLES = ('0.5 % H2O', '2 % H2O')                      # the fitted samples whose recipe holds no TMSOH
TRACE_SCAN_M = (0.0, 1e-9, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3)
TRACE_REFIT_M = {'1uM': 1e-6, '1mM': 1e-3}


def trace_options(level_M: float) -> dict:
    return {'overrides': {name: {'added_M': {'TMSOH': level_M}} for name in WATER_SAMPLES}}


def stage_traces(args) -> None:
    """Every fit starts the two water samples from the recipe alone, with no TMSOH at mixing. A TMSPA stock holds
    some hydrolysis product, so each fit of the main structure is (a) evaluated unchanged with a trace of TMSOH
    added at mixing and (b) fitted again with 1 µM and with 1 mM of it. A fit that (a) destroys rests on an
    induction time seeded by the reaction itself; (b) says whether the scenario can still be fitted."""
    rows = []

    def record(scenario, level, how, problem, theta, **extra):
        x = problem.vector(theta)
        table = problem.table(x)
        worst = table.loc[table['z'].abs().idxmax()]
        barriers = tabulate_reaction_barriers(get_structure(MAIN_STRUCTURE), theta)
        rows.append({'scenario': scenario, 'TMSOH at mixing (M)': level, 'parameters': how, 'chi2': problem.chi2(x),
                     'max |z|': float(abs(worst['z'])), 'worst': f"{worst['sample']}: {worst['quantity']}",
                     'fits': bool(abs(worst['z']) <= 3.0), **extra, **theta,
                     **{f'dG‡ {rxn}': b['dG_barrier_eV'] for rxn, b in barriers.items()}})
        pd.DataFrame(rows).to_csv(RESULTS / 'trace_sensitivity.csv', index=False)

    fits = {scenario: load_fit_result(str(fit_path(MAIN_STRUCTURE, scenario)))
            for scenario in SCENARIO_ORDER if fit_path(MAIN_STRUCTURE, scenario).exists()}
    for scenario, fit in fits.items():
        for level in TRACE_SCAN_M:
            problem = FitProblem(MAIN_STRUCTURE, sample_set(scenario, **trace_options(level)), solver=REPORT_SOLVER)
            record(scenario, level, 'as fitted without it', problem, fit['theta'])
        scan = [r for r in rows if r['scenario'] == scenario]
        log(f'trace scan {MAIN_STRUCTURE} {scenario}: chi2 ' + ', '.join(
            f"{r['chi2']:.1f} at {r['TMSOH at mixing (M)']:.0e} M" for r in scan))
    starts = [{k: v for k, v in fit['theta'].items() if not k.startswith('log10 age')} for fit in fits.values()]
    for label, level in TRACE_REFIT_M.items():
        for scenario, fit in fits.items():
            if not fit['fits']:
                continue
            new = run_one(MAIN_STRUCTURE, scenario, args, tag=f'__TMSOH_{label}', starts=[fit['theta']] + starts,
                          sample_options=trace_options(level))
            problem = FitProblem(MAIN_STRUCTURE, sample_set(scenario, **trace_options(level)), solver=REPORT_SOLVER)
            record(scenario, level, 'fitted again', problem, new['theta'], evaluations=new['n_evaluations'])


# ------------------------------------------------------------------------------
# predictions: hold-out and checks, storage at 25 °C, Peter's protocol with its ΔS‡ band
# ------------------------------------------------------------------------------
STORAGE_C0_M = {'TMSPA': 0.050, 'H2O': 0.020, 'EC': 4.5}      # the sealed electrolyte of notebook 01, Block 8
STORAGE_T_K = 298.15
DS_RANGE_J_MOL_K = (-150.0, -100.0, -50.0, 0.0, 50.0)         # declared range of the activation entropy
T_OBSERVED_K = {'hydrolysis': 295.65, 'transfer': 295.65, 'condensation': 353.15, 'solvent_attack': 353.15}
EC_DEC_DENSITY_G_ML = 0.5 * (1.321 + 0.975)                   # as in notebook 03 §5.4


def calculate_storage_times(model, tmsoh_M: float = 0.0) -> dict:
    """Hours to consume 90 % of the water, and TMSPA left after one year, in the storage case at 25 °C.
    tmsoh_M: TMSOH present at mixing on top of the recipe."""
    times = np.geomspace(1.0, 100 * 365.25 * 86400.0, 400)
    c0 = {**STORAGE_C0_M, 'TMSOH': tmsoh_M} if tmsoh_M > 0.0 else STORAGE_C0_M
    sim = simulate_history(c0, [(STORAGE_T_K, times[-1])], times, model=model, **REPORT_SOLVER)
    if not sim['success']:
        return {'t90_water_h': np.nan, 'tmspa_left_1y': np.nan}
    t_h, C, idx = times / 3600.0, sim['C_M'], sim['idx']
    year = int(np.argmin(np.abs(times - 365.25 * 86400.0)))
    return {'t90_water_h': find_crossing_time(t_h, C[idx['H2O']], 0.1 * STORAGE_C0_M['H2O']),
            'tmspa_left_1y': float(C[idx['TMSPA'], year] / STORAGE_C0_M['TMSPA'])}


def with_activation_entropy(model, dS_J_mol_K: float):
    """The model with ΔS‡ = dS for every class, each barrier anchored at the temperature where it was observed."""
    for family in list(model.family_params):
        if family == 'default':
            continue
        T_ref = next(T for key, T in T_OBSERVED_K.items() if family.startswith(key))
        model = model.with_family_params(family, dS_act_J_mol_K=float(dS_J_mol_K), T_ref_K=T_ref, name=model.name)
    return model


def accepted_thetas(structure: str, scenario: str) -> list:
    """The best fit and every profile point within Δχ² = 3.84 of it."""
    thetas = [load_fit_result(str(fit_path(structure, scenario)))['theta']]
    path = result_path('profiles', structure, scenario)
    if path.exists():
        profile = load_fit_result(str(path))['profile']
        thetas += [get_profile_theta(row) for _, row in profile[profile['delta'] <= DELTA_CHI2_95].iterrows()]
    return thetas


def check_hold_out(structure: str, scenario: str, theta: dict) -> pd.DataFrame:
    built = get_structure(structure).build(theta)
    samples = sample_set(scenario, roles=('hold_out',), paper=())
    return tabulate_residuals(built['model'], samples, ages=built['ages'], water_fraction=built['water_fraction'])


CO_PRODUCTS = ('none', 'TMSOH', 'HMDSO')


def check_tmspa_alone(structure: str, theta: dict, co_products: str = 'none') -> dict:
    """TMSPa without added water: the water content with which the model turns the first spectrum into the second.

    The first spectrum already shows hydrolysed phosphate, so silyl groups had left it by then. co_products says
    where they are at that time: 'none' (not in the solution), 'TMSOH' or 'HMDSO'. The files do not tell.
    """
    from scipy.optimize import minimize_scalar
    from kinetics.fitting.residuals import _shares_from_state
    from kinetics.data.observables import DEFAULT_ERROR_FLOOR
    record = load_lab_observables()['samples']['TMSPa alone']
    first, late = record['spectra'][0], record['spectra'][-1]
    windows = list(first['shares'])
    start = np.clip([first['shares'][w]['share'] for w in windows], 0.0, None)
    start = start / start.sum()
    total = record['c0_M']['TMSPA']
    measured = np.array([late['shares'][w]['share'] for w in windows])
    floor = DEFAULT_ERROR_FLOOR['31P']
    sigma = np.sqrt(np.array([late['shares'][w]['sigma_noise'] ** 2 + late['shares'][w]['sigma_baseline'] ** 2 for w in windows])
                    + np.maximum(floor['absolute'], floor['relative'] * np.maximum(measured, 0.0)) ** 2)
    dt_s = 3600.0 * (late['t_since_first_h'] - first['t_since_first_h'])
    model = get_structure(structure).build(theta)['model']
    lost = {'TMSPA': 0, 'BMSPA': 1, 'MMSPA': 2, 'H3PO4': 3}            # silyl groups off the phosphate, per species
    released = total * float(sum(share * lost[w] for w, share in zip(windows, start)))
    extra = {'none': {}, 'TMSOH': {'TMSOH': released}, 'HMDSO': {'HMDSO': 0.5 * released}}[co_products]

    def predict(water_M):
        c0 = {**{w: total * share for w, share in zip(windows, start)}, **extra, 'H2O': water_M, 'EC': record['c0_M']['EC']}
        sim = simulate_history(c0, [(record['T_rt_C'] + 273.15, dt_s)], [dt_s], model=model, **REPORT_SOLVER)
        return _shares_from_state(sim['C_M'], sim['idx'], '31P')[0] if sim['success'] else np.full(4, np.nan)

    def chi2(water_M):
        return float(np.nansum(((predict(water_M) - measured) / sigma) ** 2))

    grid = np.linspace(0.0, 0.5, 26)
    values = [chi2(w) for w in grid]
    j = int(np.argmin(values))
    res = minimize_scalar(chi2, bounds=(grid[max(j - 1, 0)], grid[min(j + 1, 25)]), method='bounded', options={'xatol': 1e-4})
    water = float(res.x if res.fun < values[j] else grid[j])
    predicted = predict(water)
    return {'water_mM': 1000.0 * water, 'water_ppm': water * 18.015 / (EC_DEC_DENSITY_G_ML * 1000.0) * 1e6,
            'chi2': chi2(water), 'max |z|': float(np.nanmax(np.abs((predicted - measured) / sigma))),
            **{f'{w} predicted': p for w, p in zip(windows, predicted)}, **{f'{w} measured': m for w, m in zip(windows, measured)}}


def check_paper_observations(structure: str, theta: dict) -> list:
    """Paper observations no fit uses: E4 (TMSPA + TMSOH) and the 1 and 5 vol% samples of the water series."""
    from kinetics.reactor import calculate_phosphate_fractions, calculate_water_series_c0
    experimental = load_experimental_data()
    model = get_structure(structure).build(theta)['model']
    rows = []
    e4 = next(e for e in experimental['control_experiments'] if e['id'] == 'E4')
    sim = simulate_history(e4['c0_M'], [(e4['T_K'], e4['t_s'])], [e4['t_s']], model=model, **REPORT_SOLVER)
    conversion = 1.0 - sim['C_M'][sim['idx']['TMSPA'], 0] / e4['c0_M']['TMSPA']
    rows.append({'observation': 'E4', 'age': '24 h (assumed)', 'statement': 'TMSPA conversion ≥ 0.95', 'predicted': f'{conversion:.3f}',
                 'agrees': bool(conversion >= 0.95)})
    series = experimental['water_series']
    for sample in series['samples']:
        if sample['id'] not in ('W1', 'W5'):
            continue
        c0 = calculate_water_series_c0(sample['h2o_vol_pct'], series['medium'])
        ages_s = 3600.0 * np.array([1.0, 24.0, 168.0])
        sim = simulate_history(c0, [(series['medium']['T_K'], ages_s[-1])], ages_s, model=model, **REPORT_SOLVER)
        fractions = calculate_phosphate_fractions(np.maximum(sim['C_M'], 0.0), sim['idx'])
        for k, age in enumerate(('1 h', '1 d', '7 d')):
            ok = all((low is None or sum(fractions[sp][k] for sp in key.split('+')) >= low)
                     and (high is None or sum(fractions[sp][k] for sp in key.split('+')) <= high)
                     for key, (low, high) in sample['windows'].items())
            rows.append({'observation': sample['id'], 'age': age, 'statement': sample['observation'],
                         'predicted': ', '.join(f"{sp} {fractions[sp][k]:.2f}" for sp in ('TMSPA', 'BMSPA', 'MMSPA', 'H3PO4')),
                         'agrees': bool(ok)})
    return rows


def tabulate_sample_trajectories(structure: str, scenario: str, fit: dict, trace_M: float = 1e-6) -> pd.DataFrame:
    """Phosphate shares of the two water samples against the time since mixing, at room temperature, for one fit:
    from the recipe alone and with trace_M of TMSOH at mixing. Rows of kind 'measured' place the measured shares
    at the age this fit gives the sample."""
    from kinetics.fitting.residuals import SHARE_SPECIES, _shares_from_state
    model = get_structure(structure).build(fit['theta'])['model']
    times = 3600.0 * np.geomspace(0.1, 2000.0, 240)
    windows, rows = list(SHARE_SPECIES['31P']), []
    for sample in sample_set(scenario):
        if sample['name'] not in WATER_SAMPLES:
            continue
        for level in (0.0, trace_M):
            c0 = {**sample['c0_M'], 'TMSOH': sample['c0_M'].get('TMSOH', 0.0) + level}
            sim = simulate_history(c0, [(sample['T_rt_K'], times[-1])], times, model=model, **REPORT_SOLVER)
            if not sim['success']:
                continue
            shares = _shares_from_state(sim['C_M'], sim['idx'], '31P')
            rows += [{'structure': structure, 'scenario': scenario, 'sample': sample['name'], 'kind': 'model',
                      'TMSOH at mixing (M)': level, 't_h': t / 3600.0, **dict(zip(windows, y))} for t, y in zip(times, shares)]
        age_h = fit['ages'][sample['name']]
        rows += [{'structure': structure, 'scenario': scenario, 'sample': sample['name'], 'kind': 'measured',
                  'TMSOH at mixing (M)': np.nan, 't_h': age_h + t / 3600.0, **dict(zip(windows, y))}
                 for t, y in zip(sample['t_spec_s'], sample['measured'])]
    return pd.DataFrame(rows)


def stage_predictions(args) -> None:
    fits = load_main_fits()
    selected = select_consistent(fits)
    structures = [MAIN_STRUCTURE, 'M1'] + ([selected[0]] if selected and selected[0] not in (MAIN_STRUCTURE, 'M1') else [])
    hold_out, alone, paper, storage, protocol, traces, trajectories = [], [], [], [], [], [], []
    stages, _ = build_protocol_schedule()
    recipe = calculate_recipe_molarities()
    for structure in structures:
        for scenario in SCENARIOS:
            if (structure, scenario) not in fits:
                continue
            theta = fits[(structure, scenario)]['theta']
            built = get_structure(structure).build(theta)
            hold_out.append(check_hold_out(structure, scenario, theta).assign(structure=structure, scenario=scenario))
            checked = {'the fit': theta}
            for label in TRACE_REFIT_M if structure == MAIN_STRUCTURE else ():
                if fit_path(structure, scenario, f'__TMSOH_{label}').exists():
                    checked[f'fitted with TMSOH {label} at mixing'] = load_fit_result(
                        str(fit_path(structure, scenario, f'__TMSOH_{label}')))['theta']
            alone += [{'structure': structure, 'scenario': scenario, 'parameters': which, 'co-products': co,
                       **check_tmspa_alone(structure, parameters, co)}
                      for which, parameters in checked.items() for co in CO_PRODUCTS]
            paper += [{'structure': structure, 'scenario': scenario, **row} for row in check_paper_observations(structure, theta)]
            # Storage at 25 °C: the best fit and the band over the accepted parameter sets
            thetas = accepted_thetas(structure, scenario)
            best = calculate_storage_times(built['model'])
            band = calculate_prediction_band(
                lambda th: calculate_storage_times(get_structure(structure).build(th)['model'])['t90_water_h'], thetas)
            storage.append({'structure': structure, 'scenario': scenario, 't90 water (h), best fit': best['t90_water_h'],
                            'band low (h)': band['low'], 'band high (h)': band['high'], 'accepted sets': len(thetas),
                            'sets that never reach 90 % in 100 years': int(np.isnan(band['values']).sum()),
                            'TMSPA left after 1 year, best fit': best['tmspa_left_1y']})
            log(f"storage {structure} {scenario}: 90 % of the water gone after {best['t90_water_h']:.3g} h "
                f"(accepted sets: {band['low']:.3g} to {band['high']:.3g} h, n = {len(thetas)})")
            if structure == MAIN_STRUCTURE:
                trajectories.append(tabulate_sample_trajectories(structure, scenario, fits[(structure, scenario)]))
            # The same prediction with a trace of TMSOH at mixing: the fit as it is, and the fit made with that trace
            for label, level in TRACE_REFIT_M.items() if structure == MAIN_STRUCTURE else ():
                seeded = calculate_storage_times(built['model'], level)
                row = {'structure': structure, 'scenario': scenario, 'TMSOH at mixing (M)': level,
                       't90 water (h), no trace': best['t90_water_h'], 't90 water (h), fit as it is': seeded['t90_water_h']}
                if fit_path(structure, scenario, f'__TMSOH_{label}').exists():
                    refit = load_fit_result(str(fit_path(structure, scenario, f'__TMSOH_{label}')))
                    again = calculate_storage_times(get_structure(structure).build(refit['theta'])['model'], level)
                    row.update({'t90 water (h), fitted with the trace': again['t90_water_h'],
                                'chi2 fitted with the trace': refit['chi2'], 'fits with the trace': refit['fits']})
                traces.append(row)
            # Peter's protocol: TMSPA left at each acquisition over the declared range of ΔS‡
            for dS in DS_RANGE_J_MOL_K:
                sim = simulate_protocol(stages, recipe, model=with_activation_entropy(built['model'], dS))
                left = calculate_remaining_fraction(sim)
                protocol += [{'structure': structure, 'scenario': scenario, 'dS (J/mol/K)': dS, 'hold T (°C)': T,
                              'TMSPA left': float(v)} for T, v in left.items()]
    pd.concat(hold_out, ignore_index=True).to_csv(RESULTS / 'hold_out.csv', index=False)
    pd.DataFrame(alone).to_csv(RESULTS / 'check_tmspa_alone.csv', index=False)
    pd.DataFrame(paper).to_csv(RESULTS / 'check_paper.csv', index=False)
    pd.DataFrame(storage).to_csv(RESULTS / 'prediction_storage.csv', index=False)
    pd.DataFrame(protocol).to_csv(RESULTS / 'prediction_protocol.csv', index=False)
    pd.DataFrame(traces).to_csv(RESULTS / 'prediction_storage_traces.csv', index=False)
    pd.concat(trajectories, ignore_index=True).to_csv(RESULTS / 'sample_trajectories.csv', index=False)
    worst = pd.concat(hold_out).groupby(['structure', 'scenario'])['z'].apply(lambda z: z.abs().max())
    log('hold-out tube B, largest |z| per fit: ' + '; '.join(f'{st} {sc} {v:.2f}' for (st, sc), v in worst.items()))


STAGES = {'registered': stage_registered, 'throughput': stage_throughput, 'recovery': stage_recovery,
          'baselines': stage_baselines, 'extended': stage_extended, 'refit': stage_refit, 'night': stage_night,
          'traces': stage_traces, 'predictions': stage_predictions}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('stage', choices=list(STAGES))
    parser.add_argument('--workers', type=int, default=12)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--scale', type=float, default=1.0, help='multiplies every Sobol screen (smoke tests: 0.1)')
    parser.add_argument('--profile-nfev', type=int, default=15, help='least-squares iterations per profile point')
    parser.add_argument('--results', default=None, help='results directory (default notebooks/results/05)')
    parser.add_argument('--quick', action='store_true', help='smoke test: tiny screens, two starts, few iterations')
    args = parser.parse_args()
    global RESULTS
    if args.results:
        RESULTS = Path(args.results).resolve()
    if args.quick:
        args.scale, args.profile_nfev = 0.05, 3
        for name in EFFORT:
            EFFORT[name] = (EFFORT[name][0], 2, 5)
    warnings.filterwarnings('ignore', category=RuntimeWarning)
    RESULTS.mkdir(parents=True, exist_ok=True)
    log(f'--- stage {args.stage} (workers {args.workers}) ---')
    t0 = time.perf_counter()
    STAGES[args.stage](args)
    log(f'--- stage {args.stage} done in {(time.perf_counter() - t0) / 60:.1f} min ---')


if __name__ == '__main__':
    main()
