"""Block 14F — Estimation: best fit of a candidate structure and profile intervals.

The objective is the χ² of kinetics.fitting.residuals (measured shares) plus, for
freed reaction energies, their distance from the computed values in the metric of the computed covariance.
Unknown ages are a scenario of the sample set, not a prior. An amount at mixing that a recipe does not fix is a
parameter with a declared range (kinetics.fitting.candidates).

    FitProblem                     a structure and a sample set: parameters, box, residual vector
    fit_structure                  Sobol screen of the box, then bounded least squares from the best points
    calculate_profile              χ² re-minimised over the other parameters along a grid of one parameter
    find_confidence_interval       where a profile rises by Δχ² (3.84 = 95 % for one parameter); open sides stay open
    write_fit_result, load_fit_result  a fit result as JSON (the residual table as records)

Fit rule (fixed before fitting): a structure fits a scenario if no standardised residual exceeds 3 in absolute
value. The pull of each freed reaction energy from its computed value counts as a residual.

Search and derivatives use a loose relative solver tolerance (rtol 1e-6); the reported χ², residual table and
final polish use rtol 1e-8. The absolute tolerance is tight in both (see residuals.SEARCH_SOLVER). Parameters are searched in box coordinates u = 1 + (x − lower)/(upper − lower), so that a
relative finite-difference step is the same fraction of every box (about 1 meV for a barrier).

Source: Y. Alcaraz Galván
"""

import json
import multiprocessing as mp
import os
import time

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.stats import qmc

from kinetics.fitting.candidates import get_structure, scale_energy_prior
from kinetics.fitting.residuals import REPORT_SOLVER, SEARCH_SOLVER, calculate_residuals, tabulate_residuals

FIT_RULE_MAX_Z = 3.0
DELTA_CHI2_95 = 3.84
DIFF_STEP = 1.0e-3
BARRIER_OFFSETS_EV = (-0.30, -0.20, -0.15, -0.10, -0.075, -0.05, -0.025, 0.025, 0.05, 0.075, 0.10, 0.15, 0.20, 0.30)
ENERGY_OFFSETS_SE = (-2.0, -1.5, -1.0, -0.75, -0.5, -0.25, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0)
RECIPE_FACTORS = (0.25, 0.5, 0.67, 0.8, 0.9, 1.1, 1.25, 1.5, 2.0, 4.0)      # multiples of the best amount at mixing


class FitProblem:
    """A structure and a sample set: the free parameters, their box and the residual vector."""

    def __init__(self, structure, samples: list, *, solver: dict = None, fixed: dict = None, prior_scale: float = 1.0):
        self.structure = get_structure(structure)
        self.samples = samples
        self.solver = dict(solver or SEARCH_SOLVER)
        self.prior_scale = float(prior_scale)
        self.all_parameters = self.structure.parameters(samples)
        self.fixed = dict(fixed or {})
        unknown = set(self.fixed) - {p.name for p in self.all_parameters}
        if unknown:
            raise KeyError(f'Cannot fix unknown parameter(s) {sorted(unknown)}')
        self.parameters = [p for p in self.all_parameters if p.name not in self.fixed]
        self.names = [p.name for p in self.parameters]
        self.lower = np.array([p.lower for p in self.parameters], dtype=float)
        self.upper = np.array([p.upper for p in self.parameters], dtype=float)

    @property
    def n(self) -> int:
        return len(self.parameters)

    def theta(self, x) -> dict:
        """{name: value} of every parameter, fixed ones included."""
        values = {**self.fixed, **dict(zip(self.names, np.asarray(x, dtype=float)))}
        return {p.name: float(values[p.name]) for p in self.all_parameters}

    def vector(self, theta: dict) -> np.ndarray:
        """Free-parameter vector from {name: value}, clipped to the box."""
        return np.clip([theta[name] for name in self.names], self.lower, self.upper)

    def to_unit(self, x) -> np.ndarray:
        return 1.0 + (np.asarray(x, dtype=float) - self.lower) / (self.upper - self.lower)

    def from_unit(self, u) -> np.ndarray:
        return self.lower + (np.asarray(u, dtype=float) - 1.0) * (self.upper - self.lower)

    def build(self, x) -> dict:
        theta = self.theta(x)
        return scale_energy_prior(self.structure, self.structure.build(theta), theta, self.prior_scale)

    def residuals(self, x) -> np.ndarray:
        built = self.build(x)
        data = calculate_residuals(built['model'], self.samples, added_M=built.get('added_M'), solver=self.solver)
        return np.concatenate([data, built['prior']])

    def chi2(self, x) -> float:
        r = self.residuals(x)
        return float(r @ r)

    def table(self, x) -> pd.DataFrame:
        """Standardised residuals of every measured quantity and of the freed energies, at this problem's solver."""
        built = self.build(x)
        table = tabulate_residuals(built['model'], self.samples, added_M=built.get('added_M'), solver=self.solver)
        if built['prior_rows']:
            table = pd.concat([table, pd.DataFrame(built['prior_rows'])], ignore_index=True)
        return table

    def replace(self, *, samples: list = None, solver: dict = None, fixed: dict = None, prior_scale: float = None):
        return FitProblem(self.structure, self.samples if samples is None else samples,
                          solver=self.solver if solver is None else solver,
                          fixed=self.fixed if fixed is None else fixed,
                          prior_scale=self.prior_scale if prior_scale is None else prior_scale)


# ------------------------------------------------------------------------------
# Parallel helpers: the problem is inherited by forked workers, never pickled per task
# ------------------------------------------------------------------------------
_ACTIVE = None


def _chi2_unit(u):
    return _ACTIVE.chi2(_ACTIVE.from_unit(u))


def _polish(problem: FitProblem, u0, max_nfev: int) -> dict:
    if problem.n == 0:
        return {'u': np.zeros(0), 'chi2': problem.chi2(np.zeros(0)), 'nfev': 1}
    res = least_squares(lambda u: problem.residuals(problem.from_unit(u)), np.clip(u0, 1.0, 2.0), bounds=(1.0, 2.0),
                        method='trf', diff_step=DIFF_STEP, ftol=1e-9, xtol=1e-7, gtol=1e-9, max_nfev=max_nfev)
    return {'u': res.x, 'chi2': float(2.0 * res.cost), 'nfev': int(res.nfev * (problem.n + 1))}


def _polish_unit(args):
    u0, max_nfev = args
    return _polish(_ACTIVE, u0, max_nfev)


def _map(function, items: list, problem: FitProblem, workers: int) -> list:
    """function over items with `problem` as the active problem, in `workers` forked processes."""
    global _ACTIVE
    previous, _ACTIVE = _ACTIVE, problem
    try:
        if workers <= 1 or len(items) <= 1:
            return [function(item) for item in items]
        with mp.get_context('fork').Pool(min(workers, len(items))) as pool:
            return pool.map(function, items, chunksize=1)
    finally:
        _ACTIVE = previous


# ------------------------------------------------------------------------------
# Best fit
# ------------------------------------------------------------------------------
def _select_starts(points: np.ndarray, chi2: np.ndarray, n_starts: int, min_distance: float) -> list:
    """The best points, skipping any closer than min_distance (box units) to one already chosen."""
    chosen = []
    for i in np.argsort(chi2):
        if all(np.linalg.norm(points[i] - points[j]) > min_distance for j in chosen):
            chosen.append(int(i))
        if len(chosen) == n_starts:
            break
    return chosen


def fit_structure(problem: FitProblem, *, n_screen: int = 512, n_starts: int = 16, max_nfev: int = 40, seed: int = 0,
                  workers: int = 1, start_points=(), seed_values: dict = None, report: bool = True) -> dict:
    """Best fit of a problem: Sobol screen of the box, bounded least squares from the best points, final polish.

    n_screen: number of screen points (rounded up to a power of two). n_starts: how many of the best screen
    points are polished. start_points: parameter sets {name: value} or vectors polished in addition (e.g. the
    best fit of the structure this one extends). max_nfev: least-squares iterations per start.
    seed_values: {name: value} known from a part of the data that sees those parameters alone (e.g. the
    barriers the TMSOH-only samples fix). The screen is then repeated with these coordinates held at their
    values, and the best half as many of those points are polished too; every parameter is free in the polish.
    With report: the best point is polished once more and tabulated with the report solver (rtol 1e-8).

    Returns a dict: 'theta', 'x', 'chi2', 'n_residuals', 'n_parameters', 'max_abs_z', 'worst', 'fits' (the fit
    rule), 'table' (standardised residuals), 'ages', 'local_optima', 'at_bounds', 'n_evaluations', 'seconds'.
    """
    t_start = time.perf_counter()
    n = problem.n
    starts_u, n_eval = [], 0
    if n > 0 and n_screen > 0:
        screen = 1.0 + qmc.Sobol(n, scramble=True, seed=seed).random_base2(int(np.ceil(np.log2(n_screen))))
        screen_chi2 = np.array(_map(_chi2_unit, list(screen), problem, workers))
        n_eval += len(screen)
        starts_u = [screen[i] for i in _select_starts(screen, screen_chi2, n_starts, 0.05 * np.sqrt(n))]
        held = {problem.names.index(k): v for k, v in (seed_values or {}).items() if k in problem.names}
        if held:
            seeded = screen.copy()
            for i, value in held.items():
                seeded[:, i] = 1.0 + (np.clip(value, problem.lower[i], problem.upper[i]) - problem.lower[i]) \
                    / (problem.upper[i] - problem.lower[i])
            seeded_chi2 = np.array(_map(_chi2_unit, list(seeded), problem, workers))
            n_eval += len(seeded)
            starts_u += [seeded[i] for i in _select_starts(seeded, seeded_chi2, max(n_starts // 2, 1), 0.05 * np.sqrt(n))]
    for point in start_points:
        x = problem.vector(point) if isinstance(point, dict) else np.asarray(point, dtype=float)
        starts_u.append(problem.to_unit(np.clip(x, problem.lower, problem.upper)))
    if not starts_u:
        starts_u = [np.full(n, 1.5)]
    polished = _map(_polish_unit, [(u, max_nfev) for u in starts_u], problem, workers)
    n_eval += sum(p['nfev'] for p in polished)
    polished.sort(key=lambda p: p['chi2'])
    best = polished[0]
    final = problem.replace(solver=REPORT_SOLVER) if report else problem
    if report:
        best = _polish(final, best['u'], 15)
        n_eval += best['nfev']
    x = final.from_unit(best['u'])
    table = final.table(x)
    residuals = final.residuals(x)
    worst = table.loc[table['z'].abs().idxmax()]
    ages = {name: age for name, age in zip(table['sample'], table['age_h']) if np.isfinite(age)}
    return {
        'structure': problem.structure.name, 'parameters': list(problem.names), 'fixed': dict(problem.fixed),
        'theta': final.theta(x), 'x': [float(v) for v in x], 'chi2': float(residuals @ residuals),
        'n_residuals': int(residuals.size), 'n_parameters': n, 'max_abs_z': float(abs(worst['z'])),
        'worst': f"{worst['sample']}: {worst['quantity']}", 'fits': bool(abs(worst['z']) <= FIT_RULE_MAX_Z),
        'table': table, 'ages': ages,
        'local_optima': [{'theta': problem.theta(problem.from_unit(p['u'])), 'chi2': p['chi2']} for p in polished],
        'at_bounds': [name for name, u in zip(problem.names, best['u']) if u < 1.0 + 2e-3 or u > 2.0 - 2e-3],
        'n_evaluations': int(n_eval), 'seconds': time.perf_counter() - t_start,
    }


# ------------------------------------------------------------------------------
# Profiles and intervals
# ------------------------------------------------------------------------------
def build_profile_grid(problem: FitProblem, theta: dict, name: str) -> np.ndarray:
    """Grid of one parameter around its best value: fixed offsets for a barrier (±0.30 eV), multiples of the
    standard error for a freed energy (±2 SE), multiples of the best value for an amount at mixing (the whole
    box when it is zero)."""
    parameter = next(p for p in problem.all_parameters if p.name == name)
    best = theta[name]
    if parameter.kind == 'barrier':
        grid = best + np.array(BARRIER_OFFSETS_EV)
    elif parameter.kind == 'energy':
        sigma = next(row['sigma'] for row in problem.build(problem.vector(theta))['prior_rows']
                     if parameter.target in row['quantity'])
        grid = best + sigma * np.array(ENERGY_OFFSETS_SE)
    elif parameter.kind == 'recipe' and best > 1e-6:
        grid = best * np.array(RECIPE_FACTORS)
    else:
        grid = np.linspace(parameter.lower, parameter.upper, 13)
    grid = grid[(grid >= parameter.lower - 1e-12) & (grid <= parameter.upper + 1e-12)]
    return np.unique(np.round(grid[np.abs(grid - best) > 1e-9], 9))


def _profile_chain(args) -> list:
    """Re-minimised χ² along one direction of a profile grid, each point started from the previous solution
    and from the best fit; the lower χ² is kept."""
    name, values, theta_best, max_nfev = args
    rows, theta_previous = [], dict(theta_best)
    for value in values:
        sub = _ACTIVE.replace(fixed={**_ACTIVE.fixed, name: float(value)})
        candidates = [_polish(sub, sub.to_unit(sub.vector(theta_previous)), max_nfev)]
        if theta_previous != theta_best:
            candidates.append(_polish(sub, sub.to_unit(sub.vector(theta_best)), max_nfev))
        best = min(candidates, key=lambda c: c['chi2'])
        theta_previous = sub.theta(sub.from_unit(best['u']))
        rows.append({'parameter': name, 'value': float(value), 'chi2': best['chi2'],
                     'nfev': sum(c['nfev'] for c in candidates), **{f'θ {k}': v for k, v in theta_previous.items()}})
    return rows


def calculate_profile(problem: FitProblem, fit: dict, names=None, *, grids: dict = None, max_nfev: int = 20,
                      workers: int = 1) -> pd.DataFrame:
    """Profile of χ² along each named parameter: the others are re-minimised at every grid value.

    The grid of a parameter (build_profile_grid, or grids[name]) is walked outward from the best value in both
    directions, so each sub-fit starts close to its solution. Returns one row per (parameter, value) with the
    re-minimised χ², 'delta' = χ² − χ² of the fit, and the other parameters there ('θ name').
    """
    names = list(names) if names is not None else [p.name for p in problem.parameters]
    theta = fit['theta']
    reference = problem.chi2(problem.vector(theta))            # the fit's χ² at this problem's solver tolerance
    chains = []
    for name in names:
        grid = np.asarray(grids[name], dtype=float) if grids and name in grids else build_profile_grid(problem, theta, name)
        below, above = grid[grid < theta[name] - 1e-12][::-1], grid[grid > theta[name] + 1e-12]
        chains += [(name, values, theta, max_nfev) for values in (below, above) if len(values)]
    rows = [row for chain in _map(_profile_chain, chains, problem, workers) for row in chain]
    rows += [{'parameter': name, 'value': theta[name], 'chi2': reference, 'nfev': 0,
              **{f'θ {k}': v for k, v in theta.items()}} for name in names]
    profile = pd.DataFrame(rows).sort_values(['parameter', 'value']).reset_index(drop=True)
    profile['delta'] = profile['chi2'] - reference
    return profile


def find_confidence_interval(profile: pd.DataFrame, name: str, delta: float = DELTA_CHI2_95) -> dict:
    """Interval of one parameter over which its profile stays within `delta` of its lowest χ².

    The interval runs from the lowest to the highest profile value that is within `delta`; 'separate ranges' is
    True when values between them are not (several basins), in which case the interval is their envelope.
    An edge is placed between the outermost accepted value and its rejected neighbour by interpolating
    sqrt(Δχ²), which is linear for a parabola. A side whose last grid value is still accepted is open (NaN).
    'status' is 'interval', 'upper bound only', 'lower bound only' or 'not determined'. 'best' is the value
    with the lowest χ²; 'lower than fit' flags a profile that found a χ² more than 0.5 below the fit.
    """
    rows = profile[profile['parameter'] == name]
    p = rows.groupby('value', as_index=False)['chi2'].min().sort_values('value')
    value, chi2 = p['value'].to_numpy(), p['chi2'].to_numpy()
    rise = chi2 - chi2.min()
    accepted = np.flatnonzero(rise <= delta)
    i_low, i_high = accepted[0], accepted[-1]

    def edge(i_out, i_in):
        root_in, root_out = np.sqrt(max(rise[i_in], 0.0)), np.sqrt(rise[i_out])
        f = (np.sqrt(delta) - root_in) / (root_out - root_in)
        return float(value[i_in] + f * (value[i_out] - value[i_in]))

    low = edge(i_low - 1, i_low) if i_low > 0 else np.nan
    high = edge(i_high + 1, i_high) if i_high < len(value) - 1 else np.nan
    status = {(True, True): 'interval', (False, True): 'upper bound only', (True, False): 'lower bound only',
              (False, False): 'not determined'}[(bool(np.isfinite(low)), bool(np.isfinite(high)))]
    return {'parameter': name, 'best': float(value[int(np.argmin(chi2))]), 'low': low, 'high': high, 'status': status,
            'grid_low': float(value[0]), 'grid_high': float(value[-1]),
            'separate ranges': bool((rise[i_low:i_high + 1] > delta).any()),
            'lower than fit': bool(rows['delta'].min() < -0.5)}


# ------------------------------------------------------------------------------
# Results on disk
# ------------------------------------------------------------------------------
def _plain(value):
    if isinstance(value, pd.DataFrame):
        return {'__table__': json.loads(value.to_json(orient='records', double_precision=15))}
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, np.ndarray):
        return [_plain(v) for v in value.tolist()]
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        return value.item()
    return value


def _restore(value):
    if isinstance(value, dict):
        if set(value) == {'__table__'}:
            # JSON has no NaN: a missing number comes back as None, and a column of them would stay of object type
            table = pd.DataFrame(value['__table__'])
            return table.where(table.notna(), np.nan).infer_objects()
        return {k: _restore(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_restore(v) for v in value]
    return value


def write_fit_result(result: dict, path: str) -> str:
    """Save a result dict (fit, profile bundle, ...) as JSON; tables become lists of records. Written atomically."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    temporary = f'{path}.tmp'
    with open(temporary, 'w', encoding='utf-8') as f:
        json.dump(_plain(result), f, indent=1, ensure_ascii=False)
    os.replace(temporary, path)
    return path


def load_fit_result(path: str) -> dict:
    """Result dict saved by write_fit_result, with its tables as DataFrames again."""
    with open(path, 'r', encoding='utf-8') as f:
        return _restore(json.load(f))
