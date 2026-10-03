"""Block 14B — Constraints from observations: which barrier values each observation allows.

An observation is one sample seen at one moment: a mixture, a temperature, the interval [t_min, t_max] in
which the time since mixing lies, and one or more windows on observables that must all hold at that same
moment. A model is consistent with an observation if some time in the interval puts every observable
inside its window. With t_min = t_max this is the check of the Gogoi 2024 control experiments.

Why an interval: for most lab samples the mixing time is not recorded. The file times give a lower bound
(the sample existed before its first spectrum); an upper bound has to be assumed, and the notebook shows how
the bounds move with it. A window that needs little reaction ("TMSPA survives") is decided by t_min; one that
needs enough reaction ("TMSPA gone") is decided by t_max.

    build_share_window          window on a measured share: at least / at most / both, k standard errors wide
    build_observation           one observation in the common format (checks the fields)
    observations_from_controls  the Gogoi 2024 control experiments in that format
    select_observations         the observations that probe one family
    describe_observations       display table of a list of observations
    evaluate_observation        is a model consistent with one observation (and how far off if not)
    scan_barrier                consistency of every observation along a grid of one barrier
    find_allowed_intervals      per observation, the barrier interval it allows (edges refined by bisection)
    intersect_intervals         what is left when every observation must hold
    summarize_region            one line: the interval left, or which observations exclude each other
    trace_boundaries            edges of each observation's region in the plane of two barriers

A barrier is named by the model's family ('hydrolysis', ..., or 'hydrolysis_R1' in a split network), or by
SHARED_BARRIER for one barrier on every reaction (Peter's single E0).

Classification: PUBLIC (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván
"""

from copy import deepcopy

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

from kinetics.constants import ZERO_CELSIUS_K
from kinetics.data.experimental import load_experimental_data
from kinetics.data.network import NETWORK, NETWORK_SPECIES
from kinetics.microkinetics.models import get_model
from kinetics.reactor.batch import simulate_batch_reactor
from kinetics.reactor.validation import OBSERVABLE_REACTIONS, OBSERVABLES

SHARED_BARRIER = 'all reactions'
STATUSES = ('measured', 'derived', 'assumed', 'reading')
DEFAULT_SCAN_GRID_EV = np.round(np.arange(0.60, 1.701, 0.025), 3)
# BDF gives the same observables as Radau and is far faster in scans (see feasibility.SCAN_METHOD)
SCAN_METHOD = 'BDF'


# ------------------------------------------------------------------------------
# Observations
# ------------------------------------------------------------------------------
def build_share_window(share: float, sigma: float, side: str, k_sigma: float = 2.0) -> tuple:
    """(low, high) window on a share in [0, 1]; None is an open side.

    side: 'at least' (share − kσ, open), 'at most' (open, share + kσ) or 'both'.
    """
    if side not in ('at least', 'at most', 'both'):
        raise ValueError("side must be 'at least', 'at most' or 'both'")
    low = max(0.0, share - k_sigma * sigma) if side in ('at least', 'both') else None
    high = min(1.0, share + k_sigma * sigma) if side in ('at most', 'both') else None
    return low, high


def build_observation(obs_id: str, label: str, *, c0_M: dict, T_K: float, t_min_s: float, t_max_s: float,
                      windows: dict, status: str, source: str, time_basis: str, probes=()) -> dict:
    """One observation. c0_M lists the solutes (and EC); the other network species start at zero.

    windows: {observable (kinetics.reactor.OBSERVABLES): (low, high)}, all evaluated at the same time.
    time_basis says where t_min and t_max come from (file times, paper, assumed).
    probes: the family (or families) whose barrier the observation constrains directly; scans of one family
    use the observations that probe it.
    """
    unknown = [name for name in windows if name not in OBSERVABLES]
    if unknown:
        raise KeyError(f"Unknown observable(s) {unknown}. Known: {list(OBSERVABLES)}")
    if status not in STATUSES:
        raise ValueError(f"status must be one of {STATUSES}")
    if not 0.0 <= t_min_s <= t_max_s:
        raise ValueError(f"{obs_id}: need 0 ≤ t_min_s ≤ t_max_s")
    c0 = dict.fromkeys(NETWORK_SPECIES, 0.0)
    c0.update(c0_M)
    return {'id': obs_id, 'label': label, 'c0_M': c0, 'T_K': float(T_K), 't_min_s': float(t_min_s),
            't_max_s': float(t_max_s), 'windows': dict(windows), 'status': status, 'source': source,
            'time_basis': time_basis, 'probes': (probes,) if isinstance(probes, str) else tuple(probes)}


def observations_from_controls(experiments: list = None) -> list:
    """The control experiments of data/experimental_gogoi2024.json as observations (t_min = t_max = t_s)."""
    out = []
    for exp in experiments if experiments is not None else load_experimental_data()['control_experiments']:
        rxn_id = exp.get('reaction', OBSERVABLE_REACTIONS[exp['observable']])
        out.append(build_observation(
            exp['id'], exp['label'], c0_M=exp['c0_M'], T_K=exp['T_K'], t_min_s=exp['t_s'], t_max_s=exp['t_s'],
            windows={exp['observable']: (exp['low'], exp['high'])}, status=exp['status'],
            source='Gogoi 2024 (paper)', time_basis=exp.get('time_basis', 'paper'),
            probes=NETWORK[rxn_id]['class']))
    return out


def select_observations(observations: list, family: str) -> list:
    """The observations that probe one family directly."""
    return [obs for obs in observations if family in obs['probes']]


def _format_hours(h: float) -> str:
    for unit, size in (('d', 24.0), ('h', 1.0), ('min', 1.0 / 60.0)):
        if h >= size or unit == 'min':
            v = h / size
            return f'{v:.0f} {unit}' if v >= 10 else f'{v:.2g} {unit}'


def _format_window(low, high) -> str:
    if low is not None and high is not None:
        return f'{low:.2f}–{high:.2f}'
    return f'≥ {low:.2f}' if low is not None else f'≤ {high:.2f}'


def describe_observations(observations: list) -> pd.DataFrame:
    """One row per observation: what is in the tube, temperature, time interval, windows and provenance."""
    rows = {}
    for obs in observations:
        solutes = ', '.join(f'{sp} {1000.0 * c:.0f} mM' for sp, c in obs['c0_M'].items() if c > 0 and sp != 'EC')
        t_lo, t_hi = obs['t_min_s'] / 3600.0, obs['t_max_s'] / 3600.0
        rows[obs['id']] = {
            'sample': obs['label'],
            'solutes at t = 0': solutes,
            'T (°C)': round(obs['T_K'] - ZERO_CELSIUS_K, 1),
            'time since mixing': _format_hours(t_hi) if t_lo == t_hi else f'{_format_hours(t_lo)} – {_format_hours(t_hi)}',
            'time basis': obs['time_basis'],
            'windows': '; '.join(f'{name} {_format_window(*w)}' for name, w in obs['windows'].items()),
            'status': obs['status'],
            'probes': ', '.join(obs['probes']) or '—',
            'source': obs['source'],
        }
    return pd.DataFrame.from_dict(rows, orient='index')


# ------------------------------------------------------------------------------
# Consistency of a model with an observation
# ------------------------------------------------------------------------------
def _distance_outside(value, low, high):
    """How far value lies outside [low, high] (0 inside); arrays allowed."""
    below = np.maximum(0.0, low - value) if low is not None else 0.0
    above = np.maximum(0.0, value - high) if high is not None else 0.0
    return below + above


def evaluate_observation(obs: dict, model, *, network: dict = None, species_db: dict = None,
                         method: str = SCAN_METHOD, n_times: int = 25) -> dict:
    """Consistency of a model with one observation.

    The batch reactor runs to t_max; the observables are evaluated at n_times log-spaced times in
    [t_min, t_max] (just t_max when the two are equal). Returns
    - 'consistent': some time puts every observable inside its window;
    - 'violation': the smallest summed distance outside the windows over those times (0 if consistent);
    - 't_best_s' and 'values': the time with the smallest violation and the observables there.
    """
    t_min, t_max = obs['t_min_s'], obs['t_max_s']
    single = t_max <= t_min * (1.0 + 1e-9)
    t_start = t_max / 10.0 if single else max(t_min, 1.0)
    sim = simulate_batch_reactor(obs['c0_M'], T_K=obs['T_K'], t_end_s=t_max, model=model, network=network,
                                 species_db=species_db, n_points=2 if single else n_times, t_start_s=t_start,
                                 method=method)
    keep = slice(-1, None) if single else slice(None)
    C, t_s = sim['C_M'][:, keep], sim['t_s'][keep]
    values = {name: np.atleast_1d(OBSERVABLES[name][1](C, sim['idx'], obs['c0_M'])) for name in obs['windows']}
    violation = sum(_distance_outside(values[name], low, high) for name, (low, high) in obs['windows'].items())
    best = int(np.argmin(violation))
    return {'consistent': bool(violation[best] <= 0.0), 'violation': float(violation[best]),
            't_best_s': float(t_s[best]), 'values': {name: float(v[best]) for name, v in values.items()}}


def _with_barrier(spec, parameter: str, g_eV: float):
    if parameter == SHARED_BARRIER:
        return spec.with_barrier(g_eV, name=f'{spec.name} [all g = {g_eV:.3f}]')
    return spec.with_family_params(parameter, g_eV=g_eV, name=f'{spec.name} [{parameter} = {g_eV:.3f}]')


# ------------------------------------------------------------------------------
# One barrier at a time
# ------------------------------------------------------------------------------
def scan_barrier(model, parameter: str, observations: list, *, g_grid_eV=DEFAULT_SCAN_GRID_EV,
                 network: dict = None, species_db: dict = None) -> pd.DataFrame:
    """Consistency of every observation along a grid of one barrier; the other barriers keep the model's values.

    Returns one row per (g, observation) with 'consistent' and 'violation'.
    """
    spec = get_model(model)
    rows = []
    for g in g_grid_eV:
        variant = _with_barrier(spec, parameter, float(g))
        for obs in observations:
            check = evaluate_observation(obs, variant, network=network, species_db=species_db)
            rows.append({'parameter': parameter, 'g_eV': float(g), 'observation': obs['id'],
                         'probes': ', '.join(obs['probes']), 'status': obs['status'],
                         'consistent': check['consistent'], 'violation': check['violation']})
    return pd.DataFrame(rows, columns=['parameter', 'g_eV', 'observation', 'probes', 'status', 'consistent',
                                       'violation'])


def _bisect_edge(spec, parameter, obs, g_out, g_in, tol_eV, network, species_db):
    """Barrier between g_out (inconsistent) and g_in (consistent) where consistency switches."""
    while abs(g_in - g_out) > tol_eV:
        g_mid = 0.5 * (g_in + g_out)
        ok = evaluate_observation(obs, _with_barrier(spec, parameter, g_mid), network=network,
                                  species_db=species_db)['consistent']
        g_in, g_out = (g_mid, g_out) if ok else (g_in, g_mid)
    return 0.5 * (g_in + g_out)


def _refine_near_best(spec, parameter, obs, g, violation, tol_eV, network, species_db):
    """Grid with one point added where the violation is smallest between the neighbours of its grid minimum."""
    i = int(np.argmin(violation))
    lo, hi = g[max(i - 1, 0)], g[min(i + 1, len(g) - 1)]
    res = minimize_scalar(
        lambda x: evaluate_observation(obs, _with_barrier(spec, parameter, float(x)), network=network,
                                       species_db=species_db)['violation'],
        bounds=(lo, hi), method='bounded', options={'xatol': tol_eV})
    if res.fun > 0.0:
        return g, np.zeros(len(g), dtype=bool)
    order = np.argsort(np.append(g, res.x))
    return np.append(g, res.x)[order], np.append(np.zeros(len(g), dtype=bool), True)[order]


def find_allowed_intervals(model, parameter: str, observations: list, *, scan: pd.DataFrame = None,
                           g_grid_eV=DEFAULT_SCAN_GRID_EV, tol_eV: float = 0.002, network: dict = None,
                           species_db: dict = None) -> pd.DataFrame:
    """Per observation, the interval of one barrier it allows, with the edges refined by bisection to tol_eV.

    An edge at the end of the grid is open ('scan limit'): the observation does not bound that side within
    the scanned range. 'no value' means no grid point is consistent. If the consistent points form several
    runs, the envelope is returned and flagged.
    """
    spec = get_model(model)
    by_id = {obs['id']: obs for obs in observations}
    if scan is None:
        scan = scan_barrier(spec, parameter, observations, g_grid_eV=g_grid_eV, network=network,
                            species_db=species_db)
    rows = []
    for obs_id, grp in scan[scan['parameter'] == parameter].groupby('observation', sort=False):
        grp = grp.sort_values('g_eV')
        obs = by_id[obs_id]
        g, ok = grp['g_eV'].to_numpy(), grp['consistent'].to_numpy()
        row = {'parameter': parameter, 'observation': obs_id, 'probes': ', '.join(obs['probes']),
               'status': obs['status'], 'time_basis': obs['time_basis']}
        if not ok.any():
            # A window narrower than the grid step can fall between two points: look for it near the point
            # that comes closest, by minimising the violation between its neighbours
            g, ok = _refine_near_best(spec, parameter, obs, g, grp['violation'].to_numpy(), tol_eV, network,
                                      species_db)
        if not ok.any():
            rows.append({**row, 'g_low_eV': np.nan, 'g_high_eV': np.nan, 'low_set_by': 'no value',
                         'high_set_by': 'no value', 'verdict': 'no value in the scanned range'})
            continue
        first, last = np.flatnonzero(ok)[[0, -1]]
        runs = 1 + int(np.sum(np.diff(ok[first:last + 1].astype(int)) == 1))
        lo = _bisect_edge(spec, parameter, obs, g[first - 1], g[first], tol_eV, network, species_db) \
            if first > 0 else np.nan
        hi = _bisect_edge(spec, parameter, obs, g[last + 1], g[last], tol_eV, network, species_db) \
            if last < len(g) - 1 else np.nan
        verdict = {(True, True): 'bounded on both sides', (True, False): 'lower bound only',
                   (False, True): 'upper bound only', (False, False): 'any value in the scanned range'}[
            (bool(np.isfinite(lo)), bool(np.isfinite(hi)))]
        rows.append({**row, 'g_low_eV': lo, 'g_high_eV': hi,
                     'low_set_by': 'bound' if np.isfinite(lo) else 'scan limit',
                     'high_set_by': 'bound' if np.isfinite(hi) else 'scan limit',
                     'verdict': verdict + (' (several intervals: envelope)' if runs > 1 else '')})
    return pd.DataFrame(rows, columns=['parameter', 'observation', 'probes', 'status', 'time_basis', 'g_low_eV',
                                       'g_high_eV', 'low_set_by', 'high_set_by', 'verdict'])


def intersect_intervals(intervals: pd.DataFrame, g_range_eV: tuple = (DEFAULT_SCAN_GRID_EV[0],
                                                                       DEFAULT_SCAN_GRID_EV[-1])) -> pd.DataFrame:
    """Per barrier: the interval left when every observation must hold, and the observations that set its ends.

    Open edges count as the ends of g_range_eV. 'empty' means two observations exclude each other.
    """
    rows = {}
    for parameter, grp in intervals.groupby('parameter', sort=False):
        if (grp['verdict'] == 'no value in the scanned range').any():
            killers = ', '.join(grp.loc[grp['verdict'] == 'no value in the scanned range', 'observation'])
            rows[parameter] = {'g_low_eV': np.nan, 'g_high_eV': np.nan, 'low_set_by': killers,
                               'high_set_by': killers, 'empty': True}
            continue
        lows, highs = grp['g_low_eV'].fillna(g_range_eV[0]), grp['g_high_eV'].fillna(g_range_eV[1])
        lo, hi = lows.max(), highs.min()
        rows[parameter] = {
            'g_low_eV': lo, 'g_high_eV': hi,
            'low_set_by': grp.loc[lows.idxmax(), 'observation'] if np.isfinite(grp['g_low_eV']).any() else 'scan limit',
            'high_set_by': grp.loc[highs.idxmin(), 'observation'] if np.isfinite(grp['g_high_eV']).any() else 'scan limit',
            'empty': bool(lo > hi),
        }
    return pd.DataFrame.from_dict(rows, orient='index')


def summarize_region(intervals: pd.DataFrame) -> str:
    """One line on what the observations leave of one barrier, or why nothing is left.

    Observations that allow no value are named first; among the others, the pair whose bounds cross is named.
    """
    none = intervals.loc[intervals['verdict'].str.startswith('no value'), 'observation'].tolist()
    rest = intervals[~intervals['observation'].isin(none)]
    reasons = [f"no value allowed by {', '.join(none)}"] if none else []
    if len(rest):
        inter = intersect_intervals(rest).iloc[0]
        lo, hi = inter['g_low_eV'], inter['g_high_eV']
        if inter['empty']:
            reasons.append(f"{inter['high_set_by']} (≤ {hi:.3f}) excludes {inter['low_set_by']} (≥ {lo:.3f})")
        elif not none:
            if inter['low_set_by'] == 'scan limit':
                return f"g ≤ {hi:.3f} eV ({inter['high_set_by']})"
            if inter['high_set_by'] == 'scan limit':
                return f"g ≥ {lo:.3f} eV ({inter['low_set_by']})"
            return f"{lo:.3f}–{hi:.3f} eV ({inter['low_set_by']} … {inter['high_set_by']})"
    return 'empty: ' + '; '.join(reasons)


# ------------------------------------------------------------------------------
# Two barriers at a time
# ------------------------------------------------------------------------------
def trace_boundaries(model, x: tuple, y: tuple, observations: list, *, tol_eV: float = 0.002,
                     network: dict = None, species_db: dict = None) -> pd.DataFrame:
    """Edges of the region each observation allows in the plane of two barriers.

    x = (parameter, scan grid), y = (parameter, values). At each y value the x-interval of every observation is
    found as in find_allowed_intervals (edges refined to tol_eV), so the edges are exact in x and sampled in y.
    The other barriers keep the model's values. Returns the find_allowed_intervals table with 'y_parameter' and
    'y_eV' added.
    """
    spec = get_model(model)
    (px, gx), (py, gy) = x, y
    frames = [find_allowed_intervals(_with_barrier(spec, py, float(vy)), px, observations, g_grid_eV=gx,
                                     tol_eV=tol_eV, network=network, species_db=species_db)
              .assign(y_parameter=py, y_eV=float(vy)) for vy in gy]
    return pd.concat(frames, ignore_index=True)


def replace_observation(obs: dict, **changes) -> dict:
    """Copy of an observation with some fields changed (e.g. t_max_s for a time scenario)."""
    out = deepcopy(obs)
    out.update(changes)
    return out
