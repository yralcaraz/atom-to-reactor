"""Block 14D — Residuals: the lab observables against a model, under one reading of the unknown times.

A sample is one tube: a composition at mixing, a temperature history and one or more spectra, each a set of
area shares with their uncertainties. The spectra of one nucleus form a block; a tube recorded on two nuclei
(³¹P and the ¹³C Si–CH3 region) carries two blocks, which share the composition, the history and the age and
are read off one simulation. The time from mixing to a sample's first spectrum (its age) is not
recorded, so every result is computed under a scenario:

    'short', 'middle', 'long'   every unknown age is 1 h, 1 day or 7 days
    'free'                      each unknown age lies between 1 h and 30 days and takes the value that favours
                                the model (the reading notebook 03 uses when an observation picks its best time)

A recipe can hold a species whose amount at mixing is not known (the water of the tube mixed without water):
the sample lists it under 'unknown_M' with its range, a fit takes it as a parameter, and every function that
simulates a sample takes the amounts as added_M = {sample: {species: mol/L}}, on top of the recipe.

Error model of a share: noise ⊕ baseline spread ⊕ a declared floor. Repeated room-temperature spectra of one
block (the 2 % water sample) are one composition and n − 1 measurements of "no drift": the baseline spread and
the floor are common to them, the noise (times the Birge ratio of their scatter) is independent.

The K shares of a spectrum sum to one, so their K residuals are scaled by sqrt((K − 1)/K): that is the exact
χ² of K − 1 independent values when the K errors are equal.

    build_sample_set            observables + scenario → the samples a fit compares with
    build_sample_history        temperature segments from mixing and the spectrum times of one sample
    calculate_predicted_shares  the shares a model predicts for one sample
    calculate_residuals         whitened residual vector of a model over a sample set (χ² = Σ r²)
    find_free_age               the age of the heated sample that a given model prefers (free scenario)
    tabulate_residuals          the standardised residuals, one row per measured quantity
    summarize_residuals         χ², number of values and the largest standardised residual

Source: Y. Alcaraz Galván
"""

from copy import deepcopy

import numpy as np
import pandas as pd
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import minimize_scalar

from kinetics.constants import ZERO_CELSIUS_K
from kinetics.data.experimental import load_experimental_data
from kinetics.data.observables import DEFAULT_ERROR_FLOOR, RT_MAX_C, calculate_replicate_scatter, list_spectrum_sets
from kinetics.reactor.protocol import simulate_history

SCENARIO_AGES_H = {'short': 1.0, 'middle': 24.0, 'long': 168.0}
# Range of an unknown age in the free scenario: 1 h to 30 days. DECLARED. The upper limit was 7 days in the
# first fit, where the 2 % water sample sat on it in every free-age fit
FREE_AGE_BOUNDS_H = (1.0, 720.0)
SCENARIOS = ('short', 'middle', 'long', 'free')
# In the free scenario: samples whose age is found inside each evaluation (one spectrum set, no heating), the
# sample whose age is a fit parameter (heated after its first spectrum), and the rest (age at the lower bound)
FREE_AGE_PARAMETER_SAMPLE = '2 % H2O'

# Species behind each window: {nucleus: {window: {species: nuclei of that window per molecule}}}
SHARE_SPECIES = {
    '31P': {'TMSPA': {'TMSPA': 1}, 'BMSPA': {'BMSPA': 1}, 'MMSPA': {'MMSPA': 1}, 'H3PO4': {'H3PO4': 1}},
    # Si–CH3 carbons, three per Si. The TMSOEG window also holds the second ring-opening product (notebook 03 §2.1).
    # 'P-silyl': the silyl groups still on phosphate, a window of the tubes that hold phosphate only
    '13C': {'TMSOH': {'TMSOH': 1}, 'HMDSO': {'HMDSO': 2}, 'TMSOEG': {'TMSOEG': 1, 'TMSOdiEG': 1},
            'P-silyl': {'TMSPA': 3, 'BMSPA': 2, 'MMSPA': 1}},
}
PAPER_OBSERVATIONS = ('E1',)
PAPER_SOFTNESS = 0.5          # softness of a one-sided paper window, as a fraction of its limit. DECLARED
REPLICATE_MIN_SPECTRA = 3
PENALTY = 1.0e3               # residual given to every value of a sample whose simulation failed
# The absolute tolerance has to sit far below a nanomolar: some fits explain a sample by an induction time that
# is seeded by traces of that size, and with atol = 1e-10 M their χ² came out wrong by a factor of three
SEARCH_SOLVER = {'method': 'BDF', 'rtol': 1e-6, 'atol': 1e-13, 'max_rhs_calls': 30000}
REPORT_SOLVER = {'method': 'BDF', 'rtol': 1e-8, 'atol': 1e-14, 'max_rhs_calls': 300000}
N_AGE_GRID = 53               # ages tried inside an evaluation: the step in log(age) of the first fit (41 points to 7 days)


# ------------------------------------------------------------------------------
# Samples
# ------------------------------------------------------------------------------
def _floor(share: np.ndarray, floor: dict, scale: float) -> np.ndarray:
    return scale * np.maximum(floor['absolute'], floor['relative'] * np.maximum(share, 0.0))


def _lab_block(name: str, nucleus: str, role: str, spectra: list, observables: dict, error_floor: dict,
               floor_scale: float) -> dict:
    """The spectra of one nucleus of one tube: its windows, shares and error model."""
    unknown = set(spectra[0]['shares']) - set(SHARE_SPECIES[nucleus])
    if unknown:
        raise KeyError(f'{name}: no species declared for the {nucleus} window(s) {sorted(unknown)}')
    windows = tuple(w for w in SHARE_SPECIES[nucleus] if w in spectra[0]['shares'])
    Y = np.array([[sp['shares'][w]['share'] for w in windows] for sp in spectra])
    noise = np.array([[sp['shares'][w]['sigma_noise'] for w in windows] for sp in spectra])
    baseline = np.array([[sp['shares'][w]['sigma_baseline'] for w in windows] for sp in spectra])
    all_rt = all(sp['temperature_C'] < RT_MAX_C for sp in spectra)
    replicate = all_rt and len(spectra) >= REPLICATE_MIN_SPECTRA
    floor = error_floor[nucleus]
    if replicate:
        birge = max(1.0, calculate_replicate_scatter(observables, name, nucleus)['birge'])
        sigma_ind = birge * noise
        sigma_common = np.hypot(np.median(baseline, axis=0), _floor(Y.mean(axis=0), floor, floor_scale))
    else:
        birge = 1.0
        sigma_ind = np.sqrt(noise ** 2 + baseline ** 2 + _floor(Y, floor, floor_scale) ** 2)
        sigma_common = np.zeros(len(windows))
    return {
        'nucleus': nucleus, 'role': role,
        't_spec_s': 3600.0 * np.array([sp['t_since_first_h'] for sp in spectra]),
        'windows': windows, 'measured': Y, 'sigma_ind': sigma_ind, 'sigma_common': sigma_common,
        'replicate': replicate, 'birge': birge,
        'cholesky': [cholesky(np.diag(sigma_ind[:, k] ** 2) + sigma_common[k] ** 2, lower=True)
                     for k in range(len(windows))] if replicate else None,
    }


def _lab_sample(name: str, record: dict, observables: dict, error_floor: dict, floor_scale: float, roles: tuple):
    """One tube with the blocks whose role is in `roles`; None when it has none."""
    blocks = [_lab_block(name, nucleus, role, spectra, observables, error_floor, floor_scale)
              for nucleus, role, spectra in list_spectrum_sets(record) if role in roles]
    if not blocks:
        return None
    c0 = dict(record['c0_M'])
    return {
        'name': name, 'kind': 'lab', 'c0_M': c0, 'wet': c0.get('H2O', 0.0) > 0.0,
        'T_rt_K': record['T_rt_C'] + ZERO_CELSIUS_K,
        'pre_segments': [(h['T_C'] + ZERO_CELSIUS_K, 3600.0 * h['duration_h']) for h in record['pre_history']],
        'heated': sorted((3600.0 * h['start_h'], h['T_C'] + ZERO_CELSIUS_K, 3600.0 * h['duration_h'])
                         for h in record['heated']),
        'blocks': blocks, 'age_mode': 'none', 'age_h': None,
        'unknown_M': {sp: (float(low), float(high)) for sp, (low, high) in record.get('unknown_M', {}).items()},
    }


def _paper_sample(obs_id: str) -> dict:
    exp = next(e for e in load_experimental_data()['control_experiments'] if e['id'] == obs_id)
    limit = exp['high'] if exp['high'] is not None else exp['low']
    return {'name': obs_id, 'kind': 'paper', 'c0_M': dict(exp['c0_M']), 'wet': False,
            'T_K': float(exp['T_K']), 't_s': float(exp['t_s']), 'observable': exp['observable'],
            'low': exp['low'], 'high': exp['high'], 'softness': PAPER_SOFTNESS * limit,
            'age_mode': 'none', 'age_h': None, 'label': exp['label']}


def build_sample_set(observables: dict, scenario: str, *, roles=('fit',), exclude=(), paper=PAPER_OBSERVATIONS,
                     error_floor: dict = None, floor_scale: float = 1.0, overrides: dict = None) -> list:
    """The samples a model is compared with, under one scenario of the unknown ages.

    roles: which blocks of the lab samples enter ('fit', 'hold_out', 'check'); a sample enters with the blocks
    that have one of these roles. exclude: sample names left out (leave-one-out).
    paper: one-sided observations of Gogoi 2024 added to the set (ids of data/experimental_gogoi2024.json).
    floor_scale multiplies the declared error floor (0 removes it).
    overrides: {sample: {'hold_h': hours of its known pre-history, 'heated_duration_h': hours of every heating
    episode (at most until the next one starts), 'added_M': {species: mol/L present at mixing on top of the recipe}}}, for the sensitivity to the
    times the files do not fix and to what the recipe does not list.
    """
    if scenario not in SCENARIOS:
        raise ValueError(f"scenario must be one of {SCENARIOS}")
    error_floor = error_floor or observables['_meta'].get('error_floor', DEFAULT_ERROR_FLOOR)
    samples = []
    for name, record in observables['samples'].items():
        if name in exclude:
            continue
        sample = _lab_sample(name, record, observables, error_floor, floor_scale, tuple(roles))
        if sample is None:
            continue
        change = (overrides or {}).get(name, {})
        if 'hold_h' in change:
            sample['pre_segments'] = [(T, 3600.0 * change['hold_h']) for T, _ in sample['pre_segments']]
        if 'heated_duration_h' in change:
            # an episode never runs into the next one: it is held until the next one starts at the latest
            starts = [start for start, _, _ in sample['heated']][1:] + [np.inf]
            sample['heated'] = [(start, T, min(3600.0 * change['heated_duration_h'], following - start))
                                for (start, T, _), following in zip(sample['heated'], starts)]
        for species, added in change.get('added_M', {}).items():
            sample['c0_M'][species] = sample['c0_M'].get(species, 0.0) + added
        if not sample['pre_segments']:                       # the age is unknown
            if scenario != 'free':
                sample.update(age_mode='fixed', age_h=SCENARIO_AGES_H[scenario])
            elif name == FREE_AGE_PARAMETER_SAMPLE:
                sample.update(age_mode='parameter', age_h=None)
            elif not sample['heated']:
                sample.update(age_mode='profile', age_h=None)
            else:
                sample.update(age_mode='fixed', age_h=FREE_AGE_BOUNDS_H[0])
        samples.append(sample)
    samples += [_paper_sample(obs_id) for obs_id in paper if obs_id not in exclude]
    return samples


def build_sample_history(sample: dict, age_h: float = None) -> tuple:
    """(segments ((T_K, duration_s), ...) from mixing, spectrum times since mixing [s] of each block) of one sample."""
    if sample['kind'] == 'paper':
        return [(sample['T_K'], sample['t_s'])], [np.array([sample['t_s']])]
    segments = list(sample['pre_segments'])
    if not segments:
        age_h = sample['age_h'] if age_h is None else age_h
        if age_h is None:
            raise ValueError(f"{sample['name']}: an age is needed (free scenario)")
        segments.append((sample['T_rt_K'], 3600.0 * age_h))
    t_zero = sum(d for _, d in segments)                     # the first spectrum, on the clock that starts at mixing
    clock = 0.0
    for start, T_K, duration in sample['heated']:
        if start < clock - 1e-6:
            raise ValueError(f"{sample['name']}: heating episodes overlap")
        if start > clock:
            segments.append((sample['T_rt_K'], start - clock))
        segments.append((T_K, duration))
        clock = start + duration
    last = max(float(block['t_spec_s'].max()) for block in sample['blocks'])
    if last > clock:
        segments.append((sample['T_rt_K'], last - clock))
    return segments, [t_zero + block['t_spec_s'] for block in sample['blocks']]


# ------------------------------------------------------------------------------
# Predictions
# ------------------------------------------------------------------------------
def _shares_from_state(C_M: np.ndarray, idx: dict, nucleus: str, windows: tuple = None) -> np.ndarray:
    """Shares (times × windows) from concentrations (species × times); negative concentrations count as zero.

    windows: the windows the shares are taken over (default: every window declared for the nucleus).
    """
    C = np.maximum(C_M, 0.0)
    windows = tuple(SHARE_SPECIES[nucleus]) if windows is None else windows
    amounts = np.array([sum(n * C[idx[sp]] for sp, n in SHARE_SPECIES[nucleus][w].items()) for w in windows])
    return (amounts / amounts.sum(axis=0)).T


def _block_shares(sample: dict, C_M: np.ndarray, idx: dict) -> list:
    """Shares of each block of a sample from concentrations at the spectrum times of its blocks, one after another."""
    out, start = [], 0
    for block in sample['blocks']:
        stop = start + len(block['t_spec_s'])
        out.append(_shares_from_state(C_M[:, start:stop], idx, block['nucleus'], block['windows']))
        start = stop
    return out


def _paper_value(sample: dict, C_end: np.ndarray, idx: dict) -> float:
    from kinetics.reactor.validation import OBSERVABLES
    c0 = {sp: sample['c0_M'].get(sp, 0.0) for sp in idx}
    return float(OBSERVABLES[sample['observable']][1](np.maximum(C_end, 0.0), idx, c0))


def _initial_composition(sample: dict, water_fraction: float, added: dict = None) -> dict:
    """Composition at mixing: the recipe, the available share of its water, and `added` {species: mol/L} on top
    (the amounts of the species the recipe does not fix)."""
    c0 = dict(sample['c0_M'])
    if sample['wet'] and water_fraction != 1.0:
        c0['H2O'] = water_fraction * c0['H2O']
    for species, amount in (added or {}).items():
        c0[species] = c0.get(species, 0.0) + max(float(amount), 0.0)
    return c0


def calculate_predicted_shares(sample: dict, model, *, age_h: float = None, water_fraction: float = 1.0,
                               added_M: dict = None, solver: dict = None):
    """Shares the model predicts for one lab sample, one array (spectra × windows) per block, or the value of
    a paper observable. The blocks of a sample are read off one simulation.

    added_M: {species: mol/L} present at mixing in this sample on top of its recipe.
    Returns None when the simulation fails or runs out of its call budget.
    """
    solver = solver or SEARCH_SOLVER
    segments, times = build_sample_history(sample, age_h)
    sim = simulate_history(_initial_composition(sample, water_fraction, added_M), segments, np.concatenate(times),
                           model=model, **solver)
    if not sim['success']:
        return None
    if sample['kind'] == 'paper':
        return _paper_value(sample, sim['C_M'][:, -1], sim['idx'])
    return _block_shares(sample, sim['C_M'], sim['idx'])


def _whiten(block: dict, predicted: np.ndarray) -> np.ndarray:
    """Residuals of one block with unit variance (independent values), scaled for the sum constraint."""
    K = len(block['windows'])
    scale = np.sqrt((K - 1.0) / K)
    diff = predicted - block['measured']
    if block['replicate']:
        return scale * np.concatenate([solve_triangular(block['cholesky'][k], diff[:, k], lower=True) for k in range(K)])
    return scale * (diff / block['sigma_ind']).ravel()


def _whiten_sample(sample: dict, predicted: list) -> np.ndarray:
    """Residuals of one lab sample: those of its blocks, one after another."""
    return np.concatenate([_whiten(block, p) for block, p in zip(sample['blocks'], predicted)])


def _paper_residual(sample: dict, value: float) -> float:
    excess = 0.0
    if sample['high'] is not None:
        excess = max(excess, value - sample['high'])
    if sample['low'] is not None:
        excess = max(excess, sample['low'] - value)
    return excess / sample['softness']


def _count(sample: dict) -> int:
    return 1 if sample['kind'] == 'paper' else sum(block['measured'].size for block in sample['blocks'])


def _profile_age(sample: dict, model, water_fraction: float, solver: dict, added: dict = None):
    """Age in FREE_AGE_BOUNDS_H that minimises the χ² of one unheated sample: one simulation, every age.

    The model is run once to the oldest age; the spectrum times of every candidate age lie on that trajectory.
    The best grid age is refined on shares interpolated linearly in log(age). The blocks of the sample share
    the age: the χ² is that of all of them.
    """
    ages_h = np.geomspace(*FREE_AGE_BOUNDS_H, N_AGE_GRID)
    t_spec = np.concatenate([block['t_spec_s'] for block in sample['blocks']])
    times = (3600.0 * ages_h[:, None] + t_spec[None, :]).ravel()
    sim = simulate_history(_initial_composition(sample, water_fraction, added),
                           [(sample['T_rt_K'], float(times.max()))], times, model=model, **solver)
    if not sim['success']:
        return None, None
    C = sim['C_M'].reshape(sim['C_M'].shape[0], len(ages_h), len(t_spec))
    shares = [_block_shares(sample, C[:, a, :], sim['idx']) for a in range(len(ages_h))]
    chi2 = np.array([np.sum(_whiten_sample(sample, p) ** 2) for p in shares])
    j = int(np.argmin(chi2))
    lo, hi = max(j - 1, 0), min(j + 1, len(ages_h) - 1)
    log_age = np.log(ages_h)

    def interpolate(x):
        i = int(np.clip(np.searchsorted(log_age, x) - 1, 0, len(ages_h) - 2))
        f = (x - log_age[i]) / (log_age[i + 1] - log_age[i])
        return [(1.0 - f) * p + f * q for p, q in zip(shares[i], shares[i + 1])]

    res = minimize_scalar(lambda x: np.sum(_whiten_sample(sample, interpolate(x)) ** 2), bounds=(log_age[lo], log_age[hi]),
                          method='bounded', options={'xatol': 1e-3})
    if res.fun < chi2[j]:
        return float(np.exp(res.x)), interpolate(res.x)
    return float(ages_h[j]), shares[j]


def _evaluate_sample(sample: dict, model, ages: dict, water_fraction: float, solver: dict, added_M: dict = None) -> tuple:
    """(predicted, age used) of one sample; predicted is None on failure."""
    added = (added_M or {}).get(sample['name'])
    if sample['age_mode'] == 'profile':
        age_h = (ages or {}).get(sample['name'])
        if age_h is None:
            age_h, predicted = _profile_age(sample, model, water_fraction, solver, added)
            return predicted, age_h
    elif sample['age_mode'] == 'parameter':
        age_h = ages[sample['name']]
    else:
        age_h = sample['age_h']
    return calculate_predicted_shares(sample, model, age_h=age_h, water_fraction=water_fraction, added_M=added,
                                      solver=solver), age_h


# ------------------------------------------------------------------------------
# Residuals of a model over a sample set
# ------------------------------------------------------------------------------
def calculate_residuals(model, samples: list, *, ages: dict = None, water_fraction: float = 1.0,
                        added_M: dict = None, solver: dict = None, return_details: bool = False):
    """Whitened residuals of the model over the samples, in the order of the sample set; χ² is their sum of squares.

    ages: {sample name: age in hours} for the sample whose age is a fit parameter in the free scenario (and,
    optionally, fixed ages for the samples whose age would otherwise be found inside the evaluation).
    water_fraction: share of the added water that is available, applied to the samples mixed with water.
    added_M: {sample name: {species: mol/L}} present at mixing on top of the recipe (the species a recipe does
    not fix, see 'unknown_M' of a sample). A sample without an entry starts from its recipe alone.
    A sample whose simulation fails contributes PENALTY for each of its values.
    With return_details: also {sample: {'predicted', 'age_h', 'ok'}}, 'predicted' holding one array per block.
    """
    solver = solver or SEARCH_SOLVER
    blocks, details = [], {}
    for sample in samples:
        predicted, age_h = _evaluate_sample(sample, model, ages, water_fraction, solver, added_M)
        if predicted is None:
            blocks.append(np.full(_count(sample), PENALTY))
        elif sample['kind'] == 'paper':
            blocks.append(np.array([_paper_residual(sample, predicted)]))
        else:
            blocks.append(_whiten_sample(sample, predicted))
        details[sample['name']] = {'predicted': predicted, 'age_h': age_h, 'ok': predicted is not None}
    residuals = np.concatenate(blocks)
    return (residuals, details) if return_details else residuals


def find_free_age(model, samples: list, *, water_fraction: float = 1.0, added_M: dict = None, solver: dict = None,
                  n_grid: int = 9) -> dict:
    """{sample: age in hours} for the sample whose age is a fit parameter, chosen to minimise its χ² for this model.

    Each sample's residuals depend on its own age only, so that sample is treated alone: a coarse grid in
    log(age), then a bounded search around the best point. Returns {} when no sample has a free age parameter.
    """
    solver = solver or SEARCH_SOLVER
    ages = {}
    for sample in samples:
        if sample['age_mode'] != 'parameter':
            continue

        def chi2(log_age):
            r = calculate_residuals(model, [sample], ages={sample['name']: float(np.exp(log_age))},
                                    water_fraction=water_fraction, added_M=added_M, solver=solver)
            return float(np.sum(r ** 2))

        grid = np.log(np.geomspace(*FREE_AGE_BOUNDS_H, n_grid))
        values = [chi2(x) for x in grid]
        j = int(np.argmin(values))
        res = minimize_scalar(chi2, bounds=(grid[max(j - 1, 0)], grid[min(j + 1, n_grid - 1)]), method='bounded',
                              options={'xatol': 0.02})
        ages[sample['name']] = float(np.exp(res.x if res.fun < values[j] else grid[j]))
    return ages


def _rows_of_sample(sample: dict, predicted, age_h) -> list:
    name = sample['name']
    if predicted is None:
        return [{'sample': name, 'quantity': 'simulation failed', 'kind': 'failure', 't_h': np.nan,
                 'measured': np.nan, 'predicted': np.nan, 'sigma': np.nan, 'z': PENALTY}]
    if sample['kind'] == 'paper':
        limit = sample['high'] if sample['high'] is not None else sample['low']
        sign = '≤' if sample['high'] is not None else '≥'
        return [{'sample': name, 'quantity': f"{sample['observable']} {sign} {limit:g}", 'kind': 'limit',
                 't_h': sample['t_s'] / 3600.0, 'measured': limit, 'predicted': predicted,
                 'sigma': sample['softness'], 'z': _paper_residual(sample, predicted)}]
    _, times = build_sample_history(sample, age_h)
    rows = []
    for block, shares, t_s in zip(sample['blocks'], predicted, times):
        rows += _rows_of_block(name, block, shares, t_s / 3600.0)
    return rows


def _rows_of_block(name: str, block: dict, predicted: np.ndarray, t_h: np.ndarray) -> list:
    Y, rows = block['measured'], []
    if not block['replicate']:
        for i in range(Y.shape[0]):
            for k, w in enumerate(block['windows']):
                rows.append({'sample': name, 'quantity': w, 'kind': 'share', 't_h': t_h[i], 'measured': Y[i, k],
                             'predicted': predicted[i, k], 'sigma': block['sigma_ind'][i, k],
                             'z': (predicted[i, k] - Y[i, k]) / block['sigma_ind'][i, k]})
        return rows
    weights = 1.0 / block['sigma_ind'] ** 2
    weights = weights / weights.sum(axis=0)
    y_mean, p_mean = (weights * Y).sum(axis=0), (weights * predicted).sum(axis=0)
    sigma_mean = np.sqrt(block['sigma_common'] ** 2 + 1.0 / (1.0 / block['sigma_ind'] ** 2).sum(axis=0))
    for k, w in enumerate(block['windows']):
        rows.append({'sample': name, 'quantity': f'{w} (mean of {Y.shape[0]} spectra)', 'kind': 'mean',
                     't_h': float(t_h.mean()), 'measured': y_mean[k], 'predicted': p_mean[k], 'sigma': sigma_mean[k],
                     'z': (p_mean[k] - y_mean[k]) / sigma_mean[k]})
    for i in range(Y.shape[0]):
        for k, w in enumerate(block['windows']):
            measured, model_drift = Y[i, k] - y_mean[k], predicted[i, k] - p_mean[k]
            rows.append({'sample': name, 'quantity': f'{w} drift, spectrum {i + 1}', 'kind': 'drift', 't_h': t_h[i],
                         'measured': measured, 'predicted': model_drift, 'sigma': block['sigma_ind'][i, k],
                         'z': (model_drift - measured) / block['sigma_ind'][i, k]})
    return rows


def tabulate_residuals(model, samples: list, *, ages: dict = None, water_fraction: float = 1.0,
                       added_M: dict = None, solver: dict = None) -> pd.DataFrame:
    """Standardised residuals z = (predicted − measured) / σ, one row per measured quantity.

    kind 'share': one window of one spectrum. For a block of repeated spectra: 'mean' (its composition,
    with the common error) and 'drift' (each spectrum against that mean, with its independent error).
    The window names say which nucleus a row belongs to (SHARE_SPECIES).
    kind 'limit': a one-sided paper observation (z = 0 inside the limit). 't_h' is the time since mixing.
    The fit rule reads this table: a structure fits when no |z| exceeds 3.
    """
    solver = solver or REPORT_SOLVER
    _, details = calculate_residuals(model, samples, ages=ages, water_fraction=water_fraction, added_M=added_M,
                                     solver=solver, return_details=True)
    rows = []
    for sample in samples:
        d = details[sample['name']]
        for row in _rows_of_sample(sample, d['predicted'], d['age_h']):
            rows.append({**row, 'age_h': d['age_h'] if d['age_h'] is not None else np.nan})
    return pd.DataFrame(rows, columns=['sample', 'quantity', 'kind', 't_h', 'age_h', 'measured', 'predicted', 'sigma', 'z'])


def summarize_residuals(residuals: np.ndarray, table: pd.DataFrame = None) -> dict:
    """χ², number of residuals and (with the table) the largest standardised residual and where it sits."""
    out = {'chi2': float(np.sum(residuals ** 2)), 'n_residuals': int(residuals.size)}
    if table is not None and len(table):
        worst = table.loc[table['z'].abs().idxmax()]
        out.update(max_abs_z=float(abs(worst['z'])), worst=f"{worst['sample']}: {worst['quantity']}")
    return out


def copy_sample_set(samples: list) -> list:
    """Independent copy of a sample set (e.g. to replace the measured shares by synthetic ones)."""
    return deepcopy(samples)
