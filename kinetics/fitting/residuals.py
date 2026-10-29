"""Block 14D — Residuals: the lab observables against a model, under one reading of the unknown times.

A sample is one tube: a composition at mixing, a temperature history and one or more spectra, each a set of
area shares with their uncertainties. The spectra of one nucleus form a block; a tube recorded on two nuclei
(³¹P and the ¹³C Si–CH3 region) carries two blocks, which share the composition, the history and the age and
are read off one simulation. The time from mixing to a sample's first spectrum (its age) is not
recorded, so every result is computed under a scenario (a reading): every unknown age is 10 minutes, 1 hour or
8 hours ('10 min', '1 h', '8 h'; eighth fit, M1 tube by tube). The other readings of the earlier fits (fixed ages
of 0 min to 7 days, and ages found by the fit) are in the branch M5-full.

A recipe can hold a species whose amount at mixing is not known (the water of the tube mixed without water):
the sample lists it under 'unknown_M' with its range, a fit takes it as a parameter, and every function that
simulates a sample takes the amounts as added_M = {sample: {species: mol/L}}, on top of the recipe.

Error model of a share: noise ⊕ baseline spread ⊕ a declared floor. Repeated room-temperature spectra of one
block (the 2 % water sample) are one composition and n − 1 measurements of "no drift": the baseline spread and
the floor are common to them, the noise (times the Birge ratio of their scatter) is independent.

The K shares of a spectrum sum to one, so their K residuals are scaled by sqrt((K − 1)/K): that is the exact
χ² of K − 1 independent values when the K errors are equal.

    build_sample_set            observables + scenario → the samples a fit compares with
    set_sample_overrides        changes to the samples that hold for a whole fit run (fit_config.json)
    merge_sample_overrides      the changes of a run with those of one call on top
    build_sample_history        temperature segments from mixing and the spectrum times of one sample
    calculate_predicted_shares  the shares a model predicts for one sample
    calculate_residuals         whitened residual vector of a model over a sample set (χ² = Σ r²)
    tabulate_residuals          the standardised residuals, one row per measured quantity

Source: Y. Alcaraz Galván
"""

from copy import deepcopy

import numpy as np
import pandas as pd
from scipy.linalg import cholesky, solve_triangular

from kinetics.constants import ZERO_CELSIUS_K
from kinetics.data.observables import DEFAULT_ERROR_FLOOR, RT_MAX_C, calculate_replicate_scatter, list_spectrum_sets
from kinetics.reactor.protocol import simulate_history

SCENARIO_AGES_H = {'10 min': 10.0 / 60.0, '1 h': 1.0, '8 h': 8.0}
SCENARIOS = tuple(SCENARIO_AGES_H)
DEFAULT_SCENARIOS = SCENARIOS              # the readings of a run whose fit_config.json lists none

# Species behind each window: {nucleus: {window: {species: nuclei of that window per molecule}}}
SHARE_SPECIES = {
    '31P': {'TMSPA': {'TMSPA': 1}, 'BMSPA': {'BMSPA': 1}, 'MMSPA': {'MMSPA': 1}, 'H3PO4': {'H3PO4': 1}},
    # Si–CH3 carbons, three per Si. The TMSOEG window also holds the second ring-opening product (notebook 03 §2.1).
    # 'P-silyl': the silyl groups still on phosphate, a window of the tubes that hold phosphate only
    '13C': {'TMSOH': {'TMSOH': 1}, 'HMDSO': {'HMDSO': 2}, 'TMSOEG': {'TMSOEG': 1, 'TMSOdiEG': 1},
            'P-silyl': {'TMSPA': 3, 'BMSPA': 2, 'MMSPA': 1}},
}
REPLICATE_MIN_SPECTRA = 3
PENALTY = 1.0e3               # residual given to every value of a sample whose simulation failed
# The absolute tolerance has to sit far below a nanomolar: some fits explain a sample by an induction time that
# is seeded by traces of that size, and with atol = 1e-10 M their χ² came out wrong by a factor of three
SEARCH_SOLVER = {'method': 'BDF', 'rtol': 1e-6, 'atol': 1e-13, 'max_rhs_calls': 30000}
REPORT_SOLVER = {'method': 'BDF', 'rtol': 1e-8, 'atol': 1e-14, 'max_rhs_calls': 300000}


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
        'sigma_nmr': np.hypot(noise, baseline),               # what the spectrum itself gives: noise and baseline, no floor
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
        'name': name, 'kind': 'lab', 'c0_M': c0,
        'T_rt_K': record['T_rt_C'] + ZERO_CELSIUS_K,
        'pre_segments': [(h['T_C'] + ZERO_CELSIUS_K, 3600.0 * h['duration_h']) for h in record['pre_history']],
        'heated': sorted((3600.0 * h['start_h'], h['T_C'] + ZERO_CELSIUS_K, 3600.0 * h['duration_h'])
                         for h in record['heated']),
        'blocks': blocks, 'age_h': None,
        'unknown_M': {sp: (float(low), float(high)) for sp, (low, high) in record.get('unknown_M', {}).items()},
    }


def _hold_heated(sample: dict, hold_h: float) -> list:
    """The heating episodes of a sample as holds of hold_h hours, each over before the next room-temperature spectrum.

    An episode that would run past that spectrum starts earlier, so that it ends at the spectrum; one never starts
    before the previous ends, and is cut where the next recorded one starts. Starts are on the clock of the first spectrum.
    """
    spectra = sorted({float(t) for block in sample['blocks'] for t in block['t_spec_s']})
    starts = [start for start, _, _ in sample['heated']]
    hold, clock, episodes = 3600.0 * hold_h, 0.0, []
    for k, (start, T, _) in enumerate(sample['heated']):
        following = starts[k + 1] if k + 1 < len(starts) else np.inf
        next_spectrum = min([t for t in spectra if t > start], default=np.inf)
        if start + hold > next_spectrum:
            start = next_spectrum - hold
        start = max(start, clock)
        duration = min(hold, following - start, next_spectrum - start)
        episodes.append((start, T, duration))
        clock = start + duration
    return episodes


def _cool_before_spectra(sample: dict, cooling_h: float) -> list:
    """The heating episodes of a sample, each moved earlier where needed so that it ends at least cooling_h hours
    before the next spectrum; a moved one never starts before the previous one ends (it is shortened instead)."""
    spectra = sorted({float(t) for block in sample['blocks'] for t in block['t_spec_s']})
    cooling, clock, episodes = 3600.0 * cooling_h, -np.inf, []
    for start, T, duration in sample['heated']:
        next_spectrum = min([t for t in spectra if t > start], default=np.inf)
        stop = min(start + duration, next_spectrum - cooling)
        start = max(stop - duration, clock)
        episodes.append((start, T, stop - start))
        clock = stop
    return episodes


def _episodes(rows) -> list:
    """[(start_h, T_C, duration_h), ...] on the clock of the first spectrum → [(start_s, T_K, duration_s), ...] in time order."""
    return sorted((3600.0 * float(start), float(T_C) + ZERO_CELSIUS_K, 3600.0 * float(hours)) for start, T_C, hours in rows)


# Changes to the samples that hold for a whole fit run (fit_config.json 'sample_overrides', set by
# kinetics.fitting.config): every sample set starts from them, and the overrides of a call go on top
_ACTIVE_OVERRIDES = [{}]


def get_sample_overrides() -> dict:
    return _ACTIVE_OVERRIDES[0]


def set_sample_overrides(overrides: dict = None) -> dict:
    """Make `overrides` ({sample: {key: value}}, keys as in build_sample_set) the changes every sample set of the
    run starts from (None: none) and return them."""
    _ACTIVE_OVERRIDES[0] = deepcopy(overrides or {})
    return _ACTIVE_OVERRIDES[0]


def merge_sample_overrides(base: dict, extra: dict = None) -> dict:
    """`base` with `extra` on top, key by key within a sample ('added_M' species by species)."""
    out = deepcopy(base or {})
    for name, change in (extra or {}).items():
        entry = out.setdefault(name, {})
        for key, value in change.items():
            entry[key] = {**entry.get(key, {}), **value} if key == 'added_M' else deepcopy(value)
    return out


def build_sample_set(observables: dict, scenario: str, *, roles=('fit',), exclude=(), error_floor: dict = None, floor_scale: float = 1.0, overrides: dict = None) -> list:
    """The samples a model is compared with, under one scenario of the unknown ages.

    roles: which blocks of the lab samples enter ('fit', 'hold_out', 'check'); a sample enters with the blocks
    that have one of these roles. exclude: sample names left out.
    floor_scale multiplies the declared error floor (0 removes it).
    overrides: {sample: {'hold_h': hours of its known pre-history, 'pre_history': [(T_C or None for room
    temperature, hours), ...] from mixing to the first spectrum (the age is then known), 'heated_episodes':
    [(start_h, T_C, hours), ...] on the clock of the first spectrum, replacing the recorded heating (None: the
    recorded one), 'added_episodes': the same, added to it (a negative start lies before the first spectrum and
    is cut at mixing), 'cooling_h': hours at room temperature at least before every spectrum (an episode that
    ends later is moved earlier, see _cool_before_spectra), 'heated_duration_h': hours of every heating episode (at most until the next one starts),
    'heated_hold_h': hours of every heating episode as a hold that is over before the next room-temperature
    spectrum (see _hold_heated), 'added_M': {species: mol/L present at mixing on top of the recipe},
    'unknown_M': {species: (low, high)} replacing the amounts at mixing that are fit parameters ({}: none; sixth
    fit: the water of 'TMSPa alone' declared through 'added_M' instead)}}, for the sensitivity to the times the
    files do not fix and to what the recipe does not list. They go on top of the run's own changes
    (set_sample_overrides).
    """
    if scenario not in SCENARIOS:
        raise ValueError(f"scenario must be one of {SCENARIOS}")
    error_floor = error_floor or observables['_meta'].get('error_floor', DEFAULT_ERROR_FLOOR)
    overrides = merge_sample_overrides(get_sample_overrides(), overrides)
    samples = []
    for name, record in observables['samples'].items():
        if name in exclude:
            continue
        sample = _lab_sample(name, record, observables, error_floor, floor_scale, tuple(roles))
        if sample is None:
            continue
        change = overrides.get(name, {})
        if change.get('heated_episodes') is not None:
            sample['heated'] = _episodes(change['heated_episodes'])
        if 'pre_history' in change:
            sample['pre_segments'] = [(sample['T_rt_K'] if T_C is None else T_C + ZERO_CELSIUS_K, 3600.0 * hours)
                                      for T_C, hours in change['pre_history']]
        if 'hold_h' in change:
            sample['pre_segments'] = [(T, 3600.0 * change['hold_h']) for T, _ in sample['pre_segments']]
        if 'heated_duration_h' in change:
            # an episode never runs into the next one: it is held until the next one starts at the latest
            starts = [start for start, _, _ in sample['heated']][1:] + [np.inf]
            sample['heated'] = [(start, T, min(3600.0 * change['heated_duration_h'], following - start))
                                for (start, T, _), following in zip(sample['heated'], starts)]
        if 'heated_hold_h' in change:
            sample['heated'] = _hold_heated(sample, change['heated_hold_h'])
        if change.get('added_episodes'):
            sample['heated'] = sorted(sample['heated'] + _episodes(change['added_episodes']))
        if change.get('cooling_h') is not None:
            sample['heated'] = _cool_before_spectra(sample, change['cooling_h'])
        for species, added in change.get('added_M', {}).items():
            sample['c0_M'][species] = sample['c0_M'].get(species, 0.0) + added
        if 'unknown_M' in change:
            sample['unknown_M'] = {sp: (float(low), float(high)) for sp, (low, high) in (change['unknown_M'] or {}).items()}
        if not sample['pre_segments']:                       # the age is unknown: that of the reading
            sample['age_h'] = SCENARIO_AGES_H[scenario]
        samples.append(sample)
    return samples


def build_sample_history(sample: dict, age_h: float = None) -> tuple:
    """(segments ((T_K, duration_s), ...) from mixing, spectrum times since mixing [s] of each block) of one sample.

    age_h: hours from mixing to the first spectrum of a sample whose history is not known (default: that of the reading).
    """
    segments = list(sample['pre_segments'])
    if not segments:
        age_h = sample['age_h'] if age_h is None else age_h
        segments.append((sample['T_rt_K'], 3600.0 * age_h))
    t_zero = sum(d for _, d in segments)                     # the first spectrum, on the clock that starts at mixing
    early = [episode for episode in sample['heated'] if episode[0] < 0.0]
    if early:                                                # heating before the first spectrum, inside the age
        if sample['pre_segments']:
            raise ValueError(f"{sample['name']}: heating before the first spectrum needs an unknown age, not a known history")
        segments, clock = [], -t_zero
        for start, T_K, duration in early:
            stop = start + duration
            if stop > 1e-6:
                raise ValueError(f"{sample['name']}: a heating episode before the first spectrum must end before it")
            start = max(start, -t_zero)                      # an episode longer than the age is cut at mixing
            if stop <= start:
                continue
            if start < clock - 1e-6:
                raise ValueError(f"{sample['name']}: heating episodes overlap")
            if start > clock:
                segments.append((sample['T_rt_K'], start - clock))
            segments.append((T_K, stop - start))
            clock = stop
        if clock < 0.0:
            segments.append((sample['T_rt_K'], -clock))
    clock = 0.0
    for start, T_K, duration in sample['heated']:
        if start < 0.0:
            continue
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


def _initial_composition(sample: dict, added: dict = None) -> dict:
    """Composition at mixing: the recipe, and `added` {species: mol/L} on top (the amounts of the species the
    recipe does not fix)."""
    c0 = dict(sample['c0_M'])
    for species, amount in (added or {}).items():
        c0[species] = c0.get(species, 0.0) + max(float(amount), 0.0)
    return c0


def calculate_predicted_shares(sample: dict, model, *, age_h: float = None, added_M: dict = None, solver: dict = None):
    """Shares the model predicts for one lab sample, one array (spectra × windows) per block. The blocks of a
    sample are read off one simulation.

    added_M: {species: mol/L} present at mixing in this sample on top of its recipe.
    Returns None when the simulation fails or runs out of its call budget.
    """
    solver = solver or SEARCH_SOLVER
    segments, times = build_sample_history(sample, age_h)
    sim = simulate_history(_initial_composition(sample, added_M), segments, np.concatenate(times), model=model, **solver)
    if not sim['success']:
        return None
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


def _count(sample: dict) -> int:
    return sum(block['measured'].size for block in sample['blocks'])


def _evaluate_sample(sample: dict, model, ages: dict, solver: dict, added_M: dict = None) -> tuple:
    """(predicted, age used) of one sample; predicted is None on failure."""
    age_h = (ages or {}).get(sample['name'], sample['age_h'])
    return calculate_predicted_shares(sample, model, age_h=age_h, added_M=(added_M or {}).get(sample['name']),
                                      solver=solver), age_h


# ------------------------------------------------------------------------------
# Residuals of a model over a sample set
# ------------------------------------------------------------------------------
def calculate_residuals(model, samples: list, *, ages: dict = None, added_M: dict = None, solver: dict = None,
                        return_details: bool = False):
    """Whitened residuals of the model over the samples, in the order of the sample set; χ² is their sum of squares.

    ages: {sample name: hours from mixing to its first spectrum} in place of the age of the reading.
    added_M: {sample name: {species: mol/L}} present at mixing on top of the recipe (the species a recipe does
    not fix, see 'unknown_M' of a sample). A sample without an entry starts from its recipe alone.
    A sample whose simulation fails contributes PENALTY for each of its values.
    With return_details: also {sample: {'predicted', 'age_h', 'ok'}}, 'predicted' holding one array per block.
    """
    solver = solver or SEARCH_SOLVER
    blocks, details = [], {}
    for sample in samples:
        predicted, age_h = _evaluate_sample(sample, model, ages, solver, added_M)
        if predicted is None:
            blocks.append(np.full(_count(sample), PENALTY))
        else:
            blocks.append(_whiten_sample(sample, predicted))
        details[sample['name']] = {'predicted': predicted, 'age_h': age_h, 'ok': predicted is not None}
    residuals = np.concatenate(blocks)
    return (residuals, details) if return_details else residuals


def _rows_of_sample(sample: dict, predicted, age_h) -> list:
    name = sample['name']
    if predicted is None:
        return [{'sample': name, 'quantity': 'simulation failed', 'kind': 'failure', 't_h': np.nan,
                 'measured': np.nan, 'predicted': np.nan, 'sigma': np.nan, 'z': PENALTY}]
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


def tabulate_residuals(model, samples: list, *, ages: dict = None, added_M: dict = None,
                       solver: dict = None) -> pd.DataFrame:
    """Standardised residuals z = (predicted − measured) / σ, one row per measured quantity.

    kind 'share': one window of one spectrum. For a block of repeated spectra: 'mean' (its composition,
    with the common error) and 'drift' (each spectrum against that mean, with its independent error).
    The window names say which nucleus a row belongs to (SHARE_SPECIES).
    't_h' is the time since mixing.
    The fit rule reads this table: a structure fits when no |z| exceeds 3.
    """
    solver = solver or REPORT_SOLVER
    _, details = calculate_residuals(model, samples, ages=ages, added_M=added_M, solver=solver, return_details=True)
    rows = []
    for sample in samples:
        d = details[sample['name']]
        for row in _rows_of_sample(sample, d['predicted'], d['age_h']):
            rows.append({**row, 'age_h': d['age_h'] if d['age_h'] is not None else np.nan})
    return pd.DataFrame(rows, columns=['sample', 'quantity', 'kind', 't_h', 'age_h', 'measured', 'predicted', 'sigma', 'z'])


