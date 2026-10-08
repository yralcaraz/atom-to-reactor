"""Block 14G — Reporting: a stored fit against the measurements, spectrum by spectrum and in time.

A fit result holds a parameter set and its table of standardised residuals, in which the repeated spectra of a
tube are folded into a mean and its drifts. A reader wants the same comparison in the units of the experiment:
the share of each peak in each spectrum, the concentrations of a tube since it was mixed, and the spectrum
itself. Nothing here is fitted: a stored parameter set is run forward.

    evaluate_fit                     the samples of a fit with the age and the shares its parameter set gives each
    tabulate_fit_shares              measured and predicted share of every window of every spectrum
    simulate_fit_trajectories        concentrations of each tube from mixing, along its temperature history
    tabulate_window_amounts          a trajectory as the amount of nuclei each window of one nucleus counts
    load_measured_spectra            the raw spectra behind the shares: baseline removed, area of the windows = 1
    calculate_parameter_correlation  how the parameters move together near the best fit (local, linear)

Source: Y. Alcaraz Galván
"""

import numpy as np
import pandas as pd

from kinetics.constants import ZERO_CELSIUS_K
from kinetics.data.observables import C_REGION_PPM, C_WINDOWS_PPM, P_WINDOWS_PPM, list_spectrum_sets
from kinetics.fitting.candidates import get_structure
from kinetics.fitting.estimation import FitProblem, calculate_parameter_directions
from kinetics.fitting.residuals import (
    FREE_AGE_BOUNDS_H, N_AGE_GRID, REPORT_SOLVER, SHARE_SPECIES, _block_shares, _initial_composition, _whiten_sample,
    build_sample_history, build_sample_set, calculate_residuals,
)
from kinetics.reactor.protocol import simulate_history

P_REGION_PPM = (-32.0, 6.0)              # region of the ³¹P shares (build_lab_inventory_for_observables)
PARAMETER_KINDS = ('barrier', 'energy', 'recipe', 'age')
N_EXACT_AGES = 240                       # ages tried between the two trial ages of the fit that enclose a stored age


def _find_exact_age(sample: dict, built: dict, age_h: float, solver: dict) -> float:
    """Age near age_h at which the model, run without interpolation, is closest to the spectra of one sample.

    A fit finds the age of an unheated sample on shares interpolated between trial ages (_profile_age). Where
    the composition changes within one step of those ages, the stored age is off by a part of that step; here
    the same χ² is minimised on the trajectory itself, inside the step on either side of the stored age.
    """
    step = (FREE_AGE_BOUNDS_H[1] / FREE_AGE_BOUNDS_H[0]) ** (1.0 / (N_AGE_GRID - 1))
    ages_h = np.unique(np.clip(np.geomspace(age_h / step, age_h * step, N_EXACT_AGES), *FREE_AGE_BOUNDS_H))
    t_spec = np.concatenate([block['t_spec_s'] for block in sample['blocks']])
    times = (3600.0 * ages_h[:, None] + t_spec[None, :]).ravel()
    sim = simulate_history(_initial_composition(sample, built['water_fraction'], built['added_M'].get(sample['name'])),
                           [(sample['T_rt_K'], float(times.max()))], times, model=built['model'], **solver)
    if not sim['success']:
        return age_h
    C = sim['C_M'].reshape(sim['C_M'].shape[0], len(ages_h), len(t_spec))
    chi2 = [np.sum(_whiten_sample(sample, _block_shares(sample, C[:, a, :], sim['idx'])) ** 2) for a in range(len(ages_h))]
    return float(ages_h[int(np.argmin(chi2))])


def evaluate_fit(fit: dict, observables: dict, *, overrides: dict = None, ages: dict = None, solver: dict = None) -> dict:
    """What the parameter set of a stored fit predicts for the samples it was fitted to.

    Returns 'samples' (the sample set of the fit's scenario), 'built' (FitStructure.build of its θ), 'ages'
    {sample: hours from mixing to the first spectrum, None where the history is known} and 'predicted'
    {sample: one array (spectra × windows) per block, or the value of a paper observable}.
    The ages are those of the fit, except for a sample whose age the fit finds inside each evaluation (free
    scenario): that one is found again on the exact trajectory (_find_exact_age).
    overrides: as in build_sample_set (e.g. TMSOH present at mixing).
    ages: {sample: hours} used as they are, e.g. those of another evaluation of the same fit.
    """
    solver = solver or REPORT_SOLVER
    samples = build_sample_set(observables, fit['scenario'], overrides=overrides)
    built = get_structure(fit['structure']).build(fit['theta'])
    if ages is None:
        ages = {**fit.get('ages', {}), **built['ages']}
        for sample in samples:
            if sample['age_mode'] == 'profile' and ages.get(sample['name']) is not None:
                ages[sample['name']] = _find_exact_age(sample, built, ages[sample['name']], solver)
    _, details = calculate_residuals(built['model'], samples, ages=ages, water_fraction=built['water_fraction'],
                                     added_M=built['added_M'], solver=solver, return_details=True)
    return {'samples': samples, 'built': built, 'ages': {name: d['age_h'] for name, d in details.items()},
            'predicted': {name: d['predicted'] for name, d in details.items()}}


def tabulate_fit_shares(evaluation: dict) -> pd.DataFrame:
    """Measured and predicted share of every window of every spectrum of the lab samples, one row each.

    Unlike tabulate_residuals, repeated spectra are not folded into a mean and its drifts. 'sigma' is the error
    of a single spectrum (independent and common parts in quadrature), 't_h' the time since mixing and 'total_M'
    the nuclei the windows of that spectrum count together, so share × total_M is an amount in mol/L.
    """
    rows = []
    for sample in evaluation['samples']:
        name, predicted = sample['name'], evaluation['predicted'][sample['name']]
        if sample['kind'] != 'lab' or predicted is None:
            continue
        age_h = evaluation['ages'][name]
        c0 = _initial_composition(sample, evaluation['built']['water_fraction'], evaluation['built']['added_M'].get(name))
        _, times = build_sample_history(sample, age_h)
        for block, shares, t_s in zip(sample['blocks'], predicted, times):
            species = SHARE_SPECIES[block['nucleus']]
            total = sum(n * c0.get(sp, 0.0) for w in block['windows'] for sp, n in species[w].items())
            sigma = np.hypot(block['sigma_ind'], block['sigma_common'])
            for i in range(block['measured'].shape[0]):
                for k, window in enumerate(block['windows']):
                    rows.append({'sample': name, 'nucleus': block['nucleus'], 'spectrum': i, 't_h': t_s[i] / 3600.0,
                                 'age_h': age_h if age_h is not None else np.nan, 'window': window,
                                 'measured': block['measured'][i, k], 'sigma': sigma[i, k],
                                 'predicted': shares[i, k], 'total_M': total})
    return pd.DataFrame(rows)


def simulate_fit_trajectories(evaluation: dict, *, until_h: float = None, t_first_h: float = 0.1, n_points: int = 300,
                              solver: dict = None) -> pd.DataFrame:
    """Concentrations [M] of every species in each lab sample of a fit, from mixing to its last spectrum.

    Each sample follows its own history: the age the fit gave it, then the heating episodes between its spectra.
    until_h: the sample stays at room temperature after its last spectrum until this time since mixing.
    One row per (sample, time): 't_h' since mixing, 'T_C' and one column per species.
    """
    solver = solver or REPORT_SOLVER
    built, frames = evaluation['built'], []
    for sample in evaluation['samples']:
        if sample['kind'] != 'lab':
            continue
        name = sample['name']
        segments, _ = build_sample_history(sample, evaluation['ages'][name])
        total_s = sum(duration for _, duration in segments)
        if until_h is not None and 3600.0 * until_h > total_s:
            segments.append((sample['T_rt_K'], 3600.0 * until_h - total_s))
            total_s = 3600.0 * until_h
        ends = np.cumsum([duration for _, duration in segments])
        times = np.unique(np.concatenate([np.geomspace(3600.0 * t_first_h, total_s, n_points), ends]))
        sim = simulate_history(_initial_composition(sample, built['water_fraction'], built['added_M'].get(name)),
                               segments, times, model=built['model'], **solver)
        if not sim['success']:
            continue
        T_K = np.array([T for T, _ in segments])[np.minimum(np.searchsorted(ends, times, side='left'), len(segments) - 1)]
        frames.append(pd.DataFrame({'sample': name, 't_h': times / 3600.0, 'T_C': T_K - ZERO_CELSIUS_K,
                                    **{sp: sim['C_M'][i] for sp, i in sim['idx'].items()}}))
    return pd.concat(frames, ignore_index=True)


def tabulate_window_amounts(trajectories: pd.DataFrame, nucleus: str, windows=None) -> pd.DataFrame:
    """The trajectories with one column per NMR window: the nuclei of that window [mol/L], Σ n · C over its species.

    For ³¹P a window is one phosphate species. For ¹³C (Si–CH3) it counts silyl groups: HMDSO carries two.
    """
    species = SHARE_SPECIES[nucleus]
    windows = list(windows) if windows is not None else list(species)
    amounts = {w: sum(n * trajectories[sp].clip(lower=0.0) for sp, n in species[w].items()) for w in windows}
    return pd.DataFrame({'sample': trajectories['sample'], 't_h': trajectories['t_h'], 'T_C': trajectories['T_C'], **amounts})


def load_measured_spectra(observables: dict, root=None, *, baseline_order: int = 3, exclusion_pad_ppm: float = 0.3) -> dict:
    """The raw spectra behind the shares, ready to be drawn against a model.

    {(sample, nucleus, spectrum number): {'ppm', 'intensity', 'windows_ppm', 'axis_offset_ppm', 'noise',
    'temperature_C', 't_since_first_h'}}. Each spectrum is cut to the region its shares were taken from, a
    polynomial baseline fitted outside the windows is removed (one of the three baselines of
    calculate_area_shares) and the intensity is divided by the summed area of the windows, so the area of a
    window is close to its share. Needs the raw spectra (LAB_NMR_DIR).
    """
    from kinetics.data.lab_nmr import _polynomial_baseline, build_lab_nmr_inventory, load_spectrum
    from kinetics.paths import lab_nmr_dir
    inventory = build_lab_nmr_inventory(root or lab_nmr_dir()).drop_duplicates('data_md5').set_index('data_md5')
    out = {}
    for name, record in observables['samples'].items():
        for nucleus, _, spectra in list_spectrum_sets(record):
            more = record.get('more_nuclei', {}).get(nucleus)
            if nucleus == '31P':
                windows, region, offset = P_WINDOWS_PPM, P_REGION_PPM, 0.0
            elif more is not None:
                offset = more['axis_offset_ppm']
                windows, region = more['windows_ppm'], (C_REGION_PPM[0] + offset, C_REGION_PPM[1] + offset)
            else:
                windows, region, offset = C_WINDOWS_PPM, C_REGION_PPM, 0.0
            for number, spectrum in enumerate(spectra):
                ppm, intensity = load_spectrum(inventory.loc[spectrum['data_md5'], 'path'], window_ppm=region)
                inside = (ppm >= region[0]) & (ppm <= region[1])
                x, y = ppm[inside], intensity[inside]
                signal_free = np.ones(len(x), dtype=bool)
                for lo, hi in windows.values():
                    signal_free &= ~((x >= lo - exclusion_pad_ppm) & (x <= hi + exclusion_pad_ppm))
                y = y - _polynomial_baseline(x, y, signal_free, baseline_order)
                step = abs(x[1] - x[0])
                area = sum(float(y[(x >= lo) & (x <= hi)].sum() * step) for lo, hi in windows.values())
                out[(name, nucleus, number)] = {
                    'ppm': x, 'intensity': y / area, 'windows_ppm': {w: tuple(v) for w, v in windows.items()},
                    'axis_offset_ppm': offset, 'noise': float(np.std(y[signal_free]) / area),
                    'temperature_C': spectrum['temperature_C'], 't_since_first_h': spectrum['t_since_first_h']}
    return out


def calculate_parameter_correlation(fit: dict, observables: dict, *, kinds=PARAMETER_KINDS) -> dict:
    """How the fitted parameters of a stored fit move together near its best fit.

    From J = ∂(residuals)/∂(parameters) at the best fit (calculate_parameter_directions): covariance (JᵀJ)⁻¹,
    'correlation' its normalised form and 'sigma_local' the standard error of each parameter. Local and linear:
    it describes the neighbourhood of the best fit, and a parameter the data do not determine shows up as a
    standard error wider than its search box. A parameter in a combination the data do not see at all gets NaN.
    """
    problem = FitProblem(fit['structure'], build_sample_set(observables, fit['scenario']))
    found = calculate_parameter_directions(problem, fit['theta'], kinds=kinds)
    names, sigma = found['parameters'], found['sigma_local'].to_numpy()
    covariance = np.linalg.pinv(found['J'].T @ found['J'])
    spread = np.where(np.isfinite(sigma), np.sqrt(np.diag(covariance)), np.nan)
    return {'correlation': pd.DataFrame(covariance / np.outer(spread, spread), index=names, columns=names),
            'sigma_local': pd.Series(spread, index=names)}
