"""Checks of the kinetic models against the observations of Gogoi et al. 2024.

Control experiments (data/experimental_gogoi2024.json): each is re-simulated with the batch reactor at its
composition, temperature and time, and the predicted observable is compared with the observed window.

Water series: ³¹P of TMSPA in EC/DEC with 0.5–5 vol% water at an unstated time after mixing, and the 2 vol%
sample heated in steps to 80 °C. The windows are our reading of a qualitative description. Because the time
is unknown, a model is tested on whether one common time satisfies every sample.
"""

import numpy as np
import pandas as pd

from kinetics.constants import ZERO_CELSIUS_K
from kinetics.data.experimental import get_water_series, load_experimental_data
from kinetics.data.network import NETWORK_SPECIES
from kinetics.microkinetics.models import get_model
from kinetics.reactor.batch import simulate_batch_reactor
from kinetics.reactor.protocol import DEFAULT_DENSITY_G_ML, Stage, simulate_protocol
from kinetics.data.snapshot import calculate_molar_mass

# Observable name → (label, function of final concentrations C [M], species index, initial concentrations)
OBSERVABLES = {
    'ring_opened_fraction': ('TMSOH ring-opened into TMSOEG + TMSOdiEG (fraction)',
                             lambda C, i, c0: (C[i['TMSOEG']] + C[i['TMSOdiEG']]) / c0['TMSOH']),
    'hmdso_si_fraction': ('Si of TMSOH ending in HMDSO (fraction)',
                          lambda C, i, c0: 2.0 * C[i['HMDSO']] / c0['TMSOH']),
    'tmspa_conversion': ('TMSPA converted (fraction)',
                         lambda C, i, c0: 1.0 - C[i['TMSPA']] / c0['TMSPA']),
}

# Reaction whose barrier each observable probes most directly
OBSERVABLE_REACTIONS = {'ring_opened_fraction': 'R8', 'hmdso_si_fraction': 'R4', 'tmspa_conversion': 'R5'}

PHOSPHATES = ('TMSPA', 'BMSPA', 'MMSPA', 'H3PO4')


def simulate_control_experiment(experiment: dict,   # one entry of 'control_experiments'
                                model,              # the microkinetic model to use
                                network: dict = None,       # network to simulate (default: NETWORK)
                                species_db: dict = None,    # species database to simulate (default: NETWORK)
                                method: str = 'Radau') -> float: # ODE solver. 'Radau' is the default; 'BDF' gives the same result and is faster in scans.
    
    """
    [Checked - YA]
    Simulate one control experiment and return the value of its observable at the end.

    - experiment: one entry of 'control_experiments' (initial mixture, T, time, observable).
    - method: ODE solver. 'Radau' is the default; 'BDF' gives the same result and is faster in scans.
    """

    c0_M = dict.fromkeys(NETWORK_SPECIES, 0.0)
    c0_M.update(experiment['c0_M'])
    sim = simulate_batch_reactor(c0_M, T_K=experiment['T_K'], t_end_s=experiment['t_s'], model=model,
                                 network=network, species_db=species_db, n_points=40, t_start_s=1.0, method=method)
    return float(OBSERVABLES[experiment['observable']][1](sim['C_M'][:, -1], sim['idx'], c0_M))


def is_within_window(value: float, low: float, high: float) -> bool:
    """value inside [low, high]; None or NaN means an open side."""
    low_ok = low is None or np.isnan(low) or value >= low
    high_ok = high is None or np.isnan(high) or value <= high
    return bool(low_ok and high_ok)


def evaluate_control_experiments(models, network: dict = None, 
                                 species_db: dict = None,
                                 experiments: list = None, 
                                 method: str = 'Radau') -> pd.DataFrame:
    """
    # [Checked - YA]
    Simulate every control experiment with every model and compare with the observed window.

    experiments: entries in the format of 'control_experiments' (default: the file), e.g. with a changed
    window or time to test an assumption. Returns one row per (experiment, model).
    """
    rows = []
    for exp in experiments if experiments is not None else load_experimental_data()['control_experiments']:
        label = OBSERVABLES[exp['observable']][0]
        low = np.nan if exp['low'] is None else exp['low']
        high = np.nan if exp['high'] is None else exp['high']
        for m in models:
            spec = get_model(m)
            value = simulate_control_experiment(exp, spec, network, species_db, method)
            rows.append({
                'experiment': exp['id'], 'label': exp['label'], 'observable': exp['observable'],
                'observable_label': label, 'model': spec.name,
                'predicted': value, 'low': low, 'high': high, 'status': exp['status'],
                'observation': exp['observation'], 'consistent': is_within_window(value, low, high),
            })
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------
# Water series (TMSPA + H2O in EC/DEC)
# ------------------------------------------------------------------------------
def calculate_water_series_c0(h2o_vol_pct: float, medium: dict = None) -> dict:
    """Initial concentrations [M] of a water-series sample: the H2O stock diluted by the TMSPA addition."""
    medium = medium or get_water_series()['medium']
    dilution = 1.0 - medium['tmspa_vol_frac']
    c0 = dict.fromkeys(NETWORK_SPECIES, 0.0)
    c0['H2O'] = h2o_vol_pct / 100.0 * 1000.0 * DEFAULT_DENSITY_G_ML['H2O'] / calculate_molar_mass('H2O') * dilution
    c0['TMSPA'] = medium['tmspa_vol_frac'] * 1000.0 * DEFAULT_DENSITY_G_ML['TMSPA'] / calculate_molar_mass('TMSPA')
    c0['EC'] = medium['EC_M']
    return c0


def calculate_phosphate_fractions(C_M: np.ndarray, idx: dict) -> dict:
    """Fraction of all P in each phosphate (arrays over the last axis of C_M)."""
    total = sum(C_M[idx[sp]] for sp in PHOSPHATES)
    return {sp: C_M[idx[sp]] / total for sp in PHOSPHATES}


def _window_value(fractions: dict, key: str):
    return sum(fractions[sp] for sp in key.split('+'))


def evaluate_water_series(models, times_h=None, network: dict = None, species_db: dict = None,
                          method: str = 'Radau') -> pd.DataFrame:
    """P fractions of every water-series sample on a common time grid, and the windows each time meets.

    Returns one row per (model, time): the fractions per sample, one 'W… ok' column per sample, and
    'all ok'. A model is consistent with the series only if some common time has 'all ok'.
    """
    series = get_water_series()
    times_h = np.logspace(np.log10(1.0 / 60.0), np.log10(168.0), 120) if times_h is None else np.asarray(times_h)
    T_K = series['medium']['T_K']
    frames = []
    for m in models:
        spec = get_model(m)
        table = pd.DataFrame({'model': spec.name, 't_h': times_h})
        all_ok = np.ones(len(times_h), dtype=bool)
        for sample in series['samples']:
            c0 = calculate_water_series_c0(sample['h2o_vol_pct'], series['medium'])
            sim = simulate_batch_reactor(c0, T_K=T_K, t_end_s=times_h[-1] * 3600.0, model=spec, network=network,
                                         species_db=species_db, n_points=300, t_start_s=1.0, method=method)
            fractions = {sp: np.interp(np.log(times_h), np.log(sim['t_h']), f)
                         for sp, f in calculate_phosphate_fractions(sim['C_M'], sim['idx']).items()}
            ok = np.ones(len(times_h), dtype=bool)
            for key, (low, high) in sample['windows'].items():
                value = _window_value(fractions, key)
                ok &= np.array([is_within_window(v, low, high) for v in value])
            for sp in PHOSPHATES:
                table[f"{sample['id']} {sp}"] = fractions[sp]
            table[f"{sample['id']} ok"] = ok
            all_ok &= ok
        table['all ok'] = all_ok
        frames.append(table)
    return pd.concat(frames, ignore_index=True)


def evaluate_heating_observation(models, rt_hold_h=(1.0, 6.0, 24.0), network: dict = None,
                                 species_db: dict = None) -> pd.DataFrame:
    """The 2 vol% sample: P fractions after rt_hold_h at RT and after the last heating step.

    'RT ok' applies the W2 windows to the RT spectrum; 'heating ok' requires every P fraction to change by
    at most the observed 'max_change' between the RT spectrum and the end of the 80 °C step.
    """
    series = get_water_series()
    heat = series['heating']
    w2 = next(s for s in series['samples'] if s['h2o_vol_pct'] == heat['h2o_vol_pct'])
    c0 = calculate_water_series_c0(heat['h2o_vol_pct'], series['medium'])
    rt_C = series['medium']['T_K'] - ZERO_CELSIUS_K
    recipe = {'stock': {'H2O': c0['H2O'], 'EC': c0['EC'], 'TMSPA': 0.0},
              'after': {'H2O': c0['H2O'], 'EC': c0['EC'], 'TMSPA': c0['TMSPA']}, 'dilution_factor': 1.0}
    rows = []
    for m in models:
        spec = get_model(m)
        for rt_h in rt_hold_h:
            stages = [Stage('RT', rt_C, rt_h, rt_C)]
            stages += [Stage(f'hold {T:g} °C', T, heat['hold_h'], float(T)) for T in heat['steps_C']]
            sim = simulate_protocol(stages, recipe, model=spec, network=network, species_db=species_db,
                                    injection_stage_idx=0, points_per_stage=40)
            first, last = (calculate_phosphate_fractions(a['C_M'], sim['idx'])
                           for a in (sim['acquisitions'][0], sim['acquisitions'][-1]))
            change = max(abs(last[sp] - first[sp]) for sp in PHOSPHATES)
            rows.append({
                'model': spec.name, 'RT hold (h)': rt_h,
                **{f'RT {sp}': first[sp] for sp in PHOSPHATES},
                **{f'80 °C {sp}': last[sp] for sp in PHOSPHATES},
                'max change': change,
                'RT ok': all(is_within_window(_window_value(first, k), lo, hi) for k, (lo, hi) in w2['windows'].items()),
                'heating ok': change <= heat['max_change'],
            })
    return pd.DataFrame(rows)
