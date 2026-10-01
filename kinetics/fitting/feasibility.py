"""Block 14A — Feasibility fit: which family barriers the windowed observations allow.

The control experiments of Gogoi et al. 2024 give windows ("< 2 % ring-opened after one week"), not time
series, so they bound parameters instead of fitting them. For each family the intrinsic barrier g is scanned
with the other families fixed at the model's values. Each window edge that the predicted observable crosses
becomes a bound on g, refined by bisection; the feasible interval of a family is the intersection over all
experiments.

An observable depends on a family's g only through g(T_exp), the barrier at the experiment's temperature.
A bound found with ΔS‡ = 0 is therefore a bound on g(T_exp), and project_bounds_to_entropy maps it onto the
(g(T_ref), ΔS‡) plane.
"""

import numpy as np
import pandas as pd

from kinetics.constants import EV_TO_KJ_MOL
from kinetics.data.experimental import load_experimental_data
from kinetics.data.network import NETWORK
from kinetics.microkinetics.models import get_model
from kinetics.microkinetics.parameters import resolve_family_params
from kinetics.microkinetics.rates import KINETIC_MODELS
from kinetics.reactor.validation import OBSERVABLE_REACTIONS, is_within_window, simulate_control_experiment

FAMILIES = ('hydrolysis', 'transfer', 'condensation', 'solvent_attack')
DEFAULT_G_GRID_EV = np.round(np.arange(0.60, 1.701, 0.05), 3)
# BDF gives the same observables as Radau (6 significant figures) and is far faster when fast steps act on
# near-zero species; LSODA is faster still but can stall (e.g. hydrolysis g = 0.80 eV with 2 vol% water)
SCAN_METHOD = 'BDF'


def get_family_parameter(model, family: str, key: str = 'g_eV') -> float:
    """Value of one family parameter in a model (the set's 'default' entry if the family has none)."""
    spec = get_model(model)
    fallback = KINETIC_MODELS[spec.kinetic_model][1]
    return resolve_family_params(family, spec.family_params, fallback)[key]


def _experiments(experiments):
    return experiments if experiments is not None else load_experimental_data()['control_experiments']


def scan_family_barriers(model, 
                         g_grid_eV=DEFAULT_G_GRID_EV,   # grid of g values to scan for each family
                         families=FAMILIES,             # families to scan (default: all)
                         experiments=None,              # control experiments to scan against (default: Gogoi et al. 2024)
                         network: dict = None,          # network to simulate (default: NETWORK)
                         species_db: dict = None) -> pd.DataFrame:  # species database to simulate (default: NETWORK) 
    
    """
    [Checked - YA]
    Find which barrier values g agree with the control experiments.

    - Take one family and try each g of the grid; the other families keep the model's values.
    - For each g, simulate every control experiment and get the predicted observable.
    - Check if the prediction falls inside the measured window [low, high].

    Returns a table with one row per (family, g, experiment); 'consistent' is True if inside the window.
    """

    spec = get_model(model)
    rows = []
    for family in families:
        for g in g_grid_eV:
            variant = spec.with_family_params(family, g_eV=float(g))
            for exp in _experiments(experiments):
                value = simulate_control_experiment(exp, variant, network, species_db, SCAN_METHOD)
                rows.append({'family': family, 'g_eV': float(g), 'experiment': exp['id'],
                             'observable': exp['observable'], 'T_K': exp['T_K'], 'predicted': value,
                             'low': np.nan if exp['low'] is None else exp['low'],
                             'high': np.nan if exp['high'] is None else exp['high'],
                             'status': exp['status'],
                             'consistent': is_within_window(value, exp['low'], exp['high'])})
    return pd.DataFrame(rows)


def _bisect_edge(spec, family, exp, edge_value, g_lo, g_hi, tol_eV, network, species_db):
    """g where the observable crosses edge_value between g_lo and g_hi (sign change assumed)."""
    def residual(g):
        return simulate_control_experiment(exp, spec.with_family_params(family, g_eV=g), network, species_db,
                                           SCAN_METHOD) - edge_value
    r_lo = residual(g_lo)
    while g_hi - g_lo > tol_eV:
        g_mid = 0.5 * (g_lo + g_hi)
        r_mid = residual(g_mid)
        if np.sign(r_mid) == np.sign(r_lo):
            g_lo, r_lo = g_mid, r_mid
        else:
            g_hi = g_mid
    return 0.5 * (g_lo + g_hi)


def find_barrier_bounds(model, families=FAMILIES, *, scan: pd.DataFrame = None, experiments=None,
                        g_grid_eV=DEFAULT_G_GRID_EV, tol_eV: float = 0.002, network: dict = None,
                        species_db: dict = None) -> pd.DataFrame:
    """
    [Checked - YA]
    Find the g values where a prediction crosses an edge of an experiment's window.

    - Uses the scan to spot where the prediction crosses a window edge (low or high).
    - Refines each crossing by bisection (halving the g interval) down to tol_eV.
    - Each crossing is a bound: 'g ≥' (consistent only above it) or 'g ≤' (only below it).

    Returns one row per bound: 
    - family, 
    - experiment, 
    - bound, 
    - g_eV, 
    - T_K, 
    - the reaction the experiment probes,
    - 'acts' (directly, or indirectly through another family's reaction) 
    - and, for direct bounds, the barrier ΔG‡ at the bound (NaN if indirect).
    """
    spec = get_model(model)
    exps = {e['id']: e for e in _experiments(experiments)}
    if scan is None:
        scan = scan_family_barriers(spec, g_grid_eV, families, list(exps.values()), network, species_db)
    net = network or NETWORK
    rows = []
    for (family, exp_id), grp in scan[scan['family'].isin(families)].groupby(['family', 'experiment'], sort=False):
        grp = grp.sort_values('g_eV')
        exp = exps[exp_id]
        g, pred, ok = grp['g_eV'].values, grp['predicted'].values, grp['consistent'].values
        for edge_name in ('low', 'high'):
            edge = exp[edge_name]
            if edge is None:
                continue
            crossings = np.nonzero(np.diff(np.sign(pred - edge)) != 0)[0]
            for i in crossings:
                g_bound = _bisect_edge(spec, family, exp, edge, g[i], g[i + 1], tol_eV, network, species_db)
                side = 'g ≤' if ok[i] and not ok[i + 1] else 'g ≥' if ok[i + 1] and not ok[i] else 'edge'
                rxn_id = OBSERVABLE_REACTIONS[exp['observable']]
                direct = net[rxn_id]['class'] == family
                dG = np.nan
                if direct:
                    rates = spec.with_family_params(family, g_eV=g_bound).calculate_rates(
                        exp['T_K'], network=network, species_db=species_db)
                    dG = float(rates.loc[rxn_id, 'dG_barrier_f_eV'])
                rows.append({'family': family, 'experiment': exp_id, 'window edge': edge_name, 'bound': side,
                             'g_eV': g_bound, 'T_K': exp['T_K'], 'reaction': rxn_id,
                             'acts': 'directly' if direct else 'indirectly',
                             'dG_barrier_at_bound_eV': dG, 'status': exp['status']})
    return pd.DataFrame(rows, columns=['family', 'experiment', 'window edge', 'bound', 'g_eV', 'T_K', 'reaction',
                                       'acts', 'dG_barrier_at_bound_eV', 'status'])


def summarize_feasible_intervals(bounds: pd.DataFrame, scan: pd.DataFrame, model=None,
                                 families=FAMILIES) -> pd.DataFrame:
    """
    [Checked - YA] 
    Per family: the feasible interval of g, the experiments that set each end, and a verdict.

    The interval is limited to the scanned range; 'current g' is the model's value when a model is given.
    """
    rows = {}
    for family in families:
        b = bounds[bounds['family'] == family]
        lower, upper = b[b['bound'] == 'g ≥'], b[b['bound'] == 'g ≤']
        lo = lower.loc[lower['g_eV'].idxmax()] if len(lower) else None
        hi = upper.loc[upper['g_eV'].idxmin()] if len(upper) else None
        grid = scan[scan['family'] == family].groupby('g_eV')['consistent'].all()
        g_min, g_max = grid.index.min(), grid.index.max()
        if not grid.any():
            verdict = 'no feasible value in the scanned range'
        elif lo is not None and hi is not None:
            verdict = 'bounded on both sides'
        elif lo is not None:
            verdict = 'lower bound only'
        elif hi is not None:
            verdict = 'upper bound only'
        else:
            verdict = 'unconstrained'
        row = {
            'g low (eV)': lo['g_eV'] if lo is not None else np.nan,
            'set by (low)': f"{lo['experiment']} ({lo['status']})" if lo is not None else f'scan limit {g_min:.2f}',
            'g high (eV)': hi['g_eV'] if hi is not None else np.nan,
            'set by (high)': f"{hi['experiment']} ({hi['status']})" if hi is not None else f'scan limit {g_max:.2f}',
            'verdict': verdict,
        }
        if model is not None:
            row['current g (eV)'] = get_family_parameter(model, family)
        rows[family] = row
    return pd.DataFrame.from_dict(rows, orient='index')


def project_bounds_to_entropy(bounds: pd.DataFrame, 
                              family: str, 
                              dS_grid_J_mol_K, 
                              T_ref_K: float = 298.15,
                              base_dS_J_mol_K: float = 0.0) -> pd.DataFrame:
    """
    [Checked - YA]
    Bounds of one family as lines g(T_ref) = g_bound + (T_exp − T_ref)·(ΔS‡ − ΔS‡_base) in the (ΔS‡, g) plane.

    Each bound fixes g at its experiment's temperature; g(T) = g(T_ref) − (T − T_ref)·ΔS‡. The bounds were
    found with ΔS‡ = base_dS_J_mol_K. Returns one row per (bound, ΔS‡).
    """
    rows = []
    for _, b in bounds[bounds['family'] == family].iterrows():
        for dS in dS_grid_J_mol_K:
            shift_eV = (b['T_K'] - T_ref_K) * (dS - base_dS_J_mol_K) / 1000.0 / EV_TO_KJ_MOL
            rows.append({'experiment': b['experiment'], 'bound': b['bound'], 'T_K': b['T_K'],
                         'dS_J_mol_K': float(dS), 'g_ref_eV': b['g_eV'] + shift_eV})
    return pd.DataFrame(rows)
