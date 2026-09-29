"""Block 10 — Multi-stage protocol reactor: the bench experiment of TMSPA in wet EC.

Stages run back to back; each is isothermal, its end state seeds the next, and rate constants are
recomputed at each stage temperature. TMSPA is injected as a discrete dilute-and-add event at the start
of the injection stage, which defines t = 0. Temperature changes are instantaneous.
"""

from typing import NamedTuple, Optional

import numpy as np
import pandas as pd

from kinetics.constants import ZERO_CELSIUS_K
from kinetics.microkinetics.models import get_model
from kinetics.data.network import NETWORK, NETWORK_SPECIES
from kinetics.reactor.engine import MassActionSystem, calculate_element_totals
from kinetics.data.snapshot import calculate_molar_mass

# Densities of the pure liquids at room temperature [g/mL] (TMSPA: CAS 1017-19-2)
DEFAULT_DENSITY_G_ML = {'EC': 1.321, 'H2O': 0.997, 'TMSPA': 0.940}


class Stage(NamedTuple):
    label: str
    T_C: float
    duration_h: float
    acquired_T_C: Optional[float] = None   # hold temperature reported by an NMR acquisition ending this stage


def calculate_recipe_molarities(h2o_vol_frac_stock: float = 0.02, tmspa_vol_frac: float = 0.05,
                                densities_g_mL: dict = None) -> dict:
    """Molarities of the two-step recipe: EC + H2O stock, then TMSPA added (which dilutes the stock).

    Returns 'stock' and 'after' ({species: M}), 'dilution_factor', 'h2o_tmspa_ratio' and 'summary' (DataFrame).
    """
    rho = densities_g_mL or DEFAULT_DENSITY_G_ML

    def molarity(vol_frac, sp):
        return vol_frac * 1000.0 * rho[sp] / calculate_molar_mass(sp)

    stock = {'H2O': molarity(h2o_vol_frac_stock, 'H2O'), 'EC': molarity(1.0 - h2o_vol_frac_stock, 'EC'), 'TMSPA': 0.0}
    dilution = 1.0 - tmspa_vol_frac
    after = {'H2O': stock['H2O'] * dilution, 'EC': stock['EC'] * dilution, 'TMSPA': molarity(tmspa_vol_frac, 'TMSPA')}
    summary = pd.DataFrame({
        'vol% (final)': {'EC': (1.0 - h2o_vol_frac_stock) * dilution * 100.0,
                         'H2O': h2o_vol_frac_stock * dilution * 100.0, 'TMSPA': tmspa_vol_frac * 100.0},
        'stock (M)': stock,
        'after injection (M)': after,
    })
    return {'stock': stock, 'after': after, 'dilution_factor': dilution,
            'h2o_tmspa_ratio': after['H2O'] / after['TMSPA'], 'summary': summary}


def build_protocol_schedule(rt_C: float = 20.0, stir_h: float = 24.0, equilibrate_h: float = 12.0,
                            hold_h: float = 8.0, acquire_h: float = 1.0,
                            hold_temps_C=(20, 30, 40, 50, 60, 70, 80)) -> tuple[list, pd.DataFrame]:
    """Stages of the bench protocol and their timetable (t = 0 at TMSPA addition, i.e. the first hold)."""
    stages = [Stage('EC + H2O, stirred', rt_C, stir_h), Stage('equilibrate', rt_C, equilibrate_h)]
    for T_C in hold_temps_C:
        stages.append(Stage(f'hold {T_C:g} °C', float(T_C), hold_h))
        stages.append(Stage(f'acquire ({T_C:g} °C)', rt_C, acquire_h, float(T_C)))

    rows, t_h = [], -(stir_h + equilibrate_h)
    for st in stages:
        rows.append({'stage': st.label, 'T (°C)': st.T_C, 'start (h)': round(t_h, 2),
                     'end (h)': round(t_h + st.duration_h, 2),
                     'acquisition': f'{st.acquired_T_C:.0f} °C' if st.acquired_T_C is not None else ''})
        t_h += st.duration_h
    return stages, pd.DataFrame(rows)


def simulate_protocol(stages: list = None, recipe: dict = None, *, model='peter_reference', network: dict = None,
                      species: list = None, species_db: dict = None, injection_stage_idx: int = 2,
                      buffered_species=('EC',), points_per_stage: int = 120, rtol: float = 1e-8,
                      atol: float = 1e-12, viscosity_Pa_s: float = None) -> dict:
    """Integrate the protocol stage by stage (Radau with analytic Jacobian).

    Returns time grids, concentrations, the temperature program, the state at each NMR acquisition
    ('acquisitions': [{'T_C', 't_h', 'idx', 'C_M'}]) and Si/P balances checked after injection.
    """
    spec = get_model(model)
    stages = stages if stages is not None else build_protocol_schedule()[0]
    recipe = recipe if recipe is not None else calculate_recipe_molarities()
    net = network or NETWORK
    species = list(species or NETWORK_SPECIES)
    buffered = {sp: recipe['after'][sp] for sp in buffered_species if sp in species}
    system = MassActionSystem(net, species, buffered)
    idx = system.idx

    c = np.zeros(len(species))
    for sp in ('H2O', 'EC'):
        c[idx[sp]] = recipe['stock'][sp]

    rate_cache = {}
    t_s = -sum(st.duration_h for st in stages[:injection_stage_idx]) * 3600.0
    t_parts, c_parts, T_parts, acquisitions = [], [], [], []
    for i, st in enumerate(stages):
        if i == injection_stage_idx:
            c = c * recipe['dilution_factor']
            c[idx['TMSPA']] = recipe['after']['TMSPA']
        T_K = st.T_C + ZERO_CELSIUS_K
        if T_K not in rate_cache:
            rates = spec.calculate_rates(T_K, network=net, species_db=species_db, viscosity_Pa_s=viscosity_Pa_s)
            rate_cache[T_K] = (rates['k_f'].values, rates['k_r'].values)
        duration_s = st.duration_h * 3600.0
        sol = system.integrate(c, *rate_cache[T_K], np.linspace(0.0, duration_s, points_per_stage),
                               rtol=rtol, atol=atol)
        t_parts.append(t_s + sol.t)
        c_parts.append(sol.y)
        T_parts.append(np.full(sol.t.shape, T_K))
        c = sol.y[:, -1]
        t_s += duration_s
        if st.acquired_T_C is not None:
            acquisitions.append({'T_C': st.acquired_T_C, 't_h': t_s / 3600.0,
                                 'idx': (i + 1) * points_per_stage - 1, 'C_M': c.copy()})

    t_all_s = np.concatenate(t_parts)
    C_M = np.concatenate(c_parts, axis=1)
    injection_idx = injection_stage_idx * points_per_stage
    totals = calculate_element_totals(C_M, species)
    return {
        'model': spec.name,
        't_s': t_all_s,
        't_h': t_all_s / 3600.0,
        'T_K': np.concatenate(T_parts),
        'C_M': C_M,
        'C_mM': C_M * 1000.0,
        'species': species,
        'idx': idx,
        'injection_idx': injection_idx,
        'acquisitions': acquisitions,
        'element_totals_M': totals,
        'elements_conserved': {el: bool(np.ptp(v[injection_idx:]) < 1e-5) for el, v in totals.items()},
        'recipe': recipe,
        'stages': stages,
    }
