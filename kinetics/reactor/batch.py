"""Block 7 — Isothermal batch reactor ('tank'): stiff mass-action ODEs at one temperature.

Classification: PUBLIC (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván
"""

import numpy as np

from kinetics.microkinetics.models import get_model
from kinetics.data.network import NETWORK
from kinetics.reactor.engine import MassActionSystem, calculate_element_totals


def simulate_batch_reactor(c0_M: dict, *, T_K: float = 298.15, t_end_s: float = 1e6, model='family_bep',
                           network: dict = None, species_db: dict = None, buffered_species=('EC',),
                           n_points: int = 500, t_start_s: float = 1e-3, method: str = 'Radau',
                           rtol: float = 1e-8, atol: float = 1e-12, viscosity_Pa_s: float = None) -> dict:
    """Integrate dC/dt = S·r(C, T) from t = 0 to t_end_s at constant T_K.

    c0_M: initial concentrations [M] of every tracked species (the network species at least).
    model: registry name or ModelSpec (kinetics.microkinetics.models). Buffered species stay at their initial value.
    The output grid is log-spaced from t_start_s, so fast transients after t = 0 are resolved.
    """
    spec = get_model(model)
    net = network or spec.network or NETWORK
    species = list(c0_M)
    system = MassActionSystem(net, species, {sp: c0_M[sp] for sp in buffered_species if sp in c0_M})

    rates = spec.calculate_rates(T_K, network=net, species_db=species_db, viscosity_Pa_s=viscosity_Pa_s)
    t_eval = np.logspace(np.log10(t_start_s), np.log10(t_end_s), n_points)
    sol = system.integrate([c0_M[sp] for sp in species], rates['k_f'].values, rates['k_r'].values, t_eval,
                           method=method, rtol=rtol, atol=atol)

    totals = calculate_element_totals(sol.y, species)
    return {
        'model': spec.name,
        't_s': sol.t,
        't_h': sol.t / 3600.0,
        'T_K': np.full(sol.t.shape, T_K),
        'C_M': sol.y,
        'C_mM': sol.y * 1000.0,
        'species': species,
        'idx': system.idx,
        'rates': rates,
        'element_totals_M': totals,
        'elements_conserved': {el: bool(np.max(np.abs(v - v[0])) < 1e-6) for el, v in totals.items()},
        'success': sol.success,
        'message': sol.message,
    }
