"""Tests for the network helpers, the model registry, the mass-action engine and reactor observables.

Source: Y. Alcaraz Galván
"""

import os
import sys

import numpy as np

repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

from kinetics import (
    NETWORK, MODELS, NETWORK_SPECIES, calculate_cycle_residuals, calculate_network_thermo,
    describe_models, get_model, simulate_protocol, tabulate_acquisitions,
)
from kinetics.data.experimental import get_barrier_windows, get_measured_shifts
from kinetics.reactor.validation import evaluate_control_experiments
from kinetics.microkinetics.models import tabulate_family_barriers
from kinetics.data.network import find_reaction_cycles, list_species
from kinetics.reactor import (
    MassActionSystem, calculate_remaining_fraction, calculate_worst_case_pressure_bar, find_crossing_time,
)
from kinetics.data.snapshot import calculate_molar_mass


def test_network_cycles_and_wegscheider():
    cycles = find_reaction_cycles(NETWORK)
    assert cycles.shape == (3, 9), f"Expected 3 independent cycles, got {cycles.shape[0]}"
    expected = {(0, 3, 4): (1, 1, -1), (1, 3, 5): (1, 1, -1), (2, 3, 6): (1, 1, -1)}   # R_k + R4 - R_(k+4)
    for row in cycles:
        nz = tuple(np.flatnonzero(row))
        assert nz in expected and tuple(row[list(nz)]) == expected[nz], f"Unexpected cycle {row}"
    thermo = calculate_network_thermo()
    residuals = calculate_cycle_residuals(thermo['dG_rxn_kJ_mol'].to_dict())
    assert max(abs(r) for _, r in residuals) < 1e-6, residuals
    assert set(list_species(NETWORK)) == set(NETWORK_SPECIES)


def test_model_registry():
    table = describe_models()
    assert list(table.index) == list(MODELS)
    assert not table.loc['reference_bep', 'reversal invariant']
    assert table.loc['level1', 'reversal invariant']
    swept = get_model('reference_bep').with_barrier(1.30)
    assert swept.family_params['default']['g_eV'] == 1.30
    assert MODELS['reference_bep'].family_params['default']['g_eV'] == 1.15, "with_barrier must not mutate"
    try:
        get_model('nope')
    except KeyError:
        pass
    else:
        raise AssertionError("Unknown model name must raise")


def test_mass_action_jacobian():
    system = MassActionSystem(NETWORK, NETWORK_SPECIES, {'EC': 13.97})
    rates = get_model('family_marcus').calculate_rates(293.15)
    k_f, k_r = rates['k_f'].values, rates['k_r'].values
    C = np.random.default_rng(0).uniform(0.01, 0.2, len(NETWORK_SPECIES))
    J = system.jacobian(0.0, C, k_f, k_r)
    J_fd = np.zeros_like(J)
    for i in range(len(C)):
        h = 1e-7 * C[i]
        Cp, Cm = C.copy(), C.copy()
        Cp[i] += h
        Cm[i] -= h
        J_fd[:, i] = (system.rhs(0.0, Cp, k_f, k_r) - system.rhs(0.0, Cm, k_f, k_r)) / (2 * h)
    err = np.max(np.abs(J - J_fd)) / np.max(np.abs(J_fd))
    assert err < 1e-6, f"Analytic Jacobian disagrees with finite differences: {err:.1e}"


def test_observables():
    t = np.array([0.0, 1.0, 2.0, 3.0])
    assert find_crossing_time(t, np.array([1.0, 0.8, 0.4, 0.2]), 0.5) == 1.75
    assert np.isnan(find_crossing_time(t, np.ones(4), 0.5))
    sim = simulate_protocol(points_per_stage=30)
    table = tabulate_acquisitions(sim)
    assert list(table.index) == [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]
    frac = calculate_remaining_fraction(sim)
    assert 0.0 <= frac.min() and frac.max() <= 1.0 + 1e-9
    assert abs(calculate_molar_mass('TMSPA') - 314.54) < 0.01


def test_worst_case_pressure():
    # 1 mmol of gas at 300 K in 1 mL of headspace: p = nRT/V = 1e-3 · 8.314 · 300 / 1e-6 Pa ≈ 24.9 bar
    p_bar = calculate_worst_case_pressure_bar(1.0, 1.0, 1.0, 300.0)
    assert abs(p_bar - 1e-3 * 8.314462618 * 300.0 / 1e-6 / 1e5) < 1e-9
    assert calculate_worst_case_pressure_bar(1.0, 0.5, 1.0, 300.0) == 0.5 * p_bar


def test_experimental_reference_and_checks():
    shifts = get_measured_shifts()
    assert set(shifts['nucleus']) == {'P', 'Si', 'C', 'H'} and (shifts['low_ppm'] <= shifts['high_ppm']).all()
    assert shifts['peak'].notna().all() and shifts['source'].notna().all(), 'every shift names its peak and figure'
    constraints = get_barrier_windows()
    assert set(constraints.index) == {'R4', 'R5', 'R8'}
    checks = evaluate_control_experiments(['reference_bep', 'level1'])
    assert len(checks) == 8 and checks['predicted'].between(-1e-9, 1 + 1e-9).all()
    table = tabulate_family_barriers(['reference_bep', 'level1'])
    assert table.loc['solvent_attack', 'level1'] == 1.32 and table.loc['hydrolysis', 'reference_bep'] == 1.15


if __name__ == "__main__":
    test_network_cycles_and_wegscheider()
    test_model_registry()
    test_mass_action_jacobian()
    test_observables()
    test_worst_case_pressure()
    test_experimental_reference_and_checks()
    print("All model, network and engine tests passed.")
