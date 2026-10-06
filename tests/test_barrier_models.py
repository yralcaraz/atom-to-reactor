"""Block 5 Level 1 barrier-engine tests.

Covers:
1. Axioms of every smooth barrier model (reversal invariance, bounds, anchoring, Leffler bounds).
2. Analytic limits (Blowers-Masel → Marcus, two-parabola(½) ≡ Marcus, Agmon-Levine asymptotes).
3. Bit-identity of the legacy BEP / Marcus paths against the pre-refactor baseline.
4. Network reversal: identical dynamics when a reaction is written backwards (legacy BEP fails).
5. Temperature-dependent intrinsic barrier, Level 1 calibration and the diffusion ceiling.

Source: Y. Alcaraz Galván; expected barriers computed from the Tank dataset snapshot (P. Broqvist)
"""

import json
import os
import sys
import warnings
from copy import deepcopy

import numpy as np

repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

from kinetics import NETWORK, ModelSpec, calculate_rate_constants, load_species_database, simulate_batch_reactor
from kinetics.constants import EV_TO_KJ_MOL
from kinetics.paths import data_path
from kinetics.microkinetics import (
    FAMILY_BEP_PARAMETERS, PETER_REFERENCE_PARAMETERS, calculate_barrier, invert_marcus_barrier, reversed_params,
)

X_GRID = np.linspace(-3.0, 3.0, 121)
SMOOTH_CASES = [
    ('marcus', {}),
    ('agmon_levine', {}),
    ('blowers_masel', {'w': 5.0}),
    ('two_parabola', {'alpha0': 0.3}),
    ('two_parabola', {'alpha0': 0.7}),
]
SPECIES_DB = load_species_database()


def _params_with_work(extra):
    return {'wR_eV': 0.04, 'wP_eV': -0.02, **extra}


def test_axioms_smooth_models():
    print("--- 1. Axioms A2-A6 for smooth barrier models ---")
    g = 0.8
    for shape, extra in SMOOTH_CASES:
        p = _params_with_work(extra)
        p_rev = reversed_params(p)
        for x in X_GRID:
            f, alpha = calculate_barrier(x, shape, g, **p)
            f_rev, _ = calculate_barrier(-x, shape, g, **p_rev)
            assert abs(f_rev - (f - x)) < 1e-12, f"A2 reversal violated for {shape} {extra} at x={x}"
            assert f >= max(0.0, x) - 1e-12, f"A3 bounds violated for {shape} at x={x}"
            if shape != 'marcus' or abs(x) <= 4 * g:
                assert -1e-9 <= alpha <= 1 + 1e-9, f"A6 Leffler bounds violated for {shape} at x={x}"
        f0, _ = calculate_barrier(0.0, shape, g, **extra)
        assert abs(f0 - g) < 1e-12, f"A4 anchoring F(0)=g violated for {shape}"
    print("✓ Marcus, Agmon-Levine, Blowers-Masel and two-parabola satisfy A2-A6")


def test_legacy_cap_is_direction_dependent():
    print("\n--- 2. Legacy BEP cap violates reversal invariance (documented defect) ---")
    x, g = 0.102, 0.8
    f, _ = calculate_barrier(x, 'bep_cap', g, alpha=0.5)
    f_rev, _ = calculate_barrier(-x, 'bep_cap', g, alpha=0.5)
    assert abs(f_rev - (f - x)) > 0.05, "bep_cap unexpectedly reversal invariant"
    print(f"✓ bep_cap: reversed-direction barrier differs by {abs(f_rev - (f - x)):.4f} eV")


def test_limits():
    print("\n--- 3. Analytic limits ---")
    g = 0.8
    for x in np.linspace(-3.0, 3.0, 61):
        bm, _ = calculate_barrier(x, 'blowers_masel', g, w=1000.0)
        m, _ = calculate_barrier(x, 'marcus', g)
        assert abs(bm - m) < 2e-4, f"BM(w→∞) does not approach Marcus at x={x}"
        tp, _ = calculate_barrier(x, 'two_parabola', g, alpha0=0.5)
        assert tp == m, "two_parabola(α0=½) must equal Marcus"
    assert abs(calculate_barrier(-40.0, 'agmon_levine', g)[0]) < 1e-12, "AL must tend to 0 for x → -∞"
    assert abs(calculate_barrier(40.0, 'agmon_levine', g)[0] - 40.0) < 1e-12, "AL must tend to x for x → +∞"
    for x in [-0.474, -0.047, 0.102]:
        for dg in [0.9, 1.3]:
            g_fit = invert_marcus_barrier(dg, x)
            assert abs(calculate_barrier(x, 'marcus', g_fit)[0] - dg) < 1e-12, "invert_marcus round-trip failed"
    print("✓ BM(w=1000) ≈ Marcus, two-parabola(½) ≡ Marcus, AL asymptotes, Marcus inversion")


def test_legacy_bit_identity():
    print("\n--- 4. Legacy BEP/Marcus reproduce the pre-refactor baseline ---")
    path = data_path('tests', 'block5_legacy_baseline.json')
    if not os.path.exists(path):
        print(f"  - skipped: {path} not found (set ATOM_DATA_DIR, see .env.example)")
        return
    records = json.load(open(path))['records']
    params = {'family_default': FAMILY_BEP_PARAMETERS, 'peter_reference': PETER_REFERENCE_PARAMETERS}
    worst = 0.0
    for r in records:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')   # qRRHO records fall back to the stored G (no frequencies)
            res = calculate_rate_constants(
                r['rxn'], r['T_K'], kinetic_model=r['model'], family_params=params[r['params']],
                network=NETWORK, species_db=SPECIES_DB, thermo_mode=r['mode']
            )
        for key in ['k_f', 'k_r', 'dG_barrier_f_eV']:
            worst = max(worst, abs(res[key] - r[key]) / max(abs(r[key]), 1e-300))
    assert worst < 1e-12, f"Legacy regression: worst relative deviation {worst:.2e}"
    print(f"✓ {len(records)} legacy records reproduced (worst relative deviation {worst:.1e})")


def _tank(net, kinetic_model, family_params):
    species = ['TMSPA', 'H2O', 'BMSPA', 'MMSPA', 'H3PO4', 'TMSOH', 'HMDSO', 'EC', 'TMSOEG', 'TMSOdiEG', 'CO2']
    c0_M = {s: 0.0 for s in species}
    c0_M.update({'TMSPA': 0.05, 'H2O': 0.02, 'EC': 4.5})
    model = ModelSpec('test', 'test', kinetic_model, family_params)
    return simulate_batch_reactor(
        c0_M, T_K=298.15, t_end_s=1e5, model=model, network=net, species_db=SPECIES_DB,
        rtol=1e-10, atol=1e-14, n_points=200
    )


def test_network_reversal_invariance():
    print("\n--- 5. Dynamics independent of the direction R4 is written ---")
    reversed_net = deepcopy(NETWORK)
    reversed_net['R4'] = {'reactants': {'HMDSO': 1, 'H2O': 1}, 'products': {'TMSOH': 2},
                          'class': 'condensation', 'canonical': False}
    level1_params = {
        'hydrolysis':     {'shape': 'marcus', 'g_eV': 0.80},
        'transfer':       {'shape': 'marcus', 'g_eV': 0.80},
        'condensation':   {'shape': 'two_parabola', 'g_eV': 0.80, 'alpha0': 0.3, 'wR_eV': 0.03, 'wP_eV': 0.0},
        'solvent_attack': {'shape': 'marcus', 'g_eV': 1.32},
        'default':        {'shape': 'marcus', 'g_eV': 0.80},
    }
    ref = _tank(NETWORK, 'level1', level1_params)
    rev = _tank(reversed_net, 'level1', level1_params)
    dev = np.max(np.abs(ref['C_M'] - rev['C_M']))
    assert dev < 1e-8, f"Level 1 dynamics depend on reaction direction: max |ΔC| = {dev:.2e} M"
    print(f"✓ Level 1 (asymmetric α0 + work terms): max |ΔC| = {dev:.1e} M")

    ref_bep = _tank(NETWORK, 'bep_eyring', FAMILY_BEP_PARAMETERS)
    rev_bep = _tank(reversed_net, 'bep_eyring', FAMILY_BEP_PARAMETERS)
    dev_bep = np.max(np.abs(ref_bep['C_M'] - rev_bep['C_M']))
    assert dev_bep > 1e-4, "Legacy BEP unexpectedly direction independent"
    print(f"✓ Legacy BEP (documented defect): max |ΔC| = {dev_bep:.2e} M")


def test_temperature_and_level1_calibration():
    print("\n--- 6. g(T) with activation entropy and Level 1 calibrated barriers ---")
    params = {'default': {'shape': 'marcus', 'g_eV': 1.0, 'T_ref_K': 298.15, 'dS_act_J_mol_K': -100.0}}
    r298 = calculate_rate_constants('R8', 298.15, kinetic_model='level1', family_params=params,
                                    network=NETWORK, species_db=SPECIES_DB, thermo_mode='wb97mv')
    r353 = calculate_rate_constants('R8', 353.15, kinetic_model='level1', family_params=params,
                                    network=NETWORK, species_db=SPECIES_DB, thermo_mode='wb97mv')
    expected = 55.0 * 100.0 / 1000.0 / EV_TO_KJ_MOL
    assert abs((r353['g_eV'] - r298['g_eV']) - expected) < 1e-12, "g(T) slope must equal -ΔS‡"
    assert abs(expected - 0.057) < 1e-3

    targets = {'R1': 0.580, 'R4': 1.351, 'R8': 1.297, 'R9': 1.097}
    for rxn, target in targets.items():
        res = calculate_rate_constants(rxn, 298.15, kinetic_model='level1',
                                       network=NETWORK, species_db=SPECIES_DB, thermo_mode='wb97mv')
        assert abs(res['dG_barrier_f_eV'] - target) < 2e-3, f"{rxn}: {res['dG_barrier_f_eV']:.4f} vs {target}"
        rel = abs(res['k_f'] / res['k_r'] - res['K_eq']) / res['K_eq']
        assert rel < 1e-10, f"Detailed balance violated for {rxn}"
    print(f"✓ g(353 K) - g(298 K) = {expected:.4f} eV for ΔS‡ = -100 J/mol/K; Level 1 barriers match formulation")


def test_diffusion_ceiling():
    print("\n--- 7. Collins-Kimball diffusion ceiling ---")
    params = {'default': {'shape': 'marcus', 'g_eV': 0.05}}
    res = calculate_rate_constants('R1', 298.15, kinetic_model='level1', family_params=params,
                                   network=NETWORK, species_db=SPECIES_DB,
                                   thermo_mode='wb97mv', viscosity_Pa_s=1.9e-3)
    assert 3.0e9 < res['k_D'] < 4.0e9, f"k_D out of expected range: {res['k_D']:.2e}"
    assert res['k_f'] <= res['k_D'], "k_f exceeds the diffusion limit"
    rel = abs(res['k_f'] / res['k_r'] - res['K_eq']) / res['K_eq']
    assert rel < 1e-10, "Detailed balance violated after the diffusion ceiling"
    print(f"✓ k_f = {res['k_f']:.2e} ≤ k_D = {res['k_D']:.2e} M⁻¹s⁻¹, detailed balance kept")


if __name__ == "__main__":
    print("================================================================================")
    print("RUNNING BLOCK 5 LEVEL 1 BARRIER-ENGINE TESTS")
    print("================================================================================")
    test_axioms_smooth_models()
    test_legacy_cap_is_direction_dependent()
    test_limits()
    test_legacy_bit_identity()
    test_network_reversal_invariance()
    test_temperature_and_level1_calibration()
    test_diffusion_ceiling()
    print("\n================================================================================")
    print("ALL BARRIER-ENGINE TESTS PASSED")
    print("================================================================================")
