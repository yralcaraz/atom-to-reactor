"""End-to-end pipeline integration test.

1. Solution thermochemistry and rate constants with detailed balance.
2. Multi-stage protocol reactor (Radau IIA) with Si/P conservation.

Source: Y. Alcaraz Galván
"""

import os
import sys

repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

from kinetics import (
    build_protocol_schedule, calculate_rate_constants, calculate_recipe_molarities, calculate_solution_gibbs,
    simulate_protocol,
)


def test_thermodynamics_and_rates():
    print("--- 1. Solution thermodynamics and rate constants ---")
    T = 298.15
    assert calculate_solution_gibbs("TMSPA", T)["G_sol_eV"] < 0.0, "TMSPA G_sol should be negative"
    for kinetic_model in ('bep_eyring', 'level1'):
        rates = calculate_rate_constants("R1", T, kinetic_model=kinetic_model)
        rel_error = abs(rates["k_f"] / rates["k_r"] - rates["K_eq"]) / rates["K_eq"]
        assert rel_error < 1e-4, f"Detailed balance violated for {kinetic_model}"
        print(f"✓ {kinetic_model}: ΔG‡_f(R1) = {rates['dG_barrier_f_eV']:.3f} eV, k_f/k_r = K_eq")


def test_protocol_reactor():
    print("\n--- 2. Multi-stage protocol reactor ---")
    recipe = calculate_recipe_molarities()
    assert abs(recipe["h2o_tmspa_ratio"] - 7.04) < 0.05, "Water to TMSPA ratio should be ~ 7:1"
    stages, _ = build_protocol_schedule()
    assert len(stages) == 16, "Default schedule should have 16 stages"
    sim = simulate_protocol(stages, recipe, points_per_stage=30)
    assert all(sim["elements_conserved"].values()), f"Element balance violated: {sim['elements_conserved']}"
    assert len(sim["acquisitions"]) == 7, "Should record 7 NMR acquisitions"
    print(f"✓ Protocol solved over {sim['t_h'][-1]:.1f} h; Si and P conserved")
    return sim


if __name__ == "__main__":
    test_thermodynamics_and_rates()
    test_protocol_reactor()
    print("\nALL PIPELINE TESTS PASSED")
