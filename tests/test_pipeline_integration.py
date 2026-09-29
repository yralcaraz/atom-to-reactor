"""End-to-end pipeline integration test.

1. Solution thermochemistry and rate constants with detailed balance.
2. Multi-stage protocol reactor (Radau IIA) with Si/P conservation.
3. NMR catalog from DFT shieldings and the OH mass balance for water.
4. Reaction fingerprints, identifiability and extent recovery.
"""

import os
import sys

import numpy as np

repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

from kinetics import (
    NETWORK, build_nmr_catalog, build_protocol_schedule, calculate_rate_constants,
    calculate_recipe_molarities, calculate_solution_gibbs, calculate_water_mass_balance,
    recover_reaction_extents, run_fingerprint_analysis, simulate_acquisition_spectra, simulate_protocol,
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


def test_nmr_catalog_and_water_balance(sim):
    print("\n--- 3. NMR catalog and water mass balance ---")
    catalog = build_nmr_catalog()
    assert set(catalog) == {"Si", "P", "C", "H"}, "Catalog must cover the four nuclei"
    assert "TMSPA" in catalog["Si"], "TMSPA missing in the Si catalog"
    wb = calculate_water_mass_balance(sim["C_M"], sim["idx"], sim["recipe"]["after"]["H2O"], catalog)
    post = slice(sim["injection_idx"], None)
    diff = np.max(np.abs(wb["h2o_mass_balance_M"][post] - wb["h2o_direct_M"][post]))
    assert diff < 1e-6, f"Water mass balance failed: max diff = {diff:.2e} M"
    print(f"✓ Water from the OH mass balance matches the ODE state (max diff {diff:.1e} M)")


def test_fingerprints_and_extent_recovery(sim):
    print("\n--- 4. Reaction fingerprints and extent recovery ---")
    fp = run_fingerprint_analysis(NETWORK)
    assert fp["effective_rank"] == 6, f"Rank should be 6 (found {fp['effective_rank']})"
    print(f"✓ Rank {fp['effective_rank']} of {len(NETWORK)}: {'; '.join(fp['lumped_labels'])}")

    dD = np.diff(simulate_acquisition_spectra(sim, fp), axis=0)
    xi = recover_reaction_extents(dD, fp["F_basis_norm"], fp["weights"])
    assert xi.shape == (7, 6), f"xi shape mismatch: {xi.shape}"
    err = np.max(np.abs(xi @ fp["F_basis_raw"] - dD))
    assert err < 1e-6, f"Extent reconstruction error too large: {err:.2e}"
    print(f"✓ Noise-free spectra reconstructed from recovered extents (max error {err:.1e})")


if __name__ == "__main__":
    test_thermodynamics_and_rates()
    sim = test_protocol_reactor()
    test_nmr_catalog_and_water_balance(sim)
    test_fingerprints_and_extent_recovery(sim)
    print("\nALL PIPELINE TESTS PASSED")
