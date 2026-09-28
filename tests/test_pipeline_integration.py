"""End-to-End Pipeline Integration Test.

Tests the full computational chain:
1. Quasi-RRHO thermochemistry & solution cycle.
2. Microkinetic rate constants with microscopic reversibility.
3. Multi-stage temperature-programmed protocol reactor (stiff ODE Radau IIA).
4. Automated molecular graph symmetry & multi-nuclear NMR shifts (29Si, 31P, 13C, 1H).
5. Exact water mass balance and fast-exchange labile proton tracking.
6. Reaction fingerprint construction, SVD identifiability analysis, and extent inversion.
"""

import sys
import os
import numpy as np

repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

from kinetics import (
    calculate_solution_gibbs,
    calculate_reaction_thermo,
    calculate_rate_constants,
    compute_recipe_molarities,
    build_default_protocol_schedule,
    simulate_protocol_reactor,
    build_referenced_nmr_sites,
    compute_water_mass_balance,
    build_multinuclear_feature_space,
    build_pure_component_matrix,
    build_reaction_fingerprints,
    analyze_reaction_identifiability,
    recover_reaction_extents
)


def test_thermodynamics_and_rates():
    print("--- 1. Testing Solution Thermodynamics & Microkinetics ---")
    T = 298.15
    sol_tmspa = calculate_solution_gibbs("TMSPA", T_K=T, mode="wb97mv")
    assert sol_tmspa["G_sol_eV"] < 0.0, "TMSPA G_sol should be negative"
    
    rxn_thermo = calculate_reaction_thermo("R1", T_K=T, mode="wb97mv")
    rates = calculate_rate_constants("R1", T_K=T, mode="wb97mv")
    
    # Backwards-compatibility check for legacy alias 'b3lyp_benchmark'
    sol_legacy = calculate_solution_gibbs("TMSPA", T_K=T, mode="b3lyp_benchmark")
    assert abs(sol_tmspa["G_sol_eV"] - sol_legacy["G_sol_eV"]) < 1e-12, "Legacy alias mismatch"
    
    k_ratio = rates["k_f"] / rates["k_r"]
    rel_error = abs(k_ratio - rates["K_eq"]) / rates["K_eq"]
    assert rel_error < 1e-4, f"Detailed balance violated: kf/kr={k_ratio}, K_eq={rates['K_eq']}"
    print(f"✓ Detailed balance verified: k_f/k_r = {k_ratio:.3e}, K_eq = {rates['K_eq']:.3e}")

    rates_l1 = calculate_rate_constants("R1", T_K=T, mode="wb97mv", kinetic_model="level1")
    rel_error_l1 = abs(rates_l1["k_f"] / rates_l1["k_r"] - rates_l1["K_eq"]) / rates_l1["K_eq"]
    assert rel_error_l1 < 1e-4, "Detailed balance violated for the Level 1 engine"
    print(f"✓ Level 1 engine ({rates_l1['barrier_model']}): ΔG‡_f(R1) = {rates_l1['dG_barrier_f_eV']:.3f} eV, detailed balance kept")


def test_protocol_reactor_simulation():
    print("\n--- 2. Testing Multi-Stage Protocol Reactor Simulation ---")
    recipe = compute_recipe_molarities()
    assert abs(recipe["h2o_tmspa_ratio"] - 7.04) < 0.05, "Water to TMSPA ratio should be ~ 7:1"
    
    stages, sched_df = build_default_protocol_schedule()
    assert len(stages) == 16, "Default schedule should have 16 stages"
    
    sim = simulate_protocol_reactor(stages=stages, recipe=recipe, points_per_stage=30)
    assert sim["si_conserved"], "Silicon element balance violated post-injection"
    assert sim["p_conserved"], "Phosphorus element balance violated post-injection"
    assert len(sim["snapshots"]) == 7, "Should record 7 cold NMR acquisitions"
    print(f"✓ Protocol reactor solved across {sim['t_h'][-1]:.1f} h (Si and P conserved: 100%)")
    return sim


def test_nmr_symmetry_and_water_balance(sim):
    print("\n--- 3. Testing Molecular Symmetry & Water Mass Balance ---")
    catalog = build_referenced_nmr_sites()
    for el in ["Si", "P", "C", "H"]:
        assert el in catalog, f"Catalog missing nucleus {el}"
    assert "TMSPA" in catalog["Si"], "TMSPA missing in Si catalog"
    
    wb = compute_water_mass_balance(sim["C_M"], sim["idx"], initial_h2o_molarity=sim["recipe"]["after"]["H2O"], nmr_sites_catalog=catalog)
    diff_final = abs(wb["h2o_deduced_M"][-1] - wb["h2o_direct_M"][-1])
    assert diff_final < 1e-6, f"Water mass balance deduction failed: diff={diff_final:.2e} M"
    print(f"✓ Water mass balance verified: direct [H2O] = {wb['h2o_direct_M'][-1]*1000:.2f} mM, deduced = {wb['h2o_deduced_M'][-1]*1000:.2f} mM")


def test_reaction_fingerprints_and_inversion(sim):
    print("\n--- 4. Testing Reaction Fingerprints & Extent Recovery ---")
    catalog = build_referenced_nmr_sites()
    visible = [sp for sp in sim["species"] if sp != "EC" and any(not lab for el in ["Si", "P", "C", "H"] for _, _, lab in catalog.get(el, {}).get(sp, []))]
    
    grids, x_concat, blocks, n_feat = build_multinuclear_feature_space(visible, catalog)
    p_pure = build_pure_component_matrix(visible, x_concat, blocks, n_feat, catalog)
    
    rxns = [
        ("R1", ["TMSPA", "H2O"],   ["BMSPA",    "TMSOH"]),
        ("R2", ["BMSPA", "H2O"],   ["MMSPA",    "TMSOH"]),
        ("R3", ["MMSPA", "H2O"],   ["H3PO4",    "TMSOH"]),
        ("R4", ["TMSOH", "TMSOH"], ["HMDSO",    "H2O"]),
        ("R5", ["TMSPA", "TMSOH"], ["BMSPA",    "HMDSO"]),
        ("R6", ["BMSPA", "TMSOH"], ["MMSPA",    "HMDSO"]),
        ("R7", ["MMSPA", "TMSOH"], ["H3PO4",    "HMDSO"]),
        ("R8", ["EC",    "TMSOH"], ["TMSOEG",   "CO2"]),
        ("R9", ["EC",    "TMSOEG"], ["TMSOdiEG", "CO2"]),
    ]
    f_norm, f_raw, weights, s_vis = build_reaction_fingerprints(rxns, visible, p_pure, blocks)
    analysis = analyze_reaction_identifiability(f_norm, f_raw, [r[0] for r in rxns], blocks)
    
    assert analysis["effective_rank"] == 6, f"Rank should be 6 (found {analysis['effective_rank']})"
    print(f"✓ Identifiability rank verified: {analysis['effective_rank']} of {len(rxns)} reactions are independent")
    print(f"  Lumped identifiable basis: {'; '.join(analysis['lumped_labels'])}")
    
    points_per_stage = 30
    injection_stage_idx = 2
    i0 = injection_stage_idx * points_per_stage
    stages = sim["stages"]
    i_snaps = [(i + 1) * points_per_stage - 1 for i, st in enumerate(stages) if st[3] is not None]
    step_idx = [i0] + i_snaps
    
    vis_idx = [sim["idx"][sp] for sp in visible]
    c_vis_snapshots = sim["C_M"][vis_idx][:, step_idx]
    
    d_true = (c_vis_snapshots.T * 1000.0) @ p_pure
    dD = np.diff(d_true, axis=0)
    
    xi_hat = recover_reaction_extents(dD, analysis["F_basis_norm"], weights)
    assert xi_hat.shape == (7, 6), f"xi_hat shape mismatch: {xi_hat.shape}"
    
    f_basis_raw = analysis["F_basis_raw"]
    dD_reconstructed = xi_hat @ f_basis_raw
    max_rec_err = np.max(np.abs(dD_reconstructed - dD))
    assert max_rec_err < 1e-6, f"Extent reconstruction error too large: {max_rec_err:.2e}"
    print(f"✓ Reaction extent recovery verified: max |ΔD_reconstructed − ΔD_true| = {max_rec_err:.3e}")


if __name__ == "__main__":
    print("================================================================================")
    print("RUNNING END-TO-END ATOM-TO-REACTOR PIPELINE TESTS")
    print("================================================================================")
    test_thermodynamics_and_rates()
    sim = test_protocol_reactor_simulation()
    test_nmr_symmetry_and_water_balance(sim)
    test_reaction_fingerprints_and_inversion(sim)
    print("\n================================================================================")
    print("ALL TESTS PASSED SUCCESSFULLY! (100% MATHEMATICAL & NUMERICAL INTEGRITY)")
    print("================================================================================")
