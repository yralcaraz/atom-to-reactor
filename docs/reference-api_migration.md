# API migration: harmonized `kinetics` package

| | |
|---|---|
| **Status** | Stable |
| **Source** | Y. Alcaraz Galván |
| **Scope** | `MODULES.md`, `README.md` and `docs/theory-multiscale_microkinetics.md` still use the names on the left. This document maps them to the current code. The numerical results did not change (see *Verification*). |
| **Applies to** | Package layout of 2026-09-29 (commit `36faa93`) |

## Conventions

- **Verbs.**
  - `load_` reads from disk.
  - `build_` assembles a structure (network, schedule, catalog, matrices).
  - `calculate_` evaluates a physical quantity.
  - `simulate_` integrates in time or synthesizes signals.
  - `fit_`, `analyze_`, `recover_` do regression, diagnostics and inversion.
  - `find_` searches (cycles, windows, crossing times).
  - `format_` / `plot_` exist only in `demo/`.
- **Arguments.**
  - The same argument name is used everywhere: `thermo_mode`, `kinetic_model`, `family_params`, `network`, `species_db`, `nmr_catalog`.
  - Options come after `*`, so a misspelled keyword raises an error instead of being silently swallowed by `**kwargs`.
- **Units in names.**
  - Suffixes: `_eV`, `_Eh`, `_kJ_mol`, `_J_mol_K`, `_K`, `_M`, `_mM`, `_ppm`, `_s`, `_h`.
  - Conversions read `<FROM>_TO_<TO>`.
  - All constants live in `kinetics/constants.py`.
- **No hidden state.** Functions no longer read notebook globals through stack-frame inspection; every input is an argument with an explicit default.

## Package layout: one folder per pipeline stage

| Folder | Blocks | Files |
|---|---|---|
| `kinetics/constants.py` | all | physical constants and unit conversions |
| `kinetics/data/` | 1 | `snapshot.py` (computed-data snapshot), `species.py` (species database), `network.py` (reactions, stoichiometry, Wegscheider cycles), `experimental.py` (Gogoi 2024 reference data) |
| `kinetics/thermo/` | 2–4 | `gas.py`, `standard_state.py`, `solution.py`, `uncertainty.py`, `reaction.py` |
| `kinetics/microkinetics/` | 5–6 | `barriers.py`, `parameters.py`, `rates.py`, `models.py` (`ModelSpec`, `MODELS`), `arrhenius.py` |
| `kinetics/reactor/` | 7, 10 | `engine.py` (mass action + Jacobian), `batch.py`, `protocol.py`, `observables.py`, `validation.py` (control experiments) |
| `kinetics/spectroscopy/` | 8, 11–13 | `symmetry.py`, `spectra.py`, `fingerprints.py` |
| `demo/` | notebooks | `style.py`, `tables.py`, `thermo_plots.py`, `rate_plots.py`, `reactor_plots.py`, `nmr_plots.py`, `fingerprint_plots.py` (no science) |

### Moved files

| Old path | New path |
|---|---|
| `kinetics/thermo/species_data.py` | `kinetics/data/species.py` |
| `kinetics/thermo/gas_thermo.py` | `kinetics/thermo/gas.py` |
| `kinetics/thermo/solution_gibbs.py` | `kinetics/thermo/solution.py` |
| `kinetics/thermo/reaction_thermo.py` | `kinetics/thermo/reaction.py` |
| `kinetics/thermo/solvation_uncertainty.py` | `kinetics/thermo/uncertainty.py` |
| `kinetics/microkinetics/barrier_models.py` | `kinetics/microkinetics/barriers.py` |
| `kinetics/microkinetics/kinetic_parameters.py` | `kinetics/microkinetics/parameters.py` |
| `kinetics/microkinetics/rate_constants.py` | `kinetics/microkinetics/rates.py` |
| `kinetics/reactor/batch_reactor.py` | `kinetics/reactor/batch.py` |
| `kinetics/reactor/protocol_reactor.py` | `kinetics/reactor/protocol.py` |
| `kinetics/spectroscopy/molecular_symmetry.py` | `kinetics/spectroscopy/symmetry.py` |
| `kinetics/spectroscopy/multinuclear_nmr.py` | `kinetics/spectroscopy/spectra.py` |
| `kinetics/spectroscopy/reaction_fingerprints.py` | `kinetics/spectroscopy/fingerprints.py` |

## Renamed functions

| Old | New |
|---|---|
| `load_default_species_database` | `load_species_database` (returns a copy; the default database is cached) |
| `reaction_solvation_sigma`, `reaction_solvation_sigma_naive` | `calculate_solvation_sigma`, `calculate_solvation_sigma_naive` |
| `solvation_sigma_components` | `calculate_solvation_sigma_components` |
| `barrier`, `invert_marcus` | `calculate_barrier`, `invert_marcus_barrier` |
| `generate_arrhenius_summary` (list of strings) | `build_arrhenius_table` (numeric DataFrame) |
| `simulate_tank_reactor(C0_dict, ...)` | `simulate_batch_reactor(c0_M, *, T_K, t_end_s, model, ...)` |
| `simulate_protocol_reactor` | `simulate_protocol(stages, recipe, *, model, ...)` |
| `compute_recipe_molarities` | `calculate_recipe_molarities` (`summary_df` → `summary`) |
| `build_default_protocol_schedule(step_temps_c=...)` | `build_protocol_schedule(hold_temps_C=...)`; stages are `Stage` named tuples |
| `build_stoichiometric_matrix` (reactor) | `kinetics.data.build_stoichiometric_matrix` |
| `load_dataset_snapshot` | `kinetics.data.load_snapshot` |
| `extract_nmr_sites` | `extract_shielding_sites` (returns `[]`, never `None`) |
| `build_referenced_nmr_sites` | `build_nmr_catalog` (cached) |
| `site_peaks`, `auto_regions` | `calculate_nmr_peaks`, `find_shift_windows` |
| `simulate_multinuclear_spectra(..., lw=)` | `simulate_nmr_spectra(..., fwhm_ppm=)` |
| `compute_water_mass_balance` | `calculate_water_mass_balance`: keys `h2o_mass_balance_M`, `h2o_direct_M`, `oh_total_M`, `oh_carriers_mM` |
| `build_multinuclear_feature_space` (4-tuple) | `build_feature_space` → `FeatureSpace` dataclass |
| `build_pure_component_matrix` | `build_pure_spectra` |
| `build_reaction_fingerprints(reactions_list, ...)` | `build_reaction_fingerprints(network, ...)` → dict |
| — | `run_fingerprint_analysis`, `simulate_acquisition_spectra`, `find_visible_species` |
| — | `tabulate_family_barriers` (models.py), `calculate_eyring_rate`, `tabulate_trajectory` |
| — | `load_experimental_data`, `get_measured_shifts`, `get_barrier_windows` (`kinetics.data`), `evaluate_control_experiments` (`kinetics.reactor`) |

## Renamed arguments, keys and constants

| Old | New |
|---|---|
| `mode=`, `thermo_mode=` (both accepted) | `thermo_mode=` (default `'wb97mv'`; `'b3lyp_benchmark'` alias removed) |
| `model=`, `kinetic_model=` (both accepted) | `kinetic_model=` for rate constants; `model=` (name or `ModelSpec`) for reactors |
| `bep_params=` | `family_params=` (also used by Level 1) |
| `reactions_net=` | `network=` |
| `nu_0=` | `nu0_cm1=` |
| family key `E0_eV`, `dS_act_J_molK` | `g_eV` (`E0_eV` still accepted), `dS_act_J_mol_K`; unknown keys now raise |
| `DEFAULT_FAMILY_BEP_PARAMETERS`, `LEVEL1_FAMILY_PARAMETERS` | `FAMILY_BEP_PARAMETERS`, `LEVEL1_PARAMETERS` |
| `DEFAULT_REACTIONS_NETWORK`, `DEFAULT_SPECIES`, `DEFAULT_TRACKED_SPECIES` | `NETWORK`, `NETWORK_SPECIES` (`kinetics.data`) |
| species keys `G_wb97mv_eV`, `G_B3_eV`, `E_0K_eV` | `G_gas_eV`, `G_gas_Eh`, `H_gas_eV`, `E_scf_eV` (qRRHO now starts from the SCF energy) |
| result keys `si_total_M`, `si_conserved`, `snapshots` | `element_totals_M`, `elements_conserved`, `acquisitions` (with `idx` into the time grid) |
| `HAR2EV`, `KJMOL2EV`, `EV2KCAL`, `C_CMS`, `P_REF` | `HARTREE_TO_EV`, `KJ_MOL_TO_EV`, `EV_TO_KCAL_MOL`, `C_CM_S`, `P_STD_PA` |

## Removed

- **`simulate_virtual_nmr` and `DEFAULT_NMR_29SI`.** These held B3LYP-era ²⁹Si shifts that contradict the snapshot (Finding 13). Use `simulate_nmr_spectra` with `build_nmr_catalog()`.
- **Per-model wrappers `bep_eyring`, `marcus_eyring`, `level1_eyring`, the `'bep'`/`'marcus'` aliases and `lambda_eV`.** The dispatcher reads the declarative `KINETIC_MODELS` table instead.
- **The `siloxyl` alias for HMDSO.** The pipeline uses snapshot names only.
- **Silent fallbacks, which now raise instead:**
  - species absent from the reactor were treated as 1 M;
  - missing solvation energies counted as 0 eV;
  - molar masses defaulted to 100 g/mol;
  - an unknown `mode` fell through to qRRHO.

## Behaviour changes worth knowing

- **Default `thermo_mode` is `'wb97mv'` everywhere.** qRRHO needs frequencies the snapshot does not have. For its species qRRHO still falls back to the stored G, now with a warning. Every reaction of the default network conserves the number of molecules, so both modes give identical ΔG_rxn.
- **Stored enthalpies are reported.** In `'wb97mv'` mode `calculate_gas_thermo` returns the stored enthalpy and S = (H − G)/298.15 K instead of H = G, S = 0. G, ΔG_rxn, K_eq and every rate constant are unchanged.

## Verification

A baseline of the pre-refactor API was compared with the new code. It covered solution G, 324 rate-constant records, 4 protocol models, the batch reactor, the NMR catalog, water balance, fingerprints and σ.

- **Bit-identical:** thermochemistry, rate constants, NMR catalog, selectivity and σ.
- **Solver-level differences only:** reactor trajectories agree to ≤ 5·10⁻¹⁰ M. The analytic Jacobian slightly changes the Radau steps; the fast Marcus and Level 1 runs also no longer overflow in finite differences.
- **Legacy oracle:** `block5_legacy_baseline.json` (144 records, in the data folder, `tests/`) is still reproduced to < 10⁻¹².
