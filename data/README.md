# Data inventory: `tank_api_snapshot.json`

An offline snapshot of Peter Broqvist's Tank dataset API (`http://127.0.0.1:8000/api`), captured on 2026-09-24. It is the only source of thermochemistry, solvation, NMR and geometry data in this repo. **All of it is computed; there is no experimental data in the repo.**

| Endpoint | Records | Content |
|---|---|---|
| `datasets` | 24 | Gas-phase DFT summary per species |
| `solvation` | 12 | MACE-OMol MD solvation energies in EC |
| `nmr` | 24 | Per-atom computed NMR shieldings, keyed by `dataset_uuid` |
| `structure` | 24 | Optimised geometry (elements + Cartesian coordinates in Å) |

## `datasets` (gas-phase DFT)

- **Level of theory:** ωB97M-V/def2-TZVPD with VV10 non-local dispersion, for all 24 species. Every species is neutral and closed-shell.
- **Energies (Hartree):** `energy_scf_eh`, `zpe_eh`, `enthalpy_eh`, `gibbs_eh`. G is at 298.15 K and 1 bar. The rigid-rotor/harmonic-oscillator entropy follows as S = (H − G)/T; for H₂O it is 188.6 J/(mol·K) against 188.8 J/(mol·K) experimentally.
- **Other properties:** `dipole_debye`, `polarizability_iso_ang3`, `homo_ev`, `lumo_ev`, `gap_ev`, `optical_gap_ev`.
- **Flags only, no data:**
  - `has_xps`, `has_xas`, `xas_edges`, `has_optical`, `n_xps_sites`: the spectra themselves were not captured.
  - `n_warnings`: the warning texts were not captured.
- **`n_modes` is a count, not the frequencies.** No vibrational frequencies, Hessians or moments of inertia are in the snapshot, so the `qRRHO` mode cannot be run with real data.
- **`n_modes` is below 3N−6 for five species:** TMSPA (123 of 126), BMSPA (89 of 90), HMDSO (74 of 75), TMSOdiEG (80 of 81) and c6h6f2li2o6-2 (59 of 60). All five have `n_warnings` ≥ 1. The likely cause is imaginary or discarded soft modes, which would bias `gibbs_eh`. This is not yet confirmed (TFM open question DFT-02).

**Species.**
- **Reaction network (11):** TMSPA, BMSPA, MMSPA, H3PO4, H2O, TMSOH, HMDSO (the network calls it `siloxyl`), EC, TMSOEG, TMSOdiEG, CO2.
- **Not wired into the network (13):** TMS (NMR reference), PH3, DMSO, DEC, DMC, DME, VC, propylene_carbonate, C2H4, bis_2_oxoethyl_oxalate, dilithium_4_carboxylatooxybutyl_carbonate, c6h6f2li2o6 and c6h6f2li2o6-2 (two entries with the same formula).
- **Not in the snapshot:** no fluorine/HF chemistry species (TMSF, HF, LiPF6, POF3, PF5) and no Li⁺.

## `solvation` (MACE-OMol MD in EC)

- **Species (12):** the 11 network species plus VC. The model is `MACE-OMol-extra_large` at 298.15 K. Each box has 14 EC and the solute; the pure-EC reference box has 15 EC.
- **ΔE_solv:** `delta_e_solv_kjmol` = E(solution) − (14/15)·E(pure EC) − E(gas solute), all MD averages. This is an energy, not a free energy: there is no solvation entropy.
- **`uncertainty_kjmol`:** the quadrature sum of `raw_metadata` `E_gas_std_eV`, `E_solution_std_eV` and `E_solvent_std_eV`.
  - The pure-EC term (0.202 eV) is one shared run and cancels in 2 → 2 reactions (see `kinetics/thermo/solvation_uncertainty.py`).
  - Whether these stds are per-frame deviations or standard errors of the mean is not documented (TFM open question SLV-01).
- **Also in `raw_metadata`:** `E_*_mean_eV`, `E_solvent_scaled_eV`, `solution_start_density_g_cm3`, `smiles`.

## `nmr` (computed shieldings)

- **Content:** per atom, `element`, `atom_index` (matching `structure`), absolute `isotropic_ppm` shielding and `anisotropy_ppm`. `parse_strategy` = `summary_table` for all 24, i.e. parsed from the program's shielding summary.
- **Coverage:** 211 H, 94 C, 68 O, 12 Si, 5 P and 4 F values. Li (and S in DMSO) are not included.
- **These are DFT results, not measurements.**
  - They are absolute shieldings per atom of a single static geometry. Chemically equivalent methyl H in TMS differ by about 0.01 ppm, whereas a solution spectrum would average them.
  - They include anisotropies, which solution NMR does not measure.
  - They exist for every dataset, including the Li salts.
- **Use:** `kinetics/spectroscopy/molecular_symmetry.py` averages them over equivalent nuclei and references them (TMS for ²⁹Si/¹³C/¹H, H₃PO₄ for ³¹P). They give peak positions only, with no concentrations or times, so they cannot be used to fit kinetics.
- **Accuracy against Gogoi et al. 2024:** ³¹P within about 2.6 ppm. ²⁹Si is systematically 3–6 ppm high (see the theory doc, Block 8.2).

## Experimental data: `experimental_gogoi2024.json`

The only experimental data in the repo come from Gogoi et al., *J. Phys. Chem. C* 2024, 128, 1654 (TFM BIB folder, `RXN - 02`). The measurements are in EC/DEC 1:1 without LiPF₆, and the paper tabulates no concentrations. The file holds four kinds of data:

- **`nmr_shifts`:** measured ³¹P (TMSPA, BMSPA, MMSPA, H₃PO₄) and approximate ²⁹Si shifts (xMSPA ≈ 25, TMSOH ≈ 15, HMDSO ≈ 7 ppm).
- **`barrier_constraints`:** ΔG‡ windows for R8 (1.25–1.36 eV), R4 (≥ 1.30 eV) and R5 (≤ 0.925 eV), derived from the reported observations.
- **`control_experiments`:** four simple mixtures (TMSOH in EC at RT and 80 °C; TMSPA + TMSOH at RT), each with a window on one observable. `kinetics/reactor/validation.py` re-simulates them for any model.
- **`water_series`:** ³¹P of 5 vol% TMSPA with 0.5, 1, 2 and 5 vol% water at RT (time after mixing not stated), and the 2 vol% sample heated to 80 °C. The windows on the P fractions are our reading of the text and figures (status `reading`). `evaluate_water_series` and `evaluate_heating_observation` test a model against them (notebook 03, §4).

The windows are derived in TFM *KIN - DRAFT - Level 1 formulation refined scaling relations - 260928*, §4.4. Entries marked `assumed` rest on a detection limit or reaction time the paper does not state. There is no quantitative constraint on hydrolysis (R1–R3); the water series is used only as a test of the model structure. Lab data from the Ångström group (Erik) are pending.
