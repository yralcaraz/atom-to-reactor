# Data inventory

| File | Kind | Content | Location | Loaded by |
|---|---|---|---|---|
| `experimental_gogoi2024.json` | Measured, from the literature | NMR shifts, barrier windows, control experiments, water series (Gogoi et al. 2024) | this folder | `kinetics/data/experimental.py` |
| `tank_api_snapshot.json` | Computed | Gas-phase DFT, MD solvation in EC, NMR shieldings, geometries (Tank dataset) | data folder | `kinetics/data/snapshot.py`, `species.py` |
| `lab_observables.json` | Measured | Area shares and sample histories from the experimental NMR spectra | data folder | `kinetics/data/observables.py` |
| NMR spectra (`.jdf`) | Measured | NMR data from experiments (JEOL Delta) | `LAB_NMR_DIR` | `kinetics/data/lab_nmr.py` |

## Data folder

The data files are read from folders set with two environment variables (see `.env.example` and `kinetics/paths.py`):

| Variable | Folder |
|---|---|
| `ATOM_DATA_DIR` | `data/` (the snapshot, `lab_observables.json`, `lab_sample_folders.json`), `results/` (tables written by notebooks 01 and 03) and `tests/` (the regression baseline) |
| `LAB_NMR_DIR` | the NMR spectra |

The species database, the network thermochemistry and every test are built from the snapshot, so the numerical
pipeline needs `ATOM_DATA_DIR`. The DFT and MD data are from P. Broqvist's Tank dataset.

## Snapshot: `tank_api_snapshot.json`

An offline snapshot of the Tank dataset API, captured on 2026-09-24. It is the only source of thermochemistry,
solvation, NMR and geometry data in this repo. **All of it is computed.**

| Endpoint | Content |
|---|---|
| `datasets` | Gas-phase DFT summary per species (ωB97M-V/def2-TZVPD, neutral closed-shell species) |
| `solvation` | MACE-OMol MD solvation energies in EC |
| `nmr` | Per-atom computed NMR shieldings, keyed by `dataset_uuid` |
| `structure` | Optimised geometry (elements and Cartesian coordinates in Å) |

What to keep in mind when using it:

- **Energies** are in Hartree (`energy_scf_eh`, `zpe_eh`, `enthalpy_eh`, `gibbs_eh`); G is at 298.15 K and 1 bar in the rigid-rotor/harmonic-oscillator approximation.
- **`n_modes` is a count, not the frequencies.** No vibrational frequencies or Hessians are stored, so the `qRRHO` mode cannot be run with real data.
- **Solvation** `delta_e_solv_kjmol` is an energy, not a free energy: there is no solvation entropy. `uncertainty_kjmol` combines the standard errors of the three mean MD energies, and the pure-EC reference run is shared, so its term cancels in 2 → 2 reactions (see `kinetics/thermo/uncertainty.py`).
- **Shieldings** are absolute values of a single static geometry. They are DFT results, not measurements: chemically equivalent nuclei are not averaged and anisotropies are included. `kinetics/spectroscopy/symmetry.py` averages and references them (TMS for ²⁹Si/¹³C/¹H, H₃PO₄ for ³¹P). They give peak positions only, so they cannot be used to fit kinetics.
- **Species:** 11 are wired into the reaction network (TMSPA, BMSPA, MMSPA, H3PO4, H2O, TMSOH, HMDSO, EC, TMSOEG, TMSOdiEG, CO2); the others are references and solvents. There is no fluorine/HF chemistry and no Li⁺.

## Lab observables: `lab_observables.json`

Built by `scripts/build_observables.py` from the area shares of the ³¹P and ¹³C spectra. Per sample it holds the
recipe, what is known of its temperature history, and the shares with two measured uncertainties (noise, baseline
spread). The time from mixing to the first spectrum is unknown for every sample.

## Measured data: `experimental_gogoi2024.json`

The only experimental data in the repo come from Gogoi et al., *J. Phys. Chem. C* 2024, 128, 1654 (TFM BIB folder, `RXN - 02`). The measurements are in EC/DEC 1:1 without LiPF₆, and the paper tabulates no concentrations. The file holds four kinds of data:

- **`nmr_shifts`:** measured shifts with the figure they come from. ³¹P: TMSPA, BMSPA, MMSPA, H₃PO₄. ²⁹Si: xMSPA ≈ 25 (one signal for the three silyl phosphates), TMSOH 14.26–15, TMS-EG 18.01, HMDSO ≈ 7, and an unassigned line at 19.08 ppm. ¹³C and ¹H of the TMS methyl groups: TMSOH, HMDSO and TMS-EG (¹³C); the silyl phosphates and TMSOH/TMS-EG/HMDSO as two signals (¹H). Each entry names its `site` (P, Si or Si-CH3) and the `peak` it belongs to; species reported as one signal share a peak. `kinetics/fitting/readouts.py` derives the readouts of the experiment design from these entries. TMSOdiEG and the methyl ¹³C of the silyl phosphates were never measured.
- **`barrier_constraints`:** ΔG‡ windows for R8 (1.25–1.36 eV), R4 (≥ 1.30 eV) and R5 (≤ 0.925 eV), derived from the reported observations.
- **`control_experiments`:** four simple mixtures (TMSOH in EC at RT and 80 °C; TMSPA + TMSOH at RT), each with a window on one observable. `kinetics/reactor/validation.py` re-simulates them for any model.
- **`water_series`:** ³¹P of 5 vol% TMSPA with 0.5, 1, 2 and 5 vol% water at RT (time after mixing not stated), and the 2 vol% sample heated to 80 °C. The windows on the P fractions are our reading of the text and figures (status `reading`). `evaluate_water_series` and `evaluate_heating_observation` test a model against them (notebook 03, §4).

The windows are derived in TFM *KIN - DRAFT - Level 1 formulation refined scaling relations - 260928*, §4.4. Entries marked `assumed` rest on a detection limit or reaction time the paper does not state. There is no quantitative constraint on hydrolysis (R1–R3); the water series is used only as a test of the model structure.
