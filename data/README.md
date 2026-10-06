# Data inventory

| File | Kind | Content | In the repository | Loaded by |
|---|---|---|---|---|
| `experimental_gogoi2024.json` | Measured, published | NMR shifts, barrier windows, control experiments, water series (Gogoi et al. 2024) | yes | `kinetics/data/experimental.py` |
| `tank_api_snapshot.json` | Computed, unpublished | Gas-phase DFT, MD solvation in EC, NMR shieldings, geometries (P. Broqvist, Tank dataset) | **no: private data folder** | `kinetics/data/snapshot.py`, `species.py` |
| `lab_observables.json` | Measured, unpublished | Curated area shares and sample histories of the lab NMR spectra (N. Gogoi, 2022–2023) | **no: private data folder** | `kinetics/data/observables.py` |
| raw lab spectra (`.jdf`) | Measured, unpublished | JEOL Delta spectra (N. Gogoi) | **no: `LAB_NMR_DIR`** | `kinetics/data/lab_nmr.py` |

## Private data

The Tank dataset and the lab spectra are unpublished work of other people and are not redistributed. The code
reads them from folders that each user points to with two environment variables (see `.env.example` and
`kinetics/paths.py`):

| Variable | Folder |
|---|---|
| `ATOM_PRIVATE_DIR` | `data/` (the snapshot, `lab_observables.json`, `lab_sample_folders.json`), `results/` (tables that reproduce private inputs) and `tests/` (the regression baseline) |
| `LAB_NMR_DIR` | the raw lab spectra |

Without them the package imports, but nothing that needs the snapshot runs: the species database, the network
thermochemistry and every test are built from it. What this repository shows without the private data is the code,
the fit results (`notebooks/results/`) and the notebooks with their figures. To obtain the Tank snapshot, ask
Prof. Peter Broqvist; for the lab spectra, Neha Gogoi (Ångström Laboratory, Uppsala).

**Credit.** All DFT and MD data are from Peter Broqvist's Tank dataset (unpublished). Results in this repository that
derive from them (rate constants, barriers, fitted parameters, figures) are computed from his data and say so in their
`Source` line.

## Snapshot: `tank_api_snapshot.json` (private)

An offline snapshot of the Tank dataset API, captured on 2026-09-24. It is the only source of thermochemistry,
solvation, NMR and geometry data in this repo. **All of it is computed.**

| Endpoint | Content |
|---|---|
| `datasets` | Gas-phase DFT summary per species (ωB97M-V/def2-TZVPD, neutral closed-shell species) |
| `solvation` | MACE-OMol MD solvation energies in EC |
| `nmr` | Per-atom computed NMR shieldings, keyed by `dataset_uuid` |
| `structure` | Optimised geometry (elements and Cartesian coordinates in Å) |

What to keep in mind when using it (the values are in the private file):

- **Energies** are in Hartree (`energy_scf_eh`, `zpe_eh`, `enthalpy_eh`, `gibbs_eh`); G is at 298.15 K and 1 bar in the rigid-rotor/harmonic-oscillator approximation.
- **`n_modes` is a count, not the frequencies.** No vibrational frequencies or Hessians are stored, so the `qRRHO` mode cannot be run with real data.
- **Solvation** `delta_e_solv_kjmol` is an energy, not a free energy: there is no solvation entropy. `uncertainty_kjmol` combines the standard errors of the three mean MD energies, and the pure-EC reference run is shared, so its term cancels in 2 → 2 reactions (see `kinetics/thermo/uncertainty.py`).
- **Shieldings** are absolute values of a single static geometry. They are DFT results, not measurements: chemically equivalent nuclei are not averaged and anisotropies are included. `kinetics/spectroscopy/symmetry.py` averages and references them (TMS for ²⁹Si/¹³C/¹H, H₃PO₄ for ³¹P). They give peak positions only, so they cannot be used to fit kinetics.
- **Species:** 11 are wired into the reaction network (TMSPA, BMSPA, MMSPA, H3PO4, H2O, TMSOH, HMDSO, EC, TMSOEG, TMSOdiEG, CO2); the others are references and solvents. There is no fluorine/HF chemistry and no Li⁺.

## Lab observables: `lab_observables.json` (private)

Built by `scripts/build_observables.py` from the area shares of the lab ³¹P and ¹³C spectra. Per sample it holds the
recipe, what is known of its temperature history, and the shares with two measured uncertainties (noise, baseline
spread). The time from mixing to the first spectrum is unknown for every sample.

## Measured data: `experimental_gogoi2024.json`

The only experimental data in the repo come from Gogoi et al., *J. Phys. Chem. C* 2024, 128, 1654 (TFM BIB folder, `RXN - 02`). The measurements are in EC/DEC 1:1 without LiPF₆, and the paper tabulates no concentrations. The file holds four kinds of data:

- **`nmr_shifts`:** measured shifts with the figure they come from. ³¹P: TMSPA, BMSPA, MMSPA, H₃PO₄. ²⁹Si: xMSPA ≈ 25 (one signal for the three silyl phosphates), TMSOH 14.26–15, TMS-EG 18.01, HMDSO ≈ 7, and an unassigned line at 19.08 ppm. ¹³C and ¹H of the TMS methyl groups: TMSOH, HMDSO and TMS-EG (¹³C); the silyl phosphates and TMSOH/TMS-EG/HMDSO as two signals (¹H). Each entry names its `site` (P, Si or Si-CH3) and the `peak` it belongs to; species reported as one signal share a peak. `kinetics/fitting/readouts.py` derives the readouts of the experiment design from these entries. TMSOdiEG and the methyl ¹³C of the silyl phosphates were never measured.
- **`barrier_constraints`:** ΔG‡ windows for R8 (1.25–1.36 eV), R4 (≥ 1.30 eV) and R5 (≤ 0.925 eV), derived from the reported observations.
- **`control_experiments`:** four simple mixtures (TMSOH in EC at RT and 80 °C; TMSPA + TMSOH at RT), each with a window on one observable. `kinetics/reactor/validation.py` re-simulates them for any model.
- **`water_series`:** ³¹P of 5 vol% TMSPA with 0.5, 1, 2 and 5 vol% water at RT (time after mixing not stated), and the 2 vol% sample heated to 80 °C. The windows on the P fractions are our reading of the text and figures (status `reading`). `evaluate_water_series` and `evaluate_heating_observation` test a model against them (notebook 03, §4).

The windows are derived in TFM *KIN - DRAFT - Level 1 formulation refined scaling relations - 260928*, §4.4. Entries marked `assumed` rest on a detection limit or reaction time the paper does not state. There is no quantitative constraint on hydrolysis (R1–R3); the water series is used only as a test of the model structure.
