# File classification

| | |
|---|---|
| **Status** | Stable |
| **Classification** | PUBLIC |
| **Scope** | What may be published if the repository is made public again. Every file in the repository is listed below, by path or by folder. |
| **Rule** | This registry is the reference. Where the format allows it, the file also carries a `Classification` line in its header; if the two disagree, this registry wins. |

## Levels

| Level | Meaning |
|---|---|
| **PUBLIC** | Can be published as it is. |
| **CONFIDENTIAL** | Unpublished data of third parties, results derived from them, personal paths or internal project management. Must not be published. |
| **REVIEW** | Not decided yet. Treat as CONFIDENTIAL until it is reviewed. |

**Rules used.** Data derived from the Tank snapshot [PB] or from the raw lab spectra [NG] inherit their classification: CONFIDENTIAL. This covers results tables, stored notebook outputs, regression baselines and docs that quote computed numbers. Code that reads such data is PUBLIC unless it names private files or personal paths. Model parameters calibrated against published windows [G24] are PUBLIC.

**Where the label goes.** `.py`: module docstring. `.md`: header table (README, MODULES, the theory doc and `data/README.md` are kept unchanged on purpose: registry only). `.ipynb`: first markdown cell. `.json`: `_meta.classification` (the Tank snapshot is kept unchanged: registry only). `.csv` and other formats without a header: registry only.

## Sources

| ID | Source | Status |
|---|---|---|
| **Y** | Yeray Alcaraz Galván: code, documents, notebooks, curation of [G24] | Own work |
| **PB** | Peter Broqvist: Tank dataset API snapshot (`data/tank_api_snapshot.json`: ωB97M-V gas-phase DFT, MD solvation in EC, NMR shieldings, geometries), captured 2026-09-24; the BatteryAsTank bench protocol and reference model (`tmspa_hydrolysis_protocol.ipynb`) | Unpublished |
| **NG** | Neha Gogoi (Ångström Laboratory, Uppsala): raw JEOL NMR spectra recorded in 2022–2023. Not stored in the repository; read from `LAB_NMR_DIR`. | Unpublished |
| **G24** | Gogoi et al., *J. Phys. Chem. C* 2024, 128, 1654: shifts, control experiments, water series | Published |

## Registry

### Data

| Path | Level | Source | Reason |
|---|---|---|---|
| `data/tank_api_snapshot.json` | CONFIDENTIAL | PB | Unpublished computed dataset |
| `data/experimental_gogoi2024.json` | PUBLIC | G24, Y | Values read from the published paper; windows derived by Y |
| `data/README.md` | CONFIDENTIAL | Y, PB | Describes the snapshot and quotes one of its values (pure-EC solvation term, 0.202 eV); public once that value is removed |
| `tests/data/block5_legacy_baseline.json` | CONFIDENTIAL | PB | Regression results computed from the snapshot |
| `notebooks/results/01/` (all files) | CONFIDENTIAL | PB | Tables computed from the snapshot; `tank_trajectory_peter_reference.csv` reproduces PB's reference model |
| `notebooks/results/03/` (all files) | CONFIDENTIAL | PB | Tables computed from the snapshot |
| `notebooks/results/03/lab_shares_31P.csv`, `lab_shares_13C.csv`, `observations.csv`, `m0_intervals.csv`, `m1_intervals.csv`, `m1_region.csv`, `summary.csv`, `edges_hydrolysis_transfer.csv`, `plane_hydrolysis_transfer.csv` | CONFIDENTIAL | PB, NG | Also contain or depend on the raw lab spectra |
| `degradation_timeseries.csv`, `notebooks/degradation_timeseries.csv` | CONFIDENTIAL | PB | Simulation output computed from the snapshot |

### Notebooks

All notebooks are CONFIDENTIAL because their stored outputs show snapshot values. The code cells alone would be PUBLIC (except NB04 and the lab sections of the feasible-region notebook), so the notebooks can be published once their outputs are cleared.

| Path | Level | Source | Reason |
|---|---|---|---|
| `notebooks/01_multiscale_microkinetics_theory.ipynb` | CONFIDENTIAL | Y, PB, G24 | Outputs from the snapshot |
| `notebooks/02_operando_experimental_protocol.ipynb` | CONFIDENTIAL | Y, PB | Outputs from the snapshot; PB's protocol |
| `notebooks/03_experiment_plan.ipynb` | CONFIDENTIAL | Y, PB, G24 | Outputs from the snapshot |
| `notebooks/03_feasible_region.ipynb` | CONFIDENTIAL | Y, PB, NG, G24 | Outputs from the snapshot and the lab spectra; names NG |
| `notebooks/04_lab_nmr_data_overview.ipynb` | CONFIDENTIAL | Y, NG | Lab spectra, file names, personal paths |

### Documents

| Path | Level | Source | Reason |
|---|---|---|---|
| `README.md` | PUBLIC | Y, PB | Model descriptions and parameters only |
| `MODULES.md` | PUBLIC | Y | Code reference |
| `CLASSIFICATION.md` | PUBLIC | Y | This registry |
| `docs/theory-multiscale_microkinetics.md` | PUBLIC | Y, G24 | Theory; Level 1 parameters calibrated against G24 |
| `docs/guide-operando_experimental_protocol.md` | PUBLIC | Y, PB | Figure guide; parameter values only. PB's protocol is described, not reproduced |
| `docs/reference-api_migration.md` | PUBLIC | Y | Code reference |
| `docs/guide-experiment_plan.md` | CONFIDENTIAL | Y, PB | Quotes numbers computed in NB03 from the snapshot |
| `docs/methodology-feasible_region.md` | CONFIDENTIAL | Y, PB | Quotes numbers computed from the snapshot (NB03, Appendix A) |
| `docs/improvement_plan.md` | REVIEW | Y | Internal historical audit; to be reviewed |

### Code

| Path | Level | Source | Reason |
|---|---|---|---|
| `kinetics/` (all files except below) | PUBLIC | Y | Own code; reads data but contains none |
| `kinetics/data/lab_nmr.py` | CONFIDENTIAL | Y, NG | Default data path is a personal OneDrive folder; public once only `LAB_NMR_DIR` is used |
| `demo/` (all files) | PUBLIC | Y | Plotting code |
| `tests/test_models_and_engine.py`, `test_pipeline_integration.py`, `test_solvation_uncertainty.py` | PUBLIC | Y | Need the snapshot to run, but contain no data |
| `tests/test_barrier_models.py` | REVIEW | Y, PB | Hard-codes per-reaction barriers computed from the snapshot |
| `tests/test_fitting.py` | REVIEW | Y, PB, G24 | Hard-codes barrier bounds computed from the snapshot and G24 |
| `tests/test_lab_nmr.py` | CONFIDENTIAL | Y, NG | Names private lab files and acquisition dates |
| `scripts/sync_linear.py` | CONFIDENTIAL | Y | Internal project management (Linear) |
| `requirements.txt`, `.gitignore` | PUBLIC | Y | Configuration |

### Never committed

| Path | Reason |
|---|---|
| `.env` | API key (`LINEAR_API_KEY`); copy it between machines by hand |
| `.vscode/`, `.DS_Store`, `__pycache__/` | Local editor and system files |
