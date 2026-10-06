# File classification

| | |
|---|---|
| **Status** | Stable |
| **Classification** | PUBLIC |
| **Scope** | What is stored in this repository, and what is kept out of it. Every file of the repository is listed below, by path or by folder, with the reason. |
| **Rule** | This registry is the reference. Where the format allows it, the file also carries a `Classification` line in its header; if the two disagree, this registry wins. |
| **Revised** | 2026-10-06: the private files moved out of the repository (see "Where the private files are") |

## Levels

| Level | What |
|---|---|
| **PRIVATE** | The raw data and the audit documents. Never stored in the repository. |
| **PUBLIC** | Everything else: code, results, interpretation of the data, notebooks. |

**The rule.**

- **Private: the raw data.** The database (P. Broqvist's Tank dataset snapshot) and the raw lab spectra (N. Gogoi), together with the files that are copies or listings of them: tables that dump the snapshot, per-spectrum tables with raw file names, the curated observables file.
- **Private: the audit documents.** Internal audits, plans, methodology drafts and questions to collaborators.
- **Public: results and interpretation.** Anything computed from the data, and any way of showing it: fitted parameters, barriers, rate constants, tables, figures, and notebooks that print or plot values of the data. Code is public, including code rewritten from P. Broqvist's notebook, which is credited in its header.
- **Credit.** What is computed from the Tank dataset names P. Broqvist as its source; what is computed from the lab spectra names N. Gogoi.

## Sources

| ID | Source | Status |
|---|---|---|
| **Y** | Yeray Alcaraz Galván: code, documents, notebooks, curation of [G24] | Own work |
| **PB** | Peter Broqvist: Tank dataset API snapshot (ωB97M-V gas-phase DFT, MD solvation in EC, NMR shieldings, geometries), captured 2026-09-24; the BatteryAsTank bench protocol and reference model | Unpublished |
| **NG** | Neha Gogoi (Ångström Laboratory, Uppsala): raw JEOL NMR spectra recorded in 2022–2023 | Unpublished |
| **G24** | Gogoi et al., *J. Phys. Chem. C* 2024, 128, 1654: shifts, control experiments, water series | Published |

## Where the private files are

The private data folder is a folder of the thesis project, outside the repository. The code finds it through two
environment variables, set in the shell or in a local `.env` (git-ignored; copy `.env.example`). No path to either
folder is written in the code, the notebooks or the docs (`kinetics/paths.py` is the only reader).

| Variable | Holds |
|---|---|
| `ATOM_PRIVATE_DIR` | `data/`, `results/`, `tests/`, `scripts/` and `docs-internal/`, laid out as below |
| `LAB_NMR_DIR` | the raw lab spectra (JEOL `.jdf` and text exports) |

`scripts/check_public.py` (run by the pre-commit hook in `.githooks/`) refuses a commit that contains a path to
either folder, a private file name, a raw lab file name, or a copy of the Tank dataset in a data, code or text file
(any energy, solvation energy, shielding, coordinate or identifier, in the stored unit or a converted one; this
check needs the snapshot, so it runs where `ATOM_PRIVATE_DIR` is set). Notebooks may print and plot values; there the
guard refuses only a copied geometry or a dataset identifier. Enable it once with `git config core.hooksPath .githooks`.

## Registry

### PRIVATE (not in the repository)

| Where (under `ATOM_PRIVATE_DIR`) | Source | Reason | Read by |
|---|---|---|---|
| `data/tank_api_snapshot.json` | PB | Unpublished computed dataset | `kinetics/data/snapshot.py` |
| `data/lab_observables.json` | NG, Y | Area shares and sample histories curated from the raw spectra | `kinetics/data/observables.py` |
| `data/lab_sample_folders.json` | NG, Y | Names of the raw lab folders, files and titles | `kinetics/data/observables.py`, notebooks 03 and 04 |
| `results/01/` (`snapshot_fields`, `snapshot_inputs`, `species_database`, `species_table`, `tank_trajectory_peter_reference`) | PB | Tables that reproduce the snapshot or PB's reference model | written by notebook 01 |
| `results/03/lab_shares_31P.csv`, `lab_shares_13C.csv` | NG | Per-spectrum tables with raw file names, titles and acquisition metadata | written by notebook 03; read by `build_observables.py` |
| `results/degradation_timeseries.csv` and its notebook copy | PB | Simulation output computed from the snapshot, older runs | nothing |
| `tests/block5_legacy_baseline.json` | PB | Regression results computed from the snapshot | `tests/test_barrier_models.py` (skips without it) |
| `tests/test_lab_nmr.py` | NG, Y | Names raw lab files and acquisition dates | run with `ATOM_REPO_DIR` set |
| `scripts/sync_linear.py` | Y | Internal project management (Linear) | run by hand |
| `docs-internal/` | Y, PB, NG | Guide to the experiment plan, methodology of the feasible region, plan of the first fit, improvement plan (audit of the first notebook), follow-up pipeline overview, notes on the Tank snapshot, questions for N. Gogoi | – |
| the raw lab spectra (`LAB_NMR_DIR`) | NG | Unpublished measurements | `kinetics/data/lab_nmr.py` |

Four of these are results, not raw data, and could be public by the rule: the regression baseline, the two
`degradation_timeseries` files and `tank_trajectory_peter_reference.csv`. They are kept with the private files only
because nothing public needs them.

### PUBLIC (in the repository)

#### Data and results

| Path | Source | Reason |
|---|---|---|
| `data/experimental_gogoi2024.json` | G24, Y | Values read from the published paper; windows derived by Y |
| `data/README.md` | Y | Inventory and data model; no value of the Tank dataset |
| `notebooks/results/01/` (the six remaining tables) | Y, from PB | Derived: reaction thermodynamics, family barriers, rate constants, control experiments, scorecard, sensitivity |
| `notebooks/results/03/` (all except `lab_shares_*`) | Y, from PB, NG, G24 | Derived: barrier bounds, intervals, observations, design tables |
| `notebooks/results/05/` (all files) | Y, from PB, NG | Fit results: parameters, barriers, intervals, residuals, predictions, with the observed area shares behind each residual |

#### Notebooks

All notebooks are public, with their stored outputs: figures, derived tables, and tables that show values of the
data. Two things were taken out of notebook 04 (§5) because they are listings of the raw files, not results: the
output that listed raw lab file names (replaced by a note) and the table of questions to N. Gogoi (an internal
document, kept private).

| Path | Source | Note |
|---|---|---|
| `notebooks/01_multiscale_microkinetics_theory.ipynb` | Y, PB, G24 | Writes three tables to the private folder |
| `notebooks/02_operando_experimental_protocol.ipynb` | Y, PB | PB's protocol is described, not reproduced |
| `notebooks/03_experiment_plan.ipynb` | Y, PB, G24 | |
| `notebooks/03_feasible_region.ipynb` | Y, PB, NG, G24 | Reads raw spectra from `LAB_NMR_DIR`; writes `lab_shares_*` to the private folder |
| `notebooks/04_lab_nmr_data_overview.ipynb` | Y, NG | Spectra and sample-level tables; raw names read from `lab_sample_folders.json` |
| `notebooks/05_first_fit.ipynb`, `05_Alternative_first_fit_ladder.ipynb` | Y, PB, NG | Show the results of `scripts/run_fit.py` |

#### Documents

| Path | Source | Note |
|---|---|---|
| `README.md`, `MODULES.md`, `CLASSIFICATION.md` | Y, PB | Descriptions, parameters, this registry |
| `docs/theory-multiscale_microkinetics.md` | Y, G24 | Theory; Level 1 parameters calibrated against G24 |
| `docs/guide-operando_experimental_protocol.md` | Y, PB | Figure guide; PB's protocol is described, not reproduced |
| `docs/reference-api_migration.md` | Y | Code reference |

#### Code

| Path | Source | Note |
|---|---|---|
| `kinetics/` | Y | Own code. It reads private data but contains none and no path to them; `paths.py` is the only place that reads the environment variables |
| `demo/` | Y | Plotting code |
| `kinetics/spectroscopy/symmetry.py`, `spectra.py`, `fingerprints.py`, `demo/fingerprint_plots.py` | Y, adapted from PB | Own rewrite of the method of PB's notebook `tmspa_hydrolysis_protocol.ipynb` (shift catalog from shieldings, synthetic spectra, reaction fingerprints). Credited in each file |
| `kinetics/reactor/protocol.py`, the `peter_reference` model (`kinetics/microkinetics/parameters.py`, `models.py`) | Y, from PB | PB's bench protocol and reference model, re-implemented. Credited |
| `scripts/build_observables.py`, `scripts/run_fit.py`, `scripts/check_public.py` | Y | Builder of the observables file, runner of the fit, the guard |
| `tests/` (all `test_*.py` except `test_lab_nmr.py`, `run_all.py`, `synthetic_observables.py`) | Y | Need the private snapshot to run; contain no data. Barriers and standard errors hard-coded in `test_barrier_models`, `test_fitting` and `test_estimation` are derived values, credited to PB |
| `.env.example`, `.githooks/`, `.gitignore`, `requirements.txt` | Y | Configuration |

### Never committed

| Path | Reason |
|---|---|
| `.env` | The two private paths and the Linear API key (`LINEAR_API_KEY`); copy it between machines by hand |
| `.vscode/`, `.DS_Store`, `__pycache__/` | Local editor and system files |

## History

Files that were committed before 2026-10-06 remain in the git history of this repository, which is private. The history
holds the Tank snapshot, the lab observables, the old internal documents, and a copy of the text of G24. Before the
repository is ever made public, the history must be rewritten (`git filter-repo`, removing the files listed under
PRIVATE and `gogoi_2024_text.txt`); stopping to track them is not enough.
