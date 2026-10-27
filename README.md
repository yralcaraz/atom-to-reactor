# atom-to-reactor · branch `M1-steps`

Multiscale kinetics of the TMSPA water scavenger in ethylene carbonate: from DFT energies to reactor concentrations, fitted to the lab NMR spectra.

Master's thesis project, Department of Chemistry – Ångström Laboratory, Uppsala University.

> [!WARNING]
> **Research code, under active development.** No model is calibrated yet; APIs may change without notice.

## The branches

| Branch | What it holds |
|---|---|
| **`M1-steps`** (this one, the default) | The model M1 fitted to one lab tube at a time, each step tested on the tubes not used yet (eighth fit). Notebook `05_fitting_rev08_(M1_steps)` |
| `M5-full` | The fit of every tube at once: M1 ⊂ M5 ⊂ M5-first, each with Marcus barriers and with the capped BEP rule (seventh fit). Notebook `05_fitting_rev07_(M1_all fitted)` |
| `main` | The history up to the sixth fit: notebooks 01–05 rev05, the spectroscopy and experiment-design modules, the earlier stages of the runner. Closed at the tag `historico-261027` |

This branch keeps only what the M1 steps need. Anything removed is in `main`.

## What it does

```mermaid
flowchart LR
    DATA["Data<br/>DFT · MD solvation<br/>lab NMR observables"] --> THERMO["Thermochemistry<br/>G in EC/DEC, ΔG_rxn"]
    THERMO --> RATES["Rate constants<br/>barrier models + Eyring"]
    RATES --> REACTOR["Reactor<br/>each tube's history"]
    REACTOR --> FIT["Fit<br/>one tube at a time"]
    DATA -. "³¹P and ¹³C shares" .-> FIT
```

TMSPA (tris(trimethylsilyl) phosphate) reacts with water and with the solvent through a cascade of silyl transfers: R1–R9 and the six steps N10–N15 between the same species. Energies belong to species, so every cycle closes; every rate obeys detailed balance. The barriers are not computed: each class of reaction has an intrinsic barrier, and the fit determines it from the spectra.

## Quick start

Requires Python 3.10+ (developed on 3.11).

```bash
git clone https://github.com/yralcaraz/atom-to-reactor.git
cd atom-to-reactor
pip install -r requirements.txt
cp .env.example .env              # then set the two data folders (see Data)
```

| Task | Command |
|---|---|
| Tests | `python tests/run_all.py` |
| Check that the code still reproduces a stored fit (reads only) | `python scripts/check_stored_fits.py notebooks/results/fit8_M1-steps` |
| The registered models as they are | `python scripts/run_fit.py registered --results DIR` |
| The steps of the setup in `DIR/fit_config.json` | `python scripts/run_fit.py steps --results DIR` |
| The notebook | `jupyter lab notebooks/` |

`--results` is required: every job is written there, and a folder that already holds results returns its fits as they are.

## Layout

```text
atom-to-reactor/
├── kinetics/            the package
│   ├── data/            snapshot, species, network, lab observables, solvation medium
│   ├── thermo/          G in solution, solvation σ, ΔG_rxn
│   ├── microkinetics/   barrier models, rate constants, model registry
│   ├── reactor/         mass-action engine, batch and protocol reactors, sample histories
│   └── fitting/         residuals, candidate structures, fits and profiles, stored fits run forward, run setup
├── demo/                figures and display tables of notebook rev08
├── notebooks/           05_fitting_rev08_(M1_steps)
├── scripts/             run_fit.py (registered, steps) · check_stored_fits.py · build_observables.py
├── tests/               test scripts (run_all.py runs them)
├── data/                reference data and the data inventory (data/README.md)
├── docs/                theory and references
└── MODULES.md           module reference, written for main: it also describes modules that are only there
```

## Data

| File | Content | Origin |
|---|---|---|
| `data/experimental_gogoi2024.json` | Measured shifts, barrier windows, control experiments | Gogoi et al., *J. Phys. Chem. C* 2024, 128, 1654 |
| `tank_api_snapshot_2026-09-24.json` | Gas-phase DFT (ωB97M-V/def2-TZVPD) and MACE-OMol MD solvation in EC; the default, used by the tests | Computed; offline snapshot of the Tank dataset API |
| `tank_api_snapshot_2026-10-08.json` | The same for 27 species, with MD solvation in DEC; chosen by the fits in their `fit_config.json` | Computed; later capture of the same API |
| `lab_observables.json` | Area shares, uncertainties and histories of the NMR tubes | NMR data from experiments; built by `scripts/build_observables.py` |
| NMR spectra | JEOL Delta `.jdf` files | NMR data from experiments |

The last four are not in the repository. They are read from the folders set by `ATOM_DATA_DIR` and `LAB_NMR_DIR` (see [.env.example](.env.example) and [kinetics/paths.py](kinetics/paths.py)). The fit results under `notebooks/results/` are derived from the lab data and are not in the repository either, except the early ones kept in `notebooks/results/01`, `03` and `05`.

## Context

- **Project:** Multiscale Modelling of EMC Transesterification as a Degradation Pathway in Li-ion Battery Electrolytes, Master's thesis
- **Author:** Yeray Alcaraz Galván
- **Supervision:** Prof. Peter Broqvist
- **Affiliation:** Department of Chemistry - Ångström Laboratory, Uppsala University, Sweden

No licence file is included yet.
