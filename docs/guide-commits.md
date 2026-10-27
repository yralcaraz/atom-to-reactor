# Guide to the commits

| | |
|---|---|
| **Status** | Living document: one entry per change set, newest last |
| **Source** | Y. Alcaraz Galván |
| **Scope** | What each change set modifies against the commit before it, and why, in a few lines. Reasons in depth are in the reference documents each entry links |

Each entry is titled `date - title`, and says which commit it is (or "not committed yet"), the commit it is compared with, what changed and how it was checked.

---

## 2026-10-20 - Fourth fit review and fit improvements reference

**Commit** `d1125a2` · **against** `412484f` (fourth fit setup)

- **Notebook `05_fitting_rev04_(EC_DEC)`**: section 2.2 (concentrations in time) split into the three readings (1 day, 7 days, free ages). Clock from the addition of TMSPa, with the preparation stage drawn from −36 h, one vertical line per NMR round, and a table of model against NMR in mM under each figure. The figures run the fourth fit's parameters with 8 h heating steps (`heated_hold_h`), which no fit used; the fifth fit replaces this (next entry).
- **`demo/fit_plots.py`**: `plot_tube_concentrations` and `show_tube_comparison` for those figures and tables.
- **`kinetics/fitting/residuals.py`**: the override `heated_hold_h` (each heating step held a given time and ended before the next spectrum).
- **`docs/reference-fit_improvements.md`**: the untracked `plan-fourth_fit.md`, renamed. Section 1 is the fourth fit with its outcome; section 2 is the plan of the fifth fit.

## 2026-10-20 - Fifth fit setup: 8 h heating steps of the 2 % tube, reactions N10–N15

**Commit** "Fifth fit review and sixth fit setup: one rule for every reaction, documents" (2026-10-22, together with the two entries below) · **against** `d1125a2`

Why: [reference-fitting_improvements.md](reference-fitting_improvements.md) §2 (then `reference-fit_improvements.md`). The fit runs the bench protocol for the heated 2 % tube (8 h per step, placed on the file clock) and the six further steps of P. Broqvist's generated network, without new free parameters. Nothing changes for a results folder whose `fit_config.json` does not ask for it: the fourth fit is reproduced exactly.

**Code**

| File | Change |
|---|---|
| `kinetics/data/network.py` | `EXTRA_REACTIONS` (N10–N15, each a combination of R1–R9). Network variants `R1-R9`, `R1-R9+N10-N15` and `R1-R9+N10-N15, Peter barriers`, with an active one per process (`build_network`, `get_active_network`, `set_active_network_variant`), as for the solvation medium. N12 and N15 carry `follows` (R2, R1); N10, N11, N13 and N14 are of class `exchange` |
| `kinetics/data/__init__.py` | Exports the above |
| `kinetics/fitting/candidates.py` | A structure builds on the active network; a step that follows a reaction takes its barrier class (N15 → `hydrolysis_R1` or `hydrolysis`, N12 → `hydrolysis_R23` or `hydrolysis`). `FIXED_BARRIERS_EV = {'exchange': 1.15}`: set in every structure, never a parameter. The model always carries its network. Search box of `g hydrolysis_R1` 0.80–1.80 eV (was 1.70): N15 makes g_R1 rise about 0.05 eV and the fourth fit's 7 day reading sat at 1.635. `tabulate_reaction_barriers` lists every reaction of the structure's network |
| `kinetics/fitting/residuals.py` | Run-wide changes to the samples: `set_sample_overrides`, `get_sample_overrides`, `merge_sample_overrides` (the overrides of a call go on top, key by key). New override keys `heated_episodes` (explicit steps `[start h, T °C, hours]` on the clock of the first spectrum, replacing the recorded heating; `None` keeps it), `added_episodes` (added to it) and `cooling_h` (at least that many hours at room temperature before every spectrum: `_cool_before_spectra` moves a step earlier, never over the one before). `build_sample_history` places a step that starts before the first spectrum inside the age, cut at mixing |
| `kinetics/fitting/config.py` | `fit_config.json` keys `network`, `sample_overrides`, `sensitivities` (defaults: `R1-R9`, none, the fourth fit's list). `apply_fit_config` sets the network and the changes to the samples. `check_fit_setup`: refuses to run a stored fit forward with another network, or with other run-wide changes when the fit records them |
| `kinetics/fitting/reporting.py` | `evaluate_fit` and `calculate_parameter_correlation` call `check_fit_setup` instead of `check_fit_medium` |
| `kinetics/fitting/__init__.py` | Exports `check_fit_setup` and the three sample-override functions |
| `scripts/run_fit.py` | Every fit records `network`, `sample_overrides` and `variant_overrides`. New sensitivities `heating as recorded`, `30 C before NMR 1` (8 h at 30 °C ending 1 h before NMR 1), `cooling 3 h` (every step of the 2 % tube at least 3 h before the next round) and `Peter barriers` (a variant may name its own network, used through `network_variant`). The variants of a run come from the setup (`SENSITIVITY_NAMES`). Class `exchange` anchored at 22.5 °C for the ΔS‡ band. Docstring and setup log line |
| `tests/test_medium.py` | `test_network_variants` (combinations hold, freed shifts included; classes per structure; no new parameter; fixed barrier), `test_heating_before_the_first_spectrum` (with the cooling rule), `test_fit_config_network_and_samples` |

**Documents**: `docs/reference-fitting_system.md` §2.5 and new §9; `MODULES.md` §2.3, §7.6 and §7.10; `docs/reference-fit_improvements.md` §2 (status, the checks below, cost).

**Not in git** (lab-derived): `notebooks/results/05_fit5/fit_config.json` (network `R1-R9+N10-N15`, the six steps of the 2 % tube, the eleven sensitivities, the fourth fit's recovery setup, a `why` list) and `notebooks/results/05_fit5/run_chain.sh` (the fourth fit's chain, renamed, with 14 workers instead of 12).

**Checks**

- Test suite (`python tests/run_all.py`): the nine test files pass (4.5 min).
- The fourth fit with its own setup: χ² 50.1352 at 1 day and 39.6541 with free ages, as stored.
- The fourth fit's 1 day parameters under the fifth setup (no fitting): base heating alone χ² 90.7; network alone 1495; both 1526; both with Peter barriers 1525. One added step at a time: N10, N11, N13, N14 change nothing (50.1–50.3), N12 gives 57.8, N15 gives 2204 (it starts the reaction about seven times faster than R1 at the same g).
- Heating history of the 2 % tube under the fifth setup: 30 and 40 °C from 13.3 h after the first spectrum, then 50, 60, 70 and 80 °C at 63.0, 87.0, 109.9 and 136.0 h, 8 h each; NMR 6 1.5 h after the end of the 80 °C step.
- Cost of one evaluation (search solver, one process), mean over 24 points spread over the search box: 1 day 0.38 s (fourth fit) → 0.40 s (8 h steps only) → 0.43 s (fifth fit); free ages 0.43 → 0.45 → 0.51 s. At the fourth fit's optimum: 0.62 → 0.79 s (1 day), 0.72 → 0.76 s (free ages). No simulation ran out of its step budget in any of the 150 evaluations. Throughput of the objective over 192 points of the box: 22.5 evaluations per second with 12 workers, 23.9 with 16 (each worker about 100 MB). The chain runs 14 workers; it should take about 15–17 h (fourth fit: 14.5 h with 12).
- Smoke run of the runner (tiny search, scratch folder): a main fit and the three new sensitivities ran; each fit records its network, the run's changes to the samples and its variant's overrides, and its barrier table lists the fifteen reactions. `evaluate_fit` runs the `heating as recorded` fit with that variant's overrides and refuses the `Peter barriers` fit until its network is made active (`set_active_network_variant('R1-R9+N10-N15, Peter barriers')`).
- Cooling before each NMR round, the fourth fit's 1 day parameters with the 8 h steps (R1–R9): χ² 154 with 0 h, 91 with the declared 1 h, 76 with 3 h. After a fast cool-down the ³¹P shares relax for about 2 h (R3), so the declared margin is not neutral.

**Effects on what exists**: the fourth fit and its notebooks read and evaluate as before. Notebooks rev03 and rev04 draw the search box of `g hydrolysis_R1` from `BARRIER_BOX_EV`, so, run again, that track reaches 1.80 eV.

## 2026-10-22 - Sixth fit setup: one rule for every reaction, energies within 3 σ, declared water

**Commit** "Fifth fit review and sixth fit setup: one rule for every reaction, documents" (2026-10-22, with the entries before and after) · **against** `d1125a2` plus the fifth fit's changes (entry above)

Why: [reference-fit_results.md](reference-fit_results.md) §6.3, the decisions of 2026-10-22. Nothing changes for a results folder whose `fit_config.json` does not ask for it.

**Code**

| File | Change |
|---|---|
| `kinetics/fitting/candidates.py` | Structures `M3-first` (the first step of TMSPA is its own class in hydrolysis and in transfer: `SPLIT_FIRST_STEP`, classes `transfer_R5` and `transfer_R67`, box 0.60–1.50 eV) and `M5-first` (the same, with the energies of every independent reaction free: `FREED_INDEPENDENT` = R1–R4, R8, R9). Run-wide energy box of ±k standard errors (`set_energy_box_sigma`, `get_energy_box_sigma`); the Gaussian constraint stays inside it. `build_fitted_model` names the freed reactions of the structure |
| `kinetics/fitting/residuals.py` | Override key `unknown_M`: replaces the amounts at mixing that are fit parameters (`{}`: none) |
| `kinetics/fitting/config.py` | `fit_config.json` keys `energy_box_sigma` (applied by `apply_fit_config`), and for the runner `structures`, `effort`, `start_from`, `night` |
| `scripts/run_fit.py` | Main and extended structures, search effort (per structure, or per structure and reading with the key `structure\|scenario`), an earlier results folder whose fits start the search (`earlier_run_starts`: embedded in each structure, only as starting points), and the parts of the night that run (profiles of the listed structures, optionally only in listed readings `{structure: [readings]}` and only where they fit; indistinguishability, tolerance check and trace refits optional) come from the setup. Recovery cases can be replaced by the setup (`recovery.cases`). New sensitivity `water free`. Every fit records `energy_box_sigma` |

**Which models run** (2026-10-22, after the review of §6.3.12): `structures.readings` in the setup and `READINGS` in `scripts/run_fit.py` let a structure run in some readings only. M3 is left out (it failed in the fifth fit, and every change of the sixth restricts it); M3-split and M3-first run at 1 day and free ages only (nested in M5-first at 1 h). The time saved goes to a second recovery case, `C6-1day` (M5-first at 1 day). `start_when_plugged.sh` builds its own smoke folder in `/tmp/fit6_smoke` and runs every stage before the chain.

**Names**: the two structures were first called `M3-R1R5` and `M3-R1R5-G6`; renamed `M3-first` and `M5-first` by the rule MX-suffix (X: which energies are free; suffix: how the barriers are grouped), stated in `candidates.py`, `reference-models_summary.md` §2 and `reference-fitting_system.md` §2.5a.

**Documents**: `docs/reference-fitting_system.md` revised (status, objective, parameter table, model names, selection rule, runner stages, search effort, new §10); `docs/reference-models_summary.md` §2 (naming rule, M3-first, M5-first, energy box); `MODULES.md` (structures, extended stage).

**Not in git** (lab-derived): `notebooks/results/05_fit6/fit_config.json` (readings 1 h, 1 day, free; box ±3 σ; main structure `M5-first` in every reading, `M3-first` and `M3-split` at 1 day and free ages, M3 left out, M0 to M1-split as the ladder; water of "TMSPa alone" declared at 65.6 mM; full search for M3-first and M5-first, reduced for M3-split and for M5-first at 1 h, all started also from the fifth fit; profiles of the main structure and of M3-first at 1 day; seven sensitivities: heating as recorded, covariance × 0.5 and × 2, water free, Peter barriers, floor × 0 and × 2; two recovery cases of M5-first, free ages and 1 day, with the reduced search; a `why` list), `run_chain.sh` (the fifth fit's chain) and `start_when_plugged.sh` (waits for mains power, smoke run of every stage in `/tmp/fit6_smoke`, then the chain). **Not started**: ready to launch.

**Checks**

- Test suite: `test_fitting`, `test_estimation`, `test_medium`, `test_reporting` pass.
- Under the sixth setup, the fifth fit's parameters (embedded where needed) give χ² 58.9 at 1 day (57.9 in the fifth fit: the declared water adds 1.0) and 48.2 with free ages; `M3-first` the same; `M5-first` 61.9 and 50.9 (R8 and R9 held at their computed values cost 3 through their correlation with R1–R4). One evaluation 0.74 s (1 day), 0.94 s (free ages).
- Energy boxes of the main structure: R1 ±0.311, R2 ±0.281, R3 ±0.242, R4 ±0.360, R8 ±0.310, R9 ±0.318 eV. Barrier classes: hydrolysis_R1, hydrolysis_R23, transfer_R5, transfer_R67, condensation, solvent_attack.
- Smoke run with a tiny search (scratch folder): registered, baselines and the extended fits of M1-split, M3 and M3-split at 1 day ran (M3-split χ² 58.85, fits). Stopped there because the laptop was on battery: the new structures in the runner, the night, traces and predictions stages are not smoke-tested yet. Static check: the search each structure and reading gets, the recovery case and the starts taken from the fifth fit (3 at 1 h, 1 at 1 day per structure) are as set.

**Effects on what exists**: none for fits 1–5; their configurations have no new keys and read as before.

## 2026-10-22 - Fifth fit review and documents

**Commit** "Fifth fit review and sixth fit setup: one rule for every reaction, documents" (2026-10-22, with the two entries above) · **against** `d1125a2`

The documents and notebooks of the fifth fit's review and of the sixth fit's setup. No code beyond the two entries above.

| File | Change |
|---|---|
| `docs/reference-fit_results.md` | New. Review of the fifth fit: outcome, what each tube constrains, the "all or nothing" switch of the model, the contradiction between "TMSPa alone" and the 0.5 % tube, reaction energies against the computed ones, linear and preferential mixing, experiments and computing options. §6.3: the questions, decisions and limitations of the sixth fit (the error of every computed number through the scales, per-reaction offsets, MD against COSMO-RS, one rule for every reaction, the ±3 σ box, ages, water and purity from the literature, the weight of the 2 % tube, setup, cost, open questions) |
| `docs/reference-fitting_improvements.md` | Renamed from `reference-fit_improvements.md` (Yeray), with the fifth fit's status |
| `docs/reference-models_summary.md` | Renamed from `reference-fit_models.md` (Yeray), content unchanged |
| `docs/guide-commits.md` | New: this log |
| `docs/reference-fitting_system.md`, `MODULES.md`, `README.md` | Links to the renamed documents; the README lists the new ones |
| `notebooks/05_fitting_rev04_(EC_DEC).ipynb` | Yeray's edits during the review of the fourth fit |
| `notebooks/05_fitting_rev05_(EC_DEC) copy.ipynb` | New: Yeray's working copy for the fifth fit (reads `notebooks/results/05_fit5/`); header and prose still partly those of rev04 |

**Not in git** (lab-derived, as every results folder): `notebooks/results/05_fit2/` to `05_fit6/`.

**Checks**: links to the renamed documents resolved in every document of `docs/`, `MODULES.md` and `README.md` (the notebooks still cite the older `plan-fit_improvement.md`).

## 2026-10-22 - Seventh fit setup: near-zero ages, M5, P. Broqvist's barrier rule

**Not committed yet** (branch `fit/near-zero-ages`, worktree `/home/yeray/projects/atom-to-reactor-fit7`) · **against** `ee6ab2a` (sixth fit ready to launch)

Why: meeting of 2026-10-22, "Simple M1 fit at near-zero ages", and its follow-up the same evening. Go back to a small fit, take the time from adding TMSPA to the first NMR as 0, 5 and 10 min, compare Marcus with P. Broqvist's rule, and add one model between M1 and M5-first. Nothing changes for a results folder whose `fit_config.json` does not ask for it; the sixth fit (`05_fit6`) is untouched.

**Code**

| File | Change |
|---|---|
| `kinetics/fitting/residuals.py` | Readings `0 min`, `5 min`, `10 min` in `SCENARIO_AGES_H` and `SCENARIOS`. `DEFAULT_SCENARIOS`: the four readings of a folder without `fit_config.json` |
| `kinetics/fitting/config.py` | The default setup takes `DEFAULT_SCENARIOS` (was every scenario), so the new readings enter only a run that lists them |
| `kinetics/fitting/candidates.py` | Structure `M5`: the four barriers of M1 with the energies of every independent reaction free (R1–R4, R8, R9), nested between M1 and M5-first. `FitStructure.shape_params` and `with_peter_bep`: suffix `+BEP`, capped BEP with α = 1 (`PETER_BEP_ALPHA`) in every class, the exchange steps at 1.15 eV included; `get_structure` reads the suffix like `+W` |
| `scripts/run_fit.py` | `M5` in `EFFORT`, `PARENTS` (M1 → M5 → M5-first) and `SIZE_ORDER`; `+BEP` forms nest like their Marcus forms and also start from them. The near-zero readings first in `SCENARIO_ORDER` (5 and 10 min start from 0 min) and in the common ages of the selection rule. The summary tables hold every structure of the run (`SIZE_ORDER` extended from the setup). New option `--structures`: one call of extended, recovery, night, traces and predictions limited to some structures of the setup (`ONLY`). Night, traces and predictions analyse every structure of `night.structures` (`ANALYSED`; default the main structure, as before), in the readings the setup allows (new `night` keys `refined`, `loo_fitting_only`, `loo_effort`, `sensitivity_readings`, `recovery_required`; `tolerance_check` also as `{structure: [reading, [two parameters]]}`), and skip a reading the run does not have instead of raising (the sixth fit's smoke failure). Recovery truth `"fit"`: the run's own fit, kept in `recovery/<case>_truth.json`. Covariance sensitivities skipped without freed energies. `structure` column in `sensitivities.csv`, `trace_sensitivity.csv`, `tolerance_check.csv`; `sample_trajectories.csv` from 36 s, with "TMSPa alone" and tube A. Docstring |
| `tests/test_estimation.py` | `test_seventh_fit_setup`: near-zero readings opt-in and simulated at age 0 (the first spectrum at its own offset), M5 at M1's values is M1 and M5-first at M5's is M5, `+BEP` gives max(E₀, E₀ + ΔG) forward and max(E₀, E₀ − ΔG) backward in every reaction of N10–N15, M0-BEP keeps α = 0.5 |

**Documents**: `docs/reference-fitting_system.md` (§2.2 readings, §2.5a `+BEP`, §3 and §4 structures and `--structures`, new §11); `docs/reference-models_summary.md` (Peter's BEP in §1, M5 and `+BEP` in §2, ages and α in §3); `MODULES.md` (structures, extended stage).

**Not in git** (lab-derived): `notebooks/results/05_fit7/fit_config.json` (readings 0, 5 and 10 min; structures M1, M5, M5-first and their `+BEP` forms; everything else as in the sixth fit: network R1–R9 + N10–N15, EC/DEC 1:1 linear MD, 8 h heating steps of the 2 % tube, water of "TMSPa alone" declared at 65.6 mM, energy box ±3 σ; search: M1 4096 points and 32 starts, M5 and M5-first 8192 and 32 at 0 min and 2048 and 8 at 5 and 10 min; no earlier run as starts; every structure analysed in every reading: profiles, leave-one-out with 1024 points and 8 starts, ten sensitivities, solver tolerance at 0 min, recovery with the fit as truth at 0 min, traces, predictions; a `why` list) and `run_chain.sh` (phase 1: every stage with `--structures "M1,M1+BEP"`, to finish within 10 h; phase 2: every stage for all; `run_chain.sh M1` stops after phase 1; runs the code of the worktree).

**Checks**

- `test_seventh_fit_setup` passes; `test_sample_set_and_histories`, `test_free_ages_and_tables`, `test_structures_build_models` and every test of `test_medium` still pass. The rest of the suite (the synthetic fits of `test_estimation`, `test_fitting`, `test_reporting`) was not run.
- Static check of the setup (no fit): the three readings give every unknown age 0, 5 and 10 min (the glovebox tube keeps its known 8 h at 80 °C); M1 and M1+BEP 4 parameters, M5 and M5+BEP 10, M5-first and M5-first+BEP 12; parents and search per reading as set; `run_chain.sh` passes `bash -n`.
- One evaluation of the lab set at 0 min (no fit): M1 0.23 s and M1+BEP 0.16 s, 80 residuals, none penalised.

**Effects on what exists**: none for fits 1–6. Their folders list their own readings (or none: the default is still the four), and their structures, nesting and starts are unchanged (M5 is a new parent of M5-first only where a fit of M5 exists).


## 2026-10-23 - Eighth fit setup: 10 min and free ages up to 7 days, parallel chain for the cloud

Worktree `atom-to-reactor-fit8`, branch `fit/ten-min-free`, from the seventh fit's code (uncommitted, as in its worktree) plus:

| File | Change |
|---|---|
| `kinetics/fitting/residuals.py` | `set_free_age_bounds` / `get_free_age_bounds`: the range of the free ages for a whole run (default 1 h to 30 days, `DEFAULT_FREE_AGE_BOUNDS_H`); the functions read the active range |
| `kinetics/fitting/config.py` | Setup key `free_age_bounds_h` ([low, high] in hours; None: the default), applied by `apply_fit_config` |
| `kinetics/fitting/candidates.py`, `reporting.py`, `__init__.py` | The age parameter and the exact-age search read the active range (were the constant) |
| `scripts/run_fit.py` | Option `--readings`: one call limited to some readings of the setup (`ONLY_READINGS`), also the recovery cases. Night key `recovery` (default true): false leaves the recovery out of the night, so that it runs as its own stage beside it. `--quick` also shrinks the recovery search. Each fit records `free_age_bounds_h`. Docstring |

**Not in git** (lab-derived): `notebooks/results/05_fit8/fit_config.json` (readings 10 min and free ages 10 min to 7 days; the seventh fit's structures and setup; full search everywhere, started also from the fits of 05_fit7; recovery of every structure in both readings, profiles everywhere, leave-one-out only where a structure fits, no sensitivities, no indistinguishability, no tolerance check; a `why` list), `run_chain.sh` (every structure and reading its own process: extended per reading, recovery and night per structure and reading, then the stages again without limits for the tables, traces, predictions; `smoke` runs it with `--quick` in `05_fit8_smoke`) and `cloud.sh` (the VM: n4-highcpu-64 Spot, create, setup, smoke, start, status, resume, back, delete).

**Checks**

- The free-age, config, age and scenario tests pass (5). Static check with the eighth fit's setup (no fit): free ages from 10 min (log10 −0.78) to 7 days (2.23) on the 2 % tube, the unheated tubes profiled, the heated probe tube at 10 min; M5-first 13 parameters, one evaluation 0.1 s; 3 starts from the seventh fit per reading; `--readings free` limits the readings and the recovery cases. `run_chain.sh` and `cloud.sh` pass `bash -n`. The chain itself was not run (its smoke is the first step on the VM).

**Effects on what exists**: none for fits 1–7 (the range stays 1 h to 30 days without the key; without `--readings` and the night key every call runs as before).
