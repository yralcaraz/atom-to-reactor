# Fit improvement plan: the third fit

| | |
|---|---|
| **Status** | Implemented in the code on 2026-10-07 (section 3). The third fit was launched the same day at 16:52; its results go to `notebooks/results/05_fit3/` and are not in this document |
| **Source** | Y. Alcaraz Galván |
| **Scope** | What to change in the fit after the first and second fits, why, where in the code, how to judge the result, and what the plan cannot settle. |
| **Applies to** | Branch `fit/first-fit`, commit `e3b051b` plus the changes of section 3 |
| **Evidence** | First fit (`notebooks/results/05/`), second fit (`notebooks/results/05_fit2/`, finished 2026-10-07 15:24) and screening runs of 2026-10-07 (Appendix A; not stored in the repository) |
| **Confidentiality** | Quotes numbers derived from the lab NMR spectra. Do not commit before it is classified. |

Words used here:

- **Age**: time from mixing a tube to its first spectrum. It was not recorded.
- **Reading**: an assumption about the ages. Equal ages of 1 h, 1 day or 7 days (`short`, `middle`, `long` in the code), or free ages (`free`): each tube has its own.
- **Miss**: (predicted − measured) / standard error of one measured value. A model fits when its worst miss is at most 3.
- **χ²**: the sum of the squared misses. The second fit has 77 values.

## 1. Summary

The model structure stays as it is (M3-split). Three things change around it:

| # | Change | New parameters |
|---|---|---|
| C1 | Free ages become the base reading; equal ages stay as a sensitivity run | 0 |
| C2 | The tube "TMSPa alone" enters the fit, with its water content as an unknown | 1 |
| C3 | Trace rule: the parameters of a fit, unchanged, must still fit with a trace of TMSOH at mixing | 0 |

Expected effect: the barrier of the transfer reactions (R5–R7), today only bounded from above, gets a value (about 0.7–0.8 eV in the screening), and the readings that depend on a perfectly pure start are rejected by rule.

Not changed: the reaction network, the barrier shape, the four free reaction energies and the error model.

## 2. Logic

### Step 1. What did the two fits settle?

M3-split, second fit:

| Reading | χ² | Worst miss | ΔG‡ R1 (eV) | ΔG‡ R5 (eV) |
|---|---|---|---|---|
| 1 h | 659.5 | 13.4 | 0.98 | 0.64 |
| 1 day | 32.7 | 1.93 | 1.32 | 0.59 |
| 7 days | 32.4 | 1.93 | 1.67 | 0.55 |
| Free ages | 30.4 | 1.93 | 1.09 | 0.90 |

**Result**

- Settled, the same in the three readings that fit: the shifts of the reaction energies (R1 +0.42 to +0.43, R2 +0.19 to +0.21, R3 +0.25 to +0.27, R4 −0.24 to −0.28 eV) and ΔG‡ at 80 °C of solvent attack (1.29 eV) and condensation (1.24 eV).
- Not settled: the barriers of R1 and of transfer. Three different pairs give the same χ².

--> The data cannot choose a reading. Something else has to.

### Step 2. Which reading survives a trace of TMSOH at mixing?

The 0.5 % H2O tube has 93 % TMSPA left; the 2 % tube has none. With one water molecule per reaction step, four times more water makes hydrolysis four times faster, which is not enough. At equal ages the fit explains the gap with a reaction that waits and then goes off: R1 almost closed, transfer very fast, started by the first TMSOH formed.

Second fit, M3-split, TMSOH added at mixing to the two water tubes:

| Reading | χ² without trace | χ², parameters unchanged, 1 µM | χ², parameters unchanged, 1 mM | Fitted again, 1 µM | Fitted again, 1 mM |
|---|---|---|---|---|---|
| 1 day | 32.7 | 2128 | 1449 | χ² 32.7, fits, ΔG‡ R1 1.32 → 1.67 | χ² 54.0, worst miss 3.77, fails |
| 7 days | 32.4 | 1449 (already at 1 nM) | 1450 | χ² 34.0, fits, ΔG‡ R1 1.67 → 1.55 | χ² 39.8, worst miss 2.50, fits |
| Free ages | 30.4 | 30.4 | 30.6 | χ² 30.4, same parameters | χ² 30.4, same parameters |

**Result**

- The equal-age parameter sets depend on a start that is pure to better than 1 µM. No reagent is.
- Fitted again with the trace, an equal-age reading can sometimes be rescued, but with other barriers each time.
- The free-age parameter set does not react to the trace at any level.

--> Free ages are the base reading (C1), and the trace test becomes a rule (C3).

### Step 3. What does the free-age fit get wrong?

The tube "TMSPa alone" (no water added, two spectra 243 h apart) was never fitted. Shares as TMSPA / BMSPA / MMSPA / H3PO4:

| | First spectrum | After 243 h |
|---|---|---|
| Measured | 0.81 / 0.12 / 0.05 / 0.04 | 0.23 / 0.69 / 0.04 / 0.03 |
| Free-age fit, best water (200 mM) | 0.84 / 0.14 / 0.02 / 0.00 | 0.31 / 0.38 / 0.25 / 0.06 |

**Result**

- χ² 39.7 on 8 values. The model lets BMSPA go on to MMSPA; the tube stops at BMSPA.
- It is the only tube that shows the phosphate ladder moving at room temperature. Every other tube gives one composition that does not change.

--> The tube that carries the rate information is outside the fit.

### Step 4. Does this need a new reaction?

No. A silyl group can already pass between phosphates through TMSOH and HMDSO (R5 forward, R6 backward). The fitted tubes only say that the transfer barrier `g transfer` is below about 1.06 eV. Free ages, M3-split, `g transfer` held and the rest fitted again:

| `g transfer` (eV) | χ², fitted tubes (77 values) | χ², TMSPa alone (8 values) | χ², tube B (4 values) |
|---|---|---|---|
| 1.01 (second fit) | 30.4 | 39.7 | 3.8 |
| 0.90 | 34.2 | 25.5 at 140 mM water, the end of the scan, still falling | 6.4 |
| 0.80 | 35.3 | 3.3 | 6.9 |
| 0.70 | 33.7 | 3.2 | 7.1 |

**Result**

- With `g transfer` at 0.7–0.8 eV the same model reproduces the tube, at a cost of 3 to 5 in χ² on the fitted tubes.
- Unchanged with up to 0.1 mM TMSOH at mixing (χ² 36.9 → 38.1 on the 85 values of the third-fit sample set).
- With 1 mM it fails: χ² 146, and the 0.5 % tube is predicted with 22 % BMSPA against 6 % measured. A fast transfer step multiplies a TMSOH impurity in a tube that holds water. Found on 2026-10-07 after the code was written; see section 5.
- Other ideas were tried and dropped (Appendix A): simpler energy corrections, a higher order in water, hydrolysis catalysed by P–OH, and a direct silyl exchange between phosphates. The last one works equally well but adds a reaction that the literature search did not confirm, and it is not needed.

--> Put the tube into the fit (C2). The structure stays M3-split.

### Step 5. Can the unknown water of that tube be fitted?

`g transfer` 0.70 eV, other parameters from Step 4, water at mixing scanned:

| Water at mixing | χ² | Predicted after 243 h |
|---|---|---|
| 40 mM | 29.7 | 0.49 / 0.49 / 0.02 / 0.00 |
| 50 mM | 10.0 | 0.37 / 0.60 / 0.03 / 0.00 |
| 60 mM | 3.2 | 0.26 / 0.67 / 0.06 / 0.00 |
| 70 mM | 5.1 | 0.18 / 0.72 / 0.10 / 0.00 |
| 80 mM | 11.8 | 0.11 / 0.72 / 0.16 / 0.01 |
| 100 mM | 36.1 | 0.05 / 0.62 / 0.30 / 0.03 |

**Result**

- Water is bounded on both sides: too little leaves TMSPA, too much pushes on to MMSPA. About 55–72 mM, near 900–1100 ppm.
- The barriers are shared with all tubes; the water belongs to this tube alone. That is why the two can be separated.

--> One extra parameter is enough.

## 3. Implementation

### C1. Free ages as the base reading

| File | Change |
|---|---|
| [run_fit.py](../scripts/run_fit.py) | `select_consistent`: look at the free-age fits first; report the equal-age fits as sensitivity. The rule "common age takes precedence" was fixed before the first fit and is replaced on purpose, for the reason in Step 2. |
| Notebook 05 | Lead with free ages; equal ages move to the sensitivity section. |

Nothing changes in the residuals: with free ages, a tube that was never heated already gets its age found inside each evaluation (`age_mode = 'profile'`). That applies to "TMSPa alone" too.

### C2. "TMSPa alone" in the fit, water unknown

| File | Change |
|---|---|
| [observables.py](../kinetics/data/observables.py) | `LAB_SAMPLES['TMSPa alone']`: role `check` → `fit`, and `'unknown_M': {'H2O': (0.0, 0.3)}`, the range in M a fit may give the water (DECLARED). `lab_observables.json` rebuilt from the raw spectra: the shares are identical, only the role, the new key and one flag differ. The second-fit file is kept beside it as `lab_observables_second_fit_261006.json`. |
| [residuals.py](../kinetics/fitting/residuals.py) | A sample carries `unknown_M`. `calculate_residuals`, `tabulate_residuals`, `find_free_age` and `calculate_predicted_shares` take `added_M` (amounts at mixing on top of the recipe) and pass it to `_initial_composition`, also when the age is found inside the evaluation (`_profile_age`). |
| [candidates.py](../kinetics/fitting/candidates.py) | A parameter `c0 H2O TMSPa alone` (kind `recipe`, unit M) added in `FitStructure.parameters(samples)` the way an age is. `build` returns `added_M`. `embed_parent_theta` passes it on. |
| [estimation.py](../kinetics/fitting/estimation.py) | `FitProblem` and `simulate_synthetic_shares` pass `added_M`. The profile grid of a `recipe` parameter is a set of multiples of its best value (`RECIPE_FACTORS`). |
| [run_fit.py](../scripts/run_fit.py) | `FITTED_SAMPLES` gains the tube (leave-one-out). `RECOVERY_TRUTHS`: cases A and B get a water value, and a new case C has free ages, a transfer barrier of 0.75 eV and the water profiled. `complete_start` gives a start without the parameter 50 mM. The `registered` stage leaves the tube out, because a registered model has no value for its water. |
| [test_estimation.py](../tests/test_estimation.py) | The sample set and parameter lists with the tube; `added_M` changes only the residuals of its own sample; a new test recovers 60 mM from synthetic shares inside a closed interval. |

To reproduce the second fit with this code, build the sample set with `exclude=('TMSPa alone',)`: χ² 32.706 (1 day) and 30.369 (free ages) come back exactly.

`check_tmspa_alone` is kept in the `predictions` stage as a second reading. It starts from the composition of the first spectrum. In the fit the tube starts from the recipe plus the unknown water, and its age is found like that of any other tube, so both spectra are predicted.

### C3. Trace rule

| File | Change |
|---|---|
| [run_fit.py](../scripts/run_fit.py) | `tabulate_trace_verdicts` evaluates every fit that fits, with its parameters unchanged, at 1 µM and 1 mM TMSOH in every tube whose recipe has TMSPA and no TMSOH (`TRACE_SAMPLES`, which now includes "TMSPa alone"). It writes `trace_verdict.csv`. `select_consistent` skips a fit that is rejected. `stage_traces` still scans and refits the main structure, for information. |

A fit passes when its own parameter set, unchanged, still fits (worst miss at most 3) at both levels. The two levels are DECLARED.

The first draft of this plan judged the fit made again with the trace, and asked that no determined barrier move by more than 0.05 eV. That was dropped: the second fit showed the threshold deciding by itself (in the 7-day reading the barrier of R2 and R3 moved by 0.051 eV), and the unchanged parameters already separate the readings by a factor of 50 in χ² (Step 2). The simpler rule needs no refit, so it is applied to every structure.

### Order of work

1. C1–C3 in the code, with tests. Done.
2. Recovery on synthetic data, water unknown (`recovery`). Stop if it fails.
3. The ladder M0 → M3-split in the four readings on the new sample set (`baselines`, `extended`).
4. Profiles, leave-one-out and sensitivities for M3-split (`night`).
5. Trace scan and refits (`traces`).
6. Predictions: tube B, the 13C of tube A after 10 days, the paper observations (`predictions`).

Steps 2 to 6 are the stages of `scripts/run_fit.py`, run with `--results notebooks/results/05_fit3`. Runtime: about that of the second fit, roughly a day on 12 workers.

## 4. How the result is judged (fixed before running)

| Question | Success | If not |
|---|---|---|
| Does the pipeline recover a known water value and known barriers? | Truth inside the 95 % interval in every seed | Do not run on lab data; the water is not separable with this design |
| Does M3-split fit with the tube included, free ages? | Worst miss at most 3 | Report the failure; the tube goes back to `check` |
| Is the water determined? | Interval closed on both sides | Report as not determined; quote barriers with the water profiled |
| Is `g transfer` determined? | Interval closed on both sides | The plan did not achieve its aim; say so |
| Does the fit pass the trace rule? | Its unchanged parameters fit at both levels | Not accepted. `trace_verdict.csv` and `trace_sensitivity.csv` keep χ² at every level, so the levels can be discussed without a new run |
| Tube B (hold-out) | Worst miss at most 3 | Report; do not refit to repair it |
| What rests on one tube? | Leave-one-out table | Mark every barrier that moves by more than its interval |

## 5. Limitations

**Of the evidence**

- Steps 4 and 5 are screening: one local least-squares polish per start, `g transfer` held by hand, water scanned on a grid. No joint fit, no global search, no profile, no leave-one-out.

**Of the plan**

- The two aims may pull against each other. The transfer barrier that reproduces "TMSPa alone" in the screening (0.70 eV) passes the trace rule at 1 µM and fails it at 1 mM. Whether the joint fit finds a parameter set that does both is open. If none does, the choice is between a lower upper level for the rule (the screening set holds up to 0.1 mM) and reading the 0.5 % tube as a bound on the impurity of the TMSPA stock.

- The ages stay unknown. Free ages are chosen because they survive the trace test. The data do not exclude equal ages.
- A third reading exists and is not tested: Gogoi's thesis describes the 0.5 % tube as stopped ("the reaction stops at BMSPa, while most of TMSPa remains unreacted") and the water series as an equilibrium moved by water. The thesis gives no waiting times.
- About 60 mM water (near 900 ppm) in a tube with no water added is high and cannot be checked. The model puts all of it in at mixing; water entering slowly over the 10 days would look different and is not modelled.
- `g transfer` would rest on one tube: two spectra of 16 scans, standard errors of 0.03–0.05 per share.
- Tube B gets slightly worse in the screening (χ² 3.8 → about 7 on 4 values).
- The only independent kinetic check is used up. After C2 the checks left are tube B (noisy), the 13C of tube A after 10 days and the paper statements.
- The trace rule depends on two declared levels (1 µM, 1 mM).
- The water range of the tube (0 to 0.3 M) is declared.

**Of the model, left as they are**

- One water molecule per hydrolysis step is assumed. Orders 2 and 3 were screened at equal ages only.
- The four fitted reaction energies sit 0.2–0.4 eV from the computed ones. The solvation term is an MD energy, not a free energy. R8 and R9 keep their computed energies.
- ΔS‡ = 0. Hydrolysis and transfer are seen only at room temperature, condensation and solvent attack only at 80 °C.
- The 0.5 % tube has open audit flags (lines with the pattern of PF6⁻, D2O lock).
- With free ages the 2 % tube must be at least about 9 h old at its first spectrum.

## Appendix A. Screening runs of 2026-10-07

Second-fit sample set. Local polish from hand-made starts (up to 30 least-squares iterations), then χ² at the report tolerance. χ² over 73 values, or 77 where the four energies are free.

| Idea | Reading | χ² | Worst miss | Verdict |
|---|---|---|---|---|
| One scale factor on all solvation energies, R1 split | Free | 105.3 | 7.9 | Fails |
| Same, plus a shift of H2O | Free | 76.9 | 4.3 | Fails |
| Same two | 1 day | 1305; 1238 | 23.4; 19.3 | Fails |
| Two species shifts (H2O −0.22, TMSPA −0.23 eV), R1 split | 1 day; free | 91.9; 85.7 | 7.1; 7.0 | Fails, on MMSPA of the 2 % tube |
| Hydrolysis of order 2 in water, one hydrolysis barrier | 1 day; 7 days | 215.9; 174.1 | 8.5; 8.0 | Fails |
| Hydrolysis of order 3 in water, one hydrolysis barrier | 1 day; 7 days | 121.6; 95.6 | 6.3; 5.2 | Fails |
| Hydrolysis catalysed by P–OH, one hydrolysis barrier | 1 day | 40.8 | 2.0 | Fits, but χ² 1562 with 1 µM TMSOH |
| Same | 7 days | 395 | 11.4 | Not found |
| Direct silyl exchange between phosphates added to M3-split | Free | 32.0 | — | Fits; TMSPa alone χ² 3.5; unchanged with 1 µM TMSOH |
| M3-split, `g transfer` held at 0.70 or 0.80 eV | Free | 33.7; 35.3 | — | Fits; TMSPa alone χ² 3.2; 3.3; unchanged with 1 µM TMSOH |

Literature (open web only): no report of direct silyl exchange between trimethylsilyl phosphates was found. Read: N. Gogoi, doctoral thesis, Uppsala 2024 (urn:nbn:se:uu:diva-523121), section 4.2.2.2; Gogoi et al., J. Phys. Chem. C 2024, 128, 1654.

The screening scripts are not in the repository.
