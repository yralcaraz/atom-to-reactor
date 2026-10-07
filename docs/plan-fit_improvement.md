# Fit improvement plan: the third fit

| | |
|---|---|
| **Status** | Proposal, not implemented |
| **Source** | Y. Alcaraz Galván |
| **Scope** | What to change in the fit after the first and second fits, why, where in the code, how to judge the result, and what the plan cannot settle. |
| **Applies to** | Branch `fit/first-fit`, commit `d2c5a0b` plus the uncommitted second-fit changes (2026-10-07) |
| **Evidence** | First fit (`notebooks/results/05/`), second fit (`notebooks/results/05_fit2/`, read on 2026-10-07 while its last stages were still running) and screening runs of 2026-10-07 (Appendix A; not stored in the repository) |
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
| C3 | Every fit must also hold with a trace of TMSOH at mixing (acceptance rule) | 0 |

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

| Reading | χ², parameters unchanged, 1 µM | Fitted again, 1 µM | Fitted again, 1 mM | ΔG‡ R1 (none → 1 µM) |
|---|---|---|---|---|
| 1 day | 2128 | χ² 32.7, fits | χ² 54.0, worst miss 3.77, fails | 1.32 → 1.67 |
| 7 days | 1449 (already at 1 nM) | χ² 34.0, fits | pending | 1.67 → 1.55 |
| Free ages | 30.4 | χ² 30.4, fits | pending | 1.09 → 1.09 |

In the first fit the 1 mM column was: 1 day fails (worst miss 3.86), 7 days fits with ΔG‡ R1 moved to 1.69 eV, free ages unchanged.

**Result**

- The equal-age fits depend on a start that is pure to better than 1 µM. No reagent is.
- The free-age fit does not react to the trace at any level.

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
- Unchanged with 1 µM TMSOH at mixing.
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
| [observables.py](../kinetics/data/observables.py) | `LAB_SAMPLES['TMSPa alone']`: role `check` → `fit`; a new key that declares which species of the recipe is unknown (H2O). Rebuild `lab_observables.json` with `scripts/build_observables.py`. |
| [residuals.py](../kinetics/fitting/residuals.py) | `_lab_sample` carries the declaration. `calculate_residuals`, `tabulate_residuals` and `_profile_age` take `added_M = {sample: {species: mol/L}}` and pass it to `_initial_composition`. Today this only exists at build time (`overrides`). |
| [candidates.py](../kinetics/fitting/candidates.py) | A parameter `water TMSPa alone` (kind `recipe`, box 0–0.3 M, DECLARED) added in `FitStructure.parameters(samples)` the way an age is. `build` returns `added_M`. |
| [estimation.py](../kinetics/fitting/estimation.py) | `FitProblem.residuals` and `.table` pass `built['added_M']`. `build_profile_grid` already gives a linear grid to any kind other than barrier or energy. |
| [run_fit.py](../scripts/run_fit.py) | `FITTED_SAMPLES` gains the tube (leave-one-out). `RECOVERY_TRUTHS` gains a case with free ages and a known water value. `check_tmspa_alone` stays only for fits made without the tube. Results go to a new folder (`notebooks/results/05_fit3/`). |
| [test_estimation.py](../tests/test_estimation.py), `synthetic_observables.py` | A synthetic tube with unknown water: the fit returns the water and the barriers used to make it. |

One difference from the present check: `check_tmspa_alone` starts from the composition of the first spectrum. In the fit the tube starts from the recipe plus the unknown water, and its age is found like that of any other tube. Both spectra are then predicted, not only the second.

### C3. Trace rule

| File | Change |
|---|---|
| [run_fit.py](../scripts/run_fit.py) | `stage_traces` already refits with 1 µM and 1 mM TMSOH. Add: (a) the trace goes to every tube whose recipe has TMSPA and no TMSOH, which now includes "TMSPa alone"; (b) a verdict table `trace_verdict.csv`; (c) `select_consistent` only accepts a fit that passes. |

A fit passes when, fitted again at 1 µM and at 1 mM:

1. it still fits (worst miss at most 3), and
2. no ΔG‡ that had a closed interval moves by more than 0.05 eV.

The two levels and the 0.05 eV are DECLARED.

### Order of work

1. C2 in the code, with its test.
2. Recovery on synthetic data, water unknown. Stop if it fails.
3. The ladder M0 → M3-split in the four readings on the new sample set.
4. Profiles, leave-one-out and sensitivities for M3-split.
5. Traces and the verdict (C3).
6. Predictions: tube B, the 13C of tube A after 10 days, the paper observations.

Runtime: about that of the second fit, roughly a day on 12 workers.

## 4. How the result is judged (fixed before running)

| Question | Success | If not |
|---|---|---|
| Does the pipeline recover a known water value and known barriers? | Truth inside the 95 % interval in every seed | Do not run on lab data; the water is not separable with this design |
| Does M3-split fit with the tube included, free ages? | Worst miss at most 3 | Report the failure; the tube goes back to `check` |
| Is the water determined? | Interval closed on both sides | Report as not determined; quote barriers with the water profiled |
| Is `g transfer` determined? | Interval closed on both sides | The plan did not achieve its aim; say so |
| Does the fit pass the trace rule? | Yes at both levels | Not accepted |
| Tube B (hold-out) | Worst miss at most 3 | Report; do not refit to repair it |
| What rests on one tube? | Leave-one-out table | Mark every barrier that moves by more than its interval |

## 5. Limitations

**Of the evidence**

- Steps 4 and 5 are screening: one local least-squares polish per start, `g transfer` held by hand, water scanned on a grid. No joint fit, no global search, no profile, no leave-one-out.
- The 1 mM trace refits of the second fit were not finished for 7 days and free ages when this was written.

**Of the plan**

- The ages stay unknown. Free ages are chosen because they survive the trace test. The data do not exclude equal ages.
- A third reading exists and is not tested: Gogoi's thesis describes the 0.5 % tube as stopped ("the reaction stops at BMSPa, while most of TMSPa remains unreacted") and the water series as an equilibrium moved by water. The thesis gives no waiting times.
- About 60 mM water (near 900 ppm) in a tube with no water added is high and cannot be checked. The model puts all of it in at mixing; water entering slowly over the 10 days would look different and is not modelled.
- `g transfer` would rest on one tube: two spectra of 16 scans, standard errors of 0.03–0.05 per share.
- Tube B gets slightly worse in the screening (χ² 3.8 → about 7 on 4 values).
- The only independent kinetic check is used up. After C2 the checks left are tube B (noisy), the 13C of tube A after 10 days and the paper statements.
- The trace rule depends on declared numbers (1 µM, 1 mM, 0.05 eV).

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
