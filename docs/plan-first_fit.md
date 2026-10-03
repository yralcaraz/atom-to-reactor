# Plan: first fit of the kinetic models to the lab NMR data

I read the code, both docs, the five notebooks and the result CSVs, ran the six test scripts (all pass), timed the reactor on this laptop, and made a few forward runs with existing functions. Nothing in the repository was changed.

Evidence tags:
- **[C]** read in code, a stored notebook output or a result CSV.
- **[R]** run today with existing functions at a prescribed point; nothing was optimised.
- **[I]** my inference or a hand calculation.

## Five findings that change the starting point

1. **M3 alone cannot work.** Freeing ΔG_rxn of R2 and R3 leaves the conflict between the 0.5 % and 2 % water samples at 0.041 eV, however the energies are shifted. With one hydrolysis barrier, the freed model reaches the measured 2 % composition only after 34–130 years. The smallest structure that can work also needs separate hydrolysis barriers: 9 parameters, not 6. [R]
2. **A common age never reconciles hydrolysis.** The two water samples disagree by a factor of 115, 11 and 7 in rate for common ages of 1 h, 1 d and 7 d. They agree only if the 2 % sample was about a week older than the 0.5 % one; then ΔG‡(R1) is confined to 1.07–1.08 eV. Your three scenarios need a fourth reading with free ages. [R]
3. **The only lower bound on hydrolysis comes from a suspect sample.** L1 rests on the 0.5 % sample, the one sample locked on D₂O and showing three ³¹P lines 4.4 ppm apart near −144 ppm, the pattern of PF₆⁻. If it holds LiPF₆, it is outside the network. [C, notebook 04 §5 and Appendix A]
4. **Few numbers carry information.** The fit set holds 32 share values, but only 7 say more than "zero" or "unchanged". [C]
5. **The cost fits one night.** One evaluation of the whole data set takes 0.16 s. The laptop gives 14 evaluations per second on all cores, a gain of 4, not 12. The full programme needs about 300 000 evaluations, roughly six hours. [R]

## 0. Questions

**What your answers fixed.**
- The fit is diagnostic: a best fit per structure with the misfit per observation, plus the blind spots.
- No new data and no answers are assumed. Unknown times are scenarios, never priors.
- No sampling stage; one night of computing.
- Parameters come first, storage at 25 °C is the main prediction, Peter's protocol tests the temperature blind spot, and gas is outlook.

**Still open for you** (I proceeded on the defaults):

| # | Question | What it changes | Default used |
|---|---|---|---|
| Y1 | On the Mac: do other TMSPa samples show the −144 ppm lines? What does the ¹⁹F spectrum of the glovebox TMSOH sample show? | Whether the 0.5 % sample, and the sample behind the solvent-attack barrier, are inside the model | Keep both, flag them, and read the fit without each from leave-one-sample-out |
| Y2 | Scenario ages | Every scenario-dependent result | 1 h, 1 d, 7 d for all samples, plus free ages between 1 h and 7 d |
| Y3 | Error floor on shares | How many SE a misfit is worth; notebook 03's "7 SE" uses within-spectrum errors only | ³¹P: larger of 0.01 and 5 % of the share. ¹³C (NOE on, 2 s delay): larger of 0.03 and 15 % |
| Y4 | Drop E4 in favour of tube A? | Double counting (see §1) | Drop it from the fit; keep it as a check |
| Y5 | "TMSPa alone": fit with a free water content, or check only? | One more parameter, and the transfer bound from tube A | Check only, as in notebook 03 |
| Y6 | Is "part of the water unavailable" (M4) a rate-law variant for you? | Whether the "equilibrium or exhausted water" blind spot gets numbers | It only scales [H₂O]₀, so I include it as one optional parameter |
| Y7 | New dependencies | Optimiser and parallelism | None: scipy and the standard library |
| Y8 | Rule for "the first consistent parameter set" | What gets registered | Smallest structure whose largest standardised residual is ≤ 3; a common-age scenario takes precedence over free ages |
| Y9 | Names | The thesis table | Keep M0, M1, M1-split, M3, M4; add M3-split; reserve M2 for acid catalysis, as in the guide |

**Only Neha can answer.** None is needed; each would collapse a scenario.

| Question | One answer | The other | Covered by |
|---|---|---|---|
| When was each sample mixed, above all the 2 % relative to the 0.5 % sample? | 2 % about a week older: first-order hydrolysis stands, ΔG‡(R1) = 1.07–1.08 eV | Similar ages: every structure that is first order in water fails by a factor ≥ 7 | The four age readings |
| Is there LiPF₆ in the 0.5 % sample? | No: L1 stands | Yes: hydrolysis has no lower bound | Leave-one-sample-out |
| How long was each step of the 2 % ramp held? | Minutes (files): heating says little | 8 h (paper): a strong test of any equilibrium reading | One sensitivity run |
| Was the glovebox sample held 8 h? | Yes: ΔG‡(R8) is two-sided | Otherwise the interval shifts 0.02 eV per factor 2 | One sensitivity run |
| Water content of the solvents | Known: separates M3 from M4 | Unknown: they stay indistinguishable | Profile of the available-water fraction |
| Is the ²⁹Si line at 19.08 ppm the second ring-opening product? | — | — | R9 reported as not determined |

## 1. Where we stand

**Confirmed [C].**
- The registry holds the four models with the values you list.
- `kinetics/fitting/` has no fit: the only optimiser calls are bisection and one bounded scalar search for window edges.
- Notebook 03's verdicts match its CSVs.
- The lab set is 9 samples and 89 spectra.

**One correction to your summary.** Uncertainties on area shares already exist: `calculate_area_shares` gives noise and baseline spread. What is missing is the systematic part and the curated file. [C]

**Contradictions and corrections.**

1. **ΔG_rxn(T) is not active.** In `wb97mv` mode the stored 298 K value is returned at every temperature ([gas.py:49-59](kinetics/thermo/gas.py#L49-L59)); `qRRHO` falls back to it for lack of frequencies. K_eq changes only through 1/RT ([reaction.py:48](kinetics/thermo/reaction.py#L48)). Notebook 01 Block 4 estimates the neglected change at ≤ 0.05 eV between 20 and 80 °C. [C]
2. **"M3 admissible at 0.75 SE" is a two-dimensional marginal.** [R]
   - The most probable species shifts behind it move TMSOH by +0.09 eV and water by −0.04 eV.
   - That flips condensation from +0.10 to −0.12 eV and moves R1 (−0.47 to −0.36) and R8 (−0.05 to −0.14).
   - Restricted to the phosphates, the same shift costs 2.1 SE; restricted to MMSPA and H₃PO₄, 3.2 SE.
3. **M3 inherits M1's hydrolysis conflict**, as in finding 1. By notebook 03's own rule it is rejected under the 24 h age exactly as M1 is. [R]
4. **The rejections of M1 and M1-split rest on few points.** M1 was tested at three values of g_H, M1-split at two points. The argument behind them is sound, but "at any parameter value" still needs the optimiser. [C]
5. **E4 and L3 are probably the same experiment counted twice.** Tube A shows 0.895 ± 0.043 conversion, below E4's assumed ≥ 0.95. [C]
6. **The model makes the wrong ring-opening product.** [R]
   - In the glovebox sample the model puts 0.42 of the Si in TMSOdiEG and 0.001 in TMSOEG, because k(R9)/k(R8) = 709 at 80 °C.
   - The spectra are assigned to TMS-EG.
   - The lumped observable hides this ([validation.py:36-37](kinetics/reactor/validation.py#L36-L37)), and it doubles the predicted CO₂.
7. **Conditions differ between parts of the repo.** Room temperature is 20 °C in the protocol, 21.7–22.9 °C in the lab files and 25 °C in the Gogoi file. EC is 7.1 M for lab samples and 14 M in the protocol. `T_ref` is 298.15 K, where no data sit. [C]
8. **The existing forward functions do not fit the heated samples.** [C]
   - `evaluate_observation` is isothermal ([constraints.py:151-173](kinetics/fitting/constraints.py#L151-L173)).
   - `simulate_design` is hardwired to Radau ([design.py:136](kinetics/fitting/design.py#L136)).
9. **Housekeeping.** [C]
   - `03_feasible_region.ipynb` has no stored outputs, and the raw spectra are not on this machine, so it cannot be re-run here.
   - The README still says `reference`, lists three notebooks and calls the lab data pending.
   - Notebooks 01 and 03 say the result CSVs are git-ignored, but they are tracked.
   - "M1–M4" means rate laws in the methodology doc and barrier structures in notebook 03.

## 2. The nine questions

### Q1. Data

| Sample | Spectra | Temperature | Time from files | Shares (± SE) | Carries |
|---|---|---|---|---|---|
| 0.5 % H₂O | 1 ³¹P | 22.3 °C | 18.9 h after its first spectrum | TMSPA 0.930 ± 0.019, BMSPA 0.060 ± 0.009; MMSPA and H₃PO₄ zero within 0.03 | 1 value, 2 zeros |
| 2 % H₂O | 6 ³¹P | 22.3–22.9 °C, with steps at 30–80 °C (≥ 17 min each) in between | 0.08 to 145 h | BMSPA 0.039, MMSPA 0.315, H₃PO₄ 0.662 (means); TMSPA 0.00 ± 0.01 | 2 values, 1 zero, 15 "unchanged" |
| TMSPa + TMSOH (A) | 1 ³¹P | 21.7 °C | 0.06 h | TMSPA 0.105 ± 0.043, BMSPA 0.687 ± 0.056, MMSPA 0.194 ± 0.042 | 2 values, 1 zero |
| TMSOH, probe | 3 ¹³C | 23, 80, 22 °C | at 80 °C for 8.7 min to 5.6 d | TMSOEG ≤ 0.03 ± 0.02, HMDSO 0.02 to 0.06 ± 0.03 | 6 zeros |
| TMSOH, glovebox | 1 ¹³C | after 80 °C | 8 h hold taken from the paper | TMSOH 0.418 ± 0.025, HMDSO 0.161 ± 0.056, TMSOEG 0.427 ± 0.042 | 2 values |
| E1 (paper) | — | 25 °C | 1 week | ring-opened ≤ 0.02 (our reading) | 1 one-sided |
| Held out: tube B | 1 ³¹P | 21.9 °C | 2.0 h | BMSPA 0.91 ± 0.18, MMSPA 0.01 ± 0.08 | test only |
| Check: TMSPa alone | 2 ³¹P | 21.7 °C | 242.8 h apart | TMSPA 0.81 to 0.23, BMSPA 0.12 to 0.69 (± 0.05) | water unknown |

- **Temperatures.** Everything quantitative sits at 22 °C or 80 °C. The ³¹P spectra recorded at 30–80 °C have the NOE on and cannot be read as amounts. [C]
- **Count.** There are 7 informative values against 1 to 9 structure parameters and six unknown times: three ages, two heating durations and one hold.
- **What must be built first.**
  - The observables file.
  - The time bases, as scenario readings.
  - The error model: stored SE plus a floor, with the six 2 % spectra treated as one composition and five "no drift" measurements.
- **Referenced axes** are not needed for shares, because the windows are set around the observed peaks.

What I expect each quantity to come out as:

| Quantity | Expected status |
|---|---|
| ΔG‡(R8) at 80 °C | Interval near 1.28–1.30 eV in every scenario [C, from L5] |
| ΔG‡(R4) at 80 °C | Lower bound; confounded with ΔG_rxn(R4) [I] |
| ΔG‡(R5) at 22 °C | Upper bound that moves with the age: g_T ≤ 1.06, 1.14, 1.19 eV [R] |
| ΔG‡(R1) at 22 °C | Inconsistent at common ages; 1.07–1.08 eV with free ages [R] |
| ΔG‡(R2), ΔG‡(R3) | Upper bounds only [I] |
| ΔG(R3) − ΔG(R2) | About +0.03 eV against the computed −0.02, if the plateau is an equilibrium [I] |
| Level of ΔG(R2) against available water | Not separable [I] |
| ΔS‡ | Not determined; one-sided for solvent attack, roughly ≥ −90 J mol⁻¹ K⁻¹ [I] |
| R7, R9 | Not determined |

### Q2. Model by model

| Structure | Free | Gets right | Cannot get right |
|---|---|---|---|
| **M0, capped BEP** (`peter_reference`) | 1 | Any single family: solvent attack at 1.28–1.30 eV, or TMSPA survival at 0.5 % water for E₀ ≥ 1.05 | TMSPA + TMSOH at any E₀ in 0.60–1.70 eV, because EC (7.1 M) takes the TMSOH first. Hydrolysis (≤ 1.01) and solvent attack (≥ 1.28) together: 0.27 eV apart. [C] |
| **M0, Marcus** | 1 | TMSPA + TMSOH for g = 0.67–1.14 | That together with solvent attack (≥ 1.31): 0.17 eV apart. [C] |
| **M1** (`level1` structure) | 4 | Transfer, condensation and solvent attack jointly | The two water samples at a common age (finding 2) [R]. The 2 % composition: BMSPA 0.15–0.20 and MMSPA 0.22 where 0.04 and 0.31 are measured, 7 SE [C]. A stationary composition: the run ends at 100 % H₃PO₄ [R]. |
| **M1-split** | 5 | As M1, with R1 slower than R2 and R3 | The 2 % composition while transfer stays within its bound: 6.1 SE at one tested point. [C] |
| **M3** | 4 + ladder energies | The 2 % composition: 1.3–1.5 SE at unfitted points, against 7.2–7.5 without the shift [R] | The water-sample conflict, and reaching the composition in time (finding 1). [R] |
| **M3-split** | 5 + ladder energies | Composition reached in about a day and held for weeks at one prescribed point [R] | Open: tube A's full composition is 2.4–5.3 SE off at six prescribed points [R]. The water conflict remains unless ages differ. |
| **M4** (part of the water unavailable) | +1 | A plateau by exhaustion, if 63–81 % of the added water is missing [I] | The ratio BMSPA·H₃PO₄/MMSPA²: 2.1 from computed energies, 0.26 measured [I]. |

- **Registered values.** `peter_reference` at 1.15 eV opens 62 % of the TMSOH in a week at room temperature, where none is seen. `level1` at its 0.80 eV placeholders consumes TMSPA in under a second. [C]
- **Barrier shape.** Marcus and the other smooth shapes differ by under 0.007 eV here; the data cannot tell them apart. [C]

### Q3. Peter's model

Yes, fit it as the baseline.
- Report the best single E₀ in each scenario, in both forms, with the per-observation misfit.
- Comparing BEP with Marcus at one parameter isolates the functional form. Comparing M0 with M1 isolates the number of barriers.
- Present it fairly: say what 1.15 eV was chosen for, show the best E₀ and not only 1.15, and give misfits as factors in rate.

### Q4. Level 1 or a new model

Do not overwrite `level1`. Fit its structure as the second baseline and register the survivor as a new model.

Order, simplest first: M0 (BEP, Marcus), M1, M1-split, M3, M3-split. The available-water fraction is added as one profile in the last structure.

Two design points:
- **Free the thermodynamics as reaction energies of R1–R4.** Use their computed covariance as the constraint and apply them as species shifts, so the Wegscheider cycles stay closed. Freeing "R2 and R3" alone hides what else moves (§1, item 2).
- **Tie R2 and R3 to one barrier.** The data give only upper bounds on both.

Solvent attack stays one family, with R9 declared not determined.

### Q5. Temperature

| | Where | Status |
|---|---|---|
| Eyring prefactor and Boltzmann factor | [rates.py:74-75](kinetics/microkinetics/rates.py#L74-L75) | Active |
| Rate constants recomputed per stage | [protocol.py:108-111](kinetics/reactor/protocol.py#L108-L111) | Active |
| ΔG_rxn(T) | [gas.py:49-59](kinetics/thermo/gas.py#L49-L59) | Dormant |
| ΔS‡ per family | [rates.py:58](kinetics/microkinetics/rates.py#L58), all zero in [parameters.py](kinetics/microkinetics/parameters.py) | Dormant; plumbing tested |
| Diffusion ceiling, temperature ramps | — | Dormant |
| K_eq(T) uses ΔG where van 't Hoff needs ΔH | theory doc, Block 6.4 | Inconsistent, ≤ 0.05 eV |
| `T_ref`, room temperatures, EC concentration | §1, item 7 | Inconsistent |

**Recommendation.**
- Fix ΔS‡ = 0 in every fit. Fitting ΔH‡ and ΔS‡ separately is the same problem in other coordinates.
- Set `T_ref` per family to where it is observed: 22 °C for hydrolysis and transfer, 80 °C for condensation and solvent attack.
- Report ΔG‡ at that temperature as the primary result.
- For Peter's protocol, show a band from a declared range of ΔS‡. With −150 to +50 J mol⁻¹ K⁻¹, k at 80 °C spans a factor of about 50. [I]
- What would narrow it: one quantitative spectrum of a partly reacted sample at a second temperature at least 40 K away, at a known time.

### Q6. Bounds or fit

Combine them, with one objective.
- Measured shares enter as Gaussian terms; a share measured as zero then acts as a one-sided constraint without special treatment.
- The paper statement E1 enters as a one-sided term.
- Reaction energies enter with their computed covariance.
- Unknown times enter as scenarios.

Intervals come from profile likelihood. The data are quantitative enough for this in composition, but not in time. Most barrier profiles will therefore be flat on one side, and the honest output there is a bound or "not determined".

Notebook 03's intervals stay as a regression test: with ± 2 SE windows the new objective must reproduce them.

A sampling stage adds nothing here. It would need priors on the ages, which you ruled out, and the parameter count is low enough for profiles.

### Q7. Methods

- **Fit rule, fixed before fitting.** A structure fits in a scenario if no standardised residual exceeds 3.
- **Optimiser.** A Sobol screen (256 to 2048 points), then bounded least squares from the 16 best points. The three TMSOH-only samples fix condensation and solvent attack first.
- **Uncertainty.** One-dimensional profiles with 13 points per parameter. Edges are refined in the middle and free-age scenarios.
- **Identifiability.** Singular values of the weighted sensitivity matrix at the optimum give the determined and undetermined combinations. Flat profiles confirm them.
- **Which observation drives what.** Leave-one-sample-out refits of the six fitted samples.
- **Structures that cannot be told apart.** Fit one structure to synthetic data generated from another, with the real times and errors. If it fits within noise, the data cannot separate them. Planned pairs:
  - equilibrium against exhausted water;
  - Marcus against Agmon–Levine;
  - M3 against M3-split.
- **Comparison.** The residual table comes first; a likelihood-ratio test is used for nested pairs; AICc is indicative only.
- **Validation.**
  - Tube B, which is weak: its prediction at the points I ran is within about 2 SE whatever the model. [R]
  - The paper's 1 % and 5 % water samples, which no fit uses.
  - TMSPa alone, by solving for the water it needs.

### Q8. Cost

Measured on the i7-1270P [R]:

| Item | Cost |
|---|---|
| One isothermal sample run (BDF, default tolerances) | 5–30 ms where the bounds place the barriers; 55–270 ms at the `level1` placeholders |
| Same with Radau | About twice; 2.5 s worst case |
| Whole fit set (7 samples, 21 segments), BDF, rtol 1e-6 | Median 163 ms; error below 3e-7 M |
| Same at rtol 1e-8 | 291 ms |
| Same through the existing Radau path | 206 ms |
| Same with LSODA | 57 ms, with outliers of 1.7 s |
| Throughput | 3.8 evaluations/s on one process; 12.1 on 4; 13.9 on 12; 15.2 on 16 |

- **The 2 % sample with its ramp is half the cost.**
- **The parallel limit is the laptop, not the code.** A pure-Python loop scales the same way (×4.9 at 16 processes).
- **The full 0.60–1.70 eV box is unsafe.** Two of 128 points did not finish in 20 s, both with g_SA ≤ 0.74. With Radau, one evaluation ran over 9 minutes before I stopped it.

Plan with 14 evaluations/s, about 50 000 per hour:

| Job | Evaluations | Time |
|---|---|---|
| Registered models as they are, four scenarios | 16 | seconds |
| M0 (both forms) and M1, four scenarios | 10 000 | 12 min |
| Notebook 03 consistency check | 1 300 | 2 min |
| M1-split, M3, M3-split, four scenarios | 65 000 | 80 min |
| Synthetic recovery (before real data) | 40 000 | 50 min |
| **Night:** profiles | 124 000 | 2.5 h |
| **Night:** leave-one-sample-out, four scenarios | 72 000 | 1.4 h |
| **Night:** indistinguishability tests, sensitivities, prediction bands | 101 000 | 2 h |

The night totals about 300 000 evaluations, roughly 6 h, which leaves 2 h in reserve.

- **HPC is not needed.** The Mac is needed only for the observables file.
- **Levers, by payoff:**
  - BDF in the history integrator.
  - rtol 1e-6 for the search and 1e-8 for reported numbers.
  - A search box that excludes the corners no observation allows.
  - A call budget per solve.
- **Not relied on.** LSODA for the screen and a vectorised Jacobian could each give a factor of 2–3, but neither is measured in the fit.

### Q9. Implementation

- **[kinetics/data/observables.py](kinetics/data/observables.py)** (new): `build_lab_observables`, `load_lab_observables`. It writes `data/lab_observables.json` (confidential), with domain flags per sample.
- **[kinetics/reactor/protocol.py](kinetics/reactor/protocol.py)**: `simulate_history(c0_M, segments, sample_times_s, *, model, method, rtol, atol)`.
- **[engine.py](kinetics/reactor/engine.py)**: a call budget in `integrate`.
- **[kinetics/thermo/uncertainty.py](kinetics/thermo/uncertainty.py)**: `calculate_species_shifts` and `build_shifted_species_database`.
- **[kinetics/microkinetics/models.py](kinetics/microkinetics/models.py)**: `ModelSpec` gains optional `network` and `species_shifts_eV`, plus a `fitted` status. A fitted model is then self-contained.
- **`kinetics/fitting/candidates.py`** (new): `FitStructure` and the registry of the six structures.
- **`kinetics/fitting/residuals.py`** (new): `build_sample_set` (applies a scenario), `calculate_residuals`, `tabulate_residuals`.
- **`kinetics/fitting/estimation.py`** (new): `fit_structure`, `calculate_profile`, `find_confidence_interval`, `calculate_parameter_directions`, `fit_leave_one_out`, `compare_structures`, `simulate_synthetic_shares`, `calculate_prediction_band`.
- **`scripts/run_fit.py`**: staged, resumable runner that writes to `notebooks/results/05/`.
- **`demo/`**: `plot_residual_map`, `plot_fit_across_scenarios`, `plot_profiles`, `plot_influence`, `format_structure_table`.
- **Notebook**: `05_first_fit.ipynb`, which only loads and shows results.

**Tests** (synthetic fixture, no private data):
- `simulate_history` equals the batch and design reactors.
- Species shifts keep the cycles closed and reproduce the 0.75 SE.
- The window objective reproduces `m1_intervals.csv` within 0.005 eV.
- Six identical spectra give zero drift residuals.
- Nesting never raises χ².
- Synthetic recovery puts the truth inside the profile interval and returns "not determined" for a parameter the design cannot see.
- rtol 1e-6 and 1e-8 agree within 0.1 SE.

## 3. Phases

Each phase ends in a result you can use.

| # | Goal and what gets built | Check | Computing | Result |
|---|---|---|---|---|
| 1 | Observables file, domain audit, error model. On the Mac. | File reproduces the two `lab_shares` CSVs; scatter of the six 2 % spectra matches their noise | seconds | The Q1 table with final numbers |
| 2 | `simulate_history`, scenarios, residuals; no optimiser. Depends on 1. | Regression tests; notebook 03 intervals reproduced | seconds | Misfit of the four registered models as they are |
| 3 | Structures, optimiser, profiles, runner; measure throughput and rescale. Depends on 2. | Synthetic recovery | 50 min | What data like ours can determine, before real data are touched |
| 4 | Baselines M0 and M1 in four scenarios. Depends on 3. | Best E₀ lies in notebook 03's intervals | 12 min | The baseline table for Peter |
| **5** | **M1-split, M3, M3-split.** Depends on 4. | **Nesting; fit rule applied** | **80 min** | **The fit becomes discussable here: every structure has a best fit and a per-observation misfit in four scenarios** |
| 6 | Blind spots: profiles, directions, leave-one-out, indistinguishability, sensitivities. Depends on 5. | Intervals stable from rtol 1e-6 to 1e-8 | one night | Interval or "not determined" per parameter; what drives each value |
| 7 | Storage at 25 °C with its band; Peter's protocol with the ΔS‡ band. Depends on 6. | Band contains every accepted parameter set | under 1 h | The two predictions |
| 8 | Notebook 05, registered model, README and classification registry. Depends on 7. | Notebook runs from cached results | — | The thesis material |

Phase 2 stops the work if notebook 03 is not reproduced. Phase 3 stops it if recovery fails.

## 4. Decisions that are yours

- Y1 to Y9 above.
- The declared range of ΔS‡ for the protocol band.
- Whether the registered set may depend on free ages, if no common-age scenario fits.
- Whether to commit `lab_observables.json` and `results/05/` under the same confidentiality rules as `results/03/`.

## 5. Risks and how the plan detects them

| Risk | Detection |
|---|---|
| Too few data: 7 informative values | Synthetic recovery in phase 3; flat profiles reported as "not determined" |
| Unknown times decide the result | Four scenarios; each result tagged as holding in all or in some |
| A sample outside the model (salt) | Phase 1 audit; fit without it |
| True structure not in the set | No structure fits in any scenario. That is then the result; acid catalysis is next |
| Errors too small (T₁, NOE) | Floor at zero and doubled |
| Plateau read as equilibrium when the water ran out | Profile of the available-water fraction; indistinguishability test |
| Computed covariance is statistical only | Refit with it halved and doubled |
| Solver crawl, noisy derivatives | Call budget, search box, tolerance check |
| Many choices made after seeing results | Fit rule, selection rule and scenarios are fixed in the notebook header first |

**Measurements that would close each blind spot:**

| Blind spot | Measurement |
|---|---|
| Water order | The mixing dates, or one time series at two water contents |
| Equilibrium or exhausted water | Karl Fischer on the plateau sample, or adding water to it |
| ΔS‡ | The same step at a second temperature |
| R9 | The assignment of the 19.08 ppm line |
| Transfer | Dry TMSPA + TMSOH followed in time |

## 6. Outlook, left out on purpose

- Rate-law variants: acid catalysis, EC + water.
- R4 from ¹H after referencing.
- Split transfer.
- R9 and gas.
- Sampling.

## Housekeeping

- My timing and check scripts are in the session scratchpad, outside the repository.
- I saved one note for future sessions: which conda env to use, and that the raw spectra are only on the Mac.
- The Google Drive and Microsoft 365 connectors need authorising in your claude.ai connector settings before they can be used here.
