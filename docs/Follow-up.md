# Follow-up

|                          |                                                                                     |
| ------------------------ | ----------------------------------------------------------------------------------- |
| **Status**         | Draft                                                                               |
| **Classification** | CONFIDENTIAL (see CLASSIFICATION.md)                                                |
| **Source**         | Y. Alcaraz Galván                                                                  |
| **Scope**          | The problem of the project and its pipeline, in schematic form.                     |
| **Builds on**      | `theory-multiscale_microkinetics.md` (the detailed version), NB01`01_multiscale_microkinetics_theory.ipynb`, `kinetics/fitting/residuals.py` |

---

## 1. General view

- **Goal:** predict how fast TMSPA removes water from an EC electrolyte, and which products form.
- **What we have:** computed energies of every species (DFT + MD). They give the reaction energies ΔG_rxn.
- **What we lack:** the barriers. No transition state is computed, so each reaction family has one unknown intrinsic barrier g.
- **How we get them:** simulate the NMR spectrum the model predicts and adjust g until it matches the lab spectra.
- **Why it is hard:** the lab spectra are snapshots. No sample was followed in time, and the time since mixing is not recorded.

### 1.1 General Flowsheet

```
DFT + MD  →  G of each species  →  ΔG_rxn, K_eq  →  ΔG‡, k  →  C(t)  →  NMR  ⇄  lab NMR
energies     in solution           per reaction     rates      reactor   simulated   measured
```

### 1.2 Simulation scheme

Blocks are in bold; the rows marked ↓ say what passes from one block to the next and how.

|                         | What                                                                                                                                          | How                                                                                                                                                                                                                                          |
| ----------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **DFT + MD**      | Energy of each of the 11 species, alone                                                                                                       | 1. Gas phase: ωB97M-V/def2-TZVPD, G at 298 K and 1 bar.<br />2. Solvation: MACE MD of the solute among 14 EC molecules, ΔE_solv ± standard error.                                                                                         |
| ↓                      | G_gas and ΔE_solv per species                                                                                                                | Thermodynamic cycle: G_sol = G_gas + ΔE_solv.                                                                                                                                                                                               |
| **G in solution** | Free energy of each species in EC                                                                                                             | Known at 298 K only (no frequencies stored). ΔE_solv has no entropy.                                                                                                                                                                        |
| ↓                      | G_sol per species                                                                                                                             | Stoichiometry:<br />ΔG_rxn = Σ ν G_sol; K_eq = exp(−ΔG_rxn / RT).                                                                                                                                                                       |
| **ΔG_rxn, K_eq** | Driving force of the 9 reactions, in 4 families: hydrolysis R1–R3, condensation R4, transfer R5–R7, solvent attack R8–R9                   | Standard error 0.23–0.44 eV, larger than most ΔG_rxn.<br />Transfer = hydrolysis + condensation, so the cycles close exactly.                                                                                                              |
| ↓                      | ΔG_rxn per reaction, plusg per family (unknown, fitted). <br />The model is specified here: which barrier relation, and how many separate g. | No transition state is computed, so the barrier is estimated from ΔG_rxn: more downhill, lower barrier.<br />1. Marcus: ΔG‡ = g (1 + ΔG_rxn / 4g)². <br />2. Capped BEP (Peter): ΔG‡ = E₀ if downhill, E₀ + 0.5 ΔG_rxn if uphill. |
| **ΔG‡, k**      | Barrier and rate constants of each reaction                                                                                                   | Eyring: k_f = (k_B T / h) exp(−ΔG‡ / k_B T). Detailed balance: k_r = k_f / K_eq. ΔS‡ = 0.                                                                                                                                               |
| ↓                      | k_f and k_r at each temperature, plus the recipe and the temperature history of the sample                                                    | Mass action: dC/dt = S · r(C, T).<br />S: stoichiometric matrix. r_j = k_f,j Π C_reactants − k_r,j Π C_products: net rate of reaction j.                                                                                                |
| **C(t), reactor** | Concentration of every species against time                                                                                                   | Closed, well-mixed tank. EC held constant. Stiff ODE solver, one segment per temperature.                                                                                                                                                    |
| ↓                      | C of each species at the time of each spectrum                                                                                                | One peak per species: area ∝ (nuclei per molecule) × C.                                                                                                                                                                                    |
| **Simulated NMR** | Area share of each peak, per nucleus                                                                                                          | ³¹P: TMSPA, BMSPA, MMSPA, H₃PO₄. ¹³C (Si–CH₃): TMSOH, HMDSO, TMSOEG. Share = peak area / sum of the areas.                                                                                                                           |
| ⇄                      | Predicted share against measured share                                                                                                        | Residual z = (predicted − measured) / error. The fit moves g (and, in the larger models, ΔG_rxn of R1–R4) to make the z small.                                                                                                            |
| **Lab NMR**       | Measured area shares with their errors                                                                                                        | Raw spectra → area in a ppm window around each peak → share. Error = noise ⊕ baseline spread ⊕ declared floor.                                                                                                                           |

### 3. Th1. The probleme pipeline in detail

**Positions and amounts are separate.** The DFT shieldings give the peak positions (δ = σ_ref − σ), used to assign the peaks. The comparison with the lab uses only the shares, with windows set around the observed peaks.
The


## 2. Index

- **1. DFT + MD**
  - 1.1 Gas-phase free energy of each species
  - 1.2 Solvation energy in EC
  - 1.3 Free energy in solution
- **2. Kinetics**
  - 2.1 Driving force: ΔG_rxn, K_eq, and check if closed cycles
  - 2.2 Models
    - 2.2.1 Barrier relation (Marcus, capped BEP)
    - 2.2.2 M0: one barrier for all reactions (`peter_reference`)
    - 2.2.3 M1: one barrier per family (`level1`)
    - 2.2.4 M1-split: R1 separate from R2 and R3
    - 2.2.5 M3 and M3-split: reaction energies of R1–R4 free
  - 2.3 Eyring and detailed balance
  - 2.4 Conclusions of Part 2
- **3. Mass action**
  - 3.1 Rate of one reaction
  - 3.2 Mass balance of every species: dC/dt
  - 3.3 Reactor: closed, well mixed, EC constant
  - 3.4 Inputs of a run: recipe and temperature history
  - 3.5 Numerical solution: stiffness, solver, tolerances
- **4. NMR**
  - 4.1 Peaks
    - 4.1.1 Simulated: from concentrations to area shares
    - 4.1.2 Lab data: from raw spectrum to area shares
  - 4.2 NMR fundamentals for ³¹P and ¹³C: when an area is an amount
- **5. Model fitting**
  - 5.1 Residual: predicted against measured, in units of the error
  - 5.2 Fitting variables: barriers, reaction energies, sample age, available water
  - 5.3 Fit rule: when a model counts as fitting
  - 5.4 Unknown sample ages: the four scenarios
  - 5.5 Search: how the best fit is found
  - 5.6 Limitations
- **6. Comparison**
  - 6.1 Direct comparison: predicted and measured shares, sample by sample
  - 6.2 Model quality: which models fit, what is determined, which sample decides each value
  - 6.3 Errors and noise: the error of a share and how the result depends on it
- **7. Conclusions**

---

## Part 1. DFT + MD

**What this part delivers:** one free energy in EC solution (G_sol), for the 11 species of the network. Everything after it is built from these 11 numbers. Full detail: theory, [Block 1](theory-multiscale_microkinetics.md#block-1-first-principles-quantum-chemistry--solvation-database-27-species), [Block 2](theory-multiscale_microkinetics.md#block-2-statistical-mechanics-engine--grimmes-quasi-rrho) and [Block 3](theory-multiscale_microkinetics.md#block-3-solvation-thermochemistry--solvation-driven-reaction-free-energies).

```
DFT (one molecule, gas)  ──►  G_gas  ───┐
                                        ├──►  G_sol = G_gas + ΔE_solv
MD (molecule in liquid EC) ──► ΔE_solv ─┘
```

### 1.1 Gas-phase free energy of each species

- **What:** the free energy of one isolated molecule at 298.15 K and 1 bar, G_gas.
- **How:** DFT at the ωB97M-V/def2-TZVPD level gives the electronic energy. Thermal motion (translation, rotation, vibration) is added on top ([how each term is built](theory-multiscale_microkinetics.md#1-the-finite-temperature-bridge)): 

    $G_{gas} = E_{SCF} + ZPE + H_{thermal} − T·S$

- **What we hold:** only the final G and H at 298.15 K. The vibrational frequencies are not stored.
- **Consequence:** G_gas cannot be recomputed at another temperature. The 298 K value is used at every temperature, including the 22 °C and 80 °C of the lab samples.
- **Size of the entropy:** S = (H − G)/T gives 699 J/mol/K for TMSPA and 189 for water. Large molecules have large entropies, but reactants and products move together, so reaction energies would change by at most 0.05 eV between 20 and 80 °C.

### 1.2 Solvation energy in EC

- **What:** the energy change when the gas molecule is placed in liquid EC, ΔE_solv.
- **How:** molecular dynamics with the MACE-OMol potential, one solute among 14 EC molecules: ΔE_solv = ⟨E of the solution box⟩ − (14/15)·⟨E of the pure-EC box⟩ − ⟨E of the solute in gas⟩.
- **Standard error:** each ⟨E⟩ is an average over a noisy MD run, so ΔE_solv carries an error of 0.21–0.36 eV per species ([how it carries into a reaction](theory-multiscale_microkinetics.md#3-solvation-uncertainty-provisional)).

| Species | ΔE_solv (eV) | Standard error (eV) |
|---|---|---|
| TMSPA | −1.54 | 0.25 |
| BMSPA | −1.98 | 0.23 |
| MMSPA | −2.12 | 0.24 |
| H₃PO₄ | −2.31 | 0.22 |
| H₂O | −1.41 | 0.22 |
| TMSOH | −1.59 | 0.26 |
| HMDSO | −1.53 | 0.21 |
| EC | −0.56 | 0.29 |
| TMSOEG | −1.45 | 0.21 |
| TMSOdiEG | −1.73 | 0.32 |
| CO₂ | −0.42 | 0.36 |

- **Reading:** each TMS group replaced by an OH makes the phosphate more strongly solvated (TMSPA −1.54 → H₃PO₄ −2.31 eV). This is what drives hydrolysis in solution.

**The solvent mismatch: computed in pure EC, measured in EC/DEC**

- **Evidence:** the lab file names and the ¹H spectra show EC and DEC. The 1:1 ratio is taken from the paper, not from the files.
- **Concentration: handled.** The lab samples are simulated with EC at 7.1 M, the real EC content of the mixture. DEC does not react and is left out.
- **Energies: not handled.** ΔE_solv is computed for a molecule surrounded by EC only. EC/DEC is a less polar liquid, so polar species (H₃PO₄, MMSPA, water) are really stabilised less than computed. Hydrolysis is then less downhill than the model believes. No concentration conversion corrects this.
- **Where it can show:** only where a sample is seen at rest, because the resting composition depends on the reaction energy and not on the barrier. That is the 2 % water sample, which stops at 66 % H₃PO₄ and 32 % MMSPA where the computed energies give 100 % H₃PO₄.
- **Where it does not show:** solvent attack and condensation are caught on the way, far from equilibrium. Only their barriers matter, and those are fitted to the data.
- **It is one of three explanations of that plateau.** The other two are the error of the MD energies themselves (about 0.25 eV) and part of the added water not being available. The data cannot separate the three ([NB05 §5.4](../notebooks/05_first_fit.ipynb)).
- **It does not explain** why TMSPA survives at 0.5 % water and is gone at 2 %. That is a question of rates and of the unknown sample ages.
- **What would close it:** ΔE_solv recomputed in an EC/DEC box, or lab samples in pure EC.

### 1.3 Free energy in solution

- **What:** G_sol = G_gas + ΔE_solv, per species.
- **Standard-state term:** moving the reference from a 1 bar gas to a 1 M solution adds +7.96 kJ/mol (0.08 eV) at 298 K, the same for every species. Every reaction of the network has two molecules on each side, so it cancels and is left out ([derivation](theory-multiscale_microkinetics.md#1-the-liquid-phase-thermodynamic-cycle)).
- **Solvation decides the chemistry.** For R1 (TMSPA + H₂O → BMSPA + TMSOH), ΔG is +0.14 eV in the gas phase and −0.47 eV in EC. **Seven of the nine reactions change sign** ([why solvation does this](theory-multiscale_microkinetics.md#2-the-driving-force-solvation-flips-6-out-of-9-reactions)).
- **Solvation also decides the uncertainty.** The standard errors of 1.2 carry into every reaction energy: 0.23–0.44 eV, larger than most of the reaction energies themselves (Part 2.1).

**Limits to keep in mind**

| Limit | Why it matters later |
|---|---|
| G is known at 298 K only | Reaction energies do not change with temperature |
| ΔE_solv is an energy, not a free energy | The entropy of solvation is missing |
| The solvent is pure EC | The lab samples are in EC/DEC 1:1 |
| Standard error of 0.2–0.4 eV | The fit of notebook 05 has to let reaction energies move (Part 5.2) |

---

## Part 2. Kinetics

**What this part delivers:** the forward and reverse rate constant of each of the nine reactions, at any temperature. Full detail: theory, [Block 4](theory-multiscale_microkinetics.md#block-4-reaction-network-thermodynamics--the-wegscheider-consistency) and [Block 5](theory-multiscale_microkinetics.md#block-5-microkinetics--rate-constants-engine-with-strict-detailed-balance).

```
G_sol (Part 1)  ──►  ΔG_rxn, K_eq       ──►  ΔG‡        ──►  k_f, k_r
                     2.1 driving force       2.2 model       2.3 Eyring
```

### 2.1 Driving force: ΔG_rxn, K_eq, and closed cycles

- **What:** ΔG_rxn is the free energy of the products minus that of the reactants. Negative means downhill. K_eq says where the reaction comes to rest.
- **How:** nothing new is computed. The G_sol of Part 1 are added and subtracted:

    $\Delta G_{rxn} = \sum G_{sol}(\text{products}) - \sum G_{sol}(\text{reactants}), \qquad K_{eq} = \exp(-\Delta G_{rxn}/RT)$

- **Scale:** 0.06 eV in ΔG_rxn is a factor of 10 in K_eq at room temperature.

| | Family | Reaction | ΔG_rxn (eV) | Standard error (eV) | K_eq at 25 °C | Sign certain? |
|---|---|---|---|---|---|---|
| R1 | hydrolysis | TMSPA + H₂O ⇌ BMSPA + TMSOH | −0.47 | 0.26 | 1 × 10⁸ | yes |
| R2 | hydrolysis | BMSPA + H₂O ⇌ MMSPA + TMSOH | −0.16 | 0.25 | 470 | no |
| R3 | hydrolysis | MMSPA + H₂O ⇌ H₃PO₄ + TMSOH | −0.18 | 0.24 | 1000 | no |
| R4 | condensation | 2 TMSOH ⇌ HMDSO + H₂O | +0.10 | 0.34 | 0.02 | no |
| R5 | transfer | TMSPA + TMSOH ⇌ BMSPA + HMDSO | −0.37 | 0.25 | 2 × 10⁶ | yes |
| R6 | transfer | BMSPA + TMSOH ⇌ MMSPA + HMDSO | −0.06 | 0.23 | 9 | no |
| R7 | transfer | MMSPA + TMSOH ⇌ H₃PO₄ + HMDSO | −0.08 | 0.23 | 19 | no |
| R8 | solvent attack | EC + TMSOH ⇌ TMSOEG + CO₂ | −0.05 | 0.40 | 6 | no |
| R9 | solvent attack | EC + TMSOEG ⇌ TMSOdiEG + CO₂ | −0.47 | 0.44 | 8 × 10⁷ | yes |

**The four families** ([chemistry of each](theory-multiscale_microkinetics.md#1-the-4-mechanistic-reaction-classes))

- **Hydrolysis (R1–R3):** water removes one TMS group from the phosphate and releases TMSOH. Three steps, from TMSPA down to H₃PO₄.
- **Condensation (R4):** two TMSOH join into HMDSO and give one water back.
- **Transfer (R5–R7):** TMSOH removes the TMS group in place of water. The product is HMDSO and no water is used.
- **Solvent attack (R8–R9):** TMSOH, and then its product, opens an EC ring and releases CO₂.

**Reading the table**

- **Only three signs are certain:** R1, R5 and R9, where |ΔG_rxn| is larger than its standard error. For the other six the calculation cannot say whether the reaction is downhill or uphill.
- **The error is large in K_eq.** A standard error of 0.25 eV is a factor of about 10 000 in K_eq. For R3 the computed K_eq of 1000 could as well be 0.1 or 10⁷.
- **Consequence:** the computed energies say R2 and R3 go to completion. The 2 % water sample stops at 66 % H₃PO₄ and 32 % MMSPA. An energy of R3 near zero instead of −0.18 eV would explain it, and that is within one standard error.

**Closed cycles** ([the Wegscheider condition](theory-multiscale_microkinetics.md#2-intuitive-meaning-of-the-wegscheider-consistency-condition-cyclic-detailed-balance))

- **What:** transfer gives the same products as hydrolysis followed by condensation. Adding R1 and R4 and cancelling what appears on both sides gives R5. The same holds for R6 = R2 + R4 and R7 = R3 + R4.
- **Condition:** two routes between the same molecules must have the same total energy, so ΔG(R5) = ΔG(R1) + ΔG(R4): −0.47 + 0.10 = −0.37 eV.
- **Check:** the three sums close to 10⁻¹⁰ kJ/mol, which is rounding error. They must, because every ΔG_rxn is built from the same eleven G_sol.
- **Why it matters for the fit:** the nine reaction energies are not independent. Only R1–R4 can be moved; R5–R7 follow. The fit moves them through the species energies, so the cycles stay closed.

**Limits to keep in mind**

| Limit | Why it matters later |
|---|---|
| ΔG_rxn is fixed at its 298 K value | K_eq changes with temperature only through the 1/RT in the exponent. The neglected change is at most 0.05 eV between 20 and 80 °C |
| Six of nine signs are not certain | The fit has to let the energies of R1–R4 move (Part 5.2) |
| Energies are for pure EC | The lab samples are in EC/DEC (Part 1.2) |

### 2.2 Models

- **What a model is:** the rule that turns each ΔG_rxn into a barrier ΔG‡. It has two parts: the barrier relation (2.2.1) and how many separate intrinsic barriers g there are.
- **Why one is needed:** no transition state is computed, so no barrier is known. The g are the unknowns that the lab data must fix.
- **Source of the fit results:** every fitted value and every "fits / fails" in this document comes from notebook 05, cited as NB05 §n. Its numbers are in `notebooks/results/05/`.
- **Named model or structure:** a named model has values for its g (`peter_reference`, `level1`). A structure is the same rule with the values left free for the fit (M0, M1, …).
- **The ladder:** each structure adds freedom to the one before, so a simpler one is a special case of the next.

| Structure | Free parameters | What is free | Reaction energies |
|---|---|---|---|
| M0 | 1 | one g for all nine reactions | as computed |
| M1 | 4 | one g per family | as computed |
| M1-split | 5 | M1, with R1 separate from R2 and R3 | as computed |
| M3 | 8 | M1, plus ΔG_rxn of R1–R4 | free, held near the computed values |
| M3-split | 9 | M1-split, plus ΔG_rxn of R1–R4 | free, held near the computed values |

- **Naming caution:** these are the definitions of notebook 05. `03_feasible_region` uses the same names with other contents: there M1-split has one g per hydrolysis step (6 parameters) and M3 frees only R2 and R3.

#### 2.2.1 Barrier relation (Marcus, capped BEP)

- **Idea:** within a family, a more downhill reaction has a lower barrier. g is the barrier the reaction would have at ΔG_rxn = 0. Detail: [theory, Block 5.1](theory-multiscale_microkinetics.md#51-models).
- **Marcus:** a smooth curve.

    $\Delta G^\ddagger = g\,\left(1 + \dfrac{\Delta G_{rxn}}{4g}\right)^2 \approx g + \dfrac{\Delta G_{rxn}}{2}$

- **Capped BEP (Peter's form):** a line with a floor.

    $\Delta G^\ddagger = \max(E_0,\; E_0 + 0.5\,\Delta G_{rxn})$

- **Example, both with 0.80 eV:**

| | ΔG_rxn (eV) | ΔG‡ Marcus (eV) | ΔG‡ capped BEP (eV) |
|---|---|---|---|
| R1 | −0.47 | 0.58 | 0.80 |
| R2 | −0.16 | 0.72 | 0.80 |
| R3 | −0.18 | 0.71 | 0.80 |
| R4 | +0.10 | 0.85 | 0.85 |

- Marcus ranks the reactions of a family: One g gives three different barriers for R1, R2 and R3, the most downhill being the fastest.
- Capped BEP does not: every downhill reaction gets exactly E₀. Only the uphill R4 feels its ΔG_rxn.

- **Half of an energy error goes into the barrier.** Because ΔG‡ ≈ g + ΔG_rxn/2, a standard error of 0.25 eV in ΔG_rxn is about 0.12 eV in ΔG‡ for a fixed g. This is why the fit determines ΔG‡ well and g less well.

#### 2.2.2 M0: one barrier for all reactions (`peter_reference`)

- **Assumption:** every reaction of every family has the same intrinsic barrier. One free parameter.
- **Two forms:** M0-BEP, which is the structure of `peter_reference`, and M0-Marcus.
- **`peter_reference`:** capped BEP with E₀ = 1.15 eV. All barriers are 1.15 eV, except R4 at 1.20 eV.
- **Where 1.15 eV comes from:** it spreads the consumption of TMSPA over the 20–80 °C holds of Peter's protocol. It is a design value for planning the measurement, not a fit.
- **Against the lab data, at 1.15 eV:** solvent attack is too fast and hydrolysis too slow at the same time ([NB05 §2](../notebooks/05_first_fit.ipynb)).

| Sample | Predicted | Measured |
|---|---|---|
| TMSOH in EC, one week at room temperature | 62 % has opened an EC ring | none (below 2 %) |
| 2 % water, after one day | 84 % of the TMSPA left | none left |

- **Best single value** ([NB05 §4.1](../notebooks/05_first_fit.ipynb))**:** 1.28 eV in the capped-BEP form and 1.46–1.47 eV in the Marcus form. At that value nothing reacts at room temperature, so every phosphate sample is predicted as 100 % TMSPA. The fit gives up hydrolysis to keep TMSOH from opening EC.
- **Why no value works:** solvent attack needs about 1.29 eV, so that TMSOH is left alone at room temperature. Consuming TMSPA within days at room temperature needs about 1.1 eV or less. One number cannot be both.
- **Verdict:** the structure is rejected, not the value. This says nothing against the purpose of `peter_reference`, which was never meant to describe these samples.

#### 2.2.3 M1: one barrier per family (`level1`)

- **Assumption:** each of the four families has its own g, with Marcus. Four free parameters. Inside a family the reactions differ only through their ΔG_rxn, which stays as computed.
- **`level1`:** the named model with this structure. Where its values come from: [theory, Block 5.3](theory-multiscale_microkinetics.md#53-reaction-family-barrier-calibration--finding-8-resolution).

| Family | g in `level1` (eV) | Source |
|---|---|---|
| Hydrolysis | 0.80 | placeholder, no data behind it |
| Transfer | 0.80 | placeholder, below the upper bound from the paper |
| Condensation | 1.30 | from the paper: no HMDSO from TMSOH alone at 80 °C |
| Solvent attack | 1.32 | from the paper: EC opens at 80 °C and not at room temperature |

- **What `level1` gets right as it is** ([NB05 §2](../notebooks/05_first_fit.ipynb))**:** the TMSOH-only samples. Its condensation and solvent-attack values came from the paper and also describe the lab spectra.
- **What it gets wrong as it is:** with 0.80 eV, hydrolysis and transfer finish in under a second. The 0.5 % water sample is predicted as 90 % H₃PO₄, where 93 % TMSPA is measured.
- **With the four g fitted** ([NB05 §4](../notebooks/05_first_fit.ipynb))**:** hydrolysis rises to 1.3–1.4 eV (0.94 eV if every sample is taken as 1 h old) and solvent attack stays at 1.31 eV. The structure still fails in every scenario.
- **Why it fails:** the computed energies make R2 and R3 run to the end, so the 2 % water sample is predicted as 100 % H₃PO₄. The measurement is 66 % H₃PO₄ and 32 % MMSPA, unchanged over 145 h.

--> A barrier sets how fast a reaction goes, not where it stops. A sample at rest cannot be fitted by changing barriers. The problem is in the reaction energies, which is what M3 frees (2.2.5).

#### 2.2.4 M1-split: R1 separate from R2 and R3

- **Assumption:** as M1, but the first hydrolysis (R1) has its own g. R2 and R3 share a second one. Five free parameters. Reaction energies stay as computed.
- **Why split hydrolysis:** with one g, Marcus makes R1 the fastest of the three, because it is the most downhill (at g = 0.80 eV: 0.58 eV for R1 against 0.72 eV for R2). The samples ask for the opposite. TMSPA survives in the 0.5 % water sample, so R1 is slow. The 2 % water sample is mostly H₃PO₄, so R2 and R3 are not.
- **Why R2 and R3 stay together:** only the 2 % water sample shows hydrolysis past BMSPA. One sample cannot fix two separate barriers.
**The fit** ([NB05 §4](../notebooks/05_first_fit.ipynb))

- **What was done:** an age is assumed for the samples whose mixing time is unknown. The five barriers are then adjusted until the simulated NMR shares are as close as possible to the measured ones. This is repeated for each assumed age, so each row of the table is a separate fit.

| Assumed age of the samples | ΔG‡(R1) (eV) | ΔG‡(R2) (eV) | ΔG‡(R3) (eV) | Total misfit χ² | Worst miss |
|---|---|---|---|---|---|
| 1 h | 0.91 | 0.85 | 0.84 | 1021 | 20.1 |
| 1 day | 1.23 | 0.96 | 0.95 | 301 | 16.0 |
| 7 days | 1.18 | 1.03 | 1.02 | 291 | 15.4 |
| free (each age chosen by the fit) | 1.09 | 1.04 | 1.04 | 282 | 15.0 |

**How to read the table**

- **Barriers:** the values that give the best match under that age. R1 comes out above R2 and R3 in every row, which is what the split was for.
- **A miss:** how far a prediction is from the measurement, counted in standard errors of that measurement: (predicted − measured) / standard error. Example, 1 day: MMSPA in the 2 % water sample is measured at 0.316 ± 0.020 and predicted at about 0, a miss of 16.
- **Size of a miss:** 1 or 2 is normal noise. Above 3 the model is wrong about that value.
- **Worst miss:** the largest of the 45 misses. The rule fixed before fitting: a structure fits only if it is 3 or less. Here it is 15–20, so M1-split fails at every age.
- **χ²:** all the misses in one number, the sum of their squares. It is what the fit makes as small as it can. The one structure that fits reaches about 20; here the single MMSPA miss already gives 16² ≈ 256.

**Why the table gives ΔG‡ and not g**

- **ΔG‡ is what sets the rate.** The data see how fast each reaction goes, and that depends on its actual barrier ΔG‡.
- **g is one step removed.** It is obtained from ΔG‡ through the Marcus relation and the computed ΔG_rxn, so it inherits the error of ΔG_rxn.
- **g cannot be compared between reactions.** At 1 day the fitted g are 1.46 eV for R1 and 1.04 eV for R2–R3, a gap of 0.42 eV. The barriers are 1.23 and 0.96 eV, a gap of 0.27 eV. R1 is much more downhill, so Marcus takes more off its g. Only the ΔG‡ say directly which reaction is slower, and by how much.
- **The split helps very little.** Against M1, χ² falls from 343 to 301 in the middle scenario. The largest miss stays at 15–16 standard errors.
- **Why it still fails:** the same reason as M1. The largest miss is MMSPA in the 2 % water sample. The computed energies send R2 and R3 to the end, and no barrier can stop them partway.
- **Test with part of the water unavailable** ([NB05 §5.4](../notebooks/05_first_fit.ipynb))**:** if only 21 % of the added water takes part, χ² falls to about 80. The split between BMSPA and MMSPA in the 2 % sample is still missed by 7 standard errors. Running out of water, with the computed energies, is not enough.

--> Splitting the barrier is needed, but it only works once the reaction energies can move. That combination is M3-split (2.2.5), the only structure that fits.

#### 2.2.5 M3 and M3-split: reaction energies of R1–R4 free

- **Assumption:** the computed reaction energies may be wrong by about their standard error, so the fit may move them. M3 is M1 plus the four energies of R1–R4 (8 free parameters). M3-split is M1-split plus the same four (9 free parameters).
- **Why only R1–R4:** the transfer energies follow from them through the closed cycles (2.1). R8 and R9 keep their computed values, because no sample is seen at rest for them.
- **How an energy is held:** it is not free to go anywhere. The computed value counts as one more measurement, "ΔG_rxn = computed ± standard error", and a shift is a miss like any other. Example: moving R3 by +0.21 eV, with a standard error of 0.24 eV, costs a miss of 0.9. This adds four values to the 45, so χ² is over 49.
- **What it buys:** the model can now stop a reaction partway. That is what the 2 % water sample needs.

**The fit** ([NB05 §4](../notebooks/05_first_fit.ipynb))

| Assumed age of the samples | M3: χ² | M3: worst miss | M3-split: χ² | M3-split: worst miss |
|---|---|---|---|---|
| 1 h | 536 | 11.6 | 535 | 11.6 |
| 1 day | 294 | 16.1 | **20.4** | **1.9** |
| 7 days | 60 | 4.1 | **19.8** | **1.9** |
| free | 27 | 3.1 | **18.6** | **1.9** |

- **M3-split is the only structure that fits,** and it fits at three of the four ages (bold: worst miss below 3).
- **M3 fails,** by very little with free ages (3.06 against the limit of 3). One hydrolysis barrier cannot serve R1 and R2–R3 at once.
- **Both changes are needed together:** the split alone fails (2.2.4), free energies alone fail (M3).
- **At 1 h nothing fits:** there is not enough time for the samples to react as observed.

**What M3-split finds** ([NB05 §4.2](../notebooks/05_first_fit.ipynb))

| | Computed ΔG_rxn (eV) | Fitted ΔG_rxn (eV) | ΔG‡, 1 day (eV) | ΔG‡, 7 days (eV) | ΔG‡, free (eV) |
|---|---|---|---|---|---|
| R1, first hydrolysis | −0.47 | −0.12 to −0.04 | 1.45 | 1.39 | 1.09 |
| R2, second hydrolysis | −0.16 | −0.04 to +0.05 | 0.95 | 1.01 | 1.04 |
| R3, third hydrolysis | −0.18 | +0.01 to +0.09 | 0.97 | 1.03 | 1.06 |
| R5, transfer | −0.37 | −0.22 to −0.10 | 0.58 | 0.61 | 0.89 |
| R4, condensation | +0.10 | −0.18 to +0.01 | 1.24 | 1.24 | 1.24 |
| R8, solvent attack | −0.05 | −0.05 (not moved) | 1.29 | 1.29 | 1.29 |

- **Energies:** R2 and R3 move up to about zero. A reaction with ΔG_rxn near zero stops partway, which gives the resting composition of the 2 % water sample. All shifts stay within 1.7 standard errors of the computed values.
- **The same at every age:** solvent attack (1.29 eV), condensation (1.24 eV), and the second and third hydrolysis (0.95–1.06 eV).
- **Different at every age:** the first hydrolysis (1.09 to 1.45 eV) and transfer (0.58 to 0.89 eV). These two decide how TMSPA is consumed, and the unknown ages decide them.

**Why this is not yet the model**

- **The test is weak.** Nine parameters are matched to about 17 independent compositions. A structure this flexible is hard to contradict.
- **The 1-day and 7-day fits rest on a fragile mechanism** ([NB05 §4.3](../notebooks/05_first_fit.ipynb)). In them TMSPA does not react for hours and then disappears within one. The timing depends on a trace of TMSOH at mixing that nobody measured. With 1 nM present, the 1-day fit goes from χ² 20 to 373 (`trace_sensitivity.csv`, I did´t check it). Fitted again with a trace it can still fit, but with other barriers: ΔG‡(R1) anywhere from 1.24 to 1.69 eV. The fit with free ages (should we ask) does not depend on a trace.
- **The plateau has a second explanation** ([NB05 §5.4](../notebooks/05_first_fit.ipynb)). If only a fifth of the added water was available, the fit is as good with much smaller energy shifts.
- **No model was registered.** Parts 5 and 6 cover what can and cannot be trusted.

--> M3-split shows what a working model needs: R1 slower than R2 and R3, and R2 and R3 close to thermoneutral. It does not yet say how TMSPA is consumed.

### 2.3 Eyring and detailed balance

- **What:** the last step of Part 2. Each barrier becomes a forward rate constant, and each equilibrium constant then fixes the reverse one. Detail: [theory, Block 5](theory-multiscale_microkinetics.md#block-5-microkinetics--rate-constants-engine-with-strict-detailed-balance).

**Forward: Eyring**

  $k_f = \dfrac{k_B T}{h}\,\exp\left(-\dfrac{\Delta G^\ddagger}{k_B T}\right)$

- **Reading:** k_B T / h is how often the reactants attempt the reaction, 6 × 10¹² per second at room temperature. The exponential is the fraction of attempts with enough energy to pass the barrier.
- **Units:** M⁻¹ s⁻¹, because every reaction here is between two molecules.
- **Rule:** at room temperature, +0.06 eV in the barrier makes the reaction 10 times slower (at 80 °C it takes 0.07 eV).

| ΔG‡ (eV) | k_f at 22 °C (M⁻¹ s⁻¹) | Half-life at 22 °C, partner at 1 M | Half-life at 80 °C, partner at 1 M |
|---|---|---|---|
| 0.60 | 360 | 2 ms | below 1 ms |
| 0.80 | 0.14 | 5 s | 25 ms |
| 1.00 | 5 × 10⁻⁵ | 3.5 h | 17 s |
| 1.10 | 1 × 10⁻⁶ | 7 days | 8 min |
| 1.29 | 6 × 10⁻¹⁰ | 35 years | 3 days |

- **Why the fit can find barriers at all:** the whole range from "instant" to "never" lies between 0.6 and 1.3 eV. A sample that has reacted partly at a known time pins its barrier to a few hundredths of an eV.
- **Check with solvent attack:** at 1.29 eV with EC at 7.1 M, the half-life of TMSOH is about 5 years at room temperature and about 9 h at 80 °C. This matches what is seen: no ring-opening after a week at room temperature, and about half of the TMSOH opened after 8 h at 80 °C.

**Reverse: detailed balance**

$k_r = \dfrac{k_f}{K_{eq}}, \qquad \Delta G^\ddagger_r = \Delta G^\ddagger_f - \Delta G_{rxn}$

- **Reading:** the reverse reaction climbs the same hill from the other side. If the products lie 0.18 eV below the reactants, the way back is 0.18 eV higher.
- **Why it is imposed:** it makes every reaction stop exactly at its K_eq, and it keeps the closed cycles of 2.1 closed in the rates as well.
- **Example with R3:** with the computed −0.18 eV, the reverse is 1000 times slower than the forward and R3 runs to the end. With the fitted value near zero, the two are about equal and R3 stops partway.

**Temperature**

- **What changes with temperature:** only the T in the Eyring formula. The barrier itself is kept the same at every temperature, which is the assumption ΔS‡ = 0 (no activation entropy).
- **Effect:** from 22 to 80 °C a 1.00 eV barrier becomes about 700 times faster, and a 1.29 eV barrier about 4600 times faster.
- **Where each barrier was seen:** hydrolysis and transfer at 22 °C only; condensation and solvent attack at 80 °C only. Using a barrier at the other temperature is an extrapolation.
- **Size of the risk:** an activation entropy of 100 J mol⁻¹ K⁻¹ changes a barrier by 0.06 eV over those 58 K, a factor of 10 in rate. The data cannot measure it, because no reaction was seen at two temperatures ([NB05 §7.2](../notebooks/05_first_fit.ipynb)).

### 2.4 Conclusions of Part 2

**What Part 2 does**

- **Driving force (2.1):** the energies of Part 1 give ΔG_rxn and K_eq for nine reactions. They say where each reaction stops, but only three of the nine signs are certain.
- **Model (2.2):** no barrier is computed. A model is a barrier relation plus a set of unknown intrinsic barriers g, which the lab data must fix.
- **Rate constants (2.3):** Eyring turns each barrier into k_f, and detailed balance gives k_r.

**What the ladder of models showed** ([NB05 §4](../notebooks/05_first_fit.ipynb))

| Structure | Free parameters | Fits? | Why not |
|---|---|---|---|
| M0 | 1 | no | one barrier cannot be slow for solvent attack and fast for hydrolysis |
| M1 | 4 | no | computed energies send hydrolysis to the end; the 2 % sample stops partway |
| M1-split | 5 | no | same reason; separating R1 does not help on its own |
| M3 | 8 | no, by little | one hydrolysis barrier cannot serve R1 and R2–R3 |
| M3-split | 9 | yes, at 1 day, 7 days and free ages | — |

**What can be taken as known** ([NB05 §8](../notebooks/05_first_fit.ipynb))

- **Solvent attack:** ΔG‡(R8) = 1.29 eV at 80 °C.
- **Condensation:** ΔG‡(R4) = 1.24 eV at 80 °C.
- **Second and third hydrolysis:** ΔG‡(R2) and ΔG‡(R3) between 0.92 and 1.07 eV at 22 °C.
- **Energies of R2 and R3:** above the computed values, close to zero, if all the added water was available.
- **Structure:** no single barrier, and no barrier per family with the computed energies, describes the samples.

**What is not known**

- **How TMSPA is consumed.** The first hydrolysis (R1, 1.08–1.65 eV) and transfer (R5, 0.41–0.99 eV) trade off against each other. The sample ages and the TMSOH present at mixing decide.
- **Whether the 2 % sample rests at an equilibrium or ran out of water.**
- **The temperature dependence** of hydrolysis and transfer.
- **Each barrier except solvent attack rests on one sample.**

**Carried into every number above**

- Energies computed in pure EC, samples in EC/DEC (1.2).
- ΔS‡ = 0 (2.3).

--> Part 2 hands to Part 3 one pair of rate constants, k_f and k_r, for each reaction at each temperature.

---

## Part 3. Mass action

**What this part delivers:** the concentration of every species against time, C(t), for one sample. Full detail: theory, [Block 7](theory-multiscale_microkinetics.md#block-7-homogeneous-batch-reactor-dynamics-tank-stiff-ode-integration).

```
k_f, k_r (Part 2)  ──►  rate of each reaction  ──►  dC/dt of each species  ──►  C(t)
recipe, temperature history ──────────────────────────────────────────────────┘
```

### 3.1 Rate of one reaction

- **What:** how many moles per litre react each second, r, in M/s.
- **How (mass action):** the rate constant times the concentrations of the molecules that have to meet, forward minus reverse.

    $r_j = k_{f,j}\prod C_{\text{reactants}} \;-\; k_{r,j}\prod C_{\text{products}}$

- **Example, R1** (TMSPA + H₂O ⇌ BMSPA + TMSOH):

    $r_1 = k_f\,[\text{TMSPA}]\,[\text{H}_2\text{O}] - k_r\,[\text{BMSPA}]\,[\text{TMSOH}]$

- **At mixing** there are no products, so only the forward term acts. As products build up the reverse term grows. The reaction stops when the two terms are equal, which is the equilibrium of 2.1.
- **In numbers:** with ΔG‡(R1) = 1.09 eV (the fit with free ages), k_f = 1.8 × 10⁻⁶ M⁻¹ s⁻¹ at 22 °C. In the 2 % water sample ([H₂O] = 1.05 M) TMSPA then has a half-life of about 4 days. In the 0.5 % water sample ([H₂O] = 0.26 M) it is four times longer, about 17 days.
- **First order in water:** this rate law says four times more water makes R1 exactly four times faster. The two water samples differ by much more than that, which is why their ages matter so much (Part 5.4).

### 3.2 Mass balance of every species: dC/dt

- **What:** each species gains from the reactions that make it and loses to the ones that consume it.
- **Example, TMSOH:** made by the three hydrolysis steps, consumed by condensation (two per event), by the three transfers and by solvent attack.

    $\dfrac{d[\text{TMSOH}]}{dt} = r_1 + r_2 + r_3 - 2r_4 - r_5 - r_6 - r_7 - r_8$

- **All species at once:** one such equation per species, eleven in all. Written together ([theory, Block 7 §1](theory-multiscale_microkinetics.md#1-chemical-reactor-governing-equations)):

    $\dfrac{d\mathbf{C}}{dt} = \mathbf{S}\cdot\mathbf{r}(\mathbf{C}, T)$

    **S** is the stoichiometric matrix: one row per species, one column per reaction, holding how many molecules of that species the reaction makes (+) or consumes (−).
- **The equations are coupled.** TMSOH appears in eight of the nine rates, so no reaction can be followed alone.
- **The loop that matters:** hydrolysis releases TMSOH, and TMSOH consumes TMSPA through transfer (R5). TMSPA can therefore disappear even when its own hydrolysis (R1) is slow. This is why the data cannot separate the barriers of R1 and R5.
- **Built-in check:** no reaction creates or destroys P or Si, so their totals must stay constant during a run.

### 3.3 Reactor: closed, well mixed, EC constant

The sample tube is treated as an ideal batch reactor. Each assumption, and where it is weakest:

| Assumption | Meaning | Where it is weakest |
|---|---|---|
| Closed | nothing enters or leaves; CO₂ stays dissolved | at 80 °C, where ring-opening releases CO₂ |
| Well mixed | one liquid phase, the same composition everywhere | water mixes poorly with DEC; some ¹H spectra show the solvent in two environments |
| Constant temperature per segment | temperature changes are instantaneous | short heating steps (18 min) |
| EC constant | EC enters the rates at a fixed 7.1 M and is not depleted | the glovebox sample, where up to 0.4 M of EC is consumed (6 %) |
| Ideal solution | concentrations are used in place of activities | 1 M water in an organic solvent |
| DEC inert, no salt | DEC only dilutes; there is no Li⁺, LiPF₆ or HF | the 0.5 % water sample may contain PF₆⁻ |

- **Why EC is held constant:** it is the solvent, in 16 to 50-fold excess over anything it reacts with. Fixing it removes one equation and changes the result very little.

### 3.4 Inputs of a run: recipe and temperature history

A run needs two things besides the rate constants.

- **Recipe:** the concentrations at mixing. Everything not in the recipe starts at zero.
- **Temperature history:** a list of segments, each a temperature and a duration, from mixing to the last spectrum. The rate constants are recalculated at each segment, and the run continues from the concentrations it has reached.

| Sample | Recipe (M), all with EC 7.1 | History from mixing | The simulation is read at |
|---|---|---|---|
| 0.5 % water | TMSPA 0.149, H₂O 0.26 | **unknown age**, then 19 h at 22 °C | 1 time: the end |
| 2 % water | TMSPA 0.149, H₂O 1.05 | **unknown age**, then 145 h at 22 °C with six steps of 18 min at 30, 40, … 80 °C | 6 times: the start, after the 30 and 40 °C steps, and after each later step |
| TMSPa + TMSOH, tube A | TMSPA 0.149, TMSOH 0.18 | **unknown age** at 22 °C | 1 time: the end |
| TMSOH, heated in the probe | TMSOH 0.45 | **unknown age**, 139 h at 22 °C, 9 min at 80 °C, 5 h at 22 °C | 3 times: before, during and after the heating |
| TMSOH, glovebox | TMSOH 0.45 | 8 h at 80 °C | 1 time: the end |
| E1 (paper) | TMSOH 0.45 | one week at 25 °C | 1 time: the end (compared with a limit from the paper, not a spectrum) |

- **Why the last column:** the run gives the concentrations at every instant, but the lab only has a spectrum at a few. The simulation is read at those instants and compared there. Which spectra can be used, and why, is in Part 4.2.
- **Water against TMSPA:** the 2 % sample has 7 waters per TMSPA, more than the 3 that full hydrolysis needs. The 0.5 % sample has 1.8. It can only hydrolyse fully if the released TMSOH ends as HMDSO, which gives half of the water back (1.5 waters per TMSPA are then enough).
- **The unknown age is the first segment.** No mixing time is recorded, so the time each sample spent at room temperature before its first spectrum is not known. Only the glovebox sample has a known history. The scenarios of Part 5.4 are assumptions about this one number.

**What is assumed in these inputs**

| Input | Value used | What is actually known |
|---|---|---|
| Heating steps of the 2 % sample | 18 min each | the minimum the files allow; the paper states 8 h per step |
| Hold of the glovebox sample | 8 h at 80 °C | taken from the paper, assumed to be this sample |
| Products at mixing | zero | the TMSPa-only sample already shows 12 % BMSPA in its first spectrum, so the samples were not that clean |
| Water in the solvent | zero besides what was added | never measured; TMSPa without added water still hydrolyses |

---

## Part 4. NMR

**What this part delivers:** the two numbers that are compared for every peak: the area share the model predicts and the area share the lab measured. Full detail: theory, [Block 8 §2](theory-multiscale_microkinetics.md#2-physical-nmr-parameters).

```
C(t) from the reactor (Part 3)  ──►  predicted area shares ─────┐
                                                                ├──►  compared (Part 5)
raw lab spectrum  ─────────────►  measured area shares ± error ─┘
```

### 4.1 Peaks

- **One species, one peak.** In a spectrum of one nucleus, each species that contains that nucleus gives a peak at its own position (ppm).
- **Area share:** the area of one peak divided by the sum of the areas of all the peaks of that spectrum. It is a fraction between 0 and 1, and the shares of one spectrum add up to 1.
- **Why shares and not concentrations:** a spectrum has no absolute scale. It says how the nuclei are distributed among the species, not how many moles there are.

#### 4.1.1 Simulated: from concentrations to area shares

- **How:** each species contributes (number of observed nuclei per molecule) × (its concentration). The share is that contribution divided by the total.

    $\text{share}_k = \dfrac{n_k\,C_k}{\sum_i n_i\,C_i}$

- **³¹P:** every phosphate has one P, so the share is simply its fraction of the phosphate.

| ³¹P peak | Species |
|---|---|
| TMSPA | TMSPA |
| BMSPA | BMSPA |
| MMSPA | MMSPA |
| H₃PO₄ | H₃PO₄ |

- **¹³C, Si–CH₃ region:** every Si carries three methyl carbons, so the share is a fraction of the silicon.

| ¹³C peak | Species | Si per molecule |
|---|---|---|
| TMSOH | TMSOH | 1 |
| HMDSO | HMDSO | 2 |
| TMS-EG | TMSOEG and TMSOdiEG together | 1 each |

- **Example, the glovebox sample:** concentrations of 0.19 M TMSOH, 0.036 M HMDSO and 0.19 M ring-opened product give shares of 0.42, 0.16 and 0.42. HMDSO counts twice because it holds two Si.
- **No peak shapes or positions are needed.** The comparison uses areas only. The model can also draw a full spectrum, with positions from DFT and a peak shape ([theory, Block 8 §3](theory-multiscale_microkinetics.md#3-lorentzian-convolution)), but that is for display. The DFT positions are off by up to 2.6 ppm in ³¹P and 3–6 ppm in ²⁹Si.

**What the simulated shares cannot show**

- **Water, CO₂ and EC** give no peak in these two regions.
- **³¹P is blind to TMSOH and HMDSO.** In the phosphate samples nothing measures where the released TMS groups went.
- **The two ring-opened products share one peak.** The model puts almost all of it in the second product (TMSOdiEG), where the spectra were assigned to the first (TMS-EG). The shared peak hides this disagreement.

#### 4.1.2 Lab data: from raw spectrum to area shares

**The steps**

1. **Raw signal → spectrum.** The instrument file is Fourier-transformed and phased automatically. This was checked against the spectrometer's own processed spectra for one sample only (the 0.5 % water one).
2. **Windows.** A ppm window is set around each observed peak.
3. **Area.** The signal is summed inside each window, above a baseline.
4. **Share.** Each area is divided by the sum of the areas.

| ³¹P peak | Window (ppm) | | ¹³C peak | Window (ppm) |
|---|---|---|---|---|
| TMSPA | −27.0 to −22.0 | | TMSOH | 0.05 to 0.45 |
| BMSPA | −16.5 to −12.5 | | HMDSO | 0.6 to 1.1 |
| MMSPA | −8.5 to −5.0 | | TMS-EG | −1.95 to −1.45 |
| H₃PO₄ | −1.0 to 3.5 | | | |

- **The ppm axis is not referenced.** H₃PO₄ appears at +1.6 ppm here and at −0.3 ppm in the paper. It does not matter for shares, because the windows are placed on the peaks as they appear.

**The error of a share**

- **Noise:** measured in the parts of the spectrum with no peak, and carried into the share. It falls with the number of scans: about 0.01 with 16 scans, about 0.003 with 256.
- **Baseline:** the area depends on where the baseline is drawn. Each share is computed with three different baselines; the value kept is the middle one, and half the spread between them is the baseline error.
- **Combined:** the two are added in quadrature. A third term, a declared floor, is added later (Part 6.3).
- **Negative shares:** an empty window can give a small negative share, such as −0.013. That is the size of the error, not a signal.

**The measured shares**

| Sample | Scans | Shares (± noise and baseline) |
|---|---|---|
| 0.5 % water | 16 | TMSPA 0.93 ± 0.02, BMSPA 0.06 ± 0.01, MMSPA and H₃PO₄ zero within error |
| 2 % water, six spectra | 16, then 256 | H₃PO₄ 0.66–0.67, MMSPA 0.31–0.32, BMSPA 0.03–0.05, TMSPA zero within error; the same in all six |
| TMSPa + TMSOH, tube A | 16 | BMSPA 0.69 ± 0.06, MMSPA 0.19 ± 0.04, TMSPA 0.11 ± 0.04 |
| TMSOH, probe, three spectra | 64 | TMSOH 0.90–0.96, HMDSO 0.02–0.06, TMS-EG 0.03 at most |
| TMSOH, glovebox | 64 | TMSOH 0.42 ± 0.02, TMS-EG 0.43 ± 0.04, HMDSO 0.16 ± 0.06 |

- **Source:** the shares come from the tables exported by `03_feasible_region`, not from a new pass over the raw spectra.

### 4.2 NMR fundamentals for ³¹P and ¹³C: when an area is an amount

- **The principle:** the area of a peak is proportional to the number of nuclei behind it. That is what allows a share to be read as a fraction of molecules.
- **It only holds if the spectrum is recorded properly.** Two effects break it.

| Effect | What happens | How it is avoided |
|---|---|---|
| Incomplete relaxation | a spectrum is the sum of many scans. If the next scan starts before the nuclei have recovered, slow-recovering species give less signal than they should | wait long enough between scans (the relaxation delay) |
| NOE | the protons are irradiated to simplify the spectrum. This boosts the signal of nearby nuclei, by a different amount for each species | switch the irradiation off between scans ("NOE off") |

- **Consequence:** two spectra recorded with different settings cannot be compared, and a spectrum with the NOE on shows which species are present but not reliably how much.

**What each nucleus gives in this project**

| Nucleus | How it was recorded | Used? | Reason |
|---|---|---|---|
| ³¹P at room temperature | NOE off, 5 s delay | yes | the quantitative settings; four well-separated peaks; minutes per spectrum |
| ³¹P at 30–80 °C | NOE on, 1 s delay | no | areas distorted, not comparable with the room-temperature ones |
| ¹³C, Si–CH₃ region | NOE on, 2 s delay | yes, with a larger error | not quantitative in general, but the three peaks compared are the same kind of carbon (a methyl on Si), so they are assumed to be distorted equally |
| ²⁹Si | about 9 h per spectrum, or a faster non-quantitative method | no | shows which species are present, not how much |
| ¹H | minutes per spectrum | no | the peaks overlap and are not assigned; the solvent appears twice in some spectra |

- **³¹P is the backbone.** It sees the four phosphates directly and its shares can be trusted. Its weak point is that the recovery times were never measured, so 5 s may not be quite enough.
- **¹³C rests on an assumption.** "Distorted equally" is reasonable but not checked. This is why its shares carry a larger declared error (the larger of 0.03 and 15 %) than the ³¹P ones (the larger of 0.01 and 5 %).
- **Out of 89 spectra, 12 enter the fit.** Most were recorded to see which species are present, not to count them.

**What no spectrum measures**

- **Water.** Its amount is never read, at mixing or afterwards. This is why "equilibrium" and "the water ran out" cannot be told apart for the 2 % sample.
- **TMSOH and HMDSO in the phosphate samples.** ³¹P does not see them, and the ²⁹Si spectra are not quantitative.
- **Time.** No sample was followed while it reacted. Each share is a snapshot.

--> Part 4 hands to Part 5 a list of 44 measured shares with their errors, and the rule to compute the same 44 shares from any simulated run.
