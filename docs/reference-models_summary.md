# Fit models: equations, variables and parameters

| | |
|---|---|
| **Status** | Draft |
| **Source** | Y. Alcaraz Galván |
| **Scope** | The five model structures the fit compares, the equations they share and what is free in each. Read from [candidates.py](../kinetics/fitting/candidates.py). |

## 1. Equations shared by all models

Network (R1–R9), four families:

| Family | Reactions |
|---|---|
| hydrolysis | R1: TMSPA + H2O ⇌ BMSPA + TMSOH<br>R2: BMSPA + H2O ⇌ MMSPA + TMSOH<br>R3: MMSPA + H2O ⇌ H3PO4 + TMSOH |
| condensation | R4: 2 TMSOH ⇌ HMDSO + H2O |
| transfer | R5: TMSPA + TMSOH ⇌ BMSPA + HMDSO<br>R6: BMSPA + TMSOH ⇌ MMSPA + HMDSO<br>R7: MMSPA + TMSOH ⇌ H3PO4 + HMDSO |
| solvent attack | R8: EC + TMSOH ⇌ TMSOEG + CO2<br>R9: EC + TMSOEG ⇌ TMSOdiEG + CO2 |

Reactor (mass action, piecewise isothermal):

```
dC/dt = S · (r_f − r_r)
r_f,j = k_f,j · Π C_i^ν_ij (reactants)        r_r,j = k_r,j · Π C_i^ν_ij (products)
```

Rate constants (Eyring + detailed balance):

```
k_f = (kB·T/h) · exp(−ΔG‡ / RT)
k_r = k_f / K_eq,        K_eq = exp(−ΔG_rxn / RT)
```

Barrier from the reaction energy and the intrinsic barrier `g` of the family:

```
Marcus       ΔG‡ = g · (1 + ΔG_rxn / 4g)²
capped BEP   ΔG‡ = max(g, g + α·ΔG_rxn),   α = 0.5      (M0-BEP only)
```

Objective (what the fit minimises):

```
χ² = Σ z²,   z = (predicted share − measured share) / σ
```

## 2. The five models

Each model adds one thing to the one before it.

| # | Model | What it assumes | Free parameters | Count |
|---|---|---|---|---|
| 1 | **M0** | One barrier for all nine reactions. Two variants: `M0-BEP` (capped BEP) and `M0-Marcus` | `g all reactions` | 1 |
| 2 | **M1** | One Marcus barrier per family | `g hydrolysis`, `g transfer`, `g condensation`, `g solvent_attack` | 4 |
| 3 | **M1-split** | M1 with R1 separate from R2 and R3 | `g hydrolysis_R1`, `g hydrolysis_R23`, `g transfer`, `g condensation`, `g solvent_attack` | 5 |
| 4 | **M3** | M1 with ΔG_rxn of R1–R4 free | the 4 barriers of M1 + `dG R1`, `dG R2`, `dG R3`, `dG R4` | 8 |
| 5 | **M3-split** | M1-split with ΔG_rxn of R1–R4 free | the 5 barriers of M1-split + `dG R1` … `dG R4` | 9 |

Freed energies (M3, M3-split):

```
dG Rj = ΔG_rxn,j − computed value        (j = R1…R4)
χ² = Σ z² + dᵀ · Σ_MD⁻¹ · d              (d = the four dG; Σ_MD = covariance of the MD standard errors)
```

- R5–R7 follow from them (transfer = hydrolysis + condensation); R8 and R9 keep their computed values.
- The constraint adds 4 residuals, one per freed energy.

Optional parameters, on top of any model:

| Parameter | When | Meaning |
|---|---|---|
| `water fraction` | structure name ends in `+W` | Fraction of the added water that is available |
| `c0 H2O TMSPa alone` | the tube "TMSPa alone" is in the sample set (third fit) | Water at mixing of the tube mixed without water |
| `log10 age 2 % H2O` | `free` scenario | Age of the 2 % H2O sample at its first spectrum |

The last two come from the sample set, not from the structure: every model takes them.

Rules of the third fit ([plan](plan-fit_improvement.md)):

- A model fits when no miss exceeds 3 standard errors.
- Trace rule: its parameters, unchanged, must still fit with 1 µM and with 1 mM TMSOH at mixing in the tubes that hold TMSPA and no TMSOH.
- The selected model is the smallest one that fits and passes the trace rule, with free ages before equal ages.

## 3. Variables and parameters

| Kind | Symbol | Unit | Where it comes from |
|---|---|---|---|
| State | C_i(t): TMSPA, BMSPA, MMSPA, H3PO4, H2O, TMSOH, HMDSO, EC, TMSOEG, TMSOdiEG, CO2 | M | Solved by the reactor |
| Input | Composition at mixing, temperature history, spectrum times | M, K, s | Lab records (`lab_observables.json`). The water of "TMSPa alone" is not recorded: fitted |
| Input | Sample age (mixing → first spectrum) | h | Not recorded: 1 h, 1 day or 7 days by scenario, or free |
| Input | ΔG_rxn of R1–R9 | eV | Computed (DFT + solvation) |
| Fixed | ΔS‡ = 0, α = 0.5 | | Assumed, not fitted |
| Fitted | `g <family>` | eV | Search box 0.60–1.70 (narrower for some families) |
| Fitted | `dG R1…R4` | eV | Search box ±0.75 around the computed value |
| Fitted | `water fraction` | – | Search box 0.05–1 |
| Fitted | `c0 H2O TMSPa alone` | M | Search box 0–0.3 |
| Fitted | `log10 age` | log10 h | 1 h to 30 days |
| Observed | ³¹P shares: TMSPA, BMSPA, MMSPA, H3PO4 | – | Lab NMR |
| Observed | ¹³C shares: TMSOH, HMDSO, TMSOEG (+ TMSOdiEG), P-silyl | – | Lab NMR |
