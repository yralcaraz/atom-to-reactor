# Block-by-block audit: notebook × theory × modules

| | |
|---|---|
| **Status** | Historical. It audits the first 23-cell notebook; function and module names have since changed (see [reference-api_migration.md](reference-api_migration.md)). |
| **Scope** | Systematic cross-reference of `multiscale_microkinetics.ipynb`, `theory-multiscale_microkinetics.md` and the 8 Python modules. |
| **Markers** | ✓ verified item · ⚠ error · ⚡ false assumption · 💡 improvement |

---

## Block 0 — Physical Constants & Environment (Cells 1–2)

**Theory ↔ Code cross-check:**

| Constant | Theory | Notebook (Cell 2) | Module | Match? |
|---|---|---|---|---|
| $R$ | 8.314462618 J/(mol·K) | `R_SI = 8.314462618` | `R_SI = 8.314462618` | ✓ |
| $k_B$ | 1.380649×10⁻²³ J/K | `KB_SI = 1.380649e-23` | `KB_SI = 1.380649e-23` | ✓ |
| $h$ | 6.62607015×10⁻³⁴ J·s | `H_SI = 6.62607015e-34` | `H_SI = 6.62607015e-34` | ✓ |
| $N_A$ | 6.02214076×10²³ | `NA = 6.02214076e23` | `NA = 6.02214076e23` | ✓ |
| $P°$ | 10⁵ Pa | `P_REF = 1.0e5` | `P_REF = 1.0e5` | ✓ |
| 1 eV→kJ/mol | 96.4853 | `EV_TO_KJ_MOL = 96.485332...` | `EV_TO_KJ_MOL` computed from `eV→J × NA / 1000` | ✓ |

> [!TIP]
> **No issues found.** Constants are consistent across all three sources. CODATA 2018 values used throughout.

---

## Block 1 — Species Database (Cells 3–4)

### ✓ Verified
- 27 species loaded correctly from Peter's B3LYP-D3 dataset.
- Free energies stored in Hartree, converted to eV via `HAR2EV = 27.211386245988`.
- Solvation offsets in `dE_solv_benchmark` dict match theory Section Block 1.
- NMR shifts match theory Block 8 table (TMSPA: −19.86, BMSPA: −18.01, etc.).

### ⚠ Finding 1: **Placeholder vibrational frequencies & moments of inertia**

```python
'moments_amu_A2': (150.0, 200.0, 300.0),
'frequencies_cm1': [45.0, 85.0, 150.0, 300.0, 650.0, 1050.0, 1260.0, 2960.0],
```

**Every single species** (all 27) receives **identical** placeholder moments of inertia and the **same 8 vibrational frequencies**. This means:

- In `qRRHO` mode, all species compute identical $S_{\text{trans}}$ (mass varies), identical $S_{\text{rot}}$ (wrong! TMSPA has $I \gg$ H₂O), and identical $S_{\text{vib}}$.
- **Impact:** The `qRRHO` thermodynamic mode produces **physically meaningless** $G_{\text{sol}}(T)$ values, because all species share the same vibrational spectrum. The temperature-dependent $\Delta G_{\text{rxn}}$ from qRRHO mode is unreliable.
- **Mitigated because:** The downstream blocks use `selected_mode = 'wb97mv'`, which bypasses stat-mech entirely and returns the direct DFT Gibbs energy (`G_wb97mv_eV`) directly. So **simulation results are unaffected**, but the qRRHO comparison columns in Block 4's table are nonsense.

> [!CAUTION]
> **Severity: HIGH (for qRRHO mode) / NONE (for benchmark mode).** The qRRHO mode is fundamentally broken until real per-species Hessian data (frequencies + moments of inertia) are populated. Theory doc Section Block 2 describes the physics correctly but the code lacks the input data to execute it.

### ⚡ Finding 2: **Symmetry number heuristic is incorrect**

```python
'sigma_rot': 1 if 'TMS' in name else 2,
```

- This assigns $\sigma_{\text{rot}} = 1$ to any species containing "TMS" in the name, and $\sigma_{\text{rot}} = 2$ to everything else.
- **Errors:** H₂O should have $\sigma_{\text{rot}} = 2$ (C₂ᵥ symmetry) → gets 2 ✓. CO₂ is linear with $\sigma = 2$ → gets 2 ✓. But CH₄ has $T_d$ symmetry with $\sigma = 12$ → gets 2 ✗. EC (cyclic C₂) should be 2 → gets 2 ✓. siloxyl (C₂ᵥ) should be 2 → gets 2 ✓.
- **Impact:** Only matters in qRRHO mode (which is already broken by Finding 1). In benchmark mode, $\sigma_{\text{rot}}$ is never used.

> [!NOTE]
> **Severity: LOW.** Only affects qRRHO mode, which is already non-functional due to Finding 1.

---

## Block 2 — Statistical Mechanics Engine (Cells 5–6)

### ✓ Verified (Module: `calculate_gas_thermo.py`)
- Sackur–Tetrode translational partition: $q_{\text{trans}} = V_{\text{ideal}} / \Lambda^3$ with $V = k_BT/P°$ — **correct**.
- Rigid rotor: handles linear ($I_A = 0$, 2 DoF) vs non-linear (3 DoF) — **correct**.
- ZPE: $\sum \frac{1}{2} h c \nu_k$ — **correct**.
- Quasi-RRHO: Head-Gordon weights $w = 1/(1 + (\nu_0/\nu)^4)$ with $\nu_0 = 100$ cm⁻¹ — **matches Grimme 2012, correct**.
- Free rotor entropy: $S_{\text{free rot}} = R[0.5 + \ln\sqrt{8\pi^3 I_{\text{eff}} k_B T / h^2}]$ with $I_{\text{eff}} = \hbar / (4\pi c \nu)$ — **correct**.
- Benchmark mode: returns `G_B3_eV` directly, skipping stat-mech — **correct by design**.

### ⚡ Finding 3: **Enthalpy uses $\frac{5}{2}RT$ not $\frac{5}{2}RT + PV$**

```python
H_trans = 2.5 * R_SI * T_K  # This is 5/2 RT
```

For an ideal gas, $H = U + PV = \frac{3}{2}RT + RT = \frac{5}{2}RT$. The $PV = RT$ is already embedded in the $\frac{5}{2}$ factor. **This is correct.**

### ✓ Verified: Theory ↔ Code equations match
- All partition function formulas in theory Block 2 are faithfully implemented.
- The TMSPA temperature sweep plot (Cell 6) is cosmetic only — it demonstrates the engine works, but uses placeholder data (Finding 1).

---

## Block 3 — Solvation Thermochemistry (Cells 7–8)

### ✓ Verified (Modules: `calculate_standard_state_shift.py`, `calculate_solution_gibbs.py`)

- Standard-state shift: $\Delta G^{\circ \to *} = RT \ln(C^*/C°_{\text{gas}})$ where $C° = P°/(RT)$ — **correct**.
- Solution Gibbs: $G_{\text{sol}} = G°_{\text{gas}} + \Delta E_{\text{solv}} + \Delta G^{\circ \to *}$ — **matches theory exactly**.
- Benchmark mode: `dG_std_shift_kJ_mol = 0.0` — **correct** (matches `tank_model.ipynb` which omits this shift).

### ⚡ Finding 4: **Solution entropy calculation is questionable**

In `calculate_solution_gibbs.py`:
```python
'S_sol_J_mol_K': gas_thermo['S_gas_J_mol_K'] - (dG_std_shift_kJ_mol * 1000.0 / T_K)
```

This computes $S_{\text{sol}} = S_{\text{gas}} - \Delta G^{\circ \to *} / T$. But the standard-state shift $\Delta G^{\circ \to *} = RT\ln(C^*/C°)$ is purely entropic (no enthalpy component), so this is **thermodynamically consistent**: $\Delta S^{\circ \to *} = -\Delta G^{\circ \to *}/T = -R\ln(C^*/C°)$.

However, it **ignores** the solvation entropy contribution: $\Delta E_{\text{solv}}$ is treated as a pure energy offset with **zero entropy component** ($\Delta S_{\text{solv}} = 0$). In reality, solvation has significant entropic costs (ordering of solvent molecules). 

> [!NOTE]
> **Severity: MEDIUM (for qRRHO mode).** In benchmark mode, `S_sol` is returned as 0.0 and is never used downstream. The Van 't Hoff slopes (Block 6) derive from $\Delta G(T)$ differences, which in benchmark mode are temperature-independent, so slopes reflect only $\partial(k_BT/h)/\partial T$ and barrier terms.

### ✓ Verified: Solvation dumbbell plot (Cell 8)
- Correctly computes gas-phase vs. solution-phase $\Delta G_{\text{rxn}}$ for all 9 reactions.
- The 6-out-of-9 sign-flip claim matches theory Block 3. ✓

---

### ⚠ Finding 19: **Solvation uncertainty was propagated incorrectly, and its definition is unknown**

The snapshot's `uncertainty_kjmol` is the quadrature sum of three MD stds (gas solute, solute-in-EC box, pure-EC box). The pure-EC box is shared by all species, so its term cancels in the 2 → 2 steps. Treating it as independent doubles σ(ΔG_rxn) (0.47–0.60 → 0.23–0.44 eV). Whether the stds are per-frame deviations or standard errors of the mean is not documented.

> [!WARNING]
> **Severity: HIGH (if per-frame).** **Status: open.** With per-frame σ, the signs of R2–R4 and R6–R8 are not resolved at ±1σ. Block 3 now shows provisional error bars (`kinetics/thermo/solvation_uncertainty.py`). The question is tracked in the TFM open-questions document (MET).

---

## Block 4 — Reaction Network & Wegscheider (Cells 9–10)

### ✓ Verified (Module: `calculate_reaction_thermo.py`)
- All 9 reactions correctly wired with stoichiometric coefficients matching theory Block 4.
- $\Delta G_{\text{rxn}} = \sum_p G_{\text{sol},p} - \sum_r G_{\text{sol},r}$ — **correct**.
- $K_{\text{eq}} = \exp(-\Delta G / RT)$ — **correct**.
- Wegscheider loop: $\Delta G(R_1) + \Delta G(R_4) - \Delta G(R_5) = 0$ verified to machine precision for both modes — **correct**.

### ⚡ Finding 5: **Only 1 out of 3 possible Wegscheider cycles is tested**

The theory identifies the cycle $R_1 + R_4 = R_5$. But the reaction network contains **3 independent Hess cycles**:
1. $R_1 + R_4 - R_5 = 0$ ✓ (tested)
2. $R_2 + R_4 - R_6 = 0$ ✗ (not tested)
3. $R_3 + R_4 - R_7 = 0$ ✗ (not tested)

These additional cycles are guaranteed to close because the species energies are self-consistent, but **it would strengthen the audit** to verify all three explicitly.

> [!NOTE]
> **Severity: LOW.** The untested cycles are mathematically guaranteed to close (species-level additivity). But adding them would catch any future data-entry typo.

### ⚡ Finding 6: **`qRRHO` comparison in the table is misleading**

The Block 4 table shows "benchmark vs. qRRHO" reaction energies side-by-side. But per Finding 1, the qRRHO values use placeholder frequencies, so the $\Delta(\Delta G)$ column and qRRHO $K_{\text{eq}}$ values are physically meaningless. The theory document says the two modes should give "different but self-consistent" results — this is true mathematically, but the qRRHO numbers are garbage input.

> [!IMPORTANT]
> **Action:** Add a warning note in the notebook markdown stating that qRRHO columns are computed with placeholder vibrational data and should not be interpreted physically until real Hessian data is loaded.

---

## Block 5 — Microkinetic Rate Constants (Cells 11–12)

### ✓ Verified (Module: `calculate_rate_constants.py`)
- BEP: $\Delta G^\ddagger_f = \max(E_0, E_0 + \alpha \Delta G_{\text{rxn}})$ — **matches theory Block 5 exactly**.
- Marcus: $\Delta G^\ddagger_f = (\lambda/4)(1 + \Delta G/\lambda)^2$ with $\lambda = 4E_0$ — **matches theory exactly**.
- Eyring: $k_f = (k_BT/h)\exp(-\Delta G^\ddagger / RT)$ — **correct**.
- Detailed balance: $k_r = k_f / K_{\text{eq}}$ — **correct**.
- Wegscheider kinetic ratio verified for both BEP and Marcus — **correct**.

### ⚠ Finding 7: **BEP uses `max()` on kJ/mol, not eV — consistent but fragile**

```python
dG_barrier_f_kJ_mol = max(E0_kJ_mol, E0_kJ_mol + alpha * dG_rxn_kJ_mol)
```

The `max()` operates on kJ/mol (after converting $E_0$ from eV). This is **numerically correct** because $\max$ is unit-invariant for monotonic conversions. No bug here.

### ✓ Finding 8 (RESOLVED): **Family-Specific BEP Activation Barriers ($E_{0,\text{solv}} = 1.30$ eV)**

- **Previous Flaw:** `family_bep_parameters` previously set a uniform $E_0 = 0.80$ eV across all reaction classes. At room temperature, this caused cyclic carbonate ring-opening ($R_8, R_9$) to proceed $\sim 10^8\times$ too fast, generating artificial $\text{CO}_2$ bubbling and prematurely consuming TMSOH before silyl transfer ($R_5$–$R_7$) could occur.
- **Resolution:** Following the experimental observations of Gogoi et al. (*J. Phys. Chem. C* 2024, 128, 1654; the earlier citation to *Nat. Commun.* 2024 was incorrect), `solvent_attack` is assigned $E_0 = 1.30$ eV in both `calculate_rate_constants.py` and `multiscale_microkinetics.ipynb` (Block 5). The paper reports no barrier: 1.30 eV is **derived** from its protocol (TMSOH opens EC only after 8 h holds at 80 °C → observed barrier 1.25–1.36 eV).
- **Physical Result in Block 7:**
  - Forward rate constant drops from $0.186\text{ s}^{-1}$ to $6.59 \times 10^{-10}\text{ s}^{-1}$.
  - False $\text{CO}_2$ gassing is completely suppressed ($0.00\text{ mM}$ at $10^6\text{ s}$).
  - TMSOH is preserved for productive silyl transfer: TMSPA is fully consumed down to $10.0\text{ mM}$ with $20.0\text{ mM}$ siloxyl formed, resolving Peter Broqvist's original benchmark anomaly.
- **Documented:** Added Section 5.3 and Section 4 notes to `theory-multiscale_microkinetics.md`.

### ⚠ Finding 15: **Capped BEP is direction-dependent (correctness bug)**

`max(E0, E0 + αΔG)` is not invariant under reaction reversal. Writing R4 as `HMDSO + H2O → 2 TMSOH` changes its physical rate by ×7.2 and the tank trajectories by up to 2.5 mM. The same cap puts every exergonic step exactly at E0, so α acts on 1 of 9 reactions.

> [!WARNING]
> **Severity: HIGH (correctness).** **Status: mitigated (opt-in).** The Level 1 engine (`kinetic_model='level1'`, `kinetics/microkinetics/barrier_models.py`) uses smooth, reversal-invariant forms. This is verified by `tests/test_barrier_models.py::test_network_reversal_invariance`. The legacy default is kept for Peter-reference reproduction.

### ⚠ Finding 16: **Condensation barrier contradicts Gogoi et al. 2024**

5 vol% TMSOH in EC/DEC forms no TMSOTMS beyond impurity level even after 80 °C (Gogoi et al., *J. Phys. Chem. C* 2024). With E0 = 0.80 eV the model equilibrates R4 within about 45 s (about 22 % TMSOH conversion). The observation implies an uncatalysed ΔG‡(R4) ≳ 1.30 eV. The "R4 autocatalytic water regeneration" narrative is not supported in acid-free media.

> [!CAUTION]
> **Severity: HIGH.** `LEVEL1_FAMILY_PARAMETERS` sets g(condensation) = 1.30 eV. A catalysed R4 channel (acids from R1–R3) must be added explicitly (Level 3) if needed.

### ⚡ Finding 17: **Block 4 ΔG° table in the theory doc was stale**

The theory doc listed B3LYP-era values (e.g. R4 +0.183, R6 +0.029, R8 −0.449 eV) while the code uses the ωB97M-V snapshot (R4 +0.102, R6 −0.056, R8 −0.047 eV). **Status: resolved** (table regenerated 2026-09-28).

### ⚡ Finding 18: **Level 1 hydrolysis/transfer barriers are placeholders**

Under Marcus, the 0.80 eV placeholders make R1 ~5 × 10³ faster than under the capped BEP. TMSPA is then consumed within the first 20 °C hold of the NB02 protocol, which is faster than observed. **Severity: MEDIUM. Status: open.** Recalibrate g(hydrolysis) and g(transfer) against time-resolved TMSPA decay before making `level1` the default.

---

## Block 6 — Arrhenius Fitting (Cells 14–15)

### ✓ Verified (Module: `fit_modified_arrhenius.py`)
- OLS regression of $\ln k = \ln A + \beta \ln(T/T_0) - E_a/(RT)$ — **correct**.
- $R^2$ computation — **correct**.
- Temperature grid 273.15–373.15 K with 25 points — adequate for a smooth fit.

### ⚡ Finding 9: **Arrhenius fit is exact by construction ($R^2 = 1.0000$)**

Because the underlying model is Eyring TST with a temperature-independent barrier:
$$k(T) = \frac{k_BT}{h} \exp\left(-\frac{\Delta G^\ddagger}{RT}\right)$$

This is **already** in modified Arrhenius form with $A = k_B T_0/h$, $\beta = 1.0$, $E_a = \Delta G^\ddagger$. The 3-parameter regression recovers these values exactly. The fit adds no new information — it's a tautology.

**When would this NOT be a tautology?** Only when $\Delta G^\ddagger(T)$ varies with temperature (e.g., in qRRHO mode with real Hessian data, where entropy contributions create genuine non-Arrhenius curvature).

> [!NOTE]
> **Severity: LOW.** Not a bug — the regression is correct. But the "perfect fit" should not be presented as validation; it's a consequence of the benchmark mode's temperature-independent barriers.

**Status: documented (2026-09-28).** Block 6 in the notebook and the theory doc (§ Block 6.4) now state that the fit is a consistency check, not validation. The Block 6 analysis text was also regenerated: it had stale B3LYP-era numbers (e.g. R4 barrier 0.89 eV, ln K(R4) ≈ −7).

---

## Block 7 — Batch Reactor Simulation (Cells 17–18)

### ✓ Verified (Module: `simulate_tank_reactor.py`)
- Stoichiometric matrix construction — **correct**.
- ODE: $dC/dt = S \cdot r(C,T)$ with mass-action kinetics — **correct**.
- EC buffered mode: EC concentration held constant at 4.5 M — **correct** (documented as approximation in theory).
- `np.maximum(C, 0.0)` clamp prevents negative concentrations — **standard practice**.
- Radau IIA with `rtol=1e-8, atol=1e-12` — **adequate for stiff problems**.
- Si and P conservation checks — **correct** (hardcoded atom counts match stoichiometry).

### ⚠ Finding 10: **Integration starts at $t = 1$ s, not $t = 0$ s**

```python
t_eval = np.logspace(0, np.log10(t_end_s), n_points)
```

`np.logspace(0, ...)` starts at $10^0 = 1$ s. But `solve_ivp` is called with `t_span = [t_eval[0], t_end_s]` = `[1.0, 1e6]`. This means:
- The first second of reaction dynamics (where initial fast transients occur) is **completely skipped**.
- Initial conditions at $t = 0$ are applied at $t = 1$ s, and the integrator starts from there.

**Impact:** For the benchmark case with $k_f \sim 0.2$ s⁻¹, the first second sees $\sim 20\%$ conversion of water. This initial transient is missed in the output (though the integrator itself handles it internally between $t=1$ and the first output point).

> [!WARNING]
> **Severity: MEDIUM.** The physics during $0 < t < 1$ s is lost from the output. Fix: use `np.logspace(-3, np.log10(t_end_s), n_points)` to capture dynamics from $t = 10^{-3}$ s.

### ⚡ Finding 11: **`mode='wb97mv'` is hardcoded in `simulate_tank_reactor`**

```python
def simulate_tank_reactor(..., mode: str = 'wb97mv', ...):
```

The default `mode` in the simulator is `'wb97mv'`, while all upstream Block 4/5 cells use `selected_mode` (also set to `'wb97mv'`). This is consistent **now**, but if a user changes `selected_mode = 'qRRHO'` in Block 4 without also passing `mode='qRRHO'` to the tank simulator, the results would silently use different thermodynamics.

> [!NOTE]
> **Severity: LOW.** Currently consistent, but the default should ideally inherit from the global `selected_mode` variable.

---

## Block 8 — Virtual NMR Spectrometer (Cells 19–20)

### ✓ Verified (Module: `simulate_virtual_nmr.py`)
- Lorentzian lineshape: $L(\delta) = \gamma^2 / ((\delta - \delta_0)^2 + \gamma^2)$ with $\gamma = \text{FWHM}/2$ — **correct**.
- Intensity weighting: $I = n_{\text{Si}} \times C_i(t)$ — **correct** (matches theory Block 8).
- Chemical shifts match theory and Cell 4's `nmr_29si_shifts`.

### ⚡ Finding 12: **Stacked spectra use hardcoded FWHM = 0.8 ppm, but the 1D plot uses 0.4**

In Cell 20, the `simulate_virtual_nmr` function is called with `lw=0.8`, which generates the 2D heatmap correctly. But the 1D stacked spectra loop uses:

```python
spec += conc_mM * n_si * ( (0.4)**2 / ((delta - shift)**2 + (0.4)**2) )
```

This hardcodes $\gamma = 0.4$ ppm (i.e., $\text{FWHM} = 0.8$ ppm). Since $\gamma = \text{FWHM}/2 = 0.4$, this is **actually correct and consistent**. No bug.

### ⚡ Finding 13: **BMSPA and TMSOEG share the same chemical shift (−18.01 ppm)**

**Status: obsolete (2026-09-28).** These were B3LYP-era values. With the ωB97M-V snapshot, BMSPA is at 28.5 ppm and TMSOEG at 21.6 ppm, so they no longer overlap. The xMSPA triplet (TMSPA 27.6, BMSPA 28.5, MMSPA 31.0 ppm) is the new near-degeneracy, which Gogoi et al. 2024 also could not resolve in ²⁹Si. See the theory doc, Block 8.2.

Both species appear at $\delta = -18.01$ ppm in the NMR data. In a real experiment, these peaks would overlap and be indistinguishable. The simulation correctly sums their intensities at the same position, but the theory document does not discuss this ambiguity.

> [!NOTE]
> **Severity: LOW (physical limitation, not a code error).** This is an inherent spectroscopic degeneracy. The theory should mention that the −18.01 ppm peak is composite (BMSPA + TMSOEG) and cannot be deconvolved without ²D NMR or HSQC.

---

## Block 9 — Sensitivity Sweep (Cells 21–22)

### ✓ Verified
- E₀ sweep over [0.65, 0.70, 0.75, 0.80, 0.85] eV — **correct** range.
- Each sweep creates independent `family_bep_parameters` with uniform E₀ — **correct**.
- Half-life read off at 25 mM (50% of 50 mM initial) — **correct**.
- CSV export of final timeseries — **useful**.

### ⚡ Finding 14: **Sweep range doesn't include the corrected solvent_attack barrier**

The sweep varies $E_0$ uniformly across all families. To test the effect of a higher `solvent_attack` barrier ($E_0 \geq 1.30$ eV), the sweep should allow family-specific $E_0$ values, not just a global scalar.

> [!NOTE]
> **Severity: LOW.** This is a feature request, not a bug. Block 9 demonstrates the concept; family-specific sweeps can be added later.

---

## Summary of Findings

| # | Block | Severity | Type | Description |
|---|---|---|---|---|
| 1 | 1 | **HIGH** | ⚡ False assumption | All 27 species share identical placeholder frequencies & moments of inertia → qRRHO mode is broken |
| 2 | 1 | LOW | ⚡ False assumption | Symmetry numbers assigned by string heuristic (`'TMS' in name`), wrong for CH₄, etc. |
| 3 | 2 | — | ✓ Verified | Enthalpy $\frac{5}{2}RT$ is correct (PV already included) |
| 4 | 3 | MEDIUM | ⚡ Assumption | Solvation entropy set to zero ($\Delta S_{\text{solv}} = 0$); only affects qRRHO mode |
| 5 | 4 | LOW | ⚡ Incomplete | Only 1/3 Wegscheider cycles explicitly tested |
| 6 | 4 | LOW | ⚡ Misleading | qRRHO comparison columns use garbage placeholder data |
| 7 | 5 | — | ✓ Verified | `max()` in kJ/mol is unit-invariant, no bug |
| 8 | 5 | **CRITICAL** | ✓ Resolved | Family-specific $E_0 = 1.30$ eV for `solvent_attack` (Gogoi 2024); false gassing suppressed |
| 9 | 6 | LOW | ✓ Documented | Arrhenius $R^2 = 1.0$ is exact by construction in benchmark mode; stated as a consistency check |
| 10 | 7 | **MEDIUM** | ✓ Resolved | Integration starts at $t_0 = 0$, $t_{\text{eval}}$ from $10^{-3}$ s to capture initial fast transients |
| 11 | 7 | LOW | ✓ Resolved | `mode` inherits dynamically from upstream `selected_thermo_mode` |
| 12 | 8 | — | ✓ Verified | FWHM is consistent between 2D and 1D plots |
| 13 | 8 | LOW | Obsolete | BMSPA/TMSOEG overlap was B3LYP-era; now xMSPA triplet overlaps in ²⁹Si (resolved in ³¹P) |
| 14 | 9 | LOW | Feature | Sweep uses global $E_0$; cannot test family-specific barriers |
| 15 | 5 | **HIGH** | ⚠ Correctness | Capped BEP is direction-dependent; mitigated by opt-in Level 1 engine |
| 16 | 5 | **HIGH** | ⚠ Contradiction | Condensation E0 = 0.80 eV contradicts Gogoi 2024; Level 1 uses g = 1.30 eV |
| 17 | 4 | LOW | ✓ Resolved | Stale B3LYP ΔG° table in theory doc regenerated from ωB97M-V snapshot |
| 18 | 5 | MEDIUM | Open | Level 1 hydrolysis/transfer g are placeholders; too fast under Marcus |
| 19 | 3 | **HIGH** | Open | Solvation σ: shared EC term was double-counted; std vs. SEM undocumented; provisional error bars in Block 3 |

### Critical Path Status

Both findings on the original Critical Path that materially distorted simulation results have been **RESOLVED**:

1. ✓ **Finding 8 (CRITICAL - RESOLVED):** `solvent_attack` $E_0 = 1.30$ eV (Gogoi 2024) calibrated in `calculate_rate_constants.py` and Block 5. Forward rate drops by $3 \times 10^8$, stopping false gassing ($0.0\text{ mM}$ $\text{CO}_2$) and preserving TMSOH for stoichiometric silyl transfer.
2. ✓ **Finding 10 (MEDIUM - RESOLVED):** `t_span = [0, t_end_s]` with `t_start_s = 1e-3` in `simulate_tank_reactor.py`, capturing the sub-second fast dynamics of water scavenging.

Remaining open items are low-priority refinements (Finding 5: 3-cycle Wegscheider verification, Finding 14: family sweeps, and notes on placeholder qRRHO frequencies).
