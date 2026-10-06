"""Figures for Blocks 2–4: species free energies, driving forces, Wegscheider cycles and f(T).

Source: Y. Alcaraz Galván
"""

import matplotlib.pyplot as plt
import numpy as np

from demo.style import BLUE_RAMP, FAMILY_COLORS, GRID, INK, INK_MUTED, INK_SOFT, RULE, SERIES
from kinetics.constants import EV_TO_KCAL_MOL


def plot_solvation_dumbbell(thermo, unit: str = 'kcal/mol'):
    """Gas-phase vs solution ΔG_rxn per reaction, with the ±1 standard error of the solvation term.

    thermo: DataFrame from kinetics.calculate_network_thermo.
    """
    fig, ax = plt.subplots(figsize=(8.5, 4.8), layout='constrained')
    _draw_solvation_dumbbell(ax, thermo, thermo['sigma_solv_eV'], unit, 'solution in EC, ±1 standard error')
    ax.set_title('Driving force ΔG_rxn: gas phase vs EC solution', loc='left')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=2)
    return fig


def _draw_solvation_dumbbell(ax, thermo, err_eV, unit: str, sol_label: str):
    """Gas → solution segment per reaction, solution point with ±err_eV error bars."""
    scale = EV_TO_KCAL_MOL if unit == 'kcal/mol' else 1.0
    gas = thermo['dG_gas_eV'] * scale
    sol = thermo['dG_rxn_eV'] * scale
    y = np.arange(len(thermo))
    for yi, g, s in zip(y, gas, sol):
        ax.plot([g, s], [yi, yi], color=RULE, lw=1.5, zorder=1)
    ax.scatter(gas, y, s=60, color=SERIES[0], zorder=3, label='gas phase (ωB97M-V)')
    ax.errorbar(sol, y, xerr=err_eV * scale, fmt='o', ms=8, color=SERIES[1], ecolor=SERIES[1], elinewidth=1.2,
                capsize=3, zorder=3, label=sol_label)
    for yi, g, s in zip(y, gas, sol):
        ax.annotate(f'{g:+.1f}', (g, yi), xytext=(0, 8), textcoords='offset points', ha='center', fontsize=8, color=INK_SOFT)
        ax.annotate(f'{s:+.1f}', (s, yi), xytext=(0, -13), textcoords='offset points', ha='center', fontsize=8, color=INK_SOFT)
    ax.axvline(0.0, color=INK, lw=0.8)
    ax.set_yticks(y, [f"{rid}  {eq}" for rid, eq in zip(thermo.index, thermo['equation'])], fontsize=8.5)
    if not ax.yaxis_inverted():
        ax.invert_yaxis()
    ax.set_xlabel(f'ΔG_rxn ({unit})')
    ax.grid(axis='y', visible=False)


def plot_species_temperature(entropy_J_mol_K: dict, T_K, shift_kJ_mol):
    """(A) RRHO entropy of each species, the slope of its G(T); (B) the standard-state shift ΔG°→*(T)."""
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(12.5, 4.2), layout='constrained', gridspec_kw={'width_ratios': [1.3, 1]})
    names = list(entropy_J_mol_K)
    ax_a.barh(names, [entropy_J_mol_K[n] for n in names], color=SERIES[0], height=0.65)
    for y, n in enumerate(names):
        ax_a.annotate(f'{entropy_J_mol_K[n]:.0f}', (entropy_J_mol_K[n], y), xytext=(4, 0), textcoords='offset points',
                      va='center', fontsize=8, color=INK_SOFT)
    ax_a.invert_yaxis()
    ax_a.grid(axis='y', visible=False)
    ax_a.set_xlabel('S°_gas at 298 K (J/mol/K)   —   dG/dT = −S')
    ax_a.set_title('A. Gas-phase entropy per species (RRHO, from the stored H and G)', loc='left', fontsize=10)

    T_C = np.asarray(T_K) - 273.15
    ax_b.plot(T_C, shift_kJ_mol, color=SERIES[1])
    for T_mark in (20.0, 80.0):
        y = np.interp(T_mark, T_C, shift_kJ_mol)
        ax_b.plot(T_mark, y, 'o', color=SERIES[1], ms=7)
        ax_b.annotate(f'{y:.2f} kJ/mol at {T_mark:.0f} °C', (T_mark, y), xytext=(8, -10), textcoords='offset points',
                      fontsize=8.5, color=INK_SOFT)
    ax_b.set_xlabel('T (°C)')
    ax_b.set_ylabel('ΔG°→* (kJ/mol)')
    ax_b.set_title('B. Standard-state shift, 1 bar gas → 1 M solution', loc='left', fontsize=10)
    return fig


def plot_wegscheider_cycles(dG_eV: dict, cycles=(('R1', 'R4', 'R5'), ('R2', 'R4', 'R6'), ('R3', 'R4', 'R7')),
                            equations: dict = None):
    """Energy-level diagram per cycle: two-step route (first + second) against the direct step."""
    fig, axes = plt.subplots(1, len(cycles), figsize=(4.6 * len(cycles), 4.3), sharey=True, layout='constrained')
    for ax, (first, second, direct) in zip(np.atleast_1d(axes), cycles):
        levels = [0.0, dG_eV[first], dG_eV[first] + dG_eV[second]]
        for x, y, c in zip(range(3), levels, [INK_SOFT, SERIES[0], SERIES[2]]):
            ax.hlines(y, x - 0.28, x + 0.28, color=c, lw=3)
        ax.plot([0.28, 0.72], levels[:2], ls='--', color=RULE, lw=1.2)
        ax.plot([1.28, 1.72], levels[1:], ls='--', color=RULE, lw=1.2)
        ax.annotate(f'{first}  {dG_eV[first]:+.3f} eV', (0.5, sum(levels[:2]) / 2), xytext=(-6, 0),
                    textcoords='offset points', ha='right', va='center', fontsize=8.5, color=SERIES[0])
        ax.annotate(f'{second}  {dG_eV[second]:+.3f} eV', (1.5, sum(levels[1:]) / 2), xytext=(6, 0),
                    textcoords='offset points', ha='left', va='center', fontsize=8.5, color=SERIES[0])
        ax.annotate('', xy=(1.72, dG_eV[direct]), xytext=(0.28, 0.0),
                    arrowprops=dict(arrowstyle='->', color=SERIES[2], lw=1.8))
        ax.annotate(f'{direct} direct  {dG_eV[direct]:+.3f} eV', (2.0, dG_eV[direct]), xytext=(0, 8),
                    textcoords='offset points', ha='center', fontsize=8.5, color=SERIES[2])
        residual = levels[2] - dG_eV[direct]
        ax.set_title(f'{first} + {second} − {direct}: residual {residual:+.1e} eV', loc='left', fontsize=10)
        ax.set_xticks([0, 1, 2], ['start', f'after {first}', 'end'])
        ax.set_xlim(-0.9, 2.5)
        ax.margins(y=0.15)
        ax.grid(axis='x', visible=False)
    np.atleast_1d(axes)[0].set_ylabel('Free energy relative to start (eV)')
    fig.suptitle('Wegscheider cycles: two-step and direct routes reach the same state', x=0.01, ha='left', fontsize=11)
    return fig


def plot_temperature_dependence(dG_low_eV, dG_high_eV, sigma_eV, lnK, classes: dict, T_low_C: float = 20.0,
                                T_high_C: float = 80.0):
    """(A) ΔG_rxn at the two ends of the window (RRHO entropy kept), with the solvation standard error for scale;
    (B) van 't Hoff plot ln K vs 1000/T for the stored (T-independent) ΔG_rxn.

    dG_low_eV, dG_high_eV, sigma_eV: Series indexed by reaction; lnK: DataFrame (index T_K, columns reactions).
    """
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13, 4.8), layout='constrained')
    y = np.arange(len(dG_low_eV))
    ax_a.barh(y, 2 * sigma_eV.values, left=(dG_low_eV - sigma_eV).values, height=0.6, color=GRID, zorder=0,
              label='±1 SE solvation')
    for yi, a, b in zip(y, dG_low_eV, dG_high_eV):
        ax_a.plot([a, b], [yi, yi], color=RULE, lw=1.5, zorder=1)
    ax_a.scatter(dG_low_eV, y, color=BLUE_RAMP[0], s=50, zorder=3, label=f'{T_low_C:.0f} °C')
    ax_a.scatter(dG_high_eV, y, color=BLUE_RAMP[4], s=50, zorder=3, label=f'{T_high_C:.0f} °C (ΔH, ΔS fixed at 298 K)')
    ax_a.axvline(0, color=INK, lw=0.8)
    ax_a.set_yticks(y, list(dG_low_eV.index))
    ax_a.invert_yaxis()
    ax_a.grid(axis='y', visible=False)
    ax_a.set_xlabel('ΔG_rxn (eV)')
    ax_a.set_title('A. Driving force at 20 and 80 °C', loc='left', fontsize=10)
    ax_a.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=8)

    inv_T = 1000.0 / lnK.index.values
    ends = {}
    for rxn in lnK.columns:
        ax_b.plot(inv_T, lnK[rxn], color=FAMILY_COLORS.get(classes[rxn], INK_MUTED), lw=1.6)
        ends.setdefault(round(lnK[rxn].iloc[0], 1), []).append(rxn)
    for value, members in ends.items():
        ax_b.annotate(', '.join(members), (inv_T[0], value), xytext=(4, 0), textcoords='offset points', va='center',
                      fontsize=8, color=INK_SOFT)
    ax_b.axhline(0, color=INK, lw=0.8)
    ax_b.set_xlabel('1000 / T (K⁻¹)')
    ax_b.set_ylabel('ln K_eq')
    ax_b.set_title("B. van 't Hoff: ln K = −ΔG/RT (stored ΔG, slope −ΔG/R)", loc='left', fontsize=10)
    handles = [plt.Line2D([], [], color=c, label=f.replace('_', ' ')) for f, c in FAMILY_COLORS.items()]
    ax_b.legend(handles=handles, loc='lower left', fontsize=8)
    return fig
