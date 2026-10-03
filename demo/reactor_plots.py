"""Figures for Blocks 7 and 10: concentration panels, model overlays and depletion per protocol step.

Classification: PUBLIC (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván
"""

import matplotlib.pyplot as plt
import numpy as np

from demo.style import (
    GRID, INK, INK_MUTED, INK_SOFT, RULE, SERIES, add_temperature_strip, label_line_ends, mark_acquisitions, model_color,
    model_linestyle, select_post_injection,
)
from kinetics.reactor import calculate_remaining_fraction

CONCENTRATION_PANELS = [
    ('Phosphate ester ladder (R1–R3, R5–R7)', ['TMSPA', 'BMSPA', 'MMSPA', 'H3PO4']),
    ('Silanol / siloxane (R4–R7)', ['TMSOH', 'HMDSO']),
    ('Water (R1–R4)', ['H2O']),
    ('EC ring-opening products (R8–R9)', ['TMSOEG', 'TMSOdiEG', 'CO2']),
]
_GREYS = ['#d4d3cd', '#bdbcb5', '#a6a59e', '#8f8e87', '#787771', '#61605b', '#4a4945']


def _time_axis(sim: dict):
    """Time [h] and a mask: from TMSPA addition for protocol runs, the full grid for batch runs."""
    t = sim['t_h']
    mask = select_post_injection(sim)
    label = 'Time from TMSPA addition (h)' if 'injection_idx' in sim else 'Time (h)'
    return t, mask, label


def plot_concentration_panels(sim: dict, log_time: bool = False, panels=CONCENTRATION_PANELS):
    """2 × 2 concentration panels; protocol runs get the temperature program above each column."""
    t, mask, xlabel = _time_axis(sim)
    varying_T = np.ptp(sim['T_K'][mask]) > 0
    fig = plt.figure(figsize=(12, 8.4 if varying_T else 7.4), layout='constrained')
    grid = fig.add_gridspec(3 if varying_T else 2, 2, height_ratios=[0.4, 2, 2] if varying_T else [1, 1])
    offset = 1 if varying_T else 0
    axes = [fig.add_subplot(grid[offset + k // 2, k % 2]) for k in range(4)]
    if varying_T:
        for col in range(2):
            strip = fig.add_subplot(grid[0, col], sharex=axes[col])
            add_temperature_strip(strip, t[mask], sim['T_K'][mask], [a['t_h'] for a in sim.get('acquisitions', [])])
            strip.tick_params(labelbottom=False)
    for ax, (title, species) in zip(axes, panels):
        ends = {}
        for k, sp in enumerate(species):
            y = sim['C_mM'][sim['idx'][sp]][mask]
            ax.plot(t[mask], y, color=SERIES[k], label=sp)
            ends[sp] = y[-1]
        if log_time:
            ax.set_xscale('log')
        mark_acquisitions(ax, sim)
        ax.set_title(title, loc='left')
        ax.set_ylabel('Concentration (mM)')
        ax.set_xlabel(xlabel)
        if len(species) > 1:
            ax.legend(loc='center left')
        label_line_ends(ax, t[mask][-1], ends)
    fig.suptitle(f"Concentrations — model: {sim['model']}", x=0.01, ha='left', fontsize=11)
    return fig


def plot_model_overlay(sims: dict, species: list, log_time: bool = False):
    """One panel per species, one line per model (fixed colour and style per model)."""
    n = len(species)
    fig, axes = plt.subplots(1, n, figsize=(3.6 * n, 3.7), sharex=True, layout='constrained')
    for ax, sp in zip(np.atleast_1d(axes), species):
        for name, sim in sims.items():
            t, mask, xlabel = _time_axis(sim)
            ax.plot(t[mask], sim['C_mM'][sim['idx'][sp]][mask], color=model_color(name),
                    ls=model_linestyle(name), label=name)
        if log_time:
            ax.set_xscale('log')
        ax.set_title(sp, loc='left')
        ax.set_xlabel(xlabel)
    np.atleast_1d(axes)[0].set_ylabel('Concentration (mM)')
    handles, labels = np.atleast_1d(axes)[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=len(sims))
    return fig


def plot_step_depletion(sims: dict, species: str = 'TMSPA', highlight: str = None, title: str = None, ax=None):
    """Fraction of `species` left at each acquisition vs the preceding hold temperature, one line per run.

    Registered models keep their colours; other runs (e.g. a barrier sweep) are drawn on a grey ramp,
    with `highlight` in blue.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(7.5, 4.3))
    for k, (name, sim) in enumerate(sims.items()):
        frac = calculate_remaining_fraction(sim, species)
        color = model_color(name)
        lw, ms = 2.0, 5
        if name == highlight:
            color, lw, ms = SERIES[0], 2.8, 7
        elif color == INK_MUTED:
            color = _GREYS[min(k, len(_GREYS) - 1)]
        ax.plot(frac.index, frac.values, marker='o', ms=ms, lw=lw, color=color, label=name,
                ls=model_linestyle(name), zorder=3 if name == highlight else 2)
    ax.set_ylim(-0.03, 1.05)
    ax.set_xlabel('Hold temperature before the acquisition (°C)')
    ax.set_ylabel(f'{species} remaining (fraction)')
    ax.set_title(title or f'{species} left at each acquisition', loc='left')
    ax.legend(loc='lower left', fontsize=8)
    return ax


def plot_protocol_calibration(sweep: dict, runs: dict, highlight: str = None):
    """(A) a barrier sweep of one model and (B) the registered models, as TMSPA left per acquisition."""
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13, 4.4), sharey=True, layout='constrained')
    plot_step_depletion(sweep, highlight=highlight, ax=ax_a, title='A. Barrier sweep: which value spreads depletion over the steps?')
    plot_step_depletion(runs, ax=ax_b, title='B. Registered models')
    ax_b.set_ylabel('')
    return fig


_SHORT = {
    'ring_opened_fraction': 'TMSOH ring-opened (fraction)',
    'hmdso_si_fraction': 'Si ending in HMDSO (fraction)',
    'tmspa_conversion': 'TMSPA converted (fraction)',
}


def plot_control_experiments(checks):
    """One panel per control experiment: each model's prediction against the observed window (shaded).

    checks: DataFrame from kinetics.evaluate_control_experiments.
    """
    experiments = list(dict.fromkeys(checks['experiment']))
    fig, axes = plt.subplots(1, len(experiments), figsize=(3.4 * len(experiments), 4.0), layout='constrained')
    for ax, exp_id in zip(np.atleast_1d(axes), experiments):
        rows = checks[checks['experiment'] == exp_id].reset_index(drop=True)
        low, high = rows.loc[0, 'low'], rows.loc[0, 'high']
        ax.set_yscale('log')
        ax.set_ylim(1e-4, 6.0)
        band = (1e-4 if np.isnan(low) else low, 6.0 if np.isnan(high) else high)
        ax.axhspan(*band, color=GRID, zorder=0)
        for k, r in rows.iterrows():
            value = max(r['predicted'], 1.1e-4)
            ax.scatter(k, value, s=70, color=model_color(r['model']), zorder=3, edgecolor='white', linewidth=0.8)
            ax.annotate('✓' if r['consistent'] else '✗', (k, value), xytext=(0, 9), textcoords='offset points',
                        ha='center', fontsize=11, color=INK)
        ax.set_xticks(range(len(rows)), rows['model'], rotation=35, ha='right', fontsize=8.5)
        ax.set_xlim(-0.6, len(rows) - 0.4)
        ax.grid(axis='x', visible=False)
        ax.set_title(f"{exp_id}: {rows.loc[0, 'label']}\nobserved window: {rows.loc[0, 'status']}", loc='left',
                     fontsize=9.5)
        ax.set_ylabel(_SHORT.get(rows.loc[0, 'observable'], rows.loc[0, 'observable_label']), fontsize=9)
    fig.suptitle('Models against the control experiments of Gogoi et al. 2024 (shaded: observed window)',
                 x=0.01, ha='left', fontsize=11)
    return fig


def plot_barrier_sensitivity(sweep: dict, times_h: dict, species: str = 'H2O', highlight: float = None,
                             marks: dict = None, T_K: float = 298.15):
    """(A) species(t) for each barrier of a sweep; (B) characteristic time vs barrier with the Eyring slope.

    sweep: {E0 [eV]: batch run}; times_h: {E0: characteristic time [h]}; marks: {E0: label} to annotate.
    """
    from kinetics.constants import KB_EV
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13, 4.6), layout='constrained')
    E0s = sorted(sweep)
    for k, E0 in enumerate(E0s):
        sim = sweep[E0]
        color = SERIES[0] if E0 == highlight else _GREYS[min(k, len(_GREYS) - 1)]
        ax_a.plot(sim['t_h'], sim['C_mM'][sim['idx'][species]], color=color, lw=2.6 if E0 == highlight else 1.6,
                  label=f'E0 = {E0:.2f} eV', zorder=3 if E0 == highlight else 2)
    ax_a.set_xscale('log')
    ax_a.set_xlabel('Time (h)')
    ax_a.set_ylabel(f'{species} (mM)')
    ax_a.set_title(f'A. {species} consumption for one global E0 (capped BEP)', loc='left', fontsize=10)
    ax_a.legend(loc='lower left', fontsize=7.5, ncol=2)

    x = np.array(E0s)
    y = np.array([times_h[E0] for E0 in E0s])
    ok = np.isfinite(y)
    ax_b.plot(x[ok], y[ok], 'o', color=SERIES[0], ms=7, zorder=3, label='simulated')
    if highlight in times_h and np.isfinite(times_h[highlight]):
        ref = times_h[highlight]
        grid = np.linspace(x.min(), x.max(), 50)
        ax_b.plot(grid, ref * np.exp((grid - highlight) / (KB_EV * T_K)), color=RULE, lw=1.4, zorder=1,
                  label=f'Eyring slope: ×{np.exp(0.1 / (KB_EV * T_K)):.0f} per 0.1 eV')
    for E0, text in (marks or {}).items():
        if E0 in times_h and np.isfinite(times_h[E0]):
            ax_b.annotate(text, (E0, times_h[E0]), xytext=(8, -12), textcoords='offset points', fontsize=8.5,
                          color=INK_SOFT)
    for hours, name in ((1.0, '1 h'), (24.0, '1 day'), (8766.0, '1 year')):
        ax_b.axhline(hours, color=GRID, lw=1.0, zorder=0)
        ax_b.annotate(name, (x.min(), hours), xytext=(0, 3), textcoords='offset points', fontsize=7.5, color=INK_MUTED)
    ax_b.set_yscale('log')
    ax_b.set_xlabel('E0 (eV)')
    ax_b.set_ylabel(f'time to consume 90 % of {species} (h)')
    ax_b.set_title('B. Time scale vs barrier', loc='left', fontsize=10)
    ax_b.legend(loc='upper left', fontsize=8)
    return fig
