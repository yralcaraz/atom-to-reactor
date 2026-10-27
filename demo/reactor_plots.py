"""Figures for Blocks 7 and 10: concentration panels of a simulated run.

Source: Y. Alcaraz Galván
"""

import matplotlib.pyplot as plt
import numpy as np

from demo.style import SERIES, add_temperature_strip, label_line_ends, mark_acquisitions, select_post_injection

CONCENTRATION_PANELS = [
    ('Phosphate ester ladder (R1–R3, R5–R7)', ['TMSPA', 'BMSPA', 'MMSPA', 'H3PO4']),
    ('Silanol / siloxane (R4–R7)', ['TMSOH', 'HMDSO']),
    ('Water (R1–R4)', ['H2O']),
    ('EC ring-opening products (R8–R9)', ['TMSOEG', 'TMSOdiEG', 'CO2']),
]


def _time_axis(sim: dict):
    """Time [h] and a mask: from TMSPA addition for protocol runs, the full grid for batch runs."""
    t = sim['t_h']
    mask = select_post_injection(sim)
    label = sim.get('time_label') or ('Time from TMSPA addition (h)' if 'injection_idx' in sim else 'Time (h)')
    return t, mask, label


def plot_concentration_panels(sim: dict, log_time: bool = False, panels=CONCENTRATION_PANELS, from_zero: bool = False):
    """2 × 2 concentration panels; protocol runs get the temperature program above each column.

    from_zero: every concentration axis starts at 0 and has no offset, so a change of a few µM on top of a
    constant level is drawn as the flat line it is, not stretched over the whole panel.
    """
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
        if from_zero:
            top = max(float(np.max(sim['C_mM'][[sim['idx'][sp] for sp in species]][:, mask])), 0.0)
            ax.set_ylim(0.0, 1.08 * top if top > 0 else 1.0)
            ax.ticklabel_format(axis='y', useOffset=False)
        mark_acquisitions(ax, sim)
        ax.set_title(title, loc='left')
        ax.set_ylabel('Concentration (mM)')
        ax.set_xlabel(xlabel)
        if len(species) > 1:
            ax.legend(loc='center left')
        label_line_ends(ax, t[mask][-1], ends)
    fig.suptitle(f"Concentrations — model: {sim['model']}", x=0.01, ha='left', fontsize=11)
    return fig

