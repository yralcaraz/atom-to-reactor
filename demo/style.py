"""Shared visual language: surface and ink tokens, categorical slots, model colors and axis helpers.

Categorical hues are assigned in a fixed order (validated for colour-vision deficiency) and follow the
entity: each registered model keeps the same colour and line style in every figure.
"""

import matplotlib.pyplot as plt
import numpy as np

from kinetics.microkinetics.models import MODELS

SURFACE = '#fcfcfb'
INK = '#0b0b0b'
INK_SOFT = '#52514e'
INK_MUTED = '#898781'
GRID = '#e1e0d9'
RULE = '#c3c2b7'

# Categorical slots in validated order: blue, orange, aqua, yellow, magenta, green, violet, red
SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']
# Sequential blue ramp (steps 250 → 700), for ordered traces such as successive acquisitions
BLUE_RAMP = ['#86b6ef', '#5598e7', '#2a78d6', '#1c5cab', '#104281', '#0d366b']
SEQUENTIAL_CMAP = 'Blues'

_LINESTYLES = ['-', '--', '-.', ':']

# Fixed colour per reaction family
FAMILY_COLORS = {'hydrolysis': SERIES[0], 'condensation': SERIES[1], 'transfer': SERIES[2], 'solvent_attack': SERIES[4]}


def apply_style() -> None:
    """Notebook-wide matplotlib defaults: light surface, recessive grid and axes, thin lines."""
    plt.rcParams.update({
        'figure.dpi': 110,
        'savefig.dpi': 150,
        'figure.facecolor': SURFACE,
        'axes.facecolor': SURFACE,
        'axes.edgecolor': RULE,
        'axes.labelcolor': INK_SOFT,
        'axes.titlecolor': INK,
        'axes.titlesize': 11,
        'axes.labelsize': 10,
        'axes.grid': True,
        'axes.axisbelow': True,
        'axes.prop_cycle': plt.cycler(color=SERIES),
        'grid.color': GRID,
        'grid.linewidth': 0.7,
        'xtick.color': INK_MUTED,
        'ytick.color': INK_MUTED,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 8.5,
        'legend.frameon': False,
        'lines.linewidth': 2.0,
        'lines.markersize': 6,
        'font.family': 'sans-serif',
        'font.sans-serif': ['DejaVu Sans'],
    })


def model_color(name: str) -> str:
    """Fixed colour of a registered model (registry order); unregistered variants are neutral grey."""
    names = list(MODELS)
    return SERIES[names.index(name) % len(SERIES)] if name in names else INK_MUTED


def model_linestyle(name: str) -> str:
    names = list(MODELS)
    return _LINESTYLES[names.index(name) % len(_LINESTYLES)] if name in names else '-'


def label_line_ends(ax, x_end: float, labels_y: dict, min_gap_frac: float = 0.05) -> None:
    """Direct labels at the right end of lines, nudged apart vertically so they never overlap."""
    if not labels_y:
        return
    lo, hi = ax.get_ylim()
    gap = min_gap_frac * (hi - lo)
    items = sorted(labels_y.items(), key=lambda kv: kv[1])
    placed = []
    for label, y in items:
        y_pos = max(y, placed[-1] + gap) if placed else y
        placed.append(y_pos)
        ax.annotate(label, xy=(x_end, y_pos), xytext=(4, 0), textcoords='offset points', va='center',
                    fontsize=8.5, color=INK_SOFT, annotation_clip=False)


def add_temperature_strip(ax, t_h: np.ndarray, T_K: np.ndarray, marks_h=()) -> None:
    """Step plot of the temperature program [°C] sharing the time axis of the panels below."""
    ax.step(t_h, T_K - 273.15, where='post', color=INK_SOFT, lw=1.6)
    for t in marks_h:
        ax.axvline(t, color=RULE, lw=0.8, ls=':')
    ax.set_ylabel('T (°C)')
    ax.grid(axis='x', visible=False)


def mark_acquisitions(ax, sim: dict) -> None:
    """Dotted verticals at the NMR acquisitions of a protocol run."""
    for a in sim.get('acquisitions', []):
        ax.axvline(a['t_h'], color=RULE, lw=0.8, ls=':', zorder=0)


def select_post_injection(sim: dict) -> np.ndarray:
    """Boolean mask of the time points from TMSPA addition on (all points for a batch run)."""
    mask = np.ones(len(sim['t_h']), dtype=bool)
    if 'injection_idx' in sim:
        mask[:sim['injection_idx']] = False
    return mask
