"""Figures for Block 13: nucleus selectivity and recovery of reaction extents.

Source: Y. Alcaraz Galván; adapted from P. Broqvist
"""

import matplotlib.pyplot as plt
import numpy as np

from demo.style import INK, SEQUENTIAL_CMAP, SERIES, SURFACE


def plot_selectivity(selectivity):
    """Heatmap of the NAS selectivity (share of each fingerprint no other reaction can mimic)."""
    fig, ax = plt.subplots(figsize=(7.5, 3.8), layout='constrained')
    values = selectivity.values
    im = ax.imshow(np.ma.masked_invalid(values), cmap=SEQUENTIAL_CMAP, vmin=0, vmax=1, aspect='auto')
    ax.set_xticks(range(values.shape[1]), selectivity.columns)
    ax.set_yticks(range(values.shape[0]), selectivity.index)
    ax.grid(False)
    for (i, j), v in np.ndenumerate(values):
        ax.text(j, i, '—' if np.isnan(v) else f'{v:.2f}', ha='center', va='center', fontsize=9,
                color=SURFACE if not np.isnan(v) and v > 0.6 else INK)
    fig.colorbar(im, ax=ax, label='selectivity (0 = shared, 1 = exclusive)')
    ax.set_title('Which nucleus isolates which identifiable reaction', loc='left')
    return fig


def plot_extent_recovery(step_T_C, xi_true, xi_recovered, labels, noise_rel: float):
    """Per-step extents Δξ of the identifiable reactions: model (bars) vs recovered from noisy spectra (dots)."""
    n = len(labels)
    cols = 3
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(12.5, 3.3 * rows), sharex=True, squeeze=False,
                             layout='constrained')
    for j, ax in enumerate(axes.flat):
        if j >= n:
            ax.set_visible(False)
            continue
        ax.bar(step_T_C, xi_true[:, j], width=6, color='#e1e0d9', edgecolor='#c3c2b7', label='model (true)')
        ax.plot(step_T_C, xi_recovered[:, j], 'o', ms=6, color=SERIES[0], markeredgecolor='white',
                label='recovered from spectra')
        ax.set_title(labels[j], loc='left', fontsize=10)
        ax.grid(axis='x', visible=False)
        if j % cols == 0:
            ax.set_ylabel('Δξ per step (mM)')
        if j >= n - cols:
            ax.set_xlabel('Hold temperature (°C)')
    axes.flat[0].legend(loc='upper left')
    fig.suptitle(f'Per-step extents recovered from the acquisitions ({noise_rel:.0%} noise per nucleus)',
                 x=0.01, ha='left', fontsize=11)
    return fig
