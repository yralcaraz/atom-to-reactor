"""Figures for the measured lab spectra (NB04): acquisition timeline, stacked spectra, processing check.

Source: Y. Alcaraz Galván; plots NMR data from experiments
"""

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

from demo.nmr_plots import _draw_measured
from demo.style import BLUE_RAMP, INK, INK_MUTED, INK_SOFT, RULE, SERIES

# Fixed colour and marker per nucleus
NUCLEUS_STYLE = {'1H': (SERIES[0], 'o'), '13C': (SERIES[1], 's'), '31P': (SERIES[2], 'D'),
                 '29Si': (SERIES[3], '^'), '19F': (SERIES[4], 'v')}
# Element codes used by kinetics.get_measured_shifts
_ELEMENT = {'1H': 'H', '13C': 'C', '31P': 'P', '29Si': 'Si'}


def plot_acquisition_timeline(inventory, *, heated_above_C: float = 25.0):
    """One row per sample, one marker per acquisition (colour and shape = nucleus); hollow markers were
    recorded above heated_above_C. One panel per measurement campaign (calendar year), stacked."""
    years = sorted(inventory['acquired_at'].dt.year.unique())
    fig, axes = plt.subplots(len(years), 1, figsize=(11, 0.55 * inventory['sample'].nunique() + 2.0),
                             gridspec_kw={'height_ratios': [inventory.loc[inventory['acquired_at'].dt.year == y,
                                                                          'sample'].nunique() for y in years],
                                          'hspace': 0.45}, squeeze=False)
    axes = axes.T
    for ax, year in zip(axes[0], years):
        sub = inventory[inventory['acquired_at'].dt.year == year]
        samples = list(dict.fromkeys(sub.sort_values('acquired_at')['sample']))[::-1]
        y_of = {s: i for i, s in enumerate(samples)}
        for nucleus, (color, marker) in NUCLEUS_STYLE.items():
            rows = sub[sub['nucleus'] == nucleus]
            if rows.empty:
                continue
            heated = rows['temperature_C'] > heated_above_C
            jitter = 0.12 * (list(NUCLEUS_STYLE).index(nucleus) - 2)
            y = rows['sample'].map(y_of) + jitter
            ax.scatter(rows.loc[~heated, 'acquired_at'], y[~heated], color=color, marker=marker, s=42, zorder=3,
                       label=nucleus, edgecolor='white', linewidth=0.8)
            ax.scatter(rows.loc[heated, 'acquired_at'], y[heated], facecolor='none', edgecolor=color,
                       marker=marker, s=42, linewidth=1.4, zorder=3)
        ax.set_yticks(range(len(samples)))
        ax.set_yticklabels(samples, fontsize=8.5, color=INK_SOFT)
        ax.set_ylim(-0.6, len(samples) - 0.4)
        ax.xaxis.set_major_locator(mdates.AutoDateLocator(maxticks=10))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
        ax.set_title(str(year), loc='left')
        ax.grid(axis='y', visible=False)
    handles, labels = axes[0][0].get_legend_handles_labels()
    for ax in axes[0][1:]:
        for h, l in zip(*ax.get_legend_handles_labels()):
            if l not in labels:
                handles.append(h)
                labels.append(l)
    handles.append(plt.Line2D([], [], ls='', marker='o', mfc='none', mec=INK_SOFT))
    labels.append(f'recorded above {heated_above_C:.0f} °C')
    fig.legend(handles, labels, loc='upper center', ncol=len(labels), bbox_to_anchor=(0.5, 1.02))
    return fig


def plot_spectra_stack(traces, *, window_ppm, nucleus: str, experimental=None, title: str = None,
                       ordered: bool = False, gap: float = 1.1, ax=None):
    """Spectra stacked bottom-up (first trace lowest), each scaled to its own maximum inside window_ppm.

    traces: list of (label, ppm, intensity). ordered=True colours them along a sequential ramp (a time or
    temperature series); otherwise all traces are drawn in ink. experimental: measured shifts drawn dashed.
    """
    lo, hi = sorted(window_ppm)
    if ax is None:
        _, ax = plt.subplots(figsize=(9, 0.9 + 0.75 * len(traces)))
    ramp = np.linspace(0, len(BLUE_RAMP) - 1, len(traces)).round().astype(int) if ordered else None
    for i, (label, ppm, intensity) in enumerate(traces):
        mask = (ppm >= lo) & (ppm <= hi)
        y = intensity[mask]
        scale = np.max(y) if np.max(y) > 0 else 1.0
        color = BLUE_RAMP[ramp[i]] if ordered else INK
        ax.plot(ppm[mask], y / scale + gap * i, color=color, lw=1.0)
        ax.annotate(label, xy=(1.005, gap * i + 0.15), xycoords=('axes fraction', 'data'), fontsize=8,
                    color=INK_SOFT, va='bottom', annotation_clip=False)
    if experimental is not None and nucleus in _ELEMENT:
        _draw_measured(ax, experimental, _ELEMENT[nucleus], lo, hi, vertical=True)
    ax.set_xlim(hi, lo)
    ax.set_yticks([])
    ax.grid(axis='y', visible=False)
    ax.set_xlabel(f'δ {nucleus} (ppm, as acquired)')
    if title:
        ax.set_title(title, loc='left')
    return ax


def plot_processing_check(ours, export, *, windows_ppm, nucleus: str, title: str = None):
    """Our transform of the FID against the spectrometer's own export, one panel per window.
    ours: one (ppm, intensity) per window (phase refined for it); export: (ppm, intensity).
    Each trace is scaled to its own maximum in the window."""
    fig, axes = plt.subplots(1, len(windows_ppm), figsize=(4.2 * len(windows_ppm), 3.0), squeeze=False)
    for ax, mine, (lo, hi) in zip(axes[0], ours, windows_ppm):
        m_ref = (export[0] >= lo) & (export[0] <= hi)
        m_our = (mine[0] >= lo) & (mine[0] <= hi)
        ref_max, our_max = export[1][m_ref].max(), mine[1][m_our].max()
        ax.plot(export[0][m_ref], export[1][m_ref] / ref_max, color=RULE, lw=3.0, label='spectrometer export')
        ax.plot(mine[0][m_our], mine[1][m_our] / our_max, color=SERIES[0], lw=1.0, label='this reader')
        ax.set_xlim(hi, lo)
        ax.set_yticks([])
        ax.grid(axis='y', visible=False)
        ax.set_xlabel(f'δ {nucleus} (ppm)')
    axes[0][0].legend(loc='upper left')
    if title:
        fig.suptitle(title, x=0.01, ha='left', fontsize=11, color=INK)
    fig.tight_layout()
    return fig


def plot_peak_ratio_series(table, *, x: str, ratios: dict, ylabel: str, title: str = None, note: str = None):
    """Ratios against time, one line per ratio (direct labels at the line ends).
    table: DataFrame with column x and one column per ratio; ratios: {column: label}."""
    fig, ax = plt.subplots(figsize=(8, 3.4))
    for i, (col, label) in enumerate(ratios.items()):
        ax.plot(table[x], table[col], marker='o', color=SERIES[i], lw=1.6, label=label)
        ax.annotate(label, xy=(table[x].iloc[-1], table[col].iloc[-1]), xytext=(6, 0), textcoords='offset points',
                    va='center', fontsize=8.5, color=INK_SOFT, annotation_clip=False)
    ax.legend(loc='center', ncol=len(ratios))
    ax.set_ylabel(ylabel)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
    ax.set_ylim(0, max(1.0, 1.1 * float(table[list(ratios)].max().max())))
    if title:
        ax.set_title(title, loc='left')
    if note:
        ax.annotate(note, xy=(0.99, 0.97), xycoords='axes fraction', ha='right', va='top', fontsize=7.5,
                    color=INK_MUTED)
    return fig
