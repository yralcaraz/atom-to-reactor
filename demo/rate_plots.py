"""Figures for Blocks 5–6: barrier models compared on the network, and k(T) over the protocol window."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

from demo.style import BLUE_RAMP, FAMILY_COLORS, GRID, INK, INK_MUTED, INK_SOFT, RULE, SERIES, model_color
from kinetics.constants import H_SI, KB_EV, KB_SI
from kinetics.microkinetics import calculate_barrier, calculate_eyring_rate
from kinetics.microkinetics.models import get_model

_MARKERS = ['o', 's', '^', 'D', 'v', 'P']


def plot_barrier_comparison(rates_by_model: dict, T_K: float, g_eV: float = 0.80, constraints=None):
    """(A) barrier shapes at a common g, with the network's ΔG_rxn as ticks;
    (B) forward barrier of every reaction under each model, with k_f(T_K) on the top axis and, when given,
    the experimental ΔG‡ windows (kinetics.get_barrier_windows) shaded per reaction.

    rates_by_model: {model name: DataFrame from ModelSpec.calculate_rates(T_K)}.
    """
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13, 5.2), layout='constrained',
                                     gridspec_kw={'width_ratios': [1, 1.15]})

    # A. Shapes
    first = next(iter(rates_by_model.values()))
    x = np.linspace(-1.0, 0.6, 400)
    shapes = [('bep_cap', 'capped BEP (legacy)', {'alpha': 0.5}), ('marcus', 'Marcus', {}),
              ('agmon_levine', 'Agmon-Levine', {}), ('blowers_masel', 'Blowers-Masel (w = 5 eV)', {'w': 5.0})]
    for k, (shape, label, extra) in enumerate(shapes):
        ax_a.plot(x, [calculate_barrier(xi, shape, g_eV, **extra)[0] for xi in x], color=SERIES[k],
                  ls=['-', '--', '-.', ':'][k], label=label)
    for dG in first['dG_rxn_eV']:
        ax_a.axvline(dG, ymin=0, ymax=0.04, color=INK_SOFT, lw=1.2)
    ax_a.axvline(0.0, color=RULE, lw=0.8)
    ax_a.set_xlabel('ΔG_rxn (eV)   (ticks: the 9 network reactions)')
    ax_a.set_ylabel(f'ΔG‡_f (eV) at g = {g_eV:.2f} eV')
    ax_a.set_title('A. The BEP cap is flat for every downhill step', loc='left')
    ax_a.legend(loc='upper left')

    # B. Network barriers per model
    rxn_ids = list(first.index)
    y = np.arange(len(rxn_ids))
    n = len(rates_by_model)
    for k, (name, rates) in enumerate(rates_by_model.items()):
        ax_b.scatter(rates['dG_barrier_f_eV'], y + (k - (n - 1) / 2) * 0.16, s=42, color=model_color(name),
                     marker=_MARKERS[k % len(_MARKERS)], label=name, zorder=3, edgecolor='white', linewidth=0.6)
    if constraints is not None:
        x_lo, x_hi = ax_b.get_xlim()
        x_lo, x_hi = min(x_lo, np.nanmin(constraints[['low_eV', 'high_eV']].values) - 0.05), max(x_hi, 1.45)
        for rid, c in constraints.iterrows():
            if rid in rxn_ids:
                lo = x_lo if np.isnan(c['low_eV']) else c['low_eV']
                hi = x_hi if np.isnan(c['high_eV']) else c['high_eV']
                ax_b.barh(rxn_ids.index(rid), hi - lo, left=lo, height=0.8, color=GRID, edgecolor=RULE, zorder=0)
        ax_b.set_xlim(x_lo, x_hi)
        ax_b.add_artist(ax_b.legend(handles=[Patch(color=GRID, label='experimental window (Gogoi 2024, derived)')],
                                    loc='lower left', fontsize=8))
    ax_b.set_yticks(y, [f"{rid} ({first.loc[rid, 'class'].replace('_', ' ')})" for rid in rxn_ids], fontsize=8.5)
    ax_b.invert_yaxis()
    ax_b.grid(axis='y', visible=False)
    ax_b.set_xlabel('ΔG‡_f (eV)')
    kT = KB_EV * T_K
    prefactor = KB_SI * T_K / H_SI
    top = ax_b.secondary_xaxis('top', functions=(lambda b: np.log10(prefactor) - b / (kT * np.log(10)),
                                                 lambda lk: (np.log10(prefactor) - lk) * kT * np.log(10)))
    top.set_xlabel(f'log₁₀ k_f at {T_K - 273.15:.0f} °C  (s⁻¹ or M⁻¹s⁻¹)', color=INK_MUTED)
    ax_b.set_title('B. Forward barriers of the network per model', loc='left')
    ax_b.legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=len(rates_by_model))
    return fig


def plot_rate_window(model, T_low_K: float = 293.15, T_high_K: float = 353.15, constraints=None, **options):
    """k_f of every reaction at the two ends of the protocol window (dumbbell, log scale).

    constraints: experimental ΔG‡ windows, drawn as k windows at the nearer end of the window
    (barriers taken as temperature independent).
    """
    spec = get_model(model)
    lo = spec.calculate_rates(T_low_K, **options)
    hi = spec.calculate_rates(T_high_K, **options)
    y = np.arange(len(lo))
    fig, ax = plt.subplots(figsize=(8, 4.2), layout='constrained')
    for yi, a, b in zip(y, lo['k_f'], hi['k_f']):
        ax.plot([a, b], [yi, yi], color=RULE, lw=1.5, zorder=1)
        ax.annotate(f'×{b / a:.0f}', (b, yi), xytext=(6, 0), textcoords='offset points', va='center',
                    fontsize=8, color=INK_SOFT)
    ax.scatter(lo['k_f'], y, color=BLUE_RAMP[0], s=50, zorder=3, label=f'{T_low_K - 273.15:.0f} °C')
    ax.scatter(hi['k_f'], y, color=BLUE_RAMP[4], s=50, zorder=3, label=f'{T_high_K - 273.15:.0f} °C')
    handles = []
    if constraints is not None:
        rxn_ids = list(lo.index)
        for rid, c in constraints.iterrows():
            if rid not in rxn_ids:
                continue
            T_plot = T_high_K if abs(c['T_K'] - T_high_K) < abs(c['T_K'] - T_low_K) else T_low_K
            k_max = calculate_eyring_rate(c['low_eV'], T_plot) if not np.isnan(c['low_eV']) else 1e6
            k_min = calculate_eyring_rate(c['high_eV'], T_plot) if not np.isnan(c['high_eV']) else 1e-14
            color = BLUE_RAMP[4] if T_plot == T_high_K else BLUE_RAMP[0]
            ax.plot([k_min, k_max], [rxn_ids.index(rid) + 0.28] * 2, color=color, lw=5, alpha=0.35,
                    solid_capstyle='butt', zorder=2)
        handles = [plt.Line2D([], [], color=BLUE_RAMP[2], lw=5, alpha=0.35,
                              label='experimental window (Gogoi 2024), at the marked temperature')]
    ax.set_xscale('log')
    finite = np.concatenate([lo['k_f'].values, hi['k_f'].values])
    ax.set_xlim(finite.min() / 30, finite.max() * 300)
    ax.set_yticks(y, [f"{rid}  ΔG‡ = {b:.2f} eV" for rid, b in zip(lo.index, lo['dG_barrier_f_eV'])], fontsize=8.5)
    ax.invert_yaxis()
    ax.grid(axis='y', visible=False)
    ax.set_xlabel('k_f  (s⁻¹ or M⁻¹s⁻¹)')
    ax.set_title(f'{spec.name}: forward rate constants across the protocol window', loc='left', color=INK)
    h, l = ax.get_legend_handles_labels()
    ax.legend(h + handles, l + [x.get_label() for x in handles], loc='upper center', bbox_to_anchor=(0.5, -0.14),
              ncol=3, fontsize=8)
    return fig


def plot_arrhenius(k_f, classes: dict, windows=None, window_C=(20.0, 80.0)):
    """log10 k_f vs 1000/T per reaction (colour: family), protocol window shaded, experimental k windows as bars.

    k_f: DataFrame (index T_K, columns reactions). windows: ΔG‡ windows (kinetics.data.get_barrier_windows),
    converted to k at their own temperature.
    """
    fig, ax = plt.subplots(figsize=(9.5, 5.2), layout='constrained')
    inv_T = 1000.0 / k_f.index.values
    groups = {}
    for rxn in k_f.columns:
        color = FAMILY_COLORS.get(classes[rxn], INK_MUTED)
        ax.plot(inv_T, np.log10(k_f[rxn]), color=color, lw=1.8)
        groups.setdefault(round(float(np.log10(k_f[rxn].iloc[-1])), 2), []).append(rxn)
    for value, members in groups.items():
        ax.annotate(', '.join(members), (inv_T[-1], value), xytext=(5, 0), textcoords='offset points', va='center',
                    fontsize=8, color=INK_SOFT)
    lo, hi = 1000.0 / (window_C[1] + 273.15), 1000.0 / (window_C[0] + 273.15)
    ax.axvspan(lo, hi, color=GRID, alpha=0.6, zorder=0)
    if windows is not None:
        for rid, w in windows.iterrows():
            if rid not in k_f.columns:
                continue
            k_hi = calculate_eyring_rate(w['low_eV'], w['T_K']) if not np.isnan(w['low_eV']) else None
            k_lo = calculate_eyring_rate(w['high_eV'], w['T_K']) if not np.isnan(w['high_eV']) else None
            x = 1000.0 / w['T_K']
            y_lo = np.log10(k_lo) if k_lo else ax.get_ylim()[0]
            y_hi = np.log10(k_hi) if k_hi else ax.get_ylim()[1]
            ax.plot([x, x], [y_lo, y_hi], color=FAMILY_COLORS.get(classes[rid], INK_MUTED), lw=6, alpha=0.35,
                    solid_capstyle='butt')
            ax.annotate(f'exp. {rid}', (x, y_hi if k_hi else y_lo), xytext=(6, 0), textcoords='offset points',
                        fontsize=8, color=INK_SOFT, va='center')
    top = ax.secondary_xaxis('top', functions=(lambda x: 1000.0 / x - 273.15, lambda c: 1000.0 / (c + 273.15)))
    top.set_xlabel('T (°C)', color=INK_MUTED)
    ax.set_xlabel('1000 / T (K⁻¹)')
    ax.set_ylabel('log₁₀ k_f  (s⁻¹ or M⁻¹s⁻¹)')
    handles = [plt.Line2D([], [], color=c, label=f.replace('_', ' ')) for f, c in FAMILY_COLORS.items()]
    handles.append(plt.Line2D([], [], color=INK_MUTED, lw=6, alpha=0.35, label='experimental window (Gogoi 2024)'))
    ax.legend(handles=handles, loc='lower left', fontsize=8)
    ax.set_title('Arrhenius plot: shaded band = the 20–80 °C protocol window', loc='left')
    return fig
