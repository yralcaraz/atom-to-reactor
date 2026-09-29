"""Figures for Blocks 8, 11–12: stacked spectra per acquisition, peak tracking and the water readout."""

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

from demo.style import (
    GRID, INK_MUTED, INK_SOFT, RULE, SEQUENTIAL_CMAP, SERIES, add_temperature_strip, label_line_ends, mark_acquisitions,
    select_post_injection,
)
from kinetics.spectroscopy import (
    DEFAULT_NUCLEI, EXCHANGE_PEAK, LINE_SHAPES, build_nmr_catalog, calculate_nmr_peaks, find_shift_windows,
    simulate_nmr_spectra,
)


def _acquisition_matrix(sim: dict) -> np.ndarray:
    return np.array([a['C_M'] for a in sim['acquisitions']]).T


def _group_labels(sites: list, min_sep: float, max_names: int = 2) -> list:
    """Merge (name, ppm) labels closer than min_sep ppm: [(ppm, 'A · B')]; larger groups read 'A +n'."""
    groups = []
    for name, x in sorted(sites, key=lambda item: item[1], reverse=True):
        if groups and abs(groups[-1][0][-1] - x) < min_sep:
            groups[-1][0].append(x)
            if name not in groups[-1][1]:
                groups[-1][1].append(name)
        else:
            groups.append(([x], [name]))
    return [(float(np.mean(xs)), ' · '.join(names) if len(names) <= max_names else f'{names[0]} +{len(names) - 1}')
            for xs, names in groups]


def plot_stacked_spectra(sim: dict, element: str, nmr_catalog: dict = None, experimental=None):
    """One trace per NMR acquisition (lowest = first), on broken shift axes that skip empty ranges.

    experimental: measured shifts (kinetics.get_measured_shifts) drawn as dashed lines.
    """
    catalog = nmr_catalog if nmr_catalog is not None else build_nmr_catalog()
    shape = LINE_SHAPES[element]
    reference, isotope = DEFAULT_NUCLEI[element]
    C = _acquisition_matrix(sim)
    peaks = calculate_nmr_peaks(element, C, sim['idx'], catalog)
    shifts = [v for _, s, _ in peaks for v in (s.min(), s.max())]
    windows = find_shift_windows(shifts, shape['pad_ppm'], shape['max_gap_ppm'])
    widths = [hi - lo for hi, lo in windows]
    grids = [np.linspace(hi, lo, 2500) for hi, lo in windows]
    spectra = [simulate_nmr_spectra(element, C, sim['idx'], x, nmr_catalog=catalog) for x in grids]
    step = 0.35 * max(S.max() for S in spectra)
    colors = plt.get_cmap(SEQUENTIAL_CMAP)(np.linspace(0.4, 0.95, C.shape[1]))

    fig, axes = plt.subplots(1, len(windows), figsize=(10, 5.6), sharey=True, squeeze=False, layout='constrained',
                             gridspec_kw={'width_ratios': [max(w, 0.3 * sum(widths)) for w in widths]})
    fig.get_layout_engine().set(wspace=0.0, w_pad=0.02)
    axes = axes[0]
    top = step * (C.shape[1] - 1) + max(S.max() for S in spectra)
    label_pos = [(name, float(s.mean())) for name, s, _ in peaks]
    for k, (ax, (hi, lo), x, S) in enumerate(zip(axes, windows, grids, spectra)):
        for i in range(C.shape[1]):
            ax.plot(x, S[i] + i * step, color=colors[i], lw=1.2)
        ax.set_xlim(hi, lo)
        ax.grid(axis='y', visible=False)
        in_window = [(n, p) for n, p in label_pos if lo <= p <= hi]
        for ppm, text in _group_labels(in_window, 0.04 * sum(widths)):
            ax.annotate(text, (ppm, top), xytext=(0, 4), textcoords='offset points', rotation=90, ha='center',
                        va='bottom', fontsize=7.5, color=INK_SOFT)
        if experimental is not None:
            _draw_measured(ax, experimental, element, lo, hi, vertical=True)
        if k > 0:
            ax.spines['left'].set_visible(False)
            ax.tick_params(left=False)
        if k < len(axes) - 1:
            ax.spines['right'].set_visible(False)
    if experimental is not None and len(experimental[experimental['nucleus'] == element]):
        axes[0].legend(handles=[plt.Line2D([], [], color=INK_MUTED, ls='--', lw=1, label='measured (Gogoi 2024)')],
                       loc='upper left', fontsize=8)
    for i, a in enumerate(sim['acquisitions']):
        axes[-1].annotate(f"{a['T_C']:.0f} °C", xy=(1.01, i * step), xycoords=('axes fraction', 'data'),
                          fontsize=8, color=INK_SOFT)
    axes[0].set_ylabel('Intensity (mM of nuclei, offset)')
    axes[0].set_ylim(-0.05 * top, top * 1.25)
    fig.supxlabel(f'δ {isotope} (ppm vs {reference})', fontsize=10, color=INK_SOFT)
    fig.suptitle(f"Predicted {isotope} spectra after each hold ({sim['model']})", x=0.01, ha='left', fontsize=11)
    return fig


def plot_peak_tracking(sim: dict, element: str, nmr_catalog: dict = None, log_time: bool = False,
                       experimental=None):
    """Shift-time map of the densest window (square-root intensity scale) and every species' integral.

    experimental: measured shifts (kinetics.get_measured_shifts) drawn as dashed lines on the map.
    """
    catalog = nmr_catalog if nmr_catalog is not None else build_nmr_catalog()
    shape = LINE_SHAPES[element]
    reference, isotope = DEFAULT_NUCLEI[element]
    protocol = 'injection_idx' in sim
    mask = select_post_injection(sim)
    t, C = sim['t_h'][mask], sim['C_M'][:, mask]
    peaks = calculate_nmr_peaks(element, C, sim['idx'], catalog)
    resolved = [(n, s, i) for n, s, i in peaks if n != EXCHANGE_PEAK]
    windows = find_shift_windows([s[0] for _, s, _ in resolved], shape['pad_ppm'], shape['max_gap_ppm'])
    hi, lo = max(windows, key=lambda w: sum(lo_ <= s[0] <= hi_ for _, s, _ in resolved for hi_, lo_ in [w]))
    measured = experimental[experimental['nucleus'] == element] if experimental is not None else None
    if measured is not None and len(measured):
        near = measured[(measured['mid_ppm'] < hi + shape['max_gap_ppm']) & (measured['mid_ppm'] > lo - shape['max_gap_ppm'])]
        if len(near):
            hi = max(hi, near['high_ppm'].max() + 0.5 * shape['pad_ppm'])
            lo = min(lo, near['low_ppm'].min() - 0.5 * shape['pad_ppm'])
    x = np.linspace(hi, lo, 800)
    S = simulate_nmr_spectra(element, C, sim['idx'], x, fwhm_ppm=max(shape['fwhm_ppm'], 0.02 * (hi - lo)),
                             nmr_catalog=catalog)

    varying_T = np.ptp(sim['T_K'][mask]) > 0
    ratios = [0.45, 2.0, 1.6] if varying_T else [2.0, 1.6]
    fig, axes = plt.subplots(len(ratios), 1, figsize=(10.5, 7.8 if varying_T else 6.8), sharex=True,
                             layout='constrained', gridspec_kw={'height_ratios': ratios})
    if varying_T:
        add_temperature_strip(axes[0], t, sim['T_K'][mask], [a['t_h'] for a in sim['acquisitions']])
    ax_map, ax_int = axes[-2], axes[-1]
    ax_map.pcolormesh(t, x, S.T, cmap=SEQUENTIAL_CMAP, shading='auto', norm=mcolors.PowerNorm(0.5), rasterized=True)
    ax_map.set_ylim(lo, hi)
    ax_map.set_ylabel(f'δ {isotope} (ppm)')
    ax_map.grid(False)
    if measured is not None and len(measured):
        _draw_measured(ax_map, experimental, element, lo, hi, vertical=False)
    in_window = [(name, s[0]) for name, s, _ in resolved if lo <= s[0] <= hi]
    for ppm, text in _group_labels(in_window, 0.03 * (hi - lo)):
        ax_map.annotate(text, xy=(1.005, ppm), xycoords=('axes fraction', 'data'), fontsize=7.5,
                        color=INK_SOFT, va='center')
    integrals = {}
    for name, _, inten in resolved:
        integrals[name] = integrals.get(name, 0.0) + inten
    ends = {}
    for k, (name, inten) in enumerate(integrals.items()):
        ax_int.plot(t, inten, color=SERIES[k % len(SERIES)], label=name, lw=1.8)
        ends[name] = inten[-1]
        if protocol:
            k_acq = [a['idx'] - sim['injection_idx'] for a in sim['acquisitions']]
            ax_int.plot(t[k_acq], inten[k_acq], 'o', ms=4.5, color=SERIES[k % len(SERIES)],
                        markeredgecolor='white', zorder=4)
    if protocol:
        mark_acquisitions(ax_int, sim)
        mark_acquisitions(ax_map, sim)
    if log_time:
        ax_int.set_xscale('log')
    ax_int.set_ylabel(f'Integral (mM {element})')
    ax_int.set_xlabel('Time from TMSPA addition (h)' if protocol else 'Time (h)')
    ax_int.legend(loc='center left', ncol=2 if len(integrals) > 4 else 1)
    label_line_ends(ax_int, t[-1], ends)
    fig.suptitle(f"{isotope} peaks through the run ({sim['model']})", x=0.01, ha='left', fontsize=11)
    return fig


def _draw_measured(ax, experimental, element: str, lo: float, hi: float, vertical: bool) -> None:
    """Dashed lines at the measured shifts inside [lo, hi], labelled 'exp. <species>'."""
    rows = experimental[(experimental['nucleus'] == element) & (experimental['mid_ppm'] >= lo)
                        & (experimental['mid_ppm'] <= hi)]
    for ppm in sorted(set(rows['mid_ppm'])):
        (ax.axvline if vertical else ax.axhline)(ppm, color=INK_MUTED, ls='--', lw=1.0, zorder=2)
    sites = list(zip(rows['species'], rows['mid_ppm']))
    for ppm, text in _group_labels(sites, 0.03 * (hi - lo), max_names=3):
        if vertical:
            ax.annotate(f'exp. {text}', xy=(ppm, 0.02), xycoords=('data', 'axes fraction'), rotation=90,
                        ha='right', va='bottom', fontsize=7, color=INK_MUTED)
        else:
            ax.annotate(f'exp. {text}', xy=(0.005, ppm), xycoords=('axes fraction', 'data'), va='bottom',
                        fontsize=7.5, color=INK_SOFT)


def plot_water_balance(sim: dict, water: dict):
    """(A) who holds the exchangeable OH protons; (B) water read by mass balance at each acquisition."""
    mask = select_post_injection(sim)
    t = sim['t_h'][mask]
    carriers = water['oh_carriers_mM'][mask]
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(12.5, 4.4), layout='constrained')
    water_first = sorted(carriers.columns, key=lambda c: not c.startswith('H2O'))
    colors = [GRID] + SERIES[:len(water_first) - 1]
    ax_a.stackplot(t, carriers[water_first].values.T, labels=water_first, colors=colors,
                   edgecolor='white', linewidth=0.6)
    ax_a.set_title('OH protons only move between carriers; the total stays 2·[H2O]₀', loc='left', fontsize=10)
    ax_a.set_ylabel('OH protons (mM)')
    ax_a.set_xlabel('Time from TMSPA addition (h)')
    ax_a.legend(loc='lower left', ncol=2)
    ax_a.set_xlim(t[0], t[-1])

    h2o_mM = water['h2o_mass_balance_M'][mask] * 1000.0
    ax_b.plot(t, h2o_mM, color=SERIES[0], label='[H2O] from mass balance')
    idx = [a['idx'] for a in sim['acquisitions']]
    ax_b.plot(sim['t_h'][idx], water['h2o_mass_balance_M'][idx] * 1000.0, 'o', ms=6, color=SERIES[0],
              markeredgecolor='white', label='read at the acquisitions')
    mark_acquisitions(ax_b, sim)
    ax_b.set_title('[H2O] = [H2O]₀ − ½ Σ OH on the other species', loc='left', fontsize=10)
    ax_b.set_ylabel('Water (mM)')
    ax_b.set_xlabel('Time from TMSPA addition (h)')
    ax_b.legend(loc='lower left')
    return fig


def plot_shift_parity(nmr_catalog: dict, experimental):
    """Computed vs measured shifts for ²⁹Si and ³¹P (bars: reported range), with the y = x line."""
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4), layout='constrained')
    for ax, (nucleus, isotope, color) in zip(axes, [('Si', '29Si', SERIES[0]), ('P', '31P', SERIES[1])]):
        exp = experimental[experimental['nucleus'] == nucleus]
        xs, ys, offsets = [], [], []
        for _, r in exp.iterrows():
            sites = [d for d, _, labile in nmr_catalog.get(nucleus, {}).get(r['species'], []) if not labile]
            if not sites:
                continue
            ax.errorbar(r['mid_ppm'], sites[0], xerr=[[r['mid_ppm'] - r['low_ppm']], [r['high_ppm'] - r['mid_ppm']]],
                        fmt='o', ms=7, color=color, ecolor=color, capsize=3, zorder=3)
            ax.annotate(r['species'], (r['mid_ppm'], sites[0]), xytext=(6, -3), textcoords='offset points',
                        fontsize=8.5, color=INK_SOFT)
            xs.append(r['mid_ppm'])
            ys.append(sites[0])
            offsets.append(sites[0] - r['mid_ppm'])
        lo, hi = min(xs + ys) - 3, max(xs + ys) + 3
        ax.plot([lo, hi], [lo, hi], color=RULE, lw=1.2, zorder=1)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xlabel(f'δ {isotope} measured (ppm)')
        ax.set_ylabel(f'δ {isotope} computed (ppm)')
        ax.set_title(f'{isotope}: computed − measured = {np.mean(offsets):+.1f} ppm on average', loc='left', fontsize=10)
        ax.annotate('y = x', (hi, hi), xytext=(-30, -12), textcoords='offset points', fontsize=8, color=INK_MUTED)
    fig.suptitle('Computed (gas-phase DFT) vs measured shifts in EC/DEC (Gogoi et al. 2024)', x=0.01, ha='left', fontsize=11)
    return fig
