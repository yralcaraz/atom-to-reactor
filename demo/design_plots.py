"""Figures for Block 14: feasibility of the family barriers, the water series, timing maps and experiment design."""

import matplotlib.pyplot as plt
import numpy as np

from demo.style import (
    BLUE_RAMP, FAMILY_COLORS, GRID, INK, INK_MUTED, INK_SOFT, RULE, SEQUENTIAL_CMAP, SERIES, SURFACE,
    add_temperature_strip,
)
from kinetics.fitting import calculate_readouts

_FAMILY_TITLES = {'hydrolysis': 'hydrolysis (R1–R3)', 'transfer': 'transfer (R5–R7)',
                  'condensation': 'condensation (R4)', 'solvent_attack': 'solvent attack (R8–R9)'}
_OBSERVABLE_AXIS = {'ring_opened_fraction': 'TMSOH ring-opened', 'hmdso_si_fraction': 'Si in HMDSO',
                    'tmspa_conversion': 'TMSPA remaining'}


def _format_hours(h: float) -> str:
    if not np.isfinite(h):
        return '—'
    for unit, size in (('d', 24.0), ('h', 1.0), ('min', 1.0 / 60.0), ('s', 1.0 / 3600.0)):
        if h >= size or unit == 's':
            v = h / size
            return f'{v:.0f} {unit}' if v >= 10 else f'{v:.2g} {unit}'


def plot_feasibility_scan(scan, bounds, current_g: dict = None):
    """Rows: control experiments; columns: families. Each panel: predicted observable vs g, window shaded.

    A flat line means the experiment says nothing about that family. TMSPA conversion is drawn as TMSPA
    remaining so that every row shares a log axis.
    """
    families = list(dict.fromkeys(scan['family']))
    experiments = list(dict.fromkeys(scan['experiment']))
    fig, axes = plt.subplots(len(experiments), len(families), figsize=(3.1 * len(families), 2.3 * len(experiments)),
                             sharex=True, sharey=True, layout='constrained', squeeze=False)
    for i, exp_id in enumerate(experiments):
        for j, family in enumerate(families):
            ax = axes[i, j]
            d = scan[(scan['experiment'] == exp_id) & (scan['family'] == family)].sort_values('g_eV')
            obs = d['observable'].iloc[0]
            y, low, high = d['predicted'].values, d['low'].iloc[0], d['high'].iloc[0]
            if obs == 'tmspa_conversion':
                y, low, high = 1.0 - y, 1.0 - high if np.isfinite(high) else np.nan, 1.0 - low if np.isfinite(low) else np.nan
            ax.axhspan(low if np.isfinite(low) else 1e-5, high if np.isfinite(high) else 2.0, color=GRID, zorder=0)
            ax.plot(d['g_eV'], np.clip(y, 1e-5, None), color=FAMILY_COLORS[family], lw=2.0, zorder=2)
            for _, b in bounds[(bounds['experiment'] == exp_id) & (bounds['family'] == family)].iterrows():
                ax.axvline(b['g_eV'], color=INK_SOFT, lw=1.0, ls=':', zorder=1)
                ax.annotate(f"{b['bound']} {b['g_eV']:.2f}", (b['g_eV'], 2e-5), xytext=(3, 0),
                            textcoords='offset points', fontsize=7.5, color=INK_SOFT, rotation=90, va='bottom')
            if current_g and family in current_g:
                ax.axvline(current_g[family], color=RULE, lw=1.0, zorder=1)
            ax.set_yscale('log')
            ax.set_ylim(1e-5, 2.0)
            if i == 0:
                ax.set_title(_FAMILY_TITLES.get(family, family), loc='left', fontsize=10)
            if j == 0:
                exp = d.iloc[0]
                ax.set_ylabel(f"{exp_id}: {_OBSERVABLE_AXIS[obs]}", fontsize=9)
            if i == len(experiments) - 1:
                ax.set_xlabel('intrinsic barrier g (eV)')
    fig.suptitle('Each control experiment against the barrier of each family (others fixed). '
                 'Grey band: observed window; dotted: bound; solid grey: current value',
                 x=0.01, ha='left', fontsize=10.5)
    return fig


def plot_entropy_constraints(projection, family: str, current=None):
    """Feasible region of (ΔS‡, g at 25 °C) for one family: each bound is a line; the region between is shaded."""
    fig, ax = plt.subplots(figsize=(7.2, 4.4), layout='constrained')
    dS = np.array(sorted(projection['dS_J_mol_K'].unique()))
    lower = np.full(len(dS), -np.inf)
    upper = np.full(len(dS), np.inf)
    k = 0
    for (exp_id, side, T_K), d in projection.groupby(['experiment', 'bound', 'T_K'], sort=False):
        d = d.sort_values('dS_J_mol_K')
        g = d['g_ref_eV'].values
        if side == 'g ≥':
            lower = np.maximum(lower, g)
        else:
            upper = np.minimum(upper, g)
        ax.plot(dS, g, color=SERIES[k % len(SERIES)], lw=1.8,
                label=f"{exp_id} at {T_K - 273.15:.0f} °C: {side.replace('g', 'g(T)')}")
        k += 1
    ok = lower < upper
    ax.fill_between(dS, lower, upper, where=ok, color=BLUE_RAMP[0], alpha=0.35, lw=0, label='allowed by all', zorder=0)
    if current is not None:
        ax.plot(*current, 'o', color=INK, ms=7, zorder=4, label='current model (ΔS‡ = 0)')
    ax.set_xlabel('activation entropy ΔS‡ (J/mol/K)')
    ax.set_ylabel('g at 25 °C (eV)')
    ax.set_title(f'{_FAMILY_TITLES.get(family, family)}: what two temperatures allow', loc='left')
    ax.legend(loc='upper left', fontsize=8)
    return fig


def plot_water_series_check(water, groups: dict = None, box=((0.5, 1.0), (0.0, 0.05)), marks_h=(1.0, 24.0)):
    """TMSPA left at 1 vol% water (x) against TMSPA left at 2 vol% (y) at a common time, one line per model.

    Each line runs through time; the shaded box is what Gogoi 2024 report for the two samples.
    groups: {model: group label}; models of one group share a colour and one legend entry.
    """
    fig, ax = plt.subplots(figsize=(7.2, 4.8), layout='constrained')
    (x_lo, x_hi), (y_lo, y_hi) = box
    ax.fill_between([x_lo, x_hi], max(y_lo, 1e-4), y_hi, color=BLUE_RAMP[0], alpha=0.4, lw=0, zorder=0)
    ax.annotate('observed\n(Gogoi 2024)', ((x_lo + x_hi) / 2, 3e-3), ha='center', fontsize=9, color=INK_SOFT)
    x_ref = np.linspace(0.0, 1.0, 100)
    ax.plot(x_ref, np.clip(x_ref ** 2, 1e-4, None), color=INK_MUTED, lw=1.2, ls='--', zorder=3,
            label='first order in water: y = x²')
    group_order = list(dict.fromkeys((groups or {}).get(m, m) for m in water['model'].unique()))
    labelled = set()
    for name, d in water.groupby('model', sort=False):
        group = (groups or {}).get(name, name)
        color = SERIES[group_order.index(group) % len(SERIES)]
        x, y = d['W1 TMSPA'].values, np.clip(d['W2 TMSPA'].values, 1e-4, None)
        ax.plot(x, y, color=color, lw=1.8, label=None if group in labelled else group, zorder=2)
        labelled.add(group)
        for t in marks_h:
            i = np.argmin(np.abs(d['t_h'].values - t))
            ax.plot(x[i], y[i], 'o', color=color, ms=5, mec=SURFACE, mew=1.2, zorder=4)
    ax.set_yscale('log')
    ax.set_ylim(1e-4, 1.2)
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel('TMSPA left, 1 vol% H2O (fraction of P)')
    ax.set_ylabel('TMSPA left, 2 vol% H2O (fraction of P)')
    ax.set_title('Water series at a common time (dots: 1 h and 24 h after mixing)', loc='left')
    ax.legend(loc='upper left', fontsize=8)
    return fig


def plot_half_life_maps(maps: dict, window_h=(1.0 / 6.0, 12.0), species: str = 'TMSPA'):
    """Half-life over (g, T), one panel per condition; cells inside the measurable window are outlined."""
    fig, axes = plt.subplots(1, len(maps), figsize=(5.4 * len(maps), 4.6), layout='constrained', squeeze=False)
    lo, hi = np.log10(window_h[0]) - 2.0, np.log10(window_h[1]) + 2.0
    for ax, (title, table) in zip(axes[0], maps.items()):
        values = np.log10(table.values.astype(float))
        mesh = ax.pcolormesh(np.arange(table.shape[1] + 1), np.arange(table.shape[0] + 1),
                             np.clip(values, lo, hi), cmap=SEQUENTIAL_CMAP, vmin=lo, vmax=hi, shading='flat')
        for i in range(table.shape[0]):
            for j in range(table.shape[1]):
                h = table.values[i, j]
                inside = np.isfinite(h) and window_h[0] <= h <= window_h[1]
                dark = np.isfinite(values[i, j]) and values[i, j] > (lo + hi) / 2
                ax.text(j + 0.5, i + 0.5, _format_hours(h) if np.isfinite(h) else '> max', ha='center', va='center',
                        fontsize=7.5, color=SURFACE if dark else INK, fontweight='bold' if inside else 'normal')
                if inside:
                    ax.add_patch(plt.Rectangle((j + 0.04, i + 0.04), 0.92, 0.92, fill=False, ec=INK, lw=1.2))
        ax.set_xticks(np.arange(table.shape[1]) + 0.5, [f'{c:g}' for c in table.columns])
        ax.set_yticks(np.arange(table.shape[0]) + 0.5, [f'{g:.2f}' for g in table.index])
        ax.set_xlabel('temperature (°C)')
        ax.set_ylabel(table.index.name or 'g (eV)')
        ax.set_title(title, loc='left', fontsize=10)
        ax.grid(False)
    cbar = fig.colorbar(mesh, ax=axes[0].tolist(), shrink=0.85)
    cbar.set_label(f'log10 half-life of {species} (h)')
    fig.suptitle(f'Outlined, bold: measurable ({_format_hours(window_h[0])} to {_format_hours(window_h[1])})',
                 x=0.01, ha='left', fontsize=10)
    return fig


def plot_design_precision(scans: dict, floor_meV: float = 2.0, undetermined_meV: float = 100.0):
    """σ of each parameter against its assumed true value, one line per design.

    scans: {panel title: (DataFrame from scan_design_precision filtered to one parameter, x label)}.
    """
    fig, axes = plt.subplots(1, len(scans), figsize=(4.3 * len(scans), 4.2), sharey=True, layout='constrained',
                             squeeze=False)
    designs = list(dict.fromkeys(next(iter(scans.values()))[0]['design']))
    for ax, (title, (d, x_label)) in zip(axes[0], scans.items()):
        ax.axhspan(0.01, floor_meV, color=GRID, zorder=0)
        ax.axhspan(undetermined_meV, 1e3, color=GRID, zorder=0)
        for k, design in enumerate(designs):
            dd = d[d['design'] == design].sort_values('true value')
            ax.plot(dd['true value'], np.clip(dd['sigma'] * 1000.0, 0.05, None), marker='o', ms=5,
                    color=SERIES[k % len(SERIES)], lw=1.8, label=design)
        ax.set_yscale('log')
        ax.set_ylim(0.05, 500.0)
        ax.set_xlabel(x_label)
        ax.set_title(title, loc='left', fontsize=10)
    axes[0, 0].set_ylabel('expected σ of the barrier (meV)')
    axes[0, 0].annotate('below the systematic floor', (0.02, 0.06), xycoords=('axes fraction', 'data'),
                        fontsize=7.5, color=INK_MUTED)
    axes[0, 0].annotate('not determined', (0.02, 320.0), xycoords=('axes fraction', 'data'), fontsize=7.5,
                        color=INK_MUTED)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=len(labels), fontsize=8.5)
    return fig


def plot_design_readouts(results: dict, measured: dict = None, acquisition: dict = None, nuclei=('P', 'Si')):
    """Temperature program and the readouts of each nucleus; lines per truth, dots: one synthetic measurement.

    results: {label: simulate_design result}; measured: {label: add_readout_noise output} (drawn for the first).
    """
    first = next(iter(results.values()))
    fig, axes = plt.subplots(len(nuclei) + 1, 1, figsize=(10.5, 2.0 + 3.0 * len(nuclei)), sharex=True,
                             height_ratios=[0.8] + [3.0] * len(nuclei), layout='constrained')
    add_temperature_strip(axes[0], first['t_h'], first['T_C'] + 273.15)
    styles = ['-', '--', ':', '-.']
    for ax, nucleus in zip(axes[1:], nuclei):
        labels = list(first['readouts'][nucleus])
        for k, label in enumerate(labels):
            color = SERIES[k % len(SERIES)]
            for s, (name, res) in enumerate(results.items()):
                trace = calculate_readouts(res['C_M'], res['idx'], nucleus)[label]
                ax.plot(res['t_h'], trace * 1000.0, color=color, lw=1.6, ls=styles[s % len(styles)],
                        label=label if s == 0 else None, zorder=2)
            if measured is not None:
                first_name = next(iter(results))
                y = measured[first_name][nucleus][label] * 1000.0
                t = np.asarray(first['sampling_h'][nucleus])
                ax.plot(t, y, 'o', color=color, ms=3.2, alpha=0.75, mec='none', zorder=3)
        ax.set_ylabel(f'{nucleus} in peak (mM)')
        ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1.0), fontsize=8, title=f'{nucleus} peaks',
                  title_fontsize=8)
    if len(results) > 1:
        names = list(results)
        axes[1].set_title('   '.join(f"{'─' if s == 0 else '- -' if s == 1 else '···'} {n}" for s, n in enumerate(names)),
                          loc='left', fontsize=9, color=INK_SOFT)
    axes[-1].set_xlabel('time after mixing (h)')
    return fig
