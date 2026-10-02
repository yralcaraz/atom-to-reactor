"""Figures for the feasible region (notebook 03): allowed barrier intervals, boundary maps, time scenarios and
the structure tests (composition paths, heating, equilibrium locus)."""

import textwrap

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Patch

from demo.style import FAMILY_COLORS, GRID, INK, INK_MUTED, INK_SOFT, RULE, SERIES, SURFACE
from kinetics.fitting import summarize_region

FAMILY_TITLES = {'hydrolysis': 'hydrolysis (R1–R3)', 'transfer': 'transfer (R5–R7)',
                 'condensation': 'condensation (R4)', 'solvent_attack': 'solvent attack (R8–R9)'}
_FALLBACK_COLOR = SERIES[6]


def _family_color(probes: str) -> str:
    return FAMILY_COLORS.get(probes.split(', ')[0], _FALLBACK_COLOR)


def plot_allowed_intervals(panels: dict, labels: dict, *, g_range=(0.60, 1.70), ncols: int = 2, title: str = None):
    """One panel per entry {panel title: (intervals, intersection row)}; one row per observation.

    A bar spans the barrier values the observation allows (colour = family it probes); an arrow means the
    observation leaves that side open within the scanned range. The shaded column is what every observation
    allows together; when it is empty, the panel says which observations exclude each other.
    """
    n_rows = max(len(iv) for iv, _ in panels.values())
    ncols = min(ncols, len(panels))
    nrows = int(np.ceil(len(panels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.8 * ncols, (0.5 * n_rows + 1.4) * nrows + 0.7),
                             layout='constrained', squeeze=False)
    for ax in axes.flat[len(panels):]:
        ax.set_visible(False)
    for ax, (panel_title, (iv, inter)) in zip(axes.flat, panels.items()):
        iv = iv.reset_index(drop=True)
        y = np.arange(len(iv))[::-1]
        if not inter['empty']:
            ax.axvspan(inter['g_low_eV'], inter['g_high_eV'], color=GRID, zorder=0)
        for yi, (_, r) in zip(y, iv.iterrows()):
            color = _family_color(r['probes'])
            if r['verdict'].startswith('no value'):
                ax.text(np.mean(g_range), yi, 'no value allowed', ha='center', va='center', fontsize=8.5,
                        color=INK_SOFT, style='italic')
                continue
            lo = r['g_low_eV'] if np.isfinite(r['g_low_eV']) else g_range[0]
            hi = r['g_high_eV'] if np.isfinite(r['g_high_eV']) else g_range[1]
            ax.plot([lo, hi], [yi, yi], color=color, lw=5, solid_capstyle='round', zorder=2)
            closed = [np.isfinite(r['g_low_eV']), np.isfinite(r['g_high_eV'])]
            for edge, is_closed, marker in ((lo, closed[0], '<'), (hi, closed[1], '>')):
                if not is_closed:
                    ax.plot(edge, yi, marker=marker, color=color, ms=7, zorder=3)
            if all(closed) and hi - lo < 0.12:     # narrow: one label for both edges
                labels_at = [(0.5 * (lo + hi), f'{lo:.3f}–{hi:.3f}')]
            else:
                labels_at = [(e, f'{e:.3f}') for e, c in zip((lo, hi), closed) if c]
            for x_text, text in labels_at:
                ax.annotate(text, (x_text, yi), xytext=(0, 7), textcoords='offset points', ha='center',
                            fontsize=7.5, color=INK_SOFT)
        ax.set_yticks(y)
        ax.set_yticklabels([labels.get(o, o) for o in iv['observation']], fontsize=8.5)
        ax.set_xlim(g_range[0] - 0.03, g_range[1] + 0.03)
        ax.set_ylim(-0.7, len(iv) - 0.3)
        ax.grid(axis='y', visible=False)
        ax.set_xlabel('intrinsic barrier g (eV)')
        ax.set_title(panel_title + '\n' + textwrap.fill(summarize_region(iv), 62), loc='left', fontsize=10)
    handles = [Line2D([], [], color=FAMILY_COLORS[f], lw=5, label=f'probes {FAMILY_TITLES[f]}') for f in FAMILY_COLORS]
    handles += [Line2D([], [], color=INK_SOFT, marker='>', lw=0, label='open: not bounded on that side'),
                Patch(facecolor=GRID, label='allowed by every observation')]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=8.5, frameon=False)
    if title:
        fig.suptitle(title, x=0.01, ha='left', fontsize=10.5)
    return fig


def plot_boundary_map(edges, labels: dict, *, y_bounds: dict = None, title: str = None, ax=None):
    """Plane of two barriers: each observation shades the side of its edge that it excludes (darker where
    several do); what stays white is allowed by every observation.

    edges: table of kinetics.fitting.trace_boundaries (edges in x at each y value).
    y_bounds: {label: (low, high)} for observations that bound only the y barrier (horizontal edges).
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(6.8, 5.2), layout='constrained')
    x_range = (edges['g_low_eV'].min() - 0.1 if edges['g_low_eV'].notna().any() else 0.6,
               edges['g_high_eV'].max() + 0.1 if edges['g_high_eV'].notna().any() else 1.7)
    y_values = np.sort(edges['y_eV'].unique())
    y_range = (y_values[0], y_values[-1])
    shade = dict(color=INK_MUTED, alpha=0.2, lw=0)
    allowed_lo = np.full(len(y_values), x_range[0])
    allowed_hi = np.full(len(y_values), x_range[1])
    handles = []
    for k, (obs_id, grp) in enumerate(edges.groupby('observation', sort=False)):
        grp = grp.set_index('y_eV').reindex(y_values)
        color = SERIES[k % len(SERIES)]
        lo, hi = grp['g_low_eV'].to_numpy(), grp['g_high_eV'].to_numpy()
        none = grp['verdict'].str.startswith('no value').to_numpy()
        x_lo = np.where(none, x_range[1], np.where(np.isfinite(lo), lo, x_range[0]))
        x_hi = np.where(none, x_range[1], np.where(np.isfinite(hi), hi, x_range[1]))
        ax.fill_betweenx(y_values, x_range[0], x_lo, **shade)
        ax.fill_betweenx(y_values, x_hi, x_range[1], **shade)
        for edge in (lo, hi):
            if np.isfinite(edge).any():
                ax.plot(edge, y_values, color=color, lw=2.0, marker='o', ms=3)
        allowed_lo, allowed_hi = np.maximum(allowed_lo, x_lo), np.minimum(allowed_hi, x_hi)
        handles.append(Line2D([], [], color=color, lw=2.0, label=labels.get(obs_id, obs_id)))
    in_y = np.ones(len(y_values), dtype=bool)
    for j, (label, (low, high)) in enumerate((y_bounds or {}).items()):
        color = SERIES[(len(handles) + j) % len(SERIES)]
        for value, span in ((low, (y_range[0], low)), (high, (high, y_range[1]))):
            if value is not None and np.isfinite(value):
                ax.axhspan(*span, **shade)
                ax.axhline(value, color=color, lw=2.0)
        in_y &= (y_values >= (low if low is not None and np.isfinite(low) else -np.inf)) & \
                (y_values <= (high if high is not None and np.isfinite(high) else np.inf))
        handles.append(Line2D([], [], color=color, lw=2.0, label=label))
    left = bool(np.any((allowed_lo < allowed_hi) & in_y))
    ax.set_xlim(x_range)
    ax.set_ylim(y_range)
    ax.set_xlabel(f"g of {edges['parameter'].iloc[0]} (eV)")
    ax.set_ylabel(f"g of {edges['y_parameter'].iloc[0]} (eV)")
    ax.set_title((title + '\n' if title else '') + ('a region is left (white)' if left else
                 'no point is allowed by every observation'), loc='left', fontsize=10)
    ax.legend(handles=handles + [Patch(facecolor=INK_MUTED, alpha=0.3, label='excluded (darker: by more observations)')],
              loc='upper center', bbox_to_anchor=(0.5, -0.14), fontsize=8.5, frameon=False)
    return ax


def plot_time_scenarios(panels: dict):
    """Bounds on a barrier as a function of an assumed time (log x).

    panels: {title: {'x': times_h, 'xlabel': str, 'lines': {label: (y, 'upper'|'lower')}, 'marks_h': {label: t}}}.
    Upper bounds are drawn solid, lower bounds dashed; where the lowest upper bound lies above the highest
    lower bound the band is shaded (values allowed by both).
    """
    fig, axes = plt.subplots(1, len(panels), figsize=(4.8 * len(panels), 4.0), layout='constrained', squeeze=False)
    for ax, (panel_title, spec) in zip(axes[0], panels.items()):
        x = np.asarray(spec['x'], dtype=float)
        uppers, lowers = [], []
        for k, (label, (y, kind)) in enumerate(spec['lines'].items()):
            y = np.broadcast_to(np.asarray(y, dtype=float), x.shape)
            ax.plot(x, y, color=SERIES[k], lw=2.0, ls='-' if kind == 'upper' else '--', marker='o', ms=4,
                    label=f"{label} ({'upper' if kind == 'upper' else 'lower'} bound)")
            (uppers if kind == 'upper' else lowers).append(y)
        if uppers and lowers:
            top, bottom = np.min(uppers, axis=0), np.max(lowers, axis=0)
            ax.fill_between(x, bottom, top, where=top >= bottom, interpolate=True, color=GRID, zorder=0,
                            label='allowed by both')
        for label, t in spec.get('marks_h', {}).items():
            ax.axvline(t, color=RULE, lw=1.0, zorder=0)
            ax.annotate(label, (t, 1.0), xycoords=('data', 'axes fraction'), xytext=(3, -3),
                        textcoords='offset points', fontsize=8, color=INK_SOFT, va='top', rotation=90)
        ax.set_xscale('log')
        ax.set_xticks(x, [f'{v:g}' for v in x])
        ax.minorticks_off()
        ax.set_xlabel(spec['xlabel'])
        ax.set_ylabel('intrinsic barrier g (eV)')
        ax.set_title(panel_title, loc='left', fontsize=10)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.18), fontsize=8, frameon=False)
    return fig


def plot_composition_paths(paths: dict, measured, *, x: str = 'H3PO4', ys=('BMSPA', 'MMSPA'), title: str = None):
    """Composition paths (share of all P) of model runs against measured spectra.

    paths: {label: DataFrame of shares along a run}; measured: DataFrame with '<species>' and
    '<species>_sigma' columns, one row per spectrum. A barrier moves a run along its path but does not
    change the path, so a measured point off every path is out of reach of the barrier.
    """
    fig, axes = plt.subplots(1, len(ys), figsize=(5.2 * len(ys), 4.4), layout='constrained', squeeze=False)
    for ax, y in zip(axes[0], ys):
        for k, (label, path) in enumerate(paths.items()):
            ax.plot(path[x], path[y], color=SERIES[k % len(SERIES)], lw=2.0, label=label)
        ax.errorbar(measured[x], measured[y], xerr=2 * measured[f'{x}_sigma'], yerr=2 * measured[f'{y}_sigma'],
                    fmt='o', ms=6, color=INK, mfc=SURFACE, mew=1.5, elinewidth=1.0, capsize=2, zorder=5,
                    label='measured (± 2 SE)')
        ax.set_xlabel(f'{x} share of P')
        ax.set_ylabel(f'{y} share of P')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, None)
    axes[0][0].legend(loc='upper center', bbox_to_anchor=(0.5 * len(ys) + 0.1, -0.16), ncol=2, fontsize=8.5,
                      frameon=False)
    if title:
        fig.suptitle(title, x=0.01, ha='left', fontsize=10.5)
    return fig


def plot_rt_equivalent(dG_eV, days, *, references: dict = None, title: str = None, ax=None):
    """Room-temperature days that a heating history is worth, against the barrier ΔG‡ (log y)."""
    if ax is None:
        _, ax = plt.subplots(figsize=(6.4, 4.0), layout='constrained')
    ax.plot(dG_eV, days, color=SERIES[0], lw=2.0, marker='o', ms=4, label='minimum heating from the file times')
    for k, (label, value) in enumerate((references or {}).items()):
        ax.axhline(value, color=INK_MUTED, lw=1.0, ls=['--', ':', '-.'][k % 3], label=label)
    ax.set_yscale('log')
    ax.set_xlabel('barrier ΔG‡ of the step (eV)')
    ax.set_ylabel('equivalent time at room temperature (days)')
    if title:
        ax.set_title(title, loc='left', fontsize=10)
    ax.legend(loc='upper left', fontsize=8.5)
    return ax


def _covariance_ellipse(ax, center, cov, n_sigma, **style):
    vals, vecs = np.linalg.eigh(cov)
    angle = np.degrees(np.arctan2(vecs[1, -1], vecs[0, -1]))
    w, h = 2 * n_sigma * np.sqrt(vals[-1]), 2 * n_sigma * np.sqrt(vals[0])
    ax.add_patch(Ellipse(center, w, h, angle=angle, fill=False, **style))


def plot_equilibrium_locus(locus, center: dict, cov, closest: dict, *, x: str = 'R2', y: str = 'R3',
                           title: str = None, ax=None):
    """ΔG_rxn plane: the computed value with its 1 and 2 SE ellipses, and the locus of values for which the
    measured composition would be an equilibrium (one point per possible amount of HMDSO formed)."""
    if ax is None:
        _, ax = plt.subplots(figsize=(6.0, 5.2), layout='constrained')
    c = (center[x], center[y])
    C = cov.loc[[x, y], [x, y]].to_numpy()
    for n, ls in ((1, '-'), (2, '--')):
        _covariance_ellipse(ax, c, C, n, edgecolor=SERIES[0], lw=1.5, ls=ls)
    ax.plot(*c, 'o', color=SERIES[0], ms=8, label='computed (DFT + MACE), with 1 and 2 SE')
    ax.plot(locus[f'dG_{x}_eV'], locus[f'dG_{y}_eV'], color=SERIES[1], lw=2.0,
            label='equilibrium locus of the measured plateau')
    ax.plot(locus[f'dG_{x}_eV'].iloc[0], locus[f'dG_{y}_eV'].iloc[0], 's', color=SERIES[1], ms=7)
    ax.annotate('no HMDSO\n(most TMSOH)', (locus[f'dG_{x}_eV'].iloc[0], locus[f'dG_{y}_eV'].iloc[0]),
                xytext=(-8, 8), textcoords='offset points', ha='right', fontsize=8, color=INK_SOFT)
    row = closest['row']
    ax.plot(row[f'dG_{x}_eV'], row[f'dG_{y}_eV'], 'D', color=INK, ms=6, mfc=SURFACE, mew=1.5,
            label=f"closest point: {closest['distance']:.2f} SE from computed")
    ax.axhline(0, color=RULE, lw=1.0)
    ax.axvline(0, color=RULE, lw=1.0)
    ax.set_xlabel(f'ΔG_rxn {x} (eV)')
    ax.set_ylabel(f'ΔG_rxn {y} (eV)')
    ax.set_aspect('equal', adjustable='datalim')
    if title:
        ax.set_title(title, loc='left', fontsize=10)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), fontsize=8.5, frameon=False)
    return ax
