"""Figures for the first fit (notebook 05): residual maps, parameters across scenarios, profiles, influence of
each sample, recovery of known parameters, sample trajectories and the protocol band.

Colour does one job per figure: a diverging blue–grey–red scale for signed residuals and shifts (grey = none),
one hue per age scenario in a fixed order, a single blue ramp for an ordered quantity (ΔS‡).

Classification: PUBLIC (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D

from demo.style import BLUE_RAMP, GRID, INK, INK_MUTED, INK_SOFT, RULE, SERIES, SURFACE

SCENARIO_ORDER = ('short', 'middle', 'long', 'free')
SCENARIO_LABELS = {'short': 'short (1 h)', 'middle': 'middle (1 d)', 'long': 'long (7 d)', 'free': 'free ages'}
SCENARIO_COLORS = dict(zip(SCENARIO_ORDER, SERIES[:4]))
_NEUTRAL = '#f0efec'
DIVERGING = LinearSegmentedColormap.from_list('blue_grey_red', [SERIES[0], _NEUTRAL, SERIES[7]])


def _cell_text(ax, x, y, text, strong=False):
    ax.text(x, y, text, ha='center', va='center', fontsize=7.5, color=INK if strong else INK_SOFT,
            fontweight='bold' if strong else 'normal')


def collapse_residuals(residuals):
    """Residual table with the drift rows of a repeated sample reduced to the largest one per window."""
    drift = residuals[residuals['kind'] == 'drift'].copy()
    rest = residuals[residuals['kind'] != 'drift'].copy()
    if len(drift):
        drift['quantity'] = drift['quantity'].str.replace(r', spectrum \d+', ' (largest)', regex=True)
        keys = [c for c in ('structure', 'scenario', 'model', 'sample', 'quantity') if c in drift.columns]
        drift = drift.loc[drift.groupby(keys)['z'].apply(lambda z: z.abs().idxmax())]
    return rest, drift


def plot_residual_map(residuals, *, scenario: str, columns: str = 'structure', order=None, clip: float = 6.0,
                      title: str = None):
    """Standardised residuals z of every measured quantity (rows) for each structure (columns) in one scenario.

    Blue: the model predicts less than measured; red: more; grey: within the error. The number is z; bold where
    |z| > 3, which is what the fit rule reads. Drift rows show the largest drift of each window.
    """
    data = residuals[residuals['scenario'] == scenario].copy()
    # A sample with several independent spectra has one row per spectrum
    shares = data['kind'] == 'share'
    groups = data[shares].groupby([columns, 'sample', 'quantity'])['t_h']
    number, count = groups.rank(method='first').astype(int), groups.transform('size')
    data.loc[shares, 'quantity'] = [f'{q}, spectrum {n}' if c > 1 else q
                                    for q, n, c in zip(data.loc[shares, 'quantity'], number, count)]
    rest, drift = collapse_residuals(data)
    data = pd.concat([rest, drift]) if len(drift) else rest
    samples = list(dict.fromkeys(rest['sample']))               # drift rows follow the composition of their sample
    data = data.assign(row=data['sample'] + ' · ' + data['quantity'],
                       rank=data['sample'].map(samples.index) + 0.5 * (data['kind'] == 'drift'))
    rows = list(dict.fromkeys(data.sort_values('rank', kind='stable')['row']))
    cols = [c for c in (order or list(dict.fromkeys(data[columns]))) if c in set(data[columns])]
    z = data.pivot_table(index='row', columns=columns, values='z', aggfunc='first').reindex(index=rows, columns=cols)
    fig, ax = plt.subplots(figsize=(max(4.2 + 1.05 * len(cols), 8.0), 0.9 + 0.27 * len(rows)), layout='constrained')
    image = ax.imshow(np.clip(z.to_numpy(dtype=float), -clip, clip), cmap=DIVERGING, norm=TwoSlopeNorm(0.0, -clip, clip),
                      aspect='auto')
    for i in range(len(rows)):
        for j in range(len(cols)):
            v = z.iat[i, j]
            if np.isfinite(v):
                _cell_text(ax, j, i, f'{v:+.1f}' if abs(v) < 99.5 else f'{v:+.0f}', strong=abs(v) > 3.0)
    ax.set_xticks(range(len(cols)), cols)
    ax.set_yticks(range(len(rows)), rows, fontsize=8)
    ax.tick_params(length=0, top=True, labeltop=True, bottom=False, labelbottom=False)
    ax.set_xticks(np.arange(-0.5, len(cols)), minor=True)
    ax.set_yticks(np.arange(-0.5, len(rows)), minor=True)
    ax.grid(False)
    ax.grid(which='minor', color=SURFACE, linewidth=2)
    ax.tick_params(which='minor', length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    bar = fig.colorbar(image, ax=ax, shrink=0.5, pad=0.02, extend='both')
    bar.set_label('z = (predicted − measured) / σ', color=INK_SOFT, fontsize=8.5)
    bar.outline.set_visible(False)
    fig.suptitle(title or f'Standardised residuals, {SCENARIO_LABELS.get(scenario, scenario)} scenario', x=0.01, ha='left',
                 fontsize=11)
    return ax


def plot_fit_across_scenarios(parameters, intervals=None, *, structure: str, names=None, units: dict = None,
                              boxes: dict = None, title: str = None):
    """Best-fit value of each parameter in each scenario, with its 95 % profile interval where one was computed.

    A bar is the interval; an arrow is a side the data leave open inside the profile grid; a dot alone has no
    profile. boxes: {parameter: (lower, upper)} draws the declared search box as the panel's limits.
    """
    fits = parameters[parameters['structure'] == structure].set_index('scenario')
    names = names or [c for c in fits.columns if c != 'structure' and fits[c].notna().any()]
    ncols = min(5, len(names))
    nrows = int(np.ceil(len(names) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.5 * ncols, 2.5 * nrows + 0.4), layout='constrained', squeeze=False)
    ci = None
    if intervals is not None:
        ci = intervals[(intervals['structure'] == structure) & (intervals['tag'].fillna('') == '')].set_index(['scenario', 'parameter'])
    for ax, name in zip(axes.flat, names):
        for k, scenario in enumerate(s for s in SCENARIO_ORDER if s in fits.index):
            value = fits.loc[scenario, name]
            if not np.isfinite(value):
                continue
            color = SCENARIO_COLORS[scenario]
            if ci is not None and (scenario, name) in ci.index:
                row = ci.loc[(scenario, name)]
                low = row['low'] if np.isfinite(row['low']) else row['grid_low']
                high = row['high'] if np.isfinite(row['high']) else row['grid_high']
                ax.plot([k, k], [low, high], color=color, linewidth=3, alpha=0.45, solid_capstyle='butt')
                for edge, open_side, marker in ((low, not np.isfinite(row['low']), 'v'), (high, not np.isfinite(row['high']), '^')):
                    if open_side:
                        ax.plot(k, edge, marker=marker, color=color, markersize=6, alpha=0.8)
            ax.plot(k, value, 'o', color=color, markersize=6.5, markeredgecolor=SURFACE, markeredgewidth=1.2)
        ax.set_xticks(range(len(SCENARIO_ORDER)), SCENARIO_ORDER, fontsize=8)
        ax.set_xlim(-0.6, len(SCENARIO_ORDER) - 0.4)
        ax.set_title(name, fontsize=9.5, loc='left')
        ax.set_ylabel((units or {}).get(name, ''), fontsize=8.5)
        ax.grid(axis='x', visible=False)
    for ax in axes.flat[len(names):]:
        ax.set_visible(False)
    handles = [Line2D([], [], marker='o', linestyle='', color=INK_SOFT, label='best fit'),
               Line2D([], [], color=INK_MUTED, linewidth=3, alpha=0.6, label='95 % profile interval'),
               Line2D([], [], marker='^', linestyle='', color=INK_MUTED, label='open side: not bounded by the data')]
    fig.legend(handles=handles, loc='outside lower center', ncols=3)
    fig.suptitle(title or f'{structure}: fitted parameters in the four scenarios', x=0.01, ha='left', fontsize=11)
    return axes


def plot_profiles(profiles: dict, names=None, *, delta: float = 3.84, y_max: float = 12.0, title: str = None):
    """Profile of χ² along each parameter: {scenario: profile table}. The rule marks Δχ² = 3.84 (95 %)."""
    first = next(iter(profiles.values()))
    names = names or list(dict.fromkeys(first['parameter']))
    ncols = min(4, len(names))
    nrows = int(np.ceil(len(names) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.0 * ncols, 2.4 * nrows + 0.5), layout='constrained', squeeze=False)
    for ax, name in zip(axes.flat, names):
        for scenario, profile in profiles.items():
            p = profile[profile['parameter'] == name].sort_values('value')
            if p.empty:
                continue
            rise = p['chi2'] - p['chi2'].min()
            ax.plot(p['value'], rise, color=SCENARIO_COLORS.get(scenario, INK_SOFT), linewidth=1.8, marker='o', markersize=3)
        ax.axhline(delta, color=INK_MUTED, linewidth=0.9)
        ax.set_ylim(-0.4, y_max)
        ax.set_title(name, fontsize=9.5, loc='left')
        ax.set_ylabel('Δχ²', fontsize=8.5)
    for ax in axes.flat[len(names):]:
        ax.set_visible(False)
    handles = [Line2D([], [], color=SCENARIO_COLORS[s], linewidth=1.8, label=SCENARIO_LABELS[s]) for s in SCENARIO_ORDER if s in profiles]
    handles.append(Line2D([], [], color=INK_MUTED, linewidth=0.9, label=f'Δχ² = {delta:g} (95 %)'))
    fig.legend(handles=handles, loc='outside lower center', ncols=len(handles))
    fig.suptitle(title or 'Profiles: χ² re-minimised over the other parameters', x=0.01, ha='left', fontsize=11)
    return axes


def plot_influence(loo, fit_values: dict, *, columns: list, scale: float = 1000.0, unit: str = 'meV', clip: float = 100.0,
                   title: str = None):
    """Shift of each fitted quantity (columns) when one sample is left out of the fit (rows).

    loo: leave-one-out table of one structure and scenario; fit_values: {column: value with every sample}.
    """
    rows = list(loo['left out'])
    shift = np.array([[scale * (r[c] - fit_values[c]) for c in columns] for _, r in loo.iterrows()])
    fig, ax = plt.subplots(figsize=(max(3.4 + 0.95 * len(columns), 8.0), 1.2 + 0.36 * len(rows)), layout='constrained')
    image = ax.imshow(np.clip(shift, -clip, clip), cmap=DIVERGING, norm=TwoSlopeNorm(0.0, -clip, clip), aspect='auto')
    for i in range(len(rows)):
        for j in range(len(columns)):
            if np.isfinite(shift[i, j]):
                _cell_text(ax, j, i, f'{shift[i, j]:+.0f}', strong=abs(shift[i, j]) >= 0.25 * clip)
    ax.set_xticks(range(len(columns)), [c.replace('dG‡ ', 'ΔG‡ ') for c in columns], fontsize=8.5)
    ax.set_yticks(range(len(rows)), [f'without {r}' for r in rows], fontsize=8.5)
    ax.tick_params(length=0, top=True, labeltop=True, bottom=False, labelbottom=False)
    ax.set_xticks(np.arange(-0.5, len(columns)), minor=True)
    ax.set_yticks(np.arange(-0.5, len(rows)), minor=True)
    ax.grid(False)
    ax.grid(which='minor', color=SURFACE, linewidth=2)
    ax.tick_params(which='minor', length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    bar = fig.colorbar(image, ax=ax, shrink=0.7, pad=0.02, extend='both')
    bar.set_label(f'shift of the best fit ({unit})', color=INK_SOFT, fontsize=8.5)
    bar.outline.set_visible(False)
    fig.suptitle(title or 'Which sample drives which value', x=0.01, ha='left', fontsize=11)
    return ax


def plot_recovery(summary, title: str = None):
    """Synthetic recovery: per parameter, the 95 % interval of each case against the value that made the data."""
    names = list(dict.fromkeys(summary['parameter']))
    ncols = min(4, len(names))
    nrows = int(np.ceil(len(names) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.0 * ncols, 2.3 * nrows + 0.5), layout='constrained', squeeze=False)
    for ax, name in zip(axes.flat, names):
        rows = summary[summary['parameter'] == name].reset_index(drop=True)
        for i, r in rows.iterrows():
            low = r['low'] if np.isfinite(r['low']) else min(r['fit'], r['truth']) - 0.3
            high = r['high'] if np.isfinite(r['high']) else max(r['fit'], r['truth']) + 0.3
            ax.plot([low, high], [i, i], color=SERIES[0], linewidth=3, alpha=0.45, solid_capstyle='butt')
            if not np.isfinite(r['low']):
                ax.plot(low, i, marker='<', color=SERIES[0], markersize=6)
            if not np.isfinite(r['high']):
                ax.plot(high, i, marker='>', color=SERIES[0], markersize=6)
            ax.plot(r['fit'], i, 'o', color=SERIES[0], markersize=6, markeredgecolor=SURFACE, markeredgewidth=1.2)
            ax.plot(r['truth'], i, marker='|', color=INK, markersize=13, markeredgewidth=1.6)
        ax.set_yticks(range(len(rows)), rows['case'], fontsize=8)
        ax.invert_yaxis()
        ax.set_title(name, fontsize=9.5, loc='left')
        ax.grid(axis='y', visible=False)
    for ax in axes.flat[len(names):]:
        ax.set_visible(False)
    handles = [Line2D([], [], marker='|', linestyle='', color=INK, markersize=11, markeredgewidth=1.6, label='value that made the data'),
               Line2D([], [], marker='o', linestyle='', color=SERIES[0], label='best fit'),
               Line2D([], [], color=SERIES[0], linewidth=3, alpha=0.45, label='95 % profile interval (arrow: open side)')]
    fig.legend(handles=handles, loc='outside lower center', ncols=3)
    fig.suptitle(title or 'Synthetic recovery: known parameters refitted from shares with the lab design and errors',
                 x=0.01, ha='left', fontsize=11)
    return axes


def plot_protocol_band(protocol, *, structure: str, scenarios=None, title: str = None):
    """TMSPA left at each acquisition of Peter's protocol, for each activation entropy in the declared range.

    One panel per scenario; one line per ΔS‡, light to dark from the most negative to the most positive.
    """
    data = protocol[protocol['structure'] == structure]
    scenarios = [s for s in (scenarios or SCENARIO_ORDER) if s in set(data['scenario'])]
    entropies = sorted(data['dS (J/mol/K)'].unique())
    fig, axes = plt.subplots(1, len(scenarios), figsize=(3.3 * len(scenarios), 3.3), layout='constrained', sharey=True, squeeze=False)
    for ax, scenario in zip(axes.flat, scenarios):
        panel = data[data['scenario'] == scenario]
        curves = [panel[panel['dS (J/mol/K)'] == dS].sort_values('hold T (°C)') for dS in entropies]
        ax.fill_between(curves[0]['hold T (°C)'], np.min([c['TMSPA left'].to_numpy() for c in curves], axis=0),
                        np.max([c['TMSPA left'].to_numpy() for c in curves], axis=0), color=BLUE_RAMP[0], alpha=0.25, linewidth=0)
        for color, dS, curve in zip(BLUE_RAMP[1:], entropies, curves):
            ax.plot(curve['hold T (°C)'], curve['TMSPA left'], color=color, linewidth=1.8, marker='o', markersize=3.5)
        ax.set_title(SCENARIO_LABELS.get(scenario, scenario), fontsize=9.5, loc='left')
        ax.set_xlabel('hold before the acquisition (°C)')
        ax.set_ylim(-0.03, 1.03)
    axes.flat[0].set_ylabel('TMSPA left (fraction)')
    handles = [Line2D([], [], color=color, linewidth=1.8, label=f'ΔS‡ = {dS:+.0f}') for color, dS in zip(BLUE_RAMP[1:], entropies)]
    fig.legend(handles=handles, loc='outside lower center', ncols=len(handles), title='J mol⁻¹ K⁻¹', title_fontsize=8)
    fig.suptitle(title or f"Peter's protocol with the {structure} best fit: the band of the unknown activation entropy",
                 x=0.01, ha='left', fontsize=11)
    return axes


def plot_sample_trajectories(trajectories, *, structure: str, scenarios=('middle', 'long', 'free'), window: str = 'TMSPA',
                             title: str = None):
    """Share of one ³¹P window in the two water samples against the time since mixing, one panel per scenario.

    Solid: the fit as it is, the sample started from its recipe. Dashed: the same parameters with a trace of
    TMSOH at mixing. Markers: the measured share, at the age the scenario gives the sample.
    """
    data = trajectories[trajectories['structure'] == structure]
    scenarios = [s for s in scenarios if s in set(data['scenario'])]
    samples = list(dict.fromkeys(data['sample']))
    colors = dict(zip(samples, (SERIES[0], SERIES[1])))
    trace = data['TMSOH at mixing (M)'].max()
    fig, axes = plt.subplots(1, len(scenarios), figsize=(3.7 * len(scenarios), 3.4), layout='constrained', sharey=True, squeeze=False)
    for ax, scenario in zip(axes.flat, scenarios):
        panel = data[data['scenario'] == scenario]
        for sample in samples:
            rows = panel[panel['sample'] == sample]
            for level, dashes in ((0.0, '-'), (trace, (0, (4, 2)))):
                curve = rows[(rows['kind'] == 'model') & (rows['TMSOH at mixing (M)'] == level)]
                ax.plot(curve['t_h'], curve[window], color=colors[sample], linewidth=1.8, linestyle=dashes)
            measured = rows[rows['kind'] == 'measured']
            ax.plot(measured['t_h'], measured[window], linestyle='none', marker='o', markersize=6, markerfacecolor=colors[sample],
                    markeredgecolor=SURFACE, markeredgewidth=1.2, zorder=5)
        ax.set_xscale('log')
        ax.set_xlim(0.1, 2000.0)
        ax.set_ylim(-0.05, 1.05)
        ax.set_title(SCENARIO_LABELS.get(scenario, scenario), fontsize=9.5, loc='left')
        ax.set_xlabel('time since mixing (h)')
    axes.flat[0].set_ylabel(f'{window} share of the ³¹P area')
    handles = [Line2D([], [], color=colors[sample], linewidth=1.8, label=sample) for sample in samples]
    handles += [Line2D([], [], color=INK_SOFT, linewidth=1.8, label='from the recipe alone'),
                Line2D([], [], color=INK_SOFT, linewidth=1.8, linestyle=(0, (4, 2)), label=f'with {1e6 * trace:g} µM TMSOH at mixing'),
                Line2D([], [], color=INK_SOFT, linestyle='none', marker='o', markersize=6, label='measured, at the age of the scenario')]
    fig.legend(handles=handles, loc='outside lower center', ncols=len(handles))
    fig.suptitle(title or f'{structure} best fits: {window} in the two water samples since mixing', x=0.01, ha='left', fontsize=11)
    return axes
