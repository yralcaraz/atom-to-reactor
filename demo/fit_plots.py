"""Figures for the first fit (notebook 05): residual maps, parameters across scenarios, profiles, influence of
each sample, recovery of known parameters, sample trajectories and the protocol band; and, for the step-by-step
version of the notebook, compositions against measurements, the ladder of structures, barriers on a half-life
ruler, energy shifts, the trace scan and the storage ranges.

Colour does one job per figure: a diverging blue–grey–red scale for signed residuals and shifts (grey = none),
one hue per age scenario in a fixed order, a single blue ramp for an ordered quantity (ΔS‡), red for a predicted
share that breaks the fit rule.

Source: Y. Alcaraz Galván
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D

from demo.style import BLUE_RAMP, GRID, INK, INK_MUTED, INK_SOFT, RULE, SERIES, SURFACE
from kinetics.constants import H_SI, KB_EV, KB_SI, ZERO_CELSIUS_K

SCENARIO_ORDER = ('short', 'middle', 'long', 'free')
SCENARIO_LABELS = {'short': 'short (1 h)', 'middle': 'middle (1 d)', 'long': 'long (7 d)', 'free': 'free ages'}
AGE_LABELS = {'short': '1 h', 'middle': '1 day', 'long': '7 days', 'free': 'free ages'}     # a scenario named by its age
REACTION_LABELS = {'R1': 'R1 · first hydrolysis', 'R2': 'R2 · second hydrolysis', 'R3': 'R3 · third hydrolysis',
                   'R5': 'R5 · first transfer', 'R4': 'R4 · condensation', 'R8': 'R8 · solvent attack'}
HALF_LIVES_S = {'1 s': 1.0, '1 min': 60.0, '1 h': 3600.0, '1 day': 86400.0, '1 month': 2.592e6, '10 years': 3.156e8}
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


def plot_profiles(profiles: dict, names=None, *, delta: float = 3.84, y_max: float = 12.0, title: str = None,
                  labels: dict = None):
    """Profile of χ² along each parameter: {scenario: profile table}. The rule marks Δχ² = 3.84 (95 %)."""
    labels = labels or SCENARIO_LABELS
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
    handles = [Line2D([], [], color=SCENARIO_COLORS[s], linewidth=1.8, label=labels[s]) for s in SCENARIO_ORDER if s in profiles]
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


def plot_protocol_band(protocol, *, structure: str, scenarios=None, title: str = None, labels: dict = None):
    """TMSPA left at each acquisition of Peter's protocol, for each activation entropy in the declared range.

    One panel per scenario; one line per ΔS‡, light to dark from the most negative to the most positive.
    """
    labels = labels or SCENARIO_LABELS
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
        ax.set_title(labels.get(scenario, scenario), fontsize=9.5, loc='left')
        ax.set_xlabel('hold before the acquisition (°C)')
        ax.set_ylim(-0.03, 1.03)
    axes.flat[0].set_ylabel('TMSPA left (fraction)')
    handles = [Line2D([], [], color=color, linewidth=1.8, label=f'ΔS‡ = {dS:+.0f}') for color, dS in zip(BLUE_RAMP[1:], entropies)]
    fig.legend(handles=handles, loc='outside lower center', ncols=len(handles), title='J mol⁻¹ K⁻¹', title_fontsize=8)
    fig.suptitle(title or f"Peter's protocol with the {structure} best fit: the band of the unknown activation entropy",
                 x=0.01, ha='left', fontsize=11)
    return axes


def plot_sample_trajectories(trajectories, *, structure: str, scenarios=('middle', 'long', 'free'), window: str = 'TMSPA',
                             title: str = None, labels: dict = None, with_trace: bool = True):
    """Share of one ³¹P window in the two water samples against the time since mixing, one panel per scenario.

    Solid: the fit as it is, the sample started from its recipe. Dashed: the same parameters with a trace of
    TMSOH at mixing (left out if with_trace is False). Markers: the measured share, at the age the scenario gives
    the sample.
    """
    labels = labels or SCENARIO_LABELS
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
            for level, dashes in ((0.0, '-'), (trace, (0, (4, 2))))[:2 if with_trace else 1]:
                curve = rows[(rows['kind'] == 'model') & (rows['TMSOH at mixing (M)'] == level)]
                ax.plot(curve['t_h'], curve[window], color=colors[sample], linewidth=1.8, linestyle=dashes)
            measured = rows[rows['kind'] == 'measured']
            ax.plot(measured['t_h'], measured[window], linestyle='none', marker='o', markersize=6, markerfacecolor=colors[sample],
                    markeredgecolor=SURFACE, markeredgewidth=1.2, zorder=5)
        ax.set_xscale('log')
        ax.set_xlim(0.1, 2000.0)
        ax.set_ylim(-0.05, 1.05)
        ax.set_title(labels.get(scenario, scenario), fontsize=9.5, loc='left')
        ax.set_xlabel('time since mixing (h)')
    axes.flat[0].set_ylabel(f'{window} share of the ³¹P area')
    handles = [Line2D([], [], color=colors[sample], linewidth=1.8, label=sample) for sample in samples]
    handles.append(Line2D([], [], color=INK_SOFT, linewidth=1.8, label='from the recipe alone'))
    if with_trace:
        handles.append(Line2D([], [], color=INK_SOFT, linewidth=1.8, linestyle=(0, (4, 2)), label=f'with {1e6 * trace:g} µM TMSOH at mixing'))
    handles.append(Line2D([], [], color=INK_SOFT, linestyle='none', marker='o', markersize=6, label='measured, at the age of the scenario'))
    fig.legend(handles=handles, loc='outside lower center', ncols=len(handles))
    fig.suptitle(title or f'{structure} best fits: {window} in the two water samples since mixing', x=0.01, ha='left', fontsize=11)
    return axes


def tabulate_compositions(residuals):
    """Compositions held by the residual table of one fit: a row per sample, spectrum and window.

    The repeated sample is given by its mean composition. Drift rows and the pulls of the freed energies are left out.
    """
    data = residuals[residuals['kind'].isin(('share', 'mean', 'limit'))].copy()
    data['averaged'] = data['quantity'].str.extract(r'\((mean of \d+ spectra)\)')[0].fillna('')
    data['quantity'] = data['quantity'].str.replace(r' \(mean of \d+ spectra\)', '', regex=True)
    groups = data.groupby(['sample', 'quantity'])['t_h']
    return data.assign(spectrum=groups.rank(method='first').astype(int), spectra=groups.transform('size'))


def plot_compositions(cases: dict, panels: list, *, limit: float = 3.0, predicted: bool = True, title: str = None):
    """Measured and predicted area shares: one row per case, one panel per sample.

    cases: {row label: residual table of one model in one scenario}. panels: sample names, or (sample, n) for the
    n-th spectrum of a sample with several (default: the last). Dot: measured share ± 1 standard error; triangle: a
    published upper limit. Bar: predicted share, red where its miss |z| exceeds the fit rule, with the miss above it.
    """
    specs = [(p, None) if isinstance(p, str) else tuple(p) for p in panels]
    tables = {label: tabulate_compositions(table) for label, table in cases.items()}
    first = next(iter(tables.values()))
    widths = [max(2, int((first['sample'] == sample).sum() / first.loc[first['sample'] == sample, 'spectra'].max()))
              for sample, _ in specs]
    fig, axes = plt.subplots(len(tables), len(specs), figsize=(1.9 + 0.5 * sum(widths), 0.9 + 1.5 * len(tables)),
                             layout='constrained', sharey=True, squeeze=False, gridspec_kw={'width_ratios': widths})
    any_limit = False
    for row, (label, table) in zip(axes, tables.items()):
        for ax, (sample, number) in zip(row, specs):
            rows = table[table['sample'] == sample]
            spectra = int(rows['spectra'].max())
            rows = rows[rows['spectrum'] == (number or spectra)].reset_index(drop=True)
            x = np.arange(len(rows))
            if predicted:
                missed = rows['z'].abs() > limit
                ax.bar(x, rows['predicted'], width=0.62, color=np.where(missed, SERIES[7], BLUE_RAMP[0]), linewidth=0)
                for xi, r in rows[missed].iterrows():
                    ax.text(xi, max(r['predicted'], r['measured'] + r['sigma']) + 0.04, f"{r['z']:+.1f}" if abs(r['z']) < 9.95 else f"{r['z']:+.0f}",
                            ha='center', va='bottom', fontsize=7.5, color=INK, fontweight='bold')
            else:
                ax.vlines(x, 0.0, rows['measured'].clip(lower=0.0), color=INK_MUTED, linewidth=1.0)
            is_limit = (rows['kind'] == 'limit').to_numpy()
            any_limit = any_limit or is_limit.any()
            ax.errorbar(x[~is_limit], rows['measured'][~is_limit], yerr=rows['sigma'][~is_limit], fmt='o', color=INK,
                        markersize=4.5, markeredgecolor=SURFACE, markeredgewidth=0.9, elinewidth=1.1, capsize=2.5, zorder=4)
            ax.plot(x[is_limit], rows['measured'][is_limit], linestyle='none', marker='v', markersize=6, color=INK,
                    markeredgecolor=SURFACE, markeredgewidth=0.9, zorder=4)
            names = ['ring-opened' if limited else q for q, limited in zip(rows['quantity'], is_limit)]
            ax.set_xticks(x, names, fontsize=7)
            ax.set_xlim(-0.6, len(rows) - 0.4)
            ax.set_ylim(-0.07, 1.24)
            ax.set_yticks([0.0, 0.5, 1.0])
            ax.grid(axis='x', visible=False)
            ax.tick_params(length=0)
            if ax in axes[0]:
                note = f'spectrum {number or spectra} of {spectra}' if spectra > 1 else rows['averaged'].iloc[0]
                ax.set_title(f'{sample}\n{note}', fontsize=8.5, loc='left')
        row[0].set_ylabel(label, rotation=0, ha='right', va='center', fontsize=9, color=INK)
    handles = [Line2D([], [], marker='o', linestyle='', color=INK, markersize=4.5, label='measured share ± 1 standard error')]
    if any_limit:
        handles.append(Line2D([], [], marker='v', linestyle='', color=INK, markersize=6, label='published upper limit'))
    if predicted:
        handles += [Line2D([], [], marker='s', linestyle='', color=BLUE_RAMP[0], markersize=8, label=f'predicted, miss within {limit:g}'),
                    Line2D([], [], marker='s', linestyle='', color=SERIES[7], markersize=8, label=f'predicted, miss above {limit:g} (number: the miss)')]
    fig.legend(handles=handles, loc='outside lower center', ncols=len(handles))
    fig.suptitle(title or 'Area shares: predicted against measured', x=0.01, ha='left', fontsize=11)
    return axes


def plot_ladder(summary, *, order: list, limit: float = 3.0, title: str = None):
    """Worst miss (largest |z|) of the best fit of each structure, simplest first, one line per assumed age.

    Log axis. A structure fits an age where its worst miss lies in the shaded band, at or below the fit rule.
    """
    scenarios = [s for s in SCENARIO_ORDER if s in set(summary['scenario'])]
    fig, ax = plt.subplots(figsize=(8.8, 3.9), layout='constrained')
    ax.axhspan(1.0, limit, color=GRID, alpha=0.55, linewidth=0)
    ax.axhline(limit, color=INK_SOFT, linewidth=1.0)
    ax.text(-0.42, limit * 0.93, f'fits: worst miss ≤ {limit:g}', va='top', fontsize=8.5, color=INK_SOFT)
    labelled = set()
    for k, scenario in enumerate(scenarios):
        worst = summary[summary['scenario'] == scenario].set_index('structure')['max |z|'].reindex(order)
        x = np.arange(len(order)) + 0.13 * (k - (len(scenarios) - 1) / 2)
        ax.plot(x, worst, color=SCENARIO_COLORS[scenario], linewidth=1.3, marker='o', markersize=6.5, markeredgecolor=SURFACE,
                markeredgewidth=1.2, label=AGE_LABELS[scenario])
        for xi, name, value in zip(x, order, worst):
            # The near and full fits carry their number, once per structure and value
            if value < 2.0 * limit and (name, round(value, 1)) not in labelled:
                labelled.add((name, round(value, 1)))
                ax.annotate(f'{value:.2f}' if abs(value - limit) < 0.1 else f'{value:.1f}', (xi, value), xytext=(-8, -4),
                            textcoords='offset points', ha='right', va='top', fontsize=8, color=INK)
    free = summary.groupby('structure')['parameters'].min()
    ax.set_xticks(range(len(order)), [f"{name}\n{free[name]} parameter{'s' if free[name] > 1 else ''}" for name in order])
    ax.set_xlim(-0.5, len(order) - 0.5)
    ax.set_yscale('log')
    ax.set_ylim(1.0, 100.0)
    ax.set_yticks([1, 3, 10, 30, 100], ['1', '3', '10', '30', '100'])
    ax.minorticks_off()
    ax.set_ylabel('worst miss (standard errors)')
    ax.grid(axis='x', visible=False)
    ax.legend(title='assumed age', title_fontsize=8.5, loc='upper right', ncols=len(scenarios))
    fig.suptitle(title or 'The ladder: worst miss of each structure at each assumed age', x=0.01, ha='left', fontsize=11)
    return ax


def calculate_half_life_barrier(t_half_s, T_K: float):
    """Barrier [eV] whose Eyring rate gives the half-life t_half_s [s] with the reaction partner at 1 M."""
    return KB_EV * T_K * np.log(KB_SI * T_K / H_SI * np.asarray(t_half_s, dtype=float) / np.log(2.0))


def _draw_rows(ax, values, bands, scenarios, rows, labels):
    """Dot per scenario at values[(scenario, row)] with a bar over bands[(scenario, row)], one line of the axis per row."""
    for i, name in enumerate(rows):
        for k, scenario in enumerate(scenarios):
            y = i + 0.2 * (k - (len(scenarios) - 1) / 2)
            color = SCENARIO_COLORS[scenario]
            if bands is not None and (scenario, name) in bands:
                ax.plot(bands[(scenario, name)], [y, y], color=color, linewidth=3.5, alpha=0.45, solid_capstyle='butt')
            ax.plot(values[(scenario, name)], y, 'o', color=color, markersize=6.5, markeredgecolor=SURFACE, markeredgewidth=1.2, zorder=4)
    ax.set_yticks(range(len(rows)), [labels.get(name, name) for name in rows], fontsize=9, color=INK_SOFT)
    ax.set_ylim(len(rows) - 0.45, -0.55)
    ax.set_yticks(np.arange(0.5, len(rows) - 1), minor=True)
    ax.grid(axis='y', visible=False)
    ax.grid(axis='y', which='minor', color=GRID, linewidth=0.7)
    ax.tick_params(axis='y', which='both', length=0)


def plot_barrier_ruler(barriers, bands=None, *, structure: str, scenarios=('middle', 'long', 'free'),
                       reactions=('R1', 'R2', 'R3', 'R5', 'R4', 'R8'), T_C=(22.5, 80.0), x_range=(0.35, 1.75), title: str = None):
    """ΔG‡ of each reaction for the best fit of one structure at each assumed age, on a ruler of half-lives.

    barriers: 'structure', 'scenario', 'reaction', 'dG_barrier_eV'. bands: the same keys with 'dG‡ low' and 'dG‡ high',
    drawn as the range over the accepted parameter sets. The upper scales give the Eyring half-life that a barrier
    means at each temperature of T_C with the reaction partner at 1 M (ΔS‡ = 0).
    """
    fits = barriers[barriers['structure'] == structure].set_index(['scenario', 'reaction'])['dG_barrier_eV']
    ranges = None
    if bands is not None:
        chosen = bands[bands['structure'] == structure].set_index(['scenario', 'reaction'])
        ranges = {key: (r['dG‡ low'], r['dG‡ high']) for key, r in chosen.iterrows()}
    seen = barriers[barriers['structure'] == structure].groupby('reaction')['T_K'].first() - ZERO_CELSIUS_K
    labels = {name: f'{REACTION_LABELS.get(name, name)}\nseen at {seen[name]:.0f} °C' for name in reactions}
    fig, ax = plt.subplots(figsize=(9.2, 1.5 + 0.52 * len(reactions) + 0.42 * len(T_C)), layout='constrained')
    _draw_rows(ax, fits, ranges, scenarios, reactions, labels)
    ax.set_xlim(*x_range)
    ax.set_xlabel('ΔG‡ (eV)')
    for level, T in enumerate(T_C):
        ticks = {label: calculate_half_life_barrier(t, T + ZERO_CELSIUS_K) for label, t in HALF_LIVES_S.items()}
        ticks = {label: value for label, value in ticks.items() if x_range[0] < value < x_range[1]}
        scale = ax.secondary_xaxis(1.0 + 0.19 * level)
        scale.set_xticks(list(ticks.values()), list(ticks), fontsize=8)
        scale.set_xlabel(f'half-life at {T:g} °C, partner at 1 M', fontsize=8.5, loc='left')
        scale.tick_params(length=3, color=RULE)
        scale.spines['top'].set_color(RULE)
    handles = [Line2D([], [], marker='o', linestyle='', color=SCENARIO_COLORS[s], label=f'best fit, {AGE_LABELS[s]}') for s in scenarios]
    if ranges is not None:
        handles.append(Line2D([], [], color=INK_MUTED, linewidth=3.5, alpha=0.6, label='range over the accepted sets'))
    fig.legend(handles=handles, loc='outside lower center', ncols=len(handles))
    fig.suptitle(title or f'{structure}: barrier of each reaction at each assumed age', x=0.01, ha='left', fontsize=11)
    return ax


def plot_energy_shifts(barriers, bands=None, *, structure: str, computed: dict, sigma: dict, scenarios=('middle', 'long', 'free'),
                       reactions=('R1', 'R2', 'R3', 'R4'), title: str = None):
    """ΔG_rxn of the freed reactions: the computed value ± its standard error against the fitted value at each age.

    computed, sigma: {reaction: eV}. bands: 'dG_rxn low' and 'dG_rxn high' over the accepted parameter sets.
    """
    fits = barriers[barriers['structure'] == structure].set_index(['scenario', 'reaction'])['dG_rxn_eV']
    ranges = None
    if bands is not None:
        chosen = bands[bands['structure'] == structure].set_index(['scenario', 'reaction'])
        ranges = {key: (r['dG_rxn low'], r['dG_rxn high']) for key, r in chosen.iterrows()}
    fig, ax = plt.subplots(figsize=(8.2, 1.3 + 0.62 * len(reactions)), layout='constrained')
    ax.axvline(0.0, color=INK_SOFT, linewidth=1.0)
    _draw_rows(ax, fits, ranges, scenarios, reactions, REACTION_LABELS)
    for i, name in enumerate(reactions):
        y = i - 0.38
        ax.errorbar(computed[name], y, xerr=sigma[name], fmt='|', color=INK, markersize=11, markeredgewidth=1.6, ecolor=INK_MUTED,
                    elinewidth=1.2, capsize=2.5, zorder=3)
    ax.text(0.008, len(reactions) - 0.52, 'ΔG_rxn = 0: the reaction stops partway', fontsize=8, color=INK_SOFT, va='bottom')
    ax.set_xlabel('ΔG_rxn (eV); negative = downhill')
    handles = [Line2D([], [], marker='|', linestyle='', color=INK, markersize=11, markeredgewidth=1.6, label='computed ± 1 standard error')]
    handles += [Line2D([], [], marker='o', linestyle='', color=SCENARIO_COLORS[s], label=f'fitted, {AGE_LABELS[s]}') for s in scenarios]
    if ranges is not None:
        handles.append(Line2D([], [], color=INK_MUTED, linewidth=3.5, alpha=0.6, label='range over the accepted sets'))
    fig.legend(handles=handles, loc='outside lower center', ncols=len(handles) if len(handles) <= 4 else 3)
    fig.suptitle(title or f'{structure}: reaction energies, computed and fitted', x=0.01, ha='left', fontsize=11)
    return ax


def plot_trace_scan(traces, *, scenarios=('middle', 'long', 'free'), value: str = 'max |z|', limit: float = 3.0, title: str = None):
    """Worst miss of each fit when a trace of TMSOH is present at mixing, against the size of the trace.

    Line: the fit as it is, parameters unchanged. Open circle: the fit made again with the trace in place.
    """
    levels = sorted(traces['TMSOH at mixing (M)'].unique())
    positive = [level for level in levels if level > 0]
    none = min(positive) / 10.0                                  # where 'no trace' sits on the log axis
    fig, ax = plt.subplots(figsize=(7.6, 3.7), layout='constrained')
    ax.axhspan(1.0, limit, color=GRID, alpha=0.55, linewidth=0)
    ax.axhline(limit, color=INK_SOFT, linewidth=1.0)
    for scenario in scenarios:
        rows = traces[traces['scenario'] == scenario].sort_values('TMSOH at mixing (M)')
        x = rows['TMSOH at mixing (M)'].where(rows['TMSOH at mixing (M)'] > 0, none)
        as_is = rows['parameters'] == 'as fitted without it'
        color = SCENARIO_COLORS[scenario]
        ax.plot(x[as_is], rows.loc[as_is, value], color=color, linewidth=1.8, marker='o', markersize=5.5, markeredgecolor=SURFACE,
                markeredgewidth=1.0)
        ax.plot(x[~as_is], rows.loc[~as_is, value], linestyle='none', marker='o', markersize=8, markerfacecolor=SURFACE,
                markeredgecolor=color, markeredgewidth=1.8, zorder=4)
    units = (('mM', 1e-3), ('µM', 1e-6), ('nM', 1e-9))
    names = [next(f'{level / scale:g} {unit}' for unit, scale in units if level >= 0.999 * scale) for level in positive]
    ax.set_xscale('log')
    ax.set_xticks([none] + positive, ['none'] + names)
    ax.minorticks_off()
    ax.set_yscale('log')
    ax.set_ylim(1.0, 60.0)
    ax.set_yticks([1, 3, 10, 30], ['1', '3', '10', '30'])
    ax.yaxis.set_minor_locator(plt.NullLocator())
    ax.set_xlabel('TMSOH present at mixing')
    ax.set_ylabel('worst miss (standard errors)')
    handles = [Line2D([], [], color=SCENARIO_COLORS[s], linewidth=1.8, marker='o', markersize=5.5, label=f'{AGE_LABELS[s]}, fit as it is')
               for s in scenarios]
    handles.append(Line2D([], [], linestyle='none', marker='o', markersize=8, markerfacecolor=SURFACE, markeredgecolor=INK_SOFT,
                          markeredgewidth=1.8, label='fitted again with the trace'))
    fig.legend(handles=handles, loc='outside lower center', ncols=len(handles))
    fig.suptitle(title or 'A trace of TMSOH at mixing: does the fit survive?', x=0.01, ha='left', fontsize=11)
    return ax


def plot_storage_ranges(storage, seeded=None, *, structure: str, scenarios=('middle', 'long', 'free'), title: str = None):
    """Days until 90 % of the water is consumed in the stored electrolyte, for the fit at each assumed age.

    Dot: best fit; bar: range over the accepted parameter sets. seeded: the same prediction with TMSOH at mixing,
    for the fit as it is (triangle) and for the fit made again with the trace (open circle).
    """
    data = storage[storage['structure'] == structure].set_index('scenario')
    fig, ax = plt.subplots(figsize=(8.2, 1.4 + 0.62 * len(scenarios)), layout='constrained')
    for i, scenario in enumerate(scenarios):
        row, color = data.loc[scenario], SCENARIO_COLORS[scenario]
        ax.plot([row['band low (h)'] / 24.0, row['band high (h)'] / 24.0], [i, i], color=color, linewidth=3.5, alpha=0.45, solid_capstyle='butt')
        if seeded is not None:
            traces = seeded[(seeded['structure'] == structure) & (seeded['scenario'] == scenario)]
            ax.plot(traces['t90 water (h), fit as it is'] / 24.0, [i + 0.24] * len(traces), linestyle='none', marker='^', markersize=6.5,
                    color=color, markeredgecolor=SURFACE, markeredgewidth=1.0)
            ax.plot(traces['t90 water (h), fitted with the trace'] / 24.0, [i - 0.24] * len(traces), linestyle='none', marker='o',
                    markersize=7, markerfacecolor=SURFACE, markeredgecolor=color, markeredgewidth=1.6)
        ax.plot(row['t90 water (h), best fit'] / 24.0, i, 'o', color=color, markersize=7, markeredgecolor=SURFACE, markeredgewidth=1.2, zorder=4)
    ax.set_yticks(range(len(scenarios)), [AGE_LABELS[s] for s in scenarios], fontsize=9, color=INK_SOFT)
    ax.set_ylim(len(scenarios) - 0.5, -0.5)
    ax.set_xscale('log')
    ax.set_xlim(2.0, 400.0)
    ax.set_xticks([3, 7, 30, 90, 180, 365], ['3 days', '1 week', '1 month', '3 months', '6 months', '1 year'])
    ax.minorticks_off()
    ax.grid(axis='y', visible=False)
    ax.tick_params(axis='y', length=0)
    ax.set_xlabel('time to consume 90 % of the water at 25 °C')
    ax.set_ylabel('assumed age')
    handles = [Line2D([], [], marker='o', linestyle='', color=INK_SOFT, markersize=7, label='best fit'),
               Line2D([], [], color=INK_MUTED, linewidth=3.5, alpha=0.6, label='range over the accepted sets')]
    if seeded is not None:
        handles += [Line2D([], [], marker='^', linestyle='', color=INK_SOFT, markersize=6.5, label='fit as it is, with 1 µM or 1 mM TMSOH'),
                    Line2D([], [], marker='o', linestyle='', markersize=7, markerfacecolor=SURFACE, markeredgecolor=INK_SOFT,
                           markeredgewidth=1.6, label='fitted again with that trace')]
    fig.legend(handles=handles, loc='outside lower center', ncols=2)
    fig.suptitle(title or f'{structure}: storage prediction at each assumed age', x=0.01, ha='left', fontsize=11)
    return ax
