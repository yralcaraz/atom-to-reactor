"""Figures for the fits of notebook 05: residual maps, parameters across scenarios, profiles, influence of
each sample, recovery of known parameters, sample trajectories and the protocol band; and, for the step-by-step
version of the notebook, compositions against measurements, the ladder of structures, barriers with their
ranges, energy shifts, the trace scan, the storage ranges, the check of a tube that was not fitted and the
verdict of every test at every assumed age; and, for the notebook of the third fit, measured spectra against
simulated ones, concentrations in time, intervals read off profiles, the correlation between parameters, the
recovery as offsets from the truth and the error of a share.

Colour does one job per figure: a diverging blue–grey–red scale for signed residuals and shifts (grey = none),
one hue per age scenario in a fixed order, a single blue ramp for an ordered quantity (ΔS‡), red for a predicted
share or a test that breaks the fit rule, one hue per NMR window (WINDOW_COLORS) wherever species are told apart.

Source: Y. Alcaraz Galván
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D

from demo.reactor_plots import plot_concentration_panels
from demo.style import BLUE_RAMP, GRID, INK, INK_MUTED, INK_SOFT, RULE, SERIES, SURFACE
from kinetics.constants import ZERO_CELSIUS_K
from kinetics.data.observables import PHOSPHATE_WINDOWS

SCENARIO_ORDER = ('short', 'middle', 'long', 'free')
SCENARIO_LABELS = {'short': 'short (1 h)', 'middle': 'middle (1 d)', 'long': 'long (7 d)', 'free': 'free ages'}
AGE_LABELS = {'short': '1 h', 'middle': '1 day', 'long': '7 days', 'free': 'free ages'}     # a scenario named by its age
REACTION_LABELS = {'R1': 'R1 · first hydrolysis', 'R2': 'R2 · second hydrolysis', 'R3': 'R3 · third hydrolysis',
                   'R5': 'R5 · first transfer', 'R4': 'R4 · condensation', 'R8': 'R8 · solvent attack'}
SCENARIO_COLORS = dict(zip(SCENARIO_ORDER, SERIES[:4]))
_NEUTRAL = '#f0efec'
DIVERGING = LinearSegmentedColormap.from_list('blue_grey_red', [SERIES[0], _NEUTRAL, SERIES[7]])
OVERLAY_COLORS = [SERIES[1], SERIES[2], SERIES[4], SERIES[6]]    # fits drawn over one spectrum: not blue (the main fit), not yellow (measured)


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
    """TMSPA left at each acquisition of the bench protocol, for each activation entropy in the declared range.

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
    fig.suptitle(title or f"Bench protocol with the {structure} best fit: the band of the unknown activation entropy",
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
    A sample read on two nuclei shows its ³¹P windows first and its ¹³C windows after a rule.
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
            phosphate = rows['quantity'].isin(PHOSPHATE_WINDOWS).to_numpy()
            both = phosphate.any() and (~phosphate & ~is_limit).any()
            if both:
                ax.axvline(phosphate.sum() - 0.5, color=RULE, linewidth=0.9)
            names = ['ring-opened' if limited else q for q, limited in zip(rows['quantity'], is_limit)]
            ax.set_xticks(x, names, fontsize=7)
            ax.set_xlim(-0.6, len(rows) - 0.4)
            ax.set_ylim(-0.07, 1.24)
            ax.set_yticks([0.0, 0.5, 1.0])
            ax.grid(axis='x', visible=False)
            ax.tick_params(length=0)
            if ax in axes[0]:
                note = f'spectrum {number or spectra} of {spectra}' if spectra > 1 else rows['averaged'].iloc[0]
                nuclei = '' if is_limit.all() else '³¹P | ¹³C' if both else '³¹P' if phosphate.any() else '¹³C'
                ax.set_title(f'{sample}\n' + ' · '.join(part for part in (nuclei, note) if part), fontsize=8.5, loc='left')
        row[0].set_ylabel(label, rotation=0, ha='right', va='center', fontsize=9, color=INK)
    handles = [Line2D([], [], marker='o', linestyle='', color=INK, markersize=4.5, label='measured share ± 1 standard error')]
    if any_limit:
        handles.append(Line2D([], [], marker='v', linestyle='', color=INK, markersize=6, label='published upper limit'))
    if predicted:
        handles += [Line2D([], [], marker='s', linestyle='', color=BLUE_RAMP[0], markersize=8, label=f'predicted, miss within {limit:g}'),
                    Line2D([], [], marker='s', linestyle='', color=SERIES[7], markersize=8, label=f'predicted, miss above {limit:g} (number: the miss)')]
    fig.legend(handles=handles, loc='outside lower center', ncols=len(handles) if fig.get_figwidth() >= 8.0 else 2)
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


def plot_barrier_ranges(barriers, bands=None, *, structure: str, scenarios=('middle', 'long', 'free'),
                        reactions=('R1', 'R2', 'R3', 'R5', 'R4', 'R8'), x_range=(0.35, 1.75), title: str = None):
    """ΔG‡ of each reaction for the best fit of one structure at each assumed age, with its range.

    barriers: 'structure', 'scenario', 'reaction', 'T_K', 'dG_barrier_eV'. bands: the same keys with 'dG‡ low' and
    'dG‡ high', drawn as the range over the accepted parameter sets. Each barrier is given at the temperature where
    its reaction was seen, named under the reaction.
    """
    chosen = barriers[barriers['structure'] == structure]
    fits = chosen.set_index(['scenario', 'reaction'])['dG_barrier_eV']
    ranges = None
    if bands is not None:
        rows = bands[bands['structure'] == structure].set_index(['scenario', 'reaction'])
        ranges = {key: (r['dG‡ low'], r['dG‡ high']) for key, r in rows.iterrows()}
    seen = chosen.groupby('reaction')['T_K'].first() - ZERO_CELSIUS_K
    labels = {name: f'{REACTION_LABELS.get(name, name)}\nseen at {seen[name]:.0f} °C' for name in reactions}
    fig, ax = plt.subplots(figsize=(9.2, 1.3 + 0.55 * len(reactions)), layout='constrained')
    _draw_rows(ax, fits, ranges, scenarios, reactions, labels)
    ax.set_xlim(*x_range)
    ax.set_xlabel('ΔG‡ (eV)')
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


def plot_tube_check(measured, cases: dict, *, start=None, labels: dict = None, title: str = None):
    """Shares of one spectrum of a tube that was not fitted: measured, and predicted by the fit at each assumed age.

    measured: table indexed by window with 'share' and 'sigma'. cases: {scenario: predicted share per window}, one
    bar per scenario. start: the same table for an earlier spectrum of the tube, drawn as open circles.
    """
    labels = labels or AGE_LABELS
    windows = list(measured.index)
    x = np.arange(len(windows))
    width = 0.72 / len(cases)
    fig, ax = plt.subplots(figsize=(7.6, 3.5), layout='constrained')
    for k, (scenario, predicted) in enumerate(cases.items()):
        ax.bar(x + width * (k - (len(cases) - 1) / 2), [predicted[w] for w in windows], width=0.88 * width,
               color=SCENARIO_COLORS[scenario], linewidth=0, label=f'predicted, {labels[scenario]}')
    if start is not None:
        ax.plot(x - 0.43, start['share'].reindex(windows), linestyle='none', marker='o', markersize=6.5, markerfacecolor=SURFACE,
                markeredgecolor=INK, markeredgewidth=1.3, zorder=4, label='measured, first spectrum')
    ax.errorbar(x, measured['share'], yerr=measured['sigma'], fmt='o', color=INK, markersize=5.5, markeredgecolor=SURFACE,
                markeredgewidth=1.0, elinewidth=1.2, capsize=3, zorder=5, label='measured ± 1 standard error')
    ax.set_xticks(x, windows)
    ax.set_xlim(-0.6, len(windows) - 0.4)
    ax.set_ylim(-0.05, 1.0)
    ax.set_ylabel('share of the ³¹P area')
    ax.grid(axis='x', visible=False)
    ax.tick_params(axis='x', length=0)
    fig.legend(loc='outside lower center', ncols=3)
    fig.suptitle(title or 'A tube that was not fitted: predicted against measured', x=0.01, ha='left', fontsize=11)
    return ax


def plot_verdicts(worst, *, limit: float = 3.0, title: str = None):
    """The fit rule applied to every test (rows) at each assumed age (columns).

    worst: table of the worst miss |z| of each test, NaN where the test was not run. A cell is blue where the worst
    miss is within the rule and red, with ✗, where it is above. The number is the worst miss.
    """
    values = worst.to_numpy(dtype=float)
    failed = np.ma.masked_invalid(np.where(np.isnan(values), np.nan, values > limit))
    lines = max(str(label).count('\n') + 1 for label in worst.index)
    fig, ax = plt.subplots(figsize=(5.6 + 1.15 * worst.shape[1], 1.1 + (0.26 + 0.14 * lines) * worst.shape[0]), layout='constrained')
    ax.imshow(failed, cmap=ListedColormap([BLUE_RAMP[0], SERIES[7]]), vmin=0.0, vmax=1.0, aspect='auto', alpha=0.75)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            v = values[i, j]
            if np.isnan(v):
                _cell_text(ax, j, i, 'not run')
            else:
                _cell_text(ax, j, i, (f'{v:.1f}' if v < 9.95 else f'{v:.0f}') + (' ✗' if v > limit else ' ✓'), strong=v > limit)
    ax.set_xticks(range(worst.shape[1]), worst.columns, fontsize=9)
    ax.set_yticks(range(worst.shape[0]), worst.index, fontsize=8.5)
    ax.tick_params(length=0, top=True, labeltop=True, bottom=False, labelbottom=False)
    ax.set_xticks(np.arange(-0.5, worst.shape[1]), minor=True)
    ax.set_yticks(np.arange(-0.5, worst.shape[0]), minor=True)
    ax.grid(False)
    ax.grid(which='minor', color=SURFACE, linewidth=2)
    ax.tick_params(which='minor', length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    handles = [Line2D([], [], marker='s', linestyle='', color=BLUE_RAMP[0], alpha=0.75, markersize=9, label=f'✓ worst miss ≤ {limit:g}'),
               Line2D([], [], marker='s', linestyle='', color=SERIES[7], alpha=0.75, markersize=9, label=f'✗ worst miss above {limit:g}')]
    fig.legend(handles=handles, loc='outside lower center', ncols=2)
    fig.suptitle(title or 'Every test at every assumed age: worst miss', x=0.01, ha='left', fontsize=11)
    return ax


# ------------------------------------------------------------------------------
# A fit against what was measured: spectra and concentrations (kinetics.fitting.reporting)
# ------------------------------------------------------------------------------
WINDOW_COLORS = {'TMSPA': SERIES[0], 'BMSPA': SERIES[1], 'MMSPA': SERIES[2], 'H3PO4': SERIES[3],
                 'TMSOH': SERIES[4], 'HMDSO': SERIES[5], 'TMSOEG': SERIES[6], 'P-silyl': INK_MUTED}
WINDOW_LABELS = {'P-silyl': 'silyl on phosphate', 'TMSOEG': 'TMSOEG (opened EC)'}
SMOOTH_FWHM_PPM = {'31P': 0.4, '13C': 0.06}     # display broadening: far below the distance between two lines
NUCLEUS_LABELS = {'31P': '³¹P', '13C': '¹³C'}
_GAUSS = 2.0 * np.sqrt(2.0 * np.log(2.0))       # FWHM of a Gaussian in units of its standard deviation


def _smooth(x, y, fwhm: float):
    """y on the uniform axis x, convolved with a unit-area Gaussian of the given FWHM (areas are kept)."""
    sigma = fwhm / _GAUSS / abs(x[1] - x[0])
    half = int(np.ceil(4.0 * sigma))
    kernel = np.exp(-0.5 * (np.arange(-half, half + 1) / sigma) ** 2)
    return np.convolve(y, kernel / kernel.sum(), mode='same')


def _width_at_half_height(x, y, peak: int) -> float:
    """Full width at half height of the line whose maximum is at index `peak` (linear between points)."""
    half = 0.5 * y[peak]
    edges = []
    for step in (-1, 1):
        i = peak
        while 0 < i < len(y) - 1 and y[i] > half:
            i += step
        inner = i - step
        fraction = (y[inner] - half) / (y[inner] - y[i]) if y[inner] != y[i] else 0.0
        edges.append(x[inner] + fraction * (x[i] - x[inner]))
    return abs(edges[1] - edges[0])


def _line(x, centre: float, fwhm_gauss: float, fwhm_total: float):
    """Unit-area Voigt line: a Lorentzian broadened by the display Gaussian, with the total width given."""
    from scipy.special import voigt_profile
    a, b = 0.5346, 0.2166                        # Voigt width: a·L + sqrt(b·L² + G²)
    if fwhm_total <= fwhm_gauss:
        lorentz = 0.0
    else:
        root = np.sqrt(4.0 * a ** 2 * fwhm_total ** 2 - 4.0 * (a ** 2 - b) * (fwhm_total ** 2 - fwhm_gauss ** 2))
        lorentz = (2.0 * a * fwhm_total - root) / (2.0 * (a ** 2 - b))
    return voigt_profile(x - centre, fwhm_gauss / _GAUSS, 0.5 * lorentz)


def plot_spectra_comparison(spectra: dict, shares, rows: list, *, nucleus: str, x_range=None, title: str = None,
                            clear_height: float = 6.0, overlay: bool = False, fit_label: str = 'simulated from the fit'):
    """Measured spectra with the spectrum a fit predicts drawn over them, one row per spectrum.

    spectra: kinetics.fitting.reporting.load_measured_spectra. shares: tabulate_fit_shares of the fit, or
    {label: table} to compare several fits: side by side, or with overlay as lines of different colour over the
    same spectra. rows: [(sample, spectrum number, note)], top to bottom. fit_label: legend text of a single fit.
    Both curves are broadened for display (SMOOTH_FWHM_PPM) and scaled so that a line holding the whole area has
    height 1: the height of a line is close to its share. The model predicts areas only, so each simulated line
    takes the position and the width of the measured line (the centre of the window and the display width where
    the measured window shows no line higher than clear_height times the noise). The numbers over a window are
    its share: measured / simulated (measured only with overlay). Tubes whose axis is offset are moved onto a
    common axis.
    """
    width = SMOOTH_FWHM_PPM[nucleus]
    unit = 1.0645 * width                        # 1 / height of a unit-area Gaussian of that width
    cases = shares if isinstance(shares, dict) else {None: shares}
    groups = [list(cases.items())] if overlay else [[case] for case in cases.items()]
    colors = OVERLAY_COLORS if overlay else [SERIES[0]]
    fig, axes = plt.subplots(len(rows), len(groups), figsize=(9.6 if len(groups) == 1 else 1.4 + 4.1 * len(groups), 0.9 + 1.05 * len(rows)),
                             layout='constrained', sharex=True, squeeze=False)
    lows, highs = [], []
    for position, (column, group) in enumerate(zip(axes.T, groups)):
        tables = [table[table['nucleus'] == nucleus].set_index(['sample', 'spectrum', 'window']) for _, table in group]
        for ax, (sample, number, note) in zip(column, rows):
            spectrum = spectra[(sample, nucleus, number)]
            x = spectrum['ppm'] - spectrum['axis_offset_ppm']
            measured = _smooth(x, spectrum['intensity'], width)
            simulated = [np.zeros_like(x) for _ in tables]
            edges = [(low - spectrum['axis_offset_ppm'], high - spectrum['axis_offset_ppm']) for low, high in spectrum['windows_ppm'].values()]
            noise = np.std(measured[~np.any([(x >= low - 0.3) & (x <= high + 0.3) for low, high in edges], axis=0)])
            for window, (low, high) in zip(spectrum['windows_ppm'], edges):
                found = [table.loc[(sample, number, window)] for table in tables]
                inside = np.flatnonzero((x >= low) & (x <= high))
                centre, total = 0.5 * (low + high), width
                peak = inside[np.argmax(measured[inside])]
                if measured[peak] > clear_height * noise and inside[0] < peak < inside[-1]:
                    centre, total = x[peak], _width_at_half_height(x, measured, peak)
                for curve, row in zip(simulated, found):
                    curve += row['predicted'] * _line(x, centre, width, total)
                ax.axvspan(low, high, color=GRID, alpha=0.45, linewidth=0)
                number_text = f"{found[0]['measured']:.2f}" if overlay else f"{found[0]['measured']:.2f} / {found[0]['predicted']:.2f}"
                ax.text(0.5 * (low + high), 1.24, number_text, ha='center', va='center', fontsize=7.5, color=INK_SOFT)
                if ax is column[0]:
                    ax.text(0.5 * (low + high), 1.5, WINDOW_LABELS.get(window, window), ha='center', va='bottom', fontsize=8.5, color=INK)
                lows.append(low)
                highs.append(high)
            ax.fill_between(x, 0.0, unit * measured, color=SERIES[3], alpha=0.6, linewidth=0)
            for curve, color in zip(simulated, colors):
                ax.plot(x, unit * curve, color=color, linewidth=1.3 if overlay else 1.5)
            ax.set_ylim(-0.08, 1.42)
            ax.set_yticks([])
            ax.grid(False)
            if position == 0:
                ax.set_ylabel(f'{sample}\n{note}' if note else sample, rotation=0, ha='right', va='center', fontsize=8.5, color=INK)
        if not overlay and group[0][0] is not None:
            column[0].set_title(group[0][0], fontsize=9.5, loc='left', pad=20)
        column[-1].set_xlabel(f'δ({NUCLEUS_LABELS[nucleus]}) (ppm)')
    margin = 0.08 * (max(highs) - min(lows))
    axes.flat[-1].set_xlim(*(x_range or (max(highs) + margin, min(lows) - margin)))
    if overlay:
        lines = [Line2D([], [], color=color, linewidth=1.3, label=label) for (label, _), color in zip(groups[0], colors)]
        numbers = 'numbers: measured share of the area'
    else:
        lines = [Line2D([], [], color=SERIES[0], linewidth=1.5, label=fit_label)]
        numbers = 'numbers: share of the area, measured / simulated'
    handles = [plt.Rectangle((0, 0), 1, 1, color=SERIES[3], alpha=0.6, linewidth=0, label='measured spectrum'), *lines,
               plt.Rectangle((0, 0), 1, 1, color=GRID, alpha=0.45, linewidth=0, label='window of one line: its area gives the share'),
               Line2D([], [], linestyle='none', label=numbers)]
    fig.legend(handles=handles, loc='outside lower center', ncols=3 if overlay else 2)
    fig.suptitle(title or f'{NUCLEUS_LABELS[nucleus]} spectra: measured, and simulated from the fit', x=0.01, ha='left', fontsize=11)
    return axes


def plot_tube_concentrations(trajectories, sample: str, spectrum_times, *, model: str, beyond: float = 1.15):
    """Concentrations the model gives one tube since mixing, in the four panels of the reactor figures
    (demo.reactor_plots.plot_concentration_panels). No measurement is drawn.

    trajectories: kinetics.fitting.reporting.simulate_fit_trajectories. spectrum_times: hours since mixing of the
    spectra of the tube, drawn as dotted verticals; the time axis ends at beyond × the last one. model: the name
    written in the title. A heated tube gets its temperature above the panels.
    """
    curve = trajectories[(trajectories['sample'] == sample) & (trajectories['t_h'] <= beyond * max(spectrum_times))]
    species = [c for c in curve.columns if c not in ('sample', 't_h', 'T_C')]
    sim = {'model': f'{model} · tube {sample}', 't_h': curve['t_h'].to_numpy(), 'T_K': curve['T_C'].to_numpy() + ZERO_CELSIUS_K,
           'C_mM': 1000.0 * curve[species].to_numpy().T, 'idx': {sp: i for i, sp in enumerate(species)},
           'acquisitions': [{'t_h': t_h} for t_h in spectrum_times], 'time_label': 'Time since mixing (h)'}
    return plot_concentration_panels(sim)


def plot_concentrations(amounts: dict, shares, rows: dict, *, x_range=(0.1, 1500.0), title: str = None):
    """Amounts in each tube since mixing: the model as lines, the measured spectra as points.

    amounts: {nucleus: kinetics.fitting.reporting.tabulate_window_amounts of the simulated trajectories}.
    shares: tabulate_fit_shares of the same fit; a measured amount is its share × the nuclei the windows count.
    rows: {nucleus: [samples]}, one row of panels per nucleus. A point sits at the age the fit gives its tube plus
    the time of its spectrum. A grey band marks the time a tube spent above room temperature. Where the
    trajectories run on after the last spectrum of a tube, that part is dashed.
    """
    ncols = max(len(samples) for samples in rows.values())
    fig, axes = plt.subplots(len(rows), ncols, figsize=(3.05 * ncols, 0.9 + 2.55 * len(rows)), layout='constrained',
                             sharex=True, sharey='row', squeeze=False)
    used = []
    for row_axes, (nucleus, samples) in zip(axes, rows.items()):
        for ax, sample in zip(row_axes, samples):
            points = shares[(shares['sample'] == sample) & (shares['nucleus'] == nucleus)]
            curve = amounts[nucleus][amounts[nucleus]['sample'] == sample]
            heated = curve.loc[curve['T_C'] > 25.0, 't_h']
            if len(heated):
                start, stop = max(heated.min(), x_range[0]), heated.max()
                ax.axvspan(start, stop, color=GRID, alpha=0.6, linewidth=0)
                if stop > 3.0 * start:                               # wide enough on the log axis to carry its label
                    ax.text(0.03, 0.5, f"at {curve['T_C'].max():.0f} °C", transform=ax.transAxes, fontsize=8, color=INK_SOFT, va='center')
            last = points['t_h'].max()
            for window in dict.fromkeys(points['window']):
                color = WINDOW_COLORS[window]
                used.append((nucleus, window))
                before = curve['t_h'] <= last * (1.0 + 1e-9)
                ax.plot(curve.loc[before, 't_h'], 1000.0 * curve.loc[before, window], color=color, linewidth=1.8)
                ax.plot(curve.loc[~before, 't_h'], 1000.0 * curve.loc[~before, window], color=color, linewidth=1.3, linestyle=(0, (2, 2)))
                seen = points[points['window'] == window]
                ax.errorbar(seen['t_h'], 1000.0 * seen['measured'] * seen['total_M'], yerr=1000.0 * seen['sigma'] * seen['total_M'],
                            fmt='o', color=color, markersize=5.5, markeredgecolor=SURFACE, markeredgewidth=1.0, elinewidth=1.0,
                            capsize=2, zorder=4)
            ax.set_xscale('log')
            ax.set_xlim(*x_range)
            ax.set_title(sample, fontsize=9.5, loc='left')
        row_axes[0].set_ylabel({'31P': 'phosphate (mM)\nfrom ³¹P', '13C': 'silyl groups (mM)\nfrom ¹³C'}[nucleus])
        for ax in row_axes[len(samples):]:
            ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel('time since mixing (h)')
    by_nucleus = [[Line2D([], [], color=WINDOW_COLORS[w], linewidth=1.8, label=WINDOW_LABELS.get(w, w))
                   for n, w in dict.fromkeys(used) if n == nucleus] for nucleus in rows]
    kinds = [Line2D([], [], color=INK_SOFT, linewidth=1.8, label='line: model (dashed after the last spectrum)'),
             Line2D([], [], color=INK_SOFT, linestyle='none', marker='o', markersize=5.5, label='point: measured ± 1 standard error')]
    depth = max(len(group) for group in by_nucleus)
    blank = Line2D([], [], linestyle='none', label=' ')
    grid = [group + [blank] * (depth - len(group)) for group in by_nucleus][:2]
    handles = [h for pair in zip(*grid) for h in pair] if len(grid) == 2 else grid[0]
    fig.legend(handles=handles + kinds, loc='outside lower center', ncols=depth + 1)
    fig.suptitle(title or 'Each tube since mixing: model and measurement', x=0.01, ha='left', fontsize=11)
    return axes


# ------------------------------------------------------------------------------
# How well each fitted value is known
# ------------------------------------------------------------------------------
def plot_profile_cases(profile, cases: dict, *, delta: float = 3.84, y_max: float = 12.0, title: str = None):
    """How an interval is read off a profile, one panel per parameter of `cases` {parameter: panel title}.

    Line: the rise of χ² when the parameter is held at a value and every other one is fitted again. Rule: Δχ² = delta.
    Shaded: the values below the rule, which is the interval; it runs to the edge of the panel where the data do
    not close it.
    """
    fig, axes = plt.subplots(1, len(cases), figsize=(3.3 * len(cases), 3.0), layout='constrained', sharey=True, squeeze=False)
    for ax, (name, label) in zip(axes.flat, cases.items()):
        p = profile[profile['parameter'] == name].groupby('value', as_index=False)['chi2'].min().sort_values('value')
        rise = p['chi2'] - p['chi2'].min()
        accepted = p.loc[rise <= delta, 'value']
        ax.axvspan(accepted.min(), accepted.max(), color=GRID, alpha=0.7, linewidth=0)
        ax.plot(p['value'], rise, color=SERIES[0], linewidth=1.8, marker='o', markersize=4, markeredgecolor=SURFACE, markeredgewidth=0.6)
        ax.axhline(delta, color=INK_SOFT, linewidth=1.0)
        ax.set_ylim(-0.4, y_max)
        ax.set_xlim(p['value'].min(), p['value'].max())
        ax.set_title(label, fontsize=9.5, loc='left')
        ax.set_xlabel('value the parameter is held at (eV)')
    axes.flat[0].set_ylabel('rise of χ² above the best fit')
    handles = [Line2D([], [], color=SERIES[0], linewidth=1.8, marker='o', markersize=4, label='χ² with the other parameters fitted again'),
               Line2D([], [], color=INK_SOFT, linewidth=1.0, label=f'rise of {delta:g}: the end of a 95 % interval'),
               plt.Rectangle((0, 0), 1, 1, color=GRID, alpha=0.7, linewidth=0, label='interval')]
    fig.legend(handles=handles, loc='outside lower center', ncols=3)
    fig.suptitle(title or 'From a profile to an interval', x=0.01, ha='left', fontsize=11)
    return axes


def plot_parameter_intervals(parameters, intervals, *, structure: str, names: dict, boxes: dict, scenarios=('middle', 'long', 'free'),
                             xlabel: str = 'eV', title: str = None):
    """Best fit and 95 % interval of each parameter of `names` {parameter: row label}, one dot and bar per assumed age.

    intervals: 'structure', 'scenario', 'parameter', 'low', 'high'; a side the data do not close is NaN and is drawn
    to the edge of the search box with an arrow. boxes: {parameter: (lower, upper)}, the range the search was allowed,
    drawn as a pale track behind each row.
    """
    fits = parameters[parameters['structure'] == structure].set_index('scenario')
    found = intervals[intervals['structure'] == structure].set_index(['scenario', 'parameter'])
    fig, ax = plt.subplots(figsize=(8.6, 1.3 + 0.62 * len(names)), layout='constrained')
    for i, name in enumerate(names):
        low_box, high_box = boxes[name]
        ax.plot([low_box, high_box], [i, i], color=GRID, linewidth=17, alpha=0.5, solid_capstyle='butt', zorder=1)
        for k, scenario in enumerate(scenarios):
            y, color = i + 0.2 * (k - (len(scenarios) - 1) / 2), SCENARIO_COLORS[scenario]
            if (scenario, name) in found.index:
                row = found.loc[(scenario, name)]
                low, high = (row['low'] if np.isfinite(row['low']) else low_box), (row['high'] if np.isfinite(row['high']) else high_box)
                ax.plot([low, high], [y, y], color=color, linewidth=3.5, alpha=0.45, solid_capstyle='butt', zorder=2)
                for edge, is_open, marker in ((low, not np.isfinite(row['low']), '<'), (high, not np.isfinite(row['high']), '>')):
                    if is_open:
                        ax.plot(edge, y, marker=marker, color=color, markersize=6, zorder=3)
            ax.plot(fits.loc[scenario, name], y, 'o', color=color, markersize=6.5, markeredgecolor=SURFACE, markeredgewidth=1.2, zorder=4)
    ax.set_yticks(range(len(names)), list(names.values()), fontsize=9, color=INK_SOFT)
    ax.set_ylim(len(names) - 0.45, -0.55)
    ax.grid(axis='y', visible=False)
    ax.tick_params(axis='y', length=0)
    ax.set_xlabel(xlabel)
    handles = [Line2D([], [], marker='o', linestyle='', color=SCENARIO_COLORS[s], label=f'best fit, {AGE_LABELS[s]}') for s in scenarios]
    handles += [Line2D([], [], color=INK_MUTED, linewidth=3.5, alpha=0.6, label='95 % interval'),
                Line2D([], [], marker='>', linestyle='', color=INK_MUTED, markersize=6, label='open side: the data set no limit'),
                Line2D([], [], color=GRID, linewidth=8, alpha=0.6, label='range the search was allowed')]
    fig.legend(handles=handles, loc='outside lower center', ncols=3)
    fig.suptitle(title or f'{structure}: fitted values and their 95 % intervals', x=0.01, ha='left', fontsize=11)
    return ax


def plot_correlation(correlation, *, labels: dict = None, strong: float = 0.8, title: str = None):
    """Correlation between the fitted parameters near the best fit, lower triangle.

    Red: the two move together (raise one, and the other has to rise to keep the fit). Blue: one rises as the other
    falls. Grey: independent. The number is the correlation; bold from `strong` on.
    """
    names = list(correlation.index)
    values = correlation.to_numpy(dtype=float)
    shown = np.ma.masked_where(np.triu(np.ones_like(values, dtype=bool)), values)
    fig, ax = plt.subplots(figsize=(1.9 + 0.68 * len(names), 0.9 + 0.44 * len(names)), layout='constrained')
    image = ax.imshow(shown, cmap=DIVERGING, norm=TwoSlopeNorm(0.0, -1.0, 1.0), aspect='auto')
    for i in range(len(names)):
        for j in range(i):
            if np.isfinite(values[i, j]):
                _cell_text(ax, j, i, f'{values[i, j]:+.2f}', strong=abs(values[i, j]) >= strong)
    text = [(labels or {}).get(name, name) for name in names]
    ax.set_xticks(range(len(names) - 1), text[:-1], fontsize=8.5, rotation=35, ha='right', rotation_mode='anchor')
    ax.set_yticks(range(1, len(names)), text[1:], fontsize=8.5)
    ax.set_xlim(-0.5, len(names) - 1.5)
    ax.set_ylim(len(names) - 0.5, 0.5)
    ax.tick_params(length=0)
    ax.set_xticks(np.arange(-0.5, len(names)), minor=True)
    ax.set_yticks(np.arange(-0.5, len(names)), minor=True)
    ax.grid(False)
    ax.grid(which='minor', color=SURFACE, linewidth=2)
    ax.tick_params(which='minor', length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    bar = fig.colorbar(image, ax=ax, shrink=0.6, pad=0.02)
    bar.set_label('correlation', color=INK_SOFT, fontsize=8.5)
    bar.outline.set_visible(False)
    fig.suptitle(title or 'Do the fitted parameters depend on each other?', x=0.01, ha='left', fontsize=11)
    return ax


def plot_recovery_offsets(summary, *, names: dict, cases: dict = None, clip: float = 0.35, title: str = None):
    """Synthetic recovery in one picture: fitted value minus the value that made the data, per parameter.

    summary: recovery table ('case', 'parameter', 'truth', 'fit', 'low', 'high', 'truth inside'). names:
    {parameter: row label}, all in the same unit. cases: {case letter: legend text}. One dot and bar (95 % interval)
    per invented data set, coloured by case; the line at zero is the true value. A bar that does not reach zero
    missed the truth and carries a cross. An arrow is a side the data left open; offsets beyond ±clip are drawn at
    the edge.
    """
    data = summary[summary['parameter'].isin(names)]
    letters = list(dict.fromkeys(data['case'].str.split(',').str[0]))
    colors = dict(zip(letters, SERIES))
    runs = list(dict.fromkeys(data['case']))
    fig, ax = plt.subplots(figsize=(8.6, 1.2 + 0.78 * len(names)), layout='constrained')
    ax.axvline(0.0, color=INK, linewidth=1.2, zorder=1)
    for i, name in enumerate(names):
        rows = data[data['parameter'] == name]
        for _, r in rows.iterrows():
            y = i + 0.11 * (runs.index(r['case']) - (len(runs) - 1) / 2)
            color = colors[r['case'].split(',')[0]]
            low = np.clip(r['low'] - r['truth'], -clip, clip) if np.isfinite(r['low']) else -clip
            high = np.clip(r['high'] - r['truth'], -clip, clip) if np.isfinite(r['high']) else clip
            ax.plot([low, high], [y, y], color=color, linewidth=2.6, alpha=0.45, solid_capstyle='butt', zorder=2)
            for edge, is_open, marker in ((low, not np.isfinite(r['low']), '<'), (high, not np.isfinite(r['high']), '>')):
                if is_open:
                    ax.plot(edge, y, marker=marker, color=color, markersize=5, zorder=3)
            ax.plot(np.clip(r['fit'] - r['truth'], -clip, clip), y, 'o', color=color, markersize=5.5, markeredgecolor=SURFACE,
                    markeredgewidth=1.0, zorder=4)
            if not r['truth inside']:
                ax.plot(high + 0.012, y, marker='x', color=INK, markersize=6, markeredgewidth=1.6, zorder=5)
    ax.set_yticks(range(len(names)), list(names.values()), fontsize=9, color=INK_SOFT)
    ax.set_ylim(len(names) - 0.4, -0.6)
    ax.set_yticks(np.arange(0.5, len(names) - 1), minor=True)
    ax.grid(axis='y', visible=False)
    ax.grid(axis='y', which='minor', color=GRID, linewidth=0.7)
    ax.tick_params(axis='y', which='both', length=0)
    ax.set_xlim(-clip - 0.03, clip + 0.03)
    ax.set_xlabel('fitted value − true value (eV)')
    handles = [Line2D([], [], marker='o', linestyle='', color=colors[c], label=(cases or {}).get(c, f'case {c}')) for c in letters]
    handles += [Line2D([], [], color=INK_MUTED, linewidth=2.6, alpha=0.6, label='95 % interval of one invented data set (dot: its best fit)'),
                Line2D([], [], marker='>', linestyle='', color=INK_MUTED, markersize=5, label='open side: no limit found'),
                Line2D([], [], marker='x', linestyle='', color=INK, markersize=6, markeredgewidth=1.6, label='the interval misses the true value')]
    fig.legend(handles=handles, loc='outside lower center', ncols=2)
    fig.suptitle(title or 'Invented data with a known answer: does the fit give it back?', x=0.01, ha='left', fontsize=11)
    return ax


# ------------------------------------------------------------------------------
# The error of a share: what is declared and what is checked
# ------------------------------------------------------------------------------
def plot_error_parts(shares, *, parts: dict, title: str = None):
    """Size of each part of the error of a share, one dot per measured share, one panel per nucleus.

    shares: one row per measured share with 'nucleus' and one column per part. parts: {column: label}. The tick is
    the median of the column.
    """
    nuclei = [n for n in NUCLEUS_LABELS if n in set(shares['nucleus'])]
    fig, axes = plt.subplots(1, len(nuclei), figsize=(4.3 * len(nuclei), 3.2), layout='constrained', sharey=True, squeeze=False)
    rng = np.random.default_rng(0)
    for ax, nucleus in zip(axes.flat, nuclei):
        rows = shares[shares['nucleus'] == nucleus]
        for k, column in enumerate(parts):
            ax.plot(k + rng.uniform(-0.2, 0.2, len(rows)), rows[column], linestyle='none', marker='o', markersize=4.5, color=SERIES[0],
                    alpha=0.5, markeredgewidth=0)
            ax.plot([k - 0.3, k + 0.3], [rows[column].median()] * 2, color=INK, linewidth=1.6)
        ax.set_xticks(range(len(parts)), list(parts.values()), fontsize=8.5)
        ax.set_xlim(-0.6, len(parts) - 0.4)
        ax.set_title(f'{NUCLEUS_LABELS[nucleus]}: {len(rows)} measured shares', fontsize=9.5, loc='left')
        ax.grid(axis='x', visible=False)
        ax.tick_params(axis='x', length=0)
    axes.flat[0].set_ylabel('error of a share')
    axes.flat[0].set_ylim(0.0, None)
    handles = [Line2D([], [], marker='o', linestyle='', color=SERIES[0], alpha=0.5, markersize=4.5, label='one measured share'),
               Line2D([], [], color=INK, linewidth=1.6, label='median')]
    fig.legend(handles=handles, loc='outside lower center', ncols=2)
    fig.suptitle(title or 'What the error of a share is made of', x=0.01, ha='left', fontsize=11)
    return axes


def plot_replicates(shares, *, sample: str, scatter: dict = None, title: str = None):
    """Repeated spectra of one tube: each share minus the mean of its window, with its noise as the bar.

    shares: kinetics.data.tabulate_observable_shares. One panel per nucleus, one group of dots per window, in the
    order the spectra were recorded. scatter: {nucleus: scatter of the shares / their noise}, written in the title.
    """
    data = shares[(shares['sample'] == sample) & (shares['temperature_C'] < 25.0)]
    nuclei = [n for n in NUCLEUS_LABELS if n in set(data['nucleus'])]
    fig, axes = plt.subplots(1, len(nuclei), figsize=(4.6 * len(nuclei), 3.3), layout='constrained', sharey=True, squeeze=False)
    for ax, nucleus in zip(axes.flat, nuclei):
        rows = data[data['nucleus'] == nucleus]
        windows = list(dict.fromkeys(rows['window']))
        ax.axhline(0.0, color=INK_SOFT, linewidth=1.0)
        for k, window in enumerate(windows):
            group = rows[rows['window'] == window].sort_values('spectrum')
            weight = 1.0 / group['sigma_noise'] ** 2
            mean = (weight * group['share']).sum() / weight.sum()
            x = k + np.linspace(-0.3, 0.3, len(group))
            ax.errorbar(x, group['share'] - mean, yerr=group['sigma_noise'], fmt='o', color=WINDOW_COLORS[window], markersize=5,
                        markeredgecolor=SURFACE, markeredgewidth=1.0, elinewidth=1.1, capsize=2)
        ax.set_xticks(range(len(windows)), [WINDOW_LABELS.get(w, w).replace(' (', '\n(').replace(' on ', ' on\n') for w in windows], fontsize=8)
        ax.set_xlim(-0.6, len(windows) - 0.4)
        ax.grid(axis='x', visible=False)
        ax.tick_params(axis='x', length=0)
        note = f': scatter = {scatter[nucleus]:.1f} × noise' if scatter else ''
        ax.set_title(f'{NUCLEUS_LABELS[nucleus]}{note}', fontsize=9.5, loc='left')
    axes.flat[0].set_ylabel('share − mean of the repeats')
    fig.suptitle(title or f'{sample}: the same tube measured again and again', x=0.01, ha='left', fontsize=11)
    return axes


def plot_miss_distribution(residuals, *, limit: float = 3.0, bin_width: float = 0.5, title: str = None):
    """Histogram of the misses z of one fit against what errors of exactly the stated size would give.

    residuals: the residual table of the fit ('z'). Curve: a standard normal for the same number of values.
    """
    z = residuals['z'].to_numpy(dtype=float)
    edges = np.arange(-limit - 1.0, limit + 1.0 + bin_width, bin_width)
    fig, ax = plt.subplots(figsize=(6.4, 3.3), layout='constrained')
    ax.hist(z, bins=edges, color=BLUE_RAMP[0], edgecolor=SURFACE, linewidth=1.5)
    grid = np.linspace(edges[0], edges[-1], 300)
    ax.plot(grid, len(z) * bin_width * np.exp(-0.5 * grid ** 2) / np.sqrt(2.0 * np.pi), color=INK_SOFT, linewidth=1.6)
    for edge in (-limit, limit):
        ax.axvline(edge, color=INK, linewidth=1.0)
    ax.text(limit - 0.08, ax.get_ylim()[1] * 0.96, 'fit rule', ha='right', va='top', fontsize=8.5, color=INK_SOFT)
    ax.set_xlabel('miss = (predicted − measured) / error')
    ax.set_ylabel('number of values')
    ax.grid(axis='x', visible=False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=BLUE_RAMP[0], label=f'misses of the fit ({len(z)} values)'),
               Line2D([], [], color=INK_SOFT, linewidth=1.6, label='expected if the errors are exactly as stated')]
    fig.legend(handles=handles, loc='outside lower center', ncols=2)
    fig.suptitle(title or 'Are the misses as large as the errors say?', x=0.01, ha='left', fontsize=11)
    return ax


def plot_reading_trajectories(curves: dict, points: dict, *, samples: list, window: str = 'TMSPA', nucleus: str = '31P',
                              traced: dict = None, trace_label: str = None, x_range=(0.1, 2000.0), title: str = None):
    """Share of one window in a few tubes against the time since mixing, one panel per assumed age.

    curves: {scenario: kinetics.fitting.reporting.tabulate_window_amounts of the trajectories of that fit}.
    points: {scenario: tabulate_fit_shares of that fit}; a marker is a measured share at the age the scenario gives
    its tube. traced: the same as curves for the same parameters with something added at mixing, drawn dashed.
    """
    scenarios = [s for s in SCENARIO_ORDER if s in curves]
    colors = dict(zip(samples, SERIES))
    fig, axes = plt.subplots(1, len(scenarios), figsize=(3.7 * len(scenarios), 3.5), layout='constrained', sharey=True, squeeze=False)
    for ax, scenario in zip(axes.flat, scenarios):
        for sample in samples:
            for table, dashes, weight in ((curves, '-', 1.8), (traced, (0, (4, 2)), 1.5)):
                if table is None:
                    continue
                curve = table[scenario][table[scenario]['sample'] == sample]
                total = curve.drop(columns=['sample', 't_h', 'T_C']).sum(axis=1)
                ax.plot(curve['t_h'], curve[window] / total, color=colors[sample], linewidth=weight, linestyle=dashes)
            seen = points[scenario][(points[scenario]['sample'] == sample) & (points[scenario]['nucleus'] == nucleus)
                                    & (points[scenario]['window'] == window)]
            ax.plot(seen['t_h'], seen['measured'], linestyle='none', marker='o', markersize=6.5, markerfacecolor=colors[sample],
                    markeredgecolor=SURFACE, markeredgewidth=1.2, zorder=5)
        ax.set_xscale('log')
        ax.set_xlim(*x_range)
        ax.set_ylim(-0.05, 1.05)
        ax.set_title(f'assumed age: {AGE_LABELS[scenario]}' if scenario != 'free' else 'free ages', fontsize=9.5, loc='left')
        ax.set_xlabel('time since mixing (h)')
    axes.flat[0].set_ylabel(f'{window} share of the {NUCLEUS_LABELS[nucleus]} area')
    handles = [Line2D([], [], color=colors[sample], linewidth=1.8, label=sample) for sample in samples]
    handles.append(Line2D([], [], color=INK_SOFT, linewidth=1.8, label='model, from the recipe alone'))
    if traced is not None:
        handles.append(Line2D([], [], color=INK_SOFT, linewidth=1.5, linestyle=(0, (4, 2)), label=trace_label or 'model, with a trace at mixing'))
    handles.append(Line2D([], [], color=INK_SOFT, linestyle='none', marker='o', markersize=6.5, label='measured, placed at the assumed age'))
    fig.legend(handles=handles, loc='outside lower center', ncols=3)
    fig.suptitle(title or f'{window} since mixing under each assumed age', x=0.01, ha='left', fontsize=11)
    return axes
