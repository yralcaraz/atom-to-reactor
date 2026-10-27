"""Figures of a stored fit against what was measured (notebook 05 rev08): measured spectra against simulated
ones, the shares of a tube side by side, and the concentrations of a tube in time with its spectra.

Colours come from demo/style.py: blue for the main fit, other hues (OVERLAY_COLORS) for fits drawn over it.

Source: Y. Alcaraz Galván
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from demo.reactor_plots import CONCENTRATION_PANELS, plot_concentration_panels
from demo.style import GRID, INK, INK_MUTED, INK_SOFT, SERIES, SURFACE
from kinetics.constants import ZERO_CELSIUS_K
from kinetics.fitting.residuals import SHARE_SPECIES
from kinetics.reactor.protocol import build_protocol_schedule

TMSPA_VOL_FRACTION = 0.05               # TMSPA is added as 5 vol% and dilutes the water of the stock (calculate_recipe_molarities)
TRACE_MM = 1e-3                          # below 1 µM a simulated concentration is drawn as 0 (integrator noise)
OVERLAY_COLORS = [SERIES[1], SERIES[2], SERIES[4], SERIES[6]]    # fits drawn over one spectrum: not blue (the main fit), not yellow (measured)


# ------------------------------------------------------------------------------
# A fit against what was measured: spectra and concentrations (kinetics.fitting.reporting)
# ------------------------------------------------------------------------------
WINDOW_LABELS = {'P-silyl': 'silyl on phosphate', 'TMSOEG': 'TMSOEG (opened EC)'}
SMOOTH_FWHM_PPM = {'31P': 0.4, '13C': 0.06}     # display broadening: far below the distance between two lines
RAW_FWHM_PPM = {'31P': 0.02, '13C': 0.004}      # the same for the 'paper' view: about the resolution of the spectra
SPECTRUM_VIEWS = ('paper', 'smooth', 'intensity')
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
                            clear_height: float = 6.0, overlay: bool = False, fit_label: str = 'Simulated from the fit',
                            view: str = 'paper'):
    """Measured spectra with the spectrum a fit predicts drawn over them, one row per spectrum.

    spectra: kinetics.fitting.reporting.load_measured_spectra. shares: tabulate_fit_shares of the fit, or
    {label: table} to compare several fits: side by side, or with overlay as lines of different colour over the
    same spectra. rows: [(sample, spectrum number, note)], top to bottom. fit_label: legend text of a single fit.
    Both curves are broadened for display (SMOOTH_FWHM_PPM) and scaled so that a line holding the whole area has
    height 1: the height of a line is close to its share. The model predicts areas only, so each simulated line
    takes the position and the width of the measured line (the centre of the window and the display width where
    the measured window shows no line higher than clear_height times the noise). view: 'paper' draws the spectra at
    their own resolution (RAW_FWHM_PPM, noise visible, as in a paper), each one scaled to its tallest point, 'smooth' broadens both curves for display
    (SMOOTH_FWHM_PPM), 'intensity' broadens as 'smooth' and draws every row of a column on one intensity scale
    (the area each spectrum had before it was divided by its windows), so the rows compare in height. The numbers over a window are
    its share: measured / simulated, the simulated one in the colour of the fit (measured only with overlay). Tubes whose axis is offset are moved onto a
    common axis.
    """
    if view not in SPECTRUM_VIEWS:
        raise ValueError(f'view must be one of {SPECTRUM_VIEWS}, not {view!r}')
    width = (RAW_FWHM_PPM if view == 'paper' else SMOOTH_FWHM_PPM)[nucleus]
    unit = 1.0645 * width                        # 1 / height of a unit-area Gaussian of that width
    cases = shares if isinstance(shares, dict) else {None: shares}
    groups = [list(cases.items())] if overlay else [[case] for case in cases.items()]
    colors = OVERLAY_COLORS if overlay else [SERIES[0]]
    fig, axes = plt.subplots(len(rows), len(groups), figsize=(9.6 if len(groups) == 1 else 1.4 + 4.1 * len(groups), 0.9 + 1.05 * len(rows)),
                             layout='constrained', sharex=True, squeeze=False)
    lows, highs = [], []
    for position, (column, group) in enumerate(zip(axes.T, groups)):
        tables = [table[table['nucleus'] == nucleus].set_index(['sample', 'spectrum', 'window']) for _, table in group]
        if view == 'intensity':                  # one scale for the column: the tallest measured line reaches 1
            tallest = max(unit * spectra[(sample, nucleus, number)]['area'] * _smooth(
                spectra[(sample, nucleus, number)]['ppm'], spectra[(sample, nucleus, number)]['intensity'], width).max()
                for sample, number, _ in rows)
        for ax, (sample, number, note) in zip(column, rows):
            spectrum = spectra[(sample, nucleus, number)]
            scale = 1.0
            if view == 'intensity':
                scale = spectrum['area'] / tallest
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
                if overlay:
                    ax.text(0.5 * (low + high), 1.24, f"{found[0]['measured']:.2f}", ha='center', va='center', fontsize=7.5, color=INK_SOFT)
                else:                            # measured / simulated, centred on the slash, simulated in the colour of the fit
                    ax.text(0.5 * (low + high), 1.24, f"{found[0]['measured']:.2f} /", ha='right', va='center', fontsize=7.5, color=INK_SOFT)
                    ax.text(0.5 * (low + high), 1.24, f" {found[0]['predicted']:.2f}", ha='left', va='center', fontsize=7.5, color=SERIES[0])
                if ax is column[0]:
                    ax.text(0.5 * (low + high), 1.5, WINDOW_LABELS.get(window, window), ha='center', va='bottom', fontsize=8.5, color=INK)
                lows.append(low)
                highs.append(high)
            if view == 'paper':                  # the tallest point of the measured spectrum reaches 1, as each spectrum is shown alone
                scale = 1.0 / (unit * measured.max())
                ax.plot(x, scale * unit * measured, color=SERIES[7], linewidth=0.7, zorder=3)     # red NMR over the blue fit, which follows it closely
            else:
                ax.fill_between(x, 0.0, scale * unit * measured, color=SERIES[3], alpha=0.6, linewidth=0)
            for curve, color in zip(simulated, colors):
                ax.plot(x, scale * unit * curve, color=color, linewidth=1.3 if overlay else 1.5)
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
        numbers = 'Numbers: measured share of the area'
    else:
        lines = [Line2D([], [], color=SERIES[0], linewidth=1.5, label=fit_label)]
        numbers = 'Numbers: share of the area, measured / simulated (in the colour of the fit)'
    # Three sections, filled column by column: measured | the fits, two per column | window and numbers.
    blank = Line2D([], [], linestyle='none', label=' ')
    fits = [lines[i:i + 2] + [blank] * (2 - len(lines[i:i + 2])) for i in range(0, len(lines), 2)]
    measured_key = (Line2D([], [], color=SERIES[7], linewidth=0.7, label='Measured spectrum') if view == 'paper'
                    else plt.Rectangle((0, 0), 1, 1, color=SERIES[3], alpha=0.6, linewidth=0, label='Measured spectrum'))
    handles = [measured_key, blank,
               *sum(fits, []),
               plt.Rectangle((0, 0), 1, 1, color=GRID, alpha=0.45, linewidth=0, label='Window of each species'),
               Line2D([], [], linestyle='none', label=numbers)]
    fig.legend(handles=handles, loc='outside lower center', ncols=len(fits) + 2)
    fig.suptitle(title or f'{NUCLEUS_LABELS[nucleus]} spectra: measured, and simulated from the fit', x=0.01, ha='left', fontsize=11)
    return axes


def tabulate_tube_comparison(shares: pd.DataFrame, sample: str) -> pd.DataFrame:
    """What the model gives and what the NMR gave, in mM, for every window of every spectrum of one tube.

    shares: kinetics.fitting.reporting.tabulate_fit_shares; 'spectrum at' is the time since TMSPA (or TMSOH) addition,
    the clock of the figures. A concentration is share × total_M (the amount the windows
    of that spectrum count together). 'Δy' is model − NMR; 'σ' is the error of the NMR share in mM; 'z (SD)' is Δy divided by σ.
    """
    rows = shares[shares['sample'] == sample].copy()
    scale = 1000.0 * rows['total_M']
    table = pd.DataFrame({
        'spectrum at (h)': rows['t_h'].round(2), 'nucleus': rows['nucleus'].map(NUCLEUS_LABELS),
        'species': rows['window'].map(lambda w: WINDOW_LABELS.get(w, w)),
        'NMR (mM)': scale * rows['measured'], 'model (mM)': scale * rows['predicted'],
        'Δy (mM)': scale * (rows['predicted'] - rows['measured']),
        'σ (mM)': scale * rows['sigma'],
        'z (SD)': (rows['predicted'] - rows['measured']) / rows['sigma']})
    return table.sort_values(['spectrum at (h)', 'nucleus'], kind='stable').reset_index(drop=True)


def show_tube_comparison(shares: pd.DataFrame, sample: str):
    """Displays tabulate_tube_comparison, the 'z (SD)' column coloured by its size: blue at 0, red from |z| = 1 on (above 1 or below −1)."""
    from IPython.display import display
    table = tabulate_tube_comparison(shares, sample)
    display(table.style.hide(axis='index')
            .format({'spectrum at (h)': '{:.1f}', 'NMR (mM)': '{:.1f}', 'model (mM)': '{:.1f}', 'Δy (mM)': '{:+.1f}', 'σ (mM)': '{:.2f}', 'z (SD)': '{:+.1f}'})
            .background_gradient(subset=['z (SD)'], cmap='coolwarm', gmap=table['z (SD)'].abs(), vmin=0, vmax=1)
            .set_caption(f'Tube {sample}: model against NMR at each spectrum'))


def plot_tube_concentrations(trajectories, sample: str, spectrum_times, *, model: str, beyond: float = 1.15,
                             before: str = None, margin_h: float = 24.0, solvent_stage: bool = False, measured=None):
    """Concentrations the model gives one tube, in the four panels of the reactor figures
    (demo.reactor_plots.plot_concentration_panels). The measured shares are drawn as points.

    trajectories: kinetics.fitting.reporting.simulate_fit_trajectories. spectrum_times: hours since mixing of the
    spectra of the tube. model: the age reading, written in the title after the tube. A
    heated tube gets its temperature above the panels. measured: kinetics.fitting.reporting.tabulate_fit_shares
    of the fit (all samples; the tube is picked out here). A point is the measured amount (share × total_M) at its
    spectrum time, in molecules: a silyl window counts silyl groups, so HMDSO (two per molecule) is halved. A window
    of two species (TMSOEG with TMSOdiEG) is their sum and is drawn grey; the P-silyl window has no curve and is not
    drawn, nor are the points of a panel that says 'none in this tube'. The error bar of a point is ± the error the spectrum
    itself gives (noise and baseline, 'sigma_nmr'), without the declared floor the fit adds; a point or bar below 0 lowers the axis to show it.
    The clock is the experimental one: t = 0 is the addition of TMSPA (the tubes with TMSPA) or of TMSOH (mixing).
    Each NMR event is a mark on the time axis, like solvent and TMSPA: the ³¹P and ¹³C scans taken within 1 h of each other
    are one event (their mean time). Nothing is greyed out.
    solvent_stage: the tube was prepared as the bench protocol says (kinetics.reactor.protocol.build_protocol_schedule):
    EC/DEC with water only, stirred and left to equilibrate before t = 0, with the water of the stock (before the 5 vol%
    TMSPA dilutes it). At t = 0 TMSPA jumps from 0 to its concentration and the water steps down. The fit does not
    simulate this stage; it is drawn as the constant it is (no TMSPA, nothing to react).
    before: whether the age is a reading (drawn with the model running on after the last spectrum); the only text under the panels is that. margin_h (or 15 % of the measured span if longer) is
    drawn after the last spectrum; where before is None, the curve runs to beyond × the last spectrum instead.
    Concentrations below 1 µM are drawn as 0 (what is left there is the integrator's noise, with either sign), every
    axis starts at 0 unless a measured point is below 0 (then the axis reaches it), and a panel with nothing above 1 µM says so.
    """
    rows = trajectories[trajectories['sample'] == sample]
    t0, t1 = min(spectrum_times), max(spectrum_times)
    if before is None:
        curve = rows[rows['t_h'] <= beyond * t1]
    else:
        curve = rows[rows['t_h'] <= t1 + max(margin_h, 0.15 * (t1 - t0))]
    species = [c for c in curve.columns if c not in ('sample', 't_h', 'T_C')]
    t, T_C, C_mM = curve['t_h'].to_numpy(), curve['T_C'].to_numpy(), 1000.0 * curve[species].to_numpy().T
    with_tmspa = C_mM[species.index('TMSPA')][0] > 0
    pre_h = stir_h = 0.0
    if solvent_stage:
        stages, _ = build_protocol_schedule()
        stir_h = stages[0].duration_h                                         # stirred, then equilibrated
        pre_h = sum(stage.duration_h for stage in stages[:2])
        stock = np.zeros((len(species), 1))                                   # TMSPA and TMSOH are added at t = 0, so they are absent here
        dilution = 1.0 - TMSPA_VOL_FRACTION if with_tmspa else 1.0           # only TMSPA dilutes the water of the stock
        stock[species.index('H2O')] = C_mM[species.index('H2O')][0] / dilution
        for solvent in ('EC', 'DEC'):                                         # the solvent is there before and after
            if solvent in species:
                stock[species.index(solvent)] = C_mM[species.index(solvent)][0]
        t = np.r_[-pre_h, 0.0, 0.0, t]
        T_C = np.r_[T_C[0], T_C[0], T_C[0], T_C]
        C_mM = np.hstack([stock, stock, C_mM[:, :1], C_mM])                   # the jump at t = 0 is a vertical step
    C_mM[np.abs(C_mM).max(axis=1) < TRACE_MM] = 0.0
    sim = {'model': f'{model} · tube {sample}', 't_h': t, 'T_K': T_C + ZERO_CELSIUS_K,
           'C_mM': np.clip(C_mM, 0.0, None), 'idx': {sp: i for i, sp in enumerate(species)},
           'acquisitions': [],                                  # the NMR lines are drawn below, one per event
           'time_label': 'Time from TMSPA addition (h)' if with_tmspa else 'Time since mixing (h)'}
    fig = plot_concentration_panels(sim, from_zero=True)
    fig.suptitle(f'Tube {sample}, {model}', x=0.01, ha='left', fontsize=11)
    for ax, (_, members) in zip(fig.axes, CONCENTRATION_PANELS):       # the four panels come first in fig.axes
        if not sim['C_mM'][[sim['idx'][sp] for sp in members]].any():
            for artist in [*ax.lines, *ax.texts, ax.get_legend()]:          # flat lines, their legend and end labels
                if artist is not None:
                    artist.remove()
            ax.set_yticks([])
            ax.text(0.5, 0.5, 'none in this tube (all below 1 µM)', transform=ax.transAxes, ha='center', va='center',
                    fontsize=9, color=INK_MUTED)
    points = measured[measured['sample'] == sample] if measured is not None else None
    drawn = {}                                                        # per panel: the legend labels of its points
    if points is not None:
        for ax, (_, members) in zip(fig.axes[:4], CONCENTRATION_PANELS):
            if not sim['C_mM'][[sim['idx'][sp] for sp in members]].any():       # a panel that says 'none in this tube' has no points
                continue
            for window in members:
                for nucleus, seen in points[points['window'] == window].groupby('nucleus'):
                    owners = SHARE_SPECIES[nucleus][window]                 # the species of the window, with their silyl groups
                    single = len(owners) == 1
                    drawn.setdefault(ax, set()).add('NMR (measured, ± noise and baseline)' if single else
                                                    f"NMR (measured): {' + '.join(owners)}, one {NUCLEUS_LABELS[nucleus]} window")
                    per_molecule = 1.0 / owners[next(iter(owners))] if single else 1.0
                    colour = next((line.get_color() for line in ax.get_lines() if line.get_label() in owners and single), INK_MUTED)
                    y = 1000.0 * seen['measured'] * seen['total_M'] * per_molecule
                    error = 1000.0 * seen['sigma_nmr'] * seen['total_M'] * per_molecule
                    ax.errorbar(seen['t_h'], y, yerr=error, fmt='none', ecolor=colour, elinewidth=1.0, capsize=3, clip_on=False, zorder=3)
                    ax.plot(seen['t_h'], y, 'o', color=colour, markersize=5, markeredgecolor=SURFACE, markeredgewidth=0.8,
                            clip_on=False, zorder=4)
                    bottom, top = ax.get_ylim()
                    ax.set_ylim(min(bottom, 1.08 * float((y - error).min())), max(top, 1.08 * float((y + error).max())))

    def label_x(ax, marks):
        """Adds named ticks {time: (name, time text)} to the x axis; the numbered ticks and the marks that would collide are dropped or joined."""
        lo, hi = ax.get_xlim()
        gap = 0.06 * (hi - lo)
        joined = {}
        for time, (name, when) in sorted(marks.items()):
            near = [k for k in joined if abs(k - time) <= gap]
            if near:
                (first, first_when) = joined[near[0]]                 # two marks on one spot: one label, both times
                joined[near[0]] = (first if first == name else f'{first} +\n{name}',
                                   f"({first_when.strip('()').replace(' h', '')}, {when.strip('()')})")
            else:
                joined[time] = (name, when)
        ticks = [tick for tick in ax.get_xticks() if lo <= tick <= hi and all(abs(tick - m) > 1.5 * gap for m in joined)]
        ax.set_xticks(ticks + list(joined))
        ax.set_xticklabels([f'{tick:g}' for tick in ticks] + [f'{name}\n{when}' for name, when in joined.values()])

    start, end = -pre_h, t[-1]
    nmr_times = []                                 # one event for the scans of both nuclei, which follow each other within the hour
    for time in sorted(float(x) for x in spectrum_times):
        if nmr_times and time - nmr_times[-1][0] <= 1.0:
            nmr_times[-1].append(time)
        else:
            nmr_times.append([time])
    nmr_times = [float(np.mean(group)) for group in nmr_times]
    for ax in fig.axes:
        ax.set_xlim(start, end)
    marks = {0.0: ('TMSPA' if with_tmspa else 'mixing', '(0 h)')}
    if solvent_stage:
        marks[start] = ('mixing', f'({start:.0f} h)')
        marks[start + stir_h] = ('equilibrate', f'({start + stir_h:.0f} h)')
    for t_nmr in nmr_times:
        marks[t_nmr] = ('NMR', f'({t_nmr:.0f} h)')
    for ax, (_, members) in zip(fig.axes[:4], CONCENTRATION_PANELS):
        label_x(ax, marks)
        if sim['C_mM'][[sim['idx'][sp] for sp in members]].any():          # a panel that says 'none in this tube' has no legend
            handles = [*ax.get_legend_handles_labels()[0]]
            for label in sorted(drawn.get(ax, ())):                     # the plain label first, then those of the summed windows
                handles.append(Line2D([], [], color=INK_MUTED, marker='o', linestyle='none', markersize=5, label=label))
            ax.legend(handles=handles, loc='center left' if len(handles) > 2 else 'best')
    note = 'The model runs on after the last spectrum.'
    fig.text(0.01, -0.005, note, ha='left', va='top', fontsize=9, color=INK_SOFT)
    return fig
