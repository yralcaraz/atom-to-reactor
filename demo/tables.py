"""Display tables for the notebooks (formatting only; the numbers come from kinetics).

Tables are returned as DataFrames of formatted strings, so they render without optional pandas extras.
"""

import numpy as np
import pandas as pd

from kinetics.constants import EV_TO_KCAL_MOL


def format_columns(df: pd.DataFrame, formats: dict, na_rep: str = '—') -> pd.DataFrame:
    """Copy of df with the listed columns rendered as strings ({column: format spec})."""
    out = df.copy()
    for col, spec in formats.items():
        out[col] = [na_rep if v is None or (isinstance(v, float) and np.isnan(v)) else format(v, spec) for v in df[col]]
    return out


def format_species_table(species_db: dict, species: list, nmr_catalog: dict, experimental=None) -> pd.DataFrame:
    """Per species: G_gas, ΔE_solv ± standard error, computed ²⁹Si / ³¹P shifts (and measured ones when given), frequencies."""
    def computed(el, sp):
        sites = [d for d, _, labile in nmr_catalog.get(el, {}).get(sp, []) if not labile]
        return f'{sites[0]:.1f}' if sites else '—'

    def measured(el, sp):
        if experimental is None:
            return '—'
        row = experimental[(experimental['nucleus'] == el) & (experimental['species'] == sp)]
        if row.empty:
            return '—'
        lo, hi = row.iloc[0]['low_ppm'], row.iloc[0]['high_ppm']
        return f'{lo:.1f}' if lo == hi else f'{lo:.1f} … {hi:.1f}'

    rows = {}
    for sp in species:
        r = species_db[sp]
        rows[sp] = {
            'G_gas (Eh)': f"{r['G_gas_Eh']:.5f}",
            'ΔE_solv (eV)': f"{r['dE_solv_eV']:+.3f}",
            'SE ΔE_solv (eV)': f"{r['dE_solv_sigma_eV']:.2f}",
            'δ29Si calc (ppm)': computed('Si', sp),
            'δ29Si exp (ppm)': measured('Si', sp),
            'δ31P calc (ppm)': computed('P', sp),
            'δ31P exp (ppm)': measured('P', sp),
            'frequencies': 'yes' if r.get('frequencies_cm1') else 'no',
        }
    return pd.DataFrame.from_dict(rows, orient='index')


def format_thermo_table(thermo: pd.DataFrame) -> pd.DataFrame:
    """ΔG_rxn in eV and kcal/mol, standard error of the solvation term, whether the sign is resolved at ±1 SE, and K_eq."""
    df = pd.DataFrame({
        'class': thermo['class'],
        'equation': thermo['equation'],
        'ΔG_rxn (eV)': thermo['dG_rxn_eV'],
        'ΔG_rxn (kcal/mol)': thermo['dG_rxn_eV'] * EV_TO_KCAL_MOL,
        '±1 SE solv (kcal/mol)': thermo['sigma_solv_eV'] * EV_TO_KCAL_MOL,
        'sign resolved': np.where(thermo['dG_rxn_eV'].abs() > thermo['sigma_solv_eV'], 'yes', 'no'),
        'K_eq': thermo['K_eq'],
    })
    return format_columns(df, {'ΔG_rxn (eV)': '+.3f', 'ΔG_rxn (kcal/mol)': '+.2f',
                               '±1 SE solv (kcal/mol)': '.2f', 'K_eq': '.2e'})


def format_rate_comparison(rates_by_model: dict) -> pd.DataFrame:
    """ΔG‡_f [eV] and k_f of every reaction, side by side for several models."""
    parts = {}
    for name, r in rates_by_model.items():
        parts[(name, 'ΔG‡_f (eV)')] = [f'{v:.3f}' for v in r['dG_barrier_f_eV']]
        parts[(name, 'k_f')] = [f'{v:.2e}' for v in r['k_f']]
    return pd.DataFrame(parts, index=next(iter(rates_by_model.values())).index)


def format_scorecard(scorecard: pd.DataFrame, na_rep: str = 'not reached') -> pd.DataFrame:
    """Batch scorecard (kinetics.summarize_batch_runs) with 3 significant digits."""
    return format_columns(scorecard, {col: '.3g' for col in scorecard.columns}, na_rep=na_rep)


def format_family_barriers(family_barriers: pd.DataFrame, constraints: pd.DataFrame = None) -> pd.DataFrame:
    """Barrier per family and model (kinetics.microkinetics.models.tabulate_family_barriers) next to the experimental windows."""
    out = family_barriers.copy()
    for col in out.columns[1:]:
        out[col] = [f'{v:.2f}' for v in out[col]]
    if constraints is not None:
        def window(members):
            texts = []
            for rid in members.split(', '):
                if rid in constraints.index:
                    c = constraints.loc[rid]
                    if np.isnan(c['high_eV']):
                        rng = f"≥ {c['low_eV']:.3g}"
                    elif np.isnan(c['low_eV']):
                        rng = f"≤ {c['high_eV']:.3g}"
                    else:
                        rng = f"{c['low_eV']:.3g}–{c['high_eV']:.3g}"
                    texts.append(f"ΔG‡({rid}) {rng} eV at {c['T_K'] - 273.15:.0f} °C ({c['status']})")
            return '; '.join(texts) or 'none (no time-resolved data)'
        out['experimental window (Gogoi 2024)'] = [window(m) for m in out['reactions']]
    return out


def format_control_checks(checks: pd.DataFrame) -> pd.DataFrame:
    """One row per control experiment: observation, observed window and each model's prediction (✓ / ✗)."""
    rows = {}
    for exp_id, group in checks.groupby('experiment', sort=False):
        first = group.iloc[0]
        lo = '' if np.isnan(first['low']) else f"≥ {first['low']:.2g}"
        hi = '' if np.isnan(first['high']) else f"≤ {first['high']:.2g}"
        row = {'conditions': first['label'], 'observable': first['observable_label'],
               'observed window': ' and '.join(x for x in (lo, hi) if x) + f" ({first['status']})"}
        for _, r in group.iterrows():
            row[r['model']] = f"{r['predicted']:.3g} {'✓' if r['consistent'] else '✗'}"
        rows[exp_id] = row
    return pd.DataFrame.from_dict(rows, orient='index')


def format_share_table(shares: pd.DataFrame, windows: list, index=('sample', 'acquired_at')) -> pd.DataFrame:
    """Measured area shares (kinetics.data.tabulate_area_shares) as 'share ± SE', one row per spectrum."""
    value = shares.pivot_table(index=list(index), columns='window', values='share')
    sigma = shares.pivot_table(index=list(index), columns='window', values='sigma')
    return pd.DataFrame({w: [f'{v:+.3f} ± {s:.3f}' for v, s in zip(value[w], sigma[w])] for w in windows},
                        index=value.index)
