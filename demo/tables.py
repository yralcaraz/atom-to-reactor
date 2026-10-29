"""Display tables for the notebooks (formatting only; the numbers come from kinetics).

Tables are returned as DataFrames of formatted strings, so they render without optional pandas extras.

Source: Y. Alcaraz Galván
"""

import numpy as np
import pandas as pd


def format_interval(row, spec: str = '.3f') -> str:
    """One profile interval as text: 'a–b', '≥ a', '≤ b' or 'not determined'."""
    low, high = row['low'], row['high']
    if np.isfinite(low) and np.isfinite(high):
        return f'{low:{spec}}–{high:{spec}}'
    if np.isfinite(low):
        return f'≥ {low:{spec}}'
    if np.isfinite(high):
        return f'≤ {high:{spec}}'
    return 'not determined'


def format_interval_table(intervals: pd.DataFrame, structure: str, scenarios=('10 min', '1 h', '8 h'),
                          tag: str = '') -> pd.DataFrame:
    """95 % profile intervals of one structure: one row per parameter, one column per scenario.

    A side the data leave open within the profile grid is not given a number; 'not determined' means both are open.
    """
    data = intervals[(intervals['structure'] == structure) & (intervals['tag'].fillna('') == tag)]
    rows = {}
    for name in dict.fromkeys(data['parameter']):
        rows[name] = {}
        for scenario in scenarios:
            match = data[(data['parameter'] == name) & (data['scenario'] == scenario)]
            rows[name][scenario] = format_interval(match.iloc[0], '.3f') \
                if len(match) else '—'
    return pd.DataFrame.from_dict(rows, orient='index')
