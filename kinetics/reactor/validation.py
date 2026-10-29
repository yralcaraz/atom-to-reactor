"""Observations of Gogoi et al. 2024 as functions of a simulated state.

Observables: the quantity each paper observation measures, as a function of the concentrations (read by the
paper samples of a fit, kinetics.fitting.residuals). Water series: the initial concentrations of a sample of
TMSPA in EC/DEC with 0.5–5 vol% water (read by the lab observables, kinetics.data.observables).

Source: Y. Alcaraz Galván
"""

from kinetics.data.experimental import get_water_series
from kinetics.data.network import NETWORK_SPECIES
from kinetics.reactor.protocol import DEFAULT_DENSITY_G_ML
from kinetics.data.snapshot import calculate_molar_mass

PHOSPHATES = ('TMSPA', 'BMSPA', 'MMSPA', 'H3PO4')


def _phosphate_share(species: str):
    """Observable: share of all P held by one phosphate (what a ³¹P spectrum measures)."""
    return lambda C, i, c0: C[i[species]] / sum(C[i[sp]] for sp in PHOSPHATES)


# Observable name → (label, function of concentrations C [M] (final, or species × time), species index,
# initial concentrations)
OBSERVABLES = {
    'ring_opened_fraction': ('TMSOH ring-opened into TMSOEG + TMSOdiEG (fraction)',
                             lambda C, i, c0: (C[i['TMSOEG']] + C[i['TMSOdiEG']]) / c0['TMSOH']),
    'hmdso_si_fraction': ('Si of TMSOH ending in HMDSO (fraction)',
                          lambda C, i, c0: 2.0 * C[i['HMDSO']] / c0['TMSOH']),
    'tmspa_conversion': ('TMSPA converted (fraction)',
                         lambda C, i, c0: 1.0 - C[i['TMSPA']] / c0['TMSPA']),
    **{f'{sp.lower()}_share': (f'{sp} share of all P (³¹P)', _phosphate_share(sp)) for sp in PHOSPHATES},
}


# ------------------------------------------------------------------------------
# Water series (TMSPA + H2O in EC/DEC)
# ------------------------------------------------------------------------------
def calculate_water_series_c0(h2o_vol_pct: float, medium: dict = None) -> dict:
    """Initial concentrations [M] of a water-series sample: the H2O stock diluted by the TMSPA addition."""
    medium = medium or get_water_series()['medium']
    dilution = 1.0 - medium['tmspa_vol_frac']
    c0 = dict.fromkeys(NETWORK_SPECIES, 0.0)
    c0['H2O'] = h2o_vol_pct / 100.0 * 1000.0 * DEFAULT_DENSITY_G_ML['H2O'] / calculate_molar_mass('H2O') * dilution
    c0['TMSPA'] = medium['tmspa_vol_frac'] * 1000.0 * DEFAULT_DENSITY_G_ML['TMSPA'] / calculate_molar_mass('TMSPA')
    c0['EC'] = medium['EC_M']
    return c0
