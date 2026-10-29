"""The water series of Gogoi et al. 2024: the initial concentrations of a sample of TMSPA in EC/DEC with 0.5–5 vol%
water (read by the lab observables, kinetics.data.observables). The checks of the models against the paper's
control experiments are in the branch M5-full.

Source: Y. Alcaraz Galván
"""

from kinetics.data.experimental import get_water_series
from kinetics.data.network import NETWORK_SPECIES
from kinetics.reactor.protocol import DEFAULT_DENSITY_G_ML
from kinetics.data.snapshot import calculate_molar_mass


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
