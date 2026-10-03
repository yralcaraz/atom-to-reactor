"""Physical constants and unit conversions shared by every block (CODATA 2018, exact SI where defined).

Names carry their unit (`_SI`, `_EV`) or read as `<FROM>_TO_<TO>`. Import from here; never redefine
a constant locally, since small differences between copies break the detailed-balance tests.

Classification: PUBLIC (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván
"""

import math

# Fundamental constants
R_SI = 8.314462618              # gas constant [J/(mol·K)]
KB_SI = 1.380649e-23            # Boltzmann constant [J/K]
KB_EV = 8.617333262e-5          # Boltzmann constant [eV/K]
H_SI = 6.62607015e-34           # Planck constant [J·s]
HBAR_SI = H_SI / (2.0 * math.pi)
C_CM_S = 2.99792458e10          # speed of light [cm/s]
NA = 6.02214076e23              # Avogadro constant [1/mol]

# Standard states
T_STD_K = 298.15                # temperature of the stored DFT Gibbs energies [K]
P_STD_PA = 1.0e5                # gas standard state, 1 bar [Pa]
C_STD_MOL_M3 = 1000.0           # solution standard state, 1 M [mol/m³]
ZERO_CELSIUS_K = 273.15

# Unit conversions
EV_TO_J = 1.602176634e-19
EV_TO_KJ_MOL = EV_TO_J * NA / 1000.0        # ≈ 96.4853
KJ_MOL_TO_EV = 1.0 / EV_TO_KJ_MOL
EV_TO_KCAL_MOL = 23.060547830619
HARTREE_TO_EV = 27.211386245988
AMU_TO_KG = 1.66053906660e-27
ANGSTROM_TO_M = 1.0e-10
