"""Block 3A — Standard-state shift from the 1 bar ideal gas to the 1 M solution.

Classification: PUBLIC (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván
"""

import numpy as np

from kinetics.constants import C_STD_MOL_M3, P_STD_PA, R_SI


def calculate_standard_state_shift(T_K: float, *, P_std_Pa: float = P_STD_PA,
                                   c_std_mol_m3: float = C_STD_MOL_M3) -> float:
    """ΔG°→* = RT ln(c* / c°_gas(T)) with c°_gas = P°/RT  [kJ/mol]; +7.96 kJ/mol at 298.15 K."""
    c_gas_mol_m3 = P_std_Pa / (R_SI * T_K)
    return R_SI * T_K * np.log(c_std_mol_m3 / c_gas_mol_m3) / 1000.0
