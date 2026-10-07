"""Block 3B — Solution-phase Gibbs energy: G_sol = G°_gas + ΔE_solv + ΔG°→*.

In 'wb97mv' mode the standard-state shift is omitted (as in the reference notebooks of P. Broqvist); it cancels in
every reaction of the default network because each step conserves the number of molecules.

Source: Y. Alcaraz Galván
"""

from kinetics.constants import EV_TO_KJ_MOL
from kinetics.thermo.gas import calculate_gas_thermo
from kinetics.data.species import resolve_species_database
from kinetics.thermo.standard_state import calculate_standard_state_shift


def calculate_solution_gibbs(species: str, T_K: float, *, species_db: dict = None,
                             thermo_mode: str = 'wb97mv') -> dict:
    """G, H [kJ/mol] and S [J/(mol·K)] of one solvated species in EC at T_K."""
    gas = calculate_gas_thermo(species, T_K, species_db=species_db, thermo_mode=thermo_mode)
    dE_solv_eV = resolve_species_database(species_db)[species].get('dE_solv_eV')
    if dE_solv_eV is None:
        raise KeyError(f"No solvation energy for '{species}': it cannot enter a solution-phase network")
    dE_solv_kJ_mol = dE_solv_eV * EV_TO_KJ_MOL
    shift_kJ_mol = 0.0 if thermo_mode == 'wb97mv' else calculate_standard_state_shift(T_K)

    G_sol_kJ_mol = gas['G_gas_kJ_mol'] + dE_solv_kJ_mol + shift_kJ_mol
    return {
        'G_sol_kJ_mol': G_sol_kJ_mol,
        'G_sol_eV': G_sol_kJ_mol / EV_TO_KJ_MOL,
        'H_sol_kJ_mol': gas['H_gas_kJ_mol'] + dE_solv_kJ_mol,
        'S_sol_J_mol_K': gas['S_gas_J_mol_K'] - shift_kJ_mol * 1000.0 / T_K,
    }
