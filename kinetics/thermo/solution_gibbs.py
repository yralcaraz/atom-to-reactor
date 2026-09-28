# ==============================================================================
# BLOCK 3B: SOLUTION-PHASE GIBBS FREE ENERGY CALCULATOR
# ==============================================================================
try:
    from kinetics.thermo.gas_thermo import calculate_gas_thermo, EV_TO_KJ_MOL, _get_species_database
    from kinetics.thermo.standard_state import calculate_standard_state_shift
except ImportError:
    from calculate_gas_thermo import calculate_gas_thermo, EV_TO_KJ_MOL, _get_species_database
    from calculate_standard_state_shift import calculate_standard_state_shift

def calculate_solution_gibbs(species_name: str, T_K: float, species_db: dict = None, mode: str = 'qRRHO', **kwargs) -> dict:
    """
    Combines gas-phase Quasi-RRHO free energy (or benchmark DFT G), MACE/SMD solvation energy,
    and the standard-state concentration shift to yield G_sol,i(T).
    
    G_sol,i(T) = G°_gas,i(T) + ΔE_solv,i + ΔG°→*(T)
    """
    gas_thermo = calculate_gas_thermo(species_name, T_K, species_db=species_db, mode=mode, **kwargs)
    db = _get_species_database(species_db)
    dE_solv_eV = db[species_name].get('dE_solv_eV', 0.0)
    dE_solv_kJ_mol = dE_solv_eV * EV_TO_KJ_MOL
    
    # In benchmark mode, standard state shift is omitted to match tank_model.ipynb directly
    if mode in ('wb97mv', 'b3lyp_benchmark'):
        dG_std_shift_kJ_mol = 0.0
    else:
        dG_std_shift_kJ_mol = calculate_standard_state_shift(T_K, **kwargs)
    
    G_sol_kJ_mol = gas_thermo['G_gas_kJ_mol'] + dE_solv_kJ_mol + dG_std_shift_kJ_mol
    
    return {
        'G_sol_kJ_mol': G_sol_kJ_mol,
        'G_sol_eV': G_sol_kJ_mol / EV_TO_KJ_MOL,
        'H_sol_kJ_mol': gas_thermo['H_gas_kJ_mol'] + dE_solv_kJ_mol,
        'S_sol_J_mol_K': gas_thermo['S_gas_J_mol_K'] - (dG_std_shift_kJ_mol * 1000.0 / T_K)
    }
