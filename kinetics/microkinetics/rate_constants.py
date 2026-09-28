# ==============================================================================
# BLOCK 5: KINETIC RATE CONSTANTS ENGINE (WITH STRICT DETAILED BALANCE)
# ==============================================================================
import inspect
import numpy as np
try:
    from kinetics.thermo.reaction_thermo import calculate_reaction_thermo
    from kinetics.thermo.gas_thermo import R_SI, KB_SI, H_SI, EV_TO_KJ_MOL
except ImportError:
    from calculate_reaction_thermo import calculate_reaction_thermo
    from calculate_gas_thermo import R_SI, KB_SI, H_SI, EV_TO_KJ_MOL

DEFAULT_FAMILY_BEP_PARAMETERS = {
    'hydrolysis':     {'E0_eV': 0.80, 'alpha': 0.50},
    'condensation':   {'E0_eV': 0.80, 'alpha': 0.50},
    'transfer':       {'E0_eV': 0.80, 'alpha': 0.50},
    'solvent_attack': {'E0_eV': 1.30, 'alpha': 0.50},  # Finding 8: Gogoi et al. (Nat. Commun. 2024)
    'default':        {'E0_eV': 0.80, 'alpha': 0.50}
}

def _get_bep_parameters(bep_params=None):
    """Retrieves BEP parameters from arguments, globals, calling frame, or default family parameters."""
    if bep_params is not None:
        return bep_params
    if 'family_bep_parameters' in globals():
        return globals()['family_bep_parameters']
    frame = inspect.currentframe().f_back
    while frame:
        if 'family_bep_parameters' in frame.f_globals:
            return frame.f_globals['family_bep_parameters']
        if 'family_bep_parameters' in frame.f_locals:
            return frame.f_locals['family_bep_parameters']
        frame = frame.f_back
    return DEFAULT_FAMILY_BEP_PARAMETERS.copy()

# ------------------------------------------------------------------------------
# MODEL 1: Bell-Evans-Polanyi (BEP) Linear Scaling + Eyring TST
# ------------------------------------------------------------------------------
def bep_eyring(
    rxn_id: str,
    T_K: float,
    bep_params: dict = None,
    reactions_net: dict = None,
    species_db: dict = None,
    mode: str = 'qRRHO',
    **kwargs
) -> dict:
    """
    Computes k_f and k_r using Bell-Evans-Polanyi linear barrier scaling + Eyring TST:
        ΔG‡_f = max(E0, E0 + α · ΔG_rxn)
        k_f = (kB · T / h) · exp(-ΔG‡_f / RT)
        k_r = k_f / K_eq
    """
    thermo = calculate_reaction_thermo(
        rxn_id, T_K,
        reactions_net=reactions_net,
        species_db=species_db,
        mode=mode,
        **kwargs
    )
    rxn_class = thermo['class']
    all_bep = _get_bep_parameters(bep_params)
    default_fallback = DEFAULT_FAMILY_BEP_PARAMETERS.get(rxn_class, DEFAULT_FAMILY_BEP_PARAMETERS['default'])
    bep = all_bep.get(rxn_class, all_bep.get('default', default_fallback))
    
    E0_kJ_mol = bep['E0_eV'] * EV_TO_KJ_MOL
    alpha = bep['alpha']
    dG_rxn_kJ_mol = thermo['dG_rxn_kJ_mol']
    
    # 1. Forward Barrier: ΔG‡_f = max(E0, E0 + α * ΔG_rxn)
    dG_barrier_f_kJ_mol = max(E0_kJ_mol, E0_kJ_mol + alpha * dG_rxn_kJ_mol)
    
    # 2. Eyring Forward Rate: k_f = (kB*T / h) * exp(-ΔG‡ / RT)
    eyring_prefactor = (KB_SI * T_K) / H_SI
    k_f = eyring_prefactor * np.exp(-(dG_barrier_f_kJ_mol * 1000.0) / (R_SI * T_K))
    
    # 3. Detailed Balance: k_r = k_f / K_eq
    K_eq = thermo['K_eq']
    k_r = k_f / K_eq if K_eq > 1e-300 else 0.0
    dG_barrier_r_kJ_mol = dG_barrier_f_kJ_mol - dG_rxn_kJ_mol
    
    return {
        'rxn_id': rxn_id,
        'T_K': T_K,
        'kinetic_model': 'bep_eyring',
        'thermo_mode': mode,
        'class': rxn_class,
        'dG_rxn_eV': thermo['dG_rxn_eV'],
        'dG_rxn_kJ_mol': dG_rxn_kJ_mol,
        'dG_barrier_f_kJ_mol': dG_barrier_f_kJ_mol,
        'dG_barrier_f_eV': dG_barrier_f_kJ_mol / EV_TO_KJ_MOL,
        'dG_barrier_r_kJ_mol': dG_barrier_r_kJ_mol,
        'dG_barrier_r_eV': dG_barrier_r_kJ_mol / EV_TO_KJ_MOL,
        'k_f': k_f,
        'k_r': k_r,
        'K_eq': K_eq
    }

# ------------------------------------------------------------------------------
# MODEL 2: Marcus Theory Quadratic Activation + Eyring TST
# ------------------------------------------------------------------------------
def marcus_eyring(
    rxn_id: str,
    T_K: float,
    lambda_eV: float = None,
    bep_params: dict = None,
    reactions_net: dict = None,
    species_db: dict = None,
    mode: str = 'qRRHO',
    **kwargs
) -> dict:
    """
    Computes k_f and k_r using Marcus quadratic free-energy relation:
        ΔG‡_f = (λ / 4) · (1 + ΔG_rxn / λ)²
        k_f = (kB · T / h) · exp(-ΔG‡_f / RT)
        k_r = k_f / K_eq
    
    If lambda_eV is not specified, λ is set to 4 · E0 (from bep_params, default 3.20 eV),
    ensuring that Marcus and BEP share the exact same intrinsic barrier (E0 = 0.80 eV)
    and slope (α = 0.50) at ΔG_rxn = 0.
    """
    thermo = calculate_reaction_thermo(
        rxn_id, T_K,
        reactions_net=reactions_net,
        species_db=species_db,
        mode=mode,
        **kwargs
    )
    rxn_class = thermo['class']
    dG_rxn_eV = thermo['dG_rxn_eV']
    
    if lambda_eV is not None:
        lambda_val = max(lambda_eV, 0.01)
    else:
        all_bep = _get_bep_parameters(bep_params)
        default_fallback = DEFAULT_FAMILY_BEP_PARAMETERS.get(rxn_class, DEFAULT_FAMILY_BEP_PARAMETERS['default'])
        bep = all_bep.get(rxn_class, all_bep.get('default', default_fallback))
        lambda_val = 4.0 * bep.get('E0_eV', default_fallback['E0_eV']) # Default λ = 4 * E0 for Marcus-BEP consistency
        
    ratio = dG_rxn_eV / lambda_val
    dG_barrier_f_eV = (lambda_val / 4.0) * ((1.0 + ratio) ** 2)
    dG_barrier_f_kJ_mol = dG_barrier_f_eV * EV_TO_KJ_MOL
    
    eyring_prefactor = (KB_SI * T_K) / H_SI
    k_f = eyring_prefactor * np.exp(-(dG_barrier_f_kJ_mol * 1000.0) / (R_SI * T_K))
    
    K_eq = thermo['K_eq']
    k_r = k_f / K_eq if K_eq > 1e-300 else 0.0
    dG_barrier_r_eV = dG_barrier_f_eV - dG_rxn_eV
    
    return {
        'rxn_id': rxn_id,
        'T_K': T_K,
        'kinetic_model': 'marcus_eyring',
        'thermo_mode': mode,
        'class': rxn_class,
        'lambda_eV': lambda_val,
        'dG_rxn_eV': dG_rxn_eV,
        'dG_rxn_kJ_mol': thermo['dG_rxn_kJ_mol'],
        'dG_barrier_f_kJ_mol': dG_barrier_f_kJ_mol,
        'dG_barrier_f_eV': dG_barrier_f_eV,
        'dG_barrier_r_kJ_mol': dG_barrier_r_eV * EV_TO_KJ_MOL,
        'dG_barrier_r_eV': dG_barrier_r_eV,
        'k_f': k_f,
        'k_r': k_r,
        'K_eq': K_eq
    }

# ------------------------------------------------------------------------------
# DISPATCHER: calculate_rate_constants
# ------------------------------------------------------------------------------
KIN_MODELS = {
    'bep_eyring': bep_eyring,
    'bep': bep_eyring,
    'marcus_eyring': marcus_eyring,
    'marcus': marcus_eyring,
}

def calculate_rate_constants(
    rxn_id: str,
    T_K: float,
    model: str = None,
    kinetic_model: str = 'bep_eyring',
    bep_params: dict = None,
    reactions_net: dict = None,
    species_db: dict = None,
    mode: str = 'qRRHO',
    thermo_mode: str = None,
    **kwargs
) -> dict:
    """
    Main dispatcher for microkinetic rate constants supporting pluggable kinetic models:
    - 'bep_eyring' (default): Bell-Evans-Polanyi linear scaling + Eyring TST.
    - 'marcus_eyring': Marcus theory quadratic activation relation.
    
    Thermodynamic mode:
    - 'wb97mv': Direct ωB97M-V/def2-TZVPD dataset free energies.
    - 'qRRHO': Grimme quasi-RRHO statistical mechanics.
    - 'b3lyp_benchmark': Legacy alias for 'wb97mv'.
    """
    active_kin_model = model or kinetic_model
    active_thermo_mode = thermo_mode or mode
    
    if active_kin_model not in KIN_MODELS:
        raise ValueError(
            f"Unknown kinetic model '{active_kin_model}'. "
            f"Supported models: {list(KIN_MODELS.keys())}"
        )
        
    model_fn = KIN_MODELS[active_kin_model]
    return model_fn(
        rxn_id=rxn_id,
        T_K=T_K,
        bep_params=bep_params,
        reactions_net=reactions_net,
        species_db=species_db,
        mode=active_thermo_mode,
        **kwargs
    )