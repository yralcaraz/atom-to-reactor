# ==============================================================================
# BLOCK 5: KINETIC RATE CONSTANTS ENGINE (WITH STRICT DETAILED BALANCE)
# ==============================================================================
import inspect
import numpy as np
try:
    from kinetics.thermo.reaction_thermo import calculate_reaction_thermo, _get_reactions_network
    from kinetics.thermo.gas_thermo import R_SI, KB_SI, H_SI, EV_TO_KJ_MOL
    from kinetics.microkinetics.barrier_models import barrier, reversed_params
    from kinetics.microkinetics.kinetic_parameters import (
        DEFAULT_FAMILY_BEP_PARAMETERS,
        LEVEL1_FAMILY_PARAMETERS,
        PETER_REFERENCE_PARAMETERS,
        normalize_family_params,
    )
except ImportError:
    from calculate_reaction_thermo import calculate_reaction_thermo, _get_reactions_network
    from calculate_gas_thermo import R_SI, KB_SI, H_SI, EV_TO_KJ_MOL
    from barrier_models import barrier, reversed_params
    from kinetic_parameters import (
        DEFAULT_FAMILY_BEP_PARAMETERS,
        LEVEL1_FAMILY_PARAMETERS,
        PETER_REFERENCE_PARAMETERS,
        normalize_family_params,
    )

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

def _resolve_family_parameters(rxn_class, all_params, default_set):
    """Family entry for rxn_class, falling back to the set's 'default', then to default_set."""
    fallback = default_set.get(rxn_class, default_set['default'])
    return normalize_family_params(all_params.get(rxn_class, all_params.get('default', fallback)))

# ------------------------------------------------------------------------------
# CORE: barrier model + Eyring TST + optional diffusion ceiling + detailed balance
# ------------------------------------------------------------------------------
def _barrier_eyring(
    rxn_id: str,
    T_K: float,
    kinetic_model: str,
    shape: str = None,
    bep_params: dict = None,
    default_set: dict = None,
    reactions_net: dict = None,
    species_db: dict = None,
    mode: str = 'qRRHO',
    viscosity_Pa_s: float = None,
    g_override_eV: float = None,
    **kwargs
) -> dict:
    """
    ΔG‡_f = wR + F(ΔG_rxn - wR + wP; g(T)),  g(T) = g(T_ref) - (T - T_ref)·ΔS‡
    k_f   = (kB·T / h) · exp(-ΔG‡_f / RT)       [optionally combined in series with k_D]
    k_r   = k_f / K_eq

    shape=None reads the barrier model from each family's 'shape' key (default 'marcus').
    Reactions flagged 'canonical': False in the network use the reversed family parameters.
    """
    thermo = calculate_reaction_thermo(
        rxn_id, T_K,
        reactions_net=reactions_net,
        species_db=species_db,
        mode=mode,
        **kwargs
    )
    rxn_class = thermo['class']
    fam = _resolve_family_parameters(rxn_class, bep_params, default_set)
    active_shape = shape or fam.get('shape', 'marcus')

    g_eV = g_override_eV if g_override_eV is not None else fam['g_eV']
    g_eV -= (T_K - fam['T_ref_K']) * fam['dS_act_J_molK'] / 1000.0 / EV_TO_KJ_MOL

    rxn_entry = _get_reactions_network(reactions_net)[rxn_id]
    if not rxn_entry.get('canonical', True):
        fam = reversed_params(fam)

    model_params = {'wR_eV': fam['wR_eV'], 'wP_eV': fam['wP_eV']}
    if active_shape == 'bep_cap':
        model_params['alpha'] = fam['alpha']
    elif active_shape == 'blowers_masel':
        model_params['w'] = fam['w_eV']
    elif active_shape == 'two_parabola':
        model_params['alpha0'] = fam['alpha0']

    dG_rxn_kJ_mol = thermo['dG_rxn_kJ_mol']
    dG_rxn_eV = thermo['dG_rxn_eV']
    dG_barrier_f_eV, alpha_eff = barrier(dG_rxn_eV, active_shape, g_eV, **model_params)
    dG_barrier_f_kJ_mol = dG_barrier_f_eV * EV_TO_KJ_MOL

    # Eyring forward rate: k_f = (kB*T / h) * exp(-ΔG‡ / RT)
    eyring_prefactor = (KB_SI * T_K) / H_SI
    k_f = eyring_prefactor * np.exp(-(dG_barrier_f_kJ_mol * 1000.0) / (R_SI * T_K))

    # Collins-Kimball: bimolecular steps cannot exceed the Smoluchowski limit k_D = 8RT / 3η
    k_D = None
    if viscosity_Pa_s is not None and sum(rxn_entry['reactants'].values()) == 2:
        k_D = 8.0 * R_SI * T_K / (3.0 * viscosity_Pa_s) * 1000.0   # m³/(mol·s) → M⁻¹·s⁻¹
        k_f = 1.0 / (1.0 / k_f + 1.0 / k_D)
        dG_barrier_f_kJ_mol = -R_SI * T_K * np.log(k_f / eyring_prefactor) / 1000.0
        dG_barrier_f_eV = dG_barrier_f_kJ_mol / EV_TO_KJ_MOL

    # Detailed balance, applied after any ceiling: k_r = k_f / K_eq
    K_eq = thermo['K_eq']
    k_r = k_f / K_eq if K_eq > 1e-300 else 0.0
    dG_barrier_r_kJ_mol = dG_barrier_f_kJ_mol - dG_rxn_kJ_mol

    result = {
        'rxn_id': rxn_id,
        'T_K': T_K,
        'kinetic_model': kinetic_model,
        'barrier_model': active_shape,
        'thermo_mode': mode,
        'class': rxn_class,
        'g_eV': g_eV,
        'alpha_eff': alpha_eff,
        'dG_rxn_eV': dG_rxn_eV,
        'dG_rxn_kJ_mol': dG_rxn_kJ_mol,
        'dG_barrier_f_kJ_mol': dG_barrier_f_kJ_mol,
        'dG_barrier_f_eV': dG_barrier_f_eV,
        'dG_barrier_r_kJ_mol': dG_barrier_r_kJ_mol,
        'dG_barrier_r_eV': dG_barrier_r_kJ_mol / EV_TO_KJ_MOL,
        'k_f': k_f,
        'k_r': k_r,
        'K_eq': K_eq
    }
    if active_shape == 'marcus':
        result['lambda_eV'] = 4.0 * g_eV
    if k_D is not None:
        result['k_D'] = k_D
    return result

# ------------------------------------------------------------------------------
# MODEL 1: Bell-Evans-Polanyi (BEP) Linear Scaling + Eyring TST (legacy)
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
    Legacy form: kinetics depend on the direction a reaction is written (see Level 1 formulation).
    """
    return _barrier_eyring(
        rxn_id, T_K, kinetic_model='bep_eyring', shape='bep_cap',
        bep_params=_get_bep_parameters(bep_params),
        default_set=DEFAULT_FAMILY_BEP_PARAMETERS,
        reactions_net=reactions_net, species_db=species_db, mode=mode,
        **kwargs
    )

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

    If lambda_eV is not specified, λ is set to 4 · E0 of the reaction family, so that
    Marcus and BEP share the same intrinsic barrier and slope (α = 0.50) at ΔG_rxn = 0.
    """
    g_override = max(lambda_eV, 0.01) / 4.0 if lambda_eV is not None else None
    return _barrier_eyring(
        rxn_id, T_K, kinetic_model='marcus_eyring', shape='marcus',
        bep_params=_get_bep_parameters(bep_params),
        default_set=DEFAULT_FAMILY_BEP_PARAMETERS,
        reactions_net=reactions_net, species_db=species_db, mode=mode,
        g_override_eV=g_override,
        **kwargs
    )

# ------------------------------------------------------------------------------
# LEVEL 1 MODELS: smooth, reversal-invariant barrier relations
# ------------------------------------------------------------------------------
def _smooth_model(shape: str, kinetic_model: str):
    def model_fn(rxn_id, T_K, bep_params=None, reactions_net=None, species_db=None, mode='qRRHO', **kwargs):
        return _barrier_eyring(
            rxn_id, T_K, kinetic_model=kinetic_model, shape=shape,
            bep_params=_get_bep_parameters(bep_params),
            default_set=DEFAULT_FAMILY_BEP_PARAMETERS,
            reactions_net=reactions_net, species_db=species_db, mode=mode,
            **kwargs
        )
    model_fn.__name__ = f'{shape}_eyring'
    model_fn.__doc__ = f"'{shape}' barrier model + Eyring TST with detailed balance (family g from bep_params)."
    return model_fn

agmon_levine_eyring = _smooth_model('agmon_levine', 'agmon_levine')
blowers_masel_eyring = _smooth_model('blowers_masel', 'blowers_masel')
two_parabola_eyring = _smooth_model('two_parabola', 'two_parabola')

def level1_eyring(
    rxn_id: str,
    T_K: float,
    bep_params: dict = None,
    reactions_net: dict = None,
    species_db: dict = None,
    mode: str = 'qRRHO',
    **kwargs
) -> dict:
    """
    Level 1 engine: per-family barrier model ('shape', default Marcus), intrinsic barrier g(T),
    work terms and canonical-direction handling. Without bep_params it uses
    LEVEL1_FAMILY_PARAMETERS directly (notebook globals are not inspected).
    """
    return _barrier_eyring(
        rxn_id, T_K, kinetic_model='level1', shape=None,
        bep_params=bep_params if bep_params is not None else LEVEL1_FAMILY_PARAMETERS,
        default_set=LEVEL1_FAMILY_PARAMETERS,
        reactions_net=reactions_net, species_db=species_db, mode=mode,
        **kwargs
    )

# ------------------------------------------------------------------------------
# DISPATCHER: calculate_rate_constants
# ------------------------------------------------------------------------------
KIN_MODELS = {
    'bep_eyring': bep_eyring,
    'bep': bep_eyring,
    'marcus_eyring': marcus_eyring,
    'marcus': marcus_eyring,
    'agmon_levine': agmon_levine_eyring,
    'blowers_masel': blowers_masel_eyring,
    'two_parabola': two_parabola_eyring,
    'level1': level1_eyring,
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
    - 'bep_eyring' (default, legacy): Bell-Evans-Polanyi linear scaling with floor at E0.
    - 'marcus_eyring': Marcus theory quadratic activation relation.
    - 'agmon_levine', 'blowers_masel', 'two_parabola': smooth Level 1 barrier relations.
    - 'level1': per-family shapes and LEVEL1_FAMILY_PARAMETERS (Marcus, recalibrated g).

    Optional keyword: viscosity_Pa_s enables the Collins-Kimball diffusion ceiling.

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
