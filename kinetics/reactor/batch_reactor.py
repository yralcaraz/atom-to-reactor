# ==============================================================================
# BLOCK 7: HOMOGENEOUS BATCH REACTOR SIMULATOR ("TANK" STIFF ODE ENGINE)
# ==============================================================================
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
try:
    from kinetics.microkinetics.rate_constants import calculate_rate_constants
except ImportError:
    from calculate_rate_constants import calculate_rate_constants

def build_stoichiometric_matrix(reactions_net: dict, tracked_species: list) -> np.ndarray:
    """
    Constructs the stoichiometric matrix S (N_species x N_reactions).
    S_ij > 0 for products, S_ij < 0 for reactants.
    """
    ns = len(tracked_species)
    nr = len(reactions_net)
    S = np.zeros((ns, nr))
    idx = {s: i for i, s in enumerate(tracked_species)}
    
    for j, (rxn_id, data) in enumerate(reactions_net.items()):
        for r, coeff in data['reactants'].items():
            if r in idx:
                S[idx[r], j] -= coeff
        for p, coeff in data['products'].items():
            if p in idx:
                S[idx[p], j] += coeff
    return S

def simulate_tank_reactor(
    C0_dict: dict,
    t_end_s: float = 1e6,
    T_K: float = 298.15,
    reactions_net: dict = None,
    bep_params: dict = None,
    species_db: dict = None,
    ec_buffered: bool = True,
    method: str = 'Radau',
    rtol: float = 1e-8,
    atol: float = 1e-12,
    n_points: int = 500,
    mode: str = None,
    t_start_s: float = 1e-3,
    **kwargs
) -> dict:
    """
    Integrates the stiff ODE system dC/dt = S · r(C, T) for a homogeneous batch reactor.
    
    Parameters:
    -----------
    C0_dict : dict of initial concentrations in mol/L (M)
    t_end_s : total simulation time in seconds (default 1e6 s ~ 11.5 days)
    T_K : isothermal reactor temperature in Kelvin (default 298.15 K)
    reactions_net : reaction network dictionary (inherits from calling scope if None)
    bep_params : BEP barrier and alpha parameters
    species_db : species database dictionary
    ec_buffered : if True, EC solvent concentration is held constant (reservoir)
    method : stiff ODE solver method ('Radau', 'BDF', etc.)
    rtol, atol : relative and absolute solver tolerances
    n_points : number of evaluation points
    mode : thermodynamic mode ('wb97mv', 'qRRHO', or legacy 'b3lyp_benchmark'). If None, inherits
           'selected_thermo_mode' or 'selected_mode' from calling scope (default: 'wb97mv').
    t_start_s : initial evaluation time in seconds for logarithmic time grid (default: 1e-3 s = 1 ms).
                Captures initial fast transients from t=0 while remaining strictly positive for log-scale plotting.
    """
    if reactions_net is None:
        import inspect
        frame = inspect.currentframe().f_back
        while frame:
            if 'reactions_network' in frame.f_globals:
                reactions_net = frame.f_globals['reactions_network']
                break
            if 'reactions_network' in frame.f_locals:
                reactions_net = frame.f_locals['reactions_network']
                break
            frame = frame.f_back

    if mode is None:
        import inspect
        frame = inspect.currentframe().f_back
        while frame:
            for var_name in ['selected_thermo_mode', 'selected_mode']:
                if var_name in frame.f_globals:
                    mode = frame.f_globals[var_name]
                    break
                if var_name in frame.f_locals:
                    mode = frame.f_locals[var_name]
                    break
            if mode is not None:
                break
            frame = frame.f_back
        if mode is None:
            mode = 'wb97mv'
            
    tracked_species = list(C0_dict.keys())
    idx = {s: i for i, s in enumerate(tracked_species)}
    ns = len(tracked_species)
    
    # Pre-calculate rate constants at temperature T_K
    KF = []
    KR = []
    rxn_ids = list(reactions_net.keys())
    for rxn_id in rxn_ids:
        rates = calculate_rate_constants(
            rxn_id, T_K,
            bep_params=bep_params,
            reactions_net=reactions_net,
            species_db=species_db,
            mode=mode,
            **kwargs
        )
        KF.append(rates['k_f'])
        KR.append(rates['k_r'])
        
    KF = np.array(KF)
    KR = np.array(KR)
    
    # Initial concentration vector
    C0 = np.array([C0_dict[s] for s in tracked_species], dtype=float)
    C_EC_const = C0_dict.get('EC', 4.5)
    
    def odes(t, C):
        C_pos = np.maximum(C, 0.0)
        dCdt = np.zeros(ns)
        
        def get_conc(s):
            if s == 'EC' and ec_buffered:
                return C_EC_const
            if s in idx:
                return C_pos[idx[s]]
            return 1.0
        
        for j, rxn_id in enumerate(rxn_ids):
            data = reactions_net[rxn_id]
            # Forward flux
            rf = KF[j]
            for r, coeff in data['reactants'].items():
                rf *= (get_conc(r) ** coeff)
            # Reverse flux
            rr = KR[j]
            for p, coeff in data['products'].items():
                rr *= (get_conc(p) ** coeff)
                
            net_flux = rf - rr
            
            for r, coeff in data['reactants'].items():
                if r in idx and not (r == 'EC' and ec_buffered):
                    dCdt[idx[r]] -= coeff * net_flux
            for p, coeff in data['products'].items():
                if p in idx and not (p == 'EC' and ec_buffered):
                    dCdt[idx[p]] += coeff * net_flux
                    
        return dCdt

    if 't_eval' in kwargs:
        t_eval = kwargs.pop('t_eval')
    else:
        t_eval = np.logspace(np.log10(t_start_s), np.log10(t_end_s), n_points)
    
    t_span = kwargs.pop('t_span', [0.0, t_end_s])
    
    sol = solve_ivp(
        odes,
        t_span,
        C0,
        method=method,
        t_eval=t_eval,
        rtol=rtol,
        atol=atol
    )
    
    # Invariant element conservation verification
    si_counts = {'TMSPA': 3, 'BMSPA': 2, 'MMSPA': 1, 'TMSOH': 1, 'siloxyl': 2, 'TMSOEG': 1, 'TMSOdiEG': 1}
    p_counts  = {'TMSPA': 1, 'BMSPA': 1, 'MMSPA': 1, 'H3PO4': 1}
    
    si_tot = np.zeros_like(sol.t)
    p_tot = np.zeros_like(sol.t)
    for s, n_si in si_counts.items():
        if s in idx:
            si_tot += n_si * sol.y[idx[s]]
    for s, n_p in p_counts.items():
        if s in idx:
            p_tot += n_p * sol.y[idx[s]]
            
    si_conserved = bool(np.max(np.abs(si_tot - si_tot[0])) < 1e-6)
    p_conserved = bool(np.max(np.abs(p_tot - p_tot[0])) < 1e-6)

    return {
        't_s': sol.t,
        't_h': sol.t / 3600.0,
        't_d': sol.t / 86400.0,
        'C_M': sol.y,
        'C_mM': sol.y * 1000.0,
        'species': tracked_species,
        'idx': idx,
        'si_total_M': si_tot,
        'p_total_M': p_tot,
        'si_conserved': si_conserved,
        'p_conserved': p_conserved,
        'message': sol.message,
        'success': sol.success,
        'mode': mode
    }
