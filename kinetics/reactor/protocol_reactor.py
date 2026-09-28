"""Protocol Reactor Simulator: Multi-Stage Temperature-Programmed Stiff ODE Engine.

Simulates bench-scale operando kinetics across multi-stage temperature schedules
(e.g., pre-equilibration, discrete reagent injection/dilution, stepwise temperature
ramps from 20 °C to 80 °C, and chilled NMR acquisition holds).

Coupled directly with the rigorous statistical mechanics and microkinetics engine
(Quasi-RRHO, standard-state compression shift, and Wegscheider detailed balance).
"""

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from kinetics.microkinetics.rate_constants import calculate_rate_constants
from kinetics.microkinetics.kinetic_parameters import PETER_REFERENCE_PARAMETERS

# Default physical densities (g/mL) and molecular weights (g/mol)
DEFAULT_RHO = {
    "EC": 1.321,
    "H2O": 0.997,
    "TMSPA": 0.940,  # CAS 1017-19-2
}

DEFAULT_MW = {
    "EC": 88.062,
    "H2O": 18.015,
    "TMSPA": 314.54,
}

DEFAULT_TRACKED_SPECIES = [
    "TMSPA", "BMSPA", "MMSPA", "H3PO4", "H2O", "TMSOH",
    "HMDSO", "EC", "TMSOEG", "TMSOdiEG", "CO2"
]

P_PER_SPECIES = {"TMSPA": 1, "BMSPA": 1, "MMSPA": 1, "H3PO4": 1}
SI_PER_SPECIES = {
    "TMSPA": 3, "BMSPA": 2, "MMSPA": 1, "TMSOH": 1, "HMDSO": 2,
    "TMSOEG": 1, "TMSOdiEG": 1
}


def compute_recipe_molarities(
    h2o_vol_frac_stock: float = 0.02,
    tmspa_vol_frac: float = 0.05,
    rho_dict: dict = None,
    mw_dict: dict = None
) -> dict:
    """Computes exact starting solution molarities (mol/L) from volumetric recipe.
    
    Models two-step preparation:
    1. Stock solution: EC + wet additive (e.g. 2 vol% H2O).
    2. Injection of TMSPA (e.g. 5 vol% TMSPA), which physically dilutes the stock components.
    
    Parameters
    ----------
    h2o_vol_frac_stock : float
        Volume fraction of water in the initial EC/water stock (default: 0.02 = 2 vol%).
    tmspa_vol_frac : float
        Volume fraction of TMSPA added to the stock (default: 0.05 = 5 vol%).
    rho_dict : dict, optional
        Densities in g/mL.
    mw_dict : dict, optional
        Molar masses in g/mol.
        
    Returns
    -------
    dict
        Dictionary containing:
        - 'stock': dict of molarities before TMSPA addition (mol/L)
        - 'after': dict of molarities immediately after TMSPA addition (mol/L)
        - 'dilution_factor': (1 - tmspa_vol_frac)
        - 'summary_df': pandas DataFrame summarizing vol% and concentrations.
    """
    rho = rho_dict or DEFAULT_RHO
    mw = mw_dict or DEFAULT_MW
    
    def molarity(vol_frac, sp):
        return vol_frac * 1000.0 * rho[sp] / mw[sp]
        
    stock = {
        "H2O": molarity(h2o_vol_frac_stock, "H2O"),
        "EC": molarity(1.0 - h2o_vol_frac_stock, "EC"),
        "TMSPA": 0.0
    }
    
    dilution_factor = 1.0 - tmspa_vol_frac
    after = {
        "H2O": stock["H2O"] * dilution_factor,
        "EC": stock["EC"] * dilution_factor,
        "TMSPA": molarity(tmspa_vol_frac, "TMSPA")
    }
    
    summary_df = pd.DataFrame({
        "vol% (final)": {
            "EC": (1.0 - h2o_vol_frac_stock) * dilution_factor * 100.0,
            "H2O": h2o_vol_frac_stock * dilution_factor * 100.0,
            "TMSPA": tmspa_vol_frac * 100.0
        },
        "stock (M)": stock,
        "after injection (M)": after
    })
    
    return {
        "stock": stock,
        "after": after,
        "dilution_factor": dilution_factor,
        "h2o_tmspa_ratio": after["H2O"] / after["TMSPA"],
        "summary_df": summary_df
    }


def build_default_protocol_schedule(
    rt_c: float = 20.0,
    stir_h: float = 24.0,
    equilibrate_h: float = 12.0,
    hold_h: float = 8.0,
    acquire_h: float = 1.0,
    step_temps_c: list[float] = None
) -> tuple[list[tuple], pd.DataFrame]:
    """Builds the default multi-stage experimental bench protocol schedule."""
    temps = step_temps_c if step_temps_c is not None else [20, 30, 40, 50, 60, 70, 80]
    
    stages = [
        ("EC + H2O, stirred", rt_c, stir_h, None),
        ("equilibrate", rt_c, equilibrate_h, None)
    ]
    
    for t_step in temps:
        stages.append((f"hold {t_step} C", float(t_step), hold_h, None))
        stages.append((f"acquire ({t_step} C)", rt_c, acquire_h, float(t_step)))
        
    t_start = -(stir_h + equilibrate_h)
    sched_rows = []
    current_t = t_start
    for label, t_c, dur, spec in stages:
        sched_rows.append({
            "stage": label,
            "T (°C)": t_c,
            "duration (h)": dur,
            "start (h)": round(current_t, 2),
            "end (h)": round(current_t + dur, 2),
            "acquisition": f"{spec:.0f} °C" if spec is not None else ""
        })
        current_t += dur
        
    return stages, pd.DataFrame(sched_rows)


def simulate_protocol_reactor(
    stages: list = None,
    recipe: dict = None,
    injection_stage_idx: int = 2,
    reactions_net: dict = None,
    bep_params: dict = None,
    species_list: list = None,
    ec_buffered: bool = True,
    thermo_mode: str = "qRRHO",
    points_per_stage: int = 120,
    rtol: float = 1e-8,
    atol: float = 1e-12,
    **kwargs
) -> dict:
    """Simulates the complete multi-stage experimental protocol using stiff ODE integration."""
    if stages is None:
        stages, _ = build_default_protocol_schedule()
    if recipe is None:
        recipe = compute_recipe_molarities()
        
    species = species_list or DEFAULT_TRACKED_SPECIES
    ns = len(species)
    idx = {s: i for i, s in enumerate(species)}
    
    ec_reservoir = recipe["after"]["EC"]
    
    # Initial state: stock mixture (EC + H2O)
    c_current = np.zeros(ns, dtype=float)
    if "H2O" in idx:
        c_current[idx["H2O"]] = recipe["stock"]["H2O"]
    if "EC" in idx:
        c_current[idx["EC"]] = recipe["stock"]["EC"]
        
    # Pre-equilibration time offset so t=0 corresponds to TMSPA addition
    pre_eq_h = sum(st[2] for i, st in enumerate(stages) if i < injection_stage_idx)
    current_time_s = -pre_eq_h * 3600.0
    
    t_parts = []
    c_parts = []
    snapshots = []
    events = []
    
    effective_bep = bep_params if bep_params is not None else PETER_REFERENCE_PARAMETERS
    
    for stage_i, (label, t_celsius, duration_h, spec_temp) in enumerate(stages):
        # Event: discrete injection and dilution
        if stage_i == injection_stage_idx:
            c_current = c_current * recipe["dilution_factor"]
            if "TMSPA" in idx:
                c_current[idx["TMSPA"]] = recipe["after"]["TMSPA"]
            events.append(("TMSPA added", current_time_s / 3600.0))
            
        t_kelvin = t_celsius + 273.15
        duration_s = duration_h * 3600.0
        
        kf_list, kr_list, active_rxns = [], [], []
        
        rxn_items = list(reactions_net.items()) if reactions_net is not None else [
            ("R1", {"reactants": {"TMSPA": 1, "H2O": 1},     "products": {"BMSPA": 1, "TMSOH": 1}}),
            ("R2", {"reactants": {"BMSPA": 1, "H2O": 1},     "products": {"MMSPA": 1, "TMSOH": 1}}),
            ("R3", {"reactants": {"MMSPA": 1, "H2O": 1},     "products": {"H3PO4": 1, "TMSOH": 1}}),
            ("R4", {"reactants": {"TMSOH": 2},               "products": {"HMDSO": 1, "H2O": 1}}),
            ("R5", {"reactants": {"TMSPA": 1, "TMSOH": 1},   "products": {"BMSPA": 1, "HMDSO": 1}}),
            ("R6", {"reactants": {"BMSPA": 1, "TMSOH": 1},   "products": {"MMSPA": 1, "HMDSO": 1}}),
            ("R7", {"reactants": {"MMSPA": 1, "TMSOH": 1},   "products": {"H3PO4": 1, "HMDSO": 1}}),
            ("R8", {"reactants": {"EC": 1, "TMSOH": 1},      "products": {"TMSOEG": 1, "CO2": 1}}),
            ("R9", {"reactants": {"EC": 1, "TMSOEG": 1},     "products": {"TMSOdiEG": 1, "CO2": 1}}),
        ]
        
        for rxn_id, rxn_data in rxn_items:
            rates = calculate_rate_constants(
                rxn_id,
                T_K=t_kelvin,
                bep_params=effective_bep,
                reactions_net=reactions_net,
                mode=thermo_mode,
                **kwargs
            )
            kf_list.append(rates["k_f"])
            kr_list.append(rates["k_r"])
            active_rxns.append((rxn_id, rxn_data["reactants"], rxn_data["products"]))
            
        def ode_rhs(t, c):
            c_clamped = np.maximum(c, 0.0)
            dcdt = np.zeros(ns, dtype=float)
            
            def get_conc(s):
                if s not in idx:
                    return 1.0
                if s == "EC" and ec_buffered:
                    return ec_reservoir
                return c_clamped[idx[s]]
                
            for r_idx, (_, reac_dict, prod_dict) in enumerate(active_rxns):
                rate_fwd = kf_list[r_idx]
                for r_sp, nu_r in reac_dict.items():
                    rate_fwd *= (get_conc(r_sp) ** nu_r)
                    
                rate_rev = kr_list[r_idx]
                for p_sp, nu_p in prod_dict.items():
                    rate_rev *= (get_conc(p_sp) ** nu_p)
                    
                net_rate = rate_fwd - rate_rev
                
                for r_sp, nu_r in reac_dict.items():
                    if r_sp in idx and not (r_sp == "EC" and ec_buffered):
                        dcdt[idx[r_sp]] -= nu_r * net_rate
                for p_sp, nu_p in prod_dict.items():
                    if p_sp in idx and not (p_sp == "EC" and ec_buffered):
                        dcdt[idx[p_sp]] += nu_p * net_rate
                        
            return dcdt
            
        t_eval_stage = np.linspace(0.0, duration_s, points_per_stage)
        sol = solve_ivp(
            ode_rhs,
            [0.0, duration_s],
            c_current,
            method="Radau",
            t_eval=t_eval_stage,
            rtol=rtol,
            atol=atol
        )
        
        t_parts.append(current_time_s + sol.t)
        c_parts.append(sol.y)
        c_current = sol.y[:, -1]
        current_time_s += duration_s
        
        if spec_temp is not None:
            snapshots.append({
                "T_C": spec_temp,
                "t_h": current_time_s / 3600.0,
                "C": c_current.copy()
            })
            
    t_s_total = np.concatenate(t_parts)
    t_h_total = t_s_total / 3600.0
    c_m_total = np.concatenate(c_parts, axis=1)
    c_mm_total = c_m_total * 1000.0
    
    post_inj_start = injection_stage_idx * points_per_stage
    si_total = np.zeros(len(t_h_total))
    p_total = np.zeros(len(t_h_total))
    for s_name, mult in SI_PER_SPECIES.items():
        if s_name in idx:
            si_total += mult * c_m_total[idx[s_name]]
    for s_name, mult in P_PER_SPECIES.items():
        if s_name in idx:
            p_total += mult * c_m_total[idx[s_name]]
            
    si_post = si_total[post_inj_start:]
    p_post = p_total[post_inj_start:]
    si_conserved = bool(np.max(np.abs(si_post - si_post[0])) < 1e-5) if len(si_post) > 0 else True
    p_conserved = bool(np.max(np.abs(p_post - p_post[0])) < 1e-5) if len(p_post) > 0 else True
    
    return {
        "t_s": t_s_total,
        "t_h": t_h_total,
        "C_M": c_m_total,
        "C_mM": c_mm_total,
        "species": species,
        "idx": idx,
        "snapshots": snapshots,
        "events": events,
        "si_total_M": si_total,
        "p_total_M": p_total,
        "si_conserved": si_conserved,
        "p_conserved": p_conserved,
        "recipe": recipe,
        "stages": stages
    }
