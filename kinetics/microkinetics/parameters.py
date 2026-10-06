"""Block 5 — Family kinetic parameter sets.

Parameter schema per reaction family (all keys optional except the barrier):
    'g_eV'                      intrinsic barrier at T_ref (ΔG‡ of a thermoneutral member);
                                'E0_eV' is accepted as a legacy synonym
    'alpha'                     Brønsted slope, read only by the legacy 'bep_cap' shape
    'shape'                     barrier model; used by the 'level1' kinetic model
    'T_ref_K', 'dS_act_J_mol_K' g(T) = g(T_ref) - (T - T_ref)·ΔS‡
    'alpha0'                    α(0) for 'two_parabola'
    'w_eV'                      bond-energy parameter for 'blowers_masel'
    'wR_eV', 'wP_eV'            Marcus work terms (precursor / successor complexes)
    'prior', 'source'           provenance and calibration intervals (not used in rates)
The 'default' entry applies to any family without its own entry.

Source: Y. Alcaraz Galván; barrier windows from Gogoi et al., J. Phys. Chem. C 2024, 128, 1654
"""

from copy import deepcopy

_GOGOI_2024 = "Gogoi et al., J. Phys. Chem. C 2024, 128, 1654"

# Legacy family parameters (capped BEP). Values unchanged from the pre-Level 1 engine.
FAMILY_BEP_PARAMETERS = {
    'hydrolysis':     {'g_eV': 0.80, 'alpha': 0.50},
    'condensation':   {'g_eV': 0.80, 'alpha': 0.50},
    'transfer':       {'g_eV': 0.80, 'alpha': 0.50},
    # Observed barrier derived from the 80 °C protocol of Gogoi 2024 (not a value reported in the paper)
    'solvent_attack': {'g_eV': 1.30, 'alpha': 0.50},
    'default':        {'g_eV': 0.80, 'alpha': 0.50},
}

# Peter Broqvist's protocol reference: one global E0 chosen so
# that TMSPA depletion spreads across the RT → 80 °C temperature steps.
PETER_REFERENCE_PARAMETERS = {
    'default': {'g_eV': 1.15, 'alpha': 0.50},
}

# Level 1: Marcus barriers with intrinsic barriers recalibrated for the smooth form.
# See TFM "KIN - DRAFT - Level 1 formulation refined scaling relations - 260928", §4.4-4.5.
LEVEL1_PARAMETERS = {
    'hydrolysis': {
        'shape': 'marcus', 'g_eV': 0.80, 'T_ref_K': 298.15, 'dS_act_J_mol_K': 0.0,
        'prior': {'type': 'interval', 'low_eV': 0.60, 'high_eV': 1.10},
        'source': "Placeholder from Broqvist tank_model.ipynb; no quantitative anchor yet",
    },
    'transfer': {
        'shape': 'marcus', 'g_eV': 0.80, 'T_ref_K': 298.15, 'dS_act_J_mol_K': 0.0,
        'prior': {'type': 'upper_bound', 'high_eV': 1.10},
        'source': f"Upper bound from TMSPa + TMSOH consumption at RT ({_GOGOI_2024}); derived",
    },
    'condensation': {
        'shape': 'marcus', 'g_eV': 1.30, 'T_ref_K': 298.15, 'dS_act_J_mol_K': 0.0,
        'prior': {'type': 'lower_bound', 'low_eV': 1.25},
        'source': f"No TMSOTMS from TMSOH alone after 8 h at 80 °C ({_GOGOI_2024}); "
                  "uncatalysed channel only; derived",
    },
    'solvent_attack': {
        'shape': 'marcus', 'g_eV': 1.32, 'T_ref_K': 298.15, 'dS_act_J_mol_K': 0.0,
        'prior': {'type': 'interval', 'low_eV': 1.27, 'high_eV': 1.38},
        'source': f"EC ring-opening by TMSOH absent at RT, clear after 8 h at 80 °C ({_GOGOI_2024}); derived",
    },
    'default': {
        'shape': 'marcus', 'g_eV': 0.80, 'T_ref_K': 298.15, 'dS_act_J_mol_K': 0.0,
        'source': "Fallback for unclassified reactions",
    },
}

_FAMILY_DEFAULTS = {
    'alpha': 0.50,
    'T_ref_K': 298.15,
    'dS_act_J_mol_K': 0.0,
    'alpha0': 0.50,
    'w_eV': 5.0,
    'wR_eV': 0.0,
    'wP_eV': 0.0,
}


_ALLOWED_KEYS = set(_FAMILY_DEFAULTS) | {'g_eV', 'E0_eV', 'shape', 'prior', 'source'}


def normalize_family_params(params: dict) -> dict:
    """Copy with 'g_eV' resolved (legacy 'E0_eV' accepted) and defaults filled in; unknown keys raise."""
    unknown = set(params) - _ALLOWED_KEYS
    if unknown:
        raise KeyError(f"Unknown family parameter(s) {sorted(unknown)}. Allowed: {sorted(_ALLOWED_KEYS)}")
    out = deepcopy(_FAMILY_DEFAULTS)
    out.update(params)
    if 'g_eV' not in out:
        if 'E0_eV' not in out:
            raise KeyError("Family parameters need an intrinsic barrier: 'g_eV' or 'E0_eV'")
        out['g_eV'] = out['E0_eV']
    return out


def resolve_family_params(rxn_class: str, family_params: dict, fallback_set: dict) -> dict:
    """Normalized entry for rxn_class: the family, else the set's 'default', else fallback_set."""
    fallback = fallback_set.get(rxn_class, fallback_set['default'])
    return normalize_family_params(family_params.get(rxn_class, family_params.get('default', fallback)))
