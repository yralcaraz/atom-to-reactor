"""Reaction network: the 9 elementary steps of TMSPA scavenging in wet EC, and stoichiometry helpers.

A network is a dict {rxn_id: {'reactants': {species: nu}, 'products': {species: nu}, 'class': family}}.
An optional 'canonical': False flag marks a reaction written against its family's reference direction
(Level 1 engine only). Species names are snapshot `pipeline_id`s.

The fits take one of three variants of the network (the active one, set per fit run by
kinetics.fitting.config, like the solvation medium):

    'R1-R9'                          the nine steps (NETWORK); fits 1–4
    'R1-R9+N10-N15'                  with the six further steps of P. Broqvist's generated network (fifth fit)
    'R1-R9+N10-N15, Peter barriers'  the same, the six all in the fixed-barrier class 'exchange'

    build_network                the network of a variant
    get_active_network           the network of the active variant
    set_active_network_variant   change it (None: back to 'R1-R9')

Source: Y. Alcaraz Galván; N10–N15 from P. Broqvist's tmspa_reaction_network_model notebook (2026-10-08)
"""

import numpy as np
from scipy.linalg import null_space

NETWORK = {
    'R1': {'reactants': {'TMSPA': 1, 'H2O': 1},   'products': {'BMSPA': 1, 'TMSOH': 1},    'class': 'hydrolysis'},
    'R2': {'reactants': {'BMSPA': 1, 'H2O': 1},   'products': {'MMSPA': 1, 'TMSOH': 1},    'class': 'hydrolysis'},
    'R3': {'reactants': {'MMSPA': 1, 'H2O': 1},   'products': {'H3PO4': 1, 'TMSOH': 1},    'class': 'hydrolysis'},
    'R4': {'reactants': {'TMSOH': 2},             'products': {'HMDSO': 1, 'H2O': 1},      'class': 'condensation'},
    'R5': {'reactants': {'TMSPA': 1, 'TMSOH': 1}, 'products': {'BMSPA': 1, 'HMDSO': 1},    'class': 'transfer'},
    'R6': {'reactants': {'BMSPA': 1, 'TMSOH': 1}, 'products': {'MMSPA': 1, 'HMDSO': 1},    'class': 'transfer'},
    'R7': {'reactants': {'MMSPA': 1, 'TMSOH': 1}, 'products': {'H3PO4': 1, 'HMDSO': 1},    'class': 'transfer'},
    'R8': {'reactants': {'EC': 1, 'TMSOH': 1},    'products': {'TMSOEG': 1, 'CO2': 1},     'class': 'solvent_attack'},
    'R9': {'reactants': {'EC': 1, 'TMSOEG': 1},   'products': {'TMSOdiEG': 1, 'CO2': 1},   'class': 'solvent_attack'},
}

# Tracked species in plotting order: phosphate ladder, water, silanol/siloxane, solvent and its products
NETWORK_SPECIES = [
    'TMSPA', 'BMSPA', 'MMSPA', 'H3PO4', 'H2O', 'TMSOH', 'HMDSO', 'EC', 'TMSOEG', 'TMSOdiEG', 'CO2'
]

# The six steps between the same species that P. Broqvist's generated network adds. Each is a combination of
# R1–R9 (N10 = R2 − R1, N11 = R9 − R8, N12 = R2 + R7, N13 = R2 − R3, N14 = R3 − R1, N15 = R1 + R6), so its ΔG
# follows from the species energies: they add routes, not end points. 'follows': the reaction whose barrier the
# step takes (N12 and N15 are hydrolyses that release HMDSO instead of TMSOH). 'exchange': a silyl group moving
# between two phosphates or two glycols, at a fixed barrier (kinetics.fitting.candidates.FIXED_BARRIERS_EV)
EXTRA_REACTIONS = {
    'N10': {'reactants': {'BMSPA': 2},            'products': {'MMSPA': 1, 'TMSPA': 1},    'class': 'exchange'},
    'N11': {'reactants': {'TMSOEG': 2},           'products': {'TMSOH': 1, 'TMSOdiEG': 1}, 'class': 'exchange'},
    'N12': {'reactants': {'BMSPA': 1, 'H2O': 1},  'products': {'H3PO4': 1, 'HMDSO': 1},    'class': 'hydrolysis', 'follows': 'R2'},
    'N13': {'reactants': {'BMSPA': 1, 'H3PO4': 1}, 'products': {'MMSPA': 2},               'class': 'exchange'},
    'N14': {'reactants': {'BMSPA': 1, 'MMSPA': 1}, 'products': {'H3PO4': 1, 'TMSPA': 1},   'class': 'exchange'},
    'N15': {'reactants': {'TMSPA': 1, 'H2O': 1},  'products': {'MMSPA': 1, 'HMDSO': 1},    'class': 'hydrolysis', 'follows': 'R1'},
}
NETWORK_VARIANTS = ('R1-R9', 'R1-R9+N10-N15', 'R1-R9+N10-N15, Peter barriers')
DEFAULT_NETWORK_VARIANT = 'R1-R9'
_NETWORKS = {}
_ACTIVE_VARIANT = [DEFAULT_NETWORK_VARIANT]


def build_network(variant: str = DEFAULT_NETWORK_VARIANT) -> dict:
    """The network of one variant (NETWORK_VARIANTS). Shared between calls: copy it before changing it."""
    if variant not in NETWORK_VARIANTS:
        raise ValueError(f"network variant must be one of {NETWORK_VARIANTS}, not '{variant}'")
    if variant not in _NETWORKS:
        net = dict(NETWORK)
        if variant != DEFAULT_NETWORK_VARIANT:
            for rxn_id, rxn in EXTRA_REACTIONS.items():
                net[rxn_id] = dict(rxn)
                if variant.endswith('Peter barriers'):          # as in his notebook: one barrier for all six
                    net[rxn_id]['class'] = 'exchange'
                    net[rxn_id].pop('follows', None)
        _NETWORKS[variant] = net
    return _NETWORKS[variant]


def get_active_network_variant() -> str:
    return _ACTIVE_VARIANT[0]


def set_active_network_variant(variant: str = None) -> str:
    """Make `variant` the network the fit structures use (None: 'R1-R9') and return it."""
    variant = variant or DEFAULT_NETWORK_VARIANT
    build_network(variant)                                      # refuses an unknown variant
    _ACTIVE_VARIANT[0] = variant
    return variant


def get_active_network() -> dict:
    return build_network(_ACTIVE_VARIANT[0])


def list_species(network: dict) -> list:
    """Species appearing in the network, in first-appearance order."""
    seen = {}
    for rxn in network.values():
        for side in ('reactants', 'products'):
            for sp in rxn[side]:
                seen.setdefault(sp, None)
    return list(seen)


def build_stoichiometric_matrix(network: dict, species: list) -> np.ndarray:
    """S (n_species × n_reactions): products positive, reactants negative. Untracked species are skipped."""
    idx = {sp: i for i, sp in enumerate(species)}
    S = np.zeros((len(species), len(network)))
    for j, rxn in enumerate(network.values()):
        for sp, nu in rxn['reactants'].items():
            if sp in idx:
                S[idx[sp], j] -= nu
        for sp, nu in rxn['products'].items():
            if sp in idx:
                S[idx[sp], j] += nu
    return S


def build_order_matrices(network: dict, species: list) -> tuple[np.ndarray, np.ndarray]:
    """Mass-action orders (n_species × n_reactions) of the forward (reactant) and reverse (product) rates."""
    idx = {sp: i for i, sp in enumerate(species)}
    orders_f = np.zeros((len(species), len(network)))
    orders_r = np.zeros((len(species), len(network)))
    for j, rxn in enumerate(network.values()):
        for sp, nu in rxn['reactants'].items():
            orders_f[idx[sp], j] += nu
        for sp, nu in rxn['products'].items():
            orders_r[idx[sp], j] += nu
    return orders_f, orders_r


def find_reaction_cycles(network: dict, tol: float = 1e-10) -> np.ndarray:
    """Basis of reaction combinations with zero net change (null space of S), shape (n_cycles, n_reactions).

    Wegscheider's condition requires Σ_j c_j ΔG_j = 0 along each cycle c. The basis is returned in
    reduced row-echelon form, so simple cycles such as R1 + R4 − R5 read directly.
    """
    S = build_stoichiometric_matrix(network, list_species(network))
    basis = null_space(S).T
    row = 0
    for col in range(basis.shape[1]):
        if row == basis.shape[0]:
            break
        pivot = row + int(np.argmax(np.abs(basis[row:, col])))
        if abs(basis[pivot, col]) < tol:
            continue
        basis[[row, pivot]] = basis[[pivot, row]]
        basis[row] /= basis[row, col]
        for i in range(basis.shape[0]):
            if i != row:
                basis[i] -= basis[i, col] * basis[row]
        row += 1
    return np.round(basis, 10) + 0.0   # coefficients are rational; + 0.0 turns -0.0 into 0.0


def expand_stoichiometry(side: dict) -> list:
    """{'TMSOH': 2} → ['TMSOH', 'TMSOH']."""
    return [sp for sp, nu in side.items() for _ in range(int(nu))]


def format_equation(rxn: dict, arrow: str = ' ⇌ ') -> str:
    """'2 TMSOH ⇌ HMDSO + H2O'."""
    def side(d):
        return ' + '.join(f"{nu} {sp}" if nu > 1 else sp for sp, nu in d.items())
    return side(rxn['reactants']) + arrow + side(rxn['products'])
