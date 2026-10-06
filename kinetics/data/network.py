"""Reaction network: the 9 elementary steps of TMSPA scavenging in wet EC, and stoichiometry helpers.

A network is a dict {rxn_id: {'reactants': {species: nu}, 'products': {species: nu}, 'class': family}}.
An optional 'canonical': False flag marks a reaction written against its family's reference direction
(Level 1 engine only). Species names are snapshot `pipeline_id`s.

Source: Y. Alcaraz Galván
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
