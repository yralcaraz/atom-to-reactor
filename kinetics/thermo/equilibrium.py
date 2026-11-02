"""Block 4C — Equilibrium composition of a closed tube: no reactor, no time.

The reactions that reach equilibrium are given as a set (the thermodynamic fit: R1–R4). Every reaction of the network
that is a combination of them (R5–R7; N10 and N12–N15 of the generated network) is in equilibrium with them, so
the composition is the one that minimises the Gibbs energy of the solution under the amounts that these reactions
conserve (phosphorus, silyl groups, O–H groups):

    min  Σ_i C_i (μ_i + RT ln C_i − RT)     subject to  Aᵀ C = Aᵀ c0
    μ_i = G_sol,i (kinetics.thermo.solution, 1 M standard state, no activity coefficients: the reactor's convention)

A spans the amounts conserved by the reactions (the null space of their stoichiometry). It is solved in its dual,
which is convex and smooth: C_i = exp(−μ_i/RT + (Aλ)_i), λ found by Newton on Σ C_i − λ·Aᵀc0. Every C stays positive,
and a species far below detection (TMSPA in a wet tube) is computed without dividing by zero. A species that no
reaction can make from the composition at mixing stays at zero; one outside the reactions (EC, CO₂, TMSOEG when R8 and
R9 are not in the set) keeps its amount at mixing.

    list_equilibrium_reactions  the reactions of a network that are combinations of the given ones
    calculate_equilibrium       composition at equilibrium of one tube

Source: Y. Alcaraz Galván
"""

import numpy as np
from scipy.linalg import null_space

from kinetics.constants import EV_TO_KJ_MOL, R_SI
from kinetics.data.network import build_stoichiometric_matrix, get_active_network, list_species
from kinetics.thermo.solution import calculate_solution_gibbs


def list_equilibrium_reactions(reactions, network: dict = None) -> list:
    """The reactions of `network` (default: the active one) whose stoichiometry is a combination of `reactions`:
    the given ones and those that follow them (Wegscheider)."""
    net = network or get_active_network()
    species = list_species({r: net[r] for r in reactions})
    base = build_stoichiometric_matrix({r: net[r] for r in reactions}, species)
    rank = np.linalg.matrix_rank(base)
    out = []
    for rxn_id, rxn in net.items():
        if not set(rxn['reactants']) | set(rxn['products']) <= set(species):
            continue
        column = build_stoichiometric_matrix({rxn_id: rxn}, species)
        if np.linalg.matrix_rank(np.hstack([base, column])) == rank:
            out.append(rxn_id)
    return out


def _reachable(c0: dict, network: dict, reactions: list) -> set:
    """Species present at mixing or made from them by the reactions, run either way."""
    present = {sp for sp, c in c0.items() if c > 0}
    changed = True
    while changed:
        changed = False
        for rxn_id in reactions:
            rxn = network[rxn_id]
            for side, other in (('reactants', 'products'), ('products', 'reactants')):
                if set(rxn[side]) <= present and not set(rxn[other]) <= present:
                    present |= set(rxn[other])
                    changed = True
    return present


def calculate_equilibrium(c0_M: dict, T_K: float, *, reactions=('R1', 'R2', 'R3', 'R4'), network: dict = None,
                          species_db: dict = None, species_shifts_eV: dict = None, mu_eV: dict = None,
                          tol: float = 1e-8, max_iter: int = 200) -> dict:
    """Composition at equilibrium [mol/L] of a tube whose composition at mixing is c0_M ({species: mol/L}).

    reactions: those that reach equilibrium; every reaction of the network that follows them is in equilibrium too
    (list_equilibrium_reactions). species_shifts_eV ({species: eV}) moves the species free energies (a solvation
    energy moved from its computed value); mu_eV ({species: eV}) gives them directly instead (G_sol + shift, for a
    caller that computes them once). Returns {species: mol/L} for every species of c0_M and of the reactions.
    """
    net = network or get_active_network()
    active = list_equilibrium_reactions(reactions, net)
    c0 = {sp: float(c) for sp, c in c0_M.items()}
    present = _reachable(c0, net, active)
    running = [r for r in active if set(net[r]['reactants']) | set(net[r]['products']) <= present]
    out = {sp: c0.get(sp, 0.0) for sp in set(c0) | set(list_species({r: net[r] for r in active}))}
    for sp in out:
        if sp not in present:
            out[sp] = 0.0
    if not running:
        return out
    species = list_species({r: net[r] for r in running})
    S = build_stoichiometric_matrix({r: net[r] for r in running}, species)
    A = null_space(S.T)                                         # the amounts the reactions conserve
    if mu_eV is None:
        shifts = species_shifts_eV or {}
        mu_eV = {sp: calculate_solution_gibbs(sp, T_K, species_db=species_db)['G_sol_eV'] + shifts.get(sp, 0.0)
                 for sp in species}
    RT_eV = R_SI * T_K / 1000.0 / EV_TO_KJ_MOL
    mu = np.array([mu_eV[sp] for sp in species]) / RT_eV
    mu = mu - A @ (A.T @ mu)                                    # the part a conserved amount absorbs is dropped
    b = A.T @ np.array([c0.get(sp, 0.0) for sp in species])
    guess = np.log(np.array([c0.get(sp, 0.0) for sp in species]) + 1e-3)
    lam = np.linalg.lstsq(A, guess + mu, rcond=None)[0]

    def dual(lam):
        lnC = np.minimum(-mu + A @ lam, 700.0)
        C = np.exp(lnC)
        return C.sum() - lam @ b, C

    value, C = dual(lam)
    scale = max(np.abs(b).max(), 1e-300)
    for _ in range(max_iter):
        gradient = A.T @ C - b
        if np.abs(gradient).max() <= tol * scale:
            break
        hessian = A.T @ (C[:, None] * A)
        step = np.linalg.lstsq(hessian, -gradient, rcond=None)[0]
        if -gradient @ step <= 1e-16 * scale:                  # Newton decrement at round-off (EC at 7 M sets the scale)
            break
        t = 1.0
        while t > 1e-12:
            trial, C_trial = dual(lam + t * step)
            if trial <= value + 1e-4 * t * (gradient @ step):
                break
            t *= 0.5
        if t <= 1e-12:                                         # no descent left: the gradient is at round-off
            break
        lam, value, C = lam + t * step, trial, C_trial
    else:
        raise RuntimeError(f'equilibrium not reached in {max_iter} Newton steps (|gradient| {np.abs(gradient).max():.3g} M)')
    out.update({sp: float(c) for sp, c in zip(species, C)})
    return out
