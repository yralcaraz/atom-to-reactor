"""Mass-action ODE engine shared by the batch (Block 7) and protocol (Block 10) reactors.

    dC/dt = S · (r_f - r_r),   r_f,j = k_f,j ∏_i C_i^ν_ij(reactants),   r_r,j = k_r,j ∏_i C_i^ν_ij(products)

Buffered species (the EC solvent) keep their state value but enter the rates at a fixed concentration.
The analytic Jacobian is passed to the implicit solver, which avoids finite-difference overflow for
very fast models.

Source: Y. Alcaraz Galván
"""

import numpy as np
from scipy.integrate import solve_ivp

from kinetics.data.network import build_order_matrices, build_stoichiometric_matrix, list_species
from kinetics.data.snapshot import count_elements

IMPLICIT_METHODS = ('Radau', 'BDF', 'LSODA')


class IntegrationBudgetExceeded(RuntimeError):
    """An integration needed more right-hand-side calls than its budget (see MassActionSystem.integrate)."""


class MassActionSystem:
    """Right-hand side, Jacobian and integrator of a network over a fixed list of tracked species."""

    def __init__(self, network: dict, species: list, buffered_M: dict = None):
        untracked = [sp for sp in list_species(network) if sp not in species]
        if untracked:
            raise KeyError(f"Network species not tracked by the reactor: {untracked}")
        self.species = list(species)
        self.idx = {sp: i for i, sp in enumerate(self.species)}
        self.orders_f, self.orders_r = build_order_matrices(network, self.species)
        self.buffered_M = dict(buffered_M or {})
        self._buffered_idx = np.array([self.idx[sp] for sp in self.buffered_M], dtype=int)
        self._buffered_val = np.array(list(self.buffered_M.values()), dtype=float)
        self.S = build_stoichiometric_matrix(network, self.species)
        self.S[self._buffered_idx, :] = 0.0
        # (species index, order) per reaction, for the Jacobian
        self._terms_f = [[(i, a) for i, a in enumerate(col) if a > 0] for col in self.orders_f.T]
        self._terms_r = [[(i, a) for i, a in enumerate(col) if a > 0] for col in self.orders_r.T]

    def _rate_concentrations(self, C: np.ndarray) -> np.ndarray:
        c = np.maximum(C, 0.0)
        c[self._buffered_idx] = self._buffered_val
        return c

    def calculate_net_rates(self, C: np.ndarray, k_f: np.ndarray, k_r: np.ndarray) -> np.ndarray:
        """Net rate of every reaction [M/s]."""
        c = self._rate_concentrations(C)[:, None]
        return k_f * np.prod(c ** self.orders_f, axis=0) - k_r * np.prod(c ** self.orders_r, axis=0)

    def rhs(self, t: float, C: np.ndarray, k_f: np.ndarray, k_r: np.ndarray) -> np.ndarray:
        return self.S @ self.calculate_net_rates(C, k_f, k_r)

    def jacobian(self, t: float, C: np.ndarray, k_f: np.ndarray, k_r: np.ndarray) -> np.ndarray:
        c = self._rate_concentrations(C)
        dr = np.zeros((len(k_f), len(c)))
        for j in range(len(k_f)):
            for terms, k, sign in ((self._terms_f[j], k_f[j], 1.0), (self._terms_r[j], k_r[j], -1.0)):
                for i, a in terms:
                    others = np.prod([c[m] ** b for m, b in terms if m != i])
                    dr[j, i] += sign * k * a * c[i] ** (a - 1) * others
        dr[:, self._buffered_idx] = 0.0
        dr[:, C < 0.0] = 0.0   # clamped species do not affect the rates
        return self.S @ dr

    def integrate(self, c0: np.ndarray, k_f: np.ndarray, k_r: np.ndarray, t_eval: np.ndarray, *,
                  method: str = 'Radau', rtol: float = 1e-8, atol: float = 1e-12, max_rhs_calls: int = None):
        """solve_ivp from t = 0 to t_eval[-1].

        max_rhs_calls: budget of right-hand-side calls; IntegrationBudgetExceeded is raised beyond it. It keeps a
        search over parameters from stalling where the solver crawls (None: no limit, the default).
        """
        rhs = self.rhs
        if max_rhs_calls is not None:
            calls = [0]

            def rhs(t, C, k_f, k_r):
                calls[0] += 1
                if calls[0] > max_rhs_calls:
                    raise IntegrationBudgetExceeded(f'more than {max_rhs_calls} right-hand-side calls')
                return self.rhs(t, C, k_f, k_r)
        return solve_ivp(
            rhs, (0.0, float(t_eval[-1])), np.asarray(c0, dtype=float), method=method, t_eval=t_eval,
            rtol=rtol, atol=atol, args=(np.asarray(k_f), np.asarray(k_r)),
            jac=self.jacobian if method in IMPLICIT_METHODS else None,
        )


def calculate_element_totals(C_M: np.ndarray, species: list, elements=('Si', 'P')) -> dict:
    """{element: total concentration of its atoms [M] over time}, atom counts from the snapshot geometries."""
    return {el: sum(count_elements(sp)[el] * C_M[i] for i, sp in enumerate(species)) for el in elements}
