"""Block 5 (Level 1) — Free-energy barrier models.

Pure functions mapping a reaction free energy x = ΔG°_rxn (eV) onto a forward activation free energy
ΔG‡_f (eV). No thermodynamics is evaluated here; rate_constants.py supplies x and applies Eyring TST and
detailed balance. Every model returns (dG_barrier_f_eV, alpha_eff), where alpha_eff = ∂ΔG‡/∂x is the local
Brønsted/Leffler coefficient.

Axioms (TFM KIN - DRAFT - Level 1 formulation, §2):
    A2 reversal invariance  F(-x; θ̄) = F(x; θ) - x
    A3 bounds               max(0, x) <= F(x)
    A4 anchoring            F(0) = g
    A6 Leffler bounds       0 <= α(x) <= 1
'bep_cap' is the legacy form of the reference model and violates A2, A3 and A5; it is kept for the reference model
(M0-BEP and the structures '+BEP'). The smooth shapes of Level 1 other than Marcus (Agmon-Levine, Blowers-Masel,
two-parabola) are in the branch M5-full.

Source: Y. Alcaraz Galván
"""


def _check_g(g: float) -> None:
    if not g > 0.0:
        raise ValueError(f"Intrinsic barrier g must be positive, got {g!r} eV")


# ------------------------------------------------------------------------------
# Legacy: Bell-Evans-Polanyi with floor at E0
# ------------------------------------------------------------------------------
def bep_cap(x: float, g: float, alpha: float = 0.5, **_) -> tuple:
    """ΔG‡ = max(g, g + α·x). Legacy form of the reference model; direction-dependent."""
    linear = g + alpha * x
    if linear >= g:
        return linear, alpha
    return g, 0.0


# ------------------------------------------------------------------------------
# Marcus (group-transfer form, λ = 4g)
# ------------------------------------------------------------------------------
def marcus(x: float, g: float, **_) -> tuple:
    """ΔG‡ = (λ/4)(1 + x/λ)² = g + x/2 + x²/(16g);  α = 1/2 + x/(8g)."""
    _check_g(g)
    lam = 4.0 * g
    return (lam / 4.0) * (1.0 + x / lam) ** 2, 0.5 + x / (2.0 * lam)


# ------------------------------------------------------------------------------
# Registry, reversal and work terms
# ------------------------------------------------------------------------------
BARRIER_MODELS = {
    'bep_cap': bep_cap,
    'marcus': marcus,
}


def reversed_params(params: dict) -> dict:
    """Parameters θ̄ of the reversed reaction: wR ↔ wP."""
    rev = dict(params)
    rev['wR_eV'], rev['wP_eV'] = params.get('wP_eV', 0.0), params.get('wR_eV', 0.0)
    return rev


def calculate_barrier(x: float, shape: str, g: float, wR_eV: float = 0.0, wP_eV: float = 0.0, **params) -> tuple:
    """
    Forward barrier with Marcus work terms:  ΔG‡_f = max(wR + F(x - wR + wP), 0, x).
    The floor enforces A3 when a work term is negative; being symmetric under reversal
    (max(0, -x) = max(0, x) - x) it preserves A2. Returns (dG_barrier_f_eV, alpha_eff).
    """
    if shape not in BARRIER_MODELS:
        raise ValueError(f"Unknown barrier model '{shape}'. Supported: {list(BARRIER_MODELS)}")
    y = x - wR_eV + wP_eV
    value, alpha_eff = BARRIER_MODELS[shape](y, g, **params)
    value += wR_eV
    if shape != 'bep_cap' and value < max(0.0, x):
        return max(0.0, x), (1.0 if x > 0.0 else 0.0)
    return value, alpha_eff


