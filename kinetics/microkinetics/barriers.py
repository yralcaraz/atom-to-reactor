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
'bep_cap' is the legacy Peter Broqvist form and violates A2, A3 and A5; it is kept only to reproduce
earlier results.

Source: Y. Alcaraz Galván
"""

import math

LN2 = math.log(2.0)


def _softplus(t: float) -> float:
    """Overflow-safe ln(1 + e^t)."""
    return max(t, 0.0) + math.log1p(math.exp(-abs(t)))


def _logistic(t: float) -> float:
    """Overflow-safe 1 / (1 + e^-t)."""
    return 0.5 * (1.0 + math.tanh(0.5 * t))


def _check_g(g: float) -> None:
    if not g > 0.0:
        raise ValueError(f"Intrinsic barrier g must be positive, got {g!r} eV")


# ------------------------------------------------------------------------------
# Legacy: Bell-Evans-Polanyi with floor at E0
# ------------------------------------------------------------------------------
def bep_cap(x: float, g: float, alpha: float = 0.5, **_) -> tuple:
    """ΔG‡ = max(g, g + α·x). Legacy form (tank_model.ipynb); direction-dependent."""
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
# Agmon-Levine (Chem. Phys. Lett. 1977, 52, 197)
# ------------------------------------------------------------------------------
def agmon_levine(x: float, g: float, **_) -> tuple:
    """ΔG‡ = x + (g/ln2)·ln[1 + exp(-x·ln2/g)];  α logistic in (0, 1), never inverts."""
    _check_g(g)
    c = g / LN2
    return x + c * _softplus(-x / c), _logistic(x / c)


# ------------------------------------------------------------------------------
# Blowers-Masel (AIChE J. 2000, 46, 2041; piecewise form as in Cantera)
# ------------------------------------------------------------------------------
def _blowers_masel_value(x: float, g: float, w: float) -> float:
    if x <= -4.0 * g:
        return 0.0
    if x >= 4.0 * g:
        return x
    v_p = 2.0 * w * (w + g) / (w - g)
    return (w + x / 2.0) * (v_p - 2.0 * w + x) ** 2 / (v_p ** 2 - 4.0 * w ** 2 + x ** 2)


def blowers_masel(x: float, g: float, w: float = 5.0, **_) -> tuple:
    """Blowers-Masel barrier; w = mean BDE of bonds broken/formed (eV). Tends to Marcus as w → ∞."""
    _check_g(g)
    if not w > g:
        raise ValueError(f"Blowers-Masel requires w > g (w={w!r}, g={g!r})")
    value = _blowers_masel_value(x, g, w)
    if x <= -4.0 * g:
        return value, 0.0
    if x >= 4.0 * g:
        return value, 1.0
    h = 1e-6
    slope = (_blowers_masel_value(x + h, g, w) - _blowers_masel_value(x - h, g, w)) / (2.0 * h)
    return value, slope


# ------------------------------------------------------------------------------
# Two-parabola model with unequal curvatures (α(0) = α0)
# ------------------------------------------------------------------------------
def two_parabola(x: float, g: float, alpha0: float = 0.5, **_) -> tuple:
    """
    Crossing of G_R = k_R q² and G_P = k_P (q-1)² + x with k_R = g/(1-α0)², k_P = g/α0².
    F(0) = g and α(0) = α0 exactly; α0 = 1/2 recovers Marcus. Reversal maps α0 → 1 - α0.
    Outside the real-crossing window the barrier falls back to max(0, x).
    """
    _check_g(g)
    if not 0.0 < alpha0 < 1.0:
        raise ValueError(f"alpha0 must lie in (0, 1), got {alpha0!r}")
    if abs(alpha0 - 0.5) < 1e-12:
        return marcus(x, g)
    k_r = g / (1.0 - alpha0) ** 2
    k_p = g / alpha0 ** 2
    disc = k_r * k_p + (k_r - k_p) * x
    if disc < 0.0:
        return (0.0, 0.0) if x < 0.0 else (x, 1.0)
    q = (-k_p + math.sqrt(disc)) / (k_r - k_p)
    if q <= 0.0:
        return 0.0, 0.0
    if q >= 1.0:
        return x, 1.0
    return k_r * q ** 2, k_r * q / (k_r * q + k_p * (1.0 - q))


# ------------------------------------------------------------------------------
# Registry, reversal and work terms
# ------------------------------------------------------------------------------
# Reversal-invariant shapes satisfy axiom A2; the legacy cap does not
REVERSAL_INVARIANT = {'bep_cap': False, 'marcus': True, 'agmon_levine': True, 'blowers_masel': True,
                      'two_parabola': True}

BARRIER_MODELS = {
    'bep_cap': bep_cap,
    'marcus': marcus,
    'agmon_levine': agmon_levine,
    'blowers_masel': blowers_masel,
    'two_parabola': two_parabola,
}


def reversed_params(params: dict) -> dict:
    """Parameters θ̄ of the reversed reaction: α0 → 1 - α0 and wR ↔ wP."""
    rev = dict(params)
    if 'alpha0' in rev:
        rev['alpha0'] = 1.0 - rev['alpha0']
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


def invert_marcus_barrier(dG_barrier_obs_eV: float, x: float) -> float:
    """Intrinsic barrier g that reproduces an observed barrier under Marcus (closed form)."""
    a = dG_barrier_obs_eV - x / 2.0
    disc = a * a - x * x / 4.0
    if a <= 0.0 or disc < 0.0:
        raise ValueError("Observed barrier is incompatible with the Marcus relation for this ΔG°")
    return 0.5 * (a + math.sqrt(disc))
