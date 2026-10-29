"""Block 14E — Candidate structures: what is free in each model that the first fit compares.

Names are MX-suffix. X says which reaction energies are free: M0 and M1 none (computed values), M3 those of
R1–R4, M5 those of every independent reaction (M2 is not used; M4 is the '+W' of notebook 03). The suffix says how
the barriers are grouped: none, one per family (M0: one for all); '-split', R1 apart from R2 and R3; '-first', the
first step of TMSPA apart in every family (R1 and R5).

A structure says which barriers are separate parameters, whether the reaction energies of the phosphate
ladder may move within their computed uncertainty, and whether part of the added water may be unavailable.
It turns a parameter vector into a self-contained ModelSpec (own network, own species shifts).

    M0-BEP      one barrier for all nine reactions, capped BEP (the structure of `reference_bep`)
    M0-Marcus   one barrier for all nine reactions, Marcus
    M1          one Marcus barrier per family (the structure of `level1`)
    M1-split    M1 with R1 separate from R2 and R3 (R2 and R3 share one barrier)
    M3          M1 with ΔG_rxn of R1–R4 free, constrained by their computed covariance
    M3-split    M1-split with the same freed energies
    M3-first    M3 with the first step of TMSPA (no P–OH yet) apart in every family: R1 (with N15) in
                hydrolysis and R5 in transfer; R2–R3 and R6–R7 share one class each (sixth fit)
    M5          M1 with the reaction energies of every independent reaction free: R1–R4, R8 and R9 (seventh fit)
    M5-first    the same barriers as M3-first, with the reaction energies of every independent reaction free:
                R1–R4, R8 and R9 (sixth fit)
    '+W'        any of them with a common available fraction of the added water (M4 of notebook 03)
    '+BEP'      any of them with the barrier rule of the reference model instead of Marcus (seventh fit, see with_bep_rule)

Two kinds of parameter come from the sample set and not from the structure: the age of the sample heated after
its first spectrum (free scenario), and the amount at mixing of a species a recipe does not fix ('c0 H2O TMSPa
alone': the water of the tube mixed without water, third fit).

Freed energies are parameters d_j = ΔG_rxn,j − computed value for R1–R4. They are applied as species shifts,
so R5–R7 follow (transfer = hydrolysis + condensation) and R8, R9 keep their computed values. Their
constraint is Gaussian with the covariance of the MD standard errors (kinetics.thermo).

The network is the active variant (kinetics.data.network; fifth fit: R1–R9 and N10–N15). A step that
'follows' a reaction (N12 follows R2, N15 follows R1) takes the barrier class of that reaction in each structure,
so it adds no parameter; a class of FIXED_BARRIERS_EV ('exchange': N10, N11, N13, N14) has a declared barrier
and is not fitted. Their ΔG come from the species energies, freed shifts included.

Sixth fit: the box of a freed energy can be set for a whole run to ±k standard errors of its computed value
(set_energy_box_sigma, fit_config.json 'energy_box_sigma'); the Gaussian constraint stays inside it.

Seventh fit: the barrier rule of the reference model, ΔG‡ = E₀ for a downhill step and E₀ + |ΔG_rxn| for an
uphill one, is the capped BEP of `reference_bep` with α = 1 in place of 0.5. With α = 1 it no longer depends on
the direction a step is written (the reverse barrier is E₀ when the forward one is E₀ + ΔG_rxn), so M0-BEP
(α = 0.5) is not its special case. Any structure takes it through the suffix '+BEP'.

The search box of each barrier is DECLARED: it is narrower than the 0.60–1.70 eV scanned in notebook 03 where
no observation allows the excluded part and the solver crawls there (solvent attack below 1.10 eV).

Source: Y. Alcaraz Galván
"""

from dataclasses import dataclass, field, replace

import numpy as np
from scipy.linalg import cholesky, solve_triangular

from kinetics.data.medium import get_active_medium
from kinetics.data.network import NETWORK, get_active_network
from kinetics.fitting.residuals import FREE_AGE_BOUNDS_H
from kinetics.microkinetics.models import get_model
from kinetics.thermo.reaction import calculate_network_thermo
from kinetics.thermo.uncertainty import calculate_solvation_covariance, calculate_species_shifts

# Fifth fit: the box of the first hydrolysis reaches 1.80 eV. N15 (TMSPA + H2O → MMSPA + HMDSO) shares its barrier
# and starts the reaction about seven times faster than R1 at the same g, so g_R1 has to rise by about 0.05 eV; the
# fourth fit's 7 day reading already sat at 1.635 eV
BARRIER_BOX_EV = {
    'all reactions': (0.60, 1.70),
    'hydrolysis': (0.80, 1.70), 'hydrolysis_R1': (0.80, 1.80), 'hydrolysis_R23': (0.80, 1.70),
    'transfer': (0.60, 1.50), 'transfer_R5': (0.60, 1.50), 'transfer_R67': (0.60, 1.50),
    'condensation': (1.00, 1.70), 'solvent_attack': (1.10, 1.60),
}
ENERGY_BOX_EV = (-0.75, 0.75)            # about three standard errors of the computed ΔG_rxn
WATER_BOX = (0.05, 1.0)
RECIPE_PREFIX = 'c0 '                    # 'c0 <species> <sample>': amount at mixing of a species the recipe does not fix
FREED_REACTIONS = ('R1', 'R2', 'R3', 'R4')
FREED_INDEPENDENT = ('R1', 'R2', 'R3', 'R4', 'R8', 'R9')     # every reaction of R1–R9 that is not a sum of others
SPLIT_HYDROLYSIS = {'hydrolysis_R1': ('R1',), 'hydrolysis_R23': ('R2', 'R3')}
SPLIT_FIRST_STEP = {**SPLIT_HYDROLYSIS, 'transfer_R5': ('R5',), 'transfer_R67': ('R6', 'R7')}
FAMILY_ORDER = ('hydrolysis', 'hydrolysis_R1', 'hydrolysis_R23', 'transfer', 'transfer_R5', 'transfer_R67',
                'condensation', 'solvent_attack')
# Barriers that are declared, not fitted. 'exchange': a silyl group moving between two phosphates or two glycols
# (N10, N11, N13, N14). DECLARED: the value the reference model gives every step it has no data for (its class
# "other"); at room temperature it takes about a year in our tubes
FIXED_BARRIERS_EV = {'exchange': 1.15}
BEP_SUFFIX = '+BEP'                # structure name suffix of the barrier rule of the reference model (with_bep_rule)
BEP_ALPHA = 1.0                    # its slope: an uphill step costs E₀ plus all of its ΔG_rxn


@dataclass(frozen=True)
class FitParameter:
    name: str
    kind: str               # 'barrier' | 'energy' | 'water' | 'recipe' | 'age'
    target: str             # family (or 'all reactions'), reaction id, 'H2O' or sample name
    lower: float
    upper: float
    unit: str


@dataclass(frozen=True)
class FitStructure:
    name: str
    label: str
    base_model: str                              # registry name the structure starts from
    shared_barrier: bool = False                 # one barrier for every reaction
    split: dict = field(default_factory=dict)    # {new class: reactions} taken out of their family
    freed: tuple = ()                            # reactions whose ΔG_rxn is free
    water: bool = False                          # a common available fraction of the added water
    shape: str = None                            # barrier shape that replaces the base model's
    shape_params: dict = field(default_factory=dict)   # parameters of that shape ('+BEP': {'alpha': 1.0})
    parent: str = None                           # the structure this one extends (nesting)

    # --- parameters -------------------------------------------------------------------------------
    def barrier_classes(self) -> list:
        if self.shared_barrier:
            return ['all reactions']
        families = []
        for rxn_id, rxn in self.network().items():
            if rxn['class'] not in families and rxn['class'] not in FIXED_BARRIERS_EV:
                families.append(rxn['class'])
        return sorted(families, key=FAMILY_ORDER.index)

    def parameters(self, samples: list = ()) -> list:
        """The free parameters, in the order of the parameter vector. A sample with a species its recipe does
        not fix adds that amount at mixing [M]; a sample whose age is a fit parameter (free scenario) adds
        log10(age in hours)."""
        out = [FitParameter(f'g {c}', 'barrier', c, *BARRIER_BOX_EV[c], 'eV') for c in self.barrier_classes()]
        k = get_energy_box_sigma()
        boxes = [(-k * s, k * s) for s in _energy_tools(self.freed)[3]] if k and self.freed else [ENERGY_BOX_EV] * len(self.freed)
        out += [FitParameter(f'dG {r}', 'energy', r, *box, 'eV') for r, box in zip(self.freed, boxes)]
        if self.water:
            out.append(FitParameter('water fraction', 'water', 'H2O', *WATER_BOX, '1'))
        for sample in samples:
            for species, (low, high) in sample.get('unknown_M', {}).items():
                out.append(FitParameter(f"{RECIPE_PREFIX}{species} {sample['name']}", 'recipe', sample['name'], low, high, 'M'))
        for sample in samples:
            if sample['age_mode'] == 'parameter':
                out.append(FitParameter(f"log10 age {sample['name']}", 'age', sample['name'],
                                        *np.log10(FREE_AGE_BOUNDS_H), 'log10 h'))
        return out

    @property
    def n_structure_parameters(self) -> int:
        return len(self.parameters())

    # --- model ------------------------------------------------------------------------------------
    def network(self) -> dict:
        """The active network with the split classes, and each step that follows a reaction in that reaction's class
        (the active network as it is when nothing changes)."""
        active = get_active_network()
        if not self.split and not any('follows' in rxn for rxn in active.values()):
            return active
        net = {rxn_id: dict(rxn) for rxn_id, rxn in active.items()}
        for new_class, members in self.split.items():
            for rxn_id in members:
                net[rxn_id]['class'] = new_class
        for rxn in net.values():
            if 'follows' in rxn:
                rxn['class'] = net[rxn['follows']]['class']
        return net

    def build(self, theta: dict) -> dict:
        """Everything an evaluation needs from a parameter set {name: value}.

        Returns 'model' (a ModelSpec that carries its network and species shifts), 'ages' {sample: hours},
        'water_fraction', 'added_M' {sample: {species: mol/L at mixing on top of the recipe}}, and for the freed
        energies 'prior' (whitened residuals) and 'prior_rows' (one row per reaction: shift, standard error, z).
        """
        spec = get_model(self.base_model)
        params = {fam: dict(p) for fam, p in spec.family_params.items()}
        if self.shape is not None:
            for entry in params.values():
                entry['shape'] = self.shape
                entry.update(self.shape_params)
        if self.shared_barrier:
            for entry in params.values():
                entry.pop('E0_eV', None)
                entry['g_eV'] = float(theta['g all reactions'])
        else:
            for new_class, members in self.split.items():
                family = NETWORK[members[0]]['class']
                params[new_class] = dict(params.get(family, params.get('default', {})))
            for c in self.barrier_classes():
                entry = params.setdefault(c, dict(params.get('default', {})))
                entry.pop('E0_eV', None)
                entry['g_eV'] = float(theta[f'g {c}'])
        network = self.network()
        for c in sorted({rxn['class'] for rxn in network.values()} & set(FIXED_BARRIERS_EV)):
            entry = params[c] = dict(params.get('default', {}))
            entry.pop('E0_eV', None)
            entry['g_eV'] = FIXED_BARRIERS_EV[c]
            if self.shape is not None:
                entry['shape'] = self.shape
                entry.update(self.shape_params)
        shifts, prior, prior_rows = {}, np.zeros(0), []
        if self.freed:
            d = np.array([float(theta[f'dG {r}']) for r in self.freed])
            matrix, species, chol, sigma, computed = _energy_tools(self.freed)
            shifts = {sp: float(v) for sp, v in zip(species, matrix @ d) if abs(v) > 1e-15}
            prior = solve_triangular(chol, d, lower=True)
            prior_rows = [{'sample': 'computed energies', 'quantity': f'ΔG_rxn({r}) − computed', 'kind': 'prior',
                           't_h': np.nan, 'age_h': np.nan, 'measured': 0.0, 'predicted': d[j], 'sigma': sigma[j],
                           'z': d[j] / sigma[j]} for j, r in enumerate(self.freed)]
        model = replace(spec, name=f'{self.name} (fit)', label=self.label, family_params=params,
                        network=network, species_shifts_eV=shifts)
        ages = {name[len('log10 age '):]: float(10.0 ** value) for name, value in theta.items()
                if name.startswith('log10 age ')}
        added = {}
        for name, value in theta.items():
            if name.startswith(RECIPE_PREFIX):
                species, sample = name[len(RECIPE_PREFIX):].split(' ', 1)
                added.setdefault(sample, {})[species] = float(value)
        return {'model': model, 'ages': ages, 'water_fraction': float(theta.get('water fraction', 1.0)),
                'added_M': added, 'prior': prior, 'prior_rows': prior_rows}

    def with_water(self) -> 'FitStructure':
        """The same structure with a free available fraction of the added water."""
        return replace(self, name=f'{self.name}+W', label=f'{self.label}; part of the added water unavailable',
                       water=True, parent=self.name)

    def with_shape(self, shape: str) -> 'FitStructure':
        """The same structure with another smooth barrier shape (e.g. 'agmon_levine')."""
        return replace(self, name=f'{self.name} ({shape})', label=f'{self.label}; barrier shape {shape}', shape=shape)

    def with_bep_rule(self) -> 'FitStructure':
        """The same structure with the barrier rule of the reference model in every class, the declared exchange steps included:
        capped BEP with α = 1, ΔG‡ = max(E₀, E₀ + ΔG_rxn). Its barriers are fitted as E₀ under the same names
        ('g <class>'). It extends the '+BEP' version of this structure's parent."""
        return replace(self, name=f'{self.name}{BEP_SUFFIX}',
                       label=f"{self.label}; capped BEP barriers (α = 1) instead of Marcus",
                       shape='bep_cap', shape_params={'alpha': BEP_ALPHA},
                       parent=f'{self.parent}{BEP_SUFFIX}' if self.parent else None)


_ENERGY_TOOLS = {}
_ENERGY_BOX_SIGMA = [None]


def get_energy_box_sigma():
    """k of the run's energy box (±k standard errors of the computed value), or None: the fixed ENERGY_BOX_EV."""
    return _ENERGY_BOX_SIGMA[0]


def set_energy_box_sigma(k: float = None):
    _ENERGY_BOX_SIGMA[0] = None if k is None else float(k)
    return _ENERGY_BOX_SIGMA[0]


def _energy_tools(freed: tuple, prior_scale: float = 1.0):
    """(species-shift matrix, species, Cholesky factor of the covariance, standard errors, computed ΔG) for the
    freed reactions in the active solvation medium; the covariance is multiplied by prior_scale."""
    key = (get_active_medium(), freed, prior_scale)
    if key not in _ENERGY_TOOLS:
        unit = [calculate_species_shifts({r: (1.0 if r == rxn else 0.0) for r in freed}) for rxn in freed]
        species = sorted({sp for shifts in unit for sp in shifts})
        matrix = np.array([[shifts.get(sp, 0.0) for shifts in unit] for sp in species])
        cov = prior_scale * calculate_solvation_covariance(list(freed)).to_numpy()
        computed = calculate_network_thermo(298.15)['dG_rxn_eV'].loc[list(freed)].to_numpy()
        _ENERGY_TOOLS[key] = (matrix, species, cholesky(cov, lower=True), np.sqrt(np.diag(cov)), computed)
    return _ENERGY_TOOLS[key]


def scale_energy_prior(structure: FitStructure, built: dict, theta: dict, prior_scale: float) -> dict:
    """`built` with the constraint on the freed energies recomputed for a covariance multiplied by prior_scale."""
    if not structure.freed or prior_scale == 1.0:
        return built
    d = np.array([float(theta[f'dG {r}']) for r in structure.freed])
    _, _, chol, sigma, _ = _energy_tools(structure.freed, prior_scale)
    rows = [{**row, 'sigma': sigma[j], 'z': d[j] / sigma[j]} for j, row in enumerate(built['prior_rows'])]
    return {**built, 'prior': solve_triangular(chol, d, lower=True), 'prior_rows': rows}


def tabulate_reaction_barriers(structure: FitStructure, theta: dict, temperatures_K: dict = None) -> dict:
    """ΔG‡ [eV] of every reaction at the temperature where its family is observed, and its ΔG_rxn.

    temperatures_K: {rxn_id: T}; default 353.15 K for R4, R8, R9 and N11 (seen at 80 °C), 295.65 K (22.5 °C) for
    every other reaction of the structure's network.
    With ΔS‡ = 0 the barrier does not depend on T; the temperature only says where the value is anchored.
    """
    model = structure.build(theta)['model']
    temperatures_K = temperatures_K or {r: (353.15 if r in ('R4', 'R8', 'R9', 'N11') else 295.65)
                                        for r in structure.network()}
    out = {}
    for T in sorted(set(temperatures_K.values())):
        rates = model.calculate_rates(T)
        for rxn_id, T_rxn in temperatures_K.items():
            if T_rxn == T:
                out[rxn_id] = {'T_K': T, 'dG_barrier_eV': float(rates.loc[rxn_id, 'dG_barrier_f_eV']),
                               'dG_rxn_eV': float(rates.loc[rxn_id, 'dG_rxn_eV']), 'k_f': float(rates.loc[rxn_id, 'k_f'])}
    return out


_M1 = FitStructure('M1', 'one Marcus barrier per family', 'level1', parent='M0-Marcus')
_M1_SPLIT = FitStructure('M1-split', 'Marcus barrier per family, R1 separate from R2 and R3', 'level1',
                         split=SPLIT_HYDROLYSIS, parent='M1')
FIT_STRUCTURES = {
    'M0-BEP': FitStructure('M0-BEP', 'one barrier for all reactions, capped BEP', 'reference_bep', shared_barrier=True),
    'M0-Marcus': FitStructure('M0-Marcus', 'one barrier for all reactions, Marcus', 'level1', shared_barrier=True),
    'M1': _M1,
    'M1-split': _M1_SPLIT,
    'M3': replace(_M1, name='M3', label='one Marcus barrier per family; ΔG_rxn of R1–R4 free within their covariance',
                  freed=FREED_REACTIONS, parent='M1'),
    'M3-split': replace(_M1_SPLIT, name='M3-split',
                        label='R1 separate from R2 and R3; ΔG_rxn of R1–R4 free within their covariance',
                        freed=FREED_REACTIONS, parent='M1-split'),
    'M3-first': replace(_M1, name='M3-first', split=SPLIT_FIRST_STEP, freed=FREED_REACTIONS, parent='M3-split',
                       label='first step of TMSPA apart in hydrolysis (R1) and transfer (R5); ΔG_rxn of R1–R4 free '
                             'within their covariance'),
    'M5': replace(_M1, name='M5', freed=FREED_INDEPENDENT, parent='M1',
                  label='one Marcus barrier per family; ΔG_rxn of every independent reaction (R1–R4, R8, R9) free '
                        'within their covariance'),
    'M5-first': replace(_M1, name='M5-first', split=SPLIT_FIRST_STEP, freed=FREED_INDEPENDENT, parent='M3-first',
                          label='first step of TMSPA apart in hydrolysis (R1) and transfer (R5); ΔG_rxn of every '
                                'independent reaction (R1–R4, R8, R9) free within their covariance'),
}


def get_structure(structure) -> FitStructure:
    """FitStructure from a registry name; 'NAME+W' adds the available-water fraction, 'NAME+BEP' the capped
    BEP rule (α = 1); a FitStructure passes through."""
    if isinstance(structure, FitStructure):
        return structure
    if structure.endswith(BEP_SUFFIX):
        return get_structure(structure[:-len(BEP_SUFFIX)]).with_bep_rule()
    if structure.endswith('+W'):
        return get_structure(structure[:-2]).with_water()
    try:
        return FIT_STRUCTURES[structure]
    except KeyError:
        raise KeyError(f"Unknown structure '{structure}'. Registered: {list(FIT_STRUCTURES)}") from None
