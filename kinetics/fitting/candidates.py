"""Block 14E — Candidate structures: what is free in each model that the first fit compares.

A structure says which barriers are separate parameters, whether the reaction energies of the phosphate
ladder may move within their computed uncertainty, and whether part of the added water may be unavailable.
It turns a parameter vector into a self-contained ModelSpec (own network, own species shifts).

    M0-BEP      one barrier for all nine reactions, capped BEP (the structure of `peter_reference`)
    M0-Marcus   one barrier for all nine reactions, Marcus
    M1          one Marcus barrier per family (the structure of `level1`)
    M1-split    M1 with R1 separate from R2 and R3 (R2 and R3 share one barrier)
    M3          M1 with ΔG_rxn of R1–R4 free, constrained by their computed covariance
    M3-split    M1-split with the same freed energies
    '+W'        any of them with a common available fraction of the added water (M4 of notebook 03)

Freed energies are parameters d_j = ΔG_rxn,j − computed value for R1–R4. They are applied as species shifts,
so R5–R7 follow (transfer = hydrolysis + condensation) and R8, R9 keep their computed values. Their
constraint is Gaussian with the covariance of the MD standard errors (kinetics.thermo).

The search box of each barrier is DECLARED: it is narrower than the 0.60–1.70 eV scanned in notebook 03 where
no observation allows the excluded part and the solver crawls there (solvent attack below 1.10 eV).

Classification: PUBLIC (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván
"""

from dataclasses import dataclass, field, replace

import numpy as np
from scipy.linalg import cholesky, solve_triangular

from kinetics.data.network import NETWORK
from kinetics.fitting.residuals import FREE_AGE_BOUNDS_H
from kinetics.microkinetics.models import get_model
from kinetics.thermo.reaction import calculate_network_thermo
from kinetics.thermo.uncertainty import calculate_solvation_covariance, calculate_species_shifts

BARRIER_BOX_EV = {
    'all reactions': (0.60, 1.70),
    'hydrolysis': (0.80, 1.70), 'hydrolysis_R1': (0.80, 1.70), 'hydrolysis_R23': (0.80, 1.70),
    'transfer': (0.60, 1.50), 'condensation': (1.00, 1.70), 'solvent_attack': (1.10, 1.60),
}
ENERGY_BOX_EV = (-0.75, 0.75)            # about three standard errors of the computed ΔG_rxn
WATER_BOX = (0.05, 1.0)
FREED_REACTIONS = ('R1', 'R2', 'R3', 'R4')
SPLIT_HYDROLYSIS = {'hydrolysis_R1': ('R1',), 'hydrolysis_R23': ('R2', 'R3')}
FAMILY_ORDER = ('hydrolysis', 'hydrolysis_R1', 'hydrolysis_R23', 'transfer', 'condensation', 'solvent_attack')


@dataclass(frozen=True)
class FitParameter:
    name: str
    kind: str               # 'barrier' | 'energy' | 'water' | 'age'
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
    shape: str = None                            # barrier shape that replaces the base model's (smooth models only)
    parent: str = None                           # the structure this one extends (nesting)

    # --- parameters -------------------------------------------------------------------------------
    def barrier_classes(self) -> list:
        if self.shared_barrier:
            return ['all reactions']
        families = []
        for rxn_id, rxn in self.network().items():
            if rxn['class'] not in families:
                families.append(rxn['class'])
        return sorted(families, key=FAMILY_ORDER.index)

    def parameters(self, samples: list = ()) -> list:
        """The free parameters, in the order of the parameter vector. A sample whose age is a fit parameter
        (free scenario) adds log10(age in hours)."""
        out = [FitParameter(f'g {c}', 'barrier', c, *BARRIER_BOX_EV[c], 'eV') for c in self.barrier_classes()]
        out += [FitParameter(f'dG {r}', 'energy', r, *ENERGY_BOX_EV, 'eV') for r in self.freed]
        if self.water:
            out.append(FitParameter('water fraction', 'water', 'H2O', *WATER_BOX, '1'))
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
        """The network with the split classes (the default network when nothing is split)."""
        if not self.split:
            return NETWORK
        net = {rxn_id: dict(rxn) for rxn_id, rxn in NETWORK.items()}
        for new_class, members in self.split.items():
            for rxn_id in members:
                net[rxn_id]['class'] = new_class
        return net

    def build(self, theta: dict) -> dict:
        """Everything an evaluation needs from a parameter set {name: value}.

        Returns 'model' (a ModelSpec that carries its network and species shifts), 'ages' {sample: hours},
        'water_fraction', and for the freed energies 'prior' (whitened residuals) and 'prior_rows' (one row per
        reaction: shift, standard error, z).
        """
        spec = get_model(self.base_model)
        params = {fam: dict(p) for fam, p in spec.family_params.items()}
        if self.shape is not None:
            for entry in params.values():
                entry['shape'] = self.shape
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
                        network=self.network() if self.split else None, species_shifts_eV=shifts)
        ages = {name[len('log10 age '):]: float(10.0 ** value) for name, value in theta.items()
                if name.startswith('log10 age ')}
        return {'model': model, 'ages': ages, 'water_fraction': float(theta.get('water fraction', 1.0)),
                'prior': prior, 'prior_rows': prior_rows}

    def with_water(self) -> 'FitStructure':
        """The same structure with a free available fraction of the added water."""
        return replace(self, name=f'{self.name}+W', label=f'{self.label}; part of the added water unavailable',
                       water=True, parent=self.name)

    def with_shape(self, shape: str) -> 'FitStructure':
        """The same structure with another smooth barrier shape (e.g. 'agmon_levine')."""
        return replace(self, name=f'{self.name} ({shape})', label=f'{self.label}; barrier shape {shape}', shape=shape)


_ENERGY_TOOLS = {}


def _energy_tools(freed: tuple, prior_scale: float = 1.0):
    """(species-shift matrix, species, Cholesky factor of the covariance, standard errors, computed ΔG) for the
    freed reactions; the covariance is multiplied by prior_scale."""
    key = (freed, prior_scale)
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

    temperatures_K: {rxn_id: T}; default 295.65 K (22.5 °C) for R1–R3 and R5–R7, 353.15 K for R4, R8, R9.
    With ΔS‡ = 0 the barrier does not depend on T; the temperature only says where the value is anchored.
    """
    model = structure.build(theta)['model']
    temperatures_K = temperatures_K or {r: (353.15 if r in ('R4', 'R8', 'R9') else 295.65) for r in NETWORK}
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
    'M0-BEP': FitStructure('M0-BEP', 'one barrier for all reactions, capped BEP', 'peter_reference', shared_barrier=True),
    'M0-Marcus': FitStructure('M0-Marcus', 'one barrier for all reactions, Marcus', 'level1', shared_barrier=True),
    'M1': _M1,
    'M1-split': _M1_SPLIT,
    'M3': replace(_M1, name='M3', label='one Marcus barrier per family; ΔG_rxn of R1–R4 free within their covariance',
                  freed=FREED_REACTIONS, parent='M1'),
    'M3-split': replace(_M1_SPLIT, name='M3-split',
                        label='R1 separate from R2 and R3; ΔG_rxn of R1–R4 free within their covariance',
                        freed=FREED_REACTIONS, parent='M1-split'),
}


def embed_parent_theta(structure, parent_theta: dict) -> dict:
    """Parameter set of `structure` that reproduces a fit of the structure it extends.

    A split class takes the barrier of its family, every family takes a shared barrier, freed energies start
    at their computed values (shift 0) and all the water is available. Parameters the structure does not have
    in common with its parent and that cannot be derived (ages) are left out.
    """
    structure = get_structure(structure)
    out = {}
    for p in structure.parameters():
        if p.name in parent_theta:
            out[p.name] = parent_theta[p.name]
        elif p.kind == 'barrier':
            family = NETWORK[structure.split[p.target][0]]['class'] if p.target in structure.split else p.target
            for candidate in (f'g {family}', 'g all reactions'):
                if candidate in parent_theta:
                    out[p.name] = parent_theta[candidate]
                    break
        elif p.kind == 'energy':
            out[p.name] = 0.0
        elif p.kind == 'water':
            out[p.name] = 1.0
    out.update({k: v for k, v in parent_theta.items() if k.startswith('log10 age ')})
    return out


def get_structure(structure) -> FitStructure:
    """FitStructure from a registry name; 'NAME+W' adds the available-water fraction; a FitStructure passes through."""
    if isinstance(structure, FitStructure):
        return structure
    if structure.endswith('+W'):
        return get_structure(structure[:-2]).with_water()
    try:
        return FIT_STRUCTURES[structure]
    except KeyError:
        raise KeyError(f"Unknown structure '{structure}'. Registered: {list(FIT_STRUCTURES)}") from None
