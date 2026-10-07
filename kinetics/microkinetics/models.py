"""Model registry: named, documented combinations of thermochemistry and kinetics.

A ModelSpec fixes everything a reactor needs to turn the network into rate constants
(thermo_mode, kinetic_model, family parameters) and records its status and assumptions, so that
notebooks can compare several models and then select one by name. Register new models in MODELS.
A model may carry its own network (a family split per reaction) and corrections to the computed
free energies (species shifts); the registered models carry neither.

Source: Y. Alcaraz Galván; barrier windows from Gogoi et al., J. Phys. Chem. C 2024, 128, 1654
"""

from dataclasses import dataclass, field, replace

import pandas as pd

from kinetics.microkinetics.barriers import REVERSAL_INVARIANT
from kinetics.microkinetics.parameters import (
    FAMILY_BEP_PARAMETERS, LEVEL1_PARAMETERS, REFERENCE_BEP_PARAMETERS, resolve_family_params,
)
from kinetics.microkinetics.rates import KINETIC_MODELS, calculate_network_rates
from kinetics.data.network import NETWORK
from kinetics.thermo.uncertainty import build_shifted_species_database

MODEL_STATUSES = ('reference', 'sensitivity', 'provisional', 'fitted')


@dataclass(frozen=True, eq=False)
class ModelSpec:
    name: str
    label: str
    kinetic_model: str
    family_params: dict
    thermo_mode: str = 'wb97mv'
    status: str = 'reference'
    assumptions: tuple = field(default_factory=tuple)
    network: dict = None                                   # own network (e.g. a family split per reaction)
    species_shifts_eV: dict = field(default_factory=dict)  # corrections to the free energies in EC, per species

    def __post_init__(self):
        if self.kinetic_model not in KINETIC_MODELS:
            raise ValueError(f"Unknown kinetic_model '{self.kinetic_model}'")
        if self.status not in MODEL_STATUSES:
            raise ValueError(f"status must be one of {MODEL_STATUSES}")

    @property
    def rate_options(self) -> dict:
        """Keyword arguments for calculate_rate_constants / calculate_network_rates."""
        return {'kinetic_model': self.kinetic_model, 'family_params': self.family_params,
                'thermo_mode': self.thermo_mode}

    @property
    def barrier_shapes(self) -> list:
        shape = KINETIC_MODELS[self.kinetic_model][0]
        if shape is not None:
            return [shape]
        return sorted({p.get('shape', 'marcus') for p in self.family_params.values()})

    def calculate_rates(self, T_K: float, **options) -> pd.DataFrame:
        """Rate constants of every reaction at T_K (options: network, species_db, viscosity_Pa_s).

        A model with its own network uses it unless one is passed; its species shifts are applied on top of
        the species database in use.
        """
        if options.get('network') is None and self.network is not None:
            options['network'] = self.network
        if self.species_shifts_eV:
            options['species_db'] = build_shifted_species_database(self.species_shifts_eV, options.get('species_db'))
        return calculate_network_rates(T_K, **self.rate_options, **options)

    def with_barrier(self, g_eV: float, families: list = None, name: str = None) -> 'ModelSpec':
        """Copy with the intrinsic barrier of `families` (default: every entry) set to g_eV."""
        params = {fam: dict(p) for fam, p in self.family_params.items()}
        for fam in families or list(params):
            params.setdefault(fam, {}).pop('E0_eV', None)
            params[fam]['g_eV'] = g_eV
        return replace(self, name=name or f"{self.name}@{g_eV:.2f}", family_params=params,
                       label=f"{self.label} [g = {g_eV:.2f} eV]")

    def with_family_params(self, family: str, name: str = None, **updates) -> 'ModelSpec':
        """Copy with `updates` (e.g. g_eV=1.2, dS_act_J_mol_K=-80) applied to one family.

        A family without its own entry starts from the set's 'default' entry.
        """
        params = {fam: dict(p) for fam, p in self.family_params.items()}
        entry = params.setdefault(family, dict(params.get('default', {})))
        if 'g_eV' in updates:
            entry.pop('E0_eV', None)
        entry.update(updates)
        tag = ', '.join(f'{k}={v:g}' for k, v in updates.items())
        return replace(self, name=name or f"{self.name}[{family}: {tag}]", family_params=params,
                       label=f"{self.label} [{family}: {tag}]")


MODELS = {
    'reference_bep': ModelSpec(
        name='reference_bep',
        label='Reference: capped BEP, E0 = 1.15 eV for all reactions',
        kinetic_model='bep_eyring',
        family_params=REFERENCE_BEP_PARAMETERS,
        status='reference',
        assumptions=(
            'One E0 for all 9 reactions, chosen so that TMSPA depletion spreads over the 20-80 °C holds '
            '(a design choice, not a fit to data)',
            'Capped BEP: every exergonic step has ΔG‡ = E0, so ΔG_rxn does not rank their rates',
            'Barrier depends on the direction a reaction is written (Finding 15)',
        ),
    ),
    'family_bep': ModelSpec(
        name='family_bep',
        label='Family BEP: E0 = 0.80 eV, solvent attack 1.30 eV',
        kinetic_model='bep_eyring',
        family_params=FAMILY_BEP_PARAMETERS,
        status='sensitivity',
        assumptions=(
            'E0 = 0.80 eV for hydrolysis, condensation and transfer is a placeholder (tank_model.ipynb)',
            'Solvent attack E0 = 1.30 eV is derived from the 80 °C protocol of Gogoi et al. 2024',
            'Condensation at 0.80 eV contradicts Gogoi 2024 (Finding 16)',
            'Barrier depends on the direction a reaction is written (Finding 15)',
        ),
    ),
    'family_marcus': ModelSpec(
        name='family_marcus',
        label='Family Marcus: g = 0.80 eV, solvent attack 1.30 eV',
        kinetic_model='marcus_eyring',
        family_params=FAMILY_BEP_PARAMETERS,
        status='sensitivity',
        assumptions=(
            'Same intrinsic barriers as family_bep; λ = 4g',
            'Exergonic steps get lower barriers than under the BEP cap',
        ),
    ),
    'level1': ModelSpec(
        name='level1',
        label='Level 1: Marcus with recalibrated family barriers',
        kinetic_model='level1',
        family_params=LEVEL1_PARAMETERS,
        status='provisional',
        assumptions=(
            'Condensation g = 1.30 eV and solvent attack g = 1.32 eV are derived from Gogoi et al. 2024',
            'Hydrolysis and transfer g = 0.80 eV are placeholders (Finding 18): '
            'TMSPA is consumed within the first 20 °C hold',
        ),
    ),
}


# Names under which a model appears in result files written before it was renamed
LEGACY_MODEL_NAMES = {'peter_reference': 'reference_bep'}


def get_model(model) -> ModelSpec:
    """ModelSpec from a registry name (a ModelSpec is returned unchanged; a legacy name resolves to the current one)."""
    if isinstance(model, ModelSpec):
        return model
    try:
        return MODELS[LEGACY_MODEL_NAMES.get(model, model)]
    except KeyError:
        raise KeyError(f"Unknown model '{model}'. Registered: {list(MODELS)}") from None


def describe_models(models=None) -> pd.DataFrame:
    """Overview table: kinetics, intrinsic barriers, reversal invariance, status and assumptions."""
    rows = {}
    for m in (models or list(MODELS)):
        spec = get_model(m)
        barriers = ', '.join(f"{fam} {p.get('g_eV', p.get('E0_eV')):.2f}" for fam, p in spec.family_params.items())
        rows[spec.name] = {
            'label': spec.label,
            'thermo_mode': spec.thermo_mode,
            'barrier shape': ', '.join(spec.barrier_shapes),
            'g (eV)': barriers,
            'reversal invariant': all(REVERSAL_INVARIANT[s] for s in spec.barrier_shapes),
            'status': spec.status,
            'assumptions': '; '.join(spec.assumptions),
        }
    return pd.DataFrame.from_dict(rows, orient='index')


def tabulate_family_barriers(models=None, network: dict = None) -> pd.DataFrame:
    """Intrinsic barrier g [eV] (E0 of the capped BEP) that each model applies to each reaction family.

    Rows are families in network order with their member reactions; the barrier shape and α are listed
    per model in describe_models.
    """
    net = network or NETWORK
    families = {}
    for rxn_id, rxn in net.items():
        families.setdefault(rxn['class'], []).append(rxn_id)
    rows = {}
    for family, members in families.items():
        row = {'reactions': ', '.join(members)}
        for m in (models or list(MODELS)):
            spec = get_model(m)
            fallback = KINETIC_MODELS[spec.kinetic_model][1]
            row[spec.name] = resolve_family_params(family, spec.family_params, fallback)['g_eV']
        rows[family] = row
    return pd.DataFrame.from_dict(rows, orient='index')
