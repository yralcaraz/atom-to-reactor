"""Model registry: named, documented combinations of thermochemistry and kinetics.

A ModelSpec fixes everything a reactor needs to turn the network into rate constants
(thermo_mode, kinetic_model, family parameters) and records its status and assumptions, so that
notebooks can compare several models and then select one by name. Register new models in MODELS.
A model may carry its own network (a family split per reaction) and corrections to the computed
free energies (species shifts); the registered models carry neither. The registry holds the two models the fit
structures start from: reference_bep (M0-BEP) and level1 (every other structure); family_bep and family_marcus are
in the branch M5-full.

Source: Y. Alcaraz Galván; barrier windows from Gogoi et al., J. Phys. Chem. C 2024, 128, 1654
"""

from dataclasses import dataclass, field

import pandas as pd

from kinetics.microkinetics.parameters import (
    LEVEL1_PARAMETERS, REFERENCE_BEP_PARAMETERS,
)
from kinetics.microkinetics.rates import KINETIC_MODELS, calculate_network_rates
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


def get_model(model) -> ModelSpec:
    """ModelSpec from a registry name (a ModelSpec is returned unchanged)."""
    if isinstance(model, ModelSpec):
        return model
    try:
        return MODELS[model]
    except KeyError:
        raise KeyError(f"Unknown model '{model}'. Registered: {list(MODELS)}") from None


