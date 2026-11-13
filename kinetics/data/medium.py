"""Block 1 — Solvation medium: which MD solvation energies the species database takes, and how they are combined.

A snapshot holds MD solvation energies in neat solvents only. The lab tubes are EC/DEC 1:1 by volume, so a mixture
is built from the neat values with mole fractions (from the volume fractions, densities and molar masses):

    neat       one solvent: its MD values as stored
    linear        ΔE_mix = Σ x_s ΔE_s, standard errors with the same weights (the ideal-mixture limit; no
                  preference of a solute for either solvent)
    preferential  ΔE_mix = −kT ln Σ x_s exp(−ΔE_s/kT) at 298.15 K: each solute is surrounded mostly by the solvent
                  that solvates it better, so the value lies at or below the more negative of the two; standard
                  errors with the Boltzmann weights w_s = x_s exp(−ΔE_s/kT) / Σ (the same rule as the reference model)

The two rules bracket a real mixture: neither is a measurement.

The source of the solvation energies (fit_thermo_rev01 §5.1):

    md         MACE MD ΔE_solv, its standard errors as σ (the default)
    cosmors    openCOSMO-RS ΔG_solv of the snapshot's `cosmors` table; σ = |MD − COSMO-RS| of each species in the same
               medium (it has no sampling error: the disagreement of the two methods is its uncertainty)

The medium in use is one per process (the active medium): `resolve_species_database(None)` and the solvation
standard errors of `kinetics.thermo` follow it. The default is neat EC from the 2026-09-24 snapshot, the medium of
notebooks 01–03 and of fits 1–3. A fit run sets its medium once, at the start (kinetics.fitting.config).

    SolvationMedium     snapshot file, volume fractions and mixing rule; hashable, so it keys the caches
    get_active_medium   the medium in use
    set_active_medium   change it (None: back to the default)

Source: Y. Alcaraz Galván; mixing rules as in the reference model
"""

import math
from dataclasses import dataclass

from kinetics.constants import EV_TO_J, KB_SI, T_STD_K
from kinetics.data.snapshot import DEFAULT_SNAPSHOT_FILE, calculate_molar_mass
from kinetics.paths import data_path

# Densities of the neat solvents [g/mL] at room temperature (EC as in kinetics/reactor/protocol.py; DEC as in
# notebook 03 §5.4)
SOLVENT_DENSITY_G_ML = {'EC': 1.321, 'DEC': 0.975}
MIXING_RULES = ('linear', 'preferential')
KT_MIX_EV = KB_SI * T_STD_K / EV_TO_J        # kT of the preferential rule [eV]
SOLVATION_SOURCES = ('md', 'cosmors')


@dataclass(frozen=True)
class SolvationMedium:
    """snapshot: file name in the data folder; volume_fractions: ((solvent, fraction), ...) summing to one."""
    snapshot: str = DEFAULT_SNAPSHOT_FILE
    volume_fractions: tuple = (('EC', 1.0),)
    rule: str = 'linear'
    source: str = 'md'

    def __post_init__(self):
        if abs(sum(f for _, f in self.volume_fractions) - 1.0) > 1e-9:
            raise ValueError(f'volume fractions must sum to one: {self.volume_fractions}')
        if self.rule not in MIXING_RULES:
            raise ValueError(f"mixing rule must be one of {MIXING_RULES}, not '{self.rule}'")
        if self.source not in SOLVATION_SOURCES:
            raise ValueError(f"solvation source must be one of {SOLVATION_SOURCES}, not '{self.source}'")

    @property
    def path(self) -> str:
        return data_path('data', self.snapshot)

    @property
    def solvents(self) -> tuple:
        return tuple(s for s, _ in self.volume_fractions)

    @property
    def neat(self) -> bool:
        return len(self.volume_fractions) == 1

    def mole_fractions(self) -> dict:
        """{solvent: mole fraction} from the volume fractions, the densities and the molar masses of the geometries."""
        moles = {s: f * SOLVENT_DENSITY_G_ML[s] / calculate_molar_mass(s, self.path) for s, f in self.volume_fractions}
        total = sum(moles.values())
        return {s: n / total for s, n in moles.items()}

    def label(self) -> str:
        date = self.snapshot.replace('tank_api_snapshot_', '').replace('.json', '')
        method = 'MD' if self.source == 'md' else 'openCOSMO-RS'
        if self.neat:
            return f'neat {self.solvents[0]}, {method} solvation of the {date} snapshot'
        volumes = ':'.join(f'{f:.2f}' for _, f in self.volume_fractions)
        moles = ', '.join(f'x_{s} {x:.3f}' for s, x in self.mole_fractions().items())
        return f"{':'.join(self.solvents)} {volumes} by volume ({moles}), {self.rule} mixing, {method} solvation of the {date} snapshot"

    def to_dict(self) -> dict:
        record = {'snapshot': self.snapshot, 'volume_fractions': dict(self.volume_fractions), 'rule': self.rule}
        return record if self.source == 'md' else {**record, 'source': self.source}   # MD fits keep their old record

    @classmethod
    def from_dict(cls, record: dict) -> 'SolvationMedium':
        return cls(snapshot=record['snapshot'], volume_fractions=tuple(record['volume_fractions'].items()),
                   rule=record.get('rule', 'linear'), source=record.get('source', 'md'))


def mix_solvation(values: dict, mole_fractions: dict, rule: str) -> tuple:
    """(ΔE_mix, {solvent: weight}) of one solute from its neat-solvent values {solvent: eV}; the weights carry the
    standard errors, σ_mix² = Σ (w_s σ_s)²."""
    if rule == 'linear':
        return sum(mole_fractions[s] * v for s, v in values.items()), dict(mole_fractions)
    if rule != 'preferential':
        raise ValueError(f"mixing rule must be one of {MIXING_RULES}, not '{rule}'")
    low = min(values.values())
    boltzmann = {s: mole_fractions[s] * math.exp(-(v - low) / KT_MIX_EV) for s, v in values.items()}
    total = sum(boltzmann.values())
    return low - KT_MIX_EV * math.log(total), {s: b / total for s, b in boltzmann.items()}


DEFAULT_MEDIUM = SolvationMedium()
_ACTIVE = [DEFAULT_MEDIUM]


def get_active_medium() -> SolvationMedium:
    return _ACTIVE[0]


def set_active_medium(medium: SolvationMedium = None) -> SolvationMedium:
    """Make `medium` the one in use (None: the default) and return it. Caches are keyed by the medium, so nothing
    computed in another medium is reused."""
    _ACTIVE[0] = medium or DEFAULT_MEDIUM
    return _ACTIVE[0]
