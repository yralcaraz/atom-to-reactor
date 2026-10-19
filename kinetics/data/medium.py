"""Block 1 — Solvation medium: which MD solvation energies the species database takes, and how they are combined.

A snapshot holds MD solvation energies in neat solvents only. The lab tubes are EC/DEC 1:1 by volume, so a mixture
is built from the neat values with mole fractions (from the volume fractions, densities and molar masses):

    neat       one solvent: its MD values as stored
    linear     ΔE_mix = Σ x_s ΔE_s, standard errors with the same weights (the ideal-mixture limit; no
               preference of a solute for either solvent)

The medium in use is one per process (the active medium): `resolve_species_database(None)` and the solvation
standard errors of `kinetics.thermo` follow it. The default is neat EC from the 2026-09-24 snapshot, the medium of
notebooks 01–03 and of fits 1–3. A fit run sets its medium once, at the start (kinetics.fitting.config).

    SolvationMedium     snapshot file, volume fractions and mixing rule; hashable, so it keys the caches
    get_active_medium   the medium in use
    set_active_medium   change it (None: back to the default)

Source: Y. Alcaraz Galván; mixing rules as in P. Broqvist's tmspa_reaction_network_model notebook
"""

from dataclasses import dataclass

from kinetics.data.snapshot import DEFAULT_SNAPSHOT_FILE, calculate_molar_mass
from kinetics.paths import data_path

# Densities of the neat solvents [g/mL] at room temperature (EC as in kinetics/reactor/protocol.py; DEC as in
# notebook 03 §5.4)
SOLVENT_DENSITY_G_ML = {'EC': 1.321, 'DEC': 0.975}
MIXING_RULES = ('linear',)


@dataclass(frozen=True)
class SolvationMedium:
    """snapshot: file name in the data folder; volume_fractions: ((solvent, fraction), ...) summing to one."""
    snapshot: str = DEFAULT_SNAPSHOT_FILE
    volume_fractions: tuple = (('EC', 1.0),)
    rule: str = 'linear'

    def __post_init__(self):
        if abs(sum(f for _, f in self.volume_fractions) - 1.0) > 1e-9:
            raise ValueError(f'volume fractions must sum to one: {self.volume_fractions}')
        if self.rule not in MIXING_RULES:
            raise ValueError(f"mixing rule must be one of {MIXING_RULES}, not '{self.rule}'")

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
        if self.neat:
            return f'neat {self.solvents[0]}, MD solvation of the {date} snapshot'
        volumes = ':'.join(f'{f:.2f}' for _, f in self.volume_fractions)
        moles = ', '.join(f'x_{s} {x:.3f}' for s, x in self.mole_fractions().items())
        return f"{':'.join(self.solvents)} {volumes} by volume ({moles}), {self.rule} mixing, MD solvation of the {date} snapshot"

    def to_dict(self) -> dict:
        return {'snapshot': self.snapshot, 'volume_fractions': dict(self.volume_fractions), 'rule': self.rule}

    @classmethod
    def from_dict(cls, record: dict) -> 'SolvationMedium':
        return cls(snapshot=record['snapshot'], volume_fractions=tuple(record['volume_fractions'].items()),
                   rule=record.get('rule', 'linear'))


DEFAULT_MEDIUM = SolvationMedium()
_ACTIVE = [DEFAULT_MEDIUM]


def get_active_medium() -> SolvationMedium:
    return _ACTIVE[0]


def set_active_medium(medium: SolvationMedium = None) -> SolvationMedium:
    """Make `medium` the one in use (None: the default) and return it. Caches are keyed by the medium, so nothing
    computed in another medium is reused."""
    _ACTIVE[0] = medium or DEFAULT_MEDIUM
    return _ACTIVE[0]
