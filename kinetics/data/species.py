"""Block 1 — Species database built from the dataset snapshot.

One record per snapshot species, keyed by `pipeline_id`:

    name, mass_g_mol          molar mass from the stored geometry
    G_gas_Eh, G_gas_eV        ωB97M-V/def2-TZVPD Gibbs energy (RRHO, 298.15 K, 1 bar)
    H_gas_eV                  enthalpy at 298.15 K
    E_scf_eV                  electronic SCF energy (qRRHO starts from it)
    dE_solv_eV                MACE MD solvation energy in the chosen solvent (an energy, not a free energy); None if absent.
                              In a COSMO-RS medium (fit 6.2): the COSMO-RS solvation free energy, mixed by the medium's rule
    solvent, solvation_T_K    solvent of that energy and temperature of its MD run (None if absent; for a mixture
                              the solvents joined by ':' and no temperature)
    dE_solv_sigma_eV          standard error of dE_solv (snapshot `uncertainty_kjmol`, see thermo/uncertainty.py)
    sigma_rot                 rotational symmetry number (placeholder 1, Finding 2)

The snapshot has no frequencies or moments of inertia, so `frequencies_cm1` / `moments_amu_A2` are absent
and the qRRHO route falls back to the stored Gibbs energy.

The shared default database is that of the active solvation medium (kinetics/data/medium.py): neat EC from the
2026-09-24 snapshot unless a run sets another. In a mixture a species has a solvation energy only if every
solvent of the mixture has an MD record of it.

Source: Y. Alcaraz Galván; DFT and MD data by P. Broqvist
"""

import math
from copy import deepcopy
from functools import lru_cache

from kinetics.constants import HARTREE_TO_EV, KJ_MOL_TO_EV
from kinetics.data.medium import SolvationMedium, get_active_medium
from kinetics.data.snapshot import (
    DEFAULT_SNAPSHOT_PATH, DEFAULT_SOLVENT, calculate_molar_mass, get_cosmors_records, get_dataset_records,
    get_solvation_records,
)


def load_species_database(snapshot_path: str = DEFAULT_SNAPSHOT_PATH, solvent: str = DEFAULT_SOLVENT) -> dict:
    """Species database as an independent copy (safe to modify, e.g. to add frequencies)."""
    return deepcopy(_cached_species_database(snapshot_path, solvent))


def load_medium_species_database(medium: SolvationMedium) -> dict:
    """Species database of a solvation medium (neat solvent or mixture), as an independent copy."""
    return deepcopy(_cached_medium_database(medium))


def resolve_species_database(species_db: dict = None) -> dict:
    """`species_db` if given, else the shared read-only database of the active solvation medium."""
    return species_db if species_db is not None else _cached_medium_database(get_active_medium())


@lru_cache(maxsize=4)
def _cached_medium_database(medium: SolvationMedium) -> dict:
    if medium.method == 'cosmors':
        return _cosmors_database(medium)
    if medium.neat:
        return _cached_species_database(medium.path, medium.solvents[0])
    x = medium.mole_fractions()
    parts = {s: _cached_species_database(medium.path, s) for s in medium.solvents}
    db = {}
    for name, record in parts[medium.solvents[0]].items():
        values = [(x[s], parts[s][name]['dE_solv_eV'], parts[s][name]['dE_solv_sigma_eV']) for s in medium.solvents]
        complete = all(v is not None for _, v, _ in values)
        db[name] = {
            **record,
            'dE_solv_eV': sum(w * v for w, v, _ in values) if complete else None,
            'dE_solv_sigma_eV': math.sqrt(sum((w * e) ** 2 for w, _, e in values)) if complete else None,
            'solvent': ':'.join(medium.solvents) if complete else None,
            'solvation_T_K': None,
        }
    return db


def _cosmors_database(medium: SolvationMedium) -> dict:
    """The database of a COSMO-RS medium: gas-phase records as in the MD one, solvation from the COSMO-RS table,
    mixed linearly with mole fractions. No standard error is stored for it (dE_solv_sigma_eV None)."""
    x = {medium.solvents[0]: 1.0} if medium.neat else medium.mole_fractions()
    tables = {s: get_cosmors_records(medium.path, s) for s in medium.solvents}
    db = {}
    for name, record in _cached_species_database(medium.path, medium.solvents[0]).items():
        complete = all(name in tables[s] for s in medium.solvents)
        db[name] = {
            **record,
            'dE_solv_eV': sum(x[s] * tables[s][name] for s in medium.solvents) if complete else None,
            'dE_solv_sigma_eV': None,
            'solvent': ':'.join(medium.solvents) if complete else None,
            'solvation_T_K': 298.15 if complete else None,
        }
    return db


@lru_cache(maxsize=4)
def _cached_species_database(snapshot_path: str, solvent: str) -> dict:
    solvation = get_solvation_records(snapshot_path, solvent)
    db = {}
    for name, record in get_dataset_records(snapshot_path).items():
        solv = solvation.get(name)
        db[name] = {
            'name': name,
            'mass_g_mol': calculate_molar_mass(name, snapshot_path),
            'G_gas_Eh': record['gibbs_eh'],
            'G_gas_eV': record['gibbs_eh'] * HARTREE_TO_EV,
            'H_gas_eV': record['enthalpy_eh'] * HARTREE_TO_EV,
            'E_scf_eV': record['energy_scf_eh'] * HARTREE_TO_EV,
            'dE_solv_eV': solv['delta_e_solv_kjmol'] * KJ_MOL_TO_EV if solv else None,
            'dE_solv_sigma_eV': solv['uncertainty_kjmol'] * KJ_MOL_TO_EV if solv else None,
            'solvent': solvent if solv else None,
            'solvation_T_K': solv['temperature_k'] if solv else None,
            'sigma_rot': 1,
        }
    return db
