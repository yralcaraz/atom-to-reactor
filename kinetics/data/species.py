"""Block 1 — Species database built from the dataset snapshot.

One record per snapshot species, keyed by `pipeline_id`:

    name, mass_g_mol          molar mass from the stored geometry
    G_gas_Eh, G_gas_eV        ωB97M-V/def2-TZVPD Gibbs energy (RRHO, 298.15 K, 1 bar)
    H_gas_eV                  enthalpy at 298.15 K
    E_scf_eV                  electronic SCF energy (qRRHO starts from it)
    dE_solv_eV                MACE MD solvation energy in EC (an energy, not a free energy); None if absent
    dE_solv_sigma_eV          standard error of dE_solv (snapshot `uncertainty_kjmol`, see thermo/uncertainty.py)
    sigma_rot                 rotational symmetry number (placeholder 1, Finding 2)

The snapshot has no frequencies or moments of inertia, so `frequencies_cm1` / `moments_amu_A2` are absent
and the qRRHO route falls back to the stored Gibbs energy.
"""

from copy import deepcopy
from functools import lru_cache

from kinetics.constants import HARTREE_TO_EV, KJ_MOL_TO_EV
from kinetics.data.snapshot import DEFAULT_SNAPSHOT_PATH, calculate_molar_mass, get_dataset_records, get_solvation_records


def load_species_database(snapshot_path: str = DEFAULT_SNAPSHOT_PATH) -> dict:
    """Species database as an independent copy (safe to modify, e.g. to add frequencies)."""
    return deepcopy(_cached_species_database(snapshot_path))


def resolve_species_database(species_db: dict = None) -> dict:
    """`species_db` if given, else the shared read-only default database."""
    return species_db if species_db is not None else _cached_species_database(DEFAULT_SNAPSHOT_PATH)


@lru_cache(maxsize=4)
def _cached_species_database(snapshot_path: str) -> dict:
    solvation = get_solvation_records(snapshot_path)
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
            'sigma_rot': 1,
        }
    return db
