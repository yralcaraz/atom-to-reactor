"""Block 1 — Species database built from the dataset snapshot.

One record per snapshot species, keyed by `pipeline_id`:

    name, mass_g_mol          molar mass from the stored geometry
    G_gas_Eh, G_gas_eV        ωB97M-V/def2-TZVPD Gibbs energy (RRHO, 298.15 K, 1 bar)
    H_gas_eV                  enthalpy at 298.15 K
    E_scf_eV                  electronic SCF energy (qRRHO starts from it)
    dE_solv_eV                MACE MD solvation energy in the chosen solvent (an energy, not a free energy); None if absent
    solvent, solvation_T_K    solvent of that energy and temperature of its MD run (None if absent; for a mixture
                              the solvents joined by ':' and no temperature)
    dE_solv_sigma_eV          standard error of dE_solv (snapshot `uncertainty_kjmol`, see thermo/uncertainty.py)
                              (source 'cosmors': dE_solv is the openCOSMO-RS ΔG_solv and σ = |MD − COSMO-RS|)
    sigma_rot                 rotational symmetry number (placeholder 1, Finding 2)

The snapshot has no frequencies or moments of inertia, so `frequencies_cm1` / `moments_amu_A2` are absent
and the qRRHO route falls back to the stored Gibbs energy.

The shared default database is that of the active solvation medium (kinetics/data/medium.py): neat EC from the
2026-09-24 snapshot unless a run sets another. In a mixture a species has a solvation energy only if every
solvent of the mixture has an MD record of it.

Source: Y. Alcaraz Galván; DFT and MD data from the Tank dataset
"""

import math
from functools import lru_cache

from kinetics.constants import HARTREE_TO_EV, KCAL_EV, KJ_MOL_TO_EV
from kinetics.data.medium import SolvationMedium, get_active_medium, mix_solvation
from kinetics.data.snapshot import (
    calculate_molar_mass, get_cosmors_records, get_dataset_records, get_solvation_records,
)


def resolve_species_database(species_db: dict = None) -> dict:
    """`species_db` if given, else the shared read-only database of the active solvation medium."""
    return species_db if species_db is not None else _cached_medium_database(get_active_medium())


def get_medium_database(medium: SolvationMedium) -> dict:
    """The shared read-only database of one solvation medium (not necessarily the active one)."""
    return _cached_medium_database(medium)


@lru_cache(maxsize=4)
def _cached_medium_database(medium: SolvationMedium) -> dict:
    if medium.source == 'cosmors':
        return _cosmors_database(medium)
    if medium.neat:
        return _cached_species_database(medium.path, medium.solvents[0])
    x = medium.mole_fractions()
    parts = {s: _cached_species_database(medium.path, s) for s in medium.solvents}
    db = {}
    for name, record in parts[medium.solvents[0]].items():
        values = {s: parts[s][name]['dE_solv_eV'] for s in medium.solvents}
        complete = all(v is not None for v in values.values())
        value, weights = mix_solvation(values, x, medium.rule) if complete else (None, None)
        db[name] = {
            **record,
            'dE_solv_eV': value,
            'dE_solv_sigma_eV': math.sqrt(sum((weights[s] * parts[s][name]['dE_solv_sigma_eV']) ** 2 for s in medium.solvents))
                                if complete else None,
            'solvent': ':'.join(medium.solvents) if complete else None,
            'solvation_T_K': None,
        }
    return db


def _cosmors_database(medium: SolvationMedium) -> dict:
    """The MD database of the same medium with dE_solv from openCOSMO-RS (mixed by the medium's rule) and
    σ = |MD − COSMO-RS| of each species; a species without a COSMO-RS value in every solvent has none."""
    md = _cached_medium_database(SolvationMedium(medium.snapshot, medium.volume_fractions, medium.rule))
    x = medium.mole_fractions() if not medium.neat else {medium.solvents[0]: 1.0}
    cosmo = {s: get_cosmors_records(medium.path, s) for s in medium.solvents}
    db = {}
    for name, record in md.items():
        complete = all(name in cosmo[s] for s in medium.solvents)
        value = mix_solvation({s: cosmo[s][name]['dg_solv_kcal_mol'] * KCAL_EV for s in medium.solvents}, x, medium.rule)[0] \
            if complete else None
        md_value = record['dE_solv_eV']
        db[name] = {
            **record,
            'dE_solv_eV': value,
            'dE_solv_sigma_eV': abs(md_value - value) if value is not None and md_value is not None else None,
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
