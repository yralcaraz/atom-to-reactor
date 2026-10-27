"""Read-only access to a Tank dataset API snapshot (`tank_api_snapshot_<capture date>.json`).

The snapshot is the single data source of the pipeline (see `data/README.md`): the Tank dataset
(DFT, MD solvation, computed NMR shieldings, geometries). It is read from the data folder (`ATOM_DATA_DIR`,
see `kinetics/paths.py`). It is parsed once and cached; the returned structures are shared, so callers must not mutate them.

Snapshots in use (see `data/README.md`): 2026-09-24 (MD solvation in EC at 298 K, the default; fits 1–3) and 2026-10-08
(MD solvation in EC at 338 K and in DEC at 298 K, openCOSMO-RS table). Since a snapshot may hold the MD solvation in
more than one solvent, solvation records are always read for one solvent.

Source: Y. Alcaraz Galván; DFT and MD data from the Tank dataset
"""

import json
import os
from collections import Counter
from functools import lru_cache

from kinetics.paths import DATA_DIR_VAR, data_path

DEFAULT_SNAPSHOT_FILE = "tank_api_snapshot_2026-09-24.json"
DEFAULT_SNAPSHOT_PATH = data_path("data", DEFAULT_SNAPSHOT_FILE)
SNAPSHOT_2026_10_08_PATH = data_path("data", "tank_api_snapshot_2026-10-08.json")
DEFAULT_SOLVENT = "EC"

# Solvent of an MD solvation record, for records that do not name it in `raw_metadata['solvent']` (2026-09-24 snapshot)
SOLVENT_IDS = {"1,3-dioxolan-2-one": "EC", "diethyl carbonate": "DEC"}

# Standard atomic weights [g/mol] for the elements present in the snapshot
ATOMIC_MASS_G_MOL = {
    "H": 1.008, "Li": 6.94, "C": 12.011, "N": 14.007, "O": 15.999,
    "F": 18.998, "Si": 28.085, "P": 30.974, "S": 32.06,
}


@lru_cache(maxsize=4)
def load_snapshot(path: str = DEFAULT_SNAPSHOT_PATH) -> dict:
    """Parsed snapshot JSON (cached per path)."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Snapshot not found: {path} (set {DATA_DIR_VAR}, "
                                f"see .env.example)")
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


@lru_cache(maxsize=4)
def get_dataset_records(path: str = DEFAULT_SNAPSHOT_PATH) -> dict:
    """{pipeline_id: gas-phase DFT record}."""
    return {d["pipeline_id"]: d for d in load_snapshot(path)["datasets"]}


def get_solvent(record: dict) -> str:
    """Solvent pipeline id ('EC', 'DEC') of one MD solvation record."""
    return record["raw_metadata"].get("solvent") or SOLVENT_IDS[record["solvent_name"]]


def list_solvation_sets(path: str = DEFAULT_SNAPSHOT_PATH) -> dict:
    """{solvent: {'T_K': MD temperatures, 'n_species': record count}} of the MD solvation in a snapshot."""
    sets = {}
    for record in load_snapshot(path)["solvation"]:
        entry = sets.setdefault(get_solvent(record), {"T_K": set(), "n_species": 0})
        entry["T_K"].add(record["temperature_k"])
        entry["n_species"] += 1
    return {solvent: {"T_K": sorted(e["T_K"]), "n_species": e["n_species"]} for solvent, e in sets.items()}


@lru_cache(maxsize=8)
def get_solvation_records(path: str = DEFAULT_SNAPSHOT_PATH, solvent: str = DEFAULT_SOLVENT) -> dict:
    """{species name: MACE MD solvation record} in one solvent.

    Raises if the snapshot has no record in that solvent, or two records of one species in it (e.g. two MD
    temperatures): a species must never take its energy from whichever record comes last.
    """
    records = {}
    for record in load_snapshot(path)["solvation"]:
        if get_solvent(record) != solvent:
            continue
        name = record["raw_metadata"]["name"]
        if name in records:
            raise ValueError(f"Two MD solvation records of '{name}' in {solvent} in {os.path.basename(path)}")
        records[name] = record
    if not records:
        raise KeyError(f"No MD solvation in '{solvent}' in {os.path.basename(path)}: "
                       f"it has {sorted(list_solvation_sets(path))}")
    return records


def get_structure(pipeline_id: str, path: str = DEFAULT_SNAPSHOT_PATH) -> tuple[list, list]:
    """(elements, coordinates_angstrom) of the optimised geometry."""
    uuid = _dataset_uuid(pipeline_id, path)
    record = load_snapshot(path)["structure"][uuid]
    return record["elements"], record["coordinates_angstrom"]


def get_nmr_shieldings(pipeline_id: str, path: str = DEFAULT_SNAPSHOT_PATH) -> list:
    """Per-atom shielding records [{'element', 'atom_index', 'isotropic_ppm', ...}], or [] if absent."""
    record = get_dataset_records(path).get(pipeline_id)
    if record is None or not record.get("has_nmr", False):
        return []
    return load_snapshot(path).get("nmr", {}).get(record["dataset_uuid"], {}).get("nmr_shieldings", [])


@lru_cache(maxsize=64)
def count_elements(pipeline_id: str, path: str = DEFAULT_SNAPSHOT_PATH) -> Counter:
    """Atom count per element, from the stored geometry."""
    elements, _ = get_structure(pipeline_id, path)
    return Counter(elements)


def calculate_molar_mass(pipeline_id: str, path: str = DEFAULT_SNAPSHOT_PATH) -> float:
    """Molar mass [g/mol] from the stored geometry."""
    return sum(ATOMIC_MASS_G_MOL[el] * n for el, n in count_elements(pipeline_id, path).items())


def _dataset_uuid(pipeline_id: str, path: str) -> str:
    try:
        return get_dataset_records(path)[pipeline_id]["dataset_uuid"]
    except KeyError:
        raise KeyError(f"Species '{pipeline_id}' is not in the snapshot") from None
