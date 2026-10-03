"""Read-only access to the Tank dataset API snapshot (`data/tank_api_snapshot.json`).

The snapshot is the single data source of the pipeline (see `data/README.md`). It is parsed once
and cached; the returned structures are shared, so callers must not mutate them.

Classification: PUBLIC (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván; reads P. Broqvist, Tank dataset snapshot (unpublished)
"""

import json
import os
from collections import Counter
from functools import lru_cache

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_SNAPSHOT_PATH = os.path.join(_REPO_ROOT, "data", "tank_api_snapshot.json")

# Standard atomic weights [g/mol] for the elements present in the snapshot
ATOMIC_MASS_G_MOL = {
    "H": 1.008, "Li": 6.94, "C": 12.011, "N": 14.007, "O": 15.999,
    "F": 18.998, "Si": 28.085, "P": 30.974, "S": 32.06,
}


@lru_cache(maxsize=4)
def load_snapshot(path: str = DEFAULT_SNAPSHOT_PATH) -> dict:
    """Parsed snapshot JSON (cached per path)."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Snapshot not found: {path}")
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


@lru_cache(maxsize=4)
def get_dataset_records(path: str = DEFAULT_SNAPSHOT_PATH) -> dict:
    """{pipeline_id: gas-phase DFT record}."""
    return {d["pipeline_id"]: d for d in load_snapshot(path)["datasets"]}


@lru_cache(maxsize=4)
def get_solvation_records(path: str = DEFAULT_SNAPSHOT_PATH) -> dict:
    """{species name: MACE MD solvation record}."""
    return {s["raw_metadata"]["name"]: s for s in load_snapshot(path)["solvation"]}


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
