"""Species Database Loader and Default Thermodynamic Properties.

Constructs and caches molecular properties (DFT Gibbs free energy, solvation offset,
molecular weights, and NMR references) from `data/tank_api_snapshot.json` with
graceful fallbacks to benchmark literature values.
"""

import os
import json

HAR2EV = 27.211386245988
KJMOL2EV = 1.0 / 96.48533212331
EV2KCAL = 23.060547830619

# Benchmark molecular weights (g/mol)
BENCHMARK_MASSES = {
    "TMSPA": 314.54, "BMSPA": 242.38, "MMSPA": 170.22, "H3PO4": 98.00,
    "H2O": 18.015, "TMSOH": 90.20, "HMDSO": 162.38, "siloxyl": 162.38,
    "silicone": 236.53, "CH4": 16.04, "EC": 88.062, "VC": 86.05,
    "EO": 44.05, "VO": 42.04, "CO2": 44.01, "EG": 62.07,
    "VG": 60.05, "H2CO3": 62.03, "Ethane": 30.07, "Ethene": 28.05,
    "Ethyne": 26.04, "H2": 2.016, "O2": 31.999, "TMSOEG": 134.25,
    "TMSOVG": 132.23, "TMSOCH3": 104.22, "TMS": 88.22, "TMSOdiEG": 178.30
}

_CACHED_SPECIES_DB = None


def load_default_species_database(snapshot_path: str = None) -> dict:
    """Builds and returns the comprehensive species database."""
    global _CACHED_SPECIES_DB
    if _CACHED_SPECIES_DB is not None:
        return _CACHED_SPECIES_DB

    if snapshot_path is None:
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        candidate_paths = [
            os.path.join(base_dir, "data", "tank_api_snapshot.json"),
            os.path.join(base_dir, "tank_api_snapshot.json"),
        ]
        for p in candidate_paths:
            if os.path.exists(p):
                snapshot_path = p
                break

    db = {}
    if snapshot_path and os.path.exists(snapshot_path):
        with open(snapshot_path, "r", encoding="utf-8") as f:
            snap = json.load(f)

        datasets = {d["pipeline_id"]: d for d in snap.get("datasets", [])}
        solv = {s.get("raw_metadata", {}).get("name", ""): s for s in snap.get("solvation", [])}

        for pid, d in datasets.items():
            g_hartree = d.get("gibbs_eh")
            g_ev = g_hartree * HAR2EV if g_hartree is not None else 0.0
            
            s_entry = solv.get(pid, {})
            de_solv_kjmol = s_entry.get("delta_e_solv_kjmol", 0.0)
            de_solv_ev = de_solv_kjmol * KJMOL2EV
            
            db[pid] = {
                "name": pid,
                "mass_g_mol": BENCHMARK_MASSES.get(pid, 100.0),
                "G_wb97mv_Hartree": g_hartree,
                "G_wb97mv_eV": g_ev,
                "G_B3_Hartree": g_hartree,
                "G_B3_eV": g_ev,
                "E_0K_eV": g_ev,
                "dE_solv_eV": de_solv_ev,
                "dE_solv_kcal_mol": de_solv_ev * EV2KCAL,
                "dE_solv_kJ_mol": de_solv_kjmol,
                "sigma_rot": 1,
            }
            
        # Support alias siloxyl -> HMDSO
        if "HMDSO" in db and "siloxyl" not in db:
            db["siloxyl"] = db["HMDSO"].copy()
            db["siloxyl"]["name"] = "siloxyl"

    _CACHED_SPECIES_DB = db
    return db
