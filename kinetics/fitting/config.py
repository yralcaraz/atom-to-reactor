"""Block 14H — Setup of one fit run, stored with its results (`fit_config.json` in the results directory).

Every stage of a run reads the setup at its start, and so does every notebook that reads the run's results: the
solvation medium the reaction energies are computed in, the readings of the unknown ages that are fitted, the
paper observations added to the lab samples and the changes to the synthetic recovery check ('recovery': truths
replaced per case, and the search of the fits for every case; see scripts/run_fit.recovery_cases). A results directory without the file is a run of the default setup,
that of the first three fits: neat EC from the 2026-09-24 snapshot, the four readings and E1.

    DEFAULT_FIT_CONFIG  the setup of a directory without fit_config.json
    load_fit_config     the setup stored in a results directory
    apply_fit_config    load it and make its medium the active one (kinetics.data.medium); returns the setup
    check_fit_medium    raise if a stored fit was computed in another medium than the active one

Source: Y. Alcaraz Galván
"""

import json
from pathlib import Path

from kinetics.data.medium import DEFAULT_MEDIUM, SolvationMedium, get_active_medium, set_active_medium
from kinetics.fitting.residuals import PAPER_OBSERVATIONS, SCENARIOS

FIT_CONFIG_FILE = 'fit_config.json'
DEFAULT_FIT_CONFIG = {'solvation': DEFAULT_MEDIUM.to_dict(), 'scenarios': list(SCENARIOS),
                      'paper_observations': list(PAPER_OBSERVATIONS), 'recovery': {}}


def load_fit_config(results_dir) -> dict:
    """The setup of the run whose results are in `results_dir`: fit_config.json over the default, with the medium
    as a SolvationMedium under 'medium'."""
    path = Path(results_dir) / FIT_CONFIG_FILE
    stored = json.loads(path.read_text(encoding='utf-8')) if path.exists() else {}
    config = {**DEFAULT_FIT_CONFIG, **stored}
    unknown = set(config['scenarios']) - set(SCENARIOS)
    if unknown:
        raise ValueError(f'{path}: unknown scenarios {sorted(unknown)}')
    config['medium'] = SolvationMedium.from_dict(config['solvation'])
    config['source'] = str(path) if path.exists() else 'default (no fit_config.json)'
    return config


def apply_fit_config(results_dir) -> dict:
    """load_fit_config, and its medium made the active one. Call it before computing anything from a run's fits."""
    config = load_fit_config(results_dir)
    set_active_medium(config['medium'])
    return config


def check_fit_medium(fit: dict) -> None:
    """A stored fit is run forward only in the medium it was fitted in (fits without the record: the default)."""
    stored = SolvationMedium.from_dict(fit['medium']) if fit.get('medium') else DEFAULT_MEDIUM
    if stored != get_active_medium():
        raise ValueError(f"this fit was computed in '{stored.label()}', the active medium is "
                         f"'{get_active_medium().label()}': call apply_fit_config(<its results directory>) first")
