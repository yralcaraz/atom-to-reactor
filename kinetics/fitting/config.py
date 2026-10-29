"""Block 14H — Setup of one fit run, stored with its results (`fit_config.json` in the results directory).

Every stage of a run reads the setup at its start, and so does every notebook that reads the run's results: the
solvation medium the reaction energies are computed in, the reaction network (kinetics.data.network variants),
the changes to the samples that hold for the whole run ('sample_overrides', e.g. the heating steps of the 2 %
tube in the fifth fit; keys as in residuals.build_sample_set), the readings of the unknown ages that are fitted,
the sensitivity variants of the night stage ('sensitivities':
names of scripts/run_fit.SENSITIVITIES; None runs those of the fourth fit) and the changes to the synthetic
recovery check ('recovery': truths replaced per case, the search of the fits for every case, or the cases
themselves; see scripts/run_fit.recovery_cases). Sixth fit: the box of the freed energies in standard errors
('energy_box_sigma', kinetics.fitting.candidates), and for scripts/run_fit.py the structures of the run
('structures': main and extended), their search ('effort'), the fits whose parameters start the search
('start_from': an earlier results directory) and the parts of the night and traces stages that run ('night').
Eighth fit: the readings '10 min', '1 h' and '8 h' (the only ones of this branch), and the steps of a fit made one tube at a time ('steps', 'step_variants';
see scripts/run_fit.stage_steps).
A results directory without the file is a run of the default setup: neat EC from the 2026-09-24 snapshot, R1–R9
and the three readings. A setup with paper observations (Gogoi 2024, the first fits) is refused: they are in the
branch M5-full.

    DEFAULT_FIT_CONFIG  the setup of a directory without fit_config.json
    load_fit_config     the setup stored in a results directory
    apply_fit_config    load it and make its medium, network and sample changes the active ones; returns the setup
    check_fit_medium    raise if a stored fit was computed in another medium than the active one
    check_fit_setup     check_fit_medium, and the same for the network and the changes to the samples

Source: Y. Alcaraz Galván
"""

import json
from pathlib import Path

from kinetics.data.medium import DEFAULT_MEDIUM, SolvationMedium, get_active_medium, set_active_medium
from kinetics.data.network import (
    DEFAULT_NETWORK_VARIANT, NETWORK_VARIANTS, get_active_network_variant, set_active_network_variant,
)
from kinetics.fitting.candidates import set_energy_box_sigma
from kinetics.fitting.residuals import (
    DEFAULT_SCENARIOS, SCENARIOS, get_sample_overrides, set_sample_overrides,
)

FIT_CONFIG_FILE = 'fit_config.json'
DEFAULT_FIT_CONFIG = {'solvation': DEFAULT_MEDIUM.to_dict(), 'network': DEFAULT_NETWORK_VARIANT, 'sample_overrides': {},
                      'scenarios': list(DEFAULT_SCENARIOS),
                      'sensitivities': None, 'recovery': {}, 'energy_box_sigma': None, 'structures': None,
                      'effort': {}, 'start_from': None, 'night': {}, 'steps': [], 'step_variants': {}}


def load_fit_config(results_dir) -> dict:
    """The setup of the run whose results are in `results_dir`: fit_config.json over the default, with the medium
    as a SolvationMedium under 'medium'."""
    path = Path(results_dir) / FIT_CONFIG_FILE
    stored = json.loads(path.read_text(encoding='utf-8')) if path.exists() else {}
    config = {**DEFAULT_FIT_CONFIG, **stored}
    unknown = set(config['scenarios']) - set(SCENARIOS)
    if unknown:
        raise ValueError(f'{path}: unknown scenarios {sorted(unknown)} (this branch reads {list(SCENARIOS)})')
    if config.get('paper_observations'):
        raise ValueError(f"{path}: paper observations {config['paper_observations']} are not in this branch")
    if config['network'] not in NETWORK_VARIANTS:
        raise ValueError(f"{path}: unknown network '{config['network']}' (one of {NETWORK_VARIANTS})")
    config['medium'] = SolvationMedium.from_dict(config['solvation'])
    config['source'] = str(path) if path.exists() else 'default (no fit_config.json)'
    return config


def apply_fit_config(results_dir) -> dict:
    """load_fit_config, and its medium, network and changes to the samples made the active ones. Call it before
    computing anything from a run's fits."""
    config = load_fit_config(results_dir)
    set_active_medium(config['medium'])
    set_active_network_variant(config['network'])
    set_sample_overrides(config['sample_overrides'])
    set_energy_box_sigma(config['energy_box_sigma'])
    return config


def check_fit_medium(fit: dict) -> None:
    """A stored fit is run forward only in the medium it was fitted in (fits without the record: the default)."""
    stored = SolvationMedium.from_dict(fit['medium']) if fit.get('medium') else DEFAULT_MEDIUM
    if stored != get_active_medium():
        raise ValueError(f"this fit was computed in '{stored.label()}', the active medium is "
                         f"'{get_active_medium().label()}': call apply_fit_config(<its results directory>) first")


def check_fit_setup(fit: dict) -> None:
    """check_fit_medium, and the same for the network and for the run's changes to the samples: a stored fit is run
    forward with the network it was fitted with and, where it records them (fifth fit on), on top of the same
    run-wide changes. The overrides of a variant (stored under 'variant_overrides') are passed by the caller."""
    check_fit_medium(fit)
    stored = fit.get('network', DEFAULT_NETWORK_VARIANT)
    if stored != get_active_network_variant():
        raise ValueError(f"this fit was computed with the network '{stored}', the active one is "
                         f"'{get_active_network_variant()}': call apply_fit_config(<its results directory>) first")
    if 'sample_overrides' in fit and json.loads(json.dumps(get_sample_overrides())) != fit['sample_overrides']:
        raise ValueError(f"this fit was computed with the run-wide changes to the samples {fit['sample_overrides']}, "
                         f"the active ones are {get_sample_overrides()}: call apply_fit_config(<its results directory>) first")
