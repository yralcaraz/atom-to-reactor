"""Block 14H — Setup of one fit run, stored with its results (`fit_config.json` in the results directory).

Every stage of a run, and every notebook that reads its results, reads the setup first. Keys (missing: the default):

    solvation          the solvation medium of the reaction energies (kinetics.data.medium.SolvationMedium)
    network            the reaction network variant (kinetics.data.network.NETWORK_VARIANTS)
    sample_overrides   changes to the samples for the whole run (keys as in residuals.build_sample_set)
    scenarios          the readings of the unknown ages: '10 min', '1 h', '8 h'
    energy_box_sigma   the box of the freed energies, in standard errors (kinetics.fitting.candidates)
    structures         {'main': ..., 'extended': [...]}: the structures scripts/run_fit.py fits
    effort             the search per structure: {name: [screen points, starts, iterations]}
    steps              the steps of a fit made one tube at a time (scripts/run_fit.stage_steps)
    step_variants      other readings of a tube's history, scored after every step
    thermo             the thermodynamic fit (kinetics.fitting.thermo_fit.DEFAULT_THERMO_SETUP, scripts/run_thermo_fit.py)
    species_shifts_eV  fixed corrections to the solvation energies of every model, {species: eV}
    tubes              the tubes a step can fit (default: every tube of the observables with a 'fit' block)
    observables        the lab observables file, relative to the results directory (default: the run's own
                       data_dir/data/lab_observables.json, else the data folder's)

A results directory without the file is a run of the default setup: neat EC from the 2026-09-24 snapshot, R1–R9
and the three readings.

    DEFAULT_FIT_CONFIG  the setup of a directory without fit_config.json
    load_fit_config     the setup stored in a results directory
    apply_fit_config    load it and make its medium, network and sample changes the active ones; returns the setup
    check_fit_medium    raise if a stored fit was computed in another medium than the active one
    check_fit_setup     check_fit_medium, and the same for the network and the changes to the samples

Source: Y. Alcaraz Galván
"""

import json
from pathlib import Path

from kinetics.data.observables import DEFAULT_OBSERVABLES_PATH, list_spectrum_sets
from kinetics.data.medium import DEFAULT_MEDIUM, SolvationMedium, get_active_medium, set_active_medium
from kinetics.data.network import (
    DEFAULT_NETWORK_VARIANT, NETWORK_VARIANTS, get_active_network_variant, set_active_network_variant,
)
from kinetics.fitting.candidates import set_energy_box_sigma, set_species_shifts
from kinetics.fitting.residuals import (
    DEFAULT_SCENARIOS, SCENARIOS, get_sample_overrides, set_sample_overrides,
)

FIT_CONFIG_FILE = 'fit_config.json'
DEFAULT_FIT_CONFIG = {'solvation': DEFAULT_MEDIUM.to_dict(), 'network': DEFAULT_NETWORK_VARIANT, 'sample_overrides': {},
                      'scenarios': list(DEFAULT_SCENARIOS),
                      'energy_box_sigma': None, 'structures': None,
                      'effort': {}, 'steps': [], 'step_variants': {}, 'thermo': None,
                      'species_shifts_eV': {}, 'tubes': None, 'observables': None}
RUN_OBSERVABLES = Path('data_dir') / 'data' / 'lab_observables.json'   # a run's own copy of the observables


def load_fit_config(results_dir) -> dict:
    """The setup of the run whose results are in `results_dir`: fit_config.json over the default, with the medium
    as a SolvationMedium under 'medium' and 'observables' resolved to a path (resolve_observables_path)."""
    path = Path(results_dir) / FIT_CONFIG_FILE
    stored = json.loads(path.read_text(encoding='utf-8')) if path.exists() else {}
    config = {**DEFAULT_FIT_CONFIG, **stored}
    unknown = set(config['scenarios']) - set(SCENARIOS)
    if unknown:
        raise ValueError(f'{path}: unknown scenarios {sorted(unknown)} (this branch reads {list(SCENARIOS)})')
    if config['network'] not in NETWORK_VARIANTS:
        raise ValueError(f"{path}: unknown network '{config['network']}' (one of {NETWORK_VARIANTS})")
    config['medium'] = SolvationMedium.from_dict(config['solvation'])
    config['source'] = str(path) if path.exists() else 'default (no fit_config.json)'
    config['observables'] = resolve_observables_path(results_dir, config['observables'])
    return config


def resolve_observables_path(results_dir, observables: str = None) -> str:
    """The lab observables of a run: 'observables' of its setup (relative to the results directory), else the run's own
    copy (data_dir/data/lab_observables.json), else the data folder's."""
    if observables:
        return str(Path(results_dir) / observables)
    own = Path(results_dir) / RUN_OBSERVABLES
    return str(own) if own.exists() else DEFAULT_OBSERVABLES_PATH


def list_fitted_tubes(observables: dict, config: dict = None) -> tuple:
    """The lab tubes a run can fit: 'tubes' of its setup, else every tube of the observables with a block of role
    'fit' (in their order)."""
    if config and config.get('tubes'):
        return tuple(config['tubes'])
    return tuple(name for name, record in observables['samples'].items()
                 if any(role == 'fit' for _, role, _ in list_spectrum_sets(record)))


def apply_fit_config(results_dir) -> dict:
    """load_fit_config, and its medium, network and changes to the samples made the active ones. Call it before
    computing anything from a run's fits."""
    config = load_fit_config(results_dir)
    set_active_medium(config['medium'])
    set_active_network_variant(config['network'])
    set_sample_overrides(config['sample_overrides'])
    set_energy_box_sigma(config['energy_box_sigma'])
    set_species_shifts(config['species_shifts_eV'])
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
