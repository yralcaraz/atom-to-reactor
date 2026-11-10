"""Block 14T — Thermodynamic fit: the solvation energies of the species from the spectra taken at equilibrium.

A tube whose reaction has finished no longer changes: its shares depend on the reaction energies alone, not on the
barriers and not on when the first spectrum was taken (notebook 06 §1). So the energies are fitted first, with an
algebraic model (kinetics.thermo.equilibrium, no reactor, no ages), and the kinetic fit (M1) takes them afterwards.

    G_i  = G_gas,i + E_solv,i
    E_solv,i ∈ ΔE_i ± ξ·SD_i  (SD_i: standard error of the MD run)

    ΔG_j = Σ_i ν_ij G_i,   K_j = exp(−ΔG_j / RT)      every reaction in the span of R1–R4 (R5–R7, N10, N12–N15 too)

    χ²   = Σ z²,   z = (predicted − measured) / σ     the shares of the spectra at equilibrium, error model of the
                                                       kinetic fit (residuals._lab_block, residuals._whiten)
                                                       — or, setup 'error_model': 'nmr' (fit_thermo_rev01 §2.1),
                                                       the concentrations in mM with σ_tot = noise ⊕ baseline, every
                                                       point independent (_concentration_block, _standardise)

The free parameters are the shifts δ_i of E_solv,i of the species of R1–R4, inside a box of ±ξ standard errors and
with no other constraint (notebook 06 §5); an amount at mixing that the recipe does not fix is a parameter as in the
kinetic fit. Only differences of species energies enter ΔG, so the data see at most four combinations of the seven
δ (one per reaction R1–R4): the rest is flat inside the box. The result to read is ΔG_j and K_j; a single δ_i is
determined only where a profile says so.

The file reads from the whole to its parts: the problem first, then what it gives, then where its data come from.

    0  Setup              DEFAULT_THERMO_SETUP, load_thermo_setup: the 'thermo' part of fit_config.json
    1  The problem        ThermoProblem: parameters → equilibrium → predicted shares → z → χ²
    2  What it gives      tabulate_reaction_energies, calculate_energy_ranges: ΔG and K, and their range
    3  The samples        build_equilibrium_samples: the tubes at equilibrium, with their errors (what §1 reads)
    4  Which spectra      tabulate_lab_spectra → find_equilibrium_points (the rule)
    5  Replicates, steps  compare_replicates → build_cumulative_steps: the order in which tubes enter the fit

A run goes the other way, 4 → 5 → 3 → 1 → 2 (scripts/run_thermo_fit.py). The chain of §1, from the parameters to χ²:

    chi2(x) = Σ residuals(x)²
      residuals(x)         one z per measured share: residuals._whiten(block, predicted)
        _predicted         shares of each block from the concentrations (residuals._shares_from_state)
          equilibrium      concentrations at equilibrium with μ_i = G_sol,i + δ_i (thermo.equilibrium.calculate_equilibrium)
            theta(x)       the optimiser's vector as {parameter name: value}, the held ones included

Source: Y. Alcaraz Galván; notebook 06 rev01 (the problem), docs/reference-fitting_thermo.md
"""

# ------------------------------------------------------------------------------
# IMPORTS
# -----------------------------------------------------------------------------

from types import SimpleNamespace

import numpy as np
import pandas as pd

from kinetics.constants import T_STD_K
from kinetics.data.lab_spectra import FLOOR_OF, classify_tube_files, list_tube_folders, tabulate_tube_shares
from kinetics.data.network import format_equation, get_active_network
from kinetics.data.observables import DEFAULT_ERROR_FLOOR, RT_MAX_C
from kinetics.fitting.candidates import RECIPE_PREFIX, FitParameter
from kinetics.fitting.estimation import DELTA_CHI2_95
from kinetics.fitting.residuals import (
    PENALTY, SHARE_SPECIES, _lab_block, _rows_of_block, _shares_from_state, _whiten, get_sample_overrides,
)
from kinetics.constants import EV_TO_KJ_MOL, R_SI, ZERO_CELSIUS_K
from kinetics.thermo.equilibrium import calculate_equilibrium, list_equilibrium_reactions
from kinetics.thermo.reaction import calculate_reaction_thermo
from kinetics.thermo.solution import calculate_solution_gibbs
from kinetics.thermo.uncertainty import calculate_solvation_sigma_components


# -----------------------------------------------------------------------------
# CONSTANTS AND DEFAULTS
#-----------------------------------------------------------------------------

SOLVATION_PREFIX = 'dE_solv '            # 'dE_solv <species>': shift of its solvation energy from the computed value [eV]
QUANTITATIVE_31P_D1_S = 5.0              # ³¹P: relaxation delay of a quantitative spectrum (the observables' selection)

DEFAULT_THERMO_SETUP = {
    # the box of every free solvation energy: ±ξ standard errors of its MD run; each ξ is fitted (as the readings of M1)
    'xi': [1.0, 2.0, 3.0],
    'free_species': ['TMSPA', 'BMSPA', 'MMSPA', 'H3PO4', 'H2O', 'TMSOH', 'HMDSO'],
    'reactions': ['R1', 'R2', 'R3', 'R4'],
    'nuclei': {'fit': ['31P', '13C'], 'check': ['1H', '29Si']},     
    # 'check': scored after every step, never fitted
    'rule': {'min_spectra': 2, 'min_span_h': 24.0, 'max_abs_z': 3.0, 'floor': False},
    'overrides': {},        # {tube: {file: true / false}}: decided by the setup, not the rule
    'effort': [256, 8, 40], # Sobol screen points, starts polished, iterations per start
    'profile_nfev': 20,     # number of function evaluations per profile point (the optimiser is re-run at each point)
    'order': [],            # the order in which tubes enter the fit (the rest is alphabetical)
    'steps': 'cumulative',  # 'cumulative': each step adds the next tube, 'independent': each step fits one tube alone
    # 'kinetic': shares, error model of the kinetic fit (floor, replicate blocks, √((K−1)/K)); 'nmr': fit_thermo_rev01
    # §2.1, concentrations (mM), σ_tot = √(σ_noise² + σ_baseline²), every point independent
    'error_model': 'kinetic',
    # the box of E_solv,i: 'solute' = gas ⊕ solution standard errors (the shared pure-solvent term left out);
    # 'total' = gas ⊕ solution ⊕ solvent, the snapshot's uncertainty of the solvent correction (fit_thermo_rev01 §2.1)
    'sigma_solv': 'solute',
}

def load_thermo_setup(config: dict) -> dict:
    """
    The 'thermo' part of a run's setup (kinetics.fitting.config.load_fit_config), over DEFAULT_THERMO_SETUP.
    """
    setup = {**DEFAULT_THERMO_SETUP, **(config.get('thermo') or {})}
    setup['rule'] = {**DEFAULT_THERMO_SETUP['rule'], **setup['rule']}
    setup['nuclei'] = {**DEFAULT_THERMO_SETUP['nuclei'], **setup['nuclei']}
    return setup


# ------------------------------------------------------------------------------
# 1. Parameters → equilibrium → predicted shares → z → xi^2
# ------------------------------------------------------------------------------
class ThermoProblem:
    """
    The fit problem, as a box that holds the tubes and the free energies.

        ThermoProblem(samples, free_species, xi)     create it once: it keeps the tubes and the energy limits
            .residuals(x)                            the optimiser gives energies x, the box returns the errors z
            .chi2(x)                                 the sum of z²

    Inside .residuals(x):
        x  →  equilibrium of each tube  →  predicted shares  →  z = (predicted − measured) / σ
    """


    def __init__(self, samples: list, *, free_species, xi: float, reactions=('R1', 'R2', 'R3', 'R4'), fixed: dict = None,
                 network: dict = None, solver: dict = None, sigma_solv: str = 'solute'):
        # saves tubes, species, reacions, network and xi
        self.samples = samples
        self.free_species, self.xi, self.reactions = tuple(free_species), float(xi), tuple(reactions)
        self.network = network or get_active_network()

        # Calculates sd (standard error of the solvation energy of each species): 'solute' = gas ⊕ solution box,
        # 'total' = gas ⊕ solution ⊕ pure-solvent box (the snapshot's uncertainty of the solvent correction)
        self.solver = solver
        self.sigma_solv = sigma_solv
        self.structure = SimpleNamespace(name='thermo', label=f'solvation energies within ±{self.xi:g} standard errors')
        components = calculate_solvation_sigma_components()
        if sigma_solv == 'solute':
            sd = {sp: s for sp, (s, _) in components.items()}
        elif sigma_solv == 'total':
            sd = {sp: float(np.hypot(s, ref)) for sp, (s, ref) in components.items()}
        else:
            raise ValueError(f"sigma_solv '{sigma_solv}': 'solute' or 'total'")
        
        # Creates a list of parameters
        self.all_parameters = [FitParameter(f'{SOLVATION_PREFIX}{sp}', 'solvation', sp, -self.xi * sd[sp], self.xi * sd[sp], 'eV')
                               for sp in self.free_species]
        for sample in samples:
            for species, (low, high) in sample.get('unknown_M', {}).items():
                self.all_parameters.append(FitParameter(f"{RECIPE_PREFIX}{species} {sample['name']}", 'recipe', sample['name'], low, high, 'M'))
        
        # Separates the fixed parameters from the free ones
        self.fixed = dict(fixed or {})
        known = {p.name for p in self.all_parameters}
        self.fixed = {k: v for k, v in self.fixed.items() if k in known}      # an earlier step may hold more tubes
        self.parameters = [p for p in self.all_parameters if p.name not in self.fixed]
        self.names = [p.name for p in self.parameters]
        
        # Saves the lower and upper bounds of the free parameters
        self.lower = np.array([p.lower for p in self.parameters], dtype=float)
        self.upper = np.array([p.upper for p in self.parameters], dtype=float)
        
        # Calculates the initial energy of each species in each tube (_mu0) and the row number (idx)
        species = {sp for r in list_equilibrium_reactions(self.reactions, self.network)
                   for side in ('reactants', 'products') for sp in self.network[r][side]}
        self._mu0 = {s['name']: {sp: calculate_solution_gibbs(sp, s['T_K'])['G_sol_eV'] for sp in species} for s in samples}
        self._idx = {sp: i for i, sp in enumerate(sorted({sp for n in SHARE_SPECIES.values() for w in n.values() for sp in w}))}

    # From parameters to chi^2
    def chi2(self, x) -> float:
        r = self.residuals(x)   # list of all z, one for each mesuredment
        return float(r @ r)     # z_1^2 + z_2^2 + ... + z_n^2


    def residuals(self, x, roles=('fit',)) -> np.ndarray:
        theta, out = self.theta(x), []
        for sample in self.samples:
            blocks = [b for b in sample['blocks'] if b['role'] in roles]
            if not blocks:
                continue
            predicted = self._predicted(sample, theta)
            if predicted is None:
                out.append(np.full(sum(b['measured'].size for b in blocks), PENALTY))
                continue
            out += [_standardise(b, p) for b, p in zip(sample['blocks'], predicted) if b['role'] in roles]
        return np.concatenate(out) if out else np.zeros(0)

    def _predicted(self, sample: dict, theta: dict) -> list:
        """What every block of a sample predicts (spectra × windows, the same row for every spectrum): shares, or for a
        block in mM (error model 'nmr') the concentration of each window, Σ n·C_i / divisor; None on failure."""
        try:
            C = self.equilibrium(sample, theta)
        except (RuntimeError, np.linalg.LinAlgError, FloatingPointError):
            return None
        column = np.zeros((len(self._idx), 1))
        for sp, c in C.items():
            if sp in self._idx:
                column[self._idx[sp], 0] = c
        out = []
        for b in sample['blocks']:
            if b.get('units') == 'mM':
                members = SHARE_SPECIES[b['nucleus']]
                amounts = np.array([sum(n * max(C.get(sp, 0.0), 0.0) for sp, n in members[w].items()) for w in b['windows']])
                row = (1000.0 * amounts / b['divisor'])[None, :]
            else:
                row = _shares_from_state(column, self._idx, b['nucleus'], b['windows'])
            out.append(np.repeat(row, len(b['t_spec_s']), axis=0))
        return out

    def equilibrium(self, sample: dict, theta: dict) -> dict:
        """Composition at equilibrium [mol/L] of one sample under a parameter set."""
        shifts = self.species_shifts(theta)
        c0 = dict(sample['c0_M'])
        for species in sample.get('unknown_M', {}):
            c0[species] = c0.get(species, 0.0) + max(theta[f"{RECIPE_PREFIX}{species} {sample['name']}"], 0.0)
        mu = {sp: g + shifts.get(sp, 0.0) for sp, g in self._mu0[sample['name']].items()}
        return calculate_equilibrium(c0, sample['T_K'], reactions=self.reactions, network=self.network, mu_eV=mu)

    def species_shifts(self, theta: dict) -> dict:
        """{species: eV} moved from the computed solvation energies (the hand-over to the kinetic fit)."""
        return {name[len(SOLVATION_PREFIX):]: v for name, v in theta.items() if name.startswith(SOLVATION_PREFIX)}

    def theta(self, x) -> dict:
        values = {**self.fixed, **dict(zip(self.names, np.asarray(x, dtype=float)))}
        return {p.name: float(values[p.name]) for p in self.all_parameters}

    def table(self, x, roles=('fit',)) -> pd.DataFrame:
        """Standardised residuals, one row per measured quantity (residuals._rows_of_block: 'share', or 'mean' and
        'drift' of a block of repeated spectra), with 'nucleus' and 'role'. The fit rule reads the 'fit' rows."""
        theta, rows = self.theta(x), []
        for sample in self.samples:
            predicted = self._predicted(sample, theta)
            for i, block in enumerate(sample['blocks']):
                if block['role'] not in roles:
                    continue
                if predicted is None:
                    rows.append({'sample': sample['name'], 'quantity': 'equilibrium not found', 'kind': 'failure', 'z': PENALTY,
                                 'nucleus': block['nucleus'], 'role': block['role']})
                    continue
                rows += [{**row, 'nucleus': block['nucleus'], 'role': block['role'],
                          **({'kind': 'concentration (mM)'} if block.get('units') == 'mM' else {})}
                         for row in _rows_of_block(sample['name'], block, predicted[i], block['t_spec_s'] / 3600.0)]
        table = pd.DataFrame(rows, columns=['sample', 'nucleus', 'role', 'quantity', 'kind', 't_h', 'measured', 'predicted', 'sigma', 'z'])
        table['age_h'] = np.nan
        return table

    # Support functions for the optimiser
    @property
    def n(self) -> int:
        return len(self.parameters)

    def vector(self, theta: dict) -> np.ndarray:
        return np.clip([theta.get(name, 0.0) for name in self.names], self.lower, self.upper)

    def to_unit(self, x) -> np.ndarray:
        return 1.0 + (np.asarray(x, dtype=float) - self.lower) / (self.upper - self.lower)

    def from_unit(self, u) -> np.ndarray:
        return self.lower + (np.asarray(u, dtype=float) - 1.0) * (self.upper - self.lower)

    def replace(self, *, samples: list = None, solver: dict = None, fixed: dict = None):
        return ThermoProblem(self.samples if samples is None else samples, free_species=self.free_species, xi=self.xi,
                             reactions=self.reactions, fixed=self.fixed if fixed is None else fixed, network=self.network,
                             solver=self.solver if solver is None else solver, sigma_solv=self.sigma_solv)


def _standardise(block: dict, predicted: np.ndarray) -> np.ndarray:
    """z of one block: for a block in mM (error model 'nmr', fit_thermo_rev01 §2.1) (y_model − y_NMR) / σ_tot, every
    point independent; otherwise the whitened residuals of the kinetic fit (residuals._whiten)."""
    if block.get('units') == 'mM':
        return ((predicted - block['measured']) / block['sigma_ind']).ravel()
    return _whiten(block, predicted)



# ------------------------------------------------------------------------------
# 2  What it gives: ΔG and K of every reaction, and their range over a profile
# ------------------------------------------------------------------------------

def tabulate_reaction_energies(species_shifts_eV: dict, T_K: float = T_STD_K, *, network: dict = None) -> pd.DataFrame:
    # Builds a table of the computed and shifted ΔG and K of every reaction in the network, at a temperature T_K.
    net = network or get_active_network()
    RT_eV = R_SI * T_K / 1000.0 / EV_TO_KJ_MOL
    rows = {}
    for rxn_id, rxn in net.items():
        computed = calculate_reaction_thermo(rxn_id, T_K, network=net)['dG_rxn_eV']
        shift = sum(nu * species_shifts_eV.get(sp, 0.0) for sp, nu in rxn['products'].items()) \
            - sum(nu * species_shifts_eV.get(sp, 0.0) for sp, nu in rxn['reactants'].items())
        rows[rxn_id] = {'equation': format_equation(rxn), 'class': rxn['class'], 'computed ΔG (eV)': computed,
                        'shift (eV)': shift, 'ΔG (eV)': computed + shift, 'computed K': np.exp(-computed / RT_eV),
                        'K': np.exp(-(computed + shift) / RT_eV)}
    return pd.DataFrame.from_dict(rows, orient='index')


def calculate_energy_ranges(profile: pd.DataFrame, problem: ThermoProblem, *, delta: float = DELTA_CHI2_95,
                            T_K: float = T_STD_K) -> pd.DataFrame:
    """Lowest and highest ΔG of each reaction over the profile points within `delta` of the lowest χ² (the parameters
    re-minimised at each point, 'θ …' columns of calculate_profile): an envelope of the values the data allow, not a
    profile of ΔG itself. One row per reaction: 'ΔG low (eV)', 'ΔG high (eV)', 'points'."""
    accepted = profile[profile['chi2'] <= profile['chi2'].min() + delta]
    values = []
    for _, row in accepted.iterrows():
        theta = {c[2:]: row[c] for c in profile.columns if c.startswith('θ ')}
        values.append(tabulate_reaction_energies(problem.species_shifts(theta), T_K, network=problem.network)['ΔG (eV)'])
    if not values:
        return pd.DataFrame(columns=['ΔG low (eV)', 'ΔG high (eV)', 'points'])
    values = pd.concat(values, axis=1)
    return pd.DataFrame({'ΔG low (eV)': values.min(axis=1), 'ΔG high (eV)': values.max(axis=1), 'points': values.shape[1]})


# ------------------------------------------------------------------------------
# 3  The samples: the tubes at equilibrium, with their blocks and errors (what ThermoProblem reads)
# ------------------------------------------------------------------------------
def build_equilibrium_samples(shares: pd.DataFrame, points: pd.DataFrame, observables: dict, nuclei: dict, *,
                              error_floor: dict = None, error_model: str = 'kinetic') -> list:
    """The tubes with spectra at equilibrium, as samples: 'name', 'c0_M' (the recipe, with the run's 'added_M' on
    top), 'unknown_M' (the record's, or the run's where it replaces them), 'T_K' (mean temperature of the spectra
    used) and 'blocks' (one per nucleus, the error model of residuals._lab_block: three spectra or more of one block
    share their baseline and floor). nuclei: {'fit': [...], 'check': [...]}; a block of a 'check' nucleus has the
    role 'check' and enters no χ².
    error_model 'nmr' (fit_thermo_rev01 §2.1): every block in mM instead, _concentration_block (no floor, no shared errors)."""
    if error_model not in ('kinetic', 'nmr'):
        raise ValueError(f"error_model '{error_model}': 'kinetic' or 'nmr'")
    base = error_floor or observables['_meta'].get('error_floor', DEFAULT_ERROR_FLOOR)
    floor = {nucleus: base[FLOOR_OF[nucleus]] for nucleus in FLOOR_OF}
    roles = {**{n: 'check' for n in nuclei.get('check', [])}, **{n: 'fit' for n in nuclei.get('fit', [])}}
    chosen = points[points['at equilibrium'] & points['nucleus'].isin(list(roles))]
    overrides = get_sample_overrides()
    samples = []
    for tube in dict.fromkeys(chosen['tube']):
        record, change = observables['samples'][tube], overrides.get(tube, {})
        sets = {}
        for nucleus in dict.fromkeys(chosen.loc[chosen['tube'] == tube, 'nucleus']):
            md5 = chosen[(chosen['tube'] == tube) & (chosen['nucleus'] == nucleus)].sort_values('t_h')['data_md5']
            sets[nucleus] = [_spectrum_record(shares[(shares['tube'] == tube) & (shares['data_md5'] == m)]) for m in md5]
        # residuals._lab_block reads the replicate scatter of a block from a record of this shape
        first = next(iter(sets))
        view = {'samples': {tube: {'nucleus': first, 'spectra': sets[first],
                                   'more_nuclei': {n: {'spectra': s} for n, s in sets.items() if n != first}}}}
        c0 = dict(record['c0_M'])
        for species, added in change.get('added_M', {}).items():
            c0[species] = c0.get(species, 0.0) + added
        if error_model == 'nmr':
            blocks = [_concentration_block(tube, nucleus, roles[nucleus], spectra, c0) for nucleus, spectra in sets.items()]
        else:
            blocks = [_lab_block(tube, nucleus, roles[nucleus], spectra, view, floor, 1.0) for nucleus, spectra in sets.items()]
        for block, spectra in zip(blocks, sets.values()):
            block['files'] = [sp['file'] for sp in spectra]
        unknown = change['unknown_M'] if 'unknown_M' in change else record.get('unknown_M', {})
        temperatures = [sp['temperature_C'] for spectra in sets.values() for sp in spectra]
        samples.append({'name': tube, 'kind': 'equilibrium', 'c0_M': c0, 'T_K': float(np.mean(temperatures)) + ZERO_CELSIUS_K,
                        'unknown_M': {sp: (float(low), float(high)) for sp, (low, high) in (unknown or {}).items()},
                        'blocks': blocks})
    return samples


def _concentration_block(name: str, nucleus: str, role: str, spectra: list, c0: dict) -> dict:
    """The spectra of one nucleus of one tube in mM (fit_thermo_rev01 §2.1, as its tabulate_concentrations):

        y_NMR = share · total / n,   σ_tot = √(σ_noise² + σ_baseline²) · total / n

    total: the nucleus at mixing over the windows of the spectrum, Σ_windows Σ_i n_i·c0_i; n: the count of the species of a
    one-species window (HMDSO: 2 silyl groups), 1 for a window of several species (its concentration counts groups).
    Every point is independent: no floor, no error shared between spectra."""
    members = SHARE_SPECIES[nucleus]
    unknown = set(spectra[0]['shares']) - set(members)
    if unknown:
        raise KeyError(f'{name}: no species declared for the {nucleus} window(s) {sorted(unknown)}')
    windows = tuple(w for w in members if w in spectra[0]['shares'])
    total = sum(n * c0.get(sp, 0.0) for w in windows for sp, n in members[w].items())
    divisor = np.array([next(iter(members[w].values())) if len(members[w]) == 1 else 1 for w in windows], dtype=float)
    scale = 1000.0 * total / divisor
    Y = np.array([[sp['shares'][w]['share'] for w in windows] for sp in spectra])
    sigma = np.array([[np.hypot(sp['shares'][w]['sigma_noise'], sp['shares'][w]['sigma_baseline']) for w in windows]
                      for sp in spectra])
    return {'nucleus': nucleus, 'role': role, 'units': 'mM', 'windows': windows, 'divisor': divisor,
            't_spec_s': 3600.0 * np.array([sp['t_since_first_h'] for sp in spectra]),
            'measured': Y * scale, 'sigma_ind': sigma * scale, 'sigma_common': np.zeros(len(windows)),
            'sigma_nmr': sigma * scale, 'replicate': False, 'birge': 1.0, 'cholesky': None}


def _spectrum_record(rows: pd.DataFrame) -> dict:
    first = rows.iloc[0]
    return {'t_since_first_h': float(first['t_h']), 'temperature_C': float(first['T_C']), 'file': first['file'],
            'shares': {r['window']: {'share': float(r['share']), 'sigma_noise': float(r['sigma_noise']),
                                     'sigma_baseline': float(r['sigma_baseline'])} for _, r in rows.iterrows()}}


# ------------------------------------------------------------------------------
# 4  Which spectra are at equilibrium
# ------------------------------------------------------------------------------
# Which spectra are at equilibrium (find_equilibrium_points). Every usable spectrum of every tube
# (kinetics.data.lab_spectra: not noise) is a candidate if it was taken at room temperature, is quantitative (³¹P with
# the NOE off and a relaxation delay of at least 5 s; ²⁹Si not INEPT) and has at least two windows. The candidates of one
# tube and one nucleus, in time order, are at equilibrium from the first one on which every later one agrees:
#
#     a plateau  = at least `min_spectra` spectra spanning at least `min_span_h` hours, each within `max_abs_z` of the
#                  mean of the others, window by window (σ: noise ⊕ baseline of each spectrum; with floor True also
#                  the floor, as one spectrum has it in the fit)
#
# A spectrum with no later spectrum of its nucleus cannot show that nothing changes after it, so it is never at
# equilibrium by the rule. The setup can decide a file either way ('overrides'), and the table says who decided.
# The rule cannot tell an equilibrium from a reaction that stopped (a slow step, a reagent used up): equal spectra are
# the evidence it asks for, nothing more.
def tabulate_lab_spectra(observables: dict, inventory: pd.DataFrame) -> pd.DataFrame:
    """Shares of every usable spectrum of every tube of the observables (kinetics.data.lab_spectra), one row per
    (spectrum, window). inventory: build_lab_nmr_inventory of the raw folder, copies included."""
    tables = [tabulate_tube_shares(tube, observables, classify_tube_files(tube, observables, inventory))
              for tube in list_tube_folders() if tube in observables['samples']]
    return pd.concat(tables, ignore_index=True)

def find_equilibrium_points(shares: pd.DataFrame, *, rule: dict = None, overrides: dict = None,
                            error_floor: dict = None) -> pd.DataFrame:
    """One row per spectrum of the shares table (tabulate_lab_spectra): 'tube', 'nucleus', 'file', 'data_md5', 't_h',
    'T_C', 'candidate' (room temperature, quantitative, two windows or more), 'by rule', 'at equilibrium' (the rule, or
    the setup where it decides the file), 'decided by', 'plateau from h', 'plateau spectra', 'max |z|' (of the
    plateau, or of the spectrum against the later ones when they differ) and 'reason'."""
    rule = {**DEFAULT_THERMO_SETUP['rule'], **(rule or {})}
    overrides = overrides or {}
    error_floor = error_floor or DEFAULT_ERROR_FLOOR
    spectra = (shares.groupby(['tube', 'nucleus', 'data_md5'], sort=False)
               .agg(file=('file', 'first'), t_h=('t_h', 'first'), T_C=('T_C', 'first'), noe=('noe', 'first'),
                    relaxation_delay_s=('relaxation_delay_s', 'first'), experiment=('experiment', 'first'),
                    n_windows=('window', 'size')).reset_index())
    rows = []
    for (tube, nucleus), group in spectra.groupby(['tube', 'nucleus'], sort=False):
        group = group.sort_values('t_h')
        reasons = {}
        for _, r in group.iterrows():
            why = _not_quantitative(r)
            reasons[r['data_md5']] = (f"{r['T_C']:.0f} °C" if r['T_C'] >= RT_MAX_C else f'not quantitative ({why})' if why
                                      else 'one window: no share to compare' if r['n_windows'] < 2 else '')
        candidates = group[[not reasons[m] for m in group['data_md5']]]
        blocks = [shares[(shares['tube'] == tube) & (shares['data_md5'] == m)] for m in candidates['data_md5']]
        Y = np.array([b['share'].to_numpy() for b in blocks]) if blocks else np.zeros((0, 0))
        S = np.array([_spectrum_sigma(b, error_floor, rule['floor']) for b in blocks]) if blocks else np.zeros((0, 0))
        t = candidates['t_h'].to_numpy()
        start, z_plateau = None, None
        for i in range(len(t)):
            if len(t) - i < rule['min_spectra'] or t[-1] - t[i] < rule['min_span_h']:
                break
            z = _plateau_z(Y[i:], S[i:])
            if z.max() <= rule['max_abs_z']:
                start, z_plateau = i, z
                break
        for k, (_, r) in enumerate(candidates.iterrows()):
            if start is not None and k >= start:
                reasons[r['data_md5']] = (f"plateau of {len(t) - start} spectra from {t[start]:.1f} h to {t[-1]:.1f} h, "
                                          f"max |z| {z_plateau.max():.2f}")
            elif k == len(t) - 1:
                reasons[r['data_md5']] = 'no later spectrum of this nucleus'
            elif t[-1] - t[k] < rule['min_span_h']:
                reasons[r['data_md5']] = f"later spectra within {rule['min_span_h']:g} h only"
            else:
                z = _plateau_z(Y[k:], S[k:])
                reasons[r['data_md5']] = f'differs from the later spectra (max |z| {z.max():.2f})'
        for _, r in group.iterrows():
            k = list(candidates['data_md5']).index(r['data_md5']) if r['data_md5'] in set(candidates['data_md5']) else None
            by_rule = start is not None and k is not None and k >= start
            decided = overrides.get(tube, {}).get(r['file'])
            max_z = (float(z_plateau.max()) if by_rule else float(_plateau_z(Y[k:], S[k:]).max())
                     if k is not None and len(t) - k >= 2 else np.nan)
            rows.append({'tube': tube, 'nucleus': nucleus, 'file': r['file'], 'data_md5': r['data_md5'], 't_h': r['t_h'],
                         'T_C': r['T_C'], 'candidate': k is not None, 'by rule': by_rule,
                         'at equilibrium': bool(decided) if decided is not None else by_rule,
                         'decided by': 'setup' if decided is not None else 'rule',
                         'plateau from h': t[start] if by_rule else np.nan,
                         'plateau spectra': len(t) - start if by_rule else 0, 'max |z|': max_z,
                         'reason': reasons[r['data_md5']]})
    return pd.DataFrame(rows)

def _not_quantitative(row) -> str:
    """Why a spectrum is not quantitative ('' when it is)."""
    if row['nucleus'] == '31P' and (str(row['noe']).upper() == 'TRUE' or row['relaxation_delay_s'] < QUANTITATIVE_31P_D1_S):
        return ('NOE on, ' if str(row['noe']).upper() == 'TRUE' else '') + f"d1 {row['relaxation_delay_s']:g} s"
    if row['nucleus'] == '29Si' and 'inept' in str(row['experiment']).lower():
        return 'INEPT'
    return ''


def _spectrum_sigma(rows: pd.DataFrame, error_floor: dict, floor: bool) -> np.ndarray:
    sigma = np.hypot(rows['sigma_noise'].to_numpy(), rows['sigma_baseline'].to_numpy())
    if floor:
        f = error_floor[FLOOR_OF[rows['nucleus'].iloc[0]]]
        sigma = np.hypot(sigma, np.maximum(f['absolute'], f['relative'] * np.maximum(rows['share'].to_numpy(), 0.0)))
    return sigma


def _plateau_z(Y: np.ndarray, S: np.ndarray) -> np.ndarray:
    """|z| of each spectrum (row) against the noise-weighted mean of the others, window by window."""
    z = np.zeros_like(Y)
    for k in range(len(Y)):
        others = np.arange(len(Y)) != k
        w = 1.0 / S[others] ** 2
        mean, sigma_mean = (w * Y[others]).sum(axis=0) / w.sum(axis=0), 1.0 / np.sqrt(w.sum(axis=0))
        z[k] = np.abs(Y[k] - mean) / np.hypot(S[k], sigma_mean)
    return z

# ------------------------------------------------------------------------------
# 5  Replicates and steps: the order in which the tubes enter the fit
# ------------------------------------------------------------------------------
def compare_replicates(shares: pd.DataFrame, points: pd.DataFrame, observables: dict, *, max_abs_z: float = 3.0,
                       error_floor: dict = None, floor: bool = False) -> pd.DataFrame:
    """Tubes mixed from the same recipe (same composition at mixing, the run's 'added_M' included) set against each other
    at equilibrium: per nucleus they share, the noise-weighted mean of the spectra at equilibrium of each tube, window by
    window (σ: noise ⊕ baseline, with the floor if `floor`). Equal recipes must end at the same shares; when they do not,
    one model cannot sit between them, and they are two solutions (build_cumulative_steps). One row per pair and nucleus:
    'tube', 'other', 'nucleus', 'max |z|', 'compatible'; a pair with no nucleus in common has one row with nucleus None."""
    error_floor = error_floor or observables['_meta'].get('error_floor', DEFAULT_ERROR_FLOOR)
    overrides = get_sample_overrides()

    def recipe(tube):
        c0 = dict(observables['samples'][tube]['c0_M'])
        for species, added in overrides.get(tube, {}).get('added_M', {}).items():
            c0[species] = c0.get(species, 0.0) + added
        return tuple(sorted((sp, round(c, 9)) for sp, c in c0.items()))

    def means(tube, nucleus):
        md5 = points[points['at equilibrium'] & (points['tube'] == tube) & (points['nucleus'] == nucleus)]['data_md5']
        blocks = [shares[(shares['tube'] == tube) & (shares['data_md5'] == m)] for m in md5]
        if not blocks:
            return None
        Y = np.array([b['share'].to_numpy() for b in blocks])
        w = 1.0 / np.array([_spectrum_sigma(b, error_floor, floor) for b in blocks]) ** 2
        return (w * Y).sum(axis=0) / w.sum(axis=0), 1.0 / np.sqrt(w.sum(axis=0))

    tubes = list(dict.fromkeys(points.loc[points['at equilibrium'], 'tube']))
    rows = []
    for i, tube in enumerate(tubes):
        for other in tubes[i + 1:]:
            if recipe(tube) != recipe(other):
                continue
            common = [n for n in dict.fromkeys(points['nucleus']) if means(tube, n) is not None and means(other, n) is not None]
            for nucleus in common:
                (a, sa), (b, sb) = means(tube, nucleus), means(other, nucleus)
                z = float(np.max(np.abs(a - b) / np.hypot(sa, sb)))
                rows.append({'tube': tube, 'other': other, 'nucleus': nucleus, 'max |z|': z, 'compatible': z <= max_abs_z})
            if not common:
                rows.append({'tube': tube, 'other': other, 'nucleus': None, 'max |z|': np.nan, 'compatible': True})
    return pd.DataFrame(rows, columns=['tube', 'other', 'nucleus', 'max |z|', 'compatible'])


def build_cumulative_steps(order, replicates: pd.DataFrame, tubes) -> list:
    """Steps that add one tube at a time and fit every tube added so far, each from the step before (nothing fitted is
    dropped). A tube that disagrees with a replicate already in a chain (compare_replicates) splits it: one chain goes on
    with the replicate, the other with the tube in its place, so equal recipes that end apart give two solutions.
    order: the tubes in the order they are added; tubes not in `tubes` (none of their spectra at equilibrium) are skipped.
    Step ids 'T<k>' with a letter per chain once it splits ('T2a', 'T2b')."""
    clash = {frozenset((r['tube'], r['other'])) for _, r in replicates.iterrows() if not r['compatible']}
    chains = [{'label': '', 'tubes': [], 'last': None, 'history': []}]
    steps = []
    for k, tube in enumerate([t for t in order if t in tubes], start=1):
        grown = []
        for chain in chains:
            against = [t for t in chain['tubes'] if frozenset((t, tube)) in clash]
            if against:          # the chain as it is, and a copy with the tube in place of the replicates it disagrees with
                # the copy starts from the last step of the chain that did not fit those replicates
                before = [sid for sid, fitted in chain['history'] if not set(fitted) & set(against)]
                grown += [{**chain, 'label': chain['label'] + 'a'},
                          {**chain, 'label': chain['label'] + 'b', 'tubes': [t for t in chain['tubes'] if t not in against],
                           'last': before[-1] if before else None, 'swap': tube}]
            else:
                grown.append({**chain, 'add': tube})
        chains = []
        for chain in grown:
            new = chain.pop('add', None) or chain.pop('swap', None)
            if new is None:                   # the chain that keeps the replicate: nothing new at this step
                chains.append(chain)
                continue
            fitted = chain['tubes'] + [new]
            step = {'id': f"T{k}{chain['label']}", 'fit': fitted}
            if chain['last']:
                step['from'] = chain['last']
            steps.append(step)
            chains.append({'label': chain['label'], 'tubes': fitted, 'last': step['id'],
                           'history': chain['history'] + [(step['id'], fitted)]})
    return steps
