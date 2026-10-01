"""Block 14B — Experiment design: what an NMR time series of a candidate experiment could determine.

A design is a composition at t = 0 (the moment the additive is added), a temperature program and, per
nucleus, the times at which a quantitative spectrum is recorded. The readouts are integrals of the measured
peaks (readouts.py) expressed as concentrations of that nucleus [M]; each has Gaussian noise of a fixed σ per
nucleus, a placeholder until the NMR facility gives it.

For parameters θ (family barriers g, optionally activation entropies ΔS‡) the Fisher information is

    F = Σ_designs Jᵀ J + diag(1/σ_prior²),      J_ij = (∂y_i/∂θ_j) / σ_i

and sqrt(diag F⁻¹) is the smallest standard error an unbiased fit could reach (Cramér-Rao). The prior only
regularises the inversion: a parameter the data do not inform keeps σ ≈ σ_prior. F is local, so every
result holds for one assumed truth; designs are compared over a range of plausible truths.
"""

from dataclasses import dataclass, replace

import numpy as np
import pandas as pd

from kinetics.constants import EV_TO_KJ_MOL, ZERO_CELSIUS_K
from kinetics.data.network import NETWORK, NETWORK_SPECIES
from kinetics.fitting.feasibility import get_family_parameter
from kinetics.fitting.readouts import build_nmr_readouts, find_reachable_species, find_unread_species
from kinetics.microkinetics.models import get_model
from kinetics.reactor.batch import simulate_batch_reactor
from kinetics.reactor.engine import MassActionSystem
from kinetics.reactor.observables import find_crossing_time

# Peak groups per nucleus, derived from the measured shifts (kinetics/fitting/readouts.py):
# {nucleus: {peak: {species: nuclei per molecule in that peak}}}. Species without a measured shift are not read.
NMR_READOUTS = build_nmr_readouts()

# Acquisition per nucleus: time per quantitative spectrum and standard error of one peak integral. PLACEHOLDERS:
# no source in the repository gives them (Gogoi 2024 report no acquisition parameters, integrals or noise).
# The Fisher σ of a well-determined parameter is proportional to sigma_mM, so every σ in the design analysis
# scales with these numbers; 'needs' says what replaces each one.
DEFAULT_ACQUISITION = {
    'P': {'interval_h': 5.0 / 60.0, 'sigma_mM': 1.0, 'status': 'placeholder',
          'needs': 'T1 of the silyl phosphates in EC and the S/N of one quantitative ³¹P spectrum (NMR facility)'},
    'Si': {'interval_h': 0.5, 'sigma_mM': 5.0, 'status': 'placeholder',
           'needs': 'T1 and the S/N of one inverse-gated ²⁹Si spectrum; Gogoi 2024 did not detect impurity-level '
                    'HMDSO in ²⁹Si that ¹³C showed'},
    'C': {'interval_h': 1.0, 'sigma_mM': 5.0, 'status': 'placeholder',
          'needs': 'T1 and the S/N of one quantitative ¹³C spectrum'},
}

DEFAULT_PRIOR_SIGMA = {'g': 0.30, 'dS': 100.0}      # eV, J/(mol·K): regularisation only
_PARAMETER_KEYS = {'g': ('g_eV', 0.005), 'dS': ('dS_act_J_mol_K', 5.0)}   # key, finite-difference step


@dataclass(frozen=True)
class ExperimentDesign:
    name: str
    label: str
    c0_M: dict                  # concentrations at t = 0 [M]; EC is held constant as the solvent
    segments: tuple             # ((T_C, duration_h), ...) run back to back from t = 0
    sampling_h: dict            # {nucleus: tuple of acquisition times [h]}
    purpose: str = ''

    @property
    def duration_h(self) -> float:
        return float(sum(d for _, d in self.segments))


def build_composition(EC_M: float, **c_M) -> dict:
    """Concentrations [M] of every network species: EC plus the given ones, the rest zero."""
    c0 = dict.fromkeys(NETWORK_SPECIES, 0.0)
    c0['EC'] = EC_M
    for sp, c in c_M.items():
        if sp not in c0:
            raise KeyError(f"'{sp}' is not a network species")
        c0[sp] = c
    return c0


def build_in_situ_design(name: str, label: str, c0_M: dict, *, segments, nuclei=('P', 'Si'),
                         dead_time_h: float = 1.0 / 6.0, acquisition: dict = None, purpose: str = '') -> ExperimentDesign:
    """Spectra recorded in the magnet at the reaction temperature, every acquisition interval from the dead time
    on, through a program of segments ((T_C, duration_h), ...)."""
    acq = acquisition or DEFAULT_ACQUISITION
    segments = tuple((float(T), float(d)) for T, d in segments)
    total_h = sum(d for _, d in segments)
    sampling = {n: tuple(np.arange(dead_time_h, total_h + 1e-9, acq[n]['interval_h'])) for n in nuclei}
    return ExperimentDesign(name, label, c0_M, segments, sampling, purpose)


def build_isothermal_design(name: str, label: str, c0_M: dict, *, T_C: float, duration_h: float,
                            nuclei=('P', 'Si'), dead_time_h: float = 1.0 / 6.0, acquisition: dict = None,
                            purpose: str = '') -> ExperimentDesign:
    """One temperature, spectra in situ (see build_in_situ_design)."""
    return build_in_situ_design(name, label, c0_M, segments=((T_C, duration_h),), nuclei=nuclei,
                                dead_time_h=dead_time_h, acquisition=acquisition, purpose=purpose)


def build_step_design(name: str, label: str, c0_M: dict, *, hold_temps_C, hold_h: float = 8.0,
                      acquire_T_C: float = 20.0, acquire_h: float = 1.0, nuclei=('P', 'Si'),
                      purpose: str = '') -> ExperimentDesign:
    """Holds at increasing temperature, each followed by one acquisition of every nucleus at acquire_T_C."""
    segments, times, t = [], [], 0.0
    for T_C in hold_temps_C:
        segments += [(float(T_C), float(hold_h)), (float(acquire_T_C), float(acquire_h))]
        t += hold_h + acquire_h
        times.append(t)
    return ExperimentDesign(name, label, c0_M, tuple(segments), {n: tuple(times) for n in nuclei}, purpose)


def calculate_readouts(C_M: np.ndarray, idx: dict, nucleus: str) -> dict:
    """{peak label: concentration of the nucleus in that peak [M]} over the last axis of C_M."""
    return {label: sum(n * C_M[idx[sp]] for sp, n in members.items())
            for label, members in NMR_READOUTS[nucleus].items()}


def simulate_design(design: ExperimentDesign, model, *, network: dict = None, species_db: dict = None,
                    points_per_segment: int = 80, rtol: float = 1e-8, atol: float = 1e-12) -> dict:
    """Integrate the design and return the trajectory and the readouts at every acquisition time."""
    spec = get_model(model)
    net = network or NETWORK
    species = list(design.c0_M)
    system = MassActionSystem(net, species, {'EC': design.c0_M['EC']} if 'EC' in species else {})
    all_samples = np.unique(np.concatenate([np.asarray(t, dtype=float) for t in design.sampling_h.values()]))

    c = np.array([design.c0_M[sp] for sp in species], dtype=float)
    t0, t_parts, c_parts, T_parts, rate_cache = 0.0, [], [], [], {}
    for T_C, duration_h in design.segments:
        T_K = T_C + ZERO_CELSIUS_K
        if T_K not in rate_cache:
            rates = spec.calculate_rates(T_K, network=net, species_db=species_db)
            rate_cache[T_K] = (rates['k_f'].values, rates['k_r'].values)
        inside = all_samples[(all_samples > t0) & (all_samples <= t0 + duration_h)] - t0
        t_eval_h = np.unique(np.concatenate([np.linspace(0.0, duration_h, points_per_segment), inside]))
        sol = system.integrate(c, *rate_cache[T_K], t_eval_h * 3600.0, rtol=rtol, atol=atol)
        t_parts.append(t0 + sol.t / 3600.0)
        c_parts.append(sol.y)
        T_parts.append(np.full(sol.t.shape, T_C))
        c = sol.y[:, -1]
        t0 += duration_h

    t_h = np.concatenate(t_parts)
    C_M = np.concatenate(c_parts, axis=1)
    readouts = {}
    for nucleus, times in design.sampling_h.items():
        times = np.asarray(times, dtype=float)
        readouts[nucleus] = {label: np.interp(times, t_h, trace)
                             for label, trace in calculate_readouts(C_M, system.idx, nucleus).items()}
    return {'design': design.name, 'model': spec.name, 't_h': t_h, 'T_C': np.concatenate(T_parts), 'C_M': C_M,
            'species': species, 'idx': system.idx, 'sampling_h': design.sampling_h, 'readouts': readouts}


def split_family_by_reaction(model, family: str, network: dict = None) -> tuple:
    """(model, network, new families) in which every member of `family` has its own class '<family>_<rxn>'.

    Each new class starts with the family's parameters, so the kinetics are unchanged; their barriers can
    then be varied, and fitted, one reaction at a time.
    """
    spec = get_model(model)
    net = {rxn_id: {**rxn} for rxn_id, rxn in (network or NETWORK).items()}
    members = [rxn_id for rxn_id, rxn in net.items() if rxn['class'] == family]
    if not members:
        raise KeyError(f"No reaction of family '{family}' in the network")
    params = {fam: dict(p) for fam, p in spec.family_params.items()}
    base = params.get(family, params.get('default', {}))
    for rxn_id in members:
        net[rxn_id]['class'] = f'{family}_{rxn_id}'
        params[f'{family}_{rxn_id}'] = dict(base)
    variant = replace(spec, name=f'{spec.name} [{family} per reaction]', family_params=params)
    return variant, net, [f'{family}_{rxn_id}' for rxn_id in members]


def add_readout_noise(result: dict, acquisition: dict = None, seed: int = 7) -> dict:
    """Synthetic measurement: {nucleus: {peak label: readout + Gaussian noise of that nucleus' σ}} [M]."""
    acq = acquisition or DEFAULT_ACQUISITION
    rng = np.random.default_rng(seed)
    return {nucleus: {label: values + rng.normal(0.0, acq[nucleus]['sigma_mM'] / 1000.0, values.shape)
                      for label, values in peaks.items()}
            for nucleus, peaks in result['readouts'].items()}


def describe_designs(designs, acquisition: dict = None) -> pd.DataFrame:
    """One row per design: purpose, composition, temperature program, spectra per nucleus and duration."""
    acq = acquisition or DEFAULT_ACQUISITION
    rows = {}
    for d in designs:
        solutes = ', '.join(f'{sp} {1000.0 * c:.0f} mM' for sp, c in d.c0_M.items() if c > 0 and sp != 'EC')
        program = ' → '.join(f'{T:g} °C × {h:g} h' for T, h in d.segments)
        if len(d.segments) > 6:
            holds = [T for T, _ in d.segments[::2]]
            program = f'holds {holds[0]:g}…{holds[-1]:g} °C × {d.segments[0][1]:g} h, spectra at {d.segments[1][0]:g} °C'
        rows[d.name] = {
            'experiment': d.label,
            'purpose': d.purpose,
            'solutes at t = 0': solutes,
            'temperature program': program,
            'spectra': ', '.join(f'{n}: {len(t)}' for n, t in d.sampling_h.items()),
            'first spectrum (min)': round(60.0 * min(min(t) for t in d.sampling_h.values())),
            'duration (h)': d.duration_h,
            'can form, not read': '; '.join(
                f"{n}: {', '.join(unread)}" for n in d.sampling_h
                if (unread := find_unread_species(n, find_reachable_species(d.c0_M)))) or '—',
        }
    return pd.DataFrame.from_dict(rows, orient='index')


def _parse_parameter(parameter: str):
    """'hydrolysis' → ('hydrolysis', 'g'); 'hydrolysis:dS' → ('hydrolysis', 'dS')."""
    family, _, kind = parameter.partition(':')
    kind = kind or 'g'
    if kind not in _PARAMETER_KEYS:
        raise ValueError(f"Unknown parameter kind '{kind}' in '{parameter}' (use g or dS)")
    return family, kind


def _readout_vector(result: dict, acquisition: dict):
    values, sigmas, labels = [], [], []
    for nucleus, peaks in result['readouts'].items():
        sigma_M = acquisition[nucleus]['sigma_mM'] / 1000.0
        for label, series in peaks.items():
            values.append(series)
            sigmas.append(np.full(series.shape, sigma_M))
            labels += [(nucleus, label, t) for t in result['sampling_h'][nucleus]]
    return np.concatenate(values), np.concatenate(sigmas), labels


def calculate_design_information(design: ExperimentDesign, model, parameters, *, acquisition: dict = None,
                                 network: dict = None, species_db: dict = None) -> dict:
    """Noise-weighted sensitivities of every readout to every parameter, and the Fisher information JᵀJ."""
    spec = get_model(model)
    acq = acquisition or DEFAULT_ACQUISITION
    y0, sigma, labels = _readout_vector(simulate_design(design, spec, network=network, species_db=species_db), acq)
    J = np.zeros((len(y0), len(parameters)))
    for j, parameter in enumerate(parameters):
        family, kind = _parse_parameter(parameter)
        key, step = _PARAMETER_KEYS[kind]
        value = get_family_parameter(spec, family, key)
        y = [_readout_vector(simulate_design(design, spec.with_family_params(family, **{key: value + s}),
                                             network=network, species_db=species_db), acq)[0]
             for s in (step, -step)]
        J[:, j] = (y[0] - y[1]) / (2.0 * step)
    Jw = J / sigma[:, None]
    return {'design': design.name, 'parameters': list(parameters), 'y0': y0, 'sigma': sigma, 'labels': labels,
            'J': J, 'Jw': Jw, 'F': Jw.T @ Jw}


def calculate_parameter_precision(information, parameters, prior_sigma: dict = None) -> tuple:
    """(σ per parameter, correlation matrix) from one or several information results (summed).

    prior_sigma: {'g': eV, 'dS': J/(mol·K)} regularisation; σ close to it means the data do not inform.
    """
    infos = information if isinstance(information, (list, tuple)) else [information]
    prior = {**DEFAULT_PRIOR_SIGMA, **(prior_sigma or {})}
    F = sum(info['F'] for info in infos)
    F = F + np.diag([1.0 / prior[_parse_parameter(p)[1]] ** 2 for p in parameters])
    cov = np.linalg.inv(F)
    sigma = np.sqrt(np.diag(cov))
    corr = cov / np.outer(sigma, sigma)
    return pd.Series(sigma, index=parameters), pd.DataFrame(corr, index=parameters, columns=parameters)


def compare_designs(designs, model, parameters, *, combinations: dict = None, acquisition: dict = None,
                    prior_sigma: dict = None, network: dict = None, species_db: dict = None) -> pd.DataFrame:
    """σ of every parameter for each design alone and for named combinations {name: [design names]}."""
    infos = {d.name: calculate_design_information(d, model, parameters, acquisition=acquisition,
                                                  network=network, species_db=species_db) for d in designs}
    rows = {name: calculate_parameter_precision(info, parameters, prior_sigma)[0] for name, info in infos.items()}
    for name, members in (combinations or {}).items():
        rows[name] = calculate_parameter_precision([infos[m] for m in members], parameters, prior_sigma)[0]
    return pd.DataFrame(rows).T


def scan_design_precision(designs, model, parameters, *, vary, values, acquisition: dict = None,
                          prior_sigma: dict = None, network: dict = None, species_db: dict = None) -> pd.DataFrame:
    """σ of every parameter for each design, over assumed true values of `vary`.

    vary: one parameter ('family' or 'family:dS') or a list of them, all set to each value.
    Returns one row per (design, true value, parameter).
    """
    spec = get_model(model)
    varied = [vary] if isinstance(vary, str) else list(vary)
    rows = []
    for value in values:
        truth = spec
        for parameter in varied:
            family, kind = _parse_parameter(parameter)
            truth = truth.with_family_params(family, **{_PARAMETER_KEYS[kind][0]: float(value)})
        for d in designs:
            info = calculate_design_information(d, truth, parameters, acquisition=acquisition, network=network,
                                                species_db=species_db)
            sigma = calculate_parameter_precision(info, parameters, prior_sigma)[0]
            rows += [{'design': d.name, 'true value': float(value), 'parameter': p, 'sigma': s}
                     for p, s in sigma.items()]
    return pd.DataFrame(rows)


def calculate_half_life_map(model, c0_M: dict, *, family: str, g_values_eV, T_C_values, species: str = 'TMSPA',
                            t_max_h: float = 1.0e4, network: dict = None, species_db: dict = None) -> pd.DataFrame:
    """Half-life [h] of `species` at constant T for each (g of `family`, T); NaN if longer than t_max_h."""
    spec = get_model(model)
    table = {}
    for T_C in T_C_values:
        column = []
        for g in g_values_eV:
            sim = simulate_batch_reactor(c0_M, T_K=T_C + ZERO_CELSIUS_K, t_end_s=t_max_h * 3600.0,
                                         model=spec.with_family_params(family, g_eV=float(g)), network=network,
                                         species_db=species_db, n_points=400, t_start_s=1e-3)
            trace = sim['C_M'][sim['idx'][species]]
            column.append(find_crossing_time(sim['t_h'], trace, 0.5 * trace[0]))
        table[float(T_C)] = column
    return pd.DataFrame(table, index=pd.Index([float(g) for g in g_values_eV], name='g (eV)'))


def calculate_eyring_precision(temperature_sets_C: dict, sigma_g_eV: float) -> pd.DataFrame:
    """Precision of ΔH‡ and ΔS‡ from barriers measured at several temperatures, each with error sigma_g_eV.

    Linear regression of g(T) = ΔH‡ − T·ΔS‡: σ(ΔS‡) = σ_g / sqrt(Σ(T − T̄)²), σ(ΔH‡) = σ_g·sqrt(1/n + T̄²/Σ(T − T̄)²).
    """
    sigma_J_mol = sigma_g_eV * EV_TO_KJ_MOL * 1000.0
    rows = {}
    for name, temps_C in temperature_sets_C.items():
        T = np.asarray(temps_C, dtype=float) + ZERO_CELSIUS_K
        sxx = np.sum((T - T.mean()) ** 2)
        rows[name] = {
            'temperatures (°C)': ', '.join(f'{t:g}' for t in temps_C),
            'n': len(T),
            'span (K)': float(np.ptp(T)),
            'σ ΔS‡ (J/mol/K)': sigma_J_mol / np.sqrt(sxx) if sxx > 0 else np.inf,
            'σ ΔH‡ (kJ/mol)': sigma_J_mol * np.sqrt(1.0 / len(T) + T.mean() ** 2 / sxx) / 1000.0 if sxx > 0 else np.inf,
        }
    return pd.DataFrame.from_dict(rows, orient='index')


def calculate_barrier_error_budget(dG_barrier_eV: float, T_C: float, *, dT_K: float = 0.5,
                                   rel_conc_error: float = 0.05, sigma_dG_rxn_eV: dict = None,
                                   alpha: float = 0.5) -> pd.DataFrame:
    """Systematic errors on a barrier measured at T_C, beside the noise-limited σ of the Fisher analysis.

    The experiment measures ΔG‡ of each step. A sample temperature off by dT_K and a co-reactant concentration
    off by rel_conc_error bias ΔG‡ directly. The family g follows from ΔG‡ through the barrier model, so it
    also inherits the error of ΔG_rxn with weight α (∂g/∂ΔG_rxn ≈ −α at fixed ΔG‡).
    sigma_dG_rxn_eV: {label: σ of ΔG_rxn [eV]}.
    """
    T_K = T_C + ZERO_CELSIUS_K
    kT = 8.617333262e-5 * T_K
    # Exact: the rate at T + dT read as if it were at T
    d_temp = abs(dG_barrier_eV * T_K / (T_K + dT_K) - kT * np.log1p(dT_K / T_K) - dG_barrier_eV)
    d_conc = kT * abs(np.log1p(rel_conc_error))
    rows = {
        f'sample temperature ±{dT_K:g} K': (d_temp, d_temp),
        f'initial co-reactant concentration ±{100 * rel_conc_error:g} %': (d_conc, d_conc),
    }
    for label, sigma in (sigma_dG_rxn_eV or {}).items():
        rows[f'ΔG_rxn ±{sigma:.2g} eV ({label})'] = (0.0, alpha * sigma)
    out = pd.DataFrame.from_dict(rows, orient='index', columns=['on measured ΔG‡ (meV)', 'on family g (meV)']) * 1000.0
    out['factor on k'] = np.exp(out['on measured ΔG‡ (meV)'] / 1000.0 / kT)
    return out
