"""Blocks 2–4 — Thermochemistry: species free energies in EC and reaction thermodynamics.

    gas.py             gas-phase G, H, S (stored DFT values or quasi-RRHO)
    standard_state.py  1 bar ideal gas → 1 M solution shift
    solution.py        G_sol = G_gas + ΔE_solv + ΔG°→*
    uncertainty.py     σ of the solvation contribution to a reaction
    reaction.py        ΔG_rxn, K_eq, network table, Wegscheider cycles
"""

from kinetics.thermo.gas import THERMO_MODES, calculate_gas_thermo
from kinetics.thermo.standard_state import calculate_standard_state_shift
from kinetics.thermo.solution import calculate_solution_gibbs
from kinetics.thermo.uncertainty import calculate_solvation_sigma, calculate_solvation_sigma_naive
from kinetics.thermo.reaction import calculate_cycle_residuals, calculate_network_thermo, calculate_reaction_thermo

__all__ = [
    "THERMO_MODES", "calculate_gas_thermo", "calculate_standard_state_shift", "calculate_solution_gibbs",
    "calculate_solvation_sigma", "calculate_solvation_sigma_naive",
    "calculate_reaction_thermo", "calculate_network_thermo", "calculate_cycle_residuals",
]
