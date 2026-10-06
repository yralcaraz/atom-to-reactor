"""Block 14A — NMR readouts derived from measured spectra: which peak integrals a fit can use.

A readout is the integral of one measured peak. Its value is Σ_species n · c_species, with n the number of
nuclei of that species in the peak. Peaks come only from measured shifts (`nmr_shifts` in
data/experimental_gogoi2024.json): species the paper reports as one signal share a peak, species reported as
separate signals have their own. No resolution threshold and no computed shift is used to define a peak.

The counts n follow from the stored geometries: every P or Si atom for the sites 'P' and 'Si'; the carbons
bonded to Si, or the protons on those carbons, for 'Si-CH3'. A species that carries the nucleus at that site
but has no measured shift is not read: where its peak lies is unknown. The DFT shifts are tabulated next to
the measured ones only for comparison.

Source: Y. Alcaraz Galván
"""

from functools import lru_cache

import numpy as np
import pandas as pd

from kinetics.data.experimental import get_measured_shifts
from kinetics.data.network import NETWORK, NETWORK_SPECIES
from kinetics.data.snapshot import get_nmr_shieldings, get_structure
from kinetics.spectroscopy.symmetry import DEFAULT_NUCLEI, build_bond_graph

READOUT_NUCLEI = ('P', 'Si', 'C', 'H')


def find_site_atoms(species: str, nucleus: str, site: str) -> list:
    """Indices of the atoms of `nucleus` at `site` ('P', 'Si' or 'Si-CH3') in the stored geometry."""
    elements, coords = get_structure(species)
    if site in ('P', 'Si'):
        return [i for i, e in enumerate(elements) if e == nucleus == site]
    if site != 'Si-CH3':
        raise ValueError(f"Unknown site '{site}' (use 'P', 'Si' or 'Si-CH3')")
    bonded = build_bond_graph(elements, coords)
    methyl_C = [i for i, e in enumerate(elements)
                if e == 'C' and any(elements[j] == 'Si' for j in np.flatnonzero(bonded[i]))]
    if nucleus == 'C':
        return methyl_C
    if nucleus == 'H':
        return [i for i, e in enumerate(elements) if e == 'H' and any(bonded[i, j] for j in methyl_C)]
    raise ValueError(f"Site 'Si-CH3' holds C or H, not {nucleus}")


def calculate_site_dft_shift(species: str, nucleus: str, site: str) -> float:
    """Computed shift [ppm] of the site, averaged over its atoms and referenced as in the NMR catalog; NaN without data."""
    reference = DEFAULT_NUCLEI[nucleus][0]
    sigma = {}
    for sp in (species, reference):
        records = {r['atom_index']: r['isotropic_ppm'] for r in get_nmr_shieldings(sp) if r.get('element') == nucleus}
        atoms = find_site_atoms(sp, nucleus, site)
        values = [records[i] for i in atoms if i in records]
        sigma[sp] = float(np.mean(values)) if values else np.nan
    return sigma[reference] - sigma[species]


def _site_of(nucleus: str) -> str:
    sites = set(get_measured_shifts(nucleus)['site'])
    if len(sites) != 1:
        raise ValueError(f'{nucleus}: one site per nucleus expected, found {sorted(sites)}')
    return sites.pop()


@lru_cache(maxsize=4)
def _cached_readouts(nuclei: tuple) -> dict:
    readouts = {}
    for nucleus in nuclei:
        shifts = get_measured_shifts(nucleus)
        if shifts.empty:
            continue
        site = _site_of(nucleus)
        peaks = {}
        for _, row in shifts[shifts['species'].isin(NETWORK_SPECIES)].iterrows():
            peaks.setdefault(row['peak'], {})[row['species']] = len(find_site_atoms(row['species'], nucleus, site))
        readouts[nucleus] = peaks
    return readouts


def build_nmr_readouts(nuclei=READOUT_NUCLEI) -> dict:
    """{nucleus: {peak: {species: nuclei per molecule in that peak}}} from the measured shifts (cached; do not modify)."""
    return _cached_readouts(tuple(nuclei))


def find_unread_species(nucleus: str, species=None) -> list:
    """Species that carry `nucleus` at the readout site but have no measured shift, so no peak is read for them."""
    site = _site_of(nucleus)
    read = {sp for members in build_nmr_readouts((nucleus,)).get(nucleus, {}).values() for sp in members}
    return [sp for sp in (species or NETWORK_SPECIES) if sp not in read and find_site_atoms(sp, nucleus, site)]


def find_reachable_species(c0_M: dict, network: dict = None, solvent=('EC',)) -> list:
    """Species that can be present: those with c0 > 0, the solvent, and every species a reaction can form from them.

    Reactions are reversible, so a reaction also runs from its products to its reactants.
    """
    present = {sp for sp, c in c0_M.items() if c > 0} | set(solvent)
    changed = True
    while changed:
        changed = False
        for rxn in (network or NETWORK).values():
            for a, b in (('reactants', 'products'), ('products', 'reactants')):
                if set(rxn[a]) <= present and not set(rxn[b]) <= present:
                    present |= set(rxn[b])
                    changed = True
    return [sp for sp in NETWORK_SPECIES if sp in present]


def tabulate_readout_evidence(nuclei=READOUT_NUCLEI) -> pd.DataFrame:
    """One row per measured shift: the peak it defines, nuclei per molecule, and the DFT shift for comparison.

    Network species without a measured shift are listed as 'not read'.
    """
    rows = []
    for nucleus in nuclei:
        shifts = get_measured_shifts(nucleus)
        if shifts.empty:
            continue
        site = _site_of(nucleus)
        for _, r in shifts.iterrows():
            known = r['species'] in NETWORK_SPECIES
            dft = calculate_site_dft_shift(r['species'], nucleus, site) if known else np.nan
            measured = (f"{r['low_ppm']:g}" if r['low_ppm'] == r['high_ppm']
                        else f"{r['low_ppm']:g} to {r['high_ppm']:g}")
            rows.append({'nucleus': nucleus, 'site': site, 'species': r['species'], 'measured δ (ppm)': measured,
                         'source (Gogoi 2024)': r['source'], 'peak read': r['peak'] if known else 'not read',
                         'nuclei per molecule': len(find_site_atoms(r['species'], nucleus, site)) if known else np.nan,
                         'DFT δ (ppm)': dft, 'DFT − measured (ppm)': dft - r['mid_ppm'], 'note': r['note']})
        for sp in find_unread_species(nucleus):
            rows.append({'nucleus': nucleus, 'site': site, 'species': sp, 'measured δ (ppm)': '—',
                         'source (Gogoi 2024)': '—', 'peak read': 'not read: no measured shift',
                         'nuclei per molecule': len(find_site_atoms(sp, nucleus, site)),
                         'DFT δ (ppm)': calculate_site_dft_shift(sp, nucleus, site), 'DFT − measured (ppm)': np.nan,
                         'note': ''})
    return pd.DataFrame(rows)
