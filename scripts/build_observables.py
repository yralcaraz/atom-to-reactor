#!/usr/bin/env python3
"""Build lab_observables.json, the one file a fit reads from the lab NMR spectra.

    python scripts/build_observables.py                 from the share tables exported by notebook 03 (CSV)
    python scripts/build_observables.py --from-raw      from the raw spectra (needs LAB_NMR_DIR)

The file is written to the data folder (ATOM_DATA_DIR, see kinetics/paths.py).
Both routes go through kinetics.data.build_lab_observables. The raw route recomputes the area shares and reads
the heating episodes from the acquisition times; the CSV route takes the shares from the results/03
folder and the heating episodes transcribed from notebooks 03 and 04 (RECORDED_HISTORIES), and says so in the file.
Only the raw route holds the 13C spectra of the tubes with phosphate (the second nucleus of those samples).

Source: Y. Alcaraz Galván; reads share tables derived from the experimental NMR spectra
"""

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from kinetics.data import (
    DEFAULT_LAB_NMR_DIR, DEFAULT_OBSERVABLES_PATH, build_lab_inventory_for_observables, build_lab_observables,
    describe_lab_observables, load_share_tables, write_lab_observables,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--from-raw', action='store_true', help='recompute the shares from the raw spectra')
    parser.add_argument('--out', default=DEFAULT_OBSERVABLES_PATH)
    args = parser.parse_args()

    if args.from_raw:
        shares_P, shares_C, heated = build_lab_inventory_for_observables(DEFAULT_LAB_NMR_DIR)
        observables = build_lab_observables(shares_P, shares_C, heated_windows=heated, strict=True,
                                            built_from='raw spectra (folder given by LAB_NMR_DIR)')
    else:
        shares_P, shares_C = load_share_tables()
        observables = build_lab_observables(
            shares_P, shares_C,
            built_from='results/03/lab_shares_31P.csv and lab_shares_13C.csv of the data folder (area shares exported by '
                       '03_feasible_region.ipynb); heating episodes transcribed from notebooks 03 and 04. '
                       'NOT rebuilt from the raw spectra')
    path = write_lab_observables(observables, args.out)
    print(f'written: {path}')
    print(f"built from: {observables['_meta']['built_from']}")
    print(describe_lab_observables(observables).to_string())
    print('\nopen audit:')
    for item in observables['_meta']['open_audit']:
        print('  -', item)


if __name__ == '__main__':
    main()
