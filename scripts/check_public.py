#!/usr/bin/env python3
"""Check that nothing private is about to be committed (or is tracked).

    python scripts/check_public.py              the files tracked by git, plus new files not ignored
    python scripts/check_public.py --staged     the staged files only (what the pre-commit hook runs)

It looks for five things, in every text file, including the stored outputs of the notebooks:
    1. a path to the private data folders (cloud-drive folders, home directories, the university tenant)
    2. a secret (a Linear API key)
    3. the name of a raw lab folder or file (the lab data are private, see CLASSIFICATION.md)
    4. a private file name: the Tank snapshot, the lab observables, the per-spectrum share tables, raw spectra, `.env`
    5. a value of the Tank dataset itself (an energy, a solvation energy, a shielding, a coordinate, an identifier),
       in the stored unit or a converted one, in any data, code or text file. Notebooks may show and plot values:
       there only a copied geometry or a dataset identifier is refused. This needs the snapshot, so it runs only
       where ATOM_PRIVATE_DIR is set.

Exit status 1 and a list of findings when something is found. The rules are in CLASSIFICATION.md.

Classification: PUBLIC (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván
"""

import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
SELF = Path(__file__).resolve().relative_to(REPO).as_posix()

CONTENT_RULES = {
    'private path': re.compile(r'(OneDrive|CloudStorage|Uppsalauniversitet|/Users/[A-Za-z]|/home/[a-z]+/|C:\\Users\\)'),
    'secret': re.compile(r'(lin_api_[A-Za-z0-9]{8,}|LINEAR_API_KEY\s*=\s*[\'"]?[A-Za-z0-9_]{12,})'),
    'raw lab name': re.compile(r'(WW-NG|WW-S\d|NG2305\d\d|NG2306\d\d|FRESH_Hus|Hus ?7\b|VAL - NMR raw)'),
}
FORBIDDEN_NAMES = [
    re.compile(r'(^|/)tank_api_snapshot\.json$'), re.compile(r'(^|/)lab_observables\.json$'),
    re.compile(r'(^|/)lab_sample_folders\.json$'), re.compile(r'(^|/)lab_shares_[^/]*\.csv$'),
    re.compile(r'(^|/)block5_legacy_baseline\.json$'), re.compile(r'(^|/)sync_linear\.py$'),
    re.compile(r'\.jdf$'), re.compile(r'(^|/)\.env$'), re.compile(r'(^|/)docs-internal/'),
]
BINARY_SUFFIXES = {'.png', '.jpg', '.jpeg', '.gif', '.pdf', '.pyc', '.zip', '.npz', '.npy', '.pkl'}


EH_TO_EV, EH_TO_KJ, EV_TO_KJ = 27.211386245988, 2625.4996394799, 96.48533212
NUMBER = re.compile(r'(?<![\w.])-?\d+\.\d+(?![\w])')
# fewest significant digits for a number to count as a copy of a dataset value (fewer would match by chance)
MIN_DIGITS = {'DFT energy': 6, 'DFT entropy': 5, 'DFT property': 5, 'solvation': 5, 'NMR shielding': 5, 'geometry': 6}


def load_database_index() -> tuple:
    """({(kind, decimals, 'value'): label}, {identifier: label}) of the private Tank snapshot; empty when it is absent."""
    import importlib.util                     # load kinetics/paths.py alone: the package itself needs numpy
    spec = importlib.util.spec_from_file_location('atom_paths', REPO / 'kinetics' / 'paths.py')
    paths = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(paths)
    snapshot_path = Path(paths.private_path('data', 'tank_api_snapshot.json'))
    if not snapshot_path.is_file():
        return {}, {}
    snapshot = json.loads(snapshot_path.read_text(encoding='utf-8'))
    values, identifiers = [], {}

    def add(kind, value, label):
        if isinstance(value, (int, float)) and not isinstance(value, bool) and value == value and value != 0:
            values.append((kind, float(value), label))

    for record in snapshot['datasets']:
        name = record['pipeline_id']
        for key in ('energy_scf_eh', 'zpe_eh', 'enthalpy_eh', 'gibbs_eh'):
            if record.get(key) is not None:
                for factor in (1.0, EH_TO_EV, EH_TO_KJ):
                    add('DFT energy', record[key] * factor, f'{name}.{key}')
        if record.get('enthalpy_eh') is not None and record.get('gibbs_eh') is not None:
            add('DFT entropy', (record['enthalpy_eh'] - record['gibbs_eh']) * EH_TO_KJ * 1000 / 298.15, f'{name}.S')
        for key in ('dipole_debye', 'polarizability_iso_ang3', 'homo_ev', 'lumo_ev', 'gap_ev', 'optical_gap_ev'):
            add('DFT property', record.get(key), f'{name}.{key}')
        for key in ('id', 'dataset_uuid', 'molecule_id', 'inchikey'):
            if record.get(key):
                identifiers[str(record[key])] = f'{name}.{key}'
    for record in snapshot['solvation']:
        name = record['raw_metadata'].get('name', record.get('molecule_name'))
        for key in ('delta_e_solv_kjmol', 'uncertainty_kjmol'):
            if record.get(key) is not None:
                add('solvation', record[key], f'{name}.{key}')
                add('solvation', record[key] / EV_TO_KJ, f'{name}.{key}')
        add('solvation', record.get('density_g_cm3'), f'{name}.density')
        for key, value in record['raw_metadata'].items():
            if isinstance(value, (int, float)) and not isinstance(value, bool) and (key.startswith(('E_', 'dE_')) or 'density' in key):
                add('solvation', value, f'{name}.{key}')
                add('solvation', value * EV_TO_KJ, f'{name}.{key}')
            elif key == 'smiles':
                identifiers[value] = f'{name}.smiles'
        for key in ('id', 'molecule_id', 'solvent_molecule_id'):
            if record.get(key):
                identifiers[str(record[key])] = f'{name}.{key}'
    names = {record['dataset_uuid']: record['pipeline_id'] for record in snapshot['datasets']}
    for uuid, record in snapshot['nmr'].items():
        for atom in record['nmr_shieldings']:
            add('NMR shielding', atom.get('isotropic_ppm'), f"{names.get(uuid, uuid)}.{atom['element']}{atom['atom_index']}")
            add('NMR shielding', atom.get('anisotropy_ppm'), f"{names.get(uuid, uuid)}.{atom['element']}{atom['atom_index']}")
    for uuid, record in snapshot['structure'].items():
        for position in record['coordinates_angstrom']:
            for coordinate in position:
                add('geometry', coordinate, f'{names.get(uuid, uuid)} geometry')
    index = {}
    for kind, value, label in values:
        for decimals in range(1, 13):
            index.setdefault((kind, decimals, f'{abs(value):.{decimals}f}'), label)
    return index, {text: label for text, label in identifiers.items() if len(text) >= 12}


def scan_database(path: str, text: str, index: dict, identifiers: dict) -> list:
    """Findings for numbers of `text` that equal a dataset value at their own precision."""
    in_notebook = path.endswith('.ipynb')
    if in_notebook:                                   # sources and text outputs; images are skipped
        notebook = json.loads(text)
        parts = []
        for cell in notebook.get('cells', []):
            parts.append(''.join(cell.get('source', [])))
            for output in cell.get('outputs', []):
                parts.append(''.join(output.get('text', [])) + ''.join(output.get('data', {}).get('text/plain', [])))
        text = '\n'.join(parts)
    findings = [(path, 0, 'Tank dataset identifier', label) for ident, label in identifiers.items() if ident in text]
    for match in NUMBER.finditer(text):
        token = match.group(0)
        decimals = min(len(token.split('.')[1]), 12)
        digits = len(token.lstrip('-').replace('.', '').lstrip('0'))
        value = f'{abs(float(token)):.{decimals}f}'
        for kind, need in MIN_DIGITS.items():
            if in_notebook and kind != 'geometry':    # results and interpretation of the data are public
                continue
            if digits >= need and (kind, decimals, value) in index:
                findings.append((path, 0, f'Tank dataset value ({kind})', f'{token} = {index[(kind, decimals, value)]}'))
    return findings


def git_lines(*args) -> list:
    done = subprocess.run(['git', *args], cwd=REPO, capture_output=True, text=True, check=True)
    return [line for line in done.stdout.split('\n') if line]


def candidate_files(staged: bool) -> list:
    if staged:
        return git_lines('diff', '--cached', '--name-only', '--diff-filter=ACMR')
    return sorted(set(git_lines('ls-files')) | set(git_lines('ls-files', '--others', '--exclude-standard')))


def scan(path: str, index: dict, identifiers: dict) -> list:
    findings = []
    for rule in FORBIDDEN_NAMES:
        if rule.search(path):
            findings.append((path, 0, 'private file', path))
    target = REPO / path
    if not target.is_file() or target.suffix.lower() in BINARY_SUFFIXES or path == SELF:
        return findings
    try:
        text = target.read_text(encoding='utf-8')
    except (UnicodeDecodeError, OSError):
        return findings
    if index:
        findings.extend(scan_database(path, text, index, identifiers))
    for number, line in enumerate(text.split('\n'), start=1):
        for kind, rule in CONTENT_RULES.items():
            match = rule.search(line)
            if match:
                start = max(0, match.start() - 25)
                findings.append((path, number, kind, line[start:match.end() + 25].strip()))
    return findings


def main() -> int:
    staged = '--staged' in sys.argv
    index, identifiers = load_database_index()
    findings = []
    for path in candidate_files(staged):
        findings.extend(scan(path, index, identifiers))
    database = 'Tank dataset values checked' if index else 'Tank dataset values NOT checked (ATOM_PRIVATE_DIR not set)'
    if not findings:
        print('check_public: nothing private found' + (' in the staged files' if staged else '') + f' ({database})')
        return 0
    print(f'check_public: {len(findings)} finding(s). Private data must stay out of the repository (CLASSIFICATION.md).')
    for path, number, kind, excerpt in findings[:60]:
        where = f'{path}:{number}' if number else path
        print(f'  [{kind}] {where}: {excerpt[:110]}')
    if len(findings) > 60:
        print(f'  ... and {len(findings) - 60} more')
    return 1


if __name__ == '__main__':
    sys.exit(main())
