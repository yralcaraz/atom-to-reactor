"""Where the data files are. The one place that resolves their location.

The data files (the Tank snapshot, the lab observables, the per-spectrum share tables) and the NMR spectra are
read from folders that each user sets with two environment variables, in the shell or in a local `.env` file
at the repository root (copy `.env.example`):

    ATOM_DATA_DIR   data folder, laid out as in the repository:
                      data/       tank_api_snapshot.json, lab_observables.json, lab_sample_folders.json
                      results/    tables written by notebooks 01 and 03
                      tests/      regression baseline
    LAB_NMR_DIR     folder with the NMR spectra (JEOL .jdf files)

Without the variables the code still imports: a loader that needs a data file raises FileNotFoundError, and
the tests skip.

    data_dir          the data folder, or None when ATOM_DATA_DIR is not set
    require_data_dir  the same, but raises when it is not set (before writing a file there)
    data_path         a path inside it (a placeholder path when it is not set, so the error message is clear)
    lab_nmr_dir       the NMR spectra folder, or None when LAB_NMR_DIR is not set

Source: Y. Alcaraz Galván
"""

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR_VAR = 'ATOM_DATA_DIR'
LAB_NMR_DIR_VAR = 'LAB_NMR_DIR'
_PLACEHOLDER = f'<{DATA_DIR_VAR} is not set: see .env.example>'


def _read_dotenv(path: Path) -> None:
    """Load KEY=VALUE lines of `.env` into the environment, without overriding what is already set."""
    if not path.is_file():
        return
    for line in path.read_text(encoding='utf-8').splitlines():
        line = line.strip()
        if not line or line.startswith('#') or '=' not in line:
            continue
        key, _, value = line.partition('=')
        os.environ.setdefault(key.strip(), value.strip().strip('"').strip("'"))


_read_dotenv(REPO_ROOT / '.env')


def data_dir() -> 'Path | None':
    """The data folder, or None when ATOM_DATA_DIR is not set."""
    value = os.environ.get(DATA_DIR_VAR)
    return Path(value).expanduser() if value else None


def data_path(*parts: str) -> str:
    """Path of a file inside the data folder. When the variable is not set the path is a placeholder
    that does not exist, so a loader fails with a message that says what to set."""
    base = data_dir()
    return os.path.join(str(base) if base else _PLACEHOLDER, *parts)


def require_data_dir() -> Path:
    """The data folder; raises when ATOM_DATA_DIR is not set (use before writing a file there)."""
    folder = data_dir()
    if folder is None:
        raise FileNotFoundError(f'{DATA_DIR_VAR} is not set: it is the folder with the data files (see .env.example)')
    return folder


def lab_nmr_dir() -> 'Path | None':
    """The NMR spectra folder, or None when LAB_NMR_DIR is not set."""
    value = os.environ.get(LAB_NMR_DIR_VAR)
    return Path(value).expanduser() if value else None


def require_lab_nmr_dir() -> Path:
    """The NMR spectra folder; raises when LAB_NMR_DIR is not set or the folder is missing."""
    folder = lab_nmr_dir()
    if folder is None or not folder.is_dir():
        raise FileNotFoundError(f'NMR spectra not found: set {LAB_NMR_DIR_VAR} to the folder that holds them '
                                f'(see .env.example)')
    return folder
