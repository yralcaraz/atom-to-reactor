"""Where the private data live. The one place that knows about them.

Unpublished data (the Tank snapshot, the curated lab observables, the per-spectrum share tables) and the raw
lab spectra are not stored in the repository. They sit in a folder that each user points to with two
environment variables, set in the shell or in a local `.env` file at the repository root (git-ignored; copy
`.env.example`):

    ATOM_PRIVATE_DIR   folder with the private data, laid out as in the repository:
                         data/       tank_api_snapshot.json, lab_observables.json, lab_sample_folders.json
                         results/    01/ and 03/ tables that reproduce private inputs
                         tests/      regression baseline and the tests that name private lab files
    LAB_NMR_DIR        folder with the raw lab spectra (JEOL .jdf files)

No path to either folder appears anywhere in the code, the notebooks or the docs. Without the variables
the code still imports: a loader that needs a private file raises FileNotFoundError, and the tests that need
one skip.

    private_dir        the private data folder, or None when ATOM_PRIVATE_DIR is not set
    require_private_dir  the same, but raises when it is not set (before writing a private file)
    private_path       a path inside it (a placeholder path when it is not set, so the error message is clear)
    lab_nmr_dir        the raw spectra folder, or None when LAB_NMR_DIR is not set

Classification: PUBLIC (see CLASSIFICATION.md)
Source: Y. Alcaraz Galván
"""

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
PRIVATE_DIR_VAR = 'ATOM_PRIVATE_DIR'
LAB_NMR_DIR_VAR = 'LAB_NMR_DIR'
_PLACEHOLDER = f'<{PRIVATE_DIR_VAR} is not set: see .env.example>'


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


def private_dir() -> 'Path | None':
    """The private data folder, or None when ATOM_PRIVATE_DIR is not set."""
    value = os.environ.get(PRIVATE_DIR_VAR)
    return Path(value).expanduser() if value else None


def private_path(*parts: str) -> str:
    """Path of a file inside the private data folder. When the variable is not set the path is a placeholder
    that does not exist, so a loader fails with a message that says what to set."""
    base = private_dir()
    return os.path.join(str(base) if base else _PLACEHOLDER, *parts)


def require_private_dir() -> Path:
    """The private data folder; raises when ATOM_PRIVATE_DIR is not set (use before writing a private file)."""
    folder = private_dir()
    if folder is None:
        raise FileNotFoundError(f'{PRIVATE_DIR_VAR} is not set: it holds the private data (see .env.example)')
    return folder


def lab_nmr_dir() -> 'Path | None':
    """The raw lab spectra folder, or None when LAB_NMR_DIR is not set."""
    value = os.environ.get(LAB_NMR_DIR_VAR)
    return Path(value).expanduser() if value else None


def require_lab_nmr_dir() -> Path:
    """The raw lab spectra folder; raises when LAB_NMR_DIR is not set or the folder is missing."""
    folder = lab_nmr_dir()
    if folder is None or not folder.is_dir():
        raise FileNotFoundError(f'Raw lab spectra not found: set {LAB_NMR_DIR_VAR} to the folder that holds them '
                                f'(see .env.example)')
    return folder
