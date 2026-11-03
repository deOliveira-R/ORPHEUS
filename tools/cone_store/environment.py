"""The environment key: what a stored cone is only valid under."""
from __future__ import annotations

import importlib.metadata
import json
import os
from collections.abc import Sequence
from pathlib import Path

from orpheus.numerics import traced_memo
from orpheus.numerics._traced_memo_boot import digest, file_digest

#: What the interpreter, the BLAS and pytest read from the process environment and a cone can follow, by name: the
#: memo's own list of the variables it fixes for a generation (``PYTHONHASHSEED`` and the thread counts: one list,
#: the memo's), and the ones that change what pytest runs or how numbers are summed.
VARIABLES = (*traced_memo._FIXED_ENVIRONMENT, "PYTEST_ADDOPTS", "PYTEST_PLUGINS", "PYTEST_DISABLE_PLUGIN_AUTOLOAD",
             "PYTHONWARNINGS", "MKL_CBWR")
#: The files pytest reads its configuration from, at the repository root.
INI_FILES = ("pyproject.toml", "pytest.ini", "setup.cfg", "tox.ini")
#: The repository this package belongs to.
_ROOT = Path(__file__).resolve().parents[2]
#: The modules of this package that only READ a store (``select``, ``place``): an edit to them changes no cone.
_READERS = ("selection.py", "placement.py")
#: The capture's instrument: the memo's own source (the recorder, which records nothing of its own apparatus, and
#: the pin reader; :data:`~orpheus.numerics.traced_memo._MEMO_SOURCES`, one list) and every module of this package
#: but its readers (a new module counts until it is named a reader). A cone holds no pin on them (the memo pin is
#: dropped, the apparatus never recorded), so a store is valid only under the instrument that captured it: the
#: memo's own rule for its entries, read as a key.
INSTRUMENT = (*traced_memo._MEMO_SOURCES,
              *sorted(p for p in Path(__file__).resolve().parent.glob("*.py") if p.name not in _READERS))


def distributions() -> list[list[str]]:
    """Every installed distribution as [name, version], sorted."""
    return [list(row) for row in sorted({(str(d.metadata["Name"]), d.version) for d in importlib.metadata.distributions()})]


def pytest_plugins() -> list[list[str]]:
    """Every ``pytest11`` entry point installed, as [name, value, distribution], sorted: pytest loads these before
    any recorder starts, and a distribution can gain one under an unchanged version (qa's second pass, D4)."""
    return sorted([ep.name, ep.value, str(ep.dist.metadata["Name"]) if ep.dist is not None else ""]
                  for ep in importlib.metadata.entry_points(group="pytest11"))


def environment_key(repo: Path | str | None = None, conftests: Sequence[str] = ()) -> str:
    """The digest of the interpreter (version and bytecode tag), the platform (operating system, machine, ``-O``
    level, BLAS), the installed distributions with their versions and their ``pytest11`` entry points, the
    environment variables a cone can follow (:data:`VARIABLES`), the bytes of the capture's instrument
    (:data:`INSTRUMENT`), the set of conftests a full run loads
    (``conftests``, from :mod:`tools.cone_store.conftests`, the same at capture and at select) and, given the
    repository, the bytes of its pytest configuration files. A cone captured
    under another key is not trusted (ruling 34): ``select`` lists every test of it."""
    parts = {
        "interpreter": traced_memo._python_identity(),
        "platform": traced_memo.platform_tag(),
        "distributions": distributions(),
        "pytest11": pytest_plugins(),
        "conftests": sorted(conftests),
        "instrument": {str(path.relative_to(_ROOT)): file_digest(str(path)) for path in INSTRUMENT},
        "variables": {name: os.environ.get(name) for name in VARIABLES},
        "configuration": {} if repo is None else {name: file_digest(str(Path(repo) / name)) for name in INI_FILES},
    }
    return digest(json.dumps(parts, sort_keys=True).encode())
