"""The memo's recorder, loaded by path before any ``orpheus`` code runs.

``orpheus.numerics._traced_memo_boot`` imports only the standard library, but importing it by its dotted name
runs ``orpheus/__init__`` and ``orpheus/numerics/__init__`` first, and code that ran before the recorder started
is in nobody's cone. So a capture process loads the file by path, registers it under its dotted name (the memo
imports it by that name, and a second copy would hold an empty list of active recordings; the boot module's own
``_boot`` does the same alias), and starts its ``Recording``; every ``orpheus`` import after that is recorded.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

_NAME = "orpheus.numerics._traced_memo_boot"


def boot() -> ModuleType:
    """The memo's recorder module (``Recording``, ``spawning``, ``FORK_MEASURED``, ``Witness``, ...)."""
    if _NAME in sys.modules:
        return sys.modules[_NAME]
    package = importlib.util.find_spec("orpheus")
    if package is None or not package.submodule_search_locations:
        raise ImportError("tools.cone_store: the orpheus package is not importable, so the memo's recorder cannot be found")
    path = Path(next(iter(package.submodule_search_locations))) / "numerics" / "_traced_memo_boot.py"
    spec = importlib.util.spec_from_file_location(_NAME, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"tools.cone_store: cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[_NAME] = module
    spec.loader.exec_module(module)
    return module
