"""The conftests a full run loads: every ``conftest.py`` under the session's test paths.

The store promises full-run verdicts (the orchestrator's ruling on qa's second pass, from ruling 41: the full suite
is the merge gate). A full run imports every conftest under its test paths while it collects, before any test runs,
so each one's module code and session hooks (``pytest_configure``, ``pytest_collection_modifyitems``) act on every
test. One function finds them; the grandparent (and a fresh child) loads exactly these before any test runs, and
the environment key carries their set, so a conftest that appears anywhere selects everything.
"""
from __future__ import annotations

import os
from collections.abc import Iterable
from pathlib import Path
from typing import Any

#: Directories pytest does not recurse into by default (``norecursedirs``'s defaults that can hold a conftest).
_SKIPPED = ("__pycache__", "node_modules", "build", "dist")


def under(anchors: Iterable[Path]) -> list[Path]:
    """Every ``conftest.py`` under ``anchors`` (files or directories), sorted, skipping hidden directories and
    :data:`_SKIPPED`."""
    found: set[Path] = set()
    for anchor in anchors:
        root = anchor if anchor.is_dir() else anchor.parent
        for directory, names, files in os.walk(root):
            names[:] = [n for n in names if not n.startswith(".") and n not in _SKIPPED and not n.endswith(".egg-info")]
            if "conftest.py" in files:
                found.add(Path(directory, "conftest.py").resolve())
    return sorted(found)


def of_config(config: Any) -> list[Path]:
    """The conftests a full run of ``config``'s project loads: those under its ``testpaths`` (the rootdir when it
    names none), whatever the command line selected."""
    root = Path(config.rootpath)
    paths = [root / p for p in config.getini("testpaths")] or [root]
    return under(p for p in paths if p.exists())


def relative(repo: Path | str, paths: Iterable[Path]) -> list[str]:
    real = Path(os.path.realpath(repo))
    return [os.path.relpath(p, real) for p in paths]
