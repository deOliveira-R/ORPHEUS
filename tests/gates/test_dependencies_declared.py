r"""Every third-party module ``orpheus`` imports is a declared dependency (#405 P1 step 6, S6.20).

An import at module level, inside a function, or under ``TYPE_CHECKING`` is
collected by AST over every file of ``orpheus/``; the standard library and
``orpheus`` itself are excluded; each remaining top-level module is mapped to
its distribution (``importlib.metadata.packages_distributions``) and must be
named in ``[project].dependencies`` of ``pyproject.toml``. A function-local
import counts: it runs when the function does, so an install without the
distribution fails there, not at import (`[M]` 2026-10-02: ``sympy`` was
imported at run time by ``Quadrature.gauss_legendre`` and the symbolic
phase-space functions while listed only in the test and docs extras, 1 of 8
third-party modules undeclared).
"""

from __future__ import annotations

import ast
import re
import sys
import tomllib
from importlib.metadata import packages_distributions
from pathlib import Path

import pytest

pytestmark = pytest.mark.foundation

_ROOT = Path(__file__).resolve().parents[2]


def _imported_top_level_modules() -> dict[str, list[str]]:
    sites: dict[str, list[str]] = {}
    for path in sorted((_ROOT / "orpheus").rglob("*.py")):
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                names = [node.module]
            else:
                continue
            for name in names:
                top = name.split(".")[0]
                if top == "orpheus" or top in sys.stdlib_module_names:
                    continue
                sites.setdefault(top, []).append(f"{path.relative_to(_ROOT)}:{node.lineno}")
    return sites


def _declared() -> set[str]:
    project = tomllib.loads((_ROOT / "pyproject.toml").read_text())["project"]
    return {re.split(r"[<>=!~\[; ]", spec, maxsplit=1)[0].strip().lower().replace("_", "-") for spec in project["dependencies"]}


def test_s6_20_every_third_party_import_is_declared() -> None:
    sites = _imported_top_level_modules()
    print(f"S6.20: {len(sites)} third-party top-level modules imported under orpheus/: {sorted(sites)}")
    if "numpy" not in sites:
        pytest.fail(f"activation: the AST pass did not find numpy; it found {sorted(sites)}")
    declared = _declared()
    if "numpy" not in declared:
        pytest.fail(f"activation: the dependency parse did not find numpy; it found {sorted(declared)}")
    distributions = packages_distributions()
    undeclared = {}
    for module, where in sites.items():
        names = {d.lower().replace("_", "-") for d in distributions.get(module, [module])}
        if not names & declared:
            undeclared[module] = where[:3]
    if undeclared:
        pytest.fail(f"imported but not in [project].dependencies: {undeclared}")
