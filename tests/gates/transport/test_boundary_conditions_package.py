r"""The package as a unit: its file set, its facade, its import hygiene (wave 2, slice S2a, gate 12).

Spec: ``.claude/plans/wave2_boundary_carve.md`` section 2.4 (the module table), 5.0, law L-H (import hygiene), 6.2 S2a
gate 12 and the brief's "green in isolation".

The claims (``vv-principles``), all STRUCTURE (an AST pass over the package), each with a planted positive control
(X1: a scan over a list finds nothing in a clean tree and in a broken instrument alike):

* ``test_the_package_imports_nothing_across`` (gate 12): no module of ``orpheus/transport/boundary_conditions/``
  imports ``orpheus.specification``, ``orpheus.derivations``, ``orpheus.geometry.boundary`` (the old vocabulary this
  package replaces; WF2.4 deletes it), any method package (``sn``, ``diffusion``, ``cp``, ``moc``, ``mc``,
  ``kinetics``, ``fuel``, ``thermal_hydraulics``, ``homogeneous``), nor the two modules of the red chain
  (``orpheus.transport.method`` and ``orpheus.transport.boundary_realization``, which import the retired ``BC``
  until slice S2b re-spells them). Imports at ANY nesting (function-local, ``if TYPE_CHECKING:``) and relative
  imports are resolved and read. The positive population: the package does import ``orpheus.numerics``,
  ``orpheus.values`` and the L1 geometry types (``orpheus.geometry.deck`` or ``colouring``), for real.
* ``test_the_package_has_exactly_the_tabled_modules``: the file set is the nine modules of the table, and the
  ``__init__`` holds re-exports and nothing else (no residue, no extra module).
* ``test_the_facade_exports_what_the_tables_name``: every name of the tables is exported by the package, defined in
  one of its modules (not in ``__init__``), and ``__all__`` lists them.
* ``test_stamp_boundary_role_is_copied_in``: ``stamp_boundary_role`` lives here (the old one in
  ``boundary_realization.py`` was deleted in WF2.4, with the row that held the two bodies equal) and stamps the
  ``BOUNDARY`` block role on the instance and returns it.

First red, at the base tip ``4ec95bc6``: the package directory does not exist, so the file set, the scan and the
facade read ``FileNotFoundError`` / ``ModuleNotFoundError`` at their first line; each scan's CONTROLS run first and
are green (they read planted sources), so the red is the real leg's alone.

Mutations (on a COPY of the package): add ``from orpheus.derivations.common import x`` in ``network.py``, a
function-local ``from orpheus.sn import problem`` in ``realization.py``, ``from orpheus.geometry.boundary import BC``
in ``condition.py``, ``from orpheus.transport.method import TransportMethod`` in ``binding.py`` (the scan names
each); add a module ``_helpers.py``; put a function in ``__init__``.
"""

from __future__ import annotations

import ast
import types

import pytest

from tests.gates._content_identity_helpers import require
from tests.gates.transport._boundary_conditions_helpers import (
    MODULES,
    PACKAGE,
    PACKAGE_DIR,
    bc,
    imports_of,
    node,
    violations_of_imports,
)

pytestmark = pytest.mark.foundation

_HERE = "test_boundary_conditions_package.py"
_FILES = node(_HERE, "test_the_package_has_exactly_the_tabled_modules")

#: The modules (and trees) the package may not import, with the reason of each.
FORBIDDEN = (
    "orpheus.specification",
    "orpheus.derivations",
    "orpheus.geometry.boundary",
    "orpheus.sn", "orpheus.diffusion", "orpheus.cp", "orpheus.moc", "orpheus.mc",
    "orpheus.kinetics", "orpheus.fuel", "orpheus.thermal_hydraulics", "orpheus.homogeneous",
    "orpheus.transport.method", "orpheus.transport.boundary_realization",
)

#: Every name the tables of section 2.4 name, by module (the spec's own words), plus the two table errors of the
#: realization name check and the aliases the gates read.
EXPORTS: dict[str, tuple[str, ...]] = {
    "condition": ("BoundaryCondition", "WallCondition", "WallDisplacement"),
    "network": ("MirrorPermutation", "CurrentRetraction", "IsotropicSection", "Factor", "SpecularChain", "DiffuseChain", "Chain",
                "CHAIN_ORDER", "chain_composes", "ZeroNetwork", "ResponseNetwork"),
    "inflow": ("InflowSampler", "FaceFrameInflow", "LabFrameInflow", "InflowTerms"),
    "errors": ("BoundaryError", "SubmarkovViolationError", "BoundaryResponseNotPositiveError", "BoundaryRefused", "UnboundTagError",
               "UnusedConditionError", "NonInvariantConditionError"),
    "declared": ("BoundaryConditions",),
    "binding": ("BoundaryCarrier", "Binding", "bind"),
    "admission": ("Realized", "Approximated", "Refused", "RefusalAxis", "Admission", "FaceSite", "SiteClass"),
    "realization": ("StateCone", "RealizerTable", "ResponseRealizer", "SourceRealizer", "DeckRealizer", "sample_terms",
                    "stamp_boundary_role", "TableNotSealedError", "SealedTableError"),
}

_PLANTED = {
    "absolute, in a specification": ("from orpheus.specification import exterior\n", ["orpheus.specification", "orpheus.specification.exterior"]),
    "import statement, a derivation": ("import orpheus.derivations.common.xs_library\n", ["orpheus.derivations.common.xs_library"]),
    "relative, climbing to derivations": ("from ...derivations import common\n", ["orpheus.derivations", "orpheus.derivations.common"]),
    "function-local, a method package": ("def f():\n    from orpheus.sn import problem\n    return problem\n", ["orpheus.sn", "orpheus.sn.problem"]),
    "under TYPE_CHECKING, the old vocabulary": (
        "from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n    from orpheus.geometry.boundary import BC\n",
        ["orpheus.geometry.boundary", "orpheus.geometry.boundary.BC"]),
    "the red chain": ("from orpheus.transport.method import TransportMethod\n", ["orpheus.transport.method", "orpheus.transport.method.TransportMethod"]),
}
_PLANTED_CLEAN = (
    "from orpheus.geometry.deck import MirrorDeck\n"
    "from orpheus.geometry import deck\n"
    "from orpheus.numerics.spaces.angular_trace_space import TraceRole\n"
    "from orpheus.snake import x\n"
    "from orpheus.values.content import ContentIdentity\n"
    "from .errors import BoundaryError\n"
    "import numpy as np\n"
)


def _package_sources() -> dict[str, str]:
    files = sorted(PACKAGE_DIR.glob("*.py"))
    require(PACKAGE_DIR.is_dir(), f"{PACKAGE_DIR} is not a directory")
    return {path.stem: path.read_text(encoding="utf-8") for path in files}


# ── gate 12: nothing crosses ─────────────────────────────────────────────────


def test_the_package_imports_nothing_across() -> None:
    """The package is green in isolation: it imports numerics, values and the L1 geometry types, and none of the
    forbidden trees (the module docstring lists them with their reason).

    Controls (the scanner must read each shape, X1): every planted shape of ``_PLANTED`` is reported with exactly the
    names listed (an absolute ``from``, a plain ``import``, a RELATIVE ``from ...derivations``, a function-local
    import, an import under ``TYPE_CHECKING``, the red chain); a clean planted source, including a name that merely
    STARTS with a forbidden one (``orpheus.snake``), is not reported.

    Population: the ten modules are read, and between them they import ``orpheus.numerics``, ``orpheus.values`` and
    an L1 geometry module (``deck`` or ``colouring``) (so a scan that reads nothing, or a path that resolves to an
    empty directory, cannot pass).

    Reddens on: the four planted imports of the module docstring, each in its own copy.
    """
    for label, (source, expected) in _PLANTED.items():
        found = violations_of_imports(source, FORBIDDEN)
        require(found == sorted(expected), f"control {label!r}: the scanner reports {found}, expected {sorted(expected)}")
    require(violations_of_imports(_PLANTED_CLEAN, FORBIDDEN) == [], "control: the scanner flags a clean source (orpheus.snake is not orpheus.sn)")
    require(imports_of("from .errors import X\n", f"{PACKAGE}.binding") == [f"{PACKAGE}.errors", f"{PACKAGE}.errors.X"], "control: a sibling relative import does not resolve to the package")

    sources = _package_sources()
    require(len(sources) == len(MODULES), f"the scan read {len(sources)} modules: {sorted(sources)}")
    all_imports: set[str] = set()
    problems: dict[str, list[str]] = {}
    for name, source in sources.items():
        module_name = f"{PACKAGE}.{name}" if name != "__init__" else f"{PACKAGE}.__init__"
        all_imports |= set(imports_of(source, module_name))
        found = violations_of_imports(source, FORBIDDEN, module_name)
        if found:
            problems[name] = found
    require(not problems, f"the package imports across its boundary: {problems}")
    for needed in ("orpheus.numerics", "orpheus.values"):
        require(any(i == needed or i.startswith(needed + ".") for i in all_imports), f"the package imports nothing from {needed}")
    require(any(i.startswith(("orpheus.geometry.deck", "orpheus.geometry.colouring")) for i in all_imports), "the package imports no L1 geometry type")


def test_the_package_has_exactly_the_tabled_modules() -> None:
    """The ten modules of section 2.4's table, no others (``_helpers.py`` is residue); ``__init__`` holds a docstring,
    the re-exporting ``from ... import`` statements and ``__all__``, nothing else (no function, no class, no
    assignment but ``__all__``).

    Controls: a planted ``__init__`` with a function body is reported; a planted one with only re-exports is clean.
    """
    def residue(source: str) -> list[str]:
        out = []
        for stmt in ast.parse(source).body:
            if isinstance(stmt, ast.ImportFrom):
                continue
            if isinstance(stmt, ast.Expr) and isinstance(stmt.value, ast.Constant) and isinstance(stmt.value.value, str):
                continue
            if isinstance(stmt, ast.Assign | ast.AnnAssign) and any(isinstance(t, ast.Name) and t.id == "__all__" for t in (stmt.targets if isinstance(stmt, ast.Assign) else [stmt.target])):
                continue
            out.append(type(stmt).__name__)
        return out

    require(residue("def helper():\n    return 1\n") == ["FunctionDef"], "control: a function in __init__ is not reported")
    require(residue('"""d"""\nfrom .a import B\n__all__ = ["B"]\n') == [], "control: a pure re-export __init__ is reported")
    require(residue("from .a import B\nX = 1\n") == ["Assign"], "control: an assignment in __init__ is not reported")

    present = sorted(path.stem for path in PACKAGE_DIR.glob("*.py"))
    require(present == sorted(MODULES), f"the package holds {present}, the table names {sorted(MODULES)}")
    require(not residue((PACKAGE_DIR / "__init__.py").read_text(encoding="utf-8")), f"__init__ holds more than re-exports: {residue((PACKAGE_DIR / '__init__.py').read_text(encoding='utf-8'))}")
    require(not [p for p in PACKAGE_DIR.iterdir() if p.is_dir() and p.name != "__pycache__"], "the package holds a subdirectory")


@pytest.mark.rests_on(_FILES)
def test_the_facade_exports_what_the_tables_name() -> None:
    """Every name of the tables is reachable from the package, was DEFINED in the module the table gives it
    (``__module__``; a class by identity, an alias or a function by name), and is listed in ``__all__``.

    Reddens on: a module missing a name of its table; a name defined in ``__init__``; an ``__all__`` that lags the
    facade.
    """
    package = bc  # the lazy read; every attribute access below imports the package
    names = [(module, name) for module, listed in EXPORTS.items() for name in listed]
    require(len(names) == 45, f"the export table lists {len(names)} names, not the 45 it was written with")
    exported = set(package.__all__)
    for module, name in names:
        require(hasattr(package, name), f"the package does not export {name} ({module}.py)")
        require(name in exported, f"{name} is not in __all__")
        member = getattr(package, name)
        owner = getattr(member, "__module__", None)
        if isinstance(member, type) or isinstance(member, types.FunctionType):
            require(owner == f"{PACKAGE}.{module}", f"{name} is defined in {owner}, the table places it in {module}.py")


def test_stamp_boundary_role_is_copied_in() -> None:
    """``stamp_boundary_role(op)`` stamps the ``BOUNDARY`` block role on the INSTANCE and returns it, unchanged
    behaviour of the old helper (``orpheus/transport/boundary_realization.py``, deleted in WF2.4).

    Reddens on: a copy that returns a stamped clone; one that stamps ``SOURCE``; one that does not return.
    """
    from orpheus.numerics.operator import BlockRole

    op = types.SimpleNamespace(block_role=None)
    stamped = bc.stamp_boundary_role(op)  # type: ignore[arg-type]  # a duck-typed stand-in operator
    require(stamped is op, "stamp_boundary_role returns the instance it was given")
    require(op.block_role is BlockRole.BOUNDARY, f"the block role is {op.block_role}")
    require(bc.stamp_boundary_role.__module__ == f"{PACKAGE}.realization", "stamp_boundary_role is defined in realization.py")
