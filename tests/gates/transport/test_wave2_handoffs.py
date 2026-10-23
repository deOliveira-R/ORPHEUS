r"""The hand-offs of slices S1 and S2a that slice S2b closes (wave 2; plan section 6.2, "Handed over by S1" and "by S2a").

Source-level gates: they read files (``ast``) and import no production module, so each fails at the base tip for its NAMED
reason. Each scanner is validated by a planted control before its zero is believed (``instrument-doctrine`` X1).

* S2a hand-off (a): the S2a gate files read the package through a ``sys.modules`` namespace stub and a ``_Lazy`` loader in
  ``tests/gates/transport/_boundary_conditions_helpers.py``, because ``orpheus.transport/__init__`` imported the retired
  ``BC``. Once ``transport/method.py`` stops importing it the real parent imports, the stub is dead, and it is the cause of
  the 67 pyright errors the S2a gate files carry (a class subclassing ``bc.ResponseRealizer`` through an ``Any``-typed ``bc``
  cannot take the ``table=`` and ``realizes=`` keywords). Done-when: no ``sys.modules[...]`` store, no ``_Lazy``, no ``load``
  in the helper, and ``bc`` is a plain module import.
* WF2.2 rider (a): ``tests/gates/numerics/test_ordinate_pairing.py`` no longer imports a gate of ``tests/gates/geometry/``
  (it imported the poison fixtures from ``test_bc_universal_invariants.py``, which WF2.4 deletes with the laws it tests).
  Done by the move of the fixtures into the numerics file in the gates commit of S2b; the gate keeps it so.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

pytestmark = pytest.mark.foundation

_THIS = "tests/gates/transport/test_wave2_handoffs.py"
_TESTS = Path(__file__).resolve().parent.parent  # tests/gates


def stub_findings(source: str) -> list[str]:
    """What makes a helper a namespace stub: a store into ``sys.modules[...]``, a ``_Lazy`` class, a ``load`` function, and
    a ``bc`` that is not bound by a plain ``import orpheus.transport.boundary_conditions as bc``."""
    tree = ast.parse(source)
    found = []
    plain_bc = False
    for node in ast.walk(tree):
        if isinstance(node, ast.Subscript) and isinstance(node.ctx, ast.Store) and ast.unparse(node.value) == "sys.modules":
            found.append(f"{node.lineno}: a store into sys.modules")
        elif isinstance(node, ast.ClassDef) and node.name == "_Lazy":
            found.append(f"{node.lineno}: class _Lazy")
        elif isinstance(node, ast.FunctionDef) and node.name == "load":
            found.append(f"{node.lineno}: def load")
    for node in tree.body:
        if isinstance(node, ast.Import):
            plain_bc = plain_bc or any(a.name == "orpheus.transport.boundary_conditions" and a.asname == "bc" for a in node.names)
    if not plain_bc:
        found.append("bc is not bound by `import orpheus.transport.boundary_conditions as bc`")
    return found


def geometry_gate_imports(source: str) -> list[str]:
    """Every import (at any depth) of a module under ``tests.gates.geometry``."""
    found = []
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.ImportFrom) and (node.module or "").startswith("tests.gates.geometry"):
            found.append(ast.unparse(node))
        elif isinstance(node, ast.Import):
            found += [ast.unparse(node) for a in node.names if a.name.startswith("tests.gates.geometry")]
    return found


def test_the_scanners_see_what_they_must() -> None:
    r"""The positive controls: a planted stub (store, ``_Lazy``, ``load``, a lazy ``bc``) reads four findings and a plain
    helper reads none; a planted ``from tests.gates.geometry import test_bc_universal_invariants`` is found in a function
    body, a plain ``from tests.gates.numerics import x`` is not. Passes today; reddens on a scanner that reads nothing."""
    stub = (
        "import sys, types\n"
        "class _Lazy:\n    pass\n"
        "def load():\n    sys.modules['orpheus.transport'] = types.ModuleType('x')\n"
        "bc = _Lazy()\n"
    )
    assert len(stub_findings(stub)) == 4
    assert stub_findings("import orpheus.transport.boundary_conditions as bc\nx = 1\n") == []
    assert geometry_gate_imports("def f():\n    from tests.gates.geometry import test_bc_universal_invariants as m\n")
    assert not geometry_gate_imports("from tests.gates.numerics import x\n")


@pytest.mark.rests_on(f"{_THIS}::test_the_scanners_see_what_they_must")
def test_the_boundary_conditions_helper_has_no_stub() -> None:
    r"""S2a hand-off (a): ``_boundary_conditions_helpers.py`` reads the package by a plain import. Red input at the base
    tip: the helper holds ``load()`` (a ``sys.modules`` store of a bare namespace parent), ``_Lazy`` and ``bc = _Lazy()``.
    Mutation: the stub restored in a copy."""
    path = Path(__file__).with_name("_boundary_conditions_helpers.py")
    findings = stub_findings(path.read_text())
    assert not findings, f"{path.name}: {findings}"


@pytest.mark.rests_on(f"{_THIS}::test_the_scanners_see_what_they_must")
def test_the_numerics_pairing_gate_imports_no_geometry_gate() -> None:
    r"""WF2.2 rider (a): the numerics gate stands alone. Green today (the fixtures moved with S2b's gates commit); it
    reddens on the import coming back. Mutation: the lazy ``_poison_module()`` restored."""
    path = _TESTS / "numerics" / "test_ordinate_pairing.py"
    assert path.exists()
    assert not geometry_gate_imports(path.read_text())
