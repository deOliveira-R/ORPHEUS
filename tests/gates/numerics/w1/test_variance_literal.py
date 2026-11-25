r"""G03j: every factor states its variance (N03, spec 7.4).

``Factor(axis, sk=...)`` without a variance raises ``TypeError``; over ``orpheus/``, by ``ast`` and keyed on what the name RESOLVES to (a name bound by an
import of ``numerics.factor``'s ``Factor``, not the spelling: ``transport`` has an alias of the same spelling until it is renamed), ``Factor.on`` is defined
nowhere and every call ``Factor(...)`` passes a variance (a keyword, or three positional arguments).

The ``ast`` leg is green at the base tip (no module imports ``numerics.factor.Factor`` there). The positive control: the same predicate reads a
synthetic source that calls it without a variance, and a source that states one. The mutation row applies a default to a copy of ``factor.py``'s source
and requires the gate's own reading to turn red (either at class creation or at the first assertion, ``[M]`` the two red lines differ).

Exact on any host: parse trees and a ``TypeError``. First red at ``ff305e0f``: import-only for the ``TypeError`` leg.
"""

from __future__ import annotations

import ast
import inspect
import sys
import types
from pathlib import Path

import pytest

from tests.gates.numerics.w1 import _factor_probe as P

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _resolving_names(tree: ast.AST) -> tuple[set[str], set[str]]:
    """``(names bound to numerics.factor.Factor, names bound to the module numerics.factor)`` in one parse tree."""
    classes: set[str] = set()
    modules: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            mod = node.module or ""
            if mod.endswith("numerics.factor") or (node.level > 0 and mod == "factor"):
                classes |= {a.asname or a.name for a in node.names if a.name == "Factor"}
            if mod.endswith("numerics") or (node.level > 0 and mod == ""):
                modules |= {a.asname or a.name for a in node.names if a.name == "factor"}
        elif isinstance(node, ast.Import):
            for a in node.names:
                if a.name.endswith("numerics.factor"):
                    modules.add(a.asname or a.name)
    return classes, modules


def violations_in(source: str) -> list[str]:
    """Calls of ``numerics.factor.Factor`` that state no variance, and every ``Factor.on`` read, in one source text."""
    tree = ast.parse(source)
    classes, modules = _resolving_names(tree)
    out: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            f = node.func
            is_factor = (isinstance(f, ast.Name) and f.id in classes) or (
                isinstance(f, ast.Attribute) and f.attr == "Factor" and isinstance(f.value, (ast.Name, ast.Attribute)) and ast.unparse(f.value) in modules | {"orpheus.numerics.factor"}
            )
            if is_factor and not (len(node.args) >= 3 or any(k.arg == "variance" for k in node.keywords)):
                out.append(f"line {node.lineno}: a Factor call states no variance")
        if isinstance(node, ast.Attribute) and node.attr == "on" and isinstance(node.value, ast.Name) and node.value.id in classes:
            out.append(f"line {node.lineno}: Factor.on is read")
        if isinstance(node, ast.ClassDef) and node.name == "Factor":
            for b in node.body:
                if isinstance(b, ast.FunctionDef) and b.name == "on" and any(isinstance(d, ast.Name) and d.id == "classmethod" for d in b.decorator_list):
                    out.append(f"line {b.lineno}: Factor.on is defined as a classmethod")
    return out


def _orpheus_root() -> Path:
    import orpheus

    return Path(orpheus.__file__).parent


def _gate_reading(text: str) -> str | None:
    """The gate's own reading of a ``factor.py`` source: ``None`` when ``Factor(axis, sk=...)`` raises ``TypeError``, else the red it reads."""
    mod = types.ModuleType("orpheus.numerics.factor_copy")
    mod.__package__ = "orpheus.numerics"
    sys.modules[mod.__name__] = mod
    try:
        try:
            exec(compile(text, "factor_copy.py", "exec"), mod.__dict__)
        except Exception as e:  # a default before a non-default field reds at class creation
            return f"red at class creation: {type(e).__name__}: {str(e)[:80]}"
        try:
            mod.Factor(P.diag_axis([1.0, 2.0]), sk=mod.Subquotient.whole(2))
        except TypeError:
            return None
        return "red at the assertion: Factor constructed without a variance"
    finally:
        sys.modules.pop(mod.__name__, None)


def test_every_factor_states_its_variance() -> None:
    from orpheus.numerics import factor as F

    try:
        F.Factor(P.diag_axis([1.0, 2.0]), sk=F.Subquotient.whole(2))
    except TypeError:
        pass
    else:
        pytest.fail("Factor(axis, sk=...) without a variance was constructed")
    files = sorted(p for p in _orpheus_root().rglob("*.py"))
    _require(len(files) > 200, f"the census reads {len(files)} files: the filter dropped the tree")
    found = [f"{p.relative_to(_orpheus_root())}: {v}" for p in files for v in violations_in(p.read_text())]
    _require(not found, f"{len(found)} of {len(files)} files carry a Factor with no variance or a Factor.on: {found[:5]}")


def test_the_ast_predicate_can_fail() -> None:
    """The positive control of the ``ast`` leg: a call without a variance and a ``Factor.on`` are read; a call stating one is not."""
    head = "from orpheus.numerics.factor import Factor\n"
    _require(violations_in(head + "f = Factor(axis, sk=s)\n") != [], "a Factor call with no variance is not read")
    _require(violations_in(head + "f = Factor(axis, v, s)\n") == [], "a Factor call stating its variance is read as a violation")
    _require(violations_in(head + "f = Factor(axis, variance=v, sk=s)\n") == [], "a keyword variance is read as a violation")
    _require(violations_in(head + "f = Factor.on(axis, s)\n") != [], "Factor.on is not read")
    _require(violations_in("from orpheus.transport.boundary_conditions.network import Factor\nf = Factor(a, b)\n") == [], "the transport alias is read as the numerics Factor")


def test_a_default_on_variance_reddens_the_gate() -> None:
    """Mutation: a default on ``variance`` (and on ``sk``) in a copy of ``factor.py``; the gate's own reading must be red, at class creation or at the assertion."""
    from orpheus.numerics import factor as F

    src = inspect.getsource(F)
    _require(_gate_reading(src) is None, f"the unmutated copy reads red: {_gate_reading(src)}")
    old = "    variance: FormlessVariance | Variance\n    sk: Subquotient"
    _require(old in src, "the mutation target (the two Factor fields) is not in factor.py as the spec states them")
    default = "Variance(Kind.MEASURE, VarianceChart.DENSITY)"  # a legal variance on the formed axis the gate builds
    only_variance = src.replace(old, f"    variance: FormlessVariance | Variance = {default}\n    sk: Subquotient")
    both = src.replace(old, f"    variance: FormlessVariance | Variance = {default}\n    sk: Subquotient = None  # type: ignore[assignment]")
    _require(_gate_reading(only_variance) is not None, "a default on variance alone is not read")
    _require(_gate_reading(both) is not None, "a default on variance and on sk is not read")
