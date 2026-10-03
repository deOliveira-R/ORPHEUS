r"""The observables: the closed sum, its fields, its admission and its layer (#405 P2 step 3, R3.1-R3.3, R3.6).

Specified by the test-architect (2026-10-03, ``.claude/plans/reference_p2_spec.md``
§1.3), on the user's Q2 ruling of 2026-10-03: the primitive observable is
``FluxIntegral(weight)``, a linear functional of the scalar flux whose weight
is a mesh-free function over position and group; ``Rate(cells, weight)`` is a
spelling the SPECIFICATION resolves into a ``FluxIntegral`` and is not a
member. The closed sum is ``Observable = FluxIntegral | Ratio | Eigenvalue |
PointValue``. The values are physics-free and live in ``numerics``; whether a
weight fits a problem (its groups, its regions, the coordinates a symbolic
weight may read) is the specification's admission at READ time, never a
second admission here.

Every class is resolved on its module at run time (``_cls``), never captured
at collection, so a rebinding battery arm reaches every row (lessons
``L100``). Refusal rows pin the FIELD NAME as the message fragment (the
message is the gate: lessons §1); the module's wording beyond the field name
is its own.
"""

from __future__ import annotations

import ast
import dataclasses
import importlib
import inspect
import math
import os
import subprocess
import sys
import typing
from pathlib import Path
from typing import Any, assert_never

import numpy as np
import pytest

from orpheus.numerics.mesh_free_function import RegionwiseConstant, Symbolic
from tests.gates._content_identity_helpers import require

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/numerics/test_observable.py"
_ROOT = Path(__file__).resolve().parents[3]
_MODULE = "orpheus.numerics.observable"
_KINDS = ("FluxIntegral", "Ratio", "Eigenvalue", "PointValue")


def _mod() -> Any:
    return importlib.import_module(_MODULE)


def _cls(name: str) -> Any:
    """A class read off its module now (a rebinding arm reaches it)."""
    return getattr(_mod(), name)


def _table() -> RegionwiseConstant:
    return RegionwiseConstant(np.array([[1.0, 0.5], [0.0, 2.0]]))


def _flux() -> Any:
    return _cls("FluxIntegral")(_table())


def _point() -> Any:
    return _cls("PointValue")(0.75, 1)


# ═════════════════════════════════════════════════════════════════════════════
# R3.1 — the set is closed
# ═════════════════════════════════════════════════════════════════════════════


def _describe(observable: Any) -> str:
    """A test-side exhaustive ``match`` (pyright proves exhaustiveness; at run
    time a foreign object reaches ``case _``)."""
    m = _mod()
    match observable:
        case m.FluxIntegral():
            return "FluxIntegral"
        case m.Ratio():
            return "Ratio"
        case m.Eigenvalue():
            return "Eigenvalue"
        case m.PointValue():
            return "PointValue"
        case _:
            raise TypeError(f"not an observable: {observable!r}")


def test_r3_1_the_observable_set_is_closed() -> None:
    """``get_args(Observable)`` is exactly the four kinds, by name; the
    concrete ``ContentIdentity`` classes DEFINED in the module are exactly
    those four; one value of each reaches its own ``case``; a foreign value
    (a string, the class ``FluxIntegral`` itself, a ``RegionwiseConstant``)
    reaches ``case _``; every value is an instance of the union alias."""
    from orpheus.numerics.content import ContentIdentity

    m = _mod()
    members = {c.__name__ for c in typing.get_args(m.Observable)}
    require(members == set(_KINDS), f"Observable names {sorted(members)}")
    defined = {
        name for name, obj in vars(m).items()
        if inspect.isclass(obj) and issubclass(obj, ContentIdentity) and obj.__module__ == _MODULE
        and not inspect.isabstract(obj)
    }
    require(defined == set(_KINDS), f"content classes defined in the module: {sorted(defined)}")
    values = {"FluxIntegral": _flux(), "Ratio": m.Ratio(_flux(), _point()), "Eigenvalue": m.Eigenvalue(), "PointValue": _point()}
    for kind, value in values.items():
        require(_describe(value) == kind, f"{kind} reached the wrong case")
        require(isinstance(value, m.Observable), f"{kind} is not an instance of Observable")
    for foreign in ("FluxIntegral", m.FluxIntegral, _table()):
        with pytest.raises(TypeError, match="not an observable"):
            _describe(foreign)


def _check_exhaustive(observable: "Any") -> None:
    """Static leg for pyright: an ``assert_never`` over the alias."""
    from orpheus.numerics.observable import Eigenvalue, FluxIntegral, Observable, PointValue, Ratio

    o: Observable = observable
    if isinstance(o, (FluxIntegral, Ratio, Eigenvalue, PointValue)):
        return
    assert_never(o)


def test_r3_1_rate_is_not_a_numerics_value() -> None:
    """The user's Q2 ruling: ``Rate`` is the specification's resolving
    constructor, not an observable. By AST over every ``.py`` under
    ``orpheus/numerics`` (input count printed, `[M]` 62 on 9f8d1f7f, more than 50 asserted): no
    class ``Rate`` is defined; positive control: the census finds
    ``FluxIntegral`` in the observable module. And ``Rate`` is not an
    attribute of the module."""
    files = sorted((_ROOT / "orpheus" / "numerics").rglob("*.py"))
    print(f"R3.1: {len(files)} numerics files parsed")
    require(len(files) > 50, f"activation: only {len(files)} files parsed")
    found: dict[str, list[str]] = {}
    for path in files:
        for node in ast.walk(ast.parse(path.read_text(), filename=str(path))):
            if isinstance(node, ast.ClassDef):
                found.setdefault(node.name, []).append(str(path.relative_to(_ROOT)))
    require(found.get("FluxIntegral") == ["orpheus/numerics/observable.py"], f"control: FluxIntegral at {found.get('FluxIntegral')}")
    require("Rate" not in found, f"a class Rate is defined in numerics at {found.get('Rate')}")
    require(not hasattr(_mod(), "Rate"), "the observable module exposes a Rate")


# ═════════════════════════════════════════════════════════════════════════════
# R3.2 — the fields are the roles
# ═════════════════════════════════════════════════════════════════════════════

_FIELDS = {
    "FluxIntegral": ("weight",),
    "Ratio": ("numerator", "denominator"),
    "Eigenvalue": (),
    "PointValue": ("position", "group"),
}


@pytest.mark.parametrize("kind", _KINDS)
@pytest.mark.rests_on(f"{_HERE}::test_r3_1_the_observable_set_is_closed")
def test_r3_2_the_fields_are_the_roles(kind: str) -> None:
    """By name, resolved at run time: ``dataclasses.fields`` is exactly the
    roles; no field is typed ``bool`` (a flag is a missing type); every
    constructor parameter is required (a defaulted weight, operand, position
    or group would name a datum numerics cannot know)."""
    cls = _cls(kind)
    fields = dataclasses.fields(cls)
    names = tuple(f.name for f in fields)
    require(names == _FIELDS[kind], f"{kind} fields {names}, expected {_FIELDS[kind]}")
    typed_bool = [f.name for f in fields if f.type in (bool, "bool")]
    require(not typed_bool, f"{kind} has boolean fields {typed_bool}")
    defaulted = [p.name for p in inspect.signature(cls).parameters.values() if p.default is not inspect.Parameter.empty]
    require(not defaulted, f"{kind} defaults {defaulted}")


_MISSING_ROLE = (
    ("flux-integral-no-weight", "FluxIntegral", lambda: {}),
    ("ratio-one-operand", "Ratio", lambda: {"numerator": _flux()}),
    ("point-value-no-group", "PointValue", lambda: {"position": 0.5}),
)


@pytest.mark.parametrize("kind, kwargs", [r[1:] for r in _MISSING_ROLE], ids=[r[0] for r in _MISSING_ROLE])
def test_r3_2_a_missing_role_does_not_construct(kind: str, kwargs: Any) -> None:
    """``TypeError`` from the constructor's own signature (R3.2's no-default leg, read at the call)."""
    with pytest.raises(TypeError):
        _cls(kind)(**kwargs())


# ═════════════════════════════════════════════════════════════════════════════
# R3.3 — admission, eager, each refusal naming its field
# ═════════════════════════════════════════════════════════════════════════════

_BAD_WEIGHTS = (
    ("ndarray", lambda: np.ones((2, 2))),
    ("float", lambda: 1.0),
    ("none", lambda: None),
    ("nested-list", lambda: [[1.0, 0.5], [0.0, 2.0]]),
    ("an-observable", lambda: _point()),
)


@pytest.mark.parametrize("make", [m for _, m in _BAD_WEIGHTS], ids=[n for n, _ in _BAD_WEIGHTS])
def test_r3_3_the_weight_is_a_mesh_free_function(make: Any) -> None:
    """``TypeError`` naming the weight and what it must be: an extensional
    array is mesh-bound, a number has no position or group, and an observable
    is not a function. The fragment is ``the weight is a mesh-free function``,
    not the bare field name: the encoder's own refusal of a writeable array
    names its path ``FluxIntegral.weight`` too, so a bare ``weight`` keeps the
    array rows green with the type check removed (``[M]`` battery arm A1)."""
    with pytest.raises(TypeError, match="the weight is a mesh-free function"):
        _cls("FluxIntegral")(make())


@pytest.mark.parametrize("make", [_table, lambda: Symbolic.of(1 + Symbolic.r, Symbolic.r**2)], ids=["regionwise", "symbolic"])
def test_r3_3_either_mesh_free_function_is_a_weight(make: Any) -> None:
    """Positive legs (anti-#11): both forms construct and are kept."""
    weight = make()
    flux = _cls("FluxIntegral")(weight)
    require(flux.weight == weight, "the weight is kept")


_BAD_OPERANDS = (
    ("float", lambda: 1.0),
    ("a-mesh-free-function", _table),
    ("none", lambda: None),
    ("an-enclosure", lambda: importlib.import_module("orpheus.numerics.enclosure").Enclosure(1.0, 0.0)),
)


@pytest.mark.parametrize("side", ["numerator", "denominator"])
@pytest.mark.parametrize("make", [m for _, m in _BAD_OPERANDS], ids=[n for n, _ in _BAD_OPERANDS])
def test_r3_3_a_ratio_is_of_observables(make: Any, side: str) -> None:
    """``TypeError`` naming the side: a ratio of a reading (an ``Enclosure``)
    or of a number is a quotient of answers, not an observable."""
    operands = {"numerator": _flux(), "denominator": _point()}
    operands[side] = make()
    with pytest.raises(TypeError, match=side):
        _cls("Ratio")(**operands)


def test_r3_3_ratios_nest_and_any_two_observables_divide() -> None:
    """Positive legs: a ratio of a ratio, a ratio with the eigenvalue, and the
    degenerate ``Ratio(x, x)`` (a value; whether an answer can read it is the
    answer's question) all construct."""
    R, E = _cls("Ratio"), _cls("Eigenvalue")
    R(R(_flux(), _point()), E())
    R(E(), _flux())
    R(_flux(), _flux())


_BAD_POINTS = (
    ("nan-position", (math.nan, 0), ValueError, "position"),
    ("inf-position", (math.inf, 0), ValueError, "position"),
    ("str-position", ("0.5", 0), TypeError, "position"),
    ("complex-position", (0.5j, 0), TypeError, "position"),
    ("negative-group", (0.5, -1), ValueError, "group"),
    ("float-group", (0.5, 1.0), TypeError, "group"),
    ("bool-group", (0.5, True), TypeError, "group"),
    ("none-group", (0.5, None), TypeError, "group"),
)


@pytest.mark.parametrize("args, error, fragment", [r[1:] for r in _BAD_POINTS], ids=[r[0] for r in _BAD_POINTS])
def test_r3_3_a_point_value_is_a_finite_position_and_a_group_index(args: tuple[Any, Any], error: type[BaseException], fragment: str) -> None:
    """A position is a finite real (a coordinate of the geometry, resolved by
    the specification at read time); a group is a non-negative ``int`` index
    (a ``bool`` or a float is refused, as ``parse_integer`` refuses them)."""
    with pytest.raises(error, match=fragment):
        _cls("PointValue")(*args)


def test_r3_3_a_point_value_admits_an_integer_position_and_group_zero() -> None:
    """Positive legs: an integer position is stored as a double; group 0 and
    a negative position (a slab coordinate) are admitted."""
    p = _cls("PointValue")(1, 0)
    require(type(p.position) is float and p.position == 1.0, f"position stored as {p.position!r}")
    require(type(p.group) is int and p.group == 0, f"group stored as {p.group!r}")
    _cls("PointValue")(-2.5, 3)


# ═════════════════════════════════════════════════════════════════════════════
# R3.6 — the layer
# ═════════════════════════════════════════════════════════════════════════════

_LAYER_SCRIPT = """
import sys
import numpy as np
import orpheus.numerics.observable as m
from orpheus.numerics.mesh_free_function import RegionwiseConstant
f = m.FluxIntegral(RegionwiseConstant(np.ones((2, 1))))
m.Ratio(f, m.PointValue(0.5, 0)); m.Eigenvalue()
print(m.__file__)
print(sorted({k.split('.')[1] for k in sys.modules if k.startswith('orpheus.')}))
print('sympy' in sys.modules)
"""


def test_r3_6_the_module_is_numerics() -> None:
    """(a) A fresh interpreter importing the module and building one value of
    each kind (with a ``RegionwiseConstant`` weight) loads only the
    ``orpheus`` packages a cold ``import orpheus.numerics`` loads
    (``{geometry, numerics}``) and no SymPy. (b) By AST every ``orpheus``
    import of the module (module level, local, ``TYPE_CHECKING``) is under
    ``orpheus.numerics``; activation: the import of
    ``orpheus.numerics.mesh_free_function`` is seen. (c) The module is a cold
    entry point of ``tests/gates/test_layer_imports.py``."""
    env = {**os.environ, "PYTHONPATH": str(_ROOT)}
    out = subprocess.run([sys.executable, "-O", "-c", _LAYER_SCRIPT], cwd=_ROOT, env=env, capture_output=True, text=True)
    require(out.returncode == 0, out.stderr)
    file, packages, sympy_loaded = out.stdout.strip().splitlines()
    require(file.startswith(str(_ROOT)), f"the subprocess imported {file} (L22)")
    require(packages == "['geometry', 'numerics']", f"loaded orpheus packages {packages}")
    require(sympy_loaded == "False", "building the values loads SymPy")
    module_file = _mod().__file__
    require(module_file is not None, "the module has no file")
    imported: list[str] = []
    for node in ast.walk(ast.parse(Path(str(module_file)).read_text())):
        if isinstance(node, ast.ImportFrom) and node.module:
            imported.append(node.module)
        elif isinstance(node, ast.Import):
            imported.extend(alias.name for alias in node.names)
    ours = [m for m in imported if m.startswith("orpheus")]
    require("orpheus.numerics.mesh_free_function" in ours, f"activation: the mesh-free import is not seen in {ours}")
    outside = [m for m in ours if not m.startswith("orpheus.numerics")]
    require(not outside, f"imports outside numerics: {outside}")
    from tests.gates import test_layer_imports as lint

    marks = getattr(lint.test_entry_point_imports_in_a_fresh_interpreter, "pytestmark")
    entries = next(mk for mk in marks if mk.name == "parametrize").args[1]
    require(_MODULE in entries, f"{_MODULE} is not a cold entry point")
