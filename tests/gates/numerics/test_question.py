r"""The laws of the question values (#405 P1 step 7, S7.1-S7.8 and S7.13).

DRAFT (test-architect, 2026-10-02). Lands as ``tests/gates/numerics/test_question.py``
with the module ``orpheus/numerics/question.py``.

A question is physics-free: ``Eigen(parameter, point, mode)`` asks where, on
the line through ``point`` along the direction an opaque ``parameter`` key
names, the system is singular; ``FixedSource(source, point)`` asks for
``E(point)^{-1} q``; ``Response(detector, point)`` asks for the importance of
a detector, ``E(point)^{-dagger} R``. The point is a ``FrozenMapping`` from
parameter key to a real offset from the physical value, empty by default
(the physical point). The mode is ``Fundamental()`` or ``Nearest(tau)``.
No value carries a forward/adjoint flag: the question's TYPE is its role
(ruling 3 of 2026-10-02), and the eigen adjoint belongs to the eigen answer.

The content-identity rows are ``test_content_identity_question.py`` (S7.9-S7.12).
The specification is ``.claude/plans/reference_p1_spec.md`` §1.7.
"""

from __future__ import annotations

import ast
import dataclasses
import inspect
import math
import os
import subprocess
import sys
import textwrap
import typing
from pathlib import Path
from typing import Any, assert_never

import numpy as np
import pytest

import orpheus.numerics.question as question_module
from orpheus.numerics.content import ContentIdentity, ContentlessError, FrozenMapping
from orpheus.numerics.mesh_free_function import RegionwiseConstant, Symbolic
from orpheus.numerics.question import (
    Eigen,
    FixedSource,
    Fundamental,
    Mode,
    Nearest,
    Question,
    Response,
)
from tests.gates._content_identity_helpers import require

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/numerics/test_question.py"
_MF = "tests/gates/numerics/test_mesh_free_function.py"
_CI = "tests/gates/numerics/test_content_identity.py"
_ROOT = Path(__file__).resolve().parents[3]

#: An opaque parameter key: numerics never reads what it names (ruling 5).
_KEY = "fission-emission"


def _table() -> RegionwiseConstant:
    return RegionwiseConstant(np.array([[1.0, 0.5], [0.0, 2.0]]))


def _symbolic() -> Symbolic:
    return Symbolic.of(1 + Symbolic.mu, Symbolic.r**2)


def _kind(q: Question) -> str:
    """The exhaustive dispatch over the closed set: pyright proves every case is
    named (``assert_never``), and at run time a foreign object reaches ``case _``."""
    match q:
        case Eigen():
            return "eigen"
        case FixedSource():
            return "fixed source"
        case Response():
            return "response"
        case _:
            assert_never(q)


def _mode_kind(m: Mode) -> str:
    match m:
        case Fundamental():
            return "fundamental"
        case Nearest():
            return "nearest"
        case _:
            assert_never(m)


# ── S7.1: the closed set ─────────────────────────────────────────────────────


def test_s7_1_the_question_set_is_closed() -> None:
    """``Question`` is exactly the three kinds and ``Mode`` exactly the two
    selectors, read off the union aliases (the population is the TYPE, X2); the
    module defines no other content value; the exhaustive match names each kind
    once and a foreign object (the pre-posing ``EigenPosing``, a string) falls
    through to ``case _``. Adding ``time``, ``Evolution`` or ``Enclosed`` later
    is an extension that edits this row on purpose (ruling 4, ruling 6)."""
    require(set(typing.get_args(Question)) == {Eigen, FixedSource, Response},
            f"Question = {typing.get_args(Question)}")
    require(set(typing.get_args(Mode)) == {Fundamental, Nearest}, f"Mode = {typing.get_args(Mode)}")
    defined = {
        obj for obj in vars(question_module).values()
        if inspect.isclass(obj) and issubclass(obj, ContentIdentity) and obj.__module__ == question_module.__name__
    }
    require(defined == {Eigen, FixedSource, Response, Fundamental, Nearest},
            f"the module defines the content classes {sorted(c.__qualname__ for c in defined)}")
    kinds = {_kind(Eigen(_KEY)), _kind(FixedSource(_table())), _kind(Response(_table()))}
    require(kinds == {"eigen", "fixed source", "response"}, f"the match reached {kinds}")
    require({_mode_kind(Fundamental()), _mode_kind(Nearest(0.5))} == {"fundamental", "nearest"}, "mode match")
    from orpheus.numerics.posing import EigenPosing

    for foreign in ("eigen", EigenPosing):
        with pytest.raises(AssertionError):
            _kind(foreign)  # type: ignore[arg-type]


_STRUCK = ("CriticalParameter",)
_DEFERRED = ("Enclosed", "Evolution", "Alpha", "Time", "Index", "Pseudospectrum")


def test_s7_2_the_struck_and_deferred_names_are_absent() -> None:
    """The deferred kinds (α as ``Eigen(time)``, ``Evolution``, ``Enclosed``,
    the never-minted ``Index(n)``) are not attributes of the module, and the
    struck ``CriticalParameter`` is defined nowhere under ``orpheus/`` (an AST
    census of class definitions; positive control: ``SpectralMap`` is found)."""
    present = [name for name in _DEFERRED + _STRUCK if hasattr(question_module, name)]
    require(not present, f"the module defines deferred or struck names: {present}")
    files = sorted((_ROOT / "orpheus").rglob("*.py"))
    require(len(files) > 100, f"activation: {len(files)} files under orpheus/")
    classes = {
        node.name
        for path in files
        for node in ast.walk(ast.parse(path.read_text()))
        if isinstance(node, ast.ClassDef)
    }
    require("SpectralMap" in classes, "positive control: the census did not find class SpectralMap")
    require(not (set(_STRUCK) & classes), f"a struck class is defined: {set(_STRUCK) & classes}")
    print(f"S7.2: {len(files)} files, {len(classes)} class names")


# ── S7.3: the fields are the roles ───────────────────────────────────────────

#: Keyed by NAME and resolved at run time (``getattr`` on the module), so the
#: row reads the class the module holds when it runs, not an object captured
#: at collection (a rebinding or a reload would otherwise leave it reading a
#: stale class: ``[M]`` 2026-10-02, battery arm A1 left a class-keyed row green).
_FIELDS = {
    "Eigen": ("parameter", "point", "mode"),
    "FixedSource": ("source", "point"),
    "Response": ("detector", "point"),
    "Fundamental": (),
    "Nearest": ("tau",),
}


@pytest.mark.rests_on(f"{_HERE}::test_s7_1_the_question_set_is_closed")
@pytest.mark.rests_on(f"{_MF}::test_s6_11_neither_type_carries_a_role_or_a_density")
@pytest.mark.parametrize("name", list(_FIELDS))
def test_s7_3_the_fields_are_exactly_the_roles(name: str) -> None:
    """``dataclasses.fields`` is exactly the ruled tuple: no forward/adjoint
    flag, no exclusive-or of a source and a detector, no σ (ruling 2 moved it
    into the point); and no value has an ``H``, ``adjoint`` or ``transpose``
    attribute (the eigen adjoint is part of the answer, posing ruling
    2026-09-28)."""
    cls = getattr(question_module, name)
    names = tuple(f.name for f in dataclasses.fields(cls))
    require(names == _FIELDS[name], f"{name} fields {names}, ruled {_FIELDS[name]}")
    hints = typing.get_type_hints(cls)
    flags = [n for n in names if hints.get(n) is bool]
    require(not flags, f"{cls.__qualname__} has boolean fields {flags}")
    dual = [a for a in ("H", "adjoint", "transpose", "dagger", "is_adjoint") if hasattr(cls, a)]
    require(not dual, f"{cls.__qualname__} carries {dual}")


_SPELLINGS = (
    ("fixed source given only a detector", lambda: FixedSource(detector=_table()), "detector"),  # type: ignore[call-arg]
    ("response given a source", lambda: Response(source=_table()), "source"),  # type: ignore[call-arg]
    ("response with no detector", lambda: Response(), "detector"),  # type: ignore[call-arg]
    ("fixed source with no source", lambda: FixedSource(), "source"),  # type: ignore[call-arg]
    ("eigen with no parameter", lambda: Eigen(), "parameter"),  # type: ignore[call-arg]
)


@pytest.mark.rests_on(f"{_HERE}::test_s7_3_the_fields_are_exactly_the_roles")
@pytest.mark.parametrize("build,fragment", [pytest.param(b, f, id=i.replace(" ", "_")) for i, b, f in _SPELLINGS])
def test_s7_3_a_role_spelled_in_the_wrong_value_does_not_construct(build, fragment: str) -> None:
    """The detector-only fixed source, the source-given response, the empty
    response and the parameter-less eigen are not values: Python's own
    argument binding refuses them, naming the argument. In particular numerics
    names no default parameter (k is the system's name, ruling 5)."""
    with pytest.raises(TypeError, match=fragment):
        build()


def test_s7_3_the_signatures() -> None:
    """``Response`` takes exactly one detector (and the point); the point and
    the mode have defaults, the datum and the parameter have none."""
    empty = inspect.Parameter.empty
    for cls, required, defaulted in (
        (Eigen, ("parameter",), ("point", "mode")),
        (FixedSource, ("source",), ("point",)),
        (Response, ("detector",), ("point",)),
    ):
        params = inspect.signature(cls).parameters
        require(tuple(params) == required + defaulted, f"{cls.__qualname__}{tuple(params)}")
        require(all(params[p].default is empty for p in required), f"{cls.__qualname__}: a datum has a default")
        # a default_factory field shows the ``<factory>`` sentinel, which is not ``empty``
        require(all(params[p].default is not empty for p in defaulted), f"{cls.__qualname__}: point/mode lack a default")
    with pytest.raises(TypeError, match="detector is a mesh-free function"):
        Response((_table(), _table()))  # type: ignore[arg-type]


# ── S7.4: the datum is a mesh-free function ──────────────────────────────────

_NOT_FUNCTIONS = (
    ("ndarray", lambda: np.array([[1.0, 0.5]]), "ndarray"),
    ("float", lambda: 1.0, "float"),
    ("None", lambda: None, "NoneType"),
    ("list", lambda: [[1.0, 0.5]], "list"),
)


@pytest.mark.rests_on(f"{_MF}::test_s6_8_the_table_is_a_read_only_copy")
@pytest.mark.parametrize("cls,role", [(FixedSource, "source"), (Response, "detector")], ids=["source", "detector"])
@pytest.mark.parametrize("datum,fragment", [pytest.param(b, f, id=i) for i, b, f in _NOT_FUNCTIONS])
def test_s7_4_a_datum_that_is_not_a_mesh_free_function_is_refused(cls, role: str, datum, fragment: str) -> None:
    """A bare ndarray is extensional and mesh-bound; the datum is a
    ``RegionwiseConstant`` or a ``Symbolic`` (step 6). The message names the
    ROLE and the received type."""
    with pytest.raises(TypeError, match=rf"the {role} is a mesh-free function.*{fragment}"):
        cls(datum())


@pytest.mark.parametrize("cls", [FixedSource, Response], ids=["source", "detector"])
@pytest.mark.parametrize("build", [_table, _symbolic], ids=["regionwise", "symbolic"])
def test_s7_4_both_function_types_are_admitted_in_both_roles(cls, build) -> None:
    """The positive leg (anti-#11): either step-6 type in either role."""
    value = cls(build())
    require(getattr(value, dataclasses.fields(cls)[0].name) == build(), "the datum was not kept")


# ── S7.5 and S7.6: the point and the parameter ──────────────────────────────

_BAD_OFFSETS = (
    ("nan", float("nan"), ValueError, "NaN, which is not a number"),
    ("inf", math.inf, ValueError, "infinite"),
    ("-inf", -math.inf, ValueError, "infinite"),
    ("complex", 1 + 0j, TypeError, "complex"),
    ("str", "0.1", TypeError, "str"),
)


@pytest.mark.parametrize("cls,build", [(Eigen, lambda p: Eigen(_KEY, point=p)),
                                        (FixedSource, lambda p: FixedSource(_table(), point=p)),
                                        (Response, lambda p: Response(_table(), point=p))],
                         ids=["eigen", "fixed_source", "response"])
@pytest.mark.parametrize("offset,error,fragment", [pytest.param(o, e, f, id=i) for i, o, e, f in _BAD_OFFSETS])
def test_s7_5_an_offset_that_is_not_a_finite_real_is_refused(cls, build, offset, error, fragment: str) -> None:
    """At CONSTRUCTION, the message naming the key (``'boron'``): NaN and an
    infinity are not points, a complex or a string offset is not real."""
    with pytest.raises(error, match=rf"'boron'.*{fragment}|{fragment}.*'boron'"):
        build({"boron": offset})


@dataclasses.dataclass(frozen=True, eq=False)
class _ByIdentity:
    """Hashable and contentless: a dataclass compared by identity."""


#: Keys a mapping can hold (hashable, so the refusal is production's, never
#: Python's own ``unhashable`` in the test body, ``vv-principles`` Mode 8(5)).
_BAD_KEYS = (
    ("a function", lambda: (lambda: 0), ContentlessError, "function"),
    ("an identity-compared object", _ByIdentity, ContentlessError, "identity"),
)


@pytest.mark.parametrize("key,error,fragment", [pytest.param(k, e, f, id=i.replace(" ", "_")) for i, k, e, f in _BAD_KEYS])
def test_s7_5_an_undigestable_point_key_is_refused_at_construction(key, error, fragment: str) -> None:
    """Eager: the refusal is raised by the constructor, never first by ``hash``."""
    with pytest.raises(error, match=fragment):
        Eigen(_KEY, point={key(): 0.5})


@pytest.mark.parametrize("key,error,fragment", [pytest.param(k, e, f, id=i.replace(" ", "_")) for i, k, e, f in
                                               _BAD_KEYS + (("a list", lambda: [1, 2], ContentlessError, "list"),)])
def test_s7_6_an_undigestable_parameter_is_refused_at_construction(key, error, fragment: str) -> None:
    """The parameter is an opaque, DIGESTABLE key (ruling 5): a function has
    no content; a list is a mutable part. The message carries the path
    ``Eigen.parameter``."""
    with pytest.raises(error, match=rf"Eigen\.parameter.*{fragment}"):
        Eigen(key())


@pytest.mark.parametrize("key", ["k", 7, ("cell", "nu-fission", 3), FrozenMapping({"cells": ("f",)})],
                         ids=["str", "int", "tuple", "frozen_mapping"])
def test_s7_6_a_digestable_parameter_is_admitted(key: Any) -> None:
    """The positive leg: numerics does not read the key, so any value with
    content is admitted; step 8 resolves it."""
    require(Eigen(key).parameter == key, "the key was not kept")


# ── S7.7: the mode ───────────────────────────────────────────────────────────


@pytest.mark.parametrize("tau,error,fragment", [
    pytest.param(float("nan"), ValueError, "NaN, which is not a number", id="nan"),
    pytest.param(math.inf, ValueError, "infinite", id="inf"),
    pytest.param(-math.inf, ValueError, "infinite", id="-inf"),
    pytest.param(0.5 + 0j, TypeError, "complex", id="complex"),
])
def test_s7_7_tau_is_a_finite_real(tau, error, fragment: str) -> None:
    with pytest.raises(error, match=rf"tau.*{fragment}"):
        Nearest(tau)


@pytest.mark.parametrize("mode", ["fundamental", 0.5, None], ids=["str", "bare_tau", "None"])
def test_s7_7_the_mode_is_a_mode_value(mode) -> None:
    """No stringly-typed selector, no bare τ standing for ``Nearest(τ)``."""
    with pytest.raises(TypeError, match="mode"):
        Eigen(_KEY, mode=mode)


# ── S7.8: the default point ──────────────────────────────────────────────────


@pytest.mark.rests_on(f"{_CI}::TestS54EncoderCanonicalForms")
@pytest.mark.parametrize("build", [lambda **kw: Eigen(_KEY, **kw), lambda **kw: FixedSource(_table(), **kw),
                                   lambda **kw: Response(_table(), **kw)], ids=["eigen", "fixed_source", "response"])
def test_s7_8_the_default_point_is_the_physical_point(build) -> None:
    """The default point is the empty mapping, an instance of THE
    ``orpheus.numerics.content.FrozenMapping`` (no second frozen-mapping type),
    and the explicit spellings (``{}``, ``FrozenMapping()``) are one value
    with it, digest-equal."""
    implicit = build()
    require(type(implicit.point) is FrozenMapping, f"the point is a {type(implicit.point)}")
    require(len(implicit.point) == 0, f"the default point is {implicit.point!r}")
    for explicit in (build(point={}), build(point=FrozenMapping())):
        require(explicit == implicit and hash(explicit) == hash(implicit), "explicit != default")
        require(explicit.content_digest == implicit.content_digest, "digests differ")


def test_s7_8_the_default_mode_is_fundamental() -> None:
    require(Eigen(_KEY).mode == Fundamental(), f"default mode {Eigen(_KEY).mode!r}")
    require(Eigen(_KEY) == Eigen(_KEY, mode=Fundamental()), "explicit Fundamental differs")


def test_s7_8_a_mapping_is_frozen_at_the_boundary() -> None:
    """A caller's dict is copied into a ``FrozenMapping``: writing to the dict
    afterwards moves neither the value nor its digest, and the stored point
    refuses item assignment."""
    offsets = {"boron": -0.25}
    q = Eigen(_KEY, point=offsets)
    before = q.content_digest
    offsets["boron"] = 0.75
    require(type(q.point) is FrozenMapping and q.point["boron"] == -0.25, f"the point aliased: {q.point!r}")
    require(q.content_digest == before, "the digest moved with the caller's dict")
    with pytest.raises(TypeError):
        q.point["boron"] = 1.0  # type: ignore[index]


def test_s7_8_a_zero_offset_is_not_the_empty_point() -> None:
    """``{key: 0.0}`` names the key, which step 8 must resolve, so it is a
    different question from the empty point: two cache keys, a miss and never
    a wrong hit (the orchestrator's ruling of 2026-10-02, point 3)."""
    require(Eigen(_KEY, point={"boron": 0.0}) != Eigen(_KEY), "a zero offset collapsed onto the physical point")


# ── S7.13: the layer ─────────────────────────────────────────────────────────

_LAYER_SCRIPT = textwrap.dedent(
    """
    import sys, importlib
    import numpy as np
    importlib.import_module(sys.argv[1])
    if sys.argv[1] == "orpheus.numerics.question":
        from orpheus.numerics.mesh_free_function import RegionwiseConstant
        from orpheus.numerics.question import Eigen, FixedSource, Nearest, Response
        t = RegionwiseConstant(np.ones((2, 2)))
        for q in (Eigen("k", {"b": 0.5}, Nearest(1.0)), FixedSource(t), Response(t)):
            q.content_digest
    import orpheus
    tops = sorted({k.split(".")[1] for k in sys.modules if k.startswith("orpheus.")})
    print(orpheus.__file__)
    print(" ".join(tops))
    print("sympy" in sys.modules)
    """
)


def _cold(module: str) -> tuple[str, set[str], bool]:
    env = {**os.environ, "PYTHONPATH": str(_ROOT)}
    out = subprocess.run([sys.executable, "-O", "-c", _LAYER_SCRIPT, module], capture_output=True, text=True,
                         env=env, cwd=_ROOT, timeout=300)
    require(out.returncode == 0, f"{module}: {out.stderr[-2000:]}")
    file, tops, sympy = out.stdout.splitlines()
    return file, set(tops.split()), sympy == "True"


def test_s7_13_the_module_loads_nothing_above_the_numerics_tier() -> None:
    """In a fresh interpreter, importing the module and building one value of
    each kind (a ``RegionwiseConstant`` datum) loads only ``orpheus.numerics``
    and the ``geometry`` the numerics package itself loads (``[M]``
    2026-10-02: a cold ``import orpheus.numerics`` loads exactly
    ``{geometry, numerics}``), and no SymPy. Positive control: a cold
    ``orpheus.sn.problem`` loads ``transport`` (the probe can see an upward
    load)."""
    file, tops, sympy = _cold("orpheus.numerics.question")
    require(file.startswith(str(_ROOT)), f"the subprocess imported {file} (L22)")
    require(tops <= {"numerics", "geometry"}, f"the question module loads {sorted(tops)}")
    require("numerics" in tops, "activation: numerics not loaded")
    require(not sympy, "building a question with a table loaded SymPy")
    _, control, _ = _cold("orpheus.sn.problem")
    require("transport" in control, f"positive control: a cold sn.problem loaded only {sorted(control)}")


def _orpheus_imports(source: str) -> list[str]:
    names: list[str] = []
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.ImportFrom) and node.module:
            names.append(node.module)
        elif isinstance(node, ast.Import):
            names.extend(alias.name for alias in node.names)
    return [n for n in names if n.startswith("orpheus")]


def test_s7_13_the_module_imports_only_numerics() -> None:
    """By AST, every ``orpheus`` import of the module (module-level, local and
    ``TYPE_CHECKING``) is under ``orpheus.numerics``: the question values are
    L1, below the input layer and the transport vocabulary (ruling 5). The
    activation leg: the known import of ``content`` is seen."""
    imported = _orpheus_imports(Path(question_module.__file__).read_text())
    require("orpheus.numerics.content" in imported, f"activation: the AST saw {imported}")
    above = [n for n in imported if not n.startswith("orpheus.numerics")]
    require(not above, f"the module imports above L1: {above}")
