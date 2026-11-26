r"""G04c: the subclasses of ``FunctionSpace`` are a pinned allow-list, refused at class creation (N04, spec 7.5, 6.5 "The allow-list").

A module-level mapping in ``numerics/space.py`` (the one whose keys include ``TensorSpace``; the spec names its shape, "qualified class name to the removing unit", and not its identifier)
(i) equals ``TensorSpace`` and the 9 legacy entries with their removing units (``CoupledSpace`` and ``FullFieldSpace`` N05; ``AngularTraceSpace``, ``AngularFaceTraceSpace``,
``ScalarTraceSpace``, ``_RadialCharacteristicSubSpace``, ``RadialCharacteristicInteriorSpace``, ``RadialCharacteristicBoundarySpace`` N06; ``DualSpace`` N07), and (ii) equals the runtime closure of
``FunctionSpace`` after importing every ``orpheus`` module outside ``derivations`` (the closure is read by walking ``__subclasses__`` and the import count is stated in the failure message);
(iii) a stranger subclass of ``FunctionSpace`` and a stranger subclass of an allow-listed class (a grandchild) each raise ``TypeError`` "arm is not on the FunctionSpace allow-list" at class creation;
(iv) ``ast`` over ``tests/``: no class defined there has ``FunctionSpace`` (or an allow-listed class) as a base (the two doubles ``_ToySpace`` and ``_ToyDualSpace`` of
``test_coupled_operator.py`` are re-spelled on ``TensorSpace.declared``).

First red at the base tip ``cf924aaa``: import-only for (i) to (iii) (no mapping: the closure of 10 classes, including ``TensorProductSpace``, exists today); (iv) reads the two doubles.
Mutations (described): (a) a legacy entry loses its unit reddens the unit leg; (b) an entry added to the mapping, with or without a unit, reddens the EQUALITY leg
(not the unit leg, so "without a unit" is replaced by (a)); (c) a hook that checks direct subclasses only reddens the grandchild leg. Exact on any host.
"""

from __future__ import annotations

import ast
import importlib
import pkgutil
from collections.abc import Mapping
from pathlib import Path

import pytest

pytestmark = pytest.mark.foundation

ROOT = Path(__file__).resolve().parents[4]
EXPECTED = {
    "DirectSumSpace": None,
    "AngularTraceSpace": "N06",
    "AngularFaceTraceSpace": "N06",
    "ScalarTraceSpace": "N06",
    "_RadialCharacteristicSubSpace": "N06",
    "RadialCharacteristicInteriorSpace": "N06",
    "RadialCharacteristicBoundarySpace": "N06",
    "DualSpace": "N07",
}


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _mapping() -> dict[str, object]:
    from orpheus.numerics import space

    found = [v for v in vars(space).values() if isinstance(v, Mapping) and any(str(k).rsplit(".", 1)[-1] == "TensorSpace" for k in v)]
    _require(len(found) == 1, f"numerics/space.py holds {len(found)} module-level mappings whose keys include TensorSpace (one expected)")
    return {str(k).rsplit(".", 1)[-1]: u for k, u in found[0].items()}


def _closure(root: type) -> list[type]:
    out: list[type] = []
    for s in root.__subclasses__():
        out.append(s)
        out += _closure(s)
    return out


def test_allow_list_is_the_pinned_set() -> None:
    mapping = _mapping()
    _require(set(mapping) == set(EXPECTED) | {"TensorSpace"}, f"the mapping's classes differ from the pinned set: {sorted(set(mapping) ^ (set(EXPECTED) | {'TensorSpace'}))}")
    for cls, unit in EXPECTED.items():
        if unit is None:
            _require(mapping[cls] is None, f"the arm {cls} that stays at N07 names a removing unit: {mapping[cls]!r}")
        else:
            _require(unit in str(mapping[cls]), f"the legacy entry {cls} does not name its removing unit {unit}: {mapping[cls]!r}")


def test_allow_list_equals_the_runtime_closure() -> None:
    import orpheus
    from orpheus.numerics.space import FunctionSpace

    imported = 0
    failed: list[str] = []
    for m in pkgutil.walk_packages(orpheus.__path__, "orpheus."):
        if m.name.startswith("orpheus.derivations"):
            continue
        imported += 1
        try:
            importlib.import_module(m.name)
        except Exception as e:  # noqa: BLE001  (counted and reported, never swallowed silently)
            failed.append(f"{m.name}: {type(e).__name__}")
    _require(not failed, f"{len(failed)} of {imported} orpheus modules did not import, so the closure is a lower bound: {failed[:4]}")
    closure = {c.__qualname__ for c in _closure(FunctionSpace) if c.__module__.startswith("orpheus")}
    mapping = set(_mapping())
    _require(closure == mapping, f"closure of {len(closure)} (over {imported} modules) differs from the mapping of {len(mapping)}: {sorted(closure ^ mapping)}")


def test_a_stranger_is_refused_at_any_depth() -> None:
    from orpheus.numerics.space import FunctionSpace, TensorSpace

    for base in (FunctionSpace, TensorSpace):
        with pytest.raises(TypeError, match="arm is not on the FunctionSpace allow-list"):
            type("Stranger", (base,), {})


def test_no_test_double_subclasses_a_space() -> None:
    names = {"FunctionSpace", "TensorSpace", "DualSpace", "DirectSumSpace"}

    def doubles(source: str) -> list[str]:
        return [
            n.name
            for n in ast.walk(ast.parse(source))
            if isinstance(n, ast.ClassDef) and any((isinstance(b, ast.Name) and b.id in names) or (isinstance(b, ast.Attribute) and b.attr in names) for b in n.bases)
        ]

    _require(doubles("class T(FunctionSpace): pass\nclass U(m.TensorSpace): pass") == ["T", "U"], "the census does not find planted subclasses of either spelling")
    found = {str(p.relative_to(ROOT)): d for p in sorted((ROOT / "tests").rglob("*.py")) if (d := doubles(p.read_text()))}
    _require(not found, f"test classes subclass a space: {found}")
