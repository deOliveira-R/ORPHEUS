r"""G04a: one tensor type; the product is concatenation of slot lists (N04, spec 7.5; law L1 of 2.2).

``of_axes(a, e, s)``, ``of_axes(a, e) * of_axes(s)`` and ``of_axes(a) * of_axes(e, s)`` are ONE value of ONE type (both bracketings), over the four slot pairs of
a slab-like triple and a dense head, on a positive triple and on a signed one (a signed diagonal form is a legal factor from N03b). The unit leg
``of_axes() * X == X`` is asserted when the empty product constructs (``[R]`` in the spec: what ``1`` is in ``TensorSpace``); if the empty product is refused
the refusal must be a ``ValueError`` (L1 then states the unit as refused). A second leg parses ``orpheus/``: ``TensorProductSpace`` is defined nowhere.

First red at the base tip ``cf924aaa``: import-only for the space legs (``TensorSpace`` is not in ``orpheus.numerics.space``); the ``ast`` leg is
behavioural (``TensorProductSpace`` is defined once in ``numerics/space.py``). Mutation: swap the operands in ``TensorSpace.__mul__`` reddens the ``==`` leg
(the type leg has no spellable red once the allow-list hook of G04c stands: any other subclass is refused at class creation). Exact on any host: identity
relations and a parse.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from tests.gates.numerics.w1 import _tensor_probe as P

pytestmark = pytest.mark.foundation

ROOT = Path(__file__).resolve().parents[4]


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _triples() -> dict[str, tuple]:
    from orpheus.numerics.space import TensorSpace  # noqa: F401  (import-only red at the base tip)

    return {"positive": P.triple_factors(), "signed": P.triple_factors(angular_w=(1.3, -0.7, 2.1, 0.9))}


def test_both_bracketings_are_one_value_of_one_type() -> None:
    TS = P.tensor_space()
    for name, (a, e, s) in _triples().items():
        whole = TS.of_axes(a, e, s)
        left = TS.of_axes(a, e) * TS.of_axes(s)
        right = TS.of_axes(a) * TS.of_axes(e, s)
        _require(type(whole) is type(left) is type(right), f"[{name}] the product has two types: {type(whole).__name__}, {type(left).__name__}, {type(right).__name__}")
        _require(whole == left, f"[{name}] of_axes(a, e, s) != of_axes(a, e) * of_axes(s)")
        _require(whole == right, f"[{name}] of_axes(a, e, s) != of_axes(a) * of_axes(e, s)")
        _require(hash(whole) == hash(left) == hash(right), f"[{name}] equal products hash unequal")
        _require(whole.factors == (a, e, s), f"[{name}] the slot list is not the concatenation (slot order lost)")


def test_four_slot_pairs_and_a_dense_head() -> None:
    from tests.gates.numerics.w1 import _factor_probe as FP

    TS = P.tensor_space()
    a, e, s = P.triple_factors()
    axis, kernel, _g = FP.dense_axis(4, 3, True)
    from orpheus.numerics import factor as F
    from orpheus.numerics.declaration import Coord

    head = F.Factor(axis, F.Variance(F.Kind.MEASURE, F.VarianceChart.INTEGRATED), F.Subquotient(kernel.ann(), Coord(4, ())))
    for x, y in ((a, e), (a, s), (e, s), (head, s)):
        _require(TS.of_axes(x) * TS.of_axes(y) == TS.of_axes(x, y), f"{x.axis.label} * {y.axis.label}: the product is not of_axes of the concatenation")
        _require(TS.of_axes(x) * TS.of_axes(y) != TS.of_axes(y, x), f"{x.axis.label} * {y.axis.label} equals its swap: a slot order is not part of identity")


def test_unit_leg() -> None:
    TS = P.tensor_space()
    a, e, s = P.triple_factors()
    X = TS.of_axes(a, e, s)
    try:
        unit = TS.of_axes()
    except ValueError:
        return  # the empty product is refused: L1 states the unit as refused
    _require(unit * X == X, "1 (x) X != X")
    _require(X * unit == X, "X (x) 1 != X")
    _require(unit.factors == (), "the unit has factors")


def _defined_class_names(root: Path) -> list[str]:
    hits: list[str] = []
    for path in sorted((root / "orpheus").rglob("*.py")):
        tree = ast.parse(path.read_text())
        hits += [f"{path.relative_to(root)}:{n.lineno}" for n in ast.walk(tree) if isinstance(n, ast.ClassDef) and n.name == "TensorProductSpace"]
    return hits


def test_tensor_product_space_is_defined_nowhere() -> None:
    # positive control of the instrument: a planted definition is found by the same walk
    planted = [n for n in ast.walk(ast.parse("class TensorProductSpace: pass")) if isinstance(n, ast.ClassDef) and n.name == "TensorProductSpace"]
    _require(len(planted) == 1, "the census predicate does not find a planted class")
    found = _defined_class_names(ROOT)
    _require(not found, f"TensorProductSpace is still defined: {found}")
