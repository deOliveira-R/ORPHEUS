r"""G04s, G04t, G04r: the constructor, ``with_axis`` and the one sub-product ``sub`` (N04, spec 7.5; decision ``[ON04-3]``).

G04s. ``inspect.signature(TensorSpace.__init__)`` has the parameters ``[self, factors, carrier]`` (the keyword-only carrier); ``name=``, ``shape=``, ``metric=``,
``inner_product_weights=`` each raise ``TypeError``. G04t. ``X.with_axis(f)`` appends last and equals ``X * of_axes(f)``; it refuses an ``Axis`` (``TypeError``);
``of_axes(a, a)`` and ``X.with_axis(a)`` with ``a`` already a slot raise ``ValueError`` (labels are unique); ``of_axes(Axis)`` is refused (``TypeError``, a variance is stated).
G04r. ``sub(*labels)`` is the sub-product in slot order and refuses an unknown or an out-of-order label; ``tail(k)`` and ``without_axis(label)`` equal their ``sub``; the two
production readers of ``space.factors[k]`` (``transport/fields/_bases.py``, ``harmonic_moment_flux.py``) read ``factor(label)`` and ``tail(1)`` instead (``ast``: 0 subscripts of
``<x>.space.factors`` in those two files; the opening counted 2, plus 3 more in docstrings that an ``ast`` parse does not read).

First red at the base tip ``cf924aaa``: import-only for the behavioural legs; the ``ast`` leg reads 2 subscripts. Mutations (described): ``shape`` or ``metric`` re-declared ``init=True``
reddens G04s (a fixture subclass re-declaring the four fields ``init=False``, the 3.4 spelling, is green); ``with_axis`` that prepends reddens G04t; wrapping an ``Axis``
silently reddens; removing the label-uniqueness invariant reddens; ``sub`` returning the labels reversed reddens G04r; removing the out-of-order refusal reddens. Exact on any host.
"""

from __future__ import annotations

import ast
import inspect
from pathlib import Path

import pytest

from tests.gates.numerics.w1 import _tensor_probe as P

pytestmark = pytest.mark.foundation

ROOT = Path(__file__).resolve().parents[4]


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def test_init_signature() -> None:
    TS = P.tensor_space()
    names = list(inspect.signature(TS.__init__).parameters)
    _require(names == ["self", "factors", "carrier"], f"TensorSpace.__init__ takes {names}")
    a, e, s = P.triple_factors()
    for bad in ("name", "shape", "metric", "inner_product_weights"):
        with pytest.raises(TypeError):
            TS(factors=(a, e, s), **{bad: None})


def test_with_axis_and_unique_labels() -> None:
    TS = P.tensor_space()
    a, e, s = P.triple_factors()
    X = TS.of_axes(a, e)
    _require(X.with_axis(s) == TS.of_axes(a, e, s), "with_axis does not append last")
    _require(X.with_axis(s) == X * TS.of_axes(s), "with_axis != X * of_axes(f)")
    with pytest.raises(TypeError):
        X.with_axis(s.axis)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        TS.of_axes(a.axis)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        TS.of_axes(a, a)
    with pytest.raises(ValueError):
        X.with_axis(a)


def test_sub_is_the_one_sub_product() -> None:
    TS = P.tensor_space()
    a, e, s = P.triple_factors()
    V = TS.of_axes(a, e, s)
    _require(V.sub("angular", "spatial") == TS.of_axes(a, s), "sub is not the sub-product in slot order")
    _require(V.sub("energy") == TS.of_axes(e), "sub of one label")
    _require(V.sub("angular", "energy", "spatial") == V, "sub of every label is not the space")
    with pytest.raises((KeyError, ValueError)):
        V.sub("nope")
    with pytest.raises((KeyError, ValueError)):
        V.sub("spatial", "angular")
    _require(V.tail(2) == V.sub("energy", "spatial"), "tail(2) is not its sub")
    _require(V.tail(1) == V.sub("spatial"), "tail(1) is not its sub")
    _require(V.without_axis("energy") == V.sub("angular", "spatial"), "without_axis is not its sub")


def _factor_subscripts(source: str) -> int:
    n = 0
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Subscript) and isinstance(node.value, ast.Attribute) and node.value.attr == "factors":
            owner = node.value.value
            if isinstance(owner, ast.Attribute) and owner.attr == "space":
                n += 1
    return n


def test_the_production_readers_read_by_label() -> None:
    _require(_factor_subscripts("x = self.space.factors[1]\ny = space.factors[0]") == 1, "the census does not find a planted self.space.factors[k]")
    found = {}
    for rel in ("orpheus/transport/fields/_bases.py", "orpheus/transport/fields/harmonic_moment_flux.py"):
        found[rel] = _factor_subscripts((ROOT / rel).read_text())
    _require(sum(found.values()) == 0, f"readers of space.factors[k] remain: {found}")
