r"""G04q: the metric arms compare by content (N04, spec 7.5, decision ``[ON04-6]``).

``DiagonalMetric``, ``DenseMetric``, ``DualMetric`` and ``FactoredMetric``: equal content gives ``==`` a ``bool`` (a ``numpy`` array or an exception is a red) and equal hashes
(each evaluated inside a ``try`` that fails the gate on any exception); different content differs, with different digests where the content differs; ``NoForm`` is already a
``bool`` and hashable (no first red). A ``DualMetric`` of equal primals is equal to another; the dual of a form differs from the form. ``FactoredMetric`` entries are
``(block shape, form)`` pairs with a ``NoForm`` block among them.

First red at the base tip ``cf924aaa`` (``[M]`` in the spec): ``DiagonalMetric ==`` returns an array, ``DenseMetric ==`` raises, ``DualMetric`` compares by identity, all unhashable.
Mutations (the spec's stand-in: ``ContentIdentity`` bases assigned in a probe): reverting ``DiagonalMetric``, ``DenseMetric`` or ``DualMetric`` each reddens its row; reverting
``FactoredMetric`` alone stays GREEN (its dataclass ``==`` compares ``entries``, whose arms are content-equal), so the replacement reverts ``DiagonalMetric`` and ``FactoredMetric`` together and reddens
the ``FactoredMetric`` row. Exact on any host: the content is small matrices the gate states.
"""

from __future__ import annotations

import numpy as np
import pytest

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _eq(a: object, b: object, what: str) -> bool:
    try:
        r = a == b
    except Exception as e:  # noqa: BLE001
        pytest.fail(f"{what}: == raised {type(e).__name__}: {e}")
    _require(type(r) is bool, f"{what}: == returned {type(r).__name__}, not bool")
    return r


def _hash(a: object, what: str) -> int:
    try:
        return hash(a)
    except Exception as e:  # noqa: BLE001
        pytest.fail(f"{what}: hash raised {type(e).__name__}: {e}")


def _forms() -> dict[str, tuple[object, object, object]]:
    """``{arm: (a, equal copy of a, different)}``."""
    from orpheus.numerics.declaration import Coord
    from orpheus.numerics.metric import DenseMetric, DiagonalMetric, DualMetric, Face, FactoredMetric, NoForm

    g = np.array([[2.0, 0.5], [0.5, 3.0]])
    k = Coord(2, ())
    diag = (DiagonalMetric(np.array([1.0, 2.0, 3.0])), DiagonalMetric(np.array([1.0, 2.0, 3.0])), DiagonalMetric(np.array([1.0, 2.0, 3.5])))
    dense = (DenseMetric(g, kernel=k, face=Face.METRIC), DenseMetric(g.copy(), kernel=k, face=Face.METRIC), DenseMetric(g, kernel=k, face=Face.GRAM))
    dual = (DualMetric(diag[0]), DualMetric(DiagonalMetric(np.array([1.0, 2.0, 3.0]))), DualMetric(diag[2]))
    fac = (
        FactoredMetric(entries=(((3,), diag[0]), ((2,), NoForm()))),
        FactoredMetric(entries=(((3,), DiagonalMetric(np.array([1.0, 2.0, 3.0]))), ((2,), NoForm()))),
        FactoredMetric(entries=(((3,), diag[2]), ((2,), NoForm()))),
    )
    return {"DiagonalMetric": diag, "DenseMetric": dense, "DualMetric": dual, "FactoredMetric": fac}


@pytest.mark.parametrize("arm", ["DiagonalMetric", "DenseMetric", "DualMetric", "FactoredMetric"])
def test_equal_content_is_equal_and_hashes_equal(arm: str) -> None:
    a, same, other = _forms()[arm]
    _require(_eq(a, same, arm), f"{arm}: equal content compares unequal")
    _require(_hash(a, arm) == _hash(same, arm), f"{arm}: equal content hashes unequal")
    _require(not _eq(a, other, arm), f"{arm}: different content compares equal")


def test_a_dual_differs_from_its_primal_and_no_form_is_a_bool() -> None:
    from orpheus.numerics.metric import NoForm

    a, same, _ = _forms()["DiagonalMetric"]
    from orpheus.numerics.metric import DualMetric

    _require(not _eq(DualMetric(a), a, "dual vs primal"), "the dual of a diagonal form equals the form")
    _require(_eq(NoForm(), NoForm(), "NoForm") and _hash(NoForm(), "NoForm") == _hash(NoForm(), "NoForm"), "NoForm is not a hashable value")
