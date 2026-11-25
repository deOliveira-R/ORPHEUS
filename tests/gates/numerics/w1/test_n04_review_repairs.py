r"""The repairs of the N04 review (workflow W1, unit N04, phase Amend): six rows, each with the input that reddens it.

(i) ``norm`` refuses a signed block of a composite space (``CoupledSpace``, ``FullFieldSpace``) with the issue #706 message, as it refuses a signed leaf.
Mutation: ``CoupledSpace.negative`` left at the base default 0 reddens the first row; ``FullFieldSpace.negative`` the second.
(ii) ``norm`` of a raw indefinite dense form under ``python -O`` raises (the inner product is negative) and never returns ``nan``. Mutation: the
``square < 0.0`` check removed reddens it under ``-O``.
(iii) a flip that lands on a space's own factors or on its dual returns the memoised object, a partial dual is the same space and class on every route; two routes to one genuinely
partial dual mint one member class. Mutation: ``_flip`` that rebuilds with the source's carrier reddens the route row.
(iv) the known defect that remains, declared by the ``CarrierSpec`` boundary: composed partial duals that reach ``V.dual()`` without passing
through ``V.dual()`` carry the default pair, so the member class follows the construction history; a strict xfail that W3 turns into a pass.
(v) a blank block of a ``FactoredMetric`` has one spelling: ``None`` is refused, and the metric of ``V.dual()`` is the dual of ``V``'s metric.
(vi) the allow-list is keyed on the module-qualified class name: a class of an admitted bare name defined elsewhere is refused.

Exact on any host: integer counts, identities and classes only.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.gates.numerics.w1 import _tensor_probe as T

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _signed_leaf():
    return T.tensor_space().declared(
        "s", (2,), _density(), form=_diag([1.0, -2.0])
    )


def _diag(w):
    from orpheus.numerics.metric import DiagonalMetric

    return DiagonalMetric(np.array(w, dtype=float))


def _density():
    from orpheus.numerics.factor import Kind, Variance, VarianceChart

    return Variance(Kind.FUNCTION, VarianceChart.DENSITY)


def _positive_leaf():
    return T.tensor_space().declared("p", (2,), _density(), form=_diag([1.0, 2.0]))


def test_norm_refuses_a_signed_block_of_a_coupled_space() -> None:
    from orpheus.numerics.coupled_system import CoupledSpace

    for members in ((_positive_leaf(), _signed_leaf()), (_signed_leaf(), _positive_leaf())):
        space = CoupledSpace.from_systems(members)
        _require(space.negative == 1, f"the coupled space counts {space.negative} negative directions, 1 expected")
        with pytest.raises(ValueError, match="#706"):
            space.norm(np.ones(4))
    _require(CoupledSpace.from_systems((_positive_leaf(), _positive_leaf())).negative == 0, "an all-positive coupled space counts a negative direction")


def test_norm_refuses_a_signed_block_of_a_full_field_space() -> None:
    from orpheus.numerics.spaces.full_field_space import FullFieldSpace

    space = FullFieldSpace.from_blocks(_positive_leaf(), _signed_leaf())
    _require(space.negative == 1, f"the composite counts {space.negative} negative directions, 1 expected")
    with pytest.raises(ValueError, match="#706"):
        space.norm(np.ones(4))


def test_norm_of_an_indefinite_dense_form_raises_in_every_build() -> None:
    from orpheus.numerics.declaration import Coord
    from orpheus.numerics.metric import DenseMetric, Face

    form = DenseMetric(np.array([[1.0, 2.0], [2.0, 1.0]]), kernel=Coord(2, ()), face=Face.GRAM)
    try:
        space = T.tensor_space().declared("d", (2,), _density(), form=form)
    except ValueError:
        return  # a build without -O refuses the form at its door: the formal path
    with pytest.raises(ValueError, match="#706"):
        space.norm(np.array([1.0, -1.0]))


def test_composed_partial_duals_land_on_the_memoised_object_or_agree() -> None:
    space = T.triple_space()
    first, second, third = space.labels
    twice = space.dual(on=[first]).dual(on=[first])
    column = np.ones(space.shape)
    _require(twice == space and type(twice.member(column)) is type(space.member(column)), "dual on one slot twice is not the space, with its member class")
    _require(space.riesz(on=()) is space and space.dual(on=()) is space, "an empty slot set is not the space itself")
    via_dual = space.dual().dual(on=[third])
    direct = space.dual(on=[first, second])
    _require(via_dual == direct, "two routes to one partial dual are unequal")
    _require(type(via_dual.member(column)) is type(direct.member(column)), "two routes to one partial dual mint two member classes")


@pytest.mark.xfail(strict=True, reason="the CarrierSpec boundary (issue #707, W3): the member class follows the construction history of a composed partial dual")
def test_member_class_follows_the_factors_after_composed_partial_duals() -> None:
    space = T.triple_space()
    first, second, third = space.labels
    composed = space.dual(on=[first]).dual(on=[second, third])
    _require(composed == space.dual(), "composed partial duals do not equal the all-slot dual")
    column = np.ones(space.shape)
    _require(type(composed.member(column)) is type(space.dual().member(column)), "the member class differs from that of V.dual()")


def test_a_blank_block_has_one_spelling() -> None:
    from orpheus.numerics.metric import FactoredMetric

    with pytest.raises(TypeError):
        FactoredMetric((((2,), None),))  # type: ignore[arg-type]
    space = T.triple_space()
    _require(space.dual()._resolved_metric == space._resolved_metric.dual(), "the metric of V.dual() is not the dual of V's metric")


def test_the_allow_list_is_keyed_on_the_module_qualified_name() -> None:
    from orpheus.numerics.space import TensorSpace

    namespace = {"__module__": "somewhere.else"}
    with pytest.raises(TypeError, match="allow-list"):
        type("TensorSpace", (TensorSpace,), namespace)
