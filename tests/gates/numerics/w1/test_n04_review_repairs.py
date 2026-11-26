r"""The repairs of the N04 review (workflow W1, unit N04, phase Amend): six rows, each with the input that reddens it.

(i) ``norm`` refuses a signed block of a composite space (``CoupledSpace``, ``FullFieldSpace``) with the issue #706 message, as it refuses a signed leaf.
Mutation: ``CoupledSpace.negative`` left at the base default 0 reddens the first row; ``FullFieldSpace.negative`` the second.
(ii) an indefinite dense form has two refusals, one per build: without ``-O`` the constructor refuses it (a subprocess leg, whose ``orpheus`` is asserted to be this tree's);
under ``-O`` the constructor accepts it and ``norm`` raises on the negative inner product, naming issue #706 (an in-process leg). Mutation: the ``square < 0.0`` check removed reddens the
``-O`` leg; the constructor's sign test removed reddens the subprocess leg.
(iii) a flip that lands on a space's own factors or on its dual returns the memoised object, a partial dual is the same space and class on every route; two routes to one genuinely
partial dual mint one member class. Mutation: ``_flip`` that rebuilds with the source's carrier reddens the route row.
(iv) the known defect that remains, declared by the ``CarrierSpec`` boundary: the member class follows the construction path, and the space's own verbs reach two
equal spaces that mint different classes. Producer 1: ``V.dual(on=[a]).dual(on=[b, c])`` equals ``V.dual()`` and mints ``Member`` where ``V.dual()`` mints ``Covector``.
Producer 2: the dual of that space equals ``V`` and mints ``Covector`` where ``V`` mints ``Member``. Each is a strict xfail naming W3 and issue #707, so that W3 inherits a red when it
fixes the boundary; a green row states what holds today (``V.dual().dual() is V`` and the three member classes, ``on=()`` returns ``self`` for both verbs).
(v) a blank block of a ``FactoredMetric`` has one spelling: ``None`` is refused, and the metric of ``V.dual()`` is the dual of ``V``'s metric.
(vi) the allow-list is keyed on the module-qualified class name: a class of an admitted bare name defined elsewhere is refused.

Exact on any host: integer counts, identities and classes only.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

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


_INDEFINITE = r"""
import numpy as np
import orpheus
print("file: " + orpheus.__file__)
from orpheus.numerics.declaration import Coord
from orpheus.numerics.metric import DenseMetric, Face
try:
    DenseMetric(np.array([[1.0, 2.0], [2.0, 1.0]]), kernel=Coord(2, ()), face=Face.GRAM)
except ValueError as e:
    print("refused: " + str(e))
else:
    print("accepted")
"""


def test_a_debug_build_refuses_an_indefinite_dense_form_at_the_constructor() -> None:
    """Without ``-O`` the constructor's sign test (N01) refuses the matrix: read in a child whose ``orpheus`` is asserted to be this tree's."""
    root = Path(__file__).resolve().parents[4]
    env = {k: v for k, v in os.environ.items() if k != "PYTHONOPTIMIZE"}
    env["PYTHONPATH"] = os.pathsep.join([str(root), *([env["PYTHONPATH"]] if env.get("PYTHONPATH") else [])])
    r = subprocess.run([sys.executable, "-c", _INDEFINITE], env=env, cwd=root, capture_output=True, text=True, timeout=120)
    _require(r.returncode == 0, f"the child exited {r.returncode}: {r.stderr[-400:]}")
    lines = r.stdout.strip().splitlines()
    _require(bool(lines) and lines[0].startswith("file: " + str(root) + os.sep), f"the child did not import the tree under test ({root}): {lines[:1]}")
    _require(len(lines) == 2 and lines[1].startswith("refused: ") and "kept eigenvalue is not positive" in lines[1], f"the constructor did not refuse the indefinite form: {lines[1:]}")


@pytest.mark.skipif(__debug__, reason="without -O the constructor refuses the form; its leg is the subprocess row above")
def test_norm_of_an_indefinite_dense_form_raises_under_optimize() -> None:
    """Under ``python -O`` the constructor accepts the indefinite matrix and ``norm`` refuses the negative ``<x, x>`` naming issue #706."""
    from orpheus.numerics.declaration import Coord
    from orpheus.numerics.metric import DenseMetric, Face

    form = DenseMetric(np.array([[1.0, 2.0], [2.0, 1.0]]), kernel=Coord(2, ()), face=Face.GRAM)
    space = T.tensor_space().declared("d", (2,), _density(), form=form)
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


@pytest.mark.xfail(strict=True, reason="the CarrierSpec boundary (W3, issue #707): V.dual(on=[a]).dual(on=[b, c]) equals V.dual() and mints the class V mints, not the class V.dual() mints")
def test_member_class_follows_the_factors_after_composed_partial_duals() -> None:
    space = T.triple_space()
    first, second, third = space.labels
    composed = space.dual(on=[first]).dual(on=[second, third])
    _require(composed == space.dual(), "composed partial duals do not equal the all-slot dual")
    column = np.ones(space.shape)
    _require(type(composed.member(column)) is type(space.dual().member(column)), "the member class differs from that of V.dual()")


@pytest.mark.xfail(strict=True, reason="the CarrierSpec boundary (W3, issue #707): the dual of a space equal to V.dual() but built by composed partial duals equals V and mints Covector where V mints Member")
def test_member_class_of_the_dual_of_an_equal_copy_of_the_dual() -> None:
    space = T.triple_space()
    first, second, third = space.labels
    copy = space.dual(on=[first]).dual(on=[second, third])
    _require(copy == space.dual(), "the equal copy of V.dual() is not equal to V.dual()")
    back = copy.dual()
    _require(back == space, "the dual of an equal copy of V.dual() is not equal to V")
    column = np.ones(space.shape)
    _require(type(back.member(column)) is type(space.member(column)), "the member class of the dual of the copy differs from that of V")


def test_what_holds_today_of_the_carrier() -> None:
    from orpheus.numerics.space import Covector, Member

    space = T.triple_space()
    column = np.ones(space.shape)
    dual = space.dual()
    _require(dual.dual() is space, "V.dual().dual() is not V")
    _require(type(space.member(column)) is Member, "V does not mint Member")
    _require(type(dual.member(column)) is Covector, "V.dual() does not mint Covector")
    _require(type(dual.dual().member(column)) is Member, "V.dual().dual() does not mint Member")
    _require(space.dual(on=()) is space and space.riesz(on=()) is space, "on=() does not return the space itself for both verbs")


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
