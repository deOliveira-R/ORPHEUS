r"""The adjointability law on the numerics operators: leg 1 over :data:`NUMERICS_REGISTRY`, and the inverse family.

The law, its two legs and who runs which are stated once, in :mod:`tests.gates.numerics._adjointability`. The
consumers' census is ``tests/gates/test_adjointability_is_honest.py``.
"""

from __future__ import annotations

from dataclasses import field

import numpy as np
import pytest

from orpheus.numerics.green_operator import GreenOperator
from orpheus.numerics.elementary_operators import DiagonalOperator
from orpheus.numerics.operator import IdentityOperator, LinearOperator, MissingAdjoint, posed
from tests.gates.numerics._adjointability import NUMERICS_REGISTRY, assert_doors, flat

pytestmark = [
    pytest.mark.foundation,
    pytest.mark.rests_on(
        "tests/gates/numerics/test_operator.py::test_adjoint_weighted_identity",
        "tests/gates/numerics/test_riesz_legs.py::"
        "test_lowering_realizes_the_metric_pairing",
        "tests/gates/numerics/test_riesz_legs.py::"
        "test_the_adjoint_composite_routes_through_the_legs",
        "tests/gates/numerics/test_operator_capability_predicates.py::"
        "test_inverse_adjoint_contract_keystone_v2",
    ),
]


@pytest.mark.parametrize("row", list(NUMERICS_REGISTRY))
@pytest.mark.catches("ERR-117")
def test_each_transpose_door_opens_iff_predicted(row: str) -> None:
    r"""``.dual()`` works ⟺ ``is_adjointable``; ``.H`` works ⟺ ``is_adjointable ∧ posed``.

    A refused door raises ``MissingAdjoint`` at construction. Where a door
    opens, the sanity leg: Euclidean reciprocity for ``.dual()``, and
    ``⟨A x, y⟩_W = ⟨x, A.H y⟩_V`` in the operator's own end metrics for ``.H``.
    """
    assert_doors(row, NUMERICS_REGISTRY[row]())


class _SolveOnlyLeaf(LinearOperator):
    r"""An invertible leaf with no transpose: its inverse has no transposed solve."""

    def __init__(self, c: np.ndarray) -> None:
        self._c = np.asarray(c, dtype=float)

    @property
    def domain(self):
        return None

    @property
    def codomain(self):
        return None

    def apply(self, x, /):
        return self._c * x

    def solve(self, b, /):
        return b / self._c

    @property
    def is_invertible(self) -> bool:
        return True


def _inverse_family_rows():
    from orpheus.numerics.operator import InverseOperator

    D = DiagonalOperator(np.array([2.0, -0.5, 4.0]), shape=(3,), spans=(0,))
    return [
        ("inverse_of_diagonal", D.inverse(), np.array([1.0, 3.0, -2.0])),
        ("inverse_of_solve_only_leaf", InverseOperator(_SolveOnlyLeaf(np.array([2.0, 3.0]))), np.ones(2)),
    ]


#: The SN sweep's row is the root census's (layer DAG Phase 2 refresh):
#: tests/gates/test_adjointability_is_honest.py::test_the_inverse_swaps_the_forwards_transpose_verbs_on_the_sn_sweep.
@pytest.mark.parametrize("index", range(2), ids=["inverse_of_diagonal", "inverse_of_solve_only_leaf"])
def test_the_inverse_swaps_the_forwards_transpose_verbs(index: int) -> None:
    r"""``(A⁻¹)ᵀ`` solves with ``Aᵀ``: the inverse's transposed solve IS the
    forward's transpose, and it exists iff the forward transposes."""
    from orpheus.numerics.operator import adjointable

    name, inverse, y = _inverse_family_rows()[index]
    expected = adjointable(inverse.inner)
    if inverse.is_transpose_invertible != expected:
        pytest.fail(
            f"[{name}] is_transpose_invertible={inverse.is_transpose_invertible}, "
            f"but the forward's adjointable is {expected}"
        )
    if not expected:
        return
    got, ref = flat(inverse.solve_transpose(y)), flat(inverse.inner.apply_transpose(y))
    if not np.array_equal(got, ref):
        pytest.fail(f"[{name}] solve_transpose ≠ the forward's apply_transpose (max |Δ| {np.max(np.abs(got - ref)):.2e})")


def test_the_inverse_of_a_multiplier_inherits_its_adjoint_posing() -> None:
    r"""The inverse inherits the forward's adjoint posing. A constant
    multiplier is metric-free, so ``(I⁻¹)ᴴ`` is posed with no spaces; a
    spaceless FIELD multiplier is not (S6a review, #625), so ``(D⁻¹)ᴴ``
    refuses as unbound while its transpose ``(D⁻¹).dual()`` is ``diag(1/d)``
    exactly. Until the review this row asserted ``(D⁻¹)ᴴ = diag(1/d)`` on
    unbound ends, the metric-free reading of a field multiplier."""
    from orpheus.numerics.operator import IdentityOperator

    d = np.array([2.0, -0.5, 4.0])
    y = np.array([1.0, 3.0, -2.0])
    if not np.array_equal(IdentityOperator().inverse().H.apply(y), y):
        pytest.fail("(I⁻¹)ᴴ is not the identity")
    inv = DiagonalOperator(d, shape=(3,), spans=(0,)).inverse()
    with pytest.raises(MissingAdjoint, match="UNBOUND"):
        _ = inv.H
    got = inv.dual().apply(y)
    if not np.array_equal(got, y / d):
        pytest.fail(f"(D⁻¹)ᵀ y = {got}, expected diag(1/d) y = {y / d}")


def test_a_converged_first_partial_sum_is_returned_unswept() -> None:
    r"""A cold Green start is the first Neumann partial sum, checked before it
    is driven: when it already meets the promise, one preconditioner solve is
    the whole apply, as it was when the driver started from zero (`[M]`
    2026-10-18: 1 solve old, 2 without the check)."""
    from orpheus.numerics.elementary_operators import PermutationOperator
    calls = [0]

    class _CountingDiagonal(DiagonalOperator):
        def solve(self, b_vec):
            calls[0] += 1
            return super().solve(b_vec)

    D = _CountingDiagonal(np.array([2.0, 3.0, 4.0, 5.0]), shape=(4,), spans=(0,))
    gain = 3e-12 * PermutationOperator(np.array([1, 2, 3, 0]))
    green = GreenOperator(D, gain, max_iter=10, tol=1e-10)
    green.apply(np.array([1.0, -2.0, 0.5, 3.0]))
    if calls[0] != 1:
        pytest.fail(f"a converged first partial sum cost {calls[0]} preconditioner solves, expected 1")

