r"""Tests for :class:`orpheus.numerics.operator.PermutationOperator`.

Verifies the operator's invariants:

* ``apply`` gathers entries along the tagged axis per the permutation.
* ``apply_transpose`` applies the inverse permutation (so round-trip is identity).
* ``inverse()`` returns the inverse permutation as a first-class
  :class:`PermutationOperator` (since :math:`P^{-1} = P^{T}` for
  permutation matrices) — carve P4: ``solve`` retired, the inverse is
  ALGEBRA-CLOSED, solving is ``.inverse().apply(b)``. Issue #150 closeout.
* Involution, **asked of the algebra**: ``(P @ P).apply(x) == x``. The
  ``is_involution`` attribute this file used to assert against was retired
  in G6.3 step 5 (#330) — see :ref:`bc-narrowed-involution`. The
  behavioural claim was rewired here, not deleted with the flag.
* Predicates: ``is_invertible`` and ``is_adjointable`` are both ``True``
  (a permutation always transposes and always inverts).
* Composition with ``@`` reproduces function composition.
* Axis broadcasting on non-zero axis.
* Construction-time validation of malformed permutations.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.operator import (
    IdentityOperator,
    PermutationOperator,
)
from orpheus.numerics.space import FunctionSpace


@pytest.mark.l0
def test_apply_correctness_axis0():
    """L0: known perm [2, 0, 1] on [10, 20, 30] returns [30, 10, 20]."""
    P = PermutationOperator(np.array([2, 0, 1]))
    x = np.array([10, 20, 30])
    np.testing.assert_array_equal(P.apply(x), np.array([30, 10, 20]))


@pytest.mark.l0
def test_apply_transpose_round_trip():
    """L0: ``apply_transpose(apply(x)) == x`` via the inverse permutation."""
    perm = np.array([3, 1, 4, 0, 2])
    P = PermutationOperator(perm)
    x = np.array([7.0, -1.5, 11.0, 0.25, 4.5])
    np.testing.assert_array_equal(P.apply_transpose(P.apply(x)), x)
    # And the other direction: apply(apply_transpose(x)) == x.
    np.testing.assert_array_equal(P.apply(P.apply_transpose(x)), x)


@pytest.mark.l1
def test_the_square_of_a_pair_swap_is_the_identity():
    r"""L1: :math:`P \circ P = \mathrm{id}`, **asked of the algebra**.

    This is the rewire of the retired ``is_involution`` flag (G6.3 step 5,
    #330). The flag stored ``perm[perm] == arange(n)``; the composition
    ``P @ P`` *computes* it, and — the reason for the retirement — it also
    REFUSES to exist when the operator is bound between two different
    spaces, where a stored ``bool`` was obliged to answer anyway. The
    cross-space refusal is gated in
    ``tests/sn/operators/test_specular_deck_chain.py``; this file keeps the
    positive, unbound claim the flag used to carry.

    The negative leg is the same question asked of a 3-cycle, which is a
    permutation of the same shape and NOT an involution — so a mutation
    that made ``@`` return the identity, or made ``apply`` a no-op, cannot
    pass both legs.
    """
    x = np.array([1.0, 2.0, 3.0, 4.0])

    pair_swap = PermutationOperator(np.array([1, 0, 3, 2]))
    np.testing.assert_array_equal((pair_swap @ pair_swap).apply(x), x)
    # An involution is its own inverse, so the transpose gathers the same
    # rows as the forward — the property the retired flag advertised.
    np.testing.assert_array_equal(
        pair_swap.apply(x), pair_swap.apply_transpose(x)
    )

    three_cycle = PermutationOperator(np.array([2, 0, 1]))
    y = np.array([1.0, 2.0, 3.0])
    assert not np.array_equal((three_cycle @ three_cycle).apply(y), y), (
        "a 3-cycle squared is the OTHER 3-cycle, not the identity — if this "
        "passes, `@` is not composing"
    )
    # Its cube is, which pins that the composition is genuine rather than
    # merely non-identity.
    np.testing.assert_array_equal(
        (three_cycle @ three_cycle @ three_cycle).apply(y), y
    )


@pytest.mark.l0
def test_invertible_adjointable_and_solve_retired():
    """L0: both predicates ``True``; the ``solve`` verb is retired.

    Issue #150 closeout, restated post carve P4: a permutation matrix
    is always invertible (:math:`P^{-1} = P^{T}`), so ``is_invertible``
    and ``is_adjointable`` are both ``True``. The ``solve`` verb was
    retired at carve P4 (algebra-closed inverse) — the retirement pin
    keeps a resurrected twin path from shipping silently.
    """
    P = PermutationOperator(np.array([1, 0, 2]))
    assert P.is_invertible
    assert P.is_adjointable
    assert not hasattr(P, "solve")  # retired at carve P4


@pytest.mark.l1
def test_inverse_apply_is_inverse_of_apply():
    r"""L1: ``inverse().apply(apply(x)) == x`` for any input.

    For a permutation matrix, :math:`P^{-1} = P^{T}` — both the
    transpose and the inverse are the same operation. :meth:`inverse`
    returns the inverse permutation as a first-class
    :class:`PermutationOperator` whose :meth:`apply` is the same
    ``np.take(·, inverse_perm)`` gather as :meth:`apply_transpose`.

    Rewired at carve P4: ``PermutationOperator.solve`` retired
    (algebra-closed inverse); solving is ``.inverse().apply`` —
    same values, same (exact) tolerance: the gather is integer
    indexing, bit-identical.
    """
    perm = np.array([3, 1, 4, 0, 2])
    P = PermutationOperator(perm)
    x = np.array([7.0, -1.5, 11.0, 0.25, 4.5])
    # inverse().apply(apply(x)) == x — inverse round-trip.
    np.testing.assert_array_equal(P.inverse().apply(P.apply(x)), x)
    # apply(inverse().apply(x)) == x — the other direction.
    np.testing.assert_array_equal(P.apply(P.inverse().apply(x)), x)
    # inverse().apply and apply_transpose agree (both are the inverse
    # for a permutation matrix).
    np.testing.assert_array_equal(P.inverse().apply(x), P.apply_transpose(x))


@pytest.mark.l0
def test_composition_is_function_composition():
    """L0: ``(P1 @ P2).apply(x)`` applies P2 then P1.

    Operator product is function composition: :math:`(P_1 P_2) x =
    P_1(P_2 x)`. Verify against an explicit two-step application.
    """
    perm1 = np.array([2, 0, 1])  # cycles (10, 20, 30) -> (30, 10, 20)
    perm2 = np.array([1, 2, 0])  # cycles (10, 20, 30) -> (20, 30, 10)
    P1 = PermutationOperator(perm1)
    P2 = PermutationOperator(perm2)
    x = np.array([10, 20, 30])
    composed = (P1 @ P2).apply(x)
    expected = P1.apply(P2.apply(x))
    np.testing.assert_array_equal(composed, expected)


@pytest.mark.l0
def test_axis_nonzero_apply_2d():
    """L0: shape ``(3, 4)``, axis=1, perm=[3, 0, 2, 1] reorders each row."""
    perm = np.array([3, 0, 2, 1])
    P = PermutationOperator(perm, axis=1)
    x = np.array(
        [
            [1, 2, 3, 4],
            [5, 6, 7, 8],
            [9, 10, 11, 12],
        ]
    )
    out = P.apply(x)
    expected = np.array(
        [
            [4, 1, 3, 2],
            [8, 5, 7, 6],
            [12, 9, 11, 10],
        ]
    )
    np.testing.assert_array_equal(out, expected)


@pytest.mark.l0
def test_invalid_perm_duplicate_raises():
    """L0: non-permutation (duplicate entry) is rejected at construction."""
    with pytest.raises(ValueError, match="permutation"):
        PermutationOperator(np.array([0, 0, 1]))


@pytest.mark.l0
def test_invalid_perm_out_of_range_raises():
    """L0: out-of-range entries are rejected at construction."""
    with pytest.raises(ValueError, match="permutation"):
        PermutationOperator(np.array([0, 1, 3]))


@pytest.mark.l0
def test_invalid_perm_2d_raises():
    """L0: non-1-D perm input is rejected."""
    with pytest.raises(ValueError, match="1-D"):
        PermutationOperator(np.array([[0, 1], [1, 0]]))


@pytest.mark.l0
@pytest.mark.parametrize("end", ["domain", "codomain"])
def test_a_bound_space_whose_extent_contradicts_the_perm_is_refused(end):
    """L0: the space and ``n`` state the same fact; a disagreement raises.

    Construction is the only place this can be caught — a mis-bound space is
    SILENT at apply-time, because the arrays still broadcast and the operator
    computes a plausible wrong answer.

    ⚠ The extent check is structurally unable to catch a domain/codomain
    **swap**: on every shipped SN quadrature ``|Γ₊| == |Γ₋|``, so both ends
    have the same extent. That claim needs an identity assertion instead, and
    lives in ``tests/sn/operators/test_specular_deck_chain.py``.
    """
    perm = np.array([1, 0, 3, 2])
    wrong = FunctionSpace(name="too-long", shape=(5,))
    with pytest.raises(ValueError, match="along axis 0"):
        PermutationOperator(perm, axis=0, **{end: wrong})

    # Positive leg: the matching extent constructs, so the row above is
    # about the DISAGREEMENT and not about binding being rejected at all.
    right = FunctionSpace(name="matching", shape=(4,))
    bound = PermutationOperator(perm, axis=0, **{end: right})
    assert getattr(bound, end) is right


@pytest.mark.l0
def test_compose_with_identity():
    """L0: ``(P @ Identity).apply(x) == P.apply(x)``."""
    P = PermutationOperator(np.array([2, 0, 1]))
    I = IdentityOperator()
    x = np.array([10, 20, 30])
    np.testing.assert_array_equal((P @ I).apply(x), P.apply(x))
    np.testing.assert_array_equal((I @ P).apply(x), P.apply(x))
