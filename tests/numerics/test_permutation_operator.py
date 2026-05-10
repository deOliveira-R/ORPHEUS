r"""Tests for :class:`orpheus.numerics.operator.PermutationOperator`.

Verifies the operator's invariants:

* ``apply`` gathers entries along the tagged axis per the permutation.
* ``apply_transpose`` applies the inverse permutation (so round-trip is identity).
* Involution detection: ``is_involution`` matches
  ``np.array_equal(perm[perm], np.arange(n))``.
* Capability set: ``{CAP_APPLY, CAP_APPLY_TRANSPOSE}`` — solve is NOT advertised.
* Composition with ``@`` reproduces function composition.
* Axis broadcasting on non-zero axis.
* Construction-time validation of malformed permutations.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    IdentityOperator,
    PermutationOperator,
)


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
def test_involution_detection():
    """L1: ``is_involution`` flag matches ``perm[perm] == arange(n)``.

    Pair-swap permutations are involutions; cyclic permutations of
    length > 2 are not. The flag is what downstream consumers
    (e.g. SN specular reflection) use to know that
    ``apply_transpose == apply``.
    """
    pair_swap = np.array([1, 0, 3, 2])
    P_inv = PermutationOperator(pair_swap)
    assert P_inv.is_involution is True
    assert np.array_equal(pair_swap[pair_swap], np.arange(4))
    # And ``apply_transpose == apply`` for the involution.
    x = np.array([1.0, 2.0, 3.0, 4.0])
    np.testing.assert_array_equal(P_inv.apply(x), P_inv.apply_transpose(x))

    cyclic = np.array([2, 0, 1])
    P_cyc = PermutationOperator(cyclic)
    assert P_cyc.is_involution is False
    assert not np.array_equal(cyclic[cyclic], np.arange(3))


@pytest.mark.l0
def test_capability_set():
    """L0: capabilities advertise ``{CAP_APPLY, CAP_APPLY_TRANSPOSE}``."""
    P = PermutationOperator(np.array([1, 0, 2]))
    assert P.capabilities == frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})


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
def test_compose_with_identity():
    """L0: ``(P @ Identity).apply(x) == P.apply(x)``."""
    P = PermutationOperator(np.array([2, 0, 1]))
    I = IdentityOperator()
    x = np.array([10, 20, 30])
    np.testing.assert_array_equal((P @ I).apply(x), P.apply(x))
    np.testing.assert_array_equal((I @ P).apply(x), P.apply(x))
