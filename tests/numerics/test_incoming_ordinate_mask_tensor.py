r"""Tests for :class:`orpheus.numerics.operator.IncomingOrdinateMaskTensor`.

Verifies the operator's invariants:

* ``apply`` zeroes only the indices in ``inflow_indices``; preserves the rest.
* Self-adjointness: ``apply == apply_transpose``.
* Idempotence: ``M @ M == M``.
* Empty ``inflow_indices`` is identity on the masked axis.
* Full ``inflow_indices == arange(n_ordinates)`` matches
  :class:`ZeroOperator` action.
* Composition with :class:`IdentityOperator`.
* Axis broadcasting on a non-zero axis.
* Construction-time validation (out-of-range, duplicates).
* The input is NOT mutated by :meth:`apply`.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    IdentityOperator,
    IncomingOrdinateMaskTensor,
    ZeroOperator,
)


@pytest.mark.l0
def test_apply_zeroes_only_inflow_indices():
    """L0: indices [0, 2] on [1, 2, 3, 4] zero positions 0 and 2."""
    M = IncomingOrdinateMaskTensor(np.array([0, 2]), n_ordinates=4)
    x = np.array([1.0, 2.0, 3.0, 4.0])
    np.testing.assert_array_equal(M.apply(x), np.array([0.0, 2.0, 0.0, 4.0]))


@pytest.mark.l0
def test_leaves_non_inflow_untouched():
    """L0: non-inflow positions match the input exactly under random data."""
    rng = np.random.default_rng(seed=3)
    n = 8
    x = rng.standard_normal(n)
    inflow = np.array([1, 4, 6])
    M = IncomingOrdinateMaskTensor(inflow, n_ordinates=n)
    out = M.apply(x)
    # Outflow positions {0, 2, 3, 5, 7} are unchanged.
    outflow = np.setdiff1d(np.arange(n), inflow)
    np.testing.assert_array_equal(out[outflow], x[outflow])
    # Inflow positions are zero.
    np.testing.assert_array_equal(out[inflow], np.zeros(inflow.size))


@pytest.mark.l0
def test_self_adjointness():
    """L0: ``apply(x) == apply_transpose(x)`` for random ``x``."""
    rng = np.random.default_rng(seed=5)
    M = IncomingOrdinateMaskTensor(np.array([0, 3, 5]), n_ordinates=6)
    x = rng.standard_normal(6)
    np.testing.assert_array_equal(M.apply(x), M.apply_transpose(x))


@pytest.mark.l0
def test_idempotence():
    """L0: ``M(M(x)) == M(x)`` — projection onto the outflow subspace."""
    rng = np.random.default_rng(seed=11)
    M = IncomingOrdinateMaskTensor(np.array([2, 4]), n_ordinates=5)
    x = rng.standard_normal(5)
    once = M.apply(x)
    twice = M.apply(once)
    np.testing.assert_array_equal(once, twice)


@pytest.mark.l0
def test_empty_inflow_is_identity():
    """L0: empty ``inflow_indices`` makes the mask the identity."""
    rng = np.random.default_rng(seed=13)
    M = IncomingOrdinateMaskTensor(np.array([], dtype=np.intp), n_ordinates=6)
    x = rng.standard_normal(6)
    np.testing.assert_array_equal(M.apply(x), x)
    np.testing.assert_array_equal(M.apply_transpose(x), x)


@pytest.mark.l0
def test_full_inflow_matches_zero_operator():
    """L0: ``inflow_indices = arange(n_ordinates)`` matches :class:`ZeroOperator`."""
    n = 5
    rng = np.random.default_rng(seed=17)
    M = IncomingOrdinateMaskTensor(np.arange(n), n_ordinates=n)
    Z = ZeroOperator()
    x = rng.standard_normal(n)
    np.testing.assert_array_equal(M.apply(x), Z.apply(x))
    np.testing.assert_array_equal(M.apply(x), np.zeros_like(x))


@pytest.mark.l0
def test_compose_with_identity():
    """L0: ``(M @ Identity).apply(x) == M.apply(x)``."""
    M = IncomingOrdinateMaskTensor(np.array([1, 3]), n_ordinates=4)
    I = IdentityOperator()
    x = np.array([10.0, 20.0, 30.0, 40.0])
    np.testing.assert_array_equal((M @ I).apply(x), M.apply(x))
    np.testing.assert_array_equal((I @ M).apply(x), M.apply(x))


@pytest.mark.l0
def test_axis1_apply_2d():
    """L0: shape ``(3, 5)``, axis=1, ``inflow_indices=[1, 3]`` zeroes middle columns."""
    inflow = np.array([1, 3])
    M = IncomingOrdinateMaskTensor(inflow, n_ordinates=5, axis=1)
    x = np.array(
        [
            [1.0, 2.0, 3.0, 4.0, 5.0],
            [6.0, 7.0, 8.0, 9.0, 10.0],
            [11.0, 12.0, 13.0, 14.0, 15.0],
        ]
    )
    out = M.apply(x)
    expected = np.array(
        [
            [1.0, 0.0, 3.0, 0.0, 5.0],
            [6.0, 0.0, 8.0, 0.0, 10.0],
            [11.0, 0.0, 13.0, 0.0, 15.0],
        ]
    )
    np.testing.assert_array_equal(out, expected)


@pytest.mark.l0
def test_invalid_out_of_range_raises():
    """L0: out-of-range ``inflow_indices`` rejected at construction."""
    with pytest.raises(ValueError, match="out of range"):
        IncomingOrdinateMaskTensor(np.array([0, 5]), n_ordinates=4)


@pytest.mark.l0
def test_invalid_duplicates_raise():
    """L0: duplicate ``inflow_indices`` rejected at construction."""
    with pytest.raises(ValueError, match="duplicates"):
        IncomingOrdinateMaskTensor(np.array([0, 2, 2]), n_ordinates=4)


@pytest.mark.l0
def test_invalid_indices_2d_raises():
    """L0: non-1-D ``inflow_indices`` rejected."""
    with pytest.raises(ValueError, match="1-D"):
        IncomingOrdinateMaskTensor(
            np.array([[0, 1], [1, 0]]), n_ordinates=4
        )


@pytest.mark.l0
def test_input_unmodified():
    """L0: ``apply`` returns a fresh array — the input is unchanged."""
    M = IncomingOrdinateMaskTensor(np.array([0, 2]), n_ordinates=4)
    x = np.array([1.0, 2.0, 3.0, 4.0])
    original = x.copy()
    _ = M.apply(x)
    np.testing.assert_array_equal(x, original)


@pytest.mark.l0
def test_capability_set():
    """L0: capability advertises ``{CAP_APPLY, CAP_APPLY_TRANSPOSE}``."""
    M = IncomingOrdinateMaskTensor(np.array([0]), n_ordinates=4)
    assert M.capabilities == frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})
