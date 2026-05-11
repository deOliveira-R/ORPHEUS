r"""Tests for :class:`orpheus.numerics.operator.PeriodicWrapOperator`.

Verifies the operator's invariants:

* ``apply`` is identity (returns the input by reference, matching
  the legacy zero-cost angular pass-through).
* ``apply_transpose`` is identity.
* Capability set: ``{CAP_APPLY, CAP_APPLY_TRANSPOSE}``.
* Composition with :class:`IdentityOperator` is identity.
* The input is returned by reference (not copied) — matching the
  legacy :class:`~orpheus.geometry.boundary.PeriodicBoundaryOperator`
  semantics.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    IdentityOperator,
    PeriodicWrapOperator,
)


@pytest.mark.l0
def test_apply_is_identity_various_shapes():
    """L0: ``apply(x) == x`` for various input shapes."""
    P = PeriodicWrapOperator()
    rng = np.random.default_rng(seed=2)
    for shape in [(5,), (3, 4), (2, 3, 5)]:
        x = rng.standard_normal(shape)
        np.testing.assert_array_equal(P.apply(x), x)


@pytest.mark.l0
def test_apply_transpose_is_identity():
    """L0: ``apply_transpose(x) == x`` (identity body is self-adjoint)."""
    P = PeriodicWrapOperator()
    rng = np.random.default_rng(seed=4)
    x = rng.standard_normal((3, 4))
    np.testing.assert_array_equal(P.apply_transpose(x), x)


@pytest.mark.l0
def test_capability_set():
    """L0: capability advertises ``{CAP_APPLY, CAP_APPLY_TRANSPOSE}``."""
    P = PeriodicWrapOperator()
    assert P.capabilities == frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})


@pytest.mark.l0
def test_compose_with_identity():
    """L0: ``(P @ Identity).apply(x) == x`` and the other direction."""
    P = PeriodicWrapOperator()
    I = IdentityOperator()
    rng = np.random.default_rng(seed=6)
    x = rng.standard_normal(7)
    np.testing.assert_array_equal((P @ I).apply(x), x)
    np.testing.assert_array_equal((I @ P).apply(x), x)


@pytest.mark.l0
def test_apply_returns_fresh_copy():
    """L0: ``apply`` returns a fresh copy of the input.

    The legacy ``PeriodicBoundary.apply`` body is ``psi_out.copy()``;
    this primitive matches the safe-aliasing contract so that
    callers may freely mutate the returned array without affecting
    the caller's ``psi_out`` buffer. Wave 7 updated this primitive
    to perform the copy (originally returned by reference, which
    silently violated the legacy semantics — exposed when Wave 7
    delegation routed
    ``tests/geometry/test_boundary.py::test_periodic_bc_returns_input_unchanged``
    through this operator).
    """
    P = PeriodicWrapOperator()
    x = np.array([1.0, 2.0, 3.0])
    # Output is a distinct ndarray.
    assert P.apply(x) is not x
    assert P.apply_transpose(x) is not x
    # Values agree.
    np.testing.assert_array_equal(P.apply(x), x)
    np.testing.assert_array_equal(P.apply_transpose(x), x)
    # Mutating the output does not affect the input.
    out = P.apply(x)
    out[0] = 999.0
    assert x[0] == 1.0
