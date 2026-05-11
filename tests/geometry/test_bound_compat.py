r"""Tests for the Wave-8 transitional 2-arg shim (C8.2).

The :class:`_BoundBoundaryOperator` shim wraps a realized 1-arg
:class:`LinearOperator` so the 13 production call sites at
:mod:`orpheus.sn.sweep` + :mod:`orpheus.sn.operator` (which still call
``bc.apply(psi, quad)`` with 2 args) keep working through Wave 9.

These tests pin:

* ``apply(psi, *anything, **anything)`` forwards to ``inner.apply(psi)``
  — extra args are swallowed.
* ``apply_transpose`` symmetrical for inner ops that support it.
* :attr:`capabilities` property delegates to the wrapped operator.
* The shim composes cleanly with the operator-algebra dunders inherited
  from :class:`LinearOperatorMixin` — proves it is a first-class
  :class:`LinearOperator`.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.boundary._bound_compat import _BoundBoundaryOperator
from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    IdentityOperator,
    PermutationOperator,
    ZeroOperator,
)


pytestmark = pytest.mark.l0


def test_apply_forwards_and_swallows_extra_args():
    """``apply(psi, *_extra, **_kw)`` delegates to ``inner.apply(psi)`` and
    ignores any positional / keyword extras (typically a quadrature
    argument passed by the legacy 2-arg call sites).
    """
    inner = IdentityOperator()
    shim = _BoundBoundaryOperator(inner)
    psi = np.arange(12.0).reshape(4, 3)
    # 1-arg
    np.testing.assert_array_equal(shim.apply(psi), psi)
    # 2-arg — extra "quadrature-like" sentinel is swallowed
    sentinel = object()
    np.testing.assert_array_equal(shim.apply(psi, sentinel), psi)
    # kwargs swallowed too
    np.testing.assert_array_equal(
        shim.apply(psi, quad=sentinel, extra=42), psi,
    )


def test_apply_transpose_forwards_when_inner_supports_it():
    """For an inner operator with ``apply_transpose`` (e.g. a
    permutation), the shim's :meth:`apply_transpose` returns the
    same result as ``inner.apply_transpose``.
    """
    perm = np.array([2, 0, 1, 3])
    inner = PermutationOperator(perm, axis=0)
    shim = _BoundBoundaryOperator(inner)
    psi = np.arange(8.0).reshape(4, 2)
    np.testing.assert_array_equal(
        shim.apply_transpose(psi),
        inner.apply_transpose(psi),
    )
    # And again ignoring extra args
    np.testing.assert_array_equal(
        shim.apply_transpose(psi, "ignored", k=1),
        inner.apply_transpose(psi),
    )


def test_capabilities_delegate_to_inner():
    """:attr:`capabilities` returns whatever the wrapped operator
    advertises — a permutation brings both apply + apply_transpose; a
    mask-only operator brings only apply.
    """
    perm_shim = _BoundBoundaryOperator(PermutationOperator(np.array([1, 0]), axis=0))
    assert perm_shim.capabilities == frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})

    zero_shim = _BoundBoundaryOperator(ZeroOperator())
    # ZeroOperator advertises {CAP_APPLY, CAP_APPLY_TRANSPOSE}
    assert CAP_APPLY in zero_shim.capabilities


def test_composes_with_operator_algebra():
    """The shim is a first-class :class:`LinearOperator` — it inherits
    the algebra dunders from :class:`LinearOperatorMixin`. ``shim + shim``
    builds an :class:`OperatorSum`; ``2.0 * shim`` builds a
    :class:`ScaledOperator`. This is the load-bearing piece for the
    realizer's mixed-BC recursive path (when MixedBoundaryOperator is
    composed with a wrapped realized op).
    """
    inner = IdentityOperator()
    shim = _BoundBoundaryOperator(inner)
    psi = np.ones((3, 2))
    # scalar multiplication via __mul__
    scaled = 2.0 * shim
    np.testing.assert_array_equal(scaled.apply(psi), 2.0 * psi)
    # sum via __add__ — note OperatorSum.apply forwards 1-arg only
    summed = shim + shim
    np.testing.assert_array_equal(summed.apply(psi), 2.0 * psi)


def test_shim_is_not_re_exported_from_package():
    """The shim is internal to the SN-side wiring — it MUST NOT appear
    in :mod:`orpheus.geometry.boundary`'s public surface, since no
    consumer outside ``SNMesh._resolve_bcs`` has a legitimate reason
    to wrap operators in it. Pinning the lack of re-export prevents
    accidental promotion to public API.
    """
    import orpheus.geometry.boundary as boundary_pkg

    assert "_BoundBoundaryOperator" not in getattr(boundary_pkg, "__all__", [])
    assert not hasattr(boundary_pkg, "_BoundBoundaryOperator")
