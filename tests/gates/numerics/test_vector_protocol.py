r"""foundation — the :class:`~orpheus.numerics.vector.Vector` Protocol.

Pins the structural contract that the operator algebra's vector type
speaks (``apply(x: V) -> V`` over a :data:`~orpheus.numerics.vector.V`
bound to :class:`Vector`). The campaign claim (#256) is that the SAME
Protocol is satisfied — without any shared inheritance ancestor — by

* ``np.ndarray`` (the serialization wire format + axis-primitive element),
* every :class:`~orpheus.numerics.field.Field` leaf, and
* :class:`~orpheus.transport.timed_full_field.TimedFullField` (the
  composite carrier),

while objects lacking the vector-space dunders are rejected. This module
holds the numerics half of those type-gates: the ``np.ndarray`` positive, the
negatives and the type variable's bound. The carrier positives (the transport
flux leaves and ``TimedFullField``) are transport claims and live in
``tests/gates/transport/fields/test_flux_vector_algebra.py``, section 5
(moved there with the layer-DAG re-filing, 2026-10-25).

This is a ``foundation`` gate — a software-invariant on the algebra's
type surface, with no theory-page ``:label:`` (so it carries no
``verifies(...)``).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.vector import V, Vector


pytestmark = [pytest.mark.foundation]


# ───────────────────────────────────────────────────────────────────────
# Positive — the wire format satisfies Vector for free
# ───────────────────────────────────────────────────────────────────────


def test_ndarray_is_vector() -> None:
    """The flat serialization wire format is a Vector."""
    assert isinstance(np.zeros(3), Vector)


# ───────────────────────────────────────────────────────────────────────
# Negative — missing any vector-space dunder is rejected
# ───────────────────────────────────────────────────────────────────────


def test_string_is_not_vector() -> None:
    """``str`` has ``__add__``/``__mul__`` but no ``__sub__`` — rejected."""
    assert not isinstance("not a vector", Vector)


def test_partial_dunders_rejected() -> None:
    """An object with only ``__add__`` is NOT a Vector (needs all three)."""

    class _AddOnly:
        def __add__(self, other: object) -> "_AddOnly":
            return self

    assert not isinstance(_AddOnly(), Vector)


def test_no_scalar_mul_rejected() -> None:
    """add + sub but no ``__rmul__`` is rejected — pins the scalar-mul contract.

    Guards the correctness fix that the algebra's scalar multiplication is
    ``scalar * vector`` (``__rmul__``), not ``vector * scalar``
    (``__mul__``): a carrier that cannot be left-multiplied by a scalar
    cannot flow through a scaled term of a combination (``α·A``) / ``ZeroOperator``.
    """

    class _NoScalarMul:
        def __add__(self, other: object) -> "_NoScalarMul":
            return self

        def __sub__(self, other: object) -> "_NoScalarMul":
            return self

    assert not isinstance(_NoScalarMul(), Vector)


def test_no_scalar_div_rejected() -> None:
    """add + sub + rmul but no ``__truediv__`` is rejected.

    Pins the scalar-division contract the eigenvalue drivers rely on
    (``F ψ / k`` in ``KEigenvalue``; ``ψ / p`` in ``power_iteration``): a
    carrier that cannot be divided by a scalar cannot flow through those
    renormalisations, so the step-3 ``apply(x: V) -> V`` retype of
    ``F ψ / k`` would not type-check for it.
    """

    class _NoScalarDiv:
        def __add__(self, other: object) -> "_NoScalarDiv":
            return self

        def __sub__(self, other: object) -> "_NoScalarDiv":
            return self

        def __rmul__(self, scalar: float) -> "_NoScalarDiv":
            return self

    assert not isinstance(_NoScalarDiv(), Vector)


# ───────────────────────────────────────────────────────────────────────
# The type variable
# ───────────────────────────────────────────────────────────────────────


def test_typevar_bound_to_vector() -> None:
    """``V`` is bound to ``Vector`` (the endomorphism ``apply(x: V) -> V``)."""
    assert V.__bound__ is Vector
