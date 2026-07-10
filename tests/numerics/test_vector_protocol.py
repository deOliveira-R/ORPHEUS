r"""foundation — the :class:`~orpheus.numerics.vector.Vector` Protocol.

Pins the structural contract that the operator algebra's vector type
speaks (``apply(x: V) -> V`` over a :data:`~orpheus.numerics.vector.V`
bound to :class:`Vector`). The campaign claim (#256) is that the SAME
Protocol is satisfied — without any shared inheritance ancestor — by

* ``np.ndarray`` (the serialization wire format + axis-primitive element),
* every :class:`~orpheus.numerics.field.Field` leaf, and
* :class:`~orpheus.transport.timed_full_field.TimedFullField` (the
  composite carrier),

while objects lacking the vector-space dunders are rejected. These are
the positive + negative type-gates of campaign step 1.

This is a ``foundation`` gate — a software-invariant on the algebra's
type surface, with no theory-page ``:label:`` (so it carries no
``verifies(...)``).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.vector import V, Vector
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.timed_full_field import TimedFullField

from tests.sn._test_helpers import placeholder_materials

pytestmark = [pytest.mark.foundation]


def _slab_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


# ───────────────────────────────────────────────────────────────────────
# Positive — the three carrier families satisfy Vector for free
# ───────────────────────────────────────────────────────────────────────


def test_ndarray_is_vector() -> None:
    """The flat serialization wire format is a Vector."""
    assert isinstance(np.zeros(3), Vector)


def test_field_leaves_are_vectors() -> None:
    """Every Field leaf (bulk + boundary, flux + moment) satisfies Vector.

    Covers all three storage families: angular (``AngularFlux``), scalar
    (``ScalarFlux``), moment (``HarmonicMomentFlux`` — the carrier #256
    closes the moment-state path around), and the boundary locus
    (``AngularBoundaryFlux``).
    """
    m = _slab_mesh()
    assert isinstance(AngularFlux.zeros_on(m), Vector)
    assert isinstance(ScalarFlux.zeros_on(m), Vector)
    assert isinstance(HarmonicMomentFlux.zeros_for_mesh_and_L(m, 1), Vector)
    assert isinstance(AngularBoundaryFlux.zeros_on(m), Vector)


def test_timed_full_field_is_vector() -> None:
    """The composite carrier satisfies Vector via its delegate dunders.

    This is the load-bearing positive case: ``TimedFullField`` conforms
    to a numerics-level Protocol WITHOUT ``numerics`` importing
    ``transport`` (the layering payoff of a structural Protocol).
    """
    m = _slab_mesh()
    state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=m)
    assert isinstance(state, Vector)


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
    cannot flow through ``ScaledOperator`` / ``ZeroOperator``.
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
