r"""Scalar-trace flux displacement — the iterate increment of
:class:`~orpheus.transport.fields.scalar_boundary_flux.ScalarBoundaryFlux` (#290 P2).

The scalar-trace member of the displacement column: the difference
vector :math:`\Delta J^\pm = J^{\pm,(i)} \ominus J^{\pm,(i-1)}` between
two partial-current boundary states. Stored on the same cached
:class:`~orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace`
(``mesh.scalar_trace``) as its flux sibling; a DISTINCT class — a
*step*, not a *state* (the #208 affine-axiom argument: partial-current
states form an affine space, their increments the tangent vector
space). Born ONLY from ``ScalarBoundaryFlux ⊖ ScalarBoundaryFlux``; carries the
convergence diagnostics via
:class:`~orpheus.transport.displacements._displacement.Displacement`.
"""

from __future__ import annotations

from dataclasses import dataclass

from orpheus.transport.displacements._displacement import Displacement
from orpheus.transport.fields._bases import ScalarBoundaryField


__all__ = ["ScalarBoundaryDisplacement"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class ScalarBoundaryDisplacement(Displacement, ScalarBoundaryField):
    r"""The increment :math:`\Delta J^\pm` between two partial-current states.

    Shares the ``mesh.scalar_trace`` space, face layout, and area
    metric with :class:`~orpheus.transport.fields.scalar_boundary_flux.ScalarBoundaryFlux`
    (the shared Rep :class:`~orpheus.transport.fields._bases.ScalarBoundaryField`
    keys the flux↔displacement pairing); the vector-space algebra is
    inherited from :class:`~orpheus.numerics.field.Field` and the
    convergence diagnostics (contraction ratio, true-error estimate)
    from :class:`Displacement`.
    """
