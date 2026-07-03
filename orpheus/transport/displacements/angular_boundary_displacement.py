r"""Boundary-trace flux displacement :math:`\Delta\psi|_\partial`.

The difference-space (tangent) sibling of
:class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux` — the increment
between two boundary-trace iterates (the ``boundary`` member of the SI iterate's
:class:`~orpheus.transport.timed_full_field.TimedFullField`). Stored on the
unified :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` (shared with the
boundary flux — like :class:`AngularBoundaryResidual`); a DISTINCT class. Born ONLY
from ``AngularBoundaryFlux ⊖ AngularBoundaryFlux``. See
:class:`~orpheus.transport.displacements._displacement.Displacement`.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import ANGULAR_FLUX_UNITS, Unit
from orpheus.transport.displacements._displacement import Displacement
from orpheus.transport.fields._bases import AngularBoundaryField


__all__ = ["AngularBoundaryDisplacement"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class AngularBoundaryDisplacement(Displacement, AngularBoundaryField):
    r"""L2 boundary-trace flux displacement on the unified ``AngularTraceSpace``.

    Like :class:`~orpheus.transport.residuals.angular_boundary_residual.AngularBoundaryResidual`,
    it shares the ``mesh.angular_trace`` :class:`AngularTraceSpace` with
    :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux` (the boundary
    leaves do not carry a distinct ``_SPACE_NAME``); the CLASS identity is the
    role gate.
    """

    #: Same units as :class:`AngularBoundaryFlux` (:math:`1/(\mathrm{cm^2 \cdot s \cdot
    #: sr})` — a trace flux-difference).
    UNITS: ClassVar[Unit] = ANGULAR_FLUX_UNITS
