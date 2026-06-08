r"""Boundary-trace flux displacement :math:`\Delta\psi|_\partial`.

The difference-space (tangent) sibling of
:class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux` — the increment
between two boundary-trace iterates (the ``boundary`` member of the SI iterate's
:class:`~orpheus.transport.timed_full_field.TimedFullField`). Stored on the
unified :class:`~orpheus.numerics.spaces.trace_space.TraceSpace` (shared with the
boundary flux — like :class:`BoundaryResidual`); a DISTINCT class. Born ONLY
from ``BoundaryFlux ⊖ BoundaryFlux``. See
:class:`~orpheus.transport.displacements._displacement.Displacement`.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import ANGULAR_FLUX_UNITS, Unit
from orpheus.transport.displacements._displacement import Displacement
from orpheus.transport.fields._bases import BoundaryField


__all__ = ["BoundaryDisplacement"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class BoundaryDisplacement(Displacement, BoundaryField):
    r"""L2 boundary-trace flux displacement on the unified ``TraceSpace``.

    Like :class:`~orpheus.transport.residuals.boundary_residual.BoundaryResidual`,
    it shares the ``mesh.trace`` :class:`TraceSpace` with
    :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux` (the boundary
    leaves do not carry a distinct ``_SPACE_NAME``); the CLASS identity is the
    role gate.
    """

    #: Same units as :class:`BoundaryFlux` (:math:`1/(\mathrm{cm^2 \cdot s \cdot
    #: sr})` — a trace flux-difference).
    UNITS: ClassVar[Unit] = ANGULAR_FLUX_UNITS
