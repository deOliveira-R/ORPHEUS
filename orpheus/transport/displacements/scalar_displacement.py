r"""Scalar flux displacement :math:`\Delta\phi(\vec r, g)`.

The difference-space (tangent) sibling of
:class:`~orpheus.transport.fields.scalar_flux.ScalarFlux` — the increment
between two scalar-flux iterates. Same storage shape ``(ng, nx, ny)``, units,
and space as the flux; a DISTINCT class (the affine difference vector). Born
ONLY from ``ScalarFlux ⊖ ScalarFlux``. See
:class:`~orpheus.transport.displacements._displacement.Displacement`.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import SCALAR_FLUX_UNITS, Unit
from orpheus.transport.displacements._displacement import Displacement
from orpheus.transport.fields._bases import ScalarField


__all__ = ["ScalarDisplacement"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class ScalarDisplacement(Displacement, ScalarField):
    r"""L2 scalar flux displacement on ``(ng, nx, ny)``."""

    #: Shares the flux's ``"scalar_flux"`` space identity (tangent vector of the
    #: same function space). The class identity is the role gate.
    _SPACE_NAME: ClassVar[str] = "scalar_flux"

    #: Same units as :class:`ScalarFlux` (:math:`1/(\mathrm{cm^2 \cdot s})`).
    UNITS: ClassVar[Unit] = SCALAR_FLUX_UNITS
