r"""Per-ordinate flux displacement :math:`\Delta\psi(\vec r, \hat\Omega, g)`.

The difference-space (tangent) sibling of
:class:`~orpheus.transport.fields.angular_flux.AngularFlux`: the iterate
increment :math:`\Delta\psi = \psi^{(i)} \ominus \psi^{(i-1)}` in the SN
source iteration. Same storage shape ``(N, ng, nx, ny)``, same units, same
:class:`~orpheus.numerics.space.FunctionSpace` (and metric) as the flux — but a
DISTINCT class (the affine difference vector, not a state). See
:class:`~orpheus.transport.displacements._displacement.Displacement` and
:class:`~orpheus.transport.fields._flux_role.FluxRole`.

A thin leaf — born ONLY from ``AngularFlux ⊖ AngularFlux`` (the
:class:`FluxRole` mint); all storage / algebra / factories are inherited from
:class:`~orpheus.transport.fields._bases.AngularField`. Shares the flux's
``"angular_flux"`` space identity (the tangent vector lives in the flux's own
function space) — the CLASS is the role gate, not the space name.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import ANGULAR_FLUX_UNITS, Unit
from orpheus.transport.displacements._displacement import Displacement
from orpheus.transport.fields._bases import AngularField


__all__ = ["AngularDisplacement"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class AngularDisplacement(Displacement, AngularField):
    r"""L2 per-ordinate flux displacement on ``(N, ng, nx, ny)``."""

    #: Shares the flux's space identity (``"angular_flux"``) — the displacement
    #: is the tangent vector of that same function space (same quadrature
    #: metric). The CLASS identity, not the space name, is the role gate.
    _SPACE_NAME: ClassVar[str] = "angular_flux"

    #: Same units as the flux it differences (a flux − flux is a flux-difference,
    #: :math:`1/(\mathrm{cm^2 \cdot s \cdot sr})`).
    UNITS: ClassVar[Unit] = ANGULAR_FLUX_UNITS
