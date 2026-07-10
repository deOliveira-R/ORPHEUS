r"""System B's boundary ψ½ flux — the r = R corner (Phase B split).

The *flux* role leaf of
:class:`~orpheus.transport.fields._bases.RadialCharacteristicBoundaryField`: the
boundary half of the ψ½ ray split (the coupled-block campaign poses the ray as
**System B**, its own interior ⊕ boundary composite). The per-``(level, sign)``
half-angle flux :math:`\psi_{1/2}` at the outer radius r = R — System B's BC
locus, on which
:class:`~orpheus.sn.operators.boundary.RadialCharacteristicBoundaryOperator` (B_b)
acts (inflow corner ``sign = -1`` = given BC data; outflow corner ``sign = +1`` =
the defect row, ruling R13).

Storage, validation, the ``corner(level, sign)`` view, and the ``zeros_on`` /
``from_mesh`` factories are inherited from
:class:`~orpheus.transport.fields._bases.RadialCharacteristicBoundaryField`; the
affine-point flux algebra from
:class:`~orpheus.transport.fields._flux_role.FluxRole` — its ``⊖`` mints a
:class:`~orpheus.transport.displacements.radial_characteristic_boundary_displacement.RadialCharacteristicBoundaryDisplacement`
(the Rep-keyed sibling over the shared ``RadialCharacteristicBoundaryField``).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import ANGULAR_FLUX_UNITS, Unit
from orpheus.transport.fields._bases import RadialCharacteristicBoundaryField
from orpheus.transport.fields._flux_role import FluxRole

__all__ = ["RadialCharacteristicBoundaryFlux"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class RadialCharacteristicBoundaryFlux(FluxRole, RadialCharacteristicBoundaryField):
    r"""System B's boundary ψ½ flux state — the r = R corner.

    Parameters
    ----------
    values : NDArray
        Flat backing buffer, shape ``(space.shape[0],)`` — per level and
        direction sign, the ``(ng,)`` corner.
    space : RadialCharacteristicBoundarySpace
        The R12a-keyed boundary space (canonically
        ``mesh.radial_characteristic_boundary_space``) carrying the layout and
        the ``G = V(r = R)`` corner gauge.
    mesh : SNMesh
        The SN phase-space carrier (the cross-mesh-arithmetic guard).
    """

    #: ψ½ IS an angular flux (evaluated at the μ = ±1 ray) — ``1/(cm²·s·sr)``.
    UNITS: ClassVar[Unit] = ANGULAR_FLUX_UNITS
