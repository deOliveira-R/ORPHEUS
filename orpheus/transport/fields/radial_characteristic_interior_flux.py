r"""System B's interior ψ½ flux — the marched cells (Phase B split).

The *flux* role leaf of
:class:`~orpheus.transport.fields._bases.RadialCharacteristicInteriorField`: the
interior half of the ψ½ ray split (the coupled-block campaign poses the ray as
**System B**, its own interior ⊕ boundary composite). The per-``(level, sign)``
half-angle flux :math:`\psi_{1/2}` at every radial cell — the marched interior
state that
:class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator`
(A_BB) evolves, A_AB reads (inward leg), and A_BA writes.

Storage, validation, the ``cells(level, sign)`` view, and the ``zeros_on`` /
``from_mesh`` factories are inherited from
:class:`~orpheus.transport.fields._bases.RadialCharacteristicInteriorField`; the
affine-point flux algebra (``flux ⊖ flux → displacement``, torsor
``flux ⊕ displacement``, ``flux + flux → TypeError``) from
:class:`~orpheus.transport.fields._flux_role.FluxRole` — its ``⊖`` mints a
:class:`~orpheus.transport.displacements.radial_characteristic_interior_displacement.RadialCharacteristicInteriorDisplacement`
(the Rep-keyed sibling over the shared ``RadialCharacteristicInteriorField``).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import ANGULAR_FLUX_UNITS, Unit
from orpheus.transport.fields._bases import RadialCharacteristicInteriorField
from orpheus.transport.fields._flux_role import FluxRole

__all__ = ["RadialCharacteristicInteriorFlux"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class RadialCharacteristicInteriorFlux(FluxRole, RadialCharacteristicInteriorField):
    r"""System B's interior ψ½ flux state — the marched cells.

    Parameters
    ----------
    values : NDArray
        Flat backing buffer, shape ``(space.shape[0],)`` — per level and
        direction sign, the ``(ng, nx)`` cells.
    space : RadialCharacteristicInteriorSpace
        The R12a-keyed interior space (canonically
        ``mesh.radial_characteristic_interior_space``) carrying the layout and
        the SPD ``G_sd = V_cell`` state metric.
    mesh : SNMesh
        The SN phase-space carrier (the cross-mesh-arithmetic guard).
    """

    #: ψ½ IS an angular flux (evaluated at the μ = ±1 ray) — ``1/(cm²·s·sr)``,
    #: shared with ``AngularFlux`` and the unified ψ½ leaf. Metadata, not the gate.
    UNITS: ClassVar[Unit] = ANGULAR_FLUX_UNITS
