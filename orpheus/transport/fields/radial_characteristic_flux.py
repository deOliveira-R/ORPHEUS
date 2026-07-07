r"""The ψ½ starting-direction flux state — the *flux* role leaf of
:class:`~orpheus.transport.fields._bases.RadialCharacteristicField`.

The typed per-level half-angle angular flux :math:`\psi_{1/2}` of the
curvilinear Morel–Montry thread, promoted from a lagged solver-internal
estimate to composite STATE (#282 route (a), campaign #280 phase 2.5d):
the third block of the augmented carrier
:class:`~orpheus.transport.full_field.FullField` on seed-carrying meshes
(R12a: the sphere; see the
:class:`~orpheus.numerics.spaces.radial_characteristic_space.RadialCharacteristicSpace`
module docstring for the presence trichotomy).

Per ruling R13 the leaf carries BOTH direction legs per level — the
inward :math:`\mu = \mu_{\rm start}` state and the pole-continued
outward leg (:math:`\psi_{1/2}^+(0) = \psi_{1/2}^-(0)`) — each as an
``(ng, nx)`` cells block plus the ``(ng,)`` r = R corner slot (inflow
corner = given data; outflow corner = the defect row).

Storage, validation, algebra, the ``cells(level, sign)`` /
``corner(level, sign)`` views, and the ``zeros_on`` / ``from_mesh``
factories are inherited from :class:`RadialCharacteristicField`; the
affine-point flux algebra (``flux ⊖ flux → displacement``, torsor
``flux ⊕ displacement``, ``flux + flux → TypeError``) from
:class:`~orpheus.transport.fields._flux_role.FluxRole` — identical to
every other flux leaf.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import ANGULAR_FLUX_UNITS, Unit
from orpheus.transport.fields._bases import RadialCharacteristicField
from orpheus.transport.fields._flux_role import FluxRole

__all__ = ["RadialCharacteristicFlux"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class RadialCharacteristicFlux(FluxRole, RadialCharacteristicField):
    r"""L2 starting-direction flux — the ψ½ state block of the augmented composite.

    Parameters
    ----------
    values : NDArray
        Flat backing buffer, shape ``(space.shape[0],)`` — per level and
        direction sign, cells ``(ng, nx)`` ⊕ corner ``(ng,)``.
    space : RadialCharacteristicSpace
        The R12a-keyed
        :class:`~orpheus.numerics.spaces.radial_characteristic_space.RadialCharacteristicSpace`
        (canonically ``mesh.radial_characteristic_space``) — carries the
        seed levels, the layout arithmetic, and the SPD ``G_sd = V_cell``
        state metric.
    mesh : SNMesh
        The SN phase-space carrier (the cross-mesh-arithmetic guard).
    """

    #: Dimensional identity (View-G): ψ½ IS an angular flux (evaluated at
    #: the starting-direction ray), so ``1/(cm²·s·sr)`` — shared with
    #: ``AngularFlux`` / the boundary flux leaves. Metadata, not the gate.
    UNITS: ClassVar[Unit] = ANGULAR_FLUX_UNITS
