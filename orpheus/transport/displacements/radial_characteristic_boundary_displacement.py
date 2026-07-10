r"""System B boundary ψ½ flux displacement :math:`\Delta\psi_{1/2}` (Phase B split).

The difference-space (tangent) sibling of
:class:`~orpheus.transport.fields.radial_characteristic_boundary_flux.RadialCharacteristicBoundaryFlux`
— the increment between two boundary (r = R corner) ψ½ iterates. Stored on the
shared R12a-keyed
:class:`~orpheus.numerics.spaces.radial_characteristic_space.RadialCharacteristicBoundarySpace`
(``mesh.radial_characteristic_boundary_space``); a DISTINCT class. Born ONLY from
``RadialCharacteristicBoundaryFlux ⊖ RadialCharacteristicBoundaryFlux`` (the
:class:`~orpheus.transport.fields._flux_role.FluxRole` mint, resolved through the
Rep-keyed
:meth:`~orpheus.transport.displacements._displacement.Displacement.sibling_of`
registry — the split boundary field is its own ``_carrier_rep``, so this
displacement pairs with the boundary flux and NOT with the unified ψ½ or the
interior displacement).

Its existence is FORCED by the composite torsor algebra: System B's flux
composite ``ψ₂ − ψ₁`` mints a displacement PER BLOCK — the boundary block
included. See
:class:`~orpheus.transport.displacements._displacement.Displacement`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import ANGULAR_FLUX_UNITS, Unit
from orpheus.transport.displacements._displacement import Displacement
from orpheus.transport.fields._bases import RadialCharacteristicBoundaryField

__all__ = ["RadialCharacteristicBoundaryDisplacement"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class RadialCharacteristicBoundaryDisplacement(
    Displacement, RadialCharacteristicBoundaryField,
):
    r"""System B boundary ψ½ flux displacement on the shared boundary space.

    Like every displacement leaf it shares its flux sibling's space (the tangent
    vector lives in the flux's own function space — the ``G = V(r = R)`` corner
    gauge rides along); the CLASS identity is the role gate.
    """

    #: Same units as :class:`RadialCharacteristicBoundaryFlux` (a ψ½
    #: flux-difference) — ``1/(cm²·s·sr)``.
    UNITS: ClassVar[Unit] = ANGULAR_FLUX_UNITS
