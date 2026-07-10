r"""System B interior ψ½ flux displacement :math:`\Delta\psi_{1/2}` (Phase B split).

The difference-space (tangent) sibling of
:class:`~orpheus.transport.fields.radial_characteristic_interior_flux.RadialCharacteristicInteriorFlux`
— the increment between two interior ψ½ iterates. Stored on the shared R12a-keyed
:class:`~orpheus.numerics.spaces.radial_characteristic_space.RadialCharacteristicInteriorSpace`
(``mesh.radial_characteristic_interior_space``); a DISTINCT class. Born ONLY from
``RadialCharacteristicInteriorFlux ⊖ RadialCharacteristicInteriorFlux`` (the
:class:`~orpheus.transport.fields._flux_role.FluxRole` mint, resolved through the
Rep-keyed
:meth:`~orpheus.transport.displacements._displacement.Displacement.sibling_of`
registry — the split interior field is its own ``_carrier_rep``, so this
displacement pairs with the interior flux and NOT with the unified ψ½ or the
boundary displacement).

Its existence is FORCED by the composite torsor algebra: System B's flux
composite ``ψ₂ − ψ₁`` mints a displacement PER BLOCK — the interior block
included. See
:class:`~orpheus.transport.displacements._displacement.Displacement`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import ANGULAR_FLUX_UNITS, Unit
from orpheus.transport.displacements._displacement import Displacement
from orpheus.transport.fields._bases import RadialCharacteristicInteriorField

__all__ = ["RadialCharacteristicInteriorDisplacement"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class RadialCharacteristicInteriorDisplacement(
    Displacement, RadialCharacteristicInteriorField,
):
    r"""System B interior ψ½ flux displacement on the shared interior space.

    Like every displacement leaf it shares its flux sibling's space (the tangent
    vector lives in the flux's own function space — the SPD ``G_sd = V_cell``
    state metric rides along); the CLASS identity is the role gate.
    """

    #: Same units as :class:`RadialCharacteristicInteriorFlux` (a ψ½
    #: flux-difference) — ``1/(cm²·s·sr)``.
    UNITS: ClassVar[Unit] = ANGULAR_FLUX_UNITS
