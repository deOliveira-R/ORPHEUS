r"""System B's boundary ψ½ residual block — the *residual* role leaf of
:class:`~orpheus.transport.fields._bases.RadialCharacteristicBoundaryField`.

The boundary half of the ψ½ residual split (B.2d, mirroring the B.2b
SourceSink split): the r = R corner slice of System B's typed equation
residual :math:`r_B = (A\,\psi)_B - q_B` — the defect of the corner
flux-matching rows (the inflow identity / outflow rows of the ray march).
Minted ONLY by :func:`orpheus.sn.solver.evaluate_residual`'s named balance;
its interior sibling is
:class:`~orpheus.transport.residuals.radial_characteristic_interior_residual.RadialCharacteristicInteriorResidual`.

Storage, validation, the ``corner(level, sign)`` view, and the shared
face machinery are inherited from
:class:`~orpheus.transport.fields._bases.RadialCharacteristicBoundaryField`.
Same-class arithmetic is closed by the Field Layer-1 gate; the balance that
forms a residual goes through :meth:`from_balance`, never a bare cross-class
``−``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import ANGULAR_FLUX_UNITS, Unit
from orpheus.transport.fields._bases import RadialCharacteristicBoundaryField

__all__ = ["RadialCharacteristicBoundaryResidual"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class RadialCharacteristicBoundaryResidual(RadialCharacteristicBoundaryField):
    r"""The r = R corner block of System B's typed equation residual.

    Formed ONLY by the named balance :meth:`from_balance` (the residual
    discipline).

    Parameters
    ----------
    values : NDArray
        Flat backing buffer, shape ``(space.shape[0],)`` — per level and
        direction sign, the ``(ng,)`` corner.
    space : RadialCharacteristicBoundarySpace
        The R12a-keyed boundary space (canonically
        ``mesh.radial_characteristic_boundary_space``).
    mesh : SNMesh
        The SN phase-space carrier (the cross-mesh-arithmetic guard).
    """

    #: Trace-like corner rows — ``1/(cm²·s·sr)``,
    #: :data:`~orpheus.numerics.units.ANGULAR_FLUX_UNITS` (on a trace a
    #: defect does not pick up the volumetric ``cm⁻³``; the B.2b units
    #: ruling for the SourceSink split applies identically). Metadata, not
    #: the gate.
    UNITS: ClassVar[Unit] = ANGULAR_FLUX_UNITS

    @classmethod
    def from_balance(
        cls,
        lhs: RadialCharacteristicBoundaryField,
        rhs: RadialCharacteristicBoundaryField,
    ) -> "RadialCharacteristicBoundaryResidual":
        r"""The corner ψ½ balance defect ``r½ = (A·ψ)½ − q½`` (corner slots).

        ``lhs`` is the coupled loss grid's B-row boundary member (a
        :class:`~orpheus.transport.source_sinks.radial_characteristic_boundary_source_sink.RadialCharacteristicBoundarySourceSink`),
        ``rhs`` the coupled source's q½ boundary member — the same three
        guards as every ``from_balance`` factory, via
        :meth:`~orpheus.numerics.field.Field._from_balance`.
        """
        return cls._from_balance(lhs, rhs)
