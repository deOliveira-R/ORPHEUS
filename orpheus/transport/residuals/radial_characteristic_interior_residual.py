r"""System B's interior ψ½ residual block — the *residual* role leaf of
:class:`~orpheus.transport.fields._bases.RadialCharacteristicInteriorField`.

The interior half of the ψ½ residual split (B.2d, mirroring the B.2b
SourceSink split): the marched-cells slice of System B's typed equation
residual :math:`r_B = (A\,\psi)_B - q_B`. On a carrying mesh (R12a) the
coupled loss grid's B-row emits ray rows (the μ = ∓1 DD balance residuals),
the coupled source carries the q½ member, and their difference is THIS leaf
— minted ONLY by :func:`orpheus.sn.solver.evaluate_residual`'s named balance.

Why the split leaf exists (Mode 12 (b) — the full-system residual norm)
=======================================================================

The §16.C C(i) discipline measures the equation residual over the FULL
coupled state — a System-A-only residual would be structurally blind to a
wrong seed row (vv-principles Mode 12: the restricted norm's invariance
group contains the seed-row error class). Typing System B's defect as its
own composite member keeps ``evaluate_residual`` honest by construction —
and, since the split, each locus leaf declares its own honest dimensional
identity (this cells leaf is a volumetric rate; the corner rows live on
:class:`~orpheus.transport.residuals.radial_characteristic_boundary_residual.RadialCharacteristicBoundaryResidual`
with the trace convention), dissolving the retired unified
``RadialCharacteristicResidual``'s documented corner-units deviation.

Storage, validation, the ``cells(level, sign)`` view, and the shared
face machinery are inherited from
:class:`~orpheus.transport.fields._bases.RadialCharacteristicInteriorField`.
Same-class arithmetic is closed by the Field Layer-1 gate; the balance that
forms a residual goes through :meth:`from_balance`, never a bare cross-class
``−``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import ANGULAR_RATE_UNITS, Unit
from orpheus.transport.fields._bases import RadialCharacteristicInteriorField

__all__ = ["RadialCharacteristicInteriorResidual"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class RadialCharacteristicInteriorResidual(RadialCharacteristicInteriorField):
    r"""The marched-cells block of System B's typed equation residual.

    Formed ONLY by the named balance :meth:`from_balance` (the residual
    discipline — an operator output is a source/sink; the residual arises
    when it is compared against the source and the difference is taken).

    Parameters
    ----------
    values : NDArray
        Flat backing buffer, shape ``(space.shape[0],)`` — per level and
        direction sign, the ``(ng, nx)`` cells.
    space : RadialCharacteristicInteriorSpace
        The R12a-keyed interior space (canonically
        ``mesh.radial_characteristic_interior_space``).
    mesh : SNMesh
        The SN phase-space carrier (the cross-mesh-arithmetic guard).
    """

    #: Per-cells-leg rate density ``1/(cm³·s·sr)`` — the μ = ∓1 DD balance
    #: rows share the bulk residual's dimensional signature (and the
    #: interior SourceSink's). Metadata, not the gate.
    UNITS: ClassVar[Unit] = ANGULAR_RATE_UNITS

    @classmethod
    def from_balance(
        cls,
        lhs: RadialCharacteristicInteriorField,
        rhs: RadialCharacteristicInteriorField,
    ) -> "RadialCharacteristicInteriorResidual":
        r"""The interior ψ½ balance defect ``r½ = (A·ψ)½ − q½`` (cells legs).

        ``lhs`` is the coupled loss grid's B-row interior member (a
        :class:`~orpheus.transport.source_sinks.radial_characteristic_interior_source_sink.RadialCharacteristicInteriorSourceSink`),
        ``rhs`` the coupled source's q½ interior member — the same three
        guards (same-class operands, units-exact, same space + mesh) as the
        bulk/boundary ``from_balance`` factories, via
        :meth:`~orpheus.numerics.field.Field._from_balance`.
        """
        return cls._from_balance(lhs, rhs)
