r"""Starting-direction residual field — the ψ½ block's balance defect.

The *residual* role leaf of
:class:`~orpheus.transport.fields._bases.RadialCharacteristicField` (#282
route (a), #280 Phase 2.5d): the starting-direction slice of the typed
equation residual :math:`r = A\,\psi - q`.  On a carrying mesh (R12a)
the augmented operator :math:`A = L + C - S - B` emits ψ½ rows (the
μ = ∓1 DD residuals + the corner rows), the source composite carries a
q½ block, and their difference is THIS leaf — completing the fourth
member of the role triple minted at 2.5d d1
(:class:`~orpheus.transport.fields.radial_characteristic_flux.RadialCharacteristicFlux`
/
:class:`~orpheus.transport.source_sinks.radial_characteristic_source_sink.RadialCharacteristicSourceSink`
/
:class:`~orpheus.transport.displacements.radial_characteristic_displacement.RadialCharacteristicDisplacement`),
forced by :func:`orpheus.sn.solver.evaluate_residual` exactly as the
Displacement was forced by the composite torsor algebra.

Why this leaf exists (Mode 12 (b) — the full-field residual norm)
=================================================================

The §16.C C(i) discipline measures the equation residual over the FULL
augmented field — a bulk⊕trace-only residual would be structurally
blind to a wrong seed row (vv-principles Mode 12: the restricted norm's
invariance group contains the seed-row error class).  Typing the seed
defect keeps ``evaluate_residual`` honest by construction: the block
cannot be silently dropped in the composite rebuild.

Units — the cells legs are volumetric rate densities (the μ = ∓1 DD
balance rows, same dimensional signature as the bulk residual /
``RadialCharacteristicSourceSink``); the two corner slots are trace-like
flux-matching rows (the same deviation the SourceSink leaf documents).
The gate is class identity, not units.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import ANGULAR_RATE_UNITS, Unit
from orpheus.transport.fields._bases import RadialCharacteristicField

__all__ = ["RadialCharacteristicResidual"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class RadialCharacteristicResidual(RadialCharacteristicField):
    r"""The starting-direction block of a typed equation residual.

    Formed ONLY by the named balance :meth:`from_balance` (the residual
    discipline — an operator output is a source/sink; the residual
    arises when it is compared against the source and the difference is
    taken).  Same-class arithmetic is closed by the Field Layer-1 gate.
    """

    _SPACE_NAME: ClassVar[str] = "radial_characteristic_residual"

    #: Per-cells-leg rate density (shared with the SourceSink leaf; the
    #: corner slots are trace-like flux rows — documented deviation).
    UNITS: ClassVar[Unit] = ANGULAR_RATE_UNITS

    @classmethod
    def from_balance(
        cls, lhs: RadialCharacteristicField, rhs: RadialCharacteristicField,
    ) -> "RadialCharacteristicResidual":
        r"""The ψ½-block balance defect ``r½ = (A·ψ)½ − q½``.

        ``lhs`` is the operator output's starting-direction block (a
        :class:`RadialCharacteristicSourceSink`), ``rhs`` the source
        composite's q½ block — the same three guards (same-class
        operands, units-exact, same space + mesh) as the bulk/boundary
        ``from_balance`` factories, via
        :meth:`~orpheus.numerics.field.Field._from_balance`.
        """
        return cls._from_balance(lhs, rhs)
