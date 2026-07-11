r"""System B's interior q½ source block — the *source* role leaf of
:class:`~orpheus.transport.fields._bases.RadialCharacteristicInteriorField`.

The interior half of the ψ½ source split (the coupled-block campaign poses the
ray as **System B**, its own interior ⊕ boundary composite; B.2b re-types the
coupling operators onto it). The per-``(level, sign)`` cells blocks ``(ng, nx)``
carry the within-group emission :math:`\bar q_{1/2}` folded to the starting
direction (ruling R14: the full :math:`(-1)^\ell` Legendre fold) — a volumetric
rate density, exactly like the bulk
:class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`.
This is the block ``A_BA``
(:class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicEmission`,
the ``Fold ∘ K ∘ integrate`` emission) writes.

The split DISSOLVES the unified
:class:`~orpheus.transport.source_sinks.radial_characteristic_source_sink.RadialCharacteristicSourceSink`'s
documented corner-units deviation: each locus leaf now declares its own honest
dimensional identity (cells = volumetric rate here; the corner datum lives on
:class:`~orpheus.transport.source_sinks.radial_characteristic_boundary_source_sink.RadialCharacteristicBoundarySourceSink`
with flux units, the trace convention).

Storage, validation, the ``cells(level, sign)`` view, and the ``zeros_on`` /
``from_mesh`` factories are inherited from
:class:`~orpheus.transport.fields._bases.RadialCharacteristicInteriorField`.
Plain vector-space algebra (source sums are closed) — no role mixin, like every
SourceSink leaf; the Field class-identity gate keeps source, flux, and
displacement arithmetic from silently mixing (all three share the SAME
``mesh.radial_characteristic_interior_space``, so it is the class gate, not the
space gate, that rejects cross-role sums).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import ANGULAR_RATE_UNITS, Unit
from orpheus.transport.fields._bases import RadialCharacteristicInteriorField

__all__ = ["RadialCharacteristicInteriorSourceSink"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class RadialCharacteristicInteriorSourceSink(RadialCharacteristicInteriorField):
    r"""System B's interior q½ source/sink — the folded volumetric emission cells.

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

    #: The cells blocks carry the folded volumetric emission rate
    #: ``1/(cm³·s·sr)`` — :data:`~orpheus.numerics.units.ANGULAR_RATE_UNITS`,
    #: shared with ``AngularSourceSink``. Honest per-locus identity (the split
    #: dissolved the unified leaf's corner deviation). Metadata, not the gate.
    UNITS: ClassVar[Unit] = ANGULAR_RATE_UNITS
