r"""System B's boundary q½ source block — the *source* role leaf of
:class:`~orpheus.transport.fields._bases.RadialCharacteristicBoundaryField`.

The boundary half of the ψ½ source split: the per-``(level, sign)`` r = R
corner slots ``(ng,)``. The inflow corner (``sign = -1``) carries the
prescribed r = R inflow *datum* — an angular-flux value, the affine-BC ``q``
of the seed's identity row; the outflow corner (``sign = +1``) is the defect
row's source slot (ruling R13). This is the block ``B_b``
(:class:`~orpheus.sn.operators.boundary.RadialCharacteristicBoundaryOperator`,
the ray-corner ``A_ss`` action) emits.

Units — the trace convention
============================

Like the spatial trace sibling
:class:`~orpheus.transport.source_sinks.angular_boundary_source_sink.AngularBoundarySourceSink`,
the corner datum is FLUX-valued (``1/(cm²·s·sr)``): on a trace-like locus a
source does NOT pick up the volumetric ``cm⁻³``. The split DISSOLVES the
unified
:class:`~orpheus.transport.source_sinks.radial_characteristic_source_sink.RadialCharacteristicSourceSink`'s
documented corner-units deviation — each locus leaf now declares its own
honest dimensional identity (the volumetric cells live on
:class:`~orpheus.transport.source_sinks.radial_characteristic_interior_source_sink.RadialCharacteristicInteriorSourceSink`).

Storage, validation, the ``corner(level, sign)`` view, and the ``zeros_on`` /
``from_mesh`` factories are inherited from
:class:`~orpheus.transport.fields._bases.RadialCharacteristicBoundaryField`.
Plain vector-space algebra (source sums are closed) — no role mixin, like
every SourceSink leaf; the Field class-identity gate rejects cross-role sums
(source vs flux vs displacement all share the SAME
``mesh.radial_characteristic_boundary_space``).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import ANGULAR_FLUX_UNITS, Unit
from orpheus.transport.fields._bases import RadialCharacteristicBoundaryField

__all__ = ["RadialCharacteristicBoundarySourceSink"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class RadialCharacteristicBoundarySourceSink(RadialCharacteristicBoundaryField):
    r"""System B's boundary q½ source/sink — the r = R corner data slots.

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

    #: Trace-like corner datum — ``1/(cm²·s·sr)``,
    #: :data:`~orpheus.numerics.units.ANGULAR_FLUX_UNITS`, shared with
    #: ``AngularBoundarySourceSink`` (on a trace a source does not pick up the
    #: volumetric ``cm⁻³``). Honest per-locus identity (the split dissolved the
    #: unified leaf's corner deviation). Metadata, not the gate.
    UNITS: ClassVar[Unit] = ANGULAR_FLUX_UNITS
