r"""Scalar boundary-trace source/sink — the operator-output trace leaf.

The *source/sink* role leaf of the SCALAR boundary family
(:class:`~orpheus.transport.fields._bases.ScalarBoundaryField`), minted
at #290 P4 when the first operator codomain demanded it (the
consumer-driven discipline recorded at P2: "source-sink / residual role
siblings only as operator codomains require"). It mirrors
:class:`~orpheus.transport.source_sinks.angular_boundary_source_sink.AngularBoundarySourceSink`
one-for-one on the scalar trace:

* every scalar-composite operator's ``.apply`` output boundary IS a
  :class:`ScalarBoundarySourceSink`, because the operator output is
  ``Aψ`` — a source/sink, not a state. The diffusion loss family emits
  it as follows (#290 P4):

  - :class:`~orpheus.diffusion.operators.LeakageOperator` (``L``, FULL)
    — the non-zero trace block: the outflow-definition defect on the
    :attr:`~orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace.OUTFLOW_ROW`
    plus the inflow identity on the ``INFLOW_ROW``;
  - :class:`~orpheus.diffusion.operators.DiffusionBoundaryOperator`
    (``B``, BOUNDARY) — ``𝒜·J⁺`` on the ``INFLOW_ROW`` only;
  - ``C`` / ``S`` / ``F`` (BULK) — the auto-allocated zero.

  The loss matvec ``(L + C − S − B).apply`` then composes as a CLOSED
  ``ScalarBoundarySourceSink`` sum, exactly like SN's closed
  ``AngularBoundarySourceSink`` sums (B.5.2 / #208).

Units
=====

The scalar trace is **all-current**: the ``(J⁺, J⁻)`` rows of a trace
defect / boundary source are partial currents, the same areal-rate
family :math:`1/(\mathrm{cm^2\,s})` as
:class:`~orpheus.transport.fields.scalar_boundary_flux.ScalarBoundaryFlux`
(:data:`~orpheus.numerics.units.SCALAR_FLUX_UNITS`) — the trace, unlike
the bulk, picks up no volumetric ``cm⁻³`` (the angular precedent:
"the boundary is all-flux", B.4). Same units, different ROLE — the
Field Layer-1 class gate is what keeps
``ScalarBoundarySourceSink ± ScalarBoundaryFlux`` from silently mixing.

Role grid (3-axis, #290 P2.5): locus Boundary × family Scalar × role
SourceSink. The residual sibling joins when a diffusion balance
consumer demands it (``from_balance`` — the SN
``AngularBoundaryResidual`` precedent).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import SCALAR_FLUX_UNITS, Unit
from orpheus.transport.fields._bases import ScalarBoundaryField


__all__ = ["ScalarBoundarySourceSink"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class ScalarBoundarySourceSink(ScalarBoundaryField):
    r"""Scalar boundary-trace source/sink — ``(J⁺, J⁻)``-row rate pairs
    per face per group.

    Parameters
    ----------
    values : NDArray
        Flat backing buffer, shape ``(space.layout.total_size,)``.
    space : ScalarTraceSpace
        The mesh's cached scalar trace
        (:attr:`DiffusionMesh.scalar_trace
        <orpheus.diffusion.augmented_mesh.DiffusionMesh.scalar_trace>`).
    mesh : DiffusionMesh
        The diffusion phase space (cross-mesh guard; the family's
        covariant narrowing, #290 P7a).

    Notes
    -----
    A thin role leaf: storage, validation, algebra, per-face access
    (:meth:`~orpheus.transport.fields._bases.BoundaryField.face_view`)
    and the
    :meth:`~orpheus.transport.fields._bases.BoundaryField.from_face_arrays`
    packer are all inherited from :class:`ScalarBoundaryField`. The
    leaf carries no source-specific behaviour beyond its class identity
    — all three scalar boundary leaves share the SAME
    :class:`~orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace`
    (``mesh.scalar_trace``), so it is the **class** gate (not the space
    gate) that rejects ``ScalarBoundarySourceSink ± ScalarBoundaryFlux``.
    Same-class arithmetic is closed (the loss-sum composition).
    """

    #: Dimensional identity: areal rate ``1/(cm²·s)`` — the trace is
    #: all-current, shared with ``ScalarBoundaryFlux`` (same units,
    #: different role → class gate). Metadata, not the gate.
    UNITS: ClassVar[Unit] = SCALAR_FLUX_UNITS
