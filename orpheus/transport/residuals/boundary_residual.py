r"""Boundary-trace residual field :math:`r_\Gamma(\vec r_\Gamma, \hat\Omega, g)`.

The L2 typed wrapper for the **defect of the boundary balance** — the
residual of the affine boundary law

.. math::

    r_\Gamma \;=\; \gamma_-\psi \;-\; \bigl(R\,G\,\gamma_+\psi + q\bigr),

stored as a field on the unified
:class:`~orpheus.numerics.spaces.trace_space.TraceSpace` (flat
``(layout.total_size,)`` buffer). It is the *residual* role leaf of
:class:`~orpheus.transport.fields._bases.BoundaryField`, sibling to
:class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux` (flux)
and :class:`~orpheus.transport.source_sinks.boundary_source_sink.BoundarySourceSink`
(source).

Consumer: the matvec ALREADY computes this (it is mistyped today)
=================================================================

Unlike its source sibling, the boundary residual has a **real,
already-computed consumer**. The SN matvec emits the affine-BC face
defect (``γ₊ψ − bc_estimate``, the residual of ``γ₋ψ = R·G·γ₊ψ + q``)
on its output ``.boundary`` — :class:`StreamingOperator` is the only
operator that emits a non-zero face residual; :class:`CollisionOperator`
/ :class:`ScatteringOperator` / :class:`FissionOperator` emit the
auto-allocated zero. GMRES already drives this boundary defect to zero
as part of the full flat residual vector (via
:meth:`TimedFullField.to_flat`). The relevant emission sites are
``orpheus/sn/operator.py`` (the streaming face-residual build + the
zero-boundary returns of the collision/scatter/fission peers, and the
``StreamingOperator.apply`` output at ``:1263``).

**Today those sites mistype the boundary defect as
:class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`** — the
boundary half of the "dimensional sin" (operator OUTPUTS, which are
rate-density residuals, wrapped in the FLUX type). Retyping them to
:class:`BoundaryResidual` is the boundary half of the **B.5** operator-
output carve; it is INSEPARABLE from the bulk half (the
:class:`TimedFullField` composition gate requires every leaf's output
boundary to share a class, and the typed arithmetic at
``iteration.py:455`` / ``:511`` decides all operator-output types
together). This leaf is minted here (B.3) so it is **ready** for the
B.5 wiring; the wiring itself — with its ``from_balance`` named
composition and a test-architect verification plan — lands in B.5.

Units (B.4 — declared as the ``UNITS`` class constant)
======================================================

The boundary balance is a flux-matching equation, so its defect is
flux-typed: :math:`[1/(\mathrm{cm^2 \cdot s \cdot sr})]`
(:data:`~orpheus.numerics.units.ANGULAR_FLUX_UNITS`) — the SAME as
``BoundaryFlux`` / ``BoundarySourceSink``, and NOT the volumetric
``cm⁻³`` of the *bulk* residual (the trace is all-flux). eV-free per
the binned-energy convention. Same units, different role — the gate is
class identity. See :mod:`orpheus.numerics.units`.

References
==========

* ``.claude/plans/field_role_typing_view_g.md`` — Phase B step B.3
  (mint) and B.5 (the matvec operator-output retype + ``from_balance``).
* Grand Report v3 §16A.1 (affine boundary form), §5.5 (Field hierarchy).
* ``coding-elegance`` Pattern 4 (illegal states unrepresentable).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import ANGULAR_FLUX_UNITS, Unit
from orpheus.transport.fields._bases import BoundaryField


__all__ = ["BoundaryResidual"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class BoundaryResidual(BoundaryField):
    r"""L2 boundary-balance residual — the *residual* role leaf of
    :class:`~orpheus.transport.fields._bases.BoundaryField`.

    Parameters
    ----------
    values : NDArray
        Flat backing buffer, shape ``(space.layout.total_size,)``.
    space : TraceSpace
        The unified boundary
        :class:`~orpheus.numerics.spaces.trace_space.TraceSpace`
        (canonically ``mesh.trace``), carrying the per-geometry
        :class:`~orpheus.numerics.face_layout.FaceLayout`.
    mesh : SNMesh
        The SN phase-space carrier (cross-mesh-arithmetic guard).

    Notes
    -----
    A thin role leaf: all storage, validation, algebra, per-face access,
    and the ``zeros_on`` / ``from_face_arrays`` factories are
    inherited from :class:`BoundaryField`. The leaf carries no
    residual-specific behaviour beyond its class identity. All three
    boundary leaves share the SAME ``TraceSpace`` (``mesh.trace``), so it
    is the **class** gate (not the space gate) that rejects
    ``BoundaryResidual ± BoundaryFlux`` / ``BoundaryResidual ±
    BoundarySourceSink``. Same-class arithmetic is closed.

    See the module docstring for the already-computed matvec consumer and
    why the wiring (retype of the operator-output boundary) is the
    boundary half of the B.5 carve.
    """

    #: Dimensional identity (View-G, B.4): the trace is all-flux, so
    #: ``1/(cm²·s·sr)`` — :data:`~orpheus.numerics.units.ANGULAR_FLUX_UNITS`,
    #: shared with ``BoundaryFlux`` / ``BoundarySourceSink`` (same units,
    #: different role → class gate). Metadata, not the gate.
    UNITS: ClassVar[Unit] = ANGULAR_FLUX_UNITS
