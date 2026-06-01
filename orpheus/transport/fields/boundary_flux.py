r"""Pure-Field :class:`BoundaryFlux` — boundary trace on a flat layout.

L2 typed wrapper for the angular flux at the SN domain boundary,
inheriting :class:`~orpheus.transport.fields._bases.BoundaryField`.
Storage is a SINGLE flat backing buffer; per-face access goes through a
:class:`~orpheus.numerics.face_layout.FaceLayout` descriptor (carried by
the :class:`~orpheus.numerics.spaces.trace_space.TraceSpace` since A.5)
that maps face names to ``(offset, shape)`` slice views — no copies.

Migration status (Depth B step D-G → A.5 → B.1)
================================================

Pre-D-G ``orpheus.sn.boundary_flux.BoundaryFlux`` was a MUTABLE bundle
with four per-geometry-conditional ndarray attributes (``xmin_face``,
``xmax_face``, ``xmin_xmax_buf``, ``ymin_ymax_buf``) — two of which
were always ``None`` for any given geometry, and the 2-D pair
conflated boundary face cells with interior wavefront cache.

Post-D-G this class:

* **Is a pure Field** (single flat ``values`` array + L1 ``space``
  + Field-inherited algebra) — and since **A.5** (View-G, #205/#201)
  that ``space`` IS the unified
  :class:`~orpheus.numerics.spaces.trace_space.TraceSpace`, so the
  storage is genuinely ``values + space`` with NO extra geometry
  attribute. Cross-method generality: MoC's per-track-family "boundary
  faces" (thousands) work with the same flat-buffer arithmetic as SN's
  1-4 faces.
* **Is IMMUTABLE** (``frozen=True``). The pre-D-G mutable
  write-through pattern (sweep's persistent BC state) is replaced by
  sweep-side :class:`~orpheus.sn.sweep_scratch.SweepScratch` +
  functional reconstruction at iteration boundaries.
* **Carries face geometry via** :class:`~orpheus.numerics.face_layout.FaceLayout`
  **on the** :class:`TraceSpace` **(A.5)**, not as a separate field
  attribute. ``mesh.trace`` is the cached source; the underlying
  per-geometry descriptor is still
  :attr:`~orpheus.sn.geometry.SNMesh.boundary_face_layout` (1-D slab:
  ``xmin``, ``xmax``; 1-D curvilinear: ``xmax``; 2-D: ``xmin``,
  ``xmax``, ``ymin``, ``ymax``). The :attr:`layout` read-through
  property preserves the ``boundary.layout`` access surface.
* **Excludes interior wavefront cache.** The pre-D-G 2-D
  ``(N, ng, nx+1, ny)`` / ``(N, ng, nx, ny+1)`` conflated buffers
  split: face slots live here (``(N, ng, ny)`` × 2 + ``(N, ng, nx)``
  × 2 = much smaller); interior cells live on
  :class:`~orpheus.sn.sweep_scratch.SweepScratch` (sweep-private,
  owned by the sweep operator across iterations).

**B.1 (field vocabulary):** every member that was generic to a boundary
trace field — the ``mesh`` field, the cross-mesh guard, the TraceSpace
contract, the :attr:`layout` property, :meth:`face_view`, and the
:meth:`zeros_for_sn_mesh` / :meth:`from_face_arrays` factories — moved up
to the :class:`~orpheus.transport.fields._bases.BoundaryField` storage
base. ``BoundaryFlux`` is now the *flux* role leaf; ``BoundarySource`` /
``BoundaryResidual`` will join it under ``BoundaryField`` in B.3.

References
==========

* Depth B plan §3.4 (Option Ω flat-buffer storage) and §6 step D-G;
  ``field_role_typing_view_g.md`` A.5 (TraceSpace re-home) + B.1
  (storage-base ABCs).
* `.claude/agent-memory/test-architect/dg_boundary_flux_pure_field_verification.md`.
* `.claude/agent-memory/explorer/dg_boundary_flux_consumer_audit.md`.
* ``coding-elegance`` Pattern 1 (read-as-the-math via dunder),
  Pattern 2 (single source of truth — BoundaryField), Pattern 4
  (illegal states unrepresentable — immutability).
"""

from __future__ import annotations

from dataclasses import dataclass

from orpheus.transport.fields._bases import BoundaryField


__all__ = ["BoundaryFlux"]


@dataclass(frozen=True, eq=False, kw_only=True)
class BoundaryFlux(BoundaryField):
    r"""L2 boundary trace flux — the *flux* role leaf of
    :class:`~orpheus.transport.fields._bases.BoundaryField`.

    Parameters
    ----------
    values : NDArray
        Flat backing buffer, shape ``(space.layout.total_size,)``.
    space : TraceSpace
        The unified boundary
        :class:`~orpheus.numerics.spaces.trace_space.TraceSpace` — the
        L1 space anchor (Euclidean inner product) that also carries the
        per-geometry :class:`~orpheus.numerics.face_layout.FaceLayout`.
        Canonically the mesh's cached
        :attr:`~orpheus.sn.geometry.SNMesh.trace`.
    mesh : SNMesh
        The SN phase-space carrier (the cross-mesh-arithmetic guard).

    Notes
    -----
    All storage, validation, algebra, per-face access, and construction
    machinery is inherited from
    :class:`~orpheus.transport.fields._bases.BoundaryField` (B.1). This
    leaf carries no flux-specific behaviour beyond its class identity —
    which is exactly what Field's Layer-1 gate uses to keep boundary
    flux, source, and residual arithmetic from silently mixing. Build
    via :meth:`~orpheus.transport.fields._bases.BoundaryField.zeros_for_sn_mesh`
    / :meth:`~orpheus.transport.fields._bases.BoundaryField.from_face_arrays`.
    """
