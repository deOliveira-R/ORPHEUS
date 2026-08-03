r"""Pure-Field :class:`AngularBoundaryFlux` — boundary trace on a flat layout.

L2 typed wrapper for the angular flux at the SN domain boundary,
inheriting :class:`~orpheus.transport.fields._bases.AngularBoundaryField`.
Storage is a SINGLE flat backing buffer; per-face access goes through a
:class:`~orpheus.numerics.face_layout.FaceLayout` descriptor (carried by
the :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` since A.5)
that maps face names to ``(offset, shape)`` slice views — no copies.

Migration status (Depth B step D-G → A.5 → B.1)
================================================

Pre-D-G ``orpheus.sn.boundary_flux.AngularBoundaryFlux`` was a MUTABLE bundle
with four per-geometry-conditional ndarray attributes (``xmin_face``,
``xmax_face``, ``xmin_xmax_buf``, ``ymin_ymax_buf``) — two of which
were always ``None`` for any given geometry, and the 2-D pair
conflated boundary face cells with interior wavefront cache.

Post-D-G this class:

* **Is a pure Field** (single flat ``values`` array + L1 ``space``
  + Field-inherited algebra) — and since **A.5** (View-G, #205/#201)
  that ``space`` IS the unified
  :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`, so the
  storage is genuinely ``values + space`` with NO extra geometry
  attribute. Cross-method generality: MoC's per-track-family "boundary
  faces" (thousands) work with the same flat-buffer arithmetic as SN's
  1-4 faces.
* **Is IMMUTABLE** (``frozen=True``). The pre-D-G mutable
  write-through pattern (sweep's persistent BC state) is replaced by
  sweep-side EPHEMERAL local arrays + functional reconstruction at
  iteration boundaries.
* **Carries face geometry via** :class:`~orpheus.numerics.face_layout.FaceLayout`
  **on the** :class:`AngularTraceSpace` **(A.5)**, not as a separate field
  attribute. ``mesh.angular_trace`` is the cached source; the underlying
  per-geometry descriptor is still
  :attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.boundary_face_layout` (1-D slab:
  ``xmin``, ``xmax``; 1-D curvilinear: ``xmax``; 2-D: ``xmin``,
  ``xmax``, ``ymin``, ``ymax``). The :attr:`layout` read-through
  property preserves the ``boundary.layout`` access surface.
* **Excludes interior wavefront cache.** The pre-D-G 2-D
  ``(N, ng, nx+1, ny)`` / ``(N, ng, nx, ny+1)`` conflated buffers
  split: face slots live here (``(N, ng, ny)`` × 2 + ``(N, ng, nx)``
  × 2 = much smaller); interior cells are sweep-local working
  storage, allocated inside the sweep and discarded with it (every
  sweep recomputes each interior cell from upstream propagation, so
  nothing there is authoritative state). A ``SweepScratch`` type was
  minted for them and deleted within the same step (2026-05-28,
  "Option I"): a sweep-private *persistent* type re-created the very
  boundary/interior conflation being dissolved, and the 1-D unified
  sweep had already shown the right pattern — an ephemeral local
  ``psi_face_chain``. **This trace is the sole persistent face-state
  carrier.**

**B.1 (field vocabulary):** every member that was generic to a boundary
trace field — the ``mesh`` field, the cross-mesh guard, the AngularTraceSpace
contract, the :attr:`layout` property, :meth:`face_view`, and the
:meth:`zeros_on` / :meth:`from_face_arrays` factories — moved up
to the :class:`~orpheus.transport.fields._bases.AngularBoundaryField` storage
base. ``AngularBoundaryFlux`` is the *flux* role leaf; ``AngularBoundarySourceSink``
(``orpheus.transport.source_sinks``) and ``AngularBoundaryResidual``
(``orpheus.transport.residuals``) joined it under ``AngularBoundaryField`` in B.3.

References
==========

* Depth B plan §3.4 (Option Ω flat-buffer storage) and §6 step D-G;
  ``field_role_typing_view_g.md`` A.5 (AngularTraceSpace re-home) + B.1
  (storage-base ABCs).
* `.claude/agent-memory/test-architect/dg_boundary_flux_pure_field_verification.md`.
* `.claude/agent-memory/explorer/dg_boundary_flux_consumer_audit.md`.
* ``coding-elegance`` Pattern 1 (read-as-the-math via dunder),
  Pattern 2 (single source of truth — AngularBoundaryField), Pattern 4
  (illegal states unrepresentable — immutability).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.units import ANGULAR_FLUX_UNITS, Unit
from orpheus.transport.fields._bases import AngularBoundaryField
from orpheus.transport.fields._flux_role import FluxRole


__all__ = ["AngularBoundaryFlux"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class AngularBoundaryFlux(FluxRole, AngularBoundaryField):
    r"""L2 boundary trace flux — the *flux* role leaf of
    :class:`~orpheus.transport.fields._bases.AngularBoundaryField`.

    Parameters
    ----------
    values : NDArray
        Flat backing buffer, shape ``(space.layout.total_size,)``.
    space : AngularTraceSpace
        The unified boundary
        :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` — the
        L1 space anchor (Euclidean inner product) that also carries the
        per-geometry :class:`~orpheus.numerics.face_layout.FaceLayout`.
        Canonically the mesh's cached
        :attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.angular_trace`.
    mesh : SNMesh
        The SN phase-space carrier (the cross-mesh-arithmetic guard).

    Notes
    -----
    All storage, validation, algebra, per-face access, and construction
    machinery is inherited from
    :class:`~orpheus.transport.fields._bases.AngularBoundaryField` (B.1). The
    one flux-specific behaviour this leaf carries is :meth:`net_current`
    — a *current* is an angle-weighted flux moment, so the reduction is
    flux-role semantics (mirroring
    :meth:`~orpheus.transport.fields.scalar_boundary_flux.ScalarBoundaryFlux.net_current`
    on the scalar sibling); the class identity otherwise remains what
    Field's Layer-1 gate uses to keep boundary flux, source, and
    residual arithmetic from silently mixing. Build
    via :meth:`~orpheus.transport.fields._bases.AngularBoundaryField.zeros_on`
    / :meth:`~orpheus.transport.fields._bases.AngularBoundaryField.from_face_arrays`.
    """

    #: Dimensional identity (View-G, B.4): the boundary trace stores flux
    #: values, so ``1/(cm²·s·sr)`` —
    #: :data:`~orpheus.numerics.units.ANGULAR_FLUX_UNITS`, shared with
    #: ``AngularFlux`` and the sibling boundary leaves (the boundary is
    #: all-flux). Metadata, not the arithmetic gate.
    UNITS: ClassVar[Unit] = ANGULAR_FLUX_UNITS

    def net_current(self, face: str) -> NDArray:
        r"""Net OUTWARD current on ``face`` — shape ``(ng, *face_spatial)``, a copy.

        The angular-to-scalar reduction

        .. math::

            J_g \;=\; \sum_m (\Omega_m\cdot\hat n_f)\, w_m\, \psi_{m,g}

        — the angular sibling of
        :meth:`ScalarBoundaryFlux.net_current
        <orpheus.transport.fields.scalar_boundary_flux.ScalarBoundaryFlux.net_current>`
        (``J = J⁺ − J⁻``), under the same face-local outward-normal
        convention: positive = leakage out of the domain.

        Spelled through the trace space's own atoms — the signed
        projection table
        (:attr:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.omega_dot_n`)
        and the :math:`|\Omega\cdot\hat n|\odot w` partial-current metric
        (``sign(Ω·n̂) · metric = Ω·n̂ · w``) — so no consumer re-derives
        the cosine weighting (the #291 leakage term and the future DSA
        boundary restriction both read THIS reduction). Tangential
        ordinates carry zero weight and drop out.
        """
        view = self.face_view(face)  # validates the face name
        space = self.space
        row = space.omega_dot_n[space.face_names.index(face)]
        slot = space.layout.faces[face]
        metric = space.partial_current_metric[
            slot.offset : slot.offset + slot.flat_size
        ].reshape(slot.shape)
        sign_axis0 = np.sign(row).reshape(
            (row.shape[0],) + (1,) * (len(slot.shape) - 1)
        )
        return np.einsum("m...,m...->...", sign_axis0 * metric, view)
