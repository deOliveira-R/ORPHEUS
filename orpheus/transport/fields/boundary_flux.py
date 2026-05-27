r"""Pure-Field :class:`BoundaryFlux` — boundary trace on a flat layout.

L2 typed wrapper for the angular flux at the SN domain boundary,
inheriting :class:`~orpheus.numerics.field.Field`. Storage is a SINGLE
flat backing buffer; per-face access goes through a
:class:`~orpheus.numerics.face_layout.FaceLayout` descriptor that maps
face names to ``(offset, shape)`` slice views — no copies.

Migration status (Depth B step D-G)
====================================

Pre-D-G ``orpheus.sn.boundary_flux.BoundaryFlux`` was a MUTABLE bundle
with four per-geometry-conditional ndarray attributes (``xmin_face``,
``xmax_face``, ``xmin_xmax_buf``, ``ymin_ymax_buf``) — two of which
were always ``None`` for any given geometry, and the 2-D pair
conflated boundary face cells with interior wavefront cache.

Post-D-G this class:

* **Is a pure Field** (single flat ``values`` array + L1 ``space``
  + Field-inherited algebra). Cross-method generality: MoC's
  per-track-family "boundary faces" (thousands) work with the same
  flat-buffer arithmetic as SN's 1-4 faces.
* **Is IMMUTABLE** (``frozen=True``). The pre-D-G mutable
  write-through pattern (sweep's persistent BC state) is replaced by
  sweep-side :class:`~orpheus.sn.sweep_scratch.SweepScratch` +
  functional reconstruction at iteration boundaries.
* **Carries face geometry via** :class:`~orpheus.numerics.face_layout.FaceLayout`
  on the mesh. The mesh's :attr:`~orpheus.sn.geometry.SNMesh.boundary_face_layout`
  property provides the per-geometry descriptor (1-D slab: ``xmin``,
  ``xmax``; 1-D curvilinear: ``xmax``; 2-D: ``xmin``, ``xmax``,
  ``ymin``, ``ymax``).
* **Excludes interior wavefront cache.** The pre-D-G 2-D
  ``(N, ng, nx+1, ny)`` / ``(N, ng, nx, ny+1)`` conflated buffers
  split: face slots live here (``(N, ng, ny)`` × 2 + ``(N, ng, nx)``
  × 2 = much smaller); interior cells live on
  :class:`~orpheus.sn.sweep_scratch.SweepScratch` (sweep-private,
  owned by the sweep operator across iterations).

References
==========

* Depth B plan §3.4 (Option Ω flat-buffer storage) and §6 step D-G.
* `.claude/agent-memory/test-architect/dg_boundary_flux_pure_field_verification.md`.
* `.claude/agent-memory/explorer/dg_boundary_flux_consumer_audit.md`.
* ``coding-elegance`` Pattern 1 (read-as-the-math via dunder),
  Pattern 4 (illegal states unrepresentable — immutability).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Mapping

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.face_layout import FaceLayout
from orpheus.numerics.field import Field
from orpheus.numerics.space import FunctionSpace

if TYPE_CHECKING:
    from orpheus.sn.geometry import SNMesh


__all__ = ["BoundaryFlux"]


@dataclass(frozen=True, eq=False, kw_only=True)
class BoundaryFlux(Field):
    r"""L2 boundary trace flux: pure :class:`Field` with flat-layout storage.

    Parameters
    ----------
    values : NDArray
        Flat backing buffer, shape ``(layout.total_size,)``. Field
        algebra (``__add__``, scalar ``*``, etc.) operates on this
        buffer in ONE numpy call.
    space : FunctionSpace
        The L1 space anchor (shape ``(layout.total_size,)``). Used
        by Field's dunder algebra for class+space matching.
    layout : FaceLayout
        Per-geometry face descriptor. Provides ``faces[name].slice_view``
        access into the flat buffer. Typically obtained from
        :attr:`~orpheus.sn.geometry.SNMesh.boundary_face_layout`.
    mesh : SNMesh
        The SN phase-space carrier. Mesh-binding is enforced in
        :meth:`_check_partner` — arithmetic across distinct mesh
        instances is forbidden even when the spaces match.

    Notes
    -----
    Algebra is inherited from :class:`~orpheus.numerics.field.Field`.
    Same-class arithmetic is closed. Cross-class arithmetic is rejected
    by Field's Layer 1 class-identity gate. Cross-mesh arithmetic is
    rejected by the mesh-binding override (Rank 4 of the D-G
    verification memo).
    """

    layout: FaceLayout
    mesh: "SNMesh"

    # ── Construction validation ──────────────────────────────────────

    def __post_init__(self) -> None:
        super().__post_init__()
        expected = (self.layout.total_size,)
        if self.values.shape != expected:
            raise ValueError(
                f"BoundaryFlux: values.shape {self.values.shape!r} does "
                f"not match (layout.total_size,) = {expected!r}"
            )
        if self.space.shape != expected:
            raise ValueError(
                f"BoundaryFlux: space.shape {self.space.shape!r} does "
                f"not match (layout.total_size,) = {expected!r}"
            )

    # ── Algebra extensions (over Field) ──────────────────────────────

    def _check_partner(self, other: object) -> None:
        super()._check_partner(other)
        # Mesh-binding override — two BoundaryFluxes can share a space
        # ("sn_boundary_flat", same total_size) but differ in mesh
        # identity. Rank 4 failure mode per the verification memo.
        if self.mesh is not other.mesh:  # type: ignore[attr-defined]
            raise ValueError(
                "BoundaryFlux arithmetic across distinct SNMesh instances "
                "is forbidden — the field is mesh-bound."
            )
        if self.layout is not other.layout:  # type: ignore[attr-defined]
            # Layouts may match structurally without sharing identity
            # (e.g., re-derived from the same mesh). Fall back to
            # structural equality — same set of face names + same
            # per-face shapes + same offsets.
            if self.layout != other.layout:  # type: ignore[attr-defined]
                raise ValueError(
                    "BoundaryFlux layout mismatch — operands have "
                    "structurally distinct FaceLayouts."
                )

    # ── Per-face access (slice views into the flat buffer) ───────────

    def face_view(self, name: str) -> NDArray:
        r"""Return a per-face slice view into the flat backing buffer.

        The returned ndarray shares memory with :attr:`values` — writes
        propagate to the backing buffer. This is the post-D-G analog
        of the pre-D-G ``bf.xmin``, ``bf.xmax``, ``bf.ymin``,
        ``bf.ymax`` accessor properties.

        Parameters
        ----------
        name : str
            Face name (must be a key in :attr:`layout.faces`).

        Returns
        -------
        np.ndarray
            Shaped view (per-face shape from :class:`FaceSlot.shape`).
            Memory-shared with :attr:`values`.

        Raises
        ------
        KeyError
            If ``name`` is not a face in this layout.
        """
        if name not in self.layout.faces:
            raise KeyError(
                f"BoundaryFlux: no face named {name!r} in layout; "
                f"available faces: {list(self.layout.faces)!r}"
            )
        return self.layout.faces[name].slice_view(self.values)

    @property
    def face_views(self) -> Mapping[str, NDArray]:
        r"""Mapping ``{face_name: slice_view}`` for every face in the layout.

        Equivalent to ``{name: self.face_view(name) for name in
        self.layout.faces}`` — convenience for bulk iteration. All
        views memory-shared with :attr:`values`.
        """
        return {name: self.face_view(name) for name in self.layout.faces}

    # ── Construction factories ───────────────────────────────────────

    @classmethod
    def zeros_for_sn_mesh(cls, mesh: "SNMesh") -> "BoundaryFlux":
        r"""Construct a zero :class:`BoundaryFlux` sized to ``mesh``.

        The post-D-G analog of pre-D-G :meth:`BoundaryFlux.zeros`. Uses
        ``mesh.boundary_face_layout`` (added in D-G to
        :class:`~orpheus.sn.geometry.SNMesh`) to determine the per-
        geometry flat layout, then constructs an all-zero flat buffer.

        Parameters
        ----------
        mesh : SNMesh
            The SN phase-space carrier.

        Returns
        -------
        BoundaryFlux
            All-zero boundary flux on the mesh's flat boundary layout.
        """
        layout = mesh.boundary_face_layout
        space = FunctionSpace(
            name="sn_boundary_flat",
            shape=(layout.total_size,),
        )
        return cls(
            values=np.zeros(layout.total_size),
            space=space,
            layout=layout,
            mesh=mesh,
        )

    @classmethod
    def from_face_arrays(
        cls,
        mesh: "SNMesh",
        face_arrays: Mapping[str, NDArray],
    ) -> "BoundaryFlux":
        r"""Construct from per-face ndarrays, packing into the flat layout.

        The post-D-G "builder" pattern for sweep / matvec sites that
        compute per-face values (one ndarray per face) and need to
        assemble them into a :class:`BoundaryFlux`. Replaces the pre-D-G
        per-attribute mutation pattern (``bf.xmin_face = …``,
        ``bf.xmax_face = …``).

        Parameters
        ----------
        mesh : SNMesh
            Phase-space carrier; provides :attr:`boundary_face_layout`.
        face_arrays : Mapping[str, NDArray]
            ``{face_name: per_face_values}`` for EVERY face in the mesh's
            layout. Faces missing from this mapping raise; faces in the
            mapping but absent from the layout also raise.

        Returns
        -------
        BoundaryFlux
            Flat-packed boundary flux.

        Raises
        ------
        ValueError
            If ``face_arrays`` keys differ from the mesh's layout
            faces, or any per-face ndarray's shape mismatches the
            layout's expected per-face shape.
        """
        layout = mesh.boundary_face_layout
        provided = set(face_arrays.keys())
        expected = set(layout.faces.keys())
        if provided != expected:
            raise ValueError(
                f"BoundaryFlux.from_face_arrays: face_arrays keys "
                f"{sorted(provided)!r} do not match mesh's layout faces "
                f"{sorted(expected)!r}"
            )

        flat = np.zeros(layout.total_size)
        for name, slot in layout.faces.items():
            arr = face_arrays[name]
            if arr.shape != slot.shape:
                raise ValueError(
                    f"BoundaryFlux.from_face_arrays: face {name!r} array "
                    f"shape {arr.shape!r} does not match layout slot "
                    f"shape {slot.shape!r}"
                )
            slot.slice_view(flat)[:] = arr

        space = FunctionSpace(
            name="sn_boundary_flat",
            shape=(layout.total_size,),
        )
        return cls(values=flat, space=space, layout=layout, mesh=mesh)
