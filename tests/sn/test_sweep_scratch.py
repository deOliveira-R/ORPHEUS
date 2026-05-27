r"""L0 — :class:`~orpheus.sn.sweep_scratch.SweepScratch` foundation tests.

Pins the sweep-private interior wavefront cache type that splits out of
:class:`~orpheus.sn.boundary_flux.BoundaryFlux` in Depth B step D-G.

Pins:

* :meth:`SweepScratch.for_sn_mesh` allocation policy (empty for 1-D,
  allocated for 2-D).
* :meth:`SweepScratch.ensure_2d_buffers` lazy-init (idempotent for
  2-D, raises for 1-D).
* Buffer shapes match the pre-D-G ``xmin_xmax_buf`` / ``ymin_ymax_buf``
  shapes (the post-D-G scratch carries the SAME working buffers; only
  the OWNERSHIP changes from BoundaryFlux to the sweep operator).

These are foundation tests — they pin the data-shape contract of the
new scratch type. The numerical-equivalence tests (sweep output matches
pre-D-G with SweepScratch in place) live in
``tests/sn/test_sweep_scratch_split.py``, added when the sweep is
rewired to consume SweepScratch.

References
----------

* Depth B plan §3.4 and §6 step D-G.
* `.claude/agent-memory/test-architect/dg_boundary_flux_pure_field_verification.md`.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.face_layout import FaceLayout
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.sweep_scratch import SweepScratch

from tests.sn._test_helpers import placeholder_materials


pytestmark = [pytest.mark.foundation]


# ───────────────────────────────────────────────────────────────────────
# Geometry fixtures (mirrored from tests/sn/test_typed_fields.py)
# ───────────────────────────────────────────────────────────────────────


def _slab_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _spherical_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _2d_mesh(nx: int = 3, ny: int = 2, ng: int = 2) -> SNMesh:
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, nx + 1),
        edges_y=np.linspace(0, 1, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


# ───────────────────────────────────────────────────────────────────────
# A. SweepScratch.for_sn_mesh allocation policy
# ───────────────────────────────────────────────────────────────────────


class TestForSnMeshAllocation:
    def test_1d_slab_returns_empty(self) -> None:
        """1-D slab: no 2-D wavefront buffer needed; both fields None."""
        m = _slab_mesh()
        scratch = SweepScratch.for_sn_mesh(m)
        assert scratch.psi_x_buf is None
        assert scratch.psi_y_buf is None

    def test_1d_sphere_returns_empty(self) -> None:
        """1-D sphere: same — no 2-D wavefront buffer."""
        m = _spherical_mesh()
        scratch = SweepScratch.for_sn_mesh(m)
        assert scratch.psi_x_buf is None
        assert scratch.psi_y_buf is None

    def test_2d_cartesian_allocates(self) -> None:
        """2-D Cartesian: psi_x_buf (N, ng, nx+1, ny), psi_y_buf (N, ng, nx, ny+1)."""
        m = _2d_mesh(nx=3, ny=2, ng=2)
        scratch = SweepScratch.for_sn_mesh(m)
        N = m.quad.N
        assert scratch.psi_x_buf is not None
        assert scratch.psi_y_buf is not None
        assert scratch.psi_x_buf.shape == (N, 2, 4, 2)   # (N, ng, nx+1, ny)
        assert scratch.psi_y_buf.shape == (N, 2, 3, 3)   # (N, ng, nx, ny+1)

    def test_2d_buffers_zero_initialised(self) -> None:
        """First-sweep contract: buffers start at zero."""
        m = _2d_mesh()
        scratch = SweepScratch.for_sn_mesh(m)
        np.testing.assert_array_equal(scratch.psi_x_buf, 0.0)
        np.testing.assert_array_equal(scratch.psi_y_buf, 0.0)


# ───────────────────────────────────────────────────────────────────────
# B. ensure_2d_buffers — lazy-init idempotency + 1-D rejection
# ───────────────────────────────────────────────────────────────────────


class TestEnsure2dBuffers:
    def test_ensure_allocates_on_empty_scratch(self) -> None:
        """Calling ensure on an empty 2-D scratch allocates the buffers."""
        m = _2d_mesh()
        scratch = SweepScratch()   # empty
        assert scratch.psi_x_buf is None
        scratch.ensure_2d_buffers(m)
        assert scratch.psi_x_buf is not None
        assert scratch.psi_y_buf is not None

    def test_ensure_idempotent_on_allocated(self) -> None:
        """Calling ensure twice preserves the SAME buffer instance."""
        m = _2d_mesh()
        scratch = SweepScratch.for_sn_mesh(m)
        first_x = scratch.psi_x_buf
        first_y = scratch.psi_y_buf
        scratch.ensure_2d_buffers(m)
        assert scratch.psi_x_buf is first_x
        assert scratch.psi_y_buf is first_y

    def test_ensure_raises_on_1d_mesh(self) -> None:
        """Asking for 2-D buffer allocation on a 1-D mesh is a programmer error."""
        m = _slab_mesh()
        scratch = SweepScratch()
        with pytest.raises(ValueError, match="1-D mesh"):
            scratch.ensure_2d_buffers(m)


# ───────────────────────────────────────────────────────────────────────
# C. Buffer shape matches pre-D-G BoundaryFlux 2-D conflated buffers
# ───────────────────────────────────────────────────────────────────────


class TestPreDGShapeMatch:
    """The carve preserves the buffer SHAPES; only the OWNERSHIP changes
    (BoundaryFlux pre-D-G -> SweepScratch post-D-G). The sweep's
    wavefront algorithm consumes the same shape (N, ng, nx+1, ny) for
    x-direction face fluxes and (N, ng, nx, ny+1) for y-direction —
    unchanged from pre-D-G's ``xmin_xmax_buf`` / ``ymin_ymax_buf``.
    """

    def test_shape_matches_legacy_xmin_xmax_buf(self) -> None:
        """psi_x_buf shape == pre-D-G xmin_xmax_buf shape."""
        m = _2d_mesh(nx=5, ny=3, ng=2)
        scratch = SweepScratch.for_sn_mesh(m)
        from orpheus.sn.boundary_flux import BoundaryFlux
        legacy = BoundaryFlux.zeros(m)
        assert scratch.psi_x_buf.shape == legacy.xmin_xmax_buf.shape

    def test_shape_matches_legacy_ymin_ymax_buf(self) -> None:
        """psi_y_buf shape == pre-D-G ymin_ymax_buf shape."""
        m = _2d_mesh(nx=5, ny=3, ng=2)
        scratch = SweepScratch.for_sn_mesh(m)
        from orpheus.sn.boundary_flux import BoundaryFlux
        legacy = BoundaryFlux.zeros(m)
        assert scratch.psi_y_buf.shape == legacy.ymin_ymax_buf.shape


# ───────────────────────────────────────────────────────────────────────
# D. SNMesh.boundary_face_layout (added in this commit)
# ───────────────────────────────────────────────────────────────────────


class TestBoundaryFaceLayout:
    """The mesh-side FaceLayout provider. Consumed by post-D-G
    :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux` to lay
    out its flat backing buffer."""

    def test_slab_layout_has_xmin_and_xmax(self) -> None:
        m = _slab_mesh(nx=5, ng=2)
        layout = m.boundary_face_layout
        N = m.quad.N
        assert isinstance(layout, FaceLayout)
        assert set(layout.faces) == {"xmin", "xmax"}
        assert layout.faces["xmin"].shape == (N, 2)
        assert layout.faces["xmax"].shape == (N, 2)
        assert layout.total_size == 2 * N * 2

    def test_sphere_layout_has_only_xmax(self) -> None:
        m = _spherical_mesh(nx=5, ng=2)
        layout = m.boundary_face_layout
        N = m.quad.N
        assert set(layout.faces) == {"xmax"}
        assert layout.faces["xmax"].shape == (N, 2)
        assert layout.total_size == N * 2

    def test_2d_layout_has_four_faces(self) -> None:
        m = _2d_mesh(nx=3, ny=2, ng=2)
        layout = m.boundary_face_layout
        N = m.quad.N
        assert set(layout.faces) == {"xmin", "xmax", "ymin", "ymax"}
        assert layout.faces["xmin"].shape == (N, 2, 2)
        assert layout.faces["xmax"].shape == (N, 2, 2)
        assert layout.faces["ymin"].shape == (N, 2, 3)
        assert layout.faces["ymax"].shape == (N, 2, 3)
        assert layout.total_size == N * 2 * (2 * 2 + 2 * 3)

    def test_2d_layout_excludes_interior_cells(self) -> None:
        """CRITICAL: post-D-G layout contains ONLY face slots.

        Pre-D-G ``xmin_xmax_buf`` had shape ``(N, ng, nx+1, ny)`` =
        ``N*ng*(nx+1)*ny`` cells (face slots [0, :] and [nx, :] PLUS
        interior slots [1:nx, :]). Post-D-G the FaceLayout for the 2-D
        boundary contains only [0, :] and [nx, :] per x-face — i.e.
        ``2 * N*ng*ny``, NOT ``N*ng*(nx+1)*ny``.

        Catches Rank 5 failure mode per the verification memo:
        FaceLayout botched to include interior cells (carve incomplete).
        """
        m = _2d_mesh(nx=3, ny=2, ng=2)
        N = m.quad.N
        layout = m.boundary_face_layout
        # Per-direction face-only total:
        x_faces = 2 * N * 2 * 2   # xmin + xmax, each (N, ng, ny) = N*2*2
        y_faces = 2 * N * 2 * 3   # ymin + ymax, each (N, ng, nx) = N*2*3
        assert layout.total_size == x_faces + y_faces
        # Pre-D-G conflated buffer total (NOT what we want):
        legacy_total = N * 2 * (3 + 1) * 2 + N * 2 * 3 * (2 + 1)
        assert layout.total_size < legacy_total, (
            f"FaceLayout post-D-G {layout.total_size} should be smaller "
            f"than pre-D-G conflated buffer total {legacy_total} — "
            f"interior cells live on SweepScratch."
        )

    def test_layout_is_idempotent_property(self) -> None:
        """Calling boundary_face_layout twice returns structurally equal layouts."""
        m = _slab_mesh()
        a = m.boundary_face_layout
        b = m.boundary_face_layout
        # Same set of face names, same shapes, same total_size.
        assert set(a.faces) == set(b.faces)
        for name in a.faces:
            assert a.faces[name].shape == b.faces[name].shape
            assert a.faces[name].offset == b.faces[name].offset
        assert a.total_size == b.total_size
