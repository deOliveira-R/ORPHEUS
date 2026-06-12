r"""L0 — :attr:`SNMesh.boundary_face_layout` foundation tests.

Pins the per-geometry :class:`~orpheus.numerics.face_layout.FaceLayout`
provider that supplies the flat layout for post-D-G pure-Field
:class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`.

Geometry coverage:

* 1-D slab — two faces (``xmin``, ``xmax``), each ``(N, ng)``.
* 1-D curvilinear (sphere, cylinder) — single outer face (``xmax``).
* 2-D Cartesian — four faces (``xmin``, ``xmax``, ``ymin``, ``ymax``).

CRITICAL: the layout excludes interior wavefront cache cells (the
pre-D-G ``xmin_xmax_buf`` / ``ymin_ymax_buf`` interior positions
``[1:nx, :]`` and ``[:, 1:ny]`` are NOT here). Interior cells live as
ephemeral local arrays inside ``_sweep_jacobi``, matching the
1-D blueprint pattern (``psi_face_chain_scan`` is local in
``_sweep_1d_unified``).

References
----------

* Depth B plan §3.4 (Option Ω flat-buffer storage) and §6 step D-G.
* `.claude/agent-memory/test-architect/dg_boundary_flux_pure_field_verification.md`
  §3.1 (Rank 5 failure mode — "interior cells in layout").
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.face_layout import FaceLayout
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh

from tests.sn._test_helpers import placeholder_materials


pytestmark = [pytest.mark.foundation]


# ───────────────────────────────────────────────────────────────────────
# Geometry fixtures
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
# Per-geometry layouts
# ───────────────────────────────────────────────────────────────────────


class TestBoundaryFaceLayout:
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
        ``N*ng*(nx+1)*ny`` cells (face slots ``[0, :]`` and ``[nx, :]``
        PLUS interior slots ``[1:nx, :]``). Post-D-G the boundary
        face layout contains only ``[0, :]`` and ``[nx, :]`` per
        x-face — i.e. ``2 * N*ng*ny``, NOT ``N*ng*(nx+1)*ny``.

        Catches Rank 5 failure mode per the verification memo:
        layout botched to include interior cells (carve incomplete).
        Per Option I (2026-05-28), interior cells live as ephemeral
        local arrays in ``_sweep_jacobi`` — not in BoundaryFlux,
        not in any sweep-private type.
        """
        m = _2d_mesh(nx=3, ny=2, ng=2)
        N = m.quad.N
        layout = m.boundary_face_layout
        x_faces = 2 * N * 2 * 2   # xmin + xmax, each (N, ng, ny)
        y_faces = 2 * N * 2 * 3   # ymin + ymax, each (N, ng, nx)
        assert layout.total_size == x_faces + y_faces
        legacy_total = N * 2 * (3 + 1) * 2 + N * 2 * 3 * (2 + 1)
        assert layout.total_size < legacy_total, (
            f"FaceLayout post-D-G {layout.total_size} should be smaller "
            f"than pre-D-G conflated buffer total {legacy_total} — "
            f"interior cells are ephemeral local arrays in the sweep."
        )

    def test_layout_is_idempotent_property(self) -> None:
        m = _slab_mesh()
        a = m.boundary_face_layout
        b = m.boundary_face_layout
        assert set(a.faces) == set(b.faces)
        for name in a.faces:
            assert a.faces[name].shape == b.faces[name].shape
            assert a.faces[name].offset == b.faces[name].offset
        assert a.total_size == b.total_size
