"""Foundation tests for the Issue #197 PR-TYPED-2 typed field types.

Pins the structural contract of :class:`AngularFlux`,
:class:`ScalarFlux`, and :class:`BoundaryFlux` — shape validation,
mesh-binding, dunder arithmetic, reductions, and the
:class:`SNMesh` factory methods.

These are foundation tests — they verify software invariants
(constructor shape check, frozen-ness, dunder algebra) rather than
physics claims, so they carry ``@pytest.mark.foundation`` per the V&V
harness convention.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.sn.angular_flux import AngularFlux
from orpheus.sn.boundary_flux import BoundaryFlux
from orpheus.sn.geometry import SNMesh
from orpheus.sn.quadrature import GaussLegendre1D, LevelSymmetricSN
from orpheus.sn.scalar_flux import ScalarFlux

from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation


# ── Fixtures ─────────────────────────────────────────────────────────


def _slab_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
    """Build a small slab :class:`SNMesh` for unit testing."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = GaussLegendre1D.create(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _spherical_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
    """Build a small spherical :class:`SNMesh`."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = GaussLegendre1D.create(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _2d_mesh(nx: int = 3, ny: int = 3, ng: int = 1) -> SNMesh:
    """Build a small 2-D Cartesian :class:`SNMesh`."""
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, nx + 1),
        edges_y=np.linspace(0, 1, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
    )
    quad = LevelSymmetricSN.create(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


# ════════════════════════════════════════════════════════════════════
# ScalarFlux
# ════════════════════════════════════════════════════════════════════


class TestScalarFlux:
    def test_construct_from_factory(self) -> None:
        m = _slab_mesh()
        phi = m.zeros_scalar_flux()
        assert phi.values.shape == (m.ng, m.nx, m.ny)
        assert np.all(phi.values == 0.0)
        assert phi.mesh is m

    def test_shape_validation_rejects_wrong_shape(self) -> None:
        m = _slab_mesh()
        with pytest.raises(ValueError, match="ScalarFlux expects shape"):
            ScalarFlux(np.zeros((m.ng, m.nx + 1, m.ny)), m)

    def test_addition_returns_typed(self) -> None:
        m = _slab_mesh()
        a = ScalarFlux(np.ones((m.ng, m.nx, m.ny)), m)
        b = ScalarFlux(2.0 * np.ones((m.ng, m.nx, m.ny)), m)
        c = a + b
        assert isinstance(c, ScalarFlux)
        np.testing.assert_array_equal(c.values, 3.0 * np.ones((m.ng, m.nx, m.ny)))

    def test_subtraction(self) -> None:
        m = _slab_mesh()
        a = ScalarFlux(np.full((m.ng, m.nx, m.ny), 5.0), m)
        b = ScalarFlux(np.full((m.ng, m.nx, m.ny), 2.0), m)
        c = a - b
        assert isinstance(c, ScalarFlux)
        np.testing.assert_array_equal(c.values, 3.0 * np.ones((m.ng, m.nx, m.ny)))

    def test_scalar_multiplication_left_and_right(self) -> None:
        m = _slab_mesh()
        a = ScalarFlux(np.ones((m.ng, m.nx, m.ny)), m)
        left = 3.0 * a
        right = a * 3.0
        np.testing.assert_array_equal(left.values, right.values)
        np.testing.assert_array_equal(left.values, 3.0 * np.ones((m.ng, m.nx, m.ny)))

    def test_arithmetic_rejects_different_mesh(self) -> None:
        a = _slab_mesh().zeros_scalar_flux()
        b = _slab_mesh().zeros_scalar_flux()  # different SNMesh instance
        with pytest.raises(ValueError, match="distinct SNMesh"):
            _ = a + b

    def test_arithmetic_rejects_other_type(self) -> None:
        a = _slab_mesh().zeros_scalar_flux()
        with pytest.raises(TypeError):
            _ = a + np.ones_like(a.values)  # bare ndarray rejected

    def test_at_group_returns_view(self) -> None:
        m = _slab_mesh()
        vals = np.arange(m.ng * m.nx * m.ny, dtype=float).reshape(
            m.ng, m.nx, m.ny,
        )
        phi = ScalarFlux(vals, m)
        slice_g0 = phi.at_group(0)
        assert slice_g0.shape == (m.nx, m.ny)
        np.testing.assert_array_equal(slice_g0, vals[0])

    def test_linf(self) -> None:
        m = _slab_mesh()
        vals = np.full((m.ng, m.nx, m.ny), -7.0)
        phi = ScalarFlux(vals, m)
        assert phi.linf() == 7.0

    def test_copy_owns_buffer(self) -> None:
        m = _slab_mesh()
        a = ScalarFlux(np.ones((m.ng, m.nx, m.ny)), m)
        b = a.copy()
        # Mutating the copy's array does not affect ``a``.
        b.values[...] = 99.0
        np.testing.assert_array_equal(a.values, np.ones((m.ng, m.nx, m.ny)))


# ════════════════════════════════════════════════════════════════════
# AngularFlux
# ════════════════════════════════════════════════════════════════════


class TestAngularFlux:
    def test_construct_from_factory(self) -> None:
        m = _slab_mesh()
        psi = m.zeros_angular_flux()
        assert psi.values.shape == (m.quad.N, m.ng, m.nx, m.ny)
        assert np.all(psi.values == 0.0)
        assert psi.mesh is m

    def test_shape_validation_rejects_wrong_shape(self) -> None:
        m = _slab_mesh()
        with pytest.raises(ValueError, match="AngularFlux expects shape"):
            AngularFlux(
                np.zeros((m.quad.N + 1, m.ng, m.nx, m.ny)), m,
            )

    def test_integrate_angular_returns_typed(self) -> None:
        m = _slab_mesh()
        # Set ψ = c (constant in space, ordinate, group).  Then
        # φ = Σ_n w_n · c = c · Σ_n w_n.
        c = 2.5
        psi = AngularFlux(
            c * np.ones((m.quad.N, m.ng, m.nx, m.ny)), m,
        )
        phi = psi.integrate_angular()
        assert isinstance(phi, ScalarFlux)
        expected = c * m.quad.weights.sum()
        np.testing.assert_allclose(
            phi.values, expected * np.ones((m.ng, m.nx, m.ny)),
            rtol=1e-14,
        )

    def test_addition_returns_typed(self) -> None:
        m = _slab_mesh()
        a = AngularFlux(np.ones((m.quad.N, m.ng, m.nx, m.ny)), m)
        b = AngularFlux(
            2.0 * np.ones((m.quad.N, m.ng, m.nx, m.ny)), m,
        )
        c = a + b
        assert isinstance(c, AngularFlux)
        np.testing.assert_array_equal(
            c.values, 3.0 * np.ones((m.quad.N, m.ng, m.nx, m.ny)),
        )

    def test_at_ordinate_returns_view(self) -> None:
        m = _slab_mesh()
        vals = np.arange(
            m.quad.N * m.ng * m.nx * m.ny, dtype=float,
        ).reshape(m.quad.N, m.ng, m.nx, m.ny)
        psi = AngularFlux(vals, m)
        s = psi.at_ordinate(2)
        assert s.shape == (m.ng, m.nx, m.ny)
        np.testing.assert_array_equal(s, vals[2])


# ════════════════════════════════════════════════════════════════════
# BoundaryFlux
# ════════════════════════════════════════════════════════════════════


class TestBoundaryFluxFactoryShapes:
    """``zeros`` factory must allocate the right buffers per geometry."""

    def test_slab_allocates_xmin_and_xmax_faces(self) -> None:
        m = _slab_mesh()
        bf = m.zeros_boundary_flux()
        N, ng = m.quad.N, m.ng
        assert bf.xmin_face is not None
        assert bf.xmin_face.shape == (N, ng)
        assert bf.xmax_face is not None
        assert bf.xmax_face.shape == (N, ng)
        # 2D buffers stay None for 1-D meshes.
        assert bf.xmin_xmax_buf is None
        assert bf.ymin_ymax_buf is None

    def test_spherical_allocates_only_xmax_face(self) -> None:
        m = _spherical_mesh()
        bf = m.zeros_boundary_flux()
        N, ng = m.quad.N, m.ng
        # Outer radial face only — the pole is a regularity condition,
        # not a BC face.
        assert bf.xmin_face is None
        assert bf.xmax_face is not None
        assert bf.xmax_face.shape == (N, ng)

    def test_2d_allocates_persistent_buffers(self) -> None:
        m = _2d_mesh()
        bf = m.zeros_boundary_flux()
        N, ng, nx, ny = m.quad.N, m.ng, m.nx, m.ny
        assert bf.xmin_xmax_buf is not None
        assert bf.xmin_xmax_buf.shape == (N, ng, nx + 1, ny)
        assert bf.ymin_ymax_buf is not None
        assert bf.ymin_ymax_buf.shape == (N, ng, nx, ny + 1)


class TestBoundaryFluxFaceAccessors:
    """Per-face accessors must return correctly-shaped views."""

    def test_2d_xmin_xmax_slice_views(self) -> None:
        m = _2d_mesh()
        bf = m.zeros_boundary_flux()
        N, ng, nx, ny = m.quad.N, m.ng, m.nx, m.ny

        # xmin → slot [:, :, 0, :]; xmax → slot [:, :, nx, :].
        bf.xmin_xmax_buf[:, :, 0, :] = 1.0
        bf.xmin_xmax_buf[:, :, nx, :] = 2.0
        np.testing.assert_array_equal(
            bf.xmin, np.ones((N, ng, ny)),
        )
        np.testing.assert_array_equal(
            bf.xmax, 2.0 * np.ones((N, ng, ny)),
        )

    def test_2d_ymin_ymax_slice_views(self) -> None:
        m = _2d_mesh()
        bf = m.zeros_boundary_flux()
        N, ng, nx, ny = m.quad.N, m.ng, m.nx, m.ny

        bf.ymin_ymax_buf[:, :, :, 0] = 3.0
        bf.ymin_ymax_buf[:, :, :, ny] = 4.0
        np.testing.assert_array_equal(
            bf.ymin, 3.0 * np.ones((N, ng, nx)),
        )
        np.testing.assert_array_equal(
            bf.ymax, 4.0 * np.ones((N, ng, nx)),
        )

    def test_1d_xmin_xmax_via_face_buffers(self) -> None:
        m = _slab_mesh()
        bf = m.zeros_boundary_flux()
        N, ng = m.quad.N, m.ng

        bf.xmin_face[:, :] = 5.0
        bf.xmax_face[:, :] = 6.0
        np.testing.assert_array_equal(bf.xmin, 5.0 * np.ones((N, ng)))
        np.testing.assert_array_equal(bf.xmax, 6.0 * np.ones((N, ng)))

    def test_1d_has_no_y_faces(self) -> None:
        m = _slab_mesh()
        bf = m.zeros_boundary_flux()
        assert bf.ymin is None
        assert bf.ymax is None
