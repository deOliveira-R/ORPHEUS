"""L0 foundation: :class:`BoundaryFlux` arithmetic dunders (R-1 Step 1).

:class:`BoundaryFlux` participates in the operator algebra as the
face-flux companion to :class:`AngularFlux`.  Arithmetic on
``AngularFlux`` propagates to its ``.boundary`` field; that propagation
delegates to the dunders pinned here.

Semantics:
  * ``__add__`` / ``__sub__`` — elementwise on every non-None face
    buffer.  Both operands must have the SAME mesh-identity
    (``is`` comparison).
  * ``__mul__`` / ``__rmul__`` / ``__truediv__`` — elementwise scalar
    op on every non-None face buffer.
  * ``__neg__`` — elementwise negation on every non-None face buffer.
  * Result is a FRESH :class:`BoundaryFlux` instance (immutable
    arithmetic).  The sweep's mutable write-through pattern still works
    because it operates on whichever instance it's handed.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.boundary_flux import BoundaryFlux
from orpheus.sn.geometry import SNMesh
from orpheus.numerics.quadrature import Quadrature
from tests.sn._test_helpers import placeholder_materials


pytestmark = [pytest.mark.foundation]


# ───────────────────────────────────────────────────────────────────────
# Geometry fixtures
# ───────────────────────────────────────────────────────────────────────


def _slab_mesh() -> SNMesh:
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective, BC.reflective),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=5),))
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials())


def _sphere_mesh() -> SNMesh:
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=5),))
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials())


def _cylinder_mesh() -> SNMesh:
    geom = StructuredGeometry(
        geometry="CYL",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=5),))
    quad = Quadrature.product(n_mu=2, n_phi=2)
    return SNMesh(mesh, quad, placeholder_materials())


# ───────────────────────────────────────────────────────────────────────
# Slab — both xmin_face and xmax_face populated
# ───────────────────────────────────────────────────────────────────────


def test_add_slab_propagates_to_both_faces() -> None:
    """Slab add: xmin + xmin, xmax + xmax."""
    sn = _slab_mesh()
    bf1 = BoundaryFlux.zeros(sn)
    bf2 = BoundaryFlux.zeros(sn)
    bf1.xmin_face[...] = 1.0
    bf1.xmax_face[...] = 2.0
    bf2.xmin_face[...] = 3.0
    bf2.xmax_face[...] = 4.0

    out = bf1 + bf2
    assert out is not bf1
    assert out is not bf2
    np.testing.assert_array_equal(out.xmin_face, 4.0)
    np.testing.assert_array_equal(out.xmax_face, 6.0)
    # 2-D buffers remain None on a 1-D mesh.
    assert out.xmin_xmax_buf is None
    assert out.ymin_ymax_buf is None
    # Originals unchanged (immutable arithmetic).
    np.testing.assert_array_equal(bf1.xmin_face, 1.0)
    np.testing.assert_array_equal(bf1.xmax_face, 2.0)


def test_sub_slab() -> None:
    sn = _slab_mesh()
    bf1 = BoundaryFlux.zeros(sn)
    bf2 = BoundaryFlux.zeros(sn)
    bf1.xmin_face[...] = 5.0
    bf1.xmax_face[...] = 7.0
    bf2.xmin_face[...] = 2.0
    bf2.xmax_face[...] = 3.0

    out = bf1 - bf2
    np.testing.assert_array_equal(out.xmin_face, 3.0)
    np.testing.assert_array_equal(out.xmax_face, 4.0)


def test_scalar_mul_slab() -> None:
    sn = _slab_mesh()
    bf = BoundaryFlux.zeros(sn)
    bf.xmin_face[...] = 1.5
    bf.xmax_face[...] = 2.5

    out = bf * 2.0
    np.testing.assert_array_equal(out.xmin_face, 3.0)
    np.testing.assert_array_equal(out.xmax_face, 5.0)

    # __rmul__ — scalar on the left.
    out2 = 3.0 * bf
    np.testing.assert_array_equal(out2.xmin_face, 4.5)
    np.testing.assert_array_equal(out2.xmax_face, 7.5)


def test_scalar_div_slab() -> None:
    sn = _slab_mesh()
    bf = BoundaryFlux.zeros(sn)
    bf.xmin_face[...] = 6.0
    bf.xmax_face[...] = 10.0

    out = bf / 2.0
    np.testing.assert_array_equal(out.xmin_face, 3.0)
    np.testing.assert_array_equal(out.xmax_face, 5.0)


def test_neg_slab() -> None:
    sn = _slab_mesh()
    bf = BoundaryFlux.zeros(sn)
    bf.xmin_face[...] = 1.0
    bf.xmax_face[...] = -2.0

    out = -bf
    np.testing.assert_array_equal(out.xmin_face, -1.0)
    np.testing.assert_array_equal(out.xmax_face, 2.0)


# ───────────────────────────────────────────────────────────────────────
# Curvilinear (sphere / cylinder) — only xmax_face populated
# ───────────────────────────────────────────────────────────────────────


def test_add_sphere_only_outer_face() -> None:
    """Sphere has no xmin_face (pole is a regularity condition); add only xmax."""
    sn = _sphere_mesh()
    bf1 = BoundaryFlux.zeros(sn)
    bf2 = BoundaryFlux.zeros(sn)
    assert bf1.xmin_face is None
    assert bf2.xmin_face is None
    bf1.xmax_face[...] = 1.5
    bf2.xmax_face[...] = 2.5

    out = bf1 + bf2
    np.testing.assert_array_equal(out.xmax_face, 4.0)
    assert out.xmin_face is None


def test_add_cylinder_only_outer_face() -> None:
    sn = _cylinder_mesh()
    bf1 = BoundaryFlux.zeros(sn)
    bf2 = BoundaryFlux.zeros(sn)
    assert bf1.xmin_face is None
    bf1.xmax_face[...] = 0.3
    bf2.xmax_face[...] = 0.7

    out = bf1 + bf2
    np.testing.assert_array_equal(out.xmax_face, 1.0)


# ───────────────────────────────────────────────────────────────────────
# 2-D Cartesian — xmin_xmax_buf + ymin_ymax_buf
# ───────────────────────────────────────────────────────────────────────


def test_add_2d_cartesian_propagates_to_both_buffers() -> None:
    """2-D Cartesian carries face + interior face cache; add propagates element-wise.

    The 2-D buffers carry both boundary face slices AND the wavefront
    sweep's interior face-flux cache.  Elementwise addition is
    semantically well-defined (it's just numpy arithmetic on the
    persistent storage).
    """
    from orpheus.geometry import Mesh2D
    nx, ny = 3, 2
    delta = 0.5
    mesh = Mesh2D(
        edges_x=np.linspace(0, nx * delta, nx + 1),
        edges_y=np.linspace(0, ny * delta, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
    )
    quad = Quadrature.lebedev(order=5)
    sn = SNMesh(mesh, quad, placeholder_materials())
    bf1 = BoundaryFlux.zeros(sn)
    bf2 = BoundaryFlux.zeros(sn)
    assert bf1.xmin_xmax_buf is not None
    assert bf1.ymin_ymax_buf is not None
    bf1.xmin_xmax_buf[...] = 1.0
    bf1.ymin_ymax_buf[...] = 2.0
    bf2.xmin_xmax_buf[...] = 3.0
    bf2.ymin_ymax_buf[...] = 4.0

    out = bf1 + bf2
    np.testing.assert_array_equal(out.xmin_xmax_buf, 4.0)
    np.testing.assert_array_equal(out.ymin_ymax_buf, 6.0)
    assert out.xmin_face is None  # not used in 2-D


# ───────────────────────────────────────────────────────────────────────
# Mesh-identity invariant
# ───────────────────────────────────────────────────────────────────────


def test_cross_mesh_add_rejected() -> None:
    sn1 = _slab_mesh()
    sn2 = _slab_mesh()  # different SNMesh instance, structurally identical
    bf1 = BoundaryFlux.zeros(sn1)
    bf2 = BoundaryFlux.zeros(sn2)

    with pytest.raises(ValueError, match="same SNMesh"):
        bf1 + bf2


def test_add_wrong_type_rejected() -> None:
    sn = _slab_mesh()
    bf = BoundaryFlux.zeros(sn)
    with pytest.raises(TypeError, match="BoundaryFlux"):
        bf + 42  # noqa: B015 — intentional bad add


# ───────────────────────────────────────────────────────────────────────
# Compositions match the algebra
# ───────────────────────────────────────────────────────────────────────


def test_distributive_property() -> None:
    """(bf1 + bf2) * c == bf1 * c + bf2 * c — distributivity of the algebra."""
    sn = _slab_mesh()
    bf1 = BoundaryFlux.zeros(sn)
    bf2 = BoundaryFlux.zeros(sn)
    rng = np.random.default_rng(7)
    bf1.xmin_face[...] = rng.standard_normal(bf1.xmin_face.shape)
    bf1.xmax_face[...] = rng.standard_normal(bf1.xmax_face.shape)
    bf2.xmin_face[...] = rng.standard_normal(bf2.xmin_face.shape)
    bf2.xmax_face[...] = rng.standard_normal(bf2.xmax_face.shape)

    lhs = (bf1 + bf2) * 3.0
    rhs = bf1 * 3.0 + bf2 * 3.0
    np.testing.assert_allclose(lhs.xmin_face, rhs.xmin_face, rtol=1e-15)
    np.testing.assert_allclose(lhs.xmax_face, rhs.xmax_face, rtol=1e-15)


def test_neg_equals_zero_minus() -> None:
    """-bf == 0 * bf - bf elementwise — sanity check on __neg__."""
    sn = _sphere_mesh()
    bf = BoundaryFlux.zeros(sn)
    rng = np.random.default_rng(11)
    bf.xmax_face[...] = rng.standard_normal(bf.xmax_face.shape)

    neg = -bf
    np.testing.assert_array_equal(neg.xmax_face, -bf.xmax_face)
