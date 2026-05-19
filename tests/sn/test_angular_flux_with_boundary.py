"""L0 foundation: :class:`AngularFlux` embeds :class:`BoundaryFlux` (R-1 Step 2).

Issue #197 PR-TYPED-2 step R-1.2 — :class:`AngularFlux` now carries its
own ``boundary`` field, auto-allocated to zeros via
:meth:`BoundaryFlux.zeros` when not supplied.  Arithmetic on
:class:`AngularFlux` propagates to ``.boundary`` element-wise via
:class:`BoundaryFlux`'s dunders (R-1 Step 1).

Layout adapters tested here:

* :meth:`AngularFlux.from_flat_with_traces` decodes a B1''-aware packed
  1-D vector into ``(values, boundary.xmax_face, boundary.xmin_face)``.
* :meth:`AngularFlux.to_flat_with_traces` is the inverse — round-trip
  identity at bit-exact precision.

Semantics pinned:
  * Auto-allocation on construction (``boundary=None`` → zeros).
  * Arithmetic propagation (``(ψ₁+ψ₂).boundary == ψ₁.boundary+ψ₂.boundary``).
  * Mesh-identity invariant inherits from :class:`BoundaryFlux`
    (``__add__`` on mismatched meshes raises ``ValueError``).
  * ``from_flat_with_traces / to_flat_with_traces`` is the bit-exact
    inverse pair at the scipy boundary — the ONE ravel/reshape site.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.angular_flux import AngularFlux
from orpheus.sn.boundary_flux import BoundaryFlux
from orpheus.sn.geometry import SNMesh
from orpheus.sn.quadrature import GaussLegendre1D
from tests.sn._test_helpers import placeholder_materials


pytestmark = [pytest.mark.foundation]


# ───────────────────────────────────────────────────────────────────────
# Geometry fixtures
# ───────────────────────────────────────────────────────────────────────


def _slab_mesh(nx: int = 5, n_ordinates: int = 4, ng: int = 1) -> SNMesh:
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective, BC.reflective),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = GaussLegendre1D.create(n_ordinates=n_ordinates)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _sphere_mesh(nx: int = 5, n_ordinates: int = 4, ng: int = 1) -> SNMesh:
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = GaussLegendre1D.create(n_ordinates=n_ordinates)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


# ───────────────────────────────────────────────────────────────────────
# Auto-allocation — every AngularFlux has a boundary
# ───────────────────────────────────────────────────────────────────────


def test_auto_allocates_boundary_zeros() -> None:
    """Constructing without ``boundary`` auto-allocates a zeros BoundaryFlux."""
    sn = _slab_mesh()
    N, ng, nx, ny = sn.quad.N, sn.ng, sn.nx, sn.ny
    psi = AngularFlux(np.ones((N, ng, nx, ny)), sn)

    assert psi.boundary is not None
    assert isinstance(psi.boundary, BoundaryFlux)
    assert psi.boundary.mesh is sn
    # Slab → both face buffers allocated as zeros.
    assert psi.boundary.xmin_face is not None
    assert psi.boundary.xmax_face is not None
    np.testing.assert_array_equal(psi.boundary.xmin_face, 0.0)
    np.testing.assert_array_equal(psi.boundary.xmax_face, 0.0)


def test_auto_allocates_sphere_only_outer() -> None:
    """Sphere geometry → only xmax_face allocated; xmin_face stays None."""
    sn = _sphere_mesh()
    N, ng, nx, ny = sn.quad.N, sn.ng, sn.nx, sn.ny
    psi = AngularFlux(np.zeros((N, ng, nx, ny)), sn)

    assert psi.boundary.xmax_face is not None
    assert psi.boundary.xmin_face is None


def test_supplied_boundary_is_kept_not_replaced() -> None:
    """A user-supplied ``boundary`` is stored verbatim — no auto-replace."""
    sn = _slab_mesh()
    N, ng, nx, ny = sn.quad.N, sn.ng, sn.nx, sn.ny
    bf = BoundaryFlux.zeros(sn)
    bf.xmax_face[...] = 7.5
    psi = AngularFlux(np.zeros((N, ng, nx, ny)), sn, boundary=bf)

    assert psi.boundary is bf
    np.testing.assert_array_equal(psi.boundary.xmax_face, 7.5)


# ───────────────────────────────────────────────────────────────────────
# Arithmetic propagates to boundary
# ───────────────────────────────────────────────────────────────────────


def test_add_propagates_to_boundary() -> None:
    """(ψ₁ + ψ₂).boundary == ψ₁.boundary + ψ₂.boundary."""
    sn = _slab_mesh()
    N, ng, nx, ny = sn.quad.N, sn.ng, sn.nx, sn.ny

    psi1 = AngularFlux(np.full((N, ng, nx, ny), 1.0), sn)
    psi2 = AngularFlux(np.full((N, ng, nx, ny), 2.0), sn)
    psi1.boundary.xmin_face[...] = 0.5
    psi1.boundary.xmax_face[...] = 1.5
    psi2.boundary.xmin_face[...] = 2.5
    psi2.boundary.xmax_face[...] = 3.5

    out = psi1 + psi2
    np.testing.assert_array_equal(out.values, 3.0)
    np.testing.assert_array_equal(out.boundary.xmin_face, 3.0)
    np.testing.assert_array_equal(out.boundary.xmax_face, 5.0)
    # Operands unchanged (immutable arithmetic).
    np.testing.assert_array_equal(psi1.boundary.xmin_face, 0.5)


def test_sub_propagates_to_boundary() -> None:
    sn = _slab_mesh()
    N, ng, nx, ny = sn.quad.N, sn.ng, sn.nx, sn.ny

    psi1 = AngularFlux(np.full((N, ng, nx, ny), 7.0), sn)
    psi2 = AngularFlux(np.full((N, ng, nx, ny), 2.0), sn)
    psi1.boundary.xmin_face[...] = 5.0
    psi1.boundary.xmax_face[...] = 9.0
    psi2.boundary.xmin_face[...] = 1.0
    psi2.boundary.xmax_face[...] = 3.0

    out = psi1 - psi2
    np.testing.assert_array_equal(out.values, 5.0)
    np.testing.assert_array_equal(out.boundary.xmin_face, 4.0)
    np.testing.assert_array_equal(out.boundary.xmax_face, 6.0)


def test_scalar_mul_propagates_to_boundary() -> None:
    sn = _slab_mesh()
    N, ng, nx, ny = sn.quad.N, sn.ng, sn.nx, sn.ny

    psi = AngularFlux(np.full((N, ng, nx, ny), 3.0), sn)
    psi.boundary.xmin_face[...] = 2.0
    psi.boundary.xmax_face[...] = 4.0

    out = psi * 2.5
    np.testing.assert_array_equal(out.values, 7.5)
    np.testing.assert_array_equal(out.boundary.xmin_face, 5.0)
    np.testing.assert_array_equal(out.boundary.xmax_face, 10.0)

    # __rmul__ — scalar on the left.
    out2 = 0.5 * psi
    np.testing.assert_array_equal(out2.values, 1.5)
    np.testing.assert_array_equal(out2.boundary.xmin_face, 1.0)


def test_scalar_div_propagates_to_boundary() -> None:
    sn = _slab_mesh()
    N, ng, nx, ny = sn.quad.N, sn.ng, sn.nx, sn.ny

    psi = AngularFlux(np.full((N, ng, nx, ny), 6.0), sn)
    psi.boundary.xmin_face[...] = 4.0
    psi.boundary.xmax_face[...] = 8.0

    out = psi / 2.0
    np.testing.assert_array_equal(out.values, 3.0)
    np.testing.assert_array_equal(out.boundary.xmin_face, 2.0)
    np.testing.assert_array_equal(out.boundary.xmax_face, 4.0)


def test_neg_propagates_to_boundary() -> None:
    sn = _slab_mesh()
    N, ng, nx, ny = sn.quad.N, sn.ng, sn.nx, sn.ny

    psi = AngularFlux(np.full((N, ng, nx, ny), 1.5), sn)
    psi.boundary.xmin_face[...] = -2.0
    psi.boundary.xmax_face[...] = 3.0

    out = -psi
    np.testing.assert_array_equal(out.values, -1.5)
    np.testing.assert_array_equal(out.boundary.xmin_face, 2.0)
    np.testing.assert_array_equal(out.boundary.xmax_face, -3.0)


def test_distributive_property() -> None:
    """(ψ₁ + ψ₂) * c == ψ₁ * c + ψ₂ * c — algebra distributivity."""
    sn = _slab_mesh()
    N, ng, nx, ny = sn.quad.N, sn.ng, sn.nx, sn.ny
    rng = np.random.default_rng(13)
    psi1 = AngularFlux(rng.standard_normal((N, ng, nx, ny)), sn)
    psi2 = AngularFlux(rng.standard_normal((N, ng, nx, ny)), sn)
    psi1.boundary.xmin_face[...] = rng.standard_normal(psi1.boundary.xmin_face.shape)
    psi1.boundary.xmax_face[...] = rng.standard_normal(psi1.boundary.xmax_face.shape)
    psi2.boundary.xmin_face[...] = rng.standard_normal(psi2.boundary.xmin_face.shape)
    psi2.boundary.xmax_face[...] = rng.standard_normal(psi2.boundary.xmax_face.shape)

    lhs = (psi1 + psi2) * 3.0
    rhs = psi1 * 3.0 + psi2 * 3.0
    # FP-non-associativity bound: (lhs vs rhs) differ in reduction order
    # by a single add/mul; principled bound is a few × ULP per
    # vv-principles "Bit-identity vs principled-equivalence".
    np.testing.assert_allclose(lhs.values, rhs.values, rtol=1e-14)
    np.testing.assert_allclose(
        lhs.boundary.xmin_face, rhs.boundary.xmin_face, rtol=1e-14,
    )
    np.testing.assert_allclose(
        lhs.boundary.xmax_face, rhs.boundary.xmax_face, rtol=1e-14,
    )


# ───────────────────────────────────────────────────────────────────────
# Mesh-identity invariant on .boundary
# ───────────────────────────────────────────────────────────────────────


def test_cross_mesh_add_rejected_via_boundary() -> None:
    """Cross-mesh add on AngularFlux is rejected even if .values match shape."""
    sn1 = _slab_mesh()
    sn2 = _slab_mesh()  # different SNMesh instance
    N, ng, nx, ny = sn1.quad.N, sn1.ng, sn1.nx, sn1.ny
    psi1 = AngularFlux(np.zeros((N, ng, nx, ny)), sn1)
    psi2 = AngularFlux(np.zeros((N, ng, nx, ny)), sn2)

    with pytest.raises(ValueError, match="distinct SNMesh"):
        psi1 + psi2


# ───────────────────────────────────────────────────────────────────────
# from_flat_with_traces / to_flat_with_traces — round-trip identity
# ───────────────────────────────────────────────────────────────────────


def test_round_trip_slab_bit_exact() -> None:
    """from_flat_with_traces ∘ to_flat_with_traces == identity on slab."""
    sn = _slab_mesh(nx=5, n_ordinates=4, ng=2)
    # Build a flat vector of the right size by reading n_unknowns.
    from orpheus.sn.operator import build_equation_map_with_traces
    eq_map = build_equation_map_with_traces(
        sn.nx, sn.quad, sn.ng, has_inner_bc=True,
    )
    rng = np.random.default_rng(31)
    flat_in = rng.standard_normal(eq_map.n_unknowns)

    psi = AngularFlux.from_flat_with_traces(flat_in, sn)
    flat_out = psi.to_flat_with_traces()

    np.testing.assert_array_equal(flat_in, flat_out)


def test_round_trip_sphere_bit_exact() -> None:
    """Sphere geometry → no inner-face block; round-trip stays bit-exact."""
    sn = _sphere_mesh(nx=5, n_ordinates=4, ng=2)
    from orpheus.sn.operator import build_equation_map_with_traces
    eq_map = build_equation_map_with_traces(
        sn.nx, sn.quad, sn.ng, has_inner_bc=False,
    )
    rng = np.random.default_rng(41)
    flat_in = rng.standard_normal(eq_map.n_unknowns)

    psi = AngularFlux.from_flat_with_traces(flat_in, sn)
    flat_out = psi.to_flat_with_traces()

    np.testing.assert_array_equal(flat_in, flat_out)


def test_from_flat_populates_boundary_face() -> None:
    """from_flat_with_traces puts face values in boundary.xmax_face / xmin_face."""
    sn = _slab_mesh(nx=3, n_ordinates=4, ng=1)
    from orpheus.sn.operator import build_equation_map_with_traces
    eq_map = build_equation_map_with_traces(
        sn.nx, sn.quad, sn.ng, has_inner_bc=True,
    )
    rng = np.random.default_rng(7)
    flat = rng.standard_normal(eq_map.n_unknowns)

    psi = AngularFlux.from_flat_with_traces(flat, sn)

    # Boundary must be populated at the outer-ordinate slots.
    for n_out in eq_map.face_outer_ordinate:
        assert not np.array_equal(psi.boundary.xmax_face[n_out], np.zeros(sn.ng))
    # And at the inner-ordinate slots (slab only).
    for n_in in eq_map.face_inner_ordinate:
        assert not np.array_equal(psi.boundary.xmin_face[n_in], np.zeros(sn.ng))


def test_from_flat_size_mismatch_rejected() -> None:
    """Wrong flat.size raises ValueError naming n_unknowns."""
    sn = _slab_mesh()
    with pytest.raises(ValueError, match="n_unknowns"):
        AngularFlux.from_flat_with_traces(np.zeros(3), sn)


# ───────────────────────────────────────────────────────────────────────
# copy() carries owned boundary
# ───────────────────────────────────────────────────────────────────────


def test_copy_clones_boundary_buffers() -> None:
    """copy() deep-clones both ``values`` AND ``boundary`` — no sharing."""
    sn = _slab_mesh()
    N, ng, nx, ny = sn.quad.N, sn.ng, sn.nx, sn.ny
    psi = AngularFlux(np.full((N, ng, nx, ny), 2.0), sn)
    psi.boundary.xmin_face[...] = 1.0
    psi.boundary.xmax_face[...] = 3.0

    psi2 = psi.copy()
    assert psi2 is not psi
    assert psi2.boundary is not psi.boundary
    assert psi2.boundary.xmin_face is not psi.boundary.xmin_face

    # Mutating the copy doesn't touch the original.
    psi2.boundary.xmin_face[...] = -99.0
    np.testing.assert_array_equal(psi.boundary.xmin_face, 1.0)
