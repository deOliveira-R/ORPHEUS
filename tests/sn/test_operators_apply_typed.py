r"""L0 foundation: operator ``apply(AngularFlux) → AngularFlux`` (R-1 Step 3).

The four SN operators in the algebra :math:`(L + C - S - F)` each
expose a typed :class:`~orpheus.sn.angular_flux.AngularFlux` apply
surface so the algebra reads as the math:

    (L + C - S).apply(psi)  →  AngularFlux

Per-operator boundary semantics pinned here:

* :class:`~orpheus.sn.operator.StreamingOperator` (``L``) — the ONLY
  operator that emits a non-zero face residual into
  ``result.boundary``.  The B1'' face equation
  :math:`(WDD\text{-propagated face}) − \psi_{\rm face}` lives in
  ``result.boundary.xmax_face`` (and, for slab, ``xmin_face``).
* :class:`~orpheus.sn.operator.CollisionOperator` (``C``) — volumetric;
  ``result.boundary`` is the auto-allocated zero.
* :class:`~orpheus.sn.scattering.ScatteringOperator` (``S``) —
  volumetric secondary-emission; ``result.boundary`` is zero.
* :class:`~orpheus.sn.fission.FissionOperator` (``F``) — volumetric
  rank-1 emission; ``result.boundary`` is zero.

Equivalence tests pin that the typed and packed (legacy bare-ndarray)
paths compute the SAME action — typed routes through ``boundary``
buffers while packed routes through the B1''-aware flat layout, but
both consume the same compute kernel
(:func:`~orpheus.sn.operator.transport_operator_matvec_unified` for L,
:meth:`CollisionOperator._sigma_at_unknowns` for C, etc.).
Round-trip identity via :meth:`AngularFlux.from_flat_with_traces` /
:meth:`AngularFlux.to_flat_with_traces` is the bridge.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.angular_flux import AngularFlux
from orpheus.sn.fission import FissionOperator
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import CollisionOperator, StreamingOperator
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.scalar_flux import ScalarFlux
from orpheus.sn.scattering import ScatteringOperator
from tests.sn._test_helpers import placeholder_materials


pytestmark = [pytest.mark.foundation]


# ───────────────────────────────────────────────────────────────────────
# Geometry fixtures (slab + sphere — the two B1'' face geometries)
# ───────────────────────────────────────────────────────────────────────


def _slab_mesh(nx: int = 5, n_ord: int = 4, ng: int = 2) -> SNMesh:
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective, BC.reflective),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = GaussLegendre1D.create(n_ordinates=n_ord)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _sphere_mesh(nx: int = 5, n_ord: int = 4, ng: int = 2) -> SNMesh:
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = GaussLegendre1D.create(n_ordinates=n_ord)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


GEOMETRIES = [("slab", _slab_mesh), ("sphere", _sphere_mesh)]


def _random_psi(sn: SNMesh, seed: int) -> AngularFlux:
    """Build an :class:`AngularFlux` with non-trivial values + boundary."""
    rng = np.random.default_rng(seed)
    N, ng, nx, ny = sn.quad.N, sn.ng, sn.nx, sn.ny
    psi = AngularFlux(rng.standard_normal((N, ng, nx, ny)), sn)
    # Populate the boundary face buffers with non-zero random values so
    # the typed branch exercises non-trivial face state.
    if psi.boundary.xmax_face is not None:
        psi.boundary.xmax_face[...] = rng.standard_normal(
            psi.boundary.xmax_face.shape,
        )
    if psi.boundary.xmin_face is not None:
        psi.boundary.xmin_face[...] = rng.standard_normal(
            psi.boundary.xmin_face.shape,
        )
    return psi


# ───────────────────────────────────────────────────────────────────────
# F — FissionOperator.apply(AngularFlux) lift
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("name,builder", GEOMETRIES)
def test_F_apply_angular_flux_returns_angular_flux(name, builder) -> None:
    """F.apply(AngularFlux) returns AngularFlux with zero boundary."""
    sn = builder()
    psi = _random_psi(sn, seed=1)
    F = FissionOperator.from_solver_data(mat_xs=sn.material_xs_field())

    Fpsi = F.apply(psi)
    assert isinstance(Fpsi, AngularFlux)
    assert Fpsi.values.shape == psi.values.shape
    assert Fpsi.mesh is sn
    # Boundary is the auto-allocated zero — F has no face contribution.
    if Fpsi.boundary.xmax_face is not None:
        np.testing.assert_array_equal(Fpsi.boundary.xmax_face, 0.0)
    if Fpsi.boundary.xmin_face is not None:
        np.testing.assert_array_equal(Fpsi.boundary.xmin_face, 0.0)


@pytest.mark.parametrize("name,builder", GEOMETRIES)
def test_F_typed_lift_equivalent_to_scalar(name, builder) -> None:
    """F.apply(AngularFlux) broadcasts F.apply(integrate_angular(AngularFlux))."""
    sn = builder()
    psi = _random_psi(sn, seed=2)
    F = FissionOperator.from_solver_data(mat_xs=sn.material_xs_field())

    # Typed path: F(psi) — AngularFlux
    Fpsi_typed = F.apply(psi)
    # Scalar path: F(integrate_angular(psi)) — ScalarFlux, then broadcast
    phi = psi.integrate_angular()
    F_phi = F.apply(phi)
    N = sn.quad.N
    expected_values = np.broadcast_to(
        F_phi.values[None, :, :, :], Fpsi_typed.values.shape,
    )
    np.testing.assert_allclose(Fpsi_typed.values, expected_values, rtol=1e-14)


# ───────────────────────────────────────────────────────────────────────
# S — ScatteringOperator.apply(AngularFlux) (already existed pre-R-1 Step 3)
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("name,builder", GEOMETRIES)
def test_S_apply_angular_flux_returns_angular_flux_zero_boundary(
    name, builder,
) -> None:
    """S.apply(AngularFlux) returns AngularFlux; boundary stays zero."""
    sn = builder()
    psi = _random_psi(sn, seed=3)
    S = ScatteringOperator.from_solver_data(
        mat_xs=sn.material_xs_field(),
        quadrature=sn.quad,
        scattering_order=0,
    )

    Spsi = S.apply(psi)
    assert isinstance(Spsi, AngularFlux)
    assert Spsi.values.shape == psi.values.shape
    # S is volumetric; result's boundary is auto-allocated zero.
    if Spsi.boundary.xmax_face is not None:
        np.testing.assert_array_equal(Spsi.boundary.xmax_face, 0.0)


@pytest.mark.parametrize("name,builder", GEOMETRIES)
def test_S_typed_matches_raw(name, builder) -> None:
    """S.apply(AngularFlux).values == S.apply(np.ndarray) bit-exact."""
    sn = builder()
    psi = _random_psi(sn, seed=4)
    S = ScatteringOperator.from_solver_data(
        mat_xs=sn.material_xs_field(),
        quadrature=sn.quad,
        scattering_order=0,
    )
    typed_values = S.apply(psi).values
    raw_values = S.apply(psi.values)
    np.testing.assert_array_equal(typed_values, raw_values)


# ───────────────────────────────────────────────────────────────────────
# C — CollisionOperator.apply(AngularFlux)
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("name,builder", GEOMETRIES)
def test_C_apply_angular_flux_returns_angular_flux_zero_boundary(
    name, builder,
) -> None:
    """C.apply(AngularFlux) returns AngularFlux; boundary stays zero."""
    sn = builder()
    psi = _random_psi(sn, seed=5)
    sigma = np.ones((sn.ng, sn.nx, sn.ny)) * 0.7
    C = CollisionOperator(sn, sigma)

    Cpsi = C.apply(psi)
    assert isinstance(Cpsi, AngularFlux)
    assert Cpsi.values.shape == psi.values.shape
    if Cpsi.boundary.xmax_face is not None:
        np.testing.assert_array_equal(Cpsi.boundary.xmax_face, 0.0)
    if Cpsi.boundary.xmin_face is not None:
        np.testing.assert_array_equal(Cpsi.boundary.xmin_face, 0.0)


@pytest.mark.parametrize("name,builder", GEOMETRIES)
def test_C_diagonal_action(name, builder) -> None:
    """C.apply(AngularFlux).values == σ ⊙ psi.values element-wise."""
    sn = builder()
    psi = _random_psi(sn, seed=6)
    rng = np.random.default_rng(7)
    sigma = 0.3 + 0.5 * rng.random((sn.ng, sn.nx, sn.ny))
    C = CollisionOperator(sn, sigma)

    Cpsi = C.apply(psi)
    expected = sigma[None, :, :, :] * psi.values
    np.testing.assert_array_equal(Cpsi.values, expected)


@pytest.mark.parametrize("name,builder", GEOMETRIES)
def test_C_typed_matches_packed(name, builder) -> None:
    """C.apply(AngularFlux) and C.apply(to_flat_with_traces()) agree on cells."""
    sn = builder()
    psi = _random_psi(sn, seed=8)
    rng = np.random.default_rng(9)
    sigma = 0.3 + 0.5 * rng.random((sn.ng, sn.nx, sn.ny))
    C = CollisionOperator(sn, sigma)

    Cpsi_typed = C.apply(psi)
    flat_psi = psi.to_flat_with_traces()
    Cflat = C.apply(flat_psi)
    Cpsi_recon = AngularFlux.from_flat_with_traces(Cflat, sn)

    # Cell-centre block: typed and packed agree bit-exact.
    np.testing.assert_array_equal(Cpsi_typed.values, Cpsi_recon.values)
    # Face block: typed is zero; packed routes a zero through the
    # face-block (the cell-only σ-multiply leaves face slots at 0).
    if Cpsi_recon.boundary.xmax_face is not None:
        np.testing.assert_array_equal(Cpsi_recon.boundary.xmax_face, 0.0)


# ───────────────────────────────────────────────────────────────────────
# L — StreamingOperator.apply(AngularFlux) — load-bearing
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("name,builder", GEOMETRIES)
def test_L_apply_angular_flux_returns_angular_flux(name, builder) -> None:
    """L.apply(AngularFlux) returns AngularFlux; .boundary may be non-zero."""
    sn = builder()
    psi = _random_psi(sn, seed=10)
    rng = np.random.default_rng(11)
    sigma_t = 0.4 + 0.4 * rng.random((sn.ng, sn.nx, sn.ny))
    L = StreamingOperator(sn, sigma_t)

    Lpsi = L.apply(psi)
    assert isinstance(Lpsi, AngularFlux)
    assert Lpsi.values.shape == psi.values.shape
    # L emits the face residual into .boundary; for these random
    # inputs at least one face slot is non-zero.
    assert Lpsi.boundary.xmax_face is not None
    assert not np.allclose(Lpsi.boundary.xmax_face, 0.0)


@pytest.mark.parametrize("name,builder", GEOMETRIES)
def test_L_typed_matches_packed_round_trip(name, builder) -> None:
    """L typed and packed paths agree via from_flat/to_flat round-trip.

    The packed path runs ``solution_to_angular_flux_with_traces ∘
    transport_operator_matvec_unified ∘ pack_with_traces`` and
    subtracts σ_t ⊙ ψ at cell slots only.  The typed path runs the
    matvec directly on ``psi.values`` and stores the face residual in
    ``result.boundary``.  Both must produce identical numerics — same
    compute kernel, just different I/O wrapping.
    """
    sn = builder()
    psi = _random_psi(sn, seed=12)
    rng = np.random.default_rng(13)
    sigma_t = 0.4 + 0.4 * rng.random((sn.ng, sn.nx, sn.ny))
    L = StreamingOperator(sn, sigma_t)

    # Typed path
    Lpsi_typed = L.apply(psi)

    # Packed path
    flat_psi = psi.to_flat_with_traces()
    Lflat = L.apply(flat_psi)
    Lpsi_recon = AngularFlux.from_flat_with_traces(Lflat, sn)

    # Cell-centre values agree bit-exact.
    np.testing.assert_array_equal(Lpsi_typed.values, Lpsi_recon.values)
    # Face residual agrees on the populated ordinate slots.  The typed
    # output stores residual only at ``face_outer_ordinate`` /
    # ``face_inner_ordinate`` slots; the unpopulated ordinates stay 0
    # in both paths.
    np.testing.assert_array_equal(
        Lpsi_typed.boundary.xmax_face, Lpsi_recon.boundary.xmax_face,
    )
    if sn.curvature is None:  # slab — inner face exists
        np.testing.assert_array_equal(
            Lpsi_typed.boundary.xmin_face, Lpsi_recon.boundary.xmin_face,
        )


@pytest.mark.parametrize("name,builder", GEOMETRIES)
def test_L_cross_mesh_rejected(name, builder) -> None:
    """L.apply(AngularFlux) rejects a flux bound to a different SNMesh."""
    sn1 = builder()
    sn2 = builder()  # distinct instance, structurally identical
    sigma_t = np.ones((sn1.ng, sn1.nx, sn1.ny))
    L = StreamingOperator(sn1, sigma_t)
    psi_foreign = AngularFlux(
        np.zeros((sn2.quad.N, sn2.ng, sn2.nx, sn2.ny)), sn2,
    )
    with pytest.raises(ValueError, match="SAME SNMesh"):
        L.apply(psi_foreign)


# ───────────────────────────────────────────────────────────────────────
# Compose: (L + C - S).apply(psi) reads as the algebra
# ───────────────────────────────────────────────────────────────────────


def test_operator_algebra_reads_as_math() -> None:
    """``(L + C - S).apply(psi)`` returns a typed AngularFlux end-to-end.

    Pins the load-bearing R-1 design contract: the operator algebra
    consumes :class:`AngularFlux` with one argument (no second
    boundary arg).  Validates only the type contract and the algebraic
    propagation of ``.boundary`` — numerics are checked elsewhere.
    """
    sn = _slab_mesh()
    psi = _random_psi(sn, seed=14)
    sigma_t = np.full((sn.ng, sn.nx, sn.ny), 0.7)
    L = StreamingOperator(sn, sigma_t)
    C = CollisionOperator(sn, sigma_t * 0.5)
    S = ScatteringOperator.from_solver_data(
        mat_xs=sn.material_xs_field(),
        quadrature=sn.quad,
        scattering_order=0,
    )

    # The algebra reads as the math.
    Lpsi = L.apply(psi)
    Cpsi = C.apply(psi)
    Spsi = S.apply(psi)
    rhs = Lpsi + Cpsi - Spsi

    assert isinstance(rhs, AngularFlux)
    # Boundary propagation through arithmetic — L contributes the only
    # non-zero face residual; C and S boundaries are zero, so
    # rhs.boundary == Lpsi.boundary at face slots.
    np.testing.assert_array_equal(
        rhs.boundary.xmax_face, Lpsi.boundary.xmax_face,
    )
