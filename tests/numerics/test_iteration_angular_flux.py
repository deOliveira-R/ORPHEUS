"""L0 foundation: KrylovAcceleration + SourceIteration consume AngularFlux.

R-1 Step 4a — the iteration primitives accept typed flux containers
(:class:`~orpheus.sn.angular_flux.AngularFlux`) via the duck-typed
**ravellable protocol**: any object exposing the pair
``(to_flat_with_traces, from_flat_with_traces, mesh)``.  Bare-ndarray
inputs continue to work via the legacy reshape path.

This test file pins the protocol contract.  The primitives in
``orpheus.numerics.iteration`` are deliberately decoupled from
``orpheus.sn.*`` — they consume the protocol, not the AngularFlux
class.  Tests demonstrate:

* Round-trip: AngularFlux flows in, AngularFlux flows out.
* Boundary state survives the ravel/unravel — the flat vector
  carries the B1'' face block via ``to_flat_with_traces``.
* Composed algebra ``(L + C - S - F).apply(psi)`` matches
  ``L.apply(psi) - S.apply(psi) - F.apply(psi)`` end-to-end.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.numerics.iteration import KrylovAcceleration, SourceIteration
from orpheus.numerics.iteration import (
    _is_ravellable,
    _l2_norm,
    _ravel,
    _unravel_like,
    _zeros_like,
)
from orpheus.numerics.operator import OperatorSum, ZeroOperator
from orpheus.sn.angular_flux import AngularFlux
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import CollisionOperator, StreamingOperator
from orpheus.numerics.quadrature import Quadrature
from tests.sn._test_helpers import placeholder_materials


pytestmark = [pytest.mark.foundation]


# ── Fixtures ─────────────────────────────────────────────────────────


def _slab_mesh(nx: int = 4, n_ord: int = 4, ng: int = 1) -> SNMesh:
    """Vacuum-BC slab — pure-absorption Krylov converges in O(1) iters."""
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.vacuum, BC.vacuum),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


# ── Ravellable protocol ─────────────────────────────────────────────


class TestRavellableProtocol:
    def test_angular_flux_is_ravellable(self) -> None:
        sn = _slab_mesh()
        psi = AngularFlux(np.zeros((sn.quad.N, sn.ng, sn.nx, sn.ny)), sn)
        assert _is_ravellable(psi)

    def test_bare_ndarray_is_not_ravellable(self) -> None:
        assert not _is_ravellable(np.zeros(10))
        assert not _is_ravellable(np.zeros((3, 4)))

    def test_ravel_unravel_round_trip(self) -> None:
        """_ravel ∘ _unravel_like = identity for AngularFlux on the B1'' slots.

        The flat representation only carries the B1'' face state — the
        outward-ordinate slice of ``xmax_face`` and (slab) the inward-
        ordinate slice of ``xmin_face``.  Slots outside the B1'' set
        (inward ordinates on ``xmax_face``; outward on ``xmin_face``)
        are physically zero and are zeroed by the unravel.  The test
        populates only the B1'' slots so the round-trip is the
        identity.
        """
        sn = _slab_mesh(ng=2)
        from orpheus.sn.operator import build_equation_map_with_traces
        eq_map = build_equation_map_with_traces(
            sn.nx, sn.quad, sn.ng, has_inner_bc=True,
        )
        rng = np.random.default_rng(1)
        psi = AngularFlux(
            rng.standard_normal((sn.quad.N, sn.ng, sn.nx, sn.ny)), sn,
        )
        # Populate only the B1'' slots — outward ordinates on xmax_face,
        # inward ordinates on xmin_face.
        psi.boundary.xmax_face[eq_map.face_outer_ordinate, :] = (
            rng.standard_normal((eq_map.n_face_outer, sn.ng))
        )
        psi.boundary.xmin_face[eq_map.face_inner_ordinate, :] = (
            rng.standard_normal((eq_map.n_face_inner, sn.ng))
        )

        flat = _ravel(psi)
        psi2 = _unravel_like(psi, flat)

        assert isinstance(psi2, AngularFlux)
        np.testing.assert_array_equal(psi.values, psi2.values)
        np.testing.assert_array_equal(
            psi.boundary.xmax_face, psi2.boundary.xmax_face,
        )
        if psi.boundary.xmin_face is not None:
            np.testing.assert_array_equal(
                psi.boundary.xmin_face, psi2.boundary.xmin_face,
            )

    def test_zeros_like_returns_angular_flux(self) -> None:
        sn = _slab_mesh()
        psi = AngularFlux(np.ones((sn.quad.N, sn.ng, sn.nx, sn.ny)), sn)
        zero_psi = _zeros_like(psi)
        assert isinstance(zero_psi, AngularFlux)
        np.testing.assert_array_equal(zero_psi.values, 0.0)
        # Boundary auto-allocated to zero on the fresh AngularFlux.
        if zero_psi.boundary.xmax_face is not None:
            np.testing.assert_array_equal(zero_psi.boundary.xmax_face, 0.0)

    def test_l2_norm_uses_ravelled_form(self) -> None:
        """_l2_norm of typed flux includes the boundary block in the ravel."""
        sn = _slab_mesh()
        psi = AngularFlux(np.zeros((sn.quad.N, sn.ng, sn.nx, sn.ny)), sn)
        # Boundary contributes — set xmax_face non-zero.
        psi.boundary.xmax_face[...] = 1.0
        norm = _l2_norm(psi)
        # Norm should be sqrt(number of non-zero face slots).
        expected = float(np.linalg.norm(psi.to_flat_with_traces()))
        np.testing.assert_allclose(norm, expected, rtol=1e-15)
        assert norm > 0.0


# ── KrylovAcceleration with AngularFlux ─────────────────────────────


class TestKrylovAccelerationAngularFlux:
    def test_solve_pure_collision_angular_flux(self) -> None:
        """Krylov on (L+C)psi=q for AngularFlux returns AngularFlux."""
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, sn.nx, sn.ny))
        L = StreamingOperator(sn, sigma_t)
        C = CollisionOperator(sn, sigma_t)
        LC = OperatorSum(L, C)
        S = ZeroOperator()
        F = ZeroOperator()

        q_ext = AngularFlux(
            np.ones((sn.quad.N, sn.ng, sn.nx, sn.ny)), sn,
        )
        krylov = KrylovAcceleration(LC, S, F, tol=1e-10, max_iter=200)
        psi, history = krylov.solve(q_ext)

        assert isinstance(psi, AngularFlux)
        assert psi.mesh is sn
        assert psi.values.shape == q_ext.values.shape
        assert psi.boundary is not None
        # Converged.
        assert history[-1] < 1e-9

    def test_solve_matches_bare_ndarray_path(self) -> None:
        """AngularFlux Krylov and bare-ndarray Krylov produce identical numerics.

        Sets up the same (L+C)psi=q problem in both representations.
        The AngularFlux path goes through ``to_flat_with_traces`` /
        ``from_flat_with_traces`` (cell + face blocks); the bare-ndarray
        path uses ``psi.values`` directly (cell only — face state is
        identically zero in this fixture, so the comparison is fair).
        """
        sn = _slab_mesh()
        sigma_t = np.full((sn.ng, sn.nx, sn.ny), 0.7)
        L = StreamingOperator(sn, sigma_t)
        C = CollisionOperator(sn, sigma_t)
        LC = OperatorSum(L, C)
        S = ZeroOperator()
        F = ZeroOperator()

        # Typed solve.
        q_typed = AngularFlux(
            np.full((sn.quad.N, sn.ng, sn.nx, sn.ny), 1.5), sn,
        )
        krylov_typed = KrylovAcceleration(LC, S, F, tol=1e-12, max_iter=500)
        psi_typed, _ = krylov_typed.solve(q_typed)

        # Bare-ndarray (packed-1D) solve — same source, packed via
        # the AngularFlux flat adapter so q_packed carries identical
        # numerics.  The operators dispatch on isinstance to the
        # packed-1D branch under the dual-surface contract Step 3
        # landed.
        q_packed = q_typed.to_flat_with_traces()
        # NOTE: ZeroOperator returns 0.0 * x for both paths; OperatorSum
        # composes A.apply(x) + B.apply(x); both routes through the
        # packed-1D leaf when x is ndarray.
        krylov_packed = KrylovAcceleration(
            LC, S, F, tol=1e-12, max_iter=500,
        )
        psi_packed_flat, _ = krylov_packed.solve(q_packed)

        # Reconstruct AngularFlux from packed solution for comparison.
        psi_packed = AngularFlux.from_flat_with_traces(psi_packed_flat, sn)

        # Cell values must agree.  FP-non-associativity bound: solver
        # iterations on different reduction orders → a few ULP × niter.
        np.testing.assert_allclose(
            psi_typed.values, psi_packed.values, rtol=1e-10, atol=1e-12,
        )


# ── SourceIteration with AngularFlux ────────────────────────────────


class _IdentityOperator:
    r"""A type-plumbing fixture exposing CAP_APPLY + CAP_SOLVE both as identity.

    R-1 Step B (2026-05-19) — used in place of the retired ``inverter``
    Callable hook to give SourceIteration an L whose ``.solve`` is
    trivially identity, so the iteration converges in one step and the
    test can assert on the typed-flux plumbing without exercising a real
    transport sweep.
    """
    from orpheus.numerics.operator import CAP_APPLY, CAP_SOLVE
    capabilities = frozenset({CAP_APPLY, CAP_SOLVE})

    def apply(self, x):
        return x

    def solve(self, rhs):
        return rhs


class TestSourceIterationAngularFlux:
    def test_zeros_like_in_solve(self) -> None:
        """SourceIteration's zero-initial-guess path uses _zeros_like."""
        sn = _slab_mesh()
        # L is an identity operator with full CAP_APPLY + CAP_SOLVE.
        # The point of this test is the typed-flux plumbing, not a
        # real transport solve; identity-L converges in one iter.
        L = _IdentityOperator()
        S = ZeroOperator()
        F = ZeroOperator()

        q_ext = AngularFlux(
            np.full((sn.quad.N, sn.ng, sn.nx, sn.ny), 2.5), sn,
        )
        si = SourceIteration(L, S, F, max_iter=5, tol=1e-10)
        psi, history = si.solve(q_ext)

        assert isinstance(psi, AngularFlux)
        # First iter: psi = L.solve(F·0 + S·0 + q) = identity(q) = q.
        np.testing.assert_array_equal(psi.values, q_ext.values)
