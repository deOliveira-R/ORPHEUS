"""Foundation tests for :class:`orpheus.sn.scattering.ScatteringOperator`.

Round 1.2 of Wave D of the SN reshape campaign (Issue #162). The
operator carries the same math :class:`SNSolver` used to expose under
``_add_scattering_source``, ``_build_aniso_scattering``, and
``_add_n2n_source`` (now thin delegators); these tests pin the lifted
math at the operator level so any drift would be observable here as
well as via the underscore-prefixed delegators.

The load-bearing test is the **bit-identical extraction** suite: a
synthetic ``(psi, phi, Q)`` triple is fed through both the new
:meth:`ScatteringOperator.apply` (and the in-place helpers
:meth:`add_iso_source` / :meth:`add_n2n_source` / :meth:`build_aniso_source`)
and the explicit per-cell reference implementations from
``test_solver_components.py``. The two paths must agree to round-off,
because the operator is a structural extraction, not a re-derivation.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture, make_mixture
from orpheus.geometry import Mesh2D
from orpheus.numerics.operator import CAP_APPLY, LinearOperator
from orpheus.sn.geometry import SNMesh
from orpheus.sn.quadrature import LebedevSphere
from orpheus.sn.scattering import ScatteringOperator
from orpheus.sn.solver import SNSolver

pytestmark = pytest.mark.foundation  # software-invariant tier


def _uniform_2d(nx, ny, delta, mat_map):
    """Helper: build a uniform Mesh2D."""
    return Mesh2D(
        edges_x=np.linspace(0, nx * delta, nx + 1),
        edges_y=np.linspace(0, ny * delta, ny + 1),
        mat_map=np.asarray(mat_map, dtype=int),
    )


@pytest.fixture
def solver_2g_p0():
    """SNSolver fixture, 2-group, P0 only (no anisotropic data)."""
    fuel = get_mixture("A", "2g")
    mod = get_mixture("B", "2g")
    materials = {2: fuel, 0: mod}

    nx, ny = 6, 4
    delta = 0.2
    mat = np.zeros((nx, ny), dtype=int)
    mat[:3, :] = 2
    mat[3:, :] = 0

    mesh = _uniform_2d(nx, ny, delta, mat)
    quad = LebedevSphere.create(order=17)
    sn_mesh = SNMesh(mesh, quad)
    solver = SNSolver(materials, sn_mesh)
    return solver


# ── Reference implementations (per-cell loops, known correct) ─────────


def _ref_iso_scatter_inplace(solver, Q, phi):
    """Reference per-cell P0 in-scatter (bit-identical to the legacy code)."""
    out = Q.copy()
    nx, ny = solver.sn_mesh.nx, solver.sn_mesh.ny
    for ix in range(nx):
        for iy in range(ny):
            mid = int(solver.sn_mesh.mat_map[ix, iy])
            out[ix, iy, :] += solver.sig_s0[mid].T @ phi[ix, iy, :]
    return out


def _ref_n2n_inplace(solver, Q, phi):
    """Reference per-cell (n,2n) source (bit-identical to the legacy code)."""
    out = Q.copy()
    nx, ny = solver.sn_mesh.nx, solver.sn_mesh.ny
    for ix in range(nx):
        for iy in range(ny):
            mid = int(solver.sn_mesh.mat_map[ix, iy])
            out[ix, iy, :] += 2.0 * (solver.sig2[mid].T @ phi[ix, iy, :])
    return out


# ──────────────────────────────────────────────────────────────────────
# Protocol contract
# ──────────────────────────────────────────────────────────────────────


class TestProtocolCompliance:
    """ScatteringOperator must satisfy the LinearOperator Protocol."""

    def test_implements_linear_operator(self, solver_2g_p0):
        """isinstance LinearOperator (runtime-checkable Protocol)."""
        assert isinstance(solver_2g_p0.scattering_op, LinearOperator)

    def test_capability_set_apply_only(self, solver_2g_p0):
        """capabilities = {"apply"} — no useful inverse, no transpose surface yet."""
        op = solver_2g_p0.scattering_op
        assert op.capabilities == frozenset({CAP_APPLY})
        assert "solve" not in op.capabilities
        assert "apply_transpose" not in op.capabilities

    def test_apply_accepts_psi_shape(self, solver_2g_p0):
        """apply(psi) must accept (N, nx, ny, ng) angular flux."""
        op = solver_2g_p0.scattering_op
        N = op.n_ordinates
        psi = np.ones((N, op.nx, op.ny, op.ng))
        out = op.apply(psi)
        assert out.shape == psi.shape


# ──────────────────────────────────────────────────────────────────────
# Bit-identical extraction (load-bearing)
# ──────────────────────────────────────────────────────────────────────


class TestBitIdenticalExtractionP0:
    """The lifted math must match the legacy reference per-cell code."""

    def test_add_iso_source_matches_reference(self, solver_2g_p0):
        """ScatteringOperator.add_iso_source = the per-cell reference."""
        np.random.seed(42)
        nx, ny, ng = solver_2g_p0.sn_mesh.nx, solver_2g_p0.sn_mesh.ny, solver_2g_p0.ng
        phi = np.random.rand(nx, ny, ng) + 0.1
        Q = np.random.rand(nx, ny, ng)

        expected = _ref_iso_scatter_inplace(solver_2g_p0, Q, phi)

        Q_actual = Q.copy()
        solver_2g_p0.scattering_op.add_iso_source(Q_actual, phi)

        np.testing.assert_allclose(Q_actual, expected, rtol=1e-13)

    def test_add_n2n_source_matches_reference(self, solver_2g_p0):
        """ScatteringOperator.add_n2n_source = the per-cell reference."""
        np.random.seed(123)
        nx, ny, ng = solver_2g_p0.sn_mesh.nx, solver_2g_p0.sn_mesh.ny, solver_2g_p0.ng
        phi = np.random.rand(nx, ny, ng) + 0.1
        Q = np.random.rand(nx, ny, ng)

        expected = _ref_n2n_inplace(solver_2g_p0, Q, phi)

        Q_actual = Q.copy()
        solver_2g_p0.scattering_op.add_n2n_source(Q_actual, phi)

        np.testing.assert_allclose(Q_actual, expected, rtol=1e-13)

    def test_zero_flux_zero_addition(self, solver_2g_p0):
        """φ = 0 => ScatteringOperator adds zero (linearity guard)."""
        nx, ny, ng = solver_2g_p0.sn_mesh.nx, solver_2g_p0.sn_mesh.ny, solver_2g_p0.ng
        Q = np.ones((nx, ny, ng))
        phi = np.zeros_like(Q)
        Q_before = Q.copy()
        solver_2g_p0.scattering_op.add_iso_source(Q, phi)
        np.testing.assert_array_equal(Q, Q_before)

    def test_delegators_match_operator_directly(self, solver_2g_p0):
        """SNSolver._add_scattering_source delegates to op.add_iso_source bit-identically."""
        np.random.seed(7)
        nx, ny, ng = solver_2g_p0.sn_mesh.nx, solver_2g_p0.sn_mesh.ny, solver_2g_p0.ng
        phi = np.random.rand(nx, ny, ng) + 0.1
        Q = np.random.rand(nx, ny, ng)

        Q_via_delegator = Q.copy()
        solver_2g_p0._add_scattering_source(Q_via_delegator, phi)

        Q_via_operator = Q.copy()
        solver_2g_p0.scattering_op.add_iso_source(Q_via_operator, phi)

        np.testing.assert_array_equal(Q_via_delegator, Q_via_operator)

    def test_delegator_n2n_matches_operator_directly(self, solver_2g_p0):
        """SNSolver._add_n2n_source delegates to op.add_n2n_source bit-identically."""
        np.random.seed(11)
        nx, ny, ng = solver_2g_p0.sn_mesh.nx, solver_2g_p0.sn_mesh.ny, solver_2g_p0.ng
        phi = np.random.rand(nx, ny, ng) + 0.1
        Q = np.random.rand(nx, ny, ng)

        Q_via_delegator = Q.copy()
        solver_2g_p0._add_n2n_source(Q_via_delegator, phi)

        Q_via_operator = Q.copy()
        solver_2g_p0.scattering_op.add_n2n_source(Q_via_operator, phi)

        np.testing.assert_array_equal(Q_via_delegator, Q_via_operator)


# ──────────────────────────────────────────────────────────────────────
# Pℓ Galerkin reconstruction
# ──────────────────────────────────────────────────────────────────────


class TestAnisotropicScatteringExtraction:
    """build_aniso_source must reproduce the legacy _build_aniso_scattering."""

    @pytest.fixture
    def solver_2g_p1(self):
        """4-group with P1 anisotropic scattering data."""
        # Use 421-group library which carries P1; if not available, skip.
        try:
            mix = get_mixture("A", "4g")
        except Exception:
            pytest.skip("4g library unavailable")
        if len(mix.SigS) < 2:
            pytest.skip("No P1 data in 4g library")

        mesh = _uniform_2d(2, 2, 0.5, np.zeros((2, 2), dtype=int))
        quad = LebedevSphere.create(order=17)
        return SNSolver({0: mix}, SNMesh(mesh, quad), scattering_order=1)

    def test_returns_none_for_p0(self, solver_2g_p0):
        """L=0 => Pℓ contribution is None (signal: no aniso source needed)."""
        N = solver_2g_p0.quad.N
        nx, ny, ng = solver_2g_p0.sn_mesh.nx, solver_2g_p0.sn_mesh.ny, solver_2g_p0.ng
        psi = np.ones((N, nx, ny, ng))
        out = solver_2g_p0.scattering_op.build_aniso_source(psi)
        assert out is None

    def test_returns_none_for_no_angular_flux(self, solver_2g_p0):
        """First-iteration sentinel: psi=None => return None."""
        out = solver_2g_p0.scattering_op.build_aniso_source(None)
        assert out is None

    def test_isotropic_flux_zero_aniso_source(self, solver_2g_p1):
        """Isotropic ψ_n = const for every ordinate => P1+ Galerkin moments = 0."""
        op = solver_2g_p1.scattering_op
        N = op.n_ordinates
        psi_iso = np.ones((N, op.nx, op.ny, op.ng))
        Q_aniso = op.build_aniso_source(psi_iso)
        assert Q_aniso is not None
        np.testing.assert_allclose(Q_aniso, 0, atol=1e-12)

    def test_delegator_matches_operator(self, solver_2g_p1):
        """SNSolver._build_aniso_scattering delegates bit-identically."""
        op = solver_2g_p1.scattering_op
        N = op.n_ordinates
        np.random.seed(42)
        psi = np.random.rand(N, op.nx, op.ny, op.ng) + 0.1

        out_via_delegator = solver_2g_p1._build_aniso_scattering(psi)
        out_via_operator = op.build_aniso_source(psi)
        np.testing.assert_array_equal(out_via_delegator, out_via_operator)


# ──────────────────────────────────────────────────────────────────────
# apply() — the LinearOperator surface
# ──────────────────────────────────────────────────────────────────────


class TestApplySemantics:
    """apply(psi) returns the per-ordinate scattering source.

    Combines P0 in-scatter + (n,2n) (broadcast across N) + Pℓ (genuine
    per-ordinate) into a single (N, nx, ny, ng) array.
    """

    def test_apply_isotropic_flux_p0_only(self, solver_2g_p0):
        """For P0-only solver, apply(ψ) = (P0 in-scatter + (n,2n))(φ) broadcast."""
        op = solver_2g_p0.scattering_op
        N = op.n_ordinates
        nx, ny, ng = op.nx, op.ny, op.ng

        np.random.seed(5)
        psi = np.random.rand(N, nx, ny, ng) + 0.1

        # Compute scalar flux the same way apply() does internally.
        phi = np.einsum('n,nxyg->xyg', op.weights, psi)

        # Reference: compute Q_iso explicitly
        Q_iso = np.zeros((nx, ny, ng))
        op.add_iso_source(Q_iso, phi)
        op.add_n2n_source(Q_iso, phi)
        expected = np.broadcast_to(Q_iso[None, :, :, :], psi.shape)

        actual = op.apply(psi)
        np.testing.assert_allclose(actual, expected, rtol=1e-13)

    def test_apply_zero_psi_returns_zero(self, solver_2g_p0):
        """ψ = 0 => S·ψ = 0 (linearity guard)."""
        op = solver_2g_p0.scattering_op
        N = op.n_ordinates
        psi = np.zeros((N, op.nx, op.ny, op.ng))
        out = op.apply(psi)
        np.testing.assert_array_equal(out, np.zeros_like(psi))

    def test_apply_linearity(self, solver_2g_p0):
        """S·(αψ_1 + βψ_2) = αS·ψ_1 + βS·ψ_2."""
        op = solver_2g_p0.scattering_op
        N = op.n_ordinates
        nx, ny, ng = op.nx, op.ny, op.ng

        np.random.seed(13)
        psi1 = np.random.rand(N, nx, ny, ng) + 0.1
        psi2 = np.random.rand(N, nx, ny, ng) + 0.1
        alpha, beta = 2.5, -1.7

        lhs = op.apply(alpha * psi1 + beta * psi2)
        rhs = alpha * op.apply(psi1) + beta * op.apply(psi2)
        np.testing.assert_allclose(lhs, rhs, rtol=1e-12, atol=1e-13)


# ──────────────────────────────────────────────────────────────────────
# Algebraic identities (P0 + (n,2n))
# ──────────────────────────────────────────────────────────────────────


class TestP0AlgebraicIdentities:
    """Hand-checkable cases for the P0 + (n,2n) algebra."""

    def test_p0_uniform_flux_homogeneous(self):
        """Homogeneous medium, uniform φ_g = 1: Q_iso[g] = Σ_g' Σ_s0[g'->g]."""
        mix = make_mixture(
            sig_t=np.array([0.5, 1.0]),
            sig_c=np.array([0.01, 0.02]),
            sig_f=np.array([0.01, 0.08]),
            nu=np.array([2.5, 2.5]),
            chi=np.array([1.0, 0.0]),
            sig_s=np.array([[0.38, 0.10], [0.00, 0.90]]),
        )
        nx, ny = 2, 2
        mesh = _uniform_2d(nx, ny, 0.5, np.zeros((nx, ny), dtype=int))
        quad = LebedevSphere.create(order=17)
        solver = SNSolver({0: mix}, SNMesh(mesh, quad))
        op = solver.scattering_op

        phi = np.ones((nx, ny, op.ng))
        Q = np.zeros_like(phi)
        op.add_iso_source(Q, phi)

        # Hand-computed: Q[g] = Σ_g' σ_s0[g'->g] · φ[g'] = column-sum · 1.
        sig_s0_dense = np.array(mix.SigS[0].todense())
        # Convention: ORPHEUS ``SigS[l]`` matrix entry ``[g_from, g_to]``.
        # phi @ sig_s0 sums over g_from for each g_to.
        expected_per_cell = np.ones(op.ng) @ sig_s0_dense
        for ix in range(nx):
            for iy in range(ny):
                np.testing.assert_allclose(Q[ix, iy, :], expected_per_cell, rtol=1e-14)

    def test_n2n_doubling_factor(self):
        """For a pure-(n,2n) mixture (Σ_s0 = 0), Q = 2·φ·Σ_2n."""
        # Build a synthetic mixture with zero P0 scatter and known sig2.
        mix = make_mixture(
            sig_t=np.array([0.5, 1.0]),
            sig_c=np.array([0.01, 0.02]),
            sig_f=np.array([0.0, 0.0]),  # no fission
            nu=np.array([0.0, 0.0]),
            chi=np.array([1.0, 0.0]),
            sig_s=np.array([[0.0, 0.0], [0.0, 0.0]]),  # zero P0
        )
        # Inject a non-zero (n,2n) matrix manually after construction.
        from scipy.sparse import csr_matrix
        sig2_test = np.array([[0.0, 0.05], [0.0, 0.0]])
        mix.Sig2 = csr_matrix(sig2_test)

        nx, ny = 2, 2
        mesh = _uniform_2d(nx, ny, 0.5, np.zeros((nx, ny), dtype=int))
        quad = LebedevSphere.create(order=17)
        solver = SNSolver({0: mix}, SNMesh(mesh, quad))
        op = solver.scattering_op

        np.random.seed(31)
        phi = np.random.rand(nx, ny, op.ng) + 0.1
        Q = np.zeros_like(phi)
        op.add_iso_source(Q, phi)
        # P0 contribution should be zero
        np.testing.assert_allclose(Q, 0, atol=1e-15)

        # (n,2n) contribution
        op.add_n2n_source(Q, phi)
        # Hand-computed: Q[ix, iy, g] = 2 · sum_g' phi[ix, iy, g'] · sig2[g'->g]
        for ix in range(nx):
            for iy in range(ny):
                expected = 2.0 * phi[ix, iy, :] @ sig2_test
                np.testing.assert_allclose(Q[ix, iy, :], expected, rtol=1e-14)
