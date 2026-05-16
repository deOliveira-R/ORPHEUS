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
    sn_mesh = SNMesh(mesh, quad, materials)
    solver = SNSolver(sn_mesh)
    return solver


# ── Reference implementations (per-cell loops, known correct) ─────────


def _ref_iso_scatter_inplace(solver, Q, phi):
    """Reference per-cell P0 in-scatter (bit-identical to the legacy code).

    Issue #196 PR-INDEX-4: ``Q`` / ``phi`` are principled
    ``(ng, nx, ny)``.  Per-cell update reads ``(ng,)`` slices over the
    spatial pair.
    """
    out = Q.copy()
    nx, ny = solver.sn_mesh.nx, solver.sn_mesh.ny
    for ix in range(nx):
        for iy in range(ny):
            mid = int(solver.sn_mesh.mat_map[ix, iy])
            out[:, ix, iy] += solver.sig_s0[mid].T @ phi[:, ix, iy]
    return out


def _ref_n2n_inplace(solver, Q, phi):
    """Reference per-cell (n,2n) source (bit-identical to the legacy code).

    Issue #196 PR-INDEX-4 — principled ``(ng, nx, ny)`` (see
    :func:`_ref_iso_scatter_inplace`).
    """
    out = Q.copy()
    nx, ny = solver.sn_mesh.nx, solver.sn_mesh.ny
    for ix in range(nx):
        for iy in range(ny):
            mid = int(solver.sn_mesh.mat_map[ix, iy])
            out[:, ix, iy] += 2.0 * (solver.sig2[mid].T @ phi[:, ix, iy])
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
        """apply(psi) must accept (N, ng, nx, ny) angular flux (Issue #196 PR-INDEX-4)."""
        op = solver_2g_p0.scattering_op
        N = op.n_ordinates
        psi = np.ones((N, op.ng, op.nx, op.ny))
        out = op.apply(psi)
        assert out.shape == psi.shape


# ──────────────────────────────────────────────────────────────────────
# Bit-identical extraction (load-bearing)
# ──────────────────────────────────────────────────────────────────────


class TestBitIdenticalExtractionP0:
    """The lifted math must match the legacy reference per-cell code."""

    def test_add_iso_source_matches_reference(self, solver_2g_p0):
        """ScatteringOperator.add_iso_source = the per-cell reference.

        Issue #196 PR-INDEX-4: principled ``(ng, nx, ny)`` end-to-end.
        """
        np.random.seed(42)
        nx, ny, ng = solver_2g_p0.sn_mesh.nx, solver_2g_p0.sn_mesh.ny, solver_2g_p0.ng
        phi = np.random.rand(ng, nx, ny) + 0.1
        Q = np.random.rand(ng, nx, ny)

        expected = _ref_iso_scatter_inplace(solver_2g_p0, Q, phi)

        Q_actual = Q.copy()
        solver_2g_p0.scattering_op.add_iso_source(Q_actual, phi)

        np.testing.assert_allclose(Q_actual, expected, rtol=1e-13)

    def test_add_n2n_source_matches_reference(self, solver_2g_p0):
        """ScatteringOperator.add_n2n_source = the per-cell reference."""
        np.random.seed(123)
        nx, ny, ng = solver_2g_p0.sn_mesh.nx, solver_2g_p0.sn_mesh.ny, solver_2g_p0.ng
        phi = np.random.rand(ng, nx, ny) + 0.1
        Q = np.random.rand(ng, nx, ny)

        expected = _ref_n2n_inplace(solver_2g_p0, Q, phi)

        Q_actual = Q.copy()
        solver_2g_p0.scattering_op.add_n2n_source(Q_actual, phi)

        np.testing.assert_allclose(Q_actual, expected, rtol=1e-13)

    def test_zero_flux_zero_addition(self, solver_2g_p0):
        """φ = 0 => ScatteringOperator adds zero (linearity guard)."""
        nx, ny, ng = solver_2g_p0.sn_mesh.nx, solver_2g_p0.sn_mesh.ny, solver_2g_p0.ng
        Q = np.ones((ng, nx, ny))
        phi = np.zeros_like(Q)
        Q_before = Q.copy()
        solver_2g_p0.scattering_op.add_iso_source(Q, phi)
        np.testing.assert_array_equal(Q, Q_before)

    def test_delegators_match_operator_directly(self, solver_2g_p0):
        """SNSolver._add_scattering_source delegates to op.add_iso_source bit-identically.

        Issue #196 PR-INDEX-5: delegator's PUBLIC contract is now
        principled ``(ng, nx, ny)``.  No bridge.
        """
        np.random.seed(7)
        nx, ny, ng = solver_2g_p0.sn_mesh.nx, solver_2g_p0.sn_mesh.ny, solver_2g_p0.ng
        phi = np.random.rand(ng, nx, ny) + 0.1
        Q_init = np.random.rand(ng, nx, ny)

        Q_via_delegator = Q_init.copy()
        solver_2g_p0._add_scattering_source(Q_via_delegator, phi)

        Q_via_operator = Q_init.copy()
        solver_2g_p0.scattering_op.add_iso_source(Q_via_operator, phi)

        np.testing.assert_array_equal(Q_via_delegator, Q_via_operator)

    def test_delegator_n2n_matches_operator_directly(self, solver_2g_p0):
        """SNSolver._add_n2n_source delegates to op.add_n2n_source bit-identically.

        Issue #196 PR-INDEX-5: principled ``(ng, nx, ny)`` end-to-end.
        """
        np.random.seed(11)
        nx, ny, ng = solver_2g_p0.sn_mesh.nx, solver_2g_p0.sn_mesh.ny, solver_2g_p0.ng
        phi = np.random.rand(ng, nx, ny) + 0.1
        Q_init = np.random.rand(ng, nx, ny)

        Q_via_delegator = Q_init.copy()
        solver_2g_p0._add_n2n_source(Q_via_delegator, phi)

        Q_via_operator = Q_init.copy()
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
        return SNSolver(SNMesh(mesh, quad, {0: mix}), scattering_order=1)

    def test_returns_none_for_p0(self, solver_2g_p0):
        """L=0 => Pℓ contribution is None (signal: no aniso source needed)."""
        N = solver_2g_p0.quad.N
        nx, ny, ng = solver_2g_p0.sn_mesh.nx, solver_2g_p0.sn_mesh.ny, solver_2g_p0.ng
        # Issue #196 PR-INDEX-4: principled (N, ng, nx, ny).
        psi = np.ones((N, ng, nx, ny))
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
        psi_iso = np.ones((N, op.ng, op.nx, op.ny))
        Q_aniso = op.build_aniso_source(psi_iso)
        assert Q_aniso is not None
        np.testing.assert_allclose(Q_aniso, 0, atol=1e-12)

    def test_delegator_matches_operator(self, solver_2g_p1):
        """SNSolver._build_aniso_scattering delegates bit-identically.

        Issue #196 PR-INDEX-5: delegator and operator both consume /
        return principled ``(N, ng, nx, ny)``.
        """
        op = solver_2g_p1.scattering_op
        N = op.n_ordinates
        np.random.seed(42)
        psi = np.random.rand(N, op.ng, op.nx, op.ny) + 0.1

        out_via_delegator = solver_2g_p1._build_aniso_scattering(psi)
        out_via_operator = op.build_aniso_source(psi)
        np.testing.assert_array_equal(out_via_delegator, out_via_operator)


# ──────────────────────────────────────────────────────────────────────
# apply() — the LinearOperator surface
# ──────────────────────────────────────────────────────────────────────


class TestApplySemantics:
    """apply(psi) returns the per-ordinate scattering source.

    Combines P0 in-scatter + (n,2n) (broadcast across N) + Pℓ (genuine
    per-ordinate) into a single ``(N, ng, nx, ny)`` array (principled
    storage; see :ref:`theory-sn-index-convention`).
    """

    def test_apply_isotropic_flux_p0_only(self, solver_2g_p0):
        """For P0-only solver, apply(ψ) = (P0 in-scatter + (n,2n))(φ) broadcast.

        Issue #196 PR-INDEX-4: principled ``(N, ng, nx, ny)`` ψ.
        """
        op = solver_2g_p0.scattering_op
        N = op.n_ordinates
        nx, ny, ng = op.nx, op.ny, op.ng

        np.random.seed(5)
        psi = np.random.rand(N, ng, nx, ny) + 0.1

        # Compute scalar flux the same way apply() does internally.
        phi = np.einsum('n,ngxy->gxy', op.weights, psi)

        # Reference: compute Q_iso explicitly
        Q_iso = np.zeros((ng, nx, ny))
        op.add_iso_source(Q_iso, phi)
        op.add_n2n_source(Q_iso, phi)
        expected = np.broadcast_to(Q_iso[None, :, :, :], psi.shape)

        actual = op.apply(psi)
        np.testing.assert_allclose(actual, expected, rtol=1e-13)

    def test_apply_zero_psi_returns_zero(self, solver_2g_p0):
        """ψ = 0 => S·ψ = 0 (linearity guard)."""
        op = solver_2g_p0.scattering_op
        N = op.n_ordinates
        psi = np.zeros((N, op.ng, op.nx, op.ny))
        out = op.apply(psi)
        np.testing.assert_array_equal(out, np.zeros_like(psi))

    def test_apply_linearity(self, solver_2g_p0):
        """S·(αψ_1 + βψ_2) = αS·ψ_1 + βS·ψ_2."""
        op = solver_2g_p0.scattering_op
        N = op.n_ordinates
        nx, ny, ng = op.nx, op.ny, op.ng

        np.random.seed(13)
        psi1 = np.random.rand(N, ng, nx, ny) + 0.1
        psi2 = np.random.rand(N, ng, nx, ny) + 0.1
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
        """Homogeneous medium, uniform φ_g = 1: Q_iso[g] = Σ_g' Σ_s0[g'->g].

        Issue #196 PR-INDEX-4: principled ``(ng, nx, ny)`` ψ / Q.
        """
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
        solver = SNSolver(SNMesh(mesh, quad, {0: mix}))
        op = solver.scattering_op

        phi = np.ones((op.ng, nx, ny))
        Q = np.zeros_like(phi)
        op.add_iso_source(Q, phi)

        # Hand-computed: Q[g] = Σ_g' σ_s0[g'->g] · φ[g'] = column-sum · 1.
        sig_s0_dense = np.array(mix.SigS[0].todense())
        # Convention: ORPHEUS ``SigS[l]`` matrix entry ``[g_from, g_to]``.
        # phi @ sig_s0 sums over g_from for each g_to.
        expected_per_cell = np.ones(op.ng) @ sig_s0_dense
        for ix in range(nx):
            for iy in range(ny):
                np.testing.assert_allclose(Q[:, ix, iy], expected_per_cell, rtol=1e-14)

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
        solver = SNSolver(SNMesh(mesh, quad, {0: mix}))
        op = solver.scattering_op

        np.random.seed(31)
        # Issue #196 PR-INDEX-4: principled (ng, nx, ny).
        phi = np.random.rand(op.ng, nx, ny) + 0.1
        Q = np.zeros_like(phi)
        op.add_iso_source(Q, phi)
        # P0 contribution should be zero
        np.testing.assert_allclose(Q, 0, atol=1e-15)

        # (n,2n) contribution
        op.add_n2n_source(Q, phi)
        # Hand-computed: Q[g, ix, iy] = 2 · sum_g' phi[g', ix, iy] · sig2[g'->g]
        for ix in range(nx):
            for iy in range(ny):
                expected = 2.0 * phi[:, ix, iy] @ sig2_test
                np.testing.assert_allclose(Q[:, ix, iy], expected, rtol=1e-14)


# ──────────────────────────────────────────────────────────────────────
# Foldable / residual split — Phase G Step 3+4.a (Issue #196)
# ──────────────────────────────────────────────────────────────────────


@pytest.fixture
def solver_2g_p1_n2n():
    """2G solver with non-trivial cross-group P0 AND a Pℓ≥1 channel AND
    a non-zero (n,2n) matrix. Stresses every channel of the residual."""
    # Asymmetric P0 matrix: non-trivial diagonal AND off-diagonal entries.
    p0 = np.array([[0.38, 0.10], [0.05, 0.90]])
    # Non-trivial P1 block — Pℓ≥1 is unconditionally residual.
    p1 = np.array([[0.02, 0.01], [0.00, 0.04]])
    from scipy.sparse import csr_matrix
    mix = make_mixture(
        sig_t=np.array([0.5, 1.0]),
        sig_c=np.array([0.01, 0.02]),
        sig_f=np.array([0.01, 0.08]),
        nu=np.array([2.5, 2.5]),
        chi=np.array([1.0, 0.0]),
        sig_s=p0,
    )
    # Append a P1 block manually; ORPHEUS's SigS is a list[csr_matrix]
    # indexed by Legendre order.
    mix.SigS = [csr_matrix(p0), csr_matrix(p1)]
    # Inject (n,2n) — non-zero on a cross-group entry only (the brief
    # explicitly notes diagonal sig2 entries are rare but legal).
    mix.Sig2 = csr_matrix(np.array([[0.0, 0.03], [0.01, 0.0]]))

    nx, ny = 3, 2
    mesh = _uniform_2d(nx, ny, 0.4, np.zeros((nx, ny), dtype=int))
    quad = LebedevSphere.create(order=17)
    return SNSolver(SNMesh(mesh, quad, {0: mix}), scattering_order=1)


class TestFoldablePart:
    """``foldable_part()`` carries ONLY the P0 within-group diagonal."""

    def test_returns_scattering_operator_instance(self, solver_2g_p0):
        """Mechanism criterion 1 — sibling class, not a new class."""
        S = solver_2g_p0.scattering_op
        assert isinstance(S.foldable_part(), ScatteringOperator)

    def test_scattering_order_is_zero(self, solver_2g_p0):
        """Mechanism criterion 2 — no Pℓ structure in foldable."""
        S = solver_2g_p0.scattering_op
        assert S.foldable_part().scattering_order == 0

    def test_Y_is_None(self, solver_2g_p0):
        """Mechanism criterion 3 — no spherical harmonics for ℓ=0."""
        S = solver_2g_p0.scattering_op
        assert S.foldable_part().Y is None

    def test_Y_is_None_even_for_p1_source(self, solver_2g_p1_n2n):
        """Even when S carries P1+ data, the foldable sibling has no Y."""
        S = solver_2g_p1_n2n.scattering_op
        assert S.foldable_part().Y is None

    def test_sig_s_is_diagonal_only(self, solver_2g_p1_n2n):
        """Mechanism criterion 4a — sig_s[mid][0] is diagonal-only."""
        S = solver_2g_p1_n2n.scattering_op
        S_fold = S.foldable_part()
        for mid in S.sig_s:
            mat = S_fold.sig_s[mid][0]
            expected = np.diag(np.diag(S.sig_s[mid][0]))
            np.testing.assert_array_equal(mat, expected)
            # Off-diagonal is literally zero, not just small.
            off_diag = mat - np.diag(np.diag(mat))
            assert np.all(off_diag == 0.0)

    def test_sig_s0_matches_sig_s_l0(self, solver_2g_p1_n2n):
        """Mechanism criterion 4b — sig_s0 == sig_s[mid][0]."""
        S = solver_2g_p1_n2n.scattering_op
        S_fold = S.foldable_part()
        for mid in S.sig_s:
            np.testing.assert_array_equal(
                S_fold.sig_s0[mid], S_fold.sig_s[mid][0]
            )

    def test_sig_s_has_length_one(self, solver_2g_p1_n2n):
        """Mechanism criterion 4c — no Pℓ≥1 entries in foldable."""
        S = solver_2g_p1_n2n.scattering_op
        S_fold = S.foldable_part()
        for mid in S.sig_s:
            assert len(S_fold.sig_s[mid]) == 1

    def test_sig2_is_zero_matrix(self, solver_2g_p1_n2n):
        """Mechanism criterion 4d — (n,2n) belongs to residual unconditionally."""
        S = solver_2g_p1_n2n.scattering_op
        S_fold = S.foldable_part()
        for mid in S.sig2:
            assert S_fold.sig2[mid].shape == S.sig2[mid].shape
            assert S_fold.sig2[mid].dtype == S.sig2[mid].dtype
            np.testing.assert_array_equal(
                S_fold.sig2[mid], np.zeros_like(S.sig2[mid])
            )

    def test_does_not_mutate_parent_sig_s(self, solver_2g_p1_n2n):
        """Anti-rec 4 — split returns new arrays; parent unchanged."""
        S = solver_2g_p1_n2n.scattering_op
        # Snapshot every parent array.
        before = {mid: [m.copy() for m in S.sig_s[mid]] for mid in S.sig_s}
        before_sig2 = {mid: S.sig2[mid].copy() for mid in S.sig2}
        _ = S.foldable_part()
        # Parent unchanged.
        for mid in S.sig_s:
            for l, m in enumerate(S.sig_s[mid]):
                np.testing.assert_array_equal(m, before[mid][l])
        for mid in S.sig2:
            np.testing.assert_array_equal(S.sig2[mid], before_sig2[mid])


class TestResidualPart:
    """``residual_part()`` carries everything but P0 within-group diagonal."""

    def test_returns_scattering_operator_instance(self, solver_2g_p0):
        """Mechanism criterion 1 — sibling class."""
        S = solver_2g_p0.scattering_op
        assert isinstance(S.residual_part(), ScatteringOperator)

    def test_sig_s_l0_diagonal_zeroed(self, solver_2g_p1_n2n):
        """Mechanism criterion 5a — cross-group only on P0."""
        S = solver_2g_p1_n2n.scattering_op
        S_res = S.residual_part()
        for mid in S.sig_s:
            expected = S.sig_s[mid][0] - np.diag(np.diag(S.sig_s[mid][0]))
            np.testing.assert_array_equal(S_res.sig_s[mid][0], expected)
            # The diagonal IS zero, not just close.
            diag = np.diag(S_res.sig_s[mid][0])
            assert np.all(diag == 0.0)

    def test_sig_s0_matches_diagonal_zeroed(self, solver_2g_p1_n2n):
        """Mechanism criterion 5b — sig_s0 alias of sig_s[mid][0]."""
        S = solver_2g_p1_n2n.scattering_op
        S_res = S.residual_part()
        for mid in S.sig_s:
            np.testing.assert_array_equal(
                S_res.sig_s0[mid], S_res.sig_s[mid][0]
            )

    def test_pl_ge_1_carried_verbatim(self, solver_2g_p1_n2n):
        """Mechanism criterion 5c — Pℓ≥1 blocks unchanged."""
        S = solver_2g_p1_n2n.scattering_op
        assert S.scattering_order >= 1, "fixture must carry P1+ data"
        S_res = S.residual_part()
        for mid in S.sig_s:
            for l in range(1, S.scattering_order + 1):
                np.testing.assert_array_equal(
                    S_res.sig_s[mid][l], S.sig_s[mid][l]
                )

    def test_sig2_carried_verbatim(self, solver_2g_p1_n2n):
        """Mechanism criterion 5d — (n,2n) unconditionally residual."""
        S = solver_2g_p1_n2n.scattering_op
        S_res = S.residual_part()
        for mid in S.sig2:
            np.testing.assert_array_equal(S_res.sig2[mid], S.sig2[mid])

    def test_scattering_order_preserved(self, solver_2g_p1_n2n):
        """Mechanism criterion 5e — Pℓ structure preserved."""
        S = solver_2g_p1_n2n.scattering_op
        assert S.residual_part().scattering_order == S.scattering_order

    def test_Y_is_self_Y(self, solver_2g_p1_n2n):
        """Mechanism criterion 5f — precomputed harmonics reusable."""
        S = solver_2g_p1_n2n.scattering_op
        S_res = S.residual_part()
        # Either same object (preferred) or array-equal.
        assert S_res.Y is S.Y or np.array_equal(S_res.Y, S.Y)

    def test_Y_None_for_p0_solver(self, solver_2g_p0):
        """If S has no harmonics (L=0), residual has none either."""
        S = solver_2g_p0.scattering_op
        assert S.residual_part().Y is None

    def test_does_not_mutate_parent_sig_s(self, solver_2g_p1_n2n):
        """Anti-rec 4 — split returns new arrays; parent unchanged."""
        S = solver_2g_p1_n2n.scattering_op
        before = {mid: [m.copy() for m in S.sig_s[mid]] for mid in S.sig_s}
        _ = S.residual_part()
        for mid in S.sig_s:
            for l, m in enumerate(S.sig_s[mid]):
                np.testing.assert_array_equal(m, before[mid][l])


class TestFoldableSigma:
    """``foldable_sigma()`` returns the per-material (ng,) σ_{s,0}^{g→g}."""

    def test_returns_dict_of_ndarrays(self, solver_2g_p1_n2n):
        """Mechanism criterion 6a — dict[int, ndarray]."""
        S = solver_2g_p1_n2n.scattering_op
        result = S.foldable_sigma()
        assert isinstance(result, dict)
        for mid, arr in result.items():
            assert isinstance(mid, int)
            assert isinstance(arr, np.ndarray)

    def test_shape_is_ng(self, solver_2g_p1_n2n):
        """Mechanism criterion 6b — each value is (ng,)."""
        S = solver_2g_p1_n2n.scattering_op
        result = S.foldable_sigma()
        for arr in result.values():
            assert arr.shape == (S.ng,)

    def test_values_are_diagonal_of_sig_s0(self, solver_2g_p1_n2n):
        """Mechanism criterion 6c — equals np.diag(sig_s[mid][0])."""
        S = solver_2g_p1_n2n.scattering_op
        result = S.foldable_sigma()
        for mid, arr in result.items():
            np.testing.assert_array_equal(arr, np.diag(S.sig_s[mid][0]))

    def test_returned_arrays_are_copies(self, solver_2g_p1_n2n):
        """Mutating the returned dict's values must not affect ``self``."""
        S = solver_2g_p1_n2n.scattering_op
        result = S.foldable_sigma()
        # Snapshot parent diagonal.
        before = {mid: np.diag(S.sig_s[mid][0]).copy() for mid in S.sig_s}
        # Mutate the returned arrays.
        for arr in result.values():
            arr[:] = -999.0
        # Parent unchanged.
        for mid in S.sig_s:
            np.testing.assert_array_equal(
                np.diag(S.sig_s[mid][0]), before[mid]
            )


class TestAlgebraicIdentity:
    """The load-bearing contract:
    ``S.apply(ψ) ≈ S.foldable_part().apply(ψ) + S.residual_part().apply(ψ)``
    at ``rtol=1e-14`` (FP-non-associativity precision).

    Covers P0-only, Pℓ≥1, non-zero (n,2n), and cross-group + diagonal
    coupling — the four cases enumerated in the brief's criterion 7.
    """

    def _check_identity(self, op, psi):
        full = op.apply(psi)
        split_sum = op.foldable_part().apply(psi) + op.residual_part().apply(psi)
        np.testing.assert_allclose(full, split_sum, rtol=1e-14, atol=1e-15)

    def test_identity_p0_only_random_psi(self, solver_2g_p0):
        """Case 1 — scattering_order == 0 only (no Pℓ).

        Issue #196 PR-INDEX-4: principled ``(N, ng, nx, ny)`` ψ.
        """
        op = solver_2g_p0.scattering_op
        assert op.scattering_order == 0
        N = op.n_ordinates
        np.random.seed(42)
        psi = np.random.rand(N, op.ng, op.nx, op.ny) + 0.1
        self._check_identity(op, psi)

    def test_identity_p0_only_uniform_psi(self, solver_2g_p0):
        """Case 1b — uniform ψ probes the diagonal isolation path."""
        op = solver_2g_p0.scattering_op
        N = op.n_ordinates
        psi = np.ones((N, op.ng, op.nx, op.ny))
        self._check_identity(op, psi)

    def test_identity_with_pl_ge_1(self, solver_2g_p1_n2n):
        """Case 2 — scattering_order >= 1 (with non-zero P1 block)."""
        op = solver_2g_p1_n2n.scattering_op
        assert op.scattering_order >= 1
        N = op.n_ordinates
        np.random.seed(101)
        psi = np.random.rand(N, op.ng, op.nx, op.ny) + 0.1
        self._check_identity(op, psi)

    def test_identity_with_nonzero_n2n(self, solver_2g_p1_n2n):
        """Case 3 — non-zero (n,2n) coupling."""
        op = solver_2g_p1_n2n.scattering_op
        # Fixture explicitly sets (n,2n) cross-group entries.
        any_nonzero_n2n = any(
            np.any(op.sig2[mid] != 0.0) for mid in op.sig2
        )
        assert any_nonzero_n2n, "fixture must carry non-zero (n,2n)"
        N = op.n_ordinates
        np.random.seed(202)
        psi = np.random.rand(N, op.ng, op.nx, op.ny) + 0.1
        self._check_identity(op, psi)

    def test_identity_multigroup_cross_group_plus_diagonal(self, solver_2g_p1_n2n):
        """Case 4 — non-trivial cross-group P0 + diagonal coupling."""
        op = solver_2g_p1_n2n.scattering_op
        # Fixture's P0 matrix has both diagonal AND off-diagonal entries.
        for mid in op.sig_s:
            p0 = op.sig_s[mid][0]
            diag = np.diag(p0)
            off = p0 - np.diag(diag)
            assert np.any(diag != 0.0)
            assert np.any(off != 0.0)
        N = op.n_ordinates
        np.random.seed(303)
        psi = np.random.rand(N, op.ng, op.nx, op.ny) + 0.1
        self._check_identity(op, psi)

    def test_residual_zero_when_p0_diagonal_only_no_n2n(self):
        """Pure-diagonal P0 with no (n,2n) and no Pℓ≥1 ⇒ residual.apply(ψ)=0
        and full == foldable.apply(ψ) by construction."""
        from scipy.sparse import csr_matrix
        # Strictly diagonal P0, zero (n,2n).
        mix = make_mixture(
            sig_t=np.array([0.5, 1.0]),
            sig_c=np.array([0.01, 0.02]),
            sig_f=np.array([0.0, 0.0]),
            nu=np.array([0.0, 0.0]),
            chi=np.array([1.0, 0.0]),
            sig_s=np.diag([0.3, 0.8]),
        )
        mix.Sig2 = csr_matrix(np.zeros((2, 2)))
        nx, ny = 2, 2
        mesh = _uniform_2d(nx, ny, 0.5, np.zeros((nx, ny), dtype=int))
        quad = LebedevSphere.create(order=17)
        op = SNSolver(SNMesh(mesh, quad, {0: mix})).scattering_op

        N = op.n_ordinates
        np.random.seed(404)
        # Issue #196 PR-INDEX-4: principled (N, ng, nx, ny).
        psi = np.random.rand(N, op.ng, nx, ny) + 0.1
        full = op.apply(psi)
        residual_part = op.residual_part().apply(psi)
        np.testing.assert_allclose(residual_part, 0.0, atol=1e-15)
        # And full ≡ foldable up to FP-non-associativity.
        foldable_part = op.foldable_part().apply(psi)
        np.testing.assert_allclose(full, foldable_part, rtol=1e-14, atol=1e-15)


class TestPurity:
    """``foldable_part()`` / ``residual_part()`` are pure functions —
    calling twice returns instances with equal per-material arrays
    (mechanism criterion 8)."""

    def test_foldable_part_pure(self, solver_2g_p1_n2n):
        S = solver_2g_p1_n2n.scattering_op
        a, b = S.foldable_part(), S.foldable_part()
        assert a.scattering_order == b.scattering_order
        assert a.Y is None and b.Y is None
        for mid in S.sig_s:
            np.testing.assert_array_equal(a.sig_s[mid][0], b.sig_s[mid][0])
            np.testing.assert_array_equal(a.sig_s0[mid], b.sig_s0[mid])
            np.testing.assert_array_equal(a.sig2[mid], b.sig2[mid])

    def test_residual_part_pure(self, solver_2g_p1_n2n):
        S = solver_2g_p1_n2n.scattering_op
        a, b = S.residual_part(), S.residual_part()
        assert a.scattering_order == b.scattering_order
        for mid in S.sig_s:
            for l in range(S.scattering_order + 1):
                np.testing.assert_array_equal(
                    a.sig_s[mid][l], b.sig_s[mid][l]
                )
            np.testing.assert_array_equal(a.sig2[mid], b.sig2[mid])

    def test_foldable_sigma_pure(self, solver_2g_p1_n2n):
        S = solver_2g_p1_n2n.scattering_op
        a, b = S.foldable_sigma(), S.foldable_sigma()
        assert set(a.keys()) == set(b.keys())
        for mid in a:
            np.testing.assert_array_equal(a[mid], b[mid])


# ──────────────────────────────────────────────────────────────────────
# is_foldable_into_sigma_r — Phase G Step 3+4.b.i (Issue #196)
# ──────────────────────────────────────────────────────────────────────


class TestIsFoldableIntoSigmaR:
    """``S.is_foldable_into_sigma_r()`` returns True iff S carries only
    diagonal P0 + zero sig2.

    Consumed by substep 3+4.b.ii's ``OperatorSum.solve`` fusion hook
    to detect "this S is the foldable_part — fuse into σ_r and route
    to the within-group sweep". A STRUCTURAL predicate on the
    operator's data, not an identity claim about its action.
    """

    def test_full_scattering_returns_false(self, solver_2g_p1_n2n):
        """Full S with non-zero off-diagonal P0 + non-zero P1 + non-zero
        sig2 → NOT foldable."""
        S = solver_2g_p1_n2n.scattering_op
        # Sanity: the fixture's S has all three non-foldable channels.
        assert S.scattering_order >= 1
        assert S.is_foldable_into_sigma_r() is False

    def test_foldable_part_roundtrip_is_true(self, solver_2g_p1_n2n):
        """``S.foldable_part().is_foldable_into_sigma_r() == True``.

        The load-bearing round-trip: the operator constructed by
        ``foldable_part()`` IS, by definition, the foldable part of
        itself.
        """
        S = solver_2g_p1_n2n.scattering_op
        foldable = S.foldable_part()
        assert foldable.is_foldable_into_sigma_r() is True

    def test_residual_part_returns_false(self, solver_2g_p1_n2n):
        """``S.residual_part().is_foldable_into_sigma_r() == False``.

        The residual carries the cross-group off-diagonal P0
        unconditionally (every multi-group system has at least one
        cross-group entry) — so the diagonal-only check fails.
        """
        S = solver_2g_p1_n2n.scattering_op
        residual = S.residual_part()
        assert residual.is_foldable_into_sigma_r() is False

    def test_p0_only_diagonal_returns_true(self):
        """P0-only ScatteringOperator with diagonal sig_s + zero sig2 →
        True (positive control).

        Build a synthetic ScatteringOperator directly (bypassing
        SNSolver) to isolate the predicate from any fixture setup.
        """
        ng = 2
        p0_diag = np.diag([0.38, 0.90])
        S = ScatteringOperator(
            n_ordinates=12,
            nx=2, ny=2, ng=ng,
            scattering_order=0,
            sig_s={0: [p0_diag]},
            sig2={0: np.zeros((ng, ng))},
            sig_s0={0: p0_diag},
            Y=None,
            weights=np.ones(12) / 12.0,
            cells_by_mat={0: (
                np.zeros(4, dtype=int),
                np.zeros(4, dtype=int),
            )},
        )
        assert S.is_foldable_into_sigma_r() is True

    def test_p0_with_off_diagonal_returns_false(self):
        """scattering_order=0 with non-diagonal P0 → False.

        Off-diagonal P0 is cross-group scattering — couples distinct
        energy groups and cannot collapse into a per-cell scalar.
        """
        ng = 2
        # Non-diagonal P0 — non-zero off-diagonal entry.
        p0 = np.array([[0.38, 0.10], [0.00, 0.90]])
        S = ScatteringOperator(
            n_ordinates=12,
            nx=2, ny=2, ng=ng,
            scattering_order=0,
            sig_s={0: [p0]},
            sig2={0: np.zeros((ng, ng))},
            sig_s0={0: p0},
            Y=None,
            weights=np.ones(12) / 12.0,
            cells_by_mat={0: (
                np.zeros(4, dtype=int),
                np.zeros(4, dtype=int),
            )},
        )
        assert S.is_foldable_into_sigma_r() is False

    def test_p0_diagonal_with_nonzero_sig2_returns_false(self):
        """scattering_order=0 with diagonal P0 BUT non-zero sig2 → False.

        (n,2n) doubling is unconditionally residual: folding into a
        "removal" cross-section is conceptually wrong because (n,2n)
        emits two neutrons per absorption.
        """
        ng = 2
        p0_diag = np.diag([0.38, 0.90])
        sig2 = np.array([[0.0, 0.05], [0.0, 0.0]])
        S = ScatteringOperator(
            n_ordinates=12,
            nx=2, ny=2, ng=ng,
            scattering_order=0,
            sig_s={0: [p0_diag]},
            sig2={0: sig2},
            sig_s0={0: p0_diag},
            Y=None,
            weights=np.ones(12) / 12.0,
            cells_by_mat={0: (
                np.zeros(4, dtype=int),
                np.zeros(4, dtype=int),
            )},
        )
        assert S.is_foldable_into_sigma_r() is False

    def test_scattering_order_ge_1_returns_false_even_with_diagonal_p0(
        self,
    ):
        """scattering_order >= 1 → False even if P0 is diagonal.

        Pℓ ≥ 1 is direction-dependent (Y_ℓ^m(Ω_n)) — unconditionally
        residual. The presence of ANY Pℓ ≥ 1 channel disqualifies the
        operator from foldability.
        """
        ng = 2
        p0_diag = np.diag([0.38, 0.90])
        p1 = np.array([[0.02, 0.00], [0.00, 0.04]])
        S = ScatteringOperator(
            n_ordinates=12,
            nx=2, ny=2, ng=ng,
            scattering_order=1,
            sig_s={0: [p0_diag, p1]},
            sig2={0: np.zeros((ng, ng))},
            sig_s0={0: p0_diag},
            Y=None,  # Not used by the predicate.
            weights=np.ones(12) / 12.0,
            cells_by_mat={0: (
                np.zeros(4, dtype=int),
                np.zeros(4, dtype=int),
            )},
        )
        assert S.is_foldable_into_sigma_r() is False
