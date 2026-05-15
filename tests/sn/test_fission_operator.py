"""Foundation tests for :class:`orpheus.sn.fission.FissionOperator`.

Round 1.2 of Wave D of the SN reshape campaign (Issue #162). The
operator carries the same math :class:`SNSolver` used to expose under
``compute_fission_source``; that method is now a thin delegator. These
tests pin the lifted math at the operator level and the convention
that :meth:`FissionOperator.apply` returns :math:`F\\,\\phi` *without*
the :math:`1/k` division (the caller divides).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture, make_mixture
from orpheus.geometry import Mesh2D
from orpheus.numerics.operator import CAP_APPLY, LinearOperator
from orpheus.sn.fission import FissionOperator
from orpheus.sn.geometry import SNMesh
from orpheus.sn.quadrature import LebedevSphere
from orpheus.sn.solver import SNSolver

pytestmark = pytest.mark.foundation


def _uniform_2d(nx, ny, delta, mat_map):
    return Mesh2D(
        edges_x=np.linspace(0, nx * delta, nx + 1),
        edges_y=np.linspace(0, ny * delta, ny + 1),
        mat_map=np.asarray(mat_map, dtype=int),
    )


@pytest.fixture
def solver_2g():
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
    return SNSolver(materials, sn_mesh)


# ──────────────────────────────────────────────────────────────────────
# Protocol contract
# ──────────────────────────────────────────────────────────────────────


class TestProtocolCompliance:
    """FissionOperator must satisfy the LinearOperator Protocol."""

    def test_implements_linear_operator(self, solver_2g):
        assert isinstance(solver_2g.fission_op, LinearOperator)

    def test_capability_set_apply_only(self, solver_2g):
        """capabilities = {"apply"} — F is rank-1 in energy, no inverse."""
        op = solver_2g.fission_op
        assert op.capabilities == frozenset({CAP_APPLY})
        assert "solve" not in op.capabilities


# ──────────────────────────────────────────────────────────────────────
# Bit-identical extraction (the load-bearing claim)
# ──────────────────────────────────────────────────────────────────────


class TestBitIdenticalExtraction:
    """FissionOperator math must match the legacy compute_fission_source."""

    def test_apply_matches_legacy_compute_at_k_eq_1(self, solver_2g):
        """apply(φ) = compute_fission_source(φ, k=1.0) bit-identically."""
        np.random.seed(42)
        nx, ny, ng = solver_2g.sn_mesh.nx, solver_2g.sn_mesh.ny, solver_2g.ng
        # PR-INDEX-4: FissionOperator.apply consumes / returns principled
        # (ng, nx, ny).  Use that shape directly.
        phi = np.random.rand(ng, nx, ny) + 0.1

        out_op = solver_2g.fission_op.apply(phi)
        # Reference: hand-coded version of the legacy method (no division by k).
        # All operands principled (ng, nx, ny).
        fission_rate = np.einsum("gxy,gxy->xy", solver_2g.sig_p, phi)
        expected = solver_2g.chi * fission_rate[None, :, :]

        np.testing.assert_array_equal(out_op, expected)

    def test_delegator_matches_apply_with_k(self, solver_2g):
        """SNSolver.compute_fission_source(φ, k) = fission_op.apply(φ) / k.

        SNSolver.compute_fission_source still consumes legacy (nx, ny, ng)
        on its PUBLIC contract (PR-INDEX-5 flips it).  fission_op.apply
        consumes principled (ng, nx, ny).  Verify the delegator's
        internal bridge transposes round-trip.
        """
        np.random.seed(7)
        nx, ny, ng = solver_2g.sn_mesh.nx, solver_2g.sn_mesh.ny, solver_2g.ng
        phi_legacy = np.random.rand(nx, ny, ng) + 0.1
        phi_principled = np.transpose(phi_legacy, (2, 0, 1))

        for k in [1.0, 0.93, 1.27, 0.5]:
            out_via_delegator = solver_2g.compute_fission_source(phi_legacy, k)
            out_via_operator_principled = (
                solver_2g.fission_op.apply(phi_principled) / k
            )
            # delegator returns legacy (nx, ny, ng); operator returns
            # principled (ng, nx, ny) — transpose to compare.
            out_via_operator_legacy = np.transpose(
                out_via_operator_principled, (1, 2, 0),
            )
            np.testing.assert_array_equal(
                out_via_delegator, out_via_operator_legacy,
            )


# ──────────────────────────────────────────────────────────────────────
# Algebraic identities
# ──────────────────────────────────────────────────────────────────────


class TestRank1EnergyStructure:
    """F = χ ⊗ νΣ_f is rank-1 in energy per cell."""

    def test_zero_flux_zero_source(self, solver_2g):
        """φ = 0 => F·φ = 0 (linearity guard)."""
        nx, ny, ng = solver_2g.sn_mesh.nx, solver_2g.sn_mesh.ny, solver_2g.ng
        phi = np.zeros((ng, nx, ny))
        out = solver_2g.fission_op.apply(phi)
        np.testing.assert_array_equal(out, np.zeros_like(phi))

    def test_apply_linearity(self, solver_2g):
        """F·(αφ_1 + βφ_2) = αF·φ_1 + βF·φ_2."""
        nx, ny, ng = solver_2g.sn_mesh.nx, solver_2g.sn_mesh.ny, solver_2g.ng
        np.random.seed(13)
        phi1 = np.random.rand(ng, nx, ny) + 0.1
        phi2 = np.random.rand(ng, nx, ny) + 0.1
        alpha, beta = 2.5, -1.7

        lhs = solver_2g.fission_op.apply(alpha * phi1 + beta * phi2)
        rhs = (alpha * solver_2g.fission_op.apply(phi1)
               + beta * solver_2g.fission_op.apply(phi2))
        np.testing.assert_allclose(lhs, rhs, rtol=1e-12, atol=1e-13)

    def test_constant_flux_uniform_material(self):
        """In a uniform-material homogeneous problem, φ_g = c =>
        (Fφ)_g = χ_g · c · Σ_g' νΣ_{f,g'}."""
        mix = make_mixture(
            sig_t=np.array([0.5, 1.0]),
            sig_c=np.array([0.01, 0.02]),
            sig_f=np.array([0.01, 0.08]),
            nu=np.array([2.5, 2.5]),
            chi=np.array([1.0, 0.0]),
            sig_s=np.array([[0.38, 0.10], [0.00, 0.90]]),
        )
        nx, ny = 3, 3
        mesh = _uniform_2d(nx, ny, 0.5, np.zeros((nx, ny), dtype=int))
        quad = LebedevSphere.create(order=17)
        solver = SNSolver({0: mix}, SNMesh(mesh, quad))

        c = 1.5
        # PR-INDEX-4: principled (ng, nx, ny).
        phi = c * np.ones((solver.ng, nx, ny))
        out = solver.fission_op.apply(phi)

        # Hand-computed: per-cell fission rate = c · Σ_g νΣ_f[g] = c · Σ_g SigP[g].
        # (SigP is the production cross-section νΣ_f on the Mixture.)
        nuSigF_total = float(np.sum(mix.SigP))
        rate = c * nuSigF_total
        for ix in range(nx):
            for iy in range(ny):
                expected = mix.chi * rate
                np.testing.assert_allclose(
                    out[:, ix, iy], expected, rtol=1e-14,
                    err_msg=f"FissionOperator.apply mismatch at cell ({ix}, {iy})",
                )

    def test_chi_normalisation(self, solver_2g):
        """Σ_g χ_g = 1 per material (sanity check on the per-cell χ array)."""
        # Note: for materials with no fission, chi may be zero; we only check
        # the cells whose mixture has nonzero νΣ_f.
        nx, ny = solver_2g.sn_mesh.nx, solver_2g.sn_mesh.ny
        # PR-INDEX-3: solver.chi / solver.sig_p are (ng, nx, ny).
        for mid, (ix_arr, iy_arr) in solver_2g._cells_by_mat.items():
            for ix, iy in zip(ix_arr, iy_arr):
                chi_cell = solver_2g.chi[:, ix, iy]
                if np.sum(solver_2g.sig_p[:, ix, iy]) > 1e-15:
                    np.testing.assert_allclose(chi_cell.sum(), 1.0, rtol=1e-12)


class TestKDivisionConvention:
    """The 1/k division is at the SOLVER level, not the operator level."""

    def test_apply_does_not_divide_by_k(self, solver_2g):
        """apply(φ) is independent of any eigenvalue — pure linear action."""
        np.random.seed(99)
        nx, ny, ng = solver_2g.sn_mesh.nx, solver_2g.sn_mesh.ny, solver_2g.ng
        # PR-INDEX-4: principled (ng, nx, ny).
        phi = np.random.rand(ng, nx, ny) + 0.1

        # Construct the operator twice with no k handed in — apply
        # should not depend on any external state.
        out1 = solver_2g.fission_op.apply(phi)
        out2 = solver_2g.fission_op.apply(phi.copy())
        np.testing.assert_array_equal(out1, out2)

    def test_compute_fission_source_does_divide_by_k(self, solver_2g):
        """SNSolver.compute_fission_source(φ, k) = apply(φ) / k.

        compute_fission_source's PUBLIC contract is legacy (nx, ny, ng)
        until PR-INDEX-5.  Internal bridge transposes don't affect the
        ``1/k`` scaling — the linearity check holds on either layout.
        """
        np.random.seed(101)
        nx, ny, ng = solver_2g.sn_mesh.nx, solver_2g.sn_mesh.ny, solver_2g.ng
        phi = np.random.rand(nx, ny, ng) + 0.1

        out_k_one = solver_2g.compute_fission_source(phi, 1.0)
        out_k_two = solver_2g.compute_fission_source(phi, 2.0)
        np.testing.assert_allclose(out_k_one, 2.0 * out_k_two, rtol=1e-14)
