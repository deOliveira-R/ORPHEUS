r"""Adjoint scattering — Λᵀ leaf + frame-conjugated kernel transpose (campaign #276, A2).

The SN scattering kernel is the frame conjugation :math:`R\circ\Lambda\circ M`
(Funk–Hecke: the SH angular frame is scattering's eigenbasis). Its Euclidean
transpose is :math:`M^{T}\circ\Lambda^{T}\circ R^{T}` — and since the angular
:math:`M`/:math:`R` faces already transpose (the Frame carve, Phase D), the ONLY
genuinely-new piece is :math:`\Lambda^{T}`, the per-ℓ group-axis transpose of the
block-diagonal :math:`\Sigma_{s,\ell}` matmul. Once :math:`\Lambda` advertises
``apply_transpose``, ``(R∘Λ∘M).apply_transpose`` falls out of
:meth:`OperatorProduct.apply_transpose` for free.

This file gates the A2 leaf:

* **Λᵀ** — the moment-space transpose identity ``⟨Λ m, c⟩ = ⟨m, Λᵀ c⟩`` (the
  DEFINING transpose property), a structurally-independent per-material dense
  ``sigᵀ`` reference, the group-flip mutation (Λᵀ ≠ Λ on asymmetric Σ_s), and the
  capability flip.
* **kernel = R∘Λ∘M** — the Euclidean reciprocity ``⟨kernel ψ, c⟩ = ⟨ψ, kernelᵀ c⟩``
  (the aniso transpose now live via the operator algebra) + capability propagation.

The metric-correct Hilbert adjoint ``S† = G⁻¹SᵀG`` is the ``.H`` wrapper's job
(A3/A4); these gates pin the BARE Euclidean transpose (per L27: per-group / full
tensor contraction, never a weight-summed scalar that telescopes).

vv Mode-8: ``np.testing.*`` / :func:`require` only (fire under ``python -O``).
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import Mesh2D
from orpheus.numerics.operator import CAP_APPLY_TRANSPOSE, CAP_SOLVE
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.solver import SNSolver
from orpheus.transport.operators.scattering import LegendreMomentScattering

pytestmark = pytest.mark.foundation


def require(condition: bool, message: str) -> None:
    if not condition:
        pytest.fail(message)


def _uniform_2d(nx, ny, delta, mat_map):
    return Mesh2D(
        edges_x=np.linspace(0, nx * delta, nx + 1),
        edges_y=np.linspace(0, ny * delta, ny + 1),
        mat_map=np.asarray(mat_map, dtype=int),
    )


def _mix(p0, p1):
    m = make_mixture(
        sig_t=np.array([0.5, 1.0]), sig_c=np.array([0.01, 0.02]),
        sig_f=np.array([0.0, 0.0]), nu=np.array([0.0, 0.0]),
        chi=np.zeros(2), sig_s=p0,
    )
    m.SigS = [csr_matrix(p0), csr_matrix(p1)]
    m.Sig2 = csr_matrix(np.array([[0.0, 0.03], [0.01, 0.0]]))
    return m


# ASYMMETRIC P0 + P1 blocks per material (so a group-axis transpose is detectable).
_P0_A = np.array([[0.38, 0.10], [0.05, 0.90]]); _P1_A = np.array([[0.02, 0.01], [0.00, 0.04]])
_P0_B = np.array([[0.55, 0.03], [0.12, 0.40]]); _P1_B = np.array([[0.06, 0.02], [0.01, 0.03]])


@pytest.fixture
def solver_p1_het():
    nx, ny = 4, 3
    mat = np.zeros((nx, ny), dtype=int); mat[:2, :] = 0; mat[2:, :] = 1
    sn_mesh = SNMesh(_uniform_2d(nx, ny, 0.4, mat), Quadrature.lebedev(order=17),
                     {0: _mix(_P0_A, _P1_A), 1: _mix(_P0_B, _P1_B)})
    return SNSolver(sn_mesh, scattering_order=1)


def _moment_field(op, nx, ny, seed):
    return np.random.default_rng(seed).uniform(0.05, 1.0, size=(2, 3, op.ng, nx, ny))


# ═══════════════════════════════════════════════════════════════════════
# Λᵀ — the one genuinely-new leaf (per-ℓ group-axis transpose).
# ═══════════════════════════════════════════════════════════════════════


class TestLambdaTranspose:
    def test_capability_apply_and_transpose_not_solve(self, solver_p1_het):
        lam = LegendreMomentScattering(mat_xs=solver_p1_het.scattering_op.mat_xs, L=1)
        require(CAP_APPLY_TRANSPOSE in lam.capabilities,
                f"Λ must advertise apply_transpose (campaign #276); got {lam.capabilities}.")
        require(CAP_SOLVE not in lam.capabilities,
                "Λ must NOT advertise solve (ℓ=0 block rank-deficient).")

    def test_moment_space_transpose_identity(self, solver_p1_het):
        r"""``⟨Λ m, c⟩ = ⟨m, Λᵀ c⟩`` (full moment-tensor contraction, per L27)."""
        op = solver_p1_het.scattering_op
        nx, ny = op.mat_xs.spatial_shape
        lam = LegendreMomentScattering(mat_xs=op.mat_xs, L=1, skip_l0=False)
        m = _moment_field(op, nx, ny, 1); c = _moment_field(op, nx, ny, 2)
        lhs = float((lam.apply(m) * c).sum())            # ⟨Λ m, c⟩
        rhs = float((m * lam.apply_transpose(c)).sum())  # ⟨m, Λᵀ c⟩
        np.testing.assert_allclose(
            lhs, rhs, rtol=1e-12,
            err_msg="Λ moment-space transpose identity ⟨Λm,c⟩=⟨m,Λᵀc⟩ violated.",
        )

    def test_transpose_matches_dense_per_material(self, solver_p1_het):
        r"""STRUCTURALLY-INDEPENDENT: Λᵀ block = transpose of the forward matrix.

        The forward verb (``einsum("mfc,fg->mgc")``) applies, per (cell, ℓ, m),
        the matrix :math:`A = \Sigma_{s,\ell}^{T}` to the group vector
        (``out_g = Σ_f in_f·sig[f,g]``). The transpose therefore applies
        :math:`A^{T} = \Sigma_{s,\ell}` (un-transposed) — built here by an explicit
        per-(material, ℓ, cell) ``sig @ vec`` Python loop (no einsum). A wrong group
        axis in the production transpose verb disagrees with it.
        """
        op = solver_p1_het.scattering_op
        nx, ny = op.mat_xs.spatial_shape
        lam = LegendreMomentScattering(mat_xs=op.mat_xs, L=1, skip_l0=False)
        c = _moment_field(op, nx, ny, 3)
        got = lam.apply_transpose(c)

        ref = np.zeros_like(c)
        for mid, idx in op.mat_xs.cells_by_material.items():
            sig = op.mat_xs.sig_s_legendre(mid)  # list over ℓ of (ng, ng) [g_from, g_to]
            for l in range(2):
                n_m = 2 * l + 1
                # Forward applies sigᵀ ⇒ transpose applies sig (un-transposed).
                sig_l = np.asarray(sig[l])
                for (ix, iy) in zip(*idx):
                    for mom in range(n_m):
                        ref[l, mom, :, ix, iy] = sig_l @ c[l, mom, :, ix, iy]
        np.testing.assert_allclose(
            got, ref, rtol=1e-12, atol=0.0,
            err_msg="Λᵀ disagrees with the explicit per-material sig_s[ℓ] matmul "
            "(the transpose of the forward sigᵀ) — a wrong group axis in the verb.",
        )

    def test_group_flip_is_nontrivial(self, solver_p1_het):
        r"""Discriminator: with asymmetric Σ_s, Λᵀ ≠ Λ (the transpose has teeth)."""
        op = solver_p1_het.scattering_op
        nx, ny = op.mat_xs.spatial_shape
        lam = LegendreMomentScattering(mat_xs=op.mat_xs, L=1, skip_l0=False)
        m = _moment_field(op, nx, ny, 4)
        require(
            not np.allclose(lam.apply(m), lam.apply_transpose(m)),
            "Λ and Λᵀ agreed on asymmetric Σ_s — the fixture lost its asymmetry, "
            "so the group-flip gate is blind to a transpose error.",
        )


# ═══════════════════════════════════════════════════════════════════════
# kernel = R∘Λ∘M — the aniso transpose, now free via the operator algebra.
# ═══════════════════════════════════════════════════════════════════════


class TestKernelTranspose:
    def test_kernel_advertises_apply_transpose(self, solver_p1_het):
        kernel = solver_p1_het.scattering_op.kernel
        require(
            CAP_APPLY_TRANSPOSE in kernel.capabilities,
            "kernel (R∘Λ∘M) must propagate apply_transpose once Λ has it "
            f"(OperatorProduct intersection); got {kernel.capabilities}.",
        )

    def test_kernel_euclidean_reciprocity(self, solver_p1_het):
        r"""``⟨kernel ψ, c⟩ = ⟨ψ, kernelᵀ c⟩`` — the aniso R∘Λ∘M Euclidean transpose.

        Confirms ``(R∘Λ∘M)ᵀ = Mᵀ∘Λᵀ∘Rᵀ`` composes correctly through
        ``OperatorProduct.apply_transpose`` with the Phase-D M/R face transposes.
        Full per-ordinate/per-group contraction (NOT weight-summed — L27).
        """
        op = solver_p1_het.scattering_op
        nx, ny = op.mat_xs.spatial_shape
        N = op.weights.shape[0]
        kernel = op.kernel
        rng = np.random.default_rng(5)
        psi = rng.uniform(0.05, 1.0, size=(N, op.ng, nx, ny))
        c = rng.uniform(0.05, 1.0, size=(N, op.ng, nx, ny))
        lhs = float((kernel.apply(psi) * c).sum())
        rhs = float((psi * kernel.apply_transpose(c)).sum())
        np.testing.assert_allclose(
            lhs, rhs, rtol=1e-12,
            err_msg="kernel (R∘Λ∘M) Euclidean reciprocity violated — Mᵀ∘Λᵀ∘Rᵀ "
            "did not compose correctly through OperatorProduct.apply_transpose.",
        )
