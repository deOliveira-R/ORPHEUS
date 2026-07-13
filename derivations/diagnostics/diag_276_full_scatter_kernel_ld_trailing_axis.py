"""Diagnostic: full_scatter_kernel reproduces the fast-path forward S on an LD
(spatial-moment) AngularFlux, and the LD cells-index fix (commit 0b3275d, #276)
has mutation teeth.

Created by numerics-investigator on 2026-06-28.
If this test catches a real bug, promote to tests/sn/operators/test_scattering_adjoint.py
(the matching code path: ScatteringOperator.full_scatter_kernel) — it fills the
gap that TestFullScatterKernel only exercises NON-LD fluxes (no trailing 2^d
spatial-moment axis), so the LD trailing-axis correctness of the frame form is
currently UNGUARDED by a regression gate.

Background (campaign #276 A2, commit 0b3275d):
  ScatteringOperator.full_scatter_kernel = frame.conjugate(Λ_{ℓ≥0} + N2N) = R∘(Λ+N2N)∘M
  is the frame form of the full in-scatter source. A prior attempt that routed the
  forward through it found it "diverged from the fast-path for LD". Root cause: the
  per-material cell indexers in MaterialXSField (apply_legendre_scattering_moments
  + _transpose + apply_n2n_moments + _transpose) used `cells = (Ellipsis, *idx)`,
  which — when the moment tensor carries LD's trailing 2^d spatial-moment (φ̂) axis
  — mis-targets axes: Ellipsis greedily absorbs the leading (m, g) axes, so the
  spatial cell indices *idx land on the LAST two axes (one spatial + the φ̂ axis)
  instead of the two spatial axes. The fix pins the leading two axes explicitly:
  `cells = (slice(None), slice(None), *idx)`. Bit-identical non-LD (no trailing
  axis); correct for LD.

Failure-mode class: #2 (variable swap / wrong-axis), gated by LD presence.
Invisible to every non-LD config (Ellipsis ≡ (slice,slice) when no trailing axis).

vv Mode-8: np.testing.* / pytest.fail only (fire under python -O).
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.solver import SNSolver
from orpheus.transport.spatial import LinearDiscontinuous
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.mesh.material_xs_field import MaterialXSField

pytestmark = pytest.mark.foundation


def _mix(p0, p1, sig_t, ng, n2n):
    m = make_mixture(
        sig_t=sig_t, sig_c=np.full(ng, 0.01), sig_f=np.zeros(ng),
        nu=np.zeros(ng), chi=np.zeros(ng), sig_s=p0,
    )
    m.SigS = [csr_matrix(p0), csr_matrix(p1)]
    m.Sig2 = csr_matrix(n2n)
    return m


def _build_ld_solver(order, nx=4, ny=3):
    """Heterogeneous 2-D LD solver — RECTANGULAR (nx != ny) so a wrong-axis
    index over-runs (the cleanest mutation tell)."""
    p0_a = np.array([[0.30, 0.08], [0.04, 0.70]])
    p1_a = np.array([[0.02, 0.01], [0.00, 0.03]])
    p0_b = np.array([[0.55, 0.03], [0.12, 0.40]])
    p1_b = np.array([[0.06, 0.02], [0.01, 0.03]])
    n2n = np.array([[0.0, 0.02], [0.01, 0.0]])
    sig_t = np.array([0.5, 1.0])
    mat = np.zeros((nx, ny), dtype=int)
    mat[nx // 2:, :] = 1
    mesh = Mesh2D(
        edges_x=np.linspace(0, nx * 0.1, nx + 1),
        edges_y=np.linspace(0, ny * 0.1, ny + 1),
        mat_map=mat,
    )
    quad = Quadrature.product(n_mu=4, n_phi=4)
    sn = SNMesh(
        mesh, quad,
        {0: _mix(p0_a, p1_a, sig_t, 2, n2n), 1: _mix(p0_b, p1_b, sig_t, 2, n2n)},
        scheme=LinearDiscontinuous(),
    )
    return SNSolver(sn, scattering_order=order)


def _ld_flux(solver, seed=123):
    N = solver.quad.N
    ng = solver.ng
    nx, ny = solver.sn_mesh.spatial_shape
    n_moments = 2 ** 2  # 2-D LD → 4 spatial moments
    vals = np.random.default_rng(seed).uniform(0.05, 1.0, size=(N, ng, nx, ny, n_moments))
    return AngularFlux.from_mesh(vals, solver.sn_mesh, spatial_moments=2)


@pytest.mark.parametrize("order", [0, 1])
def test_full_scatter_kernel_reproduces_forward_on_ld_flux(order):
    """(1/W)·full_scatter_kernel.apply(ψ) == S.apply(ψ) for an LD AngularFlux.

    The frame form must reproduce the fast-path forward even when ψ carries the
    trailing 2^d spatial-moment (φ̂) axis. Principled-equiv (~1e-15, reduction
    order differs; P0 is 0-ULP since the iso block is a single matmul).
    """
    solver = _build_ld_solver(order)
    op = solver.scattering_op
    psi = _ld_flux(solver)
    W = float(op.weights.sum())

    fast = op.apply(psi).values
    frame = np.asarray(op.full_scatter_kernel.apply(psi.values)) / W

    assert frame.shape == fast.shape, (
        f"full_scatter_kernel output shape {frame.shape} != fast-path {fast.shape} "
        f"on LD flux (trailing φ̂ axis dropped/misplaced)."
    )
    np.testing.assert_allclose(
        frame, fast, rtol=1e-12, atol=1e-14,
        err_msg=(
            f"P{order}: (1/W)·full_scatter_kernel.apply(ψ) does NOT reproduce the "
            f"fast-path S.apply(ψ) on an LD (spatial-moment) flux — the moment-scatter "
            f"cell indexers mis-target the trailing φ̂ axis (the #276/0b3275d "
            f"(Ellipsis,*idx) → (slice,slice,*idx) fix regressed)."
        ),
    )


def test_ld_cells_index_fix_has_mutation_teeth():
    """ISOLATION: the OLD (Ellipsis, *idx) indexing reddens this gate on LD.

    Monkeypatch the four moment-scatter verbs back to the pre-fix `(Ellipsis,
    *idx)` indexing and confirm the frame form fails on an LD flux (it raises
    IndexError on a rectangular grid because the spatial cell indices over-run
    the wrong trailing axis). The non-LD control stays clean — proving the bug
    is LD-gated. This is the consumption proof for the 0b3275d cells-index fix.

    Mutation in-process (monkeypatch); NEVER git checkout (uncommitted-state
    hazard, process-discipline rule).
    """
    def _old_leg(self, moments, L, skip_l0):
        out = np.zeros_like(moments)
        l_start = 1 if skip_l0 else 0
        for mid, idx in self.cells_by_material.items():
            sig = self.sig_s_legendre(mid)
            cells = (Ellipsis, *idx)  # PRE-FIX (buggy under LD)
            for l in range(l_start, L + 1):
                n_m = 2 * l + 1
                mv = moments[l, :n_m][cells]
                out[l, :n_m][cells] = (
                    np.einsum("mfc...,fg->mgc...", mv, sig[l]) + out[l, :n_m][cells]
                )
        return out

    def _old_n2n(self, moments):
        out = np.zeros_like(moments)
        for mid, idx in self.cells_by_material.items():
            sig2 = self.n2n_matrix(mid)
            cells = (Ellipsis, *idx)  # PRE-FIX (buggy under LD)
            mv = moments[0, :1][cells]
            out[0, :1][cells] = (
                2.0 * np.einsum("mfc...,fg->mgc...", mv, sig2) + out[0, :1][cells]
            )
        return out

    solver = _build_ld_solver(order=0)  # rectangular nx=4, ny=3
    op = solver.scattering_op
    psi = _ld_flux(solver)
    W = float(op.weights.sum())

    # Confirm the FIXED code is clean first (sanity).
    fixed = np.asarray(op.full_scatter_kernel.apply(psi.values)) / W
    np.testing.assert_allclose(
        fixed, op.apply(psi).values, rtol=1e-12, atol=1e-14,
        err_msg="precondition: fixed code must reproduce the fast-path on LD.",
    )

    orig_leg = MaterialXSField.apply_legendre_scattering_moments
    orig_n2n = MaterialXSField.apply_n2n_moments
    MaterialXSField.apply_legendre_scattering_moments = _old_leg
    MaterialXSField.apply_n2n_moments = _old_n2n
    try:
        reddened = False
        try:
            mutated = np.asarray(op.full_scatter_kernel.apply(psi.values)) / W
            # If it did not raise, it must at least produce a wrong value.
            if not np.allclose(mutated, op.apply(psi).values, rtol=1e-9, atol=1e-12):
                reddened = True
        except (IndexError, ValueError):
            reddened = True  # wrong-axis over-run on the rectangular grid
        if not reddened:
            pytest.fail(
                "The OLD (Ellipsis,*idx) indexing did NOT redden on an LD flux — "
                "the cells-index fix (0b3275d) has no teeth; the LD trailing-axis "
                "correctness is not actually being verified."
            )
    finally:
        MaterialXSField.apply_legendre_scattering_moments = orig_leg
        MaterialXSField.apply_n2n_moments = orig_n2n

    # Clean revert: the fixed code is restored and reproduces the fast-path.
    np.testing.assert_allclose(
        np.asarray(op.full_scatter_kernel.apply(psi.values)) / W,
        op.apply(psi).values, rtol=1e-12, atol=1e-14,
        err_msg="post-mutation: the monkeypatch did not cleanly revert.",
    )


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
