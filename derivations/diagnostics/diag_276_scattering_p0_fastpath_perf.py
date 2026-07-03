"""Diagnostic: the P0 iso fast-path is meaningfully cheaper than the frame form
(full_scatter_kernel) — the load-bearing perf number for campaign #276 Option-2.

Created by numerics-investigator on 2026-06-28.
This is a PERF characterization, NOT a correctness gate. Promote ONLY if a
perf-regression marker exists; otherwise LEAVE here as the design-decision
evidence of record (a wall-clock test is environment-sensitive and a poor CI
gate). Recommendation: LEAVE (print + soft-assert), reference from the #276
plan / the ScatteringOperator.full_scatter_kernel docstring.

Question answered: does keeping the scalar reaction-rate fast-path for ℓ=0 (the
Option-2 design `R₀∘K∘M₀` composed-but-fast) pay off, or could the forward route
through the frame `full_scatter_kernel` with no perf cost?

Findings (Apple M-series, .venv py3.14, min-of-N-repeats):
  * The TRUE hot path is the in-place SWEEP source (raw-ndarray P0+n2n on the
    scalar flux), used by SNSolver._add_scattering_source per SI iterate.
  * The frame form (even an ℓ=0-only PRE-BUILT R₀∘(Λ₀+N2N)∘M₀, fresh-build cost
    excluded) is 1.3–1.4x the sweep in 1-D and 2.2–2.7x in 2-D, because M/R
    materialise a moment field over all N ordinates then reconstruct, vs the
    fast-path's direct scalar reduction + group-matmul + broadcast.
  * An ℓ=0-only frame is NOT cheaper than the full kernel — the cost is the M/R
    framing over ordinates, not the moment-tensor width.
  ⟹ The fast-path is justified; Option-2 (keep scalar ℓ=0 fast-path, sum with
    the aniso frame.conjugate(Λ_{ℓ≥1})) is the right design.

vv Mode-8: pytest.fail (fires under -O); the perf assertion is a SOFT lower
bound (skipped if timing is too noisy), never a hard wall-clock equality.
"""
from __future__ import annotations

import timeit

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.solver import SNSolver
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.operators.scattering import (
    LegendreMomentScattering, N2NMomentOperator,
)
from orpheus.transport.frames.harmonic_frame import HarmonicFrame

pytestmark = pytest.mark.foundation


def _mix(p0, p1, sig_t, ng, n2n):
    m = make_mixture(
        sig_t=sig_t, sig_c=np.full(ng, 0.01), sig_f=np.zeros(ng),
        nu=np.zeros(ng), chi=np.zeros(ng), sig_s=p0,
    )
    m.SigS = [csr_matrix(p0), csr_matrix(p1)]
    m.Sig2 = csr_matrix(n2n)
    return m


def _blocks(ng, seed):
    rng = np.random.default_rng(seed)
    p0 = rng.uniform(0.02, 0.4, size=(ng, ng))
    p1 = rng.uniform(0.0, 0.05, size=(ng, ng))
    n2n = rng.uniform(0.0, 0.02, size=(ng, ng)); np.fill_diagonal(n2n, 0.0)
    sig_t = rng.uniform(0.5, 1.5, size=ng)
    return p0, p1, n2n, sig_t


def _build_2d_p0(nx, ny, ng, seed=0):
    mat = np.zeros((nx, ny), dtype=int); mat[nx // 2:, :] = 1
    mesh = Mesh2D(
        edges_x=np.linspace(0, nx * 0.1, nx + 1),
        edges_y=np.linspace(0, ny * 0.1, ny + 1), mat_map=mat,
    )
    quad = Quadrature.product(n_mu=8, n_phi=8)
    mats = {}
    for mid in (0, 1):
        p0, p1, n2n, sig_t = _blocks(ng, seed + mid)
        mats[mid] = _mix(p0, p1, sig_t, ng, n2n)
    sn = SNMesh(mesh, quad, mats)
    op = SNSolver(sn, scattering_order=0).scattering_op
    vals = np.random.default_rng(seed + 99).uniform(
        0.05, 1.0, size=(quad.N, ng, *sn.spatial_shape))
    return op, AngularFlux.from_mesh(vals, sn), sn


def _ell0_kernel(op):
    """The Option-2 ℓ=0-only PRE-BUILT composition R₀∘(Λ₀+N2N)∘M₀."""
    frame0 = HarmonicFrame.from_galerkin(op.quadrature.angular_frame(0))
    lam0 = LegendreMomentScattering(mat_xs=op.mat_xs, L=0, skip_l0=False)
    n2n0 = N2NMomentOperator(mat_xs=op.mat_xs, L=0)
    return frame0.conjugate(lam0 + n2n0)


def _t(fn, n_calls=300, n_repeat=7):
    fn()  # warm
    return min(timeit.repeat(fn, number=n_calls, repeat=n_repeat)) / n_calls


def test_p0_sweep_fastpath_beats_frame_form_2d():
    """The in-place P0 sweep source is materially faster than the (pre-built)
    ℓ=0 frame form in 2-D — the perf basis for keeping the fast-path.

    Soft gate: asserts the frame form is at least 1.5x the sweep in 2-D (the
    measured factor is ~2.2-2.7x; 1.5x leaves generous margin for noisy CI).
    If timing is too noisy (the sweep itself > 5 ms/call) the assertion is
    skipped and the numbers are reported for the record.
    """
    op, psi, sn = _build_2d_p0(nx=20, ny=20, ng=2)
    W = float(op.weights.sum())
    spatial = sn.spatial_shape
    ng = op.ng

    def sweep():
        phi = psi.integrate_angular().values
        Q = np.zeros((ng, *spatial))
        op.add_iso_source(Q, phi)
        op.add_n2n_source(Q, phi)
        return Q / W

    k0 = _ell0_kernel(op)
    k0.apply(psi.values)

    def frame():
        return k0.apply(psi.values) / W

    # correctness precondition: same value (P0 is 0-ULP)
    A = op.apply(psi).values
    B = np.asarray(k0.apply(psi.values)) / W
    np.testing.assert_allclose(
        A, B, rtol=1e-12, atol=1e-14,
        err_msg="ℓ=0 frame form must reproduce the fast-path forward (P0).",
    )

    t_sweep = _t(sweep)
    t_frame = _t(frame)
    ratio = t_frame / t_sweep
    print(
        f"\n[P0 2D 20x20 ng=2 prod8x8 N={op.n_ordinates}] "
        f"sweep={t_sweep*1e9:.0f} ns  frame(ℓ=0 prebuilt)={t_frame*1e9:.0f} ns  "
        f"ratio={ratio:.2f}x"
    )

    if t_sweep > 5e-3:  # too slow / noisy environment — report only
        pytest.skip(f"timing too noisy (sweep {t_sweep*1e3:.1f} ms/call); ratio={ratio:.2f}")
    if ratio < 1.5:
        pytest.fail(
            f"frame form only {ratio:.2f}x the fast-path in 2-D (expected >=1.5x; "
            f"measured ~2.2-2.7x). The Option-2 perf premise — keep the scalar "
            f"fast-path for ℓ=0 because the frame form is materially slower — no "
            f"longer holds; re-evaluate whether the forward should route through "
            f"the frame."
        )


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v", "-s"]))
