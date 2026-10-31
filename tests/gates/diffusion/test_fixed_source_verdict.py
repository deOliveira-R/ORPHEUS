r"""Diffusion under the Fredholm verdict (U0): the fixed-source entry, the eigen entry and the kernel bases.

Spec ``scratch/membership_covector/u0/spec_u0.md`` rung 3 (D1-D8) and the bases rows P6-P8. The loss of a closed slab
with an absorption-free energy class is singular; the verdict on ``(A, q)`` is Solve (trivial kernel; the LU's
reciprocal-condition estimate recorded), Gauge (every annihilator pairs to zero with q; the resolvent inverts the
deflated loss ``A + Pi`` and returns the minimum-G-norm member) or ``NoSolution`` before any solve.

Before U0 (`[M]` 2026-10-31, the membership worktree, ``scratch/membership_covector/u0/probes/``): the inner LU
returned max|x| 1.44e17 on a positive source over the closed pure scatterer, silently; ``solve_diffusion_1d`` on
``FIS_FEEDS`` returned k = 74.39 after ``max_outer`` ran out.

Every reference is a dense oracle the code under test does not use: the SVD of the assembled loss, and the bordered
system ``[[A, R], [(G R)^T, 0]]`` for the minimum-G-norm member.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.diffusion.augmented_mesh import DiffusionMesh
from orpheus.diffusion.solver import DiffusionSolver, solve_diffusion_1d, solve_diffusion_fixed_source
from orpheus.numerics.evidence import Measured
from orpheus.numerics.flat_operator import FlattenedOperator
from orpheus.numerics.fredholm import Gauge, NoSolution, Solve
from orpheus.transport.boundary_conditions import WallCondition

from tests.gates.transport._singular_loss_fixtures import (
    ABSORBER, BOTH_2G, FIS_FEEDS, FIS_ISOLATED, FIS_ISOLATED_KINF, FIS_REGULAR, FISSILE_NO_CAPTURE, G2_ONLY,
    G2_UPSCATTER, NEAR_1E8, PURE_1G, TWO_ZONE, bordered_member, dense_of, null, slab,
)

_FRAGMENT = "no steady solution"
_EPS = np.finfo(float).eps

pytestmark = [pytest.mark.l1, pytest.mark.rests_on(
    "tests/gates/transport/test_absorbing_chain_predictor.py",
    "tests/gates/numerics/test_fredholm_verdict.py",
    "tests/gates/numerics/test_matrix_inverse_operator.py::test_the_reciprocal_condition_is_an_estimate_from_above_within_a_factor_ten",
)]


def _dm(mixtures, left=None, right=None, bcs=None, **kw):
    return DiffusionMesh(slab(left, right, **kw), mixtures, **({} if bcs is None else {"boundary_conditions": bcs}))


def _bulk_shape(dm):
    return DiffusionSolver(dm).template.interior.values.shape


def _dense(solver):
    return dense_of(FlattenedOperator(solver.loss, solver.template).apply, solver._n_flat)


def _columns(summands, n):
    cols = []
    for s in summands:
        for j in range(np.asarray(s.generators).shape[1]):
            v = np.zeros(n)
            v[np.asarray(s.support)] = np.asarray(s.generators, float)[:, j]
            cols.append(v)
    return np.column_stack(cols) if cols else np.zeros((n, 0))


def _residual(span, vecs):
    Q, sv, _ = np.linalg.svd(span, full_matrices=False)
    Q = Q[:, sv > 1e-10 * sv[0]]
    return float(np.linalg.norm(vecs - Q @ (Q.T @ vecs)) / np.linalg.norm(vecs))


# ── P6-P8: the bases the method mints from its own loss ──────────────────────────────────────────────────

@pytest.mark.parametrize("mix", [PURE_1G, BOTH_2G, G2_ONLY], ids=["pure1g", "both2g", "g2_only"])
def test_P6_the_predicted_bases_span_the_dense_kernels(mix):
    """REFERENCE: dense SVD, both directions. `[M]` the left co-vector is NOT the integral of the class: its trace rows
    read (0, 0.8, 0, 0.8) relative to the bulk (cos(w, flat(unit)) = 0.937), and g2_only's feeder group reads 0.5."""
    s = DiffusionSolver(_dm({0: mix}))
    A = _dense(s)
    L, R = null(A)
    kernel = s.predicted_kernel
    right, left = _columns(kernel.right, A.shape[0]), _columns(kernel.left, A.shape[0])
    assert R.shape[1] == right.shape[1] == left.shape[1] == 1
    for label, a, b in (("right", right, R), ("left", left, L)):
        assert _residual(a, b) <= 1e-10 and _residual(b, a) <= 1e-10, label


def test_P7_the_left_covector_on_the_trace_is_not_the_integral():
    s = DiffusionSolver(_dm({0: PURE_1G}))
    nb = s.template.interior.values.size
    (w,) = _columns(s.predicted_kernel.left, s._n_flat).T
    np.testing.assert_allclose(w[nb:] / w[0], [0.0, 0.8, 0.0, 0.8], atol=1e-12)


@pytest.mark.parametrize("mix", [PURE_1G, BOTH_2G, G2_ONLY], ids=["pure1g", "both2g", "g2_only"])
def test_P8_the_right_basis_is_annihilated_by_the_loss(mix):
    """Object leg (a count is blind to a kernel rotation): |A r| <= 64 eps ||A||_1 ||r||_inf."""
    s = DiffusionSolver(_dm({0: mix}))
    A = _dense(s)
    right = _columns(s.predicted_kernel.right, A.shape[0])
    assert np.abs(A @ right).max() <= 64 * _EPS * np.abs(A).sum(0).max() * np.abs(right).max()


# ── D1-D5: the fixed-source entry ────────────────────────────────────────────────────────────────────────

@pytest.mark.catches("ERR-123")
@pytest.mark.parametrize("mix", [PURE_1G, BOTH_2G, G2_ONLY], ids=["pure1g", "both2g", "g2_only"])
def test_D1_a_positive_source_on_a_closed_absorption_free_class_is_refused(mix):
    """`[M]` before U0 the inner LU returned max|x| 1.44e17 (pure1g), silently. Cosines 0.929 / 0.929 / 0.881."""
    dm = _dm({0: mix})
    with pytest.raises(NoSolution, match=_FRAGMENT):
        solve_diffusion_fixed_source(dm, np.ones(_bulk_shape(dm)))


def test_D2_a_balanced_source_is_gauged_to_the_minimum_G_norm_member():
    """Two-zone mesh (bulk G = 1.0 | 1.4), a random source with its left-co-vector component removed. REFERENCE: the
    dense bordered solve with (G r)^T x = 0. `[M]` G-projection agrees to 5.8e-15; a Euclidean projector is off by
    3.35e-2 (the non-uniform G is what makes this row discriminate; on the uniform mesh with an antisymmetric source
    the two coincide, 3.8e-16 each)."""
    dm = _dm({0: PURE_1G}, **TWO_ZONE)
    s = DiffusionSolver(dm)
    n, nb = s._n_flat, s.template.interior.values.size
    A = _dense(s)
    L, R = null(A)
    q = np.zeros(n)
    q[:nb] = np.random.default_rng(7).standard_normal(nb)
    q[:nb] -= L[:nb, 0] * (L[:nb, 0] @ q[:nb]) / (L[:nb, 0] @ L[:nb, 0])
    G = np.asarray(dm.full_field_space.metric_diagonal(), dtype=float)
    x_ref = bordered_member(A, q, R, (G[:, None] * R).T)
    sol = solve_diffusion_fixed_source(dm, q[:nb].reshape(_bulk_shape(dm)))
    assert isinstance(sol.exit_report.solvability, Gauge)
    x = np.asarray(sol.flux.to_flat(), dtype=float)
    kappa = np.linalg.cond(A + (R @ (G[:, None] * R).T) / float(R[:, 0] @ (G * R[:, 0])))
    np.testing.assert_allclose(x, x_ref, rtol=0, atol=64 * _EPS * kappa * np.abs(x_ref).max())


def test_D3_the_gauge_summand_spans_bulk_and_trace():
    """The #666 composite-gauge xfail's production replacement: the class's summand crosses the FullFieldSpace block
    boundary, and its projector is idempotent and G-self-adjoint on the full field."""
    dm = _dm({0: PURE_1G})
    s = DiffusionSolver(dm)
    nb, n = s.template.interior.values.size, s._n_flat
    (summand,) = s.predicted_kernel.right
    support = np.asarray(summand.support)
    assert support.min() < nb <= support.max()
    P = s.predicted_kernel.projector
    G = np.asarray(dm.full_field_space.metric_diagonal(), dtype=float)
    x, y = np.random.default_rng(3).standard_normal((2, n))
    np.testing.assert_allclose(P.apply(P.apply(x)), P.apply(x), rtol=0, atol=64 * _EPS * np.abs(x).max())
    assert float(np.sum(G * P.apply(x) * y)) == pytest.approx(float(np.sum(G * x * P.apply(y))), rel=1e-12)


@pytest.mark.parametrize("kw", [
    dict(mixtures={0: G2_UPSCATTER}),
    dict(mixtures={0: PURE_1G}, left=0, right=0, bcs={0: WallCondition.white(0.999)}),
    dict(mixtures={0: ABSORBER, 1: PURE_1G}, bps=(0.0, 5.0, 10.0), mats=(0, 1), n=(4, 4)),
    dict(mixtures={0: FISSILE_NO_CAPTURE}),
], ids=["upscatter", "white0.999", "one_absorbing_region", "fissile_no_capture"])
def test_D4_the_admitted_controls_are_solved(kw):
    dm = _dm(**kw)
    sol = solve_diffusion_fixed_source(dm, np.ones(_bulk_shape(dm)))
    assert isinstance(sol.exit_report.solvability, Solve)


def test_D5_a_near_singular_loss_is_solved_and_its_rcond_is_on_the_exit_report():
    """phi = q / Sigma_a within 10 kappa eps (`[M]` 2.19e-9 against the band 1.11e-7); the Solve verdict carries the
    LU's rcond, at least 1/cond_1 and at most 10/cond_1; admissibility keeps its own meaning (no subcriticality)."""
    dm = _dm({0: NEAR_1E8})
    s = DiffusionSolver(dm)
    sol = solve_diffusion_fixed_source(dm, np.ones(_bulk_shape(dm)))
    kappa = 0.5 / 1e-8
    np.testing.assert_allclose(np.ravel(sol.flux.interior.values), 1.0 / 1e-8, rtol=10 * kappa * _EPS)
    verdict = sol.exit_report.solvability
    assert isinstance(verdict, Solve) and isinstance(verdict.conditioning, Measured)
    true = 1.0 / np.linalg.cond(_dense(s), 1)
    assert true * (1 - max(64 * _EPS, 10 * kappa * _EPS)) <= verdict.conditioning.value <= 10 * true
    assert type(sol.exit_report.admissibility).__name__ == "NotApplicable"


# ── D6-D8: the eigen entry (q = the fission emission) ────────────────────────────────────────────────────

@pytest.mark.catches("ERR-123")
def test_D6_the_eigen_entry_refuses_emission_into_the_class():
    """`[M]` before U0: k = 74.39 after max_outer = 200 ran out."""
    with pytest.raises(NoSolution, match=_FRAGMENT):
        solve_diffusion_1d({0: FIS_FEEDS}, slab(), max_outer=200)


def test_D7_the_eigen_entry_solves_a_regular_fissile_medium():
    """`[M]` k = 0.6250000000000001 before and after; the verdict is Solve."""
    r = solve_diffusion_1d({0: FIS_REGULAR}, slab())
    assert r.keff == pytest.approx(0.625, rel=1e-12)
    assert isinstance(r.solvability, Solve)


def test_D8_the_eigen_entry_gauges_an_isolated_class():
    """Closed form k = nu Sf1 / (St1 - Ss11) = 1.25; before U0 the same value came by pivot luck, with no verdict."""
    r = solve_diffusion_1d({0: FIS_ISOLATED}, slab())
    assert r.keff == pytest.approx(FIS_ISOLATED_KINF, rel=1e-12)
    assert isinstance(r.solvability, Gauge)
