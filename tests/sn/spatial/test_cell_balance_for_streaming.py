r"""Foundation tests for :func:`cell_balance_for_streaming`.

Issue #197 PR-TYPED-6 — pins the vectorized cell-balance helper as
Pattern 2's single source of truth.  Both
:meth:`DiamondDifference.residual` (n_mask=1 case) and (later)
:meth:`StreamingOperator.apply` (n_mask=N case) consume this helper.

Pinned invariants
=================

1. **Bit-identity at n_mask=1**: the vectorized form matches
   :func:`cell_balance_terms` exactly on the same per-cell scalar
   inputs, for both curvilinear and slab geometry.

2. **Vectorization invariance**: calling the helper at n_mask=N
   produces per-ordinate results bit-identical to calling it n_mask
   times at n_mask=1.

3. **Slab degeneracy**: with neutral curvature (dA_w=0, c_in=c_out=0)
   and ``psi_angular_upstream=None`` the helper reduces to the slab
   form ``denom = 2|μ|·A_down + σ_t·V`` / ``numer = |μ|·A_total·ψ_in``.

4. **Linearity in psi_face_in**: the helper is affine in
   ``psi_face_in`` and ``psi_angular_upstream`` (linear superposition).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.reduced_operator import StreamingTerms
from orpheus.sn.spatial.cell_balance import (
    cell_balance_for_streaming,
    cell_balance_terms,
)
from orpheus.sn.spatial.cell_update import UpstreamState


# ═══════════════════════════════════════════════════════════════════════
# Helper — build a representative curvilinear per-cell packet
# ═══════════════════════════════════════════════════════════════════════


def _curvilinear_streaming_terms() -> StreamingTerms:
    """Build a representative curvilinear (sphere-ish) StreamingTerms."""
    return StreamingTerms(
        chord_length=0.5,
        mu=0.3,
        face_area_inner=1.2,
        face_area_outer=1.5,
        delta_A_over_w=0.3,
        alpha_in=0.1,
        alpha_out=0.15,
        tau_mm=0.75,
        volume=0.6,
        abs_mu=0.3,
    )


def _slab_streaming_terms() -> StreamingTerms:
    """Slab StreamingTerms (Step 2.5 neutral curvature)."""
    return StreamingTerms(
        chord_length=0.5,
        mu=0.3,
        face_area_inner=1.0,
        face_area_outer=1.0,
        delta_A_over_w=0.0,
        alpha_in=0.0,
        alpha_out=0.0,
        tau_mm=1.0,
        volume=0.5,
        abs_mu=0.3,
    )


def _scalar_to_vectorized_inputs(
    st: StreamingTerms,
    A_downstream_scalar: float,
    upstream: UpstreamState,
):
    """Convert per-cell scalar StreamingTerms to (n_mask=1,) inputs."""
    abs_mu = np.array([st.abs_mu], dtype=float)
    A_down = np.array([A_downstream_scalar], dtype=float)
    A_total = np.array([st.face_area_inner + st.face_area_outer], dtype=float)
    dA_w = np.array([st.delta_A_over_w], dtype=float)
    tau = st.tau_mm
    c_out = np.array([st.alpha_out / tau], dtype=float)
    c_in = np.array(
        [(1.0 - tau) / tau * st.alpha_out + st.alpha_in], dtype=float,
    )
    psi_face_in = upstream.spatial_upstream[:, None]
    if upstream.angular_upstream is None:
        psi_ang = None
    else:
        psi_ang = upstream.angular_upstream[:, None]
    return abs_mu, A_down, A_total, dA_w, c_in, c_out, psi_face_in, psi_ang


# ═══════════════════════════════════════════════════════════════════════
# Test 1 — bit-identity at n_mask=1
# ═══════════════════════════════════════════════════════════════════════


def test_n_mask_1_matches_scalar_curvilinear():
    """Pattern 2: the vectorized form at n_mask=1 reproduces the
    scalar :func:`cell_balance_terms` output bit-for-bit on
    curvilinear inputs."""
    st = _curvilinear_streaming_terms()
    A_downstream = 1.5
    total_xs = np.array([1.2, 0.8])
    upstream = UpstreamState(
        spatial_upstream=np.array([0.5, 0.7]),
        angular_upstream=np.array([0.4, 0.6]),
    )

    terms = cell_balance_terms(st, A_downstream, total_xs, upstream)

    inputs = _scalar_to_vectorized_inputs(st, A_downstream, upstream)
    abs_mu, A_down, A_total, dA_w, c_in, c_out, psi_face_in, psi_ang = inputs

    denom_v, numer_v = cell_balance_for_streaming(
        abs_mu=abs_mu, A_downstream=A_down, A_total=A_total,
        dA_w=dA_w, c_in=c_in, c_out=c_out,
        total_xs=total_xs, volume=st.volume,
        psi_face_in=psi_face_in, psi_angular_upstream=psi_ang,
    )

    # Bit-identity: both helpers compute the same algebra with the
    # same operation order (named intermediates).  Differences would
    # only appear from FP-non-associativity in the multi-term sum;
    # both forms add the three contributions in the same order.
    assert denom_v.shape == (2, 1)
    assert numer_v.shape == (2, 1)
    np.testing.assert_array_equal(denom_v[:, 0], terms.denom)
    np.testing.assert_array_equal(numer_v[:, 0], terms.numer_upstream)


def test_n_mask_1_matches_scalar_slab():
    """Pattern 2: bit-identity at n_mask=1 for slab geometry
    (psi_angular_upstream is None branch)."""
    st = _slab_streaming_terms()
    A_downstream = 1.0
    total_xs = np.array([1.2, 0.8])
    upstream = UpstreamState(
        spatial_upstream=np.array([0.5, 0.7]),
        angular_upstream=None,
    )

    terms = cell_balance_terms(st, A_downstream, total_xs, upstream)

    inputs = _scalar_to_vectorized_inputs(st, A_downstream, upstream)
    abs_mu, A_down, A_total, dA_w, c_in, c_out, psi_face_in, psi_ang = inputs
    assert psi_ang is None  # slab carries no angular state

    denom_v, numer_v = cell_balance_for_streaming(
        abs_mu=abs_mu, A_downstream=A_down, A_total=A_total,
        dA_w=dA_w, c_in=c_in, c_out=c_out,
        total_xs=total_xs, volume=st.volume,
        psi_face_in=psi_face_in, psi_angular_upstream=None,
    )

    np.testing.assert_array_equal(denom_v[:, 0], terms.denom)
    np.testing.assert_array_equal(numer_v[:, 0], terms.numer_upstream)


# ═══════════════════════════════════════════════════════════════════════
# Test 2 — vectorization invariance (n_mask=1 chunks vs n_mask=N batch)
# ═══════════════════════════════════════════════════════════════════════


def test_vectorization_invariance_curvilinear():
    """Vectorizing over n_mask gives the same per-ordinate result as
    calling the helper n_mask times at n_mask=1."""
    rng = np.random.default_rng(0xC4FE_BABE)
    ng = 3
    n_mask = 5
    total_xs = rng.uniform(0.5, 1.5, size=ng)
    volume = 0.42

    abs_mu = rng.uniform(0.1, 1.0, size=n_mask)
    A_down = rng.uniform(0.5, 2.0, size=n_mask)
    A_total = rng.uniform(1.5, 3.5, size=n_mask)
    dA_w = rng.uniform(0.01, 0.5, size=n_mask)
    c_in = rng.uniform(0.0, 0.5, size=n_mask)
    c_out = rng.uniform(0.0, 0.5, size=n_mask)
    psi_face_in = rng.normal(size=(ng, n_mask))
    psi_ang = rng.normal(size=(ng, n_mask))

    # Vectorized batch call.
    denom_batch, numer_batch = cell_balance_for_streaming(
        abs_mu=abs_mu, A_downstream=A_down, A_total=A_total,
        dA_w=dA_w, c_in=c_in, c_out=c_out,
        total_xs=total_xs, volume=volume,
        psi_face_in=psi_face_in, psi_angular_upstream=psi_ang,
    )

    # Per-ordinate (n_mask=1) calls.
    for k in range(n_mask):
        denom_k, numer_k = cell_balance_for_streaming(
            abs_mu=abs_mu[k:k + 1],
            A_downstream=A_down[k:k + 1],
            A_total=A_total[k:k + 1],
            dA_w=dA_w[k:k + 1],
            c_in=c_in[k:k + 1],
            c_out=c_out[k:k + 1],
            total_xs=total_xs, volume=volume,
            psi_face_in=psi_face_in[:, k:k + 1],
            psi_angular_upstream=psi_ang[:, k:k + 1],
        )
        np.testing.assert_array_equal(denom_batch[:, k], denom_k[:, 0])
        np.testing.assert_array_equal(numer_batch[:, k], numer_k[:, 0])


# ═══════════════════════════════════════════════════════════════════════
# Test 3 — slab degeneracy: helper reduces to (2|μ|A_down + σ_t V) / (|μ|A_total ψ_in)
# ═══════════════════════════════════════════════════════════════════════


def test_slab_degenerate_form():
    """With neutral curvature (dA_w=0, c_in=c_out=0) and no angular
    upstream, the helper reduces to the pure-slab form."""
    abs_mu = np.array([0.5, 0.3])
    A_down = np.array([1.0, 1.0])
    A_total = np.array([2.0, 2.0])
    dA_w = np.zeros(2)
    c_in = np.zeros(2)
    c_out = np.zeros(2)
    total_xs = np.array([1.5, 0.8])
    volume = 0.4
    psi_face_in = np.array([[0.1, 0.2], [0.3, 0.4]])  # (ng=2, n_mask=2)

    denom, numer = cell_balance_for_streaming(
        abs_mu=abs_mu, A_downstream=A_down, A_total=A_total,
        dA_w=dA_w, c_in=c_in, c_out=c_out,
        total_xs=total_xs, volume=volume,
        psi_face_in=psi_face_in, psi_angular_upstream=None,
    )

    # Slab form: denom = 2|μ|A_down + σ_t V (no redistribution).
    expected_denom = (
        2.0 * abs_mu[None, :] * A_down[None, :]
        + (total_xs * volume)[:, None]
    )
    # Slab form: numer = |μ|A_total ψ_in (no angular contribution).
    expected_numer = abs_mu[None, :] * A_total[None, :] * psi_face_in

    np.testing.assert_array_equal(denom, expected_denom)
    np.testing.assert_array_equal(numer, expected_numer)


# ═══════════════════════════════════════════════════════════════════════
# Test 4 — linearity in psi_face_in and psi_angular_upstream
# ═══════════════════════════════════════════════════════════════════════


def test_linearity_in_psi_face_in():
    """The helper is affine in ``psi_face_in``: the numer_upstream's
    spatial-streaming contribution is linear, the denom is independent."""
    rng = np.random.default_rng(0xDECAFBAD)
    ng, n_mask = 2, 4
    total_xs = rng.uniform(0.5, 1.5, size=ng)
    volume = 0.3
    abs_mu = rng.uniform(0.1, 1.0, size=n_mask)
    A_down = rng.uniform(0.5, 2.0, size=n_mask)
    A_total = rng.uniform(1.5, 3.5, size=n_mask)
    dA_w = rng.uniform(0.01, 0.5, size=n_mask)
    c_in = rng.uniform(0.0, 0.5, size=n_mask)
    c_out = rng.uniform(0.0, 0.5, size=n_mask)
    psi_a = rng.normal(size=(ng, n_mask))
    psi_b = rng.normal(size=(ng, n_mask))
    psi_ang = rng.normal(size=(ng, n_mask))

    alpha, beta = 0.6, -0.4

    denom_a, numer_a = cell_balance_for_streaming(
        abs_mu=abs_mu, A_downstream=A_down, A_total=A_total,
        dA_w=dA_w, c_in=c_in, c_out=c_out,
        total_xs=total_xs, volume=volume,
        psi_face_in=psi_a, psi_angular_upstream=psi_ang,
    )
    denom_b, numer_b = cell_balance_for_streaming(
        abs_mu=abs_mu, A_downstream=A_down, A_total=A_total,
        dA_w=dA_w, c_in=c_in, c_out=c_out,
        total_xs=total_xs, volume=volume,
        psi_face_in=psi_b, psi_angular_upstream=psi_ang,
    )
    denom_c, numer_c = cell_balance_for_streaming(
        abs_mu=abs_mu, A_downstream=A_down, A_total=A_total,
        dA_w=dA_w, c_in=c_in, c_out=c_out,
        total_xs=total_xs, volume=volume,
        psi_face_in=alpha * psi_a + beta * psi_b,
        psi_angular_upstream=psi_ang,
    )

    # denom is independent of psi_face_in.
    np.testing.assert_allclose(denom_a, denom_b, rtol=0, atol=0)
    np.testing.assert_allclose(denom_a, denom_c, rtol=0, atol=0)
    # numer_upstream is affine — linear combo of psi_face_in lifts
    # to linear combo of numer outputs (since psi_ang is the same).
    # Both branches share the same angular contribution; the spatial
    # contribution scales linearly with psi_face_in.
    diff_angular = numer_a - abs_mu[None, :] * A_total[None, :] * psi_a
    # diff_angular is the (constant) angular contribution.
    spatial_a = abs_mu[None, :] * A_total[None, :] * psi_a
    spatial_b = abs_mu[None, :] * A_total[None, :] * psi_b
    spatial_c = abs_mu[None, :] * A_total[None, :] * (alpha * psi_a + beta * psi_b)
    np.testing.assert_allclose(
        numer_c, spatial_c + diff_angular, rtol=1e-14,
    )
    np.testing.assert_allclose(
        alpha * spatial_a + beta * spatial_b, spatial_c, rtol=1e-14,
    )


def test_angular_upstream_none_equals_zero_angular_term():
    """``psi_angular_upstream=None`` must produce the same denom and
    numer_upstream as supplying a zero-valued (ng, n_mask) array,
    modulo dtype (None branch uses zeros_like)."""
    rng = np.random.default_rng(42)
    ng, n_mask = 2, 3
    total_xs = rng.uniform(0.5, 1.5, size=ng)
    volume = 0.3
    abs_mu = rng.uniform(0.1, 1.0, size=n_mask)
    A_down = rng.uniform(0.5, 2.0, size=n_mask)
    A_total = rng.uniform(1.5, 3.5, size=n_mask)
    dA_w = rng.uniform(0.01, 0.5, size=n_mask)
    c_in = rng.uniform(0.0, 0.5, size=n_mask)
    c_out = rng.uniform(0.0, 0.5, size=n_mask)
    psi_face_in = rng.normal(size=(ng, n_mask))
    zero_ang = np.zeros((ng, n_mask))

    denom_none, numer_none = cell_balance_for_streaming(
        abs_mu=abs_mu, A_downstream=A_down, A_total=A_total,
        dA_w=dA_w, c_in=c_in, c_out=c_out,
        total_xs=total_xs, volume=volume,
        psi_face_in=psi_face_in, psi_angular_upstream=None,
    )
    denom_zero, numer_zero = cell_balance_for_streaming(
        abs_mu=abs_mu, A_downstream=A_down, A_total=A_total,
        dA_w=dA_w, c_in=c_in, c_out=c_out,
        total_xs=total_xs, volume=volume,
        psi_face_in=psi_face_in, psi_angular_upstream=zero_ang,
    )
    np.testing.assert_array_equal(denom_none, denom_zero)
    np.testing.assert_array_equal(numer_none, numer_zero)


# ═══════════════════════════════════════════════════════════════════════
# Test 5 — residual delegation (DiamondDifference.residual now uses
# cell_balance_for_streaming).  Pin the n_mask=1 round-trip identity.
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize(
    "st_factory, A_down, ang_upstream",
    [
        (_curvilinear_streaming_terms, 1.5, np.array([0.4, 0.6])),
        (_slab_streaming_terms, 1.0, None),
    ],
    ids=["curvilinear", "slab"],
)
def test_diamond_residual_consumes_cell_balance_for_streaming(
    st_factory, A_down, ang_upstream,
):
    """Verify :meth:`DiamondDifference.residual` produces the same
    output as direct delegation to :func:`cell_balance_for_streaming`."""
    from orpheus.sn.spatial.cell_update import CellVisit
    from orpheus.sn.spatial.diamond import DiamondDifference

    st = st_factory()
    visit = CellVisit(
        cell_idx=7, streaming_terms=st, face_area_downstream=A_down,
    )
    total_xs = np.array([1.2, 0.8])
    source = np.array([0.05, 0.07])
    cell_avg = np.array([0.9, 1.1])
    upstream = UpstreamState(
        spatial_upstream=np.array([0.5, 0.7]),
        angular_upstream=ang_upstream,
    )

    dd = DiamondDifference()
    residual = dd.residual(cell_avg, visit, total_xs, source, upstream)

    # Direct delegation to the vectorized helper.
    inputs = _scalar_to_vectorized_inputs(st, A_down, upstream)
    abs_mu, A_d, A_tot, dA_w, c_in, c_out, psi_face_in, psi_ang = inputs
    denom, numer = cell_balance_for_streaming(
        abs_mu=abs_mu, A_downstream=A_d, A_total=A_tot,
        dA_w=dA_w, c_in=c_in, c_out=c_out,
        total_xs=total_xs, volume=st.volume,
        psi_face_in=psi_face_in, psi_angular_upstream=psi_ang,
    )
    expected = denom[:, 0] * cell_avg - (source + numer[:, 0])
    np.testing.assert_array_equal(residual, expected)


def test_diamond_residual_round_trip_at_converged_cell_avg():
    """At the converged cell_avg (from update), the residual is zero
    to FP rounding — the canonical CellUpdate apply-solve consistency
    contract (Issue #196 Phase G Step 1)."""
    from orpheus.sn.spatial.cell_update import CellVisit
    from orpheus.sn.spatial.diamond import DiamondDifference

    st = _curvilinear_streaming_terms()
    visit = CellVisit(
        cell_idx=3, streaming_terms=st, face_area_downstream=1.5,
    )
    total_xs = np.array([1.2, 0.8])
    source = np.array([0.05, 0.07])
    upstream = UpstreamState(
        spatial_upstream=np.array([0.5, 0.7]),
        angular_upstream=np.array([0.4, 0.6]),
    )

    dd = DiamondDifference()
    result = dd.update(visit, total_xs, source, upstream)
    residual = dd.residual(
        result.cell_average_flux, visit, total_xs, source, upstream,
    )
    # The solve form produces the exact solution; residual is zero.
    np.testing.assert_allclose(residual, 0.0, atol=1e-13)
