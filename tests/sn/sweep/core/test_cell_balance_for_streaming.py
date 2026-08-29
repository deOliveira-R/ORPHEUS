r"""Foundation tests for :func:`cell_balance_for_streaming`.

Issue #197 PR-TYPED-6 — pins the vectorized cell-balance helper as
Pattern 2's single source of truth.  Both
:meth:`DiamondDifference.residual` (n_mask=1 case) and the unified
matvec (n_mask=N case) consume this helper.

PR-TYPED-6.5 Phase 2.11 — helper signature changed
---------------------------------------------------
The helper now accepts ``angular_denom_term`` (shape ``(n_mask,)``)
and ``angular_numer_upstream`` (shape ``(ng, n_mask)``) — the
operated-on contributions from the closure strategy
(``closure.cell_contribution(...)``).  Pre-Phase-2.11 the helper
accepted ``delta_A_over_w``, ``c_in``, ``c_out``, ``psi_angular_upstream``
(M-M-specific names), with the closure algebra inlined.  These tests
construct the angular contributions explicitly to verify the helper's
own algebra (no closure dependence at the test surface).

Pinned invariants
=================

1. **Hand-written literal pins at n_mask=1**: the helper's output on
   the module's synthetic curvilinear and slab packets equals literals
   computed OFFLINE by hand arithmetic (``array_equal``).  ⚠ P4.9a
   claim-class note: until the ``cell_balance_terms`` retirement these
   two rows compared TWO independently-spelled production helpers; the
   literals are the anchor that survives — computed from the balance
   equation as written, never by re-calling the helper.

2. **Vectorization invariance**: calling the helper at n_mask=N
   produces per-ordinate results bit-identical to calling it n_mask
   times at n_mask=1.

3. **Slab degeneracy**: with zero ``angular_denom_term`` and zero
   ``angular_numer_upstream`` the helper reduces to the slab form
   ``denom = 2|μ|·face_area_downstream + σ_t·V`` / ``numer = |μ|·face_area_total·ψ_in``.

4. **Linearity in psi_face_in**: the helper is affine in
   ``psi_face_in`` and ``angular_numer_upstream`` (linear
   superposition).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.transport.spatial.scheme import StreamingTerms
from orpheus.transport.spatial.cell_balance import (
    cell_balance_for_streaming,
)
from orpheus.transport.spatial.scheme import UpstreamState

# Issue #236 Phase 2 B2 Fix 3 / Step C — the ONE shared hand-transcribed
# surrogate for the M-M ``(c_in, c_out)`` constants (was a private
# byte-identical copy ``_visit_c`` here; unified with ``test_diamond.py``
# and the production-stamp catcher).  Step C: τ / α are now passed
# explicitly (the geometry-side ``StreamingTerms`` τ/α packing was
# retired); these synthetic packets carry their chosen M-M triple as
# module constants below.
from tests.sn.sweep.core._c_surrogate import c_from_constants

# Synthetic M-M (τ, α_in, α_out) triple for the representative
# curvilinear packet — chosen values that were previously baked onto
# the (now-retired) ``StreamingTerms.tau_mm`` / ``alpha_in`` / ``alpha_out``
# fields.  Slab is the neutral element (τ = 1, α = 0).
_CURVILINEAR_TAU = 0.75
_CURVILINEAR_ALPHA_IN = 0.1
_CURVILINEAR_ALPHA_OUT = 0.15
# ΔA/w for the synthetic fixtures (P4.7 — the packet no longer carries
# the fusion; the value lives beside the M-M triple it couples with).
_CURVILINEAR_DAW = 0.3
_SLAB_DAW = 0.0
_SLAB_TAU = 1.0
_SLAB_ALPHA_IN = 0.0
_SLAB_ALPHA_OUT = 0.0

# Per-cell streaming-balance algebra: structural invariants of the
# cell_balance_for_streaming primitive (no theory-page :label:).
# Foundation, not a physics equation gate. (Was a V&V orphan before
# the taxonomy reorg forced a marker.)
pytestmark = pytest.mark.foundation


# ═══════════════════════════════════════════════════════════════════════
# Helper — build a representative curvilinear per-cell packet
# ═══════════════════════════════════════════════════════════════════════


def _curvilinear_streaming_terms() -> StreamingTerms:
    """Build a representative curvilinear (sphere-ish) StreamingTerms.

    The M-M triple lives in the module constants ``_CURVILINEAR_TAU`` /
    ``_CURVILINEAR_ALPHA_{IN,OUT}`` (Step C — the geometry-side ``tau_mm``
    / ``alpha_*`` packing on ``StreamingTerms`` was retired).
    """
    return StreamingTerms(
        face_area_inner=1.2,
        face_area_outer=1.5,
        volume=0.6,
        abs_mu=0.3,
    )


def _slab_streaming_terms() -> StreamingTerms:
    """Slab StreamingTerms (Step 2.5 neutral curvature; neutral M-M triple)."""
    return StreamingTerms(
        face_area_inner=1.0,
        face_area_outer=1.0,
        volume=0.5,
        abs_mu=0.3,
    )


def _scalar_to_vectorized_inputs(
    st: StreamingTerms,
    A_downstream_scalar: float,
    psi_spatial: np.ndarray,
    psi_ang: "np.ndarray | None",
    *,
    tau: float,
    alpha_in: float,
    alpha_out: float,
    delta_A_over_w: float,
):
    """Convert per-cell scalar StreamingTerms to ``(n_mask=1,)`` inputs.

    Returns ``(abs_mu, face_area_downstream, face_area_total, psi_face_in,
    angular_denom_term, angular_numer_upstream)`` — the closure-aware
    PR-TYPED-6.5 Phase 2.11 helper argument shape.  The angular
    contributions are computed inline from the explicit M-M triple
    (Step C — no longer read off ``st``; P4.9a — nor off the visit
    family, which is purely spatial now) to verify
    ``cell_balance_for_streaming`` reproduces the pre-Phase-2.11
    behaviour under equivalent inputs.
    """
    abs_mu = np.array([st.abs_mu], dtype=float)
    face_area_downstream = np.array([A_downstream_scalar], dtype=float)
    face_area_total = np.array([st.face_area_inner + st.face_area_outer], dtype=float)
    dA_w_scalar = delta_A_over_w  # P4.7: supplied beside the M-M triple
    c_out_scalar = alpha_out / tau
    c_in_scalar = (1.0 - tau) / tau * alpha_out + alpha_in
    psi_face_in = psi_spatial[:, None]                          # (ng, 1)
    angular_denom_term = np.array(
        [dA_w_scalar * c_out_scalar], dtype=float,
    )                                                            # (1,)
    if psi_ang is None:
        angular_numer_upstream = np.zeros_like(psi_face_in)
    else:
        angular_numer_upstream = (
            dA_w_scalar * c_in_scalar * psi_ang[:, None]
        )                                                        # (ng, 1)
    return (
        abs_mu, face_area_downstream, face_area_total, psi_face_in,
        angular_denom_term, angular_numer_upstream,
    )


# ═══════════════════════════════════════════════════════════════════════
# Test 1 — bit-identity at n_mask=1
# ═══════════════════════════════════════════════════════════════════════


def test_n_mask_1_hand_written_pin_curvilinear():
    """The helper's curvilinear output equals OFFLINE hand arithmetic —
    ``array_equal`` against pasted literals.

    P4.9a rewire (was ``test_n_mask_1_matches_scalar_curvilinear``): the
    scalar twin ``cell_balance_terms`` retired, so the cross-helper
    comparison lost its subject.  The literals below were computed by
    hand in the balance equation's own operation order —
    ``denom = 2·|μ|·face_area_downstream + (ΔA/w)·c_out + Σ_t·V`` (left fold),
    ``numer = |μ|·face_area_total·ψˢ + (ΔA/w)·c_in·ψᵃ`` — on the module's
    synthetic packet (τ=0.75, α=(0.1, 0.15), face_area_downstream=1.5), [M] validated
    bit-equal against the living helper before it died (2026-08-28).
    This is a claim-class PROMOTION: the values stop being cross-checked
    between two implementations and start being anchored to literals.
    ⛔ Never recompute the reference by calling the helper.
    """
    st = _curvilinear_streaming_terms()
    A_downstream = 1.5
    total_xs = np.array([1.2, 0.8])
    (
        abs_mu, face_area_downstream, face_area_total, psi_face_in,
        angular_denom_term, angular_numer_upstream,
    ) = _scalar_to_vectorized_inputs(
        st, A_downstream, np.array([0.5, 0.7]), np.array([0.4, 0.6]),
        tau=_CURVILINEAR_TAU, alpha_in=_CURVILINEAR_ALPHA_IN,
        alpha_out=_CURVILINEAR_ALPHA_OUT, delta_A_over_w=_CURVILINEAR_DAW,
    )

    denom_v, numer_v = cell_balance_for_streaming(
        abs_mu=abs_mu, A_downstream=face_area_downstream, face_area_total=face_area_total,
        total_xs=total_xs, volume=st.volume,
        psi_face_in=psi_face_in,
        angular_denom_term=angular_denom_term,
        angular_numer_upstream=angular_numer_upstream,
    )

    assert denom_v.shape == (2, 1)
    assert numer_v.shape == (2, 1)
    np.testing.assert_array_equal(
        denom_v[:, 0], np.array([1.6799999999999997, 1.44]),
    )
    np.testing.assert_array_equal(
        numer_v[:, 0], np.array([0.42300000000000004, 0.594]),
    )


def test_n_mask_1_hand_written_pin_slab():
    """The helper's slab output equals OFFLINE hand arithmetic —
    ``array_equal`` against pasted literals.

    P4.9a rewire (was ``test_n_mask_1_matches_scalar_slab``); see the
    curvilinear pin's docstring for the claim-class note.  Slab packet:
    neutral curvature (delta_A_over_w = 0, c = 0, A = 1), so
    ``denom = 2·|μ| + Σ_t·V`` and ``numer = 2|μ|·ψˢ`` — the two
    structural zero-legs on the angular contributions survive as
    independent pins.
    """
    st = _slab_streaming_terms()
    A_downstream = 1.0
    total_xs = np.array([1.2, 0.8])
    (
        abs_mu, face_area_downstream, face_area_total, psi_face_in,
        angular_denom_term, angular_numer_upstream,
    ) = _scalar_to_vectorized_inputs(
        st, A_downstream, np.array([0.5, 0.7]), None,
        tau=_SLAB_TAU, alpha_in=_SLAB_ALPHA_IN, alpha_out=_SLAB_ALPHA_OUT,
        delta_A_over_w=_SLAB_DAW,
    )
    # Slab fixture has delta_A_over_w = 0 and c_in = c_out = 0 → both angular
    # contributions are zero arrays.
    np.testing.assert_array_equal(angular_denom_term, np.zeros(1))
    np.testing.assert_array_equal(angular_numer_upstream, np.zeros_like(psi_face_in))

    denom_v, numer_v = cell_balance_for_streaming(
        abs_mu=abs_mu, A_downstream=face_area_downstream, face_area_total=face_area_total,
        total_xs=total_xs, volume=st.volume,
        psi_face_in=psi_face_in,
        angular_denom_term=angular_denom_term,
        angular_numer_upstream=angular_numer_upstream,
    )

    np.testing.assert_array_equal(denom_v[:, 0], np.array([1.2, 1.0]))
    np.testing.assert_array_equal(numer_v[:, 0], np.array([0.3, 0.42]))


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
    face_area_downstream = rng.uniform(0.5, 2.0, size=n_mask)
    face_area_total = rng.uniform(1.5, 3.5, size=n_mask)
    # Synthetic angular contributions (post-Phase-2.11 inputs).
    angular_denom_term = rng.uniform(0.0, 0.5, size=n_mask)
    angular_numer_upstream = rng.normal(size=(ng, n_mask))
    psi_face_in = rng.normal(size=(ng, n_mask))

    # Vectorized batch call.
    denom_batch, numer_batch = cell_balance_for_streaming(
        abs_mu=abs_mu, A_downstream=face_area_downstream, face_area_total=face_area_total,
        total_xs=total_xs, volume=volume,
        psi_face_in=psi_face_in,
        angular_denom_term=angular_denom_term,
        angular_numer_upstream=angular_numer_upstream,
    )

    # Per-ordinate (n_mask=1) calls.
    for k in range(n_mask):
        denom_k, numer_k = cell_balance_for_streaming(
            abs_mu=abs_mu[k:k + 1],
            A_downstream=face_area_downstream[k:k + 1],
            face_area_total=face_area_total[k:k + 1],
            total_xs=total_xs, volume=volume,
            psi_face_in=psi_face_in[:, k:k + 1],
            angular_denom_term=angular_denom_term[k:k + 1],
            angular_numer_upstream=angular_numer_upstream[:, k:k + 1],
        )
        np.testing.assert_array_equal(denom_batch[:, k], denom_k[:, 0])
        np.testing.assert_array_equal(numer_batch[:, k], numer_k[:, 0])


# ═══════════════════════════════════════════════════════════════════════
# Test 3 — slab degeneracy: helper reduces to (2|μ|face_area_downstream + σ_t V) / (|μ|face_area_total ψ_in)
# ═══════════════════════════════════════════════════════════════════════


def test_slab_degenerate_form():
    """With zero ``angular_denom_term`` and zero ``angular_numer_upstream``
    (the IdentityAngularClosure contribution for slab), the helper
    reduces to the pure-slab form."""
    abs_mu = np.array([0.5, 0.3])
    face_area_downstream = np.array([1.0, 1.0])
    face_area_total = np.array([2.0, 2.0])
    total_xs = np.array([1.5, 0.8])
    volume = 0.4
    psi_face_in = np.array([[0.1, 0.2], [0.3, 0.4]])  # (ng=2, n_mask=2)
    angular_denom_term = np.zeros(2)
    angular_numer_upstream = np.zeros((2, 2))

    denom, numer = cell_balance_for_streaming(
        abs_mu=abs_mu, A_downstream=face_area_downstream, face_area_total=face_area_total,
        total_xs=total_xs, volume=volume,
        psi_face_in=psi_face_in,
        angular_denom_term=angular_denom_term,
        angular_numer_upstream=angular_numer_upstream,
    )

    # Slab form: denom = 2|μ|face_area_downstream + σ_t V (no redistribution).
    expected_denom = (
        2.0 * abs_mu[None, :] * face_area_downstream[None, :]
        + (total_xs * volume)[:, None]
    )
    # Slab form: numer = |μ|face_area_total ψ_in (no angular contribution).
    expected_numer = abs_mu[None, :] * face_area_total[None, :] * psi_face_in

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
    face_area_downstream = rng.uniform(0.5, 2.0, size=n_mask)
    face_area_total = rng.uniform(1.5, 3.5, size=n_mask)
    angular_denom_term = rng.uniform(0.0, 0.5, size=n_mask)
    angular_numer_upstream = rng.normal(size=(ng, n_mask))
    psi_a = rng.normal(size=(ng, n_mask))
    psi_b = rng.normal(size=(ng, n_mask))

    alpha, beta = 0.6, -0.4

    denom_a, numer_a = cell_balance_for_streaming(
        abs_mu=abs_mu, A_downstream=face_area_downstream, face_area_total=face_area_total,
        total_xs=total_xs, volume=volume,
        psi_face_in=psi_a,
        angular_denom_term=angular_denom_term,
        angular_numer_upstream=angular_numer_upstream,
    )
    denom_b, numer_b = cell_balance_for_streaming(
        abs_mu=abs_mu, A_downstream=face_area_downstream, face_area_total=face_area_total,
        total_xs=total_xs, volume=volume,
        psi_face_in=psi_b,
        angular_denom_term=angular_denom_term,
        angular_numer_upstream=angular_numer_upstream,
    )
    denom_c, numer_c = cell_balance_for_streaming(
        abs_mu=abs_mu, A_downstream=face_area_downstream, face_area_total=face_area_total,
        total_xs=total_xs, volume=volume,
        psi_face_in=alpha * psi_a + beta * psi_b,
        angular_denom_term=angular_denom_term,
        angular_numer_upstream=angular_numer_upstream,
    )

    # denom is independent of psi_face_in.
    np.testing.assert_allclose(denom_a, denom_b, rtol=0, atol=0)
    np.testing.assert_allclose(denom_a, denom_c, rtol=0, atol=0)
    # numer_upstream is affine — linear combo of psi_face_in lifts
    # to linear combo of numer outputs (since the angular contribution
    # is held fixed).
    spatial_a = abs_mu[None, :] * face_area_total[None, :] * psi_a
    spatial_b = abs_mu[None, :] * face_area_total[None, :] * psi_b
    spatial_c = abs_mu[None, :] * face_area_total[None, :] * (alpha * psi_a + beta * psi_b)
    np.testing.assert_allclose(
        numer_c, spatial_c + angular_numer_upstream, rtol=1e-14,
    )
    np.testing.assert_allclose(
        alpha * spatial_a + beta * spatial_b, spatial_c, rtol=1e-14,
    )


def test_angular_upstream_none_equals_zero_angular_term():
    """``psi_angular_upstream=None`` must produce the same denom and
    numer_upstream as supplying a zero-valued (ng, n_mask) array,
    modulo dtype (post-Phase-2.11 caller supplies the zeros)."""
    rng = np.random.default_rng(42)
    ng, n_mask = 2, 3
    total_xs = rng.uniform(0.5, 1.5, size=ng)
    volume = 0.3
    abs_mu = rng.uniform(0.1, 1.0, size=n_mask)
    face_area_downstream = rng.uniform(0.5, 2.0, size=n_mask)
    face_area_total = rng.uniform(1.5, 3.5, size=n_mask)
    angular_denom_term = rng.uniform(0.0, 0.5, size=n_mask)
    psi_face_in = rng.normal(size=(ng, n_mask))
    zero_ang = np.zeros((ng, n_mask))

    # With the post-Phase-2.11 interface, the caller supplies the
    # angular contribution explicitly.  IdentityAngularClosure returns
    # zeros; an unwitting slab caller would do the same by passing
    # ``np.zeros((ng, n_mask))`` — the helper's behaviour must be
    # invariant under either form.  This test pins the no-effect
    # property of a zero angular contribution.
    denom_zero, numer_zero = cell_balance_for_streaming(
        abs_mu=abs_mu, A_downstream=face_area_downstream, face_area_total=face_area_total,
        total_xs=total_xs, volume=volume,
        psi_face_in=psi_face_in,
        angular_denom_term=angular_denom_term,
        angular_numer_upstream=zero_ang,
    )
    # Build the expected outputs directly.
    expected_denom = (
        2.0 * abs_mu[None, :] * face_area_downstream[None, :]
        + angular_denom_term[None, :]
        + (total_xs * volume)[:, None]
    )
    expected_numer = abs_mu[None, :] * face_area_total[None, :] * psi_face_in
    np.testing.assert_array_equal(denom_zero, expected_denom)
    np.testing.assert_array_equal(numer_zero, expected_numer)


# ═══════════════════════════════════════════════════════════════════════
# Test 5 — residual delegation (DiamondDifference.residual now uses
# cell_balance_for_streaming).  Pin the n_mask=1 round-trip identity.
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize(
    "st_factory, face_area_downstream, ang_upstream, triple, daw",
    [
        (
            _curvilinear_streaming_terms, 1.5, np.array([0.4, 0.6]),
            (_CURVILINEAR_TAU, _CURVILINEAR_ALPHA_IN, _CURVILINEAR_ALPHA_OUT),
            _CURVILINEAR_DAW,
        ),
        (
            _slab_streaming_terms, 1.0, None,
            (_SLAB_TAU, _SLAB_ALPHA_IN, _SLAB_ALPHA_OUT),
            _SLAB_DAW,
        ),
    ],
    ids=["curvilinear", "slab"],
)
def test_diamond_residual_consumes_cell_balance_for_streaming(
    st_factory, face_area_downstream, ang_upstream, triple, daw,
):
    """Verify :meth:`DiamondDifference.residual` produces the same
    output as direct delegation to :func:`cell_balance_for_streaming`."""
    from orpheus.transport.spatial.scheme import CellVisit
    from orpheus.transport.spatial.diamond import DiamondDifference

    tau, alpha_in, alpha_out = triple
    st = st_factory()
    c_in, c_out = c_from_constants(tau, alpha_in, alpha_out)
    visit = CellVisit(
        cell_idx=7, streaming_terms=st, face_area_downstream=face_area_downstream,
    )
    total_xs = np.array([1.2, 0.8])
    source = np.array([0.05, 0.07])
    cell_avg = np.array([0.9, 1.1])
    psi_spatial = np.array([0.5, 0.7])
    upstream = UpstreamState(spatial_upstream=psi_spatial)

    dd = DiamondDifference()
    # P4.9a: the caller assembles the closure contributions.
    residual = dd.residual(
        cell_avg, visit, total_xs, source, upstream,
        angular_denom_term=daw * c_out,
        angular_numer_upstream=(
            None if ang_upstream is None
            else daw * c_in * ang_upstream
        ),
    )

    # Direct delegation to the vectorized helper.
    (
        abs_mu, A_d, A_tot, psi_face_in,
        angular_denom_term, angular_numer_upstream,
    ) = _scalar_to_vectorized_inputs(
        st, face_area_downstream, psi_spatial, ang_upstream,
        tau=tau, alpha_in=alpha_in, alpha_out=alpha_out, delta_A_over_w=daw,
    )
    denom, numer = cell_balance_for_streaming(
        abs_mu=abs_mu, A_downstream=A_d, face_area_total=A_tot,
        total_xs=total_xs, volume=st.volume,
        psi_face_in=psi_face_in,
        angular_denom_term=angular_denom_term,
        angular_numer_upstream=angular_numer_upstream,
    )
    expected = denom[:, 0] * cell_avg - (source + numer[:, 0])
    np.testing.assert_array_equal(residual, expected)


def test_diamond_residual_round_trip_at_converged_cell_avg():
    """At the converged cell_avg (from update), the residual is zero
    to FP rounding — the canonical DiscretizationScheme apply-solve consistency
    contract (Issue #196 Phase G Step 1)."""
    from orpheus.transport.spatial.scheme import CellVisit
    from orpheus.transport.spatial.diamond import DiamondDifference

    st = _curvilinear_streaming_terms()
    c_in, c_out = c_from_constants(
        _CURVILINEAR_TAU, _CURVILINEAR_ALPHA_IN, _CURVILINEAR_ALPHA_OUT,
    )
    visit = CellVisit(
        cell_idx=3, streaming_terms=st, face_area_downstream=1.5,
    )
    total_xs = np.array([1.2, 0.8])
    source = np.array([0.05, 0.07])
    upstream = UpstreamState(spatial_upstream=np.array([0.5, 0.7]))
    adt = _CURVILINEAR_DAW * c_out
    anu = _CURVILINEAR_DAW * c_in * np.array([0.4, 0.6])

    dd = DiamondDifference()
    result = dd.update(
        visit, total_xs, source, upstream,
        angular_denom_term=adt, angular_numer_upstream=anu,
    )
    residual = dd.residual(
        result.cell_average_flux, visit, total_xs, source, upstream,
        angular_denom_term=adt, angular_numer_upstream=anu,
    )
    # The solve form produces the exact solution; residual is zero.
    np.testing.assert_allclose(residual, 0.0, atol=1e-13)
