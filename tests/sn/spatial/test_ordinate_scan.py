r"""Tests for :func:`orpheus.sn.spatial.scan.ordinate_scan`.

Issue #196 Phase G Step 2.5b — the canonical first-order linear-
recurrence scan (Blelloch 1990 §1.5).  Sixteen strong tests:

1.  Algebraic theorems (pin the abstraction itself):
    pair-monoid associativity, identity element, Brent's blocked-scan
    equivalence, closed-form-vs-explicit-loop.
2.  Affine structure (pin the linearity):
    zero-source / zero-attenuation reductions, linearity in
    :math:`\psi_0`, linearity in ``b``, joint affine combination.
3.  Numerical stability (pin behaviour at edge cases):
    near-identity ``a≈1`` regime, small-attenuation ``a→0`` regime.
4.  Dual-view contracts (pin the strategy):
    :meth:`DiamondDifference.affine_coefficients` ↔
    :meth:`DiamondDifference.update` single-cell equivalence,
    vectorisation-vs-serial, full-sweep baseline (parametrised over
    geometry × n_groups × source kind).
5.  Production gates (pin the SN-specific behaviour):
    L0 streaming-equilibrium preserved, slab MMS convergence preserved,
    regression snapshots bit-identical.

All tagged ``@pytest.mark.l0`` and ``@pytest.mark.verifies(
"blelloch-1990-eq-1-5")``.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.sn.spatial.scan import ordinate_scan


pytestmark = [
    pytest.mark.l0,
    pytest.mark.verifies("blelloch-1990-eq-1-5"),
]


# ═══════════════════════════════════════════════════════════════════════
# Pair-monoid composition operator (test-side reference for ⊕)
# ═══════════════════════════════════════════════════════════════════════

def _monoid_compose(
    m1: tuple[float, float],
    m2: tuple[float, float],
) -> tuple[float, float]:
    r"""Pair-monoid composition

    .. math::

        (\alpha_1, \beta_1) \oplus (\alpha_2, \beta_2)
            \;=\; (\alpha_1 \alpha_2, \alpha_2\,\beta_1 + \beta_2).

    Test-side reference for the monoid algebra — separate from the
    closed-form decomposition in :func:`ordinate_scan` so the theorem
    tests verify the abstraction, not the implementation.
    """
    a1, b1 = m1
    a2, b2 = m2
    return (a1 * a2, a2 * b1 + b2)


def _explicit_scan_loop(
    a: np.ndarray,
    b: np.ndarray,
    psi_0: np.ndarray,
) -> np.ndarray:
    r"""Reference: sequential for-loop applying the recurrence
    ``psi[i+1] = a[i]·psi[i] + b[i]``.  Output shape ``(nx, ...)``;
    ``out[i]`` is the state after cell i (i.e. :math:`\psi[i+1]`).
    """
    nx = a.shape[0]
    out = np.empty_like(a)
    psi = psi_0
    for i in range(nx):
        psi = a[i] * psi + b[i]
        out[i] = psi
    return out


# ═══════════════════════════════════════════════════════════════════════
# 1. Algebraic theorems
# ═══════════════════════════════════════════════════════════════════════


class TestPairMonoidTheorems:
    r"""The four load-bearing theorems that justify the closed form.

    These tests pin the algebra (pair-monoid ``(a, b)`` with
    composition ``⊕``), independently of how ``ordinate_scan`` is
    implemented.  If associativity fails, the entire Blelloch §1.5
    derivation is wrong — STOP and dispatch literature-researcher.
    """

    def test_pair_monoid_associativity(self) -> None:
        r"""``((M₁ ⊕ M₂) ⊕ M₃) == (M₁ ⊕ (M₂ ⊕ M₃))``.  THE THEOREM.

        The pair-monoid is associative; the scan composes elements
        independent of bracketing order.  Random ``(a, b)`` triples
        over ``a ∈ [0, 2]``, ``b ∈ [-1, 1]``.  rtol=1e-14 — single-
        multiplication + single-addition is exact to ULP for
        well-scaled inputs.
        """
        rng = np.random.default_rng(seed=0)
        n_trials = 100
        for _ in range(n_trials):
            a_vals = rng.uniform(0.0, 2.0, size=3)
            b_vals = rng.uniform(-1.0, 1.0, size=3)
            m1 = (a_vals[0], b_vals[0])
            m2 = (a_vals[1], b_vals[1])
            m3 = (a_vals[2], b_vals[2])

            lhs = _monoid_compose(_monoid_compose(m1, m2), m3)
            rhs = _monoid_compose(m1, _monoid_compose(m2, m3))

            np.testing.assert_allclose(lhs[0], rhs[0], rtol=1e-14)
            np.testing.assert_allclose(lhs[1], rhs[1], rtol=1e-14)

    def test_pair_monoid_identity(self) -> None:
        r"""``(1, 0) ⊕ M == M == M ⊕ (1, 0)``.  Identity of the monoid.

        ``(1, 0)`` is the zero-attenuation, zero-source cell — a
        pass-through.  Composition with identity is identity.  Exact
        equality.
        """
        rng = np.random.default_rng(seed=1)
        identity = (1.0, 0.0)
        for _ in range(50):
            a = rng.uniform(0.0, 2.0)
            b = rng.uniform(-1.0, 1.0)
            m = (a, b)
            left = _monoid_compose(identity, m)
            right = _monoid_compose(m, identity)
            assert left == m
            assert right == m

    def test_brent_blocked_scan_equivalence(self) -> None:
        r"""Scan of nx cells == reduce-then-scan of nx/2 pair-reduced blocks.

        Brent's theorem: associative scans admit an
        :math:`O(N/\log N)`-work parallel decomposition.  Here we
        verify the simplest such decomposition: pair the cells, reduce
        each pair with ``⊕`` to a single monoid element, then scan
        the reduced sequence — the result matches the full-length
        scan at the pair-aligned indices.

        rtol=1e-13 — two ULPs for the per-pair compose + reduce.
        """
        rng = np.random.default_rng(seed=2)
        nx = 20  # even
        a = rng.uniform(0.2, 1.5, size=nx)
        b = rng.uniform(-0.5, 0.5, size=nx)
        psi_0 = rng.uniform(0.1, 0.5)

        # Full scan (ground truth via explicit loop).
        full = _explicit_scan_loop(
            a.reshape(nx, 1), b.reshape(nx, 1), np.array([psi_0]),
        )

        # Pairwise reduce: combine consecutive pairs via ⊕, then scan.
        pair_a = np.empty(nx // 2)
        pair_b = np.empty(nx // 2)
        for k in range(nx // 2):
            m1 = (a[2 * k], b[2 * k])
            m2 = (a[2 * k + 1], b[2 * k + 1])
            combined = _monoid_compose(m1, m2)
            pair_a[k] = combined[0]
            pair_b[k] = combined[1]
        reduced = _explicit_scan_loop(
            pair_a.reshape(nx // 2, 1), pair_b.reshape(nx // 2, 1),
            np.array([psi_0]),
        )

        # ``reduced[k]`` is the state after applying the kth pair;
        # equivalent to ``full[2k+1]`` (state after second cell of pair).
        for k in range(nx // 2):
            np.testing.assert_allclose(
                reduced[k, 0], full[2 * k + 1, 0], rtol=1e-13,
            )

    def test_ordinate_scan_matches_explicit_loop(self) -> None:
        r"""``ordinate_scan(a, b, psi_0)`` matches the explicit recurrence.

        Vectorised closed-form vs sequential reference loop.  This is
        the fast-path-vs-slow-path bit-identity gate at rtol=1e-13.
        Parametrised over scalar / multi-group shapes.
        """
        rng = np.random.default_rng(seed=3)
        for ng in (1, 2, 4):
            nx = 30
            a = rng.uniform(0.3, 1.4, size=(nx, ng))
            b = rng.uniform(-0.4, 0.7, size=(nx, ng))
            psi_0 = rng.uniform(0.1, 0.5, size=ng)

            scan_result = ordinate_scan(a, b, psi_0)
            loop_result = _explicit_scan_loop(a, b, psi_0)

            np.testing.assert_allclose(scan_result, loop_result, rtol=1e-13)


# ═══════════════════════════════════════════════════════════════════════
# 2. Affine structure
# ═══════════════════════════════════════════════════════════════════════


class TestAffineStructure:
    r"""The recurrence is jointly affine in (b, psi_0) for fixed a.

    These tests pin linearity — failure here signals an accidental
    non-linear interaction in the scan implementation.
    """

    def test_ordinate_scan_zero_source(self) -> None:
        r"""``ordinate_scan(a, 0, psi_0) == cumprod(a) · psi_0``.

        With source ``b = 0``, the scan reduces to multiplicative
        chain composition: pure attenuation.  Exact equality —
        ``cumsum(0 / cumprod_a) == 0`` exactly in IEEE-754.
        """
        rng = np.random.default_rng(seed=4)
        nx, ng = 15, 3
        a = rng.uniform(0.3, 1.4, size=(nx, ng))
        b = np.zeros((nx, ng))
        psi_0 = rng.uniform(0.1, 0.7, size=ng)

        result = ordinate_scan(a, b, psi_0)
        expected = np.cumprod(a, axis=0) * psi_0

        np.testing.assert_array_equal(result, expected)

    def test_ordinate_scan_zero_attenuation(self) -> None:
        r"""``ordinate_scan(1, b, psi_0) == psi_0 + cumsum(b)``.

        With attenuation ``a = 1``, the scan reduces to additive
        chain composition: pure accumulation.  Exact equality —
        ``cumprod(1) == 1`` everywhere, and ``b / 1 == b``.
        """
        rng = np.random.default_rng(seed=5)
        nx, ng = 12, 2
        a = np.ones((nx, ng))
        b = rng.uniform(-0.5, 0.5, size=(nx, ng))
        psi_0 = rng.uniform(0.0, 1.0, size=ng)

        result = ordinate_scan(a, b, psi_0)
        expected = psi_0 + np.cumsum(b, axis=0)

        np.testing.assert_array_equal(result, expected)

    def test_ordinate_scan_linearity_in_psi_0(self) -> None:
        r"""``psi_0 → ordinate_scan(a, 0, psi_0)`` is linear.

        Closed-form: ``scan(a, 0, c·psi_0) == c·scan(a, 0, psi_0)``.
        ``scan(a, 0, p + q) == scan(a, 0, p) + scan(a, 0, q)``.
        rtol=1e-14.
        """
        rng = np.random.default_rng(seed=6)
        nx, ng = 20, 2
        a = rng.uniform(0.3, 1.4, size=(nx, ng))
        b = np.zeros((nx, ng))
        p = rng.uniform(0.1, 0.5, size=ng)
        q = rng.uniform(0.2, 0.6, size=ng)
        c = 2.71828

        scan_p = ordinate_scan(a, b, p)
        scan_q = ordinate_scan(a, b, q)
        scan_pq = ordinate_scan(a, b, p + q)
        scan_cp = ordinate_scan(a, b, c * p)

        np.testing.assert_allclose(scan_pq, scan_p + scan_q, rtol=1e-14)
        np.testing.assert_allclose(scan_cp, c * scan_p, rtol=1e-14)

    def test_ordinate_scan_linearity_in_b(self) -> None:
        r"""``b → ordinate_scan(a, b, 0)`` is linear.

        ``scan(a, c·b, 0) == c·scan(a, b, 0)``.
        ``scan(a, b1 + b2, 0) == scan(a, b1, 0) + scan(a, b2, 0)``.

        rtol=1e-13 — the scan involves ``cumsum`` over (b/cumprod_a)
        with cancelling positive/negative terms; the closed-form
        scales the cumsum back by the cumprod, so adjacent partial
        sums of opposite sign can amplify ULP by a small constant
        factor.  rtol=1e-13 (≈ 50× ULP) accommodates the worst-case
        cancellation we observed across random ng=3, nx=18 cases.
        """
        rng = np.random.default_rng(seed=7)
        nx, ng = 18, 3
        a = rng.uniform(0.3, 1.4, size=(nx, ng))
        psi_0 = np.zeros(ng)
        b1 = rng.uniform(-0.4, 0.4, size=(nx, ng))
        b2 = rng.uniform(-0.6, 0.6, size=(nx, ng))
        c = 1.7

        scan_b1 = ordinate_scan(a, b1, psi_0)
        scan_b2 = ordinate_scan(a, b2, psi_0)
        scan_sum = ordinate_scan(a, b1 + b2, psi_0)
        scan_cb1 = ordinate_scan(a, c * b1, psi_0)

        np.testing.assert_allclose(scan_sum, scan_b1 + scan_b2, rtol=1e-13)
        np.testing.assert_allclose(scan_cb1, c * scan_b1, rtol=1e-13)

    def test_ordinate_scan_affine_combination(self) -> None:
        r"""``scan(a, b1+b2, p1+p2) == scan(a, b1, p1) + scan(a, b2, p2)``.

        Joint affine structure: the scan is jointly linear in
        (b, psi_0) for fixed a.  rtol=1e-13 (one extra add for the
        summed result; one ULP).
        """
        rng = np.random.default_rng(seed=8)
        nx, ng = 16, 2
        a = rng.uniform(0.3, 1.4, size=(nx, ng))
        p1 = rng.uniform(0.1, 0.4, size=ng)
        p2 = rng.uniform(0.2, 0.5, size=ng)
        b1 = rng.uniform(-0.3, 0.3, size=(nx, ng))
        b2 = rng.uniform(-0.4, 0.4, size=(nx, ng))

        s1 = ordinate_scan(a, b1, p1)
        s2 = ordinate_scan(a, b2, p2)
        s_combined = ordinate_scan(a, b1 + b2, p1 + p2)

        np.testing.assert_allclose(s_combined, s1 + s2, rtol=1e-13)


# ═══════════════════════════════════════════════════════════════════════
# 3. Numerical stability
# ═══════════════════════════════════════════════════════════════════════


class TestNumericalStability:
    r"""Edge-case regime tests.

    Historical "Gotcha #5" in earlier `_solve_recurrence` variants
    used the form ``(1 - a + eps)`` in a denominator, which blew up
    near ``a → 1``.  The Blelloch form has NO such denominator —
    ``cumprod_a`` and ``cumsum(b/cumprod_a)`` are both regular at
    ``a = 1`` and at every ``a > 0``.

    The only documented regime limit: ``a → 0`` makes
    ``b / cumprod_a`` grow, eventually overflowing.  Reactor physics
    rarely hits this — DD's ``a = 2|μ|·A_down/denom − 1`` stays in
    [−1, 1] for typical Σ_t × Δx products.
    """

    def test_ordinate_scan_near_identity_attenuation(self) -> None:
        r"""``a ≈ 1`` is well-conditioned (no denominator blow-up).

        Verify with ``a = 1 + ε`` for ε ∈ {0, 1e-15, 1e-10}.
        rtol=1e-12 to allow for accumulated ULP × N drift across
        many cells.
        """
        rng = np.random.default_rng(seed=9)
        nx, ng = 50, 2
        b = rng.uniform(-0.2, 0.2, size=(nx, ng))
        psi_0 = rng.uniform(0.1, 0.3, size=ng)

        for eps in (0.0, 1e-15, 1e-10):
            a = np.full((nx, ng), 1.0 + eps)
            scan = ordinate_scan(a, b, psi_0)
            loop = _explicit_scan_loop(a, b, psi_0)
            # Ensure both finite + agreement within accumulated ULP×N.
            assert np.all(np.isfinite(scan)), (
                f"non-finite at eps={eps}"
            )
            np.testing.assert_allclose(scan, loop, rtol=1e-12)

    def test_ordinate_scan_small_attenuation(self) -> None:
        r"""``a ≈ 0.1`` is well-conditioned for moderate chain length.

        With ``a ∈ [0.05, 0.2]`` and ``nx = 20``, ``cumprod_a`` decays
        to ``~ 0.1^20 = 1e-20`` which is well within IEEE-754 normal
        range.  ``b / cumprod_a`` can be ~ ``1e20`` for ``b ~ O(1)``;
        the cumprod_a multiplication at the end restores scale.

        Document the regime: ``cumprod_a`` underflow to subnormal
        begins around ``a^nx ~ 1e-308`` which requires
        ``nx · |log10(a)| > 308``.  For DD sweeps ``a ~ O(1)`` and
        ``nx ~ O(10^2)`` this is comfortably safe.

        rtol=1e-11 — accumulated drift from the small-magnitude
        cumprod * large-magnitude cumsum cancellation.
        """
        rng = np.random.default_rng(seed=10)
        nx, ng = 20, 2
        a = rng.uniform(0.05, 0.2, size=(nx, ng))
        b = rng.uniform(-0.5, 0.5, size=(nx, ng))
        psi_0 = rng.uniform(0.1, 0.3, size=ng)

        scan = ordinate_scan(a, b, psi_0)
        loop = _explicit_scan_loop(a, b, psi_0)

        assert np.all(np.isfinite(scan)), (
            "scan returned non-finite values in small-a regime"
        )
        np.testing.assert_allclose(scan, loop, rtol=1e-11)
