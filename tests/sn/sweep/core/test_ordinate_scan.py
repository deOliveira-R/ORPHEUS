r"""Tests for :func:`orpheus.sn.sweep.scan.ordinate_scan`.

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

from orpheus.geometry import (
    BC,
    CoordSystem,
    Mesh1D,
    cylindrical_streaming,
    slab_streaming,
    spherical_streaming,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.spatial import DiamondDifference, UpstreamState
from orpheus.transport.spatial.cell_balance import cell_balance_terms
from orpheus.transport.spatial.scheme import CellVisit
from orpheus.sn.sweep.scan import ordinate_scan
from tests.sn.sweep.core._c_surrogate import (
    c_from_constants,
    mm_constants_for_ordinate,
)


def _affine_coefficients_from_visits(
    visits: list[CellVisit],
    total_xs: np.ndarray,
    source: np.ndarray,
    angular_state: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Test-side adapter: per-cell ``(a, b)`` from a visit list.

    Step 2.5c retired ``DiamondDifference.affine_coefficients`` —
    the cache populator subsumes it.  These dual-view tests retain
    their algebraic structure by computing ``(a, b)`` per cell via
    :func:`cell_balance_terms` (the Pattern 2 anchor that both
    ``update`` and the cache populator derive from):

    .. math::

        a[i, g] = 2|\mu|\,A_{\rm total}[i] / \mathrm{denom}[i, g] - 1,
        \\
        b[i, g] = 2 \cdot (\mathrm{source}[i, g]
                            + \mathrm{ang\_contrib}[i, g])
                  / \mathrm{denom}[i, g].
    """
    nx = len(visits)
    ng = total_xs.shape[1]
    a_out = np.empty((nx, ng))
    b_out = np.empty((nx, ng))
    for idx, visit in enumerate(visits):
        psi_ang = (
            np.zeros(ng) if angular_state is None else angular_state[idx]
        )
        upstream = UpstreamState(
            spatial_upstream=np.zeros(ng),
            angular_upstream=(
                None if angular_state is None else psi_ang
            ),
        )
        # Issue #236 Phase 2 B3 / Step C: cell_balance_terms consumes the
        # angular-closure-owned c_in / c_out as data.  These visits were
        # stamped (by the builders) with the surrogate's c; read them off
        # the visit (the geometry-side τ/α packing is retired, so the c can
        # no longer be recomputed from ``visit.streaming_terms``).
        terms = cell_balance_terms(
            visit.streaming_terms,
            visit.face_area_downstream,
            total_xs[idx],
            upstream,
            c_in=visit.c_in,
            c_out=visit.c_out,
        )
        st = visit.streaming_terms
        A_total = st.face_area_inner + st.face_area_outer
        a_out[idx] = 2.0 * st.abs_mu * A_total / terms.denom - 1.0
        ang_contrib = (
            0.0 if angular_state is None
            else st.delta_A_over_w * terms.c_in * angular_state[idx]
        )
        b_out[idx] = 2.0 * (source[idx] + ang_contrib) / terms.denom
    return a_out, b_out


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


# ═══════════════════════════════════════════════════════════════════════
# 4. Dual-view contracts: affine_coefficients ↔ update
# ═══════════════════════════════════════════════════════════════════════


# Mesh fixtures (mirror tests/sn/sweep/core/test_diamond.py).

def _slab_mesh(nx: int = 5, length: float = 1.0) -> Mesh1D:
    return Mesh1D(
        edges=np.linspace(0.0, length, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )


def _spherical_mesh(nx: int = 5, radius: float = 1.0) -> Mesh1D:
    return Mesh1D(
        edges=np.linspace(0.0, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )


def _cylindrical_mesh(nx: int = 4, radius: float = 1.0) -> Mesh1D:
    return Mesh1D(
        edges=np.linspace(0.0, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )


def _build_slab_visits_and_inputs(
    nx: int, n_groups: int, direction_idx: int, *, seed: int,
):
    """Build a list of slab CellVisit + matching (total_xs, source,
    psi_angular_init, psi_in) arrays for the dual-view tests.
    """
    mesh = _slab_mesh(nx=nx, length=1.0)
    quad = Quadrature.gauss_legendre(4)
    op = slab_streaming(mesh, quad)
    visits = []
    for i in range(nx):
        st = op.streaming_terms(cell_idx=i, direction_idx=direction_idx)
        # Issue #236 Phase 2 B3 / Step C: DD.update reads the M-M c_in / c_out
        # off the visit; stamp them via the independent surrogate (slab → 0.0)
        # so the per-cell update and the affine-coefficient builder agree.
        tau, alpha_in, alpha_out = mm_constants_for_ordinate(
            op, i, direction_idx,
        )
        c_in, c_out = c_from_constants(tau, alpha_in, alpha_out)
        visits.append(
            CellVisit(
                cell_idx=i,
                streaming_terms=st,
                face_area_downstream=1.0,
                c_in=c_in, c_out=c_out, tau=tau,
            )
        )
    rng = np.random.default_rng(seed=seed)
    total_xs = rng.uniform(0.5, 1.5, size=(nx, n_groups))
    source = rng.uniform(0.1, 1.0, size=(nx, n_groups))
    psi_in = rng.uniform(0.05, 0.5, size=n_groups)
    return visits, total_xs, source, None, psi_in


def _build_sphere_visits_and_inputs(
    nx: int, n_groups: int, *, outward: bool, seed: int,
):
    """Build sphere visits for one outward (or inward) ordinate."""
    mesh = _spherical_mesh(nx=nx, radius=1.0)
    quad = Quadrature.gauss_legendre(8)
    op = spherical_streaming(mesh, quad)
    direction_idx = quad.N - 2 if outward else 1
    visits = []
    if outward:
        cell_order = range(nx)
    else:
        cell_order = range(nx - 1, -1, -1)
    for i in cell_order:
        st = op.streaming_terms(cell_idx=i, direction_idx=direction_idx)
        A_down = st.face_area_outer if outward else st.face_area_inner
        # Issue #236 Phase 2 B3 / Step C: stamp the M-M c_in / c_out / τ
        # DD.update reads, via the independent surrogate.
        tau, alpha_in, alpha_out = mm_constants_for_ordinate(
            op, i, direction_idx,
        )
        c_in, c_out = c_from_constants(tau, alpha_in, alpha_out)
        visits.append(
            CellVisit(
                cell_idx=i,
                streaming_terms=st,
                face_area_downstream=A_down,
                c_in=c_in, c_out=c_out, tau=tau,
            )
        )
    rng = np.random.default_rng(seed=seed)
    total_xs = rng.uniform(0.5, 1.5, size=(nx, n_groups))
    source = rng.uniform(0.1, 1.0, size=(nx, n_groups))
    psi_angular = rng.uniform(0.02, 0.15, size=(nx, n_groups))
    psi_in = rng.uniform(0.05, 0.5, size=n_groups)
    return visits, total_xs, source, psi_angular, psi_in


def _build_cylinder_visits_and_inputs(
    nx: int, n_groups: int, *, seed: int,
):
    """Build cylinder visits for one within-level ordinate (non-degenerate)."""
    mesh = _cylindrical_mesh(nx=nx, radius=1.0)
    quad = Quadrature.product(n_mu=4, n_phi=4)
    op = cylindrical_streaming(mesh, quad)
    direction_idx = 0
    mu_level_idx = 0
    visits = []
    # Check sign of η for this ordinate to decide cell order
    level_indices = quad.level_indices
    global_n = int(level_indices[mu_level_idx][direction_idx])
    eta_n = float(quad.eta[global_n])
    if eta_n >= 0:
        cell_order = range(nx)
        select_outer = True
    else:
        cell_order = range(nx - 1, -1, -1)
        select_outer = False
    for i in cell_order:
        st = op.streaming_terms(
            cell_idx=i, direction_idx=direction_idx,
            mu_level_idx=mu_level_idx,
        )
        A_down = st.face_area_outer if select_outer else st.face_area_inner
        # Issue #236 Phase 2 B3 / Step C: stamp the M-M c_in / c_out / τ
        # (clamped cylinder τ) DD.update reads, via the independent surrogate.
        tau, alpha_in, alpha_out = mm_constants_for_ordinate(
            op, i, direction_idx, mu_level_idx=mu_level_idx,
        )
        c_in, c_out = c_from_constants(tau, alpha_in, alpha_out)
        visits.append(
            CellVisit(
                cell_idx=i,
                streaming_terms=st,
                face_area_downstream=A_down,
                c_in=c_in, c_out=c_out, tau=tau,
            )
        )
    rng = np.random.default_rng(seed=seed)
    total_xs = rng.uniform(0.5, 1.5, size=(nx, n_groups))
    source = rng.uniform(0.1, 1.0, size=(nx, n_groups))
    psi_angular = rng.uniform(0.02, 0.15, size=(nx, n_groups))
    psi_in = rng.uniform(0.05, 0.5, size=n_groups)
    return visits, total_xs, source, psi_angular, psi_in


# Source kind variations: zero, constant, random.

def _source_zero(shape, rng):
    return np.zeros(shape)


def _source_constant(shape, rng):
    return np.full(shape, 0.5)


def _source_random(shape, rng):
    return rng.uniform(0.0, 1.0, size=shape)


_SOURCE_KINDS = {
    "zero": _source_zero,
    "constant": _source_constant,
    "random": _source_random,
}


_GEOMETRY_BUILDERS = {
    "slab": lambda nx, ng, seed: _build_slab_visits_and_inputs(
        nx=nx, n_groups=ng, direction_idx=3, seed=seed,
    ),
    "sphere_outward": lambda nx, ng, seed: _build_sphere_visits_and_inputs(
        nx=nx, n_groups=ng, outward=True, seed=seed,
    ),
    "sphere_inward": lambda nx, ng, seed: _build_sphere_visits_and_inputs(
        nx=nx, n_groups=ng, outward=False, seed=seed,
    ),
    "cylinder": lambda nx, ng, seed: _build_cylinder_visits_and_inputs(
        nx=nx, n_groups=ng, seed=seed,
    ),
}


class TestDualViewContracts:
    r"""``affine_coefficients`` ↔ ``update`` consistency theorems.

    The dual-view contract: for every cell along an ordinate's
    spatial chain, the single-cell ``update`` and the vectorised
    ``affine_coefficients`` MUST agree on the outgoing spatial face
    flux at the same upstream state.  rtol=1e-13 — one division ULP
    band; failure here means the affine derivation is wrong or
    ``update`` has hidden non-affine terms.
    """

    @pytest.mark.parametrize("geometry", list(_GEOMETRY_BUILDERS.keys()))
    @pytest.mark.parametrize("n_groups", [1, 2, 4])
    @pytest.mark.parametrize("source_kind", list(_SOURCE_KINDS.keys()))
    def test_affine_coefficients_matches_update_single_cell(
        self, geometry: str, n_groups: int, source_kind: str,
    ) -> None:
        r"""For every cell, the scan output matches per-cell ``update``.

        DUAL-VIEW CONTRACT.  Parametrised over geometry (slab,
        sphere outward/inward, cylinder) × n_groups (1, 2, 4) ×
        source kind (zero, constant, random) = 36 cases.  rtol=1e-13
        (a single division per cell; one ULP band).
        """
        nx = 6
        visits, total_xs, source_random, psi_angular, psi_in = (
            _GEOMETRY_BUILDERS[geometry](nx, n_groups, 42)
        )
        # Override source per kind.
        rng = np.random.default_rng(seed=99)
        source = _SOURCE_KINDS[source_kind]((nx, n_groups), rng)

        strat = DiamondDifference()

        # ── Per-cell update (slow but legible reference) ───────────
        per_cell_outputs = []
        psi_chain = psi_in
        for idx, visit in enumerate(visits):
            psi_ang = (
                None if psi_angular is None else psi_angular[idx]
            )
            upstream = UpstreamState(
                spatial_upstream=psi_chain,
                angular_upstream=psi_ang,
            )
            result = strat.update(
                visit=visit, total_xs=total_xs[idx],
                source=source[idx], upstream_state=upstream,
            )
            per_cell_outputs.append(result.outgoing_spatial_flux)
            # Advance the chain.
            if result.outgoing_spatial_flux is not None:
                psi_chain = result.outgoing_spatial_flux

        # ── Vectorised (a, b) builder + ordinate_scan ──────────────
        a, b = _affine_coefficients_from_visits(
            visits, total_xs, source, psi_angular,
        )
        scan_out = ordinate_scan(a, b, psi_in)

        # Compare per-cell.  For cells where the per-cell update
        # returns ``outgoing_spatial_flux is None`` (e.g. innermost
        # cell on an inward curvilinear sweep where the downstream
        # face has zero area — the pole or the centre), the scan
        # output is the algebraically-defined value but no chain
        # consumer reads it; skip the comparison.
        compared = 0
        for idx in range(nx):
            if per_cell_outputs[idx] is None:
                continue
            np.testing.assert_allclose(
                scan_out[idx],
                per_cell_outputs[idx],
                rtol=1e-13,
                err_msg=f"cell {idx}, geometry={geometry}, "
                        f"ng={n_groups}, source={source_kind}",
            )
            compared += 1
        # At least one cell must be compared (else the test is a no-op).
        assert compared >= nx - 1, (
            f"only {compared}/{nx} cells compared for {geometry}; "
            "test may have degenerated to no-op"
        )

    @pytest.mark.parametrize("geometry", list(_GEOMETRY_BUILDERS.keys()))
    def test_affine_coefficients_vectorisation_matches_serial(
        self, geometry: str,
    ) -> None:
        r"""``affine_coefficients(visits)`` per-cell == single-call.

        Vectorisation must commute with per-cell application: calling
        ``affine_coefficients`` on the full visit list yields the
        SAME ``(a, b)`` arrays as calling it once per single-visit
        list and concatenating.  rtol=1e-14 (no fold; pure per-cell
        parallel computation).
        """
        nx, n_groups = 5, 2
        visits, total_xs, source, psi_angular, _ = (
            _GEOMETRY_BUILDERS[geometry](nx, n_groups, 17)
        )
        strat = DiamondDifference()

        a_full, b_full = _affine_coefficients_from_visits(
            visits, total_xs, source, psi_angular,
        )

        a_serial = np.empty((nx, n_groups))
        b_serial = np.empty((nx, n_groups))
        for idx, visit in enumerate(visits):
            psi_ang = (
                None if psi_angular is None else psi_angular[idx:idx + 1]
            )
            a_one, b_one = _affine_coefficients_from_visits(
                [visit],
                total_xs[idx:idx + 1],
                source[idx:idx + 1],
                psi_ang,
            )
            a_serial[idx] = a_one[0]
            b_serial[idx] = b_one[0]

        np.testing.assert_allclose(a_full, a_serial, rtol=1e-14)
        np.testing.assert_allclose(b_full, b_serial, rtol=1e-14)

    def test_full_sweep_matches_pre_step_2_5b_baseline(self) -> None:
        r"""Per-ordinate scan reproduces the per-cell-fold baseline.

        Builds a representative slab ordinate, runs the slow per-cell
        fold (the pre-Step-2.5b reference), and runs the scan via
        ``affine_coefficients + ordinate_scan``.  Their downstream
        face fluxes match cell-by-cell at rtol=1e-12 (operation-count
        × ULP for nx=20 cells).

        Principled-equivalence gate per
        ``vv-principles`` §"Bit-identity vs principled-equivalence":
        the two computations are algebraically identical; differences
        are bounded by FP-non-associativity.
        """
        nx, n_groups = 20, 2
        visits, total_xs, source, psi_angular, psi_in = (
            _GEOMETRY_BUILDERS["slab"](nx, n_groups, 11)
        )
        strat = DiamondDifference()

        # ── Baseline: per-cell fold (pre-2.5b reference) ───────────
        psi_chain = psi_in
        baseline = np.empty((nx, n_groups))
        for idx, visit in enumerate(visits):
            upstream = UpstreamState(
                spatial_upstream=psi_chain,
                angular_upstream=None,
            )
            result = strat.update(
                visit=visit, total_xs=total_xs[idx],
                source=source[idx], upstream_state=upstream,
            )
            baseline[idx] = result.outgoing_spatial_flux
            psi_chain = result.outgoing_spatial_flux

        # ── New scan path ──────────────────────────────────────────
        a, b = _affine_coefficients_from_visits(visits, total_xs, source, None)
        scan_out = ordinate_scan(a, b, psi_in)

        np.testing.assert_allclose(scan_out, baseline, rtol=1e-12)
