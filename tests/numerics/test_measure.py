r"""Verification gates for :mod:`orpheus.numerics.measure`.

The polynomial-degree-of-exactness tests are tagged ``L1`` because
they verify the closed-form analytical claim made by Gauss quadrature
theorems (Stoer & Bulirsch 2002, Theorem 3.6.20 — N-point Gauss
rule is exact for every polynomial of degree :math:`\le 2N - 1`).
The reference value is the analytical integral computed in closed
form, NOT the output of another numerical integrator.

The remaining tests are tagged ``foundation`` because they verify
software invariants of the composition algebra (tensor product,
direct sum, pushforward, restriction, metadata propagation, bundle
disintegration) that have no Sphinx ``:label:`` of their own — they
guard the data structure contract that downstream consumers rely on:
:class:`~orpheus.numerics.quadrature.Quadrature`, which wraps a
``DiscreteMeasure`` and exposes the four SN rule families as named
classmethod factories, and (in Wave 2) the MoC bundle measures. The
Issue-4 design that this docstring named — an ``AngularQuadrature``
Protocol with four per-family adapter classes at
``orpheus.sn.quadrature`` — was retired into those factories in R-1
Phase A detour-C; the contract guarded here did not change.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.measure import BundleMeasure, DiscreteMeasure, equispaced, gauss_chebyshev, gauss_legendre
from orpheus.numerics.manifold import (
    COSINE_INTERVAL,
    IndexSet,
    Interval,
    ManifoldMap,
    REAL_LINE,
    SPHERE,
    UNIT_INTERVAL,
)
from orpheus.numerics.symmetry import SubgroupOfO3


# ============================================================================
# L1 — polynomial-exactness against published Gauss-quadrature theorems
# ============================================================================


@pytest.mark.l1
@pytest.mark.verifies("discrete-measure-integrate")
@pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6, 8, 10])
def test_gauss_legendre_polynomial_exactness(n: int) -> None:
    r"""Stoer-Bulirsch Theorem 3.6.20: N-point Gauss-Legendre is exact
    for every polynomial of degree :math:`\le 2N - 1` on
    :math:`[-1, 1]`.

    Reference values: closed-form analytical integrals
    :math:`\int_{-1}^{1} x^k \, dx = \frac{1 - (-1)^{k+1}}{k + 1}`,
    which is :math:`2/(k+1)` for even :math:`k` and 0 for odd.
    """
    rule = gauss_legendre(n)
    # The rule advertises its degree of exactness; the test asserts
    # the advertised claim.
    assert rule.degree_of_exactness == 2 * n - 1
    for k in range(2 * n):  # k = 0, 1, ..., 2n-1
        analytical = 2.0 / (k + 1) if k % 2 == 0 else 0.0
        numerical = rule.integrate(lambda x, kk=k: x ** kk)
        assert np.isclose(numerical, analytical, atol=1e-12, rtol=1e-12), (
            f"Gauss-Legendre n={n} should be exact for x^{k}: got "
            f"{numerical}, analytical {analytical}"
        )


@pytest.mark.l1
@pytest.mark.verifies("discrete-measure-integrate")
@pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 8])
def test_gauss_chebyshev_polynomial_exactness(n: int) -> None:
    r"""Stoer-Bulirsch §3.6: N-point Gauss-Chebyshev (first kind) is
    exact for the weighted integral
    :math:`\int_{-1}^{1} q(x) (1 - x^2)^{-1/2} \, dx` for every
    polynomial :math:`q` of degree :math:`\le 2N - 1`.

    The weighted integrals of monomials are

    .. math::
        \int_{-1}^{1} \frac{x^k}{\sqrt{1 - x^2}} \, dx
        = \begin{cases}
            \pi & k = 0 \\
            0 & k \text{ odd} \\
            \pi \cdot \frac{(k-1)!!}{k!!} & k \ge 2 \text{ even}
        \end{cases}

    where :math:`k!!` is the double factorial.
    """
    rule = gauss_chebyshev(n)
    assert rule.degree_of_exactness == 2 * n - 1

    for k in range(2 * n):
        if k == 0:
            analytical = np.pi
        elif k % 2 == 1:
            analytical = 0.0
        else:
            # (k-1)!! / k!! · π for even k >= 2.
            # Compute double factorials directly.
            num = 1
            for j in range(k - 1, 0, -2):
                num *= j
            den = 1
            for j in range(k, 0, -2):
                den *= j
            analytical = np.pi * num / den
        numerical = rule.integrate(lambda x, kk=k: x ** kk)
        assert np.isclose(numerical, analytical, atol=1e-12, rtol=1e-12), (
            f"Gauss-Chebyshev n={n} should be exact for x^{k}: got "
            f"{numerical}, analytical {analytical}"
        )


@pytest.mark.l1
@pytest.mark.verifies("discrete-measure-integrate")
def test_gauss_legendre_weight_sum_is_two() -> None:
    """:math:`\\sum_i w_i = 2` for any :math:`n` — the integral of
    the constant function 1 on :math:`[-1, 1]`."""
    for n in (1, 2, 4, 8, 16, 32):
        rule = gauss_legendre(n)
        assert np.isclose(rule.weights.sum(), 2.0, atol=1e-14)


@pytest.mark.l1
@pytest.mark.verifies("discrete-measure-integrate")
def test_gauss_chebyshev_weight_sum_is_pi() -> None:
    """:math:`\\sum_i w_i = \\pi` for any :math:`n` — the weighted
    integral of the constant function 1 on :math:`[-1, 1]` against
    :math:`(1 - x^2)^{-1/2}`."""
    for n in (1, 2, 4, 8, 16, 32):
        rule = gauss_chebyshev(n)
        assert np.isclose(rule.weights.sum(), np.pi, atol=1e-14)


# ============================================================================
# Foundation — software invariants of the composition algebra
# ============================================================================


@pytest.mark.foundation
def test_n_points_property() -> None:
    """``n_points`` reads ``len(weights)``."""
    rule = gauss_legendre(7)
    assert rule.n_points == 7
    assert rule.n_points == rule.weights.shape[0]


@pytest.mark.foundation
def test_dim_property_1d_and_2d() -> None:
    """``dim`` is 1 for ``(N,)`` nodes, 2 for ``(N, 2)`` nodes."""
    one_d = gauss_legendre(5)
    assert one_d.dim == 1

    two_d = gauss_legendre(3) * gauss_legendre(2)
    assert two_d.dim == 2


@pytest.mark.foundation
def test_integrate_vector_valued_integrand() -> None:
    """For a vector-valued ``f`` returning ``(N, k)``, the result
    contracts the leading axis and has shape ``(k,)``."""
    rule = gauss_legendre(8)
    # Vector integrand: returns (1, x, x^2) at each node.
    result = rule.integrate(lambda x: np.stack([np.ones_like(x), x, x ** 2], axis=1))
    expected = np.array([2.0, 0.0, 2.0 / 3.0])
    assert result.shape == (3,)
    assert np.allclose(result, expected, atol=1e-12)


@pytest.mark.foundation
def test_tensor_product_structure() -> None:
    """``μ * ν`` produces ``(n1*n2, d1+d2)`` nodes and
    outer-producted weights, in ``ij``-meshgrid order (outer loop μ,
    inner loop ν)."""
    a = gauss_legendre(3)
    b = gauss_legendre(4)
    p = a * b
    assert p.n_points == 12
    assert p.nodes.shape == (12, 2)
    assert p.weights.shape == (12,)
    # Cartesian-product order: column 0 advances slowly, column 1 fast.
    assert np.allclose(p.nodes[:4, 0], a.nodes[0])  # first 4: μ_a[0]
    assert np.allclose(p.nodes[:4, 1], b.nodes)
    assert np.allclose(p.nodes[4:8, 0], a.nodes[1])
    # Weight is outer product: w_ij = w_a[i] * w_b[j].
    assert np.isclose(p.weights[5], a.weights[1] * b.weights[1])


@pytest.mark.foundation
def test_tensor_product_separable_integrand_factors() -> None:
    r"""For separable :math:`f(x, y) = g(x) h(y)`,
    :math:`\int f \, d(\mu \otimes \nu) = (\int g \, d\mu)(\int h \,
    d\nu)`. This is the discrete Fubini-Tonelli identity and is the
    key test from the spec."""
    a = gauss_legendre(3)
    b = gauss_legendre(4)
    p = a * b

    # f(x, y) = x * y; ∫_-1^1 x dx = 0, so the product integral is 0.
    int_xy = p.integrate(lambda xy: xy[:, 0] * xy[:, 1])
    int_x = a.integrate(lambda x: x)
    int_y = b.integrate(lambda y: y)
    assert np.isclose(int_xy, int_x * int_y, atol=1e-13)
    assert np.isclose(int_xy, 0.0, atol=1e-13)

    # g(x) = x^2, h(y) = 1; product = ∫ x^2 dx · ∫ 1 dy = (2/3) * 2 = 4/3.
    int_gh = p.integrate(lambda xy: (xy[:, 0] ** 2) * np.ones_like(xy[:, 1]))
    int_g = a.integrate(lambda x: x ** 2)
    int_h = b.integrate(lambda y: np.ones_like(y))
    assert np.isclose(int_gh, int_g * int_h, atol=1e-13)
    assert np.isclose(int_gh, 4.0 / 3.0, atol=1e-13)


@pytest.mark.foundation
def test_tensor_product_metadata_propagation() -> None:
    """``space`` concatenates with ``×``; ``degree_of_exactness`` is
    the min of the two factors; ``invariance_group`` is dropped."""
    a = gauss_legendre(3)  # dx = 5
    b = gauss_legendre(4)  # dx = 7
    p = a * b
    assert p.support == COSINE_INTERVAL * COSINE_INTERVAL
    assert p.degree_of_exactness == 5
    assert p.invariance_group is None


@pytest.mark.foundation
def test_direct_sum_concatenates_on_shared_space() -> None:
    """``μ + ν`` requires equal ``space``; concatenates atoms."""
    a = equispaced(0.0, 0.5, 4)
    b = equispaced(0.0, 0.5, 6)
    s = a + b
    assert s.n_points == 10
    assert s.support == a.support
    assert np.allclose(s.weights[:4], a.weights)
    assert np.allclose(s.weights[4:], b.weights)


@pytest.mark.foundation
def test_direct_sum_rejects_mismatched_spaces() -> None:
    """``μ + ν`` raises ``ValueError`` on mismatched spaces."""
    a = gauss_legendre(3)
    b = equispaced(0.0, 1.0, 3)
    with pytest.raises(ValueError, match="equal spaces"):
        _ = a + b


@pytest.mark.foundation
def test_direct_sum_composite_quadrature() -> None:
    r"""Composite Gauss-Legendre on :math:`[-1, 1]` built as direct
    sum of two affinely-mapped panels integrates polynomials of
    moderate degree exactly when each panel handles its own degree.

    This stresses the direct-sum + pushforward composition: the
    'composite' rule is built explicitly from primitives.
    """
    # Build N-point GL on [-1, 0] and [0, 1] via pushforward, then
    # direct-sum on the common space "[-1,1]".
    base_left = gauss_legendre(4).pushforward(
        ManifoldMap(COSINE_INTERVAL, COSINE_INTERVAL, lambda x: 0.5 * (x - 1.0)),
    )
    # The pushforward x -> (x-1)/2 maps [-1,1] to [-1, 0]. The Jacobian
    # (1/2) MUST be applied to the weights manually — the pushforward
    # is φ-image, not Radon-Nikodym (per docstring).
    base_left = DiscreteMeasure(
        nodes=base_left.nodes,
        weights=base_left.weights * 0.5,
        support=base_left.support,
    )
    base_right = gauss_legendre(4).pushforward(
        ManifoldMap(COSINE_INTERVAL, COSINE_INTERVAL, lambda x: 0.5 * (x + 1.0)),
    )
    base_right = DiscreteMeasure(
        nodes=base_right.nodes,
        weights=base_right.weights * 0.5,
        support=base_right.support,
    )
    composite = base_left + base_right

    # Each panel has dx=7; the composite is exact for x^k with k <= 7
    # only (the panels are GL(4) on each half).
    for k in range(8):
        analytical = 2.0 / (k + 1) if k % 2 == 0 else 0.0
        numerical = composite.integrate(lambda x, kk=k: x ** kk)
        assert np.isclose(numerical, analytical, atol=1e-13), (
            f"composite integration failed at k={k}: {numerical} vs {analytical}"
        )


@pytest.mark.foundation
def test_pushforward_invertible_map_change_of_variables() -> None:
    r"""For the equispaced midpoint measure on :math:`[0, 1]` and
    :math:`\varphi(x) = x^2` (invertible on that interval), the
    image measure has nodes :math:`\varphi(x_i) = x_i^2` and the
    same weights. Integrating :math:`x^2` against the pushed
    measure should give :math:`\int_0^1 x^4 \, dx = 1/5` — this is
    the change-of-variables identity :eq:`discrete-measure-pushforward`.
    """
    n = 200  # midpoint rule needs many points to reach 1/5
    base = equispaced(0.0, 1.0, n)
    pushed = base.pushforward(
        ManifoldMap(UNIT_INTERVAL, UNIT_INTERVAL, lambda x: x ** 2)
    )
    assert pushed.n_points == n
    assert np.allclose(pushed.weights, base.weights)
    assert np.allclose(pushed.nodes, base.nodes ** 2)

    # ∫_pushed x^2 dν = ∫_base φ(x)^2 dμ = ∫_0^1 x^4 dx = 1/5.
    numerical = pushed.integrate(lambda x: x ** 2)
    expected = 1.0 / 5.0
    # Midpoint rule has O(h^2) convergence; with n=200 we expect ~3e-5.
    assert np.isclose(numerical, expected, atol=1e-3)


@pytest.mark.foundation
def test_pushforward_drops_metadata() -> None:
    """Pushforward drops ``invariance_group`` and
    ``degree_of_exactness`` (φ may not preserve either)."""
    base = gauss_legendre(5).with_metadata(invariance_group=SubgroupOfO3.O3)
    assert base.invariance_group == SubgroupOfO3.O3
    assert base.degree_of_exactness == 9
    pushed = base.pushforward(
        ManifoldMap(COSINE_INTERVAL, IndexSet(label="img"), lambda x: x ** 3)
    )
    assert pushed.invariance_group is None
    assert pushed.degree_of_exactness is None


@pytest.mark.foundation
def test_pushforward_READS_its_target_off_the_map_and_refuses_a_map_out_of_elsewhere() -> None:
    r"""``pushforward`` lands where the MAP says, and only along a map out of
    the measure's own support.

    ⭐ Three states of one verb, in campaign order (#429).  Until 2026-09-01
    ``new_space`` was optional and defaulted to a FABRICATED support named
    ``f"φ_*({self.support})"`` — and this test pinned that name.  Tracker
    2.0c made the caller NAME the target (this test then pinned the refusal
    of a bare callable).  Tracker 2.3 (2026-09-02) makes the map CARRY it —
    a :class:`~orpheus.numerics.manifold.ManifoldMap` with ``domain`` and
    ``codomain`` — so the target cannot be forged at the call site at all:
    ERR-080 is the barycentre map applied with its codomain declared as
    :math:`S^2`, and through this verb that sentence has no spelling.

    Three legs (vv-principles #11): the positive one, and two refusals — a
    map out of the WRONG point set (the same numbers on a different
    manifold), and the retired call-site spelling of a target.
    """
    base = gauss_legendre(3)  # on COSINE_INTERVAL
    shift = ManifoldMap(COSINE_INTERVAL, Interval(0.0, 2.0), lambda x: x + 1.0)

    pushed = base.pushforward(shift)
    assert pushed.support == Interval(0.0, 2.0)
    assert pushed.support is shift.codomain  # READ, not re-derived
    assert np.array_equal(pushed.nodes, base.nodes + 1.0)

    # The same arithmetic, out of a point set this measure does not live on.
    elsewhere = ManifoldMap(UNIT_INTERVAL, Interval(0.0, 2.0), lambda x: x + 1.0)
    with pytest.raises(ValueError, match="map's domain"):
        _ = base.pushforward(elsewhere)

    # The 2.0c spelling — a bare callable plus a target named at the call
    # site — is gone, not merely discouraged.
    with pytest.raises(TypeError, match="new_space"):
        _ = base.pushforward(lambda x: x + 1.0, new_space=Interval(0.0, 2.0))  # type: ignore


@pytest.mark.foundation
def test_restrict_to_positive_subset() -> None:
    """``gauss_legendre(10).restrict(x > 0)`` keeps the 5 positive
    nodes (GL with even N has nodes symmetric about zero)."""
    rule = gauss_legendre(10)
    positive = rule.restrict(lambda x: x > 0)
    assert positive.n_points == 5
    assert np.all(positive.nodes > 0)
    assert np.allclose(positive.weights, rule.weights[rule.nodes > 0])


@pytest.mark.foundation
def test_restrict_drops_metadata() -> None:
    """Restriction drops ``invariance_group`` and
    ``degree_of_exactness`` (the restricted rule is not a Gauss
    rule on the cut domain)."""
    rule = gauss_legendre(8).with_metadata(invariance_group=SubgroupOfO3.O3)
    half = rule.restrict(lambda x: x >= 0)
    assert half.invariance_group is None
    assert half.degree_of_exactness is None


@pytest.mark.foundation
def test_restrict_invalid_predicate_shape_raises() -> None:
    """Restriction predicate must return ``(N,)``-shaped boolean."""
    rule = gauss_legendre(4)
    with pytest.raises(ValueError, match="predicate"):
        _ = rule.restrict(lambda x: np.array([[True, False], [True, True]]))


@pytest.mark.foundation
def test_construction_invariants() -> None:
    """``DiscreteMeasure`` validates shapes at construction."""
    with pytest.raises(ValueError, match="weights must be 1-D"):
        DiscreteMeasure(
            nodes=np.array([1.0, 2.0]),
            weights=np.array([[1.0], [1.0]]),
            support=REAL_LINE,
        )
    with pytest.raises(ValueError, match="disagree on N"):
        DiscreteMeasure(
            nodes=np.array([1.0, 2.0, 3.0]),
            weights=np.array([1.0, 1.0]),
            support=REAL_LINE,
        )
    with pytest.raises(ValueError, match="nodes must be 1-D"):
        DiscreteMeasure(
            nodes=np.zeros((2, 2, 2)),
            weights=np.array([1.0, 1.0]),
            support=REAL_LINE,
        )


@pytest.mark.foundation
def test_with_metadata_returns_new_instance() -> None:
    """``with_metadata`` returns a new frozen instance with selected
    fields updated; the original is unchanged."""
    base = gauss_legendre(3)
    tagged = base.with_metadata(invariance_group=SubgroupOfO3.O3)
    assert tagged.invariance_group == SubgroupOfO3.O3
    assert base.invariance_group is None
    assert tagged is not base


# ============================================================================
# BundleMeasure smoke test
# ============================================================================


@pytest.mark.foundation
def test_bundle_measure_separable_smoke() -> None:
    r"""Bundle disintegration: build a base measure on ``[-1, 1]``
    and a fiber-measure factory that returns ``[0, 1]`` (independent
    of the base point — i.e. the bundle is trivially a product).
    Integrate :math:`f(b, x) = b^2 x` and confirm equality with the
    explicit nested integration.

    SN does not consume :class:`BundleMeasure` directly in Wave A,
    but this smoke test exercises the disintegration plumbing that
    MoC will rely on in Wave 2.
    """
    base = gauss_legendre(6)

    # Trivial fiber: same measure for every base point.
    fiber_template = equispaced(0.0, 1.0, 50)

    def fiber_factory(b: np.ndarray) -> DiscreteMeasure:
        # b is a scalar (base.dim == 1); independent of b in this smoke.
        return fiber_template

    bundle = BundleMeasure(base=base, fiber_factory=fiber_factory)

    # f(b, x) = b^2 * x; iterated integral
    # = (∫_-1^1 b^2 db) · (∫_0^1 x dx) = (2/3) · (1/2) = 1/3.
    result = bundle.integrate(lambda b, x: (b ** 2) * x)
    assert np.isclose(result, 1.0 / 3.0, atol=2e-3)  # midpoint floor

    # Same answer via explicit nested integration.
    inner_integrals = np.array(
        [fiber_template.integrate(lambda x, bb=b: (bb ** 2) * x) for b in base.nodes]
    )
    explicit = float(np.dot(base.weights, inner_integrals))
    assert np.isclose(result, explicit, atol=1e-13)


@pytest.mark.foundation
def test_bundle_measure_varying_fiber() -> None:
    r"""When the fiber depends non-trivially on the base point, the
    bundle integral is not separable — but it must still match the
    explicit nested computation. This pins the per-base-point
    dispatch path used by MoC ray-bundle quadratures.
    """
    base = gauss_legendre(4)

    # Fiber on [0, b+2] with 30 midpoint nodes — depends on b.
    def fiber_factory(b: np.ndarray) -> DiscreteMeasure:
        b_scalar = float(b)
        return equispaced(0.0, b_scalar + 2.0, 30)

    bundle = BundleMeasure(base=base, fiber_factory=fiber_factory)

    # f(b, x) = b * x; inner ∫_0^{b+2} b·x dx = b · (b+2)^2 / 2.
    # Outer ∫_-1^1 b · (b+2)^2 / 2 db evaluates analytically;
    # GL(4) is exact for polynomials of degree ≤ 7.
    result = bundle.integrate(lambda b, x: b * x)

    # Explicit nested integration — must match exactly (same algorithm).
    inner = np.array(
        [
            fiber_factory(b).integrate(lambda x, bb=b: bb * x)
            for b in base.nodes
        ]
    )
    explicit = float(np.dot(base.weights, inner))

    assert np.isclose(result, explicit, atol=1e-13)


# ============================================================================
# Equispaced primitive
# ============================================================================


@pytest.mark.foundation
def test_equispaced_midpoint_integrates_constants_exactly() -> None:
    """Midpoint rule integrates :math:`\\int_a^b 1 \\, dx = b - a`
    exactly for any :math:`n`."""
    for a, b, n in [(0.0, 1.0, 1), (-1.0, 1.0, 4), (0.0, np.pi, 17)]:
        rule = equispaced(a, b, n)
        result = rule.integrate(lambda x: np.ones_like(x))
        assert np.isclose(result, b - a, atol=1e-14)


@pytest.mark.foundation
def test_equispaced_midpoint_integrates_linears_exactly() -> None:
    """Midpoint rule has degree of exactness 1: integrates linears
    exactly. :math:`\\int_a^b x \\, dx = (b^2 - a^2)/2`."""
    for a, b, n in [(0.0, 1.0, 5), (-1.0, 1.0, 6), (0.0, np.pi, 12)]:
        rule = equispaced(a, b, n)
        result = rule.integrate(lambda x: x)
        expected = 0.5 * (b ** 2 - a ** 2)
        assert np.isclose(result, expected, atol=1e-13)


@pytest.mark.foundation
def test_equispaced_invalid_args_raise() -> None:
    """``equispaced`` rejects ``n < 1`` and ``a >= b``."""
    with pytest.raises(ValueError, match="n >= 1"):
        equispaced(0.0, 1.0, 0)
    with pytest.raises(ValueError, match="a < b"):
        equispaced(1.0, 0.0, 5)


@pytest.mark.foundation
def test_constructor_invalid_args_raise() -> None:
    """1-D primitive constructors reject ``n < 1``."""
    with pytest.raises(ValueError, match="n >= 1"):
        gauss_legendre(0)
    with pytest.raises(ValueError, match="n >= 1"):
        gauss_chebyshev(0)


# ─────────────────────────────────────────────────────────────────────
# Issue 9.6 dunder ergonomics + array-overload integrate
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_call_aliases_integrate() -> None:
    """``mu(f)`` is the same as ``mu.integrate(f)`` for callables."""
    rule = gauss_legendre(8)
    result_call = rule(lambda x: x ** 2)
    result_int = rule.integrate(lambda x: x ** 2)
    assert np.isclose(result_call, result_int, atol=1e-14)


@pytest.mark.foundation
def test_integrate_accepts_array_overload() -> None:
    """``mu.integrate(values)`` returns ``np.dot(weights, values)``.

    Load-bearing for ``mesh.volume_measure(scalar_flux)`` at
    production integration sites where the values are pre-computed.
    """
    rule = gauss_legendre(8)
    values = np.cos(rule.nodes)
    result = rule.integrate(values)
    expected = np.dot(rule.weights, values)
    assert np.isclose(result, expected, atol=1e-14)


@pytest.mark.foundation
def test_integrate_array_overload_matches_callable() -> None:
    """Array overload returns the same value as the callable form."""
    rule = gauss_legendre(8)
    f = lambda x: np.cos(x)
    result_callable = rule.integrate(f)
    result_array = rule.integrate(f(rule.nodes))
    assert np.isclose(result_callable, result_array, atol=1e-14)


@pytest.mark.foundation
def test_call_accepts_array() -> None:
    """``mu(values)`` (call) also accepts pre-evaluated arrays."""
    rule = gauss_legendre(8)
    values = np.cos(rule.nodes)
    result = rule(values)
    expected = np.dot(rule.weights, values)
    assert np.isclose(result, expected, atol=1e-14)


@pytest.mark.foundation
def test_iter_yields_node_weight_pairs() -> None:
    """``iter(mu)`` yields (node, weight) tuples in order."""
    rule = gauss_legendre(4)
    pairs = list(iter(rule))
    assert len(pairs) == 4
    for i, (node, weight) in enumerate(pairs):
        assert node == rule.nodes[i]
        assert weight == rule.weights[i]


@pytest.mark.foundation
def test_len_equals_n_points() -> None:
    rule = gauss_legendre(7)
    assert len(rule) == rule.n_points == 7


@pytest.mark.foundation
def test_getitem_returns_node_weight_tuple() -> None:
    rule = gauss_legendre(5)
    node, weight = rule[2]
    assert node == rule.nodes[2]
    assert weight == rule.weights[2]


@pytest.mark.foundation
def test_repr_smoke() -> None:
    rule = gauss_legendre(8)
    r = repr(rule)
    assert "DiscreteMeasure" in r
    assert "n_points=8" in r


@pytest.mark.foundation
def test_repr_with_invariance_group() -> None:
    """repr surfaces invariance_group when set."""
    rule = gauss_legendre(8).with_metadata(invariance_group=SubgroupOfO3.Mirror("z"))
    r = repr(rule)
    assert "invariance_group" in r
    assert "SubgroupOfO3.Mirror('z')" in r  # the typed group's repr, incl. its PLANE


# ============================================================================
# consolidate() — the reduction pushforward documents but does not perform
# ============================================================================


@pytest.mark.foundation
def test_consolidate_merges_coincident_atoms_preserving_mass() -> None:
    """Duplicate nodes collapse; the weights are summed, never dropped.

    Mass preservation is the in-tree discriminator between a QUOTIENT and a
    restriction — ``restrict`` drops mass, this must not.
    """
    mu = DiscreteMeasure(
        nodes=np.array([0.0, 1.0, 0.0, 2.0, 1.0]),
        weights=np.array([0.5, 0.25, 0.5, 1.0, 0.75]),
        support=COSINE_INTERVAL,
    )
    out = mu.consolidate()
    assert out.n_points == 3
    assert float(out.weights.sum()) == pytest.approx(float(mu.weights.sum()))
    # First-appearance order, weights summed per position.
    np.testing.assert_allclose(out.nodes, [0.0, 1.0, 2.0])
    np.testing.assert_allclose(out.weights, [1.0, 1.0, 1.0])


@pytest.mark.foundation
def test_consolidate_is_idempotent_and_a_no_op_when_already_reduced() -> None:
    """Consolidating twice equals consolidating once."""
    mu = DiscreteMeasure(
        nodes=np.array([0.0, 1.0, 0.0]),
        weights=np.array([1.0, 2.0, 3.0]),
        support=COSINE_INTERVAL,
    )
    once = mu.consolidate()
    twice = once.consolidate()
    np.testing.assert_array_equal(once.nodes, twice.nodes)
    np.testing.assert_array_equal(once.weights, twice.weights)

    already = DiscreteMeasure(
        nodes=np.array([0.0, 1.0]), weights=np.array([1.0, 2.0]),
        support=COSINE_INTERVAL,
    )
    assert already.consolidate() is already


@pytest.mark.foundation
def test_consolidate_preserves_the_claims_pushforward_drops() -> None:
    r"""``invariance_group`` and ``degree_of_exactness`` survive.

    Consolidation moves no node and changes no integral, so neither claim
    can be invalidated by it — unlike :meth:`pushforward`, where an
    arbitrary :math:`\varphi` can invalidate both.
    """
    from orpheus.numerics.exactness import ExactnessClaim
    from orpheus.numerics.generating_measure import LEGENDRE
    from orpheus.numerics.symmetry import SubgroupOfO3

    mu = DiscreteMeasure(
        nodes=np.array([-1.0, 0.0, 0.0, 1.0]),
        weights=np.array([1.0, 0.5, 0.5, 1.0]),
        support=COSINE_INTERVAL,
        invariance_group=SubgroupOfO3.Mirror("z"),
        exactness=ExactnessClaim(reference=LEGENDRE, degree=3),
    )
    out = mu.consolidate()
    assert out.n_points == 3
    assert out.invariance_group == SubgroupOfO3.Mirror("z")
    assert out.degree_of_exactness == 3
    assert out.support == mu.support


@pytest.mark.foundation
def test_quotient_is_pushforward_then_consolidate_with_orbifold_weights() -> None:
    r"""``quotient(G) = pushforward(rep).consolidate()``, and the weights are DERIVED.

    Folding a product rule by :math:`\langle \sigma_y \rangle` must give the
    orbit-stabilizer weights :math:`W = w \cdot |G| / |\mathrm{Stab}|` — not
    a chosen convention but a consequence. With :math:`|G| = 2`, every
    folded atom carries either its parent weight (it sits ON the mirror,
    :math:`|\mathrm{Stab}| = 2`) or twice it (a free orbit).

    The orbit count is Burnside's :math:`(N + F)/2`.
    """
    from orpheus.numerics.quadrature import Quadrature

    q = Quadrature.product(n_mu=4, n_phi=8)
    mu = q.measure

    def to_representative(nodes: np.ndarray) -> np.ndarray:
        out = nodes.copy()
        out[:, 1] = np.abs(out[:, 1])  # xi -> |xi| picks the orbit rep
        return out

    folded = mu.pushforward(
        ManifoldMap(SPHERE, SPHERE.quotient(SubgroupOfO3.Mirror("y")), to_representative),
    ).consolidate()

    n_fixed = int((np.abs(q.mu_y) < 1e-14).sum())
    assert folded.n_points == (mu.n_points + n_fixed) // 2  # Burnside
    assert float(folded.weights.sum()) == pytest.approx(
        float(mu.weights.sum()), rel=1e-12
    )

    moved = to_representative(mu.nodes)
    ratios = []
    for k in range(folded.n_points):
        members = np.flatnonzero(
            np.linalg.norm(moved - folded.nodes[k], axis=1) < 1e-12
        )
        parent = mu.weights[members]
        assert np.allclose(parent, parent[0]), "orbit members must share a weight"
        assert float(folded.weights[k]) == pytest.approx(float(parent.sum()))
        ratios.append(round(float(folded.weights[k] / parent[0]), 12))

    assert sorted(set(ratios)) == [1.0, 2.0], sorted(set(ratios))
    assert sum(1 for r in ratios if r == 1.0) == n_fixed


# ============================================================================
# quotient(G) — the verb naming the composite (Q5.1)
# ============================================================================
#
# The test above spells the fold BY HAND with the geometric section ξ → |ξ|
# (it predates the verb; it stays as the structurally-independent reference
# realization). The verb uses the only group-generic section — the orbit's
# first-appearing member — so the two produce the same quotient MEASURE
# through different representatives. The gates below are the plan's own:
# mass, orbit-stabilizer weights, idempotence, and the free-action
# certificate (Σ = ∅ ⟹ every orbit has length |G|).


def _fold_ring(n_phi: int = 8) -> DiscreteMeasure:
    """A σ_y-FREE fixture: one z-level, azimuthal MIDPOINT nodes.

    ``φ_k = (2k+1)π/n_phi`` never hits ``sin φ = 0``, so no node lies on
    the ``y = 0`` mirror — the action of ``⟨σ_y⟩`` is free (T24's
    admissibility condition, hand-built because the shipped product rule
    is equispaced-LEFT and does place nodes on the mirror).
    """
    z = 0.5
    r = np.sqrt(1.0 - z * z)
    phi = (2.0 * np.arange(n_phi) + 1.0) * np.pi / n_phi
    nodes = np.column_stack([r * np.cos(phi), r * np.sin(phi), np.full(n_phi, z)])
    return DiscreteMeasure(
        nodes=nodes,
        weights=np.full(n_phi, 4.0 * np.pi / n_phi),
        support=SPHERE,
    )


@pytest.mark.foundation
def test_quotient_agrees_with_the_geometric_section_orbit_by_orbit() -> None:
    r"""Two independent representative sections realize the SAME quotient.

    The verb picks each orbit's first-appearing member; the hand-rolled
    reference above picks the :math:`\xi \geq 0` member via
    :math:`\xi \to |\xi|`. Different atoms as embedded point sets — the
    same measure on :math:`S^2/\langle\sigma_y\rangle`: applying the
    geometric section to the verb's atoms must land each one on a
    reference atom carrying the SAME weight.
    """
    from orpheus.numerics.quadrature import Quadrature

    mu = Quadrature.product(n_mu=4, n_phi=8).measure

    folded = mu.quotient(SubgroupOfO3.Mirror("y"))

    reference_nodes = mu.nodes.copy()
    reference_nodes[:, 1] = np.abs(reference_nodes[:, 1])
    reference = mu.pushforward(
        ManifoldMap(
            SPHERE,
            SPHERE.quotient(SubgroupOfO3.Mirror("y")),
            lambda nodes: reference_nodes,
        )
    ).consolidate()

    assert folded.n_points == reference.n_points == 20  # Burnside (32+8)/2
    section_of_folded = folded.nodes.copy()
    section_of_folded[:, 1] = np.abs(section_of_folded[:, 1])
    for k in range(folded.n_points):
        (match,) = np.flatnonzero(
            np.linalg.norm(reference.nodes - section_of_folded[k], axis=1) < 1e-12
        )
        assert float(folded.weights[k]) == pytest.approx(
            float(reference.weights[match]), rel=1e-14
        )


@pytest.mark.foundation
def test_quotient_integrates_invariant_functions_like_the_parent() -> None:
    r"""The defining law: :math:`\int f \, d(\mu/G) = \int f \, d\mu` for
    :math:`G`-invariant :math:`f` — and ONLY for those.

    The odd arm is the knob-dependent positive control: a
    :math:`\xi`-odd integrand vanishes on the full sphere by
    cancellation and is decidedly nonzero on the fold (the quotient
    reporting its smaller space, plan §"the in-tree discriminator").
    If ``quotient`` degenerated to the identity, the odd arm — not the
    even one — reddens.
    """
    from orpheus.numerics.quadrature import Quadrature

    mu = Quadrature.product(n_mu=4, n_phi=8).measure
    folded = mu.quotient(SubgroupOfO3.Mirror("y"))

    def f_even(x: np.ndarray) -> np.ndarray:  # cosh is even in y
        return np.exp(x[:, 0]) * np.cosh(x[:, 1]) * (1.0 + x[:, 2] ** 2)

    def f_odd(x: np.ndarray) -> np.ndarray:  # y^3 is odd in y
        return x[:, 1] ** 3 * (1.0 + x[:, 0] ** 2)

    assert float(folded.integrate(f_even)) == pytest.approx(
        float(mu.integrate(f_even)), rel=1e-13
    )
    assert abs(float(mu.integrate(f_odd))) < 1e-13
    assert abs(float(folded.integrate(f_odd))) > 1e-2


@pytest.mark.foundation
def test_quotient_weights_are_orbit_stabilizer_derived() -> None:
    r"""Every folded weight is :math:`W = w \cdot |G|/|\mathrm{Stab}|`,
    read off the CERTIFICATE — and mass is preserved (quotient, not
    restriction)."""
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.numerics.symmetry import orbit_certificate

    mu = Quadrature.product(n_mu=4, n_phi=8).measure
    group = SubgroupOfO3.Mirror("y")

    cert = orbit_certificate(mu, group)
    assert cert is not None
    order = len(cert.operators)  # |G| = 2 for a mirror
    assert order == 2

    folded = mu.quotient(group)
    assert folded.n_points == len(cert.orbits())

    stab = cert.stabilizer_order
    for k in range(folded.n_points):
        (rep,) = np.flatnonzero(
            np.linalg.norm(mu.nodes - folded.nodes[k], axis=1) < 1e-14
        )
        predicted = float(mu.weights[rep]) * order / float(stab[rep])
        assert float(folded.weights[k]) == pytest.approx(predicted, rel=1e-14)

    np.testing.assert_array_almost_equal_nulp(
        np.array(folded.weights.sum()), np.array(mu.weights.sum()), nulp=4
    )


@pytest.mark.foundation
def test_a_free_action_folds_to_uniform_double_weights() -> None:
    r"""T24, the free-action certificate: :math:`\Sigma = \varnothing`
    :math:`\Longrightarrow` every orbit has length :math:`|G|`, and the
    orbit-stabilizer weight collapses to a uniform :math:`|G| \cdot w` —
    per-atom BIT-exact (``w + w`` never rounds), with no mixed
    :math:`w/2w` sum."""
    from orpheus.numerics.symmetry import orbit_certificate, singular_set

    mu = _fold_ring(n_phi=8)
    group = SubgroupOfO3.Mirror("y")

    assert singular_set(mu, group).size == 0  # Σ = ∅ — the action is free
    cert = orbit_certificate(mu, group)
    assert cert is not None
    assert all(orbit.size == len(cert.operators) for orbit in cert.orbits())

    folded = mu.quotient(group)
    assert folded.n_points == mu.n_points // 2
    np.testing.assert_array_equal(
        folded.weights, np.full(folded.n_points, 2.0 * (4.0 * np.pi / 8.0))
    )
    np.testing.assert_array_almost_equal_nulp(
        np.array(folded.weights.sum()), np.array(mu.weights.sum()), nulp=1
    )


@pytest.mark.foundation
def test_quotient_refuses_without_a_certificate() -> None:
    """No certificate — no quotient, loudly: a non-invariant measure and a
    continuous group both refuse (the precondition is the invariance
    PROOF, enforced by construction). Positive control: the certified
    fixture folds without raising."""
    group = SubgroupOfO3.Mirror("y")

    lopsided = DiscreteMeasure(
        nodes=np.array([[0.6, 0.8, 0.0], [0.6, 0.5, 0.62449979983984]]),
        weights=np.array([1.0, 1.0]),
        support=SPHERE,
    )
    with pytest.raises(ValueError, match="is defined only for a sigma_y-invariant"):
        lopsided.quotient(group)

    with pytest.raises(ValueError, match="is defined only for a SO3-invariant"):
        _fold_ring().quotient(SubgroupOfO3.SO3)

    assert _fold_ring().quotient(group).n_points == 4  # the honest pair


@pytest.mark.foundation
def test_the_fold_consumes_the_symmetry_idempotent_only_on_a_trivial_action() -> None:
    r"""The quotient CHANGES the space (T5), so re-folding is ill-posed.

    Free action: the folded measure keeps one member per mirror pair, so
    it is no longer :math:`\sigma_y`-invariant — a second ``quotient``
    REFUSES rather than silently halving again. Trivial action (every node
    fixed — the great circle in the :math:`xz`-plane, which :math:`\sigma_y`
    fixes pointwise): the quotient is the identity on ``(nodes, weights)`` and
    literally idempotent.

    ⭐ **The fixture changed at tracker 2.0c, and the reason is the campaign's
    own subject.** It used to be :math:`(\mu, 0, 0)` on ``"[-1,1]^slab"`` —
    described as "the slab embedding". Those points lie on **no** manifold
    carrying a :math:`\sigma_y` action: :math:`(-0.7, 0, 0)` is not on
    :math:`S^2`, and an interval has no :math:`O(3)` action at all. The string
    support could not notice — it fabricated the tag ``"[-1,1]^slab/sigma_y"``
    and the test went green. Typing ``support`` made
    :meth:`Manifold.quotient` refuse it, which is correct and is ERR-080's
    defect class exactly: a point set named by a CHART CODOMAIN rather than by
    the space it is. The honest :math:`\sigma_y`-fixed subset of the sphere is
    the :math:`xz` great circle, and it makes the claim stronger — the nodes
    are genuinely on the manifold whose quotient is being taken."""
    group = SubgroupOfO3.Mirror("y")

    folded = _fold_ring().quotient(group)
    assert folded.invariance_group is None  # the section is not equivariant
    with pytest.raises(ValueError, match="is defined only for a sigma_y-invariant"):
        folded.quotient(group)

    theta = np.array([0.3, 1.7, 2.9])
    fixed = DiscreteMeasure(
        nodes=np.column_stack([np.cos(theta), np.zeros(3), np.sin(theta)]),
        weights=np.array([0.4, 1.1, 0.5]),
        support=SPHERE,
    )
    # Every node is its own σ_y image, so the fold has nothing to merge.
    assert np.allclose(np.linalg.norm(fixed.nodes, axis=1), 1.0)
    once = fixed.quotient(group)
    np.testing.assert_array_equal(once.nodes, fixed.nodes)
    np.testing.assert_array_equal(once.weights, fixed.weights)

    # ⭐ ...and the fold is idempotent on the ATOMS ONLY. Re-folding is refused
    # even here, because the SPACE moved: `once` lives on S²/σ_y, and σ_y has
    # already been spent there. Until 2026-09-01 this assertion read
    # `twice.nodes == once.nodes` and passed — the string support happily
    # fabricated "…/sigma_y/sigma_y", so the test contradicted its own
    # docstring's "the quotient CHANGES the space, so re-folding is ill-posed"
    # and nothing could tell. Typing `support` made the two agree.
    assert once.support == SPHERE.quotient(group)
    with pytest.raises(NotImplementedError, match="S\\^2/sigma_y/sigma_y"):
        once.quotient(group)


@pytest.mark.foundation
def test_quotient_tags_the_orbifold_support_and_drops_both_claims() -> None:
    """``support`` becomes the orbifold tag; ``invariance_group`` and
    ``exactness`` are dropped (the fold consumed the symmetry; a claim
    would be against the pushforward reference, a type with no consumer
    yet — the direct-sum precedent)."""
    from orpheus.numerics.quadrature import Quadrature

    mu = Quadrature.product(n_mu=4, n_phi=8).measure
    assert mu.invariance_group is not None  # the knob the drop is about

    folded = mu.quotient(SubgroupOfO3.Mirror("y"))
    assert folded.support == SPHERE.quotient(SubgroupOfO3.Mirror("y"))
    assert folded.invariance_group is None
    assert folded.exactness is None
    assert folded.degree_of_exactness is None


# ---------------------------------------------------------------------------
# The exactness claim is relative to a measure
# ---------------------------------------------------------------------------
#
# `degree_of_exactness` states a degree; `generating_measure` states what
# integral that degree is about. Combining rules that reference DIFFERENT
# measures produces a rule exact against neither, so it may claim neither's
# degree. Before 2026-08-02 it claimed both.


@pytest.mark.l1
def test_mixed_reference_product_claims_the_PRODUCT_measure() -> None:
    r"""``gauss_legendre(4) * gauss_chebyshev(4)`` **is** exact to degree 7
    — against :math:`dx \otimes (1-y^2)^{-1/2}dy`, and this says so.

    ⚠ This test was inverted on 2026-08-02, and the correction is the
    point. It previously asserted the product must advertise **no**
    degree, on the strength of "it integrates the constant 1 to 6.2832
    where the answer is 4". `[M]` Re-measured with the reference named:
    :math:`2\pi = 6.2832` **is** the correct mass of
    ``legendre ⊗ chebyshev_t``; the expected ``4`` was the *Lebesgue*
    product, which is not this product's reference. Against its actual
    reference the rule is exact to degree 7 per axis to **4.16e-13**,
    with degree 8 missing by 1.5e-2 — so the claim is true AND tight.

    The old refusal was a conservative workaround for a missing type
    (there was no way to name :math:`\lambda_1 \otimes \lambda_2`), not
    a mathematical law. Naming the reference makes the correct claim
    representable, which is what the exactness carve is for.
    """
    from scipy.integrate import quad

    from orpheus.numerics.exactness import ProductMeasure

    gl, gc = gauss_legendre(4), gauss_chebyshev(4)
    assert gl.degree_of_exactness == 7
    assert gc.degree_of_exactness == 7
    assert gl.generating_measure != gc.generating_measure

    product = gl * gc
    claim = product.exactness
    assert claim is not None
    assert isinstance(claim.reference, ProductMeasure)
    assert claim.reference.name == "legendre ⊗ chebyshev_t"
    assert claim.degree == 7

    # The claim, verified against the reference it actually names.
    for a in range(8):
        for b in range(8):
            approx = float(np.sum(
                product.weights * product.nodes[:, 0] ** a
                * product.nodes[:, 1] ** b
            ))
            exact = (
                quad(lambda x: x ** a, -1, 1)[0]
                * quad(lambda y: y ** b / np.sqrt(1 - y * y), -1, 1)[0]
            )
            assert abs(approx - exact) < 1e-11, f"not exact at ({a}, {b})"

    # Tightness: the bound must actually bind one degree higher.
    approx8 = float(np.sum(
        product.weights * product.nodes[:, 0] ** 8 * product.nodes[:, 1] ** 8
    ))
    exact8 = (
        quad(lambda x: x ** 8, -1, 1)[0]
        * quad(lambda y: y ** 8 / np.sqrt(1 - y * y), -1, 1)[0]
    )
    assert abs(approx8 - exact8) > 1e-6, "the degree bound is not tight"

    # The DIRECT SUM is the case that genuinely has no claim: it lands on
    # the shared space, so its reference would be λ₁ + λ₂.
    assert (gl + gc).exactness is None


@pytest.mark.foundation
def test_product_keeps_a_claim_where_the_direct_sum_does_not() -> None:
    r"""The two composites diverge, and the reason is the reference.

    * **Product** — lands on the square, reference
      :math:`\lambda \otimes \lambda`, degree :math:`\min(5, 7) = 5`.
    * **Direct sum** — lands on the SAME interval, so its reference is
      :math:`\lambda + \lambda = 2\lambda`, a measure neither operand is
      exact against. Its total weight is ``4``, not ``2``.

    The predecessor asserted degree ``5`` for both. That was the
    half-claim: it kept a degree while dropping the reference, so nothing
    could notice the sum is a rule for twice the measure.
    """
    p = gauss_legendre(3) * gauss_legendre(4)
    assert p.exactness is not None
    assert p.exactness.degree == 5
    assert p.exactness.reference.name == "legendre ⊗ legendre"

    s = gauss_legendre(3) + gauss_legendre(4)
    assert s.exactness is None
    assert s.degree_of_exactness is None
    # The measurement that shows why: the sum integrates 1 to 4, not 2.
    assert float(np.sum(s.weights)) == pytest.approx(4.0)


@pytest.mark.foundation
def test_a_rule_with_no_reference_has_NO_claim_at_all() -> None:
    """The successor of ``test_untagged_rules_count_as_agreeing``.

    That test pinned the pre-2026-08-02 state: a measure could carry
    ``degree_of_exactness=3`` with no reference at all, and two such
    "equally unspecified" rules composed to a degree. **That state is now
    unrepresentable** — a degree lives inside an
    :class:`~orpheus.numerics.exactness.ExactnessClaim`, which cannot be
    built without naming what the claim is about.

    So the contract it documented is retired rather than migrated, and
    this is the replacement: no reference means no claim, and composing
    with a claimless rule yields no claim. The old behaviour ("keep the
    smaller degree") was the half-claim the carve exists to remove.
    """
    a = DiscreteMeasure(
        nodes=np.array([-0.5, 0.5]), weights=np.array([1.0, 1.0]),
        support=COSINE_INTERVAL,
    )
    b = DiscreteMeasure(
        nodes=np.array([-0.25, 0.25]), weights=np.array([1.0, 1.0]),
        support=COSINE_INTERVAL,
    )
    assert a.exactness is None and b.exactness is None
    assert a.degree_of_exactness is None
    assert a.generating_measure is None and b.generating_measure is None
    assert (a * b).exactness is None
    assert (a + b).exactness is None


@pytest.mark.foundation
def test_composites_name_no_generating_measure() -> None:
    """A composite is not a Golub-Welsch product, so ``generating_measure``
    is ``None`` — but for two different reasons, and the difference is the
    point.

    * a **product**'s reference IS named — ``legendre ⊗ legendre`` — it
      is simply a ``ProductMeasure``, which no three-term recurrence
      generates, so the narrower view reports ``None``;
    * a **direct sum** carries no claim at all (its reference would be
      :math:`\\lambda_1 + \\lambda_2`).
    """
    from orpheus.numerics.exactness import ProductMeasure

    p = gauss_legendre(3) * gauss_legendre(4)
    assert p.generating_measure is None
    assert p.exactness is not None
    assert isinstance(p.exactness.reference, ProductMeasure)
    assert p.exactness.reference.name == "legendre ⊗ legendre"

    s = gauss_legendre(3) + gauss_legendre(4)
    assert s.generating_measure is None
    assert s.exactness is None


@pytest.mark.foundation
def test_gauss_rules_carry_the_measure_that_built_them() -> None:
    """The constructors delegate to the Golub-Welsch body, so their
    output records which measure it is exact against."""
    from orpheus.numerics.generating_measure import CHEBYSHEV_T, LEGENDRE

    assert gauss_legendre(5).generating_measure == LEGENDRE
    assert gauss_chebyshev(5).generating_measure == CHEBYSHEV_T


@pytest.mark.foundation
def test_consolidate_preserves_the_generating_measure() -> None:
    """Consolidation merges coincident atoms and changes no integral,
    so the reference measure survives — same reason the degree does."""
    from orpheus.numerics.generating_measure import LEGENDRE

    rule = gauss_legendre(4)
    doubled = rule + rule  # every atom coincides with its twin
    merged = doubled.consolidate()
    assert merged.n_points == rule.n_points
    # The direct sum dropped the tag (its reference is the union), but
    # a rule that HAS one keeps it through consolidate.
    tagged = rule.consolidate()
    assert tagged.generating_measure == LEGENDRE
