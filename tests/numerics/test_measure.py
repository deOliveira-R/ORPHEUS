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
guard the data structure contract that downstream consumers
(:class:`orpheus.sn.quadrature.AngularQuadrature` adapters in Issue 4,
MoC bundle measures in Wave 2) will rely on.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.measure import (
    BundleMeasure,
    DiscreteMeasure,
    SPACE_INTERVAL_M11,
    equispaced,
    gauss_chebyshev,
    gauss_legendre,
)


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
    assert p.space == "[-1,1] × [-1,1]"
    assert p.degree_of_exactness == 5
    assert p.invariance_group is None


@pytest.mark.foundation
def test_direct_sum_concatenates_on_shared_space() -> None:
    """``μ + ν`` requires equal ``space``; concatenates atoms."""
    a = equispaced(0.0, 0.5, 4)
    b = equispaced(0.0, 0.5, 6)
    s = a + b
    assert s.n_points == 10
    assert s.space == a.space
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
        lambda x: 0.5 * (x - 1.0), new_space=SPACE_INTERVAL_M11,
    )
    # The pushforward x -> (x-1)/2 maps [-1,1] to [-1, 0]. The Jacobian
    # (1/2) MUST be applied to the weights manually — the pushforward
    # is φ-image, not Radon-Nikodym (per docstring).
    base_left = DiscreteMeasure(
        nodes=base_left.nodes,
        weights=base_left.weights * 0.5,
        space=base_left.space,
    )
    base_right = gauss_legendre(4).pushforward(
        lambda x: 0.5 * (x + 1.0), new_space=SPACE_INTERVAL_M11,
    )
    base_right = DiscreteMeasure(
        nodes=base_right.nodes,
        weights=base_right.weights * 0.5,
        space=base_right.space,
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
    pushed = base.pushforward(lambda x: x ** 2, new_space="[0,1]")
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
    base = gauss_legendre(5).with_metadata(invariance_group="O(1)")
    assert base.invariance_group == "O(1)"
    assert base.degree_of_exactness == 9
    pushed = base.pushforward(lambda x: x ** 3, new_space="img")
    assert pushed.invariance_group is None
    assert pushed.degree_of_exactness is None


@pytest.mark.foundation
def test_pushforward_default_space_tag() -> None:
    """When no ``new_space`` is provided, the target tag is
    ``f"φ_*({source.space})"``."""
    base = gauss_legendre(3)
    pushed = base.pushforward(lambda x: x + 1.0)
    assert pushed.space == "φ_*([-1,1])"


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
    rule = gauss_legendre(8).with_metadata(invariance_group="O(1)")
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
            space="R",
        )
    with pytest.raises(ValueError, match="disagree on N"):
        DiscreteMeasure(
            nodes=np.array([1.0, 2.0, 3.0]),
            weights=np.array([1.0, 1.0]),
            space="R",
        )
    with pytest.raises(ValueError, match="nodes must be 1-D"):
        DiscreteMeasure(
            nodes=np.zeros((2, 2, 2)),
            weights=np.array([1.0, 1.0]),
            space="R",
        )


@pytest.mark.foundation
def test_with_metadata_returns_new_instance() -> None:
    """``with_metadata`` returns a new frozen instance with selected
    fields updated; the original is unchanged."""
    base = gauss_legendre(3)
    tagged = base.with_metadata(invariance_group="O(1)")
    assert tagged.invariance_group == "O(1)"
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
    rule = gauss_legendre(8).with_metadata(invariance_group="Z2")
    r = repr(rule)
    assert "invariance_group" in r
    assert "'Z2'" in r
