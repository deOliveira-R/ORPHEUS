r"""Verification for ``orpheus.numerics.generating_measure``.

The module replaces several hand-rolled Gauss rules with ONE Golub-Welsch
body parameterised by a generating measure. That consolidation is only
safe if the generic body reproduces every family it absorbed, so the
tests come in three layers:

**L1 — against structurally independent references.** Each family is
compared to an external implementation (numpy's ``leggauss`` /
``laggauss`` / ``hermgauss``, scipy's ``roots_jacobi`` /
``roots_chebyu``) or to a closed form, and each rule is checked to
integrate polynomials exactly to :math:`2n-1` against *its own* measure
with the exact moments computed analytically. The moment tests are the
stronger evidence: they check the mathematical property the rule exists
for, not agreement with someone else's code.

**L1 — against itself, via a second construction.** ``jacobi(0,0)`` is
Legendre, ``jacobi(±1/2, ±1/2)`` are the two Chebyshevs. The general
recurrence and the specialised ones are different code reaching the same
measure, so their agreement is a genuine cross-check — the campaign's
"≥2 realizations *prove* an implementation" rule, with the second
realization free.

**foundation — the defining laws.** A Gauss rule is not just an array
pair: its weights are positive, its nodes are distinct, interior to the
support, and its weights sum to the measure's mass. These hold for every
family by theorem, so they are tested generically rather than per-family.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.special import gammaln, roots_chebyu, roots_jacobi

from orpheus.numerics.manifold import UNIT_INTERVAL
from orpheus.numerics.generating_measure import (
    CHEBYSHEV_T,
    CHEBYSHEV_U,
    HERMITE,
    LEGENDRE,
    GeneratingMeasure,
    jacobi,
    laguerre,
)

# Every shipped family, as (label, measure). Parameterised families are
# sampled at a non-degenerate point so the general code path runs.
ALL_FAMILIES: list[tuple[str, GeneratingMeasure]] = [
    ("legendre", LEGENDRE),
    ("chebyshev_t", CHEBYSHEV_T),
    ("chebyshev_u", CHEBYSHEV_U),
    ("hermite", HERMITE),
    ("laguerre(0)", laguerre()),
    ("laguerre(1.5)", laguerre(1.5)),
    ("jacobi(1.5, 2.25)", jacobi(1.5, 2.25)),
    ("jacobi(-0.5, 0.75)", jacobi(-0.5, 0.75)),
]


def _sorted_rule(measure: GeneratingMeasure, n: int):
    mu = measure.gauss(n)
    order = np.argsort(mu.nodes)
    return mu.nodes[order], mu.weights[order]


def _sorted_pair(nodes, weights):
    order = np.argsort(nodes)
    return np.asarray(nodes)[order], np.asarray(weights)[order]


# ===========================================================================
# L1 — each family against an independent implementation
# ===========================================================================


@pytest.mark.l1
@pytest.mark.parametrize("n", [1, 2, 5, 8, 16, 32])
def test_legendre_matches_numpy_leggauss(n: int) -> None:
    """numpy runs the same eigenproblem then Newton-refines, so it is an
    independent implementation of the same mathematics."""
    got = _sorted_rule(LEGENDRE, n)
    ref = _sorted_pair(*np.polynomial.legendre.leggauss(n))
    np.testing.assert_allclose(got[0], ref[0], atol=1e-14)
    np.testing.assert_allclose(got[1], ref[1], atol=1e-14)


@pytest.mark.l1
@pytest.mark.parametrize("n", [1, 2, 5, 8, 16])
def test_chebyshev_t_matches_the_closed_form(n: int) -> None:
    r"""Chebyshev-1 has a closed form — nodes at
    :math:`\cos((2i-1)\pi/2n)`, all weights :math:`\pi/n` — so this is a
    comparison against analysis, not against another implementation."""
    i = np.arange(1, n + 1)
    ref = _sorted_pair(
        np.cos((2.0 * i - 1.0) * np.pi / (2.0 * n)), np.full(n, np.pi / n)
    )
    got = _sorted_rule(CHEBYSHEV_T, n)
    np.testing.assert_allclose(got[0], ref[0], atol=1e-14)
    np.testing.assert_allclose(got[1], ref[1], atol=1e-14)


@pytest.mark.l1
@pytest.mark.parametrize("n", [1, 2, 5, 8, 16])
def test_chebyshev_u_matches_scipy(n: int) -> None:
    got = _sorted_rule(CHEBYSHEV_U, n)
    ref = _sorted_pair(*roots_chebyu(n))
    np.testing.assert_allclose(got[0], ref[0], atol=1e-14)
    np.testing.assert_allclose(got[1], ref[1], atol=1e-14)


@pytest.mark.l1
@pytest.mark.parametrize("n", [2, 5, 8, 16])
def test_hermite_matches_numpy_hermgauss(n: int) -> None:
    got = _sorted_rule(HERMITE, n)
    ref = _sorted_pair(*np.polynomial.hermite.hermgauss(n))
    np.testing.assert_allclose(got[0], ref[0], atol=1e-13)
    np.testing.assert_allclose(got[1], ref[1], atol=1e-14)


@pytest.mark.l1
@pytest.mark.parametrize("n", [2, 5, 8])
def test_laguerre_matches_numpy_laggauss(n: int) -> None:
    got = _sorted_rule(laguerre(0.0), n)
    ref = _sorted_pair(*np.polynomial.laguerre.laggauss(n))
    np.testing.assert_allclose(got[0], ref[0], atol=1e-13)
    np.testing.assert_allclose(got[1], ref[1], atol=1e-13)


@pytest.mark.l1
@pytest.mark.parametrize("n", [2, 5, 8])
@pytest.mark.parametrize(
    ("a", "b"),
    [(0.0, 0.0), (0.5, -0.5), (1.5, 2.25), (-0.5, -0.5), (3.0, 0.25)],
)
def test_jacobi_matches_scipy(n: int, a: float, b: float) -> None:
    got = _sorted_rule(jacobi(a, b), n)
    ref = _sorted_pair(*roots_jacobi(n, a, b))
    np.testing.assert_allclose(got[0], ref[0], atol=1e-13)
    np.testing.assert_allclose(got[1], ref[1], atol=1e-12)


# ===========================================================================
# L1 — exactness to 2n-1 against the family's OWN measure
# ===========================================================================
#
# This is the property the construction exists for, checked against
# analytically-known moments rather than another implementation.


def _legendre_moment(d: int) -> float:
    r"""\int_{-1}^{1} x^d dx."""
    return 0.0 if d % 2 else 2.0 / (d + 1.0)


def _hermite_moment(d: int) -> float:
    r"""\int_{-\infty}^{\infty} x^d e^{-x^2} dx = \Gamma((d+1)/2) for even d."""
    return 0.0 if d % 2 else float(np.exp(gammaln((d + 1) / 2.0)))


def _laguerre_moment(d: int) -> float:
    r"""\int_0^\infty x^d e^{-x} dx = d!."""
    return float(np.exp(gammaln(d + 1.0)))


def _chebyshev_t_moment(d: int) -> float:
    r"""\int_{-1}^{1} x^d (1-x^2)^{-1/2} dx.

    Substituting :math:`x = \cos\theta` gives
    :math:`\int_0^\pi \cos^d\theta\, d\theta`, which is
    :math:`\pi \binom{d}{d/2} / 2^d` for even :math:`d` and 0 for odd.
    """
    if d % 2:
        return 0.0
    from math import comb
    return np.pi * comb(d, d // 2) / 2.0**d


@pytest.mark.l1
@pytest.mark.parametrize(
    ("label", "measure", "moment"),
    [
        ("legendre", LEGENDRE, _legendre_moment),
        ("chebyshev_t", CHEBYSHEV_T, _chebyshev_t_moment),
        ("hermite", HERMITE, _hermite_moment),
        ("laguerre", laguerre(0.0), _laguerre_moment),
    ],
)
@pytest.mark.parametrize("n", [2, 4, 6])
def test_exact_to_2n_minus_1_against_its_own_measure(
    label: str, measure: GeneratingMeasure, moment, n: int
) -> None:
    r"""The rule reproduces :math:`\int x^d w\,dx` for every
    :math:`d \le 2n-1`."""
    nodes, weights = _sorted_rule(measure, n)
    for d in range(2 * n):
        got = float(np.sum(weights * nodes**d))
        want = moment(d)
        scale = max(abs(want), 1.0)
        assert abs(got - want) / scale < 1e-12, (
            f"{label} n={n}: degree {d} gave {got!r}, want {want!r}"
        )


@pytest.mark.l1
@pytest.mark.parametrize("n", [2, 3, 4, 5])
def test_not_exact_beyond_2n_minus_1(n: int) -> None:
    r"""Degree :math:`2n` must FAIL — otherwise the claim is not tight.

    Gated on a RELATIVE error and at small :math:`n` deliberately. At
    large :math:`n` the absolute error at degree :math:`2n` reads
    ~1e-16 simply because :math:`x^{2n}` on :math:`[-1,1]` collapses in
    magnitude, so an absolute-error gate here would silently stop
    biting (the T12c warning).
    """
    nodes, weights = _sorted_rule(LEGENDRE, n)
    d = 2 * n
    got = float(np.sum(weights * nodes**d))
    want = _legendre_moment(d)
    assert abs(got - want) / abs(want) > 1e-3, (
        f"n={n}: degree {d} was reproduced to "
        f"{abs(got - want) / abs(want):.2e} relative — the 2n-1 claim "
        f"is not tight, or the rule is not a Gauss rule"
    )


@pytest.mark.l1
def test_the_weighted_claim_is_not_an_unweighted_claim() -> None:
    r"""The distinction the old shared ``degree_of_exactness`` field erased.

    Gauss-Chebyshev is exact to :math:`2n-1` against
    :math:`\int q (1-x^2)^{-1/2} dx` and is **wrong** against
    :math:`\int q\,dx`. Both rules used to ship ``2n-1`` in the same
    field with nothing recording which integral was meant.
    """
    n, d = 4, 6
    gc_nodes, gc_weights = _sorted_rule(CHEBYSHEV_T, n)
    gl_nodes, gl_weights = _sorted_rule(LEGENDRE, n)

    # Exact against its OWN measure...
    weighted = float(np.sum(gc_weights * gc_nodes**d))
    assert abs(weighted - _chebyshev_t_moment(d)) < 1e-13

    # ...and decisively wrong against the unweighted one.
    unweighted_error = abs(weighted - _legendre_moment(d))
    assert unweighted_error > 0.5, (
        f"expected Gauss-Chebyshev to miss the unweighted integral badly; "
        f"error was only {unweighted_error:.3e}"
    )

    # Legendre gets the unweighted integral, which is the contrast.
    assert abs(
        float(np.sum(gl_weights * gl_nodes**d)) - _legendre_moment(d)
    ) < 1e-14

    # And the constructed measures SAY which is which, which is the
    # whole point of carrying the generating measure.
    assert CHEBYSHEV_T.gauss(n).generating_measure == CHEBYSHEV_T
    assert LEGENDRE.gauss(n).generating_measure == LEGENDRE
    assert CHEBYSHEV_T.gauss(n).degree_of_exactness == 2 * n - 1
    assert LEGENDRE.gauss(n).degree_of_exactness == 2 * n - 1


# ===========================================================================
# L1 — the internal cross-check: two constructions, one measure
# ===========================================================================


@pytest.mark.l1
@pytest.mark.parametrize("n", [2, 5, 8, 16])
@pytest.mark.parametrize(
    ("a", "b", "known"),
    [
        (0.0, 0.0, LEGENDRE),
        (-0.5, -0.5, CHEBYSHEV_T),
        (0.5, 0.5, CHEBYSHEV_U),
    ],
)
def test_jacobi_specialises_to_the_named_families(
    n: int, a: float, b: float, known: GeneratingMeasure
) -> None:
    """The general recurrence and the specialised one must agree.

    They are different code — the specialised constants use closed-form
    coefficients, the Jacobi path evaluates the general formula — so
    this is a real cross-check and not an alias. `[M]` The nodes agree
    bit-identically.
    """
    general = jacobi(a, b)
    assert general.recurrence is not known.recurrence, (
        "the two must be different code paths for this to verify anything"
    )
    gen_nodes, gen_weights = _sorted_rule(general, n)
    kn_nodes, kn_weights = _sorted_rule(known, n)
    np.testing.assert_array_equal(gen_nodes, kn_nodes)
    np.testing.assert_allclose(gen_weights, kn_weights, atol=1e-15)


@pytest.mark.foundation
@pytest.mark.parametrize(
    ("a", "b", "known"),
    [
        (0.0, 0.0, LEGENDRE),
        (-0.5, -0.5, CHEBYSHEV_T),
        (0.5, 0.5, CHEBYSHEV_U),
    ],
)
def test_specialisations_are_the_same_MEASURE(
    a: float, b: float, known: GeneratingMeasure
) -> None:
    """Equality is mathematical identity, not construction path.

    ``jacobi(0,0)`` *is* Legendre. Comparing unequal because it was
    reached through a different constructor would be a false statement
    about the mathematics.
    """
    assert jacobi(a, b) == known
    assert jacobi(a, b).name == known.name
    assert jacobi(a, b).support == known.support


@pytest.mark.foundation
def test_chebyshev_t_construction_does_not_divide_by_zero() -> None:
    r"""Regression: the general Jacobi :math:`\beta_k` formula has a
    REMOVABLE singularity at :math:`k=1` when :math:`a+b=-1`.

    Its numerator carries :math:`(k+a+b)` and its denominator
    :math:`(2k+a+b-1)`; at :math:`k=1` both equal :math:`1+a+b`, so
    they cancel — but only if the code cancels them. Evaluating the
    general form there raises ``ZeroDivisionError`` for
    ``jacobi(-1/2, -1/2)``, which is Chebyshev-1, the single most common
    member of the family.
    """
    # The singular locus a + b = -1 intersected with the integrability
    # domain a, b > -1 forces BOTH exponents into (-1, 0) — so this is
    # the complete shape of the hazard, not a sample of it.
    for a, b in ((-0.5, -0.5), (-0.25, -0.75), (-0.9, -0.1), (-0.1, -0.9)):
        assert a + b == pytest.approx(-1.0), "fixture must sit on the locus"
        rule = jacobi(a, b).gauss(6)
        assert np.all(np.isfinite(rule.nodes))
        assert np.all(np.isfinite(rule.weights))
        np.testing.assert_allclose(
            float(np.sum(rule.weights)), jacobi(a, b).mass, rtol=1e-12
        )
    # And the headline case end-to-end.
    np.testing.assert_allclose(
        _sorted_rule(jacobi(-0.5, -0.5), 8)[0],
        _sorted_rule(CHEBYSHEV_T, 8)[0],
        atol=1e-15,
    )


# ===========================================================================
# foundation — the defining laws of a Gauss rule
# ===========================================================================


@pytest.mark.foundation
@pytest.mark.parametrize(("label", "measure"), ALL_FAMILIES)
@pytest.mark.parametrize("n", [1, 2, 5, 9])
def test_weights_sum_to_the_measures_mass(
    label: str, measure: GeneratingMeasure, n: int
) -> None:
    r""":math:`\sum_i w_i = \mu_0 = \int w\,dx`, for every :math:`n`.

    Degree 0 of the exactness claim, and the one invariant that ties
    the discrete rule back to the continuous measure that made it.
    """
    rule = measure.gauss(n)
    np.testing.assert_allclose(
        float(np.sum(rule.weights)), measure.mass, rtol=1e-13
    )


@pytest.mark.foundation
@pytest.mark.parametrize(("label", "measure"), ALL_FAMILIES)
@pytest.mark.parametrize("n", [2, 5, 9])
def test_gauss_weights_are_positive(
    label: str, measure: GeneratingMeasure, n: int
) -> None:
    """Positivity is a theorem for Gauss rules, not a coincidence.

    It follows from :math:`w_i = \\mu_0 [v_i]_1^2` being a squared
    quantity times a positive mass, so a negative weight would mean the
    construction had gone wrong somewhere structural.
    """
    rule = measure.gauss(n)
    assert np.all(rule.weights > 0.0), (
        f"{label} n={n}: negative weights at "
        f"{np.flatnonzero(rule.weights <= 0.0)}"
    )


@pytest.mark.foundation
@pytest.mark.parametrize(("label", "measure"), ALL_FAMILIES)
@pytest.mark.parametrize("n", [2, 5, 9])
def test_nodes_are_distinct_and_ascending(
    label: str, measure: GeneratingMeasure, n: int
) -> None:
    """The nodes are eigenvalues of a symmetric tridiagonal matrix with
    non-zero off-diagonal, hence real and SIMPLE; ``eigh`` returns them
    ascending, and the rule relies on that ordering."""
    rule = measure.gauss(n)
    assert np.all(np.diff(rule.nodes) > 0.0), (
        f"{label} n={n}: nodes not strictly ascending: {rule.nodes}"
    )


@pytest.mark.foundation
@pytest.mark.parametrize("n", [2, 5, 9, 20])
def test_bounded_family_nodes_lie_strictly_inside_the_interval(
    n: int,
) -> None:
    """Gauss nodes are interior — the endpoints are never nodes.

    (That is precisely what a Lobatto constraint would change, which is
    why node-constraint is a separate axis from the family.)
    """
    for measure in (LEGENDRE, CHEBYSHEV_T, CHEBYSHEV_U, jacobi(1.5, 2.25)):
        rule = measure.gauss(n)
        assert np.all(rule.nodes > -1.0) and np.all(rule.nodes < 1.0), (
            f"{measure.name} n={n}: node outside (-1, 1)"
        )


@pytest.mark.foundation
@pytest.mark.parametrize("n", [2, 5, 9])
def test_half_line_family_nodes_are_positive(n: int) -> None:
    """Laguerre lives on [0, inf), so every node must be > 0."""
    for a in (0.0, 1.5):
        assert np.all(laguerre(a).gauss(n).nodes > 0.0)


@pytest.mark.foundation
@pytest.mark.parametrize(("label", "measure"), ALL_FAMILIES)
def test_mass_is_read_from_the_recurrence(
    label: str, measure: GeneratingMeasure
) -> None:
    r""":attr:`mass` is :math:`\beta_0`, not a separate stored field.

    Storing the zeroth moment beside the recurrence would let the two
    disagree; reading it out means the mass is the mass of *this*
    family by construction.
    """
    _, beta = measure.recurrence(1)
    assert measure.mass == float(beta[0])
    # And it agrees with the rule it generates, at every n.
    for n in (1, 3, 7):
        np.testing.assert_allclose(
            float(np.sum(measure.gauss(n).weights)), measure.mass, rtol=1e-13
        )


@pytest.mark.foundation
@pytest.mark.parametrize(
    ("measure", "expected"),
    [
        (LEGENDRE, 2.0),
        (CHEBYSHEV_T, np.pi),
        (CHEBYSHEV_U, np.pi / 2.0),
        (HERMITE, np.sqrt(np.pi)),
        (laguerre(0.0), 1.0),
    ],
)
def test_masses_match_their_closed_forms(
    measure: GeneratingMeasure, expected: float
) -> None:
    np.testing.assert_allclose(measure.mass, expected, rtol=1e-14)


@pytest.mark.foundation
def test_generated_rule_carries_its_measure_and_degree() -> None:
    """The construction records what produced it — which is what makes
    the exactness claim self-describing."""
    for _, measure in ALL_FAMILIES:
        for n in (1, 4, 7):
            rule = measure.gauss(n)
            assert rule.generating_measure == measure
            assert rule.degree_of_exactness == 2 * n - 1
            assert rule.support == measure.support


@pytest.mark.foundation
def test_single_point_rule_is_the_mean_carrying_all_the_mass() -> None:
    r"""The :math:`n=1` branch bypasses the eigen-solve, so it needs its
    own check: :math:`x_1 = \alpha_0`, :math:`w_1 = \mu_0`."""
    for _, measure in ALL_FAMILIES:
        rule = measure.gauss(1)
        alpha, beta = measure.recurrence(1)
        assert rule.nodes.shape == (1,)
        assert rule.nodes[0] == float(alpha[0])
        assert rule.weights[0] == float(beta[0])
        # Degree 1 exactness means it integrates constants AND linears.
        assert rule.degree_of_exactness == 1


# ===========================================================================
# foundation — the affine remap morphism
# ===========================================================================


@pytest.mark.l1
@pytest.mark.parametrize(("a", "b"), [(0.0, 1.0), (-3.0, 2.0), (2.5, 7.25)])
@pytest.mark.parametrize("n", [2, 5, 8])
def test_affine_remap_matches_the_transformed_rule(
    a: float, b: float, n: int
) -> None:
    r"""``LEGENDRE.on(a, b).gauss(n)`` equals the analytically remapped
    rule: nodes :math:`\tfrac{1}{2}[(b-a)t + (a+b)]`, weights scaled by
    the Jacobian :math:`(b-a)/2`."""
    base = LEGENDRE.gauss(n)
    remapped = LEGENDRE.on(a, b).gauss(n)
    scale, shift = (b - a) / 2.0, (a + b) / 2.0
    np.testing.assert_allclose(
        np.sort(remapped.nodes), np.sort(scale * base.nodes + shift), atol=1e-13
    )
    np.testing.assert_allclose(
        np.sum(remapped.weights), scale * np.sum(base.weights), rtol=1e-14
    )


@pytest.mark.l1
@pytest.mark.parametrize(("a", "b"), [(0.0, 1.0), (-3.0, 2.0)])
def test_remapped_rule_integrates_on_the_new_interval(
    a: float, b: float
) -> None:
    r"""The remapped Legendre rule is exact for :math:`\int_a^b x^d dx`
    to degree :math:`2n-1` — the property the remap exists for."""
    n = 5
    rule = LEGENDRE.on(a, b).gauss(n)
    for d in range(2 * n):
        want = (b ** (d + 1) - a ** (d + 1)) / (d + 1)
        got = float(np.sum(rule.weights * rule.nodes**d))
        np.testing.assert_allclose(got, want, rtol=1e-11, atol=1e-12)


@pytest.mark.foundation
def test_remap_reports_the_new_interval_and_mass() -> None:
    remapped = LEGENDRE.on(0.0, 1.0)
    # ⭐ `Interval(0.0, 1.0)`, not the string `"[0.0,1.0]"` the f-string used to
    # build. The name normalises to "[0,1]" — a principled re-baseline, and the
    # reason to compare the OBJECT: the interval is the same one whatever the
    # float repr of its endpoints, which a string comparison could not say.
    assert remapped.support == UNIT_INTERVAL
    np.testing.assert_allclose(remapped.mass, 1.0, rtol=1e-15)
    assert remapped.gauss(3).support == UNIT_INTERVAL


# ===========================================================================
# foundation — guards
# ===========================================================================


@pytest.mark.foundation
def test_gauss_rejects_non_positive_n() -> None:
    for n in (0, -1):
        with pytest.raises(ValueError, match="n >= 1"):
            LEGENDRE.gauss(n)


@pytest.mark.foundation
@pytest.mark.parametrize(("a", "b"), [(-1.0, 0.0), (0.0, -1.0), (-2.0, -2.0)])
def test_jacobi_rejects_non_integrable_exponents(a: float, b: float) -> None:
    """The weight :math:`(1-x)^a(1+x)^b` is not integrable at
    :math:`a \\le -1` or :math:`b \\le -1`, so the measure does not
    exist and no rule can be built for it."""
    with pytest.raises(ValueError, match="a > -1"):
        jacobi(a, b)


@pytest.mark.foundation
def test_laguerre_rejects_non_integrable_exponent() -> None:
    with pytest.raises(ValueError, match="a > -1"):
        laguerre(-1.0)


@pytest.mark.foundation
def test_remap_refuses_unbounded_families() -> None:
    """There is no affine map from :math:`[-1,1]` onto an unbounded
    support, so the morphism is undefined there and must say so."""
    for measure in (HERMITE, laguerre()):
        with pytest.raises(ValueError, match="defined for measures on"):
            measure.on(0.0, 1.0)


@pytest.mark.foundation
def test_remap_requires_an_ordered_interval() -> None:
    with pytest.raises(ValueError, match="a < b"):
        LEGENDRE.on(1.0, 0.0)


# ===========================================================================
# foundation — symmetry is IMPOSED, and the condition for it is DERIVED
# ===========================================================================
#
# Learned by reading numpy's `leggauss` and scipy's
# `_gen_roots_and_weights`: both end with
#     x = (x - x[::-1]) / 2 ;  w = (w + w[::-1]) / 2
# for symmetric measures, and both rescale so sum(w) == mu_0. They impose
# the structural properties rather than inheriting them to within
# round-off. scipy passes `symmetrize` as a hand-set boolean; here the
# condition is a theorem about the recurrence, so it is read off the data.


@pytest.mark.foundation
@pytest.mark.parametrize(
    ("label", "measure", "expected"),
    [
        ("legendre", LEGENDRE, True),
        ("chebyshev_t", CHEBYSHEV_T, True),
        ("chebyshev_u", CHEBYSHEV_U, True),
        ("hermite", HERMITE, True),
        ("laguerre(0)", laguerre(), False),
        ("laguerre(1.5)", laguerre(1.5), False),
        ("jacobi(2,2) a==b", jacobi(2.0, 2.0), True),
        ("jacobi(0.5,0.5) a==b", jacobi(0.5, 0.5), True),
        ("jacobi(1.5,2.25) a!=b", jacobi(1.5, 2.25), False),
        ("jacobi(-0.5,0.75) a!=b", jacobi(-0.5, 0.75), False),
    ],
)
def test_is_symmetric_is_derived_correctly(
    label: str, measure: GeneratingMeasure, expected: bool
) -> None:
    r""":math:`\alpha_k \equiv 0 \iff w` is even — a theorem, not a flag.

    The parameterised cases are the ones that matter: ``jacobi(a, b)``
    is symmetric exactly when ``a == b``, and the derivation gets that
    right without being told, which a declared flag could not.
    """
    assert measure.is_symmetric is expected


@pytest.mark.foundation
@pytest.mark.parametrize(
    ("label", "measure"),
    [
        ("legendre", LEGENDRE),
        ("chebyshev_t", CHEBYSHEV_T),
        ("chebyshev_u", CHEBYSHEV_U),
        ("hermite", HERMITE),
        ("jacobi(2,2)", jacobi(2.0, 2.0)),
    ],
)
@pytest.mark.parametrize("n", [2, 7, 8, 16, 33])
def test_symmetric_measures_give_EXACTLY_symmetric_rules(
    label: str, measure: GeneratingMeasure, n: int
) -> None:
    r"""Bit-exact, not ``allclose``: :math:`x_i = -x_{n-1-i}` and
    :math:`w_i = w_{n-1-i}` to the last bit.

    This is the property an angular quadrature's invariance group rests
    on. A reflecting boundary is exactly representable only if the node
    set is closed under the face reflection; when the rule is exactly
    symmetric that closure is integer index arithmetic rather than a
    tolerance question — which is the whole thrust of the Sigma /
    invariance-group machinery in
    :mod:`orpheus.numerics.symmetry`.
    """
    rule = measure.gauss(n)
    np.testing.assert_array_equal(rule.nodes, -rule.nodes[::-1])
    np.testing.assert_array_equal(rule.weights, rule.weights[::-1])


@pytest.mark.foundation
@pytest.mark.parametrize("n", [3, 7, 9, 33])
def test_odd_order_symmetric_rules_have_an_exact_zero_node(n: int) -> None:
    """The centre node of an odd symmetric rule is exactly ``0.0``.

    Not ``1e-17``: a downstream degeneracy test that asks ``mu == 0``
    (rather than ``abs(mu) < eps``) is then answerable, and the
    ordinate sits exactly on the mirror rather than near it.
    """
    for measure in (LEGENDRE, CHEBYSHEV_T, HERMITE, jacobi(2.0, 2.0)):
        rule = measure.gauss(n)
        assert rule.nodes[n // 2] == 0.0, (
            f"{measure.name} n={n}: centre node is "
            f"{rule.nodes[n // 2]!r}, not exactly zero"
        )


@pytest.mark.foundation
@pytest.mark.parametrize(("label", "measure"), ALL_FAMILIES)
@pytest.mark.parametrize("n", [2, 8, 32])
def test_mass_is_imposed_not_inherited(
    label: str, measure: GeneratingMeasure, n: int
) -> None:
    r""":math:`\sum w_i = \mu_0` to within one rounding of the true
    value, because the weights are rescaled to it.

    Degree-0 exactness is the one coefficient known in closed form, so
    there is no reason to let the eigensolver's round-off decide it.
    """
    rule = measure.gauss(n)
    total = float(np.sum(rule.weights))
    assert abs(total - measure.mass) <= 4.0 * np.spacing(measure.mass), (
        f"{label} n={n}: sum(w) = {total!r} vs mass {measure.mass!r}"
    )


@pytest.mark.foundation
def test_asymmetric_measures_are_not_symmetrised() -> None:
    """Laguerre lives on a half-line; mirroring it would be nonsense.

    The guard is the derived :attr:`is_symmetric`, so this checks the
    derivation actually gates the transform.
    """
    rule = laguerre().gauss(8)
    assert not np.allclose(rule.nodes, -rule.nodes[::-1])
    assert np.all(rule.nodes > 0.0)
