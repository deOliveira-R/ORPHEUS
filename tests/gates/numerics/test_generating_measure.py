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

from fractions import Fraction

import numpy as np
import pytest
from scipy.special import gammaln, roots_chebyu, roots_jacobi

from tests.gates.geometry._manifold_fixtures import UNIT_INTERVAL
from orpheus.numerics.generating_measure import (
    CHEBYSHEV_T,
    CHEBYSHEV_U,
    LEGENDRE,
    GeneratingMeasure,
    jacobi,
)
from tests.gates.numerics._measure_fixtures import HERMITE, laguerre
from orpheus.geometry.manifold import CIRCLE, COSINE_INTERVAL, SPHERE
from orpheus.geometry.symmetry import SubgroupOfO3
from tests.gates.geometry._lattice_fixtures import _g22b_AXES, _ROSTER
import orpheus.numerics.recurrence as r

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
    assert CHEBYSHEV_T.gauss(n).exactness.reference == CHEBYSHEV_T
    assert LEGENDRE.gauss(n).exactness.reference == LEGENDRE
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
    assert general.exact_recurrence is not known.exact_recurrence, (
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
            assert rule.exactness.reference == measure
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
    :mod:`orpheus.geometry.symmetry`.
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
def test_weight_sum_is_within_four_float_steps_of_the_mass(
    label: str, measure: GeneratingMeasure, n: int
) -> None:
    r""":math:`\sum w_i = \mu_0` to within a few roundings of the true
    value, because every weight is correctly rounded.

    Until the correctly rounded construction (2026-10-01) the weights were
    rescaled to the mass after an eigensolve. Nothing is imposed now; the
    bound holds because each weight is the float nearest its exact value,
    so the sum carries only the rounding of its terms and of the summation.
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


# ===========================================================================
# foundation — the correctly rounded construction refuses what it cannot certify
# ===========================================================================


@pytest.mark.foundation
def test_a_value_its_error_interval_cannot_certify_is_refused(monkeypatch) -> None:
    r"""The certification is live: a value whose error interval straddles a
    float64 rounding boundary raises instead of shipping a guess.

    The witness lowers the coarse precision to 16 digits, where the error
    radius :math:`10^{10-16}` (relative) dwarfs a float64 spacing, so no value
    of a 64-point Legendre rule can be certified. The positive control (the
    shipped precisions certify the same rule) runs first.
    """
    import orpheus.numerics.generating_measure as gm

    LEGENDRE.gauss(64)  # positive control: the shipped precisions certify
    monkeypatch.setattr(gm, "_WORKING_PRECISIONS", (r.HighPrecision(16), r.HighPrecision(60)))
    gm._correctly_rounded_rule.cache_clear()
    try:
        with pytest.raises(ArithmeticError, match="cannot be certified"):
            LEGENDRE.gauss(64)
    finally:
        gm._correctly_rounded_rule.cache_clear()


@pytest.mark.foundation
def test_a_rule_that_rounds_differently_at_the_two_precisions_is_refused(monkeypatch) -> None:
    r"""The two-precision comparison is live on its own: with the
    certification disarmed (plain rounding), a 16-digit coarse run rounds the
    64-point Legendre rule differently from the 60-digit run, and the
    comparison refuses it. Disarming isolates this refusal from the
    certification, which would otherwise fire first."""
    import orpheus.numerics.generating_measure as gm

    monkeypatch.setattr(gm, "_WORKING_PRECISIONS", (r.HighPrecision(16), r.HighPrecision(60)))
    monkeypatch.setattr(gm, "_certified_float", lambda value, precision, measure, n, **_: float(value))
    gm._correctly_rounded_rule.cache_clear()
    try:
        with pytest.raises(ArithmeticError, match="rounds differently"):
            LEGENDRE.gauss(64)
    finally:
        gm._correctly_rounded_rule.cache_clear()


@pytest.mark.foundation
def test_a_node_newton_cannot_converge_is_refused(monkeypatch) -> None:
    r"""Newton's refusal is live: with no step allowed, no seed converges,
    and the construction raises rather than rounding an unrefined seed."""
    import orpheus.numerics.generating_measure as gm

    monkeypatch.setattr(gm, "_NEWTON_STEP_LIMIT", 0)
    gm._correctly_rounded_rule.cache_clear()
    try:
        with pytest.raises(ArithmeticError, match="did not converge"):
            LEGENDRE.gauss(6)
    finally:
        gm._correctly_rounded_rule.cache_clear()


# ---------------------------------------------------------------------------
# The exact tie: an exact value on a float64 rounding tie is rounded half to even
# (unit `quadrature-ties`; the law is 2.0.3 of `.claude/plans/ontology_w0_spec.md`).
#
# Ladder.  `test_rational_polynomial.py` (the polynomial ring over Q) and the generic-builder
# rows (G6) are the foundations; the decision called directly (G4) rests on both; the
# consumer-level rows (G1, G5) rest on G4; G2 and G3 are the controls that an untied value
# is unchanged and an unequal single-boundary value is still refused.  Every private name
# is read THROUGH the module (``gm.<name>``), never imported by name, so a monkeypatch of
# ``gm`` reaches it.
# ---------------------------------------------------------------------------

_THIS = "tests/gates/numerics/test_generating_measure.py::"
_G6 = _THIS + "test_the_generic_builders_admit_the_tabulated_rationals"
_G6C = _THIS + "test_the_remap_maps_the_rational_recurrence_exactly"
_G4 = _THIS + "test_the_weight_decision_refuses_an_unequal_tie_and_has_its_domain"
_G7 = "tests/gates/numerics/test_rational_polynomial.py::test_gcd_is_the_monic_greatest_common_divisor"

_WIDTH_TIE = (0.2615032494292951, 0.7646335475579331)  # (b - a) on a tie, (a + b)/2 not: the even-n witness
_N3_OUTER = (0.29134508037323065, 1.218729331071271)  # 5 s/9 on a tie (the n = 3 outer weights)
_N3_CENTRE = (0.08349338193695595, 2.0129640274100478)  # 8 s/9 on a tie (the n = 3 centre weight)
_LOWER = (0.5118216247002567, 0.9504636963259353)  # the centre's half-even neighbour is the LOWER float
_MIRROR = (-_WIDTH_TIE[1], -_WIDTH_TIE[0])  # the panel of the only rows with negative nodes


def _centre(a: float, b: float) -> Fraction:
    return (Fraction(a) + Fraction(b)) / 2


def _half_width(a: float, b: float) -> Fraction:
    return (Fraction(b) - Fraction(a)) / 2


def _is_tie(value: Fraction) -> bool:
    """``value`` is exactly the midpoint of two adjacent float64 (so "nearest" names two floats)."""
    f = float(value)
    if Fraction(f) == value:
        return False
    other = float(np.nextafter(f, np.inf if Fraction(f) < value else -np.inf))
    return (Fraction(f) + Fraction(other)) / 2 == value


def _clear_caches() -> None:
    import orpheus.numerics.generating_measure as gm

    gm._affine_remap.cache_clear()
    gm._correctly_rounded_rule.cache_clear()


def _fresh(measure, panel):
    """``measure.on(*panel)`` built with both ``lru_cache`` s cleared: a measure compares equal to another by
    NAME, so a remap cached by an earlier row would otherwise be returned without running this row's builder."""
    _clear_caches()
    return measure.on(*panel)


# (panel, n, which array, index, exact closed form of the tied value)
_TIE_ROWS = [
    ((0.2, 0.4), 1, "nodes", 0, _centre(0.2, 0.4)),
    ((0.2, 0.4), 3, "nodes", 1, _centre(0.2, 0.4)),
    ((0.55, 0.63), 63, "nodes", 31, _centre(0.55, 0.63)),
    (_LOWER, 3, "nodes", 1, _centre(*_LOWER)),
    (_WIDTH_TIE, 1, "weights", 0, 2 * _half_width(*_WIDTH_TIE)),
    (_WIDTH_TIE, 2, "weights", 0, _half_width(*_WIDTH_TIE)),
    (_WIDTH_TIE, 2, "weights", 1, _half_width(*_WIDTH_TIE)),
    (_N3_OUTER, 3, "weights", 0, _half_width(*_N3_OUTER) * Fraction(5, 9)),
    (_N3_OUTER, 3, "weights", 2, _half_width(*_N3_OUTER) * Fraction(5, 9)),
    (_N3_CENTRE, 3, "weights", 1, _half_width(*_N3_CENTRE) * Fraction(8, 9)),
    (_MIRROR, 2, "weights", 0, _half_width(*_MIRROR)),
    (_MIRROR, 2, "weights", 1, _half_width(*_MIRROR)),
    (_MIRROR, 1, "weights", 0, 2 * _half_width(*_MIRROR)),
]


@pytest.mark.foundation
@pytest.mark.catches("ERR-127")
@pytest.mark.rests_on(_G4, _G6, _G6C, _G7)
@pytest.mark.parametrize("panel, n, field, index, exact", _TIE_ROWS, ids=lambda v: str(v))
def test_an_exact_value_on_a_rounding_tie_is_rounded_half_even(panel, n, field, index, exact) -> None:
    r"""An exact value of ``LEGENDRE.on(a, b).gauss(n)`` that sits on a float64 rounding tie is the float64
    neighbour that is even, never a refusal (the law of 2.0.3).

    The oracle is the CLOSED FORM of the tied value (``c = (a + b)/2``, ``s = (b - a)/2``, ``2 s``, ``s``,
    ``5 s/9``, ``8 s/9``) rounded by ``Fraction.__float__`` (correctly rounded, half to even), independent of the
    rational recurrence the decision reads.  Rows: centres at odd ``n`` (one whose half-even neighbour is the
    LOWER float, the round-up mutant's witness), width ties at ``n = 1`` and at EVEN ``n = 2`` (parity is not the
    law), the non-dyadic reference weights ``5/9`` and ``8/9`` at ``n = 3``, and the three MIRRORED rows, the only
    ones with negative nodes (a sign-dropping rational image of a node refuses them).  Each row first asserts that
    its closed form is a tie, so a row that stopped being one cannot pass vacuously."""
    assert _is_tie(exact), "fixture control: the closed form must sit exactly on a float64 tie"
    _clear_caches()
    try:
        rule = LEGENDRE.on(*panel).gauss(n)
        assert getattr(rule, field)[index] == float(exact)
    finally:
        _clear_caches()


@pytest.mark.foundation
def test_the_untied_sibling_is_unchanged() -> None:
    """Control (builds today and after): an untied odd panel's centre is the representable midpoint."""
    assert not _is_tie(_centre(0.0, 1.0))
    assert LEGENDRE.on(0.0, 1.0).gauss(3).nodes[1] == 0.5


@pytest.mark.foundation
def test_a_single_boundary_straddle_that_is_not_the_exact_value_is_still_refused(monkeypatch) -> None:
    r"""The anti-guesser, nodes: ``_CERTIFIED_DIGITS = 23`` widens the radius to :math:`10^{23-d}` (relative), so
    22 of the 128 values of ``LEGENDRE.on(0.0, 1.0).gauss(64)`` straddle exactly ONE boundary at 40 digits and
    none straddles several.  The first is a NODE, and none of them is a tie, so the decision must say "not equal"
    and the value stays refused.  The ``match`` is load-bearing: a decision that always says "equal" is refused
    later, by the two-precision comparison (``rounds differently``), so a bare ``pytest.raises(ArithmeticError)``
    stays green under it and only ``cannot be certified`` reddens."""
    import orpheus.numerics.generating_measure as gm

    monkeypatch.setattr(gm, "_CERTIFIED_DIGITS", 23)
    _clear_caches()
    try:
        with pytest.raises(ArithmeticError, match="cannot be certified"):
            LEGENDRE.on(0.0, 1.0).gauss(64)
    finally:
        _clear_caches()


@pytest.mark.foundation
@pytest.mark.rests_on(_G6, _G6C, _G7)
def test_the_weight_decision_refuses_an_unequal_tie_and_has_its_domain() -> None:
    r"""The decision "is this float64 midpoint the exact weight of this node?" (``gm._weight_is_exactly``),
    called directly: G3 stops on a NODE, so the unequal arm of the WEIGHT decision has no witness at consumer
    level.

    ``LEGENDRE.on(0.0, 1.0)`` at ``n = 3`` has weights ``5/18, 4/9, 5/18``.  At the centre node the tie ``4/9``
    is equal; at an outer node ``5/18`` is equal, ``4/9`` is unequal with a NON-trivial gcd (the node is not one
    of its roots) and ``1/3`` is unequal with a trivial gcd.  A point that is not a node does not isolate
    (``ArithmeticError``).  The nodes are built inside ``mpmath.workdps(40)`` and the decisions are called
    OUTSIDE it, at the ambient 15 digits (a conversion that re-rounds the node to the ambient precision fails the
    1e-20 bracket); a NEGATIVE node (the mirrored panel) decides alike (a conversion that drops the sign fails).
    The domain row: a Chebyshev T measure has an irrational mass, so its weight decision does not exist
    (``None``) while its node decision does."""
    import mpmath

    import orpheus.numerics.generating_measure as gm

    def decide(measure, node):
        return gm._weight_is_exactly(gm._rational_recurrence(measure, 3), 3, node, measure.name)

    unit, mirrored = LEGENDRE.on(0.0, 1.0), LEGENDRE.on(-1.0, 0.0)
    with mpmath.workdps(40):
        centre = mpmath.mpf(1) / 2
        outer = mpmath.mpf(1) / 2 + mpmath.sqrt(mpmath.mpf(3) / 5) / 2
        off_node = mpmath.mpf("0.3")
        negative_outer = -outer
    assert decide(unit, centre)(Fraction(4, 9)) is True
    assert decide(unit, outer)(Fraction(5, 18)) is True
    assert decide(unit, outer)(Fraction(4, 9)) is False  # gcd non-trivial, this node not a root of it
    assert decide(unit, outer)(Fraction(1, 3)) is False  # gcd trivial
    with pytest.raises(ArithmeticError, match="does not isolate"):
        decide(unit, off_node)(Fraction(4, 9))
    assert decide(mirrored, negative_outer)(Fraction(5, 18)) is True
    assert decide(mirrored, negative_outer)(Fraction(4, 9)) is False
    assert decide(mirrored, negative_outer)(Fraction(1, 3)) is False
    cheb = CHEBYSHEV_T.on(0.2, 0.4)
    rational = gm._rational_recurrence(cheb, 3)
    assert gm._weight_is_exactly(rational, 3, centre, cheb.name) is None
    assert gm._node_is_exactly(rational, 3) is not None


@pytest.mark.foundation
@pytest.mark.catches("ERR-127")
@pytest.mark.rests_on(_G4, _G6, _G6C, _G7)
@pytest.mark.parametrize(
    "label, measure, panel, n",
    [
        ("chebyshev_t n=1", CHEBYSHEV_T, (0.2, 0.4), 1),
        ("chebyshev_t n=3", CHEBYSHEV_T, (0.2, 0.4), 3),
        ("chebyshev_u n=3", CHEBYSHEV_U, (0.2, 0.4), 3),
        ("jacobi(1, 1) n=3", jacobi(1.0, 1.0), (0.2, 0.4), 3),
        ("jacobi(1.5, 1.5) n=5", jacobi(1.5, 1.5), (0.2, 0.4), 5),
        ("jacobi(0, 0) width tie n=1", jacobi(0.0, 0.0), _WIDTH_TIE, 1),
    ],
    ids=[
        "chebyshev_t n=1", "chebyshev_t n=3", "chebyshev_u n=3", "jacobi(1, 1) n=3",
        "jacobi(1.5, 1.5) n=5", "jacobi(0, 0) width tie n=1",
    ],
)
def test_every_family_rounds_its_exact_tie_half_even(label, measure, panel, n) -> None:
    r"""The law is not Legendre's: every family on ``[-1, 1]`` rounds its exact tie half to even.

    ``jacobi(1, 1)`` (mass ``4/3``) and ``jacobi(1.5, 1.5)`` (mass irrational) carry NO classical alias, so no
    Chebyshev construction can serve them (``jacobi(+-1/2, +-1/2)`` are the Chebyshevs by name and would test
    Chebyshev again); ``jacobi(0, 0)`` has the rational mass 2 and equals ``LEGENDRE`` by name, so on a width tie
    its WEIGHT must equal Legendre's.  Every row is a SYMMETRIC rule (the stabiliser of the node-weight pairing:
    reversing it would redden no row here, and an asymmetric rational weight tie is not known to exist).  Each row builds through ``_fresh`` (a remap cached for the equal constant
    would be returned without running the Jacobi builder), and the ``jacobi(0, 0)`` row builds its constant
    through ``_fresh`` too and asserts that two constructions ran (the positive control)."""
    a, b = panel
    try:
        built = _fresh(measure, panel)
        assert built.gauss(n).nodes[n // 2] == float(_centre(a, b))
        if label.startswith("jacobi(1"):
            assert built.name.startswith("jacobi(")  # no classical alias serves it
        if label.startswith("jacobi(0, 0)"):
            constant = _fresh(LEGENDRE, panel)
            jacobi_built = _fresh(measure, panel)
            assert jacobi_built is not constant
            assert jacobi_built.exact_recurrence is not constant.exact_recurrence
            assert jacobi_built.gauss(n).weights[0] == constant.gauss(n).weights[0] == float(2 * _half_width(a, b))
    finally:
        _clear_caches()


_BUILDER_ROWS = [
    # (label, builder from the recurrence module, the mass beta_0 in RATIONAL)
    ("legendre", lambda r: r.legendre_recurrence, Fraction(2)),
    ("chebyshev_t", lambda r: r.chebyshev_t_recurrence, None),
    ("chebyshev_u", lambda r: r.chebyshev_u_recurrence, None),
    ("hermite", lambda r: r.hermite_recurrence, None),
    ("jacobi(0, 0)", lambda r: r.jacobi_recurrence(0.0, 0.0), Fraction(2)),
    ("jacobi(2, 3)", lambda r: r.jacobi_recurrence(2.0, 3.0), Fraction(16, 15)),
    ("jacobi(-0.5, 0.5)", lambda r: r.jacobi_recurrence(-0.5, 0.5), None),
    ("jacobi(1.5, 2.25)", lambda r: r.jacobi_recurrence(1.5, 2.25), None),
    ("laguerre(0)", lambda r: r.laguerre_recurrence(0.0), Fraction(1)),
    ("laguerre(2)", lambda r: r.laguerre_recurrence(2.0), Fraction(2)),
    ("laguerre(1.5)", lambda r: r.laguerre_recurrence(1.5), None),
]


@pytest.mark.foundation
@pytest.mark.parametrize("label, getter, mass", _BUILDER_ROWS, ids=[row[0] for row in _BUILDER_ROWS])
def test_the_generic_builders_admit_the_tabulated_rationals(label, getter, mass) -> None:
    r"""The six closed forms are ONE definition each, evaluated in two number types (2.0.4).

    (a) The DOMAIN: in ``RATIONAL`` every ``alpha_k`` and every ``beta_k`` (``k >= 1``) is a ``Fraction`` and the
    mass ``beta_0`` is exactly the tabulated value (``None`` where it is a multiple of ``pi`` or ``sqrt(pi)``, or
    a Gamma value of a non-integer).  (b) The AGREEMENT of the two arithmetics to 1e-55 relative at 60 digits:
    both evaluate one body, so this is tautological in the body and pins the TABLE (the rational evaluation takes
    the same branch as the high-precision one).  The remap's law is the next test."""
    import mpmath

    import orpheus.numerics.recurrence as r

    builder = getter(r)
    alpha, beta = builder(6, r.RATIONAL)
    assert len(alpha) == 6 and len(beta) == 6
    assert all(isinstance(x, Fraction) for x in alpha)
    assert all(isinstance(x, Fraction) for x in beta[1:])
    assert beta[0] == mass and (beta[0] is None) == (mass is None)
    high_alpha, high_beta = builder(6, r.HighPrecision(60))
    with mpmath.workdps(60):
        for q, h in list(zip(alpha, high_alpha)) + list(zip(beta[1:], high_beta[1:])):
            exact = mpmath.mpf(q.numerator) / q.denominator
            assert (abs(exact - h) <= mpmath.mpf(10) ** -55 * abs(h)) if h != 0 else (exact == 0)


@pytest.mark.foundation
@pytest.mark.parametrize("label, getter, mass", _BUILDER_ROWS, ids=[row[0] for row in _BUILDER_ROWS])
def test_a_generator_evaluates_inside_its_arithmetic_never_at_the_ambient_precision(label, getter, mass) -> None:
    """The door is nominal: every builder is a ``RecurrenceCoefficients``, whose call enters the arithmetic's
    working context.  The positive control is the undecorated ``closed_form`` called at the ambient 15 digits,
    which carries ``1/3`` at 1e-16 only; the decorated call carries it to 1e-55 at 60 digits."""
    import mpmath

    builder = getter(r)
    assert isinstance(builder, r.RecurrenceCoefficients)
    _, beta = r.legendre_recurrence(3, r.HighPrecision(60))
    _, bare = r.legendre_recurrence.closed_form(3, r.HighPrecision(60))
    with mpmath.workdps(60):
        third = mpmath.mpf(1) / 3
        assert abs(beta[1] - third) < mpmath.mpf(10) ** -55
        assert abs(bare[1] - third) > mpmath.mpf(10) ** -30


@pytest.mark.foundation
def test_the_rational_door_admits_fractions_only() -> None:
    """``_rational_recurrence`` parses what it lets through: a recurrence whose RATIONAL evaluation leaks a float
    or an ``mpf`` is outside the domain (``None``, so the value is refused as before), never decided on.  The
    positive control is ``LEGENDRE``, which the door admits."""
    import dataclasses

    import mpmath

    import orpheus.numerics.generating_measure as gm

    assert gm._rational_recurrence(LEGENDRE, 3) is not None
    for leak in (0.5, mpmath.mpf(1) / 2):
        leaking = dataclasses.replace(
            LEGENDRE,
            name=f"legendre-leaking-{type(leak).__name__}",  # a distinct name: one name is one content (#711)
            exact_recurrence=r.over_arithmetic(
                lambda n, A, leak=leak: ([A.lift(0)] * n, [A.mass(lambda: A.lift(2))] + [leak] * (n - 1))
            ),
        )
        assert gm._rational_recurrence(leaking, 3) is None


@pytest.mark.foundation
@pytest.mark.parametrize("panel", [(0.2, 0.4), _WIDTH_TIE], ids=str)
def test_the_remap_maps_the_rational_recurrence_exactly(panel) -> None:
    r"""The affine remap's law in ``RATIONAL`` against the moments of the uniform measure on ``[a, b]`` (an
    oracle that is not a recurrence): ``alpha'_0 = (a + b)/2`` (the mean), ``beta'_0 = b - a`` (the mass),
    ``beta'_1 = (b - a)^2/12`` (the variance)."""
    from orpheus.numerics.recurrence import RATIONAL

    a, b = panel
    try:
        alpha, beta = _fresh(LEGENDRE, panel).exact_recurrence(3, RATIONAL)
    finally:
        _clear_caches()
    assert alpha[0] == _centre(a, b)
    assert beta[0] == Fraction(b) - Fraction(a)
    assert beta[1] == (Fraction(b) - Fraction(a)) ** 2 / 12


# =============================================================================
# Re-filed from tests/gates/geometry/test_manifold.py (layer DAG ruling 23, 2026-10-25)
# =============================================================================


@pytest.mark.foundation
class TestTheEntryCarriesItsQuotientMapAndItsReference:
    """The numerics half of geometry's ``TestTheEntryCarriesItsQuotientMapAndItsReference``.

    Re-filed 2026-10-25 from ``tests/gates/geometry/test_manifold.py`` (layer-DAG
    ruling 23): the row below asserts a numerics object's verb, so it lives here; the class keeps its geometry name so the node id
    changes only in its path.  The manifold rows of the class stay in geometry.
    """

    _AXES = list(_g22b_AXES)

    def test_the_axial_entries_have_the_hat_box_reference_and_the_others_None(
        self,
    ) -> None:
        r"""Archimedes: :math:`\mu_* d\Omega = 2\pi\,d\mu`, so the reference
        is Lebesgue on :math:`[-1,1]` — ``LEGENDRE``, identically, on all
        three axes.  ``None`` on the mirror (a weighted disk measure no
        shipped realization spells) and on the trivial quotient (Lebesgue
        on the BASE, whose orthogonal system the type does not carry)."""
        from orpheus.numerics.generating_measure import LEGENDRE, pushforward_reference

        for axis in self._AXES:
            assert pushforward_reference(SPHERE.quotient(SubgroupOfO3.O2(axis))) is LEGENDRE
        assert LEGENDRE.support == COSINE_INTERVAL
        for axis in self._AXES:
            assert pushforward_reference(SPHERE.quotient(SubgroupOfO3.Mirror(axis))) is None
        assert pushforward_reference(SPHERE.quotient(SubgroupOfO3.Trivial)) is None
        assert pushforward_reference(CIRCLE.quotient(SubgroupOfO3.Trivial)) is None
        # S^2/O(3) is the point (#440): the pushforward is a point mass of 4 pi, which
        # no shipped reference spells.
        assert pushforward_reference(SPHERE.quotient(SubgroupOfO3.O3)) is None
        # ⚠ A reference lives on the REALIZATION — the chart codomain, where
        # the quotient map's image lands and where 2 pi d(mu) is written —
        # never on the entry, which would demand a measure carry an axis it
        # does not know.  `[M]` 2026-09-02 (archivist): nothing in orpheus/
        # reads `reference.support`, so this is the only pin on it.
        for group in _ROSTER:
            q = SPHERE.quotient(group)
            reference = pushforward_reference(q)
            if reference is not None:
                assert reference.support == q.realization
                assert reference.support != q
