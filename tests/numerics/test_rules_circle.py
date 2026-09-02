r"""Foundation + L1 tests for ``orpheus.numerics.quadrature.rules_circle``.

Two layers:

* **Foundation** — software invariants of the rule function (shape,
  space tag, mass, the shift's fundamental domain, the guards) plus the
  two *exactness* properties that are exact-arithmetic facts rather than
  numerical ones: an on-axis node's :math:`\sin` is exactly ``0.0``, and
  the mirror is an exact index permutation.
* **L1** — the closed-form aliasing theorem
  :eq:`periodic-trapezoid-aliasing`: the :math:`n`-point periodic
  trapezoid integrates :math:`e^{ik\varphi}` exactly for every
  :math:`|k| \le n - 1` and misses by the full :math:`2\pi` at
  :math:`k = \pm n`, independently of the shift.

Test functions are built as :math:`z^k` with
:math:`z = \cos\varphi + i\sin\varphi` read straight off the nodes —
deliberately NOT as ``np.cos(k * arctan2(...))``. Recovering the angle
and re-evaluating a transcendental would route the check through the
very chart the rule avoids, and would make the test's own error floor
the thing being measured.
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest

from orpheus.numerics.exactness import OrthogonalSystem
from orpheus.numerics.measure import DiscreteMeasure, equispaced
from orpheus.numerics.manifold import CIRCLE
from orpheus.numerics.quadrature.rules_circle import (
    NODE_ALIGNED,
    STAGGERED,
    periodic_trapezoid,
)

# Shifts with no mirror symmetry at all — the generic case the
# classification excludes, kept in the sweep so "exactness is
# shift-invisible" is tested outside the two special values.
_GENERIC_SHIFTS = [Fraction(1, 3), Fraction(2, 7), Fraction(3, 5)]
_ALL_SHIFTS = [NODE_ALIGNED, STAGGERED, *_GENERIC_SHIFTS]


def _unit_complex(rule: DiscreteMeasure) -> np.ndarray:
    r""":math:`z_m = e^{i\varphi_m}`, read off the node components."""
    return rule.nodes[:, 0] + 1j * rule.nodes[:, 1]


# ---------------------------------------------------------------------------
# Foundation — software invariants
# ---------------------------------------------------------------------------


@pytest.mark.foundation
@pytest.mark.parametrize("n", [1, 2, 3, 4, 8, 16])
@pytest.mark.parametrize("shift", _ALL_SHIFTS)
def test_shape_space_and_mass(n: int, shift: Fraction) -> None:
    """A rule on :math:`S^1` with :math:`n` points on the unit circle,
    carrying the mass of the uniform measure it claims exactness
    against."""
    rule = periodic_trapezoid(n, shift=shift)
    assert isinstance(rule, DiscreteMeasure)
    assert rule.n_points == n
    assert rule.nodes.shape == (n, 2)
    assert rule.support == CIRCLE
    np.testing.assert_allclose(
        np.linalg.norm(rule.nodes, axis=1), 1.0, rtol=0, atol=4e-16,
    )
    np.testing.assert_allclose(rule.weights.sum(), 2 * np.pi, rtol=1e-15)


@pytest.mark.foundation
@pytest.mark.parametrize("n", [1, 2, 5, 8])
@pytest.mark.parametrize("shift", _ALL_SHIFTS)
def test_the_claim_names_the_circle_and_its_system(
    n: int, shift: Fraction
) -> None:
    """Degree ``n-1`` against the uniform measure on :math:`S^1`, in the
    TRIGONOMETRIC system — all three parts, since reading any one alone
    is what the claim type exists to prevent."""
    claim = periodic_trapezoid(n, shift=shift).exactness
    assert claim is not None
    assert claim.degree == n - 1
    assert claim.system is OrthogonalSystem.TRIGONOMETRIC
    assert claim.reference.support == CIRCLE


@pytest.mark.foundation
def test_a_single_point_rule_claims_only_the_constant() -> None:
    """``n = 1`` is degree 0 — a real, weak claim, not an absent one.

    Positive control for the ``degree >= 0`` guard on the claim type:
    the smallest legitimate rule must not be pushed into ``None``.
    """
    claim = periodic_trapezoid(1, shift=NODE_ALIGNED).exactness
    assert claim is not None and claim.degree == 0


@pytest.mark.foundation
def test_no_invariance_group_is_declared() -> None:
    """The tag would be a claim about an embedding this rule has not
    chosen — and this package has shipped three false symmetry
    declarations (ERR-072/073/074). Absence here is the discipline, so
    it is pinned rather than left to drift back in."""
    assert periodic_trapezoid(8, shift=STAGGERED).invariance_group is None


# ---------------------------------------------------------------------------
# Foundation — the shift's fundamental domain
# ---------------------------------------------------------------------------


@pytest.mark.foundation
@pytest.mark.parametrize(
    "equivalent", [Fraction(5, 2), Fraction(-1, 2), Fraction(101, 2)]
)
def test_a_whole_step_is_the_identity(equivalent: Fraction) -> None:
    """The shift lives in :math:`\\mathbb{Q}/\\mathbb{Z}`, so the
    reduction canonicalises a representative rather than restricting the
    caller — asserted BIT-identically, because the reduction happens in
    exact rational arithmetic before any float is produced."""
    canonical = periodic_trapezoid(6, shift=STAGGERED)
    relabelled = periodic_trapezoid(6, shift=equivalent)
    assert np.array_equal(canonical.nodes, relabelled.nodes)
    assert np.array_equal(canonical.weights, relabelled.weights)


@pytest.mark.foundation
def test_guards() -> None:
    with pytest.raises(ValueError, match="n >= 1"):
        periodic_trapezoid(0, shift=NODE_ALIGNED)
    with pytest.raises(ValueError, match="Fraction shift"):
        periodic_trapezoid(8, shift=0.5)  # type: ignore[arg-type]


@pytest.mark.foundation
def test_a_float_shift_is_refused_rather_than_converted() -> None:
    r"""``0.5`` is exactly representable and would convert cleanly — but
    ``0.1`` would convert to :math:`3602879701896397/36028797018963968`,
    and the rule would then ask for a root of unity of that order.

    So the guard is on the TYPE, not on convertibility: rejecting the
    benign float is what keeps the malignant one from arriving as a
    32-digit denominator.
    """
    assert Fraction(0.5) == STAGGERED           # the benign case...
    assert Fraction(0.1).denominator > 10**16   # ...and why it is refused
    with pytest.raises(ValueError, match="Fraction shift"):
        periodic_trapezoid(8, shift=0.1)  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# L1 — the aliasing theorem, and its tightness
# ---------------------------------------------------------------------------


@pytest.mark.l1
@pytest.mark.parametrize("n", [1, 2, 3, 4, 8, 12])
@pytest.mark.parametrize("shift", _ALL_SHIFTS)
def test_exact_on_every_fourier_mode_up_to_the_claimed_degree(
    n: int, shift: Fraction
) -> None:
    r"""The claim, measured mode by mode: :math:`\int e^{ik\varphi}` is
    :math:`2\pi` at :math:`k = 0` and :math:`0` otherwise, reproduced for
    every :math:`|k| \le n - 1`."""
    rule = periodic_trapezoid(n, shift=shift)
    z = _unit_complex(rule)
    claim = rule.exactness
    assert claim is not None

    for k in range(-claim.degree, claim.degree + 1):
        quadrature = complex(rule.weights @ (z**k))
        exact = 2 * np.pi if k == 0 else 0.0
        assert abs(quadrature - exact) < 1e-13, f"n={n} shift={shift} k={k}"


@pytest.mark.l1
@pytest.mark.parametrize("n", [2, 3, 4, 8, 12])
@pytest.mark.parametrize("shift", _ALL_SHIFTS)
def test_the_degree_is_TIGHT_the_next_mode_aliases_onto_the_constant(
    n: int, shift: Fraction
) -> None:
    r"""At :math:`k = n` the rule returns the full :math:`2\pi` where the
    exact integral is :math:`0` — the bound is attained, not estimated.

    Without this leg the degree claim would be a lower bound that any
    smaller number also satisfies; with it, ``n - 1`` is the answer.
    """
    rule = periodic_trapezoid(n, shift=shift)
    z = _unit_complex(rule)
    aliased = complex(rule.weights @ (z**n))
    np.testing.assert_allclose(abs(aliased), 2 * np.pi, rtol=1e-13)


@pytest.mark.l1
@pytest.mark.parametrize("n", [4, 8, 9])
def test_the_shift_is_exactness_INVISIBLE(n: int) -> None:
    r"""Every shift yields the same claim, because the shift enters
    :eq:`periodic-trapezoid-aliasing` only as a phase multiplying a
    factor that is already zero.

    Pinned as its own law because the *next* test is its converse: the
    parameter no exactness gate can see is the one that decides
    :math:`\Sigma`. A gate over exactness alone would license changing
    it freely.
    """
    claims = {periodic_trapezoid(n, shift=s).exactness for s in _ALL_SHIFTS}
    assert len(claims) == 1


# ---------------------------------------------------------------------------
# L1 — the mirror classification (the shift's actual job)
# ---------------------------------------------------------------------------


@pytest.mark.l1
@pytest.mark.parametrize("n", [2, 3, 4, 5, 8, 9])
@pytest.mark.parametrize(
    "shift,partner",
    [
        (NODE_ALIGNED, lambda n, m: (n - m) % n),
        (STAGGERED, lambda n, m: n - 1 - m),
    ],
)
@pytest.mark.l1
def test_the_two_mirror_symmetric_shifts_reflect_BIT_exactly(
    n: int, shift: Fraction, partner
) -> None:
    r"""For :math:`s \in \{0, \tfrac{1}{2}\}` the reflection
    :math:`\varphi \to -\varphi` is an exact index permutation of the
    node set.

    ``array_equal``, not ``allclose``: the whole reason the nodes are
    generated by the group action rather than by ``cos``/``sin`` of an
    angle is that this holds to the last bit. A tolerance here would
    pass just as well against the construction the rule replaced, and so
    would test nothing.
    """
    rule = periodic_trapezoid(n, shift=shift)
    cos, sin = rule.nodes[:, 0], rule.nodes[:, 1]
    perm = np.array([partner(n, m) for m in range(n)])
    assert np.array_equal(cos, cos[perm])
    assert np.array_equal(sin, -sin[perm])


@pytest.mark.l1
@pytest.mark.parametrize("shift", _GENERIC_SHIFTS)
def test_a_generic_shift_has_NO_mirror_at_all(shift: Fraction) -> None:
    r"""The classification is exhaustive: outside
    :math:`s \in \{0, \tfrac{1}{2}\}` the reflected set is not the
    original set — so the two named shifts are theorems, not two
    conventions someone picked.

    Negative control for the test above, which on its own would be
    satisfied by a rule that is mirror-symmetric for *every* shift.
    """
    rule = periodic_trapezoid(8, shift=shift)
    reflected = np.column_stack([rule.nodes[:, 0], -rule.nodes[:, 1]])
    # Distance from each reflected node to its nearest original node.
    gaps = np.linalg.norm(
        reflected[:, None, :] - rule.nodes[None, :, :], axis=2
    ).min(axis=1)
    assert gaps.max() > 1e-3


# ---------------------------------------------------------------------------
# L1 — Σ = {ξ = 0} is decided by an EQUALITY, not by a tolerance
# ---------------------------------------------------------------------------


@pytest.mark.l1
@pytest.mark.parametrize("n", [2, 4, 6, 8, 16])
def test_a_node_aligned_rule_puts_nodes_EXACTLY_on_the_axis(n: int) -> None:
    r""":math:`\Sigma \neq \emptyset`, and the membership is exact.

    ``sin`` is exactly ``0.0`` at :math:`m = 0` and :math:`m = n/2` —
    contrast ``np.sin(np.pi) = 1.22e-16``, asserted alongside so the
    comparison is in the test rather than only in the prose. A singular
    set defined by ``abs(sin) < eps`` would have its well-posedness
    condition set by the choice of ``eps``.
    """
    sin = periodic_trapezoid(n, shift=NODE_ALIGNED).nodes[:, 1]
    on_axis = sin == 0.0
    assert on_axis.sum() == 2
    assert on_axis[0] and on_axis[n // 2]
    assert float(np.sin(np.pi)) != 0.0  # the construction being replaced


@pytest.mark.l1
@pytest.mark.parametrize("n", [2, 4, 6, 8, 16])
def test_a_staggered_EVEN_rule_has_an_EMPTY_singular_set(n: int) -> None:
    r""":math:`\Sigma = \emptyset` — the fold's well-posedness condition,
    selected rather than flagged."""
    sin = periodic_trapezoid(n, shift=STAGGERED).nodes[:, 1]
    assert (sin != 0.0).all()


@pytest.mark.l1
@pytest.mark.parametrize("n", [3, 5, 7, 9])
def test_a_staggered_ODD_rule_still_meets_the_axis(n: int) -> None:
    r"""The parity caveat, pinned: at odd :math:`n` the node
    :math:`m = (n-1)/2` sits at :math:`\varphi = \pi`, so staggering does
    NOT empty :math:`\Sigma`.

    Without this leg the previous test would read as "staggered ⟹
    :math:`\Sigma = \emptyset`", and a consumer selecting the shift on
    that basis would be silently wrong for every odd rule.
    """
    sin = periodic_trapezoid(n, shift=STAGGERED).nodes[:, 1]
    assert (sin == 0.0).sum() == 1
    assert sin[(n - 1) // 2] == 0.0


# ---------------------------------------------------------------------------
# ⭐ The load-bearing distinction: same nodes, different objects
# ---------------------------------------------------------------------------


@pytest.mark.l1
@pytest.mark.parametrize("n", [4, 8, 16])
def test_the_SAME_nodes_on_an_interval_carry_an_INCOMPARABLE_claim(
    n: int,
) -> None:
    r"""A rule on the circle and a rule on an interval are different
    objects even when their nodes coincide.

    This is the second measured bug in the :mod:`exactness` header, made
    into a permanent gate: the staggered circle rule and
    ``equispaced(0, 2π, n)`` place their nodes at the same angles, yet
    one is exact to trigonometric degree :math:`n-1` and the other to
    algebraic degree :math:`1`. Reading the first as the second is what
    made a naive ``min()`` report degree 1 for the sphere product rule.
    """
    circle = periodic_trapezoid(n, shift=STAGGERED)
    interval = equispaced(0.0, 2 * np.pi, n)

    # Same points — the coincidence that makes the confusion possible...
    np.testing.assert_allclose(
        circle.nodes[:, 0], np.cos(interval.nodes), rtol=0, atol=1e-15
    )
    np.testing.assert_allclose(
        circle.nodes[:, 1], np.sin(interval.nodes), rtol=0, atol=1e-15
    )

    # ...and three ways the claims are not the same claim.
    circle_claim, interval_claim = circle.exactness, interval.exactness
    assert circle_claim is not None and interval_claim is not None
    assert circle_claim.system is not interval_claim.system
    assert circle_claim.reference.support != interval_claim.reference.support
    assert circle_claim.degree == n - 1 and interval_claim.degree == 1


@pytest.mark.l1
@pytest.mark.parametrize("n", [4, 8, 16])
def test_the_interval_rule_is_NOT_exact_at_the_circle_rule_degree(
    n: int,
) -> None:
    r"""The converse leg, and the one that shows the distinction has
    teeth rather than being a labelling preference.

    The interval rule's own system is algebraic, so its degree-``n-1``
    claim — were the tags naively swapped — would assert exactness on
    :math:`x^{n-1}` against :math:`dx` on :math:`[0, 2\pi]`. Measured, it
    misses by a wide margin, which is why the two integers must never be
    compared.
    """
    interval = equispaced(0.0, 2 * np.pi, n)
    quadrature = float(interval.weights @ interval.nodes ** (n - 1))
    exact = (2 * np.pi) ** n / n
    assert abs(quadrature - exact) / exact > 1e-2
