"""Foundation + L1 tests for ``orpheus.numerics.quadrature.rules_1d``.

Two layers:

* **Foundation** — software invariants of the rule function (shape,
  metadata fields, bit-identical match against the legacy
  :class:`~orpheus.sn.quadrature.GaussLegendre1D` constructor that
  the regression snapshots pin).
* **L1** — closed-form polynomial-exactness theorem of Stoer &
  Bulirsch (2002) Theorem 3.6.20: the :math:`n`-point Gauss-Legendre
  rule integrates polynomials of degree :math:`\\le 2n - 1` exactly.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.measure import SPACE_INTERVAL_M11, DiscreteMeasure
from orpheus.numerics.quadrature import gauss_legendre_on_mu
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.symmetry import SubgroupOfO3


# ---------------------------------------------------------------------------
# Foundation — software invariants
# ---------------------------------------------------------------------------


@pytest.mark.foundation
@pytest.mark.parametrize("n", [1, 2, 4, 8, 16, 32])
def test_gauss_legendre_returns_discrete_measure(n: int) -> None:
    """The rule MUST return a :class:`DiscreteMeasure` with correct
    shape, space tag, and ``invariance_group`` / ``degree_of_exactness``
    fields populated."""
    m = gauss_legendre_on_mu(n)
    assert isinstance(m, DiscreteMeasure)
    assert m.n_points == n
    assert m.nodes.shape == (n,)
    assert m.weights.shape == (n,)
    assert m.support == SPACE_INTERVAL_M11
    assert m.invariance_group == SubgroupOfO3.SO2
    assert m.degree_of_exactness == 2 * n - 1


@pytest.mark.foundation
@pytest.mark.parametrize("n", [1, 2, 4, 8])
def test_gauss_legendre_weight_sum_is_two(n: int) -> None:
    """Weight sum must equal :math:`2 = \\int_{-1}^{1} dx`."""
    m = gauss_legendre_on_mu(n)
    assert m.weights.sum() == pytest.approx(2.0, abs=1e-15)


@pytest.mark.foundation
def test_gauss_legendre_rejects_zero() -> None:
    """``n < 1`` must raise ``ValueError``."""
    with pytest.raises(ValueError):
        gauss_legendre_on_mu(0)
    with pytest.raises(ValueError):
        gauss_legendre_on_mu(-3)


@pytest.mark.foundation
@pytest.mark.parametrize("n", [2, 4, 8, 16, 32])
def test_adapter_passes_the_rule_through_unmodified(n: int) -> None:
    """The adapter must not transform the values it wraps.

    A real contract, but a NARROW one — and narrower than this test
    used to claim. ``Quadrature.gauss_legendre`` *calls*
    ``gauss_legendre_on_mu`` (``directional.py``), so comparing the two
    compares a value with itself routed through a wrapper. It verifies
    the adapter is a pass-through; it can NEVER detect the node drift
    it was previously documented as catching, because both sides move
    together by construction. The drift gate is the separate test
    below.
    """
    m = gauss_legendre_on_mu(n)
    adapter = Quadrature.gauss_legendre(n)
    assert np.array_equal(adapter.mu_x, m.nodes)
    assert np.array_equal(adapter.weights, m.weights)
    # mu_y must be identically zero (1-D slab convention).
    assert np.array_equal(adapter.mu_y, np.zeros(n))


@pytest.mark.foundation
@pytest.mark.parametrize("n", [2, 4, 8, 16, 32])
def test_nodes_are_bit_identical_to_the_reference_implementation(
    n: int,
) -> None:
    """The **load-bearing** drift contract, against an EXTERNAL source.

    Every SN regression snapshot built on a slab quadrature is pinned
    to these exact bits. A change of construction that moves them —
    swapping ``numpy.leggauss`` for the Golub-Welsch body in
    :mod:`orpheus.numerics.generating_measure`, say, which `[M]` moves
    nodes by 1-4 ULP — must therefore be a deliberate re-baselining
    decision, not a silent side effect.

    So this compares against ``numpy.polynomial.legendre.leggauss``
    directly: a source that does NOT move when ORPHEUS changes, which
    is what makes the comparison capable of failing.
    """
    m = gauss_legendre_on_mu(n)
    ref_nodes, ref_weights = np.polynomial.legendre.leggauss(n)
    assert np.array_equal(m.nodes, ref_nodes), (
        f"n={n}: slab quadrature nodes moved off the reference "
        f"implementation; max delta "
        f"{np.max(np.abs(m.nodes - ref_nodes)):.3e}. Every SN snapshot "
        f"built on this rule is pinned to the old bits — re-baseline "
        f"deliberately per vv-principles, do not relax this gate."
    )
    assert np.array_equal(m.weights, ref_weights)


# ---------------------------------------------------------------------------
# L1 — closed-form polynomial-exactness theorem
# ---------------------------------------------------------------------------


@pytest.mark.l1
@pytest.mark.parametrize("n", [2, 4, 8, 16])
def test_gauss_legendre_polynomial_exactness(n: int) -> None:
    """Stoer & Bulirsch (2002) Theorem 3.6.20: the ``n``-point GL rule
    integrates :math:`x^k` exactly for :math:`k \\le 2n - 1`.

    Closed-form reference: :math:`\\int_{-1}^{1} x^k dx = 0` for odd
    :math:`k`, :math:`= 2/(k+1)` for even :math:`k`. The rule's
    ``degree_of_exactness`` is the upper bound this test checks.
    """
    m = gauss_legendre_on_mu(n)

    for k in range(2 * n):  # 0..2n-1 inclusive
        # Closed-form ∫_{-1}^{1} x^k dx
        if k % 2 == 1:
            exact = 0.0
        else:
            exact = 2.0 / (k + 1)

        approx = m.integrate(lambda x, k=k: x**k)
        assert approx == pytest.approx(exact, abs=1e-13), (
            f"GL{n} fails at degree {k}: rule returns {approx}, "
            f"closed form is {exact}"
        )


@pytest.mark.l1
def test_gauss_legendre_fails_above_degree_of_exactness() -> None:
    """At :math:`k = 2n` the rule's exactness theorem makes no claim;
    we check that the rule is in fact NOT exact at :math:`k = 2n` for
    a small ``n`` (otherwise the test below would be vacuously
    satisfied and not actually probing the boundary)."""
    n = 4
    m = gauss_legendre_on_mu(n)
    k = 2 * n  # 8 — beyond exactness
    exact = 2.0 / (k + 1)  # 2/9
    approx = m.integrate(lambda x: x**k)
    # The rule must NOT be exact here — otherwise the boundary is wrong.
    assert abs(approx - exact) > 1e-6
