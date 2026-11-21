r"""Verification for ``orpheus.numerics.rational_polynomial``: polynomials over Q as tuples of ``Fraction``,
lowest degree first, with the Euclidean algorithm the exact-tie decision of the Gauss rules reads
(unit ``quadrature-ties``, 2.0.5 E1 of ``.claude/plans/ontology_w0_spec.md``).

**foundation.** Every function returns ``Fraction`` entries whatever it is given (``trim`` is the parse door), so an
implementation that computes in floats, or that floor-divides integers, is red.  The oracle of ``evaluate`` is
``numpy.polynomial.Polynomial`` at a dyadic point (exact in float64), which shares no code with the module.

The module is imported inside a fixture, so that each test reports ``ModuleNotFoundError`` on its own until the
module exists, and a collection error does not hide the other files of the directory.
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest

pytestmark = [pytest.mark.foundation]


@pytest.fixture(scope="module")
def rp():
    import orpheus.numerics.rational_polynomial as module

    return module


def _poly(*coefficients):
    return tuple(Fraction(c) for c in coefficients)


# p = 2 (x - 1)(x - 2) = 2 x^2 - 6 x + 4 ; q = (x - 1)(x - 3) = x^2 - 4 x + 3   (lowest degree first)
P = _poly(4, -6, 2)
Q = _poly(3, -4, 1)


def _all_fractions(poly) -> bool:
    return all(isinstance(c, Fraction) for c in poly)


def test_remainder_is_the_long_division_remainder_over_q(rp) -> None:
    r"""``p = 2 q + r`` with ``r = 2 x - 2`` ... precisely ``p - 2 q = (-6 + 8) x + (4 - 6) = 2 x - 2``, i.e.
    ``(-2, 2)``.  A second case, ``x^3 + 1`` by ``x^2 + 1``, leaves ``1 - x``; the entries are ``Fraction`` (an
    integer division would floor, a float division would not be exact)."""
    r = rp.remainder(P, Q)
    assert r == _poly(-2, 2) and _all_fractions(r)
    r2 = rp.remainder(_poly(1, 0, 0, 1), _poly(1, 0, 1))
    assert r2 == _poly(1, -1) and _all_fractions(r2)


def test_remainder_by_a_constant_is_zero_and_by_zero_is_refused(rp) -> None:
    """A degree-0 divisor divides exactly (remainder the zero polynomial); the zero polynomial does not divide."""
    assert rp.remainder(_poly(1, 0, 1), _poly(1)) == rp.ZERO
    with pytest.raises(ZeroDivisionError):
        rp.remainder(P, rp.ZERO)


def test_gcd_is_the_monic_greatest_common_divisor(rp) -> None:
    """``gcd(2 (x-1)(x-2), (x-1)(x-3)) = x - 1``: monic, so ``(-1, 1)`` (the normalisation is part of the value:
    the unnormalised Euclid ends at ``(-2, 2)``)."""
    g = rp.gcd(P, Q)
    assert g == _poly(-1, 1) and _all_fractions(g)
    assert rp.gcd(_poly(1, 0, 1), _poly(2, 1)) == _poly(1)  # coprime: the constant 1


def test_every_function_returns_fractions_whatever_it_is_given(rp) -> None:
    """``trim`` is the parse door: ``int`` and ``float`` inputs are exact rationals, leading zeros go, and every
    output entry is a ``Fraction``."""
    assert rp.trim([1, 2, 0, 0]) == _poly(1, 2) and _all_fractions(rp.trim([1, 2, 0, 0]))
    assert rp.trim([0.5, 0.25]) == (Fraction(1, 2), Fraction(1, 4))
    assert rp.trim([0, 0]) == rp.ZERO
    for out in (
        rp.add((1, 2), (3,)),
        rp.scale((1, 2), 3),
        rp.multiply((1, 1), (1, -1)),
        rp.remainder((1, 0, 0, 1), (1, 0, 1)),
        rp.gcd((4, -6, 2), (3, -4, 1)),
    ):
        assert _all_fractions(out)
    assert rp.multiply((1, 1), (-1, 1)) == _poly(-1, 0, 1)
    assert rp.add((1, 2), (-1, -2)) == rp.ZERO


def test_evaluate_agrees_with_numpy_at_a_dyadic_point(rp) -> None:
    """The independent oracle: numpy's ``Polynomial`` at ``x = 3/8`` (every product is exact in float64)."""
    x = Fraction(3, 8)
    for poly in (P, Q, _poly(1, 0, 0, 1), _poly(Fraction(1, 2), Fraction(-3, 4))):
        value = rp.evaluate(poly, x)
        assert isinstance(value, Fraction)
        assert float(value) == np.polynomial.Polynomial([float(c) for c in poly])(float(x))
