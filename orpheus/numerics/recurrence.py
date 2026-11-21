r"""The three-term recurrence of an orthogonal family: the Jacobi matrix.

The monic orthogonal polynomials :math:`\{\pi_k\}` of a positive measure
:math:`d\lambda` on the real line obey

.. math::

   \pi_{k+1}(x) = (x - \alpha_k)\,\pi_k(x) - \beta_k\,\pi_{k-1}(x),
   \qquad \pi_{-1} = 0,\ \pi_0 = 1,

and this module stores exactly those coefficients: each function returns
``(alpha, beta)``, two lists of ``n`` numbers (``mpmath.mpf`` or ``fractions.Fraction``, by the arithmetic it is given), with the convention
``beta[0]`` :math:`= \mu_0 = \int d\lambda`, the zeroth moment (Gautschi 2004
§1.3). The orthonormal family :math:`p_k = \pi_k / \|\pi_k\|` obeys the
symmetric form

.. math::

   x\, p_k(x) = a_k\, p_{k+1}(x) + b_k\, p_k(x) + a_{k-1}\, p_{k-1}(x),
   \qquad b_k = \alpha_k,\quad a_k = \sqrt{\beta_{k+1}},

whose coefficients are the entries of the symmetric tridiagonal **Jacobi
matrix** (:math:`b_k` on the diagonal, :math:`a_k` off it). The recurrence is
the measure's spectral data: Favard's theorem recovers the measure from it,
and the :math:`n`-point Gauss rule is the eigen-decomposition of the
:math:`n \times n` Jacobi matrix, with weights :math:`\mu_0 v_{0j}^2` from the
first components of its eigenvectors (Golub & Welsch 1969). It is therefore
upstream of both the continuous generating measure
(:mod:`orpheus.numerics.generating_measure`) and the discrete Gauss measure
derived from it, and imports nothing from :mod:`orpheus.numerics`.

Arithmetic. Every closed form takes the arithmetic it is evaluated in as an
argument, never mpmath's ambient setting: :class:`HighPrecision` (an explicit
number of decimal digits) or :data:`RATIONAL` (exact ``fractions.Fraction``).
The closed forms are written once, each evaluated in either; every float
coefficient and every Gauss rule is derived from them. In ``RATIONAL`` the
mass ``beta[0]`` is ``None`` when it is irrational (:math:`\pi`,
:math:`\sqrt\pi`): the one place a transcendental number appears.

References
----------
- Golub, G. H. & Welsch, J. H. (1969). Calculation of Gauss quadrature
  rules. *Math. Comp.* 23, 221-230.
- Gautschi, W. (2004). *Orthogonal Polynomials: Computation and
  Approximation*. Oxford. §1.3 (recurrence coefficients of the classical
  families) and §3.1.1 (the Jacobi matrix and Gauss rules).
"""

from __future__ import annotations

import contextlib
import math
from collections.abc import Callable
from dataclasses import dataclass
from fractions import Fraction
from typing import Any, ContextManager, Protocol, TypeVar

import mpmath

Number = mpmath.mpf | Fraction
"""A coefficient of the recurrence."""

# The type of the mass beta_0: ``Number`` under HighPrecision, ``Number | None``
# under RATIONAL (None when the mass is irrational).
M = TypeVar("M", covariant=True)


class Irrational(ArithmeticError):
    """The rational arithmetic was asked for a quantity with no rational value (a transcendental mass)."""


class Arithmetic(Protocol[M]):
    """The number type a closed form is evaluated in, with the few transcendental functions a mass needs."""

    def working(self) -> ContextManager: ...
    def lift(self, x: int | float | Fraction) -> Number: ...
    def pi(self) -> Number: ...
    def sqrt(self, x: Number) -> Number: ...
    def gamma(self, x: Number) -> Number: ...
    def beta(self, a: Number, b: Number) -> Number: ...
    def power(self, base: int, exponent: Number) -> Number: ...
    def mass(self, thunk: Callable[[], Number]) -> M: ...


@dataclass(frozen=True)
class HighPrecision:
    """``mpmath.mpf`` at an explicit number of decimal digits."""

    digits: int

    def working(self) -> ContextManager:
        return mpmath.workdps(self.digits)

    def lift(self, x: int | float | Fraction) -> Number:
        if isinstance(x, Fraction):
            return mpmath.mpf(x.numerator) / x.denominator
        return mpmath.mpf(x)

    def pi(self) -> Number:
        return +mpmath.pi

    def sqrt(self, x: Number) -> Number:
        return mpmath.sqrt(x)

    def gamma(self, x: Number) -> Number:
        return mpmath.gamma(x)

    def beta(self, a: Number, b: Number) -> Number:
        return mpmath.beta(a, b)

    def power(self, base: int, exponent: Number) -> Number:
        return mpmath.power(base, exponent)

    def mass(self, thunk: Callable[[], Number]) -> Number:
        return thunk()


class _Rational:
    """Exact ``fractions.Fraction``; a transcendental quantity raises :class:`Irrational`."""

    def working(self) -> ContextManager:
        return contextlib.nullcontext()

    def lift(self, x: int | float | Fraction) -> Fraction:
        return Fraction(x)

    def pi(self) -> Fraction:
        raise Irrational("pi")

    def sqrt(self, x: Fraction) -> Fraction:
        # no shipped builder takes the root of a rational
        raise Irrational(f"sqrt({x})")

    def _positive_integer(self, x: Fraction) -> int:
        if x.denominator != 1 or x < 1:
            raise Irrational(f"gamma or beta of {x}")
        return int(x)

    def gamma(self, x: Fraction) -> Fraction:
        return Fraction(math.factorial(self._positive_integer(x) - 1))

    def beta(self, a: Fraction, b: Fraction) -> Fraction:
        ia, ib = self._positive_integer(a), self._positive_integer(b)
        return Fraction(math.factorial(ia - 1) * math.factorial(ib - 1), math.factorial(ia + ib - 1))

    def power(self, base: int, exponent: Fraction) -> Fraction:
        if exponent.denominator != 1:
            raise Irrational(f"{base}**{exponent}")
        return Fraction(base) ** int(exponent)

    def mass(self, thunk: Callable[[], Number]) -> Fraction | None:
        try:
            return thunk()
        except Irrational:
            return None


RATIONAL = _Rational()

# A recurrence-coefficient generator: (n, arithmetic) -> (alpha, beta),
# each a list of n numbers in the arithmetic's number type, with the storage
# convention beta[0] = mu_0. The arithmetic is an ARGUMENT, not mpmath's
# ambient setting, and the generator is a NOMINAL type: the only way to make
# one is ``over_arithmetic``, which evaluates the closed form inside the
# arithmetic's working context, so a closed form whose decorator was
# forgotten is a bare function and ``GeneratingMeasure.exact_recurrence``
# (annotated with this class) refuses it under pyright. The closed forms are
# written once; the float coefficients, every Gauss rule and the exact-tie
# decision are derived from them, so they cannot drift.


@dataclass(frozen=True)
class RecurrenceCoefficients:
    """A closed form ``(n, A) -> (alpha, beta)`` bound to evaluate inside ``A``'s working context."""

    closed_form: Callable[[int, Arithmetic[Any]], tuple[list[Number], list[Any]]]

    def __call__(self, n: int, arithmetic: Arithmetic[M]) -> tuple[list[Number], list[Number | M]]:
        with arithmetic.working():
            return self.closed_form(n, arithmetic)


def over_arithmetic(
    closed_form: Callable[[int, Arithmetic[Any]], tuple[list[Number], list[Any]]],
) -> RecurrenceCoefficients:
    """Lift a closed form ``(n, A) -> (alpha, beta)`` to a :class:`RecurrenceCoefficients`: evaluated inside the arithmetic's working context, never mpmath's ambient one."""
    return RecurrenceCoefficients(closed_form)


@over_arithmetic
def legendre_recurrence(n: int, A: Arithmetic[M]) -> tuple[list[Number], list[Number | M]]:
    beta: list[Number | M] = [A.mass(lambda: A.lift(2))]
    beta += [A.lift(k * k) / (4 * k * k - 1) for k in range(1, n)]
    return [A.lift(0)] * n, beta[:n]


@over_arithmetic
def chebyshev_t_recurrence(n: int, A: Arithmetic[M]) -> tuple[list[Number], list[Number | M]]:
    # beta_1 = 1/2 breaks the otherwise-constant 1/4 pattern; this is the
    # only irregular coefficient in the family.
    beta: list[Number | M] = [A.mass(A.pi), A.lift(1) / 2]
    beta += [A.lift(1) / 4] * max(n - 2, 0)
    return [A.lift(0)] * n, beta[:n]


@over_arithmetic
def chebyshev_u_recurrence(n: int, A: Arithmetic[M]) -> tuple[list[Number], list[Number | M]]:
    beta: list[Number | M] = [A.mass(lambda: A.pi() / 2)]
    beta += [A.lift(1) / 4] * (n - 1)
    return [A.lift(0)] * n, beta


@over_arithmetic
def hermite_recurrence(n: int, A: Arithmetic[M]) -> tuple[list[Number], list[Number | M]]:
    beta: list[Number | M] = [A.mass(lambda: A.sqrt(A.pi()))]
    beta += [A.lift(k) / 2 for k in range(1, n)]
    return [A.lift(0)] * n, beta


def jacobi_recurrence(a: float, b: float) -> RecurrenceCoefficients:
    r"""The Jacobi recurrence for weight :math:`(1-x)^a (1+x)^b` on :math:`[-1,1]` (Gautschi 2004 §1.3).

    ``a`` and ``b`` enter the working precision exactly (a float64 is a binary rational)."""

    @over_arithmetic
    def recurrence(n: int, A: Arithmetic[M]) -> tuple[list[Number], list[Number | M]]:
        Aa, Bb = A.lift(a), A.lift(b)
        AB = Aa + Bb
        # mu_0 = 2^(a+b+1) B(a+1, b+1).
        alpha = [(Bb - Aa) / (AB + 2)]
        beta: list[Number | M] = [A.mass(lambda: A.power(2, AB + 1) * A.beta(Aa + 1, Bb + 1))]
        if n > 1:
            # k = 1 is a REMOVABLE singularity of the general beta
            # formula below: its numerator carries (k + ab) and its
            # denominator (2k + ab - 1), and at k = 1 both are (1 + ab).
            # They cancel. Evaluating the general form there would divide
            # by zero whenever a + b = -1 — which is exactly Chebyshev-1,
            # jacobi(-1/2, -1/2), the single most common member.
            alpha.append((Bb * Bb - Aa * Aa) / ((2 + AB) * (4 + AB)))
            beta.append(4 * (1 + Aa) * (1 + Bb) / ((2 + AB) ** 2 * (3 + AB)))
        for k in range(2, n):
            two_k_ab = 2 * k + AB
            alpha.append((Bb * Bb - Aa * Aa) / (two_k_ab * (two_k_ab + 2)))
            beta.append(
                4 * k * (k + Aa) * (k + Bb) * (k + AB)
                / (two_k_ab**2 * (two_k_ab + 1) * (two_k_ab - 1))
            )
        return alpha[:n], beta[:n]

    return recurrence


def laguerre_recurrence(a: float) -> RecurrenceCoefficients:
    r"""The generalised Laguerre recurrence for weight :math:`x^a e^{-x}` on :math:`[0,\infty)` (Gautschi 2004 §1.3)."""

    @over_arithmetic
    def recurrence(n: int, A: Arithmetic[M]) -> tuple[list[Number], list[Number | M]]:
        Aa = A.lift(a)
        beta: list[Number | M] = [A.mass(lambda: A.gamma(Aa + 1))]
        beta += [k * (k + Aa) for k in range(1, n)]
        return [2 * k + Aa + 1 for k in range(n)], beta

    return recurrence


__all__ = [
    "Arithmetic",
    "HighPrecision",
    "Irrational",
    "Number",
    "RATIONAL",
    "RecurrenceCoefficients",
    "chebyshev_t_recurrence",
    "chebyshev_u_recurrence",
    "hermite_recurrence",
    "jacobi_recurrence",
    "laguerre_recurrence",
    "legendre_recurrence",
    "over_arithmetic",
]
