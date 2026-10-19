r"""The three-term recurrence of an orthogonal family: the Jacobi matrix.

The monic orthogonal polynomials :math:`\{\pi_k\}` of a positive measure
:math:`d\lambda` on the real line obey

.. math::

   \pi_{k+1}(x) = (x - \alpha_k)\,\pi_k(x) - \beta_k\,\pi_{k-1}(x),
   \qquad \pi_{-1} = 0,\ \pi_0 = 1,

and this module stores exactly those coefficients: each function returns
``(alpha, beta)``, two lists of ``n`` ``mpmath.mpf``, with the convention
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

Precision. Every coefficient is evaluated at an explicit number of decimal
digits (the precision is an argument, never mpmath's ambient setting). The
closed forms are written once, in arbitrary precision; every float
coefficient and every Gauss rule is derived from them.

References
----------
- Golub, G. H. & Welsch, J. H. (1969). Calculation of Gauss quadrature
  rules. *Math. Comp.* 23, 221-230.
- Gautschi, W. (2004). *Orthogonal Polynomials: Computation and
  Approximation*. Oxford. §1.3 (recurrence coefficients of the classical
  families) and §3.1.1 (the Jacobi matrix and Gauss rules).
"""

from __future__ import annotations

import functools
from typing import Callable

import mpmath

# Type of a recurrence-coefficient generator: (n, digits) -> (alpha, beta),
# each a list of n ``mpmath.mpf`` evaluated at ``digits`` decimal digits, with
# the storage convention beta[0] = mu_0. The precision is an ARGUMENT, not
# mpmath's ambient setting, so a bare call cannot return 15-digit values under
# the name "exact". The closed forms are written once, in arbitrary
# precision; the float coefficients and every Gauss rule are derived from
# them, so the two cannot drift.
RecurrenceCoefficients = Callable[[int, int], "tuple[list, list]"]


def _at_digits(closed_form: "Callable[[int], tuple[list, list]]") -> RecurrenceCoefficients:
    """Lift a closed form written in mpmath arithmetic to a
    :data:`RecurrenceCoefficients`: evaluate it at the precision the caller
    names, never at mpmath's ambient one."""

    @functools.wraps(closed_form)
    def recurrence(n: int, digits: int) -> "tuple[list, list]":
        with mpmath.workdps(digits):
            return closed_form(n)

    return recurrence


@_at_digits
def legendre_recurrence(n: int) -> "tuple[list, list]":
    alpha = [mpmath.mpf(0)] * n
    beta = [mpmath.mpf(2)] + [mpmath.mpf(k * k) / (4 * k * k - 1) for k in range(1, n)]
    return alpha, beta[:n]


@_at_digits
def chebyshev_t_recurrence(n: int) -> "tuple[list, list]":
    alpha = [mpmath.mpf(0)] * n
    # beta_1 = 1/2 breaks the otherwise-constant 1/4 pattern; this is the
    # only irregular coefficient in the family.
    beta = [mpmath.pi, mpmath.mpf(1) / 2] + [mpmath.mpf(1) / 4] * max(n - 2, 0)
    return alpha, beta[:n]


@_at_digits
def chebyshev_u_recurrence(n: int) -> "tuple[list, list]":
    alpha = [mpmath.mpf(0)] * n
    beta = [mpmath.pi / 2] + [mpmath.mpf(1) / 4] * (n - 1)
    return alpha, beta


@_at_digits
def hermite_recurrence(n: int) -> "tuple[list, list]":
    alpha = [mpmath.mpf(0)] * n
    beta = [mpmath.sqrt(mpmath.pi)] + [mpmath.mpf(k) / 2 for k in range(1, n)]
    return alpha, beta


def jacobi_recurrence(a: float, b: float) -> RecurrenceCoefficients:
    r"""The Jacobi recurrence for weight :math:`(1-x)^a (1+x)^b` on :math:`[-1,1]` (Gautschi 2004 §1.3).

    ``a`` and ``b`` enter the working precision exactly (a float64 is a binary rational)."""

    @_at_digits
    def recurrence(n: int) -> "tuple[list, list]":
        A, B = mpmath.mpf(a), mpmath.mpf(b)
        AB = A + B
        # mu_0 = 2^(a+b+1) B(a+1, b+1).
        alpha = [(B - A) / (AB + 2)]
        beta = [mpmath.power(2, AB + 1) * mpmath.beta(A + 1, B + 1)]
        if n > 1:
            # k = 1 is a REMOVABLE singularity of the general beta
            # formula below: its numerator carries (k + ab) and its
            # denominator (2k + ab - 1), and at k = 1 both are (1 + ab).
            # They cancel. Evaluating the general form there would divide
            # by zero whenever a + b = -1 — which is exactly Chebyshev-1,
            # jacobi(-1/2, -1/2), the single most common member.
            alpha.append((B * B - A * A) / ((2 + AB) * (4 + AB)))
            beta.append(4 * (1 + A) * (1 + B) / ((2 + AB) ** 2 * (3 + AB)))
        for k in range(2, n):
            two_k_ab = 2 * k + AB
            alpha.append((B * B - A * A) / (two_k_ab * (two_k_ab + 2)))
            beta.append(
                4 * k * (k + A) * (k + B) * (k + AB)
                / (two_k_ab**2 * (two_k_ab + 1) * (two_k_ab - 1))
            )
        return alpha[:n], beta[:n]

    return recurrence


def laguerre_recurrence(a: float) -> RecurrenceCoefficients:
    r"""The generalised Laguerre recurrence for weight :math:`x^a e^{-x}` on :math:`[0,\infty)` (Gautschi 2004 §1.3)."""

    @_at_digits
    def recurrence(n: int) -> "tuple[list, list]":
        A = mpmath.mpf(a)
        alpha = [2 * k + A + 1 for k in range(n)]
        beta = [mpmath.gamma(A + 1)] + [k * (k + A) for k in range(1, n)]
        return alpha, beta

    return recurrence


__all__ = [
    "RecurrenceCoefficients",
    "chebyshev_t_recurrence",
    "chebyshev_u_recurrence",
    "hermite_recurrence",
    "jacobi_recurrence",
    "laguerre_recurrence",
    "legendre_recurrence",
]
