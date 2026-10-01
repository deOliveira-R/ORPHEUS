r"""Rigorous floating-point error bounds for gate design (proposal: ``tests/_harness/float_bounds.py``).

Two instruments, both from the standard model of IEEE-754 binary64 arithmetic
(Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., §2.2,
eq. (2.4)): for every basic operation ``op`` on floats ``x, y`` with no
underflow or overflow, ``fl(x op y) = (x op y)(1 + δ)``, ``|δ| ≤ u = 2⁻⁵³``.

* :func:`gamma` — Higham's :math:`\gamma_k = k u / (1 - k u)` (Lemma 3.1),
  the bound on :math:`|\theta_k|` where :math:`\prod_{i=1}^k (1+\delta_i)^{\pm 1}
  = 1 + \theta_k`. The order-free class bounds every redesigned gate uses are
  instances of it:

  - a sum of ``n`` terms in ANY association (recursive, pairwise, blocked,
    SIMD-lane): each term passes through at most ``n − 1`` additions, so
    :math:`|\hat s - s| \le \gamma_{n-1} \sum |x_i|` (Higham eq. (4.4) and the
    remark after it: the bound holds for every ordering);
  - a product of ``f`` factors in any order: relative error ``γ_{f−1}``;
  - a multilinear sum (n terms, each an f-factor product) in any bracketing:
    :math:`|\hat s - s| \le \gamma_{(f-1)+(n-1)} \sum |\text{term}_i|`.

* :class:`Tracked` — running error analysis (Higham §3.3) carried in exact
  rational arithmetic: a value is the pair (exact rational ``v``, radius ``r``)
  meaning "the float the program holds is within ``r`` of ``v``". Each
  operation propagates the radius rigorously and adds the operation's own
  rounding ``u·|result|``. Used where a gate needs the bound of ONE stated
  evaluation order (the test's own oracle, or a production loop whose order is
  part of its contract), never as a substitute for the order-free bounds.

Everything is exact (``fractions.Fraction``); no bound is evaluated in floating
point, so the bound cannot itself be the victim of the rounding it bounds.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

#: Unit roundoff of binary64, round-to-nearest.
U = Fraction(1, 2**53)


def gamma(k: int) -> Fraction:
    r""":math:`\gamma_k = k u / (1 - k u)` (Higham Lemma 3.1), exactly."""
    if k < 0:
        raise ValueError(f"gamma needs k >= 0, got {k}")
    return k * U / (1 - k * U)


def ulp(x: float) -> Fraction:
    """The spacing of binary64 at ``x`` (``math.ulp``), exactly."""
    import math

    return Fraction(math.ulp(float(x)))


@dataclass(frozen=True)
class Tracked:
    """An exact value and a rigorous radius on the float that represents it."""

    v: Fraction
    r: Fraction = Fraction(0)

    @staticmethod
    def of(x: float) -> "Tracked":
        """A float input: exact, radius 0."""
        return Tracked(Fraction(float(x)), Fraction(0))

    @staticmethod
    def const(x: int | Fraction) -> "Tracked":
        return Tracked(Fraction(x), Fraction(0))

    def _round(self) -> "Tracked":
        # The operation's own rounding: |fl(z) - z| <= u |z| <= u (|v| + r).
        return Tracked(self.v, self.r + U * (abs(self.v) + self.r))

    def __add__(self, other: "Tracked") -> "Tracked":
        return Tracked(self.v + other.v, self.r + other.r)._round()

    def __sub__(self, other: "Tracked") -> "Tracked":
        return Tracked(self.v - other.v, self.r + other.r)._round()

    def __mul__(self, other: "Tracked") -> "Tracked":
        r = abs(self.v) * other.r + abs(other.v) * self.r + self.r * other.r
        return Tracked(self.v * other.v, r)._round()

    def __truediv__(self, other: "Tracked") -> "Tracked":
        if other.r >= abs(other.v):
            raise ZeroDivisionError("divisor's radius reaches zero: no bound")
        q = self.v / other.v
        r = (self.r + abs(q) * other.r) / (abs(other.v) - other.r)
        return Tracked(q, r)._round()

    def scaled(self, power_of_two: Fraction) -> "Tracked":
        """Multiplication by a power of two is exact (no underflow): no rounding."""
        return Tracked(self.v * power_of_two, self.r * abs(power_of_two))


__all__ = ["U", "Tracked", "gamma", "ulp"]


def mp_fraction(x) -> Fraction:
    """An ``mpmath.mpf`` as the exact rational it holds (mantissa · 2^exponent)."""
    import mpmath

    # An mpf is read as it stands: ``mpf(x)`` would re-round it to the
    # CALLER's working precision, which outside a ``workdps`` block is 53 bits.
    if not isinstance(x, mpmath.mpf):
        x = mpmath.mpf(x)
    if x == 0:
        return Fraction(0)
    sign, man, exp, _ = x._mpf_
    value = Fraction(int(man)) * Fraction(2) ** int(exp)
    return -value if sign else value
