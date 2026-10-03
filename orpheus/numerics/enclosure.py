r"""An enclosure of a real number: a value and a guaranteed distance to the exact one (#405 P2 step 1).

``Enclosure(value, bound)`` claims that an exact real number lies in the
closed interval :math:`[v - b, v + b]`, taken in exact real arithmetic on the
two doubles. It is the form of a reference's claim about the exact answer of
its own equation (the user's ruling of 2026-10-02, "G3"): a guarantee, never
an estimate. An exact value is the enclosure with bound 0. A statistical
interval is not an enclosure: it is a production method's self-report, and
verification puts it on trial (``.claude/plans/reference_cache.md``, "P2
rulings so far").

**The quotient** of two enclosures encloses every quotient of their members:
:math:`x/y` for :math:`x \in [a]`, :math:`y \in [b]`. On a box that excludes
:math:`y = 0` the quotient is monotone in each argument, so its extremes are
among the four corner quotients. Every floating-point step is rounded
outward by one ``nextafter`` (an endpoint, a corner, the half-width), which
covers the half-ulp error of the rounded operation that produced it; the
bound is therefore valid, and loose by at most a few ulp (the gate R1.4 of
``.claude/plans/reference_p2_spec.md`` derives 10.5). The centre is the
quotient of the centres, so a ratio observable reads the ratio of its
readings. A denominator whose enclosure holds zero has no quotient, and an
overflow is refused rather than returned as an infinite bound. Both refusals
read the ends after the outward step, so they are conservative by one
rounding step: a denominator one subnormal from zero, or a quotient within
one step of the largest double, is refused although its exact quotient
exists.

Only an enclosure divides an enclosure: an exact number is written
``Enclosure(v, 0)`` once, at the place it enters, so that no bare float can
pass for a claim.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from orpheus.numerics.content import ContentIdentity, content_digest
from orpheus.numerics.scalars import parse_finite_real

__all__ = ["Enclosure"]


def _down(x: float) -> float:
    """The next double below ``x``: one outward step under a rounded operation."""
    return math.nextafter(x, -math.inf)


def _up(x: float) -> float:
    """The next double above ``x``."""
    return math.nextafter(x, math.inf)


@dataclass(frozen=True, eq=False)
class Enclosure(ContentIdentity):
    """The exact number lies within ``bound`` of ``value``."""

    value: float
    bound: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", parse_finite_real(self.value, "Enclosure: the value"))
        bound = parse_finite_real(self.bound, "Enclosure: the bound")
        if bound < 0.0:
            raise ValueError(f"Enclosure: the bound is a distance, non-negative, got {bound!r}")
        object.__setattr__(self, "bound", bound)
        content_digest(self)

    def ends(self) -> tuple[float, float]:
        """The interval's ends as doubles, each rounded outward: a superset of the claimed interval."""
        return _down(self.value - self.bound), _up(self.value + self.bound)

    def __truediv__(self, other: Enclosure) -> Enclosure:
        if not isinstance(other, Enclosure):
            return NotImplemented
        low_b, high_b = other.ends()
        if low_b <= 0.0 <= high_b:
            raise ZeroDivisionError(
                f"Enclosure: the denominator's enclosure, rounded outward to [{low_b!r}, {high_b!r}], contains zero, "
                f"so the quotient has no enclosure"
            )
        low_a, high_a = self.ends()
        corners = [x / y for x in (low_a, high_a) for y in (low_b, high_b)]
        low, high = _down(min(corners)), _up(max(corners))
        centre = self.value / other.value
        half_width = max(_up(high - centre), _up(centre - low))
        if not math.isfinite(half_width):
            raise ValueError(
                f"Enclosure: the quotient {self!r} / {other!r} overflows a double once rounded outward: its bound would be infinite"
            )
        return Enclosure(centre, half_width)
