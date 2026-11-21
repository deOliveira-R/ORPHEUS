"""Polynomials over the rational field: the arithmetic of the exact-tie decision.

A polynomial is a tuple of ``fractions.Fraction`` coefficients, lowest degree first. A value on a float64 rounding tie is an exact
rational, and whether a Gauss node or weight EQUALS it is a question about the roots of a polynomial over Q (the monic orthogonal
polynomial :math:`\\pi_n`) and about a greatest common divisor in the Euclidean ring :math:`\\mathbb{Q}[x]`; this module is that
arithmetic and nothing else (no space, field or operator). :func:`trim` is the parse door: every entry becomes a ``Fraction`` (an
``int`` or a float64 is an exact rational), so every function returns ``Fraction`` entries whatever it is given.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from fractions import Fraction

Polynomial = tuple[Fraction, ...]
ZERO: Polynomial = (Fraction(0),)
ONE: Polynomial = (Fraction(1),)


def trim(p: Iterable[Fraction | int | float]) -> Polynomial:
    """The parse door: every entry becomes a ``Fraction`` (a float64 is an exact binary rational) and the leading zeros go."""
    q = [Fraction(x) for x in p]
    while len(q) > 1 and q[-1] == 0:
        q.pop()
    return tuple(q) if q else ZERO


def add(p: Sequence[Fraction], q: Sequence[Fraction]) -> Polynomial:
    m = max(len(p), len(q))
    return trim([(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(m)])


def scale(p: Sequence[Fraction], c: Fraction | int) -> Polynomial:
    return trim([c * x for x in p])


def multiply(p: Sequence[Fraction], q: Sequence[Fraction]) -> Polynomial:
    p, q = trim(p), trim(q)
    out = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        if x:
            for j, y in enumerate(q):
                out[i + j] += x * y
    return trim(out)


def remainder(p: Sequence[Fraction], q: Sequence[Fraction]) -> Polynomial:
    """The remainder of ``p`` by the nonzero ``q``, by long division over Q."""
    r, d = list(trim(p)), trim(q)
    if d == ZERO:
        raise ZeroDivisionError("polynomial division by zero")
    while len(r) >= len(d) and r != [0]:
        f = r[-1] / d[-1]
        shift = len(r) - len(d)
        for i, y in enumerate(d):
            r[shift + i] -= f * y
        r = list(trim(r[:-1])) if len(r) > 1 else [Fraction(0)]
    return trim(r)


def gcd(p: Sequence[Fraction], q: Sequence[Fraction]) -> Polynomial:
    """The monic greatest common divisor in Q[x] (Euclid)."""
    a, b = trim(p), trim(q)
    while b != ZERO:
        a, b = b, remainder(a, b)
    return scale(a, 1 / a[-1]) if a != ZERO else ZERO


def evaluate(p: Sequence[Fraction], x: Fraction) -> Fraction:
    acc = Fraction(0)
    for c in reversed(p):
        acc = acc * x + c
    return acc


__all__ = ["ONE", "ZERO", "Polynomial", "add", "evaluate", "gcd", "multiply", "remainder", "scale", "trim"]
