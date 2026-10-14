r"""mpmath references for the characteristic reference's Gauss rule of the attenuation weight (P3 step 3f, route B).

No import of the code under test. Every integral of a polynomial against an exponential is a CLOSED FORM by repeated
integration by parts, evaluated at a precision raised until two precisions 40 digits apart agree to 1e-30 (the
reference's own convergence control); the Legendre moments are also given through mpmath's Bessel function,
:math:`\int_{-1}^1 P_k(s) e^{-x(1-s)} ds = 2\tilde\imath_k(x)`, the identity :func:`legendre_toward_end` pins.
"""
from __future__ import annotations

import math

import mpmath as mp
import numpy as np


def _dps(x: float, degree: int) -> int:
    lost = math.lgamma(degree + 2) / math.log(10) + degree * math.log10(4.0)
    if 0.0 < x < 1.0:
        lost += (degree + 1) * -math.log10(x)
    return 40 + int(math.ceil(lost))


def _converged(compute, x: float, degree: int):
    d = _dps(x, degree)
    a = compute(d)
    for _ in range(40):
        b = compute(d + 40)
        scale = max(abs(b), mp.mpf("1e-300"))
        if abs(a - b) <= mp.mpf("1e-30") * scale:
            return b
        a, d = b, d + 40
    raise ArithmeticError(f"the closed-form reference did not converge by {d} digits (x = {x}, degree {degree})")


def _mul(p, q):
    out = [mp.mpf(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def _deriv(p):
    return [p[k] * k for k in range(1, len(p))] or [mp.mpf(0)]


def _val(p, s):
    return mp.polyval(list(reversed(p)), s)


def _antideriv(p):
    return [mp.mpf(0)] + [c / (k + 1) for k, c in enumerate(p)]


def _reflect(p):
    return [c * (-1) ** k for k, c in enumerate(p)]


def _legendre(k: int):
    p0, p1 = [mp.mpf(1)], [mp.mpf(0), mp.mpf(1)]
    if k == 0:
        return p0
    for m in range(1, k):
        a = _mul([mp.mpf(0), mp.mpf(2 * m + 1)], p1)
        b = p0 + [mp.mpf(0)] * (len(a) - len(p0))
        p0, p1 = p1, [(ai - m * bi) / (m + 1) for ai, bi in zip(a, b)]
    return p1


def _toward_end(p, X):
    r""":math:`\int_{-1}^1 p(s) e^{-X(1-s)} ds = \sum_j (-1)^j [p^{(j)}(1) - e^{-2X} p^{(j)}(-1)] / X^{j+1}`."""
    if X == 0:
        P = _antideriv(p)
        return _val(P, 1) - _val(P, -1)
    total, d, j = mp.mpf(0), p, 0
    e2 = mp.exp(-2 * X)
    while any(c != 0 for c in d):
        total += (-1) ** j * (_val(d, 1) - e2 * _val(d, -1)) / X ** (j + 1)
        d, j = _deriv(d), j + 1
    return total


def _toward_start(p, X):
    return _toward_end(_reflect(p), X)


def legendre_toward_end(k: int, x: float) -> float:
    r""":math:`\int_{-1}^1 P_k(s) e^{-x(1-s)} ds` by parts."""
    def compute(dps):
        with mp.workdps(dps):
            return _toward_end(_legendre(k), mp.mpf(repr(float(x))))

    return float(_converged(compute, x, k))


def legendre_moment(k: int, x: float) -> float:
    r""":math:`2\tilde\imath_k(x) = 2 e^{-x}\sqrt{\pi/(2x)}\,I_{k+1/2}(x)`, the same integral through the Bessel function."""
    if x == 0.0:
        return 2.0 if k == 0 else 0.0
    with mp.workdps(40):
        X = mp.mpf(repr(float(x)))
        return float(2 * mp.exp(-X) * mp.sqrt(mp.pi / (2 * X)) * mp.besseli(k + mp.mpf(1) / 2, X))


def triangle(f, g, x: float) -> float:
    r""":math:`\int_{-1}^1 f(s) \int_{-1}^s g(s') e^{-x(s-s')} ds' ds` for polynomials ``f``, ``g`` (ascending floats), by parts."""
    def compute(dps):
        with mp.workdps(dps):
            F = [mp.mpf(repr(float(c))) for c in f]
            G = [mp.mpf(repr(float(c))) for c in g]
            X = mp.mpf(repr(float(x)))
            if X == 0:
                inner = _antideriv(G)
                inner[0] -= _val(inner, -1)
                P = _antideriv(_mul(F, inner))
                return _val(P, 1) - _val(P, -1)
            A, C, d, j = [mp.mpf(0)] * len(G), mp.mpf(0), G, 0
            while any(c != 0 for c in d):
                for i, c in enumerate(d):
                    A[i] += (-1) ** j * c / X ** (j + 1)
                C += (-1) ** j * _val(d, -1) / X ** (j + 1)
                d, j = _deriv(d), j + 1
            P = _antideriv(_mul(F, A))
            return _val(P, 1) - _val(P, -1) - C * _toward_start(F, X)

    return float(_converged(compute, x, len(f) + len(g)))


def legendre_sums(offset: np.ndarray, weight: np.ndarray, count: int) -> list[float]:
    r""":math:`\sum_q w_q P_k(1 - \delta_q)` for :math:`k <` ``count``, the double offsets and weights promoted exactly, every
    :math:`P_k` and the sum formed in 40-digit mpmath by Bonnet's recurrence written in :math:`\delta` (so :math:`1 - \delta`
    is never formed and the sum measures the rule, not its evaluation)."""
    with mp.workdps(40):
        dm = [mp.mpf(repr(float(t))) for t in offset]
        wm = [mp.mpf(repr(float(t))) for t in weight]
        prev, cur = [mp.mpf(0)] * len(dm), [mp.mpf(1)] * len(dm)
        out = []
        for j in range(count):
            out.append(float(mp.fsum(w * p for w, p in zip(wm, cur))))
            prev, cur = cur, [((2 * j + 1) * (c - t * c) - j * q) / (j + 1) for c, t, q in zip(cur, dm, prev)]
        return out
