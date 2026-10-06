"""Arbitrary-precision references for the characteristic reference's basis and line transport.

Shared by ``test_characteristic_basis.py`` and ``test_characteristic_transport.py``
(P1 step (b), second rung). NOTHING here imports the code under test: the line
geometry (crossing parameters, segment regions, the closest approach), the
Gauss-Legendre nodes, the Lagrange functions, the volume densities and every
line integral are written in mpmath from the closed forms, so a reference value
never descends from a production spelling (``instrument-doctrine`` X4).

The only inputs taken from the code under test are DATA the gates pass in
explicitly: a basis's panel ends (gated on their own by the partition rows) and,
for a source written as basis coefficients, the nodes at which a per-region
polynomial is sampled (a per-region polynomial of degree <= p is reproduced by
the panel interpolant whatever the nodes, so the reference integral of the
polynomial does not depend on them).

Conventions, the ones the rung's sketch states:

* a line is ``{foot + t Omega}``, ``t`` the 3-D arc length from the foot (the
  point of the line nearest the origin, ``Line.foot``);
* the orbit coordinate is ``|P x|`` on a radial chart (the kept columns: 2 on a
  cylinder, 3 on a sphere) and ``x_0`` on the slab;
* a transit is a maximal run of the line inside ``[r_0, r_n]``; a traversal is a
  transit read forward (increasing ``t``) or reversed;
* every integral is over 3-D arc length; no ``1/4 pi`` and no ``1/|Omega_x|``
  beyond the arc length's own.
"""
from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any
from dataclasses import dataclass
from functools import lru_cache

import mpmath as mp
import numpy as np
from mpmath.calculus.quadrature import GaussLegendre

DPS = 30

#: An mpmath real (mpmath ships no type stubs, so its value type is spelled as Any for the checker).
Mpf = Any

#: A source along the line: f(region, orbit coordinate) as an mpf.
Source = Callable[[int, "Mpf"], "Mpf"]

_KEPT = {"slab": 1, "cylinder": 2, "sphere": 3}


def mpf(x) -> Mpf:
    return mp.mpf(float(x)) if not isinstance(x, mp.mpf) else x


# ── volume densities, written by hand per chart (never from measure_constant) ──


def density(chart: str, c: Mpf) -> Mpf:
    """dV/dc per unit of the discarded columns: 1 (slab), 2 pi c (cylinder per unit height), 4 pi c^2 (sphere)."""
    match chart:
        case "slab":
            return mp.mpf(1)
        case "cylinder":
            return 2 * mp.pi * c
        case "sphere":
            return 4 * mp.pi * c * c
    raise AssertionError(chart)


def volume(chart: str, a, b) -> Mpf:
    """The measure of the shell [a, b], closed form: b - a, pi (b^2 - a^2), 4/3 pi (b^3 - a^3)."""
    a, b = mpf(a), mpf(b)
    match chart:
        case "slab":
            return b - a
        case "cylinder":
            return mp.pi * (b * b - a * a)
        case "sphere":
            return 4 * mp.pi * (b ** 3 - a ** 3) / 3
    raise AssertionError(chart)


def quad(f: Callable[[Any], Any], interval: Sequence[Any]) -> Mpf:
    """mp.quad, typed: the value only (mpmath returns a (value, error) pair only when asked)."""
    return mp.quad(f, interval)


# ── Gauss-Legendre nodes and Lagrange functions, in mpmath ─────────────────


@lru_cache(maxsize=None)
def _legendre_roots(n: int) -> tuple[Mpf, ...]:
    with mp.workdps(DPS + 10):
        guesses = np.polynomial.legendre.leggauss(n)[0]
        return tuple(mp.findroot(lambda x: mp.legendre(n, x), mp.mpf(float(g))) for g in guesses)


def gl_nodes(a, b, n: int) -> list[Mpf]:
    """The n Gauss-Legendre points of [a, b], increasing, from the roots of P_n (mpmath)."""
    a, b = mpf(a), mpf(b)
    return [a + (b - a) * (x + 1) / 2 for x in _legendre_roots(n)]


def lagrange(nodes: Sequence[Mpf], i: int, c: Mpf) -> Mpf:
    """The i-th Lagrange cardinal function of ``nodes`` at c, by the product formula."""
    value = mp.mpf(1)
    for m, x in enumerate(nodes):
        if m != i:
            value *= (c - x) / (nodes[i] - x)
    return value


# ── the line, in mpmath ────────────────────────────────────────────────────


@dataclass
class Segment:
    t_a: Mpf
    t_b: Mpf
    region: int            # 0..n-1 inside, n the inner exterior (cavity), n + 1 the outer exterior


@dataclass
class Traversal:
    t_in: Mpf           # the smaller parameter of the transit
    t_out: Mpf          # the larger
    reversed: bool
    segments: list[Segment]

    @property
    def entry(self) -> Mpf:
        return self.t_out if self.reversed else self.t_in

    @property
    def exit(self) -> Mpf:
        return self.t_in if self.reversed else self.t_out


class MpLine:
    """One line through a concentric body, every geometric quantity in mpmath.

    ``extra_radii`` are further level sets at which integrals are split (a
    basis's panel ends): they never change a segment's region.
    """

    def __init__(self, chart: str, breakpoints: Sequence[float], point, direction,
                 extra_radii: Sequence[float] = ()) -> None:
        self.chart = chart
        self.r = [mpf(v) for v in breakpoints]
        self.n = len(self.r) - 1
        with mp.workdps(DPS):
            p = [mpf(v) for v in point]
            w = [mpf(v) for v in direction]
            norm = mp.sqrt(sum(v * v for v in w))
            w = [v / norm for v in w]
            pw = sum(a * b for a, b in zip(p, w))
            self.foot = [a - pw * b for a, b in zip(p, w)]
            self.omega = w
            d = _KEPT[chart]
            self.d = d
            if chart == "slab":
                self.rate = w[0]
                self.t_star = None
                self.b = None
            else:
                fk, uk = self.foot[:d], w[:d]
                speed2 = sum(v * v for v in uk)
                self.speed = mp.sqrt(speed2)
                self.t_star = -sum(a * b for a, b in zip(fk, uk)) / speed2
                closest = [a + self.t_star * b for a, b in zip(fk, uk)]
                self.b = mp.sqrt(sum(v * v for v in closest))
            self._split = sorted(set(self._crossings(self.r) + self._crossings([mpf(v) for v in extra_radii])))

    def parameter_of(self, point) -> Mpf:
        return sum((mpf(x) - f) * w for x, f, w in zip(point, self.foot, self.omega))

    def c(self, t: Mpf) -> Mpf:
        """The orbit coordinate at parameter t."""
        x = [f + t * w for f, w in zip(self.foot, self.omega)]
        if self.chart == "slab":
            return x[0]
        return mp.sqrt(sum(v * v for v in x[: self.d]))

    def _crossings(self, radii: Sequence[Mpf]) -> list[Mpf]:
        ts = []
        for r in radii:
            if self.chart == "slab":
                ts.append((r - self.foot[0]) / self.rate)
            elif r > self.b:
                h = mp.sqrt((r - self.b) * (r + self.b)) / self.speed
                ts += [self.t_star - h, self.t_star + h]
        return ts

    def region_of(self, c: Mpf) -> int:
        if c < self.r[0]:
            return self.n
        if c > self.r[-1]:
            return self.n + 1
        for k in range(self.n):
            if c <= self.r[k + 1]:
                return k
        return self.n + 1

    def segments(self) -> list[Segment]:
        """The stretches between consecutive crossings (and the closest approach), each with its region."""
        with mp.workdps(DPS):
            cuts = list(self._split)
            if self.t_star is not None and cuts and cuts[0] < self.t_star < cuts[-1]:
                cuts = sorted(set(cuts + [self.t_star]))
            return [Segment(a, b, self.region_of(self.c((a + b) / 2))) for a, b in zip(cuts[:-1], cuts[1:])]

    def transits(self) -> list[tuple[Mpf, Mpf, list[Segment]]]:
        """The maximal runs inside [r_0, r_n], in order along the line."""
        out, run = [], []
        for seg in self.segments():
            if seg.region < self.n:
                run.append(seg)
            elif run:
                out.append(run)
                run = []
        if run:
            out.append(run)
        return [(run[0].t_a, run[-1].t_b, run) for run in out]

    def traversal(self, transit: int, reversed_: bool) -> Traversal:
        t_in, t_out, segs = self.transits()[transit]
        return Traversal(t_in, t_out, reversed_, segs)


def optical_depth(trav: Traversal, sigma: Sequence[float], a: Mpf, b: Mpf) -> Mpf:
    """The optical depth between parameters a and b inside the traversal (order-free)."""
    lo, hi = (a, b) if a <= b else (b, a)
    total = mp.mpf(0)
    for s in trav.segments:
        overlap = min(hi, s.t_b) - max(lo, s.t_a)
        if overlap > 0:
            total += mpf(sigma[s.region]) * overlap
    return total


def _pieces(trav: Traversal, lo: Mpf, hi: Mpf) -> list[tuple[Mpf, Mpf, int]]:
    out = []
    for s in trav.segments:
        a, b = max(lo, s.t_a), min(hi, s.t_b)
        if b > a:
            out.append((a, b, s.region))
    return out


def outflow(line: MpLine, trav: Traversal, f: Source, sigma: Sequence[float]) -> Mpf:
    """B: the source integral of the traversal attenuated to its EXIT."""
    with mp.workdps(DPS):
        total = mp.mpf(0)
        for a, b, k in _pieces(trav, trav.t_in, trav.t_out):
            total += mp.quad(lambda t, k=k: f(k, line.c(t)) * mp.exp(-optical_depth(trav, sigma, t, trav.exit)), [a, b])
        return total


def entry_response(line: MpLine, trav: Traversal, f: Source, sigma: Sequence[float]) -> Mpf:
    """A: the integral of f attenuated from the traversal's ENTRY."""
    with mp.workdps(DPS):
        total = mp.mpf(0)
        for a, b, k in _pieces(trav, trav.t_in, trav.t_out):
            total += mp.quad(lambda t, k=k: f(k, line.c(t)) * mp.exp(-optical_depth(trav, sigma, trav.entry, t)), [a, b])
        return total


def partial_outflow(line: MpLine, trav: Traversal, f: Source, sigma: Sequence[float], t: Mpf) -> Mpf:
    """The source integral from the traversal's entry to the parameter t, attenuated to t."""
    with mp.workdps(DPS):
        lo, hi = (t, trav.t_out) if trav.reversed else (trav.t_in, t)
        total = mp.mpf(0)
        for a, b, k in _pieces(trav, lo, hi):
            total += mp.quad(lambda s, k=k: f(k, line.c(s)) * mp.exp(-optical_depth(trav, sigma, s, t)), [a, b])
        return total


def _gl_rule(degree: int):
    rule = GaussLegendre(mp.mp)
    return rule.calc_nodes(degree, mp.mp.prec)       # [(x, w)] on [-1, 1]


def volterra_bilinear(line: MpLine, trav: Traversal, f: Source, g: Source, sigma: Sequence[float],
                      degree: int = 5) -> Mpf:
    r"""x^T V y = \int f(s) \int_{entry}^{s} g(s') e^{-tau(s' -> s)} ds' ds along the traversal.

    Written through H(s) = \int_{entry}^{s} g e^{+tau(entry -> s')} ds', so that
    the inner integral is one cumulative pass: the outer rule is mpmath's
    Gauss-Legendre of the given degree on each smooth piece (in travel order),
    and H between consecutive outer nodes is mp.quad (tanh-sinh). The caller
    compares two degrees as the reference's own convergence control.
    """
    with mp.workdps(DPS):
        pieces = _pieces(trav, trav.t_in, trav.t_out)
        if trav.reversed:
            pieces = [(b, a, k) for a, b, k in reversed(pieces)]        # travel order: from entry
        rule = _gl_rule(degree)
        total, H, s_prev, k_prev = mp.mpf(0), mp.mpf(0), trav.entry, None
        for a, b, k in pieces:
            half, mid = (b - a) / 2, (a + b) / 2
            nodes = sorted(((mid + half * x, abs(half) * w) for x, w in rule), key=lambda p: abs(p[0] - a))
            for s, w in nodes:
                lo, hi = (s_prev, s) if s >= s_prev else (s, s_prev)
                # the inner integrand on [s_prev, s], split at the segment ends it straddles
                cuts = [lo] + [e for seg in trav.segments for e in (seg.t_a, seg.t_b) if lo < e < hi] + [hi]
                inc = mp.mpf(0)
                for u, v in zip(cuts[:-1], cuts[1:]):
                    kk = line.region_of(line.c((u + v) / 2))
                    inc += mp.quad(lambda x, kk=kk: g(kk, line.c(x)) * mp.exp(optical_depth(trav, sigma, trav.entry, x)), [u, v])
                H += inc
                s_prev = s
                total += w * f(k, line.c(s)) * mp.exp(-optical_depth(trav, sigma, trav.entry, s)) * H
            # the rest of the piece, up to its far end, so the next piece starts from H at its start
            lo, hi = (s_prev, b) if b >= s_prev else (b, s_prev)
            H += mp.quad(lambda x, k=k: g(k, line.c(x)) * mp.exp(optical_depth(trav, sigma, trav.entry, x)), [lo, hi])
            s_prev, k_prev = b, k
        return total


# ── per-region polynomials (the basis's exactness family) ───────────────────


def region_polynomial(coefficients: Sequence[Sequence[float]]) -> Source:
    """f(k, c) = sum_m coefficients[k][m] c^m, one polynomial per region."""
    coeffs = [[mpf(a) for a in row] for row in coefficients]

    def f(k: int, c: Mpf) -> Mpf:
        if k >= len(coeffs):
            return mp.mpf(0)
        return mp.polyval(list(reversed(coeffs[k])), c)

    return f


def sample(f: Source, regions: np.ndarray, nodes: np.ndarray) -> np.ndarray:
    """The coefficient vector of f on a nodal basis: f at each node, in the node's region (float)."""
    return np.array([float(f(int(k), mpf(x))) for k, x in zip(regions, nodes)])


def panel_function(nodes: Sequence[Mpf], i: int, a: Mpf, b: Mpf) -> Source:
    """The i-th Lagrange function of a panel [a, b], zero outside it (region-free)."""

    def f(_k: int, c: Mpf) -> Mpf:
        return lagrange(nodes, i, c) if a <= c <= b else mp.mpf(0)

    return f
