"""Arbitrary-precision references for the characteristic reference's reading at a point (P1 step (b), rung 5b).

Used by ``test_characteristic_reading.py``. NOTHING here imports the code under test: the backward path from the
point, its passages through the regions, the reflections at the walls, the angular integrals and the escape and
transmission probabilities are written in mpmath from the closed forms (``instrument-doctrine`` X4). The inputs are
data: the breakpoints, each region's total cross section and each region's emission RATE ``q`` (per unit volume,
integrated over the directions, the code's convention), and the walls' specular or diffuse amplitudes.

Conventions. The scalar flux is ``phi(x) = int psi dOmega`` with ``psi = (1/4 pi) int q e^{-tau} ds`` along the
backward path, so a closed homogeneous body reads ``phi = q / Sigma_t`` (the code's normalisation: a line carries
``4 pi`` times a flux per steradian and its weight holds ``1/4 pi``). A wall of specular amplitude ``a`` multiplies
the path that reflects at it by ``a``; on a hollow body the inner wall is a wall like the outer one (a vacuum inner
wall absorbs, the interpretation (A) of the P1 spec).

The four routes, one per chart and law, each a one-dimensional integral over the direction:

* sphere, specular laws: the backward path from the point lies in the plane of the point and the direction, a disk
  billiard of impact parameter ``b = r sqrt(1 - mu^2)``; the first leg, then the period of transits between the walls
  (one chord, or outer-to-inner and back through a cavity), closed by the geometric sum of the period
  (``psi_sphere``); ``phi = (1/2) int psi dmu`` split at every tangency;
* sphere centre: every direction is a diameter, the elementary closed form (``phi_sphere_centre``);
* cylinder axis: the polar integral of the diameter's closure (``phi_cylinder_axis``), and under vacuum the
  closed form in Bickley's ``Ki_2`` (``phi_cylinder_axis_vacuum``); a general point under vacuum, the azimuthal
  integral of the in-plane ``Ki_2`` differences (``phi_cylinder_vacuum``);
* slab: the ``E_2`` image series of the unfolded path, any albedos (``phi_slab``);
* white walls on a homogeneous sphere or slab: the first flight plus the re-entering partial currents of the
  escape and transmission probabilities (``phi_white_sphere``, ``phi_white_slab``).

These routes were validated, before the reading existed, against the old trajectory-resolvent family at its own
floor (``scratch/characteristic_architecture/p1_gates/probe_closed_forms.py``: the centre to 1.8e-16, the
general-point route at ``r = 1e-12`` against the centre's closed form to 1.2e-21, the axis's two routes to 1.2e-21,
an artificial interface in the slab series to 3.5e-26).
"""
from __future__ import annotations

from collections.abc import Sequence
from functools import lru_cache
from typing import Any

import mpmath as mp

from tests.gates.derivations._characteristic_mp import bickley_ki1

DPS = 25

Mpf = Any

#: mp.quad, typed: the value only (its stubs admit a (value, error) pair, returned only when asked).
_quad: Any = mp.quad


def _mp(values: Sequence[float]) -> list[Mpf]:
    return [mp.mpf(repr(float(v))) if not isinstance(v, mp.mpf) else v for v in values]


def _region_of(radius: Mpf, edges: Sequence[Mpf]) -> int:
    """The region holding an orbit coordinate strictly inside the body (a midpoint, never a breakpoint)."""
    for j in range(len(edges) - 1):
        if radius <= edges[j + 1]:
            return j
    return len(edges) - 2


# ── the disk billiard: the backward path in the plane of the point and the direction ──


def _passages(edges: Sequence[Mpf], b: Mpf, s_from: Mpf, s_to: Mpf) -> list[tuple[Mpf, int]]:
    """The passages (in-plane length, region) of the straight path at impact ``b`` from position ``s_from`` to ``s_to``.

    Positions ``s`` are signed along the path from its closest approach, ``rho^2 = b^2 + s^2``. The cuts are the
    crossings ``+-sqrt(r_k^2 - b^2)`` of every interface (and the closest approach) strictly inside the stretch.
    """
    lo, hi = (s_from, s_to) if s_from <= s_to else (s_to, s_from)
    cuts = {lo, hi}
    for r in edges:
        if r > b:
            h = mp.sqrt((r - b) * (r + b))
            for s in (-h, h):
                if lo < s < hi:
                    cuts.add(s)
    if lo < 0 < hi:
        cuts.add(mp.mpf(0))
    cuts = sorted(cuts)
    out = [(t - u, _region_of(mp.sqrt(b * b + ((u + t) / 2) ** 2), edges)) for u, t in zip(cuts[:-1], cuts[1:])]
    return out if s_from <= s_to else out[::-1]


def _leg(passages: list[tuple[Mpf, int]], sigma: Sequence[Mpf], q: Sequence[Mpf], stretch: Mpf,
         tau0: Mpf) -> tuple[Mpf, Mpf]:
    """``(int q e^{-tau} ds, optical depth)`` of a leg starting at optical depth ``tau0``; lengths times ``stretch``."""
    total, tau = mp.mpf(0), tau0
    for length, j in passages:
        ell = length * stretch
        if sigma[j] == 0:
            total += q[j] * ell * mp.exp(-tau)
            continue
        total += q[j] / sigma[j] * mp.exp(-tau) * (-mp.expm1(-sigma[j] * ell))
        tau += sigma[j] * ell
    return total, tau - tau0


def _transits(edges: Sequence[Mpf], b: Mpf, wall: str) -> list[tuple[list[tuple[Mpf, int]], str]]:
    """The period of transits leaving ``wall`` ('outer' or 'inner') at impact ``b``: [(passages, the wall it ends at)]."""
    r0, R = edges[0], edges[-1]
    hR = mp.sqrt((R - b) * (R + b))
    if not (r0 > 0 and b < r0):
        return [(_passages(edges, b, -hR, hR), "outer")]
    h0 = mp.sqrt((r0 - b) * (r0 + b))
    out_in = (_passages(edges, b, -hR, -h0), "inner")
    in_out = (_passages(edges, b, h0, hR), "outer")
    return [out_in, in_out] if wall == "outer" else [in_out, out_in]


def psi_disk(edges, sigma, q, amplitude: dict[str, Mpf], r: Mpf, mu: Mpf, stretch: Mpf = mp.mpf(1)) -> Mpf:
    r"""``int q e^{-tau} ds`` along the unfolded backward path from radius ``r`` with direction cosine ``mu``.

    ``mu`` is the cosine between the direction of flight and the outward radial direction at the point;
    ``amplitude`` the specular amplitude of the ``'outer'`` and ``'inner'`` walls; ``stretch`` the 3-D length per
    in-plane length (1 on a sphere, ``1/sin(theta)`` on a cylinder). The first leg runs to the first wall; the period
    of transits that follows repeats with the factor ``Pi = prod a e^{-tau}``, summed in closed form
    (``1 - Pi`` by ``expm1``).
    """
    with mp.extradps(40):           # b = r sqrt(1 - mu^2) must keep the half-chord r |mu| at mu ~ 1e-12 (F5, 2^-40)
        return +_psi_disk(edges, sigma, q, amplitude, r, mu, stretch)


def _psi_disk(edges, sigma, q, amplitude, r, mu, stretch):
    r0, R = edges[0], edges[-1]
    b = r * mp.sqrt((1 - mu) * (1 + mu))
    s_x = -r * mu                                       # the point's position along the BACKWARD path
    hR = mp.sqrt((R - b) * (R + b))
    if r0 > 0 and b < r0:
        h0 = mp.sqrt((r0 - b) * (r0 + b))
        s_end, wall = (-h0, "inner") if s_x < 0 else (hR, "outer")
    else:
        s_end, wall = hR, "outer"
    first, tau_first = _leg(_passages(edges, b, s_x, s_end), sigma, q, stretch, mp.mpf(0))
    a1 = amplitude[wall]
    if a1 == 0:
        return first
    period, tau, gain = mp.mpf(0), mp.mpf(0), mp.mpf(1)
    for passages, ends_at in _transits(edges, b, wall):
        value, depth = _leg(passages, sigma, q, stretch, tau)
        period += gain * value
        tau += depth
        gain *= amplitude[ends_at]
    # gain is now the product of the period's amplitudes; the period's depth is tau
    log_product = mp.log(gain) - tau if gain > 0 else mp.ninf
    one_minus = -mp.expm1(log_product) if gain > 0 else mp.mpf(1)
    if one_minus == 0:
        # a lossless period of optical depth 0 (a grazing chord at a mirror, the nodes tanh-sinh puts at mu = 0 at
        # the wall): the limit of the period's integral over 1 - e^{-tau} is the outer region's q / Sigma
        return first + a1 * mp.exp(-tau_first) * q[-1] / sigma[-1]
    return first + a1 * mp.exp(-tau_first) * period / one_minus


def _mu_breaks(edges, r: Mpf) -> list[Mpf]:
    """The direction cosines at which the backward path from ``r`` is tangent to an interface or a wall below ``r``."""
    breaks = {mp.mpf(-1), mp.mpf(0), mp.mpf(1)}
    for rk in edges[:-1]:
        if 0 < rk < r:
            c = mp.sqrt((1 - rk / r) * (1 + rk / r))
            breaks |= {c, -c}
    return sorted(breaks)


@lru_cache(maxsize=None)
def phi_sphere(edges: tuple, sigma: tuple, q: tuple, outer: float, inner: float, r: float) -> float:
    """The scalar flux at radius ``r`` of a concentric sphere under specular laws: ``(1/2) int psi dmu``."""
    with mp.workdps(DPS):
        e, s, qq = _mp(edges), _mp(sigma), _mp(q)
        amp = {"outer": mp.mpf(repr(float(outer))), "inner": mp.mpf(repr(float(inner)))}
        x = mp.mpf(repr(float(r)))
        if x == 0:
            return float(phi_sphere_centre(edges, sigma, q, outer))
        value = _quad(lambda m: psi_disk(e, s, qq, amp, x, m), _mu_breaks(e, x)) / 2
        return float(value)


def phi_sphere_centre(edges: tuple, sigma: tuple, q: tuple, outer: float) -> Mpf:
    """The centre of a solid sphere: every direction a diameter of rank 1, the elementary closed form."""
    with mp.workdps(DPS):
        e, s, qq = _mp(edges), _mp(sigma), _mp(q)
        a = mp.mpf(repr(float(outer)))
        radial = [(e[j + 1] - e[j], j) for j in range(len(s))]
        F, tau = _leg(radial, s, qq, mp.mpf(1), mp.mpf(0))          # centre to the wall
        inward, _ = _leg(radial[::-1], s, qq, mp.mpf(1), mp.mpf(0))  # the wall to the centre
        B = inward + mp.exp(-tau) * F                               # the diameter, attenuated to its end
        return F + a * mp.exp(-tau) * B / (-mp.expm1(mp.log(a) - 2 * tau) if a > 0 else 1) if a > 0 else F


# ── the cylinder ─────────────────────────────────────────────────────────


def ki2(x: Mpf) -> Mpf:
    """Bickley's ``Ki_2(x) = int_0^{pi/2} cos t e^{-x / cos t} dt = x (K_1(x) - Ki_1(x))`` (Bickley's recurrence at
    n = 1, ``Ki_1`` by Struve's closed form, :func:`._characteristic_mp.bickley_ki1`); ``Ki_2(0) = 1``."""
    if x == 0:
        return mp.mpf(1)
    with mp.extradps(10):
        return +(x * (mp.besselk(1, x) - bickley_ki1(x)))


@lru_cache(maxsize=None)
def phi_cylinder_axis(edges: tuple, sigma: tuple, q: tuple, outer: float) -> float:
    """The axis of a solid cylinder: ``(1/2) int_0^pi sin(theta) psi(theta) dtheta``, the in-plane diameter's closure."""
    with mp.workdps(DPS):
        e, s, qq = _mp(edges), _mp(sigma), _mp(q)
        amp = {"outer": mp.mpf(repr(float(outer))), "inner": mp.mpf(0)}
        value = _quad(lambda t: mp.sin(t) * psi_disk(e, s, qq, amp, mp.mpf(0), mp.mpf(0), 1 / mp.sin(t)),
                        [0, mp.pi / 2, mp.pi]) / 2
        return float(value)


@lru_cache(maxsize=None)
def phi_cylinder_axis_vacuum(edges: tuple, sigma: tuple, q: tuple) -> float:
    """The axis of a solid cylinder under vacuum: ``sum_j q_j / Sigma_j (Ki_2(tau_j) - Ki_2(tau_{j+1}))``."""
    with mp.workdps(DPS):
        e, s, qq = _mp(edges), _mp(sigma), _mp(q)
        tau, total = mp.mpf(0), mp.mpf(0)
        for j in range(len(s)):
            end = tau + s[j] * (e[j + 1] - e[j])
            total += qq[j] / s[j] * (ki2(tau) - ki2(end))
            tau = end
        return float(total)


@lru_cache(maxsize=None)
def phi_cylinder_vacuum(edges: tuple, sigma: tuple, q: tuple, r: float) -> float:
    r"""A solid cylinder under vacuum at radius ``r``: ``(1/pi) int_0^pi sum q/Sigma (Ki_2(tau_s) - Ki_2(tau_e)) dalpha``.

    The polar integral of each in-plane passage is ``2 Ki_2``; ``alpha`` is the in-plane angle of the backward path,
    split at the tangencies ``sin(alpha) = r_k / r``.
    """
    with mp.workdps(DPS):
        e, s, qq = _mp(edges), _mp(sigma), _mp(q)
        x = mp.mpf(repr(float(r)))
        hR = None

        def integrand(alpha):
            mu = mp.cos(alpha)                       # the in-plane cosine of the direction of flight, outward
            b = x * mp.sin(alpha)
            h = mp.sqrt((e[-1] - b) * (e[-1] + b))
            passages = _passages(e, b, -x * mu, h)
            tau, total = mp.mpf(0), mp.mpf(0)
            for length, j in passages:
                end = tau + s[j] * length
                total += qq[j] / s[j] * (ki2(tau) - ki2(end))
                tau = end
            return total

        breaks = {mp.mpf(0), mp.pi / 2, mp.pi}
        for rk in e[1:-1]:
            if rk < x:
                a = mp.asin(rk / x)
                breaks |= {a, mp.pi - a}
        return float(_quad(integrand, sorted(breaks)) / mp.pi)


@lru_cache(maxsize=None)
def phi_cylinder(edges: tuple, sigma: tuple, q: tuple, outer: float, r: float) -> float:
    r"""A solid cylinder under a specular law at radius ``r``: ``(1/pi) int_0^pi dalpha int_0^{pi/2} sin(theta) psi dtheta``.

    The two-dimensional route (``alpha`` the in-plane angle of the direction of flight from the outward radius,
    ``theta`` the polar angle): each direction's backward path is the in-plane disk billiard of
    :func:`psi_disk`, stretched by ``1/sin(theta)``; ``alpha`` split at ``pi/2`` and at the tangencies. Costly
    (seconds per point): the slow rows only.
    """
    with mp.workdps(20):
        e, s, qq = _mp(edges), _mp(sigma), _mp(q)
        amp = {"outer": mp.mpf(repr(float(outer))), "inner": mp.mpf(0)}
        x = mp.mpf(repr(float(r)))
        breaks = {mp.mpf(0), mp.pi / 2, mp.pi}
        for rk in e[1:-1]:
            if rk < x:
                a = mp.asin(rk / x)
                breaks |= {a, mp.pi - a}
        f = lambda al, th: mp.sin(th) * psi_disk(e, s, qq, amp, x, mp.cos(al), 1 / mp.sin(th))
        return float(_quad(f, sorted(breaks), [0, mp.pi / 4, mp.pi / 2]) / mp.pi)


@lru_cache(maxsize=None)
def specular_cylinder_total(sigma: float, R: float, a: float) -> float:
    r"""``1^T K 1`` per unit height of a homogeneous cylinder behind a partial mirror ``a`` (q = 1, K the flux over 4 pi).

    The integral over lines of ``int psi ds``: per line of impact ``b`` and polar angle ``theta``, chord
    ``l = 2 sqrt(R^2 - b^2) / sin(theta)``, ``B = (1 - e^{-Sigma l}) / Sigma``, ``psi_in = a B / (1 - a e^{-Sigma l})``,
    ``int psi ds = psi_in B + (l - B) / Sigma``; the line measure over 4 pi is ``2 sin^2(theta) db dtheta`` on
    ``[0, R] x [0, pi/2]`` (a closed body gives ``pi R^2 / Sigma``).
    """
    with mp.workdps(20):
        S, Rr, aa = (mp.mpf(repr(float(v))) for v in (sigma, R, a))

        def per_line(b, th):
            length = 2 * mp.sqrt((Rr - b) * (Rr + b)) / mp.sin(th)
            B = -mp.expm1(-S * length) / S
            return aa * B * B / (1 - aa * mp.exp(-S * length)) + (length - B) / S

        inner = lambda b: _quad(lambda th: 2 * mp.sin(th) ** 2 * per_line(b, th), [0, mp.pi / 4, mp.pi / 2])
        return float(_quad(inner, [0, Rr * mp.mpf("0.9"), Rr * mp.mpf("0.999"), Rr]))


@lru_cache(maxsize=None)
def specular_sphere_layered_total(edges: tuple, sigma: tuple, q: tuple, a: float) -> float:
    r"""``1^T K 1`` of a solid layered sphere behind a partial mirror ``a``: ``int_0^R 2 pi b db int psi ds`` per line.

    The rank-1 chord at impact ``b``, its segments by region (void allowed: Sigma = 0, q = 0 there): the vacuum total
    ``V = int int_{s' < s} q e^{-tau(s', s)}``, the source integral attenuated to the exit ``B``, the entry response
    ``A = int e^{-tau(entry, s)} ds``; ``int psi ds = V + A a B / (1 - a e^{-tau})``. ``q`` is the emission rate per
    region (``1`` on the emitting regions reads the block's ``1^T K 1``). Each impact panel ``[r_k, r_{k+1}]`` is
    integrated in ``y = sqrt(r_{k+1}^2 - b^2)`` with breakpoints ``span 2^{-j}`` toward ``y = 0``, where the closure's
    pole sits under a void or thin outer shell (qa's F1).
    """
    with mp.workdps(20):
        e, s, qq = _mp(edges), _mp(sigma), _mp(q)
        aa = mp.mpf(repr(float(a)))
        R = e[-1]

        def per_line(b):
            hR = mp.sqrt((R - b) * (R + b))
            V = A = B = mp.mpf(0)
            psi, depth = mp.mpf(0), mp.mpf(0)          # the vacuum flux carried, the depth from the entry
            parts = _passages(e, b, -hR, hR)
            for length, j in parts:
                sj, qj = s[j], qq[j]
                if sj == 0:
                    V += psi * length + qj * length * length / 2
                    A += mp.exp(-depth) * length
                    psi += qj * length
                    continue
                t = -mp.expm1(-sj * length)
                V += psi * t / sj + qj / sj * (length - t / sj)
                A += mp.exp(-depth) * t / sj
                psi = psi * (1 - t) + qj / sj * t
                depth += sj * length
            B = psi                                     # the vacuum flux at the exit is the source attenuated to it
            closure = aa * B / (-mp.expm1(mp.log(aa) - depth)) if aa > 0 else mp.mpf(0)
            return V + A * closure

        total = mp.mpf(0)
        for k in range(len(e) - 1):
            lo, hi = e[k], e[k + 1]
            span = mp.sqrt((hi - lo) * (hi + lo))
            ys = sorted({mp.mpf(0), span, *[span * mp.mpf(2) ** -j for j in range(1, 40)]})
            # b db = -y dy on the panel: b = sqrt(lo^2 + (span - y)(span + y))
            total += _quad(lambda y: 2 * mp.pi * y * per_line(min(hi, mp.sqrt(lo * lo + (span - y) * (span + y)))), ys)
        return float(total)


# ── the slab ─────────────────────────────────────────────────────────────


@lru_cache(maxsize=None)
def phi_slab(edges: tuple, sigma: tuple, q: tuple, left: float, right: float, x: float,
             cutoff: float = 1e-24) -> float:
    r"""The slab's ``E_2`` image series: ``(1/2) sum`` over the unfolded path's passages of ``w q/Sigma (E_2(t_s) - E_2(t_e))``."""
    with mp.workdps(DPS):
        e, s, qq = _mp(edges), _mp(sigma), _mp(q)
        aL, aR = mp.mpf(repr(float(left))), mp.mpf(repr(float(right)))
        point = mp.mpf(repr(float(x)))

        def passage(lo, hi, forward):
            out = []
            for j in range(len(s)):
                a, b = max(e[j], lo), min(e[j + 1], hi)
                if b > a:
                    out.append((j, b - a))
            return out if forward else out[::-1]

        total = mp.mpf(0)
        for leftward in (True, False):               # the backward path goes left (mu > 0) or right
            pieces = passage(e[0], point, False) if leftward else passage(point, e[-1], True)
            wall_left = leftward
            tau, weight = mp.mpf(0), mp.mpf(1)
            while True:
                for j, dx in pieces:
                    end = tau + s[j] * dx
                    total += weight * qq[j] / s[j] * (mp.expint(2, tau) - mp.expint(2, end))
                    tau = end
                weight *= aL if wall_left else aR
                if weight == 0 or weight * mp.expint(2, tau) < cutoff:
                    break
                pieces = passage(e[0], e[-1], wall_left)  # a full transit to the other wall
                wall_left = not wall_left
        return float(total / 2)


# ── white walls on a homogeneous body ────────────────────────────────────


@lru_cache(maxsize=None)
def phi_white_sphere(sigma: float, q: float, R: float, albedo: float, r: float) -> float:
    r"""A homogeneous sphere behind a white wall of albedo ``a``, at radius ``r``.

    The first flight ``(1/2) int psi_vac dmu`` plus the re-entering current: of the emission ``Q = q V``,
    ``E = Q P_esc`` escapes; the wall returns ``j = a (E + T j)`` (``T`` the uncollided transmission of an isotropic
    current, ``P_ss``), entering isotropically with ``psi = j / (pi A)`` per steradian; at the point it is attenuated
    over the distance to the wall: ``phi_ret = (j / (pi A)) 2 pi int e^{-Sigma d(mu)} dmu``. Escape and transmission
    are Hebert's closed forms.
    """
    with mp.workdps(DPS):
        S, Q, Rr, a, x = (mp.mpf(repr(float(v))) for v in (sigma, q, R, albedo, r))
        tau = S * Rr
        pesc = 3 / (8 * tau ** 3) * (2 * tau ** 2 - 1 + (1 + 2 * tau) * mp.exp(-2 * tau))
        T = (1 - (1 + 2 * tau) * mp.exp(-2 * tau)) / (2 * tau ** 2)
        V, A = 4 * mp.pi * Rr ** 3 / 3, 4 * mp.pi * Rr ** 2
        E = Q * V * pesc
        j = a * E / (1 - a * T)

        def distance(m):                              # to the wall, backward from the point, direction cosine m
            b2 = x * x * (1 - m) * (1 + m)
            return x * m + mp.sqrt(Rr * Rr - b2)

        breaks = [mp.mpf(-1), mp.mpf(0), mp.mpf(1)]
        first = _quad(lambda m: Q / S * -mp.expm1(-S * distance(m)), breaks) / 2
        returned = j / (mp.pi * A) * 2 * mp.pi * _quad(lambda m: mp.exp(-S * distance(m)), breaks)
        return float(first + returned)


@lru_cache(maxsize=None)
def phi_white_slab(sigma: float, q: float, L: float, left: float, right: float, x: float) -> float:
    r"""A homogeneous slab of width ``L`` whose faces are white with albedos ``left``, ``right``, at ``x``.

    The first flight ``(q / 2 Sigma)(2 - E_2(Sigma x) - E_2(Sigma (L - x)))`` plus each face's re-entering current
    ``2 j E_2(Sigma d)``: the escaping current per unit area ``E = q L P_esc`` leaves half by each face, and
    ``j_1 = a_1 (E/2 + T j_2)``, ``j_2 = a_2 (E/2 + T j_1)`` with ``T = 2 E_3(Sigma L)``.
    """
    with mp.workdps(DPS):
        S, Q, W, a1, a2, xx = (mp.mpf(repr(float(v))) for v in (sigma, q, L, left, right, x))
        tau = S * W
        pesc = (1 - 2 * mp.expint(3, tau)) / (2 * tau)
        T = 2 * mp.expint(3, tau)
        E = Q * W * pesc
        j1 = a1 * (E / 2 + a2 * T * E / 2) / (1 - a1 * a2 * T * T)
        j2 = a2 * (E / 2 + T * j1)
        first = Q / (2 * S) * (2 - mp.expint(2, S * xx) - mp.expint(2, S * (W - xx)))
        return float(first + 2 * j1 * mp.expint(2, S * xx) + 2 * j2 * mp.expint(2, S * (W - xx)))
