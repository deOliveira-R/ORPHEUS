r"""The per-face primitives of a hollow body, after the #578 P2 hollow-body rulings (ERR-111, ERR-112, ERR-113).

The numerics-investigator's memo (``scratch/characteristic_architecture/p2/hollow_investigation/memo.md``) ruled the
per-face model: the volume kernel stops each ray at its FIRST wall, so the escape ``P_k`` and the response ``G_k`` of
each face count the same rays, those whose first wall is face k, and reciprocity makes ``G_k = 4 P_k`` observer by
observer. Rows here, each with its first red (``scratch/characteristic_architecture/p2/ta/battery``):

* the closure identity, the decisive rung: a pure absorber whose source is ``Sigma_t phi_0`` region by region, behind
  a white outer face and around a void cavity, has the exact solution ``phi = phi_0`` (an isotropic angular flux
  satisfies the transport equation, both faces and the cavity); the rank-2 per-face closure built from the production
  ``P``, ``G`` and ``W`` reproduces it to its discretisation's accuracy, for two groups of a pure absorber (thin and
  thick) on one region (the transmissions refuse more, so the heterogeneous rung waits for a multi-region ``W``);
* a hollow body equals the solid body whose core region is void, on the outer-response consumers (the cavity's
  tangency in the angular rules);
* the inner cylinder response against an independent 3-D surface-point reference (no Bickley function, no chord);
* ``W_oo`` of the hollow cylinder at ``Sigma_t > 0`` against the solid cylinder's Bickley ``P_ss``;
* an observer on an interior wall converging in the angular rule.
"""
from __future__ import annotations

from typing import Any

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq

from orpheus.derivations.continuous.peierls_nystrom import geometry as G
from orpheus.derivations.continuous.peierls_nystrom.geometry import CurvilinearGeometry
from tests.gates.derivations import _characteristic_mp as closed

_HERE = "tests/gates/derivations/test_peierls_hollow_faces.py::"
_RAY = "tests/gates/derivations/test_peierls_ray.py::"

_RADII = np.array([0.6, 1.3, 2.0])
_SIGMAS = {"g1": np.array([0.7, 3.1, 0.2]), "g2": np.array([1.9, 0.4, 1.1])}     # three regions, two cross-section sets
#: The closure's fixture: one region (the transmissions ``W`` refuse more, a declared scope edge of
#: ``compute_hollow_{sph,cyl}_transmission``), two groups of a pure absorber, thin and thick.
_CLOSURE_RADII = np.array([2.0])
_CLOSURE_SIGMAS = {"g1": np.array([0.8]), "g2": np.array([2.5])}
_R0 = 0.25


def _call(name: str, kind: str, r0: float, nodes, radii, sigma, mode=None, n: int = 8):
    """The consumer, resolved on the module by NAME at call time."""
    kw = {"dps": 15, ("n_surf_quad" if "G_bc" in name else "n_angular"): n}
    args = (CurvilinearGeometry(kind=kind, inner_radius=r0), np.asarray(nodes, float), np.asarray(radii, float),
            np.asarray(sigma, float))
    return np.asarray(getattr(G, name)(*args, *(() if mode is None else (mode,)), **kw), dtype=float)


# ── the closure identity ────────────────────────────────────────────────────


def _shell_nodes(r0: float, radii, per_region: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Gauss-Legendre radial nodes, ``per_region`` in each region, their weights and their region."""
    x, w = np.polynomial.legendre.leggauss(per_region)
    edges = [r0, *radii]
    nodes, weights, region = [], [], []
    for k, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        nodes.append(0.5 * (b - a) * x + 0.5 * (b + a))
        weights.append(0.5 * (b - a) * w)
        region.append(np.full(per_region, k))
    return np.concatenate(nodes), np.concatenate(weights), np.concatenate(region)


_CLOSURE = {"sphere-1d": (G.compute_hollow_sph_transmission, lambda r: 4.0 * np.pi * r * r),
            "cylinder-1d": (G.compute_hollow_cyl_transmission, lambda r: 2.0 * np.pi * r)}


def _closure_residual(kind: str, sigma: np.ndarray, *, per_region: int = 24, n: int = 16) -> float:
    """``max |phi / phi_0 - 1|`` of the rank-2 per-face closure on the pure absorber with ``Q = Sigma_t phi_0``.

    ``phi(r)/phi_0 = 1 - P_o(r) - P_i(r) + sum_k G_k(r) J_k``: the volume part is ``(1/4pi) int (1 - e^{-tau_wall})
    dOmega`` because ``int Sigma e^{-tau} ds = 1 - e^{-tau_wall}`` along each ray to its first wall, whatever the
    regions; ``J = (I - W)^{-1} S / A`` with ``S_k = int Q P_k dV`` the current leaving the volume through face k,
    ``A_k`` the face's area (per unit height on the cylinder).
    """
    transmission, area = _CLOSURE[kind]
    radii = _CLOSURE_RADII
    R = float(radii[-1])
    nodes, weights, region = _shell_nodes(_R0, radii, per_region)
    P = np.stack([_call("compute_P_esc_outer", kind, _R0, nodes, radii, sigma, n=n),
                  _call("compute_P_esc_inner", kind, _R0, nodes, radii, sigma, n=n)], axis=1)
    Gr = np.stack([_call("compute_G_bc_outer", kind, _R0, nodes, radii, sigma, n=n),
                   _call("compute_G_bc_inner", kind, _R0, nodes, radii, sigma, n=n)], axis=1)
    W = np.asarray(transmission(_R0, R, radii, sigma), dtype=float)
    q = sigma[region]                                                     # Q = Sigma_t phi_0, phi_0 = 1
    S = (area(nodes) * weights * q) @ P
    J = np.linalg.solve(np.eye(2) - W, S) / np.array([area(R), area(_R0)])
    phi = 1.0 - P.sum(axis=1) + Gr @ J
    return float(np.max(np.abs(phi - 1.0)))


#: The resolution per chart (radial Gauss nodes on [r_0, R], angular n) and the bands, 10 x `[M]` 2026-10-11 rounded
#: up. The residual is the discretisation's, falling under both refinements (sphere g1: 3.2e-5 at (24, 32), 7.2e-6 at
#: (48, 32), 4.4e-6 at (48, 64); g2: 9.5e-5, 4.4e-5, 1.4e-5); the cylinder runs coarser for cost (14 s a group at
#: (24, 16)): 6.4e-5 (g1) and 1.7e-4 (g2).
_CLOSURE_RESOLUTION = {"sphere-1d": (48, 32), "cylinder-1d": (24, 16)}
_CLOSURE_BAND = {("sphere-1d", "g1"): 1e-4, ("sphere-1d", "g2"): 5e-4,
                 ("cylinder-1d", "g1"): 7e-4, ("cylinder-1d", "g2"): 2e-3}


@pytest.mark.l1
@pytest.mark.catches("ERR-111", "ERR-112")
@pytest.mark.rests_on(_RAY + "test_each_response_is_four_times_its_escape",
                      _HERE + "test_the_hollow_cylinders_outer_transmission_is_the_solid_ones_past_the_cavity")
@pytest.mark.parametrize("group", list(_CLOSURE_SIGMAS))
@pytest.mark.parametrize("kind", ["sphere-1d", "cylinder-1d"])
def test_the_per_face_closure_reproduces_the_flat_flux_of_a_pure_absorber(kind: str, group: str) -> None:
    """The rank-2 per-face closure, built from production ``P_o, P_i, G_o, G_i`` and ``W``, gives ``phi = phi_0`` on a
    hollow pure absorber with ``Q = Sigma_t phi_0`` (white outer face, void cavity), two groups; within
    ``_CLOSURE_BAND``. The fixture is one region (r_0 = 0.25, R = 2): the
    transmissions ``W`` refuse more (``NotImplementedError``, a declared scope edge), so the heterogeneous rung waits.
    First reds: the outer response counting the rays that meet the cavity first (HEAD's pairing; the memo: 6.4e-2 on
    the homogeneous sphere); the cylinder's legacy surface-centred ``Ki_1/d`` responses with the far side (1.08); the
    closure loop cut, ``W := 0`` (7.6e-2).
    """
    per_region, n = _CLOSURE_RESOLUTION[kind]
    residual = _closure_residual(kind, _CLOSURE_SIGMAS[group], per_region=per_region, n=n)
    band = _CLOSURE_BAND[(kind, group)]
    assert band is not None and residual <= band, f"{kind} {group}: max |phi/phi_0 - 1| = {residual:.3e}"


# ── a hollow body is the solid body with a void core ───────────────────────

_VOID_CORE = [("compute_G_bc", "sphere-1d", None), ("compute_G_bc", "cylinder-1d", None),
              ("compute_G_bc_mode", "sphere-1d", 1), ("compute_G_bc_mode", "cylinder-1d", 1),
              ("compute_G_bc_cylinder_3d", "cylinder-1d", None), ("compute_P_esc_cylinder_3d_mode", "cylinder-1d", 1),
              ("compute_G_bc_cylinder_3d_mode", "cylinder-1d", 1), ("compute_P_esc", "sphere-1d", None)]


@pytest.mark.l1
@pytest.mark.catches("ERR-109")
@pytest.mark.parametrize(("name", "kind", "mode"), _VOID_CORE, ids=[f"{n}-{k}" + ("" if m is None else f"-n{m}")
                                                                      for n, k, m in _VOID_CORE])
def test_a_hollow_body_is_the_solid_body_with_a_void_core(name: str, kind: str, mode) -> None:
    """Every consumer that reads rays through the cavity to the outer surface answers the same on the hollow body
    (r_0 = 0.3, one region to R = 2, Sigma_t = 0.8) and on the solid body whose core region [0, 0.3] is void, within
    1e-13 relative (`[M]` 2026-10-11: 0 exactly; 2.6e-3 to 0.30 before the angular rules received r_0).
    First reds: r_0 dropped from the drivers' breakpoints (the cavity's tangency unresolved); the cavity filled
    with the region's Sigma_t (HEAD's one-region shortcut, ERR-109).
    """
    nodes = np.array([0.4, 0.9, 1.5, 1.95])
    hollow = _call(name, kind, 0.3, nodes, [2.0], [0.8], mode)
    core = _call(name, kind, 0.0, nodes, [0.3, 2.0], [0.0, 0.8], mode)
    np.testing.assert_allclose(hollow, core, rtol=1e-13, atol=0.0, err_msg=name)


# ── the inner cylinder response against a 3-D surface-point reference ──────


def _seg_depth(p0: np.ndarray, p1: np.ndarray, r0: float, radii, sigma) -> float:
    """The optical depth of the segment ``p0 -> p1`` (in-plane), the part inside ``r_0`` void."""
    v = p1 - p0
    L = float(np.hypot(*v))
    if L == 0.0:
        return 0.0
    u = v / L
    a, b, c = 1.0, 2.0 * float(p0 @ u), float(p0 @ p0)
    cuts = [0.0, L]
    for rk in (r0, *radii):
        disc = b * b - 4 * a * (c - rk * rk)
        if disc > 0:
            cuts += [t for t in ((-b - np.sqrt(disc)) / 2, (-b + np.sqrt(disc)) / 2) if 0.0 < t < L]
    cuts.sort()
    tau = 0.0
    for s, e in zip(cuts[:-1], cuts[1:]):
        m = p0 + 0.5 * (s + e) * u
        rho = float(np.hypot(*m))
        if r0 < rho <= radii[-1]:
            tau += float(sigma[int(np.searchsorted(radii, rho))]) * (e - s)
    return tau


def _z_fold(tau2: float) -> float:
    """``int_{-pi/2}^{pi/2} cos t e^{-tau2 / cos t} dt``, numerically (no Bickley function)."""
    return quad(lambda t: np.cos(t) * np.exp(-tau2 / np.cos(t)) if np.cos(t) > 0 else 0.0,
                -np.pi / 2, np.pi / 2, epsabs=1e-14, epsrel=1e-13, limit=200)[0]


def _root(g, a: float, b: float) -> float:
    """The bracketed root of ``g`` on ``[a, b]`` (``brentq`` without ``full_output`` returns the root alone)."""
    root: Any = brentq(g, a, b, xtol=1e-15)
    return float(root)


def _cylinder_inner_response_3d(r: float, r0: float, radii, sigma) -> float:
    """The uncollided flux at ``(r, 0, 0)`` from a unit uniform isotropic current leaving the cavity's surface into the
    annulus, by surface-point integration over ``(phi, z)`` (the memo's ``ref3d.cyl_G_face``, ``near``): only the
    points whose outward normal faces the observer emit toward it; ``z`` folded numerically."""
    obs = np.array([r, 0.0])

    def f(ph: float) -> float:
        n = np.array([np.cos(ph), np.sin(ph)])
        s = r0 * n
        v = obs - s
        d = float(np.hypot(*v))
        cos_e = float(n @ v) / d
        if cos_e <= 0.0:
            return 0.0
        return (r0 / np.pi) * cos_e / d * _z_fold(_seg_depth(s, obs, r0, radii, sigma))

    edges = [0.0, 2.0 * np.pi]
    grid = np.linspace(1e-9, 2 * np.pi - 1e-9, 2001)
    for g in [lambda ph: float(np.array([np.cos(ph), np.sin(ph)]) @ (obs - r0 * np.array([np.cos(ph), np.sin(ph)])))]:
        gv = np.array([g(a) for a in grid])
        edges += [_root(g, grid[i], grid[i + 1])
                  for i in np.nonzero(np.sign(gv[:-1]) != np.sign(gv[1:]))[0]]
    edges = sorted(edges)
    return sum(quad(f, a, b, limit=200, epsabs=1e-12, epsrel=1e-11)[0] for a, b in zip(edges[:-1], edges[1:]) if b > a)


@pytest.mark.l1
@pytest.mark.catches("ERR-112")
@pytest.mark.parametrize("radii, sigma", [((2.0,), (0.8,)), ((0.6, 1.3, 2.0), (0.7, 3.1, 0.2))], ids=["one", "three"])
def test_the_inner_cylinder_response_is_the_3d_surface_integral(radii, sigma) -> None:
    """``compute_G_bc_inner`` on a hollow cylinder (r_0 = 0.3) equals the 3-D surface-point integral of the uncollided
    flux from the cavity's surface, near side only, within ``_INNER_3D_RTOL`` at n_surf_quad = 32. `[M]` 2026-10-11:
    7.9e-6, 1.0e-6, 1.3e-7 at n = 16, 32, 64, algebraic convergence (a kink the rule does not grade toward; reported
    to main), so the band is 10 x the n = 32 reading.
    First red: the legacy surface-centred ``Ki_1/d`` form over the whole circle (the far side 46 % to 58 % of it; near
    side alone still 0.89 to 1.15 x the reference, the memo).
    """
    radii, sigma = np.asarray(radii, float), np.asarray(sigma, float)
    nodes = np.array([0.4, 1.5])
    got = _call("compute_G_bc_inner", "cylinder-1d", 0.3, nodes, radii, sigma, n=32)
    want = np.array([_cylinder_inner_response_3d(float(r), 0.3, radii, sigma) for r in nodes])
    np.testing.assert_allclose(got, want, rtol=_INNER_3D_RTOL, atol=0.0)


_INNER_3D_RTOL = 1e-5


def _cylinder_outer_response_3d(r: float, radii, sigma) -> float:
    """The uncollided flux at ``(r, 0, 0)`` inside a SOLID cylinder from a unit uniform isotropic inward current on its
    outer surface, by surface-point integration over ``(phi, z)``, ``z`` folded numerically (the memo's
    ``ref3d.cyl_G_face``, outer face; on a solid body every surface point sees the observer)."""
    R = float(radii[-1])
    obs = np.array([r, 0.0])

    def f(ph: float) -> float:
        n = np.array([np.cos(ph), np.sin(ph)])
        s = R * n
        v = obs - s
        d = float(np.hypot(*v))
        cos_e = -float(n @ v) / d
        return (R / np.pi) * cos_e / d * _z_fold(_seg_depth(s, obs, 0.0, radii, sigma))

    return quad(f, 0.0, 2.0 * np.pi, limit=200, epsabs=1e-12, epsrel=1e-11)[0]


@pytest.mark.catches("ERR-114")
@pytest.mark.l1
@pytest.mark.parametrize("radii, sigma", [((2.0,), (0.8,)), ((0.6, 1.3, 2.0), (0.7, 3.1, 0.2))], ids=["one", "three"])
def test_the_solid_cylinder_response_is_the_3d_surface_integral(radii, sigma) -> None:
    """``compute_G_bc`` on a solid cylinder equals the 3-D surface-point integral of the uncollided flux from its outer
    surface, within ``_SOLID_3D_RTOL`` at n_surf_quad = 32 (ERR-114, #112's rank-1 half; the ruling of 2026-10-11 made
    it 4 times ``compute_P_esc``). First red: the legacy surface-centred ``Ki_1/d`` form (0.586 to 1.340 times this
    reference, the memo).
    """
    radii, sigma = np.asarray(radii, float), np.asarray(sigma, float)
    nodes = np.array([0.4, 1.5])
    got = _call("compute_G_bc", "cylinder-1d", 0.0, nodes, radii, sigma, n=32)
    want = np.array([_cylinder_outer_response_3d(float(r), radii, sigma) for r in nodes])
    np.testing.assert_allclose(got, want, rtol=_SOLID_3D_RTOL, atol=0.0)


#: 10 x `[M]` 2026-10-11 at n = 32: 2e-16 on one region; on three, 3.9e-16 at r = 0.4 and 5.2e-6 at r = 1.5, which
#: falls as n^-3 (6.7e-7 at 64, 8.5e-8 at 128): the observer form's integrand has a square-root kink at each shell's
#: tangent direction, a panel end the Gauss rule does not grade toward (a measured convergence limit of P2).
_SOLID_3D_RTOL = 6e-5


# ── W_oo at Sigma_t > 0 ─────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.catches("ERR-113")
@pytest.mark.parametrize("sigma", [0.5, 0.8, 1.0, 2.0])
def test_the_hollow_cylinders_outer_transmission_is_the_solid_ones_past_the_cavity(sigma: float) -> None:
    """``W_oo`` of a hollow cylinder with a vanishing cavity (r_0 = 1e-6, R = 1) is the solid cylinder's ``P_ss``
    (Bickley's closed form, ``_characteristic_mp``) less the rays that hit the cavity: ``0 <= P_ss - W_oo <= r_0 / R``
    (those rays span ``alpha < arcsin(r_0/R)`` with ``(4/pi) cos(alpha) Ki_3 <= 1``). `[M]` 2026-10-11: 2e-7 to 7.5e-7
    relative. First red: the ``Ki_3`` argument ``2 Sigma R cos(alpha)`` unclipped, -3e-26 at the node next to pi/2
    (``W_oo = -inf`` at every Sigma_t > 0; its only earlier gate ran at Sigma_t = 0, inside the defect's stabiliser).
    """
    r0, R = 1e-6, 1.0
    W = np.asarray(G.compute_hollow_cyl_transmission(r0, R, np.array([R]), np.array([sigma])), dtype=float)
    assert np.all(np.isfinite(W)), W
    gap = float(closed.transmission_probability("cylinder", sigma)) - float(W[0, 0])
    assert 0.0 <= gap <= r0 / R, f"P_ss - W_oo = {gap:.3e}"


# ── an observer on an interior wall ─────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.parametrize("r", [0.6, 1.3])
def test_an_observer_on_an_interior_wall_converges_in_the_angular_rule(r: float) -> None:
    """``compute_P_esc_outer`` at an observer exactly on a breakpoint of the solid three-region sphere converges in
    ``n_angular``: ``|P(16) - P(64)| <= 1e-7`` (`[M]` 2026-10-11: 6.6e-9 at r = 1.3, 1.3e-12 at r = 0.6; before the
    observer's own shell joined the rule, about 1e-2 at r = 1.3 and n = 16).
    The r = 1.3 sequence still converges algebraically (4.3e-8, 5.8e-9, 7.6e-10 per doubling), the r = 0.6 one
    spectrally; reported to main. First red: the observer rule taking only the shells strictly inside the observer.
    """
    P = {n: _call("compute_P_esc_outer", "sphere-1d", 0.0, [r], _RADII, _SIGMAS["g1"], n=n)[0] for n in (16, 64)}
    assert abs(P[16] - P[64]) <= 1e-7, P
