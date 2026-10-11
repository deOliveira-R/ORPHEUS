r"""The Peierls ray layer on the geometric kernel (#578, P2 group P): :meth:`CurvilinearGeometry.rays_from` and its consumers.

Verification spec ``scratch/characteristic_architecture/p2/gate_spec.md``
§P.1-§P.3 (§P.1b for the run to a distance). Until #578 the layer was four
scalar primitives (``rho_max``, ``rho_inner_intersections``,
``optical_depth_along_ray``, ``which_annulus``) called once per direction; it
is now one batch of the kernel's :class:`~orpheus.geometry.chord.Rays` per
observer, read through ``Chord.line_integral``, ``Rays.ahead``,
``Rays.to_wall``, ``Rays.wall``, ``Rays.lengths_to`` and ``Chord.cosine_at``.
The kernel verbs themselves are gated in ``tests/gates/geometry/test_rays.py``.

**The ray oracle** is written here in mpmath (40 digits) from the ray's own
algebra: the distances at which the ray from ``r`` with cosine ``mu`` meets
the circle ``r_k`` are the roots ``-r mu +- sqrt(r_k^2 - r^2 (1 - mu^2))``
(the slab's ``(r_k - r)/mu``), and a piece between two roots is attributed to
the region holding its midpoint (X4, derivation axis; the inputs built here).

**The consumers.** Measured before the carve (`[M]` 2026-10-10, ``ta/spy``):
0 of 16 936 calls the not-slow suite made into the retired
``optical_depth_along_ray`` were multi-region, and qa's review (2026-10-11,
finding F1) measured that a row comparing each consumer with itself on the
oracle's rays cannot see the consumer's CHOICE of reading (six source
transforms moved consumers by 2e-2 to 1.8 and reddened 0 of 776 rows). That
row was retired; the consumers are pinned instead by laws that tie two
consumers together (the void limit of the escape split, reciprocity) and by
references that write the consumer's integral independently (the formula
rows), each red under its named transform.
"""
from __future__ import annotations

import functools
from typing import Any

import mpmath as mp
import numpy as np
import pytest

from orpheus.derivations.common.quadrature_recipes import observer_angular_quadrature, surface_centred_angular_quadrature
from orpheus.derivations.continuous.peierls_nystrom import geometry as G
from orpheus.derivations.continuous.peierls_nystrom.geometry import CurvilinearGeometry
from orpheus.geometry import chord as kernel_chord

_EPS = np.finfo(float).eps
_HERE = "tests/gates/derivations/test_peierls_ray.py::"
_CHORD = "tests/gates/geometry/test_chord.py::"
_BETWEEN = "tests/gates/geometry/test_chord_lengths_between.py::"
_RAYS = "tests/gates/geometry/test_rays.py::"

_RADII = np.array([0.6, 1.3, 2.0])
_SIGMA = np.array([0.7, 3.1, 0.2])          # non-monotone, distinct at every interface
_R0 = 0.25


# ── the oracle ───────────────────────────────────────────────────────────────


def _roots(kind: str, r: Any, mu: Any, b: Any) -> list:
    """The distances ``s > 0`` at which the ray from ``r`` with cosine ``mu`` meets the level ``b``."""
    if kind == "slab-polar":
        return [(b - r) / mu] if mu != 0 else []
    disc = (r * mu) ** 2 - (r * r - b * b)
    if disc < 0:
        return []
    return [-r * mu - mp.sqrt(disc), -r * mu + mp.sqrt(disc)]


@functools.lru_cache(maxsize=None)
def _ray_mp(kind: str, r0: float, radii: tuple, r: float, mu: float):
    """The oracle's reading of one ray: ``(pieces, exit, near_wall, far_wall)``; pieces are ``(start, end, region)``
    with region ``None`` in the cavity and beyond the outer surface, all mpmath."""
    mp.mp.dps = 40
    bps = [mp.mpf(r0), *(mp.mpf(x) for x in radii)]
    rr, m = mp.mpf(r), mp.mpf(mu)
    roots = sorted({s for b in bps for s in _roots(kind, rr, m, b) if s > 0})
    pos = (lambda s: rr + s * m) if kind == "slab-polar" else (lambda s: mp.sqrt(rr * rr + 2 * rr * s * m + s * s))

    def region(s):
        p = pos(s)
        if p > bps[-1] or (kind == "slab-polar" and p < bps[0]):
            return "out"
        if p < bps[0]:
            return "cavity"
        return next(j for j in range(len(radii)) if p <= bps[j + 1])

    pts = [mp.mpf(0), *roots]
    pieces, exit_, near, far = [], None, None, None
    for a, b in zip(pts, [*pts[1:], pts[-1] + 4 * bps[-1]]):
        reg = region((a + b) / 2)
        if reg == "out":
            exit_ = a
            break
        if reg == "cavity":
            near, far = (a if a > 0 else None), b
        pieces.append((a, b, reg))
    return pieces, exit_, near, far


def _depth_mp(pieces, stop, sigma) -> Any:
    tau = mp.mpf(0)
    for a, b, reg in pieces:
        if isinstance(reg, int) and a < stop:
            tau += mp.mpf(float(sigma[reg])) * (min(b, stop) - a)
    return tau


# ── P.2: the observer's rays against the oracle ─────────────────────────────

_GEOMETRIES = [("sphere-1d", 0.0), ("sphere-1d", _R0), ("cylinder-1d", 0.0), ("cylinder-1d", _R0), ("slab-polar", 0.0)]


def _directions(kind: str, r0: float, r: float, rng) -> tuple[np.ndarray, int]:
    """Cosines: generic ones (axis, perpendicular, random) and, after index ``n_generic``, the exact tangents to every
    circle inside the observer and their neighbours at +-1e-9."""
    mus = [-1.0, -0.3, 0.3, 1.0, *rng.uniform(-1.0, 1.0, 6)]
    if kind != "slab-polar":
        mus.append(0.0)
    n_generic = len(mus)
    if kind != "slab-polar":
        for b in (r0, *_RADII[:-1]):
            if 0.0 < b < r:
                t = -float(np.sqrt(1.0 - (b / r) ** 2))
                mus += [t, t + 1e-9, t - 1e-9]
    return np.array(mus), n_generic


def _observers(r0: float) -> list[float]:
    return [r0 if r0 > 0 else 0.0, 0.4, 0.6, 1.0, 1.3, 1.7, 2.0]


def _depth(rays, lengths) -> np.ndarray:
    return rays.chord.line_integral(_SIGMA, lengths)


@pytest.mark.l0
@pytest.mark.verifies("geometry-chord-segment-lengths")
@pytest.mark.rests_on(_RAYS + "test_the_run_to_the_first_wall_and_the_wall_are_the_oracles",
                      _BETWEEN + "test_the_lengths_between_two_parameters_are_the_roots_pieces")
@pytest.mark.parametrize(("kind", "r0"), _GEOMETRIES, ids=[f"{k}-r0={r}" for k, r in _GEOMETRIES])
def test_the_optical_depths_of_the_rays_are_the_oracles(kind: str, r0: float) -> None:
    """P-tau: the optical depth to the exit (``ahead``), to the first wall where it is the cavity's (``to_wall``) and to
    a random distance (``lengths_to``), against the oracle, on three regions of non-monotone Sigma_t, observers at the
    cavity wall or centre, inside each region, on each breakpoint and on the surface.

    Bands, `[M]` 2026-10-10: generic directions within 20 eps (curvilinear) and 128 eps (slab, long grazing paths)
    against ``R Sigma_max``, so ``64 eps max(tau, R Sigma_max)``; a direction tangent to a circle (to rounding) is the
    problem's own conditioning, a chord of up to ``2 sqrt(eps) r`` carried at ``Sigma_max``: ``4 sqrt(eps) R
    Sigma_max`` (measured 6.8e-8 against 3.7e-7). First reds: the cosine's sign flipped in the direction mint; the cavity
    filled with the outermost value; ``to_wall`` taken as ``ahead``.
    """
    rng = np.random.default_rng(578)
    g = CurvilinearGeometry(kind=kind, inner_radius=r0)
    R, smax = float(_RADII[-1]), float(_SIGMA.max())
    for r in _observers(r0):
        mus, n_generic = _directions(kind, r0, r, rng)
        rays = g.rays_from(r, mus, _RADII)
        d = rng.uniform(0.01, 3.0, mus.shape)
        exit_tau, wall_tau, up_tau = _depth(rays, rays.ahead), _depth(rays, rays.to_wall), _depth(rays, rays.lengths_to(d))
        for i, mu in enumerate(mus):
            pieces, exit_, near, far = _ray_mp(kind, r0, tuple(_RADII), r, float(mu))
            tangent = i >= n_generic
            checks = [(exit_tau[i], exit_, "exit"), (up_tau[i], mp.mpf(float(d[i])), "up_to")]
            if far is not None and int(rays.wall[i]) == 0:
                checks.append((wall_tau[i], near if near is not None else mp.mpf(0), "cavity"))
            for got, stop, what in checks:
                want = float(_depth_mp(pieces, stop, _SIGMA))
                band = 4 * np.sqrt(_EPS) * R * smax if tangent else 64 * _EPS * max(abs(want), R * smax)
                assert abs(got - want) <= band, f"{kind} r0={r0} r={r} mu={mu!r} {what}: {got!r} vs {want!r}"


@pytest.mark.l0
@pytest.mark.rests_on(_HERE + "test_the_optical_depths_of_the_rays_are_the_oracles")
@pytest.mark.parametrize("kind", ["sphere-1d", "cylinder-1d"])
def test_the_first_wall_is_the_cavity_exactly_when_the_oracle_says(kind: str) -> None:
    """P-presence: ``wall == 0`` (the first wall ahead is the cavity's) exactly where the oracle's ray meets ``r_0``
    ahead, on generic directions from observers off the cavity wall. First red: ``wall`` read from the transit's entry.
    """
    rng = np.random.default_rng(579)
    g = CurvilinearGeometry(kind=kind, inner_radius=_R0)
    for r in (0.4, 0.6, 1.0, 1.3, 1.7, 2.0):
        mus, n_generic = _directions(kind, _R0, r, rng)
        wall = g.rays_from(r, mus[:n_generic], _RADII).wall
        want = [_ray_mp(kind, _R0, tuple(_RADII), r, float(m))[3] is not None for m in mus[:n_generic]]
        np.testing.assert_array_equal(wall == 0, want, err_msg=f"{kind} r={r}")


@pytest.mark.l0
@pytest.mark.rests_on(_HERE + "test_the_first_wall_is_the_cavity_exactly_when_the_oracle_says",
                      _RAYS + "test_a_point_on_a_wall_sits_on_its_crossing")
@pytest.mark.parametrize("kind", ["sphere-1d", "cylinder-1d"])
def test_an_observer_on_the_cavity_wall_sees_the_cavity_only_looking_inward(kind: str) -> None:
    """From ``r = r_0`` an outward ray's first wall is the outer surface, an inward one's the cavity, at distance 0.

    Found 2026-10-10 while writing these gates (the observer's parameter and its crossing's rounded differently: a
    1.4e-17 cavity slot ahead on the sphere at mu = 0.4027, 2.8e-17 on the cylinder at mu = 0.979, 0.961); fixed the
    same day. This row reads the wall; the parameter's exactness is the kernel's row it rests on
    (``test_rays.py::test_a_point_on_a_wall_sits_on_its_crossing``), which the pre-fix spelling reds (``ta/battery`` arm
    K_start_from_line); under the kernel's ``transit`` rule (closing crossing at or ahead of the point) a 1-ulp
    disagreement no longer turns the wall here, so this row is green under that arm (`[M]` 2026-10-11).
    """
    g = CurvilinearGeometry(kind=kind, inner_radius=_R0)
    outward = np.array([0.40267903915837233, 0.9793148738097697, 0.9606546472234416, 0.1, 0.5, 0.9])
    np.testing.assert_array_equal(g.rays_from(_R0, outward, _RADII).wall, np.full(6, 3))
    rays = g.rays_from(_R0, np.array([-0.5, -0.9, -0.1]), _RADII)
    np.testing.assert_array_equal(rays.wall, np.zeros(3, dtype=int))
    np.testing.assert_array_equal(rays.to_wall.sum(axis=-1), np.zeros(3))


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_optical_depths_of_the_rays_are_the_oracles")
@pytest.mark.parametrize(("kind", "r0"), _GEOMETRIES, ids=[f"{k}-r0={r}" for k, r in _GEOMETRIES])
def test_a_batch_of_rays_is_its_rays_one_by_one(kind: str, r0: float) -> None:
    """P-batch: the batch's readings equal the one-ray readings bit for bit (every kernel step is per line)."""
    rng = np.random.default_rng(580)
    g = CurvilinearGeometry(kind=kind, inner_radius=r0)
    mus, _ = _directions(kind, r0, 1.3, rng)
    batch = g.rays_from(1.3, mus, _RADII)
    d = rng.uniform(0.01, 3.0, mus.shape)
    for i, mu in enumerate(mus):
        one = g.rays_from(1.3, np.array([mu]), _RADII)
        for a, b in ((_depth(batch, batch.ahead)[i], _depth(one, one.ahead)[0]),
                     (_depth(batch, batch.to_wall)[i], _depth(one, one.to_wall)[0]),
                     (_depth(batch, batch.lengths_to(d))[i], _depth(one, one.lengths_to(d[i:i + 1]))[0]),
                     (batch.chord.cosine_at(2.0)[i], one.chord.cosine_at(2.0)[0]), (batch.wall[i], one.wall[0])):
            assert a == b, (kind, r0, mu, a, b)


@pytest.mark.foundation
@pytest.mark.parametrize("kind", ["sphere-1d", "cylinder-1d", "slab-polar"])
def test_the_region_of_a_point_is_inner_owns(kind: str) -> None:
    """P-region: :meth:`CurvilinearGeometry.region_of` is the kernel's locator: a point on ``r_k`` (k >= 1) belongs to
    the region below it (``which_annulus`` answered the one above until #578); interior points to their region."""
    g = CurvilinearGeometry(kind=kind)
    np.testing.assert_array_equal(g.region_of(np.array([0.6, 1.3, 2.0]), _RADII), [0, 1, 2])
    np.testing.assert_array_equal(g.region_of(np.array([0.3, 1.0, 1.7]), _RADII), [0, 1, 2])


# ── P.3a: laws tying two consumers together (qa F1) ─────────────────────────

_NODES = np.array([0.3, 0.8, 1.5, 1.95])
_HOLLOW_NODES = np.array([0.4, 0.9, 1.5, 1.95])


def _run(name: str, kind: str, r0: float, mode=None, *, radii=_RADII, sigma=_SIGMA, nodes=_NODES, n: int = 8):
    """The consumer, resolved on the module by NAME at call time (a battery's rebinding reaches it)."""
    kw = {"dps": 15, ("n_surf_quad" if "G_bc" in name else "n_angular"): n}
    args = (CurvilinearGeometry(kind=kind, inner_radius=r0), nodes, np.asarray(radii, float), np.asarray(sigma, float))
    return np.asarray(getattr(G, name)(*args, *(() if mode is None else (mode,)), **kw), dtype=float)


@pytest.mark.l1
@pytest.mark.rests_on(_HERE + "test_the_first_wall_is_the_cavity_exactly_when_the_oracle_says")
@pytest.mark.parametrize("kind", ["sphere-1d", "cylinder-1d"])
def test_in_a_void_hollow_body_every_neutron_escapes_through_one_surface(kind: str) -> None:
    """The escape split's void limit: at Sigma_t = 0 every ray from an observer in the shell leaves through exactly one
    surface, the cavity's if it meets the cavity first and the outer one otherwise, so ``P_esc_outer + P_esc_inner = 1``
    per observer (within 1e-13; `[M]` 2026-10-11 exactly 1 to the printed digits).
    First reds (qa's F1 arms, re-spelled on the kernel ``Rays``): ``compute_P_esc_outer`` counting the rays that meet the
    cavity first (``wall == n`` replaced by ``ahead > 0``); ``compute_P_esc_inner`` reading the run to the exit.
    """
    void = np.zeros(3)
    total = _run("compute_P_esc_outer", kind, _R0, sigma=void) + _run("compute_P_esc_inner", kind, _R0, sigma=void)
    np.testing.assert_allclose(total, np.ones(len(_NODES)), rtol=0, atol=1e-13)


#: (response G, escape P, geometry, r0, modes): pairs whose integrals are the same rays and the same kernel, so
#: reciprocity holds observer by observer: G = 4 P (`[M]` 2026-10-11, every pair below to the printed 10 digits).
#: The hollow per-face pairs (outer and inner, sphere and cylinder) joined after the ruling of 2026-10-11 (ERR-111,
#: ERR-112): each face counts the rays whose FIRST wall it is; until then the outer response counted the rays that
#: cross the cavity (G/P = 5.23 at r = 0.3 on the sphere) and the cylinder's responses were the surface-centred
#: ``Ki_1/d`` forms over the whole circle. Since that ruling ``compute_G_bc_outer`` and ``compute_G_bc_inner`` are
#: DEFINED as 4 times their escape (one integrand per face, Pattern 2), so on those four rows the law is the witness of
#: the single source (``retirement-audit`` D.18), not an independent pin; the independent pins of the per-face values are
#: ``test_peierls_hollow_faces.py`` (the flat-flux closure identity and the 3-D surface reference).
_RECIPROCAL = [
    ("compute_G_bc", "compute_P_esc", "sphere-1d", 0.0, (None,)),
    ("compute_G_bc", "compute_P_esc", "sphere-1d", _R0, (None,)),
    ("compute_G_bc_outer", "compute_P_esc_outer", "sphere-1d", 0.0, (None,)),
    ("compute_G_bc_outer", "compute_P_esc_outer", "sphere-1d", _R0, (None,)),
    ("compute_G_bc_outer", "compute_P_esc_outer", "cylinder-1d", _R0, (None,)),
    ("compute_G_bc_inner", "compute_P_esc_inner", "sphere-1d", _R0, (None,)),
    ("compute_G_bc_inner", "compute_P_esc_inner", "cylinder-1d", _R0, (None,)),
    ("compute_G_bc_outer_mode", "compute_P_esc_outer_mode", "sphere-1d", _R0, (1, 2)),
    ("compute_G_bc_inner_mode", "compute_P_esc_inner_mode", "sphere-1d", _R0, (0, 1)),
    ("compute_G_bc_outer_mode_marshak", "compute_P_esc_outer_mode_marshak", "sphere-1d", _R0, (1,)),
    ("compute_G_bc_inner_mode_marshak", "compute_P_esc_inner_mode_marshak", "sphere-1d", _R0, (1,)),
    ("compute_G_bc_cylinder_3d_mode", "compute_P_esc_cylinder_3d_mode", "cylinder-1d", 0.0, (0, 1)),
    ("compute_G_bc_cylinder_3d_mode", "compute_P_esc_cylinder_3d_mode", "cylinder-1d", _R0, (0, 1)),
]


@pytest.mark.l1
@pytest.mark.catches("ERR-111", "ERR-112")
@pytest.mark.parametrize(("g_name", "p_name", "kind", "r0", "modes"), _RECIPROCAL,
                         ids=[f"{g}-{k}-r0={r}" for g, _, k, r, _ in _RECIPROCAL])
def test_each_response_is_four_times_its_escape(g_name: str, p_name: str, kind: str, r0: float, modes) -> None:
    """Reciprocity per surface and mode: the uncollided flux at r from a unit isotropic current entering through a surface
    is 4 times the probability that a neutron born isotropically at r leaves through it (the same rays, the same kernel,
    the same surface weight), within 1e-13 relative, on three regions of non-monotone Sigma_t.
    First reds (qa's F1 arms): ``compute_G_bc_outer_mode`` not skipping the rays that meet the cavity first;
    ``compute_P_esc_inner_mode`` reading the exit cosine at the outer surface; ``compute_G_bc_inner`` (sphere) reading
    the run to the exit instead of the run to the cavity wall; ``compute_G_bc_outer`` counting the rays that cross the
    cavity (HEAD's pairing, ERR-111).
    """
    nodes = _HOLLOW_NODES if r0 > 0 else _NODES
    for mode in modes:
        g = _run(g_name, kind, r0, mode, nodes=nodes)
        p = _run(p_name, kind, r0, mode, nodes=nodes)
        assert np.max(np.abs(p)) > 0.0, f"{p_name}: vacuous (identically 0)"
        np.testing.assert_allclose(g, 4.0 * p, rtol=1e-13, atol=0.0, err_msg=f"{g_name} n={mode}")


# ── P.3b: consumers against their integrals written independently ──────────


def _ki_mp(n: int, x) -> Any:
    """Bickley ``Ki_n(x) = int_0^{pi/2} cos^{n-1}(t) exp(-x / cos t) dt`` by quadrature (independent of ``ki_n_mp``)."""
    x = mp.mpf(x)
    return mp.quad(lambda t: mp.cos(t) ** (n - 1) * mp.exp(-x / mp.cos(t)), [0, mp.pi / 2])


def _p_tilde(n: int, mu) -> Any:
    """The shifted Legendre polynomial ``P_n(2 mu - 1)``."""
    return mp.legendre(n, 2 * mp.mpf(mu) - 1)


def _segment_depth_mp(r0: float, radii, sigma, r: float, phi: float, r_surface: float) -> tuple[Any, Any]:
    """The length ``d`` of the in-plane segment from ``(r, 0)`` to ``r_surface (cos phi, sin phi)`` and its optical depth,
    the pieces inside ``r_0`` (the cavity) void, mpmath."""
    mp.mp.dps = 40
    P = (mp.mpf(r), mp.mpf(0))
    Q = (mp.mpf(r_surface) * mp.cos(mp.mpf(phi)), mp.mpf(r_surface) * mp.sin(mp.mpf(phi)))
    d = mp.sqrt((Q[0] - P[0]) ** 2 + (Q[1] - P[1]) ** 2)
    if d == 0:
        return d, mp.mpf(0)
    u = ((Q[0] - P[0]) / d, (Q[1] - P[1]) / d)
    b = 2 * (P[0] * u[0] + P[1] * u[1])
    c0 = P[0] ** 2 + P[1] ** 2
    cuts = [mp.mpf(0), d]
    for rk in (r0, *radii):
        disc = b * b - 4 * (c0 - mp.mpf(rk) ** 2)
        if disc > 0:
            cuts += [s for s in ((-b - mp.sqrt(disc)) / 2, (-b + mp.sqrt(disc)) / 2) if 0 < s < d]
    cuts.sort()
    tau = mp.mpf(0)
    for a, e in zip(cuts[:-1], cuts[1:]):
        rho = mp.sqrt(c0 + b * (a + e) / 2 + ((a + e) / 2) ** 2)
        if r0 < rho <= radii[-1]:                       # a sliver past R (a root at the surface point) is outside
            tau += mp.mpf(float(sigma[next(j for j in range(len(radii)) if rho <= radii[j])])) * (e - a)
    return d, tau


def _rule_breakpoints(r0: float, radii) -> np.ndarray:
    """The breakpoints the drivers hand the angular rules: every positive one, a hollow body's ``r_0`` included (since
    the hollow-body rulings of 2026-10-11; the cavity's tangency is a kink of every hollow integrand)."""
    return np.array([b for b in (r0, *radii) if b > 0.0], dtype=float)


def _observer_reference(kind: str, r0: float, radii, sigma, nodes, prefactor, integrand) -> np.ndarray:
    """``prefactor * sum_j w_j integrand(omega_j, ray_j)`` on the consumer's own observer rule (a trusted recipe the carve
    does not touch), the ray read from the oracle."""
    out = []
    for r in nodes:
        q = observer_angular_quadrature(r_obs=float(r), omega_low=0.0, omega_high=float(np.pi), radii=_rule_breakpoints(r0, radii),
                                        n_per_panel=8, dps=15)
        total = mp.mpf(0)
        for w, om in zip(q.wts, q.pts):
            total += mp.mpf(float(w)) * integrand(float(r), float(om), _ray_mp(kind, r0, tuple(radii), float(r),
                                                                                  float(np.cos(om))))
        out.append(float(prefactor * total))
    return np.array(out)


def _surface_reference(r0: float, radii, sigma, nodes, r_surface: float, prefactor, weight) -> np.ndarray:
    """``prefactor * sum_j w_j weight(r, phi, d) Ki_1(tau) / d`` on the consumer's own surface-centred rule."""
    out = []
    for r in nodes:
        q = surface_centred_angular_quadrature(r_obs=float(r), r_surface=float(r_surface), radii=_rule_breakpoints(r0, radii),
                                               n_per_panel=8, phi_low=0.0, phi_high=float(np.pi), dps=15)
        total = mp.mpf(0)
        for w, phi in zip(q.wts, q.pts):
            d, tau = _segment_depth_mp(r0, radii, sigma, float(r), float(phi), r_surface)
            if d > 0:
                total += mp.mpf(float(w)) * weight(float(r), float(phi), d) * _ki_mp(1, tau) / d
        out.append(float(prefactor * total))
    return np.array(out)


def _exit_cosine(r: float, om: float, radius: float) -> Any:
    """The cosine to the normal of the circle ``radius`` where the ray from ``r`` at angle ``om`` meets it."""
    s2 = (mp.mpf(r) * mp.sin(mp.mpf(om))) ** 2
    return mp.sqrt(max(mp.mpf(radius) ** 2 - s2, 0)) / mp.mpf(radius)


def _escape_sphere(n_mode, *, jacobian: bool, R: float):
    """Integrand of the sphere's outward escape at mode n: ``sin(om) [(d/R)^2] P~_n(mu_exit) e^{-tau_exit}``."""
    def f(r, om, ray):
        pieces, exit_, near, far = ray
        if exit_ is None or exit_ <= 0:
            return mp.mpf(0)
        weight = (exit_ / R) ** 2 if jacobian else 1
        return mp.sin(om) * weight * _p_tilde(n_mode, _exit_cosine(r, om, R)) * mp.exp(-_depth_mp(pieces, exit_, _SIGMA_REF[0]))
    return f


_SIGMA_REF: list = [_SIGMA]      # the Sigma_t the integrand closures read (set per row)


def _inner_escape_sphere(n_mode, r0: float):
    """Integrand of the sphere's escape through the cavity wall at mode n: rays meeting the cavity first,
    ``sin(om) P~_n(mu_wall) e^{-tau_to_wall}`` with ``mu_wall`` the cosine at ``r_0``."""
    def f(r, om, ray):
        pieces, exit_, near, far = ray
        if far is None or near is None:
            return mp.mpf(0)
        return mp.sin(om) * _p_tilde(n_mode, _exit_cosine(r, om, r0)) * mp.exp(-_depth_mp(pieces, near, _SIGMA_REF[0]))
    return f


def _knyazev(n_mode, R: float):
    """Integrand of the cylinder-3D mode n: ``sum_k c_k mu^k Ki_{k+2}(tau_exit)``, ``P~_n(mu) = sum_k c_k mu^k`` expanded
    here (n <= 1: ``P~_0 = 1``, ``P~_1 = 2 mu - 1``)."""
    coefs = {0: (1,), 1: (-1, 2)}[n_mode]
    def f(r, om, ray):
        pieces, exit_, near, far = ray
        if exit_ is None or exit_ <= 0:
            return mp.mpf(0)
        tau, mu = _depth_mp(pieces, exit_, _SIGMA_REF[0]), _exit_cosine(r, om, R)
        return sum(c * mu ** k * _ki_mp(k + 2, tau) for k, c in enumerate(coefs) if c)
    return f


def _to_cavity_ki2():
    """Integrand of the cylinder's inner response after ERR-112: the rays whose first wall is the cavity's,
    ``Ki_2(tau_to_wall)`` (the 3-D fold of the in-plane ray). Until 2026-10-11 these two oracles wrote the legacy
    surface-centred ``Ki_1/d`` form over the whole circle, the far side included, and the rows encoded the defect."""
    def f(r, om, ray):
        pieces, exit_, near, far = ray
        if far is None or near is None:
            return mp.mpf(0)
        return _ki_mp(2, _depth_mp(pieces, near, _SIGMA_REF[0]))
    return f


#: (consumer, geometry, r0, mode, reference) on the three-region fixture: the formula rows of qa's F1 transforms that no
#: law ties to a sibling.
def _references_3region():
    R = float(_RADII[-1])
    return {
        "compute_P_esc_inner_mode-sphere": ("compute_P_esc_inner_mode", "sphere-1d", _R0, 1, _HOLLOW_NODES,
            lambda nodes: _observer_reference("sphere-1d", _R0, _RADII, _SIGMA, nodes, mp.mpf(1) / 2, _inner_escape_sphere(1, _R0))),
        "compute_P_esc_mode-hsphere": ("compute_P_esc_mode", "sphere-1d", _R0, 1, _HOLLOW_NODES,
            lambda nodes: _observer_reference("sphere-1d", _R0, _RADII, _SIGMA, nodes, mp.mpf(1) / 2, _escape_sphere(1, jacobian=True, R=R))),
        "compute_G_bc_mode-cylinder": ("compute_G_bc_mode", "cylinder-1d", 0.0, 1, _NODES,
            lambda nodes: _surface_reference(0.0, tuple(_RADII), _SIGMA, nodes, R, 2 * mp.mpf(R) / mp.pi,
                                             lambda r, phi, d: _p_tilde(1, (R - r * mp.cos(mp.mpf(phi))) / d))),
    }


_REF_3 = _references_3region()


@pytest.mark.l1
@pytest.mark.rests_on(_HERE + "test_the_optical_depths_of_the_rays_are_the_oracles")
@pytest.mark.parametrize("key", list(_REF_3))
def test_each_consumer_is_its_integral_on_three_regions(key: str) -> None:
    """The consumer against its integral written here (the docstring's formula: kernel, surface weight, cosine, Jacobian
    and the rays' presence), every ray read from the oracle, on the consumer's own angular rule, within 1e-11 relative.
    First reds (qa's F1 arms): ``compute_P_esc_inner_mode`` reading the exit cosine at the outer surface;
    ``compute_P_esc_mode`` skipping the rays that meet the cavity first; ``compute_G_bc_mode`` (cylinder) spelling
    ``mu_s`` as ``(R cos phi - r) / d``.
    """
    name, kind, r0, mode, nodes, reference = _REF_3[key]
    _SIGMA_REF[0] = _SIGMA
    got = _run(name, kind, r0, mode, nodes=nodes)
    want = reference(nodes)
    assert np.max(np.abs(want)) > 0.0
    np.testing.assert_allclose(got, want, rtol=1e-11, atol=1e-14 * np.max(np.abs(want)), err_msg=key)


# ── P.3c: hollow single-region bodies (qa F2: HEAD's one-region shortcuts filled the cavity) ────

_H_RADII, _H_SIGMA, _H_R0 = (2.0,), np.array([0.8]), 0.3


def _references_hollow():
    R = 2.0
    obs = lambda kind, pref, f: (lambda nodes: _observer_reference(kind, _H_R0, _H_RADII, _H_SIGMA, nodes, pref, f))
    surf = lambda rs, pref, wgt: (lambda nodes: _surface_reference(_H_R0, _H_RADII, _H_SIGMA, nodes, rs, pref, wgt))
    return {
        "compute_G_bc-sphere": ("compute_G_bc", "sphere-1d", None, obs("sphere-1d", 2, _escape_sphere(0, jacobian=False, R=R))),
        "compute_G_bc_mode-sphere-n1": ("compute_G_bc_mode", "sphere-1d", 1, obs("sphere-1d", 2, _escape_sphere(1, jacobian=False, R=R))),
        # since ERR-114 (2026-10-11) the cylinder's G_bc is 4 P_esc, the observer form (4/pi) int Ki_2(tau_exit); until
        # then this oracle wrote the legacy surface-centred Ki_1/d form (0.586 to 1.340 x a 3-D reference)
        "compute_G_bc-cylinder": ("compute_G_bc", "cylinder-1d", None, obs("cylinder-1d", 4 / mp.pi, _knyazev(0, R))),
        "compute_G_bc_mode-cylinder-n1": ("compute_G_bc_mode", "cylinder-1d", 1,
                                          surf(R, 2 * mp.mpf(R) / mp.pi, lambda r, phi, d: _p_tilde(1, (R - r * mp.cos(mp.mpf(phi))) / d))),
        "compute_G_bc_inner-cylinder": ("compute_G_bc_inner", "cylinder-1d", None, obs("cylinder-1d", 4 / mp.pi, _to_cavity_ki2())),
        "compute_G_bc_inner_mode-cylinder-n0": ("compute_G_bc_inner_mode", "cylinder-1d", 0,
                                                obs("cylinder-1d", 4 / mp.pi, _to_cavity_ki2())),
        "compute_G_bc_cylinder_3d": ("compute_G_bc_cylinder_3d", "cylinder-1d", None, obs("cylinder-1d", 4 / mp.pi, _knyazev(0, R))),
        "compute_P_esc_cylinder_3d_mode-n1": ("compute_P_esc_cylinder_3d_mode", "cylinder-1d", 1, obs("cylinder-1d", 1 / mp.pi, _knyazev(1, R))),
        "compute_G_bc_cylinder_3d_mode-n1": ("compute_G_bc_cylinder_3d_mode", "cylinder-1d", 1, obs("cylinder-1d", 4 / mp.pi, _knyazev(1, R))),
    }


_REF_H = _references_hollow()


@pytest.mark.l1
@pytest.mark.rests_on(_HERE + "test_the_optical_depths_of_the_rays_are_the_oracles")
@pytest.mark.parametrize("key", [pytest.param(k, marks=pytest.mark.catches("ERR-112" if "inner" in k else "ERR-109"))
                                 for k in _REF_H])
def test_a_hollow_single_region_body_has_a_void_cavity_in_each_consumer(key: str) -> None:
    """qa's F2: on a hollow single-region body (r_0 = 0.3, R = 2, Sigma_t = 0.8) each consumer HEAD got wrong equals its
    integral written here with the cavity void, within 1e-11 relative. HEAD's ``len(radii) == 1`` shortcuts took the
    optical depth as ``Sigma_t d``, the cavity counted as material: `[M]` qa 2026-10-11, 1.4 % to 22 % (``qa/q7``,
    ``q8``). First red: the cavity filled with the region's Sigma_t (``ta/battery`` arm K_li_cavity, the shortcut's
    arithmetic): 7 of 9 rows. The two inner-surface cylinder rows are blind to it by construction since ERR-112 (the
    inner face counts only the rays whose first wall is the cavity, which never cross it), so they carry ERR-112, not
    ERR-109; their first red is the far side counted (``ta/battery`` arm H_inner_farside).

    """
    name, kind, mode, reference = _REF_H[key]
    _SIGMA_REF[0] = _H_SIGMA
    got = _run(name, kind, _H_R0, mode, radii=np.array(_H_RADII), sigma=_H_SIGMA, nodes=_HOLLOW_NODES)
    want = reference(_HOLLOW_NODES)
    assert np.max(np.abs(want)) > 0.0
    np.testing.assert_allclose(got, want, rtol=1e-11, atol=1e-14 * np.max(np.abs(want)), err_msg=key)


# ── the cylinder-3D pairing and the route ───────────────────────────────────

_CONSUMERS = [
    ("compute_P_esc", "sphere", None), ("compute_P_esc", "cylinder", None),
    ("compute_G_bc", "sphere", None), ("compute_G_bc", "cylinder", None),
    ("compute_G_bc_cylinder_3d", "cylinder", None), ("compute_G_bc_cylinder_3d", "hcylinder", None),
    ("compute_P_esc_outer", "hsphere", None), ("compute_P_esc_outer", "hcylinder", None),
    ("compute_P_esc_inner", "hsphere", None), ("compute_P_esc_inner", "hcylinder", None),
    ("compute_P_esc_outer_mode", "hsphere", 1), ("compute_P_esc_inner_mode", "hsphere", 1),
    ("compute_G_bc_outer", "hsphere", None), ("compute_G_bc_outer", "hcylinder", None),
    ("compute_G_bc_inner", "hsphere", None), ("compute_G_bc_inner", "hcylinder", None),
    ("compute_G_bc_outer_mode", "hsphere", 1), ("compute_G_bc_inner_mode", "hsphere", 1),
    ("compute_P_esc_outer_mode_marshak", "hsphere", 1), ("compute_P_esc_inner_mode_marshak", "hsphere", 1),
    ("compute_G_bc_outer_mode_marshak", "hsphere", 1), ("compute_G_bc_inner_mode_marshak", "hsphere", 1),
    ("compute_P_esc_mode", "sphere", 1), ("compute_P_esc_mode", "hcylinder", 2),
    ("compute_G_bc_mode", "sphere", 1), ("compute_G_bc_mode", "cylinder", 1),
    ("compute_P_esc_cylinder_3d_mode", "cylinder", 1), ("compute_G_bc_cylinder_3d_mode", "hcylinder", 1),
]
_GEOM = {"sphere": ("sphere-1d", 0.0), "cylinder": ("cylinder-1d", 0.0),
         "hsphere": ("sphere-1d", _R0), "hcylinder": ("cylinder-1d", _R0)}

def _call(name: str, geom_key: str, mode):
    """The consumer, resolved on the module by NAME at call time, on the three-region fixture."""
    kind, r0 = _GEOM[geom_key]
    fn = getattr(G, name)
    kw = {"dps": 15, ("n_surf_quad" if "G_bc" in name else "n_angular"): 8}
    args = (CurvilinearGeometry(kind=kind, inner_radius=r0), _NODES, _RADII, _SIGMA)
    return np.asarray(fn(*args, *(() if mode is None else (mode,)), **kw), dtype=float)



@pytest.mark.l1
@pytest.mark.rests_on(_HERE + "test_the_optical_depths_of_the_rays_are_the_oracles")
@pytest.mark.catches("ERR-108")
@pytest.mark.parametrize("name", ["compute_G_bc_cylinder_3d", "compute_G_bc_cylinder_3d_mode"])
def test_the_cylinder_3d_responses_are_their_characteristic_integrals(name: str) -> None:
    r"""The cylinder-3D response at the observer is ``(4/pi) int_0^pi Ki_2(tau_2D(omega)) d omega`` (mode 0; the mode
    row ``n = 0`` of the Knyazev expansion is the same integrand), with ``tau_2D`` the optical depth of ONE in-plane ray
    from the observer to its exit, evaluated here entirely in mpmath (the oracle's ray, Bickley by quadrature) on the
    consumer's own angular rule.

    Until #578 both consumers paired the exit distance along ``+cos omega`` with the optical depth along ``-cos
    omega``, which agree only on a homogeneous cylinder: `[M]` 11.0 % on this three-region cylinder (main's
    ``p_baseline/cyl3d_oracle.py``; the new code 1.3e-16). First red: that pairing restored (the ray read from the
    mirrored cosine for the depth only).
    """
    from orpheus.derivations.common.quadrature_recipes import observer_angular_quadrature

    g = CurvilinearGeometry(kind="cylinder-1d")
    kw = {"dps": 15, "n_surf_quad": 8}
    got = np.asarray(G.compute_G_bc_cylinder_3d(g, _NODES, _RADII, _SIGMA, **kw) if name == "compute_G_bc_cylinder_3d"
                     else G.compute_G_bc_cylinder_3d_mode(g, _NODES, _RADII, _SIGMA, 0, **kw), dtype=float)
    want = []
    for r in _NODES:
        q = observer_angular_quadrature(r_obs=float(r), omega_low=0.0, omega_high=float(np.pi), radii=_RADII,
                                        n_per_panel=8, dps=15)
        total = mp.mpf(0)
        for w, om in zip(q.wts, q.pts):
            pieces, exit_, *_ = _ray_mp("cylinder-1d", 0.0, tuple(_RADII), float(r), float(np.cos(om)))
            total += mp.mpf(float(w)) * _ki_mp(2, _depth_mp(pieces, exit_, _SIGMA))
        want.append(float(4 / mp.pi * total))
    want = np.array(want)
    np.testing.assert_allclose(got, want, rtol=1e-11, atol=0.0, err_msg=f"{name}: the paired ray")


@pytest.mark.foundation
@pytest.mark.parametrize(("name", "geom_key", "mode"), _CONSUMERS,
                         ids=[f"{n}-{g}" + ("" if m is None else f"-n{m}") for n, g, m in _CONSUMERS])
def test_each_consumer_reads_one_batch_of_rays_per_observer(name: str, geom_key: str, mode,
                                                             monkeypatch: pytest.MonkeyPatch) -> None:
    """PC-route (the batching ruling of 2026-10-10, ruling 3; no timing gate): each consumer forms one chord per
    observer, ``ConcentricPartition.chord`` called exactly ``len(r_nodes)`` times per call, and reaches the kernel at
    all; the retired scalar primitives are gone (an attribute error, the retirement as a refusal).
    First red: a chord per direction before the batch (`ta/battery` arm P_one_by_one, wrapping
    `CurvilinearGeometry.rays_from`: 24 of 28 rows red; the 4 surface-centred rows mint their rays on the kernel directly
    through `_to_surface_point`, which that arm does not reach).
    """
    for retired in ("rho_max", "rho_inner_intersections", "optical_depth_along_ray", "which_annulus", "source_position"):
        assert not hasattr(CurvilinearGeometry, retired), retired
    for retired in ("ObserverRays", "Run"):
        assert not hasattr(G, retired), retired
    calls = {"chord": 0}
    honest = kernel_chord.ConcentricPartition.chord

    def counting(self, line, level=None):
        calls["chord"] += 1
        return honest(self, line, level)

    monkeypatch.setattr(kernel_chord.ConcentricPartition, "chord", counting)
    _call(name, geom_key, mode)
    assert calls["chord"] == len(_NODES), f"{name} {geom_key}: {calls['chord']} chords for {len(_NODES)} observers"
