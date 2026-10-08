"""Gates for the characteristic reference's Galerkin assembly over lines (:class:`~orpheus.derivations.continuous.characteristic.assembly.LineRule`).

P1 step (b), third rung, of the characteristic-reference campaign
(``.claude/plans/characteristic_reference_architecture.md``, "P1 step (b), third
rung: API sketch", items 5-7, the user's rulings of 2026-10-06 (the cylinder's
polar angle) and the orchestrator's (the emission support, the piece budget)).
Verification spec ``scratch/characteristic_architecture/p1_step_b3/spec.md``:
rows MD3, EB4, LR1, AS1-AS5, WC1-WC7, OP1-OP4, LA1-LA4.

The object: one group's transport block K (rows on every panel, columns on the
emission support), the scalar flux over 4 pi of an isotropic emission, built as
the line part sum_L w_L (V_L + sum_forward A_k (x) in_k) plus the diffuse walls'
update R alpha (I - T alpha)^-1 U^T (``WallCoupling``).

The ladder, bottom up (``rests_on`` on each row):

1. the kernel: the one measure and its density, the line domain
   (``tests/gates/geometry/test_measure_density.py``, ``test_line_domain.py``);
   the basis, with its even panel (``test_characteristic_basis.py``); the
   period and its inflow, with an arriving flux (``test_characteristic_closure.py``);
   the line transport (``test_characteristic_transport.py``);
2. the routes: the mass and the wall area read the one density [MD3]; the
   cylinder's rule is Gauss in theta [LR1];
3. the closed-body theorems: conservation [AS1, AS1b, AS5], symmetry with its
   stabiliser declared [AS2], the piece budget is a partition [AS3], no diffuse
   wall adds nothing [AS4];
4. the walls against closed forms: the escape probability [WC1], the wall
   transmission [WC2], reciprocity [WC3, WC7], the white and specular laws [WC4],
   the refusals [WC5, WC6];
5. operator theorems on the material blocks [OP1-OP4];
6. self-convergence ladders below the working point [EB4, LA1-LA4].

Every closed form is written in ``_characteristic_mp.py`` (mpmath, no import of
the code under test): Hebert's sphere escape probability, (1 - 2 E_3)/(2 tau) for
the slab, and the cylinder's through Bickley's Ki_3 (Struve closed form of
int K_0 and the recurrence), the balance of one re-emission chain for the white
law, the per-line integral for the specular sphere.

Levels (the labels minted 2026-10-07, ``docs/theory/references/characteristic.rst``):
the rows whose docstring names a planned level carry it, with
``verifies("characteristic-boundary-resolvent")`` on the wall-coupling rows
(WC1-WC4, WC7, AS6) and ``verifies("characteristic-galerkin-assembly")`` on the
rest; EB4 and LA1-LA4 are ``l2`` (CONV). The others are ``foundation``.
``catches`` markers (ERR-101 fixed resolutions over lines, ERR-102 I - T alpha
by subtraction, ERR-103 the normal direction ungraded) sit only on rows that
redden when their defect is re-dropped under -O (the battery,
``scratch/characteristic_architecture/p1_step_b3/gates/battery/``). Rows on the
cylinder cost 30 s or more and are ``slow`` (`[M]` 2026-10-07, the orchestrator, one process per run at the piece
budget 1024, after the line rule evaluated only live intervals: a one-region white cylinder 27.8 s per group at
tau = 30 and 8 points, 102 s at 16; at tau = 100, 38 s and 145 s; the three-region cylinder 175 s at 8 points,
792 s at 16, #586).

The line rule's gradings are derived from the group's optical scale (the user's
ruling of 2026-10-06, after qa found the fixed rules blind): the grazing
coordinate halved to the thinnest absorbing panel's optical width, the impact
parameter graded in the chord half-length toward the next radius's branch point,
toward b = 0, and exponentially at the rim. The QA rows and the escape rows at
tau = 0.01, 1e-4, 100, 1000 are their witnesses; each one's first red is the
rule it replaced (battery arms A16-A21).
"""
from __future__ import annotations

import math
from functools import lru_cache

import mpmath as mp
import numpy as np
import pytest

from orpheus.derivations.continuous.characteristic import LineRule, PanelBasis, TrappedSource, Wall, Walls
from orpheus.derivations.continuous.characteristic.lines import Lines
from orpheus.derivations.common.quadrature import composite_gauss_legendre, gauss_legendre
from orpheus.derivations.continuous.characteristic.grading import graded_ends
from orpheus.geometry.chart import Chart
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.coord import CoordSystem, MeasureCoordinate

from . import _characteristic_mp as R

pytestmark = pytest.mark.filterwarnings("error::RuntimeWarning")

_EPS = float(np.finfo(float).eps)
_HERE = "tests/gates/derivations/test_characteristic_assembly.py::"
_BASIS = "tests/gates/derivations/test_characteristic_basis.py::"
_CLOSURE = "tests/gates/derivations/test_characteristic_closure.py::"
_TRANS = "tests/gates/derivations/test_characteristic_transport.py::"
_LINES = "tests/gates/geometry/test_line_domain.py::"
_DENSITY = "tests/gates/geometry/test_measure_density.py::"

_COORD = {"sphere": CoordSystem.SPHERICAL, "cylinder": CoordSystem.CYLINDRICAL, "slab": CoordSystem.CARTESIAN}
_MR3 = (0.0, 0.5, 1.5, 2.0)
_MR3H = (0.4, 0.5, 1.5, 2.0)
_SLB3 = (0.0, 0.4, 1.5, 2.3)
#: Sigma_t per region and group (the spec's §4 data): distinct at every interface and per group.
_SIGMA = {"g0": (0.6, 1.3, 0.45), "g1": (1.7, 0.35, 2.4), "hom": (0.7, 0.7, 0.7)}
_RES = (3, 2, 0.4)                 # (degree, layers, ratio): the premises' working point
_INNER = 12                        # outer and inner arc-length points

_MIRROR, _WHITE = (1.0, 0.0), (0.0, 1.0)


def _walls(chart: str, breakpoints, laws) -> Walls:
    """Walls built directly (the law reader is gated in ``test_characteristic_walls.py``): one (specular, diffuse) per
    boundary point, inner first, or ``"P"`` for a periodic face."""
    n = len(breakpoints) - 1
    ends = (n,) if chart != "slab" and breakpoints[0] == 0.0 else (0, n)
    walls = tuple(Wall(k, 1.0, 0.0, n - k) if law == "P" else Wall(k, law[0], law[1], k) for k, law in zip(ends, laws, strict=True))
    return Walls(walls, n, Chart(_COORD[chart]))


@lru_cache(maxsize=None)
def _basis(chart: str, breakpoints, resolution=_RES) -> PanelBasis:
    return PanelBasis.of(ConcentricPartition(Chart(_COORD[chart]), tuple(breakpoints)), *resolution)


@lru_cache(maxsize=None)
def _rule(chart: str, breakpoints, laws, sigma, points: int, budget: int = 1024, chunk: int = 512) -> LineRule:
    """The group's rule: its gradings (grazing, rim, small radius) are derived from the group's Sigma_t."""
    return LineRule.of(_basis(chart, breakpoints), _walls(chart, breakpoints, laws), np.asarray(sigma), points,
                       chunk=chunk, budget=budget)


@lru_cache(maxsize=None)
def _transport(chart: str, breakpoints, laws, sigma, points: int, inner: int = _INNER, budget: int = 1024):
    return _rule(chart, breakpoints, laws, sigma, points, budget).transport(inner, inner)


def _conservation(chart, breakpoints, laws, sigma, points, inner=_INNER) -> float:
    """max |K_g sigma_g - W 1| / max |W 1| over every panel's row."""
    basis = _basis(chart, breakpoints)
    g = _transport(chart, breakpoints, laws, sigma, points, inner)
    lhs = g.block @ np.asarray(sigma)[basis.region][g.support]
    rhs = basis.mass @ np.ones(basis.size)
    return float(np.max(np.abs(lhs - rhs)) / np.max(np.abs(rhs)))


# ── MD3: the one density's route ─────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.rests_on(_DENSITY + "test_the_density_integrates_to_the_one_measure",
                      _BASIS + "test_the_volume_density_is_the_charts_by_hand")
def test_the_mass_and_the_wall_area_read_the_one_density(monkeypatch: pytest.MonkeyPatch) -> None:
    """[MD3] ROUTE: with ``MeasureCoordinate.derivative`` doubled, the mass doubles and the wall's response halves,
    bit for bit (the arriving flux is 1/D, D = A/4); the transmission and the loss, fractions of the current each
    injection enters with, and the escape do not move; and the basis keeps no density of its own
    (``PanelBasis.volume_density`` retired, ``retirement-audit`` A.1).

    A decoy by a power of two is exact through every product and sum, so the
    legs are ``array_equal``. First red: a second spelling of the density left
    on the basis or in the coupling (it does not move under the decoy).
    """
    assert not hasattr(PanelBasis, "volume_density")
    laws = (_WHITE,)

    def build():
        basis = PanelBasis.of(ConcentricPartition(Chart(CoordSystem.SPHERICAL), _MR3), *_RES)
        g = LineRule.of(basis, _walls("sphere", _MR3, laws), np.asarray(_SIGMA["g0"]), 4).transport(4, 4)
        return basis.mass, g.coupling

    mass, honest = build()
    original = MeasureCoordinate.derivative
    monkeypatch.setattr(MeasureCoordinate, "derivative", lambda self, r: 2.0 * original(self, r))
    decoy_mass, decoy = build()
    np.testing.assert_array_equal(decoy_mass, 2.0 * mass)
    np.testing.assert_array_equal(decoy.response, 0.5 * honest.response)
    np.testing.assert_array_equal(decoy.transmission, honest.transmission)
    np.testing.assert_array_equal(decoy.loss, honest.loss)
    np.testing.assert_array_equal(decoy.escape, honest.escape)


# ── LR1: the cylinder's direction rule ───────────────────────────────────


@pytest.mark.foundation
@pytest.mark.catches("ERR-101")
@pytest.mark.rests_on(_LINES + "test_the_density_is_the_beam_density_times_the_folded_directions")
def test_the_cylinders_direction_rule_is_gauss_in_the_polar_angle() -> None:
    """[LR1, structure leg] The cylinder's lines sit at the Gauss points of [0, pi/2] in theta (not in mu_z), and each
    weight is the impact rule's weight times the theta weight times the domain's density over 4 pi.

    The user's ruling of 2026-10-06 (`[M]` the orchestrator's ``m12_theta.py``:
    C9 against Bickley 2.0e-6, 2.3e-9, 1.7e-12 at 8, 16, 32 points in theta;
    5.4e-4 at 8 in mu_z). The value leg is WC1's cylinder rows. First red:
    Gauss in mu_z (the nodes arccos of the mu_z Gauss points).
    """
    # a void cylinder: no absorbing panel, so the polar angle is not graded (plain Gauss on [0, pi/2])
    rule = _rule("cylinder", (0.0, 1.3), (_WHITE,), (0.0,), 8)
    theta = np.unique(rule.lines.coordinates[:, 1])
    plain = gauss_legendre(0.0, math.pi / 2, 8)
    np.testing.assert_allclose(theta, plain.pts, rtol=2 * _EPS, atol=0.0)
    domain = Chart(CoordSystem.CYLINDRICAL).line_domain()
    w_theta = dict(zip(plain.pts, plain.wts))
    # the rule orders its lines by projected speed (2026-10-06): put them back on the (b, theta) grid first
    grid = np.lexsort((rule.lines.coordinates[:, 1], rule.lines.coordinates[:, 0]))
    coordinates, weights = rule.lines.coordinates[grid], rule.lines.weights[grid]
    per_impact = weights / (domain.density(coordinates) / (4.0 * math.pi))
    by_theta = per_impact / np.array([w_theta[min(w_theta, key=lambda t: abs(t - x))] for x in coordinates[:, 1]])
    impact = by_theta.reshape(-1, 8)
    np.testing.assert_allclose(impact, impact[:, :1] * np.ones((1, 8)), rtol=8 * _EPS, atol=0.0)
    # an absorbing cylinder: the polar angle is graded toward 0 (the grazing direction), never in mu_z
    graded = np.unique(_rule("cylinder", (0.0, 1.3), (_WHITE,), (1.0,), 8).lines.coordinates[:, 1])
    assert graded.min() < plain.pts.min() / 4, "the polar angle is not graded toward the grazing direction"


# ── AS1: closed-body conservation (the re-posed C11) ─────────────────────

#: (id, chart, breakpoints, laws, points). The slab conserves under any mu rule integrating 1 and |mu| (the spec's
#: P-3: 33 of 33 runs <= 2.5e-15), so its rows run at 4 points (declared blind); the grading is derived from Sigma_t.
_SMALL_CAVITY = (1e-3, 0.5, 1.5, 2.0)
_CLOSED = [
    ("hsphere_small_cavity_white_mirror", "sphere", _SMALL_CAVITY, (_WHITE, _MIRROR), 16),
    ("sphere_small_first_region_mirror", "sphere", (0.0, 1e-3, 0.5, 1.5), (_MIRROR,), 16),
    ("sphere_mirror", "sphere", _MR3, (_MIRROR,), 16),
    ("sphere_white", "sphere", _MR3, (_WHITE,), 16),
    ("hsphere_mirror_mirror", "sphere", _MR3H, (_MIRROR, _MIRROR), 16),
    ("hsphere_white_mirror", "sphere", _MR3H, (_WHITE, _MIRROR), 16),
    ("hsphere_mirror_white", "sphere", _MR3H, (_MIRROR, _WHITE), 16),
    ("hsphere_white_white", "sphere", _MR3H, (_WHITE, _WHITE), 16),
    ("slab_mirror_mirror", "slab", _SLB3, (_MIRROR, _MIRROR), 4),
    ("slab_white_white", "slab", _SLB3, (_WHITE, _WHITE), 4),
    ("slab_mirror_white", "slab", _SLB3, (_MIRROR, _WHITE), 4),
    ("slab_periodic", "slab", _SLB3, ("P", "P"), 4),
]
#: The cylinder at both groups and 8 points (the user's ruling of 2026-10-07): one MR3 cylinder block takes 175 s at
#: 8 points and 792 s at 16 (`[M]` the orchestrator, one process, budget 1024), so the three-region rows stay at 8.
#: `[M]` 2026-10-07 at 8 points: conservation 5.4e-13 (g0), 1.7e-13 (g1, each of the three rows); AS2's symmetry 2.2e-16 (g1).
_CLOSED_SLOW = [
    ("cylinder_mirror", "cylinder", _MR3, (_MIRROR,), 8),
    ("cylinder_white", "cylinder", _MR3, (_WHITE,), 8),
    ("hcylinder_white_mirror", "cylinder", _MR3H, (_WHITE, _MIRROR), 8),
]
#: 10 x the spec's prototype figure per chart; the built code: sphere 4.1e-13, hollow sphere 1.7e-13, slab 8.9e-15,
#: cylinder 3.2e-11 (the orchestrator, 2026-10-06).
_CONSERVATION_TOL = {"sphere": 1e-11, "slab": 3e-14, "cylinder": 3e-10}


def _closed_params(rows, marks=(), groups=("g0", "g1", "hom")):
    return [pytest.param(*r[1:], g, id=f"{r[0]}-{g}", marks=marks) for r in rows for g in groups]


@pytest.mark.l1
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.parametrize(("chart", "breakpoints", "laws", "points", "group"),
                         _closed_params(_CLOSED) + _closed_params(_CLOSED_SLOW, marks=pytest.mark.slow, groups=("g0", "g1")))
@pytest.mark.rests_on(_TRANS + "test_the_volterra_triangle_is_the_double_integral_along_the_line",
                      _CLOSURE + "test_an_arriving_flux_enters_the_cycle_at_its_traversal",
                      _LINES + "test_cauchys_formula_on_the_line_domain")
def test_a_closed_body_conserves_its_emission(chart, breakpoints, laws, points, group) -> None:
    """[AS1, the re-posed C11; planned l1, ``characteristic-galerkin-assembly``] K_g sigma_g = W 1 on a closed body.

    psi = 1 solves the closed-body problem with emission Sigma_t(x), so the
    block applied to the nodal Sigma_t (piecewise constant, in the basis
    exactly) is the volume of each basis function. Mirrors, white walls at
    alpha = 1, a periodic slab, hollow bodies with each inner and outer law;
    both groups and a homogeneous column. Independence: the left side is the
    transport (lines, closure, coupling), the right side the volume measure
    only. First reds (`[M]` the spec's arms): the arriving flux without 1/D (D
    dropped): 6.6 (sphere white), 7.5e-3 (hollow, inner white), green on the
    mirrors; the density's exponent wrong; the impact rule at 8 points
    (3.4e-6); 4 outer or inner points (6.3e-7, 6.9e-8); Gauss in mu_z on the
    white cylinder (1.4e-4). Declared blind (measured): the slab's mu rule, the
    cylinder's direction rule under a mirror (each line conserves on its own),
    a short mass rule (W 1 is exact under it), the odd basis at the working point.
    """
    tol = _CONSERVATION_TOL[chart]
    err = _conservation(chart, breakpoints, laws, _SIGMA[group], points)
    assert err <= tol, f"{err:.2e} > {tol:.0e}"


@pytest.mark.l1
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.rests_on(_HERE + "test_a_closed_body_conserves_its_emission")
def test_the_closed_sphere_conserves_to_rounding_at_high_resolution() -> None:
    """[AS1b; planned l1] At 24 impact points per piece the white sphere conserves to 1e-14: the even panel's witness
    at operator level.

    `[M]` the spec's prototype: 1.3e-15 with the even panel, 2.6e-13 with the
    odd one (Lagrange in c on the centre panel), whose b^(2m+2) log b terms
    hold the centre's impact integral to algebraic convergence.
    """
    err = _conservation("sphere", _MR3, (_WHITE,), _SIGMA["g0"], 24)
    assert err <= 1e-14, f"{err:.2e}"


@pytest.mark.l1
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.rests_on(_HERE + "test_a_closed_body_conserves_its_emission")
def test_each_group_is_coupled_through_its_own_cross_section() -> None:
    """[AS5; planned l1] Each group's rule from its own Sigma_t (``LineRule.of(..., sigma_t, ...)``, one rule per
    group since the user's ruling of 2026-10-06 derived the gradings from the group's optical scale): each block
    conserves with its own Sigma_t, and the two couplings differ (the activation leg: a group-blind coupling would be
    one object).

    First red (the sketch's B8 2G mutation, battery arm A2): the first group's
    coupling reused for every later group (`[M]` the spec's prototype: 1.4 on
    the sphere, 0.10 on the slab).
    """
    for chart, bps, laws, points in (("sphere", _MR3, (_WHITE,), 16), ("slab", _SLB3, (_WHITE, _MIRROR), 4)):
        basis = _basis(chart, bps)
        blocks = {}
        for group in ("g0", "g1"):
            g = _transport(chart, bps, laws, _SIGMA[group], points)
            blocks[group] = g
            lhs = g.block @ np.asarray(_SIGMA[group])[basis.region][g.support]
            rhs = basis.mass @ np.ones(basis.size)
            err = np.max(np.abs(lhs - rhs)) / np.max(np.abs(rhs))
            assert err <= _CONSERVATION_TOL[chart], f"{chart} {group}: {err:.2e}"
        assert not np.allclose(blocks["g0"].coupling.escape, blocks["g1"].coupling.escape, rtol=1e-3)


# ── AS2: symmetry, its stabiliser declared ───────────────────────────────


@pytest.mark.foundation
@pytest.mark.parametrize(("chart", "breakpoints", "laws", "points", "group"),
                         _closed_params(_CLOSED)
                         + _closed_params([r for r in _CLOSED_SLOW if r[0] == "cylinder_white"], marks=pytest.mark.slow,
                                          groups=("g0", "g1")))
@pytest.mark.rests_on(_HERE + "test_a_closed_body_conserves_its_emission")
def test_the_block_is_symmetric_and_blind_to_the_line_measure(chart, breakpoints, laws, points, group) -> None:
    """[AS2] K_g = K_g^T to 1e-14 of max|K| (reciprocity of the line Green's function summed over a line's two orientations).

    DECLARED BLIND (the spec's P-2, measured): each line's rule is
    reversal-symmetric, so the row cannot see the impact, mu or polar rule (the
    sphere at 4 impact points: 2.7e-16 where AS1 misses by 1.2e-2), the wall
    area D, the density or its fold, nor the group pairing (U_g0 D^-1 U_g0^T is
    symmetric). Its teeth are a mismatch between the outer and inner arc-length
    rules (`[M]` 1.5e-6 at 4 outer points, 3.8e-7 at 4 inner points).
    """
    K = _transport(chart, breakpoints, laws, _SIGMA[group], points, _INNER).block   # AS1's cached block (same key)
    assert K.shape[0] == K.shape[1]
    assert np.max(np.abs(K - K.T)) <= 1e-14 * np.max(np.abs(K)), np.max(np.abs(K - K.T)) / np.max(np.abs(K))


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_block_is_symmetric_and_blind_to_the_line_measure")
def test_the_symmetry_row_sees_the_arc_length_rules() -> None:
    """[AS2, the teeth leg] At 4 outer points against 12 inner (and 4 inner against 12 outer) the block is asymmetric
    above 1e-10: the symmetry row can see the arc-length rules (`[M]` the spec's prototype 1.5e-6, 3.8e-7)."""
    for outer, inner in ((4, 12), (12, 4)):
        g = _rule("sphere", _MR3, (_WHITE,), _SIGMA["g0"], 16).transport(outer, inner)
        K = g.block
        asym = np.max(np.abs(K - K.T)) / np.max(np.abs(K))
        assert asym > 1e-10, f"outer {outer}, inner {inner}: {asym:.1e}: the row cannot see the arc-length rules"


# ── AS3: the piece budget is a partition of the line sum ─────────────────


@pytest.mark.foundation
@pytest.mark.parametrize(("chart", "breakpoints", "laws", "points"), [
    ("sphere", _MR3, (_WHITE,), 16),
    pytest.param("slab", _SLB3, ((0.3, 0.0), (0.0, 0.8)), 4, marks=pytest.mark.slow),
], ids=["sphere_white", "slab_partial_mirror_white"])
@pytest.mark.rests_on(_HERE + "test_a_closed_body_conserves_its_emission")
def test_the_block_does_not_depend_on_the_piece_budget(chart, breakpoints, laws, points) -> None:
    """[AS3] The block at the default budget (1024), a quarter of it, a budget of 1 (one line per chunk) and an unbounded one
    (one chunk) is one block to 64 eps of max|K|; and every line lands in exactly one chunk.

    Not bitwise: the chunks re-order the sum over lines (`[M]` the spec's P-4: up
    to 12.9 eps on the slab). First reds: a line over the budget dropped instead
    of given its own chunk; a chunk boundary off by one (the line-count leg).
    """
    basis = _basis(chart, breakpoints)
    walls = _walls(chart, breakpoints, laws)
    sigma = np.asarray(_SIGMA["g0"])
    blocks = []
    for budget, chunk in ((1024, 512), (256, 512), (1, 512), (10 ** 12, 10 ** 9)):
        rule = LineRule.of(basis, walls, sigma, points, chunk=chunk, budget=budget)
        counted = sum(len(rule.lines.weights[lines]) for _rule_, lines in rule.lines.chunks(_INNER, _INNER))
        assert counted == len(rule.lines.weights), f"budget {budget}: {counted} lines in the chunks of {len(rule.lines.weights)}"
        blocks.append(rule.transport(_INNER, _INNER).block)
    scale = np.max(np.abs(blocks[0]))
    for k, other in enumerate(blocks[1:], 1):
        assert np.max(np.abs(other - blocks[0])) <= 64 * _EPS * scale, (k, np.max(np.abs(other - blocks[0])) / scale / _EPS)


# ── AS4: no diffuse wall, no update ──────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.parametrize(("chart", "breakpoints", "laws", "points"), [
    ("sphere", _MR3, (_MIRROR,), 16),
    ("sphere", _MR3, ((0.0, 0.0),), 16),
    ("sphere", _MR3H, ((0.3, 0.0), (0.8, 0.0)), 16),
    ("slab", _SLB3, ("P", "P"), 4),
], ids=["sphere_mirror", "sphere_vacuum", "hsphere_partial_mirrors", "slab_periodic"])
@pytest.mark.rests_on(_HERE + "test_a_closed_body_conserves_its_emission")
def test_with_no_diffuse_wall_the_block_is_its_line_part(chart, breakpoints, laws, points) -> None:
    """[AS4] No diffuse wall: the coupling has no wall, its update is exactly 0, and the block is the line part bitwise.

    First red: the coupling reading the specular amplitude as the diffuse one
    (a mirror then carries a full white update).
    """
    g = _transport(chart, breakpoints, laws, _SIGMA["g0"], points)
    assert g.coupling.walls.amplitude.shape == (0,)
    assert g.coupling.currents.shape == (0, g.support.size)
    np.testing.assert_array_equal(g.block, g.line)


# ── WC1, WC2: the escape and the wall transmission against closed forms ──

_HOM = 1.3          # the homogeneous body's radius (or width)


def _homogeneous(chart: str, tau: float, points: int):
    bps = (0.0, _HOM)
    laws = (_WHITE,) if chart != "slab" else (_WHITE, _WHITE)
    g = _transport(chart, bps, laws, (tau / _HOM,), points)
    volume = float(np.sum(Chart(_COORD[chart]).measure(np.asarray(bps))))
    return _basis(chart, bps), g, volume


def _escape_errors(chart: str, tau: float, points: int) -> tuple[float, float, float]:
    """(P_esc from the line part, P_esc from the escape U, T_w) relative to the closed forms."""
    basis, g, volume = _homogeneous(chart, tau, points)
    ones = np.ones(basis.size)
    sigma = tau / _HOM
    with mp.workdps(R.DPS):
        p_esc = float(R.escape_probability(chart, tau))
        p_ss = float(R.transmission_probability(chart, tau))
    by_line = 1.0 - sigma * float(ones @ g.line @ ones[g.support]) / volume
    by_escape = float(np.sum(g.coupling.escape)) / volume
    T = g.coupling.transmission
    t = float(T[-1, 0] if chart == "slab" else T[0, 0])
    if chart == "slab":
        np.testing.assert_array_equal(np.diag(T), 0.0)       # one face never sees itself
    return by_line / p_esc - 1.0, by_escape / p_esc - 1.0, t / p_ss - 1.0


#: tau = 100 (sphere) reads the rim grading; tau = 0.01 (slab, cylinder) the grazing grading derived from the
#: thinnest panel (qa's inputs, 2026-10-06: an ungraded rim 4.6e-3 at tau = 100; 12 fixed halvings 4e-2 to 5e-1 on a
#: near-void slab; plain Gauss in theta 3.2e-4 at tau = 0.01). Tolerances measured at this row's first run.
_ESCAPE = [("sphere", 0.5, 16, 1e-13), ("sphere", 2.0, 16, 1e-13), ("sphere", 8.0, 16, 1e-13),
           ("sphere", 100.0, 16, 1e-11), ("sphere", 1000.0, 16, 1e-9), ("sphere", 0.01, 16, 1e-12),
           ("slab", 0.5, 8, 1e-12), ("slab", 2.0, 8, 1e-12), ("slab", 8.0, 8, 1e-11), ("slab", 30.0, 8, 1e-11),
           ("slab", 0.01, 8, 1e-11)]
#: The rows that catch a catalogued defect, each re-dropped under -O (2026-10-07): ERR-101 at the sphere's thick rim
#: (arms A16, the rim grading removed, and A20, the old ``chord_quadrature`` rule); ERR-103 at the thick slab (arm A22,
#: the normal-direction ends dropped); ERR-101 at the thick cylinder too (arm A16 at 16 points, `[M]` 2026-10-07: T_w
#: 9.4e-7 at tau = 30, 4.2e-4 at tau = 100; arm A22 moves them by 1e-14 and stays green; arm A20's anchor is the
#: sphere's impact rule and does not reach the cylinder).
_ESCAPE_CATCHES = {("sphere", 100.0): pytest.mark.catches("ERR-101"), ("sphere", 1000.0): pytest.mark.catches("ERR-101"),
                   ("slab", 30.0): pytest.mark.catches("ERR-103"),
                   ("cylinder", 30.0): pytest.mark.catches("ERR-101"), ("cylinder", 100.0): pytest.mark.catches("ERR-101")}

#: At 8 points (the slow tier's budget, the user's ruling of 2026-10-07) every leg meets its band except the
#: transmission at tau = 2 and 8 (`[M]` 2026-10-07: T_w 9.3e-10 and 2.3e-8 at 8 points), which stay at 16 points.
_ESCAPE_SLOW = [("cylinder", 0.5, 8, 1e-11), ("cylinder", 2.0, 16, 1e-11), ("cylinder", 0.01, 8, 1e-11),
                ("cylinder", 1e-4, 8, 1e-11), ("cylinder", 8.0, 16, 1e-9),
                ("cylinder", 30.0, 16, 3e-10), ("cylinder", 100.0, 16, 3e-10)]
#: The thick cylinder legs run at 16 points (the user's ruling of 2026-10-07, restoring tau = 30 and 100 after the line
#: rule evaluated only live intervals, #586): one block takes 102 s (tau = 30) and 145 s (tau = 100) at 16 points
#: (`[M]` the orchestrator, one process, budget 1024). Their bands are 10 x `[M]` 2026-10-07 at 16 points: T_w
#: 2.6e-11 at both (P_esc 1.1e-15 and 5.0e-15 at tau = 30; 1.1e-13 and 1.3e-15 at tau = 100).


@pytest.mark.l1
@pytest.mark.verifies("characteristic-boundary-resolvent")
@pytest.mark.parametrize(("chart", "tau", "points", "tol"),
                         [pytest.param(*r, id=f"{r[0]}-tau{r[1]}", marks=_ESCAPE_CATCHES.get((r[0], r[1]), []))
                          for r in _ESCAPE]
                         + [pytest.param(*r, id=f"{r[0]}-tau{r[1]}",
                                         marks=[pytest.mark.slow, *([_ESCAPE_CATCHES[r[:2]]] if r[:2] in _ESCAPE_CATCHES else [])])
                            for r in _ESCAPE_SLOW])
@pytest.mark.rests_on(_HERE + "test_a_closed_body_conserves_its_emission")
def test_the_escape_and_transmission_probabilities_are_the_closed_forms(chart, tau, points, tol) -> None:
    """[WC1, the re-posed C9; WC2; planned l1, ``characteristic-boundary-resolvent``] A homogeneous body behind white
    walls (its line part is the vacuum part): P_esc two ways, 1 - Sigma 1^T K_line 1 / V and 1^T U / V, and the wall
    transmission T_w = P_ss (sphere, cylinder) or 2 E_3 face to face (slab, whose diagonal is 0).

    References (mpmath, ``_characteristic_mp``): Hebert's sphere, (1 - 2 E_3)/(2 tau)
    for the slab, (1 - P_ss)/(2 tau) with Bickley's P_ss for the cylinder; they
    share no line integral with the code. The tolerances are 10 x the
    measurements (`[M]` the orchestrator on the built code: sphere and slab <=
    1.1e-15 at 16 points; cylinder in theta at 32 points 2.3e-13 (tau 0.5) and
    4.2e-15 (tau 2); the prototype's slab T_w at tau = 8: 3.0e-12). First reds:
    Gauss in mu_z on the cylinder (5.4e-4 at 8 points, `[M]` the spec); the cylinder's density without
    |P Omega| (O(1)); the escape summed over the reversed traversals too (x 2 on the slab);
    the arriving flux without 1/D (T_w scales by A/4).
    """
    errors = _escape_errors(chart, tau, points)
    assert max(abs(e) for e in errors) <= tol, f"P_esc(line) {errors[0]:+.1e}, P_esc(U) {errors[1]:+.1e}, T_w {errors[2]:+.1e}"


@pytest.mark.foundation
@pytest.mark.parametrize("chart", ["sphere", "slab"])
@pytest.mark.rests_on(_HERE + "test_the_escape_and_transmission_probabilities_are_the_closed_forms")
def test_a_void_body_transmits_its_geometric_fractions(chart) -> None:
    """[WC2, void legs] With Sigma_t = 0 the transmission is geometry: a hollow sphere's T_w (exit row, entry column) is
    [[0, (r_0/R)^2], [1, 1 - (r_0/R)^2]]; a slab's [[0, 1], [1, 0]]. 4 eps (`[M]` the spec's prototype: 1.0 and 0.5 eps).

    The solid-angle fraction (r_0/R)^2 of a cosine-distributed current hitting a
    concentric sphere, by hand. The support is empty (no region emits), so no
    refusal. First reds: the arriving flux normalised by the EXIT wall's area
    (the off-diagonal moves by A_in / A_out); the arriving flux without the 1/pi
    folded into 1/D (x pi).
    """
    bps = (0.4, 0.5, 1.5, 2.0) if chart == "sphere" else _SLB3
    g = _transport(chart, bps, (_WHITE, _WHITE), (0.0, 0.0, 0.0), 8 if chart == "slab" else 16)
    assert g.support.size == 0
    want = np.array([[0.0, (0.4 / 2.0) ** 2], [1.0, 1.0 - (0.4 / 2.0) ** 2]]) if chart == "sphere" else np.array([[0.0, 1.0], [1.0, 0.0]])
    np.testing.assert_allclose(g.coupling.transmission, want, rtol=0.0, atol=4 * _EPS)


# ── WC3, WC7: reciprocity ────────────────────────────────────────────────


_RECIPROCAL = [("sphere", _MR3H, "g0", 16), ("sphere", _MR3H, "void", 16), ("slab", _SLB3, "g0", 8), ("slab", _SLB3, "void", 8)]


@pytest.mark.l1
@pytest.mark.verifies("characteristic-boundary-resolvent")
@pytest.mark.parametrize(("chart", "breakpoints", "group", "points"), _RECIPROCAL,
                         ids=[f"{r[0]}-{r[2]}" for r in _RECIPROCAL])
@pytest.mark.rests_on(_HERE + "test_a_void_body_transmits_its_geometric_fractions")
def test_the_transmission_is_reciprocal_in_the_wall_areas(chart, breakpoints, group, points) -> None:
    """[WC3; planned l1] diag(A_w)^-1 T_w is symmetric (S_i P_ij = S_j P_ji), both walls white; 8 eps of max.

    A_w is ``Chart.measure_density`` at the wall, read here. `[M]` the spec's
    prototype: 0.5 eps (hollow sphere), 1.6 eps (slab). Declared stabiliser: a
    common scale of T_w passes; WC2 pins the scale. First red: the arriving flux
    normalised by the exit wall's area.
    """
    sigma = (0.0, 0.0, 0.0) if group == "void" else _SIGMA[group]
    g = _transport(chart, breakpoints, (_WHITE, _WHITE), sigma, points)
    area = Chart(_COORD[chart]).measure_density(np.array([breakpoints[0], breakpoints[-1]]))
    S = g.coupling.transmission / area[:, None]
    assert abs(S[0, 1] - S[1, 0]) <= 8 * _EPS * np.max(np.abs(S)), (S[0, 1] - S[1, 0]) / np.max(np.abs(S)) / _EPS


@pytest.mark.l1
@pytest.mark.verifies("characteristic-boundary-resolvent")
@pytest.mark.parametrize(("chart", "breakpoints", "laws", "points"), [
    ("sphere", _MR3, (_WHITE,), 16),
    ("sphere", _MR3H, (_WHITE, (0.0, 0.5)), 16),
    ("slab", _SLB3, (_WHITE, _MIRROR), 8),
], ids=["sphere_white", "hsphere_white_white0.5", "slab_white_mirror"])
@pytest.mark.rests_on(_HERE + "test_the_transmission_is_reciprocal_in_the_wall_areas")
def test_the_walls_response_is_the_escape_over_the_quarter_area(chart, breakpoints, laws, points) -> None:
    """[WC7; planned l1] On the emission support the response of a unit current entering a wall is its escape over
    D = A_w / 4: R[support] = U diag(4 / A_w), 64 eps of max|R| (the orchestrator: <= 6.6e-16 relative).

    R is computed directly (its rows cover every panel); U only on the
    support. The two are the forward and adjoint readings of one line set, so
    the row ties the arriving flux's normalisation to the escape's. First reds:
    the arriving flux without 1/D (R scales by D); R accumulated over the
    reversed traversals too.
    """
    g = _transport(chart, breakpoints, laws, _SIGMA["g0"], points)
    walls = [w for w in _walls(chart, breakpoints, laws).walls if w.diffuse > 0.0]
    area = Chart(_COORD[chart]).measure_density(np.array([breakpoints[0] if w.breakpoint == 0 else breakpoints[-1]
                                                          for w in walls]))
    response = g.coupling.response[g.support]
    want = g.coupling.escape * (4.0 / area)[None, :]
    assert np.max(np.abs(response - want)) <= 64 * _EPS * np.max(np.abs(want)), np.max(np.abs(response - want)) / np.max(np.abs(want))


# ── WC4: the white and the specular laws (B8's operator half) ────────────


@pytest.mark.l1
@pytest.mark.verifies("characteristic-boundary-resolvent")
@pytest.mark.parametrize("alpha", [0.0, 0.5, 1.0])
@pytest.mark.rests_on(_HERE + "test_the_escape_and_transmission_probabilities_are_the_closed_forms")
def test_the_white_and_specular_laws_are_their_closed_forms_and_differ(alpha) -> None:
    """[WC4, B8's operator half; planned l1] 1^T K 1 against closed forms: the white law by the balance of one
    re-emission chain, V/Sigma [(1 - P_esc) + alpha P_esc (1 - T)/(1 - alpha T)] (sphere: T = P_ss; slab with both
    faces white: T = 2 E_3); the specular sphere by the per-line mpmath integral. 1e-13. At alpha = 0.5 the white and
    specular spheres differ by more than 100 x AS1's band.

    `[M]` the spec's prototype: 4.4e-16 to 6.7e-16; the difference at 0.5:
    3.5e-3. The inequality is the discriminator a specular-only code cannot
    pass. First red: the white law realised as specular (the diffuse amplitude
    read into the period's amplitude).
    """
    sigma, R_sphere = 0.7, 2.0
    ones = lambda b: np.ones(b.size)                      # noqa: E731
    sphere = (0.0, 0.8, R_sphere)
    white = _transport("sphere", sphere, ((0.0, alpha),), (sigma, sigma), 16)
    specular = _transport("sphere", sphere, ((alpha, 0.0),), (sigma, sigma), 16)
    b = _basis("sphere", sphere)
    with mp.workdps(R.DPS):
        want_white = float(R.white_total("sphere", sigma, R_sphere, alpha))
        want_specular = float(R.specular_sphere_total(sigma, R_sphere, alpha))
        want_slab = float(R.white_total("slab", sigma, 2.3, alpha))
    got_white = ones(b) @ white.block @ ones(b)
    got_specular = ones(b) @ specular.block @ ones(b)
    assert abs(got_white / want_white - 1.0) <= 1e-13, got_white / want_white - 1.0
    assert abs(got_specular / want_specular - 1.0) <= 1e-13, got_specular / want_specular - 1.0
    slab = (0.0, 0.4, 2.3)
    s = _transport("slab", slab, ((0.0, alpha), (0.0, alpha)), (sigma, sigma), 8)
    bs = _basis("slab", slab)
    got_slab = ones(bs) @ s.block @ ones(bs)
    assert abs(got_slab / want_slab - 1.0) <= 1e-13, got_slab / want_slab - 1.0
    if alpha == 0.5:
        assert abs(got_white - got_specular) / want_specular > 100 * _CONSERVATION_TOL["sphere"]


# ── WC5, WC6: the refusals ───────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.catches("ERR-102")
@pytest.mark.rests_on(_HERE + "test_the_escape_and_transmission_probabilities_are_the_closed_forms",
                      _CLOSURE + "test_the_least_solution_on_a_lossless_trapped_line")
def test_a_source_behind_walls_that_return_everything_in_a_lossless_body_is_refused() -> None:
    """[WC5] ``WallCoupling.currents`` refuses a source reaching walls that all return everything (alpha = 1) around a
    body that loses nothing (every loss 0): ``TrappedSource`` with its own fragment. Each condition dropped once builds.

    The refusal is the coupling's own (the user's ruling of 2026-10-06 moved it
    from an input predicate in ``LineRule.transport`` into the coupling (now ``currents``), beside the
    balance row that makes I - T alpha exactly singular there). The legs:
    alpha = 0.99 builds, finite; Sigma_t > 0 in one region builds; a void body behind a MIRROR (no diffuse wall)
    is refused by the line part's own fragment, not this one (disjoint). First
    red: the check deleted (the trapped body then returns a zero update
    silently, or, with the zero branch also gone, ``LinAlgError``: the wrong
    type).
    """
    void, bps = (0.0, 0.0, 0.0), _MR3
    support = np.array([True, True, True])
    trapped = _rule("sphere", bps, (_WHITE,), void, 8).transport(8, 8, support)
    assert np.all(trapped.coupling.loss == 0.0) and np.any(trapped.coupling.escape != 0.0)
    with pytest.raises(TrappedSource, match="absorbs nothing, behind walls that return everything"):
        trapped.coupling.currents
    assert np.all(np.isfinite(_rule("sphere", bps, ((0.0, 0.99),), void, 8).transport(8, 8, support).block))
    assert np.all(np.isfinite(_rule("sphere", bps, (_WHITE,), (0.0, 0.3, 0.0), 8).transport(8, 8, support).block))
    with pytest.raises(TrappedSource, match="lossless trapped line") as caught:
        _rule("sphere", bps, (_MIRROR,), void, 8).transport(8, 8, support)
    assert "behind walls that return everything" not in str(caught.value)


@pytest.mark.foundation
@pytest.mark.catches("ERR-102")
@pytest.mark.rests_on(_HERE + "test_a_source_behind_walls_that_return_everything_in_a_lossless_body_is_refused")
def test_a_lossless_body_with_no_source_has_an_empty_block() -> None:
    """[WC5b] With no source (an empty emission support: no region emits) the trapped walls carry nothing, and the
    block is the empty (N, 0) matrix, not a refusal and not a singular solve.

    `[M]` 2026-10-06 (this row's first run, on the built code): ``update``
    solves the exactly singular balanced system for zero right-hand columns and
    raises ``LinAlgError``, because its refusal keys on ``escape != 0`` and an
    empty escape skips it. A void body behind an alpha = 1 white wall is the
    default support's case, so any ``.block`` on it fails. Reported to the
    orchestrator, who made ``update`` return the zero update for a zero or
    empty escape on a lossless body (2026-10-06). First red: that zero branch
    falling through to the solve (``LinAlgError``).
    """
    empty = _rule("sphere", _MR3, (_WHITE,), (0.0, 0.0, 0.0), 8).transport(8, 8)
    assert empty.support.size == 0
    assert empty.block.shape == (empty.line.shape[0], 0)


_BALANCE = [
    ("sphere_white", "sphere", _MR3, (_WHITE,), 16),
    ("hsphere_white_white", "sphere", _MR3H, (_WHITE, _WHITE), 16),
    ("hsphere_mirror_white0.5", "sphere", _MR3H, (_MIRROR, (0.0, 0.5)), 16),
    ("slab_white_white", "slab", _SLB3, (_WHITE, _WHITE), 8),
    ("slab_partial0.3_white", "slab", _SLB3, ((0.3, 0.0), _WHITE), 8),
]


@pytest.mark.foundation
@pytest.mark.parametrize("scale", [1.0, 1e-6])
@pytest.mark.parametrize(("chart", "breakpoints", "laws", "points"), [r[1:] for r in _BALANCE],
                         ids=[r[0] for r in _BALANCE])
@pytest.mark.rests_on(_HERE + "test_the_escape_and_transmission_probabilities_are_the_closed_forms")
def test_each_injected_current_is_transmitted_or_lost(chart, breakpoints, laws, points, scale) -> None:
    """[WC8] The current entering at each diffuse wall is accounted for: sum_w' T_w'w + loss_w = 1, 8 eps.

    T is the fraction of the injection reaching each diffuse wall; the loss is
    the fraction absorbed (in (1 - e^{-tau})) plus the fraction leaked at a
    non-diffuse exit (e^{-tau} in (1 - a)). The identity holds per line in exact
    arithmetic (the orchestrator, 2026-10-06), so the row reads the tally's
    bookkeeping, not a quadrature. Bodies with a leaking partial mirror and a
    half-white wall; Sigma_t at the MR3 data and scaled by 1e-6 (near void).
    First reds: the leak at a non-diffuse exit dropped (the partial-mirror
    rows); the absorbed fraction attenuated twice; the tally over the forward
    traversals only (the slab's reversed traversal lost).
    """
    sigma = tuple(scale * v for v in _SIGMA["g0"])
    c = _transport(chart, breakpoints, laws, sigma, points).coupling
    total = c.transmission.sum(axis=0) + c.loss
    np.testing.assert_allclose(total, 1.0, rtol=0.0, atol=8 * _EPS)


#: The near-void sweep (the user's ruling of 2026-10-06 to fix the conditioning): Sigma_t from 1 to 1e-12 on
#: regions (0 or 0.3, 0.6, 1.0). `[M]` the orchestrator's ``band.py`` on the built code: sphere white 5.4e-13 flat;
#: hollow sphere white/white and mirror/white 1.8e-14 flat; slab white/white 1.5e-14 flat. Before the fix, with the
#: diagonal formed by subtraction (1 - alpha T_ww) and no balance row, the hollow sphere reached 1.3e-4 at 1e-12 (`[M]`
#: the archivist's re-run, 2026-10-07); with the loss-formed diagonal but no balance row, 7.7e-6 (hollow sphere) and
#: 2.2e-5 (slab) at 1e-12, and 6.4e-9 and 2.0e-8 at 1e-9.
_NEAR_VOID = [
    ("sphere_white", "sphere", (0.0, 0.6, 1.0), (_WHITE,), 16, 5e-12),
    ("hsphere_white_white", "sphere", (0.3, 0.6, 1.0), (_WHITE, _WHITE), 16, 2e-13),
    ("hsphere_mirror_white", "sphere", (0.3, 0.6, 1.0), (_MIRROR, _WHITE), 16, 2e-13),
    ("slab_white_white", "slab", (0.0, 0.6, 1.0), (_WHITE, _WHITE), 16, 2e-13),
]


@pytest.mark.l1
@pytest.mark.verifies("characteristic-boundary-resolvent")
@pytest.mark.catches("ERR-102")
@pytest.mark.parametrize("sigma", [1.0, 1e-3, 1e-6, 1e-9, 1e-12])
@pytest.mark.parametrize(("chart", "breakpoints", "laws", "points", "tol"), [r[1:] for r in _NEAR_VOID],
                         ids=[r[0] for r in _NEAR_VOID])
@pytest.mark.rests_on(_HERE + "test_each_injected_current_is_transmitted_or_lost",
                      _HERE + "test_a_closed_body_conserves_its_emission")
def test_a_nearly_void_closed_body_conserves_its_emission(chart, breakpoints, laws, points, tol, sigma) -> None:
    """[AS6; planned l1] Conservation K sigma = W 1 on a closed white-walled body whose Sigma_t (two regions, equal)
    falls from 1 to 1e-12: flat in Sigma_t, within 10 x the orchestrator's measurement per body.

    Near void the walls exchange almost every neutron, I - T alpha is nearly
    singular, and its total-current mode is fixed only by what is lost. The
    coupling forms the diagonal from the loss and replaces the last row by the
    balance (1 - alpha) + alpha loss. First reds: the balance row removed
    (battery arm N2: 6 rows, the two-wall bodies at Sigma_t <= 1e-6), and the
    balance row removed together with the diagonal formed as 1 - alpha T_ww
    (arm N12: 12 rows, every body at Sigma_t <= 1e-6). DECLARED: the diagonal
    alone (arm N1) is masked by the balance row, which fixes the total-current
    mode the subtraction would spoil; N1 reds nothing.
    """
    err = _conservation(chart, breakpoints, laws, (sigma, sigma), points)
    assert err <= tol, f"{err:.2e} > {tol:.0e}"


@pytest.mark.foundation
@pytest.mark.parametrize(("specular", "diffuse", "fragment"), [
    (0.4, 0.6, "a wall is specular or diffuse"),
    (1.0, 1.0, "not a physical wall"),
    (0.0, 1.0 + 1e-12, "not a physical wall"),
], ids=["mixed", "returns_twice", "diffuse_above_one"])
def test_a_wall_returning_both_ways_or_more_than_it_receives_is_refused(specular, diffuse, fragment) -> None:
    """[WC6, re-posed by the elegance review: a refusal row] A wall with both a specular and a diffuse part is a
    SCOPE-BOUNDARY (the user's ``LawSum`` ruling, as the SN realizer refuses it); one returning more than it receives
    is not physical. Each refused by its own fragment.

    `[M]` the spec's prototype: the coupling served a 0.4 + 0.6 wall (conservation
    1.2e-12), so the boundary is the reader's choice, not the formula's limit;
    `[M]` the orchestrator: ``Wall(1, 1.0, 1.0, 1)`` assembled a block with a
    minimum of -0.0118 before the guard. First red: either guard deleted.
    """
    with pytest.raises(NotImplementedError, match=fragment):
        Wall(3, specular, diffuse, 3)


# ── OP1-OP4: operator theorems on the material blocks ────────────────────


def _region_transfer(chart, breakpoints, laws, sigma, points=16) -> np.ndarray:
    """T[k, l] = <1_k, K 1_l> over the regions: piecewise-constant q is in every basis exactly."""
    basis = _basis(chart, breakpoints)
    g = _transport(chart, breakpoints, laws, sigma, points)
    region = basis.region
    n = basis.regions.n_regions
    Q = np.stack([(region == k).astype(float) for k in range(n)], axis=-1)
    return Q.T @ g.block @ Q[g.support]


@pytest.mark.l1
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.parametrize("split", [(0.9,), (0.8, 1.2), (0.7, 0.9, 1.1, 1.3)], ids=["2", "3", "5"])
@pytest.mark.rests_on(_HERE + "test_a_closed_body_conserves_its_emission")
def test_an_interface_between_equal_materials_is_invisible(split) -> None:
    """[OP1, C7's operator half; planned l1] Region 1 of MR3 split into 2, 3, 5 equal-material regions: the merged
    region-to-region transfer is unchanged, 1e-13 (`[M]` the spec's prototype: <= 3.7e-16).

    In this architecture an interface is a panel end with a region code, so most
    interface defects are panel-end defects AS1 also reds. First red (owed at
    build, the original spec's named defect): the carried attenuation restarted
    where the region code changes.
    """
    laws = ((0.6, 0.0),)
    base = _region_transfer("sphere", _MR3, laws, _SIGMA["g0"])
    bps = tuple(sorted(_MR3 + split))
    sigma = tuple(0.6 if r < 0.5 else (1.3 if r < 1.5 else 0.45) for r in bps[:-1])
    Tn = _region_transfer("sphere", bps, laws, sigma)
    groups = [0] + [1] * (len(split) + 1) + [2]
    merge = np.zeros((len(groups), 3))
    merge[np.arange(len(groups)), groups] = 1.0
    merged = merge.T @ Tn @ merge
    assert np.max(np.abs(merged - base)) <= 1e-13 * np.max(np.abs(base)), np.max(np.abs(merged - base)) / np.max(np.abs(base))


_CAV = (0.0, 0.4, 0.5, 1.5, 2.0)


@pytest.mark.l1
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.parametrize("a_out", [0.6, 0.0, 1.0])
@pytest.mark.rests_on(_HERE + "test_a_closed_body_conserves_its_emission")
def test_a_transparent_cavity_is_an_inner_mirror(a_out) -> None:
    """[OP2, B4b; planned l1] A solid sphere whose region 0 is void (rank 1 through it) and the hollow sphere with a
    full inner mirror (rank 2 off it) share the region transfer (1e-13) and, on the material panels, the block
    (1e-11; their impact rules differ inside [0, 0.4]).

    `[M]` the spec's prototype: <= 2.5e-16 and <= 2.3e-13. The material panels
    coincide (asserted). First red: the albedo pairing swapped (the amplitude of
    the ENTRY wall): `[M]` 2.1e-3 at a_out = 0.6, 9.0e-3 at 0; BLIND at a_out = 1
    (both walls 1), declared. The cylinder leg is ``slow``-sized and owed.
    """
    sigma_cav, sigma_hol = (0.0, 0.6, 1.3, 0.45), (0.6, 1.3, 0.45)
    cav = _transport("sphere", _CAV, ((a_out, 0.0),), sigma_cav, 16)
    hol = _transport("sphere", _MR3H, (_MIRROR, (a_out, 0.0)), sigma_hol, 16)
    bc, bh = _basis("sphere", _CAV), _basis("sphere", _MR3H)
    material = bc.region > 0
    np.testing.assert_array_equal(bc.nodes[material], bh.nodes)
    np.testing.assert_array_equal(cav.support, np.flatnonzero(material))
    Kc, Kh = cav.block[material], hol.block
    assert np.max(np.abs(Kc - Kh)) <= 1e-11 * np.max(np.abs(Kh)), np.max(np.abs(Kc - Kh)) / np.max(np.abs(Kh))
    Tc = _region_transfer("sphere", _CAV, ((a_out, 0.0),), sigma_cav)[1:, 1:]
    Th = _region_transfer("sphere", _MR3H, (_MIRROR, (a_out, 0.0)), sigma_hol)
    assert np.max(np.abs(Tc - Th)) <= 1e-13 * np.max(np.abs(Th))


_VOIDOUT = (0.0, 0.5, 1.5, 2.0, 2.6)


@pytest.mark.l1
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.parametrize("alpha", [0.6, 1.0])
@pytest.mark.rests_on(_HERE + "test_a_closed_body_conserves_its_emission")
def test_a_void_outer_layer_is_invisible_on_the_emission_support(alpha) -> None:
    """[OP3, B5a; planned l1] A void layer (2.0, 2.6) under a mirror of amplitude alpha at 2.6: on the material panels
    the block (assembled on the emission support, the three material regions) equals MR3's with the mirror at 2.0.
    1e-13 (`[M]` ``ta/m13_support.py``: 1.6e-16 at 0.6, 5.4e-17 at 1).

    At alpha = 1 the lines with b > 2 are lossless and trapped; the support
    keeps their void-panel sources out of the block (the orchestrator's
    decision, 2026-10-06). First reds: the support applied to the columns only
    (the shapes differ); the void panel's Sigma read from the adjacent region.
    """
    g = _transport("sphere", _VOIDOUT, ((alpha, 0.0),), (0.6, 1.3, 0.45, 0.0), 16)
    m = _transport("sphere", _MR3, ((alpha, 0.0),), _SIGMA["g0"], 16)
    bv, bm = _basis("sphere", _VOIDOUT), _basis("sphere", _MR3)
    material = bv.region < 3
    np.testing.assert_array_equal(bv.nodes[material], bm.nodes)
    np.testing.assert_array_equal(g.support, np.flatnonzero(material))
    Kv = g.block[material]
    assert np.max(np.abs(Kv - m.block)) <= 1e-13 * np.max(np.abs(m.block)), np.max(np.abs(Kv - m.block)) / np.max(np.abs(m.block))


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_a_void_outer_layer_is_invisible_on_the_emission_support",
                      _CLOSURE + "test_the_least_solution_on_a_lossless_trapped_line")
def test_the_full_block_of_a_void_layer_under_a_full_mirror_is_refused() -> None:
    """[OP4] The FULL block (support = every region) of the void layer under a mirror at alpha = 1 raises
    ``TrappedSource`` (the line part's fragment); at alpha = 0.6 it builds, finite.

    First red: the support restricted by default whatever the caller asks (the
    refusal never fires).
    """
    sigma = (0.6, 1.3, 0.45, 0.0)
    full = np.ones(4, dtype=bool)
    with pytest.raises(TrappedSource, match="lossless trapped line"):
        _rule("sphere", _VOIDOUT, ((1.0, 0.0),), sigma, 8).transport(8, 8, full)
    g = _rule("sphere", _VOIDOUT, ((0.6, 0.0),), sigma, 8).transport(8, 8, full)
    assert g.support.size == _basis("sphere", _VOIDOUT).size and np.all(np.isfinite(g.block))


# ── EB4, LA1-LA4: self-convergence below the working point ───────────────


def _strictly_converging(errors, floor: float, ratio: float = 10.0) -> None:
    """Monotone; each step above ``floor`` falls by more than ``ratio``."""
    for a, b in zip(errors, errors[1:]):
        assert b < a or a <= floor, f"not monotone: {errors}"
        if b > floor:
            assert a / b > ratio, f"a step falls by {a / b:.1f}: {errors}"


@pytest.mark.l2
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.rests_on(_BASIS + "test_the_even_panel_spans_the_polynomials_in_c_squared")
def test_the_centre_impact_piece_converges_geometrically() -> None:
    """[EB4; planned l2, CONV] The closed sphere's K 1 from the impact piece [0, r_1/2], at n Gauss points against 128:
    below 1e-12 at n = 8, monotone over 4, 6, 8.

    `[M]` the spec's prototype: 1.9e-6, 1.4e-10, 4.6e-14 at 4, 6, 8 with the
    even panel; the odd panel (Lagrange in c) 8.1e-6, 3.3e-7, 3.5e-8 (algebraic:
    the b^(2m+2) log b of the Abel transform of c^(2m+1)). Self-convergence: the
    claim is the rate. First red: the odd panel at the centre.
    """
    bps = (0.0, 0.5, 1.0)
    basis = PanelBasis.of(ConcentricPartition(Chart(CoordSystem.SPHERICAL), bps), 3, 0, 0.5)
    walls = _walls("sphere", bps, (_MIRROR,))
    domain = Chart(CoordSystem.SPHERICAL).line_domain()
    ones = np.ones(basis.size)

    def k_one(n: int) -> np.ndarray:
        q = gauss_legendre(0.0, 0.25, n)
        coordinates = q.pts[:, None]
        # each line chorded from its own b (the level (b, 0))
        rule = LineRule(Lines(basis, walls, np.asarray((0.7, 0.7)), coordinates, q.wts * domain.density(coordinates) / (4.0 * math.pi),
                              (q.pts, np.zeros_like(q.pts)), 512, 10 ** 9))
        return rule.transport(12, 12).block @ ones

    reference = k_one(128)
    errors = [float(np.max(np.abs(k_one(n) - reference))) for n in (4, 6, 8)]
    assert errors[-1] <= 1e-12, errors
    _strictly_converging(errors, floor=1e-14)


@pytest.mark.l2
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.rests_on(_HERE + "test_a_closed_body_conserves_its_emission")
def test_conservation_converges_in_the_impact_rule_below_the_working_point() -> None:
    """[LA1; planned l2, CONV] The white sphere's conservation at 4, 6, 8, 12, 16 impact points per piece: monotone,
    each step down by more than 10, 16 within AS1's band.

    `[M]` the spec's prototype: 1.2e-2, 1.8e-4, 3.4e-6, 1.8e-9, 1.2e-12. First red
    (owed at build): the visibility-cone substitution dropped from the impact
    pieces (algebraic).
    """
    errors = [_conservation("sphere", _MR3, (_WHITE,), _SIGMA["g0"], n) for n in (4, 6, 8, 12, 16)]
    assert errors[-1] <= _CONSERVATION_TOL["sphere"], errors
    _strictly_converging(errors, floor=1e-13)


@pytest.mark.l2
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.rests_on(_HERE + "test_a_closed_body_conserves_its_emission")
def test_conservation_converges_in_the_arc_length_rules_below_the_working_point() -> None:
    """[LA2; planned l2, CONV] Outer points 4, 6, 8 (inner 12) and inner points 4, 6, 8 (outer 12): monotone, each step
    down by more than 10 above the impact rule's floor.

    `[M]` the spec's prototype: outer 6.3e-7, 2.7e-10, 1.3e-12; inner 6.9e-8,
    1.5e-11, 1.2e-12. First red (owed at build): the inner rule's exponential
    grading removed.
    """
    rule = _rule("sphere", _MR3, (_WHITE,), _SIGMA["g0"], 16)
    basis = rule.lines.basis
    rhs = basis.mass @ np.ones(basis.size)
    sigma = np.asarray(_SIGMA["g0"])

    def err(outer, inner):
        g = rule.transport(outer, inner)
        return float(np.max(np.abs(g.block @ sigma[basis.region][g.support] - rhs)) / np.max(np.abs(rhs)))

    _strictly_converging([err(n, 12) for n in (4, 6, 8)], floor=1e-11)
    _strictly_converging([err(12, n) for n in (4, 6, 8)], floor=1e-11)


def _slab_rule_with_fixed_halvings(breakpoints, laws, sigma, points: int, halvings: int) -> LineRule:
    """The slab's rule with the cosine halved a FIXED number of times toward 0 (the rule before the 2026-10-06
    ruling), built here through the public constructor: the derived grading's comparison."""
    basis = _basis("slab", breakpoints)
    half = composite_gauss_legendre([0.0, *graded_ends(0.0, 1.0, True, False, halvings, 0.5), 1.0], points)
    mu = np.concatenate([-half.pts[::-1], half.pts])[:, None]
    w = np.concatenate([half.wts[::-1], half.wts])
    domain = Chart(CoordSystem.CARTESIAN).line_domain()
    return LineRule(Lines(basis, _walls("slab", breakpoints, laws), np.asarray(sigma), mu, w * domain.density(mu) / (4.0 * math.pi),
                          None, 512, 1024))


@pytest.mark.l2
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.catches("ERR-101")
@pytest.mark.rests_on(_HERE + "test_the_escape_and_transmission_probabilities_are_the_closed_forms")
def test_the_slabs_derived_grading_converges_and_beats_twelve_fixed_halvings() -> None:
    """[LA3; planned l2, CONV] The slab's cosine rule, graded from the thinnest panel's optical width: (i) the escape
    at tau = 0.5 and the transmission at tau = 8 fall monotonically over 2, 3, 4, 6 points, by more than 10 a step;
    (ii) on a near-void slab (Sigma = 1e-9) the collision probability from the line part meets the closed form to
    1e-11 while the rule of 12 fixed halvings (built here) misses it by more than 1e-3.

    Conservation is blind to this rule (the spec's P-3). `[M]` the spec's
    prototype for (i): P_esc 1.2e-4, 1.4e-6, 9.5e-9, 1.0e-11; T_w 5.6e-2, 5.5e-4,
    1.2e-5, 2.4e-9. For (ii) qa's ``p13`` measured 4e-2 to 5e-1 with 12 fixed
    halvings. First red: the derived depth replaced by 12 halvings (battery arm A18).
    """
    escape = [abs(_escape_errors("slab", 0.5, n)[0]) for n in (2, 3, 4, 6)]
    transmission = [abs(_escape_errors("slab", 8.0, n)[2]) for n in (2, 3, 4, 6)]
    _strictly_converging(escape, floor=1e-13)
    _strictly_converging(transmission, floor=1e-11)
    sigma, L = 1e-9, 2.3
    bps, laws = (0.0, 0.9, L), ((0.0, 0.0), (0.0, 0.0))
    with mp.workdps(60):
        p_c = float(1 - R.escape_probability("slab", sigma * L))

    def collision(rule: LineRule) -> float:
        g = rule.transport(_INNER, _INNER)
        ones = np.ones(rule.lines.basis.size)
        return sigma * float(ones @ g.line @ ones[g.support]) / L / p_c - 1.0

    derived = collision(_rule("slab", bps, laws, (sigma, sigma), 8))
    fixed = collision(_slab_rule_with_fixed_halvings(bps, laws, (sigma, sigma), 8, 12))
    assert abs(derived) <= 1e-11, f"derived grading: {derived:+.1e}"
    assert abs(fixed) > 1e-3, f"12 fixed halvings reach {fixed:+.1e}: the row cannot tell the gradings apart"


# ── QA1-QA4: the input regions qa found the fixed rules blind to (2026-10-06) ──

_NEAR_VOID_SLAB = [1e-4, 1e-6, 1e-9, 1e-12]


@pytest.mark.l1
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.parametrize("sigma", [pytest.param(v, marks=pytest.mark.catches("ERR-101")) if v <= 1e-6 else v
                                   for v in _NEAR_VOID_SLAB])  # 1e-4 is green under the defect (arm A18)
@pytest.mark.rests_on(_HERE + "test_the_slabs_derived_grading_converges_and_beats_twelve_fixed_halvings")
def test_a_near_void_slab_with_open_walls_meets_its_closed_forms(sigma) -> None:
    """[QA1; planned l1] A slab (0, 0.9, 2.3) of Sigma from 1e-4 to 1e-12, its walls NOT both returning everything:
    vacuum (the collision probability from the line part, 1 - P_esc) and white of albedos 1 and 0.5 (1^T K 1 against
    the balance of the two faces' re-emission chains, ``R.slab_two_white_total``). 1e-11.

    `[M]` 2026-10-06 on the built code (``gates/new_rows_probe.py slab_nearvoid``):
    vacuum <= 7.5e-13, white <= 5.5e-13. Near void a grazing line is still
    optically thick only for |mu| ~ tau, so the cosine must be graded to the
    body's own optical width. First red: 12 fixed halvings (arm A18; qa's ``p12``
    and ``p13``: 4e-2 to 5e-1).
    """
    L = 2.3
    with mp.workdps(60):
        p_c = float(1 - R.escape_probability("slab", sigma * L))
        total = float(R.slab_two_white_total(sigma, L, 1.0, 0.5))
    vac = _transport("slab", (0.0, 0.9, L), ((0.0, 0.0), (0.0, 0.0)), (sigma, sigma), 8)
    ones = np.ones(_basis("slab", (0.0, 0.9, L)).size)
    got_c = sigma * float(ones @ vac.line @ ones) / L
    assert abs(got_c / p_c - 1.0) <= 1e-11, got_c / p_c - 1.0
    white = _transport("slab", (0.0, 0.9, L), ((0.0, 1.0), (0.0, 0.5)), (sigma, sigma), 8)
    got = float(ones @ white.block @ ones)
    assert abs(got / total - 1.0) <= 1e-11, got / total - 1.0


@pytest.mark.slow
@pytest.mark.l1
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.parametrize("sigma", [1.0, 1e-2])
@pytest.mark.rests_on(_HERE + "test_the_slabs_derived_grading_converges_and_beats_twelve_fixed_halvings")
def test_thin_graded_panels_set_the_slabs_grazing_depth(sigma) -> None:
    """[QA2; planned l1] A slab basis graded 10 layers deep (thinnest panel 4.7e-5): the collision probability from the
    line part against the closed form, 1e-11 (`[M]` 1.3e-15 at Sigma = 1, 2.6e-13 at 1e-2).

    The thinnest absorbing PANEL, not the body, sets the grazing depth: a line
    grazing a thin panel is optically thick in it at |mu| ~ its width. First
    red: 12 fixed halvings (arm A18; qa's ``p14``).
    """
    bps, L = (0.0, 0.9, 2.3), 2.3
    basis = PanelBasis.of(ConcentricPartition(Chart(CoordSystem.CARTESIAN), bps), 3, 10, 0.4)
    assert np.diff(basis.partition.breakpoints).min() < 1e-4
    rule = LineRule.of(basis, _walls("slab", bps, ((0.0, 0.0), (0.0, 0.0))), np.asarray((sigma, sigma)), 8)
    g = rule.transport(_INNER, _INNER)
    ones = np.ones(basis.size)
    with mp.workdps(40):
        p_c = float(1 - R.escape_probability("slab", sigma * L))
    got = sigma * float(ones @ g.line @ ones) / L
    assert abs(got / p_c - 1.0) <= 1e-11, got / p_c - 1.0


@pytest.mark.slow
@pytest.mark.l1
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.rests_on(_HERE + "test_thin_graded_panels_set_the_slabs_grazing_depth")
def test_every_block_entry_of_a_thin_panel_slab_is_resolved_at_eight_points() -> None:
    """[QA3; planned l1] A slab basis graded 6 layers deep at Sigma = 0.01: every diagonal entry of the block at
    8 points against the same block with 30 more halvings toward grazing at 16 points (built here through the public
    constructor), 1e-11 relative.

    A total cannot see this: totals were exact to 5e-16 when the grazing
    halving stopped at the thinnest panel's width tau_min, while one entry was
    off by 2.0e-6 (qa's ``q3``). A panel's attenuation e^{-tau_P s} changes until
    tau_P s reaches the vanishing depth 64, so the halving reaches tau_min / 64.
    First red: the halving stopped at tau_min (battery arm A23; production's
    ``VANISHING_DEPTH`` -> 1).
    """
    bps, sigma = (0.0, 0.9, 2.3), 0.01
    basis = PanelBasis.of(ConcentricPartition(Chart(CoordSystem.CARTESIAN), bps), 3, 6, 0.4)
    walls = _walls("slab", bps, ((0.0, 0.0), (0.0, 0.0)))
    honest = LineRule.of(basis, walls, np.asarray((sigma, sigma)), 8).transport(_INNER, _INNER).block
    thinnest = float(np.min(sigma * np.diff(basis.partition.breakpoints)))
    deep = int(np.ceil(np.log2(1.0 / thinnest))) + 30
    half = composite_gauss_legendre([0.0, *graded_ends(0.0, 1.0, True, False, deep, 0.5), 1.0], 16)
    mu = np.concatenate([-half.pts[::-1], half.pts])[:, None]
    w = np.concatenate([half.wts[::-1], half.wts])
    domain = Chart(CoordSystem.CARTESIAN).line_domain()
    reference = LineRule(Lines(basis, walls, np.asarray((sigma, sigma)), mu, w * domain.density(mu) / (4.0 * math.pi),
                               None, 512, 1024)).transport(_INNER, _INNER).block
    worst = float(np.max(np.abs(np.diag(honest) - np.diag(reference)) / np.abs(np.diag(reference))))
    assert worst <= 1e-11, f"worst diagonal entry {worst:.1e}"


@pytest.mark.l1
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.rests_on(_HERE + "test_a_closed_body_conserves_its_emission")
def test_a_small_first_region_conserves_at_eight_points() -> None:
    """[QA5; planned l1] A sphere whose first region has radius 1e-3, (0, 1e-3, 0.4, 1.0), white outer wall, at 8
    points: closed-body conservation to 1e-9 (`[M]` the orchestrator on the final code: 8.8e-11).

    The wide middle impact panels see the NEXT radius's branch point, imaginary
    in the chord half-length y and close to the panel; the y-rule halves toward
    it. The row was written as that grading's witness, and on the final code it
    is NOT one: DECLARED BLIND to arm A19 (``beyond = inf``). `[M]` 2026-10-07
    (``gates/a19_probe.py``, 12/12 inner points): honest 3.5e-12 at 8 points and
    1.4e-15 at 16; with the grading removed 8.1e-12 and 1.3e-15. A tolerance
    between them would be a knife edge, not a gate. The 1.3e-4 that motivated
    the grading was measured on the rule before the y-space substitution, which
    absorbs most of the branch point. The row stays as the value witness of the
    region (the old ``chord_quadrature`` rule, arm A20, and the other line-rule
    arms are its first reds); A19's witness is
    ``test_a_wide_panel_before_a_thin_one_conserves_at_eight_points``.
    """
    err = _conservation("sphere", (0.0, 1e-3, 0.4, 1.0), (_WHITE,), (0.7, 0.9, 1.1), 8)
    assert err <= 1e-9, f"{err:.2e}"


@pytest.mark.l1
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.catches("ERR-101")
@pytest.mark.rests_on(_HERE + "test_a_closed_body_conserves_its_emission")
def test_a_wide_panel_before_a_thin_one_conserves_at_eight_points() -> None:
    """[QA6; planned l1] A hollow sphere (0.2, 0.3, 1.0) on a basis graded 10 layers deep (ratio 0.4), inner mirror,
    outer white, at 8 points: closed-body conservation to 1e-11.

    Fine interface grading puts a wide panel just inside a thin one, so the
    wide panel's chord half-length y sees the next radius's branch point close
    by. `[M]` 2026-10-07 (``gates/a19_wide.py``, 12/12 inner points, MR3 g0
    data and three other Sigma_t pairs): 1.0e-12 with the hp grading toward the
    next radius, 2.0e-9 without it (the orchestrator measured 9.2e-11 and 1.7e-3
    on his own fixture). Tolerance 10 x the honest value, so the arm clears it
    200-fold. First red: that grading removed (battery arm A19, ``beyond = inf``);
    the first-region row above is blind to it.
    """
    bps, sigma = (0.2, 0.3, 1.0), np.asarray(_SIGMA["g0"][:2])
    basis = PanelBasis.of(ConcentricPartition(Chart(CoordSystem.SPHERICAL), bps), 3, 10, 0.4)
    g = LineRule.of(basis, _walls("sphere", bps, (_MIRROR, _WHITE)), sigma, 8).transport(_INNER, _INNER)
    lhs = g.block @ sigma[basis.region][g.support]
    rhs = basis.mass @ np.ones(basis.size)
    err = float(np.max(np.abs(lhs - rhs)) / np.max(np.abs(rhs)))
    assert err <= 1e-11, f"{err:.2e}"


@pytest.mark.l1
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.catches("ERR-101")
@pytest.mark.rests_on(_HERE + "test_a_closed_body_conserves_its_emission")
def test_a_tiny_cavity_assembles_and_conserves() -> None:
    """[QA4; planned l1] A hollow sphere with a cavity of radius 1e-12 (mirror inside and out) assembles and conserves
    to AS1's band.

    Before the 2026-10-06 fix the hp distance toward b = 0 was formed as
    r_{k+1} - span, which cancels to 0 for a tiny r_k, and ``_halvings``
    refused a zero distance (qa's N2). First reds: the subtraction ``hi - span``
    (arm A24: the rule refuses); b computed from hi^2 - y^2 (arm A25).
    """
    err = _conservation("sphere", (1e-12, 0.5, 1.5, 2.0), (_MIRROR, _MIRROR), _SIGMA["g0"], 16)
    assert err <= _CONSERVATION_TOL["sphere"], f"{err:.2e}"


@pytest.mark.slow
@pytest.mark.l2
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.rests_on(_HERE + "test_the_cylinders_direction_rule_is_gauss_in_the_polar_angle")
def test_the_cylinders_escape_converges_in_the_polar_angle() -> None:
    """[LA4; planned l2, CONV] The cylinder's escape at tau = 0.5 with 4, 6, 8 points in theta (and in b): monotone,
    each step down by more than 10 (the ladder up to 8, the working point of WC1's tau = 0.5 cylinder row).

    `[M]` the orchestrator on the built code: 1.5e-6, 2.3e-9 at 8, 16 (32: 2.3e-13).
    First red: Gauss in mu_z (`[M]` 4.3e-3, 5.4e-4, 7.0e-5 at 4, 8, 16: a ratio
    of 8 per doubling).
    """
    errors = [abs(_escape_errors("cylinder", 0.5, n)[0]) for n in (4, 6, 8)]
    _strictly_converging(errors, floor=1e-12)
