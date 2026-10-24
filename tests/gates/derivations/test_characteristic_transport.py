"""Gates for the transport along one line on the panel basis (:class:`~orpheus.derivations.continuous.characteristic.transport.TraversalRule`).

P1 step (b), second rung, of the characteristic-reference campaign
(``.claude/plans/characteristic_reference_architecture.md``, "P1 step (b),
second rung: API sketch", items 2 and 3, and its "Ruled 2026-10-06" line).
Re-specification: ``scratch/characteristic_architecture/p1_step_b2/spec.md``,
rows T0-T12 (the spec's A2 B half, B5c, C2c, C3, the turning-map analogue of
C10, the reciprocity identities, the Volterra triangle and the line-level
vacuum edge).

The ladder, bottom up (``rests_on`` on each row):

1. the kernel's chord and transits; the basis (``test_characteristic_basis.py``);
   the period and its least solution (``test_characteristic_closure.py``);
2. the period on the panel chord is the period on the body's chord [T0];
3. the outflow B_k against mpmath line integrals of per-region polynomials with
   a Sigma_t jump at every interface, per traversal, both ranks [T1, the spec's
   A2 B half and C3], and per basis function [T2]; the entry response A_k
   [T3]; their design identity A_fwd == B_rev [T4];
4. the edges of the arc-length rule: a 1e3-mfp slot [T7, C2c], a void region
   [T8, B5c], a turning slot with b small against its panel [T9];
5. the Volterra triangle against an mpmath double integral [T5], the reversed
   line's triangle is the transpose [T5b], the batch's weighted sum [T6];
6. the angular flux along the line: vacuum against mpmath [T10]; with its
   closure, against the explicitly unfolded backward path [T11, B2's line
   half]; amplitude 0 adds nothing, bitwise [T12, B4a at line level].

**Scope after P3 step 3 (the cylinder's axial lift).** The cylinder rows here
gate ``TraversalRule`` on 3-D cylinder chords at speed 1, the path
``angular_flux`` keeps (the user's ruling of 2026-10-11). The cylinder's Galerkin
block and point reading no longer run on that path: they read each in-plane
line at the family of speeds 1/sin(theta), gated against the same mpmath lines
in ``test_characteristic_lift.py`` (K1a-K1d, K2, K3). A cylinder row here is
therefore evidence for ``angular_flux`` and the per-line transport, not for the
cylinder's block or reading.

Every reference is written in ``_characteristic_mp.py`` (mpmath, no import of
the code under test). Every traversal a row reads is HAND-WRITTEN per fixture
(transit, reversed), never taken from the period under test.

Levels: T1, T2, T3, T7, T8 and T9 are ``l0`` with
``verifies("characteristic-traversal-integrals")`` (tau_k and B_k as one
equation on ``docs/theory/references/characteristic.rst``); T11 is ``l0`` with
``verifies("characteristic-closure")``. Each was moved with its witness: a
mutation of the integral the label names reddens it (`[M]` 2026-10-06, battery
arm T1, every row of the six functions; arm T11, 4 of 4 rows of T11). Every
other row is ``foundation``.
"""
from __future__ import annotations

import mpmath as mp
import numpy as np
import pytest

from orpheus.derivations.continuous.characteristic import LinePeriod, Walls
from orpheus.derivations.continuous.characteristic.basis import PanelBasis
from orpheus.derivations.continuous.characteristic.transport import TraversalRule
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.line import Line

from tests.gates.derivations import _posed_walls as walls_of

from . import _characteristic_mp as R

#: A floating-point warning is an error in every row: a masked 0/0 or inf*0 still raises under -W error
#: (rung 1's lesson), and a reference must not warn on its singular stratum (b = 0) or a void.
pytestmark = pytest.mark.filterwarnings("error::RuntimeWarning")

_HERE = "tests/gates/derivations/test_characteristic_transport.py::"
_BASIS = "tests/gates/derivations/test_characteristic_basis.py::"
_CLOSURE = "tests/gates/derivations/test_characteristic_closure.py::"
_TRANSITS = "tests/gates/geometry/test_chord_transits.py::"

_EPS = float(np.finfo(float).eps)

_MR3 = (0.0, 0.5, 1.5, 2.0)
_MR3H = (0.4, 0.5, 1.5, 2.0)
_SLB3 = (0.0, 0.4, 1.5, 2.3)
#: Spec §4: Sigma_t per region and group, distinct at every interface and per group.
_SIGMA = {0: (0.6, 1.3, 0.45), 1: (1.7, 0.35, 2.4)}
#: Per-region source polynomials in the orbit coordinate, degree <= 3, jumping at every interface.
_Q = [[1.0, 0.3, -0.2, 0.05], [0.25, 1.0, 0.5, -0.1], [3.0, -1.0, 0.1, 0.02]]
_Q_LOW = [[1.0, 0.3], [0.25, 1.0], [3.0, -1.0]]
#: ``_Q`` with region 0 even (rung 3, the spec's RB4): on a solid body region 0 holds the even panel at the
#: stratum (Lagrange in c^2) beside ordinary cubic panels, so the region's exactness family is {1, c^2}.
_Q_SOLID = [[1.0, 0.0, -0.2, 0.0], _Q[1], _Q[2]]


def _q(body: str):
    """The per-region source of the T1 family on ``body``: ``_Q``, with region 0 even on a solid body."""
    return _Q_SOLID if body.endswith("_solid") else _Q
_A_IN, _A_OUT = 0.3, 0.6
_RES = (3, 2, 0.5)                 # (degree, layers, ratio): the rung's working point


def _spec(a: float):
    """The wall returning the fraction ``a`` specularly."""
    return walls_of.partial(a)


def _geometry(body: str, *, a_in: float = _A_IN, a_out: float = _A_OUT, breakpoints=None) -> walls_of.Posed:
    match body:
        case "sphere_solid":
            return walls_of.sphere(breakpoints or _MR3, (0, 1, 2), outer=_spec(a_out))
        case "cylinder_solid":
            return walls_of.cylinder(breakpoints or _MR3, (0, 1, 2), outer=_spec(a_out))
        case "sphere_hollow":
            return walls_of.sphere(breakpoints or _MR3H, (0, 1, 2), inner=_spec(a_in), outer=_spec(a_out))
        case "cylinder_hollow":
            return walls_of.cylinder(breakpoints or _MR3H, (0, 1, 2), inner=_spec(a_in), outer=_spec(a_out))
        case "slab":
            return walls_of.slab(breakpoints or _SLB3, (0, 1, 2), left=_spec(a_in), right=_spec(a_out))
        case "slab_half":                  # B6 (ERR-035): equal intermediate albedos, the symmetric partial slab
            return walls_of.slab(breakpoints or _SLB3, (0, 1, 2), left=_spec(0.5), right=_spec(0.5))
        case "slab_periodic":
            return walls_of.slab(breakpoints or _SLB3, (0, 1, 2), left=walls_of.WRAP,
                                           right=walls_of.WRAP)
    raise AssertionError(body)


def _chart(body: str) -> str:
    return body.split("_")[0]


def _in_plane(b: float, wz: float = 0.0):
    """The point (b, 0, 0) and a direction with in-plane part along y: the impact parameter is b, the foot (b, 0, 0)."""
    return (b, 0.0, 0.0), (0.0, float(np.sqrt((1.0 - wz) * (1.0 + wz))), wz)


def _slab_dir(mu: float):
    return (0.7, 0.0, 0.0), (mu, float(np.sqrt((1.0 - mu) * (1.0 + mu))), 0.0)


class _Case:
    """One line through one body, built twice: by the code under test and in mpmath."""

    def __init__(self, body: str, point, direction, *, sigma_t, resolution=_RES, geometry=None,
                 ) -> None:
        self.body = body
        posed = geometry or _geometry(body)
        self.geometry, self.exterior = posed
        self.basis = PanelBasis.of(ConcentricPartition.of(self.geometry), *resolution)
        self.line = Line.through(np.asarray(point, dtype=float), np.asarray(direction, dtype=float))
        self.walls = Walls.of(self.geometry, self.exterior)
        self.sigma_t = np.asarray(sigma_t, dtype=float)
        self.rule = TraversalRule.of(self.line, self.basis, self.walls, self.sigma_t)
        self.period = self.rule.period
        self.mp = R.MpLine(_chart(body), self.geometry.breakpoints, point, direction,
                           extra_radii=self.basis.partition.breakpoints)

    def coefficients(self, q) -> np.ndarray:
        return R.sample(R.region_polynomial(q), np.asarray(self.basis.region), np.asarray(self.basis.nodes))

    def read_at(self, t) -> tuple[TraversalRule, float]:
        """The rule on the line chorded through the basis's panels refined at the level of the point at parameter
        ``t`` (``ConcentricPartition.refined_at``: the flux is read at a crossing, since P3 step 3's target-graded reading), and
        that chord's own crossing nearest ``t``, where the point is read and the reference is evaluated."""
        level = float(self.mp.c(R.mpf(t) if not isinstance(t, mp.mpf) else t))
        chord = self.basis.partition.refined_at(level).chord(self.line)
        rule = TraversalRule.on_chord(chord, self.basis, self.walls, self.sigma_t)
        crossings = rule.period.chord.crossings
        parameters = np.asarray(crossings.parameter, dtype=float)[np.asarray(crossings.present, dtype=bool)]
        return rule, float(parameters[np.argmin(np.abs(parameters - float(t)))])


def _close(got: float, want, rtol: float, what: str) -> None:
    want = float(want)
    assert abs(got - want) <= rtol * abs(want), f"{what}: got {got!r}, mpmath {want!r}, rel {abs(got - want) / abs(want):.2e}"


# ── T0: the period on the panel chord ────────────────────────────────────

#: (id, body, point, direction, the period by hand: [(transit, reversed)], exit wall by hand: "first" / "last")
_LINES = [
    ("sphere_solid_b0", "sphere_solid", *_in_plane(0.0), [(0, False, "last")]),
    ("sphere_solid_b0.3", "sphere_solid", *_in_plane(0.3), [(0, False, "last")]),
    ("sphere_solid_b0.7", "sphere_solid", *_in_plane(0.7), [(0, False, "last")]),
    ("sphere_solid_b1.7", "sphere_solid", *_in_plane(1.7), [(0, False, "last")]),
    ("cylinder_solid_b0.7_wz0.8", "cylinder_solid", *_in_plane(0.7, 0.8), [(0, False, "last")]),
    ("cylinder_solid_b0.3_wz-0.6", "cylinder_solid", *_in_plane(0.3, -0.6), [(0, False, "last")]),
    ("sphere_hollow_b0.2_cavity", "sphere_hollow", *_in_plane(0.2), [(0, False, "first"), (1, False, "last")]),
    ("sphere_hollow_b0.45_shell", "sphere_hollow", *_in_plane(0.45), [(0, False, "last")]),
    ("cylinder_hollow_b0.2_wz0.8", "cylinder_hollow", *_in_plane(0.2, 0.8), [(0, False, "first"), (1, False, "last")]),
    ("slab_rising", "slab", *_slab_dir(0.6), [(0, False, "last"), (0, True, "first")]),
    ("slab_falling", "slab", (0.7, 0.0, 0.0), (-0.3, 0.0, float(np.sqrt(0.91))),
     [(0, False, "first"), (0, True, "last")]),
    ("slab_grazing_mu0.05", "slab", *_slab_dir(0.05), [(0, False, "last"), (0, True, "first")]),
]
_LINE_IDS = [r[0] for r in _LINES]
_LINE_ARGS = [r[1:] for r in _LINES]


@pytest.mark.foundation
@pytest.mark.parametrize(("body", "point", "direction", "expected"), _LINE_ARGS, ids=_LINE_IDS)
@pytest.mark.rests_on(_CLOSURE + "test_the_period_matches_the_hand_counted_table",
                      _BASIS + "test_the_walls_rekey_onto_the_panel_partition")
def test_the_panel_chord_unfolds_into_the_hand_counted_period(body, point, direction, expected) -> None:
    """[T0] On the panel chord the period is the hand-counted one, and tau_k is the body's.

    Rank, transit, reversal and exit wall (0, or the panel count P for the
    last) are written by hand per line; the traversals' optical depth through
    the per-panel Sigma equals the body chord's ``LinePeriod.optical_depth`` of
    the per-region Sigma_t (a refinement adds interfaces between equal
    materials, which change nothing: the spec's C7 at line level), 16 ulp.
    First reds: ``Walls.on`` keeping n as the last index (the outer wall lost,
    a ``RuntimeError`` from the period); Sigma per panel read through the panel
    index as a region index.
    """
    case = _Case(body, point, direction, sigma_t=_SIGMA[0])
    period, P = case.period, case.basis.n_panels
    m = len(expected)
    assert int(period.rank) == m
    got = [(int(period.transit[k]), bool(period.reversed[k]), int(period.exit_wall[k])) for k in range(m)]
    want = [(t, r, 0 if end == "first" else P) for t, r, end in expected]
    assert got == want
    coarse = LinePeriod.of(ConcentricPartition.of(case.geometry).chord(
        Line.through(np.asarray(point, dtype=float), np.asarray(direction, dtype=float))), Walls.of(case.geometry, case.exterior))
    np.testing.assert_allclose(case.rule.optical_depth, coarse.optical_depth(case.sigma_t), rtol=16 * _EPS, atol=0.0)


# ── T1-T3: the traversal integrals ───────────────────────────────────────


def _traversals(case: _Case, expected):
    return [(k, case.mp.traversal(transit, rev)) for k, (transit, rev, _end) in enumerate(expected)]


_T1_TOL = 1e-13


@pytest.mark.l0
@pytest.mark.verifies("characteristic-traversal-integrals", "peierls-greens-mr-regionwise-source")
@pytest.mark.catches("ERR-090")
@pytest.mark.parametrize("group", [0, 1])
@pytest.mark.parametrize(("body", "point", "direction", "expected"), _LINE_ARGS, ids=_LINE_IDS)
@pytest.mark.rests_on(_HERE + "test_the_panel_chord_unfolds_into_the_hand_counted_period",
                      _BASIS + "test_the_interpolant_reproduces_every_per_region_polynomial_of_degree_p")
def test_the_outflow_of_a_per_region_polynomial_is_its_line_integral_attenuated_to_the_exit(
        body, point, direction, expected, group) -> None:
    """[T1, the spec's A2 B half and C3] B_k . q = int_k q(c(s)) e^{-tau(s -> exit_k)} ds, against mpmath.

    q is a cubic in the orbit coordinate per region (the basis's exactness
    family at p = 3, so the basis error is null), jumping at both interfaces;
    Sigma_t jumps too, distinct per group. Lines of both ranks on every chart:
    through the centre, turning in each region, oblique cylinder lines, a
    shell's ray and a ray through the cavity, the slab both ways and grazing
    (mu = 0.05, about 70 mfp). Rank 2 reads each traversal separately, the
    slab's reversed traversal included. 1e-13 relative.

    First reds: (a) B attenuated from the closest approach or the entry
    instead of to the exit (regime 8: invisible on a palindromic homogeneous
    chord, red on every multi-region row); (b) Sigma of one panel read for the
    whole slot sequence; (c) the turning grading removed (green here; T9 is its witness); (d) the
    exponential grading removed is GREEN here, declared (`[M]` battery arm T2:
    the grazing row's middle panels are within 16 points' reach; T7 is its
    witness).

    Re-posed 2026-10-06 (rung 3, the spec's RB4): on a solid body region 0's
    cubic is even ({1, c^2}, ``_Q_SOLID``), the family the even panel at the
    centre shares with the region's other panels.

    Step (e1b) made it the successor of the old family's
    ``test_mr_oracle_first_leg_matches_the_line_integral`` (the spec's C3):
    each segment reads its own region's piece of the emission
    (``peierls-greens-mr-regionwise-source``). ERR-090's class (one piece
    across an interface) is unspellable on the panel basis: re-dropped as
    interfaces that are not panel ends (arm ``interfaces-not-panel-ends``),
    every row reds by the basis's refusal "the panel ends refine the body's
    breakpoints" (``[M]`` 2026-10-10, 24 of 24). Its value form keeps the
    panel ends and reads the first panel of each region with its
    neighbour's piece (arm ``err090-panel-region-shifted``): 22 of 24 red by
    value; the 2 green rows are the sphere line at b = 1.7, which never
    reaches the first panel of the outer region (declared).
    """
    case = _Case(body, point, direction, sigma_t=_SIGMA[group])
    x = case.coefficients(_q(body))
    B = case.rule.outflow()
    assert B.shape == (2, case.basis.size)
    f = R.region_polynomial(_q(body))
    for k, trav in _traversals(case, expected):
        _close(float(B[k] @ x), R.outflow(case.mp, trav, f, _SIGMA[group]), _T1_TOL, f"B_{k}")
    np.testing.assert_array_equal(B[len(expected):], 0.0)


_PER_FUNCTION = [r for r in _LINES if r[0] in ("sphere_solid_b0", "sphere_solid_b0.7", "sphere_hollow_b0.2_cavity",
                                               "slab_rising")]


@pytest.mark.l0
@pytest.mark.verifies("characteristic-traversal-integrals")
@pytest.mark.parametrize(("body", "point", "direction", "expected"), [r[1:] for r in _PER_FUNCTION],
                         ids=[r[0] for r in _PER_FUNCTION])
@pytest.mark.rests_on(_HERE + "test_the_outflow_of_a_per_region_polynomial_is_its_line_integral_attenuated_to_the_exit")
def test_each_basis_functions_outflow_is_its_line_integral(body, point, direction, expected) -> None:
    """[T2] Every column of B_k: the mpmath line integral of that Lagrange function (degree 1, one layer).

    The reference builds each panel's Lagrange functions on mpmath
    Gauss-Legendre roots over the basis's panel ends. A function the line does
    not reach integrates to exactly 0. 1e-13 of the largest column. First reds:
    a panel's integral scattered onto its neighbour's columns (the
    per-region-polynomial row T1 can be blind to an exchange of two equal
    panels' columns; this row is not); the nodes of one panel reversed.

    Rung 3 (the spec's RB4): ``sphere_solid_b0`` crosses the even panel at the
    centre, whose functions the reference writes in c^2 (``R.panel_function_even``,
    chosen here from the panel's lower end and the chart, never from
    ``PanelBasis.even``): the line-level witness of the even basis (the rung-2
    basis, Lagrange in c there, reds this row).
    """
    case = _Case(body, point, direction, sigma_t=_SIGMA[1], resolution=(1, 1, 0.5))
    B = case.rule.outflow()
    ends = case.basis.partition.breakpoints
    p1 = case.basis.degree + 1
    for k, trav in _traversals(case, expected):
        want = np.zeros(case.basis.size)
        for p, (a, b) in enumerate(zip(ends[:-1], ends[1:])):
            nodes = R.gl_nodes(a, b, p1)
            function = R.panel_function_even if (body.endswith("_solid") and a == 0.0) else R.panel_function
            for i in range(p1):
                want[p * p1 + i] = float(R.outflow(case.mp, trav, function(nodes, i, R.mpf(a), R.mpf(b)), _SIGMA[1]))
        scale = np.max(np.abs(want))
        assert np.max(np.abs(B[k] - want)) <= 1e-13 * scale, np.max(np.abs(B[k] - want)) / scale
        np.testing.assert_array_equal(B[k][want == 0.0], 0.0)


@pytest.mark.l0
@pytest.mark.verifies("characteristic-traversal-integrals")
@pytest.mark.parametrize(("body", "point", "direction", "expected"), _LINE_ARGS, ids=_LINE_IDS)
@pytest.mark.rests_on(_HERE + "test_the_panel_chord_unfolds_into_the_hand_counted_period")
def test_the_entry_response_is_the_integral_attenuated_from_the_entry(body, point, direction, expected) -> None:
    """[T3] A_k . q = int_k q(c(s)) e^{-tau(entry_k -> s)} ds, against mpmath; 1e-13 relative.

    Gated on its own, not through the reversal identity T4 (which the code
    builds in, so it is no evidence about A). First reds: A attenuated to the
    exit (A == B); the forward and backward sweeps exchanged. Region 0 even on a
    solid body (rung 3, RB4, as T1).
    """
    case = _Case(body, point, direction, sigma_t=_SIGMA[1])
    x = case.coefficients(_q(body))
    A = case.rule.entry_response()
    f = R.region_polynomial(_q(body))
    for k, trav in _traversals(case, expected):
        _close(float(A[k] @ x), R.entry_response(case.mp, trav, f, _SIGMA[1]), _T1_TOL, f"A_{k}")


@pytest.mark.foundation
@pytest.mark.parametrize(("body", "point", "direction", "expected"), [r[1:] for r in _LINES if r[1] == "slab"],
                         ids=[r[0] for r in _LINES if r[1] == "slab"])
@pytest.mark.rests_on(_HERE + "test_the_entry_response_is_the_integral_attenuated_from_the_entry")
def test_the_entry_response_of_a_traversal_is_the_outflow_of_its_reverse(body, point, direction, expected) -> None:
    """[T4] Reciprocity along one line: A of the forward traversal == B of the reversed one, bitwise.

    A DESIGN identity (``entry_response`` is built as the reversed sweep): this
    row is no evidence about A, which T3 gates against mpmath; it pins the
    identity so a later re-implementation of A on its own reds here first if it
    drifts by a single bit. On the slab both traversals of one transit are in
    the period.
    """
    case = _Case(body, point, direction, sigma_t=_SIGMA[0])
    A, B = case.rule.entry_response(), case.rule.outflow()
    np.testing.assert_array_equal(A[0], B[1])
    np.testing.assert_array_equal(A[1], B[0])


# ── T7-T9: the edges of the arc-length rule ──────────────────────────────


@pytest.mark.l0
@pytest.mark.verifies("characteristic-traversal-integrals")
@pytest.mark.parametrize(("mu", "sigma"), [(1.0, 1.0e3), (0.25, 250.0)], ids=["normal", "oblique"])
@pytest.mark.rests_on(_HERE + "test_the_outflow_of_a_per_region_polynomial_is_its_line_integral_attenuated_to_the_exit")
def test_a_thousand_mean_free_path_slot_integrates_to_the_closed_form(mu: float, sigma: float) -> None:
    """[T7, C2c] A one-region slab of optical width 1e3: B . 1 = A . 1 = (1 - e^{-Sigma l}) / Sigma; 4 ulp.

    Also q(x) = x against mpmath (the exit layer reads the source's value at
    the exit wall), 1e-13. First red: Gauss in arc length without the
    exponential grading (the whole exit layer, one mean free path, falls
    between two nodes of a panel 1/8 of the slot wide).
    """
    geometry = walls_of.slab((0.0, 1.0), (0,), left=walls_of.VACUUM, right=walls_of.VACUUM)
    point, direction = (0.5, 0.0, 0.0), (mu, float(np.sqrt((1.0 - mu) * (1.0 + mu))), 0.0)
    case = _Case("slab", point, direction, sigma_t=(sigma,), geometry=geometry)
    ones = np.ones(case.basis.size)
    length = mp.mpf(1) / mp.mpf(mu)
    want = float(-mp.expm1(-mp.mpf(sigma) * length) / mp.mpf(sigma))
    B, A = case.rule.outflow(), case.rule.entry_response()
    for k in range(2):
        _close(float(B[k] @ ones), want, 4 * _EPS, f"B_{k} . 1")
        _close(float(A[k] @ ones), want, 4 * _EPS, f"A_{k} . 1")
    q = [[0.0, 1.0]]
    x = case.coefficients(q)
    for k, rev in ((0, False), (1, True)):
        _close(float(B[k] @ x), R.outflow(case.mp, case.mp.traversal(0, rev), R.region_polynomial(q), (sigma,)),
               1e-13, f"B_{k} . x")


@pytest.mark.l0
@pytest.mark.verifies("characteristic-traversal-integrals")
@pytest.mark.filterwarnings("error")
@pytest.mark.rests_on(_HERE + "test_the_outflow_of_a_per_region_polynomial_is_its_line_integral_attenuated_to_the_exit")
def test_a_void_region_integrates_its_source_unattenuated() -> None:
    """[T8, B5c] Sigma_t = 0 in a region: q l there, no NaN, no warning (warnings are errors here).

    (a) a one-region void slab: B . 1 = A . 1 = l, 4 ulp; q = 0 gives exactly 0;
    (b) SLB3 with a void middle region and the MR3 sphere with a void region 1,
    against mpmath, 1e-13; every array finite. First red: the segment integral
    written q (1 - e^{-Sigma l}) / Sigma (NaN at Sigma = 0); the exponential
    grading's 1/Sigma formed on a void slot.
    """
    geometry = walls_of.slab((0.0, 1.3), (0,), left=walls_of.VACUUM, right=walls_of.VACUUM)
    case = _Case("slab", *_slab_dir(0.6), sigma_t=(0.0,), geometry=geometry)
    B, A = case.rule.outflow(), case.rule.entry_response()
    length = float(mp.mpf(1.3) / mp.mpf(0.6))
    for k in range(2):
        _close(float(B[k] @ np.ones(case.basis.size)), length, 4 * _EPS, "void B . 1")
        _close(float(A[k] @ np.ones(case.basis.size)), length, 4 * _EPS, "void A . 1")
    np.testing.assert_array_equal(B @ np.zeros(case.basis.size), 0.0)
    for body, point, direction, sigma, travs in (
        ("slab", *_slab_dir(0.6), (0.6, 0.0, 0.45), [(0, False), (0, True)]),
        ("sphere_solid", *_in_plane(0.3), (0.6, 0.0, 0.45), [(0, False)]),
    ):
        case = _Case(body, point, direction, sigma_t=sigma)
        x = case.coefficients(_Q)
        B = case.rule.outflow()
        assert np.all(np.isfinite(B)) and np.all(np.isfinite(case.rule.entry_response()))
        for k, (transit, rev) in enumerate(travs):
            _close(float(B[k] @ x), R.outflow(case.mp, case.mp.traversal(transit, rev), R.region_polynomial(_Q), sigma),
                   _T1_TOL, f"{body} void B_{k}")


#: A small cavity (rung 3, the spec's RB5): a line turning in the ordinary panel next to it, with b small against
#: that panel. On a solid body a line turning near the centre now turns in the EVEN panel, whose functions are
#: polynomials in c^2 = b^2 + |P Omega|^2 (s - s*)^2, polynomial in arc length: no branch point, so the turning
#: grading has no subject there and its witnesses moved here.
_CAVITY = (1e-5, 0.5, 1.5, 2.0)
_CAVITY_FIRST_END = 0.06250875           # the first panel end at the working resolution (3, 2, 0.5), by hand: r_0 + 0.25 / 4

_TURNING = [
    # (id, body, breakpoints, point, direction)
    ("sphere_cavity_b1e-4", "sphere_hollow", _CAVITY, *_in_plane(1e-4)),
    ("cylinder_cavity_b1e-4_wz0.8", "cylinder_hollow", _CAVITY, *_in_plane(1e-4, 0.8)),
    ("sphere_cavity_b_just_below_a_panel_end", "sphere_hollow", _CAVITY, *_in_plane(float(np.nextafter(_CAVITY_FIRST_END, 0.0)))),
    ("sphere_hollow_b_just_above_r0", "sphere_hollow", None, *_in_plane(0.4 + 1e-7)),
]


@pytest.mark.l0
@pytest.mark.verifies("characteristic-traversal-integrals")
@pytest.mark.parametrize(("body", "breakpoints", "point", "direction"), [r[1:] for r in _TURNING], ids=[r[0] for r in _TURNING])
@pytest.mark.rests_on(_HERE + "test_the_outflow_of_a_per_region_polynomial_is_its_line_integral_attenuated_to_the_exit")
def test_a_turning_slot_integrates_the_square_root_at_the_closest_approach(body, breakpoints, point, direction) -> None:
    """[T9, C10's analogue] A slot ending at the closest approach with b small against its panel, against mpmath.

    q = c (odd in c, so c(s) = sqrt(b^2 + s^2) is not polynomial in arc length
    and has branch points at s = +-i b, 1e-4 from the real segment, or
    1e-7 from a shell's tangency): pieces graded by halving toward the closest approach (down to the branch
    point's scale) make Gauss in arc length converge there. 1e-13 relative.
    First red: the turning grading removed (battery arm T3a).

    Re-posed 2026-10-06 (rung 3, RB5): the rows that turned on a solid body's
    centre panel moved to a hollow body with a cavity of radius 1e-5, where the
    line turns in an ordinary panel (on the solid body the even panel makes the
    integrand polynomial in arc length).
    """
    geometry = _geometry(body, breakpoints=breakpoints)
    case = _Case(body, point, direction, sigma_t=_SIGMA[0], geometry=geometry)
    assert int(case.period.rank) == 1, "the line misses the cavity"
    q = [[0.0, 1.0]] * 3
    x = case.coefficients(q)
    B = case.rule.outflow()
    _close(float(B[0] @ x), R.outflow(case.mp, case.mp.traversal(0, False), R.region_polynomial(q), _SIGMA[0]),
           _T1_TOL, "turning B_0")


_SMALL_RADIUS = [
    # (id, chart, breakpoints, inner wall?, b, wz): a slot whose near end is the crossing of a radius r_k much
    # smaller than the panel, so its branch points sit r_k/|P Omega| from that end (qa's probe_branch.py).
    ("sphere_cavity_1e-2_b_below_r0", "sphere", (0.01, 0.5, 2.0), True, 0.01 * (1.0 - 1e-6), 0.0),
    ("sphere_cavity_1e-3_b_half_r0", "sphere", (0.001, 0.5, 2.0), True, 0.0005, 0.0),
    ("sphere_inner_region_1e-2_b_below_r1", "sphere", (0.0, 0.01, 1.0), False, 0.01 * (1.0 - 1e-6), 0.0),
    ("cylinder_cavity_1e-2_b_below_r0", "cylinder", (0.01, 0.5, 2.0), True, 0.01 * (1.0 - 1e-6), 0.8),
]


#: T9c's rows as (resolution, row). The catching set is every row the exact ERR-099 defect reddens under
#: ``-O`` (`[M]` 2026-10-06, battery arm T17): the four rows at layers 0 (one panel per region, so no panel
#: end of the basis sits near the small radius) and the 1e-3 cavity at layers 2. The controls are the three
#: layers-2 rows it leaves green: the basis's own grading puts a panel end at the small radius and hides the
#: slot's branch point there (the archivist's re-drop, 1 of 4 rows red at layers 2).
_LAYERS0, _LAYERS2 = (3, 0, 0.5), (3, 2, 0.5)
_ERR099_CATCHERS = [pytest.param(_LAYERS0, *r[1:], id=f"{r[0]}-layers0") for r in _SMALL_RADIUS] + [
    pytest.param(_LAYERS2, *r[1:], id=f"{r[0]}-layers2") for r in _SMALL_RADIUS if r[0] == "sphere_cavity_1e-3_b_half_r0"]
_ERR099_CONTROLS = [pytest.param(_LAYERS2, *r[1:], id=f"{r[0]}-layers2") for r in _SMALL_RADIUS
                    if r[0] != "sphere_cavity_1e-3_b_half_r0"]


def _small_radius_outflow(chart, breakpoints, hollow, b, wz, resolution) -> None:
    n = len(breakpoints) - 1
    laws = {"outer": walls_of.VACUUM, **({"inner": walls_of.VACUUM} if hollow else {})}
    geometry = getattr(walls_of, chart)(breakpoints, tuple(range(n)), **laws)
    sigma = _SIGMA[0][:n]
    case = _Case(f"{chart}_x", *_in_plane(b, wz), sigma_t=sigma, geometry=geometry, resolution=resolution)
    # rung 3 (RB6): a solid body's region 0 is the even panel at the centre, so its source is even (c^2); the
    # branch point ERR-099 is about sits at the start of region 1's slot, where q = c stays
    q = [[0.0, 1.0]] * n if hollow else [[0.0, 0.0, 1.0]] + [[0.0, 1.0]] * (n - 1)
    x = case.coefficients(q)
    B = case.rule.outflow()
    transits = case.mp.transits()
    assert int(case.period.rank) == len(transits)
    for k in range(len(transits)):
        _close(float(B[k] @ x), R.outflow(case.mp, case.mp.traversal(k, False), R.region_polynomial(q), sigma),
               _T1_TOL, f"B_{k}")


@pytest.mark.foundation
@pytest.mark.catches("ERR-099")
@pytest.mark.parametrize(("resolution", "chart", "breakpoints", "hollow", "b", "wz"), _ERR099_CATCHERS)
@pytest.mark.rests_on(_HERE + "test_a_turning_slot_integrates_the_square_root_at_the_closest_approach")
def test_a_slot_starting_at_a_small_radius_integrates_its_branch_point(resolution, chart, breakpoints, hollow, b,
                                                                       wz) -> None:
    """[T9c] B_k . q for q = c on every traversal of a line meeting a radius r_k much smaller than its panel; mpmath.

    The slot beyond the crossing of r_k does not end at the closest approach,
    yet c(s) = sqrt(b^2 + |P Omega|^2 s^2) has its branch points r_k/|P Omega|
    from that end: a small cavity (r_0 = 0.01, 0.001) crossed by a line through
    it, and a solid body whose inner region has radius 0.01. 1e-13 relative,
    the T1 band. Every row here is reddened by ERR-099 (grading only toward
    the closest approach, battery arm T17) under ``-O``; the rows it leaves
    green are ``test_a_small_radius_hidden_by_a_panel_end_is_the_err099_control``.
    First reds: arm T17 (ERR-099); the grading toward the near end removed (arm T3a).
    """
    _small_radius_outflow(chart, breakpoints, hollow, b, wz, resolution)


@pytest.mark.foundation
@pytest.mark.parametrize(("resolution", "chart", "breakpoints", "hollow", "b", "wz"), _ERR099_CONTROLS)
@pytest.mark.rests_on(_HERE + "test_a_slot_starting_at_a_small_radius_integrates_its_branch_point")
def test_a_small_radius_hidden_by_a_panel_end_is_the_err099_control(resolution, chart, breakpoints, hollow, b,
                                                                     wz) -> None:
    """[T9c control] The same integrals with the basis graded (layers 2): ERR-099 leaves these green (declared).

    A panel end of the basis lands near the small radius, so the piece next to
    the branch point is already short; the row still gates the value (1e-13)
    and reds under the grading removed entirely (arm T3a), but it is not a
    catcher of ERR-099.
    """
    _small_radius_outflow(chart, breakpoints, hollow, b, wz, resolution)


_HALVING = [
    # (id, body, point, direction): turning in the ordinary panel next to a cavity of radius 1e-5 (rung 3, RB5)
    ("sphere_cavity_b1e-2", "sphere_hollow", *_in_plane(1e-2)),
    ("sphere_cavity_b3e-3", "sphere_hollow", *_in_plane(3e-3)),
    ("sphere_cavity_b1e-4_control", "sphere_hollow", *_in_plane(1e-4)),
    ("cylinder_cavity_b1e-2_wz0.8", "cylinder_hollow", *_in_plane(1e-2, 0.8)),
    ("cylinder_cavity_b3e-3_wz0.8", "cylinder_hollow", *_in_plane(3e-3, 0.8)),
    ("cylinder_cavity_b1e-4_wz0.8_control", "cylinder_hollow", *_in_plane(1e-4, 0.8)),
]


@pytest.mark.foundation
@pytest.mark.parametrize(("body", "point", "direction"), [r[1:] for r in _HALVING], ids=[r[0] for r in _HALVING])
@pytest.mark.rests_on(_HERE + "test_a_turning_slot_integrates_the_square_root_at_the_closest_approach")
def test_the_turning_grading_holds_a_turning_line_next_to_a_small_cavity(body, point, direction) -> None:
    """[T9b] On a line turning next to a cavity of radius 1e-5, B_0 . q and psi(t) . q stay within 1e-13 of mpmath
    (q = c), at the node counts the law derives per piece.

    Until P3 step 3f this row was the witness of the halving RATIO at 7 points per piece; since the node counts are
    derived (route B, ``grading.gauss_nodes``) there is no point count to coarsen, and the walk's premise, each piece
    no wider than its distance to the branch points over c, is gated geometrically and n-free (PP1,
    ``test_every_piece_is_no_wider_than_its_branch_distance_over_the_tightness``). This row is the value leg. First
    red: the hp walk disabled (`[M]` 2026-10-13 on the route-B prototype, arm ``b_hp_none``: 4 of 6 rows, the two
    ``_control`` rows green as their names declare; re-run on production, ``carve/arms_*.log``). Same 1e-13 band as T1.
    """
    geometry = _geometry(body, breakpoints=_CAVITY)
    case = _Case(body, point, direction, sigma_t=_SIGMA[0], geometry=geometry)
    assert int(case.period.rank) == 1, "the line misses the cavity"
    q = [[0.0, 1.0]] * 3
    x = case.coefficients(q)
    f = R.region_polynomial(q)
    trav = case.mp.traversal(0, False)
    _close(float(case.rule.outflow()[0] @ x), R.outflow(case.mp, trav, f, _SIGMA[0]), _T1_TOL, "B_0")
    _, ts = _points_on(case, 0)
    for i, tm in enumerate(ts):
        rule, t = case.read_at(tm)
        psi = rule.angular_flux(np.array([t]), np.zeros((2, case.basis.size)))[0] @ x
        _close(float(psi), R.partial_outflow(case.mp, trav, f, _SIGMA[0], R.mpf(t)), _T1_TOL, f"psi(t_{i})")


# ── T5, T6: the Volterra triangle ────────────────────────────────────────

#: Spec §4's group-1 Sigma_t times 5: every slot of SLB3 at mu = 0.6 is above 2 mfp, so each is cut into
#: exponential pieces and the triangle's inner rule must start at the PIECE's start (the other rows are one
#: piece per slot, where piece and slot starts coincide: declared, `[M]` battery arm T7).
_SIGMA_THICK = tuple(5.0 * v for v in _SIGMA[1])

_VOLTERRA = [
    # (id, body, point, direction, the transits read forward, Sigma_t, mpmath GL degrees)
    ("sphere_solid_b0.7", "sphere_solid", *_in_plane(0.7), [0], _SIGMA[0], (4, 5)),
    ("sphere_solid_b0.3", "sphere_solid", *_in_plane(0.3), [0], _SIGMA[0], (4, 5)),
    ("sphere_hollow_b0.2_cavity", "sphere_hollow", *_in_plane(0.2), [0, 1], _SIGMA[0], (4, 5)),
    ("slab_rising", "slab", *_slab_dir(0.6), [0], _SIGMA[0], (4, 5)),
    ("slab_thick_pieces", "slab", *_slab_dir(0.6), [0], _SIGMA_THICK, (5, 6)),
    ("cylinder_solid_b0.7_wz0.8", "cylinder_solid", *_in_plane(0.7, 0.8), [0], _SIGMA[0], (4, 5)),
]


@pytest.mark.foundation
@pytest.mark.parametrize(("body", "point", "direction", "transits", "sigma", "degrees"), [r[1:] for r in _VOLTERRA],
                         ids=[r[0] for r in _VOLTERRA])
@pytest.mark.rests_on(_HERE + "test_the_outflow_of_a_per_region_polynomial_is_its_line_integral_attenuated_to_the_exit")
def test_the_volterra_triangle_is_the_double_integral_along_the_line(body, point, direction, transits, sigma,
                                                                    degrees) -> None:
    """[T5] f . V g = sum over the line's transits read forward of int f(s) int_{s' < s} e^{-tau(s', s)} g(s') ds' ds.

    f and g distinct per-region cubics (in the exactness family). The
    reference: mpmath's Gauss-Legendre outer rule at two degrees (24 and 48
    points per smooth piece, required to agree to 1e-18 first: the reference's
    own convergence control, neither count the code's) with the inner integral
    carried cumulatively by tanh-sinh. 1e-12 relative. First reds: the inner
    rule run from the slot's start instead of the piece's (the upstream
    attenuation counted twice); the carried term attenuated by the whole piece
    instead of to each node (red on ``slab_thick_pieces`` only, the one row with
    several pieces per slot); the triangle's (i, j) transposed (red on the slab
    only: on a radial chord the inbound and outbound slots' triangles are each
    other's transpose, so transposing both is in the row's stabiliser,
    declared); the carried rows dropped; the reversed traversals summed in.
    """
    case = _Case(body, point, direction, sigma_t=sigma)
    f, g = R.region_polynomial(_Q), R.region_polynomial([[0.4, -0.1, 0.2, 0.0], [2.0, 0.3, -0.5, 0.1], [0.7, 0.2, 0.3, -0.05]])
    xf = R.sample(f, np.asarray(case.basis.region), np.asarray(case.basis.nodes))
    xg = R.sample(g, np.asarray(case.basis.region), np.asarray(case.basis.nodes))
    V = case.rule.volterra(np.array(1.0))
    assert V.shape == (case.basis.size, case.basis.size)
    want = []
    for degree in degrees:
        want.append(sum(R.volterra_bilinear(case.mp, case.mp.traversal(t, False), f, g, sigma, degree=degree)
                        for t in transits))
    assert abs(want[0] - want[1]) <= mp.mpf("1e-18") * abs(want[1]), "the mpmath reference has not converged"
    _close(float(xf @ V @ xg), want[1], 1e-12, "f . V g")


@pytest.mark.foundation
@pytest.mark.parametrize("mu", [0.6, 0.15])
@pytest.mark.rests_on(_HERE + "test_the_volterra_triangle_is_the_double_integral_along_the_line")
def test_the_reversed_lines_triangle_is_the_transpose(mu: float) -> None:
    """[T5b] V of a slab line read along -Omega is V^T of the line along +Omega (reciprocity along one line).

    Each line's triangle runs over its own transit forward only, so the two
    are two different traversals of one transit; the slab's regions differ, so
    neither triangle is symmetric (asserted: the row would be vacuous on a
    palindromic chord). Tolerance 1e-13 of the largest entry (two rules with
    their nodes in mirrored order). First reds: the reversed traversal summed
    into each line (both become V + V^T: the asymmetry leg reds); the inner
    rule graded toward the piece's start.
    """
    point, direction = _slab_dir(mu)
    plus = _Case("slab", point, direction, sigma_t=_SIGMA[1]).rule.volterra(np.array(1.0))
    minus = _Case("slab", point, tuple(-v for v in direction), sigma_t=_SIGMA[1]).rule.volterra(np.array(1.0))
    scale = np.max(np.abs(plus))
    assert np.max(np.abs(plus - plus.T)) > 1e-3 * scale, "the line's triangle is symmetric: the row is vacuous"
    assert np.max(np.abs(minus - plus.T)) <= 1e-13 * scale, np.max(np.abs(minus - plus.T)) / scale


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_volterra_triangle_is_the_double_integral_along_the_line")
def test_the_triangle_of_a_batch_is_the_weighted_sum_of_its_lines() -> None:
    """[T6] volterra(w) over a batch == sum_L w_L volterra over each line alone; a parallel line and a zero weight add 0.

    The batch mixes ranks 0 (an axial cylinder line), 1 and 2 (through the
    cavity is not available on a solid cylinder; an oblique line and a line
    through the axis). 8 ulp of the largest entry times the batch size.
    First reds: the weights applied by position after the parallel line is
    dropped (a shift); the weight broadcast over the traversal axis.
    """
    posed = _geometry("cylinder_solid")
    points = np.array([(0.7, 0.0, 0.0), (0.7, 0.0, 0.0), (0.0, 0.0, 0.0), (1.2, 0.0, 0.0)])
    directions = np.array([(0.0, 0.0, 1.0), (0.0, 0.6, 0.8), (0.0, 1.0, 0.0), (0.0, 0.8, -0.6)])
    weights = np.array([2.5, 0.3, 1.7, 0.0])
    basis = PanelBasis.of(ConcentricPartition.of(posed.geometry), *_RES)
    walls = Walls.of(*posed)

    def rule(p, d):
        return TraversalRule.of(Line.through(p, d), basis, walls, np.asarray(_SIGMA[0]))

    batch = rule(points, directions).volterra(weights)
    alone = [rule(p, d).volterra(np.array(1.0)) for p, d in zip(points, directions)]
    np.testing.assert_array_equal(alone[0], 0.0)
    want = sum(w * v for w, v in zip(weights, alone))
    scale = max(np.max(np.abs(v)) for v in alone)
    assert np.all(np.isfinite(batch))
    assert np.max(np.abs(batch - want)) <= 8 * len(weights) * _EPS * scale


# ── T5c: a thick slot read by the triangle and by psi (qa finding 3) ────


def _derivatives(coefficients) -> list[list]:
    """The polynomial (ascending mpmath coefficients) and all its derivatives."""
    out = [[R.mpf(a) for a in coefficients]]
    while len(out[-1]) > 1:
        c = out[-1]
        out.append([c[k] * k for k in range(1, len(c))])
    return out


def _poly(c, x):
    return mp.polyval(list(reversed(c)), x)


def _attenuated_to(poly, sigma, t):
    r"""int_0^t P(x) e^{-Sigma (t - x)} dx, closed form: sum_k (-1)^k [P^(k)(t) - e^{-Sigma t} P^(k)(0)] / Sigma^(k+1)."""
    return sum((-1) ** k * (_poly(d, t) - mp.exp(-sigma * t) * _poly(d, 0)) / sigma ** (k + 1)
               for k, d in enumerate(_derivatives(poly)))


def _volterra_closed(f, g, sigma, length):
    r"""int_0^L f(s) int_0^s e^{-Sigma (s - s')} g(s') ds' ds, closed form.

    The inner integral is A(s) - e^{-Sigma s} C with A the polynomial
    sum_k (-1)^k g^(k)(s) / Sigma^(k+1) and C = sum_k (-1)^k g^(k)(0) / Sigma^(k+1);
    int_0^L f A is a polynomial integral, and int_0^L f e^{-Sigma s} =
    sum_k [f^(k)(0) - e^{-Sigma L} f^(k)(L)] / Sigma^(k+1).
    """
    gd, fd = _derivatives(g), _derivatives(f)
    n = max(len(f) + len(g), 2)
    a = [mp.mpf(0)] * n
    for k, d in enumerate(gd):
        for i, c in enumerate(d):
            a[i] += (-1) ** k * c / sigma ** (k + 1)
    fa = [mp.mpf(0)] * (len(f) + n)
    for i, x in enumerate(fd[0]):
        for j, y in enumerate(a):
            fa[i + j] += x * y
    poly_part = sum(c * length ** (i + 1) / (i + 1) for i, c in enumerate(fa))
    C = sum((-1) ** k * _poly(d, 0) / sigma ** (k + 1) for k, d in enumerate(gd))
    exp_part = sum((_poly(d, 0) - mp.exp(-sigma * length) * _poly(d, length)) / sigma ** (k + 1)
                   for k, d in enumerate(fd))
    return poly_part - C * exp_part


@pytest.mark.foundation
@pytest.mark.catches("ERR-100")
@pytest.mark.parametrize("sigma", [200.0, 1000.0], ids=["200mfp", "1000mfp"])
@pytest.mark.rests_on(_HERE + "test_a_thousand_mean_free_path_slot_integrates_to_the_closed_form",
                      _HERE + "test_the_volterra_triangle_is_the_double_integral_along_the_line")
def test_a_thick_slot_is_read_by_the_triangle_and_by_psi_in_closed_form(sigma: float) -> None:
    """[T5c, qa finding 3] One panel (p = 2) over a slab of width 1 at mu = 1 and Sigma = 200 or 1000 mfp.

    psi(t) . q at the middle, just past the middle piece (64 mfp before the
    exit, then a hair beyond) and half a mean free path before the exit; and
    f . V g. Every reference is a CLOSED FORM in mpmath (a polynomial times an
    exponential integrated by parts), no quadrature on the reference side.
    1e-13 relative. First red: the piece integrals taken on the piece's own
    Gauss nodes, ungraded (battery arm T14; qa measured 5e-6 at 200 mfp and
    0.54 at 1000 mfp in the code before round 3).
    """
    geometry = walls_of.slab((0.0, 1.0), (0,), left=walls_of.VACUUM, right=walls_of.VACUUM)
    case = _Case("slab", (0.5, 0.0, 0.0), (1.0, 0.0, 0.0), sigma_t=(sigma,), geometry=geometry, resolution=(2, 0, 0.5))
    assert case.basis.n_panels == 1
    f, g = [0.3, 1.0, -0.7], [1.1, -0.4, 0.9]
    xf, xg = case.coefficients([f]), case.coefficients([g])
    S = R.mpf(sigma)
    mfp = 1.0 / sigma
    for i, target in enumerate((0.5, 1.0 - 64.0 * mfp, 1.0 - 64.0 * mfp + 0.25 * mfp, 1.0 - 0.5 * mfp)):
        rule, ti = case.read_at(target)                       # the slab refined at the point's x (P3 step 3)
        psi = rule.angular_flux(np.array([ti]), np.zeros((2, case.basis.size)))[0] @ xg
        _close(float(psi), _attenuated_to(g, S, R.mpf(ti)), _T1_TOL, f"psi at t = {ti}")
    V = case.rule.volterra(np.array(1.0))
    _close(float(xf @ V @ xg), _volterra_closed(f, g, S, mp.mpf(1)), _T1_TOL, "f . V g")


# ── T10-T12: the angular flux along the line ─────────────────────────────

#: Interior fractions of a transit; the middle (0.5) is dropped since P3 step 3: on a radial line it is the closest
#: approach, tangent to its own level, which no chord crosses.
_FRACTIONS = (0.013, 0.37, 0.63, 0.91, 1.0 - 1e-9)


def _points_on(case: _Case, transit: int) -> tuple[np.ndarray, list]:
    t_in, t_out, _ = case.mp.transits()[transit]
    ts = [t_in + (t_out - t_in) * mp.mpf(fr) for fr in _FRACTIONS]
    return np.array([float(t) for t in ts]), ts


@pytest.mark.foundation
@pytest.mark.parametrize(("body", "point", "direction", "expected"),
                         [r[1:] for r in _LINES if r[0] in ("sphere_solid_b0.7", "sphere_hollow_b0.2_cavity",
                                                             "slab_falling", "cylinder_solid_b0.3_wz-0.6")],
                         ids=[r[0] for r in _LINES if r[0] in ("sphere_solid_b0.7", "sphere_hollow_b0.2_cavity",
                                                                "slab_falling", "cylinder_solid_b0.3_wz-0.6")])
@pytest.mark.rests_on(_HERE + "test_the_outflow_of_a_per_region_polynomial_is_its_line_integral_attenuated_to_the_exit")
def test_the_vacuum_angular_flux_is_the_source_integral_since_the_entry(body, point, direction, expected) -> None:
    """[T10] With no inflow, psi(t) . q = int_{entry}^{t} q e^{-tau(s -> t)} ds along the transit holding t; mpmath.

    Points at five fractions of each transit, one a hair before its exit (where
    psi meets B_k). 1e-13 relative. The refusal: a point in the cavity or past
    the outer wall is refused by its own fragment. First reds: the carried
    integral attenuated by the whole piece instead of to the point; the
    transit's earlier pieces dropped; the inflow column of the wrong traversal
    (T11).
    """
    case = _Case(body, point, direction, sigma_t=_SIGMA[1])
    x = case.coefficients(_Q)
    zero = np.zeros((2, case.basis.size))
    f = R.region_polynomial(_Q)
    for k, (transit, rev, _end) in enumerate(expected):
        if rev:
            continue
        _, ts = _points_on(case, transit)
        trav = case.mp.traversal(transit, False)
        for i, tm in enumerate(ts):
            rule, t = case.read_at(tm)
            psi = rule.angular_flux(np.array([t]), zero)
            assert psi.shape == (1, case.basis.size)
            _close(float(psi[0] @ x), R.partial_outflow(case.mp, trav, f, _SIGMA[1], R.mpf(t)), _T1_TOL, f"psi(t_{i})")
    outside = float(case.mp.transits()[-1][1]) + 0.1
    with pytest.raises(ValueError, match="is a crossing on a transit of its line"):
        case.rule.angular_flux(np.array([outside]), zero)
    # an interior parameter of the unrefined chord is refused by the same fragment (P3 step 3: read at a crossing)
    with pytest.raises(ValueError, match="is a crossing on a transit of its line"):
        case.rule.angular_flux(np.array([float(_points_on(case, 0)[1][1])]), zero)
    # no point is a read of nothing, not an ambiguous reshape (qa, #586, 2026-10-07)
    assert case.rule.angular_flux(np.zeros(0), zero).shape == (0, case.basis.size)


_EXIT_READS = [r for r in _LINES if r[0] in ("sphere_solid_b0.7", "sphere_hollow_b0.2_cavity", "cylinder_solid_b0.7_wz0.8",
                                              "cylinder_solid_b0.3_wz-0.6", "slab_rising", "slab_falling")]


@pytest.mark.foundation
@pytest.mark.parametrize(("body", "point", "direction", "expected"), [r[1:] for r in _EXIT_READS],
                         ids=[r[0] for r in _EXIT_READS])
@pytest.mark.rests_on(_HERE + "test_the_vacuum_angular_flux_is_the_source_integral_since_the_entry")
def test_the_angular_flux_is_read_at_a_transits_exit_crossing(body, point, direction, expected) -> None:
    """[T10b] psi at the parameter of each transit's CLOSING crossing (the kernel's) is accepted and equals B there.

    The point is spelled as the kernel spells the wall, ``chord.crossings.parameter``
    at the transit's ``stop_slot``, not as start + length (which rounds an ulp
    apart and is legitimately refused). With no inflow, psi there is the
    outflow of the forward traversal: against mpmath, 1e-13. First red: the
    location test ``t - start <= length`` (qa: 64 of 119 exit reads refused;
    battery arm T16).
    """
    case = _Case(body, point, direction, sigma_t=_SIGMA[1])
    x = case.coefficients(_Q)
    chord = case.period.chord
    transits = chord.transits
    f = R.region_polynomial(_Q)
    zero = np.zeros((2, case.basis.size))
    for transit in range(2):
        if not bool(transits.present[transit]):
            continue
        t_exit = float(chord.crossings.parameter[int(transits.stop_slot[transit])])
        psi = case.rule.angular_flux(np.array([t_exit]), zero)[0] @ x
        _close(float(psi), R.outflow(case.mp, case.mp.traversal(transit, False), f, _SIGMA[1]), _T1_TOL,
               f"psi at the exit of transit {transit}")


def _unfolded_psi(case: _Case, period_by_hand, amplitudes, f, sigma, k0: int, t: R.Mpf) -> R.Mpf:
    """psi(t) on traversal k0 as the explicit backward march: the wall-by-wall sum, no geometric-series division."""
    travs = [case.mp.traversal(tr, rev) for tr, rev, _ in period_by_hand]
    m = len(travs)
    B = [R.outflow(case.mp, tr, f, sigma) for tr in travs]
    gain = [mp.mpf(amplitudes[i]) * mp.exp(-R.optical_depth(travs[i], sigma, travs[i].t_in, travs[i].t_out))
            for i in range(m)]
    inflow, weight, j = mp.mpf(0), mp.mpf(1), 1
    while True:
        i = (k0 - j) % m
        inflow += weight * mp.mpf(amplitudes[i]) * B[i]
        weight *= gain[i]
        j += 1
        if weight < mp.mpf("1e-28"):
            break
    trav = travs[k0]
    return R.partial_outflow(case.mp, trav, f, sigma, t) + mp.exp(-R.optical_depth(trav, sigma, trav.entry, t)) * inflow


_CLOSED = [
    # (id, body, point, direction, period by hand, the amplitude of each traversal's EXIT wall)
    ("slab_mirrors_rising", "slab", *_slab_dir(0.6), [(0, False, "last"), (0, True, "first")], (_A_OUT, _A_IN)),
    ("slab_mirrors_falling", "slab", (0.7, 0.0, 0.0), (-0.3, 0.0, float(np.sqrt(0.91))),
     [(0, False, "first"), (0, True, "last")], (_A_IN, _A_OUT)),
    ("sphere_hollow_cavity", "sphere_hollow", *_in_plane(0.2), [(0, False, "first"), (1, False, "last")], (_A_IN, _A_OUT)),
    ("sphere_solid_partial", "sphere_solid", *_in_plane(0.7), [(0, False, "last")], (_A_OUT,)),
    # B6, ERR-035's regime: the symmetric slab at an intermediate albedo, under the non-symmetric _Q (spec §5 B6).
    ("slab_symmetric_half", "slab_half", *_slab_dir(0.6), [(0, False, "last"), (0, True, "first")], (0.5, 0.5)),
    # The cylinder's rank-1 partial mirror on an oblique line (``peierls-greens-cylinder-T``), step (e1b).
    ("cylinder_solid_partial", "cylinder_solid", *_in_plane(0.7, 0.8), [(0, False, "last")], (_A_OUT,)),
    # The annulus's rank-2 ray through the cavity (``peierls-greens-annulus-through-rank2``), step (e1b).
    ("cylinder_hollow_cavity", "cylinder_hollow", *_in_plane(0.2, 0.8), [(0, False, "first"), (1, False, "last")],
     (_A_IN, _A_OUT)),
]


@pytest.mark.l0
@pytest.mark.verifies("characteristic-closure", "peierls-greens-cylinder-T",
                      "peierls-greens-annulus-through-rank2", "peierls-greens-hollow-sph-through-rank2")
@pytest.mark.catches("ERR-035")
@pytest.mark.parametrize(("body", "point", "direction", "period_by_hand", "amplitudes"), [r[1:] for r in _CLOSED],
                         ids=[r[0] for r in _CLOSED])
@pytest.mark.rests_on(_HERE + "test_the_vacuum_angular_flux_is_the_source_integral_since_the_entry",
                      _CLOSURE + "test_the_inflow_is_the_unfolded_wall_by_wall_sum")
def test_the_closed_angular_flux_is_the_unfolded_backward_path(body, point, direction, period_by_hand, amplitudes) -> None:
    """[T11, B2's line half] psi(t) with the period's inflow == the backward path marched wall by wall in mpmath.

    The inflow is ``LinePeriod.inflow(rule.optical_depth, rule.outflow())``;
    the reference marches back from t through each wall, multiplying by the
    amplitude of the wall the backward path reflects at (written by hand per
    fixture: distinct inner and outer amplitudes 0.3 and 0.6, so a pairing
    swap moves every reading), summing until the weight is below 1e-28.
    1e-13 relative. First reds: the inflow attenuated from the exit instead of
    the entry (the prototype's ``closure_wrong_end`` control); the inflow
    column of the reversed traversal read for the forward one (slab);
    B_k read from the entry response.

    ``slab_symmetric_half`` is the spec's B6, ERR-035's regime: the symmetric
    slab at albedo 0.5 on both walls under the non-symmetric ``_Q``. ERR-035's
    heuristic (the out-and-back integral with no amplitude at the inner
    reflection, over the same denominator) agrees with the honest cycle only at
    albedo 0 and 1, so this row and every rank-2 row red when it is re-dropped
    (``[M]`` 2026-10-10, arm ``err035`` of
    ``scratch/characteristic_architecture/p1_step_e/ta_e1b/battery``). It
    succeeds the old family's
    ``test_rank1_path_now_agrees_with_rank2_via_delegation_after_ERR035_fix``.

    Step (e1b) also added the cylinder's rank-1 partial mirror on an oblique
    line and the annulus's rank-2 ray through the cavity, so the row carries
    the old page's closure labels the deleted rows verified: the cylinder's
    rank-1 T (``peierls-greens-cylinder-T``) and the rank-2 resolvent of the
    hollow sphere and the annulus (``peierls-greens-hollow-sph-through-rank2``,
    ``peierls-greens-annulus-through-rank2``, whose 2 x 2 T is the cycle's
    least solution written in closed form). First reds of the two new rows
    (``[M]`` 2026-10-10): the cylinder's obliquity dropped (arm
    ``cylinder-no-obliquity``, both cylinder rows); ERR-035's arm (the annulus row).
    """
    case = _Case(body, point, direction, sigma_t=_SIGMA[0])
    x = case.coefficients(_Q)
    f = R.region_polynomial(_Q)
    for k0, (transit, rev, _end) in enumerate(period_by_hand):
        if rev:
            continue
        _, ts = _points_on(case, transit)
        for i, tm in enumerate(ts):
            rule, t = case.read_at(tm)
            inflow = rule.period.inflow(rule.optical_depth, rule.outflow())
            psi = rule.angular_flux(np.array([t]), inflow)[0] @ x
            _close(float(psi), _unfolded_psi(case, period_by_hand, amplitudes, f, _SIGMA[0], k0, R.mpf(t)), _T1_TOL,
                   f"psi(t_{i}) on traversal {k0}")


#: ``peierls-greens-slab-T``'s hypotheses: a homogeneous slab [0, L] (Sigma_t 0.7, L = 2.3) with ONE albedo on both
#: faces and a source symmetric about its mid-plane.
_SLAB_T_L, _SLAB_T_SIGMA, _SLAB_T_ALBEDOS = 2.3, 0.7, (0.3, 0.5, 0.9)
#: q(x) = 1 + 0.8 (x - L/2)^2, as monomial coefficients in x (degree 2: inside the cubic panels' exactness family).
_SLAB_T_Q = [[1.0 + 0.8 * (_SLAB_T_L / 2) ** 2, -0.8 * _SLAB_T_L, 0.8]]


def _slab_t_closed_form(albedo: float, mu: float, q) -> float:
    r"""``peierls-greens-slab-T``: psi_surf = alpha B / (1 - alpha e^{-tau}), tau = Sigma_t L / mu, B the one-transit
    source integral attenuated to the exit, by mpmath quadrature along the transit written here (no reference
    object, no cycle): x(s) = mu s, B = int_0^{L/mu} q(x(s)) e^{-Sigma_t (L/mu - s)} ds."""
    with mp.workdps(40):
        length = mp.mpf(_SLAB_T_L) / mp.mpf(mu)
        poly = lambda x: mp.polyval([mp.mpf(c) for c in reversed(q[0])], x)                 # noqa: E731
        b = mp.quad(lambda s: poly(mp.mpf(mu) * s) * mp.exp(-mp.mpf(_SLAB_T_SIGMA) * (length - s)), [0, length])
        tau = mp.mpf(_SLAB_T_SIGMA) * length
        return float(mp.mpf(albedo) * b / (1 - mp.mpf(albedo) * mp.exp(-tau)))


@pytest.mark.l0
@pytest.mark.verifies("peierls-greens-slab-T")
@pytest.mark.parametrize("albedo", _SLAB_T_ALBEDOS)
@pytest.mark.rests_on(_HERE + "test_the_outflow_of_a_per_region_polynomial_is_its_line_integral_attenuated_to_the_exit",
                      _CLOSURE + "test_the_inflow_is_the_unfolded_wall_by_wall_sum")
def test_a_symmetric_slabs_surface_inflow_is_the_one_transit_closed_form(albedo: float) -> None:
    r"""[``peierls-greens-slab-T``] On a homogeneous slab with one albedo alpha on both faces and a source symmetric
    about its mid-plane, the inflow at each traversal's entry is alpha B / (1 - alpha e^{-tau}), B the one-transit
    outflow, to 1e-13 relative, at alpha in {0.3, 0.5, 0.9}.

    The closed form holds only under these hypotheses (qa, step (f) F1: the
    difference from the honest cycle is E alpha^2 (B_LR - B_RL) / ((E alpha - 1)
    (E alpha + 1)), zero exactly when the two directed transits' outflows are
    equal), so the row poses them; the honest cycle for any source is
    ``test_the_closed_angular_flux_is_the_unfolded_backward_path``. The right side
    is written here (mpmath along the transit), never the reference's cycle
    (X4). First reds (``[M]`` 2026-10-10,
    ``scratch/characteristic_architecture/p1_step_f/ta/``): (a) the source made
    non-symmetric (``_SLAB_T_Q`` with a linear term, arm ``slab-t-asymmetric``):
    red, the hypothesis is load-bearing; (b) ERR-035's heuristic
    alpha B_period / (1 - alpha^2 e^{-2 tau}) put in place of the closed form
    (arm ``slab-t-heuristic``): red at every intermediate alpha; and the
    production ERR-035 arm (``err035``) reds it too.
    """
    mu = 0.6
    geometry = walls_of.slab((0.0, _SLAB_T_L), (0,), left=_spec(albedo), right=_spec(albedo))
    case = _Case("slab", *_slab_dir(mu), sigma_t=[_SLAB_T_SIGMA], geometry=geometry)
    rule = case.rule
    inflow = case.period.inflow(rule.optical_depth, rule.outflow()) @ case.coefficients(_SLAB_T_Q)
    want = _slab_t_closed_form(albedo, mu, _SLAB_T_Q)
    for k in range(2):
        assert abs(float(inflow[k]) / want - 1.0) < _T1_TOL, (k, float(inflow[k]), want)


@pytest.mark.foundation
@pytest.mark.parametrize("body", ["slab", "sphere_hollow", "sphere_solid"])
@pytest.mark.rests_on(_HERE + "test_the_vacuum_angular_flux_is_the_source_integral_since_the_entry",
                      _CLOSURE + "test_every_amplitude_zero_adds_nothing")
def test_vacuum_walls_add_nothing_on_a_line(body: str) -> None:
    """[T12, B4a at line level] Every amplitude 0: the inflow is exactly 0 and psi equals the no-inflow psi, bitwise.

    The edge where the honest answer equals a foundation exactly (lessons §4).
    First red: a default inflow substituted for a zero amplitude.
    """
    geometry = _geometry(body, a_in=0.0, a_out=0.0)
    point, direction = _slab_dir(0.6) if body == "slab" else _in_plane(0.2)
    case = _Case(body, point, direction, sigma_t=_SIGMA[0], geometry=geometry)
    rule = case.rule
    inflow = case.period.inflow(rule.optical_depth, rule.outflow())
    np.testing.assert_array_equal(inflow, 0.0)
    _, ts = _points_on(case, 0)
    for tm in ts:
        refined, t = case.read_at(tm)
        on_level = refined.period.inflow(refined.optical_depth, refined.outflow())
        np.testing.assert_array_equal(on_level, 0.0)
        np.testing.assert_array_equal(refined.angular_flux(np.array([t]), on_level),
                                      refined.angular_flux(np.array([t]), np.zeros_like(on_level)))


# ── P3 step 3e: the speed family, the refinement and the middle cut (qa's and the elegance review's witnesses) ──


@pytest.mark.foundation
def test_a_speed_is_a_scalar_or_one_family_shared_by_every_line() -> None:
    """[SP1, the elegance review's B1] ``TraversalRule``'s speed is a scalar or ONE family axis ``(F,)`` shared by every
    line: a family of length F on F lines is a family (the functionals lead with it, ``(F, F, 2, N)``, and member j
    reads every line at speed v_j: its outflow on line l equals the single-speed rule's); a speed of two axes is
    refused by its fragment, as are a non-positive and a non-finite speed. First reds: the ``ndim > 1`` refusal
    removed (a ``(F, 1)`` speed read as per-line, told apart by shape alone); F speeds on F lines read per line (the
    family leg)."""
    posed = _geometry("cylinder_solid")
    basis = PanelBasis.of(ConcentricPartition.of(posed.geometry), *_RES)
    walls = Walls.of(*posed)
    domain = basis.regions.chart.line_domain()
    lines = domain.lines(np.array([[0.3, np.pi / 2], [0.7, np.pi / 2], [1.6, np.pi / 2]]))
    speeds = np.array([1.0, 1.25, 10.0])
    family = TraversalRule.of(lines, basis, walls, np.asarray(_SIGMA[0]), speed=speeds).outflow()
    assert family.shape == (3, 3, 2, basis.size), family.shape
    for j, v in enumerate(speeds):
        alone = TraversalRule.of(lines, basis, walls, np.asarray(_SIGMA[0]), speed=float(v)).outflow()
        np.testing.assert_allclose(family[j], alone, rtol=0.0, atol=1e-13 * np.max(np.abs(alone)))
    for bad in (speeds[:, None], np.array([1.0, 0.0, 2.0]), np.array([1.0, np.inf, 2.0])):
        with pytest.raises(ValueError, match=r"a scalar or a family \(F,\) shared by every line"):
            TraversalRule.of(lines, basis, walls, np.asarray(_SIGMA[0]), speed=bad)


@pytest.mark.foundation
@pytest.mark.rests_on("tests/gates/geometry/test_chord_refinement.py::test_a_partition_coarsens_only_onto_one_it_refines")
def test_a_rule_reads_only_a_chord_through_a_refinement_of_its_basis() -> None:
    """[RF5, the refinement at the rule] ``TraversalRule.on_chord`` takes a chord through the basis's panels or any
    refinement of them (positive leg: refined at a level); a chord through a partition that drops a panel breakpoint, or
    ends elsewhere, is refused by the kernel's fragment (``coarsening``, evaluated in ``__post_init__``). First red: the
    refusal removed (the rule reads the foreign partition's regions as panels)."""
    case = _Case("sphere_solid", *_in_plane(0.7), sigma_t=_SIGMA[0])
    panels = case.basis.partition
    refined = panels.refined_at(1.1)
    TraversalRule.on_chord(refined.chord(case.line), case.basis, case.walls, case.sigma_t)
    coarser = ConcentricPartition(panels.chart, (panels.breakpoints[0], *panels.breakpoints[2:]))
    longer = ConcentricPartition(panels.chart, (*panels.breakpoints, panels.breakpoints[-1] + 1.0))
    for foreign in (coarser, longer):
        with pytest.raises(ValueError, match="coarsens only onto one it refines"):
            TraversalRule.on_chord(foreign.chord(case.line), case.basis, case.walls, case.sigma_t)


# ── PP and PL: the pieces, the law's premise, and the law's node counts (P3 step 3f, route B) ────────────────────────
#
# Spec: ``scratch/characteristic_architecture/p3/product/gate_spec.md`` (R4 and REVISION 2). Each line's slots are cut
# by the hp walk (``grading.branch_walk``) and each piece integrated by the Gauss rule of its attenuation at the node
# counts the law derives (``_PieceSet.node_counts``). The rows here read batched lines (``Line.through`` of (1, 3)).

import functools  # noqa: E402

from orpheus.derivations.continuous.characteristic import transport as _transport_module  # noqa: E402
from orpheus.derivations.continuous.characteristic.grading import ALONG_LINE_ERROR, TIGHTNESS  # noqa: E402

_PP_SIGMA = (0.6, 1.3, 0.45)
_PP_LINES = [
    # (id, chart, breakpoints, point, direction, speed): b > 0 (b = 0 has no branch point: DECLARED excluded)
    ("sphere-b0.7", "sphere", (0.0, 0.5, 1.5, 2.0), (0.7, 0.0, 0.0), (0.0, 1.0, 0.0), 1.0),
    ("sphere-b1e-4", "sphere", (0.0, 0.5, 1.5, 2.0), (1e-4, 0.0, 0.0), (0.0, 1.0, 0.0), 1.0),
    ("sphere-b1.49", "sphere", (0.0, 0.5, 1.5, 2.0), (1.49, 0.0, 0.0), (0.0, 1.0, 0.0), 1.0),
    ("hsphere-cavity-b0.005", "sphere", (0.01, 0.5, 1.5, 2.0), (0.005, 0.0, 0.0), (0.0, 1.0, 0.0), 1.0),
    ("hsphere-shell-b0.2", "sphere", (0.01, 0.5, 1.5, 2.0), (0.2, 0.0, 0.0), (0.0, 1.0, 0.0), 1.0),
    ("cyl3d-b0.3-wz0.8", "cylinder", (0.0, 0.5, 1.5, 2.0), (0.3, 0.0, 0.0), (0.0, 0.6, 0.8), 1.0),
    ("cyl-family-b1e-3", "cylinder", (0.0, 0.5, 1.5, 2.0), (1e-3, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 10.0, 1e4)),
    ("hcyl-cavity-b0.004", "cylinder", (0.01, 0.5, 1.5, 2.0), (0.004, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 1e3)),
]


def _pp_rule(chart, breakpoints, point, direction, speed, sigma=_PP_SIGMA, degree: int = 3):
    laws = {"outer": walls_of.VACUUM} | ({"inner": walls_of.VACUUM} if breakpoints[0] > 0 else {})
    geometry = getattr(walls_of, chart)(breakpoints, (0, 1, 2), **laws)
    basis = PanelBasis.of(ConcentricPartition.of(geometry.geometry), degree, 2, 0.5)
    line = Line.through(np.array([point], dtype=float), np.array([direction], dtype=float))
    return TraversalRule.of(line, basis, Walls.of(*geometry), np.asarray(sigma, dtype=float),
                            speed=np.asarray(speed, dtype=float) if np.ndim(speed) else float(speed))


def _branch_points(chart, point, direction):
    """(t*, D) in the line's OWN parameterisation from ``point``: t* = -P.Omega/|P Omega|^2 and D = b/|P Omega|, 60 digits."""
    with mp.workdps(60):
        P = [mp.mpf(repr(float(v))) for v in point]
        O = [mp.mpf(repr(float(v))) for v in direction]
        keep = 3 if chart == "sphere" else 2
        v2 = sum(O[i] ** 2 for i in range(keep))
        t_star = -sum(P[i] * O[i] for i in range(keep)) / v2
        b2 = sum(P[i] ** 2 for i in range(keep)) - sum(P[i] * O[i] for i in range(keep)) ** 2 / v2
        return t_star, mp.sqrt(b2) / mp.sqrt(v2)


@pytest.mark.foundation
@pytest.mark.parametrize(("chart", "breakpoints", "point", "direction", "speed"), [r[1:] for r in _PP_LINES],
                         ids=[r[0] for r in _PP_LINES])
def test_every_piece_is_no_wider_than_its_branch_distance_over_the_tightness(chart, breakpoints, point, direction,
                                                                             speed) -> None:
    """[PP1; foundation, THEOREM: the law's premise] Every live piece [a, c] off an even panel, in the line parameter,
    satisfies TIGHTNESS (c - a) <= dist([a, c], t* +- iD), the branch points computed here from the line's point and
    direction in mpmath (never from the rule's slots), slack 1e-12 relative. A piece's ends are signed distances from
    its slot's anchor (``_Slots.origin``), read in the line's parameter, which is measured from its foot; t* is
    translated by the point's parameter (formed from the same line).
    DECLARED excluded: even panels (a polynomial in t, entire). First red (``carve/arms_*.log``): the walk at half the
    tightness (``walk_loose``)."""
    rule = _pp_rule(chart, breakpoints, point, direction, speed)
    p = rule._pieces.pieces
    start = rule._of_slot(rule._slots.origin, p.slot)            # piece ends are signed distances from the slot's anchor
    live = p.live & ~np.asarray(rule.basis.even)[p.panel]
    a, c = (start + p.lower)[live], (start + p.upper)[live]
    assert a.size >= 3, "no pieces: the row reads nothing"
    t_star, D = _branch_points(chart, point, direction)
    foot = float(np.asarray(rule.period.chord.line.parameter_of(np.array([point], dtype=float)))[0])  # the point's t
    worst = mp.mpf(0)
    with mp.workdps(60):
        for lo, hi in zip(a, c):
            lo_, hi_ = mp.mpf(repr(float(lo))) - mp.mpf(repr(foot)), mp.mpf(repr(float(hi))) - mp.mpf(repr(foot))
            real = max(lo_ - t_star, t_star - hi_, mp.mpf(0))
            worst = max(worst, TIGHTNESS * (hi_ - lo_) / mp.sqrt(real * real + D * D))
    assert worst <= 1 + 1e-12, f"a piece {float(worst):.6g} x (its branch distance / c) wide"


@pytest.mark.foundation
@pytest.mark.parametrize("which", ["sphere-b0.7", "hsphere-cavity-b0.005", "cyl-family-b1e-3"])
@pytest.mark.rests_on(_HERE + "test_every_piece_is_no_wider_than_its_branch_distance_over_the_tightness")
def test_the_pieces_depend_on_the_geometry_alone(which: str) -> None:
    """[PP2; foundation, ROUTE: the attenuation left the piece mesh] The piece set is ``array_equal`` under every cross
    section scaled by 1e-3 and by 1e3 and under a family widened by 1e4: the attenuation is the measure of each piece's
    Gauss rule, never a cut of the line. INVERTED successor of the lift file's G1b (the grading read the family's
    extremes). First red: today's ``main`` before P3 step 3f (`[M]` 2026-10-13: 36 pieces against 272 on
    ``sphere-b0.7``, ``scratch/.../p3/product/ta/runs/pieces_today.log``)."""
    _, chart, breakpoints, point, direction, speed = next(r for r in _PP_LINES if r[0] == which)
    base = _pp_rule(chart, breakpoints, point, direction, speed)._pieces.pieces
    wide = tuple(np.atleast_1d(speed)) + (float(np.max(speed)) * 1e4,)
    for sigma, sp in ((tuple(1e-3 * v for v in _PP_SIGMA), speed), (tuple(1e3 * v for v in _PP_SIGMA), speed),
                      (_PP_SIGMA, wide)):
        other = _pp_rule(chart, breakpoints, point, direction, sp, sigma)._pieces.pieces
        for name in ("slot", "lower", "upper", "panel", "live"):
            np.testing.assert_array_equal(getattr(base, name), getattr(other, name), err_msg=name)


#: PL1's fixtures: (id, chart, breakpoints, point, direction, speed, degree, margin). The margin is the law's excess
#: over the least shift that holds the target, measured on this fixture (`[M]` 2026-10-13, ``carve/probe_pl1.log``),
#: within ``margin_B.log``'s c = 2 range (2 to 3); re-measured at the 64-eps target (``carve/probe_pl1_64.log``).
_PL1 = [
    ("sphere-turning-b1e-4-p3", "sphere", (0.0, 0.5, 1.5, 2.0), (1e-4, 0.0, 0.0), (0.0, 1.0, 0.0), 1.0, 3, 3),
    ("sphere-b0.7-p5", "sphere", (0.0, 0.5, 1.5, 2.0), (0.7, 0.0, 0.0), (0.0, 1.0, 0.0), 1.0, 5, 1),
    ("hsphere-cavity-b0.005-p3", "sphere", (0.01, 0.5, 1.5, 2.0), (0.005, 0.0, 0.0), (0.0, 1.0, 0.0), 1.0, 3, 2),
    ("cyl-family-b0.3-p1", "cylinder", (0.0, 0.5, 1.5, 2.0), (0.3, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 10.0, 1e3), 1, 0),
    ("cyl-family-b0.3-p3", "cylinder", (0.0, 0.5, 1.5, 2.0), (0.3, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 10.0, 1e3), 3, 0),
]
#: The law's target (``grading.ALONG_LINE_ERROR``, the user's ruling of 2026-10-13: the rule's own band, 64 eps).
_PL1_TARGET = ALONG_LINE_ERROR


#: The law, captured once at import (every shifted property below wraps THIS, never a property already shifted).
_HONEST_COUNTS = _transport_module._PieceSet.__dict__["node_counts"].func


def _shifted_counts(shift: int):
    """``_PieceSet.node_counts`` with ``shift`` added to every walked piece's counts (at least 1): the law moved."""
    honest = _HONEST_COUNTS

    def counts(self):
        one, two = honest(self)
        rule, p = self.rule, self.pieces
        walked = p.live & rule._of_slot(rule._walked, p.slot)
        return np.where(walked, np.maximum(one + shift, 1), one), np.where(walked, np.maximum(two + shift, 1), two)

    prop = functools.cached_property(counts)
    prop.__set_name__(_transport_module._PieceSet, "node_counts")
    return prop


def _pl1_error(monkeypatch, row, shift: int, reference) -> float:
    _, chart, breakpoints, point, direction, speed, degree, _ = row
    monkeypatch.setattr(_transport_module._PieceSet, "node_counts", _shifted_counts(shift))
    rule = _pp_rule(chart, breakpoints, point, direction, speed, degree=degree)
    weight = np.ones((np.size(speed), 1)) if np.ndim(speed) else np.ones(1)
    B, V = rule.outflow(), rule.volterra(weight)
    return max(float(np.max(np.abs(B - reference[0])) / np.max(np.abs(reference[0]))),
               float(np.max(np.abs(V - reference[1])) / np.max(np.abs(reference[1]))))


@pytest.mark.foundation
@pytest.mark.parametrize("row", _PL1, ids=[r[0] for r in _PL1])
@pytest.mark.rests_on(_HERE + "test_every_piece_is_no_wider_than_its_branch_distance_over_the_tightness",
                      "tests/gates/derivations/test_characteristic_attenuation_rule.py::"
                      "test_the_rule_is_exact_on_every_legendre_polynomial_below_degree_two_m")
def test_the_laws_node_counts_meet_the_target_and_exceed_the_least_by_at_most_the_margin(monkeypatch, row) -> None:
    """[PL1; foundation, REFERENCE by self-convergence] On each fixture line, the outflow and the Volterra block at the
    law's node counts against the same rule at the counts + 3 (on the walked pieces; an entire piece is exact at its
    degree's count): SAFE, within the law's target (``ALONG_LINE_ERROR``, 64 eps) of the largest entry; TIGHT, the
    least shift from which every shift up holds the target is at most the row's margin below the law. No mpmath line reference resolves 1e-14 (the transport rows'
    bands are 1e-12 and 1e-13), so the claim is convergence in the count, its value anchored by RB3, PT1 and the
    mpmath rows. Arms (the witness of each leg): the law minus (margin + 1) reds SAFE; the law + 3 (a law grown
    over-conservative) reds TIGHT (``carve/arms_*.log``)."""
    monkeypatch.setattr(_transport_module._PieceSet, "node_counts", _shifted_counts(3))
    _, chart, breakpoints, point, direction, speed, degree, margin = row
    rule = _pp_rule(chart, breakpoints, point, direction, speed, degree=degree)
    weight = np.ones((np.size(speed), 1)) if np.ndim(speed) else np.ones(1)
    reference = (rule.outflow(), rule.volterra(weight))
    shift = int(__import__("os").environ.get("PL1_SHIFT", "0"))         # the arms: the law moved by PL1_SHIFT
    errors = {j: _pl1_error(monkeypatch, row, shift + j, reference) for j in range(0, -(margin + 2), -1)}
    assert errors[0] <= _PL1_TARGET, ("SAFE", errors)
    least = min(j for j in errors if all(errors[k] <= _PL1_TARGET for k in errors if k >= j))
    assert -least <= margin, ("TIGHT", least, margin, errors)


@pytest.mark.foundation
@pytest.mark.parametrize("which", ["sphere-b0.7", "hsphere-shell-b0.2", "cyl3d-b0.3-wz0.8", "cyl-family-b1e-3"])
def test_an_unbatched_line_is_the_same_line_in_a_batch_of_one(which: str) -> None:
    """[UB1, qa's row; foundation, THEOREM] A line built from 1-D arrays (its image's fields numpy scalars) has the
    outflow, the entry response and the Volterra block of the same line in a batch of one, ``array_equal`` (the
    batch axis dropped), on each radial body. The carve's defect 1 (2026-10-13): ``Chord._per_point`` spread only
    ndarray fields, and an unbatched line failed in ``rays_at``. First red: ``Chord.take`` gathering ndarray fields
    only (``kg_arms.py``, ``take_arrays_only``)."""
    _, chart, breakpoints, point, direction, speed = next(r for r in _PP_LINES if r[0] == which)
    laws = {"outer": walls_of.VACUUM} | ({"inner": walls_of.VACUUM} if breakpoints[0] > 0 else {})
    geometry = getattr(walls_of, chart)(breakpoints, (0, 1, 2), **laws)
    basis = PanelBasis.of(ConcentricPartition.of(geometry.geometry), 3, 2, 0.5)
    sp = np.asarray(speed, dtype=float) if np.ndim(speed) else float(speed)

    def rule(line):
        return TraversalRule.of(line, basis, Walls.of(*geometry), np.asarray(_PP_SIGMA), speed=sp)

    single = rule(Line.through(np.array(point, dtype=float), np.array(direction, dtype=float)))
    batch = rule(Line.through(np.array([point], dtype=float), np.array([direction], dtype=float)))
    family = np.ndim(speed)
    one = (slice(None), 0) if family else (0,)
    np.testing.assert_array_equal(single.outflow(), batch.outflow()[one])
    np.testing.assert_array_equal(single.entry_response(), batch.entry_response()[one])
    weights = np.ones(np.shape(speed)) if family else np.array(1.0)
    np.testing.assert_array_equal(single.volterra(weights), batch.volterra(weights[:, None] if family else np.ones(1)))
