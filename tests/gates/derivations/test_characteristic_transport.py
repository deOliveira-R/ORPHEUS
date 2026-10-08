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
from orpheus.geometry.boundary import AlbedoBoundary, PeriodicBoundary, SpecularReturn, VacuumInflow
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.line import Line
from orpheus.geometry.structured_geometry import StructuredGeometry

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
_POINTS, _INNER = 16, 16
_RES = (3, 2, 0.5)                 # (degree, layers, ratio): the rung's working point


def _spec(a: float):
    return AlbedoBoundary(a, SpecularReturn(axis="x"))


def _geometry(body: str, *, a_in: float = _A_IN, a_out: float = _A_OUT, breakpoints=None) -> StructuredGeometry:
    match body:
        case "sphere_solid":
            return StructuredGeometry.sphere(breakpoints or _MR3, (0, 1, 2), outer=_spec(a_out))
        case "cylinder_solid":
            return StructuredGeometry.cylinder(breakpoints or _MR3, (0, 1, 2), outer=_spec(a_out))
        case "sphere_hollow":
            return StructuredGeometry.sphere(breakpoints or _MR3H, (0, 1, 2), inner=_spec(a_in), outer=_spec(a_out))
        case "cylinder_hollow":
            return StructuredGeometry.cylinder(breakpoints or _MR3H, (0, 1, 2), inner=_spec(a_in), outer=_spec(a_out))
        case "slab":
            return StructuredGeometry.slab(breakpoints or _SLB3, (0, 1, 2), left=_spec(a_in), right=_spec(a_out))
        case "slab_periodic":
            return StructuredGeometry.slab(breakpoints or _SLB3, (0, 1, 2), left=PeriodicBoundary(axis="x"),
                                           right=PeriodicBoundary(axis="x"))
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
                 points: int = _POINTS, inner: int = _INNER) -> None:
        self.body = body
        self.geometry = geometry or _geometry(body)
        self.basis = PanelBasis.of(ConcentricPartition.of(self.geometry), *resolution)
        line = Line.through(np.asarray(point, dtype=float), np.asarray(direction, dtype=float))
        self.sigma_t = np.asarray(sigma_t, dtype=float)
        self.rule = TraversalRule.of(line, self.basis, Walls.of(self.geometry), self.sigma_t, points, inner)
        self.period = self.rule.period
        self.mp = R.MpLine(_chart(body), self.geometry.breakpoints, point, direction,
                           extra_radii=self.basis.partition.breakpoints)

    def coefficients(self, q) -> np.ndarray:
        return R.sample(R.region_polynomial(q), np.asarray(self.basis.region), np.asarray(self.basis.nodes))


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
        Line.through(np.asarray(point, dtype=float), np.asarray(direction, dtype=float))), Walls.of(case.geometry))
    np.testing.assert_allclose(case.rule.optical_depth, coarse.optical_depth(case.sigma_t), rtol=16 * _EPS, atol=0.0)


# ── T1-T3: the traversal integrals ───────────────────────────────────────


def _traversals(case: _Case, expected):
    return [(k, case.mp.traversal(transit, rev)) for k, (transit, rev, _end) in enumerate(expected)]


_T1_TOL = 1e-13


@pytest.mark.l0
@pytest.mark.verifies("characteristic-traversal-integrals")
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
    geometry = StructuredGeometry.slab((0.0, 1.0), (0,), left=VacuumInflow(), right=VacuumInflow())
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
    geometry = StructuredGeometry.slab((0.0, 1.3), (0,), left=VacuumInflow(), right=VacuumInflow())
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
    laws = {"outer": VacuumInflow(), **({"inner": VacuumInflow()} if hollow else {})}
    geometry = getattr(StructuredGeometry, chart)(breakpoints, tuple(range(n)), **laws)
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
def test_the_turning_grading_halves_toward_the_closest_approach_at_seven_points(body, point, direction) -> None:
    """[T9b] At 7 points per piece, B_0 . q and psi(t) . q on a turning line stay within 1e-13 of mpmath (q = c).

    The witness for the grading's RATIO, which the 16-point rows cannot see:
    pieces shrinking by 4 instead of 2 toward the branch point leave each one
    too close to it for 7 Gauss points. Same 1e-13 band as T1. First red: the
    ratio coarsened to 1/4 (battery arm T3b). Re-posed 2026-10-06 (rung 3, RB5)
    from a solid body's centre (now the even panel, polynomial in arc length) to
    a hollow body with a cavity of radius 1e-5; which rows the arm reds there,
    and so which are controls, is in the rung-3 battery
    (``scratch/characteristic_architecture/p1_step_b3/gates/battery/``).
    """
    geometry = _geometry(body, breakpoints=_CAVITY)
    case = _Case(body, point, direction, sigma_t=_SIGMA[0], points=7, inner=7, geometry=geometry)
    assert int(case.period.rank) == 1, "the line misses the cavity"
    q = [[0.0, 1.0]] * 3
    x = case.coefficients(q)
    f = R.region_polynomial(q)
    trav = case.mp.traversal(0, False)
    _close(float(case.rule.outflow()[0] @ x), R.outflow(case.mp, trav, f, _SIGMA[0]), _T1_TOL, "B_0 at 7 points")
    t, ts = _points_on(case, 0)
    psi = case.rule.angular_flux(t, np.zeros((2, case.basis.size))) @ x
    for i, tm in enumerate(ts):
        _close(float(psi[i]), R.partial_outflow(case.mp, trav, f, _SIGMA[0], tm), _T1_TOL, f"psi(t_{i}) at 7 points")


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
    geometry = _geometry("cylinder_solid")
    points = np.array([(0.7, 0.0, 0.0), (0.7, 0.0, 0.0), (0.0, 0.0, 0.0), (1.2, 0.0, 0.0)])
    directions = np.array([(0.0, 0.0, 1.0), (0.0, 0.6, 0.8), (0.0, 1.0, 0.0), (0.0, 0.8, -0.6)])
    weights = np.array([2.5, 0.3, 1.7, 0.0])
    basis = PanelBasis.of(ConcentricPartition.of(geometry), *_RES)
    walls = Walls.of(geometry)

    def rule(p, d):
        return TraversalRule.of(Line.through(p, d), basis, walls, np.asarray(_SIGMA[0]), _POINTS, _INNER)

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
    geometry = StructuredGeometry.slab((0.0, 1.0), (0,), left=VacuumInflow(), right=VacuumInflow())
    case = _Case("slab", (0.5, 0.0, 0.0), (1.0, 0.0, 0.0), sigma_t=(sigma,), geometry=geometry, resolution=(2, 0, 0.5))
    assert case.basis.n_panels == 1
    f, g = [0.3, 1.0, -0.7], [1.1, -0.4, 0.9]
    xf, xg = case.coefficients([f]), case.coefficients([g])
    S = R.mpf(sigma)
    mfp = 1.0 / sigma
    t = np.array([0.5, 1.0 - 64.0 * mfp, 1.0 - 64.0 * mfp + 0.25 * mfp, 1.0 - 0.5 * mfp])
    psi = case.rule.angular_flux(t, np.zeros((2, case.basis.size))) @ xg
    for i, ti in enumerate(t):
        _close(float(psi[i]), _attenuated_to(g, S, R.mpf(ti)), _T1_TOL, f"psi at t = {ti}")
    V = case.rule.volterra(np.array(1.0))
    _close(float(xf @ V @ xg), _volterra_closed(f, g, S, mp.mpf(1)), _T1_TOL, "f . V g")


# ── T10-T12: the angular flux along the line ─────────────────────────────

_FRACTIONS = (0.013, 0.37, 0.5, 0.91, 1.0 - 1e-9)


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
        t, ts = _points_on(case, transit)
        psi = case.rule.angular_flux(t, zero)
        assert psi.shape == (len(t), case.basis.size)
        trav = case.mp.traversal(transit, False)
        for i, tm in enumerate(ts):
            _close(float(psi[i] @ x), R.partial_outflow(case.mp, trav, f, _SIGMA[1], tm), _T1_TOL, f"psi(t_{i})")
    outside = float(case.mp.transits()[-1][1]) + 0.1
    with pytest.raises(ValueError, match="lies on a transit of its line"):
        case.rule.angular_flux(np.array([outside]), zero)
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
]


@pytest.mark.l0
@pytest.mark.verifies("characteristic-closure")
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
    """
    case = _Case(body, point, direction, sigma_t=_SIGMA[0])
    x = case.coefficients(_Q)
    rule = case.rule
    inflow = case.period.inflow(rule.optical_depth, rule.outflow())
    f = R.region_polynomial(_Q)
    for k0, (transit, rev, _end) in enumerate(period_by_hand):
        if rev:
            continue
        t, ts = _points_on(case, transit)
        psi = rule.angular_flux(t, inflow) @ x
        for i, tm in enumerate(ts):
            _close(float(psi[i]), _unfolded_psi(case, period_by_hand, amplitudes, f, _SIGMA[0], k0, tm), _T1_TOL,
                   f"psi(t_{i}) on traversal {k0}")


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
    t, _ = _points_on(case, 0)
    np.testing.assert_array_equal(rule.angular_flux(t, inflow), rule.angular_flux(t, np.zeros_like(inflow)))
