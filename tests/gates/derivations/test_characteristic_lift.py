"""Gates for the cylinder's axial lift: the in-plane line at speed 1/sin(theta) (P3 step 3 of #405).

The carve (``.claude/plans/characteristic_reference_architecture.md``, "P3 opened (2026-10-11)", step 3; verification
spec ``scratch/characteristic_architecture/p3/step3/gate_spec.md``, the user's API ruling of 2026-10-11): a cylinder
line of impact parameter b at polar angle theta is its in-plane line (b, pi/2) traversed at 3-D speed
``ds/dt = 1/sin(theta)``. :class:`~orpheus.derivations.continuous.characteristic.transport.TraversalRule` carries the
speed as a field whose leading axes are a FAMILY over the same lines, so each in-plane line is chorded once and theta
is a broadcast axis of the exponentials and of the closure a B / (1 - a e^{-kappa T}); the along-line pieces are the
geometry's alone (the hp walk), each member integrating its attenuation by its own Gauss rule (P3 step 3f); the impact rule is
graded at the slowest polar node (#587's per-polar grading retires).

The ladder, bottom up (``rests_on`` on each row):

1. the per-line transport on 3-D chords (``test_characteristic_transport.py``), its references written in
   ``_characteristic_mp.py`` (mpmath; a cylinder line built from a 3-vector direction, no stretch assumed);
2. the family against the 3-D line [K1a outflow and entry response, K1b the vacuum angular flux, K1c the closed
   angular flux under a partial mirror, K1d the Volterra triangle]; theta a pure broadcast axis [K2]; the stretch
   carried by the rule's speed equals the stretch carried by the kernel's 3-D chord [K3];
3. the polar fold of one in-plane line against Bickley's functions [K4];
4. the dyadic grading of the rules over lines [G1, the end law; G1c, its callers]; the along-line rule's value leg on
   the closed-form readings [G2] (the rule itself: ``test_characteristic_attenuation_rule.py``);
5. the readings against frozen mpmath routes [F1, F1b, F4; F2, F3 are G2], each literal recomputed by a slow twin;
6. #587's retirement loses no accuracy, per polar node [P1];
7. the sphere and the slab, moved by the target-graded reading, accepted on the three conditions [S1].

Every reference is mpmath from ``_characteristic_mp.py`` / ``_characteristic_point_mp.py``; no reference imports the
code under test. Prototype measurements behind the bands: ``scratch/characteristic_architecture/p3/step3/ta/``.
"""
from __future__ import annotations

import functools
import math

import mpmath as mp
import numpy as np
import pytest

from orpheus.derivations.continuous.characteristic import PanelBasis, Walls
from orpheus.derivations.continuous.characteristic import reading as reading_module
from orpheus.derivations.continuous.characteristic.grading import DYADIC, VANISHING_DEPTH, ExponentialGrading, exponential_ends
from orpheus.derivations.continuous.characteristic.lines import OpticalScale, polar_rule
from orpheus.derivations.continuous.characteristic.reading import PointRule
from orpheus.derivations.continuous.characteristic.transport import TraversalRule
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.line import Line

from . import _characteristic_mp as R
from . import _characteristic_point_mp as mpref
from . import test_characteristic_assembly as AS
from . import test_characteristic_reading as RE
from . import test_characteristic_transport as T

pytestmark = pytest.mark.filterwarnings("error::RuntimeWarning")

_HERE = "tests/gates/derivations/test_characteristic_lift.py::"
_TRANS = "tests/gates/derivations/test_characteristic_transport.py::"
_ASSEMBLY = "tests/gates/derivations/test_characteristic_assembly.py::"
_READING = "tests/gates/derivations/test_characteristic_reading.py::"

#: The polar SET of the kernel rows, as sin(theta): the normal direction, two oblique ones and three toward grazing
#: (the activation of the broadcast: one member reduces every row to its per-line twin). The family's speed is
#: 1 / sin(theta), spelled from these values directly.
_SINES = (1.0, 0.6, 0.1, 1e-2, 1e-3, 1e-4)


def _direction(sine: float) -> tuple[float, float, float]:
    """The 3-D direction of polar sine ``sine`` whose in-plane part is +y (the transport file's ``_in_plane``)."""
    return 0.0, sine, math.sqrt((1.0 - sine) * (1.0 + sine))


class _Family:
    """One in-plane cylinder line (b, pi/2) at a family of speeds 1/sin(theta), and its 3-D lines in mpmath."""

    def __init__(self, body: str, b: float, sines, sigma_t, *, resolution=T._RES) -> None:
        self.body = body
        self.geometry = T._geometry(body)
        self.basis = PanelBasis.of(ConcentricPartition.of(self.geometry), *resolution)
        self.walls = Walls.of(self.geometry)
        self.sigma_t = np.asarray(sigma_t, dtype=float)
        self.b = b
        self.sines = tuple(float(s) for s in sines)
        self.line = Line.through(np.array([b, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))
        self.speed = 1.0 / np.asarray(self.sines)
        self.rule = TraversalRule.of(self.line, self.basis, self.walls, self.sigma_t, speed=self.speed)
        self.planar = R.MpLine("cylinder", self.geometry.breakpoints, (b, 0.0, 0.0), (0.0, 1.0, 0.0),
                               extra_radii=self.basis.partition.breakpoints)

    def case(self, j: int) -> T._Case:
        """The 3-D line of member ``j``: the transport file's case (its mpmath line, and the rule on the 3-D chord)."""
        return T._Case(self.body, (self.b, 0.0, 0.0), _direction(self.sines[j]), sigma_t=self.sigma_t)

    def coefficients(self, q) -> np.ndarray:
        return R.sample(R.region_polynomial(q), np.asarray(self.basis.region), np.asarray(self.basis.nodes))

    def readings(self, transit: int):
        """Per interior point of a transit (fractions :data:`_INTERIOR` of its in-plane length), the family's rule chorded
        through the basis's panels refined at the point's level (``ConcentricPartition.refined_at``, the target-graded reading) and
        the point's in-plane parameter as that chord's own crossing (a reading is read at a crossing on a transit)."""
        t_in, t_out, _ = self.planar.transits()[transit]
        for fraction in _INTERIOR:
            t = float(t_in + (t_out - t_in) * mp.mpf(fraction))
            level = math.hypot(self.b, t)
            chord = self.basis.partition.refined_at(level).chord(self.line)
            rule = TraversalRule.on_chord(chord, self.basis, self.walls, self.sigma_t, speed=self.speed)
            yield rule, _crossing_near(rule, t)


#: Interior points of a transit as fractions of its length, the closest approach (0.5, tangent to its own level) excluded.
_INTERIOR = (0.013, 0.37, 0.63, 0.91)


def _crossing_near(rule: TraversalRule, t: float) -> float:
    """The present crossing parameter of ``rule``'s chord nearest ``t``: the point is read where the kernel put it."""
    crossings = rule.period.chord.crossings
    parameters = np.asarray(crossings.parameter, dtype=float)[np.asarray(crossings.present, dtype=bool)]
    return float(parameters[np.argmin(np.abs(parameters - t))])


def _crossings_on_transit(rule: TraversalRule, transit: int) -> np.ndarray:
    """The crossing parameters of ``rule``'s (unrefined) chord from a transit's first slot start to its exit."""
    chord = rule.period.chord
    first, stop = int(chord.transits.first_slot[transit]), int(chord.transits.stop_slot[transit])
    return np.asarray(chord.crossings.parameter, dtype=float)[first:stop + 1]


#: (id, body, b, the period by hand). The transport file's cylinder lines, in plane.
_K1_LINES = [
    ("solid_b0.7", "cylinder_solid", 0.7, [(0, False, "last")]),
    ("solid_b0.3", "cylinder_solid", 0.3, [(0, False, "last")]),
    ("hollow_b0.2", "cylinder_hollow", 0.2, [(0, False, "first"), (1, False, "last")]),
]


# ── K1: the family against the 3-D lines ─────────────────────────────────


@pytest.mark.l0
@pytest.mark.verifies("characteristic-traversal-integrals")
@pytest.mark.parametrize("group", [0, 1])
@pytest.mark.parametrize(("body", "b", "expected"), [r[1:] for r in _K1_LINES], ids=[r[0] for r in _K1_LINES])
@pytest.mark.rests_on(_TRANS + "test_the_outflow_of_a_per_region_polynomial_is_its_line_integral_attenuated_to_the_exit",
                      _TRANS + "test_the_entry_response_is_the_integral_attenuated_from_the_entry")
def test_each_speed_of_the_family_is_its_3d_lines_outflow_and_entry_response(body, b, expected, group) -> None:
    """[K1a; l0] Member j of the family (the in-plane line at speed 1/sin(theta_j)) has the outflow B_k and the entry
    response A_k of the 3-D line (b, theta_j), against mpmath on that 3-D line (``MpLine`` from a 3-vector direction:
    its speed is the norm of the direction's kept components, so it assumes no stretch), per traversal, 1e-13
    relative, sin(theta) from 1 to 1e-4 in ONE call.

    `[M]` 2026-10-11 (``ta/mpline_probe.log``): today's per-(b, theta) rule meets this reference to <= 7.5e-15 at
    (16, 16) points down to sin(theta) = 1e-4. First reds: the speed dropped from the length weight (``speed * dt``
    -> ``dt``: every member but sin = 1); the speed applied to the optical depth only, or to the weight only.
    """
    family = _Family(body, b, _SINES, T._SIGMA[group])
    x = family.coefficients(T._q(body))
    f = R.region_polynomial(T._q(body))
    B, A = family.rule.outflow(), family.rule.entry_response()
    assert B.shape == A.shape == (len(_SINES), 2, family.basis.size), B.shape
    for j in range(len(_SINES)):
        case = family.case(j)
        for k, trav in T._traversals(case, expected):
            T._close(float(B[j, k] @ x), R.outflow(case.mp, trav, f, T._SIGMA[group]), T._T1_TOL, f"B_{k}, sin {_SINES[j]}")
            T._close(float(A[j, k] @ x), R.entry_response(case.mp, trav, f, T._SIGMA[group]), T._T1_TOL,
                     f"A_{k}, sin {_SINES[j]}")


@pytest.mark.l0
@pytest.mark.verifies("characteristic-traversal-integrals")
@pytest.mark.parametrize(("body", "b", "expected"), [r[1:] for r in _K1_LINES], ids=[r[0] for r in _K1_LINES])
@pytest.mark.rests_on(_HERE + "test_each_speed_of_the_family_is_its_3d_lines_outflow_and_entry_response",
                      _TRANS + "test_the_vacuum_angular_flux_is_the_source_integral_since_the_entry")
def test_each_speed_of_the_family_reads_its_3d_lines_vacuum_angular_flux(body, b, expected) -> None:
    """[K1b; l0] The family's vacuum angular flux at in-plane parameters t (four interior fractions of each forward
    transit, each read where the line is chorded through the panels refined at the point's level, the target-graded
    reading's route) is, for member j, the mpmath source integral since the entry on the 3-D line (b, theta_j) at its
    parameter t / sin(theta_j) (both lines' feet are (b, 0, 0)), 1e-13 relative.

    First reds: the speed dropped from the length weight; kappa inverted in the exponent (`[M]` on the prototype
    reading, ``ta/proto_battery.log``: +2.4e-2 to +4.2); the point parameter read in 3-D units (t3 for t2).
    """
    family = _Family(body, b, _SINES, T._SIGMA[1])
    x = family.coefficients(T._Q)
    f = R.region_polynomial(T._Q)
    cases = [family.case(j) for j in range(len(_SINES))]
    for transit, rev, _end in expected:
        if rev:
            continue
        for i, (rule, t) in enumerate(family.readings(transit)):
            psi = rule.angular_flux(np.array([t]), np.zeros_like(rule.outflow()))
            assert psi.shape == (len(_SINES), 1, family.basis.size), psi.shape
            for j, sine in enumerate(_SINES):
                trav = cases[j].mp.traversal(transit, False)
                want = R.partial_outflow(cases[j].mp, trav, f, T._SIGMA[1], mp.mpf(repr(t)) / mp.mpf(repr(sine)))
                T._close(float(psi[j, 0] @ x), want, T._T1_TOL, f"psi(t_{i}), sin {sine}")


#: (id, body, b, period by hand, the amplitude of each traversal's exit wall): the transport file's T11 cylinder rows.
_K1C = [
    ("solid_partial", "cylinder_solid", 0.7, [(0, False, "last")], (T._A_OUT,)),
    ("hollow_cavity", "cylinder_hollow", 0.2, [(0, False, "first"), (1, False, "last")], (T._A_IN, T._A_OUT)),
]


@pytest.mark.l0
@pytest.mark.verifies("characteristic-closure")
@pytest.mark.parametrize(("body", "b", "period_by_hand", "amplitudes"), [r[1:] for r in _K1C], ids=[r[0] for r in _K1C])
@pytest.mark.rests_on(_HERE + "test_each_speed_of_the_family_reads_its_3d_lines_vacuum_angular_flux",
                      _TRANS + "test_the_closed_angular_flux_is_the_unfolded_backward_path")
def test_each_speed_of_the_family_reads_its_3d_lines_closed_angular_flux(body, b, period_by_hand, amplitudes) -> None:
    """[K1c; l0, the fast catcher of a closure exact only at a in {0, 1}] Under partial mirrors (0.6 out, 0.3 in), the
    family's inflow ``period.inflow(optical_depth, outflow)`` carried to in-plane parameters equals, for member j, the
    backward path marched wall by wall in mpmath on the 3-D line (b, theta_j) (``_unfolded_psi``), 1e-13 relative.

    The closure is a B / (1 - a e^{-kappa T}) per member. The step-3 prototype spelled it a B / (1 - e^{-kappa T}),
    exact at a = 0 and a = 1 only; on a homogeneous cylinder's wall it read +1.1e-2 at a = 0.9 (`[M]` 2026-10-11,
    ``ta/floor_probe3.py``), and no fast row posed 0 < a < 1 on a cylinder (``ta/census_slowest_not_slow.json``):
    this row and F1 are its catchers. First reds: that closure; the closure's kappa taken from one member
    (``kappa.max()``: `[M]` prototype -2.7e-2); B without the speed.
    """
    family = _Family(body, b, _SINES, T._SIGMA[0])
    x = family.coefficients(T._Q)
    f = R.region_polynomial(T._Q)
    cases = [family.case(j) for j in range(len(_SINES))]
    for k0, (transit, rev, _end) in enumerate(period_by_hand):
        if rev:
            continue
        for i, (rule, t) in enumerate(family.readings(transit)):
            inflow = rule.period.inflow(rule.optical_depth, rule.outflow())
            psi = rule.angular_flux(np.array([t]), inflow) @ x
            for j, sine in enumerate(_SINES):
                want = T._unfolded_psi(cases[j], period_by_hand, amplitudes, f, T._SIGMA[0], k0,
                                       mp.mpf(repr(t)) / mp.mpf(repr(sine)))
                T._close(float(psi[j, 0]), want, T._T1_TOL, f"psi(t_{i}) on traversal {k0}, sin {sine}")


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_each_speed_of_the_family_is_its_3d_lines_outflow_and_entry_response",
                      _TRANS + "test_the_volterra_triangle_is_the_double_integral_along_the_line")
def test_the_familys_volterra_triangle_is_the_weighted_sum_of_its_3d_lines() -> None:
    """[K1d] ``volterra(w)`` with one weight per member is sum_j w_j f . V_j g, V_j the triangle of the 3-D line
    (b, theta_j) in 3-D arc length (mpmath ``volterra_bilinear``, two Gauss degrees agreeing to 1e-18 first),
    1e-12 relative. The block's polar fold of one in-plane line is this sum with w_j = w_theta 2 sin^2(theta_j).

    First reds: the speed applied once instead of twice (the double integral in 3-D length is kappa^2 times the
    in-plane one); the weights applied to the wrong axis (a member's weight on another member).
    """
    sines, weights = (1.0, 0.6, 0.1), np.array([0.3, 1.7, 0.9])
    family = _Family("cylinder_solid", 0.7, sines, T._SIGMA[0])
    f = R.region_polynomial(T._Q)
    g = R.region_polynomial([[0.4, -0.1, 0.2, 0.0], [2.0, 0.3, -0.5, 0.1], [0.7, 0.2, 0.3, -0.05]])
    xf, xg = family.coefficients(T._Q), family.coefficients([[0.4, -0.1, 0.2, 0.0], [2.0, 0.3, -0.5, 0.1], [0.7, 0.2, 0.3, -0.05]])
    V = family.rule.volterra(weights)
    assert V.shape == (family.basis.size, family.basis.size)
    want = mp.mpf(0)
    for j, w in enumerate(weights):
        case = family.case(j)
        pair = [R.volterra_bilinear(case.mp, case.mp.traversal(0, False), f, g, T._SIGMA[0], degree=d) for d in (4, 5)]
        assert abs(pair[0] - pair[1]) <= mp.mpf("1e-18") * abs(pair[1]), "the mpmath reference has not converged"
        want += mp.mpf(repr(float(w))) * pair[1]
    T._close(float(xf @ V @ xg), want, 1e-12, "f . V g over the family")


# ── K2: theta is a broadcast axis ────────────────────────────────────────


def _columns(body: str, sines) -> dict[str, np.ndarray]:
    family = _Family(body, 0.7, sines, T._SIGMA[0])
    rule = family.rule
    inflow = rule.period.inflow(rule.optical_depth, rule.outflow())
    t = _crossings_on_transit(rule, 0)
    return {"outflow": rule.outflow(), "entry_response": rule.entry_response(), "optical_depth": rule.optical_depth,
            "angular_flux": rule.angular_flux(t, inflow)}


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_each_speed_of_the_family_reads_its_3d_lines_closed_angular_flux")
def test_a_member_added_to_the_family_changes_no_other_member() -> None:
    """[K2; foundation, THEOREM] Theta is a pure broadcast axis: adding a member to the family, between its slowest and
    fastest, slower than its slowest, or first in the family's order, leaves every other member's outflow, entry response,
    optical depth and closed angular flux ``array_equal``. Since P3 step 3f the pieces are the geometry's alone and
    each member integrates its own attenuation by its own Gauss rule, so no member reads another (PROMOTED: until 3f
    the along-line grading read the family's extremes, and only an interior insertion was an invariance). The
    activation leg: two members' columns differ (the broadcast is not a copy).

    Under a partial mirror (0.6), so the closure is in every column. Declared blind (Mode 12, `[M]` qa 2026-10-12,
    arm ``closure_kmax``): an optical depth read at the family's max(speed) for every member when the insertion leaves
    the max unchanged (the interior and first-placed legs); the slower leg moves the max and sees it. First red: the
    member's rule read at the family's slowest speed (``carve/arms_*.log``, arm ``family_slowest``).
    """
    base = (1.0, 0.6, 1e-3)
    honest = _columns("cylinder_solid", base)
    assert not np.array_equal(honest["outflow"][0], honest["outflow"][1]), "two members read the same column"
    for added, keep in (((1.0, 0.6, 0.1, 1e-3), [0, 1, 3]), (base + (1e-4,), [0, 1, 2]), ((0.99999,) + base, [1, 2, 3])):
        grown = _columns("cylinder_solid", added)
        for name, array in honest.items():
            np.testing.assert_array_equal(grown[name][keep], array, err_msg=f"{name} after adding {added}")


# ── K3: the stretch, carried by the rule's speed and by the kernel's 3-D chord ──


@pytest.mark.foundation
@pytest.mark.parametrize("sine", _SINES)
@pytest.mark.parametrize(("body", "b", "expected"), [r[1:] for r in _K1_LINES], ids=[r[0] for r in _K1_LINES])
@pytest.mark.rests_on(_HERE + "test_each_speed_of_the_family_is_its_3d_lines_outflow_and_entry_response")
def test_the_in_plane_line_at_speed_kappa_is_the_3d_chord_at_speed_one(body, b, expected, sine) -> None:
    """[K3; foundation, THEOREM] The in-plane line at speed 1/sin(theta) (the stretch carried by the rule's speed)
    and the 3-D line (b, theta) chorded by the geometry kernel at speed 1 (the stretch carried by the chord: the path
    ``angular_flux`` keeps) give the same outflow, entry response and closed angular flux, 1e-12 of the largest
    entry. With one member, the family's grading in 3-D length is the 3-D chord's (depths 2^k / Sigma, 7 doublings),
    so the two differ by the chords' rounding only (`[M]` ``ta/chord_lift_probe.log``: the kernel's 3-D slot lengths
    are the in-plane ones over sin(theta) to <= 418 ulp). `[R]` the band; pin 10x the first green's reading.

    First reds: the speed dropped from the length weight or from the depth; the speed squared in the depth.
    """
    family = _Family(body, b, (sine,), T._SIGMA[0])
    case = family.case(0)
    plane, chord = family.rule, case.rule
    for name in ("outflow", "entry_response"):
        got, want = getattr(plane, name)()[0], getattr(chord, name)()
        scale = np.max(np.abs(want))
        assert np.max(np.abs(got - want)) <= 1e-12 * scale, (name, np.max(np.abs(got - want)) / scale)
    for transit, rev, _end in expected:
        if rev:
            continue
        # each rule reads at its own chord's crossings (the same levels, in the same order: t3 = t2 / sin(theta))
        t2, t3 = _crossings_on_transit(plane, transit), _crossings_on_transit(chord, transit)
        np.testing.assert_allclose(t3 * sine, t2, rtol=1e-12, atol=1e-12)
        got = plane.angular_flux(t2, plane.period.inflow(plane.optical_depth, plane.outflow()))[0]
        want = chord.angular_flux(t3, chord.period.inflow(chord.optical_depth, chord.outflow()))
        scale = np.max(np.abs(want))
        assert np.max(np.abs(got - want)) <= 1e-12 * scale, ("angular_flux", np.max(np.abs(got - want)) / scale)


# ── K4: the polar fold of one in-plane line, against Bickley ─────────────

#: A homogeneous cylinder behind a partial mirror, and one in-plane line: in-plane optical width tau = 2 Sigma
#: sqrt(R^2 - b^2) (2.383 here).
_FOLD_SIGMA, _FOLD_R, _FOLD_B, _FOLD_A = 1.3, 1.0, 0.4, 0.6


def _fold_errors(points: int) -> tuple[float, float, float]:
    """The three polar folds of one in-plane line on the polar rule at ``points``, relative to their closed forms."""
    from orpheus.geometry.boundary import AlbedoBoundary, SpecularReturn
    from orpheus.geometry.structured_geometry import StructuredGeometry

    geometry = StructuredGeometry.cylinder((0.0, _FOLD_R), (0,), outer=AlbedoBoundary(_FOLD_A, SpecularReturn(axis="x")))
    basis = PanelBasis.of(ConcentricPartition.of(geometry), *T._RES)
    polar = polar_rule(OpticalScale.of(np.array([0.0, _FOLD_R]), np.array([_FOLD_SIGMA])), points)
    sine = np.sin(polar.pts)
    line = Line.through(np.array([_FOLD_B, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))
    rule = TraversalRule.of(line, basis, Walls.of(geometry), np.array([_FOLD_SIGMA]), speed=1.0 / sine)
    one = np.ones(basis.size)                                           # q = 1: the panels reproduce a constant
    B = rule.outflow()[:, 0] @ one                                      # (V,)
    reading = float(np.sum(polar.wts * sine * B))
    block = float(one @ rule.volterra(polar.wts * 2.0 * sine**2) @ one)
    inflow = rule.period.inflow(rule.optical_depth, rule.outflow())[:, 0] @ one
    closure = float(np.sum(polar.wts * sine * np.exp(-rule.optical_depth[:, 0]) * inflow))
    want_reading, want_block, want_closure = _fold_references(_FOLD_SIGMA, _FOLD_R, _FOLD_B, _FOLD_A)
    return (abs(reading / want_reading - 1.0), abs(block / want_block - 1.0), abs(closure / want_closure - 1.0))


#: The absolute target of the closure's Neumann tail, far below a double's resolution of the sum (~0.05).
_NEUMANN_TARGET = 1e-25


@functools.cache
def _fold_references(sigma: float, radius: float, b: float, albedo: float) -> tuple[float, float, float]:
    """K4's three closed forms at 30 digits, once per fixture: (1 - Ki_2(tau))/Sigma, 2(tau - pi/4 + Ki_3(tau))/Sigma^2,
    and the closure's Neumann series sum_m a^{m+1}(Ki_2((m+1)tau) - Ki_2((m+2)tau))/Sigma in explicit terms.

    The term count: Ki_2(x) = int_0^{pi/2} cos t e^{-x/cos t} dt <= e^{-x}, so the tail from term M is at most
    sum_{m >= M} (a e^{-tau})^{m+1} = q^{M+1}/(1 - q), q = a e^{-tau}; M = ceil(ln(target (1 - q))/ln q) - 1 terms past
    the first reach ``_NEUMANN_TARGET``: 20 terms for a = 0.6, tau = 2.383 (q = 0.0554). ``mp.nsum`` cost 520 s here
    (`[M]` 2026-10-12, ``ta/confirm2/k4_tot.txt``: 99 % in mpmath's hypergeometric sums at growing arguments)."""
    with mp.workdps(30):
        S, a = mp.mpf(repr(sigma)), mp.mpf(repr(albedo))
        tau = 2 * S * mp.sqrt((mp.mpf(repr(radius)) - mp.mpf(repr(b))) * (mp.mpf(repr(radius)) + mp.mpf(repr(b))))
        q = a * mp.exp(-tau)
        terms = int(mp.ceil(mp.log(mp.mpf(_NEUMANN_TARGET) * (1 - q)) / mp.log(q)))
        reading = (1 - mpref.ki2(tau)) / S
        block = 2 / S**2 * (tau - mp.pi / 4 + R.bickley_ki3(tau))
        closure = mp.fsum(a ** (m + 1) * (mpref.ki2((m + 1) * tau) - mpref.ki2((m + 2) * tau)) for m in range(terms)) / S
        return float(reading), float(block), float(closure)


@pytest.mark.l2
@pytest.mark.verifies("characteristic-quadrature")
@pytest.mark.rests_on(_HERE + "test_each_speed_of_the_family_reads_its_3d_lines_closed_angular_flux",
                      _HERE + "test_the_familys_volterra_triangle_is_the_weighted_sum_of_its_3d_lines",
                      _ASSEMBLY + "test_the_cylinders_direction_rule_is_gauss_in_the_polar_angle")
def test_the_polar_fold_of_one_line_converges_to_bickleys_functions() -> None:
    """[K4; l2, CONV] Summed over the polar rule at 4, 6, 8 points, one in-plane line's family folds to Bickley's
    functions (the polar integral done exactly): the reading's outflow sum_j w_j sin_j B_j to (1 - Ki_2(tau)) / Sigma;
    the block's sum_j w_j 2 sin_j^2 <1, V_j 1> to 2 (tau - pi/4 + Ki_3(tau)) / Sigma^2; the exit inflow under the
    mirror a, sum_j w_j sin_j e^{-kappa_j tau} in_j, to sum_m a^{m+1} (Ki_2((m+1) tau) - Ki_2((m+2) tau)) / Sigma
    (the Neumann series, which converges for a < 1: a test reference, refuted only as production's form). The errors
    contract strictly above 1e-14 and are below 1e-13 at 8 points (10 s, `[M]` 2026-10-12, once the references were cached
    and the Neumann series summed in 20 explicit terms, :func:`_fold_references`; with ``mp.nsum`` it took 523 s).

    `[M]` 2026-10-11 on today's per-(b, theta) path (``ta/k4_formulas.log``): (1.7e-8, 7.8e-9, 3.8e-7) at 4 points,
    (1.8e-15, 1.1e-15, 2.7e-14) at 8, (0, 2.2e-16, 4.4e-16) at 16. Ki_2 and Ki_3 are Struve-based closed forms, no
    line integral of the code. First reds: the polar weight's sin dropped (`[M]` prototype +5.7e-1 on a wall
    reading); the block's speed applied once (the fold reads Ki_2's family); the closure with a dropped from its
    denominator.
    """
    errors = np.array([_fold_errors(n) for n in (4, 6, 8)])           # (rung, fold)
    for fold in range(3):
        column = errors[:, fold]
        for before, after in zip(column, column[1:]):
            assert after < before or after < 1e-14, ("not contracting", fold, column)
        assert column[-1] < 1e-13, ("at 8 points", fold, column)


# ── G1: the dyadic grading of the rules over lines (along lines the attenuation is integrated by its Gauss rule) ──
#
# P3 step 3f (route B) retired the along-line exponential grading: each piece's attenuation is the measure of its own
# Gauss rule (``attenuation_rule``, gated in ``test_characteristic_attenuation_rule.py``), and the derived law
# (``ExponentialGrading.of``) and its cap went with it, with their rows (G1, the cap, G1b). The end law survives for the
# rules over lines, whose coordinates (the polar angle, the cosine, the chord half-length) do not carry the attenuation
# linearly: ``DYADIC``, first depth 1 and ratio 2.


def _hand_ends(stop: float, start: float, steepest: float, gentlest: float, first: float, ratio: float) -> np.ndarray:
    """The end law by hand, given (first, ratio): depths first ratio^k / steepest from ``stop``,
    k < K = ceil(log_ratio(VANISHING_DEPTH steepest / (first gentlest))) + 1 (50-digit logarithm, so K is exact at a
    power of the ratio), clipped to the interval; none where steepest x width <= first."""
    width = abs(stop - start)
    ends = {start, stop}
    if steepest * width > first:
        with mp.workdps(50):
            count = int(mp.ceil(mp.log(mp.mpf(VANISHING_DEPTH) * mp.mpf(steepest) / (mp.mpf(first) * mp.mpf(gentlest)))
                                / mp.log(mp.mpf(ratio)))) + 1
        inward = math.copysign(1.0, start - stop)
        ends |= {stop + inward * min((first / steepest) * ratio**k, width) for k in range(count)}
    return np.array(sorted(ends))


_G1_ENDS = [  # (id, stop, start, steepest, gentlest, grading)
    ("dyadic-one-speed", 0.0, 3.0, 5.0, 5.0, "dyadic"),                 # 7 ends: log_2 64 = 6 exactly
    ("dyadic-reversed", 3.0, 0.0, 5.0, 5.0, "dyadic"),
    ("dyadic-kappa-1e4", 0.0, 3.0, 5.0e4, 5.0, "dyadic"),               # the reach of the gentlest member
    ("dyadic-clipped", 0.0, 0.01, 5.0e4, 5.0, "dyadic"),                # the slot ends before the reach
    ("dyadic-thin", 0.0, 0.1, 9.0, 5.0, "dyadic"),                      # optical width 0.9 <= first 1: not graded
]


@pytest.mark.foundation
@pytest.mark.parametrize(("stop", "start", "steepest", "gentlest", "which"), [r[1:] for r in _G1_ENDS],
                         ids=[r[0] for r in _G1_ENDS])
def test_the_exponential_ends_follow_the_gradings_first_depth_and_ratio(stop, start, steepest, gentlest, which) -> None:
    """[G1, the end law; foundation, THEOREM] Given a grading (first h, ratio rho), ``exponential_ends`` cuts at
    h rho^k / Sigma_steep from ``stop`` until VANISHING_DEPTH mean free paths of the GENTLEST member, clipped, none on
    a stretch of optical width at most h (:func:`_hand_ends`): the first depth is the steepest (slowest) member's, the
    reach the gentlest (fastest) member's. Rounding of rho^k: 4 eps. First reds (`[M]` on the prototype,
    ``ta/floor_probe_8.log``, ``ta/reach_probe.log``): the reach at the steepest member (7.8e-9 to 2.7e-3 on the
    closed-form readings); the first depth from the gentlest; the thin threshold at 2 instead of the grading's first."""
    grading = DYADIC
    got = np.unique(exponential_ends(np.array(stop), np.array(start), np.array(steepest), np.array(gentlest), grading))
    want = _hand_ends(stop, start, steepest, gentlest, grading.first, grading.ratio)
    assert got.shape == want.shape, (got, want)
    np.testing.assert_allclose(got, want, rtol=4 * np.finfo(float).eps, atol=4 * np.finfo(float).eps * abs(stop - start))


@pytest.mark.foundation
def test_every_rule_over_lines_is_graded_dyadically(monkeypatch: pytest.MonkeyPatch) -> None:
    """[G1c; foundation, ROUTE] ``DYADIC``'s callers, every rule over lines: the polar and cosine rules toward the
    normal direction, and the impact rule toward its rim (the user's ruling of 2026-10-12; since P3 step 3f, along
    lines the attenuation is integrated by its Gauss rule and no line is graded toward a mean free path; the line
    rules' nodes are in the polar angle, the cosine or the chord half-length, where the attenuation's exponent is not
    linear). ``DYADIC`` is (1, 2). Recorded on the ``lines`` module's binding
    while each rule is built at 8 points. First reds, on the grading leg (`[M]` 2026-10-12, ``ta/e/runs/``, the arms exec'd
    in the live module namespace, ``TA_LIVE=1``, so the spy sees the mutant's calls): the polar rule graded by the node
    count (``dyadic_by_count``), the impact rule graded by it (``impact_by_count``), each recording
    ``ExponentialGrading(1.0, 1.419)``."""
    from orpheus.derivations.continuous.characteristic import lines as lines_module
    from orpheus.derivations.continuous.characteristic.lines import cosine_rule, impact_rule

    seen: list[object] = []
    honest = lines_module.exponential_ends

    def recording(stop, start, steepest, gentlest, grading):
        seen.append(grading)
        return honest(stop, start, steepest, gentlest, grading)

    monkeypatch.setattr(lines_module, "exponential_ends", recording)
    scale = OpticalScale.of(np.array([0.0, 0.5, 1.5, 2.0]), np.array([0.6, 1.3, 0.45]))
    polar_rule(scale, 8)
    cosine_rule(scale, 8)
    assert seen and all(g is DYADIC for g in seen), seen
    assert DYADIC == ExponentialGrading(1.0, 2.0)                       # its value, moved from the retired cap row
    seen.clear()
    impact_rule(np.array([0.0, 0.5, 1.5, 2.0]), np.array([0.6, 1.3, 0.45]), 8, np.array([0.1, 0.1, 0.1]))
    assert seen and all(g is DYADIC for g in seen), seen


# ── F and G2: readings against frozen mpmath routes ──────────────────────

#: (id, breakpoints, outer specular amplitude, x, Sigma per region, q per region, the frozen route, band). Each literal
#: is ``_characteristic_point_mp``'s route (``phi_cylinder`` under a mirror, ``phi_cylinder_vacuum`` under vacuum),
#: printed by ``scratch/characteristic_architecture/p3/step3/ta/frozen_refs.py`` and ``floor_probe.py``
#: (2026-10-11) and recomputed by :func:`test_the_frozen_cylinder_routes_are_the_mpmath_routes` (slow). Bands: 10 x
#: the prototype's error with its closure corrected (``ta/floor_probe_8.log``, ``ta/level_probe.log``).
_F_ROWS = [
    ("F1-mirror0.9-wall", (0.0, 2.0), 0.9, 2.0, (2.4,), (0.7,), 0.2733761078614076, 1e-13),
    ("F1b-mirror0.999-wall", (0.0, 2.0), 0.999, 2.0, (2.4,), (0.7,), 0.2914156825428953, 1e-13),
    ("F4-mr3-vacuum-1.1", (0.0, 0.5, 1.5, 2.0), 0.0, 1.1, (0.6, 1.3, 0.45), (1.0, 0.7, 0.4), 0.49269148377062716, 1e-13),
]
#: F4 reads a three-region cylinder's block: 761 s before the carve (`[M]` 2026-10-11, ``ta/f_rows_today.log``), so it
#: is slow until the lifted block's cost is measured; every other F row took 14 to 55 s before the carve.
_F_SLOW = {"F4-mr3-vacuum-1.1"}
_G2_ROWS = [
    ("F2-thick-0.5", (0.0, 1.0), 0.0, 0.5, (50.0,), (1.0,), 0.019999999999992617, 3e-12),
    ("F2-thick-wall", (0.0, 1.0), 0.0, 1.0, (50.0,), (1.0,), 0.00994999624718132, 3e-12),
    ("F3-thin-wall", (0.0, 1.0), 0.0, 1.0, (0.01,), (1.0,), 0.9713941339198011, 5e-11),
]


def _reading(breakpoints, amplitude, x, sigma, q) -> float:
    return float(RE._given("cylinder", breakpoints, ((amplitude, 0.0),), x, RE._WORK, (tuple(sigma),), (tuple(q),))[0])


@pytest.mark.l1
@pytest.mark.verifies("characteristic-quadrature")
@pytest.mark.parametrize(("breakpoints", "amplitude", "x", "sigma", "q", "route", "band"),
                         [pytest.param(*r[1:], id=r[0], marks=pytest.mark.slow if r[0] in _F_SLOW else ()) for r in _F_ROWS])
@pytest.mark.rests_on(_HERE + "test_each_speed_of_the_family_reads_its_3d_lines_closed_angular_flux",
                      _READING + "test_the_point_rule_integrates_the_direction_moments")
def test_a_cylinders_reading_is_the_frozen_mpmath_route(breakpoints, amplitude, x, sigma, q, route, band) -> None:
    """[F1, F1b, F4; l1, REFERENCE] The lifted point reading of a cylinder (degree 3, 2 layers, 8 line points, 12 along
    each line) equals the frozen mpmath route, in the fast tier: the partial mirror on the wall at a = 0.9 (F1) and
    0.999 (F1b, ERR-106's regime: the in-plane lines chorded from their levels), and the MR3 body under vacuum
    inside (F4).

    First reds (`[M]` on the prototype, ``ta/proto_battery.log``, ``ta/level_probe.log``): the closure without a in
    its denominator (F1 +1.1e-2); the closure's kappa from one member (-2.7e-2); B without the speed (-1.1e-1); the
    impact rule at speed 1 (F1 +3.2e-9); the in-plane lines chorded from b, not their levels (F1b +5.5e-12, not
    converging); kappa inverted, V without the speed, the polar sin dropped (F4 O(1)). Declared blind: F4 to every
    closure arm (vacuum) and to the impact rule at speed 1 (`[M]` -4.1e-15: x = 1.1's lines miss the slow layer).
    """
    got = _reading(breakpoints, amplitude, x, sigma, q)
    assert abs(got / route - 1.0) <= band, (got, route, got / route - 1.0)


@pytest.mark.l1
@pytest.mark.verifies("characteristic-quadrature")
@pytest.mark.parametrize(("breakpoints", "amplitude", "x", "sigma", "q", "route", "band"), [r[1:] for r in _G2_ROWS],
                         ids=[r[0] for r in _G2_ROWS])
@pytest.mark.rests_on("tests/gates/derivations/test_characteristic_attenuation_rule.py::"
                      "test_the_rule_is_exact_on_every_legendre_polynomial_below_degree_two_m")
def test_the_along_line_grading_holds_a_thick_and_a_thin_cylinder_to_their_routes(breakpoints, amplitude, x, sigma, q,
                                                                                  route, band) -> None:
    """[G2 = F2, F3; l1, REFERENCE] The value leg of G1: a thick cylinder (Sigma R = 50) inside and on its wall, a thin
    one (Sigma R = 0.01) on its wall, under vacuum, against the frozen ``phi_cylinder_vacuum`` routes.

    `[M]` on the prototype (``ta/reach_probe.log``, ``ta/floor_probe_8.log``): the reach at 7 doublings -1.9e-3,
    -2.7e-3 (thick), -7.8e-9 (thin); the first depth at 16 x the slowest member's mean free path -1.6e-11, -6.7e-11
    (thick); at 64 x, -3.2e-7, -1.3e-6. The thin band is the polar rule's own floor in that regime (the prototype read
    +5.3e-12). First reds: those three gradings.
    """
    got = _reading(breakpoints, amplitude, x, sigma, q)
    assert abs(got / route - 1.0) <= band, (got, route, got / route - 1.0)


@pytest.mark.foundation
@pytest.mark.slow
@pytest.mark.parametrize(("breakpoints", "amplitude", "x", "sigma", "q", "route"), [r[1:7] for r in _F_ROWS + _G2_ROWS],
                         ids=[r[0] for r in _F_ROWS + _G2_ROWS])
def test_the_frozen_cylinder_routes_are_the_mpmath_routes(breakpoints, amplitude, x, sigma, q, route) -> None:
    """[F slow twin] Each frozen literal of F1-F4 is its mpmath route, to 1e-15 relative: the literal is a reference
    computed from the law, never from the code under test. The routes are stored by the traced memo, so a warm store
    serves them; a cold one computes them, 1 s to 4 min each (`[M]` 2026-10-11, the thick wall 239 s). First red: a
    literal edited to a reading of the code."""
    if amplitude == 0.0:
        want = mpref.phi_cylinder_vacuum(tuple(breakpoints), tuple(sigma), tuple(q), x)
    else:
        want = mpref.phi_cylinder(tuple(breakpoints), tuple(sigma), tuple(q), amplitude, x)
    assert abs(want / route - 1.0) <= 1e-15, (want, route)


# ── P1: #587's retirement loses no accuracy, per polar node ──────────────

_SMALL = (0.0, 0.5, 1.0)
_SMALL_SIGMA = (0.6, 2.4)
_SMALL_Q = (1.0, 0.7)


def _polar_column_errors(x: float, albedo: float) -> np.ndarray:
    """Per polar node j of the point rule at ``x`` (group 0 of Q4's small cylinder), the reading's theta-column (the
    member's lines, weighted, both sides) relative to w_j g(theta_j), g the mpmath per-theta density."""
    system = RE._absorbing("cylinder", _SMALL, ((albedo, 0.0),), (2, 0, 8), (_SMALL_SIGMA, (1.7, 0.8)))
    q = RE._emission(system, (_SMALL_Q, (0.4, 1.3)))
    group = system.groups[0]
    basis, walls = system.basis, system.walls
    sigma = np.asarray(_SMALL_SIGMA)
    coef = system.emission.split(q)[0]
    rule = PointRule.of(basis, walls, sigma, x, 8)
    lines = rule.lines
    column = np.zeros(lines.speeds.size)
    for traversals, chunk in lines.chunks():
        t = traversals.period.chord.image.parameters_at(np.array([x]), rule.side)
        stacked = reading_module._stacked_flux(traversals, t, group)                     # (V, l, q, M + W)
        psi = stacked[..., : group.support.size] @ coef
        column += np.einsum("vl,vlq->v", lines.weights[:, chunk], psi)
    ends = np.asarray(basis.partition.breakpoints)
    split = np.unique(np.append(ends, x))
    polar = polar_rule(OpticalScale.of(split, basis.on_panels(sigma)[np.searchsorted(ends, split[:-1], side="right") - 1]), 8)
    theta = np.arcsin(1.0 / lines.speeds)
    order = np.argsort(theta)
    np.testing.assert_allclose(theta[order], polar.pts, rtol=1e-13, atol=0.0)
    density = np.array(mpref.cylinder_polar_density(_SMALL, _SMALL_SIGMA, _SMALL_Q, albedo, x, tuple(map(float, polar.pts))))
    return column[order] / (polar.wts * density) - 1.0


@pytest.mark.l1
@pytest.mark.verifies("characteristic-quadrature")
@pytest.mark.parametrize(("x", "albedo"), [(1.0, 0.0), (0.5, 0.0), pytest.param(1.0, 0.9, marks=pytest.mark.slow)],
                         ids=["wall", "interface", "wall-mirror"])
@pytest.mark.rests_on(_ASSEMBLY + "test_the_impact_rule_is_graded_at_the_slowest_polar_node",
                      _READING + "test_a_small_cylinders_layer_is_graded_outside_slow")
def test_each_polar_node_of_the_reading_is_its_mpmath_density(x: float, albedo: float) -> None:
    """[P1; l1, REFERENCE] #587's retirement loses no accuracy: per polar node theta_j of the point rule on Q4's small
    cylinder (degree 2, no layers, 8 points), the reading's theta-column is w_j g(theta_j) to 1e-13, g(theta) =
    (1/pi) sin(theta) int_0^pi psi_disk(x, cos alpha, stretch 1/sin theta) dalpha in mpmath, on every node.

    `[M]` 2026-10-11 on today's tree (``ta/per_polar_probe.log``), 128 nodes each: #587's own-speed impact rule
    3.7e-14, 4.4e-14, 3.6e-14 (wall, wall-mirror, interface); every node graded at the slowest speed (the carve)
    3.7e-14, 4.4e-14, 3.5e-14, worse than the own-speed rule on 0 of 384 nodes; graded at speed 1, 2.9e-4, 2.1e-5,
    8.0e-4 per node while the TOTAL moved only 8.5e-9, 2.9e-9, 1.1e-8 (a signed cancellation over the nodes,
    ``vv-principles`` anti-#24 (d)): the per-node functional is the instrument. First red: the impact rule at speed 1.
    The a = 0.9 fixture's reference takes 20 s: slow.
    """
    errors = _polar_column_errors(x, albedo)
    assert np.max(np.abs(errors)) <= 1e-13, (np.max(np.abs(errors)), int(np.argmax(np.abs(errors))))


# ── S1: the sphere and the slab, moved by the target-graded reading, accepted on the three conditions ──

#: The pre-carve capture: ``main`` at ``4ebfa7a4``, by ``scratch/characteristic_architecture/p3/step3/ta/
#: capture_speed1.py`` and ``capture_totals.py`` run against a pristine copy of ``orpheus/`` (2026-10-11; their
#: logs). Blocks: 1^T K 1 (the carried scalar the drift is measured on). Readings: both groups.
_S1_BLOCKS = [
    ("sphere_mirror", "sphere", (0.0, 0.5, 1.5, 2.0), ((1.0, 0.0),), "g0", 8, 48.07776185525604),
    ("hsphere_white_mirror", "sphere", (0.4, 0.5, 1.5, 2.0), ((0.0, 1.0), (1.0, 0.0)), "g1", 8, 27.494611277752146),
    ("slab_mirror_white", "slab", (0.0, 0.4, 1.5, 2.3), ((1.0, 0.0), (0.0, 1.0)), "g0", 4, 2.8386528063659466),
    ("slab_white_white", "slab", (0.0, 0.4, 1.5, 2.3), ((0.0, 1.0), (0.0, 1.0)), "g1", 4, 2.4145876269732365),
]
_S1_READINGS = [
    ("sphere-1.1", "sphere", (0.0, 0.5, 1.5, 2.0), ((0.6, 0.0),), 1.1, (0.825384546564032, 2.021243639703333)),
    ("sphere-wall", "sphere", (0.0, 0.5, 1.5, 2.0), ((0.6, 0.0),), 2.0, (2.717370966981192, 0.3616090157370513)),
    ("slab-0.7", "slab", (0.0, 0.4, 1.5, 2.3), ((0.3, 0.0), (0.8, 0.0)), 0.7, (0.5703883052832583, 2.4474738857710125)),
    ("slab-right", "slab", (0.0, 0.4, 1.5, 2.3), ((0.3, 0.0), (0.8, 0.0)), 2.3, (3.580680913511321, 0.341944374583829)),
]
#: The bars: `[M]` 2026-10-12 at these rows' working points (``ta/s1_drift.log``) the drift of 1^T K 1 is <= 2.2e-16 on
#: all four blocks (conservation 1.3e-13 to 1.8e-15) and of the readings <= 2.2e-16 (route <= 1.2e-14); the block band
#: is ~ 50 x that measurement, the readings' the reading file's C1/C2 bar (``_GIVEN_BAND``, 1e-12).
_S1_BLOCK_DRIFT = 1e-14


@pytest.mark.l1
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.parametrize(("chart", "breakpoints", "laws", "group", "points", "captured"), [r[1:] for r in _S1_BLOCKS],
                         ids=[r[0] for r in _S1_BLOCKS])
@pytest.mark.rests_on(_ASSEMBLY + "test_a_closed_body_conserves_its_emission")
def test_a_sphere_or_slab_block_moved_by_the_carve_conserves_and_stays_within_its_accuracy(chart, breakpoints, laws,
                                                                                          group, points, captured) -> None:
    """[S1, blocks; l1, the three-condition acceptance of a quadrature change] The target-graded reading re-graded the
    pieces the block's response is integrated on (exit- and entry-graded one-level outflows), so a sphere's and a
    slab's block moved (`[M]` the coordinator, 2026-10-11, at 4 and 6 points: <= 3e-13; at these rows' points <= 2.2e-16)
    and is not bit-identical. Accepted because: (1) every moved number is a named quantity, the block of the same
    rule on re-graded pieces; (2) against a structurally independent reference, each of these closed bodies still
    conserves (K sigma = W 1, the volume measure only) to its AS1 band; (3) the drift of 1^T K 1 from the pre-carve
    capture is below ``_S1_BLOCK_DRIFT``, a quadrature-sized move and not an algorithmic one. First reds: a re-grading
    that loses accuracy (the conservation leg), any change moving the block by more than its accuracy (the drift
    leg: e.g. the response's outflow integrated toward the wrong end)."""
    conservation = AS._conservation(chart, breakpoints, laws, AS._SIGMA[group], points)
    block = AS._transport(chart, breakpoints, laws, AS._SIGMA[group], points).block
    assert conservation <= AS._CONSERVATION_TOL[chart], (conservation, AS._CONSERVATION_TOL[chart])
    drift = abs(float(block.sum()) / captured - 1.0)
    assert drift <= _S1_BLOCK_DRIFT, drift


@pytest.mark.l1
@pytest.mark.verifies("characteristic-quadrature")
@pytest.mark.parametrize(("chart", "breakpoints", "laws", "x", "captured"), [r[1:] for r in _S1_READINGS],
                         ids=[r[0] for r in _S1_READINGS])
@pytest.mark.rests_on(_READING + "test_the_reading_at_a_general_point_is_the_mpmath_route")
def test_a_sphere_or_slab_reading_moved_by_the_carve_stays_on_its_mpmath_route(chart, breakpoints, laws, x,
                                                                               captured) -> None:
    """[S1, readings; l1, the three-condition acceptance of a quadrature change] The target-graded reading (one-level
    sums of each upstream piece's integral on the chord refined at the point's level) moved a sphere's and a slab's
    point readings (`[M]` the coordinator, 2026-10-11: ~5e-11 relative at 4 line points and 6 along the line). Accepted
    because: (1) the reading is the same named quantity on a re-graded rule; (2) against the structurally independent
    mpmath route (``_characteristic_point_mp``: the sphere's disk billiard, the slab's E_2 image series) both the new
    reading and the pre-carve capture are within the C1/C2 bar ``_GIVEN_BAND``; (3) so the drift between them is
    bounded by twice that bar, a quadrature move. Measured here at the reading file's working point (8 line points,
    12 along). First reds: a reading that leaves the route (the new value outside the bar); a capture edited to
    the new value is caught by leg (2) only if it leaves the route, so the capture is printed by its script."""
    got = RE._given(chart, breakpoints, laws, x)
    for g in range(len(got)):
        route = RE._route(chart, breakpoints, laws, x, g)
        assert abs(captured[g] / route - 1.0) < RE._GIVEN_BAND, ("the capture", g, captured[g], route)
        assert abs(float(got[g]) / route - 1.0) < RE._GIVEN_BAND, ("the reading", g, float(got[g]), route)
        assert abs(float(got[g]) / captured[g] - 1.0) < 2.0 * RE._GIVEN_BAND, ("the drift", g, float(got[g]), captured[g])
