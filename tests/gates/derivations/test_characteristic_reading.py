"""Gates for the characteristic reference's reading at a point (:mod:`~orpheus.derivations.continuous.characteristic.reading`).

P1 step (b), rung 5b, of the characteristic-reference campaign
(``.claude/plans/characteristic_reference_architecture.md``, "P1 step (b), rung 5b: API sketch", ruled
2026-10-08): ``PointValue`` reads the transported emission ``(K q)(x)`` (iterated Galerkin), on the line rule's own
lines through the point, read at the point's parameters, weighted by ``dOmega / 4 pi``; the diffuse walls enter
through ``WallCoupling.currents``; ``psi(x, Omega)`` is exposed for gates. Verification spec
``scratch/characteristic_architecture/p1_step_b5b/spec.md`` (rows M1-M3, C1, C2, C2b, C4, C6, C7, C8, W1, W2, F1-F5,
E1, E4, D1-D7, and the re-posed D10 of the door file); P1 spec ``scratch/characteristic_architecture/p1_verification_spec.md`` §5.

The ladder, bottom up (``rests_on`` on each row):

1. rungs 1-4 and 5a, landed: the walls, the closure, the basis, the line transport, the assembly
   (``test_characteristic_assembly.py``), the system, the door;
2. the point rule's measure [M1-M3]: the direction moments, the point's parameters, one constructor;
3. the reading of a GIVEN emission [C1, C2, C2b, C4, C7]: against mpmath routes that share no line, no rule and no
   closure with the code (``_characteristic_point_mp.py``);
4. the reading against the block [C6, C8, W1, W2]: one transport, two test measures; the currents refactor;
5. the angular flux [F1-F5];
6. the door [E1, E4, D1-D7]: every question's ``PointValue``, its refusal on ``Nearest``, a ``Ratio``.

Every reference value is written in mpmath or by hand from the ``Mixture`` arrays and the geometry's numbers
(``instrument-doctrine`` X4). Bands are measured (`[M]` 2026-10-08, ``.venv/bin/python -O``,
``scratch/characteristic_architecture/p1_step_b5b/ta/``); each row's first red is a battery arm
(``scratch/characteristic_architecture/p1_step_b5b/ta/battery/``). Nothing here depends on the reading's chunking.
"""
from __future__ import annotations

from collections.abc import Callable
from functools import lru_cache
from typing import Any

import numpy as np
import pytest

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.data.materials import Materials
from orpheus.derivations.continuous.characteristic import (
    CharacteristicDerivation, GalerkinSystem, LineRule, PanelBasis, RegionCrossSections, Resolution, TransportResolution,
    characteristic_reference,
)
from orpheus.derivations.continuous.characteristic import reference as reference_module
from orpheus.derivations.continuous.characteristic.lines import Lines
from orpheus.derivations.continuous.characteristic.reading import PointRule
from orpheus.geometry import BC, StructuredGeometry
from orpheus.geometry.boundary import AlbedoBoundary, SpecularReturn
from orpheus.geometry.chart import RadialImage
from orpheus.numerics.mesh_free_function import RegionwiseConstant
from orpheus.numerics.observable import FluxIntegral, PointValue, Ratio
from orpheus.numerics.question import Eigen, FixedSource, Nearest, Response
from orpheus.numerics.traced_memo import bypass
from orpheus.reference.reading import Uncertified
from orpheus.specification.specification import GeometrySpecification
from tests.gates.derivations import _characteristic_point_mp as mpref
from tests.gates.derivations._characteristic_mp import specular_sphere_total
from tests.gates.derivations.test_characteristic_assembly import assert_each_polar_angle_is_graded_at_its_own_speed
from tests.gates.derivations.test_characteristic_system import (
    _ABS, _MIRROR, _PU2, _UP2N, _VACUUM, _WHITE, _basis, _mixture, _walls, _zero_d,
)
from tests.gates.derivations.test_peierls_greens_function_garcia2021 import (
    GARCIA_2021_CASE1_PHI, GARCIA_2021_CASE1_R, GARCIA_2021_CASE1_ROUNDING,
)

pytestmark = pytest.mark.filterwarnings("error::RuntimeWarning")

_HERE = "tests/gates/derivations/test_characteristic_reading.py::"
_ASSEMBLY = "tests/gates/derivations/test_characteristic_assembly.py::"
_SYSTEM = "tests/gates/derivations/test_characteristic_system.py::"
_REFERENCE = "tests/gates/derivations/test_characteristic_reference.py::"
_CONSERVES = _ASSEMBLY + "test_a_closed_body_conserves_its_emission"
_ESCAPE = _ASSEMBLY + "test_the_escape_and_transmission_probabilities_are_the_closed_forms"
_E2 = _SYSTEM + "test_a_closed_body_with_a_uniform_source_reads_the_infinite_medium_flux"
_D5 = _SYSTEM + "test_a_closed_body_reads_k_inf_with_a_flat_flux_in_the_infinite_mediums_group_ratio"
_PAIRING = _REFERENCE + "test_a_flux_integral_is_the_volume_integral_of_the_galerkin_flux"
_M1 = _HERE + "test_the_point_rule_integrates_the_direction_moments"
_C2 = _HERE + "test_the_reading_at_a_general_point_is_the_mpmath_route"

_MR3 = (0.0, 0.5, 1.5, 2.0)
_MR3H = (0.4, 0.5, 1.5, 2.0)
_SLB3 = (0.0, 0.4, 1.5, 2.3)
#: The P1 spec's MR3 data (§4): Sigma_t and the emission RATE q per region, per group; q and Sigma_t jump at both
#: interfaces, distinct per group.
_SIGMA = ((0.6, 1.3, 0.45), (1.7, 0.35, 2.4))
_Q = ((1.0, 0.25, 3.0), (0.4, 2.0, 0.7))
#: The working point (rung 4's ``_WORK``): degree 3, 2 layers of ratio 0.4, 8 line points, 12 along each line.
_WORK = (3, 2, 8, 12)
_GIVEN_BAND = 1e-12                 # the P1 spec's C1/C2 bar; [M] the worst sphere/slab/hollow row reads 6.7e-14


def _partial(albedo: float) -> AlbedoBoundary:
    return AlbedoBoundary(albedo, SpecularReturn(axis="x"))


# ── a given emission: pure absorbers with a source posed in every region ──


@lru_cache(maxsize=None)
def _absorbing(chart: str, breakpoints, laws, resolution=_WORK, sigma=_SIGMA) -> GalerkinSystem:
    """The system of pure absorbers of total cross sections ``sigma`` (per group, per region), a source posed everywhere.

    With no scattering and no fission the emission is the given q: the reading is ``(K q)(x)`` of that q alone.
    """
    degree, layers, line_points, points = resolution
    absorbers = [_mixture([s[j] for s in sigma], np.zeros((len(sigma), len(sigma)))) for j in range(len(sigma[0]))]
    return GalerkinSystem(
        _basis(chart, breakpoints, degree, layers, 0.4), _walls(chart, breakpoints, laws),
        RegionCrossSections.of(absorbers), TransportResolution(line_points, points, points),
        source_regions=np.ones((len(sigma[0]), len(sigma)), dtype=bool),
    )


def _emission(system: GalerkinSystem, q) -> np.ndarray:
    """A per-region table q ``(G, n)`` on the emission space, read onto the nodes exactly."""
    return system.emission.restrict(system.basis.on_nodes(np.asarray(q, dtype=float).T).T)


def _given(chart: str, breakpoints, laws, x: float, resolution=_WORK, sigma=_SIGMA, q=_Q) -> np.ndarray:
    system = _absorbing(chart, breakpoints, laws, resolution, sigma)
    return system.point_flux(x, _emission(system, q))


def _route(chart: str, breakpoints, laws, x: float, g: int, sigma=_SIGMA, q=_Q) -> float:
    """The mpmath route of the point's flux: the chart's one-dimensional integral (``_characteristic_point_mp``)."""
    s, qq = tuple(sigma[g]), tuple(q[g])
    match chart:
        case "sphere":
            inner = laws[0][0] if breakpoints[0] > 0.0 else 0.0
            if x == 0.0:
                return float(mpref.phi_sphere_centre(breakpoints, s, qq, laws[-1][0]))
            return mpref.phi_sphere(breakpoints, s, qq, laws[-1][0], inner, x)
        case "slab":
            return mpref.phi_slab(breakpoints, s, qq, laws[0][0], laws[1][0], x)
        case "cylinder":
            if x == 0.0:
                return mpref.phi_cylinder_axis(breakpoints, s, qq, laws[-1][0])
            if laws[-1][0] == 0.0:
                return mpref.phi_cylinder_vacuum(breakpoints, s, qq, x)
            return mpref.phi_cylinder(breakpoints, s, qq, laws[-1][0], x)
    raise AssertionError(chart)


def _check_route(chart, breakpoints, laws, x, band=_GIVEN_BAND, **kwargs) -> None:
    got = _given(chart, breakpoints, laws, x, **kwargs)
    data = {k: v for k, v in kwargs.items() if k in ("sigma", "q")}
    for g in range(len(got)):
        expected = _route(chart, breakpoints, laws, x, g, **data)
        assert abs(got[g] / expected - 1.0) < band, (chart, laws, x, g, got[g], expected)


# ── M. the point rule's measure ──────────────────────────────────────────

_POINTS = [
    ("sphere", _MR3, (_VACUUM,), x) for x in (0.0, 0.3, 0.5, 1.1, 2.0)
] + [
    ("sphere", _MR3H, (_VACUUM, (0.8, 0.0)), x) for x in (0.4, 1.7, 2.0)
] + [
    ("slab", _SLB3, ((0.3, 0.0), (0.8, 0.0)), x) for x in (0.0, 0.7, 2.3)
] + [
    ("cylinder", _MR3, ((0.6, 0.0),), x) for x in (0.0, 1.1, 2.0)
]
_POINT_IDS = [f"{c}{len(b) - 1}{'h' if b[0] > 0 else ''}-x{x}" for c, b, _, x in _POINTS]


def _kept(chart: str) -> int:
    return {"sphere": 3, "cylinder": 2, "slab": 1}[chart]


def _rule_geometry(chart, breakpoints, laws, x):
    """The point rule at x for group 0, its lines' directions ``(L, 3)`` and the points on them ``(L, q, 3)``."""
    basis = _basis(chart, breakpoints, 3, 2, 0.4)
    rule = PointRule.of(basis, _walls(chart, breakpoints, laws), np.asarray(_SIGMA[0]), x, 8)
    chart_object = basis.regions.chart
    lines = chart_object.line_domain().lines(rule.lines.coordinates)
    image = chart_object.image(lines)
    if rule.lines.levels is not None and isinstance(image, RadialImage):
        image = image.at_level(*rule.lines.levels)
    t = image.parameters_at(np.array([x]), rule.side)                                    # (L, q)
    points = lines.foot[:, None, :] + t[..., None] * lines.direction[:, None, :]
    return rule, lines.direction, points


@pytest.mark.foundation
@pytest.mark.parametrize(("chart", "breakpoints", "laws", "x"), _POINTS, ids=_POINT_IDS)
def test_the_point_rule_integrates_the_direction_moments(chart, breakpoints, laws, x) -> None:
    """[M1; foundation, THEOREM] The point rule is a rule for dOmega / 4 pi at x: its weights, summed over the
    lines and the sides, are 1; off the singular stratum it integrates (Omega . r_hat)^2 to 1/3 and ^4 to 1/5, with
    r_hat the outward radial direction at the point (the slab's x_hat); on a cylinder (axis included) Omega_z^2 to 1/3.

    The moments are closed forms of the uniform measure on S^2; ``Omega . r_hat`` is read from the point the rule
    puts on each line (``parameters_at`` of the line's image at its level) and the line's direction, here, so a Jacobian,
    a side or a parameter error each moves a moment. On the sphere's centre the rule is the one line b = 0 (the
    flux is isotropic there), so only the total is a law. Band 1e-13 `[M]` 2026-10-08 (worst 2.2e-15,
    ``ta/measure_m1.log``). First reds: the sphere's Jacobian b / (c y) without its b (moments move O(1)); the
    cylinder's polar weight without sin(theta) (Omega_z^2 moves); one side dropped (the total halves).
    """
    rule, direction, points = _rule_geometry(chart, breakpoints, laws, x)
    weight = np.repeat(rule.lines.weights[:, None], rule.side.size, axis=1)                    # (L, q)
    assert abs(weight.sum() - 1.0) < 1e-13, ("total", weight.sum())
    kept = _kept(chart)
    radius = np.linalg.norm(points[..., :kept], axis=-1)
    assert np.all(np.abs(radius - x) <= 8 * np.spacing(max(x, 1.0))), ("the point is on the line", np.max(np.abs(radius - x)))
    if chart == "cylinder":
        assert abs(np.sum(weight * direction[:, None, 2] ** 2) - 1 / 3) < 1e-13, "Omega_z^2"
    if chart != "slab" and x == 0.0:
        return
    outward = np.zeros_like(points)
    if chart == "slab":
        outward[..., 0] = 1.0
    else:
        outward[..., :kept] = points[..., :kept] / radius[..., None]
    mu = np.sum(outward * direction[:, None, :], axis=-1)
    assert abs(np.sum(weight * mu**2) - 1 / 3) < 1e-13, ("mu^2", np.sum(weight * mu**2))
    assert abs(np.sum(weight * mu**4) - 1 / 5) < 1e-13, ("mu^4", np.sum(weight * mu**4))


@pytest.mark.foundation
@pytest.mark.parametrize(("chart", "breakpoints", "x"), [("sphere", _MR3, 2.0), ("cylinder", _MR3, 2.0), ("slab", _SLB3, 0.0),
                                                       ("slab", _SLB3, 2.3)], ids=["sphere-wall", "cylinder-wall", "slab-left", "slab-right"])
def test_on_the_outer_wall_the_point_rule_has_the_line_rules_lines(chart, breakpoints, x) -> None:
    """[M2 = C6, second clause; foundation, THEOREM] One constructor builds both direction rules: at a point on a
    wall (no new end inserted) the point rule's line coordinates are the line rule's, the same set bit for bit
    (``np.array_equal`` after sorting rows), for one group's cross sections.

    The rows are sorted because each rule orders its lines by projected speed and the slab's two signs, the
    cylinder's grid, keep their own order. First red: a point rule built by a second constructor (a plain rule
    over the direction box, the #516 hazard); a grading the point rule computes and the line rule does not (or
    the reverse: the rim law read by one of the two).
    """
    basis = _basis(chart, breakpoints, 3, 2, 0.4)
    laws = ((0.6, 0.0),) if chart != "slab" else ((0.3, 0.0), (0.8, 0.0))
    walls = _walls(chart, breakpoints, laws)
    sigma = np.asarray(_SIGMA[1])
    point = PointRule.of(basis, walls, sigma, x, 8).lines.coordinates
    line = LineRule.of(basis, walls, sigma, 8).lines.coordinates
    assert np.array_equal(np.unique(point, axis=0), np.unique(line, axis=0))


@pytest.mark.foundation
@pytest.mark.parametrize("chart", ["sphere", "cylinder", "slab"])
def test_a_line_set_carries_levels_iff_its_chart_is_radial(chart: str) -> None:
    """[M4; foundation] ``Lines`` refuses a radial chart's lines without their exact levels, and a slab's lines with
    levels (the slab's image has no half-chord), naming both (#590). Positive leg: the line rule's own set on each
    chart constructs. First red: the check deleted (a radial set without levels chords from b, the first build's
    F3 drift; a slab's levels are silently ignored)."""
    breakpoints = _SLB3 if chart == "slab" else _MR3
    laws = ((0.3, 0.0), (0.8, 0.0)) if chart == "slab" else ((0.6, 0.0),)
    good = LineRule.of(_basis(chart, breakpoints, 1, 0, 0.4), _walls(chart, breakpoints, laws), np.asarray(_SIGMA[0]), 4).lines
    n = good.weights.size
    wrong = None if good.levels is not None else (np.ones(n), np.zeros(n))
    with pytest.raises(ValueError, match="carry their exact levels and a slab's carry none"):
        Lines(good.basis, good.walls, good.sigma_t, good.coordinates, good.weights, wrong, good.chunk, good.budget)


@pytest.mark.foundation
@pytest.mark.parametrize(("chart", "breakpoints", "x"), [("sphere", _MR3, 2.0 + 1e-9), ("slab", _SLB3, -1e-9),
                                                       ("sphere", _MR3H, 0.3)], ids=["beyond-sphere", "beyond-slab", "in-cavity"])
def test_a_point_outside_the_body_is_refused(chart, breakpoints, x) -> None:
    """[M3; foundation] ``PointRule.of`` refuses a point outside ``[r_0, r_n]`` (beyond the wall, or in a hollow
    body's cavity), naming the body's interval. ``admit_observable`` refuses such a ``PointValue`` first; this row
    is the rule's own door. First red: the guard removed (the cavity point reads lines that never reach it)."""
    basis = _basis(chart, breakpoints, 3, 2, 0.4)
    laws = ((0.0, 0.0),) * (2 if chart == "slab" or breakpoints[0] > 0 else 1)
    with pytest.raises(ValueError, match="a point is read inside the body"):
        PointRule.of(basis, _walls(chart, breakpoints, laws), np.asarray(_SIGMA[0]), x, 8)


# ── C. the reading of a given emission ───────────────────────────────────

_C1 = [("sphere", _MR3, ((a, 0.0),), 0.0) for a in (0.0, 0.6, 1.0)] + [
    ("slab", _SLB3, laws, x) for laws in (((0.3, 0.0), (0.8, 0.0)), (_MIRROR, _VACUUM)) for x in (0.0, 2.3)
]


@pytest.mark.l1
@pytest.mark.verifies("characteristic-quadrature")
@pytest.mark.parametrize(("chart", "breakpoints", "laws", "x"), _C1,
                         ids=[f"{c}-{l}-x{x}".replace(" ", "") for c, _, l, x in _C1])
@pytest.mark.rests_on(_M1)
def test_the_reading_at_a_singular_point_is_the_closed_form(chart, breakpoints, laws, x) -> None:
    """[C1; l1, REFERENCE] At the sphere's centre (every direction a diameter, the elementary closed form) and on a
    slab's faces (the E_2 image series of the unfolded path, any albedos) the reading of the MR3 emission, two
    groups with q and Sigma_t jumping at both interfaces, equals the closed form to 1e-12 relative.

    `[M]` 2026-10-08 (``ta/measure_sphere.log``, ``measure_slab.log``): centre <= 1.4e-16; faces <= 4.4e-16. First
    reds: each group's emission paired with the other group's Sigma_t (every row); the slab's albedos swapped (the
    (0.3, 0.8) rows); the centre's rule given both sides at weight 1 (x2).
    """
    _check_route(chart, breakpoints, laws, x)


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("characteristic-quadrature")
@pytest.mark.parametrize("albedo", [0.0, 0.6, 1.0])
@pytest.mark.rests_on(_M1)
def test_the_reading_on_a_cylinders_axis_is_the_polar_integral(albedo: float) -> None:
    """[C1, cylinder; l1, REFERENCE] On the axis of the MR3 cylinder every direction's line has b = 0: the reading
    equals the polar integral of the in-plane diameter's closure (and, under vacuum, Bickley's Ki_2 closed form,
    which agrees with it to 1e-21 on the reference side) to 1e-12. Slow: a cylinder's block, 19 s per group at
    layers 0. First red: the axial cosine's weight without sin(theta)."""
    _check_route("cylinder", _MR3, ((albedo, 0.0),), 0.0, resolution=(3, 0, 8, 12))


_C2_SPHERE = [("sphere", _MR3, ((a, 0.0),), x) for a in (0.0, 0.6, 1.0) for x in (0.3, 0.5, 1.1, 1.5, 2.0)]
_C2_HOLLOW = [("sphere", _MR3H, laws, x) for laws in ((_VACUUM, _VACUUM), ((0.3, 0.0), (0.8, 0.0)), (_MIRROR, _MIRROR),
                                                       ((0.9, 0.0), _VACUUM)) for x in (0.4, 0.45, 1.7, 2.0)]
_C2_SLAB = [("slab", _SLB3, laws, x) for laws in (((0.3, 0.0), (0.8, 0.0)), (_MIRROR, _VACUUM)) for x in (0.4, 0.7, 1.5)]
_C2_ROWS = _C2_SPHERE + _C2_HOLLOW + _C2_SLAB


@pytest.mark.l1
@pytest.mark.verifies("characteristic-quadrature")
@pytest.mark.parametrize(("chart", "breakpoints", "laws", "x"), _C2_ROWS,
                         ids=[f"{c}{'h' if b[0] else ''}-{l}-x{x}".replace(" ", "") for c, b, l, x in _C2_ROWS])
@pytest.mark.rests_on(_M1)
def test_the_reading_at_a_general_point_is_the_mpmath_route(chart, breakpoints, laws, x) -> None:
    """[C2; l1, REFERENCE] At general points of the MR3 sphere (inside a region, on both interfaces, on the wall),
    the MR3H hollow sphere (on the inner wall, beside it, inside, on the outer wall; the absorbing, partial,
    mirror and inner-only-mirror laws) and the SLB3 slab, the reading of the given emission equals the mpmath
    route to 1e-12 relative, both groups.

    The route integrates the backward path's closed-form passages over mu split at every tangency, the path
    built in the plane of the point and the direction (a disk billiard: the first leg, then the period of
    transits closed by its geometric sum), sharing no line, rule or closure with the code. `[M]` 2026-10-08:
    worst 6.7e-14 (sphere x = 1.5), hollow 5.4e-14 on the inner wall, slab 4.4e-16 (``ta/measure_*.log``). First
    reds: the albedo pairing (inner and outer amplitudes swapped in the period); the point read only on its inward
    side (the hollow and wall rows move O(1)); the tangency of an interface below the point not a panel end of
    the point rule (the rows at x > r_1 move, 1e-6-1e-8).
    """
    _check_route(chart, breakpoints, laws, x)


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("characteristic-quadrature")
@pytest.mark.parametrize(("albedo", "x", "line_points"), [(0.0, 1.1, 8), (0.0, 2.0, 16), (0.6, 1.1, 8)],
                         ids=["vacuum-1.1", "vacuum-wall", "partial-1.1"])
@pytest.mark.rests_on(_M1)
def test_the_reading_at_a_general_point_of_a_cylinder_is_the_mpmath_route(albedo: float, x: float, line_points: int) -> None:
    """[C2, cylinder; l1, REFERENCE] Off the axis of the MR3 cylinder: under vacuum the azimuthal integral of
    in-plane Ki_2 differences, under a partial mirror the two-dimensional route (alpha, theta) over the in-plane
    disk billiard, to 1e-12. Slow (a cylinder's block, and the reference's seconds). First red: the cylinder's
    in-plane angle Jacobian db / sqrt(c^2 - b^2) replaced by the sphere's b db / (c sqrt(c^2 - b^2)).

    The vacuum wall reads at 16 line points: a line grazing the wall at polar angle theta has an exponential layer of
    width sin(theta) / (2 Sigma_out) in y, graded by the grading law at every impact-panel top (G) for every law. `[M]` 2026-10-08 (the main
    agent's re-run of ``ta/measure_cyl_vac_wall.py``): before that law 7.9e-8 (8 points) and 1.4e-9 (16); after it
    4.9e-13 (8 points, only 2x inside the band, so not the row's point) and 2.4e-15 (16). First red: the layer graded
    only for 0 < a < 1 (battery arm ``layer-only-for-a``)."""
    _check_route("cylinder", _MR3, ((albedo, 0.0),), x, resolution=(3, 0, line_points, 12))


_C2B = [("sphere", _MR3, ((a, 0.0),), x) for a in (0.6, 0.9, 0.99) for x in (2.0, 2.0 * (1 - 1e-6))] + [
    ("sphere", _MR3H, (_VACUUM, (0.99, 0.0)), 2.0),
    ("slab", _SLB3, ((0.3, 0.0), (0.8, 0.0)), 1e-6),
    ("slab", _SLB3, ((0.3, 0.0), (0.8, 0.0)), 2.3 - 1e-6),
    ("slab", _SLB3, ((0.3, 0.0), (0.8, 0.0)), 2.3 - 1e-3),
]


@pytest.mark.l1
@pytest.mark.verifies("characteristic-quadrature")
@pytest.mark.parametrize(("chart", "breakpoints", "laws", "x"), _C2B,
                         ids=[f"{c}{'h' if b[0] else ''}-{l}-x{x}".replace(" ", "") for c, b, l, x in _C2B])
@pytest.mark.rests_on(_C2)
def test_the_reading_on_and_near_a_partial_mirror_is_the_mpmath_route(chart, breakpoints, laws, x) -> None:
    """[C2b; l1, REFERENCE] On a sphere's outer wall under a partial mirror a in {0.6, 0.9, 0.99}, a millionth of
    the radius inside it, on a hollow sphere's outer wall at 0.99, and on the slab 1e-6 and 1e-3 from a face (1e-6
    to 2.4e-6 mean free paths), the reading equals the mpmath route to 1e-12 relative.

    The rim's closure 1/(1 - a e^{-2 Sigma y}) has a pole a distance -ln(a)/(2 Sigma_out) from y = 0 (0.002 at
    a = 0.99): the rule must grade toward it (the grading law G at the outermost panel top). `[M]` 2026-10-08 before that law: on the
    wall at 8 line points 1.3e-8 (a = 0.6), 1.6e-5 (0.9), 2.8e-4 (0.99); just inside, 3.2e-15 (the point's own
    panel end resolved it by accident) (``ta/measure_wall.log``). With it (the main agent's re-run): <= 2.8e-14
    at 8 points; at a = 0.99 the miss grows with the count (5.9e-14 at 10, 2.1e-13 at 24) `[R]` as rounding
    amplified by 1/(1 - a) [REFUTED 2026-10-08 by qa's F3: the line's b loses the half-chord near grazing; the level
    carries it, Q3]. First red: the grading law G returning
    infinity (``HEAD`` ``a5492017``'s rule): the wall rows at 0.6/0.9/0.99 red, the near-wall rows stay green
    (declared: the point's own end grades them).
    """
    _check_route(chart, breakpoints, laws, x)


@pytest.mark.l1
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.parametrize("albedo", [0.6, 0.9, 0.99])
@pytest.mark.rests_on(_ESCAPE)
def test_the_blocks_total_under_a_partial_mirror_is_the_closed_form(albedo: float) -> None:
    """[C2b, the block; l1, REFERENCE] 1^T K 1 of a homogeneous sphere (Sigma 2.4, R 2) behind a partial mirror a
    equals the mpmath integral over lines (``_characteristic_mp.specular_sphere_total``) to 1e-14 relative, at the
    working resolution. The rim's pole reaches the block: `[M]` 2026-10-08 before the rim law, at 8 / 16 / 32
    line points: a = 0.9 1.0e-8 / 7.3e-12 / 2.6e-16; a = 0.99 1.6e-9 / 1.8e-10 / 2.1e-12, algebraic
    (``ta/measure_rim_block.log``); after it, <= 5.1e-16 at every count from 8 (the main agent). First red:
    the grading law G returning infinity."""
    basis = _basis("sphere", (0.0, 2.0), 3, 2, 0.4)
    block = LineRule.of(basis, _walls("sphere", (0.0, 2.0), ((albedo, 0.0),)), np.array([2.4]), 8).transport(12, 12).block
    ones = np.ones(basis.size)
    expected = float(specular_sphere_total(2.4, 2.0, albedo))
    assert abs(float(ones @ block @ ones) / expected - 1.0) < 1e-14


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.parametrize("albedo", [0.9, 0.99])
@pytest.mark.rests_on(_ESCAPE)
def test_on_a_cylinders_partial_mirror_the_block_and_the_wall_reading_are_the_closed_forms(albedo: float) -> None:
    """[C2b, the cylinder; l1, REFERENCE] A homogeneous cylinder (Sigma 2.4, R 2, q 0.7) behind a partial mirror: its
    1^T K 1 equals the mpmath integral over lines (b, theta) to 1e-14, and its reading ON the wall equals the
    two-dimensional route to 1e-12, at the working resolution.

    On the cylinder the rim's pole sits at -ln(a) sin(theta) / (2 Sigma), so each polar angle's impact rule is graded
    toward it at its own sin(theta) (#587), and below sqrt(2 R ulp(R)) in y the impact parameter rounds onto R (a
    tangent line with no crossing). `[M]` 2026-10-08, before #587, when every angle was graded at the slowest
    (``ta/measure_cylinder_rim.log``): 1^T K 1 1.6e-15 at both albedos; the wall
    4.3e-14 (0.9) and 3.4e-13 (0.99); 0 tangent lines among the point rule's. Positive control (the grading's floor
    at 0, ``measure_cylinder_rim_floor0.log``): 4 (0.9) and 52 (0.99) tangent lines, and the wall reading raises
    "lies on a transit". [The floor retires with the level (K, #590); this row's first red becomes the level
    dropped.] First reds: the level not passed by the line set (the wall rows raise or drift); G infinite (the wall
    rows move). Slow: 30 s per albedo.
    """
    basis = _basis("cylinder", (0.0, 2.0), 3, 2, 0.4)
    system = _absorbing("cylinder", (0.0, 2.0), ((albedo, 0.0),), sigma=((2.4,),))
    ones = np.ones(basis.size)
    total = float(ones @ system.groups[0].block @ ones)
    assert abs(total / mpref.specular_cylinder_total(2.4, 2.0, albedo) - 1.0) < 1e-14
    wall = float(system.point_flux(2.0, system.emission.restrict(np.full((1, basis.size), 0.7)))[0])
    assert abs(wall / mpref.phi_cylinder((0.0, 2.0), (2.4,), (0.7,), albedo, 2.0) - 1.0) < 1e-12


def _per_region(breakpoints, values: Callable[[int, np.ndarray], np.ndarray], groups: int = 2):
    """A function of the orbit coordinate, one formula per region: ``(q,) -> (G, q)``, for ``PanelBasis.project``."""
    inner = np.asarray(breakpoints[1:-1])

    def f(c: np.ndarray) -> np.ndarray:
        region = np.searchsorted(inner, c, side="right")
        return np.stack([values(g, c) * (1.0 + 0.5 * region) for g in range(groups)])
    return f


def _outside(g: int, c: np.ndarray) -> np.ndarray:
    return (1.0 + g) * (1.5 + np.cos(3.0 * c))


@pytest.mark.l2
@pytest.mark.parametrize("x", [0.3, 1.1, 2.0])
@pytest.mark.rests_on(_C2)
def test_the_reading_of_an_emission_outside_the_basis_converges_in_the_degree(x: float) -> None:
    """[C4; l2, CONV] The reading of an emission outside the basis (1.5 + cos 3c times a per-region factor, its
    L2 projection on the panels) converges as the degree rises 1 -> 5 on the MR3 sphere (a = 0.6): the change
    against degree 6 falls monotonically and by at least 1e3 from degree 1 to 5.

    The reading is exact for a q in the basis (C1, C2), so this row measures only the projection error of q, which
    the transport carries to the point. `[M]` 2026-10-08 (``ta/measure_c4.log``). First red: the projection's load
    taken at the panel's nodes only (an interpolation, not the L2 projection: the convergence stalls). Declared:
    no mpmath value leg (a nested quadrature of a non-constant q per region, minutes); C2 carries the value.
    """
    readings = []
    for degree in (1, 2, 3, 4, 5, 6):
        system = _absorbing("sphere", _MR3, ((0.6, 0.0),), (degree, 2, 8, 12))
        coefficients = system.basis.project(_per_region(_MR3, _outside), 2 * (degree + 1))
        readings.append(system.point_flux(x, system.emission.restrict(coefficients)))
    change = np.array([np.max(np.abs(r / readings[-1] - 1.0)) for r in readings[:-1]])
    assert np.all(np.diff(change) < 0.0), change
    assert change[-1] < 1e-3 * change[0], change


# ── C6, C8: one transport, two test measures ─────────────────────────────


def _volume_rule(ends: np.ndarray, chart: str, layers: int = 4, points: int = 8) -> tuple[np.ndarray, np.ndarray]:
    """Gauss-Legendre on each panel, graded geometrically (ratio 1/4) toward both ends, in the chart's measure, by hand.

    The flux of a panel function has a logarithmic derivative at the panel's ends (where the emission jumps), so
    the rule is graded there; written here, sharing nothing with the basis's own panel rule.
    """
    x, w = np.polynomial.legendre.leggauss(points)
    density = {"sphere": lambda c: 4.0 * np.pi * c**2, "cylinder": lambda c: 2.0 * np.pi * c, "slab": np.ones_like}[chart]
    pts, wts = [], []
    for a, b in zip(ends[:-1], ends[1:]):
        h = b - a
        cuts = np.unique(np.r_[[a, b, 0.5 * (a + b)], [a + 0.5 * h * 0.25**k for k in range(1, layers + 1)],
                              [b - 0.5 * h * 0.25**k for k in range(1, layers + 1)]])
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            c = 0.5 * (hi - lo) * x + 0.5 * (lo + hi)
            pts.append(c)
            wts.append(0.5 * (hi - lo) * w * density(c))
    return np.concatenate(pts), np.concatenate(wts)


@lru_cache(maxsize=None)
def _rows(chart: str, breakpoints, laws, sigma: float) -> tuple[Any, np.ndarray, np.ndarray, np.ndarray]:
    """One group's block, the volume rule's points and weights, and the point rows there ``(X, M)``, at degree 2."""
    basis = _basis(chart, breakpoints, 2, 0, 0.4)
    walls = _walls(chart, breakpoints, laws)
    sig = np.full(len(breakpoints) - 1, sigma) if np.isscalar(sigma) else np.asarray(sigma)
    transport = LineRule.of(basis, walls, sig, 8).transport(12, 12)
    xs, ws = _volume_rule(np.asarray(basis.partition.breakpoints), chart, layers=8, points=12)
    rows = np.array([PointRule.of(basis, walls, sig, x, 8).row(transport, TransportResolution(8, 12, 12)) for x in xs])
    return transport, xs, ws, rows


_C6 = [("sphere", _MR3, ((0.6, 0.0),), (1.7, 0.35, 2.4)),
       pytest.param("sphere", _MR3, (_WHITE,), (1.7, 0.35, 2.4), marks=pytest.mark.slow),
       pytest.param("sphere", _MR3, ((0.0, 0.5),), (1.7, 0.35, 2.4), marks=pytest.mark.slow)]
_C6_IDS = ["sphere3-partial", "sphere3-white", "sphere3-white-half"]


@pytest.mark.l1
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.parametrize(("chart", "breakpoints", "laws", "sigma"), _C6, ids=_C6_IDS)
@pytest.mark.rests_on(_C2, _CONSERVES)
def test_the_volume_integral_of_the_reading_is_the_galerkin_block(chart, breakpoints, laws, sigma) -> None:
    """[C6; l1, THEOREM] One transport, two test measures: the reading of each basis function u_j, integrated
    against u_i over the body by a volume rule written here, equals the block entry K_ij, to 1e-11 relative to
    max|K|, under a partial mirror and white walls (the diffuse part: the reading's r(x) currents against the
    block's R currents).

    ``W phi_h = K q`` says K_ij = int u_i (K u_j); the block integrates psi over the lines of space, the reading
    over the directions at each point, with one traversal transport. `[M]` 2026-10-08 (``ta/measure_c6b.log``): the
    volume rule (8 geometric layers of ratio 1/4 toward each panel end, 12 points per piece, 648 points) leaves
    1.2e-12 to 1.3e-12 under vacuum, a = 0.6 and a = 1; at 4 layers and 8 points it left 2e-9 to 3e-9, so the band
    is the rule's floor, not the reading's. The closed white body is slow (11 s of rows per fixture).
    First reds: the reading built on a second transport (#516's two answers); the diffuse walls' currents
    transposed or applied without alpha in the reading only (the white rows); the unit current injected as
    1/A in place of 1/D (factor 4, white rows).
    """
    transport, xs, ws, rows = _rows(chart, breakpoints, laws, sigma)
    basis = _basis(chart, breakpoints, 2, 0, 0.4)
    ends = np.asarray(basis.partition.breakpoints)
    panel = np.clip(np.searchsorted(ends, xs, side="right") - 1, 0, basis.n_panels - 1)
    tests = np.zeros((xs.size, basis.size))
    tests[np.arange(xs.size)[:, None], basis.columns(panel)] = basis.values(xs, panel)
    measured = (tests * ws[:, None]).T @ rows
    block = transport.block
    assert np.max(np.abs(measured - block)) < 1e-11 * np.max(np.abs(block))


@pytest.mark.l1
@pytest.mark.parametrize(("chart", "breakpoints", "laws", "sigma"), _C6, ids=_C6_IDS)
@pytest.mark.rests_on(_C2)
def test_the_region_transfers_read_through_points_are_reciprocal(chart, breakpoints, laws, sigma) -> None:
    """[C8; l1, THEOREM] Reciprocity through the readings: T_ij = int_{region i} (K 1_j)(x) dV, the flux in region i
    of a unit emission density in region j read at the volume rule's points, is symmetric to 1e-11 relative to
    max|T|.

    Independent of C6 and of the block's symmetry (C11): T is built from point readings alone. Declared
    stabiliser: a uniform scale of the reading (C1, C2 hold it). `[M]` 2026-10-08 (``ta/measure_c6.log``). First
    reds: the sphere's weight b / (c y) read with r^2 in place of c^2 (an x-dependent scale, asymmetric between
    regions); the diffuse response at the point taken from the escape U (the block's other face) in place of r(x).
    """
    transport, xs, ws, rows = _rows(chart, breakpoints, laws, sigma)
    region = _basis(chart, breakpoints, 2, 0, 0.4).region[transport.support]
    n = len(breakpoints) - 1
    at = np.clip(np.searchsorted(np.asarray(breakpoints), xs, side="right") - 1, 0, n - 1)
    unit = np.stack([(region == j).astype(float) for j in range(n)])                 # (n, M)
    transfer = np.array([[np.sum(ws[at == i] * (rows[at == i] @ unit[j])) for j in range(n)] for i in range(n)])
    assert np.max(np.abs(transfer - transfer.T)) < 1e-11 * np.max(np.abs(transfer)), transfer


_SPLITS = [(0.0, 0.5, 1.0, 1.5, 2.0), (0.0, 0.5, 0.8, 1.2, 1.5, 2.0)]


@pytest.mark.l1
@pytest.mark.parametrize("split", _SPLITS, ids=["middle-in-2", "middle-in-3"])
@pytest.mark.rests_on(_C2)
def test_an_interface_between_equal_materials_is_invisible_to_the_reading(split) -> None:
    """[C7, the reading's leg; l1, THEOREM] Splitting the MR3 sphere's middle region into 2 and 3 regions of the
    same material and emission leaves the reading at x in {0.3, 0.8, 1.1, 1.7, 2.0} unchanged, to 1e-13 relative
    (a = 0.6, both groups).

    The two bodies have different panels, line rules and point rules; both represent the given emission exactly.
    `[M]` 2026-10-08 (``ta/measure_c7.log``). First red: the attenuation restarted at each slot (tau not carried
    across a crossing).
    """
    pieces = len(split) - 3                                   # the middle region's pieces
    sigma = tuple(tuple(s[:1] + (s[1],) * pieces + s[2:]) for s in _SIGMA)
    q = tuple(tuple(v[:1] + (v[1],) * pieces + v[2:]) for v in _Q)
    for x in (0.3, 0.8, 1.1, 1.7, 2.0):
        whole = _given("sphere", _MR3, ((0.6, 0.0),), x)
        cut = _given("sphere", split, ((0.6, 0.0),), x, sigma=sigma, q=q)
        assert np.max(np.abs(cut / whole - 1.0)) < 1e-13, (x, cut, whole)


# ── W. the currents ──────────────────────────────────────────────────────

_WHITES = [("sphere", _MR3, (_WHITE,)), ("sphere", _MR3, ((0.0, 0.5),)), ("slab", _SLB3, ((0.0, 0.5), (0.0, 0.8))),
           ("slab", _SLB3, (_WHITE, _VACUUM))]


@pytest.mark.foundation
@pytest.mark.parametrize(("chart", "breakpoints", "laws"), _WHITES, ids=["sphere-white", "sphere-white-half", "slab-two-white", "slab-white-vacuum"])
def test_the_block_folds_its_walls_through_the_currents(chart, breakpoints, laws) -> None:
    """[W1; foundation] The block is ``line + response @ currents``, bit for bit (``WallCoupling.on_emission``, the
    one fold of the walls onto the emission), with ``currents`` of shape (W, M): the block and the reading share one
    definition of the returned currents (X4, Pattern 2). With no diffuse wall the currents are (0, M) and the block is
    its line part. (Re-posed 2026-10-08 when ``WallCoupling.update`` retired: the elegance review's second round.)

    Whether the refactor moved the block is the battery's pre-carve comparison (``ta/battery``: the block at
    ``HEAD`` ``a5492017`` against the tree, ``array_equal``), not a row: the suite holds no frozen block. First red:
    the reading's currents computed by a second solve (a different balance row): the reading's rows (C6 white)
    move, this row stays green (declared: it pins the block's spelling).
    """
    system = _absorbing(chart, breakpoints, laws)
    for group in system.groups:
        coupling = group.coupling
        assert coupling.currents.shape == (coupling.walls.breakpoint.size, group.support.size)
        assert np.array_equal(group.block, group.line + coupling.response @ coupling.currents)
    vacuum = _absorbing("sphere", _MR3, (_VACUUM,)).groups[0]
    assert vacuum.coupling.currents.shape == (0, vacuum.support.size) and np.array_equal(vacuum.block, vacuum.line)


@pytest.mark.l1
@pytest.mark.parametrize(("chart", "laws", "x"), [
    ("sphere", (_WHITE,), 1.1), ("sphere", ((0.0, 0.5),), 0.0), ("sphere", ((0.0, 0.5),), 1.1), ("sphere", ((0.0, 0.5),), 2.0),
    ("slab", ((0.0, 0.5), (0.0, 0.8)), 0.3), ("slab", ((0.0, 0.5), (0.0, 0.8)), 2.3), ("slab", (_WHITE, _WHITE), 1.1),
], ids=["sphere-closed-white", "sphere-half-centre", "sphere-half-1.1", "sphere-half-wall", "slab-0.3", "slab-face",
        "slab-closed"])
@pytest.mark.rests_on(_ESCAPE)
def test_the_reading_behind_white_walls_is_the_escape_closed_form(chart, laws, x) -> None:
    """[W2; l1, REFERENCE] Behind white walls, a homogeneous body (Sigma 0.8 and 2.4, q 1.7 and 0.9) reads the first
    flight plus each wall's re-entering current, ``j = a (E + T j)`` from Hebert's escape and transmission
    probabilities, attenuated to the point (``_characteristic_point_mp.phi_white_*``), to 1e-12; a closed white
    body reads q / Sigma flat.

    The diffuse part enters the reading only through ``currents``: a reading without it is the vacuum one, O(1)
    below. `[M]` 2026-10-08 (``ta/measure_white.log``). First reds: the diffuse term dropped from the row; the
    currents without alpha; a wall's injection read at the other wall (the two-white slab rows).
    """
    breakpoints = (0.0, 2.0) if chart == "sphere" else (0.0, 2.3)
    sigma, q = ((0.8,), (2.4,)), ((1.7,), (0.9,))
    got = _given(chart, breakpoints, laws, x, sigma=sigma, q=q)
    for g in range(2):
        if chart == "sphere":
            expected = mpref.phi_white_sphere(sigma[g][0], q[g][0], 2.0, laws[0][1], x)
        else:
            expected = mpref.phi_white_slab(sigma[g][0], q[g][0], 2.3, laws[0][1], laws[1][1], x)
        assert abs(got[g] / expected - 1.0) < 1e-12, (g, got[g], expected)


# ── Q. qa's findings on the first build (``p1_step_b5b/qa.md`` F1-F5), after the widened rung (K, G, R) ──

_VOID_OUTER = ((0.0, 1.0, 2.0), ((1.0, 0.0), (2.0, 0.0)), ((1.0, 0.0), (0.5, 0.0)))
_THIN_OUTER = ((0.0, 1.99, 2.0), ((1.0, 0.01), (2.0, 0.02)), ((1.0, 0.3), (0.5, 0.2)))
_Q_REF = "tests/gates/derivations/test_characteristic_reading.py::test_the_reading_on_and_near_a_partial_mirror_is_the_mpmath_route"


@pytest.mark.l1
@pytest.mark.verifies("characteristic-quadrature")
@pytest.mark.parametrize(("fixture", "x"), [(_VOID_OUTER, 1.0), (_VOID_OUTER, 2.0), (_THIN_OUTER, 1.99), (_THIN_OUTER, 2.0)],
                         ids=["void-outer-interface", "void-outer-wall", "thin-outer-interface", "thin-outer-wall"])
@pytest.mark.rests_on(_Q_REF)
@pytest.mark.catches("ERR-104")
def test_a_void_or_thin_outer_shell_under_a_near_one_mirror_is_read_to_the_bar(fixture, x: float) -> None:
    """[Q1 = qa F1; l1, REFERENCE] Under a = 0.99 with a void outer shell (Sigma = 0, q = 0) or an optically thin one
    (Sigma 0.01-0.02), the closure's pole moves to the interior tangency below it, a distance
    (-ln a + tau_above) |P Omega|_min / (2 Sigma_k) from y = 0: the reading at that interface and at the wall equals
    the mpmath route to 1e-12, both groups, at the working resolution.

    `[M]` qa 2026-10-08 (``qa/p1_hidden_pole.log``) under the first build's rim law (outermost panel only): void
    outer interface 5.1e-4 / 2.6e-5 / 5.8e-8 at 8 / 16 / 32 line points; wall 3.0e-6; thin outer 1.0e-6. With G
    (``ta/measure_widened.log``): 7.1e-15, 2.2e-16, 4.1e-14, 6.7e-16 at 8 line points. First red:
    the grading law G applied at the outermost panel top only (the first build's rim law).
    """
    breakpoints, sigma, q = fixture
    _check_route("sphere", breakpoints, ((0.99, 0.0),), x, sigma=sigma, q=q)


@pytest.mark.l1
@pytest.mark.verifies("characteristic-galerkin-assembly")
@pytest.mark.parametrize("fixture", [_VOID_OUTER, _THIN_OUTER], ids=["void-outer", "thin-outer"])
@pytest.mark.rests_on(_ESCAPE)
@pytest.mark.catches("ERR-104")
def test_the_blocks_total_under_a_void_or_thin_outer_shell_is_the_closed_form(fixture) -> None:
    """[Q1, the block; l1, REFERENCE] 1^T K 1 (group 0) of the layered sphere under a = 0.99 equals the per-line
    closed form integrated over b in the half-chord variable (``specular_sphere_layered_total``: the chord's vacuum
    total plus the entry response times the closure, void segments as q l) to 1e-13. The reference reproduces
    ``specular_sphere_total`` on a homogeneous sphere exactly and the point route's volume integral
    (``ta/check_layered_total.log``). `[M]` qa (``qa/p2_block.log``), first build, void outer: 2.6e-6 at 8 line
    points; with G 1.1e-16 (void) and 2.2e-16 (thin) (``ta/measure_widened.log``). First red: G at the outermost
    panel top only."""
    breakpoints, sigma, _ = fixture
    basis = _basis("sphere", breakpoints, 3, 2, 0.4)
    block = LineRule.of(basis, _walls("sphere", breakpoints, ((0.99, 0.0),)), np.asarray(sigma[0]), 8).transport(12, 12).block
    emitting = tuple(1.0 if s > 0.0 else 0.0 for s in sigma[0])
    expected = mpref.specular_sphere_layered_total(breakpoints, sigma[0], emitting, 0.99)
    total = float(np.ones(block.shape[0]) @ block @ np.ones(block.shape[1]))
    assert abs(total / expected - 1.0) < 1e-13, (total, expected)


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("characteristic-quadrature")
@pytest.mark.rests_on(_M1)
def test_a_cylinders_interface_point_is_read_to_the_bar() -> None:
    """[Q2 = qa F2; l1, REFERENCE] The MR3 cylinder under vacuum, at the interface x = 1.5: the polar-speed layer
    e^{-2 Sigma y / sin theta} sits at every panel top, not only the wall; the reading equals the Ki_2 azimuthal route
    to 1e-12. `[M]` qa (``qa/p2_cyl.log``), first build, layers 0: 1.0e-6 at 8 line points, 1.8e-8 at 16. First red:
    G's layer term at the outermost panel top only. Slow (a cylinder's block)."""
    _check_route("cylinder", _MR3, (_VACUUM,), 1.5, resolution=(3, 0, 8, 12))


@pytest.mark.l1
@pytest.mark.parametrize("albedo", [0.999, 0.99999])
@pytest.mark.rests_on(_Q_REF)
@pytest.mark.catches("ERR-106")
def test_the_wall_reading_as_the_mirror_tends_to_one_holds_the_bar(albedo: float) -> None:
    """[Q3 = qa F3; l1, REFERENCE] A homogeneous sphere (Sigma 2.4 and 0.8, q 0.7 and 1.3, R 2) read ON its wall under
    a = 0.999 and 0.99999 equals the route to 1e-12: each line through the point carries its exact half-chord (the
    level), so the closure 1/(1 - a e^{-2 Sigma y}) is read at the rule's y, not at sqrt((R - b)(R + b)) of the
    rounded b. `[M]` qa (``qa/p4_near_one.log``, ``p8_b_rounding.log``), first build: 2.8e-12 to 6.6e-12 (0.999),
    1.2e-10 to 2.8e-10 (0.99999), not converging in the line points; the rule's own sum at the exact y 2.2e-16.
    With the level (``ta/measure_widened.log``): 1.1e-15 and 4.4e-16. First red: the level dropped (the chord reads b)."""
    _check_route("sphere", (0.0, 2.0), ((albedo, 0.0),), 2.0, sigma=((2.4,), (0.8,)), q=((0.7,), (1.3,)))


@pytest.mark.l1
@pytest.mark.parametrize("x", [pytest.param(1.0, marks=pytest.mark.slow), 0.5], ids=["wall", "interface"])
@pytest.mark.rests_on(_M1)
@pytest.mark.catches("ERR-105")
def test_a_small_cylinders_layer_is_graded_outside_slow(x: float) -> None:
    """[Q4 = qa F4; l1, REFERENCE] A small cylinder (radii 0.5, 1; Sigma (0.6, 2.4) and (1.7, 0.8); q (1, 0.7) and
    (0.4, 1.3); vacuum) read at its wall and at its interface, degree 2, layers 0, 8 line points, equals the Ki_2
    route to 1e-12: the fast catcher of G's layer term (on the sphere the layer is s = 1, graded by the exponential
    ends; only the cylinder's slow polar angles make it thin, each graded at its own sin(theta) since #587). `[M]`
    2026-10-08, before #587 (``ta/measure_widened.log``):
    interface 4.5e-14 (5 s, the route 10 s), wall 4.8e-15 (37 s, slow); at degree 1 and 6 points 5.5e-11 and 4.6e-12, so
    the resolution is the row's. First red: G's layer term removed (the pole term only)."""
    _check_route("cylinder", (0.0, 0.5, 1.0), (_VACUUM,), x, resolution=(2, 0, 8, 8),
                 sigma=((0.6, 2.4), (1.7, 0.8)), q=((1.0, 0.7), (0.4, 1.3)))


@pytest.mark.foundation
@pytest.mark.parametrize(("x", "albedo"), [(1.0, 0.0), (0.5, 0.0), (1.0, 0.9)], ids=["wall", "interface", "wall-mirror"])
@pytest.mark.rests_on(_ASSEMBLY + "test_each_polar_angle_carries_the_impact_rule_graded_at_its_own_speed")
def test_each_polar_angle_through_a_point_carries_the_impact_rule_graded_at_its_own_speed(x: float, albedo: float) -> None:
    """[Q4b, #587; foundation, THEOREM of the rule's construction] The cylinder's point rule is the line rule's
    ITERATED rule below the point: at each polar node theta its lines are exactly ``impact_rule`` over the panels below
    c, graded at the line's own projected speed sin theta by the hand-written law (``grading_distances`` of the
    assembly file), with their exact levels, ``array_equal``. The weights carry the point's measure and are pinned by
    the route rows (Q2, Q4, C2).

    Activation: the small cylinder of Q4 under vacuum (the layer) and a = 0.9 (the pole); at the interface only the
    inner panel carries nodes and its layer binds below the next radius out. `[M]` 2026-10-08
    (``p1_step_c/ta/probe_small.log``): 13, 12 and 21 distinct per-angle counts; the row asserts more than one. First
    reds as the line rule's row (``p1_step_c/ta/``): the tensor rule of ``11b263a5`` (127 of 128 polar nodes
    differ), and the arms ``slowest``, ``no_layer``, ``no_pole``."""
    breakpoints, sigma = (0.0, 0.5, 1.0), (0.6, 2.4)
    basis = _basis("cylinder", breakpoints, 2, 0, 0.4)
    lines = PointRule.of(basis, _walls("cylinder", breakpoints, ((albedo, 0.0),)), np.asarray(sigma), x, 8).lines
    below = breakpoints.index(x)          # the point is a panel end, so the panels below it are the first ``below``
    counts = assert_each_polar_angle_is_graded_at_its_own_speed(lines, breakpoints, sigma, albedo, 8, panels=below)
    assert len(set(counts)) > 1, f"the fixture does not distinguish the polar speeds: {counts}"


@pytest.mark.foundation
@pytest.mark.parametrize("x", [2.0 + 1e-9, 0.3], ids=["beyond-the-wall", "in-the-cavity"])
def test_a_point_outside_the_body_is_refused_with_one_message(x: float) -> None:
    """[Q5a = qa F5; foundation] ``point_flux`` and ``angular_flux`` refuse a point outside the body (beyond the wall;
    in a hollow body's cavity) with the SAME message, naming the body's interval ("inside the body"). First red: the
    angular flux left to the transport's "lies on a transit" (the first build)."""
    breakpoints = _MR3 if x > 1.0 else _MR3H
    laws = (_VACUUM,) if breakpoints[0] == 0.0 else (_VACUUM, _VACUUM)
    system = _absorbing("sphere", breakpoints, laws)
    emission = _emission(system, _Q)
    with pytest.raises(ValueError, match="a point is read inside the body") as by_point:
        system.point_flux(x, emission)
    with pytest.raises(ValueError, match="a point is read inside the body") as by_direction:
        system.angular_flux(x, np.array([[0.5]]), emission)
    assert str(by_point.value) == str(by_direction.value)


@pytest.mark.foundation
@pytest.mark.parametrize("mu", [1e-100, 1e-300], ids=["1e-100", "1e-300"])
def test_the_slabs_grazing_direction_below_its_resolution_is_refused(mu: float) -> None:
    """[Q5b = qa F5; foundation] On the slab, a direction whose cosine is below the resolution of its line's parameter
    is refused as the sphere's is, never read as 0. `[M]` qa (``qa/p9b_slab_grazing.log``), first build: 0 returned
    silently at |mu| <= 1e-100 (true q / Sigma / 4 pi: 0.0153, 0.455), 1.6e-3 off at 1e-15; refused below sqrt(eps)
    now. First red: the refusal removed (a silent 0)."""
    system = _absorbing("slab", _SLB3, ((0.3, 0.0), (0.8, 0.0)))
    with pytest.raises(ValueError, match="too large to locate the point"):
        system.angular_flux(0.7, np.array([[mu], [-mu]]), _emission(system, _Q))


@pytest.mark.l1
@pytest.mark.parametrize("c", [1e-12, 1e-100, 1e-150], ids=["1e-12", "1e-100", "1e-150"])
@pytest.mark.rests_on(_HERE + "test_the_reading_at_a_singular_point_is_the_closed_form")
def test_a_point_near_the_centre_reads_the_centre(c: float) -> None:
    """[Q5c = qa F5; l1, THEOREM] Points 1e-12 to 1e-150 from the sphere's centre (all normal floats) read the centre's
    closed form to 1e-12 (the flux is smooth at the centre: phi(c) - phi(0) = O(c^2)); no NaN weight (the file runs
    under ``error::RuntimeWarning``). `[M]` qa (``qa/p3_near_stratum.log``), first build: 1e-12 and 1e-100 to 5e-14;
    1e-160 a NaN weight; now 5.0e-14 at 1e-12, 1e-100 and 1e-150 (``ta/measure_widened.log``). First red: the
    refusal threshold lowered below sqrt(tiny) (a NaN weight under error::RuntimeWarning)."""
    got = _given("sphere", _MR3, ((0.6, 0.0),), c)
    for g in range(2):
        expected = float(mpref.phi_sphere_centre(_MR3, _SIGMA[g], _Q[g], 0.6))
        assert abs(got[g] / expected - 1.0) < 1e-12, (c, g)


@pytest.mark.foundation
@pytest.mark.parametrize("c", [1e-160, 1e-300, 1e-310, 5e-324], ids=["1e-160", "1e-300", "subnormal", "smallest"])
def test_a_point_below_the_smallest_normal_float_is_refused(c: float) -> None:
    """[Q5d = qa F5; foundation] An orbit coordinate whose half-chords underflow (below sqrt of the smallest normal
    float, about 1.5e-154) is refused before any rule is built, not left to an underflow (a NaN weight at 1e-160,
    "breakpoints must be 1-D", "need interval[0] < interval[1]" in the first build). First red: the refusal removed."""
    system = _absorbing("sphere", _MR3, ((0.6, 0.0),))
    with pytest.raises(ValueError, match="half-chords that underflow"):
        system.point_flux(c, _emission(system, _Q))


# ── F. the angular flux ──────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.parametrize(("x", "albedo"), [(1.1, 0.6), (2.0, 0.6), (0.45, 0.0), (0.0, 1.0)], ids=["1.1", "wall", "0.45", "centre"])
@pytest.mark.rests_on(_C2)
def test_the_angular_flux_is_the_backward_paths_integral(x: float, albedo: float) -> None:
    """[F1; l1, REFERENCE] ``GalerkinSystem.angular_flux`` per steradian, at directions of the point's box (the
    outward cosine mu on the MR3 sphere: grazing +-1e-3, the interface tangencies +-sqrt(1 - r_k^2/c^2) and
    beside them, +-0.5, +-1), equals the mpmath backward path ``psi_disk / 4 pi`` to 1e-12 relative.

    One direction, one line: the row reads the transport and the side of the closest approach the direction faces,
    without the point rule. Sampled beside each tangency (+-1e-7) and away from the rim (0.05 on the wall): `[M]`
    2026-10-08, exactly AT a tangency the direction misses by 3.7e-10 (x = 1.1) and 1.1e-9 (the wall), and at
    mu = -1e-3 on the wall by 1.9e-10, because b is stored through the line's moment (1 ulp) and the half-chord is
    sqrt(2 r ulp)-sensitive there: conditioning of the line's representation, not the transport. First reds: the side taken from Omega_y (always the inward side: the outward
    directions move O(1)); the per-steradian factor 1/4 pi dropped.
    """
    system = _absorbing("sphere", _MR3, ((albedo, 0.0),))
    tangent = [np.sqrt((1 - r / x) * (1 + r / x)) for r in _MR3[1:-1] if 0.0 < r < x]
    # beside each tangency, never on it: at the tangent direction b = r_k, and a float cosine moves b by an ulp, which
    # moves the chord's half-length there by sqrt(2 r_k ulp) ~ 1e-8 on both sides (measured: 3.7e-10 and 1.1e-9 apart)
    # and grazing at 1e-3 inside, at 0.05 on the wall (on the wall y = R |mu|, and b's ulp moves y by R ulp / y: 2e-10
    # relative at mu = 1e-3, measured); a direction below 1e-8 is F2's (the box rule grades there)
    grazing = 0.05 if x == _MR3[-1] else 1e-3
    cosines = [1.0, 0.5, grazing, *[t + 1e-7 for t in tangent], *[t - 1e-7 for t in tangent]]
    mu = np.array(sorted(set(cosines + [-c for c in cosines])))[:, None] if x > 0.0 else np.zeros((1, 0))
    psi = system.angular_flux(x, mu, _emission(system, _Q))                              # (k, G)
    amplitude = {"outer": mpref.mp.mpf(repr(albedo)), "inner": mpref.mp.mpf(0)}
    for k in range(psi.shape[0]):
        cosine = float(mu[k, 0]) if x > 0.0 else 1.0
        for g in range(2):
            with mpref.mp.workdps(25):
                expected = float(mpref.psi_disk(mpref._mp(_MR3), mpref._mp(_SIGMA[g]), mpref._mp(_Q[g]), amplitude,
                                                mpref.mp.mpf(repr(x)), mpref.mp.mpf(repr(cosine))) / (4 * mpref.mp.pi))
            assert abs(psi[k, g] / expected - 1.0) < 1e-12, (x, cosine, g, psi[k, g], expected)


#: The smallest |mu| the box rule reads on the slab: above its refusal, sqrt(eps) (Q5b).
_SLAB_FLOOR = 1e-7


def _box_rule(breaks: np.ndarray, points: int) -> tuple[np.ndarray, np.ndarray]:
    """Gauss-Legendre on each interval of ``breaks`` under the endpoint substitution u -> (1 - cos(pi u)) / 2, which
    makes a square-root end smooth (the P1 spec's attacker F1: 3.3e-15 where split-only left 1.0e-6). Written here."""
    x, w = np.polynomial.legendre.leggauss(points)
    u = 0.5 * (x + 1.0)
    s, ds = 0.5 * (1.0 - np.cos(np.pi * u)), 0.5 * np.pi * np.sin(np.pi * u) * 0.5 * w
    pts = np.concatenate([a + (b - a) * s for a, b in zip(breaks[:-1], breaks[1:])])
    wts = np.concatenate([(b - a) * ds for a, b in zip(breaks[:-1], breaks[1:])])
    return pts, wts


@pytest.mark.l1
@pytest.mark.parametrize(("chart", "breakpoints", "laws", "x"), [
    ("sphere", _MR3, ((0.6, 0.0),), 1.1), ("sphere", _MR3, ((0.6, 0.0),), 2.0), pytest.param("slab", _SLB3, ((0.3, 0.0), (0.8, 0.0)), 0.7, marks=pytest.mark.slow),
], ids=["sphere-1.1", "sphere-wall", "slab-0.7"])
@pytest.mark.rests_on(_HERE + "test_the_angular_flux_is_the_backward_paths_integral")
def test_the_angular_flux_integrated_over_the_box_is_the_reading(chart, breakpoints, laws, x) -> None:
    """[F2; l1, THEOREM] The point reading equals the integral of ``angular_flux`` over the point's direction box,
    ``density * sum w psi``, by a rule written here on the box's cosine, split at mu = 0 and at the tangency to
    EVERY panel end below the point (the emission is polynomial per panel, so psi kinks at each; the main agent's
    smoke: 9e-9 with splits at mu = 0 only), each interval under the endpoint substitution; to 1e-12 relative.

    Two line constructions (the line rule's lines through the point; the lines ``Line.through`` the point in each
    direction), one transport. `[M]` 2026-10-08 after the level (``ta/measure_f2_widened.log``): 4.4e-16 (x = 1.1),
    1.2e-14 (the wall, graded to 2^-40 now that every grazing direction is read), 8.9e-16 (slab, slow: its rule
    stops at _SLAB_FLOOR and adds the sliver in closed form). Before the level: 1.1e-15, 2.7e-14 (stopped at a
    representability floor), 4.4e-16. First reds: the point rule's weights
    without the box density's 1/4 pi normalisation (x 4 pi); the point rule's lines read at the closest approach
    in place of the point (every row).
    """
    system = _absorbing(chart, breakpoints, laws)
    emission = _emission(system, _Q)
    domain = system.basis.regions.chart.directions_at(x)
    ends = np.asarray(system.basis.partition.breakpoints)
    # psi has an exponential layer e^{-tau/|mu|} at mu -> 0 for every optical depth tau to a panel end: grade toward
    # it. ON a radial wall the level now reads every direction, so the old floor goes; the SLAB refuses |mu| below sqrt(eps)
    # (Q5b), so its rule stops at _SLAB_FLOOR and the sliver |mu| < _SLAB_FLOOR is added in closed form: there the
    # line runs inside the point's region for x / |mu| >> 1 mean free paths, so psi = q / (4 pi Sigma) of that region
    floor = _SLAB_FLOOR if chart == "slab" else 0.0
    toward_grazing = 2.0 ** -np.arange(1, 41)
    toward_grazing = toward_grazing[toward_grazing > 2.0 * floor]
    breaks = np.r_[-1.0, -toward_grazing, -floor, floor, toward_grazing, 1.0]
    if chart != "slab":
        breaks = np.r_[breaks, np.concatenate([domain.tangencies(r) for r in ends if 0.0 < r < x])]
    breaks = np.unique(breaks)
    mu, w = _box_rule(breaks, 16)
    keep = np.abs(mu) > floor                                     # the sliver |mu| < floor is the slab's closed form
    mu, w = mu[keep], w[keep]
    psi = system.angular_flux(x, mu[:, None], emission)                                  # (k, G)
    integral = domain.density * (w @ psi)
    if chart == "slab":
        k = int(np.searchsorted(np.asarray(breakpoints), x, side="right") - 1)
        integral = integral + domain.density * 2.0 * floor * np.array([_Q[g][k] / _SIGMA[g][k] for g in range(2)]) / (4.0 * np.pi)
    reading = system.point_flux(x, emission)
    assert np.max(np.abs(integral / reading - 1.0)) < 1e-12, (integral, reading)


@pytest.mark.l1
@pytest.mark.parametrize("albedo", [0.6, 1.0])
@pytest.mark.rests_on(_HERE + "test_the_angular_flux_is_the_backward_paths_integral")
def test_the_grazing_limit_on_the_wall(albedo: float) -> None:
    """[F4; l1, REFERENCE] On the sphere's wall, at |mu| = 1e-6 (a few times the representability floor), both
    branches of ``angular_flux`` equal the mpmath backward path within the conditioning band 10 eps / mu^2
    (b is stored through the line's moment: an ulp of b moves y = R |mu| by R eps / y; at mu = -1e-3 that read
    1.9e-10, F1): under a = 1 the inward branch tends to the outer region's q / Sigma, under a < 1 both tend to 0.
    `[M]` 2026-10-08 (``ta/measure_f.log``): a = 1 <= 4.4e-16; a = 0.6 1.6e-4 on all four, against eps / mu^2 = 2.2e-4:
    the conditioning, not the transport. First red: the inward branch read at the far end of the chord."""
    system = _absorbing("sphere", _MR3, ((albedo, 0.0),))
    mu = np.array([[-1e-6], [1e-6]])
    psi = system.angular_flux(2.0, mu, _emission(system, _Q))
    band = 10.0 * np.finfo(float).eps / 1e-12
    amplitude = {"outer": mpref.mp.mpf(repr(albedo)), "inner": mpref.mp.mpf(0)}
    for k in range(2):
        for g in range(2):
            with mpref.mp.workdps(25):
                expected = float(mpref.psi_disk(mpref._mp(_MR3), mpref._mp(_SIGMA[g]), mpref._mp(_Q[g]), amplitude,
                                                mpref.mp.mpf(2), mpref.mp.mpf(repr(float(mu[k, 0])))) / (4 * mpref.mp.pi))
            assert abs(psi[k, g] / expected - 1.0) < band, (albedo, mu[k, 0], g, psi[k, g], expected)


@pytest.mark.l1
@pytest.mark.parametrize("albedo", [0.6, 1.0])
@pytest.mark.rests_on(_HERE + "test_the_angular_flux_is_the_backward_paths_integral")
@pytest.mark.catches("ERR-106")
def test_a_direction_grazing_the_wall_is_read_through_its_level(albedo: float) -> None:
    """[F5, inverted at the widened rung; l1, REFERENCE] On the sphere's wall, the directions mu = +-2^-40: the first
    build refused them (b rounds onto R, a tangent line); each line through the point now carries its exact level
    (c, c |Omega_x| / |P Omega|), so both branches are read and equal the mpmath backward path to 1e-12.
    `[M]` 2026-10-08 (``ta/measure_widened.log``): 4.4e-16 (a = 0.6), 7.8e-16 (a = 1). The mpmath route itself needed 40
    extra digits here (``psi_disk``): at 25 digits b = r sqrt(1 - mu^2) lost the half-chord and the route, not the code,
    was off by a factor 4 (the positive control of the reference). First red: the level dropped (the line is tangent
    and the reading raises)."""
    system = _absorbing("sphere", _MR3, ((albedo, 0.0),))
    mu = np.array([[-(2.0**-40)], [2.0**-40]])
    psi = system.angular_flux(2.0, mu, _emission(system, _Q))
    amplitude = {"outer": mpref.mp.mpf(repr(albedo)), "inner": mpref.mp.mpf(0)}
    for k in range(2):
        for g in range(2):
            with mpref.mp.workdps(25):
                expected = float(mpref.psi_disk(mpref._mp(_MR3), mpref._mp(_SIGMA[g]), mpref._mp(_Q[g]), amplitude,
                                                mpref.mp.mpf(2), mpref.mp.mpf(repr(float(mu[k, 0])))) / (4 * mpref.mp.pi))
            assert abs(psi[k, g] / expected - 1.0) < 1e-12, (albedo, mu[k, 0], g, psi[k, g], expected)


@pytest.mark.l1
@pytest.mark.rests_on(_HERE + "test_the_angular_flux_is_the_backward_paths_integral")
def test_at_the_centre_the_angular_flux_is_isotropic_and_four_pi_of_it_is_the_reading() -> None:
    """[F3; l1, THEOREM] At the sphere's centre the direction box is a point: ``angular_flux`` at the empty
    coordinate, times 4 pi, is the reading, to 1e-14 relative. First red: the per-steradian factor dropped."""
    system = _absorbing("sphere", _MR3, ((0.6, 0.0),))
    emission = _emission(system, _Q)
    psi = system.angular_flux(0.0, np.zeros((0,)), emission)
    reading = system.point_flux(0.0, emission)
    assert np.max(np.abs(4.0 * np.pi * psi / reading - 1.0)) < 1e-14, (psi, reading)


# ── the door ─────────────────────────────────────────────────────────────

_K = CellCoefficient.every(Channel.FISSION_EMISSION)
_DOOR_RES = Resolution(3, 2, 0.4, TransportResolution(8, 12, 12), 8)
#: The coarse resolution of the rows that integrate door readings over the body: the identity they gate holds at
#: any resolution, so the cheapest that solves.
_COARSE = Resolution(1, 0, 0.4, TransportResolution(6, 8, 8), 4)


def _spec(geometry: StructuredGeometry, mixtures: dict[int, Any], question: Any) -> GeometrySpecification:
    return GeometrySpecification(Materials(mixtures), geometry, question)


@lru_cache(maxsize=None)
def _door(specification: GeometrySpecification, resolution: Resolution = _DOOR_RES) -> CharacteristicDerivation:
    return CharacteristicDerivation(specification, resolution)


def _point(specification: GeometrySpecification, x: float, g: int, resolution: Resolution = _DOOR_RES) -> float:
    with bypass():
        reading = _door(specification, resolution).evaluate(PointValue(x, g))
    assert type(reading) is Uncertified
    return float(reading.value)


def _door_volume_integral(specification, resolution, weight_of_region: np.ndarray) -> float:
    """sum_g int w_g(x) phi_g(x) dV over the door's point readings, by the graded volume rule written here."""
    derivation = _door(specification, resolution)
    breakpoints = np.asarray(specification.geometry.breakpoints, dtype=float)
    xs, ws = _volume_rule(np.asarray(derivation.basis.partition.breakpoints), "sphere", layers=4, points=8)
    region = np.clip(np.searchsorted(breakpoints, xs, side="right") - 1, 0, len(breakpoints) - 2)
    total = 0.0
    with bypass():
        for x, w, k in zip(xs, ws, region):
            for g in range(weight_of_region.shape[1]):
                if weight_of_region[k, g] != 0.0:
                    total += w * weight_of_region[k, g] * float(derivation.evaluate(PointValue(float(x), g)).value)
    return total


@pytest.mark.l1
@pytest.mark.parametrize(("outer", "albedo"), [(BC.vacuum, 0.0), (_partial(0.6), 0.6)], ids=["vacuum", "partial-0.6"])
@pytest.mark.rests_on(_C2)
def test_a_pure_absorbers_point_values_are_the_mpmath_route(outer, albedo: float) -> None:
    """[E4; l1, REFERENCE] Through the door: a ``FixedSource`` (the MR3 rates per region, both groups) in pure
    absorbers (no scattering, no fission) reads ``PointValue`` at the centre, x = 1.1 and the wall equal to C1/C2's
    routes, to 1e-12. The door composes with the reading: the emission IS the source. First reds: the source
    read as a detector (x 4 pi); the solve applied as (I + K S) (no effect here: declared, S = 0; E1 and D3 see it).
    """
    absorbers = {j: _mixture([_SIGMA[0][j], _SIGMA[1][j]], np.zeros((2, 2))) for j in range(3)}
    spec = _spec(StructuredGeometry.sphere(_MR3, (0, 1, 2), outer=outer), absorbers,
                 FixedSource(RegionwiseConstant(np.array(_Q).T)))
    for x in (0.0, 1.1, 2.0):
        for g in range(2):
            expected = _route("sphere", _MR3, ((albedo, 0.0),), x, g)
            assert abs(_point(spec, x, g) / expected - 1.0) < _GIVEN_BAND, (x, g)


_GARCIA = _spec(StructuredGeometry.sphere((0.0, 3.0, 5.0, 7.0), (0, 1, 2), outer=BC.vacuum),
                {0: _mixture([1.0], [[0.99]]), 1: _mixture([0.5], [[0.3]]), 2: _mixture([2.0], [[1.9]])},
                FixedSource(RegionwiseConstant(np.array([[0.5], [1.0], [1.5]]))))
_GARCIA_RES = Resolution(4, 3, 0.4, TransportResolution(12, 16, 16), 10)
#: The reference's own step at Garcia's points, 10 x the largest change (4, 3, 12, 16) -> (5, 4, 16, 20): `[M]`
#: 2026-10-08, 5.4e-6 at r = 0.5 cm... 5.5 cm (``ta/measure_garcia.log``), so 10 x is 5.4e-5, rounded up.
_GARCIA_STEP = 6e-5


@pytest.mark.l1
@pytest.mark.catches("ERR-090")
@pytest.mark.parametrize(("r", "phi", "rounding"), list(zip(GARCIA_2021_CASE1_R, GARCIA_2021_CASE1_PHI, GARCIA_2021_CASE1_ROUNDING)),
                         ids=[f"r={r}" for r in GARCIA_2021_CASE1_R])
@pytest.mark.rests_on(_HERE + "test_a_pure_absorbers_point_values_are_the_mpmath_route")
def test_garcias_case_1_per_point(r: float, phi: float, rounding: float) -> None:
    """[E1; l1, REFERENCE] Garcia 2021 Case 1 (Williams 1991 Example 5: a three-region sphere, c = 0.99 in the core,
    vacuum) per point of Table 5, through the door: ``PointValue`` equals half the printed value (Garcia's flux is
    int Psi dmu with a per-steradian source; the door's rate source gives (2 pi / 4 pi) of it) within the table's
    half unit in its last digit plus 10 x the reference's own ladder step.

    `[M]` 2026-10-08 (``ta/measure_garcia.log``): at (4, 3, 12, 16) every point within that band; at the finest
    rung (5, 4, 16, 20) every point is inside the rounding alone (worst r = 3.0 cm, 3.0e-5 of 3.5e-5); the working
    point (3, 2, 8, 12) misses r = 6.0 cm by 1.6e-4, its own step. Today's family's bands were 3e-3 inside and
    4e-2 at the surface; this row is 30 to 600 x tighter. ``catches("ERR-090")``: a one-spline emission across the
    interfaces is unspellable here (the emission is per panel), so the marker records the row's coverage of the
    defect class and is re-dropped by the battery's interface arm. First reds: the scattering's emission omitted
    from the point (the source alone transported: O(1)); the reading's sides weighted unequally.
    """
    with bypass():
        derivation = _door(_GARCIA, _GARCIA_RES)
        value = float(derivation.evaluate(PointValue(float(r), 0)).value)
    expected = 0.5 * phi
    assert abs(value / expected - 1.0) < rounding + _GARCIA_STEP, (r, value, expected)


def _closed_mirror(mixture, question) -> GeometrySpecification:
    return _spec(StructuredGeometry.sphere((0.0, 1.0), (0,), outer=BC.reflective), {0: mixture}, question)


@pytest.mark.l1
@pytest.mark.parametrize("x", [0.0, 0.37, 1.0], ids=["centre", "inside", "wall"])
@pytest.mark.rests_on(_D5)
def test_a_closed_bodys_fundamental_reads_the_infinite_medium_flux_in_the_declared_gauge(x: float) -> None:
    """[D1; l1, THEOREM] A closed homogeneous sphere (mirror, PU2, two groups) reads at every point the infinite
    medium's flux: the dominant eigenvector of (diag Sigma_t - S)^-1 F, scaled so that its production over the body
    is 1 (the declared gauge, fission + (n,2n): PU2 has no (n,2n), asserted), written here from the Mixture arrays.

    Reads the answer's EMISSION on the gauge's scale: a door that scales the flux and not the emission is off by
    the gauge factor. `[M]` 2026-10-08 (``ta/measure_door.log``). First reds: the emission unscaled; the groups of
    the reading swapped.
    """
    assert not _PU2.Sig2[0].toarray().any()
    loss, production = _zero_d(_PU2)
    values, vectors = np.linalg.eig(np.linalg.solve(loss, production))
    v = np.abs(np.real(vectors[:, int(np.argmax(np.abs(values)))]))
    volume = 4.0 * np.pi / 3.0
    flux = v / (volume * float(np.asarray(_PU2.SigP) @ v))
    spec = _closed_mirror(_PU2, Eigen(_K))
    for g in range(2):
        assert abs(_point(spec, x, g) / flux[g] - 1.0) < 1e-12, (x, g)


@pytest.mark.l1
@pytest.mark.rests_on(_PAIRING)
def test_the_fundamentals_readings_carry_the_declared_production() -> None:
    """[D2; l1, THEOREM] On a heterogeneous vacuum sphere (PU2 | UP2N | PU2), the door's point readings, integrated
    against the production nu Sigma_f + 2 Sigma_2 per region by the volume rule written here, give 1: the declared
    gauge reaches the reading through the emission's scale. To 1e-9: `[M]` 2026-10-08 (``ta/measure_door.log``) the
    volume rule at 4 layers and 8 points leaves 1.9e-11 (3 layers, 6 points: 1.4e-9). First red: the emission left at the pencil's normalisation."""
    spec = _spec(StructuredGeometry.sphere(_MR3, (0, 1, 0), outer=BC.vacuum), {0: _PU2, 1: _UP2N}, Eigen(_K))
    production = np.array([np.asarray(_PU2.SigP), np.asarray(_UP2N.SigP), np.asarray(_PU2.SigP)])
    assert not any(m.Sig2[0].toarray().any() for m in (_PU2,)), "the gauge counts (n,2n): PU2 has none"
    n2n = 2.0 * _UP2N.Sig2[0].toarray().sum(axis=1)
    production[1] = production[1] + n2n
    assert abs(_door_volume_integral(spec, _COARSE, production) - 1.0) < 1e-9


@pytest.mark.l1
@pytest.mark.parametrize("x", [0.0, 0.6, 1.0], ids=["centre", "inside", "wall"])
@pytest.mark.rests_on(_E2)
def test_a_closed_body_with_a_uniform_source_reads_the_infinite_medium_flux_at_every_point(x: float) -> None:
    """[D3; l1, THEOREM] A closed mirror sphere of ABS (downscatter, no fission) with a uniform source Q reads
    (diag Sigma_t - S)^-1 Q at every point, to 1e-12 (the source answer's emission includes every collision's).
    First reds: the point read from the source alone (the scattered emission dropped); Sigma_s transposed."""
    source = np.array([1.0, 0.4])
    loss, _ = _zero_d(_ABS)
    flux = np.linalg.solve(loss, source)
    spec = _closed_mirror(_ABS, FixedSource(RegionwiseConstant(source[None, :])))
    for g in range(2):
        assert abs(_point(spec, x, g) / flux[g] - 1.0) < 1e-12, (x, g)


@pytest.mark.l1
@pytest.mark.rests_on(_PAIRING)
def test_a_detectors_importance_read_at_points_is_reciprocal_to_the_forward_flux() -> None:
    """[D4; l1, THEOREM] Reciprocity through the readings: on the heterogeneous vacuum sphere (ABS | UP2N | ABS, two
    groups, upscatter and downscatter, so the transposition matters), (1/4 pi) sum_g int Q_g phi^dagger_g dV, the
    ``Response`` answer's point readings integrated against the source by the volume rule written here, equals the
    forward ``FluxIntegral`` of the detector on the ``FixedSource`` answer, to 1e-9 (`[M]` the volume rule at 4
    layers and 8 points leaves 1.2e-10, ``ta/measure_door.log``). First reds: the adjoint posed with untransposed cross sections; the detector's
    lift without its 4 pi."""
    geometry = StructuredGeometry.sphere(_MR3, (0, 1, 0), outer=BC.vacuum)
    materials = {0: _ABS, 1: _UP2N}
    source = np.array([[1.0, 0.4], [0.0, 0.0], [0.5, 0.2]])
    detector = np.array([[0.3, 1.1], [2.0, 0.4], [0.7, 0.9]])
    forward = _spec(geometry, materials, FixedSource(RegionwiseConstant(source)))
    adjoint = _spec(geometry, materials, Response(RegionwiseConstant(detector)))
    with bypass():
        expected = float(_door(forward, _COARSE).evaluate(FluxIntegral(RegionwiseConstant(detector))).value)
    assert abs(_door_volume_integral(adjoint, _COARSE, source) / (4.0 * np.pi) / expected - 1.0) < 1e-9


@pytest.mark.l1
@pytest.mark.parametrize("x", [0.0, 1.1, 2.0])
@pytest.mark.rests_on(_HERE + "test_a_pure_absorbers_point_values_are_the_mpmath_route")
def test_in_one_group_the_importance_is_four_pi_times_the_forward_flux(x: float) -> None:
    """[D5; l1, THEOREM] One group: the transport is self-adjoint and the transposition the identity, so the
    ``Response`` answer's point value with detector f is 4 pi times the ``FixedSource`` answer's with source f
    (the pullback's 4 pi), to 1e-13. An edge equal to a foundation. First red: the detector lifted by the section."""
    mixtures = {j: _mixture([_SIGMA[0][j]], [[0.4 * _SIGMA[0][j]]]) for j in range(3)}
    geometry = StructuredGeometry.sphere(_MR3, (0, 1, 2), outer=_partial(0.6))
    f = RegionwiseConstant(np.array([[1.0], [0.25], [3.0]]))
    forward = _point(_spec(geometry, mixtures, FixedSource(f)), x, 0)
    importance = _point(_spec(geometry, mixtures, Response(f)), x, 0)
    assert abs(importance / (4.0 * np.pi * forward) - 1.0) < 1e-13


@pytest.mark.foundation
def test_a_point_value_of_a_mode_nearest_tau_is_refused_after_the_solve() -> None:
    """[D6; foundation] A ``PointValue`` of a ``Nearest`` answer raises NotImplementedError naming the missing flux
    scale, as its ``FluxIntegral`` does (rung 5a's N4); a ``PointValue`` of the fundamental is answered. First red:
    the mode's eigenvector read at the point on the pencil's normalisation."""
    spec = _spec(StructuredGeometry.sphere(_MR3, (0, 1, 0), outer=BC.vacuum), {0: _PU2, 1: _UP2N}, Eigen(_K, mode=Nearest(0.5)))
    with bypass(), pytest.raises(NotImplementedError, match="has no flux scale"):
        _door(spec, _COARSE).evaluate(PointValue(1.0, 0))


@pytest.mark.l1
@pytest.mark.rests_on(_HERE + "test_a_closed_bodys_fundamental_reads_the_infinite_medium_flux_in_the_declared_gauge")
def test_a_ratio_of_point_values_is_the_quotient_and_gauge_free() -> None:
    """[D7; l1, THEOREM] ``read(Ratio(PointValue(a, 0), PointValue(b, 1)))`` is the quotient of the two readings, bit
    for bit; on the closed PU2 sphere it is the infinite medium's group ratio phi_0 / phi_1 at any two points, to
    1e-12, independent of the gauge. First red: a ratio read off two different solves (the memo keyed on the
    observable alone)."""
    spec = _closed_mirror(_PU2, Eigen(_K))
    solution = characteristic_reference(spec, _DOOR_RES)
    with bypass():
        ratio = solution.read(Ratio(PointValue(0.2, 0), PointValue(0.9, 1))).value
        numerator = solution.read(PointValue(0.2, 0)).value
        denominator = solution.read(PointValue(0.9, 1)).value
    assert ratio == numerator / denominator
    loss, production = _zero_d(_PU2)
    values, vectors = np.linalg.eig(np.linalg.solve(loss, production))
    v = np.abs(np.real(vectors[:, int(np.argmax(np.abs(values)))]))
    assert abs(ratio / (v[0] / v[1]) - 1.0) < 1e-12
