"""Gates for :meth:`orpheus.geometry.chart.Chart.directions_at` and its :class:`DirectionDomain`.

Migration step (a) of P1 of the characteristic-reference campaign
(``.claude/plans/characteristic_reference_architecture.md``, "P1 API sketch"
item 7 and the ruling "Directions at a point: chart verb"; verification spec
``scratch/characteristic_architecture/p1_verification_spec.md`` row A5, whose
break set the tangency rows realise).

The object: the directions at the representative point ``x = c e_x`` modulo
the point's stabiliser ``Stab(x)`` in the chart's group, as a box in at most
two measure-uniform coordinates with ``dOmega = density * d(coordinates)``.

Independence (X4). The density is DERIVED in production as ``4 pi`` over the
box's coordinate measure, so the total-measure row is a smoke row and proves
nothing about the coordinates (the main agent's note, 2026-10-06). The gate
is row (ii): an integrand that is a function of the stabiliser's invariants,
integrated over the WHOLE sphere by a product rule in the polar angle about
``e_z`` and the azimuth (written here, sharing nothing with the domain's
coordinates), against the density times its integral over the domain through
the domain's own ``direction``. Every other expectation is typed by hand or
computed here from the definition: the stabiliser's element matrices, the
invariants, the closed-form tangencies; the kernel's chord is the cross-route
for the impact parameter and the tangencies.
"""
from __future__ import annotations

import math

import mpmath as mp
import numpy as np
import pytest

from orpheus.geometry.chart import Chart, RadialImage
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.line import Line
from orpheus.geometry.transformation import RigidMotion
from orpheus.numerics.symmetry import SubgroupOfO3

_EPS = np.finfo(float).eps
_HERE = "tests/gates/geometry/test_chart_directions.py::"
_CHART = "tests/gates/geometry/test_chart.py::"
_CHORD = "tests/gates/geometry/test_chord.py::"

_SLAB, _CYL, _SPH = (Chart(c) for c in (CoordSystem.CARTESIAN, CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL))

# (id, chart, orbit coordinate, axes, bounds, stabiliser, its name in this file): the domain table, by hand.
_TABLE = [
    ("sphere_r1.5", _SPH, 1.5, ("cosine",), ((-1.0, 1.0),), SubgroupOfO3.O2("x"), "O2x"),
    ("sphere_centre", _SPH, 0.0, (), (), SubgroupOfO3.O3, "O3"),
    ("cylinder_r1.5", _CYL, 1.5, ("angle", "axial_cosine"), ((0.0, math.pi), (0.0, 1.0)),
     SubgroupOfO3.Dnh(1), "D1h"),
    ("cylinder_axis", _CYL, 0.0, ("axial_cosine",), ((0.0, 1.0),), SubgroupOfO3.Dinfh, "Dinfh"),
    ("slab_r0.5", _SLAB, 0.5, ("cosine",), ((-1.0, 1.0),), SubgroupOfO3.O2("x"), "O2x"),
    ("slab_r-0.7", _SLAB, -0.7, ("cosine",), ((-1.0, 1.0),), SubgroupOfO3.O2("x"), "O2x"),
    ("slab_r0", _SLAB, 0.0, ("cosine",), ((-1.0, 1.0),), SubgroupOfO3.O2("x"), "O2x"),
]
_IDS = [t[0] for t in _TABLE]
_POINTS = [(t[1], t[2], t[6]) for t in _TABLE]


# ── the stabiliser's invariants and elements, written here ──────────────────


def _invariants(stabiliser_name: str, omega: np.ndarray) -> np.ndarray:
    """A complete set of invariants of each stabiliser on directions, by hand (``(..., k)``)."""
    match stabiliser_name:
        case "O2x":
            return omega[..., :1]
        case "D1h":
            return np.stack([omega[..., 0], np.abs(omega[..., 1]), np.abs(omega[..., 2])], axis=-1)
        case "Dinfh":
            return np.abs(omega[..., 2:3])
        case "O3":
            return omega[..., :0]
    raise KeyError(stabiliser_name)


def _coordinates_of(name: str, axes: tuple[str, ...], omega: np.ndarray) -> np.ndarray:
    """The domain coordinates of a direction from the axis definitions alone (never through the SUT)."""
    columns = []
    for axis in axes:
        match axis:
            case "cosine":
                columns.append(omega[..., 0])
            case "angle":
                columns.append(np.arctan2(np.abs(omega[..., 1]), omega[..., 0]))
            case "axial_cosine":
                columns.append(np.abs(omega[..., 2]))
            case _:
                raise KeyError(axis)
    return np.stack(columns, axis=-1) if columns else omega[..., :0]


def _matrix(motion: RigidMotion) -> np.ndarray:
    return np.asarray(motion.linear, dtype=float)


def _rot(axis, angle) -> np.ndarray:
    return _matrix(RigidMotion.rotation_about_axis(axis=axis, angle=angle))


def _mir(normal) -> np.ndarray:
    return _matrix(RigidMotion.reflection(normal=normal))


_ROOT2 = math.sqrt(2.0)
# The elements of each stabiliser, as matrices built here (continuous groups: a
# dense-generating rotation, a half turn and two mirrors; anti-pattern #13).
_ELEMENTS = {
    "O2x": [_rot((1.0, 0.0, 0.0), _ROOT2), _rot((1.0, 0.0, 0.0), math.pi), _mir((0.0, 1.0, 0.0)),
            _mir((0.0, 0.0, 1.0)), _mir((0.0, math.cos(0.3), math.sin(0.3)))],
    "D1h": [np.diag([1.0, -1.0, 1.0]), np.diag([1.0, 1.0, -1.0]), np.diag([1.0, -1.0, -1.0])],
    "Dinfh": [_rot((0.0, 0.0, 1.0), _ROOT2), _rot((0.0, 0.0, 1.0), math.pi), _mir((0.0, 0.0, 1.0)),
              _mir((math.cos(0.3), math.sin(0.3), 0.0)), _rot((1.0, 0.0, 0.0), math.pi)],
    "O3": [_rot((1.0, 2.0, -0.5), _ROOT2), -np.eye(3), _mir((0.3, -1.0, 0.2))],
}
# Candidates for the stabiliser-by-membership row: inside and outside every group above.
_CANDIDATES = [
    np.eye(3), _rot((1.0, 0.0, 0.0), _ROOT2), _rot((1.0, 0.0, 0.0), math.pi), _rot((0.0, 0.0, 1.0), _ROOT2),
    _rot((0.0, 0.0, 1.0), math.pi), _rot((0.0, 1.0, 0.0), math.pi), _rot((1.0, 2.0, -0.5), _ROOT2),
    _mir((1.0, 0.0, 0.0)), _mir((0.0, 1.0, 0.0)), _mir((0.0, 0.0, 1.0)), _mir((0.0, 1.0, 1.0)),
    _mir((1.0, 1.0, 0.0)), -np.eye(3),
]


# ── the table and the smoke row ─────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.rests_on(_CHART + "test_the_chart_is_derived_from_its_kept_columns_and_its_group",
                      _CHART + "test_the_singular_strata_and_their_isotropy")
@pytest.mark.parametrize(("chart", "r", "axes", "bounds", "stabiliser"), [t[1:6] for t in _TABLE], ids=_IDS)
def test_the_direction_domain_of_each_chart_and_point(chart: Chart, r: float, axes, bounds, stabiliser) -> None:
    """The domain's axes, bounds, stabiliser and point per chart, on and off the stratum, by hand.

    The smoke leg: ``density * prod(widths) == 4 pi`` (tautological by
    construction since the density is derived from the widths; kept to catch
    a broken box, never credited for the coordinates).
    First reds (run, ``[M]`` 2026-10-06, battery
    ``scratch/characteristic_architecture/p1_step_a/battery``): the cylinder off
    its axis given the sphere's domain (``cosine``, ``O2("x")``, arm D2);
    ``alpha`` on ``[0, pi/2]`` (D7); the density's ``4`` dropped (D1, the smoke
    leg, all 7 rows).
    """
    dom = chart.directions_at(r)
    assert dom.axes == axes
    assert tuple(tuple(float(v) for v in b) for b in dom.bounds) == bounds
    assert dom.stabiliser == stabiliser
    np.testing.assert_array_equal(dom.point, [r, 0.0, 0.0])
    measure = math.prod(hi - lo for lo, hi in dom.bounds)
    assert abs(dom.density * measure - 4.0 * math.pi) <= 4 * _EPS * 4.0 * math.pi


# ── (ii) the gate: an invariant integrand against a full-sphere rule ────────


def _integrand(name: str, omega: np.ndarray) -> np.ndarray:
    """A smooth function of the stabiliser's invariants, not invariant under any larger group of the table.

    D1h: odd in ``Omega_x`` and weighting ``Omega_y^2`` and ``Omega_z^2``
    differently (so neither ``O2("x")`` nor ``D2h`` leaves it invariant);
    O2x: a function of ``Omega_x`` alone, odd part included; Dinfh: of
    ``Omega_z^2``; O3: a constant. Smooth on the sphere (no ``|Omega_z|``
    kink, which a Gauss rule in ``cos(theta)`` integrates only to ``[M]`` ~1e-5).
    """
    x, y, z = omega[..., 0], omega[..., 1], omega[..., 2]
    match name:
        case "D1h":
            return np.exp(1.3 * x + 0.7 * z * z) * (1.0 + 0.5 * y * y)
        case "O2x":
            return np.exp(1.3 * x) * (1.0 + x * x)
        case "Dinfh":
            return np.exp(0.7 * z * z) * (1.0 + z ** 4)
        case "O3":
            return np.full(omega.shape[:-1], 2.5)
    raise KeyError(name)


def _full_sphere_integral(name: str) -> float:
    """The integrand over S^2 by Gauss-Legendre in cos(theta) about e_z times the trapezoid in phi."""
    t, wt = np.polynomial.legendre.leggauss(64)
    phi = 2.0 * math.pi * np.arange(128) / 128
    s = np.sqrt((1.0 - t) * (1.0 + t))
    omega = np.stack([s[:, None] * np.cos(phi), s[:, None] * np.sin(phi), np.broadcast_to(t[:, None], (64, 128))], -1)
    return float(np.sum(wt[:, None] * _integrand(name, omega)) * (2.0 * math.pi / 128))


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_direction_domain_of_each_chart_and_point")
@pytest.mark.parametrize(("chart", "r", "name"), _POINTS, ids=_IDS)
def test_an_invariant_integrand_over_the_domain_is_its_integral_over_the_sphere(chart: Chart, r: float, name: str) -> None:
    """(ii): ``int_{S^2} f dOmega == density * int_domain f(direction(q)) dq`` to 1e-12, f invariant under Stab(x).

    The reference rule (polar axis ``e_z``, ``[M]`` exact to ~1e-15 for these
    smooth integrands) shares no coordinate with the domain. The integrand is
    chosen per the HAND-TYPED stabiliser of the table row, not the SUT's.
    First reds (run): the density dropped (``pi`` for ``4 pi`` over the box: a
    factor 4 everywhere; D1, 7 of 7 rows); the cylinder taken as having stabiliser ``O2("x")`` (D2)
    (its domain the cosine to ``x``, which integrates the ``D_1h`` integrand
    over the meridian ``Omega_z = 0`` only); ``alpha`` on ``[0, pi/2]`` (half the
    orbit, the density doubled: the odd ``Omega_x`` part reads the error; D7).
    """
    dom = chart.directions_at(r)
    full = _full_sphere_integral(name)
    weight_sum = _full_sphere_integral("O3") / 2.5
    assert abs(weight_sum - 4.0 * math.pi) <= 1e-14 * 4.0 * math.pi, "the reference rule's own total measure"
    domain = _domain_integral(dom, name)
    assert abs(domain - full) <= 1e-12 * abs(full), f"domain {domain!r} vs sphere {full!r}"


def _domain_integral(dom, name: str) -> float:
    """The integrand (named by the TABLE's stabiliser, never the SUT's) over the domain through its ``direction``, times the density.

    Each cosine axis is integrated in its angle (``c = cos(theta)``, weight
    ``sin(theta)``), each angle axis directly, both by 48-point Gauss-Legendre:
    a cosine's sine, ``sqrt(1 - c^2)``, has a square-root endpoint a Gauss rule
    in ``c`` resolves only algebraically.
    """
    x, w = np.polynomial.legendre.leggauss(48)
    nodes, weights = [], []
    for axis, (lo, hi) in zip(dom.axes, dom.bounds, strict=True):
        if axis.endswith("_cosine"):
            a, b = math.acos(hi), math.acos(lo)
            theta = 0.5 * (b - a) * x + 0.5 * (a + b)
            nodes.append(np.cos(theta))
            weights.append(0.5 * (b - a) * w * np.sin(theta))
        else:
            nodes.append(0.5 * (hi - lo) * x + 0.5 * (hi + lo))
            weights.append(0.5 * (hi - lo) * w)
    if not nodes:
        return float(dom.density * _integrand(name, dom.direction(np.empty((1, 0))))[0])
    grid = np.stack(np.meshgrid(*nodes, indexing="ij"), axis=-1)
    weight = np.ones(grid.shape[:-1])
    for i, wi in enumerate(weights):
        shape = [1] * len(weights)
        shape[i] = -1
        weight = weight * wi.reshape(shape)
    return float(dom.density * np.sum(weight * _integrand(name, dom.direction(grid))))


# ── (iii) the domain is a fundamental domain of the stabiliser ─────────────


def _seeded_directions(n: int, seed: int) -> np.ndarray:
    om = np.random.default_rng(seed).normal(size=(n, 3))
    return om / np.linalg.norm(om, axis=1, keepdims=True)


def _interior_samples(dom, n: int, seed: int) -> np.ndarray:
    """Seeded coordinates strictly inside the box (the boundary is a measure-zero set of strata)."""
    lo = np.array([b[0] for b in dom.bounds])
    hi = np.array([b[1] for b in dom.bounds])
    u = np.random.default_rng(seed).uniform(0.01, 0.99, (n, len(lo)))
    return lo + (hi - lo) * u


def _in_image(dom, axes, omega: np.ndarray) -> np.ndarray:
    """Whether each direction is a domain representative: its own coordinates' representative is itself."""
    rep = dom.direction(_coordinates_of("", axes, omega))
    return np.all(np.abs(rep - omega) <= 8 * _EPS, axis=-1)


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_direction_domain_of_each_chart_and_point")
@pytest.mark.parametrize(("chart", "r", "name"), _POINTS, ids=_IDS)
def test_the_domain_meets_every_orbit_of_the_stabiliser_once(chart: Chart, r: float, name: str) -> None:
    """(iii): the stabiliser's elements (matrices built here) preserve the invariants; no element carries one representative to another; every direction's orbit meets the domain.

    Legs: (a) every element built here is in ``dom.stabiliser`` (the element
    set is the SUT's group, so the legs below speak about it); (b) each
    element maps 300 interior representatives to directions with the same
    invariants; (c) an element maps an interior representative either to
    itself (a direction the element fixes, as ``sigma_z`` fixes the sphere's
    half-meridian) or OUTSIDE the domain's image; (d) 2000 isotropic directions
    have coordinates inside the bounds whose representative has their
    invariants (the orbit meets the domain); (e) on the finite ``D_1h``, exactly
    one of the four images of a generic direction is its representative.
    First reds (run): ``alpha`` on ``[0, pi/2]`` (D7, leg d: ``Omega_x < 0``
    directions fall outside); the cylinder given ``O2("x")`` (D2, leg a: the
    hand-built ``D_1h`` elements are in ``O2("x")``, so a later leg reads it).
    Not run: a representative taking ``Omega_y`` signed (leg e).
    """
    dom = chart.directions_at(r)
    axes = dom.axes
    elements = _ELEMENTS[name]
    for g in elements:
        assert dom.stabiliser.realization.contains_element(RigidMotion(g, np.zeros(3))), (
            f"an element of {name} built here is not in the domain's stabiliser {dom.stabiliser!r}"
        )
    rep = dom.direction(_interior_samples(dom, 300, 11)) if axes else dom.direction(np.empty((1, 0)))
    for g in elements:
        moved = rep @ g.T
        np.testing.assert_allclose(_invariants(name, moved), _invariants(name, rep), rtol=0, atol=8 * _EPS)
        fixed = np.all(np.abs(moved - rep) <= 8 * _EPS, axis=-1)
        assert np.all(fixed | ~_in_image(dom, axes, moved)), "an element carries a representative to another"
    omega = _seeded_directions(2000, 12)
    q = _coordinates_of(name, axes, omega)
    for j, (lo, hi) in enumerate(dom.bounds):
        assert np.all((q[:, j] >= lo) & (q[:, j] <= hi)), f"a direction's coordinate {axes[j]} leaves its bounds"
    np.testing.assert_allclose(_invariants(name, dom.direction(q)), _invariants(name, omega), rtol=0, atol=16 * _EPS)
    if name == "D1h":
        images = np.stack([omega] + [omega @ g.T for g in elements], axis=0)
        match = np.all(np.abs(images - dom.direction(q)[None]) <= 16 * _EPS, axis=-1).sum(axis=0)
        np.testing.assert_array_equal(match, 1)


# ── (iv) unit vectors, the coordinates round-trip, the sine is conditioned ─


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_direction_domain_of_each_chart_and_point")
@pytest.mark.parametrize(("chart", "r", "name"), [p for p in _POINTS if p[2] != "O3"], ids=[i for i, p in zip(_IDS, _POINTS) if p[2] != "O3"])
def test_the_representative_is_a_unit_vector_whose_coordinates_round_trip(chart: Chart, r: float, name: str) -> None:
    """(iv): ``|direction(q)| = 1`` to 2 eps and its coordinates are ``q`` (cosines to 2 eps, the angle to 2 eps pi).

    Samples: 3000 seeded points of the box with the cosines pushed to
    ``1 - 10^-k``, ``k`` up to 15, and the box's corners; the kept angle is not
    round-tripped at ``w = 1`` (the axial direction, a stratum: every angle
    names it). ``[M]`` 2026-10-06: 1 eps for the norm, 0 for the cosines,
    1 eps pi for the angle.
    First red (run): the sine replaced by the cosine (D0, the battery's
    positive control, 6 of 6 rows).
    """
    dom = chart.directions_at(r)
    lo = np.array([b[0] for b in dom.bounds])
    hi = np.array([b[1] for b in dom.bounds])
    rng = np.random.default_rng(13)
    q = lo + (hi - lo) * rng.uniform(size=(3000, len(lo)))
    cosine_axis = [j for j, a in enumerate(dom.axes) if a.endswith("_cosine")]
    for j in cosine_axis:
        q[:100, j] = hi[j] - 10.0 ** -rng.uniform(1, 15, 100)
    corners = np.stack(np.meshgrid(*[np.array(b) for b in dom.bounds], indexing="ij"), axis=-1).reshape(-1, len(lo))
    q = np.concatenate([q, corners])
    om = dom.direction(q)
    assert np.max(np.abs(np.linalg.norm(om, axis=-1) - 1.0)) <= 2 * _EPS
    back = _coordinates_of(name, dom.axes, om)
    for j, axis in enumerate(dom.axes):
        if axis == "angle":
            keep = q[:, dom.axes.index("axial_cosine")] < 1.0
            assert np.max(np.abs(back[keep, j] - q[keep, j])) <= 2 * _EPS * math.pi
        else:
            assert np.max(np.abs(back[:, j] - q[:, j])) <= 2 * _EPS


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_representative_is_a_unit_vector_whose_coordinates_round_trip")
@pytest.mark.parametrize(("chart", "r", "column"), [(_SPH, 1.5, 1), (_SLAB, 0.5, 1), (_CYL, 1.5, 1), (_CYL, 0.0, 0)],
                         ids=["sphere", "slab", "cylinder", "cylinder_axis"])
def test_the_representative_sine_is_conditioned_near_the_pole(chart: Chart, r: float, column: int) -> None:
    """The transverse component ``sqrt(1 - c^2)`` to 2 ulp of its 50-digit value, at cosines ``c`` near 1.

    ``sqrt(1 - c*c)`` loses up to ``[M]`` ~1e-4 relative at ``c = 1 - 1e-12``;
    ``sqrt((1 - c)(1 + c))`` keeps it to an ulp. The cylinder's row reads
    ``Omega_y`` at ``alpha = pi/2`` (``sin`` of it is 1.0 exactly in binary64).
    A near-pole representative feeds every grazing chord, so a lost digit
    here is a wrong ``b`` there; the unit-norm row cannot see it (the squares
    absorb an ``eps``-sized error in ``s^2``).
    First red (run): the sine as ``sqrt(1 - c*c)`` (D4, 4 of 4 rows; no other
    row of the two files reds on it).
    """
    mp.mp.dps = 50
    dom = chart.directions_at(r)
    cosines = [1.0 - 1e-12, -(1.0 - 2.0 ** -30), 0.3, 1.0 - 2.0 ** -52, 0.999999]
    for c in cosines:
        if dom.axes == ("angle", "axial_cosine"):
            q = np.array([math.pi / 2, abs(c)])
            c = abs(c)
        elif dom.axes[0].endswith("_cosine") and dom.bounds[0][0] == 0.0:
            q, c = np.array([abs(c)]), abs(c)
        else:
            q = np.array([c])
        got = float(dom.direction(q)[column])
        exact = mp.sqrt((1 - mp.mpf(c)) * (1 + mp.mpf(c)))
        assert abs(mp.mpf(got) - exact) <= 2 * mp.mpf(np.spacing(float(exact))), f"c = {c!r}: {got!r} vs {exact}"


# ── (v) the impact parameter and the tangencies, against the kernel ────────

_RADIAL = [(_SPH, 1.5), (_SPH, 0.37), (_SPH, 0.0), (_CYL, 1.5), (_CYL, 2.0), (_CYL, 0.0)]


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_representative_is_a_unit_vector_whose_coordinates_round_trip",
                      _CHART + "test_the_image_of_a_line_and_its_shift")
@pytest.mark.parametrize(("chart", "r"), _RADIAL, ids=[f"{c.coord.name.lower()}_r{r}" for c, r in _RADIAL])
def test_the_impact_parameter_is_the_kernels_and_the_closed_form(chart: Chart, r: float) -> None:
    """(v): ``impact_parameter(q)`` equals ``Chart.image(Line.through(point, direction(q))).impact_parameter``
    to 8 eps r, and the closed form ``r sin(alpha)`` / ``r sqrt(1 - mu^2)`` evaluated here.

    Samples include 200 near-axial cylinder directions (``w = 1 - 10^-k``),
    where ``|P Omega|`` read as ``sqrt(1 - Omega_z^2)`` loses ``[M]`` ~1e-4.
    ``[M]`` 2026-10-06: the two routes agree to 3.3 eps r.
    First red (run): ``|P Omega|`` as ``sqrt(1 - Omega_z**2)`` in
    ``impact_parameter`` (D3: the two cylinder rows off the axis).
    """
    dom = chart.directions_at(r)
    lo = np.array([b[0] for b in dom.bounds])
    hi = np.array([b[1] for b in dom.bounds])
    rng = np.random.default_rng(14)
    q = lo + (hi - lo) * rng.uniform(size=(2000, len(lo)))
    if "axial_cosine" in dom.axes:
        q[:200, dom.axes.index("axial_cosine")] = 1.0 - 10.0 ** -rng.uniform(4, 15, 200)
    om = dom.direction(q) if dom.axes else dom.direction(np.empty((1, 0)))
    b = dom.impact_parameter(q if dom.axes else np.empty((1, 0)))
    image = chart.image(Line.through(np.broadcast_to(dom.point, om.shape), om))
    if not isinstance(image, RadialImage):
        pytest.fail("a radial chart's line image is radial")
    kernel = image.impact_parameter
    assert np.max(np.abs(b - kernel)) <= 8 * _EPS * max(r, 1.0)
    match dom.axes:
        case ("cosine",):
            closed = r * np.sqrt((1.0 - q[:, 0]) * (1.0 + q[:, 0]))
        case ("angle", "axial_cosine"):
            closed = r * np.sin(q[:, 0])
        case _:
            closed = np.zeros(b.shape)
    assert np.max(np.abs(b - closed)) <= 8 * _EPS * max(r, 1.0)


@pytest.mark.foundation
def test_a_slab_line_has_no_impact_parameter() -> None:
    """The slab's refusal (ruled 2026-10-06): its group fixes the kept space, so a line has no impact parameter.

    First red (run): the refusal removed (D10: numpy's empty-reduction error
    raises instead, and the message key reds).
    """
    with pytest.raises(ValueError, match="no impact parameter"):
        _SLAB.directions_at(0.5).impact_parameter(np.array([0.3]))


def _closed_tangencies(axes, r: float, level: float, radial: bool) -> list:
    """The break set by hand, in mpmath: ``mu = +-sqrt(1 - l^2/r^2)`` or ``alpha = asin(l/r), pi - asin(l/r)``."""
    mp.mp.dps = 50
    if not radial or not axes or axes[0] not in ("cosine", "angle") or not 0.0 < level <= r:
        return []
    ratio = mp.mpf(level) / mp.mpf(r)
    if axes[0] == "cosine":
        m = mp.sqrt(1 - ratio ** 2)
        return [-m, m] if m != 0 else [mp.mpf(0)]
    a = mp.asin(ratio)
    return [a, mp.pi - a] if a != mp.pi / 2 else [a]


_BREAK_CASES = [  # (chart, breakpoints, r): the point inside a region, on an interface, on the outer wall
    (_SPH, (0.0, 0.3, 1.1, 2.0), 1.5), (_SPH, (0.0, 0.3, 1.1, 2.0), 1.1), (_SPH, (0.0, 0.3, 1.1, 2.0), 2.0),
    (_SPH, (0.4, 1.1, 2.0), 0.7),
    (_CYL, (0.0, 0.3, 1.1, 2.0), 1.5), (_CYL, (0.0, 0.3, 1.1, 2.0), 1.1), (_CYL, (0.0, 0.3, 1.1, 2.0), 2.0),
    (_CYL, (0.4, 1.1, 2.0), 0.7),
]
_BREAK_IDS = [f"{c.coord.name.lower()}_{len(bp) - 1}reg{'_hollow' if bp[0] else ''}_r{r}" for c, bp, r in _BREAK_CASES]


@pytest.mark.foundation
@pytest.mark.parametrize(("chart", "breakpoints", "r"), _BREAK_CASES, ids=_BREAK_IDS)
def test_the_tangencies_are_the_closed_form_break_set(chart: Chart, breakpoints, r: float) -> None:
    """A5: ``tangencies(r_k)`` is the closed form, sorted, to 4 ulp of the axis scale; empty where none.

    Two values for ``0 < r_k < r``; one, the grazing direction (``mu = 0``,
    ``alpha = pi/2``), at ``r_k = r`` (ruled 2026-10-06: the point's own level,
    where the backward ray switches sides of the surface); none for
    ``r_k > r``, for ``r_k = 0`` (a solid body's centre is a stratum), for a
    negative level.
    First reds (run): ``asin(r_k/r)`` without its reflection ``pi - asin``
    (D5); the grazing value dropped at ``r_k = r``, the brief's first draft
    (D6); ``acos`` for ``asin`` (D9).
    """
    dom = chart.directions_at(r)
    scale = math.pi if dom.axes[0] == "angle" else 1.0
    for level in (*breakpoints, -0.3, 0.0):
        got = np.asarray(dom.tangencies(level), dtype=float)
        expected = _closed_tangencies(dom.axes, r, level, chart.acts_on_kept_space)
        assert got.shape == (len(expected),), f"level {level}: {got} vs {[mp.nstr(e, 17) for e in expected]}"
        assert np.all(np.diff(got) > 0), f"level {level}: not sorted"
        for g, e in zip(got, expected, strict=True):
            assert abs(mp.mpf(float(g)) - e) <= 4 * _EPS * scale, f"level {level}: {g!r} vs {mp.nstr(e, 20)}"


@pytest.mark.foundation
@pytest.mark.parametrize(("chart", "r"), [(_SLAB, 0.5), (_SPH, 0.0), (_CYL, 0.0)], ids=["slab", "sphere_centre", "cylinder_axis"])
def test_no_tangency_on_the_slab_or_at_a_stratum(chart: Chart, r: float) -> None:
    """The slab's only break (``u = 0``) is the reference's own; at a stratum every line has ``b = 0``.

    First red (not run): the slab's parallel direction returned as a
    tangency; the guard on the first axis is the line that refuses it.
    """
    dom = chart.directions_at(r)
    for level in (-0.7, 0.0, 0.3, 0.5, 1.1, 2.0):
        assert np.asarray(dom.tangencies(level)).size == 0


def _crossing_count(partition: ConcentricPartition, dom, q: np.ndarray, k: int) -> np.ndarray:
    om = dom.direction(q)
    ch = partition.chord(Line.through(np.broadcast_to(dom.point, om.shape), om))
    present = np.broadcast_to(ch.crossings.present, ch.crossings.parameter.shape)
    return np.sum(present & (np.broadcast_to(ch.crossings.breakpoint, present.shape) == k), axis=-1)


def _on_axis(dom, values: np.ndarray) -> np.ndarray:
    """Coordinates with the tangency axis at ``values`` and the axial cosine (if any) at 0.37."""
    if dom.axes == ("angle", "axial_cosine"):
        return np.stack([values, np.full(values.shape, 0.37)], axis=-1)
    return values[:, None]


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_tangencies_are_the_closed_form_break_set",
                      _CHORD + "test_a_tangency_is_not_a_crossing")
@pytest.mark.parametrize(("chart", "breakpoints", "r"), _BREAK_CASES, ids=_BREAK_IDS)
def test_the_tangencies_are_where_the_kernels_crossings_change(chart: Chart, breakpoints, r: float) -> None:
    """(v), the kernel cross-route (elegance E6): the chord's crossing set changes exactly at the tangencies.

    For every level ``r_k < r``: (a) bracketing each tangency at ``1e-9`` in the
    coordinate flips the number of present crossings of ``c = r_k`` between 0
    and 2 (``[M]`` the main agent, 2026-10-06: at one ulp both sides can read
    the same, since the kernel recomputes ``b`` from the line's foot);
    (b) on a 4001-node grid over the axis (offset by an irrational), the count
    changes between consecutive nodes exactly in the intervals holding a
    tangency, and every tangency's interval changes. For ``r_k = r`` (the point
    on the surface) the count is 2 on both sides of the grazing value and the
    crossing AT the point changes its sense there (outgoing for
    ``Omega . x > 0``).
    First reds (run): a missing reflection (D5: a change interval with no
    tangency); ``acos`` for ``asin`` (D9); the grazing value dropped (D6);
    ``alpha`` on ``[0, pi/2]`` (D7: the grid misses the reflected change).
    """
    partition = ConcentricPartition(chart, breakpoints)
    dom = chart.directions_at(r)
    lo, hi = dom.bounds[0]
    grid = lo + (hi - lo) * (np.arange(4001) + 1.0 / math.sqrt(5.0)) / 4001.5
    checked = 0
    for k, level in enumerate(breakpoints):
        taus = np.asarray(dom.tangencies(level), dtype=float)
        if level <= 0.0 or level > r:
            continue
        if level < r:
            for tau in taus:
                below, above = _crossing_count(partition, dom, _on_axis(dom, np.array([tau - 1e-9, tau + 1e-9])), k)
                assert {int(below), int(above)} == {0, 2}, f"level {level}: tangency {tau!r} reads {below}, {above}"
                checked += 1
            counts = _crossing_count(partition, dom, _on_axis(dom, grid), k)
            changed = np.flatnonzero(np.diff(counts) != 0)
            holds = [int(np.searchsorted(grid, tau) - 1) for tau in taus]
            assert sorted(holds) == list(changed), f"level {level}: change intervals {changed} vs tangencies {holds}"
        else:
            (graze,) = taus
            # b - r is SECOND order in the offset at the point's own level ([M] 2026-10-06: at 1e-9 the
            # chord's b rounds to r and reads a tangency on both sides), so the sides are taken at 1e-3.
            for side in (graze - 1e-3, graze + 1e-3):
                om = dom.direction(_on_axis(dom, np.array([side])))[0]
                line = Line.through(dom.point, om)
                c = partition.chord(line).crossings
                on_k = np.asarray(c.present) & (np.asarray(c.breakpoint) == k)
                assert int(on_k.sum()) == 2, f"the point's own level: {int(on_k.sum())} crossings beside the grazing value"
                t0 = float(line.parameter_of(dom.point))
                nearest = np.flatnonzero(on_k)[np.argmin(np.abs(np.asarray(c.parameter)[on_k] - t0))]
                outgoing = float(om[0]) > 0.0
                assert int(np.asarray(c.sense)[nearest]) == (1 if outgoing else -1), "the crossing at the point"
            checked += 1
    assert checked >= 2, f"only {checked} tangencies exercised"


# ── (vi) the stabiliser, by membership; (vii) refusals ─────────────────────


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_direction_domain_of_each_chart_and_point",
                      _CHART + "test_membership_in_the_symmetry_group_is_invariance_of_the_orbit_coordinate")
@pytest.mark.parametrize(("chart", "r"), [p[:2] for p in _POINTS], ids=_IDS)
def test_the_stabiliser_is_the_isotropy_of_the_point(chart: Chart, r: float) -> None:
    """(vi): ``Q in dom.stabiliser`` iff ``Q in G_c`` and ``Q x = x``, over 13 candidates in and out of every group.

    The right side is evaluated here: the chart's membership (gated in
    ``test_chart.py``) and the action on the point itself. The candidates
    include a rotation about ``x`` by ``sqrt(2)`` (in ``O(2)_x``, not in
    ``D_1h``), the half turn about ``x`` (in both), the mirror ``sigma_x`` and
    the half turn about ``z`` (in ``D_inf_h``, moving a point off the axis).
    First red (run): the cylinder's stabiliser taken as ``O2("x")`` (D2: the
    ``sqrt(2)`` rotation about ``x`` is accepted, though it tilts the axis).
    """
    dom = chart.directions_at(r)
    x = np.array([r, 0.0, 0.0])
    for Q in _CANDIDATES:
        motion = RigidMotion(Q, np.zeros(3))
        expected = chart.contains(motion) and bool(np.allclose(Q @ x, x, rtol=0, atol=1e-12))
        assert dom.stabiliser.realization.contains_element(motion) is expected, (
            f"{dom.stabiliser!r} at c = {r}: membership of\n{Q}\nis not {expected}"
        )


@pytest.mark.foundation
@pytest.mark.parametrize("chart", [_SPH, _CYL], ids=["sphere", "cylinder"])
def test_a_point_off_the_charts_range_is_refused(chart: Chart) -> None:
    """(vii): a negative or non-finite orbit coordinate on a radial chart is refused; a slab accepts ``c < 0``.

    The negative leg keys on ``non-negative``, the non-finite leg on ``finite``;
    the positive leg is the table's ``slab_r-0.7`` row.
    First red (run): the radial chart accepting ``c = -0.1`` (D8, both
    charts).
    """
    with pytest.raises(ValueError, match="non-negative"):
        chart.directions_at(-0.1)
    for bad in (math.nan, math.inf):
        with pytest.raises(ValueError, match="finite"):
            chart.directions_at(bad)
    _SLAB.directions_at(-0.1)


# ── qa review rows (2026-10-06) ─────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.parametrize("r", [0.7, 1.5])
def test_a_line_along_the_cylinder_axis_has_the_points_impact_parameter(r: float) -> None:
    """At ``w = 1`` the line is parallel to the orbit space and keeps ``c = r``: its impact parameter is ``r``.

    First red (qa, run): the parallel fallback replaced by ``0.0`` left every
    other row green, because the samples stopped at ``w = 1 - 1e-15``.
    """
    dom = _CYL.directions_at(r)
    b = dom.impact_parameter(np.array([[0.0, 1.0], [1.3, 1.0], [math.pi, 1.0]]))
    np.testing.assert_array_equal(b, [r, r, r])


@pytest.mark.foundation
@pytest.mark.parametrize(("chart", "r", "bad"), [
    (_SPH, 1.5, [1.0 + 1e-12]), (_SLAB, 0.5, [-2.0]), (_CYL, 1.5, [0.3, 1.5]),
    (_CYL, 1.5, [-1e-3, 0.4]), (_CYL, 0.0, [math.nan]), (_SPH, 1.5, [math.inf]),
], ids=["sphere_cosine_above", "slab_cosine_below", "cyl_w_above", "cyl_angle_below", "axis_nan", "sphere_inf"])
def test_a_coordinate_outside_the_box_is_refused(chart: Chart, r: float, bad: list) -> None:
    """``direction`` returns unit vectors only: a coordinate outside its interval, or not finite, is refused.

    First red (qa, run): ``q = 2`` returned ``(2, 0, 0)`` and ``q = nan`` a NaN
    direction, whose impact parameter then read ``r`` (``NaN > 0`` is False).
    """
    with pytest.raises(ValueError, match="coordinate lies in"):
        chart.directions_at(r).direction(np.array([bad]))


@pytest.mark.foundation
@pytest.mark.parametrize("scale", [1e-200, 1e-160, 1.0, 1e160, 1e300])
def test_the_tangencies_are_scale_free(scale: float) -> None:
    """The break set depends on ``level / r`` only: the same at every scale, to 2 ulp.

    First red (qa, run): ``sqrt((r - l)(r + l)) / r`` gave ``[-0.0]`` at
    ``r = 1e-200`` and ``+-inf`` at ``r >= 1e160``.
    """
    for chart in (_SPH, _CYL):
        unit = chart.directions_at(2.0).tangencies(1.0)
        scaled = chart.directions_at(2.0 * scale).tangencies(1.0 * scale)
        np.testing.assert_allclose(scaled, unit, rtol=2 * _EPS, atol=0.0)
