"""Gates for :class:`orpheus.geometry.chart.Chart`, the orbit space of a 1-D coordinate system's symmetry group.

The geometric kernel seed (``.claude/plans/characteristic_reference_architecture.md``,
P0; verification spec ``scratch/characteristic_architecture/seed_verification_spec.md``
§5 row R1 and §7b). Expected values are typed by hand from the definitions;
none is derived from the chart.

Laws gated here:

- the orbit coordinate per chart (``x_0`` signed; the distance to the axis;
  the distance to the centre), on every coordinate axis;
- ``contains(g)`` decides ``c o g = c``: the closed form is checked against
  the orbit coordinate itself, evaluated on points, over motions inside and
  outside each group (an order/predicate compatibility law, anti-pattern #15);
- the singular strata: each stratum's isotropy group, realised, fixes the
  stratum and lies in ``G_c``;
- the projected speed ``|P Omega|``, summed from components, accurate near the
  axis where ``1 - Omega_z^2`` loses every digit;
- the beam density per chart (hand values; the integrals are
  ``test_line_measure.py``).
"""
from __future__ import annotations

import mpmath as mp
import numpy as np
import pytest

from orpheus.geometry.chart import Chart
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.transformation import RigidMotion
from orpheus.numerics.symmetry import SubgroupOfO3

_EPS = np.finfo(float).eps
_SLAB, _CYL, _SPH = (Chart(c) for c in (CoordSystem.CARTESIAN, CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL))
_HERE = "tests/gates/geometry/test_chart.py::"


@pytest.mark.l0
@pytest.mark.verifies("geometry-radial-coordinate")
def test_the_orbit_coordinate_of_each_chart() -> None:
    """R1: ``c(x)`` per chart, with every coordinate axis carrying a nonzero component in some row.

    Sphere: ``|x|`` on each axis (``(0.7, 0, 0)``, ``(0, 0.7, 0)``,
    ``(0, 0, -0.7)`` all give 0.7, exactly: ``sqrt(r*r) == r`` in binary64).
    Cylinder: the distance to the z axis, blind to z (``z`` in {0, -3.7, 1e3}).
    Slab: ``x_0``, signed (``(-0.5, 9, 9)`` gives -0.5).
    First reds: the cylinder read as a sphere (z enters); the sphere reading
    ``x_0`` only; the slab taking ``|x_0|``.
    """
    for p in ((0.7, 0.0, 0.0), (0.0, 0.7, 0.0), (0.0, 0.0, -0.7)):
        assert _SPH.orbit_coordinate(np.array(p)) == 0.7
    for z in (0.0, -3.7, 1e3):
        assert _CYL.orbit_coordinate(np.array([0.7, 0.0, z])) == 0.7
        assert _CYL.orbit_coordinate(np.array([0.0, -0.7, z])) == 0.7
    assert _CYL.orbit_coordinate(np.array([0.6, 0.8, 5.0])) == pytest.approx(1.0, rel=2 * _EPS)
    assert _SLAB.orbit_coordinate(np.array([-0.5, 9.0, 9.0])) == -0.5
    assert _SPH.orbit_coordinate(np.array([0.6, 0.0, 0.8])) == pytest.approx(1.0, rel=2 * _EPS)


def _rotation(axis: tuple[float, float, float], angle: float, shift=(0.0, 0.0, 0.0)) -> RigidMotion:
    return RigidMotion(RigidMotion.rotation_about_axis(axis=axis, angle=angle).linear, np.array(shift, dtype=float))


def _mirror(normal: tuple[float, float, float], shift=(0.0, 0.0, 0.0)) -> RigidMotion:
    return RigidMotion(RigidMotion.reflection(normal=normal).linear, np.array(shift, dtype=float))


_ROOT2 = float(np.sqrt(2.0))
# (chart, motion, in G_c?). sqrt(2) rad generates a dense subgroup of the rotations
# about its axis (anti-pattern #13: four right angles generate C_4 only).
_MEMBERSHIP = [
    (_SPH, _rotation((1.0, 2.0, -0.5), _ROOT2), True),
    (_SPH, _mirror((0.3, -1.0, 0.2)), True),
    (_SPH, RigidMotion.inversion(3), True),
    (_SPH, _rotation((0.0, 0.0, 1.0), _ROOT2, (0.3, 0.0, 0.0)), False),
    (_SPH, RigidMotion.translation_by([0.0, 0.0, 1e-3]), False),
    (_CYL, _rotation((0.0, 0.0, 1.0), _ROOT2, (0.0, 0.0, 3.7)), True),
    (_CYL, _mirror((1.0, 2.0, 0.0)), True),
    (_CYL, _mirror((0.0, 0.0, 1.0), (0.0, 0.0, 2.0)), True),
    (_CYL, _rotation((1.0, 0.0, 0.0), np.pi), True),               # flips the axis end for end
    (_CYL, _rotation((1.0, 0.0, 0.0), 0.4), False),                # tilts the axis
    (_CYL, RigidMotion.translation_by([0.0, 0.2, 0.0]), False),    # moves the axis
    (_SLAB, _rotation((1.0, 0.0, 0.0), _ROOT2, (0.0, 1.3, -0.4)), True),
    (_SLAB, _mirror((0.0, 1.0, 1.0)), True),
    (_SLAB, _mirror((1.0, 0.0, 0.0)), False),                      # x -> -x
    (_SLAB, _rotation((0.0, 0.0, 1.0), 0.4), False),               # tilts the normal
    (_SLAB, RigidMotion.translation_by([0.25, 0.0, 0.0]), False),
]


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_orbit_coordinate_of_each_chart")
@pytest.mark.parametrize(
    ("chart", "motion", "member"), _MEMBERSHIP,
    ids=[f"{c.coord.name.lower()}_{i}_{'in' if m else 'out'}" for i, (c, _, m) in enumerate(_MEMBERSHIP)],
)
def test_membership_in_the_symmetry_group_is_invariance_of_the_orbit_coordinate(
    chart: Chart, motion: RigidMotion, member: bool,
) -> None:
    """``contains(g)`` (closed form in ``(Q, t)``) agrees with ``c o g = c`` (evaluated on points).

    Two derivations of one predicate, sharing no code: the closed form reads
    the motion's matrix and translation; the evaluation reads the orbit
    coordinate of 200 seeded points before and after the motion. A member
    moves no orbit coordinate by more than ``16 eps (|x| + |t|)``; a
    non-member moves some coordinate by more than ``1e-4``.
    First reds: the cylinder's axis tested against ``e_x``; the sphere ignoring
    the translation; the slab accepting ``x -> -x`` (its ``Q e_x = -e_x``).
    """
    x = np.random.default_rng(20261005).uniform(-3, 3, (200, 3))
    moved = chart.orbit_coordinate(motion.on_points(x))
    still = chart.orbit_coordinate(x)
    drift = np.abs(moved - still)
    assert chart.contains(motion) is member
    if member:
        scale = np.linalg.norm(x, axis=1) + np.linalg.norm(motion.translation) + 1.0
        assert np.all(drift <= 16 * _EPS * scale)
    else:
        assert np.max(drift) > 1e-4


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_membership_in_the_symmetry_group_is_invariance_of_the_orbit_coordinate")
def test_the_singular_strata_and_their_isotropy() -> None:
    """The sphere's centre (``O(3)``) and the cylinder's axis (``D_inf_h``); the slab has none.

    Each stratum's isotropy, realised, is checked against the chart: a motion
    the isotropy group contains fixes the stratum's point (the origin) and
    lies in ``G_c``; a motion it refuses is refused by the chart too when it
    moves the stratum. First red: the cylinder's stratum declared with
    ``O(3)`` (a tilt of the axis would then be "isotropy" but not in ``G_c``).
    """
    assert _SLAB.singular_strata == ()
    (centre,) = _SPH.singular_strata
    (axis,) = _CYL.singular_strata
    assert centre.orbit_value == 0.0 and axis.orbit_value == 0.0
    assert centre.isotropy == SubgroupOfO3.O3
    assert axis.isotropy == SubgroupOfO3.Dinfh
    candidates = (
        _rotation((0.0, 0.0, 1.0), _ROOT2), _mirror((1.0, 2.0, 0.0)), _mirror((0.0, 0.0, 1.0)),
        _rotation((1.0, 0.0, 0.0), 0.4), _rotation((1.0, 2.0, -0.5), _ROOT2),
    )
    in_dinfh = [axis.isotropy.realization.contains_element(g) for g in candidates]
    assert in_dinfh == [True, True, True, False, False]
    for g, inside in zip(candidates, in_dinfh):
        assert _CYL.contains(g) is inside
        assert _SPH.contains(g)              # every linear motion fixes the centre
        assert centre.isotropy.realization.contains_element(g)


@pytest.mark.l0
@pytest.mark.verifies("geometry-cylinder-axial-factor")
def test_the_projected_speed_is_summed_from_components() -> None:
    """``|P Omega|``: ``|Omega_x|`` (slab), ``sqrt(Omega_x^2 + Omega_y^2)`` (cylinder), 1 (sphere).

    Reference: mpmath at 50 dps on the stored float components. The near-axial
    direction ``(0, 1e-6, sqrt(1 - 1e-12))`` is the conditioning row: the sum of
    squares is within 2 ulp, while ``sqrt(1 - Omega_z^2)`` errs by ``[M]``
    8.9e-5 relative (spec §2 regime 5).
    First red: ``1 - Omega_z^2`` (the same mathematics, a different float spelling).
    """
    mp.mp.dps = 50
    rows = np.array([
        [0.0, 0.6, 0.8], [0.28, 0.0, -0.96], [0.0, 1.0, 0.0], [0.0, 0.8, -0.6],
        [0.0, 1e-6, float(np.sqrt(1.0 - 1e-12))],
    ])
    speed = _CYL.projected_speed(rows)
    for omega, got in zip(rows, speed):
        exact = mp.sqrt(mp.mpf(omega[0]) ** 2 + mp.mpf(omega[1]) ** 2)
        assert abs(mp.mpf(got) - exact) <= 2 * mp.mpf(np.spacing(float(exact)))
    np.testing.assert_array_equal(_SLAB.projected_speed(rows), np.abs(rows[:, 0]))
    np.testing.assert_array_equal(_SPH.projected_speed(rows), np.ones(len(rows)))


@pytest.mark.l0
@pytest.mark.verifies("geometry-measure-on-lines")
def test_the_beam_density_of_each_chart() -> None:
    """The parallel-beam density over ``b >= 0``: ``2 pi b`` (sphere), ``2 |P Omega|`` (cylinder, the lines at ``+-b``), ``|Omega_x|`` (slab).

    ``S_{d-2} b^{d-2} |P Omega|`` with ``S_0 = 2``, ``S_1 = 2 pi``. Hand values at
    ``b = 0.7`` and ``Omega = (0.28, 0.0, -0.96)``: sphere ``1.4 pi``; cylinder
    0.56; slab 0.28, independent of ``b``. The integrals
    these densities must reproduce (areas, volumes, Cauchy's mean chord) are the
    L1 rows of ``test_line_measure.py``.
    First reds: the cylinder's density without ``|P Omega|`` (the planar measure;
    Cauchy's cylinder mean chord then reads 2.467 R, not 2 R); without ``S_0 = 2``.
    """
    omega = np.array([0.28, 0.0, -0.96])
    b = np.array([0.0, 0.7, 1.9])
    np.testing.assert_allclose(_SPH.beam_density(b, omega), 2 * np.pi * b, rtol=2 * _EPS, atol=0)
    np.testing.assert_array_equal(_CYL.beam_density(b, omega), np.full(3, 0.56))
    np.testing.assert_array_equal(_SLAB.beam_density(b, omega), np.full(3, 0.28))


@pytest.mark.foundation
def test_the_chart_is_derived_from_its_kept_columns_and_its_group() -> None:
    """The pair the chart's verbs derive from (ruled 2026-10-05): kept columns ``d`` and the linear group ``L``.

    Slab: ``d = 1``, ``L = O(2)_x``, which fixes the kept space (no reflection of
    ``e_x``), so its coordinate is signed and it has no stratum. Cylinder:
    ``d = 2``, ``D_inf_h``; sphere: ``d = 3``, ``O(3)``; both act on the kept
    space as ``O(d)``. ``d`` is the exponent of the coordinate system's measure
    coordinate ``T = r^d``, checked against ``CoordSystem.measure`` on one shell.
    First red: ``acts_on_kept_space`` true on the slab (its coordinate becomes
    ``|x_0|``; ``test_the_orbit_coordinate_of_each_chart`` reds too).
    """
    rows = ((_SLAB, 1, SubgroupOfO3.O2("x"), False), (_CYL, 2, SubgroupOfO3.Dinfh, True), (_SPH, 3, SubgroupOfO3.O3, True))
    for chart, d, group, acts in rows:
        assert chart.kept_columns == d
        assert chart.group == group
        assert chart.acts_on_kept_space is acts
        assert (len(chart.singular_strata) == 1) is acts
    shell = np.array([0.5, 1.5])
    np.testing.assert_allclose(_CYL.measure(shell), np.pi * (1.5 ** 2 - 0.5 ** 2), rtol=2 * _EPS)
    np.testing.assert_allclose(_SPH.measure(shell), 4 / 3 * np.pi * (1.5 ** 3 - 0.5 ** 3), rtol=2 * _EPS)


@pytest.mark.l0
@pytest.mark.verifies("geometry-line-crossing-law")
@pytest.mark.rests_on(_HERE + "test_the_projected_speed_is_summed_from_components")
def test_the_image_of_a_line_and_its_shift() -> None:
    """``Chart.image``: the radial image's ``(b, t*, |P Omega|)`` and the axial image's ``(c_foot, rate)``, and ``shifted``.

    Hand rows: the cylinder line through ``(0.7, 0.4, 5)`` along ``(0, 0.6, 0.8)``
    has ``b = 0.7``, ``|P Omega| = 0.6``; the slab line through the origin along
    ``(-0.6, 0.8, 0)`` has ``rate = -0.6``. Law: ``image.shifted(s)`` reads at
    ``t + s`` what ``image`` reads at ``t`` (to ``8 eps (|t| + |s| + 4)``).
    Underflow row (qa F4): ``Omega = (1e-170, 0, 1)`` through ``(0.5, 0, 0)`` has
    ``|P Omega| = 1e-170`` exactly (never 0: not parallel) and ``b = 0`` (its
    in-plane line passes through the axis); squared components would underflow.
    First reds: ``|P Omega|`` from unscaled squares (the underflow row); ``shifted``
    moving the wrong way.
    """
    from orpheus.geometry.chart import AxialImage, RadialImage
    from orpheus.geometry.line import Line

    img = _CYL.image(Line.through(np.array([0.7, 0.4, 5.0]), np.array([0.0, 0.6, 0.8])))
    assert isinstance(img, RadialImage)
    assert img.impact_parameter == pytest.approx(0.7, rel=4 * _EPS) and img.speed == pytest.approx(0.6, rel=2 * _EPS)
    ax = _SLAB.image(Line.through(np.zeros(3), np.array([-0.6, 0.8, 0.0])))
    assert isinstance(ax, AxialImage) and ax.rate == -0.6 and ax.speed == 0.6 and not ax.parallel
    t = np.array([[-1.3, 0.2, 2.5]])
    for image in (img, ax):
        for shift in (0.75, -3.1):
            moved = image.shifted(np.asarray(shift))
            np.testing.assert_allclose(moved.orbit_coordinate_at(t + shift), image.orbit_coordinate_at(t),
                                       rtol=0, atol=8 * _EPS * (np.abs(t).max() + abs(shift) + 4.0))
    tiny = _CYL.image(Line.through(np.array([0.5, 0.0, 0.0]), np.array([1e-170, 0.0, 1.0])))
    assert isinstance(tiny, RadialImage)
    assert tiny.speed == 1e-170 and not tiny.parallel
    assert tiny.impact_parameter == 0.0
