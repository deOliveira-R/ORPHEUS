"""Foundation gates for :class:`orpheus.geometry.line.Line`, the oriented line in Plücker coordinates.

The geometric kernel seed (``.claude/plans/characteristic_reference_architecture.md``,
P0; verification spec ``scratch/characteristic_architecture/seed_verification_spec.md``
§4, rows V1-V6 re-posed onto the built API). These are software invariants of a
value type, so they carry ``foundation`` and no ``verifies``.

The laws, each written in the direction that is a floating-point theorem:

- a direction is a finite unit vector to 8 ulp, refused otherwise and never
  renormalised; the base point and the moment are finite;
- the line does not depend on the base point named (the moment is invariant
  under ``p -> p + s*Omega``), exactly when the cross products are exact;
- the parameter is measured from the foot: ``at(0)`` is the foot and
  ``parameter_of(at(t)) == t``;
- ``reversed`` keeps the foot and negates every parameter, bit for bit;
- a motion moves the direction linearly (a translation leaves it bit-identical)
  and the line through the image of every point.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from orpheus.geometry.line import Line
from orpheus.geometry.transformation import RigidMotion

pytestmark = pytest.mark.foundation

_EPS = np.finfo(float).eps
_F0 = "tests/gates/geometry/test_transformation.py::test_a3_inverse_on_both_sides_and_against_the_homogeneous_matrix"
_HERE = "tests/gates/geometry/test_line.py::"


def test_a_direction_within_8_ulp_of_unit_is_admitted_and_one_beyond_is_refused() -> None:
    """V2: the admission boundary is 8 ulp of unit length, both sides.

    Positive leg: ``(1 + 4 eps, 0, 0)`` departs by 4 eps and constructs, its
    components stored bit for bit (refused, never renormalised). Negative leg:
    ``(1 + 16 eps, 0, 0)`` and the classic non-unit ``(3, 4, 0)`` are refused
    with the type's own fragment.
    First red: a renormalising constructor (the positive leg's stored
    component changes) or a widened band (the negative leg constructs).
    """
    inside = np.array([1.0 + 4 * _EPS, 0.0, 0.0])
    line = Line.through(np.zeros(3), inside)
    np.testing.assert_array_equal(line.direction, inside)
    for outside in (np.array([1.0 + 16 * _EPS, 0.0, 0.0]), np.array([3.0, 4.0, 0.0])):
        with pytest.raises(ValueError, match="a line's direction is a finite unit vector"):
            Line.through(np.zeros(3), outside)


@pytest.mark.parametrize(
    "bad", [np.nan, np.inf, -np.inf], ids=["nan", "inf", "neg_inf"],
)
def test_a_non_finite_direction_is_refused(bad: float) -> None:
    """V2b: a non-finite direction component is refused.

    ``departure > tol`` is False for a NaN, so a refusal spelled "any beyond"
    admits it; the spelling "all within" refuses it. First red (measured):
    the pre-fix ``np.any(departure > tol)`` admitted ``(nan, 0, 0)``.

    Both doors, with warnings as errors: the constructor, and ``through``, whose
    refusal must come before any arithmetic on the direction (a cross product
    with ``inf`` emits numpy's invalid-value warning first).
    """
    bad_direction = np.array([bad, 0.0, 0.0])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(ValueError, match="a line's direction is a finite unit vector"):
            Line(direction=bad_direction, moment=np.zeros(3))
        with pytest.raises(ValueError, match="a line's direction is a finite unit vector"):
            Line.through(np.zeros(3), bad_direction)


def test_a_non_finite_base_point_is_refused_before_the_cross_product() -> None:
    """V2c: a non-finite base point is refused by ``through`` before the moment is formed.

    Warnings are errors inside the row: a refusal after the cross product would
    first emit numpy's invalid-value warning. The direct constructor refuses a
    non-finite moment with its own fragment.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(ValueError, match="a line's base point is a finite point"):
            Line.through(np.array([np.inf, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))
    with pytest.raises(ValueError, match="a line's moment is finite"):
        Line(direction=np.array([0.0, 1.0, 0.0]), moment=np.array([np.nan, 0.0, 0.0]))


@pytest.mark.rests_on(_HERE + "test_a_direction_within_8_ulp_of_unit_is_admitted_and_one_beyond_is_refused")
def test_the_line_does_not_depend_on_the_base_point_named() -> None:
    """V6: ``p`` and ``p + s*Omega`` name one line.

    Exact leg: axis-aligned data, every product exact, so the moments and the
    feet are ``array_equal`` (a far base point, 1e6 along the line, included:
    the spec's regime 6, where ``b^2 = |p|^2 - (p.Omega)^2`` loses 7e-6).
    Draw leg: 500 seeded lines, the feet agree to the cross products' rounding,
    ``8 eps (|p| + |s|)``.
    First red: a moment taken as ``p`` itself, or a foot from ``|p|^2 - (p.Omega)^2``.
    """
    omega = np.array([0.0, 1.0, 0.0])
    near = Line.through(np.array([0.7, 0.4, 0.0]), omega)
    for s in (-3.25, 1e6, -1e6):
        far = Line.through(np.array([0.7, 0.4 + s, 0.0]), omega)
        np.testing.assert_array_equal(far.moment, near.moment)
        np.testing.assert_array_equal(far.foot, near.foot)
    np.testing.assert_array_equal(near.foot, [0.7, 0.0, 0.0])

    rng = np.random.default_rng(20261005)
    om = rng.normal(size=(500, 3))
    om /= np.linalg.norm(om, axis=1, keepdims=True)
    p = rng.uniform(-3, 3, (500, 3))
    s = rng.uniform(-10, 10, (500, 1))
    a, b = Line.through(p, om), Line.through(p + s * om, om)
    scale = 8 * _EPS * (np.linalg.norm(p, axis=1) + np.abs(s[:, 0]) + 1.0)
    assert np.all(np.linalg.norm(a.foot - b.foot, axis=1) <= scale)


def test_the_parameter_is_measured_from_the_foot() -> None:
    """The foot is the point closest to the origin and the origin of ``t``.

    ``at(0) == foot``; ``foot . Omega`` is 0 to rounding (the foot is the
    closest point); ``parameter_of(at(t))`` returns ``t`` to ``4 eps (|t| + |foot|)``.
    """
    rng = np.random.default_rng(20261005)
    om = rng.normal(size=(200, 3))
    om /= np.linalg.norm(om, axis=1, keepdims=True)
    line = Line.through(rng.uniform(-3, 3, (200, 3)), om)
    np.testing.assert_array_equal(line.at(np.zeros(200)), line.foot)
    feet = np.linalg.norm(line.foot, axis=1)
    assert np.all(np.abs(np.sum(line.foot * om, axis=1)) <= 8 * _EPS * (feet + 1.0))
    t = rng.uniform(-5, 5, 200)
    assert np.all(np.abs(line.parameter_of(line.at(t)) - t) <= 4 * _EPS * (np.abs(t) + feet + 1.0))


def test_reversal_keeps_the_foot_and_negates_every_parameter() -> None:
    """``reversed`` is the same set with the opposite orientation, bit for bit.

    The foot ``Omega x m`` is unchanged when both factors change sign, and
    ``parameter_of`` negates exactly. First red: ``reversed`` negating only the
    direction (the foot moves to its mirror image).
    """
    rng = np.random.default_rng(20261005)
    om = rng.normal(size=(100, 3))
    om /= np.linalg.norm(om, axis=1, keepdims=True)
    line = Line.through(rng.uniform(-3, 3, (100, 3)), om)
    back = line.reversed()
    np.testing.assert_array_equal(back.foot, line.foot)
    x = rng.uniform(-3, 3, (100, 3))
    np.testing.assert_array_equal(back.parameter_of(x), -line.parameter_of(x))


@pytest.mark.rests_on(_F0)
def test_a_motion_moves_the_direction_linearly_and_the_points_affinely() -> None:
    """V4: under ``x -> Qx + t`` the direction becomes ``Q Omega`` and every point of the line maps onto the image line.

    A pure translation leaves the direction bit-identical (its linear part is the
    identity). For a general motion, the images of three points of the line lie
    on the image line to ``16 eps (|x| + |t|)``.
    First red: a direction moved with the translation (``on_points`` in place of
    ``on_directions``); the translated direction is not unit and is refused, or
    differs from ``Omega``.
    """
    omega = np.array([0.0, 0.6, 0.8])
    line = Line.through(np.array([0.7, 0.4, 5.0]), omega)
    shifted = line.moved_by(RigidMotion.translation_by([1.3, -0.4, 2.0]))
    np.testing.assert_array_equal(shifted.direction, line.direction)

    motion = RigidMotion(
        RigidMotion.rotation_about_axis(axis=(1.0, 2.0, -0.5), angle=np.sqrt(2.0)).linear,
        np.array([0.3, -1.1, 2.5]),
    )
    image = line.moved_by(motion)
    np.testing.assert_allclose(image.direction, motion.on_directions(omega), rtol=0, atol=4 * _EPS)
    for t in (-2.0, 0.0, 3.5):
        x = motion.on_points(line.at(np.array(t)))
        off_line = x - image.at(image.parameter_of(x))
        assert np.linalg.norm(off_line) <= 16 * _EPS * (np.linalg.norm(x) + 4.0)
