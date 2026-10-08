r"""The chart of a coordinate system: the orbit space of its symmetry group.

A one-dimensional coordinate system names a function :math:`c` on space,
the position a 1-D problem keeps: :math:`c(x) = x_0` on a slab, the
distance :math:`\sqrt{x_0^2 + x_1^2}` to the axis on a cylinder, the
distance :math:`|x|` to the centre on a sphere (each in its canonical
frame: the slab's normal is :math:`\hat e_x`, the cylinder's axis
:math:`\hat e_z`, the sphere's centre the origin; the columns
:class:`~orpheus.geometry.coord.AngularChart` declares). Its symmetry group

.. math::

    G_c = \{\, g \in E(3) : c \circ g = c \,\}

is the group of rigid motions the problem cannot tell apart, and :math:`c`
is the quotient map onto the orbit space :math:`\mathbb{R}^3 / G_c`.

**Two data determine the chart.** The kept columns, the first :math:`d` of
the canonical frame (:math:`d = 1, 2, 3`: the exponent of the coordinate
system's measure coordinate :math:`T = r^d`), and the linear group
:math:`L \subset O(3)` of :math:`G_c` (the slab's :math:`O(2)_x`, the
rotations about :math:`\hat e_x` and the mirrors through it; the
cylinder's :math:`D_{\infty h}`, the stabiliser of the axis as a line; the
sphere's :math:`O(3)`). :math:`G_c` is :math:`L` together with the
translations along the discarded columns. Every verb of the chart derives
from the pair:

- :math:`L` acts on the kept space either trivially (the slab: it contains
  no reflection of :math:`\hat e_x`) or as :math:`O(d)` (the cylinder and
  the sphere). The orbit coordinate is the kept coordinate itself in the
  first case, signed, and its norm in the second;
- the singular stratum :math:`c = 0`, with isotropy :math:`L`, exists
  exactly when the action is non-trivial (the sphere's centre, the
  cylinder's axis);
- a motion :math:`(Q, t)` is in :math:`G_c` iff :math:`Q \in L` and
  :math:`t` has no component in the kept columns;
- a line's image in the orbit space moves at :math:`|P\Omega|`, the length
  of the direction's kept components.

This is the spatial twin of
:class:`~orpheus.numerics.symmetry.SubgroupOfO3` on directions, and the
1-D part of the ``Chart`` of the boundary-law ontology's architecture E
(#551, the pair (coordinate system, kept coordinates)), pulled forward for
the geometric kernel seed by the user's rulings of 2026-10-05
(``.claude/plans/characteristic_reference_architecture.md``): the chart
derived from (kept columns, group), not matched on the coordinate system.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum
from functools import cached_property
from math import gamma, pi
from typing import TYPE_CHECKING, assert_never

import numpy as np

from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.transformation import RigidMotion
from orpheus.numerics.symmetry import SubgroupOfO3

if TYPE_CHECKING:
    from orpheus.geometry.line import Line

__all__ = [
    "AxialImage",
    "Chart",
    "DirectionDomain",
    "DirectionShape",
    "LineDomain",
    "LineShape",
    "half_chord",
    "RadialImage",
    "SingularStratum",
]

#: Tolerance on a motion's translation when deciding membership in
#: :math:`G_c`; the orthogonal part is decided by the group's realization.
_MEMBERSHIP_ATOL = 1e-12


@dataclass(frozen=True)
class SingularStratum:
    r"""A stratum of the orbit space where the isotropy of :math:`G_c` jumps.

    The orbit coordinate is not a submersion there, so a level set through
    it is a point or a line, not a surface: the sphere's centre, the
    cylinder's axis. A chord never crosses a stratum; it passes through it.

    Attributes
    ----------
    orbit_value:
        The value of the orbit coordinate on the stratum (0.0).
    isotropy:
        The linear isotropy of a point of the stratum: the chart's whole
        linear group (:math:`O(3)` at the sphere's centre,
        :math:`D_{\infty h}` on the cylinder's axis).
    """

    orbit_value: float
    isotropy: SubgroupOfO3


@dataclass(frozen=True)
class Chart:
    r"""The orbit space of a 1-D coordinate system's symmetry group :math:`G_c`.

    Attributes
    ----------
    coord:
        The coordinate system; it determines the pair (kept columns,
        linear group) every verb derives from, and its
        :meth:`~CoordSystem.measure` is the chart's measure.
    """

    coord: CoordSystem

    @property
    def kept_columns(self) -> int:
        r"""The number :math:`d` of canonical-frame columns the coordinate keeps."""
        return self.coord.measure_coordinate.exponent

    @property
    def group(self) -> SubgroupOfO3:
        r"""The linear group :math:`L` of :math:`G_c`: what the coordinate system is."""
        match self.coord:
            case CoordSystem.CARTESIAN:
                return SubgroupOfO3.O2("x")
            case CoordSystem.CYLINDRICAL:
                return SubgroupOfO3.Dinfh
            case CoordSystem.SPHERICAL:
                return SubgroupOfO3.O3

    @cached_property
    def acts_on_kept_space(self) -> bool:
        r"""Whether :math:`L` moves the kept space (it then acts as :math:`O(d)` there).

        Decided by the group's realization: :math:`L` contains the
        reflection of the first kept column :math:`\hat e_x` or it fixes the
        kept space pointwise.
        """
        return self.group.realization.contains_element(RigidMotion.reflection(normal=[1.0, 0.0, 0.0]))

    def orbit_coordinate(self, points: np.ndarray) -> np.ndarray:
        r"""The coordinate map :math:`c` on points ``(..., 3)`` of the canonical frame, ``(...,)``.

        The kept coordinate, signed, where :math:`L` fixes the kept space;
        the norm of the kept components where it acts as :math:`O(d)`.
        """
        kept = _as_vectors(points)[..., : self.kept_columns]
        if not self.acts_on_kept_space:
            return kept[..., 0].copy()
        return _norm(kept)

    def contains(self, motion: RigidMotion) -> bool:
        r"""Whether ``motion`` :math:`(Q, t)` is in :math:`G_c`: :math:`Q \in L` and :math:`t` has no kept component."""
        if motion.dimension != 3:
            raise ValueError(f"a chart's group acts on R^3; the motion acts on R^{motion.dimension}")
        translation_kept = motion.translation[: self.kept_columns]
        return self.group.realization.contains_element(motion.linear_part) and bool(
            np.all(np.abs(translation_kept) <= _MEMBERSHIP_ATOL)
        )

    @property
    def singular_strata(self) -> tuple[SingularStratum, ...]:
        r"""The stratum :math:`c = 0` with isotropy :math:`L`, where :math:`L` acts on the kept space; none otherwise."""
        if not self.acts_on_kept_space:
            return ()
        return (SingularStratum(0.0, self.group),)

    def measure(self, edges: np.ndarray) -> np.ndarray:
        r"""The measures of the cells between ``edges``: :meth:`CoordSystem.measure`, the one definition."""
        return self.coord.measure(edges)

    def measure_density(self, orbit_coordinate: np.ndarray) -> np.ndarray:
        r"""The density of :meth:`measure` in the orbit coordinate, and the area of its level set: :meth:`CoordSystem.measure_density`."""
        return self.coord.measure_density(orbit_coordinate)

    def projected_speed(self, directions: np.ndarray) -> np.ndarray:
        r"""The length :math:`|P\Omega|` of the kept components of unit directions, ``(...,)``.

        The speed of a line's image in the kept space :math:`\mathbb{R}^d`,
        before the quotient by :math:`L`: an orbit-space length :math:`\ell`
        along a chord is the 3-D length :math:`\ell / |P\Omega|` (its
        reciprocal is the line's obliquity, infinite for a line parallel to
        the orbit space, which is why the speed is the stored quantity).
        The orbit coordinate itself changes at this rate only far from the
        closest approach (on the sphere :math:`\mathrm{d}|x|/\mathrm{d}t`
        runs over :math:`[-1, 1]`). Computed with scaling, never as
        :math:`1 - \Omega_z^2`, which loses every digit near the axis.
        """
        return _norm(_as_vectors(directions)[..., : self.kept_columns])

    def beam_density(self, impact_parameter: np.ndarray, directions: np.ndarray) -> np.ndarray:
        r"""The density over the impact parameter :math:`b \ge 0` of the parallel-beam measure on lines.

        The invariant measure on lines is :math:`\mathrm{d}A_\perp\,
        \mathrm{d}\Omega`; for a fixed direction, :math:`\mathrm{d}A_\perp` is
        the area element of the plane normal to :math:`\Omega`. Per unit
        measure of the discarded columns, the lines of one direction meet
        the kept space in a beam of width :math:`|P\Omega|` per unit of the
        kept space's cross-section normal to :math:`P\Omega`; where
        :math:`L` acts as :math:`O(d)`, that cross-section at impact
        parameter :math:`b` is a :math:`(d-2)`-sphere of radius :math:`b`,
        of area :math:`S_{d-2}\,b^{d-2}` (:math:`S_0 = 2`, the two lines at
        :math:`\pm b`; :math:`S_1 = 2\pi`). So the density is
        :math:`S_{d-2}\,b^{d-2}\,|P\Omega|`: :math:`2|P\Omega|` per unit height
        on the cylinder, :math:`2\pi b` on the sphere. Where :math:`L` fixes
        the kept space (the slab) a line of one direction has no impact
        parameter, and the density per unit area of a plane
        :math:`x = \text{const}` is :math:`|\Omega_x|`.

        The quadrature over :math:`b` is the consumer's: a chord length has
        a square-root endpoint at every radius, and a rule that absorbs it
        lives above this layer (``chord_quadrature`` in
        :mod:`orpheus.derivations.common`).
        """
        speed = self.projected_speed(directions)
        if not self.acts_on_kept_space:
            return speed
        k = self.kept_columns - 2
        sphere_area = 2.0 * pi ** ((k + 1) / 2.0) / gamma((k + 1) / 2.0)
        return sphere_area * np.asarray(impact_parameter, dtype=float) ** k * speed

    def directions_at(self, orbit_coordinate: float) -> "DirectionDomain":
        r"""The directions at a point modulo the point's stabiliser: :math:`S^2 / \mathrm{Stab}(x)` (:class:`DirectionDomain`).

        The representative point is :math:`x = c\,\hat e_x` in the canonical
        frame. Refused: a non-finite :math:`c`, and a negative one where the
        chart's group acts on the kept space.
        """
        r = float(orbit_coordinate)
        if not np.isfinite(r) or (self.acts_on_kept_space and r < 0.0):
            raise ValueError(
                f"a point's orbit coordinate is finite, and non-negative where the chart's group acts on the "
                f"kept space; got {orbit_coordinate!r}"
            )
        return DirectionDomain(self, r)

    def line_domain(self) -> "LineDomain":
        r"""The oriented lines of space modulo :math:`G_c`, with their invariant measure (:class:`LineDomain`).

        The lines' counterpart of :meth:`directions_at`.
        """
        return LineDomain(self)

    def image(self, line: "Line") -> "RadialImage | AxialImage":
        r"""The image of lines of the canonical frame in the orbit space.

        Where :math:`L` acts as :math:`O(d)` the image is a straight line
        in the kept space at impact parameter :math:`b`, closest at
        :math:`t^*` (:class:`RadialImage`); where it fixes the kept space,
        the orbit coordinate is affine along the line (:class:`AxialImage`).
        """
        kept_foot = line.foot[..., : self.kept_columns]
        kept_direction = line.direction[..., : self.kept_columns]
        if not self.acts_on_kept_space:
            return AxialImage(foot_coordinate=kept_foot[..., 0].copy(), rate=kept_direction[..., 0].copy())
        speed = self.projected_speed(line.direction)
        parallel = speed == 0.0
        # Through the unit kept direction u = P Omega/|P Omega|: the closest point is
        # foot - (foot . u) u, never an infinity times a zero; the foot sits at the
        # signed position s_0 = foot . u along the image (RadialImage).
        unit = kept_direction / np.where(parallel, 1.0, speed)[..., None]
        along = np.where(parallel, 0.0, np.sum(kept_foot * unit, axis=-1))
        closest = kept_foot - along[..., None] * unit
        b = _norm(closest)
        return RadialImage(
            impact_parameter=b,
            origin_position=along,
            speed=speed,
            parameter_origin=np.zeros_like(speed),
            level=b,
            level_half_chord=np.zeros_like(b),
        )


class DirectionShape(Enum):
    r"""The shapes of :math:`S^2 / \mathrm{Stab}(x)` the three charts produce, each a box of measure-uniform coordinates.

    ``WHOLE``: the stabiliser is :math:`O(3)` (the sphere's centre), one
    orbit, no coordinate. ``COSINE``: :math:`O(2)_x`, the cosine
    :math:`\Omega_x \in [-1, 1]` (Archimedes: uniform). ``AXIAL_COSINE``:
    :math:`D_{\infty h}` on the cylinder's axis, :math:`w = |\Omega_z| \in
    [0, 1]`. ``ANGLE_AXIAL``: :math:`D_{1h} = \{e, \sigma_y, \sigma_z,
    C_2(x)\}` off the cylinder's axis, the in-plane angle :math:`\alpha \in
    [0, \pi]` between :math:`P\Omega` and :math:`\hat x` times :math:`w`, with
    :math:`\mathrm{d}\Omega = \mathrm{d}w\,\mathrm{d}\alpha`.
    """

    WHOLE = ()
    COSINE = ("cosine",)
    AXIAL_COSINE = ("axial_cosine",)
    ANGLE_AXIAL = ("angle", "axial_cosine")


#: The interval of every coordinate a direction or a line domain names; one table for both.
_AXIS_BOUNDS: dict[str, tuple[float, float]] = {
    "cosine": (-1.0, 1.0),
    "axial_cosine": (0.0, 1.0),
    "angle": (0.0, pi),
    "impact": (0.0, np.inf),
    "polar_angle": (0.0, pi / 2.0),
}


def _in_box(coordinates: np.ndarray, axes: tuple[str, ...]) -> np.ndarray:
    """``coordinates`` ``(..., len(axes))`` as floats; refused if of the wrong width, non-finite, or outside the box."""
    q = np.asarray(coordinates, dtype=float)
    if q.shape[-1:] != (len(axes),):
        raise ValueError(f"coordinates on the axes {axes} have shape (..., {len(axes)}); got {q.shape}")
    for j, axis in enumerate(axes):
        lo, hi = _AXIS_BOUNDS[axis]
        column = q[..., j]
        if not np.all((column >= lo) & (column <= hi) & np.isfinite(column)):
            raise ValueError(f"the {axis} coordinate lies in [{lo}, {hi}] and is finite; got values outside it or not finite")
    return q


@dataclass(frozen=True, eq=False)
class DirectionDomain:
    r"""The directions at the point :math:`c\,\hat e_x` of a chart modulo its stabiliser (:meth:`Chart.directions_at`).

    A box in the coordinates of its :attr:`shape` (:class:`DirectionShape`),
    on which :math:`\mathrm{d}\Omega` is :attr:`density` times Lebesgue
    measure and every orbit of the stabiliser meets the box once (its
    boundary aside).

    SCOPE-BOUNDARY[guard] machinery: the point-isotropy computation and the S^2/D_1h orbit-catalogue entry, deriving the shape table.
    ruling: the user, 2026-10-06 (#581: the catalogue's barycentre lift is not a right inverse on D_1h's chart).
    revisit: when #581 lands, this table retires onto the catalogue.

    Attributes
    ----------
    chart:
        The chart the point is on.
    orbit_coordinate:
        The point's orbit coordinate :math:`c`.
    """

    chart: "Chart"
    orbit_coordinate: float

    @property
    def point(self) -> np.ndarray:
        r"""The representative point :math:`c\,\hat e_x`, ``(3,)``."""
        return np.array([self.orbit_coordinate, 0.0, 0.0])

    @property
    def on_stratum(self) -> bool:
        """Whether the point is on a singular stratum of the chart (the sphere's centre, the cylinder's axis)."""
        return any(s.orbit_value == self.orbit_coordinate for s in self.chart.singular_strata)

    @cached_property
    def shape(self) -> DirectionShape:
        """The domain's shape, from the chart's pair and whether the point is on the singular stratum."""
        chart = self.chart
        if not chart.acts_on_kept_space or chart.kept_columns == 3:
            return DirectionShape.WHOLE if self.on_stratum else DirectionShape.COSINE
        return DirectionShape.AXIAL_COSINE if self.on_stratum else DirectionShape.ANGLE_AXIAL

    @property
    def stabiliser(self) -> SubgroupOfO3:
        """The point's linear isotropy in the chart's group: the whole group on the stratum."""
        match self.shape:
            case DirectionShape.WHOLE | DirectionShape.AXIAL_COSINE:
                return self.chart.group
            case DirectionShape.COSINE:
                return SubgroupOfO3.O2("x")
            case DirectionShape.ANGLE_AXIAL:
                return SubgroupOfO3.Dnh(1)
            case unreachable:
                assert_never(unreachable)

    @property
    def axes(self) -> tuple[str, ...]:
        """The coordinate names, in order (none at the sphere's centre)."""
        return self.shape.value

    @property
    def bounds(self) -> tuple[tuple[float, float], ...]:
        """One closed interval per axis."""
        return tuple(_AXIS_BOUNDS[axis] for axis in self.axes)

    @property
    def density(self) -> float:
        r"""The constant density of :math:`\mathrm{d}\Omega`: :math:`4\pi` over the box's coordinate measure."""
        return 4.0 * pi / float(np.prod([hi - lo for lo, hi in self.bounds]))

    def direction(self, coordinates: np.ndarray) -> np.ndarray:
        r"""Representative unit directions of the canonical frame at coordinates ``(..., len(axes))``, ``(..., 3)``.

        Refused: coordinates of the wrong width, non-finite, or outside the box.
        """
        q = _in_box(coordinates, self.axes)
        zero = np.zeros(q.shape[:-1])
        match self.shape:
            case DirectionShape.WHOLE:
                return np.stack([zero + 1.0, zero, zero], axis=-1)
            case DirectionShape.COSINE:
                cosine = q[..., 0]
                return np.stack([cosine, _sine(cosine), zero], axis=-1)
            case DirectionShape.AXIAL_COSINE:
                w = q[..., 0]
                return np.stack([_sine(w), zero, w], axis=-1)
            case DirectionShape.ANGLE_AXIAL:
                alpha, w = q[..., 0], q[..., 1]
                return np.stack([_sine(w) * np.cos(alpha), _sine(w) * np.sin(alpha), w], axis=-1)
            case unreachable:
                assert_never(unreachable)

    def impact_parameter(self, coordinates: np.ndarray) -> np.ndarray:
        r"""The impact parameter :math:`b` of the line through :attr:`point` in each direction, ``(...,)``.

        The kernel's own :meth:`Chart.image` of that line, one definition of
        :math:`b`: bit for bit equal to ``Chart.image`` of the same lines,
        and within a few ulp of a partition's chord, which first moves the
        line by its pose's inverse. A slab line has no impact parameter.
        """
        from orpheus.geometry.line import Line

        directions = self.direction(coordinates)
        image = self.chart.image(Line.through(np.broadcast_to(self.point, directions.shape), directions))
        match image:
            case RadialImage():
                return image.impact_parameter
            case AxialImage():
                raise ValueError("a slab line has no impact parameter: the slab's group fixes the kept space")

    def tangencies(self, level: float) -> np.ndarray:
        r"""The values of the first axis where the line through :attr:`point` is tangent to :math:`c` = ``level``.

        Defined where the chart's group acts on the kept space, for
        :math:`0 < \ell \le c`, where :math:`b = \ell` has solutions:
        :math:`\Omega_x = \pm\sqrt{(1 - \ell/c)(1 + \ell/c)}` on the sphere,
        :math:`\alpha = \arcsin(\ell/c)` and :math:`\pi - \arcsin(\ell/c)` on
        the cylinder, one value (the grazing direction) at :math:`\ell = c`.
        Sorted; empty elsewhere, on the slab and on a stratum.
        """
        r, l = self.orbit_coordinate, float(level)
        if not self.chart.acts_on_kept_space or not 0.0 < l <= r:
            return np.empty(0)
        ratio = l / r
        match self.shape:
            case DirectionShape.COSINE:
                m = float(_sine(np.asarray(ratio)))
                return np.unique([-m, m])
            case DirectionShape.ANGLE_AXIAL:
                a = float(np.arcsin(ratio))
                return np.unique([a, pi - a])
            case DirectionShape.WHOLE | DirectionShape.AXIAL_COSINE:
                return np.empty(0)
            case unreachable:
                assert_never(unreachable)


class LineShape(Enum):
    r"""The shapes of the orbit space of oriented lines under the three charts' groups, each a box.

    ``IMPACT`` (the sphere, :math:`O(3)`): the impact parameter
    :math:`b \in [0, \infty)`; a line and its reverse are one orbit.
    ``IMPACT_POLAR`` (the cylinder, :math:`D_{\infty h}` with the axial
    translations): :math:`b` and the polar angle :math:`\theta \in [0, \pi/2]`
    between the line and the axis; the mirror normal to the axis folds
    :math:`\theta` and :math:`\pi - \theta`. ``COSINE`` (the slab,
    :math:`O(2)_x` with the transverse translations): the cosine
    :math:`\mu = \Omega_x \in [-1, 1]`; the group fixes the kept space, so
    :math:`\mu` and :math:`-\mu` are two orbits.
    """

    IMPACT = ("impact",)
    IMPACT_POLAR = ("impact", "polar_angle")
    COSINE = ("cosine",)


@dataclass(frozen=True, eq=False)
class LineDomain:
    r"""The oriented lines of space modulo a chart's group, with the invariant measure (:meth:`Chart.line_domain`).

    A box in the coordinates of its :attr:`shape` (:class:`LineShape`); every
    orbit meets it once (its boundary aside). The invariant measure on
    oriented lines is :math:`\mathrm{d}A_\perp\,\mathrm{d}\Omega`, counted
    per unit measure of the discarded columns (per unit height on the
    cylinder, per unit transverse area on the slab); :meth:`density` is its
    density in the box's coordinates. Its integral of a line's chord length
    through a body is :math:`4\pi` times the body's measure (Cauchy).

    **Why the cylinder's angle is polar, not its cosine.** The beam's speed
    :math:`|P\Omega| = \sin\theta = \sqrt{1 - \mu_z^2}` puts a square-root
    end at :math:`\mu_z = 1` into every integrand over :math:`\mu_z`; in
    :math:`\theta` the density :math:`8\pi\sin^2\theta` is analytic (the
    user's ruling of 2026-10-06, after the escape probability missed
    Bickley's closed form by 5.4e-4 at 8 points in :math:`\mu_z` and
    reached 1.7e-12 at 32 in :math:`\theta`). :meth:`Chart.directions_at`
    keeps :math:`\mu_z`: a point's direction measure is uniform in it.

    The box is unbounded in :math:`b`: which lines meet a body is the
    body's question, so a consumer truncates at its outer radius.

    SCOPE-BOUNDARY[guard] machinery: the orbit space of oriented lines under a subgroup of E(3), deriving the shape table and its fold factors.
    ruling: the user, 2026-10-06, P1 step (b) third rung (the line domain a kernel verb on `Chart`; the cylinder's polar angle).
    revisit: when the line-orbit computation exists, the table and `density`'s folds retire onto it.

    Attributes
    ----------
    chart:
        The chart whose group the lines are taken modulo.
    """

    chart: "Chart"

    @property
    def shape(self) -> LineShape:
        """The domain's shape, from the chart's pair."""
        chart = self.chart
        if not chart.acts_on_kept_space:
            return LineShape.COSINE
        return LineShape.IMPACT if chart.kept_columns == 3 else LineShape.IMPACT_POLAR

    @property
    def axes(self) -> tuple[str, ...]:
        """The coordinate names, in order."""
        return self.shape.value

    @property
    def bounds(self) -> tuple[tuple[float, float], ...]:
        """One interval per axis; the impact parameter's is unbounded above."""
        return tuple(_AXIS_BOUNDS[axis] for axis in self.axes)

    def lines(self, coordinates: np.ndarray) -> "Line":
        r"""Representative oriented lines at coordinates ``(..., len(axes))``.

        Through :math:`b\,\hat e_y` along :math:`\hat e_x` on the sphere; through
        :math:`b\,\hat e_y` along :math:`(\sin\theta, 0, \cos\theta)` on the
        cylinder; through the origin along :math:`(\mu, \sqrt{1 - \mu^2}, 0)` on
        the slab. Their impact parameter under :meth:`Chart.image` is
        :math:`b` to an ulp: a :class:`~orpheus.geometry.line.Line` stores
        its moment and returns its foot as :math:`\Omega \times m`, which
        rounds (measured 2026-10-06: 1 ulp at :math:`b = 0.3`,
        :math:`\theta = 0.2`; exact at :math:`\theta \in \{0, \pi/2\}`).
        Refused: coordinates of the wrong width, non-finite, or outside the
        box.
        """
        from orpheus.geometry.line import Line

        q = _in_box(coordinates, self.axes)
        zero = np.zeros(q.shape[:-1])
        match self.shape:
            case LineShape.IMPACT:
                foot = np.stack([zero, q[..., 0], zero], axis=-1)
                direction = np.stack([zero + 1.0, zero, zero], axis=-1)
            case LineShape.IMPACT_POLAR:
                theta = q[..., 1]
                foot = np.stack([zero, q[..., 0], zero], axis=-1)
                direction = np.stack([np.sin(theta), zero, np.cos(theta)], axis=-1)
            case LineShape.COSINE:
                cosine = q[..., 0]
                foot = np.stack([zero, zero, zero], axis=-1)
                direction = np.stack([cosine, _sine(cosine), zero], axis=-1)
            case unreachable:
                assert_never(unreachable)
        return Line.through(foot, direction)

    def density(self, coordinates: np.ndarray) -> np.ndarray:
        r"""The density of :math:`\mathrm{d}A_\perp\,\mathrm{d}\Omega` in the box's coordinates, ``(...,)``.

        :meth:`Chart.beam_density` of the representative lines times the
        direction measure the quotient folds into one point of the box:
        :math:`4\pi` on the sphere (every direction); :math:`2\pi \cdot 2
        \sin\theta` on the cylinder (the azimuth, the two polar angles
        :math:`\theta` and :math:`\pi - \theta`, and
        :math:`\mathrm{d}\Omega = \sin\theta\,\mathrm{d}\theta\,\mathrm{d}\varphi`);
        :math:`2\pi` on the slab (the azimuth about :math:`\hat e_x`). So
        :math:`8\pi^2 b`, :math:`8\pi\sin^2\theta` and :math:`2\pi|\mu|`.
        """
        q = _in_box(coordinates, self.axes)
        lines = self.lines(q)
        impact = q[..., self.axes.index("impact")] if "impact" in self.axes else np.zeros(q.shape[:-1])
        beam = self.chart.beam_density(impact, lines.direction)
        match self.shape:
            case LineShape.IMPACT:
                folded = np.full(q.shape[:-1], 4.0 * pi)
            case LineShape.IMPACT_POLAR:
                folded = 4.0 * pi * np.sin(q[..., 1])
            case LineShape.COSINE:
                folded = np.full(q.shape[:-1], 2.0 * pi)
            case unreachable:
                assert_never(unreachable)
        return beam * folded


def _sine(cosine: np.ndarray) -> np.ndarray:
    r""":math:`\sqrt{(1 - c)(1 + c)}`, without the cancellation of :math:`1 - c^2` near :math:`|c| = 1`."""
    return np.sqrt(np.clip((1.0 - cosine) * (1.0 + cosine), 0.0, None))


#: How far, in ulp of :math:`r_*^2`, a line's level may sit from its impact parameter (:meth:`RadialImage.at_level`).
_LEVEL_AGREEMENT = 64.0


def half_chord(radius: np.ndarray, level: np.ndarray, level_half_chord: np.ndarray) -> np.ndarray:
    r"""The half-chord :math:`\sqrt{r^2 - b^2}` of a line at radius ``radius``, from its level: :math:`\sqrt{(r - r_*)(r + r_*) + y_*^2}`.

    Exact at :math:`r = r_*`, and free of the cancellation :math:`r - b` at
    every radius at or above it; 0 where the square is not positive (the
    line does not cross). Broadcast over the three arguments. The one
    spelling of a half-chord: the chord's crossings
    (:meth:`RadialImage.half_chord_at`) and the characteristic reading's
    point weights read it.
    """
    square = (radius - level) * (radius + level) + level_half_chord * level_half_chord
    return np.sqrt(np.where(square > 0.0, square, 0.0))


@dataclass(frozen=True, eq=False)
class RadialImage:
    r"""A line's image where the chart's group acts as :math:`O(d)` on the kept space.

    The image is a straight line in the kept space at impact parameter
    :math:`b`; with :math:`s` the signed position along it measured from the
    closest point,

    .. math::

        c(t)^2 = b^2 + s(t)^2, \qquad s(t) = s_0 + |P\Omega|\,(t - t_0),

    where :math:`t_0` is the parameter origin and :math:`s_0` the position
    there. The closest approach is :math:`t^* = t_0 - s_0/|P\Omega|`. Held in
    this form, every parameter the chord needs is
    :math:`t_0 + (\pm h - s_0)/|P\Omega|`: a finite numerator over one
    division, which overflows to a correctly signed infinity at a subnormal
    :math:`|P\Omega|` and never forms :math:`\infty - \infty`.

    Attributes
    ----------
    impact_parameter:
        :math:`b \ge 0`, ``(...,)``.
    origin_position:
        :math:`s_0`, ``(...,)``.
    speed:
        :math:`|P\Omega|`, ``(...,)``.
    parameter_origin:
        :math:`t_0`, ``(...,)``.
    level, level_half_chord:
        A radius :math:`r_*` and the line's exact half-chord there,
        :math:`y_* = \sqrt{r_*^2 - b^2}`, ``(...,)``: every half-chord is
        formed from them (:func:`half_chord`). :meth:`Chart.image` sets
        :math:`r_* = b` and :math:`y_* = 0`, the half-chord formed from
        :math:`b` itself. A line's :math:`b` is known to an ulp (the line
        stores its moment), so near a tangency, where :math:`y` is below
        :math:`\sqrt{2r\,\epsilon(r)}`, it rounds the half-chord away; a
        caller that knows the exact pair passes it (:meth:`at_level`, through
        :meth:`~orpheus.geometry.chord.ConcentricPartition.chord`), and the
        chord is formed without the cancellation :math:`r - b` (#590).
    """

    impact_parameter: np.ndarray
    origin_position: np.ndarray
    speed: np.ndarray
    parameter_origin: np.ndarray
    level: np.ndarray
    level_half_chord: np.ndarray

    def at_level(self, radius: np.ndarray, half_chord: np.ndarray) -> "RadialImage":
        r"""The same image with its half-chords formed from the exact pair :math:`(r_*, y_*)`, each ``(...,)``.

        The pair and the line's own :math:`b` describe one line, so they are
        held to agree: :math:`|r_*^2 - y_*^2 - b^2|` at most
        :data:`_LEVEL_AGREEMENT` ulp of :math:`r_*^2`, the squares compared
        because recomputing :math:`b` from the pair amplifies its rounding by
        :math:`r_*/b`. `[M]` 2026-10-08: at most 8 ulp over about 400 000
        lines of the characteristic reference's line and point rules (spheres
        and cylinders, solid and hollow; ``scratch/characteristic_architecture/p1_step_b5b/main/level_tolerance.py``).
        """
        radius, half_chord = np.asarray(radius, dtype=float), np.asarray(half_chord, dtype=float)
        if not (np.all(np.isfinite(half_chord)) and np.all(half_chord >= 0.0) and np.all(radius >= 0.0)):
            raise ValueError("a level is a non-negative radius with a non-negative, finite half-chord there")
        radius = np.broadcast_to(radius, self.impact_parameter.shape).copy()
        half_chord = np.broadcast_to(half_chord, self.impact_parameter.shape).copy()
        b = self.impact_parameter
        if np.any(np.abs((radius - half_chord) * (radius + half_chord) - b * b) > _LEVEL_AGREEMENT * np.spacing(radius * radius)):
            raise ValueError("a line's level disagrees with its impact parameter: r*^2 - y*^2 is not b^2")
        return replace(self, level=radius, level_half_chord=half_chord)

    def half_chord_at(self, radius: np.ndarray) -> np.ndarray:
        r"""The half-chord :math:`h = \sqrt{r^2 - b^2}` at radii ``radius`` ``(m,)`` or ``(..., m)``, ``(..., m)``; 0 where uncrossed.

        Formed as :math:`\sqrt{(r - r_*)(r + r_*) + y_*^2}` from the
        :attr:`level`: exact at :math:`r_*`, and free of the cancellation
        :math:`r - b` at every radius the line crosses at or above it. A
        radius is crossed where the square is positive and the line is not
        parallel.
        """
        h = half_chord(np.asarray(radius, dtype=float), self.level[..., None], self.level_half_chord[..., None])
        return np.where(~self.parallel[..., None], h, 0.0)

    def parameters_at(self, radius: np.ndarray, side: np.ndarray) -> np.ndarray:
        r"""The parameters at which the lines reach the radii ``radius``, on the ``side`` (:math:`\pm 1`, 0 at the closest approach) of the closest approach, ``(..., m)``."""
        return self.parameter_at(np.asarray(side, dtype=float) * self.half_chord_at(radius))

    @property
    def parallel(self) -> np.ndarray:
        """Where the line keeps one orbit coordinate (:math:`|P\\Omega| = 0`)."""
        return self.speed == 0.0

    @property
    def coordinate_of_parallel(self) -> np.ndarray:
        r"""The orbit coordinate a parallel line keeps, :math:`b`, ``(...,)``."""
        return self.impact_parameter

    @property
    def closest_approach(self) -> np.ndarray:
        r""":math:`t^* = t_0 - s_0/|P\Omega|`, ``(...,)``; :math:`t_0` on a parallel line."""
        return self.parameter_at(np.zeros_like(self.impact_parameter)[..., None])[..., 0]

    def parameter_at(self, position: np.ndarray) -> np.ndarray:
        r"""The parameters where the image is at signed positions ``position`` ``(..., m)``: :math:`t_0 + (s - s_0)/|P\Omega|`."""
        safe = np.where(self.parallel, 1.0, self.speed)[..., None]
        with np.errstate(over="ignore"):                   # a subnormal |P Omega|: a truly far parameter
            offset = (np.asarray(position, dtype=float) - self.origin_position[..., None]) / safe
        return self.parameter_origin[..., None] + np.where(self.parallel[..., None], 0.0, offset)

    def orbit_coordinate_at(self, t: np.ndarray) -> np.ndarray:
        r"""The orbit coordinate at parameters ``t`` ``(..., q)``: :math:`\sqrt{b^2 + s(t)^2}`, ``(..., q)``."""
        position = self.origin_position[..., None] + self.speed[..., None] * (
            np.asarray(t, dtype=float) - self.parameter_origin[..., None]
        )
        return np.hypot(self.impact_parameter[..., None], position)

    def shifted(self, shift: np.ndarray) -> "RadialImage":
        """The same image with every parameter increased by ``shift`` ``(...,)``."""
        return replace(self, parameter_origin=self.parameter_origin + shift)


@dataclass(frozen=True, eq=False)
class AxialImage:
    r"""A line's image where the chart's group fixes the kept space: :math:`c(t) = c_{\rm foot} + \dot c\,t`.

    Attributes
    ----------
    foot_coordinate:
        :math:`c` at the parameter origin, ``(...,)``.
    rate:
        :math:`\dot c = \Omega_x`, signed, ``(...,)``.
    """

    foot_coordinate: np.ndarray
    rate: np.ndarray

    @property
    def speed(self) -> np.ndarray:
        r""":math:`|P\Omega| = |\dot c|`."""
        return np.abs(self.rate)

    @property
    def parallel(self) -> np.ndarray:
        """Where the line keeps one orbit coordinate (:math:`\\dot c = 0`)."""
        return self.rate == 0.0

    @property
    def coordinate_of_parallel(self) -> np.ndarray:
        r"""The orbit coordinate a parallel line keeps, :math:`c_{\rm foot}`, ``(...,)``."""
        return self.foot_coordinate

    def orbit_coordinate_at(self, t: np.ndarray) -> np.ndarray:
        r"""The orbit coordinate at parameters ``t`` ``(..., q)``, ``(..., q)``."""
        return self.foot_coordinate[..., None] + self.rate[..., None] * np.asarray(t, dtype=float)

    def parameters_at(self, coordinate: np.ndarray, side: np.ndarray | None = None) -> np.ndarray:
        r"""The parameters at which the lines reach the orbit coordinates ``coordinate`` ``(m,)`` or ``(..., m)``: :math:`(c - c_{\rm foot})/\dot c`.

        A line meets each level once, so ``side`` is not read; a parallel
        line is given the rate 1 (it meets no level, and its crossings are
        absent).
        """
        rate = np.where(self.parallel, 1.0, self.rate)
        return (np.asarray(coordinate, dtype=float) - self.foot_coordinate[..., None]) / rate[..., None]

    def shifted(self, shift: np.ndarray) -> "AxialImage":
        """The same image with every parameter increased by ``shift`` ``(...,)``."""
        return AxialImage(self.foot_coordinate - self.rate * shift, self.rate)


def _norm(vectors: np.ndarray) -> np.ndarray:
    r"""The Euclidean norm over the last axis, scaled so no square underflows or overflows."""
    scale = np.max(np.abs(vectors), axis=-1)
    safe = np.where(scale > 0.0, scale, 1.0)
    return scale * np.sqrt(np.sum((vectors / safe[..., None]) ** 2, axis=-1))


def _as_vectors(points: np.ndarray) -> np.ndarray:
    x = np.asarray(points, dtype=float)
    if x.shape[-1:] != (3,):
        raise ValueError(f"points and directions live in R^3, shape (..., 3); got shape {x.shape}")
    return x

