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

from dataclasses import dataclass
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

__all__ = ["AxialImage", "Chart", "DirectionDomain", "DirectionShape", "RadialImage", "SingularStratum"]

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
        return RadialImage(
            impact_parameter=_norm(closest),
            origin_position=along,
            speed=speed,
            parameter_origin=np.zeros_like(speed),
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


_SHAPE_BOUNDS: dict[str, tuple[float, float]] = {
    "cosine": (-1.0, 1.0),
    "axial_cosine": (0.0, 1.0),
    "angle": (0.0, pi),
}


@dataclass(frozen=True, eq=False)
class DirectionDomain:
    r"""The directions at the point :math:`c\,\hat e_x` of a chart modulo its stabiliser (:meth:`Chart.directions_at`).

    A box in the coordinates of its :attr:`shape` (:class:`DirectionShape`),
    on which :math:`\mathrm{d}\Omega` is :attr:`density` times Lebesgue
    measure and every orbit of the stabiliser meets the box once (its
    boundary aside).

    SCOPE-BOUNDARY[guard] machinery: the point-isotropy computation
    :math:`L \cap \mathrm{Stab}(x)` and the orbit-space catalogue entry
    :math:`S^2/D_{1h}`, which would derive the shape table below from the
    group.
    ruling: the user, 2026-10-06 (#581: the catalogue's barycentre lift is
    not a right inverse on :math:`D_{1h}`'s non-linear chart).
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

    @cached_property
    def shape(self) -> DirectionShape:
        """The domain's shape, from the chart's pair and whether the point is on the singular stratum."""
        chart = self.chart
        on_stratum = any(s.orbit_value == self.orbit_coordinate for s in chart.singular_strata)
        if not chart.acts_on_kept_space or chart.kept_columns == 3:
            return DirectionShape.WHOLE if on_stratum else DirectionShape.COSINE
        return DirectionShape.AXIAL_COSINE if on_stratum else DirectionShape.ANGLE_AXIAL

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
        return tuple(_SHAPE_BOUNDS[axis] for axis in self.axes)

    @property
    def density(self) -> float:
        r"""The constant density of :math:`\mathrm{d}\Omega`: :math:`4\pi` over the box's coordinate measure."""
        return 4.0 * pi / float(np.prod([hi - lo for lo, hi in self.bounds]))

    def direction(self, coordinates: np.ndarray) -> np.ndarray:
        r"""Representative unit directions of the canonical frame at coordinates ``(..., len(axes))``, ``(..., 3)``.

        Refused: coordinates of the wrong width, non-finite, or outside the box.
        """
        q = np.asarray(coordinates, dtype=float)
        if q.shape[-1:] != (len(self.axes),):
            raise ValueError(f"coordinates on the axes {self.axes} have shape (..., {len(self.axes)}); got {q.shape}")
        for j, (lo, hi) in enumerate(self.bounds):
            column = q[..., j]
            if not np.all((column >= lo) & (column <= hi)):              # NaN fails both comparisons
                raise ValueError(f"the {self.axes[j]} coordinate lies in [{lo}, {hi}]; got values outside it or not finite")
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

        The kernel's own :meth:`Chart.image` of that line, so the reading and
        the chord agree bit for bit. A slab line has no impact parameter.
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


def _sine(cosine: np.ndarray) -> np.ndarray:
    r""":math:`\sqrt{(1 - c)(1 + c)}`, without the cancellation of :math:`1 - c^2` near :math:`|c| = 1`."""
    return np.sqrt(np.clip((1.0 - cosine) * (1.0 + cosine), 0.0, None))


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
    """

    impact_parameter: np.ndarray
    origin_position: np.ndarray
    speed: np.ndarray
    parameter_origin: np.ndarray

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
        return RadialImage(self.impact_parameter, self.origin_position, self.speed, self.parameter_origin + shift)


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

