r"""Lines through a concentric partition: chords, crossings and point location.

A 1-D geometry partitions space by the level sets :math:`c = r_k` of its
chart's orbit coordinate (:class:`~orpheus.geometry.chart.Chart`): planes
:math:`x = r_k` on a slab, cylinders and spheres about the axis or centre.
A line meets that partition along its **chord**.

**The chord is solved once, in the orbit space** (:meth:`Chart.image`).
Where the chart's group acts as :math:`O(d)` on the kept space (the
cylinder, the sphere), the line's image is a straight line at impact
parameter :math:`b`,

.. math::

    c(t)^2 = b^2 + \bigl(|P\Omega|\,(t - t^*)\bigr)^2,

so the line crosses :math:`c = r_k` iff :math:`b < r_k` (a tangency,
:math:`b = r_k`, is not a crossing; a solid body's :math:`r_0 = 0` is the
singular stratum, never crossed), at :math:`t = t^* \pm h_k / |P\Omega|`
with the half-chord :math:`h_k = \sqrt{(r_k - b)(r_k + b)}`. Every 3-D
length is the orbit-space length over :math:`|P\Omega|`: one factor for
the charts (1 on the sphere, :math:`\sin\theta` on the cylinder). Where the
group fixes the kept space (the slab) the image is affine,
:math:`c(t) = c_{\rm foot} + \Omega_x t`, crossing :math:`x = r_k` at
:math:`t_k = (r_k - c_{\rm foot})/\Omega_x`.

**A segment's region comes from the crossing order** (the user's ruling,
2026-10-05; :eq:`geometry-crossing-order`): crossing :math:`c = r_k` with
the coordinate decreasing enters region :math:`k - 1`, increasing enters
region :math:`k`. :attr:`Crossings.region_entered` is the one home of that
rule, and the chord's slots are the stretches between consecutive
crossings, so no point on a chord is ever located.

**Codes are out of range** (ruled 2026-10-05). The regions are
:math:`0, \dots, n - 1`; the inner exterior (a hollow body's cavity, the
slab's side below :math:`r_0`) is :math:`n` and the outer exterior
:math:`n + 1`, so indexing a per-region table with an exterior raises
instead of reading a material. A line lying in a surface carries that
breakpoint's index in :attr:`Chord.interface`, and :math:`n + 1`, out of
range for the breakpoints, where it lies in none.

**Lengths are computed without cancellation.** A half-chord as
:math:`\sqrt{(r - b)(r + b)}`; a shell's segment as
:math:`(r_{k+1} - r_k)(r_{k+1} + r_k)/(h_{k+1} + h_k)`, never as a
difference of half-chords or of crossing parameters; :math:`|P\Omega|`
with scaling.

**Point location (ruled 2026-10-05).** A bare orbit coordinate is located
inner-owns on the closed domain :math:`[r_0, r_n]`: region :math:`j` is
:math:`(r_j, r_{j+1}]`, region 0 is :math:`[r_0, r_1]`.

The design and its rulings: ``.claude/plans/characteristic_reference_architecture.md``.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields, replace
from functools import cached_property
from typing import TYPE_CHECKING, assert_never

import numpy as np

from orpheus.geometry.chart import AxialImage, Chart, RadialImage, half_chord
from orpheus.geometry.line import Line
from orpheus.geometry.transformation import RigidMotion
from orpheus.numerics.content import ContentIdentity
from orpheus.numerics.scalars import canonical_reals

if TYPE_CHECKING:
    from orpheus.geometry.structured_geometry import StructuredGeometry

__all__ = ["Chord", "ConcentricPartition", "Crossings", "Rays", "Transits", "region_half_lengths"]


@dataclass(frozen=True, eq=False)
class ConcentricPartition(ContentIdentity):
    r"""The level sets :math:`c = r_0 < r_1 < \dots < r_n` of a chart, posed in space.

    **Identity is content** (:class:`~orpheus.numerics.content.ContentIdentity`):
    the chart, the breakpoints and the pose.

    Attributes
    ----------
    chart:
        The chart whose orbit coordinate :math:`c` the breakpoints are values of.
    breakpoints:
        :math:`r_0 < \dots < r_n`, at least two; :math:`r_0 \ge 0` where the
        chart's group acts on the kept space (:math:`r_0 = 0`: solid; the
        centre or axis is a singular stratum, never a surface).
    pose:
        The rigid motion carrying the canonical frame to space; the identity
        by default.
    """

    chart: Chart
    breakpoints: tuple[float, ...]
    pose: RigidMotion = field(default_factory=lambda: RigidMotion.identity(3))

    def __post_init__(self) -> None:
        r = tuple(canonical_reals(np.asarray(self.breakpoints), "the breakpoints").tolist())  # ``-0.0`` folded: one memo key
        if len(r) < 2 or not all(np.isfinite(r)) or any(b <= a for a, b in zip(r, r[1:])):
            raise ValueError(f"breakpoints are finite and strictly increasing, at least two; got {r!r}")
        if self.chart.acts_on_kept_space and r[0] < 0.0:
            raise ValueError(f"a cylinder's or sphere's breakpoints start at r_0 >= 0; got {r[0]!r}")
        if self.pose.dimension != 3:
            raise ValueError(f"a partition is posed in R^3; the pose acts on R^{self.pose.dimension}")
        object.__setattr__(self, "breakpoints", r)

    @classmethod
    def of(cls, geometry: "StructuredGeometry", pose: RigidMotion | None = None) -> "ConcentricPartition":
        r"""The partition of a :class:`~orpheus.geometry.structured_geometry.StructuredGeometry`.

        Reads its coordinate system and breakpoints only; the materials and
        the boundary laws are the consumer's.
        """
        return cls(Chart(geometry.coord), geometry.breakpoints, pose or RigidMotion.identity(3))

    @property
    def n_regions(self) -> int:
        """The number :math:`n` of regions between the breakpoints."""
        return len(self.breakpoints) - 1

    @property
    def inner_exterior(self) -> int:
        r"""The code :math:`n` of the inner exterior: a hollow body's cavity, the slab's side below :math:`r_0`."""
        return self.n_regions

    @property
    def outer_exterior(self) -> int:
        r"""The code :math:`n + 1` of the outer exterior, above :math:`r_n`."""
        return self.n_regions + 1

    @property
    def no_interface(self) -> int:
        r"""The interface code of a line lying in no surface: :math:`n + 1`, out of range for the breakpoints."""
        return len(self.breakpoints)

    @property
    def _r(self) -> np.ndarray:
        return np.asarray(self.breakpoints)

    def region_containing(self, orbit_coordinate: np.ndarray) -> np.ndarray:
        r"""The region of each orbit-coordinate value, inner-owns on :math:`[r_0, r_n]`.

        Region :math:`j` is :math:`(r_j, r_{j+1}]`, region 0 is
        :math:`[r_0, r_1]`; below :math:`r_0` is :attr:`inner_exterior`,
        above :math:`r_n` :attr:`outer_exterior`. A non-finite value, and a
        negative one on a chart whose coordinate is a distance, is no orbit
        coordinate and is refused. For points a producer placed (quadrature
        nodes), the producer's record of their region is preferred to
        re-locating them (ruled 2026-10-05).
        """
        c = np.asarray(orbit_coordinate, dtype=float)
        if not np.all(np.isfinite(c)):
            raise ValueError("a non-finite orbit coordinate is in no region and in neither exterior; it is refused")
        if self.chart.acts_on_kept_space and np.any(c < 0.0):
            raise ValueError("an orbit coordinate that is a distance is never negative; it is refused")
        r = self._r
        region = np.searchsorted(r, c, side="left") - 1
        region = np.where(c == r[0], 0, region)
        region = np.where(c < r[0], self.inner_exterior, region)
        return np.where(c > r[-1], self.outer_exterior, region)

    def region_at(self, points: np.ndarray) -> np.ndarray:
        r"""The region of points ``(..., 3)`` in space (moved into the canonical frame by the pose's inverse)."""
        canonical = self.pose.inverse().on_points(points)
        return self.region_containing(self.chart.orbit_coordinate(canonical))

    def refined_at(self, orbit_coordinate: float) -> "ConcentricPartition":
        r"""The same partition, posed alike, with the level ``orbit_coordinate`` among its breakpoints; itself where it is one.

        The level lies on :math:`[r_0, r_n]`, so the refinement keeps the
        ends and its regions coarsen onto these (:meth:`coarsening`). A line
        chorded through it that crosses the level does so at a slot end; a
        line tangent to it (impact parameter :math:`c`) only touches it, at
        its closest approach.
        """
        c = float(orbit_coordinate)
        if not self.breakpoints[0] <= c <= self.breakpoints[-1]:
            raise ValueError(f"a refinement keeps the ends: the level {c!r} lies outside [{self.breakpoints[0]!r}, {self.breakpoints[-1]!r}]")
        if c in self.breakpoints:
            return self
        return replace(self, breakpoints=tuple(sorted((*self.breakpoints, c))))

    def coarsening(self, coarse: "ConcentricPartition") -> np.ndarray:
        r"""The region of ``coarse`` holding each region of this partition, ``(n,)``: the map a refinement carries onto what it refines.

        Refused unless this partition refines ``coarse``: the same chart, the
        same ends, and every breakpoint of ``coarse`` among its own.
        """
        if coarse.chart != self.chart or (coarse.breakpoints[0], coarse.breakpoints[-1]) != (
            self.breakpoints[0], self.breakpoints[-1]) or not set(coarse.breakpoints) <= set(self.breakpoints):
            raise ValueError(
                f"a partition coarsens only onto one it refines, on the same chart with the same ends; got "
                f"{self.chart.coord} {self.breakpoints!r} onto {coarse.chart.coord} {coarse.breakpoints!r}"
            )
        return np.searchsorted(coarse._r, self._r[:-1], side="right") - 1

    def rays_from(self, orbit_coordinate: np.ndarray, directions: np.ndarray) -> "Rays":
        r"""The lines through the points :math:`c\,\hat e_x` of the canonical frame along ``directions``, read from the points (:class:`Rays`).

        ``orbit_coordinate`` is :math:`c`, ``(...,)``, broadcast against the
        unit ``directions`` ``(..., 3)`` of the canonical frame; the pose
        carries both to space. Where the chart's group acts on the kept space
        the chord is formed from each line's exact level at the point,
        :math:`(c, c|\Omega_x|/|P\Omega|)` (#590), and the point's parameter
        from the same level, on the side of the closest approach its
        direction faces (past it where :math:`\Omega_x > 0`): a point on a
        breakpoint sits exactly on that crossing. On the slab the point's
        parameter is formed as the crossings' are, :math:`(c - c_{\rm foot})/\Omega_x`.
        """
        omega = np.asarray(directions, dtype=float)
        c = np.broadcast_to(np.asarray(orbit_coordinate, dtype=float), omega.shape[:-1])
        canonical = np.stack([c, np.zeros_like(c), np.zeros_like(c)], axis=-1)
        lines = Line.through(self.pose.on_points(canonical), self.pose.on_directions(omega))
        side = np.sign(omega[..., 0])
        if not self.chart.acts_on_kept_space:
            chord = self.chord(lines)
            return Rays(chord=chord, start=chord.image.parameters_at(c[..., None], side[..., None])[..., 0], side=side, orbit_coordinate=c)
        speed = self.chart.projected_speed(omega)
        level = (c, c * np.abs(omega[..., 0]) / np.where(speed > 0.0, speed, 1.0))
        chord = self.chord(lines, level)
        return Rays(chord=chord, start=chord.image.parameters_at(c[..., None], side[..., None])[..., 0], side=side, orbit_coordinate=c)

    def chord(self, line: Line, level: tuple[np.ndarray, np.ndarray] | None = None) -> "Chord":
        r"""The chord of each line through the partition (see the module docstring).

        Solved in the canonical frame; every parameter is then the caller's:
        a rigid motion preserves length along a line, so the canonical and
        the caller's parameters differ by the caller's parameter of the
        canonical foot's image.

        ``level`` is a radius and each line's exact half-chord there, both
        ``(...,)``, for a caller that knows them better than the line's
        stored impact parameter does (:meth:`RadialImage.at_level`; #590).
        They are orbit quantities, unchanged by the pose. Refused on a chart
        whose group fixes the kept space: an axial image's parameters are
        formed without a half-chord.
        """
        canonical = line.moved_by(self.pose.inverse())
        image = self.chart.image(canonical)
        match image, level:
            case RadialImage(), None:
                crossings, slot_length = _radial_chord(image, self._r)
            case RadialImage(), (radius, half_chord):
                image = image.at_level(radius, half_chord)
                crossings, slot_length = _radial_chord(image, self._r)
            case AxialImage(), None:
                crossings, slot_length = _axial_chord(image, self._r)
            case AxialImage(), _:
                raise ValueError("a level is a radius and a half-chord; a slab's lines have neither")
        shift = line.parameter_of(self.pose.on_points(canonical.foot))
        return Chord(
            partition=self,
            line=line,
            image=image.shifted(shift),
            crossings=replace(crossings, parameter=crossings.parameter + shift[..., None]),
            traversed_length=slot_length,
        )


@dataclass(frozen=True, eq=False)
class Crossings:
    r"""The crossings of a batch of lines with a partition's surfaces, in order along each line.

    Every array is ``(..., m)``, one column per potential crossing in the
    order a line meets them; ``present`` masks the ones it makes.

    Attributes
    ----------
    parameter:
        The line parameter :math:`t` of each crossing (for an absent one,
        where it would be: the closest approach on the cylinder and the
        sphere).
    breakpoint:
        The index :math:`k` of the surface :math:`c = r_k`.
    sense:
        +1 where the orbit coordinate increases through the crossing
        (outward on the cylinder and the sphere), -1 where it decreases.
    present:
        Whether the line makes the crossing.
    n_regions:
        The partition's :math:`n`, which fixes the exterior codes.
    """

    parameter: np.ndarray
    breakpoint: np.ndarray
    sense: np.ndarray
    present: np.ndarray
    n_regions: int

    @property
    def region_entered(self) -> np.ndarray:
        r"""The region each crossing enters, from the crossing order (:eq:`geometry-crossing-order`).

        Increasing through :math:`c = r_k` enters region :math:`k`, or the
        outer exterior :math:`n + 1` at :math:`k = n`; decreasing enters
        region :math:`k - 1`, or the inner exterior :math:`n` at
        :math:`k = 0`. The one home of the ruled rule.
        """
        k, n = self.breakpoint, self.n_regions
        rising = np.where(k < n, k, n + 1)
        falling = np.where(k >= 1, k - 1, n)
        return np.where(self.sense > 0, rising, falling)


@dataclass(frozen=True, eq=False)
class Chord:
    r"""The chord of a batch of lines through a :class:`ConcentricPartition`.

    The chord is a fixed sequence of **slots**, the stretches between
    consecutive potential crossings, each in the region the first of them
    enters; a slot the line does not traverse has length 0. On the cylinder
    and the sphere the slots are the regions inbound
    (:math:`n - 1, \dots, 0`), the inner exterior (a hollow body's cavity),
    then the regions outbound (:math:`0, \dots, n - 1`); the region holding
    the closest approach is traversed in its inbound and outbound slots,
    split at :math:`t^*`. On the slab they are the regions in the order the
    line meets them.

    Every parameter is measured on the caller's :attr:`line`, from its foot.

    Attributes
    ----------
    partition:
        The partition.
    line:
        The caller's lines.
    image:
        The lines' image in the orbit space (:class:`RadialImage` or
        :class:`AxialImage`), on the caller's parameters.
    crossings:
        The :class:`Crossings`.
    traversed_length:
        The cancellation-free length of each slot a non-parallel line
        traverses, ``(..., S)``; read through :attr:`slot_length`.
    """

    partition: ConcentricPartition
    line: Line
    image: RadialImage | AxialImage
    crossings: Crossings
    traversed_length: np.ndarray

    @property
    def parallel(self) -> np.ndarray:
        r"""Where :math:`|P\Omega| = 0`: the line keeps one orbit coordinate and never crosses."""
        return self.image.parallel

    @property
    def projected_speed(self) -> np.ndarray:
        r""":math:`|P\Omega|`, ``(...,)``."""
        return self.image.speed

    @property
    def slot_region(self) -> np.ndarray:
        """The region of each slot: the region its opening crossing enters, ``(..., S)``."""
        return self.crossings.region_entered[..., :-1]

    @property
    def slot_start(self) -> np.ndarray:
        """The parameter where each slot begins, ``(..., S)``; :math:`-\\infty` on a parallel line."""
        return np.where(self.parallel[..., None], -np.inf, self.crossings.parameter[..., :-1])

    @property
    def interface(self) -> np.ndarray:
        r"""For a parallel line lying in :math:`c = r_k`, :math:`k`; :attr:`ConcentricPartition.no_interface` otherwise.

        Decided by exact equality of the line's constant orbit coordinate
        with a breakpoint; a line rounding a hair off a surface (an axial
        line through a point computed on :math:`r = 1`; a rotated pose) is
        inside one region, which is measure-zero for every integral.
        """
        partition = self.partition
        on = self.image.coordinate_of_parallel[..., None] == partition._r
        if partition.chart.acts_on_kept_space and partition.breakpoints[0] == 0.0:
            on[..., 0] = False                              # a solid body's centre or axis is a stratum
        hit = on.any(axis=-1) & self.parallel
        return np.where(hit, np.argmax(on, axis=-1), partition.no_interface)

    @property
    def slot_length(self) -> np.ndarray:
        r"""The 3-D length of each slot, ``(..., S)``, without cancellation.

        A parallel line inside one region has one slot of infinite length
        in that region; one lying in a surface has every slot 0 (it is
        assigned to neither side; ruled 2026-10-05).
        """
        partition = self.partition
        home = partition.region_containing(self.image.coordinate_of_parallel)
        inside = self.parallel & (self.interface == partition.no_interface) & (home != partition.outer_exterior)
        candidate = self.slot_region == home[..., None]
        first = candidate & (np.cumsum(candidate, axis=-1) == 1)
        infinite = inside[..., None] & first
        return np.where(self.parallel[..., None], np.where(infinite, np.inf, 0.0), self.traversed_length)

    def lengths_beyond(self, start: np.ndarray) -> np.ndarray:
        r"""The length of each slot beyond the parameter ``start`` (a half-line), ``(..., S)``: :meth:`lengths_between` to :math:`+\infty`."""
        return self.lengths_between(start, np.inf)

    def lengths_between(self, start: np.ndarray | float, stop: np.ndarray | float) -> np.ndarray:
        r"""The length of each slot between the parameters ``start`` and ``stop``, each ``(...,)``, ``(..., S)``.

        A slot wholly inside :math:`[t_0, t_1]` keeps its cancellation-free
        length; only a slot containing an end is cut, as
        :math:`\min(t_{\rm end}, t_1) - \max(t_{\rm begin}, t_0)`, never above its
        own length; an empty interval (:math:`t_1 \le t_0`) has every slot 0.
        A parallel line keeps one orbit coordinate, so its home slot holds
        :math:`t_1 - t_0` (infinite to :math:`t_1 = +\infty`) and every other 0.
        """
        t0 = np.asarray(start, dtype=float)[..., None]
        t1 = np.asarray(stop, dtype=float)[..., None]
        begin = self.crossings.parameter[..., :-1]
        end = self.crossings.parameter[..., 1:]
        with np.errstate(invalid="ignore"):                  # inf - inf on an unbounded interval of a slot never reached
            cut = np.clip(np.minimum(end, t1) - np.maximum(begin, t0), 0.0, None)
        whole = (begin >= t0) & (end <= t1)
        between = np.where(whole, self.traversed_length, np.minimum(np.nan_to_num(cut, nan=0.0), self.traversed_length))
        home = np.isinf(self.slot_length)
        parallel = np.where(home, np.clip(t1 - t0, 0.0, None), 0.0)
        return np.where(self.parallel[..., None], parallel, between)

    def orbit_coordinate_at(self, t: np.ndarray) -> np.ndarray:
        r"""The orbit coordinate at the caller's parameters ``t`` ``(..., q)``, ``(..., q)``."""
        return self.image.orbit_coordinate_at(t)

    def rays_at(self, crossing: np.ndarray) -> "Rays":
        r"""The lines read from one of their crossings, the column ``crossing`` ``(...,)`` of :attr:`crossings` (:class:`Rays`).

        Each point sits exactly on its surface: its parameter, its orbit
        coordinate :math:`r_k` and its side are the crossing's own, and on the
        cylinder and the sphere the image is re-levelled at
        :math:`(r_k, y_k)`, the crossing's half-chord formed from the line's
        level, so :meth:`Rays.orbit_coordinate_at` reads :math:`r_k` at
        distance 0 and forms every distance from the wall, never through the
        parameter. Refused at a crossing the line does not make. Several
        points on one line are read from the chord of its lines repeated
        (:meth:`take`).
        """
        column = np.asarray(crossing)[..., None]
        if not np.all(np.take_along_axis(self.crossings.present, column, axis=-1)):
            raise ValueError("a line is read from a crossing it makes; this one it does not")
        start, breakpoint, sense = (np.take_along_axis(a, column, axis=-1)[..., 0] for a in (
            self.crossings.parameter, self.crossings.breakpoint, self.crossings.sense))
        radius = self.partition._r[breakpoint]
        match self.image:
            case RadialImage() as image:
                chord = replace(self, image=image.at_level(radius, image.half_chord_at(radius[..., None])[..., 0]))
            case AxialImage():
                chord = self
            case unreachable:
                assert_never(unreachable)
        return Rays(chord=chord, start=start, side=sense.astype(float), orbit_coordinate=radius)

    def rays_at_closest_approach(self) -> "Rays":
        r"""The lines read from their closest approach :math:`t^*` (:class:`Rays`).

        The point where a line turns in the orbit space, where a slot ends
        without a crossing: its orbit coordinate is :math:`b`, its side 0, and
        the image is levelled at :math:`(b, 0)`, so
        :meth:`Rays.orbit_coordinate_at` reads :math:`\sqrt{(|P\Omega|\ell)^2 + b^2}`
        with no subtraction, and a point advanced from it
        (:meth:`Rays.advanced`) near :math:`t^*` keeps its relative precision.
        Advanced from a crossing at :math:`r_k` instead, the same point forms
        :math:`r^2 - b^2` from :math:`r_k^2`, cancelling where :math:`r \ll r_k`.
        Refused on the slab, whose lines have no closest approach.
        """
        match self.image:
            case RadialImage() as image:
                b = image.impact_parameter
                levelled = replace(self, image=image.at_level(b, np.zeros_like(b)))
                return Rays(chord=levelled, start=image.closest_approach, side=np.zeros_like(b), orbit_coordinate=b)
            case AxialImage():
                raise ValueError("a slab line has no closest approach: its orbit coordinate is affine in the parameter")
            case unreachable:
                assert_never(unreachable)

    def take(self, lines: np.ndarray) -> "Chord":
        r"""The chord of the lines at flat indices ``lines`` ``(...)`` of the batch, in C order: a batch of ``lines``' shape.

        A line may be taken several times, to read several points on it
        (:meth:`rays_at`, :meth:`rays_at_closest_approach`). Every numeric
        field of the line, its image and the crossings whose rank reaches the
        batch's is the batch's and is gathered (a scalar of an unbatched line
        included), and refused if it does not lead with the batch's shape; a
        lower-rank field is shared, like the partition.
        """
        lines = np.asarray(lines, dtype=int)
        batch = self.line.direction.ndim - 1
        shape = self.line.direction.shape[:-1]

        def gathered(a: np.ndarray) -> np.ndarray:
            if a.shape[:batch] != shape:
                raise ValueError(f"a chord's field of rank >= its batch's leads with the batch {shape}; got {a.shape}")
            return a.reshape(-1, *a.shape[batch:])[lines]

        def taken(record):
            numbers = ((f.name, getattr(record, f.name)) for f in fields(record))
            return replace(record, **{
                name: gathered(np.asarray(a)) for name, a in numbers if isinstance(a, (np.ndarray, np.generic, float)) and np.ndim(a) >= batch
            })

        return replace(
            self,
            line=taken(self.line),
            image=taken(self.image),
            crossings=taken(self.crossings),
            traversed_length=gathered(self.traversed_length),
        )

    def line_integral(self, per_region: np.ndarray, lengths: np.ndarray) -> np.ndarray:
        r"""The integral :math:`\sum_{\rm slot} f_{\rm region}\,\ell_{\rm slot}` of a per-region field over slot lengths, ``lengths.shape[:-1]``.

        ``per_region`` holds one value per region, ``(n,)``, refused at any
        other shape (an extra value would fill the cavity, a missing one
        empty the outermost region). The exteriors are void, and a slot where
        the field is 0 adds 0 even at infinite length. ``lengths`` are slot
        lengths of this chord, ``(..., S)`` or ``(..., k, S)`` with ``k``
        runs per line (:attr:`slot_length`, :meth:`lengths_between`, a
        transit of them).
        """
        f = np.asarray(per_region, dtype=float)
        n = self.partition.n_regions
        if f.shape != (n,):
            raise ValueError(f"a per-region field holds one value per region, shape ({n},); got {f.shape}")
        region = self.slot_region
        value = np.where(region < n, f[np.minimum(region, n - 1)], 0.0)       # (..., S)
        length = np.asarray(lengths, dtype=float)
        runs = length.ndim - value.ndim
        value = np.expand_dims(value, axis=tuple(range(-1 - runs, -1))) if runs > 0 else value
        value, length = np.broadcast_arrays(value, length)
        return np.multiply(length, value, out=np.zeros(length.shape), where=value != 0.0).sum(axis=-1)

    def cosine_at(self, radius: float) -> np.ndarray:
        r"""The cosine between each line and the normal of the surface :math:`c` = ``radius`` where it crosses it, ``(...,)``.

        On the cylinder and the sphere :math:`|P\Omega|\,h/r`, with
        :math:`h` the half-chord there (0 where the line does not cross); on
        the slab :math:`|\Omega_x|`.
        """
        match self.image:
            case RadialImage() as image:
                h = image.half_chord_at(np.array([float(radius)]))[..., 0]
                return image.speed * h / float(radius)
            case AxialImage() as image:
                return image.speed
            case unreachable:
                assert_never(unreachable)

    @property
    def transits(self) -> "Transits":
        r"""The maximal runs of each line inside the domain :math:`[r_0, r_n]`, at most two (:class:`Transits`).

        A transit is a maximal run of slots in which every traversed slot
        (length :math:`> 0`) lies in a region :math:`< n`; an untraversed
        slot never breaks a run, so only a traversed exterior slot (a
        hollow body's cavity) separates two transits. Its walls are the
        breakpoints of the crossings that open its first traversed slot
        and close its last one. A parallel line, inside a region or in a
        surface, has none.
        """
        n = self.partition.n_regions
        slots = self.slot_region.shape[-1]
        traversed = (self.slot_length > 0.0) & ~self.parallel[..., None]
        run = np.cumsum(traversed & (self.slot_region >= n), axis=-1)        # traversed exterior slots so far
        inside = traversed & (self.slot_region < n)
        member = inside[..., None, :] & (run[..., None, :] == np.arange(_MAX_TRANSITS)[:, None])
        present = member.any(axis=-1)
        first = np.argmax(member, axis=-1)
        last = slots - 1 - np.argmax(member[..., ::-1], axis=-1)
        breakpoint = np.broadcast_to(self.crossings.breakpoint[..., None, :], (*present.shape, slots + 1))
        entry = np.take_along_axis(breakpoint, first[..., None], axis=-1)[..., 0]
        exit_ = np.take_along_axis(breakpoint, (last + 1)[..., None], axis=-1)[..., 0]
        absent_wall = self.partition.no_interface
        return Transits(
            first_slot=np.where(present, first, slots),
            stop_slot=np.where(present, last + 1, slots),
            entry_wall=np.where(present, entry, absent_wall),
            exit_wall=np.where(present, exit_, absent_wall),
            present=present,
        )


#: A concentric partition has one cavity at most, so a line has at most two transits.
_MAX_TRANSITS = 2


@dataclass(frozen=True, eq=False)
class Transits:
    r"""The transits of a batch of lines: the maximal runs inside the domain, in order along each line.

    Every array is ``(..., 2)``, one column per possible transit in the
    order the line meets them; ``present`` masks the ones it makes. A
    wall is named by its breakpoint index (:math:`0` or :math:`n`), never by
    position; an absent transit is the empty slot range at the end of the
    chord (``first_slot == stop_slot ==`` the slot count) and carries the
    wall code :math:`n + 1`, out of range for the breakpoints, so indexing
    a per-wall table with it raises.

    Attributes
    ----------
    first_slot:
        The first traversed slot of each transit.
    stop_slot:
        One past its last traversed slot.
    entry_wall:
        The breakpoint index of the crossing that opens ``first_slot``.
    exit_wall:
        The breakpoint index of the crossing that closes ``stop_slot - 1``.
    present:
        Whether the line makes the transit.
    """

    first_slot: np.ndarray
    stop_slot: np.ndarray
    entry_wall: np.ndarray
    exit_wall: np.ndarray
    present: np.ndarray



# ── the two images' chords ──────────────────────────────────────────────


def _radial_chord(image: RadialImage, r: np.ndarray) -> tuple[Crossings, np.ndarray]:
    n = len(r) - 1
    b = image.impact_parameter
    parallel = image.parallel
    h = image.half_chord_at(r)                            # from the image's level: exact there (#590)
    crossed = h > 0.0                                     # r_0 = 0 is never crossed: b >= 0
    with np.errstate(over="ignore"):                      # a subnormal |P Omega|: the obliquity is truly inf
        inv = 1.0 / np.where(parallel, 1.0, image.speed)

    shell = region_half_lengths(r, h)
    length = _lifted(np.concatenate([shell[..., ::-1], (2.0 * h[..., 0])[..., None], shell], axis=-1), inv)

    shape = (*b.shape, 2 * (n + 1))
    crossings = Crossings(
        parameter=image.parameter_at(np.concatenate([-h[..., ::-1], h], axis=-1)),
        breakpoint=np.broadcast_to(np.concatenate([np.arange(n, -1, -1), np.arange(n + 1)]), shape),
        sense=np.broadcast_to(np.concatenate([-np.ones(n + 1, dtype=int), np.ones(n + 1, dtype=int)]), shape),
        present=np.concatenate([crossed[..., ::-1], crossed], axis=-1),
        n_regions=n,
    )
    return crossings, length


def region_half_lengths(breakpoints: np.ndarray, half_chords: np.ndarray) -> np.ndarray:
    r"""The orbit-space length of each region on one side of a line's closest approach, ``(..., n)``.

    From the half-chords :math:`h_k` at the breakpoints :math:`r_k`,
    ``(..., n + 1)``, 0 where the line does not cross: region :math:`k` is
    :math:`(r_{k+1} - r_k)(r_{k+1} + r_k)/(h_{k+1} + h_k)`, never a
    difference of half-chords; the region of closest approach holds its
    outer half-chord :math:`h_{k+1}`; a region the line does not reach, 0.
    The chord's radial slots are these, inbound and outbound; a caller with
    the half-chords of an impact parameter and no line reads them directly.
    """
    r = np.asarray(breakpoints, dtype=float)
    h = np.asarray(half_chords, dtype=float)
    crossed = h > 0.0
    denom = h[..., 1:] + h[..., :-1]
    shell = np.where(
        crossed[..., :-1],
        (r[1:] - r[:-1]) * (r[1:] + r[:-1]) / np.where(denom > 0.0, denom, 1.0),
        h[..., 1:],
    )
    return np.where(crossed[..., 1:], shell, 0.0)


def _axial_chord(image: AxialImage, r: np.ndarray) -> tuple[Crossings, np.ndarray]:
    n = len(r) - 1
    parallel = image.parallel
    rising = image.rate > 0.0
    order = np.where(rising[..., None], np.arange(n + 1), np.arange(n, -1, -1))
    crossings = Crossings(
        parameter=image.parameters_at(r[order]),
        breakpoint=order,
        sense=np.where(rising[..., None], 1, -1) * np.ones_like(order),
        present=np.broadcast_to(~parallel[..., None], order.shape).copy(),
        n_regions=n,
    )
    with np.errstate(over="ignore"):                      # a subnormal |Omega_x|: the widths are truly inf
        widths = (r[1:] - r[:-1]) / np.where(parallel, 1.0, image.speed)[..., None]
    region = crossings.region_entered[..., :-1]              # interior regions only: 0 .. n-1
    return crossings, np.take_along_axis(np.broadcast_to(widths, region.shape), region, axis=-1)


def _lifted(orbit_length: np.ndarray, obliquity: np.ndarray) -> np.ndarray:
    r"""Orbit-space lengths ``(..., m)`` times the obliquity ``(...,)``: a zero length stays zero.

    At a subnormal :math:`|P\Omega|` the obliquity overflows to infinity,
    which is the true 3-D length of a traversed slot; an untraversed one is
    still 0, never the NaN of :math:`0 \cdot \infty`.
    """
    with np.errstate(over="ignore", invalid="ignore"):
        return np.where(orbit_length > 0.0, orbit_length * obliquity[..., None], 0.0)


@dataclass(frozen=True, eq=False)
class Rays:
    r"""The chord of the lines through points, read from the points (:meth:`ConcentricPartition.rays_from`).

    Attributes
    ----------
    chord:
        The chord of the lines.
    start:
        Each point's parameter on its line, ``(...,)``.
    side:
        The side of its line's closest approach each point is on, ``(...,)``:
        :math:`\operatorname{sign}\Omega_x`, past it (+1), before it (-1),
        at it (0).
    orbit_coordinate:
        Each point's orbit coordinate :math:`c`, ``(...,)``.
    """

    chord: Chord
    start: np.ndarray
    side: np.ndarray
    orbit_coordinate: np.ndarray

    @cached_property
    def ahead(self) -> np.ndarray:
        """The length of each slot ahead of the point, ``(..., S)``: the run to the outer surface, a cavity crossed."""
        return self.chord.lengths_beyond(self.start)

    def lengths_to(self, distance: np.ndarray) -> np.ndarray:
        """The length of each slot between the point and the distance ``distance`` ahead of it, ``(..., S)``."""
        return self.chord.lengths_between(self.start, self.start + np.asarray(distance, dtype=float))

    @cached_property
    def _transit(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """The transit each point is in: (index, present, its slots ``(..., S)``).

        The first transit the line has not left at the point, its closing
        crossing at or ahead of the point's parameter. A point on a wall
        facing out of the domain is in the transit that ends there, at
        length 0: it leaves through that wall at once.
        """
        transits = self.chord.transits
        slot = np.arange(self.chord.slot_region.shape[-1])
        closing = np.take_along_axis(self.chord.crossings.parameter, transits.stop_slot, axis=-1)   # (..., 2)
        current = transits.present & (closing >= self.start[..., None])
        index = np.asarray(np.argmax(current, axis=-1))
        present = np.asarray(current.any(axis=-1))
        first = np.take_along_axis(transits.first_slot, index[..., None], axis=-1)
        stop = np.take_along_axis(transits.stop_slot, index[..., None], axis=-1)
        slots = (slot >= first) & (slot < stop) & present[..., None]
        return index, present, slots

    @cached_property
    def to_wall(self) -> np.ndarray:
        """The length of each slot from the point to the first wall ahead of it, ``(..., S)``; 0 where it leaves the domain at once."""
        return np.where(self._transit[2], self.ahead, 0.0)

    @cached_property
    def wall(self) -> np.ndarray:
        r"""The breakpoint index of the first wall ahead of each point, ``(...,)``: :math:`n` the outer surface, 0 a hollow body's cavity (or the slab's lower face); :attr:`ConcentricPartition.no_interface` where none."""
        index, present, _ = self._transit
        exit_wall = np.take_along_axis(self.chord.transits.exit_wall, index[..., None], axis=-1)[..., 0]
        return np.where(present, exit_wall, self.chord.partition.no_interface)

    def orbit_coordinate_at(self, distance: np.ndarray) -> np.ndarray:
        r"""The orbit coordinate at the distances ``distance`` ``(..., q)`` ahead of each point, ``(..., q)``.

        Formed from the point, never through the parameter, whose round trip
        rounds: on the cylinder and the sphere from its orbit-space position,
        :math:`s = \pm y_* + |P\Omega|\,\ell`; on the slab as
        :math:`c + \Omega_x \ell`. So it is exactly the point's coordinate at
        :math:`\ell = 0`, and a point on the outer surface never reads above it.
        """
        ell = np.asarray(distance, dtype=float)
        match self.chord.image:
            case RadialImage() as image:
                position = (self.side * image.level_half_chord)[..., None] + image.speed[..., None] * ell
                return half_chord(position, image.level_half_chord[..., None], image.level[..., None])
            case AxialImage() as image:
                return self.orbit_coordinate[..., None] + image.rate[..., None] * ell
            case unreachable:
                assert_never(unreachable)

    def advanced(self, distance: np.ndarray) -> "Rays":
        r"""The same lines read from the points ``distance`` ``(...,)`` further along (behind, where negative) (:class:`Rays`).

        Each new point is anchored at its own coordinate: on the cylinder and
        the sphere the image is re-levelled at the new point's radius and
        half-chord :math:`(c', |s'|)`, :math:`s' = \pm y_* + |P\Omega|\,d`, on
        its side of the closest approach; on the slab the coordinate is
        :math:`c + \Omega_x d`. Distances from the new point are then formed
        from it. One rounding of :math:`d` is spent here, so a point on a
        surface is read from the crossing (:meth:`Chord.rays_at`) and
        advanced from there only to a point inside a slot.
        """
        d = np.asarray(distance, dtype=float)
        match self.chord.image:
            case RadialImage() as image:
                position = self.side * image.level_half_chord + image.speed * d
                radius = half_chord(position, image.level_half_chord, image.level)
                chord = replace(self.chord, image=image.at_level(radius, np.abs(position)))
                return Rays(chord=chord, start=self.start + d, side=np.sign(position), orbit_coordinate=radius)
            case AxialImage() as image:
                return Rays(chord=self.chord, start=self.start + d, side=self.side, orbit_coordinate=self.orbit_coordinate + image.rate * d)
            case unreachable:
                assert_never(unreachable)
