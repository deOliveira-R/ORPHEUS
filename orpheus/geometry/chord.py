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

from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING

import numpy as np

from orpheus.geometry.chart import AxialImage, Chart, RadialImage
from orpheus.geometry.line import Line
from orpheus.geometry.transformation import RigidMotion

if TYPE_CHECKING:
    from orpheus.geometry.structured_geometry import StructuredGeometry

__all__ = ["Chord", "ConcentricPartition", "Crossings"]


@dataclass(frozen=True, eq=False)
class ConcentricPartition:
    r"""The level sets :math:`c = r_0 < r_1 < \dots < r_n` of a chart, posed in space.

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
        r = tuple(float(v) for v in self.breakpoints)
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

    def chord(self, line: Line) -> "Chord":
        r"""The chord of each line through the partition (see the module docstring).

        Solved in the canonical frame; every parameter is then the caller's:
        a rigid motion preserves length along a line, so the canonical and
        the caller's parameters differ by the caller's parameter of the
        canonical foot's image.
        """
        canonical = line.moved_by(self.pose.inverse())
        image = self.chart.image(canonical)
        match image:
            case RadialImage():
                crossings, slot_length = _radial_chord(image, self._r)
            case AxialImage():
                crossings, slot_length = _axial_chord(image, self._r)
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
        r"""The length of each slot beyond the parameter ``start`` (a half-line), ``(..., S)``.

        A slot wholly beyond ``start`` keeps its cancellation-free length;
        only the slot containing ``start`` is cut, as the difference of its
        end and ``start``. A parallel line is unbounded both ways.
        """
        t0 = np.asarray(start, dtype=float)[..., None]
        begin = self.crossings.parameter[..., :-1]
        end = self.crossings.parameter[..., 1:]
        cut = np.clip(end - t0, 0.0, None)
        beyond = np.where(begin >= t0, self.traversed_length, np.minimum(cut, self.traversed_length))
        return np.where(self.parallel[..., None], self.slot_length, beyond)

    def orbit_coordinate_at(self, t: np.ndarray) -> np.ndarray:
        r"""The orbit coordinate at the caller's parameters ``t`` ``(..., q)``, ``(..., q)``."""
        return self.image.orbit_coordinate_at(t)



# ── the two images' chords ──────────────────────────────────────────────


def _radial_chord(image: RadialImage, r: np.ndarray) -> tuple[Crossings, np.ndarray]:
    n = len(r) - 1
    b = image.impact_parameter
    parallel = image.parallel
    crossed = (b[..., None] < r) & ~parallel[..., None]    # r_0 = 0 is never crossed: b >= 0
    h = np.sqrt(np.where(crossed, (r - b[..., None]) * (r + b[..., None]), 0.0))
    with np.errstate(over="ignore"):                      # a subnormal |P Omega|: the obliquity is truly inf
        inv = 1.0 / np.where(parallel, 1.0, image.speed)

    # Orbit-space length of region j on one side: cancellation-free; the
    # region of closest approach gets its half-chord; untraversed regions 0.
    denom = h[..., 1:] + h[..., :-1]
    shell = np.where(
        crossed[..., :-1],
        (r[1:] - r[:-1]) * (r[1:] + r[:-1]) / np.where(denom > 0.0, denom, 1.0),
        h[..., 1:],
    )
    shell = np.where(crossed[..., 1:], shell, 0.0)
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


def _axial_chord(image: AxialImage, r: np.ndarray) -> tuple[Crossings, np.ndarray]:
    n = len(r) - 1
    parallel = image.parallel
    rising = image.rate > 0.0
    rate = np.where(parallel, 1.0, image.rate)
    order = np.where(rising[..., None], np.arange(n + 1), np.arange(n, -1, -1))
    crossings = Crossings(
        parameter=(r[order] - image.foot_coordinate[..., None]) / rate[..., None],
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
