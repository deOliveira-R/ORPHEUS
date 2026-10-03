r"""The geometric coordinate of the reference questions (#405 P1 step 8).

:class:`GeometryExtent` is a DIRECTION in the system's parameter space: the
width of one interval of a :class:`~orpheus.geometry.structured_geometry.StructuredGeometry`,
in cm, with every interval outside it translated outward (the user's ruling
of 2026-10-02). On a one-interval body it is the width of that interval,
the critical radius of a solid sphere or the full width of a bare slab in the
Sood benchmarks; on a reflected body it is the core grown under a reflector
of fixed thickness. An infinite medium has no geometry and so
no extent. It is the opaque,
hashable parameter key an :class:`~orpheus.numerics.question.Eigen` holds for
a critical-extent question, resolved by the reference specification against
its geometry.

The chart (the zero of the coordinate, its physical value, its admissible
range) is read by the mode law of #529, not here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from orpheus.numerics.content import ContentIdentity
from orpheus.numerics.scalars import parse_integer

if TYPE_CHECKING:
    from orpheus.geometry.structured_geometry import StructuredGeometry

__all__ = ["GeometryExtent"]


@dataclass(frozen=True, eq=False)
class GeometryExtent(ContentIdentity):
    """The width of interval ``interval``, in cm, every interval outside it translated outward.

    ``interval`` is a non-negative index from the innermost interval. A
    negative index is refused rather than read from the outside: ``-1`` would
    be a second spelling of one coordinate, whose meaning moves with the
    interval count.
    """

    interval: int

    def __post_init__(self) -> None:
        index = parse_integer(self.interval, "GeometryExtent", "the interval")
        if index < 0:
            raise ValueError(f"GeometryExtent: the interval index is non-negative, got {index}")
        object.__setattr__(self, "interval", index)

    def resolve(self, geometry: "StructuredGeometry") -> "GeometryExtent":
        """This coordinate, once ``geometry`` is shown to have the interval it names."""
        n_intervals = len(geometry.intervals)
        if self.interval >= n_intervals:
            raise ValueError(
                f"GeometryExtent: interval {self.interval} does not exist; "
                f"the geometry has {n_intervals} interval(s)"
            )
        return self
