r"""Oriented lines in :math:`\mathbb{R}^3`, held in Plücker coordinates.

A line is the set :math:`\{p + t\,\Omega\}` of a point :math:`p` and a unit
direction :math:`\Omega`; it does not depend on which of its points is
named. Its Plücker coordinates (the direction and the moment
:math:`m = p \times \Omega`) are that independence made explicit: moving the
base point along the line, :math:`p \to p + s\Omega`, leaves :math:`m`
unchanged, because :math:`\Omega \times \Omega = 0`. The point of the line
closest to the origin is the foot :math:`\Omega \times m`, its distance from
the origin :math:`|m|`, and a line's own parameter :math:`t` is measured
from its foot.

A ray is a line with a starting parameter (a backward characteristic is the
line with :math:`-\Omega`); it is not a separate type. A direction is a
point of the unit sphere, the manifold the angular quadratures'
ordinates live on (:class:`~orpheus.numerics.manifold.Sphere`); here it is
held as the array of its Cartesian components.

Every value is batch-shaped: a :class:`Line` holds ``(..., 3)`` arrays, one
line per leading index, because every consumer (a characteristic sweep over
positions and directions) asks about many lines at once.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from orpheus.geometry.transformation import RigidMotion

__all__ = ["Line"]

#: How far a direction may be from unit length before it is refused, in
#: units of the float spacing at 1: a direction computed from angles or
#: normalised once is unit to a few roundings; a vector that is not unit
#: is a caller error, never renormalised here.
_UNIT_ULPS = 8.0


@dataclass(frozen=True, eq=False)
class Line:
    r"""A batch of oriented lines :math:`\{p + t\,\Omega\}` in Plücker coordinates.

    Build from points with :meth:`through`. Two base points on one line
    give one line (equal moments, up to the rounding of the cross product).

    Attributes
    ----------
    direction:
        The unit directions :math:`\Omega`, ``(..., 3)``.
    moment:
        The moments :math:`m = p \times \Omega`, ``(..., 3)``.
    """

    direction: np.ndarray
    moment: np.ndarray

    def __post_init__(self) -> None:
        omega = np.asarray(self.direction, dtype=float)
        m = np.asarray(self.moment, dtype=float)
        if omega.shape[-1:] != (3,) or m.shape != omega.shape:
            raise ValueError(
                f"a line's direction and moment are (..., 3) arrays of one shape; "
                f"got {omega.shape} and {m.shape}"
            )
        departure = np.abs(np.sqrt(np.sum(omega * omega, axis=-1)) - 1.0)
        # Written as "all within", not "any beyond": a NaN departure fails every
        # comparison, so only the first spelling refuses a non-finite direction.
        if not np.all(departure <= _UNIT_ULPS * np.finfo(float).eps):
            raise ValueError(
                f"a line's direction is a finite unit vector; the largest departure from "
                f"unit length is {float(np.nanmax(np.where(np.isfinite(departure), departure, np.inf))):.3e} "
                f"(a direction is refused, never renormalised)"
            )
        if not np.all(np.isfinite(m)):
            raise ValueError("a line's moment is finite (its base point is a finite point)")
        object.__setattr__(self, "direction", omega)
        object.__setattr__(self, "moment", m)

    @classmethod
    def through(cls, points: np.ndarray, directions: np.ndarray) -> "Line":
        r"""The lines through ``points`` along ``directions``, both ``(..., 3)`` (broadcast)."""
        p = np.asarray(points, dtype=float)
        omega = np.asarray(directions, dtype=float)
        p, omega = np.broadcast_arrays(p, omega)
        if not np.all(np.isfinite(p)):
            raise ValueError("a line's base point is a finite point")
        if not np.all(np.isfinite(omega)):
            raise ValueError("a line's direction is a finite unit vector; got a non-finite component")
        return cls(direction=omega.copy(), moment=np.cross(p, omega))

    @property
    def shape(self) -> tuple[int, ...]:
        """The batch shape: one line per leading index."""
        return self.direction.shape[:-1]

    @property
    def foot(self) -> np.ndarray:
        r"""The point of each line closest to the origin, :math:`\Omega \times m`, ``(..., 3)``.

        The origin of the line's own parameter :math:`t`.
        """
        return np.cross(self.direction, self.moment)

    def parameter_of(self, points: np.ndarray) -> np.ndarray:
        r"""The parameter :math:`t = (x - \text{foot})\cdot\Omega` of points on the lines, ``(...,)``."""
        x = np.asarray(points, dtype=float)
        return np.sum((x - self.foot) * self.direction, axis=-1)

    def at(self, t: np.ndarray) -> np.ndarray:
        r"""The points :math:`\text{foot} + t\,\Omega`, ``(..., 3)``."""
        return self.foot + np.asarray(t, dtype=float)[..., None] * self.direction

    def reversed(self) -> "Line":
        r"""The same lines with the opposite orientation: :math:`(-\Omega, -m)`."""
        return Line(direction=-self.direction, moment=-self.moment)

    def moved_by(self, motion: RigidMotion) -> "Line":
        r"""The image of the lines under a rigid motion :math:`x \mapsto Qx + t`.

        The direction moves linearly (:meth:`RigidMotion.on_directions`), a
        point on the line affinely (:meth:`RigidMotion.on_points`); the
        image line passes through the image of the foot.

        The parameter is not carried over: :math:`t` is measured from the
        foot, the point closest to the origin, and the image's foot is the
        point closest to the origin again, not the image of the old foot.
        So ``moved_by(g).at(t)`` is not ``g.on_points(at(t))``; the two
        parameters differ by ``moved_by(g).parameter_of(g.on_points(foot))``
        (lengths along the line are preserved, so the difference is one
        shift per line).
        """
        return Line.through(motion.on_points(self.foot), motion.on_directions(self.direction))
