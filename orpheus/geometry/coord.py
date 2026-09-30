"""Coordinate systems and their volume / surface formulas.

This module is the **single point** where coordinate-system dependence
lives.  All mesh classes delegate to these functions.

Supported coordinate systems
-----------------------------
* **Cartesian** -- flat geometry (slab, plate, box)
* **Cylindrical** -- annular geometry (pin cell, tube)
* **Spherical** -- shell geometry (pebble, sphere)

Volume formulas (1-D)
~~~~~~~~~~~~~~~~~~~~~
=========== ==========================================
Cartesian   :math:`V_i = x_{i+1} - x_i`
Cylindrical :math:`V_i = \\pi (r_{i+1}^2 - r_i^2)`
Spherical   :math:`V_i = \\tfrac{4}{3}\\pi (r_{i+1}^3 - r_i^3)`
=========== ==========================================

Area formulas (1-D)
~~~~~~~~~~~~~~~~~~~~~~
=========== ==========================================
Cartesian   :math:`A = 1` (per unit transverse area)
Cylindrical :math:`A = 2\\pi r` (per unit height)
Spherical   :math:`A = 4\\pi r^2`
=========== ==========================================
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import numpy as np


@dataclass(frozen=True)
class MeasureCoordinate:
    r"""The coordinate :math:`T(r) = r^{p}`, :math:`p \in \{1, 2, 3\}`, and its inverse.

    A measure on the position axis that is uniform in :math:`T` gives an
    interval :math:`[a, b]` the measure :math:`c\,(T(b) - T(a))`; cells of
    equal measure are equal steps of :math:`T`. Evaluated on arrays, so
    :math:`T` is numpy's correctly-rounded power (Python's scalar ``b**2``
    calls libm ``pow``, which is not correctly rounded in about 1 case in
    1000).
    """

    exponent: int

    def __post_init__(self) -> None:
        if self.exponent not in (1, 2, 3):
            raise ValueError(
                f"a measure coordinate is r, r**2 or r**3; got the exponent {self.exponent!r}"
            )

    def __call__(self, r: np.ndarray) -> np.ndarray:
        return np.asarray(r, dtype=float) ** self.exponent

    def inverse(self, t: np.ndarray) -> np.ndarray:
        t = np.asarray(t, dtype=float)
        match self.exponent:
            case 1:
                return t.copy()
            case 2:
                return np.sqrt(t)
            case _:
                return np.cbrt(t)


class CoordSystem(Enum):
    r"""Coordinate system identifier, and the measure of its position axis.

    The measure of the cells between edges :math:`r_0 < r_1 < \dots` is

    .. math::

        m_j = c\,\bigl(T(r_{j+1}) - T(r_j)\bigr), \qquad T(r) = r^{d},

    with :math:`d` the dimension the position sweeps (1 for a slab, 2 for a
    cylinder, 3 for a sphere) and :math:`c` (1, :math:`\pi`,
    :math:`\tfrac43\pi`): a slab's length per unit transverse area, a
    cylinder's area per unit height, a sphere's volume. :math:`T` is the
    :class:`MeasureCoordinate` in which the measure is uniform. This is the
    one definition of the measure: :meth:`measure` evaluates it, and an
    interval's measure is its one-cell case.
    """

    CARTESIAN = "cartesian"
    CYLINDRICAL = "cylindrical"
    SPHERICAL = "spherical"

    @property
    def measure_coordinate(self) -> MeasureCoordinate:
        r"""The coordinate :math:`T(r) = r^{d}` in which this system's measure is uniform."""
        match self:
            case CoordSystem.CARTESIAN:
                return MeasureCoordinate(1)
            case CoordSystem.CYLINDRICAL:
                return MeasureCoordinate(2)
            case CoordSystem.SPHERICAL:
                return MeasureCoordinate(3)

    @property
    def measure_constant(self) -> float:
        r"""The constant :math:`c` of :math:`m_j = c\,(T(r_{j+1}) - T(r_j))`."""
        match self:
            case CoordSystem.CARTESIAN:
                return 1.0
            case CoordSystem.CYLINDRICAL:
                return np.pi
            case CoordSystem.SPHERICAL:
                return (4.0 / 3.0) * np.pi

    def boundary_points(self, r_0: float, r_R: float) -> tuple[float, ...]:
        r"""The boundary points of the interval :math:`[r_0, r_R]` in this system, inner first.

        Both ends on a slab, and on a hollow cylinder or sphere
        (:math:`r_0 > 0`); only :math:`r_R` on a solid one, whose centre
        :math:`r = 0` is an interior point of the region and carries no law.
        """
        if self is CoordSystem.CARTESIAN or r_0 > 0.0:
            return (r_0, r_R)
        return (r_R,)

    def measure(self, edges: np.ndarray) -> np.ndarray:
        r"""The measures :math:`c\,(T(r_{j+1}) - T(r_j))` of the cells between ``edges``."""
        return self.measure_constant * np.diff(self.measure_coordinate(edges))


# ── 1-D formulas ─────────────────────────────────────────────────────

def compute_volumes_1d(coord: CoordSystem, edges: np.ndarray) -> np.ndarray:
    """Cell volumes from 1-D edge positions.

    Parameters
    ----------
    coord : CoordSystem
        Coordinate system.
    edges : ndarray, shape (N+1,)
        Monotonically increasing edge positions.

    Returns
    -------
    ndarray, shape (N,)
        Volume of each cell.
    """
    return coord.measure(edges)


def compute_areas_1d(coord: CoordSystem, edges: np.ndarray) -> np.ndarray:
    """Face areas at each 1-D edge position.

    Parameters
    ----------
    coord : CoordSystem
        Coordinate system.
    edges : ndarray, shape (N+1,)
        Edge positions.

    Returns
    -------
    ndarray, shape (N+1,)
        Face area at every edge.
    """
    match coord:
        case CoordSystem.CARTESIAN:
            return np.ones_like(edges)
        case CoordSystem.CYLINDRICAL:
            return 2.0 * np.pi * edges
        case CoordSystem.SPHERICAL:
            return 4.0 * np.pi * edges**2
        case _:
            raise ValueError(f"Unknown coordinate system: {coord}")


# ── 2-D formulas ─────────────────────────────────────────────────────

def compute_volumes_2d(
    coord: CoordSystem,
    edges_x: np.ndarray,
    edges_y: np.ndarray,
) -> np.ndarray:
    """Cell volumes from 2-D edge positions.

    Parameters
    ----------
    coord : CoordSystem
        ``CARTESIAN`` for (x, y) or ``CYLINDRICAL`` for (r, z).
    edges_x : ndarray, shape (Nx+1,)
        Edge positions in x (or radial) direction.
    edges_y : ndarray, shape (Ny+1,)
        Edge positions in y (or axial) direction.

    Returns
    -------
    ndarray, shape (Nx, Ny)
        Volume of each cell.
    """
    match coord:
        case CoordSystem.CARTESIAN:
            dx = np.diff(edges_x)
            dy = np.diff(edges_y)
            return dx[:, np.newaxis] * dy[np.newaxis, :]
        case CoordSystem.CYLINDRICAL:
            # r-z geometry: V = pi * (r_out^2 - r_in^2) * dz
            dr2 = np.diff(edges_x**2)  # (Nr,)
            dz = np.diff(edges_y)       # (Nz,)
            return np.pi * dr2[:, np.newaxis] * dz[np.newaxis, :]
        case _:
            raise ValueError(
                f"2-D volumes not defined for {coord}; "
                f"use CARTESIAN or CYLINDRICAL"
            )
