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

from enum import Enum

import numpy as np


class CoordSystem(Enum):
    r"""Coordinate system identifier, and the measure of its position axis.

    The measure of an interval :math:`[a, b]` of positions is

    .. math::

        m(a, b) = c\,(b^{d} - a^{d}),

    with the exponent :math:`d` the dimension the position sweeps (1 for a
    slab, 2 for a cylinder, 3 for a sphere) and the constant :math:`c`
    (1, :math:`\pi`, :math:`\tfrac43\pi`): a slab's length per unit
    transverse area, a cylinder's area per unit height, a sphere's volume.
    :math:`T(r) = r^{d}` is the coordinate in which the measure is uniform,
    so equal-measure cells are equal steps of :math:`T`.
    """

    CARTESIAN = "cartesian"
    CYLINDRICAL = "cylindrical"
    SPHERICAL = "spherical"

    @property
    def measure_exponent(self) -> int:
        r"""The exponent :math:`d` of :math:`m(a, b) = c\,(b^d - a^d)`."""
        match self:
            case CoordSystem.CARTESIAN:
                return 1
            case CoordSystem.CYLINDRICAL:
                return 2
            case CoordSystem.SPHERICAL:
                return 3

    @property
    def measure_constant(self) -> float:
        r"""The constant :math:`c` of :math:`m(a, b) = c\,(b^d - a^d)`."""
        match self:
            case CoordSystem.CARTESIAN:
                return 1.0
            case CoordSystem.CYLINDRICAL:
                return np.pi
            case CoordSystem.SPHERICAL:
                return (4.0 / 3.0) * np.pi

    def interval_measure(self, a: float, b: float) -> float:
        r"""The measure :math:`m(a, b) = c\,(b^d - a^d)` of the interval :math:`[a, b]`.

        Evaluated on scalars, in the order the equal-volume subdivision has
        always used (``π · (b² − a²)``), so an equal-volume cell stored as
        ``m(a, b) / n`` keeps its bits (ERR-020). The array form over cell
        edges is :func:`compute_volumes_1d`.
        """
        d = self.measure_exponent
        return self.measure_constant * (b**d - a**d)


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
    match coord:
        case CoordSystem.CARTESIAN:
            return np.diff(edges)
        case CoordSystem.CYLINDRICAL:
            return np.pi * np.diff(edges**2)
        case CoordSystem.SPHERICAL:
            return (4.0 / 3.0) * np.pi * np.diff(edges**3)
        case _:
            raise ValueError(f"Unknown coordinate system: {coord}")


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
