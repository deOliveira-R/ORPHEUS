"""Geometry module for ORPHEUS reactor physics solvers.

Provides coordinate-system-aware mesh data structures and factories
for common reactor geometries.
"""

from geometry.coord import CoordSystem, compute_surfaces_1d, compute_volumes_1d, compute_volumes_2d
from geometry.factories import Zone, mesh1d_from_zones, pwr_pin_2d, pwr_pin_equivalent, pwr_slab_half_cell
from geometry.mesh import Mesh1D, Mesh2D

__all__ = [
    "CoordSystem",
    "Mesh1D",
    "Mesh2D",
    "Zone",
    "compute_surfaces_1d",
    "compute_volumes_1d",
    "compute_volumes_2d",
    "mesh1d_from_zones",
    "pwr_pin_2d",
    "pwr_pin_equivalent",
    "pwr_slab_half_cell",
]
