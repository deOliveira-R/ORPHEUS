"""Geometry module for ORPHEUS reactor physics solvers.

Provides coordinate-system-aware mesh data structures and factories
for common reactor geometries.

Architectural layers
--------------------

* **Geometry layer** — :class:`StructuredGeometry`, :class:`Region`.
  Pure shape + BC description. No cell counts, no discretization.
  Reference solvers consume this directly.
* **Mesh layer** — :class:`Mesh1D`, :class:`Mesh2D`,
  :class:`RegionMesh`. Discrete representations and the
  per-region discretization descriptor used to build them.
  Production solvers consume :class:`Mesh1D`.

The canonical bridge is :meth:`Mesh1D.from_geometry`, which takes a
:class:`StructuredGeometry` plus a tuple of :class:`RegionMesh` (one
per region) and returns a discretized :class:`Mesh1D`.
"""

from .coord import CoordSystem, compute_surfaces_1d, compute_volumes_1d, compute_volumes_2d
from .factories import (
    Zone,
    homogeneous_1d,
    mesh1d_from_zones,
    pwr_pin_2d,
    pwr_pin_equivalent,
    pwr_slab_half_cell,
    slab_fuel_moderator,
)
from .mesh import BC, Mesh1D, Mesh2D, RegionMesh
from .structured_geometry import Region, StructuredGeometry

__all__ = [
    "BC",
    "CoordSystem",
    "Mesh1D",
    "Mesh2D",
    "Region",
    "RegionMesh",
    "StructuredGeometry",
    "Zone",
    "compute_surfaces_1d",
    "compute_volumes_1d",
    "compute_volumes_2d",
    "homogeneous_1d",
    "mesh1d_from_zones",
    "pwr_pin_2d",
    "pwr_pin_equivalent",
    "pwr_slab_half_cell",
    "slab_fuel_moderator",
]
