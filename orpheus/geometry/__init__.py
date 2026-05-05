"""Geometry module for ORPHEUS reactor physics solvers.

Provides coordinate-system-aware mesh data structures and the
1-D :class:`StructuredGeometry` → :meth:`Mesh1D.from_geometry`
pipeline.

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
per region) and returns a discretized :class:`Mesh1D`. The
:meth:`StructuredGeometry.wigner_seitz_pin_cell` and
:meth:`StructuredGeometry.pwr_slab_half_cell` classmethods cover the
standard PWR-shaped 1-D geometries; for 2-D Cartesian pin meshes use
:func:`pwr_pin_2d`.
"""

from .coord import CoordSystem, compute_surfaces_1d, compute_volumes_1d, compute_volumes_2d
from .factories import pwr_pin_2d
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
    "compute_surfaces_1d",
    "compute_volumes_1d",
    "compute_volumes_2d",
    "pwr_pin_2d",
]
