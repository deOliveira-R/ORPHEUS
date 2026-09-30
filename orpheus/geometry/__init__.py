"""The geometry layer: shapes, coordinate systems and boundary laws, with no discretisation.

* :class:`StructuredGeometry`: the 1-D layered geometry (a coordinate
  system, breakpoints, one material id per interval and one law per
  boundary point), a pure shape and boundary description with no cell
  counts. Reference solvers consume it directly.
* :class:`CoordSystem` and the volume and area formulas of each coordinate
  system.
* :mod:`~orpheus.geometry.boundary`: the boundary-condition tag
  :class:`BC` and the typed boundary laws. Boundaries are defined at the
  geometry.
* :mod:`~orpheus.geometry.transformation`: rigid motions and permutations.

The mesh is an overlay on the geometry and lives in its own package,
:mod:`orpheus.mesh`, which imports this one and never the reverse:
:meth:`~orpheus.mesh.structured.Mesh1D.from_geometry` discretises a
:class:`StructuredGeometry`.
"""

from .boundary import BC
from .coord import CoordSystem, MeasureCoordinate, compute_areas_1d, compute_volumes_1d, compute_volumes_2d
from .structured_geometry import StructuredGeometry
from .transformation import (
    NotAFinitePointGroupError,
    Permutation,
    RigidMotion,
    close_group,
)

__all__ = [
    "BC",
    "CoordSystem",
    "MeasureCoordinate",
    "NotAFinitePointGroupError",
    "Permutation",
    "RigidMotion",
    "StructuredGeometry",
    "close_group",
    "compute_areas_1d",
    "compute_volumes_1d",
    "compute_volumes_2d",
]
