r"""The mesh: the discretisation overlay on a geometry.

Posing is a chain of overlays (materials, then geometry, then mesh), and each
keeps or lowers the symmetry group. The mesh is the overlay that divides a
geometry's intervals into cells, so this package imports
:mod:`orpheus.geometry` (shapes, coordinate systems, boundary laws) and
:mod:`orpheus.numerics`, and neither of those imports it
(``tests/gates/test_layer_imports.py`` enforces both directions). Binding cross
sections to the cells is the transport layer's
:class:`~orpheus.transport.mesh.material_mesh.MaterialMesh`, one layer up.

* :mod:`~orpheus.mesh.structured`: the structured meshes :class:`Mesh1D` and
  :class:`Mesh2D`, and the per-region discretisation descriptor
  :class:`RegionMesh`.
* :mod:`~orpheus.mesh.factories`: the 2-D pin-cell factory :func:`pwr_pin_2d`
  and the equal-volume subdivision of one interval.
* :mod:`~orpheus.mesh.axis`: the per-axis primitives (:class:`AxisMesh`,
  :class:`RadialAxisMesh`, :class:`AxisCoord`) whose tensor product is a
  structured phase-space mesh, and the pure shape functions on axis tuples.

Meshes read from external meshers (unstructured) will live here too.
"""

from __future__ import annotations

from orpheus.mesh.axis import (
    Axis1D,
    AxisCoord,
    AxisMesh,
    FaceLabel,
    RadialAxisMesh,
    axes_from_legacy_mesh,
    coord_system,
    face_labels,
    face_outflow_ordinates,
    face_shape,
    legacy_mesh_from_axes,
    n_unknowns_flat,
    spatial_shape,
)
from orpheus.mesh.factories import pwr_pin_2d
from orpheus.mesh.structured import Mesh1D, Mesh2D, RegionMesh

__all__ = [
    "Axis1D",
    "AxisCoord",
    "AxisMesh",
    "FaceLabel",
    "Mesh1D",
    "Mesh2D",
    "RadialAxisMesh",
    "RegionMesh",
    "axes_from_legacy_mesh",
    "coord_system",
    "face_labels",
    "face_outflow_ordinates",
    "face_shape",
    "legacy_mesh_from_axes",
    "n_unknowns_flat",
    "pwr_pin_2d",
    "spatial_shape",
]
