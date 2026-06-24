r"""Mesh + materials — the method-agnostic geometry/materials cluster (L2).

Houses the coordinate-axis primitives
(:mod:`~orpheus.transport.mesh.axis`) and the macroscopic cross-section
field (:mod:`~orpheus.transport.mesh.material_xs_field`), promoted out of
``orpheus.sn`` because both are method-agnostic transport vocabulary: the
axis machinery is pure coordinate geometry (it depends only on
``geometry`` + ``numerics``), and the XS field reads only mesh+materials
data. :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` — the
mesh+materials carrier that every method-specific mesh subclasses (e.g.
``SNMesh(MaterialMesh)``) — lands here in the same campaign.

Layer (per ``tests/test_layer_imports.py``): L2 ``transport``. These
modules import only ``geometry`` / ``numerics`` / ``data`` (and sibling
``transport`` modules); they do NOT import any L3 method package — that
is precisely what made the promotion possible.
"""

from __future__ import annotations

from orpheus.transport.mesh.axis import (
    AXIS_NAMES,
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
from orpheus.transport.mesh.material_mesh import (
    InconsistentMaterialsError,
    MaterialMesh,
)
from orpheus.transport.mesh.material_xs_field import MaterialXSField

__all__ = [
    "AXIS_NAMES",
    "Axis1D",
    "AxisCoord",
    "AxisMesh",
    "FaceLabel",
    "RadialAxisMesh",
    "axes_from_legacy_mesh",
    "coord_system",
    "face_labels",
    "face_outflow_ordinates",
    "face_shape",
    "legacy_mesh_from_axes",
    "n_unknowns_flat",
    "spatial_shape",
    "InconsistentMaterialsError",
    "MaterialMesh",
    "MaterialXSField",
]
