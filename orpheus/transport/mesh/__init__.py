r"""Mesh plus materials — the method-agnostic materials cluster (L2).

Houses :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh`, the
mesh-plus-materials carrier that every method-specific mesh subclasses (e.g.
``SNProblem(MaterialMesh)``), and the macroscopic cross-section field
(:mod:`~orpheus.transport.mesh.material_xs_field`). The mesh itself and its
per-axis primitives are the discretisation overlay on the geometry, in
:mod:`orpheus.mesh`, one layer down.

Layer (per ``tests/gates/test_layer_imports.py``): L2 ``transport``. These
modules import only ``mesh`` / ``geometry`` / ``numerics`` / ``data`` (and
sibling ``transport`` modules); they do NOT import any L3 method package.
"""

from __future__ import annotations

from orpheus.transport.mesh.material_mesh import (
    InconsistentMaterialsError,
    MaterialMesh,
)
from orpheus.transport.mesh.material_xs_field import MaterialXSField

__all__ = [
    "InconsistentMaterialsError",
    "MaterialMesh",
    "MaterialXSField",
]
