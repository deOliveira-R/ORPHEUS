"""The meshes of the dropped-law refusals of CP, MoC and MC (#513, #514; layer-DAG re-filing, Phase 1b).

The three method files ``tests/gates/cp/test_cp_dropped_laws_are_refused.py``,
``tests/gates/moc/test_moc_dropped_laws_are_refused.py`` and ``tests/gates/mc/test_mc_dropped_laws_are_refused.py`` were one file,
``tests/gates/mesh/test_dropped_laws_are_refused.py``, until the re-filing. Each method refuses a declaration it would drop, and the
three may not import one another, so the two mesh builders they share live here, one layer below all three. The bodies are the
builders of that file, unchanged.
"""
from __future__ import annotations

from orpheus.geometry import CoordSystem, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesh1D, Mesher
from tests.gates._exterior import declare_sides


def _slab(left, right) -> tuple[Mesh1D, dict]:
    """A slab with its faces declared BY VALUE (rule R1: equal conditions share a tag): the mesh and its ``{tag: condition}``."""
    roles, conditions = declare_sides(left, right)
    geometry = StructuredGeometry.slab((0.0, 2.0), (0,), left=roles[0], right=roles[1])
    return Mesher(geometry).partition(CellsByCount.uniform_width(4)).mesh, conditions


def _hollow(coord: CoordSystem, inner: int = 1, outer: int = 0) -> Mesh1D:
    """A hollow body [0.5, 2.0]: the inner face is the tag ``inner`` and the outer face the tag ``outer`` (rule R6)."""
    geometry = StructuredGeometry(coord=coord, breakpoints=(0.5, 2.0), mat_ids=(0,), face_roles=(inner, outer))
    return Mesher(geometry).partition(CellsByCount.uniform_volume(4)).mesh
