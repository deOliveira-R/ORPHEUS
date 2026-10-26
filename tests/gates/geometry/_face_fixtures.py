"""The face tables a geometry gate and its mesh twin share (layer DAG re-filing, 2026-10-26).

Each table is written by hand from a law, never read from the code under test, and is read by a geometry row and by the
mesh rows that assert the mesh agrees with it: geometry owns it, so ``tests/gates/mesh/`` may import it and
``tests/gates/geometry/`` never imports mesh.

* ``FIRST_AXIS``: the first-axis faces of a slab, a solid and a hollow cylinder, written from the topology law
  (``CoordSystem.boundary_points`` along the first axis). Read by
  ``tests/gates/geometry/test_geometry_faces.py::TestTheOneInventoryRule`` (the rule itself) and
  ``tests/gates/mesh/test_mesh2d_face_roles.py`` (both meshes name the same faces).
* ``SLABS``, ``slab_geometry``, ``reference_group``: palindromic and asymmetric slabs with the order of their colouring
  group counted by hand, and the group by a direct call of ``colouring_symmetry`` on the geometry's raw data. Read by
  ``test_geometry_faces.py::TestTheSymmetryIsTheGeometrys`` (the geometry gives the reference group) and
  ``tests/gates/mesh/test_face_roles.py::TestTheSymmetryIsTheGeometrys`` (the mesh gives the same group).
"""
from __future__ import annotations

from typing import Any

from orpheus.geometry import CoordSystem, MirrorDeck, StructuredGeometry
from orpheus.geometry.colouring import colouring_symmetry

_CART = CoordSystem.CARTESIAN
_CYL = CoordSystem.CYLINDRICAL

#: (coord, first-axis edges, the first-axis faces written by hand).
FIRST_AXIS = {
    "slab": (_CART, (0.0, 1.0, 2.0), ("xmin", "xmax")),
    "solid-cylinder": (_CYL, (0.0, 0.5, 1.0), ("xmax",)),
    "hollow-cylinder": (_CYL, (0.25, 0.5, 1.0), ("xmin", "xmax")),
}


def slab_geometry(breakpoints: tuple[float, ...], ids: tuple[int, ...], low: object, high: object) -> StructuredGeometry:
    return StructuredGeometry.slab(breakpoints, ids, left=low, right=high)  # type: ignore[arg-type]


def reference_group(geometry: StructuredGeometry) -> Any:
    """The expected group, by a direct call on the geometry's raw data (no mesh anywhere on this route)."""
    return colouring_symmetry(
        _CART, (geometry.breakpoints[0], geometry.breakpoints[-1]), geometry.breakpoints, geometry.mat_ids,
        geometry.face_roles[0], geometry.face_roles[1],
    )


#: (label, breakpoints, ids, low role, high role, number of group elements). The counts are written by hand:
#: a palindromic colouring with equal end roles has the identity and the mirror, anything else the identity only.
SLABS = [
    ("one material, one tag", (0.0, 1.0, 2.0), (0, 0), 1, 1, 2),
    ("palindrome, equal tags", (0.0, 1.0, 2.0, 3.0), (0, 1, 0), 1, 1, 2),
    ("palindrome, two mirrors", (0.0, 1.0, 2.0, 3.0), (0, 1, 0), MirrorDeck(), MirrorDeck(), 2),
    ("the float trap (0, .1, .2, .3)", (0.0, 0.1, 0.2, 0.3), (0, 1, 0), 1, 1, 2),
    ("asymmetric materials", (0.0, 1.0, 2.0, 3.0), (0, 1, 2), 1, 1, 1),
    ("asymmetric breakpoints", (0.0, 1.0, 2.5, 3.0), (0, 1, 0), 1, 1, 1),
    ("unequal tags", (0.0, 1.0, 2.0, 3.0), (0, 1, 0), 1, 2, 1),
    ("a mirror and a tag", (0.0, 1.0, 2.0, 3.0), (0, 1, 0), MirrorDeck(), 1, 1),
]
