r"""A geometry's faces: the inventory rule, the group its end roles select, the roles of the named cells.

Three claims about a :class:`~orpheus.geometry.StructuredGeometry` and its faces, each stated in geometry's vocabulary
alone (the layer-DAG re-filing, ``.claude/plans/layer_dag.md``, 2026-10-26, candidate 2 of the Phase 3 sub-cone report:
they lived in ``tests/gates/mesh/`` beside the mesh rows that compare a mesh with them).

* ``TestTheOneInventoryRule``: :func:`~orpheus.geometry.faces.face_inventory` names the faces by the topology law in
  one, two and three dimensions; the expected faces are written by hand
  (``tests/gates/geometry/_face_fixtures.py::FIRST_AXIS``). Its two-mesh rows (both meshes name the same first-axis
  faces, both constructors ask the one rule) stay in ``tests/gates/mesh/test_mesh2d_face_roles.py::TestTheOneInventoryRule``.
* ``TestTheSymmetryIsTheGeometrys``: ``StructuredGeometry.symmetry`` is the group a direct call of ``colouring_symmetry``
  gives on the geometry's raw breakpoints, ids and end roles, with the order counted by hand (gate 6 of wave 2, law L-C
  input). The rows that a mesh gives the same group stay in ``tests/gates/mesh/test_face_roles.py::TestTheSymmetryIsTheGeometrys``
  and rest on this one.
* ``TestTheNamedCells``: ``wigner_seitz_pin_cell`` is the cylindricalised lattice cell, tag 0 on its outer surface (the
  model binds it to a white wall elsewhere); ``pwr_slab_half_cell`` is bounded by two symmetry planes, a mirror deck on
  both (the user's ruling of 2026-09-30: the role is the model's).

Every test is ``foundation``; claim kind THEOREM. Class names are kept from the source files, which keep their mesh
rows under the same class names: a citation names the file.
"""
from __future__ import annotations

import inspect

import numpy as np
import pytest

from orpheus.geometry import MirrorDeck, StructuredGeometry
from orpheus.geometry.faces import face_inventory
from tests.gates._content_identity_helpers import require
from tests.gates.geometry._face_fixtures import FIRST_AXIS as _FIRST_AXIS
from tests.gates.geometry._face_fixtures import SLABS as _SLABS
from tests.gates.geometry._face_fixtures import reference_group as _reference
from tests.gates.geometry._face_fixtures import slab_geometry as _geometry

pytestmark = pytest.mark.foundation


class TestTheOneInventoryRule:
    """The rule itself: ``face_inventory(coord, first_axis_edges, dimension)`` is ``boundary_points`` along the first
    axis and both ends along every further axis, each named by ``FaceLabel(axis_index, endpoint).face_name``. First red:
    a solid radial first axis given ``xmin`` (the pre-``FaceRoles`` ``_faces_2d`` shape)."""

    @pytest.mark.parametrize("case", list(_FIRST_AXIS))
    def test_the_rule_itself(self, case):
        coord, edges, faces = _FIRST_AXIS[case]
        assert face_inventory(coord, np.asarray(edges), 1) == faces
        assert face_inventory(coord, np.asarray(edges), 2) == faces + ("ymin", "ymax")
        assert face_inventory(coord, np.asarray(edges), 3) == faces + ("ymin", "ymax", "zmin", "zmax")


@pytest.mark.rests_on("tests/gates/geometry/test_colouring_symmetry.py::test_colouring_symmetry_reads_runs_and_roles")
class TestTheSymmetryIsTheGeometrys:
    """Gate 6, the geometry half: the geometry's group is the direct call's, of the hand-counted order, identity first."""

    @pytest.mark.parametrize("label", [row[0] for row in _SLABS])
    def test_the_geometry_gives_the_reference_group(self, label: str) -> None:
        (_, breakpoints, ids, low, high, count), = [row for row in _SLABS if row[0] == label]
        geometry = _geometry(breakpoints, ids, low, high)
        require(len(geometry.symmetry.elements) == count, f"{label}: {len(geometry.symmetry.elements)} elements, expected {count}")
        require(geometry.symmetry == _reference(geometry), f"{label}: the geometry's group is the direct call's")
        require(geometry.symmetry.elements[0].is_identity, "the identity is first")


class TestTheNamedCells:
    """``wigner_seitz_pin_cell`` is the cylindricalised lattice cell, tag 0 on
    its outer surface (the model binds it to a white wall elsewhere);
    ``pwr_slab_half_cell`` is bounded by two symmetry planes, a mirror deck on
    both. Another role is another body (``StructuredGeometry.cylinder`` /
    ``.slab``). First reds: a ``face_roles=`` or ``boundaries=`` keyword
    re-added (the signature row); the model's role changed (the role row)."""

    @pytest.mark.parametrize("name", ["wigner_seitz_pin_cell", "pwr_slab_half_cell"])
    def test_no_role_keyword(self, name):
        parameters = inspect.signature(getattr(StructuredGeometry, name)).parameters
        assert not {"boundaries", "face_roles"} & set(parameters), name
        assert inspect.Parameter.VAR_KEYWORD not in {p.kind for p in parameters.values()}, name

    def test_the_models_roles(self):
        assert StructuredGeometry.wigner_seitz_pin_cell().face_roles == (0,)
        assert StructuredGeometry.pwr_slab_half_cell().face_roles == (MirrorDeck(), MirrorDeck())
