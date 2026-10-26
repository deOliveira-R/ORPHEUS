"""The bare 2-D mesh: its derived views and its construction refusals.

``Mesh2D(edges_x, edges_y, mat_map, coord, face_roles)`` is the 2-D
discretisation: its cell widths ``dx`` and ``dy``, its cell volumes in the
Cartesian and the cylindrical ``(r, z)`` chart, its flat material ids, its cell
counts, its immutability, and the two refusals of its constructor (a material
map of the wrong shape, a spherical chart). The face roles of a 2-D mesh are
``tests/gates/mesh/test_mesh2d_face_roles.py`` and
``tests/gates/mesh/test_face_roles.py``; its volume measure is
``tests/gates/mesh/test_volume_measure.py``.

``TestMesh2D`` lived in ``tests/gates/geometry/test_geometry.py`` until the
layer-DAG re-filing (`.claude/plans/layer_dag.md`, ruling 21): a claim about an
``orpheus.mesh`` object is stated in ``mesh`` vocabulary. Every test is
``foundation``.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import CoordSystem, MirrorDeck
from orpheus.mesh import Mesh2D

pytestmark = pytest.mark.foundation


# ═══════════════════════════════════════════════════════════════════════
# Mesh2D properties and validation
# ═══════════════════════════════════════════════════════════════════════

class TestMesh2D:
    """Mesh2D derived properties and input validation."""

    def test_dx_dy(self):
        mesh = Mesh2D(
            edges_x=np.array([0.0, 1.0, 3.0]),
            edges_y=np.array([0.0, 0.5]),
            mat_map=np.array([[0], [1]]),
            face_roles={"xmin": MirrorDeck(), "xmax": MirrorDeck(), "ymin": MirrorDeck(), "ymax": MirrorDeck()},
        )
        np.testing.assert_allclose(mesh.dx, [1.0, 2.0])
        np.testing.assert_allclose(mesh.dy, [0.5])

    def test_volumes_cartesian(self):
        mesh = Mesh2D(
            edges_x=np.array([0.0, 2.0]),
            edges_y=np.array([0.0, 3.0]),
            mat_map=np.array([[0]]),
            face_roles={"xmin": MirrorDeck(), "xmax": MirrorDeck(), "ymin": MirrorDeck(), "ymax": MirrorDeck()},
        )
        np.testing.assert_allclose(mesh.volumes, [[6.0]])

    def test_volumes_cylindrical_rz(self):
        mesh = Mesh2D(
            edges_x=np.array([0.0, 1.0]),  # radial
            edges_y=np.array([0.0, 5.0]),   # axial
            mat_map=np.array([[0]]),
            coord=CoordSystem.CYLINDRICAL,
            face_roles={"xmax": 0, "ymin": MirrorDeck(), "ymax": MirrorDeck()},
        )
        np.testing.assert_allclose(mesh.volumes, [[np.pi * 5.0]])

    def test_mat_ids_flat(self):
        mat_map = np.array([[0, 1], [2, 3]])
        mesh = Mesh2D(
            edges_x=np.array([0.0, 1.0, 2.0]),
            edges_y=np.array([0.0, 1.0, 2.0]),
            mat_map=mat_map,
            face_roles={"xmin": MirrorDeck(), "xmax": MirrorDeck(), "ymin": MirrorDeck(), "ymax": MirrorDeck()},
        )
        np.testing.assert_array_equal(mesh.mat_ids, [0, 1, 2, 3])

    def test_nx_ny(self):
        mesh = Mesh2D(
            edges_x=np.linspace(0, 1, 4),
            edges_y=np.linspace(0, 1, 6),
            mat_map=np.zeros((3, 5), dtype=int),
            face_roles={"xmin": MirrorDeck(), "xmax": MirrorDeck(), "ymin": MirrorDeck(), "ymax": MirrorDeck()},
        )
        assert mesh.nx == 3
        assert mesh.ny == 5

    def test_frozen(self):
        mesh = Mesh2D(
            edges_x=np.array([0.0, 1.0]),
            edges_y=np.array([0.0, 1.0]),
            mat_map=np.array([[0]]),
            face_roles={"xmin": MirrorDeck(), "xmax": MirrorDeck(), "ymin": MirrorDeck(), "ymax": MirrorDeck()},
        )
        with pytest.raises(AttributeError):
            mesh.edges_x = np.array([0.0, 2.0])

    def test_wrong_mat_map_shape_raises(self):
        with pytest.raises(ValueError, match="mat_map shape"):
            Mesh2D(
                edges_x=np.array([0.0, 1.0, 2.0]),
                edges_y=np.array([0.0, 1.0]),
                mat_map=np.array([[0, 1]]),  # should be (2, 1)
                face_roles={"xmin": MirrorDeck(), "xmax": MirrorDeck(), "ymin": MirrorDeck(), "ymax": MirrorDeck()},
            )

    def test_spherical_2d_raises(self):
        with pytest.raises(ValueError, match="CARTESIAN or CYLINDRICAL"):
            Mesh2D(
                edges_x=np.array([0.0, 1.0]),
                edges_y=np.array([0.0, 1.0]),
                mat_map=np.array([[0]]),
                coord=CoordSystem.SPHERICAL,
                face_roles={"xmin": MirrorDeck(), "xmax": MirrorDeck(), "ymin": MirrorDeck(), "ymax": MirrorDeck()},
            )
