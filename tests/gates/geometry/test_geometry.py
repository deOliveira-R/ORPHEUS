"""The geometry's measure formulas: cell volumes and face areas.

Tests cover the closed-form volume of a cell and area of a face for every
coordinate system: :func:`compute_volumes_1d` and :func:`compute_areas_1d` in
1-D (slab, cylinder, sphere), :func:`compute_volumes_2d` in 2-D (Cartesian and
cylindrical ``(r, z)``), and the refusal of a spherical 2-D chart.

The mesh classes this file held (``TestMesh1D``, ``TestMesh2D``,
``TestFaceRoles``, ``TestPWRPin2D``, ``TestPWRPin2DEdgeCases``) moved to
``tests/gates/mesh/`` in the layer-DAG re-filing (`.claude/plans/layer_dag.md`,
ruling 21): ``test_mesh1d.py``, ``test_mesh2d.py``, ``test_face_roles.py`` and
``test_factories.py``. The equal-volume invariants of a mesh built from a
geometry are exercised by :mod:`tests.gates.mesh.test_mesher` and
:mod:`tests.gates.mesh.test_partition` via the mesher,
``Mesher(geometry).partition(rule).mesh``.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import (
    CoordSystem,
    compute_areas_1d,
    compute_volumes_1d,
    compute_volumes_2d,
)

# Every test in this file is a FOUNDATION test — it verifies a
# software invariant of orpheus.geometry (the closed-form volume and
# area formulas and one refusal) rather than a physics equation
# labelled in docs/theory/*.rst. Foundation tests are orthogonal to
# the L0..L3 physics-verification ladder (Cardinal Rule 4) and never
# carry @pytest.mark.verifies decorators — they have no theory
# label to verify.
#
# See docs/theory/verification/harness.rst ("Foundation tests — software
# invariants outside the L0..L3 ladder") for the taxonomy.
pytestmark = pytest.mark.foundation


# ═══════════════════════════════════════════════════════════════════════
# Volume formulas (1-D)
# ═══════════════════════════════════════════════════════════════════════

class TestVolumes1D:
    """Volume formulas for all 1-D coordinate systems."""

    def test_cartesian_single_cell(self):
        edges = np.array([1.0, 3.0])
        vol = compute_volumes_1d(CoordSystem.CARTESIAN, edges)
        np.testing.assert_allclose(vol, [2.0])

    def test_cartesian_multiple_cells(self):
        edges = np.array([0.0, 0.5, 1.5, 4.0])
        vol = compute_volumes_1d(CoordSystem.CARTESIAN, edges)
        np.testing.assert_allclose(vol, [0.5, 1.0, 2.5])

    def test_cylindrical_single_cell(self):
        edges = np.array([1.0, 2.0])
        vol = compute_volumes_1d(CoordSystem.CYLINDRICAL, edges)
        expected = np.pi * (4.0 - 1.0)  # pi * (r2^2 - r1^2)
        np.testing.assert_allclose(vol, [expected])

    def test_cylindrical_zero_inner_radius(self):
        edges = np.array([0.0, 1.0])
        vol = compute_volumes_1d(CoordSystem.CYLINDRICAL, edges)
        np.testing.assert_allclose(vol, [np.pi])

    def test_cylindrical_multiple_cells(self):
        edges = np.array([0.0, 0.5, 1.0, 2.0])
        vol = compute_volumes_1d(CoordSystem.CYLINDRICAL, edges)
        expected = np.pi * np.diff(edges**2)
        np.testing.assert_allclose(vol, expected)

    def test_spherical_single_cell(self):
        edges = np.array([1.0, 2.0])
        vol = compute_volumes_1d(CoordSystem.SPHERICAL, edges)
        expected = (4.0 / 3.0) * np.pi * (8.0 - 1.0)
        np.testing.assert_allclose(vol, [expected])

    def test_spherical_zero_inner_radius(self):
        edges = np.array([0.0, 1.0])
        vol = compute_volumes_1d(CoordSystem.SPHERICAL, edges)
        np.testing.assert_allclose(vol, [(4.0 / 3.0) * np.pi])

    def test_spherical_multiple_cells(self):
        edges = np.array([0.0, 0.5, 1.0, 2.0])
        vol = compute_volumes_1d(CoordSystem.SPHERICAL, edges)
        expected = (4.0 / 3.0) * np.pi * np.diff(edges**3)
        np.testing.assert_allclose(vol, expected)


# ═══════════════════════════════════════════════════════════════════════
# Surface formulas (1-D)
# ═══════════════════════════════════════════════════════════════════════

class TestSurfaces1D:
    """Surface area formulas at each edge."""

    def test_cartesian(self):
        edges = np.array([0.0, 1.0, 3.0])
        surf = compute_areas_1d(CoordSystem.CARTESIAN, edges)
        np.testing.assert_allclose(surf, [1.0, 1.0, 1.0])

    def test_cylindrical(self):
        edges = np.array([0.0, 0.5, 1.0])
        surf = compute_areas_1d(CoordSystem.CYLINDRICAL, edges)
        expected = 2.0 * np.pi * edges
        np.testing.assert_allclose(surf, expected)

    def test_spherical(self):
        edges = np.array([0.0, 0.5, 1.0])
        surf = compute_areas_1d(CoordSystem.SPHERICAL, edges)
        expected = 4.0 * np.pi * edges**2
        np.testing.assert_allclose(surf, expected)


# ═══════════════════════════════════════════════════════════════════════
# Volume formulas (2-D)
# ═══════════════════════════════════════════════════════════════════════

class TestVolumes2D:
    """Volume formulas for 2-D coordinate systems."""

    def test_cartesian_uniform(self):
        edges_x = np.array([0.0, 1.0, 2.0])
        edges_y = np.array([0.0, 0.5, 1.5])
        vol = compute_volumes_2d(CoordSystem.CARTESIAN, edges_x, edges_y)
        expected = np.array([[0.5, 1.0], [0.5, 1.0]])
        np.testing.assert_allclose(vol, expected)

    def test_cylindrical_rz(self):
        # r-z: V = pi*(r_out^2 - r_in^2) * dz
        edges_r = np.array([0.0, 1.0, 2.0])
        edges_z = np.array([0.0, 3.0])
        vol = compute_volumes_2d(CoordSystem.CYLINDRICAL, edges_r, edges_z)
        expected = np.pi * np.array([[1.0], [3.0]]) * 3.0
        np.testing.assert_allclose(vol, expected)

    def test_spherical_raises(self):
        edges = np.array([0.0, 1.0])
        with pytest.raises(ValueError, match="2-D volumes not defined"):
            compute_volumes_2d(CoordSystem.SPHERICAL, edges, edges)
