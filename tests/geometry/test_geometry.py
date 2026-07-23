"""Comprehensive tests for the geometry module.

Tests cover:
- Volume and surface formulas for all coordinate systems (1-D and 2-D)
- Mesh validation (monotonicity, shape, frozen immutability)
- :class:`BC` declaration dataclass and mesh integration
- :func:`pwr_pin_2d` (2-D Cartesian factory — kept as a non-trivial
  helper after Phase F retired the 1-D Zone/factories surface)

The 1-D ``Zone`` / ``mesh1d_from_zones`` / ``pwr_*`` / ``homogeneous_1d``
/ ``slab_fuel_moderator`` factories were retired in Phase F. Their
equal-volume invariants and discretization correctness are now
exercised by :mod:`tests.geometry.test_structured_geometry` via the
:class:`StructuredGeometry` → :meth:`Mesh1D.from_geometry` flow.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import (
    BC,
    CoordSystem,
    Mesh1D,
    Mesh2D,
    compute_areas_1d,
    compute_volumes_1d,
    compute_volumes_2d,
    pwr_pin_2d,
)

# Every test in this file is a FOUNDATION test — it verifies a
# software invariant of orpheus.geometry (volume formulas, factory
# outputs, frozen immutability, input validation, algebraic
# subdivision identities) rather than a physics equation labelled
# in docs/theory/*.rst. Foundation tests are orthogonal to the
# L0..L3 physics-verification ladder (Cardinal Rule 4) and never
# carry @pytest.mark.verifies decorators — they have no theory
# label to verify.
#
# The ERR-020 bit-exact equal-volume tests in TestZoneSubdivision
# are the canonical foundation test: they assert that every cell
# in an equal-volume zone has bit-identical volume by construction,
# which is an algebraic invariant of Mesh1D.volumes, not a physics
# claim. The catches("ERR-020") decorators on those specific tests
# stay — ERR catalogue coverage is orthogonal to the marker type.
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


# ═══════════════════════════════════════════════════════════════════════
# pwr_pin_2d — 2-D Cartesian PWR pin factory (post-Phase-F survivor)
# ═══════════════════════════════════════════════════════════════════════

class TestPWRPin2D:
    """:func:`pwr_pin_2d` 2-D factory invariants."""

    def test_pin_2d_shape(self):
        mesh = pwr_pin_2d(n_cells=10)
        assert mesh.nx == 10
        assert mesh.ny == 10
        assert mesh.mat_map.shape == (10, 10)

    def test_pin_2d_has_all_materials(self):
        mesh = pwr_pin_2d(n_cells=20)
        mats = set(mesh.mat_map.ravel())
        assert mats == {0, 1, 2}

    def test_pin_2d_mat_ids_flat(self):
        """mat_ids returns flat array for assemble_cell_xs."""
        mesh = pwr_pin_2d(n_cells=5)
        assert mesh.mat_ids.shape == (25,)


# ═══════════════════════════════════════════════════════════════════════
# Mesh1D properties and validation
# ═══════════════════════════════════════════════════════════════════════

class TestMesh1D:
    """Mesh1D derived properties and input validation."""

    def test_widths(self):
        mesh = Mesh1D(edges=np.array([0.0, 1.0, 3.0, 6.0]),
                      mat_ids=np.array([0, 1, 2]))
        np.testing.assert_allclose(mesh.widths, [1.0, 2.0, 3.0])

    def test_centers(self):
        mesh = Mesh1D(edges=np.array([0.0, 2.0, 6.0]),
                      mat_ids=np.array([0, 1]))
        np.testing.assert_allclose(mesh.centers, [1.0, 4.0])

    def test_total_width(self):
        mesh = Mesh1D(edges=np.array([1.0, 3.0, 7.0]),
                      mat_ids=np.array([0, 1]))
        assert mesh.total_width == 6.0

    def test_N(self):
        mesh = Mesh1D(edges=np.arange(6.0),
                      mat_ids=np.zeros(5, dtype=int))
        assert mesh.N == 5

    def test_frozen(self):
        mesh = Mesh1D(edges=np.array([0.0, 1.0]),
                      mat_ids=np.array([0]))
        with pytest.raises(AttributeError):
            mesh.edges = np.array([0.0, 2.0])

    def test_non_monotonic_edges_raises(self):
        with pytest.raises(ValueError, match="monotonically increasing"):
            Mesh1D(edges=np.array([0.0, 2.0, 1.0]),
                   mat_ids=np.array([0, 1]))

    def test_equal_edges_raises(self):
        with pytest.raises(ValueError, match="monotonically increasing"):
            Mesh1D(edges=np.array([0.0, 1.0, 1.0]),
                   mat_ids=np.array([0, 1]))

    def test_wrong_mat_ids_length_raises(self):
        with pytest.raises(ValueError, match="len\\(mat_ids\\)"):
            Mesh1D(edges=np.array([0.0, 1.0, 2.0]),
                   mat_ids=np.array([0, 1, 2]))

    def test_too_few_edges_raises(self):
        with pytest.raises(ValueError, match="at least 2"):
            Mesh1D(edges=np.array([0.0]), mat_ids=np.array([]))

    def test_coerces_to_float_and_int(self):
        """Accepts lists; coerces edges to float, mat_ids to int."""
        mesh = Mesh1D(edges=[0, 1, 2], mat_ids=[0, 1])
        assert mesh.edges.dtype == float
        assert mesh.mat_ids.dtype == int


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
        )
        np.testing.assert_allclose(mesh.dx, [1.0, 2.0])
        np.testing.assert_allclose(mesh.dy, [0.5])

    def test_volumes_cartesian(self):
        mesh = Mesh2D(
            edges_x=np.array([0.0, 2.0]),
            edges_y=np.array([0.0, 3.0]),
            mat_map=np.array([[0]]),
        )
        np.testing.assert_allclose(mesh.volumes, [[6.0]])

    def test_volumes_cylindrical_rz(self):
        mesh = Mesh2D(
            edges_x=np.array([0.0, 1.0]),  # radial
            edges_y=np.array([0.0, 5.0]),   # axial
            mat_map=np.array([[0]]),
            coord=CoordSystem.CYLINDRICAL,
        )
        np.testing.assert_allclose(mesh.volumes, [[np.pi * 5.0]])

    def test_mat_ids_flat(self):
        mat_map = np.array([[0, 1], [2, 3]])
        mesh = Mesh2D(
            edges_x=np.array([0.0, 1.0, 2.0]),
            edges_y=np.array([0.0, 1.0, 2.0]),
            mat_map=mat_map,
        )
        np.testing.assert_array_equal(mesh.mat_ids, [0, 1, 2, 3])

    def test_nx_ny(self):
        mesh = Mesh2D(
            edges_x=np.linspace(0, 1, 4),
            edges_y=np.linspace(0, 1, 6),
            mat_map=np.zeros((3, 5), dtype=int),
        )
        assert mesh.nx == 3
        assert mesh.ny == 5

    def test_frozen(self):
        mesh = Mesh2D(
            edges_x=np.array([0.0, 1.0]),
            edges_y=np.array([0.0, 1.0]),
            mat_map=np.array([[0]]),
        )
        with pytest.raises(AttributeError):
            mesh.edges_x = np.array([0.0, 2.0])

    def test_wrong_mat_map_shape_raises(self):
        with pytest.raises(ValueError, match="mat_map shape"):
            Mesh2D(
                edges_x=np.array([0.0, 1.0, 2.0]),
                edges_y=np.array([0.0, 1.0]),
                mat_map=np.array([[0, 1]]),  # should be (2, 1)
            )

    def test_spherical_2d_raises(self):
        with pytest.raises(ValueError, match="CARTESIAN or CYLINDRICAL"):
            Mesh2D(
                edges_x=np.array([0.0, 1.0]),
                edges_y=np.array([0.0, 1.0]),
                mat_map=np.array([[0]]),
                coord=CoordSystem.SPHERICAL,
            )


# ═══════════════════════════════════════════════════════════════════════
# BC dataclass and mesh integration
# ═══════════════════════════════════════════════════════════════════════

class TestBC:
    """BC declaration dataclass and integration with mesh fields."""

    def test_bc_creation_kind_only(self):
        bc = BC("vacuum")
        assert bc.kind == "vacuum"
        assert bc.params == {}

    def test_bc_creation_with_params(self):
        bc = BC("white", {"albedo": 0.7})
        assert bc.kind == "white"
        assert bc.params == {"albedo": 0.7}

    def test_bc_frozen(self):
        bc = BC("vacuum")
        with pytest.raises(AttributeError):
            bc.kind = "reflective"

    def test_bc_equality(self):
        assert BC("vacuum") == BC("vacuum")
        assert BC("vacuum") != BC("reflective")
        assert BC("white", {"albedo": 0.7}) == BC("white", {"albedo": 0.7})
        assert BC("white", {"albedo": 0.7}) != BC("white", {"albedo": 0.5})

    def test_bc_convenience_instances(self):
        assert BC.vacuum == BC("vacuum")
        assert BC.reflective == BC("reflective")
        assert BC.white == BC("white")

    def test_bc_repr(self):
        assert repr(BC("vacuum")) == "BC('vacuum')"
        assert repr(BC("white", {"albedo": 0.7})) == "BC('white', {'albedo': 0.7})"

    # ── Mesh1D BC fields ─────────────────────────────────────────────

    def test_mesh1d_bc_defaults_none(self):
        mesh = Mesh1D(edges=[0, 1, 2], mat_ids=[0, 1])
        assert mesh.bc_left is None
        assert mesh.bc_right is None

    def test_mesh1d_bc_explicit(self):
        mesh = Mesh1D(
            edges=[0, 1, 2], mat_ids=[0, 1],
            bc_left=BC.reflective, bc_right=BC.vacuum,
        )
        assert mesh.bc_left == BC("reflective")
        assert mesh.bc_right == BC("vacuum")

    def test_mesh1d_bc_frozen(self):
        mesh = Mesh1D(
            edges=[0, 1, 2], mat_ids=[0, 1],
            bc_left=BC.reflective,
        )
        with pytest.raises(AttributeError):
            mesh.bc_left = BC.vacuum

    def test_mesh1d_bc_invalid_type_raises(self):
        with pytest.raises(TypeError, match="bc_left must be a BC instance"):
            Mesh1D(edges=[0, 1], mat_ids=[0], bc_left="vacuum")

    # ── BC.to_alpha — production-tag → continuous-albedo bridge ─────

    def test_bc_to_alpha_vacuum(self):
        assert BC.vacuum.to_alpha() == 0.0

    def test_bc_to_alpha_reflective(self):
        assert BC.reflective.to_alpha() == 1.0

    def test_bc_to_alpha_partial_albedo(self):
        assert BC("partial", {"albedo": 0.7}).to_alpha() == 0.7
        assert BC("partial", {"albedo": 0.0}).to_alpha() == 0.0
        assert BC("partial", {"albedo": 1.0}).to_alpha() == 1.0

    def test_bc_to_alpha_partial_missing_albedo_raises(self):
        with pytest.raises(ValueError, match="missing the 'albedo' parameter"):
            BC("partial").to_alpha()

    def test_bc_to_alpha_unsupported_kind_raises(self):
        with pytest.raises(NotImplementedError, match="no specular-albedo equivalent"):
            BC.white.to_alpha()

    def test_mesh1d_backward_compat(self):
        """All existing Mesh1D constructors (no BC args) still work."""
        m1 = Mesh1D(edges=np.array([0.0, 1.0, 3.0, 6.0]),
                     mat_ids=np.array([0, 1, 2]))
        assert m1.N == 3
        m2 = Mesh1D(edges=[0, 1], mat_ids=[0],
                     coord=CoordSystem.CYLINDRICAL)
        assert m2.coord == CoordSystem.CYLINDRICAL

    # ── Mesh2D BC fields ─────────────────────────────────────────────

    def test_mesh2d_bc_defaults_none(self):
        mesh = Mesh2D(
            edges_x=[0, 1], edges_y=[0, 1], mat_map=np.array([[0]]),
        )
        assert mesh.bc_xmin is None
        assert mesh.bc_xmax is None
        assert mesh.bc_ymin is None
        assert mesh.bc_ymax is None

    def test_mesh2d_bc_explicit(self):
        mesh = Mesh2D(
            edges_x=[0, 1], edges_y=[0, 1], mat_map=np.array([[0]]),
            bc_xmin=BC.reflective, bc_xmax=BC.vacuum,
            bc_ymin=BC.reflective, bc_ymax=BC.vacuum,
        )
        assert mesh.bc_xmin == BC("reflective")
        assert mesh.bc_xmax == BC("vacuum")

    def test_mesh2d_bc_invalid_type_raises(self):
        with pytest.raises(TypeError, match="bc_xmin must be a BC instance"):
            Mesh2D(
                edges_x=[0, 1], edges_y=[0, 1], mat_map=np.array([[0]]),
                bc_xmin="reflective",
            )


# ═══════════════════════════════════════════════════════════════════════
# pwr_pin_2d edge cases
# ═══════════════════════════════════════════════════════════════════════

class TestPWRPin2DEdgeCases:
    """Edge cases for :func:`pwr_pin_2d`."""

    def test_pin_2d_wrong_mat_ids_length_raises(self):
        with pytest.raises(ValueError, match="len\\(mat_ids\\)"):
            pwr_pin_2d(radii=[1.0], mat_ids=[0, 1, 2])
