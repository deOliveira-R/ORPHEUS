"""Tests for the SN boundary condition infrastructure.

Verifies the BC_REGISTRY pattern: declaration on geometry, resolution at
SNMesh construction, and correct behavior in sweeps.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Mesh2D, CoordSystem
from orpheus.sn.geometry import SNMesh
from orpheus.sn.quadrature import GaussLegendre1D


@pytest.fixture
def quad():
    return GaussLegendre1D.create(4)


@pytest.fixture
def slab_mesh():
    return Mesh1D(edges=np.linspace(0, 5, 11), mat_ids=np.zeros(10, dtype=int))


# ═══════════════════════════════════════════════════════════════════════
# Registry tests
# ═══════════════════════════════════════════════════════════════════════

class TestSNBCRegistry:
    """BC_REGISTRY is the single source of truth: resolves, validates, advertises."""

    def test_registry_keys(self):
        assert "vacuum" in SNMesh.BC_REGISTRY
        assert "reflective" in SNMesh.BC_REGISTRY

    def test_registry_docstrings(self):
        """Every factory has a docstring (used as description for UI query)."""
        for kind, factory in SNMesh.BC_REGISTRY.items():
            assert factory.__doc__ is not None, f"BC factory '{kind}' has no docstring"

    def test_registry_programmatic_query(self):
        """Descriptions are queryable via factory docstrings."""
        descriptions = {k: v.__doc__ for k, v in SNMesh.BC_REGISTRY.items()}
        assert "vacuum" in descriptions
        assert "reflective" in descriptions


# ═══════════════════════════════════════════════════════════════════════
# Resolution tests
# ═══════════════════════════════════════════════════════════════════════

class TestSNBCResolution:
    """BC resolution at SNMesh construction time."""

    def test_default_is_reflective(self, slab_mesh, quad):
        """None on mesh resolves to 'reflective' (eigenvalue default)."""
        sn = SNMesh(slab_mesh, quad)
        assert sn.bc_left == "reflective"
        assert sn.bc_right == "reflective"

    def test_explicit_vacuum(self, slab_mesh, quad):
        mesh = Mesh1D(
            edges=slab_mesh.edges, mat_ids=slab_mesh.mat_ids,
            bc_left=BC.vacuum, bc_right=BC.vacuum,
        )
        sn = SNMesh(mesh, quad)
        assert sn.bc_left == "vacuum"
        assert sn.bc_right == "vacuum"

    def test_mixed_bcs(self, slab_mesh, quad):
        mesh = Mesh1D(
            edges=slab_mesh.edges, mat_ids=slab_mesh.mat_ids,
            bc_left=BC.reflective, bc_right=BC.vacuum,
        )
        sn = SNMesh(mesh, quad)
        assert sn.bc_left == "reflective"
        assert sn.bc_right == "vacuum"

    def test_unknown_bc_raises(self, slab_mesh, quad):
        mesh = Mesh1D(
            edges=slab_mesh.edges, mat_ids=slab_mesh.mat_ids,
            bc_left=BC("white"),
        )
        with pytest.raises(ValueError, match="does not support.*'white'"):
            SNMesh(mesh, quad)

    def test_error_lists_supported(self, slab_mesh, quad):
        mesh = Mesh1D(
            edges=slab_mesh.edges, mat_ids=slab_mesh.mat_ids,
            bc_left=BC("periodic"),
        )
        with pytest.raises(ValueError, match="'reflective'.*'vacuum'"):
            SNMesh(mesh, quad)

    def test_2d_mesh_resolution(self, quad):
        mesh = Mesh2D(
            edges_x=np.linspace(0, 2, 3),
            edges_y=np.linspace(0, 2, 3),
            mat_map=np.zeros((2, 2), dtype=int),
            bc_xmin=BC.reflective, bc_xmax=BC.vacuum,
            bc_ymin=BC.reflective, bc_ymax=BC.vacuum,
        )
        sn = SNMesh(mesh, quad)
        assert sn.bc_xmin == "reflective"
        assert sn.bc_xmax == "vacuum"
        assert sn.bc_ymin == "reflective"
        assert sn.bc_ymax == "vacuum"

    def test_curvilinear_vacuum_raises(self, quad):
        """Spherical/cylindrical only support reflective; vacuum raises at resolution."""
        mesh = Mesh1D(
            edges=np.linspace(0.1, 1.0, 6), mat_ids=np.zeros(5, dtype=int),
            coord=CoordSystem.SPHERICAL,
            bc_right=BC.vacuum,
        )
        with pytest.raises(NotImplementedError, match="spherical.*reflective"):
            SNMesh(mesh, quad)


# ═══════════════════════════════════════════════════════════════════════
# Sweep behavior tests
# ═══════════════════════════════════════════════════════════════════════

class TestSNBCSweepBehavior:
    """Verify that resolved BCs produce correct sweep behavior."""

    def test_vacuum_keff_lower_than_reflective(self, quad):
        """Vacuum BC loses neutrons → lower keff than reflective."""
        from orpheus.derivations.reference_values import get
        from orpheus.geometry import homogeneous_1d
        from orpheus.sn.solver import solve_sn

        case = get("sn_slab_2eg_1rg")
        mix = next(iter(case.materials.values()))
        materials = {0: mix}

        # Reflective (default — no BCs on mesh)
        mesh_refl = homogeneous_1d(20, 2.0, mat_id=0)
        result_refl = solve_sn(materials, mesh_refl, quad)

        # Vacuum — explicit BCs on mesh
        mesh_vac = Mesh1D(
            edges=mesh_refl.edges, mat_ids=mesh_refl.mat_ids,
            bc_left=BC.vacuum, bc_right=BC.vacuum,
        )
        result_vac = solve_sn(materials, mesh_vac, quad)

        # Reflective has higher keff (no leakage vs leakage)
        assert result_refl.keff > result_vac.keff
