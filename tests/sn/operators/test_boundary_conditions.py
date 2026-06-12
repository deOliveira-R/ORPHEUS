"""Tests for the SN boundary condition infrastructure.

Verifies the BOUNDARY_OPERATOR_REGISTRY pattern: declaration on geometry, resolution at
SNMesh construction, and correct behavior in sweeps.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Mesh2D, CoordSystem
from orpheus.sn.geometry import SNMesh
from orpheus.numerics.quadrature import Quadrature
from tests.sn._test_helpers import placeholder_materials

# SN boundary-condition infrastructure: structural invariants of the
# SNMesh BC wiring (no theory-page :label:). Foundation, not a physics
# equation gate. (Was a V&V orphan before the taxonomy reorg forced a
# marker — see .claude/plans/sn_test_taxonomy.md.)
pytestmark = pytest.mark.foundation


@pytest.fixture
def quad():
    return Quadrature.gauss_legendre(4)


@pytest.fixture
def slab_mesh():
    return Mesh1D(edges=np.linspace(0, 5, 11), mat_ids=np.zeros(10, dtype=int))


# ═══════════════════════════════════════════════════════════════════════
# Registry tests
# ═══════════════════════════════════════════════════════════════════════

class TestSNBCRegistry:
    """BOUNDARY_OPERATOR_REGISTRY is the single source of truth: resolves, validates, advertises."""

    def test_registry_keys(self):
        assert "vacuum" in SNMesh.BOUNDARY_OPERATOR_REGISTRY
        assert "reflective" in SNMesh.BOUNDARY_OPERATOR_REGISTRY

    def test_registry_docstrings(self):
        """Every factory has a docstring (used as description for UI query)."""
        for kind, factory in SNMesh.BOUNDARY_OPERATOR_REGISTRY.items():
            assert factory.__doc__ is not None, f"BC factory '{kind}' has no docstring"

    def test_registry_programmatic_query(self):
        """Descriptions are queryable via factory docstrings."""
        descriptions = {k: v.__doc__ for k, v in SNMesh.BOUNDARY_OPERATOR_REGISTRY.items()}
        assert "vacuum" in descriptions
        assert "reflective" in descriptions


# ═══════════════════════════════════════════════════════════════════════
# Resolution tests
# ═══════════════════════════════════════════════════════════════════════

class TestSNBCResolution:
    """BC resolution at SNMesh construction time."""

    def test_default_is_reflective(self, slab_mesh, quad):
        """None on mesh resolves to 'reflective' (eigenvalue default)."""
        sn = SNMesh(slab_mesh, quad, placeholder_materials())
        assert sn.bc["xmin"] == "reflective"
        assert sn.bc["xmax"] == "reflective"

    def test_explicit_vacuum(self, slab_mesh, quad):
        mesh = Mesh1D(
            edges=slab_mesh.edges, mat_ids=slab_mesh.mat_ids,
            bc_left=BC.vacuum, bc_right=BC.vacuum,
        )
        sn = SNMesh(mesh, quad, placeholder_materials())
        assert sn.bc["xmin"] == "vacuum"
        assert sn.bc["xmax"] == "vacuum"

    def test_mixed_bcs(self, slab_mesh, quad):
        mesh = Mesh1D(
            edges=slab_mesh.edges, mat_ids=slab_mesh.mat_ids,
            bc_left=BC.reflective, bc_right=BC.vacuum,
        )
        sn = SNMesh(mesh, quad, placeholder_materials())
        assert sn.bc["xmin"] == "reflective"
        assert sn.bc["xmax"] == "vacuum"

    def test_unknown_bc_raises(self, slab_mesh, quad):
        mesh = Mesh1D(
            edges=slab_mesh.edges, mat_ids=slab_mesh.mat_ids,
            bc_left=BC("white"),
        )
        with pytest.raises(ValueError, match="does not support.*'white'"):
            SNMesh(mesh, quad, placeholder_materials())

    def test_error_lists_supported(self, slab_mesh, quad):
        mesh = Mesh1D(
            edges=slab_mesh.edges, mat_ids=slab_mesh.mat_ids,
            bc_left=BC("periodic"),
        )
        with pytest.raises(ValueError, match="'reflective'.*'vacuum'"):
            SNMesh(mesh, quad, placeholder_materials())

    def test_2d_mesh_resolution(self, quad):
        mesh = Mesh2D(
            edges_x=np.linspace(0, 2, 3),
            edges_y=np.linspace(0, 2, 3),
            mat_map=np.zeros((2, 2), dtype=int),
            bc_xmin=BC.reflective, bc_xmax=BC.vacuum,
            bc_ymin=BC.reflective, bc_ymax=BC.vacuum,
        )
        sn = SNMesh(mesh, quad, placeholder_materials())
        assert sn.bc["xmin"] == "reflective"
        assert sn.bc["xmax"] == "vacuum"
        assert sn.bc["ymin"] == "reflective"
        assert sn.bc["ymax"] == "vacuum"

    def test_curvilinear_vacuum_resolves(self, quad):
        """Spherical/cylindrical accept vacuum; sweep enforces zero-incoming
        at the outer face. Vacuum support added in commits 655e3e5 / 37c5bbf
        (the curvilinear-only-reflective gate was removed and the inward
        sweep now branches on ``is_vacuum_outer``)."""
        mesh = Mesh1D(
            edges=np.linspace(0.1, 1.0, 6), mat_ids=np.zeros(5, dtype=int),
            coord=CoordSystem.SPHERICAL,
            bc_right=BC.vacuum,
        )
        sn = SNMesh(mesh, quad, placeholder_materials())
        assert sn.bc["xmax"] == "vacuum"


# ═══════════════════════════════════════════════════════════════════════
# Sweep behavior tests
# ═══════════════════════════════════════════════════════════════════════

class TestSNBCSweepBehavior:
    """Verify that resolved BCs produce correct sweep behavior."""

    @pytest.mark.catches("ERR-052")
    def test_vacuum_keff_lower_than_reflective(self, quad):
        """Vacuum BC loses neutrons → lower keff than reflective.

        Catches ERR-052: power-iteration flux non-renormalisation. Without
        the per-step ``ψ /= ||ψ||`` projection in ``power_iteration`` /
        ``KEigenvalue.solve``, subcritical cases (vacuum 2cm slab here)
        underflow to denormalised FP within ~30-60 outer iterations, and
        the keff ratio becomes numerically meaningless.
        """
        from orpheus.derivations.reference_values import get
        from orpheus.geometry import (
            Mesh1D as _Mesh1D,
            Region,
            RegionMesh,
            StructuredGeometry,
        )
        from orpheus.sn.solver import solve_sn

        case = get("sn_slab_2eg_1rg")
        mix = next(iter(case.materials.values()))
        materials = {0: mix}

        # Reflective (default — no BCs on mesh)
        mesh_refl = _Mesh1D.from_geometry(
            StructuredGeometry(
                geometry="SLB",
                regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
                bcs=(BC.reflective, BC.reflective),
            ),
            region_meshes=(RegionMesh(n_cells=20),),
        )
        result_refl = solve_sn(materials, mesh_refl, quad)

        # Vacuum — explicit BCs on mesh
        mesh_vac = Mesh1D(
            edges=mesh_refl.edges, mat_ids=mesh_refl.mat_ids,
            bc_left=BC.vacuum, bc_right=BC.vacuum,
        )
        result_vac = solve_sn(materials, mesh_vac, quad)

        # Reflective has higher keff (no leakage vs leakage)
        assert result_refl.keff > result_vac.keff
