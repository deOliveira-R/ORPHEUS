"""Foundation tests for :class:`StructuredGeometry`, :class:`Region`,
and :meth:`Mesh1D.from_geometry`.

These tests pin the geometry → mesh transition: validation rules at
the geometry layer, validation rules at the mesh layer (RegionMesh),
and the construction semantics that connect them. They are
foundation-tier (software invariants) — they verify the code-shape
contract, not a physics claim.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import (
    BC,
    CoordSystem,
    Mesh1D,
    Region,
    RegionMesh,
    StructuredGeometry,
)


pytestmark = pytest.mark.foundation


# ─────────────────────────────────────────────────────────────────────
# Region — geometry-layer per-region descriptor
# ─────────────────────────────────────────────────────────────────────


class TestRegion:
    def test_minimal_construction(self):
        r = Region(mat_id=0, outer_thickness_cm=1.5)
        assert r.mat_id == 0
        assert r.outer_thickness_cm == 1.5

    def test_frozen(self):
        r = Region(mat_id=0, outer_thickness_cm=1.5)
        with pytest.raises(AttributeError):
            r.mat_id = 1

    def test_zero_thickness_rejected(self):
        with pytest.raises(ValueError, match="must be > 0"):
            Region(mat_id=0, outer_thickness_cm=0.0)

    def test_negative_thickness_rejected(self):
        with pytest.raises(ValueError, match="must be > 0"):
            Region(mat_id=0, outer_thickness_cm=-1.0)

    def test_mat_id_must_be_int(self):
        with pytest.raises(TypeError, match="must be int"):
            Region(mat_id=0.0, outer_thickness_cm=1.0)


# ─────────────────────────────────────────────────────────────────────
# StructuredGeometry — geometry-layer
# ─────────────────────────────────────────────────────────────────────


class TestStructuredGeometryValidation:
    def test_single_region_sphere(self):
        g = StructuredGeometry(
            geometry="SPH",
            regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
            bcs=(BC.vacuum,),
        )
        assert g.geometry == "SPH"
        assert g.coord == CoordSystem.SPHERICAL
        assert g.n_endpoints == 1
        assert g.domain_extent_cm == 2.0

    def test_single_region_cylinder(self):
        g = StructuredGeometry(
            geometry="CYL",
            regions=(Region(mat_id=0, outer_thickness_cm=1.5),),
            bcs=(BC("white"),),
        )
        assert g.coord == CoordSystem.CYLINDRICAL
        assert g.n_endpoints == 1

    def test_single_region_slab_requires_two_bcs(self):
        g = StructuredGeometry(
            geometry="SLB",
            regions=(Region(mat_id=0, outer_thickness_cm=4.0),),
            bcs=(BC.vacuum, BC.vacuum),
        )
        assert g.coord == CoordSystem.CARTESIAN
        assert g.n_endpoints == 2
        assert g.domain_extent_cm == 4.0

    def test_multi_region_slab(self):
        g = StructuredGeometry(
            geometry="SLB",
            regions=(
                Region(mat_id=1, outer_thickness_cm=0.5),
                Region(mat_id=0, outer_thickness_cm=2.0),
                Region(mat_id=1, outer_thickness_cm=0.5),
            ),
            bcs=(BC.vacuum, BC.vacuum),
        )
        assert g.domain_extent_cm == 3.0
        assert len(g.regions) == 3

    def test_unknown_geometry_rejected(self):
        with pytest.raises(ValueError, match="must be one of"):
            StructuredGeometry(
                geometry="HEX",
                regions=(Region(mat_id=0, outer_thickness_cm=1.0),),
                bcs=(BC.vacuum,),
            )

    def test_lowercase_geometry_rejected(self):
        with pytest.raises(ValueError, match="must be one of"):
            StructuredGeometry(
                geometry="sph",
                regions=(Region(mat_id=0, outer_thickness_cm=1.0),),
                bcs=(BC.vacuum,),
            )

    def test_empty_regions_rejected(self):
        with pytest.raises(ValueError, match="must be non-empty"):
            StructuredGeometry(geometry="SPH", regions=(), bcs=(BC.vacuum,))

    def test_regions_must_be_tuple(self):
        with pytest.raises(TypeError, match="must be a tuple"):
            StructuredGeometry(
                geometry="SPH",
                regions=[Region(mat_id=0, outer_thickness_cm=1.0)],  # type: ignore[arg-type]
                bcs=(BC.vacuum,),
            )

    def test_bcs_must_be_tuple(self):
        with pytest.raises(TypeError, match="must be a tuple"):
            StructuredGeometry(
                geometry="SPH",
                regions=(Region(mat_id=0, outer_thickness_cm=1.0),),
                bcs=[BC.vacuum],  # type: ignore[arg-type]
            )

    def test_slab_requires_two_bcs(self):
        with pytest.raises(ValueError, match="requires 2 BC"):
            StructuredGeometry(
                geometry="SLB",
                regions=(Region(mat_id=0, outer_thickness_cm=4.0),),
                bcs=(BC.vacuum,),
            )

    def test_sphere_requires_one_bc(self):
        with pytest.raises(ValueError, match="requires 1 BC"):
            StructuredGeometry(
                geometry="SPH",
                regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
                bcs=(BC.vacuum, BC.vacuum),
            )

    def test_bcs_must_be_a_tag_or_a_typed_law(self):
        """RE-POSED from ``test_bcs_must_be_BC_instances``.

        A declaration is EITHER a ``BC`` tag or an already-typed
        ``BoundaryTraceLaw``; a bare string is still neither. The claim
        widened when the declaration channel landed — a law carrying a
        FUNCTION (a prescribed inflow whose source is a manufactured
        solution) has no tag spelling, and declaring it on the GEOMETRY is
        what makes it survive the method-mesh rebuild every public solver
        entry point performs.
        """
        with pytest.raises(TypeError, match="must be a BC tag or a"):
            StructuredGeometry(
                geometry="SPH",
                regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
                bcs=("vacuum",),  # type: ignore[arg-type]
            )

    def test_bcs_accepts_a_typed_law(self):
        """The positive leg — without it the guard could reject every law."""
        from orpheus.geometry.boundary import (
            ConstantInflowSource, PrescribedInflow,
        )

        law = PrescribedInflow(source=ConstantInflowSource(value=2.5))
        geom = StructuredGeometry(
            geometry="SPH",
            regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
            bcs=(law,),
        )
        assert geom.bcs[0] is law

    def test_regions_must_be_Region_instances(self):
        with pytest.raises(TypeError, match="must be a Region"):
            StructuredGeometry(
                geometry="SPH",
                regions=({"mat_id": 0, "outer_thickness_cm": 1.0},),  # type: ignore[arg-type]
                bcs=(BC.vacuum,),
            )

    def test_frozen(self):
        g = StructuredGeometry(
            geometry="SPH",
            regions=(Region(mat_id=0, outer_thickness_cm=1.0),),
            bcs=(BC.vacuum,),
        )
        with pytest.raises(AttributeError):
            g.geometry = "CYL"


class TestWignerSeitzPinCell:
    def test_default_construction(self):
        g = StructuredGeometry.wigner_seitz_pin_cell(
            r_fuel=0.9, r_clad=1.1, pitch=3.6,
        )
        assert g.geometry == "CYL"
        assert g.n_endpoints == 1
        assert g.bcs == (BC("white"),)
        assert len(g.regions) == 3

    def test_region_thicknesses(self):
        g = StructuredGeometry.wigner_seitz_pin_cell(
            r_fuel=0.9, r_clad=1.1, pitch=3.6,
        )
        # fuel: 0 → 0.9; clad: 0.9 → 1.1; cool: 1.1 → r_cell
        r_cell = 3.6 / np.sqrt(np.pi)
        assert g.regions[0].mat_id == 2
        assert g.regions[0].outer_thickness_cm == pytest.approx(0.9)
        assert g.regions[1].mat_id == 1
        assert g.regions[1].outer_thickness_cm == pytest.approx(0.2)
        assert g.regions[2].mat_id == 0
        assert g.regions[2].outer_thickness_cm == pytest.approx(r_cell - 1.1)

    def test_extent_matches_r_cell(self):
        g = StructuredGeometry.wigner_seitz_pin_cell(
            r_fuel=0.9, r_clad=1.1, pitch=3.6,
        )
        assert g.domain_extent_cm == pytest.approx(3.6 / np.sqrt(np.pi))

    def test_custom_outer_bc(self):
        g = StructuredGeometry.wigner_seitz_pin_cell(
            r_fuel=0.9, r_clad=1.1, pitch=3.6, bcs=(BC.vacuum,),
        )
        assert g.bcs == (BC.vacuum,)


# ─────────────────────────────────────────────────────────────────────
# RegionMesh — mesh-layer per-region descriptor
# ─────────────────────────────────────────────────────────────────────


class TestRegionMesh:
    def test_default_method(self):
        rm = RegionMesh(n_cells=10)
        assert rm.n_cells == 10
        assert rm.method == "equal-volume"

    def test_uniform_method(self):
        rm = RegionMesh(n_cells=5, method="uniform")
        assert rm.method == "uniform"

    def test_zero_cells_rejected(self):
        with pytest.raises(ValueError, match="must be ≥ 1"):
            RegionMesh(n_cells=0)

    def test_negative_cells_rejected(self):
        with pytest.raises(ValueError, match="must be ≥ 1"):
            RegionMesh(n_cells=-1)

    def test_unknown_method_rejected(self):
        with pytest.raises(ValueError, match="must be 'equal-volume' or 'uniform'"):
            RegionMesh(n_cells=5, method="invalid")  # type: ignore[arg-type]


# ─────────────────────────────────────────────────────────────────────
# Mesh1D.from_geometry — the canonical geometry → mesh bridge
# ─────────────────────────────────────────────────────────────────────


class TestMesh1DFromGeometry:
    def test_single_region_sphere_equal_volume(self):
        g = StructuredGeometry(
            geometry="SPH",
            regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
            bcs=(BC.vacuum,),
        )
        mesh = Mesh1D.from_geometry(g, region_meshes=(RegionMesh(n_cells=8),))
        assert mesh.N == 8
        assert mesh.coord == CoordSystem.SPHERICAL
        assert mesh.edges[0] == 0.0
        assert mesh.edges[-1] == pytest.approx(2.0)
        # All cells in an equal-volume zone are bit-identical by
        # construction (the precomputed_volumes invariant).
        assert np.all(mesh.volumes == mesh.volumes[0])
        # BC propagation: SPH → bc_right populated, bc_left None.
        assert mesh.bc_left is None
        assert mesh.bc_right == BC.vacuum

    def test_single_region_slab_uniform(self):
        g = StructuredGeometry(
            geometry="SLB",
            regions=(Region(mat_id=0, outer_thickness_cm=4.0),),
            bcs=(BC.vacuum, BC.reflective),
        )
        mesh = Mesh1D.from_geometry(
            g, region_meshes=(RegionMesh(n_cells=4, method="uniform"),),
        )
        assert mesh.N == 4
        assert mesh.coord == CoordSystem.CARTESIAN
        np.testing.assert_allclose(mesh.edges, [0.0, 1.0, 2.0, 3.0, 4.0])
        # SLB → both BCs propagated.
        assert mesh.bc_left == BC.vacuum
        assert mesh.bc_right == BC.reflective

    def test_multi_region_slab(self):
        g = StructuredGeometry(
            geometry="SLB",
            regions=(
                Region(mat_id=1, outer_thickness_cm=0.5),
                Region(mat_id=0, outer_thickness_cm=2.0),
                Region(mat_id=1, outer_thickness_cm=0.5),
            ),
            bcs=(BC.vacuum, BC.vacuum),
        )
        mesh = Mesh1D.from_geometry(g, region_meshes=(
            RegionMesh(n_cells=2, method="uniform"),
            RegionMesh(n_cells=4, method="uniform"),
            RegionMesh(n_cells=2, method="uniform"),
        ))
        assert mesh.N == 8
        assert mesh.edges[-1] == pytest.approx(3.0)
        # mat_id walks: 1 1 | 0 0 0 0 | 1 1
        np.testing.assert_array_equal(
            mesh.mat_ids, [1, 1, 0, 0, 0, 0, 1, 1],
        )

    def test_multi_region_cylinder_equal_volume(self):
        g = StructuredGeometry.wigner_seitz_pin_cell(
            r_fuel=0.9, r_clad=1.1, pitch=3.6,
        )
        mesh = Mesh1D.from_geometry(g, region_meshes=(
            RegionMesh(n_cells=10),
            RegionMesh(n_cells=3),
            RegionMesh(n_cells=7),
        ))
        assert mesh.N == 20
        assert mesh.coord == CoordSystem.CYLINDRICAL
        # Outer edge equals r_cell.
        r_cell = 3.6 / np.sqrt(np.pi)
        assert mesh.edges[-1] == pytest.approx(r_cell)
        # mat_id: 10 fuel, 3 clad, 7 cool
        assert (mesh.mat_ids == 2).sum() == 10
        assert (mesh.mat_ids == 1).sum() == 3
        assert (mesh.mat_ids == 0).sum() == 7
        assert mesh.bc_right == BC("white")
        assert mesh.bc_left is None

    def test_length_mismatch_raises(self):
        g = StructuredGeometry(
            geometry="SLB",
            regions=(
                Region(mat_id=0, outer_thickness_cm=1.0),
                Region(mat_id=1, outer_thickness_cm=1.0),
            ),
            bcs=(BC.vacuum, BC.vacuum),
        )
        with pytest.raises(ValueError, match="must equal"):
            Mesh1D.from_geometry(g, region_meshes=(RegionMesh(n_cells=4),))

    def test_origin_offset(self):
        g = StructuredGeometry(
            geometry="SLB",
            regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
            bcs=(BC.vacuum, BC.vacuum),
        )
        mesh = Mesh1D.from_geometry(
            g, region_meshes=(RegionMesh(n_cells=2, method="uniform"),),
            origin=5.0,
        )
        np.testing.assert_allclose(mesh.edges, [5.0, 6.0, 7.0])

    @pytest.mark.catches("ERR-020")
    def test_equal_volume_cylindrical_invariant(self):
        """Equal-volume cells in a cylindrical zone are bit-identical."""
        g = StructuredGeometry(
            geometry="CYL",
            regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
            bcs=(BC("white"),),
        )
        mesh = Mesh1D.from_geometry(g, region_meshes=(RegionMesh(n_cells=10),))
        # All cells exactly equal volume — no ULP drift.
        assert np.all(mesh.volumes == mesh.volumes[0])
        # Total volume from cells matches geometric formula.
        expected_total = np.pi * 2.0 ** 2
        np.testing.assert_allclose(mesh.volumes.sum(), expected_total, rtol=1e-14)

    @pytest.mark.catches("ERR-020")
    def test_equal_volume_spherical_invariant(self):
        g = StructuredGeometry(
            geometry="SPH",
            regions=(Region(mat_id=0, outer_thickness_cm=3.0),),
            bcs=(BC.vacuum,),
        )
        mesh = Mesh1D.from_geometry(g, region_meshes=(RegionMesh(n_cells=12),))
        assert np.all(mesh.volumes == mesh.volumes[0])
        expected_total = (4.0 / 3.0) * np.pi * 3.0 ** 3
        np.testing.assert_allclose(mesh.volumes.sum(), expected_total, rtol=1e-14)
