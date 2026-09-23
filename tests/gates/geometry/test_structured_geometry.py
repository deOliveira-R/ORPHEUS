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
    compute_volumes_1d,
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


# The power of ``r`` in the cell measure: ``V ∝ x_out - x_in`` with
# ``x = r**p`` (slab p=1, cylinder p=2, sphere p=3).
_MEASURE_POWER = {
    CoordSystem.CARTESIAN: 1,
    CoordSystem.CYLINDRICAL: 2,
    CoordSystem.SPHERICAL: 3,
}

# The volume (per unit transverse area / height) of the region between
# two radii: the closed form, written from geometry, not from the mesh.
_MEASURE_PREFACTOR = {
    CoordSystem.CARTESIAN: 1.0,
    CoordSystem.CYLINDRICAL: np.pi,
    CoordSystem.SPHERICAL: 4.0 / 3.0 * np.pi,
}


def _closed_form_region_volume(
    coord: CoordSystem, inner: float, outer: float,
) -> float:
    p = _MEASURE_POWER[coord]
    return _MEASURE_PREFACTOR[coord] * (outer**p - inner**p)


def _assert_equal_volume_regions(
    mesh: Mesh1D,
    mat_ids: tuple[int, ...],
    cells_per_region: tuple[int, ...],
    radii: tuple[float, ...],
) -> None:
    """Per region: every cell volume EXACTLY equal, the total the closed form.

    The ERR-020 invariant, asserted region by region. ``radii`` are the
    region boundaries written by the caller from the geometry (never
    read off ``mesh.edges``). Each region's slice is first asserted to
    BE that region (its cells carry its material id), so a change of
    cell ordering fails here rather than silently checking the wrong
    cells.
    """
    np.testing.assert_equal(sum(cells_per_region), mesh.N)
    start = 0
    for k, (mat_id, n) in enumerate(zip(mat_ids, cells_per_region, strict=True)):
        cells = slice(start, start + n)
        np.testing.assert_array_equal(
            mesh.mat_ids[cells], np.full(n, mat_id),
            err_msg=f"region {k}: the cell slice is not the region",
        )
        volumes = mesh.volumes[cells]
        np.testing.assert_array_equal(
            volumes, np.full(n, volumes[0]),
            err_msg=(
                f"region {k} ({mesh.coord.name}, r in "
                f"[{radii[k]}, {radii[k + 1]}]): equal-volume cells are "
                f"not bit-identical — the volumes were re-derived from "
                f"the edges (ERR-020 round trip)"
            ),
        )
        np.testing.assert_allclose(
            volumes.sum(),
            _closed_form_region_volume(mesh.coord, radii[k], radii[k + 1]),
            rtol=1e-14,
            err_msg=f"region {k} ({mesh.coord.name}): total volume",
        )
        start += n


# Three regions whose inner radius changes at every interface, meshed
# with cell counts that are NOT powers of two. The counts matter for
# the Cartesian row: on a power-of-two subdivision of these dyadic
# thicknesses every edge is exact, so ``np.diff(edges)`` reproduces
# ``(outer - inner) / n`` bit for bit and the round trip cannot show
# ([M] 2026-09-22: counts (4, 8, 4) -> 0 of 16 cells differ under the
# round trip; counts (5, 7, 11) -> 3 of 5, 3 of 7, 2 of 11).
_THREE_REGION_MAT_IDS = (0, 1, 2)
_THREE_REGION_THICKNESS = (0.5, 1.0, 0.5)
_THREE_REGION_RADII = (0.0, 0.5, 1.5, 2.0)
_THREE_REGION_CELLS = (5, 7, 11)


def _three_region_mesh(geometry: str) -> Mesh1D:
    """The three-region mesh, built the way production builds one."""
    g = StructuredGeometry(
        geometry=geometry,
        regions=tuple(
            Region(mat_id=m, outer_thickness_cm=t)
            for m, t in zip(
                _THREE_REGION_MAT_IDS, _THREE_REGION_THICKNESS, strict=True,
            )
        ),
        bcs=(BC.vacuum, BC.vacuum) if geometry == "SLB" else (BC.reflective,),
    )
    return Mesh1D.from_geometry(g, region_meshes=tuple(
        RegionMesh(n_cells=n) for n in _THREE_REGION_CELLS
    ))


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

    @pytest.mark.catches("ERR-020")
    def test_multi_region_cylinder_equal_volume(self):
        """The production pin-cell mesh: equal volume within each region.

        The Wigner-Seitz pin cell meshed 10 / 3 / 7 is the default mesh
        of the CP and MoC solvers. Beyond its cell counts, material ids
        and boundary conditions, the name's claim is asserted: within
        each of the three regions the cell volumes are EXACTLY equal and
        sum to the closed-form annulus volume (ERR-020, multi-region).
        """
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
        _assert_equal_volume_regions(
            mesh,
            mat_ids=(2, 1, 0),
            cells_per_region=(10, 3, 7),
            radii=(0.0, 0.9, 1.1, r_cell),
        )

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

    @pytest.mark.catches("ERR-020")
    @pytest.mark.parametrize("geometry", ["SLB", "CYL", "SPH"])
    def test_equal_volume_multi_region_invariant(self, geometry):
        """Equal-volume cells are bit-identical within EACH of three regions.

        A multi-region mesh puts a subdivision boundary at every region
        interface, where the inner radius of the equal-volume
        subdivision changes (issue #489). Per region the cell volumes
        are asserted EXACTLY equal (``==``, not a tolerance) and their
        sum equal to the closed-form region volume at ``rtol=1e-14``.

        The two legs catch different defects. The equality leg catches
        ERR-020 itself: volumes re-derived from the edges through the
        ``sqrt``-then-square / ``cbrt``-then-cube round trip (and, on
        the slab, through the rounding of the edge positions), which
        leaves every region's total correct, since the differences
        telescope. The total leg catches a per-cell volume that drops
        the region's inner radius, ``V_cell = V(0, outer) / n``: exact
        on the first region and on every single-region mesh, whose
        inner radius is the origin, so only a region beyond the first
        can witness it.
        """
        mesh = _three_region_mesh(geometry)
        _assert_equal_volume_regions(
            mesh,
            mat_ids=_THREE_REGION_MAT_IDS,
            cells_per_region=_THREE_REGION_CELLS,
            radii=_THREE_REGION_RADII,
        )

    @pytest.mark.parametrize("geometry", ["SLB", "CYL", "SPH"])
    def test_equal_volume_edges_bound_the_volumes(self, geometry):
        """The equal-volume edges and the precomputed volumes are one mesh.

        ``Mesh1D.from_geometry`` stores the equal-volume cell volumes
        computed from the algebraic invariant (the ERR-020 fix), so the
        volumes no longer read the edges, and no volume assertion can
        see an error in the equal-volume RADIUS formula
        ``r_k = (r_in^p + k/n (r_out^p - r_in^p))^(1/p)``. This row is
        the witness that the two spellings of the mesh agree: the volume
        each pair of edges bounds, re-derived through the round trip,
        equals the stored volume to the round trip's own conditioning.

        Tolerance, derived rather than chosen: a cell's re-derived
        volume is a difference ``x_{k+1} - x_k`` of ``x = r^p`` values
        each carrying O(eps) relative error, amplified by
        ``x / Δx <= n x_out / (x_out - x_in)``. The measured worst case
        is 3.3 of that unit over 36 configurations (three coordinate
        systems, four thickness sets, three cell-count sets up to 1000
        cells; [M] 2026-09-22), so a factor 8 leaves 2.4x headroom while
        an O(1) error in the radius formula is red by orders of
        magnitude.
        """
        mesh = _three_region_mesh(geometry)
        p = _MEASURE_POWER[mesh.coord]
        eps = np.finfo(float).eps
        rederived = compute_volumes_1d(mesh.coord, mesh.edges)
        start = 0
        for k, n in enumerate(_THREE_REGION_CELLS):
            cells = slice(start, start + n)
            x_in = _THREE_REGION_RADII[k] ** p
            x_out = _THREE_REGION_RADII[k + 1] ** p
            np.testing.assert_allclose(
                rederived[cells], mesh.volumes[cells],
                rtol=8.0 * n * x_out / (x_out - x_in) * eps, atol=0.0,
                err_msg=(
                    f"region {k} ({mesh.coord.name}): the edges do not "
                    f"bound the stored equal volumes"
                ),
            )
            start += n
