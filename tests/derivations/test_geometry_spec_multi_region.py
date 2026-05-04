r"""Foundation gates for :class:`GeometrySpec` multi-region support.

The 2026-05-04 Step-3 promotion added a :class:`Region` dataclass
and a :attr:`GeometrySpec.regions` field for first-class multi-
region geometry. These tests pin the contract of that API:

* Single-region default (``regions=None``) preserves the original
  behaviour bit-for-bit.
* Multi-region slab build emits a multi-zone :class:`Mesh1D` with
  the right ``mat_ids``, edges, and BCs.
* Multi-region sphere build is sane (inside-out radial layering,
  outer radius equals the thickness sum).
* Validation errors fire on every malformed-region path:

  - Region missing both thickness fields,
  - Region carrying ``outer_thickness_cm <= 0``,
  - GeometrySpec.regions on infinite / ISLC geometries,
  - GeometrySpec.regions empty tuple,
  - GeometrySpec.regions missing ``outer_thickness_cm``.

The mathematical content is trivial (mesh edges and material IDs);
foundation-tier tests here pin the **API contract** so future
refactors don't silently regress the multi-region path. No
``verifies(...)`` decorator: these are software invariants on a
descriptor, not theory-equation claims.

References
----------

* :class:`orpheus.derivations.common.geometry_spec.GeometrySpec`
* :class:`orpheus.derivations.common.geometry_spec.Region`
* :func:`orpheus.geometry.factories.mesh1d_from_zones`
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.geometry_spec import GeometrySpec, Region
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.mesh import BC

pytestmark = [pytest.mark.foundation]


# ── Single-region default: the regions=None path is unchanged ──────


def test_regions_none_is_default():
    """``GeometrySpec.regions`` defaults to ``None``."""
    spec = GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=1.0,
        critical_dimension_cm=1.0,
        n_groups=1,
    )
    assert spec.regions is None


def test_single_region_slab_build_unchanged():
    """``regions=None`` slab still builds a single-zone mesh.

    Pinned: domain extent is ``2 * critical_dimension_cm`` (the
    F_N full-slab convention) and every cell has ``mat_id``
    equal to the descriptor's ``mat_id`` field.
    """
    spec = GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=0.93772556,
        critical_dimension_cm=0.93772556,
        n_groups=1,
        mat_id=7,
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    )
    assert spec.domain_extent_cm == pytest.approx(2.0 * 0.93772556)
    mesh = spec.build(n_cells=16)
    assert len(mesh.mat_ids) == 16
    assert np.all(mesh.mat_ids == 7)
    assert mesh.coord == CoordSystem.CARTESIAN
    assert mesh.bc_left == BC.vacuum
    assert mesh.bc_right == BC.vacuum
    assert mesh.edges[0] == pytest.approx(0.0)
    assert mesh.edges[-1] == pytest.approx(2.0 * 0.93772556)


def test_single_region_sphere_build_unchanged():
    """``regions=None`` sphere still builds ``[0, R]`` with
    homogeneous mat_id and the centreline-reflective BC convention.
    """
    spec = GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=2.4248249802,
        critical_dimension_cm=2.4248249802,
        n_groups=1,
        mat_id=3,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    )
    mesh = spec.build(n_cells=8)
    assert len(mesh.mat_ids) == 8
    assert np.all(mesh.mat_ids == 3)
    assert mesh.coord == CoordSystem.SPHERICAL
    assert mesh.bc_left == BC.reflective
    assert mesh.bc_right == BC.vacuum
    assert mesh.edges[0] == pytest.approx(0.0)
    assert mesh.edges[-1] == pytest.approx(2.4248249802)


# ── Multi-region slab: NM 1980 / Sood Table 10 layout ───────────────


def test_multi_region_slab_extent_is_thickness_sum():
    """For multi-region slabs, ``domain_extent_cm`` is the literal
    sum of region cm thicknesses (NOT ``2 * critical_dimension_cm``).
    """
    regions = (
        Region(mat_id=1, outer_thickness_mfp=0.5, outer_thickness_cm=0.5),
        Region(mat_id=0, outer_thickness_mfp=2.0, outer_thickness_cm=2.0),
        Region(mat_id=1, outer_thickness_mfp=0.5, outer_thickness_cm=0.5),
    )
    spec = GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=1.0,  # core half-thickness, ignored by extent
        critical_dimension_cm=1.0,
        n_groups=1,
        regions=regions,
    )
    assert spec.domain_extent_cm == pytest.approx(0.5 + 2.0 + 0.5)


def test_multi_region_slab_build_zones_have_correct_mat_ids():
    """Multi-region slab build: layer ordering is left-to-right and
    each zone's cells carry the right ``mat_id``.

    Layout: [reflector(1) | core(0) | reflector(1)] with thicknesses
    0.5 / 2.0 / 0.5 cm. With ``n_cells=20`` the proportional
    allocation gives 2/16/2 cells (rounding pushes residual onto
    the largest region).
    """
    regions = (
        Region(mat_id=1, outer_thickness_cm=0.5),
        Region(mat_id=0, outer_thickness_cm=2.0),
        Region(mat_id=1, outer_thickness_cm=0.5),
    )
    spec = GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=1.0,
        critical_dimension_cm=1.0,
        n_groups=1,
        regions=regions,
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    )
    mesh = spec.build(n_cells=20)
    # Total cells = n_cells when proportional alloc + residual is clean.
    assert len(mesh.mat_ids) == 20
    # First and last segments carry mat_id=1 (reflector), middle is core.
    assert mesh.mat_ids[0] == 1
    assert mesh.mat_ids[-1] == 1
    # The largest run of identical mat_ids is the core (mat_id=0).
    core_mask = mesh.mat_ids == 0
    assert core_mask.any(), "core region should populate some cells"
    # Edges should span [0, 3.0] with internal breakpoints at 0.5 and 2.5.
    assert mesh.edges[0] == pytest.approx(0.0)
    assert mesh.edges[-1] == pytest.approx(3.0)
    # Verify the layer interfaces are present in the edge array.
    assert any(np.isclose(mesh.edges, 0.5)), "interface at 0.5 cm missing"
    assert any(np.isclose(mesh.edges, 2.5)), "interface at 2.5 cm missing"


def test_multi_region_slab_n_cells_floor_one_per_region():
    """When ``n_cells < len(regions)``, every region still gets at
    least one cell (the floor)."""
    regions = (
        Region(mat_id=2, outer_thickness_cm=0.1),
        Region(mat_id=0, outer_thickness_cm=10.0),
        Region(mat_id=2, outer_thickness_cm=0.1),
    )
    spec = GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=5.0,
        critical_dimension_cm=5.0,
        n_groups=1,
        regions=regions,
    )
    mesh = spec.build(n_cells=2)  # less than len(regions)
    # Each region has ≥1 cell → at least 3 cells total.
    assert len(mesh.mat_ids) >= 3
    # Both reflector layers must appear.
    assert (mesh.mat_ids == 2).sum() >= 2
    # Core must appear.
    assert (mesh.mat_ids == 0).sum() >= 1


# ── Multi-region sphere ────────────────────────────────────────────


def test_multi_region_sphere_outer_radius_is_thickness_sum():
    """Sphere multi-region: ``domain_extent_cm`` (= outer radius) is
    the inside-out thickness sum.

    Layout: 1 cm fuel core + 0.5 cm reflector → outer R = 1.5 cm.
    """
    regions = (
        Region(mat_id=0, outer_thickness_cm=1.0),
        Region(mat_id=1, outer_thickness_cm=0.5),
    )
    spec = GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=1.0,
        critical_dimension_cm=1.0,
        n_groups=1,
        regions=regions,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    )
    assert spec.domain_extent_cm == pytest.approx(1.5)
    mesh = spec.build(n_cells=12)
    assert mesh.coord == CoordSystem.SPHERICAL
    assert mesh.edges[0] == pytest.approx(0.0)
    assert mesh.edges[-1] == pytest.approx(1.5)
    # Interface at r=1.0 must appear.
    assert any(np.isclose(mesh.edges, 1.0))
    # Inside-out ordering: first cell is mat 0 (core), last is mat 1.
    assert mesh.mat_ids[0] == 0
    assert mesh.mat_ids[-1] == 1


# ── Region dataclass validation ────────────────────────────────────


def test_region_requires_one_thickness_field():
    """A :class:`Region` with neither ``outer_thickness_mfp`` nor
    ``outer_thickness_cm`` set must raise.
    """
    with pytest.raises(ValueError, match="outer_thickness_mfp or"):
        Region(mat_id=0)


def test_region_rejects_nonpositive_thickness():
    """Negative or zero thickness in either unit is rejected."""
    with pytest.raises(ValueError, match="outer_thickness_cm must be > 0"):
        Region(mat_id=0, outer_thickness_cm=0.0)
    with pytest.raises(ValueError, match="outer_thickness_cm must be > 0"):
        Region(mat_id=0, outer_thickness_cm=-1.0)
    with pytest.raises(ValueError, match="outer_thickness_mfp must be > 0"):
        Region(mat_id=0, outer_thickness_mfp=-0.5, outer_thickness_cm=0.5)


def test_region_mfp_only_is_allowed_at_construction():
    """A :class:`Region` carrying only the mfp form constructs OK
    (the cm form gate fires only at GeometrySpec build time)."""
    # Should not raise.
    r = Region(mat_id=0, outer_thickness_mfp=1.0)
    assert r.outer_thickness_cm is None


# ── GeometrySpec.regions validation ───────────────────────────────


def test_geometryspec_regions_empty_raises():
    """``regions=()`` (set but empty) is rejected."""
    with pytest.raises(ValueError, match="non-empty"):
        GeometrySpec(
            geometry="slab",
            critical_dimension_mfp=1.0,
            critical_dimension_cm=1.0,
            n_groups=1,
            regions=(),
        )


def test_geometryspec_regions_missing_cm_raises():
    """Region without ``outer_thickness_cm`` is rejected at the
    spec level (the mesh lives in cm)."""
    regions = (
        Region(mat_id=0, outer_thickness_mfp=1.0),  # no cm!
    )
    with pytest.raises(ValueError, match="missing outer_thickness_cm"):
        GeometrySpec(
            geometry="slab",
            critical_dimension_mfp=1.0,
            critical_dimension_cm=1.0,
            n_groups=1,
            regions=regions,
        )


def test_geometryspec_regions_rejected_for_infinite():
    """Infinite-medium specs cannot carry regions (no spatial mesh)."""
    regions = (Region(mat_id=0, outer_thickness_cm=1.0),)
    with pytest.raises(ValueError, match="infinite geometry"):
        GeometrySpec(
            geometry="infinite",
            critical_dimension_mfp=None,
            critical_dimension_cm=None,
            n_groups=1,
            regions=regions,
        )


def test_geometryspec_regions_rejected_for_islc():
    """ISLC geometry is currently outside the multi-region path."""
    regions = (Region(mat_id=0, outer_thickness_cm=1.0),)
    with pytest.raises(ValueError, match="not supported"):
        GeometrySpec(
            geometry="ISLC",
            critical_dimension_mfp=1.0,
            critical_dimension_cm=1.0,
            n_groups=1,
            regions=regions,
        )


def test_geometryspec_multi_region_n_cells_total_matches_budget():
    """With proportional allocation + residual reconciliation the
    total cell count equals the requested budget when
    ``n_cells >= len(regions)``."""
    regions = (
        Region(mat_id=1, outer_thickness_cm=0.5),
        Region(mat_id=0, outer_thickness_cm=2.0),
        Region(mat_id=1, outer_thickness_cm=0.5),
    )
    spec = GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=1.0,
        critical_dimension_cm=1.0,
        n_groups=1,
        regions=regions,
    )
    for n in [3, 7, 16, 31, 64, 100]:
        mesh = spec.build(n_cells=n)
        assert len(mesh.mat_ids) == n, (
            f"n_cells={n}: expected total {n}, got {len(mesh.mat_ids)}"
        )
