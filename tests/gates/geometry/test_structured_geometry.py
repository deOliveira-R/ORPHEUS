"""Foundation tests for :class:`StructuredGeometry` and :meth:`Mesh1D.from_geometry`.

These tests pin the geometry value and the geometry → mesh transition: the
defining refusals of the value (its breakpoints, its material ids, its
derived boundary), the thickness constructor, the retirement of the string
kind tag, and the construction semantics of the mesh built from it. They are
foundation-tier (software invariants): they verify the code-shape contract,
not a physics claim. The gate ids S2.1 to S2.5 are those of the P1
verification specification (``.claude/plans/reference_p1_spec.md`` §1.2).
"""
from __future__ import annotations

import itertools
import math
from typing import Any

import numpy as np
import pytest

import orpheus.geometry
import orpheus.geometry.structured_geometry as structured_geometry_module
from orpheus.geometry import (
    BC,
    CoordSystem,
    StructuredGeometry,
    compute_volumes_1d,
)
from orpheus.mesh import Mesh1D, RegionMesh


pytestmark = pytest.mark.foundation

_SLAB = CoordSystem.CARTESIAN
_CYLINDER = CoordSystem.CYLINDRICAL
_SPHERE = CoordSystem.SPHERICAL
_CURVILINEAR = (_CYLINDER, _SPHERE)


def _geometry(**overrides: Any) -> StructuredGeometry:
    """A valid two-interval slab, with the named fields replaced."""
    fields: dict[str, Any] = dict(
        coord=_SLAB,
        breakpoints=(0.0, 1.0, 2.0),
        mat_ids=(0, 1),
        boundaries=(BC.vacuum, BC.vacuum),
    )
    fields.update(overrides)
    return StructuredGeometry(**fields)


# ─────────────────────────────────────────────────────────────────────
# S2.1 — the breakpoint laws: the value's defining refusals
# ─────────────────────────────────────────────────────────────────────

#: (case id, overrides, error type, the fragment that keys the refusal).
#: One row per clause of the constructor; S2.1's disjointness row reads
#: every message against every other row's fragment.
_BREAKPOINT_REFUSALS = [
    ("fewer-than-two", dict(breakpoints=(1.0,), mat_ids=()),
     ValueError, "at least 2 breakpoints"),
    ("equal-pair", dict(breakpoints=(0.0, 1.0, 1.0)),
     ValueError, "strictly increasing"),
    ("decreasing-pair", dict(breakpoints=(0.0, 2.0, 1.0)),
     ValueError, "strictly increasing"),
    ("nan", dict(breakpoints=(0.0, math.nan, 2.0)),
     ValueError, "must be finite"),
    ("inf", dict(breakpoints=(0.0, 1.0, math.inf)),
     ValueError, "must be finite"),
    ("negative-radius-cylinder",
     dict(coord=_CYLINDER, breakpoints=(-0.5, 1.0), mat_ids=(0,)),
     ValueError, "starts at r_0 >= 0"),
    ("negative-radius-sphere",
     dict(coord=_SPHERE, breakpoints=(-0.5, 1.0), mat_ids=(0,)),
     ValueError, "starts at r_0 >= 0"),
    ("material-id-count", dict(mat_ids=(0,)),
     ValueError, "one material id per interval"),
    ("material-id-not-int", dict(mat_ids=(0, 1.0)),
     TypeError, "a material id is an int"),
]


class TestBreakpointLaws:
    """S2.1: the breakpoints and material ids a geometry admits."""

    @pytest.mark.parametrize(
        "overrides, error, fragment",
        [pytest.param(o, e, f, id=i) for i, o, e, f in _BREAKPOINT_REFUSALS],
    )
    def test_refusal(self, overrides, error, fragment):
        with pytest.raises(error, match=fragment):
            _geometry(**overrides)

    def test_refusal_fragments_are_disjoint(self):
        """Each refusal's message carries its own fragment and no other's."""
        messages = {}
        for case, overrides, error, _ in _BREAKPOINT_REFUSALS:
            with pytest.raises(error) as caught:
                _geometry(**overrides)
            messages[case] = str(caught.value)
        for case, _, _, fragment in _BREAKPOINT_REFUSALS:
            for other, message in messages.items():
                own = [f for c, _, _, f in _BREAKPOINT_REFUSALS if c == other]
                if fragment in own:
                    continue
                assert fragment not in message, (
                    f"the fragment of {case!r} ({fragment!r}) also appears in "
                    f"the refusal of {other!r}: {message!r}"
                )

    def test_a_slab_admits_any_origin(self):
        g = _geometry(breakpoints=(-1.0, 0.0, 2.0))
        assert g.breakpoints == (-1.0, 0.0, 2.0)

    def test_a_sphere_may_be_hollow(self):
        g = StructuredGeometry(
            coord=_SPHERE,
            breakpoints=(0.1, 1.0),
            mat_ids=(0,),
            boundaries=(BC.reflective, BC.vacuum),
        )
        assert g.breakpoints[0] == 0.1

    def test_adjacent_intervals_may_share_a_material(self):
        g = _geometry(mat_ids=(0, 0))
        assert g.mat_ids == (0, 0)

    def test_entries_are_canonicalised_to_float_and_int(self):
        """Integers and numpy scalars are parsed at the boundary, bit for bit."""
        g = _geometry(breakpoints=(0, np.float64(1.0), 2), mat_ids=(np.int64(3), 1))
        assert all(type(b) is float for b in g.breakpoints)
        assert all(type(m) is int for m in g.mat_ids)
        assert g.breakpoints == (0.0, 1.0, 2.0)
        assert g.mat_ids == (3, 1)

    def test_frozen(self):
        g = _geometry()
        with pytest.raises(AttributeError):
            g.coord = _SPHERE  # type: ignore[misc]

    def test_fields_are_keyword_only(self):
        """Two adjacent tuples of numbers cannot be swapped positionally."""
        with pytest.raises(TypeError):
            StructuredGeometry(_SLAB, (0.0, 1.0), (0,), (BC.vacuum, BC.vacuum))  # type: ignore[misc]


# ─────────────────────────────────────────────────────────────────────
# S2.2 — the boundary is derived from the coordinate system and r_0
# ─────────────────────────────────────────────────────────────────────

#: (coord, r_0) → the boundary points of the geometry [r_0, 2.0].
_BOUNDARY_TABLE = [
    (_SLAB, 0.0, (0.0, 2.0)),
    (_SLAB, 0.5, (0.5, 2.0)),
    (_CYLINDER, 0.0, (2.0,)),
    (_CYLINDER, 0.5, (0.5, 2.0)),
    (_SPHERE, 0.0, (2.0,)),
    (_SPHERE, 0.5, (0.5, 2.0)),
]


class TestTheBoundaryIsDerived:
    """S2.2: one law per point of the region's topological boundary."""

    @pytest.mark.parametrize(
        "coord, r_0, points", _BOUNDARY_TABLE,
        ids=[f"{c.name.lower()}-r0={r}" for c, r, _ in _BOUNDARY_TABLE],
    )
    def test_boundary_points(self, coord, r_0, points):
        laws = tuple(BC.vacuum for _ in points)
        g = StructuredGeometry(
            coord=coord, breakpoints=(r_0, 2.0), mat_ids=(0,), boundaries=laws,
        )
        assert g.boundary_points == points
        assert len(g.boundaries) == len(points)

    def test_a_law_at_the_centre_is_refused(self):
        """A DISCRIMINATION row: the old table refused this input too, for
        another reason ('requires 1 BC'); the new message names the centre."""
        with pytest.raises(ValueError) as caught:
            StructuredGeometry(
                coord=_SPHERE, breakpoints=(0.0, 2.0), mat_ids=(0,),
                boundaries=(BC.vacuum, BC.vacuum),
            )
        assert "the centre r = 0" in str(caught.value)
        assert "requires 1 BC" not in str(caught.value)

    @pytest.mark.parametrize("coord", _CURVILINEAR, ids=lambda c: c.name.lower())
    def test_a_hollow_body_without_an_inner_law_is_refused(self, coord):
        with pytest.raises(ValueError, match="inner surface, which needs its own law"):
            StructuredGeometry(
                coord=coord, breakpoints=(0.5, 2.0), mat_ids=(0,),
                boundaries=(BC.vacuum,),
            )

    def test_a_slab_with_one_law_is_refused(self):
        with pytest.raises(ValueError, match="a slab has two boundary points"):
            _geometry(boundaries=(BC.vacuum,))

    def test_none_is_not_a_boundary_law(self):
        with pytest.raises(TypeError, match="None is not a boundary law"):
            _geometry(boundaries=(None, BC.vacuum))

    def test_boundaries_are_parsed_like_the_other_fields(self):
        """A list is canonicalised to a tuple; a non-sequence is refused."""
        assert _geometry(boundaries=[BC.vacuum, BC.reflective]).boundaries == (
            BC.vacuum, BC.reflective,
        )
        with pytest.raises(TypeError, match="boundaries must be a sequence"):
            _geometry(boundaries=BC.vacuum)

    def test_a_law_is_a_tag_or_a_typed_law(self):
        """RE-POSED from ``test_bcs_must_be_a_tag_or_a_typed_law``.

        A declaration is EITHER a ``BC`` tag or an already-typed
        ``BoundaryTraceLaw``; a bare string is neither. A law carrying a
        FUNCTION (a prescribed inflow whose source is a manufactured
        solution) has no tag spelling, and declaring it on the GEOMETRY is
        what makes it survive the method-mesh rebuild every public solver
        entry point performs.
        """
        with pytest.raises(TypeError, match="must be a BC tag or a"):
            _geometry(boundaries=("vacuum", BC.vacuum))

    def test_a_typed_law_is_accepted(self):
        """The positive leg; without it the guard could reject every law."""
        from orpheus.geometry.boundary import (
            ConstantInflowSource, PrescribedInflow,
        )

        law = PrescribedInflow(source=ConstantInflowSource(value=2.5))
        g = StructuredGeometry(
            coord=_SPHERE, breakpoints=(0.0, 2.0), mat_ids=(0,), boundaries=(law,),
        )
        assert g.boundaries[0] is law


# ─────────────────────────────────────────────────────────────────────
# S2.3, S2.4 — breakpoints are stored; thicknesses are folded once
# ─────────────────────────────────────────────────────────────────────

#: A thickness list on which the sequential fold and a correctly rounded
#: cumulative sum differ (at the fourth breakpoint), so a constructor that
#: re-associated the sum would be seen.
_ASSOCIATION_SENSITIVE = (0.1, 0.2, 0.3, 0.4, 0.5)


class TestFromThicknesses:
    """S2.3: the thickness constructor is the left fold, bitwise."""

    def test_the_input_discriminates_the_association(self):
        """Positive control: on this input the two sums differ."""
        sequential = tuple(itertools.accumulate(_ASSOCIATION_SENSITIVE, initial=0.0))
        rounded = tuple(
            math.fsum(_ASSOCIATION_SENSITIVE[:k])
            for k in range(len(_ASSOCIATION_SENSITIVE) + 1)
        )
        assert sequential != rounded

    @pytest.mark.parametrize("r_0", [0.0, 0.25, -1.5])
    def test_breakpoints_are_the_left_fold(self, r_0):
        g = StructuredGeometry.from_thicknesses(
            coord=_SLAB,
            thicknesses=_ASSOCIATION_SENSITIVE,
            mat_ids=range(5),
            boundaries=(BC.vacuum, BC.vacuum),
            r_0=r_0,
        )
        folded = [r_0]
        for t in _ASSOCIATION_SENSITIVE:
            folded.append(folded[-1] + t)
        assert g.breakpoints == tuple(folded)
        assert g.mat_ids == (0, 1, 2, 3, 4)

    @pytest.mark.parametrize("thickness", [0.0, -1.0], ids=["zero", "negative"])
    def test_a_non_positive_thickness_is_refused(self, thickness):
        """RE-POSED from ``TestRegion::test_zero_thickness_rejected`` and
        ``test_negative_thickness_rejected``: a thickness <= 0 is a
        non-increasing breakpoint pair."""
        with pytest.raises(ValueError, match="strictly increasing"):
            StructuredGeometry.from_thicknesses(
                coord=_SLAB, thicknesses=(1.0, thickness), mat_ids=(0, 1),
                boundaries=(BC.vacuum, BC.vacuum),
            )


class TestFromThicknessesParsesLikeTheConstructor:
    """S1 of the elegance review: the thickness constructor admits exactly
    the scalars the constructor admits (no coercion of a string or a bool)."""

    @pytest.mark.parametrize(
        "thicknesses, r_0",
        [(("0.5", "2.0"), 0.0), ((True, True), 0.0), ((0.5, 2.0), "0.0")],
        ids=["string-thickness", "bool-thickness", "string-origin"],
    )
    def test_a_non_real_scalar_is_refused(self, thicknesses, r_0):
        with pytest.raises(TypeError, match="must be a real number"):
            StructuredGeometry.from_thicknesses(
                coord=_SLAB, thicknesses=thicknesses, mat_ids=(0, 1),
                boundaries=(BC.vacuum, BC.vacuum), r_0=r_0,
            )

    def test_negative_zero_is_canonicalised(self):
        """-0.0 and 0.0 are one breakpoint; the stored bits are +0.0."""
        g = _geometry(breakpoints=(-0.0, 1.0, 2.0))
        assert math.copysign(1.0, g.breakpoints[0]) == 1.0


class TestBreakpointsAreStored:
    """S2.4: a breakpoint is never re-derived from widths."""

    _E = (0.0, 0.17, 0.45, 0.62, 1.0)

    def test_the_input_discriminates_a_re_derivation(self):
        """Positive control: re-adding the widths misses two breakpoints by
        one ULP (the census's site ``test_g_adjoint_reciprocity.py:244``)."""
        widths = np.diff(self._E).tolist()
        assert tuple(itertools.accumulate(widths, initial=0.0)) != self._E

    def test_breakpoints_round_trip_bitwise(self):
        g = StructuredGeometry(
            coord=_SLAB, breakpoints=self._E, mat_ids=(0, 1, 2, 3),
            boundaries=(BC.vacuum, BC.vacuum),
        )
        assert g.breakpoints == self._E


# ─────────────────────────────────────────────────────────────────────
# S2.5 — the string kind tag retires
# ─────────────────────────────────────────────────────────────────────


class TestTheKindTagRetires:
    @pytest.mark.parametrize("tag", ["SLB", "CYL", "SPH", "HEX", "sph"])
    def test_a_string_coordinate_is_refused(self, tag):
        with pytest.raises(TypeError, match="must be a CoordSystem member"):
            _geometry(coord=tag)

    @pytest.mark.parametrize(
        "name", ["_GEOMETRY_TO_COORD", "_GEOMETRY_TO_N_ENDPOINTS", "Region"],
    )
    def test_the_tag_machinery_is_gone(self, name):
        assert not hasattr(structured_geometry_module, name)
        assert not hasattr(orpheus.geometry, name)

    @pytest.mark.parametrize("name", ["geometry", "regions", "bcs", "n_endpoints"])
    def test_the_old_fields_are_gone(self, name):
        assert not hasattr(_geometry(), name)


# ─────────────────────────────────────────────────────────────────────
# The factories
# ─────────────────────────────────────────────────────────────────────


class TestWignerSeitzPinCell:
    def test_default_construction(self):
        g = StructuredGeometry.wigner_seitz_pin_cell(
            r_fuel=0.9, r_clad=1.1, pitch=3.6,
        )
        assert g.coord is _CYLINDER
        assert g.boundary_points == (g.breakpoints[-1],)
        assert g.boundaries == (BC("white"),)
        assert g.mat_ids == (2, 1, 0)

    def test_breakpoints_are_the_radii(self):
        """The factory states radii, so the radii are the breakpoints."""
        g = StructuredGeometry.wigner_seitz_pin_cell(
            r_fuel=0.9, r_clad=1.1, pitch=3.6,
        )
        assert g.breakpoints == (0.0, 0.9, 1.1, float(3.6 / np.sqrt(np.pi)))

    def test_extent_is_the_cell_radius(self):
        g = StructuredGeometry.wigner_seitz_pin_cell(
            r_fuel=0.9, r_clad=1.1, pitch=3.6,
        )
        assert g.domain_extent_cm == float(3.6 / np.sqrt(np.pi))

    def test_custom_outer_law(self):
        g = StructuredGeometry.wigner_seitz_pin_cell(
            r_fuel=0.9, r_clad=1.1, pitch=3.6, boundaries=(BC.vacuum,),
        )
        assert g.boundaries == (BC.vacuum,)


class TestPwrSlabHalfCell:
    def test_is_the_thickness_fold(self):
        g = StructuredGeometry.pwr_slab_half_cell(
            fuel_half=0.9, clad_thick=0.2, cool_thick=0.7,
        )
        assert g == StructuredGeometry.from_thicknesses(
            coord=_SLAB, thicknesses=(0.9, 0.2, 0.7), mat_ids=(2, 1, 0),
            boundaries=(BC.reflective, BC.reflective),
        )


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


def _three_region_mesh(coord: CoordSystem) -> Mesh1D:
    """The three-region mesh, built the way production builds one."""
    g = StructuredGeometry.from_thicknesses(
        coord=coord,
        thicknesses=_THREE_REGION_THICKNESS,
        mat_ids=_THREE_REGION_MAT_IDS,
        boundaries=(
            (BC.vacuum, BC.vacuum) if coord is _SLAB else (BC.reflective,)
        ),
    )
    return Mesh1D.from_geometry(g, region_meshes=tuple(
        RegionMesh(n_cells=n) for n in _THREE_REGION_CELLS
    ))


class TestMesh1DFromGeometry:
    def test_single_region_sphere_equal_volume(self):
        g = StructuredGeometry(
            coord=CoordSystem.SPHERICAL,
            breakpoints=(0.0, 2.0),
            mat_ids=(0,),
            boundaries=(BC.vacuum,),
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
            coord=CoordSystem.CARTESIAN,
            breakpoints=(0.0, 4.0),
            mat_ids=(0,),
            boundaries=(BC.vacuum, BC.reflective),
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
        g = StructuredGeometry.from_thicknesses(
            coord=CoordSystem.CARTESIAN,
            thicknesses=(0.5, 2.0, 0.5),
            mat_ids=(1, 0, 1),
            boundaries=(BC.vacuum, BC.vacuum),
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
        g = StructuredGeometry.from_thicknesses(
            coord=CoordSystem.CARTESIAN,
            thicknesses=(1.0, 1.0),
            mat_ids=(0, 1),
            boundaries=(BC.vacuum, BC.vacuum),
        )
        with pytest.raises(ValueError, match="must equal"):
            Mesh1D.from_geometry(g, region_meshes=(RegionMesh(n_cells=4),))

    def test_the_first_breakpoint_is_the_origin(self):
        """RE-POSED from ``test_origin_offset``: ``origin=`` retired, since
        the geometry's first breakpoint states where the mesh starts."""
        g = StructuredGeometry(
            coord=CoordSystem.CARTESIAN,
            breakpoints=(5.0, 7.0),
            mat_ids=(0,),
            boundaries=(BC.vacuum, BC.vacuum),
        )
        mesh = Mesh1D.from_geometry(
            g, region_meshes=(RegionMesh(n_cells=2, method="uniform"),),
        )
        np.testing.assert_allclose(mesh.edges, [5.0, 6.0, 7.0])

    @pytest.mark.parametrize("coord", _CURVILINEAR, ids=lambda c: c.name.lower())
    def test_a_hollow_body_propagates_its_inner_law(self, coord):
        """The two laws of a hollow body go to (bc_left, bc_right), and the
        mesh starts at the inner radius."""
        g = StructuredGeometry(
            coord=coord,
            breakpoints=(0.5, 2.0),
            mat_ids=(0,),
            boundaries=(BC.reflective, BC.vacuum),
        )
        mesh = Mesh1D.from_geometry(g, region_meshes=(RegionMesh(n_cells=4),))
        assert mesh.edges[0] == 0.5
        assert mesh.edges[-1] == 2.0
        assert mesh.bc_left == BC.reflective
        assert mesh.bc_right == BC.vacuum

    @pytest.mark.catches("ERR-020")
    def test_equal_volume_cylindrical_invariant(self):
        """Equal-volume cells in a cylindrical zone are bit-identical."""
        g = StructuredGeometry(
            coord=CoordSystem.CYLINDRICAL,
            breakpoints=(0.0, 2.0),
            mat_ids=(0,),
            boundaries=(BC("white"),),
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
            coord=CoordSystem.SPHERICAL,
            breakpoints=(0.0, 3.0),
            mat_ids=(0,),
            boundaries=(BC.vacuum,),
        )
        mesh = Mesh1D.from_geometry(g, region_meshes=(RegionMesh(n_cells=12),))
        assert np.all(mesh.volumes == mesh.volumes[0])
        expected_total = (4.0 / 3.0) * np.pi * 3.0 ** 3
        np.testing.assert_allclose(mesh.volumes.sum(), expected_total, rtol=1e-14)

    @pytest.mark.catches("ERR-020")
    @pytest.mark.parametrize(
        "coord", [_SLAB, _CYLINDER, _SPHERE], ids=lambda c: c.name.lower(),
    )
    def test_equal_volume_multi_region_invariant(self, coord):
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
        mesh = _three_region_mesh(coord)
        _assert_equal_volume_regions(
            mesh,
            mat_ids=_THREE_REGION_MAT_IDS,
            cells_per_region=_THREE_REGION_CELLS,
            radii=_THREE_REGION_RADII,
        )

    @pytest.mark.parametrize(
        "coord", [_SLAB, _CYLINDER, _SPHERE], ids=lambda c: c.name.lower(),
    )
    def test_equal_volume_edges_bound_the_volumes(self, coord):
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
        mesh = _three_region_mesh(coord)
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
