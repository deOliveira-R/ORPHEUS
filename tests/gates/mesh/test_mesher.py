r"""The mesher lifts a geometry onto cells (P1 step 3b).

``Mesher(geometry).partition(rule | (rule_0, rule_1, ...)).mesh`` is the one
call site that builds a mesh from a :class:`~orpheus.geometry.StructuredGeometry`
(the user's ruling of 2026-09-29). It applies an interval rule to every
interval (or one rule per interval), concatenates the cells, gives each cell
the material of its interval and each boundary face the geometry's role at
that point, and builds the bare :class:`~orpheus.mesh.Mesh1D`. ``refine(k)``
re-partitions by ``k * rule`` for each rule.

Gate ids (``.claude/plans/reference_p1_spec.md`` §1.3a, "3b owes"): the lift
(S3.8 (b), (d)), S3.1's cross-interval legs and the posing sequence's seed 2
(every breakpoint is a cell edge, every cell lies in exactly one interval),
the mesh-level #495 law (S3.5), and the mesher's own laws (broadcast against
per-interval, refinement nests, the session's order). Every test is
``foundation``; claim kind THEOREM.

``TestMeshingAGeometry`` holds the meshes of named geometries (a solid sphere,
a slab, a three-region slab, the Wigner-Seitz pin cell, a hollow body, a
shifted origin) and four of the six ERR-020 catchers: within each
equal-volume interval the stored cell volumes are bitwise equal and sum to the
closed-form volume of the interval, and the equal-volume edges bound the
stored volumes. It lived in ``tests/gates/geometry/test_structured_geometry.py``
until the layer-DAG re-filing (`.claude/plans/layer_dag.md`, ruling 21).
"""
from __future__ import annotations

import itertools
import re

import numpy as np
import pytest

from orpheus.geometry import CoordSystem, MirrorDeck, StructuredGeometry, compute_volumes_1d
from orpheus.mesh import (
    CellEdges,
    CellsByCount,
    CellsByMaxWidth,
    EqualVolume,
    EqualWidth,
    Mesh1D,
    Mesher,
)

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/mesh/test_mesher.py"
_PART = "tests/gates/mesh/test_partition.py"
_M1D = "tests/gates/mesh/test_mesh1d.py"
_NESTING = f"{_PART}::TestNesting::test_cells_by_count_nests"
_CONSTRUCTION = f"{_M1D}::TestConstructionLaws::test_refusal"
_REFINEMENT = f"{_PART}::TestRefinement::test_the_fine_edges_contain_the_coarse"

_SLAB = CoordSystem.CARTESIAN
_CYLINDER = CoordSystem.CYLINDRICAL
_SPHERE = CoordSystem.SPHERICAL
_COORDS = (_SLAB, _CYLINDER, _SPHERE)
_CURVILINEAR = (_CYLINDER, _SPHERE)

_TAG = 3  # the outer role: a tag, so that it differs from the inner role and from a deck


def _geometry(coord: CoordSystem, breakpoints: tuple[float, ...]) -> StructuredGeometry:
    """Distinct materials (7, 3, 5, ...) so a lift that indexes by position
    rather than by interval is seen; distinct roles on the two faces."""
    mat_ids = (7, 3, 5, 1)[: len(breakpoints) - 1]
    if coord is _SLAB:
        return StructuredGeometry.slab(breakpoints, mat_ids, left=MirrorDeck(), right=_TAG)
    build = StructuredGeometry.cylinder if coord is _CYLINDER else StructuredGeometry.sphere
    if breakpoints[0] > 0.0:
        return build(breakpoints, mat_ids, inner=1, outer=_TAG)
    return build(breakpoints, mat_ids, outer=_TAG)


_BODIES = {
    "solid": (0.0, 0.5, 1.5, 2.0),
    "hollow": (0.3, 0.5, 1.5, 2.0),
}


def _rules() -> list:
    return [
        pytest.param(CellsByCount.uniform_width(5), id="count-width"),
        pytest.param(CellsByCount.uniform_volume(7), id="count-volume"),
        pytest.param(CellsByMaxWidth(0.13, EqualVolume()), id="max-width-volume"),
        pytest.param(2 * CellsByMaxWidth(0.2, EqualWidth()), id="refined-max-width"),
        pytest.param(
            (CellsByCount.uniform_width(3), CellsByMaxWidth(0.3, EqualVolume()),
             CellEdges(np.array([1.5, 1.6, 2.0]))),
            id="per-interval-mixed",
        ),
    ]


def _rule_for(rule, k: int):
    return rule[k] if isinstance(rule, tuple) else rule


# ─────────────────────────────────────────────────────────────────────
# The lift
# ─────────────────────────────────────────────────────────────────────


class TestTheLift:
    r"""The mesh is the concatenation of each interval's cells, with the
    interval's material on its cells and the geometry's roles on the faces.

    Claim kind: THEOREM. Every quantity is decided twice: from the mesh, and
    from the interval rules applied to the geometry directly.
    Mutation witnesses: a lift that repeats the materials by POSITION
    (``np.repeat(mat_ids, counts[::-1])``) reds the material leg; one that
    drops the geometry's roles for a default reds the role leg; one that keeps
    each shared breakpoint twice is refused by the constructor.
    """

    @pytest.mark.rests_on(_NESTING, _CONSTRUCTION)
    @pytest.mark.parametrize("body", _BODIES, ids=str)
    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    @pytest.mark.parametrize("rule", _rules())
    def test_the_mesh_is_the_intervals_cells(self, rule, coord, body):
        g = _geometry(coord, _BODIES[body])
        mesh = Mesher(g).partition(rule).mesh
        cells = [_rule_for(rule, k).cells(g, iv) for k, iv in enumerate(g.intervals)]
        counts = [len(v) for _, v in cells]
        np.testing.assert_array_equal(
            mesh.edges, np.concatenate([cells[0][0], *(e[1:] for e, _ in cells[1:])]),
        )
        np.testing.assert_array_equal(mesh.volumes, np.concatenate([v for _, v in cells]))
        np.testing.assert_array_equal(mesh.mat_ids, np.repeat(g.mat_ids, counts))
        assert mesh.coord is g.coord
        # the geometry's roles, one per boundary point inner first, land on the
        # faces named by the inventory rule: xmin/xmax, or xmax alone on a
        # solid radial body (written here by hand, not read from face_inventory)
        faces = ("xmin", "xmax") if len(g.boundary_points) == 2 else ("xmax",)
        assert tuple(mesh.face_roles) == faces
        assert all(mesh.face_roles[f] == role for f, role in zip(faces, g.face_roles, strict=True))
        assert mesh.boundary_points == g.boundary_points

    @pytest.mark.rests_on(f"{_HERE}::TestTheLift::test_the_mesh_is_the_intervals_cells")
    @pytest.mark.parametrize("body", _BODIES, ids=str)
    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    @pytest.mark.parametrize("rule", _rules())
    def test_every_cell_lies_in_exactly_one_interval(self, rule, coord, body):
        r"""Posing seed 2 and S3.1's cross-interval leg: every breakpoint is a
        cell edge (bitwise, at the offsets the counts give), and each cell is
        contained in exactly one interval, the one whose material it carries.
        The owner is decided by containment in the breakpoints, independently
        of the mesher's grouping."""
        g = _geometry(coord, _BODIES[body])
        mesh = Mesher(g).partition(rule).mesh
        r = np.asarray(g.breakpoints)
        left, right = mesh.edges[:-1], mesh.edges[1:]
        contained = (r[None, :-1] <= left[:, None]) & (right[:, None] <= r[None, 1:])
        np.testing.assert_array_equal(contained.sum(axis=1), np.ones(mesh.N, dtype=int))
        owner = np.argmax(contained, axis=1)
        np.testing.assert_array_equal(mesh.mat_ids, np.asarray(g.mat_ids)[owner])
        counts = np.bincount(owner, minlength=len(g.intervals))
        offsets = np.concatenate([[0], np.cumsum(counts)])
        np.testing.assert_array_equal(mesh.edges[offsets], r)

    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_one_rule_is_the_rule_on_every_interval(self, coord):
        """Broadcast and per-interval are one spelling of one mesh."""
        g = _geometry(coord, _BODIES["solid"])
        rule = CellsByCount.uniform_volume(6)
        assert Mesher(g).partition(rule).mesh == Mesher(g).partition((rule, rule, rule)).mesh

    def test_the_session(self):
        """``partition`` returns the mesher itself (a session chains); the
        mesher holds its geometry; a new ``partition`` replaces the mesh."""
        g = _geometry(_SLAB, _BODIES["solid"])
        mesher = Mesher(g)
        assert mesher.geometry is g
        assert mesher.partition(CellsByCount.uniform_width(2)) is mesher
        first = mesher.mesh
        second = mesher.partition(CellsByCount.uniform_width(3)).mesh
        assert isinstance(first, Mesh1D) and first.N == 6 and second.N == 9


# ─────────────────────────────────────────────────────────────────────
# Refinement
# ─────────────────────────────────────────────────────────────────────


class TestRefine:
    r"""``refine(k)`` re-partitions by ``k * rule``: every coarse edge is a
    fine edge (bitwise), the count is ``k`` times, and refining twice is the
    product. Mutation witness: a refine that halves a max-width bound gives
    an odd count on ``[0.5, 1.5]`` at ``h = 0.3`` and does not nest."""

    @pytest.mark.rests_on(_REFINEMENT, f"{_HERE}::TestTheLift::test_the_mesh_is_the_intervals_cells")
    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    @pytest.mark.parametrize("rule", _rules()[:4])
    def test_refine_nests(self, rule, coord):
        g = _geometry(coord, _BODIES["hollow"])
        mesher = Mesher(g).partition(rule)
        coarse = mesher.mesh
        fine = mesher.refine(2).mesh
        assert fine.N == 2 * coarse.N
        np.testing.assert_array_equal(fine.edges[::2], coarse.edges)
        assert fine == Mesher(g).partition(2 * rule).mesh
        assert mesher.refine(2).mesh == Mesher(g).partition(4 * rule).mesh

    def test_a_refinement_is_not_a_halved_width(self):
        g = _geometry(_SLAB, (0.5, 1.5))
        mesher = Mesher(g).partition(CellsByMaxWidth(0.3, EqualWidth()))
        assert mesher.mesh.N == 4
        assert mesher.refine(2).mesh.N == 8
        assert Mesher(g).partition(CellsByMaxWidth(0.15, EqualWidth())).mesh.N == 7

    def test_explicit_edges_cannot_be_refined(self):
        g = _geometry(_SLAB, (0.5, 1.5))
        mesher = Mesher(g).partition(CellEdges(np.array([0.5, 0.7, 1.5])))
        with pytest.raises(TypeError, match="no spacing rule"):
            mesher.refine(2)

    def test_the_factor_is_a_power_of_two(self):
        mesher = Mesher(_geometry(_SLAB, (0.5, 1.5))).partition(CellsByCount.uniform_width(3))
        with pytest.raises(ValueError, match="a refinement factor is a power of two"):
            mesher.refine(3)


# ─────────────────────────────────────────────────────────────────────
# The session's refusals
# ─────────────────────────────────────────────────────────────────────

_G = StructuredGeometry.slab((0.0, 1.0, 2.0), (0, 1), left=0, right=0)

_SESSION_REFUSALS = [
    ("not-a-geometry", lambda: Mesher(_G.breakpoints), TypeError,  # type: ignore[arg-type]  # a refusal input
     "Mesher loads a StructuredGeometry"),
    ("mesh-before-partition", lambda: Mesher(_G).mesh, ValueError,
     "Mesher.mesh: partition the geometry first"),
    ("refine-before-partition", lambda: Mesher(_G).refine(2), ValueError,
     "Mesher.refine: partition the geometry first"),
    ("rule-count", lambda: Mesher(_G).partition((CellsByCount.uniform_width(2),)), ValueError,
     "1 rule(s) for 2 interval(s)"),
    ("not-a-rule", lambda: Mesher(_G).partition((CellsByCount.uniform_width(2), 3)),  # type: ignore[arg-type]  # a refusal input
     TypeError, "rule 1 is an interval rule"),
    ("explicit-edges-off-the-interval",
     lambda: Mesher(_G).partition((CellEdges(np.array([0.0, 0.9])), CellsByCount.uniform_width(2))),
     ValueError, "the end edges are the breakpoints"),
    ("rule-misses-its-interval", lambda: Mesher(_G).partition((_ShiftedRule(), CellsByCount.uniform_width(2))),
     ValueError, "every breakpoint must be a cell edge"),
]


class _ShiftedRule:
    """A stub interval rule whose cells miss the interval's right end by one
    ulp: the open protocol admits it, so the mesher must check the span."""

    def cells(self, geometry, interval):
        a, b = interval
        edges = np.array([a, 0.5 * (a + b), float(np.nextafter(b, a))])
        return edges, geometry.measure(edges)

    def __rmul__(self, factor):
        raise TypeError("a stub rule has no refinement")


class _GappedRule:
    """A stub whose cells start one ulp after the interval's left end."""

    def cells(self, geometry, interval):
        a, b = interval
        edges = np.array([float(np.nextafter(a, b)), 0.5 * (a + b), b])
        return edges, geometry.measure(edges)

    def __rmul__(self, factor):
        raise TypeError("a stub rule has no refinement")


class TestEveryBreakpointIsAnEdge:
    r"""The mesher's own check that a rule's cells span their interval bit for
    bit (the elegance review of step 3b, Q3: an interval rule is an open
    protocol, so the built-in rules' pinning cannot stand in for it).

    Claim kind: THEOREM. First red: a stub rule returning shifted edges,
    which before the check built a mesh whose breakpoints were not edges (a
    right shift) or was refused later for another reason."""

    @pytest.mark.parametrize("stub", [_ShiftedRule(), _GappedRule()], ids=["right-end", "left-end"])
    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_a_rule_that_misses_its_interval_is_refused(self, stub, coord):
        g = _geometry(coord, _BODIES["hollow"])
        rules = (CellsByCount.uniform_width(2), stub, CellsByCount.uniform_width(2))
        with pytest.raises(ValueError, match="every breakpoint must be a cell edge"):
            Mesher(g).partition(rules)
        with pytest.raises(ValueError, match="rule 1 "):
            Mesher(g).partition(rules)

    def test_the_stub_is_an_interval_rule(self):
        """The stub passes the protocol door, so the span check is the one
        that refuses it (not the type check)."""
        from orpheus.mesh import IntervalRule

        assert isinstance(_ShiftedRule(), IntervalRule)


class TestSessionRefusals:
    @pytest.mark.parametrize(
        "build, error, fragment", [pytest.param(b, e, f, id=c) for c, b, e, f in _SESSION_REFUSALS],
    )
    def test_refusal(self, build, error, fragment):
        with pytest.raises(error, match=re.escape(fragment)):
            build()

    @pytest.mark.rests_on(f"{_HERE}::TestSessionRefusals::test_refusal")
    def test_refusal_fragments_are_disjoint(self):
        messages = {}
        for case, build, error, _ in _SESSION_REFUSALS:
            with pytest.raises(error) as caught:
                build()
            messages[case] = str(caught.value)
        for case, *_, fragment in _SESSION_REFUSALS:
            for other, message in messages.items():
                if other != case:
                    assert fragment not in message, (case, other, message)


# ─────────────────────────────────────────────────────────────────────
# S3.5 — #495 at the mesh
# ─────────────────────────────────────────────────────────────────────

_495_COUNTS = tuple(range(1, 65)) + (100, 127, 255, 1000)


class TestIssue495Mesh:
    r"""S3.5's mesh leg: on a slab, ``uniform_width(n)`` and
    ``uniform_volume(n)`` give EQUAL meshes (edges, volumes, materials,
    laws); on a cylinder and a sphere they differ at every ``n >= 2``.
    ``n = 5, 7, 9, 11`` are in the population by name (today's
    ``"uniform"`` stored 3, 5, 5 and 3 distinct volumes there on ``[0, 3]``)."""

    @pytest.mark.rests_on(f"{_PART}::TestIssue495::test_equal_width_is_equal_volume_on_a_slab")
    @pytest.mark.parametrize(
        "breakpoints", [(0.0, 3.0), (0.0, 2.872), (1.1, 1.8), (0.0, 0.5, 1.5, 2.0)],
        ids=lambda b: "-".join(str(x) for x in b),
    )
    def test_equal_width_is_equal_volume_on_a_slab(self, breakpoints):
        g = _geometry(_SLAB, breakpoints)
        unequal = [
            n for n in _495_COUNTS
            if Mesher(g).partition(CellsByCount.uniform_width(n)).mesh
            != Mesher(g).partition(CellsByCount.uniform_volume(n)).mesh
        ]
        assert not unequal, f"{len(unequal)} of {len(_495_COUNTS)} counts differ: {unequal[:5]}"
        for n in (5, 7, 9, 11):
            volumes = Mesher(g).partition(CellsByCount.uniform_width(n)).mesh.volumes
            assert len(set(volumes[:n].tolist())) == 1

    @pytest.mark.parametrize("coord", [_CYLINDER, _SPHERE], ids=lambda c: c.name.lower())
    def test_the_rules_differ_on_a_curved_body(self, coord):
        g = _geometry(coord, (0.0, 3.0))
        equal = [
            n for n in _495_COUNTS[1:]
            if Mesher(g).partition(CellsByCount.uniform_width(n)).mesh
            == Mesher(g).partition(CellsByCount.uniform_volume(n)).mesh
        ]
        assert not equal, f"{len(equal)} counts give equal meshes: {equal[:5]}"


def test_the_rules_cover_the_population():
    """The harness's own control: the rule table and the bodies are non-empty
    and every rule builds on every body (a filter over an empty list prints
    the same green as a clean one)."""
    built = 0
    for param, coord, body in itertools.product(_rules(), _COORDS, _BODIES):
        Mesher(_geometry(coord, _BODIES[body])).partition(param.values[0])
        built += 1
    assert built == len(_rules()) * len(_COORDS) * len(_BODIES) == 30


# ─────────────────────────────────────────────────────────────────────
# The mesher lifts a named geometry onto its cells (P1 step 3b). The class
# and test names were kept from the pre-carve node ids at P1 step 3b, so that
# the migration capture compared each mesh bitwise with the one it replaced.
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
        face_roles=(0, 0) if coord is _SLAB else (0,),
    )
    return Mesher(g).partition(tuple(
        CellsByCount.uniform_volume(n) for n in _THREE_REGION_CELLS
    )).mesh


class TestMeshingAGeometry:
    def test_single_region_sphere_equal_volume(self):
        g = StructuredGeometry(
            coord=CoordSystem.SPHERICAL,
            breakpoints=(0.0, 2.0),
            mat_ids=(0,),
            face_roles=(0,),
        )
        mesh = Mesher(g).partition(CellsByCount.uniform_volume(8)).mesh
        assert mesh.N == 8
        assert mesh.coord == CoordSystem.SPHERICAL
        assert mesh.edges[0] == 0.0
        assert mesh.edges[-1] == pytest.approx(2.0)
        # All cells in an equal-volume interval are bit-identical: the
        # stored volume is the equal share m/n (ERR-020's invariant).
        assert np.all(mesh.volumes == mesh.volumes[0])
        # One face role, the outer one: the centre carries none.
        assert dict(mesh.face_roles) == {"xmax": 0}

    def test_single_region_slab_uniform(self):
        g = StructuredGeometry(
            coord=CoordSystem.CARTESIAN,
            breakpoints=(0.0, 4.0),
            mat_ids=(0,),
            face_roles=(0, MirrorDeck()),
        )
        mesh = Mesher(g).partition(CellsByCount.uniform_width(4)).mesh
        assert mesh.N == 4
        assert mesh.coord == CoordSystem.CARTESIAN
        np.testing.assert_allclose(mesh.edges, [0.0, 1.0, 2.0, 3.0, 4.0])
        # A slab has two faces, left then right.
        assert dict(mesh.face_roles) == {"xmin": 0, "xmax": MirrorDeck()}

    def test_multi_region_slab(self):
        g = StructuredGeometry.from_thicknesses(
            coord=CoordSystem.CARTESIAN,
            thicknesses=(0.5, 2.0, 0.5),
            mat_ids=(1, 0, 1),
            face_roles=(0, 0),
        )
        mesh = Mesher(g).partition((
            CellsByCount.uniform_width(2),
            CellsByCount.uniform_width(4),
            CellsByCount.uniform_width(2),
        )).mesh
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
        mesh = Mesher(g).partition((
            CellsByCount.uniform_volume(10),
            CellsByCount.uniform_volume(3),
            CellsByCount.uniform_volume(7),
        )).mesh
        assert mesh.N == 20
        assert mesh.coord == CoordSystem.CYLINDRICAL
        # Outer edge equals r_cell.
        r_cell = 3.6 / np.sqrt(np.pi)
        assert mesh.edges[-1] == pytest.approx(r_cell)
        # mat_id: 10 fuel, 3 clad, 7 cool
        assert (mesh.mat_ids == 2).sum() == 10
        assert (mesh.mat_ids == 1).sum() == 3
        assert (mesh.mat_ids == 0).sum() == 7
        assert dict(mesh.face_roles) == {"xmax": 0}
        _assert_equal_volume_regions(
            mesh,
            mat_ids=(2, 1, 0),
            cells_per_region=(10, 3, 7),
            radii=(0.0, 0.9, 1.1, r_cell),
        )

    def test_length_mismatch_raises(self):
        """RE-POSED from the ``from_geometry`` count mismatch: a tuple of
        rules pairs one to one with the intervals."""
        g = StructuredGeometry.from_thicknesses(
            coord=CoordSystem.CARTESIAN,
            thicknesses=(1.0, 1.0),
            mat_ids=(0, 1),
            face_roles=(0, 0),
        )
        with pytest.raises(ValueError, match=r"1 rule\(s\) for 2 interval\(s\)"):
            Mesher(g).partition((CellsByCount.uniform_volume(4),))

    def test_the_first_breakpoint_is_the_origin(self):
        """RE-POSED from ``test_origin_offset``: ``origin=`` retired, since
        the geometry's first breakpoint states where the mesh starts."""
        g = StructuredGeometry(
            coord=CoordSystem.CARTESIAN,
            breakpoints=(5.0, 7.0),
            mat_ids=(0,),
            face_roles=(0, 0),
        )
        mesh = Mesher(g).partition(CellsByCount.uniform_width(2)).mesh
        np.testing.assert_allclose(mesh.edges, [5.0, 6.0, 7.0])

    @pytest.mark.parametrize("coord", _CURVILINEAR, ids=lambda c: c.name.lower())
    def test_a_hollow_body_propagates_its_inner_role(self, coord):
        """The two roles of a hollow body are its two face roles, inner
        first (both tags: no radial face is a deck), and the mesh starts at the inner radius."""
        g = StructuredGeometry(
            coord=coord,
            breakpoints=(0.5, 2.0),
            mat_ids=(0,),
            face_roles=(1, 0),
        )
        mesh = Mesher(g).partition(CellsByCount.uniform_volume(4)).mesh
        assert mesh.edges[0] == 0.5
        assert mesh.edges[-1] == 2.0
        assert dict(mesh.face_roles) == {"xmin": 1, "xmax": 0}
        assert mesh.boundary_points == (0.5, 2.0)

    @pytest.mark.catches("ERR-020")
    def test_equal_volume_cylindrical_invariant(self):
        """Equal-volume cells in a cylindrical zone are bit-identical."""
        g = StructuredGeometry(
            coord=CoordSystem.CYLINDRICAL,
            breakpoints=(0.0, 2.0),
            mat_ids=(0,),
            face_roles=(0,),
        )
        mesh = Mesher(g).partition(CellsByCount.uniform_volume(10)).mesh
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
            face_roles=(0,),
        )
        mesh = Mesher(g).partition(CellsByCount.uniform_volume(12)).mesh
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
        """The equal-volume edges and the stored equal shares are one mesh.

        The mesher stores the equal-volume cell volumes as the equal
        share of the interval's measure (the ERR-020 fix), so the
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
