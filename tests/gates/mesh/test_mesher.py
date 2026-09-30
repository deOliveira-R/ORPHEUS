r"""The mesher lifts a geometry onto cells (P1 step 3b).

``Mesher(geometry).partition(rule | (rule_0, rule_1, ...)).mesh`` is the one
call site that builds a mesh from a :class:`~orpheus.geometry.StructuredGeometry`
(the user's ruling of 2026-09-29). It applies an interval rule to every
interval (or one rule per interval), concatenates the cells, gives each cell
the material of its interval and each boundary face the geometry's law at
that point, and builds the bare :class:`~orpheus.mesh.Mesh1D`. ``refine(k)``
re-partitions by ``k * rule`` for each rule.

Gate ids (``.claude/plans/reference_p1_spec.md`` §1.3a, "3b owes"): the lift
(S3.8 (b), (d)), S3.1's cross-interval legs and the posing sequence's seed 2
(every breakpoint is a cell edge, every cell lies in exactly one interval),
the mesh-level #495 law (S3.5), and the mesher's own laws (broadcast against
per-interval, refinement nests, the session's order). Every test is
``foundation``; claim kind THEOREM.
"""
from __future__ import annotations

import itertools
import re

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, StructuredGeometry
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

_ALBEDO = BC("albedo", {"albedo": 0.3})


def _geometry(coord: CoordSystem, breakpoints: tuple[float, ...]) -> StructuredGeometry:
    """Distinct materials (7, 3, 5, ...) so a lift that indexes by position
    rather than by interval is seen; distinct laws on the two faces."""
    mat_ids = (7, 3, 5, 1)[: len(breakpoints) - 1]
    if coord is _SLAB:
        return StructuredGeometry.slab(breakpoints, mat_ids, left=BC.reflective, right=_ALBEDO)
    build = StructuredGeometry.cylinder if coord is _CYLINDER else StructuredGeometry.sphere
    if breakpoints[0] > 0.0:
        return build(breakpoints, mat_ids, inner=BC.reflective, outer=_ALBEDO)
    return build(breakpoints, mat_ids, outer=_ALBEDO)


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
    interval's material on its cells and the geometry's laws on the faces.

    Claim kind: THEOREM. Every quantity is decided twice: from the mesh, and
    from the interval rules applied to the geometry directly.
    Mutation witnesses: a lift that repeats the materials by POSITION
    (``np.repeat(mat_ids, counts[::-1])``) reds the material leg; one that
    drops the geometry's laws for a default reds the law leg; one that keeps
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
        assert len(mesh.face_laws) == len(g.boundaries)
        assert all(m is l for m, l in zip(mesh.face_laws, g.boundaries, strict=True))
        assert mesh.boundary_faces == g.boundary_points

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

_G = StructuredGeometry.slab((0.0, 1.0, 2.0), (0, 1), left=BC.vacuum, right=BC.vacuum)

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
