r"""The laws of the interval rules and of the one measure (P1 step 3a, round 2).

An interval rule divides ONE interval :math:`[a, b]` of a
:class:`~orpheus.geometry.StructuredGeometry` into cells:
``rule.cells(geometry, (a, b)) -> (edges, measures)``. A rule places the
edges; the geometry gives every measure, from the one definition owned by its
coordinate system, :math:`m_j = c\,(T(r_{j+1}) - T(r_j))` with
:math:`T(r) = r^d` evaluated by numpy (``CoordSystem.measure``). The rules are
:class:`~orpheus.mesh.CellsByCount`, :class:`~orpheus.mesh.CellsByMaxWidth`,
:class:`~orpheus.mesh.Refined` (``k * rule``, ``k`` a power of two) and
:class:`~orpheus.mesh.CellEdges`, with the spacing rules
:class:`~orpheus.mesh.EqualWidth` (steps of :math:`r`) and
:class:`~orpheus.mesh.EqualVolume` (steps of the coordinate system's
:math:`T`). Applying rules to a whole geometry is the ``Mesher``'s (step 3b),
so the whole-geometry legs (nesting across intervals, every cell in exactly
one interval, the mesh-level #495 law) are step 3b's, on the mesh.

The gate ids S3.1 to S3.7 are those of the P1 verification specification
(``.claude/plans/reference_p1_spec.md`` §1.3, placed per step in §1.3a). Every
test is ``foundation`` (no theory-page label, so no ``verifies``) and states
its claim kind: THEOREM (a law over a stated population) or RECORD (what the
code does today, designed to red when a later step changes it on purpose).

Bit identity with ``_subdivide_zone`` is NOT a law here (the ruling of
2026-09-29: the measure has one correctly-rounded definition, numpy's power,
where the predecessor used Python's scalar ``**``). ERR-020's invariant is:
equal shares are bit-identical and each is the interval's measure over ``n``,
never a re-derivation from the edges. Two rows carry ``catches("ERR-020")``,
earned by the battery of spec §1.3a (ERR-020 re-dropped in-process, both rows
red on all three coordinates under ``python -O``).
"""
from __future__ import annotations

import itertools
import math
import re
from collections.abc import Callable

import numpy as np
import pytest

import orpheus.geometry.coord as coord_module
from orpheus.geometry import BC, CoordSystem, StructuredGeometry, compute_volumes_1d
from orpheus.geometry.coord import MeasureCoordinate
from orpheus.mesh import (
    CellEdges,
    CellsByCount,
    CellsByMaxWidth,
    CountedRule,
    EqualVolume,
    EqualWidth,
    IntervalRule,
    Refined,
    Spacing,
)

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/mesh/test_partition.py"
_GEOMETRY = "tests/gates/geometry/test_structured_geometry.py"
_NAMED = "tests/gates/geometry/test_named_face_constructors.py"
_S2_1 = f"{_GEOMETRY}::TestBreakpointLaws::test_refusal"
_UNIFORM_BOUNDARY = f"{_NAMED}::TestUniformBoundary::test_equals_the_bare_constructor"
_ONE_MEASURE = f"{_HERE}::TestTheOneMeasure::test_the_measure_is_c_times_the_difference_of_T"

_SLAB = CoordSystem.CARTESIAN
_CYLINDER = CoordSystem.CYLINDRICAL
_SPHERE = CoordSystem.SPHERICAL
_COORDS = (_SLAB, _CYLINDER, _SPHERE)
_CURVILINEAR = (_CYLINDER, _SPHERE)
_SPACINGS = (EqualWidth(), EqualVolume())
_SPACING_IDS = {EqualWidth: "equal-width", EqualVolume: "equal-volume"}

#: The measure's exponent d and constant c, written from geometry here and
#: never read from ``CoordSystem``.
_EXPONENT = {_SLAB: 1, _CYLINDER: 2, _SPHERE: 3}
_CONSTANT = {_SLAB: 1.0, _CYLINDER: np.pi, _SPHERE: (4.0 / 3.0) * np.pi}


def _measure(coord: CoordSystem, edges) -> np.ndarray:
    """The ruled definition, written in the test: ``c * diff(edges**d)`` on arrays."""
    e = np.asarray(edges, dtype=float)
    return _CONSTANT[coord] * np.diff(e ** _EXPONENT[coord])


def _r2_edges(p: int, a: float, b: float, n: int) -> np.ndarray:
    """The ruled R2 body ``T^-1(T(a) + f_j (T(b) - T(a)))``, ends pinned, on arrays.

    The definition pin, not an independent reference (a rule and this copy
    agree if both are transcribed from the ruling); the independent anchors
    are the sum law, the refinement law and the ERR-020 invariant.
    """
    t_a, t_b = np.array([a, b]) ** p
    t = t_a + np.linspace(0.0, 1.0, n + 1) * (t_b - t_a)
    edges = {1: t.copy(), 2: np.sqrt(t), 3: np.cbrt(t)}[p]
    edges[0], edges[-1] = a, b
    return edges


def _power(spacing: Spacing, coord: CoordSystem) -> int:
    """The spacing's exponent, from the ruling (not read from the SUT)."""
    return 1 if isinstance(spacing, EqualWidth) else _EXPONENT[coord]


def _geometry(coord: CoordSystem, breakpoints: tuple[float, ...]) -> StructuredGeometry:
    """A geometry with one material per interval and vacuum everywhere."""
    return StructuredGeometry.uniform_boundary(
        coord, breakpoints, tuple(range(len(breakpoints) - 1)), BC.vacuum,
    )


def _cells(rule: IntervalRule, coord: CoordSystem, interval: tuple[float, float]):
    return rule.cells(_geometry(coord, interval), interval)


def _ulp(x: float) -> float:
    return float(np.spacing(abs(x)))


def _require(condition: bool, message: str) -> None:
    """A helper's assertion that survives ``python -O``."""
    if not condition:
        raise AssertionError(message)


def _report(failures: list[str], population: int, what: str) -> None:
    """Fail with ``k of N`` and the first failures; an empty population is a
    broken harness, never a pass."""
    _require(population > 0, f"{what}: the population is empty")
    _require(
        not failures,
        f"{what}: {len(failures)} of {population} cases fail; first: " + "; ".join(failures[:5]),
    )


def _coord_spacing_params() -> list:
    return [
        pytest.param(c, s, id=f"{c.name.lower()}-{_SPACING_IDS[type(s)]}")
        for c, s in itertools.product(_COORDS, _SPACINGS)
    ]


# ─────────────────────────────────────────────────────────────────────
# Populations
# ─────────────────────────────────────────────────────────────────────

_COUNTS = (1, 2, 3, 5, 8, 13, 64, 100, 1000, 4096)

#: The probe's own interval population (S3.6), hollow bodies included.
_INTERVALS: tuple[tuple[float, float], ...] = (
    (0.0, 1.0), (0.0, 3.0), (0.0, 2.872), (0.9, 1.1), (1.1, 1.8),
    (0.0, 0.41), (0.01, 2.0), (2.0, 7.0),
)


def _seeded_intervals(k: int, seed: int) -> tuple[tuple[float, float], ...]:
    """``k`` intervals, half from the origin, drawn from one fixed seed."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(k):
        a = float(rng.uniform(0.0, 3.0)) if rng.random() < 0.5 else 0.0
        out.append((a, a + float(rng.uniform(0.01, 5.0))))
    return tuple(out)


#: ``(0.3, 0.7)``: on a sphere ``c * (b**3 - a**3)`` and ``c * b**3 - c * a**3``
#: differ there, so a re-association of the measure is visible (spec §1.3a,
#: premise 6).
_LAW_INTERVALS = _INTERVALS + ((0.3, 0.7),) + _seeded_intervals(8, seed=405)

#: The S3.6 counts: the probe's ``ns`` up to 500.
_REFINE_COUNTS = tuple(range(1, 65)) + (100, 127, 128, 200, 255, 256)


# ─────────────────────────────────────────────────────────────────────
# The one measure
# ─────────────────────────────────────────────────────────────────────


class TestTheOneMeasure:
    r"""The measure has one definition, owned by the coordinate system (the
    ruling of 2026-09-29): ``coord.measure(edges) = c * diff(T(edges))``,
    ``T`` numpy's power; the geometry asks it; ``compute_volumes_1d`` is it.

    Claim kind: THEOREM. The route row swaps the definition for a decoy and
    requires every reader to move (a second spelling left anywhere stays
    unmoved and reds).
    """

    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_the_measure_is_c_times_the_difference_of_T(self, coord):
        failures = []
        for (a, b), n in itertools.product(_LAW_INTERVALS, (1, 5, 64)):
            edges = _r2_edges(1, a, b, n)
            if not np.array_equal(coord.measure(edges), _measure(coord, edges)):
                failures.append(f"[{a}, {b}] n={n}")
        _report(failures, 3 * len(_LAW_INTERVALS), f"{coord.name} measure")

    @pytest.mark.rests_on(_ONE_MEASURE)
    def test_every_reader_asks_the_one_definition(self, monkeypatch):
        """A ROUTE gate: with ``CoordSystem.measure`` replaced by a decoy
        (twice the honest measure), ``compute_volumes_1d``, the geometry's
        ``measure`` and both kinds of rule measure (equal shares and realised
        shells) all move by exactly that factor."""
        g = _geometry(_SPHERE, (0.5, 2.0))
        edges = _r2_edges(1, 0.5, 2.0, 7)
        honest = {
            "compute_volumes_1d": compute_volumes_1d(_SPHERE, edges),
            "geometry.measure": g.measure(edges),
            "equal shares": CellsByCount.uniform_volume(7).cells(g, (0.5, 2.0))[1],
            "realised shells": CellsByCount.uniform_width(7).cells(g, (0.5, 2.0))[1],
        }
        original = CoordSystem.measure
        monkeypatch.setattr(CoordSystem, "measure", lambda self, e: 2.0 * original(self, e))
        decoyed = {
            "compute_volumes_1d": compute_volumes_1d(_SPHERE, edges),
            "geometry.measure": g.measure(edges),
            "equal shares": CellsByCount.uniform_volume(7).cells(g, (0.5, 2.0))[1],
            "realised shells": CellsByCount.uniform_width(7).cells(g, (0.5, 2.0))[1],
        }
        for reader, value in honest.items():
            np.testing.assert_array_equal(
                decoyed[reader], 2.0 * value, err_msg=f"{reader} does not read CoordSystem.measure",
            )

    def test_compute_volumes_1d_is_a_delegate(self):
        """RECORD of the retirement's residue: the legacy name survives only
        as a one-line delegate (3b's migration retires it)."""
        import inspect

        body = inspect.getsource(coord_module.compute_volumes_1d)
        assert "coord.measure(edges)" in body
        assert "match" not in body

    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_the_measure_coordinate_is_the_systems(self, coord):
        T = coord.measure_coordinate
        assert T == MeasureCoordinate(_EXPONENT[coord])
        r = np.array([0.0, 0.3, 1.7, 2.0])
        np.testing.assert_array_equal(T(r), r ** _EXPONENT[coord])
        np.testing.assert_allclose(T.inverse(T(r)), r, rtol=4 * np.finfo(float).eps, atol=0.0)
        assert coord.measure_constant == _CONSTANT[coord]

    def test_the_identity_inverse_does_not_alias(self):
        t = np.array([0.0, 1.0])
        assert MeasureCoordinate(1).inverse(t) is not t

    @pytest.mark.parametrize("exponent", [0, 4, -1, 1.5], ids=["0", "4", "-1", "1.5"])
    def test_a_measure_coordinate_is_r_r2_or_r3(self, exponent):
        with pytest.raises(ValueError, match=re.escape("a measure coordinate is r, r**2 or r**3")):
            MeasureCoordinate(exponent)  # type: ignore[arg-type]  # a refusal input

    def test_the_geometry_s_intervals_and_measure(self):
        g = _geometry(_CYLINDER, (0.0, 0.5, 1.5, 2.0))
        assert g.intervals == ((0.0, 0.5), (0.5, 1.5), (1.5, 2.0))
        edges = np.array([0.0, 0.5, 1.5, 2.0])
        np.testing.assert_array_equal(g.measure(edges), _CYLINDER.measure(edges))


# ─────────────────────────────────────────────────────────────────────
# S3.1 — an interval rule's cells nest in their interval
# ─────────────────────────────────────────────────────────────────────


def _nesting_failures(edges, measures, a, b, n) -> list[str]:
    failures = []
    if edges[0] != a or edges[-1] != b:
        failures.append(f"ends ({edges[0]!r}, {edges[-1]!r}) are not ({a!r}, {b!r})")
    if not np.all(np.diff(edges) > 0):
        failures.append("edges not strictly increasing")
    if len(edges) != n + 1 or len(measures) != n:
        failures.append(f"{len(edges)} edges, {len(measures)} measures for {n} cells")
    if not np.all(measures > 0):
        failures.append("a non-positive measure")
    return failures


class TestNesting:
    r"""S3.1, per interval: the end edges ARE the interval's breakpoints, bit
    for bit; the edges strictly increase; the measures are positive; the
    realised count is the rule's.

    Claim kind: THEOREM. The whole-geometry legs are step 3b's, on the mesh.
    Mutation witness: remove the end pinning; the unpinned last edge misses
    the breakpoint in 168 cylinder and 21 sphere cases of 21 035 random
    cases per coordinate and in NONE of this population, so the pinning's
    catcher is the row that searches for its own witness.
    """

    @pytest.mark.rests_on(_S2_1, _UNIFORM_BOUNDARY, _ONE_MEASURE)
    @pytest.mark.parametrize("coord, spacing", _coord_spacing_params())
    def test_cells_by_count_nests(self, coord, spacing):
        failures = []
        cases = 0
        for (a, b), n in itertools.product(_LAW_INTERVALS, _COUNTS):
            edges, measures = _cells(CellsByCount(n, spacing), coord, (a, b))
            cases += 1
            failures += [f"[{a}, {b}] n={n}: {f}" for f in _nesting_failures(edges, measures, a, b, n)]
        _report(failures, cases, f"{coord.name} {type(spacing).__name__} nesting")

    @pytest.mark.rests_on(_S2_1, _UNIFORM_BOUNDARY)
    @pytest.mark.parametrize("coord", _CURVILINEAR, ids=lambda c: c.name.lower())
    def test_the_ends_are_pinned_where_the_body_misses_them(self, coord):
        p = _EXPONENT[coord]
        rng = np.random.default_rng(11)
        witness = None
        for _ in range(4000):
            a = float(rng.uniform(0.0, 3.0))
            b = a + float(rng.uniform(0.01, 5.0))
            n = int(rng.integers(1, 100))
            t_a, t_b = np.array([a, b]) ** p
            t_last = t_a + 1.0 * (t_b - t_a)  # the body's last step, f_n = 1
            last = float(np.sqrt(t_last) if p == 2 else np.cbrt(t_last))
            if last != b:
                witness = (a, b, n)
                break
        assert witness is not None, "no interval in the search misses its end: the row is vacuous"
        a, b, n = witness
        edges, _ = _cells(CellsByCount(n, EqualVolume()), coord, (a, b))
        assert edges[-1] == b

    @pytest.mark.parametrize("rule", [
        pytest.param(CellsByCount(3, EqualWidth()), id="count"),
        pytest.param(CellsByMaxWidth(0.3, EqualVolume()), id="max-width"),
        pytest.param(2 * CellsByCount(3, EqualVolume()), id="refined"),
    ])
    def test_the_returned_arrays_are_fresh(self, rule):
        """Two calls return equal, distinct arrays: a caller's write cannot
        reach the rule or the next call."""
        g = _geometry(_SPHERE, (0.0, 2.0))
        (e1, m1), (e2, m2) = rule.cells(g, (0.0, 2.0)), rule.cells(g, (0.0, 2.0))
        np.testing.assert_array_equal(e1, e2)
        np.testing.assert_array_equal(m1, m2)
        assert e1 is not e2 and m1 is not m2


# ─────────────────────────────────────────────────────────────────────
# S3.2 — the measures sum to the interval's measure
# ─────────────────────────────────────────────────────────────────────


class TestSumLaw:
    r"""S3.2: :math:`|\sum m_j - m(a, b)| \le (2 + \lceil\log_2 n\rceil)\,\mathrm{ulp}(m)`.

    Claim kind: THEOREM with a derived tolerance (spec §1.3 S3.2); ``m(a, b)``
    is the one-cell measure written in the test. Mutation witness: equal
    shares of ``m(0, b)`` (the inner radius dropped) red every ``a > 0``.
    """

    @pytest.mark.rests_on(f"{_HERE}::TestNesting::test_cells_by_count_nests")
    @pytest.mark.parametrize("coord, spacing", _coord_spacing_params())
    def test_measures_sum_to_the_interval_measure(self, coord, spacing):
        failures = []
        cases = 0
        for (a, b), n in itertools.product(_LAW_INTERVALS, _COUNTS):
            _, measures = _cells(CellsByCount(n, spacing), coord, (a, b))
            m = float(_measure(coord, [a, b])[0])
            total = float(np.sum(measures))
            bound = (2 + math.ceil(math.log2(n))) * _ulp(m)
            cases += 1
            if abs(total - m) > bound:
                failures.append(f"[{a}, {b}] n={n}: {abs(total - m) / _ulp(m):.1f} ulp")
        _report(failures, cases, f"{coord.name} {type(spacing).__name__} sum law")


# ─────────────────────────────────────────────────────────────────────
# S3.3 — equal volume: ERR-020's invariant
# ─────────────────────────────────────────────────────────────────────


class TestEqualVolume:
    r"""S3.3: equal-volume measures are bit-identical, and each is the
    interval's measure over ``n``; the edges are the R2 body with ``p = d``.

    Claim kind: THEOREM. Mutation witness: measures re-derived from the
    realised edges (ERR-020 itself) red on all three coordinates.
    """

    @pytest.mark.catches("ERR-020")
    @pytest.mark.rests_on(
        f"{_HERE}::TestNesting::test_cells_by_count_nests",
        f"{_HERE}::TestSumLaw::test_measures_sum_to_the_interval_measure",
    )
    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_measures_are_the_interval_measure_over_n(self, coord):
        failures = []
        cases = 0
        for (a, b), n in itertools.product(_LAW_INTERVALS, _COUNTS):
            _, measures = _cells(CellsByCount.uniform_volume(n), coord, (a, b))
            expected = _measure(coord, [a, b])[0] / n
            cases += 1
            if not np.array_equal(measures, np.full(n, expected)):
                failures.append(f"[{a}, {b}] n={n}: {len(set(measures.tolist()))} distinct")
        _report(failures, cases, f"{coord.name} equal-volume measures")

    @pytest.mark.catches("ERR-020")
    @pytest.mark.rests_on(f"{_HERE}::TestNesting::test_cells_by_count_nests")
    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_each_interval_takes_its_own_share(self, coord):
        """The three-interval body ``(0, .5, 1.5, 2)`` meshed 5 / 7 / 11: in
        each interval the shares are that interval's measure over its count,
        so a share computed from the origin reds beyond the first interval."""
        g = _geometry(coord, (0.0, 0.5, 1.5, 2.0))
        for (a, b), n in zip(g.intervals, (5, 7, 11), strict=True):
            _, measures = CellsByCount(n, EqualVolume()).cells(g, (a, b))
            np.testing.assert_array_equal(
                measures, np.full(n, _measure(coord, [a, b])[0] / n),
                err_msg=f"[{a}, {b}] ({coord.name}): not m(a, b)/n",
            )

    @pytest.mark.rests_on(f"{_HERE}::TestNesting::test_cells_by_count_nests")
    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_edges_are_the_r2_body(self, coord):
        failures = []
        for (a, b), n in itertools.product(_LAW_INTERVALS, _COUNTS):
            edges, _ = _cells(CellsByCount.uniform_volume(n), coord, (a, b))
            if not np.array_equal(edges, _r2_edges(_EXPONENT[coord], a, b, n)):
                failures.append(f"[{a}, {b}] n={n}")
        _report(failures, len(_LAW_INTERVALS) * len(_COUNTS), f"{coord.name} equal-volume edges")


# ─────────────────────────────────────────────────────────────────────
# S3.4 — equal width
# ─────────────────────────────────────────────────────────────────────


class TestEqualWidth:
    r"""S3.4: equal-width edges are the R2 body with ``p = 1``; widths within
    ``2 ulp(b)`` of ``fl((b - a)/n)``; the measures are ``fl((b - a)/n)`` on a
    slab (equal shares) and the geometry's shells between the realised edges
    on a cylinder or a sphere.

    Claim kind: THEOREM. Mutation witness: slab measures from the edges
    (today's ``"uniform"``) → the slab leg reds at ``n = 5`` on ``[0, 3]``.
    """

    @pytest.mark.rests_on(f"{_HERE}::TestNesting::test_cells_by_count_nests")
    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_edges_are_the_r2_body_and_widths_are_nominal(self, coord):
        failures = []
        for (a, b), n in itertools.product(_LAW_INTERVALS, _COUNTS):
            edges, _ = _cells(CellsByCount.uniform_width(n), coord, (a, b))
            if not np.array_equal(edges, _r2_edges(1, a, b, n)):
                failures.append(f"[{a}, {b}] n={n}: edges != R2 body")
            excess = float(np.max(np.abs(np.diff(edges) - (b - a) / n)))
            if excess > 2 * _ulp(b):
                failures.append(f"[{a}, {b}] n={n}: a width is {excess / _ulp(b):.2f} ulp(b) off")
        _report(failures, len(_LAW_INTERVALS) * len(_COUNTS), f"{coord.name} equal-width edges")

    @pytest.mark.rests_on(f"{_HERE}::TestNesting::test_cells_by_count_nests", _ONE_MEASURE)
    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_measures(self, coord):
        failures = []
        for (a, b), n in itertools.product(_LAW_INTERVALS, _COUNTS):
            edges, measures = _cells(CellsByCount.uniform_width(n), coord, (a, b))
            expected = np.full(n, (b - a) / n) if coord is _SLAB else _measure(coord, edges)
            if not np.array_equal(measures, expected):
                failures.append(f"[{a}, {b}] n={n}")
        _report(failures, len(_LAW_INTERVALS) * len(_COUNTS), f"{coord.name} equal-width measures")

    @pytest.mark.rests_on(f"{_HERE}::TestEqualWidth::test_measures")
    def test_the_slab_measure_is_not_the_edge_difference(self):
        """The discriminating input of the slab leg: on ``[0, 3]`` with 5
        cells the R2 edge differences are not all equal (2 distinct values;
        today's ``np.linspace`` edges give 3), the measures are."""
        edges, measures = _cells(CellsByCount.uniform_width(5), _SLAB, (0.0, 3.0))
        assert len(set(np.diff(edges).tolist())) > 1
        assert len(set(measures.tolist())) == 1


# ─────────────────────────────────────────────────────────────────────
# S3.5 — #495: on a slab, equal width IS equal volume
# ─────────────────────────────────────────────────────────────────────

_SLAB_495_COUNTS = tuple(range(1, 65)) + (100, 127, 255, 1000)
_SLAB_495_INTERVALS = ((0.0, 3.0), (0.0, 2.872), (1.1, 1.8), (0.5, 1.5), (-1.0, 2.0))


def _same_cells(x, y) -> bool:
    return all(np.array_equal(p, q) for p, q in zip(x, y, strict=True))


class TestIssue495:
    r"""S3.5, the #495 law per interval: on a slab the two spacing rules are
    one body, so they give the same cells (edges and measures, bitwise); on a
    cylinder and a sphere they differ at every ``n >= 2`` (the negative leg).

    Claim kind: THEOREM. ``n = 5, 7, 9, 11`` are named rows (the dyadic
    counts agreed even under today's ``"uniform"``). Mutation witness: give
    ``EqualWidth`` its own body ``np.linspace(a, b, n + 1)`` (R1). The mesh
    leg (equal meshes) is step 3b's.
    """

    @pytest.mark.rests_on(
        f"{_HERE}::TestEqualVolume::test_measures_are_the_interval_measure_over_n",
        f"{_HERE}::TestEqualWidth::test_measures",
    )
    @pytest.mark.parametrize("interval", _SLAB_495_INTERVALS, ids=lambda i: f"{i[0]}-{i[1]}")
    def test_equal_width_is_equal_volume_on_a_slab(self, interval):
        failures = [
            f"n={n}" for n in _SLAB_495_COUNTS
            if not _same_cells(
                _cells(CellsByCount.uniform_width(n), _SLAB, interval),
                _cells(CellsByCount.uniform_volume(n), _SLAB, interval),
            )
        ]
        _report(failures, len(_SLAB_495_COUNTS), f"slab {interval}: width != volume")

    @pytest.mark.rests_on(f"{_HERE}::TestIssue495::test_equal_width_is_equal_volume_on_a_slab")
    @pytest.mark.parametrize("n", (5, 7, 9, 11))
    def test_the_counts_today_s_uniform_got_wrong(self, n):
        width = _cells(CellsByCount.uniform_width(n), _SLAB, (0.0, 3.0))
        assert _same_cells(width, _cells(CellsByCount.uniform_volume(n), _SLAB, (0.0, 3.0)))
        assert len(set(width[1].tolist())) == 1

    @pytest.mark.rests_on(f"{_HERE}::TestIssue495::test_equal_width_is_equal_volume_on_a_slab")
    @pytest.mark.parametrize("coord", _CURVILINEAR, ids=lambda c: c.name.lower())
    def test_the_rules_differ_on_a_curved_body(self, coord):
        failures = []
        for interval in ((0.0, 3.0), (0.0, 2.872), (1.1, 1.8)):
            for n in _SLAB_495_COUNTS[1:]:
                if _same_cells(
                    _cells(CellsByCount.uniform_width(n), coord, interval),
                    _cells(CellsByCount.uniform_volume(n), coord, interval),
                ):
                    failures.append(f"{interval} n={n}: equal")
        _report(failures, 3 * (len(_SLAB_495_COUNTS) - 1), f"{coord.name}: rules coincide")

    @pytest.mark.rests_on(
        f"{_HERE}::TestIssue495::test_equal_width_is_equal_volume_on_a_slab",
        f"{_HERE}::TestMaxWidthCount::test_equal_width_count_is_the_least_nominal_count",
    )
    @pytest.mark.parametrize(
        "interval, h, count",
        [
            pytest.param((0.0, 3.0), 0.3, 10, id="0-3-h0.3"),
            pytest.param((0.0, 1.0), 0.1, 10, id="0-1-h0.1"),
            pytest.param((1.1, 1.8), 0.07, None, id="1.1-1.8-h0.07"),
        ],
    )
    def test_max_width_is_one_body_on_a_slab(self, interval, h, count):
        """On a slab ``CellsByMaxWidth`` counts the same for both spacings.
        ``(0, 3), h = 0.3``: the realised first cell of the stepped body is
        ``0.30000000000000004 > h`` while the nominal ``fl(3/10) = 0.3`` fits,
        so a count read off the realised first cell gives 11."""
        g = _geometry(_SLAB, interval)
        width, volume = CellsByMaxWidth(h, EqualWidth()), CellsByMaxWidth(h, EqualVolume())
        assert width.count(g, interval) == volume.count(g, interval)
        assert _same_cells(width.cells(g, interval), volume.cells(g, interval))
        if count is not None:
            assert width.count(g, interval) == count


# ─────────────────────────────────────────────────────────────────────
# S3.6 — k * d refines d
# ─────────────────────────────────────────────────────────────────────


class TestRefinement:
    r"""S3.6: ``(k * d).cells`` refines ``d.cells`` for ``k`` a power of two:
    the fine edges at indices ``0, k, 2k, ...`` ARE the coarse edges, bitwise,
    and the count is ``k`` times the coarse count. The factor law lives on
    :class:`Refined` (a type, not the operator), and a refinement of a
    refinement flattens: ``k * (m * r) == Refined(r, k m)``.

    Claim kind: THEOREM. Mutation witnesses: ``2 * CellsByMaxWidth(h)`` spelled
    ``CellsByMaxWidth(h / 2)`` reds the max-width rows; a composition by
    ``k + m`` instead of ``k m`` reds the ``4 * (2 * d)`` row (qa F3: at
    ``2 * (2 * d)`` the sum and the product coincide).
    """

    @pytest.mark.rests_on(f"{_HERE}::TestNesting::test_cells_by_count_nests")
    @pytest.mark.parametrize("coord, spacing", _coord_spacing_params())
    @pytest.mark.parametrize("factor", (2, 4, 8))
    def test_the_fine_edges_contain_the_coarse(self, coord, spacing, factor):
        failures = []
        cases = 0
        for (a, b), n in itertools.product(_INTERVALS, _REFINE_COUNTS):
            g = _geometry(coord, (a, b))
            rule = factor * CellsByCount(n, spacing)
            coarse, _ = CellsByCount(n, spacing).cells(g, (a, b))
            fine, _ = rule.cells(g, (a, b))
            cases += 1
            if rule.count(g, (a, b)) != factor * n or len(fine) != factor * n + 1:
                failures.append(f"[{a}, {b}] n={n}: count {rule.count(g, (a, b))}")
            elif not np.array_equal(fine[::factor], coarse):
                failures.append(f"[{a}, {b}] n={n}: a coarse edge is not a fine edge")
        _report(failures, cases, f"{coord.name} {type(spacing).__name__} x{factor}")

    @pytest.mark.rests_on(f"{_HERE}::TestRefinement::test_the_fine_edges_contain_the_coarse")
    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_equal_shares_are_additive_under_refinement(self, coord):
        """Where the measure is an equal share (equal volume, and equal width
        on a slab), two fine shares sum to the coarse one bit for bit:
        ``fl(m/(2n)) = fl(m/n)/2``."""
        spacings = _SPACINGS if coord is _SLAB else (EqualVolume(),)
        failures = []
        cases = 0
        for spacing, (a, b), n in itertools.product(spacings, _INTERVALS, _REFINE_COUNTS):
            _, coarse = _cells(CellsByCount(n, spacing), coord, (a, b))
            _, fine = _cells(2 * CellsByCount(n, spacing), coord, (a, b))
            cases += 1
            if not np.array_equal(fine[0::2] + fine[1::2], coarse):
                failures.append(f"{type(spacing).__name__} [{a}, {b}] n={n}")
        _report(failures, cases, f"{coord.name}: shares not additive")

    @pytest.mark.parametrize("spacing", _SPACINGS, ids=lambda s: _SPACING_IDS[type(s)])
    def test_the_refinement_type(self, spacing):
        """``k * rule`` is :class:`Refined`; refinements flatten by the
        PRODUCT of the factors; the spacing is the rule's."""
        r = CellsByCount(3, spacing)
        assert 2 * r == Refined(r, 2)
        assert 2 * (4 * r) == Refined(r, 8)
        assert 4 * (2 * r) == Refined(r, 8)
        assert Refined(Refined(r, 2), 4) == Refined(r, 8)
        assert np.int64(2) * r == Refined(r, 2)
        assert (1 * r).count(_geometry(_SLAB, (0.0, 1.0)), (0.0, 1.0)) == 3
        assert (2 * r).spacing == spacing
        assert isinstance(2 * r, CountedRule)

    @pytest.mark.rests_on(f"{_HERE}::TestRefinement::test_the_refinement_type")
    @pytest.mark.parametrize("coord, spacing", _coord_spacing_params())
    def test_a_composed_refinement_multiplies_its_factors(self, coord, spacing):
        """qa F3: ``4 * (2 * d)`` realises 8 times the coarse count and nests
        in it (a factor composed by the sum, 6, would be refused or give 6)."""
        g = _geometry(coord, (0.01, 2.0))
        d = CellsByCount(5, spacing)
        coarse, _ = d.cells(g, (0.01, 2.0))
        fine, _ = (4 * (2 * d)).cells(g, (0.01, 2.0))
        assert (4 * (2 * d)).count(g, (0.01, 2.0)) == 40
        np.testing.assert_array_equal(fine[::8], coarse)

    @pytest.mark.rests_on(f"{_HERE}::TestMaxWidthCount::test_equal_width_count_is_the_least_nominal_count")
    @pytest.mark.parametrize("coord, spacing", _coord_spacing_params())
    def test_a_refined_max_width_rule_doubles_its_realised_count(self, coord, spacing):
        failures = []
        cases = 0
        for (a, b), h in itertools.product(_INTERVALS, (0.3, 0.07, 1.0)):
            if _real_count(_power(spacing, coord), a, b, h) > 2000:
                continue
            g = _geometry(coord, (a, b))
            d = CellsByMaxWidth(h, spacing)
            n = d.count(g, (a, b))
            coarse, _ = d.cells(g, (a, b))
            fine, _ = (2 * d).cells(g, (a, b))
            cases += 1
            if not _same_cells((2 * d).cells(g, (a, b)), CellsByCount(2 * n, spacing).cells(g, (a, b))):
                failures.append(f"[{a}, {b}] h={h}: 2 * d is not the doubled count")
            elif not np.array_equal(fine[::2], coarse):
                failures.append(f"[{a}, {b}] h={h}: does not nest")
        _report(failures, cases, f"{coord.name} {type(spacing).__name__} refined max width")

    @pytest.mark.rests_on(f"{_HERE}::TestRefinement::test_a_refined_max_width_rule_doubles_its_realised_count")
    def test_halving_the_width_is_not_a_refinement(self):
        """On ``[0, 1]`` with ``h = 0.3`` the rule realises 4 cells; ``h / 2``
        realises 7, an odd count that cannot nest; ``2 * d`` realises 8."""
        g = _geometry(_SLAB, (0.0, 1.0))
        d = CellsByMaxWidth(0.3, EqualWidth())
        assert d.count(g, (0.0, 1.0)) == 4
        assert CellsByMaxWidth(0.15, EqualWidth()).count(g, (0.0, 1.0)) == 7
        assert (2 * d).count(g, (0.0, 1.0)) == 8
        assert isinstance(2 * d, Refined)


_FACTOR_REFUSALS = [
    ("zero", 0, ValueError, "a refinement factor is a power of two"),
    ("negative", -2, ValueError, "a refinement factor is a power of two"),
    ("three", 3, ValueError, "a refinement factor is a power of two"),
    ("six", 6, ValueError, "a refinement factor is a power of two"),
    ("float", 2.0, TypeError, "a refinement factor is an int"),
    ("bool", True, TypeError, "a refinement factor is an int"),
]


class TestRefinementRefusals:
    """The factor law, on the type and through the operator; the rules that
    cannot be refined. Claim kind: THEOREM (defining refusals)."""

    @pytest.mark.parametrize(
        "factor, error, fragment", [pytest.param(f, e, m, id=i) for i, f, e, m in _FACTOR_REFUSALS],
    )
    @pytest.mark.parametrize("rule", [
        pytest.param(CellsByCount(3, EqualWidth()), id="cells-by-count"),
        pytest.param(CellsByMaxWidth(0.3, EqualVolume()), id="cells-by-max-width"),
        pytest.param(Refined(CellsByCount(3, EqualWidth()), 2), id="refined"),
    ])
    def test_factor_refusal_through_the_operator(self, rule, factor, error, fragment):
        with pytest.raises(error, match=re.escape(fragment)):
            _ = factor * rule

    @pytest.mark.parametrize(
        "factor, error, fragment", [pytest.param(f, e, m, id=i) for i, f, e, m in _FACTOR_REFUSALS],
    )
    def test_factor_refusal_on_the_type(self, factor, error, fragment):
        """elegance F2: the law is the type's, so the constructor refuses too."""
        with pytest.raises(error, match=re.escape(fragment)):
            Refined(CellsByCount(3, EqualWidth()), factor)

    @pytest.mark.rests_on(f"{_HERE}::TestRefinementRefusals::test_factor_refusal_on_the_type")
    def test_factor_fragments_are_disjoint(self):
        _assert_disjoint([
            (c, (lambda f=f: Refined(CellsByCount(3, EqualWidth()), f)), e, m)
            for c, f, e, m in _FACTOR_REFUSALS
        ])

    def test_only_a_counted_rule_is_refined(self):
        with pytest.raises(TypeError, match=re.escape("Refined.rule is a counted rule")):
            Refined(CellEdges(np.array([0.0, 1.0])), 2)  # type: ignore[arg-type]  # a refusal input
        with pytest.raises(TypeError, match="no spacing rule"):
            _ = 2 * CellEdges(np.array([0.0, 1.0]))


# ─────────────────────────────────────────────────────────────────────
# S3.7 — CellsByMaxWidth: the least count whose nominal widest cell fits
# ─────────────────────────────────────────────────────────────────────

_MAX_WIDTH_INTERVALS = (
    (0.0, 1.0), (0.0, 3.0), (0.0, 2.872), (0.9, 1.1), (1.1, 1.8), (0.01, 2.0), (2.0, 7.0),
)


def _max_width_bounds() -> tuple[float, ...]:
    rng = np.random.default_rng(7)
    return (0.1, 0.2, 0.3, 0.05, 1 / 3, 0.7, 5.0) + tuple(float(h) for h in rng.uniform(0.01, 2.0, 16))


_MAX_WIDTH_BOUNDS = _max_width_bounds()


def _real_count(p: int, a: float, b: float, h: float) -> float:
    """The count whose first cell is exactly ``h`` wide, in real arithmetic:
    ``(b^p - a^p) / ((a + h)^p - a^p)``. It sizes the population: from the
    origin of a sphere it grows as ``(b/h)^3``."""
    return (b**p - a**p) / ((a + h) ** p - a**p)


class TestMaxWidthCount:
    r"""S3.7: ``n = min{n : w(n) <= h}`` with ``w`` the nominal widest cell.

    Claim kind: THEOREM. Equal width: ``w(n) = fl((b - a)/n)`` in every
    coordinate system, checked by brute force over every smaller count. Equal
    volume on a cylinder or sphere: the widest realised cell is within
    ``2 ulp(b)`` of ``h`` from below and no smaller count's first cell fits.
    Discriminating rows: ``(1, fl(1/3)) -> 3``, ``(3, fl(1/3)) -> 9`` (the
    exact-rational count gives 4, 10); ``(1, 0.1) -> 10`` (a count checked on
    realised widths gives 11).
    """

    @pytest.mark.rests_on(f"{_HERE}::TestEqualWidth::test_edges_are_the_r2_body_and_widths_are_nominal")
    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_equal_width_count_is_the_least_nominal_count(self, coord):
        failures = []
        cases = 0
        for (a, b), h in itertools.product(_MAX_WIDTH_INTERVALS, _MAX_WIDTH_BOUNDS):
            if _real_count(1, a, b, h) > 2000:
                continue
            g = _geometry(coord, (a, b))
            rule = CellsByMaxWidth(h, EqualWidth())
            n = rule.count(g, (a, b))
            cases += 1
            if (b - a) / n > h:
                failures.append(f"[{a}, {b}] h={h}: n={n} does not fit")
            elif np.any((b - a) / np.arange(1, n, dtype=float) <= h):
                failures.append(f"[{a}, {b}] h={h}: n={n} is not the least")
            edges, _ = rule.cells(g, (a, b))
            if float(np.max(np.diff(edges))) > h + 2 * _ulp(b):
                failures.append(f"[{a}, {b}] h={h}: a realised width exceeds h + 2 ulp(b)")
        _report(failures, cases, f"{coord.name} equal-width max-width count")

    @pytest.mark.rests_on(f"{_HERE}::TestEqualVolume::test_edges_are_the_r2_body")
    @pytest.mark.parametrize("coord", _CURVILINEAR, ids=lambda c: c.name.lower())
    def test_equal_volume_count_is_the_least_fitting_count(self, coord):
        failures = []
        cases = 0
        for (a, b), h in itertools.product(_MAX_WIDTH_INTERVALS, _MAX_WIDTH_BOUNDS):
            if _real_count(_EXPONENT[coord], a, b, h) > 300:
                continue
            g = _geometry(coord, (a, b))
            rule = CellsByMaxWidth(h, EqualVolume())
            n = rule.count(g, (a, b))
            edges, _ = rule.cells(g, (a, b))
            cases += 1
            if float(np.max(np.diff(edges))) > h + 2 * _ulp(b):
                failures.append(f"[{a}, {b}] h={h}: n={n}, the widest cell exceeds h")
            for m in range(1, n):
                fewer, _ = CellsByCount(m, EqualVolume()).cells(g, (a, b))
                if fewer[1] - fewer[0] <= h - 2 * _ulp(b):
                    failures.append(f"[{a}, {b}] h={h}: n={n} but {m} cells fit")
                    break
        _report(failures, cases, f"{coord.name} equal-volume max-width count")

    @pytest.mark.rests_on(f"{_HERE}::TestMaxWidthCount::test_equal_width_count_is_the_least_nominal_count")
    @pytest.mark.parametrize(
        "length, h, count",
        [
            pytest.param(1.0, 1 / 3, 3, id="1-third"),
            pytest.param(3.0, 1 / 3, 9, id="3-third"),
            pytest.param(1.0, 0.1, 10, id="1-tenth"),
            pytest.param(3.0, 0.3, 10, id="3-h0.3"),
            pytest.param(1.0, 0.3, 4, id="1-h0.3"),
            pytest.param(1.0, 5.0, 1, id="wider-than-the-interval"),
        ],
    )
    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_discriminating_rows(self, coord, length, h, count):
        g = _geometry(coord, (0.0, length))
        rule = CellsByMaxWidth(h, EqualWidth())
        assert rule.count(g, (0.0, length)) == count
        assert len(rule.cells(g, (0.0, length))[1]) == count

    @pytest.mark.parametrize(
        "coord, interval, h, fragment",
        [
            pytest.param(_SLAB, (1000.0, 2000.0), 1e-14, "below the float resolution", id="slab-below-resolution"),
            pytest.param(_SPHERE, (1000.0, 1001.0), 1e-14, "below the float resolution", id="hollow-sphere-below-resolution"),
            pytest.param(_SLAB, (0.0, 1.0), 1e-10, "more than 2147483648", id="slab-too-many-cells"),
            pytest.param(_SPHERE, (0.0, 2.0), 1e-3, "more than 2147483648", id="sphere-too-many-cells"),
            pytest.param(_SPHERE, (1000.0, 1001.0), 1e-13, "more than 2147483648", id="hollow-sphere-too-many-cells"),
        ],
    )
    def test_an_unrealisable_width_is_refused_at_once(self, coord, interval, h, fragment):
        """qa F2: a width below the interval's float resolution raised a raw
        ``ZeroDivisionError``, and one just above it walked about 1e13 counts
        (a hang); both are keyed refusals now, before any walk."""
        g = _geometry(coord, interval)
        # From a sphere's centre only equal volume places (b/h)^3 cells;
        # equal width places b/h = 2000 there and is admitted.
        spacings = (EqualVolume(),) if (coord, interval, h) == (_SPHERE, (0.0, 2.0), 1e-3) else _SPACINGS
        for spacing in spacings:
            with pytest.raises(ValueError, match=re.escape(fragment)):
                CellsByMaxWidth(h, spacing).count(g, interval)


# ─────────────────────────────────────────────────────────────────────
# The rules' refusals and values
# ─────────────────────────────────────────────────────────────────────

_COUNT_REFUSALS: list[tuple[str, Callable[[], object], type[Exception], str]] = [
    ("zero", lambda: CellsByCount(0, EqualWidth()), ValueError, "a cell count is at least 1"),
    ("negative", lambda: CellsByCount(-3, EqualWidth()), ValueError, "a cell count is at least 1"),
    ("float", lambda: CellsByCount(2.0, EqualWidth()),  # type: ignore[arg-type]  # a refusal input
     TypeError, "a cell count is an int"),
    ("bool", lambda: CellsByCount(True, EqualWidth()), TypeError, "a cell count is an int"),
    ("tuple", lambda: CellsByCount((2, 3), EqualWidth()),  # type: ignore[arg-type]  # a refusal input
     TypeError, "a cell count is an int"),
    ("no-spacing", lambda: CellsByCount(3),  # type: ignore[call-arg]  # a refusal input
     TypeError, "'spacing'"),
    ("string-spacing", lambda: CellsByCount(3, "uniform"),  # type: ignore[arg-type]  # a refusal input
     TypeError, "the spacing is EqualWidth() or EqualVolume()"),
]

_MAX_WIDTH_REFUSALS: list[tuple[str, Callable[[], object], type[Exception], str]] = [
    ("zero", lambda: CellsByMaxWidth(0.0, EqualWidth()), ValueError, "a width is positive and finite"),
    ("negative", lambda: CellsByMaxWidth(-1.0, EqualWidth()), ValueError, "a width is positive and finite"),
    ("inf", lambda: CellsByMaxWidth(math.inf, EqualWidth()), ValueError, "a width is positive and finite"),
    ("nan", lambda: CellsByMaxWidth(math.nan, EqualWidth()), ValueError, "width is NaN, which is not a number"),  # NaN: parse_real's own refusal since #405 P1 step 5
    ("string", lambda: CellsByMaxWidth("0.1", EqualWidth()),  # type: ignore[arg-type]  # a refusal input
     TypeError, "must be a real number"),
    ("bool", lambda: CellsByMaxWidth(True, EqualWidth()), TypeError, "must be a real number"),
    ("tuple", lambda: CellsByMaxWidth((0.1, 0.2), EqualWidth()),  # type: ignore[arg-type]  # a refusal input
     TypeError, "must be a real number"),
    ("no-spacing", lambda: CellsByMaxWidth(0.1),  # type: ignore[call-arg]  # a refusal input
     TypeError, "'spacing'"),
    ("string-spacing", lambda: CellsByMaxWidth(0.1, "uniform"),  # type: ignore[arg-type]  # a refusal input
     TypeError, "the spacing is EqualWidth() or EqualVolume()"),
]

_EDGES_REFUSALS: list[tuple[str, Callable[[], object], type[Exception], str]] = [
    ("string-entry", lambda: CellEdges(["0", 1.0]),  # type: ignore[arg-type]  # a refusal input
     TypeError, "must be a real number"),
    ("bool-entry", lambda: CellEdges([False, True]),  # type: ignore[arg-type]  # a refusal input
     TypeError, "must be a real number"),
    ("string", lambda: CellEdges("01"),  # type: ignore[arg-type]  # a refusal input
     TypeError, "must be a sequence of"),
    ("nan", lambda: CellEdges(np.array([0.0, math.nan, 1.0])), ValueError, "edges[1] is NaN, which is not a number"),  # NaN: parse_real's own refusal since #405 P1 step 5
    ("one-edge", lambda: CellEdges(np.array([0.0])), ValueError, "at least two strictly increasing"),
    ("equal-edges", lambda: CellEdges(np.array([0.0, 0.5, 0.5, 1.0])), ValueError, "at least two strictly increasing"),
    ("wrong-ends", lambda: CellEdges(np.array([0.0, 0.5, float(np.nextafter(1.0, 0.0))])).cells(
        _geometry(_SLAB, (0.0, 1.0)), (0.0, 1.0)),
     ValueError, "the end edges are the breakpoints, bit for bit"),
]


def _assert_disjoint(table) -> None:
    """Each row's fragment is in its own message and in no other row's."""
    messages = {}
    for case, build, error, _ in table:
        with pytest.raises(error) as caught:
            build()
        messages[case] = str(caught.value)
    for case, _, _, fragment in table:
        for other, message in messages.items():
            other_fragment = next(f for c, _, _, f in table if c == other)
            if other_fragment == fragment:
                continue
            _require(fragment not in message, f"the fragment of {case!r} ({fragment!r}) is in {other!r}: {message!r}")


def _refusal_rows(table):
    return [pytest.param(build, err, frag, id=case) for case, build, err, frag in table]


class TestRuleRefusals:
    """The rules' defining refusals (claim kind: THEOREM). There is no
    default spacing, no ``CellsByCount.uniform``, and a count or width is one
    scalar per rule (the per-interval choice is the Mesher's, step 3b)."""

    @pytest.mark.parametrize("build, error, fragment", _refusal_rows(_COUNT_REFUSALS))
    def test_cells_by_count(self, build, error, fragment):
        with pytest.raises(error, match=re.escape(fragment)):
            build()

    @pytest.mark.parametrize("build, error, fragment", _refusal_rows(_MAX_WIDTH_REFUSALS))
    def test_cells_by_max_width(self, build, error, fragment):
        with pytest.raises(error, match=re.escape(fragment)):
            build()

    @pytest.mark.parametrize("build, error, fragment", _refusal_rows(_EDGES_REFUSALS))
    def test_cell_edges(self, build, error, fragment):
        """qa F5: strings and bools are refused by the shared scalar parser."""
        with pytest.raises(error, match=re.escape(fragment)):
            build()

    @pytest.mark.parametrize("table", [_COUNT_REFUSALS, _MAX_WIDTH_REFUSALS, _EDGES_REFUSALS],
                             ids=["cells-by-count", "cells-by-max-width", "cell-edges"])
    def test_fragments_are_disjoint(self, table):
        _assert_disjoint(table)

    def test_there_is_no_default_named_uniform(self):
        assert not hasattr(CellsByCount, "uniform")

    def test_the_named_spellings(self):
        assert CellsByCount.uniform_width(4) == CellsByCount(4, EqualWidth())
        assert CellsByCount.uniform_volume(4) == CellsByCount(4, EqualVolume())
        assert CellsByCount.uniform_width(4) != CellsByCount.uniform_volume(4)
        assert CellsByCount(np.int64(4), EqualWidth()) == CellsByCount(4, EqualWidth())  # type: ignore[arg-type]  # a numpy int is admitted


class TestCellEdges:
    """The explicit rule: one interval's edges, the geometry's measures.

    Claim kind: THEOREM.
    """

    @pytest.mark.rests_on(_ONE_MEASURE)
    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_measures_are_the_geometry_s(self, coord):
        edges = [0.5, 0.9, 1.3, 1.5]
        cells_edges, measures = _cells(CellEdges(np.array(edges)), coord, (0.5, 1.5))
        np.testing.assert_array_equal(cells_edges, edges)
        np.testing.assert_array_equal(measures, _measure(coord, edges))

    @pytest.mark.rests_on(f"{_HERE}::TestEqualWidth::test_measures")
    @pytest.mark.parametrize("coord", _CURVILINEAR, ids=lambda c: c.name.lower())
    def test_equal_width_on_a_curved_body_is_its_own_edges(self, coord):
        width = _cells(CellsByCount.uniform_width(13), coord, (0.01, 2.0))
        assert _same_cells(_cells(CellEdges(width[0]), coord, (0.01, 2.0)), width)

    @pytest.mark.rests_on(f"{_HERE}::TestEqualWidth::test_the_slab_measure_is_not_the_edge_difference")
    def test_equal_width_on_a_slab_is_not_its_own_edges(self):
        width = _cells(CellsByCount.uniform_width(5), _SLAB, (0.0, 3.0))
        explicit = _cells(CellEdges(width[0]), _SLAB, (0.0, 3.0))
        np.testing.assert_array_equal(explicit[0], width[0])
        assert not np.array_equal(explicit[1], width[1])

    def test_value_equality_and_hash(self):
        """Equal edges from different objects are one value with one hash; a
        one-ULP change of an edge is another value."""
        x = CellEdges(np.array([0.0, 0.5, 1.0]))
        y = CellEdges(np.array([0.0, 0.5, 1.0]))
        z = CellEdges(np.array([0.0, float(np.nextafter(0.5, 1.0)), 1.0]))
        assert x == y and hash(x) == hash(y)
        assert x != z
        assert len({x, y, z}) == 2
        assert (x == [0.0, 0.5, 1.0]) is False

    def test_negative_zero_is_canonicalised(self):
        """qa F4: ``-0.0`` and ``+0.0`` are one position, so one value with
        one hash, and the stored bit is ``+0.0``."""
        x = CellEdges(np.array([-0.0, 1.0]))
        assert not np.signbit(x.edges[0])
        assert x == CellEdges(np.array([0.0, 1.0])) and hash(x) == hash(CellEdges(np.array([0.0, 1.0])))

    def test_the_edges_are_copied_and_read_only(self):
        source = np.array([0.0, 0.5, 1.0])
        x = CellEdges(source)
        source[1] = 0.25
        assert x.edges[1] == 0.5
        assert not x.edges.flags.writeable
        returned, _ = x.cells(_geometry(_SLAB, (0.0, 1.0)), (0.0, 1.0))
        assert returned is not x.edges


@pytest.mark.parametrize("rule", [
    pytest.param(CellsByCount(2, EqualWidth()), id="cells-by-count"),
    pytest.param(CellsByMaxWidth(0.5, EqualVolume()), id="cells-by-max-width"),
    pytest.param(2 * CellsByMaxWidth(0.5, EqualVolume()), id="refined"),
    pytest.param(CellEdges(np.array([0.0, 0.4, 1.0])), id="cell-edges"),
])
def test_every_rule_is_an_interval_rule(rule):
    """Claim kind: THEOREM. Every rule answers ``cells(geometry, interval)``."""
    assert isinstance(rule, IntervalRule)
    edges, measures = rule.cells(_geometry(_SPHERE, (0.0, 1.0)), (0.0, 1.0))
    assert len(edges) == len(measures) + 1


def test_the_spacing_rules_are_values():
    """Claim kind: THEOREM."""
    assert EqualWidth() == EqualWidth()
    assert EqualVolume() == EqualVolume()
    assert EqualWidth() != EqualVolume()
    assert isinstance(EqualWidth(), Spacing) and isinstance(EqualVolume(), Spacing)
    assert EqualWidth().coordinate(_SPHERE) == MeasureCoordinate(1)
    assert EqualVolume().coordinate(_SPHERE) == MeasureCoordinate(3)
