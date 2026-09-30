r"""The laws of a partition of a geometry's intervals into cells (P1 step 3a).

A :class:`~orpheus.mesh.Partition` refines a
:class:`~orpheus.geometry.StructuredGeometry`'s cut at its breakpoints
:math:`r_0 < \dots < r_R`: every interval :math:`[r_k, r_{k+1}]` is divided
into cells, and per interval the partition stores the cell edges and the cell
measures (lengths, areas per unit height, volumes). The rules that make one
are :class:`~orpheus.mesh.CellsByCount`, :class:`~orpheus.mesh.CellsByMaxWidth`,
:class:`~orpheus.mesh.CellEdges` and ``2 * rule``, with the spacing rules
:class:`~orpheus.mesh.EqualWidth` and :class:`~orpheus.mesh.EqualVolume`.

The gate ids S3.1 to S3.7 are those of the P1 verification specification
(``.claude/plans/reference_p1_spec.md`` §1.3, placed per step in §1.3a). Every
test here is ``foundation`` (a mathematical or software invariant with no
theory-page label, so no ``verifies``) and states its claim kind:

* THEOREM: a law true for every admissible input, asserted over a stated
  population;
* RECORD: what the code does on a given day, designed to red when a later step
  changes it on purpose.

The measure of an interval :math:`[a, b]` is :math:`m(a, b) = c\,(b^d - a^d)`
with :math:`(c, d) = (1, 1)` on a slab, :math:`(\pi, 2)` on a cylinder and
:math:`(\tfrac43\pi, 3)` on a sphere. The spacing body (R2, the ruling of
2026-09-25) places ``n`` cells at equal steps of :math:`T(r) = r^p`,
:math:`p = 1` for equal width and :math:`p = d` for equal volume, with both
end edges pinned to the breakpoints.

Two rows carry ``catches("ERR-020")``, earned rather than inherited: the
step-3a battery re-dropped ERR-020 into the production rule (every measure
re-derived from the realised edges by ``compute_volumes_1d``, installed
in-process) and both rows reddened on all three coordinates under
``python -O`` ([M] 2026-09-29, spec §1.3a). The existing catchers on
``Mesh1D.from_geometry`` move onto ``Mesh1D(g, CellsByCount.uniform_volume(n))``
in step 3b, when ``from_geometry`` retires.
"""
from __future__ import annotations

import itertools
import math
import re
from collections.abc import Callable

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, StructuredGeometry, compute_volumes_1d
from orpheus.mesh import (
    CellEdges,
    CellsByCount,
    CellsByMaxWidth,
    EqualVolume,
    EqualWidth,
    Partition,
    PartitionRule,
    Refined,
    Spacing,
)

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/mesh/test_partition.py"
_GEOMETRY = "tests/gates/geometry/test_structured_geometry.py"
_NAMED = "tests/gates/geometry/test_named_face_constructors.py"

#: S2.1, the breakpoint laws every partition's geometry rests on.
_S2_1 = f"{_GEOMETRY}::TestBreakpointLaws::test_refusal"
#: The coordinate-generic constructor every fixture geometry here is built by.
_UNIFORM_BOUNDARY = f"{_NAMED}::TestUniformBoundary::test_equals_the_bare_constructor"

_SLAB = CoordSystem.CARTESIAN
_CYLINDER = CoordSystem.CYLINDRICAL
_SPHERE = CoordSystem.SPHERICAL
_COORDS = (_SLAB, _CYLINDER, _SPHERE)
_CURVILINEAR = (_CYLINDER, _SPHERE)
_SPACINGS = (EqualWidth(), EqualVolume())

#: The measure's exponent d and constant c, written from geometry here and
#: never read from ``CoordSystem``, so the tests do not share the SUT's table.
_EXPONENT = {_SLAB: 1, _CYLINDER: 2, _SPHERE: 3}


def _measure(coord: CoordSystem, a: float, b: float) -> float:
    """The closed-form measure of [a, b], spelled in the order the
    equal-volume subdivision has always used (``c * (b**d - a**d)``)."""
    match coord:
        case CoordSystem.CARTESIAN:
            return b - a
        case CoordSystem.CYLINDRICAL:
            return np.pi * (b**2 - a**2)
        case CoordSystem.SPHERICAL:
            return (4.0 / 3.0) * np.pi * (b**3 - a**3)
    raise ValueError(coord)


def _r2_edges(p: int, a: float, b: float, n: int) -> np.ndarray:
    """The ruled R2 body: ``T^-1(T(a) + f_j (T(b) - T(a)))``, ends pinned.

    This is the LAW the spacing rules are ruled to realise (spec §0 item 2),
    written out in the test. It is not an independent reference for the
    formula (a rule and this copy agree by construction if both are
    transcribed from the ruling); its job is the definition pin. The
    independent anchors are the keystone against ``_subdivide_zone`` (the
    verified predecessor), the sum law against the closed-form measure, and
    the refinement law.
    """
    f = np.linspace(0.0, 1.0, n + 1)
    t = a**p + f * (b**p - a**p)
    edges = {1: np.asarray(t, dtype=float), 2: np.sqrt(t), 3: np.cbrt(t)}[p]
    edges[0], edges[-1] = a, b
    return edges


def _power(spacing: Spacing, coord: CoordSystem) -> int:
    """The spacing's power, from the ruling (not read from the SUT)."""
    return 1 if isinstance(spacing, EqualWidth) else _EXPONENT[coord]


def _geometry(coord: CoordSystem, breakpoints: tuple[float, ...]) -> StructuredGeometry:
    """A geometry with one material per interval and vacuum everywhere."""
    return StructuredGeometry.uniform_boundary(
        coord, breakpoints, tuple(range(len(breakpoints) - 1)), BC.vacuum,
    )


def _arrays(*rows: object) -> tuple[np.ndarray, ...]:
    """Per-interval rows as float arrays (the stored type of the fields)."""
    return tuple(np.array(row, dtype=float) for row in rows)


def _ulp(x: float) -> float:
    return float(np.spacing(abs(x)))


def _require(condition: bool, message: str) -> None:
    """A helper's assertion that survives ``python -O``."""
    if not condition:
        raise AssertionError(message)


def _report(failures: list[str], population: int, what: str) -> None:
    """Fail with ``k of N`` and the first failures, or pass; a gate whose
    population is empty is a broken harness, never a pass."""
    _require(population > 0, f"{what}: the population is empty")
    _require(
        not failures,
        f"{what}: {len(failures)} of {population} cases fail; first: "
        + "; ".join(failures[:5]),
    )


# ─────────────────────────────────────────────────────────────────────
# Populations
# ─────────────────────────────────────────────────────────────────────

#: The cell counts of S3.2-S3.4: small, prime, odd, a power of two, and large.
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


#: ``(0.3, 0.7)`` is the witness of the measure's SPELLING: on a sphere
#: ``c * (b**3 - a**3)`` and ``c * b**3 - c * a**3`` differ there and agree on
#: every other interval of this population ([M] 2026-09-29: 0 of 16 on a
#: sphere, 1 of 16 on a cylinder), so without it the bitwise ``fl(m/n)`` pin
#: could not see a re-association of the measure.
_LAW_INTERVALS = _INTERVALS + ((0.3, 0.7),) + _seeded_intervals(8, seed=405)

#: The S3.6 counts: the probe's ``ns`` up to 500.
_REFINE_COUNTS = tuple(range(1, 65)) + (100, 127, 128, 200, 255, 256)

#: A three-interval body whose inner radius changes at every interface, meshed
#: with counts that are not powers of two (ERR-020's multi-region fixture).
_THREE_BREAKPOINTS = (0.0, 0.5, 1.5, 2.0)
_THREE_COUNTS = (5, 7, 11)

_SPACING_IDS = {EqualWidth: "equal-width", EqualVolume: "equal-volume"}


def _coord_spacing_params() -> list:
    return [
        pytest.param(c, s, id=f"{c.name.lower()}-{_SPACING_IDS[type(s)]}")
        for c, s in itertools.product(_COORDS, _SPACINGS)
    ]


# ─────────────────────────────────────────────────────────────────────
# The value: a Partition's own refusals and its equality
# ─────────────────────────────────────────────────────────────────────

#: A two-interval slab [0, 1] | [1, 3] and a valid partition of it.
_G2 = _geometry(_SLAB, (0.0, 1.0, 3.0))
_OK_EDGES = _arrays([0.0, 0.5, 1.0], [1.0, 2.0, 3.0])
_OK_MEASURES = _arrays([0.5, 0.5], [1.0, 1.0])
_ONE_BELOW = float(np.nextafter(1.0, 0.0))
_THREE_BELOW = float(np.nextafter(3.0, 0.0))


#: (case id, edges, measures, error type, the fragment that keys the refusal).
#: Each row is built INSIDE ``pytest.raises`` and checked against ``_G2``, so a
#: clause may live at construction or in ``partition(geometry)``.
_PARTITION_REFUSALS = [
    ("interval-count", ([0.0, 0.5, 1.0],), ([0.5, 0.5],),
     ValueError, "refines the geometry's intervals one for one"),
    ("interior-end-edge-one-ulp-off",
     ([0.0, 0.5, _ONE_BELOW], [_ONE_BELOW, 2.0, 3.0]), _OK_MEASURES,
     ValueError, "end edges must be the breakpoints, bit for bit"),
    ("outer-end-edge-one-ulp-off",
     ([0.0, 0.5, 1.0], [1.0, 2.0, _THREE_BELOW]), _OK_MEASURES,
     ValueError, "end edges must be the breakpoints, bit for bit"),
    ("equal-interior-edges", ([0.0, 1.0, 1.0], [1.0, 2.0, 3.0]), _OK_MEASURES,
     ValueError, "must be strictly increasing"),
    ("decreasing-interior-edges", ([0.0, 0.5, 1.0], [1.0, 3.5, 3.0]), _OK_MEASURES,
     ValueError, "must be strictly increasing"),
    ("zero-measure", _OK_EDGES, ([0.5, 0.0], [1.0, 1.0]),
     ValueError, "a cell measure is positive"),
    ("negative-measure", _OK_EDGES, ([0.5, 0.5], [-1.0, 1.0]),
     ValueError, "a cell measure is positive"),
    ("measure-count", _OK_EDGES, ([0.5], [1.0, 1.0]),
     ValueError, "2 cell(s) need 2 measure(s), got 1"),
    ("measure-array-count", _OK_EDGES, ([0.5, 0.5],),
     ValueError, "one measure array per interval"),
    ("unshared-breakpoint", ([0.0, 0.5, 1.0], [1.5, 2.0, 3.0]), _OK_MEASURES,
     ValueError, "must share their breakpoint"),
    ("no-cell", ([0.0], [1.0, 2.0, 3.0]), ([], [1.0, 1.0]),
     ValueError, "needs at least one cell"),
    ("non-finite-edge", ([0.0, math.nan, 1.0], [1.0, 2.0, 3.0]), _OK_MEASURES,
     ValueError, "must be finite"),
]


class TestPartitionValue:
    r"""S3.1's explicit half: a :class:`Partition` checked against a geometry.

    Claim kind: THEOREM (the defining refusals of a new type, and its
    equality). First red: every row is defining; the input that motivates the
    end-edge row exists in the tree (a breakpoint re-derived by re-adding
    thicknesses is 1 ULP off, spec §0 item 4's census site).
    """

    @pytest.mark.rests_on(_S2_1)
    @pytest.mark.parametrize(
        "edges, measures, error, fragment",
        [pytest.param(e, m, err, f, id=i) for i, e, m, err, f in _PARTITION_REFUSALS],
    )
    def test_refusal(self, edges, measures, error, fragment):
        with pytest.raises(error, match=_literal(fragment)):
            Partition(edges=_arrays(*edges), measures=_arrays(*measures)).partition(_G2)

    @pytest.mark.rests_on(f"{_HERE}::TestPartitionValue::test_refusal")
    def test_refusal_fragments_are_disjoint(self):
        """Each refusal carries its own fragment and no other row's."""
        messages = {}
        for case, edges, measures, error, _ in _PARTITION_REFUSALS:
            with pytest.raises(error) as caught:
                Partition(edges=_arrays(*edges), measures=_arrays(*measures)).partition(_G2)
            messages[case] = str(caught.value)
        for case, *_, fragment in _PARTITION_REFUSALS:
            for other, message in messages.items():
                if fragment == _fragment_of(other):
                    continue
                assert fragment not in message, (
                    f"the fragment of {case!r} ({fragment!r}) also appears in "
                    f"the refusal of {other!r}: {message!r}"
                )

    def test_a_valid_partition_is_its_own_rule(self):
        """The positive leg: a checked partition returns itself."""
        p = Partition(edges=_OK_EDGES, measures=_OK_MEASURES)
        assert p.partition(_G2) is p

    def test_arrays_are_copied_and_read_only(self):
        """A frozen value: the caller's arrays are copied, the stored ones
        refuse writes, and the flat views are stored once."""
        edges = [np.array(e) for e in _OK_EDGES]
        measures = [np.array(m) for m in _OK_MEASURES]
        p = Partition(edges=tuple(edges), measures=tuple(measures))
        edges[0][1] = 0.25
        measures[0][0] = 9.0
        np.testing.assert_array_equal(p.edges[0], [0.0, 0.5, 1.0])
        np.testing.assert_array_equal(p.measures[0], [0.5, 0.5])
        for array in (*p.edges, *p.measures, p.all_edges, p.all_measures):
            assert not array.flags.writeable
        assert p.all_edges is p.all_edges
        assert p.all_measures is p.all_measures

    def test_derived_views(self):
        """``cell_counts``, and the flat edges with each shared breakpoint
        once, and the flat measures, in order."""
        p = Partition(edges=_OK_EDGES, measures=_OK_MEASURES)
        assert p.cell_counts == (2, 2)
        np.testing.assert_array_equal(p.all_edges, [0.0, 0.5, 1.0, 2.0, 3.0])
        np.testing.assert_array_equal(p.all_measures, [0.5, 0.5, 1.0, 1.0])

    @pytest.mark.parametrize(
        "change",
        ["interior-edge", "measure", "interval-structure"],
    )
    def test_equality_is_bitwise_over_edges_and_measures(self, change):
        """Two partitions built from equal values in different objects are
        equal; a one-ULP change of an edge or a measure, or the same flat
        cells grouped into different intervals, makes them unequal."""
        p = Partition(edges=_OK_EDGES, measures=_OK_MEASURES)
        twin = Partition(
            edges=tuple(np.array(e) for e in _OK_EDGES),
            measures=tuple(np.array(m) for m in _OK_MEASURES),
        )
        assert p == twin
        match change:
            case "interior-edge":
                other = Partition(
                    edges=_arrays([0.0, float(np.nextafter(0.5, 1.0)), 1.0], _OK_EDGES[1]),
                    measures=_OK_MEASURES,
                )
            case "measure":
                other = Partition(
                    edges=_OK_EDGES,
                    measures=_arrays(_OK_MEASURES[0], [1.0, float(np.nextafter(1.0, 2.0))]),
                )
            case "interval-structure":
                other = Partition(
                    edges=_arrays([0.0, 0.5, 1.0, 2.0, 3.0]), measures=_arrays([0.5, 0.5, 1.0, 1.0]),
                )
            case _:
                raise ValueError(change)
        assert p != other
        assert other != p

    def test_a_partition_is_not_equal_to_another_type(self):
        p = Partition(edges=_OK_EDGES, measures=_OK_MEASURES)
        assert p != _OK_EDGES
        assert (p == 3) is False

    def test_a_partition_is_unhashable_until_step_5(self):
        """RECORD: content identity (the digest) is P1 step 5's (ruling 4 of
        2026-09-29). When it lands this row reds on purpose: re-pose it as
        the eq/hash contract, never delete it."""
        with pytest.raises(TypeError, match="unhashable"):
            hash(Partition(edges=_OK_EDGES, measures=_OK_MEASURES))

    @pytest.mark.parametrize(
        "rule",
        [
            pytest.param(CellsByCount(2, EqualWidth()), id="cells-by-count"),
            pytest.param(CellsByMaxWidth(0.5, EqualVolume()), id="cells-by-max-width"),
            pytest.param(2 * CellsByMaxWidth(0.5, EqualVolume()), id="refined"),
            pytest.param(CellEdges(edges=_OK_EDGES), id="cell-edges"),
            pytest.param(Partition(edges=_OK_EDGES, measures=_OK_MEASURES), id="partition"),
        ],
    )
    def test_every_rule_is_a_partition_rule(self, rule):
        """Every rule answers ``partition(geometry) -> Partition``."""
        assert isinstance(rule, PartitionRule)
        assert isinstance(rule.partition(_G2), Partition)


def _literal(fragment: str) -> str:
    """A fragment as a regular expression matching it literally."""
    return re.escape(fragment)


def _fragment_of(case: str) -> str:
    return next(f for c, *_, f in _PARTITION_REFUSALS if c == case)


# ─────────────────────────────────────────────────────────────────────
# S3.1 — every rule's partition nests in the geometry
# ─────────────────────────────────────────────────────────────────────


def _check_nesting(p: Partition, g: StructuredGeometry, counts: tuple[int, ...]) -> list[str]:
    """The nesting law of one partition against its geometry, as failures."""
    failures = []
    if p.cell_counts != counts:
        failures.append(f"cell_counts {p.cell_counts} != {counts}")
        return failures
    for k, (e, m) in enumerate(zip(p.edges, p.measures, strict=True)):
        a, b = g.breakpoints[k], g.breakpoints[k + 1]
        if e[0] != a or e[-1] != b:
            failures.append(f"interval {k}: ends ({e[0]!r}, {e[-1]!r}) are not ({a!r}, {b!r})")
        if not np.all(np.diff(e) > 0):
            failures.append(f"interval {k}: edges not strictly increasing")
        if len(e) != counts[k] + 1 or len(m) != counts[k]:
            failures.append(f"interval {k}: {len(e)} edges, {len(m)} measures for {counts[k]} cells")
        if not np.all(m > 0):
            failures.append(f"interval {k}: a non-positive measure")
        if e.flags.writeable or m.flags.writeable:
            failures.append(f"interval {k}: a writeable array")
    return failures


class TestNesting:
    r"""S3.1: a rule's partition refines the geometry's intervals.

    Claim kind: THEOREM. In every interval the end edges ARE the breakpoints,
    bit for bit; the edges strictly increase; the measures are positive; the
    realised count is the rule's. Mutation witness: remove the end pinning
    from the spacing body; the unpinned last edge misses the breakpoint in
    168 cylinder and 21 sphere cases of 21 035 random (interval, n) cases
    per coordinate, and in NONE of this class's population ([M] 2026-09-29,
    the step-3a battery: the population rows stay green), so the pinning's
    catcher is ``test_the_ends_are_pinned_where_the_body_misses_them``,
    which searches for its own witness.
    """

    @pytest.mark.rests_on(_S2_1, _UNIFORM_BOUNDARY, f"{_HERE}::TestPartitionValue::test_refusal")
    @pytest.mark.parametrize("coord, spacing", _coord_spacing_params())
    def test_cells_by_count_nests(self, coord, spacing):
        failures = []
        cases = 0
        for (a, b), n in itertools.product(_LAW_INTERVALS, _COUNTS):
            g = _geometry(coord, (a, b))
            cases += 1
            failures += [
                f"[{a}, {b}] n={n}: {f}"
                for f in _check_nesting(CellsByCount(n, spacing).partition(g), g, (n,))
            ]
        g = _geometry(coord, _THREE_BREAKPOINTS)
        cases += 1
        failures += [
            f"three intervals: {f}"
            for f in _check_nesting(CellsByCount(_THREE_COUNTS, spacing).partition(g), g, _THREE_COUNTS)
        ]
        _report(failures, cases, f"{coord.name} {type(spacing).__name__} nesting")

    @pytest.mark.rests_on(_S2_1, _UNIFORM_BOUNDARY)
    @pytest.mark.parametrize("coord", _CURVILINEAR, ids=lambda c: c.name.lower())
    def test_the_ends_are_pinned_where_the_body_misses_them(self, coord):
        """The discriminating input of the pinning: an interval on which the
        unpinned R2 body lands its last edge off the breakpoint. Found by a
        search in the test itself, so the row states its own witness."""
        p = _EXPONENT[coord]
        rng = np.random.default_rng(11)
        witness = None
        for _ in range(4000):
            a = float(rng.uniform(0.0, 3.0))
            b = a + float(rng.uniform(0.01, 5.0))
            n = int(rng.integers(1, 100))
            f = np.linspace(0.0, 1.0, n + 1)
            t = a**p + f * (b**p - a**p)
            last = float(np.sqrt(t[-1]) if p == 2 else np.cbrt(t[-1]))
            if last != b:
                witness = (a, b, n)
                break
        assert witness is not None, "no interval in the search misses its end: the row is vacuous"
        a, b, n = witness
        partition = CellsByCount(n, EqualVolume()).partition(_geometry(coord, (a, b)))
        assert partition.edges[0][-1] == b


# ─────────────────────────────────────────────────────────────────────
# S3.2 — the measures sum to the interval's measure
# ─────────────────────────────────────────────────────────────────────


class TestSumLaw:
    r"""S3.2: per interval :math:`|\sum m_j - m(a, b)| \le (2 + \lceil\log_2 n\rceil)\,\mathrm{ulp}(m)`.

    Claim kind: THEOREM with a derived tolerance: a stored equal measure
    ``fl(m/n)`` carries half an ULP of ``m/n``, ``n`` of them about one ULP
    of ``m``, and the pairwise sum at most ``ceil(log2 n)`` half-ULPs; the
    shell measures of equal width telescope to the same bound. ``m(a, b)`` is
    the closed form written in the test. Mutation witness: a measure that
    drops the interval's inner radius, ``m(0, b)/n``, reds every interval
    with ``a > 0`` (ERR-020's second leg, generalised).
    """

    @pytest.mark.rests_on(f"{_HERE}::TestNesting::test_cells_by_count_nests")
    @pytest.mark.parametrize("coord, spacing", _coord_spacing_params())
    def test_measures_sum_to_the_interval_measure(self, coord, spacing):
        failures = []
        cases = 0
        for (a, b), n in itertools.product(_LAW_INTERVALS, _COUNTS):
            partition = CellsByCount(n, spacing).partition(_geometry(coord, (a, b)))
            m = _measure(coord, a, b)
            total = float(np.sum(partition.measures[0]))
            bound = (2 + math.ceil(math.log2(n))) * _ulp(m)
            cases += 1
            if abs(total - m) > bound:
                failures.append(
                    f"[{a}, {b}] n={n}: |sum - m| = {abs(total - m) / _ulp(m):.1f} ulp "
                    f"> {bound / _ulp(m):.0f}"
                )
        _report(failures, cases, f"{coord.name} {type(spacing).__name__} sum law")


# ─────────────────────────────────────────────────────────────────────
# S3.3 — equal volume: ERR-020's invariant, on the partition
# ─────────────────────────────────────────────────────────────────────


class TestEqualVolume:
    r"""S3.3: equal-volume measures are bit-identical and each is ``fl(m/n)``.

    Claim kind: THEOREM. Per interval every stored measure is the closed
    form ``m(a, b)`` divided by ``n``, broadcast (ERR-020's invariant: never
    re-derived from the edges, whose ``sqrt``/``cbrt`` round trip loses about
    an ULP per cell). The edges are the R2 body with ``p = d``. Mutation
    witness: measures from ``compute_volumes_1d(coord, edges)`` (ERR-020
    itself) → the equality leg reds on all three coordinates (the three-region
    fixture takes 3/3/2 distinct values per region on a slab, 4/3/3 on a
    cylinder, 5/6/5 on a sphere, spec §1.3 S3.3).
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
            partition = CellsByCount.uniform_volume(n).partition(_geometry(coord, (a, b)))
            expected = _measure(coord, a, b) / n
            cases += 1
            if not np.array_equal(partition.measures[0], np.full(n, expected)):
                distinct = len(set(partition.measures[0].tolist()))
                failures.append(f"[{a}, {b}] n={n}: {distinct} distinct measure(s), not fl(m/n)")
        _report(failures, cases, f"{coord.name} equal-volume measures")

    @pytest.mark.catches("ERR-020")
    @pytest.mark.rests_on(f"{_HERE}::TestNesting::test_cells_by_count_nests")
    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_multi_interval_measures_are_per_interval(self, coord):
        """The three-region fixture: in each interval the measures are
        bit-identical, and they are THAT interval's measure over its count
        (so a measure computed from the origin, ``m(0, b)/n``, reds beyond
        the first interval)."""
        g = _geometry(coord, _THREE_BREAKPOINTS)
        partition = CellsByCount(_THREE_COUNTS, EqualVolume()).partition(g)
        for k, n in enumerate(_THREE_COUNTS):
            a, b = _THREE_BREAKPOINTS[k], _THREE_BREAKPOINTS[k + 1]
            np.testing.assert_array_equal(
                partition.measures[k], np.full(n, _measure(coord, a, b) / n),
                err_msg=f"interval {k} ({coord.name}): not fl(m(a, b)/n)",
            )

    @pytest.mark.rests_on(f"{_HERE}::TestNesting::test_cells_by_count_nests")
    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_edges_are_the_r2_body(self, coord):
        failures = []
        cases = 0
        for (a, b), n in itertools.product(_LAW_INTERVALS, _COUNTS):
            partition = CellsByCount.uniform_volume(n).partition(_geometry(coord, (a, b)))
            cases += 1
            if not np.array_equal(partition.edges[0], _r2_edges(_EXPONENT[coord], a, b, n)):
                failures.append(f"[{a}, {b}] n={n}")
        _report(failures, cases, f"{coord.name} equal-volume edges != R2 body")

    @pytest.mark.rests_on(f"{_HERE}::TestNesting::test_cells_by_count_nests")
    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_keystone_the_predecessor_subdivision(self, coord):
        """The carve keystone (step 3a only): ``_subdivide_zone``, the
        verified equal-volume predecessor that ``Mesh1D.from_geometry`` still
        calls, gives the same interior edges and measures bit for bit; the
        end edges are the breakpoints (the predecessor does not pin them).
        This row retires with ``_subdivide_zone`` at step 3b; the capture of
        spec §3 then carries the bit-identity end to end."""
        from orpheus.mesh.factories import _subdivide_zone

        failures = []
        cases = 0
        for (a, b), n in itertools.product(_LAW_INTERVALS, _COUNTS):
            partition = CellsByCount.uniform_volume(n).partition(_geometry(coord, (a, b)))
            edges, volumes = _subdivide_zone(a, b, n, coord)
            cases += 1
            if not (
                np.array_equal(partition.edges[0][1:-1], edges[1:-1])
                and np.array_equal(partition.measures[0], volumes)
            ):
                failures.append(f"[{a}, {b}] n={n}")
        _report(failures, cases, f"{coord.name} equal volume != _subdivide_zone")


# ─────────────────────────────────────────────────────────────────────
# S3.4 — equal width
# ─────────────────────────────────────────────────────────────────────


class TestEqualWidth:
    r"""S3.4: equal-width edges are the R2 body with ``p = 1``; the measures
    are ``fl((b - a)/n)`` on a slab and the shells between the realised edges
    on a cylinder or a sphere.

    Claim kind: THEOREM. Realised widths are within ``2 ulp(b)`` of the
    nominal ``fl((b - a)/n)`` (bitwise equal widths are unrealisable, spec §0
    item 6). Mutation witness: store ``diff(edges)`` as the slab measure
    (today's ``"uniform"``) → the slab measure leg reds at ``n = 5`` on
    ``[0, 3]`` (3 distinct values, spec §2).
    """

    @pytest.mark.rests_on(f"{_HERE}::TestNesting::test_cells_by_count_nests")
    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_edges_are_the_r2_body_and_widths_are_nominal(self, coord):
        failures = []
        cases = 0
        for (a, b), n in itertools.product(_LAW_INTERVALS, _COUNTS):
            partition = CellsByCount.uniform_width(n).partition(_geometry(coord, (a, b)))
            edges = partition.edges[0]
            cases += 1
            if not np.array_equal(edges, _r2_edges(1, a, b, n)):
                failures.append(f"[{a}, {b}] n={n}: edges != R2 body")
            excess = float(np.max(np.abs(np.diff(edges) - (b - a) / n)))
            if excess > 2 * _ulp(b):
                failures.append(f"[{a}, {b}] n={n}: a width is {excess / _ulp(b):.2f} ulp(b) off")
        _report(failures, cases, f"{coord.name} equal-width edges")

    @pytest.mark.rests_on(f"{_HERE}::TestNesting::test_cells_by_count_nests")
    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_measures(self, coord):
        failures = []
        cases = 0
        for (a, b), n in itertools.product(_LAW_INTERVALS, _COUNTS):
            partition = CellsByCount.uniform_width(n).partition(_geometry(coord, (a, b)))
            if coord is _SLAB:
                expected = np.full(n, (b - a) / n)
            else:
                expected = compute_volumes_1d(coord, partition.edges[0])
            cases += 1
            if not np.array_equal(partition.measures[0], expected):
                failures.append(f"[{a}, {b}] n={n}")
        _report(failures, cases, f"{coord.name} equal-width measures")

    @pytest.mark.rests_on(f"{_HERE}::TestEqualWidth::test_measures")
    def test_the_slab_measure_is_not_the_edge_difference(self):
        """The discriminating input of the slab leg: on ``[0, 3]`` with 5
        cells the R2 edge differences are not all equal (2 distinct values;
        today's ``np.linspace`` edges give 3, spec §2), so a slab measure
        read off the edges is visible here."""
        partition = CellsByCount.uniform_width(5).partition(_geometry(_SLAB, (0.0, 3.0)))
        assert len(set(np.diff(partition.edges[0]).tolist())) > 1
        assert len(set(partition.measures[0].tolist())) == 1


# ─────────────────────────────────────────────────────────────────────
# S3.5 — #495: on a slab, equal width IS equal volume
# ─────────────────────────────────────────────────────────────────────

_SLAB_495_COUNTS = tuple(range(1, 65)) + (100, 127, 255, 1000)
_SLAB_495_BODIES = (
    (0.0, 3.0), (0.0, 2.872), (1.1, 1.8), (0.0, 0.5, 1.5, 2.0),
)


class TestIssue495:
    r"""S3.5, the #495 law: on a slab the two spacing rules are one body.

    Claim kind: THEOREM. ``uniform_width(n)`` and ``uniform_volume(n)`` give
    EQUAL partitions (edges and measures, bitwise) on a slab; on a cylinder
    and a sphere they differ at every ``n >= 2`` (the negative leg, which
    shows the comparison can fail). The dyadic counts are blind (``n = 4``
    and ``16`` agreed even under today's ``"uniform"``), so ``n = 5, 7, 9,
    11`` are named rows. Mutation witness: give ``EqualWidth`` its own body
    ``np.linspace(a, b, n + 1)`` (R1) → the edge leg reds in 576 of 852
    cases (``partition_laws.out`` (b)). The mesh-level leg (equal meshes) is
    step 3b's, once ``Mesh1D(geometry, partition)`` exists.
    """

    @pytest.mark.rests_on(
        f"{_HERE}::TestEqualVolume::test_measures_are_the_interval_measure_over_n",
        f"{_HERE}::TestEqualWidth::test_measures",
    )
    @pytest.mark.parametrize(
        "breakpoints", _SLAB_495_BODIES, ids=lambda b: "-".join(str(x) for x in b),
    )
    def test_equal_width_is_equal_volume_on_a_slab(self, breakpoints):
        g = _geometry(_SLAB, breakpoints)
        failures = [
            f"n={n}"
            for n in _SLAB_495_COUNTS
            if CellsByCount.uniform_width(n).partition(g)
            != CellsByCount.uniform_volume(n).partition(g)
        ]
        _report(failures, len(_SLAB_495_COUNTS), f"slab {breakpoints}: width != volume")

    @pytest.mark.rests_on(f"{_HERE}::TestIssue495::test_equal_width_is_equal_volume_on_a_slab")
    @pytest.mark.parametrize("n", (5, 7, 9, 11))
    def test_the_counts_today_s_uniform_got_wrong(self, n):
        """The first red of #495 on ``[0, 3]``: today's ``"uniform"`` stored
        3, 5, 5 and 3 distinct volumes at n = 5, 7, 9, 11 (spec §2)."""
        g = _geometry(_SLAB, (0.0, 3.0))
        width = CellsByCount.uniform_width(n).partition(g)
        assert width == CellsByCount.uniform_volume(n).partition(g)
        assert len(set(width.measures[0].tolist())) == 1

    @pytest.mark.rests_on(f"{_HERE}::TestIssue495::test_equal_width_is_equal_volume_on_a_slab")
    @pytest.mark.parametrize("coord", _CURVILINEAR, ids=lambda c: c.name.lower())
    def test_the_rules_differ_on_a_curved_body(self, coord):
        """The negative leg: the comparison can fail, and does, wherever the
        measure is not uniform in r."""
        failures = []
        for breakpoints in ((0.0, 3.0), (0.0, 2.872), (1.1, 1.8)):
            g = _geometry(coord, breakpoints)
            for n in _SLAB_495_COUNTS[1:]:
                if (CellsByCount.uniform_width(n).partition(g)
                        == CellsByCount.uniform_volume(n).partition(g)):
                    failures.append(f"{breakpoints} n={n}: equal")
        _report(failures, 3 * (len(_SLAB_495_COUNTS) - 1), f"{coord.name}: rules coincide")

    @pytest.mark.rests_on(
        f"{_HERE}::TestIssue495::test_equal_width_is_equal_volume_on_a_slab",
        f"{_HERE}::TestMaxWidthCount::test_equal_width_count_is_the_least_nominal_count",
    )
    @pytest.mark.parametrize(
        "breakpoints, h, count",
        [
            pytest.param((0.0, 3.0), 0.3, 10, id="0-3-h0.3"),
            pytest.param((0.0, 1.0), 0.1, 10, id="0-1-h0.1"),
            pytest.param((1.1, 1.8), 0.07, None, id="1.1-1.8-h0.07"),
            pytest.param((0.0, 0.5, 1.5, 2.0), 0.13, None, id="three-intervals-h0.13"),
        ],
    )
    def test_max_width_is_one_body_on_a_slab(self, breakpoints, h, count):
        """The one-body law reaches ``CellsByMaxWidth``: on a slab its count
        rule is the same for both spacings. ``(0, 3), h = 0.3`` is the input
        on which the first REALISED cell of the equal-volume body,
        ``fl(fl(1/10) * 3) = 0.30000000000000004``, exceeds ``h`` while the
        nominal ``fl(3/10) = 0.3`` does not ([M] scratchpad probe, 1 of 54 144
        slab pairs), so a count read off the realised first cell gives 11."""
        g = _geometry(_SLAB, breakpoints)
        width = CellsByMaxWidth(h, EqualWidth()).partition(g)
        assert width == CellsByMaxWidth(h, EqualVolume()).partition(g)
        if count is not None:
            assert width.cell_counts == (count,)


# ─────────────────────────────────────────────────────────────────────
# S3.6 — 2 * d refines d
# ─────────────────────────────────────────────────────────────────────


class TestRefinement:
    r"""S3.6: ``(k * d).partition(g)`` refines ``d.partition(g)``.

    Claim kind: THEOREM. Per interval the fine edges at the indices
    ``0, k, 2k, ...`` ARE the coarse edges, bitwise, for ``k`` a power of two
    (``fl(1/(kn)) = fl(1/n)/k``, a power-of-two scaling), and the rule
    realises ``k`` times the coarse counts. A factor that is not a power of
    two is refused: its fine fractions miss the coarse ones by up to 2 ULP
    (``[M]`` 3909 of 7176 cases at k = 3), and cells that do not nest are not
    a refinement (the ruling of 2026-09-29). Mutation witness: spell
    ``2 * CellsByMaxWidth(h)`` as ``CellsByMaxWidth(h / 2)`` → the max-width
    leg reds at ``L = 1, h = 0.3`` (counts 4 and 7).
    """

    @pytest.mark.rests_on(f"{_HERE}::TestNesting::test_cells_by_count_nests")
    @pytest.mark.parametrize("coord, spacing", _coord_spacing_params())
    @pytest.mark.parametrize("factor", (2, 4, 8))
    def test_the_fine_edges_contain_the_coarse(self, coord, spacing, factor):
        failures = []
        cases = 0
        for (a, b), n in itertools.product(_INTERVALS, _REFINE_COUNTS):
            g = _geometry(coord, (a, b))
            coarse = CellsByCount(n, spacing).partition(g)
            fine = (factor * CellsByCount(n, spacing)).partition(g)
            cases += 1
            if fine.cell_counts != (factor * n,):
                failures.append(f"[{a}, {b}] n={n}: counts {fine.cell_counts}")
            elif not np.array_equal(fine.edges[0][::factor], coarse.edges[0]):
                failures.append(f"[{a}, {b}] n={n}: a coarse edge is not a fine edge")
        _report(failures, cases, f"{coord.name} {type(spacing).__name__} x{factor}")

    @pytest.mark.rests_on(f"{_HERE}::TestRefinement::test_the_fine_edges_contain_the_coarse")
    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_equal_measures_are_additive_under_refinement(self, coord):
        """Where the measure is ``fl(m/n)`` (equal volume, and equal width on
        a slab) two fine measures sum to the coarse one bit for bit:
        ``fl(m/(2n)) = fl(m/n)/2``."""
        spacings = _SPACINGS if coord is _SLAB else (EqualVolume(),)
        failures = []
        cases = 0
        for spacing, (a, b), n in itertools.product(spacings, _INTERVALS, _REFINE_COUNTS):
            g = _geometry(coord, (a, b))
            coarse = CellsByCount(n, spacing).partition(g).measures[0]
            fine = (2 * CellsByCount(n, spacing)).partition(g).measures[0]
            cases += 1
            if not np.array_equal(fine[0::2] + fine[1::2], coarse):
                failures.append(f"{type(spacing).__name__} [{a}, {b}] n={n}")
        _report(failures, cases, f"{coord.name}: measures not additive")

    @pytest.mark.rests_on(f"{_HERE}::TestRefinement::test_the_fine_edges_contain_the_coarse")
    @pytest.mark.parametrize("spacing", _SPACINGS, ids=lambda s: _SPACING_IDS[type(s)])
    def test_a_multi_interval_rule_refines_per_interval(self, spacing):
        g = _geometry(_SPHERE, _THREE_BREAKPOINTS)
        coarse = CellsByCount(_THREE_COUNTS, spacing).partition(g)
        fine = (2 * CellsByCount(_THREE_COUNTS, spacing)).partition(g)
        assert fine.cell_counts == tuple(2 * n for n in _THREE_COUNTS)
        for k in range(3):
            np.testing.assert_array_equal(fine.edges[k][::2], coarse.edges[k])

    @pytest.mark.parametrize("spacing", _SPACINGS, ids=lambda s: _SPACING_IDS[type(s)])
    def test_the_refined_rule_is_the_rule_with_doubled_counts(self, spacing):
        assert 2 * CellsByCount(3, spacing) == CellsByCount(6, spacing)
        assert 2 * CellsByCount((3, 5), spacing) == CellsByCount((6, 10), spacing)
        assert 4 * CellsByCount(3, spacing) == 2 * (2 * CellsByCount(3, spacing))
        assert np.int64(2) * CellsByCount(3, spacing) == CellsByCount(6, spacing)
        assert 1 * CellsByCount(3, spacing) == CellsByCount(3, spacing)

    @pytest.mark.rests_on(f"{_HERE}::TestMaxWidthCount::test_equal_width_count_is_the_least_nominal_count")
    @pytest.mark.parametrize("coord, spacing", _coord_spacing_params())
    def test_a_refined_max_width_rule_doubles_its_realised_counts(self, coord, spacing):
        """``2 * CellsByMaxWidth(h, s)`` on ``g`` is ``CellsByCount(2 * <the
        counts it realises on g>, s)``, and it nests; ``2 * (2 * d)`` is the
        factor 4."""
        failures = []
        cases = 0
        for (a, b), h in itertools.product(_INTERVALS, (0.3, 0.07, 1.0)):
            if _real_count(_power(spacing, coord), a, b, h) > 2000:
                continue
            g = _geometry(coord, (a, b))
            d = CellsByMaxWidth(h, spacing)
            coarse = d.partition(g)
            fine = (2 * d).partition(g)
            quadruple = (2 * (2 * d)).partition(g)
            cases += 1
            if fine != CellsByCount(tuple(2 * n for n in d.counts(g)), spacing).partition(g):
                failures.append(f"[{a}, {b}] h={h}: 2 * d is not the doubled count")
            elif not np.array_equal(fine.edges[0][::2], coarse.edges[0]):
                failures.append(f"[{a}, {b}] h={h}: does not nest")
            elif quadruple != (4 * CellsByCount(d.counts(g), spacing)).partition(g):
                failures.append(f"[{a}, {b}] h={h}: 2 * (2 * d) is not the factor 4")
        _report(failures, cases, f"{coord.name} {type(spacing).__name__} refined max width")

    @pytest.mark.rests_on(f"{_HERE}::TestRefinement::test_a_refined_max_width_rule_doubles_its_realised_counts")
    def test_halving_the_width_is_not_a_refinement(self):
        """The refuted spelling, exhibited: on ``[0, 1]`` with ``h = 0.3``
        the rule realises 4 cells; ``h / 2`` realises 7, an odd count that
        cannot nest; ``2 * d`` realises 8."""
        g = _geometry(_SLAB, (0.0, 1.0))
        d = CellsByMaxWidth(0.3, EqualWidth())
        assert d.counts(g) == (4,)
        assert CellsByMaxWidth(0.15, EqualWidth()).counts(g) == (7,)
        assert (2 * d).partition(g).cell_counts == (8,)
        assert isinstance(2 * d, Refined)


#: (case id, factor, error type, fragment). ``factor * rule`` for each rule.
_FACTOR_REFUSALS = [
    ("zero", 0, TypeError, "a refinement factor is a positive int"),
    ("negative", -2, TypeError, "a refinement factor is a positive int"),
    ("float", 2.0, TypeError, "a refinement factor is a positive int"),
    ("bool", True, TypeError, "a refinement factor is a positive int"),
    ("three", 3, ValueError, "a refinement factor is a power of two"),
    ("six", 6, ValueError, "a refinement factor is a power of two"),
]


class TestRefinementRefusals:
    """The factor and the rules that cannot be refined.

    Claim kind: THEOREM (defining refusals). A :class:`Partition` and a
    :class:`CellEdges` carry no spacing rule, so there is nowhere to put the
    new edges.
    """

    @pytest.mark.parametrize(
        "factor, error, fragment",
        [pytest.param(f, e, m, id=i) for i, f, e, m in _FACTOR_REFUSALS],
    )
    @pytest.mark.parametrize(
        "rule",
        [
            pytest.param(CellsByCount(3, EqualWidth()), id="cells-by-count"),
            pytest.param(CellsByMaxWidth(0.3, EqualVolume()), id="cells-by-max-width"),
            pytest.param(2 * CellsByMaxWidth(0.3, EqualVolume()), id="refined"),
        ],
    )
    def test_factor_refusal(self, rule, factor, error, fragment):
        with pytest.raises(error, match=_literal(fragment)):
            _ = factor * rule

    def test_factor_fragments_are_disjoint(self):
        rule = CellsByCount(3, EqualWidth())
        messages = {}
        for case, factor, error, _ in _FACTOR_REFUSALS:
            with pytest.raises(error) as caught:
                _ = factor * rule
            messages[case] = str(caught.value)
        for case, _, _, fragment in _FACTOR_REFUSALS:
            for other, message in messages.items():
                other_fragment = next(f for c, _, _, f in _FACTOR_REFUSALS if c == other)
                if other_fragment == fragment:
                    continue
                assert fragment not in message, (case, other, message)

    @pytest.mark.parametrize(
        "rule",
        [
            pytest.param(Partition(edges=_OK_EDGES, measures=_OK_MEASURES), id="partition"),
            pytest.param(CellEdges(edges=_OK_EDGES), id="cell-edges"),
        ],
    )
    def test_a_rule_with_no_spacing_cannot_be_refined(self, rule):
        with pytest.raises(TypeError, match="no spacing rule"):
            _ = 2 * rule


# ─────────────────────────────────────────────────────────────────────
# S3.7 — CellsByMaxWidth: the least count whose nominal widest cell fits
# ─────────────────────────────────────────────────────────────────────

_MAX_WIDTH_INTERVALS = (
    (0.0, 1.0), (0.0, 3.0), (0.0, 2.872), (0.9, 1.1), (1.1, 1.8), (0.01, 2.0), (2.0, 7.0),
)


def _max_width_bounds() -> tuple[float, ...]:
    rng = np.random.default_rng(7)
    return (0.1, 0.2, 0.3, 0.05, 1 / 3, 0.7, 5.0) + tuple(
        float(h) for h in rng.uniform(0.01, 2.0, 16)
    )


_MAX_WIDTH_BOUNDS = _max_width_bounds()


def _real_count(p: int, a: float, b: float, h: float) -> float:
    """The count whose first cell is exactly ``h`` wide, in real arithmetic:
    ``(b^p - a^p) / ((a + h)^p - a^p)``. It sizes the population: from the
    origin of a sphere it grows as ``(b/h)^3``, so a bound on ``(b - a)/h``
    alone admits partitions of millions of cells."""
    return (b**p - a**p) / ((a + h) ** p - a**p)


class TestMaxWidthCount:
    r"""S3.7: per interval, ``n = min{n : w(n) <= h}`` with ``w`` the rule's
    nominal widest cell.

    Claim kind: THEOREM. For equal width ``w(n) = fl((b - a)/n)`` in every
    coordinate system, written in the test and checked by brute force over
    every smaller count (the least count, never "a count that fits"). For
    equal volume on a cylinder or sphere the nominal widest cell is the
    first cell of the spacing body (the widest: the shells thin outward); the
    law is asserted on the realised partition: its widest cell is within
    ``2 ulp(b)`` of ``h`` from below, and the partition with one cell fewer
    has a first cell wider than ``h``. Discriminating rows: ``(1, fl(1/3))
    -> 3`` and ``(3, fl(1/3)) -> 9`` (the exact-rational count gives 4 and
    10); ``(1, 0.1) -> 10`` (a count checked on realised widths gives 11,
    spec §0 item 7). Mutation witnesses: ``ceil(L/h)`` checked on realised
    widths → the ``(1, 0.1)`` row reds; ``ceil(Fraction(L)/Fraction(h))`` →
    the ``fl(1/3)`` rows red.
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
            (n,) = CellsByMaxWidth(h, EqualWidth()).counts(g)
            smaller = (b - a) / np.arange(1, n, dtype=float)
            cases += 1
            if (b - a) / n > h:
                failures.append(f"[{a}, {b}] h={h}: n={n} does not fit")
            elif np.any(smaller <= h):
                failures.append(f"[{a}, {b}] h={h}: n={n} is not the least")
            widest = float(np.max(np.diff(CellsByMaxWidth(h, EqualWidth()).partition(g).edges[0])))
            if widest > h + 2 * _ulp(b):
                failures.append(f"[{a}, {b}] h={h}: realised width {widest!r}")
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
            (n,) = rule.counts(g)
            edges = rule.partition(g).edges[0]
            cases += 1
            if float(np.max(np.diff(edges))) > h + 2 * _ulp(b):
                failures.append(f"[{a}, {b}] h={h}: n={n}, widest cell exceeds h")
            for m in range(1, n):
                fewer = CellsByCount(m, EqualVolume()).partition(g).edges[0]
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
        assert CellsByMaxWidth(h, EqualWidth()).counts(g) == (count,)
        assert CellsByMaxWidth(h, EqualWidth()).partition(g).cell_counts == (count,)

    def test_a_width_per_interval(self):
        """One bound per interval: each interval takes its own least count."""
        g = _geometry(_SLAB, _THREE_BREAKPOINTS)
        rule = CellsByMaxWidth((0.1, 0.5, 0.25), EqualWidth())
        assert rule.counts(g) == (5, 2, 2)
        assert rule.partition(g).cell_counts == (5, 2, 2)


# ─────────────────────────────────────────────────────────────────────
# The rules' refusals
# ─────────────────────────────────────────────────────────────────────

_G3 = _geometry(_SLAB, _THREE_BREAKPOINTS)

#: (case id, a callable building the rule and partitioning _G3, error, fragment).
_COUNT_REFUSALS: list[tuple[str, Callable[[], object], type[Exception], str]] = [
    ("zero", lambda: CellsByCount(0, EqualWidth()), ValueError, "a cell count is at least 1"),
    ("negative", lambda: CellsByCount(-3, EqualWidth()), ValueError, "a cell count is at least 1"),
    ("zero-in-tuple", lambda: CellsByCount((2, 0, 2), EqualWidth()),
     ValueError, "CellsByCount.counts[1]: a cell count is at least 1"),
    ("float", lambda: CellsByCount(2.0, EqualWidth()),  # type: ignore[arg-type]  # a refusal input
     TypeError, "a cell count is an int"),
    ("bool", lambda: CellsByCount(True, EqualWidth()), TypeError, "a cell count is an int"),
    ("string", lambda: CellsByCount("3", EqualWidth()),  # type: ignore[arg-type]  # a refusal input
     TypeError, "a cell count is an int"),
    ("no-spacing", lambda: CellsByCount(3),  # type: ignore[arg-type]  # a refusal input
     TypeError, "'spacing'"),
    ("string-spacing", lambda: CellsByCount(3, "uniform"),  # type: ignore[arg-type]  # a refusal input
     TypeError, "the spacing is EqualWidth() or EqualVolume()"),
    ("count-per-interval", lambda: CellsByCount((2, 2), EqualWidth()).partition(_G3),
     ValueError, "CellsByCount.counts: 2 entries for 3 interval(s)"),
]

_MAX_WIDTH_REFUSALS: list[tuple[str, Callable[[], object], type[Exception], str]] = [
    ("zero", lambda: CellsByMaxWidth(0.0, EqualWidth()), ValueError, "a width is positive and finite"),
    ("negative", lambda: CellsByMaxWidth(-1.0, EqualWidth()), ValueError, "a width is positive and finite"),
    ("inf", lambda: CellsByMaxWidth(math.inf, EqualWidth()), ValueError, "a width is positive and finite"),
    ("nan", lambda: CellsByMaxWidth(math.nan, EqualWidth()), ValueError, "a width is positive and finite"),
    ("string", lambda: CellsByMaxWidth("0.1", EqualWidth()),  # type: ignore[arg-type]  # a refusal input
     TypeError, "a width is a real number"),
    ("bool", lambda: CellsByMaxWidth(True, EqualWidth()), TypeError, "a width is a real number"),
    ("no-spacing", lambda: CellsByMaxWidth(0.1),  # type: ignore[arg-type]  # a refusal input
     TypeError, "'spacing'"),
    ("string-spacing", lambda: CellsByMaxWidth(0.1, "uniform"),  # type: ignore[arg-type]  # a refusal input
     TypeError, "the spacing is EqualWidth() or EqualVolume()"),
    ("width-per-interval", lambda: CellsByMaxWidth((0.1, 0.2), EqualWidth()).partition(_G3),
     ValueError, "CellsByMaxWidth.widths: 2 entries for 3 interval(s)"),
]


def _refusal_rows(table):
    return [pytest.param(build, err, frag, id=case) for case, build, err, frag in table]


def _assert_disjoint(table) -> None:
    messages = {}
    for case, build, error, _ in table:
        with pytest.raises(error) as caught:
            build()
        messages[case] = str(caught.value)
    for case, _, _, fragment in table:
        for other, message in messages.items():
            other_fragment = next(f for c, _, _, f in table if c == other)
            if other_fragment == fragment or fragment in other_fragment:
                continue
            _require(
                fragment not in message,
                f"the fragment of {case!r} ({fragment!r}) is in {other!r}: {message!r}",
            )


class TestCellsByCountRefusals:
    """The count rule's defining refusals and its spellings.

    Claim kind: THEOREM. There is no default spacing (the ruling of
    2026-09-25), and there is no ``CellsByCount.uniform``.
    """

    @pytest.mark.parametrize("build, error, fragment", _refusal_rows(_COUNT_REFUSALS))
    def test_refusal(self, build, error, fragment):
        with pytest.raises(error, match=_literal(fragment)):
            build()

    @pytest.mark.rests_on(f"{_HERE}::TestCellsByCountRefusals::test_refusal")
    def test_refusal_fragments_are_disjoint(self):
        _assert_disjoint(_COUNT_REFUSALS)

    def test_there_is_no_default_named_uniform(self):
        assert not hasattr(CellsByCount, "uniform")

    def test_the_named_spellings(self):
        assert CellsByCount.uniform_width(4) == CellsByCount(4, EqualWidth())
        assert CellsByCount.uniform_volume(4) == CellsByCount(4, EqualVolume())
        assert CellsByCount.uniform_width(4) != CellsByCount.uniform_volume(4)
        assert CellsByCount(np.int64(4), EqualWidth()) == CellsByCount(4, EqualWidth())  # type: ignore[arg-type]  # a numpy int is admitted


class TestCellsByMaxWidthRefusals:
    """The max-width rule's defining refusals. Claim kind: THEOREM."""

    @pytest.mark.parametrize("build, error, fragment", _refusal_rows(_MAX_WIDTH_REFUSALS))
    def test_refusal(self, build, error, fragment):
        with pytest.raises(error, match=_literal(fragment)):
            build()

    @pytest.mark.rests_on(f"{_HERE}::TestCellsByMaxWidthRefusals::test_refusal")
    def test_refusal_fragments_are_disjoint(self):
        _assert_disjoint(_MAX_WIDTH_REFUSALS)


# ─────────────────────────────────────────────────────────────────────
# CellEdges — the explicit rule
# ─────────────────────────────────────────────────────────────────────


class TestCellEdges:
    r"""The explicit rule: given edges per interval, the measures are the
    geometry's shell measures between them.

    Claim kind: THEOREM. Its nesting refusals are the partition's own (one
    check, reached through ``Partition.partition``).
    """

    @pytest.mark.rests_on(f"{_HERE}::TestPartitionValue::test_refusal")
    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    def test_measures_are_the_shells_between_the_edges(self, coord):
        edges = _arrays([0.0, 0.17, 0.5], [0.5, 0.9, 1.3, 1.5], [1.5, 2.0])
        g = _geometry(coord, _THREE_BREAKPOINTS)
        partition = CellEdges(edges=edges).partition(g)
        for k, e in enumerate(edges):
            np.testing.assert_array_equal(partition.edges[k], e)
            np.testing.assert_array_equal(
                partition.measures[k], compute_volumes_1d(coord, np.asarray(e)),
            )

    @pytest.mark.rests_on(f"{_HERE}::TestEqualWidth::test_measures")
    @pytest.mark.parametrize("coord", _CURVILINEAR, ids=lambda c: c.name.lower())
    def test_equal_width_on_a_curved_body_is_its_own_edges(self, coord):
        """On a cylinder or a sphere the equal-width measures ARE the shells
        between the realised edges, so re-stating the edges gives the same
        partition."""
        g = _geometry(coord, (0.01, 2.0))
        width = CellsByCount.uniform_width(13).partition(g)
        assert CellEdges(edges=width.edges).partition(g) == width

    @pytest.mark.rests_on(f"{_HERE}::TestEqualWidth::test_the_slab_measure_is_not_the_edge_difference")
    def test_equal_width_on_a_slab_is_not_its_own_edges(self):
        """On a slab they differ: the equal-width measure is ``fl(L/n)``, the
        explicit edges' measure is their difference (``[0, 3]``, n = 5)."""
        g = _geometry(_SLAB, (0.0, 3.0))
        width = CellsByCount.uniform_width(5).partition(g)
        explicit = CellEdges(edges=width.edges).partition(g)
        np.testing.assert_array_equal(explicit.edges[0], width.edges[0])
        assert explicit != width

    def test_the_nesting_refusal_is_the_partitions(self):
        with pytest.raises(ValueError, match="end edges must be the breakpoints, bit for bit"):
            CellEdges(edges=_arrays([0.0, 0.5, _ONE_BELOW], [_ONE_BELOW, 2.0, 3.0])).partition(_G2)
        with pytest.raises(ValueError, match="one for one"):
            CellEdges(edges=_arrays([0.0, 0.5, 1.0])).partition(_G2)


# ─────────────────────────────────────────────────────────────────────
# Posing seed 2 — every cell lies in exactly one interval
# ─────────────────────────────────────────────────────────────────────


def _rules_on_three_intervals() -> list:
    return [
        pytest.param(CellsByCount(_THREE_COUNTS, EqualWidth()), id="count-width"),
        pytest.param(CellsByCount(_THREE_COUNTS, EqualVolume()), id="count-volume"),
        pytest.param(CellsByMaxWidth((0.07, 0.3, 0.11), EqualVolume()), id="max-width-volume"),
        pytest.param(2 * CellsByMaxWidth(0.13, EqualWidth()), id="refined-max-width"),
        pytest.param(
            CellEdges(edges=_arrays([0.0, 0.17, 0.5], [0.5, 0.9, 1.3, 1.5], [1.5, 2.0])),
            id="cell-edges",
        ),
    ]


class TestEveryCellInExactlyOneInterval:
    r"""The posing sequence's seed 2, stated on the partition (step 3a): the
    partition refines the material partition, so every flat cell lies in
    exactly one interval, the one whose material it will carry. The
    mesh-level leg (``mesh.mat_ids`` is each interval's material broadcast
    over its cells) is step 3b's.

    Claim kind: THEOREM. The owner of each flat cell is decided twice,
    independently: from the partition's grouping (``cell_counts``), and from
    containment in the geometry's breakpoints. Mutation witness: a flat
    ``all_edges`` that keeps each shared breakpoint twice → a zero-width
    cell lies in two intervals and the containment count reads 2.
    """

    @pytest.mark.rests_on(f"{_HERE}::TestNesting::test_cells_by_count_nests")
    @pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
    @pytest.mark.parametrize("rule", _rules_on_three_intervals())
    def test_every_cell_in_exactly_one_interval(self, coord, rule):
        g = _geometry(coord, _THREE_BREAKPOINTS)
        partition = rule.partition(g)
        flat = partition.all_edges
        n_cells = sum(partition.cell_counts)
        assert len(flat) == n_cells + 1
        assert len(partition.all_measures) == n_cells
        np.testing.assert_array_equal(
            partition.all_measures, np.concatenate(partition.measures),
        )
        left, right = flat[:-1], flat[1:]
        r = np.asarray(g.breakpoints)
        contained = (r[None, :-1] <= left[:, None]) & (right[:, None] <= r[None, 1:])
        np.testing.assert_array_equal(contained.sum(axis=1), np.ones(n_cells, dtype=int))
        owner = np.repeat(np.arange(len(g.mat_ids)), partition.cell_counts)
        np.testing.assert_array_equal(np.argmax(contained, axis=1), owner)
        offsets = np.concatenate([[0], np.cumsum(partition.cell_counts)])
        np.testing.assert_array_equal(flat[offsets], r)


# ─────────────────────────────────────────────────────────────────────
# The spacing rules are values
# ─────────────────────────────────────────────────────────────────────


def test_the_spacing_rules_are_values():
    """Claim kind: THEOREM. Two instances of one spacing rule are equal, the
    two rules are not (the string tags they replace are refused by the
    rules' ``string-spacing`` rows)."""
    assert EqualWidth() == EqualWidth()
    assert EqualVolume() == EqualVolume()
    assert EqualWidth() != EqualVolume()
    assert isinstance(EqualWidth(), Spacing) and isinstance(EqualVolume(), Spacing)
