r"""Partitions of a 1-D geometry into cells, and the rules that make them.

A :class:`~orpheus.geometry.StructuredGeometry` cuts its interval of
positions at its breakpoints :math:`r_0 < r_1 < \dots < r_R` into material
intervals. A **partition** refines that cut: it divides every interval
:math:`[r_k, r_{k+1}]` into cells, so every cell lies in exactly one
interval and every breakpoint is a cell edge. A :class:`Partition` stores,
per interval, the cell edges and the cell measures (lengths on a slab,
areas per unit height on a cylinder, volumes on a sphere).

The measures are stored, not re-derived from the edges. An equal-volume
cell's measure is the interval's measure divided by the cell count,
:math:`m(r_k, r_{k+1}) / n`, while the same volume recomputed from the
realised edges loses about one unit in the last place through the
``sqrt``/``cbrt`` round trip; that loss is ERR-020, and storing the exact
measure is its fix.

A **partition rule** is anything with ``partition(geometry) -> Partition``:

* :class:`CellsByCount`: a number of cells per interval, placed by a
  spacing rule;
* :class:`CellsByMaxWidth`: the fewest cells per interval whose widest
  cell is no wider than a bound, placed by a spacing rule;
* :class:`CellEdges`: the cell edges written out, per interval;
* a :class:`Partition` itself, which returns itself once it is checked
  against the geometry.

``2 * rule`` is the rule's refinement: twice the cells the rule realises
in every interval, placed by the same spacing, so every coarse edge is a
fine edge.

The spacing rules
-----------------

A spacing rule places :math:`n` cells in :math:`[a, b]` at equal steps of
the power :math:`T(r) = r^p`:

.. math::

    r_j = T^{-1}\!\left(T(a) + f_j\,\bigl(T(b) - T(a)\bigr)\right),
    \qquad f_j = j / n,

with both end edges pinned to the breakpoints. :class:`EqualWidth` spaces
in :math:`p = 1` in every coordinate system; :class:`EqualVolume` spaces in
the coordinate system's own exponent :math:`d` (1, 2, 3), in which the
measure :math:`m(a, b) = c\,(b^d - a^d)` is uniform
(:class:`~orpheus.geometry.CoordSystem`). When the spacing's power equals
the measure's exponent the cells have equal measure, stored as
:math:`m(a, b)/n`; otherwise each measure is the shell between the
realised edges. On a slab both rules space in :math:`p = 1 = d`, so they
are one body and give one partition (#495, fixed at its root).
"""
from __future__ import annotations

import itertools
import math
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from numbers import Real
from typing import TYPE_CHECKING, Protocol, runtime_checkable

import numpy as np

from orpheus.geometry.coord import CoordSystem, compute_volumes_1d

if TYPE_CHECKING:
    from orpheus.geometry.structured_geometry import StructuredGeometry


# ═══════════════════════════════════════════════════════════════════════
# The spacing rules
# ═══════════════════════════════════════════════════════════════════════


class Spacing(ABC):
    r"""Where a spacing rule places the interior edges of an interval's cells.

    A spacing is the power :math:`p` in which the cells are equal steps
    (the module docstring).
    """

    @abstractmethod
    def power(self, coord: CoordSystem) -> int:
        r"""The power :math:`p` of :math:`T(r) = r^p` the cells are equal steps of."""

    def cells(
        self, coord: CoordSystem, a: float, b: float, n: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""The ``n + 1`` edges and ``n`` measures of ``n`` cells in :math:`[a, b]`."""
        p = self.power(coord)
        f = np.linspace(0.0, 1.0, n + 1)
        edges = _root(p, a**p + f * (b**p - a**p))
        edges[0], edges[-1] = a, b
        if p == coord.measure_exponent:
            measures = np.full(n, coord.interval_measure(a, b) / n)
        else:
            measures = compute_volumes_1d(coord, edges)
        return edges, measures

    def widest_cell(self, coord: CoordSystem, a: float, b: float, n: int) -> float:
        r"""The nominal width of the widest of ``n`` cells in :math:`[a, b]`.

        For :math:`p = 1` it is the nominal width :math:`\mathrm{fl}((b - a)/n)`;
        for :math:`p > 1` the cells narrow outward, and it is the first
        cell's width from the spacing body.
        """
        p = self.power(coord)
        if p == 1:
            return (b - a) / n
        return float(_root(p, a**p + (1.0 / n) * (b**p - a**p))) - a


@dataclass(frozen=True)
class EqualWidth(Spacing):
    r"""Cells of equal width: equal steps of :math:`r` in every coordinate system."""

    def power(self, coord: CoordSystem) -> int:
        return 1


@dataclass(frozen=True)
class EqualVolume(Spacing):
    r"""Cells of equal measure: equal steps of :math:`r^d`, :math:`d` the coordinate system's exponent."""

    def power(self, coord: CoordSystem) -> int:
        return coord.measure_exponent


def _root(p: int, t: "np.ndarray | float") -> np.ndarray:
    r""":math:`T^{-1}(t) = t^{1/p}` for :math:`p = 1, 2, 3`."""
    match p:
        case 1:
            return np.asarray(t, dtype=float)
        case 2:
            return np.sqrt(t)
        case 3:
            return np.cbrt(t)
    raise ValueError(f"a spacing power is 1, 2 or 3; got {p}")


# ═══════════════════════════════════════════════════════════════════════
# The partition
# ═══════════════════════════════════════════════════════════════════════


@runtime_checkable
class PartitionRule(Protocol):
    """Anything that partitions a geometry's intervals into cells."""

    def partition(self, geometry: "StructuredGeometry") -> "Partition": ...


def _read_only(values: object, where: str) -> np.ndarray:
    array = np.array(values, dtype=float) + 0.0  # -0.0 becomes +0.0
    if array.ndim != 1:
        raise ValueError(f"{where} must be 1-D; got shape {array.shape}")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{where} must be finite; got {array}")
    array.flags.writeable = False
    return array


@dataclass(frozen=True, eq=False)
class Partition:
    r"""Per interval of a geometry, its cell edges and its cell measures.

    Parameters
    ----------
    edges : sequence of 1-D arrays
        Per interval, its :math:`n_k + 1` cell edges, strictly increasing;
        interval :math:`k`'s last edge is interval :math:`k + 1`'s first.
    measures : sequence of 1-D arrays
        Per interval, its :math:`n_k` cell measures, each positive.

    Equality is bitwise over both. A partition is checked against a
    geometry by :meth:`partition`: one interval per material interval,
    and the end edges of interval :math:`k` are the breakpoints
    :math:`r_k` and :math:`r_{k+1}`, bit for bit.
    """

    edges: tuple[np.ndarray, ...]
    measures: tuple[np.ndarray, ...]
    all_edges: np.ndarray = field(init=False, repr=False)
    """Every cell edge, in order, each shared breakpoint once (read-only)."""
    all_measures: np.ndarray = field(init=False, repr=False)
    """Every cell measure, in order (read-only): the mesh's volumes are this array."""

    def __post_init__(self) -> None:
        edges = tuple(
            _read_only(e, f"Partition.edges[{k}]") for k, e in enumerate(self.edges)
        )
        measures = tuple(
            _read_only(m, f"Partition.measures[{k}]") for k, m in enumerate(self.measures)
        )
        if not edges:
            raise ValueError("a Partition has at least one interval")
        if len(measures) != len(edges):
            raise ValueError(
                f"a Partition has one measure array per interval: "
                f"{len(edges)} edge arrays, {len(measures)} measure arrays"
            )
        for k, (e, m) in enumerate(zip(edges, measures, strict=True)):
            if len(e) < 2:
                raise ValueError(
                    f"Partition interval {k} needs at least one cell (2 edges); "
                    f"got {len(e)} edge(s)"
                )
            if not np.all(np.diff(e) > 0):
                raise ValueError(
                    f"Partition interval {k}: the cell edges must be strictly "
                    f"increasing; got {e}"
                )
            if len(m) != len(e) - 1:
                raise ValueError(
                    f"Partition interval {k}: {len(e) - 1} cell(s) need "
                    f"{len(e) - 1} measure(s), got {len(m)}"
                )
            if not np.all(m > 0):
                raise ValueError(
                    f"Partition interval {k}: a cell measure is positive; got {m}"
                )
        for k, (left, right) in enumerate(itertools.pairwise(edges)):
            if left[-1] != right[0]:
                raise ValueError(
                    f"Partition intervals {k} and {k + 1} must share their "
                    f"breakpoint: {left[-1]!r} against {right[0]!r}"
                )
        object.__setattr__(self, "edges", edges)
        object.__setattr__(self, "measures", measures)
        all_edges = np.concatenate([edges[0], *(e[1:] for e in edges[1:])])
        all_measures = np.concatenate(measures)
        all_edges.flags.writeable = False
        all_measures.flags.writeable = False
        object.__setattr__(self, "all_edges", all_edges)
        object.__setattr__(self, "all_measures", all_measures)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Partition):
            return NotImplemented
        return all(
            len(mine) == len(theirs)
            and all(a.tobytes() == b.tobytes() for a, b in zip(mine, theirs, strict=True))
            for mine, theirs in ((self.edges, other.edges), (self.measures, other.measures))
        )

    __hash__ = None  # type: ignore[assignment]  # content identity is P1 step 5's

    @property
    def cell_counts(self) -> tuple[int, ...]:
        """The number of cells in each interval."""
        return tuple(len(m) for m in self.measures)

    def partition(self, geometry: "StructuredGeometry") -> "Partition":
        """This partition, once it is checked to refine ``geometry``'s intervals."""
        breakpoints = geometry.breakpoints
        if len(self.edges) != len(breakpoints) - 1:
            raise ValueError(
                f"the partition has {len(self.edges)} interval(s) and the "
                f"geometry {len(breakpoints) - 1}: a partition refines the "
                f"geometry's intervals one for one"
            )
        for k, e in enumerate(self.edges):
            if e[0] != breakpoints[k] or e[-1] != breakpoints[k + 1]:
                raise ValueError(
                    f"partition interval {k} spans [{e[0]!r}, {e[-1]!r}], not the "
                    f"geometry's interval [{breakpoints[k]!r}, {breakpoints[k + 1]!r}]: "
                    f"its end edges must be the breakpoints, bit for bit"
                )
        return self

    def __rmul__(self, factor: object) -> "Partition":
        raise TypeError(
            "a Partition has no spacing rule to place a refinement's new edges; "
            "refine the rule that produced it (2 * CellsByCount(...))"
        )


# ═══════════════════════════════════════════════════════════════════════
# The rules
# ═══════════════════════════════════════════════════════════════════════


def _per_interval[T](value: "T | tuple[T, ...]", n_intervals: int, what: str) -> tuple[T, ...]:
    """A scalar broadcast to every interval, or a tuple with one entry each."""
    if isinstance(value, tuple):
        if len(value) != n_intervals:
            raise ValueError(
                f"{what}: {len(value)} entries for {n_intervals} interval(s); "
                f"give one per interval, or one value for all"
            )
        return value
    return (value,) * n_intervals


def _parse_count(value: object, where: str) -> int:
    if not isinstance(value, (int, np.integer)) or isinstance(value, bool):
        raise TypeError(f"{where}: a cell count is an int, got {type(value).__name__}")
    if value < 1:
        raise ValueError(f"{where}: a cell count is at least 1, got {value}")
    return int(value)


def _parse_width(value: object, where: str) -> float:
    if not isinstance(value, Real) or isinstance(value, bool):
        raise TypeError(f"{where}: a width is a real number, got {type(value).__name__}")
    width = float(value)
    if not (math.isfinite(width) and width > 0):
        raise ValueError(f"{where}: a width is positive and finite, got {value}")
    return width


def _parse_spacing(spacing: object, where: str) -> Spacing:
    if not isinstance(spacing, Spacing):
        raise TypeError(
            f"{where}: the spacing is EqualWidth() or EqualVolume(), got "
            f"{type(spacing).__name__}"
        )
    return spacing


def _parse_factor(factor: object) -> int:
    """A refinement factor: a power of two, the only factors whose cells nest.

    The fine fractions are ``j * fl(1 / (k n))``; for ``k = 2^m`` that is
    ``(j / k) * fl(1 / n)`` exactly (a power-of-two scaling), so every coarse
    edge is a fine edge bit for bit. For any other ``k`` they differ by up to
    2 ulp (`[M]` 3909 of 7176 cases at k = 3), and cells that do not nest
    are not a refinement.
    """
    if not isinstance(factor, (int, np.integer)) or isinstance(factor, bool) or factor < 1:
        raise TypeError(f"a refinement factor is a positive int, got {factor!r}")
    if factor & (factor - 1):
        raise ValueError(
            f"a refinement factor is a power of two, got {factor}: only then does "
            f"every coarse edge stay a fine edge bit for bit"
        )
    return int(factor)


@dataclass(frozen=True)
class CellsByCount:
    r"""A number of cells per interval, placed by a spacing rule.

    Parameters
    ----------
    counts : int or tuple of int
        One count for every interval, or one per interval.
    spacing : Spacing
        :class:`EqualWidth` or :class:`EqualVolume`; there is no default.
    """

    counts: int | tuple[int, ...]
    spacing: Spacing

    def __post_init__(self) -> None:
        if isinstance(self.counts, tuple):
            counts: int | tuple[int, ...] = tuple(
                _parse_count(c, f"CellsByCount.counts[{k}]") for k, c in enumerate(self.counts)
            )
        else:
            counts = _parse_count(self.counts, "CellsByCount.counts")
        object.__setattr__(self, "counts", counts)
        _parse_spacing(self.spacing, "CellsByCount.spacing")

    @classmethod
    def uniform_width(cls, n: int) -> "CellsByCount":
        """``n`` cells of equal width in every interval."""
        return cls(n, EqualWidth())

    @classmethod
    def uniform_volume(cls, n: int) -> "CellsByCount":
        """``n`` cells of equal measure in every interval."""
        return cls(n, EqualVolume())

    def partition(self, geometry: "StructuredGeometry") -> Partition:
        counts = _per_interval(self.counts, len(geometry.mat_ids), "CellsByCount.counts")
        cells = [
            self.spacing.cells(geometry.coord, a, b, n)
            for (a, b), n in zip(itertools.pairwise(geometry.breakpoints), counts, strict=True)
        ]
        return Partition(
            edges=tuple(e for e, _ in cells), measures=tuple(m for _, m in cells),
        ).partition(geometry)

    def __rmul__(self, factor: object) -> "CellsByCount":
        k = _parse_factor(factor)
        if isinstance(self.counts, tuple):
            return CellsByCount(tuple(k * c for c in self.counts), self.spacing)
        return CellsByCount(k * self.counts, self.spacing)


@dataclass(frozen=True)
class CellsByMaxWidth:
    r"""The fewest cells per interval whose widest cell is no wider than a bound.

    In each interval the count is the least :math:`n` whose nominal widest
    cell (:meth:`Spacing.widest_cell`) is at most the bound :math:`h`.

    Parameters
    ----------
    widths : float or tuple of float
        The bound :math:`h`, for every interval or one per interval.
    spacing : Spacing
        :class:`EqualWidth` or :class:`EqualVolume`; there is no default.
    """

    widths: float | tuple[float, ...]
    spacing: Spacing

    def __post_init__(self) -> None:
        if isinstance(self.widths, tuple):
            widths: float | tuple[float, ...] = tuple(
                _parse_width(h, f"CellsByMaxWidth.widths[{k}]") for k, h in enumerate(self.widths)
            )
        else:
            widths = _parse_width(self.widths, "CellsByMaxWidth.widths")
        object.__setattr__(self, "widths", widths)
        _parse_spacing(self.spacing, "CellsByMaxWidth.spacing")

    def counts(self, geometry: "StructuredGeometry") -> tuple[int, ...]:
        """The cell count this rule realises in each of ``geometry``'s intervals."""
        widths = _per_interval(self.widths, len(geometry.mat_ids), "CellsByMaxWidth.widths")
        return tuple(
            self._least_count(a, b, h, geometry.coord)
            for (a, b), h in zip(itertools.pairwise(geometry.breakpoints), widths, strict=True)
        )

    def _least_count(self, a: float, b: float, h: float, coord: CoordSystem) -> int:
        p = self.spacing.power(coord)
        # The first cell is [a, T⁻¹(T(a) + ΔT/n)]; it is within h from this n on.
        n = max(1, math.ceil((b**p - a**p) / ((a + h) ** p - a**p)))
        while self.spacing.widest_cell(coord, a, b, n) > h:
            n += 1
        while n > 1 and self.spacing.widest_cell(coord, a, b, n - 1) <= h:
            n -= 1
        return n

    def partition(self, geometry: "StructuredGeometry") -> Partition:
        return CellsByCount(self.counts(geometry), self.spacing).partition(geometry)

    def __rmul__(self, factor: object) -> "Refined":
        return Refined(self, _parse_factor(factor))


@dataclass(frozen=True)
class Refined:
    r"""A rule's refinement: ``factor`` times the cells it realises in each interval.

    ``2 * CellsByMaxWidth(h, s)`` is this, not ``CellsByMaxWidth(h / 2, s)``:
    halving the bound can give an odd count, whose cells do not nest in the
    coarse ones.
    """

    rule: CellsByMaxWidth
    factor: int

    def partition(self, geometry: "StructuredGeometry") -> Partition:
        coarse = self.rule.counts(geometry)
        return CellsByCount(
            tuple(self.factor * n for n in coarse), self.rule.spacing,
        ).partition(geometry)

    def __rmul__(self, factor: object) -> "Refined":
        return Refined(self.rule, _parse_factor(factor) * self.factor)


@dataclass(frozen=True, eq=False)
class CellEdges:
    r"""The cell edges written out, per interval; the measures follow from them.

    For a partition no rule produces (an irregular grid whose irregularity
    is the point). Each measure is the shell between two given edges in the
    geometry's coordinate system.
    """

    edges: tuple[np.ndarray, ...]

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "edges",
            tuple(_read_only(e, f"CellEdges.edges[{k}]") for k, e in enumerate(self.edges)),
        )

    def partition(self, geometry: "StructuredGeometry") -> Partition:
        return Partition(
            edges=self.edges,
            measures=tuple(compute_volumes_1d(geometry.coord, e) for e in self.edges),
        ).partition(geometry)

    def __rmul__(self, factor: object) -> "CellEdges":
        raise TypeError(
            "CellEdges has no spacing rule to place a refinement's new edges"
        )


__all__ = [
    "CellEdges",
    "CellsByCount",
    "CellsByMaxWidth",
    "EqualVolume",
    "EqualWidth",
    "Partition",
    "PartitionRule",
    "Refined",
    "Spacing",
]
