r"""Interval rules: how the cells of one interval of a geometry are placed.

A :class:`~orpheus.geometry.StructuredGeometry` cuts its interval of
positions at its breakpoints :math:`r_0 < r_1 < \dots < r_R` into material
intervals. A mesh refines that cut: it divides every interval
:math:`[r_k, r_{k+1}]` into cells, so every cell lies in exactly one
interval and every breakpoint is a cell edge. An **interval rule** says how
one interval is divided; the mesher applies one rule to every interval, or
one rule per interval, and builds the mesh from the cells.

A rule places edges; the **geometry** gives every measure
(:meth:`~orpheus.geometry.StructuredGeometry.measure`: lengths on a slab,
areas per unit height on a cylinder, volumes on a sphere), so no measure
exists apart from the geometry whose coordinate system defines it.

* :class:`CellsByCount`: a number of cells, placed by a spacing rule;
* :class:`CellsByMaxWidth`: the fewest cells whose nominal widest cell is
  no wider than a bound, placed by a spacing rule;
* :class:`Refined`: a counted rule's refinement, ``k * rule`` for ``k`` a
  power of two;
* :class:`CellEdges`: the cell edges written out.

The spacing rules
-----------------

A spacing rule places :math:`n` cells in :math:`[a, b]` at equal steps of a
:class:`~orpheus.geometry.coord.MeasureCoordinate` :math:`T`:

.. math::

    r_j = T^{-1}\!\left(T(a) + f_j\,\bigl(T(b) - T(a)\bigr)\right),
    \qquad f_j = j / n,

with both end edges pinned to the breakpoints. :class:`EqualWidth` steps in
:math:`T(r) = r` in every coordinate system; :class:`EqualVolume` steps in
the coordinate system's own measure coordinate (:math:`r`, :math:`r^2`,
:math:`r^3`), in which the measure is uniform. When the spacing's
coordinate is the measure coordinate the cells are equal shares of the
interval's measure, and each is stored as that share, :math:`m/n`: the
same volume recomputed from the realised edges differs by about one unit
in the last place (the ``sqrt``/``cbrt`` round trip), which is ERR-020.
Otherwise each measure is the geometry's measure of the realised cell. On a
slab both spacings step in :math:`r`, so they are one body and give one
mesh (#495, fixed at its root).
"""
from __future__ import annotations

import math
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import TYPE_CHECKING, Protocol, runtime_checkable

import numpy as np

from orpheus.geometry.coord import CoordSystem, MeasureCoordinate
from orpheus.numerics.content import ContentIdentity
from orpheus.geometry.scalars import (
    parse_integer,
    parse_positions,
    parse_positive_integer,
    parse_positive_real,
)

if TYPE_CHECKING:
    from orpheus.geometry.structured_geometry import StructuredGeometry

Interval = tuple[float, float]
Cells = tuple[np.ndarray, np.ndarray]

_MAX_CELLS = 2**31
"""The most cells a rule places in one interval: an index beyond int32 is no mesh."""


# ═══════════════════════════════════════════════════════════════════════
# The spacing rules
# ═══════════════════════════════════════════════════════════════════════


class Spacing(ABC):
    r"""The measure coordinate whose equal steps place an interval's cells."""

    @abstractmethod
    def coordinate(self, coord: CoordSystem) -> MeasureCoordinate:
        r"""The coordinate :math:`T` the cells are equal steps of, in ``coord``."""

    def cells(self, geometry: "StructuredGeometry", interval: Interval, n: int) -> Cells:
        r"""The ``n + 1`` edges and ``n`` measures of ``n`` cells in ``interval``."""
        a, b = interval
        T = self.coordinate(geometry.coord)
        t_a, t_b = T(np.array([a, b]))
        edges = T.inverse(t_a + np.linspace(0.0, 1.0, n + 1) * (t_b - t_a))
        edges[0], edges[-1] = a, b
        if T == geometry.coord.measure_coordinate:
            share = geometry.measure(np.array([a, b]))[0] / n
            return edges, np.full(n, share)
        return edges, geometry.measure(edges)

    def widest_cell(self, coord: CoordSystem, interval: Interval, n: int) -> float:
        r"""The nominal width of the widest of ``n`` cells in ``interval``.

        The two cases define "nominal" differently. Stepping in :math:`r`,
        every cell's nominal width is :math:`\mathrm{fl}((b - a)/n)` (the
        realised widths differ from it by up to 2 ulp of :math:`b`).
        Stepping in :math:`r^2` or :math:`r^3`, the cells narrow outward, and
        the nominal widest cell is the first one as the spacing body
        realises it.
        """
        a, b = interval
        T = self.coordinate(coord)
        if T.exponent == 1:
            return (b - a) / n
        t_a, t_b = T(np.array([a, b]))
        return float(T.inverse(t_a + (1.0 / n) * (t_b - t_a))) - a


@dataclass(frozen=True)
class EqualWidth(Spacing):
    r"""Cells of equal width: equal steps of :math:`r` in every coordinate system."""

    def coordinate(self, coord: CoordSystem) -> MeasureCoordinate:
        return MeasureCoordinate(1)


@dataclass(frozen=True)
class EqualVolume(Spacing):
    r"""Cells of equal measure: equal steps of the coordinate system's measure coordinate."""

    def coordinate(self, coord: CoordSystem) -> MeasureCoordinate:
        return coord.measure_coordinate


def _parse_spacing(spacing: object, where: str) -> Spacing:
    if not isinstance(spacing, Spacing):
        raise TypeError(
            f"{where}: the spacing is EqualWidth() or EqualVolume(), got "
            f"{type(spacing).__name__}"
        )
    return spacing


# ═══════════════════════════════════════════════════════════════════════
# The interval rules
# ═══════════════════════════════════════════════════════════════════════


@runtime_checkable
class IntervalRule(Protocol):
    """How one interval of a geometry is divided into cells."""

    def cells(self, geometry: "StructuredGeometry", interval: Interval) -> Cells:
        """The edges and the measures of the cells of ``interval``."""
        ...

    def __rmul__(self, factor: int) -> "IntervalRule":
        """The rule's refinement, ``factor`` a power of two (refused by rules that cannot refine)."""
        ...


class CountedRule(ABC):
    r"""An interval rule that places a number of cells by a spacing rule.

    ``k * rule`` is its refinement (:class:`Refined`).
    """

    spacing: Spacing
    """The spacing rule that places the cells."""

    @abstractmethod
    def count(self, geometry: "StructuredGeometry", interval: Interval) -> int:
        """The number of cells this rule places in ``interval``."""

    def cells(self, geometry: "StructuredGeometry", interval: Interval) -> Cells:
        return self.spacing.cells(geometry, interval, self.count(geometry, interval))

    def __rmul__(self, factor: object) -> "Refined":
        return Refined(self, _parse_factor(factor))


@dataclass(frozen=True)
class CellsByCount(CountedRule):
    r"""``n_cells`` cells, placed by a spacing rule (there is no default spacing)."""

    n_cells: int
    spacing: Spacing

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "n_cells",
            parse_positive_integer(self.n_cells, "CellsByCount.n_cells", "a cell count"),
        )
        _parse_spacing(self.spacing, "CellsByCount.spacing")

    @classmethod
    def uniform_width(cls, n: int) -> "CellsByCount":
        """``n`` cells of equal width."""
        return cls(n, EqualWidth())

    @classmethod
    def uniform_volume(cls, n: int) -> "CellsByCount":
        """``n`` cells of equal measure."""
        return cls(n, EqualVolume())

    def count(self, geometry: "StructuredGeometry", interval: Interval) -> int:
        return self.n_cells


@dataclass(frozen=True)
class CellsByMaxWidth(CountedRule):
    r"""The fewest cells whose nominal widest cell is no wider than ``width``.

    The count is the least :math:`n` with
    :meth:`Spacing.widest_cell` :math:`\le h`. The bound is on the nominal
    width; a realised width can exceed it by up to 2 ulp of the interval's
    outer end.
    """

    width: float
    spacing: Spacing

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "width",
            parse_positive_real(self.width, "CellsByMaxWidth.width", "a width"),
        )
        _parse_spacing(self.spacing, "CellsByMaxWidth.spacing")

    def count(self, geometry: "StructuredGeometry", interval: Interval) -> int:
        a, b = interval
        h, coord = self.width, geometry.coord
        T = self.spacing.coordinate(coord)
        t_a, t_ah, t_b = T(np.array([a, a + h, b]))
        if not t_ah > t_a:
            raise ValueError(
                f"CellsByMaxWidth: the width {h!r} is below the float resolution "
                f"of the interval [{a!r}, {b!r}] (a + h == a)"
            )
        # The first cell is [a, T⁻¹(T(a) + ΔT/n)]; it fits within h from this n on.
        estimate = (t_b - t_a) / (t_ah - t_a)
        if not estimate <= _MAX_CELLS:
            raise ValueError(
                f"CellsByMaxWidth: the width {h!r} would place about {estimate:.3g} "
                f"cells in [{a!r}, {b!r}], more than {_MAX_CELLS}"
            )
        n = max(1, math.ceil(estimate))
        while self.spacing.widest_cell(coord, interval, n) > h:
            n += 1
        while n > 1 and self.spacing.widest_cell(coord, interval, n - 1) <= h:
            n -= 1
        return n


def _parse_factor(factor: object) -> int:
    r"""A refinement factor: a power of two, the only factors whose cells nest.

    The fine fractions are ``j * fl(1 / (k n))``; for ``k = 2^m`` that is
    ``(j / k) * fl(1 / n)`` exactly (a power-of-two scaling), so every coarse
    edge is a fine edge bit for bit. For any other ``k`` they differ by up to
    2 ulp (`[M]` 3909 of 7176 cases at k = 3), and cells that do not nest
    are not a refinement.
    """
    k = parse_integer(factor, "Refined.factor", "a refinement factor")
    if k < 1 or k & (k - 1):
        raise ValueError(
            f"Refined.factor: a refinement factor is a power of two, got {k}: only "
            f"then does every coarse edge stay a fine edge bit for bit"
        )
    return k


@dataclass(frozen=True)
class Refined(CountedRule):
    r"""``factor`` times the cells ``rule`` places, by the same spacing.

    ``k * rule`` spells it. Every coarse edge is a fine edge, bit for bit,
    because ``factor`` is a power of two. ``2 * CellsByMaxWidth(h, s)`` is
    this, not ``CellsByMaxWidth(h / 2, s)``: halving the bound can give an
    odd count, whose cells do not nest in the coarse ones.
    """

    rule: CountedRule
    factor: int

    def __post_init__(self) -> None:
        if not isinstance(self.rule, CountedRule):
            raise TypeError(
                f"Refined.rule is a counted rule (CellsByCount, CellsByMaxWidth), "
                f"got {type(self.rule).__name__}"
            )
        factor = _parse_factor(self.factor)
        if isinstance(self.rule, Refined):  # k * (m * rule) is (k m) * rule
            factor *= self.rule.factor
            object.__setattr__(self, "rule", self.rule.rule)
        object.__setattr__(self, "factor", factor)

    @property
    def spacing(self) -> Spacing:  # type: ignore[override]  # derived, not stored
        return self.rule.spacing

    def count(self, geometry: "StructuredGeometry", interval: Interval) -> int:
        return self.factor * self.rule.count(geometry, interval)


@dataclass(frozen=True, eq=False)
class CellEdges(ContentIdentity):
    r"""One interval's cell edges, written out; the geometry gives their measures.

    For cells no rule places (a grid whose irregularity is the point). The
    first and last edges are the interval's breakpoints, bit for bit.
    """

    edges: np.ndarray

    def __post_init__(self) -> None:
        edges = parse_positions(self.edges, "CellEdges.edges")
        if len(edges) < 2 or not np.all(np.diff(edges) > 0):
            raise ValueError(
                f"CellEdges.edges are at least two strictly increasing positions; got {edges}"
            )
        object.__setattr__(self, "edges", edges)

    def cells(self, geometry: "StructuredGeometry", interval: Interval) -> Cells:
        a, b = interval
        if self.edges[0] != a or self.edges[-1] != b:
            raise ValueError(
                f"CellEdges span [{self.edges[0]!r}, {self.edges[-1]!r}], not the "
                f"interval [{a!r}, {b!r}]: the end edges are the breakpoints, bit for bit"
            )
        return self.edges.copy(), geometry.measure(self.edges)

    def __rmul__(self, factor: object) -> "CellEdges":
        raise TypeError(
            "CellEdges has no spacing rule to place a refinement's new edges; only a "
            "counted rule (CellsByCount, CellsByMaxWidth) refines"
        )


__all__ = [
    "CellEdges",
    "CellsByCount",
    "CellsByMaxWidth",
    "CountedRule",
    "EqualVolume",
    "EqualWidth",
    "IntervalRule",
    "Refined",
    "Spacing",
]
