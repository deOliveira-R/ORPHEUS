r"""Coefficient layouts — how a basis indexes its coefficients.

A basis's coefficient space is an axis whose index set the basis fixes: the
real spherical harmonics fill a rectangular ``(L+1, 2L+1)`` table, the
Legendre basis on an orbit space :math:`S^2/O(2)_a` a flat ``(L+1,)`` row (one
coefficient per degree: the trivial isotypic component of :math:`O(2)_a` is
one-dimensional in every degree), and the indicator families one coefficient
per cell of a partition. That index structure is a VALUE the basis carries
(ruling 81 of the numerics untangling): the basis ranks below the axis in the
concept order, so it cannot mint its axis, and the one generic coefficient
axis reads the layout instead (:attr:`orpheus.numerics.axis.Axis.layout`).
No consumer dispatches on the family: each layout answers the reads its
consumers make.

The degree layouts (:class:`HarmonicLayout`, :class:`LegendreLayout`) answer
what a moment head's consumers read: the truncation order ``L``, the index
tuple of the isotropic coefficient (``isotropic_slot``), the index tuple of
the degree-:math:`\ell` block (``degree_block``), and ``unit``, the
coefficients of the constant function
:math:`1` (both families normalise :math:`Y_{00} = P_0 = 1`, so it is
:math:`e_0` at the isotropic slot; ruling 37). Until S7 these lived on two
space classes, ``SphericalHarmonicSpace`` and ``LegendreSpace``, behind the
``MomentHead`` protocol.

:class:`CellLayout` is the indicator families' index: one coefficient per cell
of a partition (no measure: a partition's cells carry volume only against
the frame's fine measure, so the raw cell head refuses to integrate, the
user's ruling, 2026-10-21). It carries the partition's per-axis EDGES, because the cells
are the atoms a partition is made of, and atoms are content (ruling 71): two
partitions with equal widths and shifted edges index two coefficient spaces.

:class:`ModeLayout` is a basis computed from data, a row of mode amplitudes
with no degree structure (the loss-kernel gauge's kernel modes); its unit is
supplied by the basis and may be zero.

The layouts are content values (:class:`~orpheus.values.content.ContentIdentity`):
two layouts are equal iff of one class with equal parts.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from orpheus.values.content import ContentIdentity

__all__ = ["CellLayout", "CoefficientLayout", "DegreeLayout", "HarmonicLayout", "LegendreLayout", "ModeLayout"]


def _e0(shape: tuple[int, ...]) -> NDArray:
    unit = np.zeros(shape)
    unit[(0,) * len(shape)] = 1.0
    unit.setflags(write=False)
    return unit


def _check_order(L: int, owner: str) -> int:
    if int(L) != L or L < 0:
        raise ValueError(f"{owner}: the truncation order must be a non-negative integer, got {L!r}")
    return int(L)


@dataclass(frozen=True, eq=False)
class HarmonicLayout(ContentIdentity):
    r"""The rectangular ``(L+1, 2L+1)`` table of the real spherical harmonics.

    Row :math:`\ell`, column :math:`\ell + m` (the addition-theorem shift), zero
    padding outside :math:`|m| \le \ell`.
    """

    L: int

    def __post_init__(self) -> None:
        object.__setattr__(self, "L", _check_order(self.L, "HarmonicLayout"))

    @property
    def label(self) -> str:
        return "harmonic"

    @property
    def shape(self) -> tuple[int, ...]:
        return (self.L + 1, 2 * self.L + 1)

    @property
    def isotropic_slot(self) -> tuple[int, ...]:
        return (0, 0)

    def degree_block(self, l: int, /) -> tuple[int | slice, ...]:
        r"""Row :math:`\ell`, its :math:`2\ell + 1` live columns."""
        if not 0 <= l <= self.L:
            raise ValueError(f"HarmonicLayout(L={self.L}): degree {l} is outside [0, {self.L}]")
        return (l, slice(0, 2 * l + 1))

    def spread(self, per_degree: NDArray) -> NDArray:
        r"""A per-degree value broadcast over each degree's live columns, zero on
        the padding (the continuum Gram :math:`4\pi/(2\ell+1)` on the table)."""
        values = np.asarray(per_degree, dtype=float)
        if values.shape != (self.L + 1,):
            raise ValueError(f"HarmonicLayout(L={self.L}): one value per degree, got shape {values.shape}")
        table = np.zeros(self.shape)
        for l in range(self.L + 1):
            table[self.degree_block(l)] = values[l]
        return table

    @property
    def unit(self) -> NDArray:
        return _e0(self.shape)


@dataclass(frozen=True, eq=False)
class LegendreLayout(ContentIdentity):
    r"""The flat ``(L+1,)`` row of the Legendre basis on an orbit space
    :math:`S^2/O(2)_a`: one coefficient per degree."""

    L: int

    def __post_init__(self) -> None:
        object.__setattr__(self, "L", _check_order(self.L, "LegendreLayout"))

    @property
    def label(self) -> str:
        return "legendre"

    @property
    def shape(self) -> tuple[int, ...]:
        return (self.L + 1,)

    @property
    def isotropic_slot(self) -> tuple[int, ...]:
        return (0,)

    def degree_block(self, l: int, /) -> tuple[int | slice, ...]:
        r"""The one coefficient of degree :math:`\ell`."""
        if not 0 <= l <= self.L:
            raise ValueError(f"LegendreLayout(L={self.L}): degree {l} is outside [0, {self.L}]")
        return (l,)

    def spread(self, per_degree: NDArray) -> NDArray:
        r"""A per-degree value on the flat row: the row itself."""
        values = np.asarray(per_degree, dtype=float)
        if values.shape != self.shape:
            raise ValueError(f"LegendreLayout(L={self.L}): one value per degree, got shape {values.shape}")
        return values.copy()

    @property
    def unit(self) -> NDArray:
        return _e0(self.shape)


@dataclass(frozen=True, eq=False)
class CellLayout(ContentIdentity):
    r"""One coefficient per cell of a tensor partition, flat, row-major.

    ``edges`` is the partition's per-axis cell edges (a tuple of strictly
    increasing 1-D arrays, stored read-only); the cell count is the product of
    the per-axis cell counts.
    """

    edges: tuple[NDArray, ...]

    def __post_init__(self) -> None:
        frozen = []
        for k, e in enumerate(self.edges):
            arr = np.array(e, dtype=float).ravel()
            if arr.size < 2:
                raise ValueError(f"CellLayout: axis {k} needs at least two edges, got {arr.size}")
            if not bool(np.isfinite(arr).all()) or not bool((np.diff(arr) > 0.0).all()):
                raise ValueError(f"CellLayout: axis {k}'s edges must be finite and strictly increasing, got {arr!r}")
            arr.setflags(write=False)
            frozen.append(arr)
        if not frozen:
            raise ValueError("CellLayout: a partition has at least one axis")
        object.__setattr__(self, "edges", tuple(frozen))

    @property
    def label(self) -> str:
        return "cells"

    @property
    def n_cells(self) -> int:
        return int(np.prod([e.size - 1 for e in self.edges]))

    @property
    def shape(self) -> tuple[int, ...]:
        return (self.n_cells,)


@dataclass(frozen=True, eq=False)
class ModeLayout(ContentIdentity):
    r"""A flat row of ``n`` mode amplitudes with no degree structure: a basis
    computed from data (the loss-kernel gauge's G-orthonormal kernel modes).

    ``unit`` is the coefficients of the constant function's projection onto
    the span, which the basis DECLARES (a computed projection would be
    round-off where the span is orthogonal to the constants, and the unit is
    content). It is ZERO for the loss-kernel gauge's modes (the user's
    rulings, 2026-10-21, at S7: a zero unit is legal, and the gauge's is
    declared), and then the marginal over the factor integrates every field
    to exactly 0, which is what each mode integrates to. ``label`` names the
    family. A frame dressing a mode basis keeps the declared unit
    (:class:`~orpheus.numerics.axis.FromFrame`).
    """

    unit: NDArray
    label: str

    def __post_init__(self) -> None:
        unit = np.array(self.unit, dtype=float)
        if unit.ndim != 1 or unit.size < 1:
            raise ValueError(f"ModeLayout: the unit is one coefficient per mode, got shape {unit.shape}")
        unit.setflags(write=False)
        object.__setattr__(self, "unit", unit)

    @property
    def shape(self) -> tuple[int, ...]:
        return (int(self.unit.size),)


#: The layouts with degree structure: the moment heads' families.
DegreeLayout = HarmonicLayout | LegendreLayout

#: Every coefficient layout a basis carries (ruling 81).
CoefficientLayout = DegreeLayout | CellLayout | ModeLayout
