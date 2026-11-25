r"""Constructors the N03 gates share (spec 7.4). Not a test module: every function RETURNS data and none asserts.

``python -O`` strips a bare ``assert`` of a helper, so the collected test files assert. The new names (``Factor``, ``Face``, ``NoForm``,
``Variance``) are read from their module AT CALL TIME, inside each function, so that on a tree without them a gate row fails in its body
(a behavioural red a reader can interpret) and not at collection.
"""

from __future__ import annotations

import itertools
from typing import Any

import numpy as np


def pairs(n: int) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    """Every ``(S, K)`` of coordinate subsets of ``range(n)`` with ``K ⊆ S``: ``3**n`` pairs."""
    return [
        (S, K)
        for r in range(n + 1)
        for S in itertools.combinations(range(n), r)
        for k in range(len(S) + 1)
        for K in itertools.combinations(S, k)
    ]


def variances() -> list[Any]:
    """The four values of ``Variance``: the Klein four-group of (kind, chart)."""
    from orpheus.numerics import factor as F

    return [F.Variance(k, c) for k in F.Kind for c in F.VarianceChart]


def diag_axis(weights: list[float] | np.ndarray, label: str = "a") -> Any:
    from orpheus.numerics.axis import Axis, BasisKind, Declared
    from orpheus.numerics.metric import DiagonalMetric

    w = np.asarray(weights, dtype=float)
    return Axis(label, (len(w),), origin=Declared(BasisKind.NODAL, DiagonalMetric(w)))


def dense_axis(n: int = 4, seed: int = 3, deficient: bool = True) -> tuple[Any, Any, np.ndarray]:
    """``(axis, kernel, G)``: a dense symmetric positive semi-definite ``G`` of rank ``n - 1`` (or full rank), stated as the metric face."""
    from orpheus.numerics.axis import Axis, BasisKind, Declared
    from orpheus.numerics.declaration import Coord, Spanned
    from orpheus.numerics.metric import DenseMetric, Face

    rng = np.random.default_rng(seed)
    c = rng.standard_normal((n, n - 1 if deficient else n))
    g = c @ c.T
    g = (g + g.T) / 2
    kernel = Spanned(n, ("test-gram-kernel", seed), 1) if deficient else Coord(n, ())
    axis = Axis("h", (n,), origin=Declared(BasisKind.NODAL, DenseMetric(g, kernel=kernel, face=Face.METRIC)))
    return axis, kernel, g


def dense_coord_axis() -> tuple[Any, Any, np.ndarray]:
    """``(axis, kernel, G)``: a dense symmetric positive semi-definite ``G`` on ``R^4`` with an exactly zero row and column at index 2, whose
    declared kernel is the non-empty ``Coord(4, (2,))`` (the exact ``Coord`` admission)."""
    from orpheus.numerics.axis import Axis, BasisKind, Declared
    from orpheus.numerics.declaration import Coord
    from orpheus.numerics.metric import DenseMetric, Face

    g = np.array([[2.0, 0.5, 0.0, 1.0], [0.5, 3.0, 0.0, 0.25], [0.0, 0.0, 0.0, 0.0], [1.0, 0.25, 0.0, 4.0]])
    kernel = Coord(4, (2,))
    return Axis("h", (4,), origin=Declared(BasisKind.NODAL, DenseMetric(g, kernel=kernel, face=Face.METRIC))), kernel, g


def density(axis: Any, sk: Any) -> Any:
    """The density factor ``(MEASURE, DENSITY)`` of ``sk`` on ``axis``."""
    from orpheus.numerics import factor as F

    return F.Factor(axis, F.Variance(F.Kind.MEASURE, F.VarianceChart.DENSITY), sk)


def sk_of(n: int, S: tuple[int, ...], K: tuple[int, ...]) -> Any:
    from orpheus.numerics import factor as F
    from orpheus.numerics.declaration import Coord

    return F.Subquotient(Coord(n, S), Coord(n, K))
