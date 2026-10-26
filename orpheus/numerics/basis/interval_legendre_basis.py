r"""The Legendre polynomials on the reference interval, and the tensor product of bases.

**The interval basis.** :class:`IntervalLegendreBasis` is
:math:`\{P_k\}_{k < n}` on the reference interval :math:`\xi \in [-1, 1]`,
normalised by :math:`P_k(1) = 1`. Its Gram is taken against the NORMALISED
measure :math:`\mathrm d\xi / 2` (total mass 1), the per-volume convention of a
cell expansion:

.. math::

   \int_{-1}^{1} P_k P_l \,\frac{\mathrm d\xi}{2} = \frac{\delta_{kl}}{2k + 1},

so the coefficient weights are :math:`1/(2k+1)` and the canonical dual factor is
:math:`2k + 1` (DLMF Eq. 18.3, Table 18.3.1, with the measure halved). It is NOT
the angular :class:`~orpheus.numerics.basis.legendre_basis.LegendreBasis`, whose
functions are :math:`P_\ell(\Omega\cdot e_a)` on the orbit space
:math:`S^2/O(2)_a` with the pushforward Gram :math:`4\pi/(2\ell+1)`; its layout is a
:class:`~orpheus.numerics.coefficient_layout.PolynomialLayout`, not a degree layout,
so a spatial axis of the same shape cannot pass for an angular moment head.

Two facts of the family that cell schemes read: the endpoint values
:math:`P_k(\pm 1) = (\pm 1)^k` (the trace of an expansion on a cell face), and
the parity :math:`P_k(-\xi) = (-1)^k P_k(\xi)` (the action of the mirror
:math:`\xi \mapsto -\xi`).

**The tensor product.** :class:`TensorProductBasis` is the product
:math:`\phi_{k_0}(\xi_0)\cdots\phi_{k_{d-1}}(\xi_{d-1})` of orthogonal factor bases,
indexed by a :class:`~orpheus.numerics.coefficient_layout.TensorLayout`: factor 0
outermost (the Kronecker order). Its weights, dual factor and unit are the Kronecker
products of the factors'; its face along factor ``a`` is the product of the other
factors (:meth:`TensorProductBasis.face`).

The linear-discontinuous cell basis is the tensor product of ``d`` copies of the
two-term interval basis (:func:`spatial_moment_basis`): coefficients
:math:`(\bar\psi, \hat\psi)` per axis, :math:`\psi \approx \bar\psi + \xi\hat\psi`
(:eq:`discretization-ld-moments`, Larsen and Morel 1989 Eq. 4.1), in the order
:eq:`spatial-moment-kronecker-order` records.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property, reduce
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from orpheus.geometry.manifold import Interval, Product
from orpheus.numerics.basis.base import Basis, GramStructure
from orpheus.numerics.basis.legendre_basis import legendre_table
from orpheus.numerics.coefficient_layout import PolynomialLayout, TensorLayout

if TYPE_CHECKING:
    from orpheus.geometry.manifold import Manifold

__all__ = [
    "AVERAGE_MOMENT",
    "IntervalLegendreBasis",
    "SPATIAL_MOMENT_AXIS_LABEL",
    "TensorProductBasis",
    "spatial_moment_basis",
]

#: The reference interval of a cell expansion, :math:`[-1, 1]`.
_REFERENCE_INTERVAL = Interval(-1.0, 1.0)

#: The label of the within-cell spatial-moment axis: the label of the cell basis's
#: tensor layout, so the axis that :meth:`~orpheus.numerics.axis.Axis.of_basis`
#: mints from it carries this label.
SPATIAL_MOMENT_AXIS_LABEL = "spatial_moment"

#: The flat slot of the constant function in a tensor Legendre layout: the
#: index of the layout's unit :math:`e_0`, since :math:`P_0 = 1` on every
#: factor and the all-zero multi-index is first in the Kronecker order. The
#: slot of the cell average (and of the face average) in every moment buffer.
AVERAGE_MOMENT = 0


def _kron_all(vectors: list[NDArray]) -> NDArray:
    return reduce(np.kron, vectors, np.array([1.0]))


@dataclass(frozen=True)
class IntervalLegendreBasis(Basis):
    r"""The Legendre polynomials :math:`P_0, \dots, P_{n-1}` on :math:`[-1, 1]`, :math:`P_k(1) = 1`.

    Parameters
    ----------
    n : int
        The number of polynomials, so the highest degree is :math:`n - 1`
        (``n = 1``: the cell average alone; ``n = 2``: the average and the
        slope).
    """

    n: int

    def __post_init__(self) -> None:
        if int(self.n) != self.n or self.n < 1:
            raise ValueError(f"IntervalLegendreBasis: n must be a positive integer, got {self.n!r}")

    # ── The family's convention-bearing values ──────────────────────
    @cached_property
    def dual_factor(self) -> NDArray:
        r"""The canonical dual factor :math:`2k + 1`, shape ``(n,)``: reconstruction is synthesis scaled by it."""
        factor = 2.0 * np.arange(self.n) + 1.0
        factor.setflags(write=False)
        return factor

    def endpoint_values(self, side: int, /) -> NDArray:
        r""":math:`P_k(\pm 1) = (\pm 1)^k`, shape ``(n,)``: the trace of an expansion at the end ``side`` (+1 or -1)."""
        if side not in (1, -1):
            raise ValueError(f"IntervalLegendreBasis.endpoint_values: side is +1 or -1, got {side!r}")
        return legendre_table(self.n - 1, np.array([float(side)]))[0]

    @cached_property
    def parity(self) -> NDArray:
        r"""The action of the mirror :math:`\xi \mapsto -\xi`: :math:`P_k(-\xi) = (-1)^k P_k(\xi)`, shape ``(n,)``."""
        parity = (-1.0) ** np.arange(self.n)
        parity.setflags(write=False)
        return parity

    # ── Tabulation ──────────────────────────────────────────────────
    def evaluate(self, points: NDArray, /) -> NDArray:
        r"""The ``(N, n)`` table :math:`P_k(\xi_m)` at ``(N,)`` or ``(N, 1)`` points of :math:`[-1, 1]`."""
        xi = np.asarray(points, dtype=float)
        if not (xi.ndim == 1 or (xi.ndim == 2 and xi.shape[1] == 1)):
            raise ValueError(f"IntervalLegendreBasis.evaluate: expected (N,) or (N, 1) points of [-1, 1], got shape {xi.shape}")
        xi = xi.reshape(-1)
        if not self.domain.contains(xi):
            raise ValueError("IntervalLegendreBasis.evaluate: the points are not on the reference interval [-1, 1]")
        return legendre_table(self.n - 1, xi)

    # ── Table contractions ──────────────────────────────────────────
    def synthesize(self, coefficients: NDArray, table: NDArray, /) -> NDArray:
        r"""Naked synthesis :math:`\sum_k P_k(\xi_m)\, c_k`."""
        return np.einsum("nk,k...->n...", table, coefficients)

    def analyze(self, values: NDArray, table: NDArray, weights: NDArray, /) -> NDArray:
        r"""Analysis :math:`\sum_m w_m P_k(\xi_m)\, f_m`."""
        return np.einsum("n,nk,n...->k...", weights, table, values)

    def analyze_transpose(self, coefficients: NDArray, table: NDArray, weights: NDArray, /) -> NDArray:
        r"""Representation transpose :math:`w_m \sum_k P_k(\xi_m)\, c_k`."""
        return np.einsum("n,nk,k...->n...", weights, table, coefficients)

    def reconstruct(self, coefficients: NDArray, table: NDArray, /) -> NDArray:
        r"""Reconstruction :math:`\sum_k (2k+1) P_k(\xi_m)\, \varphi_k`."""
        return np.einsum("nk,k,k...->n...", table, self.dual_factor, coefficients)

    def reconstruct_transpose(self, values: NDArray, table: NDArray, /) -> NDArray:
        r"""Representation transpose :math:`(2k+1) \sum_m P_k(\xi_m)\, v_m`."""
        return np.einsum("nk,k,n...->k...", table, self.dual_factor, values)

    # ── Gram, domain, layout ────────────────────────────────────────
    @property
    def gram_structure(self) -> GramStructure:
        r"""DIAGONAL: the Legendre polynomials are orthogonal on :math:`[-1, 1]`."""
        return GramStructure.DIAGONAL

    @property
    def domain(self) -> Manifold:
        return _REFERENCE_INTERVAL

    @property
    def layout(self) -> PolynomialLayout:
        return PolynomialLayout(self.n)

    @cached_property
    def gram_per_degree(self) -> NDArray:
        r"""The continuum Gram :math:`1/(2k+1)` against :math:`\mathrm d\xi/2`, shape ``(n,)``."""
        weights = 1.0 / self.dual_factor
        weights.setflags(write=False)
        return weights

    @property
    def coefficient_weights(self) -> NDArray:
        r"""The continuum Gram diagonal, :attr:`gram_per_degree`."""
        return self.gram_per_degree


@dataclass(frozen=True)
class TensorProductBasis(Basis):
    r"""The tensor product of orthogonal interval bases, factor 0 outermost.

    Parameters
    ----------
    factors : tuple[IntervalLegendreBasis, ...]
        One basis per coordinate of the product domain, in coordinate order.
    label : str
        The label of the product's coefficient axis.
    """

    factors: tuple[IntervalLegendreBasis, ...]
    label: str

    def __post_init__(self) -> None:
        if not self.factors:
            raise ValueError("TensorProductBasis: a tensor product has at least one factor")
        if not all(isinstance(f, IntervalLegendreBasis) for f in self.factors):
            raise TypeError("TensorProductBasis: every factor is an IntervalLegendreBasis")
        object.__setattr__(self, "factors", tuple(self.factors))

    @property
    def ndim(self) -> int:
        """The number of factors, the product domain's dimension."""
        return len(self.factors)

    @property
    def dimension(self) -> int:
        """The number of coefficients: the product of the factors' counts."""
        return self.layout.shape[0]

    def face(self, axis: int, /) -> TensorProductBasis:
        r"""The basis of the face normal to ``axis``: the product of the other factors, in their order.

        The trace of an expansion on that face pairs the coefficients of factor
        ``axis`` with :meth:`IntervalLegendreBasis.endpoint_values` and leaves the
        transverse coefficients in this order. Refused on a one-factor product,
        whose face is a point.
        """
        self._require_axis(axis)
        if self.ndim == 1:
            raise ValueError("TensorProductBasis.face: the face of a one-factor product is a point, which carries no basis")
        return TensorProductBasis(self.factors[:axis] + self.factors[axis + 1:], self.label)

    def face_dimension(self, axis: int, /) -> int:
        r"""The dimension of the trace space on the face normal to ``axis``: the product of the other factors' counts.

        Total where :meth:`face` is not: the face of a one-factor product is a
        point, whose function space is the constants alone, so its dimension is
        the empty product, 1.
        """
        self._require_axis(axis)
        return int(np.prod([f.n for a, f in enumerate(self.factors) if a != axis], dtype=int))

    def _require_axis(self, axis: int) -> None:
        if not 0 <= axis < self.ndim:
            raise ValueError(f"TensorProductBasis: axis {axis} is not one of the {self.ndim} factor axes")

    def reflection_signs(self, mirrored: tuple[bool, ...], /) -> NDArray:
        r"""The diagonal action of the mirrors :math:`\xi_a \mapsto -\xi_a` on the coefficients.

        ``mirrored[a]`` says whether factor ``a`` is mirrored: each mirrored factor
        contributes its :attr:`IntervalLegendreBasis.parity`, each other factor
        ones; the result is their Kronecker product, shape ``layout.shape``.
        """
        if len(mirrored) != self.ndim:
            raise ValueError(f"TensorProductBasis.reflection_signs: one flag per factor ({self.ndim}), got {len(mirrored)}")
        return _kron_all([f.parity if m else np.ones(f.n) for f, m in zip(self.factors, mirrored)])

    @cached_property
    def dual_factor(self) -> NDArray:
        r"""The Kronecker product of the factors' dual factors."""
        return _kron_all([f.dual_factor for f in self.factors])

    # ── Tabulation ──────────────────────────────────────────────────
    def evaluate(self, points: NDArray, /) -> NDArray:
        r"""The ``(N, prod n_a)`` table at ``(N, d)`` points: row-wise Kronecker products of the factor tables."""
        pts = np.asarray(points, dtype=float)
        if self.ndim == 1 and pts.ndim == 1:
            pts = pts[:, None]
        if pts.ndim != 2 or pts.shape[1] != self.ndim:
            raise ValueError(f"TensorProductBasis.evaluate: expected (N, {self.ndim}) points, got shape {pts.shape}")
        table = np.ones((pts.shape[0], 1))
        for a, factor in enumerate(self.factors):
            factor_table = factor.evaluate(pts[:, a])
            table = np.einsum("ni,nj->nij", table, factor_table).reshape(pts.shape[0], -1)
        return table

    # ── Table contractions ──────────────────────────────────────────
    def synthesize(self, coefficients: NDArray, table: NDArray, /) -> NDArray:
        return np.einsum("nk,k...->n...", table, coefficients)

    def analyze(self, values: NDArray, table: NDArray, weights: NDArray, /) -> NDArray:
        return np.einsum("n,nk,n...->k...", weights, table, values)

    def analyze_transpose(self, coefficients: NDArray, table: NDArray, weights: NDArray, /) -> NDArray:
        return np.einsum("n,nk,k...->n...", weights, table, coefficients)

    def reconstruct(self, coefficients: NDArray, table: NDArray, /) -> NDArray:
        return np.einsum("nk,k,k...->n...", table, self.dual_factor, coefficients)

    def reconstruct_transpose(self, values: NDArray, table: NDArray, /) -> NDArray:
        return np.einsum("nk,k,n...->k...", table, self.dual_factor, values)

    # ── Gram, domain, layout ────────────────────────────────────────
    @property
    def gram_structure(self) -> GramStructure:
        r"""DIAGONAL: a product of orthogonal families is orthogonal."""
        return GramStructure.DIAGONAL

    @property
    def domain(self) -> Manifold:
        return reduce(Product, [f.domain for f in self.factors])

    @property
    def layout(self) -> TensorLayout:
        return TensorLayout(tuple(f.layout for f in self.factors), self.label)

    @cached_property
    def gram(self) -> NDArray:
        r"""The Kronecker product of the factors' continuum Grams, shape ``layout.shape``."""
        weights = _kron_all([f.coefficient_weights for f in self.factors])
        weights.setflags(write=False)
        return weights

    @property
    def coefficient_weights(self) -> NDArray:
        r"""The continuum Gram diagonal, :attr:`gram`."""
        return self.gram


def spatial_moment_basis(n_per_axis: int, ndim: int, /) -> TensorProductBasis:
    r"""The within-cell spatial-moment basis of a ``ndim``-dimensional cell: ``ndim`` interval Legendre factors.

    ``n_per_axis`` is the per-axis polynomial count (1 for the cell-average
    closures, 2 for linear discontinuous). Its layout's label is
    :data:`SPATIAL_MOMENT_AXIS_LABEL`.
    """
    return TensorProductBasis((IntervalLegendreBasis(n_per_axis),) * ndim, SPATIAL_MOMENT_AXIS_LABEL)
