r"""Piecewise-constant (P0) indicator basis on a tensor-product cell partition.

The second concrete :class:`~orpheus.numerics.basis.base.Basis` — the box /
characteristic-function basis whose functions are the **cell indicators**

.. math::
   :label: indicator-basis-functions

   \phi_R(x) \;=\; \mathbf{1}_R(x) \;=\;
   \begin{cases} 1 & x \in R \\ 0 & \text{otherwise,} \end{cases}

one per cell :math:`R` of a tensor-product partition given by per-axis edge
arrays.  It is the trial (synthesis) side of the **spatial / energy
homogenisation Frame**, exactly as
:class:`~orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis`
is the trial side of the angular Frame.  Binding it to a
:class:`~orpheus.numerics.measure.DiscreteMeasure` (a coarse mesh's
``volume_measure``) in a :class:`~orpheus.numerics.frame.FrameBase` makes the
fine→coarse projection that homogenises cross sections — flux-weighted
homogenisation is the Petrov-Galerkin case
(:class:`~orpheus.numerics.frame.PetrovGalerkinFrame`); see
:meth:`orpheus.sn.solution.Solution.homogenize`.

Why a basis at all (and why the *mesh* yields it)
=================================================

Homogenisation is the :math:`L^2(\phi V)`-orthogonal (Galerkin) projection of a
fine cross-section field onto the space of functions **piecewise-constant on the
coarse cells**.  That coarse space IS the span of the cell indicators
:eq:`indicator-basis-functions` — so the coarse mesh, viewed as a basis, is the
:math:`K` of Saad's :math:`(K, L)` projection pair.  The mesh does **not inherit**
:class:`Basis` (a :class:`Basis` is the *measure-free* half of a frame, but a mesh
carries the volume measure — inheriting would conflate the two roles); instead the
mesh *yields* this view, via ``mesh.indicator_basis()``, symmetric with how it
already yields ``mesh.volume_measure``.  The basis itself is **geometry-free**: it
holds only the per-axis edge arrays, so :mod:`orpheus.numerics` stays free of any
:mod:`orpheus.geometry` dependency.

The membership table is n-D by construction
-------------------------------------------

:meth:`evaluate` tabulates the indicators at sample points — a ``(N, n_cells)``
one-hot **membership table** ``T[i, R] = 1`` iff fine point ``i`` lies in coarse
cell ``R``.  It is built **per axis** (each axis's coordinate is bucketed into that
axis's edges by ``searchsorted``) and combined with
:func:`numpy.ravel_multi_index` in C / ``"ij"`` order — the *same* flat-cell
ordering ``DiscreteMeasure`` uses for its nodes
(``np.meshgrid(*centres, indexing="ij").ravel()``, see
:meth:`orpheus.transport.mesh.material_mesh.MaterialMesh.volume_measure`).  So the
table column index and the measure's node/weight index agree by construction, in
any dimension — there is no 1-D special case.

The dual factor is the identity; the Gram-inverse lives in the space metric
=========================================================================

For spherical harmonics the canonical-dual factor :math:`2\ell+1` is **analytic**
(measure-free), so it is folded into :meth:`reconstruct`.  For the indicator basis
the canonical dual of :math:`\mathbf{1}_R` is :math:`\mathbf{1}_R / m_R` with the
region mass :math:`m_R = \langle \mathbf{1}_R, \mathbf{1}_R\rangle_\mu =
\sum_{i\in R} w_i` — which is **measure-dependent**.  A :class:`Basis` is the
measure-free half of a frame, so a measure-dependent factor cannot live here:
:meth:`reconstruct` therefore uses the identity dual factor (it equals
:meth:`synthesize`, the plain broadcast), and the Gram-inverse normalisation
:math:`G^{-1} = \mathrm{diag}(1/m_R)` is applied **separately** by the coefficient
space's :meth:`~orpheus.numerics.space.FunctionSpace.apply_inverse_metric` (its
Moore–Penrose pseudo-inverse zeroes empty regions for free).  The orthogonal
projector onto the coarse space is the composite
:math:`R \circ G^{-1} \circ M` — broadcast ∘ inverse-Gram ∘ analysis — with the
inverse-Gram supplied by the bound measure, never by the basis.

References
----------

* Brenner, S. C. and Scott, L. R. (2008). *The Mathematical Theory of Finite
  Element Methods*, 3rd ed. Springer. §3.4 (P0 / piecewise-constant elements as a
  Galerkin trial space).
* Christensen, O. (2016). *An Introduction to Frames and Riesz Bases*, 2nd ed.
  Birkhäuser. §1 (the analysis/synthesis pair; the canonical dual frame
  :math:`\{S^{-1}\phi_k\}` whose factor is the inverse Gram).
* Hébert, A. (2009). *Applied Reactor Physics*. Polytechnique. §6 (flux-weighted
  cross-section homogenisation — the reactor-physics consumer of this projection).
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.basis.base import Basis, GramStructure

if TYPE_CHECKING:
    from orpheus.numerics.measure import DiscreteMeasure
    from orpheus.numerics.space import FunctionSpace


__all__ = ["IndicatorBasis"]


@dataclass(frozen=True, eq=False)
class IndicatorBasis(Basis):
    r"""Piecewise-constant cell-indicator basis on a tensor-product partition.

    Parameters
    ----------
    edges_per_axis : tuple[NDArray, ...]
        One sorted ``(n_a + 1,)`` array of cell edges per spatial axis.  The
        cells are the tensor product of the per-axis intervals; cell count is
        :math:`\prod_a n_a`.  A 1-D partition passes a length-1 tuple
        ``(edges,)``.  Held as arrays (not a mesh) so the basis is geometry-free.

    Notes
    -----
    Frozen, ``eq=False`` (identity equality): the fields are NumPy arrays, which
    have no value equality usable by a dataclass; an :class:`IndicatorBasis` is a
    transient *view* a mesh yields, not a value compared for equality.  (Contrast
    :class:`SphericalHarmonicBasis`, whose sole field ``L`` IS its value identity.)
    """

    edges_per_axis: tuple[NDArray, ...]

    # ── Gram structure: disjoint cells ⟹ diagonal Gram ────────────────────
    @property
    def gram_structure(self) -> GramStructure:
        r"""DIAGONAL — the cell indicators have disjoint support (:meth:`mass_matrix`)."""
        return GramStructure.DIAGONAL

    # ── Derived cell counts ───────────────────────────────────────────────
    @cached_property
    def cells_per_axis(self) -> tuple[int, ...]:
        r"""Cell count per axis :math:`(n_0, n_1, \ldots)` — the ``"ij"`` grid shape."""
        return tuple(int(edges.shape[0] - 1) for edges in self.edges_per_axis)

    @property
    def n_cells(self) -> int:
        r"""Total cell count :math:`\prod_a n_a` (the coefficient-space dimension)."""
        return int(np.prod(self.cells_per_axis)) if self.cells_per_axis else 0

    @property
    def ndim(self) -> int:
        """Number of spatial axes."""
        return len(self.edges_per_axis)

    # ── Tabulation (the only points-consuming method) ─────────────────────
    def evaluate(self, points: NDArray, /) -> NDArray:
        r"""Tabulate the cell indicators → the ``(N, n_cells)`` membership table.

        ``T[i, R] = 1`` iff sample point ``i`` lies in coarse cell ``R`` (zero
        otherwise) — a one-hot row per point.  Built per axis (each axis bucketed
        into its own edges with ``side="right"``: a point on an interior edge
        joins the right cell, a point on the outer edge the last cell) and flattened
        with :func:`numpy.ravel_multi_index` in C / ``"ij"`` order, so the column
        index matches the cell ordering of a ``DiscreteMeasure`` over the same
        partition.

        Parameters
        ----------
        points : NDArray, shape ``(N,)`` or ``(N, d)``
            Sample coordinates (typically the fine-mesh cell centres from a
            measure's ``nodes``).  ``d`` must equal :attr:`ndim`.

        Returns
        -------
        NDArray, shape ``(N, n_cells)``
            The one-hot membership table.
        """
        pts = np.asarray(points, dtype=float)
        if pts.ndim == 1:
            pts = pts[:, None]
        if pts.ndim != 2 or pts.shape[1] != self.ndim:
            raise ValueError(
                f"IndicatorBasis.evaluate expects points of shape (N,) or "
                f"(N, {self.ndim}); got {np.asarray(points).shape}."
            )
        per_axis_index = tuple(
            np.clip(
                np.searchsorted(edges, pts[:, axis], side="right") - 1,
                0, n_axis - 1,
            )
            for axis, (edges, n_axis) in enumerate(
                zip(self.edges_per_axis, self.cells_per_axis)
            )
        )
        flat_cell = np.ravel_multi_index(per_axis_index, self.cells_per_axis)
        table = np.zeros((pts.shape[0], self.n_cells))
        table[np.arange(pts.shape[0]), flat_cell] = 1.0
        return table

    # ── Table contractions (the Frame caches the table, delegates here) ───
    def synthesize(self, coefficients: NDArray, table: NDArray, /) -> NDArray:
        r"""Naked synthesis :math:`S_0(c)_i = \sum_R \mathbf{1}_R(x_i)\, c_R`.

        Broadcasts each cell's coefficient onto the points it contains — the
        piecewise-constant field.  Measure-free.
        """
        return np.einsum("nR,R...->n...", table, coefficients)

    def analyze(
        self, values: NDArray, table: NDArray, weights: NDArray, /,
    ) -> NDArray:
        r"""Analysis :math:`(M f)_R = \sum_{i\in R} w_i\, f(x_i)` — the region integral.

        The W-weighted projection onto the indicators: each coefficient is the
        ``weights``-weighted sum of ``values`` over the cell's points.  With the
        volume measure ``weights = V_i`` this is :math:`\int_R f\,\mathrm{d}V`;
        trailing axes (group, …) ride the einsum unchanged.
        """
        return np.einsum("n,nR,n...->R...", weights, table, values)

    def analyze_transpose(
        self, coefficients: NDArray, table: NDArray, weights: NDArray, /,
    ) -> NDArray:
        r"""Representation transpose :math:`(M^\top c)_i = w_i \sum_R \mathbf{1}_R(x_i)\, c_R`.

        The matrix transpose of :meth:`analyze` (:math:`= w_i \cdot S_0`) — the
        quadrature-weighted broadcast.  The metric-aware adjoint machinery combines
        it with the domain/codomain Gram to form the Hilbert adjoint, so the
        Frame's analysis face gets ``.H`` for free.
        """
        return np.einsum("n,nR,R...->n...", weights, table, coefficients)

    def reconstruct(self, coefficients: NDArray, table: NDArray, /) -> NDArray:
        r"""Reconstruction :math:`R`: coefficients → values, identity dual factor.

        For the indicator basis the canonical-dual factor :math:`1/m_R` (region
        mass) is **measure-dependent**, so it cannot be folded into this
        measure-free method (see the module docstring): :meth:`reconstruct` is the
        plain broadcast (:math:`d_R = 1`), equal to :meth:`synthesize`.  The
        Gram-inverse normalisation is applied separately by the coefficient space's
        :meth:`~orpheus.numerics.space.FunctionSpace.apply_inverse_metric`; the
        coarse-space projector is the composite :math:`R \circ G^{-1} \circ M`.
        """
        return np.einsum("nR,R...->n...", table, coefficients)

    def reconstruct_transpose(self, values: NDArray, table: NDArray, /) -> NDArray:
        r"""Representation transpose :math:`(R^\top v)_R = \sum_{i\in R} v_i`, identity dual factor.

        The matrix transpose of :meth:`reconstruct` (:math:`d_R = 1`); the naked
        analysis :math:`S_0^\top`.  Measure-free, symmetric with :meth:`reconstruct`.
        """
        return np.einsum("nR,n...->R...", table, values)

    # ── The discrete Gram + the coefficient space ─────────────────────────
    def mass_matrix(self, measure: "DiscreteMeasure", /) -> NDArray:
        r"""Discrete Gram :math:`g_{RS} = \sum_n w_n\, \mathbf{1}_R(x_n)\, \mathbf{1}_S(x_n)` over ``measure``.

        **Diagonal** — the indicators have disjoint support, so
        :math:`g_{RS} = \delta_{RS}\, m_R` with the region mass
        :math:`m_R = \sum_{i\in R} w_i` (the cell's measure-weighted size).  The
        diagonal is the inverse of the canonical-dual factor; its reciprocal is the
        :meth:`~orpheus.numerics.space.FunctionSpace.apply_inverse_metric` weight a
        Frame applies to normalise the projection.

        Returns
        -------
        NDArray, shape ``(n_cells, n_cells)``
            The diagonal Gram against ``measure``.
        """
        table = self.evaluate(measure.nodes)
        return np.einsum("n,nR,nS->RS", measure.weights, table, table)

    @property
    def space(self) -> "FunctionSpace":
        r"""The coarse coefficient :class:`FunctionSpace`, shape ``(n_cells,)``.

        Euclidean (no intrinsic metric): the indicator basis is **measure-free**,
        so its Gram :math:`\mathrm{diag}(m_R)` exists only against a bound measure
        (:meth:`mass_matrix`), not as a standalone basis property — unlike the
        analytic spherical-harmonic Gram.  A consumer that needs the metric installs
        it from the bound measure (e.g. homogenisation installs the flux·volume
        region masses :math:`\Phi_R` via :func:`dataclasses.replace`).
        """
        from orpheus.numerics.space import FunctionSpace

        return FunctionSpace(
            name=f"L2[coarse_cells_R{self.ndim}]", shape=(self.n_cells,),
        )
