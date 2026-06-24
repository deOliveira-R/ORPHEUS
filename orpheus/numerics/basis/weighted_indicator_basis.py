r"""The flux-/production-weighted cell indicator — the **test side** of a
Petrov-Galerkin homogenisation frame.

The third concrete :class:`~orpheus.numerics.basis.base.Basis`: a cell-indicator
basis whose **analysis** carries a solution-dependent nodal weight :math:`w`,

.. math::
   :label: weighted-indicator-test-functional

   \chi_R(x) \;=\; w(x)\,\mathbf{1}_R(x),

so the analysis functional is the :math:`w`-reweighted region integral
:math:`\langle\chi_R, f\rangle_W = \int_R w\,f\,\mathrm{d}V`.  It is the **test**
basis of the spatial homogenisation :class:`~orpheus.numerics.frame.PetrovGalerkinFrame`,
paired with the plain :class:`~orpheus.numerics.basis.indicator_basis.IndicatorBasis`
as the **trial** side: the flux-weighted effective cross section
:math:`\Sigma_R = \int_R \varphi\Sigma\,\mathrm{d}V / \int_R \varphi\,\mathrm{d}V`
is the frame's coefficient extraction :math:`G^{-1}M` (``frame.project``) with
``test = `` :math:`\varphi\cdot\mathbf 1_R`, ``trial = `` :math:`\mathbf 1_R`,
``measure = `` the geometric volume measure.

Why the weight is the TEST side, not a measure metric
=====================================================

The campaign's load-bearing physics ruling (GitHub #48): homogenisation is
**Petrov-Galerkin**, not Galerkin.  Folding :math:`\varphi` into the measure (the
"Galerkin in :math:`L^2(\varphi V)`" reading) is a forward-flux, reaction-rate-only
convenience that breaks under eigenvalue-consistent (adjoint-weighted)
homogenisation, where the preserved functional is the **bilinear**
:math:`\langle\varphi^*,\Sigma\varphi\rangle` and ``test`` :math:`= \varphi^*\cdot
\mathbf 1_R \ne` ``trial``.  The measure carries the axis + the fixed :math:`L^2`
metric and **never** the discipline; the solution-weighting :math:`w` (:math:`\varphi`
for forward homogenisation, the production density :math:`p=\sum_g\nu\Sigma_{f,g}
\varphi_g` for the :math:`\chi` emission collapse) lives here, on the test basis.

A test-only basis (the synthesis side is the trial's)
=====================================================

A :class:`~orpheus.numerics.frame.PetrovGalerkinFrame` takes its **reconstruction**
from the **trial** basis (:meth:`~orpheus.numerics.frame.FrameBase.reconstruction`
reads ``frame.basis``), so this basis is consumed ONLY through its
:meth:`evaluate` (the test table = the geometric membership, weight-free — the
weight is an *analysis* weight, not a tabulation) and :meth:`analyze` (the
:math:`w`-reweighted projection), plus :attr:`space` (the coarse coefficient
space, shared with the trial).  The synthesis-side operations (:meth:`synthesize`,
:meth:`reconstruct`, the representation transposes, :meth:`mass_matrix`) have **no
consumer** — building the weighted synthesis / weighted Hilbert adjoint is
speculation until one arrives — so they **raise** rather than silently delegate to
the *unweighted* indicator (which would be a half-consistent basis: a landmine
where ``analysis`` carries :math:`w` but its transpose does not).  When an adjoint
consumer arrives (e.g. a least-squares frame) the weighted transposes get built
then, consistently.

References
----------

* Brenner, S. C. and Scott, L. R. (2008). *The Mathematical Theory of Finite
  Element Methods*, 3rd ed. Springer. §3.4 — Petrov-Galerkin (test :math:`\ne`
  trial space).
* Hébert, A. (2009). *Applied Reactor Physics*. Polytechnique. §6 — flux-weighted
  cross-section homogenisation (the reactor-physics consumer of this test side).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from numpy.typing import NDArray

from orpheus.numerics.basis.base import Basis
from orpheus.numerics.basis.indicator_basis import IndicatorBasis

if TYPE_CHECKING:
    from orpheus.numerics.measure import DiscreteMeasure
    from orpheus.numerics.space import FunctionSpace


__all__ = ["WeightedIndicatorBasis"]


_TEST_ONLY = (
    "WeightedIndicatorBasis is the TEST (analysis) side of a Petrov-Galerkin "
    "frame; its {op} has no consumer and is not built. The frame's "
    "reconstruction is trial-side (the plain IndicatorBasis); folding the weight "
    "into an unbuilt {op} would make a half-consistent basis."
)


@dataclass(frozen=True, eq=False)
class WeightedIndicatorBasis(Basis):
    r"""Cell indicators reweighted by a nodal weight — the Petrov-Galerkin test side.

    The analysis functional is :math:`\langle\chi_R, f\rangle_W = \int_R w f\,
    \mathrm{d}V` for the test functions :math:`\chi_R = w\cdot\mathbf 1_R`
    :eq:`weighted-indicator-test-functional`.  The geometric support (which node
    lies in which cell) comes from the wrapped trial :class:`IndicatorBasis`; the
    weight :math:`w` is a nodal field aligned with the bound measure's nodes.

    Parameters
    ----------
    indicator : IndicatorBasis
        The trial cell-indicator basis — supplies the geometric membership table
        (:meth:`evaluate`) and the coarse coefficient :attr:`space`.
    weight : NDArray
        The nodal test weight :math:`w[n, \ldots]`, indexed by measure node on its
        leading axis.  Forward spatial homogenisation passes the per-group scalar
        flux :math:`\varphi` (shape ``(n_nodes, ng)``); the :math:`\chi` emission
        collapse passes the scalar production density :math:`p` (shape
        ``(n_nodes,)``).  Trailing axes ride :meth:`analyze`'s contraction and
        broadcast against the analysed field's leading axes (so a per-group flux
        weights a vector channel elementwise and a ``[g_from, g_to]`` matrix channel
        by its source group).

    Notes
    -----
    Frozen, ``eq=False`` (identity equality): the fields are a NumPy-array-bearing
    basis and a NumPy weight array, neither value-comparable by a dataclass — and an
    instance is a transient *view* a homogenisation builds, not a value compared for
    equality (matching :class:`IndicatorBasis`).
    """

    indicator: IndicatorBasis
    weight: NDArray

    # ── Tabulation: the geometric membership (the weight is an ANALYSIS weight) ──
    def evaluate(self, points: NDArray, /) -> NDArray:
        r"""Tabulate the geometric cell membership → the ``(N, n_cells)`` table.

        Identical to the trial :meth:`IndicatorBasis.evaluate` — **weight-free**.
        The weight :math:`w` enters the *analysis* contraction (:meth:`analyze`),
        NOT the tabulation: :math:`w` reweights the inner product, it does not move
        the support.  Keeping the table the plain 2-D membership is what lets the
        per-group flux ride :meth:`analyze`'s trailing axes instead of inflating the
        table to ``(N, n_cells, ng)``.
        """
        return self.indicator.evaluate(points)

    # ── The weighted analysis (the only weight-bearing contraction) ───────
    def analyze(
        self, values: NDArray, table: NDArray, weights: NDArray, /,
    ) -> NDArray:
        r"""Weighted analysis :math:`(M f)_R = \sum_{i\in R} w_i\, \tt{weights}_i\, f_i`.

        The region integral of the :math:`w`-reweighted integrand:
        :math:`\int_R w f\,\mathrm{d}V` with the measure (volume) ``weights``.  The
        weight is folded into the integrand (:meth:`_weighted`) and the contraction
        delegates to the trial :meth:`IndicatorBasis.analyze`, so for ``w == 1`` this
        is bit-identical to the unweighted region integral.
        """
        return self.indicator.analyze(self._weighted(values), table, weights)

    def _weighted(self, values: NDArray, /) -> NDArray:
        r"""Fold the nodal weight into a field by a **leading-aligned** broadcast.

        ``weight`` and ``values`` share their leading (node) axis; the shorter rank
        is padded with trailing size-1 axes so they broadcast over the union of
        trailing shapes.  This handles every homogenisation shape with one rule:

        * per-group flux ``w (n, ng)`` × vector channel ``f (n, ng)`` → ``(n, ng)``
          (elementwise per group);
        * per-group flux ``w (n, ng)`` × matrix channel ``f (n, ng, ng')`` →
          ``(n, ng, ng')`` (``w`` aligns to the **source** group ``ng``);
        * scalar production ``w (n,)`` × :math:`\chi` channel ``f (n, ng)`` →
          ``(n, ng)`` (``w`` broadcasts over the group axis);
        * the Gram probe ``w (n, ng)`` × the constant field ``𝟙 (n,)`` → ``(n, ng)``
          (``w`` longer than ``values``: the result acquires ``w``'s trailing shape,
          giving the per-group region Gram :math:`\Phi_{R,g}`).
        """
        w = self.weight
        nd = max(w.ndim, values.ndim)
        w = w.reshape(w.shape + (1,) * (nd - w.ndim))
        v = values.reshape(values.shape + (1,) * (nd - values.ndim))
        return w * v

    # ── The coarse coefficient space (shared with the trial) ──────────────
    @property
    def space(self) -> "FunctionSpace":
        r"""The coarse coefficient :class:`FunctionSpace` — the trial's, shape ``(n_cells,)``.

        The test and trial span the SAME coarse cells (one effective material per
        coarse cell), so the coefficient space is the trial indicator's; the
        Petrov-Galerkin distinction is the *analysis* weighting, not the coefficient
        layout.
        """
        return self.indicator.space

    # ── Synthesis side: no consumer (the PG reconstruction is trial-side) ──
    def synthesize(self, coefficients: NDArray, table: NDArray, /) -> NDArray:
        raise NotImplementedError(_TEST_ONLY.format(op="synthesize"))

    def reconstruct(self, coefficients: NDArray, table: NDArray, /) -> NDArray:
        raise NotImplementedError(_TEST_ONLY.format(op="reconstruct"))

    def reconstruct_transpose(self, values: NDArray, table: NDArray, /) -> NDArray:
        raise NotImplementedError(_TEST_ONLY.format(op="reconstruct_transpose"))

    def analyze_transpose(
        self, coefficients: NDArray, table: NDArray, weights: NDArray, /,
    ) -> NDArray:
        raise NotImplementedError(_TEST_ONLY.format(op="analyze_transpose"))

    def mass_matrix(self, measure: "DiscreteMeasure", /) -> NDArray:
        raise NotImplementedError(_TEST_ONLY.format(op="mass_matrix"))
