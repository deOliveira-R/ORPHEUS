r"""The :class:`Basis` ABC — a discrete spectral basis on a measure space.

A *basis* is the **synthesis (trial) side** of a discrete frame
(:class:`orpheus.numerics.frame.FrameBase`): a collection of functions that can be
tabulated at sample points, reconstructed from coefficients, and weighed against
a quadrature to form a discrete Gram. It is the choice-free, measure-free half of
the analysis/synthesis pair — the :class:`~orpheus.numerics.measure.DiscreteMeasure`
supplies the analysis (test) side, and the :class:`FrameBase` binds the two.

Why the ABC is un-deferred on a single concrete basis
=====================================================

The package previously deferred a formal ``Basis`` ABC "until a second concrete
basis arrives" (``feedback_unify_after_two_instances``). It is promoted now — with
:class:`~orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis`
still the only concrete member — because the **forcing consumer has arrived**: the
generic :class:`~orpheus.numerics.frame.FrameBase` binds an *abstract* basis to a
measure and applies a non-identity morphism (analysis / reconstruction) across the
interface. The Frame needs the interface, not a second basis; the interface itself
is math-rigid (Grand Report v3 §5.4 lists the nine eventual bases — real spherical
harmonics, Chebyshev, Lagrange, FE shape functions, ...) and every one of them
tabulates, reconstructs, and has a discrete Gram. Deferring the ABC would force the
Frame to name a concrete basis, which is exactly the coupling the abstraction exists
to break.

The contract — tabulate, then contract a cached table
=====================================================

:meth:`evaluate` is the only method that takes sample *points*; it produces the
``Φ(point, mode)`` **table** (the layout-bearing object — ``(N, ℓ, m)`` for spherical
harmonics, ``(N, K)`` for a flat-mode basis). Every other operation **consumes that
table**, so the :class:`FrameBase` evaluates ONCE (``frame.table``) and the per-apply hot
path never re-tabulates. The basis owns these contractions because the index layout is
the basis's own; the :class:`FrameBase` stays layout-agnostic and merely delegates.

The shared naked synthesis :math:`S_0`
--------------------------------------

:meth:`synthesize` is the naked synthesis :math:`S_0(c) = \sum_k \phi_k\, c_k` (the
frame-theory *synthesis operator* :math:`T^*` — NO weights, NO dual factor). The four
weighted operations are each :math:`S_0` (or its transpose) post-multiplied by ONE
diagonal weight family:

* :meth:`analyze` :math:`M = w_n \cdot (\text{analysis contraction})` — the analysis
  operator :math:`T`;
* :meth:`analyze_transpose` :math:`M^\top = w_n \cdot S_0` — its representation transpose;
* :meth:`reconstruct` :math:`R = d_k \cdot S_0` — the canonical-dual synthesis (for the
  SH basis :math:`d_\ell = 2\ell+1`);
* :meth:`reconstruct_transpose` :math:`R^\top = d_k \cdot S_0^\top` — its representation
  transpose (measure-free, like :meth:`reconstruct`).

The two analysis-side weightings bake in the quadrature weight :math:`w_n`; the two
synthesis-side ones carry only the dual factor :math:`d_k` — synthesis is measure-free.

They are kept as **fused contractions** (not ``weight ⊙ synthesize``) because FP
non-associativity makes the factored form drift at the ULP level, and the scattering
kernel is pinned at 0 ULP. The shared :math:`S_0` is the *conceptual* unity, documented
here; the implementation keeps the fused einsum.

References
----------

* Grand Report v3 §5.4 — Basis hierarchy.
* Christensen, O. (2016). *An Introduction to Frames and Riesz Bases*, 2nd ed.
  Birkhäuser — the analysis/synthesis operator pair this ABC is the trial side of.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from enum import Enum
from typing import TYPE_CHECKING

from numpy.typing import NDArray

if TYPE_CHECKING:
    from orpheus.numerics.measure import DiscreteMeasure
    from orpheus.numerics.space import FunctionSpace


__all__ = ["Basis", "GramStructure"]


class GramStructure(Enum):
    r"""The structure of a TRIAL basis's frame-Gram — the projection-validity declaration.

    When a basis is the **trial** side of a :class:`~orpheus.numerics.frame.FrameBase`,
    the coefficient extraction :meth:`~orpheus.numerics.frame.FrameBase.project`
    (:math:`G^{-1}M`) normalises by the cross Gram :math:`G = MR`. The frame computes
    :math:`G` with a single **row-sum probe** (``analysis(reconstruction(ones))``) — but
    that probe equals the projection's required normalisation ONLY under one of two
    structural conditions. This enum is the basis's declaration of which (if any) holds,
    so the wrong combination cannot be spelled silently (a precondition the **type**
    enforces, not a docstring):

    * :attr:`DIAGONAL` — disjoint-support trial (orthogonal harmonics, nested cell /
      group indicators): :math:`MR` is diagonal, the row sum IS the diagonal.
    * :attr:`PARTITION_OF_UNITY` — overlapping trial whose membership rows sum to 1
      (the fractional :class:`~orpheus.numerics.basis.overlap_basis.OverlapBasis`):
      :math:`MR` is NOT diagonal, but :math:`R\mathbf 1 = \mathbf 1` collapses the probe
      to the per-region weight anyway.
    * :attr:`DENSE` — neither (a tapered weight, a higher-rank GEC moment — GitHub #275):
      the row-sum probe is wrong; the true projection needs the dense :math:`(MR)^{-1}M`
      solve, which is NOT built. :meth:`~orpheus.numerics.frame.FrameBase.project`
      **refuses** a DENSE-Gram trial rather than return a silently-wrong coarsening.

    The base :class:`Basis` defaults to :attr:`DENSE` — the safe refusal: a new basis
    must consciously declare it is row-sum-collapsible (having checked its Gram) before
    ``project`` will use the shortcut on it.
    """

    DIAGONAL = "diagonal"
    PARTITION_OF_UNITY = "partition_of_unity"
    DENSE = "dense"


class Basis(ABC):
    r"""Abstract discrete spectral basis — the synthesis side of a :class:`FrameBase`.

    Concrete bases (real spherical harmonics today; Legendre, Chebyshev,
    Lagrange/FE shape functions to come) implement the operations below. A basis
    is *choice-free*: it knows its functions and their convention, but not which
    quadrature samples them — that choice is the
    :class:`~orpheus.numerics.measure.DiscreteMeasure`, bound to the basis by a
    :class:`~orpheus.numerics.frame.FrameBase`.

    The contraction methods take the ``table`` from :meth:`evaluate` (positional,
    so concrete bases may name their arguments in their own domain vocabulary —
    e.g. ``directions``), so the :class:`FrameBase` tabulates once and the per-apply
    hot path never re-evaluates.
    """

    # ── Tabulation (the only points-consuming method) ─────────────────────
    @abstractmethod
    def evaluate(self, points: NDArray, /) -> NDArray:
        r"""Tabulate the basis functions at ``points`` → the ``Φ(point, mode)`` table.

        For the spherical-harmonic basis, ``Y[n, ℓ, ℓ+m] = Y_ℓ^m(Ω̂_n)`` at the
        ``(N, 3)`` direction cosines; a flat-mode basis returns ``(N, K)``.
        """
        ...

    # ── Table contractions (the Frame caches the table and delegates here) ──
    @abstractmethod
    def synthesize(self, coefficients: NDArray, table: NDArray, /) -> NDArray:
        r"""Naked synthesis :math:`S_0(c) = \sum_k \phi_k\, c_k` against a cached ``table``.

        The pure (frame-theory *synthesis operator* :math:`T^*`) reconstruction —
        NO measure weights, NO dual-frame factor. The shared kernel the three
        weighted contractions below are each one diagonal away from.
        """
        ...

    @abstractmethod
    def analyze(
        self, values: NDArray, table: NDArray, weights: NDArray, /,
    ) -> NDArray:
        r"""Analysis :math:`M = T`: sampled values → coefficients.

        :math:`(M f)_k = \sum_n w_n\, \phi_k(x_n)\, f(x_n)` — the W-weighted
        projection onto the basis. For spherical harmonics
        :math:`\phi_\ell^m = \sum_n w_n Y_\ell^m(\hat\Omega_n)\,\psi_n`. The
        ``weights`` are the measure's (analysis is the *measured* / test side).
        """
        ...

    @abstractmethod
    def analyze_transpose(
        self, coefficients: NDArray, table: NDArray, weights: NDArray, /,
    ) -> NDArray:
        r"""Representation transpose :math:`M^\top = w_n \cdot S_0`: coefficients → values.

        The matrix transpose of :meth:`analyze` (NOT its Hilbert adjoint): the
        naked synthesis weighted by the quadrature weight on each node. The
        metric-aware ``_AdjointOperator`` machinery combines this with the
        domain/codomain Gram to form the W-weighted Hilbert adjoint
        :math:`M^* = g_C \cdot S_0`, so the :class:`FrameBase`'s analysis face gets
        ``.H`` for free.
        """
        ...

    @abstractmethod
    def reconstruct(self, coefficients: NDArray, table: NDArray, /) -> NDArray:
        r"""Reconstruction :math:`R`: coefficients → values (the dual-frame synthesis).

        :math:`S_0` weighted by the canonical-dual factor (intrinsic to the basis;
        for spherical harmonics :math:`R = \sum_\ell (2\ell+1) \sum_m Y_\ell^m\,
        \phi_\ell^m` — the addition-theorem reconstruction). **Measure-free** (the
        dual factor needs no quadrature) — synthesis is the choice-free / trial side.
        """
        ...

    @abstractmethod
    def reconstruct_transpose(
        self, values: NDArray, table: NDArray, /,
    ) -> NDArray:
        r"""Representation transpose :math:`R^\top = d_k \cdot S_0^\top`: values → coefficients.

        The matrix transpose of :meth:`reconstruct` (NOT its Hilbert adjoint):
        :math:`(R^\top v)_k = d_k \sum_n \phi_k(x_n)\, v_n`, the naked
        analysis weighted by the canonical-dual factor :math:`d_k` (for spherical
        harmonics :math:`d_\ell = 2\ell+1`). **Measure-free**, symmetric with
        :meth:`reconstruct` — the quadrature weights are NOT baked in (unlike
        :meth:`analyze_transpose`, whose forward :meth:`analyze` carries them). The
        metric-aware ``_AdjointOperator`` combines this with the codomain/domain
        Gram to form the W-weighted Hilbert adjoint :math:`R^*`, so the
        :class:`FrameBase`'s reconstruction face gets ``.H`` for free — symmetric with
        the analysis face.
        """
        ...

    # ── The Gram structure (the projection-validity declaration) ──────────
    @property
    def gram_structure(self) -> GramStructure:
        r"""How this basis's frame-Gram is structured when it is the **trial** side.

        Declares whether :meth:`~orpheus.numerics.frame.FrameBase.project`'s row-sum
        probe is a valid coefficient normalisation for this basis (:class:`GramStructure`).
        Defaults to :attr:`GramStructure.DENSE` — the safe refusal: a new basis must
        consciously override this to :attr:`~GramStructure.DIAGONAL` or
        :attr:`~GramStructure.PARTITION_OF_UNITY` (having verified its Gram) before
        ``project`` will use the shortcut; a DENSE trial is refused rather than
        silently mis-projected.
        """
        return GramStructure.DENSE

    # ── The discrete Gram + the coefficient space ─────────────────────────
    @abstractmethod
    def mass_matrix(self, measure: "DiscreteMeasure", /) -> NDArray:
        r"""Discrete Gram :math:`\sum_n w_n\, \phi_j(x_n)\, \phi_k(x_n)` over ``measure``.

        The frame operator :math:`S = T^* T` in discrete form. Equals the
        continuous Gram (the basis's intrinsic metric) when the quadrature is
        exact to the basis's degree; the residual is a quadrature-exactness
        diagnostic. (A one-off diagnostic, so it is naturally measure-based and
        evaluates its own table.)
        """
        ...

    @property
    @abstractmethod
    def space(self) -> "FunctionSpace":
        r"""The coefficient :class:`FunctionSpace` this basis spans.

        Carries the basis's Gram as its ``inner_product_weights`` — the metric the
        :class:`FrameBase`'s codomain (and the Hilbert-adjoint machinery) reads. The
        basis owns exactly one space (the nodal/domain space comes from the
        measure), so the unqualified name is unambiguous — matching the
        ``Field.space`` convention. The :class:`FrameBase` re-exposes it provenance-
        qualified as ``frame.basis_space`` (beside ``frame.measure_space``).
        """
        ...
