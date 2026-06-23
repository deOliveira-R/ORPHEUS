r"""The :class:`Frame` — a discrete frame binding a :class:`Basis` to a measure.

A *frame* (harmonic-analysis frame theory) is a ``(basis, measure)`` pairing that
emits the operational analysis/reconstruction pair the rest of the algebra consumes:

* **analysis** :math:`M = T` — sampled values → coefficients (``frame.analysis``);
* **reconstruction** :math:`R` — coefficients → values, the canonical-dual synthesis
  (``frame.reconstruction``).

The naked synthesis :math:`T^* = S_0` is the shared :class:`Basis` primitive
(:meth:`~orpheus.numerics.basis.base.Basis.synthesize`) the weighted faces are each
one diagonal away from; it lives one level below the faces.

The pairing fixes BOTH spaces and the Galerkin discipline (see §"basis ⋈ measure"):

* the **basis** fixes the codomain — ``frame.basis_space`` (``= basis.space``), the
  coefficient space + its Gram (the disciplineless *trial* side);
* the **measure** fixes the domain — ``frame.measure_space``, the nodal space + the
  quadrature weights (the *measured* / test side).

So ``coefficient_space`` is never a third parameter: it is derived from the basis. The
:class:`Frame` is **layout-agnostic** — the index layout (the :math:`(\ell, m)` axes for
spherical harmonics; a flat mode axis for a Legendre basis) lives entirely in the basis,
which owns the weighted contractions; the :class:`Frame` caches the table ONCE
(``frame.table``) and the faces delegate.

Iso vs non-iso (a capability, not a separate path)
==================================================

The :class:`Frame` is the single mechanism for ALL choice-dependent change-of-basis
(GitHub #263): whether the analysis/reconstruction are mutually inverse (``R∘M = I`` —
an invertible Vandermonde, e.g. nodal-DG; the analysis face would advertise
``CAP_SOLVE``) or band-limiting (``R∘M`` = projector ≠ I, ``N > (L+1)²`` — spherical
harmonics; a section/retraction) is a *capability of the given frame*, not a reason for
a second mechanism. The spherical-harmonic frame is the non-iso case.

References
----------

* Grand Report v3 §5.4 / §19 — bases and harmonic projection.
* Christensen, O. (2016). *An Introduction to Frames and Riesz Bases*, 2nd ed. —
  the analysis operator :math:`T`, synthesis operator :math:`T^*`, frame operator
  :math:`S = T^*T`, and canonical dual.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

from numpy.typing import NDArray

from orpheus.numerics.basis.base import Basis
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    LinearOperatorMixin,
)
from orpheus.numerics.space import FunctionSpace


__all__ = ["Frame"]


@dataclass(frozen=True)
class Frame:
    r"""A discrete frame: a :class:`Basis` bound to a :class:`DiscreteMeasure`.

    Parameters
    ----------
    basis : Basis
        The synthesis (trial) side — fixes the codomain (:attr:`basis_space`).
    measure : DiscreteMeasure
        The analysis (test) side — fixes the domain (:attr:`measure_space`) and
        the quadrature weights.
    """

    basis: Basis
    measure: DiscreteMeasure

    @cached_property
    def table(self) -> NDArray:
        r"""The basis tabulated at the measure's nodes — :math:`\Phi(\text{node}, \text{mode})`.

        Evaluated ONCE and cached; the per-apply hot path of the faces reads this
        rather than re-tabulating (the L16 perf-regression guard).
        """
        return self.basis.evaluate(self.measure.nodes)

    @cached_property
    def measure_space(self) -> FunctionSpace:
        r"""The domain — the measure's induced discrete-:math:`L^2` space.

        Read straight off :attr:`DiscreteMeasure.space` (per-node values with
        the quadrature weights as the metric): the measure OWNS its domain
        space, symmetric with the basis owning its codomain
        (:attr:`basis_space` ``= basis.space``) — neither is fabricated here.
        """
        return self.measure.space

    @cached_property
    def basis_space(self) -> FunctionSpace:
        r"""The codomain — the basis's coefficient space (``= basis.space``)."""
        return self.basis.space

    @cached_property
    def analysis(self) -> "_FrameAnalysis":
        r"""The analysis face :math:`M = T` (``measure_space → basis_space``)."""
        return _FrameAnalysis(self)

    @cached_property
    def reconstruction(self) -> "_FrameReconstruction":
        r"""The reconstruction face :math:`R` (``basis_space → measure_space``)."""
        return _FrameReconstruction(self)


@dataclass(frozen=True)
class _FrameAnalysis(LinearOperatorMixin):
    r"""The analysis face :math:`M = T`: ``measure_space → basis_space``.

    A frame-backed :class:`LinearOperator` view; the math lives on the basis
    (:meth:`Basis.analyze` / :meth:`Basis.analyze_transpose`). Carries the swapped
    spaces and ``CAP_APPLY_TRANSPOSE``, so the metric-aware ``_AdjointOperator``
    gives ``.H`` (the W-weighted Hilbert adjoint :math:`g_C \cdot S_0`) for free.
    """

    frame: Frame
    # Plain unannotated class attr (NOT a dataclass field, NOT ClassVar) —
    # overrides the mixin's instance-attr ``capabilities`` the same way the
    # leaves override ``block_role`` (see LinearOperatorMixin).
    capabilities = frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})

    @property
    def domain(self) -> FunctionSpace:
        return self.frame.measure_space

    @property
    def codomain(self) -> FunctionSpace:
        return self.frame.basis_space

    def apply(self, values: NDArray, /) -> NDArray:
        return self.frame.basis.analyze(
            values, self.frame.table, self.frame.measure.weights,
        )

    def apply_transpose(self, coefficients: NDArray, /) -> NDArray:
        return self.frame.basis.analyze_transpose(
            coefficients, self.frame.table, self.frame.measure.weights,
        )


@dataclass(frozen=True)
class _FrameReconstruction(LinearOperatorMixin):
    r"""The reconstruction face :math:`R`: ``basis_space → measure_space``.

    A frame-backed :class:`LinearOperator` view delegating to
    :meth:`Basis.reconstruct` (the canonical-dual synthesis) and its representation
    transpose :meth:`Basis.reconstruct_transpose`. Carries the swapped spaces and
    ``CAP_APPLY_TRANSPOSE``, so the metric-aware ``_AdjointOperator`` gives ``R.H``
    (the W-weighted Hilbert adjoint :math:`R^* = \frac{(2\ell+1)^2}{4\pi} M`) for
    free — symmetric with the analysis face.
    """

    frame: Frame
    # Plain class attr (see _FrameAnalysis) — both faces now advertise the same
    # surface: apply + apply_transpose, so .H falls out of _AdjointOperator.
    capabilities = frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})

    @property
    def domain(self) -> FunctionSpace:
        return self.frame.basis_space

    @property
    def codomain(self) -> FunctionSpace:
        return self.frame.measure_space

    def apply(self, coefficients: NDArray, /) -> NDArray:
        return self.frame.basis.reconstruct(coefficients, self.frame.table)

    def apply_transpose(self, values: NDArray, /) -> NDArray:
        return self.frame.basis.reconstruct_transpose(values, self.frame.table)
