r"""The :class:`FrameBase` hierarchy — a discrete frame binding a :class:`Basis`
to a measure, with the **projection discipline carried by the type**.

A *frame* (harmonic-analysis frame theory) is a ``(basis, measure)`` pairing that
emits the operational analysis/reconstruction pair the rest of the algebra consumes:

* **analysis** :math:`M = T` — sampled values → coefficients (``frame.analysis``);
* **reconstruction** :math:`R` — coefficients → values, the canonical-dual synthesis
  (``frame.reconstruction``).

The naked synthesis :math:`T^* = S_0` is the shared :class:`Basis` primitive
(:meth:`~orpheus.numerics.basis.base.Basis.synthesize`) the weighted faces are each
one diagonal away from; it lives one level below the faces.

The pairing fixes BOTH spaces:

* the **basis** is the synthesis (**trial**) side — fixes the codomain
  (``frame.basis_space`` ``= basis.space``), the coefficient space + its Gram;
* the **measure** fixes the domain — ``frame.measure_space``, the nodal space + the
  quadrature weights.

So ``coefficient_space`` is never a third parameter: it is derived from the basis. The
hierarchy is **layout-agnostic** — the index layout (the :math:`(\ell, m)` axes for
spherical harmonics; a flat cell axis for the indicator basis) lives entirely in the
basis, which owns the weighted contractions; the frame caches the table ONCE
(``frame.table``) and the faces delegate.

Discipline IS a type — the Petrov-Galerkin base / Galerkin specialisation
====================================================================

Every Galerkin-style discretisation factors as an analysis/reconstruction pair
:math:`(M, R)`. The **discipline** — whether the *test* functions equal the *trial*
functions — is a genuine *kind of object*, not a flag, so it is carried by the **type**
(GitHub #268; this REVERSES the earlier ``projection.py`` discipline-ABC design, where
the discipline was a marker mixin on the operator role):

.. code-block:: text

   FrameBase                 abstract; the discipline-FREE mechanics
   │                         (table, spaces, reconstruction do NOT depend on the test side)
   └─ PetrovGalerkinFrame    explicit TEST basis (test ≠ trial in general); M = ⟨test, ·⟩_W
      └─ GalerkinFrame       test IS trial — STRENGTHENS the promise to Π* = R (a
                             symmetric, here-diagonal Gram). The angular spherical-
                             harmonic projection is the canonical pure-Galerkin frame.

``GalerkinFrame`` *is-a* ``PetrovGalerkinFrame`` with ``test is trial`` — Liskov-correct,
strengthening (never weakening) the base promises. A genuine ``test ≠ trial`` instance —
flux-weighted spatial homogenisation, spectrum-weighted energy condensation — preserves
a **bilinear** functional :math:`\langle\varphi^*, \Sigma\varphi\rangle` and so cannot be
posed as a Galerkin projection without *folding the solution into the metric*; that fold
is legitimate only for forward-flux, reaction-rate-only reduction and breaks under the
eigenvalue-consistent (adjoint-weighted) homogenisation reactor physics requires. Hence
the test side is a first-class **basis** (the test *space*), NEVER a weight smuggled onto
the measure: **the measure carries the axis + the fixed** :math:`L^2` **metric, never the
discipline**; the solution-weighting (:math:`\varphi`, :math:`\varphi^*`) lives on the
test side = the frame TYPE.

The trial side owns reconstruction; the test side owns analysis
---------------------------------------------------------------

``reconstruction`` (:math:`R`) synthesises a fine field from coefficients — it is purely
**trial**-side (``basis.reconstruct``), identical across disciplines. ``analysis``
(:math:`M`) measures a field against the **test** functions — :math:`(M f)_k =
\sum_n w_n\, \chi_k(x_n)\, f(x_n)` for test functions :math:`\chi_k` — so it reads the
test basis tabulated at the nodes (``frame.test_table``). For a :class:`GalerkinFrame`
``test is trial``, ``test_table is table``, and the analysis is bit-identical to the
single-discipline frame this hierarchy replaced (the 0-ULP scattering-kernel canary).

Iso vs non-iso (a capability, not a separate path)
==================================================

The frame is the single mechanism for ALL choice-dependent change-of-basis (GitHub
#263): whether the analysis/reconstruction are mutually inverse (``R∘M = I`` — an
invertible Vandermonde, e.g. nodal-DG; the analysis face would advertise ``CAP_SOLVE``)
or band-limiting (``R∘M`` = projector ≠ I, ``N > (L+1)²`` — spherical harmonics; a
section/retraction) is a *capability of the given frame*, not a reason for a second
mechanism. The spherical-harmonic frame is the non-iso case.

References
----------

* Grand Report v3 §5.4 / §19 — bases and harmonic projection.
* Christensen, O. (2016). *An Introduction to Frames and Riesz Bases*, 2nd ed. —
  the analysis operator :math:`T`, synthesis operator :math:`T^*`, frame operator
  :math:`S = T^*T`, and canonical dual.
* Brenner, S. C. and Scott, L. R. (2008). *The Mathematical Theory of Finite Element
  Methods*, 3rd ed. Springer. §3.4 — Galerkin vs Petrov-Galerkin (test vs trial space).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, replace
from functools import cached_property

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.basis.base import Basis
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    LinearOperator,
    OperatorProduct,
)
from orpheus.numerics.projection import AnalysisOperator, ReconstructionOperator
from orpheus.numerics.space import FunctionSpace


__all__ = ["FrameBase", "GalerkinFrame", "PetrovGalerkinFrame"]


@dataclass(frozen=True)
class FrameBase(ABC):
    r"""A discrete frame: a :class:`Basis` bound to a :class:`DiscreteMeasure`.

    The abstract base carries the **discipline-free** mechanics — the trial table, the
    two spaces, the reconstruction face, and the analysis face's wiring. The single
    abstract hook is :attr:`test` (the test basis), which the discipline subclasses
    fix: :class:`GalerkinFrame` (``test is trial``) and :class:`PetrovGalerkinFrame`
    (an explicit, generally distinct, test basis).

    Parameters
    ----------
    basis : Basis
        The synthesis (trial) side — fixes the codomain (:attr:`basis_space`).
    measure : DiscreteMeasure
        Fixes the domain (:attr:`measure_space`) and the quadrature weights.
    """

    basis: Basis
    measure: DiscreteMeasure

    # ── the discipline hook ───────────────────────────────────────────────
    @property
    @abstractmethod
    def test(self) -> Basis:
        r"""The **test** basis — the analysis (measured) side of the frame.

        :class:`GalerkinFrame` returns the trial :attr:`basis` (``test is trial``);
        :class:`PetrovGalerkinFrame` returns its explicit test basis. The analysis
        face reads this basis tabulated at the nodes (:attr:`test_table`).
        """
        ...

    # ── trial side (reconstruction + the Galerkin analysis) ───────────────
    @cached_property
    def table(self) -> NDArray:
        r"""The TRIAL basis tabulated at the measure's nodes — :math:`\Phi(\text{node}, \text{mode})`.

        Evaluated ONCE and cached; the reconstruction face and (for a Galerkin frame)
        the analysis face read this rather than re-tabulating (the L16 perf guard).
        """
        return self.basis.evaluate(self.measure.nodes)

    @cached_property
    def basis_space(self) -> FunctionSpace:
        r"""The codomain — the trial basis's coefficient space (``= basis.space``)."""
        return self.basis.space

    # ── test side (the analysis face) ─────────────────────────────────────
    @cached_property
    def test_table(self) -> NDArray:
        r"""The TEST basis tabulated at the measure's nodes (the analysis contraction).

        Defaults to the test basis evaluated at the nodes; :class:`GalerkinFrame`
        overrides it to *reuse* :attr:`table` (``test is trial`` ⟹ the same array,
        no re-evaluation, 0-ULP-identical analysis).
        """
        return self.test.evaluate(self.measure.nodes)

    @cached_property
    def test_space(self) -> FunctionSpace:
        r"""The analysis codomain — the test basis's coefficient space (``= test.space``).

        :class:`GalerkinFrame` overrides it to *reuse* :attr:`basis_space` (``test is
        trial`` ⟹ the same cached space object, preserving the analysis-codomain ``is``
        identity of the single-discipline frame this hierarchy replaced).
        """
        return self.test.space

    # ── the measure (domain) side ─────────────────────────────────────────
    @cached_property
    def measure_space(self) -> FunctionSpace:
        r"""The domain — the measure's induced discrete-:math:`L^2` space.

        Read straight off :attr:`DiscreteMeasure.space` (per-node values with the
        quadrature weights as the metric): the measure OWNS its domain space,
        symmetric with the basis owning its codomain — neither is fabricated here.
        """
        return self.measure.space

    # ── the two operator faces ────────────────────────────────────────────
    @cached_property
    def analysis(self) -> "_FrameAnalysis":
        r"""The analysis face :math:`M = T` (``measure_space → test_space``)."""
        return _FrameAnalysis(self)

    @cached_property
    def reconstruction(self) -> "_FrameReconstruction":
        r"""The reconstruction face :math:`R` (``basis_space → measure_space``)."""
        return _FrameReconstruction(self)

    # ── composed operators (the "define Frame, compose, done" production path) ──
    def conjugate(self, operator: LinearOperator, /) -> LinearOperator:
        r"""Frame-conjugate a coefficient-space operator: :math:`R \circ A \circ M`.

        THE production composition for a method whose action is "project to
        coefficients, act there, reconstruct" — e.g. SN anisotropic scattering
        :math:`S_{\ell\ge 1} = R\,\Lambda\,M` (the per-ordinate redistribution). Returns
        the typed :class:`~orpheus.numerics.operator.OperatorProduct`
        ``OperatorProduct(R, OperatorProduct(A, M))``, whose
        :meth:`~orpheus.numerics.operator.OperatorProduct.apply` is ``R(A(M·x))`` — the
        SAME numpy chain a hand-rolled ``reconstruction.apply(A.apply(analysis.apply(x)))``
        runs, now ONE named operator (Cardinal Rule 2: the composition IS the production
        path, not a parallel "semantic" reading of it).

        ``operator`` must compose between the faces — its ``domain`` is the analysis
        codomain (:attr:`test_space`) and its ``codomain`` the reconstruction domain
        (:attr:`basis_space`); the :class:`OperatorProduct` space-compatibility guard
        enforces it (an endomorphism on the coefficient space when ``test == trial``).
        """
        return OperatorProduct(
            self.reconstruction, OperatorProduct(operator, self.analysis),
        )

    def reconstruct_after(self, operator: LinearOperator, /) -> LinearOperator:
        r"""Reconstruct after a coefficient-space operator: :math:`R \circ A`.

        The :meth:`conjugate` sub-operator for inputs ALREADY in coefficient space (the
        analysis :math:`M` already applied) — e.g. the angular-windowed SN moment
        iterate, whose bulk IS :math:`M\psi`, so only :math:`R\,\Lambda` remains.
        Returns ``OperatorProduct(R, A)`` (``apply`` = ``R(A·c)``). Wiring a windowed
        consumer to :meth:`conjugate` instead would erroneously re-apply :math:`M` (a
        double-projection).
        """
        return OperatorProduct(self.reconstruction, operator)

    # ── the coefficient-extraction verb (homogenise / condense) ──────────
    @cached_property
    def gram(self) -> FunctionSpace:
        r"""The coefficient space carrying the diagonal frame Gram :math:`G_R = \langle\chi_R, \phi_R\rangle_W`.

        The cross Gram :math:`G_{kj} = \langle\chi_k, \phi_j\rangle_W = (M R)_{kj}` of
        the analysis :math:`M` and reconstruction :math:`R`. For a disjoint-support
        trial (cell / group indicators) it is **diagonal**, so the full operator
        collapses to its diagonal :math:`G_R = (M R\,\mathbf 1)_R = \int_R w\,
        \mathrm{d}V` (:math:`w` the test weight) — a single
        ``analysis ∘ reconstruction`` probe of the all-ones coefficient vector
        (off-diagonals are structurally zero, so the row sum IS the diagonal). The
        diagonal acquires the test weight's trailing (group, …) shape from the
        analysis face, and is installed as the coefficient space's metric so
        :meth:`project`'s normalisation is the reciprocal
        :meth:`~orpheus.numerics.space.FunctionSpace.apply_inverse_metric` (whose
        Moore–Penrose pseudo-inverse zeroes empty / zero-weight regions for free).
        """
        ones = np.ones(self.basis_space.shape)
        diagonal = self.analysis.apply(self.reconstruction.apply(ones))
        return replace(self.test_space, inner_product_weights=diagonal)

    def project(self, field: NDArray, /) -> NDArray:
        r"""Extract coefficients: :math:`G^{-1} M f` — the homogenise / condense verb.

        The Petrov-Galerkin coefficient extraction: analyse the field against the
        test functions (:math:`(M f)_k = \langle\chi_k, f\rangle_W`), then normalise
        by the cross :attr:`gram` :math:`G`. For flux-weighted spatial homogenisation
        (``test`` :math:`= \varphi\cdot\mathbf 1_R`, ``trial`` :math:`= \mathbf 1_R`)
        this is the rate-preserving effective cross section :math:`\Sigma_R = \int_R
        \varphi\Sigma\,\mathrm{d}V / \int_R\varphi\,\mathrm{d}V`; for a
        :class:`GalerkinFrame` (``test = trial``) it is the orthogonal projection onto
        the coarse space. The Gram is diagonal for every real consumer (disjoint
        indicators / orthogonal harmonics), so the normalisation is a reciprocal —
        the dense solve is the (unbuilt) least-squares seam only. Trailing (group, …)
        axes ride the analysis and broadcast against the diagonal Gram (so a vector
        channel divides by :math:`\Phi_{R,g}` and a ``[g_from, g_to]`` matrix channel
        by its source-group :math:`\Phi_{R,g_{\mathrm{from}}}`).
        """
        return self.gram.apply_inverse_metric(self.analysis.apply(field))


@dataclass(frozen=True)
class PetrovGalerkinFrame(FrameBase):
    r"""A frame with an explicit **test** basis (the Petrov-Galerkin discipline).

    The analysis measures against test functions :math:`\chi_k` that need NOT equal the
    trial functions :math:`\phi_k`, so the coefficient extraction :math:`G^{-1} M` (the
    ``project`` verb, built in a later phase) uses the *cross* Gram
    :math:`G_{kj} = \langle \chi_k, \phi_j\rangle` and the Hilbert adjoint
    :math:`M^* \ne R`.

    Flux-weighted spatial homogenisation and spectrum-weighted energy condensation are
    the headline consumers: the test basis is the trial cell/group indicator weighted
    by the solution (:math:`\varphi\cdot\mathbf 1_R`, or the adjoint-weighted
    :math:`\varphi^*\cdot\mathbf 1_R` for eigenvalue-consistent homogenisation) — a
    genuinely different basis, NOT a metric on the measure.

    The test basis is **required** — a Petrov-Galerkin frame is *defined* by carrying an
    explicit test space, so there is no implicit "``None`` means trial" default. Passing
    the trial :attr:`basis` itself as the test basis is the legal **degenerate**: the
    frame then behaves like a :class:`GalerkinFrame` (a useful "PG reduces to Galerkin"
    equivalence test) but does NOT advertise the strengthened :math:`M^* = R` promise.
    To get ``test = trial`` *without* naming a test basis, construct a
    :class:`GalerkinFrame` instead.

    Parameters
    ----------
    test_basis : Basis
        The explicit test basis. Generally distinct from the trial ``basis``; passing
        ``basis`` itself is the Galerkin degenerate.
    """

    test_basis: Basis

    @property
    def test(self) -> Basis:
        return self.test_basis


@dataclass(frozen=True, init=False)
class GalerkinFrame(PetrovGalerkinFrame):
    r"""The Galerkin specialisation: the test basis IS the trial basis.

    A :class:`GalerkinFrame` is constructed from ``(basis, measure)`` alone — it binds
    the inherited (required) :attr:`test_basis` to the trial ``basis``, so
    ``test is trial``. That coincidence strengthens the base Petrov-Galerkin promises:
    the Gram is symmetric (here diagonal — SH-orthogonal / disjoint indicators), so
    :math:`M^* = R` up to the basis's dual factor and the coefficient extraction is a
    reciprocal, not a solve. The angular spherical-harmonic projection
    (``quadrature.angular_frame(L)``) is the canonical pure-Galerkin frame; the in-sweep
    moment accumulation AND the §5.6 scattering kernel share THIS object.

    It is a genuine subtype of :class:`PetrovGalerkinFrame` (Liskov: a Galerkin frame
    IS a Petrov-Galerkin frame whose test basis equals its trial basis). The constructor
    takes **no** ``test_basis`` argument — that ``test ≠ trial`` freedom is exactly what
    distinguishes a :class:`PetrovGalerkinFrame`, so a distinct test basis is forbidden
    here by the constructor signature itself (it is spelled by building one of those).
    """

    def __init__(self, basis: Basis, measure: DiscreteMeasure) -> None:
        # test = trial: bind the (required) test_basis to the trial basis. A frozen
        # dataclass forbids attribute assignment, so set the three fields directly.
        object.__setattr__(self, "basis", basis)
        object.__setattr__(self, "measure", measure)
        object.__setattr__(self, "test_basis", basis)

    @cached_property
    def test_table(self) -> NDArray:
        # test is trial → reuse the trial table (same array; 0-ULP-identical analysis).
        return self.table

    @cached_property
    def test_space(self) -> FunctionSpace:
        # test is trial → reuse the cached trial space (preserves the codomain `is`).
        return self.basis_space


@dataclass(frozen=True)
class _FrameAnalysis(AnalysisOperator):
    r"""The analysis face :math:`M = T`: ``measure_space → test_space``.

    A frame-backed :class:`AnalysisOperator` view; the math lives on the TEST basis
    (:meth:`Basis.analyze` / :meth:`Basis.analyze_transpose`) tabulated at the nodes
    (``frame.test_table``). Carries the swapped spaces and ``CAP_APPLY_TRANSPOSE``, so
    the metric-aware ``_AdjointOperator`` gives ``.H`` (the W-weighted Hilbert adjoint)
    for free. For a :class:`GalerkinFrame` the test basis is the trial basis, so this is
    the :math:`Y^* W` projection bit-identical to the single-discipline frame.
    """

    frame: FrameBase
    # Plain unannotated class attr (NOT a dataclass field, NOT ClassVar) —
    # overrides the role ABC's annotated ``capabilities`` the same way the
    # leaves override ``block_role`` (see LinearOperatorMixin).
    capabilities = frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})

    @property
    def domain(self) -> FunctionSpace:
        return self.frame.measure_space

    @property
    def codomain(self) -> FunctionSpace:
        return self.frame.test_space

    def apply(self, values: NDArray, /) -> NDArray:
        return self.frame.test.analyze(
            values, self.frame.test_table, self.frame.measure.weights,
        )

    def apply_transpose(self, coefficients: NDArray, /) -> NDArray:
        return self.frame.test.analyze_transpose(
            coefficients, self.frame.test_table, self.frame.measure.weights,
        )


@dataclass(frozen=True)
class _FrameReconstruction(ReconstructionOperator):
    r"""The reconstruction face :math:`R`: ``basis_space → measure_space``.

    A frame-backed :class:`ReconstructionOperator` view delegating to the TRIAL basis's
    :meth:`Basis.reconstruct` (the canonical-dual synthesis) and its representation
    transpose :meth:`Basis.reconstruct_transpose` — reconstruction is purely trial-side,
    identical across disciplines. Carries the swapped spaces and ``CAP_APPLY_TRANSPOSE``,
    so the metric-aware ``_AdjointOperator`` gives ``R.H`` for free — symmetric with the
    analysis face.
    """

    frame: FrameBase
    # Plain class attr (see _FrameAnalysis) — both faces advertise apply +
    # apply_transpose, so .H falls out of _AdjointOperator.
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
