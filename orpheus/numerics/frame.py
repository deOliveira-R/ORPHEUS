r"""The :class:`FrameBase` hierarchy — a discrete frame binding a :class:`Basis`
to a measure, with the **projection discipline carried by the type**.

A *frame* (harmonic-analysis frame theory) is a ``(basis, measure)`` pairing that
emits the operational analysis/reconstruction pair the rest of the algebra consumes:

* **analysis** :math:`M = T` — sampled values → coefficients
  (:meth:`FrameAnalysis.of <orpheus.numerics.frame_faces.FrameAnalysis.of>`);
* **reconstruction** :math:`R` — coefficients → values, the canonical-dual synthesis
  (:meth:`FrameReconstruction.of <orpheus.numerics.frame_faces.FrameReconstruction.of>`).

The frame holds only its **data** — the binding, its sampling table and its Gram. The
faces are built FROM that data once spaces and operators exist, so they live
downstream, in :mod:`orpheus.numerics.frame_faces` (concept-order rulings 22 and 34),
as constructors of the frame and never as attributes of it.

The naked synthesis :math:`T^* = S_0` is the shared :class:`Basis` primitive
(:meth:`~orpheus.numerics.basis.base.Basis.synthesize`) the weighted faces are each
one diagonal away from; it lives one level below the faces.

The pairing fixes BOTH spaces:

* the **basis** is the synthesis (**trial**) side — fixes the codomain
  (:func:`~orpheus.numerics.frame_faces.basis_space` ``= basis.space`` dressed with
  the frame's **Parseval metric**, the inverse of the discrete trial Gram — a
  property of the *pairing*, not of the basis alone);
* the **measure** fixes the domain — :func:`~orpheus.numerics.frame_faces.measure_space`,
  the nodal space + the quadrature weights.

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
   │                         (the table and the Gram do NOT depend on the test side)
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
invertible Vandermonde, e.g. nodal-DG; the analysis face would become invertible)
or band-limiting (``R∘M`` = projector ≠ I, ``N > (L+1)²`` — spherical harmonics; a
section/retraction) is a *capability of the given frame*, not a reason for a second
mechanism. The spherical-harmonic frame is the non-iso case.

References
----------

* Grand Report v3 §5.4 / §19 — bases and harmonic projection.
* Christensen, O. (2016). *An Introduction to Frames and Riesz Bases*, 2nd ed. —
  the analysis and synthesis operators, the frame operator, and the canonical dual
  (Christensen writes :math:`T` for synthesis and :math:`S = TT^*`).
* Casazza, P. G. and Lynch, R. G. (2016). A brief introduction to Hilbert space
  frame theory and its applications. *Proc. Sympos. Appl. Math.* 73, 1–51, §4.1 —
  the letters used here: analysis operator :math:`T`, synthesis operator
  :math:`T^*`, frame operator :math:`S = T^*T` on the sampled space, Gram
  operator :math:`TT^*` on the coefficients (doi:10.1090/psapm/073/00627).
* Brenner, S. C. and Scott, L. R. (2008). *The Mathematical Theory of Finite Element
  Methods*, 3rd ed. Springer. §3.4 — Galerkin vs Petrov-Galerkin (test vs trial space).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from functools import cached_property
from typing import NamedTuple

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.axis import BasisKind, EnergyAxis
from orpheus.numerics.basis.base import Basis, GramStructure, TruncatedBasis
from orpheus.numerics.basis.indicator_basis import IndicatorBasis
from orpheus.numerics.manifold import IndexSet, ManifoldMap, quotient_onto
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.operator import AxisRetractionOperator, AxisSectionOperator
from orpheus.numerics.space import FunctionSpace


__all__ = ["FrameBase", "GalerkinFrame", "PetrovGalerkinFrame"]


#: Relative off-diagonal threshold for the :attr:`FrameBase.discrete_gram_structure`
#: diagonality verdict. The measured separation is wide open: an SH-degree-exact
#: sphere cubature's off-diagonals sit at ~1e-16 of the Cauchy–Schwarz scale
#: :math:`\sqrt{G_{jj}G_{kk}}`, while a dense frame's live Gram reaches 1.000
#: of it (`[M]` 2026-10-18, ``frame_faces.angular_frame(folded_product(4, 4), 2)``; the
#: pre-ERR-080 slab GL Gram sat at 0.93, `[M]` 2026-08-23) — any threshold in
#: (1e-12, 1e-2) draws the same verdict; 1e-10 leaves two orders of headroom
#: for accumulated roundoff at high mode counts.
_DISCRETE_GRAM_DIAGONALITY_RTOL = 1e-10


def _descent_arrow(basis: Basis, measure: DiscreteMeasure, role: str) -> ManifoldMap:
    r"""G0 — the frame's two halves must name ONE orbit space, up to a quotient map.

    A frame binds functions on ``basis.domain`` to a rule on
    ``measure.support``; that is well-posed iff the functions can be
    evaluated at the rule's nodes, i.e. iff a quotient map
    ``measure.support -> basis.domain`` EXISTS
    (:func:`~orpheus.numerics.manifold.quotient_onto`: equality, the entry's
    own map, or the induced :math:`M/K \to M/H` for :math:`K \subseteq H`).
    The frame's table is the basis pulled back along that arrow. ONE
    predicate (user-ruled 2026-09-02, #429 tracker 2.2): it refuses the Part
    I bug — the full harmonics (``Trivial``) on the slab's :math:`S^2/O(2)_x`
    — and admits the slab's Legendre basis, the fold's σ-even harmonics, the
    full harmonics on a full-sphere rule, AND the Legendre basis on a
    full-sphere rule. `[M]` 2026-09-02, before this gate, every one of the
    21 frames three real solves built passed it; after the forgery's
    retirement the slab-harmonics pairing is exactly the one it refuses.
    """
    arrow = quotient_onto(measure.support, basis.domain)
    if arrow is None:
        spent = getattr(measure, "quotient_group", None)
        has = getattr(basis, "invariance_group", None)
        raise ValueError(
            f"a frame binds {role} functions on {basis.domain.name!r} to a "
            f"rule on {measure.support.name!r}, and no quotient map "
            f"{measure.support.name} -> {basis.domain.name} exists (the rule "
            f"spent {spent.name if spent is not None else None}, the basis "
            f"has {has.name if has is not None else None}). A basis eats "
            f"points of its own orbit space or of a finer one; give the rule "
            f"the basis its orbit space carries "
            f"(frame_faces.angular_frame derives it) — #429, ERR-080."
        )
    return arrow


@dataclass(frozen=True)
class FrameBase(ABC):
    r"""A discrete frame: a :class:`Basis` bound to a :class:`DiscreteMeasure`.

    The abstract base carries the **discipline-free** data — the trial table, the
    discrete Gram and its measured structure, and the test table. The faces built from
    it (the spaces, the analysis and reconstruction operators) are
    :mod:`orpheus.numerics.frame_faces`. The single
    abstract hook is :attr:`test` (the test basis), which the discipline subclasses
    fix: :class:`GalerkinFrame` (``test is trial``) and :class:`PetrovGalerkinFrame`
    (an explicit, generally distinct, test basis).

    Parameters
    ----------
    basis : Basis
        The synthesis (trial) side — fixes the codomain
        (:func:`~orpheus.numerics.frame_faces.basis_space`).
    measure : DiscreteMeasure
        Fixes the domain (:func:`~orpheus.numerics.frame_faces.measure_space`) and the
        quadrature weights.
    """

    basis: Basis
    measure: DiscreteMeasure

    def __post_init__(self) -> None:
        # G0 at construction, on BOTH halves (the trial here; the test basis
        # through :attr:`test_descent` on first use, since the discipline
        # subclass binds it) — a frame over an incompatible pairing is
        # unspellable, loudly, early.
        _descent_arrow(self.basis, self.measure, "trial")

    # ── the G0 arrows: what the tables are pulled back along ──────────────
    @cached_property
    def descent(self) -> ManifoldMap:
        r"""The quotient map ``measure.support -> basis.domain`` the trial table is pulled back along (G0)."""
        return _descent_arrow(self.basis, self.measure, "trial")

    @cached_property
    def test_descent(self) -> ManifoldMap:
        r"""The same arrow for the TEST basis (``test is trial`` on a Galerkin frame)."""
        if self.test is self.basis:
            return self.descent
        return _descent_arrow(self.test, self.measure, "test")

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
        The nodes are pulled back along the G0 arrow first (:attr:`descent`):
        the identity for a basis on the rule's own orbit space, the quotient
        map for a Legendre basis on a full-sphere rule.
        """
        return self.basis.evaluate(self.descent(self.measure.nodes))

    @cached_property
    def discrete_gram(self) -> NDArray:
        r"""The TRIAL basis's discrete Gram over the measure — flattened ``(K, K)``.

        :math:`G_{jk} = \sum_n w_n\,\phi_j(x_n)\,\phi_k(x_n)`, computed from the
        CACHED :attr:`table` (the L16 perf guard — never re-tabulates), with the
        basis's mode layout flattened to one axis (``K =`` the coefficient count,
        row-major — the same order as ``table.reshape(N, K)``). The one home
        of the mass matrix: a basis alone does not mint it (it pairs the basis
        with a measure, which is this frame). The Parseval-metric source
        (:func:`~orpheus.numerics.frame_faces.basis_space`) and the diagonality verdict
        (:attr:`discrete_gram_structure`) both read it, and a gate that wants
        the Gram in the basis's own mode layout reshapes this one.
        O(N·K²), once per frame.
        """
        weights = self.measure.weights
        flat = self.table.reshape(weights.shape[0], -1)
        return np.einsum("n,nj,nk->jk", weights, flat, flat)

    @cached_property
    def discrete_gram_structure(self) -> GramStructure:
        r"""The MEASURED diagonality of :attr:`discrete_gram` — DIAGONAL or DENSE.

        The measured counterpart of the trial's DECLARED
        :attr:`~orpheus.numerics.basis.base.Basis.gram_structure`, and a genuinely
        different fact: the declaration states which structure the *cross* Gram
        :math:`MR` has for :func:`~orpheus.numerics.frame_faces.project`'s row-sum probe; this verdict states
        whether the *trial* Gram :math:`G` is diagonal **on this measure** — the
        precondition for the Parseval metric
        (:func:`~orpheus.numerics.frame_faces.basis_space`). The two can
        disagree: the SH basis declares DIAGONAL (continuum-orthogonal) yet
        measures DENSE on a rule that does not integrate its products exactly
        (`[M]` 2026-10-18 on ``frame_faces.angular_frame(product(4, 4), 2)``: live
        off-diagonals up to 0.293 of the Cauchy–Schwarz scale; on
        ``frame_faces.angular_frame(folded_product(4, 4), 2)``, whose mirror-even basis
        also declares DIAGONAL, 1.000; ``scratch/metric_docstring/probe_flagship.py``.
        Until the ERR-080 repair the example was the slab GL measure, `[M]`
        2026-08-23, ``scratch/probe_f1_parseval_slab.py``),
        and an :class:`~orpheus.numerics.basis.overlap_basis.OverlapBasis`
        declares PARTITION_OF_UNITY while its trial Gram measures DENSE.

        DIAGONAL iff no diagonal entry is negative (a negative-weight quadrature
        can make :math:`G` indefinite, and an indefinite form is not a metric)
        AND every live off-diagonal is below
        ``_DISCRETE_GRAM_DIAGONALITY_RTOL`` of the Cauchy–Schwarz scale
        :math:`\sqrt{G_{jj}G_{kk}}`. Structurally dead slots (:math:`G_{kk}=0` —
        layout padding, σ-odd folded columns, empty indicator regions) are
        exempt, but any coupling INTO a dead slot is DENSE.
        """
        gram = self.discrete_gram
        diag = np.diagonal(gram)
        if np.any(diag < 0.0):
            return GramStructure.DENSE
        off = np.abs(gram - np.diag(diag))
        scale = np.sqrt(np.outer(diag, diag))
        live = scale > 0.0
        if np.any(off[~live] > 0.0):
            return GramStructure.DENSE
        if np.any(live) and (
            float(np.max(off[live] / scale[live]))
            > _DISCRETE_GRAM_DIAGONALITY_RTOL
        ):
            return GramStructure.DENSE
        return GramStructure.DIAGONAL

    # ── test side (the analysis face) ─────────────────────────────────────
    @cached_property
    def test_table(self) -> NDArray:
        r"""The TEST basis tabulated at the measure's nodes (the analysis contraction).

        Defaults to the test basis evaluated at the nodes; :class:`GalerkinFrame`
        overrides it to *reuse* :attr:`table` (``test is trial`` ⟹ the same array,
        no re-evaluation, 0-ULP-identical analysis).
        """
        return self.test.evaluate(self.test_descent(self.measure.nodes))


@dataclass(frozen=True)
class PetrovGalerkinFrame(FrameBase):
    r"""A frame with an explicit **test** basis (the Petrov-Galerkin discipline).

    The analysis measures against test functions :math:`\chi_k` that need NOT equal the
    trial functions :math:`\phi_k`, so the coefficient extraction :math:`G^{-1} M` (the
    :func:`~orpheus.numerics.frame_faces.project` verb) uses the *cross* Gram
    :math:`G_{kj} = \langle \chi_k, \phi_j\rangle` and the Hilbert adjoint
    :math:`M^* \ne R`.

    Flux-weighted spatial homogenisation and spectrum-weighted energy condensation are
    the headline consumers: the test basis is the trial cell/group indicator weighted
    by the solution (:math:`\varphi\cdot\mathbf 1_R` forward, or the bilinear PAIR
    weight :math:`\varphi^*\!\odot\varphi\cdot\mathbf 1_R` for the eigenvalue-consistent
    adjoint-weighted collapse — the product, NOT a bare :math:`\varphi^*` swap; theorem
    T1 of :mod:`orpheus.derivations.common.homogenization`) — a genuinely different
    basis, NOT a metric on the measure.

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
    (``frame_faces.angular_frame(quadrature, L)``) is the canonical pure-Galerkin frame; the in-sweep
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
        _descent_arrow(basis, measure, "trial")   # G0 (the dataclass path runs __post_init__; this ctor must too)

    @cached_property
    def test_table(self) -> NDArray:
        # test is trial → reuse the trial table (same array; 0-ULP-identical analysis).
        return self.table

    def at_order(self, L_new: int, /) -> "GalerkinFrame":
        r"""THIS frame over the same measure with its truncated trial family
        cut at ``L_new`` — the verb a frame-dressed moment head truncates
        THROUGH (:func:`~orpheus.numerics.spaces.moment_head.truncated_head`,
        CS4c step 6 item 6.2c-ii): the head's metric is the frame's to
        install at every order (the discrete Gram's verdict can flip with
        :math:`L`), so a lower-order head is a lower-order FRAME's dressing,
        never a slice of the parent's. The same class (a
        :class:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame`
        re-poses as a ``HarmonicFrame``), the same measure; structurally
        equal to the quadrature's own interned frame at that order.
        """
        basis = self.basis
        if not isinstance(basis, TruncatedBasis):
            raise TypeError(
                f"{type(self).__name__}.at_order: the trial "
                f"{type(basis).__name__} carries no truncation order, so "
                f"there is no family to cut at L_new={L_new}."
            )
        lower = basis.at_order(L_new)
        assert isinstance(lower, Basis)  # narrowing only: the protocol's conformers are bases
        return type(self)(lower, self.measure)


# ── The axis collapse pair — the rank-one frame's induced output (S6.0b) ──


class _AxisCollapsePair(NamedTuple):
    """One axis's minted collapse pair (``section`` is ``None`` iff Σw = 0)."""

    retraction: AxisRetractionOperator
    section: AxisSectionOperator | None


def _collapse_pair(space: FunctionSpace, axis_label: str) -> _AxisCollapsePair:
    r"""Mint the axis collapse pair — the single-region indicator frame's output.

    THE stage-2 generator discipline, applied at rank one (user, ruled
    2026-08-24): *"A stage-2 generator induces structure on both the space
    and the operator, and the two inductions must be minted together, at
    one site. … Forgetting = retaining the induced parts; accessors are
    provenance."* Here the generator is the frame over the axis's index
    set — the single-region :class:`IndicatorBasis` covering synthetic
    index nodes :math:`\{0, \ldots, n-1\}`, bound to the axis's measure
    in a :class:`GalerkinFrame` — built eagerly HERE, read for its induced
    data, and DISCARDED (the forgetful map: neither the axis nor the
    operators keep the generator, and a frame FACE is a view holding
    ``frame:``, so no face is retained either). What is read off it:

    * the kernel weights — the frame measure's diagonal (the analysis
      face's content IS the retraction's contraction);
    * the section's divisor — the 1×1 :attr:`FrameBase.discrete_gram`
      entry :math:`\Sigma w`, i.e. the rank-one **Parseval metric**
      (F-0's inverse-discrete-Gram theorem at :math:`K = 1`):
      :math:`E = R_{\text{frame}} \circ G^{-1}`. A zero entry (a signed
      measure summing to zero) means the Gram is singular — the frame has
      no canonical dual, so no section exists: that arm is ``None`` and
      :meth:`FunctionSpace.section` refuses at access.

    The operators are constructed together at this one site (the
    two-inductions clause; the retraction's adjoint, the pullback, is
    minted inside the retraction's constructor, so all three arrows come
    from this mint); the tightness gate
    (``tests/gates/numerics/test_axis_marginal.py``) pins the minted kernels
    against the literal frame's face contents, and the gram-derivation
    gate pins the divisor against :attr:`FrameBase.discrete_gram`.

    Admission (typed refusals, all at this mint — the public path is
    :meth:`FunctionSpace.retraction` / :meth:`FunctionSpace.section`,
    which memoize this mint per axis label):

    * the space must be **axis-built** (``axes is not None``) — an
      axes-less space (a hand-named legacy space, or a product with an
      axes-less head) has no named factors to marginalise over;
    * ``axis_label`` must name **exactly one** axis;
    * the axis must be **NODAL**: a "marginal" over a MODAL axis would
      contract expansion COEFFICIENTS with the basis mass, which is not
      an integral of the represented function — the modal average is the
      coefficient at the average slot, not a weighted sum;
    * the axis must not be a typed :class:`~orpheus.numerics.axis.EnergyAxis`
      — collapse doctrine clause 2 (partition-integration of an
      :math:`L^1` class): the energy axis PERSISTS at its one-cell
      member, because the one-group limit keeps its edges and spectrum
      (:math:`\langle\bar\sigma, \phi\rangle` consumes the
      partition). The clause-2 collapse already ships as condensation
      (``EnergyGrid.overlap_to``, the Petrov-Galerkin condensation
      frames), and a drop-form marginal here would twin it (Cardinal
      Rule 2). Untyped generic axes stay admitted whatever their label —
      the clause gate reads the TYPE, never the label string (stringly
      dispatch rejected); it becomes fully structural axis-family
      polymorphism when CS2's typed axes land;
    * a single-axis space is refused — its marginal would be a bare
      scalar, which is not a :class:`FunctionSpace` (contract with the
      space's inner product instead);
    * the axis must carry no positioned form (a dense Gram the space's
      overlay places on the axis's block): the retraction integrates with
      the axis's measure, so a second measure on the same block leaves
      fibre integration undefined.

    The marginal keeps the positioned forms of the axes that remain, so
    the full space's metric is exactly the product (the collapsed axis's
    measure) ⊗ (the marginal's metric). That product condition is what
    makes the pullback :class:`~orpheus.numerics.operator.AxisPullbackOperator`
    the Hilbert adjoint of the retraction, with no metric sandwich.
    """
    # Arms 1-2 (axis-built; exactly-one label) discharge in the shared
    # resolver FunctionSpace._axis_index (un-weld arc S-1) — one home for
    # the by-label refusal vocabulary.
    k = space._axis_index(axis_label)
    axes = space.axes
    assert axes is not None  # narrowing only: _axis_index refused the None case
    axis = axes[k]
    if axis.kind is BasisKind.MODAL:
        raise TypeError(
            f"collapse pair: axis {axis_label!r} is MODAL — contracting "
            f"expansion coefficients with the basis mass is not an "
            f"integral of the represented function. The modal average is "
            f"the coefficient at the average slot (slice it), not a "
            f"weighted sum."
        )
    if isinstance(axis, EnergyAxis):
        raise TypeError(
            f"collapse pair: axis {axis_label!r} is a typed EnergyAxis, "
            f"which PERSISTS at its one-cell member (collapse doctrine "
            f"clause 2 — the one-group limit keeps its edges and "
            f"spectrum, because ⟨σ̄,φ⟩ consumes the partition). The "
            f"energy collapse is condensation: use EnergyGrid.overlap_to "
            f"/ the Petrov-Galerkin condensation frame, not a drop-form "
            f"marginal."
        )
    if len(axes) == 1:
        raise ValueError(
            f"collapse pair: {space!r} has only the {axis_label!r} axis "
            f"— its marginal would be a bare scalar, which is not a "
            f"FunctionSpace. Contract with the space's inner product "
            f"instead."
        )

    # A positioned form on the collapsed axis (a dense Gram on a
    # measure-less axis's block, the one form no axis measure can spell)
    # would put a second measure on the block the retraction integrates.
    form = space._overlay_forms()[k]
    if form is not None:
        raise ValueError(
            f"collapse pair: axis {axis_label!r} carries a positioned "
            f"{type(form).__name__} on {space!r} — the retraction "
            f"integrates with the axis's own measure while the space's "
            f"metric on that block is the form, two measures on one block, "
            f"so fibre integration over the axis is not defined."
        )

    # The ndarray dims this axis occupies: axes map to dims by
    # cumulative rank (an axis's shape may span several dims).
    start = sum(len(ax.shape) for ax in axes[:k])
    dims = tuple(range(start, start + len(axis.shape)))
    # The marginal keeps the remaining axes' forms, so the full metric is
    # exactly (the axis measure) ⊗ (the marginal metric): the product
    # condition under which the retraction's adjoint is the pullback.
    marginal_space = space.without_axis(axis_label)

    # The generator, eagerly: the axis's measure on synthetic index
    # nodes, under the single-region indicator covering them all.
    w = axis.weights
    flat_weights = (
        np.ones(int(np.prod(axis.shape)))
        if w is None
        else np.asarray(w, dtype=float).ravel()
    )
    n = int(flat_weights.shape[0])
    # ⭐ One frame, one manifold — and now literally one expression. Both
    # halves used to say it separately (the basis as an ``IndexSet``, the
    # measure as ``f"index({axis_label})"``), which is two spellings of one
    # fact and exactly what a frame's two halves must not have.
    points = IndexSet(label=axis_label, n=n)
    frame = GalerkinFrame(
        basis=IndicatorBasis(
            edges_per_axis=(np.array([-0.5, n - 0.5]),),
            partition_of=points,
        ),
        measure=DiscreteMeasure(
            nodes=np.arange(n, dtype=float),
            weights=flat_weights,
            support=points,
        ),
    )

    # The induced reads — everything the operators retain comes off the
    # frame; the frame goes out of scope at return (forgetful).
    kernel_weights = frame.measure.weights
    total_weight = float(frame.discrete_gram[0, 0])

    retraction = AxisRetractionOperator(
        full_space=space,
        marginal_space=marginal_space,
        axis_shape=axis.shape,
        dims=dims,
        flat_weights=kernel_weights,
    )
    section = (
        None
        if total_weight == 0.0
        else AxisSectionOperator(
            full_space=space,
            marginal_space=marginal_space,
            axis_shape=axis.shape,
            dims=dims,
            flat_weights=kernel_weights,
            total_weight=total_weight,
        )
    )
    return _AxisCollapsePair(retraction, section)
