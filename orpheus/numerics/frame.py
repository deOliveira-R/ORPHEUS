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
  (:func:`~orpheus.numerics.frame_faces.basis_space` ``= FunctionSpace.of_basis(basis)`` dressed with
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
=========================================================================

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

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.basis.base import Basis, GramStructure
from orpheus.numerics.declaration import Coord, Declaration, Spanned
from orpheus.numerics.metric import DenseMetric
from orpheus.values.content import content_digest
from orpheus.geometry.manifold import ManifoldMap, quotient_onto
from orpheus.numerics.measure import DiscreteMeasure


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
    (:func:`~orpheus.geometry.manifold.quotient_onto`: equality, the entry's
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
class FrameBase(ABC):  # concept: frame | role: head | kind: morphism
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
class PetrovGalerkinFrame(FrameBase):  # concept: frame | role: arm | kind: morphism
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
class GalerkinFrame(PetrovGalerkinFrame):  # concept: frame | role: arm | kind: morphism
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


#: The one rank rule of a frame's Gram: the count of its singular values at or below ``FRAME_GRAM_RANK_RTOL`` of the largest.
#: SCOPE-BOUNDARY[guard] machinery: an exactness-degree certificate on the quadrature (W8a); consumer: the harmonic head on
#: under-resolved rules; ruling: D9/U4, orchestrator, 2026-11-23; revisit: W8a's exactness certificate
FRAME_GRAM_RANK_RTOL: float = 1e-12


def _rank_and_signs(gram: NDArray) -> tuple[NDArray, int, NDArray]:
    r"""``(g_sym, dim, signs)``: the symmetrised Gram, the count of its singular values at or below ``FRAME_GRAM_RANK_RTOL`` of the
    largest, and the sign of each kept eigenvalue (``+1`` or ``-1``), all from the ONE ``svd(..., hermitian=True)`` of the
    symmetrised Gram.

    The sign of the ``i``-th eigenvalue of a symmetric matrix is ``u_i . vt_i`` of its hermitian singular value decomposition;
    the kept directions are the first ``n - dim`` (the singular values descend). A kept exact zero cannot occur (a kept value
    exceeds the cutoff), so no sign is read on a zero.
    """
    g = np.asarray(gram, dtype=float)
    if g.ndim != 2 or g.shape[0] != g.shape[1]:
        raise ValueError(f"a frame Gram is a square 2-D matrix, got shape {g.shape}")
    g_sym = DenseMetric.symmetrised(g)
    if g_sym.shape[0] == 0:
        return g_sym, 0, np.ones(0)
    u, singular, vt = np.linalg.svd(g_sym, hermitian=True)
    dim = int(np.sum(singular <= FRAME_GRAM_RANK_RTOL * singular.max()))
    return g_sym, dim, np.einsum("ji,ij->i", u, vt)[: singular.shape[0] - dim]


def _declared(g_sym: NDArray, dim: int) -> Declaration:
    r"""The kernel of dimension ``dim`` that the symmetrised Gram ``g_sym`` declares."""
    n = g_sym.shape[0]
    if dim == 0:
        return Coord(n, ())
    if dim == n:
        raise ValueError(f"a frame Gram has no value above {FRAME_GRAM_RANK_RTOL:.0e} of its largest: it is the zero form, which carries no inner product")
    g_sym.setflags(write=False)
    return Spanned(n, ("frame-gram-kernel", content_digest(g_sym)), dim)


def frame_gram_kernel(gram: NDArray) -> Declaration:
    r"""The kernel a frame declares for its discrete Gram ``gram``: the rank rule on a bare Gram.

    ``Coord(n, ())`` when no singular value of the symmetrised Gram is at or below ``FRAME_GRAM_RANK_RTOL`` of the largest,
    else a ``Spanned`` of the number of such values, keyed on the digest of the symmetrised Gram. The count reads the
    singular values of ``np.linalg.svd(g_sym, hermitian=True)`` with the vectors, the call :class:`DenseMetric` makes on the same
    matrix, so the declared count and the inverse the metric builds split one decomposition at the cutoff.
    """
    g_sym, dim, _ = _rank_and_signs(gram)
    return _declared(g_sym, dim)


def frame_kernel(frame: FrameBase) -> Declaration:
    r"""The kernel a frame declares for its Gram, through the formal door: the same declaration as :func:`frame_gram_kernel`, and on
    a measure with a negative weight (only then) the signs of the kept directions are read from the same decomposition and an
    indefinite Gram is refused.

    On a measure with no negative weight the Gram is positive semi-definite by construction and nothing is read. On a measure with a
    negative weight the Gram :math:`\sum_q w_q\,Y(\hat\Omega_q)Y(\hat\Omega_q)^{\mathsf T}` is the continuum Gram wherever the rule
    integrates a product of two basis functions exactly, and it can be indefinite where the rule does not: the basis is then
    under-resolved (a Gram block of total degree :math:`2L` exceeds the rule's exactness degree) and no dense form of the positive arm
    exists; the refusal is a real ``raise``, so it holds under ``python -O``. The diagonal verdict
    needs no counterpart: :attr:`FrameBase.discrete_gram_structure` routes a Gram with a negative diagonal entry to DENSE, so a frame
    that reaches the diagonal arm has a non-negative diagonal and a diagonal Gram is positive semi-definite.

    SCOPE-BOUNDARY[guard] a frame whose Gram is indefinite is refused at its mint.
    machinery: an exactness-degree certificate on the quadrature, which would decide this without reading the Gram (W8a)
    ruling: the user's ruling of 2026-11-25 (D0) and R1' (the orchestrator, 2026-11-24) as re-read at the opening of N04
    revisit: issue #706, the Krein arm
    """
    g_sym, dim, signs = _rank_and_signs(frame.discrete_gram)
    if _has_negative_weight(frame):
        negative = int(np.count_nonzero(signs < 0.0))
        if negative:
            raise ValueError(
                f"the frame's Gram has {negative} negative kept direction(s): the quadrature rule does not resolve the basis at this "
                f"order (a Gram block of total degree 2L exceeds the rule's exactness degree), and the Krein arm of the form is not "
                f"built (issue #706)"
            )
    return _declared(g_sym, dim)


def _has_negative_weight(frame: FrameBase) -> bool:
    r"""True when the frame's measure carries a negative weight (read exactly from the signs)."""
    return bool(np.any(np.asarray(frame.measure.weights) < 0.0))
