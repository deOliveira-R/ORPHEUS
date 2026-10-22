r"""The faces of a frame — the operators and spaces built FROM the frame's data.

A frame (:mod:`orpheus.numerics.frame`) is the join of a basis and a discrete
measure: the binding, its sampling table and its Gram. That datum sits upstream
of the axis, the space and the operators. Its **faces** are built from the frame
once spaces and operators exist, so they sit downstream, and this module is
their home (concept-order rulings 22 and 34: each face is a downstream
constructor of the frame, not an attribute of it):

* the analysis :math:`M = T`, :meth:`FrameAnalysis.of` ``(frame)``:
  ``measure_space(frame) → test_space(frame)``;
* the reconstruction :math:`R`, :meth:`FrameReconstruction.of` ``(frame)``:
  ``basis_space(frame) → measure_space(frame)``;
* the cross-Gram inverse :math:`G^{-1}`, :meth:`CrossGramInverse.of` ``(frame)``:
  ``test_space(frame) → basis_space(frame)``;
* the three dressed spaces :func:`basis_space` (the trial coefficient space
  carrying the Parseval metric), :func:`test_space` and :func:`measure_space`;
* the composed operators :func:`conjugate` (:math:`R\,A\,M`) and
  :func:`reconstruct_after` (:math:`R\,A`), and the coefficient-extraction verb
  :func:`project` (:math:`G^{-1}M f`).

**One build per frame instance.** Each face is built once and then returned as the
same object, through :func:`~orpheus.numerics.instance_cache.cached_on` under this
module's slot ``"_frame_faces"`` (on the quadrature too, for its angular frames,
keyed ``("angular_frame", L)``). Identities rest on it: a Galerkin frame's
:func:`test_space` IS its :func:`basis_space` (the analysis codomain ``is`` the
reconstruction domain), and every operator composed through a face carries the one
space object.

**The angular frame of a quadrature** is a constructor of a frame from a
quadrature, so it is homed here, downstream of both: the quadrature is a discrete
measure with a fibration and does not know which basis its orbit space admits, and
a frame is the join of a basis and a measure (concept-order ruling 21).
:func:`angular_frame` ``(quadrature, L)`` builds the Galerkin frame once per
(quadrature instance, ``L``) and caches it in the quadrature's own
``"_frame_faces"`` slot (the same slot and helper as the frame faces, keyed
``("angular_frame", L)``), so every consumer shares one frame object per rule and
order.

The frame itself carries none of these as an attribute: ``frame.analysis`` and the
other old spellings do not exist.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.axis import FromFrame
from orpheus.numerics.basis.base import Basis, GramStructure
from orpheus.numerics.basis.descent import Descent
from orpheus.numerics.basis.spherical_harmonic_basis import SphericalHarmonicBasis
from orpheus.numerics.frame import FrameBase, GalerkinFrame
from orpheus.numerics.instance_cache import cached_on
from orpheus.geometry.manifold import Manifold, Quotient, Sphere
from orpheus.numerics.metric import DenseMetric, DiagonalMetric, FactoredMetric
from orpheus.numerics.operator import LinearOperator, NotInvertible
from orpheus.numerics.projection import AnalysisOperator, ReconstructionOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace


#: This module's private slot on a frame or a quadrature (:func:`~orpheus.numerics.instance_cache.cached_on`).
_SLOT = "_frame_faces"

__all__ = [
    "CrossGramInverse",
    "FrameAnalysis",
    "FrameReconstruction",
    "angular_frame",
    "basis_space",
    "conjugate",
    "measure_space",
    "project",
    "reconstruct_after",
    "test_space",
]



def _harmonic_basis(support: Manifold, L: int) -> Basis:
    r"""The degree-:math:`L` basis a rule's function space admits, DERIVED from the point set its measure lives on.

    The measure's :attr:`~orpheus.numerics.measure.DiscreteMeasure.support`
    says which functions the rule can carry, and nothing else does
    (#429's thesis, *angular spaces are derived from symmetry*):

    * the bare sphere: the full real harmonics
      :class:`~orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis`;
    * an orbit space of the sphere: whatever
      :class:`~orpheus.numerics.basis.descent.Descent`'s discriminator
      says for that ENTRY: the Legendre basis on :math:`S^2/O(2)_a` (a
      1-D rule, THE FIX for ERR-080, 2026-09-02), the sigma-even harmonics
      on :math:`S^2/\sigma_a` (the fold; rectangular layout, odd columns
      structurally zeroed);
    * anything else refuses, naming the work.

    Until 2026-09-02 this dispatched on the ``folded_by`` TAG and
    knew nothing of the 1-D case: a 1-D rule got the full harmonics on a
    measure forged onto :math:`S^2`, which is the whole of ERR-080.
    """
    if isinstance(support, Sphere):
        return SphericalHarmonicBasis(L=L)
    if isinstance(support, Quotient) and isinstance(support.base, Sphere):
        return Descent.for_entry(support, L).frame_basis
    raise NotImplementedError(
        f"harmonic machinery on a rule whose measure lives on "
        f"{support.name}: only the sphere and its catalogued orbit "
        f"spaces carry a shipped angular basis."
    )


def angular_frame(quadrature: Quadrature, L: int) -> GalerkinFrame:
    r"""The degree-:math:`L` spherical-harmonic :class:`~orpheus.numerics.frame.GalerkinFrame` on ``quadrature``.

    Binds the basis the rule's orbit space admits (``_harmonic_basis``: the full
    harmonics on a sphere rule, the Legendre basis on a 1-D rule's
    :math:`S^2/O(2)_a`, the sigma-even harmonics on a fold) to **the
    quadrature's own measure**, always. Until 2026-09-02 a 1-D rule's measure
    was rebuilt here by padding :math:`\mu` to :math:`(\mu, 0, 0)` and calling
    the result :math:`S^2`, ERR-080's construction, retired with #429's fused
    commit; the route gate ``test_q8_4`` pins ``frame.measure is q.measure`` on
    12 of 12 shipped rules. The weights are the analysis metric. The frame is a
    pure-Galerkin frame (test IS trial: the basis is its own test basis); its
    :meth:`~orpheus.numerics.frame_faces.FrameAnalysis.of` face is the
    :math:`Y^* W` moment projection :math:`M` and its
    :meth:`~orpheus.numerics.frame_faces.FrameReconstruction.of` face the
    addition-theorem synthesis :math:`R`.

    ``frame.table`` IS the rule's angular basis table, the ONE spelling since
    2026-09-02: the pass-through accessor ``spherical_harmonics(L)`` (zero
    production consumers, and a name that became false the day a 1-D rule bound
    its Legendre basis) was retired with #429's fused commit; its nine test
    readers read this frame's table.

    **One frame per (quadrature instance, L)** (CS4c section 14.4, the hub):
    every consumer reaching the frame (S, F, the windowing method) shares ONE
    object and its cached table; the cache is the quadrature's private
    ``"_frame_faces"`` slot, a provenance that is never identity.

    The single source of the angular frame consumed by
    :class:`~orpheus.transport.operators.scattering.ScatteringOperator`: its
    section 5.6 kernel AND the in-sweep moment accumulation share THIS object,
    so the projection table is never re-evaluated or allowed to diverge.

    Naming tripwire: ``angular_frame(quadrature, L)`` names the quadrature's
    permanent physical **axis** (the :math:`S^2` directions), NOT the contingent
    spherical-harmonic basis, preserving the greppable
    ``angular_frame``/``spatial_frame``/``energy_frame`` family. The
    :class:`SphericalHarmonicBasis` hardcode is honest while the angular axis has
    exactly one basis; a *second* angular basis (a Fourier / Slepian angular
    expansion) parametrises the basis, ``angular_frame(quadrature, L, basis=...)``,
    a SIGNATURE change, and must NOT rename the function.
    """
    return cached_on(
        quadrature,
        _SLOT,
        ("angular_frame", L),
        lambda: GalerkinFrame(_harmonic_basis(quadrature.measure.support, L), quadrature.measure),
    )


class CrossGramInverse(LinearOperator):
    r"""The row-sum-collapsed inverse cross Gram :math:`G^{-1}` of a frame, as
    a typed arrow ``test_space → basis_space`` (CS4c step 6 item 6.2c-ii).

    The projection normalisation of :func:`project` and the middle
    factor of the :math:`G`-orthogonal projector
    ``conjugate(frame, CrossGramInverse.of(frame))`` :math:`= R\,G^{-1}M`. Its action is
    the Moore–Penrose reciprocal of the probe diagonal — ONE spelling, the
    :class:`~orpheus.numerics.metric.DiagonalMetric`'s own inverse
    arithmetic (zero on dead slots), broadcast against the coefficient
    array's leading block — and it reads NO space's metric, so a frame's
    Parseval dressing (diagonal or dense) can never leak into the
    normalisation. Self-adjoint as a real diagonal.

    Parameters
    ----------
    diagonal : NDArray
        The probe :math:`(M R\,\mathbf 1)` — the cross Gram's row sums, in the
        analysis codomain's layout.
    domain, codomain : FunctionSpace
        The analysis codomain (the TEST coefficient space) and the
        reconstruction domain (the TRIAL coefficient space) — one space on a
        Galerkin frame, two on a Petrov–Galerkin one.
    """

    def __init__(
        self, diagonal: NDArray, *, domain: FunctionSpace, codomain: FunctionSpace,
    ) -> None:
        self._metric = DiagonalMetric(np.asarray(diagonal, dtype=float))
        self._domain = domain
        self._codomain = codomain

    @classmethod
    def of(cls, frame: FrameBase) -> "CrossGramInverse":
        r"""The projection normalisation :math:`G^{-1}` of ``frame`` as an ARROW,
        ``test_space(frame) → basis_space(frame)``; built once per frame.

        :func:`project` normalises by the cross Gram :math:`G_{kj} = \langle\chi_k,
        \phi_j\rangle_W = (M R)_{kj}` of the analysis :math:`M` and reconstruction
        :math:`R`. This constructor takes a **single** ``analysis ∘ reconstruction``
        probe of the all-ones coefficient vector — the **row sum** of :math:`M R` —
        and returns its masked reciprocal as a typed operator between the two
        faces' coefficient ends (the Moore–Penrose pseudo-inverse zeroes empty /
        zero-weight regions for free), so :func:`project` is ``G⁻¹ ∘ M``. What
        ``conjugate(frame, CrossGramInverse.of(frame))`` then IS depends on the trial (the
        two cases below): on a disjoint-support trial it is the
        :math:`G`-orthogonal projector :math:`R\,G^{-1}M` onto ``span(basis)``;
        on a partition-of-unity trial it is NOT (see that bullet). The diagonal
        acquires the test weight's trailing (group, …) shape from the analysis
        face.

        ⭐ **Why an arrow and not a metric-twin space** (CS4c step 6 item
        6.2c-ii, 2026-09-08). Until then this face (``gram``) returned the
        TEST space re-dressed with the probe diagonal as its metric, and the
        projector was spelled ``conjugate(InverseMetricOperator(frame.gram))``
        — an endomorphism of that twin. Under structural space identity (the
        metric enters the identity of an axis-built head) the twin is a
        DIFFERENT space from ``test_space``, so the composition's own
        compatibility guard refused it on every harmonic frame; and the twin
        was always a costume: :math:`G^{-1}` maps the test space's covariant
        output to the trial space's contravariant coefficients — its ends ARE
        the faces' ends. Spelled as the arrow it composes by construction,
        and the hazard the old spelling had to guard against (a dense-dressed
        test space leaking its metric OBJECT into the probe — `[M]`
        2026-08-30, rel 1.625 on the overlap frame) is unspellable: the
        arrow's action reads the probe diagonal, never a space's metric.

        The row-sum probe is the weight :func:`project` needs under EITHER of two
        sufficient conditions — distinguish them, because the second does NOT imply a
        diagonal Gram:

        * **disjoint-support trial** (orthogonal harmonics, nested cell / group
          indicators): :math:`M R` is **diagonal**, off-diagonals structurally zero,
          so the row sum simply IS the diagonal :math:`(M R\,\mathbf 1)_R = \int_R
          w\,\mathrm{d}V` (:math:`w` the test weight).
        * **partition-of-unity trial** (the fractional
          :class:`~orpheus.numerics.basis.OverlapBasis`: a straddling fine cell lands
          in ≥2 coarse columns, so :math:`M R` is **NOT** diagonal): here the rows of
          the membership table sum to 1, i.e. :math:`R\,\mathbf 1 = \mathbf 1`, so the
          probe collapses to :math:`(M\,\mathbf 1)_R = \sum_i T_{iR}\,w_i = \Phi_R` —
          exactly the per-region weight, *even though the Gram is non-diagonal*. The
          conservative re-binning is correct because of this PoU collapse, NOT because
          the off-diagonals vanish. ⚠ What it computes is the per-region AVERAGE
          (the disintegration of the fine rate over each coarse region: conservative,
          positive, :math:`P\mathbf 1 = \mathbf 1`), which is the right object for
          condensation, and it is NOT the :math:`G`-orthogonal projector
          :math:`R\,(MR)^{-1}M`: on a straddling grid :math:`R\,G^{-1}M` is not
          idempotent (`[M]` 2026-09-28, defect 0.11 on a three-coarse-group
          fixture, 0.124 on a real straddling ``EnergyGrid``, 0 when nested) and
          differs from the orthogonal projection by 0.45 relative, while it
          preserves each coarse group's rate exactly (the design record, with
          the probe's command: ``.claude/plans/posing_sequence.md``, "the frame
          hierarchy").

        **A trial that is neither disjoint NOR a partition of unity** (a tapered
        weight, a higher-rank GEC moment — #275) makes the row-sum probe ≠ the true
        projection :math:`(M R)^{-1} M f`; such a basis needs the dense cross-Gram solve
        (the unbuilt least-squares seam), not this diagonal probe. That precondition is
        now carried by the **type**: the trial declares its
        :attr:`~orpheus.numerics.basis.base.Basis.gram_structure`, and this constructor
        **refuses** a :attr:`~orpheus.numerics.basis.base.GramStructure.DENSE` trial
        (raising :class:`~orpheus.numerics.operator.NotInvertible` — the projection
        normalisation :math:`G^{-1}` is not realizable as posed) rather than
        return a silently-wrong probe.
        """
        return cached_on(frame, _SLOT, "gram_inverse", lambda: cls._probe(frame))

    @classmethod
    def _probe(cls, frame: FrameBase) -> "CrossGramInverse":
        if frame.basis.gram_structure is GramStructure.DENSE:
            raise NotInvertible(
                f"frame_faces.project / CrossGramInverse.of needs a row-sum-collapsible trial Gram, but "
                f"the trial {type(frame.basis).__name__} declares GramStructure.DENSE: "
                f"its cross Gram MR is neither diagonal nor a partition of unity, so the "
                f"row-sum probe is NOT the projection normalisation G⁻¹. The dense "
                f"(MR)⁻¹M least-squares solve is not built (GitHub #275 — higher-rank "
                f"GEC moments / tapered weights); build it before projecting through "
                f"this basis."
            )
        ones = np.ones(basis_space(frame).shape)
        diagonal = np.asarray(FrameAnalysis.of(frame).apply(FrameReconstruction.of(frame).apply(ones)))
        return cls(
            diagonal, domain=test_space(frame), codomain=basis_space(frame),
        )

    @property
    def diagonal(self) -> NDArray:
        r"""The cross Gram's row sums :math:`(M R\,\mathbf 1)` this arrow inverts."""
        return self._metric.weights

    @property
    def domain(self) -> FunctionSpace:
        return self._domain

    @property
    def codomain(self) -> FunctionSpace:
        return self._codomain

    def apply(self, x: NDArray) -> NDArray:
        return self._metric.apply_inverse(x)

    def apply_transpose(self, x: NDArray) -> NDArray:
        # Self-adjoint: a real diagonal is its own transpose (and the two
        # ends share one shape).
        return self._metric.apply_inverse(x)

    @property
    def is_adjointable(self) -> bool:
        return True

    def __repr__(self) -> str:
        return f"CrossGramInverse({self._domain!r} → {self._codomain!r})"


def _dressed(
    space: FunctionSpace,
    *,
    frame: FrameBase,
    diagonal: NDArray | None = None,
    dense: DenseMetric | None = None,
) -> FunctionSpace:
    r"""Install a FRAME-owned metric on a basis's coefficient space — the
    one spelling of the F-0 dressing (:func:`basis_space`).

    Every coefficient space is AXIS-BUILT (the harmonic heads since CS4c
    step 6 item 6.2c-ii, the indicator families' since S7:
    :meth:`FunctionSpace.of_basis <orpheus.numerics.space.FunctionSpace.of_basis>`),
    so the metric lives where the axis doctrine put it:
    a diagonal dressing RE-WEIGHTS the head axis (the measure IS the
    metric), a dense one empties the axis's measure and POSITIONS the form
    on the space's derived metric object (item 6.2c-i) — and in both the
    axis's origin becomes :class:`~orpheus.numerics.axis.FromFrame`, which
    holds the object that can re-dress the head at another order
    (:meth:`~orpheus.numerics.frame.GalerkinFrame.at_order`) and re-derives
    the head's :attr:`~orpheus.numerics.axis.Axis.unit`: a dressed head holds
    COVARIANT moments, so the coordinates of the constant function there are
    its moments :math:`u = \Phi^{\mathsf T} w`, not the basis's coefficient
    :math:`e_0` (ruling 37). On a coarse indicator head those moments are the
    cell masses, so the dressed head is ``Modal(masses)`` (ruling 83, until
    S6c moves the moments to the dual).
    """
    if diagonal is None and dense is None:
        raise ValueError("_dressed: a dressing installs a diagonal or a dense form")
    axes = space.axes
    if axes is None or len(axes) != 1:
        raise NotImplementedError(
            f"a frame dresses a single-axis coefficient space; "
            f"{space.name!r} carries axes={axes!r}."
        )
    axis = axes[0]
    if dense is not None:
        return replace(
            space,
            axes=(replace(axis, weights=None, origin=FromFrame(frame)),),
            inner_product_weights=None,
            metric=FactoredMetric(((axis.shape, dense),)),
        )
    assert diagonal is not None  # narrowing only
    return replace(
        space,
        axes=(replace(
            axis, weights=np.asarray(diagonal).reshape(axis.shape), origin=FromFrame(frame),
        ),),
        inner_product_weights=None,
        metric=None,
    )


def basis_space(frame: FrameBase) -> FunctionSpace:
    r"""The codomain — the trial coefficient space, carrying the PARSEVAL metric.

    ``FunctionSpace.of_basis(basis)`` dressed with the inverse of the discrete trial Gram's
    diagonal (zero on dead slots). The analysis face's outputs are the
    COVARIANT moments :math:`\varphi = G c` of a band-limited field
    :math:`S_0 c`, so the inner product that makes analysis an isometry
    onto its image — Parseval, :math:`\|\varphi\|_{G^{-1}} = \|S_0 c\|_W`
    — is :math:`G^{-1}`, the **inverse of the discrete Gram**. That is a
    property of the *pairing* (basis ⊗ measure), so the frame owns it: the
    basis's own ``space`` keeps the continuum Gram (for the SH basis
    :math:`4\pi/(2\ell+1)` — the cross-Gram vocabulary of
    :class:`CrossGramInverse`/:func:`project`), and the frame REPLACES the metric with
    the measured inverse. With this metric the faces' ``.H`` is the
    physical Hilbert adjoint, and for the SH frame the square
    :math:`M^{*} = R/W`, :math:`R^{*} = W\,M` (with :math:`W = \sum_n w_n`)
    closes with one scalar exactly when the discrete trial Gram equals the
    continuum Gram on the live modes, i.e. when the rule integrates every
    product of two live basis functions exactly. `[M]` 2026-10-18
    (``scratch/numerics_untangling/frame_bridge/attack_probes/p6_tree_square.py``,
    37 rule × order rows): it closes to ≤1e-12 on 34 and fails on the 3
    whose order outruns the rule — ``folded_product(2, 4)`` at L=2
    (0.365), ``level_symmetric(4)`` at L=3 (0.067), ``level_symmetric(8)``
    at L=5 (0.039). The SN scattering adjoint does not rest on the square:
    it is standard S_N, :math:`(1/W)\,R\Lambda^{\mathsf T}M`, self-adjoint
    in the quadrature weights on every rule, and forward/adjoint
    reciprocity holds to round-off on those 3 rows
    (``scratch/adjoint_underresolved/``). A spelling that daggers the two
    faces separately and assumes the square closes is wrong there by 3 % to
    26 % (ERR-039's class). The pre-F-0 stored continuum metric was the WRONG
    side — off the physical adjoint by exactly :math:`(4\pi/(2\ell+1))^2`
    per ℓ. Live witnesses: ``tests/gates/numerics/test_frame.py``'s
    ``test_parseval_*`` suite, whose negative leg re-installs the
    pre-repair metric in-process; discovery record
    ``scratch/probe_f1_parseval.py`` — note the probe reads the frame's
    NOW-DRESSED space, so its "stored" row prints 1.000 post-repair.

    For a frame whose discrete Gram is NOT diagonal
    (:attr:`~orpheus.numerics.frame.FrameBase.discrete_gram_structure` DENSE — e.g.
    ``folded_product(4, 4).angular_frame(2)``, where NO diagonal metric
    satisfies Parseval), the metric is the
    matrix pseudo-inverse :math:`G^{+}` installed as a
    :class:`~orpheus.numerics.metric.DenseMetric` (campaign 1, P7):
    Parseval is then a THEOREM for any Gram, singular or not
    (:math:`\|Mc\|^2_{G^{+}} = c^{\mathsf T} G G^{+} G c =
    \|S_0 c\|^2_W`), and the faces' ``.H`` is the physical Hilbert
    adjoint on every frame, not only the diagonal ones. `[M]`
    2026-10-18 on ``folded_product(4, 4).angular_frame(2)``, draw-free
    over the band-limited subspace: the dense dressing's Parseval ratio
    ranges over [1.000, 1.000], the best diagonal candidate
    ``1/diag(G)``'s over [1.000, 2.000] and the undressed continuum
    metric's over [6.317, 157.9] — a diagonal metric there is not
    merely unavailable but provably insufficient
    (``scratch/metric_docstring/probe_theory.py``; on the pre-ERR-080
    slab GL8/L=2 frame, 2026-08-30, one draw read 1.000000000000,
    1.806 and 25.53). (Until P7 the DENSE arm returned the space
    UNDRESSED — *"Parseval is unavailable"* was the recorded F-0
    limitation, with the matrix-metric home deferred to the CS4c
    Riesz-leg machinery; P7 landed the metric object CS4c's legs will
    wrap, and the refusal era's record lives in this file's history
    and the error catalog's ERR-039 entry.)

    ⭐ Since CS4c step 6 item 6.2c-ii (2026-09-08) the dressing is
    VISIBLE to equality: the head is axis-built and its measure enters
    the identity, so ``basis_space(frame) != FunctionSpace.of_basis(basis)`` (`[M]` 0 of 33
    shipped (rule, L) rows equal), the frame becomes the head axis's
    generator, and the tree binds THIS space everywhere the moment
    space is spelled (ruling R-6.2c-1) — the metric-blind seam that
    let the two pass for one is gone.
    """
    return cached_on(frame, _SLOT, "basis_space", lambda: _parseval_dressed(frame))


def _parseval_dressed(frame: FrameBase) -> FunctionSpace:
    r"""The trial coefficient space dressed with the Parseval metric (see :func:`basis_space`)."""
    space = FunctionSpace.of_basis(frame.basis)
    if frame.discrete_gram_structure is not GramStructure.DIAGONAL:
        # DENSE verdict (P7): the Parseval metric exists — it is just
        # not diagonal. Install the matrix pseudo-inverse of the
        # measured Gram (DenseMetric.inverse_of keeps the exact
        # symmetrized Gram as the inverse face) and STRIP the basis's
        # continuum measure — the dressing REPLACES the metric on
        # this arm exactly as the diagonal arm overwrites it (on an
        # axis-built head the form is POSITIONED on the derived object
        # and the axis carries no measure, item 6.2c-i).
        return _dressed(
            space, dense=DenseMetric.inverse_of(frame.discrete_gram), frame=frame,
        )
    diag = np.diagonal(frame.discrete_gram).reshape(space.shape)
    live = diag > 0.0
    inverse = np.where(live, 1.0 / np.where(live, diag, 1.0), 0.0)
    return _dressed(space, diagonal=inverse, frame=frame)


def test_space(frame: FrameBase) -> FunctionSpace:
    r"""The analysis codomain — the test basis's coefficient space.

    The predicate is the frame's data, ``frame.test is frame.basis`` (ruling
    58): when the test basis IS the trial basis, the analysis codomain is the
    trial coefficient space the frame dressed, :func:`basis_space` (the same
    cached object, preserving the analysis-codomain ``is`` identity of the
    single-discipline frame the hierarchy replaced), whatever class spelled
    the frame — a :class:`~orpheus.numerics.frame.GalerkinFrame` or a
    Petrov–Galerkin frame handed its own trial basis as its test. A distinct
    test basis reads its own coefficient space,
    :meth:`FunctionSpace.of_basis(frame.test) <orpheus.numerics.space.FunctionSpace.of_basis>`.
    """
    return cached_on(frame, _SLOT, "test_space",
        lambda: basis_space(frame) if frame.test is frame.basis else FunctionSpace.of_basis(frame.test),
    )


def measure_space(frame: FrameBase) -> FunctionSpace:
    r"""The domain — the measure's induced discrete-:math:`L^2` space.

    :meth:`FunctionSpace.of_measure
    <orpheus.numerics.space.FunctionSpace.of_measure>` over the frame's
    measure, labelled with the measure's own
    :attr:`~orpheus.numerics.measure.DiscreteMeasure.axis_label` (per-node
    values with the quadrature weights as the metric). On a directional rule
    this IS the SN angular space over the same rule (rulings 63 and 71).
    """
    return cached_on(frame, _SLOT, "measure_space",
        lambda: FunctionSpace.of_measure(frame.measure, frame.measure.axis_label),
    )


@dataclass(frozen=True)
class FrameAnalysis(AnalysisOperator):
    r"""The analysis face :math:`M = T`: ``measure_space → test_space``.

    A frame-backed :class:`AnalysisOperator` view; the math lives on the TEST basis
    (:meth:`Basis.analyze` / :meth:`Basis.analyze_transpose`) tabulated at the nodes
    (``frame.test_table``). Carries the swapped spaces and a working
    ``apply_transpose`` (``is_adjointable=True``), so the metric-aware
    ``AdjointOperator`` gives ``.H`` (the W-weighted Hilbert adjoint)
    for free. For a :class:`~orpheus.numerics.frame.GalerkinFrame` the test basis is
    the trial basis, so this is the :math:`Y^* W` projection bit-identical to the
    single-discipline frame.
    """

    frame: FrameBase

    @classmethod
    def of(cls, frame: FrameBase) -> "FrameAnalysis":
        r"""The analysis face of ``frame``, built once per frame."""
        return cached_on(frame, _SLOT, "analysis", lambda: cls(frame))

    @property
    def domain(self) -> FunctionSpace:
        return measure_space(self.frame)

    @property
    def codomain(self) -> FunctionSpace:
        return test_space(self.frame)

    def apply(self, values: NDArray, /) -> NDArray:
        return self.frame.test.analyze(
            values, self.frame.test_table, self.frame.measure.weights,
        )

    def apply_transpose(self, coefficients: NDArray, /) -> NDArray:
        return self.frame.test.analyze_transpose(
            coefficients, self.frame.test_table, self.frame.measure.weights,
        )


@dataclass(frozen=True)
class FrameReconstruction(ReconstructionOperator):
    r"""The reconstruction face :math:`R`: ``basis_space → measure_space``.

    A frame-backed :class:`ReconstructionOperator` view delegating to the TRIAL basis's
    :meth:`Basis.reconstruct` (the canonical-dual synthesis) and its representation
    transpose :meth:`Basis.reconstruct_transpose` — reconstruction is purely trial-side,
    identical across disciplines. Carries the swapped spaces and a working
    ``apply_transpose`` (``is_adjointable=True``), so the metric-aware
    ``AdjointOperator`` gives ``R.H`` for free — symmetric with the
    analysis face.
    """

    frame: FrameBase

    @classmethod
    def of(cls, frame: FrameBase) -> "FrameReconstruction":
        r"""The reconstruction face of ``frame``, built once per frame."""
        return cached_on(frame, _SLOT, "reconstruction", lambda: cls(frame))

    @property
    def is_adjointable(self) -> bool:
        # The frame's reconstruction FACE adds apply_transpose on top of
        # the bare ReconstructionOperator role (apply-only), so ``R.H``
        # is free here. The override lives on the face.
        return True

    @property
    def domain(self) -> FunctionSpace:
        return basis_space(self.frame)

    @property
    def codomain(self) -> FunctionSpace:
        return measure_space(self.frame)

    def apply(self, coefficients: NDArray, /) -> NDArray:
        return self.frame.basis.reconstruct(coefficients, self.frame.table)

    def apply_transpose(self, values: NDArray, /) -> NDArray:
        return self.frame.basis.reconstruct_transpose(values, self.frame.table)


# ── composed operators (the "define Frame, compose, done" production path) ──


def conjugate(frame: FrameBase, operator: LinearOperator, /) -> LinearOperator:
    r"""Frame-conjugate a coefficient-space operator: :math:`R \circ A \circ M`.

    THE production composition for a method whose action is "project to
    coefficients, act there, reconstruct" — e.g. SN anisotropic scattering
    :math:`S_{\ell\ge 1} = R\,\Lambda\,M` (the per-ordinate redistribution). Returns
    the :class:`~orpheus.numerics.operator.OperatorProduct` word
    ``R @ A @ M`` (factors ``(R, A, M)``; a scaled ``A = cΛ`` gives the
    one-term combination ``c·(R @ Λ @ M)``, the coefficient moved out), whose
    :meth:`~orpheus.numerics.operator.OperatorProduct.apply` is ``R(A(M·x))`` — the
    SAME numpy chain a hand-rolled ``reconstruction.apply(A.apply(analysis.apply(x)))``
    runs, now ONE named operator (Cardinal Rule 2: the composition IS the production
    path, not a parallel "semantic" reading of it).

    **The double-category 2-cell — and, for an eigenbasis frame, the spectral
    theorem.** :func:`conjugate` is the 2-cell of the (Representation × Role)
    carrier grid: a vertical Role-morphism ``A`` (e.g. scattering's :math:`\Lambda`,
    the diagonal :math:`\Sigma_{s,\ell}` multiply) conjugated by the horizontal
    Representation-adjoint pair :math:`(M, R)`. When the Frame is an operator's
    EIGENBASIS — the SH angular frame is the scattering kernel's, by Funk–Hecke —
    :math:`R\circ\Lambda\circ M` IS the spectral theorem :math:`A = U\Sigma U^{*}`
    written out (:math:`M=U` analysis into the eigenbasis, :math:`\Lambda=\Sigma`
    the spectrum :math:`\Sigma_{s,\ell}`, :math:`R=U^{*}` synthesis). The frame is
    then OWNED by that operator (scattering owns its angular frame), not by the
    phase space — see :ref:`frame-eigenbasis-ownership`. The 0-ULP
    ``test_scattering_kernel_crosscheck`` is this 2-cell's interchange-law witness.

    ``operator`` must compose between the faces — its ``domain`` is the analysis
    codomain (:func:`test_space`) and its ``codomain`` the reconstruction domain
    (:func:`basis_space`); the :class:`OperatorProduct` space-compatibility guard
    enforces it (an endomorphism on the coefficient space when ``test == trial``).
    """
    return FrameReconstruction.of(frame) @ operator @ FrameAnalysis.of(frame)


def reconstruct_after(frame: FrameBase, operator: LinearOperator, /) -> LinearOperator:
    r"""Reconstruct after a coefficient-space operator: :math:`R \circ A`.

    The :func:`conjugate` sub-operator for inputs ALREADY in coefficient space (the
    analysis :math:`M` already applied) — e.g. the angular-windowed SN moment
    iterate, whose bulk IS :math:`M\psi`, so only :math:`R\,\Lambda` remains.
    Returns the word ``R @ A`` (``apply`` = ``R(A·c)``). Wiring a windowed
    consumer to :func:`conjugate` instead would erroneously re-apply :math:`M` (a
    double-projection).
    """
    return FrameReconstruction.of(frame) @ operator


# ── the coefficient-extraction verb (homogenise / condense) ──────────


def project(frame: FrameBase, field: NDArray, /) -> NDArray:
    r"""Extract coefficients: :math:`G^{-1} M f` — the homogenise / condense verb.

    The Petrov-Galerkin coefficient extraction: analyse the field against the
    test functions (:math:`(M f)_k = \langle\chi_k, f\rangle_W`), then normalise
    by the cross Gram :math:`G` (:meth:`CrossGramInverse.of`). For flux-weighted spatial homogenisation
    (``test`` :math:`= \varphi\cdot\mathbf 1_R`, ``trial`` :math:`= \mathbf 1_R`)
    this is the rate-preserving effective cross section :math:`\Sigma_R = \int_R
    \varphi\Sigma\,\mathrm{d}V / \int_R\varphi\,\mathrm{d}V`; for a
    :class:`~orpheus.numerics.frame.GalerkinFrame` (``test = trial``) it is the orthogonal projection onto
    the coarse space. The normalisation is a reciprocal (the diagonal :class:`CrossGramInverse`
    probe) for both the disjoint-indicator / orthogonal-harmonic consumers (Gram
    genuinely diagonal) AND the partition-of-unity
    :class:`~orpheus.numerics.basis.OverlapBasis` (non-diagonal Gram, but the
    row-sum probe is still the right per-region weight — see :meth:`CrossGramInverse.of`); the
    dense solve is the (unbuilt) least-squares seam only. Trailing (group, …)
    axes ride the analysis and broadcast against the diagonal Gram (so a vector
    channel divides by :math:`\Phi_{R,g}` and a ``[g_from, g_to]`` matrix channel
    by its source-group :math:`\Phi_{R,g_{\mathrm{from}}}`).
    """
    return CrossGramInverse.of(frame).apply(FrameAnalysis.of(frame).apply(field))
