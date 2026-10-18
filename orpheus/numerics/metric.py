r"""Hilbert metrics — a space's inner-product structure as an OBJECT.

The metric of a discrete Hilbert space is a symmetric positive-semi-definite
bilinear form :math:`G`: :math:`\langle x, y\rangle = y^{\mathsf T} G x`.
Until P7 (campaign 1) the tree could spell exactly two realizations, both
*multiplied* into the element: a broadcast diagonal weight array (the
Hadamard product :math:`G \odot x`) and the per-axis factor measures of an
axis-built space. Neither can express a metric with off-diagonal structure,
and the tree carries frames whose discrete Gram has one. A frame's Gram is
diagonal only when its rule integrates every product of two basis functions
exactly; on a rule that does not, the Gram couples the modes, and then **no
diagonal metric satisfies Parseval on it**. `[M]` 2026-10-18 on
``Quadrature.folded_product(4, 4).angular_frame(2)``, the mirror-even
spherical harmonics on :math:`S^2/\sigma_y`, the smallest of the three dense
angular heads the frame-bridge census found production building
(``scratch/numerics_untangling/frame_bridge/facts/report.md`` §5): the live block's largest off-diagonal is 1.000 of
the Cauchy–Schwarz scale, and over the whole band-limited subspace the
Parseval ratio ranges over ``[1.000, 2.000]`` under the best diagonal
candidate ``1/diag(G)`` and over ``[1.000, 1.000]`` under the dense
pseudo-inverse (``scratch/metric_docstring/probe_flagship.py``, a draw-free
generalized-eigenvalue range). The committed wrong-metric discriminator,
``test_no_diagonal_metric_can_satisfy_parseval_on_a_dense_frame`` in
``tests/gates/numerics/test_frame.py``, reads the same statistic on
``folded_product(2, 4).angular_frame(3)``: ``[1.000, 3.000]`` diagonal,
``[1.000, 1.000]`` dense.

This module makes the metric a thing that is **applied** rather than a thing
that is multiplied. The family owns the arithmetic; the Hadamard weight is
the diagonal special case it always was:

* :class:`DiagonalMetric` — the Hadamard realization: a weight array
  broadcast against the **leading** axes of the element (the tree's metric
  convention — :func:`_broadcast_leading`, the one home since
  ``FunctionSpace._broadcast_metric`` retired at CS4c step 6 item 6.2a).
* :class:`DenseMetric` — a dense symmetric matrix acting on the flattened
  leading block (row-major — the same flattening as
  :attr:`~orpheus.numerics.frame.FrameBase.discrete_gram`, whose
  pseudo-inverse is the realization's flagship occupant).

A :class:`~orpheus.numerics.space.FunctionSpace` resolves its metric
*source* (the ``metric`` field, the legacy ``inner_product_weights`` array,
or nothing) into one of these realizations exactly once and delegates its
three metric verbs — ``apply_metric`` / ``apply_inverse_metric`` /
``inner_product`` — to it. An axis-built space DERIVES its realization from
its axes — one :class:`FactoredMetric` with a :class:`DiagonalMetric` block
per weighted axis (CS4c step 6 item 6.2c-i; until then an inline per-axis
loop duplicated :meth:`DiagonalMetric.apply_block`): the axes ARE the
metric source there, and a measure is diagonal by nature — a Gram is a
*form*, a different concept, which is why
:class:`~orpheus.numerics.axis.Axis` never grows a matrix slot (the P7
design ruling: the generator induces, the space holds). Where a form is
needed on an axis-built space — a dense Gram on a moment head — it is a
:class:`DenseMetric` POSITIONED on that axis's block of the space's own
object (ruling R-6.2c-2, 2026-09-08), never on the axis.

The inverse face is Moore–Penrose, everywhere, by doctrine
==========================================================

``apply_inverse`` is the **pseudo**-inverse :math:`G^{+}`: the reciprocal on
the metric's range, zero on its kernel. This is not a numerical convenience:
on every dense Gram the tree builds, :math:`G` is singular, so the
pseudo-inverse is the only inverse that exists. `[M]` 2026-10-18
(``scratch/metric_docstring/probe_flagship.py``, numpy 2.4.4):

* ``folded_product(4, 4).angular_frame(2)``: :math:`15\times 15`, 6 live
  diagonal slots, **rank 5**. The other 9 slots are structurally dead
  layout. The one kernel mode inside the live block is structural, not
  small: the combination :math:`(-\sqrt3/2,\ 1/2)` of two
  :math:`\ell = 2` columns annihilates the evaluation **table** on the 8
  nodes (``‖A v‖∞ = 2.8e-16``,
  ``scratch/metric_docstring/probe_kernel.py``), so the rule cannot tell
  those two harmonics apart.
* ``folded_product(4, 4).angular_frame(4)``: :math:`45\times 45`, 15 live
  slots, rank 8.
* ``product(4, 8).angular_frame(4)``: :math:`45\times 45`, 25 live slots,
  rank 23.

On each of the three ``np.linalg.inv`` raises ``Singular matrix``, while
:func:`np.linalg.pinv` returns the exact object Parseval needs:
:math:`G G^{+} G = G` holds to a max-abs residual of ``1.8e-15`` on the
first (relative ``1.4e-16`` to :math:`\max|G| = 4\pi`), ``4.0e-15`` and
``3.0e-15`` on the other two. That identity makes
:math:`\|M\psi\|^2_{G^{+}} = c^{\mathsf T} G G^{+} G c = \|S_0 c\|^2_W` a
**theorem for any** :math:`G`, singular or not. The diagonal realization
has carried the same doctrine since the trace-metric work (zero-weight
tangential ordinates map to zero); the dense realization extends it
unchanged.

The pairing has ONE spelling: :math:`\sum (Gx) \odot y`
=======================================================

:meth:`HilbertMetric.pairing` is ``float(np.sum(self.apply(x) * y))`` for
every realization — the single home of the pairing arithmetic. Two reasons,
both load-bearing:

* **Bit-identity.** The legacy diagonal spelling ``np.sum(w * x * y)``
  evaluates left-to-right as ``(w*x)*y``, and ``apply(x)`` IS ``w_b * x``,
  so the reduction tree is preserved exactly. The matmul spelling
  ``y @ (diag(w) @ x)`` is NOT equivalent: `[M]` it differs on 60–70 %
  of draws at ``n = 15`` (the design scan's draw read 1360 of 2000),
  worst per-seed deviation banded 46–16384 ULP / rel ``9.2e-15`` to
  ``2.1e-12`` over 40 seeds (archivist census) — routing the shipped
  diagonal path through a densified matmul would move pinned numbers
  tree-wide. The bit-exact witness is
  ``tests/gates/numerics/test_dense_metric.py``'s pairing-spelling gate.
* **Single source.** ``AdjointOperator`` builds the Hilbert adjoint from
  ``apply_metric``/``apply_inverse_metric`` while the pairing that judges
  it comes from ``inner_product``; deriving the pairing FROM ``apply``
  makes the two agree by construction (the ERR-067 family — two metric
  spellings diverging silently — is unspellable here).

CS4c compatibility (the leg target)
===================================

The Hilbert adjoint is the Riesz-leg composition
``A* = A.domain.riesz_raise ∘ A.dual() ∘ A.codomain.riesz_lower``, which
:class:`~orpheus.numerics.operator.AdjointOperator` builds at construction
(CS4c R2, 2026-08-30). The legs wrap this family's two faces:
:attr:`~orpheus.numerics.space.FunctionSpace.riesz_lower` applies
``apply_metric`` and
:attr:`~orpheus.numerics.space.FunctionSpace.riesz_raise` applies
``apply_inverse_metric``, and the space delegates both to its realization
here. So the metric arithmetic lives here once, and the adjoint holds the
composition, never a second spelling of it.

History
=======

* Until the ERR-080 repair (#429, 2026-09-02) the witness of both sections
  above was the slab ``gauss_legendre(8).angular_frame(2)``. A 1-D rule
  then bound the full real spherical harmonics on a fabricated azimuth, and
  its Gram measured :math:`15\times 15` with 5 live slots and rank 4; the
  best diagonal dressing read a Parseval ratio of 1.806 there. After the
  repair that frame binds the Legendre basis and its Gram is diagonal,
  :math:`3\times 3`, rank 3, exactly :math:`2/(2\ell+1)`.
* Until CS4c R2 (2026-08-30) the section "CS4c compatibility" recorded the
  Riesz legs as debt that this family's faces would one day be wrapped by
  (the recorded-debt admonition on ``docs/theory/foundations/frame.rst``,
  now marked discharged).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from numpy.typing import NDArray

__all__ = [
    "HilbertMetric",
    "DiagonalMetric",
    "DenseMetric",
    "FactoredMetric",
]


#: Relative cutoff handed to :func:`np.linalg.pinv` by
#: :meth:`DenseMetric.inverse_of` (and by :class:`DenseMetric` when it must
#: derive its own inverse face). Pinned rather than left to numpy's default:
#: an implicit default is a silent dependency on a numpy version, and the
#: mutation battery's over-truncation arm is only meaningful against a
#: pinned value.
#:
#: **What it is FOR: structural RANK DETERMINATION, not noise suppression.**
#: A frame whose rule cannot resolve its basis has a Gram that is
#: rank-deficient by construction. `[M]` 2026-10-18 on
#: ``folded_product(4, 4).angular_frame(2)``: the dropped singular value is
#: ``2.2e-17`` of the largest, and its singular vector annihilates the
#: evaluation **table** (``‖A v‖∞ = 2.8e-16``), not merely the Gram. So the
#: dropped value is the floating-point residue of an EXACT zero, not a
#: small-but-real mode, and ``rcond`` separates a genuine range from a
#: structural kernel.
#:
#: The admissible window has TWO edges. Configuration for every number
#: below: the shipped :attr:`~orpheus.numerics.frame.FrameBase.discrete_gram`
#: of the three dense production heads, numpy 2.4.4, `[M]` 2026-10-18,
#: ``scratch/metric_docstring/probe_window.py``.
#:
#: * **Upper: over-truncation**, at the ratio of the smallest kept to the
#:   largest singular value: ``2.000000e-01`` on
#:   ``folded_product(4, 4).angular_frame(2)``, ``1.224672e-01`` on its
#:   ``angular_frame(4)``, ``6.666667e-02`` on
#:   ``product(4, 8).angular_frame(4)``. Construction cannot see it
#:   (``G⁺ G G⁺ = G⁺`` holds for a TRUNCATED pinv too, so
#:   :data:`_DENSE_METRIC_PENROSE_RTOL` is blind). The catchers are the
#:   Parseval gates on dense frames: `[M]` with the module global mutated
#:   in-process (``scratch/metric_docstring/mutate_rcond.py`` over
#:   ``tests/gates/numerics/test_frame.py`` and ``test_dense_metric.py``),
#:   ``rcond = 0.25`` reddens 5 of 82: the three dense rows of
#:   ``test_parseval_analysis_is_an_isometry_onto_its_image`` (the three
#:   heads above),
#:   ``test_the_dressing_lands_parseval_on_the_production_anisotropic_frame``
#:   (``product(4, 4).angular_frame(2)``, edge ``0.2``) and one
#:   ``test_dense_metric.py`` law; ``rcond = 0.5`` reddens 7, the three
#:   ``folded_product(2, 4).angular_frame(3)`` gates joining (edge
#:   ``1/3``); ``1e-12`` reddens 0.
#: * **Lower: a construction REFUSAL**, at the largest relative round-off
#:   residue as ``eigh`` reports it (the decomposition
#:   ``pinv(hermitian=True)`` actually cuts on): ``1.670e-16`` on
#:   ``folded_product(4, 4).angular_frame(2)``, ``1.46e-16`` on the other
#:   two. Below it the residue is inverted and
#:   :meth:`DenseMetric.__post_init__` refuses the pair: `[M]`
#:   :meth:`inverse_of` raises at half the residue and builds at twice it,
#:   on all three. So no ADMISSIBLE rcond can admit the kernel.
#:
#: ``1e-12`` sits ``11.30`` decades below the upper edge and ``3.78`` above
#: the lower on ``folded_product(4, 4).angular_frame(2)`` (``10.82`` and
#: ``3.84`` on ``product(4, 8).angular_frame(4)``). ⚠ The lower margin must
#: stay wide **because the number it is measured from is round-off with no
#: stable value**: `[M]` the same matrix's largest residue reads
#: ``1.67e-16`` (``eigh``), ``4.71e-17`` (``eigvalsh`` and
#: ``svd(hermitian=True)``) and ``2.17e-17`` (general ``svd``), a 7.7×
#: spread over four numpy routines, so it moves with the LAPACK driver, the
#: BLAS and the machine.
#:
#: **The census at this rcond.** `[M]` 2026-10-18,
#: ``scratch/probe_rcond_16_census.py`` re-run unchanged (``REFUSED`` =
#: ``basis_space`` raises, ``BREACH`` = builds but ``|Parseval − 1| >
#: 1e-12``, seed 1234): **0 of 105** slab rows (``gauss_legendre`` orders
#: ``2,3,4,5,6,7,8,9,11,12,16,17,24,32,33`` × ``L ∈ {0,1,2,3,4,5,7}``) and
#: **0 of 196** 3-D rows (``level_symmetric(2..16)``, ``product``,
#: ``folded_product``, ``lebedev(5..29)``, same ``L`` grid) are affected.
#: The census can read positive: with the global mutated to ``0.1`` it
#: reports 1 slab and 27 3-D rows in BREACH
#: (``scratch/metric_docstring/probe_census_control.py``). The dense arm is
#: taken on 59 of the 196 3-D rows and on 15 of the 105 slab rows, every one
#: of them with ``L ≥ n``, more Legendre modes than nodes (15 of the grid's
#: 16 such rows, each of rank ``n``). The 89 slab rows with ``L ≤ n − 1`` are
#: all DIAGONAL, exactly ``2/(2ℓ+1)`` with ``cond = 2L+1``, because
#: Gauss–Legendre integrates the Legendre products exactly to degree
#: ``2n − 1``. The smallest distance from the pin to a kept or a dropped
#: singular value is 4.44 decades on the slab (``GL2``, ``L = 7``) and 3.71
#: on 3-D (``product(3, 6)``, ``L = 7``).
#:
#: History. Until the ERR-080 repair (#429, 2026-09-02) this block was
#: measured on the slab ``gauss_legendre(8).angular_frame(2)``, whose Gram
#: was rank-deficient because a 1-D rule fabricated an azimuth; it first read
#: the dropped mode as a "noise floor" (refuted 2026-08-31), and the census
#: then found 31 of 105 slab rows affected and 0 of 196 3-D rows. The repair
#: retired all 31, as predicted. Its other prediction, that no
#: :class:`DenseMetric` would be built on a slab and the constant would
#: become a pure conditioning guard, did not hold: the dense arm survives on
#: under-resolved slab rows and on the 3-D heads above, where the kernel is
#: structural.
_DENSE_METRIC_RCOND: float = 1e-12

#: Symmetry admission threshold for :class:`DenseMetric`, relative to the
#: matrix's own scale. An asymmetric form is not an inner product, so an
#: asymmetric matrix is REFUSED rather than silently symmetrized — the
#: producer symmetrizes (as :meth:`DenseMetric.inverse_of` does). `[M]`
#: 2026-10-18 over six dense frames
#: (``scratch/metric_docstring/probe_hermitian.py``): ``np.linalg.pinv``
#: WITHOUT ``hermitian=True`` leaves ``max|M − Mᵀ|`` up to ``1.75e-15``
#: (``product(4, 4).angular_frame(2)``), and at most ``6.9e-17`` with it.
#: Both sit more than two decades below this threshold, so `[R]` what it
#: refuses is a form that is asymmetric in fact, not round-off. (On the
#: pre-ERR-080 slab Gram the non-hermitian spelling read ``4.74e-14``.)
_DENSE_METRIC_SYMMETRY_RTOL: float = 1e-12

#: Penrose-consistency admission threshold for an explicitly supplied
#: inverse face: ``M @ M⁺ @ M`` must reproduce ``M`` to this relative
#: tolerance, or the pair is refused — an inconsistent (matrix, inverse)
#: pair is an illegal state, not a representation choice.
_DENSE_METRIC_PENROSE_RTOL: float = 1e-10


def _broadcast_leading(w: NDArray, target_ndim: int) -> NDArray:
    r"""Pad ``w`` with trailing singleton axes up to ``target_ndim``.

    The tree's metric convention is **leading-aligned**: a metric spanning
    the first axes of an element acts on the full tensor by broadcasting
    over the trailing (element) axes. This is the arithmetic that lived in
    ``FunctionSpace._broadcast_metric`` until P7, which delegated here until
    it retired at CS4c step 6 item 6.2a (2026-09-07) — one home, one door;
    a no-op whenever ``w`` already spans every axis.
    """
    w = np.asarray(w)
    if w.ndim >= target_ndim:
        return w
    return w.reshape(w.shape + (1,) * (target_ndim - w.ndim))


class HilbertMetric(ABC):
    r"""A symmetric positive-semi-definite bilinear form, as an object.

    The contract is three verbs and one derived pairing:

    * :meth:`apply` — :math:`x \mapsto Gx` (the forward face).
    * :meth:`apply_inverse` — :math:`x \mapsto G^{+}x`, the Moore–Penrose
      pseudo-inverse face (reciprocal on the range, zero on the kernel).
    * :meth:`validate_for` — the admission check a
      :class:`~orpheus.numerics.space.FunctionSpace` runs at construction
      when this object is installed as its metric source.
    * :meth:`pairing` — :math:`\langle x, y\rangle = \sum (Gx) \odot y`,
      **final in spirit**: realizations inherit it so the pairing
      arithmetic has exactly one spelling (see the module docstring's
      bit-identity note).
    """

    @abstractmethod
    def apply(self, x: NDArray) -> NDArray:
        r"""Return :math:`G x` (leading-aligned on the element's axes)."""

    @abstractmethod
    def apply_inverse(self, x: NDArray) -> NDArray:
        r"""Return :math:`G^{+} x` — pseudo-inverse on the range, 0 on the kernel."""

    @abstractmethod
    def validate_for(self, shape: tuple[int, ...]) -> None:
        r"""Raise if this metric cannot serve a space of element ``shape``."""

    def pairing(self, x: NDArray, y: NDArray) -> float:
        r"""Return :math:`\langle x, y\rangle = \sum (Gx) \odot y` — the one spelling."""
        return float(np.sum(self.apply(x) * y))


@dataclass(frozen=True)
class DiagonalMetric(HilbertMetric):
    r"""The Hadamard realization: a diagonal weight array, leading-aligned.

    The special case every pre-P7 metric was. ``apply`` broadcasts the
    weights against the leading axes (:func:`_broadcast_leading`) and
    multiplies; ``apply_inverse`` is the masked reciprocal — the
    Moore–Penrose doctrine's diagonal form, operation-for-operation the
    arithmetic that lived in ``FunctionSpace._diagonal_apply_metric`` /
    ``_diagonal_apply_inverse_metric`` (bit-identity is a contract here,
    not an accident: the resolved legacy ``inner_product_weights`` path
    routes through this class).
    """

    weights: NDArray

    def apply(self, x: NDArray) -> NDArray:
        return _broadcast_leading(self.weights, np.ndim(x)) * x

    def apply_inverse(self, x: NDArray) -> NDArray:
        wb = _broadcast_leading(self.weights, np.ndim(x))
        nonzero = wb != 0.0
        return np.where(nonzero, x / np.where(nonzero, wb, 1.0), 0.0)

    def validate_for(self, shape: tuple[int, ...]) -> None:
        # The legacy dense-array contract: any array broadcast-compatible
        # with the element is admissible, and broadcast errors surface at
        # first application exactly as they always did. Nothing to refuse
        # at construction beyond what numpy will say more precisely later.
        return None

    def apply_block(
        self,
        x: NDArray,
        *,
        start: int,
        block_shape: tuple[int, ...],
        inverse: bool,
    ) -> NDArray:
        r"""Apply to ONE interior index block of ``x`` (``start`` leading
        ranks precede it) — the positioned form
        :class:`FactoredMetric` composes. THE per-axis arithmetic (explicit
        reshape: leading 1s, the block shape, trailing 1s) — since CS4c
        step 6 item 6.2c-i the only spelling of it: an axis-built space
        derives its metric as a ``FactoredMetric`` over its axes and
        applies it through this method (the inline
        ``FunctionSpace._apply_axes_weights`` loop it duplicated retired,
        `[M]` bit-identical)."""
        out = np.asarray(x)
        rank = len(block_shape)
        w = np.ascontiguousarray(
            np.broadcast_to(np.asarray(self.weights), block_shape)
        )
        wb = w.reshape(
            (1,) * start + block_shape + (1,) * (out.ndim - start - rank)
        )
        if inverse:
            nonzero = wb != 0.0
            return np.where(nonzero, out / np.where(nonzero, wb, 1.0), 0.0)
        return out * wb


@dataclass(frozen=True)
class DenseMetric(HilbertMetric):
    r"""A dense symmetric matrix metric on the flattened leading block.

    ``matrix`` is :math:`G`, ``(K, K)`` with ``K`` the row-major flattening
    of the metric's block — the same layout as
    :attr:`~orpheus.numerics.frame.FrameBase.discrete_gram`
    (``table.reshape(N, K)`` order), so a frame's Gram algebra and this
    realization agree on what "slot :math:`k`" means by construction.
    Elements with extra **trailing** axes are handled leading-aligned, like
    every metric in the tree: ``apply`` reshapes to ``(K, -1)``, matmuls,
    and reshapes back.

    Construction refuses two illegal states (the guard thresholds are the
    module constants — a gate on this type must quote them, per
    ``vv-principles`` #16):

    * an **asymmetric** matrix (an asymmetric form is not an inner
      product) — producers symmetrize, this type does not;
    * an explicitly supplied inverse face that fails the **Penrose
      identity** :math:`G G^{+} G = G` (an inconsistent pair is an illegal
      state, not a representation choice).

    When no inverse face is supplied it is derived here, once, as
    ``np.linalg.pinv(matrix, hermitian=True, rcond=_DENSE_METRIC_RCOND)``.
    Prefer :meth:`inverse_of` when what you hold is the GRAM whose inverse
    the metric should BE (the Parseval case): it keeps the exact Gram as
    the inverse face instead of a second pseudo-inversion.
    """

    matrix: NDArray
    inverse_matrix: Optional[NDArray] = field(default=None)

    def __post_init__(self) -> None:
        m = np.asarray(self.matrix, dtype=float)
        if m.ndim != 2 or m.shape[0] != m.shape[1]:
            raise ValueError(
                f"DenseMetric requires a square 2-D matrix, got shape "
                f"{np.shape(self.matrix)}"
            )
        scale = max(1.0, float(np.max(np.abs(m)))) if m.size else 1.0
        asym = float(np.max(np.abs(m - m.T))) if m.size else 0.0
        if asym > _DENSE_METRIC_SYMMETRY_RTOL * scale:
            raise ValueError(
                f"DenseMetric requires a symmetric matrix: "
                f"max|G - G^T| = {asym:.3e} exceeds "
                f"{_DENSE_METRIC_SYMMETRY_RTOL:.0e} of the matrix scale "
                f"{scale:.3e}. An asymmetric form is not an inner product — "
                f"symmetrize at the producer."
            )
        object.__setattr__(self, "matrix", m)
        if self.inverse_matrix is None:
            inv = np.linalg.pinv(
                m, hermitian=True, rcond=_DENSE_METRIC_RCOND
            )
            object.__setattr__(self, "inverse_matrix", inv)
        else:
            inv = np.asarray(self.inverse_matrix, dtype=float)
            if inv.shape != m.shape:
                raise ValueError(
                    f"DenseMetric inverse face shape {inv.shape} does not "
                    f"match the matrix shape {m.shape}"
                )
            residual = float(np.max(np.abs(m @ inv @ m - m)))
            if residual > _DENSE_METRIC_PENROSE_RTOL * scale:
                raise ValueError(
                    f"DenseMetric was handed an inconsistent inverse face: "
                    f"max|G G+ G - G| = {residual:.3e} exceeds "
                    f"{_DENSE_METRIC_PENROSE_RTOL:.0e} of the matrix scale "
                    f"{scale:.3e} (the Penrose identity fails)."
                )
            object.__setattr__(self, "inverse_matrix", inv)

    @classmethod
    def inverse_of(cls, gram: NDArray) -> "DenseMetric":
        r"""The metric whose MATRIX is the pseudo-inverse of ``gram``.

        The Parseval constructor: the frame's discrete Gram :math:`G` goes
        in, the metric :math:`G^{+}` comes out — with the (symmetrized)
        Gram itself installed as the inverse face, which is EXACT for a
        symmetric PSD form (:math:`(G^{+})^{+} = G` globally: both share
        the range, and :math:`G` is already zero on its own kernel) and
        strictly better conditioned than a second pseudo-inversion.
        """
        g = np.asarray(gram, dtype=float)
        g_sym = (g + g.T) / 2.0
        matrix = np.linalg.pinv(
            g_sym, hermitian=True, rcond=_DENSE_METRIC_RCOND
        )
        return cls(matrix=matrix, inverse_matrix=g_sym)

    @property
    def dim(self) -> int:
        r"""``K`` — the flattened slot count the matrix acts on."""
        return int(self.matrix.shape[0])

    def _matmul_leading(self, m: NDArray, x: NDArray) -> NDArray:
        flat = np.asarray(x).reshape(self.dim, -1)
        return (m @ flat).reshape(np.shape(x))

    def apply(self, x: NDArray) -> NDArray:
        return self._matmul_leading(self.matrix, x)

    def apply_inverse(self, x: NDArray) -> NDArray:
        inv = self.inverse_matrix
        assert inv is not None  # established by __post_init__; narrowing only
        return self._matmul_leading(inv, x)

    def validate_for(self, shape: tuple[int, ...]) -> None:
        size = int(np.prod(shape)) if shape else 1
        if size != self.dim:
            raise ValueError(
                f"DenseMetric of dimension {self.dim} cannot serve a space "
                f"of shape {shape} ({size} slots) — the matrix must span "
                f"the element's flattened leading block."
            )

    def apply_block(
        self,
        x: NDArray,
        *,
        start: int,
        block_shape: tuple[int, ...],
        inverse: bool,
    ) -> NDArray:
        r"""Apply to ONE interior index block of ``x`` — the positioned
        form :class:`FactoredMetric` composes: the element is reshaped to
        ``(lead, K, trail)`` around the block's row-major flattening and
        the matrix acts on the middle axis."""
        m = self.inverse_matrix if inverse else self.matrix
        assert m is not None  # established by __post_init__; narrowing only
        out = np.asarray(x)
        lead_shape = out.shape[:start]
        lead = int(np.prod(lead_shape)) if lead_shape else 1
        x3 = out.reshape(lead, self.dim, -1)
        return np.einsum("kl,alb->akb", m, x3).reshape(out.shape)


@dataclass(frozen=True)
class FactoredMetric(HilbertMetric):
    r"""A lazy tensor product of per-block metrics — :math:`G = G_1 \otimes G_2 \otimes \cdots`.

    Each entry is ``(block_shape, factor)``, the factor one of the
    POSITIONED leaf realizations (:class:`DiagonalMetric` /
    :class:`DenseMetric`) or ``None`` for a Euclidean block; the blocks'
    concatenation must equal the space shape (:meth:`validate_for`).
    Factors apply in sequence, each to its own index block — the
    Kronecker product is **never materialized**, the same discipline the
    per-axis path already follows ("never materializes the outer
    product"), so the pairing is exact and the cost is one pass per
    non-Euclidean factor.

    The occupant this realization exists for: a tensor product with a
    dense-metric factor — e.g. the harmonic frame's moment space, the
    Parseval-dressed spherical-harmonic coefficient block ⊗ the spatial
    cell measure — where the dense weights array ``*`` carried until CS4c
    step 6 item 6.2a could not carry the off-diagonal block and dropping
    it would silently revert the
    product to Euclidean on that factor (`[M]` 2026-08-30, the
    pre-repair behaviour: the probe pairing read 33.0 where
    :math:`G \otimes w` gives 109.0 — a value bug wearing a
    representation costume).
    """

    entries: tuple[
        tuple[tuple[int, ...], "DiagonalMetric | DenseMetric | None"], ...
    ]

    def __post_init__(self) -> None:
        for block_shape, factor in self.entries:
            if isinstance(factor, DenseMetric):
                size = int(np.prod(block_shape)) if block_shape else 1
                if size != factor.dim:
                    raise ValueError(
                        f"FactoredMetric entry: a DenseMetric of dimension "
                        f"{factor.dim} cannot serve a factor block of shape "
                        f"{block_shape} ({size} slots)."
                    )

    def _walk(self, x: NDArray, *, inverse: bool) -> NDArray:
        out = np.asarray(x)
        start = 0
        for block_shape, factor in self.entries:
            if factor is not None:
                out = factor.apply_block(
                    out, start=start, block_shape=block_shape, inverse=inverse
                )
            start += len(block_shape)
        return out

    def apply(self, x: NDArray) -> NDArray:
        return self._walk(x, inverse=False)

    def apply_inverse(self, x: NDArray) -> NDArray:
        return self._walk(x, inverse=True)

    def validate_for(self, shape: tuple[int, ...]) -> None:
        concat: tuple[int, ...] = ()
        for block_shape, _ in self.entries:
            concat = concat + block_shape
        if concat != shape:
            raise ValueError(
                f"FactoredMetric blocks {concat} do not concatenate to the "
                f"space shape {shape}."
            )
