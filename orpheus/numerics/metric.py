r"""Hilbert metrics — a space's inner-product structure as an OBJECT.

The metric of a discrete Hilbert space is a symmetric bilinear form :math:`G`
(positive semi-definite, except a signed :class:`DiagonalMetric`): :math:`\langle x, y\rangle = y^{\mathsf T} G x`. A
:class:`DiagonalMetric` may carry signed weights and reads its signature exactly;
the :class:`~orpheus.numerics.factor.Factor` requires an invertible form, not a positive one, and the verbs that need
positivity ask for it, through the one predicate :attr:`HilbertMetric.positive`. Of the 5 verbs that need positivity, 3 refuse a form that
reads a negative direction: :meth:`~orpheus.numerics.space.FunctionSpace.norm` and the weighted orthonormalisation
(:meth:`~orpheus.numerics.basis.orthonormal_span_basis.OrthonormalSpanBasis.of_generators`), both through :func:`_require_positive`, and the
frame's mint, which refuses an under-resolved Gram (:func:`~orpheus.numerics.frame.frame_kernel`). The other 2 of the 5, the dense kernel
basis and the kernel projection gauge's displacement, do not refuse (an indefinite dense form has no producer, and the gauge path is refused
upstream); the norm of a composite space, which also does not refuse, is owed to issue #706.
Until P7 (campaign 1) the tree could spell exactly two realizations, both
*multiplied* into the element: a broadcast diagonal weight array (the
Hadamard product :math:`G \odot x`) and the per-axis factor measures of an
axis-built space. Neither can express a metric with off-diagonal structure,
and the tree carries frames whose discrete Gram has one. A frame's Gram is
diagonal only when its rule integrates every product of two basis functions
exactly; on a rule that does not, the Gram couples the modes, and then **no
diagonal metric satisfies Parseval on it**. `[M]` 2026-10-18 on
``frame_faces.angular_frame(Quadrature.folded_product(4, 4), 2)``, the mirror-even
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
``frame_faces.angular_frame(folded_product(2, 4), 3)``: ``[1.000, 3.000]`` diagonal,
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

* ``frame_faces.angular_frame(folded_product(4, 4), 2)``: :math:`15\times 15`, 6 live
  diagonal slots, **rank 5**. The other 9 slots are structurally dead
  layout. The one kernel mode inside the live block is structural, not
  small: the combination :math:`(-\sqrt3/2,\ 1/2)` of two
  :math:`\ell = 2` columns annihilates the evaluation **table** on the 8
  nodes (``‖A v‖∞ = 2.8e-16``,
  ``scratch/metric_docstring/probe_kernel.py``), so the rule cannot tell
  those two harmonics apart.
* ``frame_faces.angular_frame(folded_product(4, 4), 4)``: :math:`45\times 45`, 15 live
  slots, rank 8.
* ``frame_faces.angular_frame(product(4, 8), 4)``: :math:`45\times 45`, 25 live slots,
  rank 23.

On each of the three ``np.linalg.inv`` raises ``Singular matrix``, while
the exact inverse of :func:`_exact_inverse`, zeroing the declared kernel, is the object Parseval needs:
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
``A* = RieszRaiseOperator(A.domain) ∘ A.dual() ∘ RieszLowerOperator(A.codomain)``, which
:class:`~orpheus.numerics.operator.AdjointOperator` builds at construction
(CS4c R2, 2026-08-30). The legs wrap this family's two faces:
:class:`~orpheus.numerics.operator.RieszLowerOperator` applies
``apply_metric`` and
:class:`~orpheus.numerics.operator.RieszRaiseOperator` applies
``apply_inverse_metric``, and the space delegates both to its realization
here. So the metric arithmetic lives here once, and the adjoint holds the
composition, never a second spelling of it.

History
=======

* Until the ERR-080 repair (#429, 2026-09-02) the witness of both sections
  above was the slab ``frame_faces.angular_frame(gauss_legendre(8), 2)``. A 1-D rule
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
from enum import Enum, unique
from functools import cached_property
from typing import Protocol, TypeIs, cast

import numpy as np
from numpy.typing import NDArray

from .coordinate_space import leading_coordinates
from .declaration import Ann, Coord, Declaration, Spanned

__all__ = [
    "HilbertMetric",
    "DiagonalMetric",
    "DenseMetric",
    "Face",
    "Form",
    "NoForm",
    "DualMetric",
    "FactoredMetric",
]


#: Symmetry admission threshold for :class:`DenseMetric`, relative to the matrix's own scale. An asymmetric form is not an inner
#: product, so an asymmetric matrix is refused rather than silently symmetrised; the producer symmetrises, as
#: :meth:`DenseMetric.symmetrised` does for the frame's Gram. `[M]` 2026-10-18 over six dense frames: a
#: general-purpose inverse leaves ``max|M - M^T|`` up to ``1.75e-15`` (``frame_faces.angular_frame(product(4, 4), 2)``), so a
#: producer's round-off sits more than two decades below this threshold and what is refused is a form that is asymmetric in fact.
#:
#: ELEGANCE-DEBT[guard] #753 -- retires when a Gram producer returns an exactly symmetric type, so an asymmetric matrix is
#: unspellable.
_DENSE_METRIC_SYMMETRY_RTOL: float = 1e-12


def _exact_inverse(matrix: NDArray, kernel_dim: int) -> NDArray:
    r"""The Moore-Penrose inverse of the symmetric ``matrix``, zeroing exactly ``kernel_dim`` singular values.

    One body for both faces of :class:`DenseMetric`. The hermitian singular value decomposition
    :math:`G = U \Sigma V^{\mathsf T}` is read once; the ``kernel_dim`` smallest singular values are the declared kernel and
    are mapped to zero, the others to their reciprocals, and the inverse is :math:`V \Sigma^{+} U^{\mathsf T}`. No cutoff is read:
    the declaration is the rank. The split is well defined when every kept value exceeds every dropped one strictly, which a
    frame's own rule guarantees; a raw-path declaration at a tie is the caller's.

    Under ``if __debug__:`` the sign of ``u_i . vt_i`` (the sign of the i-th eigenvalue of a symmetric matrix) is tested on the
    kept directions, and a kept value that is not positive is refused: a kept exact zero contradicts the declared rank as surely as
    a negative one. This is a misuse check on the raw path, not a guarantee (it stays, as ruled for the raw constructor, and a dense form on a
    factor is positive semi-definite by the formal path), and it is not a rank detector: it sees only the sign of a
    round-off eigenvalue, so an under-declared kernel whose smallest kept value is a round-off residue is caught only when that residue
    happens to be negative or zero (the measured rate over the deficient harness Grams is on the theory page, under the metric of a
    function space). A kept value that is tiny and positive is accepted in every build, and the inverse carries its reciprocal. Under
    ``python -O`` a kept exact zero divides by zero and the faces are NaN.
    """
    u, s, vt = np.linalg.svd(matrix, hermitian=True)
    kept = s.shape[0] - kernel_dim
    if __debug__:
        signs = np.einsum("ji,ij->i", u, vt)[:kept]
        # A kept exact zero is read on the singular value itself: the sign the decomposition gives a zero eigenvalue is the
        # library's convention (0 on one build of numpy, +1 on another), so it is not asked.
        if not np.all(signs > 0.0) or np.any(s[:kept] == 0.0):
            raise ValueError(
                f"DenseMetric: a kept eigenvalue is not positive (the matrix is not positive semi-definite of rank "
                f"{kept}): the declared kernel of dimension {kernel_dim} does not contain every non-positive direction."
            )
    sinv = np.zeros_like(s)
    sinv[:kept] = 1.0 / s[:kept]
    return vt.T @ (sinv[:, None] * u.T)


def _inverse_on_live_block(matrix: NDArray, kernel: Declaration) -> NDArray:
    r"""The exact-count inverse of ``matrix`` under the declared ``kernel``: on a ``Coord`` kernel the inverse of the live block,
    embedded with exactly zero rows and columns at the named indices; else :func:`_exact_inverse` of the whole matrix."""
    if isinstance(kernel, Coord) and kernel.indices:
        live = [i for i in range(matrix.shape[0]) if i not in kernel.indices]
        out = np.zeros_like(matrix)
        out[np.ix_(live, live)] = _exact_inverse(matrix[np.ix_(live, live)], 0)
        return out
    return _exact_inverse(matrix, kernel.dim)


def _require_positive(form: "HilbertMetric", verb: str) -> None:
    r"""Refuse, for a verb that needs a positive form, a form with a negative direction: a real ``raise``, so it survives ``python -O``."""
    if not form.positive:
        minus = form.negative
        raise ValueError(
            f"{verb} needs a positive form: this form has {minus} negative direction(s), and the Krein arm of the form is not built (issue #706)"
        )


def _masked_reciprocal(weights: NDArray) -> NDArray:
    r"""The reciprocal of the positive weights and exactly zero elsewhere: the diagonal integrated end of a frame head.

    ELEGANCE-DEBT[guard] #756 the two tagged readers of a diagonal frame head (``Factor.form``, #756, and the space's
    ``_entry_of_axis``, #755) spell this arithmetic once, here; it retires when N04 and N10 give the diagonal integrated end one
    spelling. A diagonal-verdict Gram has no negative diagonal entry (``discrete_gram_structure`` routes one to the dense arm), so the
    ``> 0.0`` here never meets a negative weight.
    """
    live = weights > 0.0
    return np.where(live, 1.0 / np.where(live, weights, 1.0), 0.0)


def _require_square(matrix: NDArray) -> None:
    r"""Refuse a matrix that is not square and 2-D: the one shape check of both construction doors of :class:`DenseMetric`."""
    shape = np.shape(matrix)
    if len(shape) != 2 or shape[0] != shape[1]:
        raise ValueError(f"DenseMetric requires a square 2-D matrix, got shape {shape}")


def _parse_kernel(kernel: object, n: int, matrix: NDArray) -> Declaration:
    r"""``kernel`` as the declared kernel of a dense form on :math:`\mathbb{R}^n`, else the refusal: the one parse of both doors.

    A dense form's kernel is ``Coord(n, ())`` (full rank), a :class:`~orpheus.numerics.declaration.Spanned` (the decomposition
    drops the smallest singular directions, so a named index set or an annihilator could disagree with them), or a non-empty
    ``Coord`` exactly when ``matrix`` has exactly zero rows and columns at those indices: the named indices are then the kernel
    by the entries themselves, and no decomposition decides them.
    """
    if not isinstance(kernel, Declaration):
        raise TypeError(f"DenseMetric kernel must be a Declaration, got {type(kernel).__name__}")
    if kernel.n != n:
        raise ValueError(f"kernel is declared on R^{kernel.n}, the matrix is {n} x {n}")
    if isinstance(kernel, Coord) and kernel.indices and kernel.dim < n:
        named = list(kernel.indices)
        if not np.any(matrix[named, :]) and not np.any(matrix[:, named]):
            return kernel
    if isinstance(kernel, Ann) or (isinstance(kernel, Coord) and kernel.indices):
        raise ValueError(
            f"kernel is not declared by a spanned subspace: got {kernel!r}; a dense form's kernel is "
            f"Coord({n}, ()) (full rank) or a Spanned"
        )
    return kernel


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


class PositionedMetric(Protocol):  # concept: fmetric | role: interface | kind: morphism
    r"""A metric that applies to ONE interior index block of an element — the
    positioned leaf form a :class:`FactoredMetric` composes
    (:class:`DiagonalMetric`, :class:`DenseMetric`, and the
    :class:`DualMetric` of either)."""

    def apply_block(
        self,
        x: NDArray,
        *,
        start: int,
        block_shape: tuple[int, ...],
        inverse: bool,
    ) -> NDArray: ...


class HilbertMetric(ABC):  # concept: fmetric | role: head | kind: morphism
    r"""A symmetric bilinear form, as an object: positive semi-definite, except a signed :class:`DiagonalMetric`.

    A :class:`DiagonalMetric` is the one realization that may carry signed weights (it reads its :attr:`signature` exactly from
    the signs); every other realization is positive semi-definite by the formal path. A factor requires an invertible form, and
    positivity is asked by the verbs that need it: :meth:`~orpheus.numerics.space.FunctionSpace.norm` refuses a form that reads a
    negative direction (exact on a diagonal form, declared as 0 on a dense one); 3 of the 5 verbs that need positivity refuse (the module
    docstring names them).

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
    * :meth:`dual` — the metric of the dual space :math:`V^*`, the inverse
      Riesz form :math:`G^{+}` (:class:`DualMetric`, numerics_untangling.md
      ruling 26), built once per metric object.
    """

    #: The declared kernel of the form on the flattened index set: ``DiagonalMetric`` reads it from the exact zero weights,
    #: ``DenseMetric.kernel`` is its field and ``DualMetric`` takes its primal's. A :class:`FactoredMetric` declares none. An
    #: annotation and not an abstract property, because the dataclass field of ``DenseMetric`` is the implementation there.
    kernel: Declaration

    @property
    def negative(self) -> int:
        r"""The number of negative directions: 0 for a positive semi-definite form, the count of negative weights for a signed
        :class:`DiagonalMetric`."""
        return 0

    @property
    def positive(self) -> bool:
        r"""Whether the form has no negative direction: the one predicate of "this form is positive".

        :func:`_require_positive`, :attr:`~orpheus.numerics.factor.Factor.positive` and the weighted orthonormalisation
        (:meth:`~orpheus.numerics.basis.orthonormal_span_basis.OrthonormalSpanBasis.of_generators`) all ask it. Exact on a diagonal form
        (the signs of the weights) and on its dual and tensor products; on a dense form it is the declared value (``negative`` is 0),
        and the scope boundary of the declared signature (:attr:`~orpheus.numerics.factor.Factor.signature`, issue #706) covers it.
        """
        return self.negative == 0

    @property
    def signature(self) -> tuple[int, int, int]:
        r"""``(n_+, n_0, n_-)`` on the form's stored index set, read once over the declared :attr:`kernel` and the :attr:`negative` count.

        This is the signature of the form on its stored index set, with the kernel in ``n_0``; the signature of a factor's induced form
        on the quotient ``S/K`` is :attr:`~orpheus.numerics.factor.Factor.signature`, whose ``n_0`` is 0.

        A diagonal form's signature is exact, read from the signs of its weights. A dense form's signature is DECLARED by its kernel
        (``negative`` is 0): the sign check of the raw constructor runs under ``__debug__`` only, so under ``python -O`` a stated
        indefinite matrix reads ``(n - dim K, dim K, 0)``; the guarantee is the frame's door, the formal path.
        """
        n, zero, minus = self.kernel.n, self.kernel.dim, self.negative
        return (n - zero - minus, zero, minus)

    def dual(self) -> "HilbertMetric":
        r"""Return the metric of :math:`V^*`: this one with its two faces swapped.

        :math:`\langle \phi, \chi\rangle_{V^*} = \langle G^{+}\phi, \chi\rangle`,
        so the dual's forward face is this metric's inverse face and its
        inverse face is this forward face (:math:`(G^{+})^{+} = G` for a
        symmetric positive semi-definite form). Nothing is re-inverted, so
        the dual's arithmetic is bit-for-bit this metric's. Cached, and an
        involution by object identity: ``m.dual().dual() is m``.
        """
        return self._dual

    @cached_property
    def _dual(self) -> "HilbertMetric":
        return DualMetric(self)

    @abstractmethod
    def apply(self, x: NDArray) -> NDArray:
        r"""Return :math:`G x` (leading-aligned on the element's axes)."""

    @abstractmethod
    def apply_inverse(self, x: NDArray) -> NDArray:
        r"""Return :math:`G^{+} x` — pseudo-inverse on the range, 0 on the kernel."""

    @abstractmethod
    def validate_for(self, shape: tuple[int, ...]) -> None:
        r"""Raise if this metric cannot serve a space of element ``shape``."""

    @abstractmethod
    def rank(self, shape: tuple[int, ...]) -> int:
        r"""The rank of the form on elements of ``shape``: the dimension of
        :math:`V / \ker G`, the space the inner product is non-degenerate on
        (numerics_untangling.md ruling 74). A degenerate metric (rank below
        the element size: the tangential trace slots, dead moment-head
        slots, a rank-deficient Parseval head) is an inner product on that
        quotient, and the adjoint is posed there."""

    def kernel_basis(self, shape: tuple[int, ...]) -> NDArray:
        r"""An orthonormal basis of :math:`\ker G` on elements of ``shape``.

        The ``(n, n − rank)`` array whose columns span the kernel, ``n`` the
        element size, in the flat C-order chart: the kernel half of
        :meth:`_block_split`. This is the basis the kernel-condition check
        applies an operator to (#705): an operator :math:`A : V \to W` has an
        adjoint on the quotient iff :math:`P_W A k = 0` for every column
        :math:`k`.
        """
        return self._block_split(shape)[1]

    @abstractmethod
    def _block_split(self, block_shape: tuple[int, ...]) -> tuple[NDArray, NDArray]:
        r"""Orthonormal bases of :math:`\operatorname{range} G` and
        :math:`\ker G` on elements of ``block_shape``, the ``(n, rank)`` and
        ``(n, n − rank)`` columns in the flat C-order chart, read off the same
        stored data :meth:`rank` reads. :class:`FactoredMetric` takes the tensor
        product of its factors' splits, so a metric positioned as a factor
        block answers through this one law."""

    def pairing(self, x: NDArray, y: NDArray) -> float:
        r"""Return :math:`\langle x, y\rangle = \sum (Gx) \odot y` — the one spelling."""
        return float(np.sum(self.apply(x) * y))

    def diagonal(self, shape: tuple[int, ...]) -> NDArray | None:
        r"""Return :math:`G\mathbf 1` over elements of ``shape`` when :math:`G` is diagonal, else ``None``.

        A diagonal form is determined by its action on the constant, so the
        diagonal is :meth:`apply` of ones, exact (a product with ``1.0``). The
        base answer is ``None``, "not known to be diagonal": a realization
        that IS diagonal says so by overriding this, and a consumer that needs
        node weights (a frame's measure on a support) refuses ``None``.
        """
        return None


def _unit_columns(mask: NDArray) -> NDArray:
    r"""The unit vectors :math:`e_j` of the coordinates where ``mask`` is true, as columns."""
    indices = np.flatnonzero(mask)
    columns = np.zeros((mask.size, indices.size))
    columns[indices, np.arange(indices.size)] = 1.0
    return columns


def _frozen_copy(array: NDArray) -> NDArray:
    r"""A read-only copy: a metric is a frozen value, and a positioned form is
    space content (``FunctionSpace.identity_key`` digests it), so its arrays
    cannot be moved after the space is keyed."""
    out = np.array(array, dtype=float)
    out.setflags(write=False)
    return out


@dataclass(frozen=True)
class DiagonalMetric(HilbertMetric):  # concept: fmetric | role: arm | kind: morphism
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

    def __post_init__(self) -> None:
        # ``+ 0.0`` maps ``-0.0`` to ``+0.0``: the two zeros are one weight, one identity and one byte pattern.
        object.__setattr__(self, "weights", _frozen_copy(np.asarray(self.weights, dtype=float) + 0.0))

    @classmethod
    def counting(cls, shape: int | tuple[int, ...]) -> "DiagonalMetric":
        r"""The counting form: all-ones weights over ``shape`` (an int ``n`` is the shape ``(n,)``)."""
        return cls(np.ones((shape,) if isinstance(shape, int) else tuple(shape)))

    @property
    def kernel(self) -> Coord:  # pyright: ignore[reportIncompatibleVariableOverride]
        r"""The exact zero set of the weights, a ``Coord`` on the flattened index set; read from the weights themselves, no cutoff."""
        flat = np.asarray(self.weights).ravel()
        return Coord(flat.size, tuple(int(i) for i in np.flatnonzero(flat == 0.0)))

    @property
    def negative(self) -> int:
        return int(np.count_nonzero(np.asarray(self.weights) < 0.0))

    def apply(self, x: NDArray) -> NDArray:
        return _broadcast_leading(self.weights, np.ndim(x)) * x

    def diagonal(self, shape: tuple[int, ...]) -> NDArray:
        r"""The weights extended to ``shape``: :math:`G\mathbf 1`."""
        return self.apply(np.ones(shape))

    def apply_inverse(self, x: NDArray) -> NDArray:
        wb = _broadcast_leading(self.weights, np.ndim(x))
        nonzero = wb != 0.0
        return np.where(nonzero, x / np.where(nonzero, wb, 1.0), 0.0)

    def rank(self, shape: tuple[int, ...]) -> int:
        w = _broadcast_leading(self.weights, len(shape))
        return int(np.count_nonzero(np.broadcast_to(w, shape)))

    def _block_split(self, block_shape: tuple[int, ...]) -> tuple[NDArray, NDArray]:
        r"""The unit vectors of the nonzero- and zero-weight coordinates, the
        weights leading-aligned as :meth:`rank` and :meth:`apply` read them."""
        w = np.broadcast_to(_broadcast_leading(self.weights, len(block_shape)), block_shape)
        zero = np.asarray(w).ravel() == 0.0
        return _unit_columns(~zero), _unit_columns(zero)

    def validate_for(self, shape: tuple[int, ...]) -> None:
        # The legacy weights array may be broadcast-shaped (``(1, 4, 1)`` on a
        # ``(3, 4, 2)`` space): the weights are leading-aligned and must
        # extend to the element's shape, never the element to theirs. The
        # check runs once, here, when the space is built. Each element is
        # parsed per call as well: by the space
        # (``FunctionSpace.leading_coordinates_of``) and, below it, by the
        # metric's own legs (``apply_block`` and ``FactoredMetric``'s walk
        # parse their block with the coordinate-space primitive).
        aligned = _broadcast_leading(self.weights, len(shape)).shape
        try:
            extended = np.broadcast_shapes(aligned, shape)
        except ValueError:
            extended = None
        if len(aligned) > len(shape) or extended != tuple(shape):
            raise ValueError(
                f"DiagonalMetric weights of shape {np.shape(self.weights)} "
                f"cannot serve a space of shape {tuple(shape)}: leading-aligned, "
                f"they must extend to the element's shape."
            )

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
        leading_coordinates(out, out.shape[:start] + block_shape, owner="DiagonalMetric.apply_block")
        rank = len(block_shape)
        w = np.ascontiguousarray(
            np.broadcast_to(_broadcast_leading(self.weights, rank), block_shape)
        )
        wb = w.reshape(
            (1,) * start + block_shape + (1,) * (out.ndim - start - rank)
        )
        if inverse:
            nonzero = wb != 0.0
            return np.where(nonzero, out / np.where(nonzero, wb, 1.0), 0.0)
        return out * wb


@unique
class Face(Enum):  # concept: fmetric | role: component | kind: set
    r"""Which face of a dense form its stated matrix is: the ``metric`` (the matrix ``apply`` multiplies by) or the ``gram``
    (the matrix ``apply_inverse`` multiplies by, the inverse's role)."""

    METRIC = "metric"
    GRAM = "gram"


@dataclass(frozen=True)
class DenseMetric(HilbertMetric):  # concept: fmetric | role: arm | kind: morphism
    r"""A dense symmetric matrix metric on the flattened leading block.

    ``stated`` is the matrix the form is handed and ``face`` says which face it is, ``(K, K)`` with ``K`` the row-major flattening
    of the metric's block, the same layout as :attr:`~orpheus.numerics.frame.FrameBase.discrete_gram` (``table.reshape(N, K)``
    order), so a frame's Gram algebra and this realization agree on what "slot :math:`k`" means by construction. Elements with
    extra **trailing** axes are handled leading-aligned, like every metric in the tree: ``apply`` reshapes to ``(K, -1)``, matmuls,
    and reshapes back.

    The form has two faces, both functions of the fields: ``matrix`` (the ``metric`` face, what ``apply`` multiplies by) is
    ``stated`` when ``face`` is :attr:`Face.METRIC` and the exact-count inverse of ``stated`` when it is :attr:`Face.GRAM`;
    ``inverse_matrix`` (the ``gram`` face) is the other way. The content is ``(stated, face, kernel)``: a form stated as
    :math:`G^{+}` at ``METRIC`` and the same form stated as :math:`G` at ``GRAM`` are two values, and
    :func:`dataclasses.replace` keeps both faces' bytes. A raw, unsymmetrised matrix is taken as stored at either face;
    :meth:`symmetrised` is the one spelling of the exact symmetrisation a producer hands in.

    The rank is declared, never read from a cutoff. ``kernel`` is the declared kernel of the matrix: ``Coord(K, ())`` for a
    full-rank form, a :class:`~orpheus.numerics.declaration.Spanned` of dimension :math:`\dim\ker`, or a non-empty ``Coord``
    when the matrix has exactly zero rows and columns at those indices (the inverse is then taken on the live block and embedded
    with zero rows and columns). An ``Ann``, and a ``Coord`` whose rows are not zero, are refused. The signature of the form is
    :math:`(K - \dim\ker,\ \dim\ker,\ 0)` by the declaration.

    Construction refuses an **asymmetric** matrix (an asymmetric form is not an inner product: producers symmetrise, this type
    does not), a ``kernel`` on another dimension, and, under ``if __debug__:``, a kept direction whose eigenvalue is not positive.
    """

    stated: NDArray
    kernel: Declaration = field(kw_only=True)
    face: Face = field(kw_only=True)
    matrix: NDArray = field(init=False, compare=False)
    inverse_matrix: NDArray = field(init=False, compare=False)

    def __post_init__(self) -> None:
        if not isinstance(self.face, Face):
            raise TypeError(f"DenseMetric face must be a Face, got {type(self.face).__name__}")
        m = np.asarray(self.stated, dtype=float)
        _require_square(m)
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
        kernel = _parse_kernel(self.kernel, m.shape[0], m)
        stated = _frozen_copy(m)
        other = _frozen_copy(_inverse_on_live_block(m, kernel))
        object.__setattr__(self, "stated", stated)
        object.__setattr__(self, "matrix", stated if self.face is Face.METRIC else other)
        object.__setattr__(self, "inverse_matrix", other if self.face is Face.METRIC else stated)

    @staticmethod
    def symmetrised(gram: NDArray) -> NDArray:
        r"""The exact symmetrisation :math:`(G + G^{\mathsf T}) / 2`, spelled once; a matrix that is not square 2-D is refused."""
        g = np.asarray(gram, dtype=float)
        _require_square(g)
        return (g + g.T) / 2.0

    def flipped(self) -> "DenseMetric":
        r"""The same stated matrix read as the other face: the chart change of a dense form. Exact, because ``stated`` is already
        symmetric and nothing is inverted again by the caller."""
        return DenseMetric(self.stated, kernel=self.kernel, face=Face.GRAM if self.face is Face.METRIC else Face.METRIC)

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
        return self._matmul_leading(self.inverse_matrix, x)

    def rank(self, shape: tuple[int, ...]) -> int:
        r"""``(K - kernel.dim)`` times the trailing multiplicity ``size // K``: the declared rank, read from no cutoff."""
        size = int(np.prod(shape)) if shape else 1
        return (self.dim - self.kernel.dim) * (size // self.dim)

    def _block_split(self, block_shape: tuple[int, ...]) -> tuple[NDArray, NDArray]:
        r"""Orthonormal live and kernel bases of the block, from :meth:`_head_split`, on the leading ``dim`` slots and repeated over each
        trailing coordinate: column ``kron(b, e_t)`` (C-order, the block leading). The vectors are the eigenvectors of the symmetrised
        ``matrix`` and the kernel is its first ``kernel.dim`` columns, the count :meth:`rank` reads."""
        size = int(np.prod(block_shape)) if block_shape else self.dim
        trailing = np.eye(size // self.dim)
        live, kernel = self._head_split()
        return np.kron(live, trailing), np.kron(kernel, trailing)

    def _head_split(self) -> tuple[NDArray, NDArray]:
        r"""The ``(dim, rank)`` and ``(dim, dim − rank)`` bases on the block itself.

        The COUNT is the declared ``kernel.dim``, the one count :meth:`rank`
        reads, so the two cannot disagree. The VECTORS are the eigenvectors of :math:`G` itself, positive
        semi-definite, so its ``dim − rank`` smallest eigenvalues span the
        kernel: an eigen-solve of the projector :math:`I − G^{+}G` instead
        loses digits as :math:`\varepsilon / \lambda_{\min}` of the live
        spectrum (``[M]`` the A1 qa review, 2026-11-12: ‖GK‖/‖G‖ up to
        5.9e-6 at a live eigenvalue of 1.2e-12, against ≤ 2e-16 here)."""
        n_kernel = self.kernel.dim
        _, vectors = np.linalg.eigh(self.symmetrised(self.matrix))
        # eigh sorts ascending: the kernel is the first n_kernel columns.
        return vectors[:, n_kernel:], vectors[:, :n_kernel]

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
        out = np.asarray(x)
        leading_coordinates(out, out.shape[:start] + block_shape, owner="DenseMetric.apply_block")
        lead_shape = out.shape[:start]
        lead = int(np.prod(lead_shape)) if lead_shape else 1
        x3 = out.reshape(lead, self.dim, -1)
        return np.einsum("kl,alb->akb", m, x3).reshape(out.shape)


@dataclass(frozen=True)
class NoForm:  # concept: fmetric | role: component | kind: element
    r"""The absence of a form on an axis: no weight, no inner product, no norm.

    Not a :class:`HilbertMetric`: it has no ``apply``, ``inner_product``, ``norm`` or ``flat``, and each consumer matches on
    :data:`Form` exhaustively. The body carries a ``__repr__`` because a field-less empty body is the fingerprint of an unrelated
    class in the concept seed.
    """

    def __repr__(self) -> str:
        return "NoForm()"


#: The closed sum of what an axis carries as its form.
Form = DiagonalMetric | DenseMetric | NoForm


def _is_blank(entry: object) -> TypeIs[NoForm | None]:
    r"""True for a block of a :class:`FactoredMetric` that carries no form (``None``, Euclidean, or :class:`NoForm`)."""
    return entry is None or isinstance(entry, NoForm)


@dataclass(frozen=True)
class FactoredMetric(HilbertMetric):  # concept: fmetric | role: arm | kind: morphism
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
        tuple[tuple[int, ...], "PositionedMetric | NoForm | None"], ...
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

    @property
    def kernel(self) -> Declaration:  # pyright: ignore[reportIncompatibleVariableOverride]
        r"""A tensor product declares no kernel on the flattened index set: the kernel of a Kronecker product is a coordinate
        subspace only when every block's is, so the question has no one answer here and the read raises ``AttributeError``."""
        raise AttributeError("a FactoredMetric declares no kernel: read the kernel of its blocks")

    @property
    def signature(self) -> tuple[int, int, int]:
        r"""``(n_+, n_0, n_-)`` of the Kronecker product, derived from the blocks' signatures: a product of two directions is
        positive when their signs agree and negative when they differ, and zero when either is zero. A Euclidean block is
        all positive."""
        total, pos, neg = 1, 1, 0
        for block_shape, factor in self.entries:
            size = int(np.prod(block_shape)) if block_shape else 1
            p, _, m = (size, 0, 0) if _is_blank(factor) else cast(HilbertMetric, factor).signature
            total, pos, neg = total * size, pos * p + neg * m, pos * m + neg * p
        return (pos, total - pos - neg, neg)

    @property
    def negative(self) -> int:
        r"""The number of negative directions of the product, read off :attr:`signature`."""
        return self.signature[2]

    def _walk(self, x: NDArray, *, inverse: bool) -> NDArray:
        # The product acts on the concatenation of its blocks: a Euclidean
        # (``None`` or :class:`NoForm`) block applies nothing, so no block alone sees its extent.
        verb = "apply_inverse" if inverse else "apply"
        out = leading_coordinates(np.asarray(x), self._block_shape, owner=f"FactoredMetric.{verb}")
        start = 0
        for block_shape, factor in self.entries:
            if not _is_blank(factor):
                out = factor.apply_block(
                    out, start=start, block_shape=block_shape, inverse=inverse
                )
            start += len(block_shape)
        return out

    @cached_property
    def _block_shape(self) -> tuple[int, ...]:
        r"""The concatenation of the blocks' shapes, the leading shape the product acts on."""
        concat: tuple[int, ...] = ()
        for block_shape, _ in self.entries:
            concat = concat + block_shape
        return concat

    def apply(self, x: NDArray) -> NDArray:
        return self._walk(x, inverse=False)

    def diagonal(self, shape: tuple[int, ...]) -> NDArray | None:
        r"""Diagonal iff every block is Euclidean or diagonal; then :math:`G\mathbf 1`."""
        for block_shape, factor in self.entries:
            if not _is_blank(factor) and cast(HilbertMetric, factor).diagonal(block_shape) is None:
                return None
        return self.apply(np.ones(shape))

    def apply_inverse(self, x: NDArray) -> NDArray:
        return self._walk(x, inverse=True)

    @cached_property
    def _dual(self) -> "HilbertMetric":
        r"""The dual of a tensor product is the product of the block duals,
        :math:`(G_1 \otimes G_2)^{+} = G_1^{+} \otimes G_2^{+}`, so a dual
        factored metric stays a positioned product that a tensor product of
        spaces can flatten into its own entries. The involution by object
        identity is seeded on the dual."""
        dual = FactoredMetric(tuple(
            (
                block_shape,
                # A positioned leaf's dual is a DualMetric, itself positioned.
                None if _is_blank(factor)
                else cast(PositionedMetric, cast(HilbertMetric, factor).dual()),
            )
            for block_shape, factor in self.entries
        ))
        vars(dual)["_dual"] = self
        return dual

    def rank(self, shape: tuple[int, ...]) -> int:
        r"""The rank of a tensor product is the product of the block ranks
        (a Euclidean block is full rank)."""
        rank = 1
        blocks = 1
        for block_shape, factor in self.entries:
            size = int(np.prod(block_shape)) if block_shape else 1
            blocks *= size
            rank *= size if _is_blank(factor) else cast(HilbertMetric, factor).rank(block_shape)
        # Trailing element axes beyond the blocks multiply the rank, as they
        # do for a leading-aligned diagonal or dense metric.
        size = int(np.prod(shape)) if shape else 1
        return rank * (size // blocks)

    def _block_split(self, block_shape: tuple[int, ...]) -> tuple[NDArray, NDArray]:
        r"""The split of a tensor product: the range is the product of the
        factors' ranges, and the kernel every product of block basis vectors
        with AT LEAST ONE kernel factor, :math:`\ker(G_1\otimes G_2) =
        \ker G_1 \otimes V_2 + \operatorname{range} G_1 \otimes \ker G_2`; the
        count is :math:`n − \prod r_i`, the complement of :meth:`rank`."""
        live: NDArray = np.ones((1, 1))
        kernel: NDArray = np.zeros((1, 0))
        for entry_shape, factor in self.entries:
            size = int(np.prod(entry_shape)) if entry_shape else 1
            if _is_blank(factor):
                block_live, block_kernel = np.eye(size), np.zeros((size, 0))
            else:
                block_live, block_kernel = cast(HilbertMetric, factor)._block_split(entry_shape)
            kernel = np.hstack([
                np.kron(kernel, np.hstack([block_live, block_kernel])),
                np.kron(live, block_kernel),
            ])
            live = np.kron(live, block_live)
        size = int(np.prod(block_shape)) if block_shape else 1
        trailing = np.eye(size // kernel.shape[0])
        return np.kron(live, trailing), np.kron(kernel, trailing)

    def validate_for(self, shape: tuple[int, ...]) -> None:
        concat = self._block_shape
        if concat != shape:
            raise ValueError(
                f"FactoredMetric blocks {concat} do not concatenate to the "
                f"space shape {shape}."
            )


class DualMetric(HilbertMetric):  # concept: fmetric | role: arm | kind: morphism
    r"""The metric of a dual space :math:`V^*`: a primal metric with its two faces swapped.

    :math:`V^*`'s inner product is the inverse Riesz form,
    :math:`\langle \phi, \chi\rangle_{V^*} = \langle G^{+}\phi, \chi\rangle`
    (numerics_untangling.md ruling 26, issue #617). Every realization
    already stores both faces of its form (a diagonal metric its weights and
    their masked reciprocal; a dense one :math:`G` and :math:`G^{+}`; a
    factored one each block's two faces), so the dual is the same object
    read the other way round: :meth:`apply` is the primal's
    :meth:`~HilbertMetric.apply_inverse`, :meth:`apply_inverse` the primal's
    :meth:`~HilbertMetric.apply`, and :meth:`apply_block` the primal's
    block with the face flipped. One definition for every realization, and
    the arithmetic is bit-for-bit the primal's: a dual never re-inverts
    (``DiagonalMetric(1/w)`` would differ from the masked reciprocal in the
    last bit, and :math:`1/(1/w) \ne w`).

    Built by :meth:`HilbertMetric.dual`, never directly; its own
    :meth:`dual` is the primal, by object identity.
    """

    def __init__(self, primal: HilbertMetric) -> None:
        #: The metric of :math:`V` whose inverse form this is.
        self.primal = primal

    @property
    def kernel(self) -> Declaration:  # pyright: ignore[reportIncompatibleVariableOverride]
        return self.primal.kernel  # G⁺ and G share their kernel

    @property
    def negative(self) -> int:
        return self.primal.negative

    def dual(self) -> HilbertMetric:
        return self.primal  # (G⁺)⁺ = G, the involution by object identity

    def apply(self, x: NDArray) -> NDArray:
        return self.primal.apply_inverse(x)

    def diagonal(self, shape: tuple[int, ...]) -> NDArray | None:
        r"""Diagonal iff the primal is: then :math:`G^{+}\mathbf 1`, the primal's masked reciprocal."""
        if self.primal.diagonal(shape) is None:
            return None
        return self.apply(np.ones(shape))

    def apply_inverse(self, x: NDArray) -> NDArray:
        return self.primal.apply(x)

    def validate_for(self, shape: tuple[int, ...]) -> None:
        self.primal.validate_for(shape)

    def rank(self, shape: tuple[int, ...]) -> int:
        return self.primal.rank(shape)  # G⁺ and G share their range

    def _block_split(self, block_shape: tuple[int, ...]) -> tuple[NDArray, NDArray]:
        return self.primal._block_split(block_shape)  # G⁺ and G share range and kernel

    def apply_block(
        self,
        x: NDArray,
        *,
        start: int,
        block_shape: tuple[int, ...],
        inverse: bool,
    ) -> NDArray:
        r"""The positioned form :class:`FactoredMetric` composes: the primal's
        block with its face flipped."""
        # A DualMetric sits in a factored entry only over a positioned leaf
        # (FactoredMetric.dual distributes over its blocks).
        return cast(PositionedMetric, self.primal).apply_block(
            x, start=start, block_shape=block_shape, inverse=not inverse,
        )

    def __repr__(self) -> str:
        return f"DualMetric({self.primal!r})"
