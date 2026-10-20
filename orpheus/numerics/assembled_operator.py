r"""The assembly functor out of the operator algebra, and its sparse carrier.

Assembly is a homomorphism out of the algebra (numerics_untangling.md ruling
24): every operator that knows its stencil emits its own raw sparse matrix,
and the composites are realised by the laws

.. math::

   [\textstyle\sum_i c_i w_i] = \sum_i c_i [w_i], \qquad
   [A_1 \cdots A_n] = [A_1] \cdots [A_n], \qquad
   [\text{atom}] = [\text{atom.expansion()}].

The laws live here, downstream of :mod:`orpheus.numerics.operator`, which
knows nothing of assembly beyond the structural flag
:attr:`~orpheus.numerics.operator.LinearOperator.is_assemblable`. A domain
composite that emits (``CoupledOperator``'s block law, ``BulkLift``'s
embedding, ``FlattenedOperator``'s pass-through) is an atom to the free
algebra: it keeps its law in its own module and recurses through the
emitter it is handed (ruling 61).

Three consumption modes of one per-cell closure algebra (the
stencil-assembly campaign, Phase 2b):

===========  =========================================  =====================
Mode         Consumer                                   Walk
===========  =========================================  =====================
``solve``    sweep = forward/triangular substitution    cells in walk order
``apply``    Krylov matvec = action                     cells in any order
assembly     emit the SAME per-cell coefficients as     cells, scatter into
             ``(row, col, value)``                      a global sparse matrix
===========  =========================================  =====================

The matrix and its algorithms (ruling 60)
=========================================

An operator's matrix representation is fixed by its spaces' coordinates (the
axes in their order): :meth:`~orpheus.numerics.operator.LinearOperator.as_matrix`
realises it by probing, column :math:`j` being :math:`A e_j`, and works for
every operator with spaces. :func:`assemble` is a faster algorithm for the
same matrix, from the operator's structure instead of :math:`n`
applications; :func:`dense` densifies through it where it exists. The two
agree in exact arithmetic, which is what the probed-against-assembled gates
compare: two algorithms, one matrix.

Why a thin carrier and NOT a new emission type (2-P0 attacker ruling)
=====================================================================

A first-class COO-builder type with its own algebra would twin the operator
algebra one layer down: every law restated on triplet buffers. The honest
carrier is :mod:`scipy.sparse` itself (its CSR algebra IS the
Mat-category realisation), wrapped just thickly enough to conform to
:class:`~orpheus.numerics.operator.LinearOperator`, as
:class:`~orpheus.numerics.flat_operator.FlattenedOperator` serialises a
structure that already exists. The duplicate-summing FEM scatter is scipy's
own ``COO → CSR`` conversion (:meth:`~scipy.sparse.coo_array.tocsr` sums
repeated ``(row, col)`` entries), so an emitter may scatter per-cell and
per-face contributions freely.

The flat carrier contract
=========================

``apply`` acts on **bare 1-D flat vectors** (the
:class:`~orpheus.numerics.flat_operator.FlattenedOperator` serialisation
convention): the global DOF numbering is the operator family's existing flat
layout (``FullField.to_flat()``, bulk C-ravel ⊕ trace buffer, for composite
carriers; the bulk C-ravel alone for per-ordinate SN spatial blocks) and is
NEVER re-derived here. ``domain`` and ``codomain`` are passthrough metadata
from the emitting operator.

``apply_transpose`` is the exact CSR transpose (``is_adjointable`` True);
``is_invertible`` stays the base ``False`` (inversion goes through an
explicit :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`,
chosen by the consumer).
"""

from __future__ import annotations

from functools import reduce
from typing import TYPE_CHECKING, Callable, Optional, Protocol, cast

import numpy as np
from scipy import sparse

from orpheus.numerics.operator import (
    FusedAtom,
    LinearOperator,
    MatrixTooLarge,
    OperatorProduct,
    OperatorSum,
    _check_dimension,
    _left_fold,
    _resolve_basis_shape,
)

if TYPE_CHECKING:
    from typing import TypeGuard

    from orpheus.numerics.space import FunctionSpace


__all__ = [
    "MissingAssembly",
    "SparseAssembledOperator",
    "SupportsAssembly",
    "assemblable",
    "assemble",
    "dense",
]


class MissingAssembly(TypeError):
    r"""Asked for the sparse assembly of an operator that has none.

    The ASSEMBLY-axis refusal, sibling of
    :class:`~orpheus.numerics.operator.NotInvertible` (inverse axis) and
    :class:`~orpheus.numerics.operator.MissingAdjoint` (adjoint axis): raised
    by :func:`assemble` when some part of the operator emits no stencil. The
    message names that part's type and the law that needed it. An operator
    without a stencil keeps the probing
    :meth:`~orpheus.numerics.operator.LinearOperator.as_matrix`, its total
    dense representation.
    """


class SupportsAssembly(Protocol):
    r"""An operator that emits its own raw sparse matrix.

    ``assembled_matrix(emit)`` returns the operator's matrix in the flat
    layout of its spaces; ``emit`` is the functor's raw emission, handed to
    the domain composites that recurse into their parts (ruling 61). A true
    leaf ignores it. Never ``isinstance`` this: narrow with
    :func:`assemblable`.
    """

    @property
    def is_assemblable(self) -> bool: ...

    def assembled_matrix(
        self, emit: "Callable[[LinearOperator], sparse.csr_array]", /,
    ) -> "sparse.sparray | sparse.spmatrix": ...


def assemblable(op: "LinearOperator") -> "TypeGuard[SupportsAssembly]":
    r"""Narrow ``op`` to :class:`SupportsAssembly` iff its flag says it emits.

    The assembly-axis member of the bridge family
    (:func:`~orpheus.numerics.operator.invertible`,
    :func:`~orpheus.numerics.operator.adjointable`): the runtime check and the
    static permission to call ``assembled_matrix`` are one construct.
    """
    return op.is_assemblable


class SparseAssembledOperator(LinearOperator):
    r"""A structurally-assembled sparse matrix, conforming to the operator algebra.

    (Deliberately unparameterized — the ndarray-carrier family convention,
    the :class:`~orpheus.numerics.flat_operator.FlattenedOperator` /
    :class:`~orpheus.numerics.elementary_operators.DiagonalOperator` precedent: the
    carrier is the bare flat vector.)

    Parameters
    ----------
    matrix : scipy sparse matrix / array
        The assembled matrix in any scipy-sparse format — canonically the
        COO an emitter scatters ``(row, col, value)`` triplets into.
        Converted once to CSR at construction; the conversion SUMS
        duplicate entries (the FEM scatter semantics an emitter relies on
        for shared-face / shared-diagonal contributions).
    domain, codomain : FunctionSpace, optional
        Passthrough space metadata from the emitting operator (the flat
        ``shape`` identity of a composite space matches this matrix's
        column/row dimension; the
        :class:`~orpheus.numerics.operator.OperatorSum` /
        ``OperatorProduct`` guards read them).

    Notes
    -----
    ``apply(x)`` accepts any array raveling to ``matrix.shape[1]`` and
    returns the 1-D CSR matvec. The matrix is exposed as :attr:`matrix`
    (CSR) — consumers that need structure (the walk-order triangularity
    gates, the DSA moment reduction ``R·A·P``) read it directly.
    """

    def __init__(
        self,
        matrix: "sparse.sparray | sparse.spmatrix",
        *,
        domain: "Optional[FunctionSpace]" = None,
        codomain: "Optional[FunctionSpace]" = None,
    ) -> None:
        #: The assembled matrix, CSR-canonicalized (duplicates summed).
        self.matrix: "sparse.csr_array" = sparse.csr_array(matrix)
        # The ONE scipy-boundary seam: scipy 1.17's inline annotations
        # (no scipy-stubs in this environment) type ``sparray.shape`` as
        # ``None``, while the runtime value is the ``(rows, cols)``
        # tuple. Freeze it here, once, behind a documented cast; every
        # dimension read below consumes the well-typed frozen tuple.
        rows, cols = cast("tuple[int, int]", self.matrix.shape)
        self._shape: tuple[int, int] = (int(rows), int(cols))
        self._domain = domain
        self._codomain = codomain

    @property
    def domain(self) -> "Optional[FunctionSpace]":
        return self._domain

    @property
    def codomain(self) -> "Optional[FunctionSpace]":
        return self._codomain

    @property
    def shape(self) -> tuple[int, int]:
        """The assembled ``(rows, cols)`` dimensions."""
        return self._shape

    def apply(self, x: "np.ndarray", /) -> np.ndarray:
        r"""Return the CSR matvec :math:`M\,x` on a flat 1-D vector."""
        flat = np.asarray(x, dtype=float).ravel()
        if flat.size != self._shape[1]:
            raise ValueError(
                f"SparseAssembledOperator.apply: input size {flat.size} "
                f"does not match the assembled column dimension "
                f"{self._shape[1]}."
            )
        return np.asarray(self.matrix @ flat)

    def apply_transpose(self, x: "np.ndarray", /) -> np.ndarray:
        r"""Return :math:`M^{\mathsf T}\,x` — the exact CSR transpose matvec."""
        flat = np.asarray(x, dtype=float).ravel()
        if flat.size != self._shape[0]:
            raise ValueError(
                f"SparseAssembledOperator.apply_transpose: input size "
                f"{flat.size} does not match the assembled row dimension "
                f"{self._shape[0]}."
            )
        return np.asarray(self.matrix.T @ flat)

    @property
    def is_adjointable(self) -> bool:
        return True  # a materialized matrix's Euclidean transpose is free

    # is_invertible inherits the base ``False``: the sparse realization
    # reads values, not structure — inversion is an explicit consumer
    # choice (MatrixInverseOperator over it, or a sparse-LU sibling).

    # ── The assembly functor is closed on its carrier ─────────────────

    @property
    def is_assemblable(self) -> bool:
        return True  # already assembled — the functor's fixed point

    def assembled_matrix(self, emit: "Callable[[LinearOperator], sparse.csr_array]", /) -> "sparse.csr_array":
        r"""Return the held matrix; :func:`assemble` returns this carrier itself."""
        return self.matrix

    # ── Densification (the Mat-functor serialization boundary) ────────

    def as_matrix(
        self,
        *,
        basis_shape: tuple[int, ...] | None = None,
        max_dimension: int = 4096,
    ) -> np.ndarray:
        r"""Densify the assembled matrix — the structure-specific override.

        The dimension is already known from assembly, so the generic
        apply-to-basis loop collapses to one ``.toarray()`` (the
        :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
        override precedent). The base contract is honored: an explicit
        ``basis_shape`` must agree with the assembled column dimension,
        and ``max_dimension`` still gates the densification.
        """
        rows, cols = self._shape
        if basis_shape is not None and int(np.prod(basis_shape)) != cols:
            raise ValueError(
                f"as_matrix on SparseAssembledOperator: basis_shape "
                f"{basis_shape} (dimension {int(np.prod(basis_shape))}) "
                f"contradicts the assembled column dimension {cols}."
            )
        if cols > max_dimension:
            raise MatrixTooLarge(
                f"as_matrix on SparseAssembledOperator: dimension {cols} "
                f"exceeds max_dimension={max_dimension}. Densifying would "
                f"commit ~{8 * rows * cols / 1e6:.0f} MB; keep the "
                f"operator sparse (read .matrix directly) or raise the "
                f"gate."
            )
        return np.asarray(self.matrix.toarray())

    def __repr__(self) -> str:
        rows, cols = self._shape
        return (
            f"SparseAssembledOperator({rows}×{cols}, "
            f"nnz={self.matrix.nnz})"
        )


# ── The functor ──────────────────────────────────────────────────────────


def assemble(op: "LinearOperator") -> SparseAssembledOperator:
    r"""Return :math:`[A]`, the sparse matrix of ``op``, through the assembly homomorphism.

    An assembled operator is its own assembly. Otherwise the raw emission
    (:func:`_emit`) is wrapped once, with ``op``'s ends as the carrier's
    passthrough metadata.

    Raises
    ------
    MissingAssembly
        Some part of ``op`` emits no stencil; the message names it.
    """
    if isinstance(op, SparseAssembledOperator):
        return op
    return SparseAssembledOperator(_emit(op), domain=op.domain, codomain=op.codomain)


def _emit(op: "LinearOperator", law: str = "assemble") -> "sparse.csr_array":
    r"""The raw emission of ``op``: the composite laws of the free algebra, else ``op``'s own.

    A leaf's emission is canonicalised to CSR as it enters the laws, which
    sums its duplicate ``(row, col)`` entries before any coefficient or
    product touches them (the FEM scatter semantics, and the rounding the
    carrier fixed before ruling 24). ``law`` names the law that asked for
    this part, so a refusal deep in an expression says both what failed and
    what needed it.
    """
    match op:
        case OperatorSum():
            return _left_fold(
                op.terms, lambda w: _emit(w, "[Σ cᵢwᵢ] = Σ cᵢ[wᵢ]"),
            )
        case OperatorProduct():
            return reduce(
                lambda m, f: m @ _emit(f, "[A₁…Aₙ] = [A₁]…[Aₙ]"),
                op.factors[1:],
                _emit(op.factors[0], "[A₁…Aₙ] = [A₁]…[Aₙ]"),
            )
        case FusedAtom():
            return _emit(op.expanded, f"[{type(op).__name__}] = [its expansion]")
    own_law = f"{type(op).__name__}'s assembly"
    if not assemblable(op):
        # A domain composite's flag is its parts': run its own recursion so
        # the refusal names the part that cannot emit, not the composite.
        emission = getattr(op, "assembled_matrix", None)
        if emission is not None:
            emission(lambda part: _emit(part, own_law))
        raise MissingAssembly(
            f"{law}: {type(op).__name__} emits no assembly "
            f"(is_assemblable is False)."
        )
    return sparse.csr_array(op.assembled_matrix(lambda part: _emit(part, own_law)))


def dense(
    op: "LinearOperator",
    *,
    basis_shape: tuple[int, ...] | None = None,
    max_dimension: int = 4096,
) -> np.ndarray:
    r"""Return the dense matrix of ``op``, by assembly where ``op`` emits, else by probing.

    The same matrix as :meth:`~orpheus.numerics.operator.LinearOperator.as_matrix`
    (ruling 60), computed from the operator's structure when it has one: the
    choice of algorithm belongs to the consumer that pays for it
    (:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`).
    The gate contract is ``as_matrix``'s: an explicit ``basis_shape`` must
    agree with the assembled column dimension, and ``max_dimension`` bounds
    the densification.
    """
    if not op.is_assemblable:
        return op.as_matrix(basis_shape=basis_shape, max_dimension=max_dimension)
    shape = _resolve_basis_shape(op, basis_shape)
    _check_dimension(op, shape, max_dimension)
    return assemble(op).as_matrix(basis_shape=shape, max_dimension=max_dimension)

