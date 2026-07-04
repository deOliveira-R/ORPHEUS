r"""The sparse assembled realization of a structurally-emitted operator.

:class:`SparseAssembledOperator` is the CARRIER of the assembly mode — the
third consumption mode of the one per-cell closure algebra (stencil-assembly
campaign, Phase 2b):

===========  =========================================  =====================
Mode         Consumer                                   Walk
===========  =========================================  =====================
``solve``    sweep = forward/triangular substitution    cells in walk order
``apply``    Krylov matvec = action                     cells in any order
``assemble``  emit the SAME per-cell coefficients as     cells, scatter into
             ``(row, col, value)``                      a global sparse matrix
===========  =========================================  =====================

``assemble()`` is the same **additive-monoidal functor** ``Op → Mat`` that
:meth:`~orpheus.numerics.operator.LinearOperator.as_matrix` realizes by
apply-to-basis probing — landed in a *sparse* carrier by *structural
emission* instead of :math:`n` operator applications. Leaves that know their
stencil override ``assemble()``; composites recurse through the homomorphism
laws realized on the composers (Sum → ``+``, Product → ``@``, Scaled →
scalar ``*``); this class closes the functor (an assembled operator's
``assemble()`` is itself).

Why a thin wrapper and NOT a new emission type (2-P0 attacker ruling)
=====================================================================

A first-class COO-builder type with its own algebra would twin the operator
algebra one layer down — every law (`+`, `@`, scaling) restated on triplet
buffers. The honest carrier is :mod:`scipy.sparse` itself (its CSR algebra IS
the Mat-category realization), wrapped just thickly enough to conform to
:class:`~orpheus.numerics.operator.LinearOperator` — exactly the
:class:`~orpheus.numerics.flat_operator.FlattenedOperator` precedent: a
*serialization* of structure that already exists, not new structure. The
duplicate-summing FEM scatter is scipy's own ``COO → CSR`` conversion
(:meth:`~scipy.sparse.coo_array.tocsr` sums repeated ``(row, col)`` entries),
so an emitter may scatter per-cell / per-face contributions freely and the
carrier performs the assembly summation.

The flat carrier contract
=========================

``apply`` acts on **bare 1-D flat vectors** (the
:class:`~orpheus.numerics.flat_operator.FlattenedOperator` serialization
convention): the global DOF numbering is the operator family's existing flat
layout — ``FullField.to_flat()`` (bulk C-ravel ⊕ trace buffer) for composite
carriers, the bulk C-ravel alone for per-ordinate SN spatial blocks — and is
NEVER re-derived here (2-P0 attacker ruling: a reified ``LocalToGlobalMap``
is deferred until an unstructured consumer arrives). ``domain`` /
``codomain`` are passthrough metadata from the emitting operator.

Structural axes: ``apply_transpose`` is the exact CSR transpose (a
materialized matrix's Euclidean transpose is free and honest —
``is_adjointable`` True); ``is_invertible`` stays the base ``False`` (the
matrix realization reads values, not structure — inversion goes through an
explicit
:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
over it, or a future sparse-LU sibling, chosen by the consumer).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional, cast

import numpy as np
from scipy import sparse

from orpheus.numerics.operator import LinearOperator, MatrixTooLarge

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace


__all__ = ["SparseAssembledOperator"]


class SparseAssembledOperator(LinearOperator):
    r"""A structurally-assembled sparse matrix, conforming to the operator algebra.

    (Deliberately unparameterized — the ndarray-carrier family convention,
    the :class:`~orpheus.numerics.flat_operator.FlattenedOperator` /
    :class:`~orpheus.numerics.operator.DiagonalOperator` precedent: the
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

    def assemble(self) -> "SparseAssembledOperator":
        r"""Return this operator itself — assembly is idempotent."""
        return self

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
