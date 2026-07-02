r"""The dense direct inverse :math:`A^{-1}` of a small materializable operator.

The fourth member of the #226 inverse family (taxonomy §12 step 5).
``A.inverse()`` is a STRUCTURE-KEYED factory — the concrete subclass is
the mathematical OBJECT the inverse is, never the algorithm's name
(taxonomy §13): a schedule-triangular forward returns the
direct-substitution :class:`~orpheus.sn.operators.sweep_operator.SweepOperator`;
a general sum returns the preconditioned-splitting
:class:`~orpheus.numerics.green_operator.GreenOperator`; a value-bearing
leaf returns the pointwise :class:`~orpheus.numerics.operator.InverseOperator`.
This class is the inverse of a **structureless small** operator — the 0-D
energy spectrum, a CP collision-probability block, any composition whose
only exploitable property is that it FITS: materialize
:math:`[A]` (:meth:`~orpheus.numerics.operator.LinearOperator.as_matrix`),
factor it once (LU), and every :meth:`apply` is a direct backsolve.

**The name is earned** (taxonomy §13, M-row):

* **M-materialise** — ``Ainv.as_matrix() @ A.as_matrix() ≈ I`` BOTH ways
  **at MACHINE·cond precision**. The precision GRAIN is the
  distinguisher, not the method's existence (spec §27.A): ``as_matrix``
  is a universal base method, so an iterative inverse also satisfies the
  identity — but only to its driver tolerance (each of its columns is an
  iterative solve), while this class's :meth:`as_matrix` override is one
  batched backsolve on the identity against the SAME LU factorization —
  no second realization, no iteration floor.
* **M-direct** — the true residual :math:`\lVert A(A^{-1}q) - q\rVert \le
  K\,\varepsilon_{\mathrm{mach}}\,\kappa(A)\,\lVert q\rVert` at
  machine-times-condition grain (a DIRECT solve, not an iteration
  tolerance), and the result is bit-identical under any
  ``initial_guess`` — an exact inverse has nothing to seed, so the
  canonical family keyword is accepted and ignored (contrast the sweep's
  Carlson closure and the Green's splitting start, which consume it).

**Values, not structure — the guard difference.** The constructor
deliberately does NOT require ``inner.is_invertible``. That predicate is
STRUCTURAL — it reads the operand tree (a sum is preconditionable only at
the canonical spelling with an invertible LEADING term). The matrix
realization reads VALUES: ``(−S) + (L+C)`` — which
:class:`~orpheus.numerics.green_operator.GreenOperator` refuses at
construction because its left-spine head is not invertible — materializes
to a perfectly invertible matrix, and THIS class inverts it. That is the
taxonomy §3 ``strategy=`` override seam realized honestly: not a flag on
``.inverse()`` but an explicit construction by a consumer who knows the
problem is small (the type IS the strategy choice). The guards that
remain are the honest ones — the size gate
(:class:`~orpheus.numerics.operator.MatrixTooLarge` propagates from the
eager materialization), squareness (a rectangular operator has no
two-sided inverse), and exact singularity (a zero LU pivot raises
:class:`numpy.linalg.LinAlgError` at CONSTRUCTION — this module family's
composition-time-not-call-time principle; NEAR-singularity is not
refused, it is priced into the M-direct :math:`\kappa(A)` bound).

**Consumers.** :func:`~orpheus.numerics.eigenvalue.direct_eigenvalue`'s
dense ``np.linalg.solve(A, F)`` IS this operator's action written as free
functions — the latent consumer that made the promotion finishable
(taxonomy §5); the engine itself stays ndarray-pure (its closed-form
verification is the point). The homogeneous solver materializes through
``as_matrix`` today and its full operator spelling
``K = MatrixInverseOperator(loss) @ F`` is the follow-on step (task
#138); CP's ``[P]`` — dense by construction, §14b — is the production
method that earns the class. Until then the M-gates construct it
directly, the same consumer ruling GreenOperator shipped under at step 4.

**Placement.** A leaf module importing :mod:`orpheus.numerics.operator`
(the ``green_operator.py`` precedent) — and unlike Green, nothing in
``operator.py`` routes back here: no automatic ``.inverse()`` returns
this type, so there is no late-import seam at all.
"""
from __future__ import annotations

import warnings
from typing import Any

import numpy as np
from numpy.linalg import LinAlgError
from scipy.linalg import LinAlgWarning, lu_factor, lu_solve

from orpheus.numerics.operator import (
    InverseWrapMixin,
    LinearOperator,
    MatrixTooLarge,
    _resolve_basis_shape,
)

__all__ = ["MatrixInverseOperator"]


class MatrixInverseOperator(InverseWrapMixin["LinearOperator"], LinearOperator):
    r"""The inverse operator :math:`A^{-1}` of a small materializable ``A``,
    realized as a one-time LU factorization of ``inner.as_matrix()`` plus a
    direct backsolve per :meth:`apply`.

    Construction is EAGER: the materialization IS the guard. Size
    (:class:`~orpheus.numerics.operator.MatrixTooLarge`), squareness
    (``ValueError``), and exact singularity
    (:class:`numpy.linalg.LinAlgError`) all surface here, never at apply
    time. ``inner.is_invertible`` is deliberately not consulted — the
    matrix realization reads values, not structure (see the module
    docstring for the ``(−S)+(L+C)`` witness).

    The wrap-delegate back-half — (``apply`` inverts,
    ``solve`` un-inverts), the domain↔codomain swap, ``solve`` = the
    forward matvec ``inner.apply``, and the object-identity involution
    ``inverse() → inner`` — is inherited from
    :class:`~orpheus.numerics.operator.InverseWrapMixin`. This class
    keeps the family's three per-sibling pieces: its ctor guard, its
    :meth:`apply` body, and ``__repr__``.

    Parameters
    ----------
    inner : LinearOperator
        The forward operator to invert. Anything with a working
        :meth:`~orpheus.numerics.operator.LinearOperator.as_matrix` —
        which is any apply-bearing operator below the size gate.
    basis_shape : tuple[int, ...], optional
        The element shape ``inner.apply`` consumes; resolution shares
        :meth:`~orpheus.numerics.operator.LinearOperator.as_matrix`'s
        single source (explicit wins, else ``inner.domain.shape``, else
        an informative ``ValueError``). Solutions are reshaped back to
        this carrier shape.
    max_dimension : int, optional
        The materialization size gate, forwarded to ``as_matrix``
        (default ``4096`` — see
        :class:`~orpheus.numerics.operator.MatrixTooLarge`).
    """

    def __init__(
        self,
        inner: "LinearOperator",
        *,
        basis_shape: tuple[int, ...] | None = None,
        max_dimension: int = 4096,
    ) -> None:
        shape = _resolve_basis_shape(inner, basis_shape)
        matrix = inner.as_matrix(basis_shape=shape, max_dimension=max_dimension)
        rows, cols = matrix.shape
        if rows != cols:
            raise ValueError(
                f"MatrixInverseOperator requires a SQUARE materialization; "
                f"{type(inner).__name__} materialized as {rows}×{cols} "
                f"(basis shape {shape}). A rectangular operator has no "
                f"two-sided inverse — M-materialise is unsatisfiable."
            )
        self._basis_shape = shape
        # One factorization for the operator's lifetime; scipy's own
        # singularity signal is a WARNING (getrf info > 0), which would
        # let a zero pivot flow into inf/nan backsolves — silence it and
        # raise the loud construction-time error from the U diagonal
        # instead (exactly the condition scipy warns on). Cardinal Rule 1:
        # fail at construction, never return a non-inverse.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", LinAlgWarning)
            # ``matrix`` is a fresh ctor-local from as_matrix — factor it
            # in place (1×n² transient instead of 2× at the gate ceiling).
            self._lu = lu_factor(matrix, overwrite_a=True)
        if np.any(np.diag(self._lu[0]) == 0.0):
            raise LinAlgError(
                f"MatrixInverseOperator: the materialization of "
                f"{type(inner).__name__} is exactly singular (zero LU "
                f"pivot) — the operator has no inverse. (Near-singular "
                f"operators are NOT refused; their conditioning is priced "
                f"into the M-direct residual bound.)"
            )
        super().__init__(inner)

    def apply(self, x: Any, /, *, initial_guess: Any | None = None) -> Any:
        r"""Return :math:`A^{-1}x` — one LU backsolve against the stored factors.

        ``initial_guess`` is the inverse family's canonical driver
        signature (taxonomy §12 step 3); an exact direct inverse has
        nothing to seed, so it is accepted and unused — which is
        precisely the M-direct seed-independence invariant.
        """
        del initial_guess  # exact direct inverse — M-direct IS seed-independence
        solution = lu_solve(self._lu, np.asarray(x).ravel())
        return solution.reshape(self._basis_shape)

    def as_matrix(
        self,
        *,
        basis_shape: tuple[int, ...] | None = None,
        max_dimension: int = 4096,
    ) -> np.ndarray:
        r"""Materialize :math:`[A^{-1}]` — one batched backsolve on the identity.

        The structure-specific override (taxonomy §11.4): the dimension
        is already known from construction and the LU factors already
        held, so the generic apply-to-basis loop collapses to
        ``lu_solve(lu, I)`` — the same math, batched in one LAPACK call.
        The base contract is honored: an explicit ``basis_shape`` must
        agree with the materialized dimension, and ``max_dimension``
        still gates (the caller may pass a tighter budget than the
        constructor's).
        """
        n = self._lu[0].shape[0]
        if basis_shape is not None and int(np.prod(basis_shape)) != n:
            raise ValueError(
                f"as_matrix on MatrixInverseOperator: basis_shape "
                f"{basis_shape} (dimension {int(np.prod(basis_shape))}) "
                f"contradicts the materialized dimension {n}."
            )
        if n > max_dimension:
            raise MatrixTooLarge(
                f"as_matrix on MatrixInverseOperator: dimension {n} "
                f"exceeds max_dimension={max_dimension}."
            )
        return lu_solve(self._lu, np.eye(n))

    def __repr__(self) -> str:
        n = self._lu[0].shape[0]
        return f"MatrixInverseOperator({self.inner!r}, n={n})"
