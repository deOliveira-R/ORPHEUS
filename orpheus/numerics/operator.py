r"""Linear-operator algebra for matrix-free transport solvers.

The neutron transport eigenvalue problem and its fixed-source cousin
both reduce to compositions of a small set of linear operators acting
on a discrete flux distribution :math:`\psi`:

.. math::

    (L - S - F)\,\psi \;=\; q
    \qquad \text{(fixed source)}

.. math::

    (L - S)\,\psi \;=\; \tfrac{1}{k}\,F\,\psi
    \qquad \text{(eigenvalue)}

where :math:`L` is the streaming + collision operator, :math:`S` is
the scattering source operator, :math:`F` is the fission source
operator, and :math:`q` is an external source (Trefethen & Bau 1997,
§3.2 frame the matrix-free Krylov view). For an SN sweep, an MoC
ray-tracer, a CP collision-probability matrix, or a diffusion BiCGSTAB
solve, the *outer* algebra is identical even though the *implementation*
of each operator differs by orders of magnitude in cost and structure.

This module installs the **algebra** as runtime-checkable Protocols.
Any object providing ``apply(x) -> Lx`` participates; objects that can
also offer ``solve(b) -> L^{-1}b`` or ``apply_transpose(x) -> L^T x``
advertise those abilities through a :pydata:`capabilities` set.
Composers (:class:`OperatorSum`, :class:`OperatorProduct`,
:class:`ScaledOperator`, :class:`IdentityOperator`,
:class:`ZeroOperator`) compute their own capability set from
constituents — capability mismatches raise :class:`MissingCapability`
at composition time, NEVER at call time, so a downstream
:class:`scipy.sparse.linalg.LinearOperator` consumer never silently
hits a broken stub.

The choice of a capability *set* (rather than subclassing or abstract
methods) is deliberate. Many transport operators have no efficient
``solve``: the scattering source S has rank in the thousands and is
never inverted directly; the fission source F is rank-deficient (it
projects onto the fission spectrum). Forcing those classes to provide
``solve`` stubs that raise ``NotImplementedError`` is harmful — the
contract should be "I can do *these* things" not "I can do everything
or fail." See :ref:`operator-algebra` for the full design rationale.

The :func:`as_scipy_linop` adapter exposes any
:class:`LinearOperator` to scipy's Krylov solvers (BiCGSTAB, GMRES) so
the existing iterative-solver call sites in
:mod:`orpheus.sn.operator` and :mod:`orpheus.diffusion` can consume
operators built from this module without rewriting their inner loops.
"""

from __future__ import annotations

from typing import Protocol, runtime_checkable

import numpy as np
import scipy.sparse.linalg as spla

__all__ = [
    "LinearOperator",
    "LinearOperatorMixin",
    "MissingCapability",
    "OperatorSum",
    "OperatorProduct",
    "ScaledOperator",
    "IdentityOperator",
    "ZeroOperator",
    "as_scipy_linop",
    "CAP_APPLY",
    "CAP_SOLVE",
    "CAP_APPLY_TRANSPOSE",
]


# Capability tag literals. Strings (rather than an enum) so user
# operators can advertise method-specific tags without subclassing.
CAP_APPLY: str = "apply"
CAP_SOLVE: str = "solve"
CAP_APPLY_TRANSPOSE: str = "apply_transpose"


class MissingCapability(TypeError):
    """A composition would require a capability the constituents lack.

    Raised at composition time so that downstream Krylov / power-iteration
    consumers never hit a broken stub at call time. The exception message
    names the missing capability and the operand that lacks it.
    """


@runtime_checkable
class LinearOperator(Protocol):
    r"""Contract for a matrix-free linear operator on a flux vector.

    Any object exposing :meth:`apply` and a :pydata:`capabilities`
    frozenset participates. :meth:`solve` and :meth:`apply_transpose`
    are *optional* — their availability is advertised through the
    capability set rather than the type system, because many transport
    operators (S, F, and the BiCGSTAB-Jacobi-preconditioner family)
    have no efficient ``solve``.

    Composition operators (:class:`OperatorSum`, :class:`OperatorProduct`,
    :class:`ScaledOperator`) are wired through ``__add__``, ``__sub__``,
    ``__mul__`` (scalar), and ``__matmul__`` (operator product) so the
    typical algebra of the Boltzmann transport equation,
    :math:`(L - S - F)`, can be built with the natural Python syntax.
    The capability set of the composition is computed by the composer
    according to the closure laws documented in :ref:`operator-algebra`.

    Attributes
    ----------
    capabilities : frozenset[str]
        Subset of ``{"apply", "solve", "apply_transpose"}`` (extra tags
        are permitted for method-specific dispatch). The capability set
        is consulted by composers; it is the **single source of truth**
        for what an operator can do.

    Notes
    -----
    Shape and dtype are deliberately not part of the protocol. numpy
    duck-typing (broadcasting + dtype promotion) handles them at
    ``apply`` call time. Imposing a static shape would forbid operators
    whose action shape depends on the input (a multi-group transport
    sweep can output a different layout than its input vector). If a
    consumer needs shape information, it can probe ``op.apply(x)`` on a
    known-size probe vector once at setup.
    """

    capabilities: frozenset[str]

    def apply(self, x: np.ndarray) -> np.ndarray:
        r"""Return :math:`L\,x`.

        Mandatory. Every :class:`LinearOperator` must implement this.
        The :pydata:`capabilities` set MUST include
        :pydata:`CAP_APPLY` whenever this method is functional.
        """
        ...


def _has(op: object, cap: str) -> bool:
    """Return True iff ``op`` advertises capability ``cap``.

    Bare protocol-check that tolerates objects without the attribute
    (returns False). Used by composers to decide which capabilities
    propagate.
    """
    caps = getattr(op, "capabilities", None)
    return bool(caps) and cap in caps


# ───────────────────────────────────────────────────────────────────────
# Composition primitives
# ───────────────────────────────────────────────────────────────────────


class LinearOperatorMixin:
    """Mixin installing operator-algebra dunders on a :class:`LinearOperator`.

    Inherit this on any user-defined operator class to gain the natural
    Python algebra (``+``, ``-``, ``*`` for scalar multiplication, ``@``
    for operator composition). The composers (:class:`OperatorSum`,
    :class:`OperatorProduct`, :class:`ScaledOperator`,
    :class:`IdentityOperator`, :class:`ZeroOperator`) already inherit
    it, so the algebra is closed under further composition without any
    extra effort from user code.

    The mixin defines no state; it relies on the user class providing
    an :pydata:`apply` method and a :pydata:`capabilities` attribute to
    satisfy the :class:`LinearOperator` Protocol. The ``+``/``-``/``*``/
    ``@`` dunders simply delegate to the composer constructors, which
    enforce the capability closure laws.
    """

    capabilities: frozenset[str]

    def __add__(self, other: "LinearOperator") -> "OperatorSum":
        return OperatorSum(self, other)  # type: ignore[arg-type]

    def __radd__(self, other: "LinearOperator") -> "OperatorSum":
        return OperatorSum(other, self)  # type: ignore[arg-type]

    def __sub__(self, other: "LinearOperator") -> "OperatorSum":
        return OperatorSum(self, ScaledOperator(-1.0, other))  # type: ignore[arg-type]

    def __rsub__(self, other: "LinearOperator") -> "OperatorSum":
        return OperatorSum(other, ScaledOperator(-1.0, self))  # type: ignore[arg-type]

    def __mul__(self, other: float) -> "ScaledOperator":
        if not isinstance(other, (int, float, np.floating, np.integer)):
            return NotImplemented
        return ScaledOperator(float(other), self)  # type: ignore[arg-type]

    def __rmul__(self, other: float) -> "ScaledOperator":
        return self.__mul__(other)

    def __matmul__(self, other: "LinearOperator") -> "OperatorProduct":
        return OperatorProduct(self, other)  # type: ignore[arg-type]


class OperatorSum(LinearOperatorMixin):
    r"""Sum of two linear operators: :math:`(A + B)\,x = A\,x + B\,x`.

    Capability closure laws:

    * ``apply`` propagates iff **both** operands have ``apply``
      (the action is well-defined only when both summands act).
    * ``solve`` does **not** propagate. There is no general algorithm
      for inverting a sum :math:`(A + B)^{-1}` from the inverses of
      the operands — Sherman-Morrison-Woodbury applies only under
      low-rank structure.
    * ``apply_transpose`` propagates iff **both** operands have it
      (transposition distributes over sums: :math:`(A + B)^T = A^T + B^T`).

    Raises
    ------
    MissingCapability
        If either operand lacks :pydata:`CAP_APPLY` at construction
        time. Catch the failure here, not at the first ``apply`` call,
        so downstream Krylov consumers don't see a stub failure
        mid-iteration.
    """

    def __init__(self, a: LinearOperator, b: LinearOperator) -> None:
        if not _has(a, CAP_APPLY):
            raise MissingCapability(
                f"OperatorSum requires apply on both operands; left "
                f"operand {type(a).__name__} lacks {CAP_APPLY!r}."
            )
        if not _has(b, CAP_APPLY):
            raise MissingCapability(
                f"OperatorSum requires apply on both operands; right "
                f"operand {type(b).__name__} lacks {CAP_APPLY!r}."
            )
        self.a = a
        self.b = b

        caps = {CAP_APPLY}
        if _has(a, CAP_APPLY_TRANSPOSE) and _has(b, CAP_APPLY_TRANSPOSE):
            caps.add(CAP_APPLY_TRANSPOSE)
        # solve does NOT propagate — see docstring.
        self.capabilities = frozenset(caps)

    def apply(self, x: np.ndarray) -> np.ndarray:
        return self.a.apply(x) + self.b.apply(x)

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        return self.a.apply_transpose(x) + self.b.apply_transpose(x)  # type: ignore[attr-defined]


class OperatorProduct(LinearOperatorMixin):
    r"""Composition of two linear operators: :math:`(A\,B)\,x = A(B\,x)`.

    Capability closure laws:

    * ``apply`` propagates iff **both** operands have ``apply``
      (function composition).
    * ``solve`` propagates iff **both** operands have ``solve``, with
      the order reversed: :math:`(A\,B)^{-1} = B^{-1}\,A^{-1}`. We
      apply ``B.solve`` first, then ``A.solve``.
    * ``apply_transpose`` propagates iff **both** operands have it,
      with the order reversed: :math:`(A\,B)^T = B^T\,A^T`.

    Raises
    ------
    MissingCapability
        If either operand lacks :pydata:`CAP_APPLY` at construction.
    """

    def __init__(self, a: LinearOperator, b: LinearOperator) -> None:
        if not _has(a, CAP_APPLY):
            raise MissingCapability(
                f"OperatorProduct requires apply on both operands; "
                f"left operand {type(a).__name__} lacks {CAP_APPLY!r}."
            )
        if not _has(b, CAP_APPLY):
            raise MissingCapability(
                f"OperatorProduct requires apply on both operands; "
                f"right operand {type(b).__name__} lacks {CAP_APPLY!r}."
            )
        self.a = a
        self.b = b

        caps = {CAP_APPLY}
        if _has(a, CAP_SOLVE) and _has(b, CAP_SOLVE):
            caps.add(CAP_SOLVE)
        if _has(a, CAP_APPLY_TRANSPOSE) and _has(b, CAP_APPLY_TRANSPOSE):
            caps.add(CAP_APPLY_TRANSPOSE)
        self.capabilities = frozenset(caps)

    def apply(self, x: np.ndarray) -> np.ndarray:
        return self.a.apply(self.b.apply(x))

    def solve(self, b_vec: np.ndarray) -> np.ndarray:
        # (AB)^{-1} = B^{-1} A^{-1}
        return self.b.solve(self.a.solve(b_vec))  # type: ignore[attr-defined]

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        # (AB)^T = B^T A^T
        return self.b.apply_transpose(self.a.apply_transpose(x))  # type: ignore[attr-defined]


class ScaledOperator(LinearOperatorMixin):
    r"""Scalar multiple of a linear operator: :math:`(\alpha L)\,x = \alpha\,(L\,x)`.

    Scalar multiplication is a unitary operation on the capability
    set: every capability of the underlying operator is preserved.
    For ``solve``, :math:`(\alpha L)^{-1} = (1/\alpha)\,L^{-1}` so we
    divide on the way out (and the scalar must be non-zero — caught
    at composition time).
    """

    def __init__(self, scalar: float, op: LinearOperator) -> None:
        if not _has(op, CAP_APPLY):
            raise MissingCapability(
                f"ScaledOperator requires apply on its operand; "
                f"{type(op).__name__} lacks {CAP_APPLY!r}."
            )
        if scalar == 0.0:
            # Zero scaling is a degenerate case: the result has the
            # same capabilities as ZeroOperator, not as the underlying
            # operator. The user should construct ZeroOperator
            # explicitly to make the intent clear.
            raise ValueError(
                "ScaledOperator with zero scalar is degenerate; "
                "use ZeroOperator explicitly."
            )
        self.scalar = float(scalar)
        self.op = op
        # All capabilities of the underlying operator survive.
        # Filter to the recognised set; users may have method-specific
        # tags that we don't know how to forward.
        survivors = {CAP_APPLY}
        if _has(op, CAP_SOLVE):
            survivors.add(CAP_SOLVE)
        if _has(op, CAP_APPLY_TRANSPOSE):
            survivors.add(CAP_APPLY_TRANSPOSE)
        self.capabilities = frozenset(survivors)

    def apply(self, x: np.ndarray) -> np.ndarray:
        return self.scalar * self.op.apply(x)

    def solve(self, b_vec: np.ndarray) -> np.ndarray:
        # (αL)^{-1} = (1/α) L^{-1}
        return (1.0 / self.scalar) * self.op.solve(b_vec)  # type: ignore[attr-defined]

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        return self.scalar * self.op.apply_transpose(x)  # type: ignore[attr-defined]


class IdentityOperator(LinearOperatorMixin):
    r"""The identity operator :math:`I\,x = x`.

    All three primitive capabilities are trivially supported:
    :math:`I^{-1} = I` and :math:`I^T = I`. ``solve`` is the same
    code path as ``apply`` — it is just relabelled.
    """

    capabilities: frozenset[str] = frozenset(
        {CAP_APPLY, CAP_SOLVE, CAP_APPLY_TRANSPOSE}
    )

    def apply(self, x: np.ndarray) -> np.ndarray:
        return x

    def solve(self, b_vec: np.ndarray) -> np.ndarray:
        return b_vec

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        return x


class ZeroOperator(LinearOperatorMixin):
    r"""The zero operator :math:`0\,x = 0`.

    Has ``apply`` (returns ``np.zeros_like(x)``) and
    ``apply_transpose`` (also zero), but **not** ``solve``: the zero
    operator is not invertible. Forcing it to advertise ``solve``
    would be the harmful-stub anti-pattern this whole module is
    designed against — composers downstream that ask for ``solve``
    on :math:`L = A + 0` would silently get a meaningless answer.
    """

    capabilities: frozenset[str] = frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})

    def apply(self, x: np.ndarray) -> np.ndarray:
        return np.zeros_like(x)

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        return np.zeros_like(x)


# ───────────────────────────────────────────────────────────────────────
# scipy interop
# ───────────────────────────────────────────────────────────────────────


def as_scipy_linop(
    op: LinearOperator,
    shape: tuple[int, int],
    dtype: type | np.dtype = float,
) -> spla.LinearOperator:
    r"""Wrap an ORPHEUS :class:`LinearOperator` for scipy Krylov solvers.

    scipy's BiCGSTAB / GMRES / etc. consume
    :class:`scipy.sparse.linalg.LinearOperator` and call ``matvec``
    (and optionally ``rmatvec``). This adapter delegates to
    :meth:`LinearOperator.apply` (and :meth:`apply_transpose` if the
    capability is advertised), keeping the capability set as the
    single source of truth.

    Parameters
    ----------
    op
        Source operator. Must advertise :pydata:`CAP_APPLY`.
    shape
        ``(n_rows, n_cols)`` of the equivalent dense matrix. scipy
        requires a static shape; the ORPHEUS protocol does not.
        Picking the wrong shape here is an error scipy will detect at
        first ``matvec``.
    dtype
        scipy LinearOperator dtype — passed through unchanged.

    Returns
    -------
    scipy.sparse.linalg.LinearOperator
        Wrapper whose ``matvec`` is ``op.apply`` and (when available)
        ``rmatvec`` is ``op.apply_transpose``.

    Raises
    ------
    MissingCapability
        If ``op`` does not advertise :pydata:`CAP_APPLY`.
    """
    if not _has(op, CAP_APPLY):
        raise MissingCapability(
            f"as_scipy_linop requires {CAP_APPLY!r}; "
            f"{type(op).__name__} advertises "
            f"{getattr(op, 'capabilities', frozenset())}."
        )

    matvec = op.apply
    rmatvec = (
        op.apply_transpose  # type: ignore[attr-defined]
        if _has(op, CAP_APPLY_TRANSPOSE)
        else None
    )
    return spla.LinearOperator(shape, matvec=matvec, rmatvec=rmatvec, dtype=dtype)
