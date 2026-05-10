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

from typing import TYPE_CHECKING, Optional, Protocol, runtime_checkable

import numpy as np
import scipy.sparse.linalg as spla

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace

__all__ = [
    "LinearOperator",
    "LinearOperatorMixin",
    "MissingCapability",
    "IncompatibleOperatorComposition",
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


class IncompatibleOperatorComposition(ValueError):
    """A composition's operands carry incompatible function spaces.

    Raised at composition time when two operators with declared
    :attr:`domain`/:attr:`range` carry shapes that cannot be combined
    (Sum: ``a.domain != b.domain`` or ``a.range != b.range``;
    Product ``A @ B``: ``A.domain != B.range``). The check is
    skipped when either operand has ``None`` for its domain or range
    — backward-compatible with operators predating Issue 9.6 that
    carry no function-space metadata.
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

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        """The function space this operator consumes, or ``None``.

        Operators that pre-date Issue 9.6 (and any operator that has
        no canonical function-space tagging) return ``None``. When
        either operand of a composition has ``None`` for ``domain``
        or ``range``, the composability check is skipped — preserving
        the legacy duck-typed behaviour for code paths that do not
        track spaces.
        """
        ...

    @property
    def range(self) -> Optional["FunctionSpace"]:
        """The function space this operator produces, or ``None``.

        See :attr:`domain` for the ``None`` semantics.
        """
        ...

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

    Issue 9.6 additions:

    * :meth:`__call__` — alias for :meth:`apply`, matching :math:`A(x)`
      math notation.
    * :meth:`__pow__` — repeated composition for non-negative integer
      powers (``A**0`` is the identity, ``A**n`` for ``n>=1`` is
      ``A @ A @ ... @ A``).
    * :meth:`adjoint` — weight-aware Hilbert adjoint for operators
      whose domain and range carry inner products.
    * :attr:`H` — property alias for ``adjoint()`` matching the Grand
      Report v3 §6.3 vocabulary.
    * Default :attr:`domain` / :attr:`range` return ``None`` —
      backward-compatible with operators predating Issue 9.6.
    * :meth:`__repr__` — uniform default reporting class name,
      domain/range, and capabilities.
    """

    capabilities: frozenset[str]

    # ------------------------------------------------------------------
    # Function-space tagging (defaults — concrete operators may override)
    # ------------------------------------------------------------------

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        return None

    @property
    def range(self) -> Optional["FunctionSpace"]:
        return None

    # ------------------------------------------------------------------
    # Algebra dunders
    # ------------------------------------------------------------------

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

    def __call__(self, *args, **kwargs):
        """Alias for :meth:`apply`. Lets user code write ``A(x)``.

        Accepts ``*args, **kwargs`` so multi-argument applies (e.g.
        :meth:`BoundaryOperator.apply` taking ``(psi_out, quadrature)``)
        compose ergonomically: ``bc(psi_out, quad)`` reads as math.
        """
        return self.apply(*args, **kwargs)  # type: ignore[attr-defined]

    def __pow__(self, n: int) -> "LinearOperator":
        r"""Return :math:`A^n` for non-negative integer ``n``.

        ``n == 0`` returns :class:`IdentityOperator`. ``n == 1``
        returns ``self`` unchanged. ``n >= 2`` builds the composition
        ``A @ A @ ... @ A`` via repeated :meth:`__matmul__`. Negative
        powers raise :class:`ValueError` — operator inverse construction
        is not part of this API; use the operator's :meth:`solve`
        capability directly when an inverse is needed.
        """
        if not isinstance(n, (int, np.integer)):
            return NotImplemented
        if n < 0:
            raise ValueError(
                "operator inverse not constructed via __pow__; "
                "consult the operator's solve() capability for inverse "
                "actions."
            )
        if n == 0:
            return IdentityOperator()
        if n == 1:
            return self  # type: ignore[return-value]
        result: "LinearOperator" = self  # type: ignore[assignment]
        for _ in range(n - 1):
            result = result @ self  # type: ignore[operator]
        return result

    # ------------------------------------------------------------------
    # Adjoint
    # ------------------------------------------------------------------

    def adjoint(self) -> "LinearOperator":
        r"""Return the Hilbert adjoint :math:`A^*`.

        For an operator :math:`A : V \to W` with diagonal inner-product
        weights :math:`w_V` (on the domain) and :math:`w_W` (on the
        range), the Hilbert adjoint satisfies

        .. math::

           \langle A x, y \rangle_W \;=\; \langle x, A^* y \rangle_V

        which gives :math:`A^* y = (1/w_V) \odot
        \mathrm{apply\_transpose}(w_W \odot y)`. When both weight
        arrays are ``None`` (Euclidean inner product on both sides)
        the adjoint reduces to the representation transpose.

        The returned wrapper preserves :meth:`apply` (= adjoint
        action) and swaps :attr:`domain` ↔ :attr:`range`. The
        :pydata:`capabilities` set advertises ``apply`` whenever the
        underlying operator advertises :pydata:`CAP_APPLY_TRANSPOSE`;
        otherwise the call to :meth:`apply` will raise at call time.
        """
        return _AdjointOperator(self)  # type: ignore[arg-type]

    @property
    def H(self) -> "LinearOperator":
        """Alias for :meth:`adjoint`. Matches the Grand Report v3 §6.3
        Hilbert-adjoint vocabulary (``A.H`` reads as "A dagger")."""
        return self.adjoint()

    # ------------------------------------------------------------------
    # Repr
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        cls = type(self).__name__
        d = getattr(self, "domain", None)
        r = getattr(self, "range", None)
        d_name = repr(d.name) if d is not None else "'?'"
        r_name = repr(r.name) if r is not None else "'?'"
        caps = sorted(getattr(self, "capabilities", frozenset()))
        return f"<{cls} domain={d_name} range={r_name} caps={caps}>"


# ---------------------------------------------------------------------------
# Adjoint wrapper
# ---------------------------------------------------------------------------


class _AdjointOperator(LinearOperatorMixin):
    """Hilbert-adjoint wrapper around a :class:`LinearOperator`.

    Constructed by :meth:`LinearOperatorMixin.adjoint` (and its alias
    ``A.H``). Domain/range are swapped relative to the inner operator;
    :meth:`apply` performs the weight-aware adjoint action.

    The capability set is derived from the inner operator: ``apply``
    on the adjoint maps to ``apply_transpose`` on the inner, so the
    adjoint advertises :pydata:`CAP_APPLY` iff the inner advertises
    :pydata:`CAP_APPLY_TRANSPOSE`. The reverse direction
    (apply_transpose on the adjoint = apply on the inner) is not
    needed by any current consumer in 9.6 and is deferred —
    :meth:`apply_transpose` raises :class:`NotImplementedError`
    until a consumer demands it.
    """

    def __init__(self, inner: "LinearOperator") -> None:
        self.inner = inner
        # Capability swap: adjoint can apply iff inner can apply_transpose.
        caps: set[str] = set()
        if _has(inner, CAP_APPLY_TRANSPOSE):
            caps.add(CAP_APPLY)
        # Adjoint can apply_transpose iff inner can apply — but defer
        # the reverse direction until a consumer requires it.
        # Solve generally does NOT propagate (the adjoint of A^{-1}
        # would need A.H.solve = (A.solve).H, additional algebra).
        self.capabilities = frozenset(caps)

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        # Adjoint of A: V → W is A.H: W → V — domain swaps with inner.range.
        return getattr(self.inner, "range", None)

    @property
    def range(self) -> Optional["FunctionSpace"]:
        return getattr(self.inner, "domain", None)

    def apply(self, y: np.ndarray) -> np.ndarray:
        if not _has(self.inner, CAP_APPLY_TRANSPOSE):
            raise MissingCapability(
                f"adjoint application requires {CAP_APPLY_TRANSPOSE!r} on "
                f"the inner operator; {type(self.inner).__name__} "
                f"advertises {getattr(self.inner, 'capabilities', frozenset())}."
            )
        # Hilbert-adjoint action:
        #   (A^* y)_V = (1/w_V) ⊙ apply_transpose(w_W ⊙ y)
        # On the adjoint wrapper, ``range`` is the inner operator's
        # domain (V) — so ``range.inner_product_weights`` is w_V; and
        # ``domain`` is the inner operator's range (W) — so
        # ``domain.inner_product_weights`` is w_W. Read with care.
        w_W = None
        w_V = None
        inner_range = getattr(self.inner, "range", None)
        inner_domain = getattr(self.inner, "domain", None)
        if inner_range is not None:
            w_W = inner_range.inner_product_weights
        if inner_domain is not None:
            w_V = inner_domain.inner_product_weights

        z = y if w_W is None else (w_W * y)
        result = self.inner.apply_transpose(z)  # type: ignore[attr-defined]
        if w_V is not None:
            result = result / w_V
        return result

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        raise NotImplementedError(
            "apply_transpose on an _AdjointOperator wrapper is not "
            "supported in 9.6; if a consumer needs it, take the adjoint "
            "of the original inner operator's transpose directly."
        )


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
        # Domain/range compatibility check (skipped when either
        # operand lacks function-space metadata — backward-compatible
        # with operators that pre-date Issue 9.6).
        a_dom, a_rng = getattr(a, "domain", None), getattr(a, "range", None)
        b_dom, b_rng = getattr(b, "domain", None), getattr(b, "range", None)
        if (
            a_dom is not None and b_dom is not None and a_dom != b_dom
        ):
            raise IncompatibleOperatorComposition(
                f"OperatorSum requires equal domains; got {a_dom!r} and "
                f"{b_dom!r}."
            )
        if (
            a_rng is not None and b_rng is not None and a_rng != b_rng
        ):
            raise IncompatibleOperatorComposition(
                f"OperatorSum requires equal ranges; got {a_rng!r} and "
                f"{b_rng!r}."
            )
        self.a = a
        self.b = b

        caps = {CAP_APPLY}
        if _has(a, CAP_APPLY_TRANSPOSE) and _has(b, CAP_APPLY_TRANSPOSE):
            caps.add(CAP_APPLY_TRANSPOSE)
        # solve does NOT propagate — see docstring.
        self.capabilities = frozenset(caps)

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        a_dom = getattr(self.a, "domain", None)
        return a_dom if a_dom is not None else getattr(self.b, "domain", None)

    @property
    def range(self) -> Optional["FunctionSpace"]:
        a_rng = getattr(self.a, "range", None)
        return a_rng if a_rng is not None else getattr(self.b, "range", None)

    def apply(self, x, *extra, **kwextra):
        # Forward extra positional/keyword args so BoundaryOperator-style
        # multi-argument apply signatures (e.g.
        # ``apply(psi_out, quadrature)``) compose under sums.
        return self.a.apply(x, *extra, **kwextra) + self.b.apply(
            x, *extra, **kwextra
        )

    def apply_transpose(self, x, *extra, **kwextra):
        return self.a.apply_transpose(x, *extra, **kwextra) + self.b.apply_transpose(  # type: ignore[attr-defined]
            x, *extra, **kwextra,
        )


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
        # Domain/range compatibility check for ``A @ B``: A.domain
        # must equal B.range. Skipped when either is None.
        a_dom = getattr(a, "domain", None)
        b_rng = getattr(b, "range", None)
        if (
            a_dom is not None and b_rng is not None and a_dom != b_rng
        ):
            raise IncompatibleOperatorComposition(
                f"OperatorProduct A @ B requires A.domain == B.range; "
                f"got A.domain={a_dom!r}, B.range={b_rng!r}."
            )
        self.a = a
        self.b = b

        caps = {CAP_APPLY}
        if _has(a, CAP_SOLVE) and _has(b, CAP_SOLVE):
            caps.add(CAP_SOLVE)
        if _has(a, CAP_APPLY_TRANSPOSE) and _has(b, CAP_APPLY_TRANSPOSE):
            caps.add(CAP_APPLY_TRANSPOSE)
        self.capabilities = frozenset(caps)

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        # A @ B: input space is B.domain.
        return getattr(self.b, "domain", None)

    @property
    def range(self) -> Optional["FunctionSpace"]:
        # A @ B: output space is A.range.
        return getattr(self.a, "range", None)

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

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        return getattr(self.op, "domain", None)

    @property
    def range(self) -> Optional["FunctionSpace"]:
        return getattr(self.op, "range", None)

    def apply(self, x, *extra, **kwextra):
        return self.scalar * self.op.apply(x, *extra, **kwextra)

    def solve(self, b_vec: np.ndarray) -> np.ndarray:
        # (αL)^{-1} = (1/α) L^{-1}
        return (1.0 / self.scalar) * self.op.solve(b_vec)  # type: ignore[attr-defined]

    def apply_transpose(self, x, *extra, **kwextra):
        return self.scalar * self.op.apply_transpose(x, *extra, **kwextra)  # type: ignore[attr-defined]


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
