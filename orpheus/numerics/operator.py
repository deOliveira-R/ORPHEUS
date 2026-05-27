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
    "PermutationOperator",
    "IncomingOrdinateMaskTensor",
    "PeriodicWrapOperator",
    "DiagonalOperator",
    "TensorProductOperator",
    "SumOfTensorProductsOperator",
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
    :attr:`domain`/:attr:`codomain` carry shapes that cannot be combined
    (Sum: ``a.domain != b.domain`` or ``a.codomain != b.codomain``;
    Product ``A @ B``: ``A.domain != B.codomain``). The check is
    skipped when either operand has ``None`` for its domain or codomain
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
        or ``codomain``, the composability check is skipped — preserving
        the legacy duck-typed behaviour for code paths that do not
        track spaces.
        """
        ...

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
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


def _broadcast_for_leading_axes(w: Optional[np.ndarray], target_ndim: int) -> Optional[np.ndarray]:
    r"""Reshape ``w`` so it broadcasts against the leading axes of a tensor of ``target_ndim``.

    When the domain / codomain of an operator carries an inner-product
    weight tensor of shape ``w.shape == (a₀, …, a_{k-1})``, and the
    data tensor crossing the operator has shape ``(a₀, …, a_{k-1}, b₀,
    …, b_{m-1})`` (i.e., the metric's axes are the LEADING axes of the
    data), numpy's right-aligned broadcasting does NOT make ``w * data``
    work directly. Padding ``w`` with trailing 1s — shape ``(a₀, …,
    a_{k-1}, 1, …, 1)`` — restores correct broadcasting.

    Called by :meth:`_AdjointOperator.apply` so the metrics defined on
    function spaces (e.g.
    :class:`~orpheus.numerics.spaces.SphericalHarmonicSpace`'s
    ``inner_product_weights`` of shape ``(L+1, 2L+1)``) broadcast
    correctly against arbitrarily-shaped data tensors (e.g. moment
    fields of shape ``(L+1, 2L+1, ng, nx, ny)`` in the production SN
    layout).

    Parameters
    ----------
    w : Optional[np.ndarray]
        The metric tensor, or ``None`` (when the space uses the
        Euclidean inner product).
    target_ndim : int
        The number of dimensions of the data tensor the metric must
        multiply.

    Returns
    -------
    Optional[np.ndarray]
        ``None`` if ``w`` is ``None``; otherwise ``w`` reshaped to
        ``w.shape + (1,) * (target_ndim - w.ndim)`` when ``w.ndim <
        target_ndim``; ``w`` unchanged otherwise (already aligned).
    """
    if w is None:
        return None
    if w.ndim >= target_ndim:
        return w
    new_shape = w.shape + (1,) * (target_ndim - w.ndim)
    return w.reshape(new_shape)


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
      whose domain and codomain carry inner products.
    * :attr:`H` — property alias for ``adjoint()`` matching the Grand
      Report v3 §6.3 vocabulary.
    * Default :attr:`domain` / :attr:`codomain` return ``None`` —
      backward-compatible with operators predating Issue 9.6.
    * :meth:`__repr__` — uniform default reporting class name,
      domain/codomain, and capabilities.
    """

    capabilities: frozenset[str]

    # ------------------------------------------------------------------
    # Function-space tagging (defaults — concrete operators may override)
    # ------------------------------------------------------------------

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        return None

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
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

    def __neg__(self) -> "ScaledOperator":
        r"""Unary minus: return :math:`-A` as ``ScaledOperator(-1.0, A)``.

        Pythonic completion of the ``__sub__`` family — when ``A - B``
        works (which ``__sub__`` already provides via the ``A +
        ScaledOperator(-1.0, B)`` rewrite) Python's arithmetic
        convention is that ``-A`` should also work. Useful for
        adjoint-flux sensitivity rewrites (the streaming term flips
        sign under the adjoint), source-iteration residual
        corrections (``correction = -L @ delta``), and Jacobi-style
        splitting (``A = D - (L + U)``, where ``-(L + U)`` is the
        off-diagonal coupling).
        """
        return ScaledOperator(-1.0, self)  # type: ignore[arg-type]

    def __truediv__(self, scalar: float) -> "ScaledOperator":
        r"""Scalar division: ``A / α`` is ``(1/α) * A``.

        Used for normalisation in eigenvalue / Krylov iterations
        (``fission_normalised = F / k_eff``), homogenisation
        averages (``avg_streaming = sum_streaming / N``), and any
        consumer that reads more naturally with ``/`` than with the
        reciprocal-multiply form ``(1.0 / α) * A``.

        Raises :class:`TypeError` if ``scalar`` is not numeric.
        Division by zero raises :class:`ZeroDivisionError` per the
        standard Python convention (handled by ``1.0 / scalar``).
        """
        if not isinstance(scalar, (int, float, np.floating, np.integer)):
            return NotImplemented
        return ScaledOperator(1.0 / float(scalar), self)  # type: ignore[arg-type]

    def __matmul__(self, other: "LinearOperator") -> "OperatorProduct":
        return OperatorProduct(self, other)  # type: ignore[arg-type]

    def __and__(self, other: "LinearOperator") -> "TensorProductOperator":
        r"""Return :math:`A \otimes B` — the per-axis tensor-product operator.

        Per Grand Report v3 §6.3 line 721 and §15.1 line 2044. For two
        operators acting on independent tensor axes, ``A & B`` produces
        the operator whose action is "apply A on its axis, apply B on
        its axis" (sequentially; commutative because axes are disjoint).

        If either operand is already a :class:`TensorProductOperator`,
        the result is flattened so ``(A & B) & C`` and ``A & (B & C)``
        produce the same instance ``TensorProductOperator((A, B, C))``.
        """
        return TensorProductOperator._build(self, other)  # type: ignore[arg-type]

    def __rand__(self, other: "LinearOperator") -> "TensorProductOperator":
        return TensorProductOperator._build(other, self)  # type: ignore[arg-type]

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
        codomain), the Hilbert adjoint satisfies

        .. math::

           \langle A x, y \rangle_W \;=\; \langle x, A^* y \rangle_V

        which gives :math:`A^* y = (1/w_V) \odot
        \mathrm{apply\_transpose}(w_W \odot y)`. When both weight
        arrays are ``None`` (Euclidean inner product on both sides)
        the adjoint reduces to the representation transpose.

        The returned wrapper preserves :meth:`apply` (= adjoint
        action) and swaps :attr:`domain` ↔ :attr:`codomain`. The
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
        c = getattr(self, "codomain", None)
        d_name = repr(d.name) if d is not None else "'?'"
        c_name = repr(c.name) if c is not None else "'?'"
        caps = sorted(getattr(self, "capabilities", frozenset()))
        return f"<{cls} domain={d_name} codomain={c_name} caps={caps}>"


# ---------------------------------------------------------------------------
# Adjoint wrapper
# ---------------------------------------------------------------------------


class _AdjointOperator(LinearOperatorMixin):
    """Hilbert-adjoint wrapper around a :class:`LinearOperator`.

    Constructed by :meth:`LinearOperatorMixin.adjoint` (and its alias
    ``A.H``). Domain/codomain are swapped relative to the inner operator;
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
        # Adjoint of A: V → W is A.H: W → V — domain swaps with inner.codomain.
        return getattr(self.inner, "codomain", None)

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
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
        # On the adjoint wrapper, ``codomain`` is the inner operator's
        # domain (V) — so ``codomain.inner_product_weights`` is w_V; and
        # ``domain`` is the inner operator's codomain (W) — so
        # ``domain.inner_product_weights`` is w_W. Read with care.
        w_W = None
        w_V = None
        inner_codomain = getattr(self.inner, "codomain", None)
        inner_domain = getattr(self.inner, "domain", None)
        if inner_codomain is not None:
            w_W = inner_codomain.inner_product_weights
        if inner_domain is not None:
            w_V = inner_domain.inner_product_weights

        # Leading-axis broadcast: if the metric has fewer dims than the
        # input array, pad it with trailing 1s so it broadcasts against
        # the input's leading axes. This is the "space carries the
        # metric shaped for axis-0 broadcast; _AdjointOperator
        # broadcasts trailing" pattern documented in the moment-space
        # + layering plan §P1.4 — applies generically to any operator
        # whose domain/codomain metric is shape-aligned with the leading
        # axes of the data tensor (e.g., MomentProjection's
        # SphericalHarmonicSpace metric on the (ℓ, m) axes).
        w_W_b = _broadcast_for_leading_axes(w_W, y.ndim) if w_W is not None else None
        z = y if w_W_b is None else (w_W_b * y)
        result = self.inner.apply_transpose(z)  # type: ignore[attr-defined]
        if w_V is not None:
            w_V_b = _broadcast_for_leading_axes(w_V, result.ndim)
            result = result / w_V_b
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
        # Domain/codomain compatibility check (skipped when either
        # operand lacks function-space metadata — backward-compatible
        # with operators that pre-date Issue 9.6).
        a_dom, a_cod = getattr(a, "domain", None), getattr(a, "codomain", None)
        b_dom, b_cod = getattr(b, "domain", None), getattr(b, "codomain", None)
        if (
            a_dom is not None and b_dom is not None and a_dom != b_dom
        ):
            raise IncompatibleOperatorComposition(
                f"OperatorSum requires equal domains; got {a_dom!r} and "
                f"{b_dom!r}."
            )
        if (
            a_cod is not None and b_cod is not None and a_cod != b_cod
        ):
            raise IncompatibleOperatorComposition(
                f"OperatorSum requires equal codomains; got {a_cod!r} and "
                f"{b_cod!r}."
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
    def codomain(self) -> Optional["FunctionSpace"]:
        a_cod = getattr(self.a, "codomain", None)
        return a_cod if a_cod is not None else getattr(self.b, "codomain", None)

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
        # Domain/codomain compatibility check for ``A @ B``: A.domain
        # must equal B.codomain. Skipped when either is None.
        a_dom = getattr(a, "domain", None)
        b_cod = getattr(b, "codomain", None)
        if (
            a_dom is not None and b_cod is not None and a_dom != b_cod
        ):
            raise IncompatibleOperatorComposition(
                f"OperatorProduct A @ B requires A.domain == B.codomain; "
                f"got A.domain={a_dom!r}, B.codomain={b_cod!r}."
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
    def codomain(self) -> Optional["FunctionSpace"]:
        # A @ B: output space is A.codomain.
        return getattr(self.a, "codomain", None)

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
    def codomain(self) -> Optional["FunctionSpace"]:
        return getattr(self.op, "codomain", None)

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

    Has ``apply`` (returns a zero of the same shape and type as ``x``)
    and ``apply_transpose`` (also zero), but **not** ``solve``: the
    zero operator is not invertible.  Forcing it to advertise
    ``solve`` would be the harmful-stub anti-pattern this whole
    module is designed against — composers downstream that ask for
    ``solve`` on :math:`L = A + 0` would silently get a meaningless
    answer.

    Type preservation
    -----------------

    The implementation routes through ``0.0 * x``, which preserves
    the input type via the right multiplication dunder
    (:meth:`__rmul__` / :meth:`__mul__`):

    * Bare ``np.ndarray`` — numpy scalar multiply returns
      ``np.zeros_like(x)`` bit-exactly.
    * Typed flux containers
      (:class:`~orpheus.sn.angular_flux.AngularFlux`,
      :class:`~orpheus.sn.scalar_flux.ScalarFlux`,
      :class:`~orpheus.sn.harmonic_moment_field.HarmonicMomentField`) —
      dunder arithmetic produces a fresh typed instance whose values
      (and ``.boundary`` for AngularFlux) are all zero.

    This is the load-bearing detail for the R-1 Step 4 typed-end-to-
    end algebra: composing ``(L + C - S - F)`` where ``F =
    ZeroOperator`` (within-group inner solve) MUST produce a typed
    AngularFlux at every step; a bare ``np.zeros_like(AngularFlux)``
    would raise.
    """

    capabilities: frozenset[str] = frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})

    def apply(self, x):
        return 0.0 * x

    def apply_transpose(self, x):
        return 0.0 * x


class PermutationOperator(LinearOperatorMixin):
    r"""Index permutation along a configurable axis: :math:`(P x)_i = x_{\pi(i)}`.

    For a permutation :math:`\pi : \{0, \ldots, N-1\} \to \{0, \ldots, N-1\}`
    represented as an integer array ``perm`` of length :math:`N`, the
    apply action gathers entries along ``axis`` according to ``perm``:

    .. math::

        (P\,x)_{i_0 \ldots i_{a-1}\,j\,i_{a+1} \ldots} \;=\;
        x_{i_0 \ldots i_{a-1}\,\pi(j)\,i_{a+1} \ldots}

    The transpose is the inverse permutation :math:`\pi^{-1}`, computed
    once at construction via ``np.argsort(perm)``. When
    ``perm[perm] == np.arange(N)`` the permutation is an involution
    (:math:`\pi \circ \pi = \mathrm{id}`) — detected at construction
    and exposed as :attr:`is_involution` for downstream consumers that
    benefit from knowing self-adjointness in the unweighted inner
    product. In particular, SN specular reflection through
    :func:`~orpheus.sn.quadrature.Quadrature.reflection_index` is an
    involution; periodic shifts and rotational reorderings are not.

    Capability set: ``{CAP_APPLY, CAP_APPLY_TRANSPOSE}``. ``solve`` is
    NOT advertised even though permutations are invertible — the
    standard solve idiom is :meth:`apply_transpose` (since
    :math:`P^{-1} = P^T` for permutations); a user who wants ``solve``
    semantics should compose with the explicit inverse.

    Parameters
    ----------
    perm
        1-D integer array of length :math:`N` whose entries are a
        permutation of :math:`\{0, \ldots, N-1\}`. Validated at
        construction; rejecting duplicates and out-of-range entries
        with :class:`ValueError`.
    axis
        Tensor axis along which the permutation acts. The operator
        broadcasts on every other axis.

    Attributes
    ----------
    perm : np.ndarray
        Forward permutation :math:`\pi`, as 1-D ``intp`` array.
    inverse_perm : np.ndarray
        Inverse permutation :math:`\pi^{-1}`, precomputed via
        :func:`numpy.argsort`.
    axis : int
        The tagged tensor axis.
    n : int
        Length of the permuted axis.
    is_involution : bool
        ``True`` iff :math:`\pi \circ \pi = \mathrm{id}`.
    """

    capabilities: frozenset[str] = frozenset(
        {CAP_APPLY, CAP_APPLY_TRANSPOSE, CAP_SOLVE}
    )

    def __init__(self, perm: np.ndarray, axis: int = 0) -> None:
        perm = np.asarray(perm, dtype=np.intp)
        if perm.ndim != 1:
            raise ValueError(
                f"PermutationOperator perm must be 1-D; got shape {perm.shape}"
            )
        n = perm.size
        # Validate: perm is a true permutation of {0, ..., n-1}.
        if n == 0 or not (
            perm.min() == 0
            and perm.max() == n - 1
            and np.unique(perm).size == n
        ):
            raise ValueError(
                "PermutationOperator perm must be a permutation of "
                f"{{0, ..., {n - 1}}}; got {perm!r}."
            )
        self.perm = perm
        self.inverse_perm = np.argsort(perm)
        self.axis = int(axis)
        self.n = n
        # Involution detection: perm[perm] == arange(n).
        self.is_involution: bool = bool(
            np.array_equal(perm[perm], np.arange(n))
        )

    def apply(self, x: np.ndarray) -> np.ndarray:
        return np.take(x, self.perm, axis=self.axis)

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        return np.take(x, self.inverse_perm, axis=self.axis)

    def solve(self, b: np.ndarray) -> np.ndarray:
        r"""Inverse permutation action: :math:`P^{-1} b`.

        For a permutation matrix, :math:`P^{-1} = P^{T}` — the
        inverse equals the transpose. Both :meth:`solve` and
        :meth:`apply_transpose` therefore route through the same
        ``np.take(b, inverse_perm, axis=axis)`` call. The pair is
        provided separately because they advertise different
        capabilities (``CAP_SOLVE`` vs ``CAP_APPLY_TRANSPOSE``)
        and consumers may filter by either.
        """
        return np.take(b, self.inverse_perm, axis=self.axis)


class IncomingOrdinateMaskTensor(LinearOperatorMixin):
    r"""Sparse inflow-ordinate mask: zeroes selected entries along an axis.

    For SN vacuum BCs at face :math:`f`, the inflow ordinate set is
    :math:`\{n : \mathrm{sign}(\Omega_n \cdot \hat n_f) < 0\}`. The
    canonical vacuum-trace representation (Grand Report v3 §16A.10
    line 3165) is the operator that ZEROES those inflow entries while
    leaving the outflow entries untouched:

    .. math::

        (M\,\psi)_n \;=\;
        \begin{cases}
            0      & n \in \mathcal{I}_{\text{in}} \\
            \psi_n & \text{otherwise}
        \end{cases}

    Distinct from :class:`ZeroOperator` (which zeroes ALL entries):
    the sparse mask preserves the outflow trace, which downstream
    consumers (sensitivity adjoints, post-BC field readers) require.

    Self-adjoint (:math:`M = M^T = M^* = M^{1/2}`). Idempotent
    (:math:`M^2 = M`) — projection onto the outflow subspace. The
    apply action returns a copy; original input is unmodified.

    Capability set: ``{CAP_APPLY, CAP_APPLY_TRANSPOSE}``. ``solve`` is
    NOT advertised — the mask is rank-deficient (it projects), so it
    is non-invertible. Forcing a ``solve`` stub would be the
    harmful-stub anti-pattern this module is designed against.

    Parameters
    ----------
    inflow_indices
        1-D integer array of ordinate indices to zero. May be empty
        (then the operator is identity on the chosen axis).
        Duplicates and out-of-range entries are rejected at
        construction.
    n_ordinates
        Length of the masked axis. Indices must satisfy
        ``0 <= idx < n_ordinates``.
    axis
        Tensor axis along which the mask acts.
    """

    capabilities: frozenset[str] = frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})

    def __init__(
        self,
        inflow_indices: np.ndarray,
        n_ordinates: int,
        axis: int = 0,
    ) -> None:
        inflow_indices = np.asarray(inflow_indices, dtype=np.intp)
        if inflow_indices.ndim != 1:
            raise ValueError(
                f"IncomingOrdinateMaskTensor inflow_indices must be 1-D; "
                f"got shape {inflow_indices.shape}"
            )
        if inflow_indices.size > 0:
            if inflow_indices.min() < 0 or inflow_indices.max() >= n_ordinates:
                raise ValueError(
                    f"IncomingOrdinateMaskTensor inflow_indices out of range "
                    f"[0, {n_ordinates}); got min={int(inflow_indices.min())}, "
                    f"max={int(inflow_indices.max())}"
                )
            if np.unique(inflow_indices).size != inflow_indices.size:
                raise ValueError(
                    "IncomingOrdinateMaskTensor inflow_indices contains duplicates"
                )
        self.inflow_indices = inflow_indices
        self.n_ordinates = int(n_ordinates)
        self.axis = int(axis)

    def apply(self, x: np.ndarray) -> np.ndarray:
        out = np.asarray(x).copy()
        if self.inflow_indices.size == 0:
            return out
        if self.axis == 0:
            out[self.inflow_indices] = 0.0
        else:
            idx: list = [slice(None)] * out.ndim
            idx[self.axis] = self.inflow_indices
            out[tuple(idx)] = 0.0
        return out

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        # Self-adjoint: same code path.
        return self.apply(x)


class PeriodicWrapOperator(LinearOperatorMixin):
    r"""Spatial-pushforward operator for periodic boundaries.

    Represents the angular trace map that connects opposite faces of
    a periodic mesh. Currently the body is angular identity — this
    matches the legacy
    :class:`~orpheus.geometry.boundary.PeriodicBoundaryOperator`
    semantics, where the SN sweep handles the spatial wrap via its
    own face-pair indexing and the BC operator only needs to pass
    the angular trace through unchanged.

    The type is reserved for a future spatial-pushforward extension
    when the periodic map needs to act on a per-cell flux field
    rather than a per-face angular trace (e.g. for coupling periodic
    BCs into a curvilinear Krylov solve where spatial wrap is not
    handled by sweep indexing). See follow-up issue
    "BC: PeriodicWrapOperator spatial-pushforward implementation".

    Capability set: ``{CAP_APPLY, CAP_APPLY_TRANSPOSE}``. The
    identity body is self-adjoint. The output is a fresh copy of
    the input — matching the legacy
    :class:`~orpheus.geometry.boundary.periodic.PeriodicBoundary.apply`
    aliasing-safety contract (``psi_out.copy()``) and the project-
    wide convention that ``op.apply(psi)`` may be mutated freely by
    the caller without affecting ``psi``.

    Wave 7 update: pre-Wave-7 this operator returned ``x`` by
    reference; the Wave 7 delegation of
    :class:`~orpheus.geometry.boundary.PeriodicBoundary` to this
    operator exposed the aliasing-safety mismatch via
    ``tests/geometry/test_boundary.py::test_periodic_bc_returns_input_unchanged``.
    The copy is now performed here so every consumer inherits the
    safe-aliasing contract uniformly.
    """

    capabilities: frozenset[str] = frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})

    def apply(self, x: np.ndarray) -> np.ndarray:
        # Wave 7: return a fresh copy. The legacy
        # ``PeriodicBoundary.apply`` body was ``psi_out.copy()``;
        # delegating through this operator must preserve the
        # caller-mutates-output safe-aliasing contract.
        return np.asarray(x).copy()

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        return np.asarray(x).copy()


class TensorProductOperator(LinearOperatorMixin):
    r"""Per-axis tensor product :math:`A \otimes B \otimes \cdots`.

    Given a tuple of linear operators :math:`A_1, A_2, \ldots, A_k`
    acting on **independent** tensor axes (i.e. each carries an
    ``axis`` attribute and broadcasts on the rest), the tensor product
    operator's action is the sequential per-axis application

    .. math::

        (A_1 \otimes A_2 \otimes \cdots \otimes A_k)\,x
        \;=\; A_k\bigl(\cdots A_2(A_1\,x) \cdots\bigr).

    Because the constituents act on disjoint axes, the order does not
    matter (the operators commute on the joint tensor). The
    :pydata:`capabilities` set is the **intersection** of the
    constituents' capabilities — the tensor product can apply iff
    every factor can apply, can apply_transpose iff every factor can,
    can solve iff every factor can.

    Algebraic laws (verified by tests):

    * **Adjoint distributivity**:
      :math:`(A \otimes B)^* = A^* \otimes B^*`.
    * **Per-axis composition**:
      :math:`(A \otimes B) \circ (C \otimes D) = (A \circ C) \otimes (B \circ D)`
      when ``A``/``C`` share an axis and ``B``/``D`` share an axis.
    * **Inverse on every axis**:
      :math:`(A \otimes B)^{-1} = A^{-1} \otimes B^{-1}` when both
      factors are invertible.

    Parameters
    ----------
    ops : tuple of LinearOperator
        The tensor-product factors. Each MUST advertise an ``axis``
        attribute (or accept an ``axis`` kwarg in :meth:`apply`) and
        broadcast on every other axis. :class:`IdentityOperator`,
        :class:`DiagonalOperator`, and any
        :class:`OperatorProduct`/:class:`OperatorSum` of such operators
        satisfy the contract.

    Notes
    -----

    Relation to numpy: :func:`numpy.kron`, :func:`numpy.tensordot`,
    :func:`numpy.einsum` are array primitives — the *implementation*
    layer. :class:`TensorProductOperator` is the *operator algebra
    type* — it carries axis tags, capability set, and the algebraic
    laws above. Its :meth:`apply` routes through each constituent's
    :meth:`apply`, which is itself typically a single ``np.einsum`` or
    broadcast-multiply. Different abstraction layers, complementary.

    Operators with non-axis-preserving signatures (e.g. an angular
    moment projection that consumes one ordinate axis and produces
    two harmonic-coefficient axes) do not fit this contract — their
    action changes tensor rank. Use them directly via their own
    :meth:`apply`; do not wrap in :class:`TensorProductOperator`.
    """

    def __init__(self, ops: tuple) -> None:
        if len(ops) < 1:
            raise ValueError("TensorProductOperator requires at least one factor")
        # Validate every factor advertises CAP_APPLY.
        for op in ops:
            if not _has(op, CAP_APPLY):
                raise MissingCapability(
                    f"TensorProductOperator factor must advertise "
                    f"{CAP_APPLY!r}; {type(op).__name__} lacks it."
                )
        self.ops: tuple = tuple(ops)
        # Capability intersection: tensor product has cap c iff EVERY
        # factor has cap c.
        recognised = {CAP_APPLY, CAP_SOLVE, CAP_APPLY_TRANSPOSE}
        intersected = {
            cap for cap in recognised if all(_has(op, cap) for op in ops)
        }
        self.capabilities = frozenset(intersected)

    @staticmethod
    def _build(a: "LinearOperator", b: "LinearOperator") -> "TensorProductOperator":
        """Construct a flattened ``A & B`` instance.

        If either operand is itself a :class:`TensorProductOperator`,
        absorb its factors so ``(A & B) & C`` and ``A & (B & C)`` both
        produce ``TensorProductOperator((A, B, C))``.
        """
        a_ops = a.ops if isinstance(a, TensorProductOperator) else (a,)
        b_ops = b.ops if isinstance(b, TensorProductOperator) else (b,)
        return TensorProductOperator(a_ops + b_ops)

    def apply(self, x: np.ndarray) -> np.ndarray:
        out = x
        for op in self.ops:
            out = op.apply(out)
        return out

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        if CAP_APPLY_TRANSPOSE not in self.capabilities:
            raise MissingCapability(
                "TensorProductOperator.apply_transpose requires every "
                "factor to advertise CAP_APPLY_TRANSPOSE."
            )
        out = x
        # Adjoint of tensor product is tensor product of adjoints.
        # Apply transposes of factors (order irrelevant for disjoint axes).
        for op in self.ops:
            out = op.apply_transpose(out)  # type: ignore[attr-defined]
        return out

    def solve(self, b_vec: np.ndarray) -> np.ndarray:
        if CAP_SOLVE not in self.capabilities:
            raise MissingCapability(
                "TensorProductOperator.solve requires every factor to "
                "advertise CAP_SOLVE."
            )
        out = b_vec
        for op in self.ops:
            out = op.solve(out)  # type: ignore[attr-defined]
        return out


class SumOfTensorProductsOperator(LinearOperatorMixin):
    r"""Sum of tensor products :math:`\sum_k A_k \otimes B_k \otimes \cdots`.

    The §15.2 / §15A.2 canonical form for scattering and streaming in
    the operator-algebra view:

    * **Streaming** (§15.1):
      :math:`L = D_x \otimes \Omega_x \otimes I_g + D_y \otimes \Omega_y \otimes I_g`.
    * **Scattering** (§15.2):
      :math:`S = \sum_\ell P_\ell \otimes \Sigma_{s,\ell}` (per-:math:`\ell`
      block-diagonal on moment space).

    Algebraically just :class:`OperatorSum` over
    :class:`TensorProductOperator` summands, but exposed as a named
    type because the structure carries V&V invariants worth checking
    explicitly:

    * Each summand IS a :class:`TensorProductOperator` —
      :meth:`assert_separable`.
    * (Future) common-axis factorisation — when many summands share
      an axis-factor, refactoring saves work.

    Parameters
    ----------
    summands : tuple of TensorProductOperator
        The tensor-product summands. Each MUST be a
        :class:`TensorProductOperator`; mixing in non-separable
        operators makes the type label misleading.

    Notes
    -----

    The implementation backs onto :class:`OperatorSum` —
    :meth:`apply` simply sums each summand's action — so all the
    algebra of :class:`OperatorSum` (composition, scaling, capability
    intersection) is inherited by delegation. The named subclass
    exists for the type signal and the assertion methods.
    """

    def __init__(self, summands: tuple) -> None:
        if len(summands) < 1:
            raise ValueError(
                "SumOfTensorProductsOperator requires at least one summand"
            )
        for s in summands:
            if not isinstance(s, TensorProductOperator):
                raise TypeError(
                    f"SumOfTensorProductsOperator summands must be "
                    f"TensorProductOperator instances; got {type(s).__name__}. "
                    f"Use OperatorSum for general operator addition."
                )
        self.summands: tuple = tuple(summands)
        # Capability intersection across summands (sum can apply iff
        # every summand can apply, etc.).
        recognised = {CAP_APPLY, CAP_APPLY_TRANSPOSE}
        # Solve does NOT propagate through sums (sum of inverses != inverse of sum).
        intersected = {
            cap for cap in recognised if all(_has(s, cap) for s in summands)
        }
        self.capabilities = frozenset(intersected)

    def apply(self, x: np.ndarray) -> np.ndarray:
        out = self.summands[0].apply(x)
        for s in self.summands[1:]:
            out = out + s.apply(x)
        return out

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        if CAP_APPLY_TRANSPOSE not in self.capabilities:
            raise MissingCapability(
                "SumOfTensorProductsOperator.apply_transpose requires "
                "every summand to advertise CAP_APPLY_TRANSPOSE."
            )
        out = self.summands[0].apply_transpose(x)
        for s in self.summands[1:]:
            out = out + s.apply_transpose(x)
        return out

    def assert_separable(self) -> None:
        """Assert every summand is a :class:`TensorProductOperator`.

        Holds by construction (the constructor enforces it), so this
        method is a no-op contract-validator. Useful as documentation
        and as a hook for subclasses or future invariant checks.
        """
        for s in self.summands:
            if not isinstance(s, TensorProductOperator):
                raise AssertionError(
                    f"SumOfTensorProductsOperator summand is not "
                    f"separable: {type(s).__name__}"
                )


class DiagonalOperator(LinearOperatorMixin):
    r"""Diagonal multiplication on a tagged tensor axis.

    For a 1-D weight vector :math:`w \in \mathbb{R}^N` and target axis
    ``axis``, the operator acts on a multi-axis tensor :math:`x` by
    elementwise multiplication along ``axis``:

    .. math::

        (D x)_{\ldots,\,n,\,\ldots} \;=\; w_n \, x_{\ldots,\,n,\,\ldots}

    All other axes broadcast through unchanged. This is the canonical
    "diagonal in some basis" operator — the abstraction the Grand Report
    v3 §9 names :math:`W` (``AngularWeightMatrix``) when the basis is
    the discrete-ordinate set of an angular cubature, and the same
    primitive any method needs for "multiply-by-weights along one axis"
    (MoC track-weight diagonal, CP region-volume weighting, MC
    importance weighting).

    Self-adjoint by construction (real-valued weights), so
    :meth:`apply_transpose` is the same code path as :meth:`apply`.
    Invertible iff all weights are non-zero, in which case
    :pydata:`CAP_SOLVE` is advertised and :meth:`solve` divides by the
    weights along ``axis``.

    Parameters
    ----------
    weights : np.ndarray
        1-D weight vector :math:`w`, shape ``(N,)``.
    axis : int, default 0
        Tensor axis on which to multiply. The operator broadcasts on
        every other axis.

    Notes
    -----

    Construction does NOT eagerly materialise an :math:`N \times N`
    diagonal matrix; the action is implemented as a single
    broadcast-multiply (``self._reshape() * x``) so memory cost is
    :math:`O(N)` regardless of the input tensor's shape.

    Use :meth:`from_measure` when the weights live on a
    :class:`~orpheus.numerics.measure.DiscreteMeasure` — common for the
    angular axis of an SN field, where the operator is built from
    ``quad.weights``.
    """

    def __init__(self, weights: np.ndarray, axis: int = 0) -> None:
        weights_arr = np.asarray(weights, dtype=float)
        if weights_arr.ndim != 1:
            raise ValueError(
                f"DiagonalOperator weights must be 1-D, got shape "
                f"{weights_arr.shape}"
            )
        self.weights = weights_arr
        self.axis = int(axis)

        # Solve is supported iff every weight is non-zero. We check
        # eagerly at construction so the capability is honest.
        invertible = bool(np.all(weights_arr != 0.0))
        caps = {CAP_APPLY, CAP_APPLY_TRANSPOSE}
        if invertible:
            caps.add(CAP_SOLVE)
        self.capabilities = frozenset(caps)

    @classmethod
    def from_measure(
        cls, measure, axis: int = 0,
    ) -> "DiagonalOperator":
        """Construct from the weights of a :class:`DiscreteMeasure`.

        Convenience constructor for the canonical case where the
        diagonal IS the discrete measure's weights — e.g.
        ``DiagonalOperator.from_measure(quad.measure, axis=0)`` is the
        Grand Report v3 §9 ``AngularWeightMatrix``.
        """
        return cls(measure.weights, axis=axis)

    def _reshape(self, ndim: int) -> np.ndarray:
        """Reshape ``self.weights`` to broadcast over an ``ndim``-tensor.

        Returns a view of shape ``(1, ..., 1, N, 1, ..., 1)`` with
        ``N`` at position ``self.axis``.
        """
        shape = [1] * ndim
        shape[self.axis] = -1
        return self.weights.reshape(shape)

    def apply(self, x: np.ndarray) -> np.ndarray:
        x_arr = np.asarray(x)
        if x_arr.shape[self.axis] != self.weights.shape[0]:
            raise ValueError(
                f"DiagonalOperator(axis={self.axis}) expects axis size "
                f"{self.weights.shape[0]}; got {x_arr.shape[self.axis]} "
                f"in input of shape {x_arr.shape}."
            )
        return self._reshape(x_arr.ndim) * x_arr

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        # Real-valued diagonal is self-adjoint.
        return self.apply(x)

    def solve(self, b_vec: np.ndarray) -> np.ndarray:
        if CAP_SOLVE not in self.capabilities:
            raise MissingCapability(
                "DiagonalOperator.solve requires non-zero weights; "
                "this operator has at least one zero weight."
            )
        b_arr = np.asarray(b_vec)
        if b_arr.shape[self.axis] != self.weights.shape[0]:
            raise ValueError(
                f"DiagonalOperator(axis={self.axis}) expects axis size "
                f"{self.weights.shape[0]}; got {b_arr.shape[self.axis]} "
                f"in input of shape {b_arr.shape}."
            )
        return b_arr / self._reshape(b_arr.ndim)


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
