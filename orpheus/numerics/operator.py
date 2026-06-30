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
"""

from __future__ import annotations

from enum import Enum
from typing import (
    TYPE_CHECKING,
    Callable,
    ClassVar,
    Generic,
    Optional,
    Protocol,
    TypeVar,
    cast,
    runtime_checkable,
)

import numpy as np

from orpheus.numerics.vector import Vector

if TYPE_CHECKING:
    from orpheus.numerics.functional import Functional
    from orpheus.numerics.space import FunctionSpace


# ── Two-parameter operator typevars (#65 / P4.5) ──────────────────────
# The honest operator type is ``LinearOperator[Domain, Codomain]``:
# ``apply`` maps an input carrier ``Domain`` to a (possibly distinct)
# output carrier ``Codomain`` — the carrier's ``(Representation, Role)``
# grid cell IS the operator's ``(Domain, Codomain)`` (the double-category
# 1-morphism between cells; see :ref:`operator-algebra`). The names are
# spelled in full (NOT ``Din``/``Cout``) because ``Domain`` already reads
# as "in" and ``Codomain`` as "out" — the abbreviation said nothing.
#
# ONE invariant pair (W-A collapse, #65): :class:`LinearOperator` is now a
# SINGLE base — a ``@runtime_checkable`` Protocol that ALSO carries the
# algebra dunders as default-method bodies — so there is no longer a
# separate variant read-Protocol and an invariant impl-Mixin to reconcile.
# ``Domain``/``Codomain`` are therefore **invariant**: the variance the old
# split needed never reached the leaves (every leaf inherits the one base
# and the static carrier collapses to ``Vector`` at the numerics layer, so
# co/contra-variance bought nothing — and a contravariant TypeVar cannot be
# passed to the invariant composer bases). PEP-696 default
# ``Codomain = Domain`` makes ``LinearOperator[V] ≡ LinearOperator[V, V]``,
# so the endomorphic majority — and every existing single-parameter
# subscript site — keeps one parameter. The native
# ``typing.TypeVar(default=…)`` requires-python ``>=3.13``.
Domain = TypeVar("Domain", bound=Vector)  # operator input carrier
Codomain = TypeVar("Codomain", bound=Vector, default=Domain)  # operator output carrier
Cmid = TypeVar("Cmid", bound=Vector)  # OperatorProduct intermediate carrier
D2 = TypeVar("D2", bound=Vector)  # __matmul__ other-operand domain

__all__ = [
    "LinearOperator",
    "SupportsInverse",
    "SupportsAdjoint",
    "BlockRole",
    "BulkOperator",
    "FullOperator",
    "BoundaryOperator",
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
    "RankOneOperator",
    "outer",
    "TensorProductOperator",
    "SumOfTensorProductsOperator",
    "CAP_APPLY",
    "CAP_SOLVE",
    "CAP_APPLY_TRANSPOSE",
]


# Capability tag literals. Strings (rather than an enum) so user
# operators can advertise method-specific tags without subclassing.
CAP_APPLY: str = "apply"
CAP_SOLVE: str = "solve"
CAP_APPLY_TRANSPOSE: str = "apply_transpose"


# ───────────────────────────────────────────────────────────────────────
# Block-role classification (Issue #208 / Wave O)
# ───────────────────────────────────────────────────────────────────────
#
# On the direct-sum transport state space ``V = V_bulk ⊕ V_boundary`` a
# linear operator is, by the biproduct theorem, a 2×2 block matrix::
#
#     A = [ A_bb  A_bs ]      A_bb : bulk → bulk        A_bs : boundary → bulk
#         [ A_sb  A_ss ]      A_sb : bulk → boundary    A_ss : boundary → boundary
#
# :class:`BlockRole` classifies a leaf by WHICH blocks its action touches —
# the single fact :meth:`OperatorSum.apply` dispatches on (Wave O step O.2)
# and the adjoint composition routes by. The classification is a partition
# (each leaf is exactly one role), and it lives on the INSTANCE (via the
# :attr:`LinearOperator.block_role` attribute), NOT the class, because
# the same generic operator class can play different roles in different
# contexts — e.g. :class:`IdentityOperator` is the bulk identity in one
# composition and a realized vacuum boundary law in another.


class BlockRole(Enum):
    r"""Which bulk/boundary blocks an operator's action touches.

    * :attr:`BULK` — only ``A_bb`` (bulk → bulk). The collision ``C``,
      scattering ``S`` and fission ``F`` operators: they read the bulk
      flux and write a bulk source/sink, with no boundary action.
    * :attr:`FULL` — has off-diagonal coupling (``A_bs`` and/or ``A_sb``).
      The streaming operator ``L``: it reads the inflow trace to seed the
      sweep and writes the outflow trace, coupling bulk ↔ boundary. The
      only irreducibly-full primitive.
    * :attr:`BOUNDARY` — only ``A_ss`` (boundary → boundary). A realized
      boundary law ``B`` (vacuum / reflective / albedo / white / periodic):
      it maps the outflow trace to the inflow trace, with no bulk action.
      The :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` stamps
      this role on its realized outputs (Wave O step O.4a.1-γ). ``B``
      becomes a first-class algebra leaf — a sibling of ``L`` in
      ``(L_full + C − S − F − B)`` — when the boundary conditions are
      extracted from the streaming sweep (Wave O step O.4a.2); until that
      wiring lands ``B`` carries the role but is still consumed inside the
      sweep. (The :class:`BulkOperator` / :class:`FullOperator` markers
      shipped in O.1; the :class:`BoundaryOperator` marker ships in
      O.4a.1-γ. The ``geometry.boundary.BoundaryOperator`` alias — a
      misnamed re-export of
      :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` — was retired
      in O.4a.1-β, freeing the name for this marker.)
    """

    BULK = "bulk"
    FULL = "full"
    BOUNDARY = "boundary"


def _join_block_roles(
    a: Optional["BlockRole"], b: Optional["BlockRole"],
) -> Optional["BlockRole"]:
    r"""The block role of a SUM ``A + B`` — the union of the blocks touched.

    Reading a role as the *set* of blocks its action touches
    (:attr:`BlockRole.BULK` = ``{bulk}``, :attr:`BlockRole.BOUNDARY` =
    ``{boundary}``, :attr:`BlockRole.FULL` = ``{bulk, boundary}``), the sum
    touches the union: ``BULK ⊔ BULK = BULK``, ``BOUNDARY ⊔ BOUNDARY =
    BOUNDARY``, and any mix (or anything with ``FULL``) is ``FULL``. So the
    join is simply *"same role stays, anything different becomes FULL"*. If
    either operand is unclassified (``None`` — a generic operator outside the
    SN bulk/boundary partition) the sum is unclassified too: ``None``
    propagates (a conservative "don't know" rather than a guessed role).

    This is the derivation that lets ``(L + C - S - F - B)`` carry its role by
    construction (Wave O / O.2b 4.5) — retiring the hand-stamped
    ``InvertibleOperator`` FULL tag. ``L`` is ``FULL``, ``C``/``S``/``F`` are
    ``BULK``, ``B`` is ``BOUNDARY`` → every within-group loss sum joins to
    ``FULL``, exactly the irreducibly bulk↔boundary-coupling streaming role.
    """
    if a is None or b is None:
        return None
    return a if a is b else BlockRole.FULL


class _BlockRoleMeta(type):
    r"""Metaclass making ``isinstance(op, BulkOperator)`` read ``op.block_role``.

    The role markers (:class:`BulkOperator`, :class:`FullOperator`,
    :class:`BoundaryOperator`) are never instantiated and carry no
    state. They exist so the block-role classification reads like the
    domain (``isinstance(L, FullOperator)``) while the single source of
    truth stays the :attr:`~LinearOperator.block_role` enum on the
    operator instance. Exclusivity is automatic: an operator carries one
    ``block_role`` and therefore satisfies at most one marker.

    A value-based check is required (not a plain ``@runtime_checkable``
    :class:`Protocol`) because Protocols can only test attribute
    *presence*, never the *value* the partition discriminates on — every
    operator has a ``block_role`` attribute, so a structural Protocol would
    match them all.
    """

    _role: "BlockRole"

    def __instancecheck__(cls, obj: object) -> bool:
        return getattr(obj, "block_role", None) is cls._role


class BulkOperator(metaclass=_BlockRoleMeta):
    r"""``isinstance`` marker for a :attr:`BlockRole.BULK` operator (``A_bb`` only)."""

    _role = BlockRole.BULK


class FullOperator(metaclass=_BlockRoleMeta):
    r"""``isinstance`` marker for a :attr:`BlockRole.FULL` operator (off-diagonal coupling)."""

    _role = BlockRole.FULL


class BoundaryOperator(metaclass=_BlockRoleMeta):
    r"""``isinstance`` marker for a :attr:`BlockRole.BOUNDARY` operator (``A_ss`` only).

    The realized boundary laws produced by
    :meth:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer.realize`
    (vacuum / reflective / white / albedo / periodic) carry
    :attr:`BlockRole.BOUNDARY`; the rank-0 affine ``PrescribedInflow``
    source does NOT — it is the boundary *source* ``q.boundary``, not a
    linear boundary operator ``B``.
    """

    _role = BlockRole.BOUNDARY


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
class LinearOperator(Protocol[Domain, Codomain]):
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

    #: Block-role classification (Issue #208 / Wave O) — see
    #: :class:`BlockRole`. A single enum value, NOT a capability tag: the
    #: role is a *partition* (an operator is exactly one of
    #: bulk/full/boundary), whereas :attr:`capabilities` is a *lattice*
    #: (apply AND solve AND …, non-exclusive). Modelling the partition as
    #: one enum makes the illegal "BULK and FULL at once" state
    #: unrepresentable; a frozenset would not.
    #:
    #: ``None`` = unclassified — the default for the generic algebra
    #: (composition operators derive their role from operands at O.2).
    #: ``None`` satisfies none of the :class:`BulkOperator` /
    #: :class:`FullOperator` markers. Concrete leaves override with a
    #: **plain (unannotated) class attribute** ``block_role = BlockRole.X``
    #: — NOT a ``ClassVar[...]`` annotation, which under ``from __future__
    #: import annotations`` is mis-detected by the ``@dataclass`` machinery
    #: as a field (it became a string and the ClassVar heuristic missed
    #: it). The annotation HERE is a **plain instance attribute** (NOT
    #: ``ClassVar``) precisely because the composers
    #: (:class:`OperatorSum` / :class:`ScaledOperator` /
    #: :class:`_AdjointOperator`) and the
    #: :func:`~orpheus.sn.boundary.realizer._as_boundary` stamp assign
    #: ``self.block_role`` per-instance (the role is DERIVED from operands,
    #: not fixed by the class). A ``ClassVar`` would (correctly) reject
    #: that instance assignment. This base is not a ``@dataclass``, so the
    #: annotation is never field-processed; the leaves' unannotated
    #: class-attr override keeps the class-level read
    #: (``ScatteringOperator.block_role``) working.
    block_role: Optional[BlockRole] = None

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        """The function space this operator consumes, or ``None``.

        Operators that pre-date Issue 9.6 (and any operator that has
        no canonical function-space tagging) return ``None`` — the
        default supplied by this base. When either operand of a
        composition has ``None`` for ``domain`` or ``codomain``, the
        composability check is skipped — preserving the legacy
        duck-typed behaviour for code paths that do not track spaces.
        """
        return None

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        """The function space this operator produces, or ``None``.

        See :attr:`domain` for the ``None`` semantics.
        """
        return None

    # ------------------------------------------------------------------
    # Per-axis structural capability predicates (#226 inverse-as-operator
    # carve) — the typed, instance-accurate successors to the stringly-
    # typed ``capabilities`` frozenset. Each is the RUNTIME advertisement
    # for one operator-returning method (:meth:`inverse` / :meth:`H`); the
    # propagation LAW lives in the composer method bodies, and these
    # predicates compute the matching "does it work?" answer recursively
    # from the operands — NOT a cached string set that can drift. The
    # STATIC contract is carried by the :class:`SupportsInverse` /
    # :class:`SupportsAdjoint` Protocols.
    # ------------------------------------------------------------------

    @property
    def is_invertible(self) -> bool:
        r"""Whether this operator can produce its inverse OPERATOR (:meth:`inverse`).

        The RUNTIME, instance-accurate successor to the ``CAP_SOLVE``
        capability tag. Unlike ``isinstance(op, SupportsInverse)`` — which
        sees only class-level method presence — this property reads the
        operator's actual structure and values, so it correctly reports a
        generic sum as non-invertible and a zero-coefficient
        :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
        as singular (``min|f| = 0``). Composites derive it recursively from
        their operands; the default is ``False`` — an operator is
        invertible only by explicit override.
        """
        return False

    @property
    def is_adjointable(self) -> bool:
        r"""Whether this operator exposes a Hilbert adjoint (:attr:`H` / transpose).

        The RUNTIME successor to the ``CAP_APPLY_TRANSPOSE`` tag. The
        transpose-of-a-sum law :math:`(A+B)^{\mathsf T} = A^{\mathsf T} +
        B^{\mathsf T}` is realised in the composer method bodies; this
        predicate is the matching *advertisement* —
        ``(A+B).is_adjointable == A.is_adjointable and B.is_adjointable`` —
        structurally computed rather than cached in a string set. Default
        ``False``; an operator with a working ``apply_transpose`` overrides.
        """
        return False

    def apply(self, x: Domain, /) -> Codomain:
        r"""Return :math:`L\,x`.

        Mandatory. Every concrete :class:`LinearOperator` must implement
        this (the body here is the Protocol contract stub). The
        :pydata:`capabilities` set MUST include :pydata:`CAP_APPLY`
        whenever this method is functional.

        The two type variables express the operator honestly: ``apply``
        maps an input carrier :data:`Domain` to a (possibly distinct)
        output carrier :data:`Codomain`. The endomorphic majority (``C``,
        the loss solve, the flat ``np.ndarray`` of the scipy
        serialization boundary) is the special case ``Codomain == Domain``,
        spelled with a single parameter via the PEP-696 default
        (``LinearOperator[V] ≡ LinearOperator[V, V]``); the
        source-producing operators ``S``/``F`` are the genuine
        ``Codomain ≠ Domain`` case (flux carrier → source/sink carrier).
        """
        ...

    # ------------------------------------------------------------------
    # Algebra dunders — default-method bodies (W-A collapse, #65)
    #
    # These carry real bodies ON this Protocol, so an explicit subclass
    # ``class Foo(LinearOperator[A, B])`` inherits BOTH the ``apply``
    # contract AND the natural Python algebra (``+``, ``-``, ``*`` scalar,
    # ``@`` composition) with no separate mixin. The dunders delegate to
    # the composer constructors, which enforce the capability-closure laws.
    # ------------------------------------------------------------------

    def __add__(
        self, other: "LinearOperator[Domain, Codomain]",
    ) -> "OperatorSum[Domain, Codomain]":
        return OperatorSum(self, other)

    def __radd__(
        self, other: "LinearOperator[Domain, Codomain]",
    ) -> "OperatorSum[Domain, Codomain]":
        return OperatorSum(other, self)

    def __sub__(
        self, other: "LinearOperator[Domain, Codomain]",
    ) -> "OperatorSum[Domain, Codomain]":
        return OperatorSum(self, ScaledOperator(-1.0, other))

    def __rsub__(
        self, other: "LinearOperator[Domain, Codomain]",
    ) -> "OperatorSum[Domain, Codomain]":
        return OperatorSum(other, ScaledOperator(-1.0, self))

    def __mul__(self, other: float) -> "ScaledOperator[Domain, Codomain]":
        if not isinstance(other, (int, float, np.floating, np.integer)):
            return NotImplemented
        return ScaledOperator(float(other), self)

    def __rmul__(self, other: float) -> "ScaledOperator[Domain, Codomain]":
        return self.__mul__(other)

    def __neg__(self) -> "ScaledOperator[Domain, Codomain]":
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
        return ScaledOperator(-1.0, self)

    def __truediv__(self, scalar: float) -> "ScaledOperator[Domain, Codomain]":
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
        return ScaledOperator(1.0 / float(scalar), self)

    def __matmul__(
        self, other: "LinearOperator[D2, Domain]",
    ) -> "OperatorProduct[D2, Codomain]":
        # ``self`` (Domain → Codomain) ∘ ``other`` (D2 → Domain): the
        # intermediate carrier is ``self``'s domain ``Domain`` =
        # ``other``'s codomain — captured as ``OperatorProduct``'s
        # ``Cmid``, giving the honest ``D2 → Codomain``.
        return OperatorProduct(self, other)

    def __and__(self, other: "LinearOperator[Domain]") -> "TensorProductOperator":
        r"""Return :math:`A \otimes B` — the per-axis tensor-product operator.

        Per Grand Report v3 §6.3 line 721 and §15.1 line 2044. For two
        operators acting on independent tensor axes, ``A & B`` produces
        the operator whose action is "apply A on its axis, apply B on
        its axis" (sequentially; commutative because axes are disjoint).

        If either operand is already a :class:`TensorProductOperator`,
        the result is flattened so ``(A & B) & C`` and ``A & (B & C)``
        produce the same instance ``TensorProductOperator((A, B, C))``.
        """
        return TensorProductOperator._build(self, other)

    def __rand__(self, other: "LinearOperator[Domain]") -> "TensorProductOperator":
        return TensorProductOperator._build(other, self)

    def __call__(self, *args, **kwargs):
        """Alias for :meth:`apply`. Lets user code write ``A(x)``.

        Accepts ``*args, **kwargs`` so any multi-argument ``apply``
        composes ergonomically (``op(x, y)`` reads as math). The
        forwarding originally served the pre-refactor 2-arg BC apply
        ``apply(psi_out, quadrature)``; that signature was retired
        with the BC descriptor cleanup (every realized BC is now a
        1-arg :class:`LinearOperator`), but the generic forwarding is
        retained for future multi-argument operators.
        """
        return self.apply(*args, **kwargs)

    def __pow__(self, n: int) -> "LinearOperator[Domain, Domain]":
        r"""Return :math:`A^n` for non-negative integer ``n``.

        Only an *endomorphic* operator is powerable (``A @ A`` requires
        ``A``'s codomain to equal its domain), so the return is the
        single-carrier ``[Domain, Domain]`` — a flux→source ``S`` has no
        ``S**2``.

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
            return cast("LinearOperator[Domain, Domain]", self)
        # __pow__ is only valid for an endomorphic operator (``A @ A``
        # needs codomain == domain); cast ``self`` to its
        # ``[Domain, Domain]`` view to express the precondition the
        # general ``[Domain, Codomain]`` method cannot carry in its
        # own signature.
        endo = cast("LinearOperator[Domain, Domain]", self)
        result: "LinearOperator[Domain, Domain]" = endo
        for _ in range(n - 1):
            result = result @ endo
        return result

    # ------------------------------------------------------------------
    # Adjoint
    # ------------------------------------------------------------------

    def adjoint(self) -> "LinearOperator[Codomain, Domain]":
        r"""Return the Hilbert adjoint :math:`A^*`.

        The adjoint SWAPS the carriers: for ``A : Domain → Codomain`` the
        adjoint is ``A^* : Codomain → Domain`` (it maps the codomain back
        to the domain), so the return type is the swapped
        ``[Codomain, Domain]``.

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
        return _AdjointOperator(self)

    @property
    def H(self) -> "LinearOperator[Codomain, Domain]":
        """Alias for :meth:`adjoint`. Matches the Grand Report v3 §6.3
        Hilbert-adjoint vocabulary (``A.H`` reads as "A dagger"). Swaps
        the carriers ``[Domain, Codomain] → [Codomain, Domain]`` (see
        :meth:`adjoint`)."""
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


# ───────────────────────────────────────────────────────────────────────
# Static capability Protocols (#226 inverse-as-operator carve)
# ───────────────────────────────────────────────────────────────────────
#
# These are STATIC contracts ONLY (pyright / annotation targets), the
# successors to the ``CAP_SOLVE`` / ``CAP_APPLY_TRANSPOSE`` string tags.
# They are deliberately NOT ``runtime_checkable``: an ``isinstance`` check
# against them reads class-level method presence, which is class-uniform on
# composites (every ``OperatorSum`` defines ``apply_transpose`` even when a
# summand cannot transpose) and blind to value-dependent leaves (a
# zero-coefficient multiplier still has an ``inverse`` method). Forbidding
# ``isinstance`` at the type level steers every runtime query to the
# instance-accurate :attr:`LinearOperator.is_invertible` /
# :attr:`LinearOperator.is_adjointable` property — the single correct
# mechanism (test-architect §1d.3).


class SupportsInverse(Protocol):
    r"""Static contract: an operator that exposes its inverse OPERATOR.

    The pyright / annotation target (``def precondition(L: SupportsInverse)``)
    for the invertibility axis. The RUNTIME, instance-accurate query is the
    :attr:`LinearOperator.is_invertible` property — NOT ``isinstance`` (this
    Protocol is intentionally not ``runtime_checkable``; see the module
    comment above). Together they replace the ``CAP_SOLVE`` string tag:
    Protocol = static contract, property = runtime truth.
    """

    def inverse(self) -> "LinearOperator": ...


class SupportsAdjoint(Protocol):
    r"""Static contract: an operator that exposes a Hilbert adjoint ``.H``.

    The static counterpart of :attr:`LinearOperator.is_adjointable`, matching
    the Grand Report v3 §3.1 ``Supports*`` family. As with
    :class:`SupportsInverse`, the instance-accurate query is the
    :attr:`~LinearOperator.is_adjointable` property, not ``isinstance``.
    Replaces the ``CAP_APPLY_TRANSPOSE`` string tag at the static layer.
    """

    @property
    def H(self) -> "LinearOperator": ...


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


# ---------------------------------------------------------------------------
# Adjoint wrapper
# ---------------------------------------------------------------------------


class _AdjointOperator(LinearOperator[Codomain, Domain], Generic[Domain, Codomain]):
    """Hilbert-adjoint wrapper around a :class:`LinearOperator`.

    Presents the SWAPPED carriers: an inner ``A : Domain → Codomain``
    becomes ``A^* : Codomain → Domain``. The explicit
    ``Generic[Domain, Codomain]`` pins the class's type-parameter order to
    ``[Domain, Codomain]`` (the non-defaulted ``Domain`` first) even though
    the base is the swapped ``LinearOperator[Codomain, Domain]`` — without
    it ``Codomain`` (which carries ``default=Domain``) would land before
    the non-defaulted ``Domain`` in the inferred parameter list, which
    PEP-696 forbids. This is the ONLY composer with ``[Codomain, Domain]``
    order.

    Constructed by :meth:`LinearOperator.adjoint` (and its alias
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

    def __init__(self, inner: "LinearOperator[Domain, Codomain]") -> None:
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
        # The G-adjoint transposes the 2×2 block matrix (A_bs ↔ A_sb^T),
        # which preserves WHICH blocks are touched — so the role is the
        # inner operator's role (Wave O / O.2b 4.5): ``L.H`` is FULL,
        # ``B.H`` is BOUNDARY, ``C.H`` is BULK.
        self.block_role = getattr(inner, "block_role", None)

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        # Adjoint of A: V → W is A.H: W → V — domain swaps with inner.codomain.
        return getattr(self.inner, "codomain", None)

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        return getattr(self.inner, "domain", None)

    def apply(self, y: Codomain) -> Domain:
        if not _has(self.inner, CAP_APPLY_TRANSPOSE):
            raise MissingCapability(
                f"adjoint application requires {CAP_APPLY_TRANSPOSE!r} on "
                f"the inner operator; {type(self.inner).__name__} "
                f"advertises {getattr(self.inner, 'capabilities', frozenset())}."
            )
        # Hilbert-adjoint action:
        #   (A^* y)_V = G_V⁺ ⊙ apply_transpose(G_W ⊙ y)
        # On the adjoint wrapper, ``codomain`` is the inner operator's
        # domain (V) and ``domain`` its codomain (W). The metric application
        # is delegated to the function space's :meth:`~FunctionSpace.apply_metric`
        # / :meth:`~FunctionSpace.apply_inverse_metric` (Wave O / O.2b) so that
        # the SAME wrapper serves BOTH a flat-ndarray metric (e.g. the
        # spherical-harmonic ``(L+1, 2L+1)`` leading-axis metric — bit-identical
        # to the former in-line leading-axis-broadcast multiply, now
        # :meth:`FunctionSpace._broadcast_metric`) AND a
        # composite bulk ⊕ trace metric on a structured ``FullField`` (the
        # direct-sum space applies a per-block metric, with a pseudo-inverse on
        # the singular partial-current trace). The space owns the metric; the
        # adjoint wrapper is metric-representation-agnostic.
        inner_codomain = getattr(self.inner, "codomain", None)
        inner_domain = getattr(self.inner, "domain", None)
        z = inner_codomain.apply_metric(y) if inner_codomain is not None else y
        result = self.inner.apply_transpose(z)  # type: ignore[attr-defined]
        if inner_domain is not None:
            result = inner_domain.apply_inverse_metric(result)
        return result

    def apply_transpose(self, x: Domain) -> Codomain:
        raise NotImplementedError(
            "apply_transpose on an _AdjointOperator wrapper is not "
            "supported in 9.6; if a consumer needs it, take the adjoint "
            "of the original inner operator's transpose directly."
        )


class OperatorSum(LinearOperator[Domain, Codomain]):
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

    def __init__(self, a: LinearOperator[Domain, Codomain], b: LinearOperator[Domain, Codomain]) -> None:
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
        # Block role DERIVED from the operands (Wave O / O.2b 4.5): the sum
        # touches the union of the blocks its summands touch. Replaces the
        # former hand-stamped ``InvertibleOperator`` FULL tag — ``(L+C)`` and
        # the whole ``(L+C-S-F-B)`` loss now carry FULL by construction.
        self.block_role = _join_block_roles(
            getattr(a, "block_role", None), getattr(b, "block_role", None),
        )

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        a_dom = getattr(self.a, "domain", None)
        return a_dom if a_dom is not None else getattr(self.b, "domain", None)

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        a_cod = getattr(self.a, "codomain", None)
        return a_cod if a_cod is not None else getattr(self.b, "codomain", None)

    def apply(self, x: Domain, /) -> Codomain:
        return self.a.apply(x) + self.b.apply(x)

    def apply_transpose(self, x: Codomain, /) -> Domain:
        return self.a.apply_transpose(x) + self.b.apply_transpose(x)  # type: ignore[attr-defined]

    @property
    def is_adjointable(self) -> bool:
        # (A+B)^T = A^T + B^T (the law in :meth:`apply_transpose`) — the
        # sum is adjointable iff BOTH summands are. Structural successor to
        # the ``CAP_APPLY_TRANSPOSE`` closure computed in :meth:`__init__`.
        return self.a.is_adjointable and self.b.is_adjointable

    # is_invertible inherits the base ``False``: there is no general
    # ``(A+B)^{-1}`` from the operand inverses (see the class docstring).
    # The sweep-invertible ``InvertibleOperator`` subclass overrides it.


class OperatorProduct(LinearOperator[Domain, Codomain]):
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

    def __init__(self, a: LinearOperator[Cmid, Codomain], b: LinearOperator[Domain, Cmid]) -> None:
        # ``A @ B``: ``B`` maps the input ``V`` to the intermediate
        # ``Cmid``, ``A`` maps ``Cmid`` to the output ``W`` — so the
        # product is honestly ``V → W`` with ``Cmid`` captured here.
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

    def apply(self, x: Domain, /) -> Codomain:
        return self.a.apply(self.b.apply(x))

    def solve(self, b_vec: Codomain) -> Domain:
        # (AB)^{-1} = B^{-1} A^{-1} — maps the codomain W back to the domain V.
        return self.b.solve(self.a.solve(b_vec))  # type: ignore[attr-defined]

    def apply_transpose(self, x: Codomain, /) -> Domain:
        # (AB)^T = B^T A^T — maps the codomain W back to the domain V.
        return self.b.apply_transpose(self.a.apply_transpose(x))  # type: ignore[attr-defined]

    @property
    def is_invertible(self) -> bool:
        # (AB)^{-1} = B^{-1} A^{-1} (the law in :meth:`solve`) — the product
        # is invertible iff BOTH factors are. Structural successor to the
        # ``CAP_SOLVE`` closure computed in :meth:`__init__`.
        return self.a.is_invertible and self.b.is_invertible

    @property
    def is_adjointable(self) -> bool:
        # (AB)^T = B^T A^T (the law in :meth:`apply_transpose`) — adjointable
        # iff BOTH factors are. Successor to the ``CAP_APPLY_TRANSPOSE`` closure.
        return self.a.is_adjointable and self.b.is_adjointable


class ScaledOperator(LinearOperator[Domain, Codomain]):
    r"""Scalar multiple of a linear operator: :math:`(\alpha L)\,x = \alpha\,(L\,x)`.

    Scalar multiplication is a unitary operation on the capability
    set: every capability of the underlying operator is preserved.
    For ``solve``, :math:`(\alpha L)^{-1} = (1/\alpha)\,L^{-1}` so we
    divide on the way out (and the scalar must be non-zero — caught
    at composition time).
    """

    def __init__(self, scalar: float, op: LinearOperator[Domain, Codomain]) -> None:
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
        # Scaling preserves which blocks the action touches (Wave O / O.2b 4.5).
        self.block_role = getattr(op, "block_role", None)

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        return getattr(self.op, "domain", None)

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        return getattr(self.op, "codomain", None)

    def apply(self, x: Domain, /, *extra, **kwextra) -> Codomain:
        return self.scalar * self.op.apply(x, *extra, **kwextra)

    def solve(self, b_vec: Codomain) -> Domain:
        # (αL)^{-1} = (1/α) L^{-1} — maps the codomain W back to the domain V.
        return (1.0 / self.scalar) * self.op.solve(b_vec)  # type: ignore[attr-defined]

    def apply_transpose(self, x: Codomain, /, *extra, **kwextra) -> Domain:
        return self.scalar * self.op.apply_transpose(x, *extra, **kwextra)  # type: ignore[attr-defined]

    @property
    def is_invertible(self) -> bool:
        # (αL)^{-1} = (1/α) L^{-1} — α ≠ 0 is enforced at construction, so
        # the scaled operator is invertible iff the operand is.
        return self.op.is_invertible

    @property
    def is_adjointable(self) -> bool:
        # (αL)^T = α L^T — scaling preserves adjointability.
        return self.op.is_adjointable


class IdentityOperator(LinearOperator[Domain]):
    r"""The identity operator :math:`I\,x = x`.

    All three primitive capabilities are trivially supported:
    :math:`I^{-1} = I` and :math:`I^T = I`. ``solve`` is the same
    code path as ``apply`` — it is just relabelled.
    """

    capabilities: frozenset[str] = frozenset(
        {CAP_APPLY, CAP_SOLVE, CAP_APPLY_TRANSPOSE}
    )

    def apply(self, x: Domain, /) -> Domain:
        return x

    def solve(self, b_vec: Domain) -> Domain:
        return b_vec

    def apply_transpose(self, x: Domain, /) -> Domain:
        return x

    @property
    def is_invertible(self) -> bool:
        return True  # I^{-1} = I

    @property
    def is_adjointable(self) -> bool:
        return True  # I^T = I


class ZeroOperator(LinearOperator[Domain, Codomain]):
    r"""The zero operator :math:`0\,x = 0`.

    Has ``apply`` (returns the zero of its CODOMAIN) and
    ``apply_transpose`` (also zero), but **not** ``solve``: the zero
    operator is not invertible.  Forcing it to advertise ``solve`` would
    be the harmful-stub anti-pattern this whole module is designed
    against — composers downstream that ask for ``solve`` on
    :math:`L = A + 0` would silently get a meaningless answer.

    A zero operator's output lives in its CODOMAIN
    ---------------------------------------------

    An operator :math:`A : \mathcal D \to \mathcal C` maps the domain
    :math:`\mathcal D` to the codomain :math:`\mathcal C`; its action —
    including the zero action — produces an element of
    :math:`\mathcal C`.  The zero map therefore yields the *zero of the
    codomain*, which is the input's type ONLY when :math:`\mathcal D =
    \mathcal C` (an endomorphism, or a bare-``np.ndarray`` operator).

    * **Default (``codomain_zero=None``) — endomorphism / bare ndarray.**
      Routes through ``0.0 * x``, echoing the INPUT type via the
      right-multiplication dunder.  Correct when domain == codomain:
      bare ``np.ndarray`` → ``np.zeros_like(x)`` bit-exact; a typed
      endomorphism → a fresh same-class zero.
    * **``codomain_zero`` supplied — genuine map between spaces.** When
      the operator changes space — e.g. a fission operator
      :math:`F : \psi \mapsto q_{\rm fis}` maps FLUX to SOURCE, so its
      zero output is a *zero SOURCE*, never a zero flux — pass
      ``codomain_zero``: a callable ``x ↦ 0_\mathcal C`` building the
      codomain's zero from ``x`` (typically reading ``x``'s mesh/shape).
      ``apply`` returns ``codomain_zero(x)``.

    The typed SN within-group inner solve wires the zero fission slot as
    ``ZeroOperator(codomain_zero=…)`` returning an
    :class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`-bulk
    zero composite, so the typed RHS ``S.apply(ψ) + F.apply(ψ) + q_ext``
    and the Krylov matvec ``L.apply − S.apply − F.apply`` stay CLOSED
    ``AngularSourceSink`` sums (the field-role-typing B.5.2 fix — a
    flux-echoing zero would hit the cross-class gate).  Formal operator
    codomain typing is issue #208; ``codomain_zero`` is the pre-#208 hook
    that keeps the zero operator honest about what space it maps into.
    ``apply_transpose`` stays an input-echo: its codomain is the domain,
    and the transpose of the zero slot is not exercised pre-#208.
    """

    capabilities: frozenset[str] = frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})

    def __init__(
        self, codomain_zero: "Callable[[Domain], Codomain] | None" = None,
    ) -> None:
        self._codomain_zero = codomain_zero

    def apply(self, x: Domain, /) -> Codomain:
        if self._codomain_zero is not None:
            return self._codomain_zero(x)
        # Endomorphic default (domain == codomain ⟹ ``W`` is ``V``): the
        # zero of the codomain IS ``0.0 * x``. ``cast`` is the PEP-484
        # bridge for the one genuinely-untypeable spot — this branch is
        # reached ONLY when no ``codomain_zero`` was supplied, i.e. when
        # the operator is endomorphic and ``W == V``.
        return cast("Codomain", 0.0 * x)

    def apply_transpose(self, x: Domain, /) -> Domain:
        # The transpose's codomain is the domain (V); a zero map's
        # transpose is the zero echo. Not exercised for the non-endo
        # (codomain_zero) case pre-#208.
        return 0.0 * x

    @property
    def is_adjointable(self) -> bool:
        return True  # 0^T = 0

    # is_invertible inherits the base ``False`` — the zero map is singular.


class PermutationOperator(LinearOperator):
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

    @property
    def is_invertible(self) -> bool:
        return True  # P^{-1} = P^T — a permutation is always invertible

    @property
    def is_adjointable(self) -> bool:
        return True


class IncomingOrdinateMaskTensor(LinearOperator):
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

    @property
    def is_adjointable(self) -> bool:
        return True  # M = M^T (self-adjoint projection)

    # is_invertible inherits the base ``False`` — the mask is rank-deficient.


class PeriodicWrapOperator(LinearOperator):
    r"""Spatial-pushforward operator for periodic boundaries.

    Represents the angular trace map that connects opposite faces of
    a periodic mesh. Currently the body is angular identity — this
    matches the
    :class:`~orpheus.geometry.boundary.PeriodicBoundary`
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

    @property
    def is_adjointable(self) -> bool:
        return True  # identity body is self-adjoint


class TensorProductOperator(LinearOperator):
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

    @property
    def is_invertible(self) -> bool:
        # (A⊗B)^{-1} = A^{-1}⊗B^{-1} — invertible iff every factor is
        # (the ``CAP_SOLVE`` intersection in :meth:`__init__`).
        return all(op.is_invertible for op in self.ops)

    @property
    def is_adjointable(self) -> bool:
        # (A⊗B)^T = A^T⊗B^T — adjointable iff every factor is.
        return all(op.is_adjointable for op in self.ops)


class SumOfTensorProductsOperator(LinearOperator):
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

    @property
    def is_adjointable(self) -> bool:
        # ∑ A_k⊗B_k transposes summand-wise — adjointable iff every summand
        # is. (Solve does not propagate through sums, so is_invertible
        # inherits the base ``False``.)
        return all(s.is_adjointable for s in self.summands)


class DiagonalOperator(LinearOperator):
    r"""Diagonal (pointwise) multiplication by a coefficient field.

    The operator multiplies a carrier tensor :math:`x` by a coefficient
    array :math:`c` that occupies a **sub-product** of the carrier's
    axes and is **constant** over the complementary ``broadcast_axes``:

    .. math::

        (D x)_{\mathbf{i}} \;=\;
        c_{\,\mathbf{i}\setminus\mathrm{bcast}} \; x_{\mathbf{i}}

    i.e. ``D.apply(x) == np.expand_dims(c, broadcast_axes) * x``. The
    coefficient's rank equals ``x.ndim - len(broadcast_axes)`` and its
    axes map, in order, onto the carrier axes NOT in ``broadcast_axes``.

    This is the canonical "diagonal in some basis" / pointwise-multiply
    operator. Two regimes it must express:

    * **1-D special case** — a 1-D coefficient on ONE carrier axis,
      broadcast over all others. This is the Grand Report v3 §9
      :math:`W` (``AngularWeightMatrix``) and the "multiply-by-weights
      along one axis" primitive (MoC track-weight diagonal, CP
      region-volume weighting, MC importance weighting). Spell it
      ``DiagonalOperator(w, axis=k)``; the action is rank-agnostic
      (broadcasts over however many other axes the carrier has).
    * **The multigroup-collision case** — a coefficient of shape
      ``(ng, *spatial)`` broadcast over the LEADING ordinate axis of a
      ``(N, ng, *spatial)`` angular flux, so that
      ``D.apply(psi) == sigma[None] * psi``. This is the broadcast
      engine the transport-layer ``MultiplicationOperator(σ_t)``
      delegates to. Spell it
      ``DiagonalOperator(sigma, broadcast_axes=(0,))``.

    Self-adjoint by construction (real-valued coefficient), so
    :meth:`apply_transpose` is the same code path as :meth:`apply`.
    Invertible iff every coefficient entry is non-zero, in which case
    :pydata:`CAP_SOLVE` is advertised and :meth:`solve` divides.

    Parameters
    ----------
    coefficient : np.ndarray
        The coefficient field :math:`c`. Its rank determines how many
        carrier axes it occupies (``x.ndim - len(broadcast_axes)``).
    broadcast_axes : tuple of int, optional
        The carrier axes over which the coefficient is constant — the
        positions :func:`numpy.expand_dims` inserts singleton dims.
        When omitted, the **1-D special case** applies: ``coefficient``
        MUST be 1-D and ``axis`` selects the single carrier axis it
        occupies (rank-agnostic broadcast over the rest).
    axis : int, default 0
        Used ONLY in the 1-D special case (``broadcast_axes is None``):
        the single carrier axis the 1-D coefficient occupies. Ignored
        when ``broadcast_axes`` is given.

    Notes
    -----

    Construction does NOT materialise a dense diagonal matrix; the
    action is a single broadcast-multiply
    (``self._broadcast(x.ndim) * x``) so memory cost is
    :math:`O(\mathrm{size}(c))` regardless of the carrier's shape.

    The two construction modes are unified through one broadcast helper
    (:meth:`_broadcast`): the 1-D ``axis`` mode is rank-agnostic (the
    carrier's rank is read at apply-time, so the same operator acts on
    a 1-D, 2-D, or N-D carrier), whereas the explicit ``broadcast_axes``
    mode pins both the coefficient rank and the complement, which a
    multi-axis coefficient requires.

    The ``weights``/``axis`` attributes remain available in the 1-D
    case (``coefficient`` is exposed for both modes) as the back-compat
    alias for ``from_measure`` ergonomics and the existing 1-D call
    sites. Composition (``&`` / :class:`TensorProductOperator`,
    ``@`` / :class:`OperatorProduct`) does NOT read them — every composer
    routes purely through ``apply`` / ``solve``.

    Use :meth:`from_measure` when a 1-D coefficient lives on a
    :class:`~orpheus.numerics.measure.DiscreteMeasure` — common for the
    angular axis of an SN field, where the operator is built from
    ``quad.weights``.
    """

    def __init__(
        self,
        coefficient: np.ndarray,
        broadcast_axes: tuple[int, ...] | None = None,
        *,
        axis: int = 0,
    ) -> None:
        coeff = np.asarray(coefficient, dtype=float)

        if broadcast_axes is None:
            # 1-D special case: a single-axis coefficient broadcasting,
            # rank-agnostically, over every other carrier axis. The
            # carrier rank is unknown at construction, so the broadcast
            # placement is deferred to apply-time (see _broadcast).
            if coeff.ndim != 1:
                raise ValueError(
                    f"DiagonalOperator without broadcast_axes is the 1-D "
                    f"special case; coefficient must be 1-D, got shape "
                    f"{coeff.shape}. For an N-D coefficient pass "
                    f"broadcast_axes=(...)."
                )
            self.broadcast_axes: tuple[int, ...] | None = None
        else:
            # General case: an N-D coefficient pinned onto an explicit
            # complement of carrier axes.
            bcast = tuple(int(a) for a in broadcast_axes)
            if len(set(bcast)) != len(bcast):
                raise ValueError(
                    f"DiagonalOperator broadcast_axes must be distinct, "
                    f"got {bcast}."
                )
            self.broadcast_axes = bcast

        # ``axis`` is consulted ONLY in the 1-D special case; storing it
        # as a plain int (default 0) keeps the attribute well-typed and
        # harmless in broadcast mode.
        self.axis = int(axis)
        self.coefficient = coeff

        # Solve is supported iff every coefficient entry is non-zero. We
        # check eagerly at construction so the capability is honest
        # (Pattern 4: an operator that cannot invert never advertises
        # CAP_SOLVE; the legacy bare-σ collision path had no such gate
        # and produced silent IEEE NaN on σ=0).
        invertible = bool(np.all(coeff != 0.0))
        caps = {CAP_APPLY, CAP_APPLY_TRANSPOSE}
        if invertible:
            caps.add(CAP_SOLVE)
        self.capabilities = frozenset(caps)

    @classmethod
    def from_measure(
        cls, measure, axis: int = 0,
    ) -> "DiagonalOperator":
        """Construct from the weights of a :class:`DiscreteMeasure`.

        Convenience constructor for the canonical 1-D case where the
        diagonal IS the discrete measure's weights — e.g.
        ``DiagonalOperator.from_measure(quad.measure, axis=0)`` is the
        Grand Report v3 §9 ``AngularWeightMatrix``.
        """
        return cls(measure.weights, axis=axis)

    @property
    def weights(self) -> np.ndarray:
        """The 1-D coefficient vector (the historical ``weights`` name).

        Available ONLY in the 1-D special case (``broadcast_axes is
        None``); it is the back-compat alias for ``from_measure`` and
        the existing 1-D call sites. Reading it on a multi-axis-
        coefficient instance is an illegal state and raises (Pattern 4)
        rather than returning an N-D array under a 1-D name.
        """
        if self.broadcast_axes is not None:
            raise AttributeError(
                "DiagonalOperator.weights is the 1-D special case's "
                "coefficient vector; this operator has an N-D coefficient "
                f"of shape {self.coefficient.shape} on broadcast_axes="
                f"{self.broadcast_axes}. Use .coefficient instead."
            )
        return self.coefficient

    def _broadcast(self, ndim: int) -> np.ndarray:
        """Reshape the coefficient to broadcast over an ``ndim`` carrier.

        Single source of truth for both construction modes:

        * 1-D ``axis`` mode — return a view of shape
          ``(1, ..., 1, N, 1, ..., 1)`` with ``N`` at ``self.axis``
          (rank-agnostic: built fresh for the carrier's actual ``ndim``).
        * explicit ``broadcast_axes`` mode — return
          ``np.expand_dims(coefficient, broadcast_axes)``, inserting a
          singleton at each broadcast axis so the coefficient occupies
          the complementary axes in order.
        """
        if self.broadcast_axes is None:
            shape = [1] * ndim
            shape[self.axis] = -1
            return self.coefficient.reshape(shape)
        return np.expand_dims(self.coefficient, self.broadcast_axes)

    def _check_shape(self, x: np.ndarray) -> None:
        """Validate the carrier's rank/axis sizes against the coefficient."""
        if self.broadcast_axes is None:
            if x.shape[self.axis] != self.coefficient.shape[0]:
                raise ValueError(
                    f"DiagonalOperator(axis={self.axis}) expects axis size "
                    f"{self.coefficient.shape[0]}; got {x.shape[self.axis]} "
                    f"in input of shape {x.shape}."
                )
            return
        expected_rank = x.ndim - len(self.broadcast_axes)
        if self.coefficient.ndim != expected_rank:
            raise ValueError(
                f"DiagonalOperator(broadcast_axes={self.broadcast_axes}) "
                f"expects a rank-{expected_rank} coefficient for a "
                f"{x.ndim}-D carrier; got rank-{self.coefficient.ndim} "
                f"coefficient of shape {self.coefficient.shape}."
            )

    def apply(self, x: np.ndarray) -> np.ndarray:
        x_arr = np.asarray(x)
        self._check_shape(x_arr)
        return self._broadcast(x_arr.ndim) * x_arr

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        # Real-valued diagonal is self-adjoint.
        return self.apply(x)

    def solve(self, b_vec: np.ndarray) -> np.ndarray:
        if CAP_SOLVE not in self.capabilities:
            raise MissingCapability(
                "DiagonalOperator.solve requires non-zero coefficient "
                "entries; this operator has at least one zero entry."
            )
        b_arr = np.asarray(b_vec)
        self._check_shape(b_arr)
        return b_arr / self._broadcast(b_arr.ndim)

    @property
    def is_invertible(self) -> bool:
        # Invertible iff every coefficient entry is non-zero (D^{-1} = 1/c) —
        # mirrors the eager ``CAP_SOLVE`` gate in :meth:`__init__`.
        return bool(np.all(self.coefficient != 0.0))

    @property
    def is_adjointable(self) -> bool:
        return True  # real diagonal is self-adjoint


class RankOneOperator(LinearOperator):
    r"""The rank-1 dyad :math:`|v\rangle\langle w|` — a reconstruction column ⊗ a functional row.

    A :class:`RankOneOperator` is the outer product of a **reconstruction**
    vector :math:`v` (the column, an output direction) and a
    :class:`~orpheus.numerics.functional.Functional` :math:`\langle w|` (the
    row, the contraction). Its action on a carrier :math:`x` is

    .. math::

        (\,|v\rangle\langle w|\,)\,x \;=\; v \,\cdot\, \langle w, x\rangle ,

    i.e. ``reconstruction * functional.evaluate(x)``: the functional contracts
    :math:`x` to the inner product :math:`\langle w, x\rangle` (with
    ``keepdims`` on the contracted axis), and the reconstruction broadcasts back
    over that length-1 axis. The functional OWNS the contraction (its weight and
    axis); the operator only broadcasts — so the matvec routes THROUGH
    ``functional.evaluate`` and there is no parallel inline reduction to drift
    from it.

    Build instances with :func:`outer` (the readable verb,
    ``outer(reconstruction, functional)``). A genuine ``M × K`` rank-1 operator
    (``v ∈ ℝ^M``, ``w ∈ ℝ^K``, ``M ≠ K``) is legal — there is no same-shape
    constraint between the column and the row (the old ``left.shape ==
    right.shape`` check was an artifact of the square-only legacy form).

    Native to the multigroup fission emission
    :math:`F = |\chi\rangle\langle\nu\Sigma_f| =
    \texttt{outer}(\chi,\ \mathrm{ReactionRateFunctional}(\nu\Sigma_f))`
    (Grand Report v3 §15): the production-rate co-vector
    :math:`\langle\nu\Sigma_f, \phi\rangle = \sum_{g'}\nu\Sigma_{f,g'}\phi_{g'}`
    is the
    :class:`~orpheus.transport.reaction_rate_functional.ReactionRateFunctional`
    row, and the emission spectrum :math:`\chi` is the reconstruction column.
    Fission is the **rank-1 (single-mode) degenerate** of the multi-mode
    scattering kernel ``R∘Λ∘M`` (a :class:`~orpheus.numerics.frame.FrameBase`
    manages the analogous *stack* of dyads); see
    :mod:`orpheus.transport.operators.integral_kernel_operator`.

    Relation to :class:`TensorProductOperator`
    -------------------------------------------

    A :class:`RankOneOperator` satisfies the TP-factor contract (it acts on the
    functional's contracted axis and broadcasts on the others), so it composes
    as a TP factor when the algebra wants the type-visible separable form:

    .. code-block:: python

        F_kernel = outer(chi, InnerProductFunctional(nu_sigma_f, axis=0)) & IdentityOperator()

    The :class:`IdentityOperator` factor advertises the spatial-axis broadcast;
    the TP fold reduces to :meth:`RankOneOperator.apply` bit-identically
    (``IdentityOperator.apply`` returns ``x``).

    Capability set: ``{CAP_APPLY, CAP_APPLY_TRANSPOSE}`` when the row is an
    :class:`~orpheus.numerics.functional.InnerProductFunctional` (the usual
    case, including its
    :class:`~orpheus.transport.reaction_rate_functional.ReactionRateFunctional`
    specialisation), else ``{CAP_APPLY}`` alone. Rank-1 operators are
    **non-invertible** by construction (``solve`` is NEVER advertised — the
    kernel is the orthogonal complement of the row along the contracted axis),
    but they DO have a **transpose**: :meth:`apply_transpose` is the dual dyad
    :math:`|w\rangle\langle v|` — swap the column :math:`v` with the row's
    weight :math:`w`, contracting :math:`\langle v,\cdot\rangle` over the same
    axis. This is the Euclidean transpose :math:`A^{T}`; the metric-correct
    Hilbert adjoint :math:`A^\dagger = G^{-1}A^{T}G` is the
    :attr:`~LinearOperator.H` wrapper's job. The fission adjoint
    :math:`F^\dagger\psi^* = \nu\Sigma_f\,(\chi\cdot\psi^*)` is exactly this
    dyad-swap (campaign #276). A nonlinear / opaque functional has no dual
    column, so such a dyad advertises ``apply`` only. See
    :ref:`operator-algebra-adjoint`.

    Parameters
    ----------
    reconstruction : Vector | numpy.ndarray
        The column :math:`v` — the output direction the inner product is
        broadcast against. Aligns with the carrier on the complement of the
        functional's contracted axis; its size on that axis is the output
        dimension ``M``.
    functional : Functional
        The row co-vector :math:`\langle w|` — contracts the carrier to
        :math:`\langle w, x\rangle` over its own axis (typically the leading
        group axis for the multigroup reaction rate). Usually an
        :class:`~orpheus.numerics.functional.InnerProductFunctional` (generic)
        or a
        :class:`~orpheus.transport.reaction_rate_functional.ReactionRateFunctional`
        (the domain-typed reaction rate).

    Notes
    -----
    Bit-identity with the legacy ``(right * x).sum(axis, keepdims) * left``
    formulation is preserved because
    :meth:`~orpheus.numerics.functional.InnerProductFunctional.evaluate` performs
    that exact ``(w * x).sum(axis, keepdims=True)`` reduction — the same numpy
    primitive, the same axis, the same order — and the reconstruction broadcast
    is elementwise. IEEE-754 pairwise-reduction order is preserved.
    """

    capabilities: frozenset[str] = frozenset({CAP_APPLY})

    def __init__(
        self,
        reconstruction: "Vector | np.ndarray",
        functional: "Functional",
    ) -> None:
        # The dyad |v⟩⟨w|: ``reconstruction`` is the column (output) vector v;
        # ``functional`` is the row co-vector ⟨w| that OWNS the contraction (its
        # weight and axis). NO same-shape guard — a genuine M×K rank-1 operator
        # (M ≠ K) is legal; the functional validates its own contraction against
        # the carrier at apply time (the old left.shape == right.shape check was
        # an artifact of the square-only legacy form, not a real constraint).
        self.reconstruction = reconstruction
        self.functional = functional
        # Capabilities are per-INSTANCE: the dual dyad |w⟩⟨v| (apply_transpose)
        # exists iff the row ⟨w| is a genuine co-vector whose weight IS the dual
        # column — i.e. an InnerProductFunctional (the ReactionRateFunctional
        # specialisation included). A nonlinear / opaque functional has no dual
        # column, so its dyad advertises CAP_APPLY only. ``solve`` is NEVER
        # advertised: rank-1 is non-invertible by construction.
        from orpheus.numerics.functional import InnerProductFunctional

        caps = {CAP_APPLY}
        if isinstance(functional, InnerProductFunctional):
            caps.add(CAP_APPLY_TRANSPOSE)
        self.capabilities = frozenset(caps)

    def apply(self, x: np.ndarray) -> np.ndarray:
        # |v⟩⟨w| x = v · ⟨w, x⟩. The functional IS the contraction — it returns
        # the inner product ⟨w, x⟩ with ``keepdims`` on the contracted axis, and
        # the reconstruction broadcasts back over that length-1 axis. Routing the
        # matvec THROUGH ``functional.evaluate`` is what makes the row-factor a
        # first-class object (no parallel inline reduction to drift from it).
        recon = np.asarray(getattr(self.reconstruction, "values", self.reconstruction))
        return recon * self.functional.evaluate(x)

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        # The dual dyad: (|v⟩⟨w|)ᵀ = |w⟩⟨v|. The transpose swaps the column and
        # the row — the new column is the old row's weight w (the dual column),
        # the new row is ⟨v| (the old reconstruction as a co-vector on the SAME
        # contracted axis). So Aᵀx = w · ⟨v, x⟩, the Euclidean transpose (the .H
        # wrapper adds the metric). Reuses the InnerProductFunctional contraction
        # primitive — single source of truth with the forward `apply` row-factor.
        from orpheus.numerics.functional import InnerProductFunctional

        if not isinstance(self.functional, InnerProductFunctional):
            raise MissingCapability(
                "RankOneOperator.apply_transpose requires the row functional to "
                "be an InnerProductFunctional (a co-vector with a dual column); "
                f"got {type(self.functional).__name__} — a nonlinear functional "
                "has no dual column."
            )
        column = np.asarray(
            getattr(self.functional.weight, "values", self.functional.weight)
        )
        dual_row = InnerProductFunctional(
            np.asarray(getattr(self.reconstruction, "values", self.reconstruction)),
            axis=self.functional.axis,
        )
        return column * dual_row.evaluate(x)


def outer(
    reconstruction: "Vector | np.ndarray",
    functional: "Functional",
) -> RankOneOperator:
    r"""Build the rank-1 dyad :math:`|v\rangle\langle w|` from a column and a co-vector.

    The universal constructor for a rank-1 :class:`LinearOperator`: a
    :class:`~orpheus.numerics.vector.Vector` (or ``ndarray``) ``reconstruction``
    :math:`v` (the column, the output direction) tensored with a
    :class:`~orpheus.numerics.functional.Functional` ``functional`` :math:`\langle w|`
    (the row, the contraction). The action is

    .. math::

        (\,|v\rangle\langle w|\,)\,x \;=\; v \,\cdot\, \langle w, x\rangle ,

    i.e. ``reconstruction * functional.evaluate(x)``. Every separable rank-1
    operator in the algebra is one of these; the multi-mode generalisation is a
    sum of dyads managed by a :class:`~orpheus.numerics.frame.FrameBase` (the
    spectral / block-term decomposition). The canonical transport instance is
    the fission emission kernel
    :math:`F = \texttt{outer}(\chi,\ \mathrm{ReactionRateFunctional}(\nu\Sigma_f))`
    (see :class:`~orpheus.transport.operators.fission.FissionOperator`).
    """
    return RankOneOperator(reconstruction, functional)

