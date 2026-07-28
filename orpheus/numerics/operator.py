r"""Linear-operator algebra for matrix-free transport solvers.

The neutron transport eigenvalue problem and its fixed-source cousin
both reduce to compositions of a small set of linear operators acting
on a discrete flux distribution :math:`\psi`:

.. math::

    \Bigl(A - \sum_i g_i\Bigr)\,\psi \;=\; q
    \qquad \text{(fixed source)}

.. math::

    \Bigl(A - \sum_i g_i\Bigr)\,\psi \;=\; \tfrac{1}{k}\,F\,\psi
    \qquad \text{(eigenvalue)}

where :math:`A` is the INVERTIBLE resolvent operand and the
:math:`g_i` are the lagged coupling gains — for SN the binding is
:math:`A = L + C`, streaming plus collision, with gains
:math:`(S,\ B)`, the honest within-group operator :math:`L+C-S-B`;
the letter matters: project-wide, ``L`` names the STREAMING LEAF
(alone not invertible) and the invertible left-hand-side operand is
``A`` — :math:`S` is the scattering source operator, :math:`F` is the
fission source operator (never a gain in the eigenvalue posing: the
outer loop scales it by :math:`1/k`), and :math:`q` is an external
source (Trefethen & Bau 1997, §3.2 frame the matrix-free Krylov
view). For an SN sweep, an MoC
ray-tracer, a CP collision-probability matrix, or a diffusion BiCGSTAB
solve, the *outer* algebra is identical even though the *implementation*
of each operator differs by orders of magnitude in cost and structure.

This module installs the **algebra** as runtime-checkable Protocols.
Any object providing ``apply(x) -> Lx`` participates. Each further
ability is a per-axis THREE-LAYER surface (#226, Design C):

* a **predicate** (:attr:`~LinearOperator.is_invertible` /
  :attr:`~LinearOperator.is_adjointable`) — the runtime,
  instance-accurate truth, reading structure AND values (a
  zero-coefficient multiplier reports ``False``; a sum reports its
  leading term), recursive on composites;
* an **operator-returning method** (``inverse()`` / :attr:`~LinearOperator.H`)
  — the canonical act; ``.H`` lives on the base (one generic wrapper
  realization exists) and refuses EAGERLY (:class:`MissingAdjoint`),
  while ``inverse()`` lives per-class: a structurally-non-invertible
  type simply does not declare it (misuse is a *static* error), and a
  value-dependent type declares it and raises :class:`NotInvertible`;
* a **realization verb** (``solve`` / ``apply_transpose``) — present
  exactly where a native realization exists (the wrap-delegate family
  delegates through ``solve``; the composer transpose laws recurse
  through ``apply_transpose``), never as an exists-but-raises stub.

The checked bridges :func:`invertible` / :func:`adjointable` (PEP-647
``TypeGuard``) convert the runtime predicate into the static permission
at guarded call sites — you cannot obtain the permission without
executing the check. Composition mismatches still fail at COMPOSITION
time, never mid-iteration: the composers guard ``apply`` eagerly
(``TypeError``), ``.H`` gates at construction, and the value-dependent
``inverse()`` guards raise before any inverse object exists — so a
downstream :class:`scipy.sparse.linalg.LinearOperator` consumer never
silently hits a broken stub. Many transport operators have no
efficient inverse action — the scattering source S is never inverted
directly; the fission source F is rank-deficient — and their honest
surface is METHOD ABSENCE, not an advertising flag. See
:ref:`operator-algebra` for the full design rationale.
"""

from __future__ import annotations

from abc import ABCMeta, abstractmethod
from enum import Enum
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Final,
    Generic,
    Optional,
    Protocol,
    TypeGuard,
    TypeVar,
    cast,
    runtime_checkable,
)

import numpy as np

from orpheus.numerics.vector import Vector

if TYPE_CHECKING:
    from orpheus.numerics.assembled_operator import SparseAssembledOperator
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
# ONE invariant pair (#65): :class:`LinearOperator` is now a
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

# Composition-leg type parameters (C4 F1 — the composition wrappers are
# generic over their LEG types, so a named composition subclass carries its
# legs' identities at the type level and its accessors need no casts:
# ``StreamingCollisionOperator = OperatorSum[FF, FF, StreamingOperator,
# MultiplicationOperator]`` reads ``self.a`` as a ``StreamingOperator``).
# COVARIANT because a pinned composition must upcast to the defaulted
# spelling (``StreamingOperator.__add__`` returns ``StreamingCollisionOperator``
# where the base dunder contract says ``OperatorSum[Domain, Codomain]``) —
# which is also why the legs are read-only properties over ``Final``
# storage: covariance is sound only without a setter. PEP-696 defaults
# keep every existing ``OperatorSum[D, C]`` / bare spelling valid.
SummandA = TypeVar(
    "SummandA", bound="LinearOperator", covariant=True,
    default="LinearOperator[Domain, Codomain]",
)
SummandB = TypeVar(
    "SummandB", bound="LinearOperator", covariant=True,
    default="LinearOperator[Domain, Codomain]",
)
FactorA = TypeVar(  # A of ``A @ B`` — maps the intermediate to the output
    "FactorA", bound="LinearOperator", covariant=True,
    default="LinearOperator[Any, Codomain]",
)
FactorB = TypeVar(  # B of ``A @ B`` — maps the input to the intermediate
    "FactorB", bound="LinearOperator", covariant=True,
    default="LinearOperator[Domain, Any]",
)
ScaledOperand = TypeVar(
    "ScaledOperand", bound="LinearOperator", covariant=True,
    default="LinearOperator[Domain, Codomain]",
)

__all__ = [
    "LinearOperator",
    "SupportsInverse",
    "SupportsAdjoint",
    "SupportsAssembly",
    "BlockRole",
    "BulkOperator",
    "FullOperator",
    "BoundaryOperator",
    "NotInvertible",
    "MissingAdjoint",
    "MissingAssembly",
    "invertible",
    "adjointable",
    "assemblable",
    "IncompatibleOperatorComposition",
    "MatrixTooLarge",
    "InverseWrapMixin",
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
]


# ───────────────────────────────────────────────────────────────────────
# Block-role classification (Issue #208)
# ───────────────────────────────────────────────────────────────────────
#
# On the direct-sum transport state space ``V = V_bulk ⊕ V_boundary`` a
# linear operator is, by the biproduct theorem, a 2×2 block matrix::
#
#     A = [ A_bb  A_bs ]      A_bb : bulk → bulk        A_bs : boundary → bulk
#         [ A_sb  A_ss ]      A_sb : bulk → boundary    A_ss : boundary → boundary
#
# :class:`BlockRole` classifies a leaf by WHICH blocks its action touches —
# the single fact :meth:`OperatorSum.apply` dispatches on, and the adjoint
# composition routes by. The classification is a partition
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
      this role on its realized outputs. ``B`` becomes a first-class
      algebra leaf — a sibling of ``L`` in ``(L_full + C − S − F − B)`` —
      when the boundary conditions are extracted from the streaming sweep;
      until that wiring lands ``B`` carries the role but is still consumed
      inside the sweep.
    """

    BULK = "bulk"
    FULL = "full"
    BOUNDARY = "boundary"


class SystemRole(Enum):
    r"""Which of the two coupled systems an operator's action maps between.

    The curvilinear-S\ :sub:`N` within-group system is a 2×2 coupled block
    operator over two systems (see
    ``docs/theory/foundations/coupled_block_operator.rst §coupled-block-operator``):

    .. math::

        \begin{bmatrix} A_{AA} & A_{AB} \\ A_{BA} & A_{BB} \end{bmatrix}
        \begin{bmatrix} \text{transport} \\ \text{ray} \end{bmatrix}

    * **System A** — the transport bulk ⊕ trace (the angular-flux
      :class:`~orpheus.transport.full_field.FullField`: a bulk field ⊕ its
      spatial boundary trace), governed by ``A_AA = L + C − S − B``.
    * **System B** — the ψ½ radial-characteristic ray (the starting-direction
      cells at each radial cell), governed by ``A_BB`` (the radial
      straight-characteristic march).

    This role is the COARSE two-system partition — orthogonal to
    :class:`BlockRole`, which refines the bulk ↔ boundary structure *within*
    System A. An operator carries at most one:

    * :attr:`A` — acts within System A only.
    * :attr:`B` — acts within System B only: the self-block ``A_BB``
      (:class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator`)
      and the ray boundary ``B_b``
      (:class:`~orpheus.sn.operators.boundary.RadialCharacteristicBoundaryOperator`).
    * :attr:`COUPLED` — maps BETWEEN the systems (an off-diagonal block, or the
      assembled 2×2): the ray→bulk seed ``A_AB``
      (:class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicSeeding`),
      the bulk→ray fold ``A_BA``, and the assembled ``CoupledOperator``.

    Reading each role as the SET of systems its action touches (``A = {A}``,
    ``B = {B}``, ``COUPLED = {A, B}``), a sum touches the union — see
    :func:`_join_system_roles`. Operators outside the two-system decomposition
    — every model-generic family (diffusion / CP / MoC) AND the model-generic
    reaction leaves ``C`` / ``S`` / ``F`` that a curvilinear-S\ :sub:`N` context
    COMPOSES into System A but that carry no intrinsic two-system membership —
    leave :attr:`~LinearOperator.system_role` at its ``None`` default: the
    conservative "not part of the ψ½ augmentation" reading, exactly as an
    unclassified :attr:`~LinearOperator.block_role` is ``None``.
    """

    A = "system_a"
    B = "system_b"
    COUPLED = "coupled"


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

    This is what lets ``(L + C - S - F - B)`` carry its role BY
    CONSTRUCTION (no hand-stamped tag): ``L`` is ``FULL``, ``C``/``S``/``F``
    are ``BULK``, ``B`` is ``BOUNDARY`` → every within-group loss sum joins
    to ``FULL``, exactly the irreducibly bulk↔boundary-coupling streaming
    role.

    Twin of :func:`_join_system_roles` (the two-system analogue — ``COUPLED``
    there plays the top-of-lattice role ``FULL`` plays here): both are the SAME
    union-lattice join (two atoms + a top + ``None`` propagation) on ORTHOGONAL
    role axes. Kept as a deliberate twin while only two axes exist — a generic
    ``RoleAxis`` join would need ``setattr``-driven dispatch that regresses the
    #226 pyright ratchet. **Collapse trigger:** a THIRD parallel role axis (a
    DSA / multiphysics role) makes the shared abstraction pay — unify then.
    """
    if a is None or b is None:
        return None
    return a if a is b else BlockRole.FULL


def _join_system_roles(
    a: Optional["SystemRole"], b: Optional["SystemRole"],
) -> Optional["SystemRole"]:
    r"""The system role of a SUM ``A + B`` — the union of the systems touched.

    Reading a role as the *set* of systems its action touches
    (:attr:`SystemRole.A` = ``{A}``, :attr:`SystemRole.B` = ``{B}``,
    :attr:`SystemRole.COUPLED` = ``{A, B}``), the sum touches the union:
    ``A ⊔ A = A``, ``B ⊔ B = B``, and any mix (or anything with ``COUPLED``) is
    ``COUPLED``. So the join is *"same role stays, anything different becomes
    COUPLED"* — the two-system analogue of :func:`_join_block_roles` with
    ``COUPLED`` playing the top-of-lattice role that ``FULL`` plays there. If
    either operand is unclassified (``None`` — an operator outside the
    two-system decomposition) the sum is unclassified too: ``None`` propagates
    (a conservative "don't know" rather than a guessed role).
    """
    if a is None or b is None:
        return None
    return a if a is b else SystemRole.COUPLED


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

    The realized boundary laws produced by the functional method
    realizers —
    :meth:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer.realize`
    (vacuum / reflective / white / albedo / periodic) and
    :meth:`~orpheus.diffusion.boundary_realizer.DiffusionBoundaryRealizer.realize`
    (the albedo family incl. zero-flux, #290) — carry
    :attr:`BlockRole.BOUNDARY` via
    :func:`~orpheus.geometry.boundary.stamp_boundary_role`; the rank-0
    affine ``PrescribedInflow`` source does NOT — it is the boundary
    *source* ``q.boundary``, not a linear boundary operator ``B``.
    """

    _role = BlockRole.BOUNDARY


class NotInvertible(TypeError):
    r"""Asked for the inverse of an operator that cannot produce one.

    The INVERSE-axis refusal: raised **eagerly** by
    :meth:`inverse` overrides (and the inverse-family constructors) when
    the operator's :attr:`~LinearOperator.is_invertible` is ``False`` —
    the VALUE-dependent arm of the two-kinds split. A zero-coefficient
    :class:`DiagonalOperator`, a sum whose leading term is not
    invertible, a product with a singular factor: the TYPE supports
    inversion, this INSTANCE refuses, at construction of the inverse and
    never mid-iteration. (The STRUCTURAL arm — :class:`ZeroOperator`,
    masks, source dyads, for which no inverse exists mathematically —
    does not declare :meth:`inverse` at all, so misuse there is a
    *static* error, not this exception.) ``TypeError`` parentage carries
    the retired ``MissingCapability``'s public contract forward — no
    ``except`` clause written against the old gate changes meaning.
    """


class MissingAdjoint(TypeError):
    r"""Asked for the Hilbert adjoint of an operator that has none.

    The ADJOINT-axis refusal: raised **eagerly**
    by :meth:`LinearOperator.adjoint` / :attr:`LinearOperator.H` when
    :attr:`~LinearOperator.is_adjointable` is ``False`` — at wrapper
    CONSTRUCTION, never lazily at the first ``.apply`` (the pre-carve
    behaviour). Also the refusal of the raw-transpose realization verb
    (``apply_transpose``) on composites whose operands cannot all
    transpose. ``TypeError`` parentage mirrors :class:`NotInvertible`.
    """


class MissingAssembly(TypeError):
    r"""Asked for the sparse assembly of an operator that has none.

    The ASSEMBLY-axis refusal — the third sibling of
    :class:`NotInvertible` (inverse axis) and :class:`MissingAdjoint`
    (adjoint axis): raised **eagerly** by the composer ``assemble()``
    bodies when an operand is not :attr:`~LinearOperator.is_assemblable`
    — a structural emission exists only where a leaf declared one, and a
    composite can recurse (Sum → ``+``, Product → ``@``, Scaled →
    scalar ``*``) only when every leg emits. Operators without a
    stencil realization simply do not declare :meth:`assemble` (misuse
    is a *static* error); the probing
    :meth:`~LinearOperator.as_matrix` remains their total (dense)
    Mat-functor. ``TypeError`` parentage mirrors the sibling refusals.
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


class MatrixTooLarge(RuntimeError):
    r"""A :meth:`LinearOperator.as_matrix` materialization exceeds its size gate.

    A **resource effect on a TOTAL functor**: every
    linear operator on a finite-dimensional space *has* a matrix — the
    functor ``Op → Mat`` is total — but materializing it commits
    :math:`O(n^2)` memory and :math:`n` applications, which this
    environment may refuse. That is why this is a ``RuntimeError``
    (a refused resource commitment), NOT a ``TypeError``/``ValueError``
    (the request is neither ill-typed nor ill-valued), and why there is
    deliberately NO ``is_materializable`` predicate alongside
    :attr:`LinearOperator.is_invertible` / ``is_adjointable``: those are
    *structural restriction* predicates (they read the operator's
    structure and values), whereas the size gate is a pure resource
    precheck that reads nothing but a dimension. Callers that want the
    fallback pattern write ``try: A.as_matrix() except MatrixTooLarge:
    <iterative path>`` — or raise the per-call ``max_dimension``.
    """


def _resolve_basis_shape(
    op: "LinearOperator",
    basis_shape: "tuple[int, ...] | None",
) -> tuple[int, ...]:
    r"""Resolve the basis shape a materialization iterates over.

    The SINGLE SOURCE for the resolution rule shared by
    :meth:`LinearOperator.as_matrix` and the eager
    :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
    constructor (which must know the resolved shape to reshape solutions
    back into carriers): an explicit ``basis_shape`` wins; otherwise the
    operator's own :attr:`~LinearOperator.domain` supplies its ``shape``;
    an operator with neither is un-materializable *as posed* and the
    caller is told both remedies.
    """
    if basis_shape is not None:
        return tuple(int(d) for d in basis_shape)
    domain = op.domain
    if domain is None:
        raise ValueError(
            f"as_matrix on {type(op).__name__}: the operator carries no "
            f"domain FunctionSpace, so the basis shape cannot be derived. "
            f"Either construct the operator with a space, or pass an "
            f"explicit basis_shape= (the element shape apply consumes, "
            f"e.g. (ng, 1) for a meshless single-cell group operator)."
        )
    return tuple(domain.shape)


@runtime_checkable
class LinearOperator(Protocol[Domain, Codomain]):
    r"""Contract for a matrix-free linear operator on a flux vector.

    Any object exposing :meth:`apply` participates. The further
    abilities are per-axis structural surfaces (the module docstring's
    three layers): the recursive predicates
    :attr:`is_invertible`/:attr:`is_adjointable` are the runtime truth,
    ``inverse()``/:attr:`H` the operator-returning acts, and
    ``solve``/``apply_transpose`` the per-class realization verbs —
    declared exactly where a native realization exists, never as
    stubs. There is no capability registry to keep in sync: the single
    source of truth for what an operator can do is the operator's own
    structure and values, read through the predicates.

    Composition operators (:class:`OperatorSum`, :class:`OperatorProduct`,
    :class:`ScaledOperator`) are wired through ``__add__``, ``__sub__``,
    ``__mul__`` (scalar), and ``__matmul__`` (operator product) so the
    typical algebra of the Boltzmann transport equation,
    :math:`(L + C - S - B)`, can be built with the natural Python syntax.
    The composites derive their predicates recursively per the closure
    laws documented in :ref:`operator-algebra`.

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

    #: Block-role classification (Issue #208) — see
    #: :class:`BlockRole`. A single enum value: the role is a
    #: *partition* (an operator is exactly one of bulk/full/boundary),
    #: so one enum makes the illegal "BULK and FULL at once" state
    #: unrepresentable; a set would not.
    #:
    #: ``None`` = unclassified — the default for the generic algebra
    #: (composition operators derive their role from operands).
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
    #: :func:`~orpheus.geometry.boundary.stamp_boundary_role` stamp assign
    #: ``self.block_role`` per-instance (the role is DERIVED from operands,
    #: not fixed by the class). A ``ClassVar`` would (correctly) reject
    #: that instance assignment. This base is not a ``@dataclass``, so the
    #: annotation is never field-processed; the leaves' unannotated
    #: class-attr override keeps the class-level read
    #: (``ScatteringOperator.block_role``) working.
    block_role: Optional[BlockRole] = None

    #: The COARSE two-system membership (:class:`SystemRole` — System A / System
    #: B / COUPLED) of a curvilinear-S_N augmented operator, orthogonal to
    #: :attr:`block_role` (which refines the bulk↔boundary structure *within*
    #: System A). ``None`` for every operator outside the ψ½ two-system
    #: decomposition. Derived through composition exactly as :attr:`block_role`
    #: is — the passthrough (:class:`_AdjointOperator`, :class:`ScaledOperator`)
    #: and the :func:`_join_system_roles` union (:class:`OperatorSum`).
    system_role: Optional[SystemRole] = None

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
    # Per-axis structural predicates (#226 inverse-as-operator carve).
    # Each is the RUNTIME advertisement for one operator-returning
    # method (:meth:`inverse` / :meth:`H`); the propagation LAW lives in
    # the composer method bodies, and these predicates compute the
    # matching "does it work?" answer recursively from the operands —
    # never a cached registry that can drift. The static bridges are
    # :func:`invertible` / :func:`adjointable` (narrowing to
    # :class:`SupportsInverse` / :class:`SupportsAdjoint`).
    # ------------------------------------------------------------------

    @property
    def is_invertible(self) -> bool:
        r"""Whether this operator can produce its inverse OPERATOR (:meth:`inverse`).

        The RUNTIME, instance-accurate predicate. Unlike
        ``isinstance(op, SupportsInverse)`` — which
        sees only class-level method presence — this property reads the
        operator's actual structure and values, so it correctly reports a
        sum with a non-invertible LEADING term as non-invertible and a
        zero-coefficient
        :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
        as singular (``min|f| = 0``). Composites derive it recursively from
        their operands; the default is ``False`` — an operator is
        invertible only by explicit override.
        """
        return False

    @property
    def is_adjointable(self) -> bool:
        r"""Whether this operator exposes a Hilbert adjoint (:attr:`H` / transpose).

        The RUNTIME predicate for the adjoint axis. The
        transpose-of-a-sum law :math:`(A+B)^{\mathsf T} = A^{\mathsf T} +
        B^{\mathsf T}` is realised in the composer method bodies; this
        predicate is the matching *advertisement* —
        ``(A+B).is_adjointable == A.is_adjointable and B.is_adjointable`` —
        structurally computed rather than cached in a string set. Default
        ``False``; an operator with a working ``apply_transpose`` overrides.
        """
        return False

    @property
    def is_assemblable(self) -> bool:
        r"""Whether this operator can emit its sparse assembly (:meth:`assemble`).

        The ASSEMBLY axis — the third structural surface beside
        :attr:`is_invertible` / :attr:`is_adjointable`: ``True`` iff a
        structural ``(row, col, value)`` emission of this operator
        exists (the stencil-assembly third consumption mode of the
        per-cell closure algebra; see
        :class:`~orpheus.numerics.assembled_operator.SparseAssembledOperator`).
        Composites derive it recursively (a sum/product assembles iff
        both legs do — the homomorphism laws in their ``assemble()``
        bodies); the default is ``False`` — an operator is assemblable
        only by explicit override, and the probing :meth:`as_matrix`
        remains the total dense Mat-functor for everything else. The
        static bridge is :func:`assemblable` (narrowing to
        :class:`SupportsAssembly`).
        """
        return False

    def apply(self, x: Domain, /) -> Codomain:
        r"""Return :math:`L\,x`.

        Mandatory. Every concrete :class:`LinearOperator` must implement
        this (the body here is the Protocol contract stub); the
        composers guard its presence eagerly at composition time.

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
    # Algebra dunders — default-method bodies (#65)
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
        convention is that ``-A`` should also work (adjoint-flux sign
        flips, residual corrections ``-L @ delta``, Jacobi splitting).
        """
        return ScaledOperator(-1.0, self)

    def __truediv__(self, scalar: float) -> "ScaledOperator[Domain, Codomain]":
        r"""Scalar division: ``A / α`` is ``(1/α) * A``.

        Reads more naturally than the reciprocal-multiply form
        ``(1.0 / α) * A`` — eigenvalue/Krylov normalisation
        (``F / k_eff``), homogenisation averages.

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

        For two operators acting on independent tensor axes, ``A & B``
        produces the operator whose action is "apply A on its axis, apply
        B on its axis" (sequentially; commutative because axes are
        disjoint). See
        ``docs/theory/foundations/operator_tensor_network.rst §tensor-network-decomposition``.

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
        composes ergonomically (``op(x, y)`` reads as math); the generic
        forwarding is retained for future multi-argument operators.
        """
        return self.apply(*args, **kwargs)

    def __pow__(
        self: "LinearOperator[Domain, Domain]", n: int,
    ) -> "LinearOperator[Domain, Domain]":
        r"""Return :math:`A^n` for non-negative integer ``n``.

        Only an *endomorphic* operator is powerable (``A @ A`` requires
        ``A``'s codomain to equal its domain) — the precondition lives in
        the ``self`` annotation, so ``S**2`` on a flux→source ``S`` is a
        call-site type error, not a runtime surprise.

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
            return self
        result: "LinearOperator[Domain, Domain]" = self
        for _ in range(n - 1):
            result = result @ self
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
        action) and swaps :attr:`domain` ↔ :attr:`codomain`.

        Raises
        ------
        MissingAdjoint
            **Eagerly, here at construction** when
            this operator is not :attr:`is_adjointable` — a wrapper that
            could only fail at its first ``.apply`` is the broken-stub
            anti-pattern this module refuses. The :func:`adjointable`
            guard doubles as the static bridge: the wrapper's
            constructor consumes the narrowed :class:`SupportsAdjoint`.
        """
        if not adjointable(self):
            raise MissingAdjoint(
                f"{type(self).__name__} is not adjointable — .H/.adjoint() "
                f"requires is_adjointable=True (a working apply_transpose "
                f"on every constituent). The Hilbert adjoint of this "
                f"operator does not exist as posed."
            )
        return _AdjointOperator(self)

    @property
    def H(self) -> "LinearOperator[Codomain, Domain]":
        """Alias for :meth:`adjoint` — the Hilbert-adjoint vocabulary
        (``A.H`` reads as "A dagger"). Swaps the carriers
        ``[Domain, Codomain] → [Codomain, Domain]`` (see
        :meth:`adjoint`)."""
        return self.adjoint()

    # ------------------------------------------------------------------
    # Materialization — the functor OUT of the operator category
    # ------------------------------------------------------------------

    def as_matrix(
        self,
        *,
        basis_shape: tuple[int, ...] | None = None,
        max_dimension: int = 4096,
    ) -> np.ndarray:
        r"""Materialize the explicit matrix :math:`[A]` of this operator.

        Where :meth:`inverse` / :attr:`H` / composition are
        *endofunctors* (``Op → Op``), ``as_matrix`` is the **functor OUT
        of the operator category** (``Op → Mat``) — the
        serialization boundary. Column :math:`j` is the operator applied
        to the :math:`j`-th basis element:

        .. math::

            [A]_{:,j} \;=\; \operatorname{ravel}\bigl(A\,e_j\bigr),
            \qquad e_j = \operatorname{unravel}(\delta_j),

        with basis elements enumerated in **C-order** over
        ``basis_shape`` and outputs raveled the same way — so for a
        group-leading ``(ng, 1)`` carrier, column ``j`` is the response
        to a unit source in group ``j``, and ``[A] @ x.ravel() ==
        A.apply(x).ravel()`` exactly. The matrix is
        ``(prod(output shape), prod(basis_shape))`` — RECTANGULAR when
        the operator is not endomorphic; the output dimension emerges
        from :meth:`apply` itself, never from declared metadata.

        This default is the apply-to-basis pattern. Structured operators
        MAY override with a direct assembly (the future per-octant
        sparse-triangular streaming assembly noted at
        ``sweep_graph.py:66`` — DEFERRED with its 3-D consumer;
        :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
        overrides with one batched LU backsolve). Until a sparse
        consumer exists, the return is a DENSE :class:`numpy.ndarray` —
        keyed by the operator's structural override, with dense the only
        realization built (defer-until-consumer).

        Parameters
        ----------
        basis_shape : tuple[int, ...], optional
            The element shape :meth:`apply` consumes. Default: derived
            from :attr:`domain` (``domain.shape``); REQUIRED explicitly
            for operators carrying no space (e.g. the meshless
            single-cell group operators, ``basis_shape=(ng, 1)``).
        max_dimension : int, optional
            The size gate: refuse (``MatrixTooLarge``) when
            ``prod(basis_shape) > max_dimension``. Default ``4096`` —
            a 4096² float64 is 134 MB and 4096 applications, generous
            for every dense-by-construction consumer (0-D energy
            spectra, CP ``[P]``) and prohibitive for none of them;
            a meshed SN full-field operator is refused by design.
            Per-call configurable — a RESOURCE knob, not structure
            (see :class:`MatrixTooLarge`).

        Raises
        ------
        ValueError
            No ``basis_shape`` given and the operator carries no
            :attr:`domain` to derive one from.
        MatrixTooLarge
            The resolved basis dimension exceeds ``max_dimension``.

        Notes
        -----
        **Honest scope**: the default serves ndarray-carrier operators
        (the energy/scattering blocks, small compositions, quadrature
        maps). Typed-carrier operators (``FullField`` SN composites)
        are not constructible from ndarray basis columns — and sit far
        above any sane gate; they stay matrix-free.

        **Assembly delegation**: when
        this operator is :func:`assemblable`, the densification routes
        through the structural sparse emission
        (``assemble().as_matrix(...)`` — same gate contract, same
        dimension checks) instead of :math:`n` probing applications.
        The probing loop is RETAINED as :meth:`_as_matrix_by_probing` —
        the fallback for assembly-less operators AND the permanent
        fuller-view oracle the probed≡assembled equivalence gates pin
        (an assembly bug must never be able to hide inside its own
        densification).
        """
        shape = _resolve_basis_shape(self, basis_shape)
        n = int(np.prod(shape)) if shape else 1
        if n > max_dimension:
            raise MatrixTooLarge(
                f"as_matrix on {type(self).__name__}: basis dimension "
                f"{n} (= prod{shape}) exceeds max_dimension="
                f"{max_dimension}. Materializing would commit ~"
                f"{8 * n * n / 1e6:.0f} MB and {n} operator "
                f"applications. Raise max_dimension= if this size is "
                f"intended, or keep the operator matrix-free."
            )
        if assemblable(self):
            # The assembly delegation: densified structural assembly, with
            # the assembled column dimension checked against the resolved
            # basis shape (SparseAssembledOperator.as_matrix enforces it).
            return self.assemble().as_matrix(
                basis_shape=shape, max_dimension=max_dimension,
            )
        return self._as_matrix_by_probing(shape)

    def _as_matrix_by_probing(self, shape: tuple[int, ...]) -> np.ndarray:
        r"""The apply-to-basis probing loop — the RETAINED fuller-view pathway.

        Column :math:`j` = ``apply(e_j)`` raveled (the pre-assembly
        ``as_matrix`` body, byte-identical). Kept as its own named
        method — NOT inlined in :meth:`as_matrix` — for two consumers:
        the delegation fallback (assembly-less operators), and the
        probed≡assembled equivalence gates, which MUST be able to force
        this pathway on an *assemblable* operator (otherwise
        ``as_matrix ≡ assemble().to_dense()`` is assembly compared with
        itself — vacuous). Size/shape gating is the caller's job
        (:meth:`as_matrix` resolves and gates before delegating here).
        """
        n = int(np.prod(shape)) if shape else 1
        columns = []
        for j in range(n):
            e_flat = np.zeros(n)
            e_flat[j] = 1.0
            column = self.apply(cast(Domain, e_flat.reshape(shape)))
            columns.append(np.asarray(column).ravel())
        return np.column_stack(columns)

    # ------------------------------------------------------------------
    # Repr
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        cls = type(self).__name__
        d = getattr(self, "domain", None)
        c = getattr(self, "codomain", None)
        d_name = repr(d.name) if d is not None else "'?'"
        c_name = repr(c.name) if c is not None else "'?'"
        # The two-axis surface, tokens present iff True.
        axes = "".join(
            f" {token}"
            for token, on in (
                ("invertible", getattr(self, "is_invertible", False)),
                ("adjointable", getattr(self, "is_adjointable", False)),
            )
            if on
        )
        return f"<{cls} domain={d_name} codomain={c_name}{axes}>"


# ───────────────────────────────────────────────────────────────────────
# Narrowing targets + checked bridges (#226 — Design C)
# ───────────────────────────────────────────────────────────────────────
#
# Each capability axis has THREE layers, and each layer carries the truth
# it alone can express:
#
#   1. the PREDICATE (``is_invertible`` / ``is_adjointable``) — runtime,
#      instance-accurate, value- and structure-aware; the polymorphic
#      override point every class defines its truth through;
#   2. the NARROWING TARGET below (``SupportsInverse``/``SupportsAdjoint``,
#      a Protocol EXTENDING :class:`LinearOperator`) — the static type a
#      guarded branch may treat the operand as;
#   3. the CHECKED BRIDGE (:func:`invertible` / :func:`adjointable`, a
#      PEP-647 ``TypeGuard``) — the ONE construct that converts the
#      runtime predicate into the static permission. You cannot obtain
#      the permission without executing the check (contrast the retired
#      ``cast(...)`` bridge, which asserted without checking).
#
# The Protocols are deliberately NOT ``runtime_checkable``: an
# ``isinstance`` check reads class-level method presence, which is
# class-uniform on composites (every ``OperatorSum`` defines
# ``apply_transpose`` even when a summand cannot transpose) and blind to
# value-dependent leaves (a zero-coefficient multiplier still has an
# ``inverse`` method). The bridge functions are the only sanctioned
# runtime→static conversion.


class SupportsInverse(LinearOperator[Domain, Codomain], Protocol):
    r"""Narrowing target: an operator whose :meth:`inverse` may be called.

    Extends :class:`LinearOperator`, so a branch narrowed by
    :func:`invertible` keeps the WHOLE algebra (``apply``, ``H``,
    composition dunders) alongside the licensed :meth:`inverse`. Never
    ``isinstance`` this (see the section comment); never annotate a
    parameter with it to *demand* invertibility — the static layer can
    certify only SPELLING (the method exists), never SOLVABILITY (the
    value-level predicate) — guard with :func:`invertible` instead.
    """

    def inverse(self) -> "LinearOperator[Codomain, Domain]": ...


class SupportsAdjoint(LinearOperator[Domain, Codomain], Protocol):
    r"""Narrowing target: an operator whose ``apply_transpose`` may be called.

    Extends :class:`LinearOperator`; the branch narrowed by
    :func:`adjointable` may call the raw Euclidean transpose verb
    ``apply_transpose`` (the realization the metric-aware
    :attr:`~LinearOperator.H` wrapper delegates to — two DIFFERENT
    objects: :math:`T^{\mathsf T}` vs :math:`G^{-1}T^{\mathsf T}G`).
    ``.H`` itself needs no narrowing — it lives on the base with an
    eager :class:`MissingAdjoint` gate.
    """

    def apply_transpose(self, x: Codomain, /) -> Domain: ...


class SupportsAssembly(LinearOperator[Domain, Codomain], Protocol):
    r"""Narrowing target: an operator whose :meth:`assemble` may be called.

    Extends :class:`LinearOperator`; the branch narrowed by
    :func:`assemblable` may call the structural sparse emission
    :meth:`assemble` — the ASSEMBLY-axis sibling of
    :class:`SupportsInverse` / :class:`SupportsAdjoint`, with the same
    discipline: never ``isinstance`` this (class-level method presence
    is class-uniform on composites — every ``OperatorSum`` defines
    ``assemble`` even when a summand cannot emit); never annotate a
    parameter with it to *demand* assemblability — guard with
    :func:`assemblable` instead.
    """

    def assemble(self) -> "SparseAssembledOperator": ...


def invertible(
    op: "LinearOperator[Domain, Codomain]",
) -> "TypeGuard[SupportsInverse[Domain, Codomain]]":
    r"""Checked bridge: narrow ``op`` to :class:`SupportsInverse` iff invertible.

    The runtime check and the static permission are ONE construct: a
    branch guarded by this function may call ``op.inverse()`` with no
    ``cast`` — and deleting the guard un-narrows the call, so CLI
    pyright REDs (the guard is type-load-bearing).

    Deliberately ``TypeGuard``, NOT ``TypeIs``: the predicate is
    VALUE-dependent (a zero-coefficient multiplier structurally *has*
    ``inverse()`` while reporting ``False``), so only the one-directional
    promise is honest — ``True`` licenses the call; ``False`` makes no
    static claim. A free function because PEP 647 narrowing applies only
    through a call expression and a method form narrows its first
    *explicit* argument, never ``self`` — no property spelling exists.
    Guard at ``LinearOperator``-typed sites only: ``TypeGuard`` REPLACES
    (does not intersect) the declared type, so guarding an
    already-concrete operand would widen it.
    """
    return op.is_invertible


def adjointable(
    op: "LinearOperator[Domain, Codomain]",
) -> "TypeGuard[SupportsAdjoint[Domain, Codomain]]":
    r"""Checked bridge: narrow ``op`` to :class:`SupportsAdjoint` iff adjointable.

    The adjoint-axis twin of :func:`invertible` — same one-directional
    ``TypeGuard`` semantics, same free-function necessity, same
    guard-at-``LinearOperator``-typed-sites discipline. Licenses the raw
    ``apply_transpose`` realization verb in composer law bodies
    (:math:`(A+B)^{\mathsf T} = A^{\mathsf T} + B^{\mathsf T}`) and
    gates the eager :attr:`~LinearOperator.H` construction.
    """
    return op.is_adjointable


def assemblable(
    op: "LinearOperator[Domain, Codomain]",
) -> "TypeGuard[SupportsAssembly[Domain, Codomain]]":
    r"""Checked bridge: narrow ``op`` to :class:`SupportsAssembly` iff assemblable.

    The assembly-axis member of the bridge family — same one-directional
    ``TypeGuard`` semantics, same free-function necessity, same
    guard-at-``LinearOperator``-typed-sites discipline as
    :func:`invertible` / :func:`adjointable`. Licenses the structural
    :meth:`~SupportsAssembly.assemble` emission in the composer
    homomorphism-law bodies (Sum → ``+``, Product → ``@``, Scaled →
    scalar ``*``) and the
    :meth:`~LinearOperator.as_matrix` densification delegation (R2).
    """
    return op.is_assemblable


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

    Construction is gated EAGERLY by :meth:`LinearOperator.adjoint`:
    only an :func:`adjointable`-narrowed operator
    reaches this constructor, so ``inner`` is statically a
    :class:`SupportsAdjoint` and :meth:`apply`'s delegation to
    ``inner.apply_transpose`` is typed — there is no lazy capability
    gate left to fail at call time. The reverse direction
    (apply_transpose on the adjoint = apply on the inner) is not
    needed by any current consumer in 9.6 and is deferred —
    :meth:`apply_transpose` raises :class:`NotImplementedError`
    until a consumer demands it.
    """

    def __init__(self, inner: "SupportsAdjoint[Domain, Codomain]") -> None:
        self.inner = inner
        # The G-adjoint transposes the 2×2 block matrix (A_bs ↔ A_sb^T),
        # which preserves WHICH blocks are touched — so the role is the
        # inner operator's role: ``L.H`` is FULL, ``B.H`` is BOUNDARY,
        # ``C.H`` is BULK.
        self.block_role = getattr(inner, "block_role", None)
        # The G-adjoint also preserves which SYSTEMS are touched: A_AB.H
        # (bulk→ray) is still COUPLED, A_BB.H still System B.
        self.system_role = getattr(inner, "system_role", None)

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        # Adjoint of A: V → W is A.H: W → V — domain swaps with inner.codomain.
        return getattr(self.inner, "codomain", None)

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        return getattr(self.inner, "domain", None)

    def apply(self, y: Codomain) -> Domain:
        # No capability gate here: the eager MissingAdjoint raise in
        # LinearOperator.adjoint() is the ONLY entrance (no direct
        # constructions exist — verified), so ``inner`` always
        # carries a working apply_transpose by the time apply runs.
        # Hilbert-adjoint action:
        #   (A^* y)_V = G_V⁺ ⊙ apply_transpose(G_W ⊙ y)
        # On the adjoint wrapper, ``codomain`` is the inner operator's
        # domain (V) and ``domain`` its codomain (W). The metric application
        # is delegated to the function space's :meth:`~FunctionSpace.apply_metric`
        # / :meth:`~FunctionSpace.apply_inverse_metric` so that
        # the SAME wrapper serves BOTH a flat-ndarray metric (e.g. the
        # spherical-harmonic ``(L+1, 2L+1)`` leading-axis metric) AND a
        # composite bulk ⊕ trace metric on a structured ``FullField`` (the
        # direct-sum space applies a per-block metric, with a pseudo-inverse on
        # the singular partial-current trace). The space owns the metric; the
        # adjoint wrapper is metric-representation-agnostic.
        inner_codomain = getattr(self.inner, "codomain", None)
        inner_domain = getattr(self.inner, "domain", None)
        z = inner_codomain.apply_metric(y) if inner_codomain is not None else y
        result = self.inner.apply_transpose(z)
        if inner_domain is not None:
            result = inner_domain.apply_inverse_metric(result)
        return result

    def apply_transpose(self, x: Domain) -> Codomain:
        raise NotImplementedError(
            "apply_transpose on an _AdjointOperator wrapper is not "
            "supported in 9.6; if a consumer needs it, take the adjoint "
            "of the original inner operator's transpose directly."
        )

    @property
    def is_invertible(self) -> bool:
        r"""Whether the adjoint's inverse operator exists — the swap law (#280).

        The inverse of the adjoint IS the adjoint of the inverse:
        :math:`(A^{*})^{-1} = (A^{-1})^{*}`. Honest
        iff the inner :math:`A` is invertible AND its inverse operator is
        adjointable (so ``.H`` on that inverse is well-posed) — spelled
        generally over the inner, no leaf specifics. :func:`invertible`
        narrows ``self.inner`` to :class:`SupportsInverse` for the RHS call;
        the ``and`` short-circuits so ``inner.inverse()`` is never built for a
        non-invertible inner.
        """
        return invertible(self.inner) and adjointable(self.inner.inverse())

    def inverse(self) -> "LinearOperator[Domain, Codomain]":
        r"""The inverse of the adjoint = the adjoint of the inverse (#280).

        :math:`(A^{*})^{-1} = (A^{-1})^{*}` — the operator-algebra swap law,
        an OBJECT IDENTITY here (not a computed numerical equivalence): this
        wrapper IS ``A.H``, so its inverse routes to ``A.inverse().H``. Gated
        by :attr:`is_invertible` (``A`` invertible and ``A.inverse()``
        adjointable). The metric adjoint-solve
        :math:`A^{-1\,*} b = G_V^{+}\,\mathrm{apply\_transpose}(G_W\,b)` then
        falls out of :meth:`apply` (which already routes
        ``inner.apply_transpose``) FOR FREE — no ``_AdjointOperator.solve`` /
        no metric code enters the sweep.
        """
        if not invertible(self.inner):
            raise NotInvertible(
                f"_AdjointOperator.inverse(): the inner "
                f"{type(self.inner).__name__} is not invertible, so the "
                f"adjoint-inverse swap law (A.H).inverse() = (A.inverse()).H "
                f"does not apply (is_invertible is False)."
            )
        inner_inverse = self.inner.inverse()
        if not adjointable(inner_inverse):
            # The swap law needs the inverse to be adjointable so ``.H`` exists
            # — matches :attr:`is_invertible`'s second clause. Raise
            # NotInvertible (NOT MissingAdjoint from the ``.H`` below): the
            # adjoint-INVERSE is what is absent here.
            raise NotInvertible(
                f"_AdjointOperator.inverse(): the inner's inverse "
                f"{type(inner_inverse).__name__} is not adjointable, so "
                f"(A.inverse()).H does not exist — the swap law needs an "
                f"adjointable inverse (e.g. an SN SweepOperator, #280 2.5c). "
                f"is_invertible is False."
            )
        return inner_inverse.H


class OperatorSum(
    LinearOperator[Domain, Codomain],
    Generic[Domain, Codomain, SummandA, SummandB],
):
    r"""Sum of two linear operators: :math:`(A + B)\,x = A\,x + B\,x`.

    Generic over its SUMMAND types: a named composition subclass
    pins them — ``StreamingCollisionOperator = OperatorSum["FullField",
    "FullField", StreamingOperator, MultiplicationOperator]`` — so its
    leg accessors are typed by construction, no casts. The PEP-696
    defaults (``LinearOperator[Domain, Codomain]``) keep every
    ``OperatorSum[D, C]`` / bare spelling valid; the legs are covariant,
    read-only properties so a pinned composition upcasts to the
    defaulted spelling (the ``__add__`` contract).

    Structural closure laws (realized in the method bodies; the
    predicates are the matching advertisements):

    * ``apply`` requires **both** operands to act — guarded eagerly at
      construction (``TypeError``), never at the first call.
    * Invertibility does NOT propagate operand-wise: there is no
      algorithm for :math:`(A + B)^{-1}` from :math:`A^{-1}` and
      :math:`B^{-1}` alone — Sherman–Morrison–Woodbury applies only
      under low-rank structure (which the boundary block B has — rank
      ≤ N/2 per face, Issue #300 — and the bulk C, L do not).  What a
      sum DOES have, when its LEADING
      (left-spine head) term is invertible, is a
      preconditioned-SPLITTING inverse: :meth:`inverse` returns a
      :class:`~orpheus.numerics.green_operator.GreenOperator` iterating
      :math:`x_{n+1} = A^{-1}(q - B\,x_n)`.  A generic sum carries no
      ``solve`` verb — solving with it IS applying that inverse OBJECT
      (the sweep-invertible ``(L+C)`` subclass overrides with its own
      direct-sweep ``solve``).  See :attr:`is_invertible` for the
      canonical-ordering contract.
    * ``apply_transpose`` requires **both** operands to transpose
      (:math:`(A + B)^T = A^T + B^T`) — guarded in the verb body with
      :class:`MissingAdjoint`; :attr:`is_adjointable` is the recursion.

    Raises
    ------
    TypeError
        If either operand has no callable ``apply`` at construction
        time. Catch the failure here, not at the first ``apply`` call,
        so downstream Krylov consumers don't see a stub failure
        mid-iteration.
    """

    def __init__(self, a: SummandA, b: SummandB) -> None:
        if not callable(getattr(a, "apply", None)):
            raise TypeError(
                f"OperatorSum requires apply on both operands; left "
                f"operand {type(a).__name__} lacks 'apply'."
            )
        if not callable(getattr(b, "apply", None)):
            raise TypeError(
                f"OperatorSum requires apply on both operands; right "
                f"operand {type(b).__name__} lacks 'apply'."
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
        self._a: Final = a
        self._b: Final = b
        # Block role DERIVED from the operands: the sum touches the union
        # of the blocks its summands touch, so ``(L+C)`` and the whole
        # ``(L+C-S-F-B)`` loss carry FULL by construction (no hand-stamp).
        self.block_role = _join_block_roles(
            getattr(a, "block_role", None), getattr(b, "block_role", None),
        )
        # System role DERIVED the same way — the sum touches the union of the
        # systems its summands touch (``A ⊔ B = COUPLED``); the two-system
        # analogue of the block-role join.
        self.system_role = _join_system_roles(
            getattr(a, "system_role", None), getattr(b, "system_role", None),
        )

    @property
    def a(self) -> SummandA:
        """The left summand (read-only — covariant leg typing)."""
        return self._a

    @property
    def b(self) -> SummandB:
        """The right summand (read-only — covariant leg typing)."""
        return self._b

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
        # (A+B)^T = A^T + B^T — the guard-narrow licenses the operand
        # calls (Design C: the runtime check IS the static permission).
        if not adjointable(self.a) or not adjointable(self.b):
            raise MissingAdjoint(
                f"OperatorSum.apply_transpose requires both summands to "
                f"transpose ((A+B)^T = A^T + B^T); got "
                f"{type(self.a).__name__} / {type(self.b).__name__} with "
                f"is_adjointable {self.a.is_adjointable} / "
                f"{self.b.is_adjointable}."
            )
        return self.a.apply_transpose(x) + self.b.apply_transpose(x)

    @property
    def is_adjointable(self) -> bool:
        # (A+B)^T = A^T + B^T (the law in :meth:`apply_transpose`) — the
        # sum is adjointable iff BOTH summands are.
        return self.a.is_adjointable and self.b.is_adjointable

    @property
    def is_invertible(self) -> bool:
        r"""``True`` iff the LEADING (left-spine head) term is invertible.

        There is no operand-wise law for a sum's inverse — but a sum
        whose leading term :math:`A` is invertible CAN produce its
        inverse OPERATOR: the preconditioned-splitting
        :class:`~orpheus.numerics.green_operator.GreenOperator`
        (:math:`x_{n+1} = A^{-1}(q - B\,x_n)`).
        The recursion ``self.a.is_invertible`` designates the left-spine
        head as the splitting's preconditioner — the CANONICAL-ORDERING
        contract: spell the invertible operator FIRST (``A - S``,
        mirroring the ``L + C`` fusion rule of
        :meth:`~orpheus.sn.operators.streaming.StreamingOperator.__add__`,
        #261).  Whether the splitting CONVERGES is a spectral
        (value-level) property no construction-time predicate can read —
        a divergent split raises
        :class:`~orpheus.numerics.green_operator.ConvergenceFailure`
        loudly at apply time, never a silent wrong answer.  (The
        sweep-invertible
        :class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator`
        subclass shadows this by MRO with its own ``True`` +
        direct-sweep :meth:`inverse` — the type-as-structure dispatch.)
        """
        return self.a.is_invertible

    def inverse(self) -> "LinearOperator[Codomain, Domain]":
        r"""Return the preconditioned-splitting inverse — a
        :class:`~orpheus.numerics.green_operator.GreenOperator`.

        The annotation is the factory's honest STATIC face — "an inverse
        operator on the swapped spaces" — because subclass overrides
        return their own structure's inverse (the sweep-invertible
        composite returns a ``SweepOperator``; type-as-structure) and the
        family members are siblings, not a hierarchy.

        Late import: ``green_operator`` is a LEAF module wrapping the
        iteration drivers, which import THIS module — the same one-way
        late-import pattern as
        :meth:`~orpheus.sn.operators.streaming.StreamingCollisionOperator.inverse`
        → ``SweepOperator``.
        """
        from orpheus.numerics.green_operator import GreenOperator

        return GreenOperator(self)

    # NO ``solve`` on a generic sum: its inverse action is DRIVER-realized
    # (the GreenOperator), not a substrate verb — solving is
    # ``.inverse().apply(b)``. The sweep-invertible ``(L+C)`` subclass
    # overrides with its own direct sweep ``solve`` (streaming.py).

    @property
    def is_assemblable(self) -> bool:
        # [A+B] = [A] + [B] (the law in :meth:`assemble`) — the sum
        # assembles iff BOTH summands emit.
        return self.a.is_assemblable and self.b.is_assemblable

    def assemble(self) -> "SparseAssembledOperator":
        r"""Return :math:`[A+B] = [A] + [B]` — the additive homomorphism law.

        The assembly functor is additive-monoidal (the ``as_matrix``
        docstring's ``Op → Mat``, sparse carrier): a sum's structural
        emission is the SPARSE SUM of its summands' emissions —
        realized by the carrier's own CSR addition, never a re-walk of
        the stencils. The guard-narrow licenses the operand calls
        (Design C) and raises the assembly-axis refusal eagerly.
        """
        if not assemblable(self.a) or not assemblable(self.b):
            raise MissingAssembly(
                f"OperatorSum.assemble requires both summands to emit "
                f"([A+B] = [A] + [B]); got {type(self.a).__name__} / "
                f"{type(self.b).__name__} with is_assemblable "
                f"{self.a.is_assemblable} / {self.b.is_assemblable}."
            )
        from orpheus.numerics.assembled_operator import SparseAssembledOperator

        return SparseAssembledOperator(
            self.a.assemble().matrix + self.b.assemble().matrix,
            domain=self.domain,
            codomain=self.codomain,
        )


class OperatorProduct(
    LinearOperator[Domain, Codomain],
    Generic[Domain, Codomain, FactorA, FactorB],
):
    r"""Composition of two linear operators: :math:`(A\,B)\,x = A(B\,x)`.

    Generic over its FACTOR types (as :class:`OperatorSum` over
    its summands): a named composition pins them — ``WindowedSweep =
    OperatorProduct["FullField", "TimedFullField", BulkAnalysisOperator,
    SweepOperator]`` — so its factor accessors are typed by
    construction. The PEP-696 defaults keep ``OperatorProduct[D, C]``
    valid; the legs are covariant read-only properties. The
    intermediate-space coupling (``A.domain == B.codomain``) is the
    RUNTIME guard below — the leg parameters do not re-encode it.

    Structural closure laws (method bodies; predicates advertise):

    * ``apply`` requires **both** operands (function composition) —
      guarded eagerly at construction (``TypeError``).
    * Invertibility propagates iff **both** factors are invertible,
      with the order reversed: :math:`(A\,B)^{-1} = B^{-1}\,A^{-1}`.
      The product IS a wrap-delegate conformer — :meth:`solve` is its
      native realization verb, re-routed through the factors' CANONICAL
      surface (``.inverse().apply``) so factor kinds whose own
      ``solve`` retired (algebra-closed permutations/scalings, Green-
      invertible sums) compose without one.
    * ``apply_transpose`` requires **both**, order reversed
      (:math:`(A\,B)^T = B^T\,A^T`) — :class:`MissingAdjoint` in the
      verb body; :attr:`is_adjointable` is the recursion.

    Raises
    ------
    TypeError
        If either operand has no callable ``apply`` at construction.
    """

    def __init__(self, a: FactorA, b: FactorB) -> None:
        # ``A @ B``: ``B`` maps the input ``V`` to the intermediate,
        # ``A`` maps the intermediate to the output ``W`` — so the
        # product is honestly ``V → W``, the intermediate guarded below.
        if not callable(getattr(a, "apply", None)):
            raise TypeError(
                f"OperatorProduct requires apply on both operands; "
                f"left operand {type(a).__name__} lacks 'apply'."
            )
        if not callable(getattr(b, "apply", None)):
            raise TypeError(
                f"OperatorProduct requires apply on both operands; "
                f"right operand {type(b).__name__} lacks 'apply'."
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
        self._a: Final = a
        self._b: Final = b

    @property
    def a(self) -> FactorA:
        """The left factor ``A`` of ``A @ B`` (read-only — covariant leg)."""
        return self._a

    @property
    def b(self) -> FactorB:
        """The right factor ``B`` of ``A @ B`` (read-only — covariant leg)."""
        return self._b

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
        r"""Solve :math:`(AB)\,x = b` — :math:`B^{-1}(A^{-1}\,b)`, factor-wise.

        The product's native realization verb (the wrap-delegate family
        wraps it: :meth:`inverse` returns ``InverseOperator(self)`` whose
        ``apply`` delegates here). The recursion goes
        through each factor's CANONICAL surface — ``.inverse().apply`` —
        not a factor ``solve``: bit-identical for every factor kind (the
        inverse objects delegate to the same realizations) and total
        over the kinds whose own ``solve`` retired (a permutation's
        inverse is a first-class forward; a Green-invertible sum's is
        the GreenOperator). The guard-narrow licenses the calls and
        raises the value-dependent refusal eagerly.
        """
        if not invertible(self.a) or not invertible(self.b):
            raise NotInvertible(
                f"OperatorProduct.solve requires both factors to be "
                f"invertible ((AB)^{{-1}} = B^{{-1}}A^{{-1}}); got "
                f"{type(self.a).__name__} / {type(self.b).__name__} with "
                f"is_invertible {self.a.is_invertible} / "
                f"{self.b.is_invertible}."
            )
        return self.b.inverse().apply(self.a.inverse().apply(b_vec))

    def apply_transpose(self, x: Codomain, /) -> Domain:
        # (AB)^T = B^T A^T — the guard-narrow licenses the operand calls.
        if not adjointable(self.a) or not adjointable(self.b):
            raise MissingAdjoint(
                f"OperatorProduct.apply_transpose requires both factors "
                f"to transpose ((AB)^T = B^T A^T); got "
                f"{type(self.a).__name__} / {type(self.b).__name__} with "
                f"is_adjointable {self.a.is_adjointable} / "
                f"{self.b.is_adjointable}."
            )
        return self.b.apply_transpose(self.a.apply_transpose(x))

    @property
    def is_invertible(self) -> bool:
        # (AB)^{-1} = B^{-1} A^{-1} (the law in :meth:`solve`) — the product
        # is invertible iff BOTH factors are.
        return self.a.is_invertible and self.b.is_invertible

    def inverse(self) -> "InverseOperator":
        r"""Return :math:`(AB)^{-1}` — the generic family member wrapping this product.

        The functoriality law :math:`(AB)^{-1} = B^{-1}A^{-1}` holds
        BEHAVIORALLY through the wrapper: ``inverse().apply(q)`` delegates
        to this product's own :meth:`solve` ``= b.solve(a.solve(q))``.
        What the family wrapper adds is the CONTRACT (#285): a raw
        ``OperatorProduct`` of inverses carries no ``initial_guess``
        keyword, so a driver seeding it raised ``TypeError`` at iteration
        time; :class:`InverseOperator` carries the family's canonical
        seeded ``apply`` (accept-and-ignore — the solve path never
        threaded seeds either, so behavior is unchanged) and every
        ``.inverse()`` in the system now returns a seeded-apply
        conformer. The involution is object identity —
        ``(A@B).inverse().inverse() is (A@B)`` via the mixin. The factors
        stay reachable as ``.inner.a`` / ``.inner.b``.

        (Contrast the ALGEBRA-CLOSED inverses — a
        :class:`PermutationOperator`'s inverse IS a permutation, an
        :class:`IdentityOperator` is self-inverse, a
        :class:`ScaledOperator`'s is a scaled inverse: those inverses
        are first-class FORWARD operators in their own closed structure,
        the other kind of inverse, and stay unwrapped.)
        """
        if not self.is_invertible:
            raise NotInvertible(
                "OperatorProduct.inverse requires both factors to be "
                "invertible ((AB)^{-1} = B^{-1}A^{-1})."
            )
        return InverseOperator(self)

    @property
    def is_adjointable(self) -> bool:
        # (AB)^T = B^T A^T (the law in :meth:`apply_transpose`) — adjointable
        # iff BOTH factors are.
        return self.a.is_adjointable and self.b.is_adjointable

    @property
    def is_assemblable(self) -> bool:
        # [AB] = [A] @ [B] (the law in :meth:`assemble`) — the product
        # assembles iff BOTH factors emit.
        return self.a.is_assemblable and self.b.is_assemblable

    def assemble(self) -> "SparseAssembledOperator":
        r"""Return :math:`[AB] = [A]\,[B]` — the multiplicative homomorphism law.

        The composition's structural emission is the SPARSE PRODUCT of
        its factors' emissions (dimension compatibility enforced by the
        carrier's own matmul). Same eager guard-narrow discipline as
        :meth:`OperatorSum.assemble`.
        """
        if not assemblable(self.a) or not assemblable(self.b):
            raise MissingAssembly(
                f"OperatorProduct.assemble requires both factors to emit "
                f"([AB] = [A][B]); got {type(self.a).__name__} / "
                f"{type(self.b).__name__} with is_assemblable "
                f"{self.a.is_assemblable} / {self.b.is_assemblable}."
            )
        from orpheus.numerics.assembled_operator import SparseAssembledOperator

        return SparseAssembledOperator(
            self.a.assemble().matrix @ self.b.assemble().matrix,
            domain=self.domain,
            codomain=self.codomain,
        )


class ScaledOperator(
    LinearOperator[Domain, Codomain],
    Generic[Domain, Codomain, ScaledOperand],
):
    r"""Scalar multiple of a linear operator: :math:`(\alpha L)\,x = \alpha\,(L\,x)`.

    Generic over its OPERAND type (as the other composition
    wrappers over their legs): ``ScaledOperator["FullField",
    "FullField", SNMaskedBoundaryOperator]`` reads ``.op`` as the masked
    boundary leaf — the ``-1·B_lower`` leg of the G-S splitting needs no
    cast. The PEP-696 default keeps ``ScaledOperator[D, C]`` valid; the
    operand is a covariant read-only property.

    Scaling (:math:`\alpha \neq 0`, caught at composition time) passes
    both structural axes through unchanged: the operand's
    invertibility/adjointability ARE the scaled operator's, and the
    algebra is closed — :meth:`inverse` is a
    :class:`ScaledOperator` (:math:`(\alpha L)^{-1} = (1/\alpha)L^{-1}`)
    and the transpose scales (:math:`(\alpha L)^T = \alpha L^T`). No
    ``solve`` verb: an algebra-closed inverse is a first-class forward,
    so solving is ``.inverse().apply(b)``.
    """

    def __init__(self, scalar: float, op: ScaledOperand) -> None:
        if not callable(getattr(op, "apply", None)):
            raise TypeError(
                f"ScaledOperator requires apply on its operand; "
                f"{type(op).__name__} lacks 'apply'."
            )
        if scalar == 0.0:
            # Zero scaling is a degenerate case: the result behaves as
            # a ZeroOperator (singular, structurally), not as the
            # underlying operator. The user should construct
            # ZeroOperator explicitly to make the intent clear.
            raise ValueError(
                "ScaledOperator with zero scalar is degenerate; "
                "use ZeroOperator explicitly."
            )
        self.scalar = float(scalar)
        self._op: Final = op
        # Scaling preserves which blocks the action touches.
        self.block_role = getattr(op, "block_role", None)
        # Scaling preserves which systems the action touches, too.
        self.system_role = getattr(op, "system_role", None)

    @property
    def op(self) -> ScaledOperand:
        """The scaled operand ``L`` (read-only — covariant leg typing)."""
        return self._op

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        return getattr(self.op, "domain", None)

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        return getattr(self.op, "codomain", None)

    def apply(self, x: Domain, /, *extra, **kwextra) -> Codomain:
        return self.scalar * self.op.apply(x, *extra, **kwextra)

    # NO ``solve``: the inverse is ALGEBRA-CLOSED —
    # :meth:`inverse` returns the first-class forward
    # ``ScaledOperator(1/α, op.inverse())`` — so there is no wrapped
    # realization verb to keep; solving is ``.inverse().apply(b)``.

    def apply_transpose(self, x: Codomain, /, *extra, **kwextra) -> Domain:
        # (αL)^T = α L^T — the guard-narrow licenses the operand call.
        if not adjointable(self.op):
            raise MissingAdjoint(
                f"ScaledOperator.apply_transpose requires an adjointable "
                f"operand ((αL)^T = αL^T); {type(self.op).__name__}."
                f"is_adjointable is False."
            )
        return self.scalar * self.op.apply_transpose(x, *extra, **kwextra)

    @property
    def is_invertible(self) -> bool:
        # (αL)^{-1} = (1/α) L^{-1} — α ≠ 0 is enforced at construction, so
        # the scaled operator is invertible iff the operand is.
        return self.op.is_invertible

    def inverse(self) -> "ScaledOperator[Codomain, Domain]":
        r"""Return :math:`(\alpha L)^{-1} = (1/\alpha)\,L^{-1}`.

        The natural structural inverse: a scaled operator's inverse IS a
        scaled operator — on the SWAPPED carriers (an inverse maps the
        forward's codomain back to its domain), so the return type is
        ``ScaledOperator[Codomain, Domain]``. ``1/α`` is exact whenever
        ``α`` is a power of two
        (the dominant −1.0 case); the action is bit-identical to
        :meth:`solve` given the operand's own ``inverse().apply ≡ solve``
        identity (both spell ``(1/α) * op_solve(b)``).
        """
        if not invertible(self.op):
            raise NotInvertible(
                "ScaledOperator.inverse requires an invertible operand "
                "((αL)^{-1} = (1/α)L^{-1})."
            )
        return ScaledOperator(1.0 / self.scalar, self.op.inverse())

    @property
    def is_adjointable(self) -> bool:
        # (αL)^T = α L^T — scaling preserves adjointability.
        return self.op.is_adjointable

    @property
    def is_assemblable(self) -> bool:
        # [αL] = α[L] (the law in :meth:`assemble`) — scaling preserves
        # assemblability.
        return self.op.is_assemblable

    def assemble(self) -> "SparseAssembledOperator":
        r"""Return :math:`[\alpha L] = \alpha\,[L]` — the scalar homomorphism law.

        The scaled emission is the carrier's own scalar multiply of the
        operand's emission. Same eager guard-narrow discipline as the
        other composer laws.
        """
        if not assemblable(self.op):
            raise MissingAssembly(
                f"ScaledOperator.assemble requires an assemblable operand "
                f"([αL] = α[L]); {type(self.op).__name__}.is_assemblable "
                f"is False."
            )
        from orpheus.numerics.assembled_operator import SparseAssembledOperator

        return SparseAssembledOperator(
            self.scalar * self.op.assemble().matrix,
            domain=self.domain,
            codomain=self.codomain,
        )


class IdentityOperator(LinearOperator[Domain]):
    r"""The identity operator :math:`I\,x = x`.

    Both axes hold trivially — :math:`I^{-1} = I` and :math:`I^T = I` —
    and both are ALGEBRA-CLOSED: :meth:`inverse` returns this very
    instance, so there is no ``solve`` verb to keep (solving with the
    identity IS applying its inverse, itself).
    """

    def apply(self, x: Domain, /) -> Domain:
        return x

    def apply_transpose(self, x: Domain, /) -> Domain:
        return x

    @property
    def is_invertible(self) -> bool:
        return True  # I^{-1} = I

    def inverse(self) -> "IdentityOperator[Domain]":
        r"""Return :math:`I^{-1} = I` — this very instance (stateless)."""
        return self

    @property
    def is_adjointable(self) -> bool:
        return True  # I^T = I


class ZeroOperator(LinearOperator[Domain, Codomain]):
    r"""The zero operator :math:`0\,x = 0`.

    Has ``apply`` (returns the zero of its CODOMAIN) and
    ``apply_transpose`` (also zero), and is STRUCTURALLY non-invertible
    — the singular map par excellence — so it declares no ``inverse()``
    at all: misuse is a static error, the honest surface for a type
    whose inverse does not exist mathematically (Design C; forcing a
    raising stub would be the harmful-stub anti-pattern this module is
    designed against).

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
    ``ZeroOperator(codomain_zero=…)`` returning a source-typed
    (:class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`)
    zero, so the typed RHS ``S.apply(ψ) + F.apply(ψ) + q_ext`` and the
    Krylov matvec ``L.apply − S.apply − F.apply`` stay CLOSED source-typed
    sums (a flux-echoing zero would hit the cross-class gate). Formal
    operator codomain typing is issue #208; ``codomain_zero`` is the
    pre-#208 hook that keeps the zero operator honest about what space it
    maps into. Since #276 A4 the hook is SYMMETRIC: ``transpose_zero``
    supplies the typed zero the TRANSPOSE emits (the domain's dual-role
    zero under duality typing) — first exercised by the coupled fission
    grid's (B, B) slot, whose whole grid the daggered eigen posing
    transposes. Without either hook both directions stay the endomorphic
    ``0.0 * x`` echo, bit-identical to the pre-A4 behaviour.
    """

    def __init__(
        self,
        codomain_zero: "Callable[[Domain], Codomain] | None" = None,
        transpose_zero: "Callable[[Codomain], Domain] | None" = None,
    ) -> None:
        self._codomain_zero = codomain_zero
        self._transpose_zero = transpose_zero

    def apply(self, x: Domain, /) -> Codomain:
        if self._codomain_zero is not None:
            return self._codomain_zero(x)
        # Endomorphic default (domain == codomain ⟹ ``W`` is ``V``): the
        # zero of the codomain IS ``0.0 * x``. ``cast`` is the PEP-484
        # bridge for the one genuinely-untypeable spot — this branch is
        # reached ONLY when no ``codomain_zero`` was supplied, i.e. when
        # the operator is endomorphic and ``W == V``.
        return cast("Codomain", 0.0 * x)

    def apply_transpose(self, x: "Codomain", /) -> Domain:
        # The transpose of a role-MAPPING zero slot emits a typed zero
        # too: under duality typing (#276 A4) the transpose consumes the
        # codomain-cotangent and emits the domain-cotangent, whose class
        # is the domain's DUAL role — supplied by ``transpose_zero``
        # exactly as ``codomain_zero`` supplies the forward's.  First
        # exercised by the coupled fission grid's (B, B) slot (#276 A4 —
        # the daggered eigen posing transposes the whole grid), which
        # ended the pre-#208 "transpose of the zero slot is not
        # exercised" era.
        if self._transpose_zero is not None:
            return self._transpose_zero(x)
        # Endomorphic default — the zero echo (see ``apply``).
        return cast("Domain", 0.0 * x)

    @property
    def is_adjointable(self) -> bool:
        return True  # 0^T = 0

    # is_invertible inherits the base ``False`` — the zero map is singular.


class _WrappedForward(Protocol):
    r"""The MINIMAL structural contract the wrap-delegate back-half consumes.

    Exactly what :class:`InverseWrapMixin` itself reads of its wrapped
    forward :math:`A`: the function-space pair the inverse SWAPS
    (``domain``/``codomain``) and the forward matvec its ``solve``
    un-inverts through (``apply``). Nothing more — the
    :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
    sibling inverts the MATERIALIZATION and never touches ``inner.solve``
    or ``inner.is_invertible`` (it reads values, not structure), so the
    minimal contract is these three members only (the tighter
    :class:`_InvertibleForward` bound fits the solve-backed siblings).

    Each sibling NARROWS ``_ForwardT`` to what its own ctor guard and
    algorithm need: :class:`InverseOperator` to
    :class:`_InvertibleForward`;
    :class:`~orpheus.sn.operators.sweep_operator.SweepOperator` to the
    schedule-triangular union;
    :class:`~orpheus.numerics.green_operator.GreenOperator` to
    :class:`OperatorSum`; ``MatrixInverseOperator`` to
    :class:`LinearOperator` (its guard needs :meth:`~LinearOperator.as_matrix`).
    """

    @property
    def domain(self) -> Optional["FunctionSpace"]: ...

    @property
    def codomain(self) -> Optional["FunctionSpace"]: ...

    def apply(self, x: Any, /) -> Any: ...


class _InvertibleForward(_WrappedForward, Protocol):
    r"""A solve-backed invertible forward — :class:`InverseOperator`'s narrowing.

    Extends the family-minimal :class:`_WrappedForward` with the two
    members the GENERIC sibling consumes: ``is_invertible`` (its ctor
    guard is the leaf's own value check) and :meth:`solve` — the
    forward's NATIVE inverse-action realization, the permanent face the
    family wrapper delegates through, so the verb lives exactly on the
    conformers the family wraps (value leaves, the product, the sweep
    composites). Delegating through one contract keeps the inverse OBJECT
    and the realization on ONE body (``coding-elegance`` Pattern 2 — no
    reciprocal twin path that could drift by a rounding).
    """

    @property
    def is_invertible(self) -> bool: ...

    def solve(self, b: Any, /) -> Any: ...


_ForwardT = TypeVar("_ForwardT", bound=_WrappedForward)


class InverseWrapMixin(Generic[_ForwardT], metaclass=ABCMeta):
    r"""The wrap-delegate back-half shared by every inverse-family sibling.

    An inverse operator in this codebase is a thin typed wrapper around
    its own FORWARD operator :math:`A` (:attr:`inner`): the wrapper's
    :meth:`apply` realizes :math:`A^{-1}` by the sibling's algorithm, and
    everything else is delegation — the byte-identical back-half shared
    by every inverse-family sibling (:class:`InverseOperator`,
    :class:`~orpheus.sn.operators.sweep_operator.SweepOperator`,
    :class:`~orpheus.numerics.green_operator.GreenOperator`), extracted
    once the third sibling appeared (defer-until-≥2):

    * :attr:`domain` / :attr:`codomain` — the SWAP of the forward's: an
      inverse maps the forward's codomain back to its domain.
    * :meth:`solve` — solving :math:`A^{-1}\,y = b` IS applying
      :math:`A`: the forward matvec ``inner.apply``, delegated.
    * ``is_invertible`` is ``True`` and :meth:`inverse` returns the
      wrapped forward ITSELF — the involution :math:`(A^{-1})^{-1} = A`
      holds by OBJECT IDENTITY, typed per sibling
      through ``_ForwardT``.

    **The canonical seeded-apply signature is part of the back-half**
    (#285, resolved STRUCTURAL): the abstract :meth:`apply` declares
    ``apply(x, /, *, initial_guess=None)`` — the
    :class:`~orpheus.numerics.iteration.SupportsSeededApply` contract
    the iteration drivers consume — so a new sibling CANNOT forget the
    keyword: pyright rejects an override that drops it (LSP), and
    ``ABCMeta`` blocks instantiating a sibling that fails to implement
    it.  Members with no use for a start accept-and-ignore, documented
    per type (an exact inverse has nothing to seed; the sweep threads it
    into the curvilinear Carlson closure; the Green threads it as its
    splitting iteration's start).

    Siblings keep exactly three things of their own: the constructor
    GUARD (what makes their ``inner`` invertible — a value check, a
    type, a derivable splitting), the :meth:`apply` body (the inversion
    algorithm), and ``__repr__``.

    The ADJOINT axis is NOT part of the back-half: ``is_adjointable`` /
    ``.H`` stay at the base defaults — the adjoint-inverse is the #280
    family, deferred (free for the iterative branch, a reverse-DAG
    ``sweep_transpose`` for the direct sweep).

    (This wrap-delegate family is ONE of two kinds of inverse in the
    algebra: ALGEBRA-CLOSED structures invert into themselves — a
    permutation's inverse IS a permutation, a scaled operator's a scaled
    operator — and stay unwrapped as first-class forwards. The canonical
    statement of the split lives on :meth:`OperatorProduct.inverse`.)
    """

    def __init__(self, inner: _ForwardT) -> None:
        #: The forward operator :math:`A` this is the inverse of.
        self.inner = inner

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        # An inverse maps the forward's CODOMAIN back to its DOMAIN.
        return self.inner.codomain

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        return self.inner.domain

    @abstractmethod
    def apply(self, x: Any, /, *, initial_guess: Any | None = None) -> Any:
        r"""Return :math:`A^{-1}\,x` by this sibling's inversion algorithm.

        ``initial_guess`` is the inverse family's canonical driver
        signature: iterative drivers thread the
        previous iterate uniformly, with no per-type signature probes.
        """
        ...

    def solve(self, b: Any, /) -> Any:
        r"""Solve :math:`A^{-1}\,y = b`, i.e. return :math:`A\,b` (the forward).

        The un-invert face: an inverse object IS invertible, and its
        realization verb is the forward matvec — keeping the involution
        web closed (``is_invertible ⟺ a working solve`` on every family
        member).
        """
        return self.inner.apply(b)

    @property
    def is_invertible(self) -> bool:
        return True  # (A^{-1})^{-1} = A — the wrapped forward itself

    def inverse(self) -> _ForwardT:
        r"""Return :math:`(A^{-1})^{-1} = A` — the wrapped forward, by identity.

        The involution law holds as an OBJECT-IDENTITY
        fact: ``A.inverse().inverse() is A``.
        """
        return self.inner


class InverseOperator(InverseWrapMixin[_InvertibleForward], LinearOperator):
    r"""The inverse operator :math:`A^{-1}` of a solve-backed leaf, in operator form.

    The GENERIC member of the #226 inverse family — the name is earned by
    exactly the universal contract and nothing more ("round-trip
    alone earns only *InverseOperator*"): :meth:`apply` inverts, the
    round-trip :math:`A^{-1}(A\,x) = x` holds to the forward's own ``solve``
    precision, and no fancier invariant (S-direct seed-independence,
    G-Neumann, M-materialise) is promised. Structured inverses with a
    distinguishing invariant get their own named types
    (:class:`~orpheus.sn.operators.sweep_operator.SweepOperator` for the
    triangular sweep;
    :class:`~orpheus.numerics.green_operator.GreenOperator` for the
    preconditioned-splitting sum;
    :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
    for the dense direct factorization) — this class serves any
    solve-backed forward with NO more specific named inverse: the
    value-bearing LEAVES (:class:`DiagonalOperator`,
    :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`),
    whose inverse action is an exact pointwise division, AND the
    invertible COMPOSITES: :meth:`OperatorProduct.inverse` returns
    ``InverseOperator(self)`` (#285), so the wrapped inverse action there
    is the product's own ``solve``, :math:`B^{-1}(A^{-1}\,q)`, not a
    division.

    **One realization, not a reciprocal twin.** :meth:`apply` DELEGATES to
    the forward's own :meth:`solve`, bit-identical to today's gated call —
    it does NOT re-derive the inverse action. For a value-bearing LEAF the
    delegation matters doubly: a reciprocal twin would (a) differ from
    ``solve`` by a rounding (:math:`(1/c)\cdot b \neq b/c` in FP), and
    (b) for a cross-section multiplier mint a units-dishonest "reciprocal
    cross-section" field (:math:`1/\Sigma` is a mean free path, a
    DIFFERENT named quantity). The division realization carries the
    inverse semantics without either lie.

    The wrap-delegate back-half (domain↔codomain swap /
    ``solve→inner.apply`` / ``is_invertible`` / ``inverse()→inner``) is
    inherited from :class:`InverseWrapMixin`. This class keeps only its
    ctor guard (the leaf's own ``is_invertible`` value check),
    :meth:`apply`, and ``__repr__``.
    """

    def __init__(self, inner: _InvertibleForward) -> None:
        if not inner.is_invertible:
            raise NotInvertible(
                f"InverseOperator requires an invertible leaf; "
                f"{type(inner).__name__}.is_invertible is False."
            )
        super().__init__(inner)

    def apply(self, x: Any, /, *, initial_guess: Any | None = None) -> Any:
        r"""Return :math:`A^{-1}\,x` — the leaf's own ``solve`` (bit-identical).

        ``initial_guess`` is the inverse family's CANONICAL driver
        signature: iterative drivers thread the previous
        iterate uniformly, with no per-type signature probes.  An EXACT
        pointwise inverse has no use for a starting point — the argument is
        accepted and unused (contrast
        :class:`~orpheus.sn.operators.sweep_operator.SweepOperator`, whose
        sweep threads it into the curvilinear Carlson closure).
        """
        del initial_guess  # exact inverse — no iterative start to seed
        return self.inner.solve(x)

    def __repr__(self) -> str:
        return f"InverseOperator({self.inner!r})"


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

    A permutation is always invertible (:math:`P^{-1} = P^T`), and its
    inverse is ALGEBRA-CLOSED: :meth:`inverse` returns the inverse
    permutation as a first-class :class:`PermutationOperator` (#226
    taxonomy step 1) whose ``apply`` is the same
    ``np.take(·, inverse_perm)`` gather as :meth:`apply_transpose`.

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

    # NO ``solve``: the inverse is ALGEBRA-CLOSED —
    # :meth:`inverse` returns the inverse permutation as a first-class
    # forward whose ``apply`` is the SAME ``np.take(·, inverse_perm)``
    # gather (P^{-1} = P^T, bit-identical) — so solving is
    # ``.inverse().apply(b)``; ``apply_transpose`` keeps the gather as
    # the Euclidean-transpose verb.

    @property
    def is_invertible(self) -> bool:
        return True  # P^{-1} = P^T — a permutation is always invertible

    def inverse(self) -> "PermutationOperator":
        r"""Return :math:`P^{-1}` as a first-class :class:`PermutationOperator`.

        Built on the precomputed :attr:`inverse_perm` — its :meth:`apply`
        is the SAME integer gather as this operator's :meth:`solve` /
        :meth:`apply_transpose` (bit-identical: no arithmetic at all), and
        ``argsort`` of a permutation is exactly involutive in integer math,
        so :math:`(P^{-1})^{-1}` reproduces :attr:`perm` EXACTLY (§13 I2).
        """
        return PermutationOperator(self.inverse_perm, axis=self.axis)

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

    STRUCTURALLY non-invertible — the mask is rank-deficient (it
    projects) — so it declares no ``inverse()``; misuse is a static
    error (Design C). Self-adjoint on the Euclidean axis, so
    ``apply_transpose`` is honest.

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

    The identity body is self-adjoint (``apply_transpose`` echoes) and
    the operator is structurally non-invertible AS POSED here (the wrap
    is a boundary pushforward, not a solvable map) — no ``inverse()``
    declared. The output is a fresh copy of
    the input — matching the legacy
    :class:`~orpheus.geometry.boundary.periodic.PeriodicBoundary.apply`
    aliasing-safety contract (``psi_out.copy()``) and the project-
    wide convention that ``op.apply(psi)`` may be mutated freely by
    the caller without affecting ``psi``.
    """

    def apply(self, x: np.ndarray) -> np.ndarray:
        # Return a fresh copy — the caller-mutates-output
        # safe-aliasing contract.
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
    matter (the operators commute on the joint tensor). Both structural
    axes are factor-wise INTERSECTIONS — invertible iff every factor
    is, adjointable iff every factor is — computed recursively by the
    predicates, and the inverse is ALGEBRA-CLOSED (a tensor product of
    the factor inverses), so there is no ``solve`` verb.

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
        # Eager apply-guard per factor (composition time, never at call).
        for op in ops:
            if not callable(getattr(op, "apply", None)):
                raise TypeError(
                    f"TensorProductOperator factor must expose 'apply'; "
                    f"{type(op).__name__} lacks it."
                )
        self.ops: tuple = tuple(ops)

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
        out = x
        # Adjoint of tensor product is tensor product of adjoints.
        # Apply transposes of factors (order irrelevant for disjoint
        # axes); the per-factor guard-narrow licenses each call.
        for op in self.ops:
            if not adjointable(op):
                raise MissingAdjoint(
                    f"TensorProductOperator.apply_transpose requires "
                    f"every factor to transpose ((A⊗B)^T = A^T⊗B^T); "
                    f"{type(op).__name__}.is_adjointable is False."
                )
            out = op.apply_transpose(out)
        return out

    # NO ``solve``: the inverse is ALGEBRA-CLOSED —
    # :meth:`inverse` returns the tensor product of the factor inverses,
    # a first-class forward — so solving is ``.inverse().apply(b)``.

    @property
    def is_invertible(self) -> bool:
        # (A⊗B)^{-1} = A^{-1}⊗B^{-1} — invertible iff every factor is
        # (recursive over the factors, like every composite predicate).
        return all(op.is_invertible for op in self.ops)

    def inverse(self) -> "TensorProductOperator":
        r"""Return :math:`(A \otimes B \otimes \cdots)^{-1} = A^{-1} \otimes B^{-1} \otimes \cdots`.

        The factor-wise structural inverse (the docstring's "inverse on
        every axis" law). Factor ORDER is preserved —
        the factors act on disjoint axes and commute, exactly as
        :meth:`solve` applies them in stored order — so the action is
        bit-identical to :meth:`solve` given each factor's own
        ``inverse().apply ≡ solve`` identity.
        """
        factor_inverses = []
        for op in self.ops:
            if not invertible(op):
                raise NotInvertible(
                    f"TensorProductOperator.inverse requires every factor "
                    f"to be invertible ((A⊗B)^{{-1}} = A^{{-1}}⊗B^{{-1}}); "
                    f"{type(op).__name__}.is_invertible is False."
                )
            factor_inverses.append(op.inverse())
        return TensorProductOperator(tuple(factor_inverses))

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

    def apply(self, x: np.ndarray) -> np.ndarray:
        out = self.summands[0].apply(x)
        for s in self.summands[1:]:
            out = out + s.apply(x)
        return out

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        for s in self.summands:
            if not adjointable(s):
                raise MissingAdjoint(
                    f"SumOfTensorProductsOperator.apply_transpose requires "
                    f"every summand to transpose; "
                    f"{type(s).__name__}.is_adjointable is False."
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
    Invertible iff every coefficient entry is non-zero — the
    VALUE-dependent arm of the split: :meth:`inverse` is declared and
    :meth:`solve` divides, both refusing eagerly
    (:class:`NotInvertible`) on a zero entry.

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
        # The value-dependent guard (Pattern 4 heritage: never a silent
        # IEEE NaN on a σ=0 division — the legacy bare-σ path had no gate).
        if not self.is_invertible:
            raise NotInvertible(
                "DiagonalOperator.solve requires non-zero coefficient "
                "entries; this operator has at least one zero entry."
            )
        b_arr = np.asarray(b_vec)
        self._check_shape(b_arr)
        return b_arr / self._broadcast(b_arr.ndim)

    @property
    def is_invertible(self) -> bool:
        # Invertible iff every coefficient entry is non-zero (D^{-1} = 1/c).
        return bool(np.all(self.coefficient != 0.0))

    def inverse(self) -> "InverseOperator":
        r"""Return :math:`D^{-1}` as an :class:`InverseOperator` over this leaf.

        Delegation, NOT a reciprocal-coefficient twin: the returned
        object's ``apply`` IS :meth:`solve` (the division ``b / c``),
        bit-identical — whereas ``DiagonalOperator(1/c)`` would multiply
        by a rounded reciprocal and drift by an ulp. The generic name is
        the honest one (round-trip alone earns exactly
        "InverseOperator"; a diagonal division carries no distinguishing
        invariant beyond it).
        """
        if not self.is_invertible:
            raise NotInvertible(
                "DiagonalOperator.inverse requires non-zero coefficient "
                "entries; this operator has at least one zero entry."
            )
        return InverseOperator(self)

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

    Adjointable exactly when the row is an
    :class:`~orpheus.numerics.functional.InnerProductFunctional` (the usual
    case, including its
    :class:`~orpheus.transport.reaction_rate_functional.ReactionRateFunctional`
    specialisation) — the VALUE-dependent arm on the adjoint axis. Rank-1
    operators are **structurally non-invertible** (no ``inverse()`` declared —
    the kernel is the orthogonal complement of the row along the contracted axis),
    but they DO have a **transpose**: :meth:`apply_transpose` is the dual dyad
    :math:`|w\rangle\langle v|` — swap the column :math:`v` with the row's
    weight :math:`w`, contracting :math:`\langle v,\cdot\rangle` over the same
    axis. This is the Euclidean transpose :math:`A^{T}`; the metric-correct
    Hilbert adjoint :math:`A^\dagger = G^{-1}A^{T}G` is the
    :attr:`~LinearOperator.H` wrapper's job. The fission adjoint
    :math:`F^\dagger\psi^* = \nu\Sigma_f\,(\chi\cdot\psi^*)` is exactly this
    dyad-swap (#276). A nonlinear / opaque functional has no dual
    column, so such a dyad advertises ``apply`` only. See
    :ref:`operator-adjoint`.

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

    @property
    def is_adjointable(self) -> bool:
        # The VALUE-dependent arm on the ADJOINT axis: the dual dyad
        # |w⟩⟨v| exists iff the row ⟨w| is a genuine co-vector whose
        # weight IS the dual column — an InnerProductFunctional (the
        # ReactionRateFunctional specialisation included). A nonlinear /
        # opaque functional has no dual column. Mirrors — and gates —
        # the apply_transpose realization below.
        from orpheus.numerics.functional import InnerProductFunctional

        return isinstance(self.functional, InnerProductFunctional)

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

        # The isinstance IS the narrowing (the same fact is_adjointable
        # advertises): the body reads the IPF-typed row's .weight/.axis.
        if not isinstance(self.functional, InnerProductFunctional):
            raise MissingAdjoint(
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

