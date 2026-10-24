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
  zero-coefficient multiplier reports ``False``; a sum of two or more
  terms reports ``False``, its inverse being named by its split),
  recursive on composites;
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
from collections.abc import Hashable, Iterable
from enum import Enum
from functools import cached_property, reduce
from typing import (
    final,
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
D2 = TypeVar("D2", bound=Vector)  # __matmul__ other-operand domain
_Word = TypeVar("_Word")  # a term of the left fold (an operator, narrowed or not)

__all__ = [
    "LinearOperator",
    "ZeroMorphism",
    "ZeroOperator",
    "SupportsInverse",
    "SupportsAdjoint",
    "SupportsTransposeSolve",
    "BlockRole",
    "BulkOperator",
    "FullOperator",
    "BoundaryOperator",
    "NotInvertible",
    "MissingAdjoint",
    "invertible",
    "adjointable",
    "transpose_invertible",
    "posed",
    "IncompatibleOperatorComposition",
    "MatrixTooLarge",
    "InverseWrapMixin",
    "OperatorSum",
    "OperatorProduct",
    "linear_combination",
    "composition",
    "FusedAtom",
    "AdjointOperator",
    "RieszLowerOperator",
    "RieszRaiseOperator",
    "IdentityOperator",
    "PointwiseOperator",
    "InverseMetricOperator",
    "TensorProductOperator",
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
      The typed realized operators of SN
      (:class:`~orpheus.sn.boundary.realized_operators.RealizedBoundary`) carry
      this role on their leaves. ``B`` becomes a first-class
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


def _agreed_space(
    ops: tuple, role: str, owner: str,
) -> Optional["FunctionSpace"]:
    r"""The space a COMMUTATIVE composite's operands agree on, or ``None``.

    The two composites whose operands commute — the sum :math:`A + B` and the
    tensor product :math:`A \otimes B` on disjoint axes — cannot resolve their
    spaces BY POSITION, because position is not part of their identity: if
    ``(A & B).domain`` were ``ops[0].domain``, the derived space of an
    order-INDEPENDENT operator would be order-DEPENDENT, contradicting the
    type's own defining law. So they resolve by **agreement** instead, and this
    is that law, written once:

    * every operand that declares the space must declare the SAME one, else
      the composite is ill-posed and this raises;
    * an operand that declares nothing contributes nothing — silence is not
      disagreement (the module-wide ``None`` semantics: spaces are optional
      while the tree migrates, and a skipped check is never a failed one);
    * all silent ⟹ ``None``, the composite is unbound like its operands.

    Contrast :class:`OperatorProduct`, the one NON-commutative composite:
    ``A @ B`` genuinely maps ``B.domain → A.codomain``, so it resolves by
    position and an unbound factor there poisons that end. The law follows the
    algebra in both cases.

    ⚠ **Why agreement is enough, and where it would stop being enough.** A
    factor's binding in this module is a WHOLE-space binding, not a per-leg
    one: a :class:`~orpheus.numerics.elementary_operators.PermutationOperator` with ``axis=0`` acting on a
    ``(4, 3)`` trace declares ``domain.shape == (4, 3)`` — both axes — because
    it broadcasts on the rest. So the factors of a tensor product are not
    describing separate legs to be multiplied together; each bound factor is
    describing the WHOLE space, and at most one of them can be non-trivial.
    Should genuine per-leg bindings ever arrive (an energy-dependent group
    kernel bound on its own axis), agreement becomes the wrong law and a
    product-space constructor is what has to be built — which is what the
    docstrings on the two tensor-product classes say, rather than silently
    picking one leg.

    ⚠ **The message keeps the phrase** ``equal <role>s`` **deliberately.** It
    predates this helper (it was :class:`OperatorSum`'s own inline wording) and
    two gates in ``tests/gates/sn/operators/test_typed_residual_evaluation.py`` pin it
    as the provenance marker that says *this* guard fired and not some
    incidental raise elsewhere. ``owner`` is prefixed so the marker still
    identifies WHICH composite refused now that the law is shared.
    """
    declared = [
        (op, space)
        for op in ops
        if (space := getattr(op, role, None)) is not None
    ]
    if not declared:
        return None
    first_op, agreed = declared[0]
    for op, space in declared[1:]:
        if space != agreed:
            raise IncompatibleOperatorComposition(
                f"{owner} requires equal {role}s; "
                f"{type(first_op).__name__} declares {agreed!r} while "
                f"{type(op).__name__} declares {space!r}."
            )
    return agreed


class _Identity:
    r"""An operator compared by OBJECT IDENTITY, never by its own ``__eq__``.

    The free vector space is the free vector space on operator OBJECTS (ruling
    23, design D2): two separately built ``L`` stay two terms, because whether
    two different objects are the same operator is a question the algebra
    cannot answer. Keying on this wrapper is what keeps an atom's ``__eq__``
    out of the algebra — so an atom with ``__hash__ = None`` (the mutable
    dataclass ``StreamingOperator``) sums, compares and collects like any
    other.
    """

    __slots__ = ("operator",)

    def __init__(self, operator: object) -> None:
        self.operator = operator

    def __eq__(self, other: object) -> bool:
        return isinstance(other, _Identity) and other.operator is self.operator

    def __hash__(self) -> int:
        return id(self.operator)


def _algebra_key(op: object) -> Hashable:
    r"""The key by which the algebra collects and compares ``op``.

    An atom is keyed by its identity (:class:`_Identity`). A formal arrow (the
    adjoint, the dual, the inverse) is keyed by its kind and its inner's key,
    so ``A.H == A.H`` across two builds although each ``.H`` mints a fresh
    wrapper (ruling 53's spec, Q5; memoisation was rejected). A combination
    or a word is keyed by its structure. The keyed classes are exactly the
    :class:`_KeyedEquality` ones, each declaring its key; every other operator,
    nominal or structural, is an atom.
    """
    return op._algebra_key() if isinstance(op, _KeyedEquality) else _Identity(op)


class _KeyedEquality:
    r"""Equality and hashing by :func:`_algebra_key` — the formal objects.

    Mixed into the classes whose identity is their structure: the
    combination, the word and the three formal arrows. An atom inherits
    ``object``'s identity equality, which is the same relation the key reads.
    """

    __slots__ = ()

    @abstractmethod
    def _algebra_key(self) -> Hashable:
        """The structural key: the kind, and the keys of the parts."""

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, _KeyedEquality):
            return NotImplemented
        return _algebra_key(self) == _algebra_key(other)

    def __hash__(self) -> int:
        return hash(_algebra_key(self))


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
    :data:`~orpheus.sn.boundary.realizer.SN_REALIZERS`
    (vacuum / reflective / white / albedo / periodic / prescribed inflow)
    and
    :data:`~orpheus.diffusion.boundary_realizer.DIFFUSION_REALIZERS`
    (the albedo family incl. zero flux, #290) — carry
    :attr:`BlockRole.BOUNDARY` via
    :func:`~orpheus.transport.boundary_realization.stamp_boundary_role`.

    **EVERY realizable law is stamped, prescribed inflow included** (P3,
    2026-08-05). Until P3 this docstring named it as the one exception,
    on the reading that the rank-0 affine source is the boundary *source*
    ``q.boundary`` rather than a linear boundary operator ``B``. The
    affine split is real and unchanged — see
    :ref:`bc-affine-source-channel` — but it does not put the law outside
    this marker: the law is ``γ₋ψ = L γ₊ψ + q``, the realizer realizes
    ``L``, and for prescribed inflow ``L = 0``. A zero morphism is an
    ordinary linear boundary operator, and it is the same one vacuum
    realizes to.

    ⛔ The exception was also not doing the job it appeared to do. An
    unstamped leaf is NOT excluded from ``B``:
    :attr:`~orpheus.sn.operators.boundary.SNBoundaryOperator._face_operators`
    collects every face's operator with no ``block_role`` filter, so the
    pre-P3 AFFINE operator reached the block regardless — measured
    ``|B(0)| = q`` and ``|B(2x) − 2B(x)| = q``, and on the Krylov path a
    raised ``ConvergenceClaimError``, because an affine map breaks
    the Arnoldi relation GMRES's residual depends on. The stamp is
    honest metadata about a leaf's role; it is not, and never was, the
    fence.
    """

    _role = BlockRole.BOUNDARY


class NotInvertible(TypeError):
    r"""Asked for the inverse of an operator that cannot produce one.

    The INVERSE-axis refusal: raised **eagerly** by
    :meth:`inverse` overrides (and the inverse-family constructors) when
    the operator's :attr:`~LinearOperator.is_invertible` is ``False`` —
    the VALUE-dependent arm of the two-kinds split. A zero-coefficient
    :class:`~orpheus.numerics.elementary_operators.DiagonalOperator`, a sum of two or more terms
    (whose inverse is named by its split), a product with a singular factor: the TYPE supports
    inversion, this INSTANCE refuses, at construction of the inverse and
    never mid-iteration. (The STRUCTURAL arm — :class:`~orpheus.numerics.operator.ZeroOperator`,
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


def _check_dimension(
    op: "LinearOperator", shape: tuple[int, ...], max_dimension: int,
) -> None:
    r"""Refuse a densification whose basis dimension exceeds ``max_dimension``.

    The one size gate of every densification of ``op``, whichever algorithm
    computes it (:meth:`LinearOperator.as_matrix` by probing,
    :func:`~orpheus.numerics.assembled_operator.dense` by assembly).
    """
    n = int(np.prod(shape)) if shape else 1
    if n > max_dimension:
        raise MatrixTooLarge(
            f"densifying {type(op).__name__}: basis dimension "
            f"{n} (= prod{shape}) exceeds max_dimension="
            f"{max_dimension}. A dense matrix would commit ~"
            f"{8 * n * n / 1e6:.0f} MB. Raise max_dimension= if this "
            f"size is intended, or keep the operator matrix-free."
        )


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
            f"e.g. (ng, 1) for an infinite-medium group operator)."
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

    The two composites are wired through the dunders: ``+``, ``-``, unary
    ``-``, ``α *`` and ``/ α`` mint the formal linear combination
    (:func:`linear_combination`, an :class:`OperatorSum`), and ``@`` the
    flattened word (:func:`composition`, an :class:`OperatorProduct`), so the
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
    #: (:class:`OperatorSum` / :class:`AdjointOperator`) and the
    #: :func:`~orpheus.transport.boundary_realization.stamp_boundary_role` stamp assign
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
    #: is — the passthrough (:class:`AdjointOperator`) and the
    #: :func:`_join_system_roles` union over the terms (:class:`OperatorSum`).
    system_role: Optional[SystemRole] = None

    @property
    @abstractmethod
    def domain(self) -> Optional["FunctionSpace"]:
        """The function space this operator consumes — DEMANDED of every
        subclass (the S4-amendment, 2026-08-22: an operator is not an
        operator without its two spaces, and this base no longer supplies
        a silent ``None`` default an implementer can inherit unawares).

        Every class must ANSWER, in one of four ways:

        * **bind** — return a space stored/threaded at construction (the
          :class:`~orpheus.transport.operators.fission.FissionOperator`
          precedent: space MANDATORY, non-Optional return);
        * **derive** — compute it from operands/held data (the
          composites' agreement/position laws, :class:`InverseWrapMixin`'s
          swap, a mesh-holding leaf reading its mesh's space);
        * **the pointwise law** — :class:`PointwiseOperator` members
          answer ``None`` BY LAW (space-polymorphic: the domain is the
          operand's space at operation time, discriminated by type);
        * **a documented Optional** — an explicit override returning
          ``FunctionSpace | None`` whose docstring names the campaign
          that owns its mandatory flip. ⛔ This arm is EMPTY since CS4c
          step 6 item 6.4 (2026-09-07): S/C/iso flipped at CS4c K2b /
          steps 2–5 (``BoundOperator``'s mandatory kw-only ends), L and
          the three boundary leaves at item 6.4 — the last documented
          Optionals in the tree. It stays listed so a future leaf that
          reaches for it knows the answer is "bind or derive", never a
          new Optional.

        The post-flip law: a BOUND leaf's ends are never ``None`` — only
        the pointwise-law members answer ``None``, BY LAW. What was never
        legal is SILENCE.
        """
        ...

    @property
    @abstractmethod
    def codomain(self) -> Optional["FunctionSpace"]:
        """The function space this operator produces — DEMANDED of every
        subclass. See :attr:`domain` for the four legitimate answers and
        the Optional-until-terminal semantics.
        """
        ...

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
        one-term combination :math:`cA` as invertible exactly when :math:`A`
        is, a sum of two or more terms as non-invertible, and a
        zero-coefficient
        :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
        as singular (``min|f| = 0``). Composites derive it recursively from
        their operands; the default is ``False`` — an operator is
        invertible only by explicit override.
        """
        return False

    @property
    def is_adjointable(self) -> bool:
        r"""Whether this operator realises its Euclidean transpose (``apply_transpose``, :meth:`dual`).

        The RUNTIME predicate for the transpose axis: ``True`` iff
        ``apply_transpose`` runs on every input of this instance, so
        :meth:`dual` constructs and applies. The Hilbert adjoint
        :attr:`H` = ``♯ ∘ dual ∘ ♭`` needs a second fact, that both ends'
        metrics give Riesz legs (:func:`posed`): ``.H`` exists iff both hold,
        and refuses with :class:`MissingAdjoint` otherwise (`[M]`
        2026-10-18, ``tests/gates/numerics/test_adjointability_is_honest.py``,
        which pins both directions). The
        transpose-of-a-sum law :math:`(A+B)^{\mathsf T} = A^{\mathsf T} +
        B^{\mathsf T}` is realised in the composer method bodies; this
        predicate is the matching *advertisement* —
        ``(A+B).is_adjointable == A.is_adjointable and B.is_adjointable`` —
        structurally computed rather than cached in a string set. Default
        ``False``; an operator with a working ``apply_transpose`` overrides.
        """
        return False

    @property
    def is_transpose_invertible(self) -> bool:
        r"""Whether this forward realises its TRANSPOSED solve :math:`(A^{\mathsf T})^{-1}`.

        The third capability of a forward, beside :attr:`is_invertible`
        (``solve``) and :attr:`is_adjointable` (``apply_transpose``):
        ``True`` iff ``solve_transpose`` exists AND runs on every input of
        this instance. It is the capability the inverse's adjoint needs, by
        the law :math:`(A^{-1})^{\mathsf T} = (A^{\mathsf T})^{-1}`: an
        inverse object is adjointable iff its forward is transpose-invertible
        (:class:`InverseOperator`, the SN sweep, the coupled substitution all
        read this one predicate). Neither of the other two flags answers it:
        a multi-D ``(L+C)`` transposes its action but not its sweep, and a
        ``callable(solve_transpose)`` test passes a method that raises (`[M]`
        2026-10-18, ``tests/gates/numerics/test_adjointability_is_honest.py``).
        Default ``False``; a forward with a working ``solve_transpose``
        overrides. Three composites derive it: :class:`OperatorProduct`
        through each factor's inverse OBJECT (whether it transposes, so a
        scaled or permuted factor qualifies though it has no transposed
        solve of its own), :class:`~orpheus.numerics.coupled_system.CoupledOperator`
        from its diagonal blocks' own transposed solves (the substitution
        calls ``block.solve_transpose``), and the inverse family
        (:class:`InverseWrapMixin`) from its forward's transpose. A sum has
        none: its inverse is a Green, whose transpose is the Green of the
        dual sum. A one-term combination :math:`cA` has none either; a native
        one is owed when a scaled block first sits on a coupled grid's
        diagonal.
        """
        return False

    @property
    def is_metric_free_adjoint(self) -> bool:
        r"""Whether this operator's Hilbert adjoint needs NO metric — i.e.
        it equals its own Euclidean transpose under EVERY inner product.

        ``True`` on the constant multipliers :math:`cI`
        (:class:`ConstantMultiplier`), which commute with every metric, and
        on the bound classes whose adjoint is a metric-free law
        (:class:`ZeroMorphism`'s swapped zero); composites derive it
        recursively (the module's predicate discipline — computed from
        operands, never cached). A field multiplier is NOT metric-free: it
        commutes only with diagonal metrics (issue #625). The adjoint door
        reads this: an UNBOUND operator may take ``.H`` iff this is
        ``True``, because then and only then is the transpose the honest
        answer rather than the R2 hazard. Default ``False``.
        """
        return False

    @property
    def is_assemblable(self) -> bool:
        r"""Whether this operator emits its own sparse matrix.

        The ASSEMBLY axis, the third structural surface beside
        :attr:`is_invertible` and :attr:`is_adjointable`: ``True`` iff a
        structural ``(row, col, value)`` emission of this operator exists
        (the stencil-assembly third consumption mode of the per-cell
        closure algebra). Composites derive it from their parts: a sum or a
        word emits iff every term or factor does. The default is ``False``:
        an operator emits only by explicit override, and the probing
        :meth:`as_matrix` stays the total representation for everything
        else. The emission itself and its laws live downstream
        (:func:`~orpheus.numerics.assembled_operator.assemble`, ruling 24);
        this module states only the flag.
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

    def __add__(self, other: "LinearOperator[Domain, Codomain]") -> "LinearOperator[Domain, Codomain]":
        if _is_scalar(other):
            return NotImplemented
        return linear_combination(((1.0, self), (1.0, other)))

    def __radd__(self, other: "LinearOperator[Domain, Codomain]") -> "LinearOperator[Domain, Codomain]":
        # The scalar 0 is the zero of ``sum(operators)``'s start value: absorbed.
        if _is_scalar(other):
            return self if other == 0 else NotImplemented
        return linear_combination(((1.0, other), (1.0, self)))

    def __sub__(self, other: "LinearOperator[Domain, Codomain]") -> "LinearOperator[Domain, Codomain]":
        if _is_scalar(other):
            return NotImplemented
        return linear_combination(((1.0, self), (-1.0, other)))

    def __rsub__(self, other: "LinearOperator[Domain, Codomain]") -> "LinearOperator[Domain, Codomain]":
        if _is_scalar(other):
            return NotImplemented
        return linear_combination(((1.0, other), (-1.0, self)))

    def __mul__(self, other: float) -> "LinearOperator[Domain, Codomain]":
        if not _is_scalar(other):
            return NotImplemented
        return linear_combination(((float(other), self),))

    def __rmul__(self, other: float) -> "LinearOperator[Domain, Codomain]":
        return self.__mul__(other)

    def __neg__(self) -> "LinearOperator[Domain, Codomain]":
        r"""Unary minus: :math:`-A`, the one-term combination :math:`(-1)A`."""
        return linear_combination(((-1.0, self),))

    def __truediv__(self, scalar: float) -> "LinearOperator[Domain, Codomain]":
        r"""Scalar division: ``A / α`` is ``(1/α) * A``.

        Reads more naturally than the reciprocal-multiply form
        ``(1.0 / α) * A`` — eigenvalue/Krylov normalisation
        (``F / k_eff``), homogenisation averages. Division by zero raises
        :class:`ZeroDivisionError` per the standard Python convention.
        """
        if not _is_scalar(scalar):
            return NotImplemented
        return linear_combination(((1.0 / float(scalar), self),))

    def __matmul__(self, other: "LinearOperator[D2, Domain]") -> "LinearOperator[D2, Codomain]":
        # ``self`` (Domain → Codomain) ∘ ``other`` (D2 → Domain) is D2 → Codomain.
        return composition(self, other)

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
            **Eagerly, here at construction** when this operator is not
            :attr:`is_adjointable`, or when its adjoint is not
            :func:`posed` (an unbound end; the
            refusal is :class:`AdjointOperator`'s) — a wrapper that
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
        return AdjointOperator(self)

    @property
    def H(self) -> "LinearOperator[Codomain, Domain]":
        """Alias for :meth:`adjoint` — the Hilbert-adjoint vocabulary
        (``A.H`` reads as "A dagger"). Swaps the carriers
        ``[Domain, Codomain] → [Codomain, Domain]`` (see
        :meth:`adjoint`)."""
        return self.adjoint()

    def dual(self) -> "LinearOperator[Codomain, Domain]":
        r"""Return the dual arrow :math:`A^{\mathsf T} : W^* \to V^*`.

        The METRIC-FREE half of the Hilbert adjoint: the representation
        transpose, carried between the dual spaces. For ``A : V → W`` the
        dual maps ``W* → V*`` — the same array arithmetic as
        :meth:`apply_transpose`, with the arrow bookkeeping made explicit
        (:meth:`~orpheus.numerics.space.FunctionSpace.dual` on each end).
        The Hilbert adjoint is then the theorem-shaped composition

        .. math::

           A^{*} \;=\; \sharp_V \circ A^{\mathsf T} \circ \flat_W
           \qquad\text{(``RieszRaiseOperator(domain) ∘ A.dual() ∘
           RieszLowerOperator(codomain)``)}

        realized by :class:`AdjointOperator` — the dual is the middle
        factor, and the metrics live entirely in the Riesz legs
        (:class:`RieszLowerOperator` / :class:`RieszRaiseOperator`).

        Involution: ``A.dual().dual() is A`` (object identity — the
        wrapper's own :meth:`_DualOperator.dual` returns the inner).

        Raises
        ------
        MissingAdjoint
            Eagerly, when this operator is not :attr:`is_adjointable`
            (no working ``apply_transpose``) — same broken-stub-refusing
            style as :meth:`adjoint`.
        """
        if not adjointable(self):
            raise MissingAdjoint(
                f"{type(self).__name__} is not adjointable — .dual() "
                f"requires is_adjointable=True (a working apply_transpose); "
                f"the dual arrow IS the transpose, carried between the "
                f"dual spaces."
            )
        return _DualOperator(self)

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

        This default is the apply-to-basis pattern. An operator that HOLDS
        its matrix reads it instead
        (:class:`~orpheus.numerics.assembled_operator.SparseAssembledOperator`
        densifies its CSR;
        :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
        runs one batched LU backsolve). An operator that only KNOWS its
        stencil does not override: its structural densification is
        :func:`~orpheus.numerics.assembled_operator.dense`, chosen by the
        consumer. The return is a DENSE :class:`numpy.ndarray`.

        Parameters
        ----------
        basis_shape : tuple[int, ...], optional
            The element shape :meth:`apply` consumes. Default: derived
            from :attr:`domain` (``domain.shape``); REQUIRED explicitly
            for operators carrying no space (bare/test-constructed
            operators on the legal-until-CS4 ``None`` path — the
            infinite-medium production operators thread the
            problem's pose since campaign 1 CS1 and derive).
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

        **The representation, by probing** (numerics_untangling.md ruling
        60): the matrix is fixed by the spaces' coordinates, the axes in
        their order, and column :math:`j` is :math:`A e_j` raveled. This is
        total over operators with spaces. Structural assembly is a faster
        algorithm for the same matrix and lives downstream
        (:func:`~orpheus.numerics.assembled_operator.dense`); the consumer
        that pays for the densification chooses it. The probed and the
        assembled matrices are two algorithms the equivalence gates compare:
        they share the leaves' coefficients, and differ in everything that
        composes them (the laws against :math:`n` applications).
        """
        shape = _resolve_basis_shape(self, basis_shape)
        _check_dimension(self, shape, max_dimension)
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


class SupportsTransposeSolve(LinearOperator[Domain, Codomain], Protocol):
    r"""Narrowing target: an operator whose ``solve_transpose`` may be called.

    The transposed-solve sibling of :class:`SupportsInverse` /
    :class:`SupportsAdjoint`: the branch narrowed by
    :func:`transpose_invertible` may call ``solve_transpose``, which returns
    :math:`(A^{\mathsf T})^{-1} y`, the realisation the inverse family's
    transpose rides by the law :math:`(A^{-1})^{\mathsf T} = (A^{\mathsf T})^{-1}`.
    Same discipline: never ``isinstance`` this, never annotate a parameter
    with it to demand the capability.
    """

    def solve_transpose(self, y: Domain, /) -> Codomain: ...


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


def transpose_invertible(
    op: "LinearOperator[Domain, Codomain]",
) -> "TypeGuard[SupportsTransposeSolve[Domain, Codomain]]":
    r"""Checked bridge: narrow ``op`` to :class:`SupportsTransposeSolve` iff it realises its transposed solve.

    The third bridge beside :func:`invertible` and :func:`adjointable`, with
    their one-directional ``TypeGuard`` semantics: ``True`` licenses
    ``op.solve_transpose`` (:attr:`~LinearOperator.is_transpose_invertible`).
    The attribute is read with a ``False`` default because an operand typed
    only by a narrower protocol need not carry it.
    """
    return bool(getattr(op, "is_transpose_invertible", False))


def posed(op: "LinearOperator[Domain, Codomain]") -> bool:
    r"""Whether ``op``'s Hilbert adjoint :math:`A^{*} = \sharp_V A^{\mathsf T} \flat_W` is POSED.

    The free-function reader beside :func:`adjointable` (which answers only
    whether a Euclidean transpose is coded, so ``.dual()`` runs): ``.H``
    exists iff both hold. The adjoint is defined by the two inner products,
    so it is posed exactly when the Riesz legs of both ends exist: both ends
    are declared spaces, primal or dual alike (a dual space carries the
    inverse Riesz form, numerics_untangling.md ruling 26), or the operator's
    adjoint needs no metric (:attr:`~LinearOperator.is_metric_free_adjoint`).
    The adjoint door (:class:`AdjointOperator`) reads this one predicate.
    """
    if getattr(op, "is_metric_free_adjoint", False):
        return True
    ends = (getattr(op, "domain", None), getattr(op, "codomain", None))
    return all(end is not None for end in ends)


def _inverse_transposes(op: "LinearOperator[Domain, Codomain]") -> bool:
    r"""Whether :math:`(A^{-1})^{\mathsf T}` exists through the canonical surface: ``op`` is
    invertible and its inverse OBJECT transposes.

    The predicate of the swap laws :math:`(A^{*})^{-1} = (A^{-1})^{*}`
    (:attr:`AdjointOperator.is_invertible`) and
    :math:`(A^{\mathsf T})^{-1} = (A^{-1})^{\mathsf T}`
    (:attr:`_DualOperator.is_invertible`), and of a product's transposed
    solve: one spelling. A leaf's NATIVE transposed solve is
    :func:`transpose_invertible`, which the inverse object's own flag reads
    (:class:`InverseWrapMixin`), so the two agree on every wrapped leaf.
    """
    return invertible(op) and adjointable(op.inverse())


# ───────────────────────────────────────────────────────────────────────
# Composition primitives
# ───────────────────────────────────────────────────────────────────────


# ---------------------------------------------------------------------------
# Adjoint wrapper
# ---------------------------------------------------------------------------


class RieszLowerOperator(LinearOperator[Domain, Domain], Generic[Domain]):
    r"""The Riesz LOWERING leg :math:`\flat : V \to V^*` — apply the metric.

    ``♭ x = G x``: the isomorphism a Hilbert space's inner product induces
    between the space and its dual (the Riesz representation theorem,
    realized). Delegates to
    :meth:`~orpheus.numerics.space.FunctionSpace.apply_metric`, so the
    metric ARITHMETIC stays single-sourced in the space's resolved
    :class:`~orpheus.numerics.metric.HilbertMetric` (or its per-axis
    path) — this class contributes only the ARROW: ``domain`` is the
    primal ``V``, ``codomain`` is ``V.dual()``, both non-Optional by
    construction (the leg is born bound).

    Constructed by :class:`~orpheus.numerics.operator.RieszLowerOperator`
    and by :class:`AdjointOperator` (whose codomain-side factor this is).

    **On a dual space** the leg is the inverse map: a
    :class:`~orpheus.numerics.space.DualSpace` carries the inverse Riesz form
    :math:`G^{+}` (numerics_untangling.md ruling 26, issue #617), so
    ``RieszLowerOperator(V*)`` : :math:`V^* \to V^{**} = V` applies
    :math:`G^{+}`, and ``lower_{V*}(lower_V(x))`` is :math:`x` on the
    metric's range. Until 2026-10-21 a dual space carried :math:`G` itself,
    so this composition applied :math:`G^2` (`[M]` 2026-08-30:
    ``[0.25, 4, 16]`` for ``w = [0.5, 2, 4]``) and the constructor refused
    dual spaces.

    ``apply_transpose`` is ``apply``: every shipped metric realization is
    symmetric (a diagonal weight; a :class:`DenseMetric` admitted only
    through its symmetry guard), so ``♭ᵀ = ♭`` under the reflexive
    identification ``V** = V``.
    """

    def __init__(self, space: "FunctionSpace") -> None:
        self.space: Final = space

    @property
    def domain(self) -> "FunctionSpace":
        return self.space

    @property
    def codomain(self) -> "FunctionSpace":
        return self.space.dual()

    def apply(self, x: Domain) -> Domain:
        return self.space.apply_metric(x)

    def apply_transpose(self, x: Domain) -> Domain:
        # ♭ᵀ = ♭ (symmetric metric, V** = V) — see the class docstring.
        return self.space.apply_metric(x)

    @property
    def is_adjointable(self) -> bool:
        # The transpose is coded (``apply_transpose``: a symmetric G is its
        # own transpose), so ``.dual()`` runs, and ``.H`` is posed on both
        # ends (the codomain V* carries G⁺, ruling 26): the law ♭ᴴ = ♯.
        return True


class RieszRaiseOperator(LinearOperator[Domain, Domain], Generic[Domain]):
    r"""The Riesz RAISING leg :math:`\sharp : V^* \to V` — apply the
    pseudo-inverse metric.

    ``♯ f = G⁺ f``: the inverse Riesz map, Moore–Penrose everywhere by
    the metric family's doctrine (the reciprocal on the metric's range,
    zero on its kernel — the tangential ``|Ω·n| = 0`` trace slots).
    Delegates to
    :meth:`~orpheus.numerics.space.FunctionSpace.apply_inverse_metric`;
    this class contributes only the arrow: ``domain`` is ``V.dual()``,
    ``codomain`` the primal ``V``.

    The round trip is the honest law of the pseudo-inverse, NOT the
    identity: ``♯ ∘ ♭ = P_range(G)`` — the identity exactly when the
    metric is strictly positive, the tangential-zeroing projector on a
    singular trace block (`[M]` 2026-08-30: a legal 2-D
    ``Quadrature.product(4,4)`` mesh has 32/64 tangential trace slots and
    round-trip trace defect ``2.87``; every strictly-positive fixture
    reads ``≤ 4.4e-16``).

    The same symmetric-transpose identity and the same construction sites as
    :class:`RieszLowerOperator` (its exact mirror — the domain-side factor of
    :class:`AdjointOperator`); on a dual space it applies :math:`G`.
    """

    def __init__(self, space: "FunctionSpace") -> None:
        self.space: Final = space

    @property
    def domain(self) -> "FunctionSpace":
        return self.space.dual()

    @property
    def codomain(self) -> "FunctionSpace":
        return self.space

    def apply(self, x: Domain) -> Domain:
        return self.space.apply_inverse_metric(x)

    def apply_transpose(self, x: Domain) -> Domain:
        # ♯ᵀ = ♯ (symmetric pseudo-inverse of a symmetric metric).
        return self.space.apply_inverse_metric(x)

    @property
    def is_adjointable(self) -> bool:
        # As for ♭ (:class:`RieszLowerOperator`): the transpose is coded,
        # and ``.H`` is posed on both ends (ruling 26): ♯ᴴ = ♭.
        return True


class _DualOperator(_KeyedEquality, LinearOperator[Codomain, Domain], Generic[Domain, Codomain]):
    r"""The dual arrow :math:`A^{\mathsf T} : W^* \to V^*` — the
    representation transpose with honest dual-space bookkeeping.

    The metric-free middle factor of the Hilbert adjoint (see
    :meth:`LinearOperator.dual`). ``apply`` IS ``inner.apply_transpose``;
    what this wrapper adds is the ARROW: domain/codomain are the DUALS of
    the inner's codomain/domain, so a composition chain can track bras vs
    kets through :class:`~orpheus.numerics.space.DualSpace` (its minted
    consumer — the §1 non-endomorphism made physical).

    Same explicit ``Generic[Domain, Codomain]`` pinning as
    :class:`AdjointOperator` (the PEP-696 parameter-order note there).

    Laws, as structure:

    * involution — ``A.dual().dual() is A`` (:meth:`dual` returns the
      inner: object identity);
    * transpose of the dual is the inner's action —
      :meth:`apply_transpose` delegates to ``inner.apply``;
    * dagger–dual commutation — ``(A.dual()).H = (A.H).dual()``, a theorem
      once each dual space carries the inverse Riesz form (ruling 26): the
      dual arrow's adjoint is the base sandwich on its dual ends, and the
      commutation is gated as a law, not taken as a second route.
    """

    def __init__(self, inner: "SupportsAdjoint[Domain, Codomain]") -> None:
        self.inner: Final = inner

    def _algebra_key(self) -> Hashable:
        # A formal arrow is its kind and its inner (ruling 53's spec, Q5).
        return ("dual", _algebra_key(self.inner))

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        c = getattr(self.inner, "codomain", None)
        return c.dual() if c is not None else None

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        d = getattr(self.inner, "domain", None)
        return d.dual() if d is not None else None

    def apply(self, y: Codomain) -> Domain:
        return self.inner.apply_transpose(y)

    def apply_transpose(self, x: Domain) -> Codomain:
        # (Aᵀ)ᵀ = A — the inner's forward action, no metric anywhere.
        return self.inner.apply(x)

    @property
    def is_adjointable(self) -> bool:
        return True

    def dual(self) -> "LinearOperator[Domain, Codomain]":
        # Involution as OBJECT IDENTITY: (Aᵀ)ᵀ = A.
        return self.inner  # type: ignore[return-value]

    @property
    def is_invertible(self) -> bool:
        r"""The swap law on the dual arrow: :math:`(A^{\mathsf T})^{-1} = (A^{-1})^{\mathsf T}`.

        The Euclidean twin of :attr:`AdjointOperator.is_invertible`: the dual
        is invertible iff the inner is AND the inner's inverse transposes.
        """
        return _inverse_transposes(self.inner)

    def inverse(self) -> "LinearOperator[Domain, Codomain]":
        r"""The inverse of the dual = the dual of the inverse, as an object.

        :math:`(A^{\mathsf T})^{-1} = (A^{-1})^{\mathsf T}`, routed as
        ``A.inverse().dual()`` (the twin of :meth:`AdjointOperator.inverse`).
        Gated by :attr:`is_invertible`.
        """
        if not (invertible(self.inner) and adjointable(self.inner.inverse())):
            raise NotInvertible(
                f"_DualOperator.inverse(): the swap law (A^T)^{{-1}} = "
                f"(A^{{-1}})^T needs the inner {type(self.inner).__name__} "
                f"invertible and its inverse adjointable."
            )
        return self.inner.inverse().dual()

class AdjointOperator(_KeyedEquality, LinearOperator[Codomain, Domain], Generic[Domain, Codomain]):
    r"""The Hilbert-adjoint ARROW :math:`A^{*} : W \to V` — first-class.

    Public since CS4c step 1 (the dagger-arrow promotion, R2 ruling —
    until then the private ``_AdjointOperator``): the realization of the
    dagger functor, constructed by :meth:`LinearOperator.adjoint` /
    ``A.H`` (the canonical door — direct construction is legal for an
    :func:`adjointable`-narrowed inner and repeats the same eager gates).

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
    :class:`SupportsAdjoint` — there is no lazy capability gate left to
    fail at call time.

    **The arrow is a composition of three first-class factors, built at
    construction** (CS4c R2 ruling, 2026-08-30): the domain/codomain
    exchange is realized in code — for ``A : V → W``,

    .. math::

       A^{*} \;=\; \underbrace{\sharp_V}_{V^* \to V}
       \circ \underbrace{A^{\mathsf T}}_{W^* \to V^*}
       \circ \underbrace{\flat_W}_{W \to W^*}

    i.e. ``RieszRaiseOperator(inner.domain) ∘ inner.dual() ∘
    RieszLowerOperator(inner.codomain)`` (the formula reserved at
    ``metric.py``'s CS4c-compatibility note). The metric arithmetic
    lives in the legs (single-sourced through the space's
    :class:`~orpheus.numerics.metric.HilbertMetric`); this class holds
    the composition plus what only the named arrow can carry: the
    eager refusals, the #280 swap law (:meth:`inverse` — an object
    identity), the role passthrough, and the dagger laws as structure:

    * involution — ``A.H.H is A`` (:meth:`adjoint` returns the inner);
    * transpose of the adjoint — a THEOREM of the legs,
      :math:`(A^{*})^{\mathsf T} = \flat_W \circ A \circ \sharp_V`
      (metrics symmetric by admission), which closed #375's dead-end:
      ``A.H`` is adjointable, and ``A.H.H`` is reachable — and IS ``A``.

    An UNBOUND end (a ``None`` space on a metric-free-exempt inner)
    contributes an :class:`IdentityOperator` leg — the Euclidean metric
    as the neutral element, bit-identical to the pre-CS4c skip branch.
    """

    def __init__(self, inner: "SupportsAdjoint[Domain, Codomain]") -> None:
        # The S4-amendment's unbound-.H refusal (user ruling, 2026-08-22):
        # the Hilbert adjoint is defined BY the two inner products, so an
        # operator that has not declared its spaces has no Hilbert adjoint
        # to take — and the pre-amendment behaviour (apply the Euclidean
        # transpose and skip the metric sandwiches) was the catalogued R2
        # hazard: "a bare Euclidean transpose wearing the Hilbert
        # adjoint's name". Eagerly, here at construction, in this class's
        # own broken-stub-refusing style. The metric-free families never
        # reach this constructor: a PointwiseOperator's adjoint() returns
        # itself, ZeroMorphism's the swapped zero map.
        # The door reads ONE predicate, :func:`posed`.
        if not posed(inner):
            raise MissingAdjoint(
                f"{type(inner).__name__} is UNBOUND (domain/codomain "
                f"None) — the Hilbert adjoint needs the two spaces' "
                f"metrics; declare both, or use apply_transpose for the "
                f"bare representation transpose. (A constant multiplier, "
                f"ConstantMultiplier, is metric-free; a field multiplier "
                f"needs the space whose metric it meets.)"
            )
        self.inner = inner
        # The G-adjoint transposes the 2×2 block matrix (A_bs ↔ A_sb^T),
        # which preserves WHICH blocks are touched — so the role is the
        # inner operator's role: ``L.H`` is FULL, ``B.H`` is BOUNDARY,
        # ``C.H`` is BULK.
        self.block_role = getattr(inner, "block_role", None)
        # The G-adjoint also preserves which SYSTEMS are touched: A_AB.H
        # (bulk→ray) is still COUPLED, A_BB.H still System B.
        self.system_role = getattr(inner, "system_role", None)
        # The three factors, built HERE — the domain/codomain exchange
        # realized at construction (R2 ruling): the codomain-side ♭ and
        # the domain-side ♯ are the swapped ends' Riesz legs; the dual
        # arrow is the metric-free middle. A None end (metric-free
        # exemption above) gets the Euclidean neutral element — the same
        # arithmetic the pre-leg skip branch performed, bit-identically.
        inner_domain = getattr(inner, "domain", None)
        inner_codomain = getattr(inner, "codomain", None)
        self._lower: "LinearOperator" = (
            RieszLowerOperator(inner_codomain)
            if inner_codomain is not None
            else IdentityOperator()
        )
        self._dual: "_DualOperator[Domain, Codomain]" = _DualOperator(inner)
        self._raise: "LinearOperator" = (
            RieszRaiseOperator(inner_domain)
            if inner_domain is not None
            else IdentityOperator()
        )

    def _algebra_key(self) -> Hashable:
        # A formal arrow is its kind and its inner (ruling 53's spec, Q5).
        return ("adjoint", _algebra_key(self.inner))

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        # Adjoint of A: V → W is A.H: W → V — domain swaps with inner.codomain.
        return getattr(self.inner, "codomain", None)

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        return getattr(self.inner, "domain", None)

    def apply(self, y: Codomain) -> Domain:
        # The Hilbert-adjoint action, as the leg composition built at
        # construction:
        #   (A^* y)_V = ♯_V(Aᵀ(♭_W y)) = G_V⁺ ⊙ apply_transpose(G_W ⊙ y)
        # Same call order, same delegation targets as the pre-leg inline
        # spelling (each leg delegates to the space's apply_metric /
        # apply_inverse_metric, so the SAME composition serves a
        # flat-ndarray metric AND a composite bulk ⊕ trace metric with
        # its pseudo-inverse on the singular partial-current trace) —
        # bit-identical by construction, gated by G-A1. The legs are the
        # individually-mutable seams the ledger's per-leg battery reads.
        return self._raise.apply(self._dual.apply(self._lower.apply(y)))

    def apply_transpose(self, x: Domain) -> Codomain:
        r"""The representation transpose of the adjoint — a THEOREM of
        the legs (#375's four-line composition, landed).

        :math:`(A^{*})^{\mathsf T} = (G_V^{+} A^{\mathsf T}
        G_W)^{\mathsf T} = G_W\, A\, G_V^{+} = \flat_W \circ A \circ
        \sharp_V` — the metrics transpose to themselves because every
        shipped realization is symmetric (diagonal weights; a
        :class:`~orpheus.numerics.metric.DenseMetric` is admitted only
        through its symmetry guard). Until CS4c this raised a stub with
        zero witnesses; the capability replaces it, and
        :attr:`is_adjointable` advertises it honestly.
        """
        return self._lower.apply(self.inner.apply(self._raise.apply(x)))

    @property
    def is_adjointable(self) -> bool:
        # True by the theorem above: apply_transpose needs only the legs
        # (always built) and the inner's forward action (always present).
        return True

    def adjoint(self) -> "LinearOperator[Domain, Codomain]":
        r"""The dagger involution as an OBJECT IDENTITY: ``A.H.H is A``.

        :math:`(A^{*})^{*} = A` — no double wrapper, no arithmetic: the
        inner IS the answer. (`[M]` pre-CS4c, ``A.H.H`` was unreachable —
        #375's headline: ``is_adjointable`` read ``False`` and the
        transpose raised.)
        """
        return self.inner  # type: ignore[return-value]

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
        return _inverse_transposes(self.inner)

    def inverse(self) -> "LinearOperator[Domain, Codomain]":
        r"""The inverse of the adjoint = the adjoint of the inverse (#280).

        :math:`(A^{*})^{-1} = (A^{-1})^{*}` — the operator-algebra swap law,
        an OBJECT IDENTITY here (not a computed numerical equivalence): this
        wrapper IS ``A.H``, so its inverse routes to ``A.inverse().H``. Gated
        by :attr:`is_invertible` (``A`` invertible and ``A.inverse()``
        adjointable). The metric adjoint-solve
        :math:`A^{-1\,*} b = G_V^{+}\,\mathrm{apply\_transpose}(G_W\,b)` then
        falls out of :meth:`apply` (which already routes
        ``inner.apply_transpose``) FOR FREE — no ``AdjointOperator.solve`` /
        no metric code enters the sweep.
        """
        if not invertible(self.inner):
            raise NotInvertible(
                f"AdjointOperator.inverse(): the inner "
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
                f"AdjointOperator.inverse(): the inner's inverse "
                f"{type(inner_inverse).__name__} is not adjointable, so "
                f"(A.inverse()).H does not exist — the swap law needs an "
                f"adjointable inverse (e.g. an SN SweepOperator, #280 2.5c). "
                f"is_invertible is False."
            )
        return inner_inverse.H


# ───────────────────────────────────────────────────────────────────────
# The free algebra: a sum is a formal linear combination, a product a word
# ───────────────────────────────────────────────────────────────────────
#
# The operators V → W form a real vector space, and under composition a
# monoid. The two composites realise exactly those structures and nothing
# more (numerics_untangling.md ruling 23, package S2):
#
#   OperatorSum     = {word ↦ coefficient}   unordered; empty = the zero map
#   OperatorProduct = (f₁, …, fₙ), n ≥ 2     ordered; no factor a product
#   a word          = an atom, or a product
#
# The laws of the free vector space (associativity, commutativity, collection
# of like terms, scalar distributivity, 0·A = 0, A − A = 0) and of the free
# monoid (associativity) are EQUALITIES OF THE REPRESENTATION: the two
# constructors below normalise, so ``(A + B) + C`` and ``A + (B + C)`` are the
# same combination and ``(AB)C`` and ``A(BC)`` the same word. Composition is
# NOT distributed over a sum — ``A(B + C)`` stays a word whose second factor
# is a sum — because distributing changes the cost of ``apply`` (two
# applications of ``A``), and ruling 23 claims only the abelian-group laws and
# scalar distributivity.


def _is_scalar(value: object) -> bool:
    """A real scalar the algebra multiplies by (``bool`` excluded)."""
    return isinstance(value, (int, float, np.floating, np.integer)) and not isinstance(value, bool)


def _capability_refusal(
    subject: str, requirement: str, law: str, parts: "Iterable[object]", predicate: str,
) -> str:
    r"""The one refusal message of a composite capability law: the law, and each part's predicate.

    A composite has a capability exactly when every term or factor has it;
    the refusal names the law that needs it and reads the predicate off each
    part, so the part that lacks it is in the message.
    """
    readings = ", ".join(f"{type(part).__name__}.{predicate}={getattr(part, predicate)}" for part in parts)
    return f"{subject} requires {requirement} ({law}); {readings}."


def _require_apply(op: object, owner: str) -> None:
    # Eager, at composition: a term that cannot act is refused HERE, never at
    # the first ``apply`` inside a Krylov iteration.
    if not callable(getattr(op, "apply", None)):
        raise TypeError(
            f"{owner} requires apply on every operand; {type(op).__name__} "
            f"lacks 'apply'."
        )


def _ends(op: object) -> "tuple[Optional[FunctionSpace], Optional[FunctionSpace]]":
    """The (domain, codomain) an operand declares; ``None`` where it declares none."""
    return getattr(op, "domain", None), getattr(op, "codomain", None)


def linear_combination(terms: "Iterable[tuple[float, LinearOperator]]") -> "LinearOperator":
    r"""Return :math:`\sum_i c_i A_i`, normalised — the one constructor of a sum.

    Every ``+``, ``-``, unary ``-``, ``α *`` and ``/ α`` routes here. The
    result is normalised, so each law of the free vector space is an
    equality of what this returns:

    * an operand that is itself an :class:`OperatorSum` is FLATTENED, its
      coefficients scaled (associativity; ``α(A + B) = αA + αB``);
    * terms with the same :func:`_algebra_key` COLLECT, their coefficients
      added (``A + A = 2A``; ``(α + β)A = αA + βA``), in first-insertion
      order;
    * a zero coefficient drops its term, and a zero operand
      (:class:`ZeroMorphism`, :class:`ZeroOperator`) is ABSORBED; each still
      contributes its spaces;
    * no surviving term is the zero map: :class:`ZeroMorphism` on the agreed
      spaces, or the pointwise :class:`ZeroOperator` when every operand was
      pointwise (spaceless by law); with an unbound non-pointwise operand
      the zero cannot be minted and the construction refuses;
    * one surviving term ``1·A`` whose own ends are the combination's is
      ``A`` itself (``(A + B) − B is A``; ``1·A is A``).

    Every operand that ENTERS contributes its spaces to the agreement law
    (:func:`_agreed_space`), whether it survives as a term, cancels, is an
    absorbed zero, or is a flattened inner sum: ``A − A`` must still know it
    is the zero map ``V → W``, or it could not mint its :class:`ZeroMorphism`.

    The coefficients are floats, so the collection is exact only where the
    float arithmetic is: ``0.1A + 0.2A − 0.3A`` leaves ``5.55e-17·A``, not the
    zero. That is recorded (``test_linear_combination.py`` c05), never hidden
    behind a tolerance in ``__eq__``.
    """
    words: list[tuple[float, LinearOperator]] = []
    carriers: list[object] = []
    entered: list[object] = []
    for coefficient, op in terms:
        _require_apply(op, "OperatorSum")
        c = float(coefficient)
        carriers.append(op)
        if isinstance(op, OperatorSum):
            entered.extend(w for _, w in op.terms)
            if c != 0.0:
                words.extend((c * ci, wi) for ci, wi in op.terms)
            continue
        entered.append(op)
        if c != 0.0 and not isinstance(op, (ZeroMorphism, ZeroOperator)):
            words.append((c, op))
    collected: dict[Hashable, list] = {}
    for c, w in words:
        slot = collected.get(key := _algebra_key(w))
        if slot is None:
            collected[key] = [c, w]
        else:
            slot[0] += c
    surviving = tuple((c, w) for c, w in collected.values() if c != 0.0)
    agreed_domain = _agreed_space(tuple(carriers), "domain", "OperatorSum")
    agreed_codomain = _agreed_space(tuple(carriers), "codomain", "OperatorSum")
    if not surviving:
        return _zero_of(entered, agreed_domain, agreed_codomain)
    if len(surviving) == 1:
        c, w = surviving[0]
        if c == 1.0 and _ends(w) == (agreed_domain, agreed_codomain):
            return w
    return OperatorSum(surviving, domain=agreed_domain, codomain=agreed_codomain)


def _zero_of(
    entered: list[object],
    domain: Optional["FunctionSpace"],
    codomain: Optional["FunctionSpace"],
) -> "LinearOperator":
    r"""The zero map of a combination whose every term cancelled or was zero.

    :class:`ZeroMorphism` on the agreed spaces; the pointwise
    :class:`ZeroOperator` when every operand was pointwise (an endomorphism at
    every space, by the family law, so it needs none); otherwise a refusal.

    SCOPE-BOUNDARY[guard] machinery: every operator end a space, so that a
    non-pointwise operator always brings both (numerics_untangling.md ruling
    26, case (1)). ruling: the S2 gate spec's Q3, the orchestrator, 2026-10-20
    (the plan's calls under ruling 51). revisit: when every end is bound, this
    refusal cannot fire.
    """
    if domain is not None and codomain is not None:
        zero: LinearOperator = ZeroMorphism(domain=domain, codomain=codomain)
    elif all(isinstance(op, PointwiseOperator) for op in entered):
        zero = ZeroOperator()
    else:
        raise ValueError(
            f"OperatorSum: the zero combination of "
            f"{', '.join(sorted({type(op).__name__ for op in entered}))} has no "
            f"bound spaces, so its ZeroMorphism cannot be minted. Bind the "
            f"operator's domain and codomain (a declared scope edge: every "
            f"operator end a space is ruling 26, not yet built)."
        )
    # The zero of a combination touches the blocks and systems its cancelled
    # terms touched: ``0·B`` of a boundary law is still a boundary-block map
    # (the role lives on the leaves and every composite derives it, the
    # user's ruling of 2026-10-20).
    if entered:
        zero.block_role = reduce(_join_block_roles, (getattr(op, "block_role", None) for op in entered))
        zero.system_role = reduce(_join_system_roles, (getattr(op, "system_role", None) for op in entered))
    return zero


class OperatorSum(_KeyedEquality, LinearOperator[Domain, Codomain]):
    r"""A formal linear combination :math:`\sum_i c_i\,w_i` of operator words.

    The element of the free vector space on operator words, :math:`V \to W`:
    :attr:`terms` holds the ``(coefficient, word)`` pairs, each word an atom
    or an :class:`OperatorProduct`, never a sum, never zero, no two with the
    same :func:`_algebra_key`. Equality and the hash read the MULTISET of
    terms and the ends (``A + B == B + A``); the insertion order is kept only
    for the floating-point fold of :meth:`apply`.

    Minted by :func:`linear_combination` — every ``+``, ``-``, ``α *`` routes
    there — which normalises; constructing the class directly is reserved to
    that function (its terms must already be collected).

    Every capability is the linear law over the terms, read once:

    * :meth:`apply` :math:`= \sum_i c_i\,w_i(x)`, a LEFT FOLD in insertion
      order, so a left-nested spelling ``((A − S) − N) − B`` applies as
      ``((Ax − Sx) − Nx) − Bx``, bit for bit what the binary sum computed.
      The fold multiplies only where :math:`c \ne \pm 1` (a −1 is a
      subtraction, bit-identical to adding :math:`(-1)\,w(x)`), and it
      accumulates in place into a buffer it allocated itself, never into an
      operand's output (``test_composition_cost.py`` holds the allocation at
      ≤ 4 vectors per degree of freedom);
    * :meth:`apply_transpose`, :meth:`dual`, :attr:`is_assemblable`,
      :attr:`is_adjointable`, :attr:`block_role`, :attr:`system_role`: the
      fold of the terms' own;
    * :attr:`is_invertible`: only a ONE-term combination with an invertible
      word, :math:`(cA)^{-1} = c^{-1}A^{-1}`. A sum of two or more terms has
      no inverse of its own (ruling 50): its inverse is named by its
      splitting, :class:`~orpheus.numerics.green_operator.GreenOperator`
      ``(M, N)`` for :math:`A = M - N`.

    The Hilbert adjoint is the base's :math:`\sharp \circ \text{dual} \circ
    \flat`, taken ONCE at the combination's ends, with the dual distributed
    over the terms (ruling 53's precise form). The literal distribution
    :math:`\sum_i c_i\,w_i^{*}` is a theorem only where every term is bound
    to a non-degenerate metric: a spaceless pointwise term on a dense metric
    misses reciprocity by 1.26e-3 under it (#625), which the once-at-the-ends
    form does not.
    """

    def __init__(
        self,
        terms: "tuple[tuple[float, LinearOperator], ...]",
        *,
        domain: Optional["FunctionSpace"],
        codomain: Optional["FunctionSpace"],
    ) -> None:
        self._terms: Final = terms
        self._domain: Final = domain
        self._codomain: Final = codomain
        # Roles DERIVED from the terms: the sum touches the union of the blocks
        # (systems) its terms touch, so ``(L + C − S − F − B)`` carries FULL by
        # construction, whatever order it was spelled in.
        self.block_role = reduce(
            _join_block_roles, (getattr(w, "block_role", None) for _, w in terms),
        )
        self.system_role = reduce(
            _join_system_roles, (getattr(w, "system_role", None) for _, w in terms),
        )

    @property
    def terms(self) -> "tuple[tuple[float, LinearOperator], ...]":
        """The ``(coefficient, word)`` pairs, in first-insertion order."""
        return self._terms

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        return self._domain

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        return self._codomain

    def _algebra_key(self) -> Hashable:
        return (
            "+",
            frozenset((_algebra_key(w), c) for c, w in self._terms),
            self._domain,
            self._codomain,
        )

    def apply(self, x: Domain, /, *args: Any, **kwargs: Any) -> Codomain:
        if len(self._terms) == 1:
            return self._apply_scaled(x, *args, **kwargs)
        if args or kwargs:
            raise TypeError(
                f"OperatorSum.apply: a combination of {len(self._terms)} terms "
                f"takes one operand; got the extra arguments "
                f"{args!r} {sorted(kwargs)!r}."
            )
        return _left_fold(self._terms, lambda w: w.apply(x))

    def _apply_scaled(self, x: Domain, /, *args: Any, **kwargs: Any) -> Codomain:
        r"""The one-term action :math:`c\,w(x)`, its extra arguments forwarded.

        A seed ``initial_guess`` lives in the OUTPUT space of :math:`c\,w`, so
        the word's own solve is seeded with it divided by :math:`c`:
        :math:`(\alpha A)^{-1} = \alpha^{-1}A^{-1}` seeds :math:`A^{-1}` with
        :math:`\alpha x_0` (the S2 spec's Q4; its witness is c30).
        """
        c, w = self._terms[0]
        seed = kwargs.get("initial_guess")
        if seed is not None:
            kwargs = {**kwargs, "initial_guess": seed / c}
        return c * w.apply(x, *args, **kwargs)

    def _adjointable_terms(self, verb: str) -> "tuple[tuple[float, SupportsAdjoint], ...]":
        # The checked bridge, once per call: (Σ cᵢwᵢ)ᵀ = Σ cᵢwᵢᵀ needs every term's.
        terms = tuple((c, w) for c, w in self._terms if adjointable(w))
        if len(terms) != len(self._terms):
            raise MissingAdjoint(_capability_refusal(
                f"OperatorSum.{verb}", "every term to transpose", "(Σ cᵢwᵢ)ᵀ = Σ cᵢwᵢᵀ",
                (w for _, w in self._terms), "is_adjointable",
            ))
        return terms

    def apply_transpose(self, y: Codomain, /, *args: Any, **kwargs: Any) -> Domain:
        terms = self._adjointable_terms("apply_transpose")
        if len(terms) == 1:
            c, w = terms[0]
            seed = kwargs.get("initial_guess")
            if seed is not None:  # the seed lives in the output of cAᵀ, as in apply
                kwargs = {**kwargs, "initial_guess": seed / c}
            return c * w.apply_transpose(y, *args, **kwargs)
        if args or kwargs:
            raise TypeError(
                "OperatorSum.apply_transpose: a combination of two or more "
                "terms takes one operand."
            )
        return _left_fold(terms, lambda w: w.apply_transpose(y))

    @property
    def is_adjointable(self) -> bool:
        return all(w.is_adjointable for _, w in self._terms)

    @property
    def is_metric_free_adjoint(self) -> bool:
        # A combination of metric-free operators is metric-free (derived).
        return all(w.is_metric_free_adjoint for _, w in self._terms)

    def dual(self) -> "LinearOperator[Codomain, Domain]":
        r"""Return :math:`(\sum_i c_i w_i)^{\mathsf T} = \sum_i c_i w_i^{\mathsf T}`.

        The dual distributes over the terms exactly, with no metric anywhere
        (ruling 53). The combination of the duals carries the dual spaces of
        this one's ends through the absorbed zero map between them, so a sum
        whose every term is spaceless keeps the spaces an absorbed zero gave
        it.
        """
        dual_terms: list[tuple[float, LinearOperator]] = [
            (c, w.dual()) for c, w in self._adjointable_terms("dual")
        ]
        if self._domain is not None and self._codomain is not None:
            dual_terms.append(
                (1.0, ZeroMorphism(domain=self._codomain.dual(), codomain=self._domain.dual()))
            )
        return linear_combination(dual_terms)

    @property
    def is_invertible(self) -> bool:
        r"""``True`` iff this is ONE term :math:`cA` with :math:`A` invertible.

        :math:`(cA)^{-1} = c^{-1}A^{-1}` is algebra-closed. A sum of two or
        more terms has no inverse of its own (ruling 50): there is no
        algorithm for :math:`(A + B)^{-1}` from the terms, and the
        preconditioned splitting that does invert it needs a SPLIT the sum
        does not carry — which term is implicit is the Strategy's choice,
        never the spelling's.
        """
        return len(self._terms) == 1 and self._terms[0][1].is_invertible

    def inverse(self) -> "LinearOperator[Codomain, Domain]":
        r"""Return :math:`(cA)^{-1} = c^{-1}A^{-1}`, or refuse a true sum."""
        if len(self._terms) != 1:
            raise NotInvertible(
                f"OperatorSum.inverse: a sum of {len(self._terms)} terms has no "
                f"inverse of its own. Name its splitting A = M − N and build "
                f"GreenOperator(M, N) (numerics_untangling.md ruling 50); a "
                f"Strategy mints it from its splitting, which in S_N is "
                f"Splitting.green (ruling 65)."
            )
        c, w = self._terms[0]
        if not invertible(w):
            raise NotInvertible(
                f"OperatorSum.inverse: (cA)^{{-1}} = c^{{-1}}A^{{-1}} needs an "
                f"invertible A; {type(w).__name__}.is_invertible is False."
            )
        return linear_combination(((1.0 / c, w.inverse()),))

    @property
    def is_assemblable(self) -> bool:
        return all(w.is_assemblable for _, w in self._terms)

    def __repr__(self) -> str:
        return "(" + " + ".join(f"{c!r}·{w!r}" for c, w in self._terms) + ")"


def _left_fold(
    terms: "tuple[tuple[float, _Word], ...]",
    act: "Callable[[_Word], Any]",
) -> Any:
    r"""The left fold :math:`((c_1 a_1 + c_2 a_2) + c_3 a_3) + \cdots`, :math:`a_i = \mathrm{act}(w_i)`.

    The rounding of every multi-term action is fixed here, once: the left
    association in insertion order; after the first term a coefficient
    :math:`\pm 1` is applied as an addition or a subtraction (no multiply, no
    negated temporary), any other coefficient multiplied in; the first term
    is multiplied by its coefficient unless that is 1. The accumulator is mutated only once it
    is a value the fold's own arithmetic produced, never an operand's output
    (which may be the operand's input, or a buffer the operand keeps), so the
    remaining terms accumulate into it in place. In place is taken only
    where it is the out-of-place result exactly (:func:`_accumulates_in_place`),
    so the fold's VALUE never depends on the optimisation.
    """
    iterator = iter(terms)
    c, w = next(iterator)
    accumulator = act(w)
    owned = c != 1.0
    if owned:
        accumulator = c * accumulator
    for c, w in iterator:
        value = act(w)
        if owned and _accumulates_in_place(accumulator, value):
            if c == -1.0:
                accumulator -= value
            else:
                accumulator += value if c == 1.0 else c * value
        else:
            accumulator = (
                accumulator - value if c == -1.0
                else accumulator + (value if c == 1.0 else c * value)
            )
            owned = True
        del value
    return accumulator


def _accumulates_in_place(accumulator: Any, value: Any) -> bool:
    r"""Whether ``accumulator += value`` is ``accumulator + value``, exactly.

    NumPy's in-place operators keep the accumulator's shape and dtype, so they
    differ from the out-of-place sum when a later term broadcasts wider or
    promotes the dtype (`[M]` the S2 qa review: a wider broadcast raises, a
    float32 accumulator silently stays float32). A typed carrier with no
    in-place operator falls back to the out-of-place one, which is exact.
    """
    if not (isinstance(accumulator, np.ndarray) and isinstance(value, np.ndarray)):
        return not isinstance(accumulator, np.ndarray)
    return (
        np.result_type(accumulator, value) == accumulator.dtype
        and np.broadcast_shapes(accumulator.shape, value.shape) == accumulator.shape
    )


def composition(*operators: "LinearOperator") -> "LinearOperator":
    r"""Return :math:`A_1 A_2 \cdots A_n`, normalised — the one constructor of a word.

    Every ``@`` routes here. An operand that is an :class:`OperatorProduct`
    is FLATTENED (associativity: ``(AB)C`` and ``A(BC)`` are one word), and
    a one-term combination :math:`cW` gives its coefficient to the result
    and its word to the factors (:math:`A(\alpha B) = \alpha(AB)`) when the
    word carries the combination's ends (a spaceless word whose spaces an
    absorbed zero supplied stays one factor, ends and all). A sum of
    two or more terms is ONE factor: composition is not distributed. An
    :class:`IdentityOperator` factor is kept, because dropping it would
    change what :attr:`OperatorProduct.domain` infers through a pointwise
    neighbour, and the unit law is not claimed.
    """
    coefficient = 1.0
    factors: list[LinearOperator] = []
    for op in operators:
        _require_apply(op, "OperatorProduct")
        # A one-term combination gives up its coefficient only when its word
        # carries the combination's own ends (the collapse test of
        # :func:`linear_combination`); a word that is spaceless where the
        # combination is bound stays inside it, so its ends are not lost.
        if isinstance(op, OperatorSum) and len(op.terms) == 1 and _ends(op.terms[0][1]) == _ends(op):
            c, op = op.terms[0]
            coefficient *= c
        if isinstance(op, OperatorProduct):
            factors.extend(op.factors)
        else:
            factors.append(op)
    word = factors[0] if len(factors) == 1 else OperatorProduct(tuple(factors))
    return word if coefficient == 1.0 else linear_combination(((coefficient, word),))


class OperatorProduct(_KeyedEquality, LinearOperator[Domain, Codomain]):
    r"""A word :math:`A_1 A_2 \cdots A_n` of operators, :math:`n \ge 2`.

    The element of the free monoid: :attr:`factors` holds the factors left
    to right, no factor a product, no factor a one-term combination (its
    coefficient moved to the enclosing combination). Equality and the hash
    read the tuple of the factors' keys. Minted by :func:`composition` —
    every ``@`` routes there.

    The ends: :attr:`domain` is the rightmost factor's that speaks, and
    :attr:`codomain` the leftmost's, a spaceless :class:`PointwiseOperator`
    passing its neighbour's through (it acts at its operand's space). The
    adjacent-ends law is checked at construction through the same
    pass-through, so ``A @ D @ B`` with ``D`` pointwise checks ``A`` against
    ``B``.

    Laws, as folds over the factors:

    * :meth:`apply`: right to left, :math:`A_1(A_2(\cdots A_n x))`;
    * :meth:`solve`: :math:`A_n^{-1} \cdots A_1^{-1}`, through each factor's
      canonical surface ``.inverse().apply``; :meth:`inverse` wraps it as
      :class:`InverseOperator`;
    * :meth:`apply_transpose`: left to right,
      :math:`A_n^{\mathsf T}(\cdots A_1^{\mathsf T} y)`; :meth:`dual` the
      reversed word of the factors' duals;
    * :attr:`is_assemblable`: every factor emits (the law
      :math:`[A_1 \cdots A_n] = [A_1] \cdots [A_n]` lives downstream, in
      :func:`~orpheus.numerics.assembled_operator.assemble`).

    The Hilbert adjoint is the base's, taken ONCE at the word's ends (ruling
    53's precise form). The reversed word of the factors' adjoints
    :math:`A_n^{*} \cdots A_1^{*}` inserts :math:`\flat \sharp = G G^{+}` at
    every middle space, a projector where the middle metric is singular
    (the S\ :sub:`N` partial-current trace until S6): 6.4e-2 off
    reciprocity on such a middle, which the once-at-the-ends form is not.
    """

    def __init__(self, factors: "tuple[LinearOperator, ...]") -> None:
        self._factors: Final = factors
        # The adjacent-ends law, right to left, a spaceless pointwise factor
        # passing its operand's space through.
        passing = getattr(factors[-1], "codomain", None)
        for left in reversed(factors[:-1]):
            if isinstance(left, PointwiseOperator):
                continue
            left_domain = getattr(left, "domain", None)
            if left_domain is not None and passing is not None and left_domain != passing:
                raise IncompatibleOperatorComposition(
                    f"OperatorProduct requires each factor's domain to equal "
                    f"the codomain of the factors to its right; "
                    f"{type(left).__name__} declares domain {left_domain!r} "
                    f"against {passing!r}."
                )
            passing = getattr(left, "codomain", None)

    @property
    def factors(self) -> "tuple[LinearOperator, ...]":
        """The factors, left to right (``A @ B`` is ``(A, B)``)."""
        return self._factors

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        for factor in reversed(self._factors):
            if not isinstance(factor, PointwiseOperator):
                return getattr(factor, "domain", None)
        return None

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        for factor in self._factors:
            if not isinstance(factor, PointwiseOperator):
                return getattr(factor, "codomain", None)
        return None

    def _algebra_key(self) -> Hashable:
        return ("@",) + tuple(_algebra_key(f) for f in self._factors)

    def apply(self, x: Domain, /) -> Codomain:
        value: Any = x
        for factor in reversed(self._factors):
            value = factor.apply(value)
        return value

    def solve(self, b: Codomain, /) -> Domain:
        r"""Solve :math:`(A_1 \cdots A_n)\,x = b` as :math:`A_n^{-1}(\cdots A_1^{-1} b)`."""
        if not self.is_invertible:
            raise NotInvertible(_capability_refusal(
                "OperatorProduct.solve", "every factor to be invertible",
                "(A₁…Aₙ)⁻¹ = Aₙ⁻¹…A₁⁻¹", self._factors, "is_invertible",
            ))
        value: Any = b
        for factor in self._factors:
            value = cast("SupportsInverse", factor).inverse().apply(value)
        return value

    def _adjointable_factors(self, verb: str) -> "tuple[SupportsAdjoint, ...]":
        # The checked bridge, once per call: (A₁…Aₙ)ᵀ needs every factor's.
        factors = tuple(f for f in self._factors if adjointable(f))
        if len(factors) != len(self._factors):
            raise MissingAdjoint(_capability_refusal(
                f"OperatorProduct.{verb}", "every factor to transpose",
                "(A₁…Aₙ)ᵀ = Aₙᵀ…A₁ᵀ", self._factors, "is_adjointable",
            ))
        return factors

    def apply_transpose(self, y: Codomain, /) -> Domain:
        value: Any = y
        for factor in self._adjointable_factors("apply_transpose"):
            value = factor.apply_transpose(value)
        return value

    def solve_transpose(self, y: Domain, /) -> Codomain:
        r"""Solve :math:`(A_1 \cdots A_n)^{\mathsf T} x = y` as :math:`A_1^{-\mathsf T}(\cdots A_n^{-\mathsf T} y)`."""
        if not self.is_transpose_invertible:
            raise NotInvertible(
                "OperatorProduct.solve_transpose requires every factor's "
                "inverse to transpose ((A₁…Aₙ)⁻ᵀ = A₁⁻ᵀ…Aₙ⁻ᵀ); "
                + ", ".join(type(f).__name__ for f in self._factors)
                + "."
            )
        value: Any = y
        for factor in reversed(self._factors):
            inverse = cast("SupportsInverse", factor).inverse()
            if not adjointable(inverse):  # excluded by is_transpose_invertible above
                raise NotInvertible(f"{type(factor).__name__}'s inverse has no transpose.")
            value = inverse.apply_transpose(value)
        return value

    @property
    def is_invertible(self) -> bool:
        return all(f.is_invertible for f in self._factors)

    @property
    def is_transpose_invertible(self) -> bool:
        return all(_inverse_transposes(f) for f in self._factors)

    def inverse(self) -> "InverseOperator":
        r"""Return :math:`(A_1 \cdots A_n)^{-1}`, the generic family member wrapping this word.

        ``inverse().apply(q)`` delegates to :meth:`solve`; the wrapper adds
        the family's seeded-apply contract (#285) and the object-identity
        involution. (The ALGEBRA-CLOSED inverses — a permutation's, the
        identity's, a one-term combination's — are first-class forward
        operators in their own structure and stay unwrapped.)
        """
        if not self.is_invertible:
            raise NotInvertible(_capability_refusal(
                "OperatorProduct.inverse", "every factor to be invertible",
                "(A₁…Aₙ)⁻¹ = Aₙ⁻¹…A₁⁻¹", self._factors, "is_invertible",
            ))
        return InverseOperator(self)

    @property
    def is_adjointable(self) -> bool:
        return all(f.is_adjointable for f in self._factors)

    @property
    def is_metric_free_adjoint(self) -> bool:
        # A word of metric-free operators is metric-free (derived).
        return all(f.is_metric_free_adjoint for f in self._factors)

    def dual(self) -> "LinearOperator[Codomain, Domain]":
        r"""Return :math:`(A_1 \cdots A_n)^{\mathsf T} = A_n^{\mathsf T} \cdots A_1^{\mathsf T}`, exact (ruling 53)."""
        return composition(*(f.dual() for f in reversed(self._adjointable_factors("dual"))))

    @property
    def is_assemblable(self) -> bool:
        return all(f.is_assemblable for f in self._factors)

    def __repr__(self) -> str:
        return "(" + " @ ".join(repr(f) for f in self._factors) + ")"


class FusedAtom(LinearOperator[Domain, Codomain]):
    r"""An atom that realises an algebraic expression with an action of its own.

    The free algebra never fuses (numerics_untangling.md ruling 52): ``+`` and
    ``@`` build combinations and words, never a realisation chosen by a type
    test on the operands. Where a method needs an expression realised by a
    faster or direct action — the WDD sweep for :math:`L + C`, the scheduled
    forward substitution for :math:`L + C - B_{\rm lower}`, the moment-emitting
    sweep for :math:`P\,A^{-1}` — its OWNER builds an atom carrying that action,
    and the atom states the expression it realises as :meth:`expansion`.

    Everything the atom does not realise itself is its expansion's, read from
    ONE cached expansion: the ends, the transpose (through the expansion's
    dual arrow), the adjointability and metric-freeness, assemblability, and the
    block and system roles (:meth:`_adopt_expansion_roles`, called by each
    atom's constructor once its legs are set). A subclass overrides a face
    only where it has a realisation of its own (the sweep atom's
    :math:`M(\sigma)^{\mathsf T}`), so a capability its expansion gains reaches
    the atom with no edit here; the assembly functor reads the expansion too.
    """

    @abstractmethod
    def expansion(self) -> "LinearOperator[Domain, Codomain]":
        """The plain combination or word this atom realises."""

    @cached_property
    def expanded(self) -> "LinearOperator[Domain, Codomain]":
        """The ONE cached :meth:`expansion` every non-realised face reads."""
        return self.expansion()

    def _adopt_expansion_roles(self) -> None:
        # The blocks and systems the realised expression touches.
        self.block_role = getattr(self.expanded, "block_role", None)
        self.system_role = getattr(self.expanded, "system_role", None)

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        return self.expanded.domain

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        return self.expanded.codomain

    def apply_transpose(self, y: Codomain, /) -> Domain:
        # The expansion's dual arrow; its construction refuses (MissingAdjoint)
        # an expression with a part that does not transpose.
        return self.expanded.dual().apply(y)

    @property
    def is_adjointable(self) -> bool:
        return self.expanded.is_adjointable

    @property
    def is_metric_free_adjoint(self) -> bool:
        return self.expanded.is_metric_free_adjoint

    @final
    @property
    def is_assemblable(self) -> bool:
        # Final: an atom realises an ACTION of its own; its assembly is its
        # expansion's, always (ruling 52, read at the S3 review).
        return self.expanded.is_assemblable


class PointwiseOperator(LinearOperator[Domain, Domain]):
    r"""The space-polymorphic POINTWISE (multiplier) family — an
    endomorphism at EVERY admissible space, acting identically.

    The natural stratum of the multiplier algebra the transport layer
    documents at
    :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
    (``M[1] = I``, ``M[0] = 0``, ``M[f]`` — that class is the BOUND,
    typed-carrier realization; this base is the polymorphic engine
    stratum). A member is *multiplication by a coefficient* — the
    constant ``1`` (:class:`IdentityOperator`), the constant ``0``
    (:class:`~orpheus.numerics.operator.ZeroOperator`), a field ``f`` (:class:`~orpheus.numerics.elementary_operators.DiagonalOperator`) —
    so each output entry depends only on the SAME input entry: nothing
    couples, nothing is transposed across.

    The two laws the type declares (the S4-amendment):

    * **No stored space pair.** ``domain``/``codomain`` answer ``None``
      BY LAW — the member acts at the operand's space, obtained at
      operation time, and ``domain == codomain == that space`` for every
      apply. This is a *typed* statement (discriminated by
      ``isinstance``), not the legacy silent default: a bound operator
      answering ``None`` is a declared migration debt; a pointwise
      operator answering ``None`` is permanently, correctly
      space-polymorphic. Composition handles it through the agreement
      law unchanged ("an operand that declares nothing contributes
      nothing" — :func:`_agreed_space`): a pointwise operand CONFORMS to
      its neighbours' space, which is precisely its naturality.
    * **The transpose is the operator:** ``A.dual()`` acts as ``A``
      (each output entry depends only on the same input entry), so every
      member is adjointable in the Euclidean sense.
    * **The Hilbert adjoint needs a metric only for a FIELD.** A constant
      multiplier :math:`cI` commutes with every metric
      (:math:`G\,(cI) = (cI)\,G`), so its Hilbert adjoint is itself on any
      space, with no spaces needed: :class:`ConstantMultiplier`
      (:class:`IdentityOperator`, :class:`ZeroOperator`) is metric-free and
      ``A.H is A``. A field multiplier
      (:class:`~orpheus.numerics.elementary_operators.DiagonalOperator`)
      commutes only with DIAGONAL metrics, so on a dense-metric space its
      Hilbert adjoint is :math:`G^{+} D G \ne D` (`[M]` issue #625: a
      1.4e-3 reciprocity miss on a 3×3 dense metric); a spaceless member has
      no metric to take it against, so its ``.H`` refuses, and it meets a
      metric only inside a composite whose bound ends pose the adjoint
      (numerics_untangling.md ruling 53). (Real coefficients throughout —
      the complexification thread revisits at Campaign 2's resolvent.)

    The boundary of the family is the commutation theorem, not
    endomorphy: :class:`~orpheus.numerics.elementary_operators.PermutationOperator` is endomorphic and acts on
    many spaces, yet permuting does NOT commute with a general diagonal
    metric, so its adjoint genuinely needs the spaces — it is correctly
    NOT a member. Bound multipliers (``MultiplicationOperator``,
    :class:`InverseMetricOperator`) are pointwise in ACTION but defined
    by a space's data — they stay bound classes.
    """

    @property
    def domain(self) -> None:
        r"""``None`` BY LAW — see the class docstring (space-polymorphic:
        the domain is the operand's space, at operation time)."""
        return None

    @property
    def codomain(self) -> None:
        r"""``None`` BY LAW — an endomorphism at the operand's space."""
        return None

    @property
    def is_adjointable(self) -> bool:
        return True  # the Euclidean transpose of a real multiplier is itself


class ConstantMultiplier(PointwiseOperator[Domain]):
    r"""Multiplication by a constant :math:`c`, the scalar multiples of the
    identity: the pointwise members whose Hilbert adjoint needs no metric.

    :math:`cI` commutes with every metric, :math:`G\,(cI) = (cI)\,G`, so the
    sandwich :math:`\sharp\,(cI)^{\mathsf T}\,\flat = G^{+}cG = cI` on the
    metric's range, on any space, and the adjoint is the operator itself:
    ``A.H is A``, algebra-closed, with no spaces needed (issue #625: the
    law holds for constants, not for a field).
    """

    @property
    def is_metric_free_adjoint(self) -> bool:
        return True  # cI commutes with every metric

    def adjoint(self) -> "ConstantMultiplier[Domain]":
        r"""``(cI)^* = cI`` — this very instance."""
        return self


class IdentityOperator(ConstantMultiplier[Domain]):
    r"""The identity operator :math:`I\,x = x` — multiplication by ``1``,
    the pointwise family's unit (:class:`PointwiseOperator`).

    Both axes hold trivially — :math:`I^{-1} = I` and :math:`I^* = I` —
    and both are ALGEBRA-CLOSED: :meth:`inverse` and the inherited
    :meth:`~ConstantMultiplier.adjoint` each return this very
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




class ZeroOperator(ConstantMultiplier[Domain]):
    r"""The ENDOMORPHIC zero :math:`0\,x = 0` — multiplication by ``0``,
    the pointwise family's absorbing element (:class:`PointwiseOperator`).

    Stateless and space-polymorphic: the action routes through ``0.0 * x``,
    echoing the operand's type — the zero OF the operand's own space,
    which is correct precisely because a pointwise member is an
    endomorphism there (``domain == codomain == the operand's space``,
    by the family law). ``apply_transpose`` is the same map
    (:math:`0^{\mathsf T} = 0`), and the inherited
    :meth:`~ConstantMultiplier.adjoint` returns ``self``.

    STRUCTURALLY non-invertible — the singular map par excellence — so it
    declares no ``inverse()`` at all: misuse is a static error, the
    honest surface for a type whose inverse does not exist mathematically
    (Design C; a raising stub would be the harmful-stub anti-pattern this
    module is designed against).

    The zero MAP between two DIFFERENT spaces is a different object —
    :class:`ZeroMorphism`, born bound to its pair. Until the
    S4-amendment (2026-08-22) this class straddled both roles through
    per-site ``codomain_zero``/``transpose_zero`` closures and optional
    space fields (#330); the un-weld split them: the closures' two
    production consumers dissolved (the fission pencil member became
    stack ∘ restriction at A2; the vacuum law binds a
    :class:`ZeroMorphism`), and the hooks retired with them.
    """

    def apply(self, x: Domain, /) -> Domain:
        # Endomorphic by the family law: the zero of the codomain IS the
        # zero of the operand's own space — ``0.0 * x`` echoes it (bare
        # ``np.ndarray`` → ``np.zeros_like(x)`` bit-exact; a typed
        # carrier → a fresh same-class zero via its scalar dunder).
        return cast("Domain", 0.0 * x)

    def apply_transpose(self, x: Domain, /) -> Domain:
        # 0^T = 0 — the same endomorphic echo.
        return cast("Domain", 0.0 * x)

    # is_invertible inherits the base ``False`` — the zero map is singular.


class ZeroMorphism(LinearOperator):
    r"""The zero MAP between two DECLARED spaces:
    :math:`0 : \mathcal D \to \mathcal C`.

    An operator :math:`A : \mathcal D \to \mathcal C` maps the domain
    to the codomain; its action — including the zero action — produces an
    element of :math:`\mathcal C`. The zero map is the ONE operator whose
    action cannot reveal which spaces it connects (every probe returns
    zero), so unlike every other operator it cannot even be *probed* into
    a pair — the pair must be DECLARED, and this class demands both at
    construction (the S4-amendment's base demand, in its sharpest
    instance). The canonical production member is the vacuum boundary
    law's :math:`R = 0 : \Gamma_+ \to \Gamma_-` — the α → 0 member of
    the albedo family ``α · (geometric link)``, realized structurally
    (the link is never built just to be zeroed).

    The action mints the codomain's zero from the BOUND shape, with the
    seam's payload convention: the space describes the STRUCTURAL axes,
    and any trailing axes the operand carries beyond the domain's rank
    (group/spatial payload at the trace seam) ride through unchanged —
    ``zeros(codomain.shape + payload)``. The transpose mints the
    domain's the same way (:math:`0^{\mathsf T} : \mathcal C \to
    \mathcal D` — under #276 A4 duality typing it consumes the
    codomain-cotangent and emits the domain-cotangent, and a zero is a
    zero on either side). This is
    :class:`~orpheus.numerics.elementary_operators.TraceRestrictionOperator`'s trailing-axes-intact convention,
    stated on a bound pair.
    The pre-amendment spelling reached these shapes through per-site
    ``codomain_zero``/``transpose_zero`` closures, which duplicated what
    the bound spaces already knew; deriving from the binding retires the
    closures AND their caller-side consistency check (the hook-vs-space
    agreement is now unspellable — single source).

    :meth:`adjoint` is algebra-closed and metric-free: the Hilbert
    adjoint of the zero map is the zero map of the swapped pair, in any
    inner product — no sandwich, no wrapper.

    Deliberately UNPARAMETERIZED (the :class:`~orpheus.numerics.elementary_operators.PermutationOperator`
    precedent): the ndarray-seam consumers (the boundary trace algebra)
    do not satisfy the ``Vector`` protocol statically, and the two
    ``FunctionSpace``\ s — not static type parameters — are this
    operator's identity.
    """

    def __init__(
        self,
        *,
        domain: "FunctionSpace",
        codomain: "FunctionSpace",
    ) -> None:
        if domain is None or codomain is None:  # a loud guard, not typing
            raise TypeError(
                "ZeroMorphism demands BOTH spaces at construction — the "
                "zero map is the one operator whose action cannot reveal "
                "its pair, so the pair must be declared. For the"
                " endomorphic zero at the operand's own space use"
                " ZeroOperator()."
            )
        self._domain = domain
        self._codomain = codomain

    @property
    def domain(self) -> "FunctionSpace":
        return self._domain

    @property
    def codomain(self) -> "FunctionSpace":
        return self._codomain

    def apply(self, x: object, /) -> "np.ndarray":
        r"""The zero of the CODOMAIN — bound structural axes + the
        operand's payload tail, in the operand's dtype."""
        arr = np.asarray(x)
        payload = arr.shape[len(self._domain.shape):]
        return np.zeros(self._codomain.shape + payload, dtype=arr.dtype)

    def apply_transpose(self, y: object, /) -> "np.ndarray":
        r"""The zero of the DOMAIN (the transpose lands there), payload
        riding through as in :meth:`apply`."""
        arr = np.asarray(y)
        payload = arr.shape[len(self._codomain.shape):]
        return np.zeros(self._domain.shape + payload, dtype=arr.dtype)

    @property
    def is_adjointable(self) -> bool:
        return True  # 0^* exists in any metric: the swapped zero map

    def adjoint(self) -> "ZeroMorphism":
        r"""``0^* : \mathcal C \to \mathcal D`` — the swapped pair's zero
        map (metric-free: the zero is self-dual up to the swap)."""
        return ZeroMorphism(domain=self._codomain, codomain=self._domain)

    # is_invertible inherits ``False``; no ``inverse()`` is declared —
    # the same Design-C surface as the endomorphic ZeroOperator.


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
    :class:`~orpheus.numerics.green_operator.GreenOperator` to the
    difference :math:`M - N` of its split; ``MatrixInverseOperator`` to
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

    * The TRANSPOSE axis, by the law
      :math:`(A^{-1})^{\mathsf T} = (A^{\mathsf T})^{-1}`:
      :meth:`apply_transpose` is the forward's transposed solve
      ``inner.solve_transpose``, and :attr:`is_adjointable` is
      :func:`transpose_invertible` of the forward (the SN sweep's reverse
      scan, a coupled block substitution, a diagonal leaf all answer
      through it). Read the other way, the inverse's transposed solve is
      the forward's transpose: :meth:`solve_transpose` is
      ``inner.apply_transpose`` and :attr:`is_transpose_invertible` is
      :func:`adjointable` of the forward.
    * :attr:`is_metric_free_adjoint` is the forward's: an operator that
      commutes with every diagonal metric has an inverse that does too, so
      the inverse of a multiplier needs no spaces for its ``.H``.

    A sibling whose transpose is not the forward's transposed solve
    overrides these: :class:`~orpheus.numerics.green_operator.GreenOperator`
    (its transpose is the Green of the dual sum), and
    ``MatrixInverseOperator`` (the ``trans=1`` backsolve on its own LU
    factors). ``.H`` itself is the base's: ``♯ ∘ dual ∘ ♭`` when
    :func:`posed`.

    (This wrap-delegate family is ONE of two kinds of inverse in the
    algebra: ALGEBRA-CLOSED structures invert into themselves — a
    permutation's inverse IS a permutation, a one-term combination's a
    one-term combination — and stay unwrapped as first-class forwards. The canonical
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

    def apply_transpose(self, x: Any, /) -> Any:
        r"""Return :math:`(A^{-1})^{\mathsf T} x = (A^{\mathsf T})^{-1} x` — the forward's transposed solve."""
        inner = cast("LinearOperator", self.inner)
        if not transpose_invertible(inner):
            raise MissingAdjoint(
                f"{type(self).__name__}.apply_transpose: the forward "
                f"{type(inner).__name__} is not transpose-invertible, so "
                f"(A^{{-1}})^T = (A^T)^{{-1}} has no realisation."
            )
        return inner.solve_transpose(x)

    @property
    def is_adjointable(self) -> bool:
        # (A^{-1})^T = (A^T)^{-1}: the inverse transposes iff its forward
        # realises the transposed solve.
        return transpose_invertible(cast("LinearOperator", self.inner))

    def solve_transpose(self, y: Any, /) -> Any:
        r"""Solve :math:`(A^{-1})^{\mathsf T} x = y`, i.e. return :math:`A^{\mathsf T} y` (the forward's transpose)."""
        inner = cast("LinearOperator", self.inner)
        if not adjointable(inner):
            raise MissingAdjoint(
                f"{type(self).__name__}.solve_transpose: the forward "
                f"{type(inner).__name__} has no transpose, so "
                f"((A^{{-1}})^T)^{{-1}} = A^T has no realisation."
            )
        return inner.apply_transpose(y)

    @property
    def is_transpose_invertible(self) -> bool:
        # ((A^{-1})^T)^{-1} = A^T: exists iff the forward transposes.
        return adjointable(cast("LinearOperator", self.inner))

    @property
    def is_metric_free_adjoint(self) -> bool:
        # An operator commuting with every diagonal metric has an inverse
        # that commutes with it too.
        return bool(getattr(self.inner, "is_metric_free_adjoint", False))


class InverseOperator(_KeyedEquality, InverseWrapMixin[_InvertibleForward], LinearOperator):
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
    :class:`~orpheus.numerics.green_operator.GreenOperator` for a named
    split :math:`M - N`;
    :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
    for the dense direct factorization) — this class serves any
    solve-backed forward with NO more specific named inverse: the
    value-bearing LEAVES (:class:`~orpheus.numerics.elementary_operators.DiagonalOperator`,
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

    def _algebra_key(self) -> Hashable:
        # A formal arrow is its kind and its inner (ruling 53's spec, Q5).
        return ("inverse", _algebra_key(self.inner))

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
    * **Spaces by AGREEMENT, not position** (:attr:`domain`): since the
      factors commute, a position-based rule would give an
      order-independent operator order-dependent spaces. Every factor
      that declares a space must declare the same one; silence
      contributes nothing; disagreement raises
      :class:`IncompatibleOperatorComposition`. This is the law
      :class:`OperatorSum` — the other commutative composite — already
      obeyed, shared as :func:`_agreed_space`.

    Parameters
    ----------
    ops : tuple of LinearOperator
        The tensor-product factors. Each MUST advertise an ``axis``
        attribute (or accept an ``axis`` kwarg in :meth:`apply`) and
        broadcast on every other axis. :class:`IdentityOperator`,
        :class:`~orpheus.numerics.elementary_operators.DiagonalOperator`, and any
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
        # Eager space agreement, for the same reason :class:`OperatorSum`
        # checks it and by the same shared law: the factors COMMUTE, so the
        # product's spaces cannot be a function of factor order.
        _agreed_space(self.ops, "domain", "TensorProductOperator")
        _agreed_space(self.ops, "codomain", "TensorProductOperator")

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        r"""The space the factors agree on — see :func:`_agreed_space`.

        A factor that declares nothing (the group-axis
        :class:`IdentityOperator` every shipped boundary law carries) leaves
        the binding to the factor that does, so ``K_ω ⊗ I`` is bound exactly
        where ``K_ω`` is. Before G6.3 step 8.0 this returned the base's
        ``None``, which meant the binding was real at the inner factor and
        INVISIBLE at the object a realizer hands out — and, because
        :class:`AdjointOperator` reads the spaces to apply the metrics, it
        also meant ``(K_ω ⊗ I).H`` silently degraded to the Euclidean
        transpose (`[M]` 87 % relative error against the weighted adjoint on
        the Lambertian, exact only for the specular mirror, whose metric
        cancels).
        """
        return _agreed_space(self.ops, "domain", "TensorProductOperator")

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        """The space the factors agree on — see :attr:`domain`."""
        return _agreed_space(self.ops, "codomain", "TensorProductOperator")

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

    @property
    def is_metric_free_adjoint(self) -> bool:
        # A tensor product of metric-free operators is metric-free (derived).
        return all(op.is_metric_free_adjoint for op in self.ops)




class InverseMetricOperator(LinearOperator):
    r"""A :class:`~orpheus.numerics.space.FunctionSpace`'s inverse metric, as an OPERATOR.

    :math:`G^{+} : V^{*} \to V` — the adapter that lets a space's metric
    enter the operator algebra instead of only being *applied* to arrays.

    **Why this exists.** A space's metric enters the operator algebra
    through this adapter — the trace metrics' inverses, a degenerate
    metric's Moore–Penrose face — wherever :math:`G^{+}` of a SPACE is the
    operator wanted. (Until CS4c step 6 item 6.2c-ii the frame's projector
    was spelled through it too, as
    ``conjugate(frame, InverseMetricOperator(frame.gram))`` over a metric-twin
    of the test space; the frame now owns that factor as a typed arrow,
    :meth:`~orpheus.numerics.frame_faces.CrossGramInverse.of` — an endomorphism
    of a twin space cannot compose with the faces once the metric enters
    space identity.)

    **The arithmetic is the SPACE's, not ours** (Cardinal Rule 2): every
    apply delegates to
    :meth:`~orpheus.numerics.space.FunctionSpace.apply_inverse_metric`, so
    the Moore–Penrose masking on a degenerate metric —
    :math:`1/G` where :math:`G \neq 0`, **0** on
    :math:`\ker G` — is single-sourced there and cannot drift.  That
    masking is not academic: `[M]` the SN trace metric
    :math:`G = |\Omega\cdot\hat n|\,w_n` is **exactly zero** on tangential
    ordinates — 50 % of rows under ``product(4,4)``, 16 % under
    ``lebedev(11)``, 0 % under ``level_symmetric``.

    ⛔ **Not invertible, and it says so by ABSENCE.** On a degenerate
    metric this is a pseudo-inverse, so ``G⁺G ≠ I`` and there is no
    ``inverse()`` method to call — the
    :class:`~orpheus.numerics.elementary_operators.TraceRestrictionOperator` spelling, not a raising stub
    (:attr:`is_invertible` stays ``False``).  ⚠ Consequently
    ``InverseMetricOperator(space)`` composed with a forward metric does
    **not** cancel; if you want the round trip, assert it on a field you
    know is off the null space.

    Self-adjoint: the metric is a symmetric positive-semi-definite form —
    a diagonal weight, or since P7 a dense
    :class:`~orpheus.numerics.metric.DenseMetric` (symmetry guarded at its
    construction) — so :math:`(G^{+})^{\mathsf T} = G^{+}` and
    ``apply_transpose`` is ``apply`` for every realization.

    Parameters
    ----------
    space :
        The space whose metric to invert.  Serves as BOTH ``domain`` and
        ``codomain`` — the metric is an endomorphism of the carrier's
        shape (it re-weights, it does not move between spaces), and
        binding both ends is what lets
        :class:`OperatorProduct`'s compatibility guard check the
        composition.
    """

    def __init__(self, space: FunctionSpace) -> None:
        self._space = space

    @property
    def space(self) -> FunctionSpace:
        """The space whose metric this inverts."""
        return self._space

    @property
    def domain(self) -> Optional[FunctionSpace]:
        return self._space

    @property
    def codomain(self) -> Optional[FunctionSpace]:
        return self._space

    def apply(self, x: np.ndarray) -> np.ndarray:
        return self._space.apply_inverse_metric(x)

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        # Self-adjoint: a real diagonal weight is its own transpose.
        return self._space.apply_inverse_metric(x)

    @property
    def is_adjointable(self) -> bool:
        return True

    def __repr__(self) -> str:
        return f"InverseMetricOperator({self._space!r})"

