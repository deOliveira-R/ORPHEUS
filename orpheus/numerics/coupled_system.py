r"""The N-system coupled block machinery: ``CoupledField`` / ``CoupledSpace`` / ``CoupledOperator``.

A *coupled system* is N sub-systems solved together: the state is the direct
sum of the systems' own fields, and the operator is the N × N grid of blocks

.. math::

    A \;=\; \begin{bmatrix} A_{11} & \cdots & A_{1N} \\
                            \vdots &        & \vdots \\
                            A_{N1} & \cdots & A_{NN} \end{bmatrix},
    \qquad A_{ij} : V_j \to W_i ,

with the diagonal blocks the systems' self-operators and the off-diagonals
the couplings. This module realizes that object **semantics-agnostically**:
nothing here knows transport, rays, or meshes — the members are any fields
satisfying the ravellable protocol, the member spaces any
:class:`~orpheus.numerics.space.FunctionSpace`. The ψ½ radial-characteristic
2×2 (System A = SN transport, System B = the curvilinear ray closure) is
**instance #1**, wired in ``orpheus.sn`` by the coupled-block campaign
(``.claude/plans/archive/coupled_block_operator_campaign.md``, Phase B).

Why a typed block operator (the rejected route)
===============================================

The alternative — a flat :class:`~orpheus.numerics.operator.OperatorSum` of
same-space operators, each padding the blocks it does not touch with
present-zeros — was REJECTED (campaign re-scope, 2026-07-10): padding keeps
**wrong multiplications representable**. A padded block accepts any composite;
nothing at the type level stops a ray operator from receiving a bulk field or
an emission from landing in the trace slot, and ``system_role`` tags are
runtime metadata, not type constraints. The honest object is a **typed block
operator over a typed block vector**, where the block matvec

.. math:: y_i \;=\; \sum_j A_{ij}\, x_j

is the *only* spelling and every term type-checks per block
(``coding-elegance`` Pattern 1 ∘ Pattern 4: the algebra is the syntax, and a
wrong pairing is unconstructable). Block existence doubles as **coupling
sparsity**: a missing block (``None``) IS the zero map, structurally — no
zero-padding arithmetic ever runs.

The three consumption modes (the ``Op → Mat`` functor at block level)
=====================================================================

=============  ===========================================================
Mode           Realization here
=============  ===========================================================
``apply``      the block matvec ``y_i = Σ_j A_ij x_j`` (Krylov action)
``assemble``   each block's sparse assembly scattered at its
               ``(row_i, col_j)`` block offset into ONE flat matrix
               (:func:`scipy.sparse.block_array`), closing over
               :class:`~orpheus.numerics.assembled_operator.SparseAssembledOperator`
``solve``      the structure-keyed DIRECT solve (step 5): a
               block-triangular grid runs back/forward SUBSTITUTION
               through the diagonal blocks' own ``solve`` verbs
               (matrix-free); a non-triangular square grid
               materializes and LU-factors —
               :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
               via :meth:`CoupledOperator.inverse`. The ITERATIVE
               splitting solve (block-Jacobi / block-G-S over
               ``A = M − N``) deliberately stays with the drivers
               (``SourceIteration``) — convergence is spectral
               (ρ(M⁻¹N) < 1), never a structural capability
=============  ===========================================================

**The offsets ARE the block structure.** ``assemble`` needs a local→global
DOF map, and :class:`CoupledSpace` provides it for free
(:attr:`CoupledSpace.system_slices`): member ``i``'s flat DOFs occupy the
contiguous slice at the cumulative offset of the members before it — the
same layout :meth:`CoupledField.to_flat` packs. This is the **scoped
realization of the deferred ``LocalToGlobalMap``** (the 2-P0 attacker ruling
at ``orpheus/numerics/assembled_operator.py`` deferred a reified map "until a
structured consumer arrives" — the block grid is that structured consumer,
and the map it needs is exactly the member-offset table, nothing more).

The Hilbert adjoint comes FREE (Mode-12 closure by construction)
================================================================

``CoupledOperator`` implements only the Euclidean :meth:`~CoupledOperator.
apply_transpose` — the transposed grid ``(Aᵀ)_{ji} = (A_{ij})ᵀ`` — and
carries :class:`CoupledSpace` domain/codomain whose metric methods dispatch
member-wise. The metric adjoint ``A.H = G⁺ Aᵀ G`` is then realized ONCE by
the existing :class:`~orpheus.numerics.operator._AdjointOperator` wrapper
(which calls ``codomain.apply_metric`` before the transpose and
``domain.apply_inverse_metric`` after). No adjoint code lives here — which is
precisely what keeps the block adjoint Mode-12-closed: a hand-rolled
"Euclidean block ``.H``" that skips the metric conjugation is the ERR-067
reopening the campaign's A2 reciprocity gates exist to red
(``vv-principles`` Mode 12; test tooth M-ADJ-metric).

Members and the ravellable protocol
===================================

A member ("one system's field") is any object satisfying the **ravellable
protocol** already established at :mod:`orpheus.numerics.iteration` —
``to_flat()`` instance method + ``from_flat(flat, template)`` class-level
factory — plus ``copy()`` and the vector dunders (``+ − · /``), i.e. exactly
the surface of the transport ``Composite`` family. :class:`CoupledField`
itself satisfies the same protocol (its :meth:`~CoupledField.to_flat`
concatenates the members' flats in system order), so the scipy-Krylov
boundary (:func:`orpheus.numerics.iteration._as_scipy_linop`), the
:class:`~orpheus.numerics.flat_operator.FlattenedOperator` adapter, and every
``restart = n_dof = template.to_flat().size`` sizing site track the coupled
dimension **automatically** (the ERR-053 GMRES-truncation family is closed by
this conformance, not by per-site edits).

Role semantics stay on the members: the machinery never adds arrays — every
``+`` it evaluates is a member ``+``, so the member family's fiber guards
(class/space/mesh) and units law apply unchanged (the same delegation
discipline as ``Composite._map_binary``; the affine torsor and displacement
minting this delegation once carried retired at campaign 1 CS3 — flux lives
in V).

References
----------

* ``.claude/plans/archive/coupled_block_operator_campaign.md`` — the campaign plan:
  the 2×2 posing, the re-scope ruling (typed grid over present-zero padding),
  Phase B.2a (this module), B.2b–d (the ψ½ instance re-type + builder +
  driver wire).
* ``.claude/agent-memory/test-architect/coupled_operator_step4_verification.md``
  — the gate spec this module's suite realizes (M1 block matvec, M2
  assemble≡probe, M3 type-safety, M4 block-``.H`` Mode-12, M5 system role).
* Trefethen & Bau (1997), *Numerical Linear Algebra*, §1 — the Hilbert
  adjoint :math:`G^{-1}A^{\mathsf T}G` vs. the representation transpose.
* Hackbusch (2016), *Iterative Solution of Large Sparse Systems*, §11 —
  block partitionings and block relaxation (the step-5 solve consumers).
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from functools import reduce
from typing import TYPE_CHECKING, Any, Callable, Optional, Protocol, Sequence, cast

import numpy as np
from scipy import sparse

from orpheus.numerics.assembled_operator import SparseAssembledOperator
from orpheus.numerics.operator import (
    IncompatibleOperatorComposition,
    InverseWrapMixin,
    LinearOperator,
    MatrixTooLarge,
    MissingAdjoint,
    MissingAssembly,
    NotInvertible,
    SystemRole,
    adjointable,
    assemblable,
)
from orpheus.numerics.space import FunctionSpace

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator

__all__ = [
    "CoupledField",
    "CoupledOperator",
    "CoupledSpace",
    "CoupledSubstitutionOperator",
    "SystemField",
]


class SystemField(Protocol):
    r"""One system's complete field — the member contract of :class:`CoupledField`.

    Structural (duck-typed): the transport ``Composite`` family (``FullField``,
    ``RadialCharacteristicField``) satisfies it without an import edge out
    of the numerics layer — the same deliberate decoupling as the ravellable
    protocol at :mod:`orpheus.numerics.iteration` (whose ``to_flat()`` /
    ``from_flat(flat, template)`` pair this protocol embeds).

    Beyond the two members declared here, a member must support the vector
    dunders (``+``, ``-``, unary ``-``, scalar ``*`` / ``/``) and the
    class-level ``from_flat(flat, template)`` factory; those stay duck-typed
    (not declared) TODAY only as inherited shape. The original reason — the
    member family's affine-torsor signatures (``flux − flux → displacement``,
    ``flux + displacement → flux``) fit neither
    :class:`~orpheus.numerics.vector.Vector`'s ``Self + Self → Self``
    spelling nor any simple declaration — was DISSOLVED at campaign 1 CS3
    (2026-08-19): every member family now carries exactly the ``Self``-typed
    vector dunders. The machinery still delegates through ``Any``-typed
    lambdas (the ``iteration.py`` posture) without inspecting them, and
    numerics still carries THREE member-contract concepts — ``Vector`` (the
    named arithmetic dunders), the ``iteration.py`` ravellable pair (ad-hoc
    private), and this protocol. **The collapse trigger is therefore LIVE**:
    mint the named ravellable-Protocol home (the ``Vector`` promotion
    precedent for the arithmetic half) and converge all three onto it —
    tracked as #391.
    """

    def to_flat(self) -> "NDArray": ...

    def copy(self) -> "SystemField": ...


def _member_from_flat(member: "SystemField", piece: "NDArray") -> "SystemField":
    r"""Rebuild one member from its flat slice — the ravellable-protocol pair.

    ``type(member).from_flat(piece, member)`` is the class-level half of the
    duck-typed ravellable protocol (see :class:`SystemField`); the cast is the
    ONE documented static seam. DELIBERATELY coextensive with (not delegated
    to) :func:`orpheus.numerics.iteration._unravel_like`'s ravellable branch —
    that helper is an ad-hoc ``_``-private the codebase is migrating away
    from; converge BOTH into the named ravellable-Protocol home when it is
    minted (the :class:`SystemField` collapse trigger), which also deletes
    this helper and its cast.
    """
    factory = cast("Any", type(member))
    return factory.from_flat(piece, member)


def _is_system_field(member: object) -> bool:
    r"""Detect the member contract (the :func:`~orpheus.numerics.iteration
    ._is_ravellable` pair + ``copy``), duck-typed.

    The first two conjuncts are deliberately coextensive with
    ``_is_ravellable`` (same migration note as :func:`_member_from_flat` —
    converge at the named ravellable-Protocol home, don't import the
    ``_``-private)."""
    return (
        hasattr(member, "to_flat")
        and hasattr(type(member), "from_flat")
        and hasattr(member, "copy")
    )


# ───────────────────────────────────────────────────────────────────────
# The block vector
# ───────────────────────────────────────────────────────────────────────


@dataclass(frozen=True, kw_only=True, eq=False)
class CoupledField:
    r"""The N-system block vector: one member field per coupled system.

    The state of a coupled solve — e.g. the ψ½ instance's
    ``[ψ_A, ψ_B] = [FullField, RadialCharacteristicField]``. The vector
    algebra (``±``, scalar ``·``, :meth:`copy`, the flat protocol) is realized
    ONCE here by member-wise delegation: every ``+`` is a member ``+``, so the
    member family's fiber guards and units law apply unchanged — this class
    adds structure, never arithmetic (the affine torsor this note once named
    retired at campaign 1 CS3).

    ``eq=False`` deliberately: an ndarray-bearing aggregate has no
    well-defined ``==`` (element-wise comparison is ambiguous as a truth
    value), so identity semantics are kept — compare members explicitly.

    Parameters
    ----------
    systems : tuple[SystemField, ...]
        The member fields, in system order (the order fixes the flat layout
        and must match the paired :class:`CoupledSpace`). At least one member;
        a 1-tuple is the legitimate degenerate (an uncoupled system — e.g. a
        non-ray-carrying mesh has no System B, so its "coupled" state is the
        1-system vector, and a System-B slot simply does not exist to be
        filled wrongly).
    """

    systems: tuple[SystemField, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.systems, tuple):
            raise TypeError(
                f"CoupledField: systems must be a tuple of member fields; got "
                f"{type(self.systems).__name__}."
            )
        if len(self.systems) == 0:
            raise ValueError(
                "CoupledField: at least one system is required — a 0-system "
                "coupled field carries no state."
            )
        for i, member in enumerate(self.systems):
            if not _is_system_field(member):
                raise TypeError(
                    f"CoupledField: system {i} ({type(member).__name__}) does "
                    f"not satisfy the member contract (to_flat / "
                    f"from_flat(flat, template) / copy — the ravellable "
                    f"protocol plus copy)."
                )

    @property
    def n_systems(self) -> int:
        r"""The number of coupled systems (the block-vector arity)."""
        return len(self.systems)

    # ── Algebra (member-wise delegation; roles live on the members) ────

    def _check_partner(self, other: object) -> "CoupledField":
        r"""Reject a partner that is not an arity-matched :class:`CoupledField`.

        Container level only (the ``Composite._check_partner`` discipline):
        member-level role/type/mesh law is the members' single source of
        truth — a member pre-check here would be a second spelling of the
        members' own law.
        """
        if not isinstance(other, CoupledField):
            raise TypeError(
                f"CoupledField arithmetic requires a same-class partner; got "
                f"{type(other).__name__}."
            )
        if other.n_systems != self.n_systems:
            raise ValueError(
                f"CoupledField arithmetic requires matching system arity; got "
                f"{self.n_systems} vs {other.n_systems} systems."
            )
        return other

    def _map_binary(
        self, other: "CoupledField", op: "Callable[[Any, Any], Any]",
    ) -> "CoupledField":
        # ``Any``-typed op: the member dunders are the duck-typed half of the
        # member contract (see SystemField). ``replace`` re-runs
        # ``__post_init__`` so the member-contract invariant re-fires for free
        # (coding-elegance Pattern 4 ∩ Pattern 2).
        return replace(
            self,
            systems=tuple(
                op(a, b) for a, b in zip(self.systems, other.systems)
            ),
        )

    def _map_unary(self, op: "Callable[[Any], Any]") -> "CoupledField":
        return replace(self, systems=tuple(op(m) for m in self.systems))

    def __add__(self, other: "CoupledField") -> "CoupledField":
        return self._map_binary(self._check_partner(other), lambda a, b: a + b)

    def __sub__(self, other: "CoupledField") -> "CoupledField":
        return self._map_binary(self._check_partner(other), lambda a, b: a - b)

    def __neg__(self) -> "CoupledField":
        return self._map_unary(lambda m: -m)

    def __mul__(self, scalar: float) -> "CoupledField":
        return self._map_unary(lambda m: m * float(scalar))

    def __rmul__(self, scalar: float) -> "CoupledField":
        return self.__mul__(scalar)

    def __truediv__(self, scalar: float) -> "CoupledField":
        return self._map_unary(lambda m: m / float(scalar))

    def copy(self) -> "CoupledField":
        r"""A member-wise deep copy (owned member buffers)."""
        return self._map_unary(lambda m: m.copy())

    # ── Flat protocol (the ravellable pair — Krylov / ERR-053 closure) ──

    def to_flat(self) -> "NDArray":
        r"""Pack the members into ONE flat vector, in system order.

        The layout is ``concat(system_0.to_flat(), …, system_{N−1}.to_flat())``
        — the same member-offset table :attr:`CoupledSpace.system_slices`
        exposes. Every ``n_dof = template.to_flat().size`` sizing site (GMRES
        ``restart``, Krylov workspace — the ERR-053 family) therefore counts
        BOTH systems automatically.
        """
        return np.concatenate(
            [np.asarray(m.to_flat(), dtype=float) for m in self.systems],
        )

    @classmethod
    def from_flat(
        cls, flat: "NDArray", template: "CoupledField",
    ) -> "CoupledField":
        r"""Reconstruct a coupled field from a flat vector + template.

        The inverse of :meth:`to_flat` (the ravellable-protocol factory half):
        the flat vector is sliced by the template members' flat sizes and each
        piece delegated to the member's own ``from_flat``. The result carries
        the template's concrete class (via :func:`dataclasses.replace`).
        """
        flat = np.asarray(flat, dtype=float).ravel()
        sizes = [np.asarray(m.to_flat()).size for m in template.systems]
        total = int(sum(sizes))
        if flat.size != total:
            raise ValueError(
                f"CoupledField.from_flat: flat.size = {flat.size} does not "
                f"match template total size "
                f"{' + '.join(str(s) for s in sizes)} = {total}"
            )
        members = []
        start = 0
        for member, size in zip(template.systems, sizes):
            members.append(_member_from_flat(member, flat[start : start + size]))
            start += size
        return replace(template, systems=tuple(members))

    def __repr__(self) -> str:
        names = ", ".join(type(m).__name__ for m in self.systems)
        return f"CoupledField[{names}]"


# ───────────────────────────────────────────────────────────────────────
# The block space
# ───────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class CoupledSpace(FunctionSpace["CoupledField"]):
    r"""The N-system direct sum :math:`V = V_1 \oplus \cdots \oplus V_N`.

    The carrier space of a :class:`CoupledOperator`: one member
    :class:`~orpheus.numerics.space.FunctionSpace` per system, in the SAME
    system order as the paired :class:`CoupledField`. Like
    :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace` (the
    2-block precedent one level down), this space owns ONLY the direct-sum
    *structure*: every metric primitive (:meth:`apply_metric` /
    :meth:`apply_inverse_metric` / :meth:`inner_product`) dispatches to the
    member space on the member field — no new metric arithmetic. That is what
    keeps the block Hilbert adjoint Mode-12-closed for free: the generic
    :class:`~orpheus.numerics.operator._AdjointOperator` conjugates through
    THESE methods, so the composite metric ``G = diag(G_1, …, G_N)`` is
    applied per member by the spaces that own each ``G_i``.

    Identity is the inherited ``(name, shape)`` tuple with the name derived
    from the members (``"coupled(a ⊕ b)"`` — the
    :class:`~orpheus.numerics.space.TensorProductSpace` naming precedent) and
    ``shape`` the flat direct-sum dimension.

    Parameters
    ----------
    name, shape, inner_product_weights
        Inherited. The composite's own ``inner_product_weights`` stays
        ``None`` — the metric is carried per member, and the member-wise
        overrides never read the base slot.
    systems : tuple[FunctionSpace, ...]
        The member spaces (``compare=False`` leaf metadata, not part of the
        identity). Each must be the space of the matching
        :class:`CoupledField` member — its carrier-generic metric surface is
        called ON that member field.
    zeros_factory : Callable[[], CoupledField], optional
        Mints this space's ZERO element — a fresh (owned-buffer)
        :class:`CoupledField` of zeros per call. ``compare=False`` leaf
        metadata like ``systems``; wired by the instance builder (which
        holds the member classes/mesh the fields are minted from — a
        FunctionSpace deliberately is not a field factory, so the zero
        exemplar must be SUPPLIED). This is the typed-carrier
        materialization seam (step 5): :meth:`CoupledOperator.as_matrix`'s
        basis probe and
        :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`'s
        typed return template both read it — a typed block operator whose
        space carries no zero exemplar stays matrix-free (the base
        ``as_matrix`` "honest scope" note).
    """

    systems: tuple[FunctionSpace, ...] = field(
        default=(), repr=False, compare=False,
    )
    zeros_factory: Optional[Callable[[], "CoupledField"]] = field(
        default=None, repr=False, compare=False,
    )

    # ── Equality / hashing inherited from FunctionSpace ───────────────
    #
    # Same rationale as TensorProductSpace / FullFieldSpace: the frozen
    # dataclass would regenerate __eq__ over every field (including the
    # ndarray weights slot); explicit delegation restores the (name, shape)
    # identity convention. ``systems`` is already compare=False.

    def __eq__(self, other: object) -> bool:
        return FunctionSpace.__eq__(self, other)

    def __hash__(self) -> int:
        return FunctionSpace.__hash__(self)

    def __repr__(self) -> str:
        return f"CoupledSpace({self.name!r}, shape={self.shape})"

    @classmethod
    def from_systems(
        cls,
        systems: "Sequence[FunctionSpace]",
        *,
        zeros: "Callable[[], CoupledField] | None" = None,
    ) -> "CoupledSpace":
        r"""Build the coupled space from the member spaces, in system order.

        Derives ``name = "coupled(" + " ⊕ ".join(member names) + ")"`` and
        ``shape = (Σ prod(member.shape),)`` (the flat direct-sum dimension —
        ``prod`` per member so multi-axis member shapes flatten honestly).
        A 1-member coupling is the legitimate degenerate (the uncoupled
        system). ``zeros`` wires the optional zero-element factory (the
        typed-carrier materialization seam — see the class docstring).
        """
        members = tuple(systems)
        if len(members) == 0:
            raise ValueError(
                "CoupledSpace.from_systems: at least one member space is "
                "required — a 0-system coupled space carries no DOFs."
            )
        name = "coupled(" + " ⊕ ".join(s.name for s in members) + ")"
        total = int(sum(int(np.prod(s.shape)) for s in members))
        return cls(
            name=name, shape=(total,), systems=members, zeros_factory=zeros,
        )

    def zeros(self) -> "CoupledField":
        r"""Mint this space's ZERO element — a fresh zero :class:`CoupledField`.

        The space's distinguished element, realized through the
        builder-wired :attr:`zeros_factory` (a FunctionSpace is not a field
        factory — the member classes live above this layer, so the exemplar
        is supplied at construction). Each call returns a FRESH field with
        owned buffers (factory semantics — callers may write into it).
        Arity-checked against the member spaces; loud when unwired.
        """
        if self.zeros_factory is None:
            raise RuntimeError(
                f"CoupledSpace {self.name!r} carries no zero-element "
                f"factory — supply zeros= at from_systems (the "
                f"typed-carrier materialization seam; without it a typed "
                f"block operator over this space stays matrix-free)."
            )
        minted = self.zeros_factory()
        self._check_arity(minted)
        return minted

    def _require_systems(self) -> tuple[FunctionSpace, ...]:
        r"""Return the member spaces, guarding the bare-constructor footgun.

        The ``systems`` field defaults to ``()`` (the ``compare=False``
        dataclass-field convention), so a space built via the bare constructor
        instead of :meth:`from_systems` would dispatch to nothing and silently
        no-op. Fail at the boundary with intent (parse-don't-validate; the
        :meth:`FullFieldSpace._require_blocks` precedent).
        """
        if len(self.systems) == 0:
            raise RuntimeError(
                "CoupledSpace has no member spaces; build it via "
                "CoupledSpace.from_systems((space_1, …, space_N)), not the "
                "bare dataclass constructor."
            )
        return self.systems

    @property
    def n_systems(self) -> int:
        r"""The number of coupled systems (the block arity)."""
        return len(self._require_systems())

    @property
    def system_slices(self) -> tuple[slice, ...]:
        r"""Member ``i``'s DOF range inside the flat direct-sum layout.

        The cumulative-offset table ``(slice(0, n_1), slice(n_1, n_1+n_2),
        …)`` — the SAME layout :meth:`CoupledField.to_flat` packs, so a flat
        Krylov iterate, an assembled matrix row/column, and a member field
        agree on where system ``i`` lives. **This IS the scoped
        ``LocalToGlobalMap``** the assembly layer deferred until a structured
        consumer arrived (2-P0 attacker ruling,
        :mod:`orpheus.numerics.assembled_operator`): the block structure
        provides the offsets, nothing more is reified.
        """
        members = self._require_systems()
        slices = []
        start = 0
        for s in members:
            size = int(np.prod(s.shape))
            slices.append(slice(start, start + size))
            start += size
        return tuple(slices)

    # ── Direct-sum metric dispatch (member-wise, on a CoupledField) ────

    def _check_arity(self, x: "CoupledField") -> tuple[FunctionSpace, ...]:
        members = self._require_systems()
        if x.n_systems != len(members):
            raise ValueError(
                f"CoupledSpace: the coupled field carries {x.n_systems} "
                f"systems but this space couples {len(members)} — the "
                f"field/space pairing is inconsistent."
            )
        return members

    def apply_metric(self, x: "CoupledField") -> "CoupledField":
        r"""Apply the block-diagonal Hilbert metric ``G = diag(G_1, …, G_N)``.

        Member ``i`` is weighted by ITS space's metric — delegated to
        ``systems[i].apply_metric`` on the member field (the carrier-generic
        surface: a composite member space applies its own per-block metrics
        recursively). No metric arithmetic lives on the coupling.
        """
        members = self._check_arity(x)
        return replace(
            x,
            systems=tuple(
                s.apply_metric(m) for s, m in zip(members, x.systems)
            ),
        )

    def apply_inverse_metric(self, x: "CoupledField") -> "CoupledField":
        r"""Apply the block-diagonal pseudo-inverse metric ``G⁺`` member-wise."""
        members = self._check_arity(x)
        return replace(
            x,
            systems=tuple(
                s.apply_inverse_metric(m) for s, m in zip(members, x.systems)
            ),
        )

    def inner_product(self, x: "CoupledField", y: "CoupledField") -> float:
        r"""The direct-sum inner product
        :math:`\langle x, y\rangle_G = \sum_i \langle x_i, y_i\rangle_{G_i}`."""
        members = self._check_arity(x)
        self._check_arity(y)
        return float(
            sum(
                s.inner_product(mx, my)
                for s, mx, my in zip(members, x.systems, y.systems)
            )
        )


# ───────────────────────────────────────────────────────────────────────
# The block operator
# ───────────────────────────────────────────────────────────────────────


class CoupledOperator(LinearOperator["CoupledField", "CoupledField"]):
    r"""The N×N typed block grid ``A_ij : System_j → System_i``.

    The block matvec ``y_i = Σ_j A_ij x_j`` is the ONLY spelling of the
    coupled action, and it type-checks per block: at construction every block
    that declares spaces is checked against the coupled domain/codomain
    members at its ``(i, j)`` position
    (:class:`~orpheus.numerics.operator.IncompatibleOperatorComposition` on
    mismatch — a mis-placed coupling is unconstructable, not a runtime shape
    crash three calls later), and at apply time the iterate must be an
    arity-matched :class:`CoupledField`. A ``None`` block IS the zero map —
    coupling sparsity is structural, no zero arithmetic runs.

    Structural axes: :attr:`is_adjointable` / :attr:`is_assemblable` derive
    recursively from the blocks (the
    :class:`~orpheus.numerics.operator.OperatorSum` closure-law discipline);
    ``A.H`` is realized by the generic
    :class:`~orpheus.numerics.operator._AdjointOperator` over THIS class's
    Euclidean :meth:`apply_transpose` (the transposed grid) and the
    :class:`CoupledSpace` member-wise metrics — see the module docstring's
    Mode-12 note. :attr:`is_invertible` is structure-keyed (step 5), two
    DIRECT routes: block-triangular grids with invertible ``solve``-bearing
    diagonals run the substitution (:meth:`solve` /
    :class:`CoupledSubstitutionOperator`); non-triangular square
    materializable grids take the LU EXTRACT
    (:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`).
    The ITERATIVE splitting solve stays with the drivers — spectral, never
    structural.

    Parameters
    ----------
    blocks : Sequence[Sequence[LinearOperator | None]]
        The row-major grid: ``blocks[i][j]`` maps domain system ``j`` to
        codomain system ``i``; ``None`` is the structural zero map. Every row
        and every column must carry at least one block — an all-``None`` row
        (or column) leaves that system's output (input) undefined; drop the
        system from the coupling instead of spelling a zero system.
    domain, codomain : CoupledSpace
        The coupled member spaces the grid is typed against (REQUIRED — the
        offsets, metrics, and per-block checks all read them). A square
        coupling passes the same space twice; rectangular grids (restriction /
        prolongation stacks) are honest and admitted.
    """

    # The grid spans systems by construction (campaign step 4a lattice).
    # Plain (unannotated) class attribute per the LinearOperator base note.
    system_role = SystemRole.COUPLED

    def __init__(
        self,
        blocks: "Sequence[Sequence[Optional[LinearOperator]]]",
        *,
        domain: "CoupledSpace",
        codomain: "CoupledSpace",
    ) -> None:
        domain_members = domain._require_systems()
        codomain_members = codomain._require_systems()
        grid = tuple(tuple(row) for row in blocks)
        n_rows, n_cols = len(codomain_members), len(domain_members)
        if len(grid) != n_rows:
            raise ValueError(
                f"CoupledOperator: the grid has {len(grid)} rows but the "
                f"codomain couples {n_rows} systems — one row per codomain "
                f"system."
            )
        for i, row in enumerate(grid):
            if len(row) != n_cols:
                raise ValueError(
                    f"CoupledOperator: row {i} has {len(row)} blocks but the "
                    f"domain couples {n_cols} systems — one column per domain "
                    f"system."
                )
        for i, row in enumerate(grid):
            if all(b is None for b in row):
                raise ValueError(
                    f"CoupledOperator: row {i} has no blocks — system {i}'s "
                    f"output would be undefined. Drop the system from the "
                    f"coupling instead of spelling an all-zero row."
                )
        for j in range(n_cols):
            if all(row[j] is None for row in grid):
                raise ValueError(
                    f"CoupledOperator: column {j} has no blocks — system "
                    f"{j}'s input would never be read. Drop the system from "
                    f"the coupling instead of spelling an all-zero column."
                )
        # Per-block space typing (Pattern 4): a block that DECLARES spaces
        # must sit at a grid position whose domain/codomain members match.
        # ``None``-spaced blocks skip the check — the house convention for
        # legacy space-less operators (LinearOperator.domain docstring).
        for i, row in enumerate(grid):
            for j, block in enumerate(row):
                if block is None:
                    continue
                block_domain = getattr(block, "domain", None)
                if block_domain is not None and block_domain != domain_members[j]:
                    raise IncompatibleOperatorComposition(
                        f"CoupledOperator: block ({i}, {j}) "
                        f"({type(block).__name__}) declares domain "
                        f"{block_domain.name!r} but column {j} carries "
                        f"system space {domain_members[j].name!r} — the "
                        f"block is placed at the wrong column (or typed "
                        f"against the wrong system)."
                    )
                block_codomain = getattr(block, "codomain", None)
                if (
                    block_codomain is not None
                    and block_codomain != codomain_members[i]
                ):
                    raise IncompatibleOperatorComposition(
                        f"CoupledOperator: block ({i}, {j}) "
                        f"({type(block).__name__}) declares codomain "
                        f"{block_codomain.name!r} but row {i} carries system "
                        f"space {codomain_members[i].name!r} — the block is "
                        f"placed at the wrong row (or typed against the "
                        f"wrong system)."
                    )
        self._blocks = grid
        self._domain = domain
        self._codomain = codomain

    # ── Structure ──────────────────────────────────────────────────────

    @property
    def blocks(self) -> tuple[tuple[Optional[LinearOperator], ...], ...]:
        r"""The row-major block grid (``None`` = structural zero map)."""
        return self._blocks

    @property
    def n_rows(self) -> int:
        return len(self._blocks)

    @property
    def n_cols(self) -> int:
        return len(self._blocks[0])

    @property
    def domain(self) -> "CoupledSpace":
        return self._domain

    @property
    def codomain(self) -> "CoupledSpace":
        return self._codomain

    def _present_blocks(self) -> list[LinearOperator]:
        return [b for row in self._blocks for b in row if b is not None]

    # ── The block matvec (apply) and its Euclidean transpose ───────────

    def _check_iterate(self, x: object, *, arity: int, side: str) -> "CoupledField":
        if not isinstance(x, CoupledField):
            raise TypeError(
                f"CoupledOperator.{side}: the block matvec is the only "
                f"spelling — expected a CoupledField with {arity} systems, "
                f"got {type(x).__name__}."
            )
        if x.n_systems != arity:
            raise ValueError(
                f"CoupledOperator.{side}: the coupled field carries "
                f"{x.n_systems} systems but this grid expects {arity}."
            )
        return x

    def apply(self, x: "CoupledField", /) -> "CoupledField":
        r"""The block matvec :math:`y_i = \sum_j A_{ij}\, x_j`.

        Each present block consumes ITS domain member and emits ITS codomain
        member; the row sum is a member ``+`` (role semantics on the member
        family). ``None`` blocks contribute nothing — structurally, not as a
        computed zero.
        """
        x = self._check_iterate(x, arity=self.n_cols, side="apply")
        outs = []
        for row in self._blocks:
            contributions = [
                block.apply(x.systems[j])
                for j, block in enumerate(row)
                if block is not None
            ]
            # Non-empty by the all-None-row construction guard.
            outs.append(reduce(lambda a, b: a + b, contributions))
        return CoupledField(systems=tuple(outs))

    def apply_transpose(self, y: "CoupledField", /) -> "CoupledField":
        r"""The Euclidean transposed grid :math:`(A^{\mathsf T})_{ji} = (A_{ij})^{\mathsf T}`.

        ``x_j = Σ_i (A_ij)ᵀ y_i`` — the representation transpose only; the
        metric conjugation that makes the Hilbert adjoint lives on
        :class:`~orpheus.numerics.operator._AdjointOperator` +
        :class:`CoupledSpace` (see the module docstring). Reachable per the
        :attr:`is_adjointable` advertisement (all present blocks adjointable).
        """
        y = self._check_iterate(y, arity=self.n_rows, side="apply_transpose")
        outs = []
        for j in range(self.n_cols):
            contributions = []
            for i in range(self.n_rows):
                block = self._blocks[i][j]
                if block is None:
                    continue
                # Guard-narrow (Design C: the runtime check IS the static
                # permission) — mirrors OperatorSum.apply_transpose.
                if not adjointable(block):
                    raise MissingAdjoint(
                        f"CoupledOperator.apply_transpose requires every "
                        f"present block to transpose ((Aᵀ)_ji = (A_ij)ᵀ); "
                        f"block ({i}, {j}) ({type(block).__name__}) has "
                        f"is_adjointable {block.is_adjointable}."
                    )
                contributions.append(block.apply_transpose(y.systems[i]))
            # Non-empty by the all-None-column construction guard.
            outs.append(reduce(lambda a, b: a + b, contributions))
        return CoupledField(systems=tuple(outs))

    @property
    def is_adjointable(self) -> bool:
        r"""``True`` iff every present block is adjointable (the sum-law
        discipline: the transposed grid exists iff every block's transpose
        does)."""
        return all(b.is_adjointable for b in self._present_blocks())

    # ── The solve mode (step 5 — the structure-keyed DIRECT inverse) ───

    def _diagonal_blocks(self) -> "list[LinearOperator] | None":
        r"""The diagonal blocks when the grid is square with a full
        diagonal, else ``None`` (the shared narrow of the solve routes)."""
        if self.n_rows != self.n_cols:
            return None
        diagonals = [self._blocks[i][i] for i in range(self.n_rows)]
        if any(d is None for d in diagonals):
            return None
        return cast("list[LinearOperator]", diagonals)

    def _triangular_orientation(self) -> "str | None":
        r"""``"upper"`` / ``"lower"`` when the grid is block-triangular,
        else ``None``.

        Structural — the ``None``-pattern IS the sparsity: a square grid
        with every diagonal block present and present off-diagonals on at
        most ONE side of the diagonal. A block-diagonal grid is both;
        reported ``"upper"`` (the canonical spelling — the substitution
        degenerates to independent diagonal solves either way).
        """
        if self._diagonal_blocks() is None:
            return None
        n = self.n_rows
        lower = any(
            self._blocks[i][j] is not None
            for i in range(n)
            for j in range(i)
        )
        upper = any(
            self._blocks[i][j] is not None
            for i in range(n)
            for j in range(i + 1, n)
        )
        if lower and upper:
            return None
        return "lower" if lower else "upper"

    def _substitution_ready(self) -> bool:
        r"""``True`` iff the DIRECT substitution route exists: triangular
        AND every diagonal block invertible with a native ``solve`` verb.

        The ``solve`` conjunct keeps the substitution honestly DIRECT: an
        invertible-by-splitting diagonal (an ``OperatorSum`` whose
        ``inverse()`` is the ITERATIVE Green splitting — no ``solve`` verb)
        must not silently smuggle an inner iteration into
        "block-triangular direct"; such a grid falls through to the
        materialized-LU route or refuses.
        """
        if self._triangular_orientation() is None:
            return False
        diagonals = self._diagonal_blocks()
        if diagonals is None:  # unreachable past the orientation check
            return False
        return all(
            d.is_invertible and callable(getattr(d, "solve", None))
            for d in diagonals
        )

    def _transpose_substitution_ready(self) -> bool:
        r"""``True`` iff the TRANSPOSED substitution route exists: the
        forward route plus every coupling block adjointable and every
        diagonal bearing BOTH the direct ``solve_transpose`` verb and an
        adjointable reverse (the #280 two-factor discipline surfaces here
        per block — an LD / multi-D diagonal that defers its reverse scan
        honestly reports ``False``)."""
        if not self._substitution_ready():
            return False
        n = self.n_rows
        for i in range(n):
            for j in range(n):
                block = self._blocks[i][j]
                if block is None:
                    continue
                if i == j:
                    if not callable(getattr(block, "solve_transpose", None)):
                        return False
                    if not block.is_adjointable:
                        return False
                elif not block.is_adjointable:
                    return False
        return True

    def _materialization_route(self) -> bool:
        r"""``True`` iff the dense direct route exists: a square grid whose
        matrix CAN be realized — a full structural assembly, or the typed
        basis probe through the domain's zero exemplar
        (:meth:`CoupledSpace.zeros`). Squareness is arity AND flat
        dimension (a two-sided inverse needs both). Size is deliberately
        NOT consulted — the gate
        (:class:`~orpheus.numerics.operator.MatrixTooLarge`) is a loud
        construction-time failure of :meth:`inverse`, not a structural
        fact."""
        if self.n_rows != self.n_cols:
            return False
        if self._domain.shape != self._codomain.shape:
            return False
        if all(b.is_assemblable for b in self._present_blocks()):
            return True
        return self._domain.zeros_factory is not None

    @property
    def is_invertible(self) -> bool:
        r"""``True`` iff a DIRECT solve route exists (step 5).

        Two structure-keyed routes, in order: **triangular substitution**
        (:meth:`_substitution_ready` — matrix-free, through the diagonal
        blocks' own direct solves) and **materialize/LU**
        (:meth:`_materialization_route` — the
        :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
        EXTRACT). The predicate advertises the ROUTE; exact singularity is
        the LU's loud construction-time :class:`numpy.linalg.LinAlgError`
        (the ``numpy.linalg.solve`` convention — structural predicate,
        value failure loud). The ITERATIVE splitting solve (SI over
        ``A = M − N``) is deliberately NOT an ``is_invertible`` route:
        convergence is spectral (ρ(M⁻¹N) < 1), never structural — the
        drivers own that choice.
        """
        return self._substitution_ready() or self._materialization_route()

    def solve(self, rhs: "CoupledField") -> "CoupledField":
        r"""The structure-keyed DIRECT solve ``A⁻¹·rhs`` (step 5).

        A triangular grid runs the block substitution matrix-free — e.g.
        the ψ½ resolvent ``M = [[LC, Seeding], [None, A_BB]]``: System B's
        march first, then System A's sweep on ``q_A − Seeding·ψ_B``. A
        non-triangular square grid delegates to a FRESH materialized LU
        per call — a ONE-SHOT convenience (hold :meth:`inverse` for
        repeated solves; the factorization is the expensive half).
        """
        rhs = self._check_iterate(rhs, arity=self.n_rows, side="solve")
        orientation = self._triangular_orientation()
        if orientation is not None and self._substitution_ready():
            return self._solve_triangular(rhs, orientation, transpose=False)
        if self._materialization_route():
            return self.inverse().apply(rhs)
        raise NotInvertible(
            f"CoupledOperator.solve: no direct route — the grid is "
            f"neither block-triangular with invertible solve-bearing "
            f"diagonals nor square-materializable (assemblable blocks or "
            f"a domain zeros factory). is_invertible is False."
        )

    def solve_transpose(self, b: "CoupledField") -> "CoupledField":
        r"""The transposed DIRECT solve ``A⁻ᵀ·b``.

        The transpose of a triangular grid is triangular the OTHER way
        (``(Aᵀ)_{ij} = (A_{ji})ᵀ``), so the same one-body substitution
        runs with the visit order flipped; a materialized grid backsolves
        the SAME LU factors under the transpose flag
        (:meth:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator.apply_transpose`).
        Per-block adjoint guards raise
        :class:`~orpheus.numerics.operator.MissingAdjoint` naming the
        offending block.
        """
        b = self._check_iterate(b, arity=self.n_cols, side="solve_transpose")
        orientation = self._triangular_orientation()
        if orientation is not None and self._substitution_ready():
            return self._solve_triangular(b, orientation, transpose=True)
        if self._materialization_route():
            inverse_op = self.inverse()
            return inverse_op.apply_transpose(b)
        raise NotInvertible(
            f"CoupledOperator.solve_transpose: no direct route — the grid "
            f"is neither block-triangular with invertible solve-bearing "
            f"diagonals nor square-materializable. is_invertible is False."
        )

    def _solve_triangular(
        self, rhs: "CoupledField", orientation: str, *, transpose: bool,
    ) -> "CoupledField":
        r"""Block back/forward substitution — ONE body for the four
        orientation × transpose combinations.

        Members are visited in the order that makes every referenced
        partner already solved: upper-triangular forward descends
        (back-substitution — member N−1 first), and each transpose flips
        the direction (upper turns lower under ``(Aᵀ)_{ij} = (A_{ji})ᵀ``).
        Every present coupling block is consumed exactly once through its
        own ``apply`` / ``apply_transpose``; each diagonal inverts its
        member through its native direct verb (``solve`` /
        ``solve_transpose`` — :meth:`_substitution_ready` vouches the
        forward verb exists; the transpose verb is guarded per block
        here). Role algebra stays on the members: the rhs update
        ``rhs_i − A_ij·x_j`` is a member ``−``.
        """
        n = self.n_rows
        descending = (orientation == "upper") != transpose
        order = range(n - 1, -1, -1) if descending else range(n)
        solved: list[Any] = [None] * n
        for i in order:
            acc = rhs.systems[i]
            for j in range(n):
                if j == i:
                    continue
                block = self._blocks[j][i] if transpose else self._blocks[i][j]
                if block is None:
                    continue
                if solved[j] is None:
                    # Unreachable under the triangularity guard — loud,
                    # never a silently-dropped coupling (Cardinal Rule 1).
                    raise RuntimeError(
                        f"CoupledOperator._solve_triangular: substitution "
                        f"ordering bug — the coupling block at grid "
                        f"position {(j, i) if transpose else (i, j)} "
                        f"references the unsolved member {j}."
                    )
                if transpose:
                    if not adjointable(block):
                        raise MissingAdjoint(
                            f"CoupledOperator.solve_transpose: coupling "
                            f"block ({j}, {i}) ({type(block).__name__}) "
                            f"does not transpose (is_adjointable False)."
                        )
                    acc = acc - block.apply_transpose(solved[j])
                else:
                    acc = acc - block.apply(solved[j])
            diagonal = self._blocks[i][i]
            if diagonal is None:  # unreachable past the orientation check
                raise RuntimeError(
                    f"CoupledOperator._solve_triangular: diagonal ({i}, "
                    f"{i}) is absent — the triangularity guard is broken."
                )
            if transpose:
                transpose_verb = getattr(diagonal, "solve_transpose", None)
                if not callable(transpose_verb):
                    raise MissingAdjoint(
                        f"CoupledOperator.solve_transpose: diagonal block "
                        f"({i}, {i}) ({type(diagonal).__name__}) has no "
                        f"solve_transpose — the transposed substitution "
                        f"needs the diagonal's direct transpose-solve."
                    )
                solved[i] = transpose_verb(acc)
            else:
                forward_verb = getattr(diagonal, "solve")
                solved[i] = forward_verb(acc)
        return CoupledField(systems=tuple(solved))

    def inverse(
        self,
    ) -> "CoupledSubstitutionOperator | MatrixInverseOperator":
        r"""The structure-keyed inverse OPERATOR (taxonomy §12/§13).

        * triangular with invertible ``solve``-bearing diagonals →
          :class:`CoupledSubstitutionOperator` (matrix-free block
          substitution — the direct inverse family's coupled sibling);
        * square materializable →
          :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
          (the dense EXTRACT: one LU of :meth:`as_matrix`, loud on exact
          singularity and the size gate);
        * neither → :class:`~orpheus.numerics.operator.NotInvertible`.
        """
        if self._substitution_ready():
            return CoupledSubstitutionOperator(self)
        if self._materialization_route():
            from orpheus.numerics.matrix_inverse_operator import (
                MatrixInverseOperator,
            )

            return MatrixInverseOperator(self)
        raise NotInvertible(
            f"CoupledOperator.inverse: no direct route — the grid is "
            f"neither block-triangular with invertible solve-bearing "
            f"diagonal blocks nor square-materializable. is_invertible "
            f"is False."
        )

    def as_matrix(
        self,
        *,
        basis_shape: tuple[int, ...] | None = None,
        max_dimension: int = 4096,
    ) -> "NDArray":
        r"""Materialize ``[A]`` — the typed-carrier ``Op → Mat`` realization.

        An all-assemblable grid routes through the base assembly
        delegation (:meth:`assemble` → dense, the structural emission).
        Otherwise the base's ndarray basis probe cannot run (a typed block
        matvec refuses flat columns — the base ``as_matrix`` "honest
        scope" note), so the probe runs over TYPED basis elements minted
        from the domain's zero exemplar: column ``j`` is
        ``apply(from_flat(δ_j, domain.zeros()))`` raveled — the same
        C-order column convention as the base. Requires the domain space's
        zeros factory (:meth:`CoupledSpace.zeros`, wired by the instance
        builder); an explicit ``basis_shape`` must match the coupled flat
        dimension, and ``max_dimension`` gates as usual.
        """
        if all(b.is_assemblable for b in self._present_blocks()):
            return super().as_matrix(
                basis_shape=basis_shape, max_dimension=max_dimension,
            )
        n = int(np.prod(self._domain.shape))
        if basis_shape is not None and int(np.prod(basis_shape)) != n:
            raise ValueError(
                f"as_matrix on CoupledOperator: basis_shape {basis_shape} "
                f"(dimension {int(np.prod(basis_shape))}) contradicts the "
                f"coupled domain dimension {n}."
            )
        if n > max_dimension:
            raise MatrixTooLarge(
                f"as_matrix on CoupledOperator: dimension {n} exceeds "
                f"max_dimension={max_dimension}."
            )
        template = self._domain.zeros()
        columns = []
        for j in range(n):
            basis = np.zeros(n)
            basis[j] = 1.0
            image = self.apply(CoupledField.from_flat(basis, template))
            columns.append(np.asarray(image.to_flat(), dtype=float))
        return np.column_stack(columns)

    # ── Assembly (the Op → Mat functor at block offsets) ───────────────

    @property
    def is_assemblable(self) -> bool:
        r"""``True`` iff every present block emits its sparse assembly."""
        return all(b.is_assemblable for b in self._present_blocks())

    def assemble(self) -> "SparseAssembledOperator":
        r"""Scatter each block's assembly at its ``(row_i, col_j)`` offset.

        The block grid → ONE flat sparse matrix over the coupled flat layout
        (:attr:`CoupledSpace.system_slices` — the scoped ``LocalToGlobalMap``):
        block ``(i, j)`` lands at row offset ``Σ_{k<i} n_k``, column offset
        ``Σ_{k<j} n_k``, realized by :func:`scipy.sparse.block_array` (the
        honest carrier — the same 2-P0 ruling that made
        :class:`~orpheus.numerics.assembled_operator.SparseAssembledOperator`
        a thin scipy wrapper). ``None`` blocks stay structural zeros. The
        all-``None`` row/column construction guard doubles as the size-inference
        precondition (``block_array`` derives each block-row/column dimension
        from a present block).
        """
        grid: list[list[Optional["sparse.csr_array"]]] = []
        for i, row in enumerate(self._blocks):
            grid_row: list[Optional["sparse.csr_array"]] = []
            for j, block in enumerate(row):
                if block is None:
                    grid_row.append(None)
                    continue
                # Guard-narrow (Design C) — mirrors OperatorSum.assemble's
                # eager assembly-axis refusal, naming the offending block.
                if not assemblable(block):
                    raise MissingAssembly(
                        f"CoupledOperator.assemble: block ({i}, {j}) "
                        f"({type(block).__name__}) emits no structural "
                        f"assembly (is_assemblable False) — the block "
                        f"scatter needs every present block's (row, col, "
                        f"value) emission; the probing as_matrix fallback "
                        f"is not available on a typed-carrier grid."
                    )
                grid_row.append(block.assemble().matrix)
            grid.append(grid_row)
        return SparseAssembledOperator(
            sparse.block_array(grid),
            domain=self._domain,
            codomain=self._codomain,
        )

    def __repr__(self) -> str:
        present = sum(1 for _ in self._present_blocks())
        return (
            f"<CoupledOperator {self.n_rows}×{self.n_cols} "
            f"blocks={present} domain={self._domain.name!r} "
            f"codomain={self._codomain.name!r}>"
        )


# ───────────────────────────────────────────────────────────────────────
# The substitution inverse (the direct inverse family's coupled sibling)
# ───────────────────────────────────────────────────────────────────────


class CoupledSubstitutionOperator(
    InverseWrapMixin["CoupledOperator"],
    LinearOperator["CoupledField", "CoupledField"],
):
    r"""``A⁻¹`` of a block-triangular :class:`CoupledOperator`, realized by
    block back/forward SUBSTITUTION.

    The inverse family's coupled sibling (taxonomy §12): the wrap-delegate
    back-half — the domain↔codomain swap, ``solve`` = the forward block
    matvec, ``is_invertible`` ``True`` with the object-identity involution
    ``inverse() → inner`` — is inherited from
    :class:`~orpheus.numerics.operator.InverseWrapMixin`. The three
    per-sibling pieces:

    * **ctor guard** — the inner grid must be block-triangular with
      invertible ``solve``-bearing diagonal blocks
      (:meth:`CoupledOperator._substitution_ready`); a full grid belongs to
      :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
      (the materialize/LU route).
    * **apply** — delegates to the grid's own
      :meth:`CoupledOperator.solve`: the substitution body lives on the
      FORWARD's native solve verb (the carve-P4 one-body contract — the
      inverse OBJECT and the realization never twin). ``initial_guess`` is
      accepted and dropped: a substitution over direct diagonal solves is
      an exact direct inverse, nothing to seed.
    * **the adjoint axis** — ``(A⁻¹)ᵀ = (Aᵀ)⁻¹``: :meth:`apply_transpose`
      delegates to :meth:`CoupledOperator.solve_transpose` (the transposed
      substitution), advertised iff every coupling block transposes and
      every diagonal carries the direct ``solve_transpose`` verb
      (:meth:`CoupledOperator._transpose_substitution_ready` — the #280
      two-factor discipline per block).
    """

    def __init__(self, inner: "CoupledOperator") -> None:
        if not inner._substitution_ready():
            raise NotInvertible(
                f"CoupledSubstitutionOperator requires a block-triangular "
                f"grid with invertible solve-bearing diagonal blocks; "
                f"{inner!r} is not (a non-triangular square grid takes the "
                f"materialized-LU route — MatrixInverseOperator)."
            )
        super().__init__(inner)

    def apply(
        self,
        rhs: "CoupledField",
        /,
        *,
        initial_guess: "CoupledField | None" = None,
    ) -> "CoupledField":
        r"""Return ``A⁻¹·rhs`` by block substitution.

        ``initial_guess`` is the inverse family's canonical driver
        signature — an exact direct substitution has nothing to seed, so
        it is accepted and unused.
        """
        del initial_guess  # exact direct substitution — nothing to seed
        return self.inner.solve(rhs)

    @property
    def is_adjointable(self) -> bool:
        r"""``True`` iff the transposed substitution is realizable
        (:meth:`CoupledOperator._transpose_substitution_ready`)."""
        return self.inner._transpose_substitution_ready()

    def apply_transpose(self, b: "CoupledField", /) -> "CoupledField":
        r"""Return ``A⁻ᵀ·b`` — the transposed block substitution."""
        return self.inner.solve_transpose(b)

    def __repr__(self) -> str:
        return f"CoupledSubstitutionOperator({self.inner!r})"
