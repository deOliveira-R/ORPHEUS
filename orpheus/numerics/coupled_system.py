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
(``.claude/plans/coupled_block_operator_campaign.md``, Phase B).

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
``solve``      NOT realized here — the block solve (block-triangular
               direct / block-Jacobi / block-G-S) is the campaign's
               step-5 deliverable, keyed by spectral character
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
``+`` it evaluates is a member ``+``, so the affine flux torsor, displacement
minting, and units law of the member family apply unchanged (the same
delegation discipline as ``Composite._map_binary``).

References
----------

* ``.claude/plans/coupled_block_operator_campaign.md`` — the campaign plan:
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
    LinearOperator,
    MissingAdjoint,
    MissingAssembly,
    SystemRole,
    adjointable,
    assemblable,
)
from orpheus.numerics.space import FunctionSpace

if TYPE_CHECKING:
    from numpy.typing import NDArray

__all__ = [
    "CoupledField",
    "CoupledOperator",
    "CoupledSpace",
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
    (not declared) because the member family's AFFINE-TORSOR signatures
    (``flux − flux → displacement``, ``flux + displacement → flux``) fit
    neither :class:`~orpheus.numerics.vector.Vector`'s ``Self + Self → Self``
    spelling nor any simple declaration — the machinery delegates through
    ``Any``-typed lambdas (the ``iteration.py`` posture) without inspecting
    them. Numerics now carries THREE member-contract concepts — ``Vector``
    (the named arithmetic dunders), the ``iteration.py`` ravellable pair
    (ad-hoc private), and this protocol; **collapse trigger**: mint the named
    ravellable-Protocol home (the ``Vector`` promotion precedent for the
    arithmetic half) and converge all three onto it.
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
    member family's role semantics (affine flux torsor, displacement minting,
    units) apply unchanged — this class adds structure, never arithmetic.

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
        truth — pre-checking member types here would block the legitimate
        member torsor (flux member + displacement member).
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
    """

    systems: tuple[FunctionSpace, ...] = field(
        default=(), repr=False, compare=False,
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
        cls, systems: "Sequence[FunctionSpace]",
    ) -> "CoupledSpace":
        r"""Build the coupled space from the member spaces, in system order.

        Derives ``name = "coupled(" + " ⊕ ".join(member names) + ")"`` and
        ``shape = (Σ prod(member.shape),)`` (the flat direct-sum dimension —
        ``prod`` per member so multi-axis member shapes flatten honestly).
        A 1-member coupling is the legitimate degenerate (the uncoupled
        system).
        """
        members = tuple(systems)
        if len(members) == 0:
            raise ValueError(
                "CoupledSpace.from_systems: at least one member space is "
                "required — a 0-system coupled space carries no DOFs."
            )
        name = "coupled(" + " ⊕ ".join(s.name for s in members) + ")"
        total = int(sum(int(np.prod(s.shape)) for s in members))
        return cls(name=name, shape=(total,), systems=members)

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
    Mode-12 note. ``is_invertible`` stays the base ``False``: the block SOLVE
    (block-triangular / block-Jacobi / block-G-S, keyed by spectral character)
    is the campaign's step-5 deliverable, not a free capability.

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
