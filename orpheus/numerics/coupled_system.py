r"""The N-system coupled block machinery: ``CoupledField`` / ``CoupledOperator`` over :class:`~orpheus.numerics.direct_sum.DirectSumSpace`.

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
an emission from landing in the trace slot, and a runtime tag is
not a type constraint. The honest object is a **typed block
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
assembly       each block's emission scattered at its ``(row_i, col_j)``
               block offset into ONE flat matrix
               (:func:`scipy.sparse.block_array`); the blocks recurse
               through the functor
               :func:`~orpheus.numerics.assembled_operator.assemble`
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

**The offsets ARE the block structure.** The assembly needs a local→global
DOF map, and the block law's offsets are those that
:func:`scipy.sparse.block_array` derives from the blocks' shapes: member ``i``'s flat DOFs occupy the
contiguous slice at the cumulative offset of the members before it — the
same layout :meth:`CoupledField.to_flat` packs. They coincide with the coupled space's
:meth:`~orpheus.numerics.direct_sum.DirectSumSpace.slice_of` (the coupled space is a
:class:`~orpheus.numerics.direct_sum.DirectSumSpace` labelled by the coupled builder), because both are cumulative sums of
the same summand sizes in label order; the assembly itself does not read ``slice_of``. This is the **scoped
realization of the deferred ``LocalToGlobalMap``** (the 2-P0 attacker ruling
at ``orpheus/numerics/assembled_operator.py`` deferred a reified map "until a
structured consumer arrives" — the block grid is that structured consumer,
and the map it needs is exactly the member-offset table, nothing more).

The Hilbert adjoint comes FREE (Mode-12 closure by construction)
================================================================

``CoupledOperator`` implements only the Euclidean :meth:`~CoupledOperator.
apply_transpose` — the transposed grid ``(Aᵀ)_{ji} = (A_{ij})ᵀ`` — and
carries coupled-sum domain/codomain whose metric methods dispatch
summand-wise. The metric adjoint ``A.H = G⁺ Aᵀ G`` is then realized ONCE by
the existing :class:`~orpheus.numerics.operator.AdjointOperator` wrapper
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
  assemble≡probe, M3 type-safety, M4 block-``.H`` Mode-12; its M5, the system-role stamp, retired with the tag).
* Trefethen & Bau (1997), *Numerical Linear Algebra*, §1 — the Hilbert
  adjoint :math:`G^{-1}A^{\mathsf T}G` vs. the representation transpose.
* Hackbusch (2016), *Iterative Solution of Large Sparse Systems*, §11 —
  block partitionings and block relaxation (the step-5 solve consumers).
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from functools import reduce
from typing import TYPE_CHECKING, Any, Callable, Optional, Protocol, Sequence, cast, overload

import numpy as np
from scipy import linalg as scipy_linalg, sparse

from orpheus.numerics.assembled_operator import MissingAssembly
from orpheus.numerics.operator import (
    IncompatibleOperatorComposition,
    InverseWrapMixin,
    LinearOperator,
    MatrixTooLarge,
    MissingAdjoint,
    NotInvertible,
    adjointable,
    transpose_invertible,
)
from orpheus.numerics.coproduct_layout import CoproductLayout
from orpheus.numerics.direct_sum import DirectSumSpace, SumCarrier
from orpheus.numerics.space import FunctionSpace

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator

__all__ = [
    "CoupledField",
    "CoupledOperator",
    "CoupledSubstitutionOperator",
    "SystemField",
    "SystemRestrictionOperator",
]


class SystemField(Protocol):  # concept: member | role: interface | kind: element.member
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

    @property
    def space(self) -> "FunctionSpace": ...


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
class CoupledField:  # concept: member | role: weld | kind: element.member | welded: member, record
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
        and must match the paired coupled sum). At least one member;
        a 1-tuple is the legitimate degenerate (an uncoupled system — e.g. a
        non-ray-carrying mesh has no System B, so its "coupled" state is the
        1-system vector, and a System-B slot simply does not exist to be
        filled wrongly).
    labels : tuple[str, ...]
        The summand labels of the coupled sum this field is an element of, one per
        member, REQUIRED (no default: a default label is a positional guess). The
        builder's constants name them; ``x.space == V`` holds because the space is
        derived from the members and these labels. The length check is a raw
        constructor's misuse check, under ``if __debug__:``, and is not a
        guarantee: the formal paths (``V.member(c)``, ``V.zero()``,
        :meth:`with_blocks`, the restriction) cannot produce a mismatch.
    """

    systems: tuple[SystemField, ...]
    labels: tuple[str, ...]

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
        if __debug__:
            if not isinstance(self.labels, tuple) or len(self.labels) != len(self.systems) or len(set(self.labels)) != len(self.labels):
                raise ValueError(
                    f"CoupledField: labels must be a tuple with one distinct label per system; got {self.labels!r} "
                    f"for {len(self.systems)} systems."
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

    @property
    def blocks(self) -> tuple[SystemField, ...]:
        r"""The summand elements in label order: the carrier's access that a direct sum reads (not the operator's ``blocks``)."""
        return self.systems

    def system(self, label: str) -> SystemField:
        r"""The member labelled ``label``; a ``KeyError`` listing the labels otherwise, as :meth:`DirectSumSpace.summand` refuses."""
        try:
            return self.systems[self.labels.index(label)]
        except ValueError:
            raise KeyError(f"{label!r} is not a member label; the labels are {list(self.labels)}") from None

    def with_blocks(self, blocks: "Sequence[SystemField]") -> "CoupledField":
        r"""The same concrete class rebuilt with ``blocks`` as its members, the labels kept."""
        return replace(self, systems=tuple(blocks))

    @property
    def space(self) -> DirectSumSpace:
        r"""The direct-sum space this block vector is an element of — DERIVED
        from the members and the labels, never stored.

        ``DirectSumSpace.of_pairs(zip(labels, member spaces))`` in system order: a
        coupled field IS its member fields, each an element of its own
        composite space (``Composite.space``, CS4b S4 — the F4 ruling), so
        the coupled space is determined by theirs and compares by content
        against the one the instance builder minted (structural identity,
        campaign 1 step 6).  Carries no zero factory: that is builder-side
        metadata (``compare=False``), not part of what the field IS.

        Until step 3 of the consumers campaign (2026-09-17) this member did
        not exist, so every ends check that asked a coupled STATE for its
        space read ``None`` — ``[M]`` ``SourcePosing.__post_init__``'s
        "source lives on the operator's codomain" law was INERT on every SN
        source, all of which are lifted to the one-system coupled state.
        """
        return DirectSumSpace.of_pairs(tuple(zip(self.labels, (member.space for member in self.systems), strict=True)))

    @property
    def principal_bulk_leaf(self) -> "Any | None":
        r"""The leaf whose space norm carries this iterate's convergence
        diagnostics — the PRINCIPAL system's own answer.

        The coupling's first member is the principal field (the same
        system-order convention that fixes the flat layout — for the ψ½
        instance, System A's transport composite), so this delegates to
        ``systems[0]``'s own ``principal_bulk_leaf``. ``None`` when the
        principal member does not expose one (the duck-typed member
        contract — #391 tracks declaring it).
        """
        return getattr(self.systems[0], "principal_bulk_leaf", None)

    # ── Algebra (member-wise delegation; roles live on the members) ────

    def _check_partner(self, other: object) -> "CoupledField":
        r"""Reject a partner that is not a :class:`CoupledField` of the same arity and the same labels.

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
        if other.labels != self.labels:
            raise ValueError(
                f"CoupledField arithmetic requires matching labels; got {list(self.labels)} vs {list(other.labels)}: "
                f"the blocks are paired by label, never by position."
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
        — the same member-offset table :class:`~orpheus.numerics.coproduct_layout.CoproductLayout`
        derives. Every ``n_dof = template.to_flat().size`` sizing site (GMRES
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
        layout = CoproductLayout(tuple((label, (size,)) for label, size in zip(template.labels, sizes)))
        pieces = layout.unpack(flat)
        members = [_member_from_flat(member, pieces[label]) for member, label in zip(template.systems, template.labels)]
        return replace(template, systems=tuple(members))

    def __repr__(self) -> str:
        names = ", ".join(type(m).__name__ for m in self.systems)
        return f"CoupledField[{names}]"


# ───────────────────────────────────────────────────────────────────────
# The block operator
# ───────────────────────────────────────────────────────────────────────


class CoupledOperator(LinearOperator["CoupledField", "CoupledField"]):  # concept: block | role: head | kind: morphism.operator | welded: block, inverse, algorithm, strategy
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
    :class:`~orpheus.numerics.operator.AdjointOperator` over THIS class's
    Euclidean :meth:`apply_transpose` (the transposed grid) and the
    coupled sum's summand-wise metrics — see the module docstring's
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
    domain, codomain : DirectSumSpace
        The coupled member spaces the grid is typed against (REQUIRED — the
        offsets, metrics, and per-block checks all read them). A square
        coupling passes the same space twice; rectangular grids (restriction /
        prolongation stacks) are honest and admitted.
    """

    def __init__(
        self,
        blocks: "Sequence[Sequence[Optional[LinearOperator]]]",
        *,
        domain: DirectSumSpace,
        codomain: DirectSumSpace,
    ) -> None:
        domain_members = tuple(space for _, space in domain.summands)
        codomain_members = tuple(space for _, space in codomain.summands)
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
    def domain(self) -> DirectSumSpace:
        return self._domain

    @property
    def codomain(self) -> DirectSumSpace:
        return self._codomain

    def _present_blocks(self) -> list[LinearOperator]:
        return [b for row in self._blocks for b in row if b is not None]

    # ── The block matvec (apply) and its Euclidean transpose ───────────

    def _check_iterate(self, x: object, *, space: DirectSumSpace, side: str) -> "CoupledField":
        arity = len(space.labels)
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
        space.admit_labels(x.labels, f"CoupledOperator.{side}")
        return x

    def apply(self, x: "CoupledField", /) -> "CoupledField":
        r"""The block matvec :math:`y_i = \sum_j A_{ij}\, x_j`.

        Each present block consumes ITS domain member and emits ITS codomain
        member; the row sum is a member ``+`` (role semantics on the member
        family). ``None`` blocks contribute nothing — structurally, not as a
        computed zero.
        """
        x = self._check_iterate(x, space=self._domain, side="apply")
        outs = []
        for row in self._blocks:
            contributions = [
                block.apply(x.systems[j])
                for j, block in enumerate(row)
                if block is not None
            ]
            # Non-empty by the all-None-row construction guard.
            outs.append(reduce(lambda a, b: a + b, contributions))
        return CoupledField(systems=tuple(outs), labels=self._codomain.labels)

    def apply_transpose(self, y: "CoupledField", /) -> "CoupledField":
        r"""The Euclidean transposed grid :math:`(A^{\mathsf T})_{ji} = (A_{ij})^{\mathsf T}`.

        ``x_j = Σ_i (A_ij)ᵀ y_i`` — the representation transpose only; the
        metric conjugation that makes the Hilbert adjoint lives on
        :class:`~orpheus.numerics.operator.AdjointOperator` +
        :class:`~orpheus.numerics.direct_sum.DirectSumSpace` (see the module docstring). Reachable per the
        :attr:`is_adjointable` advertisement (all present blocks adjointable).
        """
        y = self._check_iterate(y, space=self._codomain, side="apply_transpose")
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
        return CoupledField(systems=tuple(outs), labels=self._domain.labels)

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
        invertible diagonal with no ``solve`` verb (a one-term combination
        :math:`cA`, whose inverse is algebra-closed) must not silently
        smuggle an indirect inverse into "block-triangular direct"; such a
        grid falls through to the materialized-LU route or refuses. A
        diagonal that is a sum of two or more terms is not invertible at all
        (its inverse is the Green of a split the Strategy names,
        numerics_untangling.md ruling 65).
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
        diagonal transpose-invertible
        (:attr:`~orpheus.numerics.operator.LinearOperator.is_transpose_invertible`)
        and adjointable. A ``callable(solve_transpose)`` test, used before,
        passed a multi-D ``(L+C)`` whose reverse scan then raised (`[M]`
        2026-10-18, ``tests/gates/test_adjointability_is_honest.py``)."""
        if not self._substitution_ready():
            return False
        n = self.n_rows
        for i in range(n):
            for j in range(n):
                block = self._blocks[i][j]
                if block is None:
                    continue
                if i == j:
                    if not transpose_invertible(block):
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
        (the sum's ``has_zero`` asks whether it mints one). Squareness is arity AND flat
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
        return self._domain.has_zero

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
        rhs = self._check_iterate(rhs, space=self._codomain, side="solve")
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

    @property
    def is_transpose_invertible(self) -> bool:
        r"""``True`` iff the transposed block substitution runs (:meth:`_transpose_substitution_ready`)."""
        return self._transpose_substitution_ready()

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
        b = self._check_iterate(b, space=self._domain, side="solve_transpose")
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
        return CoupledField(systems=tuple(solved), labels=(self._codomain if transpose else self._domain).labels)

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

        The base's ndarray basis probe cannot run on a block grid (a typed
        block matvec refuses flat columns, the base ``as_matrix`` "honest
        scope" note), so the probe runs over TYPED basis elements minted
        from the domain's zero exemplar: column ``j`` is
        ``apply(from_flat(δ_j, domain.zero()))`` raveled — the same
        C-order column convention as the base. Requires the domain space's
        zero element (``zero()`` of the sum, which refuses when the instance
        builder wired no zero factory); an explicit ``basis_shape`` must match the coupled flat
        dimension, and ``max_dimension`` gates as usual.
        """
        n = self._domain.n_coordinates
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
        template = self._domain.zero()
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

    def assembled_matrix(
        self, emit: "Callable[[LinearOperator], sparse.csr_array]", /,
    ) -> "sparse.sparray | sparse.spmatrix":
        r"""Scatter each block's emission at its ``(row_i, col_j)`` offset.

        The block grid → ONE flat sparse matrix over the coupled flat layout
        (the scoped ``LocalToGlobalMap``, whose offsets coincide with the sum's
        :meth:`~orpheus.numerics.direct_sum.DirectSumSpace.slice_of`):
        block ``(i, j)`` lands at row offset ``Σ_{k<i} n_k``, column offset
        ``Σ_{k<j} n_k``, realized by :func:`scipy.sparse.block_array` (the
        honest carrier — the same 2-P0 ruling that made
        :class:`~orpheus.numerics.assembled_operator.SparseAssembledOperator`
        a thin scipy wrapper). ``None`` blocks stay structural zeros. The
        all-``None`` row/column construction guard doubles as the size-inference
        precondition (``block_array`` derives each block-row/column dimension
        from a present block).
        """
        grid: list[list[Any]] = []
        for i, row in enumerate(self._blocks):
            grid_row: list[Any] = []
            for j, block in enumerate(row):
                if block is None:
                    grid_row.append(None)
                    continue
                if not block.is_assemblable:
                    raise MissingAssembly(
                        f"CoupledOperator's assembly: block ({i}, {j}) "
                        f"({type(block).__name__}) emits no structural "
                        f"assembly (is_assemblable False): the block "
                        f"scatter needs every present block's (row, col, "
                        f"value) emission. The grid's dense matrix is its "
                        f"typed probe, as_matrix."
                    )
                grid_row.append(emit(block))
            grid.append(grid_row)
        return sparse.block_array(grid)

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
):  # concept: inverse | role: weld | kind: morphism.operator | welded: inverse, algorithm
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

    def __repr__(self) -> str:
        return f"CoupledSubstitutionOperator({self.inner!r})"

# ───────────────────────────────────────────────────────────────────────
# The member restriction (the split pair's coupled realization)
# ───────────────────────────────────────────────────────────────────────


class SystemRestrictionOperator(LinearOperator["CoupledField", "CoupledField"]):  # concept: restrict | role: arm | kind: morphism.operator
    r"""Restriction onto ONE member system of a coupling:
    :math:`r_i : \bigoplus_j V_j \to V_i` (as the legitimate 1-member
    degenerate coupling).

    The coupled realization of the tree's ONE restriction concept — a
    **split pair** ``(extension by zero e, restriction r)`` with ``r ∘ e = id`` on
    the member, ``r = e†`` under the direct sum's block-diagonal metric,
    and ``P = e ∘ r`` the orthogonal projector onto the member's copy.
    :class:`~orpheus.numerics.elementary_operators.TraceRestrictionOperator` is the
    index-subset sibling (gather / scatter along an axis); this one
    selects a member SYSTEM of a :class:`CoupledField`. The
    measure-weighted realization of a split pair (the axis marginal: the
    integrating retraction and its section, :mod:`orpheus.numerics.axis_marginal`)
    is a different object. A shared family ABC is deliberately
    NOT minted: no consumer treats "any restriction" generically yet —
    that consumer is the CS4c restriction verb, and the type-minting rule
    defers the base to it.

    The transpose is the **extension-by-zero**
    :math:`\iota_i = r_i^{\mathsf T}`: the input written into member
    ``i`` of a fresh zero of the FULL coupled space, minted through the
    space's own ``dual_zero()`` door (the builder-wired
    COTANGENT-side typed-carrier materialization point — no per-operator
    zero closures; the S4-amendment's un-weld). It follows that the
    operator is usable in a daggered composition only on a space whose
    dual zeros factory is wired — :meth:`apply_transpose` lets the
    space's own loud refusal propagate otherwise.

    Born fully bound (the S4-amendment's base demand): ``domain`` is the
    full coupled sum, ``codomain`` the derived 1-member
    coupling ("a 1-tuple is the legitimate degenerate" —
    :class:`CoupledField`'s ``systems``), whose zeros factory is DERIVED
    from the parent's when wired.

    Why it exists — the first consumer
    ----------------------------------

    The eigen pencil's fission member annihilates the ray system: fission
    is pure bulk, and the ``w = 0`` closed rays never source it. Spelled
    directly as an N×N grid that fact needs an all-``None`` column —
    which the :class:`CoupledOperator` constructor rightly refuses for a
    standalone coupling, and which its ``apply_transpose`` cannot emit
    (the empty-column reduce). The honest algebra is the composition
    ``F_posed = [[F], [E]] ∘ r_bulk``: a rectangular prolongation stack
    (already blessed) after THIS restriction, with no zero blocks
    anywhere and the transpose falling out as
    ``r† ∘ stack† = extension-by-zero ∘ [F†, E†]``.
    """

    def __init__(self, space: DirectSumSpace, *, label: str) -> None:
        member = space.summand(label)
        index = space.labels.index(label)
        self._space = space
        self._label = label
        self._index = index
        derived_zeros: "Callable[[], CoupledField] | None" = None
        derived_dual: "Callable[[], CoupledField] | None" = None
        if space.has_zero:

            def _mint_primal() -> "CoupledField":
                return CoupledField(systems=(space.zero().system(label),), labels=(label,))

            derived_zeros = _mint_primal
        if space.has_dual_zero:

            def _mint_dual() -> "CoupledField":
                return CoupledField(systems=(space.dual_zero().system(label),), labels=(label,))

            derived_dual = _mint_dual
        self._codomain = DirectSumSpace.of_pairs(
            ((label, member),), carrier=SumCarrier(zeros=derived_zeros, dual_zeros=derived_dual),
        )

    @property
    def domain(self) -> DirectSumSpace:
        r"""The full coupled space being restricted."""
        return self._space

    @property
    def codomain(self) -> DirectSumSpace:
        r"""The 1-summand coupling of the selected system, under the parent's label."""
        return self._codomain

    @property
    def label(self) -> str:
        r"""The selected summand's label."""
        return self._label

    def apply(self, x: "CoupledField", /) -> "CoupledField":
        r"""``r_i x`` — select member ``i``, as a 1-member coupling."""
        if not isinstance(x, CoupledField):
            raise TypeError(
                f"SystemRestrictionOperator.apply expects a CoupledField; "
                f"got {type(x).__name__}."
            )
        n = len(self._space.labels)
        if len(x.systems) != n:
            raise ValueError(
                f"SystemRestrictionOperator.apply: iterate couples "
                f"{len(x.systems)} systems but the domain couples {n}."
            )
        self._space.admit_labels(x.labels, f"SystemRestrictionOperator.apply (the system {self._label!r})")
        return CoupledField(systems=(x.system(self._label),), labels=(self._label,))

    def apply_transpose(self, y: "CoupledField", /) -> "CoupledField":
        r"""``\iota_i y`` — the extension-by-zero: member ``i`` carries the
        input, every other member the full space's minted DUAL-role zero.

        Dual, because a transpose is the cotangent map (#276 A4 duality
        typing — it consumes the codomain-cotangent and emits the
        domain-cotangent), so the members it extends WITH are
        cotangent-side zeros: ``dual_zero()``, the door's
        cotangent twin. [M] 2026-08-22, the class question is
        load-bearing: extending with the PRIMAL mint put a flux-classed
        ray zero into the daggered fission chain, and the cross-class
        arithmetic gate refused the sum (the same fact the retired
        per-site ``transpose_zero`` closure encoded via the then-named
        ``source_zeros_on``, today's ``source_zeros``). A primal-side
        member extension, should one
        ever be needed, is a different verb (a future primal-side
        ``extend``), not this transpose."""
        if not isinstance(y, CoupledField):
            raise TypeError(
                f"SystemRestrictionOperator.apply_transpose expects a "
                f"CoupledField; got {type(y).__name__}."
            )
        if len(y.systems) != 1:
            raise ValueError(
                f"SystemRestrictionOperator.apply_transpose: the cotangent "
                f"couples {len(y.systems)} systems but the codomain is the "
                f"1-member coupling."
            )
        self._codomain.admit_labels(y.labels, f"SystemRestrictionOperator.apply_transpose (the system {self._label!r})")
        minted = self._space.dual_zero()
        systems = list(minted.systems)
        systems[self._index] = y.system(self._label)
        return CoupledField(systems=tuple(systems), labels=self._space.labels)

    @property
    def is_adjointable(self) -> bool:
        r"""``True`` — the split pair is exact: ``r† = ι`` (the scatter)."""
        return True
