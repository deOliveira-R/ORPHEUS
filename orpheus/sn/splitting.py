r"""The Strategy's splitting of the within-group loss — ``A = M − N`` as a VALUE.

The consumers campaign's step 2 (R-cc6 (i), ruled 2026-09-13) — the operator
/strategy campaign's P5 exit criterion *"a strategy is a value;
``_select_si_splitting`` and the ``str`` flags retire"*, landed on the
LABELLED-PIECE primitive the campaign's PAIRED-CONSTRUCTION ruling
formalized:

    *"make the ASSIGNMENT the primitive and DERIVE ``M``, ``N``: … every piece
    carries exactly one label — implicit or explicit; ``M = Σ(implicit pieces)``,
    ``N = Σ(explicit pieces)``"* (user, 2026-07-29).

The Problem (the hub, through its posed record
:class:`~orpheus.sn.coupled_system.WithinGroupSystem`) owns the loss
:math:`A` and the bound FACTORS it is the signed sum of; it never chooses
which of them is inverted and which is lagged.  That choice is this module's
:class:`Splitting`: a frozen value minted FROM the factors by ONE labelling
site (:meth:`Splitting.from_schedule`), whose members ``M``
(:attr:`Splitting.implicit`) and ``N`` (:attr:`Splitting.explicit`) are
DERIVED from the labelled terms, and whose law ``M − N = A`` is a property
of the value, certifiable against the Problem it was minted from
(:meth:`Splitting.law_residual`).

The O-3 query contract (ruled at the campaign's §27.6): *the Problem answers
with the operators as POSED — original, unmodified; the Strategy answers with
the operators as it USES them — labelled, split, re-bound, lowered — and every
Strategy value certifies itself against the Problem by an algebraic law.*
``Splitting`` is the first such value.

The terms and their signs
=========================

The within-group loss on the loss-sign convention (B.2c) is the SIGNED sum

.. math::

   A \;=\; \underbrace{(L + C)}_{+} \;-\; S \;-\; N_{2n} \;-\; B_a
   \qquad\bigl(+\;A_{AB}\ \text{(seeding)},\ -\;E\ \text{(emission)},\
   +\;A_{BB}\ \text{(march)},\ -\;B_b\ \text{on a carrying mesh}\bigr),

so a term is an operator TOGETHER with the coefficient it carries in
:math:`A` — :class:`LossTerm` ``(operator, sign)``.  A splitting labels each
term implicit or explicit; then

.. math::

   M \;=\; \sum_{k \in \mathcal I} s_k T_k, \qquad
   N \;=\; -\sum_{k \in \mathcal E} s_k T_k, \qquad A = M - N,

with ``N``'s sign flip making the lagged GAINS positive on the right-hand
side (``ψ ← M⁻¹(q + N·ψ)``).  The two shipped labellings:

* **Jacobi** (every geometry; the only labelling on a carrying mesh):
  :math:`\mathcal I = \{L+C\}` (+ seeding, march), :math:`\mathcal E =
  \{S, N_{2n}, B_a\}` (+ emission, :math:`B_b`) — the whole boundary
  lagged as an external gain.
* **Boundary Gauss-Seidel** (multi-D Cartesian, seedless — RULING P1:
  gradings live on :math:`B_a` and only there): :math:`B_a` splits under
  the octant-group schedule
  (:meth:`~orpheus.sn.operators.boundary.SNBoundaryOperator.split` → the
  named pair :math:`B_{\rm lower} + B_{\rm upper} = B_a`, disjoint rows) and
  the strictly-lower half moves ACROSS the boundary: :math:`\mathcal I =
  \{L+C, B_{\rm lower}\}`, :math:`\mathcal E = \{S, N_{2n}, B_{\rm upper}\}`
  (#226 §17 W2).  The absorbed piece rides the derivation as the reified
  forward :class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`
  ``(L+C) − B_lower`` (the domain's own ``−`` dispatch — its ``solve`` is
  the octant-group forward substitution), which is RULED to retire with
  the partition phase (P5, *"retire, don't rename"*): :attr:`Splitting.implicit`
  is annotated by the CAPABILITY it must carry (an invertible operator),
  never by that class, so the retirement touches no signature here.

The collision gains ``S`` and ``N₂ₙ`` are lagged under BOTH labellings
(the sweep never re-scatters mid-sweep; only the boundary coupling gets
Gauss-Seidel).  Under the boundary-Gauss-Seidel splitting ``M − N = A``
holds BIT-EXACTLY because the split writes disjoint rows (no addition is
reordered) — measured ``0.0`` on every gated configuration; the in-class
ERR-056-family mutation (the lower half kept implicit AND lagged whole)
shifts the law by exactly ``|B_lower·x|`` — the absorbed piece's own
action (``2.778702e+00`` on the 2-D architecture box: a fixture-bound
magnitude, never a contract).  ⚠ It is a splitting, NOT a *regular*
splitting (Varga) — no comparison theorem bounds its rate (#341); the
Jacobi/Gauss-Seidel pair is the Mode-9 discriminator (same fixed point,
different rate).

Where ``M``/``N`` are derived
=============================

On a seedless mesh the derivation is a fold over the labelled terms through
the operator algebra (``+``/``−`` by the term's sign); ``M`` for the Jacobi
labelling is the bare ``L+C`` object itself (identity preserved — the
stage-separation gates compare by ``is``).  On a carrying mesh every term is
placed into the 2×2 System-A ⊕ System-B block grid BY ITS ENDS — the coupled
space's members are the addresses, so a mis-placed coupling is
unconstructable
(:class:`~orpheus.numerics.coupled_system.CoupledOperator` type-checks every
block at construction) — giving the HONEST upper-triangular
``M = [[L+C, Seeding], [∅, march]]`` (its ``solve`` the block
back-substitution: System B's march first, then the bulk sweep on
``q_A − Seeding·ψ_B``) and the gain grid
``N = [[S + N₂ₙ + B_a, ∅], [Emission, B_b]]`` (the (A,B) slot structurally
zero — Seeding lives in ``M``).  Placement by ends is the ONE mechanism the
partition phase generalizes (``A_ij = Rᵢ A Jⱼ``, P5.2 — *"one interface, two
constructors"*); the labelling above is the first constructor.

What is NOT here (forward-looking, by ruling)
=============================================

* the point in :math:`\Lambda` (a shift ``σ``): the α row's absorbed shift
  ``M[Σ_t − σ/v]`` is a further Strategy lowering of the PENCIL's
  ``at(σ)`` (the hub's terminal object) — a shifted factor set is minted
  and labelled exactly like this one;
* the SECOND constructor from a partition (P5), the SCC-derived schedule
  (#324) behind the same :class:`~orpheus.sn.loss_representation.sweep_schedule.SweepSchedule`
  call, the ρ-policy that chooses AMONG values (P6), and the Wielandt
  shift ``−(1/k)·F`` as an explicit term (P4.3) — all consume the
  assignment this value carries, not the products it derives;
* the window (the 2-D moment re-binding of the lagged lifts) — a lowering
  of the EXPLICIT pieces the SI driver performs one by one
  (:func:`orpheus.sn.solver._within_group_si`), which is why they stay
  addressable as pieces rather than one fused operator (P3's gain grid on
  both arms is a derivation from this same primitive).
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property, reduce
from typing import TYPE_CHECKING, Any, cast

import numpy as np

from orpheus.numerics.coupled_system import CoupledField, CoupledOperator
from orpheus.numerics.operator import (
    LinearOperator,
    ScaledOperator,
    SupportsInverse,
    SystemRole,
)

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace
    from orpheus.sn.coupled_system import WithinGroupSystem
    from orpheus.sn.loss_representation.sweep_schedule import SweepSchedule
    from orpheus.sn.mesh.augmented_mesh import SNMesh

__all__ = [
    "LossTerm",
    "Splitting",
    "resolve_schedule",
]


@dataclass(frozen=True)
class LossTerm:
    r"""One term of the within-group loss — an operator and the coefficient
    it carries in :math:`A` (``+1`` for the sweepable structure ``L+C``,
    seeding and march; ``−1`` for every gain).

    The sign is DATA of the term, not of the label: a gain lagged into
    ``N`` keeps its ``−1`` (and ``N`` flips it, so gains are positive on the
    right-hand side); a gain absorbed into ``M`` (``B_lower`` under
    Gauss-Seidel) is subtracted there.  Carrying the sign beside the
    operator — rather than as a ``ScaledOperator(-1, ·)`` wrapper — is what
    keeps the derivation on the domain's own algebra: ``(L+C) − B_lower``
    dispatches to the sweep-invertible scheduled composite, which a wrapped
    ``+ (−B_lower)`` would not.
    """

    operator: "LinearOperator"
    sign: int

    def __post_init__(self) -> None:
        if self.sign not in (1, -1):
            raise ValueError(
                f"LossTerm: the coefficient of a term in the loss is +1 or -1; "
                f"got {self.sign!r} for {type(self.operator).__name__}."
            )

    @property
    def is_gain(self) -> bool:
        """Whether the term is a GAIN (a ``−`` term of the loss)."""
        return self.sign < 0


def resolve_schedule(sn_mesh: "SNMesh", inner_schedule: str) -> "SweepSchedule":
    r"""The ONE site where the entry-level ``inner_schedule`` string becomes
    the schedule object a :class:`Splitting` is labelled by.

    * ``"jacobi"`` — the bare all-octants sweep, every geometry.
    * ``"gauss_seidel"`` — the octant-group order with the mesh's
      specular-reflective faces folded (:func:`~orpheus.sn.loss_representation.sweep_schedule.reflective_faces`)
      on a multi-D CARTESIAN mesh; on 1-D or curvilinear meshes it falls
      back to Jacobi — boundary G-S is a no-op on the scattering-dominated
      1-D regime, the 1-D scan is not a wavefront, and a carrying mesh's
      boundary is the ``B_a + B_b`` composite the schedule split must never
      see (RULING P1).  C5.4 (#225): the gate is the GENUINE condition
      ``is_cartesian and not is_1d`` — the pre-C5.4 ``reduced is None``
      proxy was 2-D-Cartesian by coincidence only.

    The string survives on the entries' signatures until the partition phase
    (P5) retires it; nothing downstream of this function reads it.
    """
    from orpheus.sn.loss_representation.sweep_schedule import (
        SweepSchedule,
        reflective_faces,
    )

    if inner_schedule not in ("jacobi", "gauss_seidel"):
        raise ValueError(
            f"Unknown inner_schedule: {inner_schedule!r}. "
            f"Valid choices are 'gauss_seidel' (boundary G-S, multi-D "
            f"Cartesian — auto-falls-back to Jacobi on 1-D) or 'jacobi'."
        )
    if (
        inner_schedule == "gauss_seidel"
        and sn_mesh.is_cartesian
        and not sn_mesh.is_1d
    ):
        return SweepSchedule.gauss_seidel(
            sn_mesh.ndim, sn_mesh.quad.octants, reflective_faces(sn_mesh),
        )
    return SweepSchedule.jacobi(sn_mesh.ndim, sn_mesh.quad.octants)


@dataclass(frozen=True)
class Splitting:
    r"""A splitting ``A = M − N`` of a posed within-group system, as a VALUE.

    The primitive is the ASSIGNMENT — the loss terms labelled implicit
    (:attr:`implicit_pieces`, summed into ``M``) or explicit
    (:attr:`explicit_pieces`, lagged into ``N``); ``M`` and ``N`` are
    derived (:attr:`implicit`, :attr:`explicit`).  Minted by
    :meth:`from_schedule` — the one labelling site — from the Problem's
    record and the schedule it is labelled by.

    Parameters
    ----------
    system : WithinGroupSystem
        The Problem's posed record the value was minted from — its factors
        are the terms, its ``loss`` the operator the law certifies against,
        its ``space`` the address book for block placement.
    schedule : SweepSchedule
        The octant-group schedule that labelled the boundary
        (:func:`resolve_schedule`); the Strategy datum the value records.
    implicit_pieces : tuple[LossTerm, ...]
        The terms solved IMPLICITLY — inverted every step.  The sweepable
        structure (``L+C``; seeding and march on a carrying mesh) is always
        here; a gain piece may join it (``B_lower`` under Gauss-Seidel).
    explicit_pieces : tuple[LossTerm, ...]
        The terms lagged EXPLICITLY — evaluated on the previous iterate.
    """

    system: "WithinGroupSystem"
    schedule: "SweepSchedule"
    implicit_pieces: "tuple[LossTerm, ...]"
    explicit_pieces: "tuple[LossTerm, ...]"

    def __post_init__(self) -> None:
        # The labelling is a PARTITION of the terms: a piece carries exactly
        # one label (PAIRED CONSTRUCTION).  A piece labelled twice would sit
        # in both M and N and the law M − N = A would read the same defect
        # as a piece DROPPED (measured: both 2.78 on the 2-D box), so the
        # partition is refused at construction rather than diagnosed later.
        implicit_ids = {id(term.operator) for term in self.implicit_pieces}
        twice = [
            type(term.operator).__name__
            for term in self.explicit_pieces
            if id(term.operator) in implicit_ids
        ]
        if twice:
            raise ValueError(
                f"Splitting: a term carries exactly one label — {twice} "
                f"appear in BOTH the implicit and the explicit pieces (M and "
                f"N would each contain it, and M − N = A would read the same "
                f"defect as the term dropped)."
            )

    # ── the one labelling site ─────────────────────────────────────────
    @classmethod
    def from_schedule(
        cls, system: "WithinGroupSystem", schedule: "SweepSchedule",
    ) -> "Splitting":
        r"""Label the record's factors under ``schedule``.

        A sequenced schedule (more than one octant group — the
        boundary-Gauss-Seidel order) splits :math:`B_a` and absorbs its
        strictly-lower half into ``M``; the bare all-octants schedule
        (Jacobi) lags the whole boundary.  A carrying mesh (System B
        present) admits only the Jacobi labelling: its boundary is the
        ``B_a + B_b`` composite no octant schedule may split (RULING P1),
        and :func:`resolve_schedule` never hands it a sequenced schedule —
        a direct caller that does is refused rather than silently
        re-labelled.
        """
        f = system.factors
        LC = LossTerm(f.streaming_collision, +1)
        if schedule.is_sequenced:
            if f.is_coupled:
                raise ValueError(
                    "Splitting.from_schedule: a sequenced (boundary "
                    "Gauss-Seidel) schedule labels a SEEDLESS multi-D "
                    "Cartesian system only — a carrying mesh's boundary is "
                    "the B_a + B_b composite the octant split must not "
                    "touch (RULING P1: gradings live on B_a). Resolve the "
                    "schedule through resolve_schedule(), which falls back "
                    "to Jacobi there."
                )
            parts = f.boundary.split(schedule)
            implicit: tuple[LossTerm, ...] = (LC, LossTerm(parts.lower, -1))
            explicit: tuple[LossTerm, ...] = (
                LossTerm(f.scattering, -1),
                LossTerm(f.n2n, -1),
                LossTerm(parts.upper, -1),
            )
        else:
            implicit = (LC,)
            explicit = (
                LossTerm(f.scattering, -1),
                LossTerm(f.n2n, -1),
                LossTerm(f.boundary, -1),
            )
        rc = f.radial_characteristic
        if rc is not None:
            implicit = implicit + (LossTerm(rc.seeding, +1), LossTerm(rc.march, +1))
            explicit = explicit + (LossTerm(rc.emission, -1), LossTerm(rc.boundary, -1))
        return cls(
            system=system, schedule=schedule,
            implicit_pieces=implicit, explicit_pieces=explicit,
        )

    # ── the derived members ────────────────────────────────────────────
    @cached_property
    def implicit(self) -> "SupportsInverse[Any, Any]":
        r"""``M`` — the implicit part, derived from the implicit terms: the
        signed fold on a seedless mesh (``L+C``, or ``(L+C) − B_lower``), the
        upper-triangular System-A ⊕ System-B grid on a carrying mesh.  An
        INVERTIBLE operator by construction — its ``inverse()`` is the
        driver's resolvent (the sweep / the scheduled forward substitution /
        the block back-substitution)."""
        # ``SupportsInverse`` is a Protocol over a class with MUTABLE typed
        # attributes (``system_role``), which pyright matches invariantly —
        # so the concrete grids/composites, each carrying ``inverse()``,
        # are stated to the checker rather than inferred.
        if self.system.is_coupled:
            return cast(
                "SupportsInverse[Any, Any]", self._grid(self.implicit_pieces, negate=False),
            )
        # The fold opens with the sweepable structure and only ever subtracts
        # a boundary piece through the domain's own dispatch, so the result
        # carries ``inverse()`` by construction — the CAPABILITY the
        # annotation names (never a class: the scheduled composite retires
        # with the partition phase without touching this signature).
        return cast("SupportsInverse[Any, Any]", _signed_fold(self.implicit_pieces))

    @cached_property
    def explicit(self) -> "tuple[LinearOperator, ...]":
        r"""``N`` — the lagged part, as the operators the driver applies each
        step (``rhs = q + Σ N_i·ψ``): the gains POSITIVE (a ``−`` term of the
        loss lagged is a ``+`` gain on the right-hand side; a ``+`` term
        lagged is negated).  One operator per piece on a seedless mesh —
        the window re-binds the lifts among them individually — and ONE
        gain grid ``[[S + N₂ₙ + B_a, ∅], [Emission, B_b]]`` on a carrying
        mesh (the coupled iterate is one member; the grid is the operator
        the pieces compose to)."""
        if self.system.is_coupled:
            return (self._grid(self.explicit_pieces, negate=True),)
        return tuple(_as_gain(term) for term in self.explicit_pieces)

    # ── the law ────────────────────────────────────────────────────────
    def law_residual(self, state: "CoupledField") -> np.ndarray:
        r"""The residual of the splitting law, ``(M − Σ N_i − A)·x``, flat and
        aligned with ``A·x`` — the value's certificate against the Problem
        it was minted from (exactly ``0.0`` on every shipped labelling: the
        pieces are the loss's own terms, and the Gauss-Seidel split writes
        disjoint rows).  The carrier bridge lives here: on a carrying mesh
        ``M`` and ``N`` consume the coupled state; on a seedless mesh they
        are System-A operators and the state's single member is applied
        (the asymmetry the partition phase's gain grid on both arms
        removes)."""
        loss_image = self.system.loss.apply(state).to_flat()
        if self.system.is_coupled:
            image = self.implicit.apply(state).to_flat()
            for gain in self.explicit:
                image = image - gain.apply(state).to_flat()
            return image - loss_image
        member = state.systems[0]
        image = self.implicit.apply(member).to_flat()
        for gain in self.explicit:
            image = image - gain.apply(member).to_flat()
        return image - loss_image

    # ── block placement by ends (the carrying arm) ─────────────────────
    def _grid(self, terms: "tuple[LossTerm, ...]", *, negate: bool) -> "CoupledOperator":
        r"""Assemble the labelled terms into the 2×2 block grid, each placed
        BY ITS ENDS (codomain → row, domain → column against the coupled
        space's members); terms sharing a slot are summed in labelling
        order.  ``negate`` builds the GAIN grid (``N = M − A``: every lagged
        term enters with the sign that makes it a right-hand-side gain)."""
        space = self.system.space
        members = space.systems
        slots: dict[tuple[int, int], LinearOperator] = {}
        order: list[tuple[int, int]] = []
        for term in terms:
            op = _as_gain(term) if negate else _as_loss_term(term)
            key = (_member_index(members, term.operator.codomain, "codomain", term),
                   _member_index(members, term.operator.domain, "domain", term))
            if key in slots:
                slots[key] = slots[key] + op
            else:
                slots[key] = op
                order.append(key)
        n = len(members)
        blocks: list[list[LinearOperator | None]] = [
            [slots.get((i, j)) for j in range(n)] for i in range(n)
        ]
        # C-fwd explicit stamp: System membership is the composition context's
        # fact — the model-generic members' honest None would poison the join.
        head = blocks[0][0]
        if head is not None:
            head.system_role = SystemRole.A
        return CoupledOperator(blocks, domain=space, codomain=space)


def _as_loss_term(term: LossTerm) -> "LinearOperator":
    """The term as it enters ``A`` (and ``M``): the operator, or its negation."""
    return term.operator if term.sign > 0 else ScaledOperator(-1.0, term.operator)


def _as_gain(term: LossTerm) -> "LinearOperator":
    """The term as it enters ``N = M − A`` (a right-hand-side gain): a ``−``
    term of the loss is the operator itself; a ``+`` term is negated."""
    return term.operator if term.sign < 0 else ScaledOperator(-1.0, term.operator)


def _signed_fold(terms: "tuple[LossTerm, ...]") -> "LinearOperator":
    r"""``Σ s_k T_k`` through the operator algebra — the first term starts the
    sum (so a single term IS that operator, identity preserved), every
    further term is added or subtracted per its sign (the domain's own
    ``−`` dispatch: ``(L+C) − B_lower`` is the scheduled composite)."""
    if not terms:
        raise ValueError("Splitting: the implicit part must carry at least one term.")
    first, *rest = terms
    if first.sign < 0:
        raise ValueError(
            "Splitting: the implicit part opens with the sweepable structure "
            f"(a + term); got a gain ({type(first.operator).__name__}) first."
        )
    return reduce(
        lambda acc, t: acc + t.operator if t.sign > 0 else acc - t.operator,
        rest,
        first.operator,
    )


def _member_index(
    members: "tuple[FunctionSpace, ...]", end: "FunctionSpace | None",
    which: str, term: LossTerm,
) -> int:
    """The coupled-space member a term's end names — its block row/column."""
    hits = [i for i, member in enumerate(members) if end is not None and member == end]
    if len(hits) != 1:
        raise ValueError(
            f"Splitting: the {which} of {type(term.operator).__name__} names "
            f"{len(hits)} of the coupled space's {len(members)} members — a "
            f"term is placed by its ends, so every end must be exactly one "
            f"member (unbound ends are a P1 defect, not a placement choice)."
        )
    return hits[0]
