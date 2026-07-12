r"""The ψ½ coupled system — ``build_coupled_system``, instance #1 of the block machinery.

The co-producing builder of the campaign's 2×2 coupled block operator (step
4d.2 / Phase B.2c): the within-group augmented SN system posed as a typed
block grid over its two systems,

.. math::

    \begin{bmatrix} A_{AA} & A_{AB} \\ A_{BA} & A_{BB} \end{bmatrix}
    \begin{bmatrix} \psi_A \\ \psi_B \end{bmatrix}
    \;=\;
    \begin{bmatrix} q_A \\ q_B \end{bmatrix},

with System A = the SN transport composite (``sn_mesh.full_field_space``)
and System B = the ψ½ radial-characteristic closure
(``sn_mesh.radial_characteristic_field_space``). The grid and its
:class:`~orpheus.numerics.coupled_system.CoupledSpace` are emitted TOGETHER,
aligned by construction (RULING P1: a mismatched operator/space pairing is
unconstructable — this builder is the ψ½ instance's only constructor). The
machinery itself is semantics-agnostic and N-general
(:mod:`orpheus.numerics.coupled_system`); THIS module is where the ψ½
physics binds to it.

The blocks (loss-sign convention IN the grid)
=============================================

The block matvec ``grid.apply([ψ_A, ψ_B])`` IS the within-group loss action
— the object the SI/Krylov drivers realize as ``M·ψ − N·ψ`` through the
record's splitting (since B.2d; the fused flat spelling
``(L+C)·ψ − S·ψ − A_BA·ψ − B·ψ`` of the retired triple/gain-seam pair was
the same action on the pre-eviction 3-block carrier). Signs live IN the
block slots, and they
are NOT uniform — the two off-diagonals differ, a trap worth spelling out:

* ``(A,A) = L + C − S − B_a`` — System A's self-block: streaming +
  collision minus the bulk scattering gain minus the trace boundary gain.
  Explicitly stamped ``SystemRole.A`` (the C-fwd elegance ruling: its
  model-generic members all carry the honest ``None``, so the join would
  poison to ``None`` — the SYSTEM membership is the composition context's
  fact, stamped here, not derivable from the leaves).
* ``(A,B) = +RadialCharacteristicSeeding`` — POSITIVE: its ``apply``
  already emits the seed's term of the fused bulk residual row
  (``−numer_seed/V`` — the loss sign is INTERNAL to the operator, matching
  the in-sweep placement it wraps).
* ``(B,A) = −RadialCharacteristicEmission`` — NEGATED (a
  :class:`~orpheus.numerics.operator.ScaledOperator`): the emission is a
  GAIN (the driver lags ``+A_BA·ψ`` on the rhs), so the loss row carries
  the minus — the ray equation is ``(A_BB − B_b)·ψ_B − Emission·ψ_A = q_B``.
* ``(B,B) = A_BB − B_b`` — System B's self-block: the radial
  straight-characteristic march minus the ray-corner boundary gain.

Presence is STRUCTURAL (P2, subsuming the step-6 guard retirement for the
grid arm): a seed-carrying mesh (the sphere, R12a) builds the 2×2; a
non-carrying mesh (slab / cylinder) builds the 1×1 ``[[A_AA]]`` — System B's
blocks are never constructed there (their constructors refuse seedless
meshes), so "applying a System-B block on a non-carrying mesh" is not a
runtime branch, it is an object that does not exist.

The eviction end-state (B.2d) — a live-ray ψ_A is unrepresentable
=================================================================

Since the B.2d eviction ``FullField`` is a pure 2-block composite: ψ_A
cannot carry ray state at all, so the B.2c dead-slot convention (and its
double-count hazard — the welded seed feed AND the explicit ``A_AB`` block
both firing on a live-ray ψ_A) dissolved STRUCTURALLY, exactly as ruled
(memo R3: no runtime guard — the type system is the guard). The grid's
``(A,A)`` entry acting on a 2-block member IS the ray-decoupled block
action (the fused surfaces substitute a zero seed internally and discard
the ray rows — bit-identical to the retired dead-slot arithmetic); the
explicit ``A_AB``/``A_BA`` blocks carry ALL the coupling. The coupled flat
dimension is the honest two-system sum (no dead padding — the ERR-053
``restart`` sizing reads the true count).

The named splitting ``A = M − N`` (B.2d — the driver's system record)
=====================================================================

The SI/Krylov drivers do not consume the loss grid raw: they consume its
**regular splitting** ``A = M − N`` (Hackbusch 2016 §11) — ``M`` the
sweepable part inverted every step, ``N`` the lagged coupling gains. Both
are constructed HERE, from the SAME piece objects as the grid, and shipped
together as the frozen :class:`WithinGroupSystem` record
(:func:`build_within_group_system` — since B.2d the ONE construction site
of the within-group decomposition; the former
``orpheus.sn.solver._within_group_triple`` / ``_lagged_gains`` pair retired
into it, which is what dissolved this module's tracked construction twin):

* ``M = [[L+C, +Seeding], [0, A_BB-march]]`` — realized by the FUSED
  ``(L+C)`` walk, which computes exactly this joint block action (the
  welded seed feed IS the (A,B) coupling; the in-walk ψ½ recurrence IS the
  (B,B) march), re-typed onto the coupled carrier by
  :class:`CoupledInvertibleOperator` through the B.2d explicit leaf-kwarg
  legs (``radial_characteristic_source``/``_flux``, ``seed_cot``/``_out``
  — System B's leaves cross into the walk as kwargs, never on the
  composite).
* ``N = M − A = [[S + B_a, ∅], [+Emission, B_b]]`` — ONE
  :class:`~orpheus.numerics.coupled_system.CoupledOperator` gain grid. The
  (A,B) slot is STRUCTURALLY zero (Seeding lives in M), and the signs are
  all POSITIVE here (gains on the rhs: ``rhs = q + N·ψ``) — the loss grid's
  ``−Emission``/``−B_b`` minus signs are the ``M − N`` complement, not a
  contradiction.

On a seedless mesh the record degrades structurally: ``M`` is the plain
``(L+C)`` and ``N`` the ``(S, B_a)`` tuple — the seedless driver paths
(multi-D G-S split, 2-D windowing) consume those bare pieces ZERO-TOUCH
(the B.2d DP-seedless ruling: the coupled carrier appears exactly where
System B exists).

References
==========

* ``.claude/plans/coupled_block_operator_campaign.md`` — the 2×2 posing,
  the B.2 phase rulings, the measured block algebra (``A_sb = 0`` exact,
  ``ρ(M⁻¹N) = 0.371``).
* ``.claude/agent-memory/test-architect/coupled_operator_b2c_builder_verification.md``
  — the B.2c gate spec (P1/P2, the grid≡fused centrepiece, memo F1/F2).
* Hackbusch (2016), *Iterative Solution of Large Sparse Systems*, §11 —
  block partitionings (the step-5 solve consumers).
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING

from orpheus.numerics.coupled_system import CoupledField, CoupledOperator, CoupledSpace
from orpheus.numerics.operator import InverseWrapMixin, LinearOperator, SystemRole
from orpheus.sn.operators.boundary import (
    RadialCharacteristicBoundaryOperator,
    SNBoundaryOperator,
)
from orpheus.sn.operators.radial_characteristic import (
    RadialCharacteristicEmission,
    RadialCharacteristicOperator,
    RadialCharacteristicSeeding,
)
from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.transport.full_field import FullField
from orpheus.transport.operators.multiplication_operator import (
    MultiplicationOperator,
)
from orpheus.transport.operators.scattering import ScatteringOperator
from orpheus.transport.radial_characteristic_field import (
    RadialCharacteristicField,
)

if TYPE_CHECKING:
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.sn.operators.streaming import InvertibleOperator
    from orpheus.transport.mesh.material_xs_field import MaterialXSField

__all__ = [
    "CoupledInvertibleOperator",
    "CoupledSweepOperator",
    "WithinGroupSystem",
    "build_coupled_system",
    "build_streaming_collision",
    "build_within_group_system",
]


def build_coupled_system(
    sn_mesh: "SNMesh",
    mat_xs: "MaterialXSField",
    *,
    scattering_order: int = 0,
) -> "tuple[CoupledOperator, CoupledSpace]":
    r"""Build the ψ½ coupled block operator and its space, aligned by construction.

    Since B.2d this is the loss-grid VIEW of :func:`build_within_group_system`
    (delegation — one construction site): a consumer needing only the loss
    surface ``A`` will take this pair — the PLANNED such consumers are the
    d2 ``evaluate_residual`` re-type, the assembly arm, and the DSA
    substrate (#2); today the campaign gates are the callers. The drivers
    take the full :class:`WithinGroupSystem` record (loss + splitting).

    The co-producing mechanism (P1): the typed grid and the
    :class:`~orpheus.numerics.coupled_system.CoupledSpace` it is typed
    against are emitted together — ``op.domain == op.codomain == space`` and
    every block sits at the grid position its declared spaces admit (a
    mis-placed coupling raises
    :class:`~orpheus.numerics.operator.IncompatibleOperatorComposition` at
    construction). Presence is structural (P2): the returned grid is 2×2 on
    a seed-carrying mesh, 1×1 (``[[A_AA]]``) otherwise. See the module
    docstring for the block table, the sign conventions, and the
    transitional dead-slot hazard.

    Parameters
    ----------
    sn_mesh : SNMesh
        The augmented geometry — supplies both member spaces, the
        quadrature, and the R12a presence predicate
        (``radial_characteristic_field_space is not None``).
    mat_xs : MaterialXSField
        The mesh-materialized macroscopic cross sections (the solver's
        ``sn_mesh.material_xs_field()``): σ_t feeds ``C`` AND ``A_BB`` (one
        typed field object — the mesh-identity invariant holds by
        construction), the scattering table feeds ``S``, whose isotropic
        kernel the emission block shares (single-sourced ``K_iso``).
    scattering_order : int
        Legendre truncation for ``S`` (0 = P0 — the
        :class:`~orpheus.sn.solver.SNSolver` default).

    Returns
    -------
    (CoupledOperator, CoupledSpace)
        The within-group loss grid and its carrier space. The grid carries
        ``apply``/``apply_transpose`` (and ``.H`` via the member-wise
        metrics); the block SOLVE is the campaign's step-5 deliverable
        (``is_invertible`` is ``False`` here).
    """
    system = build_within_group_system(
        sn_mesh, mat_xs, scattering_order=scattering_order,
    )
    return (system.loss, system.space)


# ───────────────────────────────────────────────────────────────────────
# The transitional fused-state bridge (B.2d d1 — retires at d2)
# ───────────────────────────────────────────────────────────────────────


def _system_a_member(state: "CoupledField | FullField") -> "FullField":
    r"""Read System A's member off a driver state — coupled or bare.

    The ONE presence-dispatch reader of the DP-seedless ruling: the solve
    sites hold a :class:`~orpheus.numerics.coupled_system.CoupledField` on
    a carrying mesh and the bare 2-block composite elsewhere (the coupled
    carrier appears exactly where System B exists), so every System-A read
    (the φ reduction, the boundary trace, the finalize threading) funnels
    through this one seam instead of branching at each call site.
    """
    if isinstance(state, CoupledField):
        member = state.systems[0]
        if not isinstance(member, FullField):
            raise TypeError(
                f"_system_a_member: the ψ½ coupled pair carries System A "
                f"(a FullField composite) at position 0; got "
                f"{type(member).__name__}."
            )
        return member
    return state


def _system_b_member(
    state: "CoupledField | FullField",
) -> "RadialCharacteristicField | None":
    r"""Read System B's member off a driver state, ``None`` where it does
    not exist (a bare seedless composite).

    The sibling of :func:`_system_a_member`. A bare :class:`FullField` is
    pure 2-block since the B.2d eviction — a "fused state carrying a live
    ψ½ ray" (the B.2c double-count hazard this reader used to refuse) is
    no longer representable, so the bare arm simply has no System B.
    """
    if isinstance(state, CoupledField):
        if state.n_systems != 2:
            raise ValueError(
                f"_system_b_member: the ψ½ coupled pair has exactly 2 "
                f"systems; got {state.n_systems}."
            )
        member = state.systems[1]
        if not isinstance(member, RadialCharacteristicField):
            raise TypeError(
                f"_system_b_member: the ψ½ coupled pair carries System B "
                f"(a RadialCharacteristicField) at position 1; got "
                f"{type(member).__name__}."
            )
        return member
    return None


def _require_coupled_pair(
    state: "CoupledField", sn_mesh: "SNMesh", context: str,
) -> "tuple[FullField, RadialCharacteristicField]":
    r"""Parse a coupled ψ½ pair ``[ψ_A, ψ_B]`` at an operator boundary.

    The shared parse of :class:`CoupledInvertibleOperator`'s four action
    surfaces (parse-don't-validate): System A through
    :func:`_system_a_member`, System B through the composite's own
    :meth:`~orpheus.transport.radial_characteristic_field.RadialCharacteristicField.require_member`
    (carrier class + mesh identity). A bare System-A composite refuses —
    ``M`` is the JOINT operator; the seedless ``(L+C)`` is its own resolvent
    (DP-seedless).
    """
    psi_a = _system_a_member(state)
    psi_b = _system_b_member(state)
    if psi_b is None:
        raise TypeError(
            f"{context}: expected the coupled pair [ψ_A, ψ_B] — a bare "
            f"System-A composite has no ray member (on a seedless mesh the "
            f"resolvent is the plain (L+C), never this joint M)."
        )
    return psi_a, RadialCharacteristicField.require_member(
        psi_b, mesh=sn_mesh, context=context,
    )


# ───────────────────────────────────────────────────────────────────────
# M — the sweepable part of the splitting, on the coupled carrier
# ───────────────────────────────────────────────────────────────────────


class CoupledInvertibleOperator(LinearOperator["CoupledField", "CoupledField"]):
    r"""``M`` — the sweepable part of the within-group splitting, re-typed
    onto the coupled pair.

    The named ``M`` of ``A = M − N`` (:class:`WithinGroupSystem`): the fused
    ``(L+C)`` walk, which computes the joint block action

    .. math::

        M \;=\; \begin{bmatrix} L+C & +\text{Seeding} \\ 0 & A_{BB}\text{-march} \end{bmatrix}

    — the welded seed feed IS the (A,B) coupling and the in-walk ψ½
    recurrence IS the (B,B) march — presented on the coupled carrier
    through the B.2d EXPLICIT LEAF LEGS: each surface splits the pair and
    hands System B's member DIRECTLY to the fused surface's composite leg
    kwargs (``radial_characteristic_source``/``_flux`` forward,
    ``seed_cot``/``_out`` transposed — the walk marches the split members
    natively since 4e; the ``to_unified``/``from_unified`` bridge and its
    buffers are retired), then packs the walk-filled composite as ψ_B.

    All four action surfaces: :meth:`apply` (the joint forward matvec —
    the Krylov action), :meth:`solve` (the joint WDD sweep — the SI step),
    :meth:`apply_transpose` / :meth:`solve_transpose` (the reverse-scan
    pair, #280, gated by the fused forward's adjointability). ``inverse()``
    returns the :class:`CoupledSweepOperator` wrap-delegate sibling.

    Lifecycle: at step 5 (#41) the block SOLVE re-poses M as an honest
    block-triangular substitution over the grid (the fused delegation
    dissolves). The M-M block split — the fused joint recurrence vs the
    seed-zeroed/bulk-zeroed block sum — is the campaign's intrinsic
    principled-equiv row (B.2c §0: ~5.5e-16), so grid-vs-M comparisons
    carry that bar, never bitwise.
    """

    # M spans both systems by construction (the step-4a lattice).
    system_role = SystemRole.COUPLED

    def __init__(
        self,
        fused: "InvertibleOperator",
        *,
        space: "CoupledSpace",
        sn_mesh: "SNMesh",
    ) -> None:
        #: The fused ``(L+C)`` forward realizing the joint action (read by
        #: the driver-routing gates' wraps-predicate).
        self.fused = fused
        self._space = space
        self._sn_mesh = sn_mesh

    @property
    def domain(self) -> "CoupledSpace":
        return self._space

    @property
    def codomain(self) -> "CoupledSpace":
        return self._space

    # ── The four bridged action surfaces ──────────────────────────────

    def apply(self, x: "CoupledField", /) -> "CoupledField":
        r"""The joint forward matvec ``M·[ψ_A, ψ_B]`` via the fused walk.

        System B's member crosses in DIRECTLY as the
        ``radial_characteristic_flux`` leg (the walk marches the split
        composite natively since 4e); the emitted ray-block rows fill the
        source composite, which IS ψ_B of the output (source-role — the
        matvec output).
        """
        x_a, x_b = _require_coupled_pair(
            x, self._sn_mesh, "CoupledInvertibleOperator.apply",
        )
        rows_buf = RadialCharacteristicField.source_zeros_on(self._sn_mesh)
        y_a = self.fused.apply(
            x_a,
            radial_characteristic_flux=x_b,
            radial_characteristic_source=rows_buf,
        )
        return CoupledField(systems=(y_a, rows_buf))

    def solve(self, rhs: "CoupledField") -> "CoupledField":
        r"""The joint WDD sweep ``M⁻¹·[rhs_A, rhs_B]`` — the SI step.

        rhs_B (the TRUE q½ source) crosses in DIRECTLY as the
        ``radial_characteristic_source`` leg; the walk marches the seed and
        bulk jointly and fills the flux composite in place — ψ_B of the
        output is that marched ray flux, ψ_A the 2-block bulk ⊕ trace.
        """
        rhs_a, rhs_b = _require_coupled_pair(
            rhs, self._sn_mesh, "CoupledInvertibleOperator.solve",
        )
        flux_buf = RadialCharacteristicField.from_mesh(self._sn_mesh)
        psi_a = self.fused.solve(
            rhs_a,
            radial_characteristic_source=rhs_b,
            radial_characteristic_flux=flux_buf,
        )
        return CoupledField(systems=(psi_a, flux_buf))

    @property
    def is_adjointable(self) -> bool:
        r"""The reverse-scan pair exists iff the fused forward's does
        (two-factor geometry gate — DD-1D yes; LD / multi-D defer)."""
        return self.fused.is_adjointable

    def apply_transpose(self, y: "CoupledField", /) -> "CoupledField":
        r"""The joint transposed matvec ``Mᵀ·[y_A, y_B]`` via the fused walk.

        y_B (the cotangent of the forward's ray rows) crosses in DIRECTLY
        as the ``seed_cot`` leg; the domain-side pullback ``Seedingᵀ·y_A +
        D_BBᵀ·y_B`` fills the ``seed_cot_out`` composite, which IS ψ_B of
        the output.
        """
        y_a, y_b = _require_coupled_pair(
            y, self._sn_mesh, "CoupledInvertibleOperator.apply_transpose",
        )
        cot_buf = RadialCharacteristicField.source_zeros_on(self._sn_mesh)
        x_a_bar = self.fused.apply_transpose(
            y_a, seed_cot=y_b, seed_cot_out=cot_buf,
        )
        return CoupledField(systems=(x_a_bar, cot_buf))

    def solve_transpose(self, b: "CoupledField") -> "CoupledField":
        r"""The joint reverse-scan ``M⁻ᵀ·[b_A, b_B]`` (#280 2.5b).

        b_B (the flux-side ray cotangent) crosses in DIRECTLY as the
        ``seed_cot`` leg; the source-side cotangent fills the
        ``seed_cot_out`` composite, which IS ψ_B of the output.
        """
        b_a, b_b = _require_coupled_pair(
            b, self._sn_mesh, "CoupledInvertibleOperator.solve_transpose",
        )
        cot_buf = RadialCharacteristicField.source_zeros_on(self._sn_mesh)
        x_a_bar = self.fused.solve_transpose(
            b_a, seed_cot=b_b, seed_cot_out=cot_buf,
        )
        return CoupledField(systems=(x_a_bar, cot_buf))

    # ── Invertibility (the schedule-triangular family) ────────────────

    @property
    def is_invertible(self) -> bool:
        r"""``True`` — joint walk-triangularity (#284: the fused sweep IS a
        forward substitution) makes M direct-invertible by construction."""
        return True

    def inverse(self) -> "CoupledSweepOperator":
        r"""``M⁻¹`` as an operator — the :class:`CoupledSweepOperator`."""
        return CoupledSweepOperator(self)

    def __repr__(self) -> str:
        return f"CoupledInvertibleOperator({self.fused!r})"


class CoupledSweepOperator(
    InverseWrapMixin["CoupledInvertibleOperator"],
    LinearOperator["CoupledField", "CoupledField"],
):
    r"""``M⁻¹`` — the joint WDD sweep as the inverse operator of
    :class:`CoupledInvertibleOperator` (the coupled sibling of
    :class:`~orpheus.sn.operators.sweep_operator.SweepOperator`).

    The wrap-delegate back-half (domain↔codomain swap, ``solve`` = the
    forward matvec, ``is_invertible``/``inverse() → inner``, the canonical
    seeded ``apply`` signature) is inherited from
    :class:`~orpheus.numerics.operator.InverseWrapMixin`. The ctor guard is
    the :class:`CoupledInvertibleOperator` TYPE itself — joint
    walk-triangularity is what makes the fused forward sweep-invertible.
    """

    def apply(
        self, rhs: "CoupledField", *, initial_guess: "CoupledField | None" = None,
    ) -> "CoupledField":
        r"""Return ``M⁻¹·rhs`` via the joint WDD sweep.

        ``initial_guess`` is the inverse family's canonical seeded keyword
        (#285) — accepted and DROPPED, exactly as
        :meth:`~orpheus.sn.operators.sweep_operator.SweepOperator.apply`
        does: the WDD sweep is a DIRECT exact inverse (the curvilinear ψ½
        starting direction is computed directly from the source since #282
        route (a) 2.5d — nothing to seed).
        """
        del initial_guess  # direct exact inverse — nothing to seed (#282/2.5d)
        return self.inner.solve(rhs)

    @property
    def is_adjointable(self) -> bool:
        r"""``(M⁻¹)ᵀ = (Mᵀ)⁻¹`` exists iff the bridged reverse-scan does."""
        return self.inner.is_adjointable

    def apply_transpose(self, b: "CoupledField", /) -> "CoupledField":
        r"""Return ``M⁻ᵀ·b`` — the bridged reverse-scan transpose-solve."""
        return self.inner.solve_transpose(b)

    def __repr__(self) -> str:
        return f"CoupledSweepOperator({self.inner!r})"


# ───────────────────────────────────────────────────────────────────────
# The system record — the loss and its splitting, from ONE construction
# ───────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class WithinGroupSystem:
    r"""The POSED within-group system: the loss ``A`` and its named
    splitting ``A = M − N``, constructed together.

    The record every within-group solve consumes (eigenvalue SI/Krylov,
    fixed-source SI/Krylov — they differ ONLY in the driver and the
    ``q_ext``, never in this decomposition): ``loss`` is the typed block
    grid (the equation), ``resolvent``/``gains`` its regular splitting
    (Hackbusch 2016 §11 — the drivers iterate ``ψ ← M⁻¹(q + N·ψ)`` / GMRES
    on ``(M − N)·ψ = q``). All four members share the SAME piece objects
    (one ``L+C``, one ``S``, one ``B_a``, one ``B_b``, …) — the single
    construction site :func:`build_within_group_system` is what retired
    the ``_within_group_triple``/``_lagged_gains`` construction twin.

    Parameters
    ----------
    loss : CoupledOperator
        ``A`` — the typed loss grid (2×2 carrying / 1×1 seedless; the
        B.2c presence-structural P2 shape).
    space : CoupledSpace
        The coupled carrier space ``loss`` is typed against (P1
        co-production).
    resolvent : CoupledInvertibleOperator | InvertibleOperator
        ``M`` — the sweepable part: the coupled-bridged fused walk on a
        carrying mesh, the plain ``(L+C)`` seedless.
    gains : tuple[LinearOperator, ...]
        ``N`` — the lagged couplings the driver applies each step: ONE
        :class:`~orpheus.numerics.coupled_system.CoupledOperator` gain grid
        ``[[S+B_a, ∅], [Emission, B_b]]`` on a carrying mesh; the
        ``(S, B_a)`` tuple seedless (``B_a`` LAST — the boundary-gain
        convention the G-S schedule arm parses).
    """

    loss: "CoupledOperator"
    space: "CoupledSpace"
    resolvent: "CoupledInvertibleOperator | InvertibleOperator"
    gains: "tuple[LinearOperator, ...]"


def build_streaming_collision(
    sn_mesh: "SNMesh", mat_xs: "MaterialXSField",
) -> "InvertibleOperator":
    r"""The fused within-group loss factor ``L + C`` — THE one LC spelling.

    ``L`` is the σ-free :class:`~orpheus.sn.operators.streaming.StreamingOperator`
    (#257 S8b) and ``C = M[σ_t]`` the collision diagonal from the typed
    ``total_cross_section_field`` accessor; the ``+`` fuses them to the
    :class:`~orpheus.sn.operators.streaming.InvertibleOperator` whose
    ``solve`` is the WDD sweep.

    Single source of truth (Cardinal Rule 2): :func:`build_within_group_system`
    and both :class:`~orpheus.sn.solver.SNSolver` binding legs (``__init__`` +
    ``rebind_cross_sections``) previously spelled this composition
    independently — the third spelling was the L-002 collapse trigger
    (B.2d d3 estate). A σ_t rebind flows through by RE-CALLING this builder
    with the mutated ``mat_xs`` (the field accessor is the live read).
    """
    return StreamingOperator(sn_mesh) + MultiplicationOperator(
        coefficient=mat_xs.total_cross_section_field,
        space=sn_mesh.full_field_space,
    )


def build_within_group_system(
    sn_mesh: "SNMesh",
    mat_xs: "MaterialXSField",
    *,
    scattering_op: "ScatteringOperator | None" = None,
    scattering_order: int = 0,
) -> "WithinGroupSystem":
    r"""Build the within-group system — loss grid + splitting — from ONE
    piece-construction pass.

    The single source of truth (Cardinal Rule 2) for the within-group
    decomposition every solve consumes. The pieces and their composition
    rules (formerly ``orpheus.sn.solver._within_group_triple``'s contract):

    * ``L`` — the σ-free :class:`~orpheus.sn.operators.streaming.StreamingOperator`
      (#257 S8b: the streaming leaf reads no σ).
    * ``C = M[σ_t]`` — the collision diagonal from the typed
      ``total_cross_section_field`` accessor (#257 S3b / #261); the
      composite ``full_field_space`` lets the ``L + C`` OperatorSum guard
      validate the build. ``L + C`` fuses to the
      :class:`~orpheus.sn.operators.streaming.InvertibleOperator` — the
      resolvent whose ``solve`` is the WDD sweep.
    * ``S`` — the bulk scattering gain (producer-side ``/W`` normalisation
      inside ``S.apply``; no consumer-side rescale). The solver's cached
      instance injects through ``scattering_op`` (a cache seam, NOT a
      configuration flag — ``scattering_order`` is consulted only when
      constructing fresh).
    * ``B_a`` — the System-A trace boundary
      (:class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`),
      a SEPARATE first-class gain (Wave O #208 O.2a): it lives on the
      trace, cannot join the ``L + C`` preconditioner, and carries the
      cosine-weighted ``|Ω·n|·w`` adjoint metric on its trace domain.
      The multi-D Cartesian G-S schedule split lives on ``B_a`` and ONLY
      ``B_a`` (RULING P1 — multi-D Cartesian is seedless by construction,
      so the schedule path always receives the plain operator).
    * System B's blocks (carrying meshes only — R12a): ``Seeding``
      (A,B), ``Emission`` (B,A — sharing ``S.isotropic_kernel``, the
      single-sourced K_iso), ``A_BB`` (the radial straight-characteristic
      march), ``B_b`` (the ray corner). Their constructors refuse seedless
      meshes, so presence is structural (P2).

    Within-group fission is zero (it enters as ``q_ext`` per the
    eigenvalue outer / within-group decomposition) — there is no fission
    piece here.

    The sign table: the LOSS grid carries ``A_AA = L+C−S−B_a``,
    ``+Seeding``, ``−Emission``, ``A_BB−B_b`` (the loss-sign convention,
    B.2c); the GAIN grid ``N = M − A`` carries everything POSITIVE
    (``[[S+B_a, ∅], [Emission, B_b]]`` — gains on the rhs). Both grids'
    (A,A) entries are stamped ``SystemRole.A`` explicitly (the C-fwd
    ruling: the model-generic members' honest ``None`` would poison the
    join).

    Parameters
    ----------
    sn_mesh : SNMesh
        The augmented geometry — supplies both member spaces, the
        quadrature, and the R12a presence predicate.
    mat_xs : MaterialXSField
        The mesh-materialized macroscopic cross sections: σ_t feeds ``C``
        AND ``A_BB`` (one typed field object — mesh-identity by
        construction); the scattering table feeds ``S``.
    scattering_op : ScatteringOperator, optional
        The already-constructed scattering operator (the solver's cached
        instance). ``None`` constructs fresh from ``mat_xs`` at
        ``scattering_order``.
    scattering_order : int
        Legendre truncation for a fresh ``S`` (0 = P0 — the solver
        default). Ignored when ``scattering_op`` is injected.
    """
    full_field_space = sn_mesh.full_field_space
    S = (
        scattering_op
        if scattering_op is not None
        else ScatteringOperator.from_solver_data(
            mat_xs=mat_xs,
            quadrature=sn_mesh.quad,
            scattering_order=scattering_order,
            full_field_space=full_field_space,
        )
    )
    # L = pure σ-free streaming; C = M[σ_t] — the ONE LC spelling.
    LC = build_streaming_collision(sn_mesh, mat_xs)
    B_a = SNBoundaryOperator(sn_mesh)
    A_AA = LC - S - B_a
    # C-fwd explicit stamp: System membership is the composition context's
    # fact — the model-generic members' honest None would poison the join.
    A_AA.system_role = SystemRole.A

    # P2 presence predicate = System B's member space itself (the grid needs
    # exactly it; presence-coextensive with the unified engine carrier, and
    # the read doubles as the Optional narrow).
    member_space = sn_mesh.radial_characteristic_field_space
    if member_space is None:
        # The non-carrying degenerate: System B does not exist — the loss
        # is the 1-system grid and the splitting the bare (L+C, (S, B_a))
        # the seedless driver paths consume zero-touch (DP-seedless).
        space = CoupledSpace.from_systems((full_field_space,))
        return WithinGroupSystem(
            loss=CoupledOperator([[A_AA]], domain=space, codomain=space),
            space=space,
            resolvent=LC,
            gains=(S, B_a),
        )

    # System B's pieces, constructed ONCE and shared between the loss grid
    # and the gain grid (single-sourced objects, two compositions).
    A_AB = RadialCharacteristicSeeding(sn_mesh)
    emission = RadialCharacteristicEmission(sn_mesh, S.isotropic_kernel)
    B_b = RadialCharacteristicBoundaryOperator(sn_mesh)
    A_BB = RadialCharacteristicOperator(
        sn_mesh, mat_xs.total_cross_section_field,
    ) - B_b
    space = CoupledSpace.from_systems((full_field_space, member_space))
    # The LOSS grid: (A,B) POSITIVE / (B,A) NEGATED — the sign asymmetry
    # documented in the module docstring (Seeding internalizes its loss
    # sign; Emission is a gain, negated into the loss row).
    loss = CoupledOperator(
        [[A_AA, A_AB], [-emission, A_BB]], domain=space, codomain=space,
    )
    # The GAIN grid N = M − A: all POSITIVE (rhs gains); the (A,B) slot is
    # STRUCTURALLY zero — Seeding lives in M (the walk's welded feed).
    N_AA = S + B_a
    N_AA.system_role = SystemRole.A  # C-fwd stamp, as on A_AA
    N = CoupledOperator(
        [[N_AA, None], [emission, B_b]], domain=space, codomain=space,
    )
    return WithinGroupSystem(
        loss=loss,
        space=space,
        resolvent=CoupledInvertibleOperator(LC, space=space, sn_mesh=sn_mesh),
        gains=(N,),
    )
