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
grid arm): a seed-carrying mesh builds the 2×2 — the GL sphere and, since
the Q5.6.3 admission flip, EVERY admitted cylinder, whose folded rule
carries on every μ-level (R12a); a non-carrying mesh builds the 1×1
``[[A_AA]]`` — the Cartesian charts and any curvilinear rule with no
carrying level (cylinder admission refuses that class since Q5.6.3; a
μ = −1-noded sphere rule reaches it, the sphere arm having no admission
gate) — System B's
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
**splitting** ``A = M − N`` (Hackbusch 2016 §11) — ``M`` the
sweepable part inverted every step, ``N`` the lagged coupling gains. Both
are constructed HERE, from the SAME piece objects as the grid, and shipped
together as the frozen :class:`WithinGroupSystem` record
(:func:`build_within_group_system` — since B.2d the ONE construction site
of the within-group decomposition; the former
``orpheus.sn.solver._within_group_triple`` / ``_lagged_gains`` pair retired
into it, which is what dissolved this module's tracked construction twin):

* ``M = [[L+C, +Seeding], [None, march]]`` — since step 5 an HONEST
  upper-triangular :class:`~orpheus.numerics.coupled_system.CoupledOperator`
  grid over the SAME piece objects: ``solve`` is the numerics block
  back-substitution (``ψ_B = march⁻¹ q_B`` — System B first, exactly the
  curvilinear sweep order — then ``ψ_A = LC⁻¹(q_A − Seeding·ψ_B)``, the
  bulk sweep on its ray-DECOUPLED channel), ``apply`` the block matvec,
  ``inverse()`` the
  :class:`~orpheus.numerics.coupled_system.CoupledSubstitutionOperator`.
  The FUSED-walk delegation (the ``CoupledInvertibleOperator`` bridge
  threading the B.2d leaf kwargs, DELETED at 5d — R-5.1/R-5.4) dissolved
  at 5b; the walk's joint legs remain a supported channel for the
  operator-free wrapper entries until the step-6 kwarg retirement.
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

.. warning::

   ⛔ **It is a splitting, NOT a "regular splitting."** ``A = M − N`` was
   described as *regular* here and at three further sites until 2026-08-09
   (#341).  "Regular splitting" is Varga's technical term — ``M`` nonsingular
   with ``M⁻¹ ≥ 0`` **and** ``N ≥ 0`` — and it is load-bearing, because its
   payoff is the comparison theorem that would force ``ρ_GS ≤ ρ_J``
   unconditionally.  The tree does not satisfy it.

   `[M]` at zero leakage the boundary operator's octant action is the
   hypercube ``Q_d``, and the intra-octant DD face-to-face transmission is
   ``Σ = (2/D)·1wᵀ − I`` with ``w_a = 2|μ_a|A_a`` and ``D = Σ_t V + Σ_b w_b``
   — rank-one minus identity.  Every ``v`` with ``wᵀv = 0`` gives
   ``Σv = −v``, so the spectrum is one absorption-damped eigenvalue
   ``1 − 2Σ_t V/D`` plus **d−1 eigenvalues exactly −1**.  Negative
   eigenvalues ⟹ not a non-negative operator ⟹ no regular splitting ⟹ no
   comparison theorem.  And the conclusion it would have licensed —
   ``ρ_GS ≤ ρ_J`` always, i.e. G-S is *never slower* — is measurably false:
   `[M]` boundary G-S takes **1631** sweeps against Jacobi's **838** on the
   d=3 all-reflective box (extents (1,2,3), cells (3,4,5), Σ_t = (0.8,1.6),
   LS4, ``inner_tol`` 1e-13).

   ⚠ Do NOT read that as "d=3 is the slow case".  The sign is **not** a
   function of dimension: on the same fixture family a d=2 box at extents
   (1,2), cells (1,1) measures **202 vs 38 — G-S 5.3× slower**, worse than
   any d=3 reading, while d=2 at (3,4) is 2.5× *faster*.  The mesh alone
   moves it past the d=3 value.  The full treatment, including which claims
   about this survived measurement and which were tombstoned, is
   :ref:`sn-boundary-gs-not-regular` (with the rate regimes at
   :ref:`sn-boundary-gs-rate-regime`).

   The ``−1`` eigenvectors are weighted-zero-average across an octant's
   faces — precisely the undamped face sawtooth (trace ratios
   ``1.074414 / 0.925586``, summing to exactly 2) that ``Σ_t V ψ_c`` cannot
   damp, which is the same object #340 measured as the slow mode costing
   1631 sweeps.  Instruments: ``scratch/probe_341_iteration_spectrum.py``
   (eigen-solves the production iteration matrix),
   ``scratch/probe_inner_budget_law.py``.

   The splitting is still perfectly valid and the fixed point is still
   splitting-invariant.  What is void is only the *guarantee* — which is
   exactly why the schedule choice had to be settled by measurement.

References
==========

* ``.claude/plans/archive/coupled_block_operator_campaign.md`` — the 2×2 posing,
  the B.2 phase rulings, the measured block algebra (``A_sb = 0`` exact,
  ``ρ(M⁻¹N) = 0.371``).
* ``.claude/agent-memory/test-architect/coupled_operator_b2c_builder_verification.md``
  — the B.2c gate spec (P1/P2, the grid≡fused centrepiece, memo F1/F2).
* Hackbusch (2016), *Iterative Solution of Large Sparse Systems*, §11 —
  block partitionings (the step-5 solve consumers).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from orpheus.numerics.coupled_system import CoupledField, CoupledOperator, CoupledSpace
from orpheus.numerics.operator import LinearOperator, SystemRole
from orpheus.sn.operators.boundary import (
    RadialCharacteristicBoundaryOperator,
    SNBoundaryOperator,
)
from orpheus.sn.operators.radial_characteristic import (
    march_start_cosines,
    RadialCharacteristicEmission,
    RadialCharacteristicOperator,
    RadialCharacteristicSeeding,
)
from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.transport.full_field import FullField
from orpheus.transport.operators.multiplication_operator import (
    MultiplicationOperator,
)
from orpheus.transport.operators.n2n import N2NOperator
from orpheus.transport.operators.scattering import ScatteringOperator
from orpheus.transport.radial_characteristic_field import (
    RadialCharacteristicField,
)

if TYPE_CHECKING:
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.sn.operators.streaming import StreamingCollisionOperator
    from orpheus.transport.mesh.material_xs_field import MaterialXSField

__all__ = [
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
        metrics) and, since step 5, the DIRECT solve surface: the full
        2×2 is ``is_invertible`` via the materialize/LU route (the space's
        zero exemplar is wired), the EXTRACT the R5/R11 swap-law gates
        ride — production solves stay the splitting iteration on the
        record's ``implicit_operator``/``explicit_gains``.
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
    grid (the equation), ``implicit_operator``/``explicit_gains`` its splitting
    (Hackbusch 2016 §11 — the drivers iterate ``ψ ← M⁻¹(q + N·ψ)`` / GMRES
    on ``(M − N)·ψ = q``). All four members share the SAME piece objects
    (one ``L+C``, one ``S``, one ``B_a``, one ``B_b``, …) — the single
    construction site :func:`build_within_group_system` is what retired
    the ``_within_group_triple``/``_lagged_gains`` construction twin.

    Parameters
    ----------
    loss : CoupledOperator
        ``A`` — the typed loss grid (2×2 carrying / 1×1 seedless; the
        B.2c presence-structural P2 shape). Since step 5 it is
        DIRECT-invertible via the materialize/LU route
        (``CoupledSpace.zeros`` is wired), the EXTRACT the swap-law
        gates ride — production solves stay the splitting iteration.
    space : CoupledSpace
        The coupled carrier space ``loss`` is typed against (P1
        co-production), carrying the zero-exemplar factory (step 5 —
        the typed-carrier materialization seam).
    implicit_operator : CoupledOperator | StreamingCollisionOperator
        ``M`` — **the sweepable part**, solved IMPLICITLY (inverted) each
        step: on a carrying mesh the HONEST upper-triangular grid
        ``[[LC, Seeding], [None, march]]`` whose
        ``solve`` is the numerics block back-substitution and whose
        ``inverse()`` the
        :class:`~orpheus.numerics.coupled_system.CoupledSubstitutionOperator`
        (step 5 — the fused ``CoupledInvertibleOperator`` bridge
        dissolved at 5b, deleted at 5d); the plain ``(L+C)`` seedless.
    explicit_gains : tuple[LinearOperator, ...]
        ``N`` — the lagged couplings, evaluated EXPLICITLY from the previous
        iterate, that the driver applies each step: ONE
        :class:`~orpheus.numerics.coupled_system.CoupledOperator` gain grid
        ``[[S+B_a, ∅], [Emission, B_b]]`` on a carrying mesh; the
        ``(S, B_a)`` tuple seedless (``B_a`` LAST — the boundary-gain
        convention the G-S schedule arm parses).
    """

    loss: "CoupledOperator"
    space: "CoupledSpace"
    implicit_operator: "CoupledOperator | StreamingCollisionOperator"
    explicit_gains: "tuple[LinearOperator, ...]"


def _zero_full_field(sn_mesh: "SNMesh") -> "FullField":
    r"""The zero System-A state — the coupled space's zero-exemplar member.

    Flux-role by construction (``AngularFlux`` bulk ⊕ ``AngularBoundaryFlux``
    trace — the DOMAIN of the within-group operators), spatial-moment-aware
    (a multi-moment closure's φ̂ tail rides the carrier's
    ``angular_trial_space`` mint, CS4b S5).  Consumed by :meth:`CoupledSpace.zeros` — the
    typed-carrier materialization seam (the loss grid's ``as_matrix``
    basis probe and ``MatrixInverseOperator``'s typed return template),
    which is what flips ``loss.is_invertible`` True (the step-5 EXTRACT
    route the swap-law gate rides).
    """
    from orpheus.transport.fields.angular_boundary_flux import (
        AngularBoundaryFlux,
    )
    from orpheus.transport.fields.angular_flux import AngularFlux

    return FullField(
        interior=AngularFlux.zeros(sn_mesh.angular_trial_space),
        boundary=AngularBoundaryFlux.zeros(sn_mesh.angular_trace),
    )


def _zero_full_field_dual(sn_mesh: "SNMesh") -> "FullField":
    r"""The zero System-A COTANGENT — the dual seam's exemplar member.

    Source-role by construction (``AngularSourceSink`` bulk ⊕
    ``AngularBoundarySourceSink`` trace): under #276 A4 duality typing a
    flux member's cotangent is source-classed, so this is what an
    extension-by-zero transpose extends WITH.  Consumed by
    :meth:`CoupledSpace.dual_zeros` — the cotangent twin of the
    :func:`_zero_full_field` seam.
    """
    from orpheus.transport.source_sinks.angular_boundary_source_sink import (
        AngularBoundarySourceSink,
    )
    from orpheus.transport.source_sinks.angular_source_sink import (
        AngularSourceSink,
    )

    return FullField(
        interior=AngularSourceSink.zeros(sn_mesh.angular_trial_space),
        boundary=AngularBoundarySourceSink.zeros(sn_mesh.angular_trace),
    )


def build_streaming_collision(
    sn_mesh: "SNMesh", mat_xs: "MaterialXSField",
) -> "StreamingCollisionOperator":
    r"""The fused within-group loss factor ``L + C`` — THE one LC spelling.

    ``L`` is the σ-free :class:`~orpheus.sn.operators.streaming.StreamingOperator`
    (#257 S8b) and ``C = M[σ_t]`` the collision diagonal from the typed
    ``total_cross_section_field`` accessor; the ``+`` fuses them to the
    :class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` whose
    ``solve`` is the WDD sweep.

    Single source of truth (Cardinal Rule 2): :func:`build_within_group_system`
    and both :class:`~orpheus.sn.solver.SNSolver` binding legs (``__init__`` +
    ``rebind_cross_sections``) previously spelled this composition
    independently — the third spelling was the L-002 collapse trigger
    (B.2d d3 estate). A σ_t rebind flows through by RE-CALLING this builder
    with the mutated ``mat_xs`` (the field accessor is the live read).
    """
    return StreamingOperator.pose(sn_mesh) + MultiplicationOperator(
        coefficient=mat_xs.total_cross_section_field,
        domain=sn_mesh.full_field_space,
        codomain=sn_mesh.full_field_space,
    )


def build_within_group_system(
    sn_mesh: "SNMesh",
    mat_xs: "MaterialXSField",
    *,
    scattering_op: "ScatteringOperator | None" = None,
    n2n_op: "N2NOperator | None" = None,
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
      :class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` — the
      invertible composite whose ``solve`` is the WDD sweep.
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
      (A,B), ``Emission`` (B,A — consuming the K_iso composed at this
      site: ``S.isotropic_energy + N2N.energy``, the §14.1 grouping), ``A_BB`` (the radial straight-characteristic
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
    n2n_op : N2NOperator, optional
        The already-constructed :math:`(n,2n)` source operator (§14.1 —
        first-class since CS4c step 3; the within-group algebra spells
        ``− S − N₂ₙ`` explicitly). ``None`` constructs fresh.
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
            scattering_order=scattering_order,
            space=full_field_space,
        )
    )
    N2N = (
        n2n_op
        if n2n_op is not None
        else N2NOperator.from_solver_data(
            mat_xs=mat_xs, space=full_field_space,
        )
    )
    # L = pure σ-free streaming; C = M[σ_t] — the ONE LC spelling.
    LC = build_streaming_collision(sn_mesh, mat_xs)
    B_a = SNBoundaryOperator(sn_mesh)
    A_AA = LC - S - N2N - B_a
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
        space = CoupledSpace.from_systems(
            (full_field_space,),
            zeros=lambda: CoupledField(
                systems=(_zero_full_field(sn_mesh),),
            ),
            dual_zeros=lambda: CoupledField(
                systems=(_zero_full_field_dual(sn_mesh),),
            ),
        )
        return WithinGroupSystem(
            loss=CoupledOperator([[A_AA]], domain=space, codomain=space),
            space=space,
            implicit_operator=LC,
            explicit_gains=(S, N2N, B_a),
        )

    # System B's pieces, constructed ONCE and shared between the loss
    # grid, the gain grid, AND the implicit-operator grid (single-sourced objects,
    # three compositions — the step-5 construction-seam collapse: the
    # walk's in-solve engine constructions retired with the fused
    # delegation, so THIS is the one march-construction site).
    A_AB = RadialCharacteristicSeeding(sn_mesh)
    reduced = sn_mesh.reduced
    assert reduced is not None  # carrying ⇒ 1-D ⇒ minted by the ctor; narrowing only
    # K_iso composed HERE (the §14.1 grouping): the P0 energy binding of
    # S's own datum + the (n,2n) energy binding — the solver-side sum
    # that replaced ``S.isotropic_kernel``.
    emission = RadialCharacteristicEmission(
        S.isotropic_energy + N2N.energy,
        field_space=member_space,
        full_field_space=full_field_space,
        angular_bulk_space=sn_mesh.angular_bulk_space,
        angular_trace=sn_mesh.angular_trace,
        quadrature=sn_mesh.quad,
        coord=sn_mesh.coord,
    )
    B_b = RadialCharacteristicBoundaryOperator(
        member_space, sn_mesh.bc["xmax"].law,
    )
    march = RadialCharacteristicOperator(
        member_space,
        mat_xs.total_cross_section_field,
        bulk_space=sn_mesh.bulk_space,
        dr=sn_mesh.axis_widths[0],
        start_cosines=march_start_cosines(
            reduced, sn_mesh.radial_characteristic_levels,
        ),
    )
    A_BB = march - B_b
    space = CoupledSpace.from_systems(
        (full_field_space, member_space),
        zeros=lambda: CoupledField(
            systems=(
                _zero_full_field(sn_mesh),
                RadialCharacteristicField.flux_zeros(sn_mesh.radial_characteristic_field_space),
            ),
        ),
        dual_zeros=lambda: CoupledField(
            systems=(
                _zero_full_field_dual(sn_mesh),
                RadialCharacteristicField.source_zeros(sn_mesh.radial_characteristic_field_space),
            ),
        ),
    )
    # The LOSS grid: (A,B) POSITIVE / (B,A) NEGATED — the sign asymmetry
    # documented in the module docstring (Seeding internalizes its loss
    # sign; Emission is a gain, negated into the loss row).
    loss = CoupledOperator(
        [[A_AA, A_AB], [-emission, A_BB]], domain=space, codomain=space,
    )
    # The GAIN grid N = M − A: all POSITIVE (rhs gains); the (A,B) slot is
    # STRUCTURALLY zero — Seeding lives in M (the (A,B) block below).
    N_AA = S + N2N + B_a
    N_AA.system_role = SystemRole.A  # C-fwd stamp, as on A_AA
    N = CoupledOperator(
        [[N_AA, None], [emission, B_b]], domain=space, codomain=space,
    )
    # M — the sweepable part, an HONEST upper-triangular grid (step 5):
    # ``[[LC, Seeding], [None, march]]``.  Its ``solve`` is the numerics
    # block back-substitution (System B's march first, then the bulk
    # sweep on ``q_A − Seeding·ψ_B`` — the ray-DECOUPLED (L+C) leg), its
    # ``apply`` the block matvec — the fused joint delegation
    # (``CoupledInvertibleOperator``, deleted at 5d) dissolved
    # (R-5.1/R-5.4).
    implicit_operator = CoupledOperator(
        [[LC, A_AB], [None, march]], domain=space, codomain=space,
    )
    return WithinGroupSystem(
        loss=loss,
        space=space,
        implicit_operator=implicit_operator,
        explicit_gains=(N,),
    )
