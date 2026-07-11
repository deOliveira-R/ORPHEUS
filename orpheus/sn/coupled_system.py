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
(``sn_mesh.radial_characteristic_composite_space``). The grid and its
:class:`~orpheus.numerics.coupled_system.CoupledSpace` are emitted TOGETHER,
aligned by construction (RULING P1: a mismatched operator/space pairing is
unconstructable — this builder is the ψ½ instance's only constructor). The
machinery itself is semantics-agnostic and N-general
(:mod:`orpheus.numerics.coupled_system`); THIS module is where the ψ½
physics binds to it.

The blocks (loss-sign convention IN the grid)
=============================================

The block matvec ``grid.apply([ψ_A, ψ_B])`` IS the within-group loss action
— the same object the SI/Krylov drivers realize today as
``(L+C)·ψ − S·ψ − A_BA·ψ − B·ψ`` on the augmented 3-block carrier
(:func:`orpheus.sn.solver._within_group_triple` /
:func:`orpheus.sn.solver._lagged_gains`). Signs therefore live IN the block
slots, and they are NOT uniform — the two off-diagonals differ, a trap worth
spelling out:

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

The transitional dead slot (B.2c → B.2d) — a DOCUMENTED hazard
==============================================================

Through B.2c System A's ``full_field_space`` is still the 3-block composite
(interior ⊕ trace ⊕ ψ½): the ray has not yet left ``FullField`` (the
eviction is B.2d, atomic per the explorer's V1). The coupled state therefore
carries ψ_A with a **present-ZERO** ``radial_characteristic`` slot and the
REAL ray state in ψ_B. Under that convention the fused ``(L+C)`` walk's
welded seed-feed arm reads zeros (contributes nothing) and the explicit
``A_AB`` block carries the coupling — the grid row is exact. **Nothing
structural enforces the convention yet**: a LIVE-ray ψ_A double-counts the
seed feed (once welded inside ``A_AA``'s walk, once through ``A_AB``). This
is a known, deliberately-unguarded transitional Pattern-4 hole (memo R3): a
runtime guard would calcify a shape B.2d dissolves structurally (the
eviction makes a live-ray ψ_A unrepresentable). The hazard-witness gate
(``TestCoupledBuilder``, G-c2.6) keeps it visible until then.

Sizing note: the coupled flat dimension counts the dead A-side ψ½ slot AND
System B (the ray twice) until the eviction — honest for Krylov workspace
(every DOF the carriers hold is real memory) but NOT the mathematical DOF
count. B.2d restores the honest count when ``full_field_space`` collapses
to 2-block.

Transient construction twin (retired at B.2d)
=============================================

This builder constructs ``L+C``/``S``/``B_a`` with the SAME spellings as
:func:`orpheus.sn.solver._within_group_triple` — two construction sites for
one decomposition, sanctioned exactly as long as the B.2b gain adapters
live: at B.2d the driver goes block-native (the iterate becomes a
:class:`~orpheus.numerics.coupled_system.CoupledField`), the triple's flat
spelling re-founds on this grid, and the twin collapses. Until then the
grid≡fused centrepiece gate (G-c2.3) pins the two sites to each other.

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

from typing import TYPE_CHECKING

from orpheus.numerics.coupled_system import CoupledOperator, CoupledSpace
from orpheus.numerics.operator import SystemRole
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
from orpheus.transport.operators.multiplication_operator import (
    MultiplicationOperator,
)
from orpheus.transport.operators.scattering import ScatteringOperator

if TYPE_CHECKING:
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.transport.mesh.material_xs_field import MaterialXSField

__all__ = ["build_coupled_system"]


def build_coupled_system(
    sn_mesh: "SNMesh",
    mat_xs: "MaterialXSField",
    *,
    scattering_order: int = 0,
) -> "tuple[CoupledOperator, CoupledSpace]":
    r"""Build the ψ½ coupled block operator and its space, aligned by construction.

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
        (``radial_characteristic_space is not None``).
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
    full_field_space = sn_mesh.full_field_space
    S = ScatteringOperator.from_solver_data(
        mat_xs=mat_xs,
        quadrature=sn_mesh.quad,
        scattering_order=scattering_order,
        full_field_space=full_field_space,
    )
    # The _within_group_triple spellings (the sanctioned transient twin —
    # module docstring): L = σ-free streaming, C = M[σ_t] from the typed
    # field accessor, B_a = the System-A trace boundary.
    LC = StreamingOperator(sn_mesh) + MultiplicationOperator(
        coefficient=mat_xs.total_cross_section_field,
        space=full_field_space,
    )
    A_AA = LC - S - SNBoundaryOperator(sn_mesh)
    # C-fwd explicit stamp: System membership is the composition context's
    # fact — the model-generic members' honest None would poison the join.
    A_AA.system_role = SystemRole.A

    # P2 presence predicate = System B's member space itself (the grid needs
    # exactly it; presence-coextensive with the unified engine carrier, and
    # the read doubles as the Optional narrow).
    member_space = sn_mesh.radial_characteristic_composite_space
    if member_space is None:
        # The non-carrying degenerate: System B does not exist, so the
        # coupling is the 1-system grid (an uncoupled System A).
        space = CoupledSpace.from_systems((full_field_space,))
        return (
            CoupledOperator([[A_AA]], domain=space, codomain=space),
            space,
        )

    # (A,B) POSITIVE / (B,A) NEGATED — the sign asymmetry documented in the
    # module docstring (Seeding internalizes its loss sign; Emission is a
    # gain, negated into the loss row).
    A_AB = RadialCharacteristicSeeding(sn_mesh)
    A_BA = -RadialCharacteristicEmission(sn_mesh, S.isotropic_kernel)
    A_BB = RadialCharacteristicOperator(
        sn_mesh, mat_xs.total_cross_section_field,
    ) - RadialCharacteristicBoundaryOperator(sn_mesh)
    space = CoupledSpace.from_systems((full_field_space, member_space))
    return (
        CoupledOperator(
            [[A_AA, A_AB], [A_BA, A_BB]], domain=space, codomain=space,
        ),
        space,
    )
