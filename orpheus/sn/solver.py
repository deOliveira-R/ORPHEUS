"""Unified SN (Discrete Ordinates) eigenvalue solver — operator-algebra form.

Wave E Round 2 (Issue #164): :class:`SNSolver` now constructs the
operator triple :math:`(L, S, F)` at ``__init__`` and consumes the
Wave E Round 1 iteration primitives.  The legacy BiCGSTAB FD-operator
path is replaced by Krylov-on-:meth:`L.apply` with the sweep as
preconditioner — the symmetric closure that closes ERR-026 for
curvilinear geometries.

Inner solver dispatch
=====================

* ``inner_solver="source_iteration"`` (default).  Sweep-driven within-
  group fixed-point iteration.  The closure is the WDD asymmetric
  closure that the curvilinear sweep ships (ERR-026 affected).  This
  path is bit-identical to the Wave A-D source iteration, by
  construction — the loop math is preserved character-for-character so
  the 11 frozen regression snapshots stay green.
* ``inner_solver="krylov"``.  GMRES on the symmetric closure carried
  by the algebraic composition
  :class:`~orpheus.sn.operators.streaming.InvertibleOperator` (= ``L + C``).
  This is the Wave E reconciliation that makes the curvilinear
  ``solve_sn_fixed_source`` path discretization-correct (closes
  ERR-026).  On Cartesian meshes it is bit-identical math to the
  legacy BiCGSTAB FD path.

Boundary conditions default to reflective (infinite lattice) but are
configurable via :class:`~orpheus.geometry.mesh.BC` on the mesh.

.. seealso:: :ref:`theory-discrete-ordinates` — Key Facts, equations, gotchas.
"""

from __future__ import annotations

import time
from collections.abc import Iterable
from dataclasses import replace
from typing import TYPE_CHECKING, Any

import numpy as np
from scipy.sparse.linalg import gmres

from orpheus.data.macro_xs.cell_xs import assemble_cell_xs
from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry import BC, Mesh1D, Mesh2D
from orpheus.numerics.eigenvalue import power_iteration
from orpheus.transport.operators.fission import FissionOperator
from orpheus.transport.reaction_rate_functional import IntegratedReactionRate
from .coupled_system import (
    CoupledInvertibleOperator,
    WithinGroupSystem,
    _system_a_member,
    _system_b_member,
    build_within_group_system,
)
from orpheus.numerics.coupled_system import CoupledField
from orpheus.transport.radial_characteristic_composite import (
    RadialCharacteristicComposite,
)
from .mesh.augmented_mesh import SNMesh
from orpheus.transport.spatial.scheme import DiscretizationSchemeBase
from .spatial.sweep_cache import CollisionCache, GeometryCoefficients
from .operators.streaming import StreamingOperator
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.numerics.moment_layout import (
    AVERAGE_MOMENT,
    cell_moment_count,
    face_moment_tail,
    is_moment_valued_by_flat_rank,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.operators.scattering import ScatteringOperator
from orpheus.transport.mesh.axis import Axis1D
from .loss_representation import transport_sweep
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.radial_characteristic_flux import RadialCharacteristicFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.source_sinks import RadialCharacteristicSourceSink
from orpheus.transport.timed_full_field import TimedFullField

if TYPE_CHECKING:
    # Annotation-only names for the late-imported operator/driver types
    # (their runtime imports stay inside the function bodies — the
    # boundary/iteration modules are one-way late imports here).
    from orpheus.numerics.iteration import KrylovAcceleration, SourceIteration
    from orpheus.numerics.operator import LinearOperator
    from .operators.boundary import SNBoundaryOperator
    from .operators.scheduled_invertible import ScheduledInvertibleOperator
    from .operators.streaming import InvertibleOperator
    from .operators.sweep_operator import SweepOperator
    from .operators.windowing import WindowedSweep


def _apply_default_bcs(
    geometry: "Mesh1D | Mesh2D | tuple[Axis1D, ...]",
    boundary_condition: str,
) -> "Mesh1D | Mesh2D | tuple[Axis1D, ...]":
    """Apply *boundary_condition* string to all faces that lack explicit BCs.

    Returns the original declaration unchanged when it already carries
    ANY explicit :class:`~orpheus.geometry.mesh.BC`, so user-set BCs
    always take precedence over the ``boundary_condition`` parameter.

    C5.5 (#225): handles BOTH entry-surface geometry declarations — a
    legacy :class:`Mesh1D` / :class:`Mesh2D` (per-face dataclass
    fields) and an axis tuple (per-endpoint ``bc`` slots on each
    :class:`~orpheus.transport.mesh.axis.AxisMesh` /
    :class:`~orpheus.transport.mesh.axis.RadialAxisMesh`). The all-or-nothing
    semantics are identical on both representations.
    """
    bc = BC(boundary_condition)
    if isinstance(geometry, Mesh1D):
        if geometry.bc_left is None and geometry.bc_right is None:
            return replace(geometry, bc_left=bc, bc_right=bc)
        return geometry
    if isinstance(geometry, Mesh2D):
        faces = ("bc_xmin", "bc_xmax", "bc_ymin", "bc_ymax")
        if all(getattr(geometry, f) is None for f in faces):
            return replace(geometry, **{f: bc for f in faces})
        return geometry
    axes = tuple(geometry)
    if any(b is not None for ax in axes for b in ax.bc.values()):
        return axes
    return tuple(ax.with_uniform_bc(bc) for ax in axes)


def _as_sn_mesh(
    geometry: "Mesh1D | Mesh2D | tuple[Axis1D, ...]",
    quadrature: "Quadrature",
    materials: "dict[int, Mixture]",
    boundary_condition: "str | None" = None,
    mat_map: "np.ndarray | None" = None,
    *,
    scheme: "DiscretizationSchemeBase | None" = None,
) -> "SNMesh":
    r"""Normalize the entry-surface geometry declaration into an SNMesh.

    The single inbound seam for both ``solve_sn`` entries (C5.5,
    #225): ``geometry`` is a legacy :class:`Mesh1D` / :class:`Mesh2D`
    (the d≤2 user-facing declaration) or an axis tuple — the
    axis-native surface and the ONLY 3-D entry
    (:meth:`SNMesh.from_axes`). ``boundary_condition`` (the
    fixed-source vacuum convention) fills faces only when the
    declaration carries no explicit BC, on either representation;
    ``None`` (the eigenvalue entry) leaves the declaration verbatim —
    unset faces then resolve to the SNMesh-level reflective default
    (the infinite-lattice eigenvalue convention). ``mat_map`` is the
    axes-entry material-assignment channel (shape ``spatial_shape``;
    defaults to single-material id 0) — a legacy mesh carries its own
    and combining the two raises.
    """
    if boundary_condition is not None:
        geometry = _apply_default_bcs(geometry, boundary_condition)
    if isinstance(geometry, (Mesh1D, Mesh2D)):
        if mat_map is not None:
            raise ValueError(
                "mat_map is the axes-entry material channel; a legacy "
                "Mesh1D/Mesh2D carries its own mat_ids/mat_map — "
                "declare the assignment on the mesh."
            )
        return SNMesh(geometry, quadrature, materials, scheme=scheme)
    return SNMesh.from_axes(
        geometry, quadrature, materials, mat_map=mat_map, scheme=scheme,
    )


# Issue #197 PR-TYPED-5: SNFixedSourceResult + SNResult RETIRED.
# Both solver entry points now return a typed :class:`Solution`
# (orpheus.sn.solution); the two legacy data bags collapse into one
# (coding-elegance Pattern 7 — single source of truth) with the
# discrimination living in Solution.is_eigenvalue() / is_fixed_source().
from .solution import IterationHistory, Solution


# ``_zero_within_group_fission`` retired 2026-07-03 (C4, zero callers):
# within-group fission enters as ``q_ext``, so no ``F = 0`` slot exists.


# The former ``_within_group_triple`` (Wave O #208 O.2a) and its gain seam
# ``_lagged_gains`` (campaign step 4c) RETIRED at B.2d d1: the within-group
# decomposition every solve consumes — the loss grid AND its regular
# splitting ``A = M − N`` — is constructed ONCE by
# :func:`orpheus.sn.coupled_system.build_within_group_system` and shipped as
# the :class:`~orpheus.sn.coupled_system.WithinGroupSystem` record (which is
# also what dissolved the builder's tracked construction twin). The triple's
# composition contract (why ``B`` is a separate first-class gain, RULING P1
# on the G-S split, the producer-side ``/W`` note) lives on the builder's
# docstring now. The B.2b transient FullField-gain adapters
# (``_RayEmissionFullFieldGain`` / ``_RayBoundaryFullFieldGain``) retired
# with them — the coupled driver consumes System B's blocks natively.


def _system_a_residual(lhs: "FullField", q_ext: "FullField") -> "FullField":
    r"""System A's typed balance defect — the named bulk/trace ``from_balance`` pair.

    The shared body of :func:`evaluate_residual`'s two arms (bare seedless /
    the coupled pair's A-member): parses the angular family at the composite
    boundary (the #289 F2-sibling role erasure) and mints the 2-block
    residual composite.
    """
    from orpheus.transport.fields._bases import AngularField, AngularBoundaryField
    from orpheus.transport.residuals import AngularResidual, AngularBoundaryResidual

    # Role parse at the composite boundary: ``AngularResidual.from_balance``
    # demands the angular family, but the ``FullField.interior`` slot erases the
    # role (the F2-sibling erasure — #289).
    # A scalar-bulk composite here is a caller error worth raising loudly —
    # on BOTH sides of the balance.
    lhs_bulk = lhs.interior
    q_bulk = q_ext.interior
    if not isinstance(lhs_bulk, AngularField) or not isinstance(
        q_bulk, AngularField
    ):
        raise TypeError(
            f"evaluate_residual: both composites must carry angular-family "
            f"per-ordinate bulks; got lhs {type(lhs_bulk).__name__}, "
            f"q_ext {type(q_bulk).__name__}."
        )
    # Same parse on the trace legs: the widened ``FullField.boundary`` slot
    # (a BoundaryField since #290 P2) erases the family; the SN residual builder
    # demands the ANGULAR trace on both sides.
    lhs_boundary = lhs.boundary
    q_boundary = q_ext.boundary
    if not isinstance(lhs_boundary, AngularBoundaryField) or not isinstance(
        q_boundary, AngularBoundaryField
    ):
        raise TypeError(
            f"evaluate_residual: both composites must carry angular "
            f"(AngularBoundaryField-family) traces; got lhs "
            f"{type(lhs_boundary).__name__}, rhs {type(q_boundary).__name__}."
        )
    # A residual is a one-shot balance defect, not an iterate — it carries no
    # history, so it is the timeless FullField (the history_depth=0 degenerate
    # of TimedFullField; W-C confines the timed type to the driver iterate).
    return FullField(
        interior=AngularResidual.from_balance(lhs=lhs_bulk, rhs=q_bulk),
        boundary=AngularBoundaryResidual.from_balance(
            lhs=lhs_boundary, rhs=q_boundary,
        ),
    )


def evaluate_residual(
    loss_op: "LinearOperator",
    psi: "TimedFullField | CoupledField",
    q_ext: "FullField | CoupledField",
) -> "FullField | CoupledField":
    r"""The typed equation residual :math:`r = A\,\psi - q`, per system.

    Evaluates the within-group balance defect

    .. math::

        r \;=\; A\,\psi \;-\; q

    via the named ``from_balance`` compositions (NOT a bare cross-class
    ``−``, which would mis-type the defect as a source). On a SEEDLESS mesh
    ``A = L + C - S - B`` and the result is the typed 2-block composite
    ``FullField(interior=AngularResidual, boundary=AngularBoundaryResidual)``.
    On a CARRYING mesh (B.2d — the coupled arm) ``loss_op`` is the 2×2
    coupled loss grid, ``psi``/``q_ext`` the coupled pairs, and the result is
    the coupled residual ``CoupledField[r_A, r_B]`` with ``r_B`` the split
    ψ½ pair (:class:`~orpheus.transport.residuals.radial_characteristic_interior_residual.RadialCharacteristicInteriorResidual`
    ⊕ :class:`~orpheus.transport.residuals.radial_characteristic_boundary_residual.RadialCharacteristicBoundaryResidual`)
    — the Mode-12 (b) closure: a System-A-only residual would be
    structurally blind to a wrong seed row, so System B's defect is a typed
    member that cannot be silently dropped.

    This is the **#208 box-7 consumer** of the residual mint — the first
    production-reachable site that types the equation residual (the mint was
    previously unconsumed). It is a diagnostic (``balance_map`` /
    :func:`boundary_vs_interior_split` / ``relative_to``) AND the substrate
    the consistent-DSA low-order correction (`#2`) will consume (``r`` is the
    transport residual the diffusion solve corrects). NOT in the convergence
    path — it is evaluated on a (typically converged) flux, additive.

    Parameters
    ----------
    loss_op : LinearOperator
        The within-group loss operator: the seedless ``L + C - S - B``
        composition, or the carrying 2×2 grid
        (:attr:`~orpheus.sn.coupled_system.WithinGroupSystem.loss`).
        ``loss_op.apply(psi)`` returns source-role members.
    psi : TimedFullField or CoupledField
        The FULL-angular flux state (``bulk`` an ``AngularFlux``; for a
        windowed 2-D solve pass the reconstructed ``Solution.angular_flux``,
        NOT the moment iterate — the operators consume per-ordinate flux).
        The coupled pair on a carrying mesh.
    q_ext : FullField or CoupledField
        The external source (source-role members), matching ``psi``'s
        carrier shape.
    """
    lhs = loss_op.apply(psi)  # A·ψ — source-role members
    if isinstance(psi, CoupledField):
        from orpheus.transport.residuals import (
            RadialCharacteristicBoundaryResidual,
            RadialCharacteristicInteriorResidual,
        )

        if not isinstance(lhs, CoupledField) or not isinstance(
            q_ext, CoupledField
        ):
            raise TypeError(
                "evaluate_residual: a coupled ψ requires the coupled loss "
                "grid and a coupled q_ext — got "
                f"lhs {type(lhs).__name__}, q_ext {type(q_ext).__name__}."
            )
        lhs_a = _system_a_member(lhs)
        q_a = _system_a_member(q_ext)
        lhs_b = _system_b_member(lhs)
        q_b = _system_b_member(q_ext)
        if lhs_b is None or q_b is None:
            raise TypeError(
                "evaluate_residual: the coupled arm requires System-B "
                "members on both the loss output and q_ext."
            )
        r_b = RadialCharacteristicComposite(
            interior=RadialCharacteristicInteriorResidual.from_balance(
                lhs=lhs_b.interior, rhs=q_b.interior,
            ),
            boundary=RadialCharacteristicBoundaryResidual.from_balance(
                lhs=lhs_b.boundary, rhs=q_b.boundary,
            ),
        )
        return CoupledField(
            systems=(_system_a_residual(lhs_a, q_a), r_b),
        )
    if not isinstance(q_ext, FullField):
        raise TypeError(
            f"evaluate_residual: a bare System-A ψ takes a bare FullField "
            f"q_ext; got {type(q_ext).__name__}."
        )
    # The eviction's signature footgun, closed loudly: on a CARRYING mesh a
    # bare System-A call would silently omit r_B — exactly the Mode-12 (b)
    # blindness the split-residual mint exists to prevent (a DSA consumer
    # feeding solution.angular_flux alone would get a residual blind to a
    # wrong seed row). The full-system residual takes the coupled pair.
    bare_mesh = getattr(q_ext.interior, "mesh", None)
    if getattr(bare_mesh, "radial_characteristic_space", None) is not None:
        raise ValueError(
            "evaluate_residual: this mesh carries starting-direction levels "
            "(R12a) — pass the coupled pair [ψ_A, ψ_B] (and the coupled loss "
            "grid) for the FULL-system residual; a bare System-A call would "
            "silently drop System B's defect (Mode-12 (b))."
        )
    return _system_a_residual(lhs, q_ext)


def boundary_vs_interior_split(
    residual: "FullField | CoupledField",
) -> tuple[float, float]:
    r"""Split a typed System-A residual into ``(boundary, interior)`` L2 norms.

    Returns the flat-L2 norm of the boundary residual and of the interior
    (bulk) residual, so :math:`\sqrt{b^2 + i^2} = \lVert r_A\rVert` (System
    A's composite flat norm — EXACT since B.2d made ``FullField`` 2-block;
    the pre-eviction 3-block silently excluded the seed rows from this
    identity, the closed diagnostic gap). Discriminates a BC-realizer /
    reflective-trace defect (large ``boundary``) from an interior-streaming
    defect (large ``interior``). A coupled residual funnels to its System-A
    member; System B's defect is its own typed member — read
    ``residual.systems[1]`` directly (its flat norm is the ray-row defect).
    """
    member = _system_a_member(residual)
    interior = float(np.linalg.norm(np.asarray(member.interior.values).ravel()))
    boundary = float(np.linalg.norm(np.asarray(member.boundary.values).ravel()))
    return boundary, interior


def _within_group_krylov(
    LC: "LinearOperator", *gains: "LinearOperator",
    n_dof: int, max_iter: int, tol: float,
) -> "KrylovAcceleration[Any]":
    r"""GMRES driver on the within-group system ``(M − N)·ψ = q``.

    Single source of truth (Cardinal Rule 2 / Phase 1 R2) for the
    :class:`~orpheus.numerics.iteration.KrylovAcceleration` construction shared
    by the eigenvalue (:meth:`SNSolver._solve_krylov`) and fixed-source
    (:func:`_solve_fixed_source_krylov`) Krylov paths — they previously built
    byte-identical instances.  The call sites pass the
    :class:`~orpheus.sn.coupled_system.WithinGroupSystem` record's splitting:
    ``LC`` is ``M`` (the fused ``(L+C)`` seedless; the
    :class:`~orpheus.sn.coupled_system.CoupledInvertibleOperator` on a
    carrying mesh) and ``*gains`` is ``N`` (the ``(S, B_a)`` pair seedless;
    the ONE coupled gain grid carrying); the matvec is
    ``M.apply − Σ Nᵢ.apply``.

    GMRES is UNPRECONDITIONED (explicit identity) per `issue #200
    <https://github.com/deOliveira-R/ORPHEUS/issues/200>`_ (the block-inverse
    face preconditioner re-enablement).  ``restart`` is sized to the FULL
    problem ``n_dof = N·ng·nx·ny`` — the legacy ``min(50, …)`` clamp left GMRES
    structurally truncated on any mesh with ``n_dof > 50`` (ERR-053).
    """
    from orpheus.numerics.iteration import KrylovAcceleration

    return KrylovAcceleration(
        LC, *gains,
        preconditioner=lambda q: q,  # explicit identity — issue #200 tracks re-enable
        tol=tol, max_iter=max_iter,
        restart=n_dof,
    )


# The former ``_GaussSeidelResolvent`` (Phase 3 sub-step 3c) dissolved in
# #226 step 2 (§17 W2): it paired ``apply = (L+C)ψ`` with
# ``solve = (L+C−B_lower)⁻¹`` — inverses of DIFFERENT operators.  The G-S
# splitting matrix is now REIFIED as
# :class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`
# (``M = (L + C) - B_lower``, an honest forward whose ``solve`` is the
# scheduled forward substitution), and the lagged complement ``B_upper``
# rides the SI driver as an ordinary external gain — structurally congruent
# with the Jacobi path's ``(S, B)``.  See :func:`_select_si_resolvent`.


# The former ``_MomentWindowedResolvent`` (a ``.solve``-shaped driver
# adapter over the windowed product) dissolved at #226 taxonomy step 3:
# :class:`~orpheus.numerics.iteration.SourceIteration` now consumes the
# inverse-application operator DIRECTLY, so :func:`_maybe_window` returns
# the typed composition ``P @ A.inverse()``
# (:class:`~orpheus.sn.operators.windowing.WindowedSweep`) itself.  The
# adapter's accepted-and-dropped ``initial_guess`` note lives on
# :meth:`WindowedSweep.apply`, where the contract actually is.


def _maybe_window(
    sweep: "SweepOperator", scattering_op: "ScatteringOperator",
    sn_mesh: "SNMesh",
) -> "tuple[WindowedSweep | SweepOperator, bool]":
    r"""Phase 5a — compose the 2-D Cartesian angular-windowing product over
    ``sweep`` (the inverse operator ``A.inverse()``), else passthrough.
    Returns ``(step, windowed)``.

    The SINGLE site of the windowing-eligibility gate AND the factory of the
    windowed product (``coding-elegance`` Pattern 7 — the convention lives in
    one place, shared by the eigenvalue and fixed-source SI drivers):
    genuinely 2-D Cartesian holds the SI iterate as harmonic moments via the
    typed composition ``P @ A.inverse()`` (#226 §17 W1) — the ``@`` dispatch
    on a :class:`~orpheus.sn.operators.sweep_operator.SweepOperator` right
    factor fuses to :class:`~orpheus.sn.operators.windowing.WindowedSweep`;
    curvilinear (1-D) stays full-angular — the Morel–Montry Carlson seed
    reads the previous per-ordinate iterate at ``μ=−1`` (lesson L21), which
    the moment tensor does not carry.  ``P`` is sourced from the scattering
    operator's own frame ⇒ the stored moments match ``S``'s internal
    projection term-for-term.

    C5.4 (#225, vv Mode 9): the gate is the GENUINE condition
    ``is_cartesian and ndim == 2`` — the pre-C5.4 ``reduced is None``
    proxy was a coincidence that is ALSO true at d=3 Cartesian and would
    have silently moment-windowed a 3-D solve (the in-sweep moment
    emission is a 2-D kernel; ``FullFieldWavefront`` refuses moment mode).
    """
    if sn_mesh.is_cartesian and sn_mesh.ndim == 2:
        from .operators.windowing import BulkAnalysisOperator

        return (
            BulkAnalysisOperator(scattering_op.frame, sn_mesh) @ sweep,
            True,
        )
    return sweep, False


def _windowed_cold_start(scattering_op, sn_mesh, *, history_depth):
    r"""Zero windowed (moment-bulk) SI cold-start iterate.

    The moment representation the windowed resolvent emits and the
    moment-consuming ``S.apply`` / ``B.apply`` expect — shared by both SI
    drivers (``coding-elegance`` Pattern 2).
    """
    from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
    from orpheus.transport.fields.harmonic_moment_flux import (
        HarmonicMomentFlux,
    )
    from orpheus.transport.timed_full_field import TimedFullField

    return TimedFullField(
        interior=HarmonicMomentFlux.zeros_for_mesh_and_L(
            sn_mesh, scattering_op.scattering_order,
            spatial_moments=sn_mesh.scheme.spatial_basis_per_axis,
        ),
        boundary=AngularBoundaryFlux.zeros_on(sn_mesh),
        _history=(),
        history_depth=history_depth,
    )


def _unwindowed_cold_start(sn_mesh, *, history_depth):
    r"""Zero un-windowed (full-angular) SI cold-start iterate.

    The full-angular ``AngularFlux`` iterate the 1-D / curvilinear SI driver
    holds.  Selects the SpatialMomentSpace factor (#240 D5b-S3) so a
    multi-moment closure (LD) carries the φ̂ axis end-to-end — composing with the
    moment-carrying ``q_ext`` + ``S.apply(ψ)`` in the SI rhs.  DD/Step
    (per_axis == 1) → no factor (byte-identical to the prior ``TimedFullField.zeros``).
    The un-windowed sibling of :func:`_windowed_cold_start` (Pattern 2).
    On a seed-carrying mesh the COUPLED cold start wraps this System-A
    frame with a zero ψ_B via :func:`_coupled_flux_state` (B.2d)."""
    from orpheus.transport.fields.angular_flux import AngularFlux
    from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
    from orpheus.transport.timed_full_field import TimedFullField

    return TimedFullField(
        interior=AngularFlux.zeros_on(
            sn_mesh, spatial_moments=sn_mesh.scheme.spatial_basis_per_axis,
        ),
        boundary=AngularBoundaryFlux.zeros_on(sn_mesh),
        _history=(),
        history_depth=history_depth,
    )


def _radial_characteristic_zeros(sn_mesh) -> "RadialCharacteristicFlux | None":
    r"""The mesh-keyed zero UNIFIED ψ½ flux buffer (``None`` on non-carrying meshes).

    The walk-buffer factory of #282 route (a): the final-reconstruction
    ``transport_sweep`` fills this carrier in place (the walk internals
    stay unified through Phase C/4e — the demote ruling). Presence is
    decided by the ONE R12a source (``SNMesh.radial_characteristic_space``),
    never by the call site. Driver STATES no longer ride it — the coupled
    pair's ψ_B is System B's own composite (B.2d,
    :func:`_coupled_flux_state`)."""
    if sn_mesh.radial_characteristic_space is None:
        return None
    return RadialCharacteristicFlux.zeros_on(sn_mesh)


def _radial_characteristic_source_from_per_ordinate(
    per_ordinate_values: "np.ndarray", sn_mesh,
) -> "RadialCharacteristicSourceSink | None":
    r"""Fold a PER-ORDINATE source to its q½ seed —
    :meth:`RadialCharacteristicSourceSink.from_angular_source` (Legendre-project
    the per-ordinate source to moments, then fold at :math:`\mu=\pm 1`).

    The PER-ORDINATE typed entry to the one fold kernel
    (:func:`~orpheus.numerics.spaces.radial_characteristic_space.fold_moments_to_radial_characteristic`),
    used where the source is genuinely per-ordinate (possibly anisotropic): the
    final total-source reconstruction, the fixed-source rhs external source, and
    the operator-free ``transport_sweep``. The eigenvalue FISSION q½ seed uses the
    MOMENTS entry :func:`_radial_characteristic_fission_seed` instead (its ℓ=0
    emission is already a moment, so the direct fold skips this factory's
    per-ordinate round-trip). Both bottom out in the SAME kernel — no twin, a
    different typed input.
    """
    return RadialCharacteristicSourceSink.from_angular_source(
        per_ordinate_values, sn_mesh,
    )


def _radial_characteristic_fission_seed(
    fission_source: "np.ndarray", sn_mesh,
) -> "RadialCharacteristicSourceSink | None":
    r"""The ψ½ FISSION ray seed — the direct ℓ=0 moments-fold of the fission
    emission (``A_BA_fission = Fold ∘ F.kernel``, factored).

    The eigenvalue outer loop computes ``fission_source = χνΣf·φ/k`` from the
    SCALAR flux (:meth:`SNSolver.compute_fission_source`), so ``F.kernel ∘
    integrate`` is already applied and only the fold remains: the isotropic ℓ=0
    emission is reconstructed at the closed :math:`\mu=\pm 1` rays through the
    migrated
    :class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicReconstruction`
    (the single fold operator). ``None`` on a seedless mesh (no ψ½ ray).

    REPLACES the ``from_isotropic → from_angular_source`` round-trip (campaign
    step 4c commit 2): the old path broadcast ``fission_source`` to a per-ordinate
    iso source then re-projected it back to the ℓ=0 moment; the direct moments-fold
    is one step. Both bottom out in the SAME
    :func:`~orpheus.numerics.spaces.radial_characteristic_space.fold_moments_to_radial_characteristic`
    kernel — principled-equivalent (~ULP), NOT bit-identical (the removed
    round-trip's per-ordinate ``·w`` reassociates). Carrying meshes are 1-D
    curvilinear (R12a) so ``fission_source`` is ``(ng, nx)``; the fold takes the
    unit-ℓ ``[None]`` axis (``RadialCharacteristicReconstruction.apply`` guards it).
    """
    if sn_mesh.radial_characteristic_space is None:
        return None
    from orpheus.sn.operators.radial_characteristic import (
        RadialCharacteristicReconstruction,
    )
    return RadialCharacteristicReconstruction(sn_mesh).apply(
        np.asarray(fission_source)[None],
    )


def _coupled_flux_state(
    psi_a: "TimedFullField", sn_mesh: "SNMesh",
) -> "CoupledField":
    r"""Pair a System-A FLUX iterate with a zero System-B flux composite.

    The coupled cold-start / iterate birth on a carrying mesh (B.2d — the
    pair is born native; there is no fused 3-block to split). ψ_B's zero
    flux composite comes from the presence-gated
    :meth:`~orpheus.transport.radial_characteristic_composite.RadialCharacteristicComposite.from_mesh`.
    """
    return CoupledField(
        systems=(psi_a, RadialCharacteristicComposite.from_mesh(sn_mesh)),
    )


def _coupled_source_state(
    q_a: "FullField", q_half: "RadialCharacteristicSourceSink | None",
    sn_mesh: "SNMesh", *, context: str,
) -> "CoupledField":
    r"""Pair a System-A SOURCE composite with its q½ System-B member.

    The coupled rhs birth on a carrying mesh: ``q_half`` is the unified q½
    fold (:func:`_radial_characteristic_source_from_per_ordinate` /
    :func:`_radial_characteristic_fission_seed`), lifted role-preserving
    onto System B's composite. ``None`` refuses loudly — a carrying mesh's
    rhs MUST carry the true q½ (the direct ψ½ solve consumes it).
    """
    if q_half is None:
        raise ValueError(
            f"{context}: a carrying mesh's coupled rhs requires the q½ fold "
            f"(got None) — the joint solve consumes System B's true source."
        )
    return CoupledField(
        systems=(q_a, RadialCharacteristicComposite.from_unified(q_half)),
    )


def _select_si_resolvent(
    LC: "InvertibleOperator", S: "ScatteringOperator",
    B: "LinearOperator[FullField, FullField]",
    sn_mesh: "SNMesh", inner_schedule: str,
) -> "tuple[InvertibleOperator | ScheduledInvertibleOperator, tuple[LinearOperator[FullField], ...]]":
    r"""Pick the ``(resolvent, gains)`` for the SEEDLESS within-group SI
    driver per ``inner_schedule`` — the single source of truth for the
    Jacobi/G-S choice.

    Seedless-only since B.2d: a seed-carrying mesh takes the coupled
    block-native arm in :func:`_within_group_si` and never reaches this
    selector (its M/N splitting is fixed by
    :func:`~orpheus.sn.coupled_system.build_within_group_system`).

    * ``"jacobi"`` (or any 1-D mesh) → ``(L+C, (S, B_a))``: the boundary
      lagged as an external gain (inter-sweep Jacobi — every geometry).
    * ``"gauss_seidel"`` on a multi-D Cartesian mesh →
      ``((L+C) - B_lower, (S, B_upper))``: the regular splitting
      ``(L+C−B) = M − B_upper`` (#226 §17 W2).  ``B`` splits under the
      octant-group schedule (:meth:`SNBoundaryOperator.split`); the
      strictly-lower half folds into the REIFIED forward
      :class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`
      (whose ``solve`` is the octant-group forward substitution), and the
      complement lags as an ordinary external gain — structurally congruent
      with the Jacobi arm, so the driver needs no case split.  ``S`` stays a
      lagged gain in BOTH (only the boundary coupling gets G-S; the sweep
      never re-scatters mid-sweep).

    1-D falls back to Jacobi: boundary G-S is a no-op on the scattering-
    dominated 1-D regime AND the 1-D scan is not a wavefront.  The converged
    fixed point is identical either way — this only selects the SI spectral
    rate.

    C5.4 (#225): the G-S gate is the GENUINE condition ``is_cartesian and
    not is_1d`` — the pre-C5.4 ``reduced is None`` proxy was 2-D-Cartesian
    by coincidence only. ``SweepSchedule.gauss_seidel`` and the scheduled
    sweep are d-generic (C3); d=3 G-S FP-invariance is value-gated by the
    C5.5 Mode-9 mixed-BC box (vv Mode 9 — never trust a splitting on a
    degenerate regime alone).
    """
    if inner_schedule not in ("jacobi", "gauss_seidel"):
        raise ValueError(
            f"Unknown inner_schedule: {inner_schedule!r}. "
            f"Valid choices are 'gauss_seidel' (boundary G-S, multi-D "
            f"Cartesian) or 'jacobi' (the splitting-invariant control)."
        )
    if (
        inner_schedule == "gauss_seidel"
        and sn_mesh.is_cartesian
        and not sn_mesh.is_1d
    ):
        from .loss_representation.sweep_schedule import SweepSchedule
        from .operators.boundary import SNBoundaryOperator

        # Multi-D Cartesian ⟹ SEEDLESS ⟹ B is the plain SNBoundaryOperator
        # (B_a alone; no ray block). The schedule split lives on B_a, never the
        # B_a + B_b composite (RULING P1 corollary) — this narrowing asserts that
        # invariant (a seed-carrying composite would be curvilinear, not
        # Cartesian, so it never reaches here).
        if not isinstance(B, SNBoundaryOperator):
            raise TypeError(
                "boundary Gauss-Seidel split requires the plain "
                "SNBoundaryOperator (a seedless multi-D Cartesian mesh); got "
                f"{type(B).__name__} — a seed-carrying composite must not reach "
                "the G-S schedule path (RULING P1: gradings live on B_a)."
            )
        parts = B.split(SweepSchedule.gauss_seidel(sn_mesh))
        return LC - parts.lower, (S, parts.upper)
    return LC, (S, B)


def _within_group_si(
    system: "WithinGroupSystem",
    sn_mesh: "SNMesh", *, inner_schedule: str, max_iter: int, tol: float,
) -> "tuple[SourceIteration[Any], CoupledInvertibleOperator | InvertibleOperator | ScheduledInvertibleOperator, tuple[LinearOperator, ...], bool]":
    r"""SourceIteration driver on the within-group system ``A = M − N``.

    Single source of truth (Cardinal Rule 2 / Phase 1 R1) for the
    :class:`~orpheus.numerics.iteration.SourceIteration` construction shared by
    the eigenvalue (:meth:`SNSolver._solve_source_iteration`) and fixed-source
    (:func:`_solve_fixed_source_si`) paths — the SI sibling of
    :func:`_within_group_krylov` (both inner methods have ONE construction
    helper consumed by their eigenvalue + fixed-source sites). Consumes the
    :class:`~orpheus.sn.coupled_system.WithinGroupSystem` record — the loss
    and its splitting from the ONE construction site
    (:func:`~orpheus.sn.coupled_system.build_within_group_system`).

    Two structurally-dispatched arms (B.2d DP-seedless — the coupled
    carrier appears exactly where System B exists):

    * **coupled** (the record's ``resolvent`` is the
      :class:`~orpheus.sn.coupled_system.CoupledInvertibleOperator` — a
      seed-carrying 1-D curvilinear mesh): the block-native driver
      ``ψ ← M⁻¹(q + N·ψ)`` on the ``CoupledField [ψ_A, ψ_B]`` iterate.
      Never windowed (carrying ⟹ 1-D, R12a) and never schedule-split
      (G-S is multi-D Cartesian ⟹ seedless, RULING P1) — both machineries
      are bypassed structurally, and ``inner_schedule`` is inert here (the
      1-D Jacobi fallback the seedless arm spells explicitly).
    * **seedless**: exactly the pre-B.2d composition — the schedule
      splitting (:func:`_select_si_resolvent`, Jacobi vs boundary-G-S),
      the INVERSE build (``base_resolvent.inverse()`` — the
      :class:`~orpheus.sn.operators.sweep_operator.SweepOperator`), and
      the Phase-5a angular-windowing composition (:func:`_maybe_window` —
      2-D Cartesian holds the iterate as harmonic moments via
      ``P @ A.inverse()``).

    Returns ``(si, base_resolvent, gains, windowed)``:

    * ``si`` — the :class:`SourceIteration` primitive;
    * ``base_resolvent`` — the un-inverted FORWARD ``M`` (the fixed-source
      path needs it for the one-shot full-angular reconstruction of
      ``Solution.angular_flux``);
    * ``gains`` — the lagged couplings ``N`` actually driven (the record's,
      except the G-S arm's re-split ``(S, B_upper)``), also needed to
      rebuild the converged source for that reconstruction;
    * ``windowed`` — whether the iterate is the moment representation (2-D
      Cartesian) vs full-angular (curvilinear / 1-D).

    Both paths forward their caller's ``inner_schedule`` (default boundary
    Gauss-Seidel on 2-D Cartesian — `#218
    <https://github.com/deOliveira-R/ORPHEUS/issues/218>`_ closed the
    eigenvalue-inner gap; the eigenvalue path reads ``SNSolver.inner_schedule``,
    the fixed-source path its ``inner_schedule`` argument).
    """
    from orpheus.numerics.iteration import SourceIteration

    if isinstance(system.resolvent, CoupledInvertibleOperator):
        # The ψ½ coupled block-native arm (B.2d): the record's splitting IS
        # the driver's — M⁻¹ = the joint sweep, N = the coupled gain grid.
        si = SourceIteration(
            system.resolvent.inverse(), *system.gains,
            max_iter=max_iter, tol=tol,
        )
        return si, system.resolvent, system.gains, False
    # Seedless: the record's gains are the (S, B_a) pair — loud on drift.
    S, B = system.gains
    if not isinstance(S, ScatteringOperator):
        raise TypeError(
            f"_within_group_si: the seedless record's first gain must be "
            f"the ScatteringOperator (the builder's (S, B_a) convention); "
            f"got {type(S).__name__}."
        )
    base_resolvent, gains = _select_si_resolvent(
        system.resolvent, S, B, sn_mesh, inner_schedule,
    )
    step, windowed = _maybe_window(base_resolvent.inverse(), S, sn_mesh)
    si = SourceIteration(step, *gains, max_iter=max_iter, tol=tol)
    return si, base_resolvent, gains, windowed


# ═══════════════════════════════════════════════════════════════════════
# Solver class (EigenvalueSolver protocol)
# ═══════════════════════════════════════════════════════════════════════

class SNSolver:
    """Unified SN eigenvalue solver satisfying the EigenvalueSolver protocol.

    Constructs the operator triple :math:`(L, S, F)` at construction
    time and routes ``solve_fixed_source`` through one of two inner-
    solver paths:

    * ``"source_iteration"`` — sweep-driven within-group fixed-point
      iteration (WDD asymmetric closure; ERR-026-affected for
      curvilinear).  Bit-identical to the Wave A-D path.
    * ``"krylov"`` — GMRES on
      :class:`~orpheus.sn.operators.streaming.InvertibleOperator` (the algebraic
      composition ``L + C``; symmetric closure) with the sweep as
      preconditioner.  Closes ERR-026 on curvilinear; bit-identical
      math to the legacy BiCGSTAB FD path on Cartesian.

    The legacy ``"bicgstab"`` value is no longer accepted — call sites
    must migrate to ``"krylov"``.

    Parameters
    ----------
    sn_mesh : SNMesh — augmented geometry (wraps Mesh1D or Mesh2D with
        precomputed streaming stencil + materials dict + ``ng``).
        Issue #197 PR-TYPED-0: ``sn_mesh.materials`` IS the single
        source of truth for cross sections and group count; the
        legacy ``materials`` / ``n_groups`` SNSolver constructor
        parameters were retired (aggressive retirement per
        ``feedback_aggressive_retirement``).
    inner_solver : "source_iteration" or "krylov".
    scattering_order : int — Legendre order for scattering (0 = P0).
    keff_tol, flux_tol : outer iteration convergence.
    max_inner, inner_tol : inner iteration parameters.
    """

    def __init__(
        self,
        sn_mesh: SNMesh,
        inner_solver: str = "source_iteration",
        scattering_order: int = 0,
        keff_tol: float = 1e-7,
        flux_tol: float = 1e-6,
        max_inner: int = 200,
        inner_tol: float = 1e-8,
        inner_schedule: str = "jacobi",
    ):
        if inner_solver not in ("source_iteration", "krylov"):
            raise ValueError(
                f"Unknown inner solver: {inner_solver!r}. "
                f"Valid choices are 'source_iteration' or 'krylov'. "
                f"(The legacy 'bicgstab' alias was retired in Wave E "
                f"Round 2; use 'krylov' which routes through "
                f"InvertibleOperator (L + C) with the sweep as "
                f"preconditioner.)"
            )
        if inner_schedule not in ("jacobi", "gauss_seidel"):
            raise ValueError(
                f"Unknown inner_schedule: {inner_schedule!r}. "
                f"Valid choices are 'gauss_seidel' (boundary G-S, multi-D "
                f"Cartesian — auto-falls-back to Jacobi on 1-D) or 'jacobi'."
            )
        self.sn_mesh = sn_mesh
        self.quad = sn_mesh.quad
        self.inner_solver = inner_solver
        # SI BOUNDARY splitting for the eigenvalue inner (#218 — the eigenvalue
        # SI now CAN reach the boundary-Gauss-Seidel accelerator the
        # fixed-source path got in Phase 3, via the shared ``_within_group_si``
        # builder; validated SI(G-S)≡Krylov≡k_inf).  DEFAULT stays ``"jacobi"``:
        # the eigenvalue inner is warm-started (the G-S rate benefit is modest
        # there), and a schedule change shifts the converged k_eff by ~inner_tol
        # (1e-10-scale — same fixed point, vv Mode 9; only the inner SI stopping
        # differs), which the keff_tol-tight regression snapshots cannot absorb.
        # ``"gauss_seidel"`` is opt-in (2-D Cartesian; ``_select_si_resolvent``
        # auto-falls-back to Jacobi on 1-D / curvilinear).
        self.inner_schedule = inner_schedule
        self.scattering_order = scattering_order
        self.keff_tol = keff_tol
        self.flux_tol = flux_tol
        self.max_inner = max_inner
        self.inner_tol = inner_tol

        # Issue #197 PR-TYPED-0: materials + ng are read from sn_mesh,
        # not from constructor parameters.  ``sn_mesh.ng`` is the
        # validated single source of truth (raises
        # ``InconsistentMaterialsError`` if materials disagree).
        materials = sn_mesh.materials
        self.ng = sn_mesh.ng

        # Issue #197 PR-TYPED-1: the canonical XS state collapses to ONE
        # attribute — ``self.mat_xs`` is the :class:`MaterialXSField`
        # wrapping both the per-material :class:`Mixture` data AND the
        # per-cell typed views.  Every operator (L, C, S, F) reads cross
        # sections through this single source of truth; the seven
        # separate per-XS attributes (sig_t, sig_a, sig_p, chi, sig_s,
        # sig2, sig_s0) plus ``_cells_by_mat`` / ``_sig2_sum`` collapse
        # into ``self.mat_xs.*`` accessors.
        #
        # The thin ``sig_t / sig_a / sig_p / chi`` read-through
        # properties below preserve the SNSolver API surface for the
        # one-PR migration window (downstream consumers that read
        # ``solver.sig_t`` etc.).  PR-TYPED-2 retires them by rewiring
        # consumers to read ``solver.mat_xs.total_cross_section`` etc.
        # directly.
        self.mat_xs = sn_mesh.material_xs_field()

        # __debug__ cell-flattening invariant pinning (formerly at
        # construction of self.sig_t — now exercised through the
        # mat_xs.total_cross_section accessor, populated lazily).
        if __debug__:
            xs_check = assemble_cell_xs(materials, sn_mesh.mat_map)
            _sig_t_old = xs_check.sig_t.reshape(*sn_mesh.spatial_shape, self.ng)
            assert np.array_equal(
                _sig_t_old,
                np.moveaxis(self.mat_xs.total_cross_section, 0, -1),
            ), "PR-INDEX-3 cell-flattening invariant broke"

        # Scattering order — clamp to the minimum Legendre count
        # available across all materials.
        L = min(
            scattering_order,
            min(len(m.SigS) - 1 for m in materials.values()),
        )
        self.scattering_order = L

        # Weight normalization (1/sum(w) — works for both GL and Lebedev)
        self.weight_norm = 1.0 / sn_mesh.quad.weights.sum()

        # Persistent boundary flux state (passed to sweep).
        # Issue #197 PR-TYPED-2: typed :class:`AngularBoundaryFlux` replaces
        # the stringly-typed ``psi_bc: dict``.  Per-face buffers
        # become named attributes; typos surface as AttributeError.
        self._boundary_flux = AngularBoundaryFlux.zeros_on(sn_mesh)

        # Phase 3 measurement seam: total inner (within-group) SI/Krylov
        # iterations consumed across the eigenvalue outer loop — the
        # measurand for the SI spectral-rate / Gauss-Seidel-recovery
        # diagnostics.  Accumulated per outer step in
        # ``_solve_source_iteration`` / ``_solve_krylov`` and read by
        # ``solve_sn`` into ``IterationHistory.total_inner_iterations``.
        # Fresh per solve (a new SNSolver is built per ``solve_sn`` call).
        self._total_inner_iterations = 0

        # Volume array for keff computation
        self.volume = sn_mesh.volumes

        # ── Operator triple — Wave E Round 2 -----------------------------
        # The (L, S, F) algebra-of-record framing.  Constructed once at
        # __init__ so downstream consumers (the iteration primitives, the
        # _solve_krylov path, future sensitivity/adjoint hooks) see a
        # consistent operator triple over the lifetime of this solver.
        # Issue #197 PR-TYPED-1: both S and F now consume the single
        # ``self.mat_xs`` — the per-material dispatch lives inside
        # :class:`MaterialXSField`'s typed verbs, not on the operators.
        self.scattering_op = ScatteringOperator.from_solver_data(
            mat_xs=self.mat_xs,
            quadrature=sn_mesh.quad,
            scattering_order=self.scattering_order,
            full_field_space=sn_mesh.full_field_space,
        )
        self.fission_op = FissionOperator.from_solver_data(
            mat_xs=self.mat_xs,
            full_field_space=sn_mesh.full_field_space,
        )
        # The full transport operator :math:`L = \Omega\cdot\nabla + \Sigma_t`
        # as the algebraic sum streaming + the collision multiplier
        # ``M[σ_t]`` = :class:`InvertibleOperator`.
        # :meth:`InvertibleOperator.apply` returns ``(L_stream + C)·ψ`` in
        # :class:`TimedFullField` form (the typed direct-sum carrier);
        # :meth:`InvertibleOperator.solve` consumes ``initial_guess`` for
        # the Carlson seed (R-1 Phase 1.2 unification).
        self.L = (
            StreamingOperator(sn_mesh)
            + MultiplicationOperator(
                coefficient=self.mat_xs.total_cross_section_field,
                space=sn_mesh.full_field_space,
            )
        )
        self.S = self.scattering_op
        self.F = self.fission_op

        # ── Sweep cache (Issue #196 Phase G Step 2.5c) ───────────────
        # Two-stratum cache: GeometryCoefficients built once at __init__
        # (geometry × quadrature only); CollisionCache built once after
        # σ_t binding.  Hot path consumes (geom, coll) without per-cell
        # StreamingTerms allocation.  Only applicable to 1-D meshes with
        # ReducedStreamingOperator — 2-D Cartesian uses the wavefront path.
        self.geom_cache: GeometryCoefficients | None = None
        self.coll_cache: CollisionCache | None = None
        # The two-stratum scan cache feeds the DAG-FREE scan strategies
        # (CumprodScan / ScanMarch) ONLY — its σ_t stratum is the closed-form
        # affine recurrence ``affine_scan_coefficients`` (the scan-family
        # triple), which a NON-affine-scannable scheme (LinearDiscontinuous,
        # #158) does not supply.  Such schemes run on the DAG wavefront
        # (FullFieldWavefront), which consumes the per-cell ``cell_kernel_batch``
        # directly — never the scan cache.  Build the cache only when the scan
        # path can actually be selected (DD keeps its bit-identical cache).
        if sn_mesh.reduced is not None and sn_mesh.scheme.is_affine_scannable:
            self.geom_cache = GeometryCoefficients.from_mesh_and_quad(sn_mesh)
            # No bridge needed: ``mat_xs.total_cross_section`` is the
            # principled ``(ng, nx)`` 1-D layout the cache expects
            # (rank-d (N, ng, *spatial); no phantom ny axis to drop).
            sig_t_1d = self.mat_xs.total_cross_section  # (ng, nx)
            self.coll_cache = CollisionCache.from_geometry(
                self.geom_cache, sig_t_1d, sn_mesh.scheme,
            )
            # Stash the caches on the SNMesh so the sweep can find them
            # without threading a solver reference through ``transport_sweep``.
            sn_mesh._geom_cache = self.geom_cache  # type: ignore[attr-defined]
            sn_mesh._coll_cache = self.coll_cache  # type: ignore[attr-defined]

    # Issue #197 PR-TYPED-2: the PR-TYPED-1 transient read-through
    # shims (``sig_t``, ``sig_a``, ``sig_p``, ``chi``, ``sig_s``,
    # ``sig2``, ``sig_s0``, ``_cells_by_mat``) have been RETIRED.
    # Every downstream consumer now reads ``self.mat_xs.*`` directly:
    # ``mat_xs.total_cross_section``,
    # ``mat_xs.absorption_cross_section``,
    # ``mat_xs.fission_production``,
    # ``mat_xs.emission_spectrum``,
    # ``mat_xs.sig_s_legendre(mid)``,
    # ``mat_xs.n2n_matrix(mid)``,
    # ``mat_xs.cells_by_material``.

    def rebind_cross_sections(self, new_sig_t: np.ndarray) -> None:
        """Rebind the total cross-section and rebuild only :class:`CollisionCache`.

        :class:`GeometryCoefficients` survives — Stratum 1 is geometry-only.
        Only the σ_t-dependent Stratum 2 rebuilds.  Used by depletion /
        thermal-feedback consumers.

        Parameters
        ----------
        new_sig_t
            New total cross-section in the principled ``(ng, nx, ny)``
            layout (Issue #196 PR-INDEX-3).

        Notes
        -----
        Issue #197 PR-TYPED-1 — ``rebind_cross_sections`` overrides
        ``self.mat_xs._sig_t_cell`` directly (without re-deriving from
        materials) because the depletion / thermal-feedback consumer
        adjusts σ_t per-cell without revisiting the per-material data.
        """
        # Override the lazy cache on mat_xs.  Force the dense
        # per-cell view to be populated first so the other cell views
        # (sig_a, sig_p, chi) match the rebind contract.
        _ = self.mat_xs.absorption_cross_section
        self.mat_xs._sig_t_cell = new_sig_t
        # Rebuild the composite so the rebound σ_t flows into the collision
        # diagonal C (the streaming leaf L is σ-free since #257 S8b — only C
        # carries σ; the composite single-sources it from C.sigma).
        self.L = (
            StreamingOperator(self.sn_mesh)
            + MultiplicationOperator(
                coefficient=self.mat_xs.total_cross_section_field,
                space=self.sn_mesh.full_field_space,
            )
        )
        if self.geom_cache is not None:
            sig_t_1d = self.mat_xs.total_cross_section
            self.coll_cache = CollisionCache.from_geometry(
                self.geom_cache, sig_t_1d, self.sn_mesh.scheme,
            )
            self.sn_mesh._coll_cache = self.coll_cache  # type: ignore[attr-defined]

    def initial_flux_distribution(self) -> np.ndarray:
        """Initial scalar flux guess: ones(ng, nx, ny).

        Issue #196 PR-INDEX-5: principled layout.
        """
        return np.ones((self.ng, *self.sn_mesh.spatial_shape))

    def compute_fission_source(
        self, flux_distribution: np.ndarray, keff: float,
    ) -> np.ndarray:
        """Fission source: χ · (νΣ_f · φ) / k.

        Thin delegator to :meth:`FissionOperator.apply` (Wave D Issue 13).
        The :math:`1/k` division stays at this level — the fission
        operator is a *linear* operator (Wave A Issue 1 Protocol);
        :meth:`FissionOperator.apply` returns :math:`F\\,\\phi` and the
        eigenvalue scaling lives here.

        Issue #196 PR-INDEX-5: ``flux_distribution`` is principled
        ``(ng, nx, ny)``.  No bridges — the PR-INDEX-4 transpose pair
        is GONE.
        """
        return self.fission_op.apply(flux_distribution) / keff

    def solve_fixed_source(
        self, fission_source: np.ndarray, flux_distribution: np.ndarray,
    ) -> np.ndarray:
        """Solve the within-group transport equation for given fission source.

        Returns updated scalar flux ``(ng, nx, ny)``.
        """
        if self.inner_solver == "source_iteration":
            return self._solve_source_iteration(fission_source, flux_distribution)
        if self.inner_solver == "krylov":
            return self._solve_krylov(fission_source, flux_distribution)
        # Should be unreachable — __init__ validated the choice.
        raise ValueError(f"Unknown inner solver: {self.inner_solver}")

    def compute_group_production_rate(
        self, flux_distribution: np.ndarray,
    ) -> np.ndarray:
        r"""Per-group volume-integrated neutron production rate, shape ``(ng,)``.

        Component :math:`r_g` is

        .. math::

            r_g \;=\; \int_V \nu \Sigma_{f,g}(\mathbf{r})\,\phi_g(\mathbf{r})\,dV
                       \;+\; 2 \int_V \sum_{g' } \Sigma_{2,g'\to g}(\mathbf{r})
                                                 \,\phi_{g'}(\mathbf{r})\,dV

        i.e. the per-group fission-neutron production plus the per-group
        ``(n, 2n)`` contribution (the factor of 2 accounts for the
        two-neutron-out yield).  Fission is integrated against
        ``mesh.volume_measure`` (Issue 9.6 wiring); ``(n, 2n)`` runs the
        existing per-material loop because the ``sig2`` matrices are
        keyed on material rather than cell.

        The output is the natural diagnostic intermediate for spectral
        analysis (per-group production rates are reactor-physics-meaningful
        quantities).  ``compute_keff`` consumes it via ``.sum()``.
        """
        ng = self.ng
        mu = self.sn_mesh.volume_measure

        # Fission production: ∫ νΣ_f · φ dV, vectorised over groups.
        # Issue #196 PR-INDEX-5: both ``mat_xs.fission_production`` and
        # ``flux_distribution`` are principled ``(ng, nx, ny)``.  The
        # named intermediate ``per_cell_per_group`` has units ``[1/s]``
        # per cell per group (a reactor-physics quantity — coding-
        # elegance Pattern 3).  ``volume_measure`` consumes a flat
        # ``(N_cells, ng)`` view (Issue 9.6 wiring); reshape via
        # ``transpose(1, 2, 0)`` then flatten the spatial axes.
        # Issue #197 PR-TYPED-2: consumes mat_xs directly (no shim).
        per_cell_per_group = np.einsum(
            "g...,g...->g...", self.mat_xs.fission_production, flux_distribution,
        )
        rate = mu(np.moveaxis(per_cell_per_group, 0, -1).reshape(-1, ng))

        # (n,2n) contribution — Issue #197 PR-TYPED-1: the per-material
        # dispatch loop lives ONLY inside
        # :meth:`MaterialXSField.add_n2n_to_group_rate`.
        self.mat_xs.add_n2n_to_group_rate(
            rate, flux_distribution, self.volume,
        )

        return rate

    def compute_group_absorption_rate(
        self, flux_distribution: np.ndarray,
    ) -> np.ndarray:
        r"""Per-group volume-integrated absorption rate, shape ``(ng,)``.

        Component :math:`a_g = \int_V \Sigma_{a,g}(\mathbf{r})\,\phi_g(\mathbf{r})\,dV`.

        Volume-integrated via ``mesh.volume_measure`` (Issue 9.6 wiring).
        ``.sum()`` of this vector is the ABSORPTION term of the
        ``compute_keff`` denominator (net removal = absorption + leakage
        − (n,2n) emission since #291/R7).

        Issue #196 PR-INDEX-5: ``flux_distribution`` is principled
        ``(ng, nx, ny)``.
        """
        ng = self.ng
        mu = self.sn_mesh.volume_measure
        # Issue #197 PR-TYPED-2: consumes mat_xs directly (no shim).
        per_cell_per_group = np.einsum(
            "g...,g...->g...", self.mat_xs.absorption_cross_section, flux_distribution,
        )
        return mu(np.moveaxis(per_cell_per_group, 0, -1).reshape(-1, ng))

    def compute_production_rate(self, flux_distribution: np.ndarray) -> float:
        r"""Total volume-integrated neutron production rate (scalar).

        :math:`P(\phi) = \int_V \sum_g \nu \Sigma_{f,g} \phi_g\,dV
        + 2 \int_V \sum_g \sum_{g'} \Sigma_{2,g'\to g} \phi_{g'} \,dV`.

        The fission term is the typed volume-integrated reaction rate
        :class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`
        over :math:`\nu\Sigma_f` — ``∫_V ⟨νΣf, φ⟩ dV`` — the single source of the
        ``Σx·φ`` contraction and its volume integral. The :math:`(n,2n)` channel
        is an **explicit additive term** (a second neutron-multiplying reaction,
        NOT a ``⟨Σx,φ⟩`` rate); it reuses the single ``add_n2n_to_group_rate``
        machinery and is exactly zero on a no-(n,2n) mixture.

        This is the canonical scale anchor for the SN eigenmode:
        :func:`orpheus.numerics.eigenvalue.power_iteration` renormalises
        :math:`\phi` to unit production rate at each outer step (ERR-052).

        Role split (R7, #259): this TOTAL physical production — fission
        plus the (n,2n) emission — is the renormalisation scale anchor
        ONLY. The k numerator in :meth:`compute_keff` is fission-only,
        because the posed eigenproblem scales only fission by
        :math:`1/k`; the (n,2n) gain sits on the net-removal side there.
        """
        fission = IntegratedReactionRate(
            self.mat_xs.fission_production_field
        ).evaluate(flux_distribution)
        n2n_rate = np.zeros(self.ng)
        self.mat_xs.add_n2n_to_group_rate(n2n_rate, flux_distribution, self.volume)
        return float(fission + n2n_rate.sum())

    def compute_keff(self, flux_distribution: np.ndarray) -> float:
        r"""k of the POSED eigenproblem: fission production over net removal.

        .. math::

            k \;=\; \frac{R_{\nu\Sigma_f}(\phi)}
                    {R_{\Sigma_a}(\phi) \;+\; L \;-\; E_{2n}(\phi)}

        Every inner solve poses the eigenproblem with ONLY fission scaled
        by :math:`1/k` — scattering and the (n,2n) emission are plain
        gains inside :meth:`solve_fixed_source` — so the reported k must
        be the eigenvalue of exactly that problem (#291 leakage omission
        + the R7 (n,2n) convention, 2026-07-03; an estimator that
        disagrees with its own posed problem converges cleanly to a
        non-eigenvalue ratio):

        * **numerator** — the typed volume-integrated FISSION production
          :math:`R_{\nu\Sigma_f}(\phi)`
          (:class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`,
          the :math:`\phi^\dagger\!=\!1` degenerate of the homogenization
          PG bilinear). The (n,2n) emission is NOT production here —
          contrast :meth:`compute_production_rate`, the ERR-052 scale
          anchor, which keeps total physical production.
        * **denominator** — net removal, assembled so no term CAN be
          forgotten: absorption :math:`R_{\Sigma_a}` (``absorption_xs``
          counts the (n,2n) COLLISION once), **plus** the net
          vacuum-boundary leakage :math:`L` (#291 — the historically
          omitted term; see :meth:`_boundary_leakage_rate`), **minus**
          the (n,2n) EMISSION :math:`E_{2n}(\phi) = \int_V \sum_{g,g'}
          2\,\Sigma_{2,g'\to g}\,\phi_{g'}\,dV` (a gain reduces net
          removal).

        Balance identity at any converged eigenpair:
        :math:`R_{\nu\Sigma_f}/k = R_{\Sigma_a} + L - E_{2n}`.

        On an all-reflective (lattice) problem :math:`L` is a STRUCTURAL
        zero, and on a Σ₂-free mixture :math:`E_{2n}` is exactly ``0.0``
        — so this reduces **bit-identically** to the historical lattice
        functional ``production / absorption``.

        The per-group :meth:`compute_group_production_rate` /
        :meth:`compute_group_absorption_rate` remain available as
        spectral diagnostics (not on the keff path).
        """
        production = IntegratedReactionRate(
            self.mat_xs.fission_production_field
        ).evaluate(flux_distribution)
        absorption = IntegratedReactionRate(
            self.mat_xs.absorption_cross_section_field
        ).evaluate(flux_distribution)
        emission_n2n = np.zeros(self.ng)
        self.mat_xs.add_n2n_to_group_rate(
            emission_n2n, flux_distribution, self.volume,
        )
        leakage = self._boundary_leakage_rate(production)
        return production / (absorption + leakage - emission_n2n.sum())

    def _boundary_leakage_rate(self, fission_production: float) -> float:
        r"""Net neutron outflow rate through the vacuum boundary faces [1/s].

        .. math::

            L \;=\; \sum_{f\,\in\,\text{vacuum}} \oint_{f} dA\,
                    \sum_g J_g(\mathbf{r}_f)
            \,, \qquad
            J_g \;=\; \sum_m (\Omega_m\cdot\hat n_f)\, w_m\, \psi_{m,g}

        — the face-area integral of the boundary trace's net outward
        current (:meth:`AngularBoundaryFlux.net_current
        <orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux.net_current>`,
        the single source of the :math:`\Omega\cdot\hat n\,w`
        contraction), read from the trace of the last inner solve
        (``self._psi_typed.boundary``). On the converged trace a vacuum
        face's inflow slots are zero, so net = outflow; the signed form
        stays honest if a prescribed-inflow law ever lands.

        Reflective faces are a **structural zero**: the reflective law
        equates inflow to the reflected outflow exactly, so their net
        current vanishes by construction — they are SKIPPED, never
        accumulated, which keeps all-reflective problems' float
        arithmetic bit-identical to the lattice functional (no
        ±cancelling angular-sum noise enters the denominator).

        Scale bridge: the stored trace belongs to the UN-renormalised
        last inner iterate, while the estimator's :math:`\phi` may be
        its renormalised multiple
        (:func:`~orpheus.numerics.eigenvalue.power_iteration` divides by
        the production rate between the solve and the k-update). Leakage
        is degree-1 homogeneous in :math:`\psi`, so it is rescaled by
        the fission-production ratio of the two — exactly ``1.0`` when
        the caller passes the returned flux itself. Contract: the flux
        handed to :meth:`compute_keff` must be a scalar multiple of the
        last inner solve's flux (true for ``power_iteration`` and for
        every manual solve-then-estimate loop).

        Raises
        ------
        RuntimeError
            If a vacuum face exists but no inner solve has stored a
            boundary trace yet — the leakage term cannot be answered
            honestly, and answering without it would silently reproduce
            the #291 omission (fail loud; never return a non-eigenvalue).
        """
        vacuum_faces = [
            name for name, op in self.sn_mesh.bc.items()
            if op.kind == "vacuum"
        ]
        if not vacuum_faces:
            return 0.0
        psi = getattr(self, "_psi_typed", None)
        phi_of_trace = getattr(self, "_phi_of_trace", None)
        if psi is None or phi_of_trace is None:
            raise RuntimeError(
                "compute_keff on a vacuum-bounded problem needs the "
                "boundary trace of an inner solve (call "
                "solve_fixed_source first): the leakage term of the k "
                "denominator is read from psi.boundary, and answering "
                "without it would silently drop leakage (#291)."
            )
        rate = 0.0
        trace = _system_a_member(psi).boundary  # System A carries the trace
        # #289-F2 role parse: the widened FullField.boundary slot erases the
        # family; the SN leakage read needs the ANGULAR trace's net_current.
        if not isinstance(trace, AngularBoundaryFlux):
            raise TypeError(
                f"compute_keff: the converged iterate's trace must be an "
                f"AngularBoundaryFlux; got {type(trace).__name__}."
            )
        for face in vacuum_faces:
            net_current = trace.net_current(face)  # (ng, *face_spatial)
            rate += float(np.sum(net_current * self._face_area_of(face)))
        reference = IntegratedReactionRate(
            self.mat_xs.fission_production_field
        ).evaluate(phi_of_trace)
        if reference <= 0.0:
            raise RuntimeError(
                "leakage scale bridge is degenerate: the last inner "
                "solve's flux carries non-positive fission production."
            )
        return rate * (fission_production / reference)

    def _face_area_of(self, face: str) -> "float | np.ndarray":
        r"""Spatial measure of a boundary face, matching ``volume_measure``.

        1-D: the scalar face area from
        :attr:`~orpheus.transport.mesh.material_mesh.MaterialMesh.areas`
        — ``1.0`` (slab, per unit cross-section), :math:`2\pi R`
        (cylinder, per unit height), :math:`4\pi R^2` (sphere) — the
        same geometric convention the cell volumes integrate under, so
        the balance identity closes. 2-D Cartesian: the per-edge-cell
        transverse widths (unit depth), broadcast against the
        ``(ng, n_edge)`` net current. 3-D vacuum leakage has no consumer
        yet and fails loud rather than guessing the transverse-area
        product's cell ordering.
        """
        mesh = self.sn_mesh
        if mesh.ndim == 1:
            areas = mesh.areas
            return float(areas[0] if face == "xmin" else areas[-1])
        if mesh.ndim == 2:
            axis_index = {"x": 0, "y": 1}[face[0]]
            transverse = mesh.axes[1 - axis_index]
            return np.diff(np.asarray(transverse.edges, dtype=float))
        raise NotImplementedError(
            f"boundary leakage for a {mesh.ndim}-D vacuum face "
            f"({face!r}): wire the transverse-area product when the "
            f"first 3-D vacuum eigenvalue consumer arrives."
        )

    def converged(
        self, keff: float, keff_old: float,
        flux_distribution: np.ndarray, flux_old: np.ndarray,
        iteration: int,
    ) -> bool:
        if iteration <= 2:
            return False
        dk = abs(keff - keff_old)
        dphi = np.linalg.norm(flux_distribution - flux_old) / \
            max(np.linalg.norm(flux_distribution), 1e-30)
        return bool(dk < self.keff_tol and dphi < self.flux_tol)

    # ── Inner solver: source iteration ────────────────────────────────

    def _solve_source_iteration(
        self, fission_source: np.ndarray, flux_distribution: np.ndarray,
    ) -> np.ndarray:
        r"""Inner within-group solve via :class:`SourceIteration` on typed AngularFlux.

        R-1 Step E (2026-05-19) — carved onto
        :class:`~orpheus.numerics.iteration.SourceIteration` consuming
        the same typed-flux operator triple as the Krylov carve:

        .. math::

            A \;=\; L + C\,, \quad
            S \;=\; \tfrac{1}{W}\,\text{full multi-group scatter}\,, \quad
            F \;=\; 0_{\rm wg}

        on :class:`~orpheus.sn.angular_flux.AngularFlux`.  The
        iteration step is

        .. math::

            \psi_{n+1} \;=\; (L + C)^{-1}\,\bigl(S\,\psi_n
                                                + F\,\psi_n
                                                + q_{\rm ext}\bigr)

        where the driver applies ``(L + C).inverse()`` — a
        :class:`~orpheus.sn.operators.sweep_operator.SweepOperator` whose
        ``apply`` IS the WDD sweep (#226 taxonomy step 3; the solver
        builds the inverse, :class:`SourceIteration` applies it).  The
        previous iterate :math:`\psi_n` travels into the sweep via the
        explicit ``initial_guess`` kwarg on the inverse's ``apply``; the
        M-M closure's ``psi_half_seed`` strategy reads it to derive the
        curvilinear Carlson coupled-pole seed (pinned by the
        seed-threading spy).

        Scope
        =====

        ALL geometries — slab, sphere, cylinder, AND 2-D Cartesian
        (Wave O "2-D SI Phase A", 2026-06-04).  The 2-D Cartesian
        eigenvalue SI inner is geometry-agnostic: it is the structural
        twin of :meth:`_solve_krylov` (the live 2-D eigenvalue path) —
        identical composite RHS, identical loss decomposition
        (``LC = StreamingOperator + MultiplicationOperator`` — the collision
        multiplier ``C = M[σ_t]`` — plus the
        scattering ``S`` and boundary ``B`` coupling gains —
        :func:`~orpheus.sn.coupled_system.build_within_group_system`, zero
        within-group fission), identical
        ``psi_typed.interior.integrate_angular()`` reduction — differing
        ONLY in the driver (:class:`SourceIteration` vs
        :class:`KrylovAcceleration`), and neither driver carries any
        geometry dependence.  The reflective coupling rides the BARE
        sweep via the ``B`` coupling gain on the 4-face
        :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
        (:class:`SNBoundaryOperator` is natively 4-face —
        xmin/xmax/ymin/ymax — and is the SAME operator the working 2-D
        Krylov path uses).

        The legacy "B1'' face block" that the 2-D path was once said to
        lack is RETIRED — it never existed as code on this branch; it
        was a 1-D boundary-closure name fully superseded by the L2
        ``AngularBoundaryFlux`` + ``SNBoundaryOperator`` bare-boundary
        architecture (Wave O O.4a.2 / O.4b).  The historical
        ``NotImplementedError`` guard predated that architecture by two
        weeks (the body's ``S + B`` fold landed 2026-06-03) and was
        never lifted; ``solve_sn`` defaults ``inner_solver="source_-
        iteration"`` for every geometry, so the stale guard broke the
        DEFAULT 2-D Cartesian eigenvalue entry point.

        Verified: 2-D SI ≡ Krylov ≡ closed-form ``k_inf`` (1g→1.5,
        2g→1.875, 4g→1.4878) to machine precision; heterogeneous
        non-flat 2-D flux SHAPE agrees SI-vs-Krylov to ~1e-9
        (``tests/sn/eigenvalue/test_keff_2d.py::TestSIKrylov2DEquivalence``).
        """
        from orpheus.transport.fields.angular_flux import (
            AngularFlux,
        )
        from orpheus.transport.fields.angular_boundary_flux import (
            AngularBoundaryFlux,
        )
        from orpheus.transport.source_sinks import (
            AngularSourceSink,
            AngularBoundarySourceSink,
        )
        from orpheus.transport.timed_full_field import TimedFullField

        # ── Build the composite RHS ─────────────────────────────────
        # R-1 Step 4 A1 — q_ext as per-ordinate density via the canonical
        # ``AngularSourceSink.from_isotropic`` factory.  The /W projection
        # lives at the factory boundary (Pattern 7 producer-side
        # normalisation); the legacy
        # ``q_per_ord = (fission_source / sum_w)[None, :, :, :]``
        # broadcast pattern is GONE.
        q_ext_per_ord = AngularSourceSink.from_isotropic(
            fission_source, self.sn_mesh,
        )
        # D-H.1c stage 2 — q_ext composite carries:
        #   * bulk = per-ordinate source values on the L2 AngularFlux.
        #   * boundary = ZERO (the EXTERNAL boundary source — zero for
        #     vacuum/reflective; non-zero only for prescribed inflow).
        # Wave O (#208) O.4a.2 — the partner-flux seeding
        # (``boundary = self._boundary_flux``) is RETIRED: the reflective
        # inflow is no longer pre-staged into the source.  It is driven by
        # the sibling ``−B`` delivered as a SEPARATE coupling gain (O.2a;
        # since B.2d the system record's ``gains``): each SI iterate
        # adds ``B·ψ.outflow`` to ``rhs.boundary``, which ``(L+C).solve``'s
        # bare sweep reads as the inflow seed (operator.py
        # ``_solve_timed_full_field`` seeds from ``rhs.boundary``).  The
        # boundary inflow is thus a live solved unknown carried in
        # ``ψ.boundary``, not an externally-recomputed partner trace.
        q_ext_composite = TimedFullField(
            # B.5.2: q_ext IS a source — carry the AngularSourceSink bulk AND
            # the AngularBoundarySourceSink inflow trace (zero for vacuum/reflective;
            # prescribed inflow otherwise). The SI rhs q_ext + S.apply + B.apply
            # closes on AngularBoundarySourceSink (operator outputs are sources).
            interior=q_ext_per_ord,
            boundary=AngularBoundarySourceSink.zeros_on(self.sn_mesh),
            _history=(),
            history_depth=2,
        )

        # ── Build the within-group system (single source of truth —
        # :func:`~orpheus.sn.coupled_system.build_within_group_system`;
        # shared with the Krylov and fixed-source paths; the solver's
        # cached scattering operator injects through the cache seam).
        # ``#218``: the eigenvalue inner honours ``self.inner_schedule``
        # (default boundary-G-S on 2-D Cartesian; the coupled arm and 1-D
        # fall to Jacobi structurally).  Phase-5a angular-windowing folds
        # in via :func:`_maybe_window` inside the SI builder. ──────────
        system = build_within_group_system(
            self.sn_mesh, self.mat_xs, scattering_op=self.scattering_op,
        )
        si, _base, _gains, windowed = _within_group_si(
            system, self.sn_mesh,
            inner_schedule=self.inner_schedule,
            max_iter=self.max_inner, tol=self.inner_tol,
        )
        # B.2d DP-seedless: the coupled pair appears exactly where System B
        # exists; the seedless paths (windowed 2-D, G-S) stay fused.
        coupled = isinstance(system.resolvent, CoupledInvertibleOperator)

        # ── Warm start (composite / coupled pair) ───────────────────
        # SourceIteration threads the previous iterate to the inverse
        # operator's ``apply`` via the explicit ``initial_guess`` kwarg
        # (accepted-and-dropped by the direct sweeps since #282/2.5d; the
        # previous iterate's boundary trace seeds the reflective-BC
        # partner-flux state — pinned by the seed-threading spy).
        # Post-B.2d ``self._psi_typed`` is a CoupledField on a carrying
        # mesh (the previous inner's converged pair — already split) and
        # a :class:`TimedFullField` elsewhere; both propagate through the
        # iteration primitive via the ravellable protocol.
        initial_guess = getattr(self, "_psi_typed", None)
        if initial_guess is None:
            # B.5.2: cold-start iterate is an all-zeros FLUX composite,
            # decoupled from q_ext's AngularSourceSink type.  Phase 5a: when
            # windowed the bulk is a zero HarmonicMomentFlux (single-sourced
            # in :func:`_windowed_cold_start`); else a zero AngularFlux —
            # paired NATIVE with a zero ψ_B on a carrying mesh (B.2d).
            if windowed:
                initial_guess = _windowed_cold_start(
                    self.scattering_op, self.sn_mesh, history_depth=2,
                )
            else:
                cold = _unwindowed_cold_start(self.sn_mesh, history_depth=2)
                initial_guess = (
                    _coupled_flux_state(cold, self.sn_mesh) if coupled else cold
                )

        # #282 route (a): on a carrying mesh the coupled rhs pairs the
        # 2-block source with the q½ fold of the (isotropic) fission
        # source — the TRUE starting-direction source the joint solve
        # consumes (System B's member, never a composite block — B.2d).
        q_driver = (
            _coupled_source_state(
                q_ext_composite,
                _radial_characteristic_fission_seed(
                    fission_source, self.sn_mesh,
                ),
                self.sn_mesh,
                context="SNSolver._solve_source_iteration",
            )
            if coupled else q_ext_composite
        )
        psi_typed, _residuals = si.solve(
            q_driver, initial_guess=initial_guess,
        )
        # Phase 3 measurement seam: accumulate this outer step's inner SI
        # iterate count (the eigenvalue path's only window onto the inner
        # spectral rate — see IterationHistory.total_inner_iterations).
        self._total_inner_iterations += len(_residuals)
        self._psi_typed = psi_typed

        # Scalar flux for the eigenvalue outer's contract.  Windowed: the
        # ℓ=0 moment IS the scalar flux (Y_0^0 = 1 ⇒ bit-identical to
        # integrate_angular) via the typed ``scalar_flux`` accessor that
        # carries the convention; un-windowed: reduce the full angular bulk.
        # The isinstance parses reify the driver-template contract (the
        # solve echoes the initial_guess representation) — the static
        # iterate type is the operators' ``FullField`` carrier, so the
        # bulk's representation re-narrows here, loudly on mismatch.
        from orpheus.transport.fields.harmonic_moment_flux import (
            HarmonicMomentFlux,
        )

        bulk = _system_a_member(psi_typed).interior
        if windowed:
            if not isinstance(bulk, HarmonicMomentFlux):
                raise TypeError(
                    f"windowed SI iterate must carry a HarmonicMomentFlux "
                    f"bulk (the moment template); got {type(bulk).__name__}."
                )
            phi = bulk.scalar_flux().values
        else:
            if not isinstance(bulk, AngularFlux):
                raise TypeError(
                    f"un-windowed SI iterate must carry an AngularFlux bulk "
                    f"(the flux template); got {type(bulk).__name__}."
                )
            phi = bulk.integrate_angular().values
        # The scale partner of the stored boundary trace (#291): the
        # leakage term of ``compute_keff`` rescales the trace's net
        # current by the fission-production ratio of the estimator's
        # flux to THIS flux (exactly 1.0 when the caller passes it back).
        self._phi_of_trace = phi
        return phi

    # ── Inner solver: Krylov on (L+C-S)·ψ = q_ext (R-1 Step D carve) ──

    def _solve_krylov(
        self, fission_source: np.ndarray, flux_distribution: np.ndarray,
    ) -> np.ndarray:
        r"""Inner solve via GMRES on :math:`(L + C - S)\,\psi = q_{\rm ext}`.

        R-1 Step D (2026-05-19) — carved onto
        :class:`~orpheus.numerics.iteration.KrylovAcceleration` consuming
        the operator triple

        .. math::

            A \;=\; L + C\,, \quad
            S \;=\; \text{full multi-group scatter}\,, \quad
            F \;=\; 0_{\rm wg}

        on typed :class:`~orpheus.sn.angular_flux.AngularFlux`.  The
        composite ``L + C`` returns an
        :class:`~orpheus.sn.operators.streaming.InvertibleOperator` (R-1 Step C);
        its ``.solve`` IS the WDD sweep but R-1 ships GMRES
        UNPRECONDITIONED (``preconditioner=None``) per user direction
        ("consolidating the foundational architecture; the block-inverse
        face preconditioner is `issue #200
        <https://github.com/deOliveira-R/ORPHEUS/issues/200>`_").  The
        sweep-as-preconditioner reactivation lives there.

        The within-group fission ``F`` is zero — the fission source
        comes in as the external ``q_{\rm ext}`` per the eigenvalue
        outer/within-group inner decomposition (Lewis & Miller §6.4).

        Scope
        =====

        D-H.2-C4d (2026-05-29) — 2-D Cartesian unblocked.  The L2
        :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
        face layout is the natural 4-face descriptor (xmin / xmax /
        ymin / ymax) that the legacy 1-D B1'' face block lacked; the
        L2-native representation ``loss_action`` 2-D matvec walk
        (S6.3 moved it off the operator; ``ScanMarch`` default since
        S6.9) operates on it directly.

        Returns the updated scalar flux ``(ng, nx, ny)``.
        """

        from orpheus.transport.fields.angular_flux import (
            AngularFlux,
        )
        from orpheus.transport.fields.angular_boundary_flux import (
            AngularBoundaryFlux,
        )
        from orpheus.transport.source_sinks import (
            AngularSourceSink,
            AngularBoundarySourceSink,
        )
        from orpheus.transport.timed_full_field import TimedFullField

        # ── Build the composite RHS ─────────────────────────────────
        # R-1 Step 4 A1 — q_ext as per-ordinate density via the canonical
        # ``AngularSourceSink.from_isotropic`` factory.  /W lives at the
        # factory boundary (Pattern 7 producer-side normalisation).
        # D-H.1c stage 2 — TimedFullField bulk + zero boundary (the
        # Krylov path does NOT pre-seed q_ext.boundary; reflective-BC
        # state threads through ``initial_guess`` per the audit §5
        # contract).
        q_ext_per_ord = AngularSourceSink.from_isotropic(
            fission_source, self.sn_mesh,
        )
        q_ext_composite = TimedFullField(
            # B.5.2: q_ext IS a source — carry the AngularSourceSink bulk AND
            # the AngularBoundarySourceSink inflow trace (zero for vacuum/reflective;
            # prescribed inflow otherwise). The SI rhs q_ext + S.apply + B.apply
            # closes on AngularBoundarySourceSink (operator outputs are sources).
            interior=q_ext_per_ord,
            boundary=AngularBoundarySourceSink.zeros_on(self.sn_mesh),
            _history=(),
            history_depth=2,
        )

        # ── Build the within-group system (single source of truth —
        # :func:`~orpheus.sn.coupled_system.build_within_group_system` /
        # ``_within_group_krylov``; shared with the SI and fixed-source
        # paths; the cached scattering operator injects through the cache
        # seam). ──────────────────────────────────────────────────────
        system = build_within_group_system(
            self.sn_mesh, self.mat_xs, scattering_op=self.scattering_op,
        )
        coupled = isinstance(system.resolvent, CoupledInvertibleOperator)

        # ── Warm start (composite / coupled pair) — built BEFORE the
        # driver so the GMRES restart is sized from the FULL ravel. ───
        # Post-B.2d ``self._psi_typed`` is a CoupledField on a carrying
        # mesh and a TimedFullField elsewhere; the Krylov ravellable
        # protocol detects either via ``to_flat`` / ``from_flat``
        # (D-H.1b.1) and threads it through the matvec / unravel cycle
        # natively.
        initial_guess = getattr(self, "_psi_typed", None)
        if initial_guess is None:
            # B.5.2: cold-start iterate is a FLUX composite, decoupled from
            # q_ext's now-AngularSourceSink type.  x0 stays all-zeros
            # (bit-identical); the flux template fixes the Krylov return
            # type — paired NATIVE with a zero ψ_B on a carrying mesh (B.2d).
            cold = TimedFullField.zeros(
                interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=self.sn_mesh,
            )
            initial_guess = (
                _coupled_flux_state(cold, self.sn_mesh) if coupled else cold
            )

        # ERR-053 (#282 route (a)): ``restart`` MUST cover the FULL ravel —
        # bulk ⊕ trace on System A plus BOTH System-B legs on the coupled
        # pair (``CoupledField.to_flat`` concatenates the systems, so the
        # sizing tracks automatically — the B.2a conformance closure; the
        # count is HONEST since B.2d, no dead padding).  A bulk-sized
        # restart re-truncates GMRES on the trace+seed DOFs (the sphere
        # Krylov stall).  Size it from the state the driver ravels.
        krylov = _within_group_krylov(
            system.resolvent, *system.gains,
            n_dof=int(initial_guess.to_flat().size),
            max_iter=self.max_inner, tol=self.inner_tol,
        )

        q_driver = (
            _coupled_source_state(
                q_ext_composite,
                _radial_characteristic_fission_seed(
                    fission_source, self.sn_mesh,
                ),
                self.sn_mesh,
                context="SNSolver._solve_krylov",
            )
            if coupled else q_ext_composite
        )
        psi_typed, _residuals = krylov.solve(
            q_driver, initial_guess=initial_guess,
        )
        # Phase 3 measurement seam: accumulate this outer step's inner
        # Krylov iterate count (see IterationHistory.total_inner_iterations).
        self._total_inner_iterations += len(_residuals)
        self._psi_typed = psi_typed

        # Reduce angular → scalar flux for the eigenvalue outer's contract.
        # The parse reifies the driver-template contract (the solve echoes
        # the flux initial_guess) loudly on mismatch.
        bulk = _system_a_member(psi_typed).interior
        if not isinstance(bulk, AngularFlux):
            raise TypeError(
                f"eigenvalue Krylov iterate must carry an AngularFlux bulk "
                f"(the flux template); got {type(bulk).__name__}."
            )
        # The scale partner of the stored boundary trace (#291) — see the
        # SI path's twin store.
        phi = bulk.integrate_angular().values
        self._phi_of_trace = phi
        return phi

    # D-J (2026-05-29): ``_make_sweep_preconditioner`` retired —
    # the method built a scipy LinearOperator wrapping the sweep as
    # the legacy packed-vector GMRES preconditioner.  Zero callers
    # remained post-D-K.1: the production Krylov path
    # (``_solve_fixed_source_krylov``) wraps the L+C composite in a
    # KrylovAcceleration consumer that handles the to_flat/from_flat
    # ravellable protocol on :class:`TimedFullField` natively.

    # ── Source computation helpers ────────────────────────────────────
    #
    # Wave D Issue 13: the math has been lifted into ScatteringOperator
    # (orpheus/sn/scattering.py). The methods below are thin delegators
    # preserved for the EigenvalueSolver Protocol surface and the
    # underscore-prefixed test probes in tests/sn/test_solver_components.py.
    # All four delegate to the same precomputed-by-construction operator
    # held on self.scattering_op, so bit-identity is by construction.

    def _add_scattering_source(self, Q: np.ndarray, phi: np.ndarray) -> None:
        """Add P0 scattering source to Q in-place (delegates to ScatteringOperator).

        Issue #196 PR-INDEX-5: every solver-public flux is principled
        ``(ng, nx, ny)``; the PR-INDEX-4 transpose pair retired.
        """
        self.scattering_op.add_iso_source(Q, phi)

    def _build_aniso_scattering(
        self, angular_flux: np.ndarray | None,
    ) -> np.ndarray | None:
        """Build per-ordinate anisotropic Pℓ scattering source (delegates to ScatteringOperator).

        Issue #196 PR-INDEX-5: ``angular_flux`` is principled
        ``(N, ng, nx, ny)``; return shape ``(N, ng, nx, ny)`` (or
        ``None`` when ``L == 0``).

        D-I.2 (2026-05-29): the underlying
        :meth:`ScatteringOperator.build_aniso_source` retired its
        bare-ndarray arm and now accepts only typed
        :class:`AngularFlux`.  This delegator wraps the inbound
        ``angular_flux`` ndarray as :class:`AngularFlux` at the helper
        boundary, then unwraps the resulting
        :class:`AngularSourceSink` via ``.values`` so the bare-ndarray
        external contract is preserved for legacy consumers in
        :func:`_solve_fixed_source_si` and the verification probes in
        :mod:`tests.sn.test_scattering_operator`.
        """
        if angular_flux is None:
            return None
        from orpheus.transport.fields.angular_flux import AngularFlux
        typed = AngularFlux.from_mesh(angular_flux, self.sn_mesh)
        result = self.scattering_op.build_aniso_source(typed)
        if result is None:
            return None
        return result.values

    def _add_n2n_source(self, Q: np.ndarray, phi: np.ndarray) -> None:
        """Add (n,2n) source to Q in-place (delegates to ScatteringOperator).

        Issue #196 PR-INDEX-5: both arguments are principled
        ``(ng, nx, ny)``.
        """
        self.scattering_op.add_n2n_source(Q, phi)


# ═══════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════

# ─────────────────────────────────────────────────────────────────────
# Retired in P1.7 (moment-space + layering plan):
#   `_build_rhs_cartesian`, `_build_rhs_spherical`,
#   `_build_rhs_cylindrical`  (the latter was an alias of spherical)
#
# All three were packed-1D RHS builders for the legacy BiCGSTAB
# FD-operator path.  Production call graph (verified by
# `tests/sn/test_fixed_source_g1.py::TestNoLegacyMachineryInCallPath::
# test_no_legacy_eq_map_or_decoder_in_g1_path`): ZERO callers — the
# G1 Krylov path goes through typed `KrylovAcceleration` +
# `(L + C).apply` instead of building a packed-1D RHS.
#
# The cartesian variant carried an inline `(2*l+1)` literal at its
# Pℓ-source build, duplicating the canonical addition-theorem factor
# now sourced from `SphericalHarmonicBasis.addition_theorem_factor`
# (the frame's reconstruction face reads it live).  Per the moment-space
# plan §P1.3 "exactly one place" claim, retiring this dead duplicate is
# required.
#
# The spherical/cylindrical variants had no (2*l+1) duplicate (they
# were P0-isotropic-only), but they share the same dead-code status
# and retire alongside cartesian for cleanliness — superseded code is
# noise per `feedback_aggressive_retirement`.
# ─────────────────────────────────────────────────────────────────────


def _is_curvilinear(mesh: Mesh1D | Mesh2D) -> bool:
    """Return ``True`` iff *mesh* is a 1-D spherical or cylindrical Mesh1D."""
    if not isinstance(mesh, Mesh1D):
        return False
    coord = getattr(mesh, "coord", None)
    if coord is None:
        return False
    name = getattr(coord, "name", str(coord)).upper()
    return name in ("SPHERICAL", "CYLINDRICAL")


def _reflect_outflow_into_inflow(
    boundary_flux, sn_mesh: SNMesh, faces: "Iterable[str] | None" = None,
    *,
    radial_characteristic: "RadialCharacteristicFlux | None" = None,
) -> None:
    r"""In-place: fill each face's inflow ordinate slots with the realized
    boundary law applied to that face's outflow trace — the ``−B`` reflective
    coupling, externalised for the bare ``transport_sweep`` (Wave O #208 O.4a.2).

    The bare sweep reads the inflow ordinate slots of its boundary buffer as
    the inflow seed; it no longer re-applies ``bc`` to the outflow internally.
    The SI driver path supplies ``B·ψ.outflow`` through ``rhs.boundary`` (as the
    ``B`` coupling gain), but the DIRECT fixed-source SI loop
    (:func:`_solve_fixed_source_si`) and the final eigenvalue reconstruction
    sweep (:func:`solve_sn`) do not route through that driver — they call this
    helper to set ``ψ.inflow = B·ψ.outflow`` on the buffer before each sweep.
    Post-un-weld ``B = B_a + B_b`` is the direct sum of the two per-system
    boundaries (RULING P1), so this helper reflects each system through its OWN
    operator: the trace via
    :meth:`~orpheus.sn.operators.boundary.SNBoundaryOperator.reflect_inflow_inplace`
    (``B_a``) and the ψ½ ray corner via
    :meth:`~orpheus.sn.operators.boundary.RadialCharacteristicBoundaryOperator.reflect_corner_inplace`
    (``B_b``) — the SAME two cores the matvec / SI driver consume as the
    ``B_a + B_b`` gain (single source of truth; the two routes cannot drift).

    For vacuum ``B = 0`` so the inflow slots stay zero (bit-identical to the
    pre-extraction ``bc.apply`` of a vacuum law); for reflective/white/albedo
    it is the same ``R·G`` reflection the pre-extraction sweep applied at entry,
    merely relocated to the caller.

    This is the TWIN of the driver route: the within-group SI/Krylov drivers
    deliver the SAME ``−B`` coupling via the SAME
    :class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`, as a first-class
    coupling gain (Wave O O.2a — since B.2d the
    :func:`~orpheus.sn.coupled_system.build_within_group_system` record's
    ``gains``).  The two differ only in plumbing: this helper writes the
    buffer's inflow slots directly; the driver's ``B`` gain rides
    ``B·ψ.outflow`` in ``rhs.boundary``.  The DRIVER route no longer needs this
    helper (O.2a collapsed it); it survives for the final eigenvalue
    reconstruction sweep (which has no driver to route through) AND Phase 3's
    octant-group Gauss-Seidel resolvent (which calls it face-restricted between
    octant-group sweeps — see ``faces``).

    ``faces`` (Phase 3 G-S): ``None`` (default) reflects EVERY boundary face —
    the whole-trace ``−B`` used by the reconstruction sweep + the SI seed. A
    face subset restricts the reflect to those faces' inflow slots, leaving the
    rest untouched: the G-S resolvent absorbs a just-swept group's outgoing
    reflective faces into ``boundary_flux``, then calls this to reflect ONLY
    those faces' outflow into inflow, so the next group reads the fresh
    current-iterate inflow. ``B`` is block-diagonal over faces ⟹ exact
    restriction.
    """
    from orpheus.sn.operators.boundary import (
        RadialCharacteristicBoundaryOperator,
        SNBoundaryOperator,
    )

    # Trace-only ``A_ss`` action — no zero-bulk probe (the bulk was only ever
    # a carrier to reach ``B``'s boundary block).  The mutating write-back is
    # ``B_a``'s own :meth:`reflect_inflow_inplace` verb (#226 step 2 moved it
    # onto the operator so the scheduled sweep's inter-group reflect and this
    # helper share ONE body); it routes through ``_reflect_trace`` with
    # ``B_a.apply``, so the helper and the matvec / SI driver cannot drift.
    SNBoundaryOperator(sn_mesh).reflect_inflow_inplace(boundary_flux, faces=faces)
    # ``radial_characteristic`` (#282 route (a)): the ψ½ carrier whose
    # inflow-corner slots get the law's corner action — the seed analogue,
    # through System B's OWN boundary ``B_b`` (RULING P1: one reflect per
    # system; ``B_b.reflect_corner_inplace`` is the ray sibling of ``B_a``'s
    # trace reflect above). Presence-guarded at the call site since B.2b:
    # a seedless mesh has no System B, so B_b is unconstructable there (the
    # old in-method ``None`` no-op moved here — a None ray ⟺ seedless).
    if radial_characteristic is not None:
        RadialCharacteristicBoundaryOperator(sn_mesh).reflect_corner_inplace(
            radial_characteristic,
        )


# ═══════════════════════════════════════════════════════════════════════
# Public API
# ═══════════════════════════════════════════════════════════════════════

def solve_sn(
    materials: dict[int, Mixture],
    mesh: "Mesh1D | Mesh2D | tuple[Axis1D, ...]",
    quadrature: Quadrature,
    inner_solver: str = "source_iteration",
    scattering_order: int = 0,
    max_outer: int = 500,
    keff_tol: float = 1e-7,
    flux_tol: float = 1e-6,
    max_inner: int = 200,
    inner_tol: float = 1e-8,
    inner_schedule: str = "jacobi",
    mat_map: "np.ndarray | None" = None,
) -> Solution:
    """Solve the multi-group SN eigenvalue problem.

    This is the **canonical entry point** for the production SN solver.
    Production callers consume ``(materials, mesh, quadrature, ...)``
    directly: materials are :class:`~orpheus.data.macro_xs.mixture.Mixture`
    objects keyed by material ID, ``mesh`` is a
    :class:`~orpheus.geometry.Mesh1D` / :class:`~orpheus.geometry.Mesh2D`
    (build via :meth:`Mesh1D.from_geometry` for multi-region 1-D cases)
    OR an axis tuple — the axis-native surface and the ONLY 3-D entry
    (C5.5, #225; per-axis BCs ride the axes, ``mat_map=`` carries the
    material assignment), and ``quadrature`` is an explicitly chosen
    :class:`~orpheus.numerics.quadrature.Quadrature` — Gauss-Legendre
    for slab, level-symmetric / product quadrature for curvilinear, or
    Lebedev for 2-D.

    The mesh's boundary conditions (``bc_left`` / ``bc_right`` for 1-D,
    ``bc_xmin`` / ``bc_xmax`` / ``bc_ymin`` / ``bc_ymax`` for 2-D) are
    honoured verbatim — the SN sweep handles ``vacuum`` and
    ``reflective``.

    Parameters
    ----------
    materials : dict mapping material ID to Mixture.
    mesh : Mesh1D or Mesh2D (base geometry).
    quadrature : Quadrature
        Explicitly chosen by the caller — Gauss-Legendre for slab,
        level-symmetric / product quadrature for curvilinear 1-D,
        Lebedev for 2-D. Mismatches between geometry and quadrature
        family are not silently coerced.
    inner_solver : "source_iteration" (default) or "krylov".
        Wave E Round 2 deviation from the campaign plan: ``solve_sn``
        keeps the ``source_iteration`` default for **all** geometries
        (Cartesian and curvilinear).  The eigenvalue is shape-
        independent (k = production / absorption is a volume-weighted
        ratio), so even on the ERR-026-affected curvilinear sweep the
        keff is correct to the eigenvalue's tolerance, even though the
        flux *shape* would have a closure-bias drift.  Preserving the
        default keeps the 6 curvilinear regression snapshots bit-
        identical.  ``solve_sn_fixed_source`` *does* auto-flip to
        ``"krylov"`` for curvilinear because shape correctness is the
        whole point of fixed-source MMS verification.
    scattering_order : Legendre order for scattering (0 = P0).
    max_outer : maximum outer (power) iterations.
    keff_tol, flux_tol : outer convergence.
    max_inner, inner_tol : inner solver parameters.
    inner_schedule : {"jacobi", "gauss_seidel"}
        Boundary splitting for the ``source_iteration`` inner (#218 — the
        eigenvalue inner CAN now reach the same boundary-G-S accelerator the
        fixed-source path got in Phase 3, via the shared ``_within_group_si``
        builder; validated SI(G-S)≡Krylov≡closed-form k_inf).  ``"jacobi"``
        (default) lags the reflective boundary ``B`` as an external gain;
        ``"gauss_seidel"`` (opt-in) folds ``B`` into an octant-group forward
        substitution on 2-D Cartesian (1-D / curvilinear auto-fall-back to
        Jacobi — G-S is 2-D-Cartesian-only and a no-op on the
        scattering-dominated 1-D regime).  The converged eigenvalue is
        identical either way to within ``inner_tol`` (vv-principles Mode 9 —
        same fixed point; only the inner SI stopping differs).  Default stays
        Jacobi: the eigenvalue inner is warm-started (modest G-S benefit) and a
        schedule change shifts k_eff by ~inner_tol, which the keff_tol-tight
        regression snapshots cannot absorb.  Ignored when
        ``inner_solver="krylov"`` (Krylov is splitting-invariant).

    Returns
    -------
    Solution
        Typed return carrying eigenvalue, typed
        :class:`~orpheus.sn.angular_flux.AngularFlux` +
        :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux` +
        :class:`~orpheus.sn.boundary_flux.AngularBoundaryFlux` fields plus an
        :class:`~orpheus.sn.solution.IterationHistory` carrying the
        eigenvalue trajectory.  Issue #197 PR-TYPED-5 — the legacy
        :class:`SNResult` data bag was retired in favour of the
        unified :class:`Solution` type that covers both eigenvalue and
        fixed-source problems.
    """
    t_start = time.perf_counter()

    # Build augmented geometry (precomputes streaming stencil).
    # Issue #197 PR-TYPED-0: materials now lives on SNMesh — the
    # phase-space-as-such object. C5.5 (#225): the declaration may be a
    # legacy mesh or an axis tuple (the only 3-D entry); unset faces
    # resolve to the SNMesh reflective default (eigenvalue convention).
    sn_mesh = _as_sn_mesh(mesh, quadrature, materials, mat_map=mat_map)

    solver = SNSolver(
        sn_mesh,
        inner_solver=inner_solver,
        scattering_order=scattering_order,
        keff_tol=keff_tol, flux_tol=flux_tol,
        max_inner=max_inner, inner_tol=inner_tol,
        inner_schedule=inner_schedule,
    )

    keff, keff_history, scalar_flux = power_iteration(solver, max_iter=max_outer)

    # Final sweep to get angular flux.  Issue #196 PR-INDEX-5: every
    # array is principled — scalar_flux ``(ng, nx, ny)``, angular_flux
    # ``(N, ng, nx, ny)``.
    # R-1 Step 4 A1: typed source via the canonical
    # :meth:`AngularSourceSink.from_isotropic` factory (Pattern 7
    # producer-side normalisation — /W projection at the factory
    # boundary).
    from orpheus.transport.source_sinks import AngularSourceSink
    from orpheus.transport.fields.angular_boundary_flux import (
        AngularBoundaryFlux as _BoundaryFlux,
    )
    Q_final = solver.compute_fission_source(scalar_flux, keff)
    solver._add_scattering_source(Q_final, scalar_flux)
    solver._add_n2n_source(Q_final, scalar_flux)
    # Wave O #208 O.4a.2 — BARE final reconstruction sweep: seed its inflow
    # from the CONVERGED boundary trace.  The inner solves drive the inflow as
    # a live unknown in ``ψ.boundary`` (carried on ``solver._psi_typed`` from
    # the last inner solve), so ``solver._boundary_flux`` is no longer the
    # partner-flux carrier (it stays all-zeros).  Reflect the converged outflow
    # into the inflow slots via −B (no-op for vacuum; idempotent here since the
    # converged inflow already equals ``B·ψ.outflow``), then sweep.
    converged = getattr(solver, "_psi_typed", None)
    # B.2d: the converged iterate is the coupled pair on a carrying mesh;
    # System A carries the trace, System B the converged ψ½ state.
    converged_a = _system_a_member(converged) if converged is not None else None
    final_boundary = (
        converged_a.boundary if converged_a is not None
        else _BoundaryFlux.zeros_on(sn_mesh)
    )
    # #289-F2 role parse: the reconstruction sweep seeds its inflow from the
    # ANGULAR trace; the widened composite slot erases the family.
    if not isinstance(final_boundary, _BoundaryFlux):
        raise TypeError(
            f"eigenvalue finalize: the converged trace must be an "
            f"AngularBoundaryFlux; got {type(final_boundary).__name__}."
        )
    # #282 route (a): the reconstruction sweep's ψ½ pair on a carrying
    # mesh — q½ = the ℓ = 0 fold of the converged total source; the
    # carrier pre-loads the converged outflow corner so the corner
    # reflect below can seed the inflow corner (vacuum ⇒ 0).
    q_final_per_ord = AngularSourceSink.from_isotropic(Q_final, sn_mesh)
    final_seed_src = _radial_characteristic_source_from_per_ordinate(
        q_final_per_ord.values, sn_mesh,
    )
    final_seed_buf = _radial_characteristic_zeros(sn_mesh)
    converged_ray = (
        _system_b_member(converged) if converged is not None else None
    )
    if final_seed_buf is not None and converged_ray is not None:
        final_seed_buf.values[...] = converged_ray.to_unified().values
    # Wave O #208 O.4b Phase E2 — the 2-D wavefront is now BARE (reads the
    # given inflow, no in-sweep bc.apply), so the reflective coupling is the
    # external -B for 2-D too.  The guard is lifted: _reflect_outflow_into_inflow
    # is geometry-agnostic (iterates boundary_flux.layout.faces via the canonical
    # SNBoundaryOperator — verified 2-D-ready) and idempotent here (the converged
    # inflow already equals B·ψ.outflow); vacuum stays a no-op (B = 0).  The
    # seed carrier's inflow corner rides the same reflect (#282 route (a)).
    _reflect_outflow_into_inflow(
        final_boundary, sn_mesh, radial_characteristic=final_seed_buf,
    )
    # The corner datum is GIVEN data for the sweep: thread it into the q½
    # source's corner slot (the sweep reads the SOURCE corner as the
    # inward march's entry — the trace's seeded-inflow discipline).
    if final_seed_src is not None and final_seed_buf is not None:
        space = final_seed_src.space
        for _p in space.levels:
            space.corner_view(final_seed_src.values, _p, -1)[...] = (
                space.corner_view(final_seed_buf.values, _p, -1)
            )
    angular_flux, _ = transport_sweep(
        q_final_per_ord,
        solver.mat_xs.total_cross_section, sn_mesh,
        final_boundary,
        radial_characteristic_source=final_seed_src,
        radial_characteristic_flux=final_seed_buf,
    )

    elapsed = time.perf_counter() - t_start

    # Issue #197 PR-TYPED-5: build typed Solution.  The bare ndarrays
    # produced by power_iteration + the final sweep get wrapped in
    # typed fields at this single boundary; downstream consumers read
    # typed fields.
    # D-H.1c stage 2 (2026-05-28): Solution.angular_flux is a
    # TimedFullField composite (bulk + boundary + history) constructed
    # DIRECTLY from L2 types.  The legacy-AngularFlux adapter wrap is
    # GONE — ``angular_flux`` is a bare ndarray from the final sweep;
    # we wrap it once into the composite carrier.
    from orpheus.transport.fields.angular_flux import (
        AngularFlux,
    )
    from orpheus.transport.fields.angular_boundary_flux import (
        AngularBoundaryFlux,
    )
    from orpheus.transport.fields.scalar_flux import ScalarFlux
    from orpheus.transport.timed_full_field import TimedFullField
    history = IterationHistory(
        keff_history=tuple(keff_history),
        n_outer=len(keff_history),
        total_inner_iterations=solver._total_inner_iterations,
        converged=True,
    )
    return Solution(
        angular_flux=TimedFullField(
            interior=AngularFlux.from_mesh(angular_flux, sn_mesh),
            # Wave O #208 O.4a.2: the converged boundary trace from the final
            # bare sweep (inflow = B·ψ.outflow, outflow = streamed).
            boundary=final_boundary,
            _history=(),
            history_depth=2,
        ),
        scalar_flux=ScalarFlux.from_mesh(scalar_flux, sn_mesh),
        mesh=sn_mesh,
        keff=float(keff_history[-1]),
        history=history,
        # B.2d DP-Solution: System B's converged state is its OWN member —
        # the final sweep's marched ψ½ carrier, lifted role-preserving onto
        # the composite (None on non-carrying meshes).
        radial_characteristic=(
            None
            if final_seed_buf is None
            else RadialCharacteristicComposite.from_unified(final_seed_buf)
        ),
    )


def _build_fixed_source_rhs(
    external_source: "np.ndarray | TimedFullField | CoupledField",
    sn_mesh: SNMesh,
) -> "TimedFullField | CoupledField":
    r"""Normalize the external source into the driver RHS.

    A fixed-source problem's RHS is the composite source
    :class:`~orpheus.transport.timed_full_field.TimedFullField` — a bulk
    :class:`~orpheus.transport.source_sinks.AngularSourceSink` paired with a
    boundary :class:`~orpheus.transport.source_sinks.AngularBoundarySourceSink`
    (the prescribed inflow :math:`q` of the affine BC
    :math:`\gamma_-\psi = R\,G\,\gamma_+\psi + q`) — paired, on a
    carrying mesh (B.2d), with System B's q½ member into the coupled
    ``CoupledField[q_A, q_B]``. This is the ONE object that represents a
    source everywhere in the solve; this helper is its single construction
    point (Cardinal Rule 2 — the SI and Krylov inner paths both consume
    what it returns, rather than each re-deriving it).

    ``external_source`` is accepted in three forms (the bulk array is a typed
    union of TWO ndarray ranks; see :func:`_lift_external_source_to_moments`):

    * **flat** ``np.ndarray`` of shape ``(N, ng, *spatial)`` — the
      per-ordinate-density BULK source only; the boundary is vacuum (all-zero).
      The original form; every pre-existing caller keeps working unchanged
      (the slope moments :math:`\hat Q` are zeroed by the lift — the honest
      default, exact for a region-uniform source).
    * **moment-resolved** ``np.ndarray`` of shape
      ``(N, ng, *spatial, per_axis**ndim)`` — a multi-moment closure (LD)
      external source whose trailing axis carries the per-cell tensor-Legendre
      moment vector (slot 0 = cell average, the slope rows = :math:`\hat Q`,
      d=2 Kronecker order ``[Q̄, Q̂_y, Q̂_x, Q̂_xy]``).  The CALLER projects
      :math:`Q^{\rm ext}` onto the moment vector (e.g. by Gauss quadrature
      ``∫q·Pₖ``); this entry threads the slope rows through unchanged so they
      join the moment-carrying scattering source ``Σ_s·φ̂`` in the SI rhs
      (#247 — the slope-SOURCE half of the LM-1989 trap).  Only meaningful for
      LD (``per_axis > 1``); a moment-resolved input on a DD/Step mesh
      (``per_axis == 1``, no moment axis) is rejected by the flat-shape check.
    * :class:`TimedFullField` — the full COMPOSITE source (bulk + a possibly
      non-vacuum prescribed-inflow boundary, e.g. from
      :meth:`AngularBoundarySourceSink.prescribed_inflow`). Its leaf values are
      re-homed onto ``sn_mesh``: the trace/grid layout is deterministic from
      ``(mesh, quadrature, materials)``, so this is an exact values-copy onto
      the solve's own mesh instance — required because the within-group
      operators are built on ``sn_mesh`` and :class:`TimedFullField` algebra
      enforces mesh identity.  Its bulk may be flat OR moment-resolved.
    * :class:`~orpheus.numerics.coupled_system.CoupledField` — the coupled
      pair with an EXPLICIT q½ member (B.2d): System A re-homes as above;
      System B's member values re-home onto the fold slot exactly (the
      caller controls the true q½ instead of the ℓ = 0 fold default).
    """
    from orpheus.transport.source_sinks import (
        AngularSourceSink,
        AngularBoundarySourceSink,
    )

    N = sn_mesh.quad.N
    ng = sn_mesh.ng
    expected = (N, ng, *sn_mesh.spatial_shape)

    explicit_seed: "RadialCharacteristicComposite | None" = None
    source_a: "np.ndarray | FullField" = (
        external_source if not isinstance(external_source, CoupledField)
        else _system_a_member(external_source)
    )
    if isinstance(external_source, CoupledField):
        # The coupled pair: System B's explicit q½ member is captured for
        # the seed slot below; System A continues through the composite arm.
        explicit_seed = _system_b_member(external_source)
    if isinstance(source_a, FullField):
        bulk_values = np.asarray(source_a.interior.values)
        trace_size = int(sn_mesh.angular_trace.layout.total_size)
        boundary_values = source_a.boundary.values
        if boundary_values.size != trace_size:
            raise ValueError(
                f"_build_fixed_source_rhs: composite boundary source has "
                f"{boundary_values.size} values, but sn_mesh.angular_trace expects "
                f"{trace_size} (layout mismatch — the composite must be built "
                f"on the same mesh / quadrature / materials)."
            )
        boundary = AngularBoundarySourceSink.from_mesh(boundary_values.copy(), sn_mesh)
    else:
        bulk_values = np.asarray(source_a)
        if bulk_values.dtype == object:
            # A stray non-array, non-TimedFullField object (e.g. a bare
            # AngularSourceSink) — np.asarray wraps it as a 0-d object array.
            # Reject explicitly rather than failing the shape check obscurely.
            raise TypeError(
                f"_build_fixed_source_rhs: external_source must be an "
                f"(N, ng, *spatial) array (bulk-only / vacuum) or a "
                f"TimedFullField composite source; got "
                f"{type(external_source).__name__}"
            )
        boundary = AngularBoundarySourceSink.zeros_on(sn_mesh)

    # Issue #196 PR-INDEX-5 + #247: the bulk source is a typed union of TWO
    # principled ndarray ranks — flat ``(N, ng, *spatial)`` (the original path)
    # OR moment-resolved ``(N, ng, *spatial, per_axis**ndim)`` (LD only; the new
    # slope-SOURCE path).  Discriminate by RANK (NOT trailing-size — a
    # coincidental spatial dim could equal 2^d): a flat bulk has exactly
    # ``len(expected)`` axes; a moment-resolved bulk has ONE more (the trailing
    # 2^d moment axis).  Reject everything else, INCLUDING a moment axis whose
    # length ≠ per_axis**ndim, and (for DD/Step where there is no moment axis) a
    # moment-resolved input outright (only flat is valid at per_axis == 1).
    n_cell_moments = cell_moment_count(
        sn_mesh.scheme.spatial_basis_per_axis, sn_mesh.ndim
    )
    moment_expected = (*expected, n_cell_moments)
    is_flat = bulk_values.shape == expected
    is_moment_resolved = (
        n_cell_moments > 1 and bulk_values.shape == moment_expected
    )
    if not (is_flat or is_moment_resolved):
        if n_cell_moments > 1 and bulk_values.shape[:-1] == expected:
            # Right rank for a moment-resolved bulk, but the trailing moment
            # axis is the wrong width — name the expected 2^d so the relaxation
            # does not swallow a real shape bug (#247 negative pin).
            raise ValueError(
                f"fixed-source moment-resolved bulk shape {bulk_values.shape} "
                f"has trailing moment axis {bulk_values.shape[-1]}, expected "
                f"per_axis**ndim = {n_cell_moments} "
                f"(moment vector {moment_expected})"
            )
        raise ValueError(
            f"fixed-source bulk shape {bulk_values.shape} does not match "
            f"(N, ng, *spatial) = {expected}"
            + (
                f" or the moment-resolved {moment_expected}"
                if n_cell_moments > 1 else ""
            )
        )
    # The external source carries the trailing 2^d spatial-moment axis at a
    # multi-moment closure (#240 D5b-S3) so it composes with the moment-carrying
    # scattering source ``S.apply(ψ)`` in the SI rhs ``q_ext + S.apply(ψ)``.  A
    # FLAT external source is flat-in-moment (Q̂ = 0 — the slope rows are zero,
    # the honest default exact for a region-uniform source): lift onto slot 0
    # (average), rest zero.  A MOMENT-RESOLVED external source already carries
    # the slope rows Q̂ (the caller projected them — #247): thread them through
    # unchanged.  DD/Step (per_axis == 1) → no lift, byte-identical.
    bulk_values, per_axis = _lift_external_source_to_moments(bulk_values, sn_mesh)
    q_a = TimedFullField(
        interior=AngularSourceSink.from_mesh(
            bulk_values, sn_mesh, spatial_moments=per_axis,
        ),
        boundary=boundary,
    )
    if sn_mesh.radial_characteristic_space is None:
        return q_a
    # #282 route (a) → B.2d: the q½ member on carrying meshes — System B's
    # OWN composite, paired with q_A.  A coupled input's explicit member is
    # re-homed values-exactly (same deterministic-layout argument as the
    # trace copy above); otherwise the ℓ = 0 fold of the per-ordinate bulk
    # populates it.  Carrying meshes are 1-D curvilinear DD (never
    # moment-resolved), so the flat bulk is the only live shape here.
    if explicit_seed is not None:
        seed_src = RadialCharacteristicSourceSink.zeros_on(sn_mesh)
        seed_src.values[...] = explicit_seed.to_unified().values
    else:
        seed_src = _radial_characteristic_source_from_per_ordinate(
            bulk_values, sn_mesh,
        )
    return _coupled_source_state(
        q_a, seed_src, sn_mesh, context="_build_fixed_source_rhs",
    )


def _lift_external_source_to_moments(
    bulk_values: "np.ndarray", sn_mesh: SNMesh,
) -> "tuple[np.ndarray, int]":
    r"""Lift / thread an external source onto the ``2^d`` cell-moment vector,
    returning ``(lifted, per_axis)``.

    Single source of the external-source moment lift for the fixed-source path
    (#240 D5b-S3 / #247 — the slope-SOURCE widening).  One production caller
    (:func:`_build_fixed_source_rhs`); kept as a single-source helper so a future
    eigenvalue external-source hook reuses the same lift/thread policy.
    ``bulk_values`` is a typed union of TWO ndarray ranks, discriminated by RANK
    (NOT trailing-size —
    :func:`~orpheus.numerics.moment_layout.is_moment_valued_by_flat_rank` against
    the flat ``(N, ng, *spatial)`` rank; a coincidental spatial dim could equal
    ``2^d``):

    * **DD/Step** (``per_axis == 1``, ``tail == ()``): no moment axis — return
      the input unchanged (byte-identical, the backward-compat negative
      control).
    * **flat** ``(N, ng, *spatial)``: zero the ``2^d`` buffer, copy the flat
      input onto slot 0 (average), leave the slope rows :math:`\hat Q` ZERO —
      the honest default (``q̂ = 0`` is exact with no sub-cell information).
    * **moment-resolved** ``(N, ng, *spatial, 2^d)``: the caller already
      projected the slope rows (#247) — thread the moment vector through
      UNCHANGED (validate its trailing axis == ``2^d``).  Joins the
      moment-carrying scattering source ``Σ_s·φ̂`` in the SI rhs; the per-octant
      slope-sign reframe (``sweep_graph._CellSolve`` /
      ``octant_moment_frame_signs``) re-signs the external slopes global→sweep
      EXACTLY as it does the scattering slopes — no new cell branch."""
    per_axis = sn_mesh.scheme.spatial_basis_per_axis
    n_cell_moments = cell_moment_count(per_axis, sn_mesh.ndim)
    tail = face_moment_tail(n_cell_moments)
    if tail == ():
        return bulk_values, per_axis
    flat_ndim = 2 + len(sn_mesh.spatial_shape)  # (N, ng, *spatial)
    if not is_moment_valued_by_flat_rank(bulk_values, flat_ndim):
        # Flat input → lift onto slot 0, slopes zero (the honest default).
        lifted = np.zeros((*bulk_values.shape, *tail), dtype=bulk_values.dtype)
        lifted[..., AVERAGE_MOMENT] = bulk_values
        return lifted, per_axis
    # Moment-resolved input → thread the slope rows through unchanged (#247).
    if bulk_values.shape[-1] != n_cell_moments:
        raise ValueError(
            f"_lift_external_source_to_moments: moment-resolved bulk has "
            f"trailing moment axis {bulk_values.shape[-1]}, expected "
            f"per_axis**ndim = {n_cell_moments}"
        )
    return bulk_values, per_axis


def _average_moment_scalar(phi: "np.ndarray", sn_mesh: SNMesh) -> "np.ndarray":
    r"""Reduce a (possibly moment-carrying) scalar flux to its cell-AVERAGE.

    The user-facing :class:`Solution` scalar flux is the cell-average moment
    (slot 0); a multi-moment closure's φ̂ slopes are internal within-cell DG
    structure (#240 D5b-S3).  ``phi`` from a multi-moment closure carries a
    trailing ``2^d`` axis — take slot ``AVERAGE_MOMENT``; DD/Step (per_axis ==
    1) → no axis → return unchanged."""
    per_axis = sn_mesh.scheme.spatial_basis_per_axis
    if face_moment_tail(cell_moment_count(per_axis, sn_mesh.ndim)) == ():
        return phi
    return phi[..., AVERAGE_MOMENT]


def solve_sn_fixed_source(
    materials: dict[int, Mixture],
    mesh: "Mesh1D | Mesh2D | tuple[Axis1D, ...]",
    quadrature: Quadrature,
    external_source: "np.ndarray | TimedFullField",
    boundary_condition: str = "vacuum",
    scattering_order: int = 0,
    max_inner: int = 1000,
    inner_tol: float = 1e-12,
    inner_solver: str | None = None,
    inner_schedule: str = "gauss_seidel",
    mat_map: "np.ndarray | None" = None,
    scheme: "DiscretizationSchemeBase | None" = None,
) -> Solution:
    r"""Solve the multi-group SN fixed-source transport problem.

    Solves

    .. math::

        \mu_n \frac{\partial\psi_n}{\partial x}+\Sigma_t\psi_n
        = \frac{1}{W}\left(\Sigma_s\phi + Q^{\text{ext}}_n\right)

    for a prescribed per-ordinate external source ``external_source``,
    with vacuum or reflective boundary conditions. The fission source
    is zero — this is the pure transport operator. Source iteration
    converges geometrically at rate :math:`c = \max\Sigma_s/\Sigma_t`.

    Parameters
    ----------
    materials, mesh, quadrature, scattering_order :
        Same as :func:`solve_sn`.
    external_source : (N, ng, nx, ny) ndarray OR TimedFullField
        The fixed-source RHS, in either of two forms (normalised by
        :func:`_build_fixed_source_rhs`):

        * ``np.ndarray`` of shape ``(N, ng, nx, ny)`` — the per-ordinate
          volumetric BULK source :math:`Q^{\text{ext}}_n(x)` in
          **per-ordinate density magnitude** (R-1 Step 4 A1 convention),
          with a **vacuum** boundary. Callers with an iso scalar source
          :math:`Q(\vec r, g)` should project to per-ordinate via
          :meth:`~orpheus.transport.source_sinks.AngularSourceSink.from_isotropic`
          before passing (the :math:`1/W` projection lives at the producer
          boundary per Pattern 7). The sweep does NOT apply ``/W``
          internally. Issue #196 PR-INDEX-5: principled layout (``g`` axis
          after ``N``).
        * :class:`~orpheus.transport.timed_full_field.TimedFullField` — the
          full **composite** source ``q = q_bulk ⊕ q_∂`` (a bulk
          :class:`~orpheus.transport.source_sinks.AngularSourceSink` paired
          with a boundary
          :class:`~orpheus.transport.source_sinks.AngularBoundarySourceSink`). This
          is how a **non-vacuum prescribed inflow** is supplied — build the
          boundary via
          :meth:`~orpheus.transport.source_sinks.AngularBoundarySourceSink.prescribed_inflow`
          (the affine-BC inhomogeneous term :math:`q`, consumed by the sweep
          as the inflow seed). The legacy array form is exactly the
          bulk-only / vacuum special case of this composite.
    boundary_condition : {"vacuum", "reflective"}
        Applied to all faces when the mesh has no explicit BC
        declarations (``bc_left`` etc. are ``None``).  When the mesh
        carries explicit :class:`~orpheus.geometry.mesh.BC` fields,
        those take precedence and this parameter is ignored.
        Vacuum is the default because the intended consumer is
        Method of Manufactured Solutions verification on a finite slab.
    max_inner, inner_tol :
        Inner solver iteration limits.
    inner_solver : {"source_iteration", "krylov", None}
        Inner-solve strategy.  When ``None`` (default), all geometries
        use ``"source_iteration"``.

        History: from Phase D (2026-05-12) through the ERR-058 fix
        (2026-06-12, Issue #195) the curvilinear ``None`` default
        resolved to ``"krylov"``, because the SWEEP's fixed point was
        wrong on non-flat fields (the ERR-026/ERR-058 closure-seed
        family) while Krylov-on-apply tracked the (then-distinct)
        matvec system.  Post-unification the sweep and matvec are ONE
        discrete system, and the ERR-058 closure-seed fix (coupled-
        pole spatial seed + angular-edge-extrapolation half-angle
        seed) makes that system O(h²)-consistent — SI and Krylov now
        converge to bit-identical fixed points on the curvilinear MMS
        ladders, with SI ~10²× faster (no GMRES restart pathology,
        ERR-053).  The curvilinear default therefore returned to
        ``"source_iteration"``; ``"krylov"`` stays available as the
        opt-in cross-check (the SI ≡ Krylov fixed-point equivalence
        is itself a standing splitting-invariance gate, vv Mode 9).
    inner_schedule : {"gauss_seidel", "jacobi"}
        Source-iteration BOUNDARY splitting (Phase 3, ``inner_solver=
        "source_iteration"`` only).  ``"gauss_seidel"`` (default) folds the
        reflective coupling ``B`` into an octant-group Gauss-Seidel resolvent
        (2-D Cartesian) — re-reflecting each octant group's outgoing reflective
        faces between group sweeps so a later group reads the fresh
        current-iterate inflow (a modest reflective-SI rate gain, ~0.86–0.92×
        on B-mixture configs).  ``"jacobi"`` lags ``B`` fully (the
        splitting-invariant control).  The converged fixed point is IDENTICAL
        for both — this selects only the SI spectral rate.  1-D meshes always
        fall back to Jacobi (boundary G-S is a no-op on the scattering-
        dominated 1-D regime, and the scan is not a wavefront).  The dominant
        within-group SCATTERING rate is unchanged either way — that is Krylov /
        consistent-DSA territory (issue #2).

    Notes
    -----
    This entry point exists for L1 verification via MMS, not for
    engineering problems — real fixed-source calculations should still
    build on :func:`solve_sn` with an appropriate external-source hook.
    See :mod:`orpheus.derivations.continuous.mms.sn` and the MMS verification
    section of the discrete-ordinates theory page.
    """
    t_start = time.perf_counter()

    # Normalize the geometry declaration (legacy mesh OR axis tuple —
    # the only 3-D entry; C5.5 #225) into the SN phase space;
    # boundary_condition fills faces only when the declaration carries
    # no explicit BC.
    sn_mesh = _as_sn_mesh(
        mesh, quadrature, materials, boundary_condition, mat_map=mat_map,
        scheme=scheme,
    )

    # Issue #168 status (Phase D ERR-026 closure, 2026-05-12):
    #
    # * Phase A (Defects 1 + 2): the ``BoundaryFaceFlux`` Protocol
    #   retired in Phase C.
    # * Phase B (Defect 3): :class:`MorelMontryAngularSweep` is now
    #   the SNMesh-default :class:`PoleAngularClosureBase`.
    # * Phase C (sweep-frame matvec): the apply matvec is one sweep
    #   iteration semantically, with WDD spatial closure.
    # * Phase D (Carlson coupled-pole seed): the M-M angular
    #   recurrence's half-angle face flux is seeded so per-ordinate
    #   flat-flux balance holds on sphere Gate 1.1 MMS.  (Since #282
    #   route (a) / #280 Phase 2.5d the seed is first-class composite
    #   STATE, marched directly from the true q½ source by
    #   :func:`~orpheus.sn.spatial.psi_half_angle_seed.carlson_inward_sweep_from_source`;
    #   the original ``CarlsonInwardSweep`` proxy-source strategy is
    #   retired.)
    #
    # **ERR-058 default restoration (2026-06-12, Issue #195)**: every
    # geometry defaults to ``"source_iteration"``.  The Phase D
    # curvilinear→Krylov flip existed because the sweep's fixed point
    # was wrong on non-flat fields; the ERR-058 closure-seed fix made
    # the (unified) sweep/matvec system O(h²)-consistent, SI ≡ Krylov
    # bit-identical, and SI ~10²× faster.  See the docstring history.
    if inner_solver is None:
        inner_solver = "source_iteration"

    solver = SNSolver(
        sn_mesh,
        inner_solver=inner_solver,
        scattering_order=scattering_order,
        max_inner=max_inner, inner_tol=inner_tol,
    )

    # Normalise the external source (raw bulk array OR composite
    # TimedFullField) into the single composite RHS ``q = q_bulk ⊕ q_∂`` the
    # inner paths consume (Cardinal Rule 2 — one construction point; shape
    # validation lives inside the helper).
    q_ext_composite = _build_fixed_source_rhs(external_source, sn_mesh)

    if inner_solver == "source_iteration":
        return _solve_fixed_source_si(
            solver, sn_mesh, q_ext_composite,
            t_start, max_inner, inner_tol, inner_schedule=inner_schedule,
        )

    # Krylov path.  We solve T·ψ = b directly via GMRES, where b carries
    # the external per-ordinate source plus any in-scatter / (n,2n) terms
    # built from the converged scalar flux.  Wrapping that in an outer
    # source iteration converges scattering self-consistently.
    return _solve_fixed_source_krylov(
        solver, sn_mesh, q_ext_composite,
        t_start, max_inner, inner_tol,
    )


def _solve_fixed_source_si(
    solver: SNSolver,
    sn_mesh: SNMesh,
    q_ext_composite: "TimedFullField | CoupledField",
    t_start: float,
    max_inner: int,
    inner_tol: float,
    inner_schedule: str = "gauss_seidel",
) -> Solution:
    r"""Fixed-source path via the :class:`SourceIteration` primitive.

    Phase 1 (R1, 2026-06-04) — carved onto the SAME
    :class:`~orpheus.numerics.iteration.SourceIteration` primitive the
    eigenvalue inner :meth:`SNSolver._solve_source_iteration` uses, on the
    identical operator triple

    .. math::

        A \;=\; L + C\,, \quad
        S \;=\; \text{full multi-group scatter} + B\,, \quad
        F \;=\; 0_{\rm wg}

    differing ONLY in ``q_ext`` (the EXTERNAL source here vs the fission
    source in the eigenvalue inner) and the returned contract (a full typed
    :class:`Solution` here vs an angular-integrated scalar flux there).

    The previous hand-rolled ``for n_inner`` loop is RETIRED.  It rebuilt the
    scattering source each iterate via ``_add_scattering_source`` +
    ``_add_n2n_source`` + ``_build_aniso_scattering`` and drove ``−B`` through
    :func:`_reflect_outflow_into_inflow` — both are now subsumed by the
    primitive's coupling gains ``S.apply(ψ_n) + B.apply(ψ_n)``:
    :meth:`ScatteringOperator.apply` (the ``TimedFullField`` branch)
    recomputes the IDENTICAL ``(P0 in-scatter + (n,2n))/W + Pℓ`` bulk source,
    and the boundary gain ``B``
    (:class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`, a first-class
    coupling in the within-group system record) delivers the reflective
    ``B·ψ.outflow`` through ``rhs.boundary`` which the bare ``(L + C).solve``
    sweep reads as the inflow seed (single source of truth — Cardinal Rule 2).  The whole-trace
    :func:`_reflect_outflow_into_inflow` route is no longer needed on this
    path; it survives for the eigenvalue reconstruction sweep + Phase 3's
    octant-restricted Gauss-Seidel variant.

    Geometry-agnostic (slab / sphere / cylinder / 2-D Cartesian): the
    within-group solve carries no geometry dependence, exactly as the
    eigenvalue SI inner (Wave O "2-D SI Phase A").

    Equivalence note (vv-principles §bit-identity): the converged fixed point
    is identical to the retired loop's (same operators, same ``S`` and ``B``
    coupling gains, same WDD sweep), but the iteration TRAJECTORY differs — the primitive
    stops on the composite ``‖ψ_{n+1} − ψ_n‖ / ‖ψ_{n+1}‖`` residual (the full
    angular + boundary iterate, the same metric the verified eigenvalue inner
    uses) rather than the scalar-flux ``‖φ − φ_prev‖ / ‖φ‖``.  Converged ``φ``
    therefore agrees to ``~inner_tol`` (principled-equivalence), and
    ``history.n_inner`` / ``flux_residuals`` reflect the composite metric.
    """
    from orpheus.transport.fields.angular_flux import (
        AngularFlux,
    )
    from orpheus.transport.fields.angular_boundary_flux import (
        AngularBoundaryFlux,
    )
    from orpheus.transport.fields.scalar_flux import ScalarFlux

    # ``q_ext_composite`` is the normalised composite RHS ``q = q_bulk ⊕ q_∂``
    # built once by :func:`_build_fixed_source_rhs` (Cardinal Rule 2 — the SI
    # and Krylov paths share it).  The bulk is the per-ordinate-density
    # external source; the boundary is the prescribed inflow (zero for
    # vacuum / reflective — the reflective inflow rides ``rhs.boundary`` via
    # the ``B`` coupling gain, NOT ``q_ext``; a NON-vacuum prescribed inflow
    # is carried in ``q_ext_composite.boundary``).  Scattering (P0 + Pℓ +
    # (n,2n)) is NOT pre-staged — the primitive's ``S`` operator recomputes
    # it each iterate.

    # ── Build the within-group system + SI via the SHARED builders (single
    # source of truth — :func:`~orpheus.sn.coupled_system.build_within_group_system`
    # / :func:`_within_group_si`; identical construction to the eigenvalue
    # inner).  ``inner_schedule`` selects Jacobi vs boundary-G-S; the SI
    # builder folds in the Phase-5a angular-windowing.  ``base_resolvent``
    # (un-wrapped) + ``gains`` are kept for the final full-angular
    # reconstruction below. ────────────────────────────────────────────
    system = build_within_group_system(
        sn_mesh, solver.mat_xs, scattering_op=solver.scattering_op,
    )
    si, base_resolvent, gains, windowed = _within_group_si(
        system, sn_mesh,
        inner_schedule=inner_schedule, max_iter=max_inner, tol=inner_tol,
    )
    coupled = isinstance(system.resolvent, CoupledInvertibleOperator)

    # Cold-start iterate (x0 = zeros).  Fixed-source is a single solve — no
    # eigenvalue outer to warm-start from (cf. the eigenvalue inner's
    # ``self._psi_typed``).  Windowed → zero moments (single-sourced in
    # :func:`_windowed_cold_start`); else a zero AngularFlux — paired
    # NATIVE with a zero ψ_B on a carrying mesh (B.2d).
    q_a_ext = _system_a_member(q_ext_composite)
    # The rhs builder emits the TIMED System-A frame (its history_depth keys
    # the cold start); the parse reifies that producer contract loudly.
    if not isinstance(q_a_ext, TimedFullField):
        raise TypeError(
            f"fixed-source SI: the rhs's System-A member must be the timed "
            f"composite; got {type(q_a_ext).__name__}."
        )
    if windowed:
        initial_guess = _windowed_cold_start(
            solver.scattering_op, sn_mesh,
            history_depth=q_a_ext.history_depth,
        )
    else:
        cold = _unwindowed_cold_start(
            sn_mesh, history_depth=q_a_ext.history_depth,
        )
        initial_guess = _coupled_flux_state(cold, sn_mesh) if coupled else cold
    # ``q_ext_composite`` is already driver-ready (the coupled pair on a
    # carrying mesh — built once by :func:`_build_fixed_source_rhs`).
    psi_typed, residuals = si.solve(
        q_ext_composite, initial_guess=initial_guess,
    )
    # System A's converged member feeds the Solution contract; System B's
    # rides ``Solution.radial_characteristic`` (B.2d DP-Solution).
    psi_full = _system_a_member(psi_typed)
    # The parse reifies the driver-template contract (the solve echoes the
    # TimedFullField initial_guess member) — the static iterate type is the
    # operators' carrier, re-narrowed here loudly.
    if not isinstance(psi_full, TimedFullField):
        raise TypeError(
            f"fixed-source SI: the converged iterate must echo the timed "
            f"template; got {type(psi_full).__name__}."
        )
    converged_flag = bool(residuals) and residuals[-1] < inner_tol
    flux_residuals = [float(r) for r in residuals]

    # Issue #197 PR-TYPED-5: build typed Solution at the boundary.
    # (The former mesh / quadrature / materials parameters retired in C4 —
    # Solution never consumed them; the typed fluxes carry the SNMesh
    # reference, which transitively exposes those handles via
    # ``.mesh.{mesh, quad, materials}``.)
    history = IterationHistory(
        flux_residuals=tuple(flux_residuals),
        n_inner=len(residuals),
        total_inner_iterations=len(residuals),
        converged=converged_flag,
    )
    # ``Solution.angular_flux`` must carry the FULL per-ordinate angular flux.
    # Un-windowed: ``psi_typed`` already IS it (return directly, exactly as the
    # fixed-source Krylov path does; the boundary trace lives on
    # ``psi_typed.boundary`` — no legacy ``solver._boundary_flux`` writeback).
    # Windowed: ``psi_typed.interior`` is the moment iterate, so reconstruct the
    # full angular with ONE final full-angular solve of the converged source
    # ``q + Σ gains·ψ`` through the UN-wrapped base resolvent (mirrors the
    # eigenvalue reconstruction sweep).  Bit-identical to the un-windowed
    # converged ψ: S/B consume the moments == the full angular's moments
    # (de-risk proven), so the source is the same, and one sweep of the
    # converged source reproduces the converged iterate by the fixed point.
    if windowed:
        # Windowed ⟹ 2-D Cartesian ⟹ seedless (R12a): the iterate is the
        # fused composite (``psi_full IS psi_typed`` here) and the resolvent
        # is never the coupled bridge — the parse states the structural
        # fact loudly instead of assuming it.
        if isinstance(base_resolvent, CoupledInvertibleOperator):
            raise TypeError(
                "windowed reconstruction reached a coupled resolvent — "
                "structurally unreachable (windowing is 2-D Cartesian, "
                "seedless; the coupled arm never windows)."
            )
        rhs_final = q_a_ext
        for gain in gains:
            rhs_final = rhs_final + gain.apply(psi_full)
        angular_out = base_resolvent.inverse().apply(
            rhs_final, initial_guess=psi_full,
        )
    else:
        # Un-windowed: the (re-fused) converged iterate IS the full
        # per-ordinate angular flux.
        angular_out = psi_full
    # Scalar flux from the RETURNED full angular flux → the Solution is exactly
    # self-consistent (``scalar == ∫ angular dΩ``), matching the un-windowed
    # contract.  (For the un-windowed path ``angular_out`` IS ``psi_typed``, so
    # this is bit-identical to the prior ``psi_typed.interior.integrate_angular``.)
    # The user-facing scalar flux is the cell-AVERAGE moment (slot 0) — a
    # multi-moment closure's φ̂ slopes are internal DG structure, not the scalar
    # flux the Solution reports (#240 D5b-S3).  The parse reifies the
    # full-angular contract of BOTH arms (the reconstruction sweep emits
    # angular; the un-windowed iterate echoes the flux template) loudly.
    angular_bulk = angular_out.interior
    if not isinstance(angular_bulk, AngularFlux):
        raise TypeError(
            f"fixed-source SI: Solution.angular_flux must carry an "
            f"AngularFlux bulk; got {type(angular_bulk).__name__}."
        )
    phi = _average_moment_scalar(
        angular_bulk.integrate_angular().values, sn_mesh,
    )
    return Solution(
        angular_flux=angular_out,
        scalar_flux=ScalarFlux.from_mesh(phi, sn_mesh),
        mesh=sn_mesh,
        keff=None,
        history=history,
        radial_characteristic=_system_b_member(psi_typed),
    )


def _solve_fixed_source_krylov(
    solver: SNSolver,
    sn_mesh: SNMesh,
    q_ext_composite: "TimedFullField | CoupledField",
    t_start: float,
    max_inner: int,
    inner_tol: float,
) -> Solution:
    r"""Curvilinear-default fixed-source path: typed :class:`KrylovAcceleration`.

    R-1 Step 4 Step G1 (2026-05-22) — carved onto
    :class:`~orpheus.numerics.iteration.KrylovAcceleration` consuming
    the operator triple

    .. math::

        A \;=\; L + C\,, \quad
        S \;=\; \text{full multi-group scatter}\,, \quad
        F \;=\; 0_{\rm wg}

    on typed :class:`~orpheus.sn.angular_flux.AngularFlux`.  The
    composite ``L + C`` returns an
    :class:`~orpheus.sn.operators.streaming.InvertibleOperator`; its ``.solve``
    IS the WDD sweep, but R-1 ships GMRES UNPRECONDITIONED
    (explicit identity) — issue #200 tracks the block-inverse face
    preconditioner re-enablement.

    The outer Picard wrap on the scattering source dissolves at G1:
    typed ``KrylovAcceleration`` solves :math:`(L+C-S)\,\psi =
    q_{\rm ext}` directly via GMRES, putting scattering
    self-consistency INSIDE the Krylov polynomial.  Fission ``F`` is
    zero — this is the pure fixed-source path; the eigenvalue
    outer/within-group decomposition lives in :func:`solve_sn`.

    The dependency audit + retirement sequence eliminated the
    packed-1D vector machinery this function previously consumed:

    * ``build_equation_map_*`` — **RETIRED** in D-J (2026-05-30).
    * ``solution_to_angular_flux_*`` — **RETIRED** in D-J (2026-05-30).
    * ``pack_with_traces`` — **RETIRED** in D-J (2026-05-30).
    * ``_build_rhs_{cartesian, spherical, cylindrical}`` —
      **RETIRED** in P1.7 of the moment-space + layering plan (the
      Cartesian variant carried an inline ``(2*l+1)`` literal that
      duplicated
      :attr:`~orpheus.numerics.spaces.SphericalHarmonicSpace.addition_theorem_factor`;
      all three were already dead code in production).
    * ``_make_sweep_preconditioner`` — **RETIRED** pre-D-J.
    * ``SNStreamingOperator`` — **RETIRED** in D-K (commit
      ``dadf4e8``).  ``self.L`` now binds to the algebraic
      composition ``StreamingOperator + MultiplicationOperator`` (the
      collision multiplier ``C = M[σ_t]``)
      (= :class:`InvertibleOperator`).

    Scope
    =====

    ALL geometries — slab, sphere, cylinder, AND 2-D Cartesian (Phase 1
    "gap" un-gate, 2026-06-04).  The 2-D Cartesian fixed-source Krylov is
    the structural twin of the live 2-D eigenvalue Krylov inner
    :meth:`SNSolver._solve_krylov`: identical within-group system
    (:func:`~orpheus.sn.coupled_system.build_within_group_system`) and identical
    :class:`~orpheus.numerics.iteration.KrylovAcceleration` driver
    (:func:`_within_group_krylov`), differing ONLY in ``q_ext`` (the
    external source here vs the fission source there) — and it is the twin
    of the geometry-agnostic 2-D fixed-source SI path.  Verified: the
    converged per-ordinate flux hits the closed-form streaming equilibrium
    ``q/Σ_t`` on a homogeneous reflective box, and the SI≡Krylov flux
    shape agrees on a heterogeneous non-flat case
    (``tests/sn/solve/test_fixed_source_2d_equivalence.py``).  The legacy
    ``NotImplementedError`` 2-D guard (and the "B1'' face block" excuse it
    carried) is RETIRED — the 4-face bare-boundary architecture (Wave O
    O.4b) makes the path geometry-agnostic by construction.
    """
    from orpheus.transport.fields.angular_flux import (
        AngularFlux,
    )
    from orpheus.transport.fields.angular_boundary_flux import (
        AngularBoundaryFlux,
    )
    from orpheus.transport.fields.scalar_flux import ScalarFlux

    ng = solver.ng
    N = sn_mesh.quad.N

    # ``q_ext_composite`` is the normalised composite RHS ``q = q_bulk ⊕ q_∂``
    # built once by :func:`_build_fixed_source_rhs` (Cardinal Rule 2). B.5.2:
    # q_ext IS a source — bulk per-ordinate-density ``AngularSourceSink`` +
    # boundary ``AngularBoundarySourceSink`` prescribed inflow (zero for vacuum /
    # reflective — the reflective inflow rides ``initial_guess`` /
    # ``rhs.boundary``, not ``q_ext``; a NON-vacuum prescribed inflow IS
    # carried in ``q_ext_composite.boundary``). The Krylov matvec composes
    # operator-output sources; q_ext is raveled type-agnostically as the
    # GMRES rhs ``b``.

    # B.5.2: the FLUX initial_guess (built FIRST so the GMRES restart is sized
    # from the FULL composite ravel) fixes the Krylov solution_template (the
    # return type); x0 stays all-zeros (bit-identical to the prior cold start).
    # The template carries the φ̂ moment axis at a multi-moment closure AND the
    # ψ½ state on a carrying mesh (the ravel keys on its ``to_flat``).
    q_a_ext = _system_a_member(q_ext_composite)
    if not isinstance(q_a_ext, TimedFullField):
        raise TypeError(
            f"fixed-source Krylov: the rhs's System-A member must be the "
            f"timed composite; got {type(q_a_ext).__name__}."
        )
    krylov_cold_start = _unwindowed_cold_start(
        sn_mesh, history_depth=q_a_ext.history_depth,
    )

    # ── Build the within-group system + Krylov driver (single source of
    # truth — :func:`~orpheus.sn.coupled_system.build_within_group_system` /
    # ``_within_group_krylov``; shared with the eigenvalue Krylov + SI
    # paths; the cached scattering operator injects through the cache
    # seam). ──────────────────────────────────────────────────────────
    # ERR-053 (#282 route (a)): ``restart`` MUST cover the FULL ravel —
    # bulk ⊕ trace ⊕ ψ½ (both systems on the coupled pair — the
    # CoupledField ``to_flat`` concatenates them, the B.2a conformance
    # closure).  A bulk-sized restart re-truncates GMRES on the trace+seed
    # DOFs.  Size it from the state the driver ravels (the multi-moment φ̂
    # axis + the trace + the ψ½ state all track automatically).
    system = build_within_group_system(
        sn_mesh, solver.mat_xs, scattering_op=solver.scattering_op,
    )
    coupled = isinstance(system.resolvent, CoupledInvertibleOperator)
    if coupled:
        # The coupled pair is born native (B.2d): the flux template pairs
        # with a zero ψ_B; ``q_ext_composite`` is already the coupled rhs.
        krylov_cold_start = _coupled_flux_state(krylov_cold_start, sn_mesh)
    krylov = _within_group_krylov(
        system.resolvent, *system.gains,
        n_dof=int(krylov_cold_start.to_flat().size),
        max_iter=max_inner, tol=inner_tol,
    )

    psi_typed, residuals = krylov.solve(
        q_ext_composite, initial_guess=krylov_cold_start,
    )
    # System A's converged member feeds the Solution contract; System B's
    # rides ``Solution.radial_characteristic`` (B.2d DP-Solution).
    psi_full = _system_a_member(psi_typed)
    # D-H.1c stage 2 — the Krylov ravellable protocol unravels back to the
    # SOLUTION TEMPLATE (the flux ``initial_guess``), so the driver's static
    # iterate type (the operators' carrier) re-narrows to the timed flux
    # composite here. The parse reifies that template contract loudly
    # instead of assuming it.
    if not isinstance(psi_full, TimedFullField):
        raise TypeError(
            f"fixed-source Krylov: the converged iterate must echo the "
            f"timed flux template; got {type(psi_full).__name__}."
        )
    bulk = psi_full.interior
    if not isinstance(bulk, AngularFlux):
        raise TypeError(
            f"fixed-source Krylov: the converged iterate must carry an "
            f"AngularFlux bulk (the flux template); got {type(bulk).__name__}."
        )
    # Read bulk for scalar reduction (cell-average moment).
    phi = _average_moment_scalar(
        bulk.integrate_angular().values, sn_mesh,
    )
    converged_flag = bool(residuals) and residuals[-1] < inner_tol
    n_outer = len(residuals)
    flux_residuals = [float(r) for r in residuals]

    # Issue #197 PR-TYPED-5: build typed Solution at the boundary.
    # R-1 Step 4 G1 — ``psi_full`` carries the Krylov-converged
    # composite with the matvec's B1'' face residual on its boundary;
    # reuse directly. (The former mesh / quadrature / materials
    # parameters retired in C4 — Solution never consumed them.)
    history = IterationHistory(
        flux_residuals=tuple(flux_residuals),
        n_inner=n_outer + 1,
        total_inner_iterations=n_outer + 1,
        converged=converged_flag,
    )
    # D-H.1c stage 2 (2026-05-28): psi_full IS already a TimedFullField;
    # no adapter wrap at the Solution boundary.
    return Solution(
        angular_flux=psi_full,
        scalar_flux=ScalarFlux.from_mesh(phi, sn_mesh),
        mesh=sn_mesh,
        keff=None,
        history=history,
        radial_characteristic=_system_b_member(psi_typed),
    )
