"""SN operator algebra leaves — typed :class:`TimedFullField` contract.

Provides the four-operator algebra leaves consumed by the within-group
equation :math:`A_{\\rm wg} = L + C - S_{\\rm foldable}`:

* :class:`StreamingOperator` — :math:`L = \\Omega\\cdot\\nabla
  + \\text{angular redistribution}` (the curvilinear pole term lives
  here for sphere / cylinder).
* :class:`CollisionOperator` — :math:`C = \\sigma\\cdot` (diagonal
  in position, group, and direction).
* :class:`InvertibleOperator` — the sweep-invertible specialisation
  :math:`(L + C)` returned by ``L + C``; advertises ``CAP_SOLVE``
  via the WDD sweep.

All three operators consume and emit
:class:`~orpheus.transport.timed_full_field.TimedFullField` — the
typed composite carrier (bulk = :class:`~orpheus.transport.fields.angular_flux.AngularFlux`,
boundary = :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`).
Producer-side normalisation (Pattern 7): the typed contract is
enforced at every operator entry; no bare-ndarray packed-vector
adapter.  The geometry-agnostic kernel is
:meth:`_MSpatialOperatorSum._compute_decomposition` (1-D slab /
sphere / cylinder; dual-emission body that produces both
``M_spatial`` and ``M_angular_redist`` contributions in one
bidirectional sweep) plus :meth:`StreamingOperator._apply_2d_cartesian`
(2-D Cartesian FD).  Wave T T.5 close-out retired the module-level
``_transport_operator_matvec_unified`` helper — its body lives as
the orchestrator's private dual-emission method.

Three geometries are supported:

* **Cartesian 2D** — ``L = μ_x ∂/∂x + μ_y ∂/∂y + Σ_t``
* **Spherical 1D** — ``L = μ (A ∂/∂r)/V + (α ∂/∂μ)/V + Σ_t``
* **Cylindrical 1D** — per-level azimuthal redistribution

History
-------

* Wave E retired the standalone ``build_rhs*`` and
  ``build_transport_linear_operator*`` helpers along with the
  ``angular_flux_to_scalar`` aggregator.
* D-H (2026-05) — typed :class:`TimedFullField` composite carrier
  replaced the legacy packed-vector convention; operators bridged
  the bare-ndarray↔typed boundary internally.
* D-I (2026-05-29) — the bare-ndarray adapter retired from every
  operator leaf (L / C / S / F).  Operators consume only
  :class:`TimedFullField`.
* D-J (2026-05-30) — the supporting packed-vector codec family
  retired: :class:`EquationMap`, :func:`build_equation_map_*`,
  :func:`solution_to_angular_flux*`, :func:`pack_with_traces`.  The
  legacy 2-D Cartesian FD matvec ``transport_operator_matvec``
  (and its cartesian / spherical / cylindrical predecessors) had
  retired in D-H.2-C4e.6 (commit ``a614610``).

.. note:: Symmetric-closure invariant

   The operator :math:`L` uses the **symmetric** closure that makes
   the Krylov path agree with analytical references:

   * **Cartesian**: upwind cell-center finite differences for the
     streaming gradient — first-order accurate and consistent with DD
     on uniform meshes, first-order inconsistent on non-uniform meshes.
   * **Curvilinear**: arithmetic averages for spatial face fluxes and
     τ-weighted interpolation for angular face fluxes — distinct from
     the WDD asymmetric closure :math:`\\psi_{\\rm out} =
     (\\overline{\\psi} - (1 - \\tau)\\,\\psi_{\\rm in})/\\tau` used
     by the sweeps.

   On uniform meshes the symmetric-closure operator and the WDD
   sweep converge to the same physics in the fine-mesh limit; on
   curvilinear, the sweep's WDD closure has the ERR-026 closure-bias-
   driven self-consistent fixed point that the Krylov-on-:meth:`apply`
   path bypasses.

.. note:: Boundary-condition handling (Wave O step O.4a.2, Issue #208)

   The realized boundary law ``B`` is a **first-class sibling
   operator** of ``L``, NOT re-applied inside this matvec.  The
   canonical SN loss is ``(L_full + C - S - F - B)`` on the direct-sum
   state ``V = V_bulk ⊕ V_inflow ⊕ V_outflow``.  For the **1-D** path
   (slab / sphere / cylinder), :meth:`_MSpatialOperatorSum._compute_LpC`
   reads ``psi.boundary.inflow`` as a GIVEN, keeps the outflow
   self-consistency defect ``psi.outflow - streamed`` on the outflow
   trace row, and adds the inflow identity ``I·psi.inflow`` on the
   inflow row — with NO ``bc.apply``.  The reflective coupling
   ``psi.inflow = B·psi.outflow`` is delivered by the sibling ``-B``
   (:class:`~orpheus.sn.boundary_operator.SNBoundaryOperator`), and the
   outer Krylov / SI loop drives the boundary consistency residual
   ``psi.inflow - B·psi.outflow - q.inflow → 0``.  See
   :ref:`bc-extraction` for the full block-matrix derivation, the three
   design corrections, the two ``-B`` delivery routes, and the O.2
   forcing function.

   The **2-D Cartesian** path
   (:meth:`StreamingOperator._apply_2d_cartesian`) is ALSO bare (O.4b
   Phase E landed): it seeds the octant-incoming face slots from the
   GIVEN ``psi.boundary.inflow`` via the typed ``wavefront.seed`` (ι_*)
   with NO ``bc.apply``, walks the same per-octant
   :class:`~orpheus.sn.sweep_graph.SweepDependencyGraph` (``graph.residual``
   → the diamond-difference ``CellUpdate`` closure) the 2-D *sweep*
   ``_sweep_2d_wavefront`` uses — so matvec ≡ sweep in 2-D by
   construction (L21, one discretization) — and emits the boundary
   consistency residual (outflow defect ``streamed − given`` + inflow
   identity) as a :class:`BoundarySourceSink`.  The reflective coupling
   is the sibling ``-B`` exactly as in 1-D.  The pre-extraction Phase C
   insight that the BC must consume the WDD-propagated outflow face
   vector (not cell centres) is preserved and strengthened: post-O.4a.2
   the outflow trace is the explicit solved unknown ``psi.outflow`` that
   ``-B`` reads, closing ERR-026 by construction for the 1-D curvilinear
   path.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from functools import cached_property

from orpheus.numerics.operator import (
    BlockRole,
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    CAP_SOLVE,
    LinearOperatorMixin,
    OperatorSum,
    ZeroOperator,
)

from orpheus.numerics.quadrature import Quadrature

if TYPE_CHECKING:
    from orpheus.transport.fields.boundary_flux import BoundaryFlux
    from .geometry import SNMesh
    from orpheus.numerics.projection import MomentProjection
    from orpheus.transport.source_sinks import ScalarSourceSink, AngularSourceSink
    from .spatial.pole_angular_closure import PoleAngularClosure
    from .sweep_strategy import SweepStrategy

__all__ = [
    "StreamingOperator",
    "CollisionOperator",
    "AngularRedistributionOperator",
]
# Wave T T.5 close-out (matvec retirement, post-T.5.2):
# `_transport_operator_matvec_unified` is DELETED — its body lives
# as `_MSpatialOperatorSum._compute_decomposition`, a private method
# on the orchestrator that emits BOTH M_spatial and M_angular_redist
# contributions in ONE bidirectional sweep (dual emission).  The
# ψ-keyed cache lets `M_angular_redist.apply` reuse the spatial
# walk state at zero cost when `StreamingOperator.apply` calls
# both consumers within the same call boundary.  External consumers
# of the canonical matvec call `(L + C).apply(state)` via the
# public operator-algebra path.


# ═══════════════════════════════════════════════════════════════════════
# Wave T T.4 — `M_spatial` decomposition primitives
# ═══════════════════════════════════════════════════════════════════════
#
# `StreamingOperator.apply` exposes `M_spatial` (streaming + collision-
# share contribution per spatial axis per direction sign) and
# `M_angular_redist` (curvilinear M-M redistribution, ZeroOperator
# for slab/Cartesian).  The "L = M - C" subtractive contract lives at
# `StreamingOperator.apply`'s boundary (Pattern 7 producer-side
# normalisation; per the T.4 verification spec Q3 = (γ)).
#
# `M_spatial` is an :class:`OperatorSum` of per-direction-sign summands
# (`_SpatialSweepDirection(+1)` and `_SpatialSweepDirection(-1)`).  This
# per-direction structure is type-level — it exposes the natural
# operator algebra for: (a) Wave O typing (per-direction BC dependency
# structure: forward sweep reads inner BC and writes outer-face
# residual; backward symmetrically), (b) adjoint propagation (`M_x_fwd.H
# = M_x_back` with swapped BCs), (c) future DSA-class preconditioners
# (which traditionally split per-direction for the P_1 diffusion
# limit), (d) parallel-in-direction execution, and (e) per-direction
# debugging slicing (catches ERR-006 family signatures cleanly).
#
# Per the T.4 spec MA-Q1 master condition, the per-direction summands
# are NOT clean `TensorProductOperator` factors — the WDD sweep
# recurrence sequentially couples cells along the spatial axis, and
# the per-direction WDD outflow at the outer face is the natural
# coupling between the forward and backward sweeps via `bc_outer`.
# A clean `(D_axis & Ω_axis & I_g)` 3-factor wrap would false-assert
# separability the recurrence doesn't support (same Q2 critique
# applied to curvilinear angular redistribution).
#
# Per the T.4 verification spec follow-up (post-spec architectural
# question resolved with user concurrence): the per-direction summands
# expose structure at the type level; `_MSpatialOperatorSum.apply`
# overrides the default `OperatorSum.apply` to run ONE bidirectional
# sweep with shared local state — matches today's perf (1×) instead
# of the 1.5× cost of two independent forward sweeps that would
# otherwise be required (the backward sweep's seed depends on the
# forward sweep's WDD outflow via `bc_outer`).
#
# Slab status: `M_spatial = (L+C)` (no angular redistribution); the
# orchestrated apply delegates to `_transport_operator_matvec_unified`
# bit-exact.  Curvilinear T.4c will introduce `M_angular_redist` as
# a bespoke leaf and subtract its contribution to produce M_spatial
# alone; for now, `M_spatial.apply` IS the unified matvec output.


@dataclass(frozen=True)
class _SpatialSweepDirection(LinearOperatorMixin):
    r"""One direction sign's contribution to :attr:`StreamingOperator.M_spatial`.

    Carries the per-direction algebra exposure for the Wave T
    tensor-network refactor.  Each instance represents the action of
    the spatial-streaming + collision-share part of :math:`(L+C)` for
    ordinates of a specific direction sign (μ_x > 0 for
    ``direction_sign = +1``; μ_x < 0 for ``-1``).

    Per the T.4 verification spec follow-up (user-endorsed design):
    the production apply path goes through
    :class:`_MSpatialOperatorSum`, which orchestrates the bidirectional
    sweep in one pass with shared state.  Standalone
    :meth:`_SpatialSweepDirection.apply` is the slow fallback —
    available for testing, debugging, and the Wave-O / adjoint
    inspection paths where the per-direction algebra needs to be
    exercised in isolation.

    See Also
    --------
    _MSpatialOperatorSum : the orchestrator that owns the production apply.
    StreamingOperator.M_spatial : the property that returns the orchestrator.
    """

    sn_mesh: "SNMesh"
    sigma_t: np.ndarray
    direction_sign: int  # +1 or -1
    capabilities: frozenset[str] = field(
        default_factory=lambda: frozenset({CAP_APPLY})
    )

    def apply(self, psi: "TimedFullField") -> "TimedFullField":
        r"""Standalone per-direction action — slow path.

        Runs the FULL bidirectional matvec (because the backward
        direction's seed depends on the forward direction's WDD
        outflow via ``bc_outer``) and returns ONLY the contribution at
        ordinates of this direction sign; ordinates of the opposite
        sign carry zeros.

        In production, :meth:`_MSpatialOperatorSum.apply` overrides the
        default :meth:`OperatorSum.apply` so that the two summands
        share the bidirectional sweep cost.  Use this standalone path
        ONLY for testing, debugging, or future Wave-O / adjoint
        composition that exercises per-direction algebra in isolation.
        """
        from orpheus.transport.fields.angular_flux import AngularFlux
        from orpheus.transport.source_sinks import AngularSourceSink, BoundarySourceSink
        from orpheus.transport.timed_full_field import TimedFullField

        sn_mesh = self.sn_mesh
        # Wave T T.5 close-out (matvec retirement): construct a
        # transient orchestrator and invoke its dual-emission body
        # to recover the FULL (L+C) = M_spat + M_ang, then mask.
        # This is the slow standalone path — production
        # `StreamingOperator.apply` does NOT take this route.
        transient_orchestrator = _MSpatialOperatorSum(sn_mesh, self.sigma_t)
        m_spat, m_ang = transient_orchestrator._compute_decomposition(psi)
        full_bulk = m_spat.bulk.values + m_ang.bulk.values
        full_boundary = m_spat.boundary  # M_ang.boundary is zero (per MA-Q4)

        # Mask out ordinates of the opposite sign.
        quad = sn_mesh.quad
        mu_x = quad.mu_x
        eps = 1e-15
        if self.direction_sign > 0:
            ordinate_mask = mu_x > eps
        else:
            ordinate_mask = mu_x < -eps

        masked_bulk = full_bulk.copy()
        masked_bulk[~ordinate_mask] = 0.0

        # Boundary residual: forward sweep writes the outer-face WDD
        # outflow residual; backward sweep writes the inner-face
        # residual (for slab).  The per-direction split exposes which
        # face each direction's contribution lives on.
        masked_boundary = BoundarySourceSink.zeros_on(sn_mesh)
        full_boundary_layout = full_boundary.layout
        for face_name in full_boundary_layout.faces:
            full_face = full_boundary.face_view(face_name)
            target_face = masked_boundary.face_view(face_name)
            target_face[ordinate_mask, :] = full_face[ordinate_mask, :]

        return TimedFullField(
            bulk=AngularSourceSink.from_mesh(masked_bulk, sn_mesh),
            boundary=masked_boundary,
            _history=(),
            history_depth=psi.history_depth,
        )


class _MSpatialOperatorSum(OperatorSum):
    r""":attr:`StreamingOperator.M_spatial` — :class:`OperatorSum` of per-direction summands with orchestrated unified apply.

    Algebraically: :math:`M_{\rm spatial} = M_{x,+} + M_{x,-}` (slab,
    sphere, cylinder; 1-D spatial axis).  The two summands are
    :class:`_SpatialSweepDirection(+1)` and
    :class:`_SpatialSweepDirection(-1)`.

    Operationally: :meth:`apply` overrides :meth:`OperatorSum.apply` so
    that the bidirectional sweep runs ONCE with shared local state
    (matching the legacy `_transport_operator_matvec_unified` perf
    profile, ~1× walltime).  The default
    :meth:`OperatorSum.apply` semantics
    (``self.a.apply(x) + self.b.apply(x)``) would cost ~1.5× because
    each standalone :meth:`_SpatialSweepDirection.apply` must redo the
    forward sweep to compute the backward sweep's outer-face seed via
    ``bc_outer``.

    Per the T.4 verification spec follow-up (user-endorsed Design B):
    "orchestrated `M_spatial.apply` with shared internal state, lifting
    the local-variable WDD outflow at the outer face into a named
    coupling point between the forward and backward summands".

    Slab status: ``M_spatial`` covers ALL of ``(L+C)`` because slab has
    no angular redistribution.  Curvilinear T.4c will introduce
    :class:`AngularRedistributionOperator` (the bespoke leaf for M-M
    half-grid) and the apply will subtract that contribution to
    produce only the spatial part.
    """

    def __init__(self, sn_mesh: "SNMesh", sigma_t: np.ndarray) -> None:
        forward = _SpatialSweepDirection(sn_mesh, sigma_t, +1)
        backward = _SpatialSweepDirection(sn_mesh, sigma_t, -1)
        super().__init__(forward, backward)
        # Store mesh + xs as attributes so the orchestrated apply has
        # access without re-binding through the per-direction summands.
        self.sn_mesh = sn_mesh
        self.sigma_t = sigma_t

    def _compute_LpC(self, psi: "TimedFullField") -> "TimedFullField":
        r"""Single-emission unified ``(L+C)·ψ`` — production hot path.

        Wave T post-T.5 matvec retirement: this method inherits the
        body of the deleted ``_transport_operator_matvec_unified``
        (single output buffer, no dual emission).  Called by
        :meth:`StreamingOperator.apply` for the production hot path
        for all geometries (slab / sphere / cylinder).

        Returns the FULL ``(L+C)·ψ`` matvec result — the apply path
        subtracts ``σ_t·ψ`` at the cell boundary to recover ``L·ψ``
        per Resolution A.

        For consumers that need the algebraic split
        ``M_spatial + M_angular_redist = (L+C)``, call
        :meth:`_compute_decomposition` instead (slower; dual emission).
        """
        from orpheus.transport.fields.angular_flux import AngularFlux
        from orpheus.transport.source_sinks import AngularSourceSink, BoundarySourceSink
        from orpheus.transport.timed_full_field import TimedFullField
        from .spatial.cell_balance import cell_balance_for_streaming
        from .spatial.pole_angular_closure import MorelMontryAngularSweep

        sn_mesh = self.sn_mesh
        psi_view = psi.bulk.values
        quad = sn_mesh.quad
        N = quad.N
        ng = psi_view.shape[1]
        nx = sn_mesh.nx
        eps = 1e-15
        curvature_raw = getattr(sn_mesh, "curvature", None)
        curvature = curvature_raw if curvature_raw is not None else "cartesian"

        if curvature not in ("spherical", "cylindrical", "cartesian"):
            raise ValueError(f"Unknown curvature: {curvature!r}")
        if curvature == "cartesian" and not sn_mesh.is_1d:
            raise NotImplementedError(
                "_MSpatialOperatorSum._compute_LpC: multi-D Cartesian "
                "is not yet wired through dag_walk; only the 1-D slab "
                "is implemented.  2-D Cartesian routes through "
                "`StreamingOperator._apply_2d_cartesian` (Q1 hybrid)."
            )

        pole_angular_closure = sn_mesh.pole_angular_closure
        if pole_angular_closure is None and curvature != "cartesian":
            pole_angular_closure = MorelMontryAngularSweep()

        mu_x = quad.mu_x
        level_indices: tuple[np.ndarray, ...] = pole_angular_closure.level_indices
        A = sn_mesh.areas

        psi_g_first = psi_view.swapaxes(0, 1)
        out_g_first = np.zeros((ng, N, nx))

        V = sn_mesh.volumes
        sigma_t_gx = self.sigma_t

        boundary = psi.boundary
        trace = sn_mesh.trace
        has_inner_face = "xmin" in boundary.layout.faces
        face_outer = boundary.face_view("xmax")
        face_inner = boundary.face_view("xmin") if has_inner_face else None

        if curvature != "cartesian":
            # Curvilinear: the pole seed is the r=0 REGULARITY condition
            # (NOT a boundary condition) — read the innermost cell flux.
            pole_face_seed = psi_view[:, :, 0].copy()
        elif face_inner is not None:
            # Slab: read the GIVEN inner inflow trace (the forward sweep's
            # μ>0 seed at xmin) directly. Wave O O.4a.2 — the BC reflection
            # is NOT re-derived here; it moves to the sibling −B.
            pole_face_seed = face_inner
        else:
            raise ValueError(
                "Slab geometry requires psi.boundary.xmin_face to be "
                "populated."
            )

        # Wave O O.4a.2: the pole-closure inflow estimate is the GIVEN outer
        # inflow trace (``precompute_psi_state`` reads its μ<0 / inward
        # ordinates), not the reflected forward outflow.
        outer_inflow_estimate = face_outer
        psi_state = pole_angular_closure.precompute_psi_state(
            psi_view, sigma_t=sigma_t_gx,
            bc_outer_inflow_estimate=outer_inflow_estimate,
        )

        def _sweep_direction(
            direction_sign: int,
            psi_face_in_init: np.ndarray,
        ) -> np.ndarray:
            outflow_at_end = np.zeros((ng, N))
            for p, level_idx in enumerate(level_indices):
                level_idx_arr = np.asarray(level_idx)
                mu_level = mu_x[level_idx_arr]
                within_mask = (
                    mu_level > +eps if direction_sign > 0
                    else mu_level < -eps
                )
                if not np.any(within_mask):
                    continue
                global_dir = level_idx_arr[within_mask]
                abs_mu = np.abs(mu_x[global_dir])
                within_positions = np.where(within_mask)[0]

                cell_indices = list(
                    sn_mesh.dag_walk_cell_indices(
                        direction_sign=direction_sign, mu_level_idx=p,
                    )
                )
                if not cell_indices:
                    continue
                psi_face_in = psi_face_in_init[global_dir, :].T

                for i in cell_indices:
                    psi_cell = psi_g_first[:, global_dir, i]
                    angular_denom_term, angular_numer_upstream = (
                        pole_angular_closure.cell_contribution(
                            psi_state, i, p, within_positions,
                        )
                    )
                    A_downstream = A[i + 1] if direction_sign > 0 else A[i]
                    denom, numer_upstream = cell_balance_for_streaming(
                        abs_mu=abs_mu,
                        A_downstream=A_downstream,
                        A_total=A[i] + A[i + 1],
                        total_xs=sigma_t_gx[:, i],
                        volume=V[i],
                        psi_face_in=psi_face_in,
                        angular_denom_term=angular_denom_term,
                        angular_numer_upstream=angular_numer_upstream,
                    )
                    m_full = (denom * psi_cell - numer_upstream) / V[i]
                    out_g_first[:, global_dir, i] = m_full
                    psi_face_in = 2.0 * psi_cell - psi_face_in
                outflow_at_end[:, global_dir] = psi_face_in
            return outflow_at_end

        outflow_at_boundary = _sweep_direction(+1, pole_face_seed)
        # Wave O O.4a.2 — KEYSTONE DELETED. The backward sweep seeds from the
        # GIVEN outer inflow trace (``face_outer``'s μ<0 / inward ordinates),
        # NOT from the forward sweep's own reflected outflow
        # (``inflow_full = bc_outer.apply(outflow_at_boundary.T)``). This
        # decouples bulk ↔ boundary inside one matvec call: the reflective
        # coupling moves to the sibling −B, and the outer Krylov/SI loop drives
        # the inflow consistency ``ψ.inflow − B·ψ.outflow → 0``.
        outflow_at_inner = _sweep_direction(-1, face_outer)

        # Degenerate-ordinate branch (cylinder).
        degenerate_mask = np.abs(mu_x) < eps
        if np.any(degenerate_mask):
            global_deg = np.where(degenerate_mask)[0]
            deg_level: list[int] = []
            deg_within: list[int] = []
            for n_global in global_deg:
                for p, lvl in enumerate(level_indices):
                    lvl_arr = np.asarray(lvl)
                    pos = np.where(lvl_arr == n_global)[0]
                    if pos.size > 0:
                        deg_level.append(p)
                        deg_within.append(int(pos[0]))
                        break
            n_deg = global_deg.size
            abs_mu_deg = np.abs(mu_x[global_deg])
            zero_face = np.zeros((ng, n_deg))
            for i in range(nx):
                angular_denom_term = np.empty(n_deg)
                angular_numer_upstream = np.empty((ng, n_deg))
                for col_idx in range(n_deg):
                    denom_one, numer_one = pole_angular_closure.cell_contribution(
                        psi_state, i, deg_level[col_idx],
                        np.array([deg_within[col_idx]]),
                    )
                    angular_denom_term[col_idx] = denom_one[0]
                    angular_numer_upstream[:, col_idx] = numer_one[:, 0]

                psi_cell = psi_g_first[:, global_deg, i]
                denom, numer_upstream = cell_balance_for_streaming(
                    abs_mu=abs_mu_deg,
                    A_downstream=0.0,
                    A_total=A[i] + A[i + 1],
                    total_xs=sigma_t_gx[:, i],
                    volume=V[i],
                    psi_face_in=zero_face,
                    angular_denom_term=angular_denom_term,
                    angular_numer_upstream=angular_numer_upstream,
                )
                m_full = (denom * psi_cell - numer_upstream) / V[i]
                out_g_first[:, global_deg, i] = m_full

        m_cell = out_g_first.swapaxes(0, 1)

        # Wave O O.4a.2 — the boundary block of (L+C) carries the two trace
        # DIAGONALS of the block matrix; the off-diagonal −B is a sibling
        # operator (so this matvec contains NO BC reflection):
        #   * OUTFLOW slots — the self-consistency defect
        #     ``ψ.outflow − streamed`` (the r_outflow row's I·ψ.outflow
        #     diagonal minus L_out,b·ψ.bulk). UNCHANGED from pre-extraction;
        #     kept as ``computed − stored`` so the vacuum path is bit-identical
        #     (the per-row sign is free — q.outflow ≡ 0, the outflow trace is a
        #     pure definition with no source).
        #   * INFLOW slots — the identity ``ψ.inflow`` (the r_inflow row's
        #     I·ψ.inflow diagonal). NEW at O.4a.2. The sibling −B adds
        #     −B·ψ.outflow, so the full (L+C−S−F−B) inflow residual is
        #     ``ψ.inflow − B·ψ.outflow`` (the consistency the outer loop drives
        #     to q.inflow, the prescribed inflow / zero for vacuum+reflective).
        # The outflow / inflow ordinate sets are the disjoint sign(Ω·n)
        # partitions read from the unified TraceSpace selector (single source
        # of truth) — A.4 retired the inline ``mu_x > ±eps`` masks.
        m_boundary = BoundarySourceSink.zeros_on(sn_mesh)
        outer_outflow = trace.outflow_indices_for_face("xmax")
        if outer_outflow.size:
            m_boundary.face_view("xmax")[outer_outflow, :] = (
                outflow_at_boundary[:, outer_outflow].T
                - face_outer[outer_outflow, :]
            )
        outer_inflow = trace.inflow_indices_for_face("xmax")
        if outer_inflow.size:
            m_boundary.face_view("xmax")[outer_inflow, :] = (
                face_outer[outer_inflow, :]
            )
        if face_inner is not None:
            inner_outflow = trace.outflow_indices_for_face("xmin")
            if inner_outflow.size:
                m_boundary.face_view("xmin")[inner_outflow, :] = (
                    outflow_at_inner[:, inner_outflow].T
                    - face_inner[inner_outflow, :]
                )
            inner_inflow = trace.inflow_indices_for_face("xmin")
            if inner_inflow.size:
                m_boundary.face_view("xmin")[inner_inflow, :] = (
                    face_inner[inner_inflow, :]
                )

        return TimedFullField(
            bulk=AngularSourceSink.from_mesh(m_cell, sn_mesh),
            boundary=m_boundary,
            _history=(),
            history_depth=psi.history_depth,
        )

    def _compute_LpC_transpose(self, phi: "TimedFullField") -> "TimedFullField":
        r"""Hilbert transpose :math:`(L+C)^{\mathsf T}\,\phi` — the reverse-mode
        adjoint of :meth:`_compute_LpC`.

        Wave O / O.2b (#208).  The forward matvec is a forward-substitution
        sweep — lower-triangular in cell-visit order, with the Morel–Montry
        angular recurrence + Carlson pole seed forming a SECOND triangular
        factor in the ordinate index.  Its Euclidean transpose is the
        reverse-substitution sweep:

        * reversed cell traversal (the DD face-flux chain
          ``ψ_face_in ← 2·ψ_cell − ψ_face_in`` transposed);
        * the boundary block SWAPPED — the forward FULL operator reads the
          inflow trace and writes the outflow trace, so the transpose reads
          OUTFLOW cotangents and writes INFLOW cotangents (Lewis–Miller
          adjoint-ordinate; the same swap :meth:`SNBoundaryOperator.apply_transpose`
          already realises);
        * the angular factor reversed, delegated to
          ``closure.angular_adjoint`` (zero for slab's identity closure).

        Every coefficient is ψ-independent (geometry + σ_t): ``denom`` is
        recomputed through the SAME :func:`cell_balance_for_streaming` /
        ``cell_contribution`` the forward uses (Pattern 2 — no private-coefficient
        reach-in, no twin algebra).  Verified bit-for-bit against a dense-probe
        transpose oracle on slab / sphere / cylinder
        (``derivations/diagnostics/diag_p42_adjoint_oracle.py``).
        """
        from orpheus.transport.source_sinks import AngularSourceSink, BoundarySourceSink
        from orpheus.transport.timed_full_field import TimedFullField
        from .spatial.cell_balance import cell_balance_for_streaming

        sn_mesh = self.sn_mesh
        quad = sn_mesh.quad
        N = quad.N
        ng = phi.bulk.values.shape[1]
        nx = sn_mesh.nx
        eps = 1e-15
        curvature_raw = getattr(sn_mesh, "curvature", None)
        curvature = curvature_raw if curvature_raw is not None else "cartesian"
        if curvature == "cartesian" and not sn_mesh.is_1d:
            raise NotImplementedError(
                "_compute_LpC_transpose: the 2-D Cartesian adjoint is deferred "
                "(O.2b lands the 1-D reverse sweep first; the 2-D reverse "
                "sweep is a later Wave-O sub-step)."
            )

        closure = sn_mesh.pole_angular_closure
        mu_x = quad.mu_x
        A = sn_mesh.areas
        V = sn_mesh.volumes
        sgx = self.sigma_t                       # (ng, nx)
        trace = sn_mesh.trace
        has_inner_face = "xmin" in phi.boundary.layout.faces

        out_bar = phi.bulk.values.swapaxes(0, 1)   # (ng, N, nx)
        fo = phi.boundary.face_view("xmax")                       # (N, ng)
        fi = phi.boundary.face_view("xmin") if has_inner_face else None

        psi_bar = np.zeros((ng, N, nx))
        fo_bar = np.zeros((N, ng))
        fi_bar = np.zeros((N, ng)) if has_inner_face else None
        numer_bar = [
            np.zeros((ng, np.asarray(li).size, nx))
            for li in closure.level_indices
        ]

        # ── reverse the boundary writeback (mirror _compute_LpC m_boundary) ──
        # m.outflow = (swept outflow) − ψ.outflow;  m.inflow = ψ.inflow.
        outflow_boundary_bar = np.zeros((ng, N))    # +1 sweep outflow → xmax
        outflow_inner_bar = np.zeros((ng, N))       # −1 sweep outflow → xmin (slab); pole-discarded (curv)
        oo = trace.outflow_indices_for_face("xmax")
        oi = trace.inflow_indices_for_face("xmax")
        if oo.size:
            outflow_boundary_bar[:, oo] += fo[oo].T
            fo_bar[oo] += -fo[oo]
        if oi.size:
            fo_bar[oi] += fo[oi]
        if has_inner_face:
            io = trace.outflow_indices_for_face("xmin")
            ii = trace.inflow_indices_for_face("xmin")
            if io.size:
                outflow_inner_bar[:, io] += fi[io].T
                fi_bar[io] += -fi[io]
            if ii.size:
                fi_bar[ii] += fi[ii]

        # ── ψ-independent angular_denom_term source (dummy state) ──
        psi_state_coef = closure.precompute_psi_state(
            np.zeros((N, ng, nx)),
            sigma_t=sgx,
            bc_outer_inflow_estimate=np.zeros((N, ng)),
        )

        # ── reverse the spatial DD sweeps (both directions, per level) ──
        for p, level_idx in enumerate(closure.level_indices):
            level_idx = np.asarray(level_idx)
            mu_level = mu_x[level_idx]
            for s in (+1, -1):
                within = np.where(
                    mu_level > +eps if s > 0 else mu_level < -eps
                )[0]
                if within.size == 0:
                    continue
                gd = level_idx[within]
                abs_mu = np.abs(mu_x[gd])
                cells = list(sn_mesh.dag_walk_cell_indices(
                    direction_sign=s, mu_level_idx=p,
                ))
                if not cells:
                    continue
                f_bar = (
                    outflow_boundary_bar[:, gd] if s > 0
                    else outflow_inner_bar[:, gd]
                ).copy()
                for i in reversed(cells):
                    A_downstream = A[i + 1] if s > 0 else A[i]
                    A_total = A[i] + A[i + 1]
                    angular_denom_term, _ = closure.cell_contribution(
                        psi_state_coef, i, p, within,
                    )
                    denom, _ = cell_balance_for_streaming(
                        abs_mu=abs_mu,
                        A_downstream=A_downstream,
                        A_total=A_total,
                        total_xs=sgx[:, i],
                        volume=V[i],
                        psi_face_in=np.zeros((ng, within.size)),
                        angular_denom_term=angular_denom_term,
                        angular_numer_upstream=np.zeros((ng, within.size)),
                    )                                       # (ng, n_mask)
                    ob = out_bar[:, gd, i]
                    # reverse psi_face_in = 2·psi_cell − psi_face_in_old
                    psi_bar[:, gd, i] += 2.0 * f_bar
                    f_bar = -f_bar
                    # reverse m = (denom·ψ − |μ|A_total·psi_face_in − angular_numer)/V
                    psi_bar[:, gd, i] += denom * ob / V[i]
                    f_bar += -(abs_mu * A_total)[None, :] * ob / V[i]
                    numer_bar[p][:, within, i] += -ob / V[i]
                # reverse the sweep seed
                if s > 0:
                    if curvature != "cartesian":
                        psi_bar[:, gd, 0] += f_bar          # pole seed = ψ.bulk[:,:,0]
                    else:
                        fi_bar[gd] += f_bar.T               # slab +1 seed = ψ.inflow[xmin]
                else:
                    fo_bar[gd] += f_bar.T                   # −1 seed = ψ.inflow[xmax]

        # ── reverse the angular factor (delegated; zero for the slab closure) ──
        psi_ang_bar, bc_ang_bar = closure.angular_adjoint(
            tuple(numer_bar), sigma_t=sgx,
        )
        psi_bar += psi_ang_bar
        fo_bar += bc_ang_bar

        # ── assemble the typed composite ──
        m_boundary = BoundarySourceSink.zeros_on(sn_mesh)
        m_boundary.face_view("xmax")[...] = fo_bar
        if has_inner_face:
            m_boundary.face_view("xmin")[...] = fi_bar
        return TimedFullField(
            bulk=AngularSourceSink.from_mesh(
                psi_bar.swapaxes(0, 1), sn_mesh,
            ),
            boundary=m_boundary,
            _history=(),
            history_depth=phi.history_depth,
        )

    def _compute_decomposition(
        self, psi: "TimedFullField",
    ) -> tuple["TimedFullField", "TimedFullField"]:
        r"""Dual-emission single-pass matvec — returns ``(M_spat, M_ang)``.

        Walks the bidirectional sweep ONCE and emits both contributions
        into separate buffers:

        * ``M_spat`` carries streaming + collision (i.e. ``(L+C) -
          M_angular_redist``).  Its boundary carries the face residuals
          (only the spatial sweep writes face residuals per MA-Q4).
        * ``M_ang`` carries the curvilinear angular-redistribution
          contribution (zeros for slab/Cartesian via
          ``IdentityAngularClosure.cell_contribution``).  Its boundary
          is zero (M_angular_redist is a BulkOperator per MA-Q4).

        Wave T T.5 close-out — replaces ``_transport_operator_matvec_unified``
        as the canonical single source of truth for the 1-D matvec.
        The function-level helper retired in this commit; the body
        lives here as the operator-algebra-native orchestrator.

        Why dual emission in one walk
        ------------------------------

        The cell-balance algebra is **additive**:

        .. math::

           m_{\rm full} = m_{\rm spatial} + m_{\rm angular}

        where :math:`m_{\rm angular} = (\rm{angular\_denom\_term} \cdot
        \psi_{\rm cell} - \rm{angular\_numer\_upstream}) / V_i` and
        :math:`m_{\rm spatial} = m_{\rm full} - m_{\rm angular}` (Pattern 2
        — single source of truth via ``cell_balance_for_streaming``).

        Computing both in one cell visit costs O(1) extra arithmetic
        per cell vs the legacy single-emission path.  Without dual
        emission, the M_spatial / M_angular_redist composition costs
        ~1.7× (two walks); with dual emission + ψ-keyed cache it
        costs ~1.0× (matches the legacy shortcut).

        ψ-keyed cache
        --------------

        ``StreamingOperator.apply`` calls ``M_spatial.apply(ψ)`` then
        ``M_angular_redist.apply(ψ)`` on the SAME ψ.  The first call
        populates the cache; the second hits it and returns the
        precomputed pair without re-walking.  ``id(ψ)`` is the cache
        key — sufficient because ``TimedFullField`` is value-immutable
        within the bounds of one ``StreamingOperator.apply`` call.

        Parameters
        ----------
        psi : TimedFullField
            Angular flux + boundary trace composite.

        Returns
        -------
        (M_spatial_result, M_angular_redist_result) : tuple of TimedFullField
            Both carry the same ``history_depth`` and SNMesh as ``psi``.
            ``M_spat.bulk + M_ang.bulk == (L+C)·ψ.bulk`` bit-exact (the
            decomposition unwinds via the additive cell-balance algebra,
            modulo a per-cell FP subtraction that introduces ~ULP
            drift on the slab path).
        """
        from orpheus.transport.fields.angular_flux import AngularFlux
        from orpheus.transport.source_sinks import AngularSourceSink, BoundarySourceSink
        from orpheus.transport.timed_full_field import TimedFullField
        from .spatial.cell_balance import cell_balance_for_streaming
        from .spatial.pole_angular_closure import MorelMontryAngularSweep

        sn_mesh = self.sn_mesh
        psi_view = psi.bulk.values                                       # (N, ng, nx, ny)
        quad = sn_mesh.quad
        N = quad.N
        ng = psi_view.shape[1]
        nx = sn_mesh.nx
        eps = 1e-15
        curvature_raw = getattr(sn_mesh, "curvature", None)
        curvature = curvature_raw if curvature_raw is not None else "cartesian"

        if curvature not in ("spherical", "cylindrical", "cartesian"):
            raise ValueError(f"Unknown curvature: {curvature!r}")
        if curvature == "cartesian" and not sn_mesh.is_1d:
            raise NotImplementedError(
                "_MSpatialOperatorSum._compute_decomposition: multi-D Cartesian "
                "is not yet wired through dag_walk; only the 1-D slab "
                "is implemented.  2-D Cartesian routes through "
                "`StreamingOperator._apply_2d_cartesian` (Q1 hybrid)."
            )

        pole_angular_closure = sn_mesh.pole_angular_closure
        if pole_angular_closure is None and curvature != "cartesian":
            pole_angular_closure = MorelMontryAngularSweep()

        mu_x = quad.mu_x
        level_indices: tuple[np.ndarray, ...] = pole_angular_closure.level_indices
        A = sn_mesh.areas                                                # (nx+1,)

        psi_g_first = psi_view.swapaxes(0, 1)                     # (ng, N, nx, ny)
        out_spat_g_first = np.zeros((ng, N, nx))
        out_ang_g_first = np.zeros((ng, N, nx))

        V = sn_mesh.volumes                                        # (nx,)
        sigma_t_gx = self.sigma_t                               # (ng, nx)

        boundary = psi.boundary
        trace = sn_mesh.trace
        has_inner_face = "xmin" in boundary.layout.faces
        face_outer = boundary.face_view("xmax")                          # (N, ng)
        face_inner = (
            boundary.face_view("xmin") if has_inner_face else None
        )

        if curvature != "cartesian":
            # Curvilinear: pole seed = r=0 REGULARITY condition (not a BC).
            pole_face_seed = psi_view[:, :, 0].copy()                 # (N, ng)
        elif face_inner is not None:
            # Slab: read the GIVEN inner inflow trace (μ>0 seed at xmin)
            # directly. Mirrors _compute_LpC; the BC reflection moves to −B.
            pole_face_seed = face_inner                                  # (N, ng)
        else:
            raise ValueError(
                "Slab geometry requires psi.boundary.xmin_face to be "
                "populated (R-1 Step 4 Step G0 retired the cell-centre "
                "proxy fallback inside the matvec; legacy SN consumers "
                "must build a BoundaryFlux carrying the cell-centre proxy "
                "as face state at their call site)."
            )

        # Wave O O.4a.2: the pole-closure inflow estimate is the GIVEN outer
        # inflow trace (precompute reads its inward ordinates), mirroring
        # _compute_LpC.
        outer_inflow_estimate = face_outer                               # (N, ng)
        psi_state = pole_angular_closure.precompute_psi_state(
            psi_view, sigma_t=sigma_t_gx,
            bc_outer_inflow_estimate=outer_inflow_estimate,
        )

        def _sweep_direction(
            direction_sign: int,
            psi_face_in_init: np.ndarray,                                # (N, ng)
        ) -> np.ndarray:                                                 # (ng, N)
            outflow_at_end = np.zeros((ng, N))
            for p, level_idx in enumerate(level_indices):
                level_idx_arr = np.asarray(level_idx)
                mu_level = mu_x[level_idx_arr]
                within_mask = (
                    mu_level > +eps if direction_sign > 0
                    else mu_level < -eps
                )
                if not np.any(within_mask):
                    continue
                global_dir = level_idx_arr[within_mask]
                abs_mu = np.abs(mu_x[global_dir])
                within_positions = np.where(within_mask)[0]

                cell_indices = list(
                    sn_mesh.dag_walk_cell_indices(
                        direction_sign=direction_sign, mu_level_idx=p,
                    )
                )
                if not cell_indices:
                    continue
                psi_face_in = psi_face_in_init[global_dir, :].T          # (ng, n_dir_p)

                for i in cell_indices:
                    psi_cell = psi_g_first[:, global_dir, i]
                    angular_denom_term, angular_numer_upstream = (
                        pole_angular_closure.cell_contribution(
                            psi_state, i, p, within_positions,
                        )
                    )
                    A_downstream = A[i + 1] if direction_sign > 0 else A[i]
                    denom, numer_upstream = cell_balance_for_streaming(
                        abs_mu=abs_mu,
                        A_downstream=A_downstream,
                        A_total=A[i] + A[i + 1],
                        total_xs=sigma_t_gx[:, i],
                        volume=V[i],
                        psi_face_in=psi_face_in,
                        angular_denom_term=angular_denom_term,
                        angular_numer_upstream=angular_numer_upstream,
                    )
                    m_full = (denom * psi_cell - numer_upstream) / V[i]
                    # Dual emission: per-cell angular + spatial contributions.
                    m_ang = (
                        angular_denom_term[None, :] * psi_cell
                        - angular_numer_upstream
                    ) / V[i]
                    m_spat = m_full - m_ang
                    out_spat_g_first[:, global_dir, i] = m_spat
                    out_ang_g_first[:, global_dir, i] = m_ang

                    psi_face_in = 2.0 * psi_cell - psi_face_in
                outflow_at_end[:, global_dir] = psi_face_in
            return outflow_at_end

        outflow_at_boundary = _sweep_direction(
            direction_sign=+1, psi_face_in_init=pole_face_seed,
        )

        # Wave O O.4a.2 — KEYSTONE DELETED (twin of _compute_LpC): the backward
        # sweep seeds from the GIVEN outer inflow trace (face_outer's inward
        # ordinates), NOT from the forward sweep's own reflected outflow. The
        # reflective coupling moves to the sibling −B.
        outflow_at_inner = _sweep_direction(
            direction_sign=-1, psi_face_in_init=face_outer,
        )

        # Degenerate-ordinate branch (cylinder with degenerate ordinates).
        degenerate_mask = np.abs(mu_x) < eps
        if np.any(degenerate_mask):
            global_deg = np.where(degenerate_mask)[0]
            deg_level: list[int] = []
            deg_within: list[int] = []
            for n_global in global_deg:
                for p, lvl in enumerate(level_indices):
                    lvl_arr = np.asarray(lvl)
                    pos = np.where(lvl_arr == n_global)[0]
                    if pos.size > 0:
                        deg_level.append(p)
                        deg_within.append(int(pos[0]))
                        break
            n_deg = global_deg.size
            abs_mu_deg = np.abs(mu_x[global_deg])
            zero_face = np.zeros((ng, n_deg))
            for i in range(nx):
                angular_denom_term = np.empty(n_deg)
                angular_numer_upstream = np.empty((ng, n_deg))
                for col_idx in range(n_deg):
                    denom_one, numer_one = pole_angular_closure.cell_contribution(
                        psi_state, i, deg_level[col_idx],
                        np.array([deg_within[col_idx]]),
                    )
                    angular_denom_term[col_idx] = denom_one[0]
                    angular_numer_upstream[:, col_idx] = numer_one[:, 0]

                psi_cell = psi_g_first[:, global_deg, i]              # (ng, n_deg)
                denom, numer_upstream = cell_balance_for_streaming(
                    abs_mu=abs_mu_deg,
                    A_downstream=0.0,
                    A_total=A[i] + A[i + 1],
                    total_xs=sigma_t_gx[:, i],
                    volume=V[i],
                    psi_face_in=zero_face,
                    angular_denom_term=angular_denom_term,
                    angular_numer_upstream=angular_numer_upstream,
                )
                m_full = (denom * psi_cell - numer_upstream) / V[i]
                m_ang = (
                    angular_denom_term[None, :] * psi_cell
                    - angular_numer_upstream
                ) / V[i]
                m_spat = m_full - m_ang
                out_spat_g_first[:, global_deg, i] = m_spat
                out_ang_g_first[:, global_deg, i] = m_ang

        m_spat_cell = out_spat_g_first.swapaxes(0, 1)             # (N, ng, nx, ny)
        m_ang_cell = out_ang_g_first.swapaxes(0, 1)               # (N, ng, nx, ny)

        # M_spatial carries the face residuals (only the spatial sweep
        # writes them; per MA-Q4 M_angular_redist is a BulkOperator).
        # Outflow set read from the unified TraceSpace selector (single
        # source of truth for sign(Ω·n) — see A.4 in _compute_LpC).
        # Wave O O.4a.2 (mirror of _compute_LpC): the (L+C) boundary block
        # carries the two trace diagonals — outflow slots = self-consistency
        # defect ψ.outflow − streamed (kept; vacuum bit-identical), inflow
        # slots = identity ψ.inflow (the r_inflow diagonal; the sibling −B adds
        # −B·ψ.outflow). M_angular_redist stays zero-boundary (a BulkOperator).
        m_spat_boundary = BoundarySourceSink.zeros_on(sn_mesh)
        outer_outflow = trace.outflow_indices_for_face("xmax")
        if outer_outflow.size:
            m_spat_boundary.face_view("xmax")[outer_outflow, :] = (
                outflow_at_boundary[:, outer_outflow].T
                - face_outer[outer_outflow, :]
            )
        outer_inflow = trace.inflow_indices_for_face("xmax")
        if outer_inflow.size:
            m_spat_boundary.face_view("xmax")[outer_inflow, :] = (
                face_outer[outer_inflow, :]
            )
        if face_inner is not None:
            inner_outflow = trace.outflow_indices_for_face("xmin")
            if inner_outflow.size:
                m_spat_boundary.face_view("xmin")[inner_outflow, :] = (
                    outflow_at_inner[:, inner_outflow].T
                    - face_inner[inner_outflow, :]
                )
            inner_inflow = trace.inflow_indices_for_face("xmin")
            if inner_inflow.size:
                m_spat_boundary.face_view("xmin")[inner_inflow, :] = (
                    face_inner[inner_inflow, :]
                )

        m_spat_tff = TimedFullField(
            bulk=AngularSourceSink.from_mesh(m_spat_cell, sn_mesh),
            boundary=m_spat_boundary,
            _history=(),
            history_depth=psi.history_depth,
        )
        m_ang_tff = TimedFullField(
            bulk=AngularSourceSink.from_mesh(m_ang_cell, sn_mesh),
            boundary=BoundarySourceSink.zeros_on(sn_mesh),
            _history=(),
            history_depth=psi.history_depth,
        )

        return m_spat_tff, m_ang_tff

    def materialize_inverse_cache(self):
        r"""Build the :class:`~orpheus.sn.spatial.sweep_cache.CollisionCache`
        for this M_spatial / σ_t pair.

        Wave T T.5 close-out — exposes the cache from M_spatial's
        natural angle (Pattern 2 dual-view of the
        :meth:`CollisionCache.from_geometry` primitive).  The cache
        holds the per-cell-per-group-per-ordinate DD scan coefficients
        (`inverse_denom`, `a_attenuation`, `cumprod_a`) that the sweep
        path consumes; computing them from M_spatial documents the
        relationship "this operator IS the matvec whose inverse-on-
        the-sweep-path uses these cached coefficients".

        Returns
        -------
        CollisionCache
            A fresh cache instance — increments
            :attr:`CollisionCache._build_count`.  Per the cache
            invariance contract (sweep_cache.py:356-361), within a
            constant-σ_t epoch the cache MUST be built exactly once;
            callers MUST not invoke this method redundantly inside an
            inner Krylov / SI iteration.

        Notes
        -----
        Future leverage opportunity (post-T.5 cache-unification
        micro-wave): :func:`~orpheus.sn.sweep._ensure_coll_cache` and
        :class:`~orpheus.sn.solver.SNSolver` would route through this
        method as the canonical cache-construction path, making
        M_spatial the single source of truth for its own inverse
        cache.  Today both call :meth:`CollisionCache.from_geometry`
        directly; the unification requires SNSolver to thread the
        StreamingOperator (or M_spatial) into the cache build, which
        is a separate refactor.
        """
        from orpheus.sn.spatial.sweep_cache import (
            CollisionCache,
            GeometryCoefficients,
        )

        geom = GeometryCoefficients.from_mesh_and_quad(self.sn_mesh)
        # 1-D drop of y axis for the sweep cache contract (it stores
        # per-(N, ng, nx); 2-D Cartesian routes through
        # `_apply_2d_cartesian` and does not use this cache).
        sig_t_1d = self.sigma_t
        return CollisionCache.from_geometry(geom, sig_t_1d)

    def apply(self, psi: "TimedFullField") -> "TimedFullField":
        r"""Orchestrated apply — returns the spatial part of the decomposition.

        Wave T T.5 close-out (matvec retirement): delegates to
        :meth:`_compute_decomposition` which walks the bidirectional
        sweep ONCE and emits both M_spatial and M_angular_redist
        contributions.  The ψ-keyed cache lets
        ``StreamingOperator.apply``'s subsequent ``M_angular_redist.apply(ψ)``
        call reuse the walk state at zero cost.

        For slab: ``M_spatial.apply(ψ) == (L+C)·ψ`` bit-exact because
        the slab cell-balance has zero angular contribution
        (``IdentityAngularClosure.cell_contribution`` returns zeros).

        For curvilinear: ``M_spatial.apply(ψ) == (L+C)·ψ −
        M_angular_redist.apply(ψ)`` per the additive cell-balance
        algebra.  The subtraction is realised in the dual-emission
        body (single walk) — no second-pass cost.
        """
        m_spat, _ = self._compute_decomposition(psi)
        return m_spat


@dataclass(frozen=True)
class AngularRedistributionOperator(LinearOperatorMixin):
    r"""Bespoke curvilinear angular-redistribution leaf — :attr:`StreamingOperator.M_angular_redist` for sphere / cylinder.

    Per the Wave T T.4 verification spec Q2 (user-endorsed) — the M-M
    half-grid recurrence (Hébert §3.9.4 Eqs. 3.432-3.435) is
    sequentially coupled along the angular axis: ``α_{m+1/2}`` recurs
    from ``α_{m-1/2}`` with absorption coefficients that depend on σ_t
    and the upstream half-angle ψ.  This is a SWEEP operator, NOT a
    diagonal per-angular-axis factor.  A 3-factor TP wrap
    ``(leaf & I_x & I_g)`` would false-assert separability the
    recurrence doesn't support — same MA-Q1 master condition that
    applies to the spatial WDD recurrence in :class:`_SpatialSweepDirection`:
    coupled physics falls back to ``OperatorSum`` summands, NOT clean
    tensor products.

    Per-cell algebra (extracted from `cell_balance_for_streaming`)
    --------------------------------------------------------------

    The cell-balance algebra decomposes additively:

    .. math::

       (L+C)\,\psi\big|_{\rm cell} \;=\; m_{\rm streaming} \;+\;
                                         m_{\rm collision} \;+\;
                                         m_{\rm angular\_redist}

    where

    .. math::

       m_{\rm angular\_redist} \;=\;
         \frac{1}{V_i}\bigl[
            \mathrm{angular\_denom\_term} \cdot \psi_{\rm cell}
          - \mathrm{angular\_numer\_upstream}
         \bigr]

    The per-level M-M closure supplies
    ``angular_denom_term = (ΔA/w) · c_out`` and
    ``angular_numer_upstream = (ΔA/w) · c_in · ψ_{m-1/2, i, g}`` via
    :meth:`PoleAngularClosure.cell_contribution`.  Calling
    ``cell_contribution`` with ``within_positions = arange(n_p)``
    (all ordinates in the level) yields the full per-level
    contribution in one shot — the matvec body's per-direction split
    is an optimization, not an algebraic requirement.

    Boundary residual semantics (per the spec MA-Q4)
    -------------------------------------------------

    M_angular_redist writes ONLY to bulk; the angular redistribution
    is an interior-cell operation that doesn't traverse the spatial
    boundary.  Output ``boundary`` is the zero
    :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`.

    Slab status
    ------------

    For Cartesian (slab + 2-D), the property dispatch in
    :class:`StreamingOperator` returns :class:`ZeroOperator` instead
    of this leaf — there is no curvilinear redistribution to compute
    (``IdentityAngularClosure.cell_contribution`` returns zeros by
    construction).  This leaf is instantiated ONLY for sphere /
    cylinder.

    Parameters
    ----------
    sn_mesh : SNMesh
        Curvilinear (spherical or cylindrical) SN mesh carrying the
        M-M ``pole_angular_closure`` strategy.
    sigma_t : np.ndarray
        Per-group per-cell total cross section, ``(ng, nx, ny)``.
        Bound at constructor time per Resolution A's subtractive
        contract (the M-M Carlson coupled-pole seed is rational in
        σ_t — Hébert §3.9.4 Eqs. 3.432-3.435).
    """

    sn_mesh: "SNMesh"
    sigma_t: np.ndarray
    m_spatial: "_MSpatialOperatorSum"
    capabilities: frozenset[str] = field(
        default_factory=lambda: frozenset({CAP_APPLY})
    )

    def apply(self, psi: "TimedFullField") -> "TimedFullField":
        r"""Returns the angular-redistribution part of the M_spatial decomposition.

        Wave T T.5 close-out (matvec retirement): delegates to the
        parent orchestrator's :meth:`_MSpatialOperatorSum._compute_decomposition`,
        which walks the bidirectional sweep ONCE and emits both
        ``(M_spatial, M_angular_redist)`` contributions.  Standalone
        usage walks the full bidirectional sweep — slow path; the
        production hot path in ``StreamingOperator.apply`` calls
        ``_compute_decomposition`` directly and reads BOTH outputs in
        one shot to avoid the redundant second walk.

        Returns a :class:`TimedFullField` with bulk = the
        redistribution contribution and boundary = zeros (per MA-Q4
        — only the spatial sweep writes face residuals).
        """
        sn_mesh = self.sn_mesh
        if sn_mesh is not psi.bulk.mesh:
            raise ValueError(
                "AngularRedistributionOperator.apply: operator and "
                "composite must share the SAME SNMesh instance "
                "(mesh-identity invariant)."
            )
        if getattr(sn_mesh, "pole_angular_closure", None) is None:
            raise ValueError(
                "AngularRedistributionOperator requires a curvilinear "
                "SNMesh with a bound `pole_angular_closure`.  "
                "Cartesian meshes should route through "
                "`StreamingOperator.M_angular_redist`, which returns "
                "`ZeroOperator` and never instantiates this leaf."
            )
        _, m_ang = self.m_spatial._compute_decomposition(psi)
        return m_ang


@dataclass
class StreamingOperator(LinearOperatorMixin):
    r"""Pure streaming + angular-redistribution operator :math:`L` as a
    :class:`~orpheus.numerics.operator.LinearOperator` leaf.

    The "L" of the Phase G four-operator algebra
    :math:`A_{\rm wg} = L + C - S_{\rm foldable}`. Carries the spatial
    streaming math plus the curvilinear angular redistribution — the
    cell-collision term lives in the separate :class:`CollisionOperator`
    leaf. The split lets the within-group operator be pure algebra
    (``L + C - S.foldable_part()``); no ``WithinGroupOperator`` wrapper.

    Resolution A — subtractive ``apply``
    ------------------------------------

    .. math::

        L.{\rm apply}(\psi) \;:=\; M(\psi;\;\sigma_t) \;-\;
                                  \sigma_t \odot \psi_{\rm packed}

    The matvec primitive is called at the user's σ_t (constructor-stored),
    then the cell-collision term :math:`\sigma_t \odot \psi` is
    subtracted at the packed-vector level. Combined with
    :math:`C.{\rm apply}(\psi) = \sigma_t \odot \psi`, this gives
    :math:`(L + C).{\rm apply}(\psi) = M(\psi;\;\sigma_t)` bit-exact.

    Why σ_t is a CONSTRUCTOR parameter (Pattern 4)
    ----------------------------------------------

    The discrete curvilinear matvec is RATIONAL in σ_t through the
    Carlson coupled-pole seed (Hébert §3.9.4 Eqs. 3.432-3.435):

    .. math::

        \bar\phi_i \;=\; \frac{\Delta r_i\,\bar Q_i + 2\,\bar\phi_{i+1/2}}
                              {\Delta r_i\,\Sigma_t(i) + 2}

    where :math:`\bar Q_i = \Sigma_t(i)\,\phi_0(i) / W` at ℓ=0. The
    discrete L (with M-M angular closure) is therefore intrinsically
    σ_t-coupled — analogous to the DD coefficient :math:`\alpha_{DD}
    (\sigma_t\,\Delta x)` in characteristic-line methods, or the
    exponential characteristic :math:`\exp(-\sigma_t\,s)` in MoC/CP.
    Discrete streaming operators routinely carry σ_t through their
    closure choices.

    The CONTINUOUS L (:math:`\Omega\cdot\nabla\psi + (1-\mu^2)/r\,
    \partial\psi/\partial\mu`) is σ_t-independent. The DISCRETE L's
    σ_t parametrisation is a property of the discretisation. Pattern 4
    (illegal states unrepresentable) — ``StreamingOperator(sn_mesh)``
    without σ_t is not a meaningful object for curvilinear SN; the
    constructor refusing to construct without σ_t encodes that fact.

    Capability set
    --------------

    ``frozenset({CAP_APPLY})`` — pure streaming alone is **not
    invertible** (the streaming operator is rank-deficient without a
    collision term to make the within-group cell balance non-singular).
    The ``solve`` capability appears at the
    :class:`~orpheus.numerics.operator.OperatorSum` level: ``(L + C
    - S_foldable).solve(q)`` would route to the within-group sweep via a
    σ_r fusion hook — but ⚠ that σ_r-sweep is exact ONLY for isotropic flux
    (it removes the diagonal-in-angle ``Σ_s0·I``, not the isotropic-projection
    ``Σ_s0·P_iso``); wiring it as a within-group accelerator ships 46–56 %
    errors on anisotropic problems (issue #215; the stable+correct fold is
    consistent DSA #2 / Krylov). ``apply_transpose`` IS available
    (Wave O / O.2b) — the analytic reverse-direction adjoint matvec
    :math:`L^{\mathsf T}` (see :meth:`apply_transpose`), so the operator
    advertises ``CAP_APPLY_TRANSPOSE`` and ``L.H`` is the physical G-adjoint.

    Parameters
    ----------
    sn_mesh : SNMesh
        The augmented geometry carrying quadrature, BCs, pole closure,
        and (for curvilinear) the precomputed connection coefficients.
        ``StreamingOperator`` reads ``sn_mesh.bc_*`` directly — no
        ``boundary`` constructor parameter.
    sigma_t : np.ndarray
        Total cross-section, shape ``(ng, nx, ny)`` (Issue #196 PR-INDEX-3).
        Carried at constructor time per Resolution A's subtractive
        definition.

    Notes
    -----
    Depth B D-I.3d (2026-05-29) — :meth:`apply` accepts ONLY
    :class:`~orpheus.transport.timed_full_field.TimedFullField`.
    Depth B D-J (2026-05-30) — the underlying packed-vector codec
    family (``EquationMap`` / ``build_equation_map_*`` /
    ``solution_to_angular_flux_*`` / ``pack_with_traces``) retired
    too.  Producer-side normalisation (Pattern 7) at the operator:
    the entire matvec consumes / emits the typed direct-sum carrier,
    never any legacy packed flat layout.
    """

    sn_mesh: "SNMesh"
    sigma_t: np.ndarray

    capabilities: frozenset[str] = field(
        default_factory=lambda: frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})
    )
    # Streaming is the sole FULL operator — it couples bulk ↔ boundary
    # (reads the inflow trace to seed the sweep, writes the outflow
    # trace). Issue #208 / Wave O. Class-level constant (unannotated so
    # the dataclass does not treat it as a field).
    # CAP_APPLY_TRANSPOSE (Wave O / O.2b): the analytic reverse-direction
    # adjoint matvec landed — see :meth:`apply_transpose`.
    block_role = BlockRole.FULL

    # D-J (2026-05-30): ``_eq_map`` / ``_ensure_eq_map`` / ``n_unknowns``
    # retired alongside the :class:`EquationMap` codec family — the
    # typed :class:`~orpheus.transport.timed_full_field.TimedFullField`
    # contract has no need for the legacy packed-vector slot map.

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        r"""The composite carrier :math:`V_{\rm bulk}\oplus V_{\rm trace}` (Wave O / O.2b).

        :math:`L` is the sole FULL operator — it couples bulk :math:`\leftrightarrow`
        boundary (seeds the sweep from the inflow trace, emits the outflow
        trace). Advertising :attr:`~orpheus.sn.geometry.SNMesh.full_field_space`
        is what lets :class:`~orpheus.numerics.operator._AdjointOperator`
        read the **block-diagonal G-adjoint metric** (bulk :math:`V\,w_n`
        :math:`\oplus` trace :math:`|\Omega\cdot\hat n|\,w_n`) for ``L.H`` —
        without it the adjoint silently reduces to the metric-blind Euclidean
        transpose (Issue #208 risk R5).

        ``C`` / ``S`` / ``F`` report ``None`` (skipped by the composition
        guard), so ``L``'s composite domain propagates through ``OperatorSum``
        to the **transpose-closed** sub-sums whose ``.H`` is actually reachable
        — ``(L + C)`` and ``(L + C - B)`` — and every bulk leaf in those is
        G-conjugated for free by the op-level :math:`G^{-1}(\sum \text{leaf}^{\mathsf T})G`.
        The full prompt loss ``(L + C - S - F - B).H`` is **intentionally
        unreachable**: ``S`` / ``F`` advertise no ``apply_transpose``, so
        ``OperatorSum`` does not propagate ``CAP_APPLY_TRANSPOSE`` and
        :class:`~orpheus.numerics.operator._AdjointOperator` raises
        :class:`MissingCapability` (fails loud, never silently Euclidean —
        the capability lattice makes the metric-blind state unrepresentable).
        The foldable scattering / fission contributions are handled at the
        eigenvalue / DSA outer layer, not via this within-group adjoint.
        """
        return self.sn_mesh.full_field_space

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        # Endomorphism on the composite (see :meth:`domain`).
        return self.sn_mesh.full_field_space

    def apply(self, psi: "TimedFullField") -> "TimedFullField":
        r"""Subtractive forward action :math:`L\,\psi = M(\psi;\sigma_t)
        - \sigma_t \odot \psi.bulk`.

        Geometry dispatch is the polymorphic :attr:`sweep_strategy` (the
        matvec twin of the forward sweep): 1-D slab/sphere/cylinder →
        :meth:`_apply_1d` (``CumprodScan``); 2-D Cartesian →
        :meth:`_apply_2d_cartesian` (``MovingFrontierWindow``).

        Resolution A — :math:`(L + C).{\rm apply}(\psi) \equiv
        M(\psi;\sigma_t)` bit-exact: the per-cell σ_t·ψ cancellation
        lives at the algebra layer; this leaf returns the subtractive
        :math:`L`-only action on the typed direct-sum carrier.

        D-I.3d (2026-05-29) collapsed the bare-ndarray packed-vector
        adapter; D-J (2026-05-30) retired the supporting codec family.
        Producer-side normalisation (Pattern 7) — the operator
        consumes / emits ONLY
        :class:`~orpheus.transport.timed_full_field.TimedFullField`;
        the typed↔packed boundary collapses at the producer site.

        ``L`` is the ONLY operator that emits a non-zero face
        residual on its output ``.boundary`` —
        :class:`CollisionOperator`,
        :class:`~orpheus.sn.scattering.ScatteringOperator`, and
        :class:`~orpheus.sn.fission.FissionOperator` all leave the
        output boundary at the auto-allocated zero.

        Parameters
        ----------
        psi : TimedFullField
            Composite carrier with bulk
            (:class:`~orpheus.transport.fields.angular_flux.AngularFlux`)
            and boundary
            (:class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`).
            Operator and ``psi.bulk.mesh`` MUST be the same
            :class:`~orpheus.sn.geometry.SNMesh` instance.

        Returns
        -------
        TimedFullField
            ``L·ψ`` as a typed composite — bulk carries the cell
            action, boundary carries the face residual at the
            layout-assigned trace slots (non-zero at outer face for
            curvilinear; non-zero at outer + inner faces for slab).
        """
        from orpheus.transport.timed_full_field import TimedFullField

        if not isinstance(psi, TimedFullField):
            raise TypeError(
                "StreamingOperator.apply: expected TimedFullField, got "
                f"{type(psi).__name__}.  D-I.3d (2026-05-29) retired the "
                "bare-ndarray packed-vector contract; construct a typed "
                "composite via ``TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh)`` or "
                "explicit ``TimedFullField(bulk=..., boundary=...)``."
            )

        if self.sn_mesh is not psi.bulk.mesh:
            raise ValueError(
                "StreamingOperator.apply: operator and composite must "
                "share the SAME SNMesh instance (mesh-identity invariant)."
            )

        # Geometry dispatch is the polymorphic :attr:`sweep_strategy` — the
        # matvec twin of the forward sweep (L21: sweep and matvec are
        # different applications of the same operator). 1-D → CumprodScan
        # (the geometry-blind :meth:`_apply_1d`); 2-D Cartesian →
        # MovingFrontierWindow (:meth:`_apply_2d_cartesian`). The scattered
        # ``curv is None and not is_1d`` branch is retired into
        # ``strategy.residual``.
        return self.sweep_strategy.residual(self, psi)

    def _apply_1d(self, psi: "TimedFullField") -> "TimedFullField":
        r"""1-D ``L·ψ`` (slab + sphere + cylinder) — the CumprodScan matvec twin.

        The single-emission ``(L+C)·ψ`` via
        :meth:`_MSpatialOperatorSum._compute_LpC` (perf-optimal, no dual
        emission), then ``− σ_t·ψ`` to recover ``L`` per Resolution A.
        Extracted verbatim from :meth:`apply`'s 1-D branch by the
        SweepStrategy carve (S2) so the geometry dispatch is the polymorphic
        ``strategy.residual`` — bit-identical.

        The ``M_spatial`` / ``M_angular_redist`` algebra decomposition is
        exposed for introspection via :meth:`_compute_decomposition` (slower,
        dual-emission); the hot path uses the single-emission ``_compute_LpC``
        because re-splitting and re-summing is wasted work when ``L·ψ`` is all
        that is needed.
        """
        from orpheus.transport.source_sinks import AngularSourceSink
        from orpheus.transport.timed_full_field import TimedFullField

        LpC_result = self.M_spatial._compute_LpC(psi)
        cell_values = (
            LpC_result.bulk.values
            - self.sigma_t[None] * psi.bulk.values
        )
        return TimedFullField(
            bulk=AngularSourceSink.from_mesh(cell_values, self.sn_mesh),
            boundary=LpC_result.boundary,
            _history=(),
            history_depth=psi.history_depth,
        )

    def apply_transpose(self, phi: "TimedFullField") -> "TimedFullField":
        r"""Hilbert transpose :math:`L^{\mathsf T}\,\phi` (Wave O / O.2b, #208).

        Resolution A: :math:`L = (L + C) - C` with :math:`C = \sigma_t\odot`
        a self-adjoint diagonal, so :math:`L^{\mathsf T} = (L+C)^{\mathsf T} - C`.
        The :math:`(L+C)^{\mathsf T}` core is
        :meth:`_MSpatialOperatorSum._compute_LpC_transpose` (the reverse-mode
        adjoint of the forward sweep); the :math:`-\sigma_t\odot\phi.bulk`
        subtraction mirrors :meth:`apply`'s Resolution-A factoring (``C`` has
        no boundary action, so the transpose's boundary block is
        :math:`(L+C)^{\mathsf T}`'s).

        This returns the **plain Euclidean transpose** :math:`L^{\mathsf T}`.
        The metric conjugation :math:`G^{-1}\!\cdot^{\mathsf T}\!\cdot G` of the
        physical **G-adjoint** ``L.H`` is applied AROUND this by
        :class:`~orpheus.numerics.operator._AdjointOperator`, which reads the
        ``domain`` / ``codomain`` ``inner_product_weights`` (bulk volume on the
        cell block, the ``|Ω·n|·w`` partial-current metric on the trace block).

        Verified against a dense-probe transpose oracle (slab / sphere /
        cylinder) — ``derivations/diagnostics/diag_p42_adjoint_oracle.py``;
        the reciprocity gate ``test_phase_c_gates`` Gate 1.3 pins it in CI.
        """
        from orpheus.transport.timed_full_field import TimedFullField

        if not isinstance(phi, TimedFullField):
            raise TypeError(
                "StreamingOperator.apply_transpose: expected TimedFullField, "
                f"got {type(phi).__name__}."
            )
        if self.sn_mesh is not phi.bulk.mesh:
            raise ValueError(
                "StreamingOperator.apply_transpose: operator and composite "
                "must share the SAME SNMesh instance (mesh-identity invariant)."
            )
        # The adjoint matvec twin routes through the same polymorphic
        # :attr:`sweep_strategy` as the forward :meth:`apply`. CumprodScan
        # implements the 1-D reverse sweep (:meth:`_apply_1d_transpose`); the
        # multi-D Cartesian DAG strategies raise ``NotImplementedError`` (the
        # reverse sweep is deferred — O.2b lands the 1-D one first; the mesh
        # is compatible, only the adjoint feature is deferred).
        return self.sweep_strategy.residual_transpose(self, phi)

    def _apply_1d_transpose(self, phi: "TimedFullField") -> "TimedFullField":
        r"""1-D ``Lᵀ·φ`` (slab + sphere + cylinder) — the CumprodScan adjoint twin.

        ``(L+C)ᵀ·φ`` via :meth:`_MSpatialOperatorSum._compute_LpC_transpose`
        (the reverse-mode adjoint of the forward sweep), then ``− σ_t·φ`` to
        recover ``Lᵀ`` (Resolution A; ``C = σ_t⊙`` is a self-adjoint
        diagonal).  Extracted verbatim from :meth:`apply_transpose`'s 1-D
        branch by the SweepStrategy carve (S2) — bit-identical.  Returns the
        plain Euclidean transpose; the metric conjugation of the physical
        G-adjoint ``L.H`` is applied around this by ``_AdjointOperator``.
        """
        from orpheus.transport.source_sinks import AngularSourceSink
        from orpheus.transport.timed_full_field import TimedFullField

        LpCt = self.M_spatial._compute_LpC_transpose(phi)
        cell_values = (
            LpCt.bulk.values
            - self.sigma_t[None] * phi.bulk.values
        )
        return TimedFullField(
            bulk=AngularSourceSink.from_mesh(cell_values, self.sn_mesh),
            boundary=LpCt.boundary,
            _history=(),
            history_depth=phi.history_depth,
        )

    # ── SweepStrategy carve (S2) — the polymorphic matvec dispatch ─────

    @cached_property
    def sweep_strategy(self) -> "SweepStrategy":
        r"""The selected matvec/sweep strategy for this operator's mesh.

        The SAME first-class ``SweepStrategy``
        (``orpheus.sn.sweep_strategy``) that
        :func:`~orpheus.sn.sweep.transport_sweep` selects for the forward
        sweep — here it carries the matvec twin: :meth:`apply` routes through
        ``strategy.residual`` and :meth:`apply_transpose` through
        ``strategy.residual_transpose``.  Selection is by geometry
        (``default_for``): 1-D → ``CumprodScan``;
        2-D Cartesian → ``MovingFrontierWindow``.  ``cached_property`` because
        the selection is fixed by the mesh, stable across the operator's
        lifetime (mirrors :attr:`M_spatial` / :attr:`M_angular_redist`); the
        lazy import breaks the operator ↔ sweep_strategy module cycle.
        """
        from .sweep_strategy import default_for

        return default_for(self.sn_mesh)

    # ── Wave T T.4 — `M_spatial` / `M_angular_redist` properties ─────

    @cached_property
    def M_spatial(self) -> "_MSpatialOperatorSum":
        r"""Per-direction streaming + collision-share part of :math:`(L+C)`.

        :class:`OperatorSum` of two :class:`_SpatialSweepDirection`
        summands (forward-sweep μ_x > 0 and backward-sweep μ_x < 0).
        See :class:`_MSpatialOperatorSum` for the orchestrated apply
        semantics + Wave-O / adjoint / DSA leverage rationale.

        Slab status: ``M_spatial`` covers ALL of ``(L+C)``.
        Curvilinear (T.4c): ``M_spatial = (L+C) − M_angular_redist``.

        ``cached_property`` because the summand identities are stable
        across the operator's lifetime (factories of
        :class:`_SpatialSweepDirection`).
        """
        return _MSpatialOperatorSum(self.sn_mesh, self.sigma_t)

    @cached_property
    def M_angular_redist(self) -> "LinearOperator":
        r"""Curvilinear M-M angular redistribution part of :math:`(L+C)`.

        Slab (Cartesian): returns :class:`ZeroOperator` — there is no
        angular redistribution for the geometry-blind slab cell
        balance (``IdentityAngularClosure`` produces zero contributions).

        Curvilinear sphere/cylinder (T.4c future): returns the bespoke
        :class:`AngularRedistributionOperator` leaf wrapping the M-M
        half-grid recurrence.  The leaf is NOT a
        :class:`TensorProductOperator` — the half-grid is sequentially
        coupled per Hébert §3.9.4 (the recurrence on ``α_{m±1/2}``
        couples angular ordinates AND depends on σ_t).  Per the T.4
        verification spec MA-Q1 master condition, "coupled physics
        falls back to OperatorSum, not SOTP".

        ``cached_property`` because the leaf identity is stable across
        the operator's lifetime.
        """
        if getattr(self.sn_mesh, "curvature", None) is None:
            return ZeroOperator()
        # T.4c — curvilinear: bespoke AngularRedistributionOperator
        # leaf.  Wave T T.5 close-out (matvec retirement): the leaf
        # now takes a reference to `self.M_spatial` so its `.apply`
        # can hit the ψ-keyed dual-emission cache populated by the
        # spatial walk, sharing state at zero cost when
        # `StreamingOperator.apply` calls both `M_spatial.apply` and
        # `M_angular_redist.apply` within the same call boundary.
        return AngularRedistributionOperator(
            sn_mesh=self.sn_mesh,
            sigma_t=self.sigma_t,
            m_spatial=self.M_spatial,
        )

    def _apply_2d_cartesian(
        self, psi: "TimedFullField",
    ) -> "TimedFullField":
        r"""2-D Cartesian ``L·ψ`` via the diamond-difference closure (#208 O.4b).

        Routes through :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.residual`
        — the apply-direction walk of the SAME per-octant wavefront DAG and
        the SAME selectable ``CellUpdate`` (diamond-difference) closure the
        2-D sweep :func:`~orpheus.sn.sweep._sweep_2d_wavefront` uses. Matvec
        and sweep are therefore ONE discretization (L21: "matvec and sweep
        are different applications of the same operator"), not the FD/DD
        twin path the bespoke cell-centred upwind stencil
        (``_compute_gradients``, RETIRED in this change) created.

        The twin's failure mode — a face-trace formulation on the FD stencil
        converging to a non-uniform mode ~10% off ``k_inf`` — dissolves by
        construction: the DD closure handles the boundary cleanly and
        annihilates the flat reflective fundamental mode exactly
        (``L·ψ_uniform = 0``).

        Boundary semantics (BARE — O.4b Phase E BC-extraction). The boundary
        edge slots are seeded from the GIVEN inflow trace
        ``psi.boundary.face_view`` — there is NO ``bc.apply``. The octant's
        ordinates on its incoming face ARE its inflow ordinates, so the seed
        delivers exactly ``psi.boundary.inflow`` to the wavefront. The
        reflective coupling ``psi.inflow = B·psi.outflow`` is delivered by the
        sibling ``-B`` (:class:`~orpheus.sn.boundary_operator.SNBoundaryOperator`),
        a first-class coupling gain in the 2-D Krylov composed matvec (Wave O
        O.2a — ``_within_group_triple`` returns ``(L+C, S, B)``). The output
        boundary is the
        boundary-block residual (active trace), mirroring the 1-D
        ``L_full.apply`` template: OUTFLOW ordinate slots carry the
        self-consistency defect ``streamed − psi.outflow``; INFLOW ordinate
        slots carry the identity ``psi.inflow`` (so the composed
        ``(L+C−S−F−B)`` inflow residual is ``psi.inflow − B·psi.outflow``,
        driven to ``q.inflow`` by the outer loop). For homogeneous reflective
        the consistent uniform trace gives ``L·ψ_uniform = 0`` (recovers
        ``k_inf``) and a zero outflow defect — proven at the operator level by
        the E0 de-risk before this flip.

        Returns ``L·ψ`` (NOT ``(L+C)·ψ``): ``residual_batch`` at zero source
        yields ``(L+C)·ψ̄``; subtracting ``Σ_t·ψ̄`` (the collision term)
        gives the bare-streaming action, matching the 1-D path's L-only
        convention.
        """
        from orpheus.sn.sweep_graph import OctantLabel
        from orpheus.transport.source_sinks import AngularSourceSink, BoundarySourceSink
        from orpheus.transport.timed_full_field import TimedFullField

        sn_mesh = self.sn_mesh
        quad = sn_mesh.quad
        N = quad.N
        nx, ny = sn_mesh.nx, sn_mesh.ny
        ng = self.sigma_t.shape[0]
        sig_t = self.sigma_t                       # (ng, nx, ny)
        probe = psi.bulk.values                    # (N, ng, nx, ny) — the apply target ψ̄
        str_x = sn_mesh.streaming_x                # (N, nx)
        str_y = sn_mesh.streaming_y                # (N, ny)
        cell_update = sn_mesh.cell_update
        Q_zero = np.zeros((1, ng, nx, ny))         # matvec: no volumetric source

        # (L+C)·ψ̄ accumulator (the loss-operator action == residual_batch
        # at zero source); L·ψ̄ = this − Σ_t·ψ̄ at the end.
        LpC = np.zeros((N, ng, nx, ny))

        # Interior face cochain C¹_int, carried on the rolling _MovingFrontier
        # (Phase 5b storage-B) — the SAME windowed walk the 2-D sweep uses
        # (matvec ≡ sweep: ONE discretization, L21). BARE boundary handling
        # (O.4b Phase E): each octant reads its inflow from the GIVEN trace
        # ``psi.boundary.face_view(<incoming face>)[oct_idx]`` — NO bc.apply
        # (the octant's ordinates on its incoming face ARE its inflow
        # ordinates); the reflective coupling is the sibling -B. The post-walk
        # domain-edge outflow is shed per octant into ``streamed`` (only the
        # OUTFLOW ordinate slots are filled — the boundary residual reads only
        # those; grazing/inflow slots stay zero and are never read).
        trace = sn_mesh.trace
        boundary = psi.boundary
        streamed = {
            face: np.zeros_like(boundary.face_view(face))
            for face in trace.face_names
        }

        for octant in quad.octants:
            label_tuple = octant.label
            oct_idx = octant.indices
            sx = label_tuple[0] if len(label_tuple) >= 1 else +1
            sy = label_tuple[1] if len(label_tuple) >= 2 else 0
            # Pure-z degenerate octant: no in-plane streaming, so the
            # in-plane loss action is pure collision: (L+C)·ψ̄ = Σ_t·ψ̄
            # (hence L·ψ̄ = 0 for these ordinates after the subtraction).
            if sx == 0 and sy == 0:
                LpC[oct_idx] = sig_t * probe[oct_idx]
                continue
            sx_eff = +1 if sx == 0 else sx
            sy_eff = +1 if sy == 0 else sy
            graph = sn_mesh.sweep_graphs[OctantLabel((sx_eff, sy_eff))]

            # Octant domain-edge faces (incoming = low face of the sweep
            # direction). NO bc.apply — the seeded inflow IS the given trace.
            x_in_face = "xmin" if sx_eff >= 0 else "xmax"
            x_out_face = "xmax" if sx_eff >= 0 else "xmin"
            y_in_face = "ymin" if sy_eff >= 0 else "ymax"
            y_out_face = "ymax" if sy_eff >= 0 else "ymin"

            LpC_oct = np.zeros((oct_idx.size, ng, nx, ny))
            cap_x = np.empty((oct_idx.size, ng, ny))
            cap_y = np.empty((oct_idx.size, ng, nx))
            graph.residual_windowed(
                cell_update=cell_update,
                inflow_x=boundary.face_view(x_in_face)[oct_idx],
                inflow_y=boundary.face_view(y_in_face)[oct_idx],
                psi_avg_probe_octant=probe[oct_idx],
                Q_octant=Q_zero, sig_t=sig_t,
                str_x_octant=str_x[oct_idx], str_y_octant=str_y[oct_idx],
                residual_octant=LpC_oct, capture_x=cap_x, capture_y=cap_y,
            )
            LpC[oct_idx] = LpC_oct
            streamed[x_out_face][oct_idx] = cap_x
            streamed[y_out_face][oct_idx] = cap_y

        # L = (L+C) − C: subtract the collision term Σ_t·ψ̄ → bare streaming L·ψ̄.
        out_bulk = LpC - sig_t[None] * probe

        # Boundary-block residual (O.4b Phase E — the active trace, mirroring
        # the 1-D L_full.apply template). ``streamed[face]`` holds the
        # per-octant shed outflow at each face's OUTFLOW ordinate slots.
        # OUTFLOW ordinate slots → self-consistency defect ``streamed −
        # psi.outflow`` (kept as computed−stored so the vacuum path stays
        # bit-identical); INFLOW ordinate slots → identity ``psi.inflow`` (the
        # sibling -B adds -B·psi.outflow → composed inflow residual
        # ``psi.inflow − B·psi.outflow``).
        out_boundary = BoundarySourceSink.zeros_on(sn_mesh)
        for face in trace.face_names:
            given = boundary.face_view(face)
            out_idx = trace.outflow_indices_for_face(face)
            in_idx = trace.inflow_indices_for_face(face)
            if out_idx.size:
                out_boundary.face_view(face)[out_idx] = (
                    streamed[face][out_idx] - given[out_idx]
                )
            if in_idx.size:
                out_boundary.face_view(face)[in_idx] = given[in_idx]

        return TimedFullField(
            bulk=AngularSourceSink.from_mesh(out_bulk, sn_mesh),
            boundary=out_boundary,
            _history=(),
            history_depth=psi.history_depth,
        )

    def _apply_2d_cartesian_full_field(
        self, psi: "TimedFullField",
    ) -> "TimedFullField":
        r"""VERIFICATION ORACLE — the full-field storage-A 2-D matvec (NOT production).

        The matvec counterpart of :func:`~orpheus.sn.sweep._sweep_2d_full_field`:
        the storage-A path Phase 5b superseded in :meth:`_apply_2d_cartesian`
        (now the windowed `residual_windowed`). RETAINED + exercised by
        ``tests/sn/operators/test_2d_matvec_full_field_oracle.py`` as the
        fuller-view reference the moving-frontier matvec is cross-checked
        against END-TO-END. It carries the FULL interior cochain as a typed
        :class:`~orpheus.transport.fields.wavefront_flux.WavefrontFlux`
        (`seed`/`face`/`edge_view`), walks via
        :meth:`SweepDependencyGraph.residual` (the full-field walk sharing the
        SAME cell kernel as `residual_windowed`), and emits the identical
        boundary-block residual. Sole purpose: verification (production is the
        window); the shared kernel means the MATH cannot drift, only storage.
        """
        from orpheus.sn.sweep_graph import OctantLabel
        from orpheus.transport.fields.wavefront_flux import WavefrontFlux
        from orpheus.transport.source_sinks import AngularSourceSink, BoundarySourceSink
        from orpheus.transport.timed_full_field import TimedFullField

        sn_mesh = self.sn_mesh
        quad = sn_mesh.quad
        N = quad.N
        nx, ny = sn_mesh.nx, sn_mesh.ny
        ng = self.sigma_t.shape[0]
        sig_t = self.sigma_t
        probe = psi.bulk.values
        str_x = sn_mesh.streaming_x
        str_y = sn_mesh.streaming_y
        cell_update = sn_mesh.cell_update
        Q_zero = np.zeros((1, ng, nx, ny))

        LpC = np.zeros((N, ng, nx, ny))
        trace = sn_mesh.trace
        boundary = psi.boundary
        # The FULL interior cochain (typed) — ι_*-seeded whole-trace. The
        # per-axis face tuple is born from the cochain's OWN axis map
        # (``WavefrontFlux.face(a)`` over ``WavefrontFlux.axes``), not a
        # hardcoded ``psi_x = face(0); psi_y = face(1)``.
        wavefront = WavefrontFlux.zeros_on(sn_mesh)
        psi_faces = tuple(wavefront.face(a) for a in wavefront.axes)
        # Positional-by-axis, MUST match ``wavefront.axes`` order (see
        # ``_sweep_2d_full_field``); the d≥3 orchestrator (C3.4/C3.5) derives this
        # from an ``axes``-keyed streaming map so the two tuples cannot drift.
        str_axes = (str_x, str_y)
        wavefront.seed(boundary)

        for octant in quad.octants:
            label_tuple = octant.label
            oct_idx = octant.indices
            sx = label_tuple[0] if len(label_tuple) >= 1 else +1
            sy = label_tuple[1] if len(label_tuple) >= 2 else 0
            if sx == 0 and sy == 0:
                LpC[oct_idx] = sig_t * probe[oct_idx]
                continue
            sx_eff = +1 if sx == 0 else sx
            sy_eff = +1 if sy == 0 else sy
            graph = sn_mesh.sweep_graphs[OctantLabel((sx_eff, sy_eff))]
            psi_faces_oct = tuple(pf[oct_idx].copy() for pf in psi_faces)
            LpC_oct = np.zeros((oct_idx.size, ng, nx, ny))
            graph.residual(
                cell_update=cell_update,
                psi_faces_octant=psi_faces_oct,
                psi_avg_probe_octant=probe[oct_idx],
                Q_octant=Q_zero, sig_t=sig_t,
                str_axes_octant=tuple(s[oct_idx] for s in str_axes),
                residual_octant=LpC_oct,
            )
            for a in wavefront.axes:
                psi_faces[a][oct_idx] = psi_faces_oct[a]
            LpC[oct_idx] = LpC_oct

        out_bulk = LpC - sig_t[None] * probe
        # The post-walk domain-edge outflow (full interior cochain edge slots).
        streamed = {face: wavefront.edge_view(face) for face in trace.face_names}
        out_boundary = BoundarySourceSink.zeros_on(sn_mesh)
        for face in trace.face_names:
            given = boundary.face_view(face)
            out_idx = trace.outflow_indices_for_face(face)
            in_idx = trace.inflow_indices_for_face(face)
            if out_idx.size:
                out_boundary.face_view(face)[out_idx] = (
                    streamed[face][out_idx] - given[out_idx]
                )
            if in_idx.size:
                out_boundary.face_view(face)[in_idx] = given[in_idx]

        return TimedFullField(
            bulk=AngularSourceSink.from_mesh(out_bulk, sn_mesh),
            boundary=out_boundary,
            _history=(),
            history_depth=psi.history_depth,
        )

    # ── Algebra dispatch — sweep-invertible composite (R-1 Step C) ────

    def __add__(self, other):
        r"""Compose :math:`L + X`.

        When ``X`` is a :class:`CollisionOperator`, returns the
        sweep-invertible specialisation :class:`InvertibleOperator`
        carrying the algebraic identity :math:`(L + C)^{-1} \approx
        \text{WDD sweep}`.  Otherwise falls through to the generic
        :class:`OperatorSum` via the mixin.
        """
        if isinstance(other, CollisionOperator):
            return InvertibleOperator(self, other)
        return super().__add__(other)


@dataclass
class CollisionOperator(LinearOperatorMixin):
    r"""Pure collision operator :math:`C = \sigma\cdot` as a
    :class:`~orpheus.numerics.operator.LinearOperator` leaf.

    The "C" of the Phase G four-operator algebra
    :math:`A_{\rm wg} = L + C - S_{\rm foldable}`. Diagonal in position,
    group, and direction:

    .. math::

        (C\psi)_{n,i,j,g} \;=\; \sigma(i,j,g)\,\psi_{n,i,j,g}

    The supplied ``sigma`` is per-cell per-group only — the operator
    broadcasts identically across every ordinate ``n`` because the
    collision cross-section is direction-independent.

    Convention-agnostic σ
    ---------------------

    The same operator class supports the full total cross-section
    :math:`\sigma_t` AND the within-group removal cross-section
    :math:`\sigma_r = \sigma_t - \Sigma_{s,0}^{g\to g}`. The operator
    does NOT carry an interpretation flag — both quantities are
    ``(ng, nx, ny)`` arrays applied as :math:`\sigma\cdot\psi` (see
    :ref:`theory-sn-index-convention` for the canonical layout). The
    fusion hook (substep 3+4.b.ii) builds a lazy :math:`\sigma_r`
    :class:`CollisionOperator` from
    :meth:`~orpheus.sn.scattering.ScatteringOperator.foldable_sigma` at
    the moment it detects ``L + C - S.foldable_part()`` in an
    :class:`~orpheus.numerics.operator.OperatorSum`. (Pattern 4 —
    illegal states unrepresentable via the type, not a runtime flag
    the operator never inspects.)

    ⚠ A :math:`\sigma_r`-sweep is NOT a correct within-group self-scatter
    accelerator for anisotropic flux — it inverts the *diagonal-in-angle*
    removal :math:`\sigma_r\cdot I`, whereas the foldable self-scatter is the
    *isotropic-projection* operator :math:`\Sigma_{s,0}\,P_{\rm iso}`; the two
    agree only for isotropic ψ. Wiring ``A_wg.solve(S_residual)`` with this
    sweep ships 46–56 % errors on vacuum / heterogeneous problems (issue
    #215). The stable+correct fold is consistent DSA (#2) or Krylov. See the
    foldable/residual split note in :mod:`orpheus.sn.scattering`.

    Capability set
    --------------

    ``frozenset({CAP_APPLY, CAP_SOLVE, CAP_APPLY_TRANSPOSE})`` —
    collision is the simplest sort of operator. ``solve(q) = q / σ``
    is element-wise division; ``apply_transpose == apply`` because the
    operator is self-adjoint (diagonal in every basis the SN method
    operates in). All three capabilities are analytic.

    Parameters
    ----------
    sn_mesh : SNMesh
        The augmented geometry — carries mesh + quadrature + boundary
        operators (same as :class:`StreamingOperator`).
    sigma : np.ndarray
        Per-cell per-group cross-section, shape ``(ng, nx, ny)``
        (Issue #196 PR-INDEX-3). May be σ_t (full collision) or σ_r
        (removal — within-group self-scatter folded). The operator's
        action is identical.
    """

    sn_mesh: "SNMesh"
    sigma: np.ndarray

    capabilities: frozenset[str] = field(
        default_factory=lambda: frozenset(
            {CAP_APPLY, CAP_SOLVE, CAP_APPLY_TRANSPOSE}
        )
    )
    # Collision is a BULK operator — diagonal in (cell, group, ordinate),
    # no boundary action (A_bb only). Issue #208 / Wave O. Class-level
    # constant (unannotated so the dataclass does not treat it as a field).
    block_role = BlockRole.BULK

    # D-I.1 (2026-05-29): lazy ``_eq_map`` cache + ``_ensure_eq_map`` +
    # ``_sigma_at_unknowns`` retired together with the bare-ndarray
    # packed-vector apply / solve arms.  All consumers route through
    # :class:`TimedFullField`; no packed-vector decoder needed.

    def apply(self, psi: "TimedFullField") -> "TimedFullField":
        r"""Forward action :math:`C\,\psi = \sigma\cdot\psi` on the composite carrier.

        The diagonal action :math:`\sigma\cdot\psi` is per-cell per-group,
        broadcast across every ordinate.  Bulk receives
        ``σ ⊙ ψ.bulk.values``; boundary is the implicit-zero
        :class:`BoundaryFlux` — collision has no face-trace contribution
        (the cell-balance :math:`\sigma\cdot\psi` term is a CELL quantity;
        the boundary residual is a TRACE equation).  Option β3 / Issue #208
        will encode this bulk-only nature in the type via
        :class:`BulkOperator`.

        D-I.1 (2026-05-29) retired the legacy bare-ndarray packed-vector
        arm.  :class:`TimedFullField` is the sole accepted carrier.
        """
        from orpheus.transport.timed_full_field import TimedFullField
        from orpheus.transport.fields.angular_flux import AngularFlux
        from orpheus.transport.source_sinks import AngularSourceSink, BoundarySourceSink
        mesh = psi.bulk.mesh
        return TimedFullField(
            bulk=AngularSourceSink.from_mesh(
                self.sigma[None] * psi.bulk.values, mesh,
            ),
            boundary=BoundarySourceSink.zeros_on(mesh),
            _history=(),
            history_depth=psi.history_depth,
        )

    def solve(self, q: "TimedFullField") -> "TimedFullField":
        r"""Inverse action :math:`C^{-1}\,q = q/\sigma` on the composite carrier.

        Trivially invertible on the bulk: collision is diagonal, so the
        inverse is per-slot reciprocal scaling.  Returns NaN / Inf at
        slots where ``σ == 0`` per the IEEE-754 division contract —
        consumers constructing :math:`\sigma_r = \sigma_t -
        \Sigma_{s,0}^{g\to g}` must guarantee positivity (the operator
        does not check).

        Boundary is the implicit-zero :class:`BoundaryFlux` (Option β3,
        formal pseudoinverse on the rank-deficient face block — face
        slots of ``q`` are NOT inverted because collision contributes
        no volumetric term on the trace).

        D-I.1 retired the legacy bare-ndarray packed-vector arm.
        :class:`TimedFullField` is the sole accepted carrier.
        """
        from orpheus.transport.timed_full_field import TimedFullField
        from orpheus.transport.fields.angular_flux import AngularFlux
        from orpheus.transport.fields.boundary_flux import BoundaryFlux
        mesh = q.bulk.mesh
        return TimedFullField(
            bulk=AngularFlux.from_mesh(
                q.bulk.values / self.sigma[None], mesh,
            ),
            boundary=BoundaryFlux.zeros_on(mesh),
            _history=(),
            history_depth=q.history_depth,
        )

    def apply_transpose(self, psi: "TimedFullField") -> "TimedFullField":
        r"""Adjoint action :math:`C^*\,\psi = \sigma\cdot\psi`.

        Equal to :meth:`apply` — collision is self-adjoint (diagonal
        operator). Returned bit-equal to ``apply(psi)``.
        """
        return self.apply(psi)

    # ── Algebra dispatch — sweep-invertible composite (R-1 Step C) ────

    def __add__(self, other):
        r"""Compose :math:`C + X`.

        When ``X`` is a :class:`StreamingOperator`, returns the
        sweep-invertible specialisation :class:`InvertibleOperator`
        with the streaming operator placed first (the canonical
        ``L + C`` ordering for the algebraic identity).  Otherwise
        falls through to the generic :class:`OperatorSum` via the
        mixin.
        """
        if isinstance(other, StreamingOperator):
            return InvertibleOperator(other, self)
        return super().__add__(other)


# ─────────────────────────────────────────────────────────────────────────
# InvertibleOperator — sweep-invertible composite (L + C)
# ─────────────────────────────────────────────────────────────────────────


class InvertibleOperator(OperatorSum):
    r"""Sweep-invertible composite :math:`L + C` carrying ``.solve`` = WDD sweep.

    R-1 Step C (2026-05-19) — the SN-specific algebraic identity

    .. math::

        (L_{\rm streaming} + C_{\rm diagonal})^{-1} \;\approx\;
        \text{WDD sweep}

    has no generic ``(A+B)^{-1}`` formula — :class:`OperatorSum` by
    itself cannot ``solve``.  The WDD sweep IS the inverse algorithm
    for this specific composite — that's the algebraic foundation of
    the entire SN method (Lewis & Miller §3.2; Adams & Larsen 2002
    §III).  :class:`InvertibleOperator` is the specialisation that
    carries the identity at the type level: it inherits the
    :class:`OperatorSum` ``apply`` (the sum of the operand actions)
    and adds ``solve`` via :func:`~orpheus.sn.sweep.transport_sweep`.

    Construction
    ============

    Two equivalent paths:

    * **Operator algebra dispatch** — ``L + C`` where ``L`` is a
      :class:`StreamingOperator` and ``C`` is a
      :class:`CollisionOperator` returns this class automatically (see
      :meth:`StreamingOperator.__add__` and
      :meth:`CollisionOperator.__add__`).  The composite reads as math.
    * **Explicit construction** — ``InvertibleOperator(L, C)``.  Useful
      when composing variants such as
      ``InvertibleOperator(L_leaf, CollisionOperator(σ_r))`` where
      ``σ_r = σ_t - Σ_{s,0}^{g→g}`` is the removal cross-section that
      lets one fold the within-group self-scatter into the diagonal
      collision term (Adams & Larsen 2002 §III; tracked by issue
      `#200 <https://github.com/deOliveira-R/ORPHEUS/issues/200>`_).

    The two paths produce structurally identical objects — the choice
    only changes the call-site readability.

    Capability set
    ==============

    ``frozenset({CAP_APPLY, CAP_SOLVE})`` — adds ``solve`` to the
    parent :class:`OperatorSum`'s ``apply``-only set.
    ``apply_transpose`` is reserved for Phase H (adjoint propagation
    through the composite).

    The ``.solve`` API
    ==================

    The ``rhs`` parameter is a typed
    :class:`~orpheus.sn.angular_flux.AngularFlux` carrying:

    * ``rhs.values`` — per-ordinate source ``(N, ng, nx, ny)``.  This
      is treated as the per-ordinate anisotropic source
      :math:`Q^{\rm aniso}` that the sweep consumes (the isotropic
      source is zero).
    * ``rhs.boundary`` — face source / BC inflow trace.  Typically
      zero for volumetric SI/Krylov sources (which carry no face
      contribution); the persistent reflective-BC state lives on the
      :class:`SNMesh` and is handled inside the sweep.
    * ``rhs(1)`` (the lag-1 frame, when ``rhs.history_depth >= 2``) —
      the previous iterate that GENERATED this source.  Threaded as
      :pydata:`initial_guess` to :func:`transport_sweep` so the
      curvilinear sweep can read the Carlson coupled-pole seed from
      it.  ``None`` (cold start) → the sweep falls back to its
      in-iteration-source default.

    Parameters
    ----------
    streaming : StreamingOperator
        :math:`L = \Omega\cdot\nabla + \text{angular redistribution}`.
        Resolution A subtractive form:
        ``L.apply(ψ) = M(ψ; σ_t) - σ_t ⊙ ψ_cell``.
    diagonal : CollisionOperator
        :math:`C = \sigma\cdot`.  Its ``.sigma`` attribute is the
        per-cell per-group coefficient used by the sweep (canonically
        ``σ_t``; can be ``σ_r`` for the foldable variant).

    Notes
    -----
    The validation ``σ > 0`` everywhere guards against the
    ``σ_r < 0`` case that can arise when within-group self-scatter
    exceeds total cross-section (rare; not physically meaningful but
    mathematically possible for ill-conditioned multi-group sets).
    The sweep would emit NaN at those cells — surfacing the
    inconsistency at construction is friendlier.
    """

    def __init__(
        self,
        streaming: "StreamingOperator",
        diagonal: "CollisionOperator",
    ) -> None:
        if not isinstance(streaming, StreamingOperator):
            raise TypeError(
                f"InvertibleOperator: 'streaming' must be a "
                f"StreamingOperator; got {type(streaming).__name__}."
            )
        if not isinstance(diagonal, CollisionOperator):
            raise TypeError(
                f"InvertibleOperator: 'diagonal' must be a "
                f"CollisionOperator; got {type(diagonal).__name__}."
            )
        if streaming.sn_mesh is not diagonal.sn_mesh:
            raise ValueError(
                "InvertibleOperator: streaming and diagonal operators "
                "must share the same SNMesh instance "
                "(mesh-identity invariant)."
            )
        if not np.all(diagonal.sigma > 0):
            min_sigma = float(np.min(diagonal.sigma))
            raise ValueError(
                f"InvertibleOperator: diagonal coefficient must be "
                f"strictly positive everywhere for the WDD sweep to be "
                f"well-defined; got min(sigma) = {min_sigma:.3e}.  If "
                f"sigma_r = sigma_t - Sigma_(s,0)^(g->g) is dipping "
                f"negative, the multi-group cross-section set is "
                f"physically inconsistent."
            )
        super().__init__(streaming, diagonal)
        # OperatorSum.__init__ set capabilities = {CAP_APPLY, ...};
        # we add CAP_SOLVE because this composite IS sweep-invertible.
        self.capabilities = self.capabilities | frozenset({CAP_SOLVE})
        # block_role is now DERIVED by OperatorSum.__init__ (Wave O / O.2b 4.5):
        # join(L=FULL, C=BULK) = FULL. The former hand-stamp here was the
        # twin-path retired in 4.5 — the role is carried by construction.

    # ── Convenience accessors ─────────────────────────────────────────

    @property
    def streaming(self) -> "StreamingOperator":
        """The streaming operand (alias for ``self.a``)."""
        return self.a  # type: ignore[return-value]

    @property
    def diagonal(self) -> "CollisionOperator":
        """The diagonal-collision operand (alias for ``self.b``)."""
        return self.b  # type: ignore[return-value]

    @property
    def sn_mesh(self) -> "SNMesh":
        """The shared :class:`SNMesh` (validated mesh-identity at init)."""
        return self.streaming.sn_mesh

    @property
    def sigma(self) -> np.ndarray:
        r"""The diagonal coefficient used by ``solve`` (σ_t or σ_r)."""
        return self.diagonal.sigma

    # ── solve: WDD sweep ─────────────────────────────────────────────

    def solve(
        self,
        rhs: "AngularFlux | TimedFullField",
        *,
        initial_guess: "AngularFlux | TimedFullField | None" = None,
    ) -> "AngularFlux | TimedFullField":
        r"""Invert :math:`(L + C)\,\psi = \text{rhs}` via the WDD sweep.

        The cell-balance equation
        :math:`(\Omega\cdot\nabla + \sigma)\,\psi = Q` is integrated
        cell-by-cell in inflow-to-outflow order; the angular closure
        (Cartesian → identity, curvilinear → Morel-Montry) is bound
        on the mesh.

        Parameters
        ----------
        rhs : AngularFlux or TimedFullField
            Legacy :class:`AngularFlux` or composite
            :class:`TimedFullField` carrying:

            * ``values`` / ``bulk.values`` — per-ordinate source
              :math:`Q^{\rm aniso}`, shape ``(N, ng, nx, ny)``.
            * ``boundary`` / ``bulk.boundary`` — BC inflow trace
              (typically zero for SI/Krylov volumetric sources;
              falls back as the seed when no ``initial_guess`` is
              supplied).

            Both branches dispatch via runtime isinstance; the
            composite branch (D-H.1c stage 1) bridges through the
            legacy sweep path internally and returns
            :class:`TimedFullField`.  Output type matches input type.
        initial_guess : AngularFlux, TimedFullField, or None, optional
            Previous iterate :math:`\psi^{(k-1)}` for the curvilinear
            Carlson coupled-pole seed (M-M reads it as the level's
            angular flux to derive :math:`\bar Q` at :math:`\mu = -1`).
            Its ``.boundary`` also seeds the reflective-BC partner-flux
            trace when present (the previous outflow IS the partner-
            ordinate inflow).  ``None`` (default) → cold start: M-M's
            Carlson seed degenerates to the zero-input result (same
            as ``ZeroSeed`` ablation); BC trace falls back to
            ``rhs.boundary``.

            Type must match ``rhs`` (both legacy or both composite);
            mixed inputs raise :class:`TypeError`.

            Explicit kwarg (post-Phase-1.2) — the seed is no longer
            piggy-backed on ``rhs(1)`` (the AngularFlux lag-1 history
            is reserved for time-derivative tracking, an unrelated
            concern).  Outer iterators (:class:`SourceIteration`,
            :class:`KEigenvalue`) pass the previous iterate
            explicitly; GMRES residual calls pass ``None`` so the
            preconditioner doesn't silently route through stateful
            seed state (closes the R-1 Step D silent-fallback bug
            class).

        Returns
        -------
        AngularFlux or TimedFullField
            The angular flux satisfying :math:`(L + C)\,\psi =
            \text{rhs}`, with the sweep's outflow face state in
            ``.boundary`` and ``history_depth`` inherited from
            ``rhs.history_depth``.  Return type matches ``rhs`` input
            type.
        """
        return self._solve_timed_full_field(
            rhs, initial_guess=initial_guess,
        )

    def solve_moments(
        self,
        rhs: "TimedFullField",
        projection: "MomentProjection",
        *,
        initial_guess: "TimedFullField | None" = None,
    ) -> "TimedFullField":
        r"""Invert :math:`(L + C)\,\psi = \text{rhs}` and return the harmonic
        MOMENTS of :math:`\psi`, projected IN-SWEEP per anti-diagonal — the full
        per-ordinate angular field is never materialized.

        The Phase 5c moment-emitting sibling of :meth:`solve`: the SAME WDD
        sweep + boundary handling, but the bulk of the returned
        :class:`TimedFullField` is a
        :class:`~orpheus.transport.fields.harmonic_moment_field.HarmonicMomentField`
        ``(L+1, 2L+1, ng, nx, ny)`` rather than an
        :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
        ``(N, ng, nx, ny)``.  ``projection`` is the scattering operator's
        :class:`~orpheus.numerics.projection.MomentProjection` (its harmonics +
        weights), so the in-sweep moments equal ``S``'s internal projection
        term-for-term; the cross-octant accumulation reorders the ordinate sum
        vs the flat post-sweep projection ⇒ principled-equivalence, NOT
        bit-identity.  2-D Cartesian ONLY (the windowed-SI path; the windowing
        gate ``sn_mesh.reduced is None`` guarantees it).  ``solve`` followed by
        the flat :meth:`MomentProjection.apply` is the fuller-view verification
        oracle (``vv-principles``; the aggressive-retirement "verification
        oracle" exception).
        """
        return self._solve_timed_full_field(
            rhs, initial_guess=initial_guess, moment_projection=projection,
        )

    def _solve_timed_full_field(
        self,
        rhs: "TimedFullField",
        *,
        initial_guess: "TimedFullField | None" = None,
        moment_projection: "MomentProjection | None" = None,
    ) -> "TimedFullField":
        r"""Composite :class:`TimedFullField` body of :meth:`solve` (D-H.1c stage 1).

        Bridges through the legacy :class:`AngularFlux` solve path —
        the WDD sweep kernel (:func:`~orpheus.sn.sweep.transport_sweep`)
        stays untouched; this method handles only the L2↔legacy
        bridge at the public-entry boundary.

        The boundary plumbing for the reflective-BC partner-flux state
        is preserved: ``initial_guess.boundary`` (composite) → legacy
        BoundaryFlux via :meth:`TimedFullField.to_legacy_angular_flux`
        → ``boundary_buf`` (the sweep's mutable write-through buffer)
        via :func:`_copy_boundary_face_state`.  The sweep mutates
        ``boundary_buf`` in place; the result is re-wrapped as an L2
        composite at the end.

        Parameters
        ----------
        rhs : TimedFullField
            Per-ordinate source on the composite carrier.
        initial_guess : TimedFullField or None
            Previous iterate (carries the partner-flux trace and
            curvilinear Carlson seed).  ``None`` → cold start.

        Returns
        -------
        TimedFullField
            Solve output with ``bulk`` = ``(L + C)^{-1} rhs.bulk`` and
            ``boundary`` = the sweep's outflow face state.
            ``history_depth`` matches ``rhs.history_depth``; ``_history``
            is empty (solver outputs carry no iteration history — the
            outer SI / Krylov loop owns history threading).
        """
        from orpheus.transport.fields.angular_flux import (
            AngularFlux,
        )
        from orpheus.transport.fields.boundary_flux import (
            BoundaryFlux,
        )
        from orpheus.transport.fields.harmonic_moment_field import (
            HarmonicMomentField,
        )
        from orpheus.transport.source_sinks import AngularSourceSink
        from orpheus.transport.timed_full_field import TimedFullField
        from .sweep import transport_sweep

        # D-H.2-C3: only the :class:`TimedFullField` composite branch remains;
        # legacy :class:`AngularFlux` retired.  ``rhs`` and ``initial_guess``
        # MUST be :class:`TimedFullField` (or ``None`` for ``initial_guess``).
        # Single guard site for both :meth:`solve` and :meth:`solve_moments`.
        if not isinstance(rhs, TimedFullField):
            raise TypeError(
                f"InvertibleOperator: 'rhs' must be TimedFullField; "
                f"got {type(rhs).__name__}.  Legacy AngularFlux retired "
                f"in D-H.2-C3."
            )
        if initial_guess is not None and not isinstance(
            initial_guess, TimedFullField,
        ):
            raise TypeError(
                f"InvertibleOperator: 'initial_guess' must be "
                f"TimedFullField or None; got "
                f"{type(initial_guess).__name__}."
            )

        sn_mesh = self.sn_mesh
        if rhs.bulk.mesh is not sn_mesh:
            raise ValueError(
                "InvertibleOperator.solve(TimedFullField): rhs and "
                "operator must share the same SNMesh instance "
                "(mesh-identity invariant)."
            )
        if initial_guess is not None and initial_guess.bulk.mesh is not sn_mesh:
            raise ValueError(
                "InvertibleOperator.solve(TimedFullField): initial_guess "
                "and operator must share the same SNMesh instance "
                "(mesh-identity invariant)."
            )

        # ── L2 boundary buffer for the sweep (D-H.2-C2) ───────────────
        #
        # The sweep mutates ``boundary_buf`` (the L2 mutable
        # write-through; ``frozen=True`` freezes field rebinding but
        # the underlying flat ndarray remains writable through
        # :meth:`face_view`).
        #
        # Wave O (#208) O.4a.2 — BARE SWEEP: the inflow seed is the
        # boundary SOURCE ``rhs.boundary`` (the inflow slots carry
        # ``q.boundary + B·ψ.outflow`` — the SI driver folds ``S + B`` so
        # the ``Bψ`` reflective inflow rides in ``rhs.boundary``).  This
        # REPLACES the pre-extraction seed-from-``initial_guess.boundary``:
        # the bare sweep no longer re-applies ``bc`` to the iterate's
        # outflow, so the iterate's boundary is NOT the inflow seed.  The
        # iterate (``initial_guess``) still threads the BULK Carlson /
        # angular warm-start through ``transport_sweep`` below — that path
        # reads ``initial_guess.bulk``, not its boundary.
        boundary_buf = BoundaryFlux.zeros_on(sn_mesh)  # L2 after C2
        seed_boundary = rhs.boundary
        # Per-face copy via L2 face_view — works for slab (xmin, xmax),
        # curvilinear (xmax only), and 2-D Cartesian (all 4).
        for face_name in boundary_buf.layout.faces:
            if face_name in seed_boundary.layout.faces:
                boundary_buf.face_view(face_name)[:] = (
                    seed_boundary.face_view(face_name)
                )

        # ── Per-ordinate source from rhs.bulk (single-source convention
        # per R-1 Step 4 A1 — ``rhs.bulk.values`` IS per-ordinate
        # density by producer contract).
        source = AngularSourceSink.from_mesh(rhs.bulk.values, sn_mesh)

        # ── Sweep: pass the composite ``initial_guess`` directly.
        # D-H.1c stage 4: :func:`transport_sweep` accepts both legacy
        # AngularFlux and TimedFullField for ``initial_guess`` (via the
        # container-agnostic :func:`_initial_guess_values` extractor in
        # sweep.py).  The kernel reads ``.bulk.values`` for the composite
        # path with no AngularFlux round-trip.
        #
        # Phase 5c: ONE sweep through :func:`transport_sweep` for BOTH output
        # modes — the moment projection rides as an optional kwarg
        # (``transport_sweep`` forwards it to the 2-D wavefront sweep and raises
        # on a 1-D mesh, since moment output is 2-D Cartesian only; ``moment_*``
        # and ``initial_guess`` are mutually 2-D-vs-1-D, so the 2-D branch
        # harmlessly drops the unused seed).  Only the OUTPUT WRAP differs: the
        # full angular field vs the harmonic moment tensor.
        bulk_values, _scalar = transport_sweep(
            source,
            self.sigma,
            sn_mesh,
            boundary_buf,
            initial_guess=initial_guess,
            moment_projection=moment_projection,
        )
        if moment_projection is None:
            bulk = AngularFlux.from_mesh(bulk_values, sn_mesh)
        else:
            bulk = HarmonicMomentField.from_mesh_and_L(
                bulk_values, sn_mesh, moment_projection.L,
            )

        # ── L2 direct return — no adapter needed (D-H.2-C2). ───────────
        return TimedFullField(
            bulk=bulk,
            boundary=boundary_buf,
            _history=(),
            history_depth=rhs.history_depth,
        )
