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
   (:meth:`StreamingOperator._apply_2d_cartesian`) is NOT yet bare — it
   still routes the incoming-at-boundary slots through the
   :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` instances on
   the :class:`~orpheus.sn.geometry.SNMesh` (``bc_xmin`` / ``bc_xmax``
   / ``bc_ymin`` / ``bc_ymax``) via ``bc.apply(outgoing)`` inside the
   sweep (deferred to O.4b).  The pre-extraction Phase C insight that
   the BC must consume the WDD-propagated outflow face vector (not
   cell centres) is preserved and strengthened: post-O.4a.2 the
   outflow trace is the explicit solved unknown ``psi.outflow`` that
   ``-B`` reads, closing ERR-026 by construction for the 1-D
   curvilinear path.
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
    from orpheus.transport.source_sinks import ScalarSourceSink, AngularSourceSink
    from .spatial.pole_angular_closure import PoleAngularClosure

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
# Finite-difference gradients (diamond scheme with reflective BCs)
# ═══════════════════════════════════════════════════════════════════════

def _compute_gradients(
    fi: np.ndarray,
    n: int, ix: int, iy: int,
    quad: AngularQuadrature,
    nx: int, ny: int,
    dx: np.ndarray, dy: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Upwind cell-center gradients with reflective BCs.

    Issue #196 PR-INDEX-7: consumes ``fi`` in principled
    ``(N, ng, nx, ny)`` layout. Per-(n, ix, iy) slice
    ``fi[n, :, ix, iy]`` is shape ``(ng,)``.

    Returns (dfi/dx, dfi/dy), each shape (ng,).

    The gradient between adjacent cell centers is divided by the
    **cell-center distance** ``(dx[i] + dx[j]) / 2`` rather than
    the local cell width ``dx[i]``. On a uniform mesh these are
    identical; on a non-uniform mesh the cell-center distance is
    the correct denominator for a first-order consistent FD stencil.
    The original MATLAB code (Mikityuk, PSI 2015) used a scalar
    ``g.delta`` (uniform mesh only), so this distinction did not arise.
    """
    ref_x = quad.reflection_index("x")
    ref_y = quad.reflection_index("y")
    mu_x, mu_y = quad.mu_x, quad.mu_y

    # X gradient
    if mu_x[n] > 1e-15:
        if ix == 0:
            dfix = fi[ref_x[n], :, ix, iy] - fi[ref_x[n], :, ix + 1, iy]
            hx = 0.5 * (dx[ix] + dx[ix + 1])
        else:
            dfix = fi[n, :, ix, iy] - fi[n, :, ix - 1, iy]
            hx = 0.5 * (dx[ix] + dx[ix - 1])
    elif mu_x[n] < -1e-15:
        if ix == nx - 1:
            dfix = fi[ref_x[n], :, ix - 1, iy] - fi[ref_x[n], :, ix, iy]
            hx = 0.5 * (dx[ix - 1] + dx[ix])
        else:
            dfix = fi[n, :, ix + 1, iy] - fi[n, :, ix, iy]
            hx = 0.5 * (dx[ix + 1] + dx[ix])
    else:
        # PR-INDEX-7: ng is fi.shape[1] (axis 1 = groups under principled layout).
        dfix = np.zeros(fi.shape[1])
        hx = 1.0

    # Y gradient
    if mu_y[n] > 1e-15:
        if iy == 0:
            dfiy = fi[ref_y[n], :, ix, iy] - fi[ref_y[n], :, ix, iy + 1]
            hy = 0.5 * (dy[iy] + dy[iy + 1])
        else:
            dfiy = fi[n, :, ix, iy] - fi[n, :, ix, iy - 1]
            hy = 0.5 * (dy[iy] + dy[iy - 1])
    elif mu_y[n] < -1e-15:
        if iy == ny - 1:
            dfiy = fi[ref_y[n], :, ix, iy - 1] - fi[ref_y[n], :, ix, iy]
            hy = 0.5 * (dy[iy - 1] + dy[iy])
        else:
            dfiy = fi[n, :, ix, iy + 1] - fi[n, :, ix, iy]
            hy = 0.5 * (dy[iy + 1] + dy[iy])
    else:
        dfiy = np.zeros(fi.shape[1])
        hy = 1.0

    return dfix / hx, dfiy / hy

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
        from orpheus.transport.source_sinks import AngularSourceSink
        from orpheus.transport.fields.boundary_flux import BoundaryFlux
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
        masked_bulk[~ordinate_mask, :, :, :] = 0.0

        # Boundary residual: forward sweep writes the outer-face WDD
        # outflow residual; backward sweep writes the inner-face
        # residual (for slab).  The per-direction split exposes which
        # face each direction's contribution lives on.
        masked_boundary = BoundaryFlux.zeros_on(sn_mesh)
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
        from orpheus.transport.source_sinks import AngularSourceSink
        from orpheus.transport.fields.boundary_flux import BoundaryFlux
        from orpheus.transport.timed_full_field import TimedFullField
        from .spatial.cell_balance import cell_balance_for_streaming
        from .spatial.pole_angular_closure import MorelMontryAngularSweep

        sn_mesh = self.sn_mesh
        psi_view = psi.bulk.values
        quad = sn_mesh.quad
        N = quad.N
        ng = psi_view.shape[1]
        nx = sn_mesh.nx
        ny = sn_mesh.ny
        eps = 1e-15
        curvature_raw = getattr(sn_mesh, "curvature", None)
        curvature = curvature_raw if curvature_raw is not None else "cartesian"

        if curvature not in ("spherical", "cylindrical", "cartesian"):
            raise ValueError(f"Unknown curvature: {curvature!r}")
        if curvature == "cartesian" and ny > 1:
            raise NotImplementedError(
                "_MSpatialOperatorSum._compute_LpC: 2-D Cartesian "
                "is not yet wired through dag_walk; only 1-D slab (ny=1) "
                "is implemented.  2-D Cartesian routes through "
                "`StreamingOperator._apply_2d_cartesian` (Q1 hybrid)."
            )

        pole_angular_closure = sn_mesh.pole_angular_closure
        if pole_angular_closure is None and curvature != "cartesian":
            pole_angular_closure = MorelMontryAngularSweep()

        mu_x = quad.mu_x
        level_indices: tuple[np.ndarray, ...] = pole_angular_closure.level_indices
        A = sn_mesh.areas

        psi_g_first = psi_view.transpose(1, 0, 2, 3)
        out_g_first = np.zeros((ng, N, nx, ny))

        V = sn_mesh.volumes[:, 0]
        sigma_t_gx = self.sigma_t[:, :, 0]

        boundary = psi.boundary
        trace = sn_mesh.trace
        has_inner_face = "xmin" in boundary.layout.faces
        face_outer = boundary.face_view("xmax")
        face_inner = boundary.face_view("xmin") if has_inner_face else None

        if curvature != "cartesian":
            # Curvilinear: the pole seed is the r=0 REGULARITY condition
            # (NOT a boundary condition) — read the innermost cell flux.
            pole_face_seed = psi_view[:, :, 0, 0].copy()
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
                    psi_cell = psi_g_first[:, global_dir, i, 0]
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
                    out_g_first[:, global_dir, i, 0] = m_full
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

                psi_cell = psi_g_first[:, global_deg, i, 0]
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
                out_g_first[:, global_deg, i, 0] = m_full

        m_cell = out_g_first.transpose(1, 0, 2, 3)

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
        m_boundary = BoundaryFlux.zeros_on(sn_mesh)
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
        from orpheus.transport.source_sinks import AngularSourceSink
        from orpheus.transport.fields.boundary_flux import BoundaryFlux
        from orpheus.transport.timed_full_field import TimedFullField
        from .spatial.cell_balance import cell_balance_for_streaming
        from .spatial.pole_angular_closure import MorelMontryAngularSweep

        sn_mesh = self.sn_mesh
        psi_view = psi.bulk.values                                       # (N, ng, nx, ny)
        quad = sn_mesh.quad
        N = quad.N
        ng = psi_view.shape[1]
        nx = sn_mesh.nx
        ny = sn_mesh.ny
        eps = 1e-15
        curvature_raw = getattr(sn_mesh, "curvature", None)
        curvature = curvature_raw if curvature_raw is not None else "cartesian"

        if curvature not in ("spherical", "cylindrical", "cartesian"):
            raise ValueError(f"Unknown curvature: {curvature!r}")
        if curvature == "cartesian" and ny > 1:
            raise NotImplementedError(
                "_MSpatialOperatorSum._compute_decomposition: 2-D Cartesian "
                "is not yet wired through dag_walk; only 1-D slab (ny=1) "
                "is implemented.  2-D Cartesian routes through "
                "`StreamingOperator._apply_2d_cartesian` (Q1 hybrid)."
            )

        pole_angular_closure = sn_mesh.pole_angular_closure
        if pole_angular_closure is None and curvature != "cartesian":
            pole_angular_closure = MorelMontryAngularSweep()

        mu_x = quad.mu_x
        level_indices: tuple[np.ndarray, ...] = pole_angular_closure.level_indices
        A = sn_mesh.areas                                                # (nx+1,)

        psi_g_first = psi_view.transpose(1, 0, 2, 3)                     # (ng, N, nx, ny)
        out_spat_g_first = np.zeros((ng, N, nx, ny))
        out_ang_g_first = np.zeros((ng, N, nx, ny))

        V = sn_mesh.volumes[:, 0]                                        # (nx,)
        sigma_t_gx = self.sigma_t[:, :, 0]                               # (ng, nx)

        boundary = psi.boundary
        trace = sn_mesh.trace
        has_inner_face = "xmin" in boundary.layout.faces
        face_outer = boundary.face_view("xmax")                          # (N, ng)
        face_inner = (
            boundary.face_view("xmin") if has_inner_face else None
        )

        if curvature != "cartesian":
            # Curvilinear: pole seed = r=0 REGULARITY condition (not a BC).
            pole_face_seed = psi_view[:, :, 0, 0].copy()                 # (N, ng)
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
                    psi_cell = psi_g_first[:, global_dir, i, 0]
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
                    out_spat_g_first[:, global_dir, i, 0] = m_spat
                    out_ang_g_first[:, global_dir, i, 0] = m_ang

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

                psi_cell = psi_g_first[:, global_deg, i, 0]              # (ng, n_deg)
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
                out_spat_g_first[:, global_deg, i, 0] = m_spat
                out_ang_g_first[:, global_deg, i, 0] = m_ang

        m_spat_cell = out_spat_g_first.transpose(1, 0, 2, 3)             # (N, ng, nx, ny)
        m_ang_cell = out_ang_g_first.transpose(1, 0, 2, 3)               # (N, ng, nx, ny)

        # M_spatial carries the face residuals (only the spatial sweep
        # writes them; per MA-Q4 M_angular_redist is a BulkOperator).
        # Outflow set read from the unified TraceSpace selector (single
        # source of truth for sign(Ω·n) — see A.4 in _compute_LpC).
        # Wave O O.4a.2 (mirror of _compute_LpC): the (L+C) boundary block
        # carries the two trace diagonals — outflow slots = self-consistency
        # defect ψ.outflow − streamed (kept; vacuum bit-identical), inflow
        # slots = identity ψ.inflow (the r_inflow diagonal; the sibling −B adds
        # −B·ψ.outflow). M_angular_redist stays zero-boundary (a BulkOperator).
        m_spat_boundary = BoundaryFlux.zeros_on(sn_mesh)
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
            boundary=BoundaryFlux.zeros_on(sn_mesh),
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
        sig_t_1d = self.sigma_t[:, :, 0]
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
    - S_foldable).solve(q)`` routes to the within-group sweep via the
    fusion hook substep 3+4.b.ii adds. No ``apply_transpose`` yet —
    the analytic reverse-direction sweep is Step 6 work.

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
        default_factory=lambda: frozenset({CAP_APPLY})
    )
    # Streaming is the sole FULL operator — it couples bulk ↔ boundary
    # (reads the inflow trace to seed the sweep, writes the outflow
    # trace). Issue #208 / Wave O. Class-level constant (unannotated so
    # the dataclass does not treat it as a field).
    block_role = BlockRole.FULL

    # D-J (2026-05-30): ``_eq_map`` / ``_ensure_eq_map`` / ``n_unknowns``
    # retired alongside the :class:`EquationMap` codec family — the
    # typed :class:`~orpheus.transport.timed_full_field.TimedFullField`
    # contract has no need for the legacy packed-vector slot map.

    def apply(self, psi: "TimedFullField") -> "TimedFullField":
        r"""Subtractive forward action :math:`L\,\psi = M(\psi;\sigma_t)
        - \sigma_t \odot \psi.bulk`.

        Routes through the geometry-agnostic
        :func:`_transport_operator_matvec_unified` for 1-D slab, sphere,
        and cylinder; through :meth:`_apply_2d_cartesian` (L2-native
        FD kernel) for 2-D Cartesian.

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
        from orpheus.transport.fields.angular_flux import AngularFlux
        from orpheus.transport.source_sinks import AngularSourceSink
        from orpheus.transport.timed_full_field import TimedFullField

        if not isinstance(psi, TimedFullField):
            raise TypeError(
                "StreamingOperator.apply: expected TimedFullField, got "
                f"{type(psi).__name__}.  D-I.3d (2026-05-29) retired the "
                "bare-ndarray packed-vector contract; construct a typed "
                "composite via ``TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh)`` or "
                "explicit ``TimedFullField(bulk=..., boundary=...)``."
            )

        sn_mesh = self.sn_mesh
        if sn_mesh is not psi.bulk.mesh:
            raise ValueError(
                "StreamingOperator.apply: operator and composite must "
                "share the SAME SNMesh instance (mesh-identity invariant)."
            )

        curv = getattr(sn_mesh, "curvature", None)
        if curv is None and sn_mesh.ny > 1:
            # 2-D Cartesian — D-H.2-C4d L2-native FD kernel.
            # T.4 Q1 = hybrid: this path stays procedural; A2D-1
            # source-hash defensive pin guards against silent edits.
            return self._apply_2d_cartesian(psi)

        # 1-D unified path (slab + sphere + cylinder) — Wave T post-T.5
        # matvec retirement.  Production hot path uses
        # :meth:`_MSpatialOperatorSum._compute_LpC` (single-emission
        # legacy matvec body, inlined as a class method) for
        # perf-optimal ``(L+C)·ψ`` computation; ``StreamingOperator.apply``
        # then subtracts ``σ_t·ψ`` to recover L per Resolution A.
        #
        # The M_spatial / M_angular_redist algebra decomposition is
        # exposed for introspection via :meth:`_compute_decomposition`
        # (slower; dual-emission), used by the per-property apply
        # paths in ``M_spatial.apply`` and ``M_angular_redist.apply``
        # — NOT on the production hot path because re-splitting and
        # re-summing is wasted work when ``StreamingOperator.apply``
        # just needs (L+C).
        LpC_result = self.M_spatial._compute_LpC(psi)
        cell_values = (
            LpC_result.bulk.values
            - self.sigma_t[None, :, :, :] * psi.bulk.values
        )
        return TimedFullField(
            bulk=AngularSourceSink.from_mesh(cell_values, sn_mesh),
            boundary=LpC_result.boundary,
            _history=(),
            history_depth=psi.history_depth,
        )

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
        r"""2-D Cartesian L2-native FD matvec (D-H.2-C4d).

        Computes :math:`L\,\psi = (\mu_x\partial_x + \mu_y\partial_y)\psi`
        with cell-centred upwind FD (the legacy stencil from
        :func:`transport_operator_matvec`'s body), wrapped in the L2
        ``TimedFullField`` carrier.

        Boundary semantics
        ------------------

        Legacy cell-centre-proxy semantics: the matvec body reads
        ``psi.bulk.values[:, :, 0, iy]`` as the outgoing trace at
        xmin (and similarly at other faces).  The BC's
        ``apply(outgoing)`` returns the incoming-direction values;
        the incoming-direction values fill the boundary cells of
        ``fi`` (overwriting the bulk value at incoming positions).
        For homogeneous reflective, the cell-centre proxy makes the
        kernel reduce to a uniform ``fi``, giving ``L·ψ_uniform = 0``
        and converging the eigenvalue to ``k_inf``.

        ``psi.boundary.face_view`` is currently passive: its values
        do NOT enter the bulk computation.  The output boundary is
        the zero L2 :class:`BoundaryFlux`.  Krylov drives the cell
        residual to zero on the bulk dimension only; the face_view
        is left at whatever value the iteration produces, but does
        not affect convergence.  This matches the legacy 2-D
        ``transport_operator_matvec`` semantics — the Krylov problem
        size is N × ng × nx × ny cells (no face unknowns).

        A more ambitious face_view-as-trace formulation (face_view
        enters the bulk computation as the boundary trace, with a
        boundary residual driving face_view ↔ bulk consistency)
        causes the eigenvalue iteration to converge to a non-uniform
        mode (~10% off from k_inf).  Deferred to a future Wave T /
        TensorProduct refactor when the BC realizers gain a proper
        composable algebra.

        Returns ``L·ψ`` (NOT ``(L+C)·ψ`` — the σ_t·ψ term subtracts
        out at the cell level, matching the 1-D path's convention
        that ``_apply_typed`` / ``_apply_timed_full_field`` return
        L-only).
        """
        from orpheus.geometry.boundary import ReflectiveBoundary
        from orpheus.transport.fields.angular_flux import (
            AngularFlux,
        )
        from orpheus.transport.source_sinks import AngularSourceSink
        from orpheus.transport.fields.boundary_flux import (
            BoundaryFlux,
        )
        from orpheus.transport.timed_full_field import TimedFullField

        sn_mesh = self.sn_mesh
        quad = sn_mesh.quad
        N = quad.N
        nx, ny = sn_mesh.nx, sn_mesh.ny
        dx, dy = sn_mesh.dx, sn_mesh.dy
        mu_x, mu_y = quad.mu_x, quad.mu_y
        trace = sn_mesh.trace

        bc_xmin = getattr(sn_mesh, "bc_xmin", None) or ReflectiveBoundary(
            axis="x", albedo=1.0,
        )
        bc_xmax = getattr(sn_mesh, "bc_xmax", None) or ReflectiveBoundary(
            axis="x", albedo=1.0,
        )
        bc_ymin = getattr(sn_mesh, "bc_ymin", None) or ReflectiveBoundary(
            axis="y", albedo=1.0,
        )
        bc_ymax = getattr(sn_mesh, "bc_ymax", None) or ReflectiveBoundary(
            axis="y", albedo=1.0,
        )

        # ── Build fi: cell-centre proxy + BC-filled incoming at boundary ─
        fi = psi.bulk.values.copy()

        # Incoming-ordinate set per face, read from the unified TraceSpace
        # selector (single source of truth for sign(Ω·n) — A.4 retired
        # the inline ``mu_{x,y} ≷ ±eps`` masks this 2-D matvec used to
        # recompute; ``inflow`` is Ω·n < −ε, i.e. the direction points
        # INTO the domain through that face).
        xmin_inflow = trace.inflow_indices_for_face("xmin")
        xmax_inflow = trace.inflow_indices_for_face("xmax")
        ymin_inflow = trace.inflow_indices_for_face("ymin")
        ymax_inflow = trace.inflow_indices_for_face("ymax")

        # xmin / xmax: outgoing trace = fi at boundary cell; BC.apply
        # returns full (N, ng) incoming; write only the incoming cells.
        for iy in range(ny):
            incoming_xmin = bc_xmin.apply(fi[:, :, 0, iy])      # (N, ng)
            fi[xmin_inflow, :, 0, iy] = incoming_xmin[xmin_inflow]
            incoming_xmax = bc_xmax.apply(fi[:, :, -1, iy])
            fi[xmax_inflow, :, -1, iy] = incoming_xmax[xmax_inflow]

        for ix in range(nx):
            incoming_ymin = bc_ymin.apply(fi[:, :, ix, 0])
            fi[ymin_inflow, :, ix, 0] = incoming_ymin[ymin_inflow]
            incoming_ymax = bc_ymax.apply(fi[:, :, ix, -1])
            fi[ymax_inflow, :, ix, -1] = incoming_ymax[ymax_inflow]

        # ── Compute M·ψ = (L+C)·ψ via cell-centred FD stencil ─────────
        out_M = np.zeros_like(psi.bulk.values)
        for n in range(N):
            for ix in range(nx):
                for iy in range(ny):
                    dfix, dfiy = _compute_gradients(
                        fi, n, ix, iy, quad, nx, ny, dx, dy,
                    )
                    out_M[n, :, ix, iy] = (
                        mu_x[n] * dfix
                        + mu_y[n] * dfiy
                        + self.sigma_t[:, ix, iy] * fi[n, :, ix, iy]
                    )

        # L = M - C: subtract σ_t · ψ at the cell-centres (using
        # ORIGINAL psi.bulk.values, not the BC-filled fi).
        out_bulk = out_M - self.sigma_t[None, :, :, :] * psi.bulk.values

        # Boundary output: zero — face_view is passive in this 2-D
        # cell-centre-proxy formulation (see method docstring).
        out_boundary = BoundaryFlux.zeros_on(sn_mesh)

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
        from orpheus.transport.source_sinks import AngularSourceSink
        from orpheus.transport.fields.boundary_flux import BoundaryFlux
        mesh = psi.bulk.mesh
        return TimedFullField(
            bulk=AngularSourceSink.from_mesh(
                self.sigma[None, :, :, :] * psi.bulk.values, mesh,
            ),
            boundary=BoundaryFlux.zeros_on(mesh),
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
                q.bulk.values / self.sigma[None, :, :, :], mesh,
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
        # (L + C): FULL inherits from the streaming operand (the bulk
        # collision C carries no boundary action). Issue #208 / Wave O.
        self.block_role = BlockRole.FULL

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
        from orpheus.transport.timed_full_field import TimedFullField

        # D-H.2-C3: only the :class:`TimedFullField` composite branch
        # remains; legacy :class:`AngularFlux` retired.  ``rhs`` and
        # ``initial_guess`` MUST be :class:`TimedFullField` (or ``None``
        # for ``initial_guess``).
        if not isinstance(rhs, TimedFullField):
            raise TypeError(
                f"InvertibleOperator.solve: 'rhs' must be TimedFullField; "
                f"got {type(rhs).__name__}.  Legacy AngularFlux retired "
                f"in D-H.2-C3."
            )
        if initial_guess is not None and not isinstance(
            initial_guess, TimedFullField,
        ):
            raise TypeError(
                f"InvertibleOperator.solve: 'initial_guess' must be "
                f"TimedFullField or None; got "
                f"{type(initial_guess).__name__}."
            )
        return self._solve_timed_full_field(
            rhs, initial_guess=initial_guess,
        )

    def _solve_timed_full_field(
        self,
        rhs: "TimedFullField",
        *,
        initial_guess: "TimedFullField | None" = None,
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
        from orpheus.transport.source_sinks import AngularSourceSink
        from orpheus.transport.timed_full_field import TimedFullField
        from .sweep import transport_sweep

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
        angular, _scalar = transport_sweep(
            source,
            self.sigma,
            sn_mesh,
            boundary_buf,
            initial_guess=initial_guess,
        )

        # ── L2 direct return — no adapter needed (D-H.2-C2). ───────────
        return TimedFullField(
            bulk=AngularFlux.from_mesh(angular, sn_mesh),
            boundary=boundary_buf,
            _history=(),
            history_depth=rhs.history_depth,
        )
