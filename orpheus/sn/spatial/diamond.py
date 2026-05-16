r"""Diamond-Difference (DD) cell-update strategy — geometry-polymorphic by data.

Issue #196 Phase G Step 2.5 collapses the historical 3-branch
:class:`DiamondDifference` (slab, curvilinear, cylindrical-degenerate)
into ONE polymorphic body.  Geometry is data carried by
:class:`~orpheus.geometry.reduced_operator.StreamingTerms` and
:class:`~orpheus.sn.spatial.cell_update.CellVisit`; the strategy does
NOT branch on geometry kind.

Architectural shift (Step 2.5)
==============================

Pre-Step-2.5 the slab branch was kept distinct for bit-identity with
the legacy :func:`orpheus.sn.sweep._sweep_1d_cumprod` operation order
(``a*ψ_in + 2q/denom``, then ``½(ψ_in + ψ_out)``).  Step 2.5 retires
the cumprod path entirely; slab now uses the same fold as sphere /
cylinder, and the per-cell algebra is the unified
``ψ_avg = (q + numer_upstream) / denom`` followed by the WDD spatial
closure ``ψ_out = 2·ψ_avg − ψ_in``.  These two operation orders are
algebraically identical (exact arithmetic) but differ at IEEE-754 ULP
— slab hand-calc tests re-baseline to ``np.allclose(rtol=1e-13)``
per the migration-endpoint clause documented in the pre-Step-2.5
docstring (see git history) and the project plan
``.claude/plans/issue_196_phase_g_replan.md`` §"Step 2.5".

The unified body
================

Three structural observations enable the collapse:

1. **Cell-balance algebra is one formula across geometries**
   when :class:`StreamingTerms` carries neutral curvature for slab
   (``α=0, ΔA/w=0, τ=1, A_in=A_out=1``).  The helper
   :func:`cell_balance_terms` (the Step-2.5 unified version)
   produces ``(denom, numer_upstream)`` for any geometry.
2. **The spatial closure** ``ψ^s_out = 2·ψ_avg − ψ^s_in`` is the
   same formula for slab and non-degenerate curvilinear.
   Cylindrical-degenerate has no downstream spatial face, signalled
   by ``visit.face_area_downstream == 0.0`` (geometric truth, not a
   numerical threshold).
3. **The angular closure** ``ψ^a_out = (ψ_avg − (1−τ)·ψ^a_in)/τ`` is
   the same formula for sphere and cylinder.  Slab has no angular
   redistribution, signalled by ``upstream_state.angular_upstream is
   None`` (the input direction does not exist in slab geometry).

The two ``if`` checks remaining inside :meth:`update`
(``face_area_downstream > 0.0`` for the spatial closure;
``angular_upstream is not None`` for the angular closure) are NOT
geometry dispatch — they test the **structural presence** of a
direction, not the geometry kind.

References
==========

* Hébert, A. (2009). *Applied Reactor Physics*.  Ch. 3 §3.9.4 —
  primary source for the curvilinear S\ :sub:`N` cell-balance + DD
  difference relations.
* Lewis, E. E., & Miller, W. F. (1984).  *Computational Methods of
  Neutron Transport*.  §4.5 (Morel–Montry angular closure feeding
  the curvilinear DD update); §5.3 (Diamond Difference, weighted-DD,
  Step, Linear Discontinuous; the negative-flux failure mode).
* Bailey, T. S., Morel, J. E., & Chang, J. H. (2010).  *Asymptotic
  Diffusion-Limit Accuracy of Sn Angular Differencing Schemes*.
  NSE 165(2):149-169 (LLNL preprint LLNL-JRNL-420356).  Auxiliary
  justification for the M-M weighted-diamond :math:`\tau` clamp via
  formal-:math:`\varepsilon` asymptotic-diffusion-limit analysis.

See also
========

* :class:`~orpheus.sn.spatial.cell_update.CellUpdate` — the
  Protocol this strategy satisfies.
* :func:`cell_balance_terms` — the unified algebra (Step 2.5).
* :doc:`/theory/discrete_ordinates`, "Cell update strategies (the
  strategy contract)" → "Diamond Difference" — the theory page.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

import numpy as np

from .cell_balance import cell_balance_for_streaming, cell_balance_terms
from .cell_update import (
    CellResult,
    CellUpdateBase,
    CellVisit,
    SweepCellSlice,
    UpstreamState,
)


@dataclass(frozen=True, slots=True)
class DiamondDifference(CellUpdateBase, key="diamond_difference"):
    r"""Diamond-Difference (DD) cell-update strategy — geometry-polymorphic.

    A **single** body handles slab, sphere, cylinder, and the
    cylindrical pure-azimuthal degenerate case.  Geometry flows
    through data on :class:`StreamingTerms` /
    :class:`CellVisit`; the strategy has no internal geometry
    dispatch.  See the module docstring for the three structural
    observations that enable the collapse.

    Notes
    -----
    Frozen + slotted: instances are immutable and lightweight, and a
    single :class:`DiamondDifference` instance can be reused across
    every cell update of a sweep without per-call allocation.
    """

    is_linear: ClassVar[bool] = True
    """DD is linear in ``source`` and ``upstream_state``: the cell
    average and downstream states are affine combinations of those
    inputs (with coefficients depending on geometry and cross
    section, not on the inputs themselves).  Lewis & Miller §5.3."""

    is_positivity_preserving: ClassVar[bool] = False
    r"""DD does **not** guarantee non-negative outputs from non-
    negative inputs.  In thin / large-source cells the WDD spatial
    closure :math:`\psi_{\rm out} = 2\overline{\psi} - \psi_{\rm in}`
    can produce negative outgoing face flux even when ``source`` and
    ``upstream_state`` are non-negative.  Lewis & Miller §5.3
    exhibits the canonical counter-example.
    """

    def update(
        self,
        visit: CellVisit,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> CellResult:
        r"""Compute the cell-average flux + downstream states.

        One body — no geometry dispatch.  See module docstring §
        "The unified body" for the three structural observations.
        """
        # ── Cell-balance solve: ONE formula, all geometries ─────────
        terms = cell_balance_terms(
            visit.streaming_terms,
            visit.face_area_downstream,
            total_xs,
            upstream_state,
        )
        psi_avg = (source + terms.numer_upstream) / terms.denom

        # ── Spatial closure (WDD) ───────────────────────────────────
        # Outputs ``None`` when there is no downstream spatial face on
        # this visit (cylindrical pure-azimuthal degenerate:
        # face_area_downstream == 0.0).  Slab and non-degenerate
        # curvilinear share the closure formula
        # ``ψ^s_out = 2·ψ_avg − ψ^s_in`` exactly.
        psi_spat_out: np.ndarray | None = None
        if visit.face_area_downstream > 0.0:
            psi_spat_out = 2.0 * psi_avg - upstream_state.spatial_upstream

        # ── Angular closure (Morel-Montry) ──────────────────────────
        # Outputs ``None`` when this geometry has no angular state to
        # propagate (slab: upstream_state.angular_upstream is None).
        # Sphere / cylinder share the closure formula
        # ``ψ^a_out = (ψ_avg − (1−τ)·ψ^a_in)/τ`` exactly.
        psi_angle_out: np.ndarray | None = None
        if upstream_state.angular_upstream is not None:
            tau = visit.streaming_terms.tau_mm
            psi_angle_out = (
                psi_avg - (1.0 - tau) * upstream_state.angular_upstream
            ) / tau

        return CellResult(
            cell_average_flux=psi_avg,
            outgoing_spatial_flux=psi_spat_out,
            outgoing_angular_state=psi_angle_out,
        )

    # ── Apply-direction residual (Issue #196 Phase G Step 1 replan) ──

    def residual(
        self,
        cell_avg: np.ndarray,
        visit: CellVisit,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> np.ndarray:
        r"""Per-cell operator residual :math:`L_{\rm cell}\,\bar\psi - q`.

        Companion to :meth:`update` — the **apply direction** of the
        same per-cell linear system.  At the converged cell average
        (i.e. when ``cell_avg == update(...).cell_average_flux``), the
        residual is zero to floating-point rounding.  See
        :meth:`CellUpdate.residual` for the full contract.

        Round-trip with :meth:`update`
        ------------------------------

        Issue #197 PR-TYPED-6: this method delegates to
        :func:`cell_balance_for_streaming` (the vectorized helper the
        SN matvec also consumes) at ``n_mask=1`` — Pattern 2, ONE
        algebra source, two consumers.  :meth:`update` still routes
        through :func:`cell_balance_terms` (the scalar form for the
        solve direction).  Both helpers compute the same algebraic
        intermediates; ``cell_balance_for_streaming`` exposes the
        vectorized broadcast-friendly shape, and at ``n_mask=1`` the
        per-ordinate result matches the scalar form bit-for-bit.

        At the converged ``cell_avg = update(...).cell_average_flux``,
        the residual is ``denom · cell_avg − (source +
        numer_upstream)`` = zero by the cell-balance equation that
        ``update`` solved.
        """
        st = visit.streaming_terms
        # Convert per-cell scalar StreamingTerms primitives to the
        # ``(n_mask=1,)`` arrays the vectorized helper consumes.
        # Pattern 2 — single algebra source via cell_balance_for_streaming.
        abs_mu_arr = np.array([st.abs_mu], dtype=float)
        A_down_arr = np.array([visit.face_area_downstream], dtype=float)
        A_total_arr = np.array(
            [st.face_area_inner + st.face_area_outer], dtype=float,
        )
        dA_w_arr = np.array([st.delta_A_over_w], dtype=float)
        # M-M closure constants (same algebra as cell_balance_terms).
        tau = st.tau_mm
        c_out_scalar = st.alpha_out / tau
        c_in_scalar = (1.0 - tau) / tau * st.alpha_out + st.alpha_in
        c_in_arr = np.array([c_in_scalar], dtype=float)
        c_out_arr = np.array([c_out_scalar], dtype=float)

        psi_face_in_mask = upstream_state.spatial_upstream[:, None]  # (ng, 1)
        psi_ang = upstream_state.angular_upstream
        psi_ang_mask = None if psi_ang is None else psi_ang[:, None]  # (ng, 1)

        denom, numer_upstream = cell_balance_for_streaming(
            abs_mu=abs_mu_arr,
            A_downstream=A_down_arr,
            A_total=A_total_arr,
            dA_w=dA_w_arr,
            c_in=c_in_arr,
            c_out=c_out_arr,
            total_xs=total_xs,
            volume=st.volume,
            psi_face_in=psi_face_in_mask,
            psi_angular_upstream=psi_ang_mask,
        )
        # Both arrays are (ng, 1); collapse to (ng,) to match the
        # scalar-form residual contract (Issue #196 Phase G Step 1).
        return denom[:, 0] * cell_avg - (source + numer_upstream[:, 0])

    # ── 2-D Cartesian batched update (Wave 2 / C2.2) ───────────────

    def update_batch(self, slice_args: SweepCellSlice) -> np.ndarray:
        r"""Vectorised DD update for one anti-diagonal level.

        Reproduces the inlined wavefront DD math at
        ``orpheus.sn.sweep._sweep_2d_wavefront`` lines 847-871
        bit-for-bit, with the ordinate axis (``N_oct``) and the
        anti-diagonal axis (``n_diag``) folded into a single batched
        call. The math is the **balance form** of WDD on a 2-D
        Cartesian cell:

        .. math::

           \overline{\psi}_{i,j}
           \;=\; \frac{Q_{i,j}
                       + s_{x,i}\,\psi^{\rm in}_{x,i,j}
                       + s_{y,j}\,\psi^{\rm in}_{y,i,j}}
                      {\Sigma_{t,i,j} + s_{x,i} + s_{y,j}},

           \qquad
           s_{x,i} = 2|\mu_x|/\Delta x_i,
           \quad s_{y,j} = 2|\mu_y|/\Delta y_j,

        with the spatial closure
        :math:`\psi^{\rm out}_x = 2\overline{\psi} - \psi^{\rm in}_x`
        (and analogously on :math:`y`). The closure values are
        scattered back into :attr:`SweepCellSlice.psi_x` /
        :attr:`SweepCellSlice.psi_y` at the outgoing-face indices in
        place — those buffers are persistent across levels.

        Bit-identity contract
        ---------------------

        The operation order of ``denom = sig_t + sx + sy``,
        ``psi_avg = (Q + sx*psi_in_x + sy*psi_in_y) / denom``, and
        ``psi_out = 2*psi_avg - psi_in`` matches the legacy inlined
        sweep (sweep.py:847-871) exactly. Per ``vv-principles``
        Bit-identity vs principled-equivalence, algebraically-
        equivalent rearrangements break the 1-ULP regression contract
        — do NOT refactor for "clarity".
        """
        s = slice_args
        ii, jj = s.ii, s.jj

        # Issue #196 PR-INDEX-5: every input is principled.
        #   sig_t       : (ng, nx, ny)
        #   Q           : (N_oct or 1, ng, nx, ny)
        #   psi_x       : (N_oct, ng, nx+1, ny)
        #   psi_y       : (N_oct, ng, nx, ny+1)
        # ``psi_avg`` (return) is principled ``(N_oct, ng, n_diag)``.

        # Gather incoming face fluxes.  Advanced indices ``jj`` (psi_in_x)
        # and ``ii`` (psi_in_y) are at the trailing position, contiguous
        # at the end of the slice expression — numpy keeps the index
        # order ``(N_oct, ng, n_diag)``.
        psi_in_x = s.psi_x[:, :, s.face_in_x_idx, jj]    # (N_oct, ng, n_diag)
        psi_in_y = s.psi_y[:, :, ii, s.face_in_y_idx]    # (N_oct, ng, n_diag)

        # Per-octant per-cell streaming coefficients, broadcast to
        # (N_oct, 1, n_diag) so they multiply across the group axis.
        sx = s.str_x[:, ii][:, None, :]                  # (N_oct, 1, n_diag)
        sy = s.str_y[:, jj][:, None, :]                  # (N_oct, 1, n_diag)

        # Total-xs slice on this level — principled ``(ng, n_diag)``.
        sigt_cells = s.sig_t[:, ii, jj]                  # (ng, n_diag)

        # denom builds as sig_t + sx + sy (operation order preserved
        # per Pattern 7 / vv-principles bit-identity discipline).
        denom = sigt_cells + sx + sy                     # (N_oct, ng, n_diag)

        # Q principled ``(N_oct or 1, ng, nx, ny)``; ``[:, :, ii, jj]``
        # gives ``(N_oct or 1, ng, n_diag)`` directly — no transpose.
        Q_cells = s.Q[:, :, ii, jj]                       # (N_oct or 1, ng, n_diag)

        psi_avg = (
            Q_cells
            + sx * psi_in_x
            + sy * psi_in_y
        ) / denom                                         # (N_oct, ng, n_diag)

        # Spatial closure — scatter outgoing face fluxes back into the
        # persistent buffers (principled layout).
        s.psi_x[:, :, s.face_out_x_idx, jj] = 2.0 * psi_avg - psi_in_x
        s.psi_y[:, :, ii, s.face_out_y_idx] = 2.0 * psi_avg - psi_in_y

        return psi_avg


__all__ = ["DiamondDifference"]
