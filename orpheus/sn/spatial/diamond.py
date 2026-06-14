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
the legacy ``_sweep_1d_cumprod`` (the dissolved ``sweep.py``) operation order
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

    is_affine_scannable: ClassVar[bool] = True
    r"""DD admits the closed-form affine recurrence
    :math:`\psi_{\rm out} = a\,\psi_{\rm in} + b` (Blelloch §1.5): the
    cell-average is an affine function of the SINGLE upstream face flux,
    so the DAG-free scan schedules (``CumprodScan``, ``ScanMarch``) can
    consume DD via :meth:`affine_scan_coefficients` +
    :meth:`cell_average_from_faces` + :meth:`outgoing_face_from_average`.
    (A Linear-Discontinuous closure couples two face moments and is
    therefore NOT affine-scannable — it would set this ``False`` and run
    on the DAG wavefront schedule.)"""

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

        psi_face_in_mask = upstream_state.spatial_upstream[:, None]  # (ng, 1)

        # PR-TYPED-6.5 Phase 2.11: the M-M-specific arguments
        # (``dA_w``, ``c_in``, ``c_out``, ``psi_angular_upstream``)
        # are gone from ``cell_balance_for_streaming``'s interface.
        # The closure now owns the M-M algebra and produces the
        # ``(angular_denom_term, angular_numer_upstream)`` contributions.
        #
        # DD.residual reads the M-M coefficients directly from
        # ``StreamingTerms`` (carried on the CellVisit), computes the
        # contributions inline, and passes them to the helper.  This
        # in-line computation is a Phase 6 cleanup target — the sweep
        # route will route through ``closure.cell_contribution(...)``
        # like the matvec does (Pattern 2: ONE strategy contract for
        # the angular contribution, two consumers).  For Phase 2.11
        # the inline path stays so we don't risk the sweep before the
        # matvec is verified.
        dA_w_scalar = st.delta_A_over_w
        tau = st.tau_mm
        c_out_scalar = st.alpha_out / tau
        c_in_scalar = (1.0 - tau) / tau * st.alpha_out + st.alpha_in
        angular_denom_term = np.array(
            [dA_w_scalar * c_out_scalar], dtype=float,
        )                                                # (1,)
        psi_ang = upstream_state.angular_upstream
        if psi_ang is None:
            angular_numer_upstream = np.zeros(
                (total_xs.size, 1), dtype=float,
            )                                            # (ng, 1)
        else:
            angular_numer_upstream = (
                dA_w_scalar * c_in_scalar * psi_ang[:, None]
            )                                            # (ng, 1)

        denom, numer_upstream = cell_balance_for_streaming(
            abs_mu=abs_mu_arr,
            A_downstream=A_down_arr,
            A_total=A_total_arr,
            total_xs=total_xs,
            volume=st.volume,
            psi_face_in=psi_face_in_mask,
            angular_denom_term=angular_denom_term,
            angular_numer_upstream=angular_numer_upstream,
        )
        # Both arrays are (ng, 1); collapse to (ng,) to match the
        # scalar-form residual contract (Issue #196 Phase G Step 1).
        return denom[:, 0] * cell_avg - (source + numer_upstream[:, 0])

    # ── 2-D Cartesian batched CELL KERNEL (storage-free; Pattern 2) ──

    def cell_kernel_batch(
        self,
        *,
        psi_in: tuple[np.ndarray, ...],   # d arrays, each (N_oct, ng, n_diag) — incoming face flux per axis
        s_axes: tuple[np.ndarray, ...],   # d arrays, each (N_oct, 1, n_diag) — streaming coeff 2|μ_a|/Δa per axis
        sigt_cells: np.ndarray,           # (ng, n_diag) — Σ_t on the level
        Q_cells: np.ndarray,              # (N_oct or 1, ng, n_diag) — weight-normalised
    ) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
        r"""Pure batched WDD cell update — dimension-generic, NO storage access.

        Given the per-axis incoming face fluxes ``psi_in`` + streaming
        coefficients ``s_axes`` (one entry per spatial axis, ``d = 1, 2, 3``) +
        source on one anti-hyperplane level, returns ``(psi_avg, psi_out)``
        where ``psi_out`` is the d-tuple of outgoing face fluxes (one per
        axis). This is the **single source of the DD cell math** (Cardinal
        Rule 2 / Pattern 2): the SOLVE arm of the ``_CellSolve`` level
        operation, consumed by BOTH storage walks — the full-cochain
        :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.walk_full` (the
        VERIFICATION ORACLE policy) AND the rolling moving-frontier
        :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.walk_windowed`
        (the storage-B PRODUCTION policy).  The two differ ONLY in how the
        incoming faces are gathered and the outgoing faces scattered — the
        cell algebra is identical, so the window walk and the oracle cannot
        drift.

        Axis convention
        ---------------

        Element ``a`` of every per-axis tuple (``psi_in``, ``s_axes``, and the
        returned ``psi_out``) is spatial axis ``a`` — the SAME axis order as
        :attr:`~orpheus.sn.sweep_graph.OctantLabel.signs` and
        :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.from_cartesian`'s
        ``shape``. The tuples are positional-by-axis; the caller MUST build
        and unpack them in axis order (``[0]`` = x, ``[1]`` = y, ``[2]`` = z).

        Operation-order discipline (and d=2 bit-identity)
        -------------------------------------------------

        The math is :math:`\psi_{\rm avg} = (Q + \sum_a s_a\,\psi^{\rm in}_a) /
        (\Sigma_t + \sum_a s_a)`, closure :math:`\psi^{\rm out}_a = 2\psi_{\rm
        avg} - \psi^{\rm in}_a`. The axis reduction is an **explicit left
        fold** (NOT ``sum()``) so the accumulation order is
        ``((sigt + s_0) + s_1) + …`` / ``((Q + s_0·in_0) + s_1·in_1) + …`` —
        which at ``d = 2`` reproduces the legacy ``sigt + sx + sy`` /
        ``Q + sx·in_x + sy·in_y`` order **bit-for-bit** (IEEE-754 addition is
        non-associative; the fold order is load-bearing). Per ``vv-principles``
        Bit-identity vs principled-equivalence, do NOT switch to ``sum()`` or
        rearrange for "clarity" — it breaks the 1-ULP regression contract.
        """
        denom = sigt_cells                                 # (ng, n_diag)
        numer = Q_cells                                    # (N_oct or 1, ng, n_diag)
        for s_a, in_a in zip(s_axes, psi_in):
            denom = denom + s_a                            # left fold → (N_oct, ng, n_diag)
            numer = numer + s_a * in_a
        psi_avg = numer / denom                            # (N_oct, ng, n_diag)
        psi_out = tuple(2.0 * psi_avg - in_a for in_a in psi_in)
        return psi_avg, psi_out

    def residual_kernel_batch(
        self,
        *,
        psi_bar: np.ndarray,             # (N_oct, ng, n_diag) — the probe cell-average
        psi_in: tuple[np.ndarray, ...],  # d arrays, each (N_oct, ng, n_diag)
        s_axes: tuple[np.ndarray, ...],  # d arrays, each (N_oct, 1, n_diag) — streaming coeff per axis
        sigt_cells: np.ndarray,          # (ng, n_diag)
        Q_cells: np.ndarray,             # (N_oct or 1, ng, n_diag)
    ) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
        r"""Pure batched DD operator residual — dimension-generic, NO storage.

        The apply-direction companion of :meth:`cell_kernel_batch` (same
        positional-by-axis tuple convention — see its "Axis convention"):
        evaluates :math:`r = (\Sigma_t + \sum_a s_a)\,\psi_{\rm bar} − (Q +
        \sum_a s_a\,\psi^{\rm in}_a)` at the PROBE cell-average, and
        reconstructs the outgoing faces from the probe
        (:math:`\psi^{\rm out}_a = 2\psi_{\rm bar} − \psi^{\rm in}_a`). Returns
        ``(residual, psi_out)`` with ``psi_out`` the d-tuple of outgoing faces.
        Single source of the matvec cell math — the APPLY arm of the
        ``_CellResidual`` level operation, shared by the same two storage
        walks as :meth:`cell_kernel_batch`
        (:meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.walk_full` /
        :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.walk_windowed`).
        Same explicit-left-fold operation-order discipline as
        :meth:`cell_kernel_batch` (d=2 bit-identical to the legacy
        ``sigt + sx + sy`` / ``Q + sx·in_x + sy·in_y`` order).
        """
        denom = sigt_cells                                 # (ng, n_diag)
        numer = Q_cells                                    # (N_oct or 1, ng, n_diag)
        for s_a, in_a in zip(s_axes, psi_in):
            denom = denom + s_a                            # left fold
            numer = numer + s_a * in_a
        residual = denom * psi_bar - numer                 # (N_oct, ng, n_diag)
        psi_out = tuple(2.0 * psi_bar - in_a for in_a in psi_in)
        return residual, psi_out

    # ── Scan-family capability (Issue #236 §2 — the DAG-free schedules) ──

    def affine_scan_coefficients(
        self,
        *,
        abs_mu: np.ndarray,    # (N,)        |μ_n|
        A_down: np.ndarray,    # (N, nx)     downstream face area (sweep-resolved)
        A_total: np.ndarray,   # (N, nx)     A_inner + A_outer
        dA_w: np.ndarray,      # (N, nx)     ΔA / w_n  (curvature redistribution)
        c_out: np.ndarray,     # (N,)        α_out / τ  (M-M outgoing closure const)
        V: np.ndarray,         # (N, nx)     cell volume per ordinate
        sig_t: np.ndarray,     # (N, ng, nx) Σ_t in the SAME cell ordering as the geometry arrays
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""DD's :math:`\Sigma_t`-epoch scan coefficients — ``(a_attenuation, inverse_denom)``.

        The single source of the DD affine-recurrence coefficients
        (Cardinal Rule 2 / Pattern 2): :class:`CollisionCache.from_geometry`
        consumes this to populate the L15 σ_t-stratum cache, and through it
        the ``CumprodScan`` / ``ScanMarch`` sweep bodies.  Before Issue #236
        §2 the same three numpy ops lived inlined in ``from_geometry``; they
        now live HERE so the cell-update scheme owns its own closure math
        and the scan schedules consume whichever scheme ``SNMesh`` selected.

        The math (Lewis & Miller §5.3; Hébert §3.9.4 for the curvilinear
        curvature/M-M terms):

        .. math::

           \mathrm{denom}[n,g,i] &= \underbrace{2|\mu_n|\,A_{\rm down}[n,i]}
                                      _{\text{streaming face}}
              + \underbrace{(\Delta A[n,i]/w_n)\,c_{\rm out}[n]}
                  _{\text{curvature redistribution}}
              + \underbrace{\Sigma_t[g,i]\,V[n,i]}_{\text{collision}}, \\
           a[n,g,i] &= \frac{2|\mu_n|\,A_{\rm total}[i]}{\mathrm{denom}[n,g,i]} - 1,
           \qquad
           \mathrm{inverse\_denom}[n,g,i] = \frac{1}{\mathrm{denom}[n,g,i]}.

        For slab the curvature fields are neutral (``dA_w = 0``,
        ``c_out = 0``, ``A_down = 1``, ``A_total = 2``) and the denominator
        collapses to the slab form :math:`2|\mu_n| + \Sigma_t V`.

        What this method does NOT compute (by design):

        * the **source-dependent** emission ``b = 2·(QV + angular_contrib)·
          inverse_denom`` — re-derived every sweep from a NEW source, so it
          stays in the sweep body, not in this σ_t-epoch cache surface;
        * the **order-dependent** ``cumprod_a = cumprod(a, cell_axis)`` — a
          prefix product along the chain (a *scan-schedule* transform), so it
          stays in :class:`CollisionCache`.

        Ordering & vectorisation contract
        ---------------------------------

        Pure broadcasting — exactly three numpy ops, NEVER a per-cell Python
        loop (that would re-introduce the cost the L15 cache eliminated;
        ``[[lessons-L16]]``).  The method is **elementwise** in
        :math:`(n, g, i)` and therefore ordering-agnostic: every ``(N, …, nx)``
        argument shares ONE per-ordinate cell ordering chosen by the caller,
        and the returned ``(N, ng, nx)`` arrays match it.  The scan caller
        passes *chain* order (so it can ``cumprod`` ``a`` along the cell axis);
        ``sig_t`` is therefore the chain-gathered ``(N, ng, nx)`` tensor, not
        the natural ``(ng, nx)`` field.

        Operation-order discipline (the bit-identity TRAP)
        --------------------------------------------------

        The three ops reproduce ``CollisionCache.from_geometry``'s factoring
        EXACTLY — ``V`` folded into the collision term, the
        ``streaming + curvature`` grouping, the ``a_numer · inverse_denom − 1``
        form.  IEEE-754 addition/multiplication is non-associative: this is
        the order the slab regression snapshots are pinned to (``rtol=1e-12``),
        so do NOT regroup for "clarity" (it is NOT the ``cell_kernel_batch``
        explicit-left-fold order — that is the batch capability's discipline).
        Per ``vv-principles`` Bit-identity vs principled-equivalence.
        """
        # streaming + curvature (no group axis) — units: dimensionless
        streaming_face_term = 2.0 * abs_mu[:, None] * A_down              # (N, nx)
        curvature_redistribution_term = dA_w * c_out[:, None]            # (N, nx)
        geometric_streaming_term = (
            streaming_face_term + curvature_redistribution_term
        )                                                                # (N, nx)
        # collision Σ_t·V (group-resolved) — units: cm² (1/cm × cm³)
        collision_volume_term = sig_t * V[:, None, :]                    # (N, ng, nx)
        denom = geometric_streaming_term[:, None, :] + collision_volume_term  # (N, ng, nx)
        inverse_denom = 1.0 / denom                                       # (N, ng, nx)
        # transmission multiplier a = 2|μ|·A_total / denom − 1
        a_numer = 2.0 * abs_mu[:, None] * A_total                        # (N, nx)
        a_attenuation = a_numer[:, None, :] * inverse_denom - 1.0        # (N, ng, nx)
        return a_attenuation, inverse_denom

    def cell_average_from_faces(
        self, face_in: np.ndarray, face_out: np.ndarray,
    ) -> np.ndarray:
        r"""Diamond cell-average :math:`\bar\psi = \tfrac12(\psi_{\rm in}+\psi_{\rm out})`.

        Closure solve-direction A (see :meth:`CellUpdate.cell_average_from_faces`):
        the scan recurrence produces the downstream face flux directly, and
        the cell-average is the diamond mean of the (incoming, outgoing) face
        pair.  Consumed by every scan body (slab ``_run_1d_unified``,
        curvilinear, ``ScanMarch`` x-row) — the single source replacing the
        ``0.5 * (psi_face_in + psi_face_out)`` written inline at each site.
        Inverse direction: :meth:`outgoing_face_from_average`.
        """
        return 0.5 * (face_in + face_out)

    def outgoing_face_from_average(
        self, cell_avg: np.ndarray, face_in: np.ndarray,
    ) -> np.ndarray:
        r"""Diamond reconstruction :math:`\psi_{\rm out} = 2\bar\psi - \psi_{\rm in}`.

        Closure solve-direction B (the inverse of
        :meth:`cell_average_from_faces`): given the cell-average and the
        incoming face, reconstruct the outgoing face via the WDD difference
        relation.  This is the SAME relation already spelled inline by the
        per-cell :meth:`update` spatial closure (:math:`2\psi_{\rm avg} -
        \psi^s_{\rm in}`) and the batched :meth:`cell_kernel_batch` /
        :meth:`residual_kernel_batch` per-axis reconstruction; the named
        method makes it a single source the ``ScanMarch`` transverse y-march
        consumes (Issue #236 §2 PR3).
        """
        return 2.0 * cell_avg - face_in

    # ── S6.4(e): the storage adapters RETIRED ───────────────────────────
    #
    # ``update_batch`` / ``residual_batch`` (the full-field gather → kernel →
    # scatter wrappers) and their ``_cell_face_selector`` /
    # ``_gather_cell_inputs`` / ``_scatter_outgoing_faces`` halves moved to
    # the WALK layer (``SweepDependencyGraph.walk_full`` — storage is the
    # walk's concern, not the discretization's).  This class is now pure
    # cell algebra in THREE capability groups, one per sweep-schedule family:
    #   1. the per-cell reference pair (``update`` / ``residual``) — the
    #      canonical contract every scheme MUST provide;
    #   2. the batched kernel pair (``cell_kernel_batch`` /
    #      ``residual_kernel_batch``) — the DAG wavefront family;
    #   3. the scan-family triple (``affine_scan_coefficients`` /
    #      ``cell_average_from_faces`` / ``outgoing_face_from_average``) — the
    #      DAG-free scan family (CumprodScan / ScanMarch), Issue #236 §2.
    # This is the ONLY direction-aware math in the SN stack, and the override
    # point for future closure strategies (Step / LD / EC supply groups 2 and,
    # if affine-scannable, 3; storage is handled once, above them).



__all__ = ["DiamondDifference"]
