r"""Diamond-Difference (DD) cell-update strategy — geometry-polymorphic by data.

Issue #196 Phase G Step 2.5 collapses the historical 3-branch
:class:`DiamondDifference` (slab, curvilinear, cylindrical-degenerate)
into ONE polymorphic body.  Geometry is data carried by
:class:`~orpheus.geometry.reduced_operator.StreamingTerms` and
:class:`~orpheus.sn.spatial.scheme.CellVisit`; the strategy does
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

* :class:`~orpheus.sn.spatial.scheme.DiscretizationScheme` — the
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
from .scheme import (
    CellResult,
    DiscretizationSchemeBase,
    CellVisit,
    UpstreamState,
)

#: Diamond-Difference cell-average blend weight ``w = ½`` (the symmetric
#: diamond mean ``ψ̄ = ½(ψ_in + ψ_out)``).  This single constant is what makes
#: a scheme "Diamond Difference": the streaming diagonal carries ``2 = 1/w_DD``
#: and every reconstruction (cell-average, outgoing face, source emission)
#: rides this ``w``.  Named so the five sites that historically hard-stamped
#: the literal ``0.5`` reference ONE source of truth (``_DD_W is exactly
#: 0.5`` — referencing it is byte-identical to the literal).
_DD_W: float = 0.5


@dataclass(frozen=True, slots=True)
class DiamondDifference(DiscretizationSchemeBase, key="diamond_difference"):
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
    so the DAG-free scan schedules (``CumprodScan``, ``ScanMarch``) consume DD
    via its per-cell coefficient triple :meth:`affine_scan_coefficients`
    ``(a, inverse_denom, w)`` and the generic base reconstruction staticmethods
    (:meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.source_emission` /
    :meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.cell_average` /
    :meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.outgoing_face_from_average`;
    #158 the coefficient model);
    DD's blend weight is ``w = ½`` (the symmetric diamond mean).
    (Linear-Discontinuous couples two face moments, but its slope is eliminated
    by the Schur complement — so LD is ALSO affine-scannable, with
    ``w = 1/(1+k)``; #158 Increment B.)"""

    transverse_coupling_is_facewise: ClassVar[bool] = True
    r"""DD's multi-D transverse coupling is **facewise** (separable): the
    coupling from a non-swept axis :math:`y` enters the x-recurrence as the
    single 0th-order face value ``s_y · ψ_{y,in}`` folded into the scan's
    affine source (the explicit per-axis left-fold in :meth:`update`), so the
    d-D DD closure factors into independent per-axis 1-D scans chained by
    scalar face traces.  DD is therefore admitted to the :math:`d \ge 2`
    scan-march (``ScanMarch``); see
    :attr:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.transverse_coupling_is_facewise`
    for the trait contract and the contrast with the (single-axis)
    ``is_affine_scannable``.  Linear-Discontinuous, whose multi-D coupling is a
    1st-order slope moment (bilinear), leaves this trait at its ``False``
    default and rides the DAG wavefront instead (#240 D5b / #38)."""

    diffusion_limit_consistent: ClassVar[bool] = True
    r"""DD's thick-diffusion limit IS a consistent diffusion discretization for
    the leading-order scalar flux (Larsen–Morel–Miller 1987 Eq. (4.24): the
    cell-average flux limits to :math:`\tfrac12(\Phi_{j+1/2}+\Phi_{j-1/2})` where
    :math:`\Phi` satisfies the edge-differenced diffusion Eq. (4.22); the
    intermediate limit Eqs. (4.33)/(4.34) hold too — LMM-1987 Table I "Diamond"
    row).  The only DD "maybe" is the thick-regime cell-EDGE flux under an
    anisotropic incident boundary (the cell-to-cell edge oscillation, Eq. 4.23),
    which is an edge-flux artefact, NOT a failure of the leading-order
    scalar-flux limit.  ⚠ This is the SPATIAL axis: DD-in-ANGLE's first-order
    :math:`\beta`-failure (the curvilinear flux dip, Bailey–Morel–Chang 2010) is
    a DISTINCT, angular result — do NOT read it as a spatial-DD deficiency.  The
    angular condition lives on the pole-angular closure; the PAIR validity is
    :func:`~orpheus.sn.spatial.pairing.pair_diffusion_limit_consistent`."""

    supports_curvilinear: ClassVar[bool] = True
    r"""DD has a curvilinear cell closure: :meth:`update` runs the Morel–Montry
    angular redistribution (the ``angular_upstream is not None`` branch) for
    sphere/cylinder, and DD rides ``CumprodScan`` on every 1-D geometry.  So a
    curvilinear mesh may select a DD scheme (the default for sphere/cylinder)."""

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
        # Issue #236 Phase 2 B3: the Morel--Montry constants c_in / c_out
        # are angular-closure-owned and arrive as DATA on the visit (stamped
        # by SNMesh._make_cell_visit from the closure's c_{in,out}_per_ordinate
        # accessors); cell_balance_terms consumes them instead of rebuilding
        # them from st.alpha_* / st.tau_mm.
        terms = cell_balance_terms(
            visit.streaming_terms,
            visit.face_area_downstream,
            total_xs,
            upstream_state,
            c_in=visit.c_in,
            c_out=visit.c_out,
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
            # DD reconstruction ``2ψ̄ − ψ_in`` is the ``w=½`` case of the
            # generic affine outflow ``(ψ̄ − (1−w)ψ_in)/w`` (byte-identical:
            # ``÷0.5`` is exact ``×2``).
            psi_spat_out = self.outgoing_face_from_average(
                psi_avg, upstream_state.spatial_upstream, _DD_W,
            )

        # ── Angular closure (Morel-Montry) ──────────────────────────
        # Outputs ``None`` when this geometry has no angular state to
        # propagate (slab: upstream_state.angular_upstream is None).
        # Sphere / cylinder share the closure formula
        # ``ψ^a_out = (ψ_avg − (1−τ)·ψ^a_in)/τ`` exactly.
        #
        # Issue #236 Phase 2 B3: τ is the angular-closure-owned weight,
        # sourced off the visit (CellVisit.tau, stamped by
        # SNMesh._make_cell_visit from the closure's tau_per_ordinate) —
        # matching the c_in / c_out provenance above.  DD no longer reads
        # the geometry-owned streaming_terms.tau_mm (which Step C retires);
        # the closure's τ is 0-ULP equal to it (Leg-1 producer-equivalence
        # gate), so this recurrence is bit-identical.
        psi_angle_out: np.ndarray | None = None
        if upstream_state.angular_upstream is not None:
            tau = visit.tau
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
        :meth:`DiscretizationScheme.residual` for the full contract.

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
        # Issue #236 Phase 2 B2: the weighted-diamond constants
        # ``c_in`` / ``c_out`` are angular-closure-owned and arrive as
        # DATA on the :class:`CellVisit` (sourced by
        # :meth:`SNMesh._make_cell_visit` from the closure's
        # ``c_{in,out}_per_ordinate`` accessors).  DD no longer rebuilds
        # them from ``st.alpha_*`` / ``st.tau_mm`` — that inline formula
        # (the former P2 duplication site) is retired.  DD stays
        # geometry-blind AND closure-blind: it consumes the angular
        # constants without seeing the closure object, preserving the
        # spatial ⊗ angular separation.  The (ΔA/w)-scaling that follows
        # is the geometry-owned redistribution factor; only the SOURCE
        # of ``c`` moved, the assembly is byte-unchanged.
        dA_w_scalar = st.delta_A_over_w
        c_out_scalar = visit.c_out
        c_in_scalar = visit.c_in
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

    # ── Shared DD Cartesian sub-primitives (single source — #240 D5a) ──
    #
    # The streaming-diagonal fold ``S = Σ_t + Σ 2 g_a`` and the ``w=½``
    # reflection coefficients are the ONE place the diamond ``2 = 1/w_DD``
    # enters the Cartesian producers.  Before #240 D5a the ``2 g_a`` fold was
    # hand-transcribed in THREE producers (``cell_kernel_batch`` /
    # ``residual_kernel_batch`` / ``cartesian_scan_coefficients``) and the
    # ``w=½`` recurrence in a fourth; these two helpers collapse them so the
    # diamond constant has exactly one home (Cardinal Rule 2 / Pattern 2).

    @staticmethod
    def _cartesian_streaming_diagonal(
        reaction_xs: np.ndarray,
        s_axes: tuple[np.ndarray, ...],
    ) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
        r"""DD's ×V streaming diagonal ``S = Σ_t + Σ_a 2 g_a`` + per-axis couplings.

        Given the local reaction-rate ``reaction_xs`` and the RAW down-face
        streaming ``g_a = |μ_a|/Δ_a`` per spatial axis (``s_axes``), returns
        ``(denom, couplings)`` where ``couplings[a] = 2 g_a`` is DD's diamond-
        scaled per-axis face term and ``denom = Σ_t + Σ_a 2 g_a`` is the cell-
        balance diagonal.  Each caller reuses ``couplings[a]`` on the upstream-
        numerator term (``Σ_a couplings[a]·ψ_in_a``) so the diamond ``2`` is
        applied ONCE per axis, never recomputed.

        Operation-order discipline (the byte-identity pin)
        --------------------------------------------------

        The diagonal is an **explicit left fold** ``((Σ_t + 2 g_0) + 2 g_1) +
        …`` over ``s_axes`` in axis order (NOT ``sum()``) — bit-identical to the
        legacy ``sigt + sx + sy`` (with ``sx = 2|μ_x|/Δx``) accumulation, because
        ``2·g_a`` equals the former pre-scaled ``2|μ_a|/Δa`` exactly (multiply by
        2 is an IEEE-754-exact power-of-2 scaling that commutes with rounding).
        The scan producer passes ``s_axes = (s_scan, *s_transverse)`` so its
        scan-first order coincides with the batch kernels' axis order — the same
        bytes from all three callers.  Per ``vv-principles`` §"Bit-identity vs
        principled-equivalence", do NOT switch to ``sum()`` or regroup.

        Scope: CARTESIAN only, by design.  The curvilinear DD diagonal
        (:meth:`affine_scan_coefficients`) is a different mathematical object —
        face-area streaming + the Morel–Montry curvature redistribution, with no
        Cartesian analogue — so it is NOT folded in here; see its "Why a separate
        diagonal" note for the structural split and the deferred unification to
        the diffusion scheme's generic advection–reaction diagonal (#242).
        """
        couplings = tuple(2.0 * s_a for s_a in s_axes)  # each: DD's 2 = 1/w_DD
        denom = reaction_xs
        for c_a in couplings:
            denom = denom + c_a                          # explicit left fold
        return denom, couplings

    @staticmethod
    def _reflection_coeffs(
        psi_bar: np.ndarray, w: "float | np.ndarray",
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""Apply-direction reflection scan coefficients ``(α = −(1−w)/w, β = ψ̄/w)``.

        The recurrence form of :meth:`outgoing_face_from_average`: with a KNOWN
        cell average ``psi_bar``, the downstream face is
        :math:`ψ_{\rm out} = (ψ̄ − (1−w)ψ_{\rm in})/w = β + α·ψ_{\rm in}`.  Pure
        in ``(ψ̄, w)``, no instance state — the shared arithmetic behind DD's
        :meth:`reflect_scan_coefficients` (``w = _DD_W`` gives ``α = −1``,
        ``β = 2ψ̄``) and the future Step twin.  Byte-identical to the inlined
        ``α = −1 / β = 2ψ̄`` at ``w = ½`` (``-(1-0.5)/0.5 == -1.0`` and
        ``ψ̄/0.5 == 2ψ̄`` exactly).
        """
        alpha = np.full_like(psi_bar, -(1.0 - w) / w)
        beta = psi_bar / w
        return alpha, beta

    # ── 2-D Cartesian batched CELL KERNEL (storage-free; Pattern 2) ──

    def cell_kernel_batch(
        self,
        *,
        psi_in: tuple[np.ndarray, ...],   # d arrays, each (N_oct, ng, n_diag) — incoming face flux per axis
        s_axes: tuple[np.ndarray, ...],   # d arrays, each (N_oct, 1, n_diag) — RAW down-face streaming g=|μ_a|/Δa per axis (DD applies its diamond 2 here)
        reaction_xs: np.ndarray,          # (ng, n_diag) — Σ_t on the level
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
        :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_full` (the
        VERIFICATION ORACLE policy) AND the rolling moving-frontier
        :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed`
        (the storage-B PRODUCTION policy).  The two differ ONLY in how the
        incoming faces are gathered and the outgoing faces scattered — the
        cell algebra is identical, so the window walk and the oracle cannot
        drift.

        Axis convention
        ---------------

        Element ``a`` of every per-axis tuple (``psi_in``, ``s_axes``, and the
        returned ``psi_out``) is spatial axis ``a`` — the SAME axis order as
        :attr:`~orpheus.sn.loss_representation.sweep_graph.OctantLabel.signs` and
        :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.from_cartesian`'s
        ``shape``. The tuples are positional-by-axis; the caller MUST build
        and unpack them in axis order (``[0]`` = x, ``[1]`` = y, ``[2]`` = z).

        Operation-order discipline (and d=2 bit-identity)
        -------------------------------------------------

        The math is :math:`\psi_{\rm avg} = (Q + \sum_a 2 g_a\,\psi^{\rm in}_a) /
        (\Sigma_t + \sum_a 2 g_a)`, closure :math:`\psi^{\rm out}_a = 2\psi_{\rm
        avg} - \psi^{\rm in}_a`, where ``s_a = g_a`` is the RAW down-face
        streaming and the **scheme** applies its diamond factor ``2 = 1/w_DD``
        (#240): both the denominator term AND the upstream-numerator term gain
        the ``2`` (a denom-only ``2`` would be a non-uniform 2× bug → wrong
        ``ψ̄``). The closure ``2ψ̄ − ψ_in`` is the diamond MEAN (also a ``2``,
        but the cell-average reconstruction, not the streaming factor).
        The streaming diagonal + the per-axis couplings ``2 g_a`` come from
        :meth:`_cartesian_streaming_diagonal` (the single source of DD's ``2 g``
        fold); the upstream-numerator term reuses those couplings.  The axis
        reduction is an **explicit left fold** (NOT ``sum()``) so the
        accumulation order is ``((sigt + 2g_0) + 2g_1) + …`` /
        ``((Q + 2g_0·in_0) + 2g_1·in_1) + …`` — bit-identical to the legacy
        ``sigt + sx + sy`` (with ``sx = 2|μ_x|/Δx``) order, because ``2·g_a``
        equals the former pre-scaled ``2|μ_a|/Δa`` exactly (multiply by 2 is an
        IEEE-754-exact power-of-2 scaling that commutes with rounding). Per
        ``vv-principles`` Bit-identity vs principled-equivalence, do NOT switch
        to ``sum()`` or rearrange for "clarity".
        """
        denom, couplings = self._cartesian_streaming_diagonal(reaction_xs, s_axes)
        numer = Q_cells                                    # (N_oct or 1, ng, n_diag)
        for c_a, in_a in zip(couplings, psi_in):
            numer = numer + c_a * in_a                     # reuse 2g_a coupling
        psi_avg = numer / denom                            # (N_oct, ng, n_diag)
        # DD diamond MEAN reconstruction = the ``w=½`` generic affine outflow
        # (byte-identical: ``÷0.5`` is exact ``×2``).
        psi_out = tuple(
            self.outgoing_face_from_average(psi_avg, in_a, _DD_W) for in_a in psi_in
        )
        return psi_avg, psi_out

    def residual_kernel_batch(
        self,
        *,
        psi_bar: np.ndarray,             # (N_oct, ng, n_diag) — the probe cell-average
        psi_in: tuple[np.ndarray, ...],  # d arrays, each (N_oct, ng, n_diag)
        s_axes: tuple[np.ndarray, ...],  # d arrays, each (N_oct, 1, n_diag) — RAW down-face streaming g per axis (DD applies its diamond 2 here)
        reaction_xs: np.ndarray,         # (ng, n_diag)
        Q_cells: np.ndarray,             # (N_oct or 1, ng, n_diag)
    ) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
        r"""Pure batched DD operator residual — dimension-generic, NO storage.

        The apply-direction companion of :meth:`cell_kernel_batch` (same
        positional-by-axis tuple convention — see its "Axis convention"):
        evaluates :math:`r = (\Sigma_t + \sum_a 2 g_a)\,\psi_{\rm bar} − (Q +
        \sum_a 2 g_a\,\psi^{\rm in}_a)` at the PROBE cell-average, and
        reconstructs the outgoing faces from the probe
        (:math:`\psi^{\rm out}_a = 2\psi_{\rm bar} − \psi^{\rm in}_a`).  ``s_a =
        g_a`` is the RAW down-face streaming; the **scheme** applies its diamond
        ``2 = 1/w_DD`` to BOTH the denom and the upstream-numerator term (#240).
        Returns ``(residual, psi_out)`` with ``psi_out`` the d-tuple of outgoing
        faces. Single source of the matvec cell math — the APPLY arm of the
        ``_CellResidual`` level operation, shared by the same two storage
        walks as :meth:`cell_kernel_batch`
        (:meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_full` /
        :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed`).
        Same explicit-left-fold operation-order discipline as
        :meth:`cell_kernel_batch`, via the shared
        :meth:`_cartesian_streaming_diagonal` (``2·g_a`` equals the former
        pre-scaled ``2|μ_a|/Δa`` bit-for-bit; d=2 bit-identical to the legacy
        fold order).
        """
        denom, couplings = self._cartesian_streaming_diagonal(reaction_xs, s_axes)
        numer = Q_cells                                    # (N_oct or 1, ng, n_diag)
        for c_a, in_a in zip(couplings, psi_in):
            numer = numer + c_a * in_a                     # reuse 2g_a coupling
        residual = denom * psi_bar - numer                 # (N_oct, ng, n_diag)
        # DD diamond MEAN reconstruction = the ``w=½`` generic affine outflow
        # (byte-identical: ``÷0.5`` is exact ``×2``).
        psi_out = tuple(
            self.outgoing_face_from_average(psi_bar, in_a, _DD_W) for in_a in psi_in
        )
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
        reaction_xs: np.ndarray,  # (N, ng, nx) Σ_t in the SAME cell ordering as the geometry arrays
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        r"""DD's :math:`\Sigma_t`-epoch scan coefficients — ``(a_attenuation, inverse_denom, face_blend_weight)``.

        The third coefficient ``face_blend_weight`` is DD's blend weight
        ``w = ½`` (broadcast; the symmetric diamond mean
        :math:`\bar\psi = \tfrac12(\psi_{\rm in}+\psi_{\rm out})`) — #158 the
        coefficient model.  The generic base reconstruction staticmethods
        (:meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.cell_average` /
        :meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.outgoing_face_from_average`
        / :meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.source_emission`)
        consume it; DD carries no cell-average / outgoing-face / source-
        emission method of its own.

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
        ``reaction_xs`` is therefore the chain-gathered ``(N, ng, nx)`` tensor,
        not the natural ``(ng, nx)`` field.

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

        Why a SEPARATE diagonal from the Cartesian helper (#242)
        -----------------------------------------------------------------------

        Both compute "the cell-balance diagonal ``S``" — a shared CONCEPT — but
        the REALIZATIONS are different mathematical objects, so they are NOT
        folded into one helper (Cardinal Rule 2 reconsidered; the Rule-of-Three
        tripwire is explicitly retired for this pair):

        * Cartesian (:meth:`_cartesian_streaming_diagonal`):
          ``S = Σ_t + Σ_a 2 g_a`` with ``g_a = |μ_a|/Δ_a`` — ÷Δ raw down-face
          streaming, NO volume weighting, NO angular coupling.
        * Curvilinear (here): ``S = Σ_t·V + 2|μ|·A_down + (ΔA/w)·c_out`` — the ×V
          collision, FACE-AREA streaming, AND the Morel–Montry curvature
          redistribution ``(ΔA/w)·c_out`` (``c_out = α_out/τ``).  That
          redistribution term couples the SPATIAL diagonal to the ANGULAR
          closure and has NO Cartesian analogue; the op-order reproduces
          ``CollisionCache.from_geometry``'s factoring — a DIFFERENT bit-identity
          pin (``rtol=1e-12`` slab snapshots) than the Cartesian explicit-left-fold.

        The genuine unification is ``S = Σ_t·V + streaming_diag`` with
        ``streaming_diag`` a geometry-parameterised contribution (Cartesian
        ``Σ 2g`` / curvilinear ``2|μ|A_down + (ΔA/w)c_out`` / a diffusion
        ``∇·D∇`` term) — the generic ADVECTION–REACTION diagonal the diffusion
        scheme (#240's next model / the consistent-DSA ``A_diff``) will build and
        consume.  That is where the merge belongs: it needs the diffusion
        consumer to validate the abstraction (``coding-elegance`` Pattern 6 —
        build the primitive once two real occupants exist), and forcing the
        curvilinear producer through the Cartesian helper NOW would re-baseline
        its snapshots for no immediate gain.  Deferred to that work (#242).
        """
        # streaming + curvature (no group axis) — units: dimensionless
        streaming_face_term = 2.0 * abs_mu[:, None] * A_down              # (N, nx)
        curvature_redistribution_term = dA_w * c_out[:, None]            # (N, nx)
        geometric_streaming_term = (
            streaming_face_term + curvature_redistribution_term
        )                                                                # (N, nx)
        # collision Σ_t·V (group-resolved) — units: cm² (1/cm × cm³)
        collision_volume_term = reaction_xs * V[:, None, :]              # (N, ng, nx)
        denom = geometric_streaming_term[:, None, :] + collision_volume_term  # (N, ng, nx)
        inverse_denom = 1.0 / denom                                       # (N, ng, nx)
        # transmission multiplier a = 2|μ|·A_total / denom − 1
        a_numer = 2.0 * abs_mu[:, None] * A_total                        # (N, nx)
        a_attenuation = a_numer[:, None, :] * inverse_denom - 1.0        # (N, ng, nx)
        # Cell-average blend weight (#158 the coefficient model): DD is the
        # symmetric diamond mean ψ̄ = ½(ψ_in+ψ_out), i.e. w = ½ everywhere.  The
        # generic base reconstruction staticmethods consume this; DD carries NO
        # cell-average / outgoing-face / source-emission method of its own.
        face_blend_weight = np.full_like(a_attenuation, _DD_W)           # (N, ng, nx)
        return a_attenuation, inverse_denom, face_blend_weight

    def cartesian_scan_coefficients(
        self,
        *,
        s_scan: np.ndarray,                # (..., n_scan) RAW down-face g on the scanned axis
        s_transverse: tuple[np.ndarray, ...],  # d−1 arrays, each broadcasting over (..., n_scan) — RAW g on the known transverse axes
        reaction_xs: np.ndarray,           # (ng, n_scan) Σ_t on the row
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[np.ndarray, ...]]:
        r"""DD's row-march scan coefficients ``(a, inverse_denom, w, transverse_couplings)``.

        The 2-D/3-D scan-march analogue of :meth:`affine_scan_coefficients`
        (#239 the coefficient-model lift): given the RAW down-face streaming
        ``g = |μ|/Δ`` on the scanned axis ``s_scan`` and the known transverse
        axes ``s_transverse``, plus :math:`\Sigma_t`, return the DD
        :math:`\Sigma_t`-epoch coefficients.  DD applies its diamond
        :math:`2 = 1/w_{\rm DD}` HERE (the scheme owns the factor — the row body
        carries NO ``2.0*`` and NO hard-coded ``0.5``):

        .. math::

           \text{diag}_{\rm scan} &= 2 g_{\rm scan},
           \qquad c_\perp = 2 g_\perp, \\
           S &= \Sigma_t + \text{diag}_{\rm scan} + \sum_\perp c_\perp,
           \qquad \mathrm{inverse\_denom} = 1/S, \\
           a &= 2\,\text{diag}_{\rm scan}\,\mathrm{inverse\_denom} - 1,
           \qquad w = \tfrac12.

        The transverse couplings :math:`c_\perp = 2 g_\perp` are scheme-owned
        and feed BOTH the diagonal (above) and the affine source (the caller
        folds :math:`\sum c_\perp \psi_\perp^{\rm in}` into ``QV`` for
        :meth:`source_emission`).  The same single source of the DD cell math
        as :meth:`cell_kernel_batch` / :meth:`residual_kernel_batch` (the
        explicit ``2 g`` left fold), specialised to one scanned + ``d−1`` known
        axes — Cardinal Rule 2 / Pattern 2.

        Operation-order note
        --------------------

        ``inverse_denom = 1/S`` is the reciprocal cell-balance diagonal — the
        ×V "denom" convention the generic reconstruction staticmethods consume
        (the ``×inverse_denom`` form, NOT the legacy ``÷S`` division).  The
        ``2 g_a`` left fold over (scan, transverse...) comes from the shared
        :meth:`_cartesian_streaming_diagonal` (passing ``s_axes =
        (s_scan, *s_transverse)`` so the scan-first order coincides with
        :meth:`cell_kernel_batch`'s axis order); per ``vv-principles``
        §"Bit-identity vs principled-equivalence" the scan SOLVE rides the
        byte-identical ``×inverse_denom`` reconstruction (the DD scan snapshots
        stay strict).
        """
        # S = Σ_t + 2g_scan + Σ 2g_⊥, with couplings[0]=2g_scan, couplings[1:]=2g_⊥.
        denom, couplings = self._cartesian_streaming_diagonal(
            reaction_xs, (s_scan, *s_transverse),
        )                                                                 # (ng, n_scan)
        scan_diag, transverse_couplings = couplings[0], couplings[1:]
        inverse_denom = 1.0 / denom
        a = 2.0 * scan_diag * inverse_denom - 1.0                         # transmission
        w = np.full_like(a, _DD_W)                                        # DD diamond mean
        return a, inverse_denom, w, transverse_couplings

    def reflect_scan_coefficients(
        self, psi_bar: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""DD apply-direction reflection scan coefficients ``(α = −1, β = 2ψ̄)``.

        The recurrence form of DD's diamond reconstruction
        :math:`\psi_{\rm out} = 2\bar\psi - \psi_{\rm in}` =
        :meth:`outgoing_face_from_average` at ``w = ½``: ``α = −(1−w)/w = −1``
        (a pure reflection, ``|α| = 1`` — no scan underflow) and ``β = ψ̄/w =
        2ψ̄``.  Delegates to the shared :meth:`_reflection_coeffs` at
        ``w = _DD_W`` — the ``w``-generic arithmetic lives once (the diamond
        ``½`` lives in the scheme; the row body carries no inline ``2.0*``), and
        the Step twin will inherit it free.  Byte-identical to the legacy inline
        ``α = −1`` / ``β = 2ψ̄`` (``-(1-0.5)/0.5 == -1.0`` and ``ψ̄/0.5 == 2ψ̄``).
        """
        return self._reflection_coeffs(psi_bar, _DD_W)

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
    #      ``residual_kernel_batch``) — the DAG wavefront family (and the matvec
    #      apply twin of the scan, since the apply direction has a concrete ψ̄);
    #   3. the scan-family coefficients ``affine_scan_coefficients`` →
    #      ``(a, inverse_denom, w)`` — the DAG-free scan SOLVE family
    #      (CumprodScan / ScanMarch), consumed by the generic base
    #      reconstruction staticmethods (``source_emission`` / ``cell_average``
    #      / ``outgoing_face_from_average`` on
    #      :class:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase`; #158 the
    #      coefficient model; the per-scheme ``cell_average_from_faces`` /
    #      ``outgoing_face_from_average`` closure methods are RETIRED — the
    #      operations are now generic in ``w``).
    # This is the ONLY direction-aware math in the SN stack, and the override
    # point for future closure strategies (Step / LD / EC supply group 2 and,
    # if affine-scannable, the group-3 coefficients; storage is handled once,
    # above them).



__all__ = ["DiamondDifference"]
