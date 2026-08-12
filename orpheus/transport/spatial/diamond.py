r"""Diamond-Difference (DD) cell-update strategy — geometry-polymorphic by data.

A **single** polymorphic body handles slab, sphere, cylinder, and the
cylindrical pure-azimuthal degenerate case.  Geometry is data carried by
:class:`~orpheus.geometry.reduced_operator.StreamingTerms` and
:class:`~orpheus.transport.spatial.scheme.CellVisit`; the strategy does
NOT branch on geometry kind.  Three structural observations enable the
collapse (derivation:
``docs/theory/foundations/discretization.rst §discretization-space-angle``):

1. **Cell-balance algebra is one formula across geometries** when
   :class:`StreamingTerms` carries neutral curvature for slab
   (``α=0, ΔA/w=0, τ=1, A_in=A_out=1``); :func:`cell_balance_terms`
   produces ``(denom, numer_upstream)`` for any geometry.
2. **The spatial closure** ``ψ^s_out = 2·ψ_avg − ψ^s_in`` is the same
   formula for slab and non-degenerate curvilinear; cylindrical-degenerate
   has no downstream spatial face, signalled by
   ``visit.face_area_downstream == 0.0`` (geometric truth, not a numerical
   threshold).
3. **The angular closure** ``ψ^a_out = (ψ_avg − (1−τ)·ψ^a_in)/τ`` is the
   same formula for sphere and cylinder; slab has no angular redistribution,
   signalled by ``upstream_state.angular_upstream is None``.

The two ``if`` checks remaining inside :meth:`update`
(``face_area_downstream > 0.0``; ``angular_upstream is not None``) are NOT
geometry dispatch — they test the **structural presence** of a direction,
not the geometry kind.

References
==========

Sources are listed by WHAT they are the authority for.  The **spatial**
DD relations, the curvilinear cell balance and the weighted **angular**
:math:`\tau` come from three different places; conflating them is the
error class this file has already paid for (theory-page record:
``docs/theory/methods/sn/curvilinear_one_group.rst
§sn-citation-corrections``).

* Hébert, A. (2009). *Applied Reactor Physics*.  Ch. 3 **§3.9.3
  (cylinder, pp. 137-141)** and **§3.9.4 (sphere, pp. 141-144)** —
  the curvilinear S\ :sub:`N` cell-balance, the :math:`\Delta A/w`
  factor, the sweep ordering and the Carlson starting direction.
  **NOT** the source of the weighted angular :math:`\tau`: he defines
  no :math:`\tau` anywhere in chapter 3 and ships the *plain* angular
  diamond (Eqs. 3.437/3.439 sphere, 3.412/3.414 cylinder, i.e.
  :math:`\tau \equiv \tfrac12`).
* Morel, J. E., & Montry, G. R. (1984).  *Analysis and Elimination of
  the Discrete-Ordinates Flux Dip*.  Transport Theory and Statistical
  Physics 13(5):615-633, doi:10.1080/00411458408211661.  **PRIMARY**
  for the weighted angular closure this file's step 3 applies.
* Bailey, T. S., Morel, J. E., & Chang, J. H. (2010).  *The Asymptotic
  Diffusion-Limit Accuracy of Sn Angular Differencing Schemes*.
  NSE 165(2):149-169 (LLNL preprint LLNL-JRNL-420356).  **Eqs. (42)/(43)
  are the form of** :math:`\tau` **implemented here** — the barycentric
  coordinate of the ordinate between its own angular cell's two edges;
  their Eq. (41) is the first-order diffusion-limit condition
  :math:`\beta = 0`, and forcing it to zero is what DETERMINES the
  weights (:math:`\tau` is derived, not chosen).  Their Eq. 53 + §I show
  the plain diamond is diffusion-limit consistent only to LEADING order
  while the weighted one is correct through FIRST order — so BMC is not
  an "auxiliary" justification, it is the scheme.
* Lewis, E. E., & Miller, W. F. (1984).  *Computational Methods of
  Neutron Transport*.  §4.5 (curvilinear angular redistribution);
  §5.3 (Diamond Difference, weighted-DD, Step, Linear Discontinuous;
  the negative-flux failure mode).  ⚠ §4.5 does **not** prescribe the
  retired :math:`[\tfrac12, 1]` clamp — that mis-citation is recorded at
  ``docs/theory/foundations/structured_geometry.rst
  §sn-tau-absorber-retirement``.

See also
========

* :class:`~orpheus.transport.spatial.scheme.DiscretizationScheme` — the
  Protocol this strategy satisfies.
* :func:`cell_balance_terms` — the unified algebra (Step 2.5).
* :doc:`/theory/methods/sn/index`, "Cell update strategies (the
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
    :math:`\psi_{\rm out} = a\,\psi_{\rm in} + b` with blend weight ``w = ½``
    (the symmetric diamond mean): the cell-average is affine in the SINGLE
    upstream face flux, so the DAG-free scan schedules (``CumprodScan``,
    ``ScanMarch``) consume DD via its coefficient triple
    :meth:`affine_scan_coefficients` and the generic base reconstruction
    staticmethods (#158 the coefficient model)."""

    transverse_coupling_is_facewise: ClassVar[bool] = True
    r"""DD's multi-D transverse coupling is **facewise** (separable): a non-swept
    axis :math:`y` enters the x-recurrence as the single 0th-order face value
    ``s_y · ψ_{y,in}`` folded into the scan's affine source (the explicit
    per-axis left-fold in :meth:`update`), so the d-D DD closure factors into
    independent per-axis 1-D scans and DD is admitted to the :math:`d \ge 2`
    scan-march (``ScanMarch``).  Contrast Linear-Discontinuous (bilinear
    slope-moment coupling → ``False``, rides the DAG wavefront); see the base
    trait for the facewise-vs-slopewise contract (#240 D5b / #38)."""

    diffusion_limit_consistent: ClassVar[bool] = True
    r"""DD's thick-diffusion limit IS a consistent diffusion discretization for
    the leading-order scalar flux (Larsen–Morel–Miller 1987 Eq. (4.24)).  ⚠ This
    is the SPATIAL axis: DD-in-ANGLE's first-order :math:`\beta`-failure (the
    curvilinear flux dip, Bailey–Morel–Chang 2010) is a DISTINCT, angular result
    — do NOT read it as a spatial-DD deficiency.  The angular condition lives on
    the pole-angular closure; the PAIR validity is
    :func:`~orpheus.sn.sweep.pairing.pair_diffusion_limit_consistent`."""

    supports_curvilinear: ClassVar[bool] = True
    r"""DD has a curvilinear cell closure: :meth:`update` runs the Morel–Montry
    angular redistribution (the ``angular_upstream is not None`` branch) for
    sphere/cylinder, and DD rides ``CumprodScan`` on every 1-D geometry.  So a
    curvilinear mesh may select a DD scheme (the default for sphere/cylinder)."""

    # ``has_transpose_kernel`` is DERIVED True — DD registers
    # ``streaming_cell_transpose`` (the relocated diamond-chain VJP, below
    # with the batch kernel pair); the trait follows the registration
    # (#310 ruling 2), never a declaration.

    def update(
        self,
        visit: CellVisit,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> CellResult:
        r"""Compute the cell-average flux + downstream states.

        One body — no geometry dispatch.  See the module docstring for the
        three structural observations enabling the collapse.
        """
        # ── Cell-balance solve: ONE formula, all geometries ─────────
        # The Morel--Montry constants c_in / c_out are angular-closure-owned
        # and arrive as DATA on the visit; cell_balance_terms consumes them —
        # it must NOT rebuild them from st.alpha_* / st.tau_mm (that would
        # re-fuse the spatial and angular closures).
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
        # ``ψ^a_out = (ψ_avg − (1−τ)·ψ^a_in)/τ`` exactly, with τ the
        # angular-closure-owned weight sourced off the visit (CellVisit.tau).
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

        This method delegates to :func:`cell_balance_for_streaming` (the
        vectorized helper the SN matvec also consumes) at ``n_mask=1``;
        :meth:`update` routes through :func:`cell_balance_terms` (the scalar
        solve-direction form).  Both compute the same intermediates, so at
        ``n_mask=1`` the per-ordinate result matches the scalar form
        bit-for-bit (Pattern 2, ONE algebra source).  At the converged
        ``cell_avg = update(...).cell_average_flux``, the residual is
        ``denom · cell_avg − (source + numer_upstream)`` = zero by the
        cell-balance equation ``update`` solved.
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

        # The M-M algebra is closure-owned: ``cell_balance_for_streaming``
        # takes only the ``(angular_denom_term, angular_numer_upstream)``
        # contributions, not raw ``c_in`` / ``c_out``.  The weighted-diamond
        # constants ``c_in`` / ``c_out`` arrive as DATA on the CellVisit; DD
        # must NOT rebuild them from ``st.alpha_*`` / ``st.tau_mm`` (it stays
        # geometry- AND closure-blind).  The (ΔA/w)-scaling below is the
        # geometry-owned redistribution factor.
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
    # enters the Cartesian producers (cell_kernel_batch / residual_kernel_batch
    # / cartesian_scan_coefficients).  Keep the diamond constant with exactly
    # one home here — never re-inline the ``2 g`` fold (Cardinal Rule 2).

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
        axis).  The **single source of the DD cell math** (Cardinal Rule 2):
        the SOLVE arm of the ``_CellSolve`` level operation, consumed
        identically by both storage walks
        (:meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_full`
        the verification oracle,
        :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed`
        the production window) — they differ ONLY in gather/scatter, so they
        cannot drift.

        Axis convention
        ---------------

        Element ``a`` of every per-axis tuple (``psi_in``, ``s_axes``, and the
        returned ``psi_out``) is spatial axis ``a`` — the SAME axis order as
        :attr:`~orpheus.sn.loss_representation.sweep_graph.OctantLabel.signs` and
        :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.from_cartesian`'s
        ``shape``. The tuples are positional-by-axis; the caller MUST build
        and unpack them in axis order (``[0]`` = x, ``[1]`` = y, ``[2]`` = z).

        The math is :math:`\psi_{\rm avg} = (Q + \sum_a 2 g_a\,\psi^{\rm in}_a) /
        (\Sigma_t + \sum_a 2 g_a)`, closure :math:`\psi^{\rm out}_a = 2\psi_{\rm
        avg} - \psi^{\rm in}_a`, where ``s_a = g_a`` is the RAW down-face
        streaming and the **scheme** applies its diamond factor ``2 = 1/w_DD``
        to BOTH the denom AND the upstream-numerator term (#240; a denom-only
        ``2`` would be a non-uniform 2× bug → wrong ``ψ̄``).  The streaming
        diagonal + the ``2 g_a`` couplings come from
        :meth:`_cartesian_streaming_diagonal` (the single source of DD's ``2 g``
        fold, carrying the explicit-left-fold bit-identity discipline — do NOT
        regroup).
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
        :meth:`_cartesian_streaming_diagonal`.
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

    def residual_kernel_batch_transpose(
        self,
        *,
        res_bar: np.ndarray,                  # (N_oct, ng, n_diag) — residual cotangent
        psi_out_bar: tuple[np.ndarray, ...],  # d arrays, each (N_oct, ng, n_diag)
        s_axes: tuple[np.ndarray, ...],       # d arrays, each (N_oct, 1, n_diag) — RAW g per axis
        reaction_xs: np.ndarray,              # (ng, n_diag)
    ) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
        r"""Batched VJP of :meth:`residual_kernel_batch` — dimension-generic.

        The exact reverse-mode pair of the DD apply kernel: with ``denom =
        \Sigma_t + \sum_a 2 g_a`` and the couplings ``2 g_a`` from the SAME
        :meth:`_cartesian_streaming_diagonal` fold the forward uses
        (Pattern 2 — no twin coefficients),

        .. math::

           \bar\psi^\dagger = \text{denom}\cdot r^\dagger
             + \sum_a \tfrac{1}{w}\,\psi_a^{\mathrm{out}\,\dagger}, \qquad
           \psi_{a,\mathrm{in}}^\dagger
             = -\tfrac{1-w}{w}\,\psi_a^{\mathrm{out}\,\dagger}
               - 2 g_a\, r^\dagger ,

        the face-chain pair riding the w-generic
        :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.outgoing_face_from_average_transpose`
        at ``w = _DD_W`` (#311 — at ``½`` the pair is the exact
        ``(2ψ_out†, −ψ_out†)``).  Consumed by the 1-D reverse loop walk's
        Cartesian arm (#310 C2) and the multi-D reverse wavefront's
        ``_CellResidualTranspose`` (#310 C3); the curvilinear reverse rides
        :meth:`streaming_cell_transpose` instead (the cell-balance arm).
        Registering this kernel is the Cartesian conjunct of the derived
        ``has_transpose_kernel``.
        """
        denom, couplings = self._cartesian_streaming_diagonal(reaction_xs, s_axes)
        psi_bar_cot = denom * res_bar
        psi_in_cots = []
        for c_a, out_bar_a in zip(couplings, psi_out_bar):
            avg_cot, in_cot = self.outgoing_face_from_average_transpose(
                out_bar_a, _DD_W,
            )
            psi_bar_cot = psi_bar_cot + avg_cot
            psi_in_cots.append(in_cot - c_a * res_bar)
        return psi_bar_cot, tuple(psi_in_cots)

    def streaming_cell_transpose(
        self,
        *,
        res_bar: np.ndarray,
        psi_out_bar: np.ndarray,
        denom: np.ndarray,
        abs_mu_A_total: np.ndarray,
        volume: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""DD's per-cell VJP — the diamond-chain transpose (relocated, bit-exact).

        The reverse-mode adjoint of the streamed DD cell relation
        {:math:`m = (\text{denom}\,\bar\psi - |\mu|A_{\rm tot}\,\psi_{\rm in}
        - \text{numer})/V`, :math:`\psi_{\rm out} = 2\bar\psi - \psi_{\rm in}`}:

        .. math::

           \bar\psi^\dagger = 2\,\psi_{\rm out}^\dagger
             + \text{denom}\cdot m^\dagger/V, \qquad
           \psi_{\rm in}^\dagger = -\psi_{\rm out}^\dagger
             - |\mu|A_{\rm tot}\cdot m^\dagger/V .

        Relocated from the 2.5a reverse-walk visit closure (#310 C1); the
        diamond-chain pair ``(2ψ_out†, −ψ_out†)`` rides the w-generic
        :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.outgoing_face_from_average_transpose`
        at ``w = _DD_W`` (#311) — byte-identical to the hand-transposed
        spelling (``÷½`` is an exact ``×2``; ``(1−½)/½ = 1`` exactly) — and
        the residual-row pullback is verbatim, so the frozen ``walk_matvec_*``
        adjoint baselines pin the whole kernel at 0 ULP.  ``denom`` arrives
        from the SAME ψ-independent ``cell_balance_for_streaming`` the forward
        uses (Pattern 2 — no twin algebra); the angular-numerator cotangent is
        the WALK's (spatial-only contract, #310 ruling 1).  This override IS
        the registration that derives ``has_transpose_kernel = True``.
        """
        psi_bar_cot, psi_in_bar = self.outgoing_face_from_average_transpose(
            psi_out_bar, _DD_W,
        )
        psi_bar_cot += denom * res_bar / volume
        psi_in_bar += -(abs_mu_A_total)[None, :] * res_bar / volume
        return psi_bar_cot, psi_in_bar

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
        (:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.cell_average` /
        :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.outgoing_face_from_average`
        / :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.source_emission`)
        consume it; DD carries no cell-average / outgoing-face / source-
        emission method of its own.

        The single source of the DD affine-recurrence coefficients (Cardinal
        Rule 2): :class:`CollisionCache.from_geometry` consumes this to populate
        the σ_t-stratum cache, and through it the ``CumprodScan`` / ``ScanMarch``
        sweep bodies.

        The math (Lewis & Miller §5.3; Hébert §3.9.3/§3.9.4 for the
        curvilinear curvature terms; the M-M weight :math:`\tau` entering
        through ``c_in`` / ``c_out`` is BMC 2010 Eqs. (42)/(43) — see the
        module References):

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

        Both compute "the cell-balance diagonal ``S``", but the REALIZATIONS are
        different mathematical objects, so they are NOT folded into one helper
        (Cardinal Rule 2 reconsidered): the Cartesian
        :meth:`_cartesian_streaming_diagonal` is ``S = Σ_t + Σ_a 2 g_a`` (÷Δ raw
        streaming, no volume, no angular coupling), while the curvilinear form
        here is ``S = Σ_t·V + 2|μ|·A_down + (ΔA/w)·c_out`` (×V collision,
        face-area streaming, AND the Morel–Montry curvature redistribution that
        couples the spatial diagonal to the angular closure — no Cartesian
        analogue, a different bit-identity pin).  The genuine unification is
        ``S = Σ_t·V + streaming_diag`` with a geometry-parameterised
        ``streaming_diag`` — the generic advection–reaction diagonal the
        diffusion scheme (#240's next model / the consistent-DSA ``A_diff``)
        will build; the merge belongs there, once two real occupants exist
        (``coding-elegance`` Pattern 6).  Deferred to #242 — forcing the
        curvilinear producer through the Cartesian helper now would re-baseline
        its snapshots for no gain.
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
        :meth:`_cartesian_streaming_diagonal` (passing
        ``s_axes = (s_scan, *s_transverse)`` so the scan-first order coincides
        with :meth:`cell_kernel_batch`'s axis order), carrying its bit-identity
        discipline — the DD scan snapshots stay strict.
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

    def reflect_scan_coefficients_transpose(
        self, res_bar: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""DD apply-transpose reflection scan coefficients ``(α = −1, β_pullback = 2)``.

        The reverse-direction pair of :meth:`reflect_scan_coefficients`
        (#310 C4 — the row-march reverse), single-sourced as the unit
        application of the w-generic #311 VJP
        :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.outgoing_face_from_average_transpose`
        at ``w = _DD_W``: ``(β_pullback, α) = VJP(1, ½) = (2, −1)`` exactly
        (``÷½`` is an exact ``×2``; ``(1−½)/½ = 1`` exactly, so ``×(−1)``
        is the exact reflection).  ``res_bar`` is a shape template.
        """
        ones = np.ones_like(res_bar)
        beta_pullback, alpha = self.outgoing_face_from_average_transpose(
            ones, _DD_W,
        )
        return alpha, beta_pullback

    # ── DiamondDifference is pure cell algebra in THREE capability groups ──
    #
    # Storage (gather → kernel → scatter) is the WALK's concern, not the
    # discretization's; this class supplies only cell math, in three groups —
    # one per sweep-schedule family:
    #   1. the per-cell reference pair (``update`` / ``residual``) — the
    #      canonical contract every scheme MUST provide;
    #   2. the batched kernel pair (``cell_kernel_batch`` /
    #      ``residual_kernel_batch``) — the DAG wavefront family (and the matvec
    #      apply twin of the scan, since the apply direction has a concrete ψ̄);
    #   3. the scan-family coefficients ``affine_scan_coefficients`` →
    #      ``(a, inverse_denom, w)`` — the DAG-free scan SOLVE family
    #      (CumprodScan / ScanMarch), consumed by the generic base
    #      reconstruction staticmethods on
    #      :class:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase`
    #      (#158 the coefficient model).
    # This is the ONLY direction-aware math in the SN stack, and the override
    # point for future closure strategies (Step / LD / EC supply group 2 and,
    # if affine-scannable, the group-3 coefficients).



__all__ = ["DiamondDifference"]
