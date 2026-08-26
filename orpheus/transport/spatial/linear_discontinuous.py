r"""Linear-Discontinuous (LD) cell-update strategy — slab / Cartesian (Issue #158).

The first higher-order, O(h²) occupant of the swappable per-cell
spatial-strategy contract (:mod:`orpheus.transport.spatial.scheme`), sibling to
the shipping :class:`~orpheus.transport.spatial.diamond.DiamondDifference` (DD).
**Full** LD (with the slope SOURCE moment :math:`\hat Q` carried) is
diffusion-limit-consistent — Larsen-Morel-Miller 1987 proves Step has **no**
valid intermediate diffusion limit while LD has all four, the load-bearing
reason a non-DD spatial scheme can lift the curvilinear pole-cell ``O(h)`` floor
(Issue #233).

.. note::

   **Diffusion-limit status (current implementation).**  The scattering slope
   source :math:`\Sigma_s\hat\phi` IS now threaded through the global
   spatial-moment iterate (Increment C / #240 D5b-S3; the per-ordinate-slope
   global-frame storage is
   ``docs/theory/methods/sn/cartesian_multid.rst §ld-ubld-sweep-global-frame``).
   Full LD therefore **recovers the thick-diffusion limit** —
   :class:`LinearDiscontinuous` declares ``diffusion_limit_consistent = True``
   and the limit is PINNED (no longer xfail) by
   ``tests/sn/verification/mms/test_mms_ld_slab.py::test_ld_thick_diffusive_limit``
   (1G) + ``::test_ld_thick_diffusive_limit_2g`` (2G).  ⚠ Remaining gap
   (Issue #247): the EXTERNAL slope source :math:`\hat Q^{\rm ext}` is still
   zeroed — the scattering channel EXERCISES the slope-source code path but does
   not CONSTRAIN its sign (vv-principles Mode 10, "activated-but-unconstrained"),
   so a manufactured solution with a non-vanishing external slope source is not
   yet supported.  This does NOT affect the diffusion limit, which is
   scattering-driven (:math:`\hat Q^{\rm ext} = 0`).

Why LD carries two moments
==========================

DD has a single per-cell unknown — the cell-average flux :math:`\bar\psi`.  LD
represents the in-cell angular flux as a **linear function** (Larsen & Morel
1989, JCP 83(1):212-236, Eqs. 4.1a-c), so each (group, ordinate) cell carries
**two** spatial moments (average :math:`\bar\psi` + slope :math:`\hat\psi`) with
an **upwind/discontinuous** face closure :math:`\psi_{\rm out} = \bar\psi +
\hat\psi` (LM-1989 Eq. 4.3c; after ORPHEUS pre-resolves the sweep direction the
strategy sees only the "downstream" face — no continuity enforced, the
"discontinuous" in LD).  The moment definitions and the 2×2 cell system
(:math:`\theta = 1/3`) are
``docs/theory/foundations/discretization.rst §discretization-ld``.

⚠ The off-diagonal / RHS signs of the **slope row** are the documented
correctness trap (LM-1989 memo §1.4/§6 — the published boxed form is internally
inconsistent).  The system ORPHEUS uses was regenerated symbolically with SymPy
and validated against the strongest oracle — **LD is exact on a linear-in-x
flux** (``ψ̄, ψ̂, ψ_out`` recovered to machine precision for any
:math:`\psi = a + bx`) — the foundation tests in
``tests/transport/spatial/test_linear_discontinuous.py``.

The Schur-complement scalar contract
====================================

The slope :math:`\hat\psi` is a **local** quantity eliminated by the Schur
complement of the 2×2, collapsing the cell update to a *scalar* relation
:math:`S\,\bar\psi = \mathrm{eff\_source} + \mathrm{eff\_numer\_upstream}`
(:math:`D_2' = \Sigma_t\theta h + |\mu|`; the closed forms are
``docs/theory/foundations/discretization.rst §discretization-ld``, label
``discretization-ld-schur``).  This is exactly the shape of DD's
``denom·ψ̄ = source + numer_upstream``, so LD fits the **existing scalar**
:meth:`residual` contract with no contract change — :math:`\hat\psi` is
reconstructed locally inside :meth:`update`.  The slope *source* :math:`\hat Q`
is carried on ``source[1]``; in a **fixed-source** problem it is supplied
directly (manufactured), which is why LD is staged on the fixed-source MMS
first.  Threading the slope through a *source-iteration* global flux iterate
(assembling the within-group :math:`\propto\Sigma_s\hat\phi` slope source) is
the deferred global moment-contract (Issue #158).

Scope (slab / Cartesian only — for now)
=======================================

This occupant implements the **slab/Cartesian** LD only.  The curvilinear
(sphere / cylinder) LD cell update — where the Morel-Montry angular weight
:math:`\tau` enters the average-moment denominator and the radial moment
:math:`(r-r_j)` weighting produces a slope-curvature coupling — is **not
published** and must be derived.  Until then, :meth:`update` / :meth:`residual`
raise :exc:`NotImplementedError` on a curvilinear visit (signalled by
``upstream_state.angular_upstream is not None``, exactly as DD signals the
presence of angular redistribution).

References
==========

* Larsen, E. W. & Morel, J. E. (1989). *Asymptotic solutions of numerical
  transport problems in optically thick, diffusive regimes II.*  JCP
  83(1):212-236.  §IV, Eqs. (4.1)-(4.3) — the slab-LD moment system + the
  full diffusion-limit proof.
* Larsen, E. W., Morel, J. E. & Miller, W. F. (1987).  *Asymptotic solutions
  of numerical transport problems in optically thick, diffusive regimes.*  JCP
  69(2):283-324.  Table I p.287 + §IX p.321 — the "LD has all four diffusion
  limits, Step has none" verdict.
* Lewis, E. E. & Miller, W. F. (1984).  *Computational Methods of Neutron
  Transport.*  §5.3 — LD, the positivity-fixup menu, the negative-flux mode.

See also
========

* :class:`~orpheus.transport.spatial.scheme.DiscretizationScheme` — the Protocol this
  strategy satisfies; :class:`~orpheus.transport.spatial.diamond.DiamondDifference`
  — the single-moment sibling whose ``update`` / ``residual`` shape this
  mirrors.
* :doc:`/theory/methods/sn/index`, "Cell update strategies" → "Linear
  Discontinuous".
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

import numpy as np

from orpheus.numerics.moment_layout import (
    AVERAGE_MOMENT,
    cell_moment_count,
    face_moment_count,
    face_moment_tail,
)

from ._ubld import (
    D1ClosedForm,
    assemble_inflow_axis,
    assemble_inflow_axis_transpose,
    assemble_ubld,
    d1_closed_form,
    mass_1d,
    per_cell_solve,
)
from .scheme import (
    CellResult,
    DiscretizationSchemeBase,
    CellVisit,
    UpstreamState,
)


def _require_slab(upstream_state: UpstreamState) -> None:
    r"""Guard: this occupant implements the slab/Cartesian LD only.

    Curvilinear geometry carries an upstream angular half-flux (the
    Morel-Montry redistribution thread); slab does not.  The presence of
    ``angular_upstream`` is therefore the geometry gate — identical to the
    signal DD uses to decide whether to run its angular closure.
    """
    if upstream_state.angular_upstream is not None:
        raise NotImplementedError(
            "LinearDiscontinuous currently implements the slab/Cartesian LD "
            "only; the curvilinear (sphere/cylinder) LD cell update is "
            "not yet implemented (Issue #158, curvilinear arm). "
            "A curvilinear visit was detected (angular_upstream is not None)."
        )


def _ld_source_moments(source: np.ndarray, ng: int) -> tuple[np.ndarray, np.ndarray]:
    r"""Split the LD moment source ``(2, ng)`` into ``(average, slope)``.

    LD's ``source`` is the cell source projected onto LD's spatial moment
    basis: ``source[0]`` is the prepared **average-moment** source
    :math:`\bar S = \bar Q\,h/(\sum_n w_n)` and ``source[1]`` the prepared
    **slope-moment** source :math:`\hat S = \hat Q\,h/(\sum_n w_n)` (both on
    the already-weight-normalised convention the contract documents for
    ``source``).  DD's single-moment ``(ng,)`` source is the ``n_moments=1``
    case; LD requires the explicit two-moment array — validated here so a
    DD-shaped source fails loudly rather than silently dropping the slope.
    """
    if source.shape != (2, ng):
        raise ValueError(
            "LinearDiscontinuous requires a two-moment source of shape "
            f"(2, ng) = (2, {ng}) — source[0]=average-moment, "
            f"source[1]=slope-moment; got shape {source.shape}. "
            "(A DiamondDifference-shaped (ng,) source drops the LD slope "
            "moment — see the LD source-moment contract.)"
        )
    return source[0], source[1]


@dataclass(frozen=True, slots=True)
class _LDCellTerms:
    r"""Per-cell LD Schur intermediates — the single source of the slab LD
    algebra for BOTH the solve (:meth:`LinearDiscontinuous.update`) and apply
    (:meth:`LinearDiscontinuous.residual`) directions.

    Eliminating the slope :math:`\hat\psi` by the Schur complement of the 2×2
    leaves a scalar relation in :math:`\bar\psi`,
    :math:`S\,\bar\psi = \mathrm{rhs}`.  The slope is reconstructed from
    :math:`\bar\psi` by :meth:`slope` — the ONE site of LD's slope-row sign
    convention (the documented correctness trap, LM-1989 memo §1.4/§6), so the
    sign lives in exactly one place rather than being re-spelled in ``update``.
    """

    eff_denom: np.ndarray            #: ``(ng,)`` Schur effective denominator S
    eff_source: np.ndarray           #: ``(ng,)`` source folded through the slope row
    eff_numer_upstream: np.ndarray   #: ``(ng,)`` inflow folded through the slope row
    _mu: np.ndarray                  #: ``|μ|·A_down`` (= ``|μ|`` for slab; slope reconstruction)
    _d2p: np.ndarray                 #: ``(ng,)`` :math:`D_2' = \Sigma_t\theta h + |\mu|`
    _slope_source: np.ndarray        #: ``(ng,)`` :math:`\theta\hat S` (slope-row RHS source)
    _psi_in: np.ndarray              #: ``(ng,)`` inflow face flux

    @property
    def rhs(self) -> np.ndarray:
        r"""Effective RHS of the Schur-reduced scalar balance (source + inflow)."""
        return self.eff_source + self.eff_numer_upstream

    def cell_average(self) -> np.ndarray:
        r"""Solve the Schur-reduced scalar balance for :math:`\bar\psi`."""
        return self.rhs / self.eff_denom

    def slope(self, cell_avg: np.ndarray) -> np.ndarray:
        r"""Reconstruct :math:`\hat\psi` from :math:`\bar\psi` via the
        slope-moment row :math:`\hat\psi = (|\mu|\bar\psi + \theta\hat S
        - |\mu|\psi_{\rm in})/D_2'` — LD's sole slope-row sign-convention site."""
        return (self._mu * cell_avg + self._slope_source
                - self._mu * self._psi_in) / self._d2p


class LinearDiscontinuous(DiscretizationSchemeBase, key="linear_discontinuous"):
    r"""Linear-Discontinuous (LD) cell-update strategy — slab / Cartesian.

    Two spatial moments per (group, ordinate) cell: the cell-average flux and
    the cell-average slope, with the upwind-discontinuous face closure.  The
    slope is eliminated locally by the Schur complement so the update fits the
    scalar :class:`~orpheus.transport.spatial.scheme.DiscretizationScheme` contract.  See
    the module docstring for the derivation, the SymPy-verified 2×2, and the
    Schur-complement scalar form.

    Not a frozen dataclass-with-fields like DD: LD carries a single immutable
    parameter :attr:`theta` (the slope-moment weight) as a class constant.
    """

    is_linear: ClassVar[bool] = True
    is_positivity_preserving: ClassVar[bool] = False
    is_affine_scannable: ClassVar[bool] = True
    r"""LD admits the single-upstream affine recurrence
    :math:`\psi_{\rm out}=a\psi_{\rm in}+b` (#158 Increment B): the slope
    :math:`\hat\psi` is Schur-eliminated locally, leaving ``a``
    source-independent, so LD supplies the scan coefficient triple via
    :meth:`affine_scan_coefficients` and rides ``CumprodScan`` (and the **1-D**
    ``ScanMarch``).  Slab/Cartesian only (raises on curvilinear).  This is the
    single-axis trait; LD leaves ``transverse_coupling_is_facewise`` at
    ``False`` because its multi-D closure is BILINEAR (an independent slope per
    axis → slope-wise coupling, not separable), so LD rides the DAG wavefront in
    :math:`d \ge 2`, not the scan-march (#38 / #240 D5b)."""

    spatial_basis_per_axis: ClassVar[int] = 2
    r"""LD carries TWO spatial moments per axis — the Legendre basis
    :math:`\{1, P_1\}` (cell-average + slope) — so the per-cell unknown is
    :math:`2^d` (d=1: ``{ψ̄, ψ̂}``; d=2: the bilinear ``{ψ̄, ψ̂_y, ψ̂_x, ψ̂_xy}``)
    and each downstream face carries :math:`2^{d-1}` transverse moments.  The
    ``per_axis > 1`` gate routes LD onto the multi-moment face-cochain + the
    moment-reducing emit (#240 D5b); DiamondDifference at ``per_axis == 1`` keeps
    the scalar path byte-identical.  The tensor-Legendre Kronecker layout (slot
    order, the diffusion-limit-load-bearing cross moment ``ψ̂_xy``) is
    ``docs/theory/methods/sn/cartesian_multid.rst §spatial-moment-space``; the
    d-generic primitive is :mod:`orpheus.transport.spatial._ubld`."""

    diffusion_limit_consistent: ClassVar[bool] = True
    r"""Full LD's thick-diffusion limit IS a consistent diffusion discretization
    (Larsen–Morel 1989 Part II §IV Eqs. (4.16)-(4.19); LMM-1987 §VII "all four
    diffusion limits") — the load-bearing reason a non-DD spatial scheme can lift
    the curvilinear pole-cell :math:`O(h)` floor (#233).  ⚠ The property is a
    statement about FULL LD — it REQUIRES the slope SOURCE moment
    :math:`\Sigma_s\hat\phi` threaded through the iterate (now landed; the limit
    is PINNED, no longer xfail, by ``test_ld_thick_diffusive_limit`` — see the
    module note).  A flat-source LD (:math:`\hat Q = 0`) would be ``O(h²)`` but
    NOT diffusion-limit-consistent.  The angular first-order condition is
    separate (the pole-angular closure's ``beta_first_order_consistent``); the
    PAIR validity is
    :func:`~orpheus.sn.sweep.pairing.pair_diffusion_limit_consistent`."""

    supports_curvilinear: ClassVar[bool] = False
    r"""LD is slab/Cartesian ONLY: the curvilinear (sphere/cylinder) LD cell
    closure is not yet implemented here (Issue #158 curvilinear arm /
    #6) — the derivation IS published: Adams-Martin 1992, NSE 111, App. A
    (the 1-D spherical LD moment balances, with the weighted-diamond
    angular closure applied per spatial moment — average AND slope;
    ⚠ the printed slope-term minus signs are a typo, the redistribution
    Gram is symmetric positive-definite, confirmed independently by
    Machorro 2007, JCP 223, and Hill 1975, LA-5990-MS §ONETRAN).
    ⛔ Do NOT implement the bare (A.1) form: Palmer-Adams 1993
    (UCRL-ID-114256) find plain spherical LD FAILS the thick diffusion
    limit (three-point removal term, unphysical boundary conditions,
    interior scalar flux low by ~2x); the fully-lumped and
    corner-balance forms pass.  :meth:`update` / :meth:`residual` raise on
    a curvilinear visit (via :func:`_require_slab`).  The conservative ``False`` (= the base default,
    declared explicitly here for the citation) makes the sweep-strategy
    selection reject a curvilinear-LD mesh at construction
    (:meth:`~orpheus.sn.loss_representation.CumprodScan.supports`) instead of
    passing on ``is_affine_scannable`` (a geometry-blind 1-D trait) and raising
    mid-sweep."""

    # ``has_transpose_kernel`` is DERIVED True (#310 C2): LD registers the
    # UBLD ``AᵀM⁻¹`` batch VJP (``residual_kernel_batch_transpose``, below),
    # which ALONE covers its slab-only forward span under the covering
    # derivation — no ``streaming_cell_transpose`` override is needed (that
    # is the CURVILINEAR arm's kernel, and LD claims no curvilinear
    # support).  The trait flips by REGISTERING the kernel, never by
    # declaration (ruling 2).  An eager ``.H`` on a 1-D LD mesh CONSTRUCTS;
    # the remaining loud backstop is the reverse entries' spatial-moment
    # tail check (a tail-mismatched cotangent raises a typed ValueError,
    # never a silent broadcast).

    theta: ClassVar[float] = 1.0 / 3.0
    r"""The slope-moment weight :math:`\theta` (LM-1989 Eq. 4.3b).  The value
    :math:`\theta = 1/3` is the SN-exact linear-discontinuous closure (LM-1989
    sets it for the consistent LD; it is a free parameter only in their
    asymptotic analysis).  Do NOT change to the mass-lumped (LLD) value — that
    is a different scheme."""

    # ── Schur-complement scalar intermediates (single source for both
    #    update and residual — Pattern 2, Cardinal Rule 2) ────────────────
    def _schur_terms(
        self,
        visit: CellVisit,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> _LDCellTerms:
        r"""Build the per-cell Schur intermediates :class:`_LDCellTerms`.

        The single algebra site consumed by BOTH :meth:`update` (solve) and
        :meth:`residual` (apply); ``_ld_source_moments`` is called here ONCE
        and the slope-row reconstruction lives on the returned bundle.
        """
        st = visit.streaming_terms
        ng = total_xs.size
        s_bar, s_hat = _ld_source_moments(source, ng)

        # Step 2.5: every StreamingTerms field is a required ``float`` (all
        # three factories populate them for every geometry), so abs_mu /
        # volume need no None-narrowing — the type guarantees they are set.
        mu: float = st.abs_mu               # |μ|, sweep-pre-resolved (slab A=1)
        h: float = st.volume                # slab cell width (V = Δx)
        theta = self.theta
        psi_in = upstream_state.spatial_upstream

        # Single-source the LD 2×2 Schur through the shared d=1 closed form
        # (slab A_down = 1, V = h, so the ÷V streaming-over-volume is g = |μ|/h).
        # ``schur_xV`` returns the ×V per-cell intermediates (S, eff_source,
        # eff_numer, θ·ŝ, |μ|A_down, D₂'); the LD 2×2 algebra lives ONCE in
        # ``_ubld.d1_closed_form``, proven == the d-generic dense primitive's
        # d=1 reduction.
        cf = d1_closed_form(mu / h, total_xs, theta)
        eff_denom, eff_source, eff_numer_upstream, slope_source, mu_Adown, d2p = (
            cf.schur_xV(h, s_bar, s_hat, psi_in)
        )
        return _LDCellTerms(
            eff_denom=eff_denom,
            eff_source=eff_source,
            eff_numer_upstream=eff_numer_upstream,
            _mu=mu_Adown,
            _d2p=d2p,
            _slope_source=slope_source,
            _psi_in=psi_in,
        )

    def update(
        self,
        visit: CellVisit,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> CellResult:
        r"""Solve the LD cell system; return the average flux + outflow face.

        Solves the Schur-reduced scalar relation for :math:`\bar\psi`,
        reconstructs the slope :math:`\hat\psi` via :meth:`_LDCellTerms.slope`
        (the sole slope-row site), and returns the downstream face
        :math:`\psi_{\rm out}=\bar\psi+\hat\psi` (the slope is recoverable from
        the (average, outflow) pair — no ``CellResult`` field change).  Slab
        only (see :func:`_require_slab`).
        """
        _require_slab(upstream_state)
        terms = self._schur_terms(visit, total_xs, source, upstream_state)
        psi_bar = terms.cell_average()
        psi_out = psi_bar + terms.slope(psi_bar)
        return CellResult(
            cell_average_flux=psi_bar,
            outgoing_spatial_flux=psi_out,
            outgoing_angular_state=None,     # slab: no angular redistribution
        )

    def residual(
        self,
        cell_avg: np.ndarray,
        visit: CellVisit,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> np.ndarray:
        r"""Per-cell operator residual :math:`S\,\bar\psi - \mathrm{rhs}`.

        The apply direction of the Schur-reduced scalar system: at the
        :math:`\bar\psi` that :meth:`update` solves for, this vanishes to
        floating-point rounding (the round-trip contract,
        :meth:`DiscretizationScheme.residual`).  Linear in ``cell_avg``; affine in
        ``source``.  Slab only.
        """
        _require_slab(upstream_state)
        terms = self._schur_terms(visit, total_xs, source, upstream_state)
        return terms.eff_denom * cell_avg - terms.rhs

    # ── DAG-family batched kernel (group 2 — the polymorphic wavefront path) ──
    #
    # The cell kernel the DAG wavefront walks call (the SAME contract
    # DiamondDifference fills): :class:`FullFieldWavefront` — the dimension-
    # generic oracle, including 1-D slab — runs LD through these exactly as it
    # runs DD.  The walk owns the storage, the scheme owns only its cell algebra.
    #
    # Convention: the kernel uses the d-generic Cartesian streaming coefficient
    # ``s = 2|μ|/Δ`` per axis (``A = 1``, ``V = Δ``) — the SAME ``s_axes`` DD's
    # kernel consumes — so ``g = s/2 = |μ|/Δ`` is LD's streaming-over-volume.
    # This is the ÷V form of the per-cell ×V ``update``/``residual`` 2×2; the two
    # are the SAME LD, pinned consistent by the group1≡group2 gate.
    #
    # A matvec is a forward APPLY, and applying the per-cell ``2^d × 2^d``
    # operator to the moment vector is intrinsically MOMENT-valued in every d, so
    # the kernel pair rides ONE d-generic dense UBLD path (``_ubld_system`` →
    # ``per_cell_solve``) for all d — theory:
    # ``docs/theory/methods/sn/cartesian_multid.rst §ld-ubld-unified-moment-matvec``.
    # The L16 closed-form fast path lives in the d=1 PRODUCTION sweep (CumprodScan
    # → ``affine_scan_coefficients`` → ``d1_closed_form``), selected by strategy
    # on geometry, NOT by an ``if d==1`` in the kernel.  The ÷V dense system is
    # SCALE-FREE (fed unit widths ``hs = [1, …]`` and ``mus = [g_0, …]``); at d=1
    # it is EXACTLY ``d1_closed_form``'s 2×2 (proven == by
    # ``tests/transport/spatial/test_ld_ubld_primitive.py``).

    def _ubld_system(
        self,
        s_axes: tuple[np.ndarray, ...],
        reaction_xs: np.ndarray,
        Q_cells: np.ndarray,
    ) -> tuple[dict, np.ndarray]:
        r"""Assemble the batched ÷V UBLD dense system + the source-moment RHS.

        Returns ``(assembled, R_source)`` where ``assembled`` is the
        :func:`~orpheus.transport.spatial._ubld.assemble_ubld` dict (``A`` of shape
        ``(N_oct, ng, n_diag, 2^d, 2^d)``) and ``R_source`` the source
        contribution ``M·S⃗_moments`` of shape ``(N_oct, ng, n_diag, 2^d)``.  The
        upwind-inflow face contributions are added by the SOLVE / APPLY arms
        (they vary by direction).

        ``Q_cells`` is the per-cell source.  It is SHAPE-AGNOSTIC over the
        spatial-moment axis (#240 D5b-S3): a genuine ``(…, 2^d)`` moment source
        — average in slot 0, the scattering-slope source ``Σ_s·φ̂`` in the slope
        rows — flows through verbatim (its trailing axis IS the ``2^d`` moment
        axis); a flat scalar ``(N_oct or 1, ng, n_diag)`` source (the matvec
        ``_MATVEC_ZERO_SOURCE`` or a flat external source) is lifted onto the
        moment vector with the average moment in slot 0 and the rest zero.  The
        average-moment slot index is the single-sourced :data:`AVERAGE_MOMENT`.
        """
        d = len(s_axes)
        size = cell_moment_count(2, d)
        # g_a = s_axes[a] (N_oct, 1, n_diag); broadcast to (N_oct, ng, n_diag)
        # against reaction_xs (ng, n_diag).
        sig_b = np.asarray(reaction_xs)
        gs = [np.asarray(s) + np.zeros_like(sig_b) for s in s_axes]
        sig_full = sig_b + np.zeros_like(gs[0])
        ones = [np.ones_like(g) for g in gs]
        assembled = assemble_ubld(ones, gs, sig_full, self.theta)
        batch = assembled["A"].shape[:-2]   # (N_oct, ng, n_diag) — the cell stack
        Q = np.asarray(Q_cells)
        # Discriminate moment vs flat by RANK against the cell-stack batch
        # (unambiguous — the flat scalar source broadcasts to ``batch``, ndim
        # ≤ len(batch); the genuine moment source is ``batch + (2^d,)``, exactly
        # one more axis).  A trailing-axis-SIZE test would be fragile (a level's
        # cell count could equal ``2^d`` at small meshes).
        if Q.ndim == len(batch) + 1:
            # Genuine moment source — its trailing axis IS the 2^d moment axis
            # (the scattering-slope source Σ_s·φ̂ rides the slope rows; #240 S3).
            S_moments = Q + np.zeros(batch + (size,))
        else:
            # Flat scalar source — lift onto slot 0 (average), rest zero.
            S_moments = np.zeros(batch + (size,))
            S_moments[..., AVERAGE_MOMENT] = Q + np.zeros(batch)
        R_source = np.einsum("...ij,...j->...i", assembled["M"], S_moments)
        return assembled, R_source

    def _ubld_inflow(
        self,
        s_axes: tuple[np.ndarray, ...],
        psi_in: tuple[np.ndarray, ...],
    ) -> np.ndarray:
        r"""Sum the per-axis upwind-inflow RHS contributions (the ÷V form).

        Each ``psi_in[a]`` is the upstream neighbour's outflow trace on axis
        ``a`` — a ``2^{d-1}``-moment transverse object (the widened face
        cochain).  :func:`~orpheus.transport.spatial._ubld.assemble_inflow_axis`
        weights it into the active-axis test functions ``B(-1)=[1,-1]`` and the
        transverse mass, times ``g_a`` (the ÷V streaming, ``mus[axis]``).

        At d=1 the transverse-moment count is ``2^{d-1} = 1`` and the walk's
        face cochain carries NO explicit moment axis (``face_moment_tail(1) ==
        ()``), so each face arrives as ``(N_oct, ng, n_diag)``.  The dense
        inflow assembler needs the transverse-moment axis EXPLICIT (it is the
        face's LAST axis), so a trailing length-1 axis is appended here at d=1
        (#240 D5b-S3 — the unified moment matvec routes d=1 LD through the dense
        primitive).  The decision keys on ``face_moment_tail(2^{d-1}) == ()``
        (deterministic on ``d``, NOT on a shape probe — a scalar face's
        trailing ``n_diag == 1`` is indistinguishable from a length-1 moment
        axis by shape alone).  This is LD-only (DD never reaches this method),
        so it does not touch DD's byte-identical scalar face convention.
        """
        d = len(s_axes)
        # The walk omits the moment axis exactly when its tail is empty (d=1).
        needs_moment_axis = face_moment_tail(face_moment_count(2, d)) == ()
        ones = [np.ones_like(g) for g in s_axes]
        gs = [np.asarray(g) for g in s_axes]
        faces = [
            np.asarray(psi_in[a])[..., None] if needs_moment_axis
            else np.asarray(psi_in[a])
            for a in range(d)
        ]
        terms = [
            assemble_inflow_axis(ones, gs, a, faces[a], self.theta)
            for a in range(d)
        ]
        R = terms[0]                                # d ≥ 1 ⇒ ≥ 1 axis, never empty
        for term in terms[1:]:
            R = R + term
        return R

    @staticmethod
    def _ubld_outgoing_faces(
        psi_moments: np.ndarray, d: int,
    ) -> tuple[np.ndarray, ...]:
        r"""Reconstruct the d downstream faces from the ``2^d`` moment vector.

        The outgoing face on axis ``a`` is the trace of the cell's
        tensor-Legendre function at the downstream node (``s_a = +1``, where
        ``P₀(+1)=P₁(+1)=1``), projected onto the ``2^{d-1}`` transverse Legendre
        moments.  In the Kronecker layout (axis 0 outer … axis d−1 inner) the
        ``2^d`` vector reshapes to ``(2,)*d`` indexed by ``(o_0, …, o_{d-1})``;
        the trace on axis ``a`` SUMS the ``o_a=0`` and ``o_a=1`` blocks (the
        ``B(+1)`` weights) and flattens the remaining transverse axes in their
        own Kronecker order — EXACTLY the order
        :meth:`_ubld_inflow`/:func:`assemble_inflow_axis` consumes on axis ``a``
        (so out-of-cell == in-of-next-cell; the closure consistency the matvec
        twin D5b.4 + the round-trip D5b.1 verify).
        """
        batch = psi_moments.shape[:-1]
        # The face tail follows the walk's cochain convention: a 2^{d-1}-moment
        # axis when > 1 (d ≥ 2), ABSENT at d=1 (face_moment_tail(1) == () — the
        # walk's scalar d=1 face cochain).  This mirrors _ubld_inflow's d=1
        # axis-append, keeping out-of-cell == in-of-next-cell consistent (#240
        # D5b-S3 — the unified moment matvec).
        tail = face_moment_tail(face_moment_count(2, d))
        tensor = psi_moments.reshape(*batch, *([2] * d))
        faces = []
        for a in range(d):
            cell_axis = len(batch) + a
            face = tensor.take(0, axis=cell_axis) + tensor.take(1, axis=cell_axis)
            faces.append(face.reshape(*batch, *tail))
        return tuple(faces)

    @staticmethod
    def _ubld_outgoing_faces_transpose(
        faces_bar: tuple[np.ndarray, ...], d: int,
    ) -> np.ndarray:
        r"""Adjoint of :meth:`_ubld_outgoing_faces` — the face-trace pullback.

        The forward trace on axis ``a`` SUMS the ``o_a = 0`` and ``o_a = 1``
        blocks (the ``B(+1)`` weights), so its adjoint BROADCASTS each face
        cotangent equally into both blocks (``np.stack`` duplicating along
        the axis-``a`` slot — the exact transpose of ``take(0) + take(1)``),
        transverse Kronecker order preserved, summed over axes.  The face
        tail follows the walk's cochain convention exactly as the forward
        (scalar at d=1 — a trailing length-1 axis is re-appended here,
        mirroring :meth:`_ubld_inflow`'s d=1 axis-append).  Returns the
        ``(..., 2^d)`` ψ⃗-cotangent contribution.  Symbolic ground: the
        face-pullback row of
        :func:`orpheus.derivations.discrete.sn.ld_ubld.derive_d1_transpose_equals_At_Minv`.
        """
        tail = face_moment_tail(face_moment_count(2, d))
        first = np.asarray(faces_bar[0], dtype=np.float64)
        batch = first.shape[: first.ndim - len(tail)]
        out = np.zeros(batch + (cell_moment_count(2, d),))
        tensor = out.reshape(*batch, *([2] * d))
        for a in range(d):
            fb = np.asarray(faces_bar[a], dtype=np.float64)
            if tail == ():
                fb = fb[..., None]                 # d=1: re-append the moment axis
            fb_t = fb.reshape(*batch, *([2] * (d - 1)))
            cell_axis = len(batch) + a
            tensor += np.stack([fb_t, fb_t], axis=cell_axis)
        return out

    def cell_kernel_batch(
        self,
        *,
        psi_in: tuple[np.ndarray, ...],
        s_axes: tuple[np.ndarray, ...],
        reaction_xs: np.ndarray,
        Q_cells: np.ndarray,
    ) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
        r"""Pure batched LD cell SOLVE — the d-generic DAG wavefront kernel.

        Assembles the batched ÷V ``2^d × 2^d`` dense UBLD system
        (``A ψ⃗ = M·S⃗ + Σ_a F_in^{(a)}``) over the level's anti-diagonal cell
        stack and solves it with one batched
        :func:`~orpheus.transport.spatial._ubld.per_cell_solve`.  Returns the full
        ``2^d``-moment ``psi_avg`` ``(N_oct, ng, n_diag, 2^d)`` — the cell-average
        in slot 0, the slope moments after — and the d-tuple of
        ``2^{d-1}``-moment downstream faces (one per axis).  ONE moment path for
        every d: at d=1 the system is the LD ``2×2`` whose Schur complement is
        :func:`~orpheus.transport.spatial._ubld.d1_closed_form` (the L16 production
        scan's fast path), proven equal by
        ``tests/transport/spatial/test_ld_ubld_primitive.py``.

        STORAGE-FREE by contract: the WALK gathers ``psi_in`` (per-axis
        ``2^{d-1}``-moment faces; a scalar at d=1) and scatters the outgoing
        faces.  The SOLVE arm of the ``_CellSolve`` level operation consumed by
        :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_full` /
        :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed`.
        """
        d = len(s_axes)
        assembled, R_source = self._ubld_system(s_axes, reaction_xs, Q_cells)
        R = R_source + self._ubld_inflow(s_axes, psi_in)
        psi_moments = per_cell_solve(assembled, R)
        return psi_moments, self._ubld_outgoing_faces(psi_moments, d)

    def moment_mass_diagonal(self, ndim: int) -> np.ndarray:
        r"""LD's unit-volume Kronecker mass diagonal ``∏_a θ^{o_a}`` — ``(2^d,)``.

        Built from the SAME 1-D factor :func:`~orpheus.transport.spatial._ubld.mass_1d`
        at unit width the UBLD assembler Kroneckers (Pattern 2): ``[1, θ]``
        at d=1, ``[1, θ, θ, θ²]`` at d=2.  This IS ``diag(M)/V`` — the
        diagonal :meth:`residual_kernel_batch` normalises by and the moment
        mass the ``.H`` bulk metric must carry (#310 C2 ruling 3).
        """
        base = np.diagonal(mass_1d(np.ones(()), self.theta))      # [1, θ]
        diag = np.array([1.0])
        for _ in range(ndim):
            diag = np.kron(diag, base)
        return diag

    def residual_kernel_batch(
        self,
        *,
        psi_bar: np.ndarray,
        psi_in: tuple[np.ndarray, ...],
        s_axes: tuple[np.ndarray, ...],
        reaction_xs: np.ndarray,
        Q_cells: np.ndarray,
    ) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
        r"""Pure batched LD operator residual — the d-generic apply twin.

        A matvec is a forward APPLY: it evaluates the bilinear UBLD residual
        :math:`r = A\,\vec\psi - R` at the PROBE ``2^d``-moment vector
        ``psi_bar`` ``(N_oct, ng, n_diag, 2^d)`` and reconstructs the d
        ``2^{d-1}``-moment downstream faces from the probe (so the matvec edge
        fluxes propagate along the wavefront exactly as the sweep's — matvec ≡
        sweep).  ONE moment path for every d — the residual is intrinsically
        moment-valued (average ROW + slope ROWS).  At the ``psi_bar``
        :meth:`cell_kernel_batch` solves for (same inputs) the residual vanishes
        to FP noise (the batched round-trip).  The APPLY arm of the
        ``_CellResidual`` level operation (the matvec walk).

        Mass-diagonal normalization (the matvec/sweep moment-source consistency).
        The UBLD RHS folds the cell source mass-weighted (``R_source = M·S⃗``, the
        test-function projection); but the operator-algebra ``A = (L+C) − S``
        subtracts the scattering moment source ``S.apply(ψ)`` RAW (per-ordinate,
        un-projected) at the ``OperatorSum`` level.  Dividing the residual by the
        (diagonal) mass ``M`` puts ``(L+C)ψ`` in RAW per-moment units so the
        row-for-row ``− S.apply(ψ)`` is consistent — the slope rows would
        otherwise disagree by the Legendre mass factor ``M_ii = θ^{|i|}`` (``|i|``
        = number of active slope axes).  The division preserves the round-trip
        (``A·ψ⃗ − R = 0 ⇒ (A·ψ⃗ − R)/M_ii = 0``) and matches the d=1 DD ÷V
        convention (``M_00 = 1`` so the average row is already raw).
        """
        d = len(s_axes)
        assembled, R_source = self._ubld_system(s_axes, reaction_xs, Q_cells)
        R = R_source + self._ubld_inflow(s_axes, psi_in)
        residual = np.einsum("...ij,...j->...i", assembled["A"], psi_bar) - R
        # M is diagonal in the Legendre basis (M = ⊗ diag(1, θ)); divide each
        # moment row by its mass so the residual is in raw per-moment units.
        mass_diag = np.diagonal(assembled["M"], axis1=-2, axis2=-1)
        residual = residual / mass_diag
        return residual, self._ubld_outgoing_faces(psi_bar, d)

    def residual_kernel_batch_transpose(
        self,
        *,
        res_bar: np.ndarray,
        psi_out_bar: tuple[np.ndarray, ...],
        s_axes: tuple[np.ndarray, ...],
        reaction_xs: np.ndarray,
    ) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
        r"""Batched VJP of :meth:`residual_kernel_batch` — the UBLD ``AᵀM⁻¹`` pullback.

        The exact reverse-mode pair of the LD apply kernel, generated from
        the SAME materialized record the forward assembles
        (:meth:`_ubld_system` — Pattern 2, zero hand-derived entries):

        * the ψ⃗ pullback is ``Aᵀ·(M⁻¹·r̄)`` — the mass inverse applies
          **BEFORE** ``Aᵀ`` (the forward normalises ``(Aψ⃗−R)/M_ii``, so
          the Jacobian is ``M⁻¹A`` and its transpose ``AᵀM⁻¹``; ``M``
          diagonal ⟹ ``M⁻ᵀ = M⁻¹``) — plus the face-trace broadcast
          :meth:`_ubld_outgoing_faces_transpose`;
        * the per-axis upstream-face pullback is
          ``−assemble_inflow_axis_transposeᵀ`` of the mass-normalised
          cotangent, mirroring :meth:`_ubld_inflow`'s factor construction
          (including its d=1 moment-axis append, stripped again on the way
          out so the returned face cotangents match the walk's scalar d=1
          cochain).

        Symbolic ground:
        :func:`orpheus.derivations.discrete.sn.ld_ubld.derive_d1_transpose_equals_At_Minv`
        (the mass-order law, with the order discriminant proven nonzero) —
        NEW-algebra (a)/(b) of the #310 C2 spec §3.3.  Source-free by the
        base contract (the source cotangent belongs to the solve path).
        **Registering this override IS the LD capability flip** — the
        covering ``has_transpose_kernel`` derivation turns ``True`` here
        (LD is slab-only, so the batch kernel alone covers its span), which
        is why this landing is atomic with the reverse-scan moment arm, the
        two-guard lift, and the rewritten deferral pins (#310 flip-safety).
        """
        d = len(s_axes)
        assembled, _ = self._ubld_system(s_axes, reaction_xs, np.zeros(()))
        mass_diag = np.diagonal(assembled["M"], axis1=-2, axis2=-1)
        r_bar = res_bar / mass_diag                       # mass-inverse FIRST
        psi_bar_cot = np.einsum("...ji,...j->...i", assembled["A"], r_bar)
        psi_bar_cot = psi_bar_cot + self._ubld_outgoing_faces_transpose(
            psi_out_bar, d,
        )
        # Mirror _ubld_inflow's factor construction (unit widths; ÷V
        # streaming as the per-axis μ) and its d=1 moment-axis convention.
        needs_moment_axis = face_moment_tail(face_moment_count(2, d)) == ()
        ones = [np.ones_like(np.asarray(g)) for g in s_axes]
        gs = [np.asarray(g) for g in s_axes]
        psi_in_cots = []
        for a in range(d):
            face_cot = assemble_inflow_axis_transpose(
                ones, gs, a, r_bar, self.theta,
            )
            if needs_moment_axis:
                face_cot = face_cot[..., 0]
            psi_in_cots.append(-face_cot)
        return psi_bar_cot, tuple(psi_in_cots)

    # ── Scan-family coefficients (group 3 — the DAG-free schedules) ──────────
    #
    # The ×V "denom" convention the scan machinery uses (CollisionCache /
    # ordinate_scan): the SAME LD as the ÷V group-2 kernel above, scaled by V
    # (S_scan = V·S_kernel, SymPy-verified), so CumprodScan-LD and
    # FullFieldWavefront-LD are principled-equivalent at nULP (the two-paths
    # gate).  Slab/Cartesian only — the curvilinear LD closure is not yet
    # implemented (#158; published: Adams-Martin 1992 App. A).

    def affine_scan_coefficients(
        self,
        *,
        abs_mu: np.ndarray,    # (N,)        |μ_n|
        A_down: np.ndarray,    # (N, nx)     downstream face area (slab: 1)
        A_total: np.ndarray,   # (N, nx)     A_inner + A_outer (unused; slab=2)
        dA_w: np.ndarray,      # (N, nx)     ΔA/w curvature redistribution (slab: 0)
        c_out: np.ndarray,     # (N, nx)     M-M outgoing closure const (slab: 0)
        V: np.ndarray,         # (N, nx)     cell volume per ordinate
        reaction_xs: np.ndarray,  # (N, ng, nx) Σ_t in the geometry's cell ordering
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        r"""LD's :math:`\Sigma_t`-epoch scan coefficients ``(a, inverse_denom, w)``.

        The Schur-reduced slab LD 2×2 (flat :math:`\hat Q = 0`) in the ×V
        convention (``m = |μ|·A_down``, ``t = Σ_t V``, ``p = m/θ``):

        .. math::

           D_2 &= t + p, \qquad k = p / D_2, \\
           S   &= (t + m) + m\,p / D_2,
                  \qquad \mathrm{inverse\_denom} = 1/S, \\
           a   &= m\,(1+k)^2 / S - k,
                  \qquad w = 1/(1+k).

        ``a`` is source-independent (the affine transmission); ``w`` is the
        cell-average blend weight (:math:`\bar\psi=(1-w)\psi_{\rm in}+w\psi_{\rm
        out}`).  All three SymPy-verified against the per-cell 2×2 + the ÷V
        :meth:`cell_kernel_batch` form (#158 Increment B).

        Slab-only guard
        ---------------

        The curvilinear (sphere/cylinder) LD scan closure is not yet
        implemented (#158).  Its
        signal is non-neutral curvature: slab carries ``dA_w == 0`` and
        ``c_out == 0`` EXACTLY (the :func:`slab_streaming` neutral element);
        curvilinear carries non-zero values.  Raising here fails fast at the
        :class:`~orpheus.sn.sweep.cache.CollisionCache` build
        (``SNSolver.__init__``) before any sweep or matvec runs — so a 1-D
        sphere/cylinder mesh carrying ``LinearDiscontinuous`` (which would match
        ``CumprodScan.supports`` via ``is_1d and is_affine_scannable``) is
        rejected loudly rather than silently running DD-shaped curvature math.
        """
        if np.any(dA_w != 0.0) or np.any(c_out != 0.0):
            raise NotImplementedError(
                "LinearDiscontinuous.affine_scan_coefficients supports the "
                "slab/Cartesian LD only; the curvilinear (sphere/cylinder) LD "
                "scan closure is not yet implemented (#158).  A non-neutral curvature "
                "was detected (dA_w / c_out are not all zero)."
            )
        # Single-source the LD 2×2 Schur through the shared d=1 closed form: the
        # ×V scan reads ``(a, inverse_denom, w)`` off the same helper the ÷V
        # kernel and ×V per-cell views use.  The ÷V streaming-over-volume is
        # g = |μ|·A_down/V; ``scan_xV`` re-applies the ×V cell-volume scaling
        # (S = V·eff_denom) to recover the transmission ``a = m(1+k)²/S − k``
        # (source-independent) and the blend weight ``w = 1/(1+k)``.  Algebra ONCE.
        V_full = V[:, None, :]                            # (N, 1, nx)
        m = abs_mu[:, None, None] * A_down[:, None, :]    # |μ|·A_down  (N, 1, nx)
        g = m / V_full                                    # ÷V streaming (N, 1, nx)
        cf = d1_closed_form(g, reaction_xs, self.theta)
        return cf.scan_xV(V_full)

    def moment_scan_closure(
        self,
        *,
        abs_mu: np.ndarray,    # (ng, nx) or (nx,)   |μ_n|·A_down/V is rebuilt below
        A_down: np.ndarray,
        V: np.ndarray,
        reaction_xs: np.ndarray,
    ) -> D1ClosedForm:
        r"""The per-cell :class:`~orpheus.transport.spatial._ubld.D1ClosedForm` for the
        1-D moment SCAN (the production LD source-iteration sweep).

        The slope-source-aware companion of :meth:`affine_scan_coefficients`: the
        flat-source scan reads only ``(a, inverse_denom, w)`` (source-independent,
        cached), but the FULL LD with the threaded scattering-slope source
        ``Σ_s·φ̂`` needs the slope fold too — the face-chain affine-source
        contribution (:meth:`~orpheus.transport.spatial._ubld.D1ClosedForm.scan_slope_face_source`)
        and the per-cell ``(ψ̄, ψ̂)`` reconstruction
        (:meth:`~orpheus.transport.spatial._ubld.D1ClosedForm.scan_reconstruct`).  Both
        ride the SAME ``D1ClosedForm`` the flat scan's ``(a, inverse_denom, w)``
        come from (Pattern 2 — ONE LD algebra site; the slope fold is shared with
        the per-cell Schur via ``_slope_fold``).  The slope-augmented affine
        source is
        ``docs/theory/methods/sn/cartesian_multid.rst §ld-ubld-moment-scan``.

        Inputs are broadcastable per-ordinate-per-cell-per-group arrays
        (``abs_mu``/``A_down``/``V`` carry NO group axis; ``reaction_xs`` is ``(ng, …)``);
        the ÷V streaming ``g = |μ|·A_down/V`` is rebuilt here from the SAME
        geometry the cached coefficients use, so the scan's flat ``b`` and its
        slope correction agree to FP.  Slab/Cartesian only (the curvilinear LD
        moment scan is not yet implemented (#158) — guarded at the cache build by
        :meth:`affine_scan_coefficients`; this method is reached only after that
        guard has passed).
        """
        g = abs_mu * A_down / V
        return d1_closed_form(g, reaction_xs, self.theta)
