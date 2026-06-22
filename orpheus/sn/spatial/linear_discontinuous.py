r"""Linear-Discontinuous (LD) cell-update strategy — slab / Cartesian (Issue #158).

The first higher-order, O(h²) occupant of the swappable per-cell
spatial-strategy contract (:mod:`orpheus.sn.spatial.scheme`), sibling to
the shipping :class:`~orpheus.sn.spatial.diamond.DiamondDifference` (DD).
**Full** LD (with the slope SOURCE moment :math:`\hat Q` carried) is
diffusion-limit-consistent — Larsen-Morel-Miller 1987 proves Step has **no**
valid intermediate diffusion limit while LD has all four, the load-bearing
reason a non-DD spatial scheme can lift the curvilinear pole-cell ``O(h)`` floor
(Issue #233).

.. warning::

   **Diffusion-limit status — read before using LD on a diffusive problem.**
   The CURRENT implementation (Issue #158 Increment A) ships LD with a **flat
   cell source** (:math:`\hat Q = 0`): the slope UNKNOWN :math:`\hat\psi` is
   always solved (that is what delivers O(h²) on smooth, streaming-dominated
   problems), but the slope SOURCE — including the scattering-source slope
   :math:`\Sigma_s\hat\phi` — is **not yet threaded** (that needs the flux slope
   :math:`\hat\phi` in the iterate, the "global moment-contract", Increment C).
   The diffusion-limit proof above is a property of the FULL (canonical) LD; the
   flat cut is **O(h²) but NOT diffusion-limit-consistent** — on an optically
   thick, scattering-dominated (``c → 1``) mesh the flat-source flux can be tens
   of percent off (it recovers only as cells are refined optically thin).  Use
   DD or a thin mesh for thick diffusive problems until Increment C lands; the
   forward tripwire is
   ``tests/sn/verification/mms/test_mms_ld_slab.py::test_ld_thick_diffusive_limit_xfail``.

Why LD carries two moments
==========================

DD has a single per-cell unknown — the cell-average flux :math:`\bar\psi`.  LD
represents the in-cell angular flux as a **linear function** (Larsen & Morel
1989, JCP 83(1):212-236, Eqs. 4.1a-c), so each (group, ordinate) cell carries
**two** spatial moments:

.. math::

    \bar\psi_{j} = \tfrac1{h_j}\!\int \psi(x)\,dx,
    \qquad
    \hat\psi_{j} = \tfrac{6}{h_j^2}\!\int (x-x_j)\,\psi(x)\,dx,
    \qquad
    \psi(x) \approx \bar\psi_j + \tfrac{2}{h_j}(x-x_j)\,\hat\psi_j .

The closure is **upwind/discontinuous** at faces (LM-1989 Eq. 4.3c): the
downstream face value :math:`\psi_{\rm out} = \bar\psi + \hat\psi` (after
ORPHEUS pre-resolves the sweep direction, so the strategy sees only the
"downstream" face) feeds the next cell; the cell's own *inflow*-face
reconstruction is discarded in favour of the upwind neighbour's outflow (no
continuity enforced — that is the "discontinuous" in LD).

The 2×2 cell system (slab, :math:`A=1`, :math:`V=h`, :math:`|\mu|`
sweep-pre-resolved, ``source`` = the full per-ordinate RHS moments) is, with
the upwind closure :math:`\psi_{\rm out}=\bar\psi+\hat\psi` substituted into
the balance + slope moment equations (LM-1989 Eqs. 4.3a-b, :math:`\theta=1/3`):

.. math::

    \begin{bmatrix} \Sigma_t h + |\mu| & |\mu| \\[2pt]
                    -\,|\mu|/\theta    & \Sigma_t h + |\mu|/\theta \end{bmatrix}
    \begin{bmatrix} \bar\psi \\ \hat\psi \end{bmatrix}
    =
    \begin{bmatrix} \bar{Q} h + |\mu|\,\psi_{\rm in} \\[2pt]
                    \hat{Q} h - (|\mu|/\theta)\,\psi_{\rm in} \end{bmatrix}.

The off-diagonal/RHS signs of the slope row are the documented correctness
trap (LM-1989 memo §1.4/§6 — the published boxed form is internally
inconsistent); the system above was regenerated symbolically with SymPy and
validated against the strongest oracle: **LD is exact on a linear-in-x flux**
(``ψ̄, ψ̂, ψ_out`` recovered to machine precision for any
:math:`\psi=a+bx`).  See the derivation in the theory page and the foundation
tests in ``tests/sn/spatial/test_linear_discontinuous.py``.

The Schur-complement scalar contract
====================================

The per-cell linear system has two unknowns, but the slope :math:`\hat\psi` is
a **local** quantity that can be eliminated by the Schur complement of the 2×2
with respect to the slope row.  This collapses the cell update to a *scalar*
relation in :math:`\bar\psi`,

.. math::

    S\,\bar\psi = \mathrm{eff\_source} + \mathrm{eff\_numer\_upstream},
    \qquad
    D_2' = \Sigma_t\theta h + |\mu|,
    \quad
    S = \frac{|\mu|^2 + (\Sigma_t h + |\mu|)\,D_2'}{D_2'},

with :math:`\mathrm{eff\_source} = \bar{S} - |\mu|\theta\,\hat{S}/D_2'` and
:math:`\mathrm{eff\_numer\_upstream} = |\mu|\,\psi_{\rm in}(D_2'+|\mu|)/D_2'`
(here :math:`\bar S, \hat S` are the prepared average/slope source moments).
This is exactly the shape of DD's ``denom·ψ̄ = source + numer_upstream``, so LD
fits the **existing scalar** :meth:`residual` contract with no contract change
— :math:`\hat\psi` is reconstructed locally inside :meth:`update`.  The slope
*source* :math:`\hat{Q}` is still required (carried on ``source[1]``); in a
**fixed-source** problem it is supplied directly (manufactured), which is why
LD is staged on the fixed-source MMS first.  Threading the slope through a
*source-iteration* global flux iterate (so the within-group scattering source
slope :math:`\propto\Sigma_s\hat\phi` is assembled) is the deferred global
moment-contract — see Issue #158 / the plan.

Scope (slab / Cartesian only — for now)
=======================================

This occupant implements the **slab/Cartesian** LD only.  The curvilinear
(sphere / cylinder) LD cell update — where the Morel-Montry angular weight
:math:`\tau` enters the average-moment denominator and the radial moment
:math:`(r-r_j)` weighting produces a slope-curvature coupling — is **not
published** and must be derived (a SymPy task with the slab-LD and
curvilinear-DD limits as the two reduction oracles).  Until then, :meth:`update`
/ :meth:`residual` raise :exc:`NotImplementedError` on a curvilinear visit
(signalled by ``upstream_state.angular_upstream is not None``, exactly as DD
signals the presence of angular redistribution).

Traits
======

* ``is_linear = True`` — the LD closure (without a positivity fixup) is linear
  in ``source`` and ``upstream_state``.  A positivity fixup (set-to-zero /
  Larsen) would make it non-linear; the first cut ships **without** a fixup
  (``is_positivity_preserving = False``), so the clean 2×2 linear solve stands.
* ``is_positivity_preserving = False`` — LD can produce negative cell-average
  AND negative edge fluxes in thin / steep-source cells (Lewis & Miller §5.3).
* ``is_affine_scannable = True`` (Issue #158 Increment B) — although LD carries
  TWO cell moments, the slope :math:`\hat\psi` is a LOCAL quantity eliminated by
  the Schur complement, leaving a *single-upstream* affine recurrence
  :math:`\psi_{\rm out}=a\,\psi_{\rm in}+b` with ``a`` source-independent (proven
  numerically to machine ε).  So LD ALSO rides the DAG-free scan schedules
  (``CumprodScan`` / ``ScanMarch``) via the coefficient model — it supplies its
  three per-cell coefficients ``(a, inverse_denom, w)`` through
  :meth:`~LinearDiscontinuous.affine_scan_coefficients` and the generic base
  reconstruction staticmethods
  (:meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.source_emission` /
  :meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.cell_average` /
  :meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.outgoing_face_from_average`)
  do the rest.  (The polymorphic
  ``FullFieldWavefront`` DAG oracle — Increment A — still backs the group-2
  kernel; the two-paths gate pins ``CumprodScan``-LD ≡ ``FullFieldWavefront``-LD.)
  **Slab/Cartesian only** — the curvilinear scan coefficients are unpublished
  (the guard in ``affine_scan_coefficients`` raises on non-neutral curvature).

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

* :class:`~orpheus.sn.spatial.scheme.DiscretizationScheme` — the Protocol this
  strategy satisfies; :class:`~orpheus.sn.spatial.diamond.DiamondDifference`
  — the single-moment sibling whose ``update`` / ``residual`` shape this
  mirrors.
* :doc:`/theory/discrete_ordinates`, "Cell update strategies" → "Linear
  Discontinuous".
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

import numpy as np

from ._ubld import (
    AVERAGE_MOMENT,
    D1ClosedForm,
    assemble_inflow_axis,
    assemble_ubld,
    d1_closed_form,
    face_moment_tail,
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
            "unpublished and must be derived (Issue #158, curvilinear arm). "
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
    scalar :class:`~orpheus.sn.spatial.scheme.DiscretizationScheme` contract.  See
    the module docstring for the derivation, the SymPy-verified 2×2, and the
    Schur-complement scalar form.

    Not a frozen dataclass-with-fields like DD: LD carries a single immutable
    parameter :attr:`theta` (the slope-moment weight) as a class constant.
    """

    is_linear: ClassVar[bool] = True
    is_positivity_preserving: ClassVar[bool] = False
    is_affine_scannable: ClassVar[bool] = True
    r"""LD admits the single-upstream affine recurrence (#158 Increment B): the
    slope :math:`\hat\psi` is eliminated by the Schur complement, leaving
    :math:`\psi_{\rm out}=a\psi_{\rm in}+b` with ``a`` source-independent.  LD
    therefore supplies the scan coefficient triple via
    :meth:`affine_scan_coefficients` and rides ``CumprodScan`` (and the **1-D**
    ``ScanMarch``).  Slab/Cartesian only — the method raises on curvilinear
    geometry.  NOTE this is the **single-axis** trait; LD leaves
    ``transverse_coupling_is_facewise`` at its ``False`` default because its
    multi-D closure is BILINEAR (an independent slope moment per axis), so the
    d-D coupling is slope-wise (1st-order), NOT a 0th-order face value — the
    d-D closure does not separate into per-axis scans, and LD rides the DAG
    wavefront in :math:`d \ge 2`, not the scan-march (#38 / #240 D5b)."""

    spatial_basis_per_axis: ClassVar[int] = 2
    r"""LD carries TWO spatial moments per axis — the Legendre basis
    :math:`\{1, P_1\}` (cell-average + slope).  The per-cell unknown is
    therefore :math:`2^d` (d=1: the slab ``{ψ̄, ψ̂}``; d=2: the bilinear
    ``{ψ̄, ψ̂_y, ψ̂_x, ψ̂_xy}`` in tensor-Legendre Kronecker order — slot 1 = the
    y-slope, slot 2 = the x-slope (axis 0 = x outer, axis 1 = y inner); the cross
    moment ``ψ̂_xy`` is the
    diffusion-limit-load-bearing term simplex-P1 lacks), and each downstream
    face carries :math:`2^{d-1}` transverse moments (d=1: a scalar; d=2:
    ``{face-bar, face-slope}``).  The ``per_axis > 1`` gate routes LD onto the
    multi-moment face-cochain + the moment-reducing emit (#240 D5b);
    DiamondDifference at ``per_axis == 1`` keeps the scalar path byte-identical.
    See :mod:`orpheus.sn.spatial._ubld` for the d-generic Kronecker primitive."""

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
        # (#240 D5b-S1 Branch 2): slab A_down = 1, V = h, so the ÷V
        # streaming-over-volume is g = |μ|/h.  ``schur_xV`` returns the ×V
        # per-cell intermediates (S, eff_source, eff_numer, θ·ŝ, |μ|A_down, D₂')
        # in the SAME shape the production closed forms used to compute inline —
        # the LD 2×2 algebra now lives ONCE (in ``_ubld.d1_closed_form``), proven
        # == the d-generic dense primitive's d=1 reduction.
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
    # generic oracle, including 1-D slab — runs LD through these, exactly as it
    # runs DD.  The discretization is polymorphic; the walk owns the storage,
    # the scheme owns only its cell algebra (#158, the basic qualification).
    #
    # Convention: the kernel uses the d-generic Cartesian streaming coefficient
    # ``s = 2|μ|/Δ`` per axis (``A = 1``, ``V = Δ``) — the SAME ``s_axes`` DD's
    # kernel consumes — so ``g = s/2 = |μ|/Δ`` is LD's streaming-over-volume.
    # This is the ÷V form of the per-cell :meth:`update`/:meth:`residual` 2×2
    # (which carry the ×V, ``∂r/∂source = −1`` contract scaling); the two are
    # the SAME LD, pinned consistent by the group1≡group2 gate.
    #
    # UNIFIED MOMENT KERNEL (#240 D5b-S3): a matvec is a forward APPLY — applying
    # the per-cell ``2^d × 2^d`` operator to the moment vector is intrinsically
    # MOMENT-valued in EVERY d.  The kernel pair therefore rides ONE d-generic
    # dense UBLD path (``_ubld_system`` → ``per_cell_solve``) for all d; the
    # former d=1 Schur-reduced SCALAR matvec arm (``_kernel_terms`` →
    # ``kernel_rhs``) was a flat-source artifact (Q̂=0 left the slope ψ̂ with no
    # global coupling) — RETIRED.  The L16 closed-form fast path lives in the
    # d=1 PRODUCTION sweep (CumprodScan → ``affine_scan_coefficients`` →
    # ``d1_closed_form``), selected by strategy on geometry, NOT by an ``if d==1``
    # in the kernel.  The ÷V dense system is SCALE-FREE: dividing the Galerkin
    # balance by the cell volume ``V = ∏_a Δ_a`` leaves a system depending ONLY
    # on the per-axis ÷V streaming ``g_a = |μ_a|/Δ_a`` (the ``s_axes`` the kernel
    # receives) and ``Σ_t``, so the dense assembler is fed unit widths
    # ``hs = [1, …]`` and ``mus = [g_0, …]``.  At d=1 this is EXACTLY
    # ``d1_closed_form``'s 2×2 (the closed form is its analytic Schur, proven ==
    # by ``tests/sn/spatial/test_ld_ubld_primitive.py``), so the moment matvec
    # and the closed-form scan agree at the solved state.

    def _ubld_system(
        self,
        s_axes: tuple[np.ndarray, ...],
        reaction_xs: np.ndarray,
        Q_cells: np.ndarray,
    ) -> tuple[dict, np.ndarray]:
        r"""Assemble the batched ÷V UBLD dense system + the source-moment RHS.

        Returns ``(assembled, R_source)`` where ``assembled`` is the
        :func:`~orpheus.sn.spatial._ubld.assemble_ubld` dict (``A`` of shape
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
        size = 2**d
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
        cochain).  :func:`~orpheus.sn.spatial._ubld.assemble_inflow_axis`
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
        needs_moment_axis = face_moment_tail(2 ** (d - 1)) == ()
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
        tail = face_moment_tail(2 ** (d - 1))
        tensor = psi_moments.reshape(*batch, *([2] * d))
        faces = []
        for a in range(d):
            cell_axis = len(batch) + a
            face = tensor.take(0, axis=cell_axis) + tensor.take(1, axis=cell_axis)
            faces.append(face.reshape(*batch, *tail))
        return tuple(faces)

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
        :func:`~orpheus.sn.spatial._ubld.per_cell_solve`.  Returns the full
        ``2^d``-moment ``psi_avg`` ``(N_oct, ng, n_diag, 2^d)`` — the cell-average
        in slot 0, the slope moments after — and the d-tuple of
        ``2^{d-1}``-moment downstream faces (one per axis).  ONE moment path for
        every d (#240 D5b-S3): at d=1 the system is the LD ``2×2`` whose Schur
        complement is :func:`~orpheus.sn.spatial._ubld.d1_closed_form` (the L16
        production scan's fast path), proven equal by
        ``tests/sn/spatial/test_ld_ubld_primitive.py``.

        STORAGE-FREE by contract: the WALK gathers ``psi_in`` (per-axis
        ``2^{d-1}``-moment faces; a scalar at d=1) and scatters the outgoing
        faces.  The SOLVE arm of the ``_CellSolve`` level operation consumed by
        :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.walk_full` /
        :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.walk_windowed`.
        """
        d = len(s_axes)
        assembled, R_source = self._ubld_system(s_axes, reaction_xs, Q_cells)
        R = R_source + self._ubld_inflow(s_axes, psi_in)
        psi_moments = per_cell_solve(assembled, R)
        return psi_moments, self._ubld_outgoing_faces(psi_moments, d)

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
        fluxes propagate along the wavefront exactly as the sweep's — L21 matvec
        ≡ sweep).  ONE moment path for every d (#240 D5b-S3): the residual is
        intrinsically moment-valued (average ROW + slope ROWS), so the former
        d=1 Schur-reduced SCALAR matvec arm — a flat-source (Q̂=0) artifact — is
        retired.

        At the ``psi_bar`` :meth:`cell_kernel_batch` solves for (same inputs)
        the residual vanishes to FP noise (the batched round-trip, D5b.1).  The
        APPLY arm of the ``_CellResidual`` level operation (the matvec walk).

        Mass-diagonal normalization (#240 D5b-S3 — the matvec/sweep moment-source
        consistency).  The UBLD RHS folds the cell source mass-weighted
        (``R_source = M·S⃗``, the test-function projection); but the operator-algebra
        ``A = (L+C) − S`` subtracts the scattering moment source ``S.apply(ψ)`` RAW
        (per-ordinate, un-projected) at the ``OperatorSum`` level.  Dividing the
        residual by the (diagonal) mass ``M`` puts ``(L+C)ψ`` in RAW per-moment
        units so the row-for-row ``− S.apply(ψ)`` is consistent — the slope rows
        would otherwise disagree by the Legendre mass factor ``M_ii = θ^{|i|}``
        (``|i|`` = number of active slope axes).  The division preserves the
        round-trip (``A·ψ⃗ − R = 0 ⇒ (A·ψ⃗ − R)/M_ii = 0``) and matches the d=1
        DD ÷V convention (where ``M_00 = 1`` so the average row is already raw).
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

    # ── Scan-family coefficients (group 3 — the DAG-free schedules) ──────────
    #
    # The ×V "denom" convention the scan machinery uses (CollisionCache /
    # ordinate_scan): the SAME LD as the ÷V group-2 kernel above, scaled by V
    # (S_scan = V·S_kernel, SymPy-verified), so CumprodScan-LD and
    # FullFieldWavefront-LD are principled-equivalent at nULP (the two-paths
    # gate).  Slab/Cartesian only — the curvilinear LD closure is unpublished.

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

        The curvilinear (sphere/cylinder) LD scan closure is unpublished.  Its
        signal is non-neutral curvature: slab carries ``dA_w == 0`` and
        ``c_out == 0`` EXACTLY (the :func:`slab_streaming` neutral element);
        curvilinear carries non-zero values.  Raising here fails fast at the
        :class:`~orpheus.sn.spatial.sweep_cache.CollisionCache` build
        (``SNSolver.__init__``) before any sweep or matvec runs — so a 1-D
        sphere/cylinder mesh carrying ``LinearDiscontinuous`` (which would match
        ``CumprodScan.supports`` via ``is_1d and is_affine_scannable``) is
        rejected loudly rather than silently running DD-shaped curvature math.
        """
        if np.any(dA_w != 0.0) or np.any(c_out != 0.0):
            raise NotImplementedError(
                "LinearDiscontinuous.affine_scan_coefficients supports the "
                "slab/Cartesian LD only; the curvilinear (sphere/cylinder) LD "
                "scan closure is unpublished (#158).  A non-neutral curvature "
                "was detected (dA_w / c_out are not all zero)."
            )
        # Single-source the LD 2×2 Schur through the shared d=1 closed form
        # (#240 D5b-S1 Branch 2): the ×V scan reads ``(a, inverse_denom, w)``
        # off the same helper the ÷V kernel and ×V per-cell views use.  The
        # ÷V streaming-over-volume is g = |μ|·A_down/V; ``scan_xV`` re-applies
        # the ×V cell-volume scaling (S = V·eff_denom) to recover the
        # transmission ``a = m(1+k)²/S − k`` (source-independent) and the
        # blend weight ``w = 1/(1+k)``.  The LD algebra lives ONCE.
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
        r"""The per-cell :class:`~orpheus.sn.spatial._ubld.D1ClosedForm` for the
        1-D moment SCAN (#240 D5b-S3 OWED-2).

        The slope-source-aware companion of :meth:`affine_scan_coefficients`: the
        flat-source scan reads only ``(a, inverse_denom, w)`` (source-independent,
        cached), but the FULL LD with the threaded scattering-slope source
        ``Σ_s·φ̂`` needs the slope fold too — the face-chain affine-source
        contribution (:meth:`~orpheus.sn.spatial._ubld.D1ClosedForm.scan_slope_face_source`)
        and the per-cell ``(ψ̄, ψ̂)`` reconstruction
        (:meth:`~orpheus.sn.spatial._ubld.D1ClosedForm.scan_reconstruct`).  Both
        ride the SAME ``D1ClosedForm`` the flat scan's ``(a, inverse_denom, w)``
        come from (Pattern 2 — ONE LD algebra site; the slope fold is shared with
        the per-cell Schur :meth:`~orpheus.sn.spatial._ubld.D1ClosedForm.schur_xV`
        via ``_slope_fold``).

        Inputs are broadcastable per-ordinate-per-cell-per-group arrays
        (``abs_mu``/``A_down``/``V`` carry NO group axis; ``reaction_xs`` is ``(ng, …)``);
        the ÷V streaming ``g = |μ|·A_down/V`` is rebuilt here from the SAME
        geometry the cached coefficients use, so the scan's flat ``b`` and its
        slope correction agree to FP.  Slab/Cartesian only (the curvilinear LD
        moment scan is unpublished — guarded at the cache build by
        :meth:`affine_scan_coefficients`; this method is reached only after that
        guard has passed).
        """
        g = abs_mu * A_down / V
        return d1_closed_form(g, reaction_xs, self.theta)
