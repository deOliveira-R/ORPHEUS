r"""Linear-Discontinuous (LD) cell-update strategy — slab / Cartesian (Issue #158).

The first higher-order, O(h²) occupant of the swappable per-cell
spatial-strategy contract (:mod:`orpheus.sn.spatial.cell_update`), sibling to
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
* ``is_affine_scannable = False`` — LD couples two face moments, so it does
  **not** admit the single-upstream affine recurrence that
  ``CumprodScan`` / ``ScanMarch`` consume; it routes through the per-cell DAG
  walk (the ``OneDimPerCellWalk`` loss-representation — Issue #158 Step 2).

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

* :class:`~orpheus.sn.spatial.cell_update.CellUpdate` — the Protocol this
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

from .cell_update import (
    CellResult,
    CellUpdateBase,
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
    _mu: float                       #: ``|μ|`` (slope reconstruction)
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


class LinearDiscontinuous(CellUpdateBase, key="linear_discontinuous"):
    r"""Linear-Discontinuous (LD) cell-update strategy — slab / Cartesian.

    Two spatial moments per (group, ordinate) cell: the cell-average flux and
    the cell-average slope, with the upwind-discontinuous face closure.  The
    slope is eliminated locally by the Schur complement so the update fits the
    scalar :class:`~orpheus.sn.spatial.cell_update.CellUpdate` contract.  See
    the module docstring for the derivation, the SymPy-verified 2×2, and the
    Schur-complement scalar form.

    Not a frozen dataclass-with-fields like DD: LD carries a single immutable
    parameter :attr:`theta` (the slope-moment weight) as a class constant.
    """

    is_linear: ClassVar[bool] = True
    is_positivity_preserving: ClassVar[bool] = False
    is_affine_scannable: ClassVar[bool] = False

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

        # Slab factories always populate abs_mu / volume; narrow explicitly
        # (-O-safe, not ``assert``) so a malformed StreamingTerms fails loudly
        # rather than tripping an opaque ``float·None`` arithmetic error.
        if st.abs_mu is None or st.volume is None:
            raise ValueError(
                "LinearDiscontinuous requires populated abs_mu and volume on "
                f"the slab StreamingTerms; got abs_mu={st.abs_mu}, "
                f"volume={st.volume}."
            )
        mu: float = st.abs_mu               # |μ|, sweep-pre-resolved (slab A=1)
        h: float = st.volume                # slab cell width (V = Δx)
        theta = self.theta
        psi_in = upstream_state.spatial_upstream

        d2p = total_xs * theta * h + mu                      # (ng,)  D₂'
        eff_denom = (mu * mu + (total_xs * h + mu) * d2p) / d2p           # S
        eff_source = s_bar - s_hat * mu * theta / d2p
        eff_numer_upstream = mu * psi_in * (d2p + mu) / d2p
        return _LDCellTerms(
            eff_denom=eff_denom,
            eff_source=eff_source,
            eff_numer_upstream=eff_numer_upstream,
            _mu=mu,
            _d2p=d2p,
            _slope_source=theta * s_hat,
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
        :meth:`CellUpdate.residual`).  Linear in ``cell_avg``; affine in
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
    # the SAME LD, pinned consistent **at the flat Q̂=0 slice** by the
    # group1≡group2 gate (mirrors how DD keeps its per-cell ``residual`` and
    # ``residual_kernel_batch`` separate).  The slope-source (Q̂≠0) equivalence
    # — which exercises the LM-1989 §1.4/§6 slope-row SIGN TRAP — is pinned when
    # the moment source is threaded (Increment C); until then only the slope
    # UNKNOWN path (Q̂=0) is cross-checked between the two forms.
    #
    # 1-D only for now: a multi-D LD is bilinear (an independent slope per
    # axis) — deferred (#158 Increment D).  The slope source ``Q̂`` is folded
    # only when a moment source is threaded (#158 Increment C); the slope
    # UNKNOWN ``ψ̂`` is ALWAYS solved (that is what delivers O(h²)).

    def _kernel_terms(
        self,
        *,
        psi_in: tuple[np.ndarray, ...],
        s_axes: tuple[np.ndarray, ...],
        sigt_cells: np.ndarray,
        Q_cells: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        r"""The ÷V Schur intermediates for the DAG kernel (1-D, flat source).

        Single source for BOTH :meth:`cell_kernel_batch` (solve) and
        :meth:`residual_kernel_batch` (apply).  Returns
        ``(eff_denom S, rhs, d2, g_over_theta, in0)`` — all broadcast to
        ``(N_oct, ng, n_diag)`` (or the inputs' broadcast).
        """
        if len(s_axes) != 1:
            raise NotImplementedError(
                "LinearDiscontinuous.cell_kernel_batch supports d=1 (slab/1-D) "
                f"only; got d={len(s_axes)} streaming axes.  A multi-D LD is "
                "bilinear (an independent slope per axis) and is deferred "
                "(#158 Increment D)."
            )
        g = 0.5 * s_axes[0]                         # |μ|/Δ  (N_oct, 1, n_diag)
        g_over_theta = g / self.theta
        in0 = psi_in[0]                             # (N_oct, ng, n_diag)
        d2 = g_over_theta + sigt_cells             # D₂ (Schur slope denom)
        eff_denom = (g + sigt_cells) + g * g_over_theta / d2          # Schur S
        rhs = Q_cells + g * in0 + g * g_over_theta * in0 / d2  # flat (Q̂ = 0)
        return eff_denom, rhs, d2, g_over_theta, in0

    def cell_kernel_batch(
        self,
        *,
        psi_in: tuple[np.ndarray, ...],
        s_axes: tuple[np.ndarray, ...],
        sigt_cells: np.ndarray,
        Q_cells: np.ndarray,
    ) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
        r"""Pure batched LD cell SOLVE — the DAG wavefront kernel (1-D, slab).

        Solves the Schur-reduced LD cell for the cell-average flux on one
        anti-hyperplane level (vectorised over ``(N_oct, ng, n_diag)``) and
        reconstructs the downstream face :math:`\psi_{\rm out}=\bar\psi+\hat\psi`.
        STORAGE-FREE by contract: the WALK gathers ``psi_in`` and scatters
        ``psi_out``.  The SOLVE arm of the ``_CellSolve`` level operation
        consumed by :meth:`SweepDependencyGraph.walk_full` /
        :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.walk_windowed`.
        """
        eff_denom, rhs, d2, g_over_theta, in0 = self._kernel_terms(
            psi_in=psi_in, s_axes=s_axes, sigt_cells=sigt_cells, Q_cells=Q_cells,
        )
        psi_avg = rhs / eff_denom
        psi_out = psi_avg + g_over_theta * (psi_avg - in0) / d2
        return psi_avg, (psi_out,)

    def residual_kernel_batch(
        self,
        *,
        psi_bar: np.ndarray,
        psi_in: tuple[np.ndarray, ...],
        s_axes: tuple[np.ndarray, ...],
        sigt_cells: np.ndarray,
        Q_cells: np.ndarray,
    ) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
        r"""Pure batched LD operator residual — the apply twin (1-D, slab).

        Evaluates :math:`r = S\,\bar\psi - \mathrm{rhs}` at the PROBE
        cell-average and reconstructs the outgoing face from the probe.  At the
        ``psi_bar`` :meth:`cell_kernel_batch` solves for (same inputs) the
        residual vanishes to FP noise (the batched round-trip).  The APPLY arm
        of the ``_CellResidual`` level operation (the matvec walk); for the
        operator action ``Q_cells`` is zero (source-free).
        """
        eff_denom, rhs, d2, g_over_theta, in0 = self._kernel_terms(
            psi_in=psi_in, s_axes=s_axes, sigt_cells=sigt_cells, Q_cells=Q_cells,
        )
        residual = eff_denom * psi_bar - rhs
        psi_out = psi_bar + g_over_theta * (psi_bar - in0) / d2
        return residual, (psi_out,)
