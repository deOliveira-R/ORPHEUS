r"""LinearOperator promotions of the SN per-cell + per-level primitives.

Phase G Step 1 of Issue #196 — **pure type-system promotion**.  This
module wraps two existing SN spatial primitives as
:class:`~orpheus.numerics.operator.LinearOperator` subclasses so the
operator algebra (``+``, ``-``, ``@``, ``.H``, capability sets, dunder
composition) becomes available at the per-cell and per-level layers
where Phase G's "four operators (L, C, S, F)" composition will fuse
streaming + collision + scattering + fission in subsequent steps.

What ships here
===============

* :class:`SNCellOperator` — wraps
  :class:`~orpheus.sn.spatial.diamond.DiamondDifference`.
  ``solve(visit, total_xs, source, upstream_state) → CellResult``
  delegates to the existing ``DiamondDifference.update`` body
  (bit-identical contract).  ``apply(cell_avg, *, visit, total_xs,
  upstream_state, source) → residual`` computes the per-cell
  discretised operator residual ``L_cell · ψ̄ − q``.  Declares
  ``capabilities = frozenset({CAP_APPLY, CAP_SOLVE})``.

* :class:`AngularRedistribution` — wraps
  :class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`.
  ``apply(psi_cells, *, alpha_half, redist_dAw, tau_mm, volume,
  level_indices, carlson_context) → R`` is the M-M angular recurrence
  transform; the Carlson seed composes as an inner
  :class:`~orpheus.sn.spatial.psi_half_angle_seed.PsiHalfAngleSeed`
  strategy bound on the underlying :class:`MorelMontryAngularSweep`
  instance.  Declares ``capabilities = frozenset({CAP_APPLY})`` at
  Step 1 — back-substitution (``CAP_SOLVE``) is deferred to a future
  step per Pattern 4 (advertise only what works).

What does NOT ship here (deferred to Step 2)
============================================

The call-site unification that closes ERR-026 manifestation #7
(residual O(h) SI-vs-Krylov WDD asymmetry on heterogeneous MR) lives
at Step 2 (``transport_operator_matvec_*`` and ``_sweep_1d_*`` both
routed through ONE :class:`SNCellOperator` instance).  Step 1 is
**wrapper-only** — no behavioural change, all 11 regression snapshots
stay bit-identical under ``rtol=1e-12``.

References
==========

* Issue #196 Phase G plan,
  ``.claude/plans/issue_196_phase_g_four_operator_unification.md``.
* Four-operator architecture reconciliation,
  ``.claude/agent-memory/general-purpose/phase_g_four_operator_architecture_reconciliation.md``
  §2.1 (SNCellOperator verdict), §2.3 (AngularRedistribution verdict).
* Step 1 verification gates,
  ``.claude/agent-memory/test-architect/issue_196_phase_g_step1_verification_gates.md``.
* Phase F closeout (failure-mode profile this defends against),
  ``.claude/agent-memory/method-implementer/issue_168_phase_f_closeout.md``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, ClassVar

import numpy as np

from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_SOLVE,
    LinearOperatorMixin,
)

from .cell_update import CellResult, CellVisit, UpstreamState
from .diamond import DiamondDifference
from .pole_angular_closure import MorelMontryAngularSweep
from .psi_half_angle_seed import CarlsonSweepContext

if TYPE_CHECKING:  # pragma: no cover
    pass


# ═══════════════════════════════════════════════════════════════════════
# SNCellOperator — LinearOperator wrapping DiamondDifference
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class SNCellOperator(LinearOperatorMixin):
    r"""Per-cell SN discrete operator promoted to :class:`LinearOperator`.

    Wraps :class:`~orpheus.sn.spatial.diamond.DiamondDifference` so the
    per-cell update can participate in the project-wide operator algebra
    (``+``, ``-``, ``@``, ``.H``, capability sets).  Step 1 of Phase G
    of Issue #196 — **pure type-system promotion**, no algorithmic
    change, bit-identical to the existing
    :meth:`DiamondDifference.update` on every parameter combination.

    Two capabilities are advertised:

    * :pydata:`~orpheus.numerics.operator.CAP_SOLVE` — the per-cell
      inverse :math:`L_{\text{cell}}^{-1} q`, delegating to the
      existing :meth:`DiamondDifference.update` body to compute the
      cell-average flux + downstream face/angular states.  Bit-identity
      with ``DiamondDifference.update`` is the Step 1 contract.
    * :pydata:`~orpheus.numerics.operator.CAP_APPLY` — the per-cell
      forward action :math:`L_{\text{cell}} \bar\psi - q`, the residual
      of the discretised per-cell balance.  Together with
      :pydata:`CAP_SOLVE` this satisfies the Protocol round-trip
      identity :math:`L_{\text{cell}}\,(L_{\text{cell}}^{-1} q) = q`
      at ``rtol=1e-12`` on non-trivial sources (Gate 2 of Step 1).

    The residual is computed from the per-cell discretised balance
    equation.  Reconstruct the outgoing face flux via the DD spatial
    closure :math:`\psi^s_{\text{out}} = 2\bar\psi - \psi^s_{\text{in}}`
    (and the M-M angular closure
    :math:`\psi_{n+1/2} = (\bar\psi - (1-\tau)\psi_{n-1/2})/\tau`
    for curvilinear); evaluate streaming + collision on
    :math:`\bar\psi`; subtract source.  At the converged cell average
    :math:`\bar\psi^* = L_{\text{cell}}^{-1} q` the residual is exactly
    zero (up to FP rounding).

    Capability rationale
    --------------------

    * No :pydata:`CAP_APPLY_TRANSPOSE` at Step 1.  The adjoint
      (`.H` propagation) is Step 5 scope; advertising a capability
      the operator cannot deliver violates Pattern 4
      (illegal-states-unrepresentable).
    * The frozen + slotted dataclass keeps the operator immutable and
      lightweight.  A single instance can be reused across every
      cell of every sweep without per-call allocation (matches the
      existing :class:`DiamondDifference` allocation pattern).

    See also
    --------
    * :class:`~orpheus.sn.spatial.diamond.DiamondDifference` — the
      wrapped per-cell strategy.
    * :class:`AngularRedistribution` — companion promotion of the
      per-level M-M angular recurrence.
    """

    capabilities: ClassVar[frozenset[str]] = frozenset({CAP_APPLY, CAP_SOLVE})

    # Single underlying DD strategy.  Frozen field with a default
    # factory so callers can simply write ``SNCellOperator()`` and
    # get a default Diamond-Difference cell operator.  An advanced
    # caller may pass a different :class:`CellUpdate` strategy
    # (e.g. a future ``LinearDiscontinuous`` or ``Step``) — but
    # Step 1's bit-identity contract is wired against the canonical
    # DD strategy specifically.
    cell_update: DiamondDifference = field(default_factory=DiamondDifference)

    def solve(
        self,
        visit: CellVisit,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> CellResult:
        r"""Per-cell solve: :math:`L_{\text{cell}}^{-1} q \to \bar\psi`.

        Delegates **verbatim** to :meth:`DiamondDifference.update`.
        Bit-identical contract (Gate 1 of Step 1): the returned
        :class:`CellResult` matches ``DiamondDifference().update(...)``
        on every parameter combination (``np.array_equal`` on every
        field).

        Parameters
        ----------
        visit
            One visit to one cell during an SN sweep.  Carries the
            pre-resolved streaming-direction view of the cell.
        total_xs
            Per-group total cross section :math:`\Sigma_t` on this
            cell, shape ``(ng,)``.
        source
            Per-group volumetric source on this cell, **already
            weight-normalised** by the sweep, shape ``(ng,)``.
            See :class:`DiamondDifference` module docstring for the
            ``Q · V · weight_norm`` convention.
        upstream_state
            Per-cell input state — spatial upstream face flux and
            (curvilinear only) angular upstream half-flux.

        Returns
        -------
        CellResult
            Cell-average flux + downstream spatial / angular states.
        """
        return self.cell_update.update(visit, total_xs, source, upstream_state)

    def apply(
        self,
        cell_avg: np.ndarray,
        *,
        visit: CellVisit,
        total_xs: np.ndarray,
        upstream_state: UpstreamState,
        source: np.ndarray,
    ) -> np.ndarray:
        r"""Per-cell forward action: residual :math:`L_{\text{cell}} \bar\psi - q`.

        Computes the per-cell discretised balance equation residual:

        .. math::

           r \;=\; \mathrm{denom} \cdot \bar\psi
                   \;-\;
                   \bigl[\,\text{upstream-streaming-contributions}\,\bigr]
                   \;-\; q .

        At the converged cell average :math:`\bar\psi^* =
        L_{\text{cell}}^{-1} q` the residual is exactly zero (up to
        FP rounding) — this IS the Protocol round-trip identity
        :meth:`apply` ∘ :meth:`solve` = id (Gate 2 of Step 1).

        Algebraic form (matches :class:`DiamondDifference`'s three
        branches)
        --------------------------------------------------------

        * **Slab**.  The balance equation
          :math:`(2|\mu| + \Delta x \cdot \Sigma_t)\,\bar\psi
          = 2|\mu|\,\psi^s_{\text{in}} + \text{source}` is derived
          from the DD recurrence + closure:

          .. math::

             2\bar\psi = \psi^s_{\text{in}} + \psi^s_{\text{out}}
             = (1+a)\,\psi^s_{\text{in}} + s
             = \bigl[1 + (2|\mu| - \Delta x \Sigma_t)/\text{denom}\bigr]
               \psi^s_{\text{in}} + 2\,\text{source}/\text{denom}

          which multiplied through by ``denom`` and divided by 2 gives
          the balance above.  The residual is
          ``2|μ|·(cell_avg − ψ_in) + chord·Σ_t·cell_avg − source``.

          Note: ``source`` here is the **already-weight-normalised**
          per-contract value (``Q · V · weight_norm`` = ``Q · Δx ·
          weight_norm`` for slab) — same convention :meth:`solve`
          consumes.

        * **Curvilinear (sphere / non-degenerate cylinder)**.  The
          balance is ``denom · cell_avg = source +
          abs_mu·(A_in + A_out)·ψ_spat_in + dA_w · c_in · ψ_angle_in``
          where ``denom = 2|μ|·A_downstream + dA_w·c_out + Σ_t·V``.
          Residual is ``denom·cell_avg − [source + ... upstream
          terms ...]``.

        * **Cylindrical pure-azimuthal degenerate** (``abs_mu < 1e-15``).
          No radial face flow.  Balance is ``denom · cell_avg = source
          + dA_w · c_in · ψ_angle_in`` with ``denom = dA_w · c_out
          + Σ_t · V``.

        Parameters
        ----------
        cell_avg
            Per-group cell average flux :math:`\bar\psi`, shape
            ``(ng,)``.  The probe point at which to evaluate the
            balance residual.
        visit
            Cell-visit packet for this cell.
        total_xs
            Per-group total cross section, shape ``(ng,)``.
        upstream_state
            Per-cell upstream face / angular state — same packet
            consumed by :meth:`solve`.
        source
            Per-group volumetric source on this cell, **already
            weight-normalised**, shape ``(ng,)``.

        Returns
        -------
        residual
            Per-group residual of the cell balance, shape ``(ng,)``.
            Zero (up to FP rounding) at the converged cell average.
        """
        st = visit.streaming_terms

        if st.alpha_in is None:
            # Slab branch.  Balance equation:
            #   (2|μ| + chord·Σ_t)·cell_avg = 2|μ|·ψ_in + source
            # (derived from the DD recurrence + closure; see docstring).
            # Residual = LHS − RHS:
            #   2|μ|·(cell_avg − ψ_in) + chord·Σ_t·cell_avg − source.
            assert st.abs_mu is not None
            abs_mu = st.abs_mu
            chord = st.chord_length
            psi_in = upstream_state.spatial_upstream
            return (
                2.0 * abs_mu * (cell_avg - psi_in)
                + chord * total_xs * cell_avg
                - source
            )

        assert st.abs_mu is not None  # populated by curvilinear factory
        if st.abs_mu < 1e-15:
            # Cylindrical pure-azimuthal degenerate: no radial face flow.
            #   denom = dA_w · c_out + Σ_t · V
            #   balance: denom · cell_avg = source + dA_w · c_in · ψ_angle_in
            return _apply_cylindrical_degenerate_residual(
                st, cell_avg, total_xs, upstream_state, source,
            )

        # Curvilinear non-degenerate branch (sphere / non-degenerate cylinder).
        assert visit.face_area_downstream is not None, (
            "Curvilinear non-degenerate cell apply requires a resolved "
            "face_area_downstream on the CellVisit packet."
        )
        return _apply_curvilinear_residual(
            st,
            visit.face_area_downstream,
            cell_avg,
            total_xs,
            upstream_state,
            source,
        )


# ─── per-branch residual helpers (private) ─────────────────────────────


def _apply_curvilinear_residual(
    st,
    A_downstream: float,
    cell_avg: np.ndarray,
    total_xs: np.ndarray,
    upstream_state: UpstreamState,
    source: np.ndarray,
) -> np.ndarray:
    r"""Residual of curvilinear non-degenerate per-cell balance.

    Mirrors the ``denom / numer`` build order in
    :meth:`DiamondDifference._update_curvilinear` so the algebraic
    expressions are the operation-order twins of the wrapped solve.
    """
    assert st.face_area_inner is not None
    assert st.face_area_outer is not None
    assert st.delta_A_over_w is not None
    assert st.alpha_in is not None
    assert st.alpha_out is not None
    assert st.tau_mm is not None
    assert st.volume is not None
    assert upstream_state.angular_upstream is not None, (
        "Curvilinear cell apply requires an upstream angular state."
    )

    abs_mu = st.abs_mu
    A_inner = st.face_area_inner
    A_outer = st.face_area_outer
    A_total = A_inner + A_outer
    dA_w = st.delta_A_over_w
    alpha_in = st.alpha_in
    alpha_out = st.alpha_out
    tau = st.tau_mm
    V = st.volume

    c_out = alpha_out / tau
    c_in = (1.0 - tau) / tau * alpha_out + alpha_in

    psi_spat_in = upstream_state.spatial_upstream
    psi_angle_in = upstream_state.angular_upstream

    denom = 2.0 * abs_mu * A_downstream + dA_w * c_out + total_xs * V
    numer_upstream = (
        abs_mu * A_total * psi_spat_in
        + dA_w * c_in * psi_angle_in
    )
    # residual = L_cell · cell_avg - q
    #          = denom·cell_avg − [source + numer_upstream]
    return denom * cell_avg - (source + numer_upstream)


def _apply_cylindrical_degenerate_residual(
    st,
    cell_avg: np.ndarray,
    total_xs: np.ndarray,
    upstream_state: UpstreamState,
    source: np.ndarray,
) -> np.ndarray:
    r"""Residual of cylindrical pure-azimuthal degenerate per-cell balance.

    Mirrors :meth:`DiamondDifference._update_cylindrical_degenerate`.
    No radial face-flow term.
    """
    assert st.delta_A_over_w is not None
    assert st.alpha_in is not None
    assert st.alpha_out is not None
    assert st.tau_mm is not None
    assert st.volume is not None
    assert upstream_state.angular_upstream is not None, (
        "Cylindrical-degenerate cell apply requires an upstream "
        "angular state."
    )

    dA_w = st.delta_A_over_w
    alpha_in = st.alpha_in
    alpha_out = st.alpha_out
    tau = st.tau_mm
    V = st.volume

    c_out = alpha_out / tau
    c_in = (1.0 - tau) / tau * alpha_out + alpha_in

    psi_angle_in = upstream_state.angular_upstream

    denom = dA_w * c_out + total_xs * V
    numer_upstream = dA_w * c_in * psi_angle_in
    return denom * cell_avg - (source + numer_upstream)


# ═══════════════════════════════════════════════════════════════════════
# AngularRedistribution — LinearOperator wrapping MorelMontryAngularSweep
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class AngularRedistribution(LinearOperatorMixin):
    r"""Curvilinear angular redistribution operator promoted to :class:`LinearOperator`.

    Wraps :class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`
    so the per-level M-M angular recurrence can participate in the
    operator algebra.  The Carlson half-angle seed
    (:class:`~orpheus.sn.spatial.psi_half_angle_seed.CarlsonInwardSweep`
    by default, supplied via ``carlson_context`` at :meth:`apply` time)
    composes as an inner strategy bound on the
    :class:`MorelMontryAngularSweep` instance.

    Step 1 of Phase G of Issue #196 — pure type-system promotion, no
    algorithmic change.  All 11 regression snapshots stay bit-identical
    under ``rtol=1e-12``.

    Capability set
    --------------

    Only :pydata:`~orpheus.numerics.operator.CAP_APPLY` is declared.

    * **No** :pydata:`~orpheus.numerics.operator.CAP_SOLVE`.  The M-M
      recurrence is upper-triangular in ordinate index, so its
      back-substitution exists algebraically — but Step 1 ships no
      working back-substitution code, and advertising a capability
      the operator cannot deliver violates Pattern 4
      (illegal-states-unrepresentable).  A future step may add
      ``CAP_SOLVE`` once back-substitution is implemented and
      verified.
    * No :pydata:`CAP_APPLY_TRANSPOSE`.  Adjoint is Step 5 scope.

    Linear by construction
    ----------------------

    The M-M recurrence and the Carlson seed (per
    :class:`~orpheus.sn.spatial.psi_half_angle_seed.CarlsonInwardSweep`)
    are both linear in the input ψ.  The 57 function-level tests at
    :file:`tests/sn/spatial/test_sweep_vs_apply_consistency.py` pin
    the linearity of the Carlson seed; the M-M recurrence's affine
    coefficients (``α``, ``ΔA/w``, ``τ``) are constants of the
    discretisation, so the redistribution output is linear in
    ``psi_cells``.  Operator-level linearity is verified at
    :file:`tests/sn/spatial/test_angular_redistribution.py` as Gate 6
    #5 of Step 1.

    See also
    --------
    * :class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`
      — the wrapped per-level recurrence strategy.
    * :class:`~orpheus.sn.spatial.psi_half_angle_seed.CarlsonInwardSweep`
      — canonical Hébert §3.9.4 (3.432)-(3.435) inward μ = −1 sweep
      bound as the M-M seed by default.
    * :class:`SNCellOperator` — companion per-cell promotion.
    """

    capabilities: ClassVar[frozenset[str]] = frozenset({CAP_APPLY})

    # The underlying M-M strategy.  Default factory constructs the
    # Phase D canonical strategy with the Carlson coupled-pole seed.
    angular_sweep: MorelMontryAngularSweep = field(
        default_factory=MorelMontryAngularSweep,
    )

    def apply(
        self,
        psi_cells: np.ndarray,
        *,
        alpha_half: "np.ndarray | list[np.ndarray]",
        redist_dAw: "np.ndarray | list[np.ndarray]",
        tau_mm: "np.ndarray | list[np.ndarray]",
        volume: np.ndarray,
        level_indices: "list[np.ndarray] | None" = None,
        carlson_context: "CarlsonSweepContext | list[CarlsonSweepContext] | None" = None,
    ) -> np.ndarray:
        r"""Apply the M-M angular redistribution to a per-cell ψ.

        Delegates **verbatim** to
        :meth:`MorelMontryAngularSweep.__call__`.  The Carlson seed is
        consumed via ``carlson_context``; when ``None`` the recurrence
        falls back to the Phase B hardcoded zero seed (the ERR-026
        anti-pattern; preserved only for regression-safety ablations).

        Parameters
        ----------
        psi_cells
            Cell-centred angular flux, shape ``(ng, M, nx)``.  The
            recurrence reads the per-ordinate slice ``psi_cells[:, m, :]``
            at each step of the M-M weighted DD chain.
        alpha_half
            :math:`\alpha_{n+1/2}` dome values.  Spherical: ``ndarray``
            of shape ``(M+1,)``.  Cylindrical: ``list`` of per-level
            arrays.
        redist_dAw
            Per-cell, per-ordinate :math:`\Delta A / w` geometry
            factor.  Spherical: ``ndarray`` of shape ``(nx, M)``.
            Cylindrical: ``list`` of per-level arrays.
        tau_mm
            Morel--Montry :math:`\tau` clamp, same per-geometry shape
            convention as ``alpha_half`` / ``redist_dAw``.
        volume
            Per-cell volume, shape ``(nx,)``.
        level_indices
            Spherical: ``None``.  Cylindrical: ``list`` of per-level
            ordinate index arrays.
        carlson_context
            Carlson inward sweep inputs (Σ_t, Δr, μ_quad, weights,
            bc_outer_value).  Spherical: single
            :class:`CarlsonSweepContext`.  Cylindrical: ``list`` of
            per-level contexts.  When ``None``, M-M falls back to the
            Phase B hardcoded zero seed.

        Returns
        -------
        redist
            Redistribution term, shape ``(ng, M, nx)``.  Same shape
            as the input ``psi_cells``.

        Notes
        -----
        Linearity in ``psi_cells`` holds when the Carlson seed is
        linear in ``psi_cells`` (it is — both
        :class:`CarlsonInwardSweep` and :class:`ZeroSeed` are linear)
        and the M-M recurrence's coefficients (``alpha_half``,
        ``redist_dAw``, ``tau_mm``) are held fixed across the linear
        combination.  See the operator algebra contract at
        :class:`~orpheus.numerics.operator.LinearOperator`.
        """
        return self.angular_sweep(
            psi_cells=psi_cells,
            alpha_half=alpha_half,
            redist_dAw=redist_dAw,
            tau_mm=tau_mm,
            volume=volume,
            level_indices=level_indices,
            carlson_context=carlson_context,
        )


__all__ = [
    "AngularRedistribution",
    "SNCellOperator",
]
