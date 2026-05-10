r"""Diamond-Difference (DD) cell-update strategy.

This module ships the first concrete
:class:`~orpheus.sn.spatial.cell_update.CellUpdate` strategy of the
SN reshape campaign — Round 2 of Wave C (Issue #158).  The Round 1
contract (:mod:`orpheus.sn.spatial.cell_update`) defined the
:class:`~orpheus.sn.spatial.cell_update.CellUpdate` Protocol and the
:class:`~orpheus.sn.spatial.cell_update.UpstreamState` /
:class:`~orpheus.sn.spatial.cell_update.CellResult` dataclasses but
shipped no concrete strategy; this module fills that hole.

Why this exists (the architectural intent)
==========================================

The historical SN sweep at :mod:`orpheus.sn.sweep` inlines the cell-
update algebra for slab, sphere, and cylinder, plus a special-case
branch for cylindrical pure-azimuthal directions.  Lifting the
closure into a strategy class makes it possible to:

* unit-test the cell-update math in isolation, without spinning up
  a full :class:`~orpheus.sn.geometry.SNMesh` + iteration loop;
* swap closures (DD :math:`\to` Linear Discontinuous, Step,
  Exponential Characteristic) without rewriting the sweep driver;
* keep the per-iteration sweep orchestration thin while the
  closure-specific algebra lives local to a strategy module.

Round 2 ships :class:`DiamondDifference` as a **bit-identical**
extraction of the existing inlined sweep math at
:mod:`orpheus.sn.sweep`.  The bit-identity contract is gated by
``np.array_equal`` hand-calc tests that mirror the sweep's scalar
formulas at the operation level — same computational order, same
float intermediates — so this module's output equals the sweep's
output bit-for-bit on synthetic per-cell inputs.  Wave D (Issue
#159) will rewrite the sweep itself to dispatch through the
:class:`CellUpdate` Protocol; until then, :class:`DiamondDifference`
lives as a parallel, testable abstraction.

Three branches, one strategy
============================

Per Wave C decision **D5** (one geometry-polymorphic class), this
module ships **one** :class:`DiamondDifference` that handles slab,
sphere, and cylinder by branching on two
:class:`~orpheus.geometry.reduced_operator.StreamingTerms` fields:

#. **Slab** (``alpha_in is None``)

   The flat / Cartesian DD recurrence — the per-cell scalar form of
   :eq:`dd-recurrence`.  The sweep's vectorised
   :func:`~orpheus.sn.sweep._sweep_1d_cumprod` path uses cumulative
   products to solve the recurrence over a whole row of cells; here
   we mirror the same algebra at the **single-cell** level so the
   strategy's per-cell output equals one row of the cumprod path's
   per-cell output bit-for-bit.

   The slab denominator is
   :math:`2|\mu| + \Delta x \cdot \Sigma_t` and the recurrence is

   .. math::

      \psi_{\rm out} \;=\; \frac{2|\mu| - \Delta x\,\Sigma_t}
                                   {2|\mu| + \Delta x\,\Sigma_t}\,
                              \psi_{\rm in}
                          \;+\; \frac{2\,Q\,\Delta x / W}
                                     {2|\mu| + \Delta x\,\Sigma_t},
      \qquad
      \overline{\psi}
        \;=\; \tfrac12\bigl(\psi_{\rm in} + \psi_{\rm out}\bigr),

   where ``source = Q · (Δx/W) · weight_norm`` arrives **already
   weight-normalised** by the sweep (``weight_norm = 1/Σ_n w_n``);
   for slab the cell volume is ``V = Δx`` so the contract's
   ``Q · V · weight_norm`` form collapses to
   ``Q · Δx · weight_norm``.  The factor of 2 in the source term
   is the slab DD source coefficient ``2 · weight_norm · Δx /
   denom`` — see :func:`orpheus.sn.sweep._solve_recurrence` lines
   208–222 (the recurrence) and the
   :func:`orpheus.sn.sweep._sweep_1d_cumprod` ``source_coeff`` /
   ``stream_coeff`` definitions at lines 117–123 for the verbatim
   reference.  The implementation here computes ``psi_out`` first
   and **then** averages, in that order, to match the sweep's
   ``0.5 * (psi_in + psi_out)`` operation order at line 222.

   For slab, ``outgoing_angular_state`` is ``None`` — slab geometry
   has no angular redistribution.

#. **Curvilinear** (``alpha_in is not None`` and
   ``abs_mu >= 1e-15``)

   Sphere or cylinder, away from the cylindrical pure-azimuthal
   degenerate case.  Mirrors :func:`orpheus.sn.sweep._sweep_1d_spherical`
   lines 350–361 (and the structurally identical cylindrical inward /
   outward branches at lines 511–531 and 548–575) verbatim:

   .. math::

      c_{\rm out} = \alpha_{n+\tfrac12} / \tau,
      \qquad
      c_{\rm in} = \frac{1 - \tau}{\tau}\,\alpha_{n+\tfrac12}
                   + \alpha_{n-\tfrac12},

   the M-M closure constants combining the Bailey 2009 dome
   :eq:`alpha-recursion` with the Morel–Montry weights
   :eq:`mm-weights`.  Then

   .. math::

      \mathrm{denom} = 2|\mu|\,A_{\rm downstream}
                       + (\Delta A / w)\,c_{\rm out}
                       + \Sigma_t\,V,

   where :math:`A_{\rm downstream}` (= ``visit.face_area_downstream``)
   is the sweep-direction-resolved outgoing face area — equal to
   ``streaming_terms.face_area_outer`` for outward sweeps
   (:math:`\mu \ge 0`) and ``streaming_terms.face_area_inner`` for
   inward sweeps (:math:`\mu < 0`).  Then

   .. math::

      \mathrm{numer} = (Q\,V/W)
                       + |\mu|\,(A_{\rm inner} + A_{\rm outer})\,
                          \psi^s_{\rm in}
                       + (\Delta A / w)\,c_{\rm in}\,\psi_{n-\tfrac12},

   the symmetric sum :math:`A_{\rm inner} + A_{\rm outer}` being
   invariant under sweep direction.

   .. math::

      \overline{\psi} = \mathrm{numer}/\mathrm{denom},

   followed by the WDD spatial closure :eq:`wdd-closure` and the
   M-M angular closure :eq:`mm-weights`:

   .. math::

      \psi^s_{\rm out} = 2\overline{\psi} - \psi^s_{\rm in},
      \qquad
      \psi_{n+\tfrac12} = (\overline{\psi}
                          - (1 - \tau)\,\psi_{n-\tfrac12}) / \tau.

   The contract guarantees ``source = Q · V · weight_norm`` arrives
   already weight-normalised by the sweep — so ``source`` directly
   plays the role of ``QV[i]`` in the sweep's lines 350–355.

#. **Cylindrical pure-azimuthal degenerate**
   (``alpha_in is not None`` and ``abs_mu < 1e-15``)

   When the level's axial direction cosine :math:`|\mu_z| \to 1` the
   radial direction cosine :math:`|\eta| \to 0` and the cell has
   **no radial face flow** — the spatial streaming term
   :math:`\mu_x \cdot \partial_r` vanishes, so the
   :math:`2|\mu| A_{\rm out}` and :math:`|\mu|(A_{\rm in} +
   A_{\rm out})\,\psi^s_{\rm in}` contributions drop out:

   .. math::

      \mathrm{denom} = (\Delta A / w)\,c_{\rm out} + \Sigma_t\,V,

   .. math::

      \mathrm{numer} = (Q\,V/W)
                       + (\Delta A / w)\,c_{\rm in}\,\psi_{n-\tfrac12},

   with the M-M angular closure still active.  Mirrors
   :func:`orpheus.sn.sweep._sweep_1d_cylindrical` lines 533–546
   verbatim.  In this branch the strategy returns
   ``CellResult(outgoing_spatial_flux=None, ...)`` to signal that
   no downstream face-flux update is meaningful — the sweep
   driver (today inlined; Wave D will dispatch via this Protocol)
   skips the face write when ``outgoing_spatial_flux is None``.

The bit-identical contract — non-negotiable
===========================================

DD's three branches reproduce the sweep's algebra at the **operation
level**.  This means: same operands in the same order, same float
intermediates, ``np.array_equal`` (not ``np.allclose``) holds against
synthetic inputs.  Achieving this requires mirroring the sweep's
``denom`` / ``numer`` build order; rearranging algebraically
equivalent expressions (e.g. computing ``denom`` as
``total_xs * V + 2.0 * abs_mu * A_out + dA_w * c_out`` instead of
``2.0 * abs_mu * A_out + dA_w * c_out + total_xs * V``) breaks
bit-equality at the ULP level even though the math is identical.

The hand-calc test gate at
``tests/sn/spatial/test_diamond.py`` constructs synthetic
:class:`~orpheus.geometry.reduced_operator.StreamingTerms` instances
from real :func:`~orpheus.geometry.reduced_operator.slab_streaming` /
:func:`~orpheus.geometry.reduced_operator.spherical_streaming` factories
and asserts ``np.array_equal`` against a per-cell scalar formula
that **also** mirrors the sweep's operation order.  If a future
edit to this module breaks the order, the hand-calc test fails by
1 ULP — the canonical signature of operation-order drift.

Linear, NOT positivity preserving
=================================

DD is :math:`\mathcal{O}(\Delta x^2)` accurate (Lewis & Miller §5.3)
and **linear** in ``upstream_state`` and ``source`` — both class-
level traits below reflect this.  However, DD is **not** positivity-
preserving: in thin cells with large source the WDD / spatial-DD
closure :math:`\psi_{\rm out} = 2\overline{\psi} - \psi_{\rm in}`
can produce negative outgoing fluxes from positive inputs (Lewis &
Miller §5.3 exhibits the canonical counter-example).  This is why
:class:`is_positivity_preserving` is ``False``.

A Wave C-extension session will ship :class:`Step`,
:class:`LinearDiscontinuous`, and :class:`ExponentialCharacteristic`
as positivity-preserving alternatives, each with its own MMS
spatial-convergence verification:

* **Step** is positivity-preserving (``is_positivity_preserving =
  True``) but only :math:`\mathcal{O}(\Delta x)` accurate.
* **Linear Discontinuous** is :math:`\mathcal{O}(\Delta x^2)` and
  has higher robustness than DD in optically-thick cells.
* **Exponential Characteristic** is positivity-preserving by
  construction (negative arguments are clipped at machine zero) —
  ``is_linear = False`` because the clip is non-linear.

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
* :class:`~orpheus.geometry.reduced_operator.StreamingTerms` —
  the per-(cell, direction) packet the strategy consumes.
* :func:`orpheus.sn.sweep._sweep_1d_cumprod` (slab),
  :func:`orpheus.sn.sweep._sweep_1d_spherical` (sphere),
  :func:`orpheus.sn.sweep._sweep_1d_cylindrical` (cylinder) — the
  inlined-math reference whose scalar form this module reproduces
  bit-for-bit.
* :doc:`/theory/discrete_ordinates`, "Cell update strategies (the
  strategy contract)" → "Diamond Difference" — the theory page.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

import numpy as np

from orpheus.geometry.reduced_operator import StreamingTerms

from .cell_update import CellResult, CellUpdateBase, CellVisit, UpstreamState


# Threshold for the cylindrical pure-azimuthal degenerate branch.
# Mirrors the threshold in :func:`orpheus.sn.sweep._sweep_1d_cylindrical`
# at line 533 (``elif abs_eta < 1e-15:``) — keeping the constant
# in sync with the sweep is what guarantees bit-equality on degenerate
# ordinates.
_DEGENERATE_ABS_MU_THRESHOLD: float = 1e-15


@dataclass(frozen=True, slots=True)
class DiamondDifference(CellUpdateBase, key="diamond_difference"):
    r"""Diamond-Difference (DD) cell-update strategy.

    A **single** geometry-polymorphic strategy that handles slab,
    sphere, and cylinder by branching on two
    :class:`~orpheus.geometry.reduced_operator.StreamingTerms`
    fields — see the module docstring for the three branches and
    the bit-identical contract against
    :mod:`orpheus.sn.sweep`.

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
    exhibits the canonical counter-example.  Wave C-extension's
    :class:`Step` / :class:`ExponentialCharacteristic` strategies
    will be positivity-preserving alternatives.
    """

    def update(
        self,
        visit: CellVisit,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> CellResult:
        r"""Compute the cell-average flux + downstream states.

        Branches on
        :attr:`visit.streaming_terms.alpha_in <orpheus.geometry.reduced_operator.StreamingTerms.alpha_in>`
        (slab vs curvilinear) and on
        :attr:`visit.streaming_terms.abs_mu <orpheus.geometry.reduced_operator.StreamingTerms.abs_mu>`
        (curvilinear vs cylindrical-degenerate).  Sweep-direction
        resolution is encapsulated in
        :attr:`visit.face_area_downstream <orpheus.sn.spatial.cell_update.CellVisit.face_area_downstream>`
        — no sign-of-:math:`\mu` branching inside this strategy.
        See the module docstring for the per-branch algebra and the
        bit-identical cross-reference into :mod:`orpheus.sn.sweep`.
        """
        st = visit.streaming_terms

        if st.alpha_in is None:
            # Branch 1: slab (Cartesian) — no curvature math, no
            # angular redistribution.  Mirrors the per-cell scalar
            # form of :func:`orpheus.sn.sweep._sweep_1d_cumprod` /
            # :func:`orpheus.sn.sweep._solve_recurrence`
            # (sweep.py:117-123 + 208-222).
            return self._update_slab(
                st, total_xs, source, upstream_state,
            )

        # Below: curvilinear branches (sphere or cylinder).
        # ``alpha_in`` populated implies the full curvature bundle is
        # populated (StreamingTerms factory contract).
        assert st.abs_mu is not None  # populated by factory

        if st.abs_mu < _DEGENERATE_ABS_MU_THRESHOLD:
            # Branch 2: cylindrical pure-azimuthal degenerate.
            # No radial face flow on this cell.  Mirrors
            # :func:`orpheus.sn.sweep._sweep_1d_cylindrical`
            # lines 533-546.  ``visit.face_area_downstream`` is
            # ``None`` here — the math doesn't read face areas.
            return self._update_cylindrical_degenerate(
                st, total_xs, source, upstream_state,
            )

        # Branch 3: curvilinear (sphere or cylinder, non-degenerate).
        # Mirrors :func:`orpheus.sn.sweep._sweep_1d_spherical`
        # lines 350-361 and the structurally identical cylindrical
        # inward / outward branches at sweep.py:511-531 / :548-575.
        # ``visit.face_area_downstream`` is the sweep-direction-
        # resolved outgoing face area.
        assert visit.face_area_downstream is not None, (
            "Curvilinear non-degenerate cell update requires a "
            "resolved face_area_downstream on the CellVisit packet."
        )
        return self._update_curvilinear(
            st,
            visit.face_area_downstream,
            total_xs,
            source,
            upstream_state,
        )

    # ── Slab ───────────────────────────────────────────────────────

    @staticmethod
    def _update_slab(
        st: StreamingTerms,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> CellResult:
        r"""Slab DD per-cell update.

        Mirrors :func:`orpheus.sn.sweep._solve_recurrence` (the
        per-cell algebra it solves over a whole row via cumprod)
        applied to a single cell.  The recurrence is

        .. math::

           \psi_{\rm out} = a\,\psi_{\rm in} + b\,Q,
           \qquad a = (2|\mu| - \Delta x\,\Sigma_t)/\mathrm{denom},
           \quad b = 2\,\Delta x\,(1/W)/\mathrm{denom},

        with ``denom = 2|μ| + Δx·Σ_t``.  Per the strategy contract,
        ``source`` arrives already weight-normalised
        (:math:`Q \cdot \Delta x / W` for slab), so ``b · Q``
        becomes ``2 · source / denom``.  The cell-average is
        ``0.5 · (psi_in + psi_out)`` — operation order matches the
        sweep's line 222.
        """
        assert st.abs_mu is not None  # populated by slab factory
        abs_mu = st.abs_mu
        chord = st.chord_length
        # ``denom`` build order matches sweep.py:119
        # (``2.0 * mu_pos[:, None, None] + dx[None, :, None] *
        # sig_t_1d[None, :, :]``).
        denom = 2.0 * abs_mu + chord * total_xs        # (ng,)
        # Stream coefficient ``a`` matches sweep.py:120-122
        # (``(2.0 * mu_pos - dx * sig_t_1d) / denom``).
        a = (2.0 * abs_mu - chord * total_xs) / denom  # (ng,)
        # Source term ``s = b · Q`` where the per-cell ``b`` for
        # slab is ``2 · Δx · weight_norm / denom``.  ``source``
        # already carries ``Q · Δx · weight_norm`` (contract;
        # for slab ``V == Δx``), so ``s`` is ``2 · source / denom``.
        # Build order matches sweep.py:123 + sweep.py:135 (``bQ =
        # source_coeff * Q_1d``) collapsed for a single cell.
        s = 2.0 * source / denom                       # (ng,)
        psi_in = upstream_state.spatial_upstream
        # Operation order matches sweep.py:221 (``psi_out = a *
        # psi_in + s``).
        psi_out = a * psi_in + s
        # Operation order matches sweep.py:222 (``return 0.5 *
        # (psi_in + psi_out)``).
        psi_avg = 0.5 * (psi_in + psi_out)
        return CellResult(
            cell_average_flux=psi_avg,
            outgoing_spatial_flux=psi_out,
            outgoing_angular_state=None,
        )

    # ── Curvilinear (sphere / cylinder, non-degenerate) ────────────

    @staticmethod
    def _update_curvilinear(
        st: StreamingTerms,
        A_downstream: float,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> CellResult:
        r"""Curvilinear DD per-cell update (sphere / non-degenerate cylinder).

        Mirrors :func:`orpheus.sn.sweep._sweep_1d_spherical` lines
        350-361 (and the structurally identical cylindrical
        inward / outward branches at sweep.py:511-531 / :548-575).
        ``c_out = α_out / τ`` and ``c_in = (1 - τ) / τ · α_out +
        α_in`` are the M-M closure constants.

        Direction resolution
        --------------------

        Reads :attr:`StreamingTerms.face_area_inner` /
        :attr:`StreamingTerms.face_area_outer` (geometric labels —
        inner = closer to :math:`r=0`) and the sweep-direction-
        resolved ``A_downstream`` (which face is downstream for this
        visit).  The DD formula uses the symmetric sum
        :math:`A_{\rm in} + A_{\rm out} \equiv A_{\rm inner} +
        A_{\rm outer}` on the :math:`|\mu| \cdot \psi^s_{\rm in}`
        term, and the asymmetric :math:`A_{\rm out}` (= downstream)
        on the :math:`2|\mu|` term in the denominator.  No sign-of-
        :math:`\mu` branching here — the sweep orchestrator already
        resolved which face is downstream before issuing the
        :class:`CellVisit`.
        """
        # Pull populated curvilinear fields off the streaming-terms
        # packet.  Asserts that the factory contract is honoured.
        assert st.abs_mu is not None
        assert st.face_area_inner is not None
        assert st.face_area_outer is not None
        assert st.delta_A_over_w is not None
        assert st.alpha_in is not None
        assert st.alpha_out is not None
        assert st.tau_mm is not None
        assert st.volume is not None
        assert upstream_state.angular_upstream is not None, (
            "Curvilinear cell update requires an upstream angular state."
        )

        abs_mu = st.abs_mu
        A_inner = st.face_area_inner
        A_outer = st.face_area_outer
        # Symmetric sum is invariant under inner/outer swap, so the
        # build order matches sweep.py:354 ``(A_in + A_out)``
        # bit-identically regardless of sweep direction.
        A_total = A_inner + A_outer
        dA_w = st.delta_A_over_w
        alpha_in = st.alpha_in
        alpha_out = st.alpha_out
        tau = st.tau_mm
        V = st.volume

        # Closure-prefactor combinations.  Operation order matches
        # sweep.py:328-329 (``c_out = alpha_out / tau_n``;
        # ``c_in = (1.0 - tau_n) / tau_n * alpha_out + alpha_in``).
        c_out = alpha_out / tau
        c_in = (1.0 - tau) / tau * alpha_out + alpha_in

        psi_spat_in = upstream_state.spatial_upstream
        psi_angle_in = upstream_state.angular_upstream

        # Denominator build order matches sweep.py:350-352 — the
        # ``A_out`` here is the sweep-direction-resolved
        # ``A_downstream`` (outgoing face area).
        # (``2.0 * abs_mu * A_out + dA_w * c_out + sig_t_1d[i] *
        # V[i]``).
        denom = 2.0 * abs_mu * A_downstream + dA_w * c_out + total_xs * V

        # Numerator build order matches sweep.py:353-355
        # (``QV[i] + abs_mu * (A_in + A_out) * psi_spatial_in +
        # dA_w * c_in * psi_angle[i]``).
        numer = (
            source
            + abs_mu * A_total * psi_spat_in
            + dA_w * c_in * psi_angle_in
        )

        # Cell-average solve (sweep.py:357 / :385 / :523 / :564).
        psi_avg = numer / denom

        # WDD spatial closure (sweep.py:360 / :388 / :525 / :566).
        psi_spat_out = 2.0 * psi_avg - psi_spat_in
        # M-M angular closure (sweep.py:361 / :389 / :526 / :567).
        psi_angle_out = (psi_avg - (1.0 - tau) * psi_angle_in) / tau

        return CellResult(
            cell_average_flux=psi_avg,
            outgoing_spatial_flux=psi_spat_out,
            outgoing_angular_state=psi_angle_out,
        )

    # ── Cylindrical pure-azimuthal degenerate ──────────────────────

    @staticmethod
    def _update_cylindrical_degenerate(
        st: StreamingTerms,
        total_xs: np.ndarray,
        source: np.ndarray,
        upstream_state: UpstreamState,
    ) -> CellResult:
        r"""Cylindrical pure-azimuthal degenerate per-cell update.

        Used when ``abs_mu < 1e-15`` — the level's axial cosine
        :math:`|\mu_z| \to 1` so the radial cosine :math:`|\eta|
        \to 0` and there is no radial face flow.  Mirrors
        :func:`orpheus.sn.sweep._sweep_1d_cylindrical` lines
        533-546.  ``outgoing_spatial_flux`` is set to ``None`` to
        signal "no face-flux write" to the sweep driver.
        """
        assert st.delta_A_over_w is not None
        assert st.alpha_in is not None
        assert st.alpha_out is not None
        assert st.tau_mm is not None
        assert st.volume is not None
        assert upstream_state.angular_upstream is not None, (
            "Cylindrical-degenerate cell update requires an upstream "
            "angular state."
        )

        dA_w = st.delta_A_over_w
        alpha_in = st.alpha_in
        alpha_out = st.alpha_out
        tau = st.tau_mm
        V = st.volume

        # Same M-M closure constants as the curvilinear branch
        # (sweep.py:328-329).
        c_out = alpha_out / tau
        c_in = (1.0 - tau) / tau * alpha_out + alpha_in

        psi_angle_in = upstream_state.angular_upstream

        # Denominator (sweep.py:538:  ``denom = dA_w * c_out +
        # sig_t_1d[i] * V[i]``).  No 2|μ| A_out term — no radial
        # face flow.
        denom = dA_w * c_out + total_xs * V

        # Numerator (sweep.py:539-540: ``numer = QV[i] + dA_w *
        # c_in * psi_angle[i]``).  No |μ|(A_in+A_out)·ψ_spatial_in
        # term — no radial face flow.
        numer = source + dA_w * c_in * psi_angle_in

        psi_avg = numer / denom

        # M-M angular closure remains active (sweep.py:543).
        psi_angle_out = (psi_avg - (1.0 - tau) * psi_angle_in) / tau

        # outgoing_spatial_flux=None signals to the sweep driver
        # that no face-flux update is meaningful — there is no
        # radial face flow on this cell.
        return CellResult(
            cell_average_flux=psi_avg,
            outgoing_spatial_flux=None,
            outgoing_angular_state=psi_angle_out,
        )


__all__ = ["DiamondDifference"]
