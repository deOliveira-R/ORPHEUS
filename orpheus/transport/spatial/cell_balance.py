r"""Unified per-cell balance algebra — geometry-blind by data.

One function — :func:`cell_balance_for_streaming` — computes the per-cell
balance ``(denom, numer_upstream)`` for slab, sphere, cylinder
(non-degenerate), and cylindrical pure-azimuthal degenerate cells,
vectorized over an ordinate mask.  Geometry is data (face areas, volume,
|μ|); the angular closure's contributions arrive ASSEMBLED
(``angular_denom_term`` / ``angular_numer_upstream``), so the helper is
blind to both the geometry kind and the closure family.  Every consumer
routes here: the vectorized matvec directly, and the per-cell
solve/apply pair through ``DiamondDifference``'s ``n_mask=1`` bridge
(P4.9a retired the scalar twin ``cell_balance_terms``, which duplicated
this algebra under M-M-specific argument names).

Mathematical content
====================

The canonical curvilinear DD cell balance (Hébert §3.9.3 cylinder /
§3.9.4 sphere; Lewis-Miller §6.1), closed with the Morel-Montry
weighted angular closure (Morel & Montry 1984; implemented form
Bailey-Morel-Chang 2010 Eqs. (42)/(43) — Hébert defines no tau and
ships the plain angular diamond):

.. math::

    \mathrm{denom} \cdot \bar\psi \;=\; q \;+\; \mathrm{numer\_upstream}

where

.. math::

    \mathrm{denom} &= 2|\mu|\,A_{\text{down}}
                      + \frac{\Delta A}{w}\,c_{\text{out}}
                      + \Sigma_t\,V,\\
    \mathrm{numer\_upstream} &= |\mu|\,(A_{\text{in}} + A_{\text{out}})
                                   \,\psi^s_{\text{in}}
                                + \frac{\Delta A}{w}\,c_{\text{in}}\,
                                  \psi^{\theta}_{\text{in}},\\
    c_{\text{out}} &= \alpha_{\text{out}} / \tau,\\
    c_{\text{in}}  &= \tfrac{1-\tau}{\tau}\,\alpha_{\text{out}}
                      + \alpha_{\text{in}}.

Pre-Step-2.5 this lived in two helpers (``cell_balance_terms`` for
the non-degenerate curvilinear branch and
``cell_balance_terms_degenerate`` for the cylindrical pure-azimuthal
case) plus an inlined slab recurrence inside ``DiamondDifference``.
Step 2.5 collapsed all three into one SCALAR helper; PR-TYPED-6.5 Phase
2.11 added this vectorized closure-blind form beside it, and P4.9a
retired the scalar one onto it (2026-08-28).  The geometry cases:

* **Slab** — neutral curvature populated by the
  :func:`~orpheus.sn.mesh.reduced_operator.slab_streaming`
  factory: ``face_area_inner = face_area_outer = 1.0``, and the
  assembled ``angular_*`` kwargs arrive as zeros (the ΔA/w coupling —
  closure-owned and cache-interned since P4.9a/P4.7 — is identically
  zero on a slab; the former ``alpha_*`` / ``tau_mm`` packing left
  the packet at #236 Step C).  ``CellVisit.face_area_downstream =
  1.0``.
  The denominator collapses to ``2|μ|·1 + 0 + Σ_t·V`` = the slab
  form; the upstream numerator collapses to ``|μ|·2·ψ^s_in`` =
  ``2|μ|·ψ^s_in``.
* **Curvilinear non-degenerate** — physical curvature values.
  ``CellVisit.face_area_downstream`` is the sweep-direction-
  resolved outgoing face area.
* **Cylindrical pure-azimuthal degenerate** —
  ``CellVisit.face_area_downstream = 0.0`` (geometric truth: no
  radial face on this ordinate).  The ``2|μ|·face_area_downstream`` term
  vanishes via ``face_area_downstream = 0`` (not via the numerical threshold
  ``|μ| < 1e-15``).

The solve branch returns ``psi_avg = (source + numer_upstream) /
denom``; the apply branch returns ``denom · cell_avg − (source +
numer_upstream)``.  Same algebra, two consumers (Pattern 2 — no
twin paths).

See :mod:`orpheus.transport.spatial.cell_balance` test gate at
``tests/sn/sweep/core/test_cell_balance_for_streaming.py``.
"""

from __future__ import annotations

import numpy as np


def cell_balance_for_streaming(
    *,
    abs_mu: np.ndarray,                # (n_mask,)  |μ_n| for ordinates in the mask
    A_downstream: "np.ndarray | float",  # (n_mask,)  sweep-resolved outgoing face
    #                                      area; broadcast-scalar 0.0 for the
    #                                      degenerate pure-z / pole cell (no
    #                                      outgoing face), like face_area_total
    face_area_total: np.ndarray,               # (n_mask,)  A_inner + A_outer  (broadcast scalar OK)
    total_xs: np.ndarray,              # (ng,)      per-group total cross section
    volume: float,                     # scalar     cell volume V_i
    psi_face_in: np.ndarray,           # (ng, n_mask) face flux entering this cell per ordinate
    angular_denom_term: np.ndarray,    # (n_mask,)  closure contribution to denom
    angular_numer_upstream: np.ndarray,  # (ng, n_mask) closure contribution to upstream numer
) -> tuple[np.ndarray, np.ndarray]:
    r"""Vectorized per-cell balance for one cell over an ORDINATE MASK.

    Issue #197 PR-TYPED-6 — the shared cell-balance algebra both
    :meth:`DiamondDifference.residual` (per-ordinate, ``n_mask=1``)
    and the unified matvec (vectorized over ordinates in a sweep
    direction) consume.  Pattern 2 — single source of truth.

    PR-TYPED-6.5 Phase 2.11 — geometry-blind interface
    ---------------------------------------------------
    Pre-Phase-2.11 the helper accepted M-M-specific arguments
    ``delta_A_over_w``, ``c_in``, ``c_out``, ``psi_angular_upstream``.  Phase 2.11
    pushes those names into the
    :class:`~orpheus.sn.angular.closure.AngularClosureBase`
    strategy: the matvec body calls
    ``closure.cell_contribution(...)`` to obtain ``(angular_denom_term,
    angular_numer_upstream)`` — the operated-on contributions in the
    shapes this helper needs.  ``cell_balance_for_streaming`` is now
    geometry-blind by interface (no M-M names; the closure produces
    zero contributions for slab via
    :class:`~orpheus.sn.angular.closure.IdentityAngularClosure`).

    Computes:

    .. math::

       \mathrm{denom}[g, n] &= 2\,|\mu_n|\,A_{\rm down,n}
                              + \mathrm{angular\_denom\_term}[n]
                              + \Sigma_{t,g}\,V,\\
       \mathrm{numer\_upstream}[g, n] &= |\mu_n|\,A_{\rm total,n}\,
                                          \psi^s_{\rm in,g,n}
                                        + \mathrm{angular\_numer\_upstream}[g, n].

    For M-M closure, ``angular_denom_term = (ΔA/w)·c_out`` and
    ``angular_numer_upstream = (ΔA/w)·c_in·ψ_{m-1/2,i,g}``.  For
    Identity closure (slab), both are zero arrays of the right shape.

    Conversion to operator residual
    -------------------------------

    At a probe value ``psi_cell_avg`` (shape ``(ng, n_mask)``), the
    per-cell **cell-balance residual** is

    .. math::

       r[g, n] \;=\; \mathrm{denom}[g, n] \cdot \overline\psi[g, n]
                     \;-\; \mathrm{numer\_upstream}[g, n].

    Converting to the matvec convention (``[neutrons/cm^2/s]`` per
    cell volume) divides by ``V``: ``L_apply[g, n] = r[g, n] / V``.
    At ``source=0`` this equals the matvec body output for the
    streaming + redistribution + collision triplet (Resolution A
    factoring at ``L+C``).

    Pattern 3 — named intermediates
    -------------------------------

    * ``streaming_denom_term`` — :math:`2|\mu|\,A_{\rm down}` (WDD
      outflow contribution to the denominator).
    * ``angular_denom_term`` (arg) — closure's contribution to the
      denominator.  M-M: :math:`(\Delta A/w)\,c_{\rm out}`.  Identity:
      zero.
    * ``collision_denom_term`` — :math:`\Sigma_t\,V` (the per-group
      collision contribution; broadcast across ``n_mask``).
    * ``spatial_upstream_term`` — :math:`|\mu|\,A_{\rm total}\,
      \psi^s_{\rm in}` (the upstream face-flux contribution).
    * ``angular_numer_upstream`` (arg) — closure's contribution to
      the upstream numerator.  M-M: :math:`(\Delta A/w)\,c_{\rm in}\,
      \psi^\theta_{\rm in}`.  Identity: zero.

    Parameters
    ----------
    abs_mu :
        Shape ``(n_mask,)``.  Absolute direction cosines for the
        ordinates in this mask.  Already sweep-direction-resolved
        (no sign flips inside the helper).
    A_downstream :
        Shape ``(n_mask,)``.  Sweep-direction-resolved outgoing face
        area for each ordinate (``face_area_outer`` for outward
        sweeps, ``face_area_inner`` for inward sweeps).  Slab carries
        ``1.0`` (neutral); cylindrical-degenerate ordinates carry
        ``0.0`` (no spatial flow).
    face_area_total :
        Shape ``(n_mask,)`` (or broadcast scalar).  Sum of the cell's
        inner and outer face areas.  Slab: ``2.0`` (1+1).
    total_xs :
        Shape ``(ng,)``.  Per-group total cross section
        :math:`\Sigma_t` on this cell.
    volume :
        Cell volume :math:`V_i`.  Scalar.
    psi_face_in :
        Shape ``(ng, n_mask)``.  Face flux entering this cell, per
        group, per ordinate-in-mask.  Already sweep-direction-resolved.
    angular_denom_term :
        Shape ``(n_mask,)``.  Closure's contribution to the denominator
        (post-Phase-2.11 strategy abstraction).  Typically the output
        of ``closure.cell_contribution(...)[0]``.
    angular_numer_upstream :
        Shape ``(ng, n_mask)``.  Closure's contribution to the
        upstream numerator.  Typically ``closure.cell_contribution(...)[1]``.

    Returns
    -------
    denom : np.ndarray
        Shape ``(ng, n_mask)``.  Per-group per-ordinate-in-mask
        cell-balance denominator.
    numer_upstream : np.ndarray
        Shape ``(ng, n_mask)``.  Per-group per-ordinate-in-mask
        cell-balance upstream numerator (excludes the volumetric
        source).
    """
    # Denominator terms — each named with its physical role.
    streaming_denom_term = 2.0 * abs_mu * A_downstream  # (n_mask,)
    collision_denom_term = total_xs * volume            # (ng,)

    # Denominator: broadcast (ng,) collision against (n_mask,) streaming +
    # angular (from closure) to produce (ng, n_mask).
    denom = (
        (streaming_denom_term + angular_denom_term)[None, :]
        + collision_denom_term[:, None]
    )

    # Upstream-numerator: ``face_area_total`` may be a scalar broadcast
    # (face area is per-cell, shared across ordinates in a sweep level).
    spatial_upstream_term = (
        abs_mu[None, :] * face_area_total * psi_face_in
    )                                                    # (ng, n_mask)
    numer_upstream = spatial_upstream_term + angular_numer_upstream

    return denom, numer_upstream


__all__ = [
    "cell_balance_for_streaming",
]
