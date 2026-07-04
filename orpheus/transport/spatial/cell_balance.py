r"""Unified per-cell balance algebra — geometry-blind by data (Step 2.5).

Issue #196 Phase G Step 2.5.  One function — :func:`cell_balance_terms` —
computes the algebraic intermediates of the per-cell DD balance for
slab, sphere, cylinder (non-degenerate), and cylindrical pure-azimuthal
degenerate cells.  Geometry is data carried by
:class:`~orpheus.geometry.reduced_operator.StreamingTerms` and
:class:`~orpheus.transport.spatial.scheme.CellVisit`; the helper does NOT
branch on geometry kind.

Mathematical content
====================

The canonical Hébert §3.9.4 / Lewis-Miller §6.1 curvilinear DD
balance with Morel-Montry angular closure:

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
Step 2.5 collapses all three into this single helper:

* **Slab** — neutral curvature populated by the
  :func:`~orpheus.geometry.reduced_operator.slab_streaming`
  factory: ``face_area_inner = face_area_outer = 1.0``,
  ``delta_A_over_w = 0.0``, ``alpha_in = alpha_out = 0.0``,
  ``tau_mm = 1.0``.  ``CellVisit.face_area_downstream = 1.0``.
  The denominator collapses to ``2|μ|·1 + 0 + Σ_t·V`` = the slab
  form; the upstream numerator collapses to ``|μ|·2·ψ^s_in`` =
  ``2|μ|·ψ^s_in``.
* **Curvilinear non-degenerate** — physical curvature values.
  ``CellVisit.face_area_downstream`` is the sweep-direction-
  resolved outgoing face area.
* **Cylindrical pure-azimuthal degenerate** —
  ``CellVisit.face_area_downstream = 0.0`` (geometric truth: no
  radial face on this ordinate).  The ``2|μ|·A_down`` term
  vanishes via ``A_down = 0`` (not via the numerical threshold
  ``|μ| < 1e-15``).

The solve branch returns ``psi_avg = (source + numer_upstream) /
denom``; the apply branch returns ``denom · cell_avg − (source +
numer_upstream)``.  Same algebra, two consumers (Pattern 2 — no
twin paths).

See :mod:`orpheus.transport.spatial.cell_balance` test gate at
``tests/sn/sweep/core/test_cell_balance_for_streaming.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from orpheus.geometry.reduced_operator import StreamingTerms

    from .scheme import UpstreamState


@dataclass(frozen=True, slots=True)
class CellBalanceTerms:
    r"""Algebraic intermediates of the unified per-cell balance.

    Constructed by :func:`cell_balance_terms`; consumed by both
    :class:`~orpheus.transport.spatial.diamond.DiamondDifference.update`
    (solve branch — divides ``(source + numer_upstream) / denom``)
    and :meth:`DiamondDifference.residual` (apply branch —
    computes ``denom · cell_avg − (source + numer_upstream)``).
    Pattern 2 — single source of truth for the cell-balance algebra
    across slab, curvilinear, and cylindrical-degenerate geometries.

    Parameters
    ----------
    denom :
        Denominator of the unified DD balance,
        :math:`2|\mu|\,A_{\text{down}} + (\Delta A / w)\,c_{\text{out}}
        + \Sigma_t\,V`.  Per-group, shape ``(ng,)``.
    numer_upstream :
        Upstream-contribution numerator (excluding the cell source),
        :math:`|\mu|\,(A_{\text{in}} + A_{\text{out}})\,\psi^s_{\text{in}}
        + (\Delta A / w)\,c_{\text{in}}\,\psi^{\theta}_{\text{in}}`.
        Per-group, shape ``(ng,)``.
    c_in :
        Morel-Montry "in" closure constant
        :math:`c_{\text{in}} = (1-\tau)/\tau \cdot \alpha_{\text{out}}
        + \alpha_{\text{in}}`.  Scalar.  ``0.0`` for slab.
    c_out :
        Morel-Montry "out" closure constant
        :math:`c_{\text{out}} = \alpha_{\text{out}} / \tau`.  Scalar.
        ``0.0`` for slab.
    """

    denom: np.ndarray
    numer_upstream: np.ndarray
    c_in: float
    c_out: float


def cell_balance_for_streaming(
    *,
    abs_mu: np.ndarray,                # (n_mask,)  |μ_n| for ordinates in the mask
    A_downstream: "np.ndarray | float",  # (n_mask,)  sweep-resolved outgoing face
    #                                      area; broadcast-scalar 0.0 for the
    #                                      degenerate pure-z / pole cell (no
    #                                      outgoing face), like A_total
    A_total: np.ndarray,               # (n_mask,)  A_inner + A_outer  (broadcast scalar OK)
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
    ``dA_w``, ``c_in``, ``c_out``, ``psi_angular_upstream``.  Phase 2.11
    pushes those names into the
    :class:`~orpheus.sn.spatial.pole_angular_closure.PoleAngularClosureBase`
    strategy: the matvec body calls
    ``closure.cell_contribution(...)`` to obtain ``(angular_denom_term,
    angular_numer_upstream)`` — the operated-on contributions in the
    shapes this helper needs.  ``cell_balance_for_streaming`` is now
    geometry-blind by interface (no M-M names; the closure produces
    zero contributions for slab via
    :class:`~orpheus.sn.spatial.pole_angular_closure.IdentityAngularClosure`).

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
    A_total :
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

    # Upstream-numerator: ``A_total`` may be a scalar broadcast
    # (face area is per-cell, shared across ordinates in a sweep level).
    spatial_upstream_term = (
        abs_mu[None, :] * A_total * psi_face_in
    )                                                    # (ng, n_mask)
    numer_upstream = spatial_upstream_term + angular_numer_upstream

    return denom, numer_upstream


def cell_balance_terms(
    st: "StreamingTerms",
    A_downstream: float,
    total_xs: np.ndarray,
    upstream_state: "UpstreamState",
    *,
    c_in: float,
    c_out: float,
) -> CellBalanceTerms:
    r"""Unified per-cell balance terms — geometry-blind by data.

    Issue #196 Phase G Step 2.5: ONE helper for slab, sphere,
    cylinder (non-degenerate), and cylindrical pure-azimuthal
    degenerate.  The geometry kind is encoded in the populated
    fields of ``st`` (neutral values for slab, physical values for
    curvilinear) and the value of ``A_downstream`` (``0.0`` for
    cylindrical-degenerate, ``1.0`` for slab, the physical outgoing
    face area for curvilinear).

    Issue #236 Phase 2 B3: the Morel--Montry weighted-diamond
    constants ``c_in`` / ``c_out`` are ANGULAR-CLOSURE-OWNED inputs,
    passed by the caller from :attr:`CellVisit.c_in` / ``c_out`` (the
    closure stamps them from its ``c_{in,out}_per_ordinate`` accessors
    via :meth:`SNMesh._make_cell_visit`).  They are NO LONGER rebuilt
    here from ``st.alpha_*`` / ``st.tau_mm`` — this helper does not read
    τ at all.  The spatial scheme receives the angular constants as
    DATA, preserving the spatial :math:`\otimes` angular separation.

    Parameters
    ----------
    st :
        Streaming-terms packet on the cell-visit.  All curvature
        fields populated (slab carries neutral values per
        :func:`~orpheus.geometry.reduced_operator.slab_streaming`).
    A_downstream :
        Sweep-direction-resolved outgoing face area.  Read from the
        :class:`CellVisit`'s ``face_area_downstream`` field.  ``1.0``
        for slab, ``0.0`` for cyl-degenerate, physical face area for
        curvilinear.
    total_xs :
        Per-group total cross section, shape ``(ng,)``.
    upstream_state :
        Per-cell upstream state.  ``spatial_upstream`` is always
        populated.  ``angular_upstream`` is ``None`` for slab and
        populated (``(ng,)``) for curvilinear.
    c_in :
        Angular-closure-owned Morel--Montry "in" constant
        :math:`(1-\tau)/\tau\,\alpha_{\rm out} + \alpha_{\rm in}` for this
        ordinate — read from :attr:`CellVisit.c_in`.  ``0.0`` for slab.
    c_out :
        Angular-closure-owned Morel--Montry "out" constant
        :math:`\alpha_{\rm out}/\tau` for this ordinate — read from
        :attr:`CellVisit.c_out`.  ``0.0`` for slab.

    Returns
    -------
    CellBalanceTerms
        Bundled algebraic intermediates.  See class docstring.
    """
    # All curvature fields populated for every geometry (slab carries
    # neutral values per Step 2.5 — see slab_streaming factory).
    abs_mu = st.abs_mu
    A_total = st.face_area_inner + st.face_area_outer
    dA_w = st.delta_A_over_w
    V = st.volume

    # The M-M closure constants c_in / c_out arrive as ANGULAR-CLOSURE-
    # OWNED data (Issue #236 Phase 2 B3) — the closure produces them from
    # its α-dome / τ and stamps them on the CellVisit; this helper consumes
    # them, it does not derive them.  Slab carries 0.0 (neutral identity
    # closure); curvilinear carries the physical weighted-diamond values.

    # Denominator: 2|μ|·A_down (vanishes for cyl-degenerate via
    # A_down=0; reduces to 2|μ|·1 for slab) + curvature redistribution
    # (vanishes for slab via dA_w=0) + collision.
    denom = 2.0 * abs_mu * A_downstream + dA_w * c_out + total_xs * V

    # Upstream numerator: spatial-streaming contribution (slab:
    # |μ|·2·ψ^s_in = 2|μ|·ψ^s_in; cyl-deg: zero by A_down=0 path —
    # but A_total is non-zero for cyl-deg, so we use A_down's
    # presence as the discriminator. Wait — cyl-deg has populated
    # A_inner+A_outer but no radial flow, so the |μ|·A_total·ψ^s_in
    # term must vanish. abs_mu < 1e-15 ensures it ≈ 0.).
    #
    # Angular-redistribution contribution: zero for slab via
    # angular_upstream=None branch below; physical for curvilinear.
    psi_ang = upstream_state.angular_upstream
    ang_contrib = (
        0.0
        if psi_ang is None
        else dA_w * c_in * psi_ang
    )
    numer_upstream = abs_mu * A_total * upstream_state.spatial_upstream + ang_contrib

    return CellBalanceTerms(
        denom=denom,
        numer_upstream=numer_upstream,
        c_in=c_in,
        c_out=c_out,
    )


__all__ = [
    "CellBalanceTerms",
    "cell_balance_for_streaming",
    "cell_balance_terms",
]
