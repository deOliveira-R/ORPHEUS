r"""Unified per-cell balance algebra — geometry-blind by data (Step 2.5).

Issue #196 Phase G Step 2.5.  One function — :func:`cell_balance_terms` —
computes the algebraic intermediates of the per-cell DD balance for
slab, sphere, cylinder (non-degenerate), and cylindrical pure-azimuthal
degenerate cells.  Geometry is data carried by
:class:`~orpheus.geometry.reduced_operator.StreamingTerms` and
:class:`~orpheus.sn.spatial.cell_update.CellVisit`; the helper does NOT
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

See :mod:`orpheus.sn.spatial.cell_balance` test gate at
``tests/sn/spatial/test_cell_balance.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from orpheus.geometry.reduced_operator import StreamingTerms

    from .cell_update import UpstreamState


@dataclass(frozen=True, slots=True)
class CellBalanceTerms:
    r"""Algebraic intermediates of the unified per-cell balance.

    Constructed by :func:`cell_balance_terms`; consumed by both
    :class:`~orpheus.sn.spatial.diamond.DiamondDifference.update`
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


def cell_balance_terms(
    st: "StreamingTerms",
    A_downstream: float,
    total_xs: np.ndarray,
    upstream_state: "UpstreamState",
) -> CellBalanceTerms:
    r"""Unified per-cell balance terms — geometry-blind by data.

    Issue #196 Phase G Step 2.5: ONE helper for slab, sphere,
    cylinder (non-degenerate), and cylindrical pure-azimuthal
    degenerate.  The geometry kind is encoded in the populated
    fields of ``st`` (neutral values for slab, physical values for
    curvilinear) and the value of ``A_downstream`` (``0.0`` for
    cylindrical-degenerate, ``1.0`` for slab, the physical outgoing
    face area for curvilinear).

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
    tau = st.tau_mm
    V = st.volume

    # M-M closure constants — degenerate to zero for slab (neutral
    # alpha values) and clamped to (½, 1] for curvilinear.
    c_out = st.alpha_out / tau
    c_in = (1.0 - tau) / tau * st.alpha_out + st.alpha_in

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
    "cell_balance_terms",
]
