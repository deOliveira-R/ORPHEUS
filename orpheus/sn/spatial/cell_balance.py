r"""Curvilinear per-cell balance terms — shared algebra for solve / apply.

Issue #196 Phase G Step 2 (qa CONCERN-2 fix).  The Diamond Difference's
``_update_curvilinear`` (which divides ``numer / denom`` to obtain the
cell-average flux) and ``SNCellOperator._apply_curvilinear_residual``
(which subtracts ``denom · cell_avg - (source + numer_upstream)`` to
obtain the residual) build the same algebraic intermediates.  Pattern 2
(no twin paths) requires the algebra to live in **one** place.  This
module factors that algebra into :func:`cell_balance_terms`, which
both consumers call.

The mathematics is the canonical Hébert §3.9.4 / Lewis-Miller §6.1
curvilinear DD balance with Morel-Montry angular closure:

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

The solve branch returns ``psi_avg = (source + numer_upstream) /
denom``; the apply branch returns ``denom · cell_avg − (source +
numer_upstream)``.  Same algebra, two consumers.

See :mod:`orpheus.sn.spatial.cell_balance` test gate at
``tests/sn/spatial/test_cell_balance.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from .cell_update import StreamingTerms, UpstreamState


@dataclass(frozen=True, slots=True)
class CellBalanceTerms:
    r"""Algebraic intermediates of the curvilinear per-cell balance.

    Constructed by :func:`cell_balance_terms`; consumed by both
    :class:`~orpheus.sn.spatial.diamond.DiamondDifference._update_curvilinear`
    (solve branch) and
    :func:`~orpheus.sn.spatial.operators._apply_curvilinear_residual`
    (apply branch).  Closes the twin-path duplication (Pattern 2)
    identified by qa CONCERN-2 in the Step 1 review.

    Parameters
    ----------
    denom :
        Denominator of the curvilinear DD balance,
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
        + \alpha_{\text{in}}`.  Scalar.
    c_out :
        Morel-Montry "out" closure constant
        :math:`c_{\text{out}} = \alpha_{\text{out}} / \tau`.  Scalar.
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
    r"""Compute the curvilinear non-degenerate per-cell balance terms.

    Algebraic core of the Hébert §3.9.4 / Lewis-Miller §6.1 curvilinear
    DD balance with Morel-Montry angular closure.  Used by

    * :meth:`DiamondDifference._update_curvilinear` — the solve branch,
      which divides ``(source + terms.numer_upstream) / terms.denom``
      to recover the cell-average flux :math:`\bar\psi`.
    * :func:`_apply_curvilinear_residual` — the apply branch, which
      computes ``terms.denom · cell_avg − (source + terms.numer_upstream)``
      as the cell-balance residual.

    Operation order matches :func:`orpheus.sn.sweep._sweep_1d_spherical`
    lines 328-355 for bit-identity inheritance from the pre-Wave-D
    inlined sweep math.

    Parameters
    ----------
    st :
        Streaming-terms packet on the cell-visit.  Must have populated
        curvilinear fields (``abs_mu``, ``face_area_inner``,
        ``face_area_outer``, ``delta_A_over_w``, ``alpha_in``,
        ``alpha_out``, ``tau_mm``, ``volume``).
    A_downstream :
        Sweep-direction-resolved outgoing face area.  Read from the
        :class:`CellVisit`'s ``face_area_downstream`` field; equal to
        ``face_area_outer`` on the outward sweep and ``face_area_inner``
        on the inward sweep.
    total_xs :
        Per-group total cross section, shape ``(ng,)``.
    upstream_state :
        Per-cell upstream state.  Must have non-``None``
        ``spatial_upstream`` (the incoming face flux) and
        ``angular_upstream`` (the incoming half-angle flux from the
        previous ordinate's M-M recurrence step).

    Returns
    -------
    CellBalanceTerms
        Bundled algebraic intermediates.  See class docstring.
    """
    assert st.abs_mu is not None
    assert st.face_area_inner is not None
    assert st.face_area_outer is not None
    assert st.delta_A_over_w is not None
    assert st.alpha_in is not None
    assert st.alpha_out is not None
    assert st.tau_mm is not None
    assert st.volume is not None
    assert upstream_state.angular_upstream is not None, (
        "Curvilinear cell balance requires an upstream angular state."
    )

    abs_mu = st.abs_mu
    A_inner = st.face_area_inner
    A_outer = st.face_area_outer
    # Symmetric sum invariant under inner/outer swap — bit-identical to
    # sweep.py:354's ``(A_in + A_out)`` regardless of sweep direction.
    A_total = A_inner + A_outer
    dA_w = st.delta_A_over_w
    alpha_in = st.alpha_in
    alpha_out = st.alpha_out
    tau = st.tau_mm
    V = st.volume

    # Closure-prefactor combinations (sweep.py:328-329).
    c_out = alpha_out / tau
    c_in = (1.0 - tau) / tau * alpha_out + alpha_in

    psi_spat_in = upstream_state.spatial_upstream
    psi_angle_in = upstream_state.angular_upstream

    # Denominator (sweep.py:350-352).
    denom = 2.0 * abs_mu * A_downstream + dA_w * c_out + total_xs * V

    # Upstream contribution to the numerator (sweep.py:353-355,
    # excluding the cell source which the consumer adds in).
    numer_upstream = (
        abs_mu * A_total * psi_spat_in + dA_w * c_in * psi_angle_in
    )

    return CellBalanceTerms(
        denom=denom,
        numer_upstream=numer_upstream,
        c_in=c_in,
        c_out=c_out,
    )


def cell_balance_terms_degenerate(
    st: "StreamingTerms",
    total_xs: np.ndarray,
    upstream_state: "UpstreamState",
) -> CellBalanceTerms:
    r"""Cylindrical pure-azimuthal degenerate variant of :func:`cell_balance_terms`.

    Used when :math:`|\mu| < 10^{-15}` — the level's axial cosine
    :math:`|\mu_z| \to 1` so the radial cosine :math:`|\eta| \to 0` and
    no radial face flow occurs.  The :math:`2|\mu|\,A_{\text{down}}` and
    :math:`|\mu|\,(A_{\text{in}}+A_{\text{out}})\,\psi^s_{\text{in}}`
    terms vanish; only the angular-redistribution + collision survive.
    """
    assert st.delta_A_over_w is not None
    assert st.alpha_in is not None
    assert st.alpha_out is not None
    assert st.tau_mm is not None
    assert st.volume is not None
    assert upstream_state.angular_upstream is not None, (
        "Cylindrical-degenerate cell balance requires an upstream "
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

    return CellBalanceTerms(
        denom=denom,
        numer_upstream=numer_upstream,
        c_in=c_in,
        c_out=c_out,
    )


__all__ = [
    "CellBalanceTerms",
    "cell_balance_terms",
    "cell_balance_terms_degenerate",
]
