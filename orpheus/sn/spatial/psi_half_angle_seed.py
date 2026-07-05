r"""The Hébert §3.9.4 starting-direction DD march — the direct ψ½ solver.

This module hosts :func:`carlson_inward_sweep_from_source`, the Hébert
(3.434)–(3.435) diamond-difference march along the starting direction
:math:`\mu = \mp 1` — the ENGINE of the #282 route-(a) direct
starting-direction treatment (#280 Phase 2.5d).

The starting-direction transport equation
=========================================

At the angular edge :math:`\mu_{\rm start}` of a curvilinear μ-level the
angular redistribution coefficient :math:`1 - \mu^2` vanishes
(:math:`\alpha_{1/2} = 0`, Hébert Eq. 3.423), so the streaming–collision
balance DECOUPLES from the α-cascade and reduces to a plain DD recurrence
in radius:

.. math::
   :label: hebert-3-434

   \bar\phi_i \;=\; \frac{\Delta r_i \cdot \bar Q_i
                            + 2 \cdot \bar\phi_{i+1/2}}
                          {\Delta r_i \cdot \Sigma_i + 2}
   \quad\text{Hébert Eq. (3.434)}

.. math::
   :label: hebert-3-435

   \bar\phi_{i-1/2} \;=\; 2 \cdot \bar\phi_i - \bar\phi_{i+1/2}
   \quad\text{Hébert Eq. (3.435)}

with the cell-averaged source at the starting direction folded from the
Legendre moments via :math:`P_\ell(\pm 1) = (\pm 1)^\ell`:

.. math::
   :label: hebert-3-432-source

   \bar Q_i \;=\; \sum_\ell \frac{2\ell+1}{2}\,Q_\ell(r_i)\,(\pm 1)^\ell

(the fold lives in
:func:`~orpheus.numerics.spaces.starting_direction_space.fold_moments_to_starting_direction`
— the R14 helper the operator seed arms and the source factories share).

Development history — the seed-strategy zoo (retired 2026-07-04, 2.5d d3)
=========================================================================

The recurrence seed :math:`\psi_{1/2,i}` was historically produced by a
swappable ``PsiHalfAngleSeed`` strategy family on the M-M closure:

* ``ZeroSeed`` — the pre-ERR-026 hardcoded zero (wrong term
  initialization, Mode 3);
* ``CarlsonInwardSweep`` — this module's march driven by the PROXY
  source :math:`\bar Q = \Sigma_t\phi_0/\!\sum w`, exact only at the
  flat-flux equilibrium (ERR-058 b — O(1) wrong off equilibrium, floored
  the curvilinear MMS at ~0.04 L2 independent of mesh, Issue #195);
* ``AngularEdgeExtrapolation`` — the operator-consistent 2-point angular
  extrapolation of the ITERATE, which fixed the matvec (#195) but left
  the SOLVE seeded from the PREVIOUS iterate — the #282 seed LAG (sphere
  cold residual 5.18e5, seed-sensitivity 4.57e-2).

Route (a) (#282, ruling R10) retired the whole strategy family: the
starting-direction flux is now first-class STATE (the
``StartingDirectionField`` block of the composite, present per level
under the R12a predicate), the SOLVE marches it directly from the TRUE
q½ source through this function, and the APPLY reads the given carrier
block.  On the non-carrying cylinder levels (R12a: product rules
τ_raw = 0, level-symmetric rules τ_raw = 1) the closure inlines the
2-point angular-edge extrapolation — see
:meth:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep.edge_extrapolated_seed`.

References
==========

* Hébert, A. (2009).  *Applied Reactor Physics*.  §3.9.4 pp. 141-144,
  Eqs. (3.432)-(3.435).  Local copy:
  ``scratch/literature/Hebert(2009)Chapter3.pdf``.
* Pomraning, G.C. (1989).  *The Transport Equation in General
  Geometry*.  NSE 101:330-340 p. 339.  Cross-reference for the
  structural singularity at r = 0 in curvilinear geometry.
* Literature memo: ``.claude/agent-memory/literature-researcher/phase_d_carlson_coupled_pole.md``.
"""

from __future__ import annotations

import numpy as np


def carlson_inward_sweep_from_source(
    Q_bar: np.ndarray,
    sigma_t: np.ndarray,
    dr: np.ndarray,
    bc_outer_value: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Run the Hébert §3.9.4 (3.434)-(3.435) starting-direction DD march.

    The direct solver of the :math:`\mu = -1` starting-direction ODE —
    marches INWARD from ``i = nx-1`` down to ``i = 0``, entering each
    cell through its OUTER face.  The same recurrence marches the
    OUTWARD :math:`\mu = +1` leg by reversing the cell axis of every
    input and un-reversing the output (orientation is carried by the
    DATA, never a flag — the 2.5a discipline):

    .. code-block:: python

       cells_rev, face_out = carlson_inward_sweep_from_source(
           Q_plus[:, ::-1], sigma_t[:, ::-1], dr[::-1], pole_face,
       )
       cells_plus = cells_rev[:, ::-1]      # outward-marched ψ½⁺ cells
       # face_out is the outer-face (r = R) value — the outflow corner.

    Parameters
    ----------
    Q_bar : np.ndarray, shape ``(ng, nx)``
        Cell-averaged TRUE source at the starting direction — the
        Legendre fold :eq:`hebert-3-432-source` of the within-group
        source moments (``fold_moments_to_starting_direction``).  For an
        isotropic source this is :math:`\tfrac12 \bar Q_{\rm iso}`.
    sigma_t : np.ndarray, shape ``(ng, nx)``
        Cell-centred total cross-section, per-group, per-cell.
    dr : np.ndarray, shape ``(nx,)``
        Radial cell widths.
    bc_outer_value : np.ndarray, shape ``(ng,)``
        The ENTRY face value — for the inward leg the outer-face
        (r = R) angular flux at :math:`\mu = -1` (the inflow corner:
        0 for vacuum, the reflected outflow corner for reflective);
        for the reversed outward call the pole-continuation value
        :math:`\psi^+_{1/2}(0) = \psi^-_{1/2}(0)`.

    Returns
    -------
    phi_cells : np.ndarray, shape ``(ng, nx)``
        Cell-centred starting-direction flux
        :math:`\bar\phi_i \equiv \psi_{1/2,i}` (the M-M recurrence seed
        on the inward leg).
    phi_face_final : np.ndarray, shape ``(ng,)``
        The EXIT face value after the last marched cell — the r = 0
        pole face on the inward leg (the pole-continuation datum the
        outward leg starts from), the r = R outer face on the reversed
        outward call (the outflow corner).

    Notes
    -----
    The recurrence is sequential in cells (each step depends on its
    upwind neighbour's face value) and vectorised across groups via
    NumPy broadcasting.
    """
    ng, nx = Q_bar.shape
    phi_aux = np.zeros((ng, nx), dtype=Q_bar.dtype)
    phi_face = bc_outer_value.copy()  # (ng,) — entry face at i = nx
    for k in range(nx - 1, -1, -1):
        denom = dr[k] * sigma_t[:, k] + 2.0       # (ng,)
        phi_cell = (dr[k] * Q_bar[:, k] + 2.0 * phi_face) / denom
        phi_aux[:, k] = phi_cell
        # Hébert (3.435): step to the next (downwind) face.
        phi_face = 2.0 * phi_cell - phi_face
    return phi_aux, phi_face


def carlson_inward_sweep_transpose(
    cells_bar: np.ndarray,
    final_face_bar: np.ndarray,
    sigma_t: np.ndarray,
    dr: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Euclidean adjoint of :func:`carlson_inward_sweep_from_source`.

    The forward march :math:`(\bar Q, f_{\rm entry}) \mapsto (\bar\phi_{\rm cells},
    f_{\rm exit})` is the linear inward DD recurrence (Hébert 3.434/3.435).  This
    is its transpose — the reverse-mode adjoint (#280 Phase 2.5b): given
    cotangents on the marched cells (:math:`\bar\phi_{\rm cells}`-cotangent
    ``cells_bar``) and the exit face (``final_face_bar``), return the cotangents
    on the source (``Q_bar``) and the entry face (``bc_outer_bar``).

    The forward marches ``k = nx-1 → 0`` (inward); the adjoint retraces
    ``k = 0 → nx-1`` (outward), threading the running face cotangent ``f_bar``
    from the exit face back to the entry face — the seed-block sibling of the
    bulk reverse-scan :func:`ordinate_scan_transpose`.

    Parameters
    ----------
    cells_bar : ndarray, shape ``(ng, nx)``
        Cotangent on the marched cell fluxes :math:`\bar\phi_i`.
    final_face_bar : ndarray, shape ``(ng,)``
        Cotangent on the exit (``phi_face_final``) face.
    sigma_t, dr : ndarray
        The SAME ``(ng, nx)`` cross-section and ``(nx,)`` widths the forward
        march consumed.

    Returns
    -------
    Q_bar : ndarray, shape ``(ng, nx)``
        Cotangent on the cell-averaged source ``Q_bar`` (the forward's input).
    bc_outer_bar : ndarray, shape ``(ng,)``
        Cotangent on the entry face ``bc_outer_value``.
    """
    ng, nx = cells_bar.shape
    Q_bar = np.zeros((ng, nx), dtype=cells_bar.dtype)
    f_bar = final_face_bar.copy()               # cotangent on f[0] = exit face
    for k in range(nx):
        denom = dr[k] * sigma_t[:, k] + 2.0     # (ng,)
        # forward: phi_cell[k] = (dr[k]·Q[k] + 2·f_in)/denom; f_out = 2·phi_cell − f_in.
        # f_bar here is the cotangent on f_out (= the face AFTER cell k).
        c_bar = cells_bar[:, k] + 2.0 * f_bar
        f_in_bar = -f_bar
        Q_bar[:, k] = (dr[k] / denom) * c_bar
        f_in_bar = f_in_bar + (2.0 / denom) * c_bar
        f_bar = f_in_bar                        # cotangent on f_in = f[k+1], next step
    return Q_bar, f_bar


__all__ = [
    "carlson_inward_sweep_from_source",
    "carlson_inward_sweep_transpose",
]
