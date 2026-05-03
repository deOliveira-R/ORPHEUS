r"""Phase-2 unification core — Variant α boundary-to-boundary scattering
resolvent and operator-application primitives shared across geometries.

This module institutionalises the **cross-domain frame** validated in
``.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md``:
sphere and cylinder Variant α solvers both close their boundary
trajectories with the same rank-1 boundary-to-boundary scattering
resolvent

.. math::

    T(\alpha, \tau_{\rm period}) \;=\; \frac{1}{1 - \alpha
        \, e^{-\tau_{\rm period}}}

(where :math:`\tau_{\rm period} = \Sigma_t\,L_{\rm period}` is the
optical depth of one bounce period). This is the rank-1 case of the
operator resolvent :math:`T = (I - S)^{-1}` with :math:`S = \alpha
R_{\rm chord}` for a 2-surface scattering operator. Sphere and cylinder
differ only in their geometry-specific evaluation of :math:`L_{\rm
period}` (sphere: :math:`2 R \mu_{\rm surf}`; cylinder:
:math:`2\sqrt{R^2 - b^2} / \sqrt{1 - \mu_{\rm axial}^2}`).

The algebraic closure is identical:

.. math::

    \psi_{\rm new} \;=\; F + e^{-\tau_{\rm first\,leg}} \cdot
        \alpha B \cdot T(\alpha, \tau_{\rm period}),

where :math:`F` and :math:`B` are the (geometry-specific) first-leg
and bounce-period source-line integrals. This module exposes the two
operations — the resolvent and the closure transform — so each
geometry's solver provides only the trajectory primitives.

Architecture
------------

The Phase-2 minimum-viable unification keeps **pure functions only** —
no ABC, no mixin, no dataclass-with-inheritance. Module-level
primitives + per-geometry call sites are sufficient for two geometric
instances and easier to refactor later if a third geometry comes
(slab, hollow sphere, annulus). See
``.claude/agent-memory/method-implementer/cylinder_variant_alpha_phase2_unification.md``
for the rationale and the ``CurvilinearGeometry`` mounting trade-off
that was rejected.

Rank-2 extension (Phase 3B)
---------------------------

For two-surface geometries with **independent reflectivities on each
surface** (asymmetric slab :math:`\alpha_L \ne \alpha_R`, hollow
sphere, annulus), the rank-1 scalar resolvent is replaced by a
:math:`2 \times 2` matrix resolvent
:func:`compute_resolvent_T_rank2` and the closure becomes
:func:`apply_variant_alpha_closure_rank2`. The rank-2 monodromy is
:math:`S = \mathrm{antidiag}(\alpha_L\,e^{-\tau}, \alpha_R\,e^{-\tau})`
with single-transit optical depth :math:`\tau = \Sigma_t L /|\mu|`
(NOT the full out+back period — the rank-2 resolvent encodes the
two-bounce pattern via matrix structure, not via :math:`\alpha^2` in a
scalar formula). At :math:`\alpha_L = \alpha_R` the rank-2 closure is
mathematically equivalent to the rank-1 symmetric closure (verified
numerically to ≤ 1e-10 in the Phase 3B test suite); the two arithmetic
forms are NOT bit-equal because IEEE-754 evaluation order differs.

What is NOT unified
-------------------

- Geometry-specific phase-space discretisation: sphere uses
  :math:`(r, \mu)`; cylinder uses :math:`(r, \mu_{\rm axial},
  \varphi_{\rm az})`. No shared abstraction would pay for itself
  here.
- Geometry-specific chord arithmetic (``_first_leg_2d_chord``,
  ``_bounce_period_2d_chord``, ``_trajectory_segments``,
  ``_chord_segments``). These are linked to the shape's parametric
  surface and don't generalise without forcing a heavyweight
  ``CurvilinearGeometry`` abstraction that breaks the rank-1
  Variant α structure.
- Multi-region piecewise-:math:`\Sigma_t` segmentation (sphere only
  for now). The resolvent and closure transform are still applied
  at the end, but the F and B integrals are computed by composing
  per-region segments — geometry-specific.

References
----------

- Sanchez, R. (1986). *Transp. Theor. Stat. Phys.* 14.
- :mod:`orpheus.derivations.continuous.trajectory_resolvent.greens_function`
  (sphere) and
  :mod:`orpheus.derivations.continuous.trajectory_resolvent.greens_function_cylinder`
  (cylinder) — call sites.
- Phase-2 unification memo:
  ``.claude/agent-memory/method-implementer/cylinder_variant_alpha_phase2_unification.md``.
- Cross-domain frame memo:
  ``.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md``.
"""
from __future__ import annotations

import numpy as np


def compute_resolvent_T(
    alpha_per_period: float, tau_period: float,
) -> float:
    r"""Rank-1 boundary-to-boundary scattering resolvent.

    .. math::

        T(\alpha_{\rm per\,period}, \tau_{\rm period}) \;=\;
            \frac{1}{1 - \alpha_{\rm per\,period}
            \, e^{-\tau_{\rm period}}}.

    This is the rank-1 case of the operator resolvent :math:`T =
    (I - S)^{-1}` with :math:`S = \alpha_{\rm per\,period}\,R_{\rm chord}`
    where :math:`R_{\rm chord} = e^{-\tau_{\rm period}}` is the chord-
    attenuation kernel along one bounce period.

    The first argument is the **per-period reflection product** —
    the cumulative reflection amplitude accumulated by ONE FULL BOUNCE
    PERIOD of the trajectory. For 1-bounce-per-period geometries
    (sphere, cylinder) this is the BC reflectivity itself,
    :math:`\alpha`. For 2-bounce-per-period geometries (slab symmetric
    rank-1), one full period samples both walls so the reflection
    product is :math:`\alpha^2` — see :func:`apply_variant_alpha_closure`
    for the back-compatible call-site contract.

    Parameters
    ----------
    alpha_per_period : float
        Per-period reflection product. Must lie in :math:`[0, 1]`.
        Sphere/cylinder: :math:`\alpha`. Slab symmetric: :math:`\alpha^2`.
    tau_period : float
        Optical depth :math:`\Sigma_t \cdot L_{\rm period}` of one
        bounce period along the trajectory. Geometry-specific:

        - Sphere: :math:`\tau_{\rm period} = \Sigma_t \cdot 2 R \mu_{\rm surf}`
        - Cylinder: :math:`\tau_{\rm period} = \Sigma_t \cdot 2\sqrt{R^2
          - b^2} / \sqrt{1 - \mu_{\rm axial}^2}`
        - Slab symmetric: :math:`\tau_{\rm period} = \Sigma_t \cdot 2L /
          |\mu|` (full transit out + reverse transit back).

    Returns
    -------
    float
        The resolvent :math:`T \ge 1`. At :math:`\alpha_{\rm per\,period}
        = 0` returns :math:`1` (no boundary feedback). At unit per-
        period reflection product and small :math:`\tau_{\rm period}`
        diverges as :math:`1/\tau_{\rm period}` (offset by :math:`B \to
        0` in the calling closure — see V_α1).
    """
    return 1.0 / (1.0 - alpha_per_period * np.exp(-tau_period))


def compute_resolvent_T_rank2(
    alpha_left: float,
    alpha_right: float,
    tau_single_transit: float | np.ndarray,
) -> np.ndarray:
    r"""Rank-2 boundary-to-boundary scattering resolvent for asymmetric
    two-surface geometries.

    Returns the :math:`2 \times 2` matrix :math:`T = (I - S)^{-1}` for
    the asymmetric-slab monodromy

    .. math::

        S(\alpha_L, \alpha_R, \tau) = \begin{pmatrix}
            0                              & \alpha_L\,e^{-\tau} \\
            \alpha_R\,e^{-\tau}            & 0
        \end{pmatrix},

    yielding

    .. math::

        T(\alpha_L, \alpha_R, \tau)
            = \frac{1}{1 - \alpha_L\,\alpha_R\,e^{-2\tau}}
              \begin{pmatrix}
                  1                       & \alpha_L\,e^{-\tau} \\
                  \alpha_R\,e^{-\tau}     & 1
              \end{pmatrix}.

    Phase-space: the 2-vector surface state is :math:`[\psi_L^+, \psi_R^-]`
    (outgoing flux at :math:`x=0` in :math:`+\mu` direction, outgoing
    flux at :math:`x=L` in :math:`-\mu` direction). One step of the
    bouncing trajectory carries :math:`\psi_L^+` to :math:`\psi_R^-`
    (transit + reflection at right wall, factor :math:`\alpha_R\,
    e^{-\tau}`) and :math:`\psi_R^-` to :math:`\psi_L^+` (transit +
    reflection at left wall, factor :math:`\alpha_L\,e^{-\tau}`). The
    diagonal of :math:`S` is zero because one step never returns to
    the same wall.

    Special cases verified by direct algebra:

    - :math:`\alpha_L = \alpha_R = \alpha` (symmetric): ``T_11 = T_22 =
      1/(1-\alpha^2 e^{-2\tau})``, ``T_12 = T_21 = \alpha\,e^{-\tau}/
      (1-\alpha^2 e^{-2\tau})``. The leading-:math:`\alpha` factor in
      the closure (which lives outside :math:`T`) reduces this to the
      rank-1 :math:`\psi_{\rm surf} = \alpha B / (1 - \alpha^2
      e^{-2\tau})` exactly, so the rank-2 framework reproduces the
      Phase-3A symmetric slab when both walls share :math:`\alpha`.
    - :math:`\alpha_L = 0` or :math:`\alpha_R = 0` (one vacuum wall):
      ``det(I - S) = 1``, ``T_11 = T_22 = 1``, ``T_12 = \alpha_L
      e^{-\tau}``, ``T_21 = \alpha_R e^{-\tau}``. The vacuum wall
      simply zeroes the corresponding column of :math:`\alpha B` in
      the closure.
    - :math:`\alpha_L = \alpha_R = 0` (vacuum-vacuum): :math:`T = I`.

    Parameters
    ----------
    alpha_left, alpha_right : float
        Per-wall reflection amplitudes in :math:`[0, 1]`.
    tau_single_transit : float or ndarray
        Optical depth :math:`\Sigma_t \cdot L / |\mu|` of one
        single-wall-to-wall transit (NOT the full out+back period).
        Vectorisable.

    Returns
    -------
    ndarray
        Shape ``(2, 2)`` for scalar input, or ``(..., 2, 2)`` for
        array input — final two axes carry the matrix.
    """
    e_tau = np.exp(-np.asarray(tau_single_transit, dtype=float))
    det = 1.0 - alpha_left * alpha_right * e_tau * e_tau
    # Build the 2x2 (or batched) matrix.
    T = np.empty(e_tau.shape + (2, 2)) if e_tau.ndim > 0 else np.empty((2, 2))
    T[..., 0, 0] = 1.0
    T[..., 0, 1] = alpha_left * e_tau
    T[..., 1, 0] = alpha_right * e_tau
    T[..., 1, 1] = 1.0
    T = T / det[..., None, None] if e_tau.ndim > 0 else T / det
    return T


def apply_variant_alpha_closure_rank2(
    F: float | np.ndarray,
    B_RL: float | np.ndarray,
    B_LR: float | np.ndarray,
    tau_first_leg: float | np.ndarray,
    tau_single_transit: float | np.ndarray,
    alpha_left: float,
    alpha_right: float,
    *,
    surface: str,
) -> float | np.ndarray:
    r"""Rank-2 Variant α surface closure for the asymmetric slab.

    Combines first-leg (:math:`F`) and **single-transit** B integrals
    (:math:`B_{RL}` from :math:`x=L` toward :math:`x=0`,
    :math:`B_{LR}` from :math:`x=0` toward :math:`x=L`) into the new
    angular flux at the interior point, by attenuating the appropriate
    surface flux back along the first-leg chord.

    The closure equations (single-transit bounce):

    .. math::

       \psi_L^+(\mu) &= \alpha_L\,B_{LR}(\mu) + \alpha_L\,e^{-\tau}\,
                          \psi_R^-(\mu), \\
       \psi_R^-(\mu) &= \alpha_R\,B_{RL}(\mu) + \alpha_R\,e^{-\tau}\,
                          \psi_L^+(\mu),

    where :math:`B_{LR}` is the single-transit chord integral
    parameterized as :math:`\int_0^{L/|\mu|} q(|\mu|\,s)\,e^{-\Sigma_t
    s}\,\mathrm d s` (source contribution AT :math:`x=0` from particles
    moving in the :math:`-\hat\mu` direction; chord traverses :math:`x
    = 0 \to L`), and :math:`B_{RL}` is the mirror :math:`\int_0^{L/|\mu|}
    q(L - |\mu|\,s)\,e^{-\Sigma_t s}\,\mathrm d s` (source contribution
    AT :math:`x=L` from particles moving in the :math:`+\hat\mu`
    direction; chord traverses :math:`x = L \to 0`).

    The intuition: :math:`\psi_L^+ = \alpha_L\,\psi(0, -\mu)`. The
    incoming :math:`\psi(0, -\mu)` integrates source emitted along
    the chord from :math:`x=0` going INTO the slab (in :math:`+x`
    direction) — that is :math:`B_{LR}`. Symmetrically, :math:`\psi_R^-
    = \alpha_R\,\psi(L, +\mu)` involves source from :math:`x=L` going
    INTO the slab (in :math:`-x` direction) which is :math:`B_{RL}`.

    The interior reconstruction is

    .. math::

       \psi(x, \mu) = F(x, \mu) + e^{-\tau_{\rm first\,leg}}\,
                       \psi_{\rm surface}(\mu),

    where :math:`\psi_{\rm surface}` is :math:`\psi_L^+(\mu)` for
    :math:`\mu > 0` (the trajectory entered from the left wall) and
    :math:`\psi_R^-(\mu)` for :math:`\mu < 0` (entered from the right
    wall). The ``surface`` argument selects which one.

    Parameters
    ----------
    F : float or ndarray
        First-leg trajectory integral.
    B_RL : float or ndarray
        Single-transit chord integral with chord traversed from
        :math:`x=L` to :math:`x=0`, parametrised as
        :math:`\int_0^{L/|\mu|} q(L - |\mu|\,s)\,e^{-\Sigma_t s}\,
        \mathrm d s`. Used in the :math:`\psi_R^-` closure equation.
    B_LR : float or ndarray
        Single-transit chord integral with chord traversed from
        :math:`x=0` to :math:`x=L`, parametrised as
        :math:`\int_0^{L/|\mu|} q(|\mu|\,s)\,e^{-\Sigma_t s}\,
        \mathrm d s`. Used in the :math:`\psi_L^+` closure equation.
    tau_first_leg : float or ndarray
        Optical depth of the backward chord from :math:`(x, \mu)` to
        the first surface arrival.
    tau_single_transit : float or ndarray
        Optical depth :math:`\Sigma_t \cdot L / |\mu|` of one
        single-wall-to-wall transit.
    alpha_left, alpha_right : float
        Per-wall reflection amplitudes in :math:`[0, 1]`.
    surface : {'left', 'right'}
        Which surface flux feeds the interior reconstruction. ``left``
        for :math:`\mu > 0` (entered from :math:`x = 0`); ``right``
        for :math:`\mu < 0` (entered from :math:`x = L`).

    Returns
    -------
    float or ndarray
        New angular flux :math:`\psi_{\rm new}` at :math:`(x, \mu)`.
    """
    if surface not in ("left", "right"):
        raise ValueError(f"surface must be 'left' or 'right'; got {surface!r}")

    e_tau = np.exp(-np.asarray(tau_single_transit, dtype=float))
    det = 1.0 - alpha_left * alpha_right * e_tau * e_tau

    # Closure: with the (α B)-vector being [α_L · B_LR, α_R · B_RL]:
    # ψ_L^+ = T_11 · α_L · B_LR + T_12 · α_R · B_RL
    #       = (1/det) · [α_L · B_LR + α_L · e^{-τ} · α_R · B_RL]
    # ψ_R^- = T_21 · α_L · B_LR + T_22 · α_R · B_RL
    #       = (1/det) · [α_R · e^{-τ} · α_L · B_LR + α_R · B_RL]
    if surface == "left":
        psi_surf = (alpha_left * B_LR
                    + alpha_left * e_tau * alpha_right * B_RL) / det
    else:  # "right"
        psi_surf = (alpha_right * e_tau * alpha_left * B_LR
                    + alpha_right * B_RL) / det

    return F + np.exp(-tau_first_leg) * psi_surf


def apply_variant_alpha_closure(
    F: float | np.ndarray,
    B: float | np.ndarray,
    tau_first_leg: float | np.ndarray,
    tau_period: float | np.ndarray,
    alpha: float,
    *,
    alpha_per_period: float | None = None,
) -> float | np.ndarray:
    r"""Variant α operator closure — combine first-leg and bounce-period
    integrals into the new angular flux.

    Implements:

    .. math::

        \psi_{\rm new} \;=\; F + e^{-\tau_{\rm first\,leg}}\,
            \alpha\,B\,T(\alpha_{\rm per\,period}, \tau_{\rm period}),

    where :math:`F`, :math:`B`, :math:`\tau_{\rm first\,leg}`, and
    :math:`\tau_{\rm period}` are computed by the calling geometry's
    trajectory machinery, :math:`\alpha` is the per-bounce reflection
    amplitude (the BC reflectivity, applied once for the FIRST surface
    arrival), and :math:`T` is :func:`compute_resolvent_T` evaluated at
    the **per-period reflection product** — the cumulative reflection
    amplitude accumulated over ONE FULL BOUNCE PERIOD.

    Per-period reflection product (:math:`\alpha_{\rm per\,period}`)
    -----------------------------------------------------------------

    The distinction between ``alpha`` and ``alpha_per_period``
    formalises a structural difference between geometries:

    - **Sphere / cylinder rank-1** — 1 bounce per period (one
      reflection on the single closed surface per trajectory cycle).
      Per-period reflection product :math:`= \alpha`. Defaults
      reproduce existing behaviour at all sphere/cylinder call sites
      (no kwarg needed).
    - **Slab symmetric rank-1** — 2 bounces per period (one
      reflection on each wall per trajectory cycle). Per-period
      reflection product :math:`= \alpha^2`. Caller must pass
      ``alpha_per_period=alpha**2``.

    The leading factor :math:`\alpha` in :math:`\psi_{\rm surf} =
    \alpha B T` is geometry-independent — it is the reflection
    amplitude for the SINGLE first surface arrival (the bouncing
    geometric series only kicks in *after* that first reflection,
    and the resolvent :math:`T` carries the whole subsequent
    geometric series).

    Vectorised over the inputs — pass arrays of identical shape and
    receive an array of the same shape.

    Parameters
    ----------
    F : float or ndarray
        First-leg source-line integral
        :math:`\int_0^{L_0} q(r(s))\,e^{-\Sigma_t s}\,\mathrm{d}s` along
        the backward chord from :math:`(r, \Omega)` to surface entry.
    B : float or ndarray
        Bounce-period source-line integral
        :math:`\int_0^{L_{\rm period}} q(r_{\rm chord}(s))\,
        e^{-\Sigma_t s}\,\mathrm{d}s` along the antipodal chord at
        the conserved impact parameter.
    tau_first_leg : float or ndarray
        Optical depth :math:`\Sigma_t \cdot L_{\rm first\,leg}` from
        :math:`(r, \Omega)` to surface entry — used to attenuate
        :math:`\psi_{\rm surf}` back to the interior point.
    tau_period : float or ndarray
        Optical depth of one full bounce period.
    alpha : float
        Reflection amplitude for the FIRST surface arrival (i.e. the
        BC reflectivity). Lies in :math:`[0, 1]`.
    alpha_per_period : float, keyword-only, optional
        Per-period reflection product feeding the geometric resolvent
        :math:`T`. **Defaults to ``alpha``** for back-compatibility
        with sphere/cylinder (1-bounce-per-period geometries). Slab
        symmetric (2-bounce-per-period) callers must pass
        ``alpha_per_period=alpha**2``.

    Returns
    -------
    float or ndarray
        The new angular flux :math:`\psi_{\rm new}` with the same
        shape as the inputs.
    """
    if alpha_per_period is None:
        alpha_per_period = alpha
    if alpha == 0.0:
        return F
    T = compute_resolvent_T(alpha_per_period, tau_period)
    psi_surf = alpha * B * T
    return F + np.exp(-tau_first_leg) * psi_surf
