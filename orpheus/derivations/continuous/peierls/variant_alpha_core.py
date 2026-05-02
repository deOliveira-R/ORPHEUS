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
- :mod:`orpheus.derivations.continuous.peierls.greens_function`
  (sphere) and
  :mod:`orpheus.derivations.continuous.peierls.greens_function_cylinder`
  (cylinder) — call sites.
- Phase-2 unification memo:
  ``.claude/agent-memory/method-implementer/cylinder_variant_alpha_phase2_unification.md``.
- Cross-domain frame memo:
  ``.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md``.
"""
from __future__ import annotations

import numpy as np


def compute_resolvent_T(alpha: float, tau_period: float) -> float:
    r"""Rank-1 boundary-to-boundary scattering resolvent.

    .. math::

        T(\alpha, \tau_{\rm period}) \;=\; \frac{1}{1 - \alpha
            \, e^{-\tau_{\rm period}}}.

    This is the rank-1 case of the operator resolvent :math:`T =
    (I - S)^{-1}` with :math:`S = \alpha\,R_{\rm chord}` where
    :math:`R_{\rm chord} = e^{-\tau_{\rm period}}` is the
    chord-attenuation kernel along one bounce period.

    Parameters
    ----------
    alpha : float
        Surface reflectivity in :math:`[0, 1]`.
    tau_period : float
        Optical depth :math:`\Sigma_t \cdot L_{\rm period}` of one
        bounce period along the trajectory. Geometry-specific:

        - Sphere: :math:`\tau_{\rm period} = \Sigma_t \cdot 2 R \mu_{\rm surf}`
        - Cylinder: :math:`\tau_{\rm period} = \Sigma_t \cdot 2\sqrt{R^2
          - b^2} / \sqrt{1 - \mu_{\rm axial}^2}`

    Returns
    -------
    float
        The resolvent :math:`T \ge 1`. At :math:`\alpha = 0` returns
        :math:`1` (no boundary feedback). At :math:`\alpha = 1` and
        small :math:`\tau_{\rm period}` diverges as :math:`1/\tau_{\rm
        period}` (offset by :math:`B \to 0` in the calling closure —
        see V_α1).
    """
    return 1.0 / (1.0 - alpha * np.exp(-tau_period))


def apply_variant_alpha_closure(
    F: float | np.ndarray,
    B: float | np.ndarray,
    tau_first_leg: float | np.ndarray,
    tau_period: float | np.ndarray,
    alpha: float,
) -> float | np.ndarray:
    r"""Variant α operator closure — combine first-leg and bounce-period
    integrals into the new angular flux.

    Implements:

    .. math::

        \psi_{\rm new} \;=\; F + e^{-\tau_{\rm first\,leg}}\,
            \alpha\,B\,T(\alpha, \tau_{\rm period}),

    where :math:`F`, :math:`B`, :math:`\tau_{\rm first\,leg}`, and
    :math:`\tau_{\rm period}` are computed by the calling geometry's
    trajectory machinery, and :math:`T` is :func:`compute_resolvent_T`.
    At :math:`\alpha = 0` the second term vanishes identically.

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
        Surface reflectivity in :math:`[0, 1]`.

    Returns
    -------
    float or ndarray
        The new angular flux :math:`\psi_{\rm new}` with the same
        shape as the inputs.
    """
    if alpha == 0.0:
        return F
    T = compute_resolvent_T(alpha, tau_period)
    psi_surf = alpha * B * T
    return F + np.exp(-tau_first_leg) * psi_surf
