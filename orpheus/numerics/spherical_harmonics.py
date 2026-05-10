r"""Real spherical harmonics evaluator on a discrete direction set.

Computes the real spherical harmonics :math:`Y_\ell^m(\hat\Omega_n)` at
every ordinate of a discrete angular direction set, returning the full
:math:`(N,\ \ell+1,\ 2\ell+1)` table consumed by the Galerkin
moment-projection / reconstruction pair (:math:`M = Y^* W`,
:math:`R = Y` with addition-theorem factor :math:`(2\ell+1)`).

Convention
----------

The project uses the **no-:math:`4\pi/(2\ell+1)`-prefactor** normalisation
of the real spherical harmonics, in which the addition theorem reads

.. math::
   :label: real-sh-addition-theorem

   \sum_{m=-\ell}^{\ell} Y_\ell^m(\hat\Omega)\,Y_\ell^m(\hat\Omega')
   \;=\; P_\ell(\hat\Omega \cdot \hat\Omega'),

so the Pℓ scattering reconstruction takes the form

.. math::

   q_n \;=\; \sum_{\ell=0}^{L} (2\ell+1) \sum_m Y_\ell^m(\hat\Omega_n)\,
            \phi^{\ell m}.

The polar axis is :math:`\mu_x` (so :math:`\cos\theta = \mu_x`,
:math:`\sin\theta = \sqrt{1-\mu_x^2}`); azimuth is measured in the
:math:`(\mu_y,\mu_z)` plane:

.. math::

   \cos\phi \;=\; \frac{\mu_y}{\sin\theta}, \qquad
   \sin\phi \;=\; \frac{\mu_z}{\sin\theta}.

The :math:`\ell\le 1` branch is hard-coded to bit-identical values for
the legacy :math:`P_0/P_1` regression tests:

.. math::

   Y_0^0 = 1, \quad Y_1^{-1} = \mu_z, \quad Y_1^0 = \mu_x,
   \quad Y_1^{+1} = \mu_y.

For :math:`\ell \ge 2` the formula uses :func:`scipy.special.lpmv` with
the Condon–Shortley phase :math:`(-1)^m` removed and norm
:math:`\sqrt{2(\ell-m)!/(\ell+m)!}` for :math:`m \ne 0`.

Cross-method consumers
----------------------

This module is generic infrastructure: every method that integrates an
angular field against the spherical-harmonic basis consumes it.

* **SN scattering** (:mod:`orpheus.sn.scattering`) — the
  :math:`Y^* W` projection that builds the per-ordinate Pℓ source.
* **PN solver** (future, §10 of the architecture report) —
  spherical-harmonic moment basis is the native space.
* **MC adjoint moments** — variance reduction with response moments
  built against :math:`Y_\ell^m`.
* **Energy-condensation diagnostics** — when within-group anisotropy
  needs an angular characterisation.

See Also
--------

* :ref:`pn-scattering` — the addition-theorem reconstruction in the
  SN solver.
* :file:`student_resources/02_spherical_harmonics.py` — pedagogical
  surface-plot visualisation of a single :math:`Y_\ell^m`. Shares the
  same :func:`scipy.special.lpmv` machinery and norm
  :math:`\sqrt{2(\ell-|m|)!/(\ell+|m|)!}`.

References
----------

* Bell, G. I. and Glasstone, S. (1970). *Nuclear Reactor Theory*.
  Van Nostrand Reinhold. §1.6 (real spherical harmonics in transport).
* Lewis, E. E. and Miller, W. F. Jr. (1993). *Computational Methods
  of Neutron Transport*. ANS. §4.7 (Pℓ scattering Galerkin
  reconstruction with the :math:`(2\ell+1)` factor).
"""

from __future__ import annotations

from math import factorial, sqrt

import numpy as np
from scipy.special import lpmv


__all__ = ["evaluate_real_sh"]


def evaluate_real_sh(
    L: int,
    mu_x: np.ndarray,
    mu_y: np.ndarray,
    mu_z: np.ndarray,
) -> np.ndarray:
    r"""Evaluate real spherical harmonics on a discrete direction set.

    Returns the full :math:`(N,\ L+1,\ 2L+1)` table of
    :math:`Y_\ell^m(\hat\Omega_n)` values used by the Galerkin
    projection / reconstruction.

    Parameters
    ----------
    L : int
        Maximum spherical-harmonic order. ``L == 0`` returns the P0
        table only; ``L < 0`` returns an all-zero array of zero shape.
    mu_x, mu_y, mu_z : np.ndarray, shape (N,)
        Direction cosines of the ordinates :math:`\hat\Omega_n =
        (\mu_x, \mu_y, \mu_z)`. The polar axis is :math:`\mu_x` (so
        :math:`\mu_x = \cos\theta`); azimuth is measured in the
        :math:`(\mu_y, \mu_z)` plane.

    Returns
    -------
    np.ndarray, shape (N, L+1, 2L+1)
        ``Y[n, l, l+m]`` is :math:`Y_\ell^m(\hat\Omega_n)` under the
        no-:math:`4\pi/(2\ell+1)`-prefactor convention; the
        :math:`m`-axis indexing is shifted by :math:`\ell` so that
        the slice ``Y[n, l, l-l : l+l+1]`` holds the :math:`2\ell+1`
        entries for order :math:`\ell`.

    Notes
    -----

    Implementation choices kept bit-identical to the legacy
    ``orpheus.sn.quadrature._build_spherical_harmonics`` (now a thin
    re-export of this function) so the regression snapshots at
    ``tests/sn/regression/snapshots/`` continue to pass.

    See the module docstring for the convention and citations.
    """
    N = len(mu_x)
    if L < 0:
        return np.zeros((N, 0, 0))
    Y = np.zeros((N, L + 1, 2 * L + 1))

    Y[:, 0, 0] = 1.0
    if L == 0:
        return Y

    Y[:, 1, 0] = mu_z   # m = -1
    Y[:, 1, 1] = mu_x   # m =  0
    Y[:, 1, 2] = mu_y   # m = +1
    if L == 1:
        return Y

    cos_theta = mu_x
    sin_theta = np.sqrt(np.maximum(1.0 - cos_theta * cos_theta, 0.0))
    on_axis = sin_theta < 1e-15
    safe_st = np.where(on_axis, 1.0, sin_theta)
    cos_phi = np.where(on_axis, 1.0, mu_y / safe_st)
    sin_phi = np.where(on_axis, 0.0, mu_z / safe_st)
    phi = np.arctan2(sin_phi, cos_phi)

    for l in range(2, L + 1):
        Y[:, l, l] = lpmv(0, l, cos_theta)  # m = 0: P_l(μ_x)
        for m in range(1, l + 1):
            P_lm = lpmv(m, l, cos_theta)
            sign = (-1.0) ** m   # remove Condon–Shortley phase
            norm = sqrt(2.0 * factorial(l - m) / factorial(l + m))
            cos_mphi = np.cos(m * phi)
            sin_mphi = np.sin(m * phi)
            Y[:, l, l + m] = sign * norm * P_lm * cos_mphi   # m > 0
            Y[:, l, l - m] = sign * norm * P_lm * sin_mphi   # m < 0
    return Y
