r"""SymPy derivation of the canonical 3-D :math:`G_{\rm bc}` kernel for
the cylinder white-BC.

The 3-D :math:`G_{\rm bc}^{\rm cyl}` kernel is obtained by integrating
the 3-D point-kernel
:math:`g(r, r') = e^{-\Sigma_t |r-r'|}/(4\pi |r-r'|^2)` over the
infinite z-extent of the cylinder lateral surface. The closed-form
result is

.. math::
   :label: peierls-cyl-Gbc-3d-derivation

   G_{\rm bc}^{\rm cyl}(r) \;=\; \frac{4 R}{\pi}\!\int_0^\pi
       \frac{R - r \cos\phi}{d^2(\phi)}\,
       \mathrm{Ki}_2\!\bigl(\Sigma_t\,d(\phi)\bigr)\,\mathrm d\phi,

with :math:`d(\phi) = \sqrt{r^2 - 2 r R \cos\phi + R^2}` the in-plane
chord and :math:`(R - r\cos\phi)/d` the Lambertian incidence-cosine
projection.

Setup
-----

Observer at interior point :math:`(r, 0, 0)`. Surface point at
:math:`(R\cos\phi, R\sin\phi, z)` for :math:`\phi \in [0, 2\pi)`,
:math:`z \in (-\infty, +\infty)`. Outward normal at the surface:
:math:`\hat n_b = (\cos\phi, \sin\phi, 0)`.

In-plane distance: :math:`d(\phi) = \sqrt{r^2 - 2 r R \cos\phi + R^2}`.
3-D distance: :math:`d_{\rm 3D} = \sqrt{d^2 + z^2}`.

For Mark closure, :math:`\psi^{-}(r_b, \Omega) = J^{-}/\pi` for all
incoming directions at any surface point. Cosine of the incidence
angle:

.. math::

   \cos\theta_{\rm inc}
   \;=\; (r_b - r_i)\cdot \hat n_b / d_{\rm 3D}
   \;=\; (R - r\cos\phi) / d_{\rm 3D}.

Therefore:

.. math::

   G_{\rm bc}^{\rm cyl}(r) \;=\; \frac{R}{\pi}\!\int_0^{2\pi}\!\!
       \mathrm d\phi\,(R - r\cos\phi)\!\int_{-\infty}^{+\infty}
       \frac{e^{-\Sigma_t\,d_{\rm 3D}}}{d_{\rm 3D}^3}\,\mathrm dz.

The :math:`z`-integral reduces to a Bickley function via the standard
substitution :math:`z = d \tan\alpha` (so :math:`d_{\rm 3D} = d/\cos\alpha`,
:math:`\mathrm dz = (d/\cos^2\alpha)\,\mathrm d\alpha`):

.. math::

   \int_{-\infty}^{+\infty}\frac{e^{-\Sigma_t\,d_{\rm 3D}}}
       {d_{\rm 3D}^3}\,\mathrm dz
   \;=\; \frac{1}{d^2}\!\int_{-\infty}^{+\infty}\!\cos\alpha\,
        e^{-\Sigma_t\,d/\cos\alpha}\,\mathrm d\alpha
   \;=\; \frac{2}{d^2}\,\mathrm{Ki}_2(\Sigma_t\,d)

(by the Bickley convention :math:`\mathrm{Ki}_n(x) = \int_0^{\pi/2}
\cos^{n-1}\theta\,e^{-x/\cos\theta}\,\mathrm d\theta` and the
:math:`\alpha \to -\alpha` symmetry that gives the factor of 2).

Symmetry around :math:`\phi = 0` collapses
:math:`\int_0^{2\pi}\to 2\!\int_0^\pi`, yielding
:eq:`peierls-cyl-Gbc-3d-derivation`.

The historical (incorrect) form
-------------------------------

Issue #112 Phase C. The legacy ``compute_G_bc`` cylinder branch used
the surface-centric form

.. math::

   G_{\rm bc}^{\rm legacy}(r) \;=\; \frac{2 R}{\pi}\!\int_0^\pi
       \frac{\mathrm{Ki}_1(\Sigma_t\,d(\phi))}{d(\phi)}\,\mathrm d\phi,

which carries three independent bugs vs the correct form: the Bickley
order (:math:`\mathrm{Ki}_1` instead of :math:`\mathrm{Ki}_2`), the
geometry factor (:math:`1/d` instead of the Lambertian
:math:`(R - r\cos\phi)/d^2`), and the leading coefficient
(:math:`2R/\pi` instead of :math:`4R/\pi`). The shipped
:func:`compute_G_bc_cylinder_3d` in :mod:`peierls_geometry` uses the
correct form (:eq:`peierls-cyl-Gbc-3d-derivation`).

Verification
------------

Pinned numerically against the production
:func:`compute_G_bc_cylinder_3d` and against the legacy 3-bug
form by ``tests/derivations/test_peierls_cylinder_g_bc_3d_symbolic.py``.
The thin-cell limits :math:`r \to 0,\,\Sigma_t R \to 0` give the
sanity values :math:`G_{\rm correct}(0) = 4\,\mathrm{Ki}_2(0) = 4`
and :math:`G_{\rm legacy}(0) = 2\,\mathrm{Ki}_1(0) = \pi`, a
:math:`\pi/4 \approx 21\,\%` overestimate at the optically-thin
end. See
``docs/theory/peierls_unified.rst:peierls-cyl-Gbc-3d-final``
(Issue #112 Phase C) for the cylinder-1G/1R row-sum probe table
and the Hébert convergence improvement that motivated the fix.
"""

from __future__ import annotations

import numpy as np

from ....common.kernels import ki_n_mp


def Ki_n(n: int, x: float, dps: int = 25) -> float:
    """Bickley-Naylor wrapper using ORPHEUS robust :func:`ki_n_mp`."""
    if x < 0:
        x = 0.0
    return float(ki_n_mp(n, x, dps))


def G_bc_cyl_correct(
    r: float,
    R: float,
    Sigma_t: float,
    n_quad: int = 64,
    dps: int = 25,
) -> float:
    r"""Reference numerical evaluator for the **correct** 3-D
    :math:`G_{\rm bc}^{\rm cyl}` kernel
    (:eq:`peierls-cyl-Gbc-3d-derivation`).

    Pure-Python Gauss-Legendre quadrature on
    :math:`(4R/\pi)\!\int_0^\pi (R - r\cos\phi)/d^2 \cdot
    \mathrm{Ki}_2(\Sigma_t\,d)\,\mathrm d\phi`. The production
    :func:`compute_G_bc_cylinder_3d` should agree with this to
    quadrature precision (verified by the symbolic test).
    """
    pts, wts = np.polynomial.legendre.leggauss(n_quad)
    phi_pts = 0.5 * (pts + 1) * np.pi
    phi_wts = wts * 0.5 * np.pi
    total = 0.0
    for k in range(n_quad):
        cf = float(np.cos(phi_pts[k]))
        d2 = r * r - 2 * r * R * cf + R * R
        d = float(np.sqrt(max(d2, 0.0)))
        if d <= 0.0:
            continue
        total += phi_wts[k] * (R - r * cf) / d2 * Ki_n(2, Sigma_t * d, dps=dps)
    return float(4.0 * R / np.pi * total)


def G_bc_cyl_legacy_three_bug(
    r: float,
    R: float,
    Sigma_t: float,
    n_quad: int = 64,
    dps: int = 25,
) -> float:
    r"""Reference numerical evaluator for the **legacy 3-bug**
    :math:`G_{\rm bc}^{\rm legacy,cyl}` form for regression-witness
    purposes.

    Reproduces the historical Issue #112 form
    :math:`(2R/\pi)\!\int_0^\pi \mathrm{Ki}_1(\Sigma_t d)/d\,\mathrm d\phi`
    so the ratio current/correct can be quantified by the test
    suite. **Do not** call this from production code — it carries the
    three independent bugs (Bickley order, geometry factor, leading
    coefficient) that the correct form fixes.
    """
    pts, wts = np.polynomial.legendre.leggauss(n_quad)
    phi_pts = 0.5 * (pts + 1) * np.pi
    phi_wts = wts * 0.5 * np.pi
    total = 0.0
    for k in range(n_quad):
        cf = float(np.cos(phi_pts[k]))
        d2 = r * r - 2 * r * R * cf + R * R
        d = float(np.sqrt(max(d2, 0.0)))
        if d <= 0.0:
            continue
        total += phi_wts[k] * Ki_n(1, Sigma_t * d, dps=dps) / d
    return float(2.0 * R / np.pi * total)
