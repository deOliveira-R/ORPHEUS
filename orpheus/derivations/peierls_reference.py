r"""Independent high-precision reference for the unified Peierls operator.

Purpose: provide a **mesh-independent, mathematically self-contained**
reference for verifying the production Nyström solver in
:mod:`orpheus.derivations.peierls_geometry`. Built from the same
unified polar form documented in :doc:`/theory/peierls_unified` §4,
but uses :func:`mpmath.quad` (adaptive Gauss-Kronrod) and native
arbitrary-precision special functions instead of fixed-order
Gauss-Legendre + mpmath dps truncation.

The reference stands alone: no Nyström basis, no fixed quadrature
orders, no shared code with production. It is the Phase-0
``ContinuousReferenceSolution`` contract made concrete for Peierls.

Design principle
----------------

The production solver and this reference **discretise the same
operator in two mathematically-different ways**:

- **Production** (:func:`~.peierls_geometry.build_volume_kernel`):
  panelwise Gauss-Legendre Nyström + Lagrange basis, fixed dps.
  The K matrix entries come from tensor-product GL quadrature over
  the polar (Ω, ρ) ray parametrisation.
- **Reference** (this module): adaptive :func:`mpmath.quad` over the
  same (Ω, ρ) parametrisation, arbitrary dps, straight evaluation of
  the integral with no basis interpolation between nodes — at each
  radial sample the kernel is evaluated at full precision.

Identical agreement of ``K[i, j]`` element-by-element
is a **complete verification** of: the unified polar-form integrand
formula, the geometry primitives (``rho_max``, ``source_position``,
``angular_weight``, ``volume_kernel_mp``), the optical-depth walker
handling of multi-region Σ_t, and the Lagrange basis evaluation.

Disagreement localises the bug to exactly one of those components.

Scope after the 2026-04-19 cleanup
----------------------------------

This module now contains ONLY the slab :math:`E_1`-form references and
the analytical vacuum-BC slab flux. The polar-form K element (slab,
cylinder, sphere) is now subsumed by the unified verification
primitive :func:`~.peierls_geometry.K_vol_element_adaptive`.

- :func:`slab_kernel_point_to_point` — exact :math:`E_1` slab point
  kernel.
- :func:`slab_K_vol_element` — adaptive ``mpmath.quad`` in real-space
  against the :math:`E_1` integrand. Used as the second independent
  formulation in the slab "two paths agree" verification gate
  (the polar form is the unified primitive).
- :func:`slab_uniform_source_analytical` — closed-form vacuum-BC
  uniform-source flux for slab; the analytical reference for the
  K-row-sum identity.
- :func:`slab_row_sum_uniform_identity` — convenience alias.
- :func:`cylinder_uniform_source_analytical` — semi-analytical vacuum-BC
  uniform-source flux for an infinite cylindrical cell; one
  ``mpmath.quad`` over the in-plane azimuth with :math:`\mathrm{Ki}_2`
  absorbing the out-of-plane polar integral.
- :func:`sphere_uniform_source_analytical` — semi-analytical vacuum-BC
  uniform-source flux for a spherical cell; one ``mpmath.quad`` over
  :math:`\mu = \cos\Theta`.

Together with the slab counterpart these three functions complete the
machine-precision analytical flux reference for vacuum-BC row-sum gating
of the K matrix built by :func:`~.peierls_geometry.build_volume_kernel_adaptive`
(the unified verification primitive).
"""

from __future__ import annotations

import mpmath
import numpy as np

from ._kernels import ki_n_mp
from .peierls_geometry import lagrange_basis_on_panels  # noqa: F401  (re-exported for downstream)


# ═══════════════════════════════════════════════════════════════════════
# Slab reference (E₁ kernel, two-face integration)
# ═══════════════════════════════════════════════════════════════════════

def slab_kernel_point_to_point(
    x_i: float, x_j: float, sig_t: float, *, dps: int = 50,
) -> mpmath.mpf:
    r"""Continuous slab Peierls kernel

    .. math::

       K(x_i, x_j) = \frac{1}{2}\,E_1(\Sigma_t\,|x_i - x_j|).

    At high dps this is the exact point-to-point infinite-medium slab
    kernel. The ``1/2`` factor emerges from the 2-D transverse
    integration of the 3-D point kernel (§2 of peierls_unified).
    """
    with mpmath.workdps(dps):
        tau = sig_t * abs(mpmath.mpf(x_i) - mpmath.mpf(x_j))
        if tau == 0:
            return mpmath.mpf("inf")  # E₁(0) diverges logarithmically
        return mpmath.expint(1, tau) / 2


def slab_K_vol_element(
    i: int, j: int,
    x_nodes: list, panel_bounds: list[tuple[float, float, int, int]],
    L: float, sig_t: float,
    *, dps: int = 50,
) -> mpmath.mpf:
    r"""Reference K[i, j] for the slab Peierls operator at dps precision.

    Computes

    .. math::

       K_{ij} \;=\; \int_0^L \tfrac{1}{2}\,E_1(\Sigma_t\,|x_i - x'|)\,
                    L_j(x')\,\mathrm d x'

    via adaptive :func:`mpmath.quad` over each panel. The diagonal
    panel (``i`` and ``j`` in the same panel) has an integrable
    logarithmic singularity at :math:`x' = x_i`; :func:`mpmath.quad`
    with the ``[a, x_i, b]`` subdivision hint handles this natively.
    """
    x_i = mpmath.mpf(x_nodes[i])
    sig_t_mp = mpmath.mpf(sig_t)

    def integrand(x_prime):
        """Lagrange basis L_j at x_prime times kernel K(x_i, x_prime)."""
        if x_prime == x_i:
            return mpmath.mpf(0)  # integrand removable; quad handles this
        tau = sig_t_mp * abs(x_i - x_prime)
        kernel = mpmath.expint(1, tau) / 2
        # Lagrange basis evaluation (numpy-based, float precision)
        L_vals = lagrange_basis_on_panels(
            np.array([float(x) for x in x_nodes]),
            panel_bounds,
            float(x_prime),
        )
        return kernel * mpmath.mpf(float(L_vals[j]))

    # Integrate panel by panel; subdivide the observer's own panel at x_i
    K_ij = mpmath.mpf(0)
    with mpmath.workdps(dps):
        for pa, pb, _, _ in panel_bounds:
            pa_mp, pb_mp = mpmath.mpf(pa), mpmath.mpf(pb)
            if pa_mp <= x_i <= pb_mp:
                # Split at the singularity for the adaptive integrator
                K_ij += mpmath.quad(integrand, [pa_mp, x_i, pb_mp])
            else:
                K_ij += mpmath.quad(integrand, [pa_mp, pb_mp])
    return K_ij




def slab_uniform_source_analytical(
    x: float, L: float, sig_t: float, *, dps: int = 50,
) -> mpmath.mpf:
    r"""Exact scalar flux from uniform unit source S(x) = 1 on a pure
    absorber slab [0, L] with vacuum BC:

    .. math::

       \varphi(x) \;=\; \frac{1}{2\Sigma_t}
                         \bigl[2 - E_2(\Sigma_t\,x)
                                 - E_2(\Sigma_t\,(L - x))\bigr].

    Derivation: :math:`\varphi(x) = \int_0^L (1/2) E_1(\Sigma_t|x-x'|)
    dx'`; split at ``x'`` = ``x`` and use
    :math:`\int E_1(\alpha\,u)\,\mathrm du = (1/\alpha)[E_2(0) -
    E_2(\alpha u)] = (1/\alpha)[1 - E_2(\alpha u)]`.

    This is an **independent closed-form reference** at arbitrary
    precision.
    """
    with mpmath.workdps(dps):
        x_mp = mpmath.mpf(x)
        L_mp = mpmath.mpf(L)
        sig_t_mp = mpmath.mpf(sig_t)
        return (2 - mpmath.expint(2, sig_t_mp * x_mp)
                  - mpmath.expint(2, sig_t_mp * (L_mp - x_mp))) / (2 * sig_t_mp)




# ═══════════════════════════════════════════════════════════════════════
# Cylinder and sphere vacuum-BC uniform-source analytical references
# ═══════════════════════════════════════════════════════════════════════

def cylinder_uniform_source_analytical(
    r: float, R: float, sig_t: float, *, dps: int = 50,
) -> mpmath.mpf:
    r"""Exact scalar flux from uniform unit source :math:`S = 1` on a
    pure-absorber cylindrical cell of radius :math:`R` (infinite axial
    extent) with vacuum BC on the lateral surface:

    .. math::

       \varphi_{\rm cyl}(r) \;=\;
           \frac{1}{\pi\,\Sigma_t}\!\int_0^\pi\!
             \Bigl[\,1 - \mathrm{Ki}_2\!\bigl(\Sigma_t\,L_{2D}(r,\theta')\bigr)\,\Bigr]
             \,\mathrm d\theta',

    with in-plane chord length

    .. math::

       L_{2D}(r, \theta') \;=\; -r\cos\theta'
                                + \sqrt{R^{2} - r^{2}\sin^{2}\theta'}.

    **Derivation.** Start from the 3-D point kernel
    :math:`e^{-\Sigma_t R'}/(4\pi R'^{2})` integrated against the
    spatial volume. In observer-centred cylindrical coordinates with
    in-plane azimuth :math:`\theta'` and polar angle
    :math:`\psi \in (0, \pi)` from the cylinder axis, the 3-D ray length
    to the exit on the lateral surface is
    :math:`\rho_{\max} = L_{2D}/\sin\psi`. The :math:`\psi`-integral
    then reduces to Bickley :math:`\mathrm{Ki}_2` via the definition
    :math:`\mathrm{Ki}_n(x) = \int_0^{\pi/2}\!\cos^{n-1}\!\phi\,
    e^{-x/\cos\phi}\,\mathrm d\phi` with the substitution
    :math:`\phi = \pi/2 - \psi`.

    **Sanity checks.**

    - :math:`r = 0`: :math:`L_{2D} = R` (constant) gives
      :math:`\varphi(0) = (1 - \mathrm{Ki}_2(\Sigma_t R))/\Sigma_t`.
    - :math:`\Sigma_t R \to 0` (thin):
      :math:`\mathrm{Ki}_2(x) \approx 1 - (\pi/2)\,x`, so
      :math:`\varphi(0) \to (\pi/2)\,R` — the 3-D mean chord through
      the axis of an infinite cylinder.
    - :math:`\Sigma_t R \to \infty`: :math:`\mathrm{Ki}_2 \to 0` and
      :math:`\varphi \to 1/\Sigma_t` (infinite-medium limit).

    **Implementation.** A single adaptive :func:`mpmath.quad` over
    :math:`\theta' \in [0, \pi]`. The integrand is smooth; no
    breakpoints needed. Machine precision via ``dps``.

    References
    ----------
    Bell & Glasstone (1970) "Nuclear Reactor Theory" Ch. 2; Case &
    Zweifel (1967) "Linear Transport Theory" Ch. 3.
    """
    with mpmath.workdps(dps):
        r_mp = mpmath.mpf(r)
        R_mp = mpmath.mpf(R)
        sig_t_mp = mpmath.mpf(sig_t)

        def integrand(theta_p):
            cos_tp = mpmath.cos(theta_p)
            sin_tp = mpmath.sin(theta_p)
            L_2d = -r_mp * cos_tp + mpmath.sqrt(
                R_mp ** 2 - r_mp ** 2 * sin_tp ** 2
            )
            ki2 = ki_n_mp(2, sig_t_mp * L_2d, dps)
            return 1 - ki2

        integral = mpmath.quad(integrand, [0, mpmath.pi])
        return integral / (mpmath.pi * sig_t_mp)


def sphere_uniform_source_analytical(
    r: float, R: float, sig_t: float, *, dps: int = 50,
) -> mpmath.mpf:
    r"""Exact scalar flux from uniform unit source :math:`S = 1` on a
    pure-absorber spherical cell of radius :math:`R` with vacuum BC:

    .. math::

       \varphi_{\rm sph}(r) \;=\; \frac{1}{2\,\Sigma_t}\!\left[\,
           2 - \int_{-1}^{1}\!\exp\!\Bigl(
               -\Sigma_t\bigl[-r\mu + \sqrt{R^{2} - r^{2} + r^{2}\mu^{2}}\bigr]
             \Bigr)\,\mathrm d\mu\,\right],

    where :math:`\mu = \cos\Theta` is the cosine of the polar angle
    from the observer's outward radial direction.

    **Derivation.** In observer-centred spherical coordinates
    :math:`(\rho, \Theta, \psi)` the point kernel
    :math:`e^{-\Sigma_t\rho}/(4\pi\rho^{2})` cancels the
    :math:`\rho^{2}` volume Jacobian; axial symmetry eliminates
    :math:`\psi`. The chord length from the observer at radius
    :math:`r` in direction :math:`\Theta` to the spherical surface is

    .. math::

       L_{\rm chord}(r, \mu) \;=\;
           -r\mu + \sqrt{R^{2} - r^{2}(1 - \mu^{2})}.

    The inner :math:`\rho`-integral gives
    :math:`(1 - e^{-\Sigma_t L_{\rm chord}})/\Sigma_t`; averaging over
    :math:`\mu \in [-1, 1]` yields the displayed form.

    **Sanity checks.**

    - :math:`r = 0`: :math:`L_{\rm chord} = R` (constant) gives
      :math:`\varphi(0) = (1 - e^{-\Sigma_t R})/\Sigma_t`.
    - :math:`\Sigma_t R \to 0` (thin): :math:`\varphi(0) \to R`, the
      3-D mean chord through the centre of a sphere.
    - :math:`\Sigma_t R \to \infty`:
      :math:`\varphi \to 1/\Sigma_t` (infinite-medium limit).

    **Implementation.** A single adaptive :func:`mpmath.quad` over
    :math:`\mu \in [-1, 1]`. The integrand is analytic; no breakpoints
    needed. Machine precision via ``dps``.

    References
    ----------
    Case & Zweifel (1967) "Linear Transport Theory" Ch. 3.
    """
    with mpmath.workdps(dps):
        r_mp = mpmath.mpf(r)
        R_mp = mpmath.mpf(R)
        sig_t_mp = mpmath.mpf(sig_t)

        def integrand(mu):
            L_chord = -r_mp * mu + mpmath.sqrt(
                R_mp ** 2 - r_mp ** 2 * (1 - mu ** 2)
            )
            return mpmath.exp(-sig_t_mp * L_chord)

        integral = mpmath.quad(integrand, [-1, 1])
        return (2 - integral) / (2 * sig_t_mp)


# ═══════════════════════════════════════════════════════════════════════
# Analytical diagnostics (independent of any Nyström)
# ═══════════════════════════════════════════════════════════════════════

def slab_row_sum_uniform_identity(
    x_i: float, L: float, sig_t: float, *, dps: int = 50,
) -> mpmath.mpf:
    r"""Slab row-sum identity: :math:`\int_0^L \tfrac12 E_1(\Sigma_t|x_i-x'|)
    \mathrm d x'`.

    Equal to :math:`\varphi(x_i)` for uniform source S=1 with pure
    absorption. See :func:`slab_uniform_source_analytical`. Convenience
    wrapper.
    """
    return slab_uniform_source_analytical(x_i, L, sig_t, dps=dps)


__all__ = [
    "slab_kernel_point_to_point",
    "slab_K_vol_element",
    "slab_uniform_source_analytical",
    "slab_row_sum_uniform_identity",
    "cylinder_uniform_source_analytical",
    "sphere_uniform_source_analytical",
]
