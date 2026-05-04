r"""Atkinson product-Nyström method for the slab Peierls integral
operator.

The 1G isotropic-scattering bare-critical slab Peierls integral
equation is

.. math::

   \phi(z) = \frac{c}{2} \int_{-a}^{a} E_1(\Sigma_t |z - z'|)\,
   \phi(z')\, dz', \qquad z \in [-a, a],

with vacuum boundaries (the ``-a``, ``+a`` limits make the
right-hand side a homogeneous Fredholm operator of the second
kind on :math:`L^\infty[-a, a]`). The eigenvalue is exactly 1 by
the criticality condition.

The kernel :math:`K(z, z') = (c/2)\, E_1(|z - z'|)` is
**logarithmically singular** at :math:`z = z'`. With the standard
small-argument expansion (Abramowitz–Stegun 5.1.11)

.. math::

   E_1(\tau) = -\gamma_E - \log\tau + \tau - \frac{\tau^2}{4}
   + \frac{\tau^3}{18} - \dots, \qquad \tau \to 0^+,

the kernel decomposes as

.. math::

   K(z, z') = \underbrace{-\frac{c}{2}\,\log|z - z'|}_{\text{singular}}
   + \underbrace{\frac{c}{2}\,\bigl[-\gamma_E + R(|z - z'|)\bigr]}
              _{\text{smooth},\;C^\infty},

where :math:`R(\tau) = E_1(\tau) + \log\tau` is the
:math:`C^\infty` Taylor remainder (analytic, with :math:`R(0) =
-\gamma_E`).

This module implements **Atkinson's product-Simpson rule** applied to
the singular log piece, and **standard Simpson** on the smooth
remainder, exactly as Atkinson 1972 [§III, Eqs. 16–17] specifies.
The singular piece is integrated **analytically** against the
piecewise-quadratic (Simpson) basis on each panel, eliminating the
log singularity from the discrete operator.

Why this exists
---------------

The earlier path
:func:`orpheus.derivations.continuous.fn_method.slab.flux_reconstruction.slab_scalar_flux_fn_projection`
discretises the kernel via plain Gauss–Legendre on the
:math:`(z, \mu)` tensor product. The :math:`\mu`-quadrature
:math:`\sum_k (w_k/\mu_k)\,\exp(-\tau/\mu_k)` is a discrete
truncation of :math:`E_1(\tau)`. Off-diagonal it converges to
machine precision; **at the diagonal** :math:`\tau = 0` it
saturates at a *finite* value (~ :math:`2 \log(n_\mu)`) instead of
diverging. The discrete kernel is therefore a *bandaged*
truncation of the true singular kernel — accurate
off-diagonal, qualitatively wrong at the diagonal.  The discrete
eigenmode picks this up as a systematic bias whose magnitude is
1–7 % depending on :math:`z/a` and :math:`c`.

Atkinson product-Nyström solves this at the operator level:
the discrete operator integrates the singular kernel
**analytically** against the Lagrange basis, never evaluating
:math:`E_1(0)`. Empirical convergence on the bare-critical slab
hits the F_N moment floor (≤ 1e-5 on flux ratios) at
:math:`n_\text{panels} \sim 128`–:math:`256`.

References
----------

* Atkinson, K.E. (1972) *Numer. Math.* **19**, 248–259.
* Atkinson, K.E. (1997) *The Numerical Solution of Integral
  Equations of the Second Kind*, §4.2 (esp. Eqs. 4.2.61–4.2.79,
  4.2.83 for the de Hoog–Weiss superconvergence rate).
* Abramowitz & Stegun (1964) §5.1.11 (the :math:`E_1` Taylor
  series used for the kernel decomposition).
"""
from __future__ import annotations

import math

import numpy as np
from scipy.special import exp1


_GAMMA_EULER = 0.5772156649015328606


def _F_k_log_primitives(t: float, s: float) -> tuple[float, float, float]:
    r"""Antiderivatives :math:`F_k(s; t) = \int s^k \log|t - s|\, ds`
    for :math:`k = 0, 1, 2`, evaluated at :math:`s`.

    .. math::

       F_0(s; t) &= (s - t)\,\log|s - t| - s, \\
       F_1(s; t) &= \tfrac{1}{2}(s^2 - t^2)\,\log|s - t|
                    - \tfrac{1}{4} s^2 - \tfrac{1}{2} t s, \\
       F_2(s; t) &= \tfrac{1}{3}(s^3 - t^3)\,\log|s - t|
                    - \tfrac{1}{9} s^3 - \tfrac{1}{6} t s^2
                    - \tfrac{1}{3} t^2 s.

    Each :math:`F_k` is finite at :math:`s = t` because
    :math:`\lim_{u \to 0} u\log|u| = 0`. The implementation uses
    that limit in the form of a guard:  when :math:`|s - t|` is
    smaller than the machine-epsilon scale, the
    :math:`\log|\cdot|`-prefactored terms are taken as zero by
    the limit, leaving the polynomial part.

    The three primitives are constructed by integration by parts:

    .. math::

       \int s^k \log|t - s|\, ds = \frac{s^{k+1}}{k+1} \log|t - s|
       - \frac{1}{k+1} \int \frac{s^{k+1}}{s - t}\, ds .

    Polynomial division of :math:`s^{k+1}/(s-t)` reduces the second
    integral to elementary terms plus :math:`t^{k+1} \log|s-t|`,
    which combines with the first term.

    Verified against SymPy: ``integrate(s**k * log(abs(t - s)), s)``.
    """
    d = s - t
    if abs(d) < 1e-300:
        # Use the limit lim_{u->0} u log|u| = 0; the
        # log-prefactored parts vanish, leaving polynomial pieces.
        log_abs_d = 0.0
        # And (s^k - t^k) -> 0, with multiplicand log -> 0, so all
        # log-prefactored terms drop.
    else:
        log_abs_d = math.log(abs(d))

    F0 = d * log_abs_d - s
    F1 = 0.5 * (s * s - t * t) * log_abs_d - 0.25 * s * s - 0.5 * t * s
    F2 = ((s ** 3 - t ** 3) / 3.0) * log_abs_d \
        - s ** 3 / 9.0 - t * s * s / 6.0 - t * t * s / 3.0
    return F0, F1, F2


def product_simpson_log_weights(t: float, a: float, b: float
                                ) -> tuple[float, float, float]:
    r"""Product-Simpson weights for :math:`\int_a^b \log|t - s|\,
    P(s)\, ds` against the unique quadratic :math:`P` interpolating
    its three node values at :math:`(a, m, b)` with :math:`m =
    (a+b)/2`.

    Returns :math:`(w_a, w_m, w_b)` such that

    .. math::

       \int_a^b \log|t - s|\, P(s)\, ds = w_a P(a) + w_m P(m) + w_b P(b)

    exactly for every quadratic :math:`P`. With

    .. math::

       L_a(s) = \frac{(s - m)(s - b)}{2 h^2}, \quad
       L_m(s) = -\frac{(s - a)(s - b)}{h^2}, \quad
       L_b(s) = \frac{(s - a)(s - m)}{2 h^2},

    where :math:`h = (b-a)/2`, the weights expand to linear
    combinations of :math:`F_k`:

    .. math::

       w_a &= \frac{1}{2 h^2}\bigl[F_2 - (m + b) F_1 + m b F_0\bigr], \\
       w_m &= -\frac{1}{h^2}\bigl[F_2 - (a + b) F_1 + a b F_0\bigr], \\
       w_b &= \frac{1}{2 h^2}\bigl[F_2 - (a + m) F_1 + a m F_0\bigr].

    Setup cost is :math:`O(1)` per panel — no quadrature is run at
    matrix-assembly time. Setup is closed-form via
    :func:`_F_k_log_primitives`.
    """
    m = 0.5 * (a + b)
    h = 0.5 * (b - a)
    F0_a, F1_a, F2_a = _F_k_log_primitives(t, a)
    F0_b, F1_b, F2_b = _F_k_log_primitives(t, b)
    F0 = F0_b - F0_a
    F1 = F1_b - F1_a
    F2 = F2_b - F2_a

    inv_2h2 = 1.0 / (2.0 * h * h)
    inv_h2 = 1.0 / (h * h)
    w_a = inv_2h2 * (F2 - (m + b) * F1 + m * b * F0)
    w_m = -inv_h2 * (F2 - (a + b) * F1 + a * b * F0)
    w_b = inv_2h2 * (F2 - (a + m) * F1 + a * m * F0)
    return w_a, w_m, w_b


def E1_smooth_remainder(tau: float) -> float:
    r"""Smooth :math:`C^\infty` remainder :math:`R(\tau) = E_1(\tau)
    + \log\tau`.

    From Abramowitz–Stegun 5.1.11:

    .. math::

       R(\tau) = -\gamma_E + \tau - \frac{\tau^2}{4}
                  + \frac{\tau^3}{18} - \dots,
                  \qquad R(0) = -\gamma_E.

    For :math:`\tau \le 10^{-15}` the leading two terms suffice
    (machine precision). For :math:`\tau > 10^{-15}` we use
    :func:`scipy.special.exp1` and add :math:`\log\tau`. The
    function is :math:`C^\infty` everywhere on :math:`[0, \infty)`.
    """
    if tau < 1e-15:
        # Leading-order Taylor: R(tau) ≈ -gamma_E + tau.
        return -_GAMMA_EULER + tau
    return float(exp1(tau)) + math.log(tau)


def _simpson_panel_smooth_weights(t: float, a: float, b: float,
                                  c_factor: float
                                  ) -> tuple[float, float, float]:
    r"""Standard Simpson weights for a single panel against the
    smooth remainder kernel :math:`(c/2)\,R(|t - s|)`.

    Returns :math:`(\tilde w_a, \tilde w_m, \tilde w_b)` such that

    .. math::

       \int_a^b \frac{c}{2}\,R(|t - s|)\, P(s)\, ds
       \approx \tilde w_a P(a) + \tilde w_m P(m) + \tilde w_b P(b),

    with the panel Simpson rule
    :math:`(h_p / 3) [f(a) + 4 f(m) + f(b)]` and
    :math:`h_p = (b - a)/2`. Here :math:`f(s) = (c/2)\,R(|t - s|)\,
    L_l(s)` for :math:`l \in \{a, m, b\}`; using the cardinal
    property :math:`L_a(a) = 1`, :math:`L_a(m) = L_a(b) = 0`, etc.,
    the panel-weights collapse to the kernel evaluated at the
    nodes:

    .. math::

       \tilde w_l = \frac{c}{2} \cdot (\text{Simpson coef})_l \cdot
       R(|t - s_l|).
    """
    m = 0.5 * (a + b)
    panel_h = (b - a) / 6.0  # so that Simpson is panel_h * (1, 4, 1)
    pref = c_factor / 2.0
    R_a = E1_smooth_remainder(abs(t - a))
    R_m = E1_smooth_remainder(abs(t - m))
    R_b = E1_smooth_remainder(abs(t - b))
    return (pref * panel_h * R_a,
            pref * panel_h * 4.0 * R_m,
            pref * panel_h * R_b)


def build_peierls_operator(c: float, a_half: float, n_panels: int
                            ) -> tuple[np.ndarray, np.ndarray]:
    r"""Build the Atkinson product-Simpson Peierls operator on a
    uniform mesh covering :math:`[-a_\text{half}, a_\text{half}]`.

    The mesh has :math:`n_\text{panels}` Simpson panels of width
    :math:`h = a_\text{half}/n_\text{panels}` (each panel covering
    width :math:`2h` from a Simpson-pair-of-intervals viewpoint, but
    here we use the convention "panel = quadratic interpolation
    interval" of width :math:`(b - a) = 2h`). The total node count
    is :math:`N = 2 n_\text{panels} + 1`.

    For each test node :math:`z_i` and each panel
    :math:`[s_{2j-2}, s_{2j}]`, the panel contributes three weighted
    rows to ``K_op``. The weighted rows are the sum of the
    product-Simpson weights for the singular log piece (negative
    sign because the kernel is :math:`-(c/2)\log`) and the
    standard Simpson weights for the smooth remainder.

    Parameters
    ----------
    c : float
        Mean number of secondaries per collision (:math:`> 1`).
    a_half : float
        Slab half-thickness in mean-free-paths.
    n_panels : int
        Number of Simpson panels covering :math:`[-a_\text{half},
        a_\text{half}]`. Must be :math:`\ge 1`. Total node count
        is :math:`2 n_\text{panels} + 1`.

    Returns
    -------
    K_op : (N, N) ndarray
        Discrete Peierls operator. Applying ``K_op @ phi`` gives
        :math:`(c/2) \int E_1(|z - z'|)\,\phi(z')\,dz'` evaluated
        at the same nodes.
    z_nodes : (N,) ndarray
        Equispaced node positions covering :math:`[-a_\text{half},
        a_\text{half}]`.

    Notes
    -----
    Setup cost is :math:`O(N^2)` (every test node looks at every
    panel). This is the same complexity as the plain-GL operator
    used in the legacy ``slab_scalar_flux_fn_projection``.  No
    runtime quadrature is performed: the singular weights are
    closed-form via :func:`_F_k_log_primitives` and the smooth
    weights are three function evaluations of :func:`exp1` per
    panel per test node.
    """
    if n_panels < 1:
        raise ValueError(f"n_panels must be >= 1, got {n_panels}")
    if c <= 0.0:
        raise ValueError(f"c must be positive, got {c}")
    if a_half <= 0.0:
        raise ValueError(f"a_half must be positive, got {a_half}")

    n_nodes = 2 * n_panels + 1
    z_nodes = np.linspace(-a_half, a_half, n_nodes)

    K_op = np.zeros((n_nodes, n_nodes))
    half_c = c / 2.0

    for j_panel in range(n_panels):
        i_lo, i_md, i_hi = 2 * j_panel, 2 * j_panel + 1, 2 * j_panel + 2
        s_lo, s_hi = z_nodes[i_lo], z_nodes[i_hi]

        for i_test, z_i in enumerate(z_nodes):
            # Singular log piece: closed-form product weights.
            w_a_log, w_m_log, w_b_log = product_simpson_log_weights(
                z_i, s_lo, s_hi
            )
            # Smooth piece: standard Simpson with kernel folded in.
            w_a_sm, w_m_sm, w_b_sm = _simpson_panel_smooth_weights(
                z_i, s_lo, s_hi, c
            )
            # K = (c/2) * [-log(|z-z'|) + R(|z-z'|)], so combine.
            K_op[i_test, i_lo] += -half_c * w_a_log + w_a_sm
            K_op[i_test, i_md] += -half_c * w_m_log + w_m_sm
            K_op[i_test, i_hi] += -half_c * w_b_log + w_b_sm

    return K_op, z_nodes


def power_iterate_dominant_eigenmode(K_op: np.ndarray,
                                     symmetric_about_index: bool = True,
                                     tol: float = 1e-13,
                                     max_iter: int = 2000
                                     ) -> tuple[float, np.ndarray]:
    r"""Extract the dominant (eigenvalue, eigenvector) pair of
    ``K_op`` by power iteration.

    For the bare-critical slab the dominant eigenvalue is exactly
    1 by the criticality condition, and the eigenvector is the
    symmetric scalar flux :math:`\phi(z)` (even about :math:`z = 0`).

    Parameters
    ----------
    K_op : (N, N) ndarray
        Discrete Peierls operator from :func:`build_peierls_operator`.
    symmetric_about_index : bool, default True
        If True, project iterates onto the even subspace at every
        step. This kills any odd-mode contamination from numerical
        noise and accelerates convergence to the dominant
        symmetric eigenmode.
    tol : float, default 1e-13
        Convergence criterion on the Rayleigh-quotient successive
        difference.
    max_iter : int, default 2000
        Maximum power-iteration steps. The slab Peierls operator
        has a moderate spectral gap so 100–500 iterations
        typically suffice.

    Returns
    -------
    eig : float
        Converged dominant eigenvalue. For a true bare-critical
        slab problem this should be 1.0 to within F_N-coefficient
        accuracy.
    phi : (N,) ndarray
        Converged eigenvector, normalized to ``np.max|phi| = 1``.
    """
    n = K_op.shape[0]
    phi_curr = np.cos(np.pi * np.arange(n) / (n - 1) - 0.5 * np.pi)
    if symmetric_about_index:
        phi_curr = 0.5 * (phi_curr + phi_curr[::-1])
    norm = np.max(np.abs(phi_curr))
    if norm > 0:
        phi_curr = phi_curr / norm

    eig_prev = 0.0
    for _ in range(max_iter):
        phi_next = K_op @ phi_curr
        # Rayleigh quotient.
        num = float(np.dot(phi_next, phi_curr))
        den = float(np.dot(phi_curr, phi_curr))
        if abs(den) < 1e-300:
            break
        eig_next = num / den
        if symmetric_about_index:
            phi_next = 0.5 * (phi_next + phi_next[::-1])
        norm = np.max(np.abs(phi_next))
        if norm < 1e-300:
            break
        phi_next = phi_next / norm
        if abs(eig_next - eig_prev) < tol:
            return float(eig_next), phi_next
        phi_curr = phi_next
        eig_prev = eig_next

    return float(eig_next), phi_curr
