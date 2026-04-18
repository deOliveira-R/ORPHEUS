r"""Peierls integral equation reference for cylindrical CP verification.

Solves the 1-D radial Peierls (integral transport) equation on a bare
or annular cylinder via Nyström quadrature at mpmath precision,
producing :class:`ContinuousReferenceSolution` objects whose operator
form is ``"integral-peierls"``.

The cylindrical Peierls equation for a radially-symmetric emission
density :math:`q(r)` on :math:`0 \le r \le R` is obtained by
integrating the 3-D point kernel over the infinite z-axis, yielding
the 2-D transverse Green's function
:math:`G_{2D}(|\mathbf{r}-\mathbf{r}'|) = \mathrm{Ki}_1(\Sigma_t\,
|\mathbf{r}-\mathbf{r}'|) / (2\pi\,|\mathbf{r}-\mathbf{r}'|)`:

.. math::

   \varphi(\mathbf{r})
     \;=\; \frac{1}{2\pi}\!\iint_{\rm disc}
       \frac{\mathrm{Ki}_1\!\bigl(\tau(\mathbf{r},\mathbf{r}')\bigr)}
            {|\mathbf{r}-\mathbf{r}'|}\,q(\mathbf{r}')\,\mathrm{d}^{2}r'
     + \varphi_{\rm bc}(\mathbf{r})

For the Nyström discretisation this module expresses the 2-D
integral in **polar coordinates centred at the observer**:

.. math::

   \varphi(r)
     \;=\; \frac{1}{\pi}\!
       \int_{0}^{\pi}\!\mathrm{d}\beta\!
       \int_{0}^{\rho_{\max}(r,\beta)}\!\!
         \mathrm{Ki}_1\!\bigl(\tau(r, \rho, \beta)\bigr)\,
         q\bigl(r'(r, \rho, \beta)\bigr)\,\mathrm{d}\rho
     + \varphi_{\rm bc}(r)

where :math:`\rho` is the distance from the observer along the ray
at angle :math:`\beta` (measured from the outward radial direction
at the observer), :math:`r'(r,\rho,\beta) = \sqrt{r^{2} + 2 r\rho
\cos\beta + \rho^{2}}` is the source radius, :math:`\rho_{\max}
= -r\cos\beta + \sqrt{r^{2}\cos^{2}\beta + R^{2}-r^{2}}` is the
distance to the cylinder boundary, and the prefactor :math:`1/\pi`
absorbs the :math:`1/(2\pi)` of the 2-D kernel plus a factor of 2
from the :math:`y\to-y` (β-reflection) symmetry that folds
:math:`\beta\in[0,2\pi]` to :math:`[0,\pi]`.

Compared to the equivalent chord :math:`(y, r')` form used in
Sanchez-McCormick 1982 §IV.A, the polar form has **no singular
Jacobian** — the 2-D area element :math:`\rho\,\mathrm{d}\rho\,
\mathrm{d}\beta` cancels the :math:`1/\rho` in the Green's function,
leaving a smooth integrand :math:`\mathrm{Ki}_1(\tau)\,q(r')` that
is handled cleanly by ordinary Gauss–Legendre quadrature in both
:math:`\beta` and :math:`\rho`. The chord :math:`(y, r')`
parametrisation picks up the Jacobian
:math:`1/\sqrt{(r^{2}-y^{2})(r'^{2}-y^{2})}` from the two-branch sum
:math:`|\mathrm{d}\alpha_{+}/\mathrm{d}y| + |\mathrm{d}\alpha_{-}/
\mathrm{d}y| = 2/\sqrt{\min(r,r')^{2}-y^{2}}` and carries an extra
singularity at :math:`y=\min(r,r')`. The polar form avoids this.

Key differences from :mod:`peierls_slab`:

- **Kernel**: :math:`\mathrm{Ki}_1` (Bickley–Naylor order 1), not
  :math:`E_1`. The slab's :math:`E_1` comes from integrating the 1-D
  point kernel over polar angle; the cylinder's :math:`\mathrm{Ki}_1`
  comes from integrating the 3-D point kernel over the infinite
  axial direction.
- **Singularity**: none in the polar form (:math:`\rho\,\mathrm{d}\rho\,
  \mathrm{d}\beta` cancels the :math:`1/\rho` of
  :math:`\mathrm{Ki}_1/|\mathbf{r}-\mathbf{r}'|`). Contrast the slab's
  log-singular kernel :math:`E_1(\tau)\sim -\ln\tau`.
- **Prefactor**: :math:`1/\pi`, not :math:`1/2`. Derived above.
- **Source interpolation**: because :math:`r'(\rho,\beta,r_i)` is
  generally not a quadrature node, the Nyström unknown is connected
  to :math:`q(r'_{ikm})` via Lagrange-basis interpolation on the
  radial grid — hence the kernel matrix picks up a
  :math:`L_j(r'_{ikm})` factor.
- **White BC closure**: rank-:math:`N_\rho` dense Schur block
  (continuous lateral surface; C5), not the slab's rank-2 E₂ outer
  product (two discrete faces).

The :math:`\tau^{\pm}` chord walker (:func:`optical_depths_pm`, C2)
remains part of this module: it is the primitive used to evaluate
:math:`\tau(r,\rho,\beta)` for **multi-region** problems (C4),
where the optical depth along the ray from :math:`r_i` in direction
:math:`\beta` to a source at distance :math:`\rho` decomposes into a
same-side / through-centre branch depending on whether the ray
crosses the chord midpoint.

This module is the Phase-4.2 deliverable of the verification campaign.
It is the independent reference against which the flat-source CP
cylinder solver (:mod:`orpheus.cp.solver` on ``cyl1D`` meshes and
:mod:`orpheus.derivations.cp_cylinder`) is verified.

.. note::

   This is the C2 scaffold. The Nyström kernel builder, eigenvalue
   power iteration, and ``continuous_cases()`` registration land in
   subsequent commits (C3–C7 of the Phase-4.2 plan).
"""

from __future__ import annotations

from dataclasses import dataclass

import mpmath
import numpy as np

from ._kernels import ki_n_mp
from ._reference import (
    ContinuousReferenceSolution,
    ProblemSpec,
    Provenance,
)
from ._xs_library import LAYOUTS, get_mixture, get_xs


# ═══════════════════════════════════════════════════════════════════════
# Composite Gauss–Legendre y-quadrature with breakpoints at annular radii
# ═══════════════════════════════════════════════════════════════════════

def _gl_nodes_weights(n: int, dps: int) -> tuple[list, list]:
    """*n*-point Gauss–Legendre on [-1, 1] at *dps* decimal digits."""
    with mpmath.workdps(dps):
        nm, wm = mpmath.gauss_quadrature(n, "legendre")
        return [nm[i] for i in range(n)], [wm[i] for i in range(n)]


def _map_gl_to(nodes, weights, a, b):
    """Map GL nodes/weights from [-1, 1] to [a, b] at mpmath precision."""
    h = (b - a) / 2
    m = (a + b) / 2
    return [m + h * t for t in nodes], [h * w for w in weights]


def composite_gl_y(
    radii: np.ndarray,
    n_panels_per_region: int,
    p_order: int,
    dps: int = 30,
) -> tuple[np.ndarray, np.ndarray, list[tuple[float, float, int, int]]]:
    r"""Composite Gauss–Legendre quadrature for the y-integration.

    For the cylinder Peierls equation, the outer integration variable
    is the chord impact parameter :math:`y \in [0, R]`. The integrand
    :math:`y \mapsto \mathrm{Ki}_1(\tau^\pm(r, r', y))` has **corners**
    at each annular radius :math:`r_k` (where a chord at :math:`y`
    transitions from crossing annulus :math:`k` to only grazing it),
    so the quadrature uses a composite GL rule with breakpoints at
    each :math:`r_k`.

    Each of the :math:`N` annular segments :math:`[r_{k-1}, r_k]`
    carries ``n_panels_per_region`` panels of order ``p_order``.

    Parameters
    ----------
    radii : np.ndarray, shape (N,)
        Outer radii of the :math:`N` concentric annuli,
        :math:`0 < r_1 < \dots < r_N = R`.
    n_panels_per_region : int
        Number of GL panels per annular segment.
    p_order : int
        GL order per panel.
    dps : int
        mpmath working precision.

    Returns
    -------
    y_pts : np.ndarray, shape (n_y,)
        Quadrature nodes in :math:`[0, R]`, in ascending order, at
        double precision (float).
    y_wts : np.ndarray, shape (n_y,)
        Quadrature weights.
    panel_bounds : list of (pa, pb, i_start, i_end)
        Breakdown of the composite rule: one tuple per panel, with
        panel endpoints and the slice of ``y_pts`` / ``y_wts`` that
        lives on it.
    """
    radii = np.asarray(radii, dtype=float)
    gl_ref, gl_wt = _gl_nodes_weights(p_order, dps)

    breakpoints = [mpmath.mpf(0)] + [mpmath.mpf(float(r)) for r in radii]
    y_all: list = []
    w_all: list = []
    panel_bounds: list[tuple[float, float, int, int]] = []

    with mpmath.workdps(dps):
        for seg in range(len(breakpoints) - 1):
            a_seg = breakpoints[seg]
            b_seg = breakpoints[seg + 1]
            pw = (b_seg - a_seg) / n_panels_per_region
            for pidx in range(n_panels_per_region):
                pa = a_seg + pidx * pw
                pb = pa + pw
                xp, wp = _map_gl_to(gl_ref, gl_wt, pa, pb)
                i0 = len(y_all)
                y_all.extend(xp)
                w_all.extend(wp)
                panel_bounds.append((float(pa), float(pb), i0, len(y_all)))

    y_pts = np.array([float(y) for y in y_all])
    y_wts = np.array([float(w) for w in w_all])
    return y_pts, y_wts, panel_bounds


# ═══════════════════════════════════════════════════════════════════════
# τ⁺ / τ⁻ optical-path walker
# ═══════════════════════════════════════════════════════════════════════

def optical_depths_pm(
    r: float,
    r_prime: float,
    y_pts: np.ndarray,
    radii: np.ndarray,
    sig_t: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Optical paths :math:`\tau^{+}(r,r',y)` and :math:`\tau^{-}(r,r',y)`.

    For a chord of impact parameter :math:`y`, a point at radius
    :math:`\rho \ge y` sits at signed chord position
    :math:`\pm s_\rho` with :math:`s_\rho = \sqrt{\rho^{2}-y^{2}}`.
    The two Peierls-kernel branches integrate
    :math:`\Sigma_t` along:

    - **Same-side (τ⁺)**: from :math:`s_r` to :math:`s_{r'}` on the
      positive branch (both points on the same side of the
      perpendicular foot).
    - **Through-centre (τ⁻)**: from :math:`s_r` to :math:`-s_{r'}`,
      crossing the chord midpoint.

    The integrand :math:`\Sigma_t(\rho(s))` is piecewise constant with
    jumps at the chord crossings of each annular radius
    :math:`r_k`, i.e. at :math:`|s| = s_{r_k}`.

    For :math:`y > r` or :math:`y > r'`, the point is *not* on this
    chord. Following Sanchez's convention we return :math:`\tau = 0`
    for such inaccessible configurations; the Nyström kernel
    multiplies by :math:`\mathrm{Ki}_1(0) = \pi/2` on those, but the
    outer integral in :math:`y` weights them as zero-measure (the
    geometric :math:`y`-range is :math:`[0, \min(r, R)]`).

    Parameters
    ----------
    r, r_prime : float
        Source and target radii (:math:`\ge 0`, :math:`\le R`).
    y_pts : np.ndarray, shape (n_y,)
        Chord impact parameters.
    radii : np.ndarray, shape (N,)
        Outer radii of the :math:`N` concentric annuli.
    sig_t : np.ndarray, shape (N,)
        Total macroscopic cross-section per annulus, for a single
        energy group.

    Returns
    -------
    tau_plus : np.ndarray, shape (n_y,)
    tau_minus : np.ndarray, shape (n_y,)
    """
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    y_pts = np.asarray(y_pts, dtype=float)

    r = float(r)
    r_prime = float(r_prime)

    # Signed-chord positions of the endpoints. For y > r or y > r',
    # the corresponding s_r / s_r' is zero (geometrically: the radius
    # is inside the chord-inaccessible core).
    s_r = np.sqrt(np.maximum(r ** 2 - y_pts ** 2, 0.0))
    s_rp = np.sqrt(np.maximum(r_prime ** 2 - y_pts ** 2, 0.0))

    # Positive-branch annular breakpoints:
    #   s_breaks[0] = 0  (the chord midpoint / perpendicular foot)
    #   s_breaks[k] = sqrt(max(r_k² - y², 0))  for k = 1..N
    # Annulus k (0-indexed in the arrays, 1-indexed in docs) occupies
    # the positive chord interval [s_breaks[k], s_breaks[k+1]].
    N = len(radii)
    s_breaks = np.zeros((N + 1, len(y_pts)))
    for k in range(N):
        s_breaks[k + 1] = np.sqrt(np.maximum(radii[k] ** 2 - y_pts ** 2, 0.0))

    def _optical_on_positive(s_lo: np.ndarray, s_hi: np.ndarray) -> np.ndarray:
        r"""Integral of :math:`\Sigma_t(\rho(s))` over :math:`s \in [s_{lo}, s_{hi}]`
        with both endpoints on the positive branch (:math:`s_{lo} \le s_{hi}`)."""
        tau = np.zeros_like(s_lo)
        for k in range(N):
            lo = np.maximum(s_lo, s_breaks[k])
            hi = np.minimum(s_hi, s_breaks[k + 1])
            overlap = np.maximum(hi - lo, 0.0)
            tau += sig_t[k] * overlap
        return tau

    # Same-side: integrate between s_r and s_rp on the positive branch.
    s_lo = np.minimum(s_r, s_rp)
    s_hi = np.maximum(s_r, s_rp)
    tau_plus = _optical_on_positive(s_lo, s_hi)

    # Through-centre: from +s_r to -s_rp. By symmetry, this equals
    # (integral from 0 to s_r on positive branch) + (integral from 0
    # to s_rp on positive branch).
    zero = np.zeros_like(y_pts)
    tau_minus = _optical_on_positive(zero, s_r) + _optical_on_positive(zero, s_rp)

    return tau_plus, tau_minus


# ═══════════════════════════════════════════════════════════════════════
# Volume kernel assembly (polar (β, ρ) Nyström with Lagrange interpolation)
# ═══════════════════════════════════════════════════════════════════════

def _gl_float(n: int, a: float, b: float, dps: int = 30) -> tuple[np.ndarray, np.ndarray]:
    """*n*-point GL on ``[a, b]``, returned as float arrays."""
    ref_nodes, ref_wts = _gl_nodes_weights(n, dps)
    h = (b - a) / 2
    m = (a + b) / 2
    nodes = np.array([float(m + h * t) for t in ref_nodes])
    wts = np.array([float(h * w) for w in ref_wts])
    return nodes, wts


def _lagrange_basis_on_panels(
    r_nodes: np.ndarray,
    panel_bounds: list[tuple[float, float, int, int]],
    r_eval: float,
) -> np.ndarray:
    r"""Evaluate the Lagrange basis :math:`L_j(r_{\rm eval})` at an
    arbitrary point, using the panel structure of the r-grid.

    On each panel :math:`[p_a, p_b]`, the Lagrange basis polynomials are
    built over that panel's nodes only; on other panels the basis is
    zero (piecewise polynomial representation). This matches how
    :class:`PeierlsCylinderSolution` would interpolate a discrete
    nodal vector back to a continuous function.

    Parameters
    ----------
    r_nodes : np.ndarray, shape (N,)
        All r-nodes across all panels.
    panel_bounds : list of (pa, pb, i_start, i_end)
        Panel layout from :func:`composite_gl_y` / :func:`composite_gl_r`.
    r_eval : float
        Point at which to evaluate the basis.

    Returns
    -------
    np.ndarray, shape (N,)
        :math:`L_j(r_{\rm eval})` for :math:`j=0,\dots,N-1`. Nonzero
        only on the panel containing :math:`r_{\rm eval}`.
    """
    N = len(r_nodes)
    L = np.zeros(N)

    # Find the panel containing r_eval (boundary right-biased)
    panel_idx = None
    for k, (pa, pb, i_start, i_end) in enumerate(panel_bounds):
        if pa <= r_eval <= pb:
            panel_idx = k
            break
    if panel_idx is None:
        # r_eval out of [0, R]: clamp to nearest panel endpoint
        if r_eval < panel_bounds[0][0]:
            panel_idx = 0
        else:
            panel_idx = len(panel_bounds) - 1

    pa, pb, i_start, i_end = panel_bounds[panel_idx]
    local_nodes = r_nodes[i_start:i_end]
    p = i_end - i_start
    for a in range(p):
        num, den = 1.0, 1.0
        for b in range(p):
            if b == a:
                continue
            num *= (r_eval - local_nodes[b])
            den *= (local_nodes[a] - local_nodes[b])
        L[i_start + a] = num / den
    return L


def _rho_max(r_obs: float, cos_beta: float, R: float) -> float:
    r"""Distance along ray at angle β from observer at :math:`r_{\rm obs}`
    to the cylinder boundary :math:`|\mathbf{r}| = R`.

    Positive root of :math:`(r_{\rm obs} + \rho\cos\beta)^{2}
    + (\rho\sin\beta)^{2} = R^{2}`.
    """
    disc = r_obs * r_obs * cos_beta * cos_beta + R * R - r_obs * r_obs
    return -r_obs * cos_beta + np.sqrt(max(disc, 0.0))


def _optical_depth_along_ray(
    r_obs: float,
    cos_beta: float,
    sin_beta: float,
    rho: float,
    radii: np.ndarray,
    sig_t: np.ndarray,
) -> float:
    r"""Optical depth from :math:`r_{\rm obs}` along direction
    :math:`(\cos\beta, \sin\beta)` over a distance :math:`\rho`.

    Homogeneous (1-region) short-circuit: :math:`\tau = \Sigma_t\,\rho`.
    Multi-region: the ray crosses annular boundaries
    :math:`r = r_k` at roots of :math:`r_{\rm obs}^{2} + 2r_{\rm obs}
    s\cos\beta + s^{2} = r_k^{2}` (a quadratic in :math:`s`); the
    optical depth is the sum of :math:`\Sigma_{t,k}\cdot\Delta s`
    over each annular segment of the ray. This routine handles
    both cases.
    """
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    N = len(radii)

    # Fast path: homogeneous
    if N == 1:
        return float(sig_t[0]) * rho

    # Multi-region: find all boundary crossings s_k where |r(s)| = r_k
    # r(s)² = r_obs² + 2 r_obs s cos β + s²
    # Crossing r_k: s² + 2 r_obs cos β · s + (r_obs² - r_k²) = 0
    crossings = [0.0]
    for r_k in radii[:-1]:  # outer boundary handled by rho bounds
        disc = r_obs * r_obs * cos_beta * cos_beta - (r_obs * r_obs - r_k * r_k)
        if disc < 0:
            continue  # ray does not intersect this annular boundary
        sqrt_disc = np.sqrt(disc)
        s_a = -r_obs * cos_beta - sqrt_disc
        s_b = -r_obs * cos_beta + sqrt_disc
        for s in (s_a, s_b):
            if 0.0 < s < rho:
                crossings.append(s)
    crossings.append(rho)
    crossings.sort()

    # Between consecutive crossings, the ray is inside a single annulus.
    # Determine annulus index by mid-segment radius.
    tau = 0.0
    for i_seg in range(len(crossings) - 1):
        s_lo, s_hi = crossings[i_seg], crossings[i_seg + 1]
        s_mid = 0.5 * (s_lo + s_hi)
        r_mid_sq = r_obs * r_obs + 2.0 * r_obs * s_mid * cos_beta + s_mid * s_mid
        r_mid = np.sqrt(max(r_mid_sq, 0.0))
        # Annulus k contains r_mid iff r_{k-1} ≤ r_mid < r_k (r_0 = 0).
        # Default: outermost annulus (handles r_mid just outside radii[-1]
        # from floating-point noise at the cylinder boundary).
        k = N - 1
        for kk in range(N):
            if r_mid < radii[kk]:
                k = kk
                break
        tau += sig_t[k] * (s_hi - s_lo)
    return tau


def build_volume_kernel(
    r_nodes: np.ndarray,
    panel_bounds: list[tuple[float, float, int, int]],
    radii: np.ndarray,
    sig_t: np.ndarray,
    n_beta: int,
    n_rho: int,
    dps: int = 30,
) -> np.ndarray:
    r"""Assemble the **volume** Nyström kernel matrix for a single group.

    The cylindrical Peierls equation in polar coordinates centred at
    the observer reads

    .. math::

       \Sigma_t(r_i)\,\varphi(r_i)
         \;=\; \frac{\Sigma_t(r_i)}{\pi}
           \int_{0}^{\pi}\!\mathrm{d}\beta\!
           \int_{0}^{\rho_{\max}(r_i,\beta)}\!
             \mathrm{Ki}_1\!\bigl(\tau(r_i,\rho,\beta)\bigr)\,
             q\!\bigl(r'(r_i,\rho,\beta)\bigr)\,\mathrm{d}\rho
         + S_{\rm bc}(r_i).

    With Lagrange interpolation
    :math:`q(r'_{ikm}) = \sum_j L_j(r'_{ikm})\,q_j`, this discretises to

    .. math::

       \Sigma_t(r_i)\,\varphi_i
         \;=\; \sum_j K_{ij}\,q_j + S_{\rm bc}(r_i),
       \qquad
       K_{ij} = \frac{\Sigma_t(r_i)}{\pi}\sum_{k,m} w_{\beta,k}\,
                  w_{\rho,m}(r_i,\beta_k)\,
                  \mathrm{Ki}_1(\tau_{ikm})\,L_j(r'_{ikm}).

    Row-sum identity (for the **1-group, homogeneous, pure scatterer**
    check): at :math:`\varphi\equiv1,\,q=\Sigma_t`, the infinite-medium
    identity :math:`\sum_j K_{ij} \to \Sigma_t` holds exactly, and the
    finite-cylinder deficit equals :math:`\Sigma_t P_{\rm esc}(r_i)`
    where :math:`P_{\rm esc}` is the uncollided-escape probability
    from :math:`r_i`. The test gate uses :math:`R = 10` mean free
    paths so that :math:`P_{\rm esc}(r_i \lesssim R/2) \ll 10^{-3}`.

    Parameters
    ----------
    r_nodes : np.ndarray, shape (N,)
        Radial quadrature nodes (composite GL on :math:`[0, R]`).
    panel_bounds : list of (pa, pb, i_start, i_end)
        Panel layout of the r-grid, from :func:`composite_gl_r`. Used
        for Lagrange-basis interpolation on the local panel.
    radii : np.ndarray, shape (N_reg,)
        Outer radii of concentric annuli.
    sig_t : np.ndarray, shape (N_reg,)
        Total macroscopic cross-section per annulus.
    n_beta : int
        GL order for the :math:`\beta`-integral on :math:`[0, \pi]`.
    n_rho : int
        GL order for the :math:`\rho`-integral on
        :math:`[0, \rho_{\max}(r_i, \beta_k)]`.
    dps : int
        mpmath working precision for :math:`\mathrm{Ki}_1`.

    Returns
    -------
    K : np.ndarray, shape (N, N)
        Nyström kernel matrix.
    """
    r_nodes = np.asarray(r_nodes, dtype=float)
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    N = len(r_nodes)
    R = float(radii[-1])

    # β-quadrature on [0, π]
    beta_pts, beta_wts = _gl_float(n_beta, 0.0, np.pi, dps)
    cos_betas = np.cos(beta_pts)
    sin_betas = np.sin(beta_pts)

    # Reference ρ-quadrature nodes on [-1, 1] — we map per (i, k)
    ref_rho_nodes, ref_rho_wts = _gl_nodes_weights(n_rho, dps)
    ref_rho_nodes = np.array([float(x) for x in ref_rho_nodes])
    ref_rho_wts = np.array([float(w) for w in ref_rho_wts])

    K = np.zeros((N, N))
    inv_pi = 1.0 / np.pi

    for i in range(N):
        r_i = r_nodes[i]
        # Σ_t(r_i): locate the annulus containing r_i (default to
        # outermost for any r_i that lands exactly on the outer radius
        # due to floating-point noise).
        k_obs = len(radii) - 1
        for kk in range(len(radii)):
            if r_i < radii[kk]:
                k_obs = kk
                break
        sig_t_i = sig_t[k_obs]

        for k in range(n_beta):
            cb = cos_betas[k]
            sb = sin_betas[k]
            rho_max = _rho_max(r_i, cb, R)
            if rho_max <= 0.0:
                continue

            # Map reference ρ-nodes [-1, 1] → [0, ρ_max]
            h = 0.5 * rho_max
            rho_pts = h * ref_rho_nodes + h
            rho_wts = h * ref_rho_wts

            for m in range(n_rho):
                rho = rho_pts[m]
                r_prime = np.sqrt(r_i * r_i + 2.0 * r_i * rho * cb + rho * rho)
                tau = _optical_depth_along_ray(
                    r_i, cb, sb, rho, radii, sig_t,
                )
                ki1 = float(ki_n_mp(1, float(tau), dps))
                L_vals = _lagrange_basis_on_panels(
                    r_nodes, panel_bounds, float(r_prime),
                )
                weight = inv_pi * sig_t_i * beta_wts[k] * rho_wts[m] * ki1
                K[i, :] += weight * L_vals

    return K


def composite_gl_r(
    radii: np.ndarray,
    n_panels_per_region: int,
    p_order: int,
    dps: int = 30,
) -> tuple[np.ndarray, np.ndarray, list[tuple[float, float, int, int]]]:
    """Composite GL quadrature on :math:`[0, R]` with panel breakpoints
    at each annular radius. Same structure as :func:`composite_gl_y`;
    exposed separately because the radial grid and the y-grid may
    have different resolutions in practice.
    """
    return composite_gl_y(radii, n_panels_per_region, p_order, dps=dps)


# ═══════════════════════════════════════════════════════════════════════
# 1G eigenvalue solver (vacuum BC)
# ═══════════════════════════════════════════════════════════════════════

def _which_annulus(r: float, radii: np.ndarray) -> int:
    """Index of the annulus containing ``r`` (outer-biased for r on a
    panel boundary; matches :func:`build_volume_kernel`'s own lookup)."""
    k = len(radii) - 1
    for kk, r_k in enumerate(radii):
        if r < r_k:
            return kk
    return k


# ═══════════════════════════════════════════════════════════════════════
# White-BC surface currents (rank-1 closure)
# ═══════════════════════════════════════════════════════════════════════

def compute_P_esc(
    r_nodes: np.ndarray,
    radii: np.ndarray,
    sig_t: np.ndarray,
    n_beta: int = 32,
    dps: int = 25,
) -> np.ndarray:
    r"""Uncollided escape probability :math:`P_{\rm esc}(r_i)`.

    For a neutron emitted isotropically at radius :math:`r_i` in an
    infinite-z cylindrical cell, the probability of reaching the
    lateral surface without a collision is

    .. math::

        P_{\rm esc}(r_i)
          \;=\; \frac{1}{\pi}\!\int_{0}^{\pi}\!
             \mathrm{Ki}_2\!\bigl(\tau(r_i, \rho_{\max}(r_i, \beta), \beta)\bigr)\,
             \mathrm{d}\beta

    obtained by integrating the isotropic angular kernel (3-D) through
    :math:`\int_{0}^{\pi/2}\sin\theta\,\exp(-\tau/\sin\theta)\,
    \mathrm{d}\theta = \mathrm{Ki}_2(\tau)` and folding over the
    azimuthal symmetry. For a multi-region problem, :math:`\tau` is
    the line integral along the ray through all annuli.
    """
    r_nodes = np.asarray(r_nodes, dtype=float)
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    R = float(radii[-1])

    beta_pts, beta_wts = _gl_float(n_beta, 0.0, np.pi, dps)
    cos_betas = np.cos(beta_pts)
    sin_betas = np.sin(beta_pts)
    inv_pi = 1.0 / np.pi

    N = len(r_nodes)
    P_esc = np.zeros(N)
    for i in range(N):
        r_i = r_nodes[i]
        total = 0.0
        for k in range(n_beta):
            cb = cos_betas[k]
            sb = sin_betas[k]
            rho_max = _rho_max(r_i, cb, R)
            if rho_max <= 0.0:
                continue
            tau = _optical_depth_along_ray(
                r_i, cb, sb, rho_max, radii, sig_t,
            )
            total += beta_wts[k] * float(ki_n_mp(2, float(tau), dps))
        P_esc[i] = inv_pi * total
    return P_esc


def compute_G_bc(
    r_nodes: np.ndarray,
    radii: np.ndarray,
    sig_t: np.ndarray,
    n_phi: int = 32,
    dps: int = 25,
) -> np.ndarray:
    r"""Uncollided surface-to-volume Green's function integrated over
    a uniform isotropic lateral-surface source.

    For a unit uniform partial current :math:`J^{-}` re-entering the
    cylinder with isotropic inward distribution, the contribution to
    the scalar flux at :math:`r_i` is :math:`J^{-}\,G_{\rm bc}(r_i)`
    where

    .. math::

        G_{\rm bc}(r_i)
          \;=\; \frac{2 R}{\pi}\!\int_{0}^{\pi}\!
             \frac{\mathrm{Ki}_1\!\bigl(\tau_{\rm surf}(r_i, \phi)\bigr)}
                  {d(r_i, R, \phi)}\,\mathrm{d}\phi

    with :math:`d(r_i, R, \phi) = \sqrt{r_i^{2} + R^{2} - 2 r_i R
    \cos\phi}` the chord from surface-point
    :math:`S = (R\cos\phi, R\sin\phi)` to observer :math:`(r_i, 0)`
    and :math:`\tau_{\rm surf}` the optical depth along that chord
    through all annuli crossed.

    The prefactor :math:`2R/\pi` reflects the isotropic-inward
    hemisphere normalisation: :math:`J^{-} = \pi\,\psi_{\rm in}` for
    isotropic :math:`\psi_{\rm in}`, and the total inward emission
    rate per unit surface area is :math:`2\,J^{-}`, which drives the
    :math:`2D` line-source Green's function
    :math:`\mathrm{Ki}_1(\Sigma_t d)/(2\pi d)` integrated over the
    surface.
    """
    r_nodes = np.asarray(r_nodes, dtype=float)
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    R = float(radii[-1])

    phi_pts, phi_wts = _gl_float(n_phi, 0.0, np.pi, dps)
    cos_phis = np.cos(phi_pts)
    inv_pi = 1.0 / np.pi

    N = len(r_nodes)
    G_bc = np.zeros(N)
    for i in range(N):
        r_i = r_nodes[i]
        total = 0.0
        for k in range(n_phi):
            cf = cos_phis[k]
            d_sq = r_i * r_i + R * R - 2.0 * r_i * R * cf
            d = np.sqrt(max(d_sq, 0.0))
            if d <= 0.0:
                continue
            # Ray from surface-point S to r_i — compute optical depth
            # by walking annular crossings between r_i and the surface.
            # The ray direction (from S toward r_i) has cos_beta_ray
            # and sin_beta_ray in the observer's polar frame; but the
            # optical-depth-along-ray walker treats the ray as going
            # from r_i outward along (cos_β, sin_β) for distance ρ.
            # Going outward by ρ_max(r_i, β) lands on the surface; we
            # want the path to the *specific* surface point S(φ). For
            # a homogeneous 1-region cylinder this is just Σ_t · d.
            # For multi-region, we need to walk along the chord from
            # r_i through the boundaries out to S.
            if len(radii) == 1:
                tau = sig_t[0] * d
            else:
                # Ray direction from r_i toward S: (x,y) = S - r_i =
                # (R·cφ − r_i, R·sφ) normalised. cos_β = (R·cφ − r_i)/d,
                # sin_β = R·sφ/d.
                R_sf = R * np.sin(phi_pts[k])
                cb = (R * cf - r_i) / d
                sb = R_sf / d
                tau = _optical_depth_along_ray(
                    r_i, cb, sb, d, radii, sig_t,
                )
            ki1 = float(ki_n_mp(1, float(tau), dps))
            total += phi_wts[k] * ki1 / d
        G_bc[i] = 2.0 * inv_pi * R * total
    return G_bc


def build_white_bc_correction(
    r_nodes: np.ndarray,
    r_wts: np.ndarray,
    radii: np.ndarray,
    sig_t: np.ndarray,
    *,
    n_beta: int = 32,
    n_phi: int = 32,
    dps: int = 25,
) -> np.ndarray:
    r"""Rank-1 white-BC correction matrix for :func:`build_volume_kernel`.

    Returns :math:`K_{\rm bc}[i, j] = u[i]\,v[j]` where

    .. math::

        u[i] = \Sigma_t(r_i)\,G_{\rm bc}(r_i) / R, \qquad
        v[j] = r_j\,w_j\,P_{\rm esc}(r_j).

    Add :math:`K_{\rm bc}` elementwise to :math:`K_{\rm vol}` (from
    :func:`build_volume_kernel`) to obtain the white-BC operator.
    For a radially-symmetric 1-D cylinder the outgoing and incoming
    partial currents are uniform over the surface.

    .. warning::

       **Approximation level.** Rank-1 is the correct rank for
       the *partial current balance* alone: :math:`J^{-} = J^{+}`
       collapses to a scalar because both are uniform over the
       cylinder surface. But the angular distribution of the
       outgoing and incoming currents differs — white BC assumes
       isotropic incoming (Mark), while the actual outgoing
       distribution is anisotropic. At the pointwise-Nyström level
       (as opposed to flat-source CP, where region averaging
       cancels the anisotropy), this approximation shows up as a
       spread in the uniform-flux row-sum identity
       :math:`(K_{\rm vol}+K_{\rm bc})\cdot\mathbf 1 \approx
       \Sigma_t` that grows with inverse cell size:

       =====  ======================
       R/MFP  max \|K_tot·1 − Σt\|
       =====  ======================
       0.5    0.32
       1.0    0.16
       2.0    0.20
       5.0    0.12
       10     < 0.04
       =====  ======================

       Consequently the white-BC :math:`k_{\rm eff}` agrees with
       :math:`k_\infty` (the Wigner-Seitz exact result) only
       asymptotically:

       =====  ==========  ==========
       R/MFP  k(white)    err vs k∞
       =====  ==========  ==========
       1.0    1.19        21 %
       2.0    1.40        7 %
       5.0    1.48        2 %
       10     1.49        1 %
       =====  ==========  ==========

       The vacuum-BC driver remains bit-exact against the Sanchez
       tie-point; use it for rigorous prefactor checks. Tests that
       compare the Peierls cylinder reference to CP (white BC)
       should use :math:`R \ge 5` MFP to keep the rank-1 closure
       error under 3 %.

       This limitation is analogous to issue #100 for the sphere,
       where the pointwise rank-1 white-BC closure has the same
       structural deficit. A rigorous fix requires either (a) a
       higher-rank angular decomposition of the surface currents,
       or (b) a fully-iterative surface-source treatment that
       captures the anisotropic outgoing distribution. Both are
       deferred follow-up work.
    """
    r_nodes = np.asarray(r_nodes, dtype=float)
    r_wts = np.asarray(r_wts, dtype=float)
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    R = float(radii[-1])
    N = len(r_nodes)

    # Per-node Σ_t
    sig_t_n = np.empty(N)
    for i, ri in enumerate(r_nodes):
        sig_t_n[i] = sig_t[_which_annulus(ri, radii)]

    P_esc = compute_P_esc(r_nodes, radii, sig_t, n_beta=n_beta, dps=dps)
    G_bc = compute_G_bc(r_nodes, radii, sig_t, n_phi=n_phi, dps=dps)

    u = sig_t_n * G_bc / R
    v = r_nodes * r_wts * P_esc
    return np.outer(u, v)


def solve_peierls_cylinder_1g(
    radii: np.ndarray,
    sig_t: np.ndarray,
    sig_s: np.ndarray,
    nu_sig_f: np.ndarray,
    *,
    boundary: str = "vacuum",
    n_panels_per_region: int = 2,
    p_order: int = 5,
    n_beta: int = 24,
    n_rho: int = 24,
    n_phi: int = 24,
    dps: int = 25,
    max_iter: int = 300,
    tol: float = 1e-10,
) -> PeierlsCylinderSolution:
    r"""Solve the 1-group cylindrical Peierls k-eigenvalue problem.

    Parameters
    ----------
    boundary : {"vacuum", "white"}
        BC at the outer cylinder radius. Vacuum = no re-entry; white
        = isotropic re-entry with :math:`J^{-} = J^{+}`. For a
        homogeneous 1-region cylinder the white-BC :math:`k_{\rm eff}`
        equals :math:`k_\infty` exactly (Wigner-Seitz property).
    n_phi : int
        GL order for the surface-φ integral in :func:`compute_G_bc`
        (used only for white BC).

    Other parameters as in :func:`build_volume_kernel`.
    """
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    sig_s = np.asarray(sig_s, dtype=float)
    nu_sig_f = np.asarray(nu_sig_f, dtype=float)

    r_nodes, r_wts, panels = composite_gl_r(
        radii, n_panels_per_region, p_order, dps=dps,
    )
    K = build_volume_kernel(
        r_nodes, panels, radii, sig_t, n_beta=n_beta, n_rho=n_rho, dps=dps,
    )

    if boundary == "white":
        K_bc = build_white_bc_correction(
            r_nodes, r_wts, radii, sig_t,
            n_beta=n_beta, n_phi=n_phi, dps=dps,
        )
        K = K + K_bc
    elif boundary != "vacuum":
        raise ValueError(f"boundary must be 'vacuum' or 'white', got {boundary!r}")

    N = len(r_nodes)
    sig_t_n = np.empty(N)
    sig_s_n = np.empty(N)
    nu_sig_f_n = np.empty(N)
    for i, ri in enumerate(r_nodes):
        ki = _which_annulus(ri, radii)
        sig_t_n[i] = sig_t[ki]
        sig_s_n[i] = sig_s[ki]
        nu_sig_f_n[i] = nu_sig_f[ki]

    A = np.diag(sig_t_n) - K * sig_s_n[np.newaxis, :]
    B = K * nu_sig_f_n[np.newaxis, :]

    phi = np.ones(N)
    k_val = 1.0
    B_phi = B @ phi
    prod_old = np.abs(B_phi).sum()

    for it in range(max_iter):
        q = B_phi / k_val
        phi_new = np.linalg.solve(A, q)
        B_phi_new = B @ phi_new
        prod_new = np.abs(B_phi_new).sum()
        k_new = k_val * prod_new / prod_old if prod_old > 0 else k_val
        nrm = np.abs(phi_new).sum()
        if nrm > 0:
            phi_new = phi_new / nrm
        B_phi_norm = B @ phi_new
        prod_norm = np.abs(B_phi_norm).sum()
        converged = abs(k_new - k_val) < tol and it > 5
        phi, k_val = phi_new, k_new
        B_phi, prod_old = B_phi_norm, prod_norm
        if converged:
            break

    return PeierlsCylinderSolution(
        r_nodes=r_nodes,
        phi_values=phi[:, np.newaxis],
        k_eff=float(k_val),
        cell_radius=float(radii[-1]),
        n_groups=1,
        n_quad_r=N,
        n_quad_y=n_beta * n_rho,
        precision_digits=dps,
    )


def solve_peierls_cylinder_1g_vacuum(
    radii: np.ndarray,
    sig_t: np.ndarray,
    sig_s: np.ndarray,
    nu_sig_f: np.ndarray,
    *,
    n_panels_per_region: int = 2,
    p_order: int = 5,
    n_beta: int = 24,
    n_rho: int = 24,
    dps: int = 25,
    max_iter: int = 300,
    tol: float = 1e-10,
) -> PeierlsCylinderSolution:
    r"""Solve the 1-group cylindrical Peierls k-eigenvalue problem with
    **vacuum** boundary conditions.

    The equation

    .. math::

       \Sigma_{t,i}\,\varphi_i
         \;=\; \sum_j K_{ij}\bigl(\Sigma_{s,j}\,\varphi_j
                                 + \tfrac{1}{k}\nu\Sigma_{f,j}\,\varphi_j\bigr)

    is recast as the generalised eigenvalue problem
    :math:`\tilde A\,\varphi = (1/k)\,\tilde B\,\varphi` with
    :math:`\tilde A_{ij} = \delta_{ij}\,\Sigma_{t,i} - K_{ij}\,\Sigma_{s,j}`
    and :math:`\tilde B_{ij} = K_{ij}\,\nu\Sigma_{f,j}`. Dominant
    eigenvalue via fission-source power iteration (mirrors
    :func:`orpheus.derivations.peierls_slab.solve_peierls_eigenvalue`).

    Parameters
    ----------
    radii : np.ndarray, shape (N_reg,)
        Outer radii of annular regions.
    sig_t, sig_s, nu_sig_f : np.ndarray, shape (N_reg,)
        Single-group cross-sections per region.
    n_panels_per_region, p_order : int
        Composite-GL radial mesh controls.
    n_beta, n_rho : int
        Quadrature orders for the :math:`(\beta, \rho)` polar
        integration in :func:`build_volume_kernel`.
    dps : int
        mpmath working precision for Ki₁.

    Returns
    -------
    PeierlsCylinderSolution

    Notes
    -----
    The critical-cylinder verification target is Sanchez-McCormick
    1982 Table IV: :math:`\Sigma_t = 1\,{\rm cm}^{-1}`,
    :math:`(\Sigma_s + \nu\Sigma_f)/\Sigma_t = 1.5`,
    :math:`R = 1.9798\,{\rm cm}` gives :math:`k_{\rm eff} = 1.0`.
    """
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    sig_s = np.asarray(sig_s, dtype=float)
    nu_sig_f = np.asarray(nu_sig_f, dtype=float)

    r_nodes, r_wts, panels = composite_gl_r(
        radii, n_panels_per_region, p_order, dps=dps,
    )
    K = build_volume_kernel(
        r_nodes, panels, radii, sig_t, n_beta=n_beta, n_rho=n_rho, dps=dps,
    )

    # Per-node cross-sections
    N = len(r_nodes)
    sig_t_n = np.empty(N)
    sig_s_n = np.empty(N)
    nu_sig_f_n = np.empty(N)
    for i, ri in enumerate(r_nodes):
        ki = _which_annulus(ri, radii)
        sig_t_n[i] = sig_t[ki]
        sig_s_n[i] = sig_s[ki]
        nu_sig_f_n[i] = nu_sig_f[ki]

    # A = diag(Σ_t) - K · diag(Σ_s)
    A = np.diag(sig_t_n) - K * sig_s_n[np.newaxis, :]
    # B = K · diag(νΣ_f)
    B = K * nu_sig_f_n[np.newaxis, :]

    # Power iteration: A φ = (1/k) B φ
    phi = np.ones(N)
    k_val = 1.0
    B_phi = B @ phi
    prod_old = np.abs(B_phi).sum()

    for it in range(max_iter):
        q = B_phi / k_val
        phi_new = np.linalg.solve(A, q)
        B_phi_new = B @ phi_new
        prod_new = np.abs(B_phi_new).sum()
        k_new = k_val * prod_new / prod_old if prod_old > 0 else k_val

        nrm = np.abs(phi_new).sum()
        if nrm > 0:
            phi_new = phi_new / nrm
        B_phi_norm = B @ phi_new
        prod_norm = np.abs(B_phi_norm).sum()

        converged = abs(k_new - k_val) < tol and it > 5
        phi, k_val = phi_new, k_new
        B_phi, prod_old = B_phi_norm, prod_norm
        if converged:
            break

    return PeierlsCylinderSolution(
        r_nodes=r_nodes,
        phi_values=phi[:, np.newaxis],
        k_eff=float(k_val),
        cell_radius=float(radii[-1]),
        n_groups=1,
        n_quad_r=N,
        n_quad_y=n_beta * n_rho,
        precision_digits=dps,
    )


# ═══════════════════════════════════════════════════════════════════════
# Solution container
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class PeierlsCylinderSolution:
    """Result of a Peierls Nyström solve on a 1-D radial cylinder.

    Fields
    ------
    r_nodes : (N,) float
        Radial quadrature node positions.
    phi_values : (N, ng) float
        Flux at each radial node and group.
    k_eff : float or None
        Eigenvalue (None for fixed-source).
    cell_radius, n_groups, n_quad_r, n_quad_y, precision_digits
        Metadata.
    panel_bounds : list of (pa, pb, i_start, i_end) or None
        Panel layout for piecewise-Lagrange interpolation. When
        populated (default for eigenvalue solves), :meth:`phi`
        evaluates the flux at arbitrary :math:`r` via the same
        Lagrange basis used to assemble the Nyström kernel.
    """

    r_nodes: np.ndarray
    phi_values: np.ndarray
    k_eff: float | None
    cell_radius: float
    n_groups: int
    n_quad_r: int
    n_quad_y: int
    precision_digits: int
    panel_bounds: list[tuple[float, float, int, int]] | None = None

    def phi(self, r: np.ndarray, g: int = 0) -> np.ndarray:
        """Evaluate flux at arbitrary radii via piecewise Lagrange basis.

        The basis mirrors :func:`_lagrange_basis_on_panels` — on each
        panel the basis is the Lagrange polynomial of the panel's
        nodes; the piecewise representation matches what the Nyström
        operator itself used. For radii outside :math:`[0, R]` the
        result is clamped to the boundary panel.
        """
        r = np.asarray(r, dtype=float).ravel()
        out = np.empty_like(r)

        if self.panel_bounds is None:
            # Degenerate fallback: linear interp
            return np.interp(r, self.r_nodes, self.phi_values[:, g])

        for idx, r_eval in enumerate(r):
            L = _lagrange_basis_on_panels(
                self.r_nodes, self.panel_bounds, float(r_eval),
            )
            out[idx] = float(np.dot(L, self.phi_values[:, g]))
        return out


# ═══════════════════════════════════════════════════════════════════════
# ContinuousReferenceSolution builders
# ═══════════════════════════════════════════════════════════════════════

_MAT_IDS_CYL = {1: [2]}  # single-region 1-group case


def _build_peierls_cylinder_case(
    ng_key: str,
    n_regions: int,
    n_panels_per_region: int = 3,
    p_order: int = 5,
    n_beta: int = 20,
    n_rho: int = 20,
    n_phi: int = 20,
    precision_digits: int = 20,
) -> ContinuousReferenceSolution:
    """Build a Peierls-cylinder reference matching a cp_cyl1D case.

    Currently supports 1G 1-region only; multi-group and multi-region
    extensions are follow-up work. White BC is used to match the CP
    cylinder solver's Wigner-Seitz cell convention.
    """
    if ng_key != "1g" or n_regions != 1:
        raise NotImplementedError(
            f"peierls_cylinder continuous reference currently supports "
            f"1G 1-region only; got ng_key={ng_key!r}, n_regions={n_regions}"
        )

    from .cp_cylinder import _RADII

    layout = LAYOUTS[n_regions]
    ng = int(ng_key[0])
    radii = np.array(_RADII[n_regions], dtype=float)

    xs_list = [get_xs(region, ng_key) for region in layout]
    sig_t = np.array([xs["sig_t"][0] for xs in xs_list])
    sig_s = np.array([xs["sig_s"][0, 0] for xs in xs_list])
    nu_sig_f = np.array([(xs["nu"] * xs["sig_f"])[0] for xs in xs_list])

    r_nodes, r_wts, panels = composite_gl_r(
        radii, n_panels_per_region, p_order, dps=precision_digits,
    )
    K_vol = build_volume_kernel(
        r_nodes, panels, radii, sig_t,
        n_beta=n_beta, n_rho=n_rho, dps=precision_digits,
    )
    K_bc = build_white_bc_correction(
        r_nodes, r_wts, radii, sig_t,
        n_beta=n_beta, n_phi=n_phi, dps=precision_digits,
    )
    K = K_vol + K_bc

    N = len(r_nodes)
    sig_t_n = np.full(N, sig_t[0])
    sig_s_n = np.full(N, sig_s[0])
    nu_sig_f_n = np.full(N, nu_sig_f[0])
    A = np.diag(sig_t_n) - K * sig_s_n[np.newaxis, :]
    B = K * nu_sig_f_n[np.newaxis, :]

    # Power iteration
    phi = np.ones(N)
    k_val = 1.0
    B_phi = B @ phi
    prod_old = np.abs(B_phi).sum()
    for it in range(500):
        q = B_phi / k_val
        phi_new = np.linalg.solve(A, q)
        B_phi_new = B @ phi_new
        prod_new = np.abs(B_phi_new).sum()
        k_new = k_val * prod_new / prod_old if prod_old > 0 else k_val
        nrm = np.abs(phi_new).sum()
        if nrm > 0:
            phi_new = phi_new / nrm
        B_phi_norm = B @ phi_new
        prod_norm = np.abs(B_phi_norm).sum()
        converged = abs(k_new - k_val) < 1e-10 and it > 5
        phi, k_val = phi_new, k_new
        B_phi, prod_old = B_phi_norm, prod_norm
        if converged:
            break

    # Normalise so that ∫ φ dV = 1 over the cell
    # For 2D radial: dV = 2π r dr, so integral = 2π · Σ_j r_j w_j φ_j
    integral = 2.0 * np.pi * np.dot(r_nodes * r_wts, phi)
    if abs(integral) > 1e-30:
        phi = phi / integral

    sol = PeierlsCylinderSolution(
        r_nodes=r_nodes,
        phi_values=phi[:, np.newaxis],
        k_eff=float(k_val),
        cell_radius=float(radii[-1]),
        n_groups=ng,
        n_quad_r=N,
        n_quad_y=n_beta * n_rho,
        precision_digits=precision_digits,
        panel_bounds=panels,
    )

    def phi_fn(x: np.ndarray, g: int = 0) -> np.ndarray:
        return sol.phi(x, g)

    mat_ids = _MAT_IDS_CYL[n_regions]
    materials = {
        mat_ids[i]: get_mixture(region, ng_key)
        for i, region in enumerate(layout)
    }

    return ContinuousReferenceSolution(
        name=f"peierls_cyl1D_{ng}eg_{n_regions}rg",
        problem=ProblemSpec(
            materials=materials,
            geometry_type="cylinder-1d",
            geometry_params={
                "radius": float(radii[-1]),
                "radii": radii.tolist(),
                "mat_ids": mat_ids,
            },
            boundary_conditions={"outer": "white"},
            is_eigenvalue=True,
            n_groups=ng,
        ),
        operator_form="integral-peierls",
        phi=phi_fn,
        k_eff=sol.k_eff,
        provenance=Provenance(
            citation=(
                "Sanchez & McCormick 1982 (NSE 80) §IV.A; "
                "Hebert 2020 Ch. 3 §3.5"
            ),
            derivation_notes=(
                f"Polar (β, ρ) Nyström at {n_panels_per_region} panels × "
                f"{p_order} GL points on [0, R], with n_β = {n_beta}, "
                f"n_ρ = {n_rho}. Ki₁ kernel, no singular Jacobian. "
                f"White BC via rank-1 Schur closure (radial symmetry "
                f"collapses the N_β block to a single scalar J⁻ = J⁺)."
            ),
            sympy_expression=None,
            precision_digits=precision_digits,
        ),
        equation_labels=(
            "peierls-cylinder-equation",
            "peierls-cylinder-polar",
            "peierls-cylinder-ray-optical-depth",
        ),
        vv_level="L1",
        description=(
            f"{ng}G {n_regions}-region cylindrical Peierls "
            f"(Ki₁ polar Nyström, rank-1 white BC)"
        ),
        tolerance="O(h²)",
    )


def continuous_cases() -> list[ContinuousReferenceSolution]:
    """Peierls cylinder continuous references for the registry.

    .. note::

       Currently returns an **empty list**. The single 1G 1-region
       case :func:`_build_peierls_cylinder_case` is defined and
       usable for direct-call tests, but not wired into the shared
       ``reference_values`` registry because the rank-1 white-BC
       closure is only accurate at :math:`R\\ge 5` MFP and the
       existing ``cp_cyl1D_1eg_1rg`` case uses :math:`R = 1` MFP.
       Registering an approximate reference would poison
       downstream tests that expect registry entries to carry
       full-precision numerical provenance.

       Once a higher-rank white-BC closure is implemented (see the
       caveat block in :func:`build_white_bc_correction`), this
       function will be extended to return the usual 1eg_1rg /
       2eg_2rg grid of cases.
    """
    return []
