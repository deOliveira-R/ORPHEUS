r"""Plan-2 follow-on A2 — Pomraning-Siewert 1982 reference solver for
homogeneous sphere with vacuum BC and isotropic scattering.

This is the **structurally-independent reference** for Variant α
verification at :math:`\alpha = 0`. PS-1982 derives the same
:math:`[E_1(|r-x|) - E_1(r+x)]` kernel as Sanchez 1986 but by a
genuinely different mathematical path (radial-µ integration + half-space
addition vs cosh-even-extension). Same final kernel, different routes —
this rules out reference contamination.

Architecture
------------

PS-1982 Eq. (21) for the homogeneous (:math:`Q = 0`) vacuum-BC case
reduces to:

.. math::

   r\,\phi(r) \;=\; \frac{c}{2}\!\int_0^R x\,\phi(x)\,
                    \big[E_1(|r-x|) - E_1(r+x)\big]\,\mathrm d x

with :math:`c = (\Sigma_s + \nu\Sigma_f/k)/\Sigma_t` the single-scattering
albedo (in optical units, :math:`r, x \in (0, R\,\Sigma_t]`).

Define the integral operator :math:`L[\phi](r) = (1/r)\int_0^R x\,\phi(x)\,
K(r,x)\,\mathrm d x` with :math:`K(r,x) = E_1(|r-x|) - E_1(r+x)`. The
eigenvalue equation :math:`\phi = (c/2)\,L[\phi]` admits a dominant
eigenvalue :math:`\lambda_{\max}(R)` (for fixed optical sphere radius
:math:`R` in MFP units) such that the critical albedo is
:math:`c_* = 2/\lambda_{\max}`. The k-eigenvalue extension: the
steady-state with fission requires :math:`c_{\rm eff} = c_*`, giving

.. math::

   k_{\rm eff} = \frac{\nu\Sigma_f}{c_*\,\Sigma_t - \Sigma_s}.

For closed sphere (:math:`R \to \infty`), :math:`\lambda_{\max} \to 2`
(no leakage) so :math:`c_* \to 1` and :math:`k_{\rm eff} \to
\nu\Sigma_f/\Sigma_a = k_\infty`. Thinner spheres have :math:`\lambda
_{\max} < 2` (leakage reduces the effective albedo), so :math:`c_* > 1`
— more multiplication is needed to reach criticality.

Numerical method
----------------

The kernel :math:`K(r,x)` has a logarithmic singularity at :math:`x = r`
(:math:`E_1(t) \sim -\ln(t) - \gamma` as :math:`t \to 0^+`), so
straight Nyström breaks down. Instead we use product-integration:

1. Maintain :math:`\phi` on a Gauss-Legendre grid :math:`\{r_i\}` in
   optical units.
2. At each iteration, build a cubic-spline interpolant
   :math:`\phi_{\rm interp}(x)`.
3. For each :math:`r_i`, compute :math:`L[\phi](r_i) = (1/r_i)\int_0^R
   x\,\phi_{\rm interp}(x)\,K(r_i, x)\,\mathrm d x` via
   :func:`scipy.integrate.quad` with the singular point :math:`x = r_i`
   as an explicit breakpoint (``points=[r_i]``). This routes the
   log-singular integrand through QUADPACK's QAGS adaptive scheme
   (Wynn epsilon table extrapolation), which is designed for exactly
   this case.
4. Power-iterate to find :math:`\lambda_{\max}` and convert to
   :math:`k_{\rm eff}`.

This is slower than direct matrix Nyström but **correct**, and
the prototype is meant for verification (used as L1 reference for
Variant α), not production speed.

References
----------

- Pomraning, G. C. & Siewert, C. E. (1982). "On the integral form of
  the equation of transfer for a homogeneous sphere". *J. Quant. Spec.
  Rad. Transfer* **28**, 503-506. DOI `10.1016/0022-4073(82)90016-4`.
- Memo: `.claude/agent-memory/literature-researcher/ps1982_and_garcia_extraction.md`
- Atkinson, K. E. (1997). *Numerical Solution of Integral Equations of
  the Second Kind*, Cambridge. §7 product integration for log-singular
  kernels.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from scipy.special import exp1


@dataclass(frozen=True)
class PS1982Result:
    """Result of PS-1982 power iteration for vacuum sphere k_eff."""

    k_eff: float
    lambda_max: float          # dominant eigenvalue of L
    c_critical: float          # critical albedo = 2/lambda_max
    phi: np.ndarray            # converged scalar flux on the GL grid (optical units)
    r_nodes_optical: np.ndarray  # r·Σ_t (MFP units)
    iterations: int
    converged: bool


def _apply_L_operator(
    phi: np.ndarray,
    r_nodes: np.ndarray,
    R_opt: float,
    *,
    quad_limit: int = 200,
    quad_epsabs: float = 1e-12,
    quad_epsrel: float = 1e-10,
) -> np.ndarray:
    r"""Apply the PS-1982 vacuum kernel operator L to a discrete
    :math:`\phi` on the GL grid.

    :math:`L[\phi](r_i) = (1/r_i) \int_0^R x\,\phi(x)\,[E_1(|r_i-x|) -
    E_1(r_i+x)]\,\mathrm d x`. The integral is computed via QUADPACK
    QAGS with the singular point :math:`x = r_i` as a breakpoint.
    """
    phi_interp = CubicSpline(r_nodes, phi, extrapolate=True)
    out = np.zeros_like(phi)
    for i, r in enumerate(r_nodes):
        def integrand(x: float, r=r) -> float:
            return float(
                x
                * phi_interp(x)
                * (exp1(abs(r - x)) - exp1(r + x))
            )
        # Breakpoint at x = r handles the log singularity.
        # Split into [0, r] and [r, R] explicitly via `points`.
        val, _ = quad(
            integrand, 0.0, R_opt,
            points=[r],
            limit=quad_limit,
            epsabs=quad_epsabs,
            epsrel=quad_epsrel,
        )
        out[i] = val / r
    return out


def solve_ps1982_vacuum_sphere(
    R: float,
    sigma_t: float,
    sigma_s: float,
    nu_sigma_f: float,
    *,
    n_quad: int = 32,
    max_iter: int = 200,
    tol: float = 1e-9,
) -> PS1982Result:
    r"""Power iteration on the PS-1982 vacuum-sphere kernel; extract
    k_eff via the homogeneous-eigenvalue trick.

    Iterates :math:`\phi^{(n+1)} = L[\phi^{(n)}]` to convergence, then
    computes the dominant eigenvalue :math:`\lambda_{\max}` via
    Rayleigh quotient and converts to :math:`k_{\rm eff}` via
    :math:`k_{\rm eff} = \nu\Sigma_f / (c_*\Sigma_t - \Sigma_s)` with
    :math:`c_* = 2/\lambda_{\max}`.

    Parameters
    ----------
    R, sigma_t, sigma_s, nu_sigma_f : float
        Sphere radius (cm) and 1G cross sections (cm⁻¹).
    n_quad : int
        Gauss-Legendre order for the radial grid.
    max_iter, tol : int, float
        Power-iteration parameters.

    Returns
    -------
    :class:`PS1982Result`
    """
    R_opt = R * sigma_t

    pts, _ = np.polynomial.legendre.leggauss(n_quad)
    r_nodes = R_opt * 0.5 * (pts + 1.0)

    # Initial guess: uniform constant.
    phi = np.ones(n_quad)
    phi /= np.linalg.norm(phi)
    lam = 1.0

    converged = False
    iterations = 0
    for it in range(max_iter):
        iterations = it + 1
        L_phi = _apply_L_operator(phi, r_nodes, R_opt)
        # Rayleigh quotient via inner product (no metric — phi already L2-normed).
        lam_new = float(np.dot(phi, L_phi))
        norm = np.linalg.norm(L_phi)
        if norm < 1e-30:
            raise RuntimeError("L_phi vanished — operator pathology?")
        phi_new = L_phi / norm
        rel = abs(lam_new - lam) / max(abs(lam_new), 1e-30)
        phi = phi_new
        lam = lam_new
        if rel < tol:
            converged = True
            break

    c_critical = 2.0 / lam
    denom = c_critical * sigma_t - sigma_s
    if denom <= 0:
        raise ValueError(
            f"PS-1982: c_*·Σ_t − Σ_s = {denom} ≤ 0 (c_* = {c_critical}, "
            f"σ_t = {sigma_t}, σ_s = {sigma_s}). For this XS the bare "
            f"scattering ratio σ_s/σ_t = {sigma_s/sigma_t:.4f} already "
            f"exceeds the critical c_* — sphere is supercritical from "
            f"scattering alone (try larger σ_a or smaller R)."
        )
    k_eff = nu_sigma_f / denom

    return PS1982Result(
        k_eff=float(k_eff),
        lambda_max=float(lam),
        c_critical=float(c_critical),
        phi=phi,
        r_nodes_optical=r_nodes,
        iterations=iterations,
        converged=converged,
    )
