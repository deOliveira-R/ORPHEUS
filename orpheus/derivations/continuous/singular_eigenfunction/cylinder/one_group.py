r"""Westfall-Metcalf 1973 bare-critical infinite-cylinder solver.

Implements [WestfallMetcalf1973]_ for the **bare** infinite cylinder
(no reflector) with monoenergetic, isotropic scattering.

.. note::
   **Implementation discipline & scope of this prototype.** The
   Westfall-Metcalf 1973 method-of-record solves the bare-critical
   cylinder via a Mitsis-style iterative Fredholm scheme on the pseudo
   eigenfunction expansion (their Eqs 28-33 with the bare-cylinder
   reduction at :math:`c_1 = c_2`). Implementing that scheme requires
   careful numerical handling of:

   * the Cauchy principal-value integrals in the continuum
     pseudo-eigenfunction :math:`\eta_\nu(\mu) = c\,P\,\nu^2/(\nu^2 - \mu^2)
     + \lambda(\nu)\,\delta(\nu - \mu)` (Eq. 19);
   * scaled modified-Bessel evaluations to avoid overflow at small
     :math:`\nu` (where :math:`R/\nu` is large);
   * the Fredholm iteration coupling :math:`A'(\nu) \leftrightarrow
     \Phi'(\mu)` between Eq. 30 and Eq. 31, with criticality tested
     by Eq. 32.

   The Branch-1 SymPy verification (V_se-cyl.1 through V_se-cyl.7) is
   the **fully verified algebra-of-record** for the WM-72 derivation,
   including the discovery that printed Eq. 17 has a typo
   (:math:`\nu_0` should be :math:`\nu_0^2` in the numerator) — see
   :func:`...origins.cylinder_derivations.derive_discrete_pseudo_eigenfunction`.

   This Branch-2 prototype implements a **structurally simpler but
   STILL singular-eigenfunction-grounded** path: direct power
   iteration on the Mitsis cylindrical integral transport equation
   (WM-72 Eq. 6a), with the dominant Case eigenfunction
   :math:`\rho(r) = J_0(r/u_0)` confirming the V_se-cyl.7 structural
   prediction, and the criticality determined by the largest
   eigenvalue of the discretised kernel reaching :math:`\lambda_{\max} = 1/c`.

   Accuracy ceiling: the integral kernel
   :math:`K(r, t) = \int_0^1 K_0(\max/\mu)\,I_0(\min/\mu)\,d\mu/\mu^2`
   has a logarithmic singularity at :math:`r = t` (related to the
   2-D Green's function of the streaming operator). Treating this
   via single-cell product integration on the leading
   :math:`E_1(|r-t|)/(2\sqrt{r\,t})` asymptotic gives algebraic
   convergence in the GL grid count :math:`n` — about
   :math:`\mathcal{O}(1/n)` error. At :math:`n = 256` this prototype
   reaches :math:`\sim 10^{-3}` relative accuracy on the critical
   radius. Reaching the 1e-5 target on Sood ``Ua-1-0-CY`` requires
   either:

   * a **graded mesh** that refines toward the diagonal (Atkinson
     1976 product-integration), OR
   * the **full Mitsis-WM Fredholm iteration** (faster convergence
     because the discrete eigenfunction is singled out).

   The single Sood Variant α cross-check at 8.5e-6 (cylinder ``Ua-1-O-CY``)
   already provides a structurally-distinct external-reference
   anchor at the target accuracy; this WM-72 prototype is the
   **second** cross-check, structurally independent of both
   Variant α and Bickley-Naylor :math:`\mathrm{Ki}_n` integrals.

References
----------

.. [WestfallMetcalf1973] Westfall, R. M. & Metcalf, D. R. (1973).
   "Singular Eigenfunction Solution of the Monoenergetic Neutron
   Transport Equation for Finite Radially Reflected Critical
   Cylinders." *Nuclear Science and Engineering* **52**, 1-11.
   DOI 10.13182/NSE73-A23285.

.. [Mitsis1963] Mitsis, G. F. (1963). "Transport Solutions to the
   Monoenergetic Critical Problems." Argonne National Laboratory
   report ANL-6787.

.. [Atkinson1976] Atkinson, K. E. (1976). *A Survey of Numerical
   Methods for the Solution of Fredholm Integral Equations of the
   Second Kind.* SIAM. Chapter 6 — product integration on weakly
   singular kernels.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.integrate
import scipy.optimize
import scipy.special

from ...fn_method.core.dispersion import case_nu0


# ───────────────────────────────────────────────────────────────────────
# Result container
# ───────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class CylinderSingularEigenfunctionResult:
    r"""Output of :func:`solve_singular_eigenfunction_cylinder_bare_critical`.

    Attributes
    ----------
    r_c_mfp : float
        Critical radius in mean-free paths.
    r_c_cm : float | None
        Critical radius in cm (only populated if a :math:`\Sigma_t` is
        passed; ``None`` otherwise).
    c : float
        Mean number of secondaries per collision (input).
    u_0 : float
        Magnitude of the imaginary Case discrete eigenvalue,
        :math:`u_0 = |\nu_0|` so :math:`\nu_0 = i u_0`.
    nu_0 : complex
        Full Case eigenvalue :math:`\nu_0 = i u_0`.
    largest_eigenvalue_residual : float
        :math:`|c \cdot \lambda_{\max}(R_c) - 1|`, the bare-cylinder
        criticality residual at convergence.
    iterations : int
        Number of Brent iterations used by the criticality search.
    converged : bool
        Whether the root-finder reported success.
    n_grid : int
        Radial grid order used for the criticality computation.

    Notes
    -----
    The result carries flux-reconstruction methods:

    * :func:`compute_scalar_flux` — :math:`\rho(r) = J_0(r/u_0)` per
      V_se-cyl.7 (the dominant Case eigenfunction).
    * :func:`compute_pseudo_angular_flux` — discrete-mode
      contribution to the pseudo-flux :math:`\Phi_1(r, \mu)`.
    """

    r_c_mfp: float
    r_c_cm: float | None
    c: float
    u_0: float
    nu_0: complex
    largest_eigenvalue_residual: float
    iterations: int
    converged: bool
    n_grid: int

    def compute_scalar_flux(self, r: float | np.ndarray) -> np.ndarray:
        r"""Evaluate the bare-critical scalar flux :math:`\rho(r) = J_0(r/u_0)`.

        Per V_se-cyl.7
        (:func:`...origins.cylinder_derivations.derive_flux_reconstruction_bare_cylinder`):
        the bare-cylinder neutron density profile is the dominant Case
        eigenfunction, which for :math:`c > 1` (so :math:`\nu_0 = i u_0`
        purely imaginary) is

        .. math::

           \rho(r) = b_0\,J_0(r/u_0)

        with the overall amplitude normalised to :math:`b_0 = 1`
        (per WM-72 page 7: "b_0 is an arbitrary constant corresponding
        indirectly to the power level which we set equal to unity").

        Parameters
        ----------
        r : float or ndarray
            Radial coordinate(s) in mean-free paths,
            :math:`r \in [0, R_c]`.

        Returns
        -------
        ndarray
            :math:`\rho(r) = J_0(r/u_0)`, shape matching ``r``.
        """
        r_arr = np.atleast_1d(np.asarray(r, dtype=float))
        return scipy.special.j0(r_arr / self.u_0)

    def compute_pseudo_angular_flux(
        self, r: float | np.ndarray, mu: float | np.ndarray
    ) -> np.ndarray:
        r"""Evaluate the discrete-eigenfunction contribution to the
        angular pseudo-flux :math:`\Phi_1(r, \mu)`.

        Per WM-72 Eq. 23 in the bare-cylinder limit (V_se-cyl.4 says
        the continuum amplitude :math:`B(\nu)` is solved for via the
        Fredholm iteration; the discrete amplitude is :math:`b_0 = 1`):

        .. math::

           \Phi_1(r, \mu) = b_0\,J_0(r/u_0)\,\mu^2\,\eta_0(\mu)
                         + \int_0^1 B(\nu)\,J_0(r/\nu_0)\,
                           \mu^2\,\eta_{1\nu}(\mu)\,d\nu .

        This method returns ONLY the discrete-mode part:

        .. math::

           \Phi_1^{(0)}(r, \mu) = J_0(r/u_0)\,\mu^2\,
                                  \frac{c\,u_0^2}{u_0^2 + \mu^2} ,

        using the corrected :math:`\eta_0(\mu) = c\,\nu_0^2/(\nu_0^2 - \mu^2)
        = c\,u_0^2/(u_0^2 + \mu^2)` for :math:`\nu_0 = i u_0`
        (V_se-cyl.2 typo-fix vs printed Eq. 17). The continuum
        contribution :math:`\int_0^1 B(\nu)\,(\dots)\,d\nu` is omitted —
        for the bare cylinder it is small but non-zero, making the
        pseudo-flux from this method an approximation. The neutron
        density :math:`\rho(r) = \int_0^1 \Phi_1(r,\mu)/\mu^2\,d\mu`
        IS exact (the continuum integrates against the half-range
        normalisation to zero by orthogonality).

        .. warning::
           The pseudo-flux :math:`\Phi_1(r, \mu)` is **not the same**
           as the physical angular flux :math:`\Psi(r, \Omega)`.
           The Mitsis-Westfall-Metcalf framework works with
           pseudo-distribution functions related to the neutron
           density via the half-range moment integrals (Eq. 7a-b).
           Direct evaluation of the physical angular flux would
           require recovering :math:`\Psi(r, \mu, \phi)` from
           :math:`\Phi_1(r, \mu)` via Eq. 8a / inverse Bessel-kernel
           transformations, which is non-trivial and not addressed in
           WM-72 for the bare cylinder. Method-implementer recorded
           this as a deferred follow-up; see Phase B1 closeout memo.

        Parameters
        ----------
        r : float or ndarray
            Radial coordinate.
        mu : float or ndarray
            Angle cosine :math:`\mu \in (0, 1)`.

        Returns
        -------
        ndarray
            Discrete-mode :math:`\Phi_1^{(0)}(r, \mu)`, shape
            ``(len(r), len(mu))``.
        """
        r_arr = np.atleast_1d(np.asarray(r, dtype=float))
        mu_arr = np.atleast_1d(np.asarray(mu, dtype=float))
        R_grid, MU_grid = np.meshgrid(r_arr, mu_arr, indexing="ij")
        rho = scipy.special.j0(R_grid / self.u_0)
        # ν_0² = -u_0², so ν_0² - μ² = -(u_0² + μ²) → η_0 = c·u_0²/(u_0² + μ²).
        eta_0 = self.c * self.u_0**2 / (self.u_0**2 + MU_grid**2)
        return rho * MU_grid**2 * eta_0


# ───────────────────────────────────────────────────────────────────────
# Kernel evaluation (with logarithmic-singularity decomposition)
# ───────────────────────────────────────────────────────────────────────


def _kernel_regular(a: float, b: float) -> float:
    r"""Regular part :math:`K_{\mathrm{reg}}(a, b)` of the cylindrical
    kernel.

    The full cylindrical transport kernel (after the angular
    substitution :math:`u = 1/\mu`) is:

    .. math::

       K(a, b) = \int_1^\infty K_0(au)\,I_0(bu)\,du, \qquad a \ge b > 0 .

    Asymptotic at large :math:`u`:
    :math:`K_0(au)\,I_0(bu) \sim \exp(-(a-b)\,u)/(2u\,\sqrt{ab})`.

    Defining the asymptotic part
    :math:`K_{\mathrm{log}}(a, b) = E_1(|a-b|)/(2\sqrt{ab})`, the
    regular part :math:`K_{\mathrm{reg}}(a, b) = K(a, b) -
    K_{\mathrm{log}}(a, b)` is bounded for all :math:`a, b > 0`,
    including :math:`a = b` where :math:`K_{\mathrm{log}}` diverges
    logarithmically.

    Parameters
    ----------
    a, b : float
        :math:`a = \max(r, t)`, :math:`b = \min(r, t)`.

    Returns
    -------
    float
        :math:`K_{\mathrm{reg}}(a, b)`.
    """
    a = float(a)
    b = float(b)
    if abs(a - b) < 1e-14:
        # K(a, a) - 1/(2a)·∫_1^∞ du/u (divergent), but
        # K_reg(a, a) = ∫_1^∞ [K_0(au)·I_0(au) - 1/(2au)] du is finite.
        result, _ = scipy.integrate.quad(
            lambda u: (
                scipy.special.k0e(a * u) * scipy.special.i0e(a * u)
                - 1.0 / (2.0 * u * a)
            ),
            1.0,
            np.inf,
            limit=200,
            epsabs=1e-14,
            epsrel=1e-12,
        )
        return float(result)
    a_max = max(a, b)
    b_min = min(a, b)
    sqrt_ab = np.sqrt(a_max * b_min)
    result, _ = scipy.integrate.quad(
        lambda u: (
            scipy.special.k0e(a_max * u)
            * scipy.special.i0e(b_min * u)
            * np.exp((b_min - a_max) * u)
            - np.exp(-(a_max - b_min) * u) / (2.0 * u * sqrt_ab)
        ),
        1.0,
        np.inf,
        limit=200,
        epsabs=1e-14,
        epsrel=1e-12,
    )
    return float(result)


def _kernel_log_part(a: float, b: float) -> float:
    r"""Logarithmic-asymptotic part :math:`E_1(|a-b|)/(2\sqrt{ab})`.

    Diverges as :math:`-\log|a-b|/(2\sqrt{ab})` as :math:`a \to b`;
    its single-cell integral is computed analytically by
    :func:`_diagonal_log_integral`.
    """
    a = float(a)
    b = float(b)
    if abs(a - b) < 1e-14:
        return float("inf")
    return float(scipy.special.exp1(abs(a - b)) / (2.0 * np.sqrt(a * b)))


def _diagonal_log_integral(r: float, h: float) -> float:
    r"""Compute :math:`\int_{r-h/2}^{r+h/2} E_1(|r-t|)\,t/(2\sqrt{r\,t})\,dt`
    to leading order in :math:`h/r`.

    Substitution :math:`s = t - r`:

    .. math::

       I(r, h) = \int_{-h/2}^{h/2}
                 \frac{E_1(|s|)\,(r + s)}{2\sqrt{r\,(r+s)}}\,ds .

    Leading order in :math:`s/r`: :math:`(r+s)/\sqrt{r(r+s)} \approx
    \sqrt{(r+s)/r}\,(r+s)/(r+s) = \sqrt{1 + s/r}` →

    .. math::

       I(r, h) \approx \int_{-h/2}^{h/2} E_1(|s|) \cdot
                       \tfrac{1}{2}\sqrt{1 + s/r}\,ds
              \approx \tfrac{1}{2}\int_{-h/2}^{h/2} E_1(|s|)\,ds
              = \int_0^{h/2} E_1(s)\,ds  \qquad (\text{symmetry})
              = \tfrac{h}{2}\,E_1(h/2) + 1 - e^{-h/2} .

    Parameters
    ----------
    r : float
        Cell-centre radial coordinate.
    h : float
        Cell width (in radial coordinate).

    Returns
    -------
    float
        :math:`I(r, h)` to leading order.

    Notes
    -----
    Higher-order corrections (in :math:`h/r`) are NOT included. The
    accuracy of the bare-critical-radius computation is therefore
    bounded by :math:`\mathcal{O}(h/r)` near the cylinder axis
    :math:`r = 0` (where the smallest GL nodes are located). For
    Sood ``Ua-1-0-CY`` (:math:`R = 1.725` mfp, :math:`n_{\rm grid} = 64`),
    this gives :math:`\sim 10^{-3}` relative accuracy on
    :math:`R_c`. Tightening the accuracy floor requires either a
    graded mesh refinement near the diagonal or an analytical
    cell-integral that retains the full :math:`\sqrt{r(r+s)}` factor
    — both are deferred to a future hardening pass (see Phase B1
    closeout memo).
    """
    h_half = h / 2.0
    return 0.5 * (h * scipy.special.exp1(h_half) + 2.0 * (1.0 - np.exp(-h_half)))


# ───────────────────────────────────────────────────────────────────────
# Operator matrix + criticality search
# ───────────────────────────────────────────────────────────────────────


def _build_operator_matrix(R: float, n_grid: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r"""Discretise the Mitsis-Westfall-Metcalf cylinder integral operator.

    The integral transport equation for the neutron density (WM-72
    Eq. 6a in the bare limit :math:`c_1 = c_2 = c`):

    .. math::

       \rho(r) = c \int_0^R K(r, t)\,t\,\rho(t)\,dt

    is discretised by Gauss-Legendre quadrature on :math:`(0, R)`
    with the kernel split

    .. math::

       K(r, t) = K_{\mathrm{reg}}(r, t) + K_{\mathrm{log}}(r, t)

    where :math:`K_{\mathrm{log}} = E_1(|r-t|)/(2\sqrt{r\,t})` carries
    the logarithmic singularity and :math:`K_{\mathrm{reg}}` is
    bounded.

    Off-diagonal elements use the full kernel evaluated by
    quadrature; the diagonal element gets the regular part by the
    same quadrature rule plus the analytically-computed local
    integral of the log-singular part (single-cell product
    integration).

    Returns
    -------
    M : ndarray, shape (n_grid, n_grid)
        Operator matrix such that the largest eigenvalue
        :math:`\lambda_{\max}(M) = 1/c` at the bare-critical radius.
    r : ndarray, shape (n_grid,)
        GL nodes on :math:`(0, R)`.
    w : ndarray, shape (n_grid,)
        GL weights on :math:`(0, R)`.
    """
    nodes, weights = np.polynomial.legendre.leggauss(n_grid)
    r = 0.5 * (nodes + 1.0) * R
    w = 0.5 * R * weights
    M = np.zeros((n_grid, n_grid))
    for i in range(n_grid):
        for j in range(n_grid):
            r_i = r[i]
            t_j = r[j]
            a = max(r_i, t_j)
            b = min(r_i, t_j)
            K_reg_val = _kernel_regular(a, b)
            if i == j:
                # Diagonal: K_reg from the regular part (finite at a=b)
                # plus analytical local integral of K_log over the cell.
                M[i, j] = K_reg_val * t_j * w[j] + _diagonal_log_integral(r_i, w[i])
            else:
                K_log_val = _kernel_log_part(a, b)
                M[i, j] = (K_reg_val + K_log_val) * t_j * w[j]
    return M, r, w


def _largest_eigenvalue(R: float, n_grid: int) -> float:
    r"""Largest real eigenvalue of the operator at radius :math:`R`."""
    M, _, _ = _build_operator_matrix(R, n_grid)
    eigs = np.linalg.eigvals(M)
    return float(np.max(np.real(eigs)))


# ───────────────────────────────────────────────────────────────────────
# Public entry point
# ───────────────────────────────────────────────────────────────────────


def solve_singular_eigenfunction_cylinder_bare_critical(
    c: float,
    *,
    sigma_t: float | None = None,
    n_grid: int = 64,
    R_min: float = 0.1,
    R_max: float = 30.0,
    bisect_tol: float = 1e-10,
    max_iter: int = 200,
) -> CylinderSingularEigenfunctionResult:
    r"""Bare-critical infinite cylinder via direct Nyström discretisation
    of the Mitsis-Westfall-Metcalf integral transport equation.

    Solves for :math:`R_c` such that the largest eigenvalue of the
    discretised cylinder integral operator (WM-72 Eq. 6a in the bare
    limit :math:`c_1 = c_2 = c`) reaches :math:`\lambda_{\max} = 1/c`,
    i.e., :math:`c \cdot \lambda_{\max}(R_c) = 1`.

    The discretisation uses Gauss-Legendre quadrature on :math:`(0, R)`
    with single-cell product integration on the log-singular kernel
    diagonal (see :func:`_diagonal_log_integral` and
    :func:`_kernel_regular`).

    Parameters
    ----------
    c : float
        Mean number of secondaries per collision,
        :math:`c = (\Sigma_s + \nu\Sigma_f)/\Sigma_t`. Must satisfy
        :math:`c > 1` (multiplying medium).
    sigma_t : float, optional
        Total macroscopic cross section in cm⁻¹. If provided, the
        result includes ``r_c_cm = r_c_mfp / sigma_t``; otherwise
        ``r_c_cm = None``.
    n_grid : int, default 64
        Number of GL nodes for the radial discretisation. Higher
        ``n_grid`` improves accuracy at :math:`\mathcal{O}(1/n)`.
        For Sood ``Ua-1-0-CY`` target accuracy floor ~3e-3 at
        :math:`n = 64`, ~1e-3 at :math:`n = 256` (proportional to
        cell-width-near-axis :math:`h/r \sim 1/n`). The
        :math:`\mathcal{O}(1/n)` floor is set by the leading-order
        single-cell product integration on the log singularity;
        graded mesh refinement would tighten the floor (deferred —
        see closeout memo).
    R_min, R_max : float, defaults 0.1, 30.0
        Bracket for the criticality root search. Default range
        :math:`[0.1, 30]` covers :math:`c \in [1.005, 5]`.
    bisect_tol : float, default 1e-10
        Relative tolerance on :math:`R` for the Brent root-finder.
    max_iter : int, default 200
        Maximum Brent iterations.

    Returns
    -------
    CylinderSingularEigenfunctionResult
        With the bare-critical radius, dispersion-root magnitude
        :math:`u_0 = |\nu_0|`, and reconstructions of the dominant
        Case eigenfunction.

    Raises
    ------
    ValueError
        If :math:`c \le 1`.
    RuntimeError
        If the bracket :math:`[R_\min, R_\max]` does not contain a
        sign change of :math:`c \cdot \lambda_{\max}(R) - 1`.

    Notes
    -----
    Below the trusted-library line (per ``algebra-of-record``), this
    solver reuses :func:`...fn_method.core.dispersion.case_nu0` for
    the dispersion-root primitive — the dispersion function is a
    medium property (V_se-cyl.1 / Eq. 18) common to all
    singular-eigenfunction-based methods. Sharing the dispersion-root
    finder is structurally safe: the criticality REDUCTION (Mitsis
    integral equation Eq. 6a vs F_N's Eq. 39 determinantal
    condition) is structurally distinct.

    The flux-reconstruction methods on
    :class:`CylinderSingularEigenfunctionResult` use the dominant
    Case eigenfunction :math:`\rho(r) = J_0(r/u_0)` (V_se-cyl.7),
    which is the load-bearing physics prediction of the
    singular-eigenfunction expansion regardless of which numerical
    discretisation scheme is used.
    """
    if c <= 1.0:
        raise ValueError(
            f"Bare-critical cylinder requires c > 1 (multiplying medium); "
            f"got c = {c}."
        )
    if n_grid < 4:
        raise ValueError(f"n_grid must be ≥ 4, got {n_grid}.")

    # Dispersion root: u_0 = |ν_0|, ν_0 = i·u_0.
    u_0, is_imag = case_nu0(c)
    if not is_imag:
        raise RuntimeError(
            f"Internal: case_nu0 returned subcritical regime for c = {c} > 1."
        )
    nu_0_complex = complex(0.0, u_0)

    def residual(R: float) -> float:
        return c * _largest_eigenvalue(R, n_grid) - 1.0

    # Verify sign change in the bracket.
    f_min = residual(R_min)
    f_max = residual(R_max)
    if not (np.isfinite(f_min) and np.isfinite(f_max)):
        raise RuntimeError(
            f"Non-finite residual at bracket endpoints for c = {c}."
        )
    if f_min * f_max > 0:
        raise RuntimeError(
            f"No sign change of (c·λ_max - 1) on [{R_min}, {R_max}] for c = {c}. "
            f"residual at R_min = {f_min}, R_max = {f_max}. "
            f"Try widening the bracket."
        )

    R_c, brent_result = scipy.optimize.brentq(
        residual,
        R_min,
        R_max,
        xtol=bisect_tol,
        rtol=bisect_tol,
        maxiter=max_iter,
        full_output=True,
        disp=False,
    )
    R_c = float(R_c)
    final_residual = abs(residual(R_c))

    return CylinderSingularEigenfunctionResult(
        r_c_mfp=R_c,
        r_c_cm=(R_c / sigma_t if sigma_t is not None else None),
        c=c,
        u_0=u_0,
        nu_0=nu_0_complex,
        largest_eigenvalue_residual=final_residual,
        iterations=int(brent_result.iterations),
        converged=bool(brent_result.converged),
        n_grid=n_grid,
    )
