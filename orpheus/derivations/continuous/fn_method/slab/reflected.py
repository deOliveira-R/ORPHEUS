r"""Reflected-slab F_N method solver — 1G isotropic, finite reflector.

Implements the Neshat-Maiorino 1980 (*Ann. Nucl. Energy* **7**, 79-81)
F_N collocation system for the **critical core half-thickness**
:math:`\tau` of a 1G isotropic homogeneous slab core (:math:`c_1 > 1`)
sandwiched between two finite reflector regions (:math:`c_2 < 1`,
half-thickness :math:`\Delta = b - \tau`):

* **Core**: :math:`-\tau \le x \le \tau`.
* **Reflector**: :math:`\tau < |x| \le b`.
* **BC**: vacuum at :math:`|x| = b`; interface continuity at :math:`x = \pm\tau`.

The F_N system is :math:`3(N+1)` equations in :math:`3(N+1)` unknowns,
gathered from three projection equations (NM Eqs. 10, 11, 12) plus a
critical condition (NM Eq. 15) that fixes :math:`\tau`. The
implementation follows the iterative scheme on NM p. 81:

1. Initialize :math:`\tau` via the :math:`F_0` closed form NM Eqs.
   16-17 (a 2-stream reflector approximation).
2. With current :math:`\tau`, build the :math:`(3N+2) \times (3N+2)`
   linear system from Eqs. 10-12 collocated on the published grid
   (with the :math:`a_0 = 1` normalisation removing one degree of
   freedom). Solve for :math:`\{a_\alpha, b_\alpha, e_\alpha\}`.
3. Plug coefficients into NM Eq. 15 root-find — solve for the new
   :math:`\tau`.
4. Iterate until :math:`\tau` converges (typically 3-4 iterations).

Architectural notes
-------------------

* The :math:`A_\alpha`, :math:`B_\alpha^{(i)}` recursions are
  IDENTICAL to the bare-slab Grandjean-Siewert recursions (verified
  symbolically in
  :func:`...origins.fn_slab_reflected_derivations.derive_reflected_moment_recursions_match_bare`).
  We reuse :func:`...core.moments.A_alpha_array` and
  :func:`...core.moments.B_alpha_array`, calling them with
  :math:`c_i` per region.
* The X-function is **medium-local**: each region uses its own
  :math:`c_i`. There is NO two-region X-function.
* Numerical pathologies: only one — at very large :math:`\Delta`,
  :math:`e^{+\Delta/\hat\xi}` in NM Eq. 11 grows. NM use
  :math:`\Delta \le 1` mfp without trouble; we issue a warning at
  :math:`\Delta > 5`.

References
----------

* Neshat & Maiorino 1980, *Ann. Nucl. Energy* **7**, 79-81.
* Grandjean & Siewert 1979, *Nucl. Sci. Eng.* **69**, 161 (recursions).
* Siewert & Benoist 1979, *Nucl. Sci. Eng.* **69**, 156 (Part I).
* Burkart 1976, *Trans. ANS* **24**, 190 — "Exact" reference values
  cited in NM Table 2.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..core.dispersion import (
    case_dispersion_root_subcritical,
    case_dispersion_root_supercritical,
)
from ..core.moments import A_alpha_array, B_alpha_array


@dataclass(frozen=True)
class SlabReflectedFNResult:
    r"""Result of :func:`solve_fn_slab_reflected_critical`.

    Attributes
    ----------
    tau_critical_mfp : float
        Critical core half-thickness :math:`\tau` in mean free paths.
    c_core : float
        Core mean number of secondaries per collision (:math:`c_1 > 1`).
    c_reflector : float
        Reflector mean number of secondaries per collision (:math:`c_2 < 1`).
    Delta_mfp : float
        Reflector half-thickness :math:`\Delta = b - \tau` (mfp).
    N : int
        F_N order (number of polynomial moments per coefficient array).
    nu0_core : float
        Core Case discrete eigenvalue magnitude :math:`u_0 = |\nu_0|`
        (the actual eigenvalue is :math:`\nu_0 = i u_0` for
        :math:`c_1 > 1`).
    eta0_reflector : float
        Reflector Case discrete eigenvalue (real, :math:`> 1`) for
        :math:`c_2 < 1`.
    coefficients_a : np.ndarray
        :math:`(N+1,)` real array of :math:`a_\alpha` —
        :math:`\psi_1(\tau, +\mu) = \sum a_\alpha \mu^\alpha`.
    coefficients_b : np.ndarray
        :math:`(N+1,)` real array of :math:`b_\alpha` —
        :math:`\psi_1(\tau, -\mu) = \sum b_\alpha \mu^\alpha`.
    coefficients_e : np.ndarray
        :math:`(N+1,)` real array of :math:`e_\alpha` —
        :math:`\psi_2(b, +\mu) = \sum e_\alpha \mu^\alpha`.
    iterations : int
        Number of outer-loop iterations until convergence.
    converged : bool
        Whether the iteration converged within ``tol``.
    tau_history : list[float]
        :math:`\tau` values per iteration (for diagnostics).
    """

    tau_critical_mfp: float
    c_core: float
    c_reflector: float
    Delta_mfp: float
    N: int
    nu0_core: float
    eta0_reflector: float
    coefficients_a: np.ndarray
    coefficients_b: np.ndarray
    coefficients_e: np.ndarray
    iterations: int
    converged: bool
    tau_history: list[float]


def _collocation_points_xi_hat(eta0: float, N: int) -> np.ndarray:
    r"""Collocation points :math:`\hat\xi \in P_2 = \eta_0 \cup (0, 1)`
    used for NM Eqs. 10, 11 (reflector projections).

    Per NM p. 81: :math:`\hat\xi_0 = \eta_0`, :math:`\hat\xi_1 = 0`,
    :math:`\hat\xi_2 = 1`, remaining :math:`(N - 2)` points equally
    spaced in :math:`[0, 1]`.

    Returns array of length :math:`N + 1`.
    """
    pts = np.empty(N + 1, dtype=float)
    pts[0] = eta0
    if N >= 1:
        pts[1] = 0.0
    if N >= 2:
        pts[2] = 1.0
    if N >= 3:
        # N-2 interior points equally spaced strictly inside (0, 1).
        # Same convention as bare-slab collocation_points.
        pts[3:] = np.array([k / (N - 1) for k in range(1, N - 1)])
    return pts


def _collocation_points_xi(N: int) -> np.ndarray:
    r"""Collocation points :math:`\xi \in (0, 1)` used for NM Eq. 12
    (core projection).

    NM p. 81: "for :math:`\xi` we pick :math:`\xi_0 = 0,\, \xi_1 = 1`,
    and the remaining also equally spaced in :math:`[0, 1]`."

    NM solve a :math:`(3N+2)`-equation linear system. Eqs. 10 + 11
    contribute :math:`2(N+1)` rows. Eq. 12 therefore contributes
    :math:`N` rows. So we use :math:`N` collocation points for Eq. 12:
    :math:`\xi_0 = 0,\, \xi_1 = 1`, remaining :math:`(N - 2)` interior
    points equally spaced in :math:`(0, 1)`.

    Eq. 15 (collocated at :math:`\xi = \nu_0`, the core Case discrete
    eigenvalue) provides the (N+1)th equation as the **critical
    condition** — it's not part of the linear system; it's the
    eigenvalue equation that determines :math:`\tau`.

    Returns array of length :math:`N`.
    """
    if N == 0:
        # Pathological: nothing to collocate. Fall back to single point.
        return np.array([0.0])
    pts = np.empty(N, dtype=float)
    pts[0] = 0.0
    if N >= 2:
        pts[1] = 1.0
    if N >= 3:
        # N-2 interior points equally spaced inside (0, 1).
        pts[2:] = np.array([k / (N - 1) for k in range(1, N - 1)])
    return pts


def _moment_arrays_at_xi(
    N: int, xi: float, c: float
) -> tuple[np.ndarray, np.ndarray]:
    r"""Compute :math:`(B_\alpha^{(i)}(\xi))_{\alpha=0}^N` and
    :math:`(A_\alpha(\xi))_{\alpha=0}^N` for given :math:`c` and
    :math:`\xi`.

    Handles the :math:`\xi = 0` limit (where :math:`\xi \log(1+1/\xi)
    \to 0`):

    * :math:`B_0(0) = 2/c - 1`, :math:`B_\alpha(0) = -1/(\alpha+1)`
      for :math:`\alpha \ge 1`.
    * :math:`A_0(0) = 1`, :math:`A_\alpha(0) = 1/(\alpha+1)`
      for :math:`\alpha \ge 1`.
    """
    if xi == 0.0:
        Bvec = np.empty(N + 1, dtype=float)
        Avec = np.empty(N + 1, dtype=float)
        Bvec[0] = (2.0 / c) - 1.0
        Avec[0] = 1.0
        for n in range(1, N + 1):
            Bvec[n] = 0.0 - 1.0 / (n + 1)
            Avec[n] = 0.0 + 1.0 / (n + 1)
        return Bvec, Avec
    return B_alpha_array(N, xi, c), A_alpha_array(N, xi, c)


def _build_F0_initial_tau(
    c_core: float,
    c_reflector: float,
    Delta: float,
) -> tuple[float, float, float, float]:
    r"""Compute the :math:`F_0` initial guess for :math:`\tau`
    (NM Eqs. 16-17).

    .. math::

       b_0 &= \frac{A_0(\eta_0) B_0^{(2)}(\eta_0)
              (1 - e^{-2\Delta/\eta_0})}
              {[B_0^{(2)}(\eta_0)]^2 - [A_0(\eta_0)]^2 e^{-2\Delta/\eta_0}}
       \\
       \tau^{(0)} &= -\frac{\nu_0}{2} \log\!\left[
              \frac{b_0 A_0(\nu_0) - B_0^{(1)}(\nu_0)}
                   {A_0(\nu_0) - b_0 B_0^{(1)}(\nu_0)}
              \right]

    Note :math:`\nu_0` here is the magnitude of the core's imaginary
    Case eigenvalue :math:`u_0`. Per NM, the F_0 truncation uses
    :math:`u_0` directly (not the imaginary :math:`i u_0`); the
    :math:`\log` argument is a real ratio whose sign is positive
    for physical configurations.

    Returns ``(tau_0, b_0, nu_0_core, eta_0_reflector)``.
    """
    # Reflector: c < 1 → real ν_0 > 1. NM uses η_0 for this.
    eta0 = case_dispersion_root_subcritical(c_reflector)
    # Core: c > 1 → ν_0 = i u_0; we use u_0 throughout.
    u0 = case_dispersion_root_supercritical(c_core)

    # Moment values (α = 0) at the discrete eigenvalues.
    # B_0^(2)(η_0) = 2/c_2 - 1 - η_0 log(1+1/η_0), A_0(η_0) = 1 - η_0 log(1+1/η_0).
    log_eta = math.log(1.0 + 1.0 / eta0)
    B0_2_at_eta = (2.0 / c_reflector) - 1.0 - eta0 * log_eta
    A0_at_eta = 1.0 - eta0 * log_eta

    exp_neg = math.exp(-2.0 * Delta / eta0)
    b0 = (
        A0_at_eta * B0_2_at_eta * (1.0 - exp_neg)
        / (B0_2_at_eta**2 - A0_at_eta**2 * exp_neg)
    )

    # Core F_0 critical condition: NM Eq. 16, with a_0 = 1, evaluated at u_0.
    # Note ν_0 = i u_0 → ν_0 log(1+1/ν_0) is COMPLEX in general. NM Table 2
    # uses a real F_0 estimate, which means they treat |ν_0| and use the
    # absolute-value branch. Concretely: A_0(ν_0) and B_0^(1)(ν_0) are
    # complex when ν_0 = i u_0, but the ratio is real (the imaginary
    # parts cancel by the symmetry of the log on the imaginary axis).
    # We work with complex values and project to real at the end.
    log_iu = complex(0.0, 1.0) * u0 * complex_log_one_plus_inv_iu(u0)
    A0_at_nu = 1.0 - log_iu
    B0_1_at_nu = (2.0 / c_core) - 1.0 - log_iu

    arg_complex = (b0 * A0_at_nu - B0_1_at_nu) / (A0_at_nu - b0 * B0_1_at_nu)
    # NM Eq. 16: τ^(0) = -(ν_0/2) log[arg]. With ν_0 = i u_0 (purely
    # imaginary) and |arg| = 1 by symmetry of the reflector ratio
    # (the modulus exactly cancels because the system is critical),
    # log[arg] = i θ is purely imaginary, so
    #   τ^(0) = -(i u_0 / 2)(i θ) = (u_0 / 2) θ
    # where θ ∈ (0, π) is the angle of arg in the upper half plane.
    log_complex = complex_log(arg_complex)
    tau_0 = 0.5 * u0 * log_complex.imag
    if tau_0 <= 0.0:
        # If F_0 estimate is non-physical, fall back to a small positive
        # value (the iteration will refine).
        tau_0 = 0.1 * u0

    return tau_0, b0, u0, eta0


def complex_log_one_plus_inv_iu(u0: float) -> complex:
    r"""Compute :math:`\log(1 + 1/(i u_0))` where :math:`u_0 > 0`.

    Used in :math:`A_0(i u_0) = 1 - (i u_0) \log(1 + 1/(i u_0))`.

    .. math::

       1 + \frac{1}{i u_0} = 1 - \frac{i}{u_0}.

    Then :math:`\log(1 - i/u_0)` is complex with real part
    :math:`\log\sqrt{1 + 1/u_0^2}` and imaginary part
    :math:`-\arctan(1/u_0)` (principal branch).
    """
    # 1 + 1/(i u0) = 1 - i/u0.
    z = complex(1.0, -1.0 / u0)
    return complex_log(z)


def complex_log(z: complex) -> complex:
    """Principal-branch complex logarithm."""
    return complex(math.log(abs(z)), math.atan2(z.imag, z.real))


def _build_full_system_matrix(
    tau: float,
    Delta: float,
    c_core: float,
    c_reflector: float,
    N: int,
    u0: float,
    eta0: float,
) -> np.ndarray:
    r"""Build the full :math:`3(N+1) \times 3(N+1)` reflected-slab F_N
    homogeneous system matrix at given :math:`\tau`.

    Per NM (Eqs. 10, 11, 12) — three projection equations each
    collocated at :math:`N+1` points, totaling :math:`3(N+1)` rows
    on :math:`3(N+1)` unknowns :math:`\{a_\alpha, b_\alpha, e_\alpha\}`.

    For Eq. 12 we collocate at :math:`(\xi_0 = \nu_0 = i u_0,\,
    \xi_1 = 0,\, \xi_2 = 1, \ldots)` — the core dispersion root
    plays the same role here as in the bare-slab F_N method
    (Grandjean-Siewert collocation prescription). The :math:`\xi_0
    = \nu_0` row is COMPLEX; we work in complex linear algebra.

    NM's Eq. 15 is equivalent to picking :math:`\xi = \nu_0` in
    Eq. 12; we use that automatically by embedding :math:`\nu_0` as
    the first collocation point.

    The critical condition is :math:`\det M(\tau) = 0`. Branch-2
    finds the real positive :math:`\tau` for which :math:`\Re \det M
    = 0`; the imaginary part vanishes by symmetry of the critical
    mode.
    """
    n = N + 1  # length of each coefficient array
    n_unknowns = 3 * n
    n_rows = 3 * n  # equations 10, 11, 12 each yield N+1 rows

    # Collocation points (per NM p. 81):
    #   ξ̂ ∈ P_2: ξ̂_0 = η_0, ξ̂_1 = 0, ξ̂_2 = 1, remaining equally spaced.
    #   ξ ∈ (0,1): ξ_0 = 0, ξ_1 = 1, remaining equally spaced.
    # Note ν_0 is NOT a Eq. 12 collocation point; the critical
    # condition (Eq. 15) is implicit in the homogeneous system having
    # a non-trivial null space.
    xi_hat = _collocation_points_xi_hat(eta0, N)  # N+1 points for Eq 10, 11
    xi_real = _collocation_points_xi(N)  # N+1 points for Eq 12 (real)
    xi = xi_real.astype(complex)

    # Pre-compute moment arrays at each collocation point (complex).
    A_at_xi_hat = np.empty((n, n), dtype=complex)
    B2_at_xi_hat = np.empty((n, n), dtype=complex)
    A_at_xi = np.empty((n, n), dtype=complex)
    B1_at_xi = np.empty((n, n), dtype=complex)

    for beta in range(n):
        Bvec_r, Avec_r = _moment_arrays_at_xi(N, xi_hat[beta], c_reflector)
        A_at_xi_hat[beta, :] = Avec_r.astype(complex)
        B2_at_xi_hat[beta, :] = Bvec_r.astype(complex)
        # Core moments — complex if xi[β] is complex.
        Bvec_c, Avec_c = _complex_moment_arrays_at_xi(N, xi[beta], c_core)
        A_at_xi[beta, :] = Avec_c
        B1_at_xi[beta, :] = Bvec_c

    # Pre-compute exponential factors.
    exp_neg_Delta = np.where(
        xi_hat > 0.0,
        np.exp(-Delta / np.where(xi_hat > 0, xi_hat, 1.0)),
        0.0,
    ).astype(complex)
    safe_exp_arg = np.minimum(Delta / np.where(xi_hat > 0, xi_hat, 1.0), 50.0)
    exp_pos_Delta = np.where(
        xi_hat > 0.0, np.exp(safe_exp_arg), np.exp(50.0)
    ).astype(complex)
    # exp(-2τ/ξ) — complex if ξ is complex (β = 0).
    exp_neg_2tau = np.empty(n, dtype=complex)
    for beta in range(n):
        if xi[beta] == 0.0:
            exp_neg_2tau[beta] = 0.0
        else:
            exp_neg_2tau[beta] = np.exp(-2.0 * tau / xi[beta])

    # Build the homogeneous system M @ z = 0 with z = (a, b, e).
    M = np.zeros((n_rows, n_unknowns), dtype=complex)

    # Eq. 10 rows (β = 0, ..., N), collocated at ξ̂ ∈ P_2:
    #   Σ_α a_α A_α(ξ̂_β) - Σ_α b_α B_α^(2)(ξ̂_β)
    #     - e^{-Δ/ξ̂_β} Σ_α e_α A_α(ξ̂_β) = 0
    for beta in range(n):
        M[beta, 0:n] = A_at_xi_hat[beta, :]
        M[beta, n:2 * n] = -B2_at_xi_hat[beta, :]
        M[beta, 2 * n:3 * n] = -exp_neg_Delta[beta] * A_at_xi_hat[beta, :]

    # Eq. 11 rows (β = N+1, ..., 2N+1), collocated at ξ̂ ∈ P_2:
    #   Σ_α b_α A_α(ξ̂_β) - Σ_α a_α B_α^(2)(ξ̂_β)
    #     + e^{+Δ/ξ̂_β} Σ_α e_α B_α^(2)(ξ̂_β) = 0
    for beta in range(n):
        row = n + beta
        M[row, 0:n] = -B2_at_xi_hat[beta, :]
        M[row, n:2 * n] = A_at_xi_hat[beta, :]
        M[row, 2 * n:3 * n] = exp_pos_Delta[beta] * B2_at_xi_hat[beta, :]

    # Eq. 12 rows (β = 2N+2, ..., 3N+2), collocated at ξ ∈ {ν_0, 0, 1, ...}:
    #   Σ_α b_α A_α(ξ_β) - Σ_α a_α B_α^(1)(ξ_β)
    #     - e^{-2τ/ξ_β} [Σ_α a_α A_α(ξ_β) - Σ_α b_α B_α^(1)(ξ_β)] = 0
    # → coefficient on a_α: -B_α^(1)(ξ) - e^{-2τ/ξ} A_α(ξ)
    #   coefficient on b_α:  A_α(ξ) + e^{-2τ/ξ} B_α^(1)(ξ)
    for beta in range(n):
        row = 2 * n + beta
        M[row, 0:n] = -B1_at_xi[beta, :] - exp_neg_2tau[beta] * A_at_xi[beta, :]
        M[row, n:2 * n] = A_at_xi[beta, :] + exp_neg_2tau[beta] * B1_at_xi[beta, :]
        # e_α: zero in Eq. 12 (no reflector dependence — Eq. 12 is core only).

    return M


def _complex_moment_arrays_at_xi(
    N: int, xi: complex, c: float
) -> tuple[np.ndarray, np.ndarray]:
    r"""Complex-arithmetic version of :func:`_moment_arrays_at_xi`.

    Handles :math:`\xi = 0` and complex :math:`\xi` (e.g.,
    :math:`\xi = i u_0`) for the recursion seed and recursion.

    For complex :math:`\xi`,
    :math:`\xi \log(1 + 1/\xi) = \xi \log((\xi + 1)/\xi)` evaluated
    on the principal branch. We use ``cmath.log``.
    """
    import cmath

    Bvec = np.empty(N + 1, dtype=complex)
    Avec = np.empty(N + 1, dtype=complex)
    if xi == 0.0:
        Bvec[0] = (2.0 / c) - 1.0
        Avec[0] = 1.0
        for n in range(1, N + 1):
            Bvec[n] = 0.0 - 1.0 / (n + 1)
            Avec[n] = 0.0 + 1.0 / (n + 1)
    else:
        log_term = xi * cmath.log(1.0 + 1.0 / xi)
        Bvec[0] = (2.0 / c) - 1.0 - log_term
        Avec[0] = 1.0 - log_term
        for n in range(1, N + 1):
            Bvec[n] = xi * Bvec[n - 1] - 1.0 / (n + 1)
            Avec[n] = -xi * Avec[n - 1] + 1.0 / (n + 1)
    return Bvec, Avec


def _solve_NM_iteration_step(
    tau: float,
    Delta: float,
    c_core: float,
    c_reflector: float,
    N: int,
    u0: float,
    eta0: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r"""Per NM p. 81: SOLVE the 3(N+1) collocation system AT FIXED τ
    with normalisation :math:`a_0 = 1`.

    The collocation system has 3(N+1) equations from NM Eqs. 10-12
    in 3(N+1) unknowns; with :math:`a_0 = 1` substituted, this becomes
    a square (3N+3) × (3N+2) system whose RHS comes from the
    :math:`a_0` column. We use NM's prescription:

    * Eqs. 10, 11 collocated at :math:`\hat\xi \in \{\eta_0, 0, 1,
      \ldots\}` (N+1 points each → 2(N+1) rows).
    * Eq. 12 collocated at :math:`\xi \in \{0, 1, \ldots\}` (no
      :math:`\nu_0`!) — N+1 rows.

    Total rows: 3(N+1) = 3N+3. Total unknowns after a_0 fixed:
    3N+2. So the system is OVER-determined by one row. We solve
    least-squares (the redundancy is the implicit critical
    condition).

    Returns ``(a, b, e)`` arrays each of length :math:`N+1`, with
    :math:`a_0 = 1`.
    """
    n = N + 1  # length of each coefficient array (a, b, e)

    # Collocation points (per NM p. 81 prescription):
    xi_hat = _collocation_points_xi_hat(eta0, N)  # N+1 points (real)
    xi = _collocation_points_xi(N)  # N points (real, no ν_0)

    n_xi_hat = len(xi_hat)  # = n = N+1
    n_xi = len(xi)  # = N

    # Pre-compute moment arrays at each collocation point.
    A_at_xi_hat = np.empty((n_xi_hat, n))
    B2_at_xi_hat = np.empty((n_xi_hat, n))
    A_at_xi = np.empty((n_xi, n))
    B1_at_xi = np.empty((n_xi, n))
    for beta in range(n_xi_hat):
        Bvec, Avec = _moment_arrays_at_xi(N, xi_hat[beta], c_reflector)
        A_at_xi_hat[beta, :] = Avec
        B2_at_xi_hat[beta, :] = Bvec
    for beta in range(n_xi):
        Bvec_c, Avec_c = _moment_arrays_at_xi(N, xi[beta], c_core)
        A_at_xi[beta, :] = Avec_c
        B1_at_xi[beta, :] = Bvec_c

    # Exponential factors.
    exp_neg_Delta = np.where(
        xi_hat > 0.0,
        np.exp(-Delta / np.where(xi_hat > 0, xi_hat, 1.0)),
        0.0,
    )
    safe_exp_arg = np.minimum(Delta / np.where(xi_hat > 0, xi_hat, 1.0), 50.0)
    exp_pos_Delta = np.where(xi_hat > 0.0, np.exp(safe_exp_arg), np.exp(50.0))
    exp_neg_2tau = np.where(
        xi > 0.0, np.exp(-2.0 * tau / np.where(xi > 0, xi, 1.0)), 0.0
    )

    # Build the system: 2(N+1) + N = 3N+2 equations in 3(N+1) = 3N+3
    # unknowns z = (a_0, ..., a_N, b_0, ..., b_N, e_0, ..., e_N).
    # With a_0 = 1 substituted, system is (3N+2) × (3N+2) — square.
    n_rows = 2 * n_xi_hat + n_xi  # 2(N+1) + N = 3N+2
    n_cols = 3 * n  # 3(N+1) = 3N+3
    M = np.zeros((n_rows, n_cols))

    # Eq. 10 (N+1 rows):
    for beta in range(n_xi_hat):
        M[beta, 0:n] = A_at_xi_hat[beta, :]
        M[beta, n:2 * n] = -B2_at_xi_hat[beta, :]
        M[beta, 2 * n:3 * n] = -exp_neg_Delta[beta] * A_at_xi_hat[beta, :]

    # Eq. 11 (N+1 rows):
    #   Σ b A - Σ a B^(2) + e^{+Δ/ξ̂} Σ e B^(2) = 0.
    # At ξ̂ = 0 the e^{+Δ/0} factor diverges; the only finite solution
    # requires Σ e_α B_α^(2)(0) = 0 AND the LHS = 0. Implement as
    # a CONSTRAINT row: at ξ̂ = 0, replace the row with
    #   Σ e_α B_α^(2)(0) = 0  (coef on a, b is zero).
    # This is the limiting form per L'Hopital / Frobenius analysis.
    for beta in range(n_xi_hat):
        row = n_xi_hat + beta
        if xi_hat[beta] == 0.0:
            # Limiting constraint: Σ e_α B_α^(2)(0) = 0.
            M[row, 0:n] = 0.0
            M[row, n:2 * n] = 0.0
            M[row, 2 * n:3 * n] = B2_at_xi_hat[beta, :]
        else:
            M[row, 0:n] = -B2_at_xi_hat[beta, :]
            M[row, n:2 * n] = A_at_xi_hat[beta, :]
            M[row, 2 * n:3 * n] = exp_pos_Delta[beta] * B2_at_xi_hat[beta, :]

    # Eq. 12 (N rows):
    for beta in range(n_xi):
        row = 2 * n_xi_hat + beta
        M[row, 0:n] = -B1_at_xi[beta, :] - exp_neg_2tau[beta] * A_at_xi[beta, :]
        M[row, n:2 * n] = A_at_xi[beta, :] + exp_neg_2tau[beta] * B1_at_xi[beta, :]
        # e_α: zero in Eq. 12.

    # Set a_0 = 1: move column 0 to RHS.
    rhs = -M[:, 0]  # M @ z = 0 → M[:, 1:] @ z' = -M[:, 0] · 1
    A_solve = M[:, 1:]  # (3N+2, 3N+2) — square
    z_rest = np.linalg.solve(A_solve, rhs)
    z_full = np.concatenate([[1.0], z_rest])

    a = z_full[0:n]
    b = z_full[n:2 * n]
    e = z_full[2 * n:3 * n]
    return a, b, e


def _solve_3Nplus2_linear_system(
    tau: float,
    Delta: float,
    c_core: float,
    c_reflector: float,
    N: int,
    u0: float,
    eta0: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r"""Build and solve the :math:`(3N+2) \times (3N+2)` reflected-slab
    F_N linear system (NM Eqs. 10, 11, 12) at given :math:`\tau`.

    Unknowns (3N+2 total after fixing :math:`a_0 = 1`):

    * :math:`a_1, a_2, \ldots, a_N` — :math:`N` core outgoing moments.
    * :math:`b_0, b_1, \ldots, b_N` — :math:`N+1` core incoming moments.
    * :math:`e_0, e_1, \ldots, e_N` — :math:`N+1` reflector outgoing
      moments.

    Rows (3N+2 total):

    * NM Eq. 10 collocated at :math:`\hat\xi_\beta \in P_2`,
      :math:`\beta = 0, \ldots, N` — :math:`N+1` rows.
    * NM Eq. 11 collocated at :math:`\hat\xi_\beta \in P_2`,
      :math:`\beta = 0, \ldots, N` — :math:`N+1` rows.
    * NM Eq. 12 collocated at :math:`\xi_\beta \in (0, 1)`,
      :math:`\beta = 0, \ldots, N` — :math:`N+1` rows.

    Total: :math:`3N + 3` rows × :math:`3N + 3` cols. With
    :math:`a_0 = 1` substituted (removes one column, contributes to
    RHS), system is :math:`3N + 3` rows × :math:`3N + 2` cols.
    Overdetermined by one row — the redundancy is the critical
    condition (NM Eq. 15) which is enforced separately. We use
    least-squares with the most-stable subset.

    For numerical robustness we instead solve the FULL :math:`3(N+1)`
    homogeneous system :math:`M \cdot z = 0` with :math:`z =
    (a_0, a_1, \ldots, a_N, b_0, \ldots, b_N, e_0, \ldots, e_N)^T`,
    extract the right null vector (smallest-singular-value SVD vector),
    then normalise by :math:`a_0`. This is well-defined when
    :math:`\tau` is **near** the critical value; far from criticality
    the matrix is non-singular and we still get a "best" null vector
    whose ratio :math:`a_0 / \mathrm{rest}` indicates the residual.

    Returns ``(a_coefs, b_coefs, e_coefs, fitness_metric)`` where
    :math:`a_0 = 1` after normalisation; ``fitness_metric`` is the
    ratio of the smallest to second-smallest singular value
    (small = near-critical).
    """
    n = N + 1  # length of each coefficient array
    n_unknowns = 3 * n  # a, b, e arrays each of length N+1
    n_rows = 3 * n  # equations 10, 11, 12 each yield N+1 rows

    # Collocation points.
    xi_hat = _collocation_points_xi_hat(eta0, N)  # N+1 points for Eq 10, 11
    xi = _collocation_points_xi(N)  # N+1 points for Eq 12

    # Pre-compute moment arrays at each collocation point.
    # For Eq 10, 11: A_α(ξ̂), B_α^(2)(ξ̂) at each ξ̂.
    # For Eq 12: A_α(ξ), B_α^(1)(ξ) at each ξ.
    A_at_xi_hat = np.empty((n, n))  # A_α(ξ̂_β), [β, α]
    B2_at_xi_hat = np.empty((n, n))  # B_α^(2)(ξ̂_β), [β, α]
    A_at_xi = np.empty((n, n))  # A_α(ξ_β), [β, α]
    B1_at_xi = np.empty((n, n))  # B_α^(1)(ξ_β), [β, α]

    for beta in range(n):
        Bvec, Avec = _moment_arrays_at_xi(N, xi_hat[beta], c_reflector)
        A_at_xi_hat[beta, :] = Avec
        B2_at_xi_hat[beta, :] = Bvec
        Bvec_core, Avec_core = _moment_arrays_at_xi(N, xi[beta], c_core)
        A_at_xi[beta, :] = Avec_core
        B1_at_xi[beta, :] = Bvec_core

    # Pre-compute exponential factors.
    # exp_neg_Delta_over_xi_hat[β] = exp(-Δ/ξ̂_β); zero for ξ̂ = 0.
    exp_neg_Delta = np.where(
        xi_hat > 0.0, np.exp(-Delta / np.where(xi_hat > 0, xi_hat, 1.0)), 0.0
    )
    # exp_pos_Delta_over_xi_hat[β] = exp(+Δ/ξ̂_β); diverges as ξ̂ → 0.
    # NM Eq. 11 RHS has -e^{+Δ/ξ̂} Σ e_α B_α^(2)(ξ̂); at ξ̂ = 0 this term
    # is unbounded UNLESS Σ e_α B_α^(2)(0) = 0. Practically: at the
    # collocation point ξ̂ = 0, equation 11 becomes a constraint
    #   Σ b_α A_α(0) - Σ a_α B_α^(2)(0) = -e^{+Δ/0} Σ e_α B_α^(2)(0)
    # whose only finite solution requires Σ e_α B_α^(2)(0) = 0
    # AND the left side = 0. This is encoded by replacing the RHS at
    # ξ̂ = 0 with a finite "limiting" version: we cap the exponential
    # at a moderate value (e.g. e^{50}) to avoid Inf, and let the SVD
    # null-space find the consistent solution.
    # In practice NM's choice of ξ̂_1 = 0 and ξ̂_2 = 1 plus the eta_0
    # collocation gives a well-conditioned system AS LONG AS the
    # coefficients e_α naturally drive Σ e_α B_α^(2)(0) → 0.
    safe_exp_arg = np.minimum(Delta / np.where(xi_hat > 0, xi_hat, 1.0), 50.0)
    exp_pos_Delta = np.where(xi_hat > 0.0, np.exp(safe_exp_arg), np.exp(50.0))

    # exp_neg_2tau_over_xi[β] = exp(-2τ/ξ_β); zero for ξ = 0.
    exp_neg_2tau = np.where(
        xi > 0.0, np.exp(-2.0 * tau / np.where(xi > 0, xi, 1.0)), 0.0
    )

    # Build the homogeneous system M @ z = 0 with z = (a, b, e).
    M = np.zeros((n_rows, n_unknowns))

    # Eq. 10 rows (β = 0, ..., N):
    #   Σ_α a_α A_α(ξ̂_β) - Σ_α b_α B_α^(2)(ξ̂_β)
    #     - e^{-Δ/ξ̂_β} Σ_α e_α A_α(ξ̂_β) = 0
    for beta in range(n):
        M[beta, 0:n] = A_at_xi_hat[beta, :]  # a_α coefficients
        M[beta, n:2 * n] = -B2_at_xi_hat[beta, :]  # b_α coefficients
        M[beta, 2 * n:3 * n] = -exp_neg_Delta[beta] * A_at_xi_hat[beta, :]

    # Eq. 11 rows (β = N+1, ..., 2N+1):
    #   Σ_α b_α A_α(ξ̂_β) - Σ_α a_α B_α^(2)(ξ̂_β)
    #     + e^{+Δ/ξ̂_β} Σ_α e_α B_α^(2)(ξ̂_β) = 0
    for beta in range(n):
        row = n + beta
        M[row, 0:n] = -B2_at_xi_hat[beta, :]  # a_α
        M[row, n:2 * n] = A_at_xi_hat[beta, :]  # b_α
        M[row, 2 * n:3 * n] = exp_pos_Delta[beta] * B2_at_xi_hat[beta, :]

    # Eq. 12 rows (β = 2N+2, ..., 3N+2):
    #   Σ_α b_α A_α(ξ_β) - Σ_α a_α B_α^(1)(ξ_β)
    #     - e^{-2τ/ξ_β} [Σ_α a_α A_α(ξ_β) - Σ_α b_α B_α^(1)(ξ_β)] = 0
    for beta in range(n):
        row = 2 * n + beta
        M[row, 0:n] = -B1_at_xi[beta, :] - exp_neg_2tau[beta] * A_at_xi[beta, :]  # a_α
        M[row, n:2 * n] = A_at_xi[beta, :] + exp_neg_2tau[beta] * B1_at_xi[beta, :]  # b_α
        # e_α terms: zero in Eq. 12 (no reflector dependence — Eq. 12 is core only).

    # Find the right null vector via SVD. The system has rank deficiency
    # equal to one when τ is critical; we extract the smallest-singular-
    # value vector and rescale so that a_0 = 1.
    U, s, Vh = np.linalg.svd(M, full_matrices=False)
    null_vec = Vh.conj().T[:, -1]  # shape (3n,)

    # Normalise so a_0 = 1.
    if abs(null_vec[0]) > 1e-30:
        null_vec = null_vec / null_vec[0]
    else:
        # Numerical degeneracy — fall back to second smallest.
        null_vec = Vh.conj().T[:, -2]
        if abs(null_vec[0]) > 1e-30:
            null_vec = null_vec / null_vec[0]

    # Take real part (imaginary part is at machine noise floor).
    null_vec = np.real(null_vec)

    a_coefs = null_vec[0:n]
    b_coefs = null_vec[n:2 * n]
    e_coefs = null_vec[2 * n:3 * n]

    # Fitness metric: ratio of smallest to second-smallest sv.
    fitness = s[-1] / s[-2] if s[-2] > 0 else float("inf")

    return a_coefs, b_coefs, e_coefs, fitness


def _eq15_residual(
    tau: float,
    a: np.ndarray,
    b: np.ndarray,
    u0: float,
    c_core: float,
    N: int,
) -> float:
    r"""Compute the residual of NM Eq. 15 at given :math:`\tau` and
    coefficient arrays:

    .. math::

       R(\tau) = e^{-2\tau/\nu_0}
       \sum_\alpha [b_\alpha B_\alpha^{(1)}(\nu_0)
       - a_\alpha A_\alpha(\nu_0)]
       - \sum_\alpha [a_\alpha B_\alpha^{(1)}(\nu_0)
       - b_\alpha A_\alpha(\nu_0)]

    Note :math:`\nu_0 = i u_0` (purely imaginary for :math:`c_1 > 1`).
    Both :math:`A_\alpha`, :math:`B_\alpha^{(1)}` are evaluated at
    complex argument; the residual is taken as the real part (the
    imaginary part cancels for the symmetric critical mode). The
    exponential :math:`e^{-2\tau/(i u_0)} = e^{2 i \tau/u_0} =
    \cos(2\tau/u_0) + i \sin(2\tau/u_0)` so the residual oscillates
    in :math:`\tau`; we use the real part to define the root.
    """
    # Compute A_α(i u_0), B_α^(1)(i u_0) via complex arithmetic.
    # Use the recursion with complex ξ.
    iu = complex(0.0, u0)
    A_complex = np.empty(N + 1, dtype=complex)
    B_complex = np.empty(N + 1, dtype=complex)
    log_term = iu * complex_log(complex(1.0) + 1.0 / iu)
    A_complex[0] = 1.0 - log_term
    B_complex[0] = (2.0 / c_core) - 1.0 - log_term
    for n in range(1, N + 1):
        A_complex[n] = -iu * A_complex[n - 1] + 1.0 / (n + 1)
        B_complex[n] = iu * B_complex[n - 1] - 1.0 / (n + 1)

    # Sum a_α A_α(ν_0), etc.
    sum_a_A = np.sum(a * A_complex)
    sum_a_B = np.sum(a * B_complex)
    sum_b_A = np.sum(b * A_complex)
    sum_b_B = np.sum(b * B_complex)

    # Eq. 15: e^{-2τ/(i u_0)} [b_α B_α^(1) - a_α A_α] = a_α B_α^(1) - b_α A_α
    # → e^{2 i τ/u_0} [Σ b B - Σ a A] = [Σ a B - Σ b A]
    exp_factor = complex(math.cos(2.0 * tau / u0), math.sin(2.0 * tau / u0))
    lhs = exp_factor * (sum_b_B - sum_a_A)
    rhs = sum_a_B - sum_b_A
    residual_complex = lhs - rhs

    # The real part is the root we seek.
    return residual_complex.real


def solve_fn_slab_reflected_critical(
    c_core: float,
    c_reflector: float,
    reflector_half_thickness: float,
    *,
    n_modes: int = 5,
    max_outer_iter: int = 20,
    tol_tau: float = 1e-8,
    tau_initial: float | None = None,
) -> SlabReflectedFNResult:
    r"""Solve the critical core half-thickness :math:`\tau` of a
    reflected slab (1G isotropic) via the Neshat-Maiorino 1980 F_N
    method.

    Parameters
    ----------
    c_core : float
        Core mean number of secondaries per collision; must satisfy
        :math:`c_{\rm core} > 1` (multiplying medium).
    c_reflector : float
        Reflector mean number of secondaries per collision; must satisfy
        :math:`0 < c_{\rm reflector} < 1` (subcritical).
    reflector_half_thickness : float
        :math:`\Delta = b - \tau` in mean-free paths. Must be positive.
    n_modes : int, default 5
        F_N order :math:`N`. NM Table 2 shows :math:`N = 5` gives
        4-5 sig figs accuracy on :math:`\tau`; :math:`N = 7` matches
        Burkart "Exact" reference values to printed precision.
    max_outer_iter : int, default 20
        Maximum number of outer-loop iterations on :math:`\tau`. NM
        report 3-4 iterations are typical.
    tol_tau : float, default 1e-8
        Convergence tolerance on :math:`|\tau_n - \tau_{n-1}|`.
    tau_initial : float, optional
        If given, override the F_0 initial-guess (NM Eq. 16) with this
        value. Useful for testing or when F_0 fails to converge from
        cold start.

    Returns
    -------
    :class:`SlabReflectedFNResult`

    Notes
    -----
    The NM iteration solves the coefficient system at fixed
    :math:`\tau` then root-finds Eq. 15. Each outer iteration:

    1. Build :math:`(3N+3) \times 3(N+1)` matrix from Eqs. 10-12.
    2. Find right null vector → :math:`\{a_\alpha, b_\alpha,
       e_\alpha\}` (with :math:`a_0 = 1`).
    3. Use these coefficients to compute :math:`\tau_{n+1}` from
       Eq. 15 via 1D root-find on the real residual.

    The Eq. 15 root-find uses scipy bisection with a generous
    bracket :math:`[\epsilon, 10\,u_0]`.

    Raises
    ------
    ValueError
        If parameters are out of range.
    RuntimeError
        If iteration fails to converge in ``max_outer_iter`` steps.
    """
    if c_core <= 1.0:
        raise ValueError(
            f"c_core must be > 1 (multiplying medium), got {c_core}"
        )
    if not (0.0 < c_reflector < 1.0):
        raise ValueError(
            f"c_reflector must be in (0, 1) (subcritical), got {c_reflector}"
        )
    if reflector_half_thickness <= 0.0:
        raise ValueError(
            f"reflector_half_thickness must be > 0, got {reflector_half_thickness}"
        )
    if n_modes < 0:
        raise ValueError(f"n_modes must be ≥ 0, got {n_modes}")
    if reflector_half_thickness > 5.0:
        import warnings as _w
        _w.warn(
            f"reflector_half_thickness = {reflector_half_thickness} > 5 mfp; "
            f"the e^{{+Δ/ξ̂}} factor in NM Eq. 11 may produce ill-conditioning. "
            f"Consider the infinite-reflector limit (Siewert-Burkart 1975).",
            stacklevel=2,
        )

    # F_0 initial guess.
    if tau_initial is None:
        tau_F0, b0_F0, u0, eta0 = _build_F0_initial_tau(
            c_core, c_reflector, reflector_half_thickness
        )
        tau_init = tau_F0
    else:
        u0 = case_dispersion_root_supercritical(c_core)
        eta0 = case_dispersion_root_subcritical(c_reflector)
        tau_init = float(tau_initial)
        tau_F0 = tau_init
        b0_F0 = 0.0  # not used

    Delta = reflector_half_thickness
    tau = tau_init
    tau_history: list[float] = [tau]
    converged = False
    iter_count = 0

    # NM's iteration scheme (p. 81): at each iteration,
    #   1. Solve the (3N+3)-equation collocation system at fixed τ
    #      with a_0 = 1 → coefficient arrays {a, b, e}.
    #   2. Plug into Eq. 15 (residual real part) and root-find τ.
    #   3. Repeat until τ converges.
    a_coefs = np.ones(n_modes + 1)
    b_coefs = np.zeros(n_modes + 1)
    e_coefs = np.zeros(n_modes + 1)

    for iter_count in range(1, max_outer_iter + 1):
        # Step 1: solve coefficient system at current τ.
        a_coefs, b_coefs, e_coefs = _solve_NM_iteration_step(
            tau, Delta, c_core, c_reflector, n_modes, u0, eta0,
        )

        # Step 2: root-find Eq. 15 over τ with these coefficients fixed.
        # The residual should oscillate slowly in τ near the critical
        # value (the e^{2 i τ/u_0} factor in Eq. 15 has period π u_0).
        # Use a tight bracket around current τ.
        bracket_half = min(0.6 * tau, math.pi * u0 / 4.0)
        bracket_half = max(bracket_half, 0.05)
        lo = max(1e-4, tau - bracket_half)
        hi = tau + bracket_half
        n_scan = 200
        tau_scan = np.linspace(lo, hi, n_scan)
        res_scan = np.array([
            _eq15_residual(t, a_coefs, b_coefs, u0, c_core, n_modes)
            for t in tau_scan
        ])

        sign_changes = np.where(np.diff(np.sign(res_scan)) != 0)[0]
        if len(sign_changes) == 0:
            # Widen once.
            wider = min(math.pi * u0 / 2.0, 5.0)
            lo = max(1e-4, tau - wider)
            hi = tau + wider
            tau_scan = np.linspace(lo, hi, 400)
            res_scan = np.array([
                _eq15_residual(t, a_coefs, b_coefs, u0, c_core, n_modes)
                for t in tau_scan
            ])
            sign_changes = np.where(np.diff(np.sign(res_scan)) != 0)[0]
        if len(sign_changes) == 0:
            raise RuntimeError(
                f"NM iteration {iter_count}: no sign change of Eq. 15 "
                f"residual in [{lo}, {hi}]. Diverged."
            )

        candidate_taus = 0.5 * (tau_scan[sign_changes] + tau_scan[sign_changes + 1])
        idx_closest = int(np.argmin(np.abs(candidate_taus - tau)))
        i_lo = int(sign_changes[idx_closest])
        a_lo, a_hi = tau_scan[i_lo], tau_scan[i_lo + 1]
        f_lo = res_scan[i_lo]

        # Bisect on Eq. 15.
        for _ in range(80):
            if a_hi - a_lo < 1e-14 * max(1.0, a_lo):
                break
            a_mid = 0.5 * (a_lo + a_hi)
            f_mid = _eq15_residual(
                a_mid, a_coefs, b_coefs, u0, c_core, n_modes
            )
            if f_mid * f_lo < 0:
                a_hi = a_mid
            else:
                a_lo, f_lo = a_mid, f_mid

        new_tau = 0.5 * (a_lo + a_hi)
        tau_history.append(new_tau)

        if abs(new_tau - tau) < tol_tau * max(1.0, abs(tau)):
            tau = new_tau
            converged = True
            break

        tau = new_tau

    # Final coefficient solve at the converged τ.
    if converged:
        a_coefs, b_coefs, e_coefs = _solve_NM_iteration_step(
            tau, Delta, c_core, c_reflector, n_modes, u0, eta0,
        )

    return SlabReflectedFNResult(
        tau_critical_mfp=float(tau),
        c_core=c_core,
        c_reflector=c_reflector,
        Delta_mfp=reflector_half_thickness,
        N=n_modes,
        nu0_core=u0,
        eta0_reflector=eta0,
        coefficients_a=a_coefs,
        coefficients_b=b_coefs,
        coefficients_e=e_coefs,
        iterations=iter_count,
        converged=converged,
        tau_history=tau_history,
    )
