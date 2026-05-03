r"""Atalay 1997 Eq 54 sphere critical-radius solver — linearly anisotropic.

The sphere is treated as the **odd-mode** of the slab problem on
:math:`[-R, R]` via the antisymmetric boundary condition
:math:`\psi(x, \mu) = -\psi(-x, -\mu)` (Atalay Eq 47). The structural
changes from slab to sphere reduce to:

* :math:`T(R, \mu) \to T_1(R, \mu)` (sign flip in second exponential).
* :math:`K_j \to L_j` (same kernel structure, T → T_1).
* sin↔cos shuffle in the LHS criticality argument:

.. math::

   \arctan\!\frac{\sin[(d+z_0)/|\nu_0|] - R \sin[(d-z_0)/|\nu_0|]}
                {\cos[(d+z_0)/|\nu_0|] + R \cos[(d-z_0)/|\nu_0|]}
   = \arctan\!\frac{\big(L_0 \bar\nu - 3 f_1 (1-c) \bar\nu
                     [L_1 \bar\nu d(\nu_0^2) - L_0 \nu_0^2 d(\bar\nu^2)]\big) |\nu_0|}
                  {(1+L_2) d(\nu_0 \bar\nu) d(-\nu_0 \bar\nu)
                   + L_1 \bar\nu d(\nu_0^2) - L_0 \nu_0^2 d(\bar\nu^2)}

where :math:`d` here is the sphere **radius** (interpretation: the
half-thickness of the equivalent slab equals the sphere radius).

References
----------

* Atalay 1997, *Prog. Nucl. Energy* **31**(3), 229-252.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..core.dispersion import case_atalay_u0, nu_bar_atalay
from ..core.extrapolated_endpoint import atalay_z0
from ..core.half_range import atalay_L_moments


@dataclass(frozen=True)
class CaseMethodSphereResult:
    r"""Result of :func:`solve_case_method_sphere_critical`.

    Attributes
    ----------
    R_critical_mfp : float
        Critical sphere radius in mean free paths.
    c, R_refl, f1 : float
        Input material/BC parameters (``R_refl`` is the reflection
        coefficient at the outer surface, distinct from the sphere
        radius which is what we solve for).
    u0 : float
        :math:`u_0 = |\nu_0|`.
    nu_bar : float
        :math:`\bar\nu`.
    z0 : float
        Extrapolated endpoint.
    L_moments : dict[int, float]
        :math:`L_0, L_1, L_2`.
    eq54_residual : float
        Eq 54 residual at the converged radius.
    n_bracket_points : int
        Bracket-scan size.
    n_bisect_iters : int
        Bisection iterations.
    """

    R_critical_mfp: float
    c: float
    R_refl: float
    f1: float
    u0: float
    nu_bar: float
    z0: float
    L_moments: dict[int, float]
    eq54_residual: float
    n_bracket_points: int
    n_bisect_iters: int


def _eq54_residual(
    radius: float,
    *,
    c: float,
    R_refl: float,
    f1: float,
    u0: float,
    nu_bar: float,
    z0: float,
) -> float:
    r"""Atalay Eq 54 residual: :math:`\sin(\theta_{LHS} - \theta_{RHS})`.

    Smooth and crosses zero at the critical radius.
    """
    Ls = atalay_L_moments(
        c=c, f1=f1, u0=u0, nu_bar=nu_bar, R=R_refl, d_thick=radius, j_max=2,
    )
    L0, L1, L2 = Ls[0], Ls[1], Ls[2]

    nu0_sq = -u0**2
    nu_bar_sq = nu_bar**2

    def d_func(ab):
        return 1.0 + 3.0 * f1 * (1.0 - c) * ab

    d_nu0_nubar_mod_sq = 1.0 + (3.0 * f1 * (1.0 - c) * u0 * nu_bar)**2
    L_factor = (1.0 + L2) * d_nu0_nubar_mod_sq

    d_nu0_sq = d_func(nu0_sq)
    d_nubar_sq = d_func(nu_bar_sq)

    eq54_rhs_den = L_factor + L1 * nu_bar * d_nu0_sq - L0 * nu0_sq * d_nubar_sq

    inner_bracket = L1 * nu_bar * d_nu0_sq - L0 * nu0_sq * d_nubar_sq
    eq54_rhs_num = (L0 * nu_bar - 3.0 * f1 * (1.0 - c) * nu_bar * inner_bracket) * u0

    theta_rhs = math.atan2(eq54_rhs_num, eq54_rhs_den)

    # Eq 54 LHS — sphere parity flip relative to slab Eq 46.
    # Numerator: sin[(d+z_0)/u_0] - R sin[(d-z_0)/u_0]
    # Denominator: cos[(d+z_0)/u_0] + R cos[(d-z_0)/u_0]
    a1 = (radius - z0) / u0   # (d - z_0)/u_0
    a2 = (radius + z0) / u0   # (d + z_0)/u_0
    lhs_num = math.sin(a2) - R_refl * math.sin(a1)
    lhs_den = math.cos(a2) + R_refl * math.cos(a1)
    theta_lhs = math.atan2(lhs_num, lhs_den)

    # Residual: tan(θ_LHS) - tan(θ_RHS), via the cross-product form
    #   sin(θ_LHS) cos(θ_RHS) - cos(θ_LHS) sin(θ_RHS) = sin(θ_LHS - θ_RHS).
    # No modulo wrapping — that introduces spurious zero crossings. The true
    # criticality has θ_LHS = θ_RHS, so sin(θ_LHS - θ_RHS) = 0.
    return math.sin(theta_lhs - theta_rhs)


def solve_case_method_sphere_critical(
    c: float,
    R_refl: float = 0.0,
    f1: float = 0.0,
    *,
    radius_min: float = 0.10,
    radius_max: float = 80.0,
    n_bracket: int = 200,
    bisect_tol: float = 1e-10,
    max_bisect: int = 80,
    mode: int = 1,
) -> CaseMethodSphereResult:
    r"""Solve Atalay 1997 sphere critical radius via the Case singular-
    eigenfunction method (linearly anisotropic, odd-mode parity-flip
    of the slab problem).

    Parameters
    ----------
    c : float
        Mean number of secondaries per collision (``c > 1``).
    R_refl : float, default 0.0
        Reflection coefficient at the outer sphere surface.
        ``R_refl = 0`` is bare; ``R_refl ∈ (0, 1)`` is partial reflection;
        ``R_refl = 1`` is rejected (degenerate case).
    f1 : float, default 0.0
        Linearly-anisotropic mean cosine.
    radius_min, radius_max : float, default 0.10, 80.0
        Bracket for the radius scan in mfp.
    n_bracket : int, default 200
        Bracket-scan resolution.
    bisect_tol : float, default 1e-10
        Tolerance on radius.
    max_bisect : int, default 80
        Max bisection iterations.
    mode : int, default 1
        Eigenvalue mode (1 = fundamental).

    Returns
    -------
    :class:`CaseMethodSphereResult`
    """
    if c <= 1.0:
        raise ValueError(f"solve_case_method_sphere_critical requires c > 1, got c={c}")
    if R_refl >= 1.0:
        raise ValueError(
            f"R_refl = {R_refl} ≥ 1 makes the sphere radius drop out "
            f"of Eq 54 (perfect reflection at the surface)."
        )
    if R_refl < 0:
        raise ValueError(f"R_refl must be in [0, 1), got {R_refl}")
    if mode < 1:
        raise ValueError(f"mode must be ≥ 1, got {mode}")

    disp = case_atalay_u0(c=c, f1=f1)
    u0 = disp.u0
    nu_bar = nu_bar_atalay(c=c, f1=f1, u0=u0)
    z0 = atalay_z0(c=c, f1=f1, u0=u0, nu_bar=nu_bar)

    radius_scan = np.linspace(radius_min, radius_max, n_bracket)
    residuals = np.empty(n_bracket)
    for i, r in enumerate(radius_scan):
        try:
            residuals[i] = _eq54_residual(
                r, c=c, R_refl=R_refl, f1=f1, u0=u0, nu_bar=nu_bar, z0=z0,
            )
        except Exception:
            residuals[i] = float("nan")

    valid = ~np.isnan(residuals)
    valid_idx = np.where(valid)[0]
    if len(valid_idx) < 2:
        raise RuntimeError(
            f"Eq 54 residual evaluation failed across the entire bracket "
            f"[{radius_min}, {radius_max}] for (c={c}, R_refl={R_refl}, f_1={f1})."
        )
    sgn = np.sign(residuals[valid_idx])
    sign_change_indices = np.where(np.diff(sgn) != 0)[0]
    if len(sign_change_indices) == 0:
        raise RuntimeError(
            f"No sign change of Eq 54 residual in [{radius_min}, {radius_max}] "
            f"for (c={c}, R_refl={R_refl}, f_1={f1})."
        )

    if mode > len(sign_change_indices):
        raise RuntimeError(
            f"Asked for mode={mode} but only {len(sign_change_indices)} sign "
            f"changes in [{radius_min}, {radius_max}]."
        )

    idx_in_valid = sign_change_indices[mode - 1]
    i_lo, i_hi = valid_idx[idx_in_valid], valid_idx[idx_in_valid + 1]
    r_lo, r_hi = radius_scan[i_lo], radius_scan[i_hi]
    f_lo = _eq54_residual(r_lo, c=c, R_refl=R_refl, f1=f1, u0=u0, nu_bar=nu_bar, z0=z0)
    f_hi = _eq54_residual(r_hi, c=c, R_refl=R_refl, f1=f1, u0=u0, nu_bar=nu_bar, z0=z0)

    iters = 0
    while iters < max_bisect:
        if r_hi - r_lo < bisect_tol * max(1.0, r_lo):
            break
        r_mid = 0.5 * (r_lo + r_hi)
        f_mid = _eq54_residual(r_mid, c=c, R_refl=R_refl, f1=f1, u0=u0, nu_bar=nu_bar, z0=z0)
        if f_mid * f_lo < 0:
            r_hi, f_hi = r_mid, f_mid
        else:
            r_lo, f_lo = r_mid, f_mid
        iters += 1

    r_crit = 0.5 * (r_lo + r_hi)
    Ls_final = atalay_L_moments(
        c=c, f1=f1, u0=u0, nu_bar=nu_bar, R=R_refl, d_thick=r_crit,
    )
    final_resid = _eq54_residual(
        r_crit, c=c, R_refl=R_refl, f1=f1, u0=u0, nu_bar=nu_bar, z0=z0,
    )

    return CaseMethodSphereResult(
        R_critical_mfp=float(r_crit),
        c=float(c), R_refl=float(R_refl), f1=float(f1),
        u0=float(u0), nu_bar=float(nu_bar), z0=float(z0),
        L_moments=Ls_final,
        eq54_residual=float(final_resid),
        n_bracket_points=n_bracket,
        n_bisect_iters=iters,
    )
