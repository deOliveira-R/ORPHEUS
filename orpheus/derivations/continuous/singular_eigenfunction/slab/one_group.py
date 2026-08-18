r"""Atalay 1997 Eq 46 reflected-slab critical-thickness solver.

For a homogeneous slab of half-thickness :math:`d` (mfp), with linearly
anisotropic scattering (mean cosine :math:`f_1`) and specular-style
reflection at both ends with reflection coefficient :math:`R`, the
**even-mode** criticality condition is (Atalay Eq 46):

.. math::

   \tan\!\Big(\pm\frac{\pi}{2} - \theta_{LHS}\Big) = \theta_{RHS}

where

.. math::

   \theta_{LHS} = \arctan\!\frac{R \sin[(d-z_0)/|\nu_0|] + \sin[(d+z_0)/|\nu_0|]}
                                  {R \cos[(d-z_0)/|\nu_0|] - \cos[(d+z_0)/|\nu_0|]}

.. math::

   \theta_{RHS} = \arctan\!\frac{\big(K_0 \bar\nu - 3 f_1 (1-c) \bar\nu
                                  [K_1 \bar\nu d(\nu_0^2) - K_0 \nu_0^2 d(\bar\nu^2)]\big) |\nu_0|}
                                  {K + K_1 \bar\nu d(\nu_0^2) - K_0 \nu_0^2 d(\bar\nu^2)}

with :math:`K = (1+K_2) d(\nu_0 \bar\nu) d(-\nu_0 \bar\nu)`,
:math:`\nu_0^2 = -u_0^2`, and :math:`K_j(c, R, d)` from
:func:`...core.half_range.atalay_K_moments`.

The :math:`\pm \pi/2` ambiguity in the LHS reflects the multiple
critical thicknesses (the second and third eigenvalues — first three
modes) that Atalay reports per :math:`(c, R, f_1)` cell in Tables 2-5.
We bracket-scan and return the **fundamental mode** (smallest
:math:`d`).

Validity
--------

Atalay Eq 5: :math:`c \le 1 + 1/(3 f_1)`. Outside this range complex
eigenvalues appear and first-order Fredholm iteration breaks.

Precision floor at small slab thicknesses (ERR-038)
---------------------------------------------------

Atalay Eq 46 is itself a **first-order approximation** to the full
Fredholm equation Eq 32 — Atalay (p.236) explicitly omits the integral
term in Eq 32 to derive Eq 46. Atalay (p.246) further notes "we expect
some improvement in the accuracy especially for the small slab
thicknesses" if higher-order iteration were performed (which the
paper does NOT do). Empirically, our solver matches Atalay's Tables
to:

* Machine precision (1e-5 relative) at :math:`2d \ge 20` mfp.
* 4e-4 relative at :math:`2d = 2` mfp.
* 3.2% at :math:`2d = 0.20` mfp.
* 5% at :math:`2d \approx 0.015` mfp (Table 2 R=0.99 column).

This 1/d-scaling error is the **published reference's first-order
approximation precision floor** at small slab thicknesses — NOT a bug
in our implementation of Eq 46. See ERR-038 in
``docs/theory/verification/error_catalog.rst`` for the cascade
evidence and structurally-independent grounding (Atalay's own text,
the 1/d_crit scaling fingerprint, and self-consistency at moderate d).

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
from ..core.half_range import atalay_K_moments


@dataclass(frozen=True)
class CaseMethodSlabResult:
    r"""Result of :func:`solve_case_method_slab_critical`.

    Attributes
    ----------
    d_critical_mfp : float
        Critical half-thickness in mean free paths. Slab full
        thickness is :math:`2 d`.
    c, R, f1 : float
        Input material/BC parameters.
    u0 : float
        Magnitude of the discrete Case eigenvalue
        :math:`\nu_0 = i u_0`.
    nu_bar : float
        :math:`\bar\nu = \gamma^{(1)}/\gamma^{(0)}` (Atalay Eq 24).
    z0 : float
        Extrapolated endpoint (Atalay Eq 42).
    K_moments : dict[int, float]
        :math:`K_0, K_1, K_2` evaluated at the converged
        :math:`d_{\rm critical}`.
    eq46_residual : float
        :math:`\tan(\theta_{LHS}^{eq46}) - \theta_{RHS}^{eq46}` at the
        converged :math:`d`. Should be small.
    n_bracket_points : int
        Number of points used in the initial bracket scan.
    n_bisect_iters : int
        Bisection iterations.
    converged : bool
        Whether the bisection met ``bisect_tol`` — **required, no
        default** (#340).

        The bracket width is the witness, NOT ``eq46_residual``: the
        residual's documented spec is *"Should be small"* with no
        threshold stated anywhere, so deriving convergence from it
        would mean inventing a tolerance — the exact "assert a number
        nobody measured" move #340 exists to remove.  ``bisect_tol``
        is a real, caller-supplied parameter, so the honest predicate
        is the loop's own: ``d_hi - d_lo < bisect_tol * max(1, d_lo)``,
        re-evaluated on the final bracket.

        A field whose safe value is the optimistic one lies by
        omission, so this one has no default and every producer must
        state it.
    """

    d_critical_mfp: float
    c: float
    R: float
    f1: float
    u0: float
    nu_bar: float
    z0: float
    K_moments: dict[int, float]
    eq46_residual: float
    n_bracket_points: int
    n_bisect_iters: int
    converged: bool


def _eq46_residual(
    d_thick: float,
    *,
    c: float,
    R: float,
    f1: float,
    u0: float,
    nu_bar: float,
    z0: float,
    dps: int,
    maxdegree: int,
    branch_sign: int,
) -> float:
    r"""Atalay Eq 46 residual: :math:`\tan(\theta_{LHS}) - \tan(\theta_{RHS})`,
    expressed in equivalent angle-difference form to avoid tan-singularity.

    Returns :math:`\sin(\theta_{LHS} - \theta_{RHS})`, which crosses zero
    at the critical thickness without tan blowing up.
    """
    # K_j(c, R, d) for j = 0, 1, 2.
    Ks = atalay_K_moments(
        c=c, f1=f1, u0=u0, nu_bar=nu_bar, R=R, d_thick=d_thick,
        j_max=2, dps=dps, maxdegree=maxdegree,
    )
    K0, K1, K2 = Ks[0], Ks[1], Ks[2]

    # ν_0² = -u_0², ν̄² = nu_bar².
    nu0_sq = -u0**2  # ν_0 = i u_0 ⇒ ν_0² = -u_0²
    nu_bar_sq = nu_bar**2

    # d(ab) = 1 + 3 f_1 (1-c) ab. Use the actual values.
    def d_func(ab):
        return 1.0 + 3.0 * f1 * (1.0 - c) * ab

    # In Eq 46 denominator: K + K_1 ν̄ d(ν_0²) - K_0 ν_0² d(ν̄²).
    # K = (1 + K_2) d(ν_0 ν̄) d(-ν_0 ν̄). For ν_0 = i u_0 (purely imag):
    #   ν_0 ν̄ = i u_0 ν̄ (purely imag).
    # d(ν_0 ν̄) = 1 + 3 f_1 (1-c) (i u_0 ν̄)  — complex.
    # d(-ν_0 ν̄) = 1 - 3 f_1 (1-c) (i u_0 ν̄) — its conjugate.
    # Their product = |1 + 3f_1(1-c)(iu_0 ν̄)|^2 = 1 + [3f_1(1-c) u_0 ν̄]^2.
    # (real positive)
    d_nu0_nubar_mod_sq = 1.0 + (3.0 * f1 * (1.0 - c) * u0 * nu_bar)**2
    K_factor = (1.0 + K2) * d_nu0_nubar_mod_sq

    # d(ν_0²) = d(-u_0²) = 1 + 3 f_1 (1-c) (-u_0²) = 1 - 3 f_1 (1-c) u_0² (real).
    # d(ν̄²) = 1 + 3 f_1 (1-c) ν̄² (real).
    d_nu0_sq = d_func(nu0_sq)
    d_nubar_sq = d_func(nu_bar_sq)

    # Eq 46 RHS denominator (real): K + K_1 ν̄ d(ν_0²) - K_0 ν_0² d(ν̄²).
    eq46_rhs_den = K_factor + K1 * nu_bar * d_nu0_sq - K0 * nu0_sq * d_nubar_sq

    # Eq 46 RHS numerator (real, gets multiplied by |ν_0|=u_0):
    #   [K_0 ν̄ - 3 f_1 (1-c) ν̄ {K_1 ν̄ d(ν_0²) - K_0 ν_0² d(ν̄²)}] · u_0
    inner_bracket = K1 * nu_bar * d_nu0_sq - K0 * nu0_sq * d_nubar_sq
    eq46_rhs_num = (K0 * nu_bar - 3.0 * f1 * (1.0 - c) * nu_bar * inner_bracket) * u0

    # θ_RHS = atan2(num, den). atan2 gives angle ∈ (-π, π].
    theta_rhs = math.atan2(eq46_rhs_num, eq46_rhs_den)

    # Eq 46 LHS:
    #   ±π/2 - atan{(R sin a1 + sin a2) / (R cos a1 - cos a2)}
    # where a1 = (d - z_0)/u_0, a2 = (d + z_0)/u_0.
    a1 = (d_thick - z0) / u0
    a2 = (d_thick + z0) / u0
    # Use atan2 with explicit signs to handle quadrant.
    lhs_num = R * math.sin(a1) + math.sin(a2)
    lhs_den = R * math.cos(a1) - math.cos(a2)
    theta_lhs_inner = math.atan2(lhs_num, lhs_den)

    theta_lhs = branch_sign * (math.pi / 2) - theta_lhs_inner

    # Residual = sin(theta_lhs - theta_rhs) (smooth, zero at criticality).
    diff = theta_lhs - theta_rhs
    # Wrap to nearest π (the criticality is mod π, since tan is π-periodic).
    diff_wrapped = (diff + math.pi / 2) % math.pi - math.pi / 2
    return math.sin(diff_wrapped)


def solve_case_method_slab_critical(
    c: float,
    R: float,
    f1: float = 0.0,
    *,
    d_min: float = 0.05,
    d_max: float = 60.0,
    n_bracket: int = 200,
    bisect_tol: float = 1e-10,
    max_bisect: int = 80,
    dps: int = 25,
    maxdegree: int = 8,
    mode: int = 1,
) -> CaseMethodSlabResult:
    r"""Solve Atalay 1997 reflected-slab critical half-thickness :math:`d`
    via the Case singular-eigenfunction method (linearly anisotropic).

    Parameters
    ----------
    c : float
        Mean number of secondaries per collision (``c > 1``).
    R : float
        Reflection coefficient at both slab faces.
        ``R = 0`` gives the bare slab (Sood ``Ua-1-0-SL`` family);
        ``R ∈ (0, 1)`` is partial reflection;
        ``R = 1`` (perfect) is rejected — slab thickness drops out
        of the criticality condition.
    f1 : float, default 0.0
        Linearly-anisotropic mean cosine. Must satisfy
        :math:`c \le 1 + 1/(3 f_1)` (Atalay Eq 5).
    d_min, d_max : float, default 0.05, 60.0
        Bracket for the half-thickness scan in mfp. Must contain at
        least one criticality.
    n_bracket : int, default 200
        Number of bracket-scan points.
    bisect_tol : float, default 1e-10
        Tolerance on :math:`d`.
    max_bisect : int, default 80
        Max bisection iterations.
    dps : int, default 25
        mpmath decimal precision.
    maxdegree : int, default 8
        mpmath quadrature max order for K_j integrals.
    mode : int, default 1
        Which mode to return (1 = fundamental, 2 = first overtone, etc).
        Atalay Tables 2-5 list up to 3 modes per cell.

    Returns
    -------
    :class:`CaseMethodSlabResult`
    """
    if c <= 1.0:
        raise ValueError(f"solve_case_method_slab_critical requires c > 1, got c={c}")
    if R >= 1.0:
        raise ValueError(
            f"R = {R} ≥ 1 makes the slab thickness drop out of Eq 46 "
            f"(perfect reflection). Atalay omits the R = 1 column from "
            f"Tables 2-5 for this reason."
        )
    if R < 0:
        raise ValueError(f"R must be in [0, 1), got R={R}")
    if mode < 1:
        raise ValueError(f"mode must be ≥ 1, got {mode}")

    # 1. Discrete eigenvalue u_0.
    disp = case_atalay_u0(c=c, f1=f1)
    u0 = disp.u0

    # 2. Moment ratio ν̄.
    nu_bar = nu_bar_atalay(c=c, f1=f1, u0=u0)

    # 3. Extrapolated endpoint z_0.
    z0 = atalay_z0(c=c, f1=f1, u0=u0, nu_bar=nu_bar, dps=dps, maxdegree=maxdegree)

    # 4. Bracket scan on d_thick using both branch signs (±π/2 in LHS).
    # The fundamental mode usually appears for branch_sign = +1.
    d_scan = np.linspace(d_min, d_max, n_bracket)
    crossings: list[tuple[float, float, int]] = []  # (d_lo, d_hi, branch_sign)

    for branch_sign in (+1, -1):
        residuals = np.empty(n_bracket)
        for i, d in enumerate(d_scan):
            try:
                residuals[i] = _eq46_residual(
                    d, c=c, R=R, f1=f1, u0=u0, nu_bar=nu_bar, z0=z0,
                    dps=dps, maxdegree=maxdegree, branch_sign=branch_sign,
                )
            except Exception:
                residuals[i] = float("nan")

        valid = ~np.isnan(residuals)
        valid_idx = np.where(valid)[0]
        if len(valid_idx) < 2:
            continue
        sgn = np.sign(residuals[valid_idx])
        diffs = np.diff(sgn)
        sign_change_indices = np.where(diffs != 0)[0]
        for idx in sign_change_indices:
            i_lo = valid_idx[idx]
            i_hi = valid_idx[idx + 1]
            crossings.append((d_scan[i_lo], d_scan[i_hi], branch_sign))

    if not crossings:
        raise RuntimeError(
            f"No sign change of Eq 46 residual in d ∈ [{d_min}, {d_max}] "
            f"with c={c}, R={R}, f_1={f1}. Try widening the bracket or "
            f"increasing n_bracket."
        )

    # Sort crossings by d_lo (smallest first = fundamental mode).
    crossings.sort(key=lambda x: x[0])

    if mode > len(crossings):
        raise RuntimeError(
            f"Asked for mode={mode} but only {len(crossings)} sign changes "
            f"found in [{d_min}, {d_max}] for (c={c}, R={R}, f_1={f1})."
        )

    d_lo, d_hi, branch_sign = crossings[mode - 1]

    # 5. Bisection.
    f_lo = _eq46_residual(
        d_lo, c=c, R=R, f1=f1, u0=u0, nu_bar=nu_bar, z0=z0,
        dps=dps, maxdegree=maxdegree, branch_sign=branch_sign,
    )
    f_hi = _eq46_residual(
        d_hi, c=c, R=R, f1=f1, u0=u0, nu_bar=nu_bar, z0=z0,
        dps=dps, maxdegree=maxdegree, branch_sign=branch_sign,
    )
    iters = 0
    while iters < max_bisect:
        if d_hi - d_lo < bisect_tol * max(1.0, d_lo):
            break
        d_mid = 0.5 * (d_lo + d_hi)
        f_mid = _eq46_residual(
            d_mid, c=c, R=R, f1=f1, u0=u0, nu_bar=nu_bar, z0=z0,
            dps=dps, maxdegree=maxdegree, branch_sign=branch_sign,
        )
        if f_mid * f_lo < 0:
            d_hi, f_hi = d_mid, f_mid
        else:
            d_lo, f_lo = d_mid, f_mid
        iters += 1

    d_crit = 0.5 * (d_lo + d_hi)

    # 6. Final K moments at converged d for diagnostic.
    Ks_final = atalay_K_moments(
        c=c, f1=f1, u0=u0, nu_bar=nu_bar, R=R, d_thick=d_crit,
        j_max=2, dps=dps, maxdegree=maxdegree,
    )
    final_residual = _eq46_residual(
        d_crit, c=c, R=R, f1=f1, u0=u0, nu_bar=nu_bar, z0=z0,
        dps=dps, maxdegree=maxdegree, branch_sign=branch_sign,
    )

    # The loop's OWN stop predicate, re-evaluated on the final bracket —
    # not a transcription of which way control flow left the loop, and not
    # `iters < max_bisect` (which misreports a bisection that tightens
    # exactly on its last allowed step).  Same discipline as
    # ``orpheus.numerics.convergence.IterationRecord.converged``: derive the
    # fact from the state, so there is nothing to transcribe wrongly.
    # (Named ``orpheus.sn.solver._claims_convergence`` until 2026-08-09; that
    # predicate retired when the drivers began returning records.)
    converged = bool(d_hi - d_lo < bisect_tol * max(1.0, d_lo))

    return CaseMethodSlabResult(
        d_critical_mfp=float(d_crit),
        c=float(c), R=float(R), f1=float(f1),
        u0=float(u0), nu_bar=float(nu_bar), z0=float(z0),
        K_moments=Ks_final,
        eq46_residual=float(final_residual),
        n_bracket_points=n_bracket,
        n_bisect_iters=iters,
        converged=converged,
    )
