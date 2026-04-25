"""Verify the sphere 1G/2R Hébert +10 % overshoot is structural, not BASE noise.

Per Issue #132 partial-fix evidence, Hébert closure on sphere 1G/2R
gives k_eff = 0.715 vs cp k_inf 0.648 (+10.33 %) at BASE quadrature.
If this is BASE-preset Issue #114 ρ-quadrature noise, RICH should
drop the residual significantly. If it stays at +10 %, the residual
is structural (the Mark uniformity assumption amplified by the
geometric series).

Also tests at higher P_ss n_quad to rule out P_ss numerical
underconvergence as the source.
"""

from __future__ import annotations

import time

import numpy as np

from orpheus.derivations import cp_sphere
from orpheus.derivations._xs_library import LAYOUTS, get_xs
from orpheus.derivations.peierls_geometry import SPHERE_1D, solve_peierls_mg


_QUAD_BASE = dict(
    n_panels_per_region=2, p_order=3,
    n_angular=24, n_rho=24, n_surf_quad=24, dps=15,
)
_QUAD_RICH = dict(
    n_panels_per_region=4, p_order=5,
    n_angular=64, n_rho=48, n_surf_quad=64, dps=20,
)


def _solve(ng_key, n_regions, *, quad):
    layout = LAYOUTS[n_regions]
    xs_list = [get_xs(r, ng_key) for r in layout]
    sig_t = np.vstack([xs["sig_t"] for xs in xs_list])
    sig_s = np.stack([xs["sig_s"] for xs in xs_list], axis=0)
    nu_sig_f = np.vstack([xs["nu"] * xs["sig_f"] for xs in xs_list])
    chi = np.vstack([xs["chi"] for xs in xs_list])
    radii = np.array(cp_sphere._RADII[n_regions])
    sol = solve_peierls_mg(
        SPHERE_1D, radii=radii, sig_t=sig_t, sig_s=sig_s,
        nu_sig_f=nu_sig_f, chi=chi,
        boundary="white_hebert", n_bc_modes=1,
        **quad,
    )
    return sol.k_eff


def main():
    print("=" * 76)
    print("Hébert closure: BASE vs RICH stability test (Issue #132 1G/2R)")
    print("=" * 76)
    cases = [("1g", 1), ("1g", 2), ("2g", 1), ("2g", 2)]
    print(f"\n  {'Config':<14} {'k_inf':>10} {'k_BASE':>14} {'err_BASE':>10} "
          f"{'k_RICH':>14} {'err_RICH':>10} {'Δ BR':>8}")
    print("  " + "-" * 92)
    for ng_key, n_regions in cases:
        kinf = cp_sphere._build_case(ng_key, n_regions).k_inf
        t0 = time.time()
        k_base = _solve(ng_key, n_regions, quad=_QUAD_BASE)
        t_base = time.time() - t0
        err_base = (k_base - kinf) / kinf * 100
        t0 = time.time()
        k_rich = _solve(ng_key, n_regions, quad=_QUAD_RICH)
        t_rich = time.time() - t0
        err_rich = (k_rich - kinf) / kinf * 100
        delta_br = abs(k_rich - k_base) / max(abs(k_rich), 1e-30) * 100
        print(
            f"  {ng_key} {n_regions}r          {kinf:>10.6f} "
            f"{k_base:>14.10f} {err_base:>+9.3f}% "
            f"{k_rich:>14.10f} {err_rich:>+9.3f}% {delta_br:>7.3f}%   "
            f"({t_base:.1f}s|{t_rich:.1f}s)"
        )

    print()
    print("=" * 76)
    print("Interpretation")
    print("=" * 76)
    print()
    print("  If sph 1G/2R RICH err ≈ 0 %: BASE quadrature noise was the cause")
    print("                                → 1G/2R 'limitation' is fictitious,")
    print("                                  Hébert is exact at high quadrature")
    print("  If sph 1G/2R RICH err ≈ +10 %: structural Mark limit confirmed")
    print("                                 → fix needs Davison / augmented Nyström")


if __name__ == "__main__":
    main()
