"""Diagnostic: Issue #131 Probe F — 2G 2-region white_f4 (the shipped benchmark).

Created by numerics-investigator on 2026-04-23.

The actual Issue #131 benchmark: `peierls_slab_2eg_2rg` with
boundary="white_f4" on the unified path. Target: after the
P_esc/G_bc closed-form fix for multi-region slab, the rel_diff
should drop from 1.5e-2 to ~1e-8 or better (limited by remaining
volume-kernel adaptive quadrature).
"""
from __future__ import annotations

import numpy as np
import pytest
import time


@pytest.mark.slow
def test_issue131_probe_f_2eg_2rg_white_f4():
    from orpheus.derivations._xs_library import LAYOUTS, get_xs
    from orpheus.derivations.cp_slab import _THICKNESSES
    from orpheus.derivations.peierls_geometry import (
        SLAB_POLAR_1D,
        solve_peierls_mg,
    )
    from orpheus.derivations.peierls_slab import solve_peierls_eigenvalue

    n_regions = 2
    ng_key = "2g"
    layout = LAYOUTS[n_regions]
    thicknesses = _THICKNESSES[n_regions]
    xs_list = [get_xs(region, ng_key) for region in layout]

    sig_t_list = [xs["sig_t"] for xs in xs_list]
    sig_s_list = [xs["sig_s"] for xs in xs_list]
    nu_list = [xs["nu"] * xs["sig_f"] for xs in xs_list]
    chi_list = [xs["chi"] for xs in xs_list]

    n_panels, p_order, dps = 2, 3, 20

    t0 = time.time()
    sol_nat = solve_peierls_eigenvalue(
        sig_t_regions=sig_t_list,
        sig_s_matrices=sig_s_list,
        nu_sig_f_all=nu_list,
        chi_all=chi_list,
        thicknesses=thicknesses,
        n_panels_per_region=n_panels, p_order=p_order,
        precision_digits=dps, boundary="white",
    )
    t_nat = time.time() - t0

    sig_t_mg = np.stack(sig_t_list)
    sig_s_mg = np.stack(sig_s_list)
    nu_sig_f_mg = np.stack(nu_list)
    chi_mg = np.stack(chi_list)
    radii_mg = np.cumsum(np.asarray(thicknesses, dtype=float))

    t0 = time.time()
    sol_mg = solve_peierls_mg(
        SLAB_POLAR_1D, radii_mg,
        sig_t=sig_t_mg, sig_s=sig_s_mg,
        nu_sig_f=nu_sig_f_mg, chi=chi_mg,
        boundary="white_f4",
        n_panels_per_region=n_panels, p_order=p_order,
        dps=dps, tol=1e-12,
    )
    t_mg = time.time() - t0

    rel = abs(sol_mg.k_eff - sol_nat.k_eff) / abs(sol_nat.k_eff)
    print(f"\nProbe F (2G 2-region white_f4):")
    print(f"  native  k_eff = {sol_nat.k_eff:.12f}  ({t_nat:.1f}s)")
    print(f"  unified k_eff = {sol_mg.k_eff:.12f}  ({t_mg:.1f}s)")
    print(f"  rel_diff      = {rel:.3e}")


if __name__ == "__main__":
    test_issue131_probe_f_2eg_2rg_white_f4()
