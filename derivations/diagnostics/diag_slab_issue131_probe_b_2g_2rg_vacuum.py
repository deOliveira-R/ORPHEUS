"""Diagnostic: Issue #131 Probe B — 2G 2-region vacuum isolation.

Created by numerics-investigator on 2026-04-23.

Probe A (1G 2-region vacuum) passed at 1e-13 — multi-region volume
kernel is fine. This Probe B runs the SHIPPED 2G 2-region fixture with
vacuum BC on both paths, isolating the MG × multi-region interaction
from the F.4 closure.

Expected outcome
----------------
- If ~1e-8 agreement: F.4 closure is the source of the ~1.5% gap.
- If ~1% or larger: MG × multi-region interaction itself has a bug
  (assembly, χ indexing, sig_s transpose on multi-region boundary).
"""
from __future__ import annotations

import numpy as np
import pytest


@pytest.mark.slow
def test_issue131_probe_b_2g_2rg_vacuum():
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

    # Native path — VACUUM
    sol_nat = solve_peierls_eigenvalue(
        sig_t_regions=sig_t_list,
        sig_s_matrices=sig_s_list,
        nu_sig_f_all=nu_list,
        chi_all=chi_list,
        thicknesses=thicknesses,
        n_panels_per_region=2,
        p_order=3,
        precision_digits=20,
        boundary="vacuum",
    )

    # Unified MG path — VACUUM
    sig_t_mg = np.stack(sig_t_list)
    sig_s_mg = np.stack(sig_s_list)
    nu_sig_f_mg = np.stack(nu_list)
    chi_mg = np.stack(chi_list)
    radii_mg = np.cumsum(np.asarray(thicknesses, dtype=float))

    sol_mg = solve_peierls_mg(
        SLAB_POLAR_1D,
        radii_mg,
        sig_t=sig_t_mg,
        sig_s=sig_s_mg,
        nu_sig_f=nu_sig_f_mg,
        chi=chi_mg,
        boundary="vacuum",
        n_panels_per_region=2,
        p_order=3,
        dps=20,
        tol=1e-12,
    )

    rel = abs(sol_mg.k_eff - sol_nat.k_eff) / abs(sol_nat.k_eff)
    print(f"\nProbe B (2G 2-region VACUUM):")
    print(f"  native  k_eff = {sol_nat.k_eff:.12f}")
    print(f"  unified k_eff = {sol_mg.k_eff:.12f}")
    print(f"  rel_diff      = {rel:.3e}")

    assert rel < 1e-6, (
        f"Probe B 2G 2-region vacuum: rel_diff={rel:.3e} "
        f"(native={sol_nat.k_eff:.12f}, unified={sol_mg.k_eff:.12f})"
    )


if __name__ == "__main__":
    test_issue131_probe_b_2g_2rg_vacuum()
