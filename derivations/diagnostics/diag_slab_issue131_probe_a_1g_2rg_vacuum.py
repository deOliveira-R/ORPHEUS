"""Diagnostic: Issue #131 Probe A — 1G 2-region vacuum isolation.

Created by numerics-investigator on 2026-04-23.

Tests whether `solve_peierls_mg(SLAB_POLAR_1D, ...)` agrees with
`peierls_slab.solve_peierls_eigenvalue` on the SIMPLEST 2-region
slab: 1 group, vacuum BC. Isolates multi-region handling from the
MG × closure interaction.

Expected outcome
----------------
- If ~1e-8 agreement: multi-region volume kernel is fine; bug sits
  in the MG × multi-region × F.4 closure interaction.
- If ~1% or larger: multi-region volume kernel itself has a subtle
  bug (ray walker, material-interface breakpoints).

If this test catches a real bug, promote to tests/derivations/
test_peierls_multigroup.py or test_peierls_rank2_bc.py.
"""
from __future__ import annotations

import numpy as np
import pytest


@pytest.mark.slow
def test_issue131_probe_a_1g_2rg_vacuum():
    from orpheus.derivations.peierls_geometry import (
        SLAB_POLAR_1D,
        solve_peierls_mg,
    )
    from orpheus.derivations.peierls_slab import solve_peierls_eigenvalue

    # Synthetic 1G, 2-region slab. Region A: dense. Region B: moderator.
    # Matches the shipped 2eg_2rg fixture thicknesses.
    thicknesses = [0.5, 0.5]
    sig_t_A = np.array([0.5])
    sig_t_B = np.array([0.6])
    nu_sf_A = np.array([0.3])
    nu_sf_B = np.array([0.0])
    sig_s_A = np.array([[0.3]])  # sig_s[src, dst]
    sig_s_B = np.array([[0.5]])
    chi_A = np.array([1.0])
    chi_B = np.array([1.0])

    # Native path (E1 Nyström)
    sol_nat = solve_peierls_eigenvalue(
        sig_t_regions=[sig_t_A, sig_t_B],
        sig_s_matrices=[sig_s_A, sig_s_B],
        nu_sig_f_all=[nu_sf_A, nu_sf_B],
        chi_all=[chi_A, chi_B],
        thicknesses=thicknesses,
        n_panels_per_region=2,
        p_order=3,
        precision_digits=20,
        boundary="vacuum",
    )

    # Unified MG path — ng=1 with vacuum
    sig_t_mg = np.stack([sig_t_A, sig_t_B])  # (2, 1)
    sig_s_mg = np.stack([sig_s_A, sig_s_B])  # (2, 1, 1)
    nu_sig_f_mg = np.stack([nu_sf_A, nu_sf_B])  # (2, 1)
    chi_mg = np.stack([chi_A, chi_B])  # (2, 1)
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
    print(f"\nProbe A (1G 2-region vacuum):")
    print(f"  native  k_eff = {sol_nat.k_eff:.12f}")
    print(f"  unified k_eff = {sol_mg.k_eff:.12f}")
    print(f"  rel_diff      = {rel:.3e}")

    assert rel < 1e-6, (
        f"Probe A 1G 2-region vacuum: rel_diff={rel:.3e} "
        f"(native={sol_nat.k_eff:.12f}, unified={sol_mg.k_eff:.12f})"
    )


if __name__ == "__main__":
    test_issue131_probe_a_1g_2rg_vacuum()
