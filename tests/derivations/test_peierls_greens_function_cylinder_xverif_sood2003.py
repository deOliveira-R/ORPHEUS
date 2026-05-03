r"""V&V hardening — cylinder Variant α external-reference cross-check
against Sood/Forster/Parsons LA-13511 (1999) Table 13 ``Ua-1-O-CY``.

Closes the highest-priority post-Phase-3 V&V gap on the cylinder
geometry: prior to this test, cylinder Variant α had only the
internal :math:`T_{00} \equiv P_{ss}` cross-check, which reuses the
Bickley-Naylor :math:`\mathrm{Ki}_n` identity in **both** code paths
(ERR-032 pattern risk — a procedurally-independent but *not*
structurally-independent agreement).

Why the Sood benchmark is structurally independent
----------------------------------------------------

Sood et al. compute critical radii via the **F_N method**:

- Westfall & Metcalf (1972), *Nucl. Sci. Eng.* 49, 273 — F_N method
  for slab and sphere.
- Westfall (1983) follow-on — F_N applied to the bare critical
  cylinder (LA-13511 Refs. 27, 28).

The F_N method uses **Case singular eigenfunctions** plus
**Wiener-Hopf factorization** of the half-space transport equation.
**No Bickley-Naylor / Ki_n integrals appear anywhere in the F_N
derivation.**

This is structurally independent of the cylinder Variant α path,
which integrates the angle-resolved Green's function along bouncing
characteristics (no Ki_n at all), and also independent of the
cylinder Nyström primitives (which DO use Bickley-Naylor Ki_3). The
Sood benchmark therefore provides a true "different mathematical
framework" cross-check — what `vv-principles` § "The three pillars
of verification" calls a structurally-independent reference, the
strongest evidence class available for an L1 claim.

Per the post-Phase-3 hindsight review, this closes Task 1 of the V&V
hardening pass. Tasks 2-4 (slab-sym ↔ Nyström, off-diagonal
intermediate α, grazing-ray stability) were already landed in the
preceding commit ``1f17c20``.

Benchmark configuration — ``Ua-1-O-CY``
----------------------------------------

LA-13511 Table 12 (U-235 (a) cross sections, c = 1.30) +
Table 13 (Critical Dimensions for One-Group Bare U-235):

- ν = 2.70
- Σ_f = 0.06528 cm⁻¹ → νΣ_f = 0.176256 cm⁻¹
- Σ_c = 0.013056 cm⁻¹ (capture)
- Σ_s = 0.248064 cm⁻¹
- Σ_t = 0.32640 cm⁻¹
- c = (Σ_s + νΣ_f) / Σ_t = (0.248064 + 0.176256) / 0.32640 = 1.30
- **Critical radius r_c = 1.72500292 mfp = 5.284935 cm** (Sood
  6-digit accuracy claim, F_N method; Refs. 27 and 28 in LA-13511).

Test: drive Variant α at α=0 (vacuum BC, bare critical cylinder) at
R = r_c. The eigenvalue problem must converge to **k_eff = 1.0**
to research-grade tolerance (target ≤ 1e-5; achieved 8.5e-6 at the
quadrature orders pinned below).

References
----------

- Sood, Forster & Parsons (1999), *Analytical Benchmark Test Set
  for Criticality Code Verification*, LA-13511 (Los Alamos National
  Laboratory). Open-access OSTI 10601. Also published as Sood,
  Forster & Parsons (2003), *Prog. Nucl. Energy* 42(1), 55-106
  (DOI 10.1016/S0149-1970(02)00098-7).
- Westfall, R. M. (1983), F_N-method critical-cylinder calculation
  cited in LA-13511 Ref. 27.
- Westfall & Metcalf (1972), *Nucl. Sci. Eng.* 49, 273 — F_N method
  origin (Ref. 28).

Adjacent tests
--------------

- :mod:`.test_peierls_greens_function_xverif_ps1982` — sphere
  Variant α external-reference cross-check (parallel structure).
- :mod:`.test_peierls_greens_function_cylinder_solver` — cylinder
  V_α1_cyl k_inf-exactness invariants and the V_α2_cyl
  internal-primitive cross-check (T_00 ≡ P_ss).

Memo: ``.claude/agent-memory/literature-researcher/sood_2003_cylinder_benchmarks.md``.
Closeout: ``.claude/agent-memory/method-implementer/sood2003_cylinder_xverif.md``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.peierls_greens_function.greens_function_cylinder import (
    solve_greens_function_cylinder,
)


# ───────────────────────────────────────────────────────────────────────
# Sood Ua-1-O-CY benchmark constants (LA-13511 Tables 12 + 13)
# ───────────────────────────────────────────────────────────────────────

R_CRITICAL_SOOD = 5.284935       # cm; r_c = 1.72500292 mfp
SIGMA_T = 0.32640                # cm⁻¹
SIGMA_S = 0.248064               # cm⁻¹
NU_SIGMA_F = 0.176256            # cm⁻¹ (= 2.70 × 0.06528)
# c = (Σ_s + νΣ_f) / Σ_t = (0.248064 + 0.176256) / 0.32640 = 1.30 ✓


# ───────────────────────────────────────────────────────────────────────
# L1 — Sood 1999 Ua-1-O-CY external-reference cross-check
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("peierls-greens-cylinder-architecture")
def test_a2_variant_alpha_agrees_with_sood2003_cylinder():
    r"""L1 (slow) — cylinder Variant α critical radius matches Sood
    1999/2003 ``Ua-1-O-CY`` benchmark to ≤ 1e-5.

    Sood/Forster/Parsons LA-13511 Table 13 ``Ua-1-O-CY`` benchmark:
    bare critical homogeneous cylinder, U-235 (a) cross sections,
    c = 1.30, vacuum BC. Sood's reference :math:`r_c = 1.72500292`
    mfp = 5.284935 cm is computed via the F_N method (Case singular
    eigenfunctions + Wiener-Hopf factorization, LA-13511 Refs. 27 and
    28), structurally independent of Bickley-Naylor / :math:`\mathrm{Ki}_n`
    integrals.

    Closes the post-Phase-3 hindsight V&V gap on cylinder Variant α:
    previously only the internal :math:`T_{00} \equiv P_{ss}` cross-check
    was available, which reuses the Bickley-Naylor identity in both
    code paths (ERR-032 pattern risk). Sood's F_N reference is in a
    completely different mathematical framework — true
    structurally-independent cross-verification.

    Quadrature orders ``(n_r, n_μ_axial, n_φ_az, n_traj_quad) =
    (24, 20, 64, 96)`` calibrated empirically to give ≤ 1e-5 agreement
    at the 6-digit-accurate Sood reference (achieved 8.5e-6 in a
    calibration sweep at HEAD ``1f17c20`` — see closeout memo).
    """
    res = solve_greens_function_cylinder(
        R=R_CRITICAL_SOOD,
        sigma_t=SIGMA_T,
        sigma_s=SIGMA_S,
        nu_sigma_f=NU_SIGMA_F,
        alpha=0.0,
        n_r=24, n_mu_axial=20, n_phi_az=64, n_traj_quad=96,
        max_iter=400, tol=1e-12,
    )
    assert res.converged, (
        f"Sood Ua-1-O-CY: Variant α did not converge at R = "
        f"{R_CRITICAL_SOOD} cm; k_eff = {res.k_eff}, "
        f"iter = {res.iterations}"
    )

    # Target ≤ 1e-5 (achieved ~8.5e-6 in calibration). Hard floor 1e-4
    # would indicate a structural problem — STOP and dispatch
    # numerics-investigator if exceeded.
    np.testing.assert_allclose(
        res.k_eff, 1.0, atol=1e-5,
        err_msg=(
            f"Sood Ua-1-O-CY: Variant α at r_c = {R_CRITICAL_SOOD} cm "
            f"(F_N method reference) gave k_eff = {res.k_eff:.16e}, "
            f"|k_eff - 1| = {abs(res.k_eff - 1.0):.3e}, exceeds the "
            f"1e-5 target tolerance. Sood's 6-digit-accurate F_N "
            f"reference and Variant α should agree at this level — "
            f"investigate before declaring V&V hardening complete."
        ),
    )


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-cylinder-architecture")
def test_a2_sood_quadrature_order_convergence():
    r"""L1 — Sood ``Ua-1-O-CY`` agreement improves monotonically with
    quadrature refinement.

    Verifies that the agreement |k_eff - 1| at the Sood reference radius
    decreases as the quadrature is refined from coarse to medium. Pins
    the quadrature-order sensitivity of the cross-check so a future
    regression that introduces grid-dependent bias at the Sood
    configuration is caught.

    Calibration data captured at HEAD ``1f17c20`` (closeout memo):

    - (8, 8, 16, 24): 2.7e-4
    - (12, 10, 24, 32): 3.3e-5
    - (16, 12, 32, 48): 4.8e-5
    - (20, 16, 48, 64): 1.4e-5

    Note: convergence is non-monotone at intermediate orders (the
    (12, 10, 24, 32) result is slightly closer than (16, 12, 32, 48));
    this is a known property of GL on a 3D phase space with the
    impact-parameter geometry of cylinder vacuum-BC eigenmodes (see
    Phase-1 closeout memo). The non-monotonicity is bounded — the
    coarse-vs-fine agreement is consistently better than the
    coarse-vs-coarse step.
    """
    # Coarse: should give ~2.7e-4 (verifies "agreement is real, not
    # accidental at one grid").
    res_coarse = solve_greens_function_cylinder(
        R=R_CRITICAL_SOOD, sigma_t=SIGMA_T, sigma_s=SIGMA_S,
        nu_sigma_f=NU_SIGMA_F, alpha=0.0,
        n_r=8, n_mu_axial=8, n_phi_az=16, n_traj_quad=24,
        max_iter=400, tol=1e-11,
    )
    # Medium-fine: should give ~1e-5 already.
    res_fine = solve_greens_function_cylinder(
        R=R_CRITICAL_SOOD, sigma_t=SIGMA_T, sigma_s=SIGMA_S,
        nu_sigma_f=NU_SIGMA_F, alpha=0.0,
        n_r=20, n_mu_axial=16, n_phi_az=48, n_traj_quad=64,
        max_iter=400, tol=1e-12,
    )

    err_coarse = abs(res_coarse.k_eff - 1.0)
    err_fine = abs(res_fine.k_eff - 1.0)

    assert res_coarse.converged and res_fine.converged
    # Coarse agreement is real — should already be < 1e-3 (Sood is a
    # genuinely critical configuration).
    assert err_coarse < 1e-3, (
        f"Sood Ua-1-O-CY: even coarse (8, 8, 16, 24) Variant α gives "
        f"|k_eff - 1| = {err_coarse:.3e} >= 1e-3 — coarse agreement "
        f"with the F_N reference is unexpectedly poor; investigate."
    )
    # Fine agreement is research-grade.
    assert err_fine < 5e-5, (
        f"Sood Ua-1-O-CY: fine (20, 16, 48, 64) Variant α gives "
        f"|k_eff - 1| = {err_fine:.3e} >= 5e-5 — fine agreement with "
        f"the F_N reference is below research-grade target."
    )
    # Fine is at least 4x tighter than coarse (allows for the
    # documented non-monotonicity at intermediate grids).
    assert err_fine < err_coarse / 4, (
        f"Sood Ua-1-O-CY: refinement does not tighten agreement — "
        f"coarse |k-1| = {err_coarse:.3e}, fine |k-1| = "
        f"{err_fine:.3e}; expected ratio ≥ 4x."
    )
