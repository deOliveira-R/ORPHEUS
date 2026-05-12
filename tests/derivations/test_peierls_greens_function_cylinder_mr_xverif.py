r"""Phase 1b cylinder Variant α MR — Gate 2 WM-72 vacuum-BC cross-check.

Structurally-independent L1 reference for cylinder MR Variant α. The
Westfall-Metcalf 1972 / 1973 F_N method (singular eigenfunctions +
Wiener-Hopf factorisation) is in a completely different mathematical
framework from the angle-resolved Green's function along bouncing
characteristics that Variant α uses. Their agreement at the
single-region MR vacuum-BC critical radius is the strongest
structurally-independent L1 evidence available for cylinder MR.

Gate 2 separated from
``test_peierls_greens_function_cylinder_mr.py`` because:

- Longest runtime (Brent root-search inside WM-72 + cylinder MR power
  iteration at adequate quadrature).
- L1 anchor — filename signals "this is the structurally-independent
  L1 evidence" per the conventional ``_xverif_`` suffix used by
  ``test_peierls_greens_function_cylinder_xverif_sood2003.py`` and
  ``test_peierls_greens_function_xverif_ps1982.py``.

Gate 5 (literature MR benchmark) status: per the verification plan
§3 Gate 5 search summary, no published cylinder MR benchmark of the
appropriate scope (concentric annular regions + isotropic scattering
+ specular / vacuum BC + tabulated k_eff or φ(r)) was identified
during the Phase 1b verification work. The chain of structural
independence for cylinder MR Variant α terminates at:

  - **Gate 2** (this file) — WM-72 single-region vacuum critical
    radius
  - **Gate 3** (closed-form analytical k_∞)
  - **Gate 4** (interface continuity invariant)
  - **Gate 6** (quadrature refinement)

This is documented as a known V&V limitation; a future search of
Russian and Japanese journal archives, plus the Sanchez-Pomraning
textbook (1989) and Calvik 1965/1967 papers, might surface a
shippable cylinder MR benchmark.

References
----------

- Westfall, R. M. & Metcalf, D. R. (1972). *Nucl. Sci. Eng.* 49, 273
  — F_N method derivation (the Wiener-Hopf factorization that
  underlies :func:`solve_singular_eigenfunction_cylinder_bare_critical`).
- Sood, Forster & Parsons (1999/2003), LA-13511 Table 13 Ua-1-O-CY
  benchmark (the same critical-radius value WM-72 produces).
- Verification plan:
  ``.claude/plans/cylinder_mr_variant_alpha_verification.md``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.singular_eigenfunction import (
    solve_singular_eigenfunction_cylinder_bare_critical,
)
from orpheus.derivations.continuous.trajectory_resolvent.greens_function_cylinder import (
    solve_greens_function_cylinder_mr,
)


# ───────────────────────────────────────────────────────────────────────
# WM-72 / Sood Ua-1-O-CY constants — c = 1.30 bare-critical cylinder
# ───────────────────────────────────────────────────────────────────────

SIGMA_T = 0.32640        # cm⁻¹
SIGMA_S = 0.248064       # cm⁻¹
NU_SIGMA_F = 0.176256    # cm⁻¹
# c = (Σ_s + νΣ_f) / Σ_t = 1.30


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("peierls-greens-cylinder-mr-wm72-vacuum")
def test_mr_single_region_vacuum_matches_wm72():
    r"""Gate 2 — cylinder MR (K=1, α=0) at WM-72 critical radius gives
    :math:`k_{\rm eff} = 1` to ≤ 1e-5.

    Branch independence: WM-72's F_N method uses Case singular
    eigenfunctions + Wiener-Hopf factorization of the half-space
    transport equation. **No Bickley-Naylor / Ki_n integrals appear
    anywhere in the F_N derivation.** Cylinder Variant α MR uses
    angle-resolved Green's-function integration along bouncing
    characteristics. These are mathematically distinct frameworks
    (different structural ground per ``vv-principles`` § "Structural
    independence"); agreement implies both are converging to the same
    true critical radius.

    Anti-pattern protection per plan §3 Gate 2: do NOT compare against
    :func:`solve_greens_function_cylinder_mg` with α=0 — that is
    procedurally independent but structurally identical (same Variant
    α algorithm). The WM-72 cross IS the structurally-independent
    pillar.

    Quadrature ``(n_r, n_mu_axial, n_phi_az, n_traj_quad) =
    (24, 20, 48, 64)`` calibrated to give ≤ 1e-5 agreement at the
    6-digit-accurate WM-72 reference (empirically 4.8e-6 in
    pre-implementation sanity check).
    """
    wm = solve_singular_eigenfunction_cylinder_bare_critical(
        c=1.30, sigma_t=SIGMA_T, n_grid=24,
    )
    assert wm.converged, (
        f"Gate 2 (WM-72): bare-critical solver did not converge; "
        f"residual = {wm.criticality_residual}"
    )
    r_c = wm.r_c_cm
    # Sanity: WM-72 should land on the Sood 6-digit reference 5.284935
    np.testing.assert_allclose(
        r_c, 5.284935, rtol=2e-6,
        err_msg=(
            f"Gate 2 (WM-72): r_c = {r_c} cm differs from Sood reference "
            f"5.284935 cm by more than 2e-6 — WM-72 solver regression."
        ),
    )

    # Drive MR (K=1) at WM-72 r_c with α=0 (vacuum)
    res = solve_greens_function_cylinder_mr(
        radii=np.array([r_c]),
        sigma_t=np.array([SIGMA_T])[None, :],
        sigma_s=np.array([SIGMA_S])[None, None, :],
        nu_sigma_f=np.array([NU_SIGMA_F])[None, :],
        alpha=0.0,
        n_r=24, n_mu_axial=20, n_phi_az=48, n_traj_quad=64,
        max_iter=400, tol=1e-12,
    )
    assert res.converged, (
        f"Gate 2 (WM-72 cross): MR at K=1, α=0, R = {r_c:.6f} cm did "
        f"not converge; k_eff = {res.k_eff}, iters = {res.iterations}"
    )
    # Target ≤ 1e-5 per plan §3 Gate 2.
    np.testing.assert_allclose(
        res.k_eff, 1.0, atol=1e-5,
        err_msg=(
            f"Gate 2 (WM-72 cross): MR at K=1, α=0, R = {r_c:.6f} cm "
            f"(WM-72 critical radius reference) gave k_eff = "
            f"{res.k_eff:.16e}; |k_eff - 1| = {abs(res.k_eff - 1.0):.3e}, "
            f"exceeds the 1e-5 target tolerance. A divergence here "
            f"after Gate 1 passes would indicate a vacuum-branch bug "
            f"specific to the MR α=0 codepath (suspect first-leg-only "
            f"integration or piecewise τ accumulation in the α=0 "
            f"early-return)."
        ),
    )
