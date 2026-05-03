r"""L1 cross-check: sphere F_N method (Siewert-Thomas 1986) vs Variant α
sphere bouncing-trajectory operator at Sood ``Ua-1-0-SP``.

Structural-independence pillar
------------------------------

Both branches reach the same critical-sphere prediction by genuinely
disjoint mathematical paths above the trusted-library line:

* **Sphere F_N (this branch's reference)**: Wiener-Hopf factorisation
  of the Case dispersion function plus half-range moment expansion.
  The collocation system :math:`\sum a_\alpha M_{\beta,\alpha}(R) = 0`
  with :math:`M_{\beta,\alpha} = B_\alpha(\xi_\beta) - e^{-2R/\xi_\beta}
  A_\alpha(\xi_\beta)`. Provides the critical condition via
  :math:`\det M(R_c) = 0`. Works in the **Case singular-eigenfunction
  representation** of the angular flux.

* **Variant α sphere (this branch's system-under-test)**: angle-
  resolved Green's function with bouncing-trajectory closure
  (rank-1 boundary-to-boundary scattering operator). Reduces to a
  fixed-point iteration on the surface inflow. Works in the
  **bouncing-characteristic representation** of the angular flux.

These representations agree on the same physics (the 1G isotropic-
scattering bare-critical sphere) but use no shared in-house code —
only ``numpy``, ``scipy.special``, and ``scipy.linalg`` at the
trusted-library level. This makes their agreement at :math:`R_c` a
genuine **L1 cross-check**, NOT cross-implementation agreement (L4).

The previous test in this file compared PS-1982 wrapper (Peierls
integral equation) vs Variant α (bouncing-trajectory). Those two
methods are NOT genuinely structurally independent — both reduce
the same Peierls integral equation by different algebraic paths.
The F_N method is genuinely independent because it works in the
Case eigenfunction basis, never reducing to an integral equation
in :math:`r`.

References
----------

* Siewert & Thomas 1986, *Nucl. Sci. Eng.* **94**, 264.
* Sanchez 1986, *J. Quant. Spec. Rad. Transfer* **35**, 121 (Variant α
  specular leg).
* Sood, Forster & Parsons 1999, LA-13511 (Table 6).
"""
from __future__ import annotations

import pytest

# Suppress numpy divide-by-zero warnings from the F_N bracket scan.
pytestmark = [
    pytest.mark.filterwarnings(
        "ignore:divide by zero encountered in det:RuntimeWarning"
    ),
    pytest.mark.filterwarnings(
        "ignore:invalid value encountered in det:RuntimeWarning"
    ),
]

from orpheus.derivations.continuous.fn_method.benchmarks import UA_1_0_SP_STUB
from orpheus.derivations.continuous.fn_method.sphere import (
    solve_fn_sphere_bare_critical,
)
from orpheus.derivations.continuous.peierls_greens_function.greens_function import (
    solve_greens_function_sphere,
)


@pytest.mark.l1
def test_fn_sphere_vs_variant_alpha_sphere_at_sood_ua_1_0_sp():
    r"""L1: F_N sphere (Wiener-Hopf via Case eigenfunctions) and Variant
    α sphere (bouncing-trajectory operator) agree on the bare-critical
    radius at Sood ``Ua-1-0-SP`` (c=1.30) to ≤ 1e-5 absolute.

    The two methods are structurally independent above the trusted-
    library line: F_N works in the Case singular-eigenfunction
    representation; Variant α works in the bouncing-characteristic
    representation. Their agreement is genuine L1 evidence for the
    correctness of both — neither is L4 cross-implementation.

    Direct comparison: F_N predicts :math:`R_c^{\rm FN}`; Variant α
    at :math:`R = R_c^{\rm FN}` should give :math:`k_{\rm eff} = 1`.
    Equivalently: solve both for c = 1.30 critical radius
    independently, compare.
    """
    # F_N method at moderate N (N=10 already reaches ~1e-7 absolute
    # vs Sood truth, well below the 1e-5 cross-check tolerance).
    res_fn = solve_fn_sphere_bare_critical(c=1.30, n_modes=10)
    R_c_fn = res_fn.R_critical_mfp

    # Variant α sphere at the F_N predicted radius. We expect k_eff = 1
    # to within 1e-5 (the structural-independence cross-check pillar).
    case = UA_1_0_SP_STUB
    sigma_t = float(case.sigma_t[0])
    sigma_s = float(case.sigma_s[0, 0])
    nu_sigma_f = float(case.nu_sigma_f[0])
    # Convert F_N R_c (mfp) to cm using the case's σ_t.
    R_c_cm = R_c_fn / sigma_t

    res_va = solve_greens_function_sphere(
        R=R_c_cm,
        sigma_t=sigma_t,
        sigma_s=sigma_s,
        nu_sigma_f=nu_sigma_f,
        alpha=0.0,
        n_r=32,
        n_mu=32,
        n_traj_quad=64,
        max_iter=400,
        tol=1e-10,
    )
    err = abs(res_va.k_eff - 1.0)
    assert err < 1e-5, (
        f"F_N sphere R_c = {R_c_fn:.10f} mfp = {R_c_cm:.6f} cm; "
        f"Variant α at this R gives k_eff = {res_va.k_eff:.8f}, "
        f"|k - 1| = {err:.3e} (target 1e-5). "
        f"F_N vs Sood Ua-1-0-SP truth ({case.critical_dimension_mfp}): "
        f"err = {abs(R_c_fn - case.critical_dimension_mfp):.3e}"
    )


@pytest.mark.l1
def test_fn_sphere_matches_sood_ua_1_0_sp_directly():
    r"""L1: F_N sphere reaches Sood truth at ≤ 1e-5 directly (no
    cross-check, just the published reference).

    This duplicates the foundation gate in test_fn_la13511_sphere
    but at the L1 level — for the V&V audit harness it documents
    "F_N method has L1 cross-check evidence at Sood Ua-1-0-SP",
    which is the structurally-independent pillar.
    """
    res = solve_fn_sphere_bare_critical(c=1.30, n_modes=10)
    truth = UA_1_0_SP_STUB.critical_dimension_mfp
    err = abs(res.R_critical_mfp - truth)
    assert err < 1e-5, (
        f"F_N sphere vs Sood Ua-1-0-SP truth: err = {err:.3e}"
    )


@pytest.mark.l1
def test_variant_alpha_sphere_at_sood_truth_radius():
    r"""L1: Variant α sphere at the Sood Ua-1-0-SP TRUTH radius
    :math:`R = 7.428998` cm gives :math:`k_{\rm eff} = 1` to ≤ 5e-5.

    Held over from the previous (PS-1982 wrapper) cross-check — the
    Variant α prediction at the truth radius is the SAME claim as
    before (Variant α has not changed). The structural-independence
    upgrade applies to the reference, not the system-under-test.
    """
    case = UA_1_0_SP_STUB
    sigma_t = float(case.sigma_t[0])
    sigma_s = float(case.sigma_s[0, 0])
    nu_sigma_f = float(case.nu_sigma_f[0])
    R_truth_cm = case.critical_dimension_cm  # 7.428998

    res_va = solve_greens_function_sphere(
        R=R_truth_cm,
        sigma_t=sigma_t,
        sigma_s=sigma_s,
        nu_sigma_f=nu_sigma_f,
        alpha=0.0,
        n_r=32,
        n_mu=32,
        n_traj_quad=64,
        max_iter=400,
        tol=1e-10,
    )
    err = abs(res_va.k_eff - 1.0)
    assert err < 5e-5, (
        f"Variant α sphere at Sood truth R={R_truth_cm} cm: "
        f"k_eff={res_va.k_eff:.8f}, err={err:.3e}"
    )
