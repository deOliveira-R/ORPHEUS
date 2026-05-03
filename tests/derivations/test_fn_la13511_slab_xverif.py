r"""L1 cross-check: F_N slab vs Variant α slab at Sood ``Ua-1-0-SL``.

This is the load-bearing **structural-independence** verification
between two methods that share only ``numpy``/``scipy`` (above the
trusted-library line):

* **F_N method** (this slice's Branch-2 solver
  :func:`solve_fn_slab_bare_critical`) — Case singular eigenfunctions
  + Wiener-Hopf factorization (slab F_N collocation per Siewert-
  Benoist Part I + Grandjean-Siewert Part II).
* **Variant α method** (existing
  :func:`solve_greens_function_slab` at :math:`\alpha = 0`) —
  bouncing-trajectory angle-resolved Green's-function operator
  (Sanchez 1986 specular-leg machinery).

Both methods produce the same Sood/KLL truth :math:`a_c =
0.93772556` mfp via mathematically disjoint paths; their agreement is
the strongest available L1 evidence for either solver.

The agreement test uses the Sood Ua-1-0-SL XS at the published
critical half-thickness :math:`a = 2.872934` cm (so the slab full
width is :math:`L = 5.745868` cm), and verifies that BOTH methods
return :math:`k_{\rm eff} \approx 1.0` to within their respective
quadrature floors.

Quadrature notes
----------------

* F_N at :math:`N = 10` reaches :math:`a_c` truth to ~5e-6.
* Variant α slab at :math:`(n_x, n_\mu) = (48, 128)` reaches
  :math:`k_{\rm eff} = 1` to ~1e-5 at the F_N truth thickness. The
  angular-quadrature floor is the bottleneck — the slab Variant α
  uses GL on :math:`[-1, 1]` for :math:`\mu`, and slab vacuum BC
  has a near-cusp at :math:`\mu = 0` that takes ~128 nodes to
  resolve to 1e-5.

Cross-check tolerance is set to 5e-5 to accommodate both floors.
"""
from __future__ import annotations

import pytest

from orpheus.derivations.continuous.sood_registry import UA_1_0_SL_STUB
from orpheus.derivations.continuous.fn_method.slab import (
    solve_fn_slab_bare_critical,
)
from orpheus.derivations.continuous.peierls_greens_function.greens_function_slab import (
    solve_greens_function_slab,
)


# Suppress the bracket-scan divide-by-zero warnings (see slab tests).
pytestmark = [
    pytest.mark.filterwarnings(
        "ignore:divide by zero encountered in det:RuntimeWarning"
    ),
    pytest.mark.filterwarnings(
        "ignore:invalid value encountered in det:RuntimeWarning"
    ),
]


@pytest.mark.l1
@pytest.mark.slow
def test_fn_slab_vs_variant_alpha_at_sood_ua_1_0_sl():
    r"""Both methods return :math:`k_{\rm eff} = 1` at the Sood Ua-1-0-SL
    truth dimension :math:`L = 5.745868` cm.

    F_N runs at :math:`N = 10` (its inherent convergence floor for this
    case is ~5e-6 in :math:`a_c`); Variant α runs at
    :math:`(n_x, n_\mu, n_{\rm traj}) = (48, 128, 96)` (~1e-5
    floor). The agreement is required to ≤ 5e-5 — the larger of
    the two floors.
    """
    case = UA_1_0_SL_STUB
    sigma_t = float(case.sigma_t[0])
    sigma_s = float(case.sigma_s[0, 0])
    nu_sigma_f = float(case.nu_sigma_f[0])

    # F_N: returns a_c in mfp; convert to cm via /sigma_t.
    res_fn = solve_fn_slab_bare_critical(c=(sigma_s + nu_sigma_f) / sigma_t,
                                          n_modes=10)
    a_fn_mfp = res_fn.a_critical_mfp
    a_fn_cm = a_fn_mfp / sigma_t  # 2.872xx cm
    L_full_cm = 2.0 * a_fn_cm

    # F_N self-consistency: should match Sood truth.
    truth_mfp = case.critical_dimension_mfp
    err_fn = abs(a_fn_mfp - truth_mfp)
    assert err_fn < 5e-5, (
        f"F_N internal: a_c={a_fn_mfp:.10f} vs truth {truth_mfp}, "
        f"err={err_fn:.2e}"
    )

    # Variant α at the F_N-truth thickness should give k_eff = 1.
    res_va = solve_greens_function_slab(
        L=L_full_cm,
        sigma_t=sigma_t,
        sigma_s=sigma_s,
        nu_sigma_f=nu_sigma_f,
        alpha=0.0,                  # vacuum BC
        n_x=48,
        n_mu=128,
        n_traj_quad=96,
        max_iter=500,
        tol=1e-9,
    )
    k_eff_at_truth = res_va.k_eff
    err_va = abs(k_eff_at_truth - 1.0)
    assert err_va < 5e-5, (
        f"Variant α at F_N truth thickness: k_eff={k_eff_at_truth:.8f}, "
        f"err={err_va:.2e}"
    )

    # The two-method cross-check: both report critical-thickness
    # agreement at the level of their respective floors. If one method
    # is biased by O(1e-3) the assertion above will fail and the bug
    # source can be pinned to whichever method moved.


@pytest.mark.l1
@pytest.mark.slow
def test_fn_slab_vs_variant_alpha_at_grandjean_siewert_c150():
    r"""Cross-check at a different :math:`c` value — Grandjean-Siewert
    Table XI :math:`c = 1.50` row. Both methods should agree at the
    same critical thickness to ≤ 1e-4 (looser at higher :math:`c` due
    to thinner slab + steeper boundary gradient).
    """
    c = 1.50
    # Solve via F_N (unit XS).
    res_fn = solve_fn_slab_bare_critical(c=c, n_modes=10)
    a_fn = res_fn.a_critical_mfp  # in mfp = cm at sigma_t=1.

    # Variant α at L = 2*a_fn using unit XS.
    res_va = solve_greens_function_slab(
        L=2.0 * a_fn,
        sigma_t=1.0,
        sigma_s=0.0,
        nu_sigma_f=c,             # so c_input = c
        alpha=0.0,
        n_x=48,
        n_mu=128,
        n_traj_quad=96,
        max_iter=500,
        tol=1e-9,
    )
    err = abs(res_va.k_eff - 1.0)
    assert err < 1e-4, (
        f"GS c=1.50: F_N a={a_fn:.10f}, Variant α k_eff at that "
        f"thickness = {res_va.k_eff:.8f}, err={err:.2e}"
    )
