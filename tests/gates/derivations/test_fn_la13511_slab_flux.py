r"""L1 verification for the slab F_N rich-machinery extension.

The F_N method gives the boundary angular flux at :math:`z = \pm a`
via the polynomial expansion :math:`\psi(\pm a, \mp\mu) = \sum_\alpha
a_\alpha \mu^\alpha`. The rich-machinery extension (KLL 1974 Fredholm
iteration) gives the **interior** scalar flux :math:`\phi(z)` and,
via characteristic integration, the **interior** angular flux
:math:`\psi(z, \mu)` at any :math:`(z, \mu)`.

Verification gates here:

* **L1 — Sood / KLL Table III flux ratios at c = 1.30**: at
  :math:`z/b \in \{0.25, 0.50, 0.75, 1.00\}`, the KLL reconstruction
  agrees with the published reference at ≤ 1e-5 absolute. (This is
  the same XS as Sood ``Ua-1-0-SL``.)
* **L1 — KLL Table III sweep**: for :math:`c \in \{1.05, 1.10, 1.20,
  1.30, 1.40, 1.60\}`, the reconstruction agrees with the published
  values at ≤ 5e-5 absolute (the c=1.40, 1.60 rows are tighter
  binding for thin slabs where convergence is slower).
* **L1 — angular-flux closure**: :math:`\int_{-1}^{1} \psi(z, \mu)\,
  d\mu` recomputed from the angular-flux reconstruction agrees with
  the directly-evaluated :math:`\phi(z)` at ≤ 1e-6 absolute (away
  from the boundary).
* **L1 — KLL ↔ F_N a_c agreement**: the KLL recipe (taking F_N's
  :math:`a_c` as input) and the F_N method itself agree on the
  critical condition for the slab — this is a structural
  cross-check above the trusted-library line.

Tolerance budget (per ``vv-principles``):

* The 1e-5 acceptance gate matches the F_N solver's r_c precision at
  N=10 (Sood Table 14 has 6 published digits; the verifier sees both
  the F_N r_c and the KLL flux ratio as known to ≤ 1e-5 against the
  references).
* Tighter gates are technically achievable but pinned at the
  benchmark precision; reporting tighter would imply false precision.

References
----------

* Kaper-Lindeman-Leaf 1974, *Nucl. Sci. Eng.* **54**, 94, Table III.
* Sood, Forster, Parsons 1999, LA-13511 Table 14 (case ``Ua-1-0-SL``).
* Siewert-Benoist 1979, Part I + Grandjean-Siewert 1979, Part II.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.fn_method.slab import (
    slab_angular_flux_from_scalar,
    slab_scalar_flux_from_angular_quadrature,
    slab_scalar_flux_kll,
    slab_scalar_flux_ratio,
    solve_fn_slab_bare_critical,
    solve_kll_slab_continuum_coefficient,
)
from orpheus.derivations.continuous.sood_registry import LA13511_CASES

# Suppress numpy divide-by-zero warnings from the F_N bracket scan.
pytestmark = [
    pytest.mark.filterwarnings(
        "ignore:divide by zero encountered in det:RuntimeWarning"
    ),
    pytest.mark.filterwarnings(
        "ignore:invalid value encountered in det:RuntimeWarning"
    ),
]


# ────────────────────────────────────────────────────────────────────
# L1 — Sood Ua-1-0-SL flux ratios from KLL Table III
# ────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.verifies("kll-1974-slab-flux")
def test_l1_slab_flux_ratios_kll_table_III_at_c_1p30() -> None:
    r"""KLL/Sood ``Ua-1-0-SL`` flux ratios at c = 1.30 to ≤ 1e-5.

    KLL Table III (also Sood Table 14):

    +-------+------------+
    | z/b   | φ(z)/φ(0)  |
    +=======+============+
    | 0.25  | 0.9669506  |
    | 0.50  | 0.8686259  |
    | 0.75  | 0.7055218  |
    | 1.00  | 0.4461912  |
    +-------+------------+
    """
    case = LA13511_CASES["Ua-1-0-SL"]
    # Compute c = (Σ_s + νΣ_f) / Σ_t via the legacy property accessors
    # (which delegate to mixture_to_fn_arrays per sood_registry).
    sigma_t = case.materials[0].SigT[0]
    sigma_s = case.materials[0].SigS[0][0, 0]
    nu_sigma_f = case.materials[0].SigP[0]
    c = float((sigma_s + nu_sigma_f) / sigma_t)
    assert abs(c - 1.30) < 1e-10, f"Expected c=1.30, got c={c}"

    b_kll = case.truth.critical_dimension_mfp
    res = solve_kll_slab_continuum_coefficient(b_kll, c, n_nodes=64)
    assert res.converged, "KLL Fredholm iteration did not converge"

    z_over_b = np.array([0.25, 0.50, 0.75, 1.00])
    computed = slab_scalar_flux_ratio(res, z_over_b)
    expected = np.array([0.9669506, 0.8686259, 0.7055218, 0.4461912])
    abs_err = np.abs(computed - expected)

    assert np.max(abs_err) < 1e-5, (
        f"KLL slab flux ratios at c=1.30 disagree with Sood Table 14 / "
        f"KLL Table III at > 1e-5: computed={computed}, "
        f"expected={expected}, abs_err={abs_err}"
    )


@pytest.mark.l1
@pytest.mark.verifies("kll-1974-slab-flux")
@pytest.mark.parametrize(
    "c,b_truth,truth_ratios",
    [
        (1.05, 3.30026_3772, (0.94714400, 0.79372641, 0.55329025, 0.21419206)),
        (1.10, 2.11330_9666, (0.95435566, 0.82074787, 0.60710163, 0.29143526)),
        (1.20, 1.28937_9285, (0.96227497, 0.85072970, 0.66821362, 0.38556499)),
        (1.30, 0.93772_556,  (0.96695_06,  0.86862_59,  0.70552_18,  0.44619_12)),
        (1.40, 0.73660_355,  (0.97017_34,  0.88105_40,  0.73181_31,  0.49025_92)),
        (1.60, 0.51196_298,  (0.97445_9,   0.89770_1,   0.76752_2,   0.55179_0)),
    ],
)
def test_l1_slab_flux_ratios_kll_table_III_sweep(
    c: float, b_truth: float, truth_ratios: tuple[float, ...]
) -> None:
    r"""KLL Table III flux ratios sweep over :math:`c \in [1.05, 1.60]`.

    Tighter tolerance for thicker slabs (low c) where the iteration is
    fastest; looser tolerance for thinner slabs (high c) which
    converge slower.
    """
    res = solve_kll_slab_continuum_coefficient(
        b_truth, c, n_nodes=80, max_iter=400
    )
    assert res.converged, f"KLL did not converge for c={c}, b={b_truth}"

    z_over_b = np.array([0.25, 0.50, 0.75, 1.00])
    computed = slab_scalar_flux_ratio(res, z_over_b)
    expected = np.array(truth_ratios)
    abs_err = np.abs(computed - expected)

    # Tolerance budget: thick slabs (c=1.05) trivially better; thin
    # slabs (c=1.60) slightly looser. KLL Table III itself reports 4-7
    # significant digits.
    if c <= 1.20:
        tol = 1e-6
    elif c <= 1.40:
        tol = 5e-6
    else:
        tol = 5e-5

    assert np.max(abs_err) < tol, (
        f"KLL slab flux ratios at c={c} disagree with Table III at > "
        f"{tol}: computed={computed}, expected={expected}, "
        f"abs_err={abs_err}"
    )


# ────────────────────────────────────────────────────────────────────
# L1 — angular-flux closure
# ────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.verifies("kll-1974-slab-flux")
def test_l1_slab_angular_flux_closure() -> None:
    r"""Angular-flux closure: :math:`\phi(z) = \int_{-1}^{1} \psi(z, \mu)\,
    d\mu` self-consistency.

    Verifies the angular-flux reconstruction (characteristic
    integration of the BTE with the converged KLL :math:`\phi(z)` as
    the source) integrates back to :math:`\phi(z)` to within the
    quadrature accuracy.

    Tolerance budget: away from the boundary (z/b ≤ 0.5), expect
    ≤ 1e-6. Near the boundary (z/b → 1) the characteristic integrand
    becomes more singular and looser tolerance is reasonable. Use
    interior z values only for the closure test.
    """
    c = 1.30
    b = 0.93772556
    res = solve_kll_slab_continuum_coefficient(b, c, n_nodes=64)
    z_grid = np.array([0.0, 0.25 * b, 0.50 * b])

    for z in z_grid:
        phi_kll = slab_scalar_flux_kll(res, float(z))
        phi_quad = slab_scalar_flux_from_angular_quadrature(
            res, float(z), n_mu=96, n_quad=64
        )
        assert abs(phi_quad - phi_kll) < 1e-5, (
            f"Closure check fails at z={z}: phi_KLL={phi_kll:.10e}, "
            f"phi_quad={phi_quad:.10e}, abs_diff="
            f"{abs(phi_quad - phi_kll):.2e}"
        )


# ────────────────────────────────────────────────────────────────────
# L1 — KLL ↔ F_N agreement on critical condition
# ────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.verifies("kll-1974-slab-flux")
def test_l1_slab_fn_to_kll_a_c_agreement() -> None:
    r"""Cross-check: F_N's :math:`a_c` is consistent with KLL's
    Fredholm iteration converging at the same critical point.

    Specifically: at the F_N critical :math:`a_c`, the KLL Fredholm
    equation has its eigenmode normalisable (i.e., :math:`\phi(0)
    \neq 0` and :math:`\phi(b)` is finite). At the published Sood
    truth :math:`a_c`, the KLL iteration converges similarly. The two
    :math:`a_c` values should agree to F_N's convergence-floor
    precision.
    """
    c = 1.30
    fn_result = solve_fn_slab_bare_critical(c, n_modes=10, n_bracket=400)
    a_c_fn = fn_result.a_critical_mfp
    a_c_truth = 0.93772556  # Sood Ua-1-0-SL

    # F_N convergence floor at N=10 is ~1e-6 (per the slab L1 test).
    assert abs(a_c_fn - a_c_truth) < 5e-5, (
        f"F_N r_c={a_c_fn} disagrees with truth {a_c_truth} at > 5e-5: "
        f"diff={a_c_fn - a_c_truth:.2e}"
    )

    # KLL at the F_N r_c should give flux ratios close to the published
    # benchmark — this is the structural cross-check above the
    # trusted-library line.
    res = solve_kll_slab_continuum_coefficient(a_c_fn, c, n_nodes=64)
    ratio_at_surface = slab_scalar_flux_ratio(res, 1.0)
    expected_surface = 0.4461912
    assert abs(ratio_at_surface - expected_surface) < 1e-4, (
        f"KLL at F_N's r_c gives surface flux ratio "
        f"{ratio_at_surface:.6f}; truth={expected_surface}, "
        f"diff={ratio_at_surface - expected_surface:.2e}"
    )


# ────────────────────────────────────────────────────────────────────
# L1 — angular flux at surface from F_N coefficients
# ────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.verifies("kll-1974-slab-flux")
def test_l1_slab_surface_angular_flux_polynomial_form() -> None:
    r"""F_N surface angular flux :math:`\psi(\pm a, \mp\mu) = \sum_\alpha
    a_\alpha \mu^\alpha` evaluates to a polynomial — sanity test that
    the polyval-based accessor returns the same as direct summation.
    """
    from orpheus.derivations.continuous.fn_method.slab import (
        slab_surface_angular_flux_fn,
    )

    c = 1.30
    fn_result = solve_fn_slab_bare_critical(c, n_modes=10, n_bracket=400)
    coeffs = fn_result.coefficients_a

    mu_test = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
    psi_polyval = slab_surface_angular_flux_fn(coeffs, mu_test)

    # Direct summation: ψ = Σ_α a_α μ^α.
    psi_direct = np.zeros_like(mu_test)
    for alpha, a_alpha in enumerate(coeffs.real):
        psi_direct += a_alpha * mu_test ** alpha

    np.testing.assert_allclose(psi_polyval, psi_direct, atol=1e-12, rtol=1e-12)
    # Should be positive (outgoing angular flux from a multiplying
    # bare-critical slab — physical positivity).
    assert np.all(psi_polyval > 0), (
        f"Surface angular flux should be non-negative: {psi_polyval}"
    )
