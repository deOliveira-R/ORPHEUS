r"""L1 verification for the sphere F_N rich-machinery extension.

Sphere analog of :mod:`.test_fn_la13511_slab_flux`. Verification gates:

* **L1 — Sood / KLL Table VII flux ratios at c = 1.30**: at
  :math:`r/R \in \{0.25, 0.50, 0.75, 1.00\}`, the KLL reconstruction
  agrees with the published reference at ≤ 1e-5 absolute. This is
  the same XS as Sood ``Ua-1-0-SP`` (U-235 (a) bare critical sphere).
* **L1 — KLL Table VII sweep**: for :math:`c \in \{1.05, 1.10, 1.20,
  1.30, 1.40, 1.60\}`, the reconstruction agrees with the published
  values at ≤ 5e-5 absolute.
* **L1 — angular-flux closure**: :math:`\int_{-1}^{1} \psi(r, \mu)\,
  d\mu = \phi(r)` self-consistency at ≤ 1e-5 absolute.

References
----------

* Kaper-Lindeman-Leaf 1974, *Nucl. Sci. Eng.* **54**, 94, Table VII.
* Sood, Forster, Parsons 1999, LA-13511 Table 6 (case ``Ua-1-0-SP``;
  no flux-ratio table for this case in LA-13511, but the XS match
  KLL's c=1.30 row exactly).
* Siewert-Thomas 1986, *Nucl. Sci. Eng.* **94**, 264.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.fn_method.sphere import (
    solve_fn_sphere_bare_critical,
    solve_kll_sphere_continuum_coefficient,
    sphere_angular_flux_from_scalar,
    sphere_scalar_flux_from_angular_quadrature,
    sphere_scalar_flux_kll,
    sphere_scalar_flux_ratio,
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
# L1 — Sood Ua-1-0-SP flux ratios from KLL Table VII
# ────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.verifies("kll-1974-sphere-flux")
def test_l1_sphere_flux_ratios_kll_table_VII_at_c_1p30() -> None:
    r"""KLL Table VII flux ratios at c = 1.30 to ≤ 1e-5.

    Same XS as Sood ``Ua-1-0-SP`` (U-235 (a) bare-critical sphere,
    c=1.30, R_c=2.4248249802 mfp).

    +-------+--------------+
    | r/R   | φ(r)/φ(0)    |
    +=======+==============+
    | 0.25  | 0.93244907   |
    | 0.50  | 0.74553332   |
    | 0.75  | 0.48095413   |
    | 1.00  | 0.17177706   |
    +-------+--------------+
    """
    case = LA13511_CASES["Ua-1-0-SP"]
    sigma_t = case.sigma_t[0]
    sigma_s = case.sigma_s[0, 0]
    nu_sigma_f = case.nu_sigma_f[0]
    c = float((sigma_s + nu_sigma_f) / sigma_t)
    assert abs(c - 1.30) < 1e-10, f"Expected c=1.30, got c={c}"

    R_kll = case.critical_dimension_mfp
    res = solve_kll_sphere_continuum_coefficient(R_kll, c, n_nodes=64)
    assert res.converged, "KLL Fredholm iteration did not converge"

    r_over_R = np.array([0.25, 0.50, 0.75, 1.00])
    computed = sphere_scalar_flux_ratio(res, r_over_R)
    # The truth values populated in sood_registry should match KLL Table VII.
    flux_ratios_in_registry = case.truth.flux_ratios
    assert flux_ratios_in_registry is not None, (
        "Ua-1-0-SP truth.flux_ratios must be populated"
    )
    expected = np.array([
        flux_ratios_in_registry[0.25],
        flux_ratios_in_registry[0.50],
        flux_ratios_in_registry[0.75],
        flux_ratios_in_registry[1.00],
    ])
    abs_err = np.abs(computed - expected)

    assert np.max(abs_err) < 1e-5, (
        f"KLL sphere flux ratios at c=1.30 disagree with KLL Table VII "
        f"at > 1e-5: computed={computed}, expected={expected}, "
        f"abs_err={abs_err}"
    )


@pytest.mark.l1
@pytest.mark.verifies("kll-1974-sphere-flux")
@pytest.mark.parametrize(
    "c,R_truth,truth_ratios",
    [
        (1.05, 7.27718_17945, (0.91612699, 0.68954766, 0.38621118, 0.07449726)),
        (1.10, 4.87271_42665, (0.92169358, 0.70849371, 0.41777283, 0.10487685)),
        (1.20, 3.17207_25129, (0.92831440, 0.73121497, 0.45624613, 0.14462575)),
        (1.30, 2.42482_49802, (0.93244907, 0.74553332, 0.48095413, 0.17177706)),
        (1.40, 1.98534_34324, (0.93538006, 0.75575352, 0.49884364, 0.19222603)),
        (1.60, 1.47609_85891, (0.93935552, 0.76971786, 0.52365345, 0.22166498)),
    ],
)
def test_l1_sphere_flux_ratios_kll_table_VII_sweep(
    c: float, R_truth: float, truth_ratios: tuple[float, ...]
) -> None:
    r"""KLL Table VII flux ratios sweep over :math:`c \in [1.05, 1.60]`."""
    res = solve_kll_sphere_continuum_coefficient(
        R_truth, c, n_nodes=80, max_iter=400
    )
    assert res.converged, f"KLL did not converge for c={c}, R={R_truth}"

    r_over_R = np.array([0.25, 0.50, 0.75, 1.00])
    computed = sphere_scalar_flux_ratio(res, r_over_R)
    expected = np.array(truth_ratios)
    abs_err = np.abs(computed - expected)

    if c <= 1.20:
        tol = 1e-7
    elif c <= 1.40:
        tol = 1e-7
    else:
        tol = 5e-7

    assert np.max(abs_err) < tol, (
        f"KLL sphere flux ratios at c={c} disagree with Table VII at > "
        f"{tol}: computed={computed}, expected={expected}, "
        f"abs_err={abs_err}"
    )


# ────────────────────────────────────────────────────────────────────
# L1 — angular-flux closure
# ────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.verifies("kll-1974-sphere-flux")
def test_l1_sphere_angular_flux_closure() -> None:
    r"""Closure: :math:`\phi(r) = \int_{-1}^{1} \psi(r, \mu)\,d\mu`.

    Tests that the angular flux reconstructed from the converged
    :math:`\phi(r)` integrates back to :math:`\phi(r)` to within the
    quadrature accuracy.
    """
    c = 1.30
    R = 2.4248249802
    res = solve_kll_sphere_continuum_coefficient(R, c, n_nodes=64)
    r_grid = np.array([0.25 * R, 0.50 * R, 0.75 * R])

    for r in r_grid:
        phi_kll = sphere_scalar_flux_kll(res, float(r))
        phi_quad = sphere_scalar_flux_from_angular_quadrature(
            res, float(r), n_mu=96, n_quad=64
        )
        assert abs(phi_quad - phi_kll) < 1e-5, (
            f"Closure check fails at r={r}: phi_KLL={phi_kll:.10e}, "
            f"phi_quad={phi_quad:.10e}, abs_diff="
            f"{abs(phi_quad - phi_kll):.2e}"
        )


# ────────────────────────────────────────────────────────────────────
# L1 — KLL ↔ F_N agreement on critical condition
# ────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.verifies("kll-1974-sphere-flux")
def test_l1_sphere_fn_to_kll_R_c_agreement() -> None:
    r"""Cross-check: F_N's :math:`R_c` is consistent with KLL's
    Fredholm iteration."""
    c = 1.30
    fn_result = solve_fn_sphere_bare_critical(c, n_modes=10, n_bracket=800)
    R_c_fn = fn_result.R_critical_mfp
    R_c_truth = 2.4248249802

    assert abs(R_c_fn - R_c_truth) < 5e-5, (
        f"F_N R_c={R_c_fn} disagrees with truth {R_c_truth} at > 5e-5: "
        f"diff={R_c_fn - R_c_truth:.2e}"
    )

    res = solve_kll_sphere_continuum_coefficient(R_c_fn, c, n_nodes=64)
    ratio_at_surface = sphere_scalar_flux_ratio(res, 1.0)
    expected_surface = 0.17177706
    assert abs(ratio_at_surface - expected_surface) < 1e-4, (
        f"KLL at F_N's R_c gives surface flux ratio "
        f"{ratio_at_surface:.6f}; truth={expected_surface}, "
        f"diff={ratio_at_surface - expected_surface:.2e}"
    )


# ────────────────────────────────────────────────────────────────────
# L1 — angular flux at surface from F_N coefficients
# ────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.verifies("kll-1974-sphere-flux")
def test_l1_sphere_surface_angular_flux_polynomial_form() -> None:
    r"""F_N surface angular flux is a polynomial in :math:`\mu`."""
    from orpheus.derivations.continuous.fn_method.sphere import (
        sphere_surface_angular_flux_fn,
    )

    c = 1.30
    fn_result = solve_fn_sphere_bare_critical(c, n_modes=10, n_bracket=800)
    coeffs = fn_result.coefficients_a

    mu_test = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
    psi_polyval = sphere_surface_angular_flux_fn(coeffs, mu_test)

    psi_direct = np.zeros_like(mu_test)
    for alpha, a_alpha in enumerate(coeffs.real):
        psi_direct += a_alpha * mu_test ** alpha

    np.testing.assert_allclose(psi_polyval, psi_direct, atol=1e-12, rtol=1e-12)
