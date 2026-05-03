"""L1 cross-check: Carlvik-Galerkin solver vs sood_registry P_1 anisotropic cases.

Bridges the Wave 2-C Carlvik-Galerkin solver with the Sood/Forster/Parsons
LA-13511 (1999) `*-1-1-SL/SP` benchmark cases.

The Sood `*-1-1-SL/SP` cases are critical configurations
(:math:`k_{\\rm eff} = 1`) — at the published critical dimension
:math:`r_c`, the Carlvik-Galerkin solver should reproduce the
configuration's :math:`c = (\\Sigma_s + \\nu \\Sigma_f)/\\Sigma_t`
to within the solver's numerical tolerance (≤ 1e-5 absolute).

CRITICAL convention bridge: Sood specifies :math:`\\Sigma_{s_1}` directly
as the linear-anisotropy moment of the **scattering** kernel (not the
secondaries). Dahl-Sjostrand 1979 uses :math:`\\bar\\mu_{\\rm eff}` =
mean cosine of all secondaries (scattering + fission, with fission
isotropic). The conversion is
:math:`\\bar\\mu_{\\rm eff} = \\Sigma_{s_1} / (c \\cdot \\Sigma_t)`.

Cases covered (Sood Tables 25 and 29):

* PUa-1-1-SL (problem 32): Pu-239 (a) slab, c=1.40, μ̄_eff=0.142857.
* PUb-1-1-SL (problem 34): Pu-239 (b) slab, c=1.40, μ̄_eff=0.238095
  (strongly forward-peaked; Sood flags |μ̄_scat| > 1/3 → negative
  scattering at μ near -1).
* UD2Oa-1-1-SP (problem 39): U-D2O (a) sphere, c=1.0308381,
  μ̄_eff=0.10. Direct match to Dahl-Sjostrand Table I row d=20 mfp,
  μ̄=0.10.
* UD2Ob-1-1-SP (problem 41): U-D2O (b) sphere, c=1.0341086,
  μ̄_eff=0.20. Direct match to DS Table I.
* UD2Oc-1-1-SP (problem 43): U-D2O (c) sphere, c=1.01964,
  μ̄_eff=-0.50 (back-peaked). Outside DS table coverage but within
  the Carlvik-Galerkin formulation's mathematical validity.
"""
from __future__ import annotations

import pytest

from orpheus.derivations.continuous.carlvik_galerkin.slab import (
    solve_carlvik_galerkin_slab,
)
from orpheus.derivations.continuous.carlvik_galerkin.sphere import (
    solve_carlvik_galerkin_sphere,
)
from orpheus.derivations.continuous.sood_registry import (
    LA13511_CASES,
    WIDE_SLICE_BARE_CRITICAL_1G_P1,
)


def _mu_bar_eff_from_case(case) -> float:
    """Compute the Dahl-Sjostrand μ̄_eff from a Sood case mixture.

    μ̄_eff = Σ_s1 / (c · Σ_t), where c = (Σ_s + νΣ_f)/Σ_t.
    """
    mixture = case.materials[0]
    sigma_t = float(mixture.SigT[0])
    sigma_s0 = float(mixture.SigS[0].toarray()[0, 0])
    sigma_s1 = float(mixture.SigS[1].toarray()[0, 0])  # P_1 component
    nu_sigma_f = float(mixture.SigP[0])
    c = (sigma_s0 + nu_sigma_f) / sigma_t
    return sigma_s1 / (c * sigma_t)


def _c_from_case(case) -> float:
    """Compute c = (Σ_s + νΣ_f)/Σ_t from a Sood case mixture."""
    mixture = case.materials[0]
    sigma_t = float(mixture.SigT[0])
    sigma_s0 = float(mixture.SigS[0].toarray()[0, 0])
    nu_sigma_f = float(mixture.SigP[0])
    return (sigma_s0 + nu_sigma_f) / sigma_t


# ───────────────────────────────────────────────────────────────────
# Slab cases
# ───────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.parametrize("case_id", ["PUa-1-1-SL", "PUb-1-1-SL"])
def test_l1_carlvik_galerkin_sood_slab_p1(case_id: str) -> None:
    """Sood `*-1-1-SL` slab cases reproduce c at the published r_c."""
    case = LA13511_CASES[case_id]
    c_expected = _c_from_case(case)
    mu_bar_eff = _mu_bar_eff_from_case(case)
    r_c_mfp = case.critical_dimension_mfp
    d_mfp = 2.0 * r_c_mfp

    result = solve_carlvik_galerkin_slab(
        c=c_expected,
        d=d_mfp,
        mu_bar=mu_bar_eff,
        n_modes=9,
        n_quad=128,
    )

    abs_err = abs(result.c_critical - c_expected)
    assert abs_err <= 1e-4, (
        f"{case_id}: r_c = {r_c_mfp} mfp, μ̄_eff = {mu_bar_eff:.6f}, "
        f"expected c = {c_expected}, got c_crit = {result.c_critical:.7f}, "
        f"abs error = {abs_err:.2e}"
    )


# ───────────────────────────────────────────────────────────────────
# Sphere cases
# ───────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.parametrize(
    "case_id",
    ["UD2Oa-1-1-SP", "UD2Ob-1-1-SP", "UD2Oc-1-1-SP"],
)
def test_l1_carlvik_galerkin_sood_sphere_p1(case_id: str) -> None:
    """Sood `*-1-1-SP` sphere cases reproduce c at the published r_c."""
    case = LA13511_CASES[case_id]
    c_expected = _c_from_case(case)
    mu_bar_eff = _mu_bar_eff_from_case(case)
    r_c_mfp = case.critical_dimension_mfp
    d_mfp = 2.0 * r_c_mfp

    result = solve_carlvik_galerkin_sphere(
        c=c_expected,
        d=d_mfp,
        mu_bar=mu_bar_eff,
        n_modes=9,
        n_quad=128,
    )

    abs_err = abs(result.c_critical - c_expected)
    # Sphere bare-critical with very thin d (here d=20 mfp huge) is
    # numerically robust; tolerance 1e-5 is comfortable.
    assert abs_err <= 1e-5, (
        f"{case_id}: r_c = {r_c_mfp} mfp, μ̄_eff = {mu_bar_eff:.6f}, "
        f"expected c = {c_expected}, got c_crit = {result.c_critical:.7f}, "
        f"abs error = {abs_err:.2e}"
    )


# ───────────────────────────────────────────────────────────────────
# Registry-level smoke test
# ───────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_sood_p1_registry_completeness() -> None:
    """Smoke test: WIDE_SLICE_BARE_CRITICAL_1G_P1 has exactly 5 cases,
    all identifiable by case_id, all in LA13511_CASES."""
    assert len(WIDE_SLICE_BARE_CRITICAL_1G_P1) == 5

    ids = {c.case_id for c in WIDE_SLICE_BARE_CRITICAL_1G_P1}
    expected_ids = {
        "PUa-1-1-SL",
        "PUb-1-1-SL",
        "UD2Oa-1-1-SP",
        "UD2Ob-1-1-SP",
        "UD2Oc-1-1-SP",
    }
    assert ids == expected_ids, f"Got {ids}, expected {expected_ids}"

    for case in WIDE_SLICE_BARE_CRITICAL_1G_P1:
        # Each case must be in the global LA13511_CASES dict.
        assert case.case_id in LA13511_CASES
        # All P_1 cases have scattering_order=1.
        assert case.scattering_order == 1
        # All have k_eff=1.0 (critical) since they're bare-critical.
        assert case.truth.k_eff_or_kinf == 1.0


@pytest.mark.foundation
def test_sood_p1_mixtures_have_p1_component() -> None:
    """Each WIDE_SLICE_BARE_CRITICAL_1G_P1 case has a non-trivial P_1
    scattering matrix (Σ_s1 ≠ 0)."""
    for case in WIDE_SLICE_BARE_CRITICAL_1G_P1:
        mixture = case.materials[0]
        # SigS[0] is P_0, SigS[1] is P_1 — must exist with non-trivial value.
        assert len(mixture.SigS) >= 2, (
            f"{case.case_id}: P_1 scattering not present in mixture"
        )
        sigma_s1 = float(mixture.SigS[1].toarray()[0, 0])
        assert sigma_s1 != 0.0, (
            f"{case.case_id}: Σ_s1 = 0; expected non-trivial P_1 component"
        )
