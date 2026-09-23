r"""Foundation tests for the Atalay 1997 Case singular-eigenfunction method.

One ``@pytest.mark.foundation`` test per ``derive_*()`` function in
:mod:`orpheus.derivations.continuous.singular_eigenfunction.origins.slab_sphere_derivations`.
These pin the **algebra-of-record**: the published Atalay equations
follow from the reduction chain documented in the paper.

The test count equals the V_n claim count (per algebra-of-record skill).

References
----------

* Atalay, M.A. (1997), *Prog. Nucl. Energy* **31**(3), 229-252.
"""
from __future__ import annotations

import pytest

from orpheus.derivations.continuous.singular_eigenfunction.origins.slab_sphere_derivations import (
    derive_atalay_critical_slab_eq46,
    derive_atalay_critical_sphere_eq54_via_parity_flip,
    derive_atalay_dispersion_linear_anisotropic,
    derive_atalay_extrapolated_endpoint_eq42,
    derive_atalay_fredholm_form_eq27_to_eq32,
    derive_atalay_half_range_eqs28_to_31,
    derive_atalay_symmetry_conditions_eq13_14_47_to_49,
    derive_atalay_validity_bound_eq5,
    derive_atalay_x_function_eq40,
)


@pytest.mark.foundation
def test_v_case_1_dispersion_linear_anisotropic():
    """V_case.1: Atalay Eq 11 ↔ Eq 12 dispersion quadratic in c."""
    r = derive_atalay_dispersion_linear_anisotropic()
    assert r["pass"], f"V_case.1 FAIL: {r}"


@pytest.mark.foundation
def test_v_case_2_symmetry_conditions():
    """V_case.2: slab Eqs 13-14 (s=+1) and sphere Eqs 48-49 (s=-1) parity."""
    r = derive_atalay_symmetry_conditions_eq13_14_47_to_49()
    assert r["pass"], f"V_case.2 FAIL: {r}"


@pytest.mark.foundation
def test_v_case_3_half_range_eqs28_to_31():
    """V_case.3: Atalay Eqs 28-31 share Wiener-Hopf X-function structure."""
    r = derive_atalay_half_range_eqs28_to_31()
    assert r["pass"], f"V_case.3 FAIL: {r}"


@pytest.mark.foundation
def test_v_case_4_fredholm_form_eq27_to_eq32():
    """V_case.4: Atalay Eq 32 first-order Fredholm form (prefactor + T limits)."""
    r = derive_atalay_fredholm_form_eq27_to_eq32()
    assert r["pass"], f"V_case.4 FAIL: {r}"


@pytest.mark.foundation
def test_v_case_5_critical_slab_eq46():
    """V_case.5: Atalay Eq 43 → Eq 46 via Re/Im of complex-conjugate quotient."""
    r = derive_atalay_critical_slab_eq46()
    assert r["pass"], f"V_case.5 FAIL: {r}"


@pytest.mark.foundation
def test_v_case_6_critical_sphere_eq54_via_parity_flip():
    """V_case.6: Atalay Eq 54 (sphere) ↔ Eq 46 (slab) via 2nd-term sign flip in z_L."""
    r = derive_atalay_critical_sphere_eq54_via_parity_flip()
    assert r["pass"], f"V_case.6 FAIL: {r}"


@pytest.mark.foundation
def test_v_case_7_extrapolated_endpoint_eq42():
    """V_case.7: Atalay Eq 42 z_0 integrand bracket matches Eq 40 X bracket."""
    r = derive_atalay_extrapolated_endpoint_eq42()
    assert r["pass"], f"V_case.7 FAIL: {r}"


@pytest.mark.foundation
def test_v_case_8_x_function_eq40():
    """V_case.8: Atalay Eq 40 X-function integrand structure + isotropic limit."""
    r = derive_atalay_x_function_eq40()
    assert r["pass"], f"V_case.8 FAIL: {r}"


@pytest.mark.foundation
def test_v_case_9_validity_bound_eq5():
    """V_case.9: Atalay Eq 5 validity bound c ≤ 1 + 1/(3 f_1)."""
    r = derive_atalay_validity_bound_eq5()
    assert r["pass"], f"V_case.9 FAIL: {r}"
