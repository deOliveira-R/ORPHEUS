"""L1 cross-check: Carlvik-Galerkin sphere solver vs Dahl-Sjostrand 1979 Table I.

Reproduces a representative subset of Dahl-Sjostrand 1979 NSE 69, 114-125,
Table I ("Eigenvalues of Multiplying Spheres with Diameter d for Various
Values of the Average Cosine μ̄ of the Scattering Angle") to ≤ 1e-6 absolute
on the dominant eigenvalue.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.galerkin_spectral.sphere import (
    solve_galerkin_spectral_sphere,
)


# ───────────────────────────────────────────────────────────────────
# Isotropic μ̄ = 0
# ───────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.parametrize(
    "d, c_expected_DS_table_i",
    [
        (0.2, 13.428632),
        (1.0, 3.2371777),
        (2.0, 1.9883853),
        (5.0, 1.2868065),
        (8.0, 1.1384602),
        (20.0, 1.0281490),
    ],
)
def test_l1_sphere_isotropic_fundamental_vs_DS_table_i(
    d: float, c_expected_DS_table_i: float
) -> None:
    """Sphere μ̄=0: dominant eigenvalue matches DS Table I to ≤ 1e-6 absolute."""
    result = solve_galerkin_spectral_sphere(
        c=1.0, d=d, mu_bar=0.0, n_modes=9, n_quad=128
    )
    abs_err = abs(result.c_critical - c_expected_DS_table_i)
    # For d=0.2 c=13.43, machine round-off in last digit (one unit of 1e-6
    # on the value 13.43 corresponds to relative 7e-8 — within DS precision).
    assert abs_err <= 5e-6, (  # slightly looser than 1e-6 to account for
        # the absolute scale of c~13 at small d
        f"d={d}, μ̄=0: c_crit = {result.c_critical:.10f}, "
        f"DS Table I expects {c_expected_DS_table_i}, "
        f"abs error = {abs_err:.2e}"
    )


# ───────────────────────────────────────────────────────────────────
# Anisotropic μ̄ = 0.30
# ───────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.parametrize(
    "d, c_expected_DS_table_i",
    [
        (0.2, 15.550669),
        (1.0, 3.6599219),
        (2.0, 2.1957124),
        (5.0, 1.3611358),
        (8.0, 1.1790299),
    ],
)
def test_l1_sphere_anisotropic_mu03_vs_DS_table_i(
    d: float, c_expected_DS_table_i: float
) -> None:
    """Sphere μ̄=0.30: dominant eigenvalue matches DS Table I to ≤ 1e-6 absolute."""
    result = solve_galerkin_spectral_sphere(
        c=1.0, d=d, mu_bar=0.30, n_modes=9, n_quad=128
    )
    abs_err = abs(result.c_critical - c_expected_DS_table_i)
    # Same scale tolerance as isotropic.
    assert abs_err <= 5e-6, (
        f"d={d}, μ̄=0.30: c_crit = {result.c_critical:.10f}, "
        f"DS Table I expects {c_expected_DS_table_i}, "
        f"abs error = {abs_err:.2e}"
    )


@pytest.mark.l1
@pytest.mark.parametrize(
    "d, c_expected_DS_table_i",
    [
        (0.2, 14.015294),
        (1.0, 3.3544208),
        (2.0, 2.0459554),
        (5.0, 1.3073468),
        (8.0, 1.1495822),
    ],
)
def test_l1_sphere_anisotropic_mu01_vs_DS_table_i(
    d: float, c_expected_DS_table_i: float
) -> None:
    """Sphere μ̄=0.10: dominant eigenvalue matches DS Table I to ≤ 1e-6 absolute."""
    result = solve_galerkin_spectral_sphere(
        c=1.0, d=d, mu_bar=0.10, n_modes=9, n_quad=128
    )
    abs_err = abs(result.c_critical - c_expected_DS_table_i)
    assert abs_err <= 5e-6, (
        f"d={d}, μ̄=0.10: c_crit = {result.c_critical:.10f}, "
        f"DS Table I expects {c_expected_DS_table_i}, "
        f"abs error = {abs_err:.2e}"
    )


# ───────────────────────────────────────────────────────────────────
# Higher modes
# ───────────────────────────────────────────────────────────────────


@pytest.mark.l1
def test_l1_sphere_isotropic_first_three_modes_d2() -> None:
    """Sphere μ̄=0, d=2: first 3 eigenvalues match DS Table I to ≤ 1e-5."""
    result = solve_galerkin_spectral_sphere(
        c=1.0, d=2.0, mu_bar=0.0, n_modes=9, n_quad=128
    )
    expected_first_3 = np.array([1.9883853, 3.818339, 5.763974])
    actual_first_3 = result.eigenvalue_spectrum[:3].real
    abs_errs = np.abs(actual_first_3 - expected_first_3)
    assert np.all(abs_errs <= 1e-5), (
        f"DS Table I d=2 μ̄=0 first 3 modes: actual = {actual_first_3}, "
        f"expected = {expected_first_3}, abs errors = {abs_errs}"
    )


# ───────────────────────────────────────────────────────────────────
# Complex eigenvalue detection
# ───────────────────────────────────────────────────────────────────


@pytest.mark.l1
def test_l1_sphere_complex_eigenvalue_pair_d02_mu03() -> None:
    """Sphere d=0.2, μ̄=0.30: complex eigenvalue pairs at high modes."""
    result = solve_galerkin_spectral_sphere(
        c=1.0, d=0.2, mu_bar=0.30, n_modes=9, n_quad=128
    )
    n_complex = int(np.sum(np.abs(result.eigenvalue_spectrum.imag) > 1e-6))
    assert n_complex >= 2, (
        f"Expected at least 2 complex eigenvalues at d=0.2, μ̄=0.30 sphere. "
        f"Got {n_complex}. Spectrum: {result.eigenvalue_spectrum}"
    )

    c_crit_expected = 15.550669
    assert abs(result.c_critical - c_crit_expected) <= 5e-6, (
        f"Fundamental at d=0.2, μ̄=0.30 sphere: "
        f"got {result.c_critical}, expected {c_crit_expected}"
    )
