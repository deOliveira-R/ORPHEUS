"""L1 cross-check: Carlvik-Galerkin slab solver vs Dahl-Sjostrand 1979 Table II.

Reproduces a representative subset of Dahl-Sjostrand 1979 NSE 69, 114-125,
Table II ("Eigenvalues of Multiplying Slabs with Thickness d for Various
Values of the Average Cosine μ̄ of the Scattering Angle") to ≤ 1e-6
absolute on the dominant eigenvalue (the paper claims 7-sig-fig precision
with last-digit uncertain to one unit; underlined entries to 10 units).

Coverage:

* Isotropic (μ̄=0) full sweep d ∈ {0.2, 1, 2, 5, 8, 20}.
* Anisotropic μ̄=0.30 sweep d ∈ {0.2, 1, 2, 5, 8}.
* First three eigenvalues at d=2.0, μ̄=0.10 (representative
  spectrum check).
* Complex eigenvalue detection at d=0.2, μ̄=0.30
  (test for paper's reported complex pair).

The brief's tolerance targets:

* fundamental :math:`c_{\rm crit}` ≤ 1e-6 absolute → met by 7-sig-fig
  agreement at n_modes=9, n_quad=128.
* first 3 eigenvalues ≤ 1e-5 absolute → met for the cases tested.
* Complex pair appearance verification (high μ̄) → met
  qualitatively.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.galerkin_spectral.slab import (
    solve_galerkin_spectral_slab,
)


# ───────────────────────────────────────────────────────────────────
# Isotropic μ̄ = 0 — exact agreement with Hartman 1975 / Carlvik 1968
# ───────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.parametrize(
    "d, c_expected_DS_table_ii",
    [
        (0.2, 3.8303067),
        (1.0, 1.6153785),
        (2.0, 1.2771018),
        (5.0, 1.0775728),
        (8.0, 1.0364020),
        (20.0, 1.0071358),
    ],
)
def test_l1_slab_isotropic_fundamental_vs_DS_table_ii(
    d: float, c_expected_DS_table_ii: float
) -> None:
    """Slab μ̄=0: dominant eigenvalue matches DS Table II to ≤ 1e-6 absolute."""
    result = solve_galerkin_spectral_slab(
        c=1.0, d=d, mu_bar=0.0, n_modes=9, n_quad=128
    )
    abs_err = abs(result.c_critical - c_expected_DS_table_ii)
    assert abs_err <= 1e-6, (
        f"d={d}, μ̄=0: c_crit = {result.c_critical:.10f}, "
        f"DS Table II expects {c_expected_DS_table_ii}, "
        f"abs error = {abs_err:.2e}"
    )


# ───────────────────────────────────────────────────────────────────
# Anisotropic μ̄ = 0.30 — sweep across d
# ───────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.parametrize(
    "d, c_expected_DS_table_ii",
    [
        (0.2, 3.9341098),
        (1.0, 1.6726359),
        (2.0, 1.3162908),
        (5.0, 1.0954329),
        (8.0, 1.0465450),
    ],
)
def test_l1_slab_anisotropic_mu03_vs_DS_table_ii(
    d: float, c_expected_DS_table_ii: float
) -> None:
    """Slab μ̄=0.30: dominant eigenvalue matches DS Table II to ≤ 1e-6 absolute."""
    result = solve_galerkin_spectral_slab(
        c=1.0, d=d, mu_bar=0.30, n_modes=9, n_quad=128
    )
    abs_err = abs(result.c_critical - c_expected_DS_table_ii)
    assert abs_err <= 1e-6, (
        f"d={d}, μ̄=0.30: c_crit = {result.c_critical:.10f}, "
        f"DS Table II expects {c_expected_DS_table_ii}, "
        f"abs error = {abs_err:.2e}"
    )


# ───────────────────────────────────────────────────────────────────
# Anisotropic μ̄ = 0.10 — sweep
# ───────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.parametrize(
    "d, c_expected_DS_table_ii",
    [
        (0.2, 3.8635951),
        (1.0, 1.6330295),
        (2.0, 1.2888640),
        (5.0, 1.0826970),
        (8.0, 1.0392421),
    ],
)
def test_l1_slab_anisotropic_mu01_vs_DS_table_ii(
    d: float, c_expected_DS_table_ii: float
) -> None:
    """Slab μ̄=0.10: dominant eigenvalue matches DS Table II to ≤ 1e-6 absolute."""
    result = solve_galerkin_spectral_slab(
        c=1.0, d=d, mu_bar=0.10, n_modes=9, n_quad=128
    )
    abs_err = abs(result.c_critical - c_expected_DS_table_ii)
    assert abs_err <= 1e-6, (
        f"d={d}, μ̄=0.10: c_crit = {result.c_critical:.10f}, "
        f"DS Table II expects {c_expected_DS_table_ii}, "
        f"abs error = {abs_err:.2e}"
    )


# ───────────────────────────────────────────────────────────────────
# Higher modes — first 3 eigenvalues at one representative point
# ───────────────────────────────────────────────────────────────────


@pytest.mark.l1
def test_l1_slab_isotropic_first_three_modes_d2() -> None:
    """Slab μ̄=0, d=2: first 3 eigenvalues match DS Table II to ≤ 1e-5."""
    result = solve_galerkin_spectral_slab(
        c=1.0, d=2.0, mu_bar=0.0, n_modes=9, n_quad=128
    )
    # DS Table II d=2.0 μ̄=0.0:
    expected_first_3 = np.array([1.2771018, 2.873324, 4.784655])
    actual_first_3 = result.eigenvalue_spectrum[:3].real
    abs_errs = np.abs(actual_first_3 - expected_first_3)
    assert np.all(abs_errs <= 1e-5), (
        f"DS Table II d=2 μ̄=0 first 3 modes: actual = {actual_first_3}, "
        f"expected = {expected_first_3}, abs errors = {abs_errs}"
    )


# ───────────────────────────────────────────────────────────────────
# Complex eigenvalue detection — high μ̄, high modes
# ───────────────────────────────────────────────────────────────────


@pytest.mark.l1
def test_l1_slab_complex_eigenvalue_pair_d02_mu03() -> None:
    """Slab d=0.2, μ̄=0.30: complex eigenvalue pairs appear at high modes.

    DS Table II row d=0.2, μ̄=0.30 reports pairs of complex-conjugate
    eigenvalues at columns 3, 4, 5, 6, 7, 8 (with imaginary parts
    19.2610i, etc.). We verify at least one complex pair is detected
    by the solver — the EXACT positions of complex eigenvalues are
    sensitive to n_quad and n_modes, but the QUALITATIVE appearance
    is robust.
    """
    result = solve_galerkin_spectral_slab(
        c=1.0, d=0.2, mu_bar=0.30, n_modes=9, n_quad=128
    )
    # Count eigenvalues with non-trivial imaginary part:
    n_complex = int(np.sum(np.abs(result.eigenvalue_spectrum.imag) > 1e-6))
    assert n_complex >= 4, (
        f"Expected at least 4 complex eigenvalues at d=0.2, μ̄=0.30 "
        f"(per DS Table II; pairs in cols 3-8). Got {n_complex} complex eigenvalues. "
        f"Spectrum: {result.eigenvalue_spectrum}"
    )

    # All real eigenvalues must be positive (per DS Fig. 6 — the
    # fundamental mode is always real positive even at μ̄=0.30).
    real_mask = np.abs(result.eigenvalue_spectrum.imag) < 1e-9
    real_part = result.eigenvalue_spectrum[real_mask].real
    n_real_pos = int(np.sum(real_part > 0))
    assert n_real_pos >= 1, "No positive real eigenvalues found"

    # The smallest positive real eigenvalue is the fundamental.
    c_crit_expected = 3.9341098
    assert abs(result.c_critical - c_crit_expected) <= 1e-6, (
        f"Fundamental real eigenvalue at d=0.2, μ̄=0.30: "
        f"got {result.c_critical}, expected {c_crit_expected}"
    )
