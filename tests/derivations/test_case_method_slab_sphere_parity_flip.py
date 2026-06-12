r"""Verify Atalay 1997 sphere/slab parity-flip equivalence numerically.

The structural identity proved symbolically in
:func:`...origins.derivations.derive_atalay_critical_sphere_eq54_via_parity_flip`
(V_case.6) is also enforced numerically here:

* Slab and sphere share the SAME core machinery
  (:mod:`...core.dispersion`, :mod:`...core.x_function`,
  :mod:`...core.extrapolated_endpoint`).
* The K_j (slab Eq 38) and L_j (sphere Eq 51) differ ONLY through
  T(R, μ) vs T_1(R, μ).
* For vacuum (R = 0), T(0, μ) = T_1(0, μ) = e^{-2d/μ}, hence
  K_j(R=0, d) = L_j(R=0, d) **bit-for-bit**.
"""
from __future__ import annotations

import math

import pytest

from orpheus.derivations.continuous.singular_eigenfunction.core.dispersion import (
    case_atalay_u0,
    nu_bar_atalay,
)
from orpheus.derivations.continuous.singular_eigenfunction.core.half_range import (
    T1_function,
    T_function,
    atalay_K_moments,
    atalay_L_moments,
    clear_X_cache,
)


@pytest.mark.l1
@pytest.mark.verifies(
    "singular-eigenfunction-eq46", "singular-eigenfunction-eq54",
)
def test_K_equals_L_at_vacuum_BC():
    """At R = 0 (vacuum), slab K_j == sphere L_j to machine precision.

    This is a numerical *consequence* of the parity-flip structural
    identity: T(0, μ) = T_1(0, μ) = e^{-2d/μ} — the only difference
    between the two integrands vanishes at vacuum BC.
    """
    clear_X_cache()
    c, R, f1 = 1.30, 0.0, 0.0
    res = case_atalay_u0(c=c, f1=f1)
    u0 = res.u0
    nu_bar = nu_bar_atalay(c=c, f1=f1, u0=u0)

    Ks = atalay_K_moments(c=c, f1=f1, u0=u0, nu_bar=nu_bar, R=R, d_thick=1.0)
    Ls = atalay_L_moments(c=c, f1=f1, u0=u0, nu_bar=nu_bar, R=R, d_thick=1.0)

    for j in (0, 1, 2):
        assert math.isclose(Ks[j], Ls[j], rel_tol=1e-13), (
            f"K_{j} != L_{j} at R=0 (vacuum): K = {Ks[j]}, L = {Ls[j]}. "
            f"Parity-flip equivalence broken — T(R=0,μ) and T_1(R=0,μ) "
            f"should both equal e^{{-2d/μ}}."
        )


@pytest.mark.l1
@pytest.mark.verifies(
    "singular-eigenfunction-eq46", "singular-eigenfunction-eq54",
)
def test_T_functions_at_vacuum():
    """T(R=0, μ) = T_1(R=0, μ) = e^{-2d/μ} (both reduce to the same vacuum form)."""
    for mu in (0.1, 0.3, 0.5, 0.7, 0.9):
        for d_thick in (0.5, 1.0, 2.0, 5.0):
            T_slab = T_function(0.0, mu, d_thick)
            T_sphere = T1_function(0.0, mu, d_thick)
            expected = math.exp(-2.0 * d_thick / mu)
            assert math.isclose(T_slab, T_sphere, rel_tol=1e-15)
            assert math.isclose(T_slab, expected, rel_tol=1e-12)


@pytest.mark.l1
@pytest.mark.verifies(
    "singular-eigenfunction-eq46", "singular-eigenfunction-eq54",
)
def test_T_functions_at_perfect_reflector():
    """T(R=1, μ) = -1, T_1(R=1, μ) = +1 (the parity-flip surfaces here)."""
    for mu in (0.1, 0.3, 0.7):
        for d_thick in (0.5, 1.0, 5.0):
            T_slab = T_function(1.0, mu, d_thick)
            T_sphere = T1_function(1.0, mu, d_thick)
            assert math.isclose(T_slab, -1.0, abs_tol=1e-12)
            assert math.isclose(T_sphere, 1.0, abs_tol=1e-12)


@pytest.mark.l1
@pytest.mark.verifies(
    "singular-eigenfunction-eq46", "singular-eigenfunction-eq54",
)
def test_T_functions_partial_reflection_signs():
    """For 0 < R < 1, T(R,μ) ∈ (-1, e^{-2d/μ}) and T_1(R,μ) ∈ (e^{-2d/μ}, 1).

    The parity flip changes the signed range of T entirely.
    """
    # At R=0.5, large d: T(R=0.5, μ=0.5, d=2) = (0.5 - e^-8)/(0.5 e^-8 - 1) ≈ 0.5/-1 ≈ -0.5
    T_slab_05 = T_function(0.5, 0.5, 2.0)
    T_sphere_05 = T1_function(0.5, 0.5, 2.0)
    expmf = math.exp(-2.0 * 2.0 / 0.5)
    assert -1.0 < T_slab_05 < expmf
    assert expmf < T_sphere_05 < 1.0
