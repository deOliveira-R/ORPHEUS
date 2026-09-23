"""Branch-1 SymPy foundation tests for the Carlvik-Galerkin method.

Each test pins one ``derive_*()`` claim from
:mod:`orpheus.derivations.continuous.galerkin_spectral.origins.derivations`.
The full V_cg.N verification list:

* V_cg.1 — Galerkin LHS = 2 F_m via Legendre orthogonality.
* V_cg.2 — Eq. (3) matrix eigenvalue form structure.
* V_cg.3 — A_{0,0}(a) closed form (slab even basis, SymPy ↔ hand-derived).
* V_cg.4 — Even / odd Legendre basis parity (slab vs sphere).
* V_cg.5 — B_{m,n} boundary-chord rank-1 structure.
* V_cg.6 — Eq. (4) block-matrix linearization equivalent to Eq. (3).
* V_cg.7 — Carlvik 1968 Eq. (4b) sign correction documented.
* V_cg.8 — μ̄ = 0 isotropic limit reduces to Carlvik 1968 isotropic eq.

These are the algebra-of-record bedrock for the Branch-2 production
solvers in :mod:`...galerkin_spectral.slab` and :mod:`...galerkin_spectral.sphere`.
"""
from __future__ import annotations

import pytest

from orpheus.derivations.continuous.galerkin_spectral.origins.derivations import (
    derive_B_mn_boundary_chord_form,
    derive_carlvik_eq4b_corrected_form,
    derive_eq3_matrix_eigenvalue,
    derive_eq4_block_linearization,
    derive_galerkin_lhs_identity,
    derive_isotropic_limit,
    derive_low_order_A_mn_slab,
    derive_low_order_A_mn_sphere,
)


@pytest.mark.foundation
def test_v_cg_1_galerkin_lhs_identity() -> None:
    """V_cg.1: Galerkin LHS = 2 F_m via Legendre orthogonality."""
    r = derive_galerkin_lhs_identity()
    assert r["pass"], f"V_cg.1 failed: {r}"


@pytest.mark.foundation
def test_v_cg_2_eq3_matrix_eigenvalue_form() -> None:
    """V_cg.2: Eq. (3) matrix eigenvalue structure A - 3μ̄(c-1)B = 1/(cd)."""
    r = derive_eq3_matrix_eigenvalue()
    assert r["pass"], f"V_cg.2 failed: {r}"


@pytest.mark.foundation
def test_v_cg_3_A_00_closed_form_slab() -> None:
    """V_cg.3: A_{0,0}(a) closed form, SymPy and hand-derived agree.

    NOTE: this test takes ~5-7s because SymPy must integrate
    E_1(a|x-y|) symbolically. Within budget for a foundation test
    but slowest in the V_cg suite.
    """
    r = derive_low_order_A_mn_slab()
    assert r["pass"], f"V_cg.3 failed: {r}"


@pytest.mark.foundation
def test_v_cg_4_basis_parity_slab_vs_sphere() -> None:
    """V_cg.4: Even-Legendre basis sums to even-x function (slab);
    odd-Legendre basis sums to odd-x function (sphere)."""
    r = derive_low_order_A_mn_sphere()
    assert r["pass"], f"V_cg.4 failed: {r}"


@pytest.mark.foundation
def test_v_cg_5_B_mn_boundary_chord_rank_one() -> None:
    """V_cg.5: B_{m,n} boundary-chord term has rank 1 for both slab (q=0)
    and sphere (q=1)."""
    r = derive_B_mn_boundary_chord_form()
    assert r["pass"], f"V_cg.5 failed: {r}"


@pytest.mark.foundation
def test_v_cg_6_eq4_block_linearization() -> None:
    """V_cg.6: Eq. (4) block-matrix form is algebraically equivalent to Eq. (3)."""
    r = derive_eq4_block_linearization()
    assert r["pass"], f"V_cg.6 failed: {r}"


@pytest.mark.foundation
def test_v_cg_7_carlvik_eq4b_correction_documented() -> None:
    """V_cg.7: Documents the Carlvik 1968 Eq. (4b) sign correction."""
    r = derive_carlvik_eq4b_corrected_form()
    assert r["pass"], f"V_cg.7 failed: {r}"


@pytest.mark.foundation
def test_v_cg_8_isotropic_limit() -> None:
    """V_cg.8: μ̄ = 0 reduces Eq. (3) to Carlvik 1968 isotropic
    eigenvalue equation A F = (1/(cd)) F."""
    r = derive_isotropic_limit()
    assert r["pass"], f"V_cg.8 failed: {r}"
