r"""Foundation symbolic gates for the F_N flux-reconstruction
extension (slab + sphere shared identities).

One ``@pytest.mark.foundation`` test per ``derive_*()`` function in
:mod:`orpheus.derivations.continuous.fn_method.origins.fn_flux_reconstruction_derivations`.
These verify the algebraic skeleton of the rich-machinery extension —
KLL Eqs. 7 and 15 structure, characteristic-integration BTE
satisfaction, and the universal angular-flux closure
:math:`\phi = \int_{-1}^{1} \psi\,d\mu`.

References
----------

* Kaper-Lindeman-Leaf 1974, Eqs. 5-7 (slab), 13-15 (sphere).
"""
from __future__ import annotations

import pytest

from orpheus.derivations.continuous.fn_method.origins.fn_flux_reconstruction_derivations import (
    derive_scalar_flux_angular_integral,
    derive_slab_kll_phi_eq7_structure,
    derive_slab_phi_endpoint_normalization,
    derive_slab_psi_from_phi_characteristic,
    derive_sphere_kll_phi_eq15_structure,
    derive_sphere_psi_from_phi_characteristic,
)


@pytest.mark.foundation
def test_v_fn_flux_slab_1_kll_eq7_structure() -> None:
    """V_fn-flux-slab.1 — KLL Eq. 7 cosine + cosh-integral structure."""
    result = derive_slab_kll_phi_eq7_structure()
    assert result["pass"], f"V_fn-flux-slab.1 failed: {result['name']}"


@pytest.mark.foundation
def test_v_fn_flux_slab_2_endpoint_normalization() -> None:
    """V_fn-flux-slab.2 — :math:`\\phi(z)/\\phi(0)` is normalisation-free."""
    result = derive_slab_phi_endpoint_normalization()
    assert result["pass"], f"V_fn-flux-slab.2 failed: {result['name']}"


@pytest.mark.foundation
def test_v_fn_flux_slab_3_psi_from_phi_characteristic() -> None:
    """V_fn-flux-slab.3 — interior :math:`\\psi(z, \\mu)` via characteristic
    integration satisfies the BTE + vacuum BC."""
    result = derive_slab_psi_from_phi_characteristic()
    assert result["pass"], f"V_fn-flux-slab.3 failed: {result['name']}"


@pytest.mark.foundation
def test_v_fn_flux_sphere_1_kll_eq15_structure() -> None:
    """V_fn-flux-sphere.1 — KLL Eq. 15 sin-c + sinh-integral structure."""
    result = derive_sphere_kll_phi_eq15_structure()
    assert result["pass"], f"V_fn-flux-sphere.1 failed: {result['name']}"


@pytest.mark.foundation
def test_v_fn_flux_sphere_2_psi_from_phi_characteristic() -> None:
    """V_fn-flux-sphere.2 — sphere chord-length :math:`s_{\\rm in}` formula
    via characteristic geometry."""
    result = derive_sphere_psi_from_phi_characteristic()
    assert result["pass"], f"V_fn-flux-sphere.2 failed: {result['name']}"


@pytest.mark.foundation
def test_v_fn_flux_shared_1_scalar_from_angular() -> None:
    """V_fn-flux-shared.1 — universal angular-flux closure."""
    result = derive_scalar_flux_angular_integral()
    assert result["pass"], f"V_fn-flux-shared.1 failed: {result['name']}"
