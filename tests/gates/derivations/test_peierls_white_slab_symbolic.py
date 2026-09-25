r"""Foundation rows for the SymPy origin of the slab white-boundary closed
form, :mod:`orpheus.derivations.continuous.peierls_nystrom.origins.white_slab`.

Two rows, one per ``derive_*`` function:

- V_ws1: :math:`\int_0^\tau E_2 = \tfrac12 - E_3(\tau)` and
  :math:`\int_0^\tau E_1 = 1 - E_2(\tau)`, by the fundamental theorem of
  calculus (derivative residual zero, value at zero zero). Its negative
  leg is the ERR-032 candidate :math:`1 - E_3`, which passes the
  derivative condition and fails the anchor (value :math:`\tfrac12`).
- V_ws2: the Peierls-plus-Mark-balance flux assembled with the true
  antiderivative is identically :math:`1/\Sigma_t`
  (:eq:`peierls-white-bc-slab`), and with the ERR-032 antiderivative it
  is not (residual 0.45097 at the arm's probe point).

These rows prove the ALGEBRA and evaluate no production code, so no
code defect can redden them: they carry no ``catches`` marker (X3). The
code-side ERR-032 catchers are ``TestSlabWhiteBCInfiniteMediumIdentity``
and ``TestSlabWhiteBCPartialCurrentBalance`` in
``tests/gates/derivations/test_peierls_reference.py``. Each row's tooth is
its own negative leg. V_ws2 carries ``verifies`` under the
algebra-of-record Branch-1 exception to the foundation rule
(``vv-principles``, the level taxonomy): it is a SymPy identity pinned
against a theory-page label.
"""
from __future__ import annotations

import mpmath
import pytest
import sympy as sp

from orpheus.derivations.continuous.peierls_nystrom.origins.white_slab import (
    derive_en_antiderivatives,
    derive_slab_white_bc_flux,
)

_V_WS1 = (
    "tests/gates/derivations/test_peierls_white_slab_symbolic.py"
    "::test_v_ws1_en_antiderivatives_anchored_at_zero"
)


@pytest.mark.foundation
def test_v_ws1_en_antiderivatives_anchored_at_zero() -> None:
    r"""The two antiderivatives hold, and the ERR-032 candidate differs
    from the true one by exactly the integration constant ½."""
    result = derive_en_antiderivatives()
    assert result["E2_derivative_residual"] == 0, result["E2_derivative_residual"]
    assert result["E2_value_at_zero"] == 0, result["E2_value_at_zero"]
    assert result["E1_derivative_residual"] == 0, result["E1_derivative_residual"]
    assert result["E1_value_at_zero"] == 0, result["E1_value_at_zero"]
    # Negative leg: same derivative, wrong anchor. A derivative-only check
    # would pass 1 - E_3; the anchor is the discriminating condition.
    assert result["err032_derivative_residual"] == 0, (
        "the ERR-032 candidate 1 - E_3 should share the derivative E_2; "
        f"got residual {result['err032_derivative_residual']}"
    )
    assert result["err032_value_at_zero"] == sp.Rational(1, 2), (
        "the ERR-032 candidate 1 - E_3 should be anchored at 1/2, not 0; "
        f"got {result['err032_value_at_zero']}"
    )
    assert result["pass"] is True


@pytest.mark.foundation
@pytest.mark.verifies("peierls-white-bc-slab")
@pytest.mark.rests_on(_V_WS1)
def test_v_ws2_white_bc_slab_flux_is_one_over_sigma_t() -> None:
    r"""The assembled flux minus :math:`1/\Sigma_t` simplifies to 0; the
    same assembly with :math:`1 - E_3` leaves 0.45097 at
    :math:`(x, L, \Sigma_t) = (0.3, 1, 1)` (the catalogue's wrong closed
    form, 1.45097, minus the true 1)."""
    result = derive_slab_white_bc_flux()
    assert result["residual"] == 0, (
        f"phi - 1/Sigma_t did not simplify to 0: {result['residual']}"
    )
    # The catalogue's wrong closed form, spelled directly from the entry
    # (phi_wrong = (1/(2 Sigma_t)) [2 + (2 beta - 1)(E_2(Sigma_t x) +
    # E_2(Sigma_t (L - x)))], beta = (1 - E_3)/(1 - 2 E_3)), evaluated by
    # mpmath: the SymPy assembly with 1 - E_3 must reproduce it.
    with mpmath.workdps(30):
        x, L, sig_t = mpmath.mpf("0.3"), mpmath.mpf(1), mpmath.mpf(1)
        e3 = mpmath.expint(3, sig_t * L)
        beta = (1 - e3) / (1 - 2 * e3)
        e2_sum = mpmath.expint(2, sig_t * x) + mpmath.expint(2, sig_t * (L - x))
        phi_wrong = (2 + (2 * beta - 1) * e2_sum) / (2 * sig_t)
        catalogue_residual = phi_wrong - 1 / sig_t
        sympy_residual = mpmath.mpf(str(result["err032_residual_at_probe"]))
        gap = abs(sympy_residual - catalogue_residual)
    assert gap < mpmath.mpf("1e-25"), (
        f"the SymPy ERR-032 assembly {sympy_residual} does not reproduce the "
        f"catalogue's wrong closed form {catalogue_residual} (gap {gap})"
    )
    assert catalogue_residual > mpmath.mpf("0.1"), catalogue_residual
    assert result["pass"] is True
