r"""Regression tests for ERR-037: Atalay Eq 42 z_0 quadrature endpoint.

Promoted from ``derivations/diagnostics/diag_05_z0_regression_atalay_table1.py``
on 2026-05-03 (numerics-investigator).

Bug history
-----------

Wave 2-B (2026-05-02) shipped the singular-eigenfunction slab/sphere
package (the former ``case_method`` folder, since consolidated into
:mod:`orpheus.derivations.continuous.singular_eigenfunction`) reproducing
Atalay Table 1 z_0 to ~1.5-2% relative across all c values. The
Wave 2-B closeout misdiagnosed this as a "Case-Zweifel
completeness-sum normalisation discrepancy" with the published
Atalay form, deferring the fix as "multi-day work."

Diagnosis (2026-05-03): the actual cause was quadrature endpoint
convergence at μ=1, where the integrand factor ``1 + c·μ²/(1-μ²)``
carries an integrable but slowly-cancelling pole. The fix is a
single-line variable substitution μ = tanh(t) that maps (0, 1) →
(0, ∞) with Jacobian sech²(t) cancelling the pole exactly.
Implementation lives in
``orpheus/derivations/continuous/singular_eigenfunction/core/extrapolated_endpoint.py``.

Reference
---------

* Atalay, M.A. (1997). *Prog. Nucl. Energy* **31**(3), 229-252.
  Table 1: ν̄ and z_0 vs c at f_1 = 0.

Cross-check
-----------

z_0 must satisfy the Davison/Milne-form identity
``z_0 = (π/2)·u_0 - a_c`` where ``a_c`` is the bare-slab critical
half-thickness. We verify this against Atalay's own Table 2 row
at c=1.30, R=0.
"""
from __future__ import annotations

import pytest

from orpheus.derivations.continuous.singular_eigenfunction.core.dispersion import (
    case_atalay_u0,
)
from orpheus.derivations.continuous.singular_eigenfunction.core.extrapolated_endpoint import (
    atalay_z0,
)


# Atalay 1997 Prog. Nucl. Energy Vol. 31 No. 3, Table 1 (f_1 = 0):
# Mean cosine ν̄ and extrapolated endpoint z_0 for the Milne problem.
ATALAY_TABLE_1 = [
    # (c, ν̄, z_0) — published 6 digits.
    (1.10, 0.694172, 0.645971),
    (1.20, 0.679619, 0.592392),
    (1.30, 0.666526, 0.547144),
    (1.40, 0.654682, 0.508410),
    (1.50, 0.643911, 0.474869),
    (1.60, 0.634071, 0.445536),
    (1.70, 0.625042, 0.419659),
    (1.80, 0.616723, 0.396659),
    (1.90, 0.609033, 0.376078),
    (2.00, 0.601898, 0.357551),
]


@pytest.mark.l1
@pytest.mark.parametrize("c, nu_bar_ref, z0_ref", ATALAY_TABLE_1)
@pytest.mark.catches("ERR-037")
@pytest.mark.verifies("atalay-eq42-extrapolated-endpoint")
def test_atalay_z0_table1_isotropic(c, nu_bar_ref, z0_ref):
    """ERR-037 regression: z_0 reproduces Atalay Table 1 to ≤ 1e-5
    after the μ=tanh(t) endpoint regularisation.

    Pre-fix (Wave 2-B baseline): z_0 missed Atalay Table 1 by ~1.5-2%
    relative across all c (e.g. c=1.30: computed 0.535568 vs Atalay
    0.547144). The Wave 2-B closeout misdiagnosed this as a
    completeness-sum normalisation issue.

    Post-fix (2026-05-03): z_0 matches Table 1 to 6-7 digits at all c
    (3.78e-7 worst-case at c=2.00).
    """
    u0 = case_atalay_u0(c, f1=0.0).u0
    z0 = atalay_z0(c=c, f1=0.0, u0=u0, nu_bar=nu_bar_ref)
    abs_err = abs(z0 - z0_ref)
    rel_err = abs_err / abs(z0_ref)
    assert abs_err < 1e-5, (
        f"c={c}: atalay_z0 = {z0:.7f}, Atalay Table 1 = {z0_ref:.6f}, "
        f"abs_err = {abs_err:.4e}, rel_err = {rel_err:.4e}"
    )


@pytest.mark.l1
@pytest.mark.catches("ERR-037")
def test_atalay_z0_consistent_with_milne_form_at_c130():
    """Cross-check: at c=1.30, z_0 from Atalay Eq 42 must equal
    ``(π/2)·u_0 - a_c`` where ``a_c`` is the bare-slab critical
    half-thickness from Atalay Table 2 (R=0).

    This is the **Milne-form definition** of z_0 (Davison's classical
    formula): z_0 = (π/2)·|ν_0| - a_c. Atalay's Table 1 z_0 is
    constructed precisely from Table 2's a_c via this relation.
    """
    import math

    from orpheus.derivations.continuous.singular_eigenfunction.core.half_range import (
        clear_X_cache,
    )

    c = 1.30
    clear_X_cache()
    u0 = case_atalay_u0(c, f1=0.0).u0
    z0 = atalay_z0(c=c, f1=0.0, u0=u0, nu_bar=0.666526)

    # Atalay Table 2 R=0: 2d_c = 1.87766 → a_c = 0.93883
    atalay_table2_2d = 1.87766
    a_c = atalay_table2_2d / 2
    z0_milne = (math.pi / 2) * u0 - a_c

    abs_diff = abs(z0 - z0_milne)
    print(
        f"\n  c=1.30 cross-check:\n"
        f"  z_0 from Eq 42:    {z0:.7f}\n"
        f"  (π/2)·u_0 - a_c:   {z0_milne:.7f}\n"
        f"  diff:              {abs_diff:.4e}"
    )
    # Atalay's own self-consistency between Tables 1 and 2 is at the
    # last printed digit (~5e-6); allow 1e-4 to absorb any mode-1
    # truncation in the Table 2 critical-thickness solve.
    assert abs_diff < 1e-4, (
        f"z_0 from Eq 42 ({z0}) inconsistent with Milne form ({z0_milne})"
    )
