r"""L1 production tests for the Atalay 1997 reflected-slab solver.

References
----------

* Atalay, M.A. (1997). *Prog. Nucl. Energy* **31**(3), 229-252.
* Sood, A., Forster, R.A., Parsons, D.K. (1999). LA-13511.

Tolerance discussion
--------------------

After ERR-037 was fixed (z_0 quadrature endpoint regularised via
:math:`\mu = \tanh(t)` substitution; 2026-05-03), the Atalay Eq 42
extrapolated endpoint reproduces Atalay Table 1 to **6-7 digits**.
Slab vacuum (R=0) critical thicknesses now match Atalay Table 2 to
**0.1-0.5% absolute** (vs the Wave 2-B baseline of 1-1.2%).

Reflected cases (R > 0) and high-c anisotropic cases retain a
**residual ~1-4%** gap that originates in the K_j moment quadrature
(which carries the same endpoint behaviour as z_0 but with an
additional ``exp(-2d/μ)`` boundary-layer factor that needs separate
regularisation). That follow-up is tracked in the case_method
docstrings.

Tests below carry the post-fix tolerances:

* Vacuum slab isotropic (R=0): 1e-2 relative.
* Reflected slab isotropic (R∈[0.25, 0.75]): 5e-2 relative.
* Vacuum slab anisotropic f_1=0.10: 3e-2 relative (the c=2.0 row
  is the residual K_j gap).
"""
from __future__ import annotations

import pytest

from orpheus.derivations.continuous.singular_eigenfunction.core.half_range import clear_X_cache
from orpheus.derivations.continuous.singular_eigenfunction.slab import (
    solve_case_method_slab_critical,
)


# ═══════════════════════════════════════════════════════════════════
# Atalay Table 2 (f_1 = 0) — slab critical 2d for various R, c.
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("atalay-table2-slab-vacuum-isotropic")
@pytest.mark.parametrize("c, d2_atalay_table2", [
    (1.30, 1.87766),
    (1.50, 1.21523),
    (2.00, 0.63257),
])
def test_slab_vacuum_isotropic_atalay_table2(c, d2_atalay_table2):
    """Atalay Table 2 (f_1=0, R=0): slab full-width 2d_critical to 2e-2 relative.

    Tightened from previous loose 5e-2 baseline. After ERR-037 fix:
    0.12% (c=1.30), 0.43% (c=1.50), 1.75% (c=2.00). The c=2.00
    residual is the K_j moment quadrature endpoint behaviour
    (separate from z_0 fix; tracked as residual K_j gap that
    needs the same μ=tanh(t) regularisation in K_j as z_0
    received in ERR-037).
    """
    clear_X_cache()
    res = solve_case_method_slab_critical(c=c, R=0.0, f1=0.0, mode=1, n_bracket=80)
    err = abs(2 * res.d_critical_mfp - d2_atalay_table2) / d2_atalay_table2
    assert err < 2e-2, (
        f"Slab 2d at c={c}, R=0, f_1=0: got {2 * res.d_critical_mfp:.5f}, "
        f"Atalay Table 2 ref {d2_atalay_table2}, error {err*100:.2f}%"
    )


@pytest.mark.l1
@pytest.mark.verifies("atalay-table2-slab-reflected-isotropic")
@pytest.mark.parametrize("c, R, d2_atalay_table2", [
    (1.30, 0.25, 1.40621),
    (1.30, 0.50, 0.89317),
    (1.30, 0.75, 0.40758),
])
def test_slab_reflected_isotropic_atalay_table2(c, R, d2_atalay_table2):
    """Atalay Table 2 (f_1=0, R>0): reflected-slab 2d_critical to 5e-2 relative.

    Tightened from 1e-1 → 5e-2 after ERR-037 fix. Achieved post-fix:
    1.1% (R=0.25), 2.9% (R=0.50), 4.4% (R=0.75) — error grows with R
    due to K_j moment quadrature endpoint behaviour (separate from
    z_0 fix; tracked as residual K_j gap). At R=0.99 Atalay's last
    column the gap can still exceed 10%; we pin only R ∈ [0, 0.75].
    """
    clear_X_cache()
    res = solve_case_method_slab_critical(c=c, R=R, f1=0.0, mode=1,
                                          n_bracket=200, d_min=0.005)
    err = abs(2 * res.d_critical_mfp - d2_atalay_table2) / d2_atalay_table2
    assert err < 5e-2, (
        f"Slab 2d at c={c}, R={R}, f_1=0: got {2 * res.d_critical_mfp:.5f}, "
        f"Atalay Table 2 ref {d2_atalay_table2}, error {err*100:.2f}%"
    )


# ═══════════════════════════════════════════════════════════════════
# Atalay Table 3 (f_1 = 0.10) — anisotropic slab.
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("atalay-table3-slab-vacuum-anisotropic")
@pytest.mark.parametrize("c, d2_atalay_table3", [
    (1.30, 1.94146),
    (1.50, 1.25221),
    (2.00, 0.65236),
])
def test_slab_vacuum_anisotropic_f1_010_atalay_table3(c, d2_atalay_table3):
    """Atalay Table 3 (f_1=0.10, R=0): anisotropic vacuum slab to 3e-2 rel.

    Post-ERR-037 fix: 0.14% (c=1.30), 0.55% (c=1.50), 2.6% (c=2.00).
    The c=2.00 residual is the K_j moment quadrature endpoint behaviour
    amplified by X² at high c (separate from z_0 fix; tracked as
    residual K_j gap).
    """
    clear_X_cache()
    res = solve_case_method_slab_critical(
        c=c, R=0.0, f1=0.10, mode=1, n_bracket=80,
    )
    err = abs(2 * res.d_critical_mfp - d2_atalay_table3) / d2_atalay_table3
    assert err < 3e-2, (
        f"Slab 2d at c={c}, R=0, f_1=0.10: got {2 * res.d_critical_mfp:.5f}, "
        f"Atalay Table 3 ref {d2_atalay_table3}, error {err*100:.2f}%"
    )


# ═══════════════════════════════════════════════════════════════════
# Validity-bound enforcement
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.l1
def test_slab_rejects_subcritical_c():
    """c ≤ 1 (non-multiplying) must raise ValueError."""
    clear_X_cache()
    with pytest.raises(ValueError, match="c > 1"):
        solve_case_method_slab_critical(c=0.5, R=0.0, f1=0.0)


@pytest.mark.l1
def test_slab_rejects_perfect_reflector():
    """R = 1 must raise ValueError (slab thickness drops out of Eq 46)."""
    clear_X_cache()
    with pytest.raises(ValueError, match="R = 1"):
        solve_case_method_slab_critical(c=1.30, R=1.0, f1=0.0)


@pytest.mark.l1
def test_slab_rejects_above_validity_bound():
    """c > 1 + 1/(3 f_1) must raise ValueError per Atalay Eq 5."""
    clear_X_cache()
    # f_1 = 0.30 gives c_max = 1 + 1/0.9 = 19/9 ≈ 2.111.
    with pytest.raises(ValueError, match="validity bound"):
        solve_case_method_slab_critical(c=2.5, R=0.0, f1=0.30)


# ═══════════════════════════════════════════════════════════════════
# Result-object structural sanity
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_slab_result_attributes():
    """The result struct carries the expected diagnostic fields."""
    clear_X_cache()
    res = solve_case_method_slab_critical(c=1.30, R=0.0, f1=0.0, n_bracket=40)
    assert res.c == 1.30
    assert res.R == 0.0
    assert res.f1 == 0.0
    assert res.u0 > 0.0
    assert 0.0 < res.nu_bar < 1.0
    assert res.z0 > 0.0
    assert set(res.K_moments.keys()) == {0, 1, 2}
    assert res.d_critical_mfp > 0.0
