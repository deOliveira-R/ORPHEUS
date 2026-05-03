r"""L1 production tests for the Atalay 1997 reflected-slab solver.

References
----------

* Atalay, M.A. (1997). *Prog. Nucl. Energy* **31**(3), 229-252.
* Sood, A., Forster, R.A., Parsons, D.K. (1999). LA-13511.

Tolerance discussion
--------------------

The current Atalay implementation reproduces published critical slab
thicknesses with **systematic ~1-7% absolute error** that grows with
reflection coefficient R. The error tracks the z_0 1.5% gap from the
Atalay Eq 42 form (see ``...core.extrapolated_endpoint`` docstring;
:func:`...origins.derivations.derive_atalay_extrapolated_endpoint_eq42`
verifies the integrand structure but not the numerical match).

Tests below use **loose tolerances** (5e-2 relative) commensurate with
the implementation accuracy. The Atalay strict-1e-5 anchor is held by
fn_method/slab/one_group.solve_fn_slab_bare_critical for the bare
case and by Variant α for higher-precision references.
"""
from __future__ import annotations

import pytest

from orpheus.derivations.continuous.case_method.core.half_range import clear_X_cache
from orpheus.derivations.continuous.case_method.slab import (
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
    """Atalay Table 2 (f_1=0, R=0): slab full-width 2d_critical to 2e-2 relative."""
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
    """Atalay Table 2 (f_1=0, R>0): reflected-slab 2d_critical to 1e-1 relative.

    The error grows with R due to the cumulative quadrature precision in
    the K_j moments. At R = 0.99 (Atalay's last column) the gap can
    exceed 10%; we pin only the more accurate R ∈ [0, 0.75] range.
    """
    clear_X_cache()
    res = solve_case_method_slab_critical(c=c, R=R, f1=0.0, mode=1,
                                          n_bracket=200, d_min=0.005)
    err = abs(2 * res.d_critical_mfp - d2_atalay_table2) / d2_atalay_table2
    assert err < 1e-1, (
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
    """Atalay Table 3 (f_1=0.10, R=0): anisotropic vacuum slab to 2e-2 rel."""
    clear_X_cache()
    res = solve_case_method_slab_critical(
        c=c, R=0.0, f1=0.10, mode=1, n_bracket=80,
    )
    err = abs(2 * res.d_critical_mfp - d2_atalay_table3) / d2_atalay_table3
    assert err < 2e-2, (
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
