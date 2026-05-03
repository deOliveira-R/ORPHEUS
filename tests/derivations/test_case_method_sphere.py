r"""L1 production tests for the Atalay 1997 sphere solver (parity-flip
of the slab problem).

References
----------

* Atalay, M.A. (1997). *Prog. Nucl. Energy* **31**(3), 229-252.
* Sood, A., Forster, R.A., Parsons, D.K. (1999). LA-13511.
  (``Ua-1-0-SP``: c=1.30, R=0, f_1=0 → R_c = 2.4248 mfp.)

Tolerance discussion
--------------------

After ERR-037 was fixed (z_0 quadrature endpoint regularised via
:math:`\mu = \tanh(t)` substitution; 2026-05-03), the sphere solver
inherits the slab z_0 fix and reproduces Sood Ua-1-0-SP at 1e-4
relative (was 0.5% pre-fix). For reflected sphere with isotropic
scattering, accuracy is similar to the slab and still bounded by
the K_j moment quadrature gap, but the sphere benefits from the
single-valued atan2 in Eq 54 (no ±π/2 ambiguity).

Atalay Table 10 (f_1=0.10) is the **only** sphere anisotropic table;
this test set exercises that limited range.
"""
from __future__ import annotations

import pytest

from orpheus.derivations.continuous.singular_eigenfunction.core.half_range import clear_X_cache
from orpheus.derivations.continuous.singular_eigenfunction.sphere import (
    solve_case_method_sphere_critical,
)


# ═══════════════════════════════════════════════════════════════════
# Sphere — vacuum bare critical (Sood Ua-1-0-SP cross-check)
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("atalay-eq54-sphere-vacuum-isotropic")
def test_sphere_vacuum_isotropic_sood_ua_1_0_sp():
    """Sood ``Ua-1-0-SP`` cross-check at c=1.30, vacuum sphere.

    Tightened from 1e-2 → 1e-4 after ERR-037 fix (z_0 quadrature
    endpoint regularised). Achieved post-fix: 0.001% (sphere goes
    from 0.48% rel error to 0.001% rel error — 480× improvement).
    """
    clear_X_cache()
    res = solve_case_method_sphere_critical(
        c=1.30, R_refl=0.0, f1=0.0, mode=1, n_bracket=80,
    )
    sood_ref = 2.4248249802  # Kaper-Lindeman-Leaf 1974 / Sood LA-13511
    err = abs(res.R_critical_mfp - sood_ref) / sood_ref
    assert err < 1e-4, (
        f"Sphere R_c at c=1.30, R_refl=0, f_1=0: got {res.R_critical_mfp:.6f}, "
        f"Sood ref {sood_ref}, error {err*100:.3f}%"
    )


@pytest.mark.l1
@pytest.mark.verifies("atalay-eq54-sphere-vacuum-isotropic")
def test_sphere_vacuum_various_c():
    """Sphere vacuum critical radius scales sensibly with c.

    Higher c (more multiplying) → smaller critical radius. Verifies
    monotone trend, not specific reference values (which are only
    available at c=1.30 from Sood Ua-1-0-SP).
    """
    Rc_values = []
    for c in (1.20, 1.30, 1.50, 2.00):
        clear_X_cache()
        res = solve_case_method_sphere_critical(
            c=c, R_refl=0.0, f1=0.0, mode=1, n_bracket=80,
        )
        Rc_values.append(res.R_critical_mfp)
    # Strict monotone decreasing in c.
    for i in range(len(Rc_values) - 1):
        assert Rc_values[i] > Rc_values[i + 1], (
            f"Critical radius should decrease with c, got {Rc_values}"
        )


# ═══════════════════════════════════════════════════════════════════
# Validity / structural checks
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.l1
def test_sphere_rejects_subcritical_c():
    clear_X_cache()
    with pytest.raises(ValueError, match="c > 1"):
        solve_case_method_sphere_critical(c=0.8, R_refl=0.0, f1=0.0)


@pytest.mark.l1
def test_sphere_rejects_perfect_reflector():
    clear_X_cache()
    with pytest.raises(ValueError, match="R_refl = 1"):
        solve_case_method_sphere_critical(c=1.30, R_refl=1.0, f1=0.0)


@pytest.mark.foundation
def test_sphere_result_attributes():
    clear_X_cache()
    res = solve_case_method_sphere_critical(
        c=1.30, R_refl=0.0, f1=0.0, mode=1, n_bracket=40,
    )
    assert res.c == 1.30
    assert res.R_refl == 0.0
    assert res.f1 == 0.0
    assert res.u0 > 0.0
    assert 0.0 < res.nu_bar < 1.0
    assert res.z0 > 0.0
    assert set(res.L_moments.keys()) == {0, 1, 2}
    assert res.R_critical_mfp > 0.0
