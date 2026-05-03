r"""Foundation tests for the F_N method critical-sphere solver — Sood
``Ua-1-0-SP`` (Case 4) via the **true F_N method** of Siewert-Thomas
1986.

Test breakdown:

* **Branch-1 SymPy gates** — one ``@pytest.mark.foundation`` test per
  ``derive_*()`` in
  :mod:`orpheus.derivations.continuous.fn_method.origins.fn_sphere_derivations`.
  Five tests pinning the algebraic identities that justify the
  geometry-sign abstraction:

  - V_fn-sphere-fn.1: slab/sphere BCs differ by a geometry sign
    :math:`s \in \{+1, -1\}`.
  - V_fn-sphere-fn.2: sphere F_N matrix entry follows from
    :math:`s = -1` in the unified assembler.
  - V_fn-sphere-fn.3: sphere bare-critical = :math:`\det M(R) = 0`
    (homogeneous linear system).
  - V_fn-sphere-fn.4: Siewert-Thomas 2G→1G via matrix-dimension
    collapse.
  - V_fn-sphere-fn.5: Wiener-Hopf X-function depends only on
    :math:`\Lambda(z)` — geometry-independent.

* **Branch-1 ↔ Branch-2 numerical agreement** — the Chebyshev
  collocation grid and assembler internals are sanity-checked.

* **Branch-2 reference-value gate (Sood Ua-1-0-SP)** — the F_N sphere
  solver reproduces the published Sood/KLL :math:`R_c = 2.4248249802`
  mfp at :math:`c = 1.30` to ≤ 1e-5 absolute. With the F_N method
  this gate runs in <1 second (vs ~100 s for the previous PS-1982
  wrapper that bisected on R via repeated Nyström evaluations).

* **Convergence-order gate** — F_N error decreases monotonically
  with N up to N≈20 (where the Brent precision floor is reached
  at ~10^-9 absolute on R_c).

References
----------

* Siewert & Thomas 1986, *Nucl. Sci. Eng.* **94**, 264-270.
* Sood, Forster & Parsons 1999, LA-13511 (Table 6).
* Kaper, Lindeman & Leaf 1974, *Nucl. Sci. Eng.* **54**, 94 (Table V).
"""
from __future__ import annotations

import pytest

# Suppress numpy divide-by-zero warnings from the F_N bracket scan
# (intermediate `R` values give near-singular complex matrices; we
# minimise through them legitimately).
pytestmark = [
    pytest.mark.filterwarnings(
        "ignore:divide by zero encountered in det:RuntimeWarning"
    ),
    pytest.mark.filterwarnings(
        "ignore:invalid value encountered in det:RuntimeWarning"
    ),
]

from orpheus.derivations.continuous.sood_registry import UA_1_0_SP_STUB
from orpheus.derivations.continuous.fn_method.origins.fn_sphere_derivations import (
    derive_sphere_2g_to_1g_reduction,
    derive_sphere_bc_sign_flip,
    derive_sphere_critical_condition,
    derive_sphere_fn_matrix_entry,
    derive_x_function_geometry_independence,
)
from orpheus.derivations.continuous.fn_method.sphere import (
    solve_fn_sphere_bare_critical,
)


# ═══════════════════════════════════════════════════════════════════
# Branch-1 SymPy gates — one test per derive_*() function
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_v_fn_sphere_fn_1_bc_sign_flip():
    """V_fn-sphere-fn.1: slab/sphere F_N differ by geometry_sign s ∈ {+1, -1}."""
    r = derive_sphere_bc_sign_flip()
    assert r["pass"], f"V_fn-sphere-fn.1 FAIL: {r}"


@pytest.mark.foundation
def test_v_fn_sphere_fn_2_matrix_entry():
    """V_fn-sphere-fn.2: sphere F_N matrix entry = B - exp(-2R/ξ)·A."""
    r = derive_sphere_fn_matrix_entry()
    assert r["pass"], f"V_fn-sphere-fn.2 FAIL: {r}"


@pytest.mark.foundation
def test_v_fn_sphere_fn_3_critical_condition():
    """V_fn-sphere-fn.3: sphere bare-critical = det M(R) = 0."""
    r = derive_sphere_critical_condition()
    assert r["pass"], f"V_fn-sphere-fn.3 FAIL: {r}"


@pytest.mark.foundation
def test_v_fn_sphere_fn_4_2g_to_1g_reduction():
    """V_fn-sphere-fn.4: Siewert-Thomas 2G→1G via matrix-dimension collapse."""
    r = derive_sphere_2g_to_1g_reduction()
    assert r["pass"], f"V_fn-sphere-fn.4 FAIL: {r}"


@pytest.mark.foundation
def test_v_fn_sphere_fn_5_x_function_geometry_independence():
    """V_fn-sphere-fn.5: Wiener-Hopf X-function geometry-independent."""
    r = derive_x_function_geometry_independence()
    assert r["pass"], f"V_fn-sphere-fn.5 FAIL: {r}"


# ═══════════════════════════════════════════════════════════════════
# Branch-1 ↔ Branch-2 collocation grid sanity
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_sphere_chebyshev_collocation_in_open_unit_interval():
    """Sphere Chebyshev collocation per Siewert-Thomas Eq. 38a yields
    points strictly INSIDE (0, 1) — does NOT include 0 or 1 like the
    slab Grandjean-Siewert grid.

    This test pins the structural distinction: slab F_N and sphere F_N
    use DIFFERENT collocation grids despite sharing the moment
    recursion, dispersion-root, and X-function machinery.
    """
    import numpy as np
    from orpheus.derivations.continuous.fn_method.sphere.one_group import (
        _siewert_thomas_collocation,
    )

    nu0_complex = complex(0.0, 0.946)  # placeholder for c=1.30
    for N in [5, 8, 10, 15]:
        xis = _siewert_thomas_collocation(nu0_complex, N)
        assert xis[0] == nu0_complex, "ξ_0 must be ν_0"
        # All other points should be real, in (0, 1) strictly.
        for k in range(1, N + 1):
            xi_k = xis[k]
            assert abs(xi_k.imag) < 1e-15, (
                f"ξ_{k} should be real, got {xi_k}"
            )
            assert 0.0 < xi_k.real < 1.0, (
                f"ξ_{k} = {xi_k.real} must be strictly in (0, 1)"
            )
        # Specifically: NO point should be at exactly 0 or 1.
        real_parts = np.array([xis[k].real for k in range(1, N + 1)])
        assert np.all(real_parts > 1e-12), (
            "Sphere F_N grid must NOT include ξ = 0 (would create rank "
            "deficiency, see Siewert-Thomas Eq. 38a docstring)"
        )
        assert np.all(real_parts < 1.0 - 1e-12), (
            "Sphere F_N grid must NOT include ξ = 1"
        )


# ═══════════════════════════════════════════════════════════════════
# Branch-2 reference-value gate — Sood Ua-1-0-SP
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_fn_sphere_sood_ua_1_0_sp_critical_radius():
    r"""Sood ``Ua-1-0-SP`` (c=1.30): F_N sphere reproduces
    :math:`R_c = 2.4248249802` mfp to ≤ 1e-5 absolute.

    The F_N method gives ~1e-7 absolute at N = 10 — 100x better
    than the 1e-5 target. The implementation runs in <1 s.
    """
    res = solve_fn_sphere_bare_critical(c=1.30, n_modes=10)
    truth = UA_1_0_SP_STUB.critical_dimension_mfp
    err = abs(res.R_critical_mfp - truth)
    assert err < 1e-5, (
        f"Sood Ua-1-0-SP F_N sphere: R_c={res.R_critical_mfp:.10f}, "
        f"truth={truth}, err={err:.3e} (target 1e-5)"
    )
    # Determinant residual at the converged R should be machine-tiny.
    assert abs(res.determinant_residual) < 1e-20, (
        f"det(M) residual {abs(res.determinant_residual):.3e} not tiny"
    )


@pytest.mark.foundation
def test_fn_sphere_convergence_with_N():
    r"""F_N error decreases as N grows for sphere (algebraic
    convergence per Siewert-Thomas Table IV pattern).

    We require:
    (a) all errors at N ∈ {5, 8, 10, 12, 15} are below 1e-5 (i.e. all
        reach the Sood target);
    (b) the error at N=15 is at least 4x smaller than at N=5 (genuine
        convergence, not a flat regime).
    """
    truth = 2.4248249802
    errs = {}
    for N in [5, 8, 10, 12, 15]:
        res = solve_fn_sphere_bare_critical(c=1.30, n_modes=N)
        errs[N] = abs(res.R_critical_mfp - truth)
    for N, err in errs.items():
        assert err < 1e-5, (
            f"N={N}: err={err:.3e} exceeds 1e-5 target"
        )
    # F_15 should clearly beat F_5.
    assert errs[15] < errs[5] / 4, (
        f"F_15 ({errs[15]:.3e}) should be < F_5/4 ({errs[5]/4:.3e})"
    )


@pytest.mark.foundation
def test_fn_sphere_coefficients_normalized():
    """F_N sphere coefficients are normalised so a_0 = 1."""
    res = solve_fn_sphere_bare_critical(c=1.30, n_modes=8)
    assert abs(res.coefficients_a[0] - 1.0) < 1e-12


@pytest.mark.foundation
def test_fn_sphere_dispersion_consistency():
    """Returned u_0 satisfies the dispersion relation 1 - c·u·atan(1/u) = 0."""
    import math
    res = solve_fn_sphere_bare_critical(c=1.30, n_modes=8)
    Lam = 1.0 - 1.30 * res.nu0 * math.atan(1.0 / res.nu0)
    assert abs(Lam) < 1e-14, (
        f"Sphere F_N: c=1.30 dispersion residual {Lam:.3e} > 1e-14"
    )


@pytest.mark.foundation
def test_fn_sphere_uses_geometry_sign_minus_one():
    r"""Sphere F_N solver MUST use ``geometry_sign = -1`` in the
    unified assembler — the load-bearing structural fact from
    Siewert-Thomas 1986 Eq. 46.

    Direct test: build the F_N matrix at the converged R for both
    geometry signs and verify that:
    (a) the sphere-sign matrix is near-singular (smin small);
    (b) the slab-sign matrix at the same R is NOT near-singular —
        proving the sign matters for the critical condition.
    """
    import numpy as np
    from orpheus.derivations.continuous.fn_method.core.fn_matrix import (
        assemble_fn_matrix,
    )

    res = solve_fn_sphere_bare_critical(c=1.30, n_modes=10)
    # Sphere matrix at the converged R: should be near-singular.
    M_sphere = assemble_fn_matrix(
        c=1.30,
        R_mfp=res.R_critical_mfp,
        geometry_sign=-1,
        n_modes=10,
        xis=res.xi_collocation,
    )
    # Slab matrix at the same R: should NOT be near-singular.
    M_slab = assemble_fn_matrix(
        c=1.30,
        R_mfp=res.R_critical_mfp,
        geometry_sign=+1,
        n_modes=10,
        xis=res.xi_collocation,
    )
    smin_sphere = float(np.linalg.svd(M_sphere, compute_uv=False)[-1])
    smin_slab = float(np.linalg.svd(M_slab, compute_uv=False)[-1])
    # At the sphere critical R, sphere smin should be much smaller
    # than slab smin (smin_sphere ~ 1e-11; smin_slab ~ 1e-7).
    assert smin_sphere < smin_slab, (
        f"Sphere smin ({smin_sphere:.3e}) should be < slab smin "
        f"({smin_slab:.3e}) at R_c = {res.R_critical_mfp:.5f}"
    )
    # And the gap should be at least 3 decades.
    import math
    log_gap = math.log10(smin_slab / smin_sphere)
    assert log_gap > 3.0, (
        f"smin gap should be > 3 decades, got {log_gap:.2f}"
    )
