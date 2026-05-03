r"""Foundation tests for the F_N method critical-slab solver (Sood
``Ua-1-0-SL`` / Case 2).

Test breakdown:

* **Branch-1 SymPy gates** — one ``@pytest.mark.foundation`` test per
  ``derive_*()`` in
  :mod:`orpheus.derivations.continuous.fn_method.origins.fn_slab_derivations`.
  These are pure-algebra: moment-integral recursion identities + the
  critical-determinant matrix structure.

* **Branch-1 ↔ Branch-2 numerical agreement** — the SymPy-verified
  closed-form moment recursion in numpy (:mod:`...core.moments`)
  satisfies the recursion to machine precision.

* **Branch-2 reference-value gates** — the F_N solver
  :func:`solve_fn_slab_bare_critical` reproduces the published
  Sood/KLL :math:`a_c = 0.93772556` mfp at :math:`c = 1.30` (Sood
  Case 2 ``Ua-1-0-SL``) to ≤ 1e-5 absolute, plus the
  Grandjean-Siewert Table XI critical thicknesses for
  :math:`c \in \{1.10, 1.30, 1.50, 1.70, 1.90\}`.

* **Convergence-order gate** — F_N error decreases with N (the
  expected algebraic convergence of the F_N method per Siewert-
  Benoist Part I).

References
----------

* Sood, Forster & Parsons 1999, LA-13511 (Table 4).
* Kaper, Lindeman & Leaf 1974, *Nucl. Sci. Eng.* **54**, 94 (Table I).
* Siewert & Benoist 1979, *Nucl. Sci. Eng.* **69**, 156-160.
* Grandjean & Siewert 1979, *Nucl. Sci. Eng.* **69**, 161-168 (Table XI).
"""
from __future__ import annotations

import pytest

# Suppress numpy divide-by-zero warnings from the F_N bracket scan
# (intermediate `a` values give near-singular complex matrices; we
# bisect through them legitimately).
pytestmark = [
    pytest.mark.filterwarnings(
        "ignore:divide by zero encountered in det:RuntimeWarning"
    ),
    pytest.mark.filterwarnings(
        "ignore:invalid value encountered in det:RuntimeWarning"
    ),
]

from orpheus.derivations.continuous.sood_registry import UA_1_0_SL_STUB
from orpheus.derivations.continuous.fn_method.core import (
    A_alpha,
    A_alpha_array,
    B_alpha,
    B_alpha_array,
    case_dispersion_root_supercritical,
    case_nu0,
    collocation_points,
)
from orpheus.derivations.continuous.fn_method.origins.fn_slab_derivations import (
    derive_A0_seed,
    derive_A_recursion,
    derive_B0_seed,
    derive_B_recursion,
    derive_critical_determinant_structure,
)
from orpheus.derivations.continuous.fn_method.slab import (
    solve_fn_slab_bare_critical,
)


# ═══════════════════════════════════════════════════════════════════
# Branch-1 SymPy gates — one test per derive_*() function
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_v_fn_slab_1_B_recursion():
    """V_fn-slab.1: B_α(ξ) recursion follows from algebraic long division."""
    r = derive_B_recursion()
    assert r["pass"], f"V_fn-slab.1 FAIL: {r}"


@pytest.mark.foundation
def test_v_fn_slab_2_A_recursion():
    """V_fn-slab.2: A_α(ξ) recursion follows from algebraic long division."""
    r = derive_A_recursion()
    assert r["pass"], f"V_fn-slab.2 FAIL: {r}"


@pytest.mark.foundation
def test_v_fn_slab_3_B0_seed():
    """V_fn-slab.3: B_0 long-division identity μ/(ξ-μ) = -1 + ξ/(ξ-μ)."""
    r = derive_B0_seed()
    assert r["pass"], f"V_fn-slab.3 FAIL: {r}"


@pytest.mark.foundation
def test_v_fn_slab_4_A0_seed():
    """V_fn-slab.4: A_0 seed integral closes."""
    r = derive_A0_seed()
    assert r["pass"], f"V_fn-slab.4 FAIL: {r}"


@pytest.mark.foundation
def test_v_fn_slab_5_critical_determinant():
    """V_fn-slab.5: critical-slab F_N system reduces to det M(a) = 0."""
    r = derive_critical_determinant_structure()
    assert r["pass"], f"V_fn-slab.5 FAIL: {r}"


# ═══════════════════════════════════════════════════════════════════
# Branch-1 SymPy ↔ Branch-2 numpy moment-recursion agreement
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_branch2_B_alpha_recursion_self_consistent():
    """Branch-2 B_alpha array satisfies the algebraic recursion to machine
    precision (no floating-point drift in the recursion loop)."""
    c = 1.30
    for xi in [0.05, 0.3, 0.5, 0.7, 0.95, 2.0, 5.0]:
        N = 10
        Bvec = B_alpha_array(N, xi, c)
        for n in range(1, N + 1):
            expected = xi * Bvec[n - 1] - 1.0 / (n + 1)
            assert abs(Bvec[n] - expected) < 1e-12, (
                f"B recursion drift at xi={xi}, n={n}: "
                f"got {Bvec[n]}, expected {expected}"
            )


@pytest.mark.foundation
def test_branch2_A_alpha_recursion_self_consistent():
    """Branch-2 A_alpha array satisfies the algebraic recursion to machine
    precision."""
    c = 1.30
    for xi in [0.05, 0.3, 0.5, 0.7, 0.95, 2.0, 5.0]:
        N = 10
        Avec = A_alpha_array(N, xi, c)
        for n in range(1, N + 1):
            expected = -xi * Avec[n - 1] + 1.0 / (n + 1)
            assert abs(Avec[n] - expected) < 1e-12, (
                f"A recursion drift at xi={xi}, n={n}: "
                f"got {Avec[n]}, expected {expected}"
            )


@pytest.mark.foundation
def test_branch2_B_alpha_array_matches_scalar():
    """Vectorized B_alpha_array equals the per-alpha B_alpha calls."""
    c = 1.20
    xi = 0.4
    N = 6
    arr = B_alpha_array(N, xi, c)
    for alpha in range(N + 1):
        scalar = B_alpha(alpha, xi, c)
        assert abs(arr[alpha] - scalar) < 1e-14


@pytest.mark.foundation
def test_branch2_A_alpha_array_matches_scalar():
    """Vectorized A_alpha_array equals the per-alpha A_alpha calls."""
    c = 1.20
    xi = 0.4
    N = 6
    arr = A_alpha_array(N, xi, c)
    for alpha in range(N + 1):
        scalar = A_alpha(alpha, xi, c)
        assert abs(arr[alpha] - scalar) < 1e-14


# ═══════════════════════════════════════════════════════════════════
# Dispersion-root sanity checks
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_case_dispersion_root_supercritical_1g():
    """Case dispersion u_0 satisfies 1 - c·u·atan(1/u) = 0 to machine precision."""
    import math
    for c in [1.10, 1.30, 1.50, 1.70, 1.90]:
        u0 = case_dispersion_root_supercritical(c)
        Lam = 1.0 - c * u0 * math.atan(1.0 / u0)
        assert abs(Lam) < 1e-14, f"c={c}: Λ(u_0={u0}) = {Lam}, expected 0."


@pytest.mark.foundation
def test_case_nu0_dispatch_supercritical():
    """case_nu0 dispatcher returns is_imag=True for c > 1."""
    u0, is_imag = case_nu0(1.30)
    assert is_imag is True
    # Sanity-check magnitude (Mitsis tabulated at c=1.30 gives u_0 ~ 0.946).
    assert 0.94 < u0 < 0.95


# ═══════════════════════════════════════════════════════════════════
# Collocation prescription
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_collocation_grid_grandjean_siewert():
    """Collocation grid {ν_0, 0, 1, equally-spaced (0,1)} matches GS prescription."""
    nu0 = 0.946
    # N = 5 (F_5): expected (ν_0, 0, 1, 0.25, 0.50, 0.75) per GS Table III implication.
    pts = collocation_points(nu0, 5)
    assert len(pts) == 6
    assert pts[0] == pytest.approx(nu0)
    assert pts[1] == 0.0
    assert pts[2] == 1.0
    # Interior 3 points: should be 0.25, 0.50, 0.75 (1/4 spacing on [0,1]).
    interior = sorted(pts[3:])
    assert interior == pytest.approx([0.25, 0.50, 0.75])


# ═══════════════════════════════════════════════════════════════════
# Branch-2 reference-value gates — Sood Ua-1-0-SL + Grandjean-Siewert XI
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_fn_slab_sood_ua_1_0_sl_critical_radius():
    """Sood ``Ua-1-0-SL`` (c=1.30): F_N reproduces a_c = 0.93772556 mfp to ≤ 1e-5."""
    res = solve_fn_slab_bare_critical(c=1.30, n_modes=10)
    truth = UA_1_0_SL_STUB.critical_dimension_mfp
    assert abs(res.a_critical_mfp - truth) < 1e-5, (
        f"Sood Ua-1-0-SL F_N: a={res.a_critical_mfp}, truth={truth}, "
        f"err={abs(res.a_critical_mfp - truth):.2e}"
    )
    # Also verify the determinant residual is small (should be ~ machine eps).
    assert abs(res.determinant_residual) < 1e-25


@pytest.mark.foundation
@pytest.mark.parametrize(
    ("c", "two_tau_truth"),
    [
        (1.10, 4.2266),
        (1.30, 1.8755),
        (1.50, 1.2101),
        (1.70, 0.8851),
        (1.90, 0.6919),
    ],
)
def test_fn_slab_grandjean_siewert_table_xi(c: float, two_tau_truth: float):
    """Grandjean-Siewert Table XI: F_5 reproduces 2τ to ≤ 5e-4 relative."""
    res = solve_fn_slab_bare_critical(c=c, n_modes=5)
    two_tau = 2.0 * res.a_critical_mfp
    rel_err = abs(two_tau - two_tau_truth) / two_tau_truth
    assert rel_err < 5e-4, (
        f"GS Table XI c={c}: 2τ={two_tau:.5f}, truth={two_tau_truth}, "
        f"rel_err={rel_err:.2e}"
    )


@pytest.mark.foundation
def test_fn_slab_convergence_with_N():
    """F_N error decreases as N grows (algebraic convergence)."""
    truth = 0.93772556
    errs = []
    for N in [3, 5, 8, 12]:
        res = solve_fn_slab_bare_critical(c=1.30, n_modes=N)
        errs.append(abs(res.a_critical_mfp - truth))
    # F_N error should generally decrease with N. Allow for the
    # well-known F_N collocation-quality oscillation (errors don't
    # decrease monotonically across all N because different ξ_β
    # placements affect the convergence pattern). We require:
    # (a) all errors ≤ 1e-3 absolute (loose envelope)
    # (b) the error at N=12 is smaller than at N=3.
    for N, err in zip([3, 5, 8, 12], errs):
        assert err < 1e-3, f"N={N}: err={err:.2e} exceeds 1e-3 envelope"
    assert errs[-1] < errs[0], (
        f"F_12 ({errs[-1]:.2e}) should beat F_3 ({errs[0]:.2e})"
    )


@pytest.mark.foundation
def test_fn_slab_coefficients_normalized():
    """F_N coefficients are normalised so a_0 = 1."""
    res = solve_fn_slab_bare_critical(c=1.30, n_modes=5)
    assert abs(res.coefficients_a[0] - 1.0) < 1e-12


@pytest.mark.foundation
def test_fn_slab_dispersion_consistency():
    """Returned u_0 satisfies the dispersion relation."""
    import math
    res = solve_fn_slab_bare_critical(c=1.30, n_modes=5)
    Lam = 1.0 - 1.30 * res.nu0 * math.atan(1.0 / res.nu0)
    assert abs(Lam) < 1e-14
