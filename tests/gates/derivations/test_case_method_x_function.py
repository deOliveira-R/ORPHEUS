"""Regression pins for Atalay 1997 Eq 40 X-function quadrature.

Documents and pins the **logarithmically-divergent integrand** that was
characterised on 2026-05-03 as the source of the residual 1-4 % gap on
Atalay Table 2 reflected-slab critical thicknesses (post-ERR-037).

ERR-037 fixed the z_0 evaluator (Atalay Eq 42) by replacing plain
``mp.quad`` over (0,1) with the ``μ = tanh(t)`` substitution that
exactly cancels the bracket's ``1/(1-μ²)`` algebraic pole. The natural
follow-up hypothesis was that the same substitution should fix the K_j
moment quadrature, but a diagnostic cascade (diag_kj_phase1_*) found:

* The K_j outer ``scipy.quad`` converges fine (rel drift ~1e-11
  between baseline and tight tolerance settings).
* The bottleneck is X(-ν), which propagates ~1e-3 rel error into K_j
  via the X² factor.
* Atalay Eq 40's X-function integrand is **logarithmically divergent**
  at ν → 1 (primitive ∝ -1/log(1-ν), bounded but truncation-dependent).
* Both ``(0,1)`` direct quadrature and the ``μ=tanh(t)`` substituted
  ``(0,∞)`` form give finite but **algorithm-dependent** values.
* Implementing the substitution in both X-function backends shifted
  the slab gap from 4.40 % → 4.44 % (marginal degradation). Reverted.

These tests **pin the post-ERR-037 baseline** so that:

1. If a future investigation closes the X-function gap, the pin fires
   and the tolerance must be tightened.
2. If a regression introduces a different X-function bug, the pin
   fires with a different sign.

References:

* Atalay 1997 Prog. Nucl. Energy 31(3) Eq 40 — X-function calculation
  formula (the divergent representation).
* Atalay 1997 Eq 26 — implicit definition of X via ``∫γ d(ν²)/(ν-μ) dν``.
* Case-Plazcek-Hofmann 1961 — closed-form X-function for isotropic
  scattering (alternative formulation, not currently in ORPHEUS).
* ``error_catalog.md`` ERR-037 entry "extended-fix hypothesis FALSIFIED"
  for the full reasoning trail.
"""
from __future__ import annotations

import mpmath as mp
import pytest


@pytest.mark.l0
@pytest.mark.verifies("singular-eigenfunction-eq40")
@pytest.mark.parametrize("c, nu_eval", [(1.30, 0.55), (1.30, 0.85)])
def test_x_integrand_logarithmic_divergence_fingerprint(c, nu_eval):
    """The Atalay Eq 40 integrand drifts >1e-4 across mpmath dps 15→60.

    For a CONVERGENT integrand mpmath would give bit-stable values
    across dps. The drift > 1e-4 confirms the integrand is divergent
    near ν = 1 — i.e. the X-function value is regularisation-dependent.

    This is a **necessary precondition** for the residual K_j gap pin
    below to be meaningful: the gap exists because production code
    uses one regularisation (scipy float-precision endpoint sampling)
    while Atalay's published Table 2 values implicitly use a different
    one.

    If a future fix makes mpmath bit-stable across dps (i.e. the
    integrand becomes convergent), this test FIRES — at which point
    the regularisation question has been resolved and the tolerance
    on ``test_residual_kj_gap_pinned`` should be tightened.
    """
    f1 = 0.0
    mu = -nu_eval

    def integrand(nu, c_val):
        bracket = 1 + c_val * nu * nu / (1 - nu * nu)
        lam = 1 - c_val * nu * mp.atanh(nu)
        pi_term = mp.pi * c_val * nu / 2
        g1 = 1 / (lam * lam + pi_term * pi_term)
        return g1 * bracket * mp.log(nu - mu)

    values = {}
    for dps in (15, 60):
        with mp.workdps(dps):
            c_mp = mp.mpf(c)
            integral = mp.quad(
                lambda nu: integrand(nu, c_mp), [0, 1], maxdegree=15,
            )
            values[dps] = float(mp.exp(-c_mp / 2 * integral))

    drift = abs(values[60] - values[15])
    assert drift > 1e-4, (
        f"c={c}, ν={nu_eval}: expected divergent-integrand drift > 1e-4 "
        f"between mpmath dps=15 and dps=60, got {drift:.3e}. "
        f"If the integrand has become convergent (drift → 0 as dps → ∞), "
        f"the Atalay Eq 40 regularisation question has been resolved — "
        f"tighten test_residual_kj_gap_pinned tolerances."
    )


@pytest.mark.l1
@pytest.mark.verifies("singular-eigenfunction-eq40", "singular-eigenfunction-eq46")
@pytest.mark.parametrize("c, R, d2_atalay_table2, expected_pct", [
    (1.30, 0.25, 1.40621, 1.13),  # 1.13% post-ERR-037 baseline
    (1.30, 0.50, 0.89317, 2.87),  # 2.87% post-ERR-037 baseline
    (1.30, 0.75, 0.40758, 4.40),  # 4.40% post-ERR-037 baseline
])
def test_residual_kj_gap_pinned(c, R, d2_atalay_table2, expected_pct):
    """Pin post-ERR-037 residual gap on Atalay Table 2 reflected slabs.

    These three cases retain a 1-4 % gap after ERR-037 fixed z_0. The
    gap originates in the X-function regularisation (Atalay Eq 40 is
    divergent — see the fingerprint test above). The μ=tanh(t)
    substitution does NOT close this gap — when applied to the X-function
    quadrature, the gap shifts marginally (4.40 % → 4.44 % on R=0.75).

    Pin tolerance: gap must be within ±50 % of the expected baseline.
    A wider miss means either (a) a regression in the production X-function
    quadrature, or (b) a fix has been applied — in either case investigate.
    """
    from orpheus.derivations.continuous.singular_eigenfunction.core.half_range import clear_X_cache
    from orpheus.derivations.continuous.singular_eigenfunction.slab.one_group import (
        solve_case_method_slab_critical,
    )

    clear_X_cache()
    res = solve_case_method_slab_critical(
        c=c, R=R, f1=0.0, mode=1, n_bracket=200, d_min=0.005,
    )
    rel_err_pct = 100.0 * abs(2 * res.d_critical_mfp - d2_atalay_table2) / d2_atalay_table2

    assert 0.5 * expected_pct < rel_err_pct < 1.5 * expected_pct, (
        f"c={c}, R={R}: rel_err {rel_err_pct:.3f}% outside expected band "
        f"[{0.5*expected_pct:.2f}%, {1.5*expected_pct:.2f}%]. "
        f"Production: 2d = {2*res.d_critical_mfp:.6f}, "
        f"Atalay Table 2: 2d = {d2_atalay_table2}. "
        f"If the gap CLOSED, the X-function regularisation has been resolved; "
        f"tighten this tolerance and update the docstring. "
        f"If the gap GREW, look for a regression in the X-function or "
        f"K_j moment quadrature."
    )
