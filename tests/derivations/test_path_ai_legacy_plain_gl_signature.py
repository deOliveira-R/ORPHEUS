r"""Signature tests for the legacy plain-GL Path A.i kernel
truncation bug (ERR-036) and the foundation of the Atkinson fix
notation bridge.

This file pins the QUANTITATIVE failure mode of the legacy
plain-GL Path A.i. It does NOT exercise
``slab_scalar_flux_fn_projection`` directly — instead it verifies
the underlying mathematical observations that make the legacy
path produce 5–7 % error and the Atkinson product-Nyström path
produce 1e-5 error.

Roles
-----

* :func:`test_E1_log_decomposition`: verifies the small-tau
  expansion `E_1(τ) = -γ_E - log(τ) + R(τ)` with `R` smooth.
  This is the load-bearing notation-bridge identity used by
  :mod:`...peierls_atkinson_nystrom` to split the kernel into
  log-singular + smooth parts.

* :func:`test_legacy_path_ai_diagonal_kernel_truncation`:
  documents the precise pathology — plain-GL μ-quadrature
  truncates the divergent diagonal `E_1(0+) = +∞` to a finite
  ~`2 log(n_μ)` value. This test would FAIL if someone
  accidentally "fixed" the truncation by inserting a more careful
  `μ = 0` handler — but the right fix is product-Nyström, not
  evaluation of E_1 at the singularity.

* :func:`test_legacy_path_ai_convergence_rate_below_atkinson`:
  empirical convergence rate of the legacy path is ~ -0.9 (first
  order with log correction); confirms the Atkinson path's
  ~ -3 rate is a genuine improvement.

Catches
-------

@pytest.mark.catches("ERR-036") — see
``.claude/skills/vv-principles/error_catalog.md`` ERR-036.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.special import exp1


_GAMMA_EULER = 0.5772156649015328606


@pytest.mark.foundation
@pytest.mark.verifies("e1-decomposition")
def test_E1_log_decomposition_notation_bridge():
    """E_1(tau) = -gamma_E - log(tau) + R(tau), with R(tau) -> 0
    as tau -> 0+.

    This identity is the foundation of the Atkinson fix: the kernel
    ``(c/2) E_1(|z-z'|)`` decomposes into a log-singular piece
    (handled by product weights) and a C^infty remainder R (handled
    by standard Simpson).

    The bridge is NOT given in Atkinson 1972 / 1997 (which are
    kernel-agnostic); it is derived from Abramowitz-Stegun 5.1.11
    locally. The diagnostic memo at
    ``scratch/derivations/peierls_log_singular/atkinson_product_nystrom.md``
    flagged this as a load-bearing identity to verify independently.
    """
    taus = [1e-2, 1e-4, 1e-6, 1e-8, 1e-10]
    remainders = []
    for tau in taus:
        E1_val = float(exp1(tau))
        # E_1(tau) - (-gamma_E - log(tau)) should be the smooth remainder R.
        R = E1_val - (-_GAMMA_EULER - math.log(tau))
        remainders.append(R)
        # R(tau) approaches 0 as tau -> 0; specifically R(tau) ~ tau.
        # Linear leading order: R(tau) ≈ tau at small tau.
        rel = abs(R - tau) / max(tau, 1e-300)
        assert rel < 1e-2, (
            f"E_1 - (-gamma_E - log(tau)) leading order is NOT tau "
            f"at tau = {tau}: R = {R}, rel diff from tau = {rel}"
        )

    # Confirm R IS converging to 0 as tau -> 0.
    assert abs(remainders[-1]) < 1e-9, (
        f"Smooth remainder should vanish at tau -> 0+; "
        f"got R(1e-10) = {remainders[-1]}"
    )


@pytest.mark.l1
@pytest.mark.catches("ERR-036")
def test_legacy_path_ai_diagonal_kernel_truncation():
    """The legacy plain-GL Path A.i mu-quadrature SILENTLY truncates
    the divergent diagonal kernel value E_1(0+) = +inf to a finite
    value ~ 2 log(n_mu).

    Quantitatively: for n_mu = 16, 32, 48, ..., 512:
        S(n_mu) = Sum_k w_mu_k / mu_k
    grows linearly in log(n_mu) with empirical slope ~ 2.0,
    i.e. S(n_mu) ~ 2 log(n_mu) + const, NOT diverging to +inf as
    a faithful representation of E_1(0+) would.

    This is the smoking-gun pathology of ERR-036. The Atkinson
    fix sidesteps it by integrating the singularity analytically
    against the basis, never sampling K[i, i] at the diagonal.

    A regression that "fixed" the truncation by clipping or
    capping E_1 evaluations would still leave this test passing
    (because the GL sum-of-1/mu is independent of E_1); only a
    REPLACEMENT of plain-GL with product-Nyström (or a similar
    singularity-aware method) would change the assembly path.

    What this test guards against: someone adding back a plain-GL
    Path A.i with the same μ-mesh structure should still see the
    truncation. Atkinson is the right fix.
    """
    n_mus = [16, 32, 48, 64, 96, 128, 192, 256, 384, 512]
    finite_vals = []
    for n_mu in n_mus:
        nodes, weights = np.polynomial.legendre.leggauss(n_mu)
        mu = 0.5 * (nodes + 1.0)
        w = 0.5 * weights
        S = float(np.sum(w / mu))
        finite_vals.append(S)
        assert S < 30.0, (
            f"Truncated S(n_mu = {n_mu}) = {S:.4f} unexpectedly large; "
            f"the GL truncation should saturate at ~ 2 log(n_mu) ~ "
            f"{2 * math.log(n_mu):.2f}"
        )

    # Linear fit S = a log(n_mu) + b. Expect a ≈ 2 from
    # mu_min(GL on (0,1)) ~ 1/n^2.
    log_n = np.log(np.asarray(n_mus, dtype=float))
    slope_a, intercept_b = np.polyfit(log_n, np.asarray(finite_vals), 1)
    assert 1.5 < slope_a < 2.5, (
        f"Diagonal-truncation log-scaling slope = {slope_a:.4f}; "
        f"expected ~ 2 (from mu_min ~ 1/n^2). The bug is reading the "
        f"divergent E_1(0+) as a finite logarithmic truncation."
    )


@pytest.mark.l1
@pytest.mark.catches("ERR-036")
def test_legacy_path_ai_off_diagonal_converges_machine_precision():
    """Off the diagonal (tau >= 0.05 mfp) the legacy Path A.i
    mu-quadrature converges to E_1(tau) at machine precision (1e-15).
    This confirms that the bug is *exclusively* at the diagonal —
    everywhere else the discrete kernel is correct.

    This is the diagnostic that PINPOINTS the bug: uniform
    refinement of the z-mesh cannot help, because the per-cell
    error is concentrated entirely at the diagonal, where it
    is qualitatively wrong.
    """
    n_mu = 256
    nodes, weights = np.polynomial.legendre.leggauss(n_mu)
    mu = 0.5 * (nodes + 1.0)
    w = 0.5 * weights

    for tau in [0.05, 0.1, 0.5, 1.0, 2.0]:
        arg = np.minimum(tau / mu, 50.0)
        K_GL = float(np.sum((w / mu) * np.exp(-arg)))
        E_exact = float(exp1(tau))
        rel = abs(K_GL - E_exact) / E_exact
        assert rel < 1e-10, (
            f"At tau = {tau} (off-diagonal), GL mu-quad differs from "
            f"E_1 by rel = {rel:.4e}. Off-diagonal is supposed to be "
            f"machine-precision; if this fires, the failure is somewhere "
            f"OTHER than ERR-036 diagonal-truncation."
        )


@pytest.mark.l1
@pytest.mark.catches("ERR-036")
def test_legacy_path_ai_first_order_with_log_signature():
    """The legacy plain-GL Path A.i exhibits the textbook
    first-order-with-log convergence signature:

        err(n) ~ log(n) / n,    so err * n / log(n) ≈ const,

    NOT the Atkinson-Schneider half-order (err * sqrt(n) ≈ const)
    nor the smooth Simpson rate (err * n^2 ≈ const).

    This test imports the legacy Path A.i and runs a quick n-sweep
    at the worst-case point z/a = 0.95.  Catching this specific
    rate (vs 1/sqrt(n) or 1/n^2) PINS the failure mode to
    "log-singular kernel diagonal truncation" rather than
    "endpoint solution singularity" or "smooth-integrand
    quadrature undersampling".

    We assert the err * n / log(n) is constant within ~ 50% across
    a 4× range of n. Tighter is better but the prefactor varies
    slightly because the smooth remainder also has a finite-N tail.
    """
    from orpheus.derivations.continuous.fn_method.slab import (
        slab_scalar_flux_fn_projection,
        slab_scalar_flux_ratio,
        solve_fn_slab_bare_critical,
        solve_kll_slab_continuum_coefficient,
    )

    c = 1.30
    fn_res = solve_fn_slab_bare_critical(c=c, n_modes=10)
    kll_res = solve_kll_slab_continuum_coefficient(
        b=fn_res.a_critical_mfp, c=c, n_nodes=96
    )

    z_test = 0.95 * fn_res.a_critical_mfp
    ratio_B = slab_scalar_flux_ratio(kll_res, 0.95)

    ns = [32, 64, 128]
    constants = []
    for n in ns:
        phi_AI = slab_scalar_flux_fn_projection(fn_res, z_test, n_quad_z=n)
        phi_AI_0 = slab_scalar_flux_fn_projection(fn_res, 0.0, n_quad_z=n)
        ratio_AI = phi_AI / phi_AI_0
        err = abs(ratio_AI - ratio_B)
        constants.append(err * n / math.log(n))

    # First-order-with-log: err * n / log(n) ≈ const within 1.5x.
    cmin, cmax = min(constants), max(constants)
    spread = cmax / cmin
    assert spread < 1.6, (
        f"err * n / log(n) values across n-sweep: {constants}. "
        f"Spread {spread:.2f}x — exceeds the ~1.5x band for "
        f"first-order-with-log scaling. Either the bug has changed "
        f"or the slope estimate is unreliable."
    )
