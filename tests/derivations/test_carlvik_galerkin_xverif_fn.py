"""L1 cross-check: Carlvik-Galerkin vs F_N method at isotropic limit.

Two structurally-independent Pillar-2 methods solving the same physical
problem: bare-critical multiplying slab/sphere with isotropic scattering.

* :mod:`...carlvik_galerkin` — Galerkin spectral expansion of Carlvik's
  integral equation (Dahl-Sjostrand 1979).
* :mod:`...fn_method` — F_N boundary collocation in Case singular
  eigenfunctions (Siewert-Benoist 1979 / Grandjean-Siewert 1979 /
  Siewert-Thomas 1986).

Cross-check: at :math:`\\bar\\mu = 0`, both methods must agree on
the critical secondaries-per-collision ratio :math:`c_{\\rm crit}(d)`.
This is the canonical Pillar-2 redundancy test.

Procedure:

1. Pick a target :math:`c \\in (1, 1.5]` (multiplying medium).
2. Run the F_N solver, which returns the critical half-thickness
   :math:`a_c` (slab) or radius :math:`R_c` (sphere) for which the
   bare medium is critical at this :math:`c`.
3. Set :math:`d = 2 a_c` (slab) or :math:`d = 2 R_c` (sphere) and
   run the Carlvik-Galerkin solver at :math:`(\\bar\\mu = 0, d)`.
4. The Carlvik-Galerkin dominant eigenvalue :math:`c_{\\rm crit}^{CG}`
   must equal the input :math:`c` to within the looser of the two
   methods' truncation tolerances.

Tolerance: ≤ 1e-5 absolute (the brief's target). Empirically the
methods agree at 1e-7 (sphere, weak quadrature dependence) to 1e-5
(slab at high c, where F_N convergence is slower).

This test does NOT prove correctness — both could be wrong. But it
is a strong consistency check across two structurally independent
mathematical frames (Galerkin spectral on integral form vs Case
singular eigenfunctions on the boundary). Combined with each method's
agreement vs published reference tables (DS Tables I, II for
Carlvik-Galerkin; KLL Table I and Sood LA-13511 for F_N), this
provides a structurally independent verification chain.
"""
from __future__ import annotations

import warnings

import pytest

from orpheus.derivations.continuous.carlvik_galerkin.slab import (
    solve_carlvik_galerkin_slab,
)
from orpheus.derivations.continuous.carlvik_galerkin.sphere import (
    solve_carlvik_galerkin_sphere,
)
from orpheus.derivations.continuous.fn_method.slab.one_group import (
    solve_fn_slab_bare_critical,
)
from orpheus.derivations.continuous.fn_method.sphere.one_group import (
    solve_fn_sphere_bare_critical,
)


# ───────────────────────────────────────────────────────────────────
# Slab cross-check
# ───────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.parametrize("c_target", [1.05, 1.20, 1.30, 1.40, 1.50])
def test_l1_slab_xverif_fn_vs_carlvik_galerkin(c_target: float) -> None:
    """Slab μ̄=0: F_N's a_c(c) ↔ Carlvik-Galerkin's c_crit(d=2a_c) agreement.

    Both methods solve the same physics. At isotropic limit, the
    F_N output a_c(c) and the Carlvik-Galerkin c_crit(d=2a_c) must
    cycle back to the same c.
    """
    # Suppress the F_N solver's numpy det() warnings (expected at
    # eigenvalue zero crossings).
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        fn_result = solve_fn_slab_bare_critical(c=c_target, n_modes=10)
    a_c = fn_result.a_critical_mfp
    d_c = 2 * a_c

    cg_result = solve_carlvik_galerkin_slab(
        c=c_target, d=d_c, mu_bar=0.0, n_modes=9, n_quad=128
    )

    abs_err = abs(cg_result.c_critical - c_target)
    assert abs_err <= 1e-5, (
        f"c={c_target}: F_N a_c = {a_c}, CG c_crit at d=2a_c = "
        f"{cg_result.c_critical}, abs error = {abs_err:.2e}"
    )


# ───────────────────────────────────────────────────────────────────
# Sphere cross-check
# ───────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.parametrize("c_target", [1.05, 1.20, 1.30, 1.40, 1.50])
def test_l1_sphere_xverif_fn_vs_carlvik_galerkin(c_target: float) -> None:
    """Sphere μ̄=0: F_N's R_c(c) ↔ Carlvik-Galerkin's c_crit(d=2R_c) agreement."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        fn_result = solve_fn_sphere_bare_critical(c=c_target, n_modes=10)
    R_c = fn_result.R_critical_mfp
    d_c = 2 * R_c

    cg_result = solve_carlvik_galerkin_sphere(
        c=c_target, d=d_c, mu_bar=0.0, n_modes=9, n_quad=128
    )

    abs_err = abs(cg_result.c_critical - c_target)
    assert abs_err <= 1e-5, (
        f"c={c_target}: F_N R_c = {R_c}, CG c_crit at d=2R_c = "
        f"{cg_result.c_critical}, abs error = {abs_err:.2e}"
    )
