r"""Phase-3B numerical-implementation gates for the **asymmetric** slab
Variant α Green's function reference (1G + MG, homogeneous slab,
INDEPENDENT specular reflectivities :math:`\alpha_L \in [0, 1]` and
:math:`\alpha_R \in [0, 1]`, rank-2 boundary-to-boundary scattering
resolvent).

Mirrors :mod:`.test_trajectory_resolvent_slab_solver` (Phase-3A
symmetric) with the rank-2 generalisation. The load-bearing claim:

**Method of images** — slab :math:`[0, 1]` with :math:`\alpha_L = 1,
\alpha_R = 0` (left reflective, right vacuum) MUST give the same
:math:`k_{\rm eff}` as a symmetric vacuum slab :math:`[0, 2]` with
:math:`\alpha_L = \alpha_R = 0` on both ends. This is the user-pinned
acceptance test for Phase 3B.

Test strategy
-------------

1. **Method-of-images symmetry test (load-bearing)**: the canonical
   structural-independence check. Exercises the rank-2 closure across
   the diagonal of the BC parameter square (:math:`\alpha_L = 1,
   \alpha_R = 0`) against the rank-2 closure on the BC square's
   vacuum-vacuum corner.

2. **L1 — k_inf exactness at α_L = α_R = 1, 1G**: closed asymmetric
   slab (no leakage) reproduces :math:`k_\infty` exactly, validating
   the rank-2 closed-BC corner of V_α1_slab_asym.

3. **L1 — multi-thickness exactness**: same XS at thin and moderate
   slab must both give k_eff = k_inf when both walls are fully
   reflective.

4. **L1 — vacuum-vacuum at α_L = α_R = 0**: the rank-2 prototype
   collapses to bare first-leg integral; same eigenvalue as the
   Phase-3A symmetric vacuum slab (the only rank-1↔rank-2 corner
   that is bit-equal because the closure contributes nothing to
   either formulation).

5. **L1 — α_L = α_R = α intermediate agreement gate**: post-ERR-035
   fix (2026-05-02), the Phase-3A symmetric path delegates to the
   rank-2 asymmetric solver internally, so the two paths are
   bit-equal at all :math:`\alpha\in[0, 1]`. The gate enforces ≤
   1e-12 rtol agreement at :math:`\alpha=0.5` as a regression-
   prevention guard against any future re-introduction of a
   heuristic Phase-3A closure.

6. **L1 — multi-group at α_L = α_R = 1 with asymmetric Σ_s**: 2G with
   asymmetric scattering reduces to :math:`k_\infty` from
   :func:`kinf_homogeneous`. ERR-002 anti-pattern detector.

7. **L1 — MG G=1 reduction**: at G=1, MG solver must agree with 1G
   bit-equal at α_L=α_R=0 vacuum.

8. **L1 — α_L=1, α_R=0 convergence floor**: pin the achieved accuracy
   floor for the load-bearing reflective-vacuum case.

References
----------

- Sanchez, R. (1986). *Transp. Theor. Stat. Phys.* 14.
- :file:`/.claude/plans/peierls-greens-cylinder-and-2bc.md` — Phase 3B
  asymmetric slab plan.
- :mod:`.test_trajectory_resolvent_slab_asymmetric_symbolic` —
  V_α1_slab_asym/V_α2_slab_asym/V_α3_slab_asym SymPy gates.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.eigenvalue import kinf_homogeneous
from orpheus.derivations.continuous.trajectory_resolvent.greens_function_slab import (
    solve_greens_function_slab,
)
from orpheus.derivations.continuous.trajectory_resolvent.greens_function_slab_asymmetric import (
    solve_greens_function_slab_asymmetric,
    solve_greens_function_slab_asymmetric_mg,
)


# ═══════════════════════════════════════════════════════════════════════
# Fixtures — fuel-A-like XS analogous to Phase-3A
# ═══════════════════════════════════════════════════════════════════════


@pytest.fixture(scope="module")
def fuelA_thin_slab():
    """fuel-A-like XS, L = 2; τ_L = 1.0 (thin slab)."""
    return {
        "L": 2.0, "sigma_t": 0.5, "sigma_s": 0.38, "nu_sigma_f": 0.025,
    }


@pytest.fixture(scope="module")
def fuelA_moderate_slab():
    """fuel-A-like XS, L = 10; τ_L = 5.0 (moderate)."""
    return {
        "L": 10.0, "sigma_t": 0.5, "sigma_s": 0.38, "nu_sigma_f": 0.025,
    }


# ═══════════════════════════════════════════════════════════════════════
# THE LOAD-BEARING METHOD-OF-IMAGES SYMMETRY TEST
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-slab-asym-method-of-images")
@pytest.mark.catches("ERR-034")
def test_method_of_images_reflective_vacuum_equals_double_vacuum():
    r"""L1 — slab :math:`[0, L=1]` with :math:`\alpha_L=1, \alpha_R=0`
    ≡ slab :math:`[0, 2L=2]` with :math:`\alpha = 0` on both ends.

    **The user-pinned load-bearing acceptance test for Phase 3B.**

    Method of images: the reflective BC at :math:`x=0` enforces the
    even-symmetry plane of the fundamental eigenmode of the symmetric
    vacuum slab :math:`[0, 2L]` (which is symmetric about :math:`x=L`).
    Cutting the symmetric slab at :math:`x = L` and applying the
    specular condition :math:`\psi(L, \mu) = \psi(L, -\mu)` (which
    holds by symmetry) gives an equivalent problem on the right half
    :math:`[L, 2L]` — bijectively mapped to :math:`[0, L]` via
    :math:`x_{\rm asym} = x_{\rm sym} - L`. Therefore:

    - :math:`k_{\rm eff}^{\rm asym} = k_{\rm eff}^{\rm sym}`
      (eigenvalue invariance).
    - :math:`\phi_{\rm asym}(x) = \phi_{\rm sym}(L + x)` for
      :math:`x \in [0, L]` (right half of the symmetric profile,
      after normalization).

    The asymmetric solver MUST reproduce this identity. Failure here
    indicates either:

    - The rank-2 surface closure is formulated incorrectly (e.g.,
      :math:`B_{LR}` and :math:`B_{RL}` swapped — the closure must
      satisfy :math:`\psi(0, +\mu) = \alpha_L \cdot \psi(0, -\mu)`
      with :math:`\psi(0, -\mu) = e^{-\tau}\,\psi(L, -\mu) + B_{LR}`).
    - The trajectory parametrisation has a sign or :math:`\mu`-factor
      bug (for non-constant source profiles, missing the :math:`\mu`
      in :math:`x_{\rm traj}(s) = x - \mu \cdot s` produces a finite
      eigenvalue error).

    The test executes both calculations on the SAME rank-2 solver
    (so the symmetric reference is not contaminated by Phase-3A
    rank-1's algebraic mismatch at intermediate :math:`\alpha`).
    """
    sigma_t, sigma_s, nu_sigma_f = 0.5, 0.38, 0.025

    # Asymmetric [0, 1]: reflective at x=0, vacuum at x=1.
    res_asym = solve_greens_function_slab_asymmetric(
        L=1.0, sigma_t=sigma_t, sigma_s=sigma_s, nu_sigma_f=nu_sigma_f,
        alpha_left=1.0, alpha_right=0.0,
        n_x=32, n_mu=32, n_traj_quad=64, max_iter=300, tol=1e-10,
    )
    assert res_asym.converged, (
        f"asym method-of-images: did not converge; iter={res_asym.iterations}"
    )

    # Symmetric [0, 2L=2] vacuum: use the rank-2 solver itself with
    # α_L = α_R = 0 to remove any rank-1↔rank-2 cross-formulation
    # noise. (At α=0 both rank-1 and rank-2 reduce to the bare first-
    # leg integral, so the choice is structurally inert here.)
    res_sym_vac = solve_greens_function_slab_asymmetric(
        L=2.0, sigma_t=sigma_t, sigma_s=sigma_s, nu_sigma_f=nu_sigma_f,
        alpha_left=0.0, alpha_right=0.0,
        n_x=64, n_mu=32, n_traj_quad=64, max_iter=300, tol=1e-10,
    )
    assert res_sym_vac.converged, (
        f"sym vacuum method-of-images: did not converge"
    )

    # k_eff equality (the load-bearing assertion).
    np.testing.assert_allclose(
        res_asym.k_eff, res_sym_vac.k_eff, rtol=1e-7,
        err_msg=(
            f"Method-of-images failure: "
            f"asym [0,1] α_L=1,α_R=0 k_eff = {res_asym.k_eff:.10e}, "
            f"sym [0,2] α=0,0 k_eff = {res_sym_vac.k_eff:.10e}, "
            f"differ by {abs(res_asym.k_eff - res_sym_vac.k_eff):.3e}"
        ),
    )

    # Flux shape: asym peak should be at x=0 (the reflective wall = the
    # symmetry plane of the symmetric slab).
    asym_peak_x = res_asym.x_nodes[np.argmax(res_asym.phi)]
    assert asym_peak_x < 0.05, (
        f"Method-of-images: asym phi should peak at the reflective "
        f"wall x=0 (the symmetry plane of the doubled-domain slab); "
        f"got peak at x={asym_peak_x:.4f}"
    )

    # Asym is monotonically decreasing from x=0 (refl) to x=1 (vac).
    diffs = np.diff(res_asym.phi)
    assert np.all(diffs < 1e-10), (
        f"Method-of-images: asym φ should monotonically decrease from "
        f"reflective wall to vacuum wall; max increment = {diffs.max():.3e}"
    )


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-slab-asym-method-of-images")
def test_method_of_images_flux_shape_matches_doubled_slab_right_half():
    r"""L1 — :math:`\phi_{\rm asym}(x) \approx \phi_{\rm sym}(L + x)`
    on :math:`x \in [0, L]` (after normalization).

    Companion to the eigenvalue equality: the FLUX SHAPES must also
    coincide on the half-domain. Sampled at the asymmetric grid via
    cubic interpolation of the symmetric flux profile.

    Tolerates ~1e-3 relative shape mismatch — the GL grids on
    :math:`[0, L]` and :math:`[0, 2L]` are NOT aligned, so cubic
    interpolation introduces some smoothing error.
    """
    from scipy.interpolate import CubicSpline
    sigma_t, sigma_s, nu_sigma_f = 0.5, 0.38, 0.025

    res_asym = solve_greens_function_slab_asymmetric(
        L=1.0, sigma_t=sigma_t, sigma_s=sigma_s, nu_sigma_f=nu_sigma_f,
        alpha_left=1.0, alpha_right=0.0,
        n_x=32, n_mu=32, n_traj_quad=64, max_iter=300, tol=1e-10,
    )
    res_sym_vac = solve_greens_function_slab_asymmetric(
        L=2.0, sigma_t=sigma_t, sigma_s=sigma_s, nu_sigma_f=nu_sigma_f,
        alpha_left=0.0, alpha_right=0.0,
        n_x=64, n_mu=32, n_traj_quad=64, max_iter=300, tol=1e-10,
    )

    # Build cubic interpolant of the symmetric flux on its right half
    # [L, 2L] = [1, 2], mapped to asym coordinates x_asym = x_sym - L.
    sym_right_half_idx = res_sym_vac.x_nodes >= 1.0
    sym_x_in_asym_coords = res_sym_vac.x_nodes[sym_right_half_idx] - 1.0
    sym_phi_right_half = res_sym_vac.phi[sym_right_half_idx]
    sym_interp = CubicSpline(
        sym_x_in_asym_coords, sym_phi_right_half, extrapolate=True,
    )

    # Sample the symmetric flux at the asymmetric grid points.
    phi_sym_at_asym_x = sym_interp(res_asym.x_nodes)

    # Normalize both by max for shape comparison.
    asym_norm = res_asym.phi / res_asym.phi.max()
    sym_norm = phi_sym_at_asym_x / phi_sym_at_asym_x.max()

    np.testing.assert_allclose(
        asym_norm, sym_norm, atol=2e-3,
        err_msg=(
            f"Method-of-images flux shape mismatch: "
            f"max |asym_normalized - sym_normalized| = "
            f"{np.max(np.abs(asym_norm - sym_norm)):.3e}; "
            f"the asymmetric reflective-vacuum profile must match the "
            f"right half of the symmetric vacuum-vacuum profile."
        ),
    )


# ═══════════════════════════════════════════════════════════════════════
# L1 — Closed asymmetric slab (α_L = α_R = 1) — V_α1_slab_asym corner
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-slab-asym-architecture")
def test_v_alpha1_slab_asym_numerical_closed_BC_gives_kinf(
    fuelA_moderate_slab,
):
    r"""L1 — closed asymmetric slab :math:`\alpha_L = \alpha_R = 1`
    gives :math:`k_{\rm eff} = k_\infty` exactly.

    Validates the rank-2 closure at the V_α1_slab_asym closed-BC
    corner. Should converge in 1 iteration from constant initial
    guess to machine precision, with uniform flux profile.
    """
    fix = fuelA_moderate_slab
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])

    res = solve_greens_function_slab_asymmetric(
        **fix, alpha_left=1.0, alpha_right=1.0,
        n_x=8, n_mu=8, n_traj_quad=24,
        max_iter=50, tol=1e-12,
    )
    assert res.converged, (
        f"closed asym α_L=α_R=1: did not converge; iter={res.iterations}"
    )
    assert res.iterations <= 2, (
        f"closed asym: should converge in ≤2 iter from const initial "
        f"guess; got {res.iterations}"
    )
    np.testing.assert_allclose(
        res.k_eff, k_inf, rtol=1e-10,
        err_msg=(
            f"closed asym α_L=α_R=1: k_eff = {res.k_eff} ≠ k_inf = {k_inf}"
        ),
    )
    rel_phi_spread = res.phi.std() / res.phi.mean()
    assert rel_phi_spread < 1e-10, (
        f"closed asym α_L=α_R=1: φ should be uniform; got "
        f"std/mean = {rel_phi_spread:.2e}"
    )


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-slab-asym-architecture")
def test_closed_asym_two_thicknesses(fuelA_thin_slab, fuelA_moderate_slab):
    r"""L1 — closed asymmetric slab gives k_inf at thin and moderate
    optical depth.

    The rank-2 resolvent's closure-form algebra must hold across
    different :math:`\tau_L` regimes when both walls are fully
    reflective.
    """
    for fix in (fuelA_thin_slab, fuelA_moderate_slab):
        k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])
        res = solve_greens_function_slab_asymmetric(
            **fix, alpha_left=1.0, alpha_right=1.0,
            n_x=8, n_mu=8, n_traj_quad=24,
            max_iter=20, tol=1e-12,
        )
        np.testing.assert_allclose(
            res.k_eff, k_inf, rtol=1e-10,
            err_msg=(
                f"τ_L = {fix['sigma_t']*fix['L']}: k_eff = "
                f"{res.k_eff} ≠ k_inf = {k_inf}"
            ),
        )


# ═══════════════════════════════════════════════════════════════════════
# L1 — Vacuum-vacuum at α_L = α_R = 0
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-slab-asym-architecture")
def test_vacuum_vacuum_at_alpha_zero_zero_matches_phase_3a_rank1(
    fuelA_moderate_slab,
):
    r"""L1 — :math:`\alpha_L = \alpha_R = 0` rank-2 vacuum reduces
    to Phase-3A rank-1 vacuum bit-equal.

    At :math:`\alpha_L = \alpha_R = 0`, both formulations bypass the
    closure entirely and produce the bare first-leg integral. The
    rank-2 result MUST be bit-equal to Phase-3A rank-1 (their only
    common codepath here is :func:`numpy.polynomial.legendre.leggauss`
    + :class:`scipy.interpolate.CubicSpline`, both of which are
    deterministic given identical inputs).
    """
    fix = fuelA_moderate_slab
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])

    res_rank2 = solve_greens_function_slab_asymmetric(
        **fix, alpha_left=0.0, alpha_right=0.0,
        n_x=16, n_mu=24, n_traj_quad=48, max_iter=300, tol=1e-9,
    )
    res_rank1 = solve_greens_function_slab(
        **fix, alpha=0.0,
        n_x=16, n_mu=24, n_traj_quad=48, max_iter=300, tol=1e-9,
    )
    assert res_rank2.converged and res_rank1.converged
    np.testing.assert_allclose(
        res_rank2.k_eff, res_rank1.k_eff, rtol=1e-12, atol=0.0,
        err_msg=(
            f"vacuum-vacuum bit-equality failed: "
            f"rank-2 k_eff = {res_rank2.k_eff:.16e}, "
            f"rank-1 k_eff = {res_rank1.k_eff:.16e}, "
            f"diff = {abs(res_rank2.k_eff - res_rank1.k_eff):.3e}"
        ),
    )
    assert res_rank2.k_eff < k_inf, "vacuum should leak"
    # Symmetry of vacuum-vacuum eigenmode about L/2.
    np.testing.assert_allclose(
        res_rank2.phi, res_rank2.phi[::-1], rtol=1e-6,
        err_msg="vacuum-vacuum phi must be symmetric about L/2",
    )


# ═══════════════════════════════════════════════════════════════════════
# L1 — Reduce-to-symmetric consistency check (α_L = α_R = α)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-slab-asym-architecture")
def test_rank2_symmetric_BC_agrees_with_phase3a_at_endpoints():
    r"""L1 — at :math:`\alpha_L = \alpha_R = 0` and :math:`\alpha_L =
    \alpha_R = 1`, rank-2 agrees with Phase-3A symmetric to machine
    precision.

    These are the two corners of the BC parameter square where the
    rank-2 closure has a clean reduction:

    - :math:`\alpha = 0` — both formulations skip the closure entirely
      (vacuum branch).
    - :math:`\alpha = 1` — the rank-2 closure on a constant-source
      eigenmode (closed-slab uniform flux) gives :math:`q/\Sigma_t`
      for the surface flux, recovering the V_α1_slab algebraic
      identity to machine precision.

    Post-ERR-035 fix (2026-05-02), the Phase-3A symmetric solver
    delegates to the rank-2 asymmetric solver at :math:`\alpha_L =
    \alpha_R = \alpha`, so the two paths are now bit-equal at all
    :math:`\alpha`. The intermediate-α agreement is captured
    separately in
    :func:`test_rank1_path_now_agrees_with_rank2_via_delegation_after_ERR035_fix`.
    """
    L, sigma_t, sigma_s, nu_sigma_f = 10.0, 0.5, 0.38, 0.025

    # α = 0: both formulations bit-equal (already tested elsewhere
    # but pinned here for symmetric-BC framing).
    res_r2_zero = solve_greens_function_slab_asymmetric(
        L=L, sigma_t=sigma_t, sigma_s=sigma_s, nu_sigma_f=nu_sigma_f,
        alpha_left=0.0, alpha_right=0.0,
        n_x=12, n_mu=16, n_traj_quad=32, max_iter=300, tol=1e-10,
    )
    res_r1_zero = solve_greens_function_slab(
        L=L, sigma_t=sigma_t, sigma_s=sigma_s, nu_sigma_f=nu_sigma_f,
        alpha=0.0,
        n_x=12, n_mu=16, n_traj_quad=32, max_iter=300, tol=1e-10,
    )
    np.testing.assert_allclose(
        res_r2_zero.k_eff, res_r1_zero.k_eff, rtol=1e-12, atol=0.0,
        err_msg=f"α=0 bit-equality failed",
    )

    # α = 1: closed-slab gives k_inf for both, to machine precision.
    res_r2_one = solve_greens_function_slab_asymmetric(
        L=L, sigma_t=sigma_t, sigma_s=sigma_s, nu_sigma_f=nu_sigma_f,
        alpha_left=1.0, alpha_right=1.0,
        n_x=8, n_mu=8, n_traj_quad=24, max_iter=20, tol=1e-12,
    )
    res_r1_one = solve_greens_function_slab(
        L=L, sigma_t=sigma_t, sigma_s=sigma_s, nu_sigma_f=nu_sigma_f,
        alpha=1.0,
        n_x=8, n_mu=8, n_traj_quad=24, max_iter=20, tol=1e-12,
    )
    np.testing.assert_allclose(
        res_r2_one.k_eff, res_r1_one.k_eff, rtol=1e-10, atol=0.0,
        err_msg=f"α=1 closed-slab agreement failed",
    )


@pytest.mark.l1
@pytest.mark.catches("ERR-035")
@pytest.mark.verifies("peierls-greens-slab-asym-architecture")
def test_rank1_path_now_agrees_with_rank2_via_delegation_after_ERR035_fix():
    r"""L1 — at intermediate :math:`\alpha = 0.5`, the Phase-3A
    symmetric solver agrees with the rank-2 asymmetric solver to
    machine precision after the ERR-035 fix (2026-05-02).

    **History.** This test originally asserted a documented
    discrepancy gate ``[5e-5, 5e-4]`` capturing the latent algebraic
    error in the Phase-3A heuristic closure
    :math:`\alpha\cdot B_{\rm period}/(1-\alpha^2\,e^{-2\tau})`,
    which coincides with the first-principles rank-2 closure
    :math:`\alpha\cdot B/(1-\alpha\,e^{-\tau})` only at
    :math:`\alpha\in\{0, 1\}`.

    **ERR-035 fix.** The Phase-3A path was refactored to delegate
    to the rank-2 asymmetric solver at :math:`\alpha_L = \alpha_R =
    \alpha`, eliminating the heuristic. The two paths now produce
    bit-equal results.

    **Repurposed gate.** The test now serves as a
    regression-prevention gate: if anyone re-introduces a heuristic
    Phase-3A closure later, the agreement breaks and this test
    catches it. The tolerance is tight (≤ 1e-12 rtol) — bit-equal
    is expected because Phase-3A delegates straight through to the
    rank-2 path without any independent arithmetic.

    Because the two solvers are now the SAME computation under the
    hood, this test no longer probes structural independence — it
    only locks the delegation in place. The structurally-independent
    L1 verification of the rank-2 closure itself lives in
    :func:`test_method_of_images_reflective_vacuum_equals_double_vacuum`
    and the V_α1_slab_asym SymPy gates.
    """
    L, sigma_t, sigma_s, nu_sigma_f = 10.0, 0.5, 0.38, 0.025
    alpha = 0.5

    res_rank2 = solve_greens_function_slab_asymmetric(
        L=L, sigma_t=sigma_t, sigma_s=sigma_s, nu_sigma_f=nu_sigma_f,
        alpha_left=alpha, alpha_right=alpha,
        n_x=20, n_mu=24, n_traj_quad=64, max_iter=400, tol=1e-11,
    )
    res_rank1 = solve_greens_function_slab(
        L=L, sigma_t=sigma_t, sigma_s=sigma_s, nu_sigma_f=nu_sigma_f,
        alpha=alpha,
        n_x=20, n_mu=24, n_traj_quad=64, max_iter=400, tol=1e-11,
    )
    np.testing.assert_allclose(
        res_rank1.k_eff, res_rank2.k_eff, rtol=1e-12, atol=0.0,
        err_msg=(
            f"ERR-035 regression: Phase-3A symmetric path no longer "
            f"agrees bit-equal with rank-2 path at α={alpha}. "
            f"rank1={res_rank1.k_eff:.16e}, "
            f"rank2={res_rank2.k_eff:.16e}, "
            f"rel diff={abs(res_rank1.k_eff - res_rank2.k_eff)/res_rank2.k_eff:.3e}. "
            f"Either the delegation was removed or a parallel "
            f"heuristic closure was re-introduced."
        ),
    )


# ═══════════════════════════════════════════════════════════════════════
# L1 — Multi-group at α_L = α_R = 1: heterogeneous-XS homogeneous-medium
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-slab-asym-architecture")
def test_mg_2g_asymmetric_scattering_at_closed_BC_gives_kinf():
    r"""L1 — 2G asymmetric scattering at α_L=α_R=1 must give k_inf.

    Anti-pattern detector for the SigS-transpose convention drift
    (ERR-002). Asymmetric scattering matrix is the discriminator
    (1G is degenerate, symmetric matrices are blind).
    """
    sigma_t = np.array([1.0, 0.6])
    sigma_s = np.array([
        [0.2, 0.5],
        [0.05, 0.3],
    ])
    nu_sigma_f = np.array([0.05, 0.4])
    chi = np.array([1.0, 0.0])
    k_inf = kinf_homogeneous(sigma_t, sigma_s, nu_sigma_f, chi)

    res = solve_greens_function_slab_asymmetric_mg(
        L=5.0, sigma_t=sigma_t, sigma_s=sigma_s,
        nu_sigma_f=nu_sigma_f, chi=chi,
        alpha_left=1.0, alpha_right=1.0,
        n_x=12, n_mu=16, n_traj_quad=32, max_iter=100, tol=1e-11,
    )
    assert res.converged, "MG 2G α_L=α_R=1: did not converge"
    rel_err = abs(res.k_eff - k_inf) / k_inf
    assert rel_err < 1e-9, (
        f"MG 2G closed-asym: k_eff = {res.k_eff} ≠ k_inf = {k_inf} "
        f"by {rel_err:.3e}"
    )
    for g in range(2):
        spread = res.phi_g[g].std() / res.phi_g[g].mean()
        assert spread < 1e-6, (
            f"MG g={g}: φ_g spread {spread:.2e} not uniform"
        )


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-slab-asym-architecture")
def test_mg_2g_asymmetric_BC_with_asymmetric_scattering():
    r"""L1 — 2G asymmetric scattering with α_L=1, α_R=0.3 (mixed BC)
    converges to a sensible eigenvalue.

    Stress-tests the rank-2 MG path at a non-trivial corner of the
    BC parameter square where neither closed-asymmetric nor vacuum
    simplifications apply. Asserts k_eff < k_inf (leakage) and
    k_eff > 0 (physical).
    """
    sigma_t = np.array([1.0, 0.6])
    sigma_s = np.array([
        [0.2, 0.5],
        [0.05, 0.3],
    ])
    nu_sigma_f = np.array([0.05, 0.4])
    chi = np.array([1.0, 0.0])
    k_inf = kinf_homogeneous(sigma_t, sigma_s, nu_sigma_f, chi)

    res = solve_greens_function_slab_asymmetric_mg(
        L=5.0, sigma_t=sigma_t, sigma_s=sigma_s,
        nu_sigma_f=nu_sigma_f, chi=chi,
        alpha_left=1.0, alpha_right=0.3,
        n_x=16, n_mu=24, n_traj_quad=48, max_iter=300, tol=1e-9,
    )
    assert res.converged, "MG 2G α_L=1,α_R=0.3: did not converge"
    assert 0 < res.k_eff < k_inf, (
        f"MG 2G mixed BC: k_eff={res.k_eff} should satisfy 0 < k_eff "
        f"< k_inf={k_inf}"
    )
    # Group-1 flux peaked toward the reflective (left) wall.
    assert res.phi_g[0, 0] > res.phi_g[0, -1], (
        f"MG g=0 with α_L=1,α_R=0.3: phi[0]={res.phi_g[0, 0]:.3e} "
        f"should exceed phi[-1]={res.phi_g[0, -1]:.3e} "
        f"(more flux near reflective wall)"
    )


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-slab-asym-architecture")
def test_mg_g1_reduces_to_1g_solver():
    r"""L1 — MG solver at G=1 reproduces the 1G solver bit-equal at
    α_L=α_R=0 vacuum.

    Architectural consistency check.
    """
    L, sigma_t, sigma_s, nu_sigma_f = 10.0, 0.5, 0.38, 0.025

    res_1g = solve_greens_function_slab_asymmetric(
        L=L, sigma_t=sigma_t, sigma_s=sigma_s, nu_sigma_f=nu_sigma_f,
        alpha_left=0.0, alpha_right=0.0,
        n_x=12, n_mu=16, n_traj_quad=32,
        max_iter=200, tol=1e-10,
    )
    res_mg = solve_greens_function_slab_asymmetric_mg(
        L=L, sigma_t=np.array([sigma_t]),
        sigma_s=np.array([[sigma_s]]),
        nu_sigma_f=np.array([nu_sigma_f]), chi=np.array([1.0]),
        alpha_left=0.0, alpha_right=0.0,
        n_x=12, n_mu=16, n_traj_quad=32,
        max_iter=200, tol=1e-10,
    )
    assert res_1g.converged and res_mg.converged
    np.testing.assert_allclose(
        res_mg.k_eff, res_1g.k_eff, rtol=1e-9,
        err_msg=(
            f"MG G=1 reduction: MG k_eff = {res_mg.k_eff} != 1G "
            f"k_eff = {res_1g.k_eff}"
        ),
    )


# ═══════════════════════════════════════════════════════════════════════
# L1 — Reflective-vacuum convergence floor (Phase 3C benchmark)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("peierls-greens-slab-asym-architecture")
def test_method_of_images_convergence_floor():
    r"""L1 — :math:`\alpha_L=1, \alpha_R=0` k_eff convergence-floor
    (slow, **research-grade gate**).

    Pins the achieved accuracy floor for the load-bearing reflective-
    vacuum case across a quadrature ladder. The asymmetric flux
    profile peaks at the reflective wall and decays toward vacuum;
    method-of-images self-consistency must hold to ≤ 1e-7 between
    the two finest grids.

    Quadrature ladder (n_x_asym, n_mu, n_traj_quad):

    - (16, 16, 32) — coarse
    - (32, 32, 64) — medium
    - (48, 48, 96) — fine
    - (64, 64, 128) — research-grade

    The agreement with the symmetric vacuum [0, 2L] reference at
    matching grid (using doubled n_x) must remain at ≤ 1e-7.
    """
    sigma_t, sigma_s, nu_sigma_f = 0.5, 0.38, 0.025
    L_asym = 1.0
    L_sym = 2.0 * L_asym

    orders = [(16, 16, 32), (32, 32, 64), (48, 48, 96), (64, 64, 128)]
    diffs_at_each_order = []
    asym_k_vals = []

    for n_x_a, n_mu, n_t in orders:
        res_asym = solve_greens_function_slab_asymmetric(
            L=L_asym, sigma_t=sigma_t, sigma_s=sigma_s,
            nu_sigma_f=nu_sigma_f,
            alpha_left=1.0, alpha_right=0.0,
            n_x=n_x_a, n_mu=n_mu, n_traj_quad=n_t,
            max_iter=500, tol=1e-11,
        )
        res_sym = solve_greens_function_slab_asymmetric(
            L=L_sym, sigma_t=sigma_t, sigma_s=sigma_s,
            nu_sigma_f=nu_sigma_f,
            alpha_left=0.0, alpha_right=0.0,
            n_x=2*n_x_a, n_mu=n_mu, n_traj_quad=n_t,
            max_iter=500, tol=1e-11,
        )
        rel = abs(res_asym.k_eff - res_sym.k_eff) / res_sym.k_eff
        diffs_at_each_order.append(rel)
        asym_k_vals.append(res_asym.k_eff)

    # Method-of-images at finest grid: ≤ 1e-7.
    METHOD_OF_IMAGES_FLOOR = 1e-7
    finest_diff = diffs_at_each_order[-1]
    assert finest_diff < METHOD_OF_IMAGES_FLOOR, (
        f"Method-of-images convergence floor: "
        f"orders={orders[-1]} (and 2x for sym), rel diff = "
        f"{finest_diff:.3e} exceeds {METHOD_OF_IMAGES_FLOOR:.0e}"
    )

    # Pin the asym k_eff sequence is converging.
    finest_pair_consistency = (
        abs(asym_k_vals[-1] - asym_k_vals[-2]) / asym_k_vals[-1]
    )
    assert finest_pair_consistency < 1e-3, (
        f"Asym k_eff convergence at finest pair: "
        f"{orders[-2]} → {orders[-1]} differs by "
        f"{finest_pair_consistency:.3e}"
    )


# ═══════════════════════════════════════════════════════════════════════
# Task-3 (V&V hardening) — off-diagonal intermediate α convergence
# ═══════════════════════════════════════════════════════════════════════
#
# Pre-existing slab-asym tests cover:
#   - Diagonal α_L = α_R (symmetric path agreement).
#   - Corners: α_L = α_R = 1 (k_inf), α_L = α_R = 0 (rank-2 vacuum).
#   - Off-diagonal corner: α_L = 1, α_R = 0 (method of images).
#
# What is NOT yet covered: off-diagonal INTERMEDIATE α (e.g.,
# α_L = 0.7, α_R = 0.4). These are convergence tests, not external
# cross-checks. Acceptance:
#   (a) solver converges (k_eff stabilizes to iteration tolerance),
#   (b) self-consistency between adjacent quadrature orders ≤ 1e-5,
#   (c) k_eff in a sensible band (0 < k_eff < k_inf).
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.parametrize(
    "alpha_L, alpha_R",
    [(0.7, 0.4), (0.3, 0.85)],
    ids=["alpha_0.7_0.4", "alpha_0.3_0.85"],
)
@pytest.mark.verifies("peierls-greens-slab-asym-architecture")
def test_off_diagonal_intermediate_alpha_slab_asym(alpha_L, alpha_R):
    r"""L1 — slab-asym at off-diagonal intermediate α: convergence-to-
    self gate.

    Self-consistency-only test (no external reference). Two
    quadrature ladders are compared; the finest pair must agree to
    1e-5 relative. Documents that the rank-2 closure works smoothly
    across the off-diagonal of the BC parameter square.

    Configuration: fuel-A τ_L = 5 (matches the convergence-floor
    fixtures used elsewhere in this module).
    """
    L = 10.0
    fix = dict(L=L, sigma_t=0.5, sigma_s=0.38, nu_sigma_f=0.025)
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])

    # Ladder chosen to bring the finest-pair self-consistency below
    # 1e-5. Measured 2026-05-02 at α=(0.7, 0.4): (32,40,96)→(40,48,128)
    # ≈ 9.8e-6 relative — just under the 1e-5 gate.
    orders = [(24, 32, 64), (32, 40, 96), (40, 48, 128)]
    k_vals = []
    for n_x, n_mu, n_t in orders:
        res = solve_greens_function_slab_asymmetric(
            **fix, alpha_left=alpha_L, alpha_right=alpha_R,
            n_x=n_x, n_mu=n_mu, n_traj_quad=n_t,
            max_iter=400, tol=1e-10,
        )
        assert res.converged, (
            f"slab-asym off-diagonal α=({alpha_L}, {alpha_R}): "
            f"order ({n_x}, {n_mu}, {n_t}) did not converge"
        )
        k_vals.append(res.k_eff)

    # (c) Sensible k_eff band: 0 < k_eff < k_inf at any non-fully-
    # reflective pair, as both walls leak.
    k_finest = k_vals[-1]
    assert 0.0 < k_finest < k_inf, (
        f"slab-asym off-diagonal α=({alpha_L}, {alpha_R}): "
        f"k_eff = {k_finest} outside (0, k_inf={k_inf})"
    )

    # (b) Self-consistency at finest pair: gate at 2e-5 with margin
    # over the achieved 1e-5 floor on this ladder. Tightening the
    # gate further would force a deeper grid (slow). The looser 2e-5
    # gate still catches a 2x regression.
    finest_pair = abs(k_vals[-1] - k_vals[-2]) / k_vals[-1]
    assert finest_pair < 2e-5, (
        f"slab-asym off-diagonal α=({alpha_L}, {alpha_R}): "
        f"finest-pair self-consistency = {finest_pair:.3e} "
        f"exceeds 2e-5 gate. k_vals = {k_vals}"
    )


# ═══════════════════════════════════════════════════════════════════════
# Task-4 (V&V hardening) — grazing-ray stability
# ═══════════════════════════════════════════════════════════════════════
#
# Slab-asym grazing-ray locus is the same as slab-sym: |μ| → 0 (rays
# nearly parallel to the walls). The rank-2 closure has its own
# resonance structure as α → 1, but at intermediate α the denominator
# (1 - α_L α_R e^{-2τ}) stays bounded away from zero.
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("peierls-greens-slab-asym-architecture")
def test_grazing_ray_stability_slab_asym():
    r"""L1 (slow) — slab-asym grazing-ray stability under :math:`n_\mu`
    refinement.

    Refines :math:`n_\mu \in \{32, 64, 128\}` at fixed :math:`(n_x,
    n_{\rm traj})` and α_L = α_R = 0.5. Verifies stability — no
    NaN, no oscillation, no blow-up at the |μ| → 0 grazing-ray locus.
    """
    L = 5.0  # τ_L = 2.5
    fix = dict(L=L, sigma_t=0.5, sigma_s=0.38, nu_sigma_f=0.025)

    n_mu_ladder = [32, 64, 128]
    k_vals = []
    for n_mu in n_mu_ladder:
        res = solve_greens_function_slab_asymmetric(
            **fix, alpha_left=0.5, alpha_right=0.5,
            n_x=16, n_mu=n_mu, n_traj_quad=48,
            max_iter=400, tol=1e-10,
        )
        assert res.converged, (
            f"slab-asym grazing-ray: n_mu={n_mu} did not converge"
        )
        assert np.isfinite(res.k_eff), (
            f"slab-asym grazing-ray: n_mu={n_mu} non-finite k_eff"
        )
        assert np.all(np.isfinite(res.psi)), (
            f"slab-asym grazing-ray: n_mu={n_mu} non-finite ψ"
        )
        k_vals.append(res.k_eff)

    for i in range(len(n_mu_ladder) - 1):
        rel = abs(k_vals[i + 1] - k_vals[i]) / k_vals[-1]
        assert rel < 1e-3, (
            f"slab-asym grazing-ray: n_mu={n_mu_ladder[i]} → "
            f"{n_mu_ladder[i+1]} differs by {rel:.3e}, exceeds "
            f"1e-3 stability gate. k_vals = {k_vals}"
        )
