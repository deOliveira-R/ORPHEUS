r"""Phase-1 numerical-implementation gates for the cylinder Variant α
Green's function reference (1G + MG, homogeneous infinite cylinder).

Mirrors :mod:`.test_peierls_greens_function_solver` (sphere) and
:mod:`.test_peierls_greens_function_xverif` (sphere cross-verification)
with cylinder geometry. The load-bearing claim:

**Variant α gives the exact :math:`k_{\rm eff} = k_\infty` for the closed
homogeneous cylinder with specular BC** (V_α1_cyl algebraic identity).
This is the cylinder analogue of the sphere's Variant α exactness
property and the canonical Phase-1 acceptance test.

Test strategy
-------------

1. **L1 — k_inf exactness at α=1, 1G**: closed cylinder with perfect
   specular has no leakage and rank-1 isotropic eigenmode with
   :math:`k_{\rm eff} = k_\infty = \nu\Sigma_f/\Sigma_a` exactly. Must
   reproduce to ≤ 1e-10 within iteration tolerance.

2. **L1 — multi-thickness exactness**: same XS at thin (τ_R = 2.5) and
   moderate (τ_R = 5) cylinder must both yield k_eff = k_inf — the
   bounce-sum geometric series is wired correctly across τ_R range.

3. **L1 — non-uniform initial guess at α=1**: power iteration starting
   from a strongly non-uniform initial state must converge to the rank-1
   isotropic eigenmode and recover k_inf within the iteration
   tolerance.

4. **L1 — multi-group asymmetric scattering**: a 2G case with
   asymmetric :math:`\Sigma_s` (:math:`\Sigma_{s,01} \ne \Sigma_{s,10}`)
   must reduce to :math:`k_\infty` from
   :func:`orpheus.derivations.common.eigenvalue.kinf_homogeneous` at
   :math:`\alpha = 1`. Asymmetric XS detect the
   ``SigS^T`` vs ``SigS`` convention drift (ERR-002 anti-pattern).

5. **L1 — vacuum BC (α=0) sanity**: at α=0 the spatial flux profile
   must be peaked at center and depleted at the surface (positive
   leakage); :math:`k_{\rm eff} < k_\infty`. This stress-tests the
   trajectory machinery — V_α1_cyl's algebraic cancellation does NOT
   protect this case (the surface fixed-point closure is bypassed).

6. **L1 — α=0 cross-check between MG and 1G reductions**: at G=1 the
   MG solver must agree with the 1G solver bit-equal.

Anti-cross-check prohibition (Issue #129)
------------------------------------------

We do NOT cross-check cylinder Variant α against:
- Thin-slab limit (cylinder is 3D angle-resolved; planar Ki_n
  pre-integrates).
- Cylinder Nyström vacuum (Issue #129: axial pre-integration introduces
  a documented ~1e-4 mismatch — both formulations are correct but
  describe slightly different approximations of the same equation).

The L1 cross-check is **k_inf-exactness at α=1**, which is an
analytical invariant from V_α1_cyl, structurally independent of any
other ORPHEUS solver.

Predecessor / sibling tests
---------------------------

- :mod:`.test_peierls_greens_function_cylinder_symbolic` — V_α1_cyl,
  V_α2_cyl, V_α3_cyl SymPy gates.
- :mod:`.test_peierls_greens_function_solver` — sphere V_α1.numerical
  (the structural template this test file mirrors).

References
----------

- Sanchez, R. (1986). *Transp. Theor. Stat. Phys.* 14.
- :file:`/.claude/plans/peierls-greens-cylinder-and-2bc.md` — Phase 1
  cylinder plan.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.eigenvalue import kinf_homogeneous
from orpheus.derivations.continuous.peierls_nystrom.geometry import (
    compute_P_ss_cylinder,
    compute_T_specular_cylinder_3d,
)
from orpheus.derivations.continuous.peierls_greens_function.greens_function_cylinder import (
    solve_greens_function_cylinder,
    solve_greens_function_cylinder_mg,
)


# ═══════════════════════════════════════════════════════════════════════
# Fixtures — matching the sphere fuel-A-like XS
# ═══════════════════════════════════════════════════════════════════════


@pytest.fixture(scope="module")
def fuelA_thin_cylinder():
    """fuel-A-like XS, R = 5; τ_R = 2.5."""
    return {
        "R": 5.0, "sigma_t": 0.5, "sigma_s": 0.38, "nu_sigma_f": 0.025,
    }


@pytest.fixture(scope="module")
def fuelA_moderate_cylinder():
    """fuel-A-like XS, R = 10; τ_R = 5."""
    return {
        "R": 10.0, "sigma_t": 0.5, "sigma_s": 0.38, "nu_sigma_f": 0.025,
    }


# ═══════════════════════════════════════════════════════════════════════
# L1 — k_inf exactness gate (the load-bearing Phase-1 acceptance test)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-cylinder-architecture")
def test_v_alpha1_cyl_numerical_constant_initial_guess(
    fuelA_thin_cylinder,
):
    r"""L1 — V_α1_cyl numerical: constant ψ₀ lands on eigenmode in
    1 iteration.

    Validates the cylinder Variant α machinery (in-plane chord +
    axial-cosine + bounce-period closure) against V_α1_cyl algebraic
    identity :math:`(K \cdot 1) = \omega_0`. Uniform initial state →
    1 iteration → k_eff = k_inf to machine precision and φ uniform to
    machine precision.

    Cylinder analog of
    :func:`tests.derivations.test_peierls_greens_function_solver.test_v_alpha1_numerical_constant_initial_guess`.
    """
    fix = fuelA_thin_cylinder
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])

    res = solve_greens_function_cylinder(
        **fix, alpha=1.0,
        n_r=8, n_mu_axial=8, n_phi_az=16, n_traj_quad=24,
        max_iter=50, tol=1e-12,
    )

    assert res.converged, (
        f"V_α1_cyl numerical: did not converge in 50 iters; "
        f"k_eff = {res.k_eff:.8f}"
    )
    assert res.iterations <= 2, (
        f"V_α1_cyl numerical: should converge in ≤2 iter from const "
        f"initial guess; got {res.iterations}"
    )
    np.testing.assert_allclose(
        res.k_eff, k_inf, rtol=1e-10,
        err_msg=(
            f"V_α1_cyl numerical: k_eff = {res.k_eff} ≠ k_inf = {k_inf} "
            f"to 1e-10 tolerance"
        ),
    )
    rel_phi_spread = res.phi.std() / res.phi.mean()
    assert rel_phi_spread < 1e-10, (
        f"V_α1_cyl numerical: φ should be uniform; got "
        f"std/mean = {rel_phi_spread:.2e}"
    )


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-cylinder-architecture")
def test_v_alpha1_cyl_numerical_two_thicknesses(
    fuelA_thin_cylinder, fuelA_moderate_cylinder,
):
    r"""L1 — k_eff = k_inf at thin (τ_R = 2.5) and moderate (τ_R = 5)
    cylinders.

    Closed cylinder has no leakage regardless of thickness. The
    bounce-sum geometric :math:`T(b, \mu_{\rm axial}) = 1/(1 -
    e^{-\Sigma_t L_{\rm period}})` factor must produce the same k_inf
    at both optical depths.
    """
    for fix in (fuelA_thin_cylinder, fuelA_moderate_cylinder):
        k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])
        res = solve_greens_function_cylinder(
            **fix, alpha=1.0,
            n_r=8, n_mu_axial=8, n_phi_az=16, n_traj_quad=24,
            max_iter=20, tol=1e-12,
        )
        np.testing.assert_allclose(
            res.k_eff, k_inf, rtol=1e-10,
            err_msg=(
                f"τ_R = {fix['sigma_t']*fix['R']}: k_eff = "
                f"{res.k_eff} ≠ k_inf = {k_inf}"
            ),
        )


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-cylinder-architecture")
def test_v_alpha1_cyl_numerical_nonuniform_initial_guess(
    fuelA_thin_cylinder,
):
    r"""L1 — non-uniform ψ₀ converges to rank-1 isotropic eigenmode and
    recovers k_inf.

    Power iteration starting from an aggressive radial perturbation
    :math:`\psi_0(r, \mu_{\rm axial}, \varphi_{\rm az}) = 1 + 2(r/R)^2`
    must converge to the constant rank-1 eigenmode. The rank-1
    isotropic eigenmode is the unique closed-cylinder solution by
    V_α1_cyl, so any non-uniform start must collapse to it.
    """
    fix = fuelA_thin_cylinder
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])

    n_r, n_mu, n_phi = 12, 10, 24
    R = fix["R"]
    r_pts, _ = np.polynomial.legendre.leggauss(n_r)
    r_nodes = R * 0.5 * (r_pts + 1.0)

    # Strong radial perturbation, axisymmetric in (μ_axial, φ_az).
    psi0 = np.ones((n_r, n_mu, n_phi))
    for i, r in enumerate(r_nodes):
        psi0[i] *= (1.0 + 2.0 * (r / R) ** 2)

    res = solve_greens_function_cylinder(
        **fix, alpha=1.0, n_r=n_r, n_mu_axial=n_mu, n_phi_az=n_phi,
        n_traj_quad=32, max_iter=100, tol=1e-11, initial_psi=psi0,
    )

    assert res.converged, (
        f"V_α1_cyl nonuniform: did not converge in 100 iters; "
        f"k_eff = {res.k_eff}, iter = {res.iterations}"
    )
    rel_err = abs(res.k_eff - k_inf) / k_inf
    assert rel_err < 1e-9, (
        f"V_α1_cyl nonuniform: k_eff = {res.k_eff} ≠ k_inf = {k_inf} "
        f"by {rel_err:.3e}"
    )
    rel_phi_spread = res.phi.std() / res.phi.mean()
    assert rel_phi_spread < 1e-6, (
        f"V_α1_cyl nonuniform: φ should converge to uniform; got "
        f"std/mean = {rel_phi_spread:.2e}"
    )


# ═══════════════════════════════════════════════════════════════════════
# L1 — α=0 vacuum: physics-plausible flux profile + leakage
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-cylinder-architecture")
def test_alpha_zero_vacuum_produces_leaky_eigenmode(fuelA_thin_cylinder):
    r"""L1 — at α=0 vacuum BC, k_eff < k_inf and φ(r) is peaked at
    center.

    The vacuum case stress-tests the trajectory machinery — V_α1_cyl's
    algebraic cancellation does NOT apply when the surface fixed-point
    closure is zeroed (α=0 → ψ_surf = 0). The flux profile must be
    physically plausible: monotonically decreasing from center to
    surface (positive leakage at the boundary).
    """
    fix = fuelA_thin_cylinder
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])

    res = solve_greens_function_cylinder(
        **fix, alpha=0.0,
        n_r=16, n_mu_axial=12, n_phi_az=32, n_traj_quad=48,
        max_iter=300, tol=1e-9,
    )
    assert res.converged, f"vacuum: did not converge"
    # k_eff < k_inf (leakage)
    assert res.k_eff < k_inf, (
        f"vacuum: k_eff = {res.k_eff} should be < k_inf = {k_inf} "
        f"(leakage)"
    )
    # Sanity bound — fuel-A-like at R=5 has substantial leakage
    assert 0.4 < res.k_eff / k_inf < 0.7, (
        f"vacuum: k_eff/k_inf = {res.k_eff/k_inf} outside expected band "
        f"[0.4, 0.7]"
    )
    # Flux profile peaked at center (interior > surface)
    assert res.phi[0] > res.phi[-1], (
        f"vacuum: φ(r=0)={res.phi[0]:.4e} should exceed "
        f"φ(r=R)={res.phi[-1]:.4e} (peaked at center)"
    )
    # Monotone decrease (allow small numerical noise)
    diffs = np.diff(res.phi)
    assert np.all(diffs <= 1e-3), (
        f"vacuum: φ(r) not monotone-decreasing; diffs[max] = "
        f"{diffs.max():.3e}"
    )


# ═══════════════════════════════════════════════════════════════════════
# L1 — Multi-group at α=1: heterogeneous-XS homogeneous-medium
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-cylinder-architecture")
def test_mg_2g_asymmetric_scattering_at_alpha_one_gives_kinf():
    r"""L1 — 2G asymmetric scattering at α=1 must give k_inf exactly.

    Anti-pattern detector for the SigS-transpose convention drift
    (ERR-002 — sphere/cylinder MG solvers using ``SigS`` instead of
    ``SigS^T``). Asymmetric scattering matrix
    (:math:`\Sigma_{s,01} \ne \Sigma_{s,10}`) is the discriminator —
    symmetric matrices are blind to the convention error.

    Per :mod:`vv-principles` H1: 1-group eigenvalue tests are
    degenerate (k = νΣ_f/Σ_a is shape-independent). Multi-group with
    asymmetric XS is mandatory.
    """
    sigma_t = np.array([1.0, 0.6])
    sigma_s = np.array([
        [0.2, 0.5],   # g=0 → g=0: 0.2, g=0 → g=1: 0.5  (downscatter)
        [0.05, 0.3],  # g=1 → g=0: 0.05 (upscatter), g=1 → g=1: 0.3
    ])
    nu_sigma_f = np.array([0.05, 0.4])
    chi = np.array([1.0, 0.0])
    k_inf = kinf_homogeneous(sigma_t, sigma_s, nu_sigma_f, chi)

    res = solve_greens_function_cylinder_mg(
        R=5.0, sigma_t=sigma_t, sigma_s=sigma_s,
        nu_sigma_f=nu_sigma_f, chi=chi,
        alpha=1.0, n_r=12, n_mu_axial=10, n_phi_az=24,
        n_traj_quad=32, max_iter=100, tol=1e-11,
    )
    assert res.converged, (
        f"MG 2G α=1: did not converge in 100 iters"
    )
    rel_err = abs(res.k_eff - k_inf) / k_inf
    assert rel_err < 1e-9, (
        f"MG 2G α=1: k_eff = {res.k_eff} ≠ k_inf = {k_inf} "
        f"by {rel_err:.3e}"
    )
    # Each group's flux uniform (rank-1 isotropic eigenmode in r)
    for g in range(2):
        spread = res.phi_g[g].std() / res.phi_g[g].mean()
        assert spread < 1e-6, (
            f"MG 2G α=1, g={g}: φ_g spread {spread:.2e} not uniform"
        )


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-cylinder-architecture")
def test_mg_g1_reduces_to_1g_solver():
    r"""L1 — MG solver at G=1 reproduces the 1G solver bit-equal at
    α=0 vacuum.

    Architectural consistency: solving a G=1 problem via the MG path
    must give exactly the same answer as the 1G path (same XS, same
    grid, same iteration). Detects accidental MG-only special cases.
    """
    R, sigma_t, sigma_s, nu_sigma_f = 5.0, 0.5, 0.38, 0.025

    res_1g = solve_greens_function_cylinder(
        R=R, sigma_t=sigma_t, sigma_s=sigma_s, nu_sigma_f=nu_sigma_f,
        alpha=0.0,
        n_r=12, n_mu_axial=10, n_phi_az=24, n_traj_quad=32,
        max_iter=200, tol=1e-10,
    )
    res_mg = solve_greens_function_cylinder_mg(
        R=R, sigma_t=np.array([sigma_t]), sigma_s=np.array([[sigma_s]]),
        nu_sigma_f=np.array([nu_sigma_f]), chi=np.array([1.0]),
        alpha=0.0,
        n_r=12, n_mu_axial=10, n_phi_az=24, n_traj_quad=32,
        max_iter=200, tol=1e-10,
    )
    # Both should converge.
    assert res_1g.converged and res_mg.converged
    np.testing.assert_allclose(
        res_mg.k_eff, res_1g.k_eff, rtol=1e-9,
        err_msg=(
            f"MG G=1 reduction: MG k_eff = {res_mg.k_eff} != 1G k_eff = "
            f"{res_1g.k_eff}"
        ),
    )


# ═══════════════════════════════════════════════════════════════════════
# Convergence study (tabulated for posterity — not enforced as L1)
# ═══════════════════════════════════════════════════════════════════════


# ═══════════════════════════════════════════════════════════════════════
# V_α2_cyl numerical primitive cross-check — the structurally-independent
# L1 evidence for T_00^cyl = P_ss^cyl (no closed-form sphere-equivalent
# available because Ki_3 has no elementary closed form, so the SymPy
# proof can only get to the integrand level via the Bickley-Naylor
# identity bridge — see derive_T00_equals_P_ss_cylinder docstring).
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.parametrize(
    "tau_R",
    [0.5, 1.0, 2.5, 5.0, 10.0],
    ids=["tauR_0.5", "tauR_1.0", "tauR_2.5", "tauR_5.0", "tauR_10.0"],
)
@pytest.mark.verifies("peierls-greens-cylinder-T")
def test_v_alpha2_cyl_T00_equals_Pss_via_production_primitives(tau_R):
    r"""V_α2_cyl — **structurally-independent L1 cross-check**:
    :math:`T_{00}^{\rm cyl} = P_{ss}^{\rm cyl}` at the production-
    primitive level over a :math:`\tau_R` sweep.

    Why this test is load-bearing for cylinder V_α2: unlike sphere,
    cylinder has no elementary closed form (Ki_3 obstruction), so the
    SymPy V_α2_cyl proof
    (:func:`derive_T00_equals_P_ss_cylinder`) can only reach the
    integrand-identity level via the Bickley-Naylor bridge
    :math:`\mathrm{Ki}_3(x) = \int_0^{\pi/2}\sin^2\beta e^{-x/\sin\beta}\mathrm d\beta`
    — a known mathematical fact applied as a substitution, not a
    SymPy-derived equality. The rigorous V&V evidence for the rank-1
    Knyazev ≡ Hébert white-BC theorem must therefore come from the
    **numerical level**.

    This test pins the identity at the production-code level: the
    matrix function :func:`compute_T_specular_cylinder_3d` (Knyazev
    expansion + scattering kernel + matrix construction) and the
    scalar function :func:`compute_P_ss_cylinder` (slanted-chord polar
    integration → Ki_3 → in-plane :math:`\alpha` integration) are
    independent code paths — different functions, different
    intermediate quantities. At rank-1 (``n_modes = 1``,
    :math:`m = n = 0`, all Knyazev shifted-Legendre coefficients
    :math:`k_m = k_n = 0`) the matrix [0, 0] element MUST equal the
    scalar :math:`P_{ss}^{\rm cyl}`.

    The L1 evidence pillar for cylinder V_α2 is therefore this
    numerical cross-check, not the SymPy integrand identity (which
    must rely on the Bickley-Naylor bridge as a known identity).
    """
    R = 5.0
    sig_t = np.array([tau_R / R])
    radii = np.array([R])

    P_ss_scalar = compute_P_ss_cylinder(radii, sig_t, n_quad=64, dps=25)
    T_matrix = compute_T_specular_cylinder_3d(
        radii, sig_t, n_modes=1, n_quad=64,
    )
    T_00 = float(T_matrix[0, 0])

    np.testing.assert_allclose(
        T_00, P_ss_scalar, rtol=1e-10, atol=1e-12,
        err_msg=(
            f"V_α2_cyl production-primitive cross-check failed at "
            f"τ_R = {tau_R}: T_00 = {T_00:.16e}, P_ss = "
            f"{P_ss_scalar:.16e}, diff = {abs(T_00 - P_ss_scalar):.3e}"
        ),
    )


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("peierls-greens-cylinder-architecture")
def test_alpha_zero_convergence_floor():
    r"""L1 — α=0 vacuum k_eff convergence-floor (slow, **research-grade
    gate**).

    Pins the achieved accuracy floor for the α=0 vacuum case at
    successive quadrature orders, **asserting research-grade
    self-consistency between the two finest grids**. Phase 2
    unification must preserve this floor — any unified architecture
    that regresses the finest-pair self-consistency above
    ``RESEARCH_GRADE_FLOOR_CYL`` halts the unification per
    ``feedback_unify_after_two_instances.md``.

    Quadrature ladder (n_r, n_mu_axial, n_phi_az, n_traj_quad):

    - (12, 10, 24, 32) — coarse
    - (16, 12, 32, 48) — medium
    - (20, 16, 48, 64) — fine
    - (24, 20, 64, 96) — research-grade

    The Phase-1 closeout memo recorded :math:`8.24 \times 10^{-8}`
    self-consistency between the fine and research-grade grids. The
    floor assertion below is set at ``5e-7`` to allow for slight
    iteration-tolerance variation while still catching any genuine
    regression — anything above 5e-7 means the (24, 20, 64, 96)
    quadrature is no longer at research-grade convergence.

    fuel-A-like XS, τ_R = 2.5: k_eff ≈ 0.12045 (not k_inf; physics
    says :math:`0.5 < k_{\rm eff}/k_\infty < 0.7`).
    """
    fix = {
        "R": 5.0, "sigma_t": 0.5, "sigma_s": 0.38, "nu_sigma_f": 0.025,
    }
    orders = [
        (12, 10, 24, 32),
        (16, 12, 32, 48),
        (20, 16, 48, 64),
        (24, 20, 64, 96),  # research-grade — the Phase-2 benchmark
    ]
    k_values = []
    for n_r, n_mu, n_phi, n_t in orders:
        res = solve_greens_function_cylinder(
            **fix, alpha=0.0,
            n_r=n_r, n_mu_axial=n_mu, n_phi_az=n_phi, n_traj_quad=n_t,
            max_iter=400, tol=1e-10,
        )
        k_values.append(res.k_eff)

    # Coarse-to-fine sanity: every refinement step should at least
    # halve the error (roughly). Loose 1e-3 catastrophic-failure bound.
    for i in range(len(k_values) - 1):
        diff = abs(k_values[i + 1] - k_values[i]) / k_values[-1]
        assert diff < 1e-3, (
            f"Convergence-floor: order {orders[i]} → {orders[i+1]} "
            f"differs by {diff:.3e} relative — exceeds 1e-3 catastrophic "
            f"bound"
        )

    # ── Research-grade gate (the Phase-2 acceptance benchmark) ───────
    # Self-consistency between the two FINEST grids must be tight.
    # Closeout memo recorded ~8.24e-8 here; gate at 5e-7 to allow
    # iteration-tolerance fluctuation while catching regressions.
    RESEARCH_GRADE_FLOOR_CYL = 5e-7
    finest_pair_diff = (
        abs(k_values[-1] - k_values[-2]) / k_values[-1]
    )
    assert finest_pair_diff < RESEARCH_GRADE_FLOOR_CYL, (
        f"Research-grade floor: {orders[-2]} → {orders[-1]} "
        f"self-consistency = {finest_pair_diff:.3e}, exceeds "
        f"{RESEARCH_GRADE_FLOOR_CYL:.0e}. Phase-2 unification "
        f"benchmark broken — investigate before proceeding "
        f"(see feedback_unify_after_two_instances.md)."
    )

    # Pin the achieved value at the finest tested grid for posterity.
    assert 0.117 < k_values[-1] < 0.125, (
        f"Convergence-floor: finest k_eff = {k_values[-1]} outside "
        f"expected physics band [0.117, 0.125]"
    )
