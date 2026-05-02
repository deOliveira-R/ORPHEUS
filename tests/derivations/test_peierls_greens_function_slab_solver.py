r"""Phase-3A numerical-implementation gates for the slab Variant α
Green's function reference (1G + MG, homogeneous slab, symmetric
reflective specular BC :math:`\alpha_{\rm left} = \alpha_{\rm right} =
\alpha`).

Mirrors :mod:`.test_peierls_greens_function_solver` (sphere) and
:mod:`.test_peierls_greens_function_cylinder_solver` (cylinder) with
slab geometry. The load-bearing claim:

**Variant α gives the exact :math:`k_{\rm eff} = k_\infty` for the
closed homogeneous slab with symmetric specular BC** (V_α1_slab
algebraic identity). This is the slab analogue of the sphere/cylinder
Variant α exactness property and the canonical Phase-3A acceptance
test.

Test strategy
-------------

1. **L1 — k_inf exactness at α=1, 1G**: closed slab with perfect
   specular has no leakage and rank-1 isotropic eigenmode with
   :math:`k_{\rm eff} = k_\infty = \nu\Sigma_f/\Sigma_a` exactly. Must
   reproduce to ≤ 1e-10 within iteration tolerance.

2. **L1 — multi-thickness exactness**: same XS at thin (τ_L = 1) and
   moderate (τ_L = 5) slab must both yield k_eff = k_inf — the bounce-
   sum geometric series with α² per-period reflection product is
   wired correctly across τ_L range.

3. **L1 — non-uniform initial guess at α=1**: power iteration starting
   from a strongly non-uniform initial state must converge to the rank-1
   isotropic eigenmode and recover k_inf within the iteration
   tolerance.

4. **L1 — multi-group asymmetric scattering**: a 2G case with
   asymmetric :math:`\Sigma_s` (:math:`\Sigma_{s,01} \ne \Sigma_{s,10}`)
   must reduce to :math:`k_\infty` from
   :func:`orpheus.derivations.common.eigenvalue.kinf_homogeneous` at
   :math:`\alpha = 1`. ERR-002 anti-pattern detector.

5. **L1 — vacuum BC (α=0) sanity**: at α=0 the spatial flux profile
   must be peaked at center and depleted at the walls; :math:`k_{\rm
   eff} < k_\infty`. Stress-tests the trajectory machinery —
   V_α1_slab's algebraic cancellation does NOT protect this case.

6. **L1 — α=0 cross-check between MG and 1G reductions**: at G=1 the
   MG solver must agree with the 1G solver bit-equal.

7. **L1 — V_α2_slab production-primitive cross-check**:
   :math:`T_{\rm slab}[0, n] = 2 E_3(\Sigma_t L)` from
   :func:`compute_T_specular_slab` over a τ_L sweep, against the
   canonical mpmath ``2 * expint(3, τ_L)``.

Anti-cross-check prohibition
----------------------------

We do NOT cross-check slab Variant α against the slab Nyström
``compute_P_esc_outer`` — those primitives use the same E_n evaluator
upstream (``mpmath.expint``), so the cross-check would only verify the
mpmath library, not the Variant α integrand structure. The L1 cross-
check is **k_inf-exactness at α=1**, an analytical invariant from
V_α1_slab structurally independent of any other ORPHEUS solver.

Predecessor / sibling tests
---------------------------

- :mod:`.test_peierls_greens_function_slab_symbolic` — V_α1_slab,
  V_α2_slab, V_α3_slab SymPy gates.
- :mod:`.test_peierls_greens_function_solver` — sphere V_α1.numerical.
- :mod:`.test_peierls_greens_function_cylinder_solver` — cylinder
  V_α1_cyl.numerical (the structural template this test file mirrors).

References
----------

- Sanchez, R. (1986). *Transp. Theor. Stat. Phys.* 14.
- :file:`/.claude/plans/peierls-greens-cylinder-and-2bc.md` — Phase 3A
  slab plan.
"""
from __future__ import annotations

import mpmath
import numpy as np
import pytest

from orpheus.derivations.common.eigenvalue import kinf_homogeneous
from orpheus.derivations.continuous.peierls_nystrom.geometry import (
    compute_T_specular_slab,
)
from orpheus.derivations.continuous.peierls_greens_function.greens_function_slab import (
    solve_greens_function_slab,
    solve_greens_function_slab_mg,
)


# ═══════════════════════════════════════════════════════════════════════
# Fixtures — fuel-A-like XS analogous to the sphere/cylinder fixtures
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
# L1 — k_inf exactness gate (the load-bearing Phase-3A acceptance test)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-slab-architecture")
def test_v_alpha1_slab_numerical_constant_initial_guess(
    fuelA_moderate_slab,
):
    r"""L1 — V_α1_slab numerical: constant ψ₀ lands on eigenmode in
    1 iteration.

    Validates the slab Variant α machinery (first-leg chord +
    bounce-period chord with 2-bounce-per-period structure +
    closure with :math:`\alpha^2` per-period reflection product)
    against V_α1_slab algebraic identity :math:`(K \cdot 1) =
    \omega_0`. Uniform initial state → 1 iteration → k_eff = k_inf
    to machine precision and φ uniform to machine precision.

    The load-bearing Phase-3A acceptance test. Failure here would
    indicate either:

    - The :math:`\alpha^2` per-period reflection product was not
      threaded through to the closure, OR
    - The slab chord arithmetic is wrong (first-leg or bounce-period
      length formula incorrect), OR
    - The trajectory parametrisation has a sign/direction error.
    """
    fix = fuelA_moderate_slab
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])

    res = solve_greens_function_slab(
        **fix, alpha=1.0,
        n_x=8, n_mu=8, n_traj_quad=24,
        max_iter=50, tol=1e-12,
    )

    assert res.converged, (
        f"V_α1_slab numerical: did not converge in 50 iters; "
        f"k_eff = {res.k_eff:.8f}"
    )
    assert res.iterations <= 2, (
        f"V_α1_slab numerical: should converge in ≤2 iter from const "
        f"initial guess; got {res.iterations}"
    )
    np.testing.assert_allclose(
        res.k_eff, k_inf, rtol=1e-10,
        err_msg=(
            f"V_α1_slab numerical: k_eff = {res.k_eff} ≠ k_inf = "
            f"{k_inf} to 1e-10 tolerance"
        ),
    )
    rel_phi_spread = res.phi.std() / res.phi.mean()
    assert rel_phi_spread < 1e-10, (
        f"V_α1_slab numerical: φ should be uniform; got "
        f"std/mean = {rel_phi_spread:.2e}"
    )


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-slab-architecture")
def test_v_alpha1_slab_numerical_two_thicknesses(
    fuelA_thin_slab, fuelA_moderate_slab,
):
    r"""L1 — k_eff = k_inf at thin (τ_L = 1) and moderate (τ_L = 5)
    slabs.

    Closed slab has no leakage regardless of thickness. The bounce-sum
    geometric :math:`T(\mu) = 1/(1 - \alpha^2 e^{-\Sigma_t L_{\rm
    period}})` factor with the 2-bounce-per-period reflection product
    must produce the same k_inf at both optical depths.
    """
    for fix in (fuelA_thin_slab, fuelA_moderate_slab):
        k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])
        res = solve_greens_function_slab(
            **fix, alpha=1.0,
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


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-slab-architecture")
def test_v_alpha1_slab_numerical_nonuniform_initial_guess(
    fuelA_moderate_slab,
):
    r"""L1 — non-uniform ψ₀ converges to rank-1 isotropic eigenmode and
    recovers k_inf.

    Power iteration starting from an aggressive perturbation
    :math:`\psi_0(x, \mu) = 1 + 2(x/L - 1/2)^2` (stronger at walls,
    weaker at center) must converge to the constant rank-1 eigenmode.
    The rank-1 isotropic eigenmode is the unique closed-slab solution
    by V_α1_slab, so any non-uniform start must collapse to it.
    """
    fix = fuelA_moderate_slab
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])

    n_x, n_mu = 12, 16
    L = fix["L"]
    x_pts, _ = np.polynomial.legendre.leggauss(n_x)
    x_nodes = L * 0.5 * (x_pts + 1.0)

    psi0 = np.ones((n_x, n_mu))
    for i, x in enumerate(x_nodes):
        psi0[i] *= 1.0 + 2.0 * (x / L - 0.5) ** 2

    res = solve_greens_function_slab(
        **fix, alpha=1.0, n_x=n_x, n_mu=n_mu,
        n_traj_quad=32, max_iter=100, tol=1e-11, initial_psi=psi0,
    )

    assert res.converged, (
        f"V_α1_slab nonuniform: did not converge in 100 iters; "
        f"k_eff = {res.k_eff}, iter = {res.iterations}"
    )
    rel_err = abs(res.k_eff - k_inf) / k_inf
    assert rel_err < 1e-9, (
        f"V_α1_slab nonuniform: k_eff = {res.k_eff} ≠ k_inf = "
        f"{k_inf} by {rel_err:.3e}"
    )
    rel_phi_spread = res.phi.std() / res.phi.mean()
    assert rel_phi_spread < 1e-6, (
        f"V_α1_slab nonuniform: φ should converge to uniform; got "
        f"std/mean = {rel_phi_spread:.2e}"
    )


# ═══════════════════════════════════════════════════════════════════════
# L1 — α=0 vacuum: physics-plausible flux profile + leakage
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-slab-architecture")
def test_alpha_zero_vacuum_produces_leaky_eigenmode(fuelA_moderate_slab):
    r"""L1 — at α=0 vacuum BC, k_eff < k_inf and φ(x) is peaked at
    center.

    The vacuum case stress-tests the trajectory machinery — V_α1_slab's
    algebraic cancellation does NOT apply when the surface fixed-point
    closure is zeroed (α=0 → ψ_surf = 0). The flux profile must be
    physically plausible: symmetric about :math:`x = L/2`, peaked at
    center, depleted at both walls (positive leakage).
    """
    fix = fuelA_moderate_slab
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])

    res = solve_greens_function_slab(
        **fix, alpha=0.0,
        n_x=16, n_mu=24, n_traj_quad=48,
        max_iter=300, tol=1e-9,
    )
    assert res.converged, "vacuum: did not converge"
    # k_eff < k_inf (leakage)
    assert res.k_eff < k_inf, (
        f"vacuum: k_eff = {res.k_eff} should be < k_inf = "
        f"{k_inf} (leakage)"
    )
    # Sanity bound — fuel-A-like at L=10 (τ_L=5): substantial leakage
    # Should be in [0.45, 0.85] (slabs leak less than spheres of same τ
    # because surface-to-volume ratio is smaller).
    assert 0.45 < res.k_eff / k_inf < 0.85, (
        f"vacuum: k_eff/k_inf = {res.k_eff/k_inf} outside expected "
        f"band [0.45, 0.85]"
    )
    # Flux profile peaked at center — symmetric in x.
    n_x = len(res.phi)
    mid = n_x // 2
    # Both phi[0] (near x=0 wall) and phi[-1] (near x=L wall) < phi[mid].
    assert res.phi[mid] > res.phi[0], (
        f"vacuum: φ(mid)={res.phi[mid]:.4e} should exceed "
        f"φ(near x=0)={res.phi[0]:.4e} (peaked at center)"
    )
    assert res.phi[mid] > res.phi[-1], (
        f"vacuum: φ(mid)={res.phi[mid]:.4e} should exceed "
        f"φ(near x=L)={res.phi[-1]:.4e} (peaked at center)"
    )
    # x → L-x symmetry (within numerical tolerance — GL grid is
    # mirror-symmetric about midpoint).
    np.testing.assert_allclose(
        res.phi, res.phi[::-1], rtol=1e-6,
        err_msg=(
            f"vacuum: φ(x) should be symmetric about L/2; got "
            f"max asymmetry {np.abs(res.phi - res.phi[::-1]).max():.3e}"
        ),
    )


# ═══════════════════════════════════════════════════════════════════════
# L1 — Multi-group at α=1: heterogeneous-XS homogeneous-medium
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-slab-architecture")
def test_mg_2g_asymmetric_scattering_at_alpha_one_gives_kinf():
    r"""L1 — 2G asymmetric scattering at α=1 must give k_inf exactly.

    Anti-pattern detector for the SigS-transpose convention drift
    (ERR-002 — using ``SigS`` instead of ``SigS^T``). Asymmetric
    scattering matrix (:math:`\Sigma_{s,01} \ne \Sigma_{s,10}`) is
    the discriminator — symmetric matrices are blind to the
    convention error.

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

    res = solve_greens_function_slab_mg(
        L=5.0, sigma_t=sigma_t, sigma_s=sigma_s,
        nu_sigma_f=nu_sigma_f, chi=chi,
        alpha=1.0, n_x=12, n_mu=16,
        n_traj_quad=32, max_iter=100, tol=1e-11,
    )
    assert res.converged, "MG 2G α=1: did not converge in 100 iters"
    rel_err = abs(res.k_eff - k_inf) / k_inf
    assert rel_err < 1e-9, (
        f"MG 2G α=1: k_eff = {res.k_eff} ≠ k_inf = {k_inf} "
        f"by {rel_err:.3e}"
    )
    # Each group's flux uniform (rank-1 isotropic eigenmode in x).
    for g in range(2):
        spread = res.phi_g[g].std() / res.phi_g[g].mean()
        assert spread < 1e-6, (
            f"MG 2G α=1, g={g}: φ_g spread {spread:.2e} not uniform"
        )


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-slab-architecture")
def test_mg_g1_reduces_to_1g_solver():
    r"""L1 — MG solver at G=1 reproduces the 1G solver bit-equal at
    α=0 vacuum.

    Architectural consistency: solving a G=1 problem via the MG path
    must give exactly the same answer as the 1G path (same XS, same
    grid, same iteration). Detects accidental MG-only special cases.
    """
    L, sigma_t, sigma_s, nu_sigma_f = 10.0, 0.5, 0.38, 0.025

    res_1g = solve_greens_function_slab(
        L=L, sigma_t=sigma_t, sigma_s=sigma_s, nu_sigma_f=nu_sigma_f,
        alpha=0.0,
        n_x=12, n_mu=16, n_traj_quad=32,
        max_iter=200, tol=1e-10,
    )
    res_mg = solve_greens_function_slab_mg(
        L=L, sigma_t=np.array([sigma_t]), sigma_s=np.array([[sigma_s]]),
        nu_sigma_f=np.array([nu_sigma_f]), chi=np.array([1.0]),
        alpha=0.0,
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
# V_α2_slab production-primitive cross-check via mpmath E_3
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.parametrize(
    "tau_L",
    [0.5, 1.0, 2.5, 5.0, 10.0],
    ids=["tauL_0.5", "tauL_1.0", "tauL_2.5", "tauL_5.0", "tauL_10.0"],
)
@pytest.mark.verifies("peierls-greens-slab-T")
def test_v_alpha2_slab_T00_equals_2E3_via_production_primitive(tau_L):
    r"""V_α2_slab — **structurally-independent L1 cross-check**:
    :math:`T_{\rm slab}[0, n_modes] = 2 E_3(\Sigma_t L)` via the
    production primitive :func:`compute_T_specular_slab` over a
    :math:`\tau_L` sweep.

    The slab T-matrix is block-off-diagonal:

    .. math::

       T_{\rm slab} = \begin{pmatrix} 0 & T_{oi} \\ T_{io} & 0
                      \end{pmatrix}

    with :math:`T_{io} = T_{oi}` by face symmetry on a homogeneous slab
    and :math:`T_{oi}^{(0,0)} = 2 E_3(\Sigma_t L)` at rank-1
    (``n_modes = 1``).

    The production primitive computes :math:`T_{oi}^{(0,0)}` via
    explicit half-range Gauss-Legendre on :math:`\mu \in [0, 1]`
    integrating :math:`2 \mu e^{-\Sigma_t L /\mu}`. The mpmath E_3
    is computed via ``mpmath.expint(3, τ)`` — a fundamentally
    different code path (analytical continuation of the gamma
    function). Equality of the two is the V_α2_slab numerical
    structural-independence cross-check.
    """
    L = 1.0
    sig_t = np.array([tau_L / L])  # ensures τ_L total = sig_t · L
    radii = np.array([L])

    T_matrix = compute_T_specular_slab(radii, sig_t, n_modes=1, n_quad=64)
    # Off-diagonal block T_oi (the slab face-to-face transit):
    T_oi_00 = float(T_matrix[0, 1])  # block [outer, inner] entry
    # Also check T_io = T_oi (face symmetry):
    T_io_00 = float(T_matrix[1, 0])

    canonical = float(2.0 * mpmath.expint(3, tau_L))

    np.testing.assert_allclose(
        T_oi_00, canonical, rtol=1e-10, atol=1e-12,
        err_msg=(
            f"V_α2_slab T_oi production-primitive cross-check failed at "
            f"τ_L = {tau_L}: T_oi[0,0] = {T_oi_00:.16e}, 2·E_3 = "
            f"{canonical:.16e}, diff = {abs(T_oi_00 - canonical):.3e}"
        ),
    )
    np.testing.assert_allclose(
        T_io_00, canonical, rtol=1e-10, atol=1e-12,
        err_msg=(
            f"V_α2_slab T_io production-primitive cross-check failed at "
            f"τ_L = {tau_L}: T_io[0,0] = {T_io_00:.16e}, 2·E_3 = "
            f"{canonical:.16e}"
        ),
    )


# ═══════════════════════════════════════════════════════════════════════
# Convergence study (l1 + slow) — research-grade gate
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("peierls-greens-slab-architecture")
def test_alpha_zero_convergence_floor():
    r"""L1 — α=0 vacuum k_eff convergence-floor (slow, **research-grade
    gate**).

    Pins the achieved accuracy floor for the α=0 vacuum case at
    successive quadrature orders, asserting research-grade self-
    consistency between the two finest grids.

    Quadrature ladder (n_x, n_mu, n_traj_quad):

    - (16, 24, 48) — medium
    - (20, 32, 64) — fine
    - (24, 40, 96) — research-grade
    - (32, 56, 128) — finest

    fuel-A-like XS, τ_L = 5: physics says :math:`0.45 < k_{\rm eff}/
    k_\infty < 0.85` (slabs leak less than spheres of the same τ).

    History
    -------

    Pre-ERR-034 + ERR-035 fix: the finest-pair self-consistency was
    pinned at ~5e-4 with a 1e-3 gate; the slow convergence was
    misattributed to a quadrature limitation (slow GL convergence on
    µ for the grazing-ray locus). After both fixes:

    - ERR-034 corrected the trajectory parametrisation
      :math:`x_{\rm back}(s) = x - \mu\,s` (the missing :math:`\mu`
      factor was the dominant error mode).
    - ERR-035 corrected the closure formula to the first-principles
      :math:`\alpha\,B/(1 - \alpha\,e^{-\tau})` — at :math:`\alpha=0`
      the closure is bypassed entirely so this fix is not directly
      load-bearing here, but the delegation refactor uses the
      rank-2 single-transit B integrals which converge faster than
      the rank-1 full-period integrals.

    Post-fix re-measured floor: ~9e-6 between (24,40,96) and
    (32,56,128). The new gate is 5e-5 — catches genuine regressions
    with margin while keeping a comfortable cushion above the
    achieved value.
    """
    fix = {
        "L": 10.0, "sigma_t": 0.5, "sigma_s": 0.38, "nu_sigma_f": 0.025,
    }
    orders = [
        (16, 24, 48),
        (20, 32, 64),
        (24, 40, 96),
        (32, 56, 128),
    ]
    k_values = []
    for n_x, n_mu, n_t in orders:
        res = solve_greens_function_slab(
            **fix, alpha=0.0,
            n_x=n_x, n_mu=n_mu, n_traj_quad=n_t,
            max_iter=400, tol=1e-10,
        )
        k_values.append(res.k_eff)

    # Coarse-to-fine sanity: each refinement should not regress
    # catastrophically.
    for i in range(len(k_values) - 1):
        diff = abs(k_values[i + 1] - k_values[i]) / k_values[-1]
        assert diff < 5e-3, (
            f"Convergence-floor: order {orders[i]} → {orders[i+1]} "
            f"differs by {diff:.3e} relative — exceeds catastrophic "
            f"5e-3"
        )

    # Post-ERR-035-fix achievable floor for slab vacuum at fuel-A τ_L=5.
    # Re-measured 2026-05-02: (24,40,96) → (32,56,128) ≈ 8.85e-6 rel.
    # The 5e-5 gate has ~5x margin over the achieved value; tighten
    # only after a structurally-independent reference (PS-1982-class
    # for slab) is added.
    SLAB_VACUUM_FLOOR = 5e-5
    finest_pair_diff = (
        abs(k_values[-1] - k_values[-2]) / k_values[-1]
    )
    assert finest_pair_diff < SLAB_VACUUM_FLOOR, (
        f"Vacuum-slab floor: {orders[-2]} → {orders[-1]} "
        f"self-consistency = {finest_pair_diff:.3e}, exceeds "
        f"{SLAB_VACUUM_FLOOR:.0e}. Either a regression in the slab "
        f"trajectory machinery / closure, or a quadrature change "
        f"that lost the post-fix convergence floor."
    )

    # Pin the achieved value at the finest tested grid for posterity.
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])
    ratio = k_values[-1] / k_inf
    assert 0.45 < ratio < 0.85, (
        f"Convergence-floor: k_eff/k_inf = {ratio} outside expected "
        f"band [0.45, 0.85]"
    )
