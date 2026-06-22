r"""Phase 1b cylinder Variant α — multi-region verification gates.

Closes the ``cyl_2g_3reg`` ERR-026 gap identified at
``.claude/plans/issue_168_phase_c.md`` by pinning the cylinder MR
Variant α solver (:func:`solve_greens_function_cylinder_mr`) against
the six mandatory gates defined in
``.claude/plans/cylinder_mr_variant_alpha_verification.md``.

Gate coverage in this file:

- **Gate 1** (foundation): Homogeneous-limit reducibility.
  All-K-regions-same-material reproduces
  :func:`solve_greens_function_cylinder_mg` at FP-non-associativity
  floor. K=3 and K=5; 1G and 2G.
- **Gate 3** (L1): Single-region MR specular k_∞ recovery with
  asymmetric 2G SigS to expose Mode #6 convention drift.
- **Gate 4** (foundation): φ(r) continuity across material interfaces
  with 10× asymmetric contrast.
- **Gate 6** (L1): Per-axis quadrature refinement convergence on a
  stress configuration different from any other gate's setup
  (anti-tautology).
- **Gate 7** (foundation, optional): Reflective symmetry probe — the
  closed-cylinder MR operator preserves the V_α1_cyl_mr.q symbolic
  identity when material is uniform (algebraic-ancestor cross-check).

Gate 2 (L1 WM-72 vacuum cross-check) lives in the separate
``..._cylinder_mr_xverif.py`` file. Gate 5 is documented in the
module docstring there.

Branch-1 SymPy verifications backing these gates live at
:mod:`...origins.specular.greens_function_cylinder` (three
``derive_*_cylinder_mr*`` functions pinning the algebra). Their
foundation tests live in
``test_peierls_greens_function_cylinder_symbolic.py``.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.interpolate import CubicSpline

from orpheus.derivations.common.eigenvalue import kinf_homogeneous
from orpheus.derivations.common.xs_library import get_xs
from orpheus.derivations.continuous.trajectory_resolvent.greens_function_cylinder import (
    solve_greens_function_cylinder_mg,
    solve_greens_function_cylinder_mr,
)


# ════════════════════════════════════════════════════════════════════════
# Fixtures
# ════════════════════════════════════════════════════════════════════════


@pytest.fixture(scope="module")
def fuelA_1g_uniform_K3():
    """Gate 1, K=3, 1G — all 3 regions same fuel-A material."""
    xs = get_xs("A", "1g")
    R = 5.0
    radii = np.array([0.4 * R, 0.7 * R, R])
    n_regions = 3
    return {
        "R": R,
        "radii": radii,
        "sigma_t": np.tile(xs["sig_t"], (n_regions, 1)),
        "sigma_s": np.tile(xs["sig_s"][None, :, :], (n_regions, 1, 1)),
        "nu_sigma_f": np.tile((xs["nu"] * xs["sig_f"]), (n_regions, 1)),
        "chi": np.tile(xs["chi"], (n_regions, 1)),
        "scalar_xs": xs,
    }


@pytest.fixture(scope="module")
def fuelA_1g_uniform_K5():
    """Gate 1, K=5, 1G — all 5 regions same fuel-A material."""
    xs = get_xs("A", "1g")
    R = 5.0
    radii = np.array([0.2 * R, 0.4 * R, 0.6 * R, 0.8 * R, R])
    n_regions = 5
    return {
        "R": R,
        "radii": radii,
        "sigma_t": np.tile(xs["sig_t"], (n_regions, 1)),
        "sigma_s": np.tile(xs["sig_s"][None, :, :], (n_regions, 1, 1)),
        "nu_sigma_f": np.tile((xs["nu"] * xs["sig_f"]), (n_regions, 1)),
        "chi": np.tile(xs["chi"], (n_regions, 1)),
        "scalar_xs": xs,
    }


@pytest.fixture(scope="module")
def fuelA_2g_uniform_K3():
    """Gate 1, K=3, 2G — all regions same fuel-A 2G material."""
    xs = get_xs("A", "2g")
    R = 5.0
    radii = np.array([0.4 * R, 0.7 * R, R])
    n_regions = 3
    return {
        "R": R,
        "radii": radii,
        "sigma_t": np.tile(xs["sig_t"], (n_regions, 1)),
        "sigma_s": np.tile(xs["sig_s"][None, :, :], (n_regions, 1, 1)),
        "nu_sigma_f": np.tile((xs["nu"] * xs["sig_f"]), (n_regions, 1)),
        "chi": np.tile(xs["chi"], (n_regions, 1)),
        "scalar_xs": xs,
    }


@pytest.fixture(scope="module")
def asymmetric_2g_xs():
    """Gate 3 — asymmetric 2G SigS designed to expose Mode #6 convention
    drift (per ``cylinder_mr_variant_alpha_verification.md`` §3 Gate 3).
    """
    return {
        "sigma_t": np.array([1.0, 1.5]),
        # SigS[g_from, g_to]; strongly asymmetric (downscatter 0.3 vs
        # upscatter 0.1) — guards against transpose bugs.
        "sigma_s": np.array([[0.4, 0.3], [0.1, 0.8]]),
        "nu_sigma_f": np.array([0.6, 0.9]),
        "chi": np.array([1.0, 0.0]),
    }


@pytest.fixture(scope="module")
def interface_continuity_xs():
    """Gate 4 — asymmetric 3-region 1G XS profile (10× contrast)."""
    R = 4.0
    # Plan §3 Gate 4 reads "radii = [0.5·R, 1.0·R, R]"; the second entry
    # equals R, which collapses K=3 into K=2 (the second interface IS
    # the outer surface). Interpreting as a typo for K=3 with a
    # genuine interior interface at 0.75·R:
    #   interface_inner @ 0.5·R = 2.0
    #   interface_mid   @ 0.75·R = 3.0
    #   outer surface   @ R = 4.0
    radii = np.array([0.5 * R, 0.75 * R, R])
    # 10× σ_t contrast at the inner interface, asymmetric profile.
    sigma_t = np.array([2.0, 0.5, 1.5])
    sigma_s = np.array([1.8, 0.3, 1.0])
    nu_sigma_f = np.array([0.1, 0.0, 0.05])
    return {
        "R": R,
        "radii": radii,
        "sigma_t": sigma_t[:, None],            # (3, 1)
        "sigma_s": sigma_s.reshape(3, 1, 1),    # (3, 1, 1)
        "nu_sigma_f": nu_sigma_f[:, None],      # (3, 1)
        "interface_radii": radii[:-1],          # interior boundaries
    }


# ════════════════════════════════════════════════════════════════════════
# Gate 1 — Homogeneous-limit reducibility (foundation)
# ════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies(
    "peierls-greens-cylinder-mr-homogeneous-reduction",
    "peierls-greens-cylinder-mr-trajectory-segments",
)
def test_mr_K3_uniform_reduces_to_mg_1g(fuelA_1g_uniform_K3):
    r"""Gate 1.K3.1G — K=3 uniform-material MR reproduces MG bit-tight.

    Cardinal Bit-Identity criterion 1 (``vv-principles`` § "Bit-identity
    vs principled-equivalence"): with all 3 regions assigned identical
    fuel-A material, the MR solver MUST reproduce
    :func:`solve_greens_function_cylinder_mg` to within
    :math:`K \cdot n_{\rm traj\_quad} \cdot {\rm ULP}` —
    the FP-non-associativity floor for the piecewise integration tree.

    This is the load-bearing inheritance gate. If it passes, ALL
    homogeneous-cylinder verification (Sood ``Ua-1-O-CY`` at 8.5e-6,
    V_α1_cyl k_∞ at machine precision, V_α2_cyl T_00 = P_ss) inherits
    into the MR code path. Every subsequent gate then tests only the
    piecewise-specific machinery.
    """
    fix = fuelA_1g_uniform_K3
    xs = fix["scalar_xs"]

    res_mg = solve_greens_function_cylinder_mg(
        R=fix["R"],
        sigma_t=xs["sig_t"],
        sigma_s=xs["sig_s"],
        nu_sigma_f=xs["nu"] * xs["sig_f"],
        chi=xs["chi"],
        alpha=1.0,
        n_r=12, n_mu_axial=12, n_phi_az=24, n_traj_quad=32,
        max_iter=20, tol=1e-12,
    )
    res_mr = solve_greens_function_cylinder_mr(
        radii=fix["radii"],
        sigma_t=fix["sigma_t"],
        sigma_s=fix["sigma_s"],
        nu_sigma_f=fix["nu_sigma_f"],
        chi=fix["chi"],
        alpha=1.0,
        n_r=12, n_mu_axial=12, n_phi_az=24, n_traj_quad=32,
        max_iter=20, tol=1e-12,
    )
    np.testing.assert_allclose(
        res_mr.k_eff, res_mg.k_eff, rtol=1e-10,
        err_msg=(
            f"Gate 1 K=3 1G uniform: MR k_eff = {res_mr.k_eff:.16e} "
            f"differs from MG k_eff = {res_mg.k_eff:.16e} by more than "
            f"the K·n_traj·ULP FP-non-associativity floor (~3e-12 for "
            "this configuration)."
        ),
    )


@pytest.mark.foundation
@pytest.mark.verifies("peierls-greens-cylinder-mr-homogeneous-reduction")
def test_mr_K5_uniform_reduces_to_mg_1g(fuelA_1g_uniform_K5):
    """Gate 1.K5.1G — K=5 uniform-material MR reproduces MG.

    Exercises a longer accumulation chain than K=3. A bug in segment
    ordering would show K=3 passing while K=5 fails (or vice versa).
    """
    fix = fuelA_1g_uniform_K5
    xs = fix["scalar_xs"]

    res_mg = solve_greens_function_cylinder_mg(
        R=fix["R"],
        sigma_t=xs["sig_t"], sigma_s=xs["sig_s"],
        nu_sigma_f=xs["nu"] * xs["sig_f"], chi=xs["chi"],
        alpha=1.0,
        n_r=12, n_mu_axial=12, n_phi_az=24, n_traj_quad=32,
        max_iter=20, tol=1e-12,
    )
    res_mr = solve_greens_function_cylinder_mr(
        radii=fix["radii"], sigma_t=fix["sigma_t"],
        sigma_s=fix["sigma_s"], nu_sigma_f=fix["nu_sigma_f"],
        chi=fix["chi"], alpha=1.0,
        n_r=12, n_mu_axial=12, n_phi_az=24, n_traj_quad=32,
        max_iter=20, tol=1e-12,
    )
    np.testing.assert_allclose(
        res_mr.k_eff, res_mg.k_eff, rtol=1e-10,
        err_msg=(
            f"Gate 1 K=5 1G uniform: MR k_eff = {res_mr.k_eff:.16e} "
            f"differs from MG k_eff = {res_mg.k_eff:.16e}"
        ),
    )


@pytest.mark.foundation
@pytest.mark.verifies("peierls-greens-cylinder-mr-homogeneous-reduction")
def test_mr_K3_uniform_reduces_to_mg_2g(fuelA_2g_uniform_K3):
    """Gate 1.K3.2G — K=3 uniform 2G MR reproduces MG.

    The 2G sub-gate is **mandatory** per the verification plan §3 Gate
    1 anti-pattern flag: a 1G-only test would be blind to the SigS-
    transpose-on-MR-broadcast variant of Mode #6 (convention drift) —
    the symmetric ``SigS[k, g, g]`` 1G case cannot distinguish
    ``SigS`` from ``SigS^T``.
    """
    fix = fuelA_2g_uniform_K3
    xs = fix["scalar_xs"]

    res_mg = solve_greens_function_cylinder_mg(
        R=fix["R"],
        sigma_t=xs["sig_t"], sigma_s=xs["sig_s"],
        nu_sigma_f=xs["nu"] * xs["sig_f"], chi=xs["chi"],
        alpha=1.0,
        n_r=10, n_mu_axial=8, n_phi_az=16, n_traj_quad=24,
        max_iter=50, tol=1e-10,
    )
    res_mr = solve_greens_function_cylinder_mr(
        radii=fix["radii"], sigma_t=fix["sigma_t"],
        sigma_s=fix["sigma_s"], nu_sigma_f=fix["nu_sigma_f"],
        chi=fix["chi"], alpha=1.0,
        n_r=10, n_mu_axial=8, n_phi_az=16, n_traj_quad=24,
        max_iter=50, tol=1e-10,
    )
    np.testing.assert_allclose(
        res_mr.k_eff, res_mg.k_eff, rtol=1e-9,
        err_msg=(
            f"Gate 1 K=3 2G uniform: MR k_eff = {res_mr.k_eff:.16e} "
            f"differs from MG k_eff = {res_mg.k_eff:.16e}; "
            f"a divergence here likely signals an SigS-transpose "
            f"convention drift in the MR broadcast (Mode #6)."
        ),
    )


# ════════════════════════════════════════════════════════════════════════
# Gate 3 — Single-region MR specular-BC k_∞ recovery (L1 closed-form)
# ════════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-cylinder-mr-kinf")
def test_mr_single_region_kinf_2g_asymmetric_sigs(asymmetric_2g_xs):
    r"""Gate 3 — single-region MR specular BC recovers k_∞ exactly.

    For ``radii=[R]`` with α=1, the V_α1_cyl algebraic identity
    (Branch-1 :func:`derive_operator_constant_trial_closed_cylinder`)
    ensures :math:`k_{\rm eff} = k_\infty = \max{\rm eig}({\bf A}^{-1}{\bf F})`
    exactly — the closed-cylinder bounce-sum collapses to flat trial
    flux regardless of geometry partitioning (V_α1_cyl_mr.q
    homogeneous-limit reducibility).

    The 2G asymmetric ``SigS`` deliberately picks ``SigS[0,1] = 0.3``
    versus ``SigS[1,0] = 0.1`` so that ``SigS ≠ SigS^T``. This guards
    against the L0-SN-009 (ERR-002) Mode #6 convention-drift pattern:
    a symmetric ``SigS`` would let any ``SigS``/``SigS^T`` broadcast
    swap pass undetected.

    Target: ``|k_eff - k_∞| / k_∞ ≤ 1e-6`` (the actual algebraic
    identity is exact; the floor is power-iteration tolerance +
    cubic-spline interpolation residual for a constant trial — both
    near machine precision at adequate quadrature).
    """
    xs = asymmetric_2g_xs

    k_inf = kinf_homogeneous(
        xs["sigma_t"], xs["sigma_s"], xs["nu_sigma_f"], xs["chi"],
    )

    res = solve_greens_function_cylinder_mr(
        radii=np.array([5.0]),
        sigma_t=xs["sigma_t"][None, :],
        sigma_s=xs["sigma_s"][None, :, :],
        nu_sigma_f=xs["nu_sigma_f"][None, :],
        chi=xs["chi"][None, :],
        alpha=1.0,
        n_r=16, n_mu_axial=12, n_phi_az=24, n_traj_quad=32,
        max_iter=80, tol=1e-12,
    )
    assert res.converged, (
        f"Gate 3 single-region 2G asymmetric: power iteration did not "
        f"converge; k_eff = {res.k_eff}, iters = {res.iterations}"
    )
    rel_err = abs(res.k_eff - k_inf) / abs(k_inf)
    # 1e-6 is the tight target; plan §3 Gate 3 brief sets a ceiling of
    # 5e-4. We expect machine-precision agreement modulo power-iter
    # tolerance.
    assert rel_err < 1e-6, (
        f"Gate 3 single-region 2G asymmetric: k_eff = {res.k_eff:.16e} "
        f"vs k_inf = {k_inf:.16e}, rel_err = {rel_err:.3e} > 1e-6. "
        f"Suspect SigS convention drift (Mode #6) or chi/nu_sigma_f "
        f"broadcast bug."
    )


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-cylinder-mr-kinf")
def test_mr_single_region_kinf_1g_fuelA(fuelA_2g_uniform_K3):
    """Gate 3 sanity — 1G fuel-A single-region MR α=1 recovers k_∞.

    Companion check at 1G: a simpler config to verify the V_α1_cyl
    invariance survives at K=1 for any standard XS. This is NOT a
    substitute for the 2G asymmetric test (1G is degenerate per
    the 1-group-degeneracy rule); it's an additional safety net.
    """
    xs = get_xs("A", "1g")
    # 1G XS from get_xs("A", "1g"): sig_t (1,), sig_s (1, 1), nu*sig_f (1,)
    k_inf_anal = float(
        (xs["nu"] * xs["sig_f"])[0]
        / (xs["sig_t"][0] - xs["sig_s"][0, 0])
    )

    res = solve_greens_function_cylinder_mr(
        radii=np.array([5.0]),
        sigma_t=xs["sig_t"][None, :],
        sigma_s=xs["sig_s"][None, :, :],
        nu_sigma_f=(xs["nu"] * xs["sig_f"])[None, :],
        chi=xs["chi"][None, :],
        alpha=1.0,
        n_r=12, n_mu_axial=12, n_phi_az=24, n_traj_quad=32,
        max_iter=30, tol=1e-12,
    )
    np.testing.assert_allclose(
        res.k_eff, k_inf_anal, rtol=1e-10,
        err_msg=(
            f"Gate 3 sanity 1G: k_eff = {res.k_eff} vs k_inf = "
            f"{k_inf_anal}"
        ),
    )


# ════════════════════════════════════════════════════════════════════════
# Gate 4 — Interface φ(r) continuity (foundation)
# ════════════════════════════════════════════════════════════════════════


@pytest.mark.slow
@pytest.mark.foundation
@pytest.mark.verifies("peierls-greens-cylinder-mr-interface-continuity")
def test_mr_interface_continuity_3region(interface_continuity_xs):
    r"""Gate 4 — φ(r) continuous across material interfaces.

    For a piecewise-homogeneous medium, the angular flux ψ has a
    discontinuous derivative at material interfaces (the source profile
    jumps with σ_t), but the scalar flux φ = ∫ψ dμ dφ_az is
    **continuous** because the angular integration averages out the
    radial-direction discontinuity.

    Test: 3-region 1G cylinder with 10× σ_t contrast at the inner
    interface, asymmetric profile [2.0, 0.5, 1.5]. Build per-region
    splines of φ vs r restricted to nodes within each region. Evaluate
    the spline limit at the exact interface radius from each side. The
    relative jump must be bounded by a multiple of the angular-
    quadrature floor (the test uses 1e-2 abs — a permissive ceiling
    given the prototype's single-domain GL grid that interpolates
    across the source jump; the SHAPE result of the test is more
    important than the exact bound, but the bound rules out gross
    region-indexing bugs).

    Anti-pattern protection per the plan §3 Gate 4: asymmetric profile
    (10× contrast) prevents sign-flip / symmetric-error bugs from
    averaging out across symmetrically-placed interfaces.

    Post-Phase-E (2026-05-12, Issue #168 Phase D Step 4b follow-up):
    the radial quadrature is now composite per-region GL (the
    canonical treatment for σ_t interface kinks), so the GL nodes
    stay strictly inside each region.  This shifts the spline error
    profile from "spline-across-jump" to "spline-extrapolation-to-
    interface".  Each region's spline must extrapolate from its
    interior GL nodes to the material boundary; the extrapolation
    error from each side is now the dominant term in ``rel_jump``.
    The ceiling was relaxed from 1e-2 → 5e-2 to accommodate this
    trade-off; the gain is monotonic k_eff convergence under
    refinement (see Gate 6 quadrature-convergence tests) and
    accurate flux profiles inside each region.
    """
    fix = interface_continuity_xs

    res = solve_greens_function_cylinder_mr(
        radii=fix["radii"],
        sigma_t=fix["sigma_t"],
        sigma_s=fix["sigma_s"],
        nu_sigma_f=fix["nu_sigma_f"],
        alpha=1.0,
        n_r=36, n_mu_axial=16, n_phi_az=32, n_traj_quad=48,
        max_iter=300, tol=1e-9,
    )
    assert res.converged, (
        f"Gate 4: did not converge; k_eff = {res.k_eff}, "
        f"iters = {res.iterations}"
    )

    phi = res.phi_g[0]
    r_nodes = res.r_nodes
    interfaces = fix["interface_radii"]  # interior boundaries

    for r_iface in interfaces:
        # Locate left/right region indices
        lreg = int(np.searchsorted(fix["radii"], r_iface, side="left"))
        rreg = lreg + 1
        l_mask = res.region_at_node == lreg
        r_mask = res.region_at_node == rreg
        assert l_mask.sum() >= 4 and r_mask.sum() >= 4, (
            f"Gate 4: not enough nodes to build per-region splines at "
            f"interface r = {r_iface}; left = {l_mask.sum()}, "
            f"right = {r_mask.sum()}. Increase n_r."
        )
        spl_l = CubicSpline(r_nodes[l_mask], phi[l_mask], extrapolate=True)
        spl_r = CubicSpline(r_nodes[r_mask], phi[r_mask], extrapolate=True)
        phi_l = float(spl_l(r_iface))
        phi_r = float(spl_r(r_iface))
        rel_jump = abs(phi_l - phi_r) / max(abs(phi_l), abs(phi_r))
        # 5e-2 ceiling, post-Phase-E composite-GL correction (2026-05-12):
        # with single-domain GL on (0, R) the nodes straddled the
        # material interfaces, so per-region cubic splines could
        # interpolate-across-jump and the rel_jump was ~3e-3.  With
        # composite per-region GL (the canonical correction for σ_t
        # interface kinks — see ``_composite_per_region_gl`` in
        # ``orpheus.derivations.continuous.trajectory_resolvent.greens_function``),
        # the GL nodes stay STRICTLY INSIDE each region; the spline
        # must EXTRAPOLATE from interior nodes to the interface, and
        # the extrapolation error from each side is now the dominant
        # term in ``rel_jump`` (not the source-jump smoothing artifact
        # the original 1e-2 ceiling was calibrated for).  Empirical
        # post-Phase-E value: ~2.85e-2 on the asymmetric 3-region
        # config; the 5e-2 ceiling rules out off-by-one region_at_node
        # bugs (≥ 10% jumps) while accommodating the extrapolation
        # trade-off.  Quadrature-convergence stress tests
        # (``test_mr_quadrature_convergence_*``) pin the MONOTONIC
        # k_eff convergence the composite-GL correction delivers.
        assert rel_jump < 5e-2, (
            f"Gate 4 interface r = {r_iface}: φ left limit = {phi_l:.6e}, "
            f"φ right limit = {phi_r:.6e}, relative jump = {rel_jump:.3e}. "
            f"Exceeds 5e-2 ceiling — investigate region_at_node or "
            f"segment-arithmetic bug."
        )


# ════════════════════════════════════════════════════════════════════════
# Gate 6 — Per-axis quadrature refinement convergence (L1)
# ════════════════════════════════════════════════════════════════════════


# Gate 6 stress config: a 2-region 1G cylinder different from any other
# gate's config (anti-tautology). Picked to have moderate-to-strong
# multiplication so power iteration converges in O(50) iters.
#
# Plan §3 Gate 6 specifies a 3-region 2G config; we use a 1G analogue
# with strong σ_t contrast because the cylinder MR oracle is Python-
# triple-nested (no JIT) — full-spec 3R/2G at n_phi_az=64 / n_mu=32
# would push the test runtime above CI budget. The 1G config still
# exercises every piece of MR machinery (segmentation, piecewise τ,
# per-region source).
#
# Stress configuration:
#   R = 4.0; radii = [0.4·R, R] = [1.6, 4.0] (K=2)
#   inner region: fuel-A (σ_t=1, σ_s=0.5, νσ_f=0.75)
#   outer region: artificial strong-scatterer (σ_t=1.0, σ_s=0.95, νσ_f=0)
#   α = 0.7 — partial reflection, neither limit
GATE6_R = 4.0
GATE6_RADII = np.array([0.4 * GATE6_R, GATE6_R])
GATE6_SIGMA_T = np.array([1.0, 1.0])[:, None]
GATE6_SIGMA_S = np.array([[[0.5]], [[0.95]]])
GATE6_NU_SIGMA_F = np.array([0.75, 0.0])[:, None]


def _run_gate6(n_r, n_mu, n_phi, n_chord, max_iter=200, tol=1e-9):
    return solve_greens_function_cylinder_mr(
        radii=GATE6_RADII,
        sigma_t=GATE6_SIGMA_T,
        sigma_s=GATE6_SIGMA_S,
        nu_sigma_f=GATE6_NU_SIGMA_F,
        alpha=0.7,
        n_r=n_r, n_mu_axial=n_mu, n_phi_az=n_phi, n_traj_quad=n_chord,
        max_iter=max_iter, tol=tol,
    )


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("peierls-greens-cylinder-mr-quadrature-convergence")
def test_mr_quadrature_convergence_chord():
    r"""Gate 6 — chord-quadrature refinement on K=2 1G α=0.7 stress.

    Refine ``n_traj_quad`` ∈ {12, 24, 48} at fixed
    ``(n_r, n_mu_axial, n_phi_az) = (12, 12, 24)``. Assert successive
    |Δk_eff| are decreasing (improving rates) and the ratio of
    successive deltas is ≤ 0.5 (super-algebraic — Gauss-Legendre on
    a smooth integrand gives geometric convergence).

    Note: plan §3 Gate 6 specifies refinement sequences up to
    n_traj_quad=128 at n_phi_az=64. The cylinder MR oracle's Python
    triple loop makes that intractable for a CI test budget; we use a
    coarser sequence at the same refinement-rate-check principle.
    Convergence rate is a **necessary** condition for L1 evidence
    (``numerical-bug-signatures`` H4) — failure of monotonicity here
    signals non-smooth integrand contributions (typically a per-
    segment piece dropped or boundary handling).
    """
    res_coarse = _run_gate6(12, 12, 24, 12)
    res_medium = _run_gate6(12, 12, 24, 24)
    res_fine = _run_gate6(12, 12, 24, 48)

    dk_1 = abs(res_medium.k_eff - res_coarse.k_eff)
    dk_2 = abs(res_fine.k_eff - res_medium.k_eff)

    assert dk_2 < dk_1, (
        f"Gate 6 chord: |Δk_eff| not monotonically decreasing under "
        f"chord refinement; dk(12→24)={dk_1:.3e}, dk(24→48)={dk_2:.3e}. "
        f"Indicates non-smooth integrand bug."
    )
    # 0.5 ratio bound (geometric convergence expected for smooth GL).
    assert dk_2 / max(dk_1, 1e-30) < 0.5, (
        f"Gate 6 chord: convergence rate too shallow; "
        f"dk_2/dk_1 = {dk_2/dk_1:.3f} ≥ 0.5"
    )


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("peierls-greens-cylinder-mr-quadrature-convergence")
def test_mr_quadrature_convergence_phi_az():
    """Gate 6 — azimuthal-quadrature refinement."""
    res_coarse = _run_gate6(12, 12, 12, 32)
    res_medium = _run_gate6(12, 12, 24, 32)
    res_fine = _run_gate6(12, 12, 48, 32)

    dk_1 = abs(res_medium.k_eff - res_coarse.k_eff)
    dk_2 = abs(res_fine.k_eff - res_medium.k_eff)

    assert dk_2 < dk_1, (
        f"Gate 6 phi_az: dk monotonicity violation: "
        f"dk(12→24)={dk_1:.3e}, dk(24→48)={dk_2:.3e}"
    )


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("peierls-greens-cylinder-mr-quadrature-convergence")
def test_mr_quadrature_convergence_mu_axial():
    """Gate 6 — axial-cosine quadrature refinement."""
    res_coarse = _run_gate6(12, 6, 24, 32)
    res_medium = _run_gate6(12, 12, 24, 32)
    res_fine = _run_gate6(12, 18, 24, 32)

    dk_1 = abs(res_medium.k_eff - res_coarse.k_eff)
    dk_2 = abs(res_fine.k_eff - res_medium.k_eff)

    assert dk_2 < dk_1, (
        f"Gate 6 mu_axial: dk monotonicity violation: "
        f"dk(6→12)={dk_1:.3e}, dk(12→18)={dk_2:.3e}"
    )


# ════════════════════════════════════════════════════════════════════════
# Gate 7 (optional, foundation) — Algebraic ancestor cross-check
# ════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies("peierls-greens-cylinder-mr-homogeneous-reduction")
def test_mr_branch_1_branch_2_algebraic_ancestor():
    r"""L1 cross-check — Branch 1 (SymPy V_α1_cyl_mr.q) and Branch 2
    (numpy ``solve_greens_function_cylinder_mr``) agree at the
    homogeneous limit.

    Branch 1 source: :func:`derive_two_region_constant_source_consistency_cylinder_mr`
    proves that ``psi_surf = q/Σ_t`` exactly when Σ_t,1 = Σ_t,2 = Σ_t
    and α = 1 — the closed-cylinder constant-trial fixed point lifts
    verbatim into the 2-region MR machinery.

    Branch 2 manifestation: a 2-region cylinder MR with both regions
    set to identical fuel-A material and α=1 must produce
    ``k_eff = k_∞ = ν Σ_f / Σ_a`` exactly (the operator on a constant
    flat trial is multiplicatively closed — V_α1_cyl_mr.q says
    F + e^{-τ_first} ψ_surf = q/Σ_t = ψ_trial).

    Above the trusted-library line, Branch 1 (SymPy symbolic algebra)
    and Branch 2 (numpy + scipy CubicSpline + numpy.linalg.eigvals)
    are STRUCTURALLY INDEPENDENT: no shared in-house primitives. They
    share only ``np.exp``, ``scipy.interpolate.CubicSpline`` and
    ``numpy.polynomial.legendre.leggauss`` which are trusted-library
    primitives (``algebra-of-record`` § "Structural independence
    applies above the trusted-library line").

    Acceptance: rel_err ≤ 1e-10 (machine-precision modulo power-
    iteration tolerance; the V_α1_cyl_mr.q identity says equality is
    algebraically exact).
    """
    xs = get_xs("A", "1g")
    R = 5.0
    radii = np.array([0.6 * R, R])  # K=2, different from K=3 K=5 fixtures
    sigma_t = np.tile(xs["sig_t"], (2, 1))
    sigma_s = np.tile(xs["sig_s"][None, :, :], (2, 1, 1))
    nuf = np.tile((xs["nu"] * xs["sig_f"]), (2, 1))
    chi = np.tile(xs["chi"], (2, 1))

    # k_inf analytical (Branch-1 algebraic ancestor result).
    # 1G XS arrays from get_xs("A", "1g"): sig_t shape (1,), sig_s
    # shape (1, 1), nu*sig_f shape (1,).
    k_inf = float(
        (xs["nu"] * xs["sig_f"])[0]
        / (xs["sig_t"][0] - xs["sig_s"][0, 0])
    )

    # Branch-2 numerical with K=2 uniform material
    res = solve_greens_function_cylinder_mr(
        radii=radii, sigma_t=sigma_t, sigma_s=sigma_s,
        nu_sigma_f=nuf, chi=chi, alpha=1.0,
        n_r=12, n_mu_axial=12, n_phi_az=24, n_traj_quad=32,
        max_iter=30, tol=1e-13,
    )
    rel_err = abs(res.k_eff - k_inf) / k_inf
    assert rel_err < 1e-10, (
        f"L1 Branch-1/Branch-2 cross-check: V_α1_cyl_mr.q says K=2 "
        f"uniform → k_eff = k_inf exactly; got k_eff = {res.k_eff:.16e}, "
        f"k_inf = {k_inf:.16e}, rel_err = {rel_err:.3e}"
    )
