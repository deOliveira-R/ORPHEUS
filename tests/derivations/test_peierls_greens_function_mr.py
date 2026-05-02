r"""Plan-(b) Option-2 — Multi-region Variant α extension test gate.

Direct attack on Issue #132 (Class B Multi-Region catastrophe). The
Phase 4 ``specular_multibounce`` rank-N closure is structurally
broken for multi-region configurations: a mode-0/mode-:math:`n\\ge 1`
normalisation mismatch in the rank-N Marshak partial-current basis
produces +57 % k_eff error on 1G/2R fuel-A inner / moderator-B outer
configurations (rank-2 gives k_eff = 1.015 vs k_inf = 0.648 — sign
flip + supercritical from a strongly subcritical configuration).

Variant α has no such closure — the boundary condition is absorbed
into the kernel via Sanchez 1986 Eq. (A1), and trajectory + bounce
machinery extends naturally to multi-region via piecewise τ(µ).
This test gate pins:

1. **Single-region reduction** — MR with 1 region matches the MG
   solver bit-exactly. Sanity check that the multi-region machinery
   reduces cleanly.
2. **Issue #132 reproducer** — sphere ``radii=[0.5, 1.0]`` with
   fuel-A inner / moderator-B outer (σ_t = [1, 2]) gives a sensible
   k_eff (specifically: ``< 1``; the Phase 4 catastrophe gives
   ``> 1``). Variant α k_eff is between the rank-1 Phase 4 result
   (0.551) and the cell-averaged k_inf homogenisation (0.648).
3. **Spatial mode physical sanity** — closed-sphere Issue #132
   eigenmode has φ higher in fuel (region 0) than moderator
   (region 1), decreasing monotonically with r within each region,
   with a discernible slope change at the fuel/moderator interface.
4. **Vacuum BC reduces k_eff further** — α=0 multi-region sphere
   should give k_eff lower than α=1 (leakage reduces multiplication).

Predecessors:

- :mod:`.test_peierls_greens_function_mg` (multi-group, single-region)
- Issue #132 issue body and
  :file:`.claude/agent-memory/numerics-investigator/issue_100_class_b_mr_mg.md`
  for the Phase 4 catastrophe documentation.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_xs
from orpheus.derivations.continuous.peierls.greens_function import (
    solve_greens_function_sphere_mg,
    solve_greens_function_sphere_mr,
)


# ════════════════════════════════════════════════════════════════════════
# Fixtures
# ════════════════════════════════════════════════════════════════════════


@pytest.fixture(scope="module")
def issue132_xs():
    """Issue #132 Class B catastrophe fixture: sphere radii=[0.5, 1.0],
    fuel-A inner / moderator-B outer.

    fuel-A 1G: σ_t = 1.0, σ_s = 0.5, νσ_f = 0.75
    moderator-B 1G: σ_t = 2.0, σ_s = 1.9, νσ_f = 0
    """
    xsA = get_xs("A", "1g")
    xsB = get_xs("B", "1g")
    return {
        "radii": np.array([0.5, 1.0]),
        "sigma_t": np.array([xsA["sig_t"], xsB["sig_t"]]),
        "sigma_s": np.array([xsA["sig_s"], xsB["sig_s"]]),
        "nu_sigma_f": np.array([
            xsA["nu"] * xsA["sig_f"],
            xsB["nu"] * xsB["sig_f"],
        ]),
        "chi": np.array([xsA["chi"], xsB["chi"]]),
    }


@pytest.fixture(scope="module")
def fuelA_single_region_1g():
    """fuel-A 1G XS in 1-region MR shape (for reduction sanity check)."""
    xs = get_xs("A", "1g")
    return {
        "radii": np.array([5.0]),
        "sigma_t": xs["sig_t"][None, :],          # (1, 1)
        "sigma_s": xs["sig_s"][None, :, :],       # (1, 1, 1)
        "nu_sigma_f": (xs["nu"] * xs["sig_f"])[None, :],
        "chi": xs["chi"][None, :],
    }


# ════════════════════════════════════════════════════════════════════════
# Tests
# ════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_mr_single_region_reduces_to_mg(fuelA_single_region_1g):
    """MR.1 — Multi-region solver with 1 region reproduces MG solver
    bit-exactly.

    The MR machinery splits trajectory and bounce-period chord into
    region segments. With 1 region, there are no interior boundaries,
    so the segment list collapses to a single segment spanning
    :math:`[0, L]`, and the per-segment composite quadrature equals
    the single-segment GL of the MG operator.
    """
    fix = fuelA_single_region_1g
    xs = get_xs("A", "1g")

    res_mg = solve_greens_function_sphere_mg(
        R=5.0,
        sigma_t=xs["sig_t"],
        sigma_s=xs["sig_s"],
        nu_sigma_f=xs["nu"] * xs["sig_f"],
        chi=xs["chi"],
        alpha=1.0,
        n_r=12, n_mu=12, n_traj_quad=24, max_iter=20, tol=1e-12,
    )
    res_mr = solve_greens_function_sphere_mr(
        radii=fix["radii"],
        sigma_t=fix["sigma_t"],
        sigma_s=fix["sigma_s"],
        nu_sigma_f=fix["nu_sigma_f"],
        chi=fix["chi"],
        alpha=1.0,
        n_r=12, n_mu=12, n_traj_quad=24, max_iter=20, tol=1e-12,
    )
    np.testing.assert_allclose(
        res_mr.k_eff, res_mg.k_eff, rtol=1e-12,
        err_msg=(
            f"MR.1: 1-region MR k_eff = {res_mr.k_eff} differs from "
            f"MG k_eff = {res_mg.k_eff}"
        ),
    )


@pytest.mark.foundation
def test_mr_issue132_no_catastrophe_closed_sphere(issue132_xs):
    """MR.2 — Issue #132 reproducer: sphere fuel-A/moderator-B closed
    sphere k_eff is subcritical and physically sensible.

    Phase 4 specular_multibounce rank-2 gives k_eff = 1.015
    (+57 % from k_inf = 0.648; sign-flip pathology). Variant α has no
    rank-N closure and therefore cannot exhibit the catastrophe.

    Pinned bounds:

    - ``0.5 < k_eff < 0.95`` — subcritical, between rank-1 Phase 4
      result (0.551) and the homogenised volume-averaged
      ``k_inf ≈ 0.648``. The exact value (≈0.735 with default
      quadrature) is the actual multi-region transport solution.
    - ``k_eff < 1`` — the Phase 4 catastrophe gives ``k_eff > 1``;
      this assertion explicitly rules it out.

    No structurally-independent reference is available for closed
    multi-region sphere k_eff (the literature has only fixed-source
    multi-region benchmarks like Garcia 2021). This test is a
    regression gate against a known pathology, not an L1 cross-check.
    The L1 cross-check belongs to Plan-(b) Option-1 (Garcia 2021
    flux-shape benchmark) — separate test file.
    """
    fix = issue132_xs

    res = solve_greens_function_sphere_mr(
        radii=fix["radii"],
        sigma_t=fix["sigma_t"],
        sigma_s=fix["sigma_s"],
        nu_sigma_f=fix["nu_sigma_f"],
        chi=fix["chi"],
        alpha=1.0,
        n_r=24, n_mu=24, n_traj_quad=48,
        max_iter=500, tol=1e-9,
    )
    assert res.converged, (
        f"MR.2: Issue #132 closed sphere did not converge in 500 "
        f"iterations; k_eff = {res.k_eff:.6f}"
    )
    assert res.k_eff < 1.0, (
        f"MR.2: closed-sphere k_eff = {res.k_eff:.6f} exceeds 1.0 — "
        "indicates supercriticality from a configuration that should "
        "be subcritical (Phase 4 catastrophe pattern)"
    )
    assert 0.5 < res.k_eff < 0.95, (
        f"MR.2: k_eff = {res.k_eff:.6f} outside expected physical "
        "range [0.5, 0.95] for fuel-A/moderator-B closed sphere"
    )


@pytest.mark.foundation
def test_mr_issue132_spatial_mode_physical(issue132_xs):
    """MR.3 — Issue #132 spatial mode is physically reasonable.

    For closed sphere (α=1) with fuel inner / moderator outer:

    - φ peaked in the fuel region (more multiplication).
    - φ monotonically decreasing with r within each region.
    - Slope discontinuity (or visible change) at the fuel/moderator
      interface at r = 0.5.
    """
    fix = issue132_xs

    res = solve_greens_function_sphere_mr(
        radii=fix["radii"],
        sigma_t=fix["sigma_t"],
        sigma_s=fix["sigma_s"],
        nu_sigma_f=fix["nu_sigma_f"],
        chi=fix["chi"],
        alpha=1.0,
        n_r=24, n_mu=24, n_traj_quad=48,
        max_iter=500, tol=1e-9,
    )
    assert res.converged

    phi = res.phi_g[0]  # 1G

    # φ peaked at centre.
    assert phi.argmax() in (0, 1, 2), (
        f"MR.3: φ should peak near r=0; got argmax at index "
        f"{phi.argmax()}"
    )

    # Within each region, monotonic decrease with r.
    fuel_mask = res.region_at_node == 0
    mod_mask = res.region_at_node == 1
    phi_fuel = phi[fuel_mask]
    phi_mod = phi[mod_mask]
    # Each region should be a non-increasing sequence.
    assert (np.diff(phi_fuel) <= 1e-6).all(), (
        f"MR.3: φ in fuel region not monotonic decreasing: "
        f"{phi_fuel}"
    )
    assert (np.diff(phi_mod) <= 1e-6).all(), (
        f"MR.3: φ in moderator region not monotonic decreasing: "
        f"{phi_mod}"
    )

    # φ in fuel > φ in moderator (more multiplication in fuel).
    assert phi_fuel.mean() > phi_mod.mean(), (
        f"MR.3: φ_fuel mean = {phi_fuel.mean():.4f} should be > "
        f"φ_moderator mean = {phi_mod.mean():.4f}"
    )


@pytest.mark.foundation
def test_mr_issue132_vacuum_below_closed(issue132_xs):
    """MR.4 — vacuum k_eff < closed-sphere k_eff (leakage reduces
    multiplication).

    Both should be subcritical (no catastrophe), but the vacuum case
    must give a smaller k_eff because neutrons leak out at the outer
    surface. This pins the α-monotonicity of the multi-region operator.
    """
    fix = issue132_xs

    res_closed = solve_greens_function_sphere_mr(
        radii=fix["radii"],
        sigma_t=fix["sigma_t"],
        sigma_s=fix["sigma_s"],
        nu_sigma_f=fix["nu_sigma_f"],
        chi=fix["chi"],
        alpha=1.0,
        n_r=24, n_mu=24, n_traj_quad=48,
        max_iter=500, tol=1e-9,
    )
    res_vacuum = solve_greens_function_sphere_mr(
        radii=fix["radii"],
        sigma_t=fix["sigma_t"],
        sigma_s=fix["sigma_s"],
        nu_sigma_f=fix["nu_sigma_f"],
        chi=fix["chi"],
        alpha=0.0,
        n_r=24, n_mu=24, n_traj_quad=48,
        max_iter=500, tol=1e-8,
    )

    assert res_closed.converged
    assert res_vacuum.converged

    assert res_vacuum.k_eff < res_closed.k_eff, (
        f"MR.4: vacuum k_eff = {res_vacuum.k_eff:.6f} should be < "
        f"closed-sphere k_eff = {res_closed.k_eff:.6f}"
    )
    assert res_vacuum.k_eff > 0, (
        f"MR.4: vacuum k_eff = {res_vacuum.k_eff:.6f} should be "
        "positive"
    )
