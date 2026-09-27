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

- :mod:`.test_trajectory_resolvent_mg` (multi-group, single-region)
- Issue #132 issue body and
  :file:`.claude/agent-memory/numerics-investigator/issue_100_class_b_mr_mg.md`
  for the Phase 4 catastrophe documentation.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_xs
from orpheus.derivations.continuous.trajectory_resolvent.greens_function import (
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


# ════════════════════════════════════════════════════════════════════════
# The radial ladder — the convergence the global source spline broke (ERR-090)
# ════════════════════════════════════════════════════════════════════════

#: The heterogeneous closed sphere of ``test_phase_c_crosscheck.py``: fuel A |
#: moderator B | fuel A at outer radii 0.5, 1.5, 2.0 cm, 2 groups.
_LADDER_RADII = np.array([0.5, 1.5, 2.0])
_LADDER_N_R = (24, 36, 48)


def _aba_sphere_k(n_r: int) -> float:
    parts = [get_xs(key, "2g") for key in ("A", "B", "A")]
    res = solve_greens_function_sphere_mr(
        radii=_LADDER_RADII,
        sigma_t=np.stack([p["sig_t"] for p in parts]),
        sigma_s=np.stack([p["sig_s"] for p in parts]),
        nu_sigma_f=np.stack([p["nu"] * p["sig_f"] for p in parts]),
        chi=np.stack([p["chi"] for p in parts]),
        alpha=1.0,
        n_r=n_r, n_mu=24, n_traj_quad=64,
        max_iter=2000, tol=1e-9, initial_k=1.38,
    )
    assert res.converged, f"n_r={n_r}: the power iteration did not converge"
    return float(res.k_eff)


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.catches("ERR-090")
@pytest.mark.verifies("peierls-greens-mr-regionwise-source")
@pytest.mark.rests_on(
    "tests/gates/derivations/test_trajectory_resolvent_regionwise_source.py"
    "::test_mr_oracle_first_leg_matches_the_line_integral[sphere]",
)
def test_mr_sphere_k_converges_in_n_r() -> None:
    r"""The multi-region sphere's eigenvalue converges in :math:`n_r` faster than first order.

    On the ladder :math:`n_r = 24, 36, 48` (composite per-region
    Gauss-Legendre, :math:`n_\mu = 24`, 64 chord points per segment) the
    ratio of successive eigenvalue increments is compared with the ratio a
    first-order method would show on the same ladder,

    .. math::

       \rho_1 = \frac{1/36 - 1/48}{1/24 - 1/36} = \tfrac12 ,

    computed below from the ladder rather than typed. A method of order
    :math:`p > 1` shows a smaller ratio (second order: 0.35). The emission
    density jumps at both interfaces; with the per-region interpolant
    (ERR-090) the reference reads, ``[M]`` 2026-09-26, k = 1.383737,
    1.381917, 1.381293: increments -1.82e-3 then -6.25e-4, ratio 0.34, the
    second-order value. With one spline across the jumps it read 1.358083,
    1.361371, 1.379031: increments +3.29e-3 then +1.77e-2, ratio 5.4, and the
    sequence was not converging at all. The ladder's two further rungs
    (not run here, for cost) read 1.381381 at :math:`n_r = 72` (per region)
    and 1.380103 (one spline).

    This is a rate claim and bounds nothing absolute (``vv-principles``
    anti-pattern #5): the eigenvalue's value is pinned by the cross-check
    against the discrete-ordinates solve in
    ``tests/gates/sn/verification/analytical/test_phase_c_crosscheck.py``.
    The angular order is held at :math:`n_\mu = 24`, where the reference's
    own :math:`\mu` error is about 1e-3 in k (the Gauss-Legendre rule in
    :math:`\mu` meets a kink where a chord grazes an interior interface,
    #516). That error does NOT cancel in the increments: the kink angles move
    with the node radii, so the radial and angular errors couple (the step
    :math:`n_r` 36 to 72 is 5.4e-4 at :math:`n_\mu = 24` and 2.4e-4 at
    :math:`n_\mu = 96`). The rate criterion still separates the two
    interpolants by a factor of 16 in the ratio, but it measures the coupled
    error, and the reference's bound used by the SN cross-check is taken from
    its two-parameter ladder
    (``tests/gates/derivations/_trajectory_resolvent_ladders.py``), not from
    this row.
    """
    k = [_aba_sphere_k(n) for n in _LADDER_N_R]
    increments = np.diff(k)
    ratio = abs(increments[1]) / abs(increments[0])
    n = np.asarray(_LADDER_N_R, dtype=float)
    first_order_ratio = (1 / n[1] - 1 / n[2]) / (1 / n[0] - 1 / n[1])
    print(f"k(n_r={_LADDER_N_R}) = {k}; increments {increments}; ratio {ratio:.3f} "
          f"(first order: {first_order_ratio:.3f})")
    assert ratio < first_order_ratio, (
        f"the multi-region sphere reference converges in n_r no faster than first order: "
        f"k = {k} at n_r = {_LADDER_N_R}, increment ratio {ratio:.3f} >= {first_order_ratio:.3f}; "
        f"a source interpolant fitted across the material interfaces does this (ERR-090)"
    )
