r"""D5b — multi-D (2-D Cartesian) Linear-Discontinuous on the DAG wavefront.

Sub-step **D5b-S2** of #240 / #38 / #37: the bilinear (UBLD) LD cell kernel now
runs in d≥2 on the wavefront (``FullFieldWavefront`` / ``MovingFrontierWindow``),
consuming the verified d-generic primitive
(:mod:`orpheus.transport.spatial._ubld`).  This module carries the END-TO-END 2-D LD
gates that need a full solve:

* **D5b.2 (smoke)** — the 2-D LD kernel produces a CONVERGENT O(h²) solution on
  the existing heterogeneous 2-D MMS (vacuum edges → no ``AngularBoundaryFlux`` change,
  the S2-permitted optional smoke-check).  The full strengthened stress-ansatz
  MMS (the vv Mode-7 override — ``SN2DCartesianLDStressMMSCase`` + the
  ``@verifies("ld-cartesian-2d")`` claim) is **S4** (it needs the non-vanishing
  boundary trace); this smoke-check verifies the LD kernel CONVERGES, not the
  full flux-shape claim.
* **D5b.3 (foundation)** — ``FullFieldWavefront`` ≡ ``MovingFrontierWindow``:
  the SAME 2-D LD via two DAG schedules.  Het, 2G-asymmetric, μ-non-trivial,
  NON-SQUARE — each leg forced to its INTENDED rep.
* **D5b.5 (foundation)** — DD ≠ LD: the same mesh/materials/source converges to
  DIFFERENT fluxes under ``DiamondDifference`` vs ``LinearDiscontinuous`` — the
  catcher that LD is GENUINELY LD, not a silent DD collapse (the D5-0 misroute).

D5b.0 (the INVERTED raise pin), D5b.1 (the kernel round-trip), and the
D5b.4 kernel-level matvec twin (``residual_kernel_batch`` faces reconstruct
identically to ``cell_kernel_batch``, both directions per axis) live with the
kernel in ``tests/transport/spatial/test_linear_discontinuous.py``
(``TestLDKernel``).  The FULL ``loss_action`` / Krylov 2-D LD matvec needs the
``2^d``-moment spatial iterate (S3); at S2 the production driver is the
within-group source-iteration sweep.

``-O``-safe (Mode 8): the load-bearing checks are ``np.testing.assert_*`` /
``pytest.fail`` (function calls that fire under ``python -O``), never bare
``assert``.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.derivations.continuous.mms.sn import (
    build_2d_cartesian_heterogeneous_mms_case,
    build_2d_cartesian_ld_stress_mms_case,
    build_nonvacuum_fixed_source,
)
from orpheus.geometry import BC, CoordSystem, Mesh2D
from orpheus.numerics.moment_layout import AVERAGE_MOMENT
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import solve_sn_fixed_source
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.loss_representation import (
    FullFieldWavefront,
    MovingFrontierWindow,
    default_for,
)
from orpheus.transport.spatial import DiamondDifference, LinearDiscontinuous
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux

from tests.sn._test_helpers import volume_weighted_l2


def _nonsquare_het_2g_mesh(nx: int = 5, ny: int = 4):
    r"""A NON-SQUARE, vacuum, 2-group het SNMesh on a level-symmetric quad.

    NON-SQUARE (``nx ≠ ny``) is the x↔y-swap defence; ``level_symmetric``
    supplies genuine ``mu_y`` (#214-safe); vacuum edges keep the domain inflow
    zero (no S4 boundary-trace widening needed)."""
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 1.3, nx + 1),
        edges_y=np.linspace(0.0, 0.9, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("vacuum"), bc_ymax=BC("vacuum"),
    )
    return SNMesh(
        mesh, Quadrature.level_symmetric(4), {0: get_mixture("A", "2g")},
        scheme=LinearDiscontinuous(),
    )


# ═══════════════════════════════════════════════════════════════════════
# D5b.2 (smoke) — 2-D LD converges O(h²) on the het MMS (vacuum edges)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.slow
def test_ld_2d_converges_second_order_smoke():
    r"""D5b.2 (smoke): the 2-D LD kernel produces a CONVERGENT O(h²) solution.

    Runs LD on the existing heterogeneous 2-D MMS (its ansatz VANISHES on all
    four edges → vacuum BCs automatic → no ``AngularBoundaryFlux`` change, the
    S2-permitted optional smoke-check).  Asserts the observed convergence order
    climbs toward 2 (last > 1.85, all > 1.7) AND the converged value lands in a
    sane band against ``phi_exact`` (vv §5: rate ≠ correctness).

    NOT the full ``@verifies("ld-cartesian-2d")`` claim: that needs the vv
    Mode-7 stress ansatz (``SN2DCartesianLDStressMMSCase``, μ-non-trivial slope
    drivers) which has a NON-vanishing boundary trace — S4.  The isotropic
    sin·sin ansatz under-stresses LD's bilinear slope coupling, so this is a
    CONVERGENCE smoke-check (the LD kernel runs end-to-end and is second-order),
    not a flux-shape verification.  Hence ``@l1`` but no ``verifies`` label."""
    case = build_2d_cartesian_heterogeneous_mms_case()
    n_cells = [20, 40, 80]
    errors = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        materials = case.build_materials(mesh)
        Q = case.external_source(mesh)
        result = solve_sn_fixed_source(
            materials, mesh, case.quadrature, Q,
            boundary_condition="vacuum", max_inner=500, inner_tol=1e-12,
            scheme=LinearDiscontinuous(),
        )
        phi = result.scalar_flux.values                     # (ng, nx, ny)
        sq = 0.0
        for g in range(phi.shape[0]):
            ref = case.phi_exact(mesh.centers_x, mesh.centers_y, g)
            sq += volume_weighted_l2(phi[g], ref, mesh.volumes) ** 2
        errors.append(np.sqrt(sq))

    errors = np.asarray(errors)
    orders = np.log2(errors[:-1] / errors[1:])
    if not (orders[-1] > 1.85 and np.all(orders > 1.7)):
        pytest.fail(
            f"2-D LD did not show O(h²) convergence: orders={orders} "
            f"from errors={errors}"
        )
    # Value band (rate ≠ correctness — the manufactured solution is the ref).
    if not (1e-8 < errors[-1] < 1e-2):
        pytest.fail(f"2-D LD finest error out of band: {errors[-1]}")


# ═══════════════════════════════════════════════════════════════════════
# D5b.3 (foundation) — FullFieldWavefront ≡ MovingFrontierWindow (two paths)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_ld_2d_two_paths_ffw_equals_mfw():
    r"""D5b.3: the 2-D LD discretization via the two DAG schedules agrees.

    ``FullFieldWavefront`` (anti-diagonal full cochain) ≡ ``MovingFrontierWindow``
    (rolling frontier) — the SAME LD ``cell_kernel_batch`` via two storage
    policies; the agreement proves the WALK storage doesn't perturb LD's cell
    math (the headline software invariant).  Mode-9 degeneracy-break config: het
    2G, level_symmetric (genuine ``mu_y``), NON-SQUARE (``nx ≠ ny``), vacuum
    edges.  Each leg is FORCED to its INTENDED rep (a two-paths gate where both
    legs secretly run the same rep is a silent false green)."""
    nx, ny = 5, 4
    sn = _nonsquare_het_2g_mesh(nx, ny)
    ng = 2
    N = sn.quad.N
    rng = np.random.default_rng(583)
    sig_t = rng.uniform(0.3, 3.0, size=(ng, nx, ny))
    Q = rng.uniform(0.0, 2.0, size=(N, ng, nx, ny))

    # Each leg builds its OWN rep instance — VERIFY each ran its intended rep.
    mfw = MovingFrontierWindow.pose(sn)
    ffw = FullFieldWavefront.pose(sn)
    if not isinstance(default_for(sn, sn.scheme, sn.angular_closure), MovingFrontierWindow):
        pytest.fail("expected the 2-D LD default rep to be MovingFrontierWindow")

    bf_win = AngularBoundaryFlux.zeros(sn.angular_trace)              # VACUUM (zero domain inflow)
    bf_full = AngularBoundaryFlux.zeros(sn.angular_trace)
    ang_w, scal_w = mfw.sweep(Q, sig_t, bf_win)
    ang_f, scal_f = ffw.sweep(Q, sig_t, bf_full)

    np.testing.assert_allclose(
        ang_w, ang_f, rtol=1e-9, atol=1e-12, err_msg="angular flux FFW≠MFW",
    )
    np.testing.assert_allclose(
        scal_w, scal_f, rtol=1e-9, atol=1e-12, err_msg="scalar flux FFW≠MFW",
    )
    # Single-sweep sharper pin (fixed reduction depth → nulp-class).
    np.testing.assert_array_almost_equal_nulp(scal_w, scal_f, nulp=nx * ny)
    np.testing.assert_allclose(
        bf_win.values, bf_full.values, atol=1e-12,
        err_msg="post-sweep boundary trace FFW≠MFW",
    )


# ═══════════════════════════════════════════════════════════════════════
# D5b.5 (foundation) — DD ≠ LD routing-flip (the cross-thread guard)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_dd_and_ld_2d_converge_to_different_values():
    r"""D5b.5: DD and LD on the SAME 2-D het 2G config converge to DIFFERENT
    fluxes — the catcher that LD is GENUINELY computing LD, not silently
    collapsing to DD (the D5-0 misroute would give DD ≡ LD).

    Drives ``solve_sn_fixed_source`` twice (same mesh/materials/source) at a
    COARSE mesh, where the O(h²) closures diverge, and asserts the converged
    scalar fluxes differ by more than the solver tolerance.  (A test asserting
    DD ≡ LD would be WRONG; a test that only checks "LD runs" is blind to the
    misroute.)"""
    case = build_2d_cartesian_heterogeneous_mms_case()
    mesh = case.build_mesh(12)                      # coarse — closures diverge
    materials = case.build_materials(mesh)
    Q = case.external_source(mesh)
    kw = dict(boundary_condition="vacuum", max_inner=300, inner_tol=1e-11)
    phi_dd = solve_sn_fixed_source(
        materials, mesh, case.quadrature, Q, scheme=DiamondDifference(), **kw,
    ).scalar_flux.values
    phi_ld = solve_sn_fixed_source(
        materials, mesh, case.quadrature, Q, scheme=LinearDiscontinuous(), **kw,
    ).scalar_flux.values

    if not np.all(np.isfinite(phi_ld)):
        pytest.fail("2-D LD produced non-finite flux")
    if np.allclose(phi_dd, phi_ld, rtol=1e-3):
        pytest.fail(
            "DD ≈ LD at coarse mesh — LD is silently computing DD (the D5-0 "
            "misroute regressed, or the bilinear closure was not exercised)"
        )


# ═══════════════════════════════════════════════════════════════════════
# D5b.6 (foundation) — 2-D LD Krylov ≡ SI on a PURE-Z-bearing quadrature
# (ERR-062: the matvec pure-z arm's missing moment-broadcast guard)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.catches("ERR-062")
@pytest.mark.verifies("ld-ubld-pure-z-collision")
def test_ld_2d_krylov_equals_si_pure_z_quadrature():
    r"""D5b.6: 2-D LD ``inner_solver='krylov'`` ≡ ``'source_iteration'`` on a
    quadrature WITH pure-z ordinates (μ_x = μ_y = 0) — the EXACT habitat of the
    qa Concern A blocker (ERR-062).

    The matvec ``loss_action`` pure-z arm ``(L+C)ψ̄ = σ·ψ̄`` lacked the
    moment-axis broadcast guard its sweep twin ``ψ̄ = Q/σ_t`` already had: when
    the LD probe carries the trailing ``2^d`` spatial-moment axis, ``σ·probe``
    broadcast-FAILED (``ValueError: operands could not be broadcast together
    with shapes (ng,nx,ny) (1,ng,nx,ny,2^d)``).  The bug hid on
    ``level_symmetric`` (NO pure-z ordinate → the pure-z arm never runs); it
    fires on a Lebedev / product quadrature that DOES carry the ±z poles.

    The fix single-sources the moment-broadcast through
    ``_moment_broadcast_sigma`` so the two collision-only twins (sweep ≡ matvec)
    cannot diverge (Pattern 2 / L21).  Krylov solves ``(L+C−S_full)`` via the
    matvec; SI via the sweep — both must reach the SAME within-group fixed point,
    so their converged scalar fluxes agree to solver tolerance.

    Mode-9 degeneracy-break config (the genuine cross-check, not a degenerate
    coincidence): a pure-z-bearing **Lebedev** quadrature (genuine ``mu_y``, the
    ±z poles), **heterogeneous** 2-material map (A/B), **2-group asymmetric**
    cross sections with **non-zero self-scatter** (``c > 0`` — the
    ``(L+C−S)`` scattering source is live), **NON-SQUARE** ``nx ≠ ny``, vacuum
    edges (zero domain inflow → no S4 boundary-trace widening).  WITHOUT the fix
    the Krylov solve raises ``ValueError`` at the first matvec; WITH it the two
    paths agree.  ``-O``-safe (``np.testing.*`` / ``pytest.fail``)."""
    nx, ny = 5, 4                                   # NON-SQUARE (x↔y-swap defence)
    mat_map = np.zeros((nx, ny), dtype=int)
    mat_map[nx // 2:, :] = 1                         # heterogeneous: A | B split
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 1.3, nx + 1),
        edges_y=np.linspace(0.0, 0.9, ny + 1),
        mat_map=mat_map,
        coord=CoordSystem.CARTESIAN,
        bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("vacuum"), bc_ymax=BC("vacuum"),
    )
    # Lebedev order 5 (N=14) carries 2 pure-z ordinates (the ±z poles) AND
    # genuine mu_y — the cheap pure-z habitat (the production 2-D Cartesian MMS
    # uses Lebedev order 17 / N=110 with the same 2 pure-z ordinates).
    quad = Quadrature.lebedev(order=5)
    is_pure_z = (np.abs(quad.mu_x) < 1e-12) & (np.abs(quad.mu_y) < 1e-12)
    if not np.any(is_pure_z):
        pytest.fail(
            "test premise broken: the quadrature has NO pure-z ordinate, so the "
            "pure-z matvec arm (ERR-062's habitat) is never exercised"
        )

    materials = {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
    ng = 2
    N = quad.N
    rng = np.random.default_rng(583)
    Q = rng.uniform(0.0, 2.0, size=(N, ng, nx, ny))

    kw = dict(
        boundary_condition="vacuum", max_inner=500, inner_tol=1e-12,
        scheme=LinearDiscontinuous(),
    )
    # WITHOUT the fix this raises ValueError on the first Krylov matvec.
    phi_kr = solve_sn_fixed_source(
        materials, mesh, quad, Q, inner_solver="krylov", **kw,
    ).scalar_flux.values
    phi_si = solve_sn_fixed_source(
        materials, mesh, quad, Q, inner_solver="source_iteration", **kw,
    ).scalar_flux.values

    if not np.all(np.isfinite(phi_kr)):
        pytest.fail("2-D LD Krylov produced non-finite flux")
    # Same (L+C−S_full) fixed point reached by two structurally-distinct inners.
    np.testing.assert_allclose(
        phi_kr, phi_si, rtol=1e-9, atol=1e-10,
        err_msg="2-D LD Krylov ≠ SI on the pure-z quadrature (ERR-062)",
    )


# ═══════════════════════════════════════════════════════════════════════
# D5b.2 (HEADLINE) — 2-D LD STRESS MMS O(h²) + value band, on FullFieldWavefront
# ═══════════════════════════════════════════════════════════════════════
#
# The vv-Mode-7 strengthened ansatz (SN2DCartesianLDStressMMSCase): μ-bilinear
# (μ_x·B + μ_y·C activates the per-axis bilinear slope rows), NON-vanishing
# boundary trace (a0>0 → prescribed inflow), NON-SQUARE het 2G.  Unlike the
# vacuum-edge smoke gate above, this is the flux-shape VERIFICATION claim
# (@verifies("ld-cartesian-2d")), not just a convergence check.  See the
# SN2DCartesianLDStressMMSCase docstring for the HONEST SCOPE: this closes the
# slope-UNKNOWN half of the LM-1989 trap; the slope-SOURCE half is deferred
# (the production lifts only the AVERAGE-moment external source — solver.py
# `_lift_external_source_to_moments`).


def _ld_stress_l2_errors(case, n_cells):
    r"""Volume-weighted L2 error of the LD stress solve at each mesh, run on the
    NON-vacuum prescribed-inflow source (the FullFieldWavefront oracle leg)."""
    errors = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        materials = case.build_materials(mesh)
        sn = SNMesh(mesh, case.quadrature, materials, scheme=LinearDiscontinuous())
        rhs = build_nonvacuum_fixed_source(case, sn)
        result = solve_sn_fixed_source(
            materials, mesh, case.quadrature, rhs,
            max_inner=500, inner_tol=1e-12, scheme=LinearDiscontinuous(),
        )
        phi = result.scalar_flux.values                 # (ng, nx, ny)
        sq = 0.0
        for g in range(phi.shape[0]):
            ref = case.phi_exact(mesh.centers_x, mesh.centers_y, g)
            sq += volume_weighted_l2(phi[g], ref, mesh.volumes) ** 2
        errors.append(np.sqrt(sq))
    return np.asarray(errors)


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("ld-cartesian-2d", "transport-cartesian-2d")
def test_ld_2d_stress_converges_second_order():
    r"""D5b.2: the 2-D bilinear LD closure converges O(h²) to the manufactured
    STRESS solution :math:`A_g(x,y)`, with a non-vanishing prescribed-inflow
    boundary trace — the headline flux-shape verification of the multi-D slope
    rows.

    The strengthened ansatz :math:`\psi = (A + \mu_x B + \mu_y C)/W` activates
    the per-axis bilinear slope rows (DD cannot represent the μ-weighted slope);
    :math:`a_0>0` stresses the boundary closure at the face-average moment;
    NON-SQUARE het 2G is the x↔y-swap + multigroup discrimination.  Asserts the
    observed order climbs to 2 (finest > 1.95, all > 1.85) AND the converged
    VALUE lands in a band against ``phi_exact`` (vv §5: rate ≠ correctness — the
    manufactured solution IS the structurally-independent reference).

    HONEST SCOPE: this gate verifies the slope-UNKNOWN sign + the AVERAGE-moment
    boundary closure on the FLAT external source (the slope rows are zero here,
    so it is BLIND to the slope-SOURCE sign — the Mode-10 gap).  The slope-SOURCE
    sign half of the LM-1989 trap is now VERIFIED separately (Leg A, #247) by the
    per-moment structural gate + mutation controls in the #247 block below; the
    BOUNDARY transverse-face-slope (Leg B) is tracked in #251.  See the case
    docstring.  ``-O``-safe (``np.testing.*`` / ``pytest.fail``)."""
    case = build_2d_cartesian_ld_stress_mms_case()
    errors = _ld_stress_l2_errors(case, [16, 32, 64])
    orders = np.log2(errors[:-1] / errors[1:])
    if not (orders[-1] > 1.95 and np.all(orders > 1.85)):
        pytest.fail(
            f"2-D LD stress MMS did not show O(h²) convergence: orders={orders} "
            f"from errors={errors}"
        )
    # Value band (rate ≠ correctness — the manufactured solution is the ref).
    if not (1e-9 < errors[-1] < 1e-2):
        pytest.fail(f"2-D LD stress finest error out of band: {errors[-1]}")


# ═══════════════════════════════════════════════════════════════════════
# D5b.4 — 2-D LD Krylov ≡ SI on the STRESS config (the matvec-twin leg, L14)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("ld-cartesian-2d")
def test_ld_2d_stress_krylov_equals_si():
    r"""D5b.4: 2-D LD ``inner_solver='krylov'`` ≡ ``'source_iteration'`` on the
    STRESS config — the L14 matvec-twin leg on the slope-active habitat.

    Krylov drives ``(L+C−S_full)`` via ``loss_action`` (``residual_kernel_batch``,
    the apply twin); SI via the sweep (``cell_kernel_batch``).  Both must reach
    the SAME within-group fixed point on the strengthened μ-bilinear, het-2G,
    non-square, NON-vacuum config (the existing pure-z gate covers a different
    habitat — collision-only ordinates; THIS gate covers the bilinear-slope +
    non-vanishing-trace habitat).  Their converged scalar fluxes agree to solver
    tolerance.  Mode-9 degeneracy-break: het 2G with live self-scatter,
    μ-non-trivial, non-square, non-vanishing inflow.  ``-O``-safe."""
    case = build_2d_cartesian_ld_stress_mms_case()
    mesh = case.build_mesh(16)
    materials = case.build_materials(mesh)
    sn = SNMesh(mesh, case.quadrature, materials, scheme=LinearDiscontinuous())
    rhs = build_nonvacuum_fixed_source(case, sn)
    kw = dict(max_inner=500, inner_tol=1e-12, scheme=LinearDiscontinuous())

    phi_kr = solve_sn_fixed_source(
        materials, mesh, case.quadrature, rhs, inner_solver="krylov", **kw,
    ).scalar_flux.values
    phi_si = solve_sn_fixed_source(
        materials, mesh, case.quadrature, rhs,
        inner_solver="source_iteration", **kw,
    ).scalar_flux.values

    if not np.all(np.isfinite(phi_kr)):
        pytest.fail("2-D LD stress Krylov produced non-finite flux")
    np.testing.assert_allclose(
        phi_kr, phi_si, rtol=1e-9, atol=1e-11,
        err_msg="2-D LD stress Krylov ≠ SI (matvec twin disagrees with sweep)",
    )


# ═══════════════════════════════════════════════════════════════════════
# D5b.3 (STRESS upgrade) — FullFieldWavefront ≡ MovingFrontierWindow on the
# strengthened ansatz (alongside the existing random-source vacuum pin)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("ld-cartesian-2d")
def test_ld_2d_stress_two_paths_ffw_equals_mfw():
    r"""D5b.3 (stress upgrade): the 2-D LD discretization via the two DAG
    schedules agrees on the STRESS source.

    ``FullFieldWavefront`` (anti-diagonal full cochain) ≡ ``MovingFrontierWindow``
    (rolling frontier) — the SAME LD ``cell_kernel_batch`` via two storage
    policies, now driven by the manufactured strengthened-ansatz source (not the
    random vacuum source of ``test_ld_2d_two_paths_ffw_equals_mfw``).  The
    agreement proves the WALK storage doesn't perturb LD's cell math on the
    real slope-active habitat.  Each leg is FORCED to its INTENDED rep (a
    two-paths gate where both legs secretly run the same rep is a silent false
    green).  The non-vanishing boundary trace is seeded identically on both.

    NOTE this drives the LOW-LEVEL ``rep.sweep`` directly (not the full solve)
    so each schedule is pinned explicitly; the source is the manufactured
    per-ordinate field on a coarse mesh (one sweep, fixed reduction depth →
    nulp-class agreement)."""
    case = build_2d_cartesian_ld_stress_mms_case()
    mesh = case.build_mesh(10)                       # coarse, NON-SQUARE (10×7)
    materials = case.build_materials(mesh)
    sn = SNMesh(mesh, case.quadrature, materials, scheme=LinearDiscontinuous())
    nx, ny = mesh.mat_map.shape
    ng = case.n_groups
    Q = case.external_source(mesh)                   # (N, ng, nx, ny)
    sig_t = np.stack([
        case.sigma_t_fn(*np.meshgrid(mesh.centers_x, mesh.centers_y,
                                     indexing="ij"), g)
        for g in range(ng)
    ])                                               # (ng, nx, ny)

    mfw = MovingFrontierWindow.pose(sn)
    ffw = FullFieldWavefront.pose(sn)
    if not isinstance(default_for(sn, sn.scheme, sn.angular_closure), MovingFrontierWindow):
        pytest.fail("expected the 2-D LD default rep to be MovingFrontierWindow")

    # Non-vanishing prescribed inflow → a seeded AngularBoundaryFlux on both legs.
    bss = case.prescribed_inflow(sn)
    bf_win = AngularBoundaryFlux.zeros(sn.angular_trace)
    bf_full = AngularBoundaryFlux.zeros(sn.angular_trace)
    bf_win.values[...] = bss.values
    bf_full.values[...] = bss.values

    ang_w, scal_w = mfw.sweep(Q, sig_t, bf_win)
    ang_f, scal_f = ffw.sweep(Q, sig_t, bf_full)

    np.testing.assert_allclose(
        ang_w, ang_f, rtol=1e-9, atol=1e-12, err_msg="angular flux FFW≠MFW",
    )
    np.testing.assert_allclose(
        scal_w, scal_f, rtol=1e-9, atol=1e-12, err_msg="scalar flux FFW≠MFW",
    )
    np.testing.assert_array_almost_equal_nulp(scal_w, scal_f, nulp=nx * ny)


# ═══════════════════════════════════════════════════════════════════════
# #247 — slope-SOURCE half of the LM-1989 trap (the Mode-10 closeout).
#
# LIVE gates — were xfail-strict / skip until the production widening LANDED
# (solve_sn_fixed_source accepts a moment-resolved bulk; the lift passes
# the slope rows; a moment-resolved boundary trace).  See
# `.claude/agent-memory/test-architect/issue_247_moment_source_gate_spec.md`
# for the full design + the PROBED Mode-10 evidence (the slope-SOURCE sign
# flip is O(h²)-small → the teeth are STRUCTURAL, not the converged flux).
#
# Projection convention (the CRUX, locked): the projected moment is the
# BARE per-volume tensor-Legendre coefficient (q̄ = cell average, slot 0;
# q̂_a = ⟨q,P₁⟩/⟨P₁,P₁⟩ on axis a — NO θ, NO h, NO V — the kernel's mass
# M = diag(h, θh) adds them).  d=2 Kronecker order [ψ̄, ψ̂_y, ψ̂_x, ψ̂_xy]
# (x slope = SLOT 2).  GLOBAL frame (production reframes per-octant).  The
# projector calls leggauss directly — NEVER _lift_external_source_to_moments
# nor any LD cell op (L11 structural independence / non-tautology).
# ═══════════════════════════════════════════════════════════════════════

_D2_MOMENT_SLOTS = {(0, 0): 0, (0, 1): 1, (1, 0): 2, (1, 1): 3}  # [bar,y,x,xy]


def _project_scalar_to_tensor_legendre(mesh, fn, *, ng, q_nodes=6, per_ord=None):
    r"""Project a scalar field ``fn(x_grid, y_grid, g) -> (len(x),len(y))`` (or
    ``(per_ord, len(x), len(y))``) onto the per-cell tensor-Legendre moments.

    Structurally INDEPENDENT of production: uses only
    :func:`numpy.polynomial.legendre.leggauss`.  Returns the BARE per-volume
    Legendre coefficients ``q̄`` (slot 0 = cell average), ``q̂_a`` (slot a =
    ``⟨q,P₁⟩/⟨P₁,P₁⟩``), in the d=2 Kronecker order ``[bar, y, x, xy]``.  This
    is exactly the ``S_moments`` the LD cell op consumes (the kernel's mass
    ``M=diag(h,θh)`` applies the volume/θ weighting — the projection must NOT).
    """
    from numpy.polynomial.legendre import leggauss

    xi, wq = leggauss(q_nodes)
    W2 = wq.sum()
    leg = {0: np.ones_like(xi), 1: xi}
    mean_p2 = {0: 1.0, 1: float((wq * xi * xi).sum() / W2)}  # mean(P1²)=1/3
    ex, ey = mesh.edges_x, mesh.edges_y
    cx, cy = mesh.centers_x, mesh.centers_y
    nx, ny = len(cx), len(cy)
    lead = () if per_ord is None else (per_ord,)
    out = np.zeros((*lead, ng, nx, ny, 4))
    for i in range(nx):
        xq = (ex[i] + ex[i + 1]) / 2 + (ex[i + 1] - ex[i]) / 2 * xi
        for j in range(ny):
            yq = (ey[j] + ey[j + 1]) / 2 + (ey[j + 1] - ey[j]) / 2 * xi
            for g in range(ng):
                F = fn(xq, yq, g)  # (q,q) or (per_ord,q,q)
                for ox in (0, 1):
                    for oy in (0, 1):
                        wgt = np.outer(wq * leg[ox], wq * leg[oy]) / (W2 * W2)
                        coeff = np.einsum("...ij,ij->...", F, wgt) / (
                            mean_p2[ox] * mean_p2[oy]
                        )
                        out[..., g, i, j, _D2_MOMENT_SLOTS[(ox, oy)]] = coeff
    return out


# ───────────────────────────────────────────────────────────────────────
# DELIVERABLE 3 — projection-correctness sub-gate (L11, foundation)
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.verifies("ld-cartesian-2d-bilinear-coeffs")
def test_tensor_legendre_projection_matches_hand_polynomial():
    r"""#247 sub-gate: the quadrature projector reproduces the EXACT
    tensor-Legendre coefficients of a KNOWN bilinear polynomial.

    ``q = a00 + a10 x + a01 y + a11 xy`` on a cell has the hand-derivable
    coefficients (Legendre ``{1,ξ}``, ``ξ∈[-1,1]``, the SAME basis as the UBLD
    mass ``diag(h,θh)``):

        q̄    = a00 + a10·xc + a01·yc + a11·xc·yc      (slot 0)
        q̂_y  = (hy/2)(a01 + a11·xc)                   (slot 1)
        q̂_x  = (hx/2)(a10 + a11·yc)                   (slot 2)
        q̂_xy = (hx/2)(hy/2)·a11                       (slot 3)

    A bilinear integrand is integrated EXACTLY by ``q_nodes>=2``, so the match
    is machine precision.  Non-tautology (L11): the projector uses only
    ``leggauss`` + the hand-derived integral — NEVER
    ``_lift_external_source_to_moments`` nor any LD cell op; the reference is
    hand-laid polynomial algebra, not a production echo.  ``-O``-safe."""
    from numpy.polynomial.legendre import leggauss

    xL, xR, yL, yR = 1.0, 1.4, 2.0, 2.6
    a00, a10, a01, a11 = 0.7, 1.3, -0.5, 2.1
    hx, hy = xR - xL, yR - yL
    xc, yc = (xL + xR) / 2, (yL + yR) / 2

    # Inline single-cell projector (the helper above is mesh-wide; here we want
    # one cell to keep the hand reference trivial).
    xi, wq = leggauss(4)
    W2 = wq.sum()
    leg = {0: np.ones_like(xi), 1: xi}
    mean_p2 = {0: 1.0, 1: float((wq * xi * xi).sum() / W2)}
    xq = xc + (hx / 2) * xi
    yq = yc + (hy / 2) * xi
    XX, YY = np.meshgrid(xq, yq, indexing="ij")
    F = a00 + a10 * XX + a01 * YY + a11 * XX * YY
    got = {}
    for ox in (0, 1):
        for oy in (0, 1):
            wgt = np.outer(wq * leg[ox], wq * leg[oy]) / (W2 * W2)
            got[(ox, oy)] = float(np.sum(F * wgt) / (mean_p2[ox] * mean_p2[oy]))

    expect = {
        (0, 0): a00 + a10 * xc + a01 * yc + a11 * xc * yc,
        (0, 1): (hy / 2) * (a01 + a11 * xc),
        (1, 0): (hx / 2) * (a10 + a11 * yc),
        (1, 1): (hx / 2) * (hy / 2) * a11,
    }
    for key, slot_name in [
        ((0, 0), "average"), ((1, 0), "x-slope"),
        ((0, 1), "y-slope"), ((1, 1), "xy"),
    ]:
        np.testing.assert_allclose(
            got[key], expect[key], rtol=0, atol=1e-13,
            err_msg=f"projector {slot_name} coeff wrong (normalization drift)",
        )


@pytest.mark.foundation
@pytest.mark.verifies("ld-cartesian-2d-projection-coeff")
def test_projection_slot0_is_cell_average_not_centre():
    r"""#247 sub-gate (2nd leg): the projector's slot-0 IS the cell AVERAGE,
    not the cell-CENTRE value (the §1 normalization caveat).

    The existing flat ``case.external_source`` evaluates Q at the cell centre;
    the projector cell-averages.  They differ by O(h²).  Pin that slot-0
    equals an INDEPENDENT fine-quadrature (q_nodes=12) cell average to ~1e-8
    (so the projector is doing genuine averaging) — and do NOT cross-check
    slot-0 against ``external_source`` (which would falsely fail by O(h²)).
    ``-O``-safe."""
    case = build_2d_cartesian_ld_stress_mms_case()
    mesh = case.build_mesh(8)
    a_only = lambda x, y, g: case._drivers(x, y, g)[0]  # noqa: E731
    coarse = _project_scalar_to_tensor_legendre(mesh, a_only, ng=case.n_groups, q_nodes=6)
    fine = _project_scalar_to_tensor_legendre(mesh, a_only, ng=case.n_groups, q_nodes=12)
    np.testing.assert_allclose(
        coarse[..., 0], fine[..., 0], rtol=0, atol=1e-7,
        err_msg="projector slot-0 is not the (quadrature-converged) cell average",
    )


# ───────────────────────────────────────────────────────────────────────
# Shared #247 plumbing — the manufactured moment source + the driver that runs
# it through the WIDENED public solve (``_solve_moment_resolved``; the #247
# production path has landed).
# ───────────────────────────────────────────────────────────────────────


def _manufactured_Q_pointwise(case, x, y, g, quad):
    r"""Closed-form ``Q^ext_{n,g}(x,y)`` per ordinate, on a matched (x,y) grid
    — the SAME residual ``case.external_source`` emits, but pointwise (for
    quadrature projection).  Per-ordinate density (÷ sum_w)."""
    mu_x, mu_y = quad.mu_x, quad.mu_y
    sum_w = float(quad.weights.sum())
    N, ng = len(mu_x), case.n_groups
    A, dAx, dAy, B, dBx, dBy, C, dCx, dCy = case._drivers(x, y, g)
    sig = case.sigma_t_fn(np.repeat(x, len(y)), np.tile(y, len(x)), g).reshape(
        len(x), len(y)
    )
    in_sc = np.zeros((len(x), len(y)))
    for gf in range(ng):
        ss = case.sigma_s_fn(
            np.repeat(x, len(y)), np.tile(y, len(x)), gf, g
        ).reshape(len(x), len(y))
        in_sc += ss * case._drivers(x, y, gf)[0]
    Q = np.zeros((N, len(x), len(y)))
    for n in range(N):
        stream = (
            mu_x[n] * dAx + mu_y[n] * dAy
            + mu_x[n] ** 2 * dBx + mu_x[n] * mu_y[n] * dBy
            + mu_x[n] * mu_y[n] * dCx + mu_y[n] ** 2 * dCy
        )
        Q[n] = (stream + sig * (A + mu_x[n] * B + mu_y[n] * C) - in_sc) / sum_w
    return Q


def _project_external_source(case, mesh):
    r"""Project the manufactured per-ordinate ``Q^ext`` onto the GLOBAL-frame
    tensor-Legendre moment vector ``(N, ng, nx, ny, 4)`` — the moment-resolved
    external source the production widening consumes (#247)."""
    quad = case.quadrature
    fn = lambda x, y, g: _manufactured_Q_pointwise(case, x, y, g, quad)  # noqa: E731
    return _project_scalar_to_tensor_legendre(
        mesh, fn, ng=case.n_groups, q_nodes=6, per_ord=quad.N,
    )


def _moment_resolved_rhs(case, sn, mesh, *, moment_source=None):
    r"""Build the composite RHS whose BULK is the moment-resolved external
    source (``(N, ng, nx, ny, 4)`` — slot 0 = average, slope rows = projected
    :math:`\hat Q`) and whose boundary is the manufactured prescribed inflow.

    Drives the WIDENED public ``solve_sn_fixed_source`` (the production #247
    path): ``_build_fixed_source_rhs`` accepts the moment-resolved bulk and
    ``_lift_external_source_to_moments`` threads the slope rows through.  The
    moment source is computed by the structurally-independent projector
    (:func:`_project_external_source` → ``leggauss`` only); pass an explicit
    ``moment_source`` to override (the mutation controls flip a slope row)."""
    from orpheus.transport.source_sinks import AngularSourceSink
    from orpheus.transport.timed_full_field import TimedFullField

    Qm = _project_external_source(case, mesh) if moment_source is None else moment_source
    return TimedFullField(
        interior=AngularSourceSink(values=Qm, space=sn.angular_trial_space),
        boundary=case.prescribed_inflow(sn),
    )


def _solve_moment_resolved(case, nc, *, moment_source=None):
    r"""Run the full public solve on the moment-resolved external source at mesh
    ``nc`` and return ``(scalar_flux_values, mesh)`` — the WIDENED #247 path."""
    mesh = case.build_mesh(nc)
    materials = case.build_materials(mesh)
    sn = SNMesh(mesh, case.quadrature, materials, scheme=LinearDiscontinuous())
    rhs = _moment_resolved_rhs(case, sn, mesh, moment_source=moment_source)
    result = solve_sn_fixed_source(
        materials, mesh, case.quadrature, rhs,
        max_inner=500, inner_tol=1e-12, scheme=LinearDiscontinuous(),
    )
    return result.scalar_flux.values, mesh


def _A_l2_err(case, phi, mesh):
    r"""Volume-weighted L2 error of a scalar flux ``phi`` (``(ng, nx, ny)``) vs
    the manufactured average :math:`\phi = A` over all groups."""
    sq = 0.0
    for g in range(phi.shape[0]):
        ref = case.phi_exact(mesh.centers_x, mesh.centers_y, g)
        sq += volume_weighted_l2(phi[g], ref, mesh.volumes) ** 2
    return float(np.sqrt(sq))


# ───────────────────────────────────────────────────────────────────────
# DELIVERABLE — the STRUCTURAL teeth: the production lift threads the
# moment-resolved slope rows through UNCHANGED (machine precision).  This is
# the sharpest, non-tautological proof of the actual #247 production change
# (the slope rows are no longer zeroed); a regression that re-zeroes them is
# caught at O(1) here, where the converged flux is only sub-floor sensitive
# (spec §0).  L11: the projected reference is built by leggauss only — it never
# calls _lift_external_source_to_moments nor any LD cell op.
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_ld_2d_external_slope_source_threaded_through_lift():
    r"""#247 STRUCTURAL teeth — the production lift threads a moment-resolved
    external source through UNCHANGED (the slope rows :math:`\hat Q` are no
    longer zeroed).

    The actual production change is ``_lift_external_source_to_moments`` going
    from "zero the 2^d buffer + slot-0 ← flat (slopes ZERO)" to "thread a
    moment-resolved bulk through".  This pins it directly + at machine precision:
    feed the projected moment vector to the lift and assert the returned moment
    source equals the projection EXACTLY (every slope slot, not just slot 0).  A
    regression that re-zeroes the slope rows breaks this O(1) (the spec §0 trap:
    the converged flux would NOT catch it — its slope-source sensitivity is
    sub-floor).  Non-tautology (L11): the projection is built by ``leggauss``
    only; the SECOND leg (the flat path lifts onto slot-0 with slopes ZERO) is
    the negative control that pins the honest default is preserved.  ``-O``-safe.
    """
    from orpheus.sn.solver import _lift_external_source_to_moments

    case = build_2d_cartesian_ld_stress_mms_case()
    mesh = case.build_mesh(8)
    materials = case.build_materials(mesh)
    sn = SNMesh(mesh, case.quadrature, materials, scheme=LinearDiscontinuous())
    Qm = _project_external_source(case, mesh)        # (N, ng, nx, ny, 4)

    lifted, per_axis = _lift_external_source_to_moments(Qm, sn)
    if per_axis != 2:
        pytest.fail(f"expected LD per_axis == 2, got {per_axis}")
    # The whole moment vector — slot 0 AND the slope rows — threads UNCHANGED.
    np.testing.assert_array_equal(
        lifted, Qm,
        err_msg="moment-resolved external source not threaded through the lift "
        "unchanged (the #247 slope rows were dropped/zeroed)",
    )
    # NEGATIVE CONTROL: a FLAT bulk still lifts onto slot 0 with slopes ZERO
    # (the honest default — the backward-compat invariant must not move).
    flat = case.external_source(mesh)                 # (N, ng, nx, ny)
    flat_lifted, _ = _lift_external_source_to_moments(flat, sn)
    np.testing.assert_array_equal(
        flat_lifted[..., 0], flat,
        err_msg="flat external source slot-0 not preserved (regression)",
    )
    np.testing.assert_array_equal(
        flat_lifted[..., 1:], 0.0,
        err_msg="flat external source slope rows are no longer ZERO (the honest "
        "default Q̂=0 regressed)",
    )


# ───────────────────────────────────────────────────────────────────────
# DELIVERABLE 1 — per-moment / converged O(h²) convergence with the
# slope-SOURCE consumed (NECESSARY but per §0 NOT SUFFICIENT for the SIGN).
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("ld-cartesian-2d")
def test_ld_2d_external_slope_source_converges_second_order():
    r"""#247 Deliverable 1 — the WIDENED public solve, fed the moment-resolved
    external source (the projected slope rows :math:`\hat Q` consumed), converges
    O(h²) to the manufactured :math:`\phi = A`.

    Proves the slope rows are CONSISTENT: the slope-UNKNOWN + the slope-SOURCE
    together produce a 2nd-order solution.  NECESSARY (a wrong-MAGNITUDE slope
    source would break the order) but per §0 NOT SUFFICIENT for the SIGN — a sign
    flip leaves the order at ~2 (the slope-SOURCE error is O(h²)-small, the
    Mode-10 trap); the SIGN teeth are the mutation control below.  ``-O``-safe.
    """
    case = build_2d_cartesian_ld_stress_mms_case()
    n_cells = [16, 32, 64]
    errors = np.array([
        _A_l2_err(case, *_solve_moment_resolved(case, nc)) for nc in n_cells
    ])
    orders = np.log2(errors[:-1] / errors[1:])
    if not (orders[-1] > 1.9 and np.all(orders > 1.8)):
        pytest.fail(
            f"moment-resolved external source did not converge O(h²): "
            f"orders={orders} from errors={errors}"
        )
    if not (1e-9 < errors[-1] < 1e-2):
        pytest.fail(f"moment-resolved finest error out of band: {errors[-1]}")


@pytest.mark.l1
@pytest.mark.verifies("ld-cartesian-2d")
def test_ld_2d_external_slope_source_improves_on_flat():
    r"""#247 — the projected slope-SOURCE rows IMPROVE the converged flux vs the
    flat (slope-zeroed) external source.

    A POSITIVE verification that the consumed :math:`\hat Q` carries real
    sub-cell information: at a fixed mesh the moment-resolved solve lands strictly
    closer to the manufactured :math:`\phi = A` than the flat-in-moment solve
    (probed ≈ 3.4e-3 vs 5.9e-3 at nc=24).  Distinct from a value-band: it
    compares two production runs that differ ONLY in whether the slope rows are
    threaded, so the discretization floor cancels and the slope-source
    contribution is the signal.  ``-O``-safe."""
    case = build_2d_cartesian_ld_stress_mms_case()
    nc = 24
    phi_mom, mesh = _solve_moment_resolved(case, nc)
    # The flat (slope-zeroed) sibling — same boundary trace, same average source.
    materials = case.build_materials(mesh)
    sn = SNMesh(mesh, case.quadrature, materials, scheme=LinearDiscontinuous())
    rhs_flat = build_nonvacuum_fixed_source(case, sn)
    phi_flat = solve_sn_fixed_source(
        materials, mesh, case.quadrature, rhs_flat,
        max_inner=500, inner_tol=1e-12, scheme=LinearDiscontinuous(),
    ).scalar_flux.values

    err_mom = _A_l2_err(case, phi_mom, mesh)
    err_flat = _A_l2_err(case, phi_flat, mesh)
    if not (err_mom < err_flat):
        pytest.fail(
            f"moment-resolved external source ({err_mom:.3e}) did NOT improve on "
            f"the flat source ({err_flat:.3e}) — the slope rows carry no signal"
        )


# ───────────────────────────────────────────────────────────────────────
# DELIVERABLE 2 — the mutation control (anti-pattern #11, the PRIMARY
# sign-catcher).  M1–M3 (EXTERNAL Q̂) + M4 (SCATTERING Σ_s·φ̂).  The teeth:
# flipping a CONSUMED source slope row changes the converged flux by ≫ solver
# tol (the consumption proof — O(1) above the 1e-12 fixed point), while the
# FLAT scalar gate (test_ld_2d_stress_converges_second_order) stays GREEN
# because it feeds an already-zero slope row (the Mode-10 asymmetry).
# ───────────────────────────────────────────────────────────────────────

#: A converged-flux change this far above the inner tolerance proves the flipped
#: source row was CONSUMED (not carried-and-ignored).  The inner solve converges
#: to 1e-12; the smallest probed slope flip moves the flux ~6e-5 (the xy slot) —
#: ~5e7× above tol.  The band is well clear of FP / iteration noise yet far below
#: the spec §0 trap (a value-band against A would need a tolerance tighter than
#: the ~1.5× discretization gap, which the O(h²) floor eats).
_CONSUMPTION_TOL = 1e-8


@pytest.mark.l1
@pytest.mark.verifies("ld-cartesian-2d")
@pytest.mark.parametrize("slot,name", [(2, "x-slope"), (1, "y-slope"), (3, "xy")])
def test_ld_2d_external_slope_source_sign_mutation_reddens(slot, name):
    r"""#247 Deliverable 2 (M1–M3) — flipping the EXTERNAL :math:`\hat Q`
    slope-row sign (slot 2=x, 1=y, 3=xy) CHANGES the converged flux (the
    consumption proof) while the FLAT scalar-flux gate stays GREEN.

    TEETH (the consumption proof, O(1) above the 1e-12 fixed point — NOT the
    sub-floor converged-flux value-band of spec §0): a sign flip on the
    CONSUMED slope-source row moves the converged scalar flux by ≫
    ``_CONSUMPTION_TOL`` (probed ~3e-3 x-slope, ~1e-2 y-slope, ~6e-5 xy at
    nc=24).  The asymmetry that closes the Mode-10 gap: the FLAT scalar gate
    ``test_ld_2d_stress_converges_second_order`` feeds an already-zero slope row,
    so a flip there is a no-op (verified directly below by the FLAT no-op leg) —
    the flat gate stays GREEN under the same flip.  ``-O``-safe."""
    case = build_2d_cartesian_ld_stress_mms_case()
    nc = 24
    phi_ok, mesh = _solve_moment_resolved(case, nc)

    Qm = _project_external_source(case, mesh)
    Qm_flip = Qm.copy()
    Qm_flip[..., slot] *= -1.0
    phi_flip, _ = _solve_moment_resolved(case, nc, moment_source=Qm_flip)

    change = np.abs(phi_flip - phi_ok).max() / np.abs(phi_ok).max()
    if not (change > _CONSUMPTION_TOL):
        pytest.fail(
            f"flipping the EXTERNAL {name} slope-source row did NOT change the "
            f"converged flux (|Δφ|/|φ|={change:.3e} ≤ {_CONSUMPTION_TOL:.0e}) — "
            f"the slope row is NOT consumed (the #247 lift regressed to zeroing)"
        )

    # The Mode-10 asymmetry: the SAME flip on the FLAT external source is a
    # no-op (the flat path zeroes the slope row, so flipping zero changes
    # nothing) — the flat scalar gate is correctly BLIND to the slope-source
    # sign, which is exactly why #247 adds the moment-resolved gate.
    materials = case.build_materials(mesh)
    sn = SNMesh(mesh, case.quadrature, materials, scheme=LinearDiscontinuous())
    from orpheus.sn.solver import _lift_external_source_to_moments
    flat = case.external_source(mesh)
    flat_lift, _ = _lift_external_source_to_moments(flat, sn)
    np.testing.assert_array_equal(
        flat_lift[..., slot], 0.0,
        err_msg=f"the FLAT external source {name} slope row is non-zero — the "
        "flat gate would NOT be blind to the flip (the Mode-10 asymmetry breaks)",
    )


@pytest.mark.l1
@pytest.mark.verifies("ld-cartesian-2d")
def test_ld_2d_scattering_slope_source_sign_mutation_reddens(monkeypatch):
    r"""#247 Deliverable 2 (M4) — flipping the SCATTERING :math:`\Sigma_s\cdot
    \hat\phi` slope-row sign CHANGES the converged flux.

    Distinct from M1–M3: this verifies the EXISTING (S3) scattering consumption
    was never sign-blind, whereas M1–M3 verify the NEW external consumption.
    Monkeypatch the per-ordinate scattering-source combine
    (:meth:`AngularLift._combine <orpheus.transport.operators.angular_lift.AngularLift._combine>`
    — the producer-side ``(iso/W) + aniso`` over every spatial moment, the
    moment-carrying source that joins the external :math:`\hat Q`; until CS4c
    step 5 spelled ``ScatteringOperator._assemble_per_ordinate_source``) to
    negate the iso slope rows (slots 1:), reverted by the fixture.  TEETH =
    the consumption proof (probed |Δφ|/|φ| ≈ 2.6e-3 at nc=24, ≫
    ``_CONSUMPTION_TOL``).  ``-O``-safe."""
    from orpheus.transport.operators.angular_lift import AngularLift

    case = build_2d_cartesian_ld_stress_mms_case()
    nc = 24
    mesh = case.build_mesh(nc)
    materials = case.build_materials(mesh)
    sn = SNMesh(mesh, case.quadrature, materials, scheme=LinearDiscontinuous())
    rhs = build_nonvacuum_fixed_source(case, sn)
    kw = dict(max_inner=500, inner_tol=1e-12, scheme=LinearDiscontinuous())

    phi_ok = solve_sn_fixed_source(
        materials, mesh, case.quadrature, rhs, **kw,
    ).scalar_flux.values

    orig = AngularLift._combine

    def _flip_iso_slopes(self, iso, aniso):
        out = orig(self, iso, aniso)
        # Moment-valued per-ordinate source → flip every slope row (slots 1:).
        if out.values.ndim >= 5 and out.values.shape[-1] > 1:
            out.values[..., 1:] *= -1.0
        return out

    monkeypatch.setattr(AngularLift, "_combine", _flip_iso_slopes)
    phi_flip = solve_sn_fixed_source(
        materials, mesh, case.quadrature, rhs, **kw,
    ).scalar_flux.values
    monkeypatch.undo()

    change = np.abs(phi_flip - phi_ok).max() / np.abs(phi_ok).max()
    if not (change > _CONSUMPTION_TOL):
        pytest.fail(
            f"flipping the SCATTERING slope-source rows did NOT change the "
            f"converged flux (|Δφ|/|φ|={change:.3e} ≤ {_CONSUMPTION_TOL:.0e}) — "
            f"the scattering slope source Σ_s·φ̂ is sign-blind"
        )


# ───────────────────────────────────────────────────────────────────────
# DELIVERABLE 5 — the negative pin (the relaxed bulk-shape check still rejects
# a wrong trailing-moment axis, so the typed-union relaxation does not swallow
# real shape bugs).
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_moment_resolved_bulk_still_rejects_wrong_trailing_axis():
    r"""#247 Deliverable 5 — the widened bulk-shape check (typed union: flat OR
    moment-resolved bulk) STILL rejects a trailing-moment axis ≠ ``per_axis**ndim
    = 4`` (the moment-vector contract), so the relaxation does not swallow real
    shape bugs.  Also pins that DD/Step (``per_axis == 1``, no moment axis)
    rejects ANY moment-resolved bulk (only flat is valid there).  ``-O``-safe
    (``pytest.raises``)."""
    from orpheus.transport.spatial import DiamondDifference

    case = build_2d_cartesian_ld_stress_mms_case()
    mesh = case.build_mesh(8)
    materials = case.build_materials(mesh)
    N = case.quadrature.N
    ng = case.n_groups
    nx, ny = mesh.mat_map.shape

    # LD: a trailing axis of length 5 (≠ 2^d = 4) is a real shape bug.
    bad_width = np.zeros((N, ng, nx, ny, 5))
    with pytest.raises(ValueError, match=r"per_axis\*\*ndim = 4"):
        solve_sn_fixed_source(
            materials, mesh, case.quadrature, bad_width,
            scheme=LinearDiscontinuous(),
        )

    # DD/Step: no moment axis exists → a moment-resolved (4-wide) bulk is
    # rejected outright (only flat (N, ng, *spatial) is valid).
    bad_dd = np.zeros((N, ng, nx, ny, 4))
    with pytest.raises(ValueError, match=r"does not match"):
        solve_sn_fixed_source(
            materials, mesh, case.quadrature, bad_dd,
            scheme=DiamondDifference(),
        )


# ═══════════════════════════════════════════════════════════════════════
# #251 — Leg B: the BOUNDARY transverse face-slope (the OTHER half of the
# LM-1989 trap).  SEPARATE production path from Leg A (#247, the bulk lift):
# the manufactured prescribed inflow varies ALONG each face (transverse), but
# the boundary trace carries one SCALAR per face/ordinate/group and
# `_inflow_to_moments` (loss_representation.py:357-378) ZEROS the 2^{d-1}
# transverse face-slope.  Leg B closes that.
#
# THE PROBED MODE-10 EVIDENCE (the design driver, see
# `.claude/agent-memory/test-architect/issue_251_legB_boundary_gate_spec.md`):
#   * The transverse face-slope is ~12% of the face-average and decays O(h)
#     (intra-cell slope: |face_slope|/|face_bar| median 0.16/0.09/0.05 at
#     nc=8/16/32).
#   * ⚠ "IMPROVES ON FLAT" IS NOT ACHIEVABLE (sharper than Leg A): seeding the
#     REAL transverse slope makes the converged near-boundary A-error SLIGHTLY
#     WORSE (1.015e-2 → 1.030e-2 at nc=16; flipped is slightly BETTER) — the
#     boundary correction is sub-floor (below the bulk O(h²) A-error).  So the
#     teeth are STRUCTURAL ONLY (no converged-value improvement leg).
#   * The CONSUMPTION signal IS detectable: a +slope/−slope flip moves the
#     converged NEAR-BOUNDARY flux 4.1e-3 (nc=16), ≫ _CONSUMPTION_TOL (1e-8),
#     LINEAR in the slope magnitude (proves genuine consumption).
#   * The scalar-inflow no-op is BYTE-IDENTICAL (slope=0 → array_equal to today).
#
# Face-moment NORMALIZATION (the CRUX, locked — coordinate with the
# method-implementer's AngularTraceSpace widening): the trace must carry the BARE
# per-transverse Legendre coefficients [face-bar=slot 0, transverse-slope=slot
# 1].  The cochain's transverse mass `mass_1d(h_t, θ)=diag(h_t, θh_t)`
# (assemble_inflow_axis, _ubld.py:270) adds the h_t/θ weighting — the projection
# must NOT (exactly Leg A's M=diag(h,θh), transposed to the transverse axis).
# ⚠ today's scalar trace carries the cell-CENTRE value, not the average; the
# Deliverable-1 structural gate compares SLOT-1 only (slot-0 may differ
# centre-vs-average).  The face projector calls leggauss directly — NEVER
# _inflow_to_moments nor any LD cell op (L11 structural independence).
# ═══════════════════════════════════════════════════════════════════════


def _face_transverse_legendre(t_edges, j, fn_t, *, q_nodes=6):
    r"""Project a 1-D transverse field ``fn_t(t_quad) -> (q,)`` on face cell ``j``
    (``[t_edges[j], t_edges[j+1]]``) onto the BARE face Legendre coefficients.

    Returns ``(face_bar, face_slope)`` where ``face_bar = ⟨ψ,P₀⟩/⟨P₀,P₀⟩`` (the
    transverse cell AVERAGE) and ``face_slope = ⟨ψ,P₁⟩/⟨P₁,P₁⟩`` (the BARE
    transverse-slope coefficient, NO θ, NO h_t — the cochain's transverse mass
    ``diag(h_t, θh_t)`` adds them).  Structurally INDEPENDENT of production
    (``leggauss`` only — never ``_inflow_to_moments`` nor ``assemble_inflow_axis``
    nor any LD cell op; L11).  The 1-D-transverse factor of the tensor projector
    :func:`_project_scalar_to_tensor_legendre`."""
    from numpy.polynomial.legendre import leggauss

    xi, wq = leggauss(q_nodes)
    W2 = wq.sum()
    mean_p1sq = float((wq * xi * xi).sum() / W2)  # mean(P1²) = 1/3
    tL, tR = t_edges[j], t_edges[j + 1]
    tq = (tL + tR) / 2 + (tR - tL) / 2 * xi
    psi = fn_t(tq)
    face_bar = float((wq * psi).sum() / W2)
    face_slope = float((wq * xi * psi).sum() / (W2 * mean_p1sq))
    return face_bar, face_slope


# ───────────────────────────────────────────────────────────────────────
# DELIVERABLE 3 — the face-projection-correctness foundation sub-gate (L11).
# GREEN NOW (it tests only the test-side face projector; no production change).
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_face_transverse_legendre_projection_matches_hand_polynomial():
    r"""#251 sub-gate (L11): the transverse face projector reproduces the EXACT
    face Legendre coefficients of a KNOWN 1-D polynomial.

    ``ψ_face(t) = c0 + c1 t`` on a face cell ``[tL,tR]`` has the hand-derivable
    BARE coefficients (Legendre ``{1,ξ}``, ``ξ∈[-1,1]`` — the SAME basis as the
    cochain's transverse mass ``diag(h_t, θh_t)``):

        face_bar   = c0 + c1·tc          (slot 0 — transverse cell AVERAGE)
        face_slope = (h_t/2)·c1          (slot 1 — bare transverse P₁ coeff)

    A linear integrand is integrated EXACTLY by ``q_nodes>=2``, so the match is
    machine precision.  Non-tautology (L11): the projector uses only ``leggauss``
    + the hand-derived integral — NEVER ``_inflow_to_moments`` /
    ``assemble_inflow_axis`` / any LD cell op; the reference is hand-laid 1-D
    polynomial algebra.  The 1-D-transverse analog of
    ``test_tensor_legendre_projection_matches_hand_polynomial``.  ``-O``-safe."""
    tL, tR = 2.0, 2.6
    c0, c1 = -0.5, 2.1
    ht = tR - tL
    tc = (tL + tR) / 2

    face_bar, face_slope = _face_transverse_legendre(
        np.array([tL, tR]), 0, lambda t: c0 + c1 * t, q_nodes=4,
    )
    np.testing.assert_allclose(
        face_bar, c0 + c1 * tc, rtol=0, atol=1e-13,
        err_msg="face projector bar (transverse average) coeff wrong",
    )
    np.testing.assert_allclose(
        face_slope, (ht / 2) * c1, rtol=0, atol=1e-13,
        err_msg="face projector transverse-slope coeff wrong (normalization "
        "drift — must be the BARE (h_t/2)·c1, NOT θ/h_t-weighted)",
    )


@pytest.mark.foundation
def test_face_projection_slot0_is_transverse_cell_average():
    r"""#251 sub-gate (2nd leg): the face projector's slot-0 IS the transverse
    cell AVERAGE, not the cell-CENTRE value (the §1 normalization caveat).

    Today's scalar boundary trace carries the cell-CENTRE eval of the
    manufactured inflow; the face projector cell-averages along the transverse
    coordinate.  They differ by O(h²).  Pin that slot-0 equals an INDEPENDENT
    fine-quadrature (q_nodes=12) transverse cell average to ~1e-8 (so the
    projector is doing genuine averaging) — and do NOT cross-check slot-0
    against ``case.prescribed_inflow`` (centre eval, which differs by O(h²)).
    ``-O``-safe."""
    case = build_2d_cartesian_ld_stress_mms_case()
    mesh = case.build_mesh(8)
    ey = mesh.edges_y
    # Manufactured inflow on face xmin (x=0, transverse y), ordinate 0, group 0.
    quad = case.quadrature
    W = float(quad.weights.sum())
    mu_x0, mu_y0 = quad.mu_x[0], quad.mu_y[0]

    def fn_t(tq):  # ψ_{0,0}(x=0, y=tq) / W
        A, _, _, B, _, _, C, _, _ = case._drivers(np.array([0.0]), tq, 0)
        return (A.reshape(-1) + mu_x0 * B.reshape(-1) + mu_y0 * C.reshape(-1)) / W

    for j in range(len(mesh.centers_y)):
        bar_coarse, _ = _face_transverse_legendre(ey, j, fn_t, q_nodes=6)
        bar_fine, _ = _face_transverse_legendre(ey, j, fn_t, q_nodes=12)
        np.testing.assert_allclose(
            bar_coarse, bar_fine, rtol=0, atol=1e-7,
            err_msg=f"face projector slot-0 (cell {j}) is not the "
            "(quadrature-converged) transverse cell average",
        )


@pytest.mark.foundation
def test_case_projector_agrees_with_test_face_projector():
    r"""#257 S9 GATE C — single-source: the PRODUCTION MMS face projector
    ``case._project_inflow_to_face_moments`` agrees with the test-side hand
    projector ``_face_transverse_legendre`` on the manufactured inflow trace at
    MACHINE PRECISION (SLOT-1), per transverse cell, per ordinate, per group.

    The two projectors are deliberately INDEPENDENT implementations of the SAME
    bare transverse Legendre projection (Cardinal Rule 2 / L11): the production
    projector lives on the case and descends from ``case._drivers`` + leggauss;
    the test projector ``_face_transverse_legendre`` is the structurally-distinct
    leggauss reference (NEVER a production op).  Both project the SAME field
    :math:`\psi_{n,g}(\text{face},t)/W` onto the bare transverse :math:`P_1`
    coefficient (NO θ/h_t), so they MUST agree to machine precision — a drift
    here is a normalization bug (a double-applied transverse mass) that GATE B's
    threading would catch downstream; this pins it at the projector level.
    ``-O``-safe."""
    case = build_2d_cartesian_ld_stress_mms_case()
    mesh = case.build_mesh(8)
    quad = case.quadrature
    W = float(quad.weights.sum())
    mu_x, mu_y = quad.mu_x, quad.mu_y
    ng = case.n_groups
    # Face xmin (x=0, transverse y) — exercise the production projector directly.
    prod = case._project_inflow_to_face_moments(
        "x", 0.0, mesh.edges_y, 2,
    )                                              # (N, ng, n_t, 2)
    ey = mesh.edges_y
    for g in range(ng):
        for n in range(len(mu_x)):
            def fn_t(tq, _g=g, _n=n):  # ψ_{n,g}(x=0, y=tq)/W
                A, _, _, B, _, _, C, _, _ = case._drivers(np.array([0.0]), tq, _g)
                return (A.reshape(-1) + mu_x[_n] * B.reshape(-1)
                        + mu_y[_n] * C.reshape(-1)) / W
            for j in range(len(mesh.centers_y)):
                _, slope_ref = _face_transverse_legendre(ey, j, fn_t, q_nodes=6)
                np.testing.assert_array_equal(
                    prod[n, g, j, 1], slope_ref,
                    err_msg=f"the production case projector slot-1 disagrees with "
                    f"the leggauss hand reference (ordinate {n}, group {g}, cell "
                    f"{j}) — a transverse-slope normalization drift (the bare "
                    f"P₁ coeff must NOT be θ/h_t-weighted)",
                )


# ───────────────────────────────────────────────────────────────────────
# Shared #251 plumbing — the per-face transverse face-slope buffers (the
# REAL projected manufactured slope) + a faithful SURROGATE that widens
# `_inflow_to_moments` to seed slot-1 from them (mirrors what the production
# moment-resolved boundary trace will do).  The surrogate keys on the inflow
# face's transverse length + matches the per-octant ordinate slice against the
# centre-eval the trace carries (PROBED: face_view ≡ centre-eval, maxdiff 0).
# ───────────────────────────────────────────────────────────────────────


def _face_transverse_buffers(case, mesh):
    r"""Per-face ``(centre, slope)`` buffers, shape ``(N, ng, n_transverse)``.

    ``centre`` = the cell-CENTRE eval the trace carries today (the surrogate's
    matching key); ``slope`` = the BARE projected transverse face-slope (slot 1).
    Structurally independent (``_face_transverse_legendre`` → ``leggauss``)."""
    quad = case.quadrature
    W = float(quad.weights.sum())
    mu_x, mu_y = quad.mu_x, quad.mu_y
    ng = case.n_groups
    N = len(mu_x)
    ex, ey = mesh.edges_x, mesh.edges_y
    cx, cy = mesh.centers_x, mesh.centers_y
    Lx, Ly = case.length_x, case.length_y
    out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    specs = {
        "xmin": ("x", 0.0, ey, cy), "xmax": ("x", Lx, ey, cy),
        "ymin": ("y", 0.0, ex, cx), "ymax": ("y", Ly, ex, cx),
    }
    for face, (const_axis, const_val, te, tc) in specs.items():
        n_t = len(tc)
        centre = np.zeros((N, ng, n_t))
        slope = np.zeros((N, ng, n_t))
        for g in range(ng):
            # Centre eval (the trace value): drivers at the transverse cell centres.
            if const_axis == "x":
                A, _, _, B, _, _, C, _, _ = case._drivers(np.array([const_val]), tc, g)
            else:
                A, _, _, B, _, _, C, _, _ = case._drivers(tc, np.array([const_val]), g)
            A, B, C = A.reshape(-1), B.reshape(-1), C.reshape(-1)
            centre[:, g, :] = (
                A[None, :] + mu_x[:, None] * B[None, :] + mu_y[:, None] * C[None, :]
            ) / W
            for j in range(n_t):
                for n in range(N):
                    def fn_t(tq, _g=g, _n=n, _ax=const_axis, _cv=const_val):
                        if _ax == "x":
                            Aa, _, _, Bb, _, _, Cc, _, _ = case._drivers(
                                np.array([_cv]), tq, _g)
                        else:
                            Aa, _, _, Bb, _, _, Cc, _, _ = case._drivers(
                                tq, np.array([_cv]), _g)
                        return (Aa.reshape(-1) + mu_x[_n] * Bb.reshape(-1)
                                + mu_y[_n] * Cc.reshape(-1)) / W
                    _, sl = _face_transverse_legendre(te, j, fn_t, q_nodes=6)
                    slope[n, g, j] = sl
        out[face] = (centre, slope)
    return out


def _solve_with_boundary_slope(case, nc, *, slope_sign):
    r"""Full PUBLIC solve with a moment-resolved transverse face-slope on the
    prescribed-inflow boundary trace (the REAL #251 production path — no
    monkeypatch).

    ``slope_sign``: ``None`` → the scalar avg/centre-only prescribed inflow (no
    transverse moment supplied — the average-only baseline the verdict pins
    compare against); ``+1.0`` → moment-resolved trace with the REAL projected
    transverse slope on slot-1; ``-1.0`` → the flipped slope; ``0.0`` → a
    moment-resolved trace with slot-1 = 0 (the scalar no-op control, must be
    byte-identical to ``None``).  In ALL branches the boundary is built TEST-SIDE
    from ``_face_transverse_buffers`` (slot 0 = the trace scalar/centre, slot 1 =
    ``slope_sign × projected slope``) so the diagnostic's average-only vs
    average+slope distinction is preserved INDEPENDENT of production's
    ``case.prescribed_inflow`` (which, since #257 S9, emits the real slope on an
    LD mesh — this helper must NOT inherit that, it is the controlled toggle).
    L11: the slope is built by ``_face_transverse_buffers`` → ``leggauss`` only
    (never ``_inflow_to_moments``).  The moment-resolved boundary is materialised
    by the public
    :meth:`~orpheus.transport.source_sinks.AngularBoundarySourceSink.prescribed_inflow`,
    bundled with the manufactured bulk source, and solved by the public
    ``solve_sn_fixed_source`` — so it drives the production moment-resolved trace
    END-TO-END (the sweep consumes slot-1, the capture stores the outflow
    moments).  Returns ``(scalar_flux_values, mesh)``."""
    from orpheus.transport.source_sinks import (
        AngularSourceSink, AngularBoundarySourceSink,
    )
    from orpheus.transport.timed_full_field import TimedFullField

    mesh = case.build_mesh(nc)
    materials = case.build_materials(mesh)
    sn = SNMesh(mesh, case.quadrature, materials, scheme=LinearDiscontinuous())
    bufs = _face_transverse_buffers(case, mesh)
    face_values = {}
    for face, (centre, slope) in bufs.items():
        if slope_sign is None:
            # Average-only: a SCALAR inflow (N, ng, n_t) — the producer seeds
            # slot 0, slope stays zero.  Built test-side so it stays average-only
            # regardless of the production prescribed_inflow's S9 slope honesty.
            face_values[face] = centre
        else:
            # Moment-resolved (N, ng, n_t, 2): slot 0 = centre, slot 1 =
            # slope_sign × projected transverse face-slope.
            slot = np.zeros((*centre.shape, 2))
            slot[..., AVERAGE_MOMENT] = centre
            slot[..., 1] = slope_sign * slope
            face_values[face] = slot
    rhs = TimedFullField(
        interior=AngularSourceSink(values=case.external_source(mesh), space=sn.angular_bulk_space),
        boundary=AngularBoundarySourceSink.prescribed_inflow(sn, face_values),
    )
    result = solve_sn_fixed_source(
        materials, mesh, case.quadrature, rhs,
        max_inner=500, inner_tol=1e-12, scheme=LinearDiscontinuous(),
    )
    return result.scalar_flux.values, mesh


def _edge_cell_mask(mesh):
    r"""Boolean ``(nx, ny)`` mask of cells adjacent to any domain edge."""
    nx, ny = mesh.mat_map.shape
    mask = np.zeros((nx, ny), dtype=bool)
    mask[0, :] = mask[-1, :] = mask[:, 0] = mask[:, -1] = True
    return mask


# ───────────────────────────────────────────────────────────────────────
# DELIVERABLE 1 (PRIMARY) — the STRUCTURAL teeth: the widened boundary trace
# threads the projected transverse slope into slot-1 at machine precision.
# (was xfail-strict until the moment-resolved boundary trace landed, #251).
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_ld_2d_boundary_slope_threaded_through_inflow_to_moments():
    r"""#251 STRUCTURAL teeth — the production widening threads the projected
    transverse face-slope into the face cochain's slot-1 at MACHINE PRECISION
    (the analog of Leg A's lift pass-through; the production-change proof).

    The sharpest non-tautological proof of the actual #251 change (the transverse
    slope is no longer zeroed).  The spec §0 converged flux is sub-floor
    sensitive to the boundary slope, so this O(1) structural check is the right
    teeth — a regression that re-zeroes slot-1 is caught here.  L11: the projected
    reference is built by ``leggauss`` only (``_face_transverse_legendre``) — it
    never calls ``_inflow_to_moments`` nor ``assemble_inflow_axis``.

    The widened ``_inflow_to_moments`` (loss_representation.py) rank-discriminates
    the moment-resolved inflow and threads slot-1 through unchanged (this gate was
    xfail-strict until that landed; #251).  ``-O``-safe."""
    case = build_2d_cartesian_ld_stress_mms_case()
    mesh = case.build_mesh(8)
    materials = case.build_materials(mesh)
    sn = SNMesh(mesh, case.quadrature, materials, scheme=LinearDiscontinuous())
    rep = default_for(sn, sn.scheme, sn.angular_closure)

    # A moment-resolved boundary inflow: per face, (N_oct, ng, n_t, 2) with
    # slot 0 = the trace scalar (centre), slot 1 = the projected transverse slope.
    # POSITIVE: the widened _inflow_to_moments must thread slot-1 through.  We
    # build the moment-resolved inflow the production trace will carry and assert
    # the producer preserves slot-1.  (Targets _inflow_to_moments directly — the
    # stable single producer site — independent of the eventual AngularTraceSpace shape.)
    bufs = _face_transverse_buffers(case, mesh)
    # Take face xmin's inflow ordinates as a representative moment-resolved face.
    centre, slope = bufs["xmin"]                       # (N, ng, n_t)
    mu_x = case.quadrature.mu_x
    inflow_ord = np.where(mu_x > 0)[0]
    scalar_face = centre[inflow_ord]                   # (N_in, ng, n_t)
    moment_face = np.zeros((*scalar_face.shape, 2))
    moment_face[..., 0] = scalar_face
    moment_face[..., 1] = slope[inflow_ord]

    # The widened producer must (a) RECOGNISE the moment-resolved input (NOT
    # append a second moment axis) and (b) thread slot-1 through unchanged.
    # TODAY: _inflow_to_moments treats the (N,ng,n_t,2) moment face as a SCALAR
    # and appends ANOTHER (...,2) axis → output rank is wrong AND slot-1 zeroed,
    # so BOTH assertions below fail (the #251 gap).  When production widens it to
    # rank-discriminate (like Leg A's lift) and pass slot-1 through, BOTH xpass.
    widened = rep._inflow_to_moments((moment_face,))[0]
    if widened.shape != moment_face.shape:
        pytest.fail(
            f"_inflow_to_moments did not RECOGNISE the moment-resolved inflow: "
            f"output shape {widened.shape} != input {moment_face.shape} (it "
            f"appended a spurious moment axis — the #251 rank discriminator is "
            f"not yet implemented)"
        )
    np.testing.assert_array_equal(
        widened[..., 1], moment_face[..., 1],
        err_msg="the production _inflow_to_moments did NOT thread the transverse "
        "face-slope (slot-1) through — it re-zeroed it (the #251 gap is open)",
    )
    # slot-0 (the trace scalar / bar) is threaded unchanged too.
    np.testing.assert_array_equal(
        widened[..., AVERAGE_MOMENT], moment_face[..., AVERAGE_MOMENT],
        err_msg="the production _inflow_to_moments altered slot-0 of a "
        "moment-resolved inflow (the bar/average must pass through unchanged)",
    )

    # ── #257 S9 GATE B — the PRODUCTION producer leg (closes the Mode-11
    # producer-blindness the #251 surrogate had).  The MMS case's own
    # `prescribed_inflow` now EMITS a moment-resolved slot; its SLOT-1 (on the
    # inflow ordinates) must match the structurally-independent leggauss
    # reference at MACHINE PRECISION (the producer stamp).  The #251 gate above
    # pinned the CONSUMER (`_inflow_to_moments`) against a test-built slot; this
    # leg pins the PRODUCER (`case.prescribed_inflow`) — together they close the
    # producer→consumer chain.  L11: the reference is `_face_transverse_buffers`
    # → leggauss only, NEVER `case._project_inflow_to_face_moments`.
    bss = case.prescribed_inflow(sn)              # the PRODUCTION boundary source
    for face in ("xmin", "xmax", "ymin", "ymax"):
        view = bss.face_view(face)                # (N, ng, n_t, 2)
        _, slope_ref = bufs[face]                 # (N, ng, n_t), leggauss
        inflow_face = sn.angular_trace.inflow_indices_for_face(face)
        np.testing.assert_array_equal(
            view[inflow_face, ..., 1], slope_ref[inflow_face],
            err_msg=f"the PRODUCTION case.prescribed_inflow did NOT emit the "
            f"projected transverse face-slope (slot-1) on face {face!r} — the "
            f"producer dropped the slope (the #257 S9 production change is not "
            f"in effect, or the projector drifted from the leggauss reference)",
        )


@pytest.mark.foundation
def test_ld_2d_boundary_scalar_inflow_no_op_negative_control():
    r"""#251 NEGATIVE control — a SCALAR inflow (no transverse moment) keeps the
    transverse slope EXACTLY ZERO through ``_inflow_to_moments`` (byte-identical
    to today).

    The Leg-B asymmetry: the scalar default is correctly BLIND to the transverse
    slope.  This is the foundation invariant the bit-identity guard (Deliverable
    4) rests on — and it holds TODAY (the scalar widening seeds only slot 0), so
    this gate is GREEN now and MUST stay green after the production widening (the
    scalar path must remain byte-identical).  ``-O``-safe."""
    case = build_2d_cartesian_ld_stress_mms_case()
    mesh = case.build_mesh(8)
    materials = case.build_materials(mesh)
    sn = SNMesh(mesh, case.quadrature, materials, scheme=LinearDiscontinuous())
    rep = default_for(sn, sn.scheme, sn.angular_closure)

    bufs = _face_transverse_buffers(case, mesh)
    centre, _ = bufs["xmin"]
    mu_x = case.quadrature.mu_x
    scalar_face = centre[np.where(mu_x > 0)[0]]        # (N_in, ng, n_t) — scalar

    widened = rep._inflow_to_moments((scalar_face,))[0]
    # slot-0 ← the scalar; slot-1 EXACTLY zero (a scalar inflow has no slope).
    np.testing.assert_array_equal(
        widened[..., AVERAGE_MOMENT], scalar_face,
        err_msg="scalar inflow did not land on slot-0 unchanged",
    )
    np.testing.assert_array_equal(
        widened[..., 1], np.zeros_like(widened[..., 1]),
        err_msg="scalar inflow produced a non-zero transverse slope (the "
        "scalar-default no-op / Leg-B asymmetry broke)",
    )


# ───────────────────────────────────────────────────────────────────────
# DELIVERABLE 2 — the consumption-proof mutation control (anti-pattern #11).
# Flipping the CONSUMED transverse-face-slope sign moves the converged
# near-boundary flux ≫ _CONSUMPTION_TOL; the scalar-inflow gate stays blind.
# (was xfail-strict until the moment-resolved boundary trace landed, #251).
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.verifies("ld-cartesian-2d")
def test_ld_2d_boundary_slope_sign_mutation_reddens():
    r"""#251 Deliverable 2 (the PRIMARY sign-catcher) — flipping the CONSUMED
    transverse-face-slope sign CHANGES the converged near-boundary flux (the
    consumption proof) while the SCALAR-inflow path stays blind (the no-op).

    TEETH (the consumption proof, O(1) above the 1e-12 fixed point — NOT a
    value-band against A, which the spec §0b shows is sub-floor): a sign flip on
    the CONSUMED transverse face-slope moves the converged NEAR-BOUNDARY scalar
    flux by ≫ ``_CONSUMPTION_TOL`` (probed |Δφ|/|φ| ≈ 4.1e-3 near-bdy at nc=16,
    ~5.6 orders above tol; LINEAR in the slope magnitude → genuine consumption).
    The Mode-10 asymmetry: with a SCALAR inflow (slot-1 ≡ 0) the flip is a no-op
    (verified by the scalar no-op leg below) — the scalar path is correctly BLIND
    to the transverse slope, which is exactly why #251 widens the trace.

    PUBLIC-API-DRIVEN (the #251 production trace landed): the helper
    ``_solve_with_boundary_slope`` builds a moment-resolved
    :meth:`~orpheus.transport.source_sinks.AngularBoundarySourceSink.prescribed_inflow`
    (slot-1 = ±projected slope) and solves via the public
    ``solve_sn_fixed_source`` — so it exercises the production moment-resolved
    boundary trace END-TO-END (the producer stamp AND the cochain consumer).  This
    closed the prior surrogate's vv Mode-11 blindness (the monkeypatch pinned only
    the consumer; the public moment trace stamp is now exercised too — the
    threading gate ``test_ld_2d_boundary_slope_threaded_through_inflow_to_moments``
    is the dedicated stamp catcher).  ``-O``-safe."""
    case = build_2d_cartesian_ld_stress_mms_case()
    nc = 16
    phi_pos, mesh = _solve_with_boundary_slope(case, nc, slope_sign=+1.0)
    phi_neg, _ = _solve_with_boundary_slope(case, nc, slope_sign=-1.0)

    edge = _edge_cell_mask(mesh)
    d = (phi_pos - phi_neg)[:, edge]
    base = phi_pos[:, edge]
    change = float(np.sqrt(np.sum(d * d)) / np.sqrt(np.sum(base * base)))
    if not (change > _CONSUMPTION_TOL):
        pytest.fail(
            f"flipping the CONSUMED transverse face-slope did NOT change the "
            f"converged near-boundary flux (|Δφ|/|φ|={change:.3e} ≤ "
            f"{_CONSUMPTION_TOL:.0e}) — the transverse slope is NOT consumed "
            f"(the #251 boundary trace is still zeroing it)"
        )

    # The Mode-10 asymmetry: a SCALAR inflow (slot-1 ≡ 0) is unaffected by the
    # flip — its slope row is zero, so flipping zero changes nothing.  The scalar
    # no-op solve is byte-identical to today's avg/centre-only solve.
    phi_today, _ = _solve_with_boundary_slope(case, nc, slope_sign=None)
    phi_zero, _ = _solve_with_boundary_slope(case, nc, slope_sign=0.0)
    np.testing.assert_array_equal(
        phi_zero, phi_today,
        err_msg="the scalar no-op (slot-1 = 0) solve is NOT byte-identical to "
        "today's avg/centre-only solve — the scalar default / Leg-B asymmetry "
        "broke (the bit-identity guard would fail)",
    )


# ───────────────────────────────────────────────────────────────────────
# DELIVERABLE 6 — the negative pin: reject a wrong transverse-moment width.
# The boundary moment width is 2^{d-1} (= 2 for d=2 LD); the widened trace
# must reject a moment-resolved inflow with a different trailing width, and
# DD/Step must reject any moment-resolved inflow.
# (was xfail-strict until the moment-resolved boundary trace landed, #251).
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_ld_2d_boundary_trace_rejects_wrong_transverse_width():
    r"""#251 Deliverable 6 — the widened boundary trace STILL rejects a
    transverse-moment width ≠ ``2^{d-1} = 2`` (the face-moment contract), so the
    moment-resolved relaxation does not swallow real shape bugs.

    The widened ``_inflow_to_moments`` (loss_representation.py) ACCEPTS a
    moment-resolved inflow (rank-discriminated against the flat face rank) and
    rejects a wrong trailing width (e.g. 3 ≠ ``2^{d-1}``) with a clear ValueError
    naming ``2^(d-1)`` (``coding-elegance`` Pattern 4 — the relaxation does not
    swallow real shape bugs; this gate was xfail-strict until that landed, #251).
    ``-O``-safe (``pytest.raises``)."""
    case = build_2d_cartesian_ld_stress_mms_case()
    mesh = case.build_mesh(8)
    materials = case.build_materials(mesh)
    sn = SNMesh(mesh, case.quadrature, materials, scheme=LinearDiscontinuous())
    rep = default_for(sn, sn.scheme, sn.angular_closure)

    bufs = _face_transverse_buffers(case, mesh)
    centre, _ = bufs["xmin"]
    scalar_face = centre[np.where(case.quadrature.mu_x > 0)[0]]  # (N_in, ng, n_t)
    # A moment-resolved inflow with a WRONG trailing width (3 ≠ 2^{d-1} = 2).
    bad_width = np.zeros((*scalar_face.shape, 3))
    bad_width[..., 0] = scalar_face
    with pytest.raises(ValueError):
        rep._inflow_to_moments((bad_width,))
