r"""D5b — multi-D (2-D Cartesian) Linear-Discontinuous on the DAG wavefront.

Sub-step **D5b-S2** of #240 / #38 / #37: the bilinear (UBLD) LD cell kernel now
runs in d≥2 on the wavefront (``FullFieldWavefront`` / ``MovingFrontierWindow``),
consuming the verified d-generic primitive
(:mod:`orpheus.sn.spatial._ubld`).  This module carries the END-TO-END 2-D LD
gates that need a full solve:

* **D5b.2 (smoke)** — the 2-D LD kernel produces a CONVERGENT O(h²) solution on
  the existing heterogeneous 2-D MMS (vacuum edges → no ``BoundaryFlux`` change,
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
kernel in ``tests/sn/spatial/test_linear_discontinuous.py``
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
)
from orpheus.geometry import BC, CoordSystem, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import solve_sn_fixed_source
from orpheus.sn.geometry import SNMesh
from orpheus.sn.loss_representation import (
    FullFieldWavefront,
    MovingFrontierWindow,
    default_for,
)
from orpheus.sn.spatial import DiamondDifference, LinearDiscontinuous
from orpheus.transport.fields.boundary_flux import BoundaryFlux


def _l2_2d(err: np.ndarray, volumes: np.ndarray) -> float:
    r"""Volume-weighted :math:`L^2` norm on a 2-D mesh."""
    return float(np.sqrt(np.sum(volumes * err * err)))


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
    four edges → vacuum BCs automatic → no ``BoundaryFlux`` change, the
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
            sq += _l2_2d(phi[g] - ref, mesh.volumes) ** 2
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
    mfw = MovingFrontierWindow(sn)
    ffw = FullFieldWavefront(sn)
    if not isinstance(default_for(sn), MovingFrontierWindow):
        pytest.fail("expected the 2-D LD default rep to be MovingFrontierWindow")

    bf_win = BoundaryFlux.zeros_on(sn)              # VACUUM (zero domain inflow)
    bf_full = BoundaryFlux.zeros_on(sn)
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
