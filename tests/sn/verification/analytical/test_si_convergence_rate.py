"""L1 SI convergence-RATE verification: Gauss-Seidel recovery on
reflective coupling.

Closes the iterations-to-converge measurement gap (Wave O O.4a.2/O.4b
G-S->Jacobi rate regression was invisible — nothing measured n_inner).

Measurand: ``solve_sn_fixed_source(...).history.n_inner`` — the SI sweep
count to reach ``inner_tol``.  Target rate: the analytic SI spectral
radius ``rho_J = c`` (Jacobi, fully-lagged -B) vs the Gauss-Seidel
schedule (``rho_GS ~ c^2`` for the symmetric reflective model problem).

Reflective-coupling background
------------------------------
The Wave O O.4a.2 / O.4b BC extraction made the transport sweep BARE — it
reads ``psi.boundary.inflow`` as GIVEN for the whole sweep, and the
reflective coupling is applied EXTERNALLY via the sibling ``-B``
(``SNBoundaryOperator``, Wave O.2a ``_within_group_triple`` -> ``(L+C, S,
B)``).  This converted the reflective coupling from intra-sweep
Gauss-Seidel (the retired ``bc.apply``-inside-the-sweep read the LIVE
boundary buffer) to inter-sweep Jacobi (``B`` fully lagged): same
converged fixed point, slower SI spectral rate.  The planned recovery
interleaves the external ``-B`` reflect at octant-group granularity in
the SI driver (a G-S SCHEDULE of the same ``-B`` = forward-substitution
``(L+C-B_lower)^{-1}``), keeping the sweep bare and ``-B`` single-sourced.
The Krylov/GMRES path is unaffected (splitting-invariant).

Spec memo (FULL design rationale, regime selection, pillar gate):
``.claude/agent-memory/test-architect/si_gauss_seidel_rate_recovery_verification_spec.md``
(REFRESHED at HEAD 7d85222, post Phase 1+2).

Structural independence (vv-principles section 1)
-------------------------------------------------
The rate target terminates in the analytic SI spectral radius
``rho_J = c`` (the scattering ratio Sigma_s/Sigma_t), a closed-form
property of the cross sections — NOT another ORPHEUS solver and NOT a
procedurally-different derivation of the sweep.  The G-S vs Jacobi ratio
is the textbook Gauss-Seidel acceleration of a consistently-ordered
splitting.  Krylov is used ONLY as a rate-optimal LOWER-BOUND anchor
(necessary, not sufficient).  NO eigenvalue claim is paired with an MMS
reference (MMS does not prove eigenvalues).
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.derivations.continuous.analytical.homogeneous import (
    derive_2g_continuous,
)
from orpheus.geometry import (
    BC,
    Mesh1D,
    Mesh2D,
    Region,
    RegionMesh,
    StructuredGeometry,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.solver import solve_sn, solve_sn_fixed_source
from orpheus.transport.source_sinks import AngularSourceSink


# NOTE: NO module-level ``pytestmark = pytest.mark.l1``.  This file mixes
# L1 rate claims (the SI spectral-rate equation ``rho=c``) with ONE
# ``foundation`` software-invariant (G-3 flat-flux limit, no theory
# ``:label:``).  Per the harness precedence (explicit > class > inherited)
# and the ``test_l1_standoff_slab_cylinder.py`` precedent, each test
# carries its own level marker explicitly so the foundation row is not
# clobbered by an inherited ``l1`` (which raised a conflicting-marker
# warning).


# ── Jacobi baselines captured at HEAD 7d85222 (2026-06-05, post Phase 1+2) ──
# B-2g/4g reflective slab, nx=20, GL N=8, inner_tol=1e-8.  These are the
# TODAY (Jacobi) SI counts; the recovery must come in BELOW
# RECOVERY_FACTOR x baseline.  Deterministic (no RNG) — safe to pin.
#
# Re-run verbatim 2026-06-05 (this worktree, post Phase 1 R5 +
# Phase 2 O.2a; the L-(S+B) -> (L-S)-B FP reassociation shifted the slab
# counts by +1 vs the 266fcf5 spec, exactly as anticipated):
#   SLAB_2G  Jacobi n_inner = 655   c_max=0.9750  log(tol)/log(c)=728  ratio=0.90
#   SLAB_4G  Jacobi n_inner = 523   c_max=0.9722  log(tol)/log(c)=654  ratio=0.80
#   VACUUM_2G Jacobi n_inner = 128  (genuine vacuum mesh — see G-4 note)
#   BOX_2G   Jacobi n_inner = 697   (2-D 8x8, product n_mu=2 n_phi=4)
_N_JACOBI_SLAB_2G = 655
_N_JACOBI_SLAB_4G = 523
_N_JACOBI_BOX_2G = 697  # 2-D 8x8, product(n_mu=2, n_phi=4)
_N_JACOBI_VACUUM_2G = 128  # negative control (coupling inactive)
_RECOVERY_FACTOR = 0.75  # conservative; ~0.5 expected (rho_GS ~ c^2)
_TOL = 1e-8


def _reflective_slab(mat: str, ng_key: str, nx: int = 20, length: float = 2.0):
    """B-mixture reflective slab fixture: (mixture, materials, mesh, quad)."""
    m = get_mixture(mat, ng_key)
    mats = {0: m}
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=length),),
        bcs=(BC.reflective, BC.reflective),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    return m, mats, mesh, quad


def _vacuum_slab(mat: str, ng_key: str, nx: int = 20, length: float = 2.0):
    """B-mixture VACUUM slab fixture (G-4 negative control).

    IMPORTANT: the BC must be baked into the mesh as ``BC.vacuum`` — the
    ``boundary_condition=`` kwarg of ``solve_sn_fixed_source`` is IGNORED
    when the mesh carries explicit BC fields (per its docstring: "When the
    mesh carries explicit BC fields, those take precedence").  A reflective
    mesh + ``boundary_condition="vacuum"`` silently re-solves REFLECTIVE
    (655, not 128) — the latent bug this fixture exists to avoid.
    """
    m = get_mixture(mat, ng_key)
    mats = {0: m}
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=length),),
        bcs=(BC.vacuum, BC.vacuum),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    return m, mats, mesh, quad


def _iso_source(m, mesh, quad, mats, nx, ny=1):
    """Uniform iso scalar source -> per-ordinate (N, ng, nx, ny) array."""
    sn = SNMesh(mesh, quad, mats)
    return AngularSourceSink.from_isotropic(np.ones((m.ng, nx, ny)), sn).values


def _si_count(mats, mesh, quad, source, bc="reflective", tol=_TOL):
    """SI sweep count to reach ``tol`` (asserts convergence first)."""
    sol = solve_sn_fixed_source(
        mats, mesh, quad, source, boundary_condition=bc,
        inner_solver="source_iteration", max_inner=20000, inner_tol=tol,
    )
    assert sol.history.converged, "SI did not converge — count is meaningless"
    return sol.history.n_inner


# ─── section 4.1 BEFORE/AFTER RATE GATE (load-bearing) ───────────────

@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("si-spectral-rate")
@pytest.mark.xfail(
    strict=True,
    reason="G-S rate recovery not yet landed — SI is Jacobi (fully "
           "lagged -B).  This gate is RED today (documents the "
           "regression) and FLIPS green when the octant-group G-S "
           "schedule lands.  Remove the xfail in the recovery PR.",
)
def test_si_gs_rate_recovery_reflective_slab_2g():
    """L1 [iteration-cost] G-S recovery on the PRIMARY reflective config.

    Catches the Jacobi regression: on today's code n_inner ~655, which
    exceeds 0.75 x baseline (=491).  On the recovered G-S code n_inner
    drops to ~327-490 (rho_GS ~ c^2 ~ halves the count).  This is the
    true before/after gate.
    """
    m, mats, mesh, quad = _reflective_slab("B", "2g")
    source = _iso_source(m, mesh, quad, mats, nx=20)
    n_inner = _si_count(mats, mesh, quad, source)
    target = _RECOVERY_FACTOR * _N_JACOBI_SLAB_2G
    np.testing.assert_array_less(
        n_inner, target,
        err_msg=(
            f"SI count {n_inner} not below {target:.0f} "
            f"(={_RECOVERY_FACTOR} x Jacobi baseline {_N_JACOBI_SLAB_2G}). "
            f"G-S recovery either absent (still Jacobi) or insufficient."
        ),
    )


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("si-spectral-rate")
@pytest.mark.xfail(strict=True, reason="G-S recovery not yet landed (see slab_2g).")
def test_si_gs_rate_recovery_reflective_slab_4g():
    """L1 [iteration-cost] G-S recovery, 4-group (catches a group-coupling
    schedule bug invisible at 2g)."""
    m, mats, mesh, quad = _reflective_slab("B", "4g")
    source = _iso_source(m, mesh, quad, mats, nx=20)
    n_inner = _si_count(mats, mesh, quad, source)
    target = _RECOVERY_FACTOR * _N_JACOBI_SLAB_4G
    np.testing.assert_array_less(n_inner, target)


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("si-spectral-rate")
@pytest.mark.xfail(strict=True, reason="G-S recovery not yet landed (see slab_2g).")
def test_si_gs_rate_recovery_reflective_box_2g():
    """L1 [iteration-cost] G-S recovery on a 2-D fully-reflective box
    (the coupling acts on all 4 faces; catches a 1-D-only schedule)."""
    m = get_mixture("B", "2g")
    mats = {0: m}
    nx = ny = 8
    edges = np.linspace(0.0, 2.0, nx + 1)
    mesh = Mesh2D(
        edges_x=edges, edges_y=edges,
        mat_map=np.zeros((nx, ny), dtype=int),
        bc_xmin=BC.reflective, bc_xmax=BC.reflective,
        bc_ymin=BC.reflective, bc_ymax=BC.reflective,
    )
    quad = Quadrature.product(n_mu=2, n_phi=4)
    source = _iso_source(m, mesh, quad, mats, nx=nx, ny=ny)
    n_inner = _si_count(mats, mesh, quad, source)
    target = _RECOVERY_FACTOR * _N_JACOBI_BOX_2G
    np.testing.assert_array_less(n_inner, target)


# ─── section 2.1 ANALYTIC RATE pin (structurally-independent target) ──

@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("si-spectral-rate")
def test_si_jacobi_rate_matches_scattering_ratio():
    """L1 [iteration-cost] The TODAY (Jacobi) SI count matches the analytic
    spectral radius rho_J = c.

    This PASSES on today's code (it pins the Jacobi rate against the
    closed-form c) and documents the structurally-independent target the
    recovery improves upon.  n_Jacobi ~ log(tol)/log(c_max), within ~25%
    (finite-slab + multigroup correction).  NOT xfail — it is the anchor.
    """
    m, mats, mesh, quad = _reflective_slab("B", "2g")
    source = _iso_source(m, mesh, quad, mats, nx=20)
    n_inner = _si_count(mats, mesh, quad, source)
    c_max = float(np.max(np.asarray(m.scattering_ratio)))
    n_predicted = math.log(_TOL) / math.log(c_max)
    ratio = n_inner / n_predicted
    assert 0.6 <= ratio <= 1.2, (
        f"Jacobi SI count {n_inner} does not track analytic "
        f"log(tol)/log(c)={n_predicted:.0f} (ratio {ratio:.2f}); "
        f"the SI spectral rate is not rho_J=c={c_max:.4f} as expected."
    )


# ─── section 2.2 Krylov lower-bound bracket (necessary, not sufficient) ─

@pytest.mark.l1
@pytest.mark.slow
def test_jacobi_si_far_above_krylov_lower_bound():
    """L1 [iteration-cost] The Jacobi SI count is FAR above the rate-optimal
    Krylov count (documents the poor splitting; Krylov is the floor any
    SI recovery sits above).  PASSES today (SI ~860 vs Krylov ~302 at
    1e-10, ratio ~2.85)."""
    m, mats, mesh, quad = _reflective_slab("B", "2g")
    source = _iso_source(m, mesh, quad, mats, nx=20)
    n_si = _si_count(mats, mesh, quad, source, tol=1e-10)
    sol_k = solve_sn_fixed_source(
        mats, mesh, quad, source, boundary_condition="reflective",
        inner_solver="krylov", max_inner=20000, inner_tol=1e-10,
    )
    n_krylov = sol_k.history.n_inner
    assert n_si >= 1.5 * n_krylov, (
        f"Jacobi SI {n_si} not >> Krylov {n_krylov} — the splitting is "
        f"not as poor as expected (regime mischosen?)."
    )


# ─── section 4.2 CORRECTNESS guards (PASS on both Jacobi and G-S) ─────

@pytest.mark.l1
@pytest.mark.verifies("matrix-eigenvalue", "fission-matrix", "removal-matrix")
def test_recovery_preserves_kinf_2g():
    """G-1 L1 [eigenvalue] recovery must not change k_inf (2g reflective
    homogeneous = 1.875, closed-form analytical reference).

    The recovery changes the inner-SI RATE only; the converged eigenvalue
    MUST be identical.  k_inf is the dominant eigenvalue of A^{-1}F solved
    by pure linear algebra on the XS matrices — structurally independent
    of the SN sweep AND of the inner splitting schedule.
    """
    case = derive_2g_continuous()
    mat_id = next(iter(case.problem.materials.keys()))
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=mat_id, outer_thickness_cm=2.0),),
        bcs=(BC.reflective, BC.reflective),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=10),))
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    sol = solve_sn(
        case.problem.materials, mesh, quad,
        inner_solver="source_iteration",
        max_outer=1000, keff_tol=1e-12, flux_tol=1e-10,
        max_inner=5000, inner_tol=1e-10,
    )
    np.testing.assert_allclose(
        sol.keff, case.k_eff, rtol=1e-10,
        err_msg=f"keff={sol.keff:.10f} vs expected={case.k_eff:.10f}",
    )


@pytest.mark.l1
@pytest.mark.slow
def test_recovery_si_matches_krylov_fixed_point_2g():
    """G-2 L1 [flux-shape] recovered SI converges to the SAME fixed point
    as Krylov (structurally-independent splitting).  rel-Linf < 1e-8.

    Krylov (GMRES) and SI are different algorithms for the SAME
    ``(L+C-S-B)`` fixed point.  Any recovery bug that changes WHICH fixed
    point SI converges to is caught here (the G-S schedule must converge
    to the same phi, only faster).

    NOTE: ``Solution.compare`` is a same-mesh refactor-consistency tool
    (it requires ``self.mesh is other.mesh``); two separate solver calls
    build distinct internal meshes, so we compare ``scalar_flux.values``
    arrays directly.  Confirmed equivalent at this HEAD: rel-Linf ~4.4e-9.
    """
    m, mats, mesh, quad = _reflective_slab("B", "2g")
    source = _iso_source(m, mesh, quad, mats, nx=20)
    sol_si = solve_sn_fixed_source(
        mats, mesh, quad, source, boundary_condition="reflective",
        inner_solver="source_iteration", max_inner=20000, inner_tol=1e-10,
    )
    sol_k = solve_sn_fixed_source(
        mats, mesh, quad, source, boundary_condition="reflective",
        inner_solver="krylov", max_inner=20000, inner_tol=1e-10,
    )
    phi_si = sol_si.scalar_flux.values
    phi_k = sol_k.scalar_flux.values
    linf = float(np.abs(phi_si - phi_k).max())
    norm = float(np.abs(phi_k).max())
    rel_linf = linf / norm
    assert rel_linf < 1e-8, (
        f"SI and Krylov disagree on the fixed point: "
        f"scalar_flux rel-Linf={rel_linf:.2e}"
    )


@pytest.mark.foundation
def test_recovery_flat_flux_balance_limit():
    """G-3 foundation: uniform Q + uniform Sigma_t + fully reflective ->
    spatially flat infinite-medium balance phi_g.  Catches a G-S-schedule
    bug that breaks the flat fixed point.  Uses B (small absorption ->
    finite limit); a pure scatterer would be singular under full
    reflection."""
    m, mats, mesh, quad = _reflective_slab("B", "2g")
    source = _iso_source(m, mesh, quad, mats, nx=20)
    sol = solve_sn_fixed_source(
        mats, mesh, quad, source, boundary_condition="reflective",
        inner_solver="source_iteration", max_inner=20000, inner_tol=1e-12,
    )
    phi = sol.scalar_flux.values  # (ng, nx, ny)
    # Flat across space (reflective + uniform -> spatially flat to round-off).
    per_group_spread = phi.std(axis=(1, 2)) / np.abs(phi.mean(axis=(1, 2)))
    np.testing.assert_array_less(
        per_group_spread, 1e-8,
        err_msg=f"flux not spatially flat: per-group rel-spread={per_group_spread}",
    )


# ─── section 4.3 G-4 VACUUM negative control ─────────────────────────

@pytest.mark.l1
@pytest.mark.slow
def test_recovery_vacuum_count_unchanged():
    """G-4 [iteration-cost negative control] on VACUUM the coupling is
    inactive (B=0) -> the G-S schedule is a no-op -> n_inner UNCHANGED
    (~128, +/-2).  Proves the recovery touches ONLY reflective coupling.
    PASSES on both Jacobi and G-S (the count is identical).

    The mesh is built with EXPLICIT ``BC.vacuum`` (see ``_vacuum_slab``) —
    NOT a reflective mesh + ``boundary_condition="vacuum"`` (which is
    silently ignored and would re-solve reflective at 655)."""
    m, mats, mesh, quad = _vacuum_slab("B", "2g")
    source = _iso_source(m, mesh, quad, mats, nx=20)
    n_inner = _si_count(mats, mesh, quad, source, bc="vacuum")
    assert abs(n_inner - _N_JACOBI_VACUUM_2G) <= 2, (
        f"vacuum SI count {n_inner} drifted from baseline "
        f"{_N_JACOBI_VACUUM_2G} — the recovery wrongly altered the bare "
        f"(uncoupled) sweep."
    )


# ─── EIGENVALUE-path measurement seam (Phase 3 sub-step 1) ───────────

@pytest.mark.foundation
def test_eigenvalue_path_surfaces_total_inner_iterations():
    """Phase 3 sub-step 1 (measurement seam): the EIGENVALUE solve surfaces
    the total inner-SI iteration count via
    ``IterationHistory.total_inner_iterations``.

    Before the seam this was invisible: ``solve_sn`` routes through
    ``power_iteration``, whose ``IterationHistory`` carried ``n_inner=None``
    and only the (splitting-INVARIANT) OUTER count.  The seam accumulates
    each outer step's inner iterate count on the ``SNSolver`` and reads it
    into ``total_inner_iterations``, un-blocking the eigenvalue-path rate
    measurement.

    Measured Jacobi baselines (A-2g reflective, n=10, GL N=8, inner_tol=1e-8,
    post-seam): **SI total_inner=371, Krylov total_inner=310** (n_outer=3 for
    both — the outer count IS splitting-invariant; the inner SI count is
    where the G-S recovery shows).  These anchor the FUTURE eigenvalue rate
    gate (n_GS < 0.75 x 371), which lands with the eigenvalue-path G-S
    sub-step (deferred after the fixed-source G-S per the plan).  The seam is
    measurement-only — keff is unperturbed (guarded by
    :func:`test_recovery_preserves_kinf_2g`)."""
    case = derive_2g_continuous()
    mat_id = next(iter(case.problem.materials.keys()))
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=mat_id, outer_thickness_cm=2.0),),
        bcs=(BC.reflective, BC.reflective),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=10),))
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    sol = solve_sn(
        case.problem.materials, mesh, quad, inner_solver="source_iteration",
        max_outer=1000, max_inner=5000, inner_tol=1e-8,
    )
    total_inner = sol.history.total_inner_iterations
    assert total_inner is not None, (
        "eigenvalue IterationHistory.total_inner_iterations is None — the "
        "Phase 3 measurement seam did not populate it."
    )
    # At least one inner iterate per outer step (and the outer count, which
    # the eigenvalue path already surfaced, is positive).
    assert total_inner >= sol.history.n_outer > 0, (
        f"total_inner={total_inner} not >= n_outer={sol.history.n_outer} > 0"
    )
