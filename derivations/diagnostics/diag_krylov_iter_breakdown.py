"""Diagnostic: Krylov per-iteration breakdown for cylinder heavy-scattering.

Created by numerics-investigator on 2026-05-14.
Pin for the Krylov inner-solver structural-improvement memo
(.claude/agent-memory/numerics-investigator/krylov_inner_solver_profile_2026_05_14.md).

The cProfile run of the slow case (cylinder Krylov nx=20, n_mu=4,
n_phi=4) attributes 51 % of cumulative time to
``transport_operator_matvec_cylindrical`` (the GMRES matvec) and 47 %
to ``_run_1d_sweep`` (the sweep wrapped as left preconditioner).  This
diagnostic isolates the GMRES iteration count and per-matvec wall
clock so the memo can:

    1. Quote the GMRES iteration count to converge at the production
       inner_tol = 1e-8.
    2. Quote mean per-iteration wall clock.
    3. Partition the per-iteration time into matvec (apply path) vs
       preconditioner (sweep path) vs orthogonalisation residual.

This is a one-shot profiling probe, NOT a regression test.  Do NOT
promote — the wall-clock thresholds are machine-dependent.  If the
Krylov structural improvement lands later, a follow-up regression
test belongs in ``tests/sn/`` and should pin GMRES iteration count
(machine-independent) and a relative speedup ratio (also machine-
independent).
"""
from __future__ import annotations

import time
import numpy as np
from scipy.sparse.linalg import LinearOperator as ScipyLinearOperator
from scipy.sparse.linalg import gmres

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.quadrature import ProductQuadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.solver import (
    SNSolver,
    _build_rhs_spherical,
)
from orpheus.sn.operator import (
    solution_to_angular_flux_cylindrical,
)


def build_solver():
    """Reproduce the slow test case's SNSolver setup."""
    fuel = get_mixture("B", "1g")
    geom = StructuredGeometry(
        geometry="CYL",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=20),))
    quad = ProductQuadrature.create(n_mu=4, n_phi=4)
    sn_mesh = SNMesh(mesh, quad)
    solver = SNSolver(
        materials={0: fuel},
        sn_mesh=sn_mesh,
        inner_solver="krylov",
        scattering_order=0,
        keff_tol=1e-7,
        flux_tol=1e-6,
        max_inner=200,
        inner_tol=1e-8,
    )
    return solver, quad, mesh


def run_breakdown():
    solver, quad, mesh = build_solver()
    nx, ny, ng = solver.sn_mesh.nx, solver.sn_mesh.ny, solver.ng

    # Drive Krylov manually so we can attach a callback that counts
    # iterations and measures wall clock per matvec.
    sum_w = float(quad.weights.sum())
    Q = np.ones((quad.N, nx, ny, ng))   # uniform isotropic source
    # The fixed-source code path scales Q internally; do the same.
    Q_iso_scalar = Q.mean(axis=0)  # (nx, ny, ng) — already uniform
    fission_source = np.zeros((nx, ny, ng))
    flux_distribution = np.zeros((nx, ny, ng))

    eq_map = solver.L._ensure_eq_map()
    n_unknowns = eq_map.n_unknowns

    rhs = _build_rhs_spherical(  # cylindrical alias
        fission_source, flux_distribution, eq_map, quad,
        solver.sig_s, solver.sig2, solver.sn_mesh.mat_map,
        nx, ng,
    )
    # Add the external isotropic source (mirroring _solve_fixed_source_krylov).
    # The RHS of L psi = b has Q/W folded in at the same place fission lives.
    Q_iso_packed = np.zeros_like(rhs)
    for k in range(eq_map.n_eq):
        Q_iso_packed_idx = slice(k * ng, (k + 1) * ng)  # Fortran-order packing
    # For this microbench we use rhs as-is (zero everywhere); the
    # streaming-equilibrium test only exercises the homogeneous scatter+source
    # convergence.  We instead build the RHS by injecting Q via the manufactured
    # source-iteration RHS.
    Q_external_packed = np.full(n_unknowns, 1.0 / sum_w)  # uniform Q per unknown
    rhs = rhs + Q_external_packed

    # Hot-path matvec + preconditioner.
    matvec_count = [0]
    matvec_wall = []

    def timed_matvec(x):
        t0 = time.perf_counter()
        y = solver.L.apply(x)
        matvec_wall.append(time.perf_counter() - t0)
        matvec_count[0] += 1
        return y

    L_scipy = ScipyLinearOperator(
        (n_unknowns, n_unknowns), matvec=timed_matvec, dtype=float,
    )
    precond = solver._make_sweep_preconditioner(eq_map, n_unknowns, sum_w)
    precond_wall = []
    orig_precond_matvec = precond.matvec

    def timed_precond_matvec(x):
        t0 = time.perf_counter()
        y = orig_precond_matvec(x)
        precond_wall.append(time.perf_counter() - t0)
        return y

    precond_timed = ScipyLinearOperator(
        (n_unknowns, n_unknowns), matvec=timed_precond_matvec, dtype=float,
    )

    iter_count = [0]
    def callback(xk):
        iter_count[0] += 1

    x0 = np.ones(n_unknowns)
    t_total = time.perf_counter()
    try:
        sol, info = gmres(
            L_scipy, rhs, x0=x0, M=precond_timed,
            rtol=1e-8, atol=0.0, maxiter=200,
            restart=min(50, n_unknowns), callback=callback,
        )
    except TypeError:
        sol, info = gmres(
            L_scipy, rhs, x0=x0, M=precond_timed,
            tol=1e-8, maxiter=200,
            restart=min(50, n_unknowns), callback=callback,
        )
    wall_total = time.perf_counter() - t_total

    return {
        "n_unknowns": n_unknowns,
        "info": info,
        "n_gmres_iterations": iter_count[0],
        "n_matvecs": matvec_count[0],
        "n_precond_matvecs": len(precond_wall),
        "wall_total": wall_total,
        "matvec_wall_mean_ms": 1000.0 * float(np.mean(matvec_wall)),
        "matvec_wall_std_ms": 1000.0 * float(np.std(matvec_wall)),
        "matvec_wall_total_s": float(np.sum(matvec_wall)),
        "precond_wall_mean_ms": 1000.0 * float(np.mean(precond_wall)),
        "precond_wall_std_ms": 1000.0 * float(np.std(precond_wall)),
        "precond_wall_total_s": float(np.sum(precond_wall)),
    }


if __name__ == "__main__":
    out = run_breakdown()
    for k, v in out.items():
        if isinstance(v, float):
            print(f"  {k:30s} = {v:.4f}")
        else:
            print(f"  {k:30s} = {v}")
