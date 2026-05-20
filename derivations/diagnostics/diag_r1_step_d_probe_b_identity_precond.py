"""R-1 Step D probe B — disable the sweep preconditioner.

Hypothesis 5 says: the sweep used as a Krylov preconditioner receives
GMRES residual VECTORS, not iterates.  These residuals have no
``rhs(1)`` history, so the curvilinear Carlson seed inside
:func:`transport_sweep` falls back to its in-iteration default.  For
slab that's fine (no Carlson seed needed).  For sphere/cylinder the
sweep is then NOT a good preconditioner — possibly making the
M⁻¹·A·δ Krylov polynomial non-identity in a way that destabilises
GMRES.

This probe forces the preconditioner to identity ``M = I`` and reruns
the 5-outer-iteration sphere k_eff trace.  If the oscillation
disappears, Hypothesis 5 is correct: the SWEEP-AS-PRECONDITIONER, not
the matvec, is the offender.

Created 2026-05-19.
"""
import sys
import time
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "tests/sn/l1_analytical"))


def _run(coord: str, preconditioner_kind: str, n_outer: int = 5, n_cells: int = 10):
    warnings.simplefilter("ignore")
    from test_kinf_homogeneous import (
        _get_continuous_case,
        _homogeneous_mesh,
        _quadrature_for,
    )
    from orpheus.sn.angular_flux import AngularFlux
    from orpheus.sn.geometry import SNMesh
    from orpheus.sn.operator import CollisionOperator, StreamingOperator
    from orpheus.sn.solver import SNSolver
    from orpheus.numerics.iteration import KrylovAcceleration
    from orpheus.numerics.operator import ZeroOperator

    case = _get_continuous_case("2eg")
    mat_id = next(iter(case.problem.materials.keys()))
    mesh = _homogeneous_mesh(
        coord=coord, n_cells=n_cells, length=2.0, mat_id=mat_id,
    )
    quad_key = "slab" if coord == "slab" else coord
    quad = _quadrature_for(quad_key)
    sn_mesh = SNMesh(mesh, quad, case.problem.materials)

    solver = SNSolver(
        sn_mesh=sn_mesh, scattering_order=0,
        max_inner=300, inner_tol=1e-12, inner_solver="krylov",
    )

    # Replicate the solver's _solve_krylov but override preconditioner.
    N = sn_mesh.quad.N
    nx = sn_mesh.nx
    ng = solver.ng
    sum_w = float(sn_mesh.quad.weights.sum())

    L_leaf = StreamingOperator(sn_mesh, solver.mat_xs.total_cross_section)
    C_t = CollisionOperator(sn_mesh, solver.mat_xs.total_cross_section)
    LC = L_leaf + C_t

    S_normalised = solver.scattering_op / sum_w

    if preconditioner_kind == "default":
        preconditioner = None  # → falls back to LC.solve (sweep)
    elif preconditioner_kind == "identity":
        preconditioner = lambda q: q
    else:
        raise ValueError(preconditioner_kind)

    krylov = KrylovAcceleration(
        LC, S_normalised, ZeroOperator(),
        preconditioner=preconditioner,
        tol=1e-12, max_iter=300,
        restart=min(50, N * ng * nx * 1),
    )

    # Outer power iteration replicating the brief.
    phi = solver.initial_flux_distribution()
    keff = 1.0
    keffs = []
    iters_per_outer = []
    psi_typed_warm = None
    for n in range(n_outer):
        fis = solver.compute_fission_source(phi, keff)
        # Build typed q_ext.
        q_per_ord = np.broadcast_to(
            (fis / sum_w)[None, :, :, :], (N, ng, nx, 1),
        ).copy()
        q_ext = AngularFlux(q_per_ord, sn_mesh)
        t0 = time.perf_counter()
        psi_typed, residuals = krylov.solve(q_ext, initial_guess=psi_typed_warm)
        wall = time.perf_counter() - t0
        psi_typed_warm = psi_typed
        phi = psi_typed.integrate_angular().values
        keff = solver.compute_keff(phi)
        keffs.append(keff)
        iters_per_outer.append(len(residuals))
        print(f"  [{coord:<6} {preconditioner_kind:<8}] n={n}: keff={keff:.6f}  "
              f"inner={wall*1000:6.1f}ms  iters={len(residuals):4d}")
    return keffs


if __name__ == "__main__":
    print("=" * 70)
    print("Probe B: preconditioner choice on Krylov outer iteration")
    print("=" * 70)
    print("\n# slab: reference k_inf = 1.875")
    print("# sphere: reference k_inf = 1.875")
    print()
    print("--- slab, default sweep precond (should converge in 1) ---")
    _run("slab", "default", n_outer=3)
    print("\n--- slab, identity precond ---")
    _run("slab", "identity", n_outer=3)
    print("\n--- sphere, default sweep precond (the failing case) ---")
    _run("sphere", "default", n_outer=5)
    print("\n--- sphere, identity precond (the test of Hypothesis 5) ---")
    _run("sphere", "identity", n_outer=5)
