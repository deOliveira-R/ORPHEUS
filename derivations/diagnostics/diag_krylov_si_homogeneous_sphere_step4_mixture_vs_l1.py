"""Diagnostic Step 4: cross-test the L1 anchor flow with mixture A.

Created by numerics-investigator on 2026-05-28.

Step 3 finding: the L1 anchor flow with mixture A (n_cells=20, GL n_ord=8,
length=2.0) also produces wrong keff (1.49) — NOT k_inf=1.875.  But the
L1 anchor with ``derive_2g_continuous`` and n_cells=10, length=2.0 PASSES.

This step bisects: is the wrongness from
  (i)  the material (mixture A vs derive_2g_continuous), or
  (ii) the mesh / geometry parameters (n_cells=20 vs 10)?

Critical insight: The brief states SI converges to 1.875 with mixture A
under solve_sn() — and Step 1 confirmed this.  So mixture A is a valid
test material; the BUG is somewhere in the Krylov inner-solver
path's interaction with it (NOT just a "loose tolerance" issue, since
Step 1 showed tightening inner_tol does not help).
"""
from __future__ import annotations

import warnings

import numpy as np


def _run_L1_flow_with_material(
    *, materials_dict, n_cells, length, max_outer=120, keff_tol=1e-10,
):
    from orpheus.geometry import (
        BC, Mesh1D, Region, RegionMesh, StructuredGeometry,
    )
    from orpheus.numerics.iteration import KrylovAcceleration
    from orpheus.numerics.operator import ZeroOperator
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.angular_flux import AngularFlux as L1AngularFlux
    from orpheus.sn.geometry import SNMesh
    from orpheus.sn.operator import CollisionOperator, StreamingOperator
    from orpheus.sn.solver import SNSolver
    from orpheus.transport.sources import PerOrdinateSource

    mat_id = next(iter(materials_dict.keys()))
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=mat_id, outer_thickness_cm=length),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(
        geom, region_meshes=(RegionMesh(n_cells=n_cells),),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    sn_mesh = SNMesh(mesh, quad, materials_dict)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        solver = SNSolver(
            sn_mesh=sn_mesh, scattering_order=0,
            inner_solver="krylov",
            max_inner=1000, inner_tol=1e-12,
        )

    L_leaf = StreamingOperator(sn_mesh, solver.mat_xs.total_cross_section)
    C_t = CollisionOperator(sn_mesh, solver.mat_xs.total_cross_section)
    LC = L_leaf + C_t
    N = sn_mesh.quad.N
    nx, ny, ng = sn_mesh.nx, sn_mesh.ny, solver.ng
    krylov = KrylovAcceleration(
        LC, solver.scattering_op, ZeroOperator(),
        preconditioner=lambda q: q,
        tol=1e-12, max_iter=1000,
        restart=min(50, N * ng * nx * ny),
    )

    phi = solver.initial_flux_distribution()
    keff = 1.0
    psi_typed_warm: L1AngularFlux | None = None
    for n_outer in range(max_outer):
        fis = solver.compute_fission_source(phi, keff)
        q_ext_per_ord = PerOrdinateSource.from_isotropic(fis, sn_mesh)
        q_ext_typed = L1AngularFlux(q_ext_per_ord.values, sn_mesh)
        psi_typed, _residuals = krylov.solve(
            q_ext_typed, initial_guess=psi_typed_warm,
        )
        psi_typed_warm = psi_typed
        phi = psi_typed.integrate_angular().values
        keff_new = solver.compute_keff(phi)
        if abs(keff_new - keff) < keff_tol:
            return keff_new, n_outer + 1
        keff = keff_new
    return keff, max_outer


def _kinf_homog_reflective(materials_dict):
    """k_inf via the 2x2 multi-group transfer matrix on a homogeneous medium.

    Solves the eigenvalue problem A·phi = (1/k)·F·phi where A is the
    multi-group removal+downscatter matrix and F is the fission matrix.
    Single-cell, no spatial dependence.
    """
    from orpheus.data.macro_xs.mixture import Mixture
    import scipy.linalg as la

    mat = next(iter(materials_dict.values()))
    ng = mat.ng
    # Removal matrix: diag(Σ_a) − Σ_s (off-diagonals)
    SigS0 = mat.SigS[0]  # (ng, ng), [g_from, g_to]
    sig_t = mat.sig_t
    # SigS rows sum to sig_s_total for that source group
    A = np.diag(sig_t) - SigS0.T  # transport on row index = sink group
    # Fission: F[g_sink, g_source] = chi[g_sink] * nu_sig_f[g_source]
    F = np.outer(mat.chi, mat.nu_sig_f)
    # Solve generalised eigenvalue: A phi = (1/k) F phi → F phi = k A phi
    w, _ = la.eig(F, A)
    return float(np.max(np.abs(w.real)))


def test_step4_kinf_mixture_a():
    from orpheus.derivations.common.xs_library import get_mixture
    fuel = get_mixture("A", "2g")
    kinf = _kinf_homog_reflective({0: fuel})
    print(f"\n[step4] k_inf (analytical, mixture A 2g) = {kinf:.10f}")


def test_step4_l1_flow_mixture_a_failing_config():
    """L1 anchor flow with mixture A at n_cells=20, length=2.0 (failing test config)."""
    from orpheus.derivations.common.xs_library import get_mixture
    fuel = get_mixture("A", "2g")
    keff, n = _run_L1_flow_with_material(
        materials_dict={0: fuel}, n_cells=20, length=2.0,
    )
    print(f"\n[step4] mixture A, n_cells=20, length=2.0: keff={keff:.10f}  n={n}")


def test_step4_l1_flow_mixture_a_l1_anchor_config():
    """Same flow with mixture A at n_cells=10, length=2.0 (L1 anchor mesh)."""
    from orpheus.derivations.common.xs_library import get_mixture
    fuel = get_mixture("A", "2g")
    keff, n = _run_L1_flow_with_material(
        materials_dict={0: fuel}, n_cells=10, length=2.0,
    )
    print(f"\n[step4] mixture A, n_cells=10, length=2.0: keff={keff:.10f}  n={n}")


if __name__ == "__main__":
    test_step4_kinf_mixture_a()
    test_step4_l1_flow_mixture_a_failing_config()
    test_step4_l1_flow_mixture_a_l1_anchor_config()
