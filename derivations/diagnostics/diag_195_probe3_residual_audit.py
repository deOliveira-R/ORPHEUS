"""Diagnostic: Issue #195 Probe 3 — discrete residual of the manufactured solution.

Created by numerics-investigator on 2026-06-12.

The L2 error PLATEAUS (Probe 1) at ~0.041 (sphere) / ~0.049 (cylinder) with a
pole spike that GROWS under refinement (Probe 2: err[0] 0.16→0.44 across
nx 40→640).  That is a mesh-independent / DIVERGING-at-pole error — root
cause (2), a real discretisation bug, NOT a pre-asymptotic transient.

This probe pins WHERE: build the typed within-group operator triple
``(L+C, S, B)`` exactly as ``solve_sn_fixed_source`` does, feed it the EXACT
manufactured (isotropic reference) angular flux ψ_ref,n = A(r)/W, and compute
the discrete balance residual

    r = (L + C - S - B)·ψ_ref - Q_ext   (per cell, per ordinate).

If the manufactured solution SATISFIES the discrete operator (r small,
O(h^2), distributed) the source is consistent and the solver converged to
the WRONG fixed point elsewhere.  If r is LARGE and pole-concentrated, the
discrete operator's pole closure does NOT admit the manufactured solution
→ MMS source vs operator MISMATCH at the pole.

If this catches a real bug, promote to tests/sn/verification/mms/.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import (
    build_spherical_mms_case,
    build_cylindrical_mms_case,
)
from orpheus.sn import solve_sn_fixed_source
from orpheus.sn.solver import (
    SNSolver,
    _within_group_triple,
    _build_fixed_source_rhs,
    _as_sn_mesh,
)
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField


def _build_solver_and_mesh(case, nc, inner_solver="source_iteration"):
    mesh = case.build_mesh(nc)
    Q = case.external_source(mesh)
    sn_mesh = _as_sn_mesh(
        mesh, case.quadrature, case.materials,
        "vacuum", mat_map=None,
    )
    solver = SNSolver(
        sn_mesh, inner_solver=inner_solver,
        scattering_order=0, max_inner=2000, inner_tol=1e-13,
    )
    q_ext = _build_fixed_source_rhs(Q, sn_mesh)
    return mesh, sn_mesh, solver, q_ext


def _reference_angular_flux(case, sn_mesh) -> AngularFlux:
    """ψ_ref,n(r) = A(r)/W, isotropic in ordinate (Y_0^0 = 1 convention)."""
    r = sn_mesh.spatial_centers if hasattr(sn_mesh, "spatial_centers") else None
    # robust: pull radial centres from the mesh the case built
    nx = sn_mesh.spatial_shape[0]
    centers = case.build_mesh(nx).centers
    A = case.phi_exact(centers)                 # (nx,)
    W = float(sn_mesh.quad.weights.sum())
    N = sn_mesh.quad.N
    ng = sn_mesh.ng
    space_shape = (N, ng, *sn_mesh.spatial_shape)  # principled layout
    vals = np.zeros(space_shape)
    # isotropic per ordinate; spatial axes broadcast A/W over the radial axis
    iso = (A / W)
    # build an index that places iso along the radial (first spatial) axis
    bshape = [1, 1] + list(sn_mesh.spatial_shape)
    rad = np.ones(bshape)
    rad_slice = [0, 0] + [slice(None)] + [0] * (len(sn_mesh.spatial_shape) - 1)
    # Simpler: spatial_shape is (nx,) for 1-D
    vals[:, 0, :] = iso[None, :]
    return AngularFlux.from_mesh(vals, sn_mesh)


def _residual_profile(case, nc):
    mesh, sn_mesh, solver, q_ext = _build_solver_and_mesh(case, nc)
    LC, S, B = _within_group_triple(solver)
    psi_ref = _reference_angular_flux(case, sn_mesh)
    rhs = TimedFullField.zeros(
        bulk=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh,
    )
    # Wrap ψ_ref into a TimedFullField (bulk only; boundary zero for vacuum)
    psi_tff = TimedFullField(bulk=psi_ref, boundary=rhs.boundary)
    # Discrete matvec:  A·ψ = (L+C)·ψ - S·ψ - B·ψ
    Apsi = LC.apply(psi_tff)
    Spsi = S.apply(psi_tff)
    Bpsi = B.apply(psi_tff)
    # residual bulk = (LC - S - B)·ψ - q_ext  (per ordinate, per cell)
    res_vals = (
        Apsi.bulk.values - Spsi.bulk.values - Bpsi.bulk.values
        - q_ext.bulk.values
    )
    # scalar residual: angular-integrate (Σ_n w_n r_n)
    w = sn_mesh.quad.weights
    res_scalar = np.einsum("n,ng...->g...", w, res_vals)[0]  # (nx, [ny])
    res_scalar = np.asarray(res_scalar).reshape(-1)          # flatten to (nx,)
    return mesh.centers, res_scalar, res_vals


def test_probe3_sphere_residual_pole_concentration():
    """Discrete residual of ψ_ref: pole-concentrated ⟹ operator pole-closure bug."""
    case = build_spherical_mms_case()
    for nc in [40, 80, 160]:
        rc, res, _ = _residual_profile(case, nc)
        ae = np.abs(res)
        print(f"\nSPHERE nc={nc}: max|scalar res|={ae.max():.4e} at "
              f"cell {ae.argmax()} (r={rc[ae.argmax()]:.4f}); "
              f"res[0]={res[0]:+.4e}  res[mid]={res[len(res)//2]:+.4e}  "
              f"res[-1]={res[-1]:+.4e}")
        n5 = max(1, nc // 20)
        mass = res * res
        print(f"  res^2 mass: pole-5%={mass[:n5].sum()/mass.sum()*100:.1f}%  "
              f"bulk={mass[n5:-n5].sum()/mass.sum()*100:.1f}%  "
              f"outer-5%={mass[-n5:].sum()/mass.sum()*100:.1f}%")
    # The verdict assertion: if the residual is pole-concentrated AND O(1),
    # the operator does not admit the manufactured solution at the pole.
    rc, res, _ = _residual_profile(case, 160)
    assert np.abs(res[0]) < 1e-2, (
        f"Pole-cell residual {res[0]:.4e} is O(1) — the discrete operator "
        f"does NOT admit the manufactured solution at the pole (root cause 2)."
    )


def test_probe3_cylinder_residual_pole_concentration():
    case = build_cylindrical_mms_case()
    for nc in [40, 80, 160]:
        rc, res, _ = _residual_profile(case, nc)
        ae = np.abs(res)
        print(f"\nCYLINDER nc={nc}: max|scalar res|={ae.max():.4e} at "
              f"cell {ae.argmax()} (r={rc[ae.argmax()]:.4f}); "
              f"res[0]={res[0]:+.4e}  res[-1]={res[-1]:+.4e}")
    rc, res, _ = _residual_profile(case, 160)
    assert np.abs(res[0]) < 1e-2, (
        f"Pole-cell residual {res[0]:.4e} is O(1)."
    )


if __name__ == "__main__":
    for builder, name in [
        (build_spherical_mms_case, "SPHERE"),
        (build_cylindrical_mms_case, "CYLINDER"),
    ]:
        case = builder()
        print(f"\n========== {name} ==========")
        for nc in [40, 80, 160]:
            rc, res, _ = _residual_profile(case, nc)
            ae = np.abs(res)
            n5 = max(1, nc // 20)
            mass = res * res
            print(f"nc={nc:4d} max|res|={ae.max():.4e}@cell{ae.argmax()}(r={rc[ae.argmax()]:.4f})  "
                  f"res[0]={res[0]:+.4e}  res[mid]={res[len(res)//2]:+.4e}  res[-1]={res[-1]:+.4e}  "
                  f"pole5%mass={mass[:n5].sum()/mass.sum()*100:.0f}%")
            print(f"     first 6 scalar residuals: {np.array2string(res[:6], precision=4)}")
