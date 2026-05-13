"""Diagnostic: homogeneous-medium streaming-equilibrium fixed-source.

Created by numerics-investigator on 2026-05-12 for Phase G Step 2 pre-step.

Homogeneous medium + uniform isotropic Q + reflective BC → analytical
streaming equilibrium: φ = Q / Σ_a (since Σ_a = Σ_t − Σ_s for pure
absorption + scattering balance, and reflective BC removes leakage).

For 1G pure-scatterer (Σ_s = Σ_t): infinite multiplication, RHS = LHS
identity. For 1G with c < 1: φ = Q / Σ_a everywhere.

This is the canonical L0-SN-001 streaming-equilibrium test.  Pomraning
1989 isotropy at the pole IS exact here: ψ_n(r=0) = Q / (Σ_t·Σw) for
every ordinate.

Run both inner solvers and compare against the analytical limit.
"""
from __future__ import annotations

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.solver import solve_sn_fixed_source


def test_homogeneous_streaming_equilibrium():
    """Pure-absorber homogeneous sphere with uniform Q.

    Pure absorber (Σ_s = 0) at Σ_t = 0.5 gives φ = Q / Σ_t = 1/0.5 = 2.0
    everywhere (reflective BC removes leakage; with no scattering,
    no within-group iteration is needed).
    """
    # Use a 1G mixture with no scattering.
    # Or use the 2G mixture but the homogeneous result is by group.
    # Mixture B: SigT=2, SigS=1.9 (so SigA=0.1), no fission. c = 0.95.
    fuel = get_mixture("B", "1g")
    print(f"Mixture: SigT={fuel.SigT}, SigS_0={fuel.SigS[0].toarray() if hasattr(fuel.SigS[0], 'toarray') else fuel.SigS[0]}")
    print(f"  SigP (nu*SigF)={fuel.SigP}, chi={fuel.chi}")

    n_cells = 40
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(
        geom,
        region_meshes=(RegionMesh(n_cells=n_cells),),
    )
    quad = GaussLegendre1D.create(n_ordinates=8)
    N = quad.N
    nx = mesh.N
    ng = 1

    Q = np.ones((N, nx, 1, ng))

    result_si = solve_sn_fixed_source(
        materials={0: fuel},
        mesh=mesh,
        quadrature=quad,
        external_source=Q,
        boundary_condition="reflective",
        inner_solver="source_iteration",
    )
    result_kr = solve_sn_fixed_source(
        materials={0: fuel},
        mesh=mesh,
        quadrature=quad,
        external_source=Q,
        boundary_condition="reflective",
        inner_solver="krylov",
    )

    psi_si = result_si.angular_flux[:, :, 0, :]
    psi_kr = result_kr.angular_flux[:, :, 0, :]
    sf_si = result_si.scalar_flux[:, 0, :]
    sf_kr = result_kr.scalar_flux[:, 0, :]

    # Analytical: for fully-scattering medium with reflective BC and
    # uniform Q, infinite cancellation — the system has no unique
    # solution. For non-pure-scattering (some absorption), φ = Q / Σ_a
    # where the scattering source iteration converges geometrically.
    #
    # For mixture A 1G: SigT, SigS_self, SigF
    SigT = float(fuel.SigT[0])
    SigS_self = float(fuel.SigS[0].toarray()[0, 0])
    NuSigF = float(fuel.SigP[0])  # nu*Sigma_f
    # For k-eigenvalue, the steady state has no external source.
    # For our fixed-source: φ steady = Q / (SigT - SigS_self - chi*NuSigF/k)
    # but we're not doing eigenvalue — pure source iteration converges
    # geometrically at c = (SigS + NuSigF) / SigT < 1 only.
    c_eff = (SigS_self + NuSigF) / SigT
    if c_eff < 1:
        phi_analytic = 1.0 / (SigT - SigS_self - NuSigF)
        print(f"  c_eff = {c_eff:.3f} (subcritical), phi_analytic = {phi_analytic:.4f}")
    else:
        print(f"  c_eff = {c_eff:.3f} (supercritical) — fixed source diverges, skip")
        return

    print(f"\nFixed-source homogeneous sphere n={n_cells}, Q=1")
    print(f"  sf_si analytic-rel-error:")
    rel_si = np.abs(sf_si[:, 0] - phi_analytic) / phi_analytic
    print(f"    i=0: {rel_si[0]:.6f}, i=20: {rel_si[20]:.6f}, i={nx-1}: {rel_si[nx-1]:.6f}")
    rel_kr = np.abs(sf_kr[:, 0] - phi_analytic) / phi_analytic
    print(f"  sf_kr analytic-rel-error:")
    print(f"    i=0: {rel_kr[0]:.6f}, i=20: {rel_kr[20]:.6f}, i={nx-1}: {rel_kr[nx-1]:.6f}")

    print(f"\nPer-ordinate ψ at pole (i=0): should be phi_analytic / Σw =")
    psi_analytic = phi_analytic / quad.weights.sum()
    print(f"  expected: {psi_analytic:.4f}")
    print(f"  ord  mu       psi_si        psi_kr")
    for n in range(N):
        print(f"  {n:3d}  {quad.mu_x[n]:+.4f}  {psi_si[n,0,0]:.4e}  {psi_kr[n,0,0]:.4e}")


if __name__ == "__main__":
    test_homogeneous_streaming_equilibrium()
