"""Minimal reproducer for Phase G Step 2 SI-vs-Krylov discrepancy.

Find the smallest possible curvilinear case showing the discrepancy.
Start at the analytical streaming-equilibrium target φ = Q / Σ_a:
mixture B 1G (Σ_t=2, Σ_s=1.9, Σ_a=0.1, c=0.95), Q=1 isotropic,
sphere radius R=2, reflective outer BC.

Analytical answer: φ = Q/Σ_a = 1/0.1 = 10.0 EVERYWHERE; Pomraning
isotropy at the pole → ψ_n(r=0) = φ/Σ_w = 10/Σ_w for every ordinate.

Sweeps over (n_cells, n_ord) until we find the smallest combination
where SI gives >1% error against this analytical answer.
"""
from __future__ import annotations

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.solver import solve_sn_fixed_source


def _run(n_cells: int, n_ord: int, bc):
    fuel = get_mixture("B", "1g")
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(bc,),
    )
    mesh = Mesh1D.from_geometry(
        geom,
        region_meshes=(RegionMesh(n_cells=n_cells),),
    )
    quad = GaussLegendre1D.create(n_ordinates=n_ord)
    N = quad.N
    nx = mesh.N
    ng = 1
    Q = np.ones((N, nx, 1, ng))

    res_si = solve_sn_fixed_source(
        materials={0: fuel}, mesh=mesh, quadrature=quad,
        external_source=Q, boundary_condition=bc,
        inner_solver="source_iteration",
        max_inner=5000, inner_tol=1e-13,
    )
    res_kr = solve_sn_fixed_source(
        materials={0: fuel}, mesh=mesh, quadrature=quad,
        external_source=Q, boundary_condition=bc,
        inner_solver="krylov",
        max_inner=2000, inner_tol=1e-13,
    )
    psi_si = res_si.angular_flux[:, :, 0, 0]   # (N, nx)
    psi_kr = res_kr.angular_flux[:, :, 0, 0]
    return psi_si, psi_kr, quad


def test_2x2():
    """2-cell × 2-ordinate (GL-2) is the SMALLEST configuration that
    demonstrates the SI-vs-Krylov ~38% L0 discrepancy. Reflective BC.
    """
    psi_si, psi_kr, quad = _run(2, 2, BC.reflective)
    print(f"\n=== n_cells=2, N_ord=2 (GL-2), reflective ===")
    print(f"  ψ_SI: {psi_si}")
    print(f"  ψ_K:  {psi_kr}")
    print(f"  |Δψ|_inf: {float(np.max(np.abs(psi_kr-psi_si))):.4e}")
    print(f"  μ: {quad.mu_x}, w: {quad.weights}")
    # Analytic isotropic: ψ = 10 / Σ_w = 10 / 2 = 5
    psi_analytic = 10.0 / quad.weights.sum()
    print(f"  expected ψ_iso = 10/Σw = {psi_analytic:.4f}")

    # Pin: Krylov recovers analytical to machine precision
    assert float(np.max(np.abs(psi_kr - psi_analytic))) < 1e-9, (
        "Krylov failed analytical streaming-equilibrium gate"
    )
    # Pin: SI baseline produces O(1) discrepancy (the bug fingerprint)
    assert float(np.max(np.abs(psi_si - psi_analytic))) > 1.0, (
        "SI lost its expected O(1) L0-failure fingerprint; the bug "
        "may have been silently fixed elsewhere — re-investigate."
    )


def test_2x4():
    psi_si, psi_kr, quad = _run(2, 4, BC.reflective)
    print(f"\n=== n_cells=2, N_ord=4 ===")
    print(f"  ψ_SI: {psi_si}")
    print(f"  ψ_K:  {psi_kr}")
    print(f"  |Δψ|_inf: {float(np.max(np.abs(psi_kr-psi_si))):.4e}")
    print(f"  expected ψ_iso = 10/Σw = {10.0/quad.weights.sum():.4f}")


def test_3x2():
    psi_si, psi_kr, quad = _run(3, 2, BC.reflective)
    print(f"\n=== n_cells=3, N_ord=2 ===")
    print(f"  ψ_SI: {psi_si}")
    print(f"  ψ_K:  {psi_kr}")
    print(f"  |Δψ|_inf: {float(np.max(np.abs(psi_kr-psi_si))):.4e}")


def test_3x4():
    psi_si, psi_kr, quad = _run(3, 4, BC.reflective)
    print(f"\n=== n_cells=3, N_ord=4 ===")
    print(f"  ψ_SI: {psi_si}")
    print(f"  ψ_K:  {psi_kr}")
    print(f"  |Δψ|_inf: {float(np.max(np.abs(psi_kr-psi_si))):.4e}")


if __name__ == "__main__":
    test_2x2()
    test_2x4()
    test_3x2()
    test_3x4()
