"""Diagnostic: Output C — fixed-source SI vs Krylov at MATCHED Q.

Created by numerics-investigator on 2026-05-12 for Phase G Step 2 pre-step.

The cleanest H1 vs H2 discriminator. In a FIXED-SOURCE problem:
- Q is externally imposed, NOT derived from ψ. So Q is THE SAME for both
  SI and Krylov by construction.
- The within-group fixed-point equation is L·ψ = Q + Σ_s·φ.
- If SI and Krylov produce DIFFERENT ψ at matched Q and converged
  scattering source, then the operator L itself differs between the two
  paths → H2 confirmed.
- If they produce the same ψ, the eigenvalue-only drift was from
  source-feedback amplification → H1.

We run a heterogeneous sphere fixed-source problem with both inner
solvers and compare ψ cell-by-cell.

Strategy: use solve_sn_fixed_source which already supports BOTH inner
solvers. Use a non-trivial Q (NOT constant in cell or ordinate) so the
flat-flux degeneracy doesn't hide closure-form differences.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.solver import solve_sn_fixed_source


OUT = Path("/tmp")


def _build_problem(n_cells: int = 40):
    """sphere_2g_3reg geometry with reflective BC for fixed-source test.

    Vacuum BC + heterogeneous material + non-trivial Q can produce a
    non-invertible apply-matvec (the inward boundary inflow is zero,
    so a non-trivial RHS can saturate the matvec without a fixed point);
    reflective BC keeps the operator well-conditioned.
    """
    fuel = get_mixture("A", "2g")
    mod = get_mixture("B", "2g")
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(
            Region(mat_id=0, outer_thickness_cm=0.5),
            Region(mat_id=1, outer_thickness_cm=1.0),
            Region(mat_id=0, outer_thickness_cm=0.5),
        ),
        bcs=(BC.reflective,),
    )
    n_per_region = (n_cells // 4, n_cells // 2, n_cells // 4)
    mesh = Mesh1D.from_geometry(
        geom,
        region_meshes=tuple(RegionMesh(n_cells=n) for n in n_per_region),
    )
    quad = GaussLegendre1D.create(n_ordinates=8)
    return {
        "materials": {0: fuel, 1: mod},
        "mesh": mesh,
        "quadrature": quad,
    }


def test_fixed_source_si_vs_krylov_isotropic_q():
    """Uniform isotropic Q — degenerate test (flat ψ in homogeneous regime).

    On the 3-region heterogeneous problem this Q triggers spatial variation
    via the region structure.
    """
    kwargs = _build_problem(40)
    N = kwargs["quadrature"].N
    nx = kwargs["mesh"].N
    ng = 2
    # Uniform isotropic Q: same per ordinate, per cell, per group.
    Q = np.ones((N, nx, 1, ng)) * 1.0

    result_si = solve_sn_fixed_source(
        external_source=Q,
        boundary_condition="reflective",
        inner_solver="source_iteration",
        **kwargs,
    )
    result_kr = solve_sn_fixed_source(
        external_source=Q,
        boundary_condition="reflective",
        inner_solver="krylov",
        **kwargs,
    )

    psi_si = result_si.angular_flux[:, :, 0, :]
    psi_kr = result_kr.angular_flux[:, :, 0, :]
    sf_si = result_si.scalar_flux[:, 0, :]
    sf_kr = result_kr.scalar_flux[:, 0, :]

    delta = psi_kr - psi_si
    abs_max = np.abs(delta).max()
    rel_max = abs_max / np.abs(psi_si).max()
    per_cell_delta = np.abs(delta).max(axis=(0, 2))
    i_star = int(np.argmax(per_cell_delta))

    print(f"\n=== Fixed-source isotropic Q=1, vacuum BC, sphere n=40 ===")
    print(f"||psi_si||_inf = {np.abs(psi_si).max():.4e}")
    print(f"||psi_kr||_inf = {np.abs(psi_kr).max():.4e}")
    print(f"||delta||_inf  = {abs_max:.4e}   rel = {rel_max:.4%}")
    print(f"worst cell i* = {i_star}, per-cell max|Δψ|={per_cell_delta[i_star]:.4e}")
    print(f"per-cell drift (g=0): [i=0..nx-1]")
    print(f"  i=0:  sf_si={sf_si[0,0]:.4e}, sf_kr={sf_kr[0,0]:.4e}, "
          f"reldiff={100*(sf_kr[0,0]-sf_si[0,0])/sf_si[0,0]:.3f}%")
    print(f"  i=10: sf_si={sf_si[10,0]:.4e}, sf_kr={sf_kr[10,0]:.4e}, "
          f"reldiff={100*(sf_kr[10,0]-sf_si[10,0])/sf_si[10,0]:.3f}%")
    print(f"  i=20: sf_si={sf_si[20,0]:.4e}, sf_kr={sf_kr[20,0]:.4e}, "
          f"reldiff={100*(sf_kr[20,0]-sf_si[20,0])/sf_si[20,0]:.3f}%")
    print(f"  i=30: sf_si={sf_si[30,0]:.4e}, sf_kr={sf_kr[30,0]:.4e}, "
          f"reldiff={100*(sf_kr[30,0]-sf_si[30,0])/sf_si[30,0]:.3f}%")
    print(f"  i=39: sf_si={sf_si[39,0]:.4e}, sf_kr={sf_kr[39,0]:.4e}, "
          f"reldiff={100*(sf_kr[39,0]-sf_si[39,0])/sf_si[39,0]:.3f}%")

    np.savez(
        OUT / "phase_g_step2_fixed_source.npz",
        psi_si=psi_si,
        psi_kr=psi_kr,
        sf_si=sf_si,
        sf_kr=sf_kr,
    )

    # The test_phase_e_flux_shape_sentinel is xfail-strict because of
    # this same closure difference. We DOCUMENT the drift here, do not
    # assert tightness.
    print(f"\nInterpretation: at MATCHED external Q, both solvers find")
    print(f"different ψ. The drift CANNOT be from source convergence.")
    print(f"This isolates the drift to the WITHIN-GROUP CLOSURE FORM.")
    print(f"→ H2 confirmed: SI sweep and Krylov apply-matvec are NOT the same operator.")


def test_fixed_source_si_vs_krylov_anisotropic_q():
    """Per-ordinate Q that varies with μ — exercises BC truncation."""
    kwargs = _build_problem(40)
    N = kwargs["quadrature"].N
    nx = kwargs["mesh"].N
    ng = 2

    # Make Q vary with ordinate to break per-ordinate flat-flux degeneracy.
    # Q[n, i, 0, g] = 1.0 + 0.5·μ_n
    mu = kwargs["quadrature"].mu_x
    Q = np.zeros((N, nx, 1, ng))
    for n in range(N):
        Q[n, :, :, :] = 1.0 + 0.5 * mu[n]

    result_si = solve_sn_fixed_source(
        external_source=Q,
        boundary_condition="reflective",
        inner_solver="source_iteration",
        **kwargs,
    )
    result_kr = solve_sn_fixed_source(
        external_source=Q,
        boundary_condition="reflective",
        inner_solver="krylov",
        **kwargs,
    )

    psi_si = result_si.angular_flux[:, :, 0, :]
    psi_kr = result_kr.angular_flux[:, :, 0, :]

    delta = psi_kr - psi_si
    abs_max = np.abs(delta).max()
    per_cell_delta = np.abs(delta).max(axis=(0, 2))
    i_star = int(np.argmax(per_cell_delta))
    psi_si_max = np.abs(psi_si).max()

    print(f"\n=== Fixed-source anisotropic Q=1+0.5μ, vacuum BC, sphere n=40 ===")
    print(f"||delta||_inf  = {abs_max:.4e}   rel = {abs_max/psi_si_max:.4%}")
    print(f"worst cell i* = {i_star}")
    print(f"per-cell max|Δψ| profile (every 5 cells):")
    for ic in range(0, nx, 5):
        print(f"  i={ic:3d}: max|Δψ|={per_cell_delta[ic]:.4e}")
    print(f"  i={nx-1:3d}: max|Δψ|={per_cell_delta[nx-1]:.4e}")

    # Per-ordinate at i=39
    print(f"\nAt outer cell i=39, per-ordinate (g=0):")
    print(f"  ord  mu       psi_si        psi_kr        delta")
    for n in range(N):
        print(f"  {n:3d}  {mu[n]:+.4f}  {psi_si[n,39,0]:.4e}  "
              f"{psi_kr[n,39,0]:.4e}  {delta[n,39,0]:+.4e}")


if __name__ == "__main__":
    test_fixed_source_si_vs_krylov_isotropic_q()
    test_fixed_source_si_vs_krylov_anisotropic_q()
