"""Diagnostic: Output A — per-cell-per-ordinate ψ_si vs ψ_kr.

Created by numerics-investigator on 2026-05-12 for Phase G Step 2 pre-step.

Question: where do SI fixed-point and Krylov fixed-point disagree at the
per-cell-per-ordinate level, given they share the Phase F Carlson seed?

Discriminates H1 (BC timing / source convergence) from H2 (closure
algebra) by attributing the drift to cells.

Saves:
  /tmp/phase_g_step2_psi_si.npz — angular_flux from SI path
  /tmp/phase_g_step2_psi_kr.npz — angular_flux from Krylov path
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.solver import solve_sn


OUT = Path("/tmp")


def _build_sphere_2g_3reg(n_cells: int = 40):
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
    return dict(
        materials={0: fuel, 1: mod},
        mesh=mesh,
        quadrature=GaussLegendre1D.create(n_ordinates=8),
        scattering_order=0,
    )


def test_output_a_per_cell_per_ordinate_psi_comparison():
    kwargs = _build_sphere_2g_3reg(40)
    result_si = solve_sn(inner_solver="source_iteration", **kwargs)
    result_kr = solve_sn(inner_solver="krylov", **kwargs)

    # angular_flux: (N_ord, nx, ny=1, ng) → squeeze ny.
    psi_si = result_si.angular_flux[:, :, 0, :]   # (N, nx, ng)
    psi_kr = result_kr.angular_flux[:, :, 0, :]
    sf_si = result_si.scalar_flux[:, 0, :]        # (nx, ng)
    sf_kr = result_kr.scalar_flux[:, 0, :]

    # Eigenvector normalisation differs between the two power iterations
    # (the matvec form has a different fixed point because of the residual
    # WDD asymmetry). Normalise by sum-scalar-flux so per-cell shapes
    # are comparable.
    norm_si = float(sf_si.sum())
    norm_kr = float(sf_kr.sum())
    psi_si_n = psi_si / norm_si
    psi_kr_n = psi_kr / norm_kr
    sf_si_n = sf_si / norm_si
    sf_kr_n = sf_kr / norm_kr

    delta_psi = psi_kr_n - psi_si_n  # (N, nx, ng)

    # Per-cell summary: max abs delta across ordinate × group.
    delta_per_cell = np.abs(delta_psi).max(axis=(0, 2))  # (nx,)
    # L2 norm of per-cell ψ drift (combined over ordinate × group).
    delta_l2_per_cell = np.sqrt((delta_psi ** 2).sum(axis=(0, 2)))

    i_star = int(np.argmax(delta_per_cell))
    psi_si_max = float(np.abs(psi_si_n).max())
    rel_drift_inf = delta_per_cell.max() / psi_si_max

    print(f"\nN_ord={psi_si.shape[0]}, nx={psi_si.shape[1]}, ng={psi_si.shape[2]}")
    print(f"keff SI={result_si.keff:.8f}   Kr={result_kr.keff:.8f}   "
          f"diff={100*(result_kr.keff-result_si.keff)/result_si.keff:.3f}%")
    print(f"\nWorst cell i*={i_star}")
    print(f"  max |Δψ| at i*  = {delta_per_cell[i_star]:.6e}")
    print(f"  ||Δψ||_2 at i*  = {delta_l2_per_cell[i_star]:.6e}")
    print(f"  max |ψ_si_n|    = {psi_si_max:.6e}")
    print(f"  rel drift ∞    = {rel_drift_inf:.4%}")

    # Per-ordinate at i*: which ordinate dominates?
    quad = result_si.quadrature
    print(f"\nAt cell i*={i_star}: per-ordinate Δψ contribution (g=0)")
    print(f"  {'n':>3s}  {'mu_n':>8s}  {'psi_si':>12s}  {'psi_kr':>12s}  "
          f"{'delta_psi':>12s}")
    for n in range(psi_si.shape[0]):
        print(
            f"  {n:>3d}  {quad.mu_x[n]:+8.4f}  "
            f"{psi_si_n[n, i_star, 0]:+12.6e}  "
            f"{psi_kr_n[n, i_star, 0]:+12.6e}  "
            f"{delta_psi[n, i_star, 0]:+12.6e}"
        )

    # Drift vs cell index for both groups.
    print(f"\nMax |Δψ| per cell (top 10 by magnitude):")
    sort_cells = np.argsort(delta_per_cell)[::-1][:10]
    for ic in sort_cells:
        print(f"  i={ic:3d}  max|Δψ|={delta_per_cell[ic]:.4e}  "
              f"||Δψ||_2={delta_l2_per_cell[ic]:.4e}  "
              f"sf_si={sf_si_n[ic, 0]:.4e}  sf_kr={sf_kr_n[ic, 0]:.4e}")

    np.savez(
        OUT / "phase_g_step2_psi.npz",
        psi_si=psi_si,
        psi_kr=psi_kr,
        sf_si=sf_si,
        sf_kr=sf_kr,
        mu=quad.mu_x,
        weights=quad.weights,
        keff_si=result_si.keff,
        keff_kr=result_kr.keff,
        i_star=i_star,
        norm_si=norm_si,
        norm_kr=norm_kr,
    )
    print(f"\nSaved to {OUT/'phase_g_step2_psi.npz'}")

    # Sanity assertion — the documented Phase F drift on this case.
    # rel drift ~ 6-10% at the worst cell after sf-normalisation
    # (the keff differs by 0.286%, but eigenvector shapes differ more).
    assert delta_per_cell[i_star] > 0, "No drift detected (Step 1 broke?)"


if __name__ == "__main__":
    test_output_a_per_cell_per_ordinate_psi_comparison()
