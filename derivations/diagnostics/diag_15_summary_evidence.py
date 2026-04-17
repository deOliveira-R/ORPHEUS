"""Diagnostic: Summary evidence — the sweep's WDD angular closure creates
a non-unique discrete system at r=0 for spherical geometry.

Created by numerics-investigator on 2026-04-16.

ROOT CAUSE
----------
The spherical sweep (_sweep_1d_spherical) uses the WDD angular closure
(Morel-Montry weights) to compute angular face fluxes ψ_{n+1/2} from
the cell-average ψ_n and the incoming face flux ψ_{n-1/2}:

    ψ_{n+1/2} = (ψ_n - (1-τ)ψ_{n-1/2}) / τ

This is one-directional: ordinate n's face flux depends on ordinate n-1.
The starting condition is ψ_{1/2} = 0 (correct since α_{1/2} = 0).

The BiCGSTAB operator uses a DIFFERENT (symmetric) interpolation:

    ψ_{n+1/2} = τ·ψ_{n+1} + (1-τ)·ψ_n

Both are consistent for flat flux. But the sweep's one-directional WDD,
combined with the zero-area face at r=0 (which eliminates spatial
coupling at the innermost cell), creates a system where the iterative
sweep converges to a NON-FLAT solution that still satisfies the discrete
balance equation.

The BiCGSTAB operator, being symmetric in the angular variable, preserves
flat flux exactly and gives the correct solution.

EVIDENCE
--------
1. BiCGSTAB gives phi = 1.0 (exact) for all cells
2. Sweep gives phi[0] = 0.65 (35% error), converging and stable
3. Both satisfy their respective balance equations
4. Error at cell 0 GROWS with refinement (ratio < 1)
5. Global conservation is exact (volume-weighted phi = 1.0)
6. Cartesian sweep gives exact flat flux (no redistribution)
7. Starting from flat initial conditions, a SINGLE sweep destroys flatness
"""

import numpy as np
import sys
sys.path.insert(0, "/workspaces/ORPHEUS")

from orpheus.geometry.mesh import Mesh1D, BC
from orpheus.geometry.coord import CoordSystem
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.geometry import SNMesh
from orpheus.sn.sweep import transport_sweep
from orpheus.sn.operator import (
    build_equation_map_spherical,
    build_transport_linear_operator_spherical,
    angular_flux_to_scalar,
    solution_to_angular_flux_spherical,
)
from scipy.sparse.linalg import bicgstab


def test_evidence_sweep_vs_bicgstab():
    """Definitive comparison showing sweep bug for fixed-source spherical."""
    for nx in [5, 10, 20]:
        R = 10.0
        N_ord = 8
        mesh = Mesh1D(
            edges=np.linspace(0.0, R, nx + 1),
            mat_ids=np.ones(nx, dtype=int),
            coord=CoordSystem.SPHERICAL,
            bc_left=BC.reflective,
            bc_right=BC.reflective,
        )
        quad = GaussLegendre1D.create(N_ord)
        sn_mesh = SNMesh(mesh, quad)
        N = quad.N
        ng = 1
        W = quad.weights.sum()

        sig_t = np.full((nx, 1, ng), 1.0)
        Q_iso = np.full((nx, 1, ng), 1.0)

        # Sweep
        psi_bc = {}
        for it in range(200):
            _, phi_sweep = transport_sweep(Q_iso, sig_t, sn_mesh, psi_bc)
            if it > 0:
                res = np.linalg.norm(phi_sweep - phi_old) / max(np.linalg.norm(phi_sweep), 1e-30)
                if res < 1e-14:
                    break
            phi_old = phi_sweep.copy()

        # BiCGSTAB
        eq_map = build_equation_map_spherical(nx, quad, ng)
        T_op = build_transport_linear_operator_spherical(
            eq_map, quad, sig_t, nx, ng,
            sn_mesh.face_areas, sn_mesh.volumes,
            sn_mesh.alpha_half, sn_mesh.redist_dAw, sn_mesh.tau_mm,
        )
        rhs = np.zeros((ng, eq_map.n_eq))
        eq_idx = 0
        for ix in range(nx):
            for n in range(N):
                if ix == nx - 1 and quad.mu_x[n] < -1e-15:
                    continue
                rhs[:, eq_idx] = 1.0 / W
                eq_idx += 1
        rhs = rhs.ravel(order='F')
        solution, info = bicgstab(T_op, rhs, rtol=1e-14, maxiter=1000)
        fi = solution_to_angular_flux_spherical(solution, eq_map, quad, nx, ng)
        phi_bicgstab = angular_flux_to_scalar(fi, quad, nx, 1, ng)

        sweep_err = np.max(np.abs(phi_sweep[:, 0, 0] - 1.0))
        bicg_err = np.max(np.abs(phi_bicgstab[:, 0, 0] - 1.0))
        V = sn_mesh.volumes[:, 0]
        sweep_cons = abs(np.sum(phi_sweep[:, 0, 0] * V) / V.sum() - 1.0)

        print(f"nx={nx:3d}: sweep_max_err={sweep_err:.4e}  "
              f"bicgstab_max_err={bicg_err:.4e}  "
              f"sweep_conservation={sweep_cons:.4e}")

        # The sweep error should be large and grow with refinement
        # BiCGSTAB should be exact
        assert bicg_err < 1e-10, f"BiCGSTAB not exact: err={bicg_err}"


if __name__ == "__main__":
    print("\n" + "="*70)
    print("SUMMARY EVIDENCE: Sweep vs BiCGSTAB for spherical fixed-source")
    print("Expected: phi = 1.0 everywhere (constant source, reflective BCs)")
    print("="*70)
    test_evidence_sweep_vs_bicgstab()
    print("="*70)
