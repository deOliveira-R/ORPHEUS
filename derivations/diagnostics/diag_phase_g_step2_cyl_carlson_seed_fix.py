"""Diagnostic: cylinder Carlson seed Q_bar normalisation fix.

Created by numerics-investigator on 2026-05-13.

HYPOTHESIS: the MorelMontryAngularSweep Carlson seed uses
``Q_bar = 0.5 · Σ_t · phi_0_level`` which is correct for sphere
(where Σw = 2, so 0.5 = 1/Σw) but wrong for cylinder where the
per-level azimuthal weights sum to 2π — the correct convention is
``Q_bar = Σ_t · phi_0_level / Σ w_level``.

This diagnostic patches the convention IN-MEMORY (via a custom closure
class) and re-runs the L·ψ residual check on flat ψ.
"""
from __future__ import annotations

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.quadrature import ProductQuadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import (
    build_equation_map_cylindrical,
    transport_operator_matvec_cylindrical,
)
from orpheus.sn.spatial.pole_angular_closure import (
    MorelMontryAngularSweep,
    _mm_weighted_angular_recurrence_single_level,
)
from orpheus.sn.spatial.psi_half_angle_seed import (
    CarlsonSweepContext, carlson_inward_sweep_from_source,
)


class FixedCarlsonSeed:
    """Carlson seed with Q_bar normalisation fixed for general Σw_level."""

    is_linear = True

    def __call__(self, psi_level, context):
        sigma_t = context.sigma_t
        dr = context.dr
        weights = context.weights
        bc_outer = context.bc_outer_value
        Sw_level = weights.sum()

        phi_0 = np.einsum("m,gmx->gx", weights, psi_level)
        # FIX: Q_bar = Σ_t · phi_0 / Sw_level (not 0.5 · Σ_t · phi_0).
        Q_bar = sigma_t * phi_0 / Sw_level

        return carlson_inward_sweep_from_source(
            Q_bar=Q_bar, sigma_t=sigma_t, dr=dr, bc_outer_value=bc_outer,
        )


from dataclasses import dataclass
from orpheus.sn.spatial.pole_angular_closure import PoleAngularClosureBase


@dataclass(frozen=True, slots=True)
class FixedCylinderMM(
    PoleAngularClosureBase, key="fixed_cylinder_mm",
):
    """MorelMontryAngularSweep with FixedCarlsonSeed."""

    is_linear = True

    def __call__(
        self, psi_cells, alpha_half, redist_dAw, tau_mm, volume,
        level_indices=None, carlson_context=None,
    ):
        if level_indices is None:
            psi_half_seed_arr = None
            if carlson_context is not None:
                psi_half_seed_arr = FixedCarlsonSeed()(
                    psi_cells, carlson_context,
                )
            return _mm_weighted_angular_recurrence_single_level(
                psi_cells, alpha_half, redist_dAw, tau_mm, volume,
                psi_half_seed=psi_half_seed_arr,
            )
        # Cylindrical
        ng, N, nx = psi_cells.shape
        redist = np.zeros((ng, N, nx), dtype=psi_cells.dtype)
        for p, level_idx in enumerate(level_indices):
            psi_level = psi_cells[:, level_idx, :]
            psi_half_seed_arr = None
            if carlson_context is not None:
                psi_half_seed_arr = FixedCarlsonSeed()(
                    psi_level, carlson_context[p],
                )
            redist_level = _mm_weighted_angular_recurrence_single_level(
                psi_level, alpha_half[p], redist_dAw[p], tau_mm[p],
                volume, psi_half_seed=psi_half_seed_arr,
            )
            redist[:, level_idx, :] = redist_level
        return redist


fuel = get_mixture("B", "1g")
R = 2.0
geom = StructuredGeometry(
    geometry="CYL",
    regions=(Region(mat_id=0, outer_thickness_cm=R),),
    bcs=(BC.reflective,),
)


def test_fixed_seed_recovers_flat():
    print("\n=== Cylinder apply-matvec L·ψ_flat residual with patched Carlson seed ===")
    for n_mu in [2, 4]:
        for n_phi in [2, 4]:
            print(f"\n  ProductQuadrature(n_mu={n_mu}, n_phi={n_phi})")
            for n_cells in [3, 5, 10, 20, 40]:
                mesh = Mesh1D.from_geometry(
                    geom, region_meshes=(RegionMesh(n_cells=n_cells),),
                )
                quad = ProductQuadrature.create(n_mu=n_mu, n_phi=n_phi)
                sn_mesh = SNMesh(
                    mesh, quad, pole_angular_closure=FixedCylinderMM(),
                )
                ng = 1
                sig_t = np.full((n_cells, 1, ng), 2.0)
                psi_flat = 10.0 / quad.weights.sum()
                eq_map = build_equation_map_cylindrical(n_cells, quad, ng)
                sol_flat = np.full(eq_map.n_unknowns, psi_flat)

                lhs_fixed = transport_operator_matvec_cylindrical(
                    sol_flat, eq_map, quad, sig_t, n_cells, ng,
                    sn_mesh.reduced.face_areas, sn_mesh.volumes,
                    sn_mesh.reduced.alpha_per_level,
                    sn_mesh.reduced.redist_dAw_per_level,
                    sn_mesh.reduced.tau_mm_per_level,
                    sn_mesh=sn_mesh,
                    pole_angular_closure=FixedCylinderMM(),
                )
                res = np.max(np.abs(lhs_fixed - 2.0 * psi_flat))
                print(f"    n_cells={n_cells:3d}, Σw={quad.weights.sum():.4f}: "
                      f"||L·ψ − Σ_t·ψ||_inf = {res:.3e}")


if __name__ == "__main__":
    test_fixed_seed_recovers_flat()
