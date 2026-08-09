"""Diagnostic: instrument the cylindrical apply-matvec internals.

Created by numerics-investigator on 2026-05-13.

Print per-step intermediates of transport_operator_matvec_cylindrical
when called on flat ψ.  Reveal exactly where the operator's per-ordinate
balance breaks vs the flat-flux invariant.
"""
from __future__ import annotations

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.quadrature import ProductQuadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operator import (
    build_equation_map_cylindrical,
    solution_to_angular_flux_cylindrical,
)


fuel = get_mixture("B", "1g")
R = 2.0
n_cells = 3
geom = StructuredGeometry(
    geometry="CYL",
    regions=(Region(mat_id=0, outer_thickness_cm=R),),
    bcs=(BC.reflective,),
)
mesh = Mesh1D.from_geometry(
    geom, region_meshes=(RegionMesh(n_cells=n_cells),),
)
quad = ProductQuadrature.create(n_mu=2, n_phi=2)
sn_mesh = SNMesh(mesh, quad)

ng = 1
N = quad.N
nx = n_cells
sig_t = np.full((nx, 1, ng), 2.0)
sigma_a = 0.1
Q_ext = 1.0
Sw = quad.weights.sum()
psi_flat_val = (Q_ext / sigma_a) / Sw

print(f"ψ_flat = {psi_flat_val:.10f}, Σw = {Sw:.6f}, 4π = {4*np.pi:.6f}")
print(f"face_areas A = {sn_mesh.reduced.face_areas}")
print(f"volumes V = {sn_mesh.volumes[:,0]}")
print(f"dr = {sn_mesh.dx}")
print(f"alpha_per_level = {sn_mesh.reduced.alpha_per_level}")
print(f"tau_mm_per_level = {sn_mesh.reduced.tau_mm_per_level}")
print(f"redist_dAw[0] = {sn_mesh.reduced.redist_dAw_per_level[0]}")
print()
print(f"quad.mu_x = {quad.mu_x}")
print(f"quad.weights = {quad.weights}")
print(f"level_indices = {quad.level_indices}")
print()

# Build the eq_map and apply-matvec by manually invoking the building blocks
eq_map = build_equation_map_cylindrical(nx, quad, ng)
sol_flat = np.full(eq_map.n_unknowns, psi_flat_val)
fi = solution_to_angular_flux_cylindrical(sol_flat, eq_map, quad, nx, ng)
print(f"fi shape: {fi.shape}")  # (ng, N, nx, 1)
print(f"fi flat? max-min: {fi[0].max() - fi[0].min():.2e}")  # Should be 0 (flat)
print()

# Step 1: outer_inflow_estimate
from orpheus.geometry.boundary import ReflectiveBoundary
from orpheus.sn.boundary_realizer import SNBoundaryRealizer, SNMethodSpace

spec_law = ReflectiveBoundary(axis="x", albedo=1.0)
bc_outer = SNBoundaryRealizer().realize(spec_law, SNMethodSpace.minimal(quad))

outer_face_psi = fi[:, :, -1, 0].T   # (N, ng)
outer_inflow_estimate = bc_outer.apply(outer_face_psi)
print(f"outer_face_psi (input to BC): shape {outer_face_psi.shape}")
print(f"  flat? values: {outer_face_psi.flatten()}")
print(f"outer_inflow_estimate (output of BC.apply): shape {outer_inflow_estimate.shape}")
print(f"  values: {outer_inflow_estimate.flatten()}")
print()

# Carlson per-level sweep — print phi_aux
sigma_t_gx = sig_t[:, 0, :].T
dr = sn_mesh.dx
mu_full = quad.mu_x

from orpheus.sn.sweep.psi_half_angle_seed import carlson_inward_sweep_from_source

for p, level_idx in enumerate(quad.level_indices):
    level_idx_arr = np.asarray(level_idx)
    mu_level = mu_full[level_idx_arr]
    w_level = quad.weights[level_idx_arr]
    print(f"--- Level p={p}: level_indices={level_idx_arr}, "
          f"mu_level={mu_level}, w_level={w_level} ---")
    within_idx_most_inward = int(np.argmin(mu_level))
    global_idx_most_inward = int(level_idx_arr[within_idx_most_inward])
    bc_outer_value_level = outer_inflow_estimate[global_idx_most_inward, :]
    print(f"  global_idx_most_inward = {global_idx_most_inward}, "
          f"bc_outer_value_level = {bc_outer_value_level}")

    # phi_0 = Σ_n w_n · ψ_n
    psi_level = fi[:, level_idx_arr, :, 0]  # (ng, M, nx)
    print(f"  psi_level shape = {psi_level.shape}, flat at {psi_level[0,0,0]:.6f}")
    phi_0 = (w_level[:, None] * psi_level[0]).sum(axis=0)  # (nx,)
    print(f"  phi_0 (from this level alone) = {phi_0}")
    print(f"  NOTE: full-quadrature phi_0 = Σw_n · ψ_n = {Sw*psi_flat_val:.6f}")

    Q_bar = 0.5 * sigma_t_gx * phi_0[None, :]  # but uses FULL phi_0
    # Wait — the actual code uses the FULL ψ Legendre fold, not per-level.
    # Let me re-check the apply-matvec call.

# Now run the actual operator (capturing intermediates would require code modification).
# Instead, run end-to-end and check.
print()
print("Now run the full transport_operator_matvec_cylindrical and "
      "print L·ψ_flat per-cell per-ordinate:")

from orpheus.sn.operator import transport_operator_matvec_cylindrical
lhs = transport_operator_matvec_cylindrical(
    sol_flat, eq_map, quad, sig_t, nx, ng,
    sn_mesh.reduced.face_areas, sn_mesh.volumes,
    sn_mesh.reduced.alpha_per_level, sn_mesh.reduced.redist_dAw_per_level,
    sn_mesh.reduced.tau_mm_per_level,
    sn_mesh=sn_mesh,
)
fi_lhs = solution_to_angular_flux_cylindrical(lhs, eq_map, quad, nx, ng)
print(f"L·ψ per (n, i):")
print(fi_lhs[0, :, :, 0])
print(f"Expected (Σ_t·ψ): {2.0*psi_flat_val:.6f}")
print(f"residual:")
print(fi_lhs[0, :, :, 0] - 2.0*psi_flat_val)
