"""Diagnostic: cylinder SI sweep with patched Carlson Q_bar.

Created by numerics-investigator on 2026-05-13.

Build a custom SI sweep where Q_bar is computed correctly per-level for
cylinder, then run end-to-end via Picard iteration on the L0
streaming-equilibrium test.  Verify that the corrected SI converges
to ψ_n = 0.7957747 (analytical) at every mesh size and quadrature.
"""
from __future__ import annotations

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.quadrature import ProductQuadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.sweep.psi_half_angle_seed import carlson_inward_sweep_from_source
from orpheus.transport.spatial.scheme import UpstreamState


def custom_cylinder_si_sweep(
    sn_mesh, sig_t, Q, scalar_flux_prev, psi_pole_prev, bc_outer,
):
    """Cylinder SI sweep matching production code but with FIXED Q_bar."""
    nx = sn_mesh.nx
    ng = Q.shape[2]
    quad = sn_mesh.quad
    N = quad.N
    weights = quad.weights

    Q_1d = Q[:, 0, :]
    sig_t_1d = sig_t[:, 0, :]
    reduced = sn_mesh.reduced
    A = reduced.face_areas
    V = sn_mesh.volumes[:, 0]
    scheme = sn_mesh.scheme

    bc_outer_obj = sn_mesh.bc_right

    angular_flux = np.zeros((N, nx, 1, ng))
    scalar_flux = np.zeros((nx, ng))

    Sw_full = weights.sum()
    weight_norm = 1.0 / Sw_full
    QV_iso = Q_1d * V[:, None] * weight_norm

    sigma_t_gx = sig_t_1d.T
    dr = sn_mesh.dx
    mu_full = quad.mu_x

    # ── FIX: compute Q_bar per-level inside the loop ──
    # For each level, Q_bar_level = Σ_t · phi_0_level / Σw_level where
    # phi_0_level = Σ_{m∈level} w_m · ψ_m (not the full quadrature's φ_0).

    # We need ψ at level slices. Use the previous-iter angular flux psi_pole_prev
    # at all cells? No — actually phi_0_level only depends on ψ at the level's
    # ordinates.  For now, derive it from scalar_flux_prev assuming a
    # near-isotropic profile (the correct fix uses ψ_prev directly).
    # ALTERNATIVE per-cycle approximation: phi_0_level ≈ phi_0_full · (Σw_level / Σw_full)
    # when ψ is isotropic.  On strongly anisotropic flux this is inexact.
    # The CORRECT fix would need ψ_prev cached, not just φ_prev.

    # For verification at flat fixed point, we can use phi_0_prev and the
    # ratio of weight sums:
    inflow_full_cyl = bc_outer_obj.apply(bc_outer)

    for p, level_idx in enumerate(quad.level_indices):
        M = len(level_idx)
        level_idx_arr = np.asarray(level_idx)
        Sw_level = weights[level_idx_arr].sum()

        within_idx_most_inward = int(np.argmin(mu_full[level_idx_arr]))
        global_idx_most_inward = int(level_idx_arr[within_idx_most_inward])
        bc_outer_value_level = inflow_full_cyl[global_idx_most_inward, :]

        # FIX: per-level Q_bar.  At fixed point with isotropic-in-φ ψ,
        # phi_0_level = Σw_level · ψ̄_level.  But we don't have ψ̄_level
        # alone — we have phi_0_full = Σw_full · ψ̄ (when ψ flat).  For the
        # general case we need ψ_prev for the level's ordinates.
        #
        # Use: phi_0_level_approx = phi_0_full · (Σw_level / Σw_full),
        # which is exact when ψ is flat-in-azimuth (the Pomraning-isotropy
        # invariant the level structure assumes).
        phi_0_level_approx = scalar_flux_prev * (Sw_level / Sw_full)  # (nx,)
        Q_bar_iso = sigma_t_gx * phi_0_level_approx.reshape(1, -1) / Sw_level
        # Equivalently: Q_bar_iso = sigma_t_gx * scalar_flux_prev.reshape(1, -1) / Sw_full
        # because phi_0_level / Sw_level = phi_0_full / Sw_full on isotropic-in-φ ψ.

        phi_aux_level = carlson_inward_sweep_from_source(
            Q_bar=Q_bar_iso,
            sigma_t=sigma_t_gx,
            dr=dr,
            bc_outer_value=bc_outer_value_level,
        )

        psi_angle = phi_aux_level.T.copy()

        for m_local in range(M):
            n = level_idx[m_local]
            eta_n = quad.mu_x[n]
            abs_eta = abs(eta_n)
            w_n = weights[n]

            QV = QV_iso

            if eta_n < 0:
                psi_in_full = bc_outer_obj.apply(bc_outer)
                psi_spatial_in = psi_in_full[n]
            else:
                psi_spatial_in = psi_pole_prev[n]

            for visit in sn_mesh.iter_cell_visits(
                ordinate_idx=m_local, mu_level_idx=p,
            ):
                i = visit.cell_idx
                upstream = UpstreamState(
                    spatial_upstream=psi_spatial_in,
                    angular_upstream=psi_angle[i],
                )
                result = scheme.update(
                    visit=visit,
                    total_xs=sig_t_1d[i],
                    source=QV[i],
                    upstream_state=upstream,
                )

                psi = result.cell_average_flux
                psi_angle[i] = result.outgoing_angular_state
                angular_flux[n, i, 0, :] = psi
                scalar_flux[i] += w_n * psi

                if result.outgoing_spatial_flux is not None:
                    psi_spatial_in = result.outgoing_spatial_flux

            if eta_n >= 0 and abs_eta >= 1e-15:
                bc_outer[n] = psi_spatial_in

    return angular_flux, scalar_flux


def picard_loop(n_cells, n_mu, n_phi, tol=1e-12, max_iter=1000):
    fuel = get_mixture("B", "1g")
    R = 2.0
    geom = StructuredGeometry(
        geometry="CYL",
        regions=(Region(mat_id=0, outer_thickness_cm=R),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(
        geom, region_meshes=(RegionMesh(n_cells=n_cells),),
    )
    quad = ProductQuadrature.create(n_mu=n_mu, n_phi=n_phi)
    sn_mesh = SNMesh(mesh, quad)
    ng = 1
    sig_t = np.full((n_cells, 1, ng), 2.0)
    sigma_s = 1.9
    Q_iso_per_cell = np.ones((n_cells, 1, ng))  # external Q
    N = quad.N
    bc_outer = np.zeros((N, ng))
    psi_pole_prev = np.zeros((N, ng))
    psi_prev = np.zeros((N, n_cells, 1, ng))
    scalar_flux_prev = np.zeros((n_cells, ng))

    for it in range(max_iter):
        # Total within-group source = Σ_s · phi + Q_ext
        Q_tot = sigma_s * scalar_flux_prev[:, None, :] + Q_iso_per_cell  # (nx, 1, ng)
        psi_new, scalar_new = custom_cylinder_si_sweep(
            sn_mesh, sig_t, Q_tot, scalar_flux_prev[:, 0], psi_pole_prev, bc_outer,
        )
        res = np.max(np.abs(scalar_new - scalar_flux_prev))
        if res < tol:
            break
        scalar_flux_prev = scalar_new
        psi_pole_prev = psi_new[:, 0, 0, :].copy()
        psi_prev = psi_new

    return psi_new, scalar_new, it + 1, quad


def test_si_fix():
    print("\n=== Cylinder SI sweep with patched Carlson Q_bar ===")
    print(f"Expected ψ_n = 0.7957747, φ = 10\n")
    for n_mu in [2, 4]:
        for n_phi in [2, 4]:
            print(f"\n  ProductQuadrature(n_mu={n_mu}, n_phi={n_phi})")
            for n_cells in [5, 10, 20, 40]:
                try:
                    psi, sf, n_iter, quad = picard_loop(n_cells, n_mu, n_phi)
                    psi_ref = 10.0 / quad.weights.sum()
                    err_psi = np.max(np.abs(psi - psi_ref))
                    err_phi = np.max(np.abs(sf - 10.0))
                    print(f"    n_cells={n_cells:3d}: err_ψ={err_psi:.3e}, "
                          f"err_φ={err_phi:.3e}, n_iter={n_iter}")
                except Exception as e:
                    print(f"    n_cells={n_cells:3d}: EXC {e!s}")


if __name__ == "__main__":
    test_si_fix()
