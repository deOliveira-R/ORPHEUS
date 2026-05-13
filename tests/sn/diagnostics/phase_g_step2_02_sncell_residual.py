"""Diagnostic: Output B — per-cell residual at the unified SNCellOperator.

Created by numerics-investigator on 2026-05-12 for Phase G Step 2 pre-step.

The KEY discrimination of H1 vs H2:
  - Plug ψ_si (the SI fixed point) into SNCellOperator.apply at each cell.
    The residual should be ZERO if SI's fixed point is also a fixed point
    of the operator SNCellOperator wraps.
  - Plug ψ_kr (the Krylov fixed point) into SNCellOperator.apply.

If both residuals are zero (within FP), then the operators are the same,
and the drift is in inputs (Q, BC, normalisation) — supports H1.

If r_si != 0 when ψ_si is plugged in (i.e., SI's fixed point is NOT a
fixed point of the operator that Krylov uses), then the closures DIFFER
at the operator level — supports H2.

NOTE: SNCellOperator wraps DiamondDifference (the SI sweep's cell update).
So at the SI's converged ψ, the residual under SNCellOperator MUST be
zero — this is the round-trip Gate 2 of Step 1.

The interesting test is plugging ψ_kr into SNCellOperator and seeing if
that's a fixed point of the SI's cell-update math. If not, the algebraic
form of the apply-matvec residual is DIFFERENT from the algebraic form
of the SI cell-update residual.
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


def test_output_b_sncell_residual():
    kwargs = _build_sphere_2g_3reg(40)
    result_si = solve_sn(inner_solver="source_iteration", **kwargs)
    result_kr = solve_sn(inner_solver="krylov", **kwargs)

    print(f"\nN=8, nx=40, ng=2")
    print(f"keff SI={result_si.keff:.8f}  Kr={result_kr.keff:.8f}")

    # Get SN mesh + materials internals.
    from orpheus.sn.geometry import SNMesh
    sn_mesh = SNMesh(kwargs["mesh"], kwargs["quadrature"])

    print(f"  mesh: nx={sn_mesh.nx}, N_ord={sn_mesh.quad.N}, ng={result_si.scalar_flux.shape[-1]}")

    # Build sig_t per cell, per group
    # The SN solver stores per-cell material properties.
    mat_ids = sn_mesh.mat_map[:, 0]  # (nx,)
    materials = kwargs["materials"]
    ng = result_si.scalar_flux.shape[-1]
    nx = sn_mesh.nx
    sigma_t = np.zeros((nx, ng))
    for i in range(nx):
        mix = materials[mat_ids[i]]
        sigma_t[i, :] = mix.SigT

    print(f"  sigma_t shape: {sigma_t.shape}")
    print(f"  sigma_t[0]: {sigma_t[0]}, sigma_t[20]: {sigma_t[20]}, sigma_t[-1]: {sigma_t[-1]}")

    # Build within-group source at SI's converged state.
    # Q_g = q_ext (= 0 for k-eigenvalue) + Σ_s @ φ + (1/k)·χ·νΣ_f·φ
    def build_Q(scalar_flux_full, keff_val):
        """scalar_flux_full shape (nx, ny=1, ng) → returns (nx, ng)."""
        sf = scalar_flux_full[:, 0, :]  # (nx, ng)
        Q = np.zeros((nx, ng))
        for i in range(nx):
            mix = materials[mat_ids[i]]
            sigs0 = mix.SigS[0].toarray()  # (ng, ng) — Σ_s[g_from, g_to]
            chi = mix.chi  # (ng,)
            sigp = mix.SigP  # nu*SigF, shape (ng,)
            # Σ_s @ φ:  Q_g = Σ_{g'} Σ_s[g',g] · φ_{g'} = (Σ_s.T @ φ)_g
            scat = sigs0.T @ sf[i, :]
            # Fission: χ_g · (1/k) · Σ_{g'} νΣ_f,g' · φ_{g'}
            fiss = chi * (1.0 / keff_val) * (sigp @ sf[i, :])
            Q[i, :] = scat + fiss
        return Q

    Q_si = build_Q(result_si.scalar_flux, result_si.keff)
    Q_kr = build_Q(result_kr.scalar_flux, result_kr.keff)

    print(f"\nQ_si[20] = {Q_si[20]}, Q_kr[20] = {Q_kr[20]}")
    print(f"Q_si[0]  = {Q_si[0]},   Q_kr[0]  = {Q_kr[0]}")
    print(f"Q_si[-1] = {Q_si[-1]}, Q_kr[-1] = {Q_kr[-1]}")

    # ratios under different normalisations
    norm_si = result_si.scalar_flux.sum()
    norm_kr = result_kr.scalar_flux.sum()
    Q_si_n = Q_si / norm_si
    Q_kr_n = Q_kr / norm_kr
    print(f"\nNormalised:")
    print(f"  Q_si_n[20] = {Q_si_n[20]}, Q_kr_n[20] = {Q_kr_n[20]}, "
          f"rel = {100*(Q_kr_n[20]/Q_si_n[20]-1)}%")

    # ---------------------------------------------------------------
    # Apply matvec to ψ at each fixed point and check residual.
    # ---------------------------------------------------------------

    psi_si = result_si.angular_flux[:, :, 0, :]
    psi_kr = result_kr.angular_flux[:, :, 0, :]
    weights = sn_mesh.quad.weights
    weight_norm = 1.0 / weights.sum()

    # Run the actual transport_operator_matvec_spherical at psi_si
    # and psi_kr.  Both must be applied to ψ folded as scalar flux
    # (per the matvec's contract takes a "solution" vector).
    from orpheus.sn.operator import (
        transport_operator_matvec_spherical,
        build_equation_map_spherical,
    )
    print("\n=== Output B-2: apply-matvec residual ψ - source(ψ)/(1/k) ===")
    print("The matvec computes T·ψ where T = L_streaming + redist + Σ_t.")
    print("Fixed point: T·ψ = (1/k)·F·ψ + Σ_s·ψ  →  T·ψ - Σ_s·ψ = (1/k)·F·ψ")
    print("Equivalently: (T - Σ_s)·ψ - (1/k)·F·ψ = 0")

    # Build equation map + convert psi to solution vector
    eq_map = build_equation_map_spherical(nx, sn_mesh.quad, ng)
    print(f"  eq_map.n_eq = {eq_map.n_eq}")

    # Build solution vector from psi by gather at (ordinate, ix) slots.
    # eq_map skips inward-at-outer-boundary slots — those are BC-resolved.
    # psi_si shape (N, nx, ng) → gather at (ordinate[k], ix[k]) per k.
    def psi_to_solution(psi):
        sol = np.zeros((ng, eq_map.n_eq))
        for g in range(ng):
            for k in range(eq_map.n_eq):
                sol[g, k] = psi[eq_map.ordinate[k], eq_map.ix[k], g]
        return sol.ravel(order='F')

    sol_si = psi_to_solution(psi_si)
    sol_kr = psi_to_solution(psi_kr)

    print(f"  sol_si shape={sol_si.shape}, sol_kr shape={sol_kr.shape}")

    # Build face areas, volumes, alpha_half, etc.
    reduced = sn_mesh.reduced
    face_areas = reduced.face_areas
    volumes = sn_mesh.volumes
    alpha_half = sn_mesh.reduced.alpha_half
    redist_dAw = sn_mesh.reduced.redist_dAw
    tau_mm = sn_mesh.reduced.tau_mm
    bc_outer = sn_mesh.bc_right
    pole_closure = sn_mesh.pole_angular_closure

    sig_t_3d = np.zeros((nx, 1, ng))
    sig_t_3d[:, 0, :] = sigma_t

    # T·ψ via the apply matvec.
    Tpsi_si = transport_operator_matvec_spherical(
        sol_si, eq_map, sn_mesh.quad, sig_t_3d, nx, ng,
        face_areas, volumes, alpha_half, redist_dAw, tau_mm,
        sn_mesh=sn_mesh, bc_outer=bc_outer, pole_angular_closure=pole_closure,
    )
    Tpsi_kr = transport_operator_matvec_spherical(
        sol_kr, eq_map, sn_mesh.quad, sig_t_3d, nx, ng,
        face_areas, volumes, alpha_half, redist_dAw, tau_mm,
        sn_mesh=sn_mesh, bc_outer=bc_outer, pole_angular_closure=pole_closure,
    )

    print(f"  T·ψ_si: shape={Tpsi_si.shape}, ||T·ψ_si||_∞ = {np.abs(Tpsi_si).max():.4e}")
    print(f"  T·ψ_kr: shape={Tpsi_kr.shape}, ||T·ψ_kr||_∞ = {np.abs(Tpsi_kr).max():.4e}")

    # T includes ONLY the LHS of the transport equation: L_streaming +
    # redist + Σ_t · ψ.  At the fixed point, T·ψ = Σ_s·ψ + (1/k)·F·ψ.
    # So the residual we want is: T·ψ - Σ_s·ψ - (1/k)·F·ψ ≡ 0?
    # Build Σ_s·ψ + (1/k)·F·ψ in (eq_map.n_eq,) form:
    #   for each ordinate n, each cell i, each group g: Q_iso[i,g]/Σw

    def build_rhs(scalar_flux_full, keff_val):
        """RHS of T·ψ = RHS where RHS = (Σ_s·φ + (1/k)·F·φ) at ordinate n.
        For isotropic scatter/fission, RHS is the within-group isotropic
        source divided by Σw (because ψ_n on RHS appears via T·ψ_n =
        (Q_iso/Σw + L·ψ_n) → RHS_n = Q_iso/Σw)."""
        Q = build_Q(scalar_flux_full, keff_val)  # (nx, ng)
        rhs_per_ord = Q * weight_norm  # (nx, ng)
        # eq_map flattens to (ng, n_eq).  Each ordinate n at cell i contributes
        # one entry; for isotropic source rhs_per_ord[i,g] is the SAME for
        # every ordinate.  Use a broadcast.
        rhs = np.zeros((ng, eq_map.n_eq))
        # For each unknown index k: rhs[g, k] = rhs_per_ord[ix[k], g]
        # The eq_map stores (ordinate, ix, iy) per unknown index.
        for g in range(ng):
            rhs[g, :] = rhs_per_ord[eq_map.ix, g]
        return rhs.ravel(order='F')

    rhs_si = build_rhs(result_si.scalar_flux, result_si.keff)
    rhs_kr = build_rhs(result_kr.scalar_flux, result_kr.keff)

    resid_si = Tpsi_si - rhs_si
    resid_kr = Tpsi_kr - rhs_kr
    print(f"\n  resid_si = T·ψ_si − RHS(ψ_si, k_si):")
    print(f"     ||resid_si||_∞ = {np.abs(resid_si).max():.4e}")
    print(f"     ||resid_si||_2 = {np.linalg.norm(resid_si):.4e}")
    print(f"     relative to ||rhs_si||_∞: {np.abs(resid_si).max() / np.abs(rhs_si).max():.4e}")
    print(f"\n  resid_kr = T·ψ_kr − RHS(ψ_kr, k_kr):")
    print(f"     ||resid_kr||_∞ = {np.abs(resid_kr).max():.4e}")
    print(f"     ||resid_kr||_2 = {np.linalg.norm(resid_kr):.4e}")
    print(f"     relative to ||rhs_kr||_∞: {np.abs(resid_kr).max() / np.abs(rhs_kr).max():.4e}")
    print()
    print("INTERPRETATION:")
    print("  Both should be near machine zero IF both methods solve the same operator.")
    print("  resid_si large → SI's fixed point is NOT a fixed point of the apply matvec → H2")
    print("  resid_kr small AND resid_si small → both solve the same operator → H1")
    print("  resid_si large AND resid_kr small → ONLY apply is consistent → H2 (apply-form is canonical)")

    # Save for further analysis
    np.savez(
        OUT / "phase_g_step2_residual.npz",
        resid_si=resid_si,
        resid_kr=resid_kr,
        Tpsi_si=Tpsi_si,
        Tpsi_kr=Tpsi_kr,
        rhs_si=rhs_si,
        rhs_kr=rhs_kr,
        Q_si=Q_si,
        Q_kr=Q_kr,
    )
    print(f"\nSaved to {OUT/'phase_g_step2_residual.npz'}")


if __name__ == "__main__":
    test_output_b_sncell_residual()
