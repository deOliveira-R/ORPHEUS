"""Isolate the two algebraic differences between SI and apply-matvec.

This script implements four variants of the SI sweep with toggleable
fixes and compares each to ψ_K (the apply-matvec converged answer):

  variant 0: production SI (baseline; should match ψ_SI*)
  variant 1: SI + pole-face IC fix (apply matvec's psi_face_in[i=0] = psi_cell)
  variant 2: SI + Carlson seed Q_bar = 0.5 · Σ_t · φ_0  (not 0.5 · Q_scatt)
  variant 3: SI + both fixes simultaneously

If variant 3 matches ψ_K to machine precision, the SI's two algebraic
defects (pole-face IC + Q_bar source mismatch) together fully account
for the 22% L0 error.

If variant 3 still doesn't match, there's a third effect.
"""
from __future__ import annotations

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.geometry import SNMesh
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.solver import solve_sn_fixed_source
from orpheus.sn.spatial.psi_half_angle_seed import carlson_inward_sweep_from_source


def custom_si_sweep(
    Q_scatt: np.ndarray,         # (nx,) — scattering source per cell
    Q_ext: np.ndarray,           # (nx,) — external source per cell
    sigma_t: float,
    sn_mesh: SNMesh,
    quad,
    psi_state: np.ndarray,       # (N, nx) — previous-iter angular flux
    bc_outer: np.ndarray,        # (N, ng) buffer
    *,
    fix_pole_face_ic: bool = False,
    fix_carlson_seed_source: bool = False,
):
    """SI sweep with optional fixes for the two algebraic defects.

    `fix_pole_face_ic`: when True, set psi_spatial_in at the pole face
        (μ≥0, i=0) to psi_state[n, 0] instead of 0.
    `fix_carlson_seed_source`: when True, build Q_bar from
        0.5 · σ_t · φ_0 (apply-matvec style) instead of 0.5 · Q_scatt.
    """
    N = quad.N
    nx = sn_mesh.nx
    mu = quad.mu_x
    w = quad.weights
    sum_w = w.sum()
    red = sn_mesh.reduced
    A = red.face_areas
    V = sn_mesh.volumes[:, 0]
    alpha_h = red.alpha_half
    tau_mm = red.tau_mm
    dAw = red.redist_dAw
    bc_outer_obj = sn_mesh.bc_right

    # Carlson seed
    inflow_full = bc_outer_obj.apply(bc_outer)
    most_inward = int(np.argmin(mu))
    bc_outer_value = inflow_full[most_inward, :]
    if fix_carlson_seed_source:
        # Apply-matvec style: 0.5 · σ_t · φ_0
        phi_0 = (w[:, None] * psi_state).sum(axis=0)  # (nx,)
        Q_bar = 0.5 * sigma_t * phi_0[None, :]  # (ng=1, nx)
    else:
        # SI style: 0.5 · Q_scatt (no external)
        Q_bar = 0.5 * Q_scatt[None, :]
    sigma_t_gx = np.full((1, nx), sigma_t)
    dr = sn_mesh.dx
    phi_aux = carlson_inward_sweep_from_source(
        Q_bar=Q_bar, sigma_t=sigma_t_gx, dr=dr, bc_outer_value=bc_outer_value,
    )  # (ng, nx)
    psi_angle = phi_aux[0, :].copy()  # (nx,)

    # Volumetric source per cell, per ordinate
    QV = (Q_scatt + Q_ext) * V / sum_w  # (nx,)

    psi_out = np.zeros((N, nx))
    bc_out = bc_outer.copy()

    for n in range(N):
        mu_n = mu[n]
        abs_mu = abs(mu_n)
        if mu_n < 0:
            psi_spatial_in = bc_outer_obj.apply(bc_out)[n, 0]
            cell_order = range(nx - 1, -1, -1)
        else:
            # SI default: 0; Pole-face IC fix: psi_state[n, 0]
            psi_spatial_in = psi_state[n, 0] if fix_pole_face_ic else 0.0
            cell_order = range(nx)
        for i in cell_order:
            A_inner = A[i]
            A_outer = A[i + 1]
            A_total = A_inner + A_outer
            A_downstream = A_outer if mu_n >= 0 else A_inner
            tau = tau_mm[n]
            c_in = (1.0 - tau) / tau * alpha_h[n + 1] + alpha_h[n]
            c_out = alpha_h[n + 1] / tau
            denom = 2.0 * abs_mu * A_downstream + dAw[i, n] * c_out + sigma_t * V[i]
            numer = (QV[i]
                     + abs_mu * A_total * psi_spatial_in
                     + dAw[i, n] * c_in * psi_angle[i])
            psi_avg = numer / denom
            psi_out[n, i] = psi_avg
            # WDD spatial update
            psi_spatial_out = 2.0 * psi_avg - psi_spatial_in
            # M-M angular update at THIS cell
            psi_angle[i] = (psi_avg - (1.0 - tau) * psi_angle[i]) / tau
            psi_spatial_in = psi_spatial_out
        if mu_n >= 0:
            bc_out[n, 0] = psi_spatial_in

    return psi_out, bc_out


def picard_to_fp(fix_pole, fix_seed, n_cells=2, n_ord=2, n_picard=2000, tol=1e-14):
    fuel = get_mixture("B", "1g")
    geom = StructuredGeometry(geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),), bcs=(BC.reflective,))
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=n_cells),))
    quad = GaussLegendre1D.create(n_ordinates=n_ord)
    sn_mesh = SNMesh(mesh, quad)
    sigma_t = float(fuel.SigT[0])
    sigma_s = float(fuel.SigS[0].toarray()[0, 0])
    N = quad.N
    nx = sn_mesh.nx

    psi = np.zeros((N, nx))
    bc_outer = np.zeros((N, 1))
    Q_ext = np.ones(nx)
    prev = 0.0
    for k in range(n_picard):
        phi_0 = (quad.weights[:, None] * psi).sum(axis=0)
        Q_scatt = sigma_s * phi_0
        psi_new, bc_outer = custom_si_sweep(
            Q_scatt, Q_ext, sigma_t, sn_mesh, quad, psi, bc_outer,
            fix_pole_face_ic=fix_pole, fix_carlson_seed_source=fix_seed,
        )
        norm = np.linalg.norm(psi_new)
        if norm > 0:
            r = np.linalg.norm(psi_new - psi) / norm
            psi = psi_new
            if r < tol:
                return psi, k + 1, r
        psi = psi_new
        prev = norm
    return psi, n_picard, r


def test_isolate_two_bugs():
    """Pin the two algebraic differences between SI sweep and apply-matvec.

    The custom SI sweep with both fixes applied MUST recover ψ_K = 5
    (the analytical streaming-equilibrium answer) to machine precision.
    If either fix is missing, the sweep converges to a wrong fixed point.
    """
    print("Variants of SI sweep on 2-cell 2-ordinate sphere, mixture B 1G:")
    psi_ref = np.full((2, 2), 5.0)  # ψ_K analytical
    errors = {}
    for fix_pole, fix_seed, label in [
        (False, False, "v0: production SI (baseline)"),
        (True,  False, "v1: + pole-face IC fix"),
        (False, True,  "v2: + Carlson seed source fix"),
        (True,  True,  "v3: + BOTH fixes"),
    ]:
        psi, n_iter, r = picard_to_fp(fix_pole, fix_seed)
        err = np.max(np.abs(psi - psi_ref))
        errors[(fix_pole, fix_seed)] = err
        print(f"\n  {label}:")
        print(f"    n_picard={n_iter}, residual={r:.3e}")
        print(f"    ψ:\n{psi}")
        print(f"    max |ψ − 5.0| = {err:.4e}")

    # Pin the algebraic-difference fingerprints
    # v0: production SI — ~1.88 (mesh-independent O(1) error)
    assert errors[(False, False)] > 1.0, (
        f"v0 baseline lost its expected O(1) error fingerprint: "
        f"{errors[(False, False)]:.4e}; this diagnostic assumes the "
        f"production SI's two-defect signature."
    )
    # v3: both fixes — machine precision
    assert errors[(True, True)] < 1e-10, (
        f"v3 (both fixes) failed to recover ψ_K to machine precision: "
        f"{errors[(True, True)]:.4e}. Either the algebraic-difference "
        f"identification is incomplete, or one of the custom-sweep "
        f"fixes drifted from the apply-matvec semantics."
    )
    # v1: pole-face IC fix alone substantially improves (1.88 → ~0.06)
    assert errors[(True, False)] < 0.1, (
        f"v1 (pole-face IC fix only) max err = {errors[(True, False)]:.4e}; "
        f"expected significant reduction relative to v0."
    )
    # Pole-face IC defect is DOMINANT over Carlson-seed defect
    assert errors[(True, False)] < errors[(False, True)], (
        "Pole-face IC defect is expected to dominate the Carlson-seed "
        "defect at this configuration. If this assertion fails, the "
        "dominance relationship has reversed and the analysis must be "
        "re-evaluated."
    )


if __name__ == "__main__":
    test_isolate_two_bugs()
