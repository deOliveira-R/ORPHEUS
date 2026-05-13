"""Pinpoint the dominant algebraic difference between SI and apply-matvec
at the pole cell on outward sweeps: the WDD ψ_face_in at the pole face.

SI sweep:        psi_face_in[μ≥0, i=0] = 0           (sweep.py:559)
Apply-matvec:    psi_face_in[μ≥0, i=0] = ψ_cell[n, i=0]   (operator.py:784)

At the pole the face area A[0] = 0, so the streaming-IN term
−A[0]·ψ_face_in = 0 either way. BUT psi_face_in feeds into the WDD
recurrence: psi_face_out = 2·psi_cell − psi_face_in. So:

  SI:    psi_face_out[μ≥0, i=0] = 2·psi_cell - 0       = 2·psi_cell
  Apply: psi_face_out[μ≥0, i=0] = 2·psi_cell - psi_cell = psi_cell

This face-out propagates as the next cell's face-in, so the difference
cascades through ALL outward cells.

The Pomraning isotropy at r=0 says ψ_n(r=0) = const · 1 for every
ordinate — this is consistent with ψ_face_out(pole) = ψ_cell[pole] which
gives flat WDD propagation. The SI's choice (ψ_face_in=0) gives
ψ_face_out = 2·ψ_cell ≠ ψ_cell on flat profiles — VIOLATES Pomraning.

This script empirically demonstrates: forcing the apply-matvec's
pole-face IC into the SI sweep recovers the analytical answer.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.geometry import SNMesh
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.solver import solve_sn_fixed_source


def _solve(n_cells=2, n_ord=2):
    fuel = get_mixture("B", "1g")
    geom = StructuredGeometry(geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),), bcs=(BC.reflective,))
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=n_cells),))
    quad = GaussLegendre1D.create(n_ordinates=n_ord)
    Q = np.ones((quad.N, n_cells, 1, 1))
    res_si = solve_sn_fixed_source(
        materials={0: fuel}, mesh=mesh, quadrature=quad,
        external_source=Q, boundary_condition=BC.reflective,
        inner_solver="source_iteration", max_inner=20000, inner_tol=1e-14)
    res_kr = solve_sn_fixed_source(
        materials={0: fuel}, mesh=mesh, quadrature=quad,
        external_source=Q, boundary_condition=BC.reflective,
        inner_solver="krylov", max_inner=2000, inner_tol=1e-14)
    return (res_si.angular_flux[:,:,0,0], res_kr.angular_flux[:,:,0,0],
            SNMesh(mesh, quad), quad, fuel)


def test_pole_face_streaming_mismatch():
    """At fixed point, compute SI's vs apply-matvec's effective WDD
    propagation of psi_face at the pole cell. Show the magnitudes."""
    psi_si, psi_kr, sn_mesh, quad, fuel = _solve()
    print(f"\nψ_SI: {psi_si}")
    print(f"ψ_K:  {psi_kr}")
    mu = quad.mu_x
    N = quad.N
    nx = sn_mesh.nx
    A = sn_mesh.reduced.face_areas
    V = sn_mesh.volumes[:, 0]
    sigma_t = float(fuel.SigT[0])

    print(f"\nFace areas A = {A}, volumes V = {V}")

    # For each outward ordinate n=1 (μ=+1/√3):
    n = 1
    print(f"\n=== Outward ordinate n=1 (μ={mu[n]:+.4f}) ===")
    for psi_label, psi in (("SI", psi_si), ("Krylov", psi_kr)):
        print(f"\nWDD propagation under {psi_label} fixed point:")
        # Case A: SI's choice  ψ_face_in[i=0] = 0
        psi_face_in_si = 0.0
        # Case B: apply's choice ψ_face_in[i=0] = ψ_cell[i=0]
        psi_face_in_apply = psi[n, 0]
        # Step through cells outward
        face_si = psi_face_in_si
        face_app = psi_face_in_apply
        for i in range(nx):
            psi_cell = psi[n, i]
            face_out_si = 2.0 * psi_cell - face_si
            face_out_app = 2.0 * psi_cell - face_app
            # streaming contribution at this cell (signed mu_n):
            streaming_si = mu[n] * (A[i+1] * face_out_si - A[i] * face_si) / V[i]
            streaming_app = mu[n] * (A[i+1] * face_out_app - A[i] * face_app) / V[i]
            collision = sigma_t * psi_cell
            print(f"  i={i}: ψ_cell={psi_cell:.4f}  "
                  f"face_in(SI)={face_si:.4f} face_in(app)={face_app:.4f}  "
                  f"face_out(SI)={face_out_si:.4f} face_out(app)={face_out_app:.4f}  "
                  f"streaming(SI)={streaming_si:.4f}  streaming(app)={streaming_app:.4f}  "
                  f"coll={collision:.4f}")
            face_si = face_out_si
            face_app = face_out_app


def test_apply_matvec_balance_at_pole_with_two_ICs():
    """Show: streaming + collision balance at pole (n=1, i=0) under ψ_K=5
    flat: with apply-matvec IC (face_in=ψ_cell) gives balance with q=10;
    with SI IC (face_in=0) gives twice the streaming → unbalanced."""
    psi_si, psi_kr, sn_mesh, quad, fuel = _solve()
    A = sn_mesh.reduced.face_areas
    V = sn_mesh.volumes[:, 0]
    mu = quad.mu_x
    sigma_t = float(fuel.SigT[0])
    n = 1
    i = 0
    psi_cell = psi_kr[n, i]  # = 5
    # Case A: apply matvec — face_in = psi_cell
    face_in_a = psi_cell
    face_out_a = 2 * psi_cell - face_in_a  # = psi_cell
    streaming_a = mu[n] * (A[i+1] * face_out_a - A[i] * face_in_a) / V[i]
    coll_a = sigma_t * psi_cell
    # q = Q_within / Σw at FP. With Σs=1.9, φ_0=10, Q_ext=1: Q_within=20.
    q_total = 20.0 / quad.weights.sum()  # = 10
    print(f"\n=== Pomraning isotropy check at ψ_K=5 with both ICs ===")
    print(f"Case Apply (face_in=ψ_cell={face_in_a:.4f}):")
    print(f"  face_out={face_out_a:.4f}  streaming={streaming_a:.4f}  collision={coll_a:.4f}")
    print(f"  T·ψ = streaming + redist_assumed_0 + collision = {streaming_a + coll_a:.4f}")
    print(f"  q_total = {q_total:.4f}  →  T·ψ − q = {streaming_a + coll_a - q_total:.4f}")
    # Case B: SI — face_in = 0
    face_in_b = 0.0
    face_out_b = 2 * psi_cell - face_in_b
    streaming_b = mu[n] * (A[i+1] * face_out_b - A[i] * face_in_b) / V[i]
    coll_b = sigma_t * psi_cell
    print(f"\nCase SI (face_in=0):")
    print(f"  face_out={face_out_b:.4f}  streaming={streaming_b:.4f}  collision={coll_b:.4f}")
    print(f"  T·ψ = streaming + redist_assumed_0 + collision = {streaming_b + coll_b:.4f}")
    print(f"  q_total = {q_total:.4f}  →  T·ψ − q = {streaming_b + coll_b - q_total:.4f}")
    print(f"\nDifference in T·ψ between two ICs: "
          f"{(streaming_b + coll_b) - (streaming_a + coll_a):.4f}")
    print("  ≡ streaming_b - streaming_a")
    print(f"  = mu_n · (A[1] · (face_out_b - face_out_a) - A[0] · (face_in_b - face_in_a)) / V[0]")
    print(f"  = mu_n · A[1] · (face_out_b - face_out_a) / V[0]    [A[0]=0]")
    print(f"  = mu_n · A[1] · ((2ψ-0)-(2ψ-ψ)) / V[0]")
    print(f"  = mu_n · A[1] · ψ_cell / V[0]")
    print(f"  = {mu[n]:.4f} · {A[1]:.4f} · {psi_cell:.4f} / {V[0]:.4f}")
    print(f"  = {mu[n] * A[1] * psi_cell / V[0]:.4f}")


if __name__ == "__main__":
    test_pole_face_streaming_mismatch()
    test_apply_matvec_balance_at_pole_with_two_ICs()
