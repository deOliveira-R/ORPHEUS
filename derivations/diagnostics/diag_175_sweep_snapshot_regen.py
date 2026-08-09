"""Issue #175 — regenerate ``sweep_ref_2g.npy`` with independent cross-check.

The frozen snapshot consumed by
``tests/sn/operators/test_solver_components.py::TestTransportSweep::
test_matches_saved_reference`` was generated under the pre-Wave-2
sweep and drifted past ``rtol=1e-14`` through the Wave-T / Phase-5
reduction-tree refactors (xfail'd since the taxonomy reorg,
``b06bc49``). Per vv-principles §bit-identity, regenerating a snapshot
requires verifying the NEW values against a structurally-independent
reference — old-vs-new distance proves nothing when the algorithm
deliberately changed.

The independent reference here is a from-scratch per-cell 2-D
diamond-difference sweep written below with PLAIN PYTHON LOOPS:

* no shared code with the production windowed-frontier kernel
  (``_MovingFrontier`` / ``cell_kernel_batch`` /
  ``sweep_octant_group``) — different traversal order (per-ordinate
  row-column loops vs batched anti-diagonal wavefront), different
  recurrence transcription, no shared helpers;
* same discrete model by construction: DD closure
  ``psi_out = 2 psi_c − psi_in``, vacuum inflow, per-ordinate source
  ``Q/W`` (producer-side Pattern-7 projection, lesson L18).

Agreement at float-precision between the two establishes L14 leg 1
(production ≡ structurally-independent reference) for the snapshot
values; the broader correctness of the production sweep family is
anchored by the streaming-equilibrium L0 gates and the 2-D MMS L1
convergence suite.

Run:  .venv/bin/python derivations/diagnostics/diag_175_sweep_snapshot_regen.py [--write]
"""

from __future__ import annotations

import sys

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.loss_representation import transport_sweep
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.source_sinks import AngularSourceSink


def build_fixture() -> tuple[SNMesh, np.ndarray, np.ndarray]:
    """Replicate the ``solver_2g`` fixture + seeded source of the test."""
    fuel = get_mixture("A", "2g")
    mod = get_mixture("B", "2g")
    materials = {2: fuel, 0: mod}

    nx, ny, delta = 6, 4, 0.2
    mat = np.zeros((nx, ny), dtype=int)
    mat[:3, :] = 2
    mat[3:, :] = 0
    mesh = Mesh2D(
        edges_x=np.linspace(0, nx * delta, nx + 1),
        edges_y=np.linspace(0, ny * delta, ny + 1),
        mat_map=mat,
    )
    quad = Quadrature.lebedev(order=17)
    sn_mesh = SNMesh(mesh, quad, materials)

    # σ_t exactly as the test consumes it (solver.mat_xs.total_cross_section).
    from orpheus.sn.solver import SNSolver

    solver = SNSolver(sn_mesh)
    sig_t = solver.mat_xs.total_cross_section
    ng = sig_t.shape[0]

    np.random.seed(7)
    Q = np.random.rand(ng, nx, ny) + 0.01
    return sn_mesh, sig_t, Q


def hand_sweep(
    sn_mesh: SNMesh, sig_t: np.ndarray, Q: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Per-cell-loop 2-D DD sweep, vacuum inflow on every face.

    Solves (Ω·∇ + σ_t) ψ_n = Q/W per ordinate with the DD closure.
    Returns (psi (N, ng, nx, ny), phi (ng, nx, ny) = Σ_n w_n ψ_n).
    """
    quad = sn_mesh.quad
    N = quad.N
    ng, nx, ny = Q.shape
    W = float(quad.weights.sum())
    dx, dy = (np.asarray(w) for w in sn_mesh.axis_widths)

    qn = Q / W  # per-ordinate magnitude, identical across n
    psi = np.zeros((N, ng, nx, ny))

    for n in range(N):
        mu = float(quad.mu_x[n])
        eta = float(quad.mu_y[n])
        xs = range(nx) if mu >= 0 else range(nx - 1, -1, -1)
        ys = range(ny) if eta >= 0 else range(ny - 1, -1, -1)
        for g in range(ng):
            # Vacuum: zero inflow at the upstream domain faces.
            inflow_y = np.zeros(nx)  # ψ on the upstream y-face, per x-column
            for j in ys:
                inflow_x = 0.0  # ψ on the upstream x-face for this row
                for i in xs:
                    cx = 2.0 * abs(mu) / dx[i]
                    cy = 2.0 * abs(eta) / dy[j]
                    psi_c = (qn[g, i, j] + cx * inflow_x + cy * inflow_y[i]) / (
                        sig_t[g, i, j] + cx + cy
                    )
                    psi[n, g, i, j] = psi_c
                    inflow_x = 2.0 * psi_c - inflow_x
                    inflow_y[i] = 2.0 * psi_c - inflow_y[i]

    phi = np.einsum("n,ngxy->gxy", quad.weights, psi)
    return psi, phi


def main() -> int:
    sn_mesh, sig_t, Q = build_fixture()

    ang, phi_prod = transport_sweep(
        AngularSourceSink.from_isotropic(Q, sn_mesh),
        sig_t,
        sn_mesh,
        AngularBoundaryFlux.zeros_on(sn_mesh),
    )
    psi_hand, phi_hand = hand_sweep(sn_mesh, sig_t, Q)

    psi_prod = ang if isinstance(ang, np.ndarray) else ang.values
    dpsi = np.abs(psi_prod - psi_hand).max()
    dphi_abs = np.abs(phi_prod - phi_hand).max()
    dphi_rel = (np.abs(phi_prod - phi_hand) / np.abs(phi_hand)).max()

    print(f"max |ψ_prod − ψ_hand|        = {dpsi:.3e}")
    print(f"max |φ_prod − φ_hand|        = {dphi_abs:.3e}")
    print(f"max rel |φ_prod − φ_hand|/φ  = {dphi_rel:.3e}")

    ok = dphi_rel < 1e-12
    print(f"cross-check: {'PASS' if ok else 'FAIL'} (gate: rel < 1e-12)")

    if "--write" in sys.argv:
        if not ok:
            print("REFUSING to write snapshot: cross-check failed.")
            return 1
        from pathlib import Path

        out = Path(__file__).parents[2] / "tests" / "sn" / "sweep_ref_2g.npy"
        np.save(out, phi_prod)
        print(f"wrote {out}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
