"""Classify the spherical CP conservation break: quadrature artifact vs real.

Sweeps total optical thickness and n_quad_y for a solid-ball spherical mesh,
reporting the white-BC conservation residual max|rowsum - 1|.
"""
import numpy as np
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.cp.solver import CPMesh, CPParams

def sphere_cons(N, sig_t, n_quad_y, cell_width=1.0):
    edges = np.arange(N + 1, dtype=float) * cell_width
    mesh = Mesh1D(edges=edges, mat_ids=np.zeros(N, dtype=int),
                  coord=CoordSystem.SPHERICAL, bc_right=BC("white"))
    cp = CPMesh(mesh, CPParams(n_quad_y=n_quad_y))
    P = cp.compute_pinf_group(np.full(N, float(sig_t)))
    rs = P.sum(axis=1)
    tau = N * sig_t * cell_width
    return tau, float(np.max(np.abs(rs - 1.0))), float(rs.min()), float(rs.max())

print("=== SPHERE conservation: tau sweep at n_quad_y=64, N=25, width=1 ===")
print(f"{'tau_tot':>8s} {'cons_resid':>12s} {'rowsum_min':>11s} {'rowsum_max':>11s}")
for sig_t in [0.01, 0.04, 0.12, 0.2, 0.4, 1.0]:
    tau, resid, rmin, rmax = sphere_cons(25, sig_t, 64)
    print(f"{tau:8.3g} {resid:12.3e} {rmin:11.4f} {rmax:11.4f}")

print("\n=== SPHERE conservation: n_quad_y sweep at tau_total=25 (N=25,sig_t=1) ===")
print(f"{'n_quad_y':>8s} {'cons_resid':>12s}")
for nq in [64, 128, 256, 512, 1024]:
    tau, resid, rmin, rmax = sphere_cons(25, 1.0, nq)
    print(f"{nq:8d} {resid:12.3e}")

print("\n=== SPHERE conservation: n_quad_y sweep at tau_total=5 (N=25,sig_t=0.2) ===")
print(f"{'n_quad_y':>8s} {'cons_resid':>12s}")
for nq in [64, 128, 256, 512]:
    tau, resid, rmin, rmax = sphere_cons(25, 0.2, nq)
    print(f"{nq:8d} {resid:12.3e}")

# Also: refine radial mesh at fixed tau_total=25 (smaller cells => smaller cell tau)
print("\n=== SPHERE conservation: radial refinement at tau_total=25, n_quad_y=256 ===")
print(f"{'N':>5s} {'cell_w':>7s} {'tau_cell':>9s} {'cons_resid':>12s}")
for N in [25, 50, 100]:
    w = 25.0 / N
    tau, resid, rmin, rmax = sphere_cons(N, 1.0, 256, cell_width=w)
    print(f"{N:5d} {w:7.3f} {1.0*w:9.3f} {resid:12.3e}")
