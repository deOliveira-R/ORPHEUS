"""Isolate the spherical CP rowsum>1 excess: P_cell (pre-BC) vs diagonal/off-diag."""
import numpy as np
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.cp.solver import CPMesh, CPParams

def probe(N, sig_t, coord, w=1.0, nq=256):
    edges = np.arange(N + 1, dtype=float) * w
    mesh = Mesh1D(edges=edges, mat_ids=np.zeros(N, dtype=int),
                  coord=coord, bc_right=BC("vacuum"))  # vacuum => returns pure P_cell
    cp = CPMesh(mesh, CPParams(n_quad_y=nq))
    P = cp.compute_pinf_group(np.full(N, float(sig_t)))
    rs = P.sum(axis=1)
    diag = np.diag(P).copy()
    offrs = rs - diag
    return P, rs, diag, offrs

print("=== PURE P_cell (vacuum) rowsums — should be <= 1 (rowsum + escape = 1) ===")
for coord, name in [(CoordSystem.SPHERICAL, "sph"), (CoordSystem.CYLINDRICAL, "cyl"),
                    (CoordSystem.CARTESIAN, "slab")]:
    for sig_t in [0.2, 1.0]:
        P, rs, diag, offrs = probe(25, sig_t, coord)
        tau = 25 * sig_t
        print(f"  {name} tau_tot={tau:5.1f}: rowsum[min,max]=[{rs.min():.4f},{rs.max():.4f}]  "
              f"diag[min,max]=[{diag.min():.4f},{diag.max():.4f}]  n(rs>1)={int(np.sum(rs>1+1e-9))}")

print("\n=== SPHERE: which rows exceed 1? (tau_total=25, N=25) ===")
P, rs, diag, offrs = probe(25, 1.0, CoordSystem.SPHERICAL)
print("  row  radius_out  P_diag    off_rowsum  total_rowsum")
for i in [0, 1, 2, 5, 12, 20, 23, 24]:
    print(f"  {i:3d}  {float(i+1):9.1f}  {diag[i]:.5f}   {offrs[i]:.5f}    {rs[i]:.5f}")

print("\n=== SPHERE self-collision diagonal alone: is P_ii > 1 possible? ===")
# P_ii is the self-collision probability; must be <= 1 physically.
P, rs, diag, offrs = probe(25, 1.0, CoordSystem.SPHERICAL)
print(f"  max P_ii = {diag.max():.5f}  (rows where P_ii>1: {int(np.sum(diag>1+1e-9))})")
