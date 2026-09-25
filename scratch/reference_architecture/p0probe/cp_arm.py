import io, contextlib, numpy as np
import orpheus.cp.solver as cps
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import CoordSystem, Mesh1D
m = get_mixture("A", "2g")
def run(tag):
    for coord in (CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL):
        mesh = Mesh1D(edges=np.linspace(0, 2.0, 11), mat_ids=np.full(10, 2), coord=coord)
        with contextlib.redirect_stdout(io.StringIO()):
            r = cps.solve_cp({2: m}, mesh, cps.CPParams(keff_tol=1e-12, flux_tol=1e-12, max_outer=2000))
        phi = r.flux; flat = np.max(np.abs(phi/phi.mean(axis=0)-1))
        print(f"{tag:14} {coord.name[:3]} k/kinf-1={r.keff/1.875-1:+.3e} flat={flat:.2e}")
run("honest")
orig = cps.CPMesh._apply_white_bc
cps.CPMesh._apply_white_bc = lambda self, P, s: P.copy()          # arm A: white re-entry dropped (vacuum)
run("arm A vacuum")
def half(self, P, s):                                               # arm B: re-entry halved (albedo 0.5)
    Pi = orig(self, P, s); return P + 0.5*(Pi - P)
cps.CPMesh._apply_white_bc = half; run("arm B albedo.5")
def clamp(self, P, s):                                              # arm C: row-sum defect in P_cell (1e-3 extra on the outer row) 
    P2 = P.copy(); P2[-1, :] *= 1.001; return orig(self, P2, s)
cps.CPMesh._apply_white_bc = clamp; run("arm C rowsum")
cps.CPMesh._apply_white_bc = orig
