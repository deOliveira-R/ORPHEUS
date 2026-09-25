import time, numpy as np
from orpheus.cp.solver import CPParams, solve_cp
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import CoordSystem, Mesh1D
from orpheus.homogeneous.solver import *  # noqa
def kinf_dense(m):
    # 0-D pencil built independently from the Mixture's raw arrays: A = diag(SigT) - SigS0^T (from->to stored [from,to]), F = chi (x) SigP
    S0 = np.asarray(m.SigS[0].todense() if hasattr(m.SigS[0], "todense") else m.SigS[0])
    A = np.diag(m.SigT) - S0.T
    F = np.outer(m.chi, m.SigP)
    ev, vec = np.linalg.eig(np.linalg.solve(A, F)); i = np.argmax(ev.real)
    v = np.abs(vec[:, i].real); return ev[i].real, v / v.sum()
for ng in ("1g", "2g", "4g"):
    m = get_mixture("A", ng); ki, spec = kinf_dense(m)
    for coord in (CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL):
        for R, n in ((0.5, 5), (2.0, 10), (10.0, 20)):
            edges = np.linspace(0, R, n + 1)
            mesh = Mesh1D(edges=edges, mat_ids=np.full(n, 2), coord=coord)
            t0 = time.perf_counter()
            r = solve_cp({2: m}, mesh, CPParams(keff_tol=1e-12, flux_tol=1e-12, max_outer=2000))
            dt = time.perf_counter() - t0
            phi = r.flux  # (n, ng)
            flat = np.max(np.abs(phi / phi.mean(axis=0) - 1))
            sp = phi.mean(axis=0) / phi.mean(axis=0).sum()
            print(f"{ng} {coord.name[:3]} R={R:5} k-kinf rel={abs(r.keff-ki)/ki:.2e} flat={flat:.2e} spec={np.max(np.abs(sp-spec)):.2e} kinf={ki:.10f} {dt:.2f}s it={len(r.keff_history)}")
