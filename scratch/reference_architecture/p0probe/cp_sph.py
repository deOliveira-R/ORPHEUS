import io, contextlib, numpy as np
from orpheus.cp.solver import CPParams, solve_cp
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import CoordSystem, Mesh1D
m = get_mixture("A", "1g")
for R, n in ((3.0,10),(4.0,10),(5.0,10),(6.0,10),(8.0,10),(10.0,10),(10.0,5),(10.0,40),(5.0,20),(5.0,40),(20.0,20)):
    for coord in (CoordSystem.SPHERICAL,):
        mesh = Mesh1D(edges=np.linspace(0, R, n+1), mat_ids=np.full(n, 2), coord=coord)
        with contextlib.redirect_stdout(io.StringIO()):
            r = solve_cp({2: m}, mesh, CPParams(keff_tol=1e-12, flux_tol=1e-12, max_outer=2000))
        phi = r.flux[:, 0]; print(f"SPH R={R:5} n={n:3} dr={R/n:.2f} flat={np.max(np.abs(phi/phi.mean()-1)):.2e} k={r.keff:.12f} argmax_dev_cell={np.argmax(np.abs(phi/phi.mean()-1))}")
