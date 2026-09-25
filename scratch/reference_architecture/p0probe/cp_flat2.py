import io, contextlib, numpy as np, inspect
from orpheus.cp.solver import CPParams, solve_cp
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import CoordSystem, Mesh1D
from orpheus.geometry.mesh import BC
print(inspect.signature(Mesh1D))
m = get_mixture("A", "2g")
S0 = np.asarray(m.SigS[0].todense() if hasattr(m.SigS[0], "todense") else m.SigS[0])
ev = np.linalg.eigvals(np.linalg.solve(np.diag(m.SigT) - S0.T, np.outer(m.chi, m.SigP))); ki = ev.real.max()
for coord in (CoordSystem.CARTESIAN, CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL):
    for R, kind in ((0.5, "geom"), (4.0, "geom"), (2.0, "one-cell")):
        if kind == "geom":
            edges = R * (1 - np.geomspace(1, 1e-2, 12)) / (1 - 1e-2); edges[0] = 0.0; edges[-1] = R
        else:
            edges = np.array([0.0, R])
        n = len(edges) - 1
        kw = dict(edges=edges, mat_ids=np.full(n, 2), coord=coord)
        if coord == CoordSystem.CARTESIAN: kw.update(bc_left=BC("white"), bc_right=BC("white"))
        try:
            mesh = Mesh1D(**kw)
            with contextlib.redirect_stdout(io.StringIO()):
                r = solve_cp({2: m}, mesh, CPParams(keff_tol=1e-12, flux_tol=1e-12, max_outer=2000))
            phi = r.flux; print(f"{coord.name[:3]} R={R} {kind:8} n={n} k/kinf-1={r.keff/ki-1:+.2e} flat={np.max(np.abs(phi/phi.mean(axis=0)-1)):.2e}")
        except Exception as e:
            print(coord.name[:3], R, kind, "ERR", type(e).__name__, str(e)[:120])
