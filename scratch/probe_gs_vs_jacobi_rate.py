"""#341 — is the boundary-G-S rate INVERSION at d=3 general, or an artifact of
the zero-leakage pure-absorber extreme it was first measured on?

The docstring at solver.py:3056 claims a "~0.86-0.92x rate gain" for
gauss_seidel, and solve_sn_fixed_source DEFAULTS to it.  The d3 diagnosis
measured G-S 2x SLOWER than Jacobi at d=3 -- but only on the all-reflective
pure absorber (zero leakage, zero scattering), which is the pathological
corner.  Changing a production default on one extreme point would be
over-fitting.  So: sweep dimension x leakage x scattering.
"""
from __future__ import annotations

import numpy as np

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn_fixed_source
from orpheus.transport.mesh.axis import AxisMesh


def axes(extents, cells, bcs):
    return tuple(
        AxisMesh(edges=np.linspace(0.0, e, n + 1), bc_low=lo, bc_high=hi)
        for e, n, (lo, hi) in zip(extents, cells, bcs)
    )


def absorber():
    return make_mixture(
        sig_t=np.array([0.8, 1.6]), sig_c=np.array([0.8, 1.6]),
        sig_f=np.array([0.0, 0.0]), nu=np.array([0.0, 0.0]),
        chi=np.zeros(2), sig_s=np.zeros((2, 2)),
    )


def scatterer(c=0.5):
    # c = Sigma_s / Sigma_t per group; downscatter in the (0->1) slot so the
    # groups are genuinely coupled (a Mode-6 convention catcher, and the
    # realistic regime the default is meant to serve).
    sig_t = np.array([0.8, 1.6])
    s = np.zeros((2, 2))
    s[0, 0] = c * sig_t[0] * 0.7
    s[1, 0] = c * sig_t[0] * 0.3          # 0 -> 1 downscatter
    s[1, 1] = c * sig_t[1]
    sig_c = sig_t - s.sum(axis=0)
    return make_mixture(
        sig_t=sig_t, sig_c=sig_c, sig_f=np.array([0.0, 0.0]),
        nu=np.array([0.0, 0.0]), chi=np.zeros(2), sig_s=s,
    )


V, R = BC("vacuum"), BC("reflective")
CONFIGS = [
    # (label, extents, cells, bcs, mixture)
    # CONTROL FIRST (vv anti-#17): this row must reproduce the diagnosis's
    # 1631 G-S sweeps, or the instrument is not measuring what I think.
    ("d3 all-reflective absorber [CTRL]", (1.0, 2.0, 3.0), (3, 4, 5),
     [(R, R)] * 3, absorber()),
    ("d2 all-reflective absorber", (1.0, 2.0), (3, 4),
     [(R, R)] * 2, absorber()),
    ("d2 all-reflective c=0.5", (1.0, 2.0), (3, 4),
     [(R, R)] * 2, scatterer()),
    ("d2 x-VACUUM c=0.5", (1.0, 2.0), (3, 4),
     [(V, V), (R, R)], scatterer()),
    ("d3 all-reflective c=0.5", (1.0, 2.0, 3.0), (3, 4, 5),
     [(R, R)] * 3, scatterer()),
    ("d3 x-VACUUM c=0.5", (1.0, 2.0, 3.0), (3, 4, 5),
     [(V, V), (R, R), (R, R)], scatterer()),
    ("d3 x,y-VACUUM c=0.5", (1.0, 2.0, 3.0), (3, 4, 5),
     [(V, V), (V, V), (R, R)], scatterer()),
    ("d3 ALL-vacuum c=0.5", (1.0, 2.0, 3.0), (3, 4, 5),
     [(V, V)] * 3, scatterer()),
]

quad = Quadrature.level_symmetric(sn_order=4)
W = float(np.sum(quad.weights))
Q_g = np.array([1.0, 0.5])

print(f"{'config':<32} {'G-S':>7} {'Jacobi':>7} {'ratio':>8}  verdict")
print("-" * 74)
for label, extents, cells, bcs, mix in CONFIGS:
    shape = (quad.N, 2, *cells)
    q = np.broadcast_to(
        (Q_g / W)[(slice(None) if False else None), :, *([None] * len(cells))],
        shape,
    ).copy()
    counts = {}
    for sched in ("gauss_seidel", "jacobi"):
        sol = solve_sn_fixed_source(
            {0: mix}, axes(extents, cells, bcs), quad,
            external_source=q, inner_solver="source_iteration",
            inner_schedule=sched, inner_tol=1e-13, max_inner=20000,
        )
        counts[sched] = (sol.history.n_inner, bool(sol.history.converged))
    gs, ja = counts["gauss_seidel"], counts["jacobi"]
    if not (gs[1] and ja[1]):
        print(f"{label:<32} {gs[0]:>7} {ja[0]:>7} {'---':>8}  "
              f"NOT CONVERGED gs={gs[1]} ja={ja[1]}")
        continue
    ratio = gs[0] / ja[0]
    verdict = "G-S faster" if ratio < 0.95 else (
        "G-S SLOWER" if ratio > 1.05 else "~tie")
    print(f"{label:<32} {gs[0]:>7} {ja[0]:>7} {ratio:>8.3f}  {verdict}")
