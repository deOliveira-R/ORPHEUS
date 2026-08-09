"""#341 — is `ndim` the discriminating variable, or a proxy?

Sweeps mesh / optical thickness / quadrature at BOTH dimensions on the
zero-leakage (all-reflective) fixture and reports rho_GS vs rho_J from the
verified spectral instrument (`probe_341_iteration_spectrum`).  Two falsifiers:

  * ANY d=2 all-reflective fixture with rho_GS > rho_J  ->  `ndim` refuted.
  * ANY d=3 all-reflective fixture with rho_GS < rho_J  ->  `ndim` refuted.

Row 1 is the CONTROL: the 2-group d=3 (3,4,5) point, which must reproduce
rho_GS = 0.98552 / rho_J = 0.97541 (the diagnosis's fitted 0.985348 / 0.975014).
Row 2 repeats it in ONE group to license the 1-group economy used below (a pure
absorber decouples the groups exactly, so the 2-group spectrum is the union).

Usage:  .venv/bin/python -O scratch/probe_341_predicate_sweep.py
"""
from __future__ import annotations

import time

import numpy as np

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC
from orpheus.numerics.quadrature import Quadrature
import scratch.probe_341_iteration_spectrum as P

R = BC("reflective")


def absorber_1g(sig_t):
    st = np.array([float(sig_t)])
    return make_mixture(
        sig_t=st, sig_c=st.copy(), sig_f=np.zeros(1), nu=np.zeros(1),
        chi=np.zeros(1), sig_s=np.zeros((1, 1)),
    )


def undamped_diag(extents, cells, sig_t, quad, ndim):
    """Report the DD single-cell transmission spectrum stats for this fixture.

    Sigma = (2/D) 1 w^T - I with w_a = 2|mu_a| A_a, s = Sigma_t V:
    one eigenvalue 1 - 2s/D (damped) and (d-1) eigenvalues exactly -1.
    Returns (max |1-2s/D| over ordinates, d-1).
    """
    h = np.array([e / n for e, n in zip(extents, cells)], dtype=float)
    Vc = float(np.prod(h))
    A = np.array([Vc / hh for hh in h])
    mu = np.abs(np.asarray(quad.nodes)[:, :ndim])
    w = 2.0 * mu * A[None, :]
    s = sig_t * Vc
    D = s + w.sum(axis=1)
    return float(np.max(np.abs(1.0 - 2.0 * s / D))), ndim - 1


def run(label, extents, cells, sig_t, quad, mixture=None, k=16):
    ndim = len(cells)
    bcs = [(R, R)] * ndim
    mix = mixture if mixture is not None else absorber_1g(sig_t)
    out = {}
    ndof = None
    for sched in ("gauss_seidel", "jacobi"):
        G, template, meta = P.build_iteration_operator(
            extents, cells, bcs, mix, quad, sched,
        )
        ndof = meta["n_dof"]
        vals, _ = P.spectrum(G, k=k)
        out[sched] = P._rho_below_one(vals)
    smooth, n_undamped = undamped_diag(extents, cells, sig_t, quad, ndim)
    rg, rj = out["gauss_seidel"], out["jacobi"]
    ratio = np.log(rj) / np.log(rg)
    print(f"{label:<34} {ndof:>6} {rg:>9.6f} {rj:>9.6f} {ratio:>7.3f} "
          f"{n_undamped:>4} {smooth:>8.4f}  "
          f"{'G-S LOSES' if ratio > 1.02 else ('G-S wins' if ratio < 0.98 else 'tie')}")
    return rg, rj, ratio


def main():
    ls4 = Quadrature.level_symmetric(sn_order=4)
    ls2 = Quadrature.level_symmetric(sn_order=2)
    ls6 = Quadrature.level_symmetric(sn_order=6)
    print(f"{'fixture':<34} {'ndof':>6} {'rho_GS':>9} {'rho_J':>9} "
          f"{'n_GS/nJ':>7} {'d-1':>4} {'|1-2s/D|':>8}")
    print("-" * 96)

    t0 = time.perf_counter()
    print("--- CONTROLS")
    run("d3 (3,4,5) 2g abs LS4 [CTRL]", (1., 2., 3.), (3, 4, 5), 0.8, ls4,
        mixture=P.absorber())
    run("d3 (3,4,5) 1g St=0.8 LS4", (1., 2., 3.), (3, 4, 5), 0.8, ls4)
    run("d2 (3,4)   1g St=0.8 LS4", (1., 2.), (3, 4), 0.8, ls4)

    print("--- d=2: mesh")
    for cells, ext in [((1, 1), (1., 2.)), ((2, 2), (1., 2.)),
                       ((3, 3), (1., 1.)), ((6, 6), (1., 1.)),
                       ((10, 10), (1., 1.)), ((4, 16), (1., 4.))]:
        run(f"d2 {cells} ext{ext} St=0.8 LS4", ext, cells, 0.8, ls4)

    print("--- d=2: optical thickness (cube 3x3, extent 1)")
    for st in (0.05, 0.2, 0.8, 3.2, 12.8, 51.2):
        run(f"d2 (3,3) ext(1,1) St={st} LS4", (1., 1.), (3, 3), st, ls4)

    print("--- d=2: COARSE + THICK (both DD axes strongly negative)")
    for ext, cells, st in [((6., 6.), (3, 3), 0.8), ((12., 12.), (3, 3), 0.8),
                           ((6., 6.), (2, 2), 2.0), ((20., 20.), (4, 4), 1.0)]:
        run(f"d2 {cells} ext{ext} St={st} LS4", ext, cells, st, ls4)

    print("--- d=2: quadrature")
    for nm, q in (("LS2", ls2), ("LS4", ls4), ("LS6", ls6)):
        run(f"d2 (3,4) ext(1,2) St=0.8 {nm}", (1., 2.), (3, 4), 0.8, q)

    print("--- d=3: mesh")
    for cells, ext in [((1, 1, 1), (1., 2., 3.)), ((2, 2, 2), (1., 1., 1.)),
                       ((3, 3, 3), (1., 1., 1.)), ((2, 4, 8), (1., 2., 4.)),
                       ((5, 5, 5), (1., 1., 1.))]:
        run(f"d3 {cells} ext{ext} St=0.8 LS4", ext, cells, 0.8, ls4)

    print("--- d=3: optical thickness (cube 3x3x3, extent 1)")
    for st in (0.05, 0.2, 0.8, 3.2, 12.8):
        run(f"d3 (3,3,3) ext(1,1,1) St={st} LS4", (1., 1., 1.), (3, 3, 3),
            st, ls4)

    print("--- d=3: quadrature")
    for nm, q in (("LS2", ls2), ("LS4", ls4), ("LS6", ls6)):
        run(f"d3 (3,4,5) ext(1,2,3) St=0.8 {nm}", (1., 2., 3.), (3, 4, 5),
            0.8, q)

    print(f"\ntotal {time.perf_counter() - t0:.0f} s")


if __name__ == "__main__":
    main()
