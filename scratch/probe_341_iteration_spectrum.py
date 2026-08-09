"""#341 — the SI iteration operator's SPECTRUM for the two boundary splittings.

Settles: is the d=2-vs-d=3 G-S/Jacobi sweep-count inversion a genuine
inversion of the two iteration operators' spectral radii, and if so what
does the dominant eigenvector look like in each case?

The SI driver (``orpheus.numerics.iteration.SourceIteration.solve``) is

    psi <- A_inv.apply(q + sum_i g_i.apply(psi), initial_guess=psi_prev)

so with q = 0 the iteration matrix is  G = M^{-1} N  with

    M = base_implicit  (Jacobi: L+C ; boundary-G-S: (L+C) - B_lower)
    N = sum of gains   (Jacobi: S + B ; boundary-G-S: S + B_upper)

built through the PRODUCTION builders (`build_within_group_system`,
`_select_si_splitting`, `_maybe_window`) so the operator measured IS the
operator the solver runs.

Usage:
    .venv/bin/python scratch/probe_341_iteration_spectrum.py [case ...]
"""
from __future__ import annotations

import sys
import time

import numpy as np
import scipy.sparse.linalg as spla

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.coupled_system import build_within_group_system
from orpheus.sn.solver import (
    SNSolver,
    _as_sn_mesh,
    _maybe_window,
    _select_si_splitting,
    _unwindowed_cold_start,
    _windowed_cold_start,
)
from orpheus.transport.mesh.axis import AxisMesh

V, R = BC("vacuum"), BC("reflective")


def axes(extents, cells, bcs):
    return tuple(
        AxisMesh(edges=np.linspace(0.0, e, n + 1), bc_low=lo, bc_high=hi)
        for e, n, (lo, hi) in zip(extents, cells, bcs)
    )


def absorber(sig_t=(0.8, 1.6)):
    sig_t = np.asarray(sig_t, dtype=float)
    return make_mixture(
        sig_t=sig_t, sig_c=sig_t.copy(),
        sig_f=np.zeros(2), nu=np.zeros(2),
        chi=np.zeros(2), sig_s=np.zeros((2, 2)),
    )


def scatterer(c=0.5, sig_t=(0.8, 1.6)):
    sig_t = np.asarray(sig_t, dtype=float)
    s = np.zeros((2, 2))
    s[0, 0] = c * sig_t[0] * 0.7
    s[1, 0] = c * sig_t[0] * 0.3
    s[1, 1] = c * sig_t[1]
    return make_mixture(
        sig_t=sig_t, sig_c=sig_t - s.sum(axis=0), sig_f=np.zeros(2),
        nu=np.zeros(2), chi=np.zeros(2), sig_s=s,
    )


def build_iteration_operator(extents, cells, bcs, mixture, quad, schedule):
    """Return (LinearOperator G, template, windowed, meta) for one splitting."""
    sn_mesh = _as_sn_mesh(axes(extents, cells, bcs), quad, {0: mixture})
    solver = SNSolver(sn_mesh, inner_solver="source_iteration")
    system = build_within_group_system(
        sn_mesh, solver.mat_xs, scattering_op=solver.scattering_op,
    )
    S, B = system.explicit_gains
    base_implicit, gains = _select_si_splitting(
        system.implicit_operator, S, B, sn_mesh, schedule,
    )
    step, windowed = _maybe_window(base_implicit.inverse(), S, sn_mesh)
    template = (
        _windowed_cold_start(solver.scattering_op, sn_mesh, history_depth=0)
        if windowed else
        _unwindowed_cold_start(sn_mesh, history_depth=0)
    )
    n = template.to_flat().size

    calls = {"n": 0}

    def matvec(flat):
        calls["n"] += 1
        x = type(template).from_flat(np.asarray(flat, dtype=float), template)
        rhs = None
        for g in gains:
            gx = g.apply(x)
            rhs = gx if rhs is None else rhs + gx
        out = step.apply(rhs, initial_guess=x)
        return out.to_flat()

    meta = dict(
        n_dof=n, windowed=windowed, sn_mesh=sn_mesh, calls=calls,
        n_groups_gs=len(
            _gs_groups(sn_mesh) if schedule == "gauss_seidel" else ()
        ),
    )
    return spla.LinearOperator((n, n), matvec=matvec, dtype=float), template, meta


def _gs_groups(sn_mesh):
    from orpheus.sn.loss_representation.sweep_schedule import SweepSchedule
    return SweepSchedule.gauss_seidel(sn_mesh).groups


def spectrum(G, k=6, tol=1e-10, maxiter=None):
    n = G.shape[0]
    kk = min(k, n - 2)
    vals, vecs = spla.eigs(
        G, k=kk, which="LM", tol=tol, maxiter=maxiter or 10 * n,
        ncv=min(n - 1, max(4 * kk, 30)),
    )
    order = np.argsort(-np.abs(vals))
    return vals[order], vecs[:, order]


def report_modes(vals, vecs, template, n_show=10):
    """Print |lambda|, arg, and the bulk/trace mass split of each mode."""
    n_int = template.interior.values.size
    print(f"    {'#':>2} {'|lambda|':>12} {'arg deg':>9} "
          f"{'bulk frac':>10} {'trace frac':>11}")
    for i in range(min(n_show, len(vals))):
        v = vecs[:, i]
        nv = np.linalg.norm(v)
        fb = np.linalg.norm(v[:n_int]) / nv
        ft = np.linalg.norm(v[n_int:]) / nv
        print(f"    {i:>2} {abs(vals[i]):>12.8f} "
              f"{np.degrees(np.angle(vals[i])):>9.2f} "
              f"{fb:>10.4f} {ft:>11.4f}")


CASES = {
    # label: (extents, cells, bcs, mixture)
    "d2-refl-abs": ((1.0, 2.0), (3, 4), [(R, R)] * 2, absorber()),
    "d3-refl-abs": ((1.0, 2.0, 3.0), (3, 4, 5), [(R, R)] * 3, absorber()),
    "d2-refl-c05": ((1.0, 2.0), (3, 4), [(R, R)] * 2, scatterer()),
    "d3-refl-c05": ((1.0, 2.0, 3.0), (3, 4, 5), [(R, R)] * 3, scatterer()),
    "d3-xvac-c05": ((1.0, 2.0, 3.0), (3, 4, 5),
                    [(V, V), (R, R), (R, R)], scatterer()),
}


def main(argv):
    quad = Quadrature.level_symmetric(sn_order=4)
    wanted = argv[1:] or ["d2-refl-abs", "d3-refl-abs"]
    store = {}
    for label in wanted:
        extents, cells, bcs, mix = CASES[label]
        for sched in ("gauss_seidel", "jacobi"):
            t0 = time.perf_counter()
            G, template, meta = build_iteration_operator(
                extents, cells, bcs, mix, quad, sched,
            )
            # linearity control: G is linear (G(0)=0, G(2x)=2G(x))
            rng = np.random.default_rng(0)
            xr = rng.standard_normal(meta["n_dof"])
            lin = np.max(np.abs(G @ (2 * xr) - 2 * (G @ xr)))
            zero = np.max(np.abs(G @ np.zeros(meta["n_dof"])))
            vals, vecs = spectrum(G, k=12)
            dt = time.perf_counter() - t0
            print(f"\n=== {label} / {sched} : ndof={meta['n_dof']} "
                  f"windowed={meta['windowed']} "
                  f"gs_groups={meta['n_groups_gs']} "
                  f"[linearity {lin:.2e}, G(0) {zero:.2e}] {dt:.1f}s")
            report_modes(vals, vecs, template)
            store[(label, sched)] = (vals, vecs, template, meta)
    print()
    for label in wanted:
        try:
            vg = store[(label, "gauss_seidel")][0]
            vj = store[(label, "jacobi")][0]
        except KeyError:
            continue
        rg = _rho_below_one(vg)
        rj = _rho_below_one(vj)
        print(f"{label}: rho_GS={rg:.6f}  rho_J={rj:.6f}  "
              f"(excluding |lambda|>=1-1e-8)   "
              f"predicted n_GS/n_J = {np.log(rj) / np.log(rg):.3f}")
    return store


def _rho_below_one(vals, cut=1.0 - 1e-8):
    below = [abs(v) for v in vals if abs(v) < cut]
    return max(below) if below else float("nan")


if __name__ == "__main__":
    main(sys.argv)
