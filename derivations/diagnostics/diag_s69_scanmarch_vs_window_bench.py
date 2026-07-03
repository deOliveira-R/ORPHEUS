"""S6.9 / S5.3 — the window-fate measurement (#222).

Does ``ScanMarch`` beat ``MovingFrontierWindow`` end-to-end?  Times the
SWEEP and the MATVEC (loss_action) for all three 2-D representations
(window / scan-march / full-field oracle) across a size × angular-order ×
group grid, plus tracemalloc peak for the sweep, plus one end-to-end
``solve_sn_fixed_source`` (window default vs ScanMarch forced).

Run:  .venv/bin/python -O derivations/diagnostics/diag_s69_scanmarch_vs_window_bench.py
"""
from __future__ import annotations

import time
import tracemalloc
from contextlib import contextmanager
from unittest.mock import patch

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.loss_representation import (
    FullFieldWavefront,
    MovingFrontierWindow,
    ScanMarch,
)
from orpheus.sn.operator import StreamingOperator
from orpheus.sn.solver import solve_sn_fixed_source
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField

REPS = {
    "window": MovingFrontierWindow,
    "scanmarch": ScanMarch,
    "fullfield": FullFieldWavefront,
}


def build(nx, ny, lvl, ng):
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 2.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_xmin=BC("reflective"), bc_xmax=BC("reflective"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    sn = SNMesh(mesh, Quadrature.level_symmetric(lvl), {0: get_mixture("A", f"{ng}g")})
    rng = np.random.default_rng(7)
    N = sn.quad.N
    sig_t = rng.uniform(0.3, 3.0, size=(ng, nx, ny))
    Q = rng.uniform(0.0, 2.0, size=(N, ng, nx, ny))
    return sn, sig_t, Q


def median_time(fn, n=9):
    ts = []
    for _ in range(n):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    return float(np.median(ts))


def peak_mem(fn):
    tracemalloc.start()
    fn()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return peak


def bench_kernels():
    print(f"{'config':<22}{'rep':<11}{'sweep_ms':>10}{'matvec_ms':>11}{'sweep_peakMB':>14}")
    grid = [
        (24, 24, 4, 2),
        (48, 48, 4, 2),
        (96, 96, 4, 2),
        (48, 48, 8, 2),
        (48, 48, 16, 2),
        (48, 48, 8, 4),
    ]
    summary = {}
    for nx, ny, lvl, ng in grid:
        sn, sig_t, Q = build(nx, ny, lvl, ng)
        L = StreamingOperator(sn)          # pure σ-free streaming (#257 S8b)
        psi = TimedFullField.zeros(bulk=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn)
        psi.bulk.values[...] = np.random.default_rng(3).standard_normal(
            psi.bulk.values.shape)
        rows = {}
        for name, cls in REPS.items():
            rep = cls(sn)
            bf = AngularBoundaryFlux.zeros_on(sn)
            t_sweep = median_time(lambda: rep.sweep(Q, sig_t, bf))
            t_mv = median_time(lambda: rep.loss_action(L, psi))
            m_sweep = peak_mem(lambda: rep.sweep(Q, sig_t, AngularBoundaryFlux.zeros_on(sn)))
            rows[name] = (t_sweep, t_mv, m_sweep)
            print(f"{f'{nx}x{ny} LS{lvl} {ng}g':<22}{name:<11}"
                  f"{t_sweep*1e3:>10.2f}{t_mv*1e3:>11.2f}{m_sweep/1e6:>14.2f}")
        summary[(nx, ny, lvl, ng)] = rows
        w, s = rows["window"], rows["scanmarch"]
        print(f"{'':<22}{'sm/win':<11}{s[0]/w[0]:>10.2f}{s[1]/w[1]:>11.2f}"
              f"{s[2]/w[2]:>14.2f}")
    return summary


@contextmanager
def scanmarch_forced():
    from orpheus.sn import loss_representation as lr
    real = lr.default_for

    def forced(mesh):
        rep = real(mesh)
        return lr.ScanMarch(mesh) if isinstance(rep, lr.MovingFrontierWindow) else rep

    with patch.object(lr, "default_for", forced):
        yield


def bench_end_to_end():
    fuel = get_mixture("A", "2g")
    mod = get_mixture("B", "2g")
    nx, ny = 48, 48
    mat = np.zeros((nx, ny), dtype=int)
    mat[: nx // 2, :] = 2
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 12.0, nx + 1),
        edges_y=np.linspace(0.0, 12.0, ny + 1),
        mat_map=mat,
        bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    quad = Quadrature.level_symmetric(sn_order=8)
    sum_w = float(quad.weights.sum())
    q = np.full((quad.N, 2, nx, ny), 1.0 / sum_w)
    mats = {2: fuel, 0: mod}

    def run():
        return solve_sn_fixed_source(
            mats, mesh, quad, q, scattering_order=1,
            inner_solver="source_iteration", max_inner=2000, inner_tol=1e-10,
        )

    # warm both paths once (caches: DAG per shape, xs prep)
    run()
    t_win = median_time(run, n=3)
    with scanmarch_forced():
        run()
        t_sm = median_time(run, n=3)
    print(f"\nend-to-end solve_sn_fixed_source 48x48 LS8 2G het (median of 3):")
    print(f"  window    : {t_win:.3f} s")
    print(f"  scanmarch : {t_sm:.3f} s   (sm/win = {t_sm/t_win:.2f})")


if __name__ == "__main__":
    bench_kernels()
    bench_end_to_end()
