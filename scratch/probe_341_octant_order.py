"""#341 experiment C — does the OCTANT SWEEP ORDER control the boundary-G-S rate?

Settles: the boundary-G-S splitting folds `Sum_a L_a` of the `d*2^d` directed
octant-coupling edges, where `L_a` is the length of the maximal SUFFIX of the
octant sweep order over which the sign on axis `a` is constant
(derivation in `scratch/issue_341_gs_jacobi_mechanism.md` §2.1).  `L_a` is a
pure property of the ORDER, so permuting `Quadrature.octants` (an in-process
monkeypatch of a `cached_property`; no production change) varies the folded
fraction at fixed physics.

Rows C1/C2 hold `Sum L_a = 7` fixed and vary WHICH axis is folded; C3/C4 drop
`Sum L_a` to its minimum `d`.  Every row re-verifies the fixed point is
splitting-invariant, so the permutation is a rate experiment only.

Usage:  .venv/bin/python -O scratch/probe_341_octant_order.py
"""
from __future__ import annotations

import time

import numpy as np

from orpheus.geometry import BC
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.loss_representation.sweep_schedule import SweepSchedule
from orpheus.sn.solver import _as_sn_mesh
import scratch.probe_341_iteration_spectrum as P

R = BC("reflective")


def permute_octants(quad, key):
    """Return a NEW Quadrature-like view whose cached `octants` is reordered.

    `key(label3) -> sort key`, label3 the 3-sign tuple of the partition entry.
    Mutates the instance's cached_property slot; caller owns the lifetime.
    """
    import copy
    q = copy.copy(quad)
    q.__dict__.pop("octants", None)
    base = quad.octants
    order = sorted(range(len(base)),
                   key=lambda i: key(tuple(int(v) for v in base[i].label)))
    q.__dict__["octants"] = tuple(base[i] for i in order)
    return q


def suffix_runs(group_labels):
    """`L_a` per axis: length of the maximal constant-sign suffix."""
    d = len(group_labels[0])
    out = []
    for a in range(d):
        n = 1
        for k in range(len(group_labels) - 2, -1, -1):
            if group_labels[k][a] == group_labels[-1][a]:
                n += 1
            else:
                break
        out.append(n)
    return out


def measure(extents, cells, mixture, quad, ndim, k=16):
    bcs = [(R, R)] * ndim
    sn_mesh = _as_sn_mesh(P.axes(extents, cells, bcs), quad, {0: mixture})
    sch = SweepSchedule.gauss_seidel(sn_mesh)
    labels = [sch.groups[i].sweeps[0].label.signs for i in range(len(sch.groups))]
    L = suffix_runs(labels)
    rows = sch.lower_inflow_rows(sn_mesh)
    tr = sn_mesh.angular_trace
    n_low = sum(len(v) for v in rows.values())
    n_all = sum(len(tr.inflow_indices_for_face(f)) for f in tr.layout.faces)
    out = {}
    for sched in ("gauss_seidel", "jacobi"):
        G, template, meta = P.build_iteration_operator(
            extents, cells, bcs, mixture, quad, sched,
        )
        vals, _ = P.spectrum(G, k=k)
        out[sched] = P._rho_below_one(vals)
    return dict(labels=labels, L=L, sumL=sum(L),
                frac=n_low / n_all, n_low=n_low, n_all=n_all, **out)


KEYS = {
    # name -> (key function on the 3-sign label, human note)
    "lex-x-major   (BASELINE)": lambda s: (s[0], s[1], s[2]),
    "reversed             C1": lambda s: (-s[0], -s[1], -s[2]),
    "lex-z-major          C2": lambda s: (s[2], s[1], s[0]),
    "antipodal-last       C3": lambda s: (
        {(-1, -1, -1): 7, (1, 1, 1): 6, (-1, -1, 1): 0, (-1, 1, -1): 1,
         (1, -1, -1): 2, (-1, 1, 1): 3, (1, -1, 1): 4, (1, 1, -1): 5}[s],
    ),
}

KEYS_D2 = {
    "lex-x-major   (BASELINE)": lambda s: (s[0], s[1], s[2]),
    # in-plane group order becomes [(+,-), (-,+), (+,+), (-,-)]  -> Sum L = 2
    "antipodal-last       C4": lambda s: (
        {(1, -1): 0, (-1, 1): 1, (1, 1): 2, (-1, -1): 3}[(s[0], s[1])], s[2],
    ),
}


def main():
    quad0 = Quadrature.level_symmetric(sn_order=4)
    mix = P.absorber()
    print("d=3  extents (1,2,3) cells (3,4,5) all-reflective absorber "
          "Sigma_t=(0.8,1.6) LS4 DD")
    print(f"{'order':<26} {'L_a':<12} {'sumL':>5} {'frac':>7} "
          f"{'rho_GS':>10} {'rho_J':>10} {'ln ratio':>9} {'s':>5}")
    print("-" * 96)
    for name, key in KEYS.items():
        t0 = time.perf_counter()
        q = permute_octants(quad0, key)
        r = measure((1.0, 2.0, 3.0), (3, 4, 5), mix, q, 3)
        print(f"{name:<26} {str(r['L']):<12} {r['sumL']:>5} {r['frac']:>7.4f} "
              f"{r['gauss_seidel']:>10.6f} {r['jacobi']:>10.6f} "
              f"{np.log(r['jacobi']) / np.log(r['gauss_seidel']):>9.3f} "
              f"{time.perf_counter() - t0:>5.0f}")
    print()
    print("d=2  extents (1,2) cells (3,4) all-reflective absorber LS4 DD")
    print(f"{'order':<26} {'L_a':<12} {'sumL':>5} {'frac':>7} "
          f"{'rho_GS':>10} {'rho_J':>10} {'ln ratio':>9} {'s':>5}")
    print("-" * 96)
    for name, key in KEYS_D2.items():
        t0 = time.perf_counter()
        q = permute_octants(quad0, key)
        r = measure((1.0, 2.0), (3, 4), mix, q, 2)
        print(f"{name:<26} {str(r['L']):<12} {r['sumL']:>5} {r['frac']:>7.4f} "
              f"{r['gauss_seidel']:>10.6f} {r['jacobi']:>10.6f} "
              f"{np.log(r['jacobi']) / np.log(r['gauss_seidel']):>9.3f} "
              f"{time.perf_counter() - t0:>5.0f}")


if __name__ == "__main__":
    main()
