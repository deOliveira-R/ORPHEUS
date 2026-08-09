"""#341 — rho_GS over EVERY achievable fold pattern (L_x, L_y, L_z).

The boundary-G-S fold is fully described by the per-axis suffix run
`L_a` of the octant sweep order (see `issue_341_gs_jacobi_mechanism.md` §2.1):
axis `a` has `L_a` of its `2^(d-1)` inflow octants reading a FRESH reflection.
Enumerating all `(2^d)!` orders gives a small set of achievable `(L_a)` tuples;
this probe measures rho_GS for one representative of each, so the fold pattern
-> rate map is complete rather than sampled.

Usage:  .venv/bin/python -O scratch/probe_341_fold_pattern.py [d2|d3]
"""
from __future__ import annotations

import itertools
import sys
import time

import numpy as np

from orpheus.geometry import BC
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.loss_representation.sweep_schedule import SweepSchedule
from orpheus.sn.solver import _as_sn_mesh
import scratch.probe_341_iteration_spectrum as P
from scratch.probe_341_octant_order import permute_octants, suffix_runs

R = BC("reflective")


def achievable_patterns(d):
    """{(L_a tuple): representative order (list of sign tuples)}."""
    labs = [tuple(s) for s in itertools.product((-1, 1), repeat=d)]
    out = {}
    for order in itertools.permutations(labs):
        key = tuple(suffix_runs(list(order)))
        out.setdefault(key, list(order))
    return out


def key_from_order(order, d):
    """A sort key over the 3-sign quadrature labels realising `order` in-plane."""
    pos = {s: i for i, s in enumerate(order)}
    return lambda s3: (pos[tuple(s3[:d])], s3[d:])


def main(argv):
    which = argv[1] if len(argv) > 1 else "d3"
    quad0 = Quadrature.level_symmetric(sn_order=4)
    mix = P.absorber()
    if which == "d2":
        d, extents, cells = 2, (1.0, 2.0), (3, 4)
    else:
        d, extents, cells = 3, (1.0, 2.0, 3.0), (3, 4, 5)
    pats = achievable_patterns(d)
    print(f"{which}: extents {extents} cells {cells} all-reflective absorber "
          f"Sigma_t=(0.8,1.6) LS4 DD; {len(pats)} achievable fold patterns")
    print(f"{'L_a':<14} {'sumL':>5} {'rep order':<44} {'rho_GS':>9} "
          f"{'rho_J':>9} {'nGS/nJ':>8} {'s':>4}")
    print("-" * 100)
    rows = []
    for key in sorted(pats, key=lambda k: (sum(k), k)):
        order = pats[key]
        t0 = time.perf_counter()
        q = permute_octants(quad0, key_from_order(order, d))
        bcs = [(R, R)] * d
        rho = {}
        for sched in ("gauss_seidel", "jacobi"):
            G, _, _ = P.build_iteration_operator(
                extents, cells, bcs, mix, q, sched,
            )
            vals, _ = P.spectrum(G, k=16)
            rho[sched] = P._rho_below_one(vals)
        # confirm the schedule really realises the intended pattern
        sn_mesh = _as_sn_mesh(P.axes(extents, cells, bcs), q, {0: mix})
        sch = SweepSchedule.gauss_seidel(sn_mesh)
        got = tuple(suffix_runs([g.sweeps[0].label.signs for g in sch.groups]))
        assert got == key, f"pattern mismatch {got} != {key}"
        rg, rj = rho["gauss_seidel"], rho["jacobi"]
        ratio = np.log(rj) / np.log(rg)
        rows.append((key, rg, rj, ratio))
        short = " ".join("".join("+" if v > 0 else "-" for v in s)
                         for s in order)
        print(f"{str(key):<14} {sum(key):>5} {short:<44} {rg:>9.6f} "
              f"{rj:>9.6f} {ratio:>8.3f} {time.perf_counter() - t0:>4.0f}")
    best = min(rows, key=lambda r: r[1])
    worst = max(rows, key=lambda r: r[1])
    print(f"\nbest  fold pattern {best[0]}: rho_GS={best[1]:.6f} "
          f"(nGS/nJ={best[3]:.3f})")
    print(f"worst fold pattern {worst[0]}: rho_GS={worst[1]:.6f} "
          f"(nGS/nJ={worst[3]:.3f})")


if __name__ == "__main__":
    main(sys.argv)
