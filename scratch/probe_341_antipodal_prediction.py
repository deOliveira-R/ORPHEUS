"""#341 PREDICTION D — does an antipodal-last octant order remove the inversion?

Claim under test (`issue_341_gs_jacobi_mechanism.md` §7.1): the amplification that
makes boundary G-S LOSE to Jacobi is carried by the n >= 2 terms of
(I - K_l)^{-1} = sum_n K_l^n, so an octant order whose fold has nilpotency
index 2 (K_l^2 = 0 -- a STAR into the last octant, not a CHAIN) cannot invert
the comparison:  rho_GS <= rho_J on every fixture.

Re-runs the §6 sweep fixtures under both the shipped lexicographic order and the
antipodal-last order.  A single row with rho_GS(antipodal) > rho_J REFUTES the claim.

Usage:  .venv/bin/python -O scratch/probe_341_antipodal_prediction.py
"""
from __future__ import annotations

import itertools
import time

import numpy as np

from orpheus.geometry import BC
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.loss_representation.sweep_schedule import SweepSchedule
from orpheus.sn.solver import _as_sn_mesh
import scratch.probe_341_iteration_spectrum as P
from scratch.probe_341_octant_order import permute_octants, suffix_runs
from scratch.probe_341_predicate_sweep import absorber_1g

R = BC("reflective")


def antipodal_key(d):
    """Order key: last two in-plane octants differ on EVERY axis -> Sum L_a = d."""
    labs = [tuple(s) for s in itertools.product((-1, 1), repeat=d)]
    last, prev = tuple([-1] * d), tuple([1] * d)
    rest = [s for s in labs if s not in (last, prev)]
    pos = {s: i for i, s in enumerate(rest + [prev, last])}
    return lambda s3: (pos[tuple(s3[:d])], s3[d:])


FIXTURES = [
    # (label, extents, cells, sig_t, quad_name)  -- the seven PREDICTED rows first
    ("d2 (1,1)   ext(1,2)  St=0.8  LS4", (1., 2.), (1, 1), 0.8, "LS4"),
    ("d2 (2,2)   ext(6,6)  St=2.0  LS4", (6., 6.), (2, 2), 2.0, "LS4"),
    ("d2 (3,3)   ext(1,1)  St=51.2 LS4", (1., 1.), (3, 3), 51.2, "LS4"),
    ("d3 (2,2,2) ext(1,1,1) St=0.8  LS4", (1., 1., 1.), (2, 2, 2), 0.8, "LS4"),
    ("d3 (3,3,3) ext(1,1,1) St=0.8  LS4", (1., 1., 1.), (3, 3, 3), 0.8, "LS4"),
    ("d3 (3,3,3) ext(1,1,1) St=12.8 LS4", (1., 1., 1.), (3, 3, 3), 12.8, "LS4"),
    ("d3 (5,5,5) ext(1,1,1) St=0.8  LS4", (1., 1., 1.), (5, 5, 5), 0.8, "LS4"),
    # cost rows: lex already wins here
    ("d2 (3,4)   ext(1,2)  St=0.8  LS4", (1., 2.), (3, 4), 0.8, "LS4"),
    ("d2 (10,10) ext(1,1)  St=0.8  LS4", (1., 1.), (10, 10), 0.8, "LS4"),
    ("d3 (2,4,8) ext(1,2,4) St=0.8  LS4", (1., 2., 4.), (2, 4, 8), 0.8, "LS4"),
    ("d3 (3,4,5) ext(1,2,3) St=0.8  LS6", (1., 2., 3.), (3, 4, 5), 0.8, "LS6"),
    # the row the claim was BUILT on (not a test)
    ("d3 (3,4,5) ext(1,2,3) St=0.8  LS4", (1., 2., 3.), (3, 4, 5), 0.8, "LS4"),
]

QUADS = {"LS2": 2, "LS4": 4, "LS6": 6}


def main():
    quads = {k: Quadrature.level_symmetric(sn_order=v) for k, v in QUADS.items()}
    print(f"{'fixture':<36} {'ord':<6} {'L_a':<10} {'rho_GS':>9} {'rho_J':>9} "
          f"{'nGS/nJ':>7}  verdict")
    print("-" * 96)
    fails = []
    for label, extents, cells, st, qn in FIXTURES:
        d = len(cells)
        bcs = [(R, R)] * d
        mix = absorber_1g(st)
        base = quads[qn]
        rho = {}
        for oname, q in (("lex", base),
                         ("anti", permute_octants(base, antipodal_key(d)))):
            t0 = time.perf_counter()
            got = {}
            for sched in ("gauss_seidel", "jacobi"):
                G, _, _ = P.build_iteration_operator(
                    extents, cells, bcs, mix, q, sched,
                )
                vals, _ = P.spectrum(G, k=16)
                got[sched] = P._rho_below_one(vals)
            sn_mesh = _as_sn_mesh(P.axes(extents, cells, bcs), q, {0: mix})
            L = suffix_runs([g.sweeps[0].label.signs
                             for g in SweepSchedule.gauss_seidel(sn_mesh).groups])
            rg, rj = got["gauss_seidel"], got["jacobi"]
            ratio = np.log(rj) / np.log(rg)
            v = "G-S LOSES" if ratio > 1.0 else "G-S wins"
            rho[oname] = (rg, rj, ratio)
            print(f"{label:<36} {oname:<6} {str(tuple(L)):<10} {rg:>9.6f} "
                  f"{rj:>9.6f} {ratio:>7.3f}  {v}  ({time.perf_counter()-t0:.0f}s)")
            if oname == "anti" and rg > rj:
                fails.append((label, rg, rj))
        print()
    print("=" * 96)
    if fails:
        print("PREDICTION D REFUTED on:")
        for lab, rg, rj in fails:
            print(f"   {lab}: rho_GS(anti)={rg:.6f} > rho_J={rj:.6f}")
    else:
        print("PREDICTION D HELD on every row: rho_GS(antipodal) <= rho_J.")


if __name__ == "__main__":
    main()
