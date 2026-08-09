"""#341 — does the boundary-G-S / Jacobi inversion survive c -> 1?

The whole #341 evidence base is a PURE ABSORBER (c = 0), and the one scattering
row in the issue is c = 0.5.  The production regime a user reaches through
`solve_sn` on an all-reflective lattice is c ~ 0.9-0.99, where source iteration's
own scattering mode (rho ~ c) competes with the boundary mode.  This probe fills
that gap: rho_GS vs rho_J across c at both dimensions.

Usage:  .venv/bin/python -O scratch/probe_341_scattering_ratio.py
"""
from __future__ import annotations

import numpy as np

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC
from orpheus.numerics.quadrature import Quadrature
import scratch.probe_341_iteration_spectrum as P

R = BC("reflective")


def mix_1g(sig_t, c):
    st = np.array([float(sig_t)])
    ss = np.array([[c * float(sig_t)]])
    return make_mixture(
        sig_t=st, sig_c=st - ss.sum(axis=0), sig_f=np.zeros(1),
        nu=np.zeros(1), chi=np.zeros(1), sig_s=ss,
    )


def main():
    quad = Quadrature.level_symmetric(sn_order=4)
    print("all-reflective, Sigma_t = 0.8, LS4, DD, shipped lex octant order")
    print(f"{'fixture':<30} {'c':>6} {'rho_GS':>9} {'rho_J':>9} "
          f"{'nGS/nJ':>7}  verdict")
    print("-" * 76)
    for label, extents, cells in [
        ("d2 (3,4) ext(1,2)", (1., 2.), (3, 4)),
        ("d3 (3,4,5) ext(1,2,3)", (1., 2., 3.), (3, 4, 5)),
    ]:
        d = len(cells)
        bcs = [(R, R)] * d
        for c in (0.0, 0.5, 0.9, 0.99, 0.999):
            mix = mix_1g(0.8, c)
            got = {}
            for sched in ("gauss_seidel", "jacobi"):
                G, _, _ = P.build_iteration_operator(
                    extents, cells, bcs, mix, quad, sched,
                )
                vals, _ = P.spectrum(G, k=16)
                got[sched] = P._rho_below_one(vals)
            rg, rj = got["gauss_seidel"], got["jacobi"]
            ratio = np.log(rj) / np.log(rg)
            print(f"{label:<30} {c:>6.3f} {rg:>9.6f} {rj:>9.6f} {ratio:>7.3f}"
                  f"  {'G-S LOSES' if ratio > 1.02 else ('G-S wins' if ratio < 0.98 else 'tie')}")
        print()


if __name__ == "__main__":
    main()
