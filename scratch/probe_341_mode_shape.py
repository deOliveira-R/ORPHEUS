"""#341 — WHICH mode is the slow one, for each splitting, at d=2 and d=3?

Extracts the dominant (|lambda| < 1) eigenvector of the SI iteration matrix and
reports where its mass lives: bulk vs trace, per face, per ordinate class, and
the spatial sign pattern on the heaviest face (a sawtooth shows up as a sign
alternation with near-zero mean).

The point of the probe: the issue asks whether boundary G-S and Jacobi are
racing the SAME slow mode.  If they are not, the "flip" is two different
comparisons, not one quantity changing sign.

Usage:  .venv/bin/python -O scratch/probe_341_mode_shape.py
"""
from __future__ import annotations

import numpy as np

from orpheus.geometry import BC
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import _as_sn_mesh
import scratch.probe_341_iteration_spectrum as P

R = BC("reflective")


def describe(label, extents, cells, mix, quad, sched):
    ndim = len(cells)
    bcs = [(R, R)] * ndim
    sn_mesh = _as_sn_mesh(P.axes(extents, cells, bcs), quad, {0: mix})
    G, template, meta = P.build_iteration_operator(
        extents, cells, bcs, mix, quad, sched,
    )
    vals, vecs = P.spectrum(G, k=16)
    idx = next(i for i in range(len(vals)) if abs(vals[i]) < 1 - 1e-8)
    lam, v = vals[idx], vecs[:, idx]
    v = v / np.linalg.norm(v)
    n_int = template.interior.values.size
    bulk, trace = v[:n_int], v[n_int:]
    print(f"\n### {label} / {sched}: lambda = {abs(lam):.6f} "
          f"exp(i {np.degrees(np.angle(lam)):+.2f} deg)   "
          f"windowed={meta['windowed']}  ndof={meta['n_dof']}")
    print(f"    mass: bulk {np.linalg.norm(bulk):.4f}  "
          f"trace {np.linalg.norm(trace):.4f}")

    tr = sn_mesh.angular_trace
    slots = tr.layout.faces
    views = {name: trace[s.offset:s.offset + s.flat_size].reshape(s.shape)
             for name, s in slots.items()}
    masses = {name: float(np.linalg.norm(v)) for name, v in views.items()}
    tot = float(np.linalg.norm(trace))
    print("    per-face trace mass fraction: " + "  ".join(
        f"{f}={m / tot:.3f}" for f, m in
        sorted(masses.items(), key=lambda kv: -kv[1])))
    heavy = max(masses, key=masses.get)
    arr = views[heavy]                            # (N, ng, *face_cells)
    nodes = np.abs(np.asarray(quad.nodes))
    m_ord = np.linalg.norm(arr.reshape(arr.shape[0], -1), axis=1)
    top = np.argsort(-m_ord)[:3]
    print(f"    heaviest face {heavy} shape {arr.shape}; top ordinates "
          f"(|mu_x|,|mu_y|,|mu_z|)=mass: " + "  ".join(
              f"({nodes[i][0]:.3f},{nodes[i][1]:.3f},{nodes[i][2]:.3f})"
              f"={m_ord[i]:.3f}" for i in top))
    k = int(top[0])
    g = int(np.argmax(np.linalg.norm(arr[k].reshape(arr.shape[1], -1), axis=1)))
    sub = arr[k, g]
    ref = sub.flat[int(np.argmax(np.abs(sub)))]
    real = np.real(sub / ref)
    print(f"    spatial pattern, ordinate {k} group {g} on {heavy} "
          f"(normalised): mean={real.mean():+.4f}  "
          f"min={real.min():+.3f} max={real.max():+.3f}")
    for line in np.array2string(real, precision=3, max_line_width=100,
                                suppress_small=True).splitlines():
        print("      " + line)
    return lam


def main():
    quad = Quadrature.level_symmetric(sn_order=4)
    mix = P.absorber()
    for label, extents, cells in [
        ("d2 (3,4) ext(1,2)", (1.0, 2.0), (3, 4)),
        ("d3 (3,4,5) ext(1,2,3)", (1.0, 2.0, 3.0), (3, 4, 5)),
    ]:
        for sched in ("gauss_seidel", "jacobi"):
            describe(label, extents, cells, mix, quad, sched)


if __name__ == "__main__":
    main()
