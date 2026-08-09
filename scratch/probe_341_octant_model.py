"""#341 — the 2^d octant-amplitude MODEL of the reflective SI iteration.

Settles the mechanism analytically on the smallest object that carries it.

Model.  On an all-reflective box every ordinate's boundary coupling is a walk on
the octant hypercube Q_d: reflection at face (a, +/-) flips the sign of the a-th
direction cosine, so octant `s` feeds octant `s XOR a`.  Collapse each octant to
ONE scalar amplitude and give axis `a` a single round-trip gain `g_a`:

    Jacobi   K = sum_a g_a T_a       (T_a = flip-bit-a permutation on {+-1}^d)
    G-S      G = (I - K_l)^{-1} K_u,  K_l + K_u = K

with `K_l` the exact ORPHEUS fold: face f = (a, sign) is reflected after its LAST
outflowing octant group, so the edge (s -> s XOR a) through face (a, s_a) is
implicit iff  sigma(s XOR a) > max{sigma(s') : s'_a = s_a}.

Two facts fall straight out:

  * The T_a commute, so K is diagonalised by the characters of (Z_2)^d and
        spec(K) = { sum_a g_a eps_a : eps in {+-1}^d },   rho_J = sum_a |g_a|.
    ** rho_J is blind to the SIGNS of g_a. **
  * If every g_a > 0 then I-K_l is a regular splitting of I-K with
    K_u <= K elementwise, so Varga's comparison theorem forces rho_GS <= rho_J
    at EVERY dimension and for EVERY order.  An inversion therefore REQUIRES a
    sign-indefinite gain vector.

Physically the sign is the multi-D diamond-difference self-transmission

    coeff_a = (w_a - s - sum_{b!=a} w_b) / (s + sum_b w_b),
    w_a = 2|mu_a| A_a,  s = Sigma_t V

which is negative unless one axis's `w_a` exceeds everything else combined --
so AT MOST ONE axis per ordinate is positive, and >= d-1 are negative.

Usage:  .venv/bin/python -O scratch/probe_341_octant_model.py
"""
from __future__ import annotations

import itertools

import numpy as np


def octant_labels(d):
    return [tuple(s) for s in itertools.product((-1, 1), repeat=d)]


def lex_order(d, major=0):
    """Lexicographic on the signs with axis `major` most significant."""
    labs = octant_labels(d)
    axes = [major] + [a for a in range(d) if a != major]
    return sorted(labs, key=lambda s: tuple(s[a] for a in axes))


def antipodal_last_order(d):
    """An order whose last two octants differ on EVERY axis -> Sum L_a = d."""
    labs = octant_labels(d)
    last = tuple([-1] * d)
    prev = tuple([1] * d)
    rest = [s for s in labs if s not in (last, prev)]
    return rest + [prev, last]


def suffix_runs(order):
    d = len(order[0])
    out = []
    for a in range(d):
        n = 1
        for k in range(len(order) - 2, -1, -1):
            if order[k][a] == order[-1][a]:
                n += 1
            else:
                break
        out.append(n)
    return out


def build(order, g):
    """Return (K, K_l) for the ORPHEUS deferred-reflect fold under `order`."""
    d = len(g)
    pos = {s: i for i, s in enumerate(order)}
    n = len(order)
    K = np.zeros((n, n))
    Kl = np.zeros((n, n))
    # T_f = last position of an octant OUTFLOWING through face (a, sign)
    T = {}
    for a in range(d):
        for sg in (-1, 1):
            T[(a, sg)] = max(pos[s] for s in order if s[a] == sg)
    for s in order:
        for a in range(d):
            t = tuple(-v if b == a else v for b, v in enumerate(s))
            K[pos[t], pos[s]] += g[a]
            if pos[t] > T[(a, s[a])]:
                Kl[pos[t], pos[s]] += g[a]
    return K, Kl


def rates(order, g):
    K, Kl = build(order, g)
    Ku = K - Kl
    G = np.linalg.solve(np.eye(len(K)) - Kl, Ku)
    rho_j = max(abs(np.linalg.eigvals(K)))
    rho_gs = max(abs(np.linalg.eigvals(G)))
    return rho_j, rho_gs, int(round(Kl.astype(bool).sum()))


def main():
    print("A. Jacobi is SIGN-BLIND, boundary G-S is not "
          "(|g| = (0.45, 0.30, 0.20), lex x-major order)")
    print(f"{'signs':<14} {'rho_J':>9} {'rho_GS':>9} {'GS/J sweeps':>12}  verdict")
    print("-" * 62)
    mag = (0.45, 0.30, 0.20)
    for sgn in itertools.product((1, -1), repeat=3):
        g = tuple(m * s for m, s in zip(mag, sgn))
        rj, rg, nl = rates(lex_order(3), g)
        v = ("G-S faster" if rg < rj * 0.999 else
             "G-S SLOWER" if rg > rj * 1.001 else "tie")
        print(f"{str(sgn):<14} {rj:>9.5f} {rg:>9.5f} "
              f"{np.log(rj) / np.log(rg):>12.3f}  {v}")

    print("\nB. dimension x sign-count, lex order, |g_a| = 0.30 each")
    print(f"{'d':>2} {'n_neg':>6} {'sumL':>5} {'rho_J':>9} {'rho_GS':>9} "
          f"{'GS/J sweeps':>12}  verdict")
    print("-" * 62)
    for d in (1, 2, 3, 4):
        order = lex_order(d)
        L = suffix_runs(order)
        for n_neg in range(d + 1):
            g = tuple(-0.30 if a < n_neg else 0.30 for a in range(d))
            rj, rg, _ = rates(order, g)
            v = ("G-S faster" if rg < rj * 0.999 else
                 "G-S SLOWER" if rg > rj * 1.001 else "tie")
            print(f"{d:>2} {n_neg:>6} {sum(L):>5} {rj:>9.5f} {rg:>9.5f} "
                  f"{np.log(rj) / np.log(rg):>12.3f}  {v}")

    print("\nC. the ORDER, at d=3, |g| = 0.30, two axes negative "
          "(the physical DD sign pattern)")
    g = (0.30, -0.30, -0.30)
    print(f"{'order':<22} {'L_a':<12} {'sumL':>5} {'rho_J':>9} {'rho_GS':>9} "
          f"{'GS/J':>8}")
    print("-" * 72)
    for name, order in [
        ("lex x-major", lex_order(3, 0)),
        ("lex z-major", lex_order(3, 2)),
        ("antipodal-last", antipodal_last_order(3)),
    ]:
        rj, rg, nl = rates(order, g)
        print(f"{name:<22} {str(suffix_runs(order)):<12} "
              f"{sum(suffix_runs(order)):>5} {rj:>9.5f} {rg:>9.5f} "
              f"{np.log(rj) / np.log(rg):>8.3f}")

    print("\nD. all-POSITIVE gains: Varga's comparison theorem must hold "
          "(rho_GS <= rho_J) at every d and every order")
    bad = 0
    rng = np.random.default_rng(0)
    for d in (2, 3, 4):
        for _ in range(200):
            g = tuple(rng.uniform(0.02, 0.9 / d, size=d))
            for order in (lex_order(d), antipodal_last_order(d)):
                rj, rg, _ = rates(order, g)
                if rg > rj + 1e-12:
                    bad += 1
    print(f"    violations over 1200 random positive-gain cases: {bad}")

    print("\nE. mixed-sign gains: how often does G-S LOSE? (random |g|, d=2 vs 3)")
    for d in (2, 3, 4):
        n_lose = {k: 0 for k in range(d + 1)}
        n_tot = {k: 0 for k in range(d + 1)}
        for _ in range(400):
            mags = rng.uniform(0.02, 0.9 / d, size=d)
            for sgn in itertools.product((1, -1), repeat=d):
                g = tuple(m * s for m, s in zip(mags, sgn))
                k = sum(1 for s in sgn if s < 0)
                rj, rg, _ = rates(lex_order(d), g)
                n_tot[k] += 1
                if rg > rj * 1.001:
                    n_lose[k] += 1
        frac = {k: n_lose[k] / max(n_tot[k], 1) for k in n_tot}
        print(f"    d={d}: P(G-S slower | n_neg) = "
              + "  ".join(f"{k}:{frac[k]:.2f}" for k in sorted(frac)))


if __name__ == "__main__":
    main()
