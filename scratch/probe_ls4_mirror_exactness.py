"""Did #337's moment-matched LS seed break the BIT-EXACT mirror structure of S4?

The d3 pure-absorber bias lives on exactly the 8 |mu_x|-dominant ordinates.
A reflective BC needs an EXACT involution on each axis; ERR-073/#326 is the
established defect class (a match loop that returns a many-to-one relation
rather than a permutation).  Check, bit-exactly:
  1. is each cosine set closed under negation, node-for-node?
  2. is every ordinate ON the unit sphere?
  3. is the x-mirror ordinate permutation an involutive BIJECTION?
  4. how far did the S4 cosines actually move at #337?
"""
from __future__ import annotations

import numpy as np

from orpheus.numerics.quadrature import Quadrature

quad = Quadrature.level_symmetric(sn_order=4)
mx = np.asarray(quad.mu_x)
my = np.asarray(quad.mu_y)
mz = np.asarray(quad.mu_z)

print("=== 1. distinct |cosine| values, full precision ===")
vals = np.unique(np.abs(np.concatenate([mx, my, mz])))
for v in vals:
    print(f"   {v!r}")

print("\n=== 2. unit-sphere defect per ordinate (bit-exact) ===")
norm = mx * mx + my * my + mz * mz
print(f"   max |sum mu^2 - 1| = {np.max(np.abs(norm - 1.0)):.6e}")
print(f"   exactly 1.0 for {int(np.sum(norm == 1.0))} / {quad.N} ordinates")

print("\n=== 3. per-axis negation closure (bit-exact set equality) ===")
for name, arr in (("mu_x", mx), ("mu_y", my), ("mu_z", mz)):
    ok = np.array_equal(np.sort(arr), np.sort(-arr))
    print(f"   {name}: set closed under negation, bit-exactly? {ok}")

print("\n=== 4. the x-mirror as a PERMUTATION (the reflective-BC map) ===")
# partner of n = the ordinate with (-mu_x, +mu_y, +mu_z), matched bit-exactly
partner = np.full(quad.N, -1, dtype=int)
for n in range(quad.N):
    hit = np.flatnonzero(
        (mx == -mx[n]) & (my == my[n]) & (mz == mz[n]),
    )
    if hit.size == 1:
        partner[n] = hit[0]
    else:
        print(f"   ordinate {n}: {hit.size} bit-exact partners {hit.tolist()}")
unmatched = int(np.sum(partner < 0))
print(f"   unmatched (no unique bit-exact partner): {unmatched}")
if unmatched == 0:
    inj = len(set(partner.tolist())) == quad.N
    invol = bool(np.array_equal(partner[partner], np.arange(quad.N)))
    print(f"   injective (a bijection)? {inj}")
    print(f"   involutive (P o P = id)? {invol}")

print("\n=== 5. what the production reflection map says ===")
try:
    from orpheus.geometry.transformation import reflection
    perm = quad.ordinate_permutation(reflection(axis=0, dim=3))
    p = np.asarray(perm)
    print(f"   ordinate_permutation(reflection x) = {p.tolist()}")
    print(f"   bijection? {len(set(p.tolist())) == quad.N}")
    print(f"   involutive? {bool(np.array_equal(p[p], np.arange(quad.N)))}")
    print(f"   agrees with the bit-exact match? "
          f"{bool(np.array_equal(p, partner)) if unmatched == 0 else 'n/a'}")
except Exception as exc:                                  # noqa: BLE001
    print(f"   (could not build the production map: {exc!r})")

print("\n=== 6. weights: are mirror partners equal-weighted, bit-exactly? ===")
w = np.asarray(quad.weights)
if unmatched == 0:
    print(f"   max |w[n] - w[partner[n]]| = {np.max(np.abs(w - w[partner])):.6e}")
