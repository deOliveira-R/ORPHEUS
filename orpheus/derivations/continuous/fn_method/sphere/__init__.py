"""Sphere F_N benchmarks — placeholder.

Will house :func:`solve_critical_sphere_1g` (Kaper-Lindeman-Leaf 1974
NSE 54, 94 — same paper as slab) and :func:`solve_critical_sphere_2g`
(Siewert-Thomas 1986 NSE 94, 264).

Targets:
  * `Ua-1-0-SP` (problem 14) — :math:`r_c = 2.4248249802` mfp.
  * `PUb-1-0-SP` (problem 8) — flux ratios at 4 sample points.
  * `PU-2-0-SP` (problem 46) — 2G bare Pu sphere :math:`r_c = 1.15513` mfp.

The sphere F_N reference solver is the natural cross-check for the
existing **Variant α** Green's-function family
(:mod:`orpheus.derivations.continuous.peierls_greens_function`), which
already cross-checks cylinder against Sood/Westfall-Metcalf at 8.5e-6.
Once this is populated, the Variant α sphere prototype will have an
independent published-method anchor that does NOT inherit any of
ORPHEUS's Bickley-Naylor or rank-1/rank-2 closure machinery.

TODO: populate after literature-researcher delivers KLL 1974 +
Siewert-Thomas 1986.
"""
from __future__ import annotations
