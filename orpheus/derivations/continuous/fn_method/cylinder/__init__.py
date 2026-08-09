"""Cylinder F_N benchmarks — placeholder.

Will house :func:`solve_critical_cylinder_1g` (Westfall-Metcalf 1973
NSE 52, 1 — singular eigenfunction, not strictly F_N for cylinder)
and :func:`solve_critical_cylinder_1g_anisotropic` (Sanchez-Ganapol
1983 NSE 64, 61 — F_N method for cylinder with linear anisotropy).

Targets:
  * `Ua-1-0-CY` (problem 13) — :math:`r_c = 1.72500292` mfp. Already
    cross-checked by Variant α at 8.5e-6 (see
    :mod:`tests.derivations.test_peierls_greens_function_cylinder_xverif_sood2003`).
  * `Ua-1-1-CY` (problem 36) — 1G linearly anisotropic.

TODO: populate after literature-researcher delivers Westfall-Metcalf
1973 + Sanchez-Ganapol 1983.
"""
from __future__ import annotations
