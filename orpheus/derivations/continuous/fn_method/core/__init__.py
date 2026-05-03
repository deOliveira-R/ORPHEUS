"""F_N method shared primitives — placeholder.

This subpackage will house the F_N method machinery shared across
slab/sphere/cylinder geometries:

* ``case_eigenvalues`` — :math:`\\nu_0, \\nu_1` transcendental roots of
  the dispersion law :math:`1 - c\\nu \\, \\mathrm{atanh}(1/\\nu) = 0`
  for :math:`c > 1`; SymPy preamble + mpmath findroot.
* ``moments`` — :math:`\\int \\mu^k / (\\nu \\pm \\mu) \\, \\mathrm d\\mu`
  integrals (slab moment library — closed-form rational + log).
* ``h_function`` — Chandrasekhar H/X half-range factorisation for
  Wiener-Hopf BC handling.
* ``bessel_kernels`` — :math:`\\mathrm{Ki}_n, K_n` integrals for
  cylinder F_N (Westfall-Metcalf 1973, Sanchez-Ganapol 1983).

TODO: populate when Kaper-Lindeman-Leaf 1974 NSE 54, 94 and
Westfall-Metcalf 1973 NSE 52, 1 are acquired (literature-researcher
dispatched in parallel with the first-slice k_inf work).
"""
from __future__ import annotations
