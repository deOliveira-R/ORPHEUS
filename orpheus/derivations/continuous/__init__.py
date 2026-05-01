"""Path 2 — continuous reference derivations.

Each sub-package commits to a *continuous, mesh-independent* form of
the transport equation that the production solvers in
:mod:`orpheus.derivations.discrete` are verified against. The
references here use analytical, semi-analytical, or transcendental
techniques (closed-form algebra, Peierls integral form via Nyström
collocation, Method of Manufactured Solutions, transfer-matrix
diffusion, etc.).

Sub-packages:

- :mod:`~orpheus.derivations.continuous.analytical` — closed-form
  homogeneous-medium and infinite-medium spectral references.
- :mod:`~orpheus.derivations.continuous.flat_source_cp` —
  flat-source collision-probability eigenvalues for slab, cylinder,
  sphere.
- :mod:`~orpheus.derivations.continuous.peierls` — Peierls integral
  form via high-precision Nyström quadrature; the active reference
  for heterogeneous CP verification.
- :mod:`~orpheus.derivations.continuous.greens_function` — placeholder
  for analytical Green's-function references (Plan-2).
- :mod:`~orpheus.derivations.continuous.mms` — Method of Manufactured
  Solutions for S\\ :sub:`N` and MOC.
- :mod:`~orpheus.derivations.continuous.cases` — per-method case
  registries that compose the derivations above.
"""
