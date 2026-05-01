"""Peierls integral form — high-precision Nyström references.

The Peierls form rewrites the linear Boltzmann equation as an
integral equation in the *scalar flux* alone (the angular flux is
eliminated by analytic integration of the streaming operator with
the appropriate single-event escape kernel). The eigenvalue and
flux-shape references in this sub-package are produced by Nyström
collocation of the integral form on a high-precision quadrature
grid; they are the reference of choice for heterogeneous CP
verification.

Sub-modules:

- :mod:`~orpheus.derivations.continuous.peierls.geometry` —
  shared geometry primitives (chord half-lengths, panel quadratures,
  closure operators).
- :mod:`~orpheus.derivations.continuous.peierls.slab` /
  :mod:`~.cylinder` /
  :mod:`~.sphere` — geometry-specific solvers and registries.
- :mod:`~orpheus.derivations.continuous.peierls.cases` — multi-region
  / multi-group case manifest.
- :mod:`~orpheus.derivations.continuous.peierls.reference` — entry
  point that exposes the Peierls-form references to the registry.
- :mod:`~orpheus.derivations.continuous.peierls.origins` — symbolic
  *origins* (specular BC R-matrix, cylindrical 3-D G-BC, Knyazev
  shifted-Legendre identities) that are imported here without
  having a continuous reference of their own.
"""
