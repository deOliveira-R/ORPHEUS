"""Peierls integral form — high-precision Nyström references.

The Peierls form rewrites the linear Boltzmann equation as an
integral equation in the *scalar flux* alone (the angular flux is
eliminated by analytic integration of the streaming operator with
the appropriate single-event escape kernel). The eigenvalue and
flux-shape references in this sub-package are produced by Nyström
collocation of the integral form on a high-precision quadrature
grid; they are the reference of choice for heterogeneous CP
verification.

The companion sub-package
:mod:`orpheus.derivations.continuous.peierls_greens_function` hosts
the angle-resolved Variant α Green's-function references that share
the Peierls integral ancestry but discretize the angle-resolved
Green's function rather than the scalar-flux integral equation.

Sub-modules:

- :mod:`~orpheus.derivations.continuous.peierls_nystrom.geometry` —
  shared geometry primitives (chord half-lengths, panel quadratures,
  closure operators).
- :mod:`~orpheus.derivations.continuous.peierls_nystrom.slab` /
  :mod:`~.cylinder` /
  :mod:`~.sphere` — geometry-specific solvers and registries.
- :mod:`~orpheus.derivations.continuous.peierls_nystrom.cases` — multi-region
  / multi-group case manifest.
- :mod:`~orpheus.derivations.continuous.peierls_nystrom.reference` — entry
  point that exposes the Peierls-form references to the registry.
- :mod:`~orpheus.derivations.continuous.peierls_nystrom.origins` — symbolic
  *origins* (specular BC R-matrix, cylindrical 3-D G-BC, Knyazev
  shifted-Legendre identities) that are imported here without
  having a continuous reference of their own.
"""
