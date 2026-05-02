"""Peierls Variant α Green's-function references.

This sub-package hosts the angle-resolved Green's-function references
produced via the Peierls Variant α formulation (Plan-2). The Variant α
form interpolates between vacuum (:math:`\\alpha = 0`) and specular
(:math:`\\alpha = 1`) boundary conditions on a single rank-1 resolvent,
sharing the closure operator across geometry kernels.

The package is the continuous-reference twin of
:mod:`orpheus.derivations.continuous.peierls_nystrom`: both descend from
the Peierls integral form, but the Nyström family discretizes the
*scalar-flux* integral equation while the Variant α family discretizes
the *angle-resolved* Green's function and reduces it analytically via
the rank-1 trick.

Sub-modules:

- :mod:`~orpheus.derivations.continuous.peierls_greens_function.variant_alpha_core` —
  the shared rank-1 Variant α resolvent and surface-flux closure
  primitives consumed by both sphere and cylinder solvers.
- :mod:`~orpheus.derivations.continuous.peierls_greens_function.greens_function` —
  sphere production solver (fixed-source Green's function and
  multi-region eigenvalue/k-inf evaluation).
- :mod:`~orpheus.derivations.continuous.peierls_greens_function.greens_function_cylinder` —
  cylinder production solver (mirror of the sphere driver, mounted on
  the same Variant α core).
- :mod:`~orpheus.derivations.continuous.peierls_greens_function.origins` —
  symbolic *origins* (SymPy V_α1..V_α3 operator-level identities) for
  sphere and cylinder.
"""
