"""Shared mathematical foundations for the derivations package.

This sub-package houses the **mathematical root** that both the
``discrete`` and ``continuous`` paths branch from, plus the cross-cutting
utilities (kernels, quadrature, cross-section library, eigenvalue
helpers, verification-case dataclass) that no single path owns.

The four utilities most consumers import directly:

- :mod:`~orpheus.derivations.common.kernels` —
  :math:`E_n`, :math:`\\mathrm{Ki}_n`, chord primitives.
- :mod:`~orpheus.derivations.common.quadrature` — 1-D quadrature
  contract (:class:`Quadrature1D`, :class:`AdaptiveQuadrature1D`).
- :mod:`~orpheus.derivations.common.quadrature_recipes` —
  geometry-aware recipes (chord, observer-angular).
- :mod:`~orpheus.derivations.common.xs_library` — synthetic
  cross-section library shared between derivations and tests.

The entry point for the math root is
:mod:`~orpheus.derivations.common.transport_equation`, which lays
out the symbolic transport equation that both Path 1 (discrete
production solvers) and Path 2 (continuous reference derivations)
discretise.
"""
