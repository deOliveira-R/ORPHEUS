"""Cross-method regression infrastructure.

Lightweight shared protocol for verifying that multiple solvers (each
of which is independently L1-verified against analytical truth) agree
on a single physical case to within their respective quadrature
floors.

The package contains:

* :mod:`tests.cross_method.protocol` — :class:`CrossMethodCase`,
  :class:`SolverAdapter` Protocol, :class:`ScalarResult`, and the
  ``agree`` tolerance helper.
* :mod:`tests.cross_method.adapters` — concrete adapters for each
  shipped continuous reference solver (fn_method, trajectory_resolvent,
  ...). New solvers register here.
* :mod:`tests.cross_method.cases` — populated case sets keyed by
  geometry (bare-critical slab / sphere; reflected slab; multi-region
  fixed-source). Each case carries truth + per-solver tolerances +
  pillar tag.
* ``test_eigenvalue.py`` — parametrised cross-method gates.
* ``test_fixed_source.py`` — fixed-source one-sided coverage cases.

See :doc:`/testing/cross_method` for the architecture rationale.
"""
