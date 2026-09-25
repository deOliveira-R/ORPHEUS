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
:mod:`orpheus.derivations.continuous.trajectory_resolvent` hosts
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

Withdrawn (#506)
----------------

The Nyström solver half of this package — the volume kernel, the
boundary-closure operators, ``solve_peierls_1g`` / ``solve_peierls_mg``,
the slab eigen-solve, the ``_build_peierls_*_case`` builders and the
case builders the lazy registry calls — is withdrawn by the maintainer's
ruling of 2026-09-24: it is not research grade, so it is neither cited as
evidence nor run by default. Each of those symbols is locked with
:data:`PEIERLS_NYSTROM_WITHDRAWAL` through
:func:`~orpheus.derivations.common.withdrawal.withdrawn_generator`; a call
raises :class:`~orpheus.derivations.common.withdrawal.GeneratorWithdrawn`
unless ``ORPHEUS_RUN_WITHDRAWN=506`` is set. The escape, transmission and
boundary primitives (``compute_P_*``, ``compute_G_*``, ``compute_T_*``),
their angular-assembly drivers, :mod:`.reference` and :mod:`.ps1982_reference`
stay in service.
"""

from orpheus.derivations.common.withdrawal import Withdrawal

#: The ruling that withdraws the Nyström solver half of this package; the
#: one value every lock in the package is keyed on, and the value the one
#: test-side mark (``tests/_harness/withdrawals.py``) is minted from.
PEIERLS_NYSTROM_WITHDRAWAL = Withdrawal(
    reason=(
        "Peierls Nyström reference is not research grade "
        "(maintainer ruling 2026-09-24)"
    ),
    issue=506,
)
