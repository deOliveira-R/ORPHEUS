"""Reference solutions and symbolic discretisations for ORPHEUS verification.

This package is organised into three top-level paths:

- :mod:`~orpheus.derivations.common` — shared math root and utilities
  (kernels, quadrature, cross-section library, eigenvalue helpers,
  the :class:`VerificationCase` /
  :class:`ContinuousReferenceSolution` types). The symbolic transport
  equation that both verification paths branch from lives in
  :mod:`~orpheus.derivations.common.transport_equation`.

- :mod:`~orpheus.derivations.discrete` — *Path 1*: symbolic
  discretisations of the production solvers
  (S\\ :sub:`N` balance / contamination, MOC characteristics, …).
  These describe the equations the solvers commit to satisfying.

- :mod:`~orpheus.derivations.continuous` — *Path 2*: continuous,
  mesh-independent reference solutions
  (closed-form analytical, flat-source CP, Peierls integral form,
  Method of Manufactured Solutions, transfer-matrix diffusion,
  Case/Placzek MC analytics, …). These are what the production
  solvers are verified against.

Two retrieval registries live at the top of this package:

- :func:`get` / :func:`all_names` — legacy scalar ``k_inf``
  :class:`VerificationCase` lookup.
- :func:`continuous_get` / :func:`continuous_all_names` —
  :class:`ContinuousReferenceSolution` lookup, the Phase-0
  contract that the verification campaign is migrating to.

See :doc:`/verification/reference_solutions` for the contract these
references commit to and the verification-campaign migration plan,
and ``orpheus/derivations/README.md`` for the architectural
treatment of the two paths.
"""

from .common.continuous_reference import (
    ContinuousReferenceSolution,
    OperatorForm,
    ProblemSpec,
    Provenance,
)
from .reference_values import (
    all_names,
    continuous_all,
    continuous_all_names,
    continuous_by_operator_form,
    continuous_get,
    continuous_register,
    get,
)
