"""Analytical and semi-analytical reference solutions for verifying ORPHEUS solvers.

SymPy and mpmath derivations that produce **continuous, mesh-independent
reference solutions** for verifying numerical solvers against the exact
form of the equation each solver discretises:

- Differential transport (SN, MOC) → analytical :math:`\\psi(x)` /
  :math:`\\phi(x)` via Method of Manufactured Solutions.
- Diffusion → closed-form eigenfunctions (sine, Bessel) or
  piecewise-smooth transfer-matrix solutions.
- Integral transport (CP) → semi-analytical Peierls solutions via
  high-precision Nyström quadrature of the exact
  :math:`E_n` / :math:`\\mathrm{Ki}_n` kernels.
- Stochastic transport (MC) → Case/Placzek, Milne, Sood analytical
  reference solutions from the transport literature.

Legacy scalar ``k_inf``-only :class:`VerificationCase` objects are
retrieved via :func:`get` / :func:`all_names`. The new Phase-0
:class:`ContinuousReferenceSolution` objects are retrieved via
:func:`continuous_get` / :func:`continuous_all_names`.

See :doc:`/verification/reference_solutions` for the contract these
references commit to, and the verification-campaign migration plan.
"""

from ._reference import (
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
