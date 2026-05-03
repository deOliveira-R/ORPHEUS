r"""Westfall-Metcalf 1973 singular-eigenfunction expansion of the
infinite-cylinder bare-critical problem.

Public API
----------

* :class:`CylinderSingularEigenfunctionResult` — container for the
  bare-critical solver output (critical radius, dominant Case
  eigenvalue, root-find diagnostics, callable flux reconstructions).
* :func:`solve_singular_eigenfunction_cylinder_bare_critical` — solves
  the criticality condition (V_se-cyl.5) for the bare-critical infinite
  cylinder of given material parameter :math:`c > 1`.
"""
from __future__ import annotations

from .one_group import (
    CylinderSingularEigenfunctionResult,
    solve_singular_eigenfunction_cylinder_bare_critical,
)

__all__ = [
    "CylinderSingularEigenfunctionResult",
    "solve_singular_eigenfunction_cylinder_bare_critical",
]
