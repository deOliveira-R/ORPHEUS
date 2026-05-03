r"""Atalay 1997 sphere criticality solver — linearly anisotropic, reflected.

Implements the criticality condition Atalay Eq 54 — the **odd-mode**
sphere criticality from the antisymmetric BC :math:`\psi(x, \mu) =
-\psi(-x, -\mu)` (Eq 47). Reuses 100% of slab kernel machinery; the
only structural change is :math:`T \to T_1` (sign flip in the second
exponential of the T-function numerator and denominator) and the
sin↔cos shuffle in the LHS form.

Public API
----------

* :class:`CaseMethodSphereResult` — result container.
* :func:`solve_case_method_sphere_critical` — given :math:`(c, R, f_1)`,
  finds critical sphere radius :math:`R_s` in mfp.
"""
from __future__ import annotations

from .one_group import (
    CaseMethodSphereResult,
    solve_case_method_sphere_critical,
)

__all__ = [
    "CaseMethodSphereResult",
    "solve_case_method_sphere_critical",
]
