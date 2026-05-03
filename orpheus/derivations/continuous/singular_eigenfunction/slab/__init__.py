r"""Atalay 1997 slab criticality solver — linearly anisotropic, reflected.

Implements the criticality condition Atalay Eq 46 for the
**even-mode** symmetric slab with linearly anisotropic scattering and
specular-like reflection coefficient :math:`R \in [0, 1)` (perfect
reflection :math:`R = 1` makes thickness drop out).

Public API
----------

* :class:`CaseMethodSlabResult` — result container.
* :func:`solve_case_method_slab_critical` — given :math:`(c, R, f_1)`,
  finds critical half-thickness :math:`d` in mfp.
"""
from __future__ import annotations

from .one_group import (
    CaseMethodSlabResult,
    solve_case_method_slab_critical,
)

__all__ = [
    "CaseMethodSlabResult",
    "solve_case_method_slab_critical",
]
