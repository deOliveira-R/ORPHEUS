r"""Linear-anisotropic-scattering utilities for the Case method.

The Atalay 1997 method handles linearly anisotropic scattering only
(P_1 expansion: :math:`\Sigma_s(\mu \to \mu') = (\Sigma_s/2)(1 + 3 f_1
\mu \mu')`, mean cosine :math:`f_1`).
"""
from __future__ import annotations

from .linear import (
    atalay_validity_max_c,
    check_atalay_validity,
)

__all__ = [
    "atalay_validity_max_c",
    "check_atalay_validity",
]
