"""Branch-2 numpy/scipy implementations for the LA-13511 :math:`k_\\infty` cases.

Public entry points:

* :func:`compute_kinf_1g` — Sood Eq 19/20 (1G infinite medium).
* :func:`compute_kinf_2g_no_upscatter` — Sood Eq 29 (2G, no upscatter).
* :func:`compute_kinf_2g_general` — Sood Eq 28 corrected form
  (general 2G, with upscatter possible).
* :func:`compute_kinf_mg` — Sood Eq 76 (general G).
* :func:`compute_flux_ratio_2g_no_upscatter` — Sood Eq 32
  (:math:`\\phi_2/\\phi_1` for no upscatter).

These are pure-numpy, structurally-independent of the production
:func:`orpheus.derivations.common.eigenvalue.kinf_homogeneous` solver:
the F_N-method module evaluates the published closed forms directly
(Sood Eq 28/29/32), while ``kinf_homogeneous`` solves the dominant-
eigenvalue problem of :math:`A^{-1} F` via :func:`numpy.linalg.eig`.
The two paths share only ``numpy`` itself (above the trusted-library
line) and must agree to floating-point precision on every case where
both apply.
"""
from __future__ import annotations

from .k_inf import (
    compute_kinf_1g,
    compute_kinf_2g_no_upscatter,
    compute_kinf_2g_general,
    compute_kinf_mg,
    compute_flux_ratio_2g_no_upscatter,
)

__all__ = [
    "compute_kinf_1g",
    "compute_kinf_2g_no_upscatter",
    "compute_kinf_2g_general",
    "compute_kinf_mg",
    "compute_flux_ratio_2g_no_upscatter",
]
