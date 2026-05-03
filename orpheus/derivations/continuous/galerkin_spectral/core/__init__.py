r"""Shared Carlvik-Galerkin numerical primitives.

* :mod:`.carlvik_recurrences` — closed-form / quadrature evaluation
  of :math:`A_{m,n}` and :math:`B_{m,n}` matrix elements.
* :mod:`.galerkin_matrix` — assembly of the Eq. (3) :math:`(A - 3\bar\mu(c-1)B)`
  matrix and Eq. (4) block-matrix form, plus eigenvalue extraction.
"""
from __future__ import annotations
