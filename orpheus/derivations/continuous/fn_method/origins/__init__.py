"""Branch-1 SymPy derivations for F_N method benchmarks.

This subpackage carries the canonical algebra-of-record for the
LA-13511 closed-form k_inf cases. Each ``derive_*()`` function
returns a dict with PASS flags whose foundations are verified by
:mod:`tests.derivations.test_fn_la13511_kinf`.

The naming convention mirrors the existing peierls_greens_function
:mod:`origins` layout: identities are ``derive_<topic>_<form>()``,
each pinned by exactly one ``@pytest.mark.foundation`` test.

References
----------

* Sood, Forster, Parsons (1999) LA-13511. Sections V + Appendix A.
"""
from __future__ import annotations

from .k_inf_derivations import (
    derive_kinf_1g_eq_19,
    derive_kinf_1g_eq_20_simplifies_to_eq_19,
    derive_kinf_2g_general_from_matrix,
    derive_kinf_2g_no_upscatter,
    derive_phi_ratio_2g_no_upscatter,
    derive_kinf_mg_matrix_form,
    derive_kinf_mg_reduces_to_1g,
)

__all__ = [
    "derive_kinf_1g_eq_19",
    "derive_kinf_1g_eq_20_simplifies_to_eq_19",
    "derive_kinf_2g_general_from_matrix",
    "derive_kinf_2g_no_upscatter",
    "derive_phi_ratio_2g_no_upscatter",
    "derive_kinf_mg_matrix_form",
    "derive_kinf_mg_reduces_to_1g",
]
