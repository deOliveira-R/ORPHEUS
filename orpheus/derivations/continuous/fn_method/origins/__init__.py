"""Branch-1 SymPy derivations for F_N method benchmarks.

This subpackage carries the canonical algebra-of-record for the
LA-13511 closed-form k_inf cases and the F_N method's published
recursion + critical-condition identities.

Naming convention: ``derive_<topic>_<form>()`` returns a dict with a
PASS flag whose foundation is verified by a 1:1 corresponding test in
:mod:`tests.derivations`.

Modules
-------

* :mod:`.k_inf_derivations` — LA-13511 Eqs 18-32, 72-76 (k_inf).
* :mod:`.fn_slab_derivations` — Slab F_N moment recursions
  (Siewert-Benoist Part I + Grandjean-Siewert Part II).
* :mod:`.fn_sphere_derivations` — Sphere F_N specialisation (Siewert-Thomas
  1986: slab + sphere unified by geometry_sign ∈ {+1, -1}).

References
----------

* Sood, Forster, Parsons (1999) LA-13511. Sections V + Appendix A.
* Siewert & Benoist 1979, *Nucl. Sci. Eng.* **69**, 156.
* Grandjean & Siewert 1979, *Nucl. Sci. Eng.* **69**, 161.
* Siewert & Thomas 1986, *Nucl. Sci. Eng.* **94**, 264.
* Kaper, Lindeman & Leaf 1974, *Nucl. Sci. Eng.* **54**, 94.
"""
from __future__ import annotations

from .fn_slab_derivations import (
    derive_A0_seed,
    derive_A_recursion,
    derive_B0_seed,
    derive_B_recursion,
    derive_critical_determinant_structure,
)
from .fn_sphere_derivations import (
    derive_sphere_2g_to_1g_reduction,
    derive_sphere_bc_sign_flip,
    derive_sphere_critical_condition,
    derive_sphere_fn_matrix_entry,
    derive_x_function_geometry_independence,
)
from .k_inf_derivations import (
    derive_kinf_1g_eq_19,
    derive_kinf_1g_eq_20_simplifies_to_eq_19,
    derive_kinf_2g_general_from_matrix,
    derive_kinf_2g_no_upscatter,
    derive_kinf_mg_matrix_form,
    derive_kinf_mg_reduces_to_1g,
    derive_phi_ratio_2g_no_upscatter,
)

__all__ = [
    # k_inf
    "derive_kinf_1g_eq_19",
    "derive_kinf_1g_eq_20_simplifies_to_eq_19",
    "derive_kinf_2g_general_from_matrix",
    "derive_kinf_2g_no_upscatter",
    "derive_phi_ratio_2g_no_upscatter",
    "derive_kinf_mg_matrix_form",
    "derive_kinf_mg_reduces_to_1g",
    # F_N slab
    "derive_B_recursion",
    "derive_A_recursion",
    "derive_B0_seed",
    "derive_A0_seed",
    "derive_critical_determinant_structure",
    # F_N sphere
    "derive_sphere_bc_sign_flip",
    "derive_sphere_fn_matrix_entry",
    "derive_sphere_critical_condition",
    "derive_sphere_2g_to_1g_reduction",
    "derive_x_function_geometry_independence",
]
