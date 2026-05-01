"""Specular boundary-condition derivations for the Peierls path.

Sub-package shell for the SymPy math-origin functions that produce
the curvilinear specular reflection operator and the slab per-face
companion. Public symbols are re-exported here so existing call sites
that imported from the flat ``orpheus.derivations.continuous.peierls.origins.specular``
module continue to resolve.

Modules
-------

- :mod:`.r_matrix` — rank-:math:`N` curvilinear reflection operator
  :math:`R_{\\rm spec} = \\tfrac{1}{2}\\,M^{-1}` (the load-bearing
  identity).
- :mod:`.slab` — per-face slab primitives (block-diagonal
  :math:`R_{\\rm slab} = {\\rm blockdiag}(R_{\\rm spec}, R_{\\rm spec})`
  and the :math:`E_{k+2}` closed-form mode primitives).
- :mod:`.continuous_mu` — Phase-5 Sanchez 1986 Eq. (A6)
  continuous-:math:`\\mu` multi-bounce kernel verifications (V1..V4
  identities for the textbook reference implementation).
- :mod:`.greens_function` — Plan-2 Variant α operator-level identities
  (V_α1..V_α3) for the angle-resolved Green's function reference.
"""

from .r_matrix import (
    build_M_closed_form,
    build_M_symbolic,
    build_R_specular_symbolic,
)
from .slab import (
    build_R_slab_blockdiag,
    derive_slab_g_outer_n0,
    derive_slab_p_outer_n0,
    slab_g_outer_mode_n_no_mu,
    slab_p_outer_mode_n_no_mu,
)
from .continuous_mu import (
    derive_diagonal_singularity,
    derive_m1_equivalence,
    derive_multi_bounce_factor,
    derive_vacuum_reduction,
)
from .greens_function import (
    derive_T00_equals_P_ss_sphere,
    derive_alpha_zero_kernel_reduction,
    derive_operator_constant_trial_closed_sphere,
)

__all__ = [
    # r_matrix
    "build_M_closed_form",
    "build_M_symbolic",
    "build_R_specular_symbolic",
    # slab
    "build_R_slab_blockdiag",
    "derive_slab_g_outer_n0",
    "derive_slab_p_outer_n0",
    "slab_g_outer_mode_n_no_mu",
    "slab_p_outer_mode_n_no_mu",
    # continuous_mu
    "derive_diagonal_singularity",
    "derive_m1_equivalence",
    "derive_multi_bounce_factor",
    "derive_vacuum_reduction",
    # greens_function
    "derive_T00_equals_P_ss_sphere",
    "derive_alpha_zero_kernel_reduction",
    "derive_operator_constant_trial_closed_sphere",
]
