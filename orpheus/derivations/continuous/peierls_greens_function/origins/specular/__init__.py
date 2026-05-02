"""Specular boundary-condition derivations for the Peierls Variant α
Green's-function family.

Sub-package shell for the SymPy math-origin functions that produce the
operator-level identities (V_α1..V_α3) for the angle-resolved Green's
function reference. Public symbols are re-exported here so call sites
can import directly from
``orpheus.derivations.continuous.peierls_greens_function.origins.specular``.

Modules
-------

- :mod:`.greens_function` — sphere V_α identities (closed-sphere
  fixed-point, :math:`T_{00}^{\\rm sphere} = P_{ss}^{\\rm sphere}`,
  :math:`\\alpha = 0` kernel reduction).
- :mod:`.greens_function_cylinder` — cylinder V_α identities (closed-
  cylinder fixed-point, :math:`T_{00}^{\\rm cyl} = P_{ss}^{\\rm cyl}`,
  :math:`\\alpha = 0` kernel reduction, bounce-period chord identity).
- :mod:`.greens_function_slab` — slab V_α identities (closed-slab
  symmetric-specular fixed-point, :math:`T_{00}^{\\rm slab} =
  P_{ss}^{\\rm slab} = 2 E_3(\\Sigma_t L)`, :math:`\\alpha = 0` kernel
  reduction). Slab is the first 2-bounce-per-period geometry.
- :mod:`.greens_function_slab_asymmetric` — asymmetric slab V_α
  identities for the **rank-2 boundary-to-boundary scattering
  resolvent** :math:`T = (I - S)^{-1}` (closed asymmetric-slab fixed-
  point at :math:`\\alpha_L = \\alpha_R = 1`, rank-2 inversion algebra
  with reductions to symmetric / one-vacuum-wall / vacuum-vacuum
  corners, :math:`\\alpha_L = \\alpha_R = 0` kernel reduction). First
  validated instance of the cross-domain frame's rank-2 prediction.
- :mod:`.greens_function_hollow_sphere` — hollow sphere V_α identities
  with rank-2 closure on the through-ray subset (:math:`b \\le R_{\\rm
  in}`) and rank-1 closure on the outer-only subset (:math:`b > R_{\\rm
  in}`). First curvilinear 2-surface instance of the rank-2 frame.
"""

from .greens_function import (
    derive_T00_equals_P_ss_sphere,
    derive_alpha_zero_kernel_reduction,
    derive_operator_constant_trial_closed_sphere,
)
from .greens_function_cylinder import (
    derive_T00_equals_P_ss_cylinder,
    derive_alpha_zero_kernel_reduction_cylinder,
    derive_bounce_period_chord_cylinder,
    derive_operator_constant_trial_closed_cylinder,
)
from .greens_function_hollow_sphere import (
    derive_alpha_zero_kernel_reduction_hollow_sphere,
    derive_operator_constant_trial_closed_hollow_sphere,
    derive_rank2_resolvent_hollow_sphere,
)
from .greens_function_slab import (
    derive_T00_equals_P_ss_slab,
    derive_alpha_zero_kernel_reduction_slab,
    derive_operator_constant_trial_closed_slab,
)
from .greens_function_slab_asymmetric import (
    derive_alpha_zero_kernel_reduction_slab_asymmetric,
    derive_operator_constant_trial_closed_slab_asymmetric,
    derive_rank2_resolvent_slab_asymmetric,
)

__all__ = [
    # greens_function (sphere)
    "derive_T00_equals_P_ss_sphere",
    "derive_alpha_zero_kernel_reduction",
    "derive_operator_constant_trial_closed_sphere",
    # greens_function_cylinder
    "derive_T00_equals_P_ss_cylinder",
    "derive_alpha_zero_kernel_reduction_cylinder",
    "derive_bounce_period_chord_cylinder",
    "derive_operator_constant_trial_closed_cylinder",
    # greens_function_slab
    "derive_T00_equals_P_ss_slab",
    "derive_alpha_zero_kernel_reduction_slab",
    "derive_operator_constant_trial_closed_slab",
    # greens_function_slab_asymmetric (rank-2)
    "derive_alpha_zero_kernel_reduction_slab_asymmetric",
    "derive_operator_constant_trial_closed_slab_asymmetric",
    "derive_rank2_resolvent_slab_asymmetric",
    # greens_function_hollow_sphere (rank-2 + impact-parameter partition)
    "derive_alpha_zero_kernel_reduction_hollow_sphere",
    "derive_operator_constant_trial_closed_hollow_sphere",
    "derive_rank2_resolvent_hollow_sphere",
]
