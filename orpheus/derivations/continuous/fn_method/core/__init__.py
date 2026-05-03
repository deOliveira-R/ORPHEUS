"""F_N method shared primitives.

Public functions:

* :func:`case_dispersion_root_subcritical` — solve :math:`1 - c \\nu \\,
  \\mathrm{atanh}(1/\\nu) = 0` for :math:`\\nu_0 > 1` at :math:`c < 1`.
* :func:`case_dispersion_root_supercritical` — solve :math:`1 - c u \\,
  \\mathrm{atan}(1/u) = 0` for :math:`u_0 > 0` at :math:`c > 1`.
* :func:`case_nu0` — unified dispatcher returning ``(|\\nu_0|, is_imag)``.
* :func:`B_alpha`, :func:`A_alpha` — slab F_N moment integrals (single).
* :func:`B_alpha_array`, :func:`A_alpha_array` — vectorised over alpha.
* :func:`collocation_points` — Grandjean-Siewert prescription
  :math:`\\xi_0 = \\nu_0,\\ \\xi_1 = 0,\\ \\xi_2 = 1,\\ \\xi_{3..}` equally
  spaced in :math:`(0, 1)`.
* :func:`assemble_fn_matrix` — slab/sphere F_N matrix assembler
  parametrised by ``geometry_sign`` (Siewert-Thomas 1986 unification).
* :func:`fn_determinant` — convenience ``det(assemble_fn_matrix(...))``.

References
----------

* Siewert & Benoist 1979, *Nucl. Sci. Eng.* **69**, 156-160.
* Grandjean & Siewert 1979, *Nucl. Sci. Eng.* **69**, 161-168.
* Siewert & Thomas 1986, *Nucl. Sci. Eng.* **94**, 264.
"""
from __future__ import annotations

from .dispersion import (
    case_dispersion_root_subcritical,
    case_dispersion_root_supercritical,
    case_nu0,
)
from .fn_matrix import (
    assemble_fn_matrix,
    fn_determinant,
)
from .moments import (
    A_alpha,
    A_alpha_array,
    B_alpha,
    B_alpha_array,
    collocation_points,
)

__all__ = [
    "case_dispersion_root_subcritical",
    "case_dispersion_root_supercritical",
    "case_nu0",
    "A_alpha",
    "A_alpha_array",
    "B_alpha",
    "B_alpha_array",
    "collocation_points",
    "assemble_fn_matrix",
    "fn_determinant",
]
