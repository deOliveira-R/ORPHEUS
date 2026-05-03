r"""Case singular-eigenfunction method for one-speed transport with
linearly anisotropic scattering — Atalay 1997 family.

This package is a structurally-independent verification ground for
ORPHEUS transport solvers in the **reflected slab + sphere with
linearly anisotropic scattering** corner of the Sood/Forster/Parsons
LA-13511 catalogue. It mounts on the published Atalay 1997 method
[Atalay1997]_, which extends the McCormick-Kušcer 1965 half-range
bi-orthogonality machinery with **four new parallel half-range
relations** (Atalay Eqs. 28-31) that close the deficit blocking
previous reflected-slab attempts.

Mathematical pillar
-------------------

This is **NOT** the F_N method (collocation in the matrix algebra).
It is the **Case singular-eigenfunction expansion** + **first-order
Fredholm iteration** branch — a different mathematical pillar
sharing only ``numpy`` / ``scipy.special`` / ``mpmath`` (above the
trusted-library line per the algebra-of-record discipline).

Slab + sphere unification
-------------------------

Atalay reuses the slab kernels for the sphere via the standard
Case-Zweifel **antisymmetric trick**: the sphere flux corresponds to
the **odd-mode** slab solution. ψ(x,μ) = -ψ(-x,-μ) gives a₀₊ = -a₀₋
and A(ν) = -A(-ν). The structural change from slab to sphere is

* T(R,μ) → T₁(R,μ) (sign flip in the :math:`Re^{d/μ} ± e^{-d/μ}` numerator).
* K_j → L_j (same kernel structure, different T → T₁).
* sin↔cos shuffle in the criticality condition.

This is the same structural pattern as Siewert-Thomas 1986 sphere F_N
sharing slab kernels via a ``geometry_sign`` parameter — see
:mod:`orpheus.derivations.continuous.fn_method.core.fn_matrix`.

Validity range
--------------

The bi-orthogonality relations are valid for :math:`c \le 1 + 1/(3 f_1)`
(Atalay Eq. 5), the range where the transport operator has only one
pair of discrete modes :math:`±\nu_0`. For :math:`f_1 = 0` (isotropic
limit) the bound is trivial (all :math:`c`); for :math:`f_1 = 0.30` the
upper limit is :math:`c \le 2.111`. Outside this band complex
eigenvalues appear (Dahl-Sjöstrand 1979, Kohut 1993) — first-order
Atalay does not detect them.

References
----------

.. [Atalay1997] Atalay, M.A. (1997).
   "The reflected slab and sphere criticality problem with anisotropic
   scattering in one-speed neutron transport theory."
   *Progress in Nuclear Energy* **31**(3), 229-252.
   DOI: 10.1016/0149-1970(95)00094-1.

.. [BurkartIshiguroSiewert1976] Burkart, A.R., Ishiguro, Y., Siewert, C.E. (1976).
   "Neutron transport in two dissimilar media with anisotropic scattering."
   *Nuclear Science and Engineering* **61**, 72-81.
   (Used as cross-check for Atalay vacuum-slab + anisotropic results.)

.. [McCormickKuscer1965] McCormick, N.J., Kušcer, I. (1965).
   "Bi-orthogonality relations for solving half-space transport problems."
   *J. Math. Phys.* **6**, 1939.
"""
from __future__ import annotations
